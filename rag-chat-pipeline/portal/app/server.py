#!/usr/bin/env python3
"""
RAG Chat Pipeline Web Portal
Flask backend server for managing the RAG chat workflow
"""

import os
import json
import subprocess
import threading
import time
import urllib.request
from pathlib import Path
from datetime import datetime
from flask import Flask, render_template, jsonify, request, Response, stream_with_context
from flask_cors import CORS
import psutil
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
try:
    import pynvml
    NVML_AVAILABLE = True
except ImportError:
    NVML_AVAILABLE = False

app = Flask(__name__,
            template_folder='../templates',
            static_folder='../static')
CORS(app)

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent.parent.absolute()
SCRIPTS_DIR = PROJECT_ROOT / 'scripts'
DATA_DIR = PROJECT_ROOT / 'data'
CONFIG_FILE = PROJECT_ROOT / '.env'
LOG_DIR = DATA_DIR / 'logs'
VCF_PREVIEW_CACHE = {}

# Global state
pipeline_state = {
    'current_step': None,
    'status': 'idle',  # idle, running, success, error
    'progress': 0,
    'logs': [],
    'start_time': None,
    'end_time': None,
    'error_message': None
}

# LLM Performance metrics
llm_metrics = {
    'ttft': None,  # Time to first token (ms)
    'tokens_per_sec': None,
    'cache_hits': 0,
    'cache_misses': 0,
    'total_queries': 0,
    'avg_latency': None,
}

running_processes = {}


def load_config():
    """Load pipeline configuration from .env"""
    config = {}
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    config[key.strip()] = value.strip()
    return config


def save_config(config):
    """Save configuration to .env"""
    lines = []
    existing_keys = set()

    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    lines.append(line)
                elif '=' in stripped:
                    key = stripped.split('=', 1)[0].strip()
                    existing_keys.add(key)
                    if key in config:
                        lines.append(f"{key}={config[key]}\n")
                    else:
                        lines.append(line)

    # Add any new config values
    for key, value in config.items():
        if key not in existing_keys:
            lines.append(f"{key}={value}\n")

    with open(CONFIG_FILE, 'w') as f:
        f.writelines(lines)


def check_file_exists(filepath):
    """Check if a file exists"""
    return Path(filepath).exists()


def get_file_size(filepath):
    """Get file size in human-readable format"""
    path = Path(filepath)
    if not path.exists():
        return "N/A"

    size = path.stat().st_size
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} PB"


def get_disk_space():
    """Get available disk space"""
    try:
        result = subprocess.run(
            ['df', '-h', str(PROJECT_ROOT)],
            capture_output=True,
            text=True
        )
        lines = result.stdout.strip().split('\n')
        if len(lines) > 1:
            parts = lines[1].split()
            return {
                'total': parts[1],
                'used': parts[2],
                'available': parts[3],
                'use_percent': parts[4]
            }
    except Exception as e:
        print(f"Error getting disk space: {e}")

    return {'total': 'N/A', 'used': 'N/A', 'available': 'N/A', 'use_percent': 'N/A'}


def get_vcf_preview(vcf_path: str, limit: int = 100):
    """Get VCF preview data (cached)"""
    global VCF_PREVIEW_CACHE

    cache_key = f"{vcf_path}:{limit}"
    if cache_key in VCF_PREVIEW_CACHE:
        return VCF_PREVIEW_CACHE[cache_key]

    variants = []
    try:
        # Try using cyvcf2 if available
        try:
            from cyvcf2 import VCF
            vcf = VCF(vcf_path)
            for i, record in enumerate(vcf):
                if i >= limit:
                    break
                variants.append({
                    'chrom': record.CHROM,
                    'pos': record.POS,
                    'ref': record.REF[:20] + ('...' if len(record.REF) > 20 else ''),
                    'alt': ','.join([str(a)[:20] + ('...' if len(str(a)) > 20 else '') for a in record.ALT]),
                    'qual': f"{record.QUAL:.1f}" if record.QUAL else 'N/A',
                    'filter': record.FILTER or 'PASS',
                })
            vcf.close()
        except ImportError:
            # Fallback to gzip/pysam or raw parsing
            import gzip
            opener = gzip.open if vcf_path.endswith('.gz') else open
            with opener(vcf_path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if len(variants) >= limit:
                        break
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        variants.append({
                            'chrom': parts[0],
                            'pos': int(parts[1]),
                            'ref': parts[3][:20] + ('...' if len(parts[3]) > 20 else ''),
                            'alt': parts[4][:20] + ('...' if len(parts[4]) > 20 else ''),
                            'qual': parts[5] if parts[5] != '.' else 'N/A',
                            'filter': parts[6] if parts[6] != '.' else 'PASS',
                        })

        VCF_PREVIEW_CACHE[cache_key] = variants
    except Exception as e:
        print(f"Error reading VCF: {e}")
        return []

    return variants


def run_script(script_name, step_name):
    """Run a pipeline script"""
    global pipeline_state, running_processes

    script_path = SCRIPTS_DIR / script_name

    if not script_path.exists():
        pipeline_state['status'] = 'error'
        pipeline_state['error_message'] = f"Script not found: {script_name}"
        return

    pipeline_state['current_step'] = step_name
    pipeline_state['status'] = 'running'
    pipeline_state['start_time'] = datetime.now().isoformat()
    pipeline_state['error_message'] = None

    try:
        # Use stdbuf to force unbuffered output
        process = subprocess.Popen(
            ['stdbuf', '-oL', '-eL', 'bash', str(script_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=0,
            cwd=str(PROJECT_ROOT),
            env={**os.environ, 'PYTHONUNBUFFERED': '1'}
        )

        running_processes[step_name] = process

        # Stream output line by line
        for line in iter(process.stdout.readline, ''):
            if line:
                pipeline_state['logs'].append({
                    'timestamp': datetime.now().isoformat(),
                    'message': line.rstrip()
                })

        process.wait()

        if process.returncode == 0:
            pipeline_state['status'] = 'success'
        else:
            pipeline_state['status'] = 'error'
            pipeline_state['error_message'] = f"Script exited with code {process.returncode}"

    except Exception as e:
        pipeline_state['status'] = 'error'
        pipeline_state['error_message'] = str(e)

    finally:
        pipeline_state['end_time'] = datetime.now().isoformat()
        if step_name in running_processes:
            del running_processes[step_name]


def run_command(command, step_name, cwd=None):
    """Run a shell command"""
    global pipeline_state, running_processes

    pipeline_state['current_step'] = step_name
    pipeline_state['status'] = 'running'
    pipeline_state['start_time'] = datetime.now().isoformat()
    pipeline_state['error_message'] = None

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=0,
            cwd=cwd or str(PROJECT_ROOT),
            env={**os.environ, 'PYTHONUNBUFFERED': '1'},
            shell=True
        )

        running_processes[step_name] = process

        for line in iter(process.stdout.readline, ''):
            if line:
                pipeline_state['logs'].append({
                    'timestamp': datetime.now().isoformat(),
                    'message': line.rstrip()
                })

        process.wait()

        if process.returncode == 0:
            pipeline_state['status'] = 'success'
        else:
            pipeline_state['status'] = 'error'
            pipeline_state['error_message'] = f"Command exited with code {process.returncode}"

    except Exception as e:
        pipeline_state['status'] = 'error'
        pipeline_state['error_message'] = str(e)

    finally:
        pipeline_state['end_time'] = datetime.now().isoformat()
        if step_name in running_processes:
            del running_processes[step_name]


# Routes
@app.route('/health')
@app.route('/healthz')
def health():
    """Health check endpoint for monitoring systems"""
    return jsonify({'status': 'healthy', 'service': 'rag-portal', 'port': 5001})


@app.route('/')
def index():
    """Main page"""
    return render_template('index.html')


@app.route('/api/status')
def get_status():
    """Get current pipeline status"""
    config = load_config()

    # Check service status
    # Milvus health on port 9091, Ollama on 11434 (replaces vLLM for ARM64)
    services = {
        'milvus': check_service_health('milvus', 'localhost', 9091),
        'ollama': check_service_health('ollama', 'localhost', 11434),
    }

    # Check data files
    vcf_path = config.get('VCF_INPUT_PATH', str(PROJECT_ROOT.parent / 'genomics-pipeline' / 'data' / 'output' / 'HG002.genome.vcf.gz'))
    data_status = {
        'vcf_exists': check_file_exists(vcf_path),
        'vcf_path': vcf_path,
        'vcf_size': get_file_size(vcf_path),
    }

    # Check if collection exists in Milvus
    collection_status = check_milvus_collection()

    # System status
    system_status = {
        'docker_installed': check_docker(),
        'gpu_available': check_gpu(),
        'disk_space': get_disk_space(),
        'cpu_utilization': get_cpu_utilization(),
        'gpu_utilization': get_gpu_utilization(),
        'memory_utilization': get_memory_utilization()
    }

    return jsonify({
        'pipeline': pipeline_state,
        'services': services,
        'data': data_status,
        'collection': collection_status,
        'system': system_status,
        'llm_metrics': llm_metrics,
    })


@app.route('/api/vcf-preview')
def vcf_preview():
    """Get VCF file preview"""
    config = load_config()
    vcf_path = config.get('VCF_INPUT_PATH', str(PROJECT_ROOT.parent / 'genomics-pipeline' / 'data' / 'output' / 'HG002.genome.vcf.gz'))
    limit = request.args.get('limit', 100, type=int)

    if not check_file_exists(vcf_path):
        return jsonify({'error': 'VCF file not found', 'path': vcf_path}), 404

    variants = get_vcf_preview(vcf_path, limit)
    return jsonify({
        'path': vcf_path,
        'count': len(variants),
        'variants': variants
    })


@app.route('/api/config', methods=['GET', 'POST'])
def config():
    """Get or update configuration"""
    if request.method == 'GET':
        return jsonify(load_config())
    else:
        config_data = request.json
        save_config(config_data)
        return jsonify({'success': True})


@app.route('/api/run/<step>')
def run_step(step):
    """Run a specific workflow step"""
    global pipeline_state

    if pipeline_state['status'] == 'running':
        return jsonify({'error': 'A step is already running'}), 400

    # Reset state
    pipeline_state['logs'] = []
    pipeline_state['progress'] = 0

    step_mapping = {
        'setup': ('source venv/bin/activate 2>/dev/null || python3 -m venv venv && source venv/bin/activate && pip install -r requirements.txt', 'Python Setup'),
        'start-milvus': ('docker-compose up -d milvus attu', 'Start Milvus'),
        'start-vllm': ('docker-compose up -d vllm', 'Start vLLM'),
        'start-all': ('docker-compose up -d', 'Start All Services'),
        'stop': ('docker-compose down', 'Stop Services'),
        'ingest': ('source venv/bin/activate && python scripts/ingest_vcf.py', 'Ingest VCF'),
        'ingest-limit': ('source venv/bin/activate && python scripts/ingest_vcf.py --limit 10000', 'Ingest VCF (10k)'),
        'chat': ('source venv/bin/activate && python scripts/run_chat.py', 'Launch Chat'),
    }

    if step not in step_mapping:
        return jsonify({'error': 'Invalid step'}), 400

    command, step_name = step_mapping[step]

    # Run in background thread
    thread = threading.Thread(target=run_command, args=(command, step_name))
    thread.daemon = True
    thread.start()

    return jsonify({'success': True, 'step': step_name})


@app.route('/api/stop')
def stop():
    """Stop running process"""
    global pipeline_state, running_processes

    for step_name, process in running_processes.items():
        try:
            process.terminate()
            process.wait(timeout=5)
        except (subprocess.TimeoutExpired, ProcessLookupError, OSError):
            process.kill()

    running_processes.clear()
    pipeline_state['status'] = 'idle'
    pipeline_state['current_step'] = None

    return jsonify({'success': True})


@app.route('/api/stream')
def stream():
    """Server-Sent Events stream for real-time updates"""
    def event_stream():
        last_log_count = 0
        while True:
            # Send status update
            status_data = {
                'type': 'status',
                'data': {
                    'status': pipeline_state['status'],
                    'current_step': pipeline_state['current_step'],
                    'start_time': pipeline_state['start_time'],
                    'end_time': pipeline_state['end_time'],
                    'error_message': pipeline_state['error_message']
                }
            }
            yield f"data: {json.dumps(status_data)}\n\n"

            # Send new logs
            current_log_count = len(pipeline_state['logs'])
            if current_log_count > last_log_count:
                new_logs = pipeline_state['logs'][last_log_count:]
                for log_entry in new_logs:
                    log_data = {
                        'type': 'log',
                        'data': log_entry
                    }
                    yield f"data: {json.dumps(log_data)}\n\n"
                last_log_count = current_log_count

            time.sleep(0.2)

    return Response(stream_with_context(event_stream()),
                   mimetype='text/event-stream',
                   headers={
                       'Cache-Control': 'no-cache',
                       'X-Accel-Buffering': 'no'
                   })


@app.route('/api/targets')
def get_targets():
    """Get saved target hypotheses"""
    try:
        from src.target_hypothesis import TargetHypothesisManager
        manager = TargetHypothesisManager(storage_dir=DATA_DIR / "targets")
        targets = manager.list_all()
        summary = manager.get_summary()

        return jsonify({
            'success': True,
            'summary': summary,
            'targets': [t.to_dict() for t in targets]
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/api/targets/export')
def export_targets():
    """Export targets for Phase 5"""
    try:
        from src.target_hypothesis import TargetHypothesisManager
        manager = TargetHypothesisManager(storage_dir=DATA_DIR / "targets")
        output_file = manager.export_for_phase5()

        with open(output_file, 'r') as f:
            data = json.load(f)

        return jsonify({
            'success': True,
            'file': str(output_file),
            'data': data
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/api/model', methods=['GET', 'POST'])
def model_api():
    """Get or set the current LLM model"""
    model_file = DATA_DIR / 'current_model.json'

    # Model to provider mapping
    model_providers = {
        'llama3.1:8b-fast-ctx2k': 'ollama',
        'llama3.1:8b-instruct-q4_0': 'ollama',
        'llama3.1:8b': 'ollama',
        'llama3.1:70b': 'ollama',
        'claude-sonnet-4-20250514': 'anthropic',
        'claude-opus-4-20250514': 'anthropic',
    }

    if request.method == 'GET':
        # Return current model
        try:
            if model_file.exists():
                with open(model_file, 'r') as f:
                    data = json.load(f)
                    return jsonify({
                        'success': True,
                        'model': data.get('model', 'llama3.1:70b'),
                        'provider': data.get('provider', 'ollama')
                    })
        except (json.JSONDecodeError, OSError):
            pass
        return jsonify({'success': True, 'model': 'llama3.1:70b', 'provider': 'ollama'})

    else:  # POST
        try:
            data = request.json
            model = data.get('model', 'llama3.1:70b')

            # Validate model
            valid_models = list(model_providers.keys())
            if model not in valid_models:
                return jsonify({'success': False, 'error': f'Invalid model. Must be one of: {valid_models}'}), 400

            # Determine provider
            provider = model_providers[model]

            # Check API key for cloud providers
            if provider == 'anthropic':
                api_key = os.getenv('ANTHROPIC_API_KEY')
                if not api_key:
                    return jsonify({
                        'success': False,
                        'error': 'ANTHROPIC_API_KEY not set. Add it to your .env file or environment.'
                    }), 400

            # Save to file for chat UI to read
            with open(model_file, 'w') as f:
                json.dump({
                    'model': model,
                    'provider': provider,
                    'updated_at': time.strftime('%Y-%m-%d %H:%M:%S')
                }, f)

            # Also update .env config
            config = load_config()
            config['LLM_MODEL'] = model
            config['LLM_PROVIDER'] = provider
            save_config(config)

            return jsonify({
                'success': True,
                'model': model,
                'provider': provider,
                'message': f'Model switched to {model} ({provider})'
            })

        except Exception as e:
            return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/api/llm-metrics')
def get_llm_metrics():
    """Get LLM performance metrics"""
    # Try to get metrics from Ollama
    ollama_stats = get_ollama_stats()

    # Try to load shared metrics from chat interface
    shared_metrics = load_shared_metrics()

    # Get KV cache usage from shared metrics or estimate from GPU
    kv_cache_usage = shared_metrics.get('kv_cache_usage')
    if kv_cache_usage is None:
        gpu_stats = get_gpu_utilization()
        if gpu_stats.get('available') and gpu_stats.get('devices'):
            gpu = gpu_stats['devices'][0]
            kv_cache_usage = gpu.get('memory_percent', 0) * 0.6

    return jsonify({
        'ttft': shared_metrics.get('ttft') or llm_metrics['ttft'],
        'tokens_per_sec': shared_metrics.get('tokens_per_sec') or llm_metrics['tokens_per_sec'],
        'cache_hits': shared_metrics.get('cache_hits', 0) or llm_metrics['cache_hits'],
        'cache_misses': shared_metrics.get('cache_misses', 0) or llm_metrics['cache_misses'],
        'total_queries': shared_metrics.get('total_queries', 0) or llm_metrics['total_queries'],
        'avg_latency': shared_metrics.get('avg_latency') or llm_metrics['avg_latency'],
        'kv_cache_usage': kv_cache_usage,
        'prefix_cache_hits': shared_metrics.get('prefix_cache_hits', 0),
        'ollama': ollama_stats,
    })


def load_shared_metrics():
    """Load shared metrics from the metrics file written by chat interface"""
    metrics_file = DATA_DIR / 'llm_metrics.json'
    try:
        if metrics_file.exists():
            with open(metrics_file, 'r') as f:
                return json.load(f)
    except (json.JSONDecodeError, OSError, IOError):
        pass
    return {}


def get_ollama_stats():
    """Get Ollama server statistics"""
    try:
        # Check Ollama health and get model info
        with urllib.request.urlopen('http://localhost:11434/api/tags', timeout=2) as response:
            data = json.loads(response.read().decode('utf-8'))
            models = data.get('models', [])

            # Find loaded models and their sizes
            model_info = []
            for model in models:
                model_info.append({
                    'name': model.get('name', 'unknown'),
                    'size': model.get('size', 0),
                    'modified': model.get('modified_at', '')
                })

            return {
                'available': True,
                'models_count': len(models),
                'models': model_info[:5],  # Top 5 models
            }
    except Exception as e:
        return {'available': False, 'error': str(e)}


def check_service_health(service_name, host, port):
    """Check if a service is healthy"""
    try:
        import socket
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(2)
        result = sock.connect_ex((host, port))
        sock.close()
        return result == 0
    except (OSError, socket.error, socket.timeout):
        return False


def check_milvus_collection():
    """Check Milvus collection status"""
    try:
        from pymilvus import connections, utility
        connections.connect('default', host='localhost', port=19530, timeout=5)
        collections = utility.list_collections()
        genomic_exists = 'genomic_evidence' in collections

        if genomic_exists:
            from pymilvus import Collection
            coll = Collection('genomic_evidence')
            stats = coll.num_entities

            connections.disconnect('default')
            return {
                'exists': True,
                'name': 'genomic_evidence',
                'num_entities': stats
            }
        else:
            connections.disconnect('default')
            return {'exists': False}
    except Exception as e:
        return {'exists': False, 'error': str(e)}


def check_docker():
    """Check if Docker is available"""
    try:
        result = subprocess.run(
            ['docker', '--version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return False


def check_gpu():
    """Check if GPU is available"""
    try:
        result = subprocess.run(
            ['nvidia-smi'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return False


def get_cpu_utilization():
    """Get CPU utilization percentage"""
    try:
        cpu_percent = psutil.cpu_percent(interval=0.1)
        cpu_per_core = psutil.cpu_percent(interval=0.1, percpu=True)
        cpu_freq = psutil.cpu_freq()
        cpu_count = psutil.cpu_count()

        return {
            'percent': round(cpu_percent, 1),
            'per_core': [round(p, 1) for p in cpu_per_core] if cpu_per_core else [],
            'frequency': round(cpu_freq.current, 0) if cpu_freq else 0,
            'count': cpu_count
        }
    except Exception as e:
        print(f"Error getting CPU utilization: {e}")
        return {'percent': 0, 'per_core': [], 'frequency': 0, 'count': 0}


def get_gpu_utilization():
    """Get GPU utilization percentage and memory usage"""
    if not NVML_AVAILABLE:
        return {'available': False, 'devices': []}

    try:
        pynvml.nvmlInit()
        device_count = pynvml.nvmlDeviceGetCount()

        devices = []
        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)

            name = pynvml.nvmlDeviceGetName(handle)
            if isinstance(name, bytes):
                name = name.decode('utf-8')

            util = pynvml.nvmlDeviceGetUtilizationRates(handle)
            gpu_util = util.gpu

            try:
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                mem_used_gb = mem_info.used / (1024**3)
                mem_total_gb = mem_info.total / (1024**3)
                mem_percent = (mem_info.used / mem_info.total) * 100
            except (pynvml.NVMLError, AttributeError, ZeroDivisionError):
                mem_used_gb = 0
                mem_total_gb = 0
                mem_percent = 0

            try:
                temp = pynvml.nvmlDeviceGetTemperature(handle, pynvml.NVML_TEMPERATURE_GPU)
            except (pynvml.NVMLError, AttributeError):
                temp = 0

            try:
                power = pynvml.nvmlDeviceGetPowerUsage(handle) / 1000.0
            except (pynvml.NVMLError, AttributeError):
                power = 0

            devices.append({
                'index': i,
                'name': name,
                'utilization': gpu_util,
                'memory_utilization': util.memory,
                'memory_used_gb': round(mem_used_gb, 2),
                'memory_total_gb': round(mem_total_gb, 2),
                'memory_percent': round(mem_percent, 1),
                'temperature': temp,
                'power_watts': round(power, 1)
            })

        pynvml.nvmlShutdown()

        return {'available': True, 'devices': devices}

    except Exception as e:
        print(f"Error getting GPU utilization: {e}")
        return {'available': False, 'devices': [], 'error': str(e)}


def get_memory_utilization():
    """Get system memory utilization"""
    try:
        mem = psutil.virtual_memory()
        swap = psutil.swap_memory()

        return {
            'percent': round(mem.percent, 1),
            'used_gb': round(mem.used / (1024**3), 2),
            'total_gb': round(mem.total / (1024**3), 2),
            'available_gb': round(mem.available / (1024**3), 2),
            'swap_percent': round(swap.percent, 1) if swap else 0
        }
    except Exception as e:
        print(f"Error getting memory utilization: {e}")
        return {
            'percent': 0,
            'used_gb': 0,
            'total_gb': 0,
            'available_gb': 0,
            'swap_percent': 0
        }


if __name__ == '__main__':
    # Create required directories
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    # Run server
    print("=" * 60)
    print("RAG Chat Pipeline Web Portal")
    print("=" * 60)
    print(f"Project Root: {PROJECT_ROOT}")
    print(f"Opening portal at: http://localhost:5001")
    print("=" * 60)

    app.run(host='0.0.0.0', port=5001, debug=True, threaded=True)
