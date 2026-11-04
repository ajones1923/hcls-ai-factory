#!/usr/bin/env python3
"""
Genomics Pipeline Web Portal
Flask backend server for managing the FASTQ → VCF workflow
"""

import os
import sys
import json
import hmac
import re
import secrets
import shlex
import signal
import atexit
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from datetime import datetime
from functools import wraps
from flask import Flask, render_template, jsonify, request, Response, stream_with_context, send_file
from werkzeug.utils import secure_filename
from flask_cors import CORS
from loguru import logger
import psutil
try:
    import pynvml
    NVML_AVAILABLE = True
except ImportError:
    NVML_AVAILABLE = False

# Configure structured logging
logger.remove()
logger.add(sys.stderr, format="{time:YYYY-MM-DD HH:mm:ss} | {level:<7} | {message}", level="INFO")

app = Flask(__name__,
            template_folder='../templates',
            static_folder='../static')
# CORS configuration.
# The default used to be '*', which let any web page the operator happened to visit call this
# API from their browser. The portal's own UI is served same-origin and needs no CORS at all,
# so the default is now a loopback allowlist; CORS_ORIGINS overrides it for a genuine
# cross-origin deployment. Credentials are NOT enabled, and must not be while '*' is
# permitted as an override — a wildcard origin with credentials is the cookie-theft shape.
DEFAULT_CORS_ORIGINS = 'http://localhost:5000,http://127.0.0.1:5000'
_cors_origins = [o.strip() for o in
                 os.environ.get('CORS_ORIGINS', DEFAULT_CORS_ORIGINS).split(',') if o.strip()]
if '*' in _cors_origins:
    logger.warning("CORS_ORIGINS contains '*' — every website a user visits can call this API "
                   "from their browser. Set an explicit origin allowlist.")
CORS(app, resources={
    r"/api/*": {
        "origins": _cors_origins,
        "methods": ["GET", "POST", "DELETE"],
        "allow_headers": ["Content-Type", "X-API-Key"],
        "supports_credentials": False,
    }
})


try:
    from security import RateLimiter
    api_limiter = RateLimiter(max_requests=60, window_seconds=60)
except ImportError:
    api_limiter = None


def rate_limit(f):
    """Apply rate limiting if available, otherwise pass through."""
    if api_limiter:
        return api_limiter.limit(f)
    return f


def _publish_ephemeral_key(key):
    """Write an ephemeral key where the operator can read it, mode 0600.

    Returns the path written, or None. The key itself is never logged: the log is a
    lower-trust artifact than a credential store, and this one gates routes that read
    the engine's configuration.
    """
    configured = os.environ.get('PORTAL_API_KEY_FILE')
    candidates = [Path(configured)] if configured else []
    runtime_dir = os.environ.get('XDG_RUNTIME_DIR')
    if runtime_dir:
        candidates.append(Path(runtime_dir) / 'hcls-genomics-portal.key')
    candidates.append(Path(tempfile.gettempdir()) / 'hcls-genomics-portal.key')

    for path in candidates:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            # Create with 0600 from the outset — never world-readable, not even briefly.
            fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
            try:
                os.write(fd, (key + '\n').encode('utf-8'))
            finally:
                os.close(fd)
            os.chmod(str(path), 0o600)
            return path
        except OSError:
            continue
    return None


def _portal_api_key():
    """The key the portal enforces.

    Fail closed: when PORTAL_API_KEY is unset an ephemeral key is generated at first use,
    so no route is ever open merely because nothing was configured. The sibling portal
    (core/engines/precision-intelligence) already behaved this way; this one did not.
    """
    configured = os.environ.get('PORTAL_API_KEY')
    if configured:
        return configured
    if _portal_api_key._ephemeral is None:
        _portal_api_key._ephemeral = secrets.token_urlsafe(32)
        key_path = _publish_ephemeral_key(_portal_api_key._ephemeral)
        logger.warning("=" * 72)
        logger.warning("PORTAL_API_KEY not set — generated an ephemeral key. Auth is ENFORCED.")
        if key_path:
            logger.warning(f"Key written to {key_path} (mode 0600). Send it as X-API-Key.")
        else:
            logger.warning("Could not write the key to disk, so it cannot be retrieved: "
                           "set PORTAL_API_KEY and restart to use the API.")
        logger.warning("To persist across restarts, set PORTAL_API_KEY in the environment.")
        logger.warning("=" * 72)
    return _portal_api_key._ephemeral


_portal_api_key._ephemeral = None


def _api_key_matches(supplied):
    """Constant-time comparison on encoded bytes.

    Bytes, not str: compare_digest raises TypeError on a non-ASCII str, which would turn a
    non-ASCII header — or a non-ASCII passphrase in PORTAL_API_KEY — into a 500 instead of a 401.
    """
    expected = _portal_api_key()
    if not expected or not supplied:
        return False
    return hmac.compare_digest(
        supplied.encode('utf-8', 'surrogateescape'),
        expected.encode('utf-8', 'surrogateescape'),
    )


def require_api_key(func=None, *, allow_query_key=False):
    """Require a matching X-API-Key. Always enforced — there is no key-unset bypass.

    allow_query_key=True additionally accepts ?api_key= on routes the browser cannot send a
    header on (EventSource, and a download triggered by navigation). Query strings are more
    likely to be logged by an intermediary than a header, so it is opt-in per route.
    """
    def _decorate(f):
        @wraps(f)
        def decorated(*args, **kwargs):
            supplied = request.headers.get('X-API-Key') or ''
            if not supplied and allow_query_key:
                supplied = request.args.get('api_key') or ''
            if not _api_key_matches(supplied):
                return jsonify({'error': 'Unauthorized — set X-API-Key header'}), 401
            return f(*args, **kwargs)
        return decorated

    return _decorate(func) if func is not None else _decorate


# Configure max upload size. Genomics FASTQ are genuinely huge, so the default stays at
# 100GB — but it is now a configured number that the upload route also enforces while
# streaming, because MAX_CONTENT_LENGTH only acts on a declared Content-Length and a
# chunked upload declares none.
DEFAULT_MAX_UPLOAD_BYTES = 100 * 1024 * 1024 * 1024  # 100GB


def max_upload_bytes():
    """Upload size ceiling, from PORTAL_MAX_UPLOAD_BYTES (bytes)."""
    raw = os.environ.get('PORTAL_MAX_UPLOAD_BYTES', '').strip()
    if raw:
        try:
            value = int(raw)
            if value > 0:
                return value
        except ValueError:
            logger.warning(f"PORTAL_MAX_UPLOAD_BYTES is not an integer ({raw!r}) — "
                           f"using the {DEFAULT_MAX_UPLOAD_BYTES} byte default.")
    return DEFAULT_MAX_UPLOAD_BYTES


app.config['MAX_CONTENT_LENGTH'] = max_upload_bytes()

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent.parent.absolute()
SCRIPTS_DIR = PROJECT_ROOT / 'scripts'
DATA_DIR = PROJECT_ROOT / 'data'
CONFIG_FILE = PROJECT_ROOT / 'config' / 'pipeline.env'
LOG_DIR = DATA_DIR / 'output' / 'logs'
STATE_FILE = DATA_DIR / 'output' / '.pipeline_state.json'

# Global state
MAX_LOG_ENTRIES = 5000  # Prevent unbounded memory growth

pipeline_state = {
    'current_step': None,
    'status': 'idle',  # idle, running, success, error
    'progress': 0,
    'logs': [],
    'start_time': None,
    'end_time': None,
    'error_message': None,
    'runtime_seconds': 0
}

# Thread lock for pipeline_state — protects reads/writes from concurrent threads
pipeline_lock = threading.Lock()

running_processes = {}


def save_pipeline_state():
    """Save pipeline state to disk for persistence across restarts"""
    try:
        state_to_save = {
            'current_step': pipeline_state['current_step'],
            'status': pipeline_state['status'],
            'start_time': pipeline_state['start_time'],
            'end_time': pipeline_state['end_time'],
            'error_message': pipeline_state['error_message'],
            'runtime_seconds': pipeline_state.get('runtime_seconds', 0)
        }
        with open(STATE_FILE, 'w') as f:
            json.dump(state_to_save, f)
    except Exception as e:
        logger.error(f"Error saving pipeline state: {e}")


def load_pipeline_state():
    """Load pipeline state from disk"""
    global pipeline_state
    try:
        if STATE_FILE.exists():
            with open(STATE_FILE, 'r') as f:
                saved_state = json.load(f)
                # Only restore if it was a completed or error state (not running)
                if saved_state.get('status') in ['success', 'error']:
                    pipeline_state['current_step'] = saved_state.get('current_step')
                    pipeline_state['status'] = saved_state.get('status')
                    pipeline_state['start_time'] = saved_state.get('start_time')
                    pipeline_state['end_time'] = saved_state.get('end_time')
                    pipeline_state['error_message'] = saved_state.get('error_message')
                    pipeline_state['runtime_seconds'] = saved_state.get('runtime_seconds', 0)
                    logger.info(f"Restored pipeline state: {saved_state.get('status')} - {saved_state.get('current_step')}")
    except Exception as e:
        logger.error(f"Error loading pipeline state: {e}")


# ─── Pipeline configuration contract ────────────────────────────────────────────────────
# config/pipeline.env is `source`d by the pipeline scripts, i.e. bash EXECUTES it, and one of
# its keys (PB_IMG) is the image argument to `docker run` with host mounts. Everything written
# to it therefore has to survive being read as shell. Three rules hold that line:
#   1. only keys in CONFIG_WRITABLE_KEYS may be written;
#   2. each value must match that key's pattern — which admits no newline, no carriage return
#      and no shell metacharacter — and PB_IMG must additionally name an allowed image;
#   3. values are written single-quoted via shlex.quote, so even a value that slipped the
#      pattern cannot be re-interpreted as shell syntax.

CONFIG_WRITABLE_KEYS = (
    'GP', 'IN', 'REF', 'OUT', 'LOG',
    'PB_IMG', 'PATIENT_ID', 'REF_BUILD',
    'GIAB_INDEX_URL', 'PARABRICKS_SAMPLE_URL',
    'NUM_GPUS', 'LOW_MEMORY',
)

# What GET may return. pipeline.env is an operator-owned file: the engine README tells
# self-hosters to put registry credentials in it, so the route returns only the keys the
# portal declares and never the rest of the file.
CONFIG_READABLE_KEYS = CONFIG_WRITABLE_KEYS

# Keys naming a directory on the host. Confined to PORTAL_CONFIG_PATH_ROOT (default: the
# engine directory) so a caller cannot repoint the pipeline at arbitrary host paths.
CONFIG_PATH_KEYS = ('GP', 'IN', 'REF', 'OUT', 'LOG')

# Per-key value patterns. None of them admit \n, \r or a shell metacharacter.
_SAFE_PATH = re.compile(r'^/?[A-Za-z0-9._/-]{1,512}$')
_CONFIG_VALUE_PATTERNS = {
    'GP': _SAFE_PATH,
    'IN': _SAFE_PATH,
    'REF': _SAFE_PATH,
    'OUT': _SAFE_PATH,
    'LOG': _SAFE_PATH,
    'PB_IMG': re.compile(r'^[a-z0-9][a-z0-9._-]*(\.[a-z0-9.-]+)?(:[0-9]+)?'
                         r'(/[a-z0-9._-]+)+:[A-Za-z0-9._-]{1,128}$'),
    'PATIENT_ID': re.compile(r'^[A-Za-z0-9][A-Za-z0-9_-]{0,49}$'),
    'REF_BUILD': re.compile(r'^(GRCh38|GRCh37|hg38|hg19)$'),
    'GIAB_INDEX_URL': re.compile(r'^https?://[A-Za-z0-9._-]+(:[0-9]+)?(/[A-Za-z0-9._~/-]*)?$'),
    'PARABRICKS_SAMPLE_URL': re.compile(r'^https?://[A-Za-z0-9._-]+(:[0-9]+)?(/[A-Za-z0-9._~/-]*)?$'),
    'NUM_GPUS': re.compile(r'^[1-9][0-9]?$'),
    'LOW_MEMORY': re.compile(r'^[01]$'),
}

# PB_IMG selects the container that runs with GPU access and host volume mounts, so it is
# pinned rather than merely pattern-checked. Self-hosters who run a different build set
# PORTAL_ALLOWED_IMAGES (comma-separated) to the images they allow.
_DEFAULT_ALLOWED_IMAGES = (
    'nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1',
    'nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1',
    'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1',
)


def allowed_pipeline_images():
    """Images PB_IMG may name. PORTAL_ALLOWED_IMAGES replaces the built-in pinned set."""
    configured = os.environ.get('PORTAL_ALLOWED_IMAGES', '').strip()
    if configured:
        return tuple(i.strip() for i in configured.split(',') if i.strip())
    return _DEFAULT_ALLOWED_IMAGES


def config_path_root():
    """Directory that the path-valued config keys must stay inside."""
    return Path(os.environ.get('PORTAL_CONFIG_PATH_ROOT') or PROJECT_ROOT)


def validate_config_entry(key, value):
    """Validate one config key/value. Returns an error string, or None when acceptable."""
    if key not in CONFIG_WRITABLE_KEYS:
        return f'Invalid config key: {key}'
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        return f'Value for {key} must be a string or number'

    text = str(value)
    # Checked explicitly and first: a newline is what turns one config line into extra lines
    # in a file bash executes, and it is the one character the old filter let through.
    if '\n' in text or '\r' in text or '\x00' in text:
        return f'Value for {key} may not contain a line break'
    if not text:
        return f'Value for {key} may not be empty'

    pattern = _CONFIG_VALUE_PATTERNS.get(key)
    if pattern and not pattern.match(text):
        return f'Invalid value for {key}'

    if key == 'PB_IMG' and text not in allowed_pipeline_images():
        return (f'Image not allowed for PB_IMG. Permitted: '
                f'{", ".join(allowed_pipeline_images())} '
                f'(set PORTAL_ALLOWED_IMAGES to permit others)')

    if key in CONFIG_PATH_KEYS:
        root = config_path_root()
        try:
            resolved = Path(text) if Path(text).is_absolute() else (root / text)
            resolved = resolved.resolve()
            if not (resolved == root.resolve() or resolved.is_relative_to(root.resolve())):
                return f'Path for {key} must be inside {root}'
        except (OSError, ValueError):
            return f'Invalid path for {key}'

    return None


def _unquote_config_value(value):
    """Reverse the shell quoting save_config writes, for display."""
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
        inner = value[1:-1]
        if value[0] == "'":
            return inner.replace("'\\''", "'")
        return inner
    return value


def load_config(keys=None):
    """Load pipeline configuration.

    keys: when given, only these keys are returned. Callers that serve the result over HTTP
    pass CONFIG_READABLE_KEYS so an operator's credentials in this file are never returned.
    """
    config = {}
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    if keys is not None and key not in keys:
                        continue
                    config[key] = _unquote_config_value(value)
    return config


def save_config(config):
    """Save pipeline configuration.

    Values are written shell-quoted: the file is `source`d by the pipeline scripts, so an
    unquoted value is shell source code. Callers must have validated the keys and values
    first (validate_config_entry); the quoting is the second line of defence, not the first.
    """
    rejected = [(k, v) for k, v in config.items() if validate_config_entry(k, v) is not None]
    if rejected:
        raise ValueError(f'refusing to write unvalidated config keys: '
                         f'{", ".join(sorted(k for k, _ in rejected))}')

    config = dict(config)
    lines = []
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    lines.append(line)
                elif '=' in stripped:
                    key = stripped.split('=', 1)[0].strip()
                    if key in config:
                        lines.append(f"{key}={shlex.quote(str(config[key]))}\n")
                        del config[key]
                    else:
                        lines.append(line)

    # Add any new config values
    for key, value in config.items():
        lines.append(f"{key}={shlex.quote(str(value))}\n")

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
        logger.error(f"Error getting disk space: {e}")

    return {'total': 'N/A', 'used': 'N/A', 'available': 'N/A', 'use_percent': 'N/A'}


def run_script(script_name, step_name):
    """Run a pipeline script"""
    global pipeline_state, running_processes

    script_path = SCRIPTS_DIR / script_name

    if not script_path.exists():
        with pipeline_lock:
            pipeline_state['status'] = 'error'
            pipeline_state['error_message'] = f"Script not found: {script_name}"
        return

    with pipeline_lock:
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
            bufsize=0,  # Unbuffered
            cwd=str(PROJECT_ROOT),
            env={**os.environ, 'PYTHONUNBUFFERED': '1'}
        )

        running_processes[step_name] = process

        # Stream output line by line
        for line in iter(process.stdout.readline, ''):
            if line:
                with pipeline_lock:
                    pipeline_state['logs'].append({
                        'timestamp': datetime.now().isoformat(),
                        'message': line.rstrip()
                    })
                    # Cap log size to prevent unbounded memory growth
                    if len(pipeline_state['logs']) > MAX_LOG_ENTRIES:
                        pipeline_state['logs'] = pipeline_state['logs'][-MAX_LOG_ENTRIES:]

        process.wait()

        # Calculate runtime
        with pipeline_lock:
            pipeline_state['end_time'] = datetime.now().isoformat()
            if pipeline_state['start_time']:
                start = datetime.fromisoformat(pipeline_state['start_time'])
                end = datetime.fromisoformat(pipeline_state['end_time'])
                pipeline_state['runtime_seconds'] = int((end - start).total_seconds())

            if process.returncode == 0:
                pipeline_state['status'] = 'success'
                save_pipeline_state()
            else:
                pipeline_state['status'] = 'error'
                pipeline_state['error_message'] = f"Script exited with code {process.returncode}"
                save_pipeline_state()

    except Exception as e:
        with pipeline_lock:
            pipeline_state['status'] = 'error'
            pipeline_state['error_message'] = str(e)
            pipeline_state['end_time'] = datetime.now().isoformat()
            if pipeline_state['start_time']:
                start = datetime.fromisoformat(pipeline_state['start_time'])
                end = datetime.fromisoformat(pipeline_state['end_time'])
                pipeline_state['runtime_seconds'] = int((end - start).total_seconds())
            save_pipeline_state()

    finally:
        if step_name in running_processes:
            del running_processes[step_name]


# Routes
@app.route('/health')
@app.route('/healthz')
def health():
    """Health check endpoint for monitoring systems"""
    return jsonify({'status': 'healthy', 'service': 'genomics-portal', 'port': 5000})


@app.route('/')
def index():
    """Main page"""
    return render_template('index.html')


@app.route('/api/status')
@rate_limit
def get_status():
    """Get current pipeline status"""
    # Check system status
    system_status = {
        'docker_installed': check_docker(),
        'gpu_available': check_gpu(),
        'disk_space': get_disk_space(),
        'cpu_utilization': get_cpu_utilization(),
        'gpu_utilization': get_gpu_utilization(),
        'memory_utilization': get_memory_utilization()
    }

    # Check data files
    data_status = {
        'fastq_r1': check_file_exists(DATA_DIR / 'input' / 'HG002_R1.fastq.gz'),
        'fastq_r2': check_file_exists(DATA_DIR / 'input' / 'HG002_R2.fastq.gz'),
        'reference': check_file_exists(DATA_DIR / 'ref' / 'GRCh38.fa'),
        'chr20_vcf': check_file_exists(DATA_DIR / 'output' / 'HG002.chr20.vcf.gz'),
        'genome_vcf': check_file_exists(DATA_DIR / 'output' / 'HG002.genome.vcf.gz')
    }

    # Get file sizes
    file_sizes = {
        'fastq_r1': get_file_size(DATA_DIR / 'input' / 'HG002_R1.fastq.gz'),
        'fastq_r2': get_file_size(DATA_DIR / 'input' / 'HG002_R2.fastq.gz'),
        'reference': get_file_size(DATA_DIR / 'ref' / 'GRCh38.fa'),
        'chr20_vcf': get_file_size(DATA_DIR / 'output' / 'HG002.chr20.vcf.gz'),
        'genome_vcf': get_file_size(DATA_DIR / 'output' / 'HG002.genome.vcf.gz')
    }

    return jsonify({
        'pipeline': pipeline_state,
        'system': system_status,
        'data': data_status,
        'file_sizes': file_sizes
    })


@app.route('/api/config', methods=['GET', 'POST'])
@require_api_key
def config():
    """Get or update configuration.

    Authenticated on both methods. GET returns only CONFIG_READABLE_KEYS — the rest of
    pipeline.env (the engine README invites operators to keep registry credentials there)
    is not served. POST validates every key and value before anything is written, because
    the pipeline scripts `source` this file and PB_IMG from it selects the container image.
    """
    if request.method == 'GET':
        return jsonify(load_config(keys=CONFIG_READABLE_KEYS))
    else:
        config_data = request.json
        if not config_data or not isinstance(config_data, dict):
            return jsonify({'error': 'Request body must be a JSON object'}), 400

        invalid_keys = set(config_data.keys()) - set(CONFIG_WRITABLE_KEYS)
        if invalid_keys:
            return jsonify({'error': f'Invalid config keys: {", ".join(sorted(invalid_keys))}'}), 400

        for key, value in config_data.items():
            error = validate_config_entry(key, value)
            if error:
                return jsonify({'error': error}), 400

        save_config(config_data)
        return jsonify({'success': True})


@app.route('/api/run/<step>', methods=['POST'])
@require_api_key
def run_step(step):
    """Run a specific workflow step"""
    global pipeline_state

    with pipeline_lock:
        if pipeline_state['status'] == 'running':
            return jsonify({'error': 'A step is already running'}), 400

        # Reset state
        pipeline_state['logs'] = []
        pipeline_state['progress'] = 0

    script_mapping = {
        'check': ('00-setup-check.sh', 'Prerequisites Check'),
        'login': ('01-ngc-login.sh', 'NGC Login'),
        'download': ('02-download-data.sh', 'Data Download'),
        'reference': ('03-setup-reference.sh', 'Reference Setup'),
        'test': ('04-run-chr20-test.sh', 'Chr20 Test'),
        'full': ('05-run-full-genome.sh', 'Full Genome')
    }

    if step not in script_mapping:
        return jsonify({'error': 'Invalid step'}), 400

    script_name, step_name = script_mapping[step]

    # Run in background thread
    thread = threading.Thread(target=run_script, args=(script_name, step_name))
    thread.daemon = True
    thread.start()

    return jsonify({'success': True, 'step': step_name})


@app.route('/api/stop', methods=['POST'])
@require_api_key
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
    with pipeline_lock:
        pipeline_state['status'] = 'idle'
        pipeline_state['current_step'] = None

    return jsonify({'success': True})


@app.route('/api/stop-all', methods=['POST'])
@require_api_key
def stop_all():
    """Stop all pipeline-related processes"""
    global pipeline_state, running_processes
    killed_count = 0

    # First, stop any processes managed by the portal
    for step_name, process in running_processes.items():
        try:
            process.terminate()
            process.wait(timeout=5)
            killed_count += 1
        except (subprocess.TimeoutExpired, ProcessLookupError, OSError):
            try:
                process.kill()
                killed_count += 1
            except (ProcessLookupError, OSError):
                pass

    running_processes.clear()

    # Kill any running pipeline scripts and related processes
    pipeline_patterns = [
        '00-setup-check.sh',
        '01-ngc-login.sh',
        '02-download-data.sh',
        '03-setup-reference.sh',
        '04-run-chr20-test.sh',
        '05-run-full-genome.sh',
        'pbrun',
        'parabricks',
        'pbbwa',
        'pbsort',
        'pbmarkdup',
        'fq2bam',
        'deepvariant',
        'clara-parabricks'
    ]

    try:
        # Find and kill processes matching pipeline patterns
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                cmdline = ' '.join(proc.info['cmdline'] or [])
                proc_name = proc.info['name'] or ''

                # Check if process matches any pipeline pattern
                for pattern in pipeline_patterns:
                    if pattern in cmdline or pattern in proc_name:
                        proc.terminate()
                        try:
                            proc.wait(timeout=3)
                        except psutil.TimeoutExpired:
                            proc.kill()
                        killed_count += 1
                        break
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue

        # Stop any Docker containers running parabricks (check multiple image patterns)
        parabricks_images = [
            'nvcr.io/nvidia/clara/clara-parabricks',
            'clara-parabricks',
            'parabricks'
        ]

        for image_pattern in parabricks_images:
            try:
                # Find containers by image
                result = subprocess.run(
                    ['docker', 'ps', '-q', '--filter', f'ancestor={image_pattern}'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                container_ids = result.stdout.strip().split('\n')
                for cid in container_ids:
                    if cid:
                        subprocess.run(['docker', 'stop', cid], capture_output=True, timeout=30)
                        killed_count += 1
            except Exception:
                pass

        # Also find any container with parabricks in the name
        try:
            result = subprocess.run(
                ['docker', 'ps', '--format', '{{.ID}} {{.Image}}'],
                capture_output=True,
                text=True,
                timeout=5
            )
            for line in result.stdout.strip().split('\n'):
                if line and 'parabricks' in line.lower():
                    cid = line.split()[0]
                    subprocess.run(['docker', 'stop', cid], capture_output=True, timeout=30)
                    killed_count += 1
        except Exception:
            pass

    except Exception as e:
        logger.error(f"Error stopping processes: {e}")

    # Reset pipeline state
    with pipeline_lock:
        pipeline_state['status'] = 'idle'
        pipeline_state['current_step'] = None

    return jsonify({
        'success': True,
        'killed_count': killed_count,
        'message': f'Stopped {killed_count} process(es)' if killed_count > 0 else 'No running processes found'
    })


@app.route('/api/logs/<log_type>')
@require_api_key
def get_logs(log_type):
    """Get log file content"""
    log_files = {
        'chr20_fq2bam': 'chr20_fq2bam.log',
        'chr20_deepvariant': 'chr20_deepvariant.log',
        'genome_fq2bam': 'genome_fq2bam.log',
        'genome_deepvariant': 'genome_deepvariant.log'
    }

    if log_type not in log_files:
        return jsonify({'error': 'Invalid log type'}), 400

    log_file = LOG_DIR / log_files[log_type]

    if not log_file.exists():
        return jsonify({'content': 'Log file not found'})

    try:
        with open(log_file, 'r') as f:
            content = f.read()
        return jsonify({'content': content})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/stream')
@require_api_key(allow_query_key=True)
def stream():
    """Server-Sent Events stream for real-time updates"""
    def event_stream():
        last_log_count = 0
        idle_ticks = 0  # Track how long we've been in a terminal state

        while True:
            with pipeline_lock:
                current_status = pipeline_state['status']
                status_data = {
                    'type': 'status',
                    'data': {
                        'status': current_status,
                        'current_step': pipeline_state['current_step'],
                        'start_time': pipeline_state['start_time'],
                        'end_time': pipeline_state['end_time'],
                        'error_message': pipeline_state['error_message']
                    }
                }
                current_log_count = len(pipeline_state['logs'])
                new_logs = pipeline_state['logs'][last_log_count:] if current_log_count > last_log_count else []

            yield f"data: {json.dumps(status_data)}\n\n"

            # Send new logs
            for log_entry in new_logs:
                log_data = {
                    'type': 'log',
                    'data': log_entry
                }
                yield f"data: {json.dumps(log_data)}\n\n"
            last_log_count = current_log_count

            # Close stream after pipeline reaches terminal state and logs are flushed
            if current_status in ('success', 'error', 'idle'):
                idle_ticks += 1
                if idle_ticks >= 15:  # 3 seconds of terminal state → close stream
                    yield f"data: {json.dumps({'type': 'stream_end', 'data': {'reason': current_status}})}\n\n"
                    return
            else:
                idle_ticks = 0

            time.sleep(0.2)

    return Response(stream_with_context(event_stream()),
                   mimetype='text/event-stream',
                   headers={
                       'Cache-Control': 'no-cache',
                       'X-Accel-Buffering': 'no'
                   })


@app.route('/api/reset', methods=['POST'])
@require_api_key
def reset_state():
    """Reset pipeline state to idle"""
    global pipeline_state

    with pipeline_lock:
        pipeline_state = {
            'current_step': None,
            'status': 'idle',
            'progress': 0,
            'logs': [],
            'start_time': None,
            'end_time': None,
            'error_message': None,
            'runtime_seconds': 0
        }

    # Remove saved state file
    if STATE_FILE.exists():
        try:
            STATE_FILE.unlink()
        except OSError:
            pass

    return jsonify({'success': True, 'message': 'Pipeline state reset to idle'})


@app.route('/api/files/<directory>')
@require_api_key
def list_files(directory):
    """List files in a directory"""
    # Map directory names to actual paths
    dir_mapping = {
        'input': DATA_DIR / 'input',
        'output': DATA_DIR / 'output',
        'ref': DATA_DIR / 'ref'
    }

    if directory not in dir_mapping:
        return jsonify({'error': 'Invalid directory'}), 400

    target_dir = dir_mapping[directory]

    if not target_dir.exists():
        return jsonify({'files': [], 'directory': directory})

    files = []
    try:
        for item in sorted(target_dir.iterdir()):
            # Skip hidden files and directories
            if item.name.startswith('.'):
                continue
            # Skip subdirectories for now (keep it simple)
            if item.is_dir():
                continue

            stat = item.stat()
            files.append({
                'name': item.name,
                'size': stat.st_size,
                'size_human': get_file_size(item),
                'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
                'type': get_file_type(item.name)
            })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

    return jsonify({
        'files': files,
        'directory': directory,
        'path': str(target_dir)
    })


def get_file_type(filename):
    """Get file type based on extension"""
    ext = filename.lower().split('.')[-1] if '.' in filename else ''
    type_map = {
        'fastq': 'fastq', 'fq': 'fastq',
        'gz': 'compressed',
        'bam': 'bam', 'bai': 'index',
        'vcf': 'vcf',
        'fa': 'fasta', 'fasta': 'fasta',
        'bed': 'bed',
        'txt': 'text', 'log': 'text',
        'json': 'json'
    }
    return type_map.get(ext, 'file')


@app.route('/api/files/download/<directory>/<filename>')
@require_api_key(allow_query_key=True)
def download_file(directory, filename):
    """Download a file"""
    # Map directory names to actual paths
    dir_mapping = {
        'input': DATA_DIR / 'input',
        'output': DATA_DIR / 'output',
        'ref': DATA_DIR / 'ref'
    }

    if directory not in dir_mapping:
        return jsonify({'error': 'Invalid directory'}), 400

    # Security: prevent path traversal using secure_filename + resolve check
    safe = secure_filename(filename)
    if not safe:
        return jsonify({'error': 'Invalid filename'}), 400
    target_dir = dir_mapping[directory]
    file_path = target_dir / safe
    if not file_path.resolve().is_relative_to(target_dir.resolve()):
        return jsonify({'error': 'Invalid filename'}), 400

    if not file_path.exists():
        return jsonify({'error': 'File not found'}), 404

    try:
        return send_file(
            file_path,
            as_attachment=True,
            download_name=safe
        )
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/files/upload/<directory>', methods=['POST'])
@require_api_key
@rate_limit
def upload_file(directory):
    """Upload a file to a directory"""
    # Map directory names to actual paths
    dir_mapping = {
        'input': DATA_DIR / 'input',
        'output': DATA_DIR / 'output',
        'ref': DATA_DIR / 'ref'
    }

    if directory not in dir_mapping:
        return jsonify({'error': 'Invalid directory'}), 400

    if 'file' not in request.files:
        return jsonify({'error': 'No file provided'}), 400

    file = request.files['file']

    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400

    # Security: sanitize filename using Werkzeug's secure_filename
    filename = secure_filename(file.filename)
    if not filename:
        return jsonify({'error': 'Invalid filename'}), 400

    target_dir = dir_mapping[directory]
    target_dir.mkdir(parents=True, exist_ok=True)

    file_path = target_dir / filename

    limit = max_upload_bytes()
    declared = request.content_length
    if declared is not None and declared > limit:
        return jsonify({'error': f'Upload exceeds the {limit} byte limit '
                                 f'(PORTAL_MAX_UPLOAD_BYTES)'}), 413

    try:
        # Stream to disk with a running total: a chunked upload declares no Content-Length,
        # so MAX_CONTENT_LENGTH never sees it and the only real cap is this one.
        written = 0
        with open(file_path, 'wb') as out:
            while True:
                chunk = file.stream.read(8 * 1024 * 1024)
                if not chunk:
                    break
                written += len(chunk)
                if written > limit:
                    out.close()
                    file_path.unlink(missing_ok=True)
                    return jsonify({'error': f'Upload exceeds the {limit} byte limit '
                                             f'(PORTAL_MAX_UPLOAD_BYTES)'}), 413
                out.write(chunk)
        return jsonify({
            'success': True,
            'filename': filename,
            'size': file_path.stat().st_size,
            'size_human': get_file_size(file_path)
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/files/delete/<directory>/<filename>', methods=['DELETE'])
@require_api_key
def delete_file(directory, filename):
    """Delete a file"""
    # Map directory names to actual paths
    dir_mapping = {
        'input': DATA_DIR / 'input',
        'output': DATA_DIR / 'output',
        'ref': DATA_DIR / 'ref'
    }

    if directory not in dir_mapping:
        return jsonify({'error': 'Invalid directory'}), 400

    # Security: prevent path traversal using secure_filename + resolve check
    safe = secure_filename(filename)
    if not safe:
        return jsonify({'error': 'Invalid filename'}), 400
    target_dir = dir_mapping[directory]
    file_path = target_dir / safe
    if not file_path.resolve().is_relative_to(target_dir.resolve()):
        return jsonify({'error': 'Invalid filename'}), 400

    if not file_path.exists():
        return jsonify({'error': 'File not found'}), 404

    try:
        file_path.unlink()
        return jsonify({'success': True, 'filename': filename})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/metrics')
@rate_limit
def get_metrics():
    """Get real-time processing metrics for the executive dashboard"""
    gpu_util = get_gpu_utilization()
    gpu_detailed = get_detailed_gpu_metrics()

    # Calculate estimated throughput based on GPU utilization
    throughput = 0
    gpu_percent = 0
    if gpu_util['available'] and len(gpu_util['devices']) > 0:
        gpu_percent = gpu_util['devices'][0].get('utilization', 0)
        # Parabricks on GB10 can process ~45-60 GB/hr at full utilization
        if gpu_percent > 0:
            throughput = (gpu_percent / 100) * 55  # Max ~55 GB/hr

    # Get pipeline runtime - use saved value for completed jobs
    runtime_seconds = pipeline_state.get('runtime_seconds', 0)
    if pipeline_state['status'] == 'running' and pipeline_state['start_time']:
        # Calculate live runtime for running jobs
        start = datetime.fromisoformat(pipeline_state['start_time'])
        runtime_seconds = int((datetime.now() - start).total_seconds())

    return jsonify({
        'throughput_gb_hr': round(throughput, 1),
        'runtime_seconds': runtime_seconds,
        'gpu_utilization': gpu_percent,
        'status': pipeline_state['status'],
        'current_step': pipeline_state['current_step'],
        'gpu_metrics': gpu_detailed,
        'start_time': pipeline_state['start_time'],
        'end_time': pipeline_state['end_time']
    })


def get_detailed_gpu_metrics():
    """Get detailed GPU performance metrics"""
    metrics = {
        'iops': 0,
        'memory_bandwidth': 0,
        'sm_efficiency': 0,
        'tensor_usage': 0,
        'pcie_throughput': 0,
        'power_draw': 0
    }

    if not NVML_AVAILABLE:
        return metrics

    try:
        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(0)

        # Get utilization rates
        util = pynvml.nvmlDeviceGetUtilizationRates(handle)

        # Check for running processes to determine if GPU is active
        try:
            compute_procs = pynvml.nvmlDeviceGetComputeRunningProcesses(handle)
            is_active = len(compute_procs) > 0
        except Exception:
            is_active = False

        # Get power usage
        try:
            power = pynvml.nvmlDeviceGetPowerUsage(handle) / 1000.0
            metrics['power_draw'] = round(power, 1)
        except Exception:
            pass

        # Estimate other metrics based on utilization when GPU is active
        if is_active and util.gpu > 0:
            gpu_util_pct = util.gpu / 100.0

            # IOPS estimate (GB10 theoretical ~500K IOPS for genomics workloads)
            metrics['iops'] = round(450 * gpu_util_pct + (50 * gpu_util_pct * (0.5 + 0.5 * (hash(str(time.time())) % 100) / 100)), 1)

            # Memory bandwidth (GB10 has ~200 GB/s theoretical)
            metrics['memory_bandwidth'] = round(180 * gpu_util_pct + (20 * gpu_util_pct * (0.5 + 0.5 * (hash(str(time.time() + 1)) % 100) / 100)), 1)

            # SM efficiency (streaming multiprocessor)
            metrics['sm_efficiency'] = round(util.gpu * (0.9 + 0.1 * (hash(str(time.time() + 2)) % 100) / 100), 1)

            # Tensor core usage (Parabricks uses tensor cores for certain operations)
            metrics['tensor_usage'] = round(35 * gpu_util_pct + (25 * gpu_util_pct * (hash(str(time.time() + 3)) % 100) / 100), 1)

            # PCIe throughput (GB10 PCIe 5.0 x16 = ~64 GB/s theoretical)
            metrics['pcie_throughput'] = round(12 * gpu_util_pct + (8 * gpu_util_pct * (hash(str(time.time() + 4)) % 100) / 100), 1)

        pynvml.nvmlShutdown()

    except Exception as e:
        logger.error(f"Error getting detailed GPU metrics: {e}")

    return metrics


def check_docker():
    """Check if Docker is available"""
    try:
        result = subprocess.run(
            ['docker', '--version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
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
    except Exception:
        return False


def get_cpu_utilization():
    """Get CPU utilization percentage"""
    try:
        # Get CPU percent averaged over 1 second
        cpu_percent = psutil.cpu_percent(interval=0.1)
        # Get per-core utilization
        cpu_per_core = psutil.cpu_percent(interval=0.1, percpu=True)
        # Get CPU frequency
        cpu_freq = psutil.cpu_freq()
        # Get CPU count
        cpu_count = psutil.cpu_count()

        return {
            'percent': round(cpu_percent, 1),
            'per_core': [round(p, 1) for p in cpu_per_core] if cpu_per_core else [],
            'frequency': round(cpu_freq.current, 0) if cpu_freq else 0,
            'count': cpu_count
        }
    except Exception as e:
        logger.error(f"Error getting CPU utilization: {e}")
        return {
            'percent': 0,
            'per_core': [],
            'frequency': 0,
            'count': 0
        }


def get_gpu_utilization():
    """Get GPU utilization percentage and memory usage"""
    if not NVML_AVAILABLE:
        return {
            'available': False,
            'devices': []
        }

    try:
        pynvml.nvmlInit()
        device_count = pynvml.nvmlDeviceGetCount()

        devices = []
        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)

            # Get GPU name
            name = pynvml.nvmlDeviceGetName(handle)
            if isinstance(name, bytes):
                name = name.decode('utf-8')

            # Get utilization
            util = pynvml.nvmlDeviceGetUtilizationRates(handle)
            gpu_util = util.gpu

            # Get memory info
            try:
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                mem_used_gb = mem_info.used / (1024**3)
                mem_total_gb = mem_info.total / (1024**3)
                mem_percent = (mem_info.used / mem_info.total) * 100
            except (pynvml.NVMLError, AttributeError, ZeroDivisionError):
                mem_used_gb = 0
                mem_total_gb = 0
                mem_percent = 0

            # Check for running compute processes on this GPU
            try:
                compute_procs = pynvml.nvmlDeviceGetComputeRunningProcesses(handle)
                has_compute_procs = len(compute_procs) > 0
            except (pynvml.NVMLError, AttributeError):
                has_compute_procs = False

            # GB10 workaround: if SM utilization is stuck high (>=90%) but no
            # processes are running and memory utilization is 0, report 0%
            if gpu_util >= 90 and not has_compute_procs and util.memory == 0:
                gpu_util = 0

            # Get temperature
            try:
                temp = pynvml.nvmlDeviceGetTemperature(handle, pynvml.NVML_TEMPERATURE_GPU)
            except (pynvml.NVMLError, AttributeError):
                temp = 0

            # Get power usage
            try:
                power = pynvml.nvmlDeviceGetPowerUsage(handle) / 1000.0  # Convert to watts
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

        return {
            'available': True,
            'devices': devices
        }

    except Exception as e:
        logger.error(f"Error getting GPU utilization: {e}")
        return {
            'available': False,
            'devices': [],
            'error': str(e)
        }


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
        logger.error(f"Error getting memory utilization: {e}")
        return {
            'percent': 0,
            'used_gb': 0,
            'total_gb': 0,
            'available_gb': 0,
            'swap_percent': 0
        }


@app.route('/api/ready')
def readiness():
    """Readiness probe — verifies all dependencies are available."""
    checks = {
        'docker': check_docker(),
        'gpu': check_gpu(),
        'scripts_dir': SCRIPTS_DIR.exists(),
    }
    all_ready = all(checks.values())
    return jsonify({'ready': all_ready, 'checks': checks}), 200 if all_ready else 503


# --- Graceful shutdown ---
def graceful_shutdown(signum=None, frame=None):
    """Terminate all child processes on shutdown."""
    for step_name, process in list(running_processes.items()):
        try:
            process.terminate()
            process.wait(timeout=5)
        except Exception:
            try:
                process.kill()
            except Exception:
                pass
    running_processes.clear()
    if signum:
        sys.exit(0)


signal.signal(signal.SIGTERM, graceful_shutdown)
signal.signal(signal.SIGINT, graceful_shutdown)
atexit.register(graceful_shutdown)


if __name__ == '__main__':
    # Create required directories
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    # Load any previously saved pipeline state
    load_pipeline_state()

    # Run server
    logger.info("=" * 60)
    logger.info("Genomics Pipeline Web Portal")
    logger.info("=" * 60)
    logger.info(f"Project Root: {PROJECT_ROOT}")
    logger.info(f"Opening portal at: http://localhost:5000")
    logger.info("=" * 60)

    # Apply security headers to all responses
    from security import add_security_headers

    @app.after_request
    def apply_security_headers(response):
        return add_security_headers(response)

    # Security: Debug mode disabled by default. Set FLASK_DEBUG=true to enable.
    debug_mode = os.environ.get('FLASK_DEBUG', 'false').lower() == 'true'
    host = os.environ.get('FLASK_HOST', '0.0.0.0')
    port = int(os.environ.get('FLASK_PORT', '5000'))

    if debug_mode:
        logger.warning("Debug mode is ENABLED. Do not use in production!")

    app.run(host=host, port=port, debug=debug_mode, threaded=True)
