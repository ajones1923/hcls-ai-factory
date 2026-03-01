#!/usr/bin/env python3
"""
AI Factory Landing Page Server

A Flask server that serves the landing page for the Precision Medicine
to Drug Discovery AI Factory.
"""

import atexit
import json
import os
import signal
import socket
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import requests
from flask import Flask, jsonify, render_template, send_from_directory
from flask_cors import CORS
from loguru import logger

# Configure structured logging
logger.remove()
logger.add(sys.stderr, format="{time:YYYY-MM-DD HH:mm:ss} | {level:<7} | {message}", level="INFO")

app = Flask(__name__)
CORS(app, resources={
    r"/api/*": {
        "origins": os.environ.get('CORS_ORIGINS', 'http://localhost:8080').split(','),
        "methods": ["GET"],
        "allow_headers": ["Content-Type"]
    }
})

# Base directory for all pipelines
TRANSFER_DIR = Path(__file__).parent.parent

# Configuration
PORT = int(os.environ.get('LANDING_PORT', 8080))

# Service definitions with health check endpoints
SERVICES = {
    # Pipeline Interfaces
    'genomics': {'port': 5000, 'name': 'Genomics Pipeline', 'path': '/health'},
    'rag-portal': {'port': 5001, 'name': 'RAG/Chat API', 'path': '/health'},
    'rag-chat': {'port': 8501, 'name': 'RAG Chat Interface', 'path': '/healthz'},
    'drug-main': {'port': 8505, 'name': 'Drug Discovery', 'path': '/healthz'},
    'drug-portal': {'port': 8510, 'name': 'Discovery Portal', 'path': '/healthz'},
    # Monitoring & Infrastructure
    'grafana': {'port': 3000, 'name': 'Grafana', 'path': '/api/health'},
    'prometheus': {'port': 9099, 'name': 'Prometheus', 'path': '/-/healthy'},
    'node-exporter': {'port': 9100, 'name': 'Node Exporter', 'path': '/metrics'},
    'dcgm': {'port': 9400, 'name': 'DCGM Exporter', 'path': '/metrics'},
    # Intelligence Agents
    'cart-agent': {'port': 8521, 'name': 'CAR-T Intelligence Agent', 'path': '/healthz'},
    'imaging-agent': {'port': 8525, 'name': 'Imaging Intelligence Agent', 'path': '/healthz'},
    'onco-agent': {'port': 8526, 'name': 'Precision Oncology Agent', 'path': '/healthz'},
    # Vector Database
    'milvus': {'port': 19530, 'name': 'Milvus', 'path': None},  # TCP check only
}

# Services to auto-start when landing page starts
AUTO_START_SERVICES = [
    {
        'name': 'Genomics Portal',
        'port': 5000,
        'dir': TRANSFER_DIR / 'genomics-pipeline' / 'web-portal',
        'cmd': ['./venv/bin/python', 'app/server.py']
    },
    {
        'name': 'RAG/Chat API',
        'port': 5001,
        'dir': TRANSFER_DIR / 'rag-chat-pipeline',
        'cmd': ['./venv/bin/python', 'portal/app/server.py']
    },
]


def is_port_open(port, host='localhost', timeout=1):
    """Check if a port is open."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(timeout)
    result = sock.connect_ex((host, port))
    sock.close()
    return result == 0


def auto_start_services():
    """Auto-start dependent services if they're not running."""
    for service in AUTO_START_SERVICES:
        if not is_port_open(service['port']):
            logger.info(f"[AUTO-START] Starting {service['name']} (port {service['port']})...")
            try:
                log_path = f"/tmp/{service['name'].lower().replace('/', '-').replace(' ', '-')}.log"
                with open(log_path, 'w') as log_file:
                    subprocess.Popen(
                        service['cmd'],
                        cwd=service['dir'],
                        stdout=log_file,
                        stderr=subprocess.STDOUT,
                        start_new_session=True
                    )
                # Wait a moment for service to start
                time.sleep(2)
                if is_port_open(service['port']):
                    logger.info(f"[AUTO-START] {service['name']} started successfully")
                else:
                    logger.warning(f"[AUTO-START] {service['name']} may still be starting...")
            except Exception as e:
                logger.error(f"[AUTO-START] Failed to start {service['name']}: {e}")
        else:
            logger.info(f"[OK] {service['name']} (port {service['port']}) already running")


def get_host_ip():
    """Get the host IP address for service URLs."""
    try:
        # Try to get the IP that would be used to connect to an external host
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "localhost"


def check_service_health(service_id, service_info, host='localhost', timeout=2):
    """Check if a single service is online."""
    port = service_info['port']
    path = service_info.get('path')

    result = {
        'id': service_id,
        'name': service_info['name'],
        'port': port,
        'online': False,
        'message': ''
    }

    try:
        if path:
            # HTTP health check
            url = f"http://{host}:{port}{path}"
            response = requests.get(url, timeout=timeout)
            if response.status_code < 500:
                result['online'] = True
                result['message'] = 'Service responding'
            else:
                result['message'] = f'HTTP {response.status_code}'
        else:
            # TCP port check (for services like Milvus)
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(timeout)
            connection_result = sock.connect_ex((host, port))
            sock.close()
            if connection_result == 0:
                result['online'] = True
                result['message'] = 'Port open'
            else:
                result['message'] = 'Port closed'
    except requests.exceptions.Timeout:
        result['message'] = 'Timeout'
    except requests.exceptions.ConnectionError:
        result['message'] = 'Connection refused'
    except Exception as e:
        result['message'] = str(e)[:50]

    return result


@app.route('/')
def index():
    """Serve the main landing page."""
    host = get_host_ip()
    return render_template('index.html', host=host)


@app.route('/static/<path:filename>')
def static_files(filename):
    """Serve static files."""
    return send_from_directory('static', filename)


@app.route('/health')
def health():
    """Health check endpoint."""
    return {'status': 'healthy', 'service': 'landing-page'}


@app.route('/api/report-status')
def report_status():
    """Return pipeline report freshness metadata."""
    meta_path = Path(__file__).parent / 'static' / 'report_meta.json'
    if meta_path.exists():
        try:
            with open(meta_path) as f:
                return jsonify(json.load(f))
        except (json.JSONDecodeError, OSError):
            pass
    # Fallback to file mtime
    pdf_path = Path(__file__).parent / 'static' / 'VCP_Drug_Candidate_Report.pdf'
    if pdf_path.exists():
        mtime = datetime.fromtimestamp(pdf_path.stat().st_mtime)
        return jsonify({'generated_at': mtime.strftime('%Y-%m-%d %H:%M:%S'), 'has_context': False})
    return jsonify({'generated_at': None})


@app.route('/api/check-service/<service_id>')
def check_single_service(service_id):
    """Check health of a single service."""
    if service_id not in SERVICES:
        return jsonify({'error': f'Unknown service: {service_id}'}), 404

    host = get_host_ip()
    result = check_service_health(service_id, SERVICES[service_id], host)
    return jsonify(result)


@app.route('/api/check-services')
def check_all_services():
    """Check health of all services in parallel."""
    host = get_host_ip()
    results = {}

    # Use ThreadPoolExecutor for parallel health checks
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {
            executor.submit(check_service_health, sid, sinfo, host): sid
            for sid, sinfo in SERVICES.items()
        }

        for future in as_completed(futures):
            service_id = futures[future]
            try:
                result = future.result()
                results[service_id] = result
            except Exception as e:
                results[service_id] = {
                    'id': service_id,
                    'name': SERVICES[service_id]['name'],
                    'port': SERVICES[service_id]['port'],
                    'online': False,
                    'message': str(e)[:50]
                }

    # Calculate summary
    online_count = sum(1 for r in results.values() if r['online'])
    total_count = len(results)

    return jsonify({
        'services': results,
        'summary': {
            'online': online_count,
            'total': total_count,
            'all_healthy': online_count == total_count
        }
    })


@app.route('/api/ready')
def readiness():
    """Readiness probe — verifies critical pipeline services are available."""
    critical_services = ['genomics', 'rag-portal', 'drug-main']
    checks = {}
    for svc_id in critical_services:
        svc_info = SERVICES.get(svc_id, {})
        checks[svc_id] = is_port_open(svc_info.get('port', 0))
    all_ready = all(checks.values())
    return jsonify({'ready': all_ready, 'checks': checks}), 200 if all_ready else 503


# --- Graceful shutdown ---
def graceful_shutdown(signum=None, frame=None):
    """Clean shutdown handler."""
    if signum:
        sys.exit(0)


signal.signal(signal.SIGTERM, graceful_shutdown)
signal.signal(signal.SIGINT, graceful_shutdown)
atexit.register(graceful_shutdown)


if __name__ == '__main__':
    host_ip = get_host_ip()
    logger.info(f"HCLS AI FACTORY LANDING PAGE | Local: http://localhost:{PORT} | Network: http://{host_ip}:{PORT}")
    # Auto-start dependent services
    logger.info("Checking dependent services...")
    auto_start_services()

    app.run(host='0.0.0.0', port=PORT, debug=False)
