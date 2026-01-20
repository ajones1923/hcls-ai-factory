# HCLS AI Factory Landing Page

[![Flask](https://img.shields.io/badge/Flask-3.0-lightgrey?logo=flask)](https://flask.palletsprojects.com/)
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat&logo=python)](https://www.python.org/)

**Central Hub for the Healthcare & Life Sciences AI Factory**

> The landing page serves as the main entry point for all pipeline interfaces, providing real-time service health monitoring and quick access to all platform components.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Quick Start](#quick-start)
- [Service Health Monitoring](#service-health-monitoring)
- [Auto-Start Feature](#auto-start-feature)
- [API Endpoints](#api-endpoints)
- [Configuration](#configuration)
- [Directory Structure](#directory-structure)
- [Customization](#customization)
- [Troubleshooting](#troubleshooting)

---

## Overview

The Landing Page is a Flask-based web application that provides:

- Central navigation hub for all HCLS AI Factory services
- Real-time health monitoring for all pipeline interfaces
- Automatic startup of dependent services
- Platform statistics and workflow visualization

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         LANDING PAGE (Port 8080)                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │                      PIPELINE INTERFACES                             │  │
│   │                                                                      │  │
│   │   ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐           │  │
│   │   │ Genomics │  │ RAG/Chat │  │ RAG Chat │  │   Drug   │           │  │
│   │   │ Pipeline │  │   API    │  │ Interface│  │ Discovery│           │  │
│   │   │  :5000   │  │  :5001   │  │  :8501   │  │  :8505   │           │  │
│   │   └──────────┘  └──────────┘  └──────────┘  └──────────┘           │  │
│   └─────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │                   MONITORING & INFRASTRUCTURE                        │  │
│   │                                                                      │  │
│   │   ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐           │  │
│   │   │ Grafana  │  │Prometheus│  │  Node    │  │   DCGM   │           │  │
│   │   │  :3000   │  │  :9099   │  │ Exporter │  │ Exporter │           │  │
│   │   │          │  │          │  │  :9100   │  │  :9400   │           │  │
│   │   └──────────┘  └──────────┘  └──────────┘  └──────────┘           │  │
│   └─────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Features

### Central Navigation
- Quick-launch cards for all pipeline interfaces
- Color-coded service categories
- Direct links to monitoring dashboards

### Real-Time Health Monitoring
- Automatic health checks every 30 seconds
- Visual online/offline indicators
- Toast notifications for offline services
- Keyboard shortcut (R) to refresh status

### Auto-Start Dependencies
- Automatically starts Genomics Portal (port 5000) if offline
- Automatically starts RAG/Chat API (port 5001) if offline
- Runs health checks before service launch

### Platform Statistics
- Animated counters for key metrics
- Lines of code, variant counts, gene coverage
- End-to-end processing time display

### Responsive Design
- Modern gradient UI with animations
- Works on desktop and tablet
- Floating particle background effects

---

## Quick Start

### Option 1: Using Virtual Environment (Recommended)

```bash
cd /home/adam/transfer/landing-page

# Create virtual environment (first time only)
python3 -m venv venv
./venv/bin/pip install flask flask-cors requests

# Start the server
./venv/bin/python server.py
```

### Option 2: Start All Services

```bash
cd /home/adam/transfer/landing-page
./start-all.sh
```

This starts all 6 pipeline services plus the landing page.

### Option 3: Direct Python

```bash
cd /home/adam/transfer/landing-page
python3 server.py
```

**Access:** http://localhost:8080 or http://<your-ip>:8080

---

## Service Health Monitoring

The landing page monitors 10 services:

### Pipeline Interfaces

| Service | Port | Health Endpoint | Check Type |
|---------|------|-----------------|------------|
| Genomics Pipeline | 5000 | `/health` | HTTP |
| RAG/Chat API | 5001 | `/health` | HTTP |
| RAG Chat Interface | 8501 | `/healthz` | HTTP |
| Drug Discovery | 8505 | `/healthz` | HTTP |
| Discovery Portal | 8510 | `/healthz` | HTTP |

### Infrastructure

| Service | Port | Health Endpoint | Check Type |
|---------|------|-----------------|------------|
| Grafana | 3000 | `/api/health` | HTTP |
| Prometheus | 9099 | `/-/healthy` | HTTP |
| Node Exporter | 9100 | `/metrics` | HTTP |
| DCGM Exporter | 9400 | `/metrics` | HTTP |
| Milvus | 19530 | N/A | TCP |

### Health Check Behavior

1. **On Page Load:** All services checked immediately
2. **Every 30 Seconds:** Automatic refresh of all service statuses
3. **On Service Click:** Individual service check before navigation
4. **Manual Refresh:** Press `R` key to refresh all statuses

### Status Indicators

- **Green (Online):** Service responding normally
- **Red (Offline):** Service not responding or error
- **Toast Notification:** Appears when clicking offline service

---

## Auto-Start Feature

When the landing page server starts, it automatically checks and starts these services:

### Services Auto-Started

| Service | Port | Directory | Command |
|---------|------|-----------|---------|
| Genomics Portal | 5000 | `genomics-pipeline/web-portal` | `./venv/bin/python app/server.py` |
| RAG/Chat API | 5001 | `rag-chat-pipeline` | `./venv/bin/python portal/app/server.py` |

### How It Works

```python
# On server startup:
1. Check if port 5000 is open
2. If not, start Genomics Portal
3. Wait 2 seconds for startup
4. Check if port 5001 is open
5. If not, start RAG/Chat API
6. Wait 2 seconds for startup
7. Start Flask server on port 8080
```

### Startup Output

```
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║    HCLS AI FACTORY LANDING PAGE                                              ║
║                                                                              ║
║    Local:    http://localhost:8080                                           ║
║    Network:  http://192.168.68.107:8080                                      ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Checking dependent services...
[AUTO-START] Starting Genomics Portal (port 5000)...
[AUTO-START] Genomics Portal started successfully
[AUTO-START] Starting RAG/Chat API (port 5001)...
[AUTO-START] RAG/Chat API started successfully

 * Serving Flask app 'server'
 * Running on all addresses (0.0.0.0)
 * Running on http://127.0.0.1:8080
 * Running on http://192.168.68.107:8080
```

---

## API Endpoints

### Health Check Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Landing page health status |
| `/api/check-services` | GET | Check all services in parallel |
| `/api/check-service/<id>` | GET | Check single service by ID |

### Example Responses

**GET /health**
```json
{
  "status": "healthy",
  "service": "landing-page"
}
```

**GET /api/check-services**
```json
{
  "services": {
    "genomics": {
      "id": "genomics",
      "name": "Genomics Pipeline",
      "port": 5000,
      "online": true,
      "message": "Service responding"
    },
    "rag-portal": {
      "id": "rag-portal",
      "name": "RAG/Chat API",
      "port": 5001,
      "online": true,
      "message": "Service responding"
    }
    // ... other services
  },
  "summary": {
    "online": 10,
    "total": 10,
    "all_healthy": true
  }
}
```

**GET /api/check-service/genomics**
```json
{
  "id": "genomics",
  "name": "Genomics Pipeline",
  "port": 5000,
  "online": true,
  "message": "Service responding"
}
```

---

## Configuration

### Server Configuration

Edit `server.py` to modify:

```python
# Port configuration
PORT = int(os.environ.get('LANDING_PORT', 8080))

# Services to monitor
SERVICES = {
    'genomics': {'port': 5000, 'name': 'Genomics Pipeline', 'path': '/health'},
    'rag-portal': {'port': 5001, 'name': 'RAG/Chat API', 'path': '/health'},
    # ... add or modify services
}

# Services to auto-start
AUTO_START_SERVICES = [
    {
        'name': 'Genomics Portal',
        'port': 5000,
        'dir': TRANSFER_DIR / 'genomics-pipeline' / 'web-portal',
        'cmd': ['./venv/bin/python', 'app/server.py']
    },
    # ... add services to auto-start
]
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `LANDING_PORT` | 8080 | Port for landing page |

### UI Configuration

Edit `templates/index.html` to modify:
- Service cards and links
- Platform statistics
- Branding and styling

Edit `static/css/style.css` for visual customization.

Edit `static/js/main.js` for JavaScript behavior.

---

## Directory Structure

```
landing-page/
├── server.py                # Flask server with auto-start
├── start-all.sh            # Start all services script
├── requirements.txt        # Python dependencies
├── venv/                   # Python virtual environment
├── templates/
│   └── index.html         # Main landing page template
├── static/
│   ├── css/
│   │   └── style.css      # Custom styles (1,041 lines)
│   └── js/
│       └── main.js        # Frontend JavaScript (750 lines)
└── README.md              # This file
```

---

## Customization

### Adding a New Service Card

1. **Add to SERVICES in server.py:**
```python
SERVICES = {
    # ... existing services
    'new-service': {'port': 9999, 'name': 'New Service', 'path': '/health'},
}
```

2. **Add card in index.html:**
```html
<a href="http://{{ host }}:9999" target="_blank" class="pipeline-card">
    <div class="card-header">
        <div class="card-status online">
            <span class="status-indicator"></span>
            <span>Online</span>
        </div>
    </div>
    <div class="card-body">
        <h3>New Service</h3>
        <p>Description of the service</p>
    </div>
    <div class="card-footer">
        <span class="port-label">Port 9999</span>
        <span class="launch-text">Launch →</span>
    </div>
</a>
```

3. **Add to main.js serviceConfig:**
```javascript
const serviceConfig = {
    // ... existing services
    'new-service': { cardSelector: '.pipeline-card.new-service', port: 9999 }
};
```

### Changing the Footer

Edit line 429-431 in `templates/index.html`:
```html
<div class="footer-copy">
    HCLS AI Factory Demo
</div>
```

### Updating Statistics

Edit the stats section in `templates/index.html`:
```html
<div class="stat-value" data-value="35992">0</div>
<div class="stat-label">Lines of Code</div>
```

---

## Troubleshooting

### Landing Page Won't Start

```bash
# Check if port is in use
lsof -i :8080

# Kill existing process
fuser -k 8080/tcp

# Check Python version
./venv/bin/python --version

# Reinstall dependencies
./venv/bin/pip install -r requirements.txt
```

### Services Show Offline

```bash
# Check service status via API
curl http://localhost:8080/api/check-services

# Manually check service
curl http://localhost:5000/health

# Restart landing page (will auto-start dependencies)
fuser -k 8080/tcp
./venv/bin/python server.py
```

### Auto-Start Not Working

```bash
# Check if venvs exist
ls -la ../genomics-pipeline/web-portal/venv/
ls -la ../rag-chat-pipeline/venv/

# Check server.py logs
./venv/bin/python server.py 2>&1 | head -50

# Start services manually
cd ../genomics-pipeline/web-portal && ./venv/bin/python app/server.py &
cd ../rag-chat-pipeline && ./venv/bin/python portal/app/server.py &
```

### JavaScript Not Working

```bash
# Clear browser cache (hard refresh)
# Chrome/Firefox: Ctrl+Shift+R (Windows/Linux) or Cmd+Shift+R (Mac)

# Check browser console for errors
# Press F12 → Console tab
```

### Health Checks Timing Out

```bash
# Increase timeout in server.py
def check_service_health(service_id, service_info, host='localhost', timeout=5):
    # Changed from timeout=2 to timeout=5
```

---

## Dependencies

### Python Packages

```
flask>=2.0.0
flask-cors>=3.0.0
requests>=2.28.0
```

### Installation

```bash
cd /home/adam/transfer/landing-page
python3 -m venv venv
./venv/bin/pip install -r requirements.txt
```

---

## Related Documentation

- **Main README:** `/home/adam/transfer/hcls-ai-factory/README.md`
- **Genomics Pipeline:** `/home/adam/transfer/genomics-pipeline/README.md`
- **RAG/Chat Pipeline:** `/home/adam/transfer/rag-chat-pipeline/README.md`
- **Drug Discovery:** `/home/adam/transfer/drug-discovery-pipeline/README.md`

---

**Part of the HCLS AI Factory Platform**

*Built with Flask, JavaScript ES6, and CSS3*
