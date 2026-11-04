# Genomics Pipeline Web Portal

A modern web interface for managing the complete FASTQ → VCF genomics workflow using NVIDIA Parabricks.

![Portal Interface](https://img.shields.io/badge/Interface-Web%20Portal-blue)
![Python](https://img.shields.io/badge/Python-3.11+-green)
![Flask](https://img.shields.io/badge/Flask-3.0-lightgrey)

## 🌟 Features

- **Visual Workflow Management**: Click-to-run interface for all pipeline steps
- **Real-Time Monitoring**: Live console output and status updates
- **System Status Dashboard**: Monitor Docker, GPU, and disk space
- **Configuration Management**: Edit pipeline settings through web UI
- **Log Viewer**: Browse and search pipeline logs
- **Progress Tracking**: Visual indicators for each workflow step
- **Responsive Design**: Works on desktop, tablet, and mobile

## 🚀 Quick Start

### Option 1: Direct Python (Recommended)

```bash
cd core/engines/genomic-foundation/web-portal
./start-portal.sh
```

Then open: **http://localhost:5000**

### Option 2: Docker Compose

```bash
cd core/engines/genomic-foundation/web-portal
docker-compose up -d
```

Then open: **http://localhost:5000**

## 📋 Prerequisites

- Python 3.11+ (for direct installation)
- Docker (if using containerized deployment)
- Same system requirements as the pipeline

## 🎯 Using the Portal

### 1. Check Prerequisites
- Click "Prerequisites Check" to verify Docker, GPU, and disk space
- Green badges indicate everything is ready

### 2. Run Workflow Steps
- Click each step in order from the left sidebar:
  1. **Prerequisites Check** - Verify system
  2. **NGC Login** - Authenticate (interactive in console)
  3. **Download Data** - Get GIAB HG002 FASTQ files
  4. **Setup Reference** - Prepare GRCh38 genome
  5. **Chr20 Test** - Quick validation run
  6. **Full Genome** - Complete analysis

### 3. Monitor Progress
- **Console Output Tab**: Real-time command output
- **System Status Panel**: Check resources
- **Data Files Panel**: Track file availability
- **Current Step Card**: Shows active operation

### 4. View Logs
- Switch to "Pipeline Logs" tab
- Select a log file from dropdown
- Browse detailed execution logs

### 5. Configure Settings
- Go to "Configuration" tab
- Adjust GPU count, memory mode, etc.
- Click "Save Configuration"

## 📊 Interface Overview

### Left Sidebar
- **Workflow Steps**: Clickable list of pipeline steps
- **System Status**: Docker, GPU, disk space indicators
- **Data Files**: Track FASTQ, reference, and VCF files

### Main Area
- **Console Output**: Live terminal output
- **Configuration**: Pipeline settings
- **Pipeline Logs**: Historical log files
- **Help**: Quick reference guide

### Status Indicators
- 🟢 **Green**: Available/Complete
- 🟡 **Yellow**: Running/Warning
- 🔴 **Red**: Error/Not Found
- ⚫ **Gray**: Not Started/Unknown

## 🔧 Configuration Options

Access via the "Configuration" tab:

| Setting | Default | Description |
|---------|---------|-------------|
| NUM_GPUS | 1 | Number of GPUs to use |
| LOW_MEMORY | 0 | Enable if GPU memory issues occur |
| PATIENT_ID | HG002 | Sample identifier |
| PB_IMG | parabricks:4.6.0-1 | Container image version |

## 📁 Directory Structure

```
web-portal/
├── app/
│   └── server.py           # Flask backend
├── templates/
│   └── index.html          # Main UI template
├── static/
│   ├── css/
│   │   └── style.css       # Custom styles
│   └── js/
│       └── app.js          # Frontend logic
├── start-portal.sh         # Startup script
├── requirements.txt        # Python dependencies
├── Dockerfile              # Container build
├── docker-compose.yml      # Container orchestration
└── README.md              # This file
```

## 🌐 API Endpoints

The portal provides a RESTful API:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Main portal page |
| `/api/status` | GET | Get pipeline and system status |
| `/api/config` | GET/POST | Get or update configuration |
| `/api/run/<step>` | GET | Execute workflow step |
| `/api/stop` | GET | Stop running process |
| `/api/logs/<type>` | GET | Get log file content |
| `/api/stream` | GET | Server-Sent Events stream |

## 🔒 Security Notes

**Authentication is always enforced.** Every endpoint except `/health`, `/healthz`,
`/api/ready`, `/api/status` and `/api/metrics` requires an `X-API-Key` header.

- Set `PORTAL_API_KEY` to a strong value. If you do not, the portal generates an
  **ephemeral key at startup**, writes it to a mode-0600 file (`PORTAL_API_KEY_FILE`, else
  `$XDG_RUNTIME_DIR/hcls-genomics-portal.key`, else the same name under the temp dir) and
  names that file in its log. The key itself is never written to the log. The portal does
  **not** fall back to an open API when no key is configured — it used to, and that was a
  remote-code-execution path through `/api/config`.
- The web UI asks for the key on first use and keeps it in `sessionStorage` for the tab.
- `/api/stream` and `/api/files/download/...` also accept `?api_key=` because `EventSource`
  and a browser-navigation download cannot send a header. Query strings are likelier to end
  up in proxy logs than headers; prefer the header for anything scripted.
- **CORS now defaults to a loopback allowlist** (`http://localhost:5000,http://127.0.0.1:5000`),
  not `*`. Set `CORS_ORIGINS` for a genuine cross-origin deployment. Do not set it to `*`.
- `/api/config` writes `../config/pipeline.env`, which the pipeline scripts `source` — bash
  executes it — and `PB_IMG` from it is the image `docker run` is given with host mounts and
  GPU access. The route therefore accepts only allowlisted keys, rejects any value with a
  line break or a shell metacharacter, requires `PB_IMG` to name an image in the portal's
  allowed set (`PORTAL_ALLOWED_IMAGES` overrides it), requires `GP`/`IN`/`REF`/`OUT`/`LOG` to
  resolve inside `PORTAL_CONFIG_PATH_ROOT` (default: the engine directory), and writes every
  value shell-quoted.
- `GET /api/config` returns only the declared pipeline keys. If you keep a registry
  credential in `pipeline.env`, it is not served — but it is still a plaintext secret on
  disk, so set its file mode accordingly.
- Uploads are capped at `PORTAL_MAX_UPLOAD_BYTES` (default 100 GB) and the cap is enforced
  while streaming, so a chunked upload cannot bypass it.
- Portal binds to `0.0.0.0:5000` (accessible from network). Bind it to `127.0.0.1` or
  firewall the port unless you intend network access.
- For production, use a reverse proxy (nginx) with SSL
- Requires access to Docker socket for container management — which is root-equivalent on
  the host. Treat the portal's API key as a host credential.

## 🐛 Troubleshooting

### Portal Won't Start
```bash
# Check Python version
python3 --version  # Should be 3.11+

# Check if port 5000 is available
lsof -i :5000

# View startup errors
./start-portal.sh
```

### Can't Connect to Portal
```bash
# Check if portal is running
ps aux | grep server.py

# Check firewall
sudo ufw status
sudo ufw allow 5000/tcp

# Try localhost
curl http://localhost:5000
```

### Steps Don't Execute
- Ensure scripts are executable: `chmod +x ../scripts/*.sh`
- Check Docker is running: `docker info`
- Verify file paths in `server.py`

### Real-Time Updates Not Working
- Check browser console for errors
- Verify Server-Sent Events support
- Try a different browser (Chrome/Firefox recommended)

## 🔄 Development

### Running in Debug Mode

```python
# In app/server.py, ensure:
app.run(host='0.0.0.0', port=5000, debug=True)
```

### Making Changes

1. **Backend**: Edit `app/server.py`
2. **Frontend**: Edit `templates/index.html`
3. **Styles**: Edit `static/css/style.css`
4. **JavaScript**: Edit `static/js/app.js`

Changes are automatically reloaded in debug mode.

### Adding New Steps

1. Add entry to `script_mapping` in `server.py`
2. Add list item to workflow steps in `index.html`
3. Create corresponding script in `../scripts/`

## 📦 Deployment Options

### Production Deployment

1. **With nginx reverse proxy:**
```nginx
server {
    listen 80;
    server_name genomics-portal.example.com;

    location / {
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

2. **With systemd service:**
```ini
[Unit]
Description=Genomics Pipeline Portal
After=network.target

[Service]
Type=simple
User=your-user
WorkingDirectory=/path/to/web-portal
ExecStart=/path/to/web-portal/venv/bin/python app/server.py
Restart=always

[Install]
WantedBy=multi-user.target
```

### Environment Variables

```bash
export FLASK_ENV=production
export FLASK_DEBUG=0
export PORT=5000

# Security (see "Security Notes" above)
export PORTAL_API_KEY=<strong secret>              # required in practice; auth is always on
export PORTAL_API_KEY_FILE=/run/portal/api.key     # where an ephemeral key is written
export CORS_ORIGINS=http://localhost:5000          # allowlist; never '*'
export PORTAL_ALLOWED_IMAGES=nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
export PORTAL_CONFIG_PATH_ROOT=/path/to/genomic-foundation
export PORTAL_MAX_UPLOAD_BYTES=107374182400        # 100 GB
```

## 🤝 Integration

The portal integrates with:
- Pipeline scripts in `../scripts/`
- Configuration in `../config/pipeline.env`
- Data directories in `../data/`
- Docker daemon for container management

## 📝 License

Part of the Genomics Pipeline project. See main README for details.

## 🆘 Support

For issues:
1. Check the Help tab in the portal
2. Review logs in Console Output tab
3. Check `../README.md` for pipeline documentation
4. Verify system prerequisites

---

**Built with:** Flask, Bootstrap 5, JavaScript ES6, Server-Sent Events
