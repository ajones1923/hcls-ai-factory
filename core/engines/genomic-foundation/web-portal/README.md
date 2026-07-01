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
cd genomics-pipeline/web-portal
./start-portal.sh
```

Then open: **http://localhost:5000**

### Option 2: Docker Compose

```bash
cd genomics-pipeline/web-portal
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

- Portal binds to `0.0.0.0:5000` (accessible from network)
- For production, use a reverse proxy (nginx) with SSL
- Requires access to Docker socket for container management
- No authentication by default - add if needed

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
