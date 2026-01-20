#!/bin/bash
#
# HCLS AI Factory - Start All Services
# This script starts all pipeline services and the landing page
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRANSFER_DIR="$(dirname "$SCRIPT_DIR")"

echo "=============================================="
echo "  HCLS AI Factory - Starting All Services"
echo "=============================================="

# Function to check if a port is in use
port_in_use() {
    nc -z localhost $1 2>/dev/null
    return $?
}

# Function to start a service if not running
start_service() {
    local name="$1"
    local port="$2"
    local dir="$3"
    local cmd="$4"

    if port_in_use $port; then
        echo "[OK] $name (port $port) - Already running"
    else
        echo "[STARTING] $name (port $port)..."
        cd "$dir"
        nohup $cmd > /tmp/${name}.log 2>&1 &
        sleep 2
        if port_in_use $port; then
            echo "[OK] $name started successfully"
        else
            echo "[WARN] $name may still be starting..."
        fi
    fi
}

# Start Genomics Pipeline Portal (Port 5000)
start_service "genomics-portal" 5000 \
    "$TRANSFER_DIR/genomics-pipeline/web-portal" \
    "./venv/bin/python app/server.py"

# Start RAG/Chat API Portal (Port 5001)
start_service "rag-portal" 5001 \
    "$TRANSFER_DIR/rag-chat-pipeline" \
    "./venv/bin/python portal/app/server.py"

# Start RAG Chat Interface (Port 8501)
start_service "rag-chat" 8501 \
    "$TRANSFER_DIR/rag-chat-pipeline" \
    "./venv/bin/streamlit run app/chat_ui.py --server.port 8501 --server.address 0.0.0.0 --server.headless true"

# Start Drug Discovery Main (Port 8505)
start_service "drug-discovery" 8505 \
    "$TRANSFER_DIR/drug-discovery-pipeline" \
    "./venv/bin/streamlit run app/discovery_ui.py --server.port 8505"

# Start Drug Discovery Portal (Port 8510)
start_service "drug-portal" 8510 \
    "$TRANSFER_DIR/drug-discovery-pipeline" \
    "./venv/bin/streamlit run portal/app.py --server.port 8510 --server.headless true"

# Start Landing Page (Port 8080)
start_service "landing-page" 8080 \
    "$SCRIPT_DIR" \
    "./venv/bin/python server.py"

echo ""
echo "=============================================="
echo "  All services started!"
echo "=============================================="
echo ""
echo "  Landing Page:      http://localhost:8080"
echo "  Genomics Portal:   http://localhost:5000"
echo "  RAG/Chat API:      http://localhost:5001"
echo "  RAG Chat:          http://localhost:8501"
echo "  Drug Discovery:    http://localhost:8505"
echo "  Discovery Portal:  http://localhost:8510"
echo ""
