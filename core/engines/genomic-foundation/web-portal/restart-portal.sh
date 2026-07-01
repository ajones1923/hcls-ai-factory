#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Genomics Pipeline Web Portal - Restart"
echo "=========================================="
echo ""

# Stop any existing portal process on port 5000
echo "Stopping existing portal..."
if lsof -ti :5000 > /dev/null 2>&1; then
    lsof -ti :5000 | xargs kill -9 2>/dev/null
    echo "  Killed existing process(es) on port 5000"
    sleep 1
else
    echo "  No existing process found on port 5000"
fi

# Verify port is free
if lsof -ti :5000 > /dev/null 2>&1; then
    echo "Error: Port 5000 is still in use"
    exit 1
fi
echo "  Port 5000 is free"
echo ""

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed"
    exit 1
fi

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install/upgrade dependencies
echo "Installing dependencies..."
pip install --upgrade pip > /dev/null 2>&1
pip install -r requirements.txt > /dev/null 2>&1

echo ""
echo "=========================================="
echo "Starting Web Portal..."
echo "=========================================="
echo ""
echo "Portal URL: http://localhost:5000"
echo "Press Ctrl+C to stop"
echo ""

# Start the server
python3 app/server.py
