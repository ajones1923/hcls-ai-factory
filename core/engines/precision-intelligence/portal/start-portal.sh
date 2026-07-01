#!/bin/bash
# Start RAG Chat Pipeline Portal

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo ""
echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}  RAG Chat Pipeline Portal${NC}"
echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
echo ""

# Check if venv exists
if [ ! -d "venv" ]; then
    echo -e "${YELLOW}Creating virtual environment...${NC}"
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Install dependencies if needed
if ! python -c "import flask" 2>/dev/null; then
    echo -e "${YELLOW}Installing portal dependencies...${NC}"
    pip install flask flask-cors psutil pynvml
fi

echo ""
echo -e "${GREEN}Starting portal at http://localhost:5001${NC}"
echo ""

# Start the portal
python portal/app/server.py
