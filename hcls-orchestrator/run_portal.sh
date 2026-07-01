#!/bin/bash
# HLS Pipeline Portal Launcher

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORTAL_DIR="${SCRIPT_DIR}/portal"

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║  Healthcare & Life Sciences Pipeline Portal                   ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo ""

# Check for virtual environment
if [ -d "${SCRIPT_DIR}/venv" ]; then
    source "${SCRIPT_DIR}/venv/bin/activate"
elif [ -d "/home/adam/projects/hcls-ai-factory/core/engines/therapeutic-discovery/venv" ]; then
    source "/home/adam/projects/hcls-ai-factory/core/engines/therapeutic-discovery/venv/bin/activate"
fi

# Install requirements if needed
pip install -q streamlit requests plotly pandas 2>/dev/null || true

# Start the portal
echo "Starting portal on http://localhost:8510"
echo ""

cd "${PORTAL_DIR}"
streamlit run app.py \
    --server.port 8510 \
    --server.headless true \
    --browser.gatherUsageStats false \
    --theme.primaryColor "#76b900" \
    --theme.backgroundColor "#ffffff" \
    --theme.secondaryBackgroundColor "#f8f9fa"
