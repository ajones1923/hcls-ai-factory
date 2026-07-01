#!/bin/bash
###############################################################################
# Precision Autoimmune Intelligence Agent — Run Script
#
# Usage:
#   ./run.sh              # Start Streamlit UI only (port 8531)
#   ./run.sh --api        # Start FastAPI only (port 8532)
#   ./run.sh --both       # Start both UI and API
#   ./run.sh --setup      # Create collections + seed knowledge
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Load .env if present
if [ -f .env ]; then
    set -a
    source .env
    set +a
fi

# Create venv if needed
if [ ! -d venv ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    ./venv/bin/pip install --upgrade pip
    ./venv/bin/pip install -r requirements.txt
fi

PYTHON="./venv/bin/python"
STREAMLIT="./venv/bin/streamlit"

case "${1:-}" in
    --api)
        echo "Starting Autoimmune Agent API on port ${AUTO_API_PORT:-8532}..."
        $PYTHON -m uvicorn api.main:app \
            --host 0.0.0.0 \
            --port "${AUTO_API_PORT:-8532}" \
            --workers 2 \
            --log-level info
        ;;
    --both)
        echo "Starting Autoimmune Agent (UI: ${AUTO_STREAMLIT_PORT:-8531}, API: ${AUTO_API_PORT:-8532})..."

        # Track child PIDs for cleanup
        CHILD_PIDS=()

        cleanup() {
            echo ""
            echo "Shutting down Autoimmune Agent..."
            for pid in "${CHILD_PIDS[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    echo "Sending SIGTERM to PID $pid..."
                    kill -TERM "$pid" 2>/dev/null || true
                fi
            done

            # Wait up to 5 seconds for graceful shutdown
            for i in 1 2 3 4 5; do
                all_stopped=true
                for pid in "${CHILD_PIDS[@]}"; do
                    if kill -0 "$pid" 2>/dev/null; then
                        all_stopped=false
                        break
                    fi
                done
                if $all_stopped; then
                    break
                fi
                sleep 1
            done

            # Force-kill any remaining processes
            for pid in "${CHILD_PIDS[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    echo "Force-killing PID $pid..."
                    kill -KILL "$pid" 2>/dev/null || true
                fi
            done

            echo "Autoimmune Agent shutdown complete."
        }

        trap cleanup SIGTERM SIGINT EXIT

        $PYTHON -m uvicorn api.main:app \
            --host 0.0.0.0 \
            --port "${AUTO_API_PORT:-8532}" \
            --workers 2 \
            --log-level info &
        API_PID=$!
        CHILD_PIDS+=("$API_PID")
        echo "API started (PID: $API_PID)"

        $STREAMLIT run app/autoimmune_ui.py \
            --server.port "${AUTO_STREAMLIT_PORT:-8531}" \
            --server.address 0.0.0.0 \
            --server.headless true &
        UI_PID=$!
        CHILD_PIDS+=("$UI_PID")
        echo "UI started (PID: $UI_PID)"

        # Wait for any child to exit, then clean up
        wait -n 2>/dev/null || wait
        ;;
    --setup)
        echo "Setting up collections and seeding knowledge..."
        $PYTHON scripts/setup_collections.py --seed "$@"
        ;;
    *)
        echo "Starting Autoimmune Agent UI on port ${AUTO_STREAMLIT_PORT:-8531}..."
        $STREAMLIT run app/autoimmune_ui.py \
            --server.port "${AUTO_STREAMLIT_PORT:-8531}" \
            --server.address 0.0.0.0 \
            --server.headless true
        ;;
esac
