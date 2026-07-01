#!/usr/bin/env bash
# Stand up the TSC Intelligence Engine on the host (portal + 3 surfaces), persistent via nohup.
#   scripts/serve.sh            # offline reasoning (no API spend) — default
#   TSC_OFFLINE=0 scripts/serve.sh   # real Claude reasoning (needs key in .env; bills the key)
set -e
cd "$(dirname "$0")/.."
PY="venv/bin/python"; [ -x "$PY" ] || PY="python3"
export TSC_OFFLINE="${TSC_OFFLINE:-1}"
export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
mkdir -p data/logs data/state

start() { nohup "$PY" $2 > "data/logs/$1.log" 2>&1 & echo $! > "data/state/$1.pid"; echo "  $1  -> pid $(cat data/state/$1.pid)  (log: data/logs/$1.log)"; }

echo "Starting TSC Intelligence Engine (TSC_OFFLINE=$TSC_OFFLINE)..."
start api       "-m uvicorn api.main:app --host 0.0.0.0 --port 8560"
start briefing  "-m streamlit run app/briefing_app.py  --server.port 8561 --server.address 0.0.0.0 --server.headless true"
start dashboard "-m streamlit run app/dashboard_app.py --server.port 8562 --server.address 0.0.0.0 --server.headless true"
start alerts    "-m streamlit run app/alerts_app.py    --server.port 8563 --server.address 0.0.0.0 --server.headless true"
echo ""
echo "  Portal:    http://localhost:8560   (API docs at /docs)"
echo "  Briefing:  http://localhost:8561"
echo "  Dashboard: http://localhost:8562"
echo "  Alerts:    http://localhost:8563"
echo "Stop with: scripts/stop.sh"
