#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
#  HCLS AI Factory — Master Startup Script
#  Starts all services with correct environment, validates health, reports status
# ═══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

TRANSFER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="/tmp/hcls-factory"
mkdir -p "$LOG_DIR"

# ─── Configuration ────────────────────────────────────────────────────────────
# Set your API key here (or export ANTHROPIC_API_KEY before running)
ANTHROPIC_API_KEY="${ANTHROPIC_API_KEY:-}"
NIM_ALLOW_MOCK_FALLBACK="${NIM_ALLOW_MOCK_FALLBACK:-true}"
GRAFANA_PASSWORD="${GRAFANA_PASSWORD:-admin}"
LLM_MODEL="${LLM_MODEL:-claude-sonnet-4-20250514}"
LLM_PROVIDER="${LLM_PROVIDER:-anthropic}"

# ─── Colors ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

ok()   { echo -e "  ${GREEN}[OK]${NC}    $1"; }
fail() { echo -e "  ${RED}[FAIL]${NC}  $1"; }
warn() { echo -e "  ${YELLOW}[WARN]${NC}  $1"; }
info() { echo -e "  ${CYAN}[....]${NC}  $1"; }

# ─── Helpers ──────────────────────────────────────────────────────────────────
port_open() {
    nc -z localhost "$1" 2>/dev/null
}

wait_for_port() {
    local port="$1" name="$2" max_wait="${3:-15}"
    local waited=0
    while [ $waited -lt $max_wait ]; do
        if port_open "$port"; then
            return 0
        fi
        sleep 1
        waited=$((waited + 1))
    done
    return 1
}

check_venv() {
    local venv_dir="$1"
    if [ ! -f "$venv_dir/bin/python3" ] || [ ! -L "$venv_dir/bin/python3" ]; then
        return 1
    fi
    # Check if the symlink target exists
    if [ ! -e "$venv_dir/bin/python3" ]; then
        return 1
    fi
    return 0
}

rebuild_venv() {
    local dir="$1" name="$2"
    warn "$name venv broken — rebuilding..."
    (
        cd "$dir"
        python3 -m venv venv --clear
        ./venv/bin/pip install --quiet -r requirements.txt 2>/dev/null
        # Streamlit is needed for UI services
        if grep -q streamlit requirements.txt 2>/dev/null; then
            true  # already in requirements
        else
            ./venv/bin/pip install --quiet streamlit 2>/dev/null
        fi
    )
    if check_venv "$dir/venv"; then
        ok "$name venv rebuilt"
    else
        fail "$name venv rebuild failed"
        return 1
    fi
}

start_python_service() {
    local name="$1" port="$2" dir="$3" cmd="$4"

    if port_open "$port"; then
        ok "$name (:$port) already running"
        return 0
    fi

    info "Starting $name (:$port)..."

    # Validate venv
    if ! check_venv "$dir/venv"; then
        rebuild_venv "$dir" "$name" || return 1
    fi

    (
        cd "$dir"
        ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
        NIM_ALLOW_MOCK_FALLBACK="$NIM_ALLOW_MOCK_FALLBACK" \
        GRAFANA_PASSWORD="$GRAFANA_PASSWORD" \
        nohup $cmd > "$LOG_DIR/${name}.log" 2>&1 &
    )

    if wait_for_port "$port" "$name" 15; then
        ok "$name (:$port) started"
    else
        warn "$name (:$port) still starting — check $LOG_DIR/${name}.log"
    fi
}

start_streamlit_service() {
    local name="$1" port="$2" dir="$3" app_file="$4"

    if port_open "$port"; then
        ok "$name (:$port) already running"
        return 0
    fi

    info "Starting $name (:$port)..."

    if ! check_venv "$dir/venv"; then
        rebuild_venv "$dir" "$name" || return 1
    fi

    (
        cd "$dir"
        ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
        NIM_ALLOW_MOCK_FALLBACK="$NIM_ALLOW_MOCK_FALLBACK" \
        nohup ./venv/bin/streamlit run "$app_file" \
            --server.port "$port" \
            --server.address 0.0.0.0 \
            --server.headless true \
            > "$LOG_DIR/${name}.log" 2>&1 &
    )

    if wait_for_port "$port" "$name" 15; then
        ok "$name (:$port) started"
    else
        warn "$name (:$port) still starting — check $LOG_DIR/${name}.log"
    fi
}

# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD}  HCLS AI Factory — Precision Medicine to Drug Discovery${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo -e "  $(date '+%Y-%m-%d %H:%M:%S')  |  Logs: $LOG_DIR/"
echo ""

# ─── Step 0: Environment Check ───────────────────────────────────────────────
echo -e "${BOLD}  Step 0: Environment${NC}"

if [ -n "$ANTHROPIC_API_KEY" ]; then
    ok "ANTHROPIC_API_KEY set (${ANTHROPIC_API_KEY:0:12}...)"
    LLM_PROVIDER="anthropic"
    LLM_MODEL="claude-sonnet-4-20250514"
else
    warn "ANTHROPIC_API_KEY not set — RAG Chat will use Ollama (local)"
    LLM_PROVIDER="ollama"
    LLM_MODEL="llama3.1:8b"
fi

if [ "$NIM_ALLOW_MOCK_FALLBACK" = "true" ]; then
    info "NIM_ALLOW_MOCK_FALLBACK=true (drug discovery uses mock data)"
fi
echo ""

# ─── Step 1: Docker Services ─────────────────────────────────────────────────
echo -e "${BOLD}  Step 1: Docker Services${NC}"

# Milvus
if port_open 19530; then
    ok "Milvus (:19530) already running"
else
    info "Starting Milvus..."
    (cd "$TRANSFER_DIR/core/engines/precision-intelligence" && docker compose up -d milvus 2>/dev/null) || true
    wait_for_port 19530 "Milvus" 30 && ok "Milvus (:19530) started" || warn "Milvus may need manual start"
fi

# Ollama
if port_open 11434; then
    ok "Ollama (:11434) already running"
else
    info "Starting Ollama..."
    docker start ollama-compose 2>/dev/null || true
    wait_for_port 11434 "Ollama" 15 && ok "Ollama (:11434) started" || warn "Ollama not available"
fi

# Monitoring stack (Grafana, Prometheus, Node Exporter, DCGM)
for svc_info in "grafana:3000" "prometheus:9099" "node-exporter:9100" "dcgm-exporter:9400"; do
    svc="${svc_info%%:*}"
    port="${svc_info##*:}"
    if port_open "$port"; then
        ok "$svc (:$port) already running"
    else
        info "Starting $svc..."
        (cd "$TRANSFER_DIR/core/engines/therapeutic-discovery/monitoring" && \
            GRAFANA_PASSWORD="$GRAFANA_PASSWORD" docker compose up -d "$svc" 2>/dev/null) || true
        wait_for_port "$port" "$svc" 15 && ok "$svc (:$port) started" || warn "$svc may need manual start"
    fi
done

# Attu (Milvus UI)
if port_open 8000; then
    ok "Attu (:8000) already running"
else
    docker start attu 2>/dev/null && \
        wait_for_port 8000 "Attu" 10 && ok "Attu (:8000) started" || true
fi
echo ""

# ─── Step 2: Flask API Services ──────────────────────────────────────────────
echo -e "${BOLD}  Step 2: API Services${NC}"

start_python_service "genomics-portal" 5000 \
    "$TRANSFER_DIR/core/engines/genomic-foundation/web-portal" \
    "./venv/bin/python app/server.py"

start_python_service "rag-api" 5001 \
    "$TRANSFER_DIR/core/engines/precision-intelligence" \
    "./venv/bin/python portal/app/server.py"

# Set the LLM model once RAG API is up
if port_open 5001; then
    curl -s -X POST http://localhost:5001/api/model \
        -H "Content-Type: application/json" \
        -d "{\"model\":\"$LLM_MODEL\",\"provider\":\"$LLM_PROVIDER\"}" > /dev/null 2>&1
    ok "LLM model set: $LLM_MODEL ($LLM_PROVIDER)"
fi

start_python_service "landing-page" 8080 \
    "$TRANSFER_DIR/landing-page" \
    "./venv/bin/python server.py"

echo ""

# ─── Step 3: Streamlit UIs ───────────────────────────────────────────────────
echo -e "${BOLD}  Step 3: Streamlit UIs${NC}"

start_streamlit_service "rag-chat-ui" 8501 \
    "$TRANSFER_DIR/core/engines/precision-intelligence" \
    "app/chat_ui.py"

start_streamlit_service "drug-discovery-ui" 8505 \
    "$TRANSFER_DIR/core/engines/therapeutic-discovery" \
    "app/discovery_ui.py"

start_streamlit_service "discovery-portal" 8510 \
    "$TRANSFER_DIR/core/engines/therapeutic-discovery" \
    "portal/app.py"

echo ""

# ─── Step 4: Health Verification ─────────────────────────────────────────────
echo -e "${BOLD}  Step 4: Health Verification${NC}"
sleep 2

TOTAL=0
HEALTHY=0

check_health() {
    local name="$1" port="$2" path="$3"
    TOTAL=$((TOTAL + 1))
    local code
    code=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:${port}${path}" --max-time 5 2>/dev/null || echo "000")
    if [ "$code" = "200" ]; then
        ok "$name (:$port) — HTTP $code"
        HEALTHY=$((HEALTHY + 1))
    else
        fail "$name (:$port) — HTTP $code"
    fi
}

check_health "Landing Page"      8080  "/health"
check_health "Genomics Portal"   5000  "/health"
check_health "RAG/Chat API"      5001  "/health"
check_health "RAG Chat UI"       8501  "/_stcore/health"
check_health "Drug Discovery UI" 8505  "/_stcore/health"
check_health "Discovery Portal"  8510  "/_stcore/health"
check_health "Grafana"           3000  "/api/health"
check_health "Prometheus"        9099  "/-/healthy"
check_health "Node Exporter"     9100  "/metrics"

# Milvus is TCP-only
TOTAL=$((TOTAL + 1))
if port_open 19530; then
    ok "Milvus (:19530) — TCP open"
    HEALTHY=$((HEALTHY + 1))
else
    fail "Milvus (:19530) — not responding"
fi

echo ""

# ─── Summary ─────────────────────────────────────────────────────────────────
if [ "$HEALTHY" -eq "$TOTAL" ]; then
    COLOR="$GREEN"
else
    COLOR="$YELLOW"
fi

echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo -e "  ${COLOR}${HEALTHY}/${TOTAL} services healthy${NC}"
echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}"
echo ""
echo -e "  ${BOLD}Demo URLs:${NC}"
echo -e "  Landing Page:      ${CYAN}http://$(hostname -I | awk '{print $1}'):8080${NC}"
echo -e "  RAG Chat:          ${CYAN}http://$(hostname -I | awk '{print $1}'):8501${NC}"
echo -e "  Drug Discovery:    ${CYAN}http://$(hostname -I | awk '{print $1}'):8505${NC}"
echo -e "  Grafana:           ${CYAN}http://$(hostname -I | awk '{print $1}'):3000${NC}  (admin/$GRAFANA_PASSWORD)"
echo ""
echo -e "  ${BOLD}LLM:${NC} $LLM_MODEL ($LLM_PROVIDER)"
echo -e "  ${BOLD}Logs:${NC} $LOG_DIR/"
echo ""
