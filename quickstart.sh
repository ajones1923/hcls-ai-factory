#!/bin/bash
#
# HCLS AI Factory -- One-Command Quickstart
# ==========================================
# Sets up the entire precision medicine platform from scratch.
# Idempotent: safe to run multiple times.
#
# Usage: ./quickstart.sh [options]
#
# Options:
#   --help      Show this help message
#   --skip-venv Skip virtual environment creation and pip installs
#   --skip-data Skip the data status check
#
# Requirements:
#   - git, python3 (3.10+), docker, nvidia-smi, pip
#   - NVIDIA DGX Spark (ARM64) or x86_64 system with NVIDIA GPU
#
# Apache 2.0 Licensed | https://github.com/NVIDIA/hcls-ai-factory
#

set -eo pipefail

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STARTTIME=$(date +%s)

# Pipeline directories
RAG_CHAT_DIR="$SCRIPT_DIR/rag-chat-pipeline"
DRUG_DISCOVERY_DIR="$SCRIPT_DIR/drug-discovery-pipeline"
LANDING_PAGE_DIR="$SCRIPT_DIR/landing-page"
GENOMICS_PORTAL_DIR="$SCRIPT_DIR/genomics-pipeline/web-portal"

# Service ports
LANDING_PORT=8080
GENOMICS_PORT=5000
RAG_API_PORT=5001
RAG_UI_PORT=8501
DRUG_DISCOVERY_PORT=8505
DRUG_PORTAL_PORT=8510
MILVUS_PORT=19530
GRAFANA_PORT=3000
PROMETHEUS_PORT=9099

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m'
BOLD='\033[1m'
DIM='\033[2m'
NVIDIA_GREEN='\033[38;5;118m'

# Parse arguments
SKIP_VENV=false
SKIP_DATA=false
for arg in "$@"; do
    case "$arg" in
        --help|-h) show_help=true ;;
        --skip-venv) SKIP_VENV=true ;;
        --skip-data) SKIP_DATA=true ;;
    esac
done

# Step counter
STEP=0
TOTAL_STEPS=9
PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

# ============================================================================
# Helper Functions
# ============================================================================

print_banner() {
    echo ""
    echo -e "${NVIDIA_GREEN}"
    echo '    ╔══════════════════════════════════════════════════════════════════╗'
    echo '    ║                                                                ║'
    echo '    ║   ██╗  ██╗ ██████╗██╗     ███████╗                             ║'
    echo '    ║   ██║  ██║██╔════╝██║     ██╔════╝                             ║'
    echo '    ║   ███████║██║     ██║     ███████╗                             ║'
    echo '    ║   ██╔══██║██║     ██║     ╚════██║                             ║'
    echo '    ║   ██║  ██║╚██████╗███████╗███████║                             ║'
    echo '    ║   ╚═╝  ╚═╝ ╚═════╝╚══════╝╚══════╝                             ║'
    echo '    ║                                                                ║'
    echo -e "    ║   ${WHITE}${BOLD}AI Factory -- Quickstart${NC}${NVIDIA_GREEN}                                      ║"
    echo '    ║                                                                ║'
    echo -e "    ║   ${CYAN}From Patient DNA to Drug Candidates in < 5 Hours${NC}${NVIDIA_GREEN}               ║"
    echo -e "    ║   ${MAGENTA}Powered by NVIDIA DGX Spark (\$3,999)${NC}${NVIDIA_GREEN}                          ║"
    echo '    ║                                                                ║'
    echo '    ╚══════════════════════════════════════════════════════════════════╝'
    echo -e "${NC}"
}

print_section() {
    STEP=$((STEP + 1))
    echo ""
    echo -e "${NVIDIA_GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${WHITE}${BOLD}  [$STEP/$TOTAL_STEPS] $1${NC}"
    echo -e "${NVIDIA_GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
}

pass() {
    echo -e "  ${GREEN}[PASS]${NC} $1"
    PASS_COUNT=$((PASS_COUNT + 1))
}

fail() {
    echo -e "  ${RED}[FAIL]${NC} $1"
    FAIL_COUNT=$((FAIL_COUNT + 1))
}

warn() {
    echo -e "  ${YELLOW}[WARN]${NC} $1"
    WARN_COUNT=$((WARN_COUNT + 1))
}

info() {
    echo -e "  ${CYAN}[INFO]${NC} $1"
}

skip() {
    echo -e "  ${DIM}[SKIP]${NC} $1"
}

check_port() {
    local port=$1
    if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

show_help() {
    print_banner
    echo "Usage: ./quickstart.sh [options]"
    echo ""
    echo "Options:"
    echo "  --help        Show this help message"
    echo "  --skip-venv   Skip virtual environment creation and pip installs"
    echo "  --skip-data   Skip the data status check"
    echo ""
    echo "This script will:"
    echo "  1. Check prerequisites (git, python3, docker, nvidia-smi, pip)"
    echo "  2. Set up .env configuration file"
    echo "  3. Create Python virtual environments for each pipeline"
    echo "  4. Install Python dependencies"
    echo "  5. Show data download status"
    echo "  6. Start Milvus vector database"
    echo "  7. Start all pipeline services"
    echo "  8. Print service URLs and summary"
    echo ""
    echo "Runtime: < 5 minutes (with dependencies pre-installed)"
    echo ""
    exit 0
}

if [ "${show_help:-}" = "true" ]; then
    show_help
fi

# ============================================================================
# Main Quickstart Flow
# ============================================================================

print_banner

ARCH=$(uname -m)
info "Architecture: ${BOLD}$ARCH${NC} $([ "$ARCH" = "aarch64" ] && echo "(DGX Spark / ARM64)" || echo "(x86_64)")"
info "Working directory: ${BOLD}$SCRIPT_DIR${NC}"
echo ""

# --------------------------------------------------------------------------
# Step 1: Check Prerequisites
# --------------------------------------------------------------------------

print_section "Checking Prerequisites"

PREREQ_CRITICAL=0

# git
if command -v git &>/dev/null; then
    GIT_VER=$(git --version 2>&1 | awk '{print $3}')
    pass "git $GIT_VER"
else
    fail "git not found -- install with: sudo apt install git"
    PREREQ_CRITICAL=1
fi

# python3 (3.10+)
if command -v python3 &>/dev/null; then
    PY_VER=$(python3 --version 2>&1 | awk '{print $2}')
    PY_MAJOR=$(echo "$PY_VER" | cut -d. -f1)
    PY_MINOR=$(echo "$PY_VER" | cut -d. -f2)
    if [ "$PY_MAJOR" -ge 3 ] && [ "$PY_MINOR" -ge 10 ]; then
        pass "python3 $PY_VER (>= 3.10)"
    else
        fail "python3 $PY_VER found but >= 3.10 required"
        PREREQ_CRITICAL=1
    fi
else
    fail "python3 not found -- install with: sudo apt install python3"
    PREREQ_CRITICAL=1
fi

# pip
if command -v pip3 &>/dev/null || python3 -m pip --version &>/dev/null 2>&1; then
    PIP_VER=$(python3 -m pip --version 2>&1 | awk '{print $2}' || pip3 --version 2>&1 | awk '{print $2}')
    pass "pip $PIP_VER"
else
    fail "pip not found -- install with: sudo apt install python3-pip"
    PREREQ_CRITICAL=1
fi

# docker
if command -v docker &>/dev/null; then
    DOCKER_VER=$(docker --version 2>&1 | awk '{print $3}' | tr -d ',')
    pass "docker $DOCKER_VER"
    # Check docker compose (v2 plugin or standalone)
    if docker compose version &>/dev/null 2>&1; then
        COMPOSE_VER=$(docker compose version 2>&1 | awk '{print $NF}')
        pass "docker compose $COMPOSE_VER"
    elif command -v docker-compose &>/dev/null; then
        COMPOSE_VER=$(docker-compose --version 2>&1 | awk '{print $NF}')
        pass "docker-compose $COMPOSE_VER"
    else
        warn "docker compose not found -- Milvus requires docker compose"
    fi
else
    fail "docker not found -- install from https://docs.docker.com/get-docker/"
    PREREQ_CRITICAL=1
fi

# nvidia-smi
if command -v nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "unknown")
    DRIVER_VER=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1 || echo "unknown")
    pass "nvidia-smi -- GPU: $GPU_NAME, Driver: $DRIVER_VER"
else
    warn "nvidia-smi not found -- GPU acceleration will not be available"
fi

if [ $PREREQ_CRITICAL -gt 0 ]; then
    echo ""
    fail "Critical prerequisites missing. Please install them and re-run."
    exit 1
fi

echo ""
info "All critical prerequisites satisfied."

# --------------------------------------------------------------------------
# Step 2: Check .env Configuration
# --------------------------------------------------------------------------

print_section "Checking Environment Configuration"

if [ -f "$SCRIPT_DIR/.env" ]; then
    pass ".env file exists"

    # Check for placeholder values
    ENV_WARNINGS=0
    if grep -q "your-.*-here" "$SCRIPT_DIR/.env" 2>/dev/null; then
        ENV_WARNINGS=1
    fi

    # Check individual keys
    if grep -q "^ANTHROPIC_API_KEY=" "$SCRIPT_DIR/.env" && ! grep -q "^ANTHROPIC_API_KEY=your-" "$SCRIPT_DIR/.env"; then
        pass "ANTHROPIC_API_KEY is set"
    else
        warn "ANTHROPIC_API_KEY not configured -- RAG/Chat pipeline will not work"
        echo -e "       ${DIM}Get your key at: https://console.anthropic.com/${NC}"
    fi

    if grep -q "^NGC_API_KEY=" "$SCRIPT_DIR/.env" && ! grep -q "^NGC_API_KEY=your-" "$SCRIPT_DIR/.env"; then
        pass "NGC_API_KEY is set"
    else
        warn "NGC_API_KEY not configured -- Parabricks and BioNeMo will not work"
        echo -e "       ${DIM}Get your key at: https://ngc.nvidia.com/setup/api-key${NC}"
    fi

    if grep -q "^NVIDIA_API_KEY=" "$SCRIPT_DIR/.env" && ! grep -q "^NVIDIA_API_KEY=your-" "$SCRIPT_DIR/.env"; then
        pass "NVIDIA_API_KEY is set"
    else
        warn "NVIDIA_API_KEY not configured -- Cloud NIM services will not work"
        echo -e "       ${DIM}Get your key at: https://build.nvidia.com/${NC}"
    fi
else
    warn ".env file not found -- creating from .env.example"
    if [ -f "$SCRIPT_DIR/.env.example" ]; then
        cp "$SCRIPT_DIR/.env.example" "$SCRIPT_DIR/.env"
        pass "Created .env from .env.example"
        echo ""
        echo -e "  ${YELLOW}${BOLD}ACTION REQUIRED:${NC} Edit ${BOLD}.env${NC} and fill in your API keys:"
        echo -e "       ${WHITE}ANTHROPIC_API_KEY${NC}  -- https://console.anthropic.com/"
        echo -e "       ${WHITE}NGC_API_KEY${NC}        -- https://ngc.nvidia.com/setup/api-key"
        echo -e "       ${WHITE}NVIDIA_API_KEY${NC}     -- https://build.nvidia.com/"
        echo ""
        echo -e "       Then re-run: ${CYAN}./quickstart.sh${NC}"
    else
        fail ".env.example not found -- cannot create .env"
    fi
fi

# --------------------------------------------------------------------------
# Step 3: Create Python Virtual Environments
# --------------------------------------------------------------------------

print_section "Creating Python Virtual Environments"

if [ "$SKIP_VENV" = true ]; then
    skip "Skipping virtual environment creation (--skip-venv)"
else
    VENV_DIRS=(
        "$RAG_CHAT_DIR"
        "$DRUG_DISCOVERY_DIR"
        "$LANDING_PAGE_DIR"
        "$GENOMICS_PORTAL_DIR"
    )

    VENV_NAMES=(
        "rag-chat-pipeline"
        "drug-discovery-pipeline"
        "landing-page"
        "genomics-pipeline/web-portal"
    )

    for i in "${!VENV_DIRS[@]}"; do
        DIR="${VENV_DIRS[$i]}"
        NAME="${VENV_NAMES[$i]}"

        if [ ! -d "$DIR" ]; then
            warn "$NAME directory not found at $DIR -- skipping"
            continue
        fi

        if [ -d "$DIR/venv" ]; then
            pass "$NAME/venv already exists"
        else
            info "Creating virtual environment for $NAME..."
            if python3 -m venv "$DIR/venv" 2>/dev/null; then
                pass "$NAME/venv created"
            else
                fail "Failed to create venv for $NAME"
                echo -e "       ${DIM}Try: sudo apt install python3-venv${NC}"
            fi
        fi
    done
fi

# --------------------------------------------------------------------------
# Step 4: Install Python Dependencies
# --------------------------------------------------------------------------

print_section "Installing Python Dependencies"

if [ "$SKIP_VENV" = true ]; then
    skip "Skipping dependency installation (--skip-venv)"
else
    VENV_DIRS=(
        "$RAG_CHAT_DIR"
        "$DRUG_DISCOVERY_DIR"
        "$LANDING_PAGE_DIR"
        "$GENOMICS_PORTAL_DIR"
    )

    VENV_NAMES=(
        "rag-chat-pipeline"
        "drug-discovery-pipeline"
        "landing-page"
        "genomics-pipeline/web-portal"
    )

    for i in "${!VENV_DIRS[@]}"; do
        DIR="${VENV_DIRS[$i]}"
        NAME="${VENV_NAMES[$i]}"

        if [ ! -d "$DIR/venv" ]; then
            warn "$NAME/venv not found -- cannot install dependencies"
            continue
        fi

        if [ ! -f "$DIR/requirements.txt" ]; then
            warn "$NAME/requirements.txt not found -- skipping"
            continue
        fi

        info "Installing dependencies for $NAME..."
        set +e
        (
            source "$DIR/venv/bin/activate"
            pip install --upgrade pip -q 2>/dev/null
            pip install -r "$DIR/requirements.txt" -q 2>/dev/null
            deactivate
        )
        PIP_EXIT=$?
        set -e

        if [ $PIP_EXIT -eq 0 ]; then
            PKG_COUNT=$(source "$DIR/venv/bin/activate" && pip list 2>/dev/null | wc -l && deactivate)
            pass "$NAME -- dependencies installed ($PKG_COUNT packages)"
        else
            fail "$NAME -- pip install failed (check $DIR/requirements.txt)"
            echo -e "       ${DIM}Try manually: cd $NAME && source venv/bin/activate && pip install -r requirements.txt${NC}"
        fi
    done
fi

# --------------------------------------------------------------------------
# Step 5: Check Data Status
# --------------------------------------------------------------------------

print_section "Checking Pipeline Data Status"

if [ "$SKIP_DATA" = true ]; then
    skip "Skipping data status check (--skip-data)"
else
    if [ -x "$SCRIPT_DIR/setup-data.sh" ]; then
        info "Running setup-data.sh --status ..."
        echo ""
        set +e
        "$SCRIPT_DIR/setup-data.sh" --status 2>&1 | while IFS= read -r line; do
            echo "    $line"
        done
        DATA_EXIT=${PIPESTATUS[0]}
        set -e
        echo ""
        if [ "${DATA_EXIT:-0}" -eq 0 ]; then
            pass "Data status check complete"
        else
            warn "Data status check returned non-zero exit"
        fi
        echo ""
        info "To download all pipeline data, run: ${BOLD}./setup-data.sh --all${NC}"
        echo -e "       ${DIM}Warning: Full data download is ~856 GB (Stage 1 genomics data is large)${NC}"
        echo -e "       ${DIM}For a quick start, try: ./setup-data.sh --stage2 --stage3${NC}"
    else
        warn "setup-data.sh not found or not executable"
        echo -e "       ${DIM}Data must be downloaded manually for each pipeline.${NC}"
    fi
fi

# --------------------------------------------------------------------------
# Step 6: Start Milvus Vector Database
# --------------------------------------------------------------------------

print_section "Starting Milvus Vector Database"

if check_port $MILVUS_PORT; then
    pass "Milvus is already running on port $MILVUS_PORT"
else
    if [ -f "$RAG_CHAT_DIR/docker-compose.yml" ]; then
        info "Starting Milvus via docker compose..."
        set +e
        cd "$RAG_CHAT_DIR"
        if docker compose version &>/dev/null 2>&1; then
            docker compose up -d milvus 2>&1 | while IFS= read -r line; do echo "    $line"; done
        else
            docker-compose up -d milvus 2>&1 | while IFS= read -r line; do echo "    $line"; done
        fi
        cd "$SCRIPT_DIR"
        set -e

        # Wait for Milvus to be ready
        info "Waiting for Milvus to be ready..."
        MILVUS_WAIT=0
        MILVUS_MAX=60
        while ! check_port $MILVUS_PORT; do
            sleep 2
            MILVUS_WAIT=$((MILVUS_WAIT + 2))
            if [ $MILVUS_WAIT -ge $MILVUS_MAX ]; then
                fail "Milvus did not start within ${MILVUS_MAX}s"
                echo -e "       ${DIM}Check logs: docker logs milvus-standalone${NC}"
                break
            fi
            echo -ne "\r  ${CYAN}[INFO]${NC} Waiting for Milvus... ${MILVUS_WAIT}s / ${MILVUS_MAX}s"
        done

        if check_port $MILVUS_PORT; then
            echo ""
            pass "Milvus is running on port $MILVUS_PORT"
        fi
    else
        fail "docker-compose.yml not found in rag-chat-pipeline/"
        echo -e "       ${DIM}Milvus is required for the RAG/Chat pipeline.${NC}"
    fi
fi

# --------------------------------------------------------------------------
# Step 7: Start All Services
# --------------------------------------------------------------------------

print_section "Starting All Services"

if [ -x "$SCRIPT_DIR/start-services.sh" ]; then
    info "Running start-services.sh --all ..."
    echo ""
    set +e
    "$SCRIPT_DIR/start-services.sh" --all 2>&1 | while IFS= read -r line; do
        echo "    $line"
    done
    SERVICES_EXIT=${PIPESTATUS[0]}
    set -e
    echo ""
    if [ "${SERVICES_EXIT:-0}" -eq 0 ]; then
        pass "All services started successfully"
    else
        warn "start-services.sh exited with non-zero status (some services may have failed)"
        echo -e "       ${DIM}Check individual logs in /tmp/ for details.${NC}"
    fi
else
    fail "start-services.sh not found or not executable"
    echo -e "       ${DIM}Services must be started manually.${NC}"
fi

# --------------------------------------------------------------------------
# Step 8: Print Service Summary
# --------------------------------------------------------------------------

print_section "Service Summary"

HOST_IP=$(hostname -I 2>/dev/null | awk '{print $1}' || echo "localhost")

echo -e "  ${NVIDIA_GREEN}╔═══════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "  ${NVIDIA_GREEN}║${NC}                                                                       ${NVIDIA_GREEN}║${NC}"
echo -e "  ${NVIDIA_GREEN}║${NC}   ${WHITE}${BOLD}HCLS AI Factory -- Service URLs${NC}                                     ${NVIDIA_GREEN}║${NC}"
echo -e "  ${NVIDIA_GREEN}║${NC}                                                                       ${NVIDIA_GREEN}║${NC}"
echo -e "  ${NVIDIA_GREEN}╚═══════════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "  ┌───────────────────────────┬───────────────────────────────┬──────────┐"
echo -e "  │ ${CYAN}Service${NC}                     │ ${CYAN}URL${NC}                           │ ${CYAN}Status${NC}   │"
echo -e "  ├───────────────────────────┼───────────────────────────────┼──────────┤"

print_service_row() {
    local name=$1
    local port=$2
    local pad_name=$3

    if check_port $port; then
        STATUS="${GREEN}Running${NC} "
    else
        STATUS="${RED}Stopped${NC} "
    fi

    printf "  │ %-27s│ %-31s│ %b│\n" \
        "$name" "http://$HOST_IP:$port" "$STATUS"
}

print_service_row "AI Factory Landing Page" $LANDING_PORT 27
print_service_row "Genomics Pipeline Portal" $GENOMICS_PORT 27
print_service_row "RAG/Chat API" $RAG_API_PORT 27
print_service_row "RAG/Chat UI (Streamlit)" $RAG_UI_PORT 27
print_service_row "Drug Discovery UI" $DRUG_DISCOVERY_PORT 27
print_service_row "Drug Discovery Portal" $DRUG_PORTAL_PORT 27
echo -e "  ├───────────────────────────┼───────────────────────────────┼──────────┤"
print_service_row "Milvus Vector DB" $MILVUS_PORT 27
print_service_row "Grafana Monitoring" $GRAFANA_PORT 27
echo -e "  └───────────────────────────┴───────────────────────────────┴──────────┘"
echo ""

# --------------------------------------------------------------------------
# Step 9: Final Summary
# --------------------------------------------------------------------------

print_section "Quickstart Complete"

ENDTIME=$(date +%s)
ELAPSED=$((ENDTIME - STARTTIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo -e "  ${GREEN}${BOLD}Results:${NC}  ${GREEN}$PASS_COUNT passed${NC}  ${YELLOW}$WARN_COUNT warnings${NC}  ${RED}$FAIL_COUNT failed${NC}"
echo -e "  ${DIM}Completed in ${MINUTES}m ${SECONDS}s${NC}"
echo ""

if [ $FAIL_COUNT -gt 0 ]; then
    echo -e "  ${YELLOW}${BOLD}Some steps had failures.${NC} Review the output above for details."
    echo -e "  ${DIM}Most services should still be functional.${NC}"
    echo ""
fi

echo -e "  ${WHITE}${BOLD}Logs:${NC}"
echo -e "  ${DIM}Landing Page:     /tmp/landing-page.log${NC}"
echo -e "  ${DIM}RAG Chat:         /tmp/rag-chat.log${NC}"
echo -e "  ${DIM}RAG Portal:       /tmp/rag-portal.log${NC}"
echo -e "  ${DIM}Drug Discovery:   /tmp/drug-discovery.log${NC}"
echo -e "  ${DIM}Drug Portal:      /tmp/drug-portal.log${NC}"
echo -e "  ${DIM}Genomics Portal:  /tmp/genomics-portal.log${NC}"
echo ""
echo -e "  ${WHITE}${BOLD}Next Steps:${NC}"
echo -e "  ${CYAN}1.${NC} Download pipeline data:  ${BOLD}./setup-data.sh --stage2 --stage3${NC}"
echo -e "  ${CYAN}2.${NC} Run the full demo:       ${BOLD}./demo.sh${NC}"
echo -e "  ${CYAN}3.${NC} Check service status:    ${BOLD}./start-services.sh --status${NC}"
echo -e "  ${CYAN}4.${NC} Stop everything:         ${BOLD}./start-services.sh --stop${NC}"
echo ""
echo -e "${NVIDIA_GREEN}══════════════════════════════════════════════════════════════════════════════${NC}"
echo ""
echo -e "  ${NVIDIA_GREEN}${BOLD}Ready! Open http://localhost:$LANDING_PORT to begin.${NC}"
echo ""
echo -e "${NVIDIA_GREEN}══════════════════════════════════════════════════════════════════════════════${NC}"
echo ""
