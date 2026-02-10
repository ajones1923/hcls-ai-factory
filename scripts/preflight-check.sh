#!/bin/bash
# ================================================================
# HCLS AI Factory - Environment Preflight Check
# ================================================================
# Validates all prerequisites before first deployment.
#
# Usage: ./scripts/preflight-check.sh
# ================================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

ERRORS=0
WARNINGS=0
PASSES=0

pass() {
    echo -e "  ${GREEN}PASS${NC} $1"
    PASSES=$((PASSES + 1))
}

fail() {
    echo -e "  ${RED}FAIL${NC} $1"
    ERRORS=$((ERRORS + 1))
}

warn() {
    echo -e "  ${YELLOW}WARN${NC} $1"
    WARNINGS=$((WARNINGS + 1))
}

echo "=========================================="
echo "HCLS AI Factory - Preflight Check"
echo "=========================================="
echo ""

# --------------------------------------------------------
# [1/12] Docker
# --------------------------------------------------------
echo -e "${CYAN}[1/12] Docker${NC}"
if command -v docker &>/dev/null; then
    DOCKER_VERSION=$(docker --version | grep -oP '\d+\.\d+' | head -1)
    DOCKER_MAJOR=$(echo "$DOCKER_VERSION" | cut -d. -f1)
    if [ "$DOCKER_MAJOR" -ge 24 ] 2>/dev/null; then
        pass "Docker $DOCKER_VERSION installed"
    else
        warn "Docker $DOCKER_VERSION found (recommended >= 24.0)"
    fi
else
    fail "Docker not installed"
fi

# --------------------------------------------------------
# [2/12] Docker Compose v2
# --------------------------------------------------------
echo -e "${CYAN}[2/12] Docker Compose${NC}"
if docker compose version &>/dev/null; then
    COMPOSE_VERSION=$(docker compose version --short 2>/dev/null || docker compose version | grep -oP '\d+\.\d+\.\d+')
    COMPOSE_MINOR=$(echo "$COMPOSE_VERSION" | cut -d. -f2)
    if [ "$(echo "$COMPOSE_VERSION" | cut -d. -f1)" -ge 2 ] && [ "$COMPOSE_MINOR" -ge 20 ] 2>/dev/null; then
        pass "Docker Compose $COMPOSE_VERSION (supports 'include')"
    else
        warn "Docker Compose $COMPOSE_VERSION (>= 2.20 recommended for 'include' support)"
    fi
else
    fail "Docker Compose v2 not available (install docker-compose-plugin)"
fi

# --------------------------------------------------------
# [3/12] Docker Daemon
# --------------------------------------------------------
echo -e "${CYAN}[3/12] Docker Daemon${NC}"
if docker info &>/dev/null; then
    pass "Docker daemon running"
else
    fail "Docker daemon not running (start with: sudo systemctl start docker)"
fi

# --------------------------------------------------------
# [4/12] NVIDIA Driver
# --------------------------------------------------------
echo -e "${CYAN}[4/12] NVIDIA Driver${NC}"
if command -v nvidia-smi &>/dev/null; then
    DRIVER_VERSION=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
    if [ -n "$DRIVER_VERSION" ]; then
        pass "NVIDIA driver $DRIVER_VERSION ($GPU_NAME)"
    else
        fail "nvidia-smi found but no GPU detected"
    fi
else
    fail "nvidia-smi not found (install NVIDIA drivers)"
fi

# --------------------------------------------------------
# [5/12] NVIDIA Container Toolkit
# --------------------------------------------------------
echo -e "${CYAN}[5/12] NVIDIA Container Toolkit${NC}"
if docker info 2>/dev/null | grep -q "nvidia"; then
    pass "NVIDIA Container Toolkit available"
elif command -v nvidia-ctk &>/dev/null; then
    pass "nvidia-ctk installed"
else
    warn "NVIDIA Container Toolkit not detected (required for GPU containers)"
fi

# --------------------------------------------------------
# [6/12] Disk Space
# --------------------------------------------------------
echo -e "${CYAN}[6/12] Disk Space${NC}"
AVAIL_GB=$(df -BG . 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G')
if [ -n "$AVAIL_GB" ]; then
    if [ "$AVAIL_GB" -lt 50 ] 2>/dev/null; then
        fail "Only ${AVAIL_GB}GB available (minimum 50GB required)"
    elif [ "$AVAIL_GB" -lt 100 ] 2>/dev/null; then
        warn "${AVAIL_GB}GB available (100GB+ recommended for full pipeline)"
    else
        pass "${AVAIL_GB}GB available"
    fi
else
    warn "Could not determine available disk space"
fi

# --------------------------------------------------------
# [7/12] System RAM
# --------------------------------------------------------
echo -e "${CYAN}[7/12] System RAM${NC}"
if command -v free &>/dev/null; then
    TOTAL_GB=$(free -g | awk '/^Mem:/{print $2}')
    if [ "$TOTAL_GB" -lt 32 ] 2>/dev/null; then
        fail "Only ${TOTAL_GB}GB RAM (minimum 32GB required)"
    elif [ "$TOTAL_GB" -lt 64 ] 2>/dev/null; then
        warn "${TOTAL_GB}GB RAM (DGX Spark has 128GB)"
    else
        pass "${TOTAL_GB}GB RAM"
    fi
else
    warn "Could not determine system memory"
fi

# --------------------------------------------------------
# [8/12] Python
# --------------------------------------------------------
echo -e "${CYAN}[8/12] Python${NC}"
if command -v python3 &>/dev/null; then
    PY_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    PY_MINOR=$(echo "$PY_VERSION" | cut -d. -f2)
    if [ "$PY_MINOR" -ge 10 ] 2>/dev/null; then
        pass "Python $PY_VERSION"
    else
        warn "Python $PY_VERSION found (3.10+ recommended)"
    fi
else
    fail "Python 3 not installed"
fi

# --------------------------------------------------------
# [9/12] Port Availability
# --------------------------------------------------------
echo -e "${CYAN}[9/12] Port Availability${NC}"
PORTS="8080 5000 5001 8501 8505 8510 19530 3000 9099 9100 9400 8001 8002"
BUSY_PORTS=""
for PORT in $PORTS; do
    if ss -tlnp 2>/dev/null | grep -q ":${PORT} " || \
       netstat -tlnp 2>/dev/null | grep -q ":${PORT} "; then
        BUSY_PORTS="$BUSY_PORTS $PORT"
    fi
done
if [ -z "$BUSY_PORTS" ]; then
    pass "All 13 service ports available"
else
    warn "Ports in use:${BUSY_PORTS} (may be existing services)"
fi

# --------------------------------------------------------
# [10/12] Environment Variables
# --------------------------------------------------------
echo -e "${CYAN}[10/12] Environment Variables${NC}"
ENV_OK=true
if [ -z "$ANTHROPIC_API_KEY" ]; then
    warn "ANTHROPIC_API_KEY not set (required for RAG Chat pipeline)"
    ENV_OK=false
fi
if [ -z "$NGC_API_KEY" ]; then
    warn "NGC_API_KEY not set (required for BioNeMo NIM services)"
    ENV_OK=false
fi
if [ "$ENV_OK" = true ]; then
    pass "Required API keys set"
fi

# --------------------------------------------------------
# [11/12] Docker Images
# --------------------------------------------------------
echo -e "${CYAN}[11/12] Key Docker Images${NC}"
IMAGES_MISSING=0
for IMAGE in \
    "nvcr.io/nvidia/k8s/dcgm-exporter:3.3.5-3.4.0-ubuntu22.04" \
    "grafana/grafana-oss:11.0.0" \
    "prom/prometheus:v2.52.0"; do
    if ! docker image inspect "$IMAGE" &>/dev/null; then
        IMAGES_MISSING=$((IMAGES_MISSING + 1))
    fi
done
if [ "$IMAGES_MISSING" -eq 0 ]; then
    pass "Key monitoring images pulled"
else
    warn "${IMAGES_MISSING} monitoring images not yet pulled (docker compose will pull them)"
fi

# --------------------------------------------------------
# [12/12] Git
# --------------------------------------------------------
echo -e "${CYAN}[12/12] Git${NC}"
if command -v git &>/dev/null; then
    GIT_VERSION=$(git --version | grep -oP '\d+\.\d+\.\d+')
    pass "Git $GIT_VERSION"
else
    warn "Git not installed (needed for development workflows)"
fi

# --------------------------------------------------------
# Summary
# --------------------------------------------------------
echo ""
echo "=========================================="
TOTAL=$((PASSES + WARNINGS + ERRORS))
echo -e "Results: ${GREEN}${PASSES} passed${NC}, ${YELLOW}${WARNINGS} warnings${NC}, ${RED}${ERRORS} errors${NC} (${TOTAL} checks)"
echo "=========================================="

if [ "$ERRORS" -gt 0 ]; then
    echo -e "${RED}Preflight check FAILED — fix errors above before deploying.${NC}"
    exit 1
else
    if [ "$WARNINGS" -gt 0 ]; then
        echo -e "${YELLOW}Preflight check PASSED with warnings.${NC}"
    else
        echo -e "${GREEN}Preflight check PASSED — ready to deploy!${NC}"
    fi
    exit 0
fi
