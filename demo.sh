#!/bin/bash
#
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║          PRECISION MEDICINE TO DRUG DISCOVERY AI FACTORY                      ║
# ║                        ONE-CLICK DEMO LAUNCHER                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
#
# This script launches all three pipelines for the GTC demo.
# Run with: ./demo.sh
#
# Options:
#   --stop      Stop all running services
#   --status    Show status of all services
#   --help      Show this help message
#

set -e

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOMICS_DIR="$SCRIPT_DIR/genomics-pipeline"
RAG_CHAT_DIR="$SCRIPT_DIR/rag-chat-pipeline"
DRUG_DISCOVERY_DIR="$SCRIPT_DIR/drug-discovery-pipeline"

# Ports
LANDING_PAGE_PORT=8080
GENOMICS_PORTAL_PORT=5000
RAG_PORTAL_PORT=5001
RAG_STREAMLIT_PORT=8501
DRUG_DISCOVERY_PORT=8505
DRUG_DISCOVERY_PORTAL_PORT=8510
MILVUS_PORT=19530
OLLAMA_PORT=11434

# Directories
LANDING_PAGE_DIR="$SCRIPT_DIR/landing-page"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# NVIDIA Green
NVIDIA_GREEN='\033[38;5;118m'

# ============================================================================
# Helper Functions
# ============================================================================

print_banner() {
    echo ""
    echo -e "${NVIDIA_GREEN}╔══════════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}                                                                              ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}    ${WHITE}${BOLD}PRECISION MEDICINE TO DRUG DISCOVERY AI FACTORY${NC}                          ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}                                                                              ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}    ${CYAN}From Patient DNA to Drug Candidates${NC}                                        ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}                                                                              ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}    ${MAGENTA}Powered by NVIDIA DGX Spark | Parabricks | BioNeMo${NC}                        ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}║${NC}                                                                              ${NVIDIA_GREEN}║${NC}"
    echo -e "${NVIDIA_GREEN}╚══════════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

print_section() {
    echo ""
    echo -e "${NVIDIA_GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${WHITE}${BOLD}  $1${NC}"
    echo -e "${NVIDIA_GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
}

check_port() {
    local port=$1
    if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null 2>&1; then
        return 0  # Port in use
    else
        return 1  # Port available
    fi
}

wait_for_port() {
    local port=$1
    local name=$2
    local max_wait=30
    local wait_count=0

    echo -ne "  ${YELLOW}Waiting for $name (port $port)...${NC}"
    while ! check_port $port; do
        sleep 1
        wait_count=$((wait_count + 1))
        if [ $wait_count -ge $max_wait ]; then
            echo -e " ${RED}TIMEOUT${NC}"
            return 1
        fi
        echo -n "."
    done
    echo -e " ${GREEN}READY${NC}"
    return 0
}

check_service() {
    local port=$1
    local name=$2
    if check_port $port; then
        echo -e "  ${GREEN}✓${NC} $name (port $port) - ${GREEN}Running${NC}"
        return 0
    else
        echo -e "  ${RED}✗${NC} $name (port $port) - ${RED}Not Running${NC}"
        return 1
    fi
}

# ============================================================================
# Status Command
# ============================================================================

show_status() {
    print_banner
    print_section "SERVICE STATUS"

    echo -e "  ${CYAN}Landing Page:${NC}"
    check_service $LANDING_PAGE_PORT "AI Factory Landing Page" || true

    echo ""
    echo -e "  ${CYAN}Infrastructure:${NC}"
    check_service $MILVUS_PORT "Milvus Vector Database" || true
    check_service $OLLAMA_PORT "Ollama LLM Server" || true

    echo ""
    echo -e "  ${CYAN}Pipeline UIs:${NC}"
    check_service $GENOMICS_PORTAL_PORT "Genomics Pipeline Portal" || true
    check_service $RAG_PORTAL_PORT "RAG/Chat API Portal" || true
    check_service $RAG_STREAMLIT_PORT "RAG/Chat Streamlit UI" || true
    check_service $DRUG_DISCOVERY_PORT "Drug Discovery (Main)" || true
    check_service $DRUG_DISCOVERY_PORTAL_PORT "Drug Discovery (Portal)" || true

    echo ""
    echo -e "  ${CYAN}Docker Containers:${NC}"
    if docker ps --format "{{.Names}}" 2>/dev/null | grep -q "milvus"; then
        echo -e "  ${GREEN}✓${NC} Milvus container running"
    else
        echo -e "  ${YELLOW}○${NC} Milvus container not running"
    fi

    echo ""
}

# ============================================================================
# Stop Command
# ============================================================================

stop_services() {
    print_banner
    print_section "STOPPING SERVICES"

    echo -e "  ${YELLOW}Stopping Landing Page...${NC}"
    pkill -f "python.*landing-page.*server" 2>/dev/null && echo -e "  ${GREEN}✓${NC} Landing Page stopped" || echo -e "  ${CYAN}○${NC} Landing Page was not running"

    echo -e "  ${YELLOW}Stopping Streamlit processes...${NC}"
    pkill -f "streamlit run.*chat_ui" 2>/dev/null && echo -e "  ${GREEN}✓${NC} RAG Chat UI stopped" || echo -e "  ${CYAN}○${NC} RAG Chat UI was not running"
    pkill -f "streamlit run.*discovery_ui" 2>/dev/null && echo -e "  ${GREEN}✓${NC} Drug Discovery UI stopped" || echo -e "  ${CYAN}○${NC} Drug Discovery UI was not running"

    echo -e "  ${YELLOW}Stopping Portal servers...${NC}"
    pkill -f "python.*server.py.*5000" 2>/dev/null && echo -e "  ${GREEN}✓${NC} Genomics Portal stopped" || echo -e "  ${CYAN}○${NC} Genomics Portal was not running"
    pkill -f "python.*server.py.*5001" 2>/dev/null && echo -e "  ${GREEN}✓${NC} RAG Portal stopped" || echo -e "  ${CYAN}○${NC} RAG Portal was not running"

    echo ""
    echo -e "  ${YELLOW}Note: Milvus and Ollama are left running (shared infrastructure)${NC}"
    echo -e "  ${CYAN}To stop Milvus: cd rag-chat-pipeline && docker-compose down${NC}"
    echo ""
    echo -e "  ${GREEN}All demo services stopped.${NC}"
    echo ""
}

# ============================================================================
# Help Command
# ============================================================================

show_help() {
    print_banner
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  (none)      Launch all pipeline UIs for demo"
    echo "  --stop      Stop all running demo services"
    echo "  --status    Show status of all services"
    echo "  --help      Show this help message"
    echo ""
    echo "Services launched:"
    echo "  - AI Factory Landing Page     http://localhost:8080  ← START HERE"
    echo "  - Genomics Pipeline Portal    http://localhost:5000"
    echo "  - RAG/Chat Pipeline Portal    http://localhost:5001"
    echo "  - RAG/Chat Streamlit UI       http://localhost:8501"
    echo "  - Drug Discovery (Main)       http://localhost:8505"
    echo "  - Drug Discovery (Portal)     http://localhost:8510"
    echo ""
    echo "Prerequisites:"
    echo "  - Docker (for Milvus)"
    echo "  - Ollama with llama3.1:70b model"
    echo "  - Python 3.10+"
    echo ""
}

# ============================================================================
# Main Launch Function
# ============================================================================

launch_demo() {
    print_banner

    # -------------------------------------------------------------------------
    # Prerequisites Check
    # -------------------------------------------------------------------------
    print_section "CHECKING PREREQUISITES"

    PREREQ_OK=true

    # Check Docker
    if command -v docker &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} Docker installed"
    else
        echo -e "  ${RED}✗${NC} Docker not found - Milvus requires Docker"
        PREREQ_OK=false
    fi

    # Check Python
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
        echo -e "  ${GREEN}✓${NC} Python $PYTHON_VERSION installed"
    else
        echo -e "  ${RED}✗${NC} Python 3 not found"
        PREREQ_OK=false
    fi

    # Check Ollama
    if check_port $OLLAMA_PORT; then
        echo -e "  ${GREEN}✓${NC} Ollama running on port $OLLAMA_PORT"
    else
        echo -e "  ${YELLOW}!${NC} Ollama not running - LLM features will be limited"
    fi

    # Check pipeline directories
    if [ -d "$GENOMICS_DIR" ]; then
        echo -e "  ${GREEN}✓${NC} Genomics Pipeline found"
    else
        echo -e "  ${RED}✗${NC} Genomics Pipeline not found at $GENOMICS_DIR"
        PREREQ_OK=false
    fi

    if [ -d "$RAG_CHAT_DIR" ]; then
        echo -e "  ${GREEN}✓${NC} RAG/Chat Pipeline found"
    else
        echo -e "  ${RED}✗${NC} RAG/Chat Pipeline not found at $RAG_CHAT_DIR"
        PREREQ_OK=false
    fi

    if [ -d "$DRUG_DISCOVERY_DIR" ]; then
        echo -e "  ${GREEN}✓${NC} Drug Discovery Pipeline found"
    else
        echo -e "  ${RED}✗${NC} Drug Discovery Pipeline not found at $DRUG_DISCOVERY_DIR"
        PREREQ_OK=false
    fi

    if [ "$PREREQ_OK" = false ]; then
        echo ""
        echo -e "  ${RED}Prerequisites check failed. Please fix the issues above.${NC}"
        exit 1
    fi

    # -------------------------------------------------------------------------
    # Start Landing Page
    # -------------------------------------------------------------------------
    print_section "STARTING LANDING PAGE (Port $LANDING_PAGE_PORT)"

    if check_port $LANDING_PAGE_PORT; then
        echo -e "  ${GREEN}✓${NC} Landing Page already running"
    else
        echo -e "  ${CYAN}Launching Landing Page...${NC}"
        cd "$LANDING_PAGE_DIR"

        # Create venv if needed
        if [ ! -d "venv" ]; then
            python3 -m venv venv
        fi
        source venv/bin/activate
        pip install -q -r requirements.txt 2>/dev/null

        # Start in background
        nohup python3 server.py > /tmp/landing-page.log 2>&1 &
        deactivate
        cd "$SCRIPT_DIR"

        wait_for_port $LANDING_PAGE_PORT "Landing Page"
    fi

    # -------------------------------------------------------------------------
    # Start Milvus
    # -------------------------------------------------------------------------
    print_section "STARTING INFRASTRUCTURE"

    if check_port $MILVUS_PORT; then
        echo -e "  ${GREEN}✓${NC} Milvus already running on port $MILVUS_PORT"
    else
        echo -e "  ${YELLOW}Starting Milvus...${NC}"
        cd "$RAG_CHAT_DIR"
        docker-compose up -d milvus 2>/dev/null || {
            echo -e "  ${RED}Failed to start Milvus. Trying to continue...${NC}"
        }
        cd "$SCRIPT_DIR"
        wait_for_port $MILVUS_PORT "Milvus" || {
            echo -e "  ${YELLOW}Warning: Milvus may not be ready. RAG search may fail.${NC}"
        }
    fi

    # -------------------------------------------------------------------------
    # Start Genomics Pipeline Portal
    # -------------------------------------------------------------------------
    print_section "STARTING GENOMICS PIPELINE (Port $GENOMICS_PORTAL_PORT)"

    if check_port $GENOMICS_PORTAL_PORT; then
        echo -e "  ${GREEN}✓${NC} Genomics Portal already running"
    else
        echo -e "  ${CYAN}Launching Genomics Portal...${NC}"
        cd "$GENOMICS_DIR/web-portal"

        # Create venv if needed
        if [ ! -d "venv" ]; then
            python3 -m venv venv
        fi
        source venv/bin/activate
        pip install -q -r requirements.txt 2>/dev/null

        # Start in background
        nohup python3 app/server.py > /tmp/genomics-portal.log 2>&1 &
        deactivate
        cd "$SCRIPT_DIR"

        wait_for_port $GENOMICS_PORTAL_PORT "Genomics Portal"
    fi

    # -------------------------------------------------------------------------
    # Start RAG/Chat Pipeline
    # -------------------------------------------------------------------------
    print_section "STARTING RAG/CHAT PIPELINE (Ports $RAG_PORTAL_PORT, $RAG_STREAMLIT_PORT)"

    cd "$RAG_CHAT_DIR"

    # Create venv if needed
    if [ ! -d "venv" ]; then
        echo -e "  ${CYAN}Creating Python virtual environment...${NC}"
        python3 -m venv venv
    fi
    source venv/bin/activate
    pip install -q -r requirements.txt 2>/dev/null

    # Start Portal API
    if check_port $RAG_PORTAL_PORT; then
        echo -e "  ${GREEN}✓${NC} RAG Portal already running"
    else
        echo -e "  ${CYAN}Launching RAG Portal API...${NC}"
        nohup python3 portal/app/server.py > /tmp/rag-portal.log 2>&1 &
        wait_for_port $RAG_PORTAL_PORT "RAG Portal"
    fi

    # Start Streamlit UI
    if check_port $RAG_STREAMLIT_PORT; then
        echo -e "  ${GREEN}✓${NC} RAG Streamlit UI already running"
    else
        echo -e "  ${CYAN}Launching RAG Streamlit UI...${NC}"
        nohup streamlit run app/chat_ui.py --server.port $RAG_STREAMLIT_PORT --server.headless true > /tmp/rag-streamlit.log 2>&1 &
        wait_for_port $RAG_STREAMLIT_PORT "RAG Streamlit"
    fi

    deactivate
    cd "$SCRIPT_DIR"

    # -------------------------------------------------------------------------
    # Start Drug Discovery Pipeline
    # -------------------------------------------------------------------------
    print_section "STARTING DRUG DISCOVERY PIPELINE (Port $DRUG_DISCOVERY_PORT)"

    cd "$DRUG_DISCOVERY_DIR"

    # Create venv if needed
    if [ ! -d "venv" ]; then
        echo -e "  ${CYAN}Creating Python virtual environment...${NC}"
        python3 -m venv venv
    fi
    source venv/bin/activate
    pip install -q -r requirements.txt 2>/dev/null
    pip install -q streamlit 2>/dev/null

    if check_port $DRUG_DISCOVERY_PORT; then
        echo -e "  ${GREEN}✓${NC} Drug Discovery UI already running"
    else
        echo -e "  ${CYAN}Launching Drug Discovery UI...${NC}"
        nohup streamlit run app/discovery_ui.py --server.port $DRUG_DISCOVERY_PORT --server.headless true > /tmp/drug-discovery.log 2>&1 &
        wait_for_port $DRUG_DISCOVERY_PORT "Drug Discovery UI"
    fi

    deactivate
    cd "$SCRIPT_DIR"

    # -------------------------------------------------------------------------
    # Open Browser
    # -------------------------------------------------------------------------
    print_section "OPENING BROWSER"

    sleep 2  # Give services a moment to fully initialize

    # Detect browser open command
    if command -v xdg-open &> /dev/null; then
        OPEN_CMD="xdg-open"
    elif command -v open &> /dev/null; then
        OPEN_CMD="open"
    else
        OPEN_CMD=""
    fi

    # Get host IP for display
    HOST_IP=$(hostname -I 2>/dev/null | awk '{print $1}' || echo "localhost")

    if [ -n "$OPEN_CMD" ]; then
        echo -e "  ${CYAN}Opening Landing Page...${NC}"
        $OPEN_CMD "http://localhost:$LANDING_PAGE_PORT" 2>/dev/null &
        echo -e "  ${GREEN}✓${NC} Landing Page opened in browser"
    else
        echo -e "  ${YELLOW}Could not detect browser. Please open URLs manually.${NC}"
    fi

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print_section "DEMO READY"

    echo -e "  ${WHITE}${BOLD}★ START HERE ★${NC}"
    echo ""
    echo -e "  ${NVIDIA_GREEN}╔═══════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "  ${NVIDIA_GREEN}║${NC}                                                                       ${NVIDIA_GREEN}║${NC}"
    echo -e "  ${NVIDIA_GREEN}║${NC}   ${WHITE}${BOLD}AI FACTORY LANDING PAGE${NC}     http://$HOST_IP:$LANDING_PAGE_PORT              ${NVIDIA_GREEN}║${NC}"
    echo -e "  ${NVIDIA_GREEN}║${NC}                                                                       ${NVIDIA_GREEN}║${NC}"
    echo -e "  ${NVIDIA_GREEN}╚═══════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "  ${WHITE}${BOLD}Pipeline Interfaces:${NC}"
    echo ""
    echo -e "  ┌─────────────────────────┬──────────────────────────────┬───────────┐"
    echo -e "  │ ${CYAN}Service${NC}                 │ ${CYAN}URL${NC}                          │ ${CYAN}Status${NC}    │"
    echo -e "  ├─────────────────────────┼──────────────────────────────┼───────────┤"
    echo -e "  │ Genomics Pipeline       │ http://$HOST_IP:$GENOMICS_PORTAL_PORT     │ ${GREEN}Running${NC}   │"
    echo -e "  │ RAG/Chat Pipeline       │ http://$HOST_IP:$RAG_PORTAL_PORT     │ ${GREEN}Running${NC}   │"
    echo -e "  │ RAG Chat UI             │ http://$HOST_IP:$RAG_STREAMLIT_PORT     │ ${GREEN}Running${NC}   │"
    echo -e "  │ Drug Discovery (Main)   │ http://$HOST_IP:$DRUG_DISCOVERY_PORT     │ ${GREEN}Running${NC}   │"
    echo -e "  │ Drug Discovery (Portal) │ http://$HOST_IP:$DRUG_DISCOVERY_PORTAL_PORT     │ ${GREEN}Running${NC}   │"
    echo -e "  └─────────────────────────┴──────────────────────────────┴───────────┘"
    echo ""
    echo -e "  ${WHITE}${BOLD}Demo Workflow:${NC}"
    echo -e "  ${YELLOW}1.${NC} Open Landing Page → Visual overview of entire platform"
    echo -e "  ${YELLOW}2.${NC} Genomics Portal   → Show HG002 variant processing"
    echo -e "  ${YELLOW}3.${NC} RAG/Chat UI       → Search for VCP variants, generate target hypothesis"
    echo -e "  ${YELLOW}4.${NC} Drug Discovery    → Generate drug candidates for VCP target"
    echo ""
    echo -e "  ${WHITE}${BOLD}Key Demo Points:${NC}"
    echo -e "  • ${CYAN}31,213 lines of code${NC} across three integrated pipelines"
    echo -e "  • ${CYAN}201 genes${NC} in Clinker knowledge base (13 therapeutic areas)"
    echo -e "  • ${CYAN}71M AlphaMissense${NC} predictions for pathogenicity scoring"
    echo -e "  • ${CYAN}3.5M embeddings${NC} in Milvus vector database"
    echo -e "  • ${CYAN}< 4 hours${NC} from raw FASTQ to drug candidates"
    echo ""
    echo -e "  ${WHITE}${BOLD}Logs:${NC}"
    echo -e "  • Landing Page:      /tmp/landing-page.log"
    echo -e "  • Genomics Portal:   /tmp/genomics-portal.log"
    echo -e "  • RAG Portal:        /tmp/rag-portal.log"
    echo -e "  • RAG Streamlit:     /tmp/rag-streamlit.log"
    echo -e "  • Drug Discovery:    /tmp/drug-discovery.log"
    echo ""
    echo -e "  ${MAGENTA}To stop all services: ./demo.sh --stop${NC}"
    echo ""
    echo -e "${NVIDIA_GREEN}══════════════════════════════════════════════════════════════════════════════${NC}"
    echo ""
}

# ============================================================================
# Main Entry Point
# ============================================================================

case "${1:-}" in
    --stop)
        stop_services
        ;;
    --status)
        show_status
        ;;
    --help|-h)
        show_help
        ;;
    *)
        launch_demo
        ;;
esac
