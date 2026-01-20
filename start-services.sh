#!/bin/bash
#
# HCLS AI Factory - Service Startup Script
# =========================================
# Brings all pipeline services online after a DGX Spark reboot.
#
# Usage: ./start-services.sh [options]
#
# Options:
#   --all           Start all services (default)
#   --landing       Start only landing page
#   --rag           Start only RAG/Chat pipeline
#   --drug          Start only Drug Discovery pipeline
#   --status        Show status of all services
#   --stop          Stop all services
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRANSFER_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# Service ports
LANDING_PORT=8080
RAG_CHAT_PORT=8501
DRUG_DISCOVERY_PORT=8505
DRUG_PORTAL_PORT=8510
MILVUS_PORT=19530

print_banner() {
    echo -e "${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════════════════════╗"
    echo "║                                                                              ║"
    echo "║   ${BOLD}HCLS AI FACTORY${NC}${CYAN}                                                          ║"
    echo "║   Precision Medicine to Drug Discovery                                       ║"
    echo "║                                                                              ║"
    echo "╚══════════════════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_header() {
    echo ""
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}  $1${NC}"
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
}

print_info() {
    echo -e "${YELLOW}  ➤ $1${NC}"
}

print_success() {
    echo -e "${GREEN}  ✓ $1${NC}"
}

print_error() {
    echo -e "${RED}  ✗ $1${NC}"
}

print_status() {
    local name=$1
    local port=$2
    local status=$3
    if [ "$status" = "running" ]; then
        echo -e "  ${GREEN}●${NC} ${BOLD}$name${NC} - Port $port ${GREEN}[RUNNING]${NC}"
    else
        echo -e "  ${RED}●${NC} ${BOLD}$name${NC} - Port $port ${RED}[STOPPED]${NC}"
    fi
}

get_host_ip() {
    hostname -I | awk '{print $1}'
}

check_port() {
    local port=$1
    lsof -i :$port > /dev/null 2>&1
}

wait_for_service() {
    local name=$1
    local port=$2
    local max_wait=${3:-30}

    for i in $(seq 1 $max_wait); do
        if check_port $port; then
            print_success "$name is ready on port $port"
            return 0
        fi
        sleep 1
    done
    print_error "$name failed to start on port $port"
    return 1
}

# ============================================================================
# SERVICE START FUNCTIONS
# ============================================================================

start_milvus() {
    print_header "Starting Milvus Vector Database"

    cd "$TRANSFER_DIR/rag-chat-pipeline"

    if docker ps | grep -q "milvus-standalone"; then
        print_success "Milvus is already running"
    else
        print_info "Starting Milvus container..."
        docker-compose up -d milvus
        sleep 5

        if docker ps | grep -q "milvus-standalone"; then
            print_success "Milvus started successfully"
        else
            print_error "Failed to start Milvus"
            return 1
        fi
    fi
}

start_landing_page() {
    print_header "Starting Landing Page"

    cd "$TRANSFER_DIR/landing-page"

    if check_port $LANDING_PORT; then
        print_success "Landing page is already running on port $LANDING_PORT"
    else
        print_info "Starting Flask server..."
        source venv/bin/activate 2>/dev/null || python3 -m venv venv && source venv/bin/activate
        pip install -q flask flask-cors 2>/dev/null

        nohup python3 server.py > /tmp/landing-page.log 2>&1 &

        wait_for_service "Landing Page" $LANDING_PORT 15
    fi
}

start_rag_chat() {
    print_header "Starting RAG/Chat Pipeline"

    cd "$TRANSFER_DIR/rag-chat-pipeline"

    # Start Milvus first
    start_milvus

    # Start Streamlit Chat UI
    if check_port $RAG_CHAT_PORT; then
        print_success "RAG Chat UI is already running on port $RAG_CHAT_PORT"
    else
        print_info "Starting Streamlit Chat UI..."
        source venv/bin/activate

        nohup streamlit run app/chat_ui.py \
            --server.port $RAG_CHAT_PORT \
            --server.address 0.0.0.0 \
            --server.headless true \
            > /tmp/rag-chat.log 2>&1 &

        wait_for_service "RAG Chat UI" $RAG_CHAT_PORT 20
    fi
}

start_drug_discovery() {
    print_header "Starting Drug Discovery Pipeline"

    cd "$TRANSFER_DIR/drug-discovery-pipeline"

    # Start main Discovery UI (port 8505)
    if check_port $DRUG_DISCOVERY_PORT; then
        print_success "Drug Discovery UI is already running on port $DRUG_DISCOVERY_PORT"
    else
        print_info "Starting Drug Discovery UI..."
        source venv/bin/activate

        nohup streamlit run app/discovery_ui.py \
            --server.port $DRUG_DISCOVERY_PORT \
            --server.address 0.0.0.0 \
            --server.headless true \
            > /tmp/drug-discovery.log 2>&1 &

        wait_for_service "Drug Discovery UI" $DRUG_DISCOVERY_PORT 20
    fi

    # Start Pipeline Portal (port 8510) - located in hls-orchestrator
    if check_port $DRUG_PORTAL_PORT; then
        print_success "Drug Discovery Portal is already running on port $DRUG_PORTAL_PORT"
    else
        print_info "Starting Drug Discovery Portal..."

        cd "$TRANSFER_DIR/hls-orchestrator"
        source "$TRANSFER_DIR/drug-discovery-pipeline/venv/bin/activate"

        nohup streamlit run portal/app.py \
            --server.port $DRUG_PORTAL_PORT \
            --server.address 0.0.0.0 \
            --server.headless true \
            > /tmp/drug-portal.log 2>&1 &

        wait_for_service "Drug Discovery Portal" $DRUG_PORTAL_PORT 20
        cd "$TRANSFER_DIR/drug-discovery-pipeline"
    fi
}

# ============================================================================
# STATUS AND STOP FUNCTIONS
# ============================================================================

show_status() {
    print_banner
    print_header "Service Status"

    HOST_IP=$(get_host_ip)

    echo ""
    echo -e "  ${BOLD}INFRASTRUCTURE${NC}"
    if docker ps | grep -q "milvus-standalone"; then
        print_status "Milvus Vector DB" $MILVUS_PORT "running"
    else
        print_status "Milvus Vector DB" $MILVUS_PORT "stopped"
    fi

    echo ""
    echo -e "  ${BOLD}WEB SERVICES${NC}"
    if check_port $LANDING_PORT; then
        print_status "Landing Page" $LANDING_PORT "running"
    else
        print_status "Landing Page" $LANDING_PORT "stopped"
    fi

    if check_port $RAG_CHAT_PORT; then
        print_status "RAG Chat UI" $RAG_CHAT_PORT "running"
    else
        print_status "RAG Chat UI" $RAG_CHAT_PORT "stopped"
    fi

    if check_port $DRUG_DISCOVERY_PORT; then
        print_status "Drug Discovery UI" $DRUG_DISCOVERY_PORT "running"
    else
        print_status "Drug Discovery UI" $DRUG_DISCOVERY_PORT "stopped"
    fi

    if check_port $DRUG_PORTAL_PORT; then
        print_status "Drug Discovery Portal" $DRUG_PORTAL_PORT "running"
    else
        print_status "Drug Discovery Portal" $DRUG_PORTAL_PORT "stopped"
    fi

    echo ""
    echo -e "  ${BOLD}ACCESS URLS${NC}"
    echo -e "  Landing Page:          ${CYAN}http://$HOST_IP:$LANDING_PORT${NC}"
    echo -e "  RAG Chat:              ${CYAN}http://$HOST_IP:$RAG_CHAT_PORT${NC}"
    echo -e "  Drug Discovery:        ${CYAN}http://$HOST_IP:$DRUG_DISCOVERY_PORT${NC}"
    echo -e "  Drug Discovery Portal: ${CYAN}http://$HOST_IP:$DRUG_PORTAL_PORT${NC}"
    echo ""
}

stop_services() {
    print_header "Stopping All Services"

    # Stop Streamlit processes
    print_info "Stopping Streamlit services..."
    pkill -f "streamlit run app/chat_ui.py" 2>/dev/null || true
    pkill -f "streamlit run app/discovery_ui.py" 2>/dev/null || true
    pkill -f "streamlit run portal/app.py" 2>/dev/null || true

    # Stop Landing Page
    print_info "Stopping Landing Page..."
    pkill -f "python3 server.py" 2>/dev/null || true

    # Optionally stop Milvus (commented out to preserve data)
    # print_info "Stopping Milvus..."
    # cd "$TRANSFER_DIR/rag-chat-pipeline" && docker-compose down

    print_success "Services stopped (Milvus kept running to preserve data)"
}

start_all() {
    print_banner

    HOST_IP=$(get_host_ip)

    start_landing_page
    start_rag_chat
    start_drug_discovery

    print_header "All Services Started!"

    echo ""
    echo -e "  ${BOLD}ACCESS URLS${NC}"
    echo ""
    echo -e "  ${GREEN}Landing Page${NC}          http://$HOST_IP:$LANDING_PORT"
    echo -e "  ${GREEN}RAG Chat${NC}              http://$HOST_IP:$RAG_CHAT_PORT"
    echo -e "  ${GREEN}Drug Discovery${NC}        http://$HOST_IP:$DRUG_DISCOVERY_PORT"
    echo -e "  ${GREEN}Drug Discovery Portal${NC} http://$HOST_IP:$DRUG_PORTAL_PORT"
    echo ""
    echo -e "  ${CYAN}Ready for demo!${NC}"
    echo ""
}

# ============================================================================
# MAIN
# ============================================================================

show_help() {
    echo "HCLS AI Factory - Service Manager"
    echo ""
    echo "Usage: ./start-services.sh [command]"
    echo ""
    echo "Commands:"
    echo "  (none)    Start all services (default)"
    echo "  --all     Start all services"
    echo "  --landing Start only landing page"
    echo "  --rag     Start only RAG/Chat pipeline"
    echo "  --drug    Start only Drug Discovery pipeline"
    echo "  --status  Show status of all services"
    echo "  --stop    Stop all services"
    echo "  --help    Show this help message"
    echo ""
}

case "${1:---all}" in
    --all|"")
        start_all
        ;;
    --landing)
        print_banner
        start_landing_page
        ;;
    --rag)
        print_banner
        start_rag_chat
        ;;
    --drug)
        print_banner
        start_drug_discovery
        ;;
    --status)
        show_status
        ;;
    --stop)
        print_banner
        stop_services
        ;;
    --help|-h)
        show_help
        ;;
    *)
        print_error "Unknown command: $1"
        show_help
        exit 1
        ;;
esac
