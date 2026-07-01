#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}  $1${NC}"
    echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo ""
}

print_info() {
    echo -e "${YELLOW}ℹ️  $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

show_help() {
    echo "RAG Chat Pipeline - Genomic Evidence Search"
    echo ""
    echo "Usage: ./run.sh <command>"
    echo ""
    echo "Commands:"
    echo "  setup       Install Python dependencies"
    echo "  start       Start Milvus database (Docker)"
    echo "  stop        Stop Milvus database"
    echo "  status      Check service status"
    echo "  ingest      Ingest VCF into vector database"
    echo "  chat        Start Streamlit chat interface"
    echo "  portal      Start web portal (port 5001)"
    echo "  demo-setup  Setup demo targets for GTC"
    echo "  all         Start services and chat interface"
    echo "  help        Show this help message"
    echo ""
    echo "Examples:"
    echo "  ./run.sh setup              # First-time setup"
    echo "  ./run.sh start              # Start Milvus"
    echo "  ./run.sh ingest --limit 10000  # Ingest first 10k variants"
    echo "  ./run.sh chat               # Start chat UI"
    echo "  ./run.sh portal             # Start web portal"
    echo ""
}

cmd_setup() {
    print_header "Setting up RAG Chat Pipeline"

    # Check Python
    if ! command -v python3 &> /dev/null; then
        print_error "Python 3 is required but not installed"
        exit 1
    fi

    # Create virtual environment if not exists
    if [ ! -d "venv" ]; then
        print_info "Creating virtual environment..."
        python3 -m venv venv
    fi

    # Activate and install dependencies
    print_info "Installing dependencies..."
    source venv/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt

    # Create .env if not exists
    if [ ! -f ".env" ]; then
        print_info "Creating .env from template..."
        cp .env.example .env
        print_info "Please edit .env and add your API keys"
    fi

    # Create data directories
    mkdir -p data/input data/output data/cache

    print_success "Setup complete!"
    echo ""
    echo "Next steps:"
    echo "  1. Edit .env and add your HF_TOKEN (Hugging Face)"
    echo "     Get token at: https://huggingface.co/settings/tokens"
    echo "  2. Accept Llama license at: https://huggingface.co/meta-llama/Llama-3.1-8B-Instruct"
    echo "  3. Run: ./run.sh start"
    echo "  4. Run: ./run.sh ingest"
    echo "  5. Run: ./run.sh chat"
}

cmd_start() {
    print_header "Starting RAG Chat Services"

    if ! command -v docker &> /dev/null; then
        print_error "Docker is required but not installed"
        exit 1
    fi

    # Check for .env
    if [ ! -f ".env" ]; then
        print_error ".env file not found. Run ./run.sh setup first"
        exit 1
    fi

    # Load .env for HF_TOKEN
    source .env

    print_info "Starting Milvus vector database..."
    docker-compose up -d milvus

    print_info "Waiting for Milvus to be ready..."
    sleep 10

    # Check if Milvus is healthy
    if docker-compose ps | grep -q "milvus.*Up"; then
        print_success "Milvus is running"
    else
        print_error "Milvus failed to start"
        docker-compose logs milvus
        exit 1
    fi

    # Check LLM provider - only start vLLM if configured
    LLM_PROVIDER="${LLM_PROVIDER:-anthropic}"

    if [ "$LLM_PROVIDER" = "vllm" ]; then
        print_info "Starting vLLM local LLM server..."
        print_info "Model: ${VLLM_MODEL:-meta-llama/Llama-3.1-8B-Instruct}"
        print_info "This may take a few minutes on first run (downloading model)..."
        docker-compose up -d vllm

        print_info "Waiting for vLLM to load model..."
        # Wait up to 5 minutes for vLLM to be ready
        for i in {1..30}; do
            if curl -s http://localhost:8080/health > /dev/null 2>&1; then
                print_success "vLLM is ready!"
                break
            fi
            echo "  Waiting for model to load... ($i/30)"
            sleep 10
        done
    elif [ "$LLM_PROVIDER" = "ollama" ]; then
        print_info "Using Ollama LLM (ensure it's running externally)"
        print_info "Ollama host: ${OLLAMA_HOST:-http://localhost:11434}"
    else
        print_info "Using Anthropic Claude API"
        if [ -z "$ANTHROPIC_API_KEY" ]; then
            print_error "ANTHROPIC_API_KEY environment variable not set!"
            print_info "Set it with: export ANTHROPIC_API_KEY='your-key-here'"
        fi
    fi

    # Final status check
    echo ""
    print_success "All services started!"
    echo ""
    echo "Services:"
    echo "  • Milvus: localhost:19530"
    echo "  • LLM:    $LLM_PROVIDER"
    echo ""
    echo "To start the chat UI: ./run.sh chat"
}

cmd_stop() {
    print_header "Stopping Services"
    docker-compose down
    print_success "Services stopped"
}

cmd_status() {
    print_header "Service Status"

    echo "Docker containers:"
    docker-compose ps

    echo ""
    echo "Service health checks:"

    # Milvus
    source venv/bin/activate 2>/dev/null || true
    python3 -c "
from pymilvus import connections
try:
    connections.connect('default', host='localhost', port='19530')
    print('  ✅ Milvus: connected')
    connections.disconnect('default')
except Exception as e:
    print(f'  ❌ Milvus: {e}')
" 2>/dev/null || echo "  ⚠️  Milvus: Python environment not set up"

    # vLLM
    if curl -s http://localhost:8080/health > /dev/null 2>&1; then
        MODEL=$(curl -s http://localhost:8080/v1/models | python3 -c "import sys,json; print(json.load(sys.stdin)['data'][0]['id'])" 2>/dev/null || echo "unknown")
        echo "  ✅ vLLM: running (model: $MODEL)"
    else
        echo "  ❌ vLLM: not responding"
    fi

    # Attu
    if curl -s http://localhost:8000 > /dev/null 2>&1; then
        echo "  ✅ Attu UI: http://localhost:8000"
    else
        echo "  ❌ Attu UI: not responding"
    fi
}

cmd_ingest() {
    print_header "Ingesting VCF into Vector Database"

    source venv/bin/activate 2>/dev/null || {
        print_error "Virtual environment not found. Run ./run.sh setup first"
        exit 1
    }

    # Pass all arguments to the ingestion script
    python scripts/ingest_vcf.py "$@"
}

cmd_chat() {
    print_header "Starting RAG Chat Interface"

    source venv/bin/activate 2>/dev/null || {
        print_error "Virtual environment not found. Run ./run.sh setup first"
        exit 1
    }

    python scripts/run_chat.py "$@"
}

cmd_portal() {
    print_header "Starting RAG Chat Pipeline Portal"

    source venv/bin/activate 2>/dev/null || {
        print_error "Virtual environment not found. Run ./run.sh setup first"
        exit 1
    }

    print_info "Portal available at: http://localhost:5001"
    python portal/app/server.py
}

cmd_demo_setup() {
    print_header "Setting Up Demo Targets"

    source venv/bin/activate 2>/dev/null || {
        print_error "Virtual environment not found. Run ./run.sh setup first"
        exit 1
    }

    python scripts/setup_demo_targets.py
}

cmd_all() {
    cmd_start
    echo ""
    print_info "Starting chat interface..."
    cmd_chat
}

# Main command router
case "${1:-help}" in
    setup)
        cmd_setup
        ;;
    start)
        cmd_start
        ;;
    stop)
        cmd_stop
        ;;
    status)
        cmd_status
        ;;
    ingest)
        shift
        cmd_ingest "$@"
        ;;
    chat)
        shift
        cmd_chat "$@"
        ;;
    portal)
        cmd_portal
        ;;
    demo-setup)
        cmd_demo_setup
        ;;
    all)
        cmd_all
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        print_error "Unknown command: $1"
        show_help
        exit 1
        ;;
esac
