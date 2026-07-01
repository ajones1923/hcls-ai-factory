#!/bin/bash

# Immediate output test
echo "SCRIPT STARTED - $(date)"
echo "Setting up GRCh38 reference genome..."

set -e

# Disable output buffering completely
export PYTHONUNBUFFERED=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "Project root: $PROJECT_ROOT"

# Source environment
echo "Loading configuration..."
source "${PROJECT_ROOT}/config/pipeline.env"
echo "Configuration loaded"

# Use local paths
REF_DIR="${PROJECT_ROOT}/data/ref"
echo "Reference directory: $REF_DIR"

# Helper functions for progress messages
log_progress() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "📊 $1"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
}

log_info() {
    echo "ℹ️  $1"
}

log_success() {
    echo "✅ $1"
}

# Function to show periodic status updates during long operations
status_heartbeat() {
    local message="$1"
    while true; do
        sleep 15
        echo "⏳ ${message} (still running...)"
    done
}

# Function to run command with periodic status updates
run_with_status() {
    local status_message="$1"
    shift

    # Start heartbeat in background
    status_heartbeat "$status_message" &
    local heartbeat_pid=$!

    # Run the actual command
    "$@"
    local exit_code=$?

    # Stop heartbeat
    kill $heartbeat_pid 2>/dev/null
    wait $heartbeat_pid 2>/dev/null

    return $exit_code
}

echo "=========================================="
echo "📚 GRCh38 Reference Genome Setup"
echo "=========================================="
echo ""
log_info "This script will download and prepare the GRCh38 human reference genome"
log_info "Source: NVIDIA Parabricks sample bundle (~11GB download)"
log_info "GPU acceleration will be used for indexing"
echo ""

echo "Creating reference directory..."
mkdir -p "${REF_DIR}"
echo "Changing to reference directory..."
cd "${REF_DIR}"
echo "Current directory: $(pwd)"

log_progress "STEP 1/4: Downloading Parabricks Sample Bundle"
if [ ! -f "parabricks_sample.tar.gz" ]; then
    log_info "Downloading from: ${PARABRICKS_SAMPLE_URL}"
    log_info "Size: ~11GB - this will take 5-20 minutes depending on connection"
    echo ""

    # Start heartbeat in background
    status_heartbeat "Downloading Parabricks sample bundle" &
    HEARTBEAT_PID=$!

    # Run download
    wget --progress=dot:giga -O parabricks_sample.tar.gz "${PARABRICKS_SAMPLE_URL}" 2>&1 | \
        grep --line-buffered -E '^\s+[0-9]+K' || true

    # Stop heartbeat
    kill $HEARTBEAT_PID 2>/dev/null
    wait $HEARTBEAT_PID 2>/dev/null

    echo ""
    log_success "Download complete"
else
    log_success "Already downloaded (parabricks_sample.tar.gz exists)"
    log_info "Skipping download step"
fi

log_progress "STEP 2/4: Extracting Reference Genome"
if [ ! -f "GRCh38.fa" ]; then
    log_info "Extracting parabricks_sample.tar.gz..."
    run_with_status "Extracting reference genome" tar -xzf parabricks_sample.tar.gz
    log_info "Copying Homo_sapiens_assembly38.fasta to GRCh38.fa..."
    cp parabricks_sample/Ref/Homo_sapiens_assembly38.fasta GRCh38.fa
    FASTA_SIZE=$(ls -lh GRCh38.fa | awk '{print $5}')
    log_success "Reference genome extracted: GRCh38.fa (${FASTA_SIZE})"
else
    log_success "Already extracted (GRCh38.fa exists)"
    log_info "Skipping extraction step"
fi

log_progress "STEP 3/5: Building FASTA Index"
log_info "Using NVIDIA Parabricks container"
log_info "Creating FASTA index (.fai) using samtools"
echo ""

if [ ! -f "${REF_DIR}/GRCh38.fa.fai" ]; then
    echo "Running samtools faidx..."
    run_with_status "Creating FASTA index" \
        docker run --rm -i \
            -v "${REF_DIR}:/ref" \
            "${PB_IMG}" \
            bash -c "samtools faidx /ref/GRCh38.fa" 2>&1
    echo "FASTA index created"
    echo ""
    log_success "FASTA index created"
else
    log_success "FASTA index already exists"
    log_info "Skipping FASTA indexing step"
fi

log_progress "STEP 4/5: Copying BWA Index Files"
log_info "Using pre-built BWA index files from Parabricks sample bundle"
log_info "BWA index is required for Parabricks fq2bam alignment"
log_info "Index files: .bwt (3.0GB), .pac (768MB), .sa (1.5GB), .amb, .ann"
echo ""

if [ ! -f "${REF_DIR}/GRCh38.fa.bwt" ]; then
    echo "Copying BWA index files..."
    cp "${REF_DIR}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.amb" "${REF_DIR}/GRCh38.fa.amb"
    cp "${REF_DIR}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.ann" "${REF_DIR}/GRCh38.fa.ann"
    echo "Copying .bwt file (3.0GB)..."
    cp "${REF_DIR}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.bwt" "${REF_DIR}/GRCh38.fa.bwt"
    echo "Copying .pac file (768MB)..."
    cp "${REF_DIR}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.pac" "${REF_DIR}/GRCh38.fa.pac"
    echo "Copying .sa file (1.5GB)..."
    cp "${REF_DIR}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.sa" "${REF_DIR}/GRCh38.fa.sa"
    echo ""
    log_success "BWA index files copied successfully"
else
    log_success "BWA index already exists"
    log_info "Skipping BWA index copy step"
fi

log_progress "STEP 5/5: Creating Sequence Dictionary"
log_info "Using samtools dict"
log_info "This creates a .dict file needed for variant calling"
echo ""

if [ ! -f "${REF_DIR}/GRCh38.dict" ]; then
    echo "Running samtools dict..."
    run_with_status "Creating sequence dictionary" \
        docker run --rm -i \
            -v "${REF_DIR}:/ref" \
            "${PB_IMG}" \
            bash -c "samtools dict /ref/GRCh38.fa -o /ref/GRCh38.dict" 2>&1
    echo "Sequence dictionary created"
    echo ""
    log_success "Sequence dictionary created"
else
    log_success "Sequence dictionary already exists"
    log_info "Skipping dictionary creation step"
fi

echo ""
echo "=========================================="
log_success "✅ REFERENCE GENOME SETUP COMPLETE!"
echo "=========================================="
echo ""
log_info "Verifying reference files..."
echo ""
echo "📁 Reference Files Created:"
ls -lh "${REF_DIR}"/GRCh38.fa* "${REF_DIR}"/GRCh38.dict 2>/dev/null | awk '{printf "   %s  %s\n", $5, $9}' || true
echo ""
log_info "Reference genome is ready for pipeline processing"
echo ""
log_info "Next steps:"
log_info "  - Quick test: ./scripts/04-run-chr20-test.sh (5-20 minutes)"
log_info "  - Full genome: ./scripts/05-run-full-genome.sh (30-110 minutes)"
echo ""
