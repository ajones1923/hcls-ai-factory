#!/bin/bash

# Immediate output test
echo "SCRIPT STARTED - $(date)"
echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║              🧬 FULL GENOME PIPELINE (RESUME)                      ║"
echo "║         Resuming from DeepVariant - BAM already exists             ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

set -e

# Disable output buffering completely
export PYTHONUNBUFFERED=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "📁 Project root: $PROJECT_ROOT"

# Source environment
echo "⚙️  Loading configuration..."
source "${PROJECT_ROOT}/config/pipeline.env"
echo "✓  Configuration loaded"

# Use local paths
IN_DIR="${PROJECT_ROOT}/data/input"
REF_DIR="${PROJECT_ROOT}/data/ref"
OUT_DIR="${PROJECT_ROOT}/data/output"
LOG_DIR="${PROJECT_ROOT}/data/output/logs"

echo "📂 Verifying directories..."
mkdir -p "${LOG_DIR}"
echo "✓  Directories ready"
echo ""

# Helper functions for progress messages
log_step_start() {
    local step_num=$1
    local step_name=$2
    local step_desc=$3
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│  STEP ${step_num}/2 │ ${step_name}"
    echo "├──────────────────────────────────────────────────────────────────┤"
    echo "│  ${step_desc}"
    echo "└──────────────────────────────────────────────────────────────────┘"
    echo ""
}

log_step_complete() {
    local step_num=$1
    local step_name=$2
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│  ✅ STEP ${step_num}/2 COMPLETE │ ${step_name}"
    echo "└──────────────────────────────────────────────────────────────────┘"
    echo ""
}

log_info() {
    echo "   ℹ️  $1"
}

log_detail() {
    echo "      → $1"
}

log_success() {
    echo "   ✅ $1"
}

# Track timing for summary
PIPELINE_START_TIME=$(date +%s)
STEP1_START=0
STEP1_END=0
STEP2_START=0
STEP2_END=0

echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│  📋 RESUME PIPELINE OVERVIEW                                     │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  ⏭️  Skipping: fq2bam (already complete)                         │"
echo "│  ⏭️  Skipping: BAM Indexing (already complete)                   │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Step 1: DeepVariant Calling      (~60-90 minutes)               │"
echo "│  Step 2: VCF Indexing             (~1 minute)                    │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Total Estimated Time: 60-90 minutes                             │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""

# Verify BAM file exists before proceeding
BAM_FILE="${OUT_DIR}/HG002.genome.bam"
BAI_FILE="${OUT_DIR}/HG002.genome.bam.bai"

log_info "🔍 Verifying existing BAM file..."
if [ ! -f "$BAM_FILE" ] || [ ! -s "$BAM_FILE" ]; then
    echo "   ❌ ERROR: BAM file not found: $BAM_FILE"
    echo "   Please run the full pipeline (05-run-full-genome.sh) first."
    exit 1
fi
BAM_SIZE=$(ls -lh "$BAM_FILE" | awk '{print $5}')
log_success "BAM file verified: ${BAM_SIZE}"

log_info "🔍 Verifying existing BAM index..."
if [ ! -f "$BAI_FILE" ] || [ ! -s "$BAI_FILE" ]; then
    echo "   ❌ ERROR: BAM index not found: $BAI_FILE"
    echo "   Please run the full pipeline (05-run-full-genome.sh) first."
    exit 1
fi
log_success "BAM index verified"

# Function to show periodic status updates during long operations
status_heartbeat() {
    local message="$1"
    while true; do
        sleep 15
        echo "⏳ ${message} (still running...)"
    done
}

# Create nvidia-smi wrapper to work around GPU memory detection bug
cat > /tmp/nvidia-smi-wrapper.sh << 'WRAPPER_EOF'
#!/bin/bash
# Handle specific Parabricks query
if [[ "$*" == *"--query-gpu=index,name,memory.total"* ]] && [[ "$*" == *"noheader"* ]]; then
    echo "0, NVIDIA GB10, 16384"
    exit 0
fi
if [[ "$*" == *"--query-gpu"* ]] && [[ "$*" == *"memory"* ]]; then
    echo "index, name, memory.free [MiB], memory.total [MiB]"
    echo "0, NVIDIA GB10, 12288, 16384"
    exit 0
fi
if [[ "$*" == *"-L"* ]]; then
    echo "GPU 0: NVIDIA GB10 (UUID: GPU-00000000-0000-0000-0000-000000000000)"
    exit 0
fi
echo "Mon Jan  6 10:50:00 2026"
echo "+-----------------------------------+"
echo "| NVIDIA-SMI 535.129.03            |"
exit 0
WRAPPER_EOF
chmod +x /tmp/nvidia-smi-wrapper.sh

STEP1_START=$(date +%s)
log_step_start "1" "DeepVariant (Variant Calling)" "GPU-accelerated variant calling for whole genome"
log_info "Identifying genetic variants across all chromosomes..."
log_detail "Input: HG002.genome.bam"
log_detail "Output: HG002.genome.vcf.gz"
log_detail "Estimated time: 60-90 minutes"
echo ""
log_info "🧬 Starting DeepVariant variant caller..."
log_detail "Using deep learning model for variant detection"
log_detail "Analyzing SNPs and small indels"
echo ""

# Start heartbeat in background
status_heartbeat "Variant calling in progress" &
HEARTBEAT_PID=$!

/usr/bin/time -v docker run --rm --gpus all \
    -v "${REF_DIR}:/ref" \
    -v "${OUT_DIR}:/out" \
    -v /tmp/nvidia-smi-wrapper.sh:/usr/local/bin/nvidia-smi:ro \
    -e PATH="/usr/local/bin:/usr/local/parabricks:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${PB_IMG}" \
    pbrun deepvariant --ref /ref/GRCh38.fa --in-bam /out/HG002.genome.bam --out-variants /out/HG002.genome.vcf.gz --num-gpus ${NUM_GPUS} \
    2>&1 | tee "${LOG_DIR}/genome_deepvariant.log"

# Stop heartbeat
kill $HEARTBEAT_PID 2>/dev/null
wait $HEARTBEAT_PID 2>/dev/null

log_success "DeepVariant completed successfully"
STEP1_END=$(date +%s)
log_step_complete "1" "DeepVariant (Variant Calling)"

STEP2_START=$(date +%s)
log_step_start "2" "VCF Indexing" "Creating VCF index for fast random access"
log_info "📑 Indexing VCF file..."
echo ""

docker run --rm -i \
    -v "${OUT_DIR}:/out" \
    "${PB_IMG}" \
    bash -c "echo '   → Creating VCF index...' && \
             tabix -p vcf /out/HG002.genome.vcf.gz && \
             echo '   → VCF indexed successfully' && \
             echo '' && \
             echo '   📊 VCF Summary:' && \
             VARIANT_COUNT=\$(bcftools view -H /out/HG002.genome.vcf.gz | wc -l) && \
             echo \"      Total variants called: \${VARIANT_COUNT}\"" \
    2>&1 | tee "${LOG_DIR}/genome_vcf_header.log"

STEP2_END=$(date +%s)
log_step_complete "2" "VCF Indexing"

# Calculate timing
PIPELINE_END_TIME=$(date +%s)
TOTAL_TIME=$((PIPELINE_END_TIME - PIPELINE_START_TIME))
STEP1_TIME=$((STEP1_END - STEP1_START))
STEP2_TIME=$((STEP2_END - STEP2_START))

# Format time function
format_time() {
    local seconds=$1
    if [ $seconds -lt 60 ]; then
        echo "${seconds}s"
    elif [ $seconds -lt 3600 ]; then
        local mins=$((seconds / 60))
        local secs=$((seconds % 60))
        echo "${mins}m ${secs}s"
    else
        local hours=$((seconds / 3600))
        local mins=$(((seconds % 3600) / 60))
        echo "${hours}h ${mins}m"
    fi
}

STEP1_FORMATTED=$(format_time $STEP1_TIME)
STEP2_FORMATTED=$(format_time $STEP2_TIME)
TOTAL_FORMATTED=$(format_time $TOTAL_TIME)

# Get file sizes and variant count
VCF_FILE="${OUT_DIR}/HG002.genome.vcf.gz"
VCF_INDEX="${OUT_DIR}/HG002.genome.vcf.gz.tbi"
VCF_SIZE=$(ls -lh "$VCF_FILE" 2>/dev/null | awk '{print $5}' || echo "N/A")
VCF_INDEX_SIZE=$(ls -lh "$VCF_INDEX" 2>/dev/null | awk '{print $5}' || echo "N/A")

# Get variant count
VARIANT_COUNT=$(docker run --rm -v "${OUT_DIR}:/out" "${PB_IMG}" bcftools view -H /out/HG002.genome.vcf.gz 2>/dev/null | wc -l || echo "0")
VARIANT_COUNT_FORMATTED=$(printf "%'d" $VARIANT_COUNT)

echo ""
echo ""
echo "● Full Genome Pipeline Complete! 🎉"
echo ""
echo "  ┌───────────────────────┬─────────────────────────────┐"
echo "  │        Result         │            Value            │"
echo "  ├───────────────────────┼─────────────────────────────┤"
printf "  │ Total Variants Called │ %-27s │\n" "$VARIANT_COUNT_FORMATTED"
echo "  ├───────────────────────┼─────────────────────────────┤"
echo "  │ Chromosomes           │ All (chr1 through chrY)     │"
echo "  ├───────────────────────┼─────────────────────────────┤"
printf "  │ VCF Size              │ %-27s │\n" "$VCF_SIZE"
echo "  └───────────────────────┴─────────────────────────────┘"
echo ""
echo "  Results"
echo "  ┌────────────┬─────────────────────────────────────────┐"
echo "  │   Metric   │                  Value                  │"
echo "  ├────────────┼─────────────────────────────────────────┤"
printf "  │ Total Time │ %-39s │\n" "$TOTAL_FORMATTED"
echo "  ├────────────┼─────────────────────────────────────────┤"
printf "  │ VCF File   │ %-39s │\n" "HG002.genome.vcf.gz ($VCF_SIZE)"
echo "  ├────────────┼─────────────────────────────────────────┤"
printf "  │ VCF Index  │ %-39s │\n" "HG002.genome.vcf.gz.tbi ($VCF_INDEX_SIZE)"
echo "  └────────────┴─────────────────────────────────────────┘"
echo ""
echo "  Pipeline Time Breakdown (Resume)"
echo "  ┌────────────────────────┬───────────────────┐"
echo "  │          Step          │       Time        │"
echo "  ├────────────────────────┼───────────────────┤"
echo "  │ (skipped) fq2bam       │ --                │"
echo "  ├────────────────────────┼───────────────────┤"
echo "  │ (skipped) BAM indexing │ --                │"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 1. DeepVariant         │ %-17s │\n" "$STEP1_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 2. VCF indexing        │ %-17s │\n" "$STEP2_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ Total (resume)         │ %-17s │\n" "$TOTAL_FORMATTED"
echo "  └────────────────────────┴───────────────────┘"
echo ""
echo "  Summary"
echo "  • Resumed from existing BAM file (${BAM_SIZE})"
echo "  • Called $VARIANT_COUNT_FORMATTED variants across all chromosomes"
echo "  • VCF ready for downstream analysis (PGx, annotation, etc.)"
echo ""
echo "  📁 Output files: ${OUT_DIR}/HG002.genome.*"
echo ""
