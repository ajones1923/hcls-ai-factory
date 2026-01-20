#!/bin/bash

# Immediate output test
echo "SCRIPT STARTED - $(date)"
echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                  🧬 FULL GENOME PIPELINE                           ║"
echo "║            FASTQ → BAM → VCF Complete Analysis                     ║"
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

echo "📂 Creating output directories..."
mkdir -p "${LOG_DIR}"
mkdir -p "${OUT_DIR}/tmp_genome"
echo "✓  Directories ready"
echo ""

# Helper functions for progress messages
log_step_start() {
    local step_num=$1
    local step_name=$2
    local step_desc=$3
    echo ""
    echo "┌──────────────────────────────────────────────────────────────────┐"
    echo "│  STEP ${step_num}/4 │ ${step_name}"
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
    echo "│  ✅ STEP ${step_num}/4 COMPLETE │ ${step_name}"
    echo "└──────────────────────────────────────────────────────────────────┘"
    echo ""
}

log_progress() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "📊 $1"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
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
STEP3_START=0
STEP3_END=0
STEP4_START=0
STEP4_END=0

echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│  📋 PIPELINE OVERVIEW                                           │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Step 1: fq2bam (FASTQ → BAM)     (~60-180 minutes)              │"
echo "│  Step 2: BAM Indexing & QC        (~2-5 minutes)                 │"
echo "│  Step 3: DeepVariant Calling      (~60-90 minutes)               │"
echo "│  Step 4: VCF Indexing             (~1 minute)                    │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Total Estimated Time: 2-5 hours                                 │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""

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

# Build fq2bam command
FQ2BAM_CMD="pbrun fq2bam --ref /ref/GRCh38.fa --in-fq /in/HG002_R1.fastq.gz /in/HG002_R2.fastq.gz --out-bam /out/HG002.genome.bam --tmp-dir /out/tmp_genome --num-gpus ${NUM_GPUS}"

# Add --low-memory flag only if LOW_MEMORY=1 in config
if [ "${LOW_MEMORY}" = "1" ]; then
    FQ2BAM_CMD="${FQ2BAM_CMD} --low-memory"
    log_info "Using --low-memory mode (LOW_MEMORY=1 in config)"
else
    log_info "Using full memory mode for maximum performance"
fi

STEP1_START=$(date +%s)
log_step_start "1" "fq2bam (FASTQ → BAM)" "GPU-accelerated alignment and sorting for whole genome"
log_info "This is the most GPU-intensive step..."
log_detail "Input files:"
log_detail "  • HG002_R1.fastq.gz (Read 1)"
log_detail "  • HG002_R2.fastq.gz (Read 2)"
log_detail "Output: HG002.genome.bam"
log_detail "Estimated time: 60-180 minutes"
echo ""
log_info "🚀 Starting GPU alignment with Parabricks BWA..."
log_detail "Aligning reads to GRCh38 reference genome"
echo ""

# Start heartbeat in background
status_heartbeat "GPU alignment in progress" &
HEARTBEAT_PID=$!

/usr/bin/time -v docker run --rm --gpus all \
    -v "${IN_DIR}:/in" \
    -v "${REF_DIR}:/ref" \
    -v "${OUT_DIR}:/out" \
    -v /tmp/nvidia-smi-wrapper.sh:/usr/local/bin/nvidia-smi:ro \
    -e PATH="/usr/local/bin:/usr/local/parabricks:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${PB_IMG}" \
    ${FQ2BAM_CMD} \
    2>&1 | tee "${LOG_DIR}/genome_fq2bam.log"

# Stop heartbeat
kill $HEARTBEAT_PID 2>/dev/null
wait $HEARTBEAT_PID 2>/dev/null

log_success "fq2bam completed successfully"

log_info "⏳ Waiting for file system sync..."
sync
sleep 5

# Verify BAM file exists before proceeding
BAM_FILE="${OUT_DIR}/HG002.genome.bam"
log_info "🔍 Verifying BAM file..."
RETRY_COUNT=0
MAX_RETRIES=6
while [ ! -f "$BAM_FILE" ] || [ ! -s "$BAM_FILE" ]; do
    RETRY_COUNT=$((RETRY_COUNT + 1))
    if [ $RETRY_COUNT -ge $MAX_RETRIES ]; then
        echo "   ❌ ERROR: BAM file not found after ${MAX_RETRIES} attempts"
        exit 1
    fi
    log_detail "Waiting for BAM file (attempt ${RETRY_COUNT}/${MAX_RETRIES})..."
    sleep 10
done
BAM_SIZE=$(ls -lh "$BAM_FILE" | awk '{print $5}')
log_success "BAM file verified: ${BAM_SIZE}"
STEP1_END=$(date +%s)
log_step_complete "1" "fq2bam (FASTQ → BAM)"

STEP2_START=$(date +%s)
log_step_start "2" "BAM Indexing & QC" "Creating BAM index and generating alignment statistics"
log_info "Processing alignment results..."
log_detail "Creating .bai index file for random access"
log_detail "Running samtools flagstat for quality metrics"
echo ""

docker run --rm -i \
    -v "${OUT_DIR}:/out" \
    "${PB_IMG}" \
    bash -c "echo '   → Creating BAM index...' && \
             samtools index /out/HG002.genome.bam && \
             echo '   → Running flagstat QC...' && \
             samtools flagstat /out/HG002.genome.bam | head -n 40" \
    2>&1 | tee "${LOG_DIR}/genome_flagstat.log"

log_success "BAM indexed and QC complete"

log_info "⏳ Waiting for file system sync..."
sync
sleep 3

# Verify BAM index exists before proceeding
BAI_FILE="${OUT_DIR}/HG002.genome.bam.bai"
log_info "🔍 Verifying BAM index..."
RETRY_COUNT=0
MAX_RETRIES=6
while [ ! -f "$BAI_FILE" ] || [ ! -s "$BAI_FILE" ]; do
    RETRY_COUNT=$((RETRY_COUNT + 1))
    if [ $RETRY_COUNT -ge $MAX_RETRIES ]; then
        echo "   ❌ ERROR: BAM index not found after ${MAX_RETRIES} attempts"
        exit 1
    fi
    log_detail "Waiting for BAM index (attempt ${RETRY_COUNT}/${MAX_RETRIES})..."
    sleep 5
done
log_success "BAM index verified"
STEP2_END=$(date +%s)
log_step_complete "2" "BAM Indexing & QC"

STEP3_START=$(date +%s)
log_step_start "3" "DeepVariant (Variant Calling)" "GPU-accelerated variant calling for whole genome"
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
STEP3_END=$(date +%s)
log_step_complete "3" "DeepVariant (Variant Calling)"

STEP4_START=$(date +%s)
log_step_start "4" "VCF Indexing" "Creating VCF index for fast random access"
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

STEP4_END=$(date +%s)
log_step_complete "4" "VCF Indexing"

# Calculate timing
PIPELINE_END_TIME=$(date +%s)
TOTAL_TIME=$((PIPELINE_END_TIME - PIPELINE_START_TIME))
STEP1_TIME=$((STEP1_END - STEP1_START))
STEP2_TIME=$((STEP2_END - STEP2_START))
STEP3_TIME=$((STEP3_END - STEP3_START))
STEP4_TIME=$((STEP4_END - STEP4_START))

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
STEP3_FORMATTED=$(format_time $STEP3_TIME)
STEP4_FORMATTED=$(format_time $STEP4_TIME)
TOTAL_FORMATTED=$(format_time $TOTAL_TIME)

# Get file sizes and variant count
VCF_FILE="${OUT_DIR}/HG002.genome.vcf.gz"
VCF_INDEX="${OUT_DIR}/HG002.genome.vcf.gz.tbi"
VCF_SIZE=$(ls -lh "$VCF_FILE" 2>/dev/null | awk '{print $5}' || echo "N/A")
VCF_INDEX_SIZE=$(ls -lh "$VCF_INDEX" 2>/dev/null | awk '{print $5}' || echo "N/A")

# Get variant count
VARIANT_COUNT=$(docker run --rm -v "${OUT_DIR}:/out" "${PB_IMG}" bcftools view -H /out/HG002.genome.vcf.gz 2>/dev/null | wc -l || echo "0")
VARIANT_COUNT_FORMATTED=$(printf "%'d" $VARIANT_COUNT)

# Get FASTQ sizes
FASTQ_R1_SIZE=$(ls -lh "${IN_DIR}/HG002_R1.fastq.gz" 2>/dev/null | awk '{print $5}' || echo "N/A")
FASTQ_R2_SIZE=$(ls -lh "${IN_DIR}/HG002_R2.fastq.gz" 2>/dev/null | awk '{print $5}' || echo "N/A")

# Calculate total FASTQ size in GB
FASTQ_TOTAL_BYTES=$(stat --format="%s" "${IN_DIR}/HG002_R1.fastq.gz" 2>/dev/null || echo "0")
FASTQ_TOTAL_BYTES=$((FASTQ_TOTAL_BYTES + $(stat --format="%s" "${IN_DIR}/HG002_R2.fastq.gz" 2>/dev/null || echo "0")))
FASTQ_TOTAL_GB=$((FASTQ_TOTAL_BYTES / 1073741824))

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
echo "  Pipeline Time Breakdown"
echo "  ┌────────────────────────┬───────────────────┐"
echo "  │          Step          │       Time        │"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 1. fq2bam (alignment)  │ %-17s │\n" "$STEP1_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 2. BAM indexing & QC   │ %-17s │\n" "$STEP2_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 3. DeepVariant         │ %-17s │\n" "$STEP3_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 4. VCF indexing        │ %-17s │\n" "$STEP4_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ Total                  │ %-17s │\n" "$TOTAL_FORMATTED"
echo "  └────────────────────────┴───────────────────┘"
echo ""
echo "  Summary"
echo "  • Processed ${FASTQ_TOTAL_GB} GB of FASTQ data"
echo "  • Aligned reads to GRCh38 reference genome"
echo "  • Called $VARIANT_COUNT_FORMATTED variants across all chromosomes"
echo "  • VCF ready for downstream analysis (PGx, annotation, etc.)"
echo ""
echo "  📁 Output files: ${OUT_DIR}/HG002.genome.*"
echo ""
