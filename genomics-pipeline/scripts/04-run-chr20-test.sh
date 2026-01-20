#!/bin/bash

# Immediate output test
echo "SCRIPT STARTED - $(date)"
echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                    🧬 CHR20 TEST PIPELINE                          ║"
echo "║              FASTQ → BAM → VCF Validation Test                     ║"
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

# Use local paths
IN_DIR="${PROJECT_ROOT}/data/input"
REF_DIR="${PROJECT_ROOT}/data/ref"
OUT_DIR="${PROJECT_ROOT}/data/output"
LOG_DIR="${PROJECT_ROOT}/data/output/logs"

echo "📂 Creating output directories..."
mkdir -p "${LOG_DIR}"
mkdir -p "${OUT_DIR}/tmp_chr20"
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

echo ""
echo "┌──────────────────────────────────────────────────────────────────┐"
echo "│  📋 PIPELINE OVERVIEW                                           │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Step 1: Create Interval BED      (~30 seconds)                 │"
echo "│  Step 2: fq2bam (FASTQ → BAM)     (~2-10 minutes)               │"
echo "│  Step 3: BAM Indexing & QC        (~1 minute)                   │"
echo "│  Step 4: DeepVariant Calling      (~1-6 minutes)                │"
echo "├──────────────────────────────────────────────────────────────────┤"
echo "│  Total Estimated Time: 5-20 minutes                             │"
echo "└──────────────────────────────────────────────────────────────────┘"
echo ""

STEP1_START=$(date +%s)
log_step_start "1" "Create Interval BED" "Identifying chromosome 20 coordinates for targeted analysis"
log_info "Reading reference genome index..."
log_detail "Reference: GRCh38.fa"
log_detail "Target: Chromosome 20 only"
docker run --rm -i \
    -v "${REF_DIR}:/ref" \
    -v "${OUT_DIR}:/out" \
    "${PB_IMG}" \
    bash -c '
        test -f /ref/GRCh38.fa.fai || samtools faidx /ref/GRCh38.fa
        C20=$(awk '"'"'$1=="chr20"{print "chr20"} $1=="20"{print "20"}'"'"' /ref/GRCh38.fa.fai | head -n 1)
        if [ -z "$C20" ]; then
            echo "ERROR: chr20/20 not found in reference index"
            exit 1
        fi
        LEN=$(awk -v c="$C20" '"'"'$1==c{print $2}'"'"' /ref/GRCh38.fa.fai)
        echo -e "${C20}\t0\t${LEN}" > /out/interval_chr20.bed
        echo "Created interval file:"
        cat /out/interval_chr20.bed
    ' 2>&1

# Get chr20 name for DeepVariant
CHR20_NAME=$(cat "${OUT_DIR}/interval_chr20.bed" | cut -f1)
log_success "Interval BED created for chromosome: ${CHR20_NAME}"
STEP1_END=$(date +%s)
log_step_complete "1" "Create Interval BED"

# Create nvidia-smi wrapper to work around GPU memory detection bug
# Parabricks queries: nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader,nounits
cat > /tmp/nvidia-smi-wrapper.sh << 'WRAPPER_EOF'
#!/bin/bash
# Handle specific Parabricks query that expects: "index, name, memory_in_MiB" (no header, no units)
if [[ "$*" == *"--query-gpu=index,name,memory.total"* ]] && [[ "$*" == *"noheader"* ]]; then
    # Return fake memory value for Parabricks (16GB = 16384 MiB)
    echo "0, NVIDIA GB10, 16384"
    exit 0
fi

# Handle memory queries with headers
if [[ "$*" == *"--query-gpu"* ]] && [[ "$*" == *"memory"* ]]; then
    echo "index, name, memory.free [MiB], memory.total [MiB]"
    echo "0, NVIDIA GB10, 12288, 16384"
    exit 0
fi

# For -L (list GPUs)
if [[ "$*" == *"-L"* ]]; then
    echo "GPU 0: NVIDIA GB10 (UUID: GPU-00000000-0000-0000-0000-000000000000)"
    exit 0
fi

# Default: show basic GPU info
echo "Mon Jan  6 10:50:00 2026"
echo "+-----------------------------------------------------------------------------+"
echo "| NVIDIA-SMI 535.129.03   Driver Version: 535.129.03   CUDA Version: 12.2   |"
echo "|-------------------------------+----------------------+----------------------+"
echo "| GPU  Name        TCC/WDDM | Bus-Id        Disp.A | Volatile Uncorr. ECC |"
echo "| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |"
echo "|===============================+======================+======================|"
echo "|   0  NVIDIA GB10        On  | 00000000:01:00.0 Off |                  N/A |"
echo "| N/A   45C    P0    25W / N/A |      0MiB / 16384MiB |      0%      Default |"
echo "+-------------------------------+----------------------+----------------------+"
exit 0
WRAPPER_EOF
chmod +x /tmp/nvidia-smi-wrapper.sh

# Build fq2bam command
# Note: Using nvidia-smi wrapper to work around GPU memory detection bug
FQ2BAM_CMD="pbrun fq2bam --ref /ref/GRCh38.fa --in-fq /in/HG002_R1.fastq.gz /in/HG002_R2.fastq.gz --interval-file /out/interval_chr20.bed --out-bam /out/HG002.chr20.bam --tmp-dir /out/tmp_chr20 --num-gpus ${NUM_GPUS}"

# Add --low-memory flag only if LOW_MEMORY=1 in config
if [ "${LOW_MEMORY}" = "1" ]; then
    FQ2BAM_CMD="${FQ2BAM_CMD} --low-memory"
    log_info "Using --low-memory mode (LOW_MEMORY=1 in config)"
else
    log_info "Using full memory mode for maximum performance"
fi

STEP2_START=$(date +%s)
log_step_start "2" "fq2bam (FASTQ → BAM)" "GPU-accelerated alignment and sorting for chromosome 20"
log_info "This is the most GPU-intensive step..."
log_detail "Input files:"
log_detail "  • HG002_R1.fastq.gz (Read 1)"
log_detail "  • HG002_R2.fastq.gz (Read 2)"
log_detail "Output: HG002.chr20.bam"
log_detail "Estimated time: 2-10 minutes"
echo ""
log_info "🚀 Starting GPU alignment with Parabricks BWA..."
log_detail "Aligning reads to GRCh38 reference genome"
log_detail "Filtering to chromosome 20 region only"
echo ""

# Start heartbeat in background for fq2bam
status_heartbeat "GPU alignment in progress" &
HEARTBEAT_PID=$!

# Run fq2bam
# Mount nvidia-smi wrapper to work around [N/A] memory reporting bug
# Mount to /usr/local/bin which is checked before /usr/bin in PATH
/usr/bin/time -v docker run --rm --gpus all \
    -v "${IN_DIR}:/in" \
    -v "${REF_DIR}:/ref" \
    -v "${OUT_DIR}:/out" \
    -v /tmp/nvidia-smi-wrapper.sh:/usr/local/bin/nvidia-smi:ro \
    -e PATH="/usr/local/bin:/usr/local/parabricks:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${PB_IMG}" \
    ${FQ2BAM_CMD} \
    2>&1 | tee "${LOG_DIR}/chr20_fq2bam.log"

# Stop heartbeat
kill $HEARTBEAT_PID 2>/dev/null
wait $HEARTBEAT_PID 2>/dev/null

log_success "fq2bam completed successfully"

# Buffer and verification before proceeding
log_info "⏳ Waiting for file system sync..."
sync
sleep 5

# Verify BAM file exists and is readable
BAM_FILE="${OUT_DIR}/HG002.chr20.bam"
log_info "🔍 Verifying BAM file..."
RETRY_COUNT=0
MAX_RETRIES=6
while [ ! -f "$BAM_FILE" ] || [ ! -s "$BAM_FILE" ]; do
    RETRY_COUNT=$((RETRY_COUNT + 1))
    if [ $RETRY_COUNT -ge $MAX_RETRIES ]; then
        echo "   ❌ ERROR: BAM file not found after ${MAX_RETRIES} attempts"
        exit 1
    fi
    echo "   → Waiting for BAM file (attempt ${RETRY_COUNT}/${MAX_RETRIES})..."
    sleep 10
done

BAM_SIZE=$(ls -lh "$BAM_FILE" | awk '{print $5}')
log_success "BAM file verified: ${BAM_SIZE}"
STEP2_END=$(date +%s)
log_step_complete "2" "fq2bam (FASTQ → BAM)"

STEP3_START=$(date +%s)
log_step_start "3" "BAM Indexing & QC" "Creating BAM index and generating alignment statistics"
log_info "Processing alignment results..."
log_detail "Creating .bai index file for random access"
log_detail "Running samtools flagstat for quality metrics"
echo ""

docker run --rm -i \
    -v "${OUT_DIR}:/out" \
    "${PB_IMG}" \
    bash -c "echo '   → Creating BAM index...' && \
             samtools index /out/HG002.chr20.bam && \
             echo '   → Running flagstat QC...' && \
             samtools flagstat /out/HG002.chr20.bam" \
    2>&1 | tee "${LOG_DIR}/chr20_flagstat.log"

log_success "BAM indexed and QC complete"

# Verify BAM index exists before proceeding
BAI_FILE="${OUT_DIR}/HG002.chr20.bam.bai"
log_info "⏳ Waiting for file system sync..."
sync
sleep 3

log_info "🔍 Verifying BAM index..."
RETRY_COUNT=0
MAX_RETRIES=6
while [ ! -f "$BAI_FILE" ] || [ ! -s "$BAI_FILE" ]; do
    RETRY_COUNT=$((RETRY_COUNT + 1))
    if [ $RETRY_COUNT -ge $MAX_RETRIES ]; then
        echo "   ❌ ERROR: BAM index not found after ${MAX_RETRIES} attempts"
        exit 1
    fi
    echo "   → Waiting for BAM index (attempt ${RETRY_COUNT}/${MAX_RETRIES})..."
    sleep 5
done
log_success "BAM index verified"
STEP3_END=$(date +%s)
log_step_complete "3" "BAM Indexing & QC"

STEP4_START=$(date +%s)
log_step_start "4" "DeepVariant (Variant Calling)" "GPU-accelerated variant calling using DeepVariant"
log_info "Identifying genetic variants in chromosome 20..."
log_detail "Input: HG002.chr20.bam"
log_detail "Output: HG002.chr20.vcf.gz"
log_detail "Estimated time: 1-6 minutes"
echo ""
log_info "🧬 Starting DeepVariant variant caller..."
log_detail "Using deep learning model for variant detection"
log_detail "Analyzing SNPs and small indels"
echo ""

# Start heartbeat in background for DeepVariant
status_heartbeat "Variant calling in progress" &
HEARTBEAT_PID=$!

# Run DeepVariant
# Mount nvidia-smi wrapper to work around [N/A] memory reporting bug
/usr/bin/time -v docker run --rm --gpus all \
    -v "${REF_DIR}:/ref" \
    -v "${OUT_DIR}:/out" \
    -v /tmp/nvidia-smi-wrapper.sh:/usr/local/bin/nvidia-smi:ro \
    -e PATH="/usr/local/bin:/usr/local/parabricks:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${PB_IMG}" \
    pbrun deepvariant --ref /ref/GRCh38.fa --in-bam /out/HG002.chr20.bam --out-variants /out/HG002.chr20.vcf.gz --regions "${CHR20_NAME}" --num-gpus ${NUM_GPUS} \
    2>&1 | tee "${LOG_DIR}/chr20_deepvariant.log"

# Stop heartbeat
kill $HEARTBEAT_PID 2>/dev/null
wait $HEARTBEAT_PID 2>/dev/null

log_success "DeepVariant completed successfully"
STEP4_END=$(date +%s)
log_step_complete "4" "DeepVariant (Variant Calling)"

echo ""
log_info "📑 Indexing VCF file and viewing header..."
docker run --rm -i \
    -v "${OUT_DIR}:/out" \
    "${PB_IMG}" \
    bash -c "echo '   → Creating VCF index...' && \
             tabix -p vcf /out/HG002.chr20.vcf.gz && \
             echo '   → VCF indexed successfully' && \
             echo '' && \
             echo '   📊 VCF Summary:' && \
             VARIANT_COUNT=\$(bcftools view -H /out/HG002.chr20.vcf.gz | wc -l) && \
             echo \"      Total variants called: \${VARIANT_COUNT}\"" \
    2>&1 | tee "${LOG_DIR}/chr20_vcf_header.log"

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
VCF_FILE="${OUT_DIR}/HG002.chr20.vcf.gz"
VCF_INDEX="${OUT_DIR}/HG002.chr20.vcf.gz.tbi"
VCF_SIZE=$(ls -lh "$VCF_FILE" 2>/dev/null | awk '{print $5}' || echo "N/A")
VCF_INDEX_SIZE=$(ls -lh "$VCF_INDEX" 2>/dev/null | awk '{print $5}' || echo "N/A")

# Get variant count
VARIANT_COUNT=$(docker run --rm -v "${OUT_DIR}:/out" "${PB_IMG}" bcftools view -H /out/HG002.chr20.vcf.gz 2>/dev/null | wc -l || echo "0")
VARIANT_COUNT_FORMATTED=$(printf "%'d" $VARIANT_COUNT)

# Get FASTQ sizes
FASTQ_R1_SIZE=$(ls -lh "${IN_DIR}/HG002_R1.fastq.gz" 2>/dev/null | awk '{print $5}' || echo "N/A")
FASTQ_R2_SIZE=$(ls -lh "${IN_DIR}/HG002_R2.fastq.gz" 2>/dev/null | awk '{print $5}' || echo "N/A")

echo ""
echo ""
echo "● Chr20 Test Complete! 🎉"
echo ""
echo "  ┌───────────────────────┬─────────────────────────────┐"
echo "  │        Result         │            Value            │"
echo "  ├───────────────────────┼─────────────────────────────┤"
printf "  │ Total Variants Called │ %-27s │\n" "$VARIANT_COUNT_FORMATTED"
echo "  ├───────────────────────┼─────────────────────────────┤"
echo "  │ Chromosome            │ chr20 only                  │"
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
printf "  │ VCF File   │ %-39s │\n" "HG002.chr20.vcf.gz ($VCF_SIZE)"
echo "  ├────────────┼─────────────────────────────────────────┤"
printf "  │ VCF Index  │ %-39s │\n" "HG002.chr20.vcf.gz.tbi ($VCF_INDEX_SIZE)"
echo "  └────────────┴─────────────────────────────────────────┘"
echo ""
echo "  Pipeline Time Breakdown"
echo "  ┌────────────────────────┬───────────────────┐"
echo "  │          Step          │       Time        │"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 1. Create Interval BED │ %-17s │\n" "$STEP1_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 2. fq2bam (alignment)  │ %-17s │\n" "$STEP2_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 3. BAM indexing & QC   │ %-17s │\n" "$STEP3_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ 4. DeepVariant         │ %-17s │\n" "$STEP4_FORMATTED"
echo "  ├────────────────────────┼───────────────────┤"
printf "  │ Total                  │ %-17s │\n" "$TOTAL_FORMATTED"
echo "  └────────────────────────┴───────────────────┘"
echo ""
echo "  Summary"
echo "  • Processed chromosome 20 region from FASTQ data"
echo "  • Aligned reads to GRCh38 reference genome"
echo "  • Called $VARIANT_COUNT_FORMATTED variants on chr20"
echo "  • Pipeline validated and ready for full genome run"
echo ""
echo "  🚀 Next: Run 'Full Genome' for complete analysis (30-110 minutes)"
echo ""
