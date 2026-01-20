#!/bin/bash
set -e

echo "=========================================="
echo "GIAB HG002 FASTQ Download (Verified)"
echo "=========================================="
echo ""
echo "This script downloads with MD5 verification"
echo "to ensure data integrity"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Source environment
source "${PROJECT_ROOT}/config/pipeline.env"

# Use local paths
DATA_DIR="${PROJECT_ROOT}/data/input/giab_hg002"
IN_DIR="${PROJECT_ROOT}/data/input"

# Helper functions
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

log_error() {
    echo "❌ $1"
}

echo "ℹ️  Dataset: Illumina 2x250bp paired-end sequencing"
echo "ℹ️  Total size: ~300GB (will take several hours)"
echo "ℹ️  Source: NCBI GIAB FTP server"
echo ""

log_progress "STEP 1/8: Preparing Download Directory"
mkdir -p "${DATA_DIR}/reads"
cd "${DATA_DIR}"
log_success "Directory prepared: ${DATA_DIR}"

log_progress "STEP 2/8: Checking aria2c"
if ! command -v aria2c &> /dev/null; then
    log_info "Installing aria2c..."
    sudo apt-get update -qq
    sudo apt-get install -y aria2
    log_success "aria2c installed"
else
    log_success "aria2c already installed"
fi

log_progress "STEP 3/8: Fetching File Index with MD5 Checksums"
if [ ! -f hg002_illumina_2x250.index.tsv ]; then
    log_info "Downloading index from NCBI..."
    wget -q --show-progress -O hg002_illumina_2x250.index.tsv "${GIAB_INDEX_URL}"
    log_success "Index downloaded"
else
    log_success "Index already exists"
fi

FILE_COUNT=$(awk -F'\t' 'NR>1 && NF>3' hg002_illumina_2x250.index.tsv | wc -l)
log_info "Found ${FILE_COUNT} file pairs in index"

log_progress "STEP 4/8: Extracting URLs and MD5 Checksums"
# Extract R1 URLs and MD5s
awk -F'\t' 'NR>1 && NF>3{print $1}' hg002_illumina_2x250.index.tsv > urls_R1.txt
awk -F'\t' 'NR>1 && NF>3{print $2}' hg002_illumina_2x250.index.tsv > md5_R1.txt

# Extract R2 URLs and MD5s
awk -F'\t' 'NR>1 && NF>3{print $3}' hg002_illumina_2x250.index.tsv > urls_R2.txt
awk -F'\t' 'NR>1 && NF>3{print $4}' hg002_illumina_2x250.index.tsv > md5_R2.txt

log_success "Extracted ${FILE_COUNT} R1 URLs and checksums"
log_success "Extracted ${FILE_COUNT} R2 URLs and checksums"

log_progress "STEP 5/8: Downloading R1 Files with Verification"
log_info "Using aria2c with:"
log_info "  - 16 connections per file for speed"
log_info "  - Auto-resume on connection failures"
log_info "  - Checksum verification after download"
echo ""

cd reads

# Download R1 files
log_info "Starting R1 downloads..."
aria2c \
    --console-log-level=notice \
    --summary-interval=10 \
    -x 16 \
    -s 16 \
    -j 4 \
    --max-tries=10 \
    --retry-wait=5 \
    --max-connection-per-server=16 \
    --split=16 \
    --file-allocation=none \
    -i ../urls_R1.txt \
    2>&1 | tee -a /tmp/aria2c_r1.log

log_success "R1 downloads complete"

log_progress "STEP 6/8: Verifying R1 Files"
log_info "Checking MD5 checksums for all R1 files..."
FAILED=0
while IFS= read -r url && IFS= read -r md5 <&3; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        echo -n "Verifying $filename... "
        calculated_md5=$(md5sum "$filename" | awk '{print $1}')
        if [ "$calculated_md5" == "$md5" ]; then
            echo "✓ OK"
        else
            echo "✗ FAILED (expected: $md5, got: $calculated_md5)"
            FAILED=$((FAILED + 1))
            rm -f "$filename"
        fi
    else
        echo "✗ MISSING: $filename"
        FAILED=$((FAILED + 1))
    fi
done < ../urls_R1.txt 3< ../md5_R1.txt

if [ $FAILED -gt 0 ]; then
    log_error "$FAILED R1 files failed verification - will retry"
    exit 1
fi
log_success "All R1 files verified successfully"

log_progress "STEP 7/8: Downloading R2 Files with Verification"
log_info "Starting R2 downloads..."
aria2c \
    --console-log-level=notice \
    --summary-interval=10 \
    -x 16 \
    -s 16 \
    -j 4 \
    --max-tries=10 \
    --retry-wait=5 \
    --max-connection-per-server=16 \
    --split=16 \
    --file-allocation=none \
    -i ../urls_R2.txt \
    2>&1 | tee -a /tmp/aria2c_r2.log

log_success "R2 downloads complete"

log_progress "STEP 8/8: Verifying R2 Files"
log_info "Checking MD5 checksums for all R2 files..."
FAILED=0
while IFS= read -r url && IFS= read -r md5 <&3; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        echo -n "Verifying $filename... "
        calculated_md5=$(md5sum "$filename" | awk '{print $1}')
        if [ "$calculated_md5" == "$md5" ]; then
            echo "✓ OK"
        else
            echo "✗ FAILED (expected: $md5, got: $calculated_md5)"
            FAILED=$((FAILED + 1))
            rm -f "$filename"
        fi
    else
        echo "✗ MISSING: $filename"
        FAILED=$((FAILED + 1))
    fi
done < ../urls_R2.txt 3< ../md5_R2.txt

if [ $FAILED -gt 0 ]; then
    log_error "$FAILED R2 files failed verification - will retry"
    exit 1
fi
log_success "All R2 files verified successfully"

cd ..

echo ""
echo "=========================================="
log_success "✅ DOWNLOAD COMPLETE AND VERIFIED!"
echo "=========================================="
echo ""
log_info "All files downloaded and checksums verified"
log_info "Total files: $((FILE_COUNT * 2))"
ls -lh reads/*.fastq.gz | wc -l | awk '{print "Files on disk: " $1}'
du -sh reads | awk '{print "Total size: " $1}'
echo ""
log_info "Next step: Run merge script to combine files"
echo ""
