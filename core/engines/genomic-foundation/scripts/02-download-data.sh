#!/bin/bash

# Immediate output test - this should appear right away
echo "SCRIPT STARTED - $(date)"
echo "Testing output streaming..."

set -e

# Disable output buffering completely
export PYTHONUNBUFFERED=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "Project root: $PROJECT_ROOT"
echo "Script directory: $SCRIPT_DIR"

# Source environment
source "${PROJECT_ROOT}/config/pipeline.env"

# Use local paths
DATA_DIR="${PROJECT_ROOT}/data/input/giab_hg002"
IN_DIR="${PROJECT_ROOT}/data/input"

# Helper function for progress messages
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

echo "=========================================="
echo "📥 GIAB HG002 FASTQ Data Download"
echo "=========================================="
echo ""
echo "ℹ️  This script will download the Genome in a Bottle (GIAB) HG002 reference sample"
echo "ℹ️  Dataset: Illumina 2x250bp paired-end sequencing data"
echo "ℹ️  Total size: ~200GB (will take several hours)"
echo "ℹ️  Source: NCBI GIAB FTP server"
echo ""
echo "⚠️  Warning: This will download ~200GB of data!"
echo "ℹ️  Auto-proceeding with download in 3 seconds..."
sleep 3
echo "✅ Starting download process..."
echo ""

log_progress "STEP 1/7: Preparing Download Directory"
log_info "Creating directory structure..."
mkdir -p "${DATA_DIR}/reads"
cd "${DATA_DIR}"
log_success "Directory prepared: ${DATA_DIR}"

# Check for aria2
log_progress "STEP 2/7: Installing Download Tools"
if ! command -v aria2c &> /dev/null; then
    log_info "aria2c not found - installing for high-speed parallel downloads..."
    sudo apt-get update -qq
    sudo apt-get install -y aria2
    log_success "aria2c installed successfully"
else
    log_success "aria2c already installed"
fi

log_progress "STEP 3/7: Fetching File Index"
log_info "Downloading sequence index from NCBI..."
log_info "URL: ${GIAB_INDEX_URL}"
wget -q --show-progress -O hg002_illumina_2x250.index.tsv "${GIAB_INDEX_URL}" 2>&1 | \
    grep --line-buffered -oP '\d+%|\d+[KMG]'
if [ ! -s hg002_illumina_2x250.index.tsv ]; then
    echo "Error: Index file download failed or is empty"
    exit 1
fi
log_success "Index file downloaded"

FILE_COUNT=$(awk -F'\t' 'NF>3' hg002_illumina_2x250.index.tsv | wc -l)
log_info "Found ${FILE_COUNT} paired-end FASTQ files (${FILE_COUNT} R1 + ${FILE_COUNT} R2)"

log_progress "STEP 4/7: Parsing Download URLs"
log_info "Extracting R1 (forward read) URLs..."
awk -F'\t' 'NF>3{print $1}' hg002_illumina_2x250.index.tsv > urls_R1.txt
log_success "Extracted ${FILE_COUNT} R1 URLs"

log_info "Extracting R2 (reverse read) URLs..."
awk -F'\t' 'NF>3{print $3}' hg002_illumina_2x250.index.tsv > urls_R2.txt
log_success "Extracted ${FILE_COUNT} R2 URLs"

log_progress "STEP 5/7: Downloading R1 Files (Forward Reads)"
log_info "Downloading ${FILE_COUNT} R1 FASTQ files in parallel..."
log_info "Using 16 connections per file for maximum speed"
log_info "This will take 1-3 hours depending on your internet connection"
echo ""
echo "⏳ Download in progress - aria2c will show periodic updates..."
echo ""

mkdir -p reads
aria2c \
    --console-log-level=notice \
    --summary-interval=5 \
    -x 8 \
    -s 8 \
    -j 4 \
    --file-allocation=none \
    -d reads \
    -i urls_R1.txt

ARIA2_EXIT=$?
if [ $ARIA2_EXIT -ne 0 ]; then
    echo "❌ Error: R1 download failed (aria2c exit code: $ARIA2_EXIT)"
    echo "ℹ️  Check network connectivity and disk space, then re-run this script."
    exit 1
fi

echo ""
log_success "All R1 files downloaded successfully"

log_progress "STEP 6/7: Downloading R2 Files (Reverse Reads)"
log_info "Downloading ${FILE_COUNT} R2 FASTQ files in parallel..."
log_info "This will take 1-3 hours depending on your internet connection"
echo ""
echo "⏳ Download in progress - aria2c will show periodic updates..."
echo ""

aria2c \
    --console-log-level=notice \
    --summary-interval=5 \
    -x 8 \
    -s 8 \
    -j 4 \
    --file-allocation=none \
    -d reads \
    -i urls_R2.txt

ARIA2_EXIT=$?
if [ $ARIA2_EXIT -ne 0 ]; then
    echo "❌ Error: R2 download failed (aria2c exit code: $ARIA2_EXIT)"
    echo "ℹ️  Check network connectivity and disk space, then re-run this script."
    exit 1
fi

echo ""
log_success "All R2 files downloaded successfully"

log_progress "STEP 7/7: Merging Lane FASTQs"
log_info "Combining multiple sequencing lanes into single R1 and R2 files..."

cd reads
ls -1 *R1*.f*q.gz 2>/dev/null | sort > R1_files.txt
R1_FILE_COUNT=$(wc -l < R1_files.txt)
if [ "$R1_FILE_COUNT" -eq 0 ]; then
    echo "❌ Error: No R1 FASTQ files found after download"
    exit 1
fi

ls -1 *R2*.f*q.gz 2>/dev/null | sort > R2_files.txt
R2_FILE_COUNT=$(wc -l < R2_files.txt)
if [ "$R2_FILE_COUNT" -eq 0 ]; then
    echo "❌ Error: No R2 FASTQ files found after download"
    exit 1
fi

log_info "Merging ${R1_FILE_COUNT} R1 lane files..."
echo "🔄 Concatenating R1 gzip files (valid per gzip spec, much faster than decompress/recompress)..."
cat $(cat R1_files.txt) > ../HG002_R1.fastq.gz
log_success "R1 files merged into HG002_R1.fastq.gz"

log_info "Merging ${R2_FILE_COUNT} R2 lane files..."
echo "🔄 Concatenating R2 gzip files..."
cat $(cat R2_files.txt) > ../HG002_R2.fastq.gz
log_success "R2 files merged into HG002_R2.fastq.gz"

cd ..

echo ""
log_info "Validating merged FASTQ files..."
echo "🔍 Counting reads in R1 file..."
R1_COUNT=$(zcat HG002_R1.fastq.gz | wc -l | awk '{print $1/4}')
echo "🔍 Counting reads in R2 file..."
R2_COUNT=$(zcat HG002_R2.fastq.gz | wc -l | awk '{print $1/4}')

log_info "R1 reads: $(printf "%'d" $R1_COUNT)"
log_info "R2 reads: $(printf "%'d" $R2_COUNT)"

if [ "$R1_COUNT" != "$R2_COUNT" ]; then
    echo "❌ Error: R1 and R2 read counts do not match!"
    exit 1
fi
log_success "Read counts match - paired-end data is valid!"

echo ""
log_info "Copying merged files to pipeline input directory..."
cp -v HG002_R1.fastq.gz "${IN_DIR}/HG002_R1.fastq.gz"
cp -v HG002_R2.fastq.gz "${IN_DIR}/HG002_R2.fastq.gz"

echo ""
echo "=========================================="
log_success "✅ DOWNLOAD COMPLETE!"
echo "=========================================="
echo ""
log_info "Downloaded and merged GIAB HG002 FASTQ data"
log_info "Total reads: $(printf "%'d" $R1_COUNT) paired-end reads"
echo ""
echo "📁 Output Files:"
ls -lh "${IN_DIR}"/HG002_*.fastq.gz | awk '{printf "   %s  %s\n", $5, $9}'
echo ""
log_info "Next step: ./scripts/03-setup-reference.sh (Setup Reference Genome)"
echo ""
