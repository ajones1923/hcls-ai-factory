#!/bin/bash
set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${PIPELINE_DIR}/data"

echo "=========================================="
echo "Fast Parallel FASTQ Re-merge"
echo "=========================================="
echo ""
echo "This uses pigz (parallel gzip) to utilize all CPU cores"
echo ""

# Install pigz if not present
if ! command -v pigz &> /dev/null; then
    echo "Installing pigz..."
    sudo apt-get update -qq
    sudo apt-get install -y pigz
    echo "✓ pigz installed"
fi

CORES=$(nproc)
echo "System: $CORES CPU cores detected"
echo "Will use all cores for parallel compression"
echo ""

cd "${DATA_DIR}/input/giab_hg002/reads"

echo "Cleaning up partial files..."
rm -f ../HG002_R1.fastq.gz ../HG002_R2.fastq.gz
echo "✓ Cleaned"
echo ""

echo "Creating file lists..."
ls -1 D1_S1_L001_R1_*.fastq.gz | sort > R1_files.txt
ls -1 D1_S1_L001_R2_*.fastq.gz | sort > R2_files.txt

R1_COUNT=$(wc -l < R1_files.txt)
R2_COUNT=$(wc -l < R2_files.txt)
echo "✓ Found $R1_COUNT R1 files and $R2_COUNT R2 files"
echo ""

echo "=========================================="
echo "Merging R1 files with pigz (parallel)"
echo "=========================================="
echo "⏳ Decompressing and parallel recompressing..."
echo "Started at: $(date)"
echo ""

# Use unpigz for decompression (parallel) and pigz for compression (parallel)
time unpigz -c $(cat R1_files.txt) | pigz -p $CORES > ../HG002_R1.fastq.gz

echo ""
echo "✓ R1 merge complete at: $(date)"
echo ""

echo "=========================================="
echo "Merging R2 files with pigz (parallel)"
echo "=========================================="
echo "⏳ Decompressing and parallel recompressing..."
echo "Started at: $(date)"
echo ""

time unpigz -c $(cat R2_files.txt) | pigz -p $CORES > ../HG002_R2.fastq.gz

echo ""
echo "✓ R2 merge complete at: $(date)"
echo ""

cd ..
echo "Verifying merged files..."
echo "Counting reads in R1..."
R1_READS=$(timeout 120 pigz -dc HG002_R1.fastq.gz | wc -l | awk '{print $1/4}')
echo "Counting reads in R2..."
R2_READS=$(timeout 120 pigz -dc HG002_R2.fastq.gz | wc -l | awk '{print $1/4}')

echo ""
echo "=========================================="
echo "✅ Fast Parallel Merge Complete!"
echo "=========================================="
echo "R1 reads: $(printf "%'d" $R1_READS)"
echo "R2 reads: $(printf "%'d" $R2_READS)"
echo ""
ls -lh HG002_*.fastq.gz
echo ""
echo "Copying to input directory..."
cp -v HG002_R1.fastq.gz "${DATA_DIR}/input/"
cp -v HG002_R2.fastq.gz "${DATA_DIR}/input/"
echo ""
echo "✅ Done! Files are ready for pipeline."
echo "   Processing used all $CORES cores for maximum speed!"
