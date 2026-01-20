#!/bin/bash
set -e

echo "=========================================="
echo "Fast Parallel FASTQ Merge (with pigz)"
echo "=========================================="
echo ""

CORES=$(nproc)
echo "System: $CORES CPU cores"
echo "Using pigz (parallel gzip) for maximum speed"
echo ""

cd /home/adam/transfer/genomics-pipeline/data/input/giab_hg002/reads

echo "Creating file lists..."
ls -1 D1_S1_L001_R1_*.fastq.gz D1_S1_L002_R1_*.fastq.gz | sort > R1_files.txt
ls -1 D1_S1_L001_R2_*.fastq.gz D1_S1_L002_R2_*.fastq.gz | sort > R2_files.txt

R1_COUNT=$(wc -l < R1_files.txt)
R2_COUNT=$(wc -l < R2_files.txt)
echo "✓ Found $R1_COUNT R1 files and $R2_COUNT R2 files"
echo ""

echo "=========================================="
echo "Merging R1 files (Lane 001 + Lane 002)"
echo "=========================================="
echo "Started at: $(date '+%I:%M %p')"
echo ""

# Use unpigz for parallel decompression, pigz for parallel compression
time unpigz -c $(cat R1_files.txt) | pigz -p $CORES > ../HG002_R1.fastq.gz

echo ""
echo "✓ R1 merge complete at: $(date '+%I:%M %p')"
R1_SIZE=$(du -h ../HG002_R1.fastq.gz | awk '{print $1}')
echo "   File size: $R1_SIZE"
echo ""

echo "=========================================="
echo "Merging R2 files (Lane 001 + Lane 002)"
echo "=========================================="
echo "Started at: $(date '+%I:%M %p')"
echo ""

time unpigz -c $(cat R2_files.txt) | pigz -p $CORES > ../HG002_R2.fastq.gz

echo ""
echo "✓ R2 merge complete at: $(date '+%I:%M %p')"
R2_SIZE=$(du -h ../HG002_R2.fastq.gz | awk '{print $1}')
echo "   File size: $R2_SIZE"
echo ""

cd ..

echo "Verifying files..."
echo "Counting reads (this may take a minute)..."
R1_READS=$(timeout 120 pigz -dc HG002_R1.fastq.gz | wc -l | awk '{print $1/4}')
R2_READS=$(timeout 120 pigz -dc HG002_R2.fastq.gz | wc -l | awk '{print $1/4}')

echo ""
echo "=========================================="
echo "✅ Parallel Merge Complete!"
echo "=========================================="
echo "R1: $(printf "%'d" $R1_READS) reads ($R1_SIZE)"
echo "R2: $(printf "%'d" $R2_READS) reads ($R2_SIZE)"
echo ""
echo "Copying to pipeline input directory..."
cp -v HG002_R1.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
cp -v HG002_R2.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
echo ""
echo "✅ Files ready for Chr20 Test and Full Genome Pipeline!"
echo "   Used all $CORES CPU cores for processing"
echo ""
