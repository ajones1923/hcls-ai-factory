#!/bin/bash
set -e

echo "=========================================="
echo "Fast Parallel FASTQ Merge (zcat + pigz)"
echo "=========================================="
echo ""
echo "Using:"
echo "  - zcat for decompression (handles quirky gzip files)"
echo "  - pigz for parallel compression (uses all CPU cores)"
echo ""

CORES=$(nproc)
echo "System: $CORES CPU cores"
echo ""

cd /home/adam/transfer/genomics-pipeline/data/input/giab_hg002

# Clean up previous attempts
rm -f HG002_R1.fastq.gz HG002_R2.fastq.gz
cd reads

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

# Use regular zcat for decompression, pigz for compression
time zcat $(cat R1_files.txt) | pigz -p $CORES > ../HG002_R1.fastq.gz

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

time zcat $(cat R2_files.txt) | pigz -p $CORES > ../HG002_R2.fastq.gz

echo ""
echo "✓ R2 merge complete at: $(date '+%I:%M %p')"
R2_SIZE=$(du -h ../HG002_R2.fastq.gz | awk '{print $1}')
echo "   File size: $R2_SIZE"
echo ""

cd ..

echo "Quick verification..."
R1_LINES=$(timeout 60 zcat HG002_R1.fastq.gz | head -1000000 | wc -l)
R2_LINES=$(timeout 60 zcat HG002_R2.fastq.gz | head -1000000 | wc -l)

echo ""
echo "=========================================="
echo "✅ Parallel Merge Complete!"
echo "=========================================="
echo "R1: $R1_SIZE ($((R1_LINES / 4)) reads in first million lines)"
echo "R2: $R2_SIZE ($((R2_LINES / 4)) reads in first million lines)"
echo ""
echo "Copying to pipeline input directory..."
cp -v HG002_R1.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
cp -v HG002_R2.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
echo ""
echo "✅ Files ready for Chr20 Test and Full Genome Pipeline!"
echo "   Compression used all $CORES CPU cores"
echo ""
