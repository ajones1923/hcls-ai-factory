#!/bin/bash
set -e

echo "=========================================="
echo "Re-merging FASTQ Files"
echo "=========================================="
echo ""
echo "This will properly merge the downloaded FASTQ chunks"
echo "using decompression and recompression instead of simple cat"
echo ""

cd /home/adam/transfer/genomics-pipeline/data/input/giab_hg002/reads

echo "Deleting old corrupted merged files..."
rm -f ../HG002_R1.fastq.gz ../HG002_R2.fastq.gz
echo "✓ Old files deleted"
echo ""

echo "Creating file lists..."
ls -1 D1_S1_L001_R1_*.fastq.gz | sort > R1_files.txt
ls -1 D1_S1_L001_R2_*.fastq.gz | sort > R2_files.txt

R1_COUNT=$(wc -l < R1_files.txt)
R2_COUNT=$(wc -l < R2_files.txt)
echo "✓ Found $R1_COUNT R1 files and $R2_COUNT R2 files"
echo ""

echo "Merging R1 files..."
echo "⏳ Decompressing and recompressing (30-60 minutes)..."
echo "Started at: $(date)"
time zcat $(cat R1_files.txt) | gzip -c > ../HG002_R1.fastq.gz
echo "✓ R1 merge complete at: $(date)"
echo ""

echo "Merging R2 files..."
echo "⏳ Decompressing and recompressing (30-60 minutes)..."
echo "Started at: $(date)"
time zcat $(cat R2_files.txt) | gzip -c > ../HG002_R2.fastq.gz
echo "✓ R2 merge complete at: $(date)"
echo ""

cd ..
echo "Verifying merged files..."
echo "Counting reads in R1..."
R1_READS=$(timeout 120 zcat HG002_R1.fastq.gz | wc -l | awk '{print $1/4}')
echo "Counting reads in R2..."
R2_READS=$(timeout 120 zcat HG002_R2.fastq.gz | wc -l | awk '{print $1/4}')

echo ""
echo "=========================================="
echo "✅ Merge Complete!"
echo "=========================================="
echo "R1 reads: $(printf "%'d" $R1_READS)"
echo "R2 reads: $(printf "%'d" $R2_READS)"
ls -lh HG002_*.fastq.gz
echo ""
echo "Copying to input directory..."
cp -v HG002_R1.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
cp -v HG002_R2.fastq.gz /home/adam/transfer/genomics-pipeline/data/input/
echo "✅ Done! Files are ready for pipeline."
