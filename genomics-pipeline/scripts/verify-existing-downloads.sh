#!/bin/bash

echo "=========================================="
echo "Verifying Existing Downloaded Files"
echo "=========================================="
echo ""

cd /home/adam/transfer/genomics-pipeline/data/input/giab_hg002

if [ ! -f hg002_illumina_2x250.index.tsv ]; then
    echo "❌ Index file not found"
    exit 1
fi

cd reads

echo "Extracting checksums from index..."
awk -F'\t' 'NR>1 && NF>3{print $2}' ../hg002_illumina_2x250.index.tsv > ../md5_R1.txt
awk -F'\t' 'NR>1 && NF>3{print $4}' ../hg002_illumina_2x250.index.tsv > ../md5_R2.txt
awk -F'\t' 'NR>1 && NF>3{print $1}' ../hg002_illumina_2x250.index.tsv > ../urls_R1.txt
awk -F'\t' 'NR>1 && NF>3{print $3}' ../hg002_illumina_2x250.index.tsv > ../urls_R2.txt

echo ""
echo "Verifying R1 files..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

GOOD_R1=0
BAD_R1=0
MISSING_R1=0

while IFS= read -r url && IFS= read -r md5 <&3; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        calculated_md5=$(md5sum "$filename" 2>/dev/null | awk '{print $1}')
        if [ "$calculated_md5" == "$md5" ]; then
            echo "✓ $filename"
            GOOD_R1=$((GOOD_R1 + 1))
        else
            echo "✗ $filename (MD5 mismatch)"
            BAD_R1=$((BAD_R1 + 1))
        fi
    else
        echo "⊘ $filename (missing)"
        MISSING_R1=$((MISSING_R1 + 1))
    fi
done < ../urls_R1.txt 3< ../md5_R1.txt

echo ""
echo "Verifying R2 files..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

GOOD_R2=0
BAD_R2=0
MISSING_R2=0

while IFS= read -r url && IFS= read -r md5 <&3; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        calculated_md5=$(md5sum "$filename" 2>/dev/null | awk '{print $1}')
        if [ "$calculated_md5" == "$md5" ]; then
            echo "✓ $filename"
            GOOD_R2=$((GOOD_R2 + 1))
        else
            echo "✗ $filename (MD5 mismatch)"
            BAD_R2=$((BAD_R2 + 1))
        fi
    else
        echo "⊘ $filename (missing)"
        MISSING_R2=$((MISSING_R2 + 1))
    fi
done < ../urls_R2.txt 3< ../md5_R2.txt

echo ""
echo "=========================================="
echo "Verification Summary"
echo "=========================================="
echo ""
echo "R1 Files:"
echo "  ✓ Good: $GOOD_R1"
echo "  ✗ Bad:  $BAD_R1"
echo "  ⊘ Missing: $MISSING_R1"
echo ""
echo "R2 Files:"
echo "  ✓ Good: $GOOD_R2"
echo "  ✗ Bad:  $BAD_R2"
echo "  ⊘ Missing: $MISSING_R2"
echo ""

TOTAL_EXPECTED=$(($(wc -l < ../urls_R1.txt) * 2))
TOTAL_GOOD=$((GOOD_R1 + GOOD_R2))
TOTAL_BAD=$((BAD_R1 + BAD_R2 + MISSING_R1 + MISSING_R2))

echo "Overall:"
echo "  Total Expected: $TOTAL_EXPECTED"
echo "  Valid Files: $TOTAL_GOOD"
echo "  Invalid/Missing: $TOTAL_BAD"
echo ""

if [ $TOTAL_BAD -eq 0 ]; then
    echo "✅ All files are valid! Ready to merge."
else
    echo "⚠️  $TOTAL_BAD files need to be re-downloaded"
    echo ""
    echo "Options:"
    echo "  1. Re-download all files (safest)"
    echo "  2. Keep good files and only download bad ones"
fi
echo ""
