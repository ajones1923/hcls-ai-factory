#!/bin/bash
set -e

echo "=========================================="
echo "GIAB HG002 Conservative Download"
echo "=========================================="
echo ""
echo "Using conservative settings for reliability:"
echo "  - One file at a time"
echo "  - 8 connections per file (reduced from 16)"
echo "  - Verify MD5 after each download"
echo "  - Keep verified files from previous attempts"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "${PROJECT_ROOT}/config/pipeline.env"

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

cd "${DATA_DIR}"

if [ ! -f hg002_illumina_2x250.index.tsv ]; then
    log_error "Index file not found. Run download script first."
    exit 1
fi

# Extract checksums if not already done
if [ ! -f md5_R1.txt ]; then
    awk -F'\t' 'NR>1 && NF>3{print $2}' hg002_illumina_2x250.index.tsv > md5_R1.txt
    awk -F'\t' 'NR>1 && NF>3{print $4}' hg002_illumina_2x250.index.tsv > md5_R2.txt
fi

if [ ! -f urls_R1.txt ]; then
    awk -F'\t' 'NR>1 && NF>3{print $1}' hg002_illumina_2x250.index.tsv > urls_R1.txt
    awk -F'\t' 'NR>1 && NF>3{print $3}' hg002_illumina_2x250.index.tsv > urls_R2.txt
fi

mkdir -p reads
cd reads

# Function to download and verify a single file
download_and_verify() {
    local url="$1"
    local expected_md5="$2"
    local filename=$(basename "$url")

    # Check if file already exists and is valid
    if [ -f "$filename" ]; then
        local current_md5=$(md5sum "$filename" 2>/dev/null | awk '{print $1}')
        if [ "$current_md5" == "$expected_md5" ]; then
            echo "✓ $filename already verified, skipping"
            return 0
        else
            echo "⚠ $filename exists but invalid, re-downloading"
            rm -f "$filename"
        fi
    fi

    # Download with conservative settings
    echo "📥 Downloading: $filename"

    local attempt=1
    local max_attempts=5

    while [ $attempt -le $max_attempts ]; do
        echo "  Attempt $attempt of $max_attempts..."

        # Download with aria2c - conservative settings
        aria2c \
            --console-log-level=warn \
            --summary-interval=0 \
            -x 8 \
            -s 8 \
            -j 1 \
            --max-tries=5 \
            --retry-wait=3 \
            --max-connection-per-server=8 \
            --split=8 \
            --file-allocation=none \
            --allow-overwrite=true \
            "$url" \
            2>&1 | grep -v "^$" || true

        # Verify MD5
        if [ -f "$filename" ]; then
            local downloaded_md5=$(md5sum "$filename" 2>/dev/null | awk '{print $1}')
            if [ "$downloaded_md5" == "$expected_md5" ]; then
                echo "✅ $filename verified successfully"
                return 0
            else
                log_error "MD5 mismatch on attempt $attempt"
                echo "  Expected: $expected_md5"
                echo "  Got:      $downloaded_md5"
                rm -f "$filename"
                attempt=$((attempt + 1))
                if [ $attempt -le $max_attempts ]; then
                    echo "  Waiting 5 seconds before retry..."
                    sleep 5
                fi
            fi
        else
            log_error "Download failed on attempt $attempt"
            attempt=$((attempt + 1))
            if [ $attempt -le $max_attempts ]; then
                sleep 5
            fi
        fi
    done

    log_error "$filename failed after $max_attempts attempts"
    return 1
}

log_progress "Downloading and Verifying R1 Files"

TOTAL_R1=$(wc -l < ../urls_R1.txt)
CURRENT=0
FAILED_R1=()

while IFS= read -r url && IFS= read -r md5 <&3; do
    CURRENT=$((CURRENT + 1))
    echo ""
    echo "[$CURRENT/$TOTAL_R1] Processing $(basename "$url")"

    if ! download_and_verify "$url" "$md5"; then
        FAILED_R1+=("$(basename "$url")")
    fi
done < ../urls_R1.txt 3< ../md5_R1.txt

echo ""
log_progress "Downloading and Verifying R2 Files"

TOTAL_R2=$(wc -l < ../urls_R2.txt)
CURRENT=0
FAILED_R2=()

while IFS= read -r url && IFS= read -r md5 <&3; do
    CURRENT=$((CURRENT + 1))
    echo ""
    echo "[$CURRENT/$TOTAL_R2] Processing $(basename "$url")"

    if ! download_and_verify "$url" "$md5"; then
        FAILED_R2+=("$(basename "$url")")
    fi
done < ../urls_R2.txt 3< ../md5_R2.txt

cd ..

echo ""
echo "=========================================="
echo "Download Summary"
echo "=========================================="
echo ""

R1_SUCCESS=$((TOTAL_R1 - ${#FAILED_R1[@]}))
R2_SUCCESS=$((TOTAL_R2 - ${#FAILED_R2[@]}))

echo "R1 Files:"
echo "  ✅ Successful: $R1_SUCCESS / $TOTAL_R1"
echo "  ❌ Failed: ${#FAILED_R1[@]}"

if [ ${#FAILED_R1[@]} -gt 0 ]; then
    echo "  Failed files:"
    for f in "${FAILED_R1[@]}"; do
        echo "    - $f"
    done
fi

echo ""
echo "R2 Files:"
echo "  ✅ Successful: $R2_SUCCESS / $TOTAL_R2"
echo "  ❌ Failed: ${#FAILED_R2[@]}"

if [ ${#FAILED_R2[@]} -gt 0 ]; then
    echo "  Failed files:"
    for f in "${FAILED_R2[@]}"; do
        echo "    - $f"
    done
fi

echo ""
TOTAL_SIZE=$(du -sh reads 2>/dev/null | awk '{print $1}')
echo "Total downloaded: $TOTAL_SIZE"
echo ""

if [ ${#FAILED_R1[@]} -eq 0 ] && [ ${#FAILED_R2[@]} -eq 0 ]; then
    log_success "✅ All files downloaded and verified successfully!"
    echo ""
    log_info "Next step: Merge files with pigz"
    echo "  ./scripts/fast-merge-fixed.sh"
else
    log_error "Some files failed to download"
    echo ""
    echo "You can:"
    echo "  1. Run this script again to retry failed files"
    echo "  2. Check network connection and try later"
    exit 1
fi
