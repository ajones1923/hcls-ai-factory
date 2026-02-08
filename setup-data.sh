#!/bin/bash
#
# HCLS AI Factory — Stage 0: Data Acquisition
# =============================================
# Downloads and verifies ALL required data files for all three pipeline stages.
#
# Usage: ./setup-data.sh [command] [options]
#
# Commands:
#   --all           Download all data for all stages (default)
#   --stage1        Download Stage 1 (Genomics) data only
#   --stage2        Download Stage 2 (RAG/Chat) data only
#   --stage3        Download Stage 3 (Drug Discovery) data only
#   --verify        Verify all downloaded data without downloading
#   --status        Show download status dashboard
#   --help          Show help message
#
# Options:
#   --skip-merge    Skip FASTQ merge step (Stage 1)
#   --dry-run       Show what would be downloaded without downloading
#   --retry N       Max retry attempts per file (default: 3)
#   --connections N aria2c connections per file (default: 8)
#

set -eo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Pipeline directories
GENOMICS_DIR="${SCRIPT_DIR}/genomics-pipeline"
RAG_CHAT_DIR="${SCRIPT_DIR}/rag-chat-pipeline"
DRUG_DISCOVERY_DIR="${SCRIPT_DIR}/drug-discovery-pipeline"

# Data URLs
GIAB_INDEX_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina_2x250bps_06012016_updated.HG002"
PARABRICKS_SAMPLE_URL="https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"
CLINVAR_SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_VCF_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_VCF_TBI_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
ALPHAMISSENSE_URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
PARABRICKS_IMG="nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

# VCP demo PDB IDs for optional pre-fetch
VCP_PDB_IDS="5FTK 7K56 8OOI 9DIL"

# State tracking
STATE_FILE="${SCRIPT_DIR}/.data-setup-state"

# Defaults
MAX_RETRIES=3
CONNECTIONS=8
SKIP_MERGE=false
DRY_RUN=false

# ============================================================================
# COLORS AND VISUAL HELPERS
# ============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m'
BOLD='\033[1m'
NVIDIA_GREEN='\033[38;5;118m'

print_banner() {
    echo ""
    echo -e "${NVIDIA_GREEN}+==============================================================================+${NC}"
    echo -e "${NVIDIA_GREEN}|${NC}                                                                              ${NVIDIA_GREEN}|${NC}"
    echo -e "${NVIDIA_GREEN}|${NC}    ${WHITE}${BOLD}HCLS AI FACTORY -- Stage 0: Data Acquisition${NC}                              ${NVIDIA_GREEN}|${NC}"
    echo -e "${NVIDIA_GREEN}|${NC}                                                                              ${NVIDIA_GREEN}|${NC}"
    echo -e "${NVIDIA_GREEN}|${NC}    ${CYAN}Download and verify all data required for the pipeline${NC}                    ${NVIDIA_GREEN}|${NC}"
    echo -e "${NVIDIA_GREEN}|${NC}                                                                              ${NVIDIA_GREEN}|${NC}"
    echo -e "${NVIDIA_GREEN}+==============================================================================+${NC}"
    echo ""
}

print_section() {
    echo ""
    echo -e "${NVIDIA_GREEN}================================================================================${NC}"
    echo -e "${WHITE}${BOLD}  $1${NC}"
    echo -e "${NVIDIA_GREEN}================================================================================${NC}"
    echo ""
}

log_info() {
    echo -e "  ${YELLOW}-->  $1${NC}"
}

log_success() {
    echo -e "  ${GREEN}[OK] $1${NC}"
}

log_error() {
    echo -e "  ${RED}[!!] $1${NC}"
}

log_skip() {
    echo -e "  ${CYAN}[--] $1${NC}"
}

log_warn() {
    echo -e "  ${YELLOW}[!!] $1${NC}"
}

print_progress() {
    local current=$1
    local total=$2
    local label="${3:-Progress}"
    local width=40

    if [ "$total" -eq 0 ]; then return; fi

    local percentage=$((current * 100 / total))
    local filled=$((current * width / total))
    local empty=$((width - filled))

    local bar=""
    for ((i = 0; i < filled; i++)); do bar+="="; done
    for ((i = 0; i < empty; i++)); do bar+="."; done

    printf "\r  ${CYAN}%s: [%s] %d/%d (%d%%)${NC}" "$label" "$bar" "$current" "$total" "$percentage"
}

format_size() {
    local size_bytes=$1
    if [ "$size_bytes" -ge 1073741824 ]; then
        echo "$(echo "scale=1; $size_bytes / 1073741824" | bc 2>/dev/null || echo "?") GB"
    elif [ "$size_bytes" -ge 1048576 ]; then
        echo "$(echo "scale=1; $size_bytes / 1048576" | bc 2>/dev/null || echo "?") MB"
    else
        echo "${size_bytes} bytes"
    fi
}

# ============================================================================
# STATE MANAGEMENT
# ============================================================================

is_verified() {
    local filepath="$1"
    local expected_md5="$2"

    if [ ! -f "$STATE_FILE" ]; then return 1; fi
    if [ -z "$expected_md5" ]; then
        grep -q "^${filepath}|" "$STATE_FILE" 2>/dev/null && return 0
    else
        grep -q "^${filepath}|${expected_md5}|" "$STATE_FILE" 2>/dev/null && return 0
    fi
    return 1
}

mark_verified() {
    local filepath="$1"
    local md5="$2"
    local timestamp
    timestamp=$(date -Iseconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S')

    # Remove old entry if exists
    if [ -f "$STATE_FILE" ]; then
        grep -v "^${filepath}|" "$STATE_FILE" > "${STATE_FILE}.tmp" 2>/dev/null || true
        mv "${STATE_FILE}.tmp" "$STATE_FILE"
    fi

    echo "${filepath}|${md5}|${timestamp}" >> "$STATE_FILE"
}

count_verified() {
    local pattern="${1:-}"
    if [ ! -f "$STATE_FILE" ]; then echo 0; return; fi
    if [ -z "$pattern" ]; then
        wc -l < "$STATE_FILE" | tr -d ' '
    else
        grep -c "$pattern" "$STATE_FILE" 2>/dev/null || echo 0
    fi
}

# ============================================================================
# PRE-FLIGHT CHECKS
# ============================================================================

preflight_checks() {
    local stages="$1"

    print_section "PRE-FLIGHT CHECKS"

    local all_ok=true

    # Check required tools
    echo -e "  ${BOLD}Required Tools${NC}"
    for tool in wget md5sum awk zcat; do
        if command -v "$tool" &>/dev/null; then
            log_success "$tool found"
        else
            log_error "$tool not found"
            all_ok=false
        fi
    done

    # aria2c (required for Stage 1)
    if echo "$stages" | grep -q "1"; then
        if command -v aria2c &>/dev/null; then
            log_success "aria2c found (parallel downloads)"
        else
            log_error "aria2c not found -- required for FASTQ downloads"
            echo ""
            echo -e "    ${YELLOW}Install with:${NC}"
            echo "      Ubuntu/Debian: sudo apt-get install -y aria2"
            echo "      CentOS/RHEL:   sudo yum install -y aria2"
            all_ok=false
        fi
    fi

    # pigz (required for merge)
    if echo "$stages" | grep -q "1" && [ "$SKIP_MERGE" = false ]; then
        if command -v pigz &>/dev/null; then
            log_success "pigz found (parallel compression)"
        else
            log_warn "pigz not found -- FASTQ merge will use gzip (slower)"
            echo ""
            echo -e "    ${YELLOW}Install with:${NC}"
            echo "      Ubuntu/Debian: sudo apt-get install -y pigz"
        fi
    fi

    # Check platform
    echo ""
    echo -e "  ${BOLD}Platform${NC}"
    local arch
    arch=$(uname -m)
    log_info "Architecture: $arch"
    log_info "OS: $(uname -s) $(uname -r)"

    # Check disk space
    echo ""
    echo -e "  ${BOLD}Disk Space${NC}"
    local available_kb
    available_kb=$(df -k "$SCRIPT_DIR" 2>/dev/null | tail -1 | awk '{print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    log_info "Available: ${available_gb} GB"

    local required_gb=0
    if echo "$stages" | grep -q "1"; then required_gb=$((required_gb + 350)); fi
    if echo "$stages" | grep -q "2"; then required_gb=$((required_gb + 3)); fi
    if echo "$stages" | grep -q "3"; then required_gb=$((required_gb + 1)); fi

    log_info "Required: ~${required_gb} GB"

    if [ "$available_gb" -lt "$required_gb" ]; then
        log_error "Insufficient disk space (need ${required_gb} GB, have ${available_gb} GB)"
        echo ""
        echo -e "    ${YELLOW}Options:${NC}"
        echo "      1. Free space: sudo apt-get clean && docker system prune"
        echo "      2. Download stages individually: ./setup-data.sh --stage2 (only ~2GB)"
        echo "      3. Symlink data directories to a larger volume"
        all_ok=false
    else
        log_success "Sufficient disk space"
    fi

    # Check network connectivity
    echo ""
    echo -e "  ${BOLD}Network Connectivity${NC}"

    if echo "$stages" | grep -q "1"; then
        if wget --spider -q --timeout=10 "https://ftp-trace.ncbi.nlm.nih.gov/" 2>/dev/null; then
            log_success "NCBI FTP (FASTQ source): reachable"
        else
            log_error "NCBI FTP: unreachable"
            echo -e "    ${YELLOW}Check: https://www.ncbi.nlm.nih.gov/Status/${NC}"
            all_ok=false
        fi

        if wget --spider -q --timeout=10 "$PARABRICKS_SAMPLE_URL" 2>/dev/null; then
            log_success "S3 (reference genome): reachable"
        else
            log_warn "S3 (reference genome): unreachable -- may be transient"
        fi
    fi

    if echo "$stages" | grep -q "2"; then
        if wget --spider -q --timeout=10 "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/" 2>/dev/null; then
            log_success "NCBI ClinVar: reachable"
        else
            log_error "NCBI ClinVar: unreachable"
            all_ok=false
        fi

        if wget --spider -q --timeout=10 "$ALPHAMISSENSE_URL" 2>/dev/null; then
            log_success "Google Cloud (AlphaMissense): reachable"
        else
            log_error "Google Cloud (AlphaMissense): unreachable"
            all_ok=false
        fi
    fi

    # Check pipeline directories exist
    echo ""
    echo -e "  ${BOLD}Pipeline Directories${NC}"
    for dir_name in "genomics-pipeline:1" "rag-chat-pipeline:2" "drug-discovery-pipeline:3"; do
        local dir="${dir_name%%:*}"
        local stage="${dir_name##*:}"
        if echo "$stages" | grep -q "$stage"; then
            if [ -d "${SCRIPT_DIR}/${dir}" ]; then
                log_success "$dir found"
            else
                log_error "$dir not found at ${SCRIPT_DIR}/${dir}"
                all_ok=false
            fi
        fi
    done

    echo ""
    if [ "$all_ok" = false ]; then
        log_error "Pre-flight checks failed. Fix the issues above and re-run."
        echo ""
        echo -e "  ${CYAN}For help: See docs/DATA_SETUP.md${NC}"
        echo -e "  ${CYAN}Report issues: https://github.com/NVIDIA/hcls-ai-factory/issues${NC}"
        exit 1
    fi

    log_success "All pre-flight checks passed"
}

# ============================================================================
# DOWNLOAD ENGINE
# ============================================================================

download_with_retry() {
    local url="$1"
    local target_path="$2"
    local expected_md5="${3:-}"
    local max_retries="${4:-$MAX_RETRIES}"
    local filename
    filename=$(basename "$target_path")
    local target_dir
    target_dir=$(dirname "$target_path")

    # Check state file for previous verification
    if [ -n "$expected_md5" ] && is_verified "$target_path" "$expected_md5"; then
        if [ -f "$target_path" ]; then
            log_skip "$filename (previously verified)"
            return 0
        fi
    fi

    # Check if file exists and verify
    if [ -f "$target_path" ]; then
        if [ -n "$expected_md5" ]; then
            local current_md5
            current_md5=$(md5sum "$target_path" 2>/dev/null | awk '{print $1}')
            if [ "$current_md5" = "$expected_md5" ]; then
                mark_verified "$target_path" "$expected_md5"
                log_skip "$filename (already verified)"
                return 0
            else
                log_warn "$filename exists but checksum mismatch, re-downloading"
                rm -f "$target_path"
            fi
        else
            # No MD5 available -- check gzip integrity
            if gzip -t "$target_path" 2>/dev/null; then
                local file_md5
                file_md5=$(md5sum "$target_path" 2>/dev/null | awk '{print $1}')
                mark_verified "$target_path" "$file_md5"
                log_skip "$filename (already present, gzip valid)"
                return 0
            else
                log_warn "$filename exists but corrupted, re-downloading"
                rm -f "$target_path"
            fi
        fi
    fi

    if [ "$DRY_RUN" = true ]; then
        log_info "Would download: $filename"
        return 0
    fi

    mkdir -p "$target_dir"

    for attempt in $(seq 1 "$max_retries"); do
        log_info "Downloading $filename (attempt $attempt/$max_retries)..."

        local dl_ok=false

        if [ "$attempt" -le 2 ] && command -v aria2c &>/dev/null; then
            # aria2c: fast parallel download
            local conn=$((CONNECTIONS / attempt))
            [ "$conn" -lt 2 ] && conn=2

            aria2c \
                --console-log-level=warn \
                --summary-interval=0 \
                -x "$conn" \
                -s "$conn" \
                -j 1 \
                --max-tries=3 \
                --retry-wait=5 \
                --max-connection-per-server="$conn" \
                --file-allocation=none \
                --allow-overwrite=true \
                -d "$target_dir" \
                -o "$filename" \
                "$url" 2>&1 | grep -v "^$" || true

            [ -f "$target_path" ] && dl_ok=true
        else
            # wget fallback: most reliable for FTP
            wget --continue --tries=3 --timeout=60 --waitretry=5 \
                -O "$target_path" "$url" 2>&1 || true

            [ -f "$target_path" ] && dl_ok=true
        fi

        if [ "$dl_ok" = true ]; then
            # Verify
            if [ -n "$expected_md5" ]; then
                local actual_md5
                actual_md5=$(md5sum "$target_path" 2>/dev/null | awk '{print $1}')
                if [ "$actual_md5" = "$expected_md5" ]; then
                    mark_verified "$target_path" "$expected_md5"
                    log_success "$filename verified (MD5 match)"
                    return 0
                else
                    log_error "$filename MD5 mismatch (expected: ${expected_md5:0:12}..., got: ${actual_md5:0:12}...)"
                    rm -f "$target_path"
                fi
            else
                if gzip -t "$target_path" 2>/dev/null; then
                    local file_md5
                    file_md5=$(md5sum "$target_path" 2>/dev/null | awk '{print $1}')
                    mark_verified "$target_path" "$file_md5"
                    log_success "$filename verified (gzip integrity)"
                    return 0
                else
                    log_error "$filename failed gzip integrity check"
                    rm -f "$target_path"
                fi
            fi
        else
            log_error "$filename download failed"
        fi

        # Exponential backoff before retry
        if [ "$attempt" -lt "$max_retries" ]; then
            local wait_time=$((attempt * 5))
            log_info "Waiting ${wait_time}s before retry..."
            sleep "$wait_time"
        fi
    done

    log_error "$filename FAILED after $max_retries attempts"
    return 1
}

# ============================================================================
# STAGE 1: GENOMICS DATA
# ============================================================================

download_stage1() {
    print_section "STAGE 1: GENOMICS DATA"

    echo -e "  ${BOLD}Components:${NC}"
    echo "    1a. HG002 FASTQ files (68 files, ~200 GB)"
    echo "    1b. GRCh38 reference genome (11 GB)"
    if [ "$SKIP_MERGE" = false ]; then
        echo "    1c. FASTQ merge (68 chunks -> 2 files)"
    fi
    echo ""

    # 1a. FASTQ download
    download_fastq_files

    # 1b. Reference genome
    download_reference

    # 1c. FASTQ merge
    if [ "$SKIP_MERGE" = false ]; then
        merge_fastq_files
    fi
}

download_fastq_files() {
    print_section "STAGE 1a: HG002 FASTQ FILES"

    local data_dir="${GENOMICS_DIR}/data/input/giab_hg002"
    local reads_dir="${data_dir}/reads"

    if [ "$DRY_RUN" = true ]; then
        log_info "Dry run: would download 68 FASTQ files (~200 GB) from GIAB"
        log_info "  Source: $GIAB_INDEX_URL"
        log_info "  Target: $reads_dir"
        return 0
    fi

    mkdir -p "$reads_dir" || {
        log_error "Cannot create directory: $reads_dir"
        log_error "Check permissions on $(dirname "$reads_dir")"
        return 1
    }

    # Download index TSV
    local index_file="${data_dir}/hg002_illumina_2x250.index.tsv"
    if [ ! -f "$index_file" ]; then
        log_info "Downloading GIAB HG002 file index..."
        wget -q --show-progress -O "$index_file" "$GIAB_INDEX_URL" 2>&1 || {
            log_error "Failed to download GIAB index file"
            diagnose_failure "GIAB index download" "connection_refused" ""
            return 1
        }
        log_success "Index downloaded"
    else
        log_skip "Index file already present"
    fi

    # Extract URLs and MD5s
    local urls_r1="${data_dir}/urls_R1.txt"
    local urls_r2="${data_dir}/urls_R2.txt"
    local md5_r1="${data_dir}/md5_R1.txt"
    local md5_r2="${data_dir}/md5_R2.txt"

    awk -F'\t' 'NR>1 && NF>3{print $1}' "$index_file" > "$urls_r1"
    awk -F'\t' 'NR>1 && NF>3{print $2}' "$index_file" > "$md5_r1"
    awk -F'\t' 'NR>1 && NF>3{print $3}' "$index_file" > "$urls_r2"
    awk -F'\t' 'NR>1 && NF>3{print $4}' "$index_file" > "$md5_r2"

    local total_r1
    total_r1=$(wc -l < "$urls_r1" | tr -d ' ')
    local total_r2
    total_r2=$(wc -l < "$urls_r2" | tr -d ' ')
    local total_files=$((total_r1 + total_r2))

    log_info "Found $total_r1 R1 files and $total_r2 R2 files ($total_files total)"
    echo ""

    # Download R1 files
    echo -e "  ${BOLD}Downloading R1 files (forward reads)${NC}"
    local current=0
    local failed_r1=0

    while IFS= read -r url && IFS= read -r md5 <&3; do
        current=$((current + 1))
        local filename
        filename=$(basename "$url")
        local target="${reads_dir}/${filename}"

        printf "\n  [%d/%d] %s\n" "$current" "$total_r1" "$filename"

        if ! download_with_retry "$url" "$target" "$md5"; then
            failed_r1=$((failed_r1 + 1))
        fi
    done < "$urls_r1" 3< "$md5_r1"

    echo ""

    # Download R2 files
    echo -e "  ${BOLD}Downloading R2 files (reverse reads)${NC}"
    current=0
    local failed_r2=0

    while IFS= read -r url && IFS= read -r md5 <&3; do
        current=$((current + 1))
        local filename
        filename=$(basename "$url")
        local target="${reads_dir}/${filename}"

        printf "\n  [%d/%d] %s\n" "$current" "$total_r2" "$filename"

        if ! download_with_retry "$url" "$target" "$md5"; then
            failed_r2=$((failed_r2 + 1))
        fi
    done < "$urls_r2" 3< "$md5_r2"

    echo ""

    # Summary
    local ok_r1=$((total_r1 - failed_r1))
    local ok_r2=$((total_r2 - failed_r2))

    echo -e "  ${BOLD}FASTQ Download Summary${NC}"
    echo "    R1: ${ok_r1}/${total_r1} verified"
    echo "    R2: ${ok_r2}/${total_r2} verified"

    if [ "$failed_r1" -gt 0 ] || [ "$failed_r2" -gt 0 ]; then
        local total_failed=$((failed_r1 + failed_r2))
        log_error "$total_failed files failed after all retries"
        echo ""
        diagnose_failure "FASTQ download" "md5_mismatch" "$total_failed files"
        return 1
    fi

    local total_size
    total_size=$(du -sh "$reads_dir" 2>/dev/null | awk '{print $1}')
    log_success "All $total_files FASTQ files verified ($total_size total)"
}

download_reference() {
    print_section "STAGE 1b: REFERENCE GENOME (GRCh38)"

    local ref_dir="${GENOMICS_DIR}/data/ref"

    # Check if reference is already complete (before mkdir to avoid permission issues)
    if [ -d "$ref_dir" ]; then
        local ref_complete=true
        for f in GRCh38.fa GRCh38.fa.fai GRCh38.fa.bwt GRCh38.fa.pac GRCh38.fa.sa GRCh38.fa.amb GRCh38.fa.ann; do
            if [ ! -f "${ref_dir}/${f}" ]; then
                ref_complete=false
                break
            fi
        done

        if [ "$ref_complete" = true ]; then
            log_skip "Reference genome already complete (all index files present)"
            return 0
        fi
    fi

    if [ "$DRY_RUN" = true ]; then
        log_info "Would download: parabricks_sample.tar.gz (~11 GB)"
        log_info "Would extract: GRCh38.fa + BWA index files"
        return 0
    fi

    mkdir -p "$ref_dir" || {
        log_error "Cannot create directory: $ref_dir"
        log_error "Check permissions on $(dirname "$ref_dir")"
        return 1
    }

    # Download Parabricks sample bundle
    local bundle="${ref_dir}/parabricks_sample.tar.gz"
    if [ ! -f "$bundle" ]; then
        log_info "Downloading Parabricks sample bundle (~11 GB)..."
        log_info "This will take 5-20 minutes depending on connection speed"

        wget --progress=dot:giga -O "$bundle" "$PARABRICKS_SAMPLE_URL" 2>&1 || {
            log_error "Failed to download reference genome bundle"
            diagnose_failure "Reference genome" "connection_refused" ""
            return 1
        }
        log_success "Bundle downloaded"
    else
        log_skip "Bundle already downloaded"
    fi

    # Extract
    if [ ! -f "${ref_dir}/GRCh38.fa" ]; then
        log_info "Extracting reference genome..."
        tar -xzf "$bundle" -C "$ref_dir" 2>&1
        cp "${ref_dir}/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta" "${ref_dir}/GRCh38.fa"
        log_success "Reference genome extracted"
    fi

    # Copy BWA index files
    if [ ! -f "${ref_dir}/GRCh38.fa.bwt" ]; then
        log_info "Copying BWA index files (5.3 GB)..."
        local pb_ref="${ref_dir}/parabricks_sample/Ref"
        for ext in amb ann bwt pac sa; do
            if [ -f "${pb_ref}/Homo_sapiens_assembly38.fasta.${ext}" ]; then
                cp "${pb_ref}/Homo_sapiens_assembly38.fasta.${ext}" "${ref_dir}/GRCh38.fa.${ext}"
            fi
        done
        log_success "BWA index files copied"
    else
        log_skip "BWA index files already present"
    fi

    # Create FASTA index
    if [ ! -f "${ref_dir}/GRCh38.fa.fai" ]; then
        if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
            log_info "Creating FASTA index via samtools..."
            docker run --rm \
                -v "${ref_dir}:/ref" \
                "$PARABRICKS_IMG" \
                bash -c "samtools faidx /ref/GRCh38.fa" 2>&1 || {
                log_warn "Docker samtools indexing failed. You can create it later with:"
                echo "    samtools faidx ${ref_dir}/GRCh38.fa"
            }
            if [ -f "${ref_dir}/GRCh38.fa.fai" ]; then
                log_success "FASTA index created"
            fi
        else
            log_warn "Docker not available. Create FASTA index manually:"
            echo "    samtools faidx ${ref_dir}/GRCh38.fa"
        fi
    else
        log_skip "FASTA index already exists"
    fi

    # Create sequence dictionary
    if [ ! -f "${ref_dir}/GRCh38.dict" ]; then
        if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
            log_info "Creating sequence dictionary..."
            docker run --rm \
                -v "${ref_dir}:/ref" \
                "$PARABRICKS_IMG" \
                bash -c "samtools dict /ref/GRCh38.fa -o /ref/GRCh38.dict" 2>&1 || {
                log_warn "Docker samtools dict failed. You can create it later with:"
                echo "    samtools dict ${ref_dir}/GRCh38.fa -o ${ref_dir}/GRCh38.dict"
            }
            if [ -f "${ref_dir}/GRCh38.dict" ]; then
                log_success "Sequence dictionary created"
            fi
        else
            log_warn "Docker not available. Create dictionary manually:"
            echo "    samtools dict ${ref_dir}/GRCh38.fa -o ${ref_dir}/GRCh38.dict"
        fi
    else
        log_skip "Sequence dictionary already exists"
    fi

    log_success "Reference genome setup complete"
}

merge_fastq_files() {
    print_section "STAGE 1c: FASTQ MERGE"

    local data_dir="${GENOMICS_DIR}/data/input/giab_hg002"
    local reads_dir="${data_dir}/reads"
    local input_dir="${GENOMICS_DIR}/data/input"

    # Check if merged files already exist
    if [ -f "${data_dir}/HG002_R1.fastq.gz" ] && [ -f "${data_dir}/HG002_R2.fastq.gz" ]; then
        local r1_size
        r1_size=$(stat -c%s "${data_dir}/HG002_R1.fastq.gz" 2>/dev/null || stat -f%z "${data_dir}/HG002_R1.fastq.gz" 2>/dev/null || echo 0)
        if [ "$r1_size" -gt 1000000 ]; then
            log_skip "Merged FASTQ files already exist ($(du -sh "${data_dir}/HG002_R1.fastq.gz" 2>/dev/null | awk '{print $1}') R1, $(du -sh "${data_dir}/HG002_R2.fastq.gz" 2>/dev/null | awk '{print $1}') R2)"
            return 0
        fi
    fi

    if [ "$DRY_RUN" = true ]; then
        log_info "Would merge 68 FASTQ chunks into 2 files"
        return 0
    fi

    # Check reads directory has files
    local file_count
    file_count=$(ls -1 "${reads_dir}"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
    if [ "$file_count" -lt 68 ]; then
        log_error "Expected 68 FASTQ files in ${reads_dir}, found $file_count"
        log_info "Run FASTQ download first: ./setup-data.sh --stage1"
        return 1
    fi

    local cores
    cores=$(nproc 2>/dev/null || echo 4)
    local compress_cmd="pigz -p $cores"
    if ! command -v pigz &>/dev/null; then
        compress_cmd="gzip"
        log_warn "pigz not found, using gzip (slower). Install pigz for faster merge."
    fi

    log_info "Merging FASTQ files using $cores CPU cores..."
    log_info "This takes 30-60 minutes depending on CPU speed"

    cd "$reads_dir"

    # Create sorted file lists
    ls -1 D1_S1_L001_R1_*.fastq.gz D1_S1_L002_R1_*.fastq.gz 2>/dev/null | sort > R1_files.txt
    ls -1 D1_S1_L001_R2_*.fastq.gz D1_S1_L002_R2_*.fastq.gz 2>/dev/null | sort > R2_files.txt

    local r1_count
    r1_count=$(wc -l < R1_files.txt | tr -d ' ')
    local r2_count
    r2_count=$(wc -l < R2_files.txt | tr -d ' ')
    log_info "Found $r1_count R1 files and $r2_count R2 files"

    # Merge R1
    log_info "Merging R1 files (this is the longest step)..."
    zcat $(cat R1_files.txt) | $compress_cmd > "${data_dir}/HG002_R1.fastq.gz"
    local r1_size_h
    r1_size_h=$(du -h "${data_dir}/HG002_R1.fastq.gz" | awk '{print $1}')
    log_success "R1 merge complete ($r1_size_h)"

    # Merge R2
    log_info "Merging R2 files..."
    zcat $(cat R2_files.txt) | $compress_cmd > "${data_dir}/HG002_R2.fastq.gz"
    local r2_size_h
    r2_size_h=$(du -h "${data_dir}/HG002_R2.fastq.gz" | awk '{print $1}')
    log_success "R2 merge complete ($r2_size_h)"

    cd "$SCRIPT_DIR"

    # Copy to pipeline input
    log_info "Copying merged files to pipeline input directory..."
    cp "${data_dir}/HG002_R1.fastq.gz" "${input_dir}/"
    cp "${data_dir}/HG002_R2.fastq.gz" "${input_dir}/"
    log_success "Merged FASTQ files ready in ${input_dir}/"
}

# ============================================================================
# STAGE 2: RAG/CHAT DATA
# ============================================================================

download_stage2() {
    print_section "STAGE 2: RAG/CHAT DATA"

    echo -e "  ${BOLD}Components:${NC}"
    echo "    - ClinVar variant_summary.txt.gz (~394 MB)"
    echo "    - ClinVar VCF + tabix index (~85 MB)"
    echo "    - AlphaMissense_hg38.tsv.gz (~614 MB)"
    echo ""

    local annot_dir="${RAG_CHAT_DIR}/data/annotations"

    if [ "$DRY_RUN" != true ]; then
        mkdir -p "$annot_dir" || {
            log_error "Cannot create directory: $annot_dir"
            return 1
        }
    fi

    local failed=0

    # ClinVar variant summary (used by ClinVarAnnotator)
    echo -e "  ${BOLD}ClinVar Variant Summary${NC}"
    if ! download_with_retry \
        "$CLINVAR_SUMMARY_URL" \
        "${annot_dir}/clinvar_variant_summary.txt.gz" \
        ""; then
        failed=$((failed + 1))
    fi

    echo ""

    # ClinVar VCF
    echo -e "  ${BOLD}ClinVar VCF (GRCh38)${NC}"
    if ! download_with_retry \
        "$CLINVAR_VCF_URL" \
        "${annot_dir}/clinvar.vcf.gz" \
        ""; then
        failed=$((failed + 1))
    fi

    # ClinVar VCF tabix index
    if ! download_with_retry \
        "$CLINVAR_VCF_TBI_URL" \
        "${annot_dir}/clinvar.vcf.gz.tbi" \
        ""; then
        # Non-critical: tabix index can be regenerated
        log_warn "ClinVar tabix index failed (non-critical, can be regenerated)"
    fi

    echo ""

    # AlphaMissense
    echo -e "  ${BOLD}AlphaMissense Predictions${NC}"
    if ! download_with_retry \
        "$ALPHAMISSENSE_URL" \
        "${annot_dir}/AlphaMissense_hg38.tsv.gz" \
        ""; then
        failed=$((failed + 1))
    fi

    echo ""

    # Embedding model info
    echo -e "  ${BOLD}Embedding Model${NC}"
    log_info "BGE-small-en-v1.5 (~130 MB) will auto-download on first pipeline run"
    log_info "No manual download required"

    echo ""

    if [ "$failed" -gt 0 ]; then
        log_error "$failed Stage 2 downloads failed"
        return 1
    fi

    local total_size
    total_size=$(du -sh "$annot_dir" 2>/dev/null | awk '{print $1}')
    log_success "Stage 2 data complete ($total_size)"
}

# ============================================================================
# STAGE 3: DRUG DISCOVERY DATA
# ============================================================================

download_stage3() {
    print_section "STAGE 3: DRUG DISCOVERY DATA"

    echo -e "  ${BOLD}Stage 3 data is primarily fetched on-demand during pipeline execution.${NC}"
    echo ""

    # Create directory structure
    local dirs=(
        "${DRUG_DISCOVERY_DIR}/data/structures/pdb_cache"
        "${DRUG_DISCOVERY_DIR}/data/structures/image_cache"
        "${DRUG_DISCOVERY_DIR}/data/molecules"
        "${DRUG_DISCOVERY_DIR}/data/targets"
    )

    if [ "$DRY_RUN" != true ]; then
        for dir in "${dirs[@]}"; do
            mkdir -p "$dir" 2>/dev/null || true
        done
        log_success "Directory structure created"
    else
        log_info "Would create directory structure under ${DRUG_DISCOVERY_DIR}/data/"
    fi

    # Optionally pre-fetch VCP demo PDB structures
    echo ""
    echo -e "  ${BOLD}VCP Demo Structures (optional)${NC}"

    local pdb_cache="${DRUG_DISCOVERY_DIR}/data/structures/pdb_cache"
    local fetched=0

    for pdb_id in $VCP_PDB_IDS; do
        local pdb_lower
        pdb_lower=$(echo "$pdb_id" | tr '[:upper:]' '[:lower:]')
        local pdb_file="${pdb_cache}/${pdb_lower}.pdb"

        if [ -f "$pdb_file" ]; then
            log_skip "${pdb_id}.pdb already cached"
        elif [ "$DRY_RUN" = true ]; then
            log_info "Would download: ${pdb_id}.pdb from RCSB"
        else
            log_info "Fetching ${pdb_id} from RCSB PDB..."
            wget -q -O "$pdb_file" "https://files.rcsb.org/download/${pdb_id}.pdb" 2>/dev/null || {
                log_warn "Could not fetch ${pdb_id}.pdb (non-critical, auto-fetched at runtime)"
                rm -f "$pdb_file"
                continue
            }
            fetched=$((fetched + 1))
            log_success "${pdb_id}.pdb cached"
        fi
    done

    echo ""
    echo -e "  ${BOLD}BioNeMo NIM Services${NC}"
    log_info "MolMIM and DiffDock containers are pulled on-demand via Docker"
    log_info "Requires NGC API key (set in .env)"
    log_info "Mock fallback available: set NIM_ALLOW_MOCK_FALLBACK=true"

    echo ""
    log_success "Stage 3 setup complete"
}

# ============================================================================
# STATUS DASHBOARD
# ============================================================================

show_status() {
    print_banner

    print_section "DATA STATUS DASHBOARD"

    local total_components=0
    local complete_components=0

    # Stage 1: Genomics
    echo -e "  ${BOLD}STAGE 1: GENOMICS${NC}"
    echo ""

    local reads_dir="${GENOMICS_DIR}/data/input/giab_hg002/reads"
    local ref_dir="${GENOMICS_DIR}/data/ref"
    local merge_dir="${GENOMICS_DIR}/data/input/giab_hg002"

    # R1 files
    total_components=$((total_components + 1))
    local r1_count=0
    if [ -d "$reads_dir" ]; then
        r1_count=$(ls -1 "${reads_dir}"/D1_S1_L*_R1_*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
    fi
    local r1_size
    r1_size=$(du -sh "$reads_dir" 2>/dev/null | awk '{print $1}' || echo "0")
    if [ "$r1_count" -ge 34 ]; then
        echo -e "    ${GREEN}[OK]${NC} FASTQ R1 files        ${r1_count}/34 files"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} FASTQ R1 files        ${r1_count}/34 files"
    fi

    # R2 files
    total_components=$((total_components + 1))
    local r2_count=0
    if [ -d "$reads_dir" ]; then
        r2_count=$(ls -1 "${reads_dir}"/D1_S1_L*_R2_*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [ "$r2_count" -ge 34 ]; then
        echo -e "    ${GREEN}[OK]${NC} FASTQ R2 files        ${r2_count}/34 files"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} FASTQ R2 files        ${r2_count}/34 files"
    fi

    # Merged files
    total_components=$((total_components + 1))
    if [ -f "${merge_dir}/HG002_R1.fastq.gz" ] && [ -f "${merge_dir}/HG002_R2.fastq.gz" ]; then
        local merged_size
        merged_size=$(du -sh "${merge_dir}/HG002_R1.fastq.gz" 2>/dev/null | awk '{print $1}')
        echo -e "    ${GREEN}[OK]${NC} Merged FASTQ           R1: ${merged_size}, R2: $(du -sh "${merge_dir}/HG002_R2.fastq.gz" 2>/dev/null | awk '{print $1}')"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} Merged FASTQ           Not created"
    fi

    # Reference genome
    total_components=$((total_components + 1))
    local ref_files=0
    for f in GRCh38.fa GRCh38.fa.fai GRCh38.fa.bwt GRCh38.fa.pac GRCh38.fa.sa; do
        [ -f "${ref_dir}/${f}" ] && ref_files=$((ref_files + 1))
    done
    if [ "$ref_files" -ge 5 ]; then
        local ref_size
        ref_size=$(du -sh "${ref_dir}/GRCh38.fa" 2>/dev/null | awk '{print $1}')
        echo -e "    ${GREEN}[OK]${NC} Reference genome       GRCh38 + BWA index (${ref_size})"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} Reference genome       ${ref_files}/5 index files"
    fi

    echo ""

    # Stage 2: RAG/Chat
    echo -e "  ${BOLD}STAGE 2: RAG/CHAT${NC}"
    echo ""

    local annot_dir="${RAG_CHAT_DIR}/data/annotations"

    # ClinVar
    total_components=$((total_components + 1))
    if [ -f "${annot_dir}/clinvar_variant_summary.txt.gz" ]; then
        local cv_size
        cv_size=$(du -sh "${annot_dir}/clinvar_variant_summary.txt.gz" 2>/dev/null | awk '{print $1}')
        echo -e "    ${GREEN}[OK]${NC} ClinVar variant_sum   ${cv_size}"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} ClinVar variant_sum   Not downloaded"
    fi

    # AlphaMissense
    total_components=$((total_components + 1))
    if [ -f "${annot_dir}/AlphaMissense_hg38.tsv.gz" ]; then
        local am_size
        am_size=$(du -sh "${annot_dir}/AlphaMissense_hg38.tsv.gz" 2>/dev/null | awk '{print $1}')
        echo -e "    ${GREEN}[OK]${NC} AlphaMissense          ${am_size}"
        complete_components=$((complete_components + 1))
    else
        echo -e "    ${RED}[!!]${NC} AlphaMissense          Not downloaded"
    fi

    echo ""

    # Stage 3: Drug Discovery
    echo -e "  ${BOLD}STAGE 3: DRUG DISCOVERY${NC}"
    echo ""

    local pdb_cache="${DRUG_DISCOVERY_DIR}/data/structures/pdb_cache"
    total_components=$((total_components + 1))
    local pdb_count=0
    if [ -d "$pdb_cache" ]; then
        pdb_count=$(ls -1 "${pdb_cache}"/*.pdb 2>/dev/null | wc -l | tr -d ' ')
    fi
    echo -e "    ${GREEN}[OK]${NC} PDB structures         ${pdb_count} cached (auto-fetched on demand)"
    complete_components=$((complete_components + 1))

    echo -e "    ${CYAN}[--]${NC} BioNeMo NIMs           Pulled via Docker on demand"

    echo ""

    # Overall
    echo -e "  ${NVIDIA_GREEN}================================================================================${NC}"
    local pct=$((complete_components * 100 / total_components))
    if [ "$complete_components" -eq "$total_components" ]; then
        echo -e "  ${GREEN}${BOLD}  OVERALL: ${complete_components}/${total_components} components complete (${pct}%) -- READY TO RUN${NC}"
    else
        local missing=$((total_components - complete_components))
        echo -e "  ${YELLOW}${BOLD}  OVERALL: ${complete_components}/${total_components} components complete (${pct}%) -- ${missing} remaining${NC}"
    fi
    echo -e "  ${NVIDIA_GREEN}================================================================================${NC}"
    echo ""
}

# ============================================================================
# VERIFY ALL
# ============================================================================

verify_all() {
    print_banner
    print_section "FULL DATA VERIFICATION"

    log_info "Re-verifying all downloaded files (this may take a while for large FASTQ files)..."
    echo ""

    local total_checked=0
    local total_passed=0
    local total_failed=0

    # Verify FASTQ files
    local reads_dir="${GENOMICS_DIR}/data/input/giab_hg002/reads"
    local index_file="${GENOMICS_DIR}/data/input/giab_hg002/hg002_illumina_2x250.index.tsv"

    if [ -f "$index_file" ] && [ -d "$reads_dir" ]; then
        echo -e "  ${BOLD}Verifying FASTQ files (MD5 checksums)${NC}"

        local urls_r1
        urls_r1=$(awk -F'\t' 'NR>1 && NF>3{print $1}' "$index_file")
        local md5_r1
        md5_r1=$(awk -F'\t' 'NR>1 && NF>3{print $2}' "$index_file")

        while IFS= read -r url && IFS= read -r md5 <&3; do
            local filename
            filename=$(basename "$url")
            local filepath="${reads_dir}/${filename}"
            total_checked=$((total_checked + 1))

            if [ -f "$filepath" ]; then
                local actual_md5
                actual_md5=$(md5sum "$filepath" 2>/dev/null | awk '{print $1}')
                if [ "$actual_md5" = "$md5" ]; then
                    total_passed=$((total_passed + 1))
                else
                    log_error "$filename: MD5 mismatch"
                    total_failed=$((total_failed + 1))
                fi
            else
                log_error "$filename: MISSING"
                total_failed=$((total_failed + 1))
            fi
            print_progress "$total_checked" 68 "FASTQ"
        done <<< "$urls_r1" 3<<< "$md5_r1"

        # Also verify R2
        local urls_r2
        urls_r2=$(awk -F'\t' 'NR>1 && NF>3{print $3}' "$index_file")
        local md5_r2
        md5_r2=$(awk -F'\t' 'NR>1 && NF>3{print $4}' "$index_file")

        while IFS= read -r url && IFS= read -r md5 <&3; do
            local filename
            filename=$(basename "$url")
            local filepath="${reads_dir}/${filename}"
            total_checked=$((total_checked + 1))

            if [ -f "$filepath" ]; then
                local actual_md5
                actual_md5=$(md5sum "$filepath" 2>/dev/null | awk '{print $1}')
                if [ "$actual_md5" = "$md5" ]; then
                    total_passed=$((total_passed + 1))
                else
                    log_error "$filename: MD5 mismatch"
                    total_failed=$((total_failed + 1))
                fi
            else
                log_error "$filename: MISSING"
                total_failed=$((total_failed + 1))
            fi
            print_progress "$total_checked" 68 "FASTQ"
        done <<< "$urls_r2" 3<<< "$md5_r2"

        echo ""
    fi

    # Verify Stage 2 files (gzip integrity)
    echo ""
    echo -e "  ${BOLD}Verifying Stage 2 files (gzip integrity)${NC}"
    local annot_dir="${RAG_CHAT_DIR}/data/annotations"

    for f in clinvar_variant_summary.txt.gz AlphaMissense_hg38.tsv.gz clinvar.vcf.gz; do
        total_checked=$((total_checked + 1))
        if [ -f "${annot_dir}/${f}" ]; then
            if gzip -t "${annot_dir}/${f}" 2>/dev/null; then
                log_success "$f: valid"
                total_passed=$((total_passed + 1))
            else
                log_error "$f: corrupt"
                total_failed=$((total_failed + 1))
            fi
        else
            log_error "$f: MISSING"
            total_failed=$((total_failed + 1))
        fi
    done

    echo ""
    echo -e "  ${BOLD}Verification Summary${NC}"
    echo "    Checked: $total_checked"
    echo "    Passed:  $total_passed"
    echo "    Failed:  $total_failed"
    echo ""

    if [ "$total_failed" -eq 0 ]; then
        log_success "All files verified successfully"
    else
        log_error "$total_failed files failed verification"
        echo ""
        echo -e "  ${YELLOW}Re-download failed files:${NC}"
        echo "    ./setup-data.sh --all"
    fi
}

# ============================================================================
# ERROR DIAGNOSIS
# ============================================================================

diagnose_failure() {
    local component="$1"
    local error_type="$2"
    local details="$3"

    echo ""
    echo -e "${RED}================================================================================${NC}"
    echo -e "${RED}  TROUBLESHOOTING: $component${NC}"
    echo -e "${RED}================================================================================${NC}"
    echo ""

    case "$error_type" in
        md5_mismatch)
            echo "  Problem: Files downloaded but checksum verification failed ($details)."
            echo "  Likely cause: Partial download, network corruption, or FTP mirror issue."
            echo ""
            echo "  Try:"
            echo "    1. Re-run this script (idempotent -- retries only failed files)"
            echo "    2. Use fewer connections: ./setup-data.sh --stage1 --connections 4"
            echo "    3. The NCBI FTP server sometimes throttles aggressive downloaders."
            echo "       Wait 15-30 minutes and try again."
            echo "    4. If failures persist, try at off-peak hours (US Eastern: 2-6 AM)"
            ;;
        connection_refused)
            echo "  Problem: Cannot connect to download server."
            echo "  Likely cause: Firewall, DNS resolution, or server maintenance."
            echo ""
            echo "  Try:"
            echo "    1. Check DNS:     nslookup ftp-trace.ncbi.nlm.nih.gov"
            echo "    2. Check HTTPS:   curl -v https://ftp-trace.ncbi.nlm.nih.gov/"
            echo "    3. NCBI status:   https://www.ncbi.nlm.nih.gov/Status/"
            echo "    4. Wait 15 minutes and retry (NCBI rate-limits aggressive downloaders)"
            ;;
        disk_full)
            local avail
            avail=$(df -h "$SCRIPT_DIR" 2>/dev/null | tail -1 | awk '{print $4}')
            echo "  Problem: Insufficient disk space."
            echo "  Available: $avail"
            echo ""
            echo "  Try:"
            echo "    1. Free space: sudo apt-get clean && docker system prune"
            echo "    2. Download stages individually: ./setup-data.sh --stage2 (only ~2 GB)"
            echo "    3. Symlink data directories to a larger volume:"
            echo "       ln -s /mnt/large-drive/genomics-data genomics-pipeline/data"
            ;;
    esac

    echo ""
    echo "  For more help: See docs/DATA_SETUP.md"
    echo "  Report issues: https://github.com/NVIDIA/hcls-ai-factory/issues"
    echo -e "${RED}================================================================================${NC}"
    echo ""
}

# ============================================================================
# HELP
# ============================================================================

show_help() {
    print_banner

    echo "Usage: ./setup-data.sh [command] [options]"
    echo ""
    echo "Downloads and verifies all data required for the HCLS AI Factory."
    echo "Total data: ~500 GB across three pipeline stages."
    echo ""
    echo "Commands:"
    echo "  --all           Download all data for all stages (default)"
    echo "  --stage1        Stage 1 only: HG002 FASTQ + reference genome (~300 GB)"
    echo "  --stage2        Stage 2 only: ClinVar + AlphaMissense (~2 GB)"
    echo "  --stage3        Stage 3 only: PDB structures + directory setup"
    echo "  --verify        Verify all downloaded data (no downloads)"
    echo "  --status        Show download status dashboard"
    echo "  --help          Show this help message"
    echo ""
    echo "Options:"
    echo "  --skip-merge    Skip FASTQ merge step (Stage 1)"
    echo "  --dry-run       Show what would be downloaded without downloading"
    echo "  --retry N       Max retry attempts per file (default: 3)"
    echo "  --connections N Parallel connections per file (default: 8)"
    echo ""
    echo "Examples:"
    echo "  ./setup-data.sh --all                    # Download everything"
    echo "  ./setup-data.sh --stage2                  # Just ClinVar + AlphaMissense (quick)"
    echo "  ./setup-data.sh --stage1 --connections 4  # Conservative FASTQ download"
    echo "  ./setup-data.sh --status                  # Check what's downloaded"
    echo "  ./setup-data.sh --verify                  # Re-verify all checksums"
    echo "  ./setup-data.sh --all --dry-run           # Preview without downloading"
    echo ""
    echo "Data breakdown:"
    echo "  Stage 1 (Genomics):       ~300 GB  (HG002 FASTQ + GRCh38 reference)"
    echo "  Stage 2 (RAG/Chat):         ~2 GB  (ClinVar + AlphaMissense)"
    echo "  Stage 3 (Drug Discovery): On-demand (PDB auto-fetched at runtime)"
    echo ""
    echo "Time estimates (depends on connection speed):"
    echo "  Stage 1 FASTQ download:   2-6 hours"
    echo "  Stage 1 reference genome: 5-20 minutes"
    echo "  Stage 1 FASTQ merge:      30-60 minutes"
    echo "  Stage 2:                  5-15 minutes"
    echo ""
    echo "For troubleshooting: See docs/DATA_SETUP.md"
    echo ""
}

# ============================================================================
# MAIN
# ============================================================================

main() {
    local command="--all"
    local stages=""

    # Parse arguments
    while [ $# -gt 0 ]; do
        case "$1" in
            --all)       command="--all"; shift ;;
            --stage1)    command="--stage"; stages="${stages}1"; shift ;;
            --stage2)    command="--stage"; stages="${stages}2"; shift ;;
            --stage3)    command="--stage"; stages="${stages}3"; shift ;;
            --verify)    command="--verify"; shift ;;
            --status)    command="--status"; shift ;;
            --help|-h)   command="--help"; shift ;;
            --skip-merge) SKIP_MERGE=true; shift ;;
            --dry-run)   DRY_RUN=true; shift ;;
            --retry)     MAX_RETRIES="$2"; shift 2 ;;
            --connections) CONNECTIONS="$2"; shift 2 ;;
            *)
                echo "Unknown option: $1"
                echo "Run ./setup-data.sh --help for usage"
                exit 1
                ;;
        esac
    done

    case "$command" in
        --help)
            show_help
            ;;
        --status)
            show_status
            ;;
        --verify)
            print_banner
            verify_all
            ;;
        --all)
            print_banner
            stages="123"
            preflight_checks "$stages"
            download_stage1
            download_stage2
            download_stage3
            echo ""
            show_status
            ;;
        --stage)
            print_banner

            if [ -z "$stages" ]; then
                echo "No stage specified. Use --stage1, --stage2, or --stage3."
                exit 1
            fi

            preflight_checks "$stages"

            if echo "$stages" | grep -q "1"; then download_stage1; fi
            if echo "$stages" | grep -q "2"; then download_stage2; fi
            if echo "$stages" | grep -q "3"; then download_stage3; fi

            echo ""
            show_status
            ;;
    esac
}

main "$@"
