#!/bin/bash

# Genomics Pipeline - Main Runner Script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

show_help() {
    cat << EOF
Genomics Pipeline - FASTQ → VCF using NVIDIA Parabricks

Usage: ./run.sh <command>

Commands:
  check          Check prerequisites (Docker, GPU, disk space)
  login          Authenticate with NGC
  download       Download GIAB HG002 data (~200GB)
  reference      Setup GRCh38 reference genome
  test           Run chr20 fast test (~5-20 min)
  full           Run full genome pipeline (~30-110 min)
  clean          Clean output files (keeps input data)
  clean-all      Clean everything including downloaded data
  help           Show this help message

Quick Start:
  ./run.sh check       # Verify prerequisites
  ./run.sh login       # NGC authentication
  ./run.sh download    # Download data (takes hours)
  ./run.sh reference   # Setup reference
  ./run.sh test        # Quick test on chr20
  ./run.sh full        # Full genome run

Examples:
  ./run.sh check && ./run.sh login
  ./run.sh test | tee test-run.log
  ./run.sh full

EOF
}

case "$1" in
    check)
        ./scripts/00-setup-check.sh
        ;;
    login)
        ./scripts/01-ngc-login.sh
        ;;
    download)
        ./scripts/02-download-data.sh
        ;;
    reference)
        ./scripts/03-setup-reference.sh
        ;;
    test)
        ./scripts/04-run-chr20-test.sh
        ;;
    full)
        ./scripts/05-run-full-genome.sh
        ;;
    clean)
        echo "Cleaning output files..."
        rm -rf data/output/*
        mkdir -p data/output/logs
        echo "✅ Output files cleaned"
        ;;
    clean-all)
        echo "⚠️  This will delete ALL data including downloads!"
        echo "Press Ctrl+C to cancel, or Enter to continue..."
        read
        rm -rf data/input/* data/ref/* data/output/*
        mkdir -p data/input data/ref data/output/logs
        echo "✅ All data cleaned"
        ;;
    help|"")
        show_help
        ;;
    *)
        echo "Unknown command: $1"
        echo ""
        show_help
        exit 1
        ;;
esac
