#!/bin/bash
set -e

echo "=========================================="
echo "Genomics Pipeline - Prerequisites Check"
echo "=========================================="
echo ""

# Check Docker
echo "[1/5] Checking Docker..."
if ! command -v docker &> /dev/null; then
    echo "❌ Docker not found. Please install Docker first."
    exit 1
fi
echo "✅ Docker version: $(docker --version)"
echo ""

# Check Docker daemon
echo "[2/5] Checking Docker daemon..."
if ! docker info &> /dev/null; then
    echo "❌ Docker daemon not running. Please start Docker."
    exit 1
fi
echo "✅ Docker daemon is running"
echo ""

# Check NVIDIA Container Runtime
echo "[3/5] Checking NVIDIA Container Runtime..."
if ! docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu22.04 nvidia-smi &> /dev/null; then
    echo "❌ NVIDIA Container Runtime not working."
    echo "   Please ensure nvidia-container-toolkit is installed and configured."
    exit 1
fi
echo "✅ NVIDIA Container Runtime is working"
echo ""

# Check GPU
echo "[4/5] Checking GPU availability..."
docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu22.04 nvidia-smi
echo ""

# Check disk space
echo "[5/5] Checking disk space..."
REQUIRED_GB=500
AVAILABLE_KB=$(df -k . | tail -1 | awk '{print $4}')
AVAILABLE_GB=$((AVAILABLE_KB / 1024 / 1024))

echo "Available disk space: ${AVAILABLE_GB} GB"
if [ "$AVAILABLE_GB" -lt "$REQUIRED_GB" ]; then
    echo "⚠️  Warning: Less than ${REQUIRED_GB} GB available."
    echo "   Full genome processing requires significant disk space."
    echo "   You may proceed for chr20 test only."
else
    echo "✅ Sufficient disk space available"
fi
echo ""

echo "=========================================="
echo "✅ Prerequisites check complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Authenticate with NGC: ./scripts/01-ngc-login.sh"
echo "2. Download data: ./scripts/02-download-data.sh"
echo "3. Setup reference: ./scripts/03-setup-reference.sh"
echo "4. Run chr20 test: ./scripts/04-run-chr20-test.sh"
echo "5. Run full genome: ./scripts/05-run-full-genome.sh"
