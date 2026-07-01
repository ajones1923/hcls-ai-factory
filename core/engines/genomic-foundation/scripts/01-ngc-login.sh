#!/bin/bash
set -e

echo "=========================================="
echo "NGC Authentication Setup"
echo "=========================================="
echo ""
echo "You need an NGC API key to pull NVIDIA containers."
echo "Get your API key from: https://ngc.nvidia.com/setup/api-key"
echo ""
echo "When prompted:"
echo "  Username: \$oauthtoken"
echo "  Password: <your NGC API key>"
echo ""

docker login nvcr.io

echo ""
echo "✅ NGC authentication successful!"
echo ""
echo "Next step: ./scripts/02-download-data.sh"
