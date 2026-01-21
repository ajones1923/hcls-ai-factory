#!/bin/bash
# ================================================================
# HCLS AI Factory - Pre-commit Secret Scanner
# ================================================================
# Run this before pushing to public repo to catch accidental secrets
#
# Usage: ./scripts/check-secrets.sh
# ================================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "HCLS AI Factory - Secret Scanner"
echo "=========================================="
echo ""

FOUND_SECRETS=0

# Placeholder patterns to ignore (documentation examples)
PLACEHOLDER_PATTERNS="\.\.\.|\.\.\.|xxx|XXX|your[-_]|here|example|placeholder|REPLACE|change[-_]this"

# Function to filter out placeholder matches
filter_placeholders() {
    grep -v -i -E "$PLACEHOLDER_PATTERNS" || true
}

# Patterns to search for (real secrets only)
echo "Scanning for potential secrets..."
echo ""

# 1. Check for full Anthropic API keys (sk-ant-api03- followed by actual key chars)
echo -n "  Checking Anthropic API keys... "
results=$(grep -r -n -E "sk-ant-api[0-9]{2}-[A-Za-z0-9]{20,}" . \
    --include="*.py" --include="*.sh" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="*.example" --exclude="check-secrets.sh" \
    2>/dev/null | filter_placeholders || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

# 2. Check for NGC API keys (nvapi- followed by actual key)
echo -n "  Checking NGC API keys... "
results=$(grep -r -n -E "nvapi-[A-Za-z0-9_-]{20,}" . \
    --include="*.py" --include="*.sh" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="*.example" --exclude="check-secrets.sh" \
    2>/dev/null | filter_placeholders || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

# 3. Check for OpenAI-style keys
echo -n "  Checking OpenAI-style keys... "
results=$(grep -r -n -E "sk-[A-Za-z0-9]{32,}" . \
    --include="*.py" --include="*.sh" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="*.example" --exclude="check-secrets.sh" \
    2>/dev/null | filter_placeholders || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

# 4. Check for HuggingFace tokens (hf_ followed by actual chars)
echo -n "  Checking HuggingFace tokens... "
results=$(grep -r -n -E "hf_[A-Za-z0-9]{20,}" . \
    --include="*.py" --include="*.sh" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="*.example" --exclude="check-secrets.sh" \
    2>/dev/null | filter_placeholders || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

# 5. Check for AWS keys
echo -n "  Checking AWS keys... "
results=$(grep -r -n -E "AKIA[0-9A-Z]{16}" . \
    --include="*.py" --include="*.sh" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="*.example" --exclude="check-secrets.sh" \
    2>/dev/null | filter_placeholders || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

# 6. Check for private keys
echo -n "  Checking private keys... "
results=$(grep -r -l "PRIVATE KEY" . \
    --include="*.py" --include="*.yml" --include="*.yaml" \
    --include="*.json" --include="*.env" --include="*.txt" --include="*.pem" \
    --exclude-dir=".git" --exclude-dir="venv" --exclude-dir="node_modules" \
    --exclude="check-secrets.sh" \
    2>/dev/null || true)
if [ -n "$results" ]; then
    echo -e "${RED}FOUND${NC}"
    echo "$results" | head -5
    FOUND_SECRETS=1
else
    echo -e "${GREEN}OK${NC}"
fi

echo ""

# Check for .env files that shouldn't exist
echo "Checking for .env files..."
env_files=$(find . -name ".env" -not -name ".env.example" -not -path "./.git/*" -not -path "./venv/*" 2>/dev/null || true)
if [ -n "$env_files" ]; then
    echo -e "${YELLOW}WARNING:${NC} Found .env files (ensure these are gitignored):"
    echo "$env_files"
    echo ""
fi

# Final result
echo "=========================================="
if [ $FOUND_SECRETS -eq 1 ]; then
    echo -e "${RED}SCAN FAILED: Potential secrets found!${NC}"
    echo "Review the files above and remove any real credentials."
    echo "Use environment variables or .env files (gitignored) instead."
    exit 1
else
    echo -e "${GREEN}SCAN PASSED: No secrets detected.${NC}"
    echo ""
    echo "Note: This scanner ignores placeholder examples like:"
    echo "  - sk-ant-...  nvapi-...  your-key-here  etc."
    echo ""
    echo "Always manually review before pushing to public repo."
    exit 0
fi
