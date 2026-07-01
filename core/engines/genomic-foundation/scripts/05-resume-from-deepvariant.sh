#!/bin/bash

# ============================================================================
# DEPRECATED: This script is no longer needed.
#
# The main pipeline script (05-run-full-genome.sh) now automatically detects
# existing BAM/BAI files and resumes from DeepVariant when appropriate.
#
# Simply run:   ./run.sh full
#
# If BAM+BAI exist but no VCF, it will skip Steps 1-2 and resume from
# DeepVariant with retry logic (up to 3 attempts).
# ============================================================================

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  NOTE: This resume script is deprecated.                          ║"
echo "║                                                                    ║"
echo "║  The main pipeline now auto-detects existing BAM files and        ║"
echo "║  resumes from DeepVariant automatically.                          ║"
echo "║                                                                    ║"
echo "║  Forwarding to: ./run.sh full                                     ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/05-run-full-genome.sh" "$@"
