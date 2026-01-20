#!/bin/bash
# =============================================================================
# Cache Cryo-EM Structure Images for Offline Demo
# =============================================================================
# Downloads all structure images from RCSB/EMDB for offline use

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="${SCRIPT_DIR}/../data/structures/images"

echo "Creating image cache directory: $CACHE_DIR"
mkdir -p "$CACHE_DIR"

# VCP Structures from vcp_structures.json
declare -A IMAGES=(
    # RCSB PDB images
    ["8ooi_assembly"]="https://cdn.rcsb.org/images/structures/oo/8ooi/8ooi_assembly-1.jpeg"
    ["9dil_assembly"]="https://cdn.rcsb.org/images/structures/di/9dil/9dil_assembly-1.jpeg"
    ["7k56_assembly"]="https://cdn.rcsb.org/images/structures/k5/7k56/7k56_assembly-1.jpeg"
    ["5ftk_assembly"]="https://cdn.rcsb.org/images/structures/ft/5ftk/5ftk_assembly-1.jpeg"
    # EMDB images
    ["EMD-36540"]="https://www.ebi.ac.uk/emdb/images/entry/EMD-36540/400_36540.gif"
    ["EMD-19622"]="https://www.ebi.ac.uk/emdb/images/entry/EMD-19622/400_19622.gif"
    ["EMD-22890"]="https://www.ebi.ac.uk/emdb/images/entry/EMD-22890/400_22890.gif"
)

echo ""
echo "Downloading structure images..."
echo ""

for name in "${!IMAGES[@]}"; do
    url="${IMAGES[$name]}"
    ext="${url##*.}"
    output="$CACHE_DIR/${name}.${ext}"

    if [ -f "$output" ]; then
        echo "  [SKIP] $name (already cached)"
    else
        echo "  [GET]  $name"
        if curl -sL --fail "$url" -o "$output" 2>/dev/null; then
            echo "         -> Saved to $output"
        else
            echo "         -> FAILED to download"
        fi
    fi
done

echo ""
echo "=============================================="
echo "Image cache complete!"
echo "Location: $CACHE_DIR"
echo ""
ls -la "$CACHE_DIR"
echo "=============================================="
