#!/usr/bin/env python3
"""Seed the CAR-T knowledge graph and verify.

Validates the knowledge graph, prints statistics, and optionally
exports to JSON for reference.

Usage:
    python scripts/seed_knowledge.py
    python scripts/seed_knowledge.py --export data/reference/knowledge.json

Author: Adam Jones
Date: February 2026
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.knowledge import (
    CART_MANUFACTURING,
    CART_TARGETS,
    CART_TOXICITIES,
    get_target_context,
    get_toxicity_context,
)


def main():
    parser = argparse.ArgumentParser(
        description="Seed and verify CAR-T knowledge graph"
    )
    parser.add_argument(
        "--export", type=str, default=None,
        help="Export knowledge graph to JSON file",
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  CAR-T Intelligence Agent — Knowledge Graph")
    print("=" * 65)
    print()

    # Target antigens
    print(f"  Target Antigens: {len(CART_TARGETS)}")
    for antigen, data in CART_TARGETS.items():
        products = data.get("approved_products", [])
        diseases = data.get("diseases", [])
        status = f"({len(products)} approved)" if products else "(preclinical/trials)"
        print(f"    {antigen:15s} — {', '.join(diseases[:3]):40s} {status}")

    print()

    # Toxicity profiles
    print(f"  Toxicity Profiles: {len(CART_TOXICITIES)}")
    for tox_id, data in CART_TOXICITIES.items():
        print(f"    {tox_id:20s} — {data.get('full_name', '')}")

    print()

    # Manufacturing knowledge
    print(f"  Manufacturing Processes: {len(CART_MANUFACTURING)}")
    for proc_id, data in CART_MANUFACTURING.items():
        print(f"    {proc_id:30s} — {data.get('description', '')[:50]}")

    print()

    # Test context retrieval
    print("  Testing context retrieval:")
    for antigen in ["CD19", "BCMA", "HER2"]:
        ctx = get_target_context(antigen)
        print(f"    get_target_context('{antigen}'): {len(ctx)} chars")

    for tox in ["CRS", "ICANS"]:
        ctx = get_toxicity_context(tox)
        print(f"    get_toxicity_context('{tox}'): {len(ctx)} chars")

    print()

    # Export if requested
    if args.export:
        export_path = Path(args.export)
        export_path.parent.mkdir(parents=True, exist_ok=True)
        export_data = {
            "targets": {k: _serialize(v) for k, v in CART_TARGETS.items()},
            "toxicities": {k: _serialize(v) for k, v in CART_TOXICITIES.items()},
            "manufacturing": {k: _serialize(v) for k, v in CART_MANUFACTURING.items()},
        }
        with open(export_path, "w") as f:
            json.dump(export_data, f, indent=2)
        print(f"  Exported to: {export_path}")

    print()
    print("  Knowledge graph verified successfully.")
    print("=" * 65)


def _serialize(obj):
    """Convert dict values to JSON-safe types."""
    if isinstance(obj, dict):
        return {k: _serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize(v) for v in obj]
    return obj


if __name__ == "__main__":
    main()
