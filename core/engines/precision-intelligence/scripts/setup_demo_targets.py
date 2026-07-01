#!/usr/bin/env python3
"""
Setup demo targets for GTC presentation.
Pre-populates target hypotheses with known druggable targets.
"""
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.target_hypothesis import TargetHypothesis, TargetHypothesisManager
from config.settings import settings


# Demo targets - well-known druggable targets with PDB structures
DEMO_TARGETS = [
    {
        "gene": "EGFR",
        "protein": "Epidermal growth factor receptor",
        "uniprot_id": "P00533",
        "rationale": "Receptor tyrosine kinase with variants affecting the ATP-binding pocket. Well-established cancer target with multiple approved inhibitors (erlotinib, gefitinib, osimertinib).",
        "therapeutic_area": "Oncology",
        "mechanism": "Tyrosine kinase inhibition",
        "confidence": "high",
        "priority": 5,
        "pdb_ids": ["7SYE", "1M17", "4HJO"],
        "status": "validated",
    },
    {
        "gene": "BRAF",
        "protein": "Serine/threonine-protein kinase B-raf",
        "uniprot_id": "P15056",
        "rationale": "Key driver in MAPK signaling pathway. V600E mutation is highly druggable with vemurafenib and dabrafenib.",
        "therapeutic_area": "Oncology",
        "mechanism": "Serine/threonine kinase inhibition",
        "confidence": "high",
        "priority": 5,
        "pdb_ids": ["6P3D", "4MNE", "5CSW"],
        "status": "validated",
    },
    {
        "gene": "PIK3CA",
        "protein": "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha",
        "uniprot_id": "P42336",
        "rationale": "Lipid kinase in PI3K/AKT pathway. Hotspot mutations (E545K, H1047R) are common in breast cancer. Alpelisib approved for PIK3CA-mutant breast cancer.",
        "therapeutic_area": "Oncology",
        "mechanism": "PI3K inhibition",
        "confidence": "high",
        "priority": 4,
        "pdb_ids": ["4JPS", "7L1C"],
        "status": "hypothesis",
    },
    {
        "gene": "KRAS",
        "protein": "GTPase KRas",
        "uniprot_id": "P01116",
        "rationale": "Previously considered undruggable, now targetable with covalent inhibitors (sotorasib for G12C). Critical oncogene in multiple cancers.",
        "therapeutic_area": "Oncology",
        "mechanism": "Covalent GTPase inhibition",
        "confidence": "medium",
        "priority": 4,
        "pdb_ids": ["6OIM", "4EPR"],
        "status": "hypothesis",
    },
    {
        "gene": "JAK2",
        "protein": "Tyrosine-protein kinase JAK2",
        "uniprot_id": "O60674",
        "rationale": "Non-receptor tyrosine kinase involved in cytokine signaling. V617F mutation drives myeloproliferative neoplasms. Ruxolitinib approved for JAK2-mutant conditions.",
        "therapeutic_area": "Hematology",
        "mechanism": "JAK inhibition",
        "confidence": "high",
        "priority": 4,
        "pdb_ids": ["4IVA", "6VGL"],
        "status": "hypothesis",
    },
]


def setup_demo_targets():
    """Create demo target hypotheses."""
    print("Setting up demo targets for GTC presentation...")
    print("=" * 60)

    manager = TargetHypothesisManager(storage_dir=settings.DATA_DIR / "targets")

    # Clear existing targets if desired
    existing = manager.list_all()
    if existing:
        print(f"Found {len(existing)} existing targets.")
        response = input("Clear existing targets? (y/N): ")
        if response.lower() == 'y':
            for t in existing:
                manager.delete(t.id)
            print("Cleared existing targets.")

    # Add demo targets
    for target_data in DEMO_TARGETS:
        target = TargetHypothesis(**target_data)
        manager.add(target)
        print(f"  Added: {target.gene} (Priority: {target.priority}, PDB: {', '.join(target.pdb_ids)})")

    print("=" * 60)
    print(f"Created {len(DEMO_TARGETS)} demo targets.")

    # Show summary
    summary = manager.get_summary()
    print(f"\nSummary:")
    print(f"  Total: {summary['total']}")
    print(f"  High Priority: {summary['high_priority']}")
    print(f"  By Status: {summary['by_status']}")

    # Export for Phase 5
    output_file = manager.export_for_phase5()
    print(f"\nExported to: {output_file}")

    print("\nDemo targets ready!")
    print("These will appear in both the Streamlit Chat and Portal interfaces.")


if __name__ == "__main__":
    setup_demo_targets()
