#!/usr/bin/env python3
"""Seed Milvus collections with data from the knowledge graph module.

Extracts structured knowledge from src/knowledge.py (ACTIONABLE_TARGETS,
THERAPY_MAP, RESISTANCE_MECHANISMS, PATHWAY_MAP, BIOMARKERS) and seeds
the corresponding Milvus collections with embedded records.

Usage: python3 scripts/seed_knowledge.py
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager
from src.knowledge import (
    ACTIONABLE_TARGETS,
    THERAPY_MAP,
    PATHWAY_MAP,
)


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        if isinstance(texts, str):
            texts = [texts]
        return self.model.encode(texts).tolist()


def seed_variants_from_knowledge(manager, embedder):
    """Seed onco_variants with ACTIONABLE_TARGETS entries."""
    records = []
    for gene, info in ACTIONABLE_TARGETS.items():
        description = info.get("description", "")
        cancer_types = ", ".join(info.get("cancer_types", []))
        key_variants = ", ".join(info.get("key_variants", []))
        therapies = ", ".join(info.get("targeted_therapies", []))

        text = (
            f"{gene} ({info.get('full_name', '')}). "
            f"Actionable in: {cancer_types}. "
            f"Key variants: {key_variants}. "
            f"Targeted therapies: {therapies}. "
            f"{description}"
        )
        embedding = embedder.encode(text)[0]

        records.append({
            "id": f"kg_target_{gene.lower()}",
            "embedding": embedding,
            "gene": gene,
            "variant_name": key_variants[:100] if key_variants else "",
            "variant_type": "snv",
            "cancer_type": cancer_types[:50] if cancer_types else "",
            "evidence_level": info.get("evidence_level", "B"),
            "drugs": therapies[:500] if therapies else "",
            "civic_id": "",
            "vrs_id": "",
            "text_summary": text[:3000],
            "clinical_significance": f"Actionable target in {cancer_types}"[:200],
            "allele_frequency": 0.0,
        })

    if records:
        count = manager.insert_batch("onco_variants", records)
        print(f"  Inserted {count} variant records from knowledge graph")
    return len(records)


def seed_therapies_from_knowledge(manager, embedder):
    """Seed onco_therapies with THERAPY_MAP entries."""
    records = []
    for drug, info in THERAPY_MAP.items():
        description = info.get("description", "")
        targets = ", ".join(info.get("targets", []))
        indications = ", ".join(info.get("indications", []))
        moa = info.get("mechanism", "")

        text = (
            f"{drug} ({info.get('category', 'targeted')} therapy). "
            f"Targets: {targets}. "
            f"Indications: {indications}. "
            f"Mechanism: {moa}. "
            f"{description}"
        )
        embedding = embedder.encode(text)[0]

        records.append({
            "id": f"kg_therapy_{drug.lower().replace(' ', '_').replace('+', '_')}",
            "embedding": embedding,
            "drug_name": drug[:200],
            "category": info.get("category", "targeted")[:30],
            "targets": targets[:200],
            "approved_indications": indications[:500],
            "resistance_mechanisms": "",
            "evidence_level": info.get("evidence_level", "A")[:20],
            "text_summary": text[:3000],
            "mechanism_of_action": moa[:500] if moa else "",
        })

    if records:
        count = manager.insert_batch("onco_therapies", records)
        print(f"  Inserted {count} therapy records from knowledge graph")
    return len(records)


def seed_pathways_from_knowledge(manager, embedder):
    """Seed onco_pathways with PATHWAY_MAP entries."""
    records = []
    for pathway_name, info in PATHWAY_MAP.items():
        description = info.get("description", "")
        genes = ", ".join(info.get("key_genes", []))
        targets = ", ".join(info.get("therapeutic_targets", []))
        cross_talk = ", ".join(info.get("cross_talk", []))

        text = (
            f"{pathway_name} signaling pathway. "
            f"Key genes: {genes}. "
            f"Therapeutic targets: {targets}. "
            f"Cross-talk: {cross_talk}. "
            f"{description}"
        )
        embedding = embedder.encode(text)[0]

        records.append({
            "id": f"kg_pathway_{pathway_name.lower().replace('/', '_').replace(' ', '_')}",
            "embedding": embedding,
            "name": pathway_name[:100],
            "key_genes": genes[:500],
            "therapeutic_targets": targets[:300],
            "cross_talk": cross_talk[:500],
            "text_summary": text[:3000],
        })

    if records:
        count = manager.insert_batch("onco_pathways", records)
        print(f"  Inserted {count} pathway records from knowledge graph")
    return len(records)


def main():
    print("=" * 60)
    print("Oncology Intelligence -- Knowledge Graph Seeder")
    print("=" * 60)

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()

    # Ensure collections exist
    for col in ["onco_variants", "onco_therapies", "onco_pathways"]:
        try:
            manager.get_collection(col)
        except Exception:
            manager.create_collection(col)

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Seeding from knowledge graph...")
    total = 0

    print("\n  --- Actionable Targets -> onco_variants ---")
    total += seed_variants_from_knowledge(manager, embedder)

    print("\n  --- Therapy Profiles -> onco_therapies ---")
    total += seed_therapies_from_knowledge(manager, embedder)

    print("\n  --- Pathways -> onco_pathways ---")
    total += seed_pathways_from_knowledge(manager, embedder)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {total} total records from knowledge graph")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
