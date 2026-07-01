#!/usr/bin/env python3
"""
Setup Milvus collections for the Precision Autoimmune Agent.

Usage:
    python scripts/setup_collections.py                     # Create collections
    python scripts/setup_collections.py --drop-existing     # Recreate (drop first)
    python scripts/setup_collections.py --seed              # Create + seed knowledge
    python scripts/setup_collections.py --host milvus --port 19530
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from loguru import logger

from config.settings import settings
from src.collections import AutoimmuneCollectionManager


def seed_knowledge(cm: AutoimmuneCollectionManager, embedder):
    """Seed collections with static knowledge base data."""
    from src.knowledge import (
        AUTOANTIBODY_DISEASE_MAP,
        BIOLOGIC_THERAPIES,
        DISEASE_ACTIVITY_THRESHOLDS,
        FLARE_BIOMARKER_PATTERNS,
        HLA_DISEASE_ASSOCIATIONS,
    )

    # ── HLA Associations ──
    records = []
    for allele, assocs in HLA_DISEASE_ASSOCIATIONS.items():
        for i, a in enumerate(assocs):
            text = (
                f"{allele} is associated with {a['disease'].replace('_', ' ')} "
                f"with an odds ratio of {a['odds_ratio']}. "
                f"{a.get('note', '')}"
            )
            records.append({
                "id": f"hla_{allele}_{i}".replace("*", "").replace(":", ""),
                "text_chunk": text,
                "allele": allele,
                "disease": a["disease"],
                "odds_ratio": a["odds_ratio"],
                "population": "European",
                "pmid": a.get("pmid", ""),
                "mechanism": a.get("note", ""),
                "clinical_implication": f"Genetic susceptibility marker for {a['disease'].replace('_', ' ')}",
            })
    if records:
        logger.info(f"Seeding {len(records)} HLA associations...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_HLA_ASSOCIATIONS, records)

    # ── Autoantibody Panels ──
    records = []
    for ab_name, assocs in AUTOANTIBODY_DISEASE_MAP.items():
        diseases = ", ".join(a["disease"].replace("_", " ") for a in assocs)
        sens_spec = "; ".join(
            f"{a['disease']}: sens={a['sensitivity']}, spec={a['specificity']}"
            for a in assocs
        )
        text = (
            f"{ab_name} autoantibody is associated with: {diseases}. "
            f"Sensitivity/specificity: {sens_spec}. "
            f"{assocs[0].get('note', '')}"
        )
        records.append({
            "id": f"ab_{ab_name.replace('-', '_').replace('/', '_').replace(' ', '_')}",
            "text_chunk": text,
            "antibody_name": ab_name,
            "associated_diseases": diseases,
            "sensitivity": assocs[0]["sensitivity"],
            "specificity": assocs[0]["specificity"],
            "pattern": "",
            "clinical_significance": f"Key autoantibody marker for {diseases}",
            "interpretation_guide": sens_spec,
        })
    if records:
        logger.info(f"Seeding {len(records)} autoantibody panels...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_AUTOANTIBODY_PANELS, records)

    # ── Biologic Therapies ──
    records = []
    for t in BIOLOGIC_THERAPIES:
        diseases = ", ".join(d.replace("_", " ") for d in t["indicated_diseases"])
        pgx = "; ".join(t.get("pgx_considerations", []))
        text = (
            f"{t['drug_name']} ({t['drug_class']}): {t.get('mechanism', '')}. "
            f"Indicated for: {diseases}. "
            f"PGx: {pgx}. "
            f"Contraindications: {', '.join(t.get('contraindications', []))}."
        )
        records.append({
            "id": f"bio_{t['drug_name'].lower().replace(' ', '_')}",
            "text_chunk": text,
            "drug_name": t["drug_name"],
            "drug_class": t["drug_class"],
            "mechanism": t.get("mechanism", ""),
            "indicated_diseases": diseases,
            "pgx_considerations": pgx,
            "contraindications": ", ".join(t.get("contraindications", [])),
            "monitoring": ", ".join(t.get("monitoring_requirements", [])),
            "dosing": "",
            "evidence_level": "Level A",
        })
    if records:
        logger.info(f"Seeding {len(records)} biologic therapies...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_BIOLOGIC_THERAPIES, records)

    # ── Disease Activity Scoring ──
    records = []
    for score_name, info in DISEASE_ACTIVITY_THRESHOLDS.items():
        components = ", ".join(info.get("components", []))
        thresholds = str(info.get("thresholds", {}))
        text = (
            f"{score_name} is a disease activity score for "
            f"{info['disease'].replace('_', ' ')}. "
            f"Components: {components}. "
            f"Thresholds: {thresholds}. "
            f"Reference: {info.get('reference', 'N/A')}."
        )
        records.append({
            "id": f"act_{score_name.lower().replace('-', '_')}",
            "text_chunk": text,
            "score_name": score_name,
            "disease": info["disease"],
            "components": components,
            "thresholds": thresholds,
            "interpretation": "Higher scores indicate more active disease",
            "monitoring_frequency": "Every 3-6 months in active disease",
        })
    if records:
        logger.info(f"Seeding {len(records)} disease activity scores...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_DISEASE_ACTIVITY, records)

    # ── Flare Patterns ──
    records = []
    for disease, pattern in FLARE_BIOMARKER_PATTERNS.items():
        markers = ", ".join(pattern["early_warning_biomarkers"])
        protective = ", ".join(pattern.get("protective_signals", []))
        text = (
            f"Flare prediction for {disease.replace('_', ' ')}: "
            f"Early warning biomarkers: {markers}. "
            f"Protective signals: {protective}."
        )
        records.append({
            "id": f"flare_{disease}",
            "text_chunk": text,
            "disease": disease,
            "biomarker_pattern": markers,
            "early_warning_signs": markers,
            "typical_timeline": "2-4 weeks before clinical flare",
            "protective_factors": protective,
            "intervention_triggers": str(pattern.get("thresholds", {})),
        })
    if records:
        logger.info(f"Seeding {len(records)} flare patterns...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_FLARE_PATTERNS, records)

    # ── Disease Classification Criteria ──
    from src.diagnostic_engine import CLASSIFICATION_CRITERIA
    records = []
    for disease, criteria_set in CLASSIFICATION_CRITERIA.items():
        criteria_items = str(criteria_set.get("criteria", {}))
        text = (
            f"{criteria_set['name']} for {disease.value.replace('_', ' ')}: "
            f"Threshold: {criteria_set['threshold']} points. "
            f"{criteria_set.get('entry_criterion', '')}. "
            f"Criteria: {criteria_items[:2000]}"
        )
        records.append({
            "id": f"crit_{disease.value}",
            "text_chunk": text[:3000],
            "disease": disease.value,
            "criteria_set": criteria_set["name"],
            "criteria_type": "classification",
            "year": 2019,
            "required_score": str(criteria_set["threshold"]),
            "criteria_items": criteria_items[:3000],
            "sensitivity_specificity": "",
        })
    if records:
        logger.info(f"Seeding {len(records)} classification criteria...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_DISEASE_CRITERIA, records)

    # ── Cross-Disease Patterns ──
    from src.diagnostic_engine import OVERLAP_SYNDROMES
    records = []
    for syndrome, info in OVERLAP_SYNDROMES.items():
        diseases = [d.value if hasattr(d, "value") else str(d)
                    for d in info.get("diseases", info.get("features_from", []))]
        shared = info.get("shared_markers", info.get("required", []))
        text = (
            f"{syndrome.replace('_', ' ').title()}: "
            f"Involves {', '.join(d.replace('_', ' ') for d in diseases)}. "
            f"Shared markers: {', '.join(shared)}."
        )
        records.append({
            "id": f"cross_{syndrome}",
            "text_chunk": text,
            "primary_disease": diseases[0] if diseases else "",
            "associated_conditions": ", ".join(diseases),
            "shared_pathways": "",
            "shared_biomarkers": ", ".join(shared),
            "overlap_criteria": str(info),
            "co_occurrence_rate": 0.0,
        })
    if records:
        logger.info(f"Seeding {len(records)} cross-disease patterns...")
        texts = [r["text_chunk"] for r in records]
        embeddings = embedder.encode(texts, batch_size=32, show_progress_bar=False)
        for rec, emb in zip(records, embeddings):
            rec["embedding"] = emb.tolist()
        cm.insert_batch(settings.COLL_CROSS_DISEASE, records)

    logger.info("Knowledge seeding complete!")


def main():
    parser = argparse.ArgumentParser(description="Setup Autoimmune Agent Milvus collections")
    parser.add_argument("--host", default=settings.MILVUS_HOST)
    parser.add_argument("--port", type=int, default=settings.MILVUS_PORT)
    parser.add_argument("--drop-existing", action="store_true")
    parser.add_argument("--seed", action="store_true", help="Seed with knowledge base data")
    args = parser.parse_args()

    cm = AutoimmuneCollectionManager(host=args.host, port=args.port)
    cm.connect()

    logger.info(f"Creating collections (drop_existing={args.drop_existing})...")
    collections = cm.create_all_collections(drop_existing=args.drop_existing)
    logger.info(f"Created {len(collections)} collections: {list(collections.keys())}")

    if args.seed:
        logger.info("Loading embedding model for seeding...")
        from sentence_transformers import SentenceTransformer
        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
        seed_knowledge(cm, embedder)

    # Print stats
    stats = cm.get_collection_stats()
    for name, count in sorted(stats.items()):
        logger.info(f"  {name}: {count} records")

    cm.disconnect()
    logger.info("Setup complete!")


if __name__ == "__main__":
    main()
