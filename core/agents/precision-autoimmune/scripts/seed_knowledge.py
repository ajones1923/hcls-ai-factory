#!/usr/bin/env python3
"""Seed Milvus collections from the curated autoimmune knowledge base.

Why this file exists: this agent was the only one of the eight with no seeder. Its 13
collections were created empty and stayed that way, so it answered requests from a corpus of
zero -- the same failure the clinical-trial agent had, where a service returns 200 and nothing
in it is grounded.

The knowledge is already curated in src/knowledge.py (101 entries across five domains). This
maps each domain onto its collection, embeds the prose, and projects every row onto the
declared schema -- collections here are structured (allele/disease/odds_ratio, not a generic
text field) and dynamic fields are off, so an unexpected key aborts the whole batch.

    python scripts/seed_knowledge.py            # seed
    python scripts/seed_knowledge.py --dry-run  # show what would be written
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

EMBED_MODEL = "BAAI/bge-small-en-v1.5"


def _rows_hla(k) -> list[dict]:
    """HLA_DISEASE_ASSOCIATIONS is keyed by ALLELE; each entry names the disease."""
    out = []
    for allele, entries in k.HLA_DISEASE_ASSOCIATIONS.items():
        for e in entries:
            disease = str(e.get("disease", "")).replace("_", " ")
            out.append({
                "text_chunk": (f"HLA association: {allele} with {disease}. "
                               f"Odds ratio {e.get('odds_ratio', 'NR')}"
                               + (f" in {e['population']} populations" if e.get("population") else "")
                               + ". " + str(e.get("note", e.get("mechanism", ""))).strip()),
                "allele": allele, "disease": disease,
                "odds_ratio": float(e.get("odds_ratio") or 0.0),
                "population": str(e.get("population", "")), "pmid": str(e.get("pmid", "")),
                "mechanism": str(e.get("mechanism", e.get("note", ""))),
                "clinical_implication": str(e.get("clinical_implication", e.get("note", ""))),
            })
    return out


def _rows_autoantibody(k) -> list[dict]:
    out = []
    # keyed by ANTIBODY; each entry names the disease it associates with
    for ab, entries in k.AUTOANTIBODY_DISEASE_MAP.items():
        for e in entries:
            dis = str(e.get("disease", "")).replace("_", " ")
            out.append({
                "text_chunk": (f"Autoantibody {ab} in {dis}. "
                               f"Sensitivity {e.get('sensitivity', 'NR')}, "
                               f"specificity {e.get('specificity', 'NR')}. "
                               f"{e.get('clinical_significance', '')}"),
                "antibody_name": ab, "associated_diseases": dis,
                "sensitivity": float(e.get("sensitivity") or 0.0),
                "specificity": float(e.get("specificity") or 0.0),
                "pattern": e.get("pattern", ""),
                "clinical_significance": e.get("clinical_significance", ""),
                "interpretation_guide": e.get("interpretation", ""),
            })
    return out


def _rows_biologics(k) -> list[dict]:
    out = []
    for d in k.BIOLOGIC_THERAPIES:
        name = d.get("drug_name") or d.get("name", "")
        out.append({
            "text_chunk": (f"Biologic therapy {name} ({d.get('drug_class', '')}). "
                           f"Mechanism: {d.get('mechanism', '')}. "
                           f"Indications: {d.get('indications', d.get('indicated_diseases', ''))}. "
                           f"Monitoring: {d.get('monitoring_requirements', d.get('monitoring', ''))}"),
            "drug_name": name, "drug_class": d.get("drug_class", ""),
            "mechanism": d.get("mechanism", ""),
            "indicated_diseases": str(d.get("indications", d.get("indicated_diseases", ""))),
            "pgx_considerations": str(d.get("pgx", d.get("pgx_considerations", ""))),
            "contraindications": str(d.get("contraindications", "")),
            "monitoring": str(d.get("monitoring_requirements", d.get("monitoring", ""))),
            "evidence_level": str(d.get("evidence_level", "")),
        })
    return out


def _rows_activity(k) -> list[dict]:
    out = []
    for disease, d in k.DISEASE_ACTIVITY_THRESHOLDS.items():
        out.append({"text_chunk": f"Disease activity thresholds for {disease}: {d}",
                    "disease": disease, "instrument": str(d.get("instrument", "")),
                    "thresholds": str(d)})
    return out


def _rows_flares(k) -> list[dict]:
    out = []
    for disease, d in k.FLARE_BIOMARKER_PATTERNS.items():
        out.append({"text_chunk": f"Flare biomarker pattern for {disease}: {d}",
                    "disease": disease, "pattern": str(d)})
    return out


DOMAINS = [
    ("autoimmune_hla_associations", _rows_hla),
    ("autoimmune_autoantibody_panels", _rows_autoantibody),
    ("autoimmune_biologic_therapies", _rows_biologics),
    ("autoimmune_disease_activity", _rows_activity),
    ("autoimmune_flare_patterns", _rows_flares),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    import src.knowledge as k
    from config.settings import settings

    total = 0
    if args.dry_run:
        for coll, fn in DOMAINS:
            logger.info("  %-38s %d rows", coll, len(fn(k)))
        return 0

    from pymilvus import MilvusClient
    from sentence_transformers import SentenceTransformer

    client = MilvusClient(uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}")
    model = SentenceTransformer(EMBED_MODEL)

    for coll, fn in DOMAINS:
        rows = fn(k)
        if not rows:
            continue
        try:
            fields = client.describe_collection(coll).get("fields", [])
        except Exception as exc:
            logger.warning("  %s: not found (%s) — run setup_collections.py first", coll, exc)
            continue
        valid = {f["name"] for f in fields}
        numeric = {f["name"] for f in fields
                   if any(t in getattr(f.get("type", ""), "name", str(f.get("type", ""))).upper()
                          for t in ("INT", "FLOAT", "DOUBLE"))} - {"id"}

        vecs = model.encode([r["text_chunk"] for r in rows], show_progress_bar=False).tolist()
        # `id` here is a VARCHAR primary key WITHOUT auto_id, so it must be supplied.
        # A stable, content-derived id also makes reseeding idempotent (upsert by key)
        # rather than duplicating the corpus on every run.
        import hashlib
        payload = []
        for i, (r, v) in enumerate(zip(rows, vecs)):
            row = {kk: vv for kk, vv in r.items() if kk in valid}
            if "id" in valid:
                digest = hashlib.sha1(r["text_chunk"].encode()).hexdigest()[:24]
                row["id"] = f"{coll}-{digest}"
            for name in valid:
                if name in ("id", "embedding"):
                    continue    # id set above; embedding set below
                row.setdefault(name, 0 if name in numeric else "")
            for name in numeric:
                try:
                    row[name] = float(row.get(name) or 0.0)
                except (TypeError, ValueError):
                    row[name] = 0.0
            row["embedding"] = v
            payload.append(row)
        try:
            client.upsert(collection_name=coll, data=payload)
            client.flush(coll)
            client.load_collection(collection_name=coll)
            logger.info("  inserted %d into %s", len(payload), coll)
            total += len(payload)
        except Exception as exc:
            logger.warning("  %s insert failed: %s", coll, exc)

    logger.info("Seed complete: %d records", total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
