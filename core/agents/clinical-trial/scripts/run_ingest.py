"""CLI for running ingest pipelines for the Clinical Trial Intelligence Agent.

Usage:
    python scripts/run_ingest.py --source clinicaltrials --max-results 100
    python scripts/run_ingest.py --source pubmed --query "breast cancer clinical trial"
    python scripts/run_ingest.py --source regulatory --drugs pembrolizumab osimertinib
    python scripts/run_ingest.py --source all

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.ingest.clinicaltrials_parser import ClinicalTrialsParser
from src.ingest.pubmed_parser import PubMedTrialParser
from src.ingest.regulatory_parser import RegulatoryParser

logger = logging.getLogger(__name__)


def _persist(records, default_collection: str) -> int:
    """Write IngestRecords to Milvus. Returns rows written.

    This script fetched, parsed and validated — and then dropped everything on the floor
    unless --output was given. It logged "60 records validated" and wrote nothing, which is
    why the trial corpus sat at 59 rows while a working ClinicalTrials.gov ingest existed.

    Each record carries its own `collection_name`; rows are grouped by it, projected onto the
    declared schema (these collections are structured and dynamic fields are off, so one
    unexpected key aborts the batch), and given a content-derived id so re-running the ingest
    updates rather than duplicates.
    """
    import hashlib
    from collections import defaultdict

    try:
        from pymilvus import MilvusClient
        from sentence_transformers import SentenceTransformer
        from config.settings import settings
    except Exception as exc:
        logger.warning("Cannot persist (%s) — install pymilvus/sentence-transformers", exc)
        return 0

    client = MilvusClient(uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}")
    model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    groups = defaultdict(list)
    for r in records:
        groups[getattr(r, "collection_name", None) or default_collection].append(r)

    written = 0
    for coll, rows in groups.items():
        try:
            fields = client.describe_collection(coll).get("fields", [])
        except Exception:
            logger.warning("  %s: no such collection — run setup_collections.py", coll)
            continue
        valid = {fl["name"] for fl in fields}
        def _tname(fl):
            t = fl.get("type", "")
            return getattr(t, "name", str(t)).upper()

        numeric = {fl["name"] for fl in fields
                   if ("INT" in _tname(fl) or "FLOAT" in _tname(fl))
                   and "VECTOR" not in _tname(fl)}
        # Only genuine ARRAY columns may receive a list. `phase` is a VARCHAR here and a
        # list value aborts the whole batch.
        arrays = {fl["name"] for fl in fields if "ARRAY" in _tname(fl)}
        text_field = next((c for c in ("text", "text_chunk", "text_content", "description")
                           if c in valid), None)

        vecs = model.encode([r.text for r in rows], show_progress_bar=False).tolist()
        payload = []
        for r, v in zip(rows, vecs):
            row = {}
            meta = getattr(r, "metadata", None) or {}
            if isinstance(meta, dict):
                row.update({k: val for k, val in meta.items() if k in valid})
            if text_field:
                row[text_field] = r.text[:8192]
            for name in valid:
                if name in ("id", "embedding"):
                    continue
                row.setdefault(name, 0 if name in numeric else "")
            for name in numeric:
                if name == "id":
                    continue          # set from the content hash below
                try:
                    row[name] = int(float(row.get(name) or 0))
                except (TypeError, ValueError):
                    row[name] = 0
            for k2, v2 in list(row.items()):
                if k2 in numeric:
                    continue
                if isinstance(v2, list) and k2 not in arrays:
                    row[k2] = ", ".join(str(x) for x in v2)[:4096]
                elif not isinstance(v2, (str, list)):
                    row[k2] = str(v2)[:4096]
                elif isinstance(v2, str):
                    row[k2] = v2[:4096]
            # The primary key type differs per collection: VARCHAR here, int64 there.
            # A content hash keeps re-ingest idempotent either way.
            if "id" in valid:
                digest = hashlib.sha1(f"{coll}:{r.text}".encode()).hexdigest()
                row["id"] = (int(digest[:15], 16) if "id" in numeric
                             else f"{coll}-{digest[:24]}")
            row["embedding"] = v
            payload.append(row)
        try:
            client.upsert(collection_name=coll, data=payload)
            client.flush(coll)
            client.load_collection(collection_name=coll)
            logger.info("  persisted %d rows into %s", len(payload), coll)
            written += len(payload)
        except Exception as exc:
            logger.warning("  %s: insert failed: %s", coll, str(exc)[:160])
    return written


def run_clinicaltrials(args: argparse.Namespace) -> None:
    """Run the ClinicalTrials.gov ingest pipeline."""
    parser = ClinicalTrialsParser(api_key=args.api_key)
    conditions = args.conditions or ["cancer", "heart failure", "alzheimer"]

    logger.info("Running ClinicalTrials.gov ingest for conditions: %s", conditions)
    records, stats = parser.run(
        conditions=conditions,
        max_results=args.max_results,
    )

    logger.info(
        "ClinicalTrials.gov ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)
    if not args.dry_run:
        logger.info("Persisted %d rows", _persist(records, "trial_protocols"))


def run_pubmed(args: argparse.Namespace) -> None:
    """Run the PubMed ingest pipeline."""
    parser = PubMedTrialParser(api_key=args.api_key)
    query = args.query or '"clinical trial"[Publication Type]'

    logger.info("Running PubMed ingest for query: %s", query[:100])
    records, stats = parser.run(
        query=query,
        max_results=args.max_results,
    )

    logger.info(
        "PubMed ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)


def run_regulatory(args: argparse.Namespace) -> None:
    """Run the regulatory ingest pipeline."""
    parser = RegulatoryParser(api_key=args.api_key)
    drug_names = args.drugs or []

    logger.info("Running regulatory ingest for drugs: %s", drug_names)
    records, stats = parser.run(
        drug_names=drug_names,
        max_results=args.max_results,
        include_milestones=True,
    )

    logger.info(
        "Regulatory ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)


def _write_output(records: list, output_path: str) -> None:
    """Write ingest records to a JSON file."""
    data = [r.to_dict() for r in records]
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    logger.info("Wrote %d records to %s", len(data), output_path)


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Run Clinical Trial Intelligence Agent ingest pipelines"
    )
    parser.add_argument(
        "--source",
        choices=["clinicaltrials", "pubmed", "regulatory", "all"],
        required=True,
        help="Data source to ingest from",
    )
    parser.add_argument("--max-results", type=int, default=100, help="Max results to fetch")
    parser.add_argument("--api-key", default=None, help="API key for the data source")
    parser.add_argument("--output", default=None, help="Output JSON file path")
    parser.add_argument("--dry-run", action="store_true",
                        help="Fetch and validate without writing to Milvus (the old behaviour)")
    parser.add_argument("--query", default=None, help="PubMed search query")
    parser.add_argument("--conditions", nargs="*", help="ClinicalTrials.gov conditions")
    parser.add_argument("--drugs", nargs="*", help="Drug names for regulatory search")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.source == "clinicaltrials":
        run_clinicaltrials(args)
    elif args.source == "pubmed":
        run_pubmed(args)
    elif args.source == "regulatory":
        run_regulatory(args)
    elif args.source == "all":
        run_clinicaltrials(args)
        run_pubmed(args)
        run_regulatory(args)


if __name__ == "__main__":
    main()
