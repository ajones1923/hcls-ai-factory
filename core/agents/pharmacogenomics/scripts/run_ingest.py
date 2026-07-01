#!/usr/bin/env python3
"""Run all PGx data ingest parsers with progress tracking.

Executes each ingest pipeline (CPIC, PharmVar, PharmGKB, FDA labels, PubMed,
ClinicalTrials) in sequence, fetching live data from upstream sources and
storing embedded records into the corresponding Milvus collections.

Assumes collections have already been created by setup_collections.py.

Usage:
    python scripts/run_ingest.py
    python scripts/run_ingest.py --sources cpic pharmgkb
    python scripts/run_ingest.py --host milvus-standalone --max-results 100

Author: Adam Jones
Date: March 2026
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Dict

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger
from tqdm import tqdm

from src.collections import PGxCollectionManager


# ═══════════════════════════════════════════════════════════════════════
# INGEST PIPELINE REGISTRY
# ═══════════════════════════════════════════════════════════════════════

INGEST_SOURCES = {
    "cpic": {
        "label": "CPIC Guidelines",
        "module": "src.ingest.cpic_parser",
        "class": "CPICGuidelineParser",
        "collection": "pgx_drug_guidelines",
        "description": "Clinical Pharmacogenetics Implementation Consortium guidelines",
    },
    "pharmvar": {
        "label": "PharmVar Alleles",
        "module": "src.ingest.pharmvar_parser",
        "class": "PharmVarParser",
        "collection": "pgx_gene_reference",
        "description": "Pharmacogene Variation Consortium star allele definitions",
    },
    "pharmgkb": {
        "label": "PharmGKB Annotations",
        "module": "src.ingest.pharmgkb_parser",
        "class": "PharmGKBParser",
        "collection": "pgx_drug_interactions",
        "description": "PharmGKB drug-gene interaction annotations",
    },
    "fda": {
        "label": "FDA Drug Labels",
        "module": "src.ingest.fda_label_parser",
        "class": "FDALabelParser",
        "collection": "pgx_fda_labels",
        "description": "FDA pharmacogenomic labeling information",
    },
    "pubmed": {
        "label": "PubMed Literature",
        "module": "src.ingest.pubmed_parser",
        "class": "PubMedPGxParser",
        "collection": "pgx_clinical_evidence",
        "description": "PubMed pharmacogenomics publications",
    },
    "population": {
        "label": "Population Frequencies",
        "module": "src.ingest.population_parser",
        "class": "PopulationFrequencyParser",
        "collection": "pgx_population_data",
        "description": "Population-specific allele frequency data",
    },
}


def _import_pipeline_class(module_path: str, class_name: str):
    """Dynamically import an ingest pipeline class."""
    import importlib

    module = importlib.import_module(module_path)
    return getattr(module, class_name)


def run_single_source(
    source_key: str,
    source_config: dict,
    manager: PGxCollectionManager,
    model,
    max_results: int,
) -> Dict[str, int]:
    """Run a single ingest pipeline and return results.

    Returns:
        Dict with keys: fetched, parsed, inserted, errors.
    """
    result = {"fetched": 0, "parsed": 0, "inserted": 0, "errors": 0}

    try:
        PipelineClass = _import_pipeline_class(
            source_config["module"], source_config["class"]
        )
    except (ImportError, AttributeError) as e:
        logger.warning(f"  [{source_key}] Pipeline not yet implemented: {e}")
        result["errors"] = 1
        return result

    try:
        pipeline = PipelineClass(manager, model)
        count = pipeline.run(
            collection_name=source_config["collection"],
            max_results=max_results,  # passed as fetch_kwargs
        )
        result["inserted"] = count
        logger.info(f"  [{source_key}] Inserted {count} records")
    except Exception as e:
        logger.error(f"  [{source_key}] Pipeline failed: {e}")
        result["errors"] = 1

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run PGx data ingest pipelines"
    )
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    parser.add_argument(
        "--sources",
        nargs="+",
        default=None,
        choices=list(INGEST_SOURCES.keys()),
        help="Specific source(s) to ingest (default: all)",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=500,
        help="Max records to fetch per source (default: 500)",
    )
    args = parser.parse_args()

    logger.info("=" * 65)
    logger.info("  PGx Intelligence Agent — Data Ingest")
    logger.info("=" * 65)

    # ── Load embedding model ──
    try:
        from sentence_transformers import SentenceTransformer

        logger.info("Loading embedding model: BAAI/bge-small-en-v1.5")

        class SimpleEmbedder:
            """Wrapper around SentenceTransformer for pipeline compatibility."""

            def __init__(self, st_model):
                self._model = st_model

            def encode(self, texts):
                return self._model.encode(texts).tolist()

        st_model = SentenceTransformer("BAAI/bge-small-en-v1.5")
        embedder = SimpleEmbedder(st_model)
        logger.info("Embedding model loaded (384-dim)")
    except ImportError as e:
        logger.error(f"Cannot load embedding model: {e}")
        sys.exit(1)

    # ── Connect to Milvus ──
    manager = PGxCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # ── Determine sources ──
    sources = args.sources or list(INGEST_SOURCES.keys())
    logger.info(f"Running {len(sources)} ingest pipeline(s): {', '.join(sources)}")

    # ── Run pipelines ──
    start = time.time()
    all_results: Dict[str, Dict[str, int]] = {}

    progress = tqdm(sources, desc="Ingest pipelines", unit="source")
    for source_key in progress:
        config = INGEST_SOURCES[source_key]
        progress.set_postfix_str(config["label"])
        logger.info(f"  Starting: {config['label']} ({config['description']})")

        result = run_single_source(
            source_key, config, manager, embedder, args.max_results
        )
        all_results[source_key] = result

    elapsed = time.time() - start

    # ── Summary ──
    logger.info("")
    logger.info("=" * 65)
    logger.info("  Ingest Summary")
    logger.info("=" * 65)
    logger.info(f"  Time elapsed: {elapsed:.1f}s")
    logger.info("")
    logger.info(f"  {'Source':<20s} {'Inserted':>10s} {'Errors':>8s}")
    logger.info(f"  {'-' * 20} {'-' * 10} {'-' * 8}")

    total_inserted = 0
    total_errors = 0
    for source_key, result in all_results.items():
        label = INGEST_SOURCES[source_key]["label"]
        inserted = result["inserted"]
        errors = result["errors"]
        total_inserted += inserted
        total_errors += errors
        logger.info(f"  {label:<20s} {inserted:>10,} {errors:>8}")

    logger.info(f"  {'-' * 20} {'-' * 10} {'-' * 8}")
    logger.info(f"  {'TOTAL':<20s} {total_inserted:>10,} {total_errors:>8}")
    logger.info("=" * 65)

    # ── Final stats ──
    stats = manager.get_collection_stats()
    logger.info("")
    logger.info("Milvus collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    manager.disconnect()
    logger.info("Ingest complete.")


if __name__ == "__main__":
    main()
