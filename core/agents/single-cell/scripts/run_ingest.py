#!/usr/bin/env python3
"""Run ingest pipelines for the Single-Cell Intelligence Agent.

Supports running individual or all ingest parsers with optional
Milvus insertion.

Usage:
    python scripts/run_ingest.py --source cellxgene
    python scripts/run_ingest.py --source markers
    python scripts/run_ingest.py --source tme
    python scripts/run_ingest.py --source all
"""

import argparse
import logging
import sys
from pathlib import Path

# Ensure project root on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

from hcls_common.ingest_persist import persist_records



DRY_RUN = "--dry-run" in sys.argv          # fetch + validate without writing (the old behaviour)


def _milvus_uri() -> str:
    try:
        from config.settings import settings
        return f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}"
    except Exception:
        import os
        return (f"http://{os.getenv('MILVUS_HOST', 'localhost')}:"
                f"{os.getenv('MILVUS_PORT', '19530')}")


def run_parser(source: str):
    """Run a specific ingest parser by source name."""
    from src.ingest.cellxgene_parser import CellxGeneParser
    from src.ingest.marker_parser import MarkerParser
    from src.ingest.tme_parser import TMEParser

    parser_map = {
        "cellxgene": ("CellxGene Cell Types", CellxGeneParser),
        "markers": ("Marker Genes", MarkerParser),
        "tme": ("TME Profiles", TMEParser),
    }

    if source == "all":
        sources = list(parser_map.keys())
    elif source in parser_map:
        sources = [source]
    else:
        logger.error("Unknown source: %s. Valid: %s, all", source, ", ".join(parser_map.keys()))
        sys.exit(1)

    total_records = 0
    total_written = 0
    for src in sources:
        name, parser_cls = parser_map[src]
        logger.info("Running %s ingest ...", name)
        parser = parser_cls()
        records, stats = parser.run()
        total_records += len(records)
        logger.info(
            "  %s: %d validated records in %.1fs",
            name, stats.total_validated, stats.duration_seconds,
        )
        # This used to stop here: fetched, validated, and wrote nothing. The log said
        # "N validated records" and the corpus never changed.
        if not DRY_RUN:
            total_written += persist_records(records, "sc_literature",
                                             milvus_uri=_milvus_uri())

    logger.info("Ingest complete: %d validated, %d persisted", total_records, total_written)


def main():
    """Parse arguments and run ingest."""
    parser = argparse.ArgumentParser(description="Run single-cell ingest pipelines")
    parser.add_argument(
        "--source",
        default="all",
        choices=["cellxgene", "markers", "tme", "all"],
        help="Ingest source to run (default: all)",
    )
    args = parser.parse_args()
    run_parser(args.source)


if __name__ == "__main__":
    main()
