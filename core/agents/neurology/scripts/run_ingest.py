#!/usr/bin/env python3
"""Run ingest pipelines for the Neurology Intelligence Agent.

Supports running individual or all ingest parsers with optional
Milvus insertion.

Usage:
    python scripts/run_ingest.py --source pubmed
    python scripts/run_ingest.py --source neuroimaging
    python scripts/run_ingest.py --source eeg
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


def run_parser(source: str):
    """Run a specific ingest parser by source name."""
    from src.ingest.pubmed_neuro_parser import PubMedNeuroParser
    from src.ingest.neuroimaging_parser import NeuroimagingParser
    from src.ingest.eeg_parser import EEGParser

    parser_map = {
        "pubmed": ("PubMed Neurology", PubMedNeuroParser),
        "neuroimaging": ("Neuroimaging Protocols", NeuroimagingParser),
        "eeg": ("EEG Patterns", EEGParser),
    }

    if source == "all":
        sources = list(parser_map.keys())
    elif source in parser_map:
        sources = [source]
    else:
        logger.error("Unknown source: %s. Valid: %s, all", source, ", ".join(parser_map.keys()))
        sys.exit(1)

    total_records = 0
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

    logger.info("Ingest complete: %d total validated records", total_records)


def main():
    """Parse arguments and run ingest."""
    parser = argparse.ArgumentParser(description="Run neurology ingest pipelines")
    parser.add_argument(
        "--source",
        default="all",
        choices=["pubmed", "neuroimaging", "eeg", "all"],
        help="Ingest source to run (default: all)",
    )
    args = parser.parse_args()
    run_parser(args.source)


if __name__ == "__main__":
    main()
