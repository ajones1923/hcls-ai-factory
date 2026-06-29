"""CLI for running ingest pipelines for the Rare Disease Diagnostic Agent.

Usage:
    python scripts/run_ingest.py --source omim
    python scripts/run_ingest.py --source hpo
    python scripts/run_ingest.py --source orphanet
    python scripts/run_ingest.py --source gene_therapy
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

from src.ingest.omim_parser import OMIMParser
from src.ingest.hpo_parser import HPOParser
from src.ingest.orphanet_parser import OrphanetParser
from src.ingest.gene_therapy_parser import GeneTherapyParser

logger = logging.getLogger(__name__)


def run_omim(args: argparse.Namespace) -> None:
    """Run the OMIM ingest pipeline."""
    parser = OMIMParser(api_key=args.api_key)
    logger.info("Running OMIM ingest")
    records, stats = parser.run()

    logger.info(
        "OMIM ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)


def run_hpo(args: argparse.Namespace) -> None:
    """Run the HPO ingest pipeline."""
    parser = HPOParser()
    logger.info("Running HPO ingest")
    records, stats = parser.run()

    logger.info(
        "HPO ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)


def run_orphanet(args: argparse.Namespace) -> None:
    """Run the Orphanet ingest pipeline."""
    parser = OrphanetParser(api_key=args.api_key)
    logger.info("Running Orphanet ingest")
    records, stats = parser.run()

    logger.info(
        "Orphanet ingest complete: %d records validated "
        "(fetched=%d, parsed=%d, errors=%d, duration=%.1fs)",
        stats.total_validated, stats.total_fetched, stats.total_parsed,
        stats.total_errors, stats.duration_seconds,
    )

    if args.output:
        _write_output(records, args.output)


def run_gene_therapy(args: argparse.Namespace) -> None:
    """Run the Gene Therapy ingest pipeline."""
    parser = GeneTherapyParser()
    logger.info("Running Gene Therapy ingest")
    records, stats = parser.run()

    logger.info(
        "Gene Therapy ingest complete: %d records validated "
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
        description="Run Rare Disease Diagnostic Agent ingest pipelines"
    )
    parser.add_argument(
        "--source",
        choices=["omim", "hpo", "orphanet", "gene_therapy", "all"],
        required=True,
        help="Data source to ingest from",
    )
    parser.add_argument("--api-key", default=None, help="API key for the data source")
    parser.add_argument("--output", default=None, help="Output JSON file path")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.source == "omim":
        run_omim(args)
    elif args.source == "hpo":
        run_hpo(args)
    elif args.source == "orphanet":
        run_orphanet(args)
    elif args.source == "gene_therapy":
        run_gene_therapy(args)
    elif args.source == "all":
        run_omim(args)
        run_hpo(args)
        run_orphanet(args)
        run_gene_therapy(args)


if __name__ == "__main__":
    main()
