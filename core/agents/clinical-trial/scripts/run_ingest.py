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
