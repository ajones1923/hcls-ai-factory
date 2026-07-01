#!/usr/bin/env python3
# =============================================================================
# Cardiology Intelligence Agent - Live Data Ingestion
# Author: Adam Jones
# Date: March 2026
# =============================================================================
"""Run live data ingestion from external sources into Milvus collections.

This script fetches, parses, embeds, and inserts cardiovascular research data
from PubMed and ClinicalTrials.gov into the agent's vector store, enabling
retrieval-augmented generation over the latest evidence base.

Data sources:
    PubMed          Cardiovascular research abstracts via NCBI E-utilities API.
                    Default query covers heart failure, coronary artery disease,
                    cardiomyopathy, arrhythmia, valvular disease, and structural
                    heart disease publications from the last 24 months.

    ClinicalTrials  Active and recruiting cardiovascular trials from
                    ClinicalTrials.gov v2 API, including NCTID, phase, endpoints,
                    enrollment, and sponsor information.

Usage:
    python scripts/run_ingest.py                           # ingest all sources
    python scripts/run_ingest.py --source pubmed           # PubMed only
    python scripts/run_ingest.py --source trials           # trials only
    python scripts/run_ingest.py --max-results 500         # limit result count
    python scripts/run_ingest.py --source pubmed --query "SGLT2 inhibitor heart failure"
    python scripts/run_ingest.py --dry-run                 # preview without insert
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.ingest.pubmed_parser import PubMedCardioParser
from src.ingest.clinical_trials_parser import ClinicalTrialsCardioParser

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Default search queries
# ---------------------------------------------------------------------------
DEFAULT_PUBMED_QUERIES = [
    "cardiovascular disease treatment",
    "heart failure management",
    "coronary artery disease",
    "atrial fibrillation",
    "cardiomyopathy genetics",
    "valvular heart disease",
    "cardiac imaging",
    "interventional cardiology",
    "SGLT2 inhibitor cardiovascular",
    "GLP-1 receptor agonist cardiovascular",
]

DEFAULT_TRIAL_CONDITIONS = [
    "Heart Failure",
    "Coronary Artery Disease",
    "Atrial Fibrillation",
    "Hypertrophic Cardiomyopathy",
    "Aortic Stenosis",
    "Pulmonary Hypertension",
]


# ---------------------------------------------------------------------------
# Ingest functions
# ---------------------------------------------------------------------------
def run_pubmed_ingest(
    max_results: int = 3000,
    query: Optional[str] = None,
    dry_run: bool = False,
) -> List[dict]:
    """Fetch and ingest PubMed cardiovascular literature.

    Parameters
    ----------
    max_results : int
        Maximum number of abstracts to retrieve per query.
    query : str, optional
        Custom PubMed query string. If None, uses DEFAULT_PUBMED_QUERIES.
    dry_run : bool
        If True, fetch and parse but do not insert into Milvus.

    Returns
    -------
    list[dict]
        Parsed PubMed records.
    """
    parser = PubMedCardioParser()
    queries = [query] if query else DEFAULT_PUBMED_QUERIES
    all_records = []

    for q in queries:
        logger.info("PubMed query: %s (max=%d)", q, max_results)
        records = parser.run(query=q, max_results=max_results, dry_run=dry_run)
        all_records.extend(records)
        logger.info("  Retrieved %d records for query: %s", len(records), q)

    # Deduplicate by PMID
    seen = set()
    unique = []
    for r in all_records:
        pmid = r.get("pmid", "")
        if pmid and pmid not in seen:
            seen.add(pmid)
            unique.append(r)

    logger.info(
        "PubMed ingest complete: %d total, %d unique records (dry_run=%s)",
        len(all_records),
        len(unique),
        dry_run,
    )
    return unique


def run_trials_ingest(
    max_results: int = 500,
    condition: Optional[str] = None,
    dry_run: bool = False,
) -> List[dict]:
    """Fetch and ingest ClinicalTrials.gov cardiovascular trials.

    Parameters
    ----------
    max_results : int
        Maximum number of trials to retrieve per condition.
    condition : str, optional
        Specific condition to query. If None, uses DEFAULT_TRIAL_CONDITIONS.
    dry_run : bool
        If True, fetch and parse but do not insert into Milvus.

    Returns
    -------
    list[dict]
        Parsed clinical trial records.
    """
    parser = ClinicalTrialsCardioParser()
    conditions = [condition] if condition else DEFAULT_TRIAL_CONDITIONS
    all_records = []

    for cond in conditions:
        logger.info("ClinicalTrials query: %s (max=%d)", cond, max_results)
        records = parser.run(condition=cond, max_results=max_results, dry_run=dry_run)
        all_records.extend(records)
        logger.info("  Retrieved %d trials for condition: %s", len(records), cond)

    # Deduplicate by NCT ID
    seen = set()
    unique = []
    for r in all_records:
        nct_id = r.get("nct_id", "")
        if nct_id and nct_id not in seen:
            seen.add(nct_id)
            unique.append(r)

    logger.info(
        "Clinical trials ingest complete: %d total, %d unique trials (dry_run=%s)",
        len(all_records),
        len(unique),
        dry_run,
    )
    return unique


def run_all(
    max_results: int = 3000,
    dry_run: bool = False,
) -> dict:
    """Run full data ingest from all external sources.

    Returns
    -------
    dict
        Summary with counts per source.
    """
    logger.info("Starting full data ingest...")
    start = time.time()

    pubmed_records = run_pubmed_ingest(max_results=max_results, dry_run=dry_run)
    trial_records = run_trials_ingest(
        max_results=min(max_results, 500), dry_run=dry_run
    )

    elapsed = time.time() - start
    summary = {
        "pubmed": len(pubmed_records),
        "trials": len(trial_records),
        "total": len(pubmed_records) + len(trial_records),
        "elapsed_seconds": round(elapsed, 1),
    }
    logger.info(
        "Full data ingest complete: %d PubMed + %d trials = %d total (%.1fs)",
        summary["pubmed"],
        summary["trials"],
        summary["total"],
        summary["elapsed_seconds"],
    )
    return summary


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Run Cardiology Intelligence Agent data ingestion",
    )
    parser.add_argument(
        "--source",
        choices=["pubmed", "trials", "all"],
        default="all",
        help="Data source to ingest (default: all)",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=3000,
        help="Maximum results per query (default: 3000)",
    )
    parser.add_argument(
        "--query",
        type=str,
        default=None,
        help="Custom PubMed search query (pubmed source only)",
    )
    parser.add_argument(
        "--condition",
        type=str,
        default=None,
        help="Specific condition to search (trials source only)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch and parse data without inserting into Milvus",
    )
    args = parser.parse_args()

    if args.source == "pubmed":
        run_pubmed_ingest(
            max_results=args.max_results,
            query=args.query,
            dry_run=args.dry_run,
        )
    elif args.source == "trials":
        run_trials_ingest(
            max_results=args.max_results,
            condition=args.condition,
            dry_run=args.dry_run,
        )
    else:
        run_all(max_results=args.max_results, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
