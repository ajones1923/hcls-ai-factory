#!/usr/bin/env python3
"""End-to-end validation for the Oncology Intelligence Agent.

Performs comprehensive validation across all agent subsystems:
  1. Service health checks (Milvus, API)
  2. Collection stats and integrity
  3. Sample vector search queries
  4. Case creation workflow
  5. MTB packet generation
  6. Knowledge graph lookups

Usage: python3 scripts/validate_e2e.py [--verbose]
"""

import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import argparse


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

class ValidationResult:
    """Container for a single validation step result."""

    def __init__(self, name: str, passed: bool, message: str = "", duration: float = 0):
        self.name = name
        self.passed = passed
        self.message = message
        self.duration = duration

    def __str__(self):
        status = "PASS" if self.passed else "FAIL"
        time_str = f" ({self.duration:.2f}s)" if self.duration > 0 else ""
        msg = f" -- {self.message}" if self.message else ""
        return f"  [{status}] {self.name}{time_str}{msg}"


def timed_check(name, func):
    """Run a validation function and capture timing + exceptions."""
    t0 = time.time()
    try:
        passed, message = func()
        duration = time.time() - t0
        return ValidationResult(name, passed, message, duration)
    except Exception as e:
        duration = time.time() - t0
        return ValidationResult(name, False, f"Exception: {e}", duration)


# ---------------------------------------------------------------------------
# Validation steps
# ---------------------------------------------------------------------------

def check_milvus_connection():
    """Step 1: Verify Milvus connectivity."""
    from src.collections import OncoCollectionManager

    manager = OncoCollectionManager()
    manager.connect()
    manager.disconnect()
    return True, "Milvus connection successful"


def check_collection_stats():
    """Step 2: Check all collection stats."""
    from src.collections import COLLECTION_SCHEMAS, OncoCollectionManager

    manager = OncoCollectionManager()
    manager.connect()

    stats_summary = []
    total = 0
    empty_collections = []

    for name in sorted(COLLECTION_SCHEMAS.keys()):
        try:
            stats = manager.get_collection_stats(name)
            count = stats.get("num_entities", 0)
            total += count
            stats_summary.append(f"{name}={count}")
            if count == 0:
                empty_collections.append(name)
        except Exception:
            stats_summary.append(f"{name}=N/A")
            empty_collections.append(name)

    manager.disconnect()

    msg = f"{total} total records across {len(COLLECTION_SCHEMAS)} collections"
    if empty_collections:
        msg += f" (empty: {', '.join(empty_collections)})"

    # Pass if at least some collections have data
    return total > 0, msg


def check_embedding_model():
    """Step 3: Verify embedding model loads and produces correct dimensions."""
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    embedding = model.encode(["test oncology query"])
    dim = len(embedding[0])

    if dim != 384:
        return False, f"Expected dim=384, got dim={dim}"
    return True, f"BGE-small-en-v1.5 loaded, dim={dim}"


def check_vector_search():
    """Step 4: Run sample vector searches across collections."""
    from sentence_transformers import SentenceTransformer
    from src.collections import OncoCollectionManager

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    query_vector = model.encode(
        ["Represent this sentence for searching relevant passages: "
         "EGFR mutation treatment options in NSCLC"]
    ).tolist()[0]

    manager = OncoCollectionManager()
    manager.connect()

    collections_searched = 0
    total_hits = 0
    tested_collections = [
        "onco_variants", "onco_therapies", "onco_guidelines",
        "onco_biomarkers", "onco_trials",
    ]

    for col in tested_collections:
        try:
            hits = manager.search(name=col, query_vector=query_vector, top_k=3)
            total_hits += len(hits)
            collections_searched += 1
        except Exception:
            pass

    manager.disconnect()

    if collections_searched == 0:
        return False, "No collections could be searched"

    return total_hits > 0, f"{total_hits} hits from {collections_searched} collections"


def check_knowledge_graph():
    """Step 5: Verify knowledge graph lookups work."""
    from src.knowledge import (
        ACTIONABLE_TARGETS,
        THERAPY_PROFILES,
        RESISTANCE_MECHANISMS,
        PATHWAYS,
        BIOMARKERS,
        get_target_context,
    )

    # Check data structures exist and have content
    counts = {
        "targets": len(ACTIONABLE_TARGETS),
        "therapies": len(THERAPY_PROFILES),
        "resistance": len(RESISTANCE_MECHANISMS),
        "pathways": len(PATHWAYS),
        "biomarkers": len(BIOMARKERS),
    }

    total = sum(counts.values())
    if total == 0:
        return False, "Knowledge graph is empty"

    # Test a specific lookup
    egfr = ACTIONABLE_TARGETS.get("EGFR", {})
    if not egfr:
        return False, "EGFR target not found in knowledge graph"

    # Test context generation
    ctx = get_target_context("EGFR")
    if not ctx:
        return False, "get_target_context('EGFR') returned empty"

    counts_str = ", ".join(f"{k}={v}" for k, v in counts.items())
    return True, f"Knowledge graph OK ({counts_str})"


def check_case_creation():
    """Step 6: Verify case creation workflow (without Milvus insert)."""
    from src.models import CaseSnapshot

    # Create a synthetic case
    case = CaseSnapshot(
        id="test-case-001",
        patient_id="TEST-001",
        cancer_type="nsclc",
        stage="IV",
        variants="EGFR L858R, TP53 R175H",
        biomarkers="PD-L1 TPS 50%, TMB 8 mut/Mb",
        prior_therapies="none",
        text_summary="Stage IV NSCLC with EGFR L858R and TP53 R175H mutations. PD-L1 TPS 50%.",
    )

    if not case.id or not case.cancer_type:
        return False, "Case creation failed"

    return True, f"Case created: {case.id}, cancer_type={case.cancer_type}"


def check_seed_data_files():
    """Step 7: Verify all seed data JSON files exist and are valid."""
    seed_dir = PROJECT_ROOT / "data" / "reference"
    expected_files = [
        "variant_seed_data.json",
        "biomarker_seed_data.json",
        "therapy_seed_data.json",
        "pathway_seed_data.json",
        "guideline_seed_data.json",
        "trial_seed_data.json",
        "resistance_seed_data.json",
        "outcome_seed_data.json",
    ]

    missing = []
    record_counts = {}
    for fname in expected_files:
        fpath = seed_dir / fname
        if not fpath.exists():
            missing.append(fname)
            continue
        try:
            with open(fpath, "r") as f:
                data = json.load(f)
            record_counts[fname] = len(data) if isinstance(data, list) else 0
        except (json.JSONDecodeError, OSError):
            missing.append(f"{fname} (invalid JSON)")

    if missing:
        return False, f"Missing or invalid: {', '.join(missing)}"

    total = sum(record_counts.values())
    return True, f"All {len(expected_files)} seed files valid ({total} total records)"


def check_mtb_packet():
    """Step 8: Verify MTB packet model can be instantiated."""
    from src.models import MTBPacket

    packet = MTBPacket(
        case_id="test-case-001",
        patient_id="TEST-001",
        cancer_type="nsclc",
        stage="IV",
        actionable_variants=[],
        therapy_recommendations=[],
        trial_matches=[],
        resistance_alerts=[],
        open_questions=[],
    )

    if not packet.case_id:
        return False, "MTB packet creation failed"

    return True, f"MTB packet created: case_id={packet.case_id}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="End-to-end validation for Oncology Intelligence Agent"
    )
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    parser.parse_args()

    print("=" * 60)
    print("Oncology Intelligence -- End-to-End Validation")
    print("=" * 60)

    validations = [
        ("Milvus Connection", check_milvus_connection),
        ("Collection Stats", check_collection_stats),
        ("Embedding Model", check_embedding_model),
        ("Vector Search", check_vector_search),
        ("Knowledge Graph", check_knowledge_graph),
        ("Case Creation", check_case_creation),
        ("Seed Data Files", check_seed_data_files),
        ("MTB Packet", check_mtb_packet),
    ]

    results = []
    for name, func in validations:
        print(f"\n  Checking: {name}...")
        result = timed_check(name, func)
        results.append(result)
        print(result)

    # Summary
    passed = sum(1 for r in results if r.passed)
    failed = sum(1 for r in results if not r.passed)
    total_time = sum(r.duration for r in results)

    print(f"\n{'=' * 60}")
    print(f"Validation Summary: {passed}/{len(results)} checks passed")
    print(f"  Total time: {total_time:.2f}s")

    if failed > 0:
        print("\n  Failed checks:")
        for r in results:
            if not r.passed:
                print(f"    - {r.name}: {r.message}")

    status = "ALL PASSED" if failed == 0 else f"{failed} FAILED"
    print(f"\n  STATUS: {status}")
    print(f"{'=' * 60}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
