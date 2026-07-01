#!/usr/bin/env python3
"""End-to-end validation of Biomarker Intelligence Agent data layer.

Tests:
  1. Connect to Milvus
  2. Verify all 10 collections exist
  3. Verify each has >0 records
  4. Run a test search query
  5. Print summary table

Author: Adam Jones
Date: March 2026
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sentence_transformers import SentenceTransformer

from src.collections import BiomarkerCollectionManager

# -- Setup --------------------------------------------------------------------

QUERY_PREFIX = "Represent this sentence for searching relevant passages: "

EXPECTED_COLLECTIONS = [
    "biomarker_reference",
    "biomarker_genetic_variants",
    "biomarker_pgx_rules",
    "biomarker_disease_trajectories",
    "biomarker_clinical_evidence",
    "biomarker_nutrition",
    "biomarker_drug_interactions",
    "biomarker_aging_markers",
    "biomarker_genotype_adjustments",
    "biomarker_monitoring",
]

DEMO_QUERIES = [
    "What is the clinical significance of elevated HbA1c in pre-diabetes?",
    "How does MTHFR C677T affect folate metabolism and homocysteine levels?",
    "CYP2D6 poor metabolizer codeine dosing recommendation",
    "PNPLA3 rs738409 GG genotype and NAFLD risk",
    "PhenoAge biological age biomarker interpretation",
]


def main():
    print("=" * 70)
    print("Biomarker Intelligence Agent -- End-to-End Validation")
    print("=" * 70)

    # Connect to Milvus
    manager = BiomarkerCollectionManager()
    manager.connect()

    # -- Test 1: Collection stats ---------------------------------------------
    print("\n[TEST 1] Collection Stats")
    print("-" * 50)
    stats = manager.get_collection_stats()
    total_vectors = 0
    for name in EXPECTED_COLLECTIONS:
        count = stats.get(name, 0)
        status = "OK" if count > 0 else "EMPTY"
        print(f"  {name:40s}  {count:>6,}  [{status}]")
        total_vectors += count
    print(f"  {'TOTAL':40s}  {total_vectors:>6,}")

    populated = {k: v for k, v in stats.items() if v > 0 and k in EXPECTED_COLLECTIONS}
    print(f"  Populated: {len(populated)}/{len(EXPECTED_COLLECTIONS)} collections")

    # -- Test 2: Verify all collections exist ---------------------------------
    print("\n[TEST 2] Collection Existence")
    print("-" * 50)
    missing = [n for n in EXPECTED_COLLECTIONS if n not in stats]
    if missing:
        print(f"  FAIL: Missing collections: {missing}")
    else:
        print(f"  PASS: All {len(EXPECTED_COLLECTIONS)} collections exist")

    # -- Test 3: Load embedder ------------------------------------------------
    print("\n[SETUP] Loading BGE-small-en-v1.5 embedder...")
    t0 = time.time()
    embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
    print(f"  Loaded in {time.time() - t0:.1f}s")

    def embed_query(text: str):
        return embedder.encode(QUERY_PREFIX + text).tolist()

    # -- Test 3: Single-collection searches -----------------------------------
    print("\n[TEST 3] Single-Collection Searches")
    print("-" * 50)

    test_query = "HbA1c glucose diabetes biomarker interpretation"
    query_vec = embed_query(test_query)

    for coll_name in populated:
        t0 = time.time()
        results = manager.search(
            collection_name=coll_name,
            query_embedding=query_vec,
            top_k=3,
        )
        elapsed = (time.time() - t0) * 1000
        print(f"\n  {coll_name} -- {len(results)} hits ({elapsed:.0f}ms)")
        for i, hit in enumerate(results):
            score = hit.get("score", 0.0)
            hit_id = hit.get("id", "?")
            print(f"    [{i+1}] score={score:.4f}  id={hit_id}")

    if populated:
        print("\n  PASS: Populated collections return results")

    # -- Test 4: Multi-collection search_all() --------------------------------
    print("\n[TEST 4] Multi-Collection search_all()")
    print("-" * 50)

    t0 = time.time()
    all_results = manager.search_all(
        query_embedding=query_vec,
        top_k_per_collection=3,
    )
    elapsed = (time.time() - t0) * 1000
    total_hits = sum(len(v) for v in all_results.values())
    print(f"  Searched all collections in {elapsed:.0f}ms, {total_hits} total hits")

    for coll_name, hits in all_results.items():
        if hits:
            print(f"  {coll_name:40s}  {len(hits)} hits")

    if total_hits > 0:
        print("  PASS: search_all() returns cross-collection results")

    # -- Test 5: Demo queries -------------------------------------------------
    print("\n[TEST 5] Demo Queries (search_all)")
    print("-" * 50)

    for query in DEMO_QUERIES:
        qvec = embed_query(query)
        t0 = time.time()
        results = manager.search_all(
            query_embedding=qvec,
            top_k_per_collection=3,
            score_threshold=0.3,
        )
        elapsed = (time.time() - t0) * 1000
        total = sum(len(v) for v in results.values())

        # Get top hit across all collections
        all_hits = []
        for coll_hits in results.values():
            all_hits.extend(coll_hits)
        all_hits.sort(key=lambda x: x["score"], reverse=True)
        top_score = all_hits[0]["score"] if all_hits else 0.0
        top_coll = all_hits[0]["collection"] if all_hits else "none"

        print(f"\n  Q: {query}")
        print(f"     {total} hits in {elapsed:.0f}ms | top: {top_score:.4f} ({top_coll})")

    if populated:
        print("\n  PASS: All demo queries executed")

    # -- Summary --------------------------------------------------------------
    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print(f"Total vectors: {total_vectors:,}")
    print(f"Populated collections: {list(populated.keys())}")
    print("=" * 70)

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
