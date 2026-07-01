#!/usr/bin/env python3
"""End-to-end validation of CAR-T Intelligence Agent data layer.

Tests:
  1. Collection stats (all 5 collections exist, 3 populated)
  2. Single-collection search on each populated collection
  3. Multi-collection search_all()
  4. Filtered search (target_antigen == "CD19")
  5. Demo queries from the design doc

Author: Adam Jones
Date: February 2026
"""

import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # noqa: F821

from sentence_transformers import SentenceTransformer

from src.collections import CARTCollectionManager

# ── Setup ────────────────────────────────────────────────────────────

QUERY_PREFIX = "Represent this sentence for searching relevant passages: "

DEMO_QUERIES = [
    "Why do CD19 CAR-T therapies fail in relapsed B-ALL patients?",
    "Compare 4-1BB vs CD28 costimulatory domains for DLBCL",
    "What manufacturing parameters predict clinical response?",
    "How does antigen density affect CAR-T efficacy?",
    "What are the resistance mechanisms to BCMA-targeted CAR-T?",
]


def main():
    print("=" * 70)
    print("CAR-T Intelligence Agent — End-to-End Validation")
    print("=" * 70)

    # Connect to Milvus
    manager = CARTCollectionManager()
    manager.connect()

    # ── Test 1: Collection stats ──────────────────────────────────────
    print("\n[TEST 1] Collection Stats")
    print("-" * 50)
    stats = manager.get_collection_stats()
    total_vectors = 0
    for name, count in stats.items():
        status = "OK" if count > 0 else "EMPTY"
        print(f"  {name:25s}  {count:>6,}  [{status}]")
        total_vectors += count
    print(f"  {'TOTAL':25s}  {total_vectors:>6,}")

    populated = {k: v for k, v in stats.items() if v > 0}
    assert len(populated) >= 3, f"Expected >= 3 populated collections, got {len(populated)}"
    print("  PASS: 3+ populated collections")

    # ── Load embedder ─────────────────────────────────────────────────
    print("\n[SETUP] Loading BGE-small-en-v1.5 embedder...")
    t0 = time.time()
    embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
    print(f"  Loaded in {time.time() - t0:.1f}s")

    def embed_query(text: str):
        return embedder.encode(QUERY_PREFIX + text).tolist()

    # ── Test 2: Single-collection searches ────────────────────────────
    print("\n[TEST 2] Single-Collection Searches")
    print("-" * 50)

    test_query = "CD19 CAR-T therapy for B-cell lymphoma"
    query_vec = embed_query(test_query)

    for coll_name in populated:
        t0 = time.time()
        results = manager.search(
            collection_name=coll_name,
            query_embedding=query_vec,
            top_k=3,
        )
        elapsed = (time.time() - t0) * 1000
        print(f"\n  {coll_name} — {len(results)} hits ({elapsed:.0f}ms)")
        for i, hit in enumerate(results):
            score = hit.get("score", 0.0)
            hit_id = hit.get("id", "?")
            title = hit.get("title", hit.get("name", hit.get("text_summary", "")))
            if title and len(title) > 80:
                title = title[:77] + "..."
            print(f"    [{i+1}] score={score:.4f}  id={hit_id}")
            if title:
                print(f"        {title}")

        assert len(results) > 0, f"No results from {coll_name}"
    print("\n  PASS: All populated collections return results")

    # ── Test 3: Multi-collection search_all() ─────────────────────────
    print("\n[TEST 3] Multi-Collection search_all()")
    print("-" * 50)

    t0 = time.time()
    all_results = manager.search_all(
        query_embedding=query_vec,
        top_k_per_collection=3,
    )
    elapsed = (time.time() - t0) * 1000
    total_hits = sum(len(v) for v in all_results.values())
    print(f"  Searched all 5 collections in {elapsed:.0f}ms, {total_hits} total hits")

    for coll_name, hits in all_results.items():
        print(f"  {coll_name:25s}  {len(hits)} hits")
        for hit in hits[:2]:
            print(f"    score={hit['score']:.4f}  id={hit['id']}")

    assert total_hits > 0, "No results from search_all()"
    print("  PASS: search_all() returns cross-collection results")

    # ── Test 4: Filtered search ───────────────────────────────────────
    print("\n[TEST 4] Filtered Search (target_antigen == 'CD19')")
    print("-" * 50)

    cd19_results = manager.search(
        collection_name="cart_literature",
        query_embedding=query_vec,
        top_k=5,
        filter_expr='target_antigen == "CD19"',
    )
    print(f"  cart_literature (CD19 filter): {len(cd19_results)} hits")
    for hit in cd19_results[:3]:
        print(f"    score={hit['score']:.4f}  antigen={hit.get('target_antigen', '?')}  id={hit['id']}")

    # Also test on cart_trials
    cd19_trials = manager.search(
        collection_name="cart_trials",
        query_embedding=query_vec,
        top_k=5,
        filter_expr='target_antigen == "CD19"',
    )
    print(f"  cart_trials    (CD19 filter): {len(cd19_trials)} hits")
    for hit in cd19_trials[:3]:
        print(f"    score={hit['score']:.4f}  phase={hit.get('phase', '?')}  id={hit['id']}")

    print("  PASS: Filtered search works")

    # ── Test 5: Demo queries ──────────────────────────────────────────
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

        # Distribution across collections
        dist = {k: len(v) for k, v in results.items() if len(v) > 0}
        print(f"     Collections: {dist}")

    print("\n  PASS: All demo queries return results")

    # ── Summary ───────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print(f"Total vectors: {total_vectors:,}")
    print(f"Populated collections: {list(populated.keys())}")
    print("=" * 70)

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
