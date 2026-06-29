#!/usr/bin/env python3
"""Test the full Oncology Intelligence RAG pipeline.

Connects to Milvus, loads the embedding model, creates the RAG engine,
and runs 3 sample queries to validate the end-to-end pipeline:
  1. EGFR-mutant NSCLC treatment options
  2. BRAF V600E resistance mechanisms
  3. Biomarker testing for colorectal cancer

Usage: python3 scripts/test_rag_pipeline.py
"""

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, text):
        if isinstance(text, str):
            return self.model.encode([text]).tolist()[0]
        return self.model.encode(text).tolist()


SAMPLE_QUERIES = [
    {
        "question": "What are the first-line treatment options for EGFR L858R-mutant advanced NSCLC?",
        "collections": ["onco_variants", "onco_therapies", "onco_guidelines"],
        "expected_terms": ["osimertinib", "EGFR", "NSCLC"],
    },
    {
        "question": "What resistance mechanisms develop after BRAF V600E targeted therapy in melanoma?",
        "collections": ["onco_resistance", "onco_variants", "onco_pathways"],
        "expected_terms": ["BRAF", "resistance", "MEK"],
    },
    {
        "question": "What molecular biomarkers should be tested in metastatic colorectal cancer?",
        "collections": ["onco_biomarkers", "onco_guidelines", "onco_variants"],
        "expected_terms": ["KRAS", "MSI", "BRAF"],
    },
]


def run_query(manager, embedder, query_info, query_num):
    """Run a single RAG query and display results."""
    question = query_info["question"]
    collections = query_info["collections"]
    expected_terms = query_info["expected_terms"]

    print(f"\n  Query {query_num}: {question}")
    print(f"  Collections: {', '.join(collections)}")
    print(f"  {'─' * 50}")

    # Embed the query
    t0 = time.time()
    query_vector = embedder.encode(f"Represent this sentence for searching relevant passages: {question}")
    embed_time = time.time() - t0

    # Search each collection
    all_hits = []
    t0 = time.time()
    for col_name in collections:
        try:
            hits = manager.search(
                name=col_name,
                query_vector=query_vector,
                top_k=3,
            )
            all_hits.extend(hits)
        except Exception as e:
            print(f"  WARNING: Search failed for {col_name}: {e}")

    search_time = time.time() - t0

    # Sort by distance (lower = better for cosine)
    all_hits.sort(key=lambda h: h.get("_distance", float("inf")))

    # Display results
    print(f"  Embed time: {embed_time:.3f}s | Search time: {search_time:.3f}s")
    print(f"  Results: {len(all_hits)} hits across {len(collections)} collections\n")

    for i, hit in enumerate(all_hits[:5], 1):
        distance = hit.get("_distance", 0)
        collection = hit.get("_collection", "unknown")
        record_id = hit.get("id", "N/A")

        # Get text content
        text = (
            hit.get("text_summary", "")
            or hit.get("text_chunk", "")
            or hit.get("text", "")
            or ""
        )
        text_preview = text[:150].replace("\n", " ") + "..." if len(text) > 150 else text

        print(f"    {i}. [{collection}] (dist={distance:.4f}) id={record_id}")
        print(f"       {text_preview}")

    # Check for expected terms in results
    all_text = " ".join(
        str(hit.get("text_summary", "")) + str(hit.get("text", ""))
        for hit in all_hits[:5]
    ).lower()

    found = [t for t in expected_terms if t.lower() in all_text]
    missing = [t for t in expected_terms if t.lower() not in all_text]

    if found:
        print(f"\n  Expected terms found: {', '.join(found)}")
    if missing:
        print(f"  Expected terms MISSING: {', '.join(missing)}")

    return len(all_hits) > 0


def main():
    print("=" * 60)
    print("Oncology Intelligence -- RAG Pipeline Test")
    print("=" * 60)

    print("\n[1/4] Connecting to Milvus...")
    manager = OncoCollectionManager()
    try:
        manager.connect()
    except Exception as e:
        print(f"  ERROR: Cannot connect to Milvus: {e}")
        print("  Make sure Milvus is running on localhost:19530")
        return 1
    print("  Connected successfully.")

    print("\n[2/4] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()
    print("  Embedder loaded.")

    print("\n[3/4] Checking collection availability...")
    all_collections = set()
    for q in SAMPLE_QUERIES:
        all_collections.update(q["collections"])

    available = []
    for col in sorted(all_collections):
        try:
            stats = manager.get_collection_stats(col)
            count = stats.get("num_entities", 0)
            available.append(col)
            print(f"  {col:25s}  {count:>6,} records")
        except Exception:
            print(f"  {col:25s}  NOT FOUND (will skip)")

    if not available:
        print("\n  ERROR: No collections available. Run setup_collections.py --seed first.")
        manager.disconnect()
        return 1

    print(f"\n[4/4] Running {len(SAMPLE_QUERIES)} sample queries...")
    results = []
    for i, query_info in enumerate(SAMPLE_QUERIES, 1):
        # Only search available collections
        query_info["collections"] = [
            c for c in query_info["collections"] if c in available
        ]
        if not query_info["collections"]:
            print(f"\n  Query {i}: SKIPPED (no available collections)")
            results.append(False)
            continue

        success = run_query(manager, embedder, query_info, i)
        results.append(success)

    # Summary
    passed = sum(results)
    total = len(results)
    print(f"\n{'=' * 60}")
    print(f"RAG Pipeline Test Results: {passed}/{total} queries returned results")
    if passed == total:
        print("STATUS: ALL PASSED")
    elif passed > 0:
        print("STATUS: PARTIAL PASS (some queries returned no results)")
    else:
        print("STATUS: FAILED (no queries returned results)")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
