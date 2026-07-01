#!/usr/bin/env python3
"""Smoke test the Imaging Intelligence Agent RAG pipeline.

Creates the full RAG stack (ImagingCollectionManager + SentenceTransformer
+ NIMServiceManager + ImagingRAGEngine) and runs 3 test queries to verify
the retrieval and synthesis pipeline works end-to-end.

Usage:
    python scripts/test_rag_pipeline.py
    python scripts/test_rag_pipeline.py --host 10.0.0.1 --port 19530

Author: Adam Jones
Date: February 2026
"""

import argparse
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from config.settings import settings
from src.collections import ImagingCollectionManager
from src.nim.service_manager import NIMServiceManager
from src.rag_engine import ImagingRAGEngine

TEST_QUERIES = [
    "What is Lung-RADS classification?",
    "CT head hemorrhage detection AI",
    "Compare CT vs MRI for brain imaging",
]


def main():
    parser = argparse.ArgumentParser(
        description="Smoke test the Imaging RAG pipeline"
    )
    parser.add_argument(
        "--host", default=None, help="Milvus host (default: localhost)"
    )
    parser.add_argument(
        "--port", type=int, default=None, help="Milvus port (default: 19530)"
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Top-K results per collection",
    )
    parser.add_argument(
        "--retrieve-only",
        action="store_true",
        help="Only test retrieval, skip LLM synthesis",
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  Imaging Intelligence Agent — RAG Pipeline Smoke Test")
    print("=" * 65)
    print()

    # --- Step 1: Initialize components ---
    logger.info("Loading embedding model...")
    from sentence_transformers import SentenceTransformer

    embedder = SentenceTransformer(settings.EMBEDDING_MODEL)

    logger.info("Connecting to Milvus...")
    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Show collection stats
    stats = manager.get_collection_stats()
    total_records = sum(stats.values())
    logger.info(f"Milvus connected: {total_records:,} total records across {len(stats)} collections")
    for name, count in stats.items():
        if count > 0:
            logger.info(f"  {name}: {count:,}")

    if total_records == 0:
        logger.warning(
            "No records in any collection. Run seed scripts first. "
            "Retrieval will return empty results."
        )

    logger.info("Initializing NIM Service Manager (mock fallback enabled)...")
    nim_manager = NIMServiceManager(settings)

    # Create LLM client for RAG synthesis
    llm_client = nim_manager.llm

    logger.info("Creating ImagingRAGEngine...")
    rag_engine = ImagingRAGEngine(
        collection_manager=manager,
        embedder=embedder,
        llm_client=llm_client,
        nim_service_manager=nim_manager,
    )

    # --- Step 2: Run test queries ---
    results_summary = []

    for i, query in enumerate(TEST_QUERIES, 1):
        print()
        print(f"  Query {i}/{len(TEST_QUERIES)}: {query}")
        print("  " + "-" * 50)

        start = time.time()

        try:
            # Retrieve evidence
            evidence = rag_engine.retrieve(query, top_k_per_collection=args.top_k)
            retrieval_time = (time.time() - start) * 1000

            hit_count = evidence.hit_count
            collections_with_hits = len(
                [c for c, hits in evidence.hits_by_collection().items() if hits]
            )

            logger.info(
                f"  Retrieved {hit_count} hits from {collections_with_hits} collections "
                f"in {retrieval_time:.0f}ms"
            )

            # Show top hits
            for hit in evidence.hits[:3]:
                logger.info(
                    f"    [{hit.collection}] score={hit.score:.3f}: "
                    f"{hit.text[:100]}..."
                )

            # Optionally run LLM synthesis
            answer = ""
            if not args.retrieve_only:
                try:
                    answer = rag_engine.query(query)
                    total_time = (time.time() - start) * 1000
                    logger.info(f"  LLM answer ({total_time:.0f}ms): {answer[:200]}...")
                except Exception as e:
                    logger.warning(f"  LLM synthesis failed (retrieval OK): {e}")
                    answer = f"[LLM error: {e}]"
                    total_time = (time.time() - start) * 1000
            else:
                total_time = retrieval_time

            results_summary.append({
                "query": query,
                "hits": hit_count,
                "collections": collections_with_hits,
                "time_ms": total_time,
                "status": "PASS" if hit_count > 0 or total_records == 0 else "WARN",
            })

        except Exception as e:
            logger.error(f"  Query failed: {e}")
            results_summary.append({
                "query": query,
                "hits": 0,
                "collections": 0,
                "time_ms": (time.time() - start) * 1000,
                "status": "FAIL",
            })

    # --- Step 3: Print summary ---
    print()
    print("=" * 65)
    print("  RAG Pipeline Smoke Test Results")
    print("=" * 65)
    for r in results_summary:
        print(
            f"  [{r['status']:4s}] {r['query'][:45]:45s} "
            f"hits={r['hits']:3d}  colls={r['collections']}  "
            f"time={r['time_ms']:.0f}ms"
        )

    passed = sum(1 for r in results_summary if r["status"] == "PASS")
    total = len(results_summary)
    print()
    print(f"  Result: {passed}/{total} queries passed")
    print("=" * 65)

    manager.disconnect()


if __name__ == "__main__":
    main()
