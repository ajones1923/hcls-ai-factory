#!/usr/bin/env python3
"""Integration test: Full CAR-T RAG pipeline with Claude LLM.

Tests the complete flow:
  1. Embed query (BGE-small-en-v1.5)
  2. Search all Milvus collections
  3. Knowledge graph augmentation
  4. Query expansion
  5. Claude LLM response generation

Requires:
  - Milvus running on localhost:19530 with populated collections
  - ANTHROPIC_API_KEY set (or in rag-chat-pipeline/.env)

Author: Adam Jones
Date: February 2026
"""

import os
import sys
import time
from pathlib import Path

# Load API key from rag-chat-pipeline .env if not already set
if not os.environ.get("ANTHROPIC_API_KEY"):
    env_path = Path(os.environ.get("CART_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")) / ".env"
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            if line.startswith("ANTHROPIC_API_KEY="):
                os.environ["ANTHROPIC_API_KEY"] = line.split("=", 1)[1].strip().strip('"')
                break

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import anthropic
from sentence_transformers import SentenceTransformer

from src.collections import CARTCollectionManager
from src.rag_engine import CARTRAGEngine
from src import knowledge as kg
from src import query_expansion as qe


# ── Lightweight wrappers ──────────────────────────────────────────────

class SimpleEmbedder:
    """Wrapper around SentenceTransformer with embed_text() API."""
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def embed_text(self, text: str):
        return self.model.encode(text).tolist()

    def encode(self, texts):
        return self.model.encode(texts).tolist()


class SimpleLLMClient:
    """Wrapper around Anthropic SDK with generate() / generate_stream()."""
    def __init__(self, model: str = "claude-sonnet-4-20250514"):
        self.client = anthropic.Anthropic()
        self.model = model

    def generate(self, prompt, system_prompt="", max_tokens=2048, temperature=0.7):
        msg = self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}],
        )
        return msg.content[0].text

    def generate_stream(self, prompt, system_prompt="", max_tokens=2048, temperature=0.7):
        with self.client.messages.stream(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}],
        ) as stream:
            for text in stream.text_stream:
                yield text


DEMO_QUERIES = [
    "Why do CD19 CAR-T therapies fail in relapsed B-ALL patients?",
    "Compare 4-1BB vs CD28 costimulatory domains for DLBCL",
    "What are the resistance mechanisms to BCMA-targeted CAR-T?",
]


def main():
    print("=" * 70)
    print("CAR-T Intelligence Agent — Full RAG Pipeline Integration Test")
    print("=" * 70)

    # ── Step 1: Initialize components ─────────────────────────────────
    print("\n[1/5] Initializing components...")

    # Check API key
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("  ERROR: ANTHROPIC_API_KEY not set")
        return 1
    print("  ANTHROPIC_API_KEY: set")

    # Milvus
    manager = CARTCollectionManager()
    manager.connect()
    stats = manager.get_collection_stats()
    total = sum(stats.values())
    print(f"  Milvus: {total:,} vectors across {sum(1 for v in stats.values() if v > 0)} collections")

    # Embedder
    t0 = time.time()
    embedder = SimpleEmbedder()
    print(f"  Embedder: BGE-small-en-v1.5 loaded ({time.time()-t0:.1f}s)")

    # LLM
    llm = SimpleLLMClient()
    print(f"  LLM: {llm.model}")

    # ── Step 2: Build RAG engine ──────────────────────────────────────
    print("\n[2/5] Building CARTRAGEngine...")
    engine = CARTRAGEngine(
        collection_manager=manager,
        embedder=embedder,
        llm_client=llm,
        knowledge=kg,
        query_expander=qe,
    )
    print("  Engine ready with knowledge graph + query expansion")

    # ── Step 3: Test retrieve() only ──────────────────────────────────
    print("\n[3/5] Testing retrieve() (search + knowledge + expansion)...")
    from src.models import AgentQuery

    query = AgentQuery(question=DEMO_QUERIES[0])
    t0 = time.time()
    evidence = engine.retrieve(query)
    elapsed = time.time() - t0

    print(f"  Query: {query.question}")
    print(f"  Results: {evidence.hit_count} hits in {elapsed*1000:.0f}ms")
    print(f"  Collections searched: {evidence.total_collections_searched}")

    by_coll = evidence.hits_by_collection()
    for coll, hits in by_coll.items():
        top = hits[0] if hits else None
        print(f"    {coll}: {len(hits)} hits" + (f" (top: {top.score:.3f})" if top else ""))

    if evidence.knowledge_context:
        ctx_lines = evidence.knowledge_context.count("\n") + 1
        print(f"  Knowledge context: {ctx_lines} lines")
    else:
        print("  Knowledge context: none")

    assert evidence.hit_count > 0, "retrieve() returned no hits"
    print("  PASS: retrieve() works")

    # ── Step 4: Test full query() with Claude ─────────────────────────
    print("\n[4/5] Testing full query() -> Claude LLM response...")
    print(f"  Sending to {llm.model}...")

    t0 = time.time()
    answer = engine.query(DEMO_QUERIES[0])
    elapsed = time.time() - t0

    print(f"  Response received in {elapsed:.1f}s ({len(answer)} chars)")
    print("  First 500 chars:")
    print("-" * 50)
    print(answer[:500])
    print("-" * 50)

    assert len(answer) > 100, "LLM response too short"
    print("  PASS: Full RAG query works")

    # ── Step 5: Test streaming with a second query ────────────────────
    print(f"\n[5/5] Testing query_stream() with: {DEMO_QUERIES[1][:60]}...")

    t0 = time.time()
    full_answer = ""
    chunk_count = 0
    evidence_received = False

    for event in engine.query_stream(DEMO_QUERIES[1]):
        if event["type"] == "evidence":
            evidence_received = True
            ev = event["content"]
            print(f"  Evidence: {ev.hit_count} hits ({ev.search_time_ms:.0f}ms)")
        elif event["type"] == "token":
            full_answer += event["content"]
            chunk_count += 1
        elif event["type"] == "done":
            pass

    elapsed = time.time() - t0
    print(f"  Streamed {chunk_count} chunks in {elapsed:.1f}s ({len(full_answer)} chars)")
    print("  First 300 chars:")
    print("-" * 50)
    print(full_answer[:300])
    print("-" * 50)

    assert evidence_received, "No evidence event received"
    assert len(full_answer) > 100, "Streamed response too short"
    print("  PASS: Streaming RAG works")

    # ── Summary ───────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("ALL INTEGRATION TESTS PASSED")
    print(f"  Data layer: {total:,} vectors")
    print(f"  Knowledge graph: {kg.get_knowledge_stats()}")
    print(f"  Query expansion: {len(qe.ALL_EXPANSION_MAPS)} categories")
    print(f"  LLM: {llm.model}")
    print("=" * 70)

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
