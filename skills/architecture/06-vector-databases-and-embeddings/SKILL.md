---
name: 06-vector-databases-and-embeddings
description: >-
  Best-practice standards for Pillar 06 (Vector Databases & Embeddings) of the HCLS AI Factory. Use when
  designing, building, operating, or reviewing semantic retrieval, similarity search, or embedding
  pipelines. Triggers: adding a RAG collection, choosing an index/metric, pinning an embedding model,
  building a protein/text similarity index, or wiring retrieve-then-rerank.
---

# Pillar 06 — Vector Databases & Embeddings

The semantic-retrieval tier: one shared vector database backs every engine and agent's RAG and every
similarity index, fed by pinned, domain-appropriate embedding models. On one box this is deliberately
a single, reused store — not one vector database per service.

## In the HCLS AI Factory
- **Milvus 2.4 — the one shared vector database** (`vector-db`, :19530, `serving: container`, tags
  `milvus`, `infra`). ~3.5M vectors, IVF_FLAT + COSINE. Backs all RAG retrieval (the Precision
  Intelligence Engine and the eight agents) and every similarity index. Accessed through the shared
  `UnifiedMilvusClient` (`milvus_client.py`) — connection pool, TTL/LRU cache, Prometheus latency
  metrics, sanitized filter expressions, and `parallel_search` for multi-collection (CAR-T-style)
  fan-out.
- **Text embeddings — BGE-small-en-v1.5** (`embedder.py`, `LocalEmbedder`, 384-dim, on-device
  sentence-transformers) with the BGE query instruction prefix, an in-memory + on-disk cache, and
  `embed_with_variant_id` for genomic evidence. Snowflake arctic-embed is the alternate backend some
  services use; the client is a swappable ABC so the backend is pluggable behind one interface.
- **Protein embeddings — ESM-2 650M** (`esm2-search`, :8571, 1280-dim, GPU on the GB10) → Milvus ANN,
  reusing the same shared vector database. Verified: ubiquitin self-retrieves cosine=1.0, a variant
  0.997. **Retrieve-then-rerank:** an optional exact Smith-Waterman re-rank (Biopython local
  alignment) over the ANN top-k adds %identity + alignment score for precision. `POST /embed /index
  /search`.

## Best-practice standards
- **Reuse the ONE shared vector database.** New RAG or similarity features are new *collections* in
  Milvus, not a new vector store. Go through `UnifiedMilvusClient`, not a raw pymilvus connection.
- **Pin the embedding model per collection.** The model + dimension are part of the collection's
  contract; write them into config and the collection metadata. Never mix embedders (or versions)
  within one collection — a re-embed is a versioned rebuild, not an in-place swap.
- **Match dimension exactly.** Text=384 (BGE), protein=1280 (ESM-2). A dimension mismatch on upsert is
  a silent-garbage bug — assert it at the client boundary.
- **Choose the index by scale.** IVF_FLAT + COSINE is the factory default and correct at ~3.5M
  vectors on one box; tune `nlist`/`nprobe` for the recall/latency trade-off. Reserve HNSW/heavier
  indexes for cases that measurably need them (memory is shared with GPU inference).
- **Retrieve-then-rerank for precision.** Use ANN to get a cheap top-k, then an exact re-rank
  (Smith-Waterman for proteins, a cross-encoder/exact metric for text) when precision matters —
  don't push k arbitrarily high on the ANN pass.
- **Normalize for cosine.** COSINE assumes normalized vectors; normalize on write and query, and keep
  the query path (BGE instruction prefix) identical between index and search.
- **Sanitize filter expressions.** Gene/chrom/scalar filters flow through the client's validation
  (`sanitize_gene_name`, `validate_milvus_filter`) — never string-concatenate a Milvus `expr` from
  user input.
- **Embeddings are recomputable; index the source of truth elsewhere.** Keep the canonical documents
  and their ids so a collection can be rebuilt; don't treat the vector database as the primary store.

## Do / Don't
**Do:** create a new Milvus collection per corpus; pin `{model, dimension, metric, index}` in config;
go through `UnifiedMilvusClient`; normalize vectors; ANN-then-exact-rerank; cache embeddings; reuse
:8571 for protein search.
**Don't:** spin up a second vector database or a per-agent store; mix embedding models/dims in one
collection; hand-roll a pymilvus client that bypasses caching/metrics/sanitization; interpolate raw
input into a filter `expr`; store raw weights/corpora in git.

## Wiring it in
- Text RAG through the shared client:
  ```python
  from hcls_common.milvus_client import UnifiedMilvusClient
  from hcls_common.embedder import LocalEmbedder          # BAAI/bge-small-en-v1.5, 384-dim
  emb = LocalEmbedder()
  cli = UnifiedMilvusClient(host="localhost", port=19530, collection_name="genomic_evidence",
                            embedding_dim=384)
  hits = cli.search(emb.embed_query("VCP pathogenic variants"), top_k=10)
  ```
- Protein similarity + exact re-rank via the service:
  ```bash
  curl -s -X POST localhost:8571/search -d '{"sequence":"MSGA...","top_k":25,"rerank":true}'
  ```
- Register a new retrieval capability in `lib/hcls_common/capabilities.json` (reference `esm2-search`),
  map the directory in `scripts/validate_registry.py`, then `python scripts/validate_registry.py`.
- Milvus runs as a container in `docker-compose.dgx-spark.yml` at 19530 with a healthcheck; connect
  via `config.py`/env, never a hardcoded host/IP.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **Milvus memory + the 128 GB unified pool.** Loaded collections and GPU inference share LPDDR5x —
  ~3.5M vectors is fine, but loading many large collections at once can pressure ESMFold/DiffDock.
  Release idle collections; prefer IVF over memory-hungry indexes.
- **Milvus has had crash-loop / recovery issues on this box** — run it on a backed, named volume, add
  a healthcheck, and never point a `live` retrieval route at a store that's mid-recovery.
- **ESM-2 embedding is GPU-served** on the GB10; if the GPU is saturated, batch and queue rather than
  falling back silently — and never label a mock/random-vector fallback as real (the honesty rule).
- **ARM wheels for sentence-transformers/torch** must be the GB10 (cu130) builds; a wrong wheel yields
  CPU-only or import failure — pin and verify.
- **Dimension/model drift is silent.** Re-embedding with a new model into an old collection returns
  plausible-but-wrong neighbors; rebuild the collection and bump its version instead.

## Related
- Pillars: 05-structured-databases, 08-inference-serving, 07-message-bus-and-async
- build-housekeeping-standards
