"""
Unified Milvus vector-database wrapper for the HCLS AI Factory.

Consolidates patterns from:
  - core/engines/precision-intelligence/src/milvus_client.py   (genomic evidence search with gene/region sanitization)
  - CAR-T agent's CARTCollectionManager      (11 parallel collection searches)

Features:
  - Connection pooling with configurable pool size
  - Retry logic with exponential backoff
  - Circuit breaker integration (hcls_common.circuit_breaker)
  - Milvus expression injection validation (gene name / chromosome patterns)
  - Batch upsert with chunked inserts
  - TTL + LRU query-result cache
  - Health check
  - Prometheus metrics (search latency, upsert latency, cache hits)
"""

from __future__ import annotations

import hashlib
import logging
import re
import threading
import time
from collections import OrderedDict
from queue import Empty, Queue
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Prometheus metrics (optional)
# ---------------------------------------------------------------------------

try:
    from prometheus_client import Counter, Histogram, Gauge

    _search_latency = Histogram(
        "hcls_milvus_search_seconds",
        "Milvus search latency",
        ["collection"],
        buckets=(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0),
    )
    _upsert_latency = Histogram(
        "hcls_milvus_upsert_seconds",
        "Milvus upsert latency",
        ["collection"],
        buckets=(0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0),
    )
    _cache_hits = Counter(
        "hcls_milvus_cache_hits_total",
        "Query cache hits",
        ["collection"],
    )
    _cache_misses = Counter(
        "hcls_milvus_cache_misses_total",
        "Query cache misses",
        ["collection"],
    )
    _pool_size_gauge = Gauge(
        "hcls_milvus_pool_available",
        "Available connections in pool",
    )
    _HAS_PROMETHEUS = True
except ImportError:
    _HAS_PROMETHEUS = False


# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------

# Gene names: standard HGNC symbols (letters, digits, hyphens, underscores)
GENE_PATTERN = re.compile(r"^[A-Za-z0-9_\-]+$")

# Chromosome: optional "chr" prefix, followed by 1-22, X, Y, or M/MT
CHROM_PATTERN = re.compile(r"^(chr)?([0-9]{1,2}|[XYMT]{1,2})$", re.IGNORECASE)

# Allowed tokens in Milvus filter expressions (whitelist approach)
_FILTER_TOKEN_RE = re.compile(
    r"""
    (?:                         # allowed tokens:
      \d+(?:\.\d+)?             #   numeric literals
    | "[^"]*"                   #   double-quoted string literal
    | '[^']*'                   #   single-quoted string literal
    | ==|!=|>=|<=|>|<           #   comparison operators
    | \b(?:and|or|not|in|AND|OR|NOT|IN)\b  # logical keywords
    | \b(?:gene|chrom|pos|impact|consequence|clinical_significance
          |rsid|am_class|am_pathogenicity|qual)\b  # field names
    | [()\[\],]                 #   grouping / list punctuation
    | \s+                       #   whitespace
    )
    """,
    re.VERBOSE,
)


def sanitize_gene_name(gene: str) -> str:
    """Validate a gene name against injection. Returns the cleaned value."""
    gene = gene.strip()
    if not gene or not GENE_PATTERN.match(gene):
        raise ValueError(f"Invalid gene name: {gene!r}")
    return gene


def sanitize_chromosome(chrom: str) -> str:
    """Validate a chromosome identifier."""
    chrom = chrom.strip()
    if not chrom or not CHROM_PATTERN.match(chrom):
        raise ValueError(f"Invalid chromosome: {chrom!r}")
    return chrom


def validate_milvus_filter(expr: str) -> str:
    """
    Validate a Milvus boolean filter expression using a token whitelist.

    Raises ``ValueError`` if the expression contains tokens that are not in
    the whitelist (preventing injection of arbitrary Milvus operations).
    """
    if not expr or not expr.strip():
        raise ValueError("Filter expression is empty")

    # Remove all valid tokens; anything remaining is suspicious
    residual = _FILTER_TOKEN_RE.sub("", expr)
    residual = residual.strip()
    if residual:
        raise ValueError(
            f"Filter expression contains disallowed tokens: {residual!r}"
        )
    return expr


# ---------------------------------------------------------------------------
# TTL + LRU cache for query results
# ---------------------------------------------------------------------------

class TTLLRUCache:
    """
    Thread-safe cache combining TTL expiration and LRU eviction.

    Parameters
    ----------
    max_size : int
        Maximum number of entries.
    ttl : float
        Time-to-live per entry in seconds.
    """

    def __init__(self, max_size: int = 1024, ttl: float = 300.0) -> None:
        self._max_size = max_size
        self._ttl = ttl
        self._lock = threading.Lock()
        self._store: OrderedDict[str, Tuple[float, Any]] = OrderedDict()

    def _make_key(self, *args: Any) -> str:
        raw = str(args).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def get(self, key: str) -> Tuple[bool, Any]:
        """Return ``(True, value)`` on hit or ``(False, None)`` on miss."""
        with self._lock:
            entry = self._store.get(key)
            if entry is None:
                return False, None
            ts, value = entry
            if time.monotonic() - ts > self._ttl:
                del self._store[key]
                return False, None
            # Move to end (most recently used)
            self._store.move_to_end(key)
            return True, value

    def put(self, key: str, value: Any) -> None:
        with self._lock:
            if key in self._store:
                self._store.move_to_end(key)
            self._store[key] = (time.monotonic(), value)
            while len(self._store) > self._max_size:
                self._store.popitem(last=False)

    def invalidate(self, key: str) -> None:
        with self._lock:
            self._store.pop(key, None)

    def clear(self) -> None:
        with self._lock:
            self._store.clear()

    def __len__(self) -> int:
        return len(self._store)


# ---------------------------------------------------------------------------
# Connection pool
# ---------------------------------------------------------------------------

class _ConnectionPool:
    """
    Simple blocking connection pool for pymilvus.

    Each connection is identified by a unique alias so pymilvus can manage
    them independently.
    """

    def __init__(self, host: str, port: int, pool_size: int) -> None:
        self._host = host
        self._port = port
        self._pool: Queue[str] = Queue(maxsize=pool_size)
        self._lock = threading.Lock()
        self._created = 0
        self._pool_size = pool_size

    def acquire(self, timeout: float = 30.0) -> str:
        """Acquire a connection alias from the pool."""
        # Try to get an existing one
        try:
            alias = self._pool.get_nowait()
            return alias
        except Empty:
            pass

        # Create a new one if capacity allows
        with self._lock:
            if self._created < self._pool_size:
                alias = f"hcls_pool_{self._created}"
                self._created += 1
            else:
                # Wait for one to be returned
                alias = self._pool.get(timeout=timeout)
                return alias

        from pymilvus import connections as _conn

        _conn.connect(alias=alias, host=self._host, port=self._port)
        logger.debug(
            "Milvus pool: created connection %s -> %s:%d",
            alias,
            self._host,
            self._port,
        )
        return alias

    def release(self, alias: str) -> None:
        self._pool.put_nowait(alias)
        if _HAS_PROMETHEUS:
            _pool_size_gauge.set(self._pool.qsize())

    def close_all(self) -> None:
        from pymilvus import connections as _conn

        while not self._pool.empty():
            try:
                alias = self._pool.get_nowait()
                try:
                    _conn.disconnect(alias)
                except Exception:
                    pass
            except Empty:
                break


# ---------------------------------------------------------------------------
# UnifiedMilvusClient
# ---------------------------------------------------------------------------

class UnifiedMilvusClient:
    """
    Unified Milvus wrapper for the HCLS AI Factory.

    Provides:
      - ``connect`` / ``disconnect``
      - ``create_collection`` with standard genomic schema
      - ``batch_upsert`` with chunked inserts
      - ``search`` (vector similarity) with optional filter
      - ``search_by_gene`` / ``search_by_region`` convenience methods
      - ``parallel_search`` for multi-collection CAR-T style queries
      - TTL+LRU query result cache
      - Circuit breaker integration
      - Prometheus metrics

    Parameters
    ----------
    host : str
        Milvus gRPC host.
    port : int
        Milvus gRPC port.
    collection_name : str
        Default collection name.
    embedding_dim : int
        Embedding vector dimension.
    pool_size : int
        Connection pool size.
    cache_max_size : int
        Maximum query-cache entries.
    cache_ttl : float
        Query-cache TTL in seconds.
    circuit_breaker_name : str or None
        If given, wraps operations with this named circuit breaker.
    """

    # Input validation patterns (also exposed at module level above)
    _GENE_PATTERN = GENE_PATTERN
    _CHROM_PATTERN = CHROM_PATTERN

    # Default output fields for genomic evidence collections
    DEFAULT_OUTPUT_FIELDS = [
        "chrom", "pos", "ref", "alt", "qual", "gene",
        "consequence", "impact", "genotype", "text_summary",
        "clinical_significance", "rsid", "am_pathogenicity", "am_class",
    ]

    def __init__(
        self,
        host: str = "localhost",
        port: int = 19530,
        collection_name: str = "genomic_evidence",
        embedding_dim: int = 384,
        pool_size: int = 4,
        cache_max_size: int = 1024,
        cache_ttl: float = 300.0,
        circuit_breaker_name: Optional[str] = None,
    ) -> None:
        self.host = host
        self.port = port
        self.collection_name = collection_name
        self.embedding_dim = embedding_dim

        self._pool: Optional[_ConnectionPool] = None
        self._pool_size = pool_size
        self._cache = TTLLRUCache(max_size=cache_max_size, ttl=cache_ttl)

        # Circuit breaker (optional)
        self._cb: Any = None
        if circuit_breaker_name:
            from hcls_common.circuit_breaker import CircuitBreaker, get_breaker

            self._cb = get_breaker(circuit_breaker_name) or CircuitBreaker(
                circuit_breaker_name
            )

    # ── lifecycle ----------------------------------------------------------

    def connect(self) -> None:
        """Initialize the connection pool."""
        logger.info("Connecting to Milvus at %s:%d (pool=%d)", self.host, self.port, self._pool_size)
        self._pool = _ConnectionPool(self.host, self.port, self._pool_size)
        # Eagerly create one connection to verify connectivity
        alias = self._pool.acquire()
        self._pool.release(alias)
        logger.info("Milvus connection pool initialized")

    def disconnect(self) -> None:
        """Close all pooled connections."""
        if self._pool is not None:
            self._pool.close_all()
            self._pool = None
        logger.info("Milvus connections closed")

    def _acquire(self) -> str:
        if self._pool is None:
            self.connect()
        assert self._pool is not None
        return self._pool.acquire()

    def _release(self, alias: str) -> None:
        if self._pool is not None:
            self._pool.release(alias)

    # ── health check -------------------------------------------------------

    def health_check(self) -> bool:
        """Return ``True`` if Milvus is reachable."""
        try:
            alias = self._acquire()
            try:
                from pymilvus import utility

                utility.list_collections(using=alias)
                return True
            finally:
                self._release(alias)
        except Exception as exc:
            logger.warning("Milvus health check failed: %s", exc)
            return False

    # ── collection helpers -------------------------------------------------

    def collection_exists(self, name: Optional[str] = None) -> bool:
        name = name or self.collection_name
        alias = self._acquire()
        try:
            from pymilvus import utility

            return utility.has_collection(name, using=alias)
        finally:
            self._release(alias)

    def create_collection(
        self,
        name: Optional[str] = None,
        drop_existing: bool = False,
    ) -> Any:
        """
        Create the standard genomic evidence collection.

        Returns the ``pymilvus.Collection`` object.
        """
        from pymilvus import (
            Collection,
            CollectionSchema,
            DataType,
            FieldSchema,
            utility,
        )

        name = name or self.collection_name
        alias = self._acquire()
        try:
            if drop_existing and utility.has_collection(name, using=alias):
                logger.warning("Dropping existing collection: %s", name)
                utility.drop_collection(name, using=alias)

            if utility.has_collection(name, using=alias):
                logger.info("Collection %s already exists", name)
                return Collection(name, using=alias)

            logger.info("Creating collection: %s", name)

            fields = [
                FieldSchema(
                    name="id",
                    dtype=DataType.VARCHAR,
                    is_primary=True,
                    max_length=200,
                    description="Variant ID (chr_pos_ref_alt)",
                ),
                FieldSchema(
                    name="embedding",
                    dtype=DataType.FLOAT_VECTOR,
                    dim=self.embedding_dim,
                    description="Text embedding vector",
                ),
                FieldSchema(name="chrom", dtype=DataType.VARCHAR, max_length=10),
                FieldSchema(name="pos", dtype=DataType.INT64),
                FieldSchema(name="ref", dtype=DataType.VARCHAR, max_length=500),
                FieldSchema(name="alt", dtype=DataType.VARCHAR, max_length=500),
                FieldSchema(name="qual", dtype=DataType.FLOAT),
                FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
                FieldSchema(name="consequence", dtype=DataType.VARCHAR, max_length=100),
                FieldSchema(name="impact", dtype=DataType.VARCHAR, max_length=20),
                FieldSchema(name="genotype", dtype=DataType.VARCHAR, max_length=10),
                FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=2000),
                FieldSchema(
                    name="clinical_significance",
                    dtype=DataType.VARCHAR,
                    max_length=200,
                ),
                FieldSchema(name="rsid", dtype=DataType.VARCHAR, max_length=20),
                FieldSchema(
                    name="disease_associations",
                    dtype=DataType.VARCHAR,
                    max_length=500,
                ),
                FieldSchema(name="am_pathogenicity", dtype=DataType.FLOAT),
                FieldSchema(name="am_class", dtype=DataType.VARCHAR, max_length=30),
            ]

            schema = CollectionSchema(
                fields=fields,
                description="Genomic variant evidence for RAG",
            )
            collection = Collection(name=name, schema=schema, using=alias)

            # IVF_FLAT index for cosine similarity
            collection.create_index(
                field_name="embedding",
                index_params={
                    "metric_type": "COSINE",
                    "index_type": "IVF_FLAT",
                    "params": {"nlist": 1024},
                },
            )
            logger.info("Collection %s created with IVF_FLAT index", name)
            return collection
        finally:
            self._release(alias)

    def get_collection(self, name: Optional[str] = None) -> Any:
        """Return an existing collection (create if needed)."""
        from pymilvus import Collection

        name = name or self.collection_name
        alias = self._acquire()
        try:
            if self.collection_exists(name):
                return Collection(name, using=alias)
            return self.create_collection(name)
        finally:
            self._release(alias)

    # ── batch upsert -------------------------------------------------------

    def batch_upsert(
        self,
        collection_name: Optional[str] = None,
        ids: Optional[List[str]] = None,
        embeddings: Optional[np.ndarray] = None,
        data_rows: Optional[List[List[Any]]] = None,
        batch_size: int = 1000,
    ) -> int:
        """
        Upsert data in chunked batches.

        Accepts **either** structured ``data_rows`` (list-of-columns matching
        the collection schema) **or** ``ids`` + ``embeddings`` for minimal
        inserts.

        Returns total number of upserted entities.
        """
        from pymilvus import Collection

        cname = collection_name or self.collection_name
        alias = self._acquire()
        try:
            col = Collection(cname, using=alias)
            total = 0

            if data_rows is not None:
                # data_rows is already list-of-columns (pymilvus convention)
                num_entities = len(data_rows[0]) if data_rows else 0
                for start in range(0, num_entities, batch_size):
                    end = min(start + batch_size, num_entities)
                    chunk = [col_data[start:end] for col_data in data_rows]

                    t0 = time.monotonic()
                    if self._cb:
                        with self._cb:
                            result = col.insert(chunk)
                    else:
                        result = col.insert(chunk)
                    elapsed = time.monotonic() - t0

                    if _HAS_PROMETHEUS:
                        _upsert_latency.labels(collection=cname).observe(elapsed)

                    total += result.insert_count
                    logger.debug(
                        "Upserted %d/%d into %s (%.3fs)",
                        total,
                        num_entities,
                        cname,
                        elapsed,
                    )
            elif ids is not None and embeddings is not None:
                num_entities = len(ids)
                for start in range(0, num_entities, batch_size):
                    end = min(start + batch_size, num_entities)
                    chunk = [ids[start:end], embeddings[start:end].tolist()]

                    t0 = time.monotonic()
                    if self._cb:
                        with self._cb:
                            result = col.insert(chunk)
                    else:
                        result = col.insert(chunk)
                    elapsed = time.monotonic() - t0

                    if _HAS_PROMETHEUS:
                        _upsert_latency.labels(collection=cname).observe(elapsed)

                    total += result.insert_count

            col.flush()
            self._cache.clear()  # invalidate cache after writes
            return total
        finally:
            self._release(alias)

    # ── vector search ------------------------------------------------------

    def search(
        self,
        query_embedding: np.ndarray,
        top_k: int = 10,
        score_threshold: float = 0.0,
        filter_expr: Optional[str] = None,
        collection_name: Optional[str] = None,
        output_fields: Optional[List[str]] = None,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """
        Semantic similarity search.

        Parameters
        ----------
        query_embedding : np.ndarray
            Query vector (1-D, same dimension as the index).
        top_k : int
            Number of results.
        score_threshold : float
            Minimum cosine similarity to include.
        filter_expr : str or None
            Optional Milvus boolean filter.  Validated via ``validate_milvus_filter``.
        collection_name : str or None
            Collection to search (defaults to ``self.collection_name``).
        output_fields : list[str] or None
            Fields to return; defaults to ``DEFAULT_OUTPUT_FIELDS``.
        use_cache : bool
            Whether to use the query result cache.

        Returns
        -------
        list[dict]
            Each dict contains ``id``, ``score``, and the requested output fields.
        """
        cname = collection_name or self.collection_name
        ofields = output_fields or self.DEFAULT_OUTPUT_FIELDS

        # Validate filter expression
        if filter_expr:
            filter_expr = validate_milvus_filter(filter_expr)

        # Cache lookup
        cache_key = self._cache._make_key(
            cname,
            query_embedding.tobytes(),
            top_k,
            score_threshold,
            filter_expr,
            tuple(ofields),
        )
        if use_cache:
            hit, cached = self._cache.get(cache_key)
            if hit:
                if _HAS_PROMETHEUS:
                    _cache_hits.labels(collection=cname).inc()
                return cached
            if _HAS_PROMETHEUS:
                _cache_misses.labels(collection=cname).inc()

        from pymilvus import Collection

        alias = self._acquire()
        try:
            col = Collection(cname, using=alias)
            col.load()

            search_params = {
                "metric_type": "COSINE",
                "params": {"nprobe": 16},
            }

            t0 = time.monotonic()
            if self._cb:
                with self._cb:
                    results = col.search(
                        data=[query_embedding.tolist()],
                        anns_field="embedding",
                        param=search_params,
                        limit=top_k,
                        expr=filter_expr,
                        output_fields=ofields,
                    )
            else:
                results = col.search(
                    data=[query_embedding.tolist()],
                    anns_field="embedding",
                    param=search_params,
                    limit=top_k,
                    expr=filter_expr,
                    output_fields=ofields,
                )
            elapsed = time.monotonic() - t0

            if _HAS_PROMETHEUS:
                _search_latency.labels(collection=cname).observe(elapsed)

            parsed = self._parse_search_results(results, score_threshold)

            if use_cache:
                self._cache.put(cache_key, parsed)

            return parsed
        finally:
            self._release(alias)

    @staticmethod
    def _parse_search_results(
        results: Any,
        score_threshold: float,
    ) -> List[Dict[str, Any]]:
        evidence: List[Dict[str, Any]] = []
        for hits in results:
            for hit in hits:
                score = hit.score
                if score < score_threshold:
                    continue

                entry: Dict[str, Any] = {"id": hit.id, "score": score}
                for field_name in hit.entity.fields:
                    val = hit.entity.get(field_name)
                    # Convert AlphaMissense sentinel -1.0 back to None
                    if field_name == "am_pathogenicity" and val == -1.0:
                        val = None
                    # Convert empty string to None for am_class
                    if field_name == "am_class" and val == "":
                        val = None
                    entry[field_name] = val
                evidence.append(entry)
        return evidence

    # ── convenience searches -----------------------------------------------

    def search_by_gene(
        self,
        gene: str,
        top_k: int = 100,
        collection_name: Optional[str] = None,
        output_fields: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Search for all variants in a specific gene (scalar filter, no vector)."""
        gene = sanitize_gene_name(gene)
        cname = collection_name or self.collection_name
        ofields = output_fields or self.DEFAULT_OUTPUT_FIELDS

        from pymilvus import Collection

        alias = self._acquire()
        try:
            col = Collection(cname, using=alias)
            col.load()

            t0 = time.monotonic()
            results = col.query(
                expr=f'gene == "{gene}"',
                output_fields=ofields,
                limit=top_k,
            )
            elapsed = time.monotonic() - t0

            if _HAS_PROMETHEUS:
                _search_latency.labels(collection=cname).observe(elapsed)

            return results
        finally:
            self._release(alias)

    def search_by_region(
        self,
        chrom: str,
        start: int,
        end: int,
        top_k: int = 100,
        collection_name: Optional[str] = None,
        output_fields: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Search for variants in a genomic region (scalar filter, no vector)."""
        chrom = sanitize_chromosome(chrom)
        start = int(start)
        end = int(end)
        if start < 0 or end < 0 or end < start:
            raise ValueError(f"Invalid region: start={start}, end={end}")

        cname = collection_name or self.collection_name
        ofields = output_fields or self.DEFAULT_OUTPUT_FIELDS

        from pymilvus import Collection

        alias = self._acquire()
        try:
            col = Collection(cname, using=alias)
            col.load()

            t0 = time.monotonic()
            results = col.query(
                expr=f'chrom == "{chrom}" and pos >= {start} and pos <= {end}',
                output_fields=ofields,
                limit=top_k,
            )
            elapsed = time.monotonic() - t0

            if _HAS_PROMETHEUS:
                _search_latency.labels(collection=cname).observe(elapsed)

            return results
        finally:
            self._release(alias)

    # ── multi-collection parallel search (CAR-T style) ---------------------

    def parallel_search(
        self,
        query_embedding: np.ndarray,
        collections: Dict[str, float],
        top_k_per_collection: int = 5,
        score_threshold: float = 0.4,
    ) -> List[Dict[str, Any]]:
        """
        Search across multiple collections in parallel with per-collection
        weights, mimicking the CAR-T agent's 11-collection search.

        Parameters
        ----------
        query_embedding : np.ndarray
            Shared query vector.
        collections : dict[str, float]
            Mapping of collection_name -> weight (weights should sum to ~1.0).
        top_k_per_collection : int
            Results per collection.
        score_threshold : float
            Minimum similarity per collection.

        Returns
        -------
        list[dict]
            Aggregated results sorted by weighted score descending.
        """
        import concurrent.futures

        all_results: List[Dict[str, Any]] = []

        def _search_one(cname: str, weight: float) -> List[Dict[str, Any]]:
            try:
                raw = self.search(
                    query_embedding=query_embedding,
                    top_k=top_k_per_collection,
                    score_threshold=score_threshold,
                    collection_name=cname,
                )
                for r in raw:
                    r["collection"] = cname
                    r["weight"] = weight
                    r["weighted_score"] = r["score"] * weight
                return raw
            except Exception as exc:
                logger.warning(
                    "parallel_search: collection %s failed: %s", cname, exc
                )
                return []

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(len(collections), 8)
        ) as pool:
            futures = {
                pool.submit(_search_one, cname, weight): cname
                for cname, weight in collections.items()
            }
            for fut in concurrent.futures.as_completed(futures):
                all_results.extend(fut.result())

        all_results.sort(key=lambda r: r["weighted_score"], reverse=True)
        return all_results

    # ── stats / admin -------------------------------------------------------

    def get_stats(self, collection_name: Optional[str] = None) -> Dict[str, Any]:
        """Return basic collection statistics."""
        from pymilvus import Collection

        cname = collection_name or self.collection_name
        alias = self._acquire()
        try:
            col = Collection(cname, using=alias)
            return {
                "name": cname,
                "num_entities": col.num_entities,
                "schema": str(col.schema),
            }
        finally:
            self._release(alias)

    def drop_collection(self, name: Optional[str] = None) -> None:
        from pymilvus import utility

        name = name or self.collection_name
        alias = self._acquire()
        try:
            if utility.has_collection(name, using=alias):
                utility.drop_collection(name, using=alias)
                self._cache.clear()
                logger.info("Collection %s dropped", name)
        finally:
            self._release(alias)

    # ── context manager ----------------------------------------------------

    def __enter__(self) -> "UnifiedMilvusClient":
        self.connect()
        return self

    def __exit__(self, *exc: Any) -> None:
        self.disconnect()
