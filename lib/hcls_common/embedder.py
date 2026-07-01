"""
Unified embedding provider for the HCLS AI Factory.

Consolidates patterns from:
  - core/engines/precision-intelligence/src/embedder.py  (EvidenceEmbedder / CachedEmbedder with
    sentence-transformers BGE, disk cache)
  - CAR-T agent embedding usage

Providers:
  - LocalEmbedder   -- sentence-transformers on-device (BAAI/bge-small-en-v1.5)
  - TEIEmbedder     -- NVIDIA / HuggingFace Text Embeddings Inference HTTP API

Features:
  - ABC base class for swappable backends
  - LRU memory cache (10K entries, SHA-256 key hashing)
  - Disk cache for variant embeddings (npz)
  - BGE instruction prefix for queries
  - Prometheus metrics (encode latency, cache hits)
  - Singleton ``get_embedder()``
"""

from __future__ import annotations

import hashlib
import logging
import threading
import time
from abc import ABC, abstractmethod
from collections import OrderedDict
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Prometheus metrics (optional)
# ---------------------------------------------------------------------------

try:
    from prometheus_client import Counter, Histogram

    _encode_latency = Histogram(
        "hcls_embedder_encode_seconds",
        "Time to encode a batch of texts",
        ["provider"],
        buckets=(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5),
    )
    _cache_hits = Counter(
        "hcls_embedder_cache_hits_total",
        "Embedding cache hits",
        ["level"],  # memory | disk
    )
    _cache_misses = Counter(
        "hcls_embedder_cache_misses_total",
        "Embedding cache misses",
    )
    _HAS_PROMETHEUS = True
except ImportError:
    _HAS_PROMETHEUS = False


# ---------------------------------------------------------------------------
# LRU memory cache
# ---------------------------------------------------------------------------

class _EmbeddingMemoryCache:
    """
    Thread-safe LRU memory cache for embedding vectors keyed by SHA-256
    hash of the input text.
    """

    def __init__(self, max_size: int = 10_000) -> None:
        self._max_size = max_size
        self._lock = threading.Lock()
        self._store: OrderedDict[str, np.ndarray] = OrderedDict()

    @staticmethod
    def _key(text: str) -> str:
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    def get(self, text: str) -> Optional[np.ndarray]:
        key = self._key(text)
        with self._lock:
            vec = self._store.get(key)
            if vec is not None:
                self._store.move_to_end(key)
                return vec
        return None

    def put(self, text: str, vec: np.ndarray) -> None:
        key = self._key(text)
        with self._lock:
            if key in self._store:
                self._store.move_to_end(key)
            self._store[key] = vec
            while len(self._store) > self._max_size:
                self._store.popitem(last=False)

    def get_many(self, texts: List[str]) -> Dict[int, np.ndarray]:
        """Return {index: vector} for cache hits."""
        found: Dict[int, np.ndarray] = {}
        with self._lock:
            for i, text in enumerate(texts):
                key = self._key(text)
                vec = self._store.get(key)
                if vec is not None:
                    self._store.move_to_end(key)
                    found[i] = vec
        return found

    def put_many(self, texts: List[str], vecs: np.ndarray) -> None:
        with self._lock:
            for text, vec in zip(texts, vecs):
                key = self._key(text)
                if key in self._store:
                    self._store.move_to_end(key)
                self._store[key] = vec
                while len(self._store) > self._max_size:
                    self._store.popitem(last=False)

    def clear(self) -> None:
        with self._lock:
            self._store.clear()

    def __len__(self) -> int:
        return len(self._store)


# ---------------------------------------------------------------------------
# Disk cache for variant embeddings
# ---------------------------------------------------------------------------

class _DiskCache:
    """
    Simple npz-backed disk cache for variant embeddings, compatible with
    the CachedEmbedder pattern from core/engines/precision-intelligence.
    """

    def __init__(self, cache_dir: Path, model_name: str) -> None:
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        model_slug = model_name.replace("/", "_")
        self._file = self._cache_dir / f"embeddings_{model_slug}.npz"
        self._data: Dict[str, np.ndarray] = {}
        self._dirty = False
        self._load()

    def _load(self) -> None:
        if self._file.exists():
            try:
                npz = np.load(self._file, allow_pickle=True)
                self._data = dict(npz["cache"].item())
                logger.info("Disk cache: loaded %d entries from %s", len(self._data), self._file)
            except Exception as exc:
                logger.warning("Disk cache: could not load %s: %s", self._file, exc)
                self._data = {}

    def get(self, key: str) -> Optional[np.ndarray]:
        return self._data.get(key)

    def put(self, key: str, vec: np.ndarray) -> None:
        self._data[key] = vec
        self._dirty = True

    def flush(self) -> None:
        if self._dirty:
            try:
                np.savez_compressed(self._file, cache=self._data)
                self._dirty = False
                logger.debug("Disk cache: saved %d entries", len(self._data))
            except Exception as exc:
                logger.warning("Disk cache: save failed: %s", exc)

    def clear(self) -> None:
        self._data.clear()
        self._dirty = False
        if self._file.exists():
            self._file.unlink()

    def __len__(self) -> int:
        return len(self._data)


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

class BaseEmbedder(ABC):
    """
    Abstract embedding provider.

    Subclasses must implement ``_encode_batch`` which receives a list of raw
    text strings and returns a 2-D ``np.ndarray`` of shape
    ``(len(texts), dimension)``.
    """

    def __init__(
        self,
        model_name: str = "BAAI/bge-small-en-v1.5",
        dimension: int = 384,
        batch_size: int = 32,
        cache_dir: Optional[Path] = None,
        memory_cache_size: int = 10_000,
    ) -> None:
        self.model_name = model_name
        self.dimension = dimension
        self.batch_size = batch_size

        self._mem_cache = _EmbeddingMemoryCache(max_size=memory_cache_size)
        self._disk_cache: Optional[_DiskCache] = None
        if cache_dir:
            self._disk_cache = _DiskCache(cache_dir, model_name)

    @abstractmethod
    def _encode_batch(self, texts: List[str]) -> np.ndarray:
        """Encode raw text strings (no caching)."""
        ...

    @property
    @abstractmethod
    def provider_name(self) -> str:
        """Label for metrics (e.g. ``local_bge``, ``tei``)."""
        ...

    # ── public interface ---------------------------------------------------

    def embed_text(self, text: str) -> np.ndarray:
        """Embed a single text string (uses memory cache)."""
        cached = self._mem_cache.get(text)
        if cached is not None:
            if _HAS_PROMETHEUS:
                _cache_hits.labels(level="memory").inc()
            return cached
        if _HAS_PROMETHEUS:
            _cache_misses.inc()

        vec = self._encode_batch([text])[0]
        self._mem_cache.put(text, vec)
        return vec

    def embed_texts(self, texts: List[str]) -> np.ndarray:
        """Embed multiple texts with batching and memory caching."""
        if not texts:
            return np.empty((0, self.dimension), dtype=np.float32)

        # Check memory cache
        cached = self._mem_cache.get_many(texts)
        missing_indices = [i for i in range(len(texts)) if i not in cached]

        if _HAS_PROMETHEUS:
            _cache_hits.labels(level="memory").inc(len(cached))
            _cache_misses.inc(len(missing_indices))

        if not missing_indices:
            # All cached
            return np.array([cached[i] for i in range(len(texts))])

        # Encode missing in batches
        missing_texts = [texts[i] for i in missing_indices]
        all_new: List[np.ndarray] = []
        for start in range(0, len(missing_texts), self.batch_size):
            chunk = missing_texts[start : start + self.batch_size]

            t0 = time.monotonic()
            batch_vecs = self._encode_batch(chunk)
            elapsed = time.monotonic() - t0

            if _HAS_PROMETHEUS:
                _encode_latency.labels(provider=self.provider_name).observe(elapsed)

            all_new.append(batch_vecs)

        new_vecs = np.concatenate(all_new, axis=0)

        # Update cache
        self._mem_cache.put_many(missing_texts, new_vecs)

        # Assemble full result
        result = np.empty((len(texts), self.dimension), dtype=np.float32)
        for i, vec in cached.items():
            result[i] = vec
        for j, idx in enumerate(missing_indices):
            result[idx] = new_vecs[j]

        return result

    def embed_query(self, query: str) -> np.ndarray:
        """
        Embed a search query.  For BGE models the instruction prefix is
        prepended automatically.
        """
        if "bge" in self.model_name.lower():
            query = (
                f"Represent this sentence for searching relevant passages: {query}"
            )
        return self.embed_text(query)

    def embed_with_variant_id(
        self,
        variant_id: str,
        text: str,
    ) -> np.ndarray:
        """
        Embed with disk-cache keyed by ``variant_id``.

        Useful for large-scale ingestion where the same variant may be
        re-processed across runs.
        """
        if self._disk_cache is not None:
            cached = self._disk_cache.get(variant_id)
            if cached is not None:
                if _HAS_PROMETHEUS:
                    _cache_hits.labels(level="disk").inc()
                return cached

        vec = self.embed_text(text)

        if self._disk_cache is not None:
            self._disk_cache.put(variant_id, vec)

        return vec

    def flush_disk_cache(self) -> None:
        """Persist the disk cache to disk."""
        if self._disk_cache is not None:
            self._disk_cache.flush()

    def clear_caches(self) -> None:
        """Clear both memory and disk caches."""
        self._mem_cache.clear()
        if self._disk_cache is not None:
            self._disk_cache.clear()


# ---------------------------------------------------------------------------
# LocalEmbedder (sentence-transformers)
# ---------------------------------------------------------------------------

class LocalEmbedder(BaseEmbedder):
    """
    On-device embedding via ``sentence-transformers``.

    Default model: BAAI/bge-small-en-v1.5 (384 dim, fast, good quality).
    """

    def __init__(
        self,
        model_name: str = "BAAI/bge-small-en-v1.5",
        dimension: int = 384,
        batch_size: int = 32,
        device: Optional[str] = None,
        cache_dir: Optional[Path] = None,
        model_cache_folder: Optional[str] = None,
        memory_cache_size: int = 10_000,
    ) -> None:
        super().__init__(
            model_name=model_name,
            dimension=dimension,
            batch_size=batch_size,
            cache_dir=cache_dir,
            memory_cache_size=memory_cache_size,
        )

        try:
            from sentence_transformers import SentenceTransformer
        except ImportError as exc:
            raise ImportError(
                "sentence-transformers is required for LocalEmbedder. "
                "Install with: pip install sentence-transformers"
            ) from exc

        logger.info("Loading local embedding model: %s", model_name)
        self._model = SentenceTransformer(
            model_name,
            cache_folder=model_cache_folder,
            device=device,
        )
        actual_dim = self._model.get_sentence_embedding_dimension()
        if actual_dim != dimension:
            logger.warning(
                "Model dimension (%d) differs from configured dimension (%d). "
                "Using actual model dimension.",
                actual_dim,
                dimension,
            )
            self.dimension = actual_dim
        logger.info(
            "Local embedder ready: %s  dim=%d  device=%s",
            model_name,
            self.dimension,
            self._model.device,
        )

    @property
    def provider_name(self) -> str:
        return "local_bge"

    def _encode_batch(self, texts: List[str]) -> np.ndarray:
        return self._model.encode(
            texts,
            normalize_embeddings=True,
            batch_size=self.batch_size,
            show_progress_bar=False,
        )


# ---------------------------------------------------------------------------
# TEIEmbedder (Text Embeddings Inference HTTP API)
# ---------------------------------------------------------------------------

class TEIEmbedder(BaseEmbedder):
    """
    Embedding via a remote Text Embeddings Inference (TEI) server.

    Compatible with NVIDIA NIM TEI and HuggingFace TEI containers.
    """

    def __init__(
        self,
        endpoint: str = "http://localhost:8081",
        model_name: str = "BAAI/bge-small-en-v1.5",
        dimension: int = 384,
        batch_size: int = 32,
        timeout: float = 30.0,
        cache_dir: Optional[Path] = None,
        memory_cache_size: int = 10_000,
    ) -> None:
        super().__init__(
            model_name=model_name,
            dimension=dimension,
            batch_size=batch_size,
            cache_dir=cache_dir,
            memory_cache_size=memory_cache_size,
        )
        self._endpoint = endpoint.rstrip("/")
        self._timeout = timeout

        try:
            import requests  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "requests is required for TEIEmbedder. "
                "Install with: pip install requests"
            ) from exc

        logger.info(
            "TEI embedder configured: endpoint=%s  model=%s  dim=%d",
            self._endpoint,
            model_name,
            dimension,
        )

    @property
    def provider_name(self) -> str:
        return "tei"

    def _encode_batch(self, texts: List[str]) -> np.ndarray:
        import requests

        payload = {"inputs": texts, "normalize": True}
        resp = requests.post(
            f"{self._endpoint}/embed",
            json=payload,
            timeout=self._timeout,
        )
        resp.raise_for_status()
        embeddings = resp.json()

        # TEI returns list of lists
        return np.array(embeddings, dtype=np.float32)

    def health_check(self) -> bool:
        """Check if TEI endpoint is reachable."""
        import requests

        try:
            resp = requests.get(f"{self._endpoint}/health", timeout=5)
            return resp.status_code == 200
        except requests.RequestException:
            return False


# ---------------------------------------------------------------------------
# Singleton accessor
# ---------------------------------------------------------------------------

_embedder_instance: Optional[BaseEmbedder] = None
_embedder_lock = threading.Lock()


def get_embedder(settings: Any = None) -> BaseEmbedder:
    """
    Return a singleton embedder instance based on ``HCLSSettings``.

    Uses ``embedding_provider`` to choose between ``LocalEmbedder`` and
    ``TEIEmbedder``.  Pass ``settings=None`` to auto-load via
    ``get_settings()``.
    """
    global _embedder_instance
    if _embedder_instance is not None:
        return _embedder_instance

    with _embedder_lock:
        # Double-check after acquiring lock
        if _embedder_instance is not None:
            return _embedder_instance

        if settings is None:
            from hcls_common.config import get_settings

            settings = get_settings()

        provider = getattr(settings, "embedding_provider", "local_bge").lower()
        model = getattr(settings, "embedding_model", "BAAI/bge-small-en-v1.5")
        dim = getattr(settings, "embedding_dimension", 384)
        batch = getattr(settings, "embedding_batch_size", 32)
        cache_dir = getattr(settings, "embedding_cache_dir", None)
        device = getattr(settings, "embedding_device", None)

        if provider == "tei":
            endpoint = getattr(settings, "tei_endpoint", None)
            if not endpoint:
                raise ValueError(
                    "tei_endpoint must be set when embedding_provider=tei"
                )
            _embedder_instance = TEIEmbedder(
                endpoint=endpoint,
                model_name=model,
                dimension=dim,
                batch_size=batch,
                cache_dir=cache_dir,
            )
        else:
            # Default: local_bge
            _embedder_instance = LocalEmbedder(
                model_name=model,
                dimension=dim,
                batch_size=batch,
                device=device,
                cache_dir=cache_dir,
            )

        return _embedder_instance


def reset_embedder() -> None:
    """Clear the singleton (useful in tests)."""
    global _embedder_instance
    with _embedder_lock:
        _embedder_instance = None
