"""
Embeddings (PRD §2.5.4). Production uses BAAI/bge-large-en-v1.5 + a BiomedBERT-derived
clinical model via sentence-transformers. When those weights aren't available (offline /
no download), a deterministic hashing embedder keeps the RAG path runnable and tests
reproducible. Same interface either way: embed(list[str]) -> np.ndarray (L2-normalized).
"""
from __future__ import annotations

import hashlib
import os

import numpy as np

_DIM = 384
_model = None
_mode = "deterministic"


def _try_real_model():
    global _model, _mode
    if _model is not None:
        return _model
    # Real embeddings are the DEFAULT for production retrieval. The deterministic hash
    # embedder is reserved for tests / offline runs, selected EXPLICITLY via TSC_OFFLINE=1
    # (set by conftest) or TSC_USE_REAL_EMBEDDINGS=0. Any other context attempts the real
    # model and, if weights/network are unavailable, falls back to the hash embedder with a
    # loud warning — a non-semantic path must never masquerade silently as real (honesty register).
    _offline = os.environ.get("TSC_OFFLINE") == "1"
    _explicit_off = os.environ.get("TSC_USE_REAL_EMBEDDINGS") == "0"
    if _offline or _explicit_off:
        _mode = "deterministic"
        return None
    try:  # pragma: no cover - depends on local weights / network
        from sentence_transformers import SentenceTransformer

        _model = SentenceTransformer("BAAI/bge-large-en-v1.5")
        _mode = "bge-large-en-v1.5"
    except Exception as e:  # graceful, LOUD fallback — never a silent non-semantic path
        import sys
        print(
            "[TSC][WARN] real embeddings (BAAI/bge-large-en-v1.5) unavailable "
            f"({type(e).__name__}); falling back to NON-SEMANTIC deterministic hash embedder. "
            "Retrieval quality will be degraded. Install sentence-transformers + weights, "
            "or set TSC_USE_REAL_EMBEDDINGS=0 to silence.",
            file=sys.stderr,
        )
        _model = None
        _mode = "deterministic"
    return _model


def _hash_embed(text: str) -> np.ndarray:
    vec = np.zeros(_DIM, dtype=np.float32)
    for tok in text.lower().split():
        h = int(hashlib.sha1(tok.encode()).hexdigest(), 16)
        vec[h % _DIM] += 1.0
    n = np.linalg.norm(vec)
    return vec / n if n else vec


def embed(texts: list[str]) -> np.ndarray:
    model = _try_real_model()
    if model is not None:  # pragma: no cover
        emb = model.encode(texts, normalize_embeddings=True)
        return np.asarray(emb, dtype=np.float32)
    return np.vstack([_hash_embed(t) for t in texts])


def mode() -> str:
    _try_real_model()
    return _mode
