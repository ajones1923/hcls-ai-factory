"""
ESM-2 protein embeddings (B2) — turns a protein sequence into a fixed-length vector
for similarity search.

Model: facebook/esm2_t33_650M_UR50D (1280-dim). Embedding is a single forward pass
(no folding trunk) so it is fast even on CPU. Vectors are mean-pooled over residues
(excluding the BOS/EOS special tokens) and L2-normalized, so inner-product == cosine.

Pure logic (pooling, normalization) is separated from the heavy model load and the
model is injectable, so this is unit-testable without the 2.5 GB download.
"""
from __future__ import annotations

import math
from typing import Callable, Sequence

from esmfold_service import validate_sequence  # reuse the shared AA validator

ESM2_MODEL = "facebook/esm2_t33_650M_UR50D"
ESM2_DIM = 1280


def l2_normalize(vec: Sequence[float]) -> list[float]:
    norm = math.sqrt(sum(x * x for x in vec)) or 1.0
    return [x / norm for x in vec]


class ESM2Embedder:
    """Lazy ESM-2 embedder. `_embed_fn` (callable(seq)->raw vector) overrides the real model for tests."""

    def __init__(
        self,
        model_name: str = ESM2_MODEL,
        dim: int = ESM2_DIM,
        max_len: int = 1024,
        _embed_fn: Callable[[str], Sequence[float]] | None = None,
    ) -> None:
        self.model_name = model_name
        self.dim = dim
        self.max_len = max_len
        self._embed_fn = _embed_fn
        self._tok = None
        self._model = None
        self._torch = None

    def _load(self) -> None:
        if self._model is None:
            import torch
            from transformers import AutoTokenizer, EsmModel
            self._torch = torch
            self.device = "cuda" if torch.cuda.is_available() else "cpu"   # GB10 when available
            self._tok = AutoTokenizer.from_pretrained(self.model_name)
            self._model = EsmModel.from_pretrained(self.model_name).eval().to(self.device)

    def _raw_embed(self, seq: str) -> Sequence[float]:
        if self._embed_fn is not None:
            return self._embed_fn(seq)
        self._load()
        enc = self._tok(seq, return_tensors="pt", add_special_tokens=True).to(self.device)
        with self._torch.no_grad():
            hidden = self._model(**enc).last_hidden_state[0]   # [L+2, dim]
        # mean-pool residues only (drop BOS at 0 and EOS at -1)
        return hidden[1:-1].mean(0).cpu().tolist()

    def embed(self, sequence: str) -> list[float]:
        """Validated, mean-pooled, L2-normalized embedding (cosine-ready)."""
        seq = validate_sequence(sequence, self.max_len)
        vec = self._raw_embed(seq)
        if len(vec) != self.dim:
            raise ValueError(f"embedding dim {len(vec)} != expected {self.dim}")
        return l2_normalize(vec)

    def embed_batch(self, sequences: list[str]) -> list[list[float]]:
        return [self.embed(s) for s in sequences]
