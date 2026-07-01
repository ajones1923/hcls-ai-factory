"""
Exact Smith-Waterman re-rank stage for protein search (B3).

Upgrades sequence search from pure-ANN recall to the retrieve-then-rerank pattern: ANN
recalls candidates fast (cosine over ESM-2 embeddings), then exact local alignment re-scores
the top candidates for precision (%identity + alignment score). Uses Biopython's
PairwiseAligner (BLOSUM62 local SW) — pure, ARM-clean, permissively licensed. (parasail's SIMD
SW has no aarch64 wheel on this hardware; Biopython gives the same exact alignment, slower,
which is fine for re-ranking a small candidate set.)
"""
from __future__ import annotations

from typing import Any

_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        a = PairwiseAligner()
        a.mode = "local"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
        _ALIGNER = a
    return _ALIGNER


def sw_align(query: str, target: str) -> dict[str, Any]:
    """Exact local (Smith-Waterman) alignment: score + %identity over aligned columns."""
    try:
        aln = _aligner().align(query, target)[0]  # best local alignment
    except (IndexError, StopIteration):           # no positive-scoring local region
        return {"sw_score": 0.0, "pct_identity": 0.0, "aln_length": 0, "matches": 0}
    top, bottom = str(aln[0]), str(aln[1])        # gapped aligned strings
    cols = max(len(top), 1)
    matches = sum(1 for a, b in zip(top, bottom) if a == b and a != "-")
    return {
        "sw_score": round(float(aln.score), 1),
        "pct_identity": round(matches / cols, 4),
        "aln_length": cols,
        "matches": matches,
    }


def rerank(query: str, hits: list[dict[str, Any]], seq_lookup: dict[str, str]) -> list[dict[str, Any]]:
    """Re-score ANN hits by exact SW alignment and re-sort by alignment score (desc)."""
    out = []
    for h in hits:
        seq = seq_lookup.get(h.get("name"))
        out.append({**h, **sw_align(query, seq)} if seq else {**h, "sw_score": None})
    out.sort(key=lambda h: (h.get("sw_score") if h.get("sw_score") is not None else -1.0), reverse=True)
    return out
