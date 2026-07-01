"""
pMHC-I immunogenicity (B10) via MHCflurry 2.x (Apache-2.0).

Scores peptide presentation by class-I HLA alleles — a developability axis for designed
proteins/biologics (immunogenic-epitope burden) and a neoantigen pre-screen. Real model
(MHCflurry presentation predictor); a sliding-window scan turns a full sequence into a
per-allele epitope-burden score.
"""
from __future__ import annotations

from typing import Any

_PRED = None


def _predictor():
    global _PRED
    if _PRED is None:
        from mhcflurry import Class1PresentationPredictor
        _PRED = Class1PresentationPredictor.load()
    return _PRED


def predict(peptides: list[str], alleles: list[str]) -> list[dict[str, Any]]:
    """Presentation scores for explicit peptides against a set of HLA alleles."""
    df = _predictor().predict(peptides=peptides, alleles=alleles, verbose=0)
    return df.to_dict("records")


def scan_sequence(sequence: str, alleles: list[str], length: int = 9,
                  threshold: float = 0.5) -> dict[str, Any]:
    """Sliding-window epitope scan → immunogenic burden (fraction of strong-presented k-mers)."""
    seq = "".join(sequence.split()).upper()
    peptides = [seq[i:i + length] for i in range(max(0, len(seq) - length + 1))]
    if not peptides:
        return {"sequence_len": len(seq), "n_peptides": 0, "alleles": alleles,
                "n_strong_binders": 0, "immunogenic_burden": 0.0, "top": []}
    recs = predict(peptides, alleles)
    strong = [r for r in recs if r.get("presentation_score", 0.0) >= threshold]
    top = sorted(recs, key=lambda r: -r.get("presentation_score", 0.0))[:5]
    return {
        "sequence_len": len(seq),
        "n_peptides": len(peptides),
        "alleles": alleles,
        "n_strong_binders": len(strong),
        "immunogenic_burden": round(len(strong) / len(peptides), 3),
        "top": [{"peptide": r.get("peptide"), "best_allele": r.get("best_allele"),
                 "presentation_score": round(r.get("presentation_score", 0.0), 4)} for r in top],
    }
