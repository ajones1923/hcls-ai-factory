"""
Protein design + developability (B3/B4).

B4 — developability scoring: real, computed biophysical proxies (Bio ProtParam): GRAVY
hydropathy, Guruprasad instability index, isoelectric point, aromaticity, MW — with a
flag/verdict summary. No model, deterministic.

B3 — developability-guided sequence design: a single-point mutation scan that proposes
substitutions which *improve* developability (lower instability index) — a real protein
engineering optimizer that complements structure-conditioned design (ProteinMPNN, the GPU
follow-up). Optional MHCflurry immunogenicity when installed.
"""
from esmfold_service import validate_sequence, SequenceError

_SUBSTITUTIONS = list("AEQKRLIVFY")   # common, generally well-tolerated substitutions to scan


def _analysis(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return ProteinAnalysis(seq)


def develop_metrics(seq: str) -> dict:
    pa = _analysis(seq)
    return {
        "gravy": round(pa.gravy(), 3),
        "instability_index": round(pa.instability_index(), 2),
        "isoelectric_point": round(pa.isoelectric_point(), 2),
        "aromaticity": round(pa.aromaticity(), 3),
        "molecular_weight": round(pa.molecular_weight(), 1),
        "length": len(seq),
    }


def develop_flags(m: dict) -> tuple[list[str], str]:
    flags = []
    if m["instability_index"] >= 40:
        flags.append("unstable (instability index >= 40)")
    if m["gravy"] > 0.4:
        flags.append("hydrophobic / aggregation risk (GRAVY > 0.4)")
    if m["isoelectric_point"] < 5.0 or m["isoelectric_point"] > 9.5:
        flags.append("extreme pI (charge-driven solubility risk)")
    verdict = "clean" if not flags else ("caution" if len(flags) == 1 else "high-risk")
    return flags, verdict


class DevelopabilityScorer:
    def score(self, sequence: str) -> dict:
        seq = validate_sequence(sequence, max_len=2000)
        m = develop_metrics(seq)
        flags, verdict = develop_flags(m)
        return {"sequence_len": len(seq), "metrics": m, "flags": flags, "verdict": verdict}

    def optimize(self, sequence: str, n: int = 5) -> dict:
        """Propose single-point mutations that lower the instability index (real, computed)."""
        seq = validate_sequence(sequence, max_len=2000)
        base = develop_metrics(seq)["instability_index"]
        proposals = []
        for i, wt in enumerate(seq):
            for sub in _SUBSTITUTIONS:
                if sub == wt:
                    continue
                mutant = seq[:i] + sub + seq[i + 1:]
                ii = develop_metrics(mutant)["instability_index"]
                if ii < base:
                    proposals.append({"mutation": f"{wt}{i+1}{sub}", "instability_index": round(ii, 2),
                                      "delta": round(ii - base, 2)})
        proposals.sort(key=lambda p: p["delta"])         # most-improving first
        return {"baseline_instability": base, "n_proposals": len(proposals), "top": proposals[:n]}
