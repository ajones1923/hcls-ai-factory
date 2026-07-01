"""
ACMG-AMP combinatorial classifier (PRD §3 FR-VC-5/6; master paper §8).

Implements the Richards et al. 2015 (Genet Med) Table 5 combining rules. Criteria are
assigned with an explicit strength BUCKET (PVS/PS/PM/PP for pathogenic evidence;
BA/BS/BP for benign), and this validator maps the assembled buckets to one of the five
classifications, returning the rule that fired so the call is auditable. The Opus
synthesis step proposes per-criterion reasoning; THIS function is authoritative for the
classification, exactly as the PRD requires ("validated against the standard ACMG-AMP
combinatorial rules; if not satisfied, the classification is rejected").
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

PATH_BUCKETS = {"PVS", "PS", "PM", "PP"}
BENIGN_BUCKETS = {"BA", "BS", "BP"}


@dataclass
class Criterion:
    code: str          # e.g. "PVS1", "PM2", "PP4", "PS1", "BA1"
    bucket: str        # PVS | PS | PM | PP | BA | BS | BP
    rationale: str = ""


def _pathogenic(pvs: int, ps: int, pm: int, pp: int) -> str | None:
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
        return "1 Very Strong (PVS1) + supporting combination"
    if ps >= 2:
        return "≥2 Strong"
    if ps >= 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)):
        return "1 Strong + Moderate/Supporting combination"
    return None


def _likely_pathogenic(pvs: int, ps: int, pm: int, pp: int) -> str | None:
    if pvs >= 1 and pm >= 1:
        return "1 Very Strong + 1 Moderate"
    if ps >= 1 and 1 <= pm <= 2:
        return "1 Strong + 1-2 Moderate"
    if ps >= 1 and pp >= 2:
        return "1 Strong + ≥2 Supporting"
    if pm >= 3:
        return "≥3 Moderate"
    if pm >= 2 and pp >= 2:
        return "2 Moderate + ≥2 Supporting"
    if pm >= 1 and pp >= 4:
        return "1 Moderate + ≥4 Supporting"
    return None


def _benign(ba: int, bs: int) -> str | None:
    if ba >= 1:
        return "BA1 stand-alone"
    if bs >= 2:
        return "≥2 Strong (benign)"
    return None


def _likely_benign(bs: int, bp: int) -> str | None:
    if bs >= 1 and bp >= 1:
        return "1 Strong + 1 Supporting (benign)"
    if bp >= 2:
        return "≥2 Supporting (benign)"
    return None


def classify(criteria: list[Criterion]) -> tuple[str, str]:
    """Return (classification, rule_that_fired)."""
    c = Counter(cr.bucket for cr in criteria)
    pvs, ps, pm, pp = c["PVS"], c["PS"], c["PM"], c["PP"]
    ba, bs, bp = c["BA"], c["BS"], c["BP"]

    path_rule = _pathogenic(pvs, ps, pm, pp)
    lp_rule = _likely_pathogenic(pvs, ps, pm, pp)
    ben_rule = _benign(ba, bs)
    lb_rule = _likely_benign(bs, bp)

    path_side = path_rule or lp_rule
    benign_side = ben_rule or lb_rule
    if path_side and benign_side:
        return "Variant of Uncertain Significance", "conflicting pathogenic and benign evidence"
    if path_rule:
        return "Pathogenic", path_rule
    if lp_rule:
        return "Likely Pathogenic", lp_rule
    if ben_rule:
        return "Benign", ben_rule
    if lb_rule:
        return "Likely Benign", lb_rule
    return "Variant of Uncertain Significance", "criteria insufficient for a confident call"
