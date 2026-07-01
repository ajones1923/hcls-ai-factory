"""
Deterministic TAND aggregation and scoring (PRD §3 FR-TS-2/4; master paper §11).

A tunable, transparent scoring layer (NOT an LLM), configured in config/tand_scoring.yaml:
per-marker weights, recency decay, a recurrence bonus, a negation filter, and a flag
threshold. It scores from the note's ground-truth discourse SPANS when present — honoring
polarity so a DENIED concern ("mother denies…") is never counted — and falls back to
substring detection for spanless notes (e.g. external or LLM-added text). Flags clusters
whose weighted longitudinal score crosses the configured threshold.
"""
from __future__ import annotations

import yaml

from config.settings import settings
from src.agents.tand_surveillance.taxonomy import (
    TAND_CLUSTERS, detect_clusters, detect_markers,
)

_CFG: dict | None = None


def _cfg() -> dict:
    global _CFG
    if _CFG is None:
        with open(settings.TAND_SCORING_PATH) as f:
            _CFG = yaml.safe_load(f)
    return _CFG


def _recency_factor(date: str | None, cfg: dict) -> float:
    rc = cfg.get("recency", {})
    if not rc.get("enabled") or not date:
        return 1.0
    try:
        ny, nm = int(date[:4]), int(date[5:7])
        ry, rm = int(rc["reference_date"][:4]), int(rc["reference_date"][5:7])
        months_back = max(0, (ry - ny) * 12 + (rm - nm))
        return 0.5 ** (months_back / float(rc.get("half_life_months", 18)))
    except Exception:
        return 1.0


def _marker_weight(markers, cfg) -> float:
    mw = cfg["marker_weights"]
    return sum(mw.get(m, 0.0) for m in markers) or cfg["base_mention_weight"]


def score(signals: list[dict], notes: list[dict]) -> dict:
    """Return per-cluster weighted score + flagged clusters + highlighted passages."""
    cfg = _cfg()
    negate = cfg.get("negation", True)
    clusters = {c: {"score": 0.0, "mentions": 0, "sources": [], "markers": set()}
                for c in TAND_CLUSTERS}
    highlights = []

    # 1) explicit cluster-tagged structured signals (each a base-weight mention)
    for s in signals:
        c = s.get("cluster")
        if c in clusters:
            clusters[c]["score"] += cfg["base_mention_weight"]
            clusters[c]["mentions"] += 1
            clusters[c]["sources"].append(s.get("source"))

    # 2) notes: prefer ground-truth spans (polarity-aware), else substring fallback
    for note in notes:
        src = f"{note.get('specialty')} {note.get('date')}"
        rf = _recency_factor(note.get("date"), cfg)
        tand_spans = [sp for sp in (note.get("spans") or []) if sp.get("kind") == "tand"]
        if tand_spans:
            per_cluster: dict[str, dict] = {}
            for sp in tand_spans:
                if negate and sp.get("polarity") == "absent":
                    continue                                  # negation: a denied concern is not a signal
                c = sp.get("cluster")
                if c not in clusters:
                    continue
                pc = per_cluster.setdefault(c, {"w": 0.0, "markers": set(), "spans": []})
                pc["w"] += _marker_weight(sp.get("markers", []), cfg)
                pc["markers"].update(sp.get("markers", []))
                pc["spans"].append(sp)
            for c, pc in per_cluster.items():
                clusters[c]["score"] += pc["w"] * rf
                clusters[c]["mentions"] += 1
                clusters[c]["sources"].append(src)
                clusters[c]["markers"].update(pc["markers"])
                highlights.append({
                    "cluster": c, "markers": sorted(pc["markers"]),
                    "passage": note.get("text", ""), "source": src,
                    "spans": [{"start": s["start"], "end": s["end"], "quote": s["quote"]} for s in pc["spans"]],
                })
        else:
            text = note.get("text", "")
            markers = detect_markers(text)
            for c in detect_clusters(text):
                clusters[c]["score"] += _marker_weight(markers, cfg) * rf
                clusters[c]["mentions"] += 1
                clusters[c]["sources"].append(src)
                clusters[c]["markers"].update(markers)
                if markers:
                    highlights.append({"cluster": c, "markers": markers, "passage": text, "source": src})

    # 3) recurrence bonus + flagging
    rb = cfg.get("recurrence_bonus", 0.0)
    thr = cfg["cluster_flag_threshold"]
    flagged = []
    for c, v in clusters.items():
        if v["mentions"] > 1:
            v["score"] += rb * (v["mentions"] - 1)
        v["score"] = round(v["score"], 3)
        v["markers"] = sorted(v["markers"])
        v["signal_density"] = v["mentions"]          # back-compat key
        if v["score"] >= thr:
            flagged.append(c)
    return {"clusters": clusters, "flagged_clusters": flagged, "highlights": highlights}
