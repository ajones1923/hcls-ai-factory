"""
The Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy and the six
TAND clusters (PRD §3 FR-TS-3; master paper §11). Versioned configuration — a direct
extension of Dr. Hagedorn's published clinical-NLP methodology to TSC neuropsychiatry.

The marker patterns are deliberately conservative substrings; the LLM (Sonnet) layer
adds nuance over these when keyed. This taxonomy is the versioned artifact a future
validation study would refine (a stated W5 / institutional deliverable).
"""
from __future__ import annotations

TAXONOMY_VERSION = "marshall-hagedorn/v0.1"

TAND_CLUSTERS = [
    "behavioral", "psychiatric", "intellectual",
    "academic", "neuropsychological", "psychosocial",
]

# diagnostic-uncertainty discourse markers (the methodology's core categories)
DISCOURSE_MARKERS = {
    "hedging": ["some difficulty", "mild", "possible", "seems", "appears", "reportedly", "a bit"],
    "deferral": ["will reassess", "follow up", "next visit", "monitor", "revisit", "watch"],
    "third_party_attribution": ["mother reports", "father reports", "teacher concerned",
                                "parent notes", "per mother", "per father", "family reports"],
    "conditional": ["if it persists", "may need", "consider", "would warrant"],
    "follow_up_without_formalization": ["will reassess next visit", "plan to revisit",
                                        "reassess at next", "no formal"],
}

# keyword -> TAND cluster (for deterministic per-note cluster assignment offline)
CLUSTER_KEYWORDS = {
    "academic": ["school", "focus", "grade", "learning", "academic", "homework"],
    "behavioral": ["aggression", "tantrum", "behavior", "sleep", "self-injur"],
    "psychiatric": ["anxiety", "anxious", "depress", "mood", "panic"],
    "intellectual": ["developmental delay", "intellectual disability", "cognitive delay"],
    "neuropsychological": ["attention", "memory", "executive function", "concentration"],
    "psychosocial": ["friend", "social", "isolation", "peer", "withdraw"],
}


def detect_markers(text: str) -> list[str]:
    t = text.lower()
    return [cat for cat, pats in DISCOURSE_MARKERS.items() if any(p in t for p in pats)]


def detect_clusters(text: str) -> list[str]:
    t = text.lower()
    return [cl for cl, kws in CLUSTER_KEYWORDS.items() if any(k in t for k in kws)]
