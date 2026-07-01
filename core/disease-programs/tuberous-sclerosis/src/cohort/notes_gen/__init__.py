"""
Synthetic clinical-note generator with ground-truth spans (PRD §3 FR-CG-3; master paper §15).

Produces a deterministic set of realistic, multi-specialty notes per patient and — crucially
— records the exact character-offset span of every embedded HPO phenotype and TAND discourse
marker, with polarity (present/absent) and temporality (onset/current/historical). Those
spans are the ground truth the Phenome Mapper (FR-PM-2: verbatim-span extraction) and the
TAND agent (FR-TS-1/3: marker offsets + negation handling) consume — so the offline build has
real, span-grounded NLP signal rather than one template sentence. Frontier-model elaboration
(when keyed) layers on top; these spans remain the checkable substrate.

Deterministic by patient id (stable cohort hash). Every `quote` equals `text[start:end]`.
"""
from __future__ import annotations

import random

from src.agents.tand_surveillance.taxonomy import DISCOURSE_MARKERS

# organ routing: which specialty note an HPO feature is mentioned in
_SPECIALTY_FOR = {
    "HP:0009717": "neurology",    # cortical tubers
    "HP:0009718": "neurology",    # SEGA
    "HP:0001250": "neurology",    # seizure
    "HP:0012469": "neurology",    # infantile spasms
    "HP:0006772": "nephrology",   # renal AML
    "HP:0009719": "dermatology",  # hypomelanotic macule
    "HP:0009720": "dermatology",  # adenoma sebaceum
    "HP:0001249": "developmental",  # intellectual disability
    "HP:0000717": "developmental",  # autism
}
_DATES = ["2024-03-12", "2024-09-08", "2025-02-19", "2025-06-03", "2025-09-15", "2025-11-20"]


class _NoteBuilder:
    """Accumulates note text while recording verbatim spans at their true offsets."""

    def __init__(self) -> None:
        self._parts: list[str] = []
        self._spans: list[dict] = []
        self._pos = 0

    def add(self, text: str) -> None:
        self._parts.append(text)
        self._pos += len(text)

    def mark(self, text: str, **span) -> None:
        start = self._pos
        self.add(text)
        self._spans.append({"start": start, "end": self._pos, "quote": text, **span})

    def text(self) -> str:
        return "".join(self._parts)

    def spans(self) -> list[dict]:
        return self._spans


def _hpo_phrase(label: str) -> str:
    return label[0].lower() + label[1:] if label else label


def _temporality(onset: str | None) -> str:
    return "onset" if onset in ("infancy", "childhood") else "current"


def _note(specialty: str, date: str, builder: _NoteBuilder) -> dict:
    return {"note_id": f"{specialty}:{date}", "specialty": specialty, "date": date,
            "text": builder.text(), "spans": builder.spans()}


def build_notes(spec, phenos: list[dict], tand_clusters: dict | None) -> list[dict]:
    rng = random.Random(f"notes:{spec.patient_id}")
    tand_clusters = tand_clusters or {}
    by_specialty: dict[str, list[dict]] = {}
    for p in phenos:
        by_specialty.setdefault(_SPECIALTY_FOR.get(p["hpo_id"], "neurology"), []).append(p)

    notes: list[dict] = []
    di = 0

    # 1) genetics summary (molecular finding; no HPO span needed)
    b = _NoteBuilder()
    b.add(f"[SYNTHETIC] Genetics — {spec.age}yo {spec.sex}. ")
    if spec.gene:
        b.add(f"Confirmed {spec.gene} pathogenic variant "
              f"({'mosaic, affected-tissue' if spec.zygosity == 'mosaic' else 'germline'}). ")
    else:
        b.add("No mutation identified on blood; tissue sequencing advised. ")
    b.add("Meets clinical diagnostic criteria for TSC.")
    notes.append(_note("genetics", _DATES[di], b)); di += 1

    # 2..n) one note per involved specialty, embedding each HPO feature as a verbatim span
    for specialty in ("neurology", "nephrology", "dermatology", "developmental"):
        feats = by_specialty.get(specialty)
        if not feats:
            continue
        b = _NoteBuilder()
        b.add(f"[SYNTHETIC] {specialty.title()} — findings include ")
        for k, p in enumerate(feats):
            if k:
                b.add(", ")
            b.mark(_hpo_phrase(p["label"]), kind="hpo", hpo_id=p["hpo_id"], label=p["label"],
                   polarity="present", temporality=_temporality(p.get("onset")))
        b.add(". Continue ITSC surveillance.")
        if specialty == "nephrology" and spec.extras.get("egfr_series"):
            b.add(" Renal function trending down on serial panels.")
        notes.append(_note(specialty, _DATES[min(di, len(_DATES) - 1)], b)); di += 1

    # longitudinal discordance: a feature present at onset but resolved now (FR-PM-3). Infantile
    # spasms in infancy, absent currently — a present/absent conflict the Phenome Mapper logs.
    if any(p["hpo_id"] == "HP:0012469" for p in phenos):
        b = _NoteBuilder()
        b.add("[SYNTHETIC] Neurology follow-up — ")
        b.mark("no current infantile spasms", kind="hpo", hpo_id="HP:0012469",
               label="Infantile spasms", polarity="absent", temporality="current")
        b.add(" since infancy; other seizures stable.")
        notes.append(_note("neurology", _DATES[min(di, len(_DATES) - 1)], b)); di += 1

    # primary-care note: under-formalized TAND signal(s) with real discourse markers + a
    # NEGATION example (must NOT be counted) — the FR-TS-3 specificity case.
    b = _NoteBuilder()
    b.add(f"[SYNTHETIC] Primary care — {spec.age}yo for routine follow-up. ")
    for cluster, n in ((c, k) for c, k in tand_clusters.items() if k):
        for _ in range(n):
            tp = rng.choice(DISCOURSE_MARKERS["third_party_attribution"])
            hedge = rng.choice(DISCOURSE_MARKERS["hedging"])
            kw = {"academic": "focus at school", "behavioral": "behavior at home",
                  "psychiatric": "anxiety", "neuropsychological": "attention",
                  "intellectual": "developmental delay", "psychosocial": "social withdrawal",
                  }.get(cluster, "focus at school")
            defer = rng.choice(DISCOURSE_MARKERS["deferral"])
            b.mark(f"{tp.title()} {hedge} with {kw}", kind="tand", cluster=cluster,
                   markers=["third_party_attribution", "hedging"], polarity="present")
            b.add("; ")
            b.mark(defer if defer.startswith("will") else f"will {defer}",
                   kind="tand", cluster=cluster, markers=["deferral"], polarity="present")
            b.add(". ")
    # negation: a behavioral concern explicitly DENIED (polarity absent -> not counted)
    b.mark("Mother denies new behavioral concerns", kind="tand", cluster="behavioral",
           markers=["third_party_attribution"], polarity="absent")
    b.add(" at home.")
    notes.append(_note("primary care", _DATES[min(di, len(_DATES) - 1)], b)); di += 1

    return notes
