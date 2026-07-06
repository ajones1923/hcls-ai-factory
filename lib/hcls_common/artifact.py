"""
Artifact envelope (PF-3) — the universal, typed wrapper every payload wears when it
crosses a capability boundary, so **provenance and honesty can never be dropped at a seam**.

Today each of the seventeen capabilities carries its own ad-hoc `honesty_label` / `_meta` /
`honesty` field, and a payload crossing a hop loses its patient thread, its lineage, and its
honesty label. This module is the single object they collapse into:

    Artifact(shape, payload, patient_id, run_id, honesty, provenance)

- ``shape``      — the *semantic* ``ArtifactShape`` (PF-2): what the payload IS.
- ``honesty``   — a typed ``Honesty`` (maturity + labels + required reviews) that travels WITH
                  the data and is checked at every gate. The maturity ladder is **ordered**, so a
                  composite can be held to the *weakest leg* it derives from (the non-inflation law,
                  PF-4). This module provides the ordering + the combining primitive
                  (``combine_honesty``); PF-4 wires the enforcing gate into governance + composer.
- ``provenance``— the lineage edge (producer, where it ACTUALLY ran, and the input artifact ids),
                  so a composed run is replayable and a claim can cite the exact calls beneath it.

Stdlib-only, declarative dataclasses — consistent with the rest of hcls_common.
"""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from hcls_common.capability_registry import ArtifactShape


# --------------------------------------------------------------------------- #
# Honesty — a typed, ordered maturity + labels that ride with every artifact
# --------------------------------------------------------------------------- #
class Maturity(str, Enum):
    """How much clinical weight an artifact may carry — **ordered strongest → weakest**.

    The order is load-bearing: the non-inflation law (PF-4) forbids a downstream artifact from
    being *stronger* (more clinically weighty) than the weakest input it derives from — a
    ``hypothesis_only`` cell state may not become a ``live`` treatment recommendation.
    """
    live = "live"                      # real, validated for use as decision support
    preclinical = "preclinical"        # design/analysis bench, not a treatment today
    research_use = "research_use"      # research / trial-use input, not routine diagnosis
    hypothesis_only = "hypothesis_only"  # a hypothesis to test (e.g. drug-response prediction)
    roadmap = "roadmap"                # not built / not real yet

    @property
    def rank(self) -> int:
        """0 = strongest (live) … 4 = weakest (roadmap). Higher = more cautious."""
        return _MATURITY_ORDER.index(self)


_MATURITY_ORDER: tuple[Maturity, ...] = (
    Maturity.live,
    Maturity.preclinical,
    Maturity.research_use,
    Maturity.hypothesis_only,
    Maturity.roadmap,
)


@dataclass
class Honesty:
    """The honesty envelope that travels with an artifact and is checked at every gate."""
    maturity: Maturity
    labels: list[str] = field(default_factory=list)    # e.g. decision-support, pediatric-caution, gated, rna_inferred
    requires: list[str] = field(default_factory=list)  # e.g. clinician-review, tumor-board, orthogonal-confirmation
    rationale: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "maturity": self.maturity.value,
            "labels": list(self.labels),
            "requires": list(self.requires),
            "rationale": self.rationale,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Honesty":
        return cls(
            maturity=Maturity(d["maturity"]),
            labels=list(d.get("labels", [])),
            requires=list(d.get("requires", [])),
            rationale=d.get("rationale", ""),
        )


def weakest_maturity(maturities: list[Maturity]) -> Maturity:
    """The most cautious maturity in a set (the 'weakest leg'). Empty → live."""
    return max(maturities, key=lambda m: m.rank) if maturities else Maturity.live


def combine_honesty(inputs: list[Honesty], *, rationale: str = "") -> Honesty:
    """Non-inflating combine (PF-4's core rule): a composite's maturity is the *weakest* of its
    inputs, and its labels/requires are the **union** of all inputs — a claim is never more
    confident, and never less caveated, than the evidence it stands on. Empty → live."""
    if not inputs:
        return Honesty(maturity=Maturity.live, rationale=rationale)
    labels = sorted({l for h in inputs for l in h.labels})
    requires = sorted({r for h in inputs for r in h.requires})
    return Honesty(
        maturity=weakest_maturity([h.maturity for h in inputs]),
        labels=labels,
        requires=requires,
        rationale=rationale or f"non-inflated from {len(inputs)} input(s)",
    )


# --------------------------------------------------------------------------- #
# Provenance — the lineage edge
# --------------------------------------------------------------------------- #
@dataclass
class Provenance:
    producer_id: str               # registry capability id that emitted this
    invoke_path: str = "/"         # the endpoint path used
    serving: str = "native"        # where it ACTUALLY ran (native | cloud_nim | remote_gpu ...)
    inputs: list[str] = field(default_factory=list)  # artifact ids this derived from (lineage edge)
    knowledge_version: str = ""    # guideline/atlas/db snapshot tag
    as_of: str = ""                # data-currency stamp
    ts: str = ""                   # emission timestamp (ISO-8601)

    def to_dict(self) -> dict[str, Any]:
        return {
            "producer_id": self.producer_id,
            "invoke_path": self.invoke_path,
            "serving": self.serving,
            "inputs": list(self.inputs),
            "knowledge_version": self.knowledge_version,
            "as_of": self.as_of,
            "ts": self.ts,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Provenance":
        return cls(
            producer_id=d["producer_id"],
            invoke_path=d.get("invoke_path", "/"),
            serving=d.get("serving", "native"),
            inputs=list(d.get("inputs", [])),
            knowledge_version=d.get("knowledge_version", ""),
            as_of=d.get("as_of", ""),
            ts=d.get("ts", ""),
        )


# --------------------------------------------------------------------------- #
# Artifact — the envelope itself
# --------------------------------------------------------------------------- #
@dataclass
class Artifact:
    id: str
    shape: ArtifactShape
    payload: dict[str, Any]
    honesty: Honesty
    provenance: Provenance
    patient_id: str | None = None   # the multi-omics join key (PF-5 resolves it)
    run_id: str = ""                # the composed-run thread

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "shape": self.shape.value,
            "payload": self.payload,
            "honesty": self.honesty.to_dict(),
            "provenance": self.provenance.to_dict(),
            "patient_id": self.patient_id,
            "run_id": self.run_id,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Artifact":
        return cls(
            id=d["id"],
            shape=ArtifactShape(d["shape"]),
            payload=dict(d.get("payload", {})),
            honesty=Honesty.from_dict(d["honesty"]),
            provenance=Provenance.from_dict(d["provenance"]),
            patient_id=d.get("patient_id"),
            run_id=d.get("run_id", ""),
        )


def new_artifact(
    shape: ArtifactShape,
    payload: dict[str, Any],
    *,
    producer_id: str,
    honesty: Honesty | None = None,
    maturity: Maturity = Maturity.live,
    labels: list[str] | None = None,
    requires: list[str] | None = None,
    inputs: list[str] | None = None,
    invoke_path: str = "/",
    serving: str = "native",
    knowledge_version: str = "",
    as_of: str = "",
    patient_id: str | None = None,
    run_id: str = "",
    id: str | None = None,
    ts: str | None = None,
) -> Artifact:
    """Convenience constructor — fills a uuid id and an ISO timestamp, builds the Provenance,
    and either takes a supplied ``honesty`` or builds one from ``maturity``/``labels``/``requires``."""
    h = honesty or Honesty(
        maturity=maturity, labels=list(labels or []), requires=list(requires or [])
    )
    prov = Provenance(
        producer_id=producer_id,
        invoke_path=invoke_path,
        serving=serving,
        inputs=list(inputs or []),
        knowledge_version=knowledge_version,
        as_of=as_of,
        ts=ts if ts is not None else datetime.now(timezone.utc).isoformat(),
    )
    return Artifact(
        id=id or uuid.uuid4().hex,
        shape=shape,
        payload=payload,
        honesty=h,
        provenance=prov,
        patient_id=patient_id,
        run_id=run_id,
    )


# --------------------------------------------------------------------------- #
# PF-4 — the non-inflation gate (the load-bearing wall at the seams)
# --------------------------------------------------------------------------- #
def non_inflation_issues(artifact: Artifact, inputs: list[Artifact]) -> list[str]:
    """The non-inflation gate: a downstream artifact may not be **more clinically weighty than
    the weakest input it derives from**, and it may not **drop** any caution its inputs carry.

    Returns a list of blocking ``ERROR:`` issues (empty = clean) — the honesty label is computed
    from the evidence, never inflated as it crosses a seam. Wired into the composer/output gate in
    the integration wave; usable standalone now.
    """
    issues: list[str] = []
    if not inputs:
        return issues
    weakest = weakest_maturity([a.honesty.maturity for a in inputs])
    if artifact.honesty.maturity.rank < weakest.rank:
        issues.append(
            f"ERROR: honesty inflation — artifact {artifact.id} claims "
            f"'{artifact.honesty.maturity.value}' but derives from weaker '{weakest.value}' "
            f"input(s); a claim may not be more confident than its evidence"
        )
    dropped_labels = {l for a in inputs for l in a.honesty.labels} - set(artifact.honesty.labels)
    dropped_requires = {r for a in inputs for r in a.honesty.requires} - set(artifact.honesty.requires)
    if dropped_labels:
        issues.append(f"ERROR: dropped honesty labels {sorted(dropped_labels)} carried by inputs")
    if dropped_requires:
        issues.append(f"ERROR: dropped required reviews {sorted(dropped_requires)} carried by inputs")
    return issues


def derive_artifact(
    shape: ArtifactShape,
    payload: dict[str, Any],
    *,
    producer_id: str,
    inputs: list[Artifact],
    own_maturity: Maturity = Maturity.live,
    extra_labels: list[str] | None = None,
    extra_requires: list[str] | None = None,
    invoke_path: str = "/",
    serving: str = "native",
    knowledge_version: str = "",
    as_of: str = "",
    run_id: str = "",
    id: str | None = None,
    ts: str | None = None,
) -> Artifact:
    """Build a downstream artifact whose honesty is **non-inflated by construction**: its maturity
    is capped at the weakest input (a producer may add its own caution with ``own_maturity`` but can
    never claim more confidence), input labels/requires are carried forward and unioned, its
    ``provenance.inputs`` chains the input artifact ids, and ``patient_id`` is inherited if the
    inputs agree. The result always passes ``non_inflation_issues``."""
    base = combine_honesty([a.honesty for a in inputs])
    maturity = weakest_maturity([base.maturity, own_maturity])   # own_maturity may only weaken
    labels = sorted(set(base.labels) | set(extra_labels or []))
    requires = sorted(set(base.requires) | set(extra_requires or []))
    pids = {a.patient_id for a in inputs if a.patient_id is not None}
    patient_id = next(iter(pids)) if len(pids) == 1 else None
    return new_artifact(
        shape, payload, producer_id=producer_id,
        honesty=Honesty(maturity=maturity, labels=labels, requires=requires,
                        rationale=f"non-inflated from {len(inputs)} input(s)"),
        inputs=[a.id for a in inputs],
        invoke_path=invoke_path, serving=serving,
        knowledge_version=knowledge_version, as_of=as_of,
        patient_id=patient_id, run_id=run_id or (inputs[0].run_id if inputs else ""),
        id=id, ts=ts,
    )
