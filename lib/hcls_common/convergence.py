"""
Multi-omics ConvergenceReasoner (connective-tissue §7.2) — the engine behind the factory's
"same signal across omics" story.

It reads a ``PatientContext`` (PF-1), resolves entities to canonical ``BioKey``s (PF-5) so aliases
converge, finds the entities that light up across more than one omics layer, and ranks each
convergence by **breadth × clinical weight × honesty floor**. The load-bearing rule is the honesty
floor: a convergence is **never more confident than its weakest contributing layer** — its label is
*computed, not copied*. A convergence whose weakest layer is ``hypothesis_only`` is presented as a
**hypothesis**, never a finding.

Scope note (honest): the fully-principled score is a *combination of likelihood ratios* per layer.
Calibrated LRs are clinical data (not fabricated here); this provides the breadth × honesty-floor
scaffold with an optional per-signal LR hook (``layer_lrs``) that slots real LRs in later. The
honesty floor — the part that must never be wrong — is exact.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from hcls_common.artifact import Honesty, Maturity, weakest_maturity
from hcls_common.biokey import BioKeyResolver, DEFAULT_RESOLVER
from hcls_common.multiomics import PatientContext

# how much a convergence's score is discounted by its honesty floor (a hypothesis weighs less
# than a validated finding). Not a likelihood ratio — a transparent floor weight.
_FLOOR_WEIGHT = {
    Maturity.live: 1.0,
    Maturity.preclinical: 0.7,
    Maturity.research_use: 0.5,
    Maturity.hypothesis_only: 0.3,
    Maturity.roadmap: 0.1,
}


@dataclass
class ConvergenceSignal:
    entity: str                 # the canonical entity (gene/BioKey id) that converges
    layers: list[str]           # the omics layers it appears in
    breadth: int                # number of layers (the multi-modal breadth)
    honesty: Honesty            # floor = weakest contributing layer (computed, not copied)
    score: float                # breadth × floor-weight (× LR product if supplied)

    @property
    def is_finding(self) -> bool:
        """True only if every contributing layer is at least ``live``; otherwise this is presented
        as a *hypothesis* (the weakest layer carries caution)."""
        return self.honesty.maturity is Maturity.live

    @property
    def presentation(self) -> str:
        return "finding" if self.is_finding else "hypothesis"

    def to_dict(self) -> dict[str, Any]:
        return {"entity": self.entity, "layers": list(self.layers), "breadth": self.breadth,
                "honesty": self.honesty.to_dict(), "score": self.score,
                "presentation": self.presentation}


class ConvergenceReasoner:
    """Ranks cross-omics convergence for one patient, honesty-floored.

    ``layer_maturity`` declares each layer's clinical weight (e.g. a single-cell drug-response
    layer is ``hypothesis_only``); layers absent from the map default to ``live``. Pass a PF-5
    ``resolver`` so entity aliases (``tuberin`` ≡ ``TSC2``) converge instead of missing.
    """

    def __init__(self, resolver: BioKeyResolver | None = None,
                 layer_maturity: dict[str, Maturity] | None = None) -> None:
        self.resolver = resolver or DEFAULT_RESOLVER
        self.layer_maturity = dict(layer_maturity or {})

    def _floor(self, layers: list[str]) -> Maturity:
        return weakest_maturity([self.layer_maturity.get(l, Maturity.live) for l in layers])

    def signals(self, ctx: PatientContext,
                layer_lrs: dict[str, dict[str, float]] | None = None) -> list[ConvergenceSignal]:
        """Ranked convergence signals for the patient. ``layer_lrs[entity][layer]`` optionally
        supplies a calibrated likelihood ratio per (entity, layer); when present the score is
        ``floor_weight × Π(LRs)`` (a principled combination), else ``breadth × floor_weight``."""
        out: list[ConvergenceSignal] = []
        for link in ctx.cross_omics_links(resolver=self.resolver):
            entity, layers, breadth = link["gene"], link["layers"], link["n_layers"]
            floor = self._floor(layers)
            honesty = Honesty(maturity=floor, labels=["decision-support"], requires=["clinician-review"])
            fw = _FLOOR_WEIGHT[floor]
            if layer_lrs and entity in layer_lrs:
                lr_product = 1.0
                for l in layers:
                    lr_product *= layer_lrs[entity].get(l, 1.0)
                score = round(fw * lr_product, 4)
            else:
                score = round(breadth * fw, 4)
            out.append(ConvergenceSignal(entity, layers, breadth, honesty, score))
        out.sort(key=lambda s: (-s.score, -s.breadth, s.entity))
        return out
