"""
Composed-chain lineage / reproducibility manifest (connective-tissue §8).

The existing ``ReproducibilityManifest`` captures a single node/run (hardware, pinned packages,
seeds, IO). But an MTB packet stands on a *composed* chain touching E1+E5+E4+E8+A2+A6+A3 — and there
is no single manifest covering every hop (and the RunPod leg). This walks the ``Artifact``
provenance DAG that ``WorkflowComposer.run(governed=True)`` produces into one manifest for the whole
chain:

- the **lineage graph** (nodes + edges) chained via ``provenance.inputs``,
- **where each hop ACTUALLY ran** (``serving``: native / cloud_nim / remote_gpu) — the elastic-burst
  transparency, honest across a RunPod burst,
- each hop's **honesty** (maturity + labels) and ``knowledge_version`` / ``as_of``,
- the **chain-level honesty floor** (the whole chain is never more confident than its weakest hop),
- a deterministic **replay hash**: the same chain *structure* (producers × shapes × input-producers)
  yields the same hash, independent of volatile uuids/timestamps — so a re-run is verifiable.

This is the mechanism behind a packet citing "the exact E1 variant call and E4 image read it stands on."
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from hcls_common.artifact import Artifact, combine_honesty


def _as_artifacts(artifacts: Any) -> list[Artifact]:
    """Accept a dict of {node: Artifact|dict}, a list of Artifacts, or their to_dict() forms."""
    items = list(artifacts.values()) if isinstance(artifacts, dict) else list(artifacts)
    return [a if isinstance(a, Artifact) else Artifact.from_dict(a) for a in items]


def chain_lineage(artifacts: Any) -> dict[str, Any]:
    """Walk a governed run's Artifacts into a composed-chain lineage/reproducibility manifest.

    ``artifacts`` is what ``WorkflowComposer.run(governed=True)`` returns under ``artifacts``
    (a dict of node -> Artifact or its dict form), or any iterable of Artifacts.
    """
    arts = _as_artifacts(artifacts)
    by_id = {a.id: a for a in arts}

    nodes = [{
        "id": a.id,
        "shape": a.shape.value,
        "producer": a.provenance.producer_id,
        "serving": a.provenance.serving,                 # where it ACTUALLY ran
        "knowledge_version": a.provenance.knowledge_version,
        "as_of": a.provenance.as_of,
        "maturity": a.honesty.maturity.value,
        "labels": list(a.honesty.labels),
        "patient_id": a.patient_id,
    } for a in arts]

    edges = [{"from": inp, "to": a.id}
             for a in arts for inp in a.provenance.inputs if inp in by_id]

    chain_honesty = combine_honesty([a.honesty for a in arts])

    # replay hash: structure-only (producer × shape × the producers of its inputs), so it is stable
    # across uuids/timestamps — an identical chain re-run hashes identically.
    struct = sorted(
        [a.provenance.producer_id, a.shape.value,
         sorted(by_id[i].provenance.producer_id for i in a.provenance.inputs if i in by_id)]
        for a in arts
    )
    lineage_hash = hashlib.sha256(json.dumps(struct, sort_keys=True).encode()).hexdigest()[:16]

    run_ids = sorted({a.run_id for a in arts if a.run_id})
    return {
        "run_id": run_ids[0] if len(run_ids) == 1 else (run_ids or [""])[0],
        "n_hops": len(arts),
        "nodes": nodes,
        "edges": edges,
        "chain_honesty": chain_honesty.to_dict(),        # the whole chain's honesty floor
        "lineage_hash": lineage_hash,
    }
