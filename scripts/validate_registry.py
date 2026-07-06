#!/usr/bin/env python3
"""Validate the capability registry and its coverage of the engine/agent tree.

Run in CI to prevent drift as new engines/agents are added:
  * the manifest parses and every capability passes the registry's rules
    (valid enums, unique IDs, and the honesty rule: a LIVE capability may
    never be MOCK-served);
  * every directory under core/engines/ and core/agents/ is represented by
    at least one registered capability.

Exit code 0 = clean, 1 = problems (prints them). Requires `hcls_common`
importable (pip install -e lib/hcls_common).
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Known component -> capability-id coverage. When you add an engine/agent, add
# its directory here pointing at the capability id(s) that represent it (or give
# the capability a tag/id/name that contains the directory name).
COVERAGE: dict[str, list[str]] = {
    # engines
    "core/engines/genomic-foundation": ["genomics-engine", "variant-store", "acmg-secondary-findings", "gwas-association"],
    "core/engines/precision-intelligence": ["precision-intelligence-engine"],
    "core/engines/therapeutic-discovery": ["therapeutic-discovery-engine", "molmim-nim", "diffdock-nim", "genmol-nim", "chemprop-admet", "molecule-generator"],
    "core/engines/clinical-imaging": ["imaging-intelligence-agent"],
    "core/engines/precision-oncology": ["precision-oncology-agent"],
    "core/engines/cardiology": ["cardiology-intelligence-agent"],
    "core/engines/structural-biology": ["esmfold-model", "esm2-search", "protein-developability", "mhcflurry-immunogenicity", "proteinmpnn-design", "esm2-finetune", "chai1-structure"],
    "core/engines/single-cell": ["singlecell-compute"],
    # agents
    "core/agents/cart": ["cart-intelligence-agent"],
    "core/agents/precision-biomarker": ["precision-biomarker-agent"],
    "core/agents/pharmacogenomics": ["pharmacogenomics-intelligence-agent"],
    "core/agents/precision-autoimmune": ["precision-autoimmune-agent"],
    "core/agents/neurology": ["neurology-intelligence-agent"],
    "core/agents/clinical-trial": ["clinical-trial-intelligence-agent"],
    "core/agents/rare-disease-diagnostic": ["rare-disease-diagnostic-agent"],
    "core/agents/single-cell": ["single-cell-intelligence-agent"],
}


def main() -> int:
    errors: list[str] = []

    # 1) Manifest parses + registry rules (enums, unique IDs, honesty rule).
    try:
        from hcls_common.capability_registry import CapabilityRegistry
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: cannot import hcls_common ({e}). Run: pip install -e lib/hcls_common")
        return 1
    try:
        reg = CapabilityRegistry().load_manifest()
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: manifest failed to load/validate: {e}")
        return 1

    ids = reg.ids()
    if len(ids) != len(set(ids)):
        errors.append("duplicate capability ids in manifest")
    id_set = set(ids)

    # 2) Coverage: every engine/agent directory maps to a registered capability.
    def dirs(kind: str) -> list[str]:
        base = REPO / "core" / kind
        return sorted(f"core/{kind}/{p.name}" for p in base.iterdir() if p.is_dir() and not p.name.startswith("."))

    for comp in dirs("engines") + dirs("agents"):
        expected = COVERAGE.get(comp)
        if expected is None:
            errors.append(f"{comp} has no registry coverage declared — add it to COVERAGE in scripts/validate_registry.py and register a capability")
            continue
        if not any(cid in id_set for cid in expected):
            errors.append(f"{comp} is declared to map to {expected} but none are registered in capabilities.json")

    # 3) Stale coverage entries (a mapped id no longer exists).
    for comp, cids in COVERAGE.items():
        for cid in cids:
            if cid not in id_set:
                errors.append(f"COVERAGE[{comp}] references unknown capability id '{cid}'")

    # 4) PF-12 drift-guard: no two live engine/agent capabilities may share a port
    #    (an active routing break the composer/MCP cannot resolve).
    for port_, cids in reg.port_collisions().items():
        errors.append(f"port collision on :{port_} — {cids} (each service needs a unique port)")

    n_eng = len(reg.by_type("engine")) if hasattr(reg, "by_type") else sum(1 for c in reg.all() if c.type.value == "engine")
    print(f"registry: {len(ids)} capabilities, {n_eng} typed 'engine'")
    if errors:
        print(f"\nFAILED — {len(errors)} problem(s):")
        for e in errors:
            print(f"  - {e}")
        return 1
    print("OK — manifest valid and every engine/agent directory is registered.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
