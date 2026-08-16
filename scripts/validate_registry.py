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
    "core/engines/genomic-foundation": ["genomics-engine", "variant-store", "acmg-secondary-findings", "gwas-association", "mosaicism-vaf"],
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

    # 4a) PORT CONVENTION drift-guard (adopted 2026-08-15).
    #
    #     The registry endpoint is the UI port; the API is UI + 1.
    #
    # Before this guard existed, the registry and health-monitor.sh (the cron supervisor that is
    # the ACTUAL deployment mechanism) disagreed on 8 of 13 subjects, and precision-biomarker and
    # neurology each claimed BOTH 8528 and 8529 -- they could not run together. reg.port_collisions()
    # could not see it, because it only compares registry entries with each other; the supervisor's
    # ports are invisible to the registry. So parse the supervisor and cross-check.
    hm = REPO / "health-monitor.sh"
    if hm.exists():
        import re as _re
        supervised: dict[str, int] = {}
        for line in hm.read_text().splitlines():
            m = _re.match(r'\s*"([a-z0-9-]+)\|(\d+)\|', line)
            if m:
                supervised[m.group(1)] = int(m.group(2))

        # a port may be claimed once and once only, across BOTH sources
        claimed: dict[int, list[str]] = {}
        SINGLE = {"genomics-engine", "precision-intelligence-engine",
                  "therapeutic-discovery-engine", "singlecell-compute"}
        for c in reg.all():
            if c.type.value not in ("engine", "agent") or not c.endpoint:
                continue
            try:
                ui_port = int(str(c.endpoint).rsplit(":", 1)[-1])
            except ValueError:
                continue
            claimed.setdefault(ui_port, []).append(f"{c.id}(ui)")
            if c.id not in SINGLE:
                claimed.setdefault(ui_port + 1, []).append(f"{c.id}(api)")
        for port_, who in sorted(claimed.items()):
            if len(who) > 1:
                errors.append(
                    f"port convention violated on :{port_} — claimed by {who}. "
                    f"The registry endpoint is the UI and the API is UI+1, so two capabilities "
                    f"whose UI ports are adjacent will always collide; re-seat one of them."
                )
        # every supervised port must be a port the convention actually allocates
        allocated = set(claimed)
        for name, port_ in sorted(supervised.items()):
            if port_ not in allocated and port_ not in (3000, 8080, 8501, 8510, 9099, 9100, 9400, 19530):
                errors.append(
                    f"health-monitor.sh supervises '{name}' on :{port_}, which the registry does "
                    f"not allocate under UI/UI+1 — the supervisor has drifted from the registry"
                )

    # 4b) taxonomy drift-guard: tags must not contradict type (engine tagged 'agent' or vice versa)
    for cid in reg.type_tag_conflicts():
        errors.append(f"{cid}: tags contradict type (an engine tagged 'agent' or an agent tagged 'engine')")

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
