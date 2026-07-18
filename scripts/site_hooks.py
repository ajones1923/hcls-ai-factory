"""Native mkdocs build hook — injects registry-sourced values into templates.

Keeps the designed home hero *honest by construction* (INV-1): the maturity board it renders is
computed from the capability registry at build time, never hand-typed in the template.
Wired via `hooks:` in mkdocs.yml (no extra plugin).
"""
import json
import pathlib
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"

# A curated, recognizable spread for the hero board (id -> short display label). The BADGE for each
# is looked up live from the registry (INV-1 — never hand-typed); only the label + selection are
# curated, and the selection is an honest range (verified / live / gated / planned), full matrix a
# click away. If an id ever leaves the registry it is simply skipped.
_HERO_PICK = [
    ("variant-store", "Variant Store"),
    ("singlecell-compute", "Single-Cell"),
    ("esmfold-model", "ESMFold"),
    ("esm2-search", "ESM-2 Search"),
    ("genomics-engine", "Genomic Foundation"),
    ("precision-intelligence-engine", "Precision Intelligence"),
    ("chai2-binder-design", "Chai-2 Binder Design"),
    ("chai1-structure", "Chai-1 Co-fold"),
]


def _effective(c: dict) -> str:
    return c.get("maturity") or c.get("status", "unknown")


def _board_html() -> str:
    caps = json.loads(REGISTRY.read_text())["capabilities"]
    by_id = {c["id"]: c for c in caps}
    counts = Counter(_effective(c) for c in caps)
    cells = "".join(
        f'<li class="hcls-mm__cell"><span class="hcls-mm__name">{label}</span>'
        f'<span class="cap-badge cap-{_effective(by_id[cid])}">{_effective(by_id[cid])}</span></li>'
        for cid, label in _HERO_PICK if cid in by_id
    )
    summary = f'{counts.get("verified", 0)} verified · {counts.get("live", 0)} live'
    return (
        '<aside class="hcls-mm" aria-label="Capability maturity — generated from the registry">'
        f'<div class="hcls-mm__head"><span>Capability Maturity</span>'
        f'<span class="hcls-mm__sum">{summary}</span></div>'
        f'<ul class="hcls-mm__grid">{cells}</ul>'
        '<div class="hcls-mm__cap">the site goes green as the factory matures</div></aside>'
    )


def on_config(config, **kwargs):
    config["extra"]["maturity_board"] = _board_html()
    return config
