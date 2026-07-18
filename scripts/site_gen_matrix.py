"""Build-time generation of the Capability Maturity Matrix for the public docs site.

Reads the capability registry (the single source of truth) and emits a Markdown page, so the
site can never assert a status the registry does not declare — honesty by construction (Pillar 16).
Run by mkdocs-gen-files during `mkdocs build`. Not published itself.
"""
import json
import pathlib
from collections import Counter

import mkdocs_gen_files

ROOT = pathlib.Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"

data = json.loads(REGISTRY.read_text())
caps = data["capabilities"]
version = data.get("version", "?")

# Display order for the honest labels (most-proven -> not-built -> access-gated).
ORDER = ["verified", "live", "research-use", "preclinical", "planned", "roadmap", "gated", "unknown"]


def effective(c: dict) -> str:
    """The single most honest badge: the evidence/access tier if set, else the serving status."""
    return c.get("maturity") or c.get("status", "unknown")


def badge(label: str) -> str:
    known = {"live", "planned", "verified", "preclinical", "research-use", "roadmap", "gated"}
    cls = label if label in known else "unknown"
    return f'<span class="cap-badge cap-{cls}">{label}</span>'


def table(rows) -> str:
    out = [
        '<table class="cap-matrix" markdown="0">',
        "<thead><tr><th>Capability</th><th>Domain</th><th>Status</th></tr></thead>",
        "<tbody>",
    ]
    for c in sorted(rows, key=lambda c: c["name"]):
        out.append(
            f'<tr><td><strong>{c["name"]}</strong></td>'
            f'<td>{c.get("domain", "")}</td>'
            f'<td>{badge(effective(c))}</td></tr>'
        )
    out.append("</tbody></table>")
    return "\n".join(out)


GROUPS = [
    ("Engines", lambda c: c["type"] == "engine"),
    ("Intelligence Agents", lambda c: c["type"] == "agent"),
    ("Models & NIMs", lambda c: c["type"] in ("model", "nim")),
    ("Platform Services & Pipeline Stages", lambda c: c["type"] in ("service", "stage")),
]

counts = Counter(effective(c) for c in caps)
summary = " · ".join(f"{counts[s]} {s}" for s in ORDER if counts.get(s))

lines = [
    "# Capability Maturity Matrix",
    "",
    "> Generated at build time from the capability registry "
    f"(`lib/hcls_common/capabilities.json`, v{version}). "
    "The site cannot show a status the registry does not declare — **honesty by construction.**",
    "",
    f"**{summary}** — across {len(caps)} registered capabilities.",
    "",
    "Each badge shows the single most honest label: a capability's **evidence/access tier** "
    "(`verified` · `research-use` · `preclinical` · `roadmap` · `gated`) where the registry "
    "documents one, otherwise its **serving status** (`live` · `planned`).",
    "",
]
for title, pred in GROUPS:
    rows = [c for c in caps if pred(c)]
    if rows:
        lines += [f"## {title}", "", table(rows), ""]

lines += [
    '!!! note "How to read the tiers"',
    "    **`live`** — served by a real model against real input (never mock-served). "
    "**`verified`** — additionally proven against real, recorded input (e.g. on real HG002). "
    "**`planned`** — on the roadmap, not yet running. **`gated`** — real but partnership-/"
    "license-gated. **`preclinical` · `research-use` · `roadmap`** — honest evidence tiers carried "
    "where documented. The tier is **orthogonal to serving** — a capability can be `live` *and* "
    "`research-use`. Tiers are populated conservatively, only where the registry documents a basis; "
    "feature-level caveats (e.g. TSC gene therapy is preclinical) are on the "
    "[Honesty & Governance](../honesty/index.md) ledger.",
    "",
]

with mkdocs_gen_files.open("honesty/maturity-matrix.md", "w") as f:
    f.write("\n".join(lines))

# "Edit this page" points at the registry — the real source of these statuses.
mkdocs_gen_files.set_edit_path("honesty/maturity-matrix.md", "lib/hcls_common/capabilities.json")
