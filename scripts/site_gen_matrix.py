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


def badge(status: str) -> str:
    known = {"live", "planned", "verified", "preclinical", "research-use", "roadmap", "gated"}
    cls = status if status in known else "unknown"
    return f'<span class="cap-badge cap-{cls}">{status}</span>'


def table(rows) -> str:
    out = [
        '<table class="cap-matrix" markdown="0">',
        "<thead><tr><th>Capability</th><th>Domain</th><th>Status</th></tr></thead>",
        "<tbody>",
    ]
    for c in rows:
        out.append(
            f'<tr><td><strong>{c["name"]}</strong></td>'
            f'<td>{c.get("domain", "")}</td>'
            f'<td>{badge(c.get("status", "unknown"))}</td></tr>'
        )
    out.append("</tbody></table>")
    return "\n".join(out)


engines = sorted((c for c in caps if c["type"] == "engine"), key=lambda c: c["name"])
agents = sorted((c for c in caps if c["type"] == "agent"), key=lambda c: c["name"])
counts = Counter(c.get("status") for c in caps)

lines = [
    "# Capability Maturity Matrix",
    "",
    "> Generated at build time from the capability registry "
    f"(`lib/hcls_common/capabilities.json`, v{version}). "
    "The site cannot show a status the registry does not declare — **honesty by construction.**",
    "",
    f"**{counts.get('live', 0)} live · {counts.get('planned', 0)} planned** "
    f"across {len(caps)} registered capabilities.",
    "",
    "## Engines",
    "",
    table(engines),
    "",
    "## Intelligence Agents",
    "",
    table(agents),
    "",
    '!!! note "On the status vocabulary"',
    "    The registry currently declares two states — **live** (served by a real model against real",
    "    input; a `live` capability is never mock-served) and **planned**. A richer maturity",
    "    vocabulary (*verified · preclinical · research-use · roadmap · gated*) is a tracked",
    "    enhancement — see the [Roadmap](../roadmap.md). When the registry schema carries it, this",
    "    matrix upgrades automatically; the badge styling already supports all seven states.",
    "",
]

with mkdocs_gen_files.open("honesty/maturity-matrix.md", "w") as f:
    f.write("\n".join(lines))

# "Edit this page" points at the registry — the real source of these statuses.
mkdocs_gen_files.set_edit_path("honesty/maturity-matrix.md", "lib/hcls_common/capabilities.json")
