"""Build-time generation of per-engine and per-agent detail pages for the docs site.

Each page is sourced from the capability registry (the single source of truth), so the depth
pages stay honest and in-sync automatically. Run by mkdocs-gen-files during `mkdocs build`.
Index pages go in the nav; the per-capability detail pages are linked from them.
"""
import json
import pathlib

import mkdocs_gen_files

ROOT = pathlib.Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"

data = json.loads(REGISTRY.read_text())
caps = data["capabilities"]
by_id = {c["id"]: c for c in caps}


def effective(c: dict) -> str:
    return c.get("maturity") or c.get("status", "unknown")


def badge(label: str) -> str:
    known = {"live", "planned", "verified", "preclinical", "research-use", "roadmap", "gated"}
    cls = label if label in known else "unknown"
    return f'<span class="cap-badge cap-{cls}">{label}</span>'


def ports_table(ports: list) -> str:
    if not ports:
        return "_None._"
    rows = ["| Name | Shape | Semantic | Notes |", "|---|---|---|---|"]
    for p in ports:
        rows.append(
            f'| `{p["name"]}` | {p.get("shape", "")} | {p.get("semantic", "") or "—"} '
            f'| {p.get("description", "") or ""} |'
        )
    return "\n".join(rows)


def detail_page(c: dict, kind: str) -> str:
    ep = c.get("endpoint") or "—"
    lines = [
        f"# {c['name']}",
        "",
        f'{badge(effective(c))} &nbsp; **Domain:** {c.get("domain", "—")} &nbsp; '
        f'**Type:** {c["type"]}',
        "",
        c.get("description", ""),
        "",
        "## Interface",
        "",
        f"- **Endpoint:** `{ep}`  ·  **Invoke path:** `{c.get('invoke_path', '/')}`",
        f"- **Serving:** {c.get('serving', '—')}  ·  **GPU:** {'yes' if c.get('gpu') else 'no'}"
        f"  ·  **Cost class:** {c.get('cost_class', '—')}",
        "",
        "**Inputs**",
        "",
        ports_table(c.get("inputs", [])),
        "",
        "**Outputs**",
        "",
        ports_table(c.get("outputs", [])),
        "",
    ]
    if c.get("tags"):
        lines += ["**Tags:** " + " · ".join(f"`{t}`" for t in c["tags"]), ""]
    if kind == "Intelligence Agents":
        lines += [
            '!!! note "Runtime dependency"',
            "    This agent is a Retrieval-Augmented Generation service: at runtime it needs a "
            "populated vector database and an LLM API key. When those are absent it returns an "
            "honest degraded response (e.g. HTTP 503) and never fabricates clinical content.",
            "",
        ]
    lines += [
        "---",
        "",
        f"↩ Back to the [{kind} index](index.md) · the "
        "[Capability Maturity Matrix](../../honesty/maturity-matrix.md) · the "
        "[Capability Brief](../../brief/README.md).",
        "",
        "!!! note",
        "    Status and interface are generated from the capability registry "
        "(`lib/hcls_common/capabilities.json`) — the site cannot claim ahead of the code. "
        "All clinical output is decision support for a qualified clinician, never autonomous diagnosis.",
        "",
    ]
    return "\n".join(lines)


def index_page(title: str, rows: list, folder: str, grid_img: str, blurb: str) -> str:
    lines = [
        f"# {title}",
        "",
        blurb,
        "",
        f"![{title}](../../assets/infographics/{grid_img})",
        "",
        "| Capability | Domain | Status |",
        "|---|---|---|",
    ]
    for c in sorted(rows, key=lambda c: c["name"]):
        lines.append(
            f'| [{c["name"]}]({c["id"]}.md) | {c.get("domain", "")} | {badge(effective(c))} |'
        )
    lines += [
        "",
        "Every row is generated from the capability registry; open a capability for its interface, "
        "ports, and honest status. See the full [Capability Maturity Matrix]"
        "(../../honesty/maturity-matrix.md).",
        "",
    ]
    return "\n".join(lines)


engines = [c for c in caps if c["type"] == "engine"]
agents = [c for c in caps if c["type"] == "agent"]

with mkdocs_gen_files.open("factory/engines/index.md", "w") as f:
    f.write(index_page(
        "Engines", engines, "engines", "eight-engines.png",
        "The horizontal compute muscle of the factory. Each engine is a registered capability with a "
        "typed interface and an honest status.",
    ))
for c in engines:
    with mkdocs_gen_files.open(f"factory/engines/{c['id']}.md", "w") as f:
        f.write(detail_page(c, "Engines"))
    mkdocs_gen_files.set_edit_path(f"factory/engines/{c['id']}.md", "lib/hcls_common/capabilities.json")

with mkdocs_gen_files.open("factory/agents/index.md", "w") as f:
    f.write(index_page(
        "Intelligence Agents", agents, "agents", "eight-agents.png",
        "The clinical reasoning layers over the engines. Each agent is **decision support for a "
        "qualified clinician — never autonomous diagnosis or prescribing.**",
    ))
for c in agents:
    with mkdocs_gen_files.open(f"factory/agents/{c['id']}.md", "w") as f:
        f.write(detail_page(c, "Intelligence Agents"))
    mkdocs_gen_files.set_edit_path(f"factory/agents/{c['id']}.md", "lib/hcls_common/capabilities.json")
