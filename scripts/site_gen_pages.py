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


def _narrative(c: dict, kind: str) -> str:
    """Inject the hand-authored, expanded explainer for a capability if one has been written.

    Content lives OUTSIDE docs/ (in site_content/{engines,agents}/{id}.md) so mkdocs never serves
    it as a standalone page; the generator splices it into the registry-driven detail page. Returns
    the markdown, or "" if no explainer exists yet — so pages degrade cleanly during the rollout.
    """
    folder = "engines" if kind == "Engines" else "agents"
    src = ROOT / "site_content" / folder / f"{c['id']}.md"
    return src.read_text().rstrip() + "\n\n" if src.exists() else ""


def _video(c: dict) -> str:
    """Embed the capability's narrated 60-second explainer video if one has been produced.

    Uses the hero infographic as the poster frame; degrades to nothing if no video exists yet.
    Video files live in docs/assets/videos/<id>.mp4 (web-optimized, faststart).
    """
    vid = ROOT / "docs" / "assets" / "videos" / f"{c['id']}.mp4"
    if not vid.exists():
        return ""
    poster = ROOT / "docs" / "assets" / "infographics" / "heros" / f"{c['id']}.png"
    poster_attr = f' poster="../../assets/infographics/heros/{c["id"]}.png"' if poster.exists() else ""
    return (
        f'<video class="cap-video" controls preload="none" playsinline{poster_attr}>\n'
        f'  <source src="../../assets/videos/{c["id"]}.mp4" type="video/mp4">\n'
        f'  Your browser can\'t play embedded video — '
        f'<a href="../../assets/videos/{c["id"]}.mp4">download the explainer</a>.\n'
        f'</video>\n'
        '/// caption\nA narrated, captioned explainer. Decision support for a qualified clinician.\n///\n\n'
    )


def _hero(c: dict) -> str:
    """Embed the capability's flow-diagram hero if one has been drawn for it (keyed by id).

    Returns markdown (image + illustrative caption) or an empty string if no hero exists, so the
    detail page degrades cleanly for any capability that hasn't been illustrated yet.
    """
    img = ROOT / "docs" / "assets" / "infographics" / "heros" / f"{c['id']}.png"
    if not img.exists():
        return ""
    return (
        f"![{c['name']} — what it takes in, what it computes, what it returns]"
        f"(../../assets/infographics/heros/{c['id']}.png)\n"
        "/// caption\nIllustrative. Decision support for a qualified clinician — never autonomous "
        "diagnosis or prescribing.\n///\n"
    )


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
        # Narrated 60s explainer video (if produced), then the registry-sourced hero infographic.
        _video(c),
        _hero(c),
        # Hand-authored expanded explainer (site_content/), injected before the generated interface.
        _narrative(c, kind),
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
        "All clinical output is [decision support](../../honesty/decision-support.md) for a qualified "
        "clinician, never autonomous diagnosis.",
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


# ── Capability Brief — Technical Cut (registry-generated; retires the hard-coded HTML, SC-4) ──
def brief_row(c: dict) -> str:
    ep = c.get("endpoint") or ("library" if c.get("serving") == "none" else "—")
    return f'| **{c["name"]}** | {c.get("domain", "")} | `{ep}` | {badge(effective(c))} |'


def brief_group(title: str, rows: list) -> str:
    out = [f"### {title}", "", "| Capability | Domain | Endpoint | Status |", "|---|---|---|---|"]
    out += [brief_row(c) for c in sorted(rows, key=lambda c: c["name"])]
    return "\n".join(out) + "\n"


brief = [
    "# Capability Brief — Technical Cut",
    "",
    "The whole factory in one read: **8 engines · 8 intelligence agents · 1 disease program (TSC)** on "
    "one platform. Every endpoint and status below is **generated at build time from the capability "
    "registry** — no hand-typed drift. For the human *why*, see the [Mission Cut](mission.md).",
    "",
    "![HCLS AI Factory architecture](architecture.svg)",
    "/// caption",
    "Illustrative architecture diagram — teal engines, indigo agents, ember TSC.",
    "///",
    "",
    "## The roster",
    "",
    brief_group("Engines", [c for c in caps if c["type"] == "engine"]),
    brief_group("Intelligence Agents", [c for c in caps if c["type"] == "agent"]),
    brief_group("Models & NIMs", [c for c in caps if c["type"] in ("model", "nim")]),
    brief_group("Platform Services & Pipeline Stages", [c for c in caps if c["type"] in ("service", "stage")]),
    "## Flagship disease program",
    "",
    "[**Tuberous Sclerosis Complex**](../programs/tsc.md) composes the horizontals for one child in one "
    "governed afternoon — TSC1/TSC2 → mTORC1 → everolimus. Gene therapy is **preclinical**; everolimus "
    "is real and approved.",
    "",
    "## Honesty",
    "",
    "All clinical output is [decision support](../honesty/decision-support.md) for a qualified "
    "clinician, never autonomous diagnosis. See the live "
    "[Capability Maturity Matrix](../honesty/maturity-matrix.md) and the "
    "[Honesty Ledger](../honesty/ledger.md). A `live` capability is never mock-served.",
    "",
    '!!! note "Generated from the registry"',
    "    This page's roster, endpoints, and statuses are generated from "
    "`lib/hcls_common/capabilities.json` at build time — the site cannot claim ahead of the code.",
    "",
]
with mkdocs_gen_files.open("brief/README.md", "w") as f:
    f.write("\n".join(brief))
mkdocs_gen_files.set_edit_path("brief/README.md", "lib/hcls_common/capabilities.json")
