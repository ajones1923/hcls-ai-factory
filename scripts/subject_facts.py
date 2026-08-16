#!/usr/bin/env python3
"""Extract real, per-subject facts from the repo.

This exists so the 85 documentation assets are built from what the code actually is — module names,
ports, test counts, registered capabilities, dependencies — rather than from a template with the
subject's name substituted in. Every number a generated guide states comes from here.

Usage:
    .venv/bin/python scripts/subject_facts.py            # table
    .venv/bin/python scripts/subject_facts.py --json out.json
    .venv/bin/python scripts/subject_facts.py cardiology # one subject, full detail
"""
from __future__ import annotations
import ast, json, pathlib, re, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SKIP = {".venv", "venv", "site-packages", "node_modules", "__pycache__", ".git",
        ".pytest_cache", "vendor_rfdiffusion", "vendor_proteinmpnn"}
keep = lambda p: not any(x in p.parts for x in SKIP)

ROSTER = [
    ("engine", "genomic-foundation", "Genomic Foundation"),
    ("engine", "precision-intelligence", "Precision Intelligence"),
    ("engine", "therapeutic-discovery", "Therapeutic Discovery"),
    ("engine", "clinical-imaging", "Clinical Imaging"),
    ("engine", "precision-oncology", "Precision Oncology"),
    ("engine", "cardiology", "Cardiology"),
    ("engine", "structural-biology", "Structural Biology"),
    ("engine", "single-cell", "Single-Cell Analysis"),
    ("agent", "cart", "CAR-T Intelligence"),
    ("agent", "precision-biomarker", "Precision Biomarker"),
    ("agent", "pharmacogenomics", "Pharmacogenomics"),
    ("agent", "precision-autoimmune", "Precision Autoimmune"),
    ("agent", "neurology", "Neurology Intelligence"),
    ("agent", "clinical-trial", "Clinical Trial Intelligence"),
    ("agent", "rare-disease-diagnostic", "Rare Disease Diagnostic"),
    ("agent", "single-cell", "Single-Cell Intelligence"),
    ("program", "tuberous-sclerosis", "Tuberous Sclerosis Complex"),
]
BASE = {"engine": "core/engines", "agent": "core/agents", "program": "core/disease-programs"}

# demo keys from docs/demos/DEMO_CATALOG.md (E/A/P axis, NOT the D1-D7 portfolio)
DEMO_KEY = {("engine", i): f"E{i+1}" for i in range(8)}
DEMO_KEY.update({("agent", i): f"A{i+1}" for i in range(8)})
DEMO_KEY[("program", 0)] = "P1"


def registry():
    m = json.loads((ROOT / "lib/hcls_common/capabilities.json").read_text())
    return m["capabilities"]


def coverage_map():
    """component path -> capability ids, from validate_registry.COVERAGE."""
    src = (ROOT / "scripts/validate_registry.py").read_text()
    m = re.search(r"COVERAGE\s*[:=][^{]*\{(.*?)\n\}", src, re.S)
    out = {}
    if m:
        for line in m.group(1).splitlines():
            km = re.match(r'\s*"([^"]+)":\s*\[(.*?)\]', line)
            if km:
                out[km.group(1)] = re.findall(r'"([^"]+)"', km.group(2))
    return out


def public_symbols(py: pathlib.Path, limit=14):
    """Top-level classes/functions — what a reader would actually call."""
    try:
        tree = ast.parse(py.read_text(errors="ignore"))
    except Exception:
        return []
    out = []
    for n in tree.body:
        if isinstance(n, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if n.name.startswith("_"):
                continue
            doc = (ast.get_docstring(n) or "").strip().split("\n")[0]
            out.append({"kind": "class" if isinstance(n, ast.ClassDef) else "func",
                        "name": n.name, "doc": doc[:150]})
    return out[:limit]


def facts_for(kind, name, title, idx):
    d = ROOT / BASE[kind] / name
    root = d / "agent" if (d / "agent").is_dir() else d
    py = [f for f in root.rglob("*.py") if keep(f)]
    tests = [f for f in py if f.name.startswith("test_") or "tests" in f.parts]
    loc = 0
    for f in py:
        try:
            loc += sum(1 for _ in f.open(errors="ignore"))
        except Exception:
            pass

    caps = [c for c in registry()
            if c["id"] in coverage_map().get(f"{BASE[kind]}/{name}", [])]
    # Single-service portals have NO UI/API split, so "API = UI + 1" does not apply to them --
    # asserting it would publish genomic-foundation's API as 5001, which is precision-intelligence's
    # UI. See docs/build/PORT_MAP.md.
    SINGLE = {"genomics-engine", "precision-intelligence-engine",
              "therapeutic-discovery-engine", "singlecell-compute"}
    ui, single = None, False
    for c in caps:
        if c.get("type") in ("engine", "agent") and c.get("endpoint"):
            ui = int(str(c["endpoint"]).rsplit(":", 1)[-1])
            single = c["id"] in SINGLE
            break

    # the modules a reader should start from: biggest non-test sources
    mods = sorted((f for f in py if f not in tests),
                  key=lambda f: -f.stat().st_size)[:6]
    modules = [{"path": str(f.relative_to(root)), "symbols": public_symbols(f)} for f in mods]

    reqs = []
    for r in root.rglob("requirements*.txt"):
        if keep(r):
            reqs += [l.split("#")[0].strip() for l in r.read_text().splitlines()
                     if l.strip() and not l.startswith("#")]

    readme = root / "README.md"
    summary = ""
    if readme.exists():
        for line in readme.read_text(errors="ignore").splitlines():
            s = line.strip()
            if s and not s.startswith("#") and not s.startswith("!["):
                summary = s
                break

    return {
        "kind": kind, "slug": name, "title": title,
        "demo_key": DEMO_KEY[(kind, idx)],
        "path": str(root.relative_to(ROOT)),
        "py_files": len(py), "loc": loc, "test_files": len(tests),
        "capabilities": [{"id": c["id"], "type": c.get("type"), "status": c.get("status"),
                          "endpoint": c.get("endpoint"),
                          "description": (c.get("description") or "")[:400]} for c in caps],
        "ui_port": ui, "api_port": (None if (single or not ui) else ui + 1),
        "single_service": single,
        "modules": modules,
        "dependencies": sorted(set(reqs))[:20],
        "summary": summary,
        "has_docker": any(keep(f) for f in root.rglob("Dockerfile")),
    }


def build():
    out, ei, ai = [], 0, 0
    for kind, name, title in ROSTER:
        if kind == "engine":
            idx, ei = ei, ei + 1
        elif kind == "agent":
            idx, ai = ai, ai + 1
        else:
            idx = 0
        out.append(facts_for(kind, name, title, idx))
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    data = build()
    if "--json" in sys.argv:
        p = sys.argv[sys.argv.index("--json") + 1]
        json.dump(data, open(p, "w"), indent=1)
        print(f"wrote {p} ({len(data)} subjects)")
        return 0
    if args:
        for f in data:
            if args[0] in f["slug"]:
                print(json.dumps(f, indent=1))
        return 0
    print(f"{'subject':26s}{'demo':6s}{'py':>5s}{'LOC':>8s}{'tests':>6s}{'UI':>6s}{'API':>5s}  caps")
    for f in data:
        print(f"  {f['slug']:24s}{f['demo_key']:6s}{f['py_files']:5d}{f['loc']:8d}"
              f"{f['test_files']:6d}{str(f['ui_port'] or '-'):>6s}{str(f['api_port'] or '-'):>5s}"
              f"  {len(f['capabilities'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
