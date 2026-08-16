#!/usr/bin/env python3
"""Run every engine / agent / disease-program test suite and print one table.

Why this exists: each subject needed a DIFFERENT pytest invocation, so there was no single
command that told you the health of the platform. Two traps made that worse:

  1. Several subjects ship `src/collections.py`, which SHADOWS the Python standard library.
     Putting their `src/` on PYTHONPATH kills the interpreter before collection
     ("cannot import name 'namedtuple' from partially initialized module 'collections'").
     So `src/` is added only for subjects that do not shadow a stdlib module.

  2. `core/engines/structural-biology/vendor_rfdiffusion/` is VENDORED third-party code whose
     own tests need `rfdiffusion` + `dgl` (gated, GPU-only). Those are not our tests and are
     excluded; the engine's own suite passes without them.

Usage:
    .venv/bin/python scripts/run_all_tests.py            # all subjects
    .venv/bin/python scripts/run_all_tests.py cardio     # substring filter
    .venv/bin/python scripts/run_all_tests.py --json out.json
"""
from __future__ import annotations
import json, os, pathlib, re, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY = str(ROOT / ".venv" / "bin" / "python")

ROSTER = [
    ("engine", "genomic-foundation"), ("engine", "precision-intelligence"),
    ("engine", "therapeutic-discovery"), ("engine", "clinical-imaging"),
    ("engine", "precision-oncology"), ("engine", "cardiology"),
    ("engine", "structural-biology"), ("engine", "single-cell"),
    ("agent", "cart"), ("agent", "precision-biomarker"), ("agent", "pharmacogenomics"),
    ("agent", "precision-autoimmune"), ("agent", "neurology"), ("agent", "clinical-trial"),
    ("agent", "rare-disease-diagnostic"), ("agent", "single-cell"),
    ("program", "tuberous-sclerosis"),
]
BASE = {"engine": "core/engines", "agent": "core/agents", "program": "core/disease-programs"}
IGNORE = {"structural-biology": ["vendor_rfdiffusion"]}
STDLIB = set(sys.stdlib_module_names)


def shadows_stdlib(src: pathlib.Path) -> list[str]:
    """Module files under src/ whose name collides with the standard library."""
    if not src.is_dir():
        return []
    return sorted(f.stem for f in src.glob("*.py") if f.stem in STDLIB)


def run(kind: str, name: str) -> dict:
    d = ROOT / BASE[kind] / name
    src = d / "src"
    clash = shadows_stdlib(src)
    path = [d, ROOT / "lib"] + ([] if clash else ([src] if src.is_dir() else []))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(str(p) for p in path))
    # Installing scanpy/anndata brings zarr and fast-array-utils, which register pytest plugins
    # that import an application Settings model at startup and kill COLLECTION for three subjects.
    # Disable those two specifically -- do NOT set PYTEST_DISABLE_PLUGIN_AUTOLOAD, because the
    # suites genuinely need pytest-asyncio (38 async tests fail without it).
    cmd = [PY, "-m", "pytest", "-q", "--no-header", "-p", "no:cacheprovider",
           "-p", "no:zarr", "-p", "no:fast_array_utils"]
    for ig in IGNORE.get(name, []):
        cmd += [f"--ignore={ig}"]
    try:
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=1800, env=env)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        out = "TIMEOUT"
    g = lambda p: int(m.group(1)) if (m := re.search(p, out)) else 0
    res = dict(kind=kind, name=name,
               passed=g(r"(\d+) passed"), failed=g(r"(\d+) failed"),
               errors=g(r"(\d+) error"), skipped=g(r"(\d+) skipped"),
               shadowed=clash, note="")
    if out == "TIMEOUT":
        res["note"] = "TIMED OUT"
    elif sum(res[k] for k in ("passed", "failed", "errors", "skipped")) == 0:
        res["note"] = "no tests collected"
    miss = sorted(set(re.findall(r"No module named '([A-Za-z0-9_]+)'", out)))
    if miss:
        res["note"] = (res["note"] + " · " if res["note"] else "") + "missing: " + ",".join(miss[:5])
    return res


def main() -> int:
    argv = sys.argv[1:]
    jsonout = None
    if "--json" in argv:
        i = argv.index("--json")
        jsonout = argv[i + 1] if i + 1 < len(argv) else None
        del argv[i:i + 2]                      # drop BOTH flag and value, or the value
    args = [a for a in argv if not a.startswith("--")]   # would be read as a name filter
    roster = [(k, n) for k, n in ROSTER if not args or any(a in n for a in args)]

    rows = []
    print(f"{'subject':28s}{'kind':9s}{'pass':>6s}{'fail':>6s}{'err':>5s}{'skip':>6s}  note")
    for kind, name in roster:
        r = run(kind, name)
        rows.append(r)
        flag = "  [src/ withheld: shadows " + ",".join(r["shadowed"]) + "]" if r["shadowed"] else ""
        print(f"  {name:26s}{kind:9s}{r['passed']:6d}{r['failed']:6d}{r['errors']:5d}"
              f"{r['skipped']:6d}  {r['note']}{flag}")
    t = lambda k: sum(r[k] for r in rows)
    print(f"\n  {len(rows)} subjects · passed {t('passed')} · failed {t('failed')} "
          f"· errors {t('errors')} · skipped {t('skipped')}")
    if jsonout:
        json.dump(rows, open(jsonout, "w"), indent=1)
        print(f"  wrote {jsonout}")
    return 1 if (t("failed") or t("errors")) else 0


if __name__ == "__main__":
    raise SystemExit(main())
