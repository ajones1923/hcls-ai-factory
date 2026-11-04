#!/usr/bin/env python3
"""How much of each workflow's declared retrieval budget points at an empty collection.

`check_corpus.py` answers "is there a corpus". This answers the question that decides whether an
answer can be trusted: **of the places this workflow was built to look, how many hold anything?**

The two come apart badly. clinical-trial holds 328 vectors across 14 collections, which reads as
thin but working. Ten of those collections are empty, and they are not the cheap ones: patient
matching puts 0.25 of its retrieval weight on `trial_eligibility` and 0.12 on `trial_sites`, both
empty, and **73% of that workflow's budget lands on nothing**. The query succeeds, the answer
arrives, and `collections_searched` reports 14.

Measured 2026-11-03: **56 of 70 workflows below 75% coverage.** Worst per subject:

    cardiology               valvular_disease           2%    12 of 13 collections NEVER CREATED
    clinical-trial           recruitment_optimization  24%
    rare-disease-diagnostic  prenatal_diagnosis        37%
    cart                     car_design                38%
    single-cell              drug_response             46%
    neurology                neuromuscular_eval        70%
    pharmacogenomics         (all six)                100%

Weight, not collection count, because ten trivial collections going quiet matters less than one
heavy one doing it, and a count cannot tell those apart.

Usage:
    .venv/bin/python scripts/check_corpus_weights.py
    .venv/bin/python scripts/check_corpus_weights.py --enforce     # non-zero below the floor
    .venv/bin/python scripts/check_corpus_weights.py --min 0.5

Three subjects (precision-autoimmune, precision-biomarker, clinical-imaging) declare no weight
table — they search a flat collection set — so they have no coverage to report and are listed as
such rather than silently omitted.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

# The platform venv when it is here, otherwise whatever is running this -- a git worktree has no
# .venv of its own, and this script is most useful on a branch before it is deployed.
_VENV_PY = ROOT / ".venv" / "bin" / "python"
PY = str(_VENV_PY) if _VENV_PY.exists() else sys.executable

# subject -> (dir under core/, module, attribute). The attribute is whichever table the engine
# actually READS at retrieval time, which is not the one that reads like the declared budget.
#
# Five subjects carry two tables each: `WORKFLOW_COLLECTION_WEIGHTS` in vector_collections.py,
# which looks authoritative (14 collections summing to 1.00, a `get_search_weights()` accessor,
# tests asserting the sum), and `WORKFLOW_COLLECTION_BOOST` in agent.py, which is what
# `_get_boosted_weights` multiplies against the base weights and ships to Milvus. In
# clinical-trial, neurology, rare-disease-diagnostic, single-cell and cardiology alike,
# `get_search_weights()` has no caller outside tests. This map follows the live path; the other
# table is recorded in docs/build/EMPTY_COLLECTIONS.md as the duplication it is.
#
# cart and pharmacogenomics declare theirs as {workflow: [collection, ...]} in rag_engine.py --
# a boost list with no weights. pharmacogenomics also has a weighted table of the same name in
# agent.py which its retrieval does not read.
WEIGHT_TABLES = {
    "clinical-trial": ("agents/clinical-trial", "src.agent", "WORKFLOW_COLLECTION_BOOST"),
    "cart": ("agents/cart", "src.rag_engine", "STAGE_COLLECTION_BOOST"),
    "pharmacogenomics": ("agents/pharmacogenomics", "src.rag_engine",
                         "WORKFLOW_COLLECTION_BOOST"),
    "neurology": ("agents/neurology", "src.agent", "WORKFLOW_COLLECTION_BOOST"),
    "rare-disease-diagnostic": ("agents/rare-disease-diagnostic", "src.agent",
                                "WORKFLOW_COLLECTION_BOOST"),
    "single-cell": ("agents/single-cell", "src.agent", "WORKFLOW_COLLECTION_BOOST"),
    "cardiology": ("engines/cardiology", "src.agent", "WORKFLOW_COLLECTION_BOOST"),
}

NO_WEIGHT_TABLE = ("precision-autoimmune", "precision-biomarker", "clinical-imaging")

# A workflow below this share of its declared budget is retrieving from a corpus too partial for
# its own confidence score to mean much. Same judgment call as corpus_coverage.DEGRADED_BELOW.
MIN_COVERAGE = 0.75

# The effective weight map, computed the way the engine computes it -- a boost table is NOT a
# weight map. `_get_boosted_weights` multiplies each collection's base weight from
# COLLECTION_CONFIG by its boost (1.0 when the table does not name it, so every collection stays
# in play) and renormalizes. Reading the boost table as weights instead drops the unboosted
# collections entirely and overstates coverage; clinical-trial's patient_matching reads 36% that
# way and 29% correctly.
#
# cart and pharmacogenomics declare {workflow: [collection, ...]} and their engines use it to
# pick which collections to search at all, so for those the list IS the collection set, at base
# weight. `test_retrieval_coverage.py` pins this arithmetic against the engine's own method.
_EXTRACT = """
import json
import {mod} as m
import src.rag_engine as R

t = getattr(m, {attr!r})
base = {{str(n): float(cfg.get("weight", 0.05))
        for n, cfg in R.COLLECTION_CONFIG.items()}}


def effective(v):
    if not isinstance(v, dict):
        return {{str(c): base.get(str(c), 0.05) for c in v}}
    boosted = {{n: w * float(v.get(n, 1.0)) for n, w in base.items()}}
    total = sum(boosted.values())
    return {{n: w / total for n, w in boosted.items()}} if total > 0 else boosted


print("@@" + json.dumps({{str(k).split(".")[-1].lower(): effective(v) for k, v in t.items()}}))
"""


def weight_maps(subject: str) -> tuple[dict[str, dict[str, float]], str]:
    """Import one subject's weight table out of process.

    Out of process because these modules import pymilvus and an application Settings model at
    module scope, and two subjects shadow each other's `src` package. Same cwd + PYTHONPATH
    convention as scripts/run_all_tests.py.
    """
    rel, mod, attr = WEIGHT_TABLES[subject]
    d = ROOT / "core" / rel
    path = os.pathsep.join([str(d), str(ROOT / "lib"), str(d / "src")])
    try:
        r = subprocess.run(
            [PY, "-c", _EXTRACT.format(mod=mod, attr=attr)],
            cwd=d, capture_output=True, text=True, timeout=300,
            env=dict(os.environ, PYTHONPATH=path),
        )
    except subprocess.TimeoutExpired:
        return {}, "timed out importing the weight table"
    for line in (r.stdout or "").splitlines():
        if line.startswith("@@"):
            return json.loads(line[2:]), ""
    err = (r.stderr or "").strip().splitlines()
    return {}, (err[-1][:90] if err else "no table emitted")


def live_counts() -> tuple[dict[str, int], str]:
    """Row counts for every collection in the live instance."""
    try:
        from pymilvus import Collection, connections, utility
    except ImportError:
        return {}, "pymilvus not installed"

    host = os.getenv("MILVUS_HOST", "localhost")
    port = os.getenv("MILVUS_PORT", "19530")
    try:
        connections.connect(alias="weightcheck", host=host, port=port)
    except Exception as exc:
        return {}, f"cannot reach Milvus at {host}:{port} -- {exc}"

    counts: dict[str, int] = {}
    for name in utility.list_collections(using="weightcheck"):
        try:
            col = Collection(name, using="weightcheck")
            col.load()
            # count(*), not num_entities: the latter counts soft-deleted rows until compaction.
            counts[name] = int(col.query(expr="", output_fields=["count(*)"])[0]["count(*)"])
        except Exception:
            counts[name] = 0
    return counts, ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--enforce", action="store_true",
                    help="exit non-zero if any workflow is below the floor")
    ap.add_argument("--min", type=float, default=MIN_COVERAGE, dest="floor")
    ap.add_argument("--json", metavar="PATH", help="write the full result as JSON")
    a = ap.parse_args()

    from hcls_common.corpus_coverage import coverage

    counts, err = live_counts()
    if err:
        print(f"  {err}")
        return 1 if a.enforce else 0

    rows = []
    notes = []
    for subject in sorted(WEIGHT_TABLES):
        maps, err = weight_maps(subject)
        if err:
            notes.append(f"{subject}: {err}")
            continue
        for workflow, weights in maps.items():
            v = coverage(weights, {c: counts.get(c, 0) for c in weights})
            rows.append((subject, workflow, v))

    rows.sort(key=lambda r: r[2].coverage)

    print(f"  {'subject':26s}{'workflow':26s}{'covered':>9s}{'sources':>10s}  heaviest empty")
    for subject, workflow, v in rows:
        heaviest = v.silent[0] if v.silent else "-"
        # An absent collection and an empty one both retrieve nothing, but one needs
        # setup_collections.py and the other needs a seeder. Say which.
        if heaviest != "-" and heaviest not in counts:
            heaviest += "  (never created)"
        print(f"  {subject:26s}{workflow[:25]:26s}{v.coverage:8.0%}"
              f"{len(v.contributing):>5d}/{len(v.declared):<4d}  {heaviest}")

    # A collection that was never created is invisible to check_corpus.py, which can only
    # iterate what Milvus holds. cardiology was exactly this on 2026-11-03: twelve workflows,
    # twelve cardio_* collections, not one of them created, and no row anywhere to show it. The
    # only reason its coverage is not 0% is `genomic_evidence`, the shared collection it borrows.
    for subject in sorted(WEIGHT_TABLES):
        declared = {c for s, _, v in rows if s == subject for c in v.declared}
        absent = sorted(declared - counts.keys())
        if absent:
            print(f"\n  !! {subject}: {len(absent)} of {len(declared)} declared collections are "
                  f"ABSENT from Milvus, not merely empty -- they were never created.\n"
                  f"     Run its scripts/setup_collections.py, then seed. Missing: "
                  f"{', '.join(absent[:6])}" + (" ..." if len(absent) > 6 else ""))

    for n in notes:
        print(f"  ! {n}")
    for subject in NO_WEIGHT_TABLE:
        print(f"  - {subject}: no weight table declared; searches a flat collection set")

    below = [r for r in rows if r[2].coverage < a.floor]
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(
            [{"subject": s, "workflow": w, **v.as_dict()} for s, w, v in rows], indent=1))

    if below:
        worst = below[0]
        print(f"\n  {len(below)} of {len(rows)} workflows below {a.floor:.0%} coverage."
              f" Worst: {worst[0]} / {worst[1]} at {worst[2].coverage:.0%}.")
        print("  These answers are drawn from a fraction of the evidence the workflow declares,")
        print("  and nothing in the retrieval path fails. Seed the collections --")
        print("  docs/build/CORPUS_SEEDING.md, and docs/build/EMPTY_COLLECTIONS.md for what can")
        print("  fill each one.")
        return 1 if a.enforce else 0

    print("\n  OK -- every workflow retrieves from the sources it declares.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
