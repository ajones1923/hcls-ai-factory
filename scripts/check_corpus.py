#!/usr/bin/env python3
"""Report what each agent's corpus actually contains.

An agent with an empty collection answers HTTP 200, passes its whole test suite, and retrieves
nothing from it. The service is up, the collection is loaded, the query succeeds, and the result
set is empty — there is no error anywhere in that chain. The model then answers from its own
knowledge and the reply looks exactly like a sourced one.

Measured on 2026-09-16: **44 of 113 collections were empty**, including 11 of 13 for
clinical-trial, 10 of 13 for rare-disease and 8 of 13 for precision-autoimmune. Every one of those
agents was passing its clinical eval at the time, because the eval graded the ANSWER and never
asked whether retrieval contributed. `scripts/run_clinical_eval.py` now reports an UNGROUNDED
verdict for that; this script shows the corpus side of the same question.

Usage:
    .venv/bin/python scripts/check_corpus.py              # report
    .venv/bin/python scripts/check_corpus.py --enforce    # non-zero if any agent is below floor
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

# An agent below this many vectors is retrieving from almost nothing, whatever it returns.
MIN_VECTORS_PER_AGENT = 500

# Milvus collection prefix -> the subject that searches it.
PREFIXES = {
    "cart": "cart", "biomarker": "precision-biomarker", "pgx": "pharmacogenomics",
    "autoimmune": "precision-autoimmune", "neuro": "neurology", "trial": "clinical-trial",
    "rd": "rare-disease-diagnostic", "sc": "single-cell", "imaging": "clinical-imaging",
    "genomic": "shared (genomic evidence)", "tsc": "tuberous-sclerosis",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--enforce", action="store_true")
    ap.add_argument("--min-vectors", type=int, default=MIN_VECTORS_PER_AGENT)
    a = ap.parse_args()

    try:
        from pymilvus import Collection, connections, utility
    except ImportError:
        print("  pymilvus not installed"); return 0

    host = os.getenv("MILVUS_HOST", "localhost")
    port = os.getenv("MILVUS_PORT", "19530")
    try:
        connections.connect(alias="corpuscheck", host=host, port=port)
    except Exception as exc:
        print(f"  cannot reach Milvus at {host}:{port} — {exc}")
        return 1 if a.enforce else 0

    stats = collections.defaultdict(lambda: {"cols": 0, "empty": 0, "vectors": 0, "names": []})
    total_empty = []
    for name in utility.list_collections(using="corpuscheck"):
        try:
            col = Collection(name, using="corpuscheck")
            # `num_entities` counts soft-deleted rows until compaction, so a collection that has
            # been re-ingested reports MORE than it holds: tsc_literature read 18 while holding
            # 6 live rows. count(*) is the number an operator actually wants.
            try:
                col.load()
                n = int(col.query(expr="", output_fields=["count(*)"])[0]["count(*)"])
            except Exception:
                n = col.num_entities
        except Exception:
            n = 0
        subject = PREFIXES.get(name.split("_")[0], name.split("_")[0])
        s = stats[subject]
        s["cols"] += 1
        s["vectors"] += n
        if n == 0:
            s["empty"] += 1
            s["names"].append(name)
            total_empty.append(name)

    print(f"  {'subject':28s}{'collections':>12s}{'empty':>7s}{'vectors':>10s}")
    thin = []
    for subject, s in sorted(stats.items(), key=lambda kv: kv[1]["vectors"]):
        mark = ""
        if s["vectors"] < a.min_vectors:
            thin.append((subject, s["vectors"]))
            mark = f"   << under {a.min_vectors}"
        print(f"  {subject:28s}{s['cols']:12d}{s['empty']:7d}{s['vectors']:10,d}{mark}")

    ncols = sum(s["cols"] for s in stats.values())
    print(f"\n  {len(total_empty)} of {ncols} collections are empty")
    if total_empty:
        print("  empty: " + ", ".join(sorted(total_empty)[:8])
              + (f" … +{len(total_empty)-8} more" if len(total_empty) > 8 else ""))
    if thin:
        print(f"\n  {len(thin)} subject(s) below {a.min_vectors} vectors: "
              + ", ".join(f"{n} ({v:,})" for n, v in thin))
        print("  An empty collection is not an error anywhere in the stack: the query succeeds and\n"
              "  returns nothing, and the model answers from memory. Seed it — docs/build/CORPUS_SEEDING.md")
        return 1 if a.enforce else 0
    print("\n  OK — every subject has a corpus above the floor.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
