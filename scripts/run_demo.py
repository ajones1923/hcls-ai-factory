#!/usr/bin/env python3
"""Run one of the 17 per-subject demonstrations and record a transcript.

Keys are E1-E8 (engines), A1-A8 (agents), P1 (TSC program) — deliberately NOT D1-D7, which is the
established patient-story portfolio in docs/demos/index.md. See DEMO_CATALOG.md.

The rule this script exists to enforce (PRD DR-3):

    A demo declared LIVE whose service is unreachable FAILS. It does not quietly degrade to a
    canned result.

That is the same discipline that caught `singlecell-compute` and `structural-biology-engine`
registered `live` in the capability registry while nothing bound their ports.

Usage:
    .venv/bin/python scripts/run_demo.py --list
    .venv/bin/python scripts/run_demo.py E8
    .venv/bin/python scripts/run_demo.py --check-all      # prerequisites only, runs nothing
"""
from __future__ import annotations
import argparse, importlib.util, json, pathlib, socket, subprocess, sys, urllib.request
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parent.parent
TRANSCRIPTS = ROOT / "demo" / "transcripts"

LIVE, REPRESENTATIVE, BURST = "LIVE", "REPRESENTATIVE", "BURST"


class Demo:
    def __init__(self, key, subject, title, label, *, port=None, packages=(),
                 payload=None, gated=(), runner=None):
        self.key, self.subject, self.title, self.label = key, subject, title, label
        self.port, self.packages, self.payload = port, packages, payload
        self.gated, self.runner = gated, runner

    def missing_packages(self):
        return [p for p in self.packages if importlib.util.find_spec(p) is None]

    def port_open(self):
        if not self.port:
            return None
        with socket.socket() as s:
            s.settimeout(1.5)
            return s.connect_ex(("127.0.0.1", self.port)) == 0

    def check(self):
        """-> (ok, [reasons]). A REPRESENTATIVE demo may run without its service."""
        reasons = []
        miss = self.missing_packages()
        if miss:
            reasons.append(f"missing packages: {', '.join(miss)}")
        if self.payload and not (ROOT / self.payload).exists():
            reasons.append(f"missing payload: {self.payload}")
        if self.label == LIVE and self.port is not None and not self.port_open():
            reasons.append(f"declared LIVE but nothing is listening on :{self.port}")
        if self.gated:
            reasons.append(f"gated (informational): {', '.join(self.gated)}")
        blocking = [r for r in reasons if not r.startswith("gated (informational)")]
        return (not blocking), reasons


DEMOS = [
    Demo("E1", "genomic-foundation", "The variant that was always there", REPRESENTATIVE,
         packages=("duckdb", "statsmodels"), gated=("Parabricks (G2)",)),
    Demo("E2", "precision-intelligence", "Ask the evidence layer a question", LIVE, port=5001),
    Demo("E3", "therapeutic-discovery", "From one protein to a hundred candidates", REPRESENTATIVE,
         gated=("MolMIM (G3)", "DiffDock (G4)")),
    Demo("E4", "clinical-imaging", "The scan that had already answered", LIVE, port=8524,
         payload="demo/requests/imaging_query.json"),
    Demo("E5", "precision-oncology", "The molecular tumour board", LIVE, port=8527,
         payload="demo/requests/oncology_query.json"),
    Demo("E6", "cardiology", "Risk that changes management", LIVE, port=8127,
         payload="demo/requests/cardiology_risk.json"),
    Demo("E7", "structural-biology", "The shape you have to fit", REPRESENTATIVE,
         gated=("CUDA torch (G1)", "ESMFold (G5)")),
    Demo("E8", "single-cell", "Nine populations from one sample", LIVE,
         packages=("scanpy", "anndata"), runner="single_cell"),
    Demo("A1", "cart", "Why this construct, for this patient", LIVE, port=8522,
         payload="demo/requests/cart_query.json"),
    Demo("A2", "precision-biomarker", "The marker that changes the decision", LIVE, port=8529,
         payload="demo/requests/biomarker_query.json"),
    Demo("A3", "pharmacogenomics", "Two patients, same dose, different outcome", LIVE, port=8508,
         payload="demo/requests/pharmacogenomics_query.json"),
    Demo("A4", "precision-autoimmune", "Before the third flare", LIVE, port=8532,
         payload="demo/requests/autoimmune_query.json"),
    Demo("A5", "neurology", "NIHSS, and what comes next", LIVE, port=8536,
         payload="demo/requests/neurology_nihss.json"),
    Demo("A6", "clinical-trial", "The trial that was open all along", LIVE, port=8539,
         payload="demo/requests/trial_match.json"),
    Demo("A7", "rare-disease-diagnostic", "Ending the odyssey", LIVE, port=8545,
         payload="demo/requests/rare_disease_diagnose.json"),
    Demo("A8", "single-cell", "The engine computes, the agent interprets", LIVE, port=8541,
         payload="demo/requests/single_cell_query.json"),
    Demo("P1", "tuberous-sclerosis", "The whole factory, one child", REPRESENTATIVE, port=8561),
]
BY_KEY = {d.key: d for d in DEMOS}


def run_single_cell(log):
    """E8 — the only demo that is real end-to-end today with no gated software and no GPU."""
    import scanpy as sc  # noqa: F401
    src = ROOT / "core/engines/single-cell/src"
    sys.path.insert(0, str(src))
    from single_cell_compute import SingleCellAnalysis, PBMC_MARKERS

    h5ad = ROOT / "core/engines/single-cell/data/pbmc3k_raw.h5ad"
    log(f"dataset       {h5ad.name} ({h5ad.stat().st_size/1e6:.1f} MB)")
    log(f"marker panel  {len(PBMC_MARKERS)} cell types")
    adata = sc.read_h5ad(h5ad)
    log(f"loaded        {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log("running       QC -> normalize -> HVG -> PCA -> neighbors -> Leiden -> marker DE")
    result = SingleCellAnalysis().run(adata, resolution=1.0)
    log(f"clusters      {result.get('n_clusters', '?')}")
    for ct in result.get("cell_types", []):
        log(f"  cell type   {ct}")
    return result


RUNNERS = {"single_cell": run_single_cell}


def execute(demo, verbose=True):
    TRANSCRIPTS.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(msg):
        lines.append(msg)
        if verbose:
            print(f"    {msg}")

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    log(f"demo          {demo.key} · {demo.subject}")
    log(f"title         {demo.title}")
    log(f"label         {demo.label}")
    log(f"started       {stamp}")

    ok, reasons = demo.check()
    for r in reasons:
        log(f"prerequisite  {r}")
    if not ok:
        log("RESULT        BLOCKED — prerequisites not met")
        (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
        return False

    if demo.runner:
        try:
            RUNNERS[demo.runner](log)
            log("RESULT        PASS — ran on real input")
        except Exception as e:  # noqa: BLE001
            log(f"RESULT        FAIL — {type(e).__name__}: {e}")
            (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
            return False
    else:
        log("RESULT        NOT IMPLEMENTED — spec only, see docs/demos/DEMO_CATALOG.md")
        (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
        return False

    (TRANSCRIPTS / f"{demo.key}.txt").write_text("\n".join(lines) + "\n")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("key", nargs="?")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--check-all", action="store_true")
    a = ap.parse_args()

    if a.list or (not a.key and not a.check_all):
        print(f"{'key':5s}{'subject':26s}{'label':16s}title")
        for d in DEMOS:
            print(f"  {d.key:3s}{d.subject:26s}{d.label:16s}{d.title}")
        return 0

    if a.check_all:
        print(f"{'key':5s}{'label':16s}{'ready':7s}why")
        ready = 0
        for d in DEMOS:
            ok, reasons = d.check()
            ready += ok
            why = "; ".join(r for r in reasons if not r.startswith("gated")) or "-"
            print(f"  {d.key:3s}{d.label:16s}{'yes' if ok else 'NO':7s}{why}")
        print(f"\n  {ready}/{len(DEMOS)} demos have their prerequisites met")
        return 0

    d = BY_KEY.get(a.key.upper())
    if not d:
        print(f"unknown demo '{a.key}' — try --list")
        return 2
    print(f"\n{d.key} · {d.title}\n")
    return 0 if execute(d) else 1


if __name__ == "__main__":
    raise SystemExit(main())
