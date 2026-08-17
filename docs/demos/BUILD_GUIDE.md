# Build Guide — the 17 Demonstrations

Companion to [`DEMO_CATALOG.md`](DEMO_CATALOG.md) (specs) and [`PRD.md`](PRD.md) (requirements).
Keys are **E1–E8 / A1–A8 / P1** — not D1–D7, which is the patient-story portfolio.

---

## 0. Check what is ready

```bash
.venv/bin/python scripts/run_demo.py --list        # all 17
.venv/bin/python scripts/run_demo.py --check-all   # prerequisites, runs nothing
```

Measured 2026-08-15: **7 of 17 prerequisites met**, one (**E8**) proven end-to-end.

The runner enforces PRD DR-3: **a demo declared LIVE whose service is unreachable fails.** It does
not fall back to a canned result. That is deliberate — the same discipline caught two capabilities
registered `live` in the registry with nothing bound to their ports.

## 1. Run the one that already works

```bash
.venv/bin/python scripts/run_demo.py E8
```

Real output, no gated software, no GPU:

```
loaded        2,700 cells x 32,738 genes
running       QC -> normalize -> HVG -> PCA -> neighbors -> Leiden -> marker DE
clusters      9
  cell type   B cells / CD14+ Monocytes / CD4 T cells / Dendritic cells
              FCGR3A+ Monocytes / Megakaryocytes / NK cells
RESULT        PASS — ran on real input
```

Transcript: `demo/transcripts/E8.txt`. Diff it after any change — a silent regression in the
clustering shows up here.

Prerequisite (ungated): `.venv/bin/pip install scanpy anndata leidenalg igraph`

## 2. Unblock the eleven agent/engine demos — cheapest win available

E2, E5, E6 and A1–A8 need only two things:

```bash
# 1. an API key
cp .env.example .env && ${EDITOR:-nano} .env      # set ANTHROPIC_API_KEY

# 2. their service running
docker compose -f docker-compose.dgx-spark.yml up -d milvus-etcd milvus-minio milvus
docker compose -f docker-compose.dgx-spark.yml up -d
.venv/bin/python scripts/run_demo.py --check-all
```

**Seed Milvus before demoing an agent.** An empty collection makes a working agent look broken. If
no seeded corpus exists, `core/*/scripts/seed_*.py` and `ingest_*.py` are the entry points.

Ports follow the convention in [`../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md) — registry endpoint is
the **UI**, API is **UI + 1**. The runner probes the API port.

## 3. The gated demos — E1, E3, E7

These stay **REPRESENTATIVE** until `../build/GATED_SOFTWARE_BUILD_GUIDE.md` is executed:

| Demo | Needs | Becomes |
|---|---|---|
| E1 Genomic Foundation | Parabricks (G2) | LIVE — real FASTQ→BAM→VCF |
| E3 Therapeutic Discovery | MolMIM + DiffDock (G3/G4) | LIVE, or **BURST** if the NIMs are x86-only |
| E7 Structural Biology | CUDA torch + ESMFold (G1/G5) | LIVE — prediction in the room |

**Do not relabel any of these by editing the catalogue.** Change the label only after the runner
passes with the service actually answering.

## 4. Adding or changing a demo

1. Add a `Demo(...)` entry in `scripts/run_demo.py` with its key, label, port, packages and payload.
2. If it can execute, write a runner function and register it in `RUNNERS`.
3. Run it; commit the transcript.
4. Update `DEMO_CATALOG.md` — story first, honesty flags inline, not in a footer.

**The label is a claim about reality, not a category.** LIVE means it ran, now, on real input, in
front of the audience.

## 5. Running a session

The catalogue's recommended arc for a room:

```
P1  (weight — one child, the whole factory)
 └─► E4  (proof — imaging, already working end-to-end on this box)
      └─► E8  (something completely real, no caveats)
           └─► the audience's agent  (A3 pharmacogenomics · A7 rare disease · A1 CAR-T)
                └─► A6  (close the loop — hand a molecule to the trial agent)
```

**Strongest three today:** E4 (proven), E8 (entirely real), A3 (44 verified CPIC pairs).
**Handle with most care:** E3 — the flagship, the most watched, and currently the least live.
