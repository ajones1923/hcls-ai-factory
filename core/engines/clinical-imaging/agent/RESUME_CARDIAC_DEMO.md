# Resume — Clinical Imaging Engine, cardiac coronary demo

State as of **2026-08-13**, written before a reboot. Everything below is done and verified unless
it says otherwise. Nothing is committed — 73 files are modified/untracked on `main`.

## Restart the services first

Milvus runs under Docker with `restart: unless-stopped`, so it comes back by itself. The API and
the portal do not — start those two by hand.

```bash
cd /home/adam/projects/hcls-ai-factory/core/engines/clinical-imaging/agent

# Milvus — only if `docker ps` does not already show imaging-milvus-standalone healthy.
docker compose up -d milvus-etcd milvus-minio milvus-standalone

# API — port 8524. Starts in ~3 s when Milvus is up.
setsid nohup python3 -m uvicorn api.main:app --host 0.0.0.0 --port 8524 \
  > /tmp/imaging-api.log 2>&1 < /dev/null &

# Portal — Vite dev server, port 3002
cd portal && nohup npx vite --host 0.0.0.0 > /tmp/imaging-portal.log 2>&1 &
```

**Check the log line, not just the clock.** With Milvus down the API still returns `healthy` and
still serves the cardiac demo — it just reads every collection as 0. Want:

```
INFO  Milvus connected
INFO  Milvus collections ensured
```

**Warm the index before demoing.** The first `/search` after an API restart takes ~15 s while the
collections load into memory; every one after is ~100 ms. Fire one throwaway query first:

```bash
curl -s -X POST http://127.0.0.1:8524/search -H 'Content-Type: application/json' \
  -d '{"question":"warmup","top_k":1}' > /dev/null
```

Then open **http://192.168.68.107:3002/workflows** → run *Cardiac Workup: Coronary Artery Disease
Assessment* → **Advanced Imagery** tab.

Health check:

```bash
curl -s -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8524/health          # want 200
curl -s -X POST http://127.0.0.1:8524/demo-cases/DEMO-003/run | python3 -m json.tool | head -40
```

## What the demo shows (four panels, all agreeing)

| Panel | Asset | Content |
|---|---|---|
| A | `coronary_tree_annotated.png` | Coronary tree + reference heart, ring at the measured lesion, strip reads `4A/P3/HRP · 77.9% · 385` |
| B | `coronary_stenosis_detail.png` | 20 mm close-up, solid/dashed caliper key in the free corner |
| C | `coronary_tree_spin.gif` | 60 frames × 280 ms = **16.8 s/revolution** |
| D | `cardiac_pathway_dark.png` | Nano Banana Pro dark infographic: heart → calcium score 385 → LDLR → statin |

All four live in `data/demo/coronary/` and are served at `/coronary/<file>` from the API.

## Regenerate everything with one command

```bash
python3 scripts/precompute_coronary_analysis.py
```

That writes `coronary_analysis.json`, re-renders A/B/C, re-runs the illustration restyle, and
**checks panel D's baked-in numbers against the analysis**, printing a loud warning if they drift
(D is generated art — its numbers are pixels and cannot follow the data).

The prompt that produced D is kept at `data/demo/coronary/cardiac_pathway_dark.prompt.txt`; the
values painted into it are recorded in `cardiac_pathway_dark.json`.

## Key numbers (all measured, not asserted)

- Case **V1P2**, CoronariesNC6, `manualGT_rater1.stl`
- Lumen 0.60 mm / reference 2.70 mm → **77.9% diameter**, **95.1% area** stenosis
- **CAD-RADS 4A/P3/HRP** — 4A from 77.9% (70–99% band), P3 from Agatston 301–999, HRP from
  high-risk plaque features
- Agatston **385 (92nd percentile)** is REPRESENTATIVE — it needs Hounsfield units a surface mesh
  does not carry, and is labelled "(repr.)" on the panel
- Rater agreement: A 77.9% vs B 73.2%, grades agree

## Files changed in this pass

**Python**
- `src/workflows/ct_coronary_angiography.py` — CAD-RADS 4B fixed to 3-vessel **≥70%** (was ≥50%),
  left main inclusive at 50%, vessel counting folded onto parent arteries, `cadrads_report()` with
  P/HRP modifiers, `statin_intensity()` + therapy-gap finding, metric names fixed
- `src/workflows/base.py` — `apply_measurements()` hook so demo-case overrides can't clobber a real
  measurement (this is why the findings used to say 72% next to panels saying 77.9%)
- `scripts/render_coronary_mesh.py` — `TITLE_PX` shared by panels A and B, two-line title, measured
  lesion marker, legend key, 20 mm FOV, badges removed
- `scripts/restyle_cardiac_story.py` — illustration corrections (light version; the dark variant is
  superseded by the Nano Banana asset, `to_dark()` kept as fallback)
- `scripts/precompute_coronary_analysis.py` — one-command regeneration + panel D staleness guard
- `tests/test_cadrads_grading.py` — **41 new tests** pinning every fix to the guideline

**Portal**
- `portal/src/components/ChromosomeViewer.tsx` — gene cards rebuilt: MANE transcripts, full HGVS,
  GRCh38, dbSNP, evidence lines, ClinVar verify links, OMIM `*`/`#` sigils, PGx vs therapeutic
  split, LPA direction fixed, "Exemplar Variant · representative"
- `portal/src/pages/Workflows.tsx` — panels 760 px, cache-busted asset URLs, honesty banner,
  model chips marked "not run" in mock mode, latency badge → "cached result"
- `Dashboard.tsx` / `Sidebar.tsx` / `Benchmarks.tsx` — test count 1,324 → **1,365**

**Data**
- `data/reference/demo_cases.json` — DEMO-003: LAD 77.9%, RCA plaque `none` → `calcified`,
  classification `CAD-RADS 4A/P3/HRP`, medications wired in, talking points rewritten

## Verification commands

```bash
# repo gate
cd /home/adam/projects/hcls-ai-factory
ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && python3 -m pytest -q )                 # 372 pass
PYTHONPATH="$PWD/lib" python3 scripts/validate_registry.py     # note the PYTHONPATH

# engine tests
cd core/engines/clinical-imaging/agent && python3 -m pytest tests -q   # 1365 pass
```

## Open items

1. **Nothing is committed.** ~85 modified/untracked files on `main`. Branch before committing.
   `scripts/flip_chai1_live.py` is unrelated to this work — keep it out of the cardiac commit.
2. Slides 15–19 of the deck still need re-shooting from the updated portal.
4. The findings name vessels (LAD proximal, LCx, RCA) but `coronary_analysis.json` states branch
   identity is **not** resolved from the mesh — those labels come from `workflow_overrides` in
   `demo_cases.json`. Only the 77.9% magnitude is measured. Have the answer ready; a cardiologist
   will ask how LAD was derived from an unlabelled surface.

## Closed since (2026-08-13)

- **Milvus was down for two weeks**, so every collection read 0. Root cause was in the compose
  file, not the container: three `ETCD_*` env vars duplicated their own command-line flags, which
  etcd 3.5 treats as fatal. Fixed here and in **10 other compose files** carrying the same defect
  (cardiology, precision-oncology, 7 agents, the lite compose). 440 vectors now indexed across
  10 collections.
- `api/main.py` called `manager.ensure_collections()`, which does not exist — the AttributeError
  was caught and logged as "Milvus connection deferred", making a code bug look like a network
  fault. Now calls `create_all_collections()` and logs connect/setup separately.
- 9p21 card cited OMIM **611082**, which is MIAT on 22q12 — wrong locus, wrong chromosome, and a
  gene entry labelled as a phenotype. Now `*613149` (CDKN2B-AS1/ANRIL) + `#611139` (CHDS8).
- SLCO1B1 band 12p12.2 → **12p12.1**. LPA `*152200` and SLCO1B1 21.13 Mb were already correct.
- The "92nd percentile" CAC claim now reads "above the 90th percentile", noting that an exact MESA
  percentile also needs race/ethnicity, which this case does not specify.
- `npm run build` **passes**. Two of the 13 errors were a corrupt `@types/three` install (missing
  `CanvasTexture.d.ts`); the rest were unused bindings and `as Record<string, unknown>` casts that
  typed rendered values `unknown`. `DemoCaseResult` now declares `patient_demographics` and
  `clinical_scenario` properly.
- `generate_annotated_images.py` still baked "LAD Proximal: 72% Stenosis" — corrected to 77.9%.
  It generates no coronary asset today, so this was latent, not on screen.
- Panel D regenerated (Nano Banana Pro, 21:9 4K, downsampled to 2560×1086 to match the original).
  Its caption now reads `cardiac · CAD-RADS 4A/P3/HRP · statin`. The staleness guard compared only
  the bare grade, so it stayed silent while the caption sat a modifier behind every other panel —
  it now compares `cadrads_report` as well.

## Traps worth remembering

- `pkill -f "uvicorn api.main:app"` **matches its own shell** and kills the command that ran it.
  Kill by PID.
- The API takes ~15 s to start (Milvus); a health check loop that gives up early reports failure on
  a server that is actually fine.
- Assets are regenerated in place, so a browser holding the old copy shows a panel that disagrees
  with the data. The portal now cache-busts per page load; if you ever see stale panels, hard
  refresh first before believing the screen.
