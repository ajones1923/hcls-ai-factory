# TSC Anatomical Digital Twin — Spark → RunPod → Omniverse render guide

The demo's "wow": the engine's **real forecast for a child's lesion**, rendered as a
volumetric, time-scrubbed anatomical twin with the uncertainty made visible. The DGX Spark
runs the engine and **authors** the scene (free, CPU); a RunPod RTX instance **renders** it.

## The split (why it's cheap and reproducible)

| Stage | Where | GPU | Output |
| --- | --- | --- | --- |
| Engine + analytics | DGX Spark | no | projections (validated, provenance-stamped) |
| Scene authoring | DGX Spark | **no** | deterministic `.usda` (this layer) |
| RTX render | RunPod RTX pod | yes | real-time navigation · path-traced stills/film |

Nothing here depends on Omniverse running on the Spark — authoring is pure Python.

## Two scenes

- **Lesion trajectory twin** (`<pid>_lesion_trajectory.usda`) — the growing lesion + uncertainty
  envelope + threshold membrane, scrubbable across the forecast. Default: Patient B (SEGA), C (AML).
- **Mosaic "powers-of-ten"** (`<pid>_mosaic.usda`) — a `PointInstancer` cell field where exactly
  the recovered VAF fraction of cells (e.g. **12 of 144 ≈ 8.3%** for Patient A) carries the glowing
  variant prototype; the variant call (gene/HGVS/ACMG/recovered-from-tissue) rides as `customData`.
  This is the Act-One showpiece: the signal a blood test calls negative, made countable.
- **Whole-child atlas** (`<pid>_atlas.usda`) — a stylized body with the TSC-affected organ systems
  (brain/kidney/skin/heart/lung) **lit by the patient's HPO profile**; the multisystem burden at a glance.
- **Population** (`cohort_population.usda`) — all **50 patients** as a grid of whole-child figures,
  body colour = ACMG class, organs lit by phenome, and the **7 recovered mosaics ringed in a gold
  halo**. The scale story (`GET /viz/population`): one Spark, fifty children, a whole disease.

## 1. Author scenes on the Spark

```bash
# featured A/B/C  (or --all for the whole cohort)
venv/bin/python scripts/export_usd.py
# scenes land in data/usd/<patient>_{lesion_trajectory,mosaic}.usda
```

Or pull one over the API: `GET /viz/lesion/<pid>?inline=true` or `GET /viz/mosaic/<pid>?inline=true`.

## 2. Render on RunPod (RTX)

1. Launch a RunPod pod with an RTX GPU and the **NVIDIA Omniverse / USD Composer** (or a
   Kit-based) image.
2. Copy the `.usda` (e.g. `TSC-0001_lesion_trajectory.usda`) to the pod.
3. Open it in USD Composer. Set the renderer to **RTX – Interactive (Path Tracing)** for
   reference quality, or **RTX – Real-Time** to navigate.
4. **Scrub the timeline** from frame −24 (two years ago) through 0 (now) to +18 (the
   forecast). Watch the SEGA grow along the engine's trajectory, the **uncertainty envelope**
   (the two translucent shells) widen into the forecast, and the lesion's colour cross from
   teal → amber → vermillion as it approaches the threshold membrane.

## What you are looking at (and what is honest about it)

- **The lesion** (`Lesion` sphere) — radius time-sampled to the engine's measured sizes
  (−24/−12/−6 mo) and forecast means (+6/+12/+18 mo). Colour = threshold state.
- **The uncertainty envelope** (`Envelope50` / `Envelope90`) — the shell radii **are** the
  engine's 50%/90% prediction-interval bounds. The cloud is exactly as wide as the engine is
  uncertain; it cannot imply false precision. (Asserted by `tests/viz`.)
- **The threshold membrane** (`ThresholdShell`) — the SEGA discussion line (1.8 cm) / AML
  bleeding-risk line (4.0 cm).
- **Provenance + watermark** — stage/prim `customData`; every element traces to a validated
  projection. The scene is labelled SYNTHETIC, atlas-anchored (not patient imaging), decision
  support, not FDA-cleared.

## Polish (later phases, per the Digital Twin PRD)

- **V1** — MDL glass/subsurface materials for true volumetric translucency; the mosaic
  "powers-of-ten" and whole-child/population scenes; an optional Omniverse Kit panel.
- **V2** — USDZ export for tablet AR in clinic; an in-portal degraded web preview.
- **V3 (institutional)** — replace the atlas sublayer with **segmented patient anatomy**
  (DICOM → MONAI → USD). The scene schema does not change — only the source of the geometry.
