# TSC Anatomical Digital Twin — Visualization Layer PRD (v0.1)

### HCLS AI Factory, Engine 7 · Surface (d) · Omniverse / OpenUSD Volumetric Twin

**Owner:** Adam M. Jones (architect, sole engineer for the v0.1 build)
**Target:** A deterministic OpenUSD authoring layer on the existing TSC Intelligence Engine that renders the engine's projections as a volumetric, time-resolved, uncertainty-aware anatomical scene — authored CPU-side on the DGX Spark, rendered with RTX in NVIDIA Omniverse on an RTX workstation or burst RunPod instance.
**License:** Apache 2.0 · **Companion to:** *The TSC Anatomical Digital Twin* (design research) and *The TSC Intelligence Engine* (master volume) · 2026

> **v0.1 = synthetic-data, atlas-anchored visualization build.** No patient imaging, no segmentation, no IRB, no institutional commitment. The anatomy is stylized reference geometry scaled by the engine's real measurements; the DICOM→segmentation→USD pipeline for patient-specific anatomy is institutional Phase-3 and explicitly out of scope. Acceptance is defined by the criteria in §6.

---

## Implementation Status — Specification, Not Yet Built (June 2026)

Unlike the engine PRD it accompanies — which has been executed and verified — **this PRD describes a build that has not yet started.** It is written to the same standard so the build is well-specified, but the honest status is:

- **The engine substrate is built.** Every input this layer consumes — the trajectory forecasts with prediction intervals, the recovered mosaic variants, the ontology-grounded HPO profiles, the threshold flags, the provenance records — exists today in the engine's projections and is covered by the engine's 58-test suite and validation harness.
- **Phase V0 is buildable on the current DGX Spark with no new hardware.** Authoring OpenUSD is a CPU-side data transformation (`usd-core`/`pxr`, or an ASCII-USD fallback with zero third-party dependency). No GPU, no Omniverse install, and no new clinical data are required to produce and test the first scene.
- **Phases V1–V3 require an RTX render target and, for V3, institutional data.** RTX rendering (Omniverse/USD Composer) runs on an RTX workstation or RunPod; patient-specific anatomy (V3) requires real imaging and segmentation and is institutional Phase-1/3 work.
- **Nothing here depends on unverified Grace-Blackwell Omniverse Kit support.** The authoring/render split (§2.4) is chosen precisely so the high-value first build does not block on it.

---

## Table of Contents

1. Overview, Goals, Scope & Users
2. System Architecture, Repo Layout & Deployment
3. Functional Requirements — SceneSpec, the OpenUSD Exporter & the Scene Schema
4. Functional Requirements — The Three Scenes, Surface (d) & Portal Preview
5. Non-Functional Requirements, Data Specs, Render Targets & Validation Harness
6. Build Plan, Test Plan, Risk Register, Dependencies & Roadmap

---


## 1. Overview, Goals, Scope & Users

### 1.1 Product summary

The **TSC Anatomical Digital Twin** is Surface (d) of the TSC Intelligence Engine: a layer that authors the engine's existing per-patient projections as OpenUSD scenes for volumetric, time-resolved viewing in NVIDIA Omniverse. It adds **no new analytics**. It is a deterministic transformation `projection → SceneSpec → .usd`, plus the three scene templates the design research specifies (lesion-trajectory twin, mosaic "powers of ten," whole-child/population atlas), plus a thin surface and API to launch and preview them.

The layer is built around one architectural commitment — **separate authoring from rendering** — so that the differentiating artifact (a correct, provenance-stamped, uncertainty-aware scene) is produced for free on the Spark's CPU, while RTX rendering happens wherever an RTX GPU is available. The layer reuses the engine's projection access, provenance model, featured-patient map, watermark discipline, and offline-first defaults.

### 1.2 Goals

1. **Render the engine's forecast with its uncertainty, volumetrically and honestly.** Author Scene 1 so the lesion grows along the engine's trajectory and is wrapped in nested 50%/90% shells whose radii equal the engine's prediction intervals — the render's crispness equals the engine's certainty.
2. **Author entirely CPU-side on the Spark, deterministically.** Same projection in → byte-identical `.usda` out, with no GPU and (in the fallback path) no third-party dependency.
3. **Carry provenance and the SYNTHETIC watermark into the scene**, so the twin is as auditable and as clearly-labeled as every other surface.
4. **Make the render a finish, not a foundation.** Value (a correct, testable scene) exists before any GPU touches it; RTX rendering and interactivity are additive phases.
5. **Generalize by configuration.** The SceneSpec and scene schema are disease-agnostic; a new disease swaps an atlas and a lesion vocabulary, not the exporter.

### 1.3 Non-goals (v0.1)

- **Not** patient-specific anatomy. No DICOM, no segmentation, no per-patient meshes; anatomy is stylized atlas geometry scaled by measurement. (Phase V3 / institutional.)
- **Not** a physiological or biomechanical growth simulation. The twin displays the engine's forecast; it does not compute one.
- **Not** a new source of clinical claims. No quantity may appear in a scene that is absent from the projection.
- **Not** dependent on Omniverse for the first build. Phase V0 authors and tests USD with no Omniverse install.
- **Not** an autonomous or diagnostic surface. Decision support only, clinician-review framing, never a diagnosis from a render.
- **Not** a live EHR/PACS integration. As with the engine, integrations are Phase-1 institutional work.

### 1.4 Users

The primary user is the TSC clinician already using the in-visit dashboard, who flies into the twin for a selected lesion. Secondary users are the multidisciplinary team and families (the twin as a shared communication object), trainees (the twin as a teaching object), and Adam (the twin as the demonstration's most memorable artifact). A technical user — a reviewer or collaborator — consumes the authored `.usda` directly and audits it as text.

### 1.5 Guiding principle

The principle inherited from the engine holds: **the build is the argument, and honesty is the build.** A render that implies more than the engine knows would betray the whole program. The uncertainty envelope, the provenance metadata, the watermark, and the determinism are therefore not features bolted on — they are the load-bearing requirements, and the rest of the scene is built around them.

---

## 2. System Architecture, Repo Layout & Deployment

### 2.1 Where the layer sits

The layer is a new subsystem alongside the engine's agents and surfaces. It reads projections through the engine's existing store and writes scenes to the artifact store; it never writes events or mutates engine state. Data flows in one direction:

```
engine projection (existing)
        │
        ▼
  SceneSpec  (engine-agnostic intermediate representation)
        │
        ├──► UsdExporter (pxr)        ──►  patient.usd / .usda   ──►  Omniverse RTX render
        └──► AsciiUsdExporter (no dep) ──►  patient.usda          ──►  (same)
```

The SceneSpec is the seam. It captures everything a scene needs (lesions with sizes/positions/time-samples/intervals/thresholds, organ illumination, variant annotation, provenance, watermark) as plain data, with no USD or rendering concept in it. Two exporters consume it: a primary `pxr`-based exporter (full USD authoring API) and an ASCII-USD fallback that emits `.usda` text directly when `usd-core` is unavailable. Both are required so the layer runs on a bare Spark.

### 2.2 Repo layout (new files only)

```
src/viz/
  __init__.py
  scene_spec.py        # SceneSpec dataclasses (engine-agnostic IR) — FR-VZ-1
  build_spec.py        # projection -> SceneSpec for each scene template — FR-VZ-2..7, FR-VZ-12..14
  atlas.py             # stylized anatomical anchors + reference geometry — FR-VZ-6
  palette.py           # colorblind-safe threshold-state palette + opacity maps — FR-VZ-5
  usd_pxr.py           # pxr-based exporter (time-samples, customData, materials) — FR-VZ-10a
  usd_ascii.py         # dependency-free .usda text exporter — FR-VZ-10b
  export.py            # façade: pick exporter, write scene, return path — FR-VZ-10
  provenance_usd.py    # projection provenance -> USD customData + /Provenance scope — FR-VZ-8
  watermark.py         # in-scene SYNTHETIC watermark prim + HUD spec — FR-VZ-9
api/routes/viz.py      # GET /viz/{scene}/{patient_id} -> authored scene (path/inline) — FR-VZ-15
app/twin_app.py        # Surface (d): launch/preview, patient+lesion picker — §4
scripts/export_usd.py  # CLI: author scenes for featured patients / whole cohort — §6
omniverse/             # (Phase V1+) Kit extension + USD Composer project + MDL materials
  ext_tsc_twin/        # optional Omniverse Kit extension (TSC panel, scrub, toggles)
tests/viz/
  test_scene_spec.py   # IR invariants
  test_usd_export.py   # determinism, geometry scale, envelope==interval, provenance, watermark
  test_viz_api.py      # endpoint contract
```

This mirrors the engine's conventions: an engine-agnostic core (`scene_spec`), deterministic builders, a runs-here fallback (`usd_ascii`) for every production path (`usd_pxr`), tests under a dedicated package, and a CLI in `scripts/`.

### 2.3 Dependencies

| Dependency | Role | Required for | Fallback |
| --- | --- | --- | --- |
| `usd-core` (`pxr`) | Full USD authoring (time-samples, customData, materials, USDZ) | Phase V0 primary path, USDZ (V2) | `usd_ascii.py` ASCII `.usda` writer (zero-dep) |
| NumPy (already present) | Interval→radius math, positions | All phases | — |
| NVIDIA Omniverse / USD Composer | RTX rendering, MDL materials, interactive scrub | V1+ (render only) | Author-only; view in any USD viewer |
| Omniverse Kit | Custom TSC panel extension | V1 polish (optional) | Native USD timeline |
| MONAI + imaging stack | DICOM→segmentation→USD | V3 (institutional) | Atlas geometry |

The only new hard dependency for the first build is `usd-core`, and even it is optional via the ASCII fallback. Everything GPU-bound is render-side and out of the authoring path.

### 2.4 Deployment — the authoring/render split

- **Author on the Spark.** `python scripts/export_usd.py` (or `GET /viz/...`) produces `.usda`/`.usd` scenes from projections, CPU-only, offline-capable, deterministic. This is the whole of Phase V0.
- **Render off-box.** Open the authored scene in Omniverse/USD Composer on an RTX workstation or a RunPod RTX instance for real-time navigation, path-traced stills/film, and (V2) USDZ export for tablet AR.
- **Settings.** New `TSC_VIZ_*` settings (output dir, scene-unit scale, exporter preference `auto|pxr|ascii`, watermark text) extend the engine's pydantic settings with the `TSC_` env prefix. Authoring honors `TSC_OFFLINE` (no network ever needed to author).

### 2.5 Reuse of the engine

The layer imports the engine's projection access (`orch.store.projection`), featured-patient map, provenance records (already on every event), artifact store (for writing scenes), and watermark/settings conventions. It adds no agent, no event type, and no change to engine state — it is a pure consumer, which keeps the engine's 58-test suite untouched and green.

---

## 3. Functional Requirements — SceneSpec, the OpenUSD Exporter & the Scene Schema

> Requirements are labeled **FR-VZ-n**. "Real in v0.1" means implemented and tested CPU-side; "render-side" means the requirement is satisfied by the authored scene and realized visually only when rendered.

### FR-VZ-1 — SceneSpec intermediate representation (real)

A small set of dataclasses captures a scene as engine-agnostic data: `SceneSpec(patient_id, scene_kind, watermark, provenance, time_range, lesions[], organs[], annotations[])`; `LesionSpec(label, anchor, observed[(month,cm)], forecast{m:(mean,pi50_lo,pi50_hi,pi90_lo,pi90_hi)}, threshold_cm, time_to_threshold, state_by_time)`; `OrganSpec(name, anchor, illuminated, source_hpo[])`; `Annotation(text, target, provenance_ref)`. The SceneSpec contains no USD or rendering concept and is fully serializable (for golden-file tests and for the API).

### FR-VZ-2 — lesion geometry from the trajectory (real)

For each lesion in the projection's trajectory, the builder emits a `LesionSpec` whose geometry is a sphere (default) or ellipsoid with radius = (measured/forecast diameter ÷ 2) × `scene_unit_scale`, positioned at the lesion's anatomical anchor (FR-VZ-6). Observed sizes and forecast means become the radius time-samples (FR-VZ-3). No size appears that is not in the projection.

### FR-VZ-3 — time-sampled growth animation (real)

The lesion radius and color are authored as **USD time-samples** across a timeline whose codes map to months (observed −24…0 and forecast +6/+12/+18). The stage declares `startTimeCode`, `endTimeCode`, and `timeCodesPerSecond`. Scrubbing the timeline interpolates the lesion through its history into the forecast. The mapping is deterministic and documented in the scene's `customData`.

### FR-VZ-4 — the uncertainty envelope (real; the load-bearing requirement)

Every forecast lesion is wrapped in two concentric translucent shells: an inner shell at the 50% prediction-interval radius and an outer shell at the 90% radius, with their radii **time-sampled to equal the engine's interval bounds at each forecast code**. Opacity maps inversely to interval width within configured bounds (wider interval → softer, larger cloud). **A forecast lesion may not be authored without this envelope** — the exporter raises if asked to. This requirement encodes the program's core safety inversion: the render is exactly as crisp as the engine is certain.

### FR-VZ-5 — threshold surface and color-state machine (real)

Each lesion with a threshold gets a translucent threshold shell at the threshold radius and a color state in {`below`, `approaching`, `at_or_above`} computed from the engine's threshold flags (`crosses_in_12_18mo_window`, `at_or_above_threshold`, `months_to_threshold`) at each time code. Colors come from a colorblind-safe palette (FR-VZ via `palette.py`). The state and color are time-sampled so the lesion visibly changes as it nears the boundary.

### FR-VZ-6 — anatomical atlas anchors (real)

`atlas.py` provides canonical anchors and stylized reference geometry: a low-poly brain with ventricular landmarks including the **foramen of Monro** (SEGA anchor), stylized **kidneys** (AML anchor), **cortex** (tuber anchor), **skin/body shell**, and **heart** (rhabdomyoma anchor). Anchors are positions+orientations in scene space; geometry is reference meshes loaded from `omniverse/atlas/` or generated procedurally. The atlas is a **sublayer** so Phase V3 can replace it with segmented anatomy without touching the lesion/envelope/provenance logic.

### FR-VZ-7 — HPO→organ illumination (real)

For the whole-child scene, each organ's `illuminated` flag and intensity derive from the patient's ontology-grounded HPO profile (e.g., a renal-angiomyolipoma term illuminates the kidney). The mapping table is explicit and disease-swappable. Only organs implied by present HPO terms are illuminated.

### FR-VZ-8 — provenance travels into the scene (real)

Each renderable prim carries, as USD `customData`, the provenance of the projection field that produced it (event type, model id, tier, latency, prompt-template version, input hash). A dedicated `/Provenance` scope in the scene mirrors the engine's audit trail so a viewer can select an element and read its source. Provenance completeness in the scene is asserted by test (FR-VZ parity with the engine's 100% provenance metric).

### FR-VZ-9 — in-scene SYNTHETIC watermark (real)

A persistent watermark prim (textual + a HUD overlay spec) marks the scene as synthetic demonstration data, with the engine's standard "decision support, clinician review required, not FDA-cleared" framing. The watermark is authored into the base layer and is part of every scene; removing it is not a supported operation.

### FR-VZ-10 — dual exporter with runs-here fallback (real)

`export.py` selects an exporter by `TSC_VIZ_EXPORTER` (`auto` default): **FR-VZ-10a** `usd_pxr.py` uses the `pxr` API (UsdGeom prims, `Set(value, time)` time-samples, `customData`, UsdShade/MDL material bindings, USDZ packaging) when `usd-core` is importable; **FR-VZ-10b** `usd_ascii.py` emits well-formed `.usda` text (prims, attributes, time-samples, `customData`) with zero third-party dependency. Both satisfy FR-VZ-2..9 and produce equivalent scenes; the ASCII path guarantees the layer authors on a bare Spark.

### FR-VZ-11 — deterministic, diffable output (real)

Given a fixed projection and settings, the exporter produces a byte-identical scene (stable prim ordering, fixed float formatting, no timestamps or randomness in the authored file). Scenes are therefore reviewable and testable like source code via golden files.

### FR-VZ-12 — variant annotation prim (real)

For the mosaic scene, the recovered variant is authored as an annotation prim carrying gene, HGVS, VAF, mosaic status, ACMG criteria and classification, and the recovery narrative — all from the projection — with provenance. The cellular-field VAF visualization (FR-VZ-13) references this prim.

### FR-VZ-13 — VAF cellular field (real geometry, render-side effect)

The mosaic scene authors a field of instanced cell prims of which exactly the recovered VAF fraction are flagged "variant-carrying" (e.g., ~1 in 12 for an 8.3% VAF), via USD point-instancing. The fraction equals the projection's VAF to instancing resolution; the glow is a render-side material effect.

### FR-VZ-14 — population array (real)

The population scene authors the cohort as an array of ghosted body prims keyed to each patient's classification and recovery status, with the recovered mosaics flagged — the engine's `/cohort` aggregation, rendered. Layout is deterministic.

### FR-VZ-15 — viz API (real)

`GET /viz/{scene}/{patient_id}` authors (or returns a cached) scene for the patient and returns the artifact path plus the SceneSpec JSON; `scene ∈ {lesion, mosaic, atlas}` and a `population` variant ignores `patient_id`. The route reads engine state via `request.app.state` exactly as the cohort/eval routes do. Errors are explicit for unknown scene/patient.

---

## 4. Functional Requirements — The Three Scenes, Surface (d) & Portal Preview

### FR-VZ-16 — Scene 1: lesion-trajectory twin (real authoring; RTX render V1)

Compose FR-VZ-2..9 for a selected patient+lesion into the centerpiece scene: atlas-anchored lesion, time-sampled growth across observed→forecast, the 50%/90% uncertainty envelope, the threshold shell and color state, a time-to-threshold marker, and the 12–18-month window highlighted on the timeline. Default targets: Patient B (SEGA, foramen of Monro) and Patient C (AML, ~4 cm). This is the Phase V0 deliverable.

### FR-VZ-17 — Scene 2: mosaic "powers of ten" (real authoring; render V1)

Compose the organ→cellular-field→genome→variant scale traversal using FR-VZ-12/13, with camera keyframes authored as USD and the ACMG annotation surfaced. Default target: Patient A (8.3% VAF *TSC2* frameshift in tuber tissue).

### FR-VZ-18 — Scene 3: whole-child atlas + population (real authoring; render V1)

Compose the whole-child organ illumination (FR-VZ-7) with the ITSC surveillance overlay (from the engine's gap analysis) and the population array (FR-VZ-14). Default: any featured patient for the body view; the full cohort for the population view.

### FR-VZ-19 — Surface (d): the twin app (real)

`app/twin_app.py` is a thin surface that (a) lets the user pick a patient, lesion, and scene; (b) authors the scene on demand via the viz layer; (c) provides the **fly-in from the dashboard** (a link/handoff from the in-visit dashboard's 2D trajectory to the twin of the same lesion); and (d) offers download of the `.usda`/`.usdz` and a preview (FR-VZ-20). It carries the engine's watermark and offline/live indicator and never computes clinical conclusions.

### FR-VZ-20 — portal preview with honest degradation (real)

When Omniverse is not present, the surface still shows *something* truthful: a pre-rendered turntable image/clip of the scene if one exists, or a constrained web-3D preview (a lightweight viewer of the authored geometry), each clearly labeled as a preview of an atlas-anchored synthetic scene. The preview never implies RTX fidelity it is not showing, and never implies patient-specific anatomy.

### FR-VZ-21 — Omniverse Kit extension (optional, V1 polish)

An optional Kit extension (`omniverse/ext_tsc_twin/`) adds a purpose-built TSC panel — patient/lesion picker, timeline scrub with the threshold window marked, organ toggles, and a provenance inspector that reads the `/Provenance` scope. The extension is additive; all scenes are fully usable with the native USD timeline without it.

---

## 5. Non-Functional Requirements, Data Specs, Render Targets & Validation

### 5.1 Non-functional requirements

| ID | Requirement |
| --- | --- |
| NFR-V-1 | **Author CPU-side, offline.** Scene authoring requires no GPU and no network; honors `TSC_OFFLINE`. |
| NFR-V-2 | **Deterministic.** Same projection+settings → byte-identical scene (FR-VZ-11). |
| NFR-V-3 | **Performance.** Authoring a single-lesion scene completes in well under a second; the whole featured set and the population scene in a few seconds, on the Spark CPU. |
| NFR-V-4 | **Prim budget.** Scenes stay within a documented prim/instance budget so RTX real-time navigation is smooth on a single RTX GPU. |
| NFR-V-5 | **Accessibility.** The threshold palette is colorblind-safe; state is also encoded by a non-color channel (label/opacity) so color is never the sole carrier. |
| NFR-V-6 | **Auditability.** 100% of renderable elements carry provenance metadata (FR-VZ-8); the watermark is present in every scene (FR-VZ-9). |
| NFR-V-7 | **No-new-claim invariant.** No authored quantity may be absent from the projection; enforced by the build_spec layer and asserted by test. |
| NFR-V-8 | **Portability.** Authored `.usda` opens in any compliant USD viewer; `.usdz` (V2) opens in tablet AR. |

### 5.2 Data specs — scene units and mappings

- **Scene-unit scale.** Centimeters map to scene units by `TSC_VIZ_SCALE` (default chosen so a few-centimeter lesion is comfortably navigable); the scale is recorded in `customData`.
- **Time mapping.** Month *m* maps to time code *f(m)*; `f` is linear and documented in the scene; the forecast window 12–18 months is tagged as a named interval for the UI.
- **Interval→radius.** The 50%/90% shell radii equal (interval bound ÷ 2) × scale at each forecast code; opacity is a monotone function of interval width within `[opacity_min, opacity_max]`.
- **Threshold state.** Derived solely from the engine's flags; no independent thresholding in the viz layer.

### 5.3 Render targets

| Target | Phase | Use |
| --- | --- | --- |
| Any USD viewer (author check) | V0 | Validate geometry/animation without RTX |
| Omniverse / USD Composer (RTX, MDL) | V1 | Real-time navigation; path-traced stills and film |
| RunPod RTX instance | V1 | Burst rendering when no local RTX |
| USDZ on tablet (ARKit) | V2 | In-clinic AR; family communication |
| In-portal web preview | V2 | Honest degraded preview when Omniverse absent |

### 5.4 Validation harness (extends the engine's harness)

Because the twin is a deterministic function of the projection, it is tested as software. `tests/viz/` asserts:

- **Determinism** (NFR-V-2): re-authoring yields a byte-identical file (golden-file diff).
- **Geometry-scale fidelity** (FR-VZ-2): rendered lesion radius equals measured size under scale, to tolerance.
- **Envelope-equals-interval** (FR-VZ-4): the 50%/90% shell radii at each forecast code equal the engine's prediction-interval bounds — the visual uncertainty *is* the numeric uncertainty. This is the headline visualization test.
- **Threshold-state correctness** (FR-VZ-5): color state matches the engine's threshold flags at every time code.
- **Provenance completeness** (FR-VZ-8): every renderable element carries provenance; the `/Provenance` scope mirrors the audit trail.
- **Watermark presence** (FR-VZ-9): every authored scene contains the watermark prim.
- **No-new-claim** (NFR-V-7): every authored numeric equals a projection field.
- **Exporter parity** (FR-VZ-10): the pxr and ASCII exporters produce equivalent scenes for the same SceneSpec.
- **API contract** (FR-VZ-15): the route returns a valid scene + SceneSpec and errors explicitly.

Clinical-impact validation (does the twin change decision confidence/time/communication?) is a forward usability study, distinguished — as in the engine — from construct-level correctness, which the tests above establish now.

---

## 6. Build Plan, Test Plan, Risk Register, Dependencies & Roadmap

### 6.1 Phased build plan

**Phase V0 — CPU-side authoring (buildable now, no new hardware).**
1. `scene_spec.py` (FR-VZ-1) and `build_spec.py` for the lesion scene (FR-VZ-2..5, 8, 9, 11).
2. `atlas.py` anchors (FR-VZ-6, SEGA + AML first), `palette.py` (FR-VZ-5), `watermark.py` (FR-VZ-9), `provenance_usd.py` (FR-VZ-8).
3. `usd_ascii.py` (FR-VZ-10b) first — it has no dependency and unblocks tests — then `usd_pxr.py` (FR-VZ-10a).
4. `export.py` façade (FR-VZ-10), `api/routes/viz.py` (FR-VZ-15), `scripts/export_usd.py`.
5. `tests/viz/` (§5.4) including the envelope-equals-interval test. **Exit criterion:** Scene 1 authors for Patients B and C on the Spark, opens in a USD viewer, and passes the full viz test suite.

**Phase V1 — RTX rendering and the remaining scenes.** Omniverse/USD Composer on an RTX target; MDL materials for the envelope/threshold/glow; Scenes 2 (FR-VZ-17) and 3 (FR-VZ-18); the timeline scrub as a polished experience; optional Kit extension (FR-VZ-21).

**Phase V2 — AR and portal preview.** USDZ packaging (FR-VZ-20/NFR-V-8); the in-portal degraded preview; surface (d) polish (FR-VZ-19/20).

**Phase V3 — patient-specific anatomy (institutional).** DICOM → MONAI segmentation → USD; replace the atlas sublayer with segmented meshes; the scene schema is unchanged.

### 6.2 Acceptance criteria

| AC | Criterion | Phase |
| --- | --- | --- |
| AC-V1 | Scene 1 authors deterministically for Patients B and C, CPU-side, offline | V0 |
| AC-V2 | The uncertainty envelope's 50%/90% radii equal the engine's prediction intervals (tested) | V0 |
| AC-V3 | Every authored scene carries the watermark and 100% element provenance (tested) | V0 |
| AC-V4 | The ASCII and pxr exporters produce equivalent scenes (tested) | V0 |
| AC-V5 | `GET /viz/...` returns a valid scene + SceneSpec; CLI authors the featured set | V0 |
| AC-V6 | Scenes render with RTX in Omniverse with the envelope/threshold reading correctly | V1 |
| AC-V7 | Scenes 2 and 3 author and render | V1 |
| AC-V8 | USDZ opens in tablet AR; portal preview degrades honestly | V2 |

AC-V1..V5 are buildable now; AC-V6..V8 require an RTX target / tablet; patient-specific anatomy is explicitly beyond these criteria.

### 6.3 Risk register

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| Grace-Blackwell Omniverse Kit support unverified | Medium | Low | Authoring is CPU-side; rendering is off-box by design — no core dependency on Spark Kit |
| Render implies false precision | Medium | **High** | Uncertainty envelope is mandatory on every forecast lesion (FR-VZ-4); palette calibrated to numbers; tested |
| Atlas mistaken for patient anatomy | Medium | High | Stylized geometry visually distinct; scenes labeled atlas-anchored; watermark always present |
| `usd-core` unavailable on target | Low | Low | ASCII-USD fallback (FR-VZ-10b) has zero dependency |
| Scope creep into physiology simulation | Medium | Medium | Non-goal stated; twin defined as a pure function of the projection |
| Prim/instance budget hurts real-time perf | Low | Medium | NFR-V-4 budget; point-instancing for cellular/population fields |
| Visualization drifts from analytics | Low | High | Determinism + envelope-equals-interval + no-new-claim tests bind the render to the math |

### 6.4 Dependencies

`usd-core` (optional; ASCII fallback otherwise); NumPy (present); the engine's projection/provenance/featured/artifact APIs (present). Render-side: NVIDIA Omniverse / USD Composer on an RTX GPU (V1); a RunPod RTX instance (optional burst); a tablet for AR (V2). Institutional (V3): imaging access and a MONAI-based segmentation pipeline.

### 6.5 Roadmap and generalization

After V0–V2, the layer generalizes by configuration: a new disease provides an atlas sublayer and a lesion/organ vocabulary; the SceneSpec, exporters, envelope, threshold, provenance, and watermark are reused unchanged — extending the engine's NF1/NF2/Rett/Williams/mTORopathy replication thesis from analytics to visualization. The long-term direction is the convergence of V3 patient-specific anatomy with the engine's analytic forecasts, so that a clinician sees *this child's* segmented lesion grow along *this engine's* calibrated, uncertainty-bounded forecast — the data-driven anatomical twin in full, still honest, still auditable, still running first on a single Spark.

### 6.6 Summary

This PRD specifies a visualization layer that adds no analytics, depends on no unverified hardware for its first build, authors deterministically and freely on the existing Spark, binds every rendered element to a validated projection and its provenance, and makes its single most important property — that the render is exactly as confident as the engine — a tested, load-bearing requirement. It is the surface that turns the engine's rigor into something a clinician, a team, and a family can see together, without ever letting the picture outrun the proof.
