# TSC Intelligence Engine — Cincinnati Children's Demo Runbook

### A three-act walkthrough · one DGX Spark + RunPod (Omniverse RTX) · SYNTHETIC data

**Owner:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Companion to:** the master research volume, the engine PRD, and the Digital Twin docs
**Live driver:** `venv/bin/python scripts/demo_runbook.py` (prints this with *live* numbers + authors every scene)
**License:** Apache 2.0 · 2026

> **What this is.** The script you walk Cincinnati Children's through. Three acts — mosaic
> recovery, the longitudinal twin, scale & infrastructure — each with the story, exactly what
> to show, the measured number, and the line to say. Everything runs on the Spark on synthetic
> data ($0); the film-quality render is a RunPod RTX pod, off-box. Decision support, clinician
> review required, not FDA-cleared.

---

## Pre-flight (before they arrive)

| Step | Command / action |
| --- | --- |
| Start the engine + surfaces | `bash scripts/serve.sh` → portal `http://<spark>:8560` |
| Author every USD scene | `venv/bin/python scripts/export_usd.py --all` → `data/usd/*.usda` |
| Stage the render | open the scenes in Omniverse / USD Composer on a RunPod RTX pod (checklist below) |
| Sanity check | 50 patients enrolled · reasoning offline ($0) · all four ports HTTP 200 |
| Print the live runbook | `venv/bin/python scripts/demo_runbook.py` (numbers pulled live) |

Keep one browser on the **dashboard** (`:8562`) and one on the **portal** (`:8560`); keep the
Omniverse viewport ready on the RunPod pod.

---

## Act One — Mosaic recovery: ending the diagnostic odyssey (Patient A)

**The story.** Patient A is four years old. Standard blood testing came back *no mutation
identified* — the ~15% of TSC children who never get a molecular diagnosis. We sequence the
resected tuber tissue instead.

**What to show.**

- **Dashboard → Patient A** (`:8562`): the Variant Curator recovers a *TSC2* frameshift at
  ~8.3% VAF (mosaic), classifies it **Likely Pathogenic** (PVS1 applied Strong for the low-VAF
  mosaic + PM2 + PP4), and recommends ddPCR. **Open the audit trail** — model, tier, latency,
  prompt version on every step.
- **Omniverse → mosaic scene** (`TSC-0043_mosaic.usda`): a field of cells in which exactly
  **~1 in 12 glows** — 8.3% VAF made *countable*. The variant cells are emissive (MDL); the
  call rides along as metadata. This is the signal a blood test reports as negative.

**The measured number (this cohort).**

| Metric | Value |
| --- | --- |
| Variant detection (standard → engine) | 86% → **98%** (+12 pts) |
| Sub-threshold mosaics recovered | **6 / 6** (sensitivity 100%) |
| Truncating variant called benign | **0** (a strand-biased artifact is rejected every sample) |

**The line to say.** *"This is the capability your biobank unlocks — banked surgical tissue
becomes a molecular diagnosis the blood test missed. On real specimens this is a retrospective
NMI re-analysis, on tissue you already have in the freezer."*

---

## Act Two — The longitudinal twin: seeing the future of a lesion (B & C)

**The story.** Patient B (12) has a SEGA near the foramen of Monro; the question every visit is
*when* it reaches the size the team discusses intervening — and *how sure* are we. Patient C
(18) has a growing renal AML, declining eGFR, and refractory seizures — three trajectories at
once.

**What to show.**

- **Dashboard → Patient B** (`:8562`): the four-quadrant view (variant · phenotype · trajectory
  · TAND). The 2-D chart shows the SEGA crossing the discussion line in the **12–18-month
  window** with its 50/90% bands.
- **Omniverse → lesion twin** (`TSC-0001_lesion_trajectory.usda`): **scrub the timeline −24 → +18
  months.** The lesion grows along the forecast inside a **glass uncertainty envelope whose radii
  *are* the 50/90% prediction intervals** (MDL OmniGlass), changing colour as it crosses the red
  threshold membrane. The render is exactly as confident as the engine.
- **Dashboard → Patient C**: eGFR and seizure-frequency tabs, surveillance cadence **tightened
  vs. the ITSC floor**, and the six-section therapeutics brief (every claim source-attributed).

**The line to say.** *"The envelope can't lie — it widens exactly as the model is uncertain.
That's what makes a forecast safe to put in front of a family."*

---

## Act Three — Scale & infrastructure: one Spark, a whole disease

**What to show.**

- **Portal** (`:8560`): the **population command center** (all 50 patients), the **validation
  scorecard**, and **`/system/usage`** (the cost ledger — **$0** offline).
- **Omniverse → population scene** (`cohort_population.usda`): 50 whole-child figures, body
  colour = ACMG class, organs lit by the phenome, and the **7 recovered mosaics ringed in a gold
  emissive halo**. The scale story, in 3-D.

**The numbers (validation scorecard — construct validity on synthetic ground truth).**

| Metric | Value |
| --- | --- |
| ACMG classification accuracy | 100% (50/50) |
| Diagnostic detection uplift · mosaic recovery | +12 pts · 100% |
| Truncating variant called benign · provenance complete | 0 · 100% |
| Classification distribution | 38 Pathogenic · 10 VUS · 1 LP · 1 NMI |

**The thesis.** *"Everything you just saw ran on a single $4,699 DGX Spark, on synthetic data,
for $0 — and the film-quality render is a RunPod RTX pod, off-box. This is the engine. Phase 2
is your institution: the **Discover Together Biobank** as the substrate, the **Winslow Pavilion**
as the envelope, **Dr. Hagedorn's BMI** as the methodology, the **TSC clinic** as the patients,
**Epic + the LIMS** as the plumbing. Point us at the biobank and we'll show whether the recovery
yield holds on your specimens."*

**The honest caveat (say it out loud).** *"These are construct-validity metrics on synthetic
data with known ground truth — they prove the logic recovers the planted signal, not prospective
clinical accuracy. The real test is running the Variant Curator on actual Cincinnati specimens.
That's exactly the Phase-1 study we're proposing."*

---

## RunPod RTX render checklist (Omniverse)

1. Launch a RunPod pod with an RTX GPU and **NVIDIA Omniverse / USD Composer**.
2. Copy `data/usd/*.usda` to the pod (four scene kinds: `mosaic`, `lesion_trajectory`, `atlas`,
   `population`).
3. Open a scene; set the renderer to **RTX – Interactive (Path Tracing)**.
4. MDL is already wired — envelopes render as **glass**, variant cells and recovery halos as
   **emissive glow**. (Modules `OmniGlass.mdl` / `OmniPBR.mdl` resolve inside Omniverse.)
5. For the lesion twin, **scrub the timeline** (frame −24 → +18). For the population, fly the grid
   and land on a gold-haloed recovery.

**Authoring is CPU-side on the Spark; rendering is the RunPod RTX step. Nothing depends on
Omniverse running on the Spark.**

---

## One-line cheat sheet

- **Act One:** lead with the *tissue* — mosaic recovery, +12-pt detection, the glowing cell field.
- **Act Two:** scrub the *twin* — the glass envelope is the prediction interval; it can't lie.
- **Act Three:** show the *scale* — 50 children on one Spark, 7 gold halos, $0, and the honest caveat.
- **Close:** *"Lead with the tissue. Show the glass. Tell the truth."*
