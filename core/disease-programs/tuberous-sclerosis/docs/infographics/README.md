# TSC Engine (Engine 7) — Nano Banana Pro infographic prompts

Copy-paste prompts for generating infographics of the TSC Intelligence Engine with Nano Banana
Pro. The **filenames are numbered `01`–`18` in walkthrough order** (orientation → disease →
problem → engine → capabilities → visuals → proof → demo → strategy), so they sort and read as
a sequence. Each generated image saves alongside its prompt with the same number (e.g.,
`05_engine7_portrait.md` → `05_engine7_portrait.png`).

**All 18 are now a cohesive storyboard:** every prompt is **16:9**, embeds the same **HOUSE STYLE
block** (deep-navy theme, sky-blue accents, gold highlights — see [STYLE_PREAMBLE.md](STYLE_PREAMBLE.md)),
carries a small gold **"NN / 18" corner index**, and uses only **truthful** facts from the build
(5 agents, deterministic orchestrator, 3 surfaces + the twin, 19,389-term HPO grounding, 100% ACMG
accuracy, +12-pt detection, 6/6 mosaics, 13 event types, Haiku/Sonnet/Opus tiers, DGX Spark + RunPod).
Frames **10–12** are a darker, cinematic "render-showcase" act. For the best batch, generate one
anchor frame, then feed it back as a **style-reference image** for the rest.

## The 18 prompts, in order

| # | File | Family | What it shows |
| --- | --- | --- | --- |
| 01 | [01_summary_executive_onepager.md](01_summary_executive_onepager.md) | Summary | The whole engine on one page (problem · what · proof · where) — **start here** |
| 02 | [02_biology_mtor_mechanism.md](02_biology_mtor_mechanism.md) | Foundations | Why TSC happens — the mTOR-brake mechanism, organs, mosaicism |
| 03 | [03_story_problem_solution.md](03_story_problem_solution.md) | Story | Six clinical struggles → six engine responses ("why this exists") |
| 04 | [04_engine7_hero.md](04_engine7_hero.md) | Architecture | The engine, minimal one-glance title slide |
| 05 | [05_engine7_portrait.md](05_engine7_portrait.md) | Architecture | The engine, full "shows + explains everything" poster |
| 06 | [06_engine7_landscape.md](06_engine7_landscape.md) | Architecture | The engine, deck-ready left-to-right pipeline |
| 07 | [07_architecture_ai_tiering.md](07_architecture_ai_tiering.md) | Architecture | AI under the hood — deterministic core vs. LLM model tiers |
| 08 | [08_story_patient_journey.md](08_story_patient_journey.md) | Story | One patient end to end — how the data flows through the engine |
| 09 | [09_story_mosaic_recovery.md](09_story_mosaic_recovery.md) | Story | Ending the diagnostic odyssey — the mosaic-recovery story |
| 10 | [10_omniverse_four_scenes.md](10_omniverse_four_scenes.md) | Omniverse | The four digital-twin scenes: what each looks like + shows |
| 11 | [11_omniverse_pipeline.md](11_omniverse_pipeline.md) | Omniverse | The Spark→RunPod authoring/render pipeline |
| 12 | [12_omniverse_hero.md](12_omniverse_hero.md) | Omniverse | Cinematic hero: "the glass cloud IS the uncertainty" |
| 13 | [13_story_validation_scorecard.md](13_story_validation_scorecard.md) | Story | Measured results + the honest construct-validity caveat |
| 14 | [14_story_trust_layer.md](14_story_trust_layer.md) | Story | Why you can trust it — the five honesty disciplines |
| 15 | [15_engine7_demo_flow.md](15_engine7_demo_flow.md) | Demo | The 3-act live-demo process/steps (mirrors the runbook) |
| 16 | [16_story_cincinnati_roadmap.md](16_story_cincinnati_roadmap.md) | Story | Phase 1 → biobank study → Phase 2 (the 5 CCHMC touchpoints) |
| 17 | [17_story_generalization.md](17_story_generalization.md) | Story | One engine, many diseases ("swap the box labels") |
| 18 | [18_story_small_cheap_open.md](18_story_small_cheap_open.md) | Story | Small · Cheap · Open — the positioning thesis |

## Six families (for picking by purpose)

- **Foundations** — the disease biology: `02`
- **Architecture** — what the engine is: `04`, `05`, `06`, `07`
- **Demo flow** — the live demonstration: `15`
- **Omniverse** — the 3-D digital-twin visuals: `10`, `11`, `12`
- **Story / strategy** — why it matters, why trust it, where it goes: `03`, `08`, `09`, `13`, `14`, `16`, `17`, `18`
- **Summary** — the one-page at-a-glance: `01`

## How to use

- See **[STORYBOARD.md](STORYBOARD.md)** for deck recipes (which numbers to show, in what order,
  for which audience — Hagedorn/Cincinnati, the platform execs, keynote, families, researchers).
- See **`../TSC_ENGINE_MASTER_STUDY_GUIDE`** to study them as a learning path (each module pairs
  a numbered infographic with in-depth reading).
- **Keep the labels** — every prompt carries SYNTHETIC / decision-support / not-FDA-cleared
  framing, which is what makes the graphics credible to a clinical audience.
- If text renders crowded or misspelled (the hard case for any image model), regenerate with:
  *"keep the exact layout but regenerate only the text so every label is sharp and correctly
  spelled."*
