# Nano Banana Pro prompt — TSC Engine OMNIVERSE: the Spark → RunPod render pipeline

Explains the authoring/render split — the DGX Spark authors the OpenUSD scenes (CPU, free),
RunPod renders them with RTX in Omniverse. Copy the block below into Nano Banana Pro.

> **Frame 11 / 18 · 16:9.** The Spark→RunPod pipeline frame (Omniverse act). See `STYLE_PREAMBLE.md`.

---

## Prompt

**HOUSE STYLE — keep identical across all 18 frames (this is frame 11 of 18).** 16:9 landscape. Background deep navy #0F1830; panels #16213F with thin #26345C borders; primary accent sky-blue #7FB0FF; gold #FFD166 highlights; off-white text #E8EEF9. Clean modern geometric sans-serif, strong hierarchy, generous spacing, premium NVIDIA / scientific-keynote feel — credible enterprise diagram, never a marketing cartoon. Small pink "SYNTHETIC · not FDA-cleared" pill top-right; a thin footer band; a small gold **"11 / 18"** index in the bottom-right corner. Render all text crisply and spelled exactly as written.

Create **FRAME 11 — FROM ENGINE RESULTS TO A FILM-QUALITY TWIN**: a clean four-stage **left-to-right pipeline** of the authoring/render split — the "DGX Spark" stages (results → author the scene, CPU / free) grouped on the left, the "RunPod · Omniverse RTX" render stage on the right ending in a cinematic mini-render of the twin; gold ribbons mark which machine does what; the honest principle is prominent at the bottom.

**TOP HEADER (full width):** Title **"From Engine Results to a Film-Quality Twin"**, subtitle **"Author on the DGX Spark · render on RunPod (NVIDIA Omniverse RTX)"**; top-right pink pill **"SYNTHETIC · not FDA-cleared"**.

**STAGE 1 (far left) — "ENGINE RESULTS"** (icon: data/cards), label "on the DGX Spark": *"validated, provenance-stamped outputs — lesion sizes, forecasts with 50/90% prediction intervals, thresholds, phenotype, recovered-mosaic flags."*

**Arrow →**

**STAGE 2 — "AUTHOR THE SCENE"** (icon: CPU chip + a small wireframe), a panel with a gold ribbon **"on the DGX Spark · CPU only · no GPU · free · deterministic"**: *"the engine writes an OpenUSD (.usda) scene — lesion radius and color as time-samples, the uncertainty envelope radii set equal to the prediction intervals, the threshold membrane, watermark + provenance as metadata."*

**Arrow → (label the arrow "copy the .usda files")**

**STAGE 3 — "RENDER WITH RTX"** (icon: GPU + glowing cube), a panel with a ribbon **"on RunPod · RTX GPU · NVIDIA Omniverse / USD Composer"**: *"open the scene, set RTX Path Tracing, apply MDL materials — OmniGlass for the uncertainty envelope, OmniPBR emissive for variant cells and recovery halos. Scrub the timeline −24 → +18 months."*

**Arrow →**

**STAGE 4 (far right) — "THE TWIN"** (a glowing mini-render: a teal sphere inside translucent blue glass shells near a red membrane): *"a film-quality, scrubbable anatomical twin — the lesion grows along its forecast inside a glass uncertainty cloud, crossing the threshold."*

**BOTTOM CALLOUT BAND (full width, gold-bordered):** a single bold honest principle: **"The render is exactly as confident as the engine — the glass envelope's size IS the forecast's 50% / 90% prediction interval."**

**FOOTER (thin):** *"OpenUSD · Apache 2.0 · Engine 7 of the HCLS AI Factory · authoring is CPU-side and dependency-free; only the RTX render needs a GPU · SYNTHETIC demonstration data · stylized atlas, not patient imaging · not FDA-cleared."*

Composition: a clean four-stage left-to-right flow, the "Spark" stages clearly grouped on the left and the "RunPod / Omniverse" render stage on the right, gold ribbons marking which machine does what, the honest principle prominent at the bottom. Credible enterprise pipeline look. Every label sharp and readable.

---

## Tips

- The point of this graphic is the **authoring/render split** — make Stages 1–2 visually read as
  "Spark / CPU / free" and Stage 3 as "RunPod / GPU / render." If it blurs, add: *"visually
  group the first two stages under a 'DGX Spark' bracket and the render stage under a 'RunPod'
  bracket."*
- Truthful throughout (CPU-side authoring, OmniGlass/OmniPBR, envelope = prediction interval).
  Keep the SYNTHETIC / not-FDA-cleared labels.
- Second pass for text: *"keep the exact pipeline layout but regenerate only the text so every
  label is sharp and correctly spelled."*
