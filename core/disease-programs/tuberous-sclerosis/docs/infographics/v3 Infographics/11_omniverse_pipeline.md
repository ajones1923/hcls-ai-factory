# Slide 11 — How the visuals are built

## 🎙️ Narration script (~70 sec)

You might assume rendering film-quality 3-D needs a data center. It doesn't — and the trick is a
clean split. [beat]

*Building* the scene is just turning the engine's results into a structured 3-D format called
OpenUSD. That's pure data work — it runs on the little DGX Spark, on the CPU, for free, no graphics
card needed. The lesion's size and color become an animation over time; the uncertainty envelope's
radii are set *equal* to the prediction intervals; the threshold and the provenance ride along as
metadata. And notice — the twin never invents a number. It's a *lens*, not a source: it only shows
what the engine already computed and checked.

Then the heavy, beautiful part — the RTX path tracing, the glass, the glow — happens separately, on a
rented cloud GPU, in Omniverse. In short: the DGX Spark does the thinking and building; the cloud GPU does the beautiful picture.

And because the envelope's size is literally the prediction interval, the render is exactly as
confident as the engine. [beat] It can't lie. Down the road, that same twin could sit on a clinic
table in augmented reality — and one day be built from a child's own scan.

---

## 🖼️ Nano Banana Pro prompt

HOUSE STYLE — keep identical across all 18 frames (this is frame 11 of 18). 16:9 landscape. Background deep navy #0F1830; panels #16213F with thin #26345C borders; primary accent sky-blue #7FB0FF; gold #FFD166 highlights; off-white text #E8EEF9. Clean modern geometric sans-serif, strong hierarchy, generous spacing, premium scientific-keynote feel — credible enterprise diagram, never a marketing cartoon. Small pink "SYNTHETIC · not FDA-cleared" pill top-right; a thin footer band; a small gold "11 / 18" index in the bottom-right corner. Render all text crisply and spelled exactly as written. (Everything in this HOUSE STYLE note is design direction — do NOT print any of its words on the image, including "premium scientific-keynote", "frame", "frame index", "house style", or "index"; the only thing it adds to the canvas is the small gold "11 / 18" number in the bottom-right corner. Structural words below such as "footer", "callout band" are layout directions, not text to draw — render on the image only the words shown inside quotation marks. The stage labels "STAGE 1"…"STAGE 4" ARE meant to appear.)

Create FRAME 11 — FROM ENGINE RESULTS TO A FILM-QUALITY TWIN: a clean four-stage left-to-right pipeline of the authoring/render split — the "DGX Spark" stages (results → author the scene, CPU / free) grouped on the left; the "cloud GPU · Omniverse RTX" render stage on the right ending in a cinematic mini-render of the twin; gold ribbons mark which machine does what; the honest principle prominent at the bottom.

Across the top, a header — "From Engine Results to a Film-Quality Twin", subtitle "Author on the DGX Spark · render on a cloud GPU (Omniverse RTX)".
STAGE 1 (far left) "ENGINE RESULTS" (data icon), label "on the DGX Spark": "validated, provenance-stamped outputs — lesion sizes, forecasts with 50/90% prediction intervals, thresholds, phenotype, recovered-mosaic flags."
Arrow → STAGE 2 "AUTHOR THE SCENE" (CPU chip + wireframe), gold ribbon "on the DGX Spark · CPU only · no GPU · free · deterministic": "the engine writes an OpenUSD (.usda) scene — lesion radius and color as time-samples, the uncertainty envelope radii set equal to the prediction intervals, the threshold membrane, watermark + provenance as metadata."
Arrow (labeled "copy the .usda files") → STAGE 3 "RENDER WITH RTX" (GPU + glowing cube), ribbon "on a cloud GPU · RTX · Omniverse / USD Composer": "open the scene, set RTX Path Tracing, apply MDL materials — OmniGlass for the uncertainty envelope, OmniPBR emissive for variant cells and recovery halos. Scrub the timeline −24 → +18 months."
Arrow → STAGE 4 (far right) "THE TWIN" (a glowing mini-render: teal sphere inside blue glass shells near a red membrane): "a film-quality, scrubbable anatomical twin — the lesion grows along its forecast inside a glass uncertainty cloud, crossing the threshold."
Along the bottom, a gold-bordered callout band reading: "The render is exactly as confident as the engine — the glass envelope's size IS the forecast's 50% / 90% prediction interval."
Along the very bottom, a thin footer strip reading: "OpenUSD · Apache 2.0 · Engine 7 of the HCLS AI Factory · authoring is CPU-side and dependency-free; only the RTX render needs a GPU · SYNTHETIC demonstration data · stylized atlas, not patient imaging · not FDA-cleared."

A clean four-stage flow, Spark stages grouped left and the cloud-GPU / Omniverse render right, the honest principle prominent at the bottom. Every label sharp.
