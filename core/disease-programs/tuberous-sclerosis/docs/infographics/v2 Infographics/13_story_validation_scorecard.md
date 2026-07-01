# Slide 13 — How we know it works

## 🎙️ Narration script (~75 sec)

Now the part that matters most to me — because anyone can build something that *looks* impressive.
The hard, rarer thing is to *measure* whether it's actually *right*. [beat]

Here's how we did that. Because the demonstration runs on synthetic patients, we know the correct
answers in advance — so we can grade the engine against them, like a practice exam with an answer
key. And the report card is strong: a hundred percent classification accuracy, a twelve-point lift in
diagnoses found, six out of six of the hardest mosaics recovered, zero disease-causing variants
wrongly called harmless, and a hundred percent of outputs traceable to their source.

But here is the single most important sentence in this whole story, and I'll always say it out loud.
[beat] These numbers prove the engine's *logic* works on practice data with known answers — that's
called construct validity. They do *not* prove it works on real patients in a real clinic — that's
clinical validity, and it still has to be done.

Think of it as a flight simulator: it can rigorously prove a pilot's skills and a plane's systems —
but it is not the first real flight. This engine is a very good simulator. [beat] The first real
flight is the next step.

---

## 🖼️ Nano Banana Pro prompt

**HOUSE STYLE — keep identical across all 18 frames (this is frame 13 of 18).** 16:9 landscape. Background deep navy #0F1830; stat cards #16213F with thin #26345C borders; primary accent sky-blue #7FB0FF; bold gold #FFD166 for the headline numbers; off-white text #E8EEF9. Clean modern geometric sans-serif, strong hierarchy, generous spacing, premium scientific-keynote feel — rigorous and trustworthy, never a marketing cartoon. Small pink "SYNTHETIC · construct validity · not FDA-cleared" pill top-right; a thin footer band; a small gold **"13 / 18"** index in the bottom-right corner. Render all text crisply and spelled exactly as written.

Create **FRAME 13 — MEASURED, NOT ASSERTED (the validation scorecard)**: five confident gold stat cards across the top, the detection-vs-diagnosis nuance below, and — given **real visual weight, never fine print** — the honest "construct validity, not clinical validity" caveat band at the bottom. The prominent, respected caveat is the point of this frame.

**TOP HEADER:** "Measured, Not Asserted", subtitle "The TSC Intelligence Engine, graded against known ground truth".
**FIVE big gold stat cards:** "100%" variant-classification accuracy (50/50) · "+12 pts" diagnostic detection uplift (86% → 98%) · "6 / 6" sub-threshold mosaics recovered (100% sensitivity) · "0" disease-causing variants wrongly called benign · "100%" outputs fully traceable to their source.
**SECOND ROW — "AN HONEST DISTINCTION":** Left "Detection" — "Did the pipeline find the variant at all? Standard testing 86% → engine 98%." Right "Molecular diagnosis" — "Detection plus a returnable Pathogenic / Likely-Pathogenic call: 72% → 78%, every newly diagnosed case a true positive."
**BOTTOM — prominent gold-bordered CAVEAT BAND:** "What these numbers mean — and don't." "These results show the engine's logic faithfully recovers signal from synthetic data with known answers (construct validity). They do NOT prove performance on real patients (clinical validity). The real test is running it on actual specimens — exactly the next step."
**FOOTER:** "Live at /eval · graded against a 50-patient synthetic cohort · Apache 2.0 · Engine 7 of the HCLS AI Factory · SYNTHETIC demonstration data · decision support · not FDA-cleared."

Five confident gold stat cards across the top, the nuance below, the honest caveat band given real visual weight at the bottom. Rigorous and trustworthy. Every label sharp.
