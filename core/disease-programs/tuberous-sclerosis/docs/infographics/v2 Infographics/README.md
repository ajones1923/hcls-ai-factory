# TSC AI Engine — v2 storyboard (slides + narration)

The launch version. Each of the 18 frames is a **slide**, and each has two things in its file:

1. **🎙️ Narration script** — spoken, first-person, ready for an ElevenLabs voice clone. Written to
   *walk through everything*, grounded in `TSC_ENGINE_LEARN_GUIDE` (the plain-language story, with its
   analogies) and `TSC_ENGINE_OMNIVERSE_VISUALS_GUIDE` (the 3-D visuals). Combine the 18 in order and
   you have the full voiceover (also assembled in `FULL_NARRATION_SCRIPT.md`).
2. **🖼️ Nano Banana Pro prompt** — the updated, genericized image prompt (16:9, shared house style,
   `NN / 18` corner index).

## Changed from v1

- **Storytelling is now guide-based** — the narration carries the learn-guide framings and analogies
  (mosaicism = a *misprint in only some copies of a book*; the uncertainty envelope = a *hurricane
  forecast cone*; the orchestrator = an *air-traffic controller*; construct validity = a *flight
  simulator*; "swap the box labels, keep the wiring").
- **No institution or person named in the content.** Cincinnati Children's and the CHIO are *not*
  named anywhere in the v2 slides or scripts — a partner children's hospital, a tissue biobank, a
  clinical-informatics sponsor, and "a diagnostic-uncertainty NLP method" are used instead. (A
  Cincinnati-specific variant of the roadmap slide remains in the parent folder for the direct pitch.)

## How to produce the video

1. Generate each slide image from its Nano Banana Pro prompt (16:9). Use one anchor frame as a
   style-reference for the rest so the look stays locked.
2. Run each slide's narration script through ElevenLabs with your voice clone (one clip per slide, or
   the full script in one pass).
3. Assemble: slide image + its narration clip, in order 01 → 18. Total runtime ≈ 16–20 minutes.

## Script style notes (for ElevenLabs)

- First person, warm and confident, but **honest** — the construct-validity / synthetic / not-FDA
  caveat is said plainly where it belongs (slides 13 and the close).
- Plain spoken prose; punctuation drives the pacing. Bracketed `[beat]` marks an optional pause.
- Numbers are spoken naturally (e.g., "a hundred percent," "six out of six," "plus twelve points").
- Every figure is truthful to the build.

| # | Slide | ~secs |
| --- | --- | --- |
| 01 | The one-line story (open) | 70 |
| 02 | What TSC is — the biology | 70 |
| 03 | Why it's hard — six struggles | 65 |
| 04 | Meet the engine | 55 |
| 05 | The full architecture | 70 |
| 06 | How the data flows | 55 |
| 07 | How the AI is used responsibly | 60 |
| 08 | One patient, end to end | 60 |
| 09 | The flagship — mosaic recovery | 75 |
| 10 | The digital twin — four scenes | 65 |
| 11 | How the visuals are built | 55 |
| 12 | The honest-geometry hero | 55 |
| 13 | How we know it works | 75 |
| 14 | Why you can trust it | 60 |
| 15 | The live demo, in three acts | 65 |
| 16 | Where it goes | 65 |
| 17 | One engine, many diseases | 50 |
| 18 | Small, affordable, open (close) | 60 |

*SYNTHETIC demonstration data · decision support, clinician review required · not FDA-cleared · Apache 2.0.*
