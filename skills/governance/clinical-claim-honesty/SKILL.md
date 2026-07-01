---
name: clinical-claim-honesty
description: >-
  The active honesty discipline applied to every clinical or scientific claim before it ships
  ANYWHERE — code docstrings, docs, demo scripts, UI copy, READMEs, outreach, slides. Consult
  whenever you write or review a sentence that asserts something about a disease, variant, drug,
  biomarker, therapy, or capability. It is the load-bearing wall of the mission: overstatement
  betrays the project more than any missing feature ever could.
---

# Clinical-Claim Honesty — the five checks every claim must pass

The power of this platform is that **what people see is real**. So before any clinical or
scientific claim leaves your hands — in a docstring, a doc, a demo line, a UI string, or an
outreach draft — run it through the five checks below. A claim that fails any one of them is not
"done"; it is cut or qualified until it passes. This is not review theater — the honesty register
in `lib/hcls_common/verify_gate.py` enforces the sharp edges deterministically, but most claims
live in prose the gate never sees. That prose is your responsibility.

## The five checks (run on every claim, everywhere)

1. **Decision support, never autonomy.** Is the output framed as *decision support for a qualified
   clinician* — never autonomous diagnosis, prescribing, or a treatment directive to a patient?
   "Flags a likely-pathogenic variant for clinician review" — yes. "Diagnoses," "confirms the
   diagnosis," "you should stop taking," "prescribe X" — no. (These exact shapes are `block`/`warn`
   patterns in `verify_gate.honesty_check`.)
2. **Maturity is labeled.** Is the claim *proven / preclinical / research-use / roadmap / gated* —
   and does the text **say which, out loud**? An unlabeled clinical claim reads as "shipping and
   validated." When in doubt, label down, not up. See the standing ledger below.
3. **Traceable to a source.** Is it traceable to a primary or authoritative source — **PMID / DOI /
   guideline** — and is that citation *present*? Clinical calculators, variant calls, biomarker
   thresholds, and dosing logic must carry their source with them. If a clinical claim can't be
   cited, it doesn't go on screen.
4. **Live ≠ mock (the registry honesty rule).** Does it honor the rule that a `live` capability is
   **never mock-served**? The registry rejects `status=live` with `serving=mock`
   (`capability_registry.py`). Never advertise a simulated or fallback result as real; label
   graceful fallbacks and representative data plainly ("representative output," "simulated for the
   demo").
5. **Elastic burst, not "all local."** If the claim touches a heavy or ARM-incompatible model
   (Chai-1/2, RFdiffusion, Evo 2), does it say **"elastic burst"** rather than imply everything runs
   on the one box? These burst to RunPod over a private Tailscale mesh — say so.

## The standing honesty ledger (state plainly wherever relevant)

Carry these forward verbatim in intent; do not soften a label upward:

- **Gene therapy for TSC1/TSC2 (and gene-correction generally): preclinical.** The factory is the
  open design/analysis bench, not a cure today.
- **MAISI synthetic imaging: research / augmentation / QA use — never a diagnostic source.**
- **Single-cell atlas similarity + foundation-model cell embeddings: roadmap.**
- **Chai-2 de novo binder / antibody design: gated / partnership access; integration contingent.**
- **α-synuclein SAA, plasma p-tau217, NSD-ISS / SynNeurGe staging: research- / trial-use biomarker
  inputs and frameworks the agents reason over — not routine clinical diagnostics.**
- **Pediatric RNA fusion detection (Arriba / STAR-Fusion): recommended near-term addition** to make
  the pediatric molecular tumor board airtight for a St. Jude / Cincinnati audience.
- **All clinical outputs: decision support only** — never autonomous diagnosis or prescribing.

## Where this applies (there is no exempt surface)

Docstrings and code comments · service READMEs and `docs/` · demo scripts and the words spoken over
a live demo · UI labels, tooltips, and result banners · report/PDF templates · commit messages ·
outreach, decks, social, and release notes. If a human or an LLM will read the sentence and take it
as a factual clinical claim, it is in scope.

## Do / Don't

**Do:** frame every output as clinician decision support; label maturity (preclinical / research /
roadmap / gated) out loud; cite the PMID/DOI/guideline next to the claim; say "elastic burst" when a
model leaves the box; label fallbacks and representative data as such; label *down* when unsure; run
`honesty_flags(text)` / `assert_publishable(text, llm=...)` on generated clinical text before it
ships.

**Don't:** write "diagnoses," "confirms/establishes the diagnosis," "cures," "FDA-approved/cleared,"
"100% / guaranteed," or "prescribe / you should stop taking"; present a preclinical, roadmap, or
gated item as if it ships; show a mock or fallback as a real `live` result; imply "all on one box"
when a model bursts to RunPod; put an uncited clinical claim on screen; round a maturity label up to
sound more impressive.

## Checklist (before any claim ships)

- [ ] Framed as decision support for a qualified clinician — no autonomy, no patient directive.
- [ ] Maturity labeled: proven / preclinical / research-use / roadmap / gated — and stated in the text.
- [ ] Cited to a primary/authoritative source (PMID / DOI / guideline) present alongside the claim.
- [ ] Honors live ≠ mock; any fallback or representative data is labeled, not passed off as real.
- [ ] Says "elastic burst" (RunPod over Tailscale) if a heavy/ARM-incompatible model is involved.
- [ ] Ledger-consistent (above); when unsure, labeled down, not up.
- [ ] If generated text: passed `verify_gate` (`is_publishable`) with no `block`, no unlabeled
      clinical content.

## Related
- `hcls-core-vision-mission` — why honesty is load-bearing (the north star it protects).
- `regulatory-compliance-posture` — the conservative regulatory/compliance framing this defends.
- `demo-foundation-alignment` — the honesty ledger and "earn a demo home" companion discipline.
- `11-security-and-secrets` — the governance/honesty gates (`api_gate`, `verify_gate`) that enforce this.
