---
name: persona-domain-experts
description: >-
  Deep drill-down for reaching domain experts — clinicians, geneticists, oncologists, and pharma
  R&D. Consult when an outreach artifact primarily targets specialists: technical talks, methods
  sections, the "technical cut" of a video, expert-facing docs, or anything a working scientist
  will scrutinize. Pairs with broad-general-persona and the clinical honesty discipline.
---

# Domain Experts — the respect that decides the project's fate

Clinicians, geneticists, oncologists, and pharma R&D scientists — the people whose respect the
project needs and whose public pushback carries real weight. They read the seams, not the sizzle.
They convert from skeptic to advocate **on honesty alone**, and losing one is expensive and hard to
reverse: a single credible expert calling it overhyped travels further than any launch post.

## What wins them / What loses them

**Wins them:** precision, real tools shown actually running, and limits stated plainly *before they
have to ask*. They trust the person who labels their own caveats — honest "preclinical" vs.
"proven," named tools (Parabricks, DeepVariant, MolMIM, DiffDock, ClinVar, AlphaMissense), real
provenance, reproducibility they can check. Understatement reads as competence.

**Loses them / the risk:** any overclaim, and especially any blur between **decision support and
diagnosis**. Spectacle outrunning evidence is disqualifying — one unlabeled roadmap item shown as
if it ships, one "diagnoses" where you meant "flags for a clinician," and they write the whole
thing off. Their wince is the leading indicator that the artifact went too far.

## How to speak to them

- Use the **technical cut**: show tools running, not just polished results. "You can run this
  yourself, end to end, on one box" is the claim they cannot dismiss.
- **Label the honesty status out loud**: live vs. planned vs. preclinical vs. research-use vs.
  roadmap. Name the gated/burst pieces honestly (heavy models burst to remote GPUs — say "elastic
  burst," never imply all-local).
- Say **"decision support for a qualified clinician — never autonomous diagnosis or prescribing"**
  and mean it in every framing. This is the line they are watching for.
- **Cite clinical claims** to primary or authoritative sources; if it can't be cited, keep it off
  the slide. Distinguish what the platform *proves* (structure-validated candidates) from what it
  *supports* (a clinician's reasoning).
- Frame openness as an invitation to *audit*: fork it, run it, find the cracks, improve it. Experts
  respect a project that hands them the tools to check its own claims.

## Do / Don't

**Do:** state limits before they ask; show it running; cite sources; keep decision-support vs.
diagnosis crisp; welcome scrutiny. **Don't:** dramatize past the evidence; call anything a cure or
a diagnosis; hide the preclinical/roadmap label; imply everything runs locally when it bursts.

## Related
- `broad-general-persona` — the five-persona map and the technical-cut principle.
- `clinical-claim-honesty` — the decision-support-not-diagnosis discipline in depth.
- `demo-foundation-alignment` — the honesty ledger and citation rule to check every claim against.
