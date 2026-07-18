## In plain terms

The Clinical Trial Intelligence Agent matches a patient to the **right clinical trial** — and helps
design better trials. It reasons over a patient's profile and biomarkers to surface eligible trials,
grounded and cited, as decision support for the care team.

## Why it matters

The right trial can be the best available option, but eligibility criteria are dense and trials are
scattered. For a rare-disease or pediatric-oncology patient especially, a grounded matching layer can
turn a needle-in-a-haystack search into a short, explained shortlist.

## How it works

![How the Clinical Trial agent reasons — patient profile to matched, ranked trials](../../assets/infographics/pages/clinical-trial-intelligence-agent-how.png)
/// caption
Patient profile and biomarkers to a ranked, cited trial shortlist. Decision support. Illustrative.
///

1. **Ingest** — the patient profile and biomarkers.
2. **Match** — candidate trials retrieved and checked against eligibility.
3. **Optimize** — trial optimization, adaptive-design and biomarker-strategy reasoning.
4. **Ground the answer** — a ranked, cited trial shortlist that refuses to fabricate where evidence is
   thin.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (profile, biomarkers).
- **Out:** a grounded, cited set of **ranked trial matches**.

## Where it fits

![Where the Clinical Trial agent sits — the closing step for oncology and rare disease](../../assets/infographics/pages/clinical-trial-intelligence-agent-fits.png)
/// caption
It takes a patient or a discovery and finds the trial; the closing step of many workflows. Illustrative.
///

It is often the **closing step**: the [Precision Oncology Engine](../engines/precision-oncology-agent.md)
and [Rare Disease agent](rare-disease-diagnostic-agent.md) hand off to it to turn findings into trial
options.

## Honest limits

- **Decision support, never enrollment.** It surfaces and explains options for the care team; it does
  not enroll or decide.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
- **Eligibility must be confirmed.** Matches are a starting shortlist; formal eligibility is always
  confirmed with the trial.
