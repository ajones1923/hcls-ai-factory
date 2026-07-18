## In plain terms

The Precision Autoimmune Intelligence Agent reasons about autoimmune disease: it interprets
**autoantibodies**, analyzes **HLA** type, and reasons about **flare** risk — grounding each answer in
the literature as decision support for a clinician.

## Why it matters

Autoimmune conditions are heterogeneous and relapsing: the same diagnosis behaves differently across
patients, autoantibody patterns are subtle, and HLA associations carry real risk information. Bringing
those threads together into a grounded, cited picture helps a clinician see the whole patient.

## How it works

![How the Autoimmune agent reasons — autoantibodies and HLA to flare-risk reasoning](../../assets/infographics/pages/precision-autoimmune-agent-how.png)
/// caption
Autoantibody interpretation, HLA analysis, and flare reasoning. Decision support. Illustrative.
///

1. **Ingest** — autoantibody results and HLA typing, with patient context.
2. **Interpret autoantibodies** — pattern interpretation grounded in curated evidence.
3. **Analyze HLA** — HLA associations relevant to the condition and to drug risk.
4. **Reason about flares** — flare-risk reasoning, returned as a cited answer that refuses to fabricate.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (autoantibodies, HLA).
- **Out:** a grounded, cited autoimmune **answer** — interpretation and flare reasoning.

## Where it fits

![Where the Autoimmune agent sits — joining HLA and pharmacogenomics](../../assets/infographics/pages/precision-autoimmune-agent-fits.png)
/// caption
It joins HLA analysis with pharmacogenomics; anchor of the autoimmune program. Illustrative.
///

It draws on HLA analysis and hands off to the
[Pharmacogenomics agent](pharmacogenomics-intelligence-agent.md) for drug safety, anchoring the
autoimmune program.

## Honest limits

- **Decision support, never diagnosis.** It supports a qualified clinician; it does not diagnose.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
- **Flare reasoning is supportive.** Flare-risk reasoning is an aid to clinical judgment, not a
  prediction to act on unaccompanied.
