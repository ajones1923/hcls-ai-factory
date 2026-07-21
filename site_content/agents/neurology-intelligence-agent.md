## In plain terms

The Neurology Intelligence Agent is a broad neurology reasoning layer — decision support across **ten
domains**, from stroke triage to dementia evaluation to multiple sclerosis, with **Parkinson's**
featured. It grounds each answer in the literature and returns it, cited, for a clinician to weigh.

## Why it matters

Neurology spans very different problems that share a need for careful, staged reasoning: a stroke needs
minutes-fast triage, a neurodegenerative disease needs longitudinal staging. One grounded agent that
can reason across those domains — and cite why — is a genuine aid to a busy service.

*For a patient: one consistent, evidence-grounded read across very different neurological problems, from stroke to slow-moving disease.*

## How it works

![How the Neurology agent reasons — presentation and biomarkers across ten domains to a cited summary](../../assets/infographics/pages/neurology-intelligence-agent-how.png)
/// caption
Ten neurology domains, with Parkinson's S+N+G staging. Research/trial-use frameworks, decision support. Illustrative.
///

1. **Ingest** — the neurological presentation and any biomarkers, with patient context.
2. **Route across domains** — reasoning across ten domains (Parkinson's, Alzheimer's, dementia, ALS,
   MS, stroke, and more).
3. **Stage** — domain-specific staging such as **Parkinson's S+N+G**, **EDSS** for MS, and dementia
   evaluation.
4. **Ground the answer** — a cited neurology summary that refuses to fabricate where evidence is thin.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (presentation, biomarkers).
- **Out:** a grounded, cited neurology **summary**.

## Where it fits

![Where the Neurology agent sits — the neurology suite across the factory](../../assets/infographics/pages/neurology-intelligence-agent-fits.png)
/// caption
The neurology suite, Parkinson's featured; it draws on genomics and biomarker context. Illustrative.
///

It is the neurology suite of the factory — Parkinson's featured — drawing on genomics and biomarker
context and anchoring the neurology demonstration.

## Honest limits

- **Research/trial-use staging.** Frameworks such as Parkinson's **S+N+G** staging and related
  biomarker inputs are research- or trial-use, not routine clinical diagnostics.
- **Decision support, never diagnosis.** It supports a qualified clinician; it does not diagnose.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
