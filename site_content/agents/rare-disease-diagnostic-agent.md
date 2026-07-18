## In plain terms

The Rare Disease Intelligence Agent helps narrow a **diagnostic odyssey**. It reasons from a patient's
phenotype (coded as **HPO** terms) and family (trio) variants toward a ranked shortlist of candidate
diagnoses, classifies variants by **ACMG** criteria, and does it all as grounded, cited decision
support.

## Why it matters

Rare-disease patients often wait years for a diagnosis, cycling through specialists. Systematically
matching phenotype to gene, classifying variants rigorously, and using trio (parent-child) data to
resolve inheritance is exactly the kind of patient, evidence-heavy reasoning an intelligence layer can
accelerate.

## How it works

![How the Rare Disease agent reasons — HPO phenotype and trio variants to ranked diagnoses](../../assets/infographics/pages/rare-disease-diagnostic-agent-how.png)
/// caption
Phenotype and trio variants to ACMG-classified candidate diagnoses. Decision support. Illustrative.
///

1. **Ingest** — the **HPO** phenotype and trio (parent-child) variants.
2. **Match phenotype** — HPO-to-gene matching against curated knowledge.
3. **Classify variants** — **ACMG** classification, with trio analysis to resolve inheritance.
4. **Ground the answer** — a ranked candidate-diagnosis shortlist (with gene-therapy tracking) that
   refuses to fabricate where evidence is thin.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (HPO phenotype, trio variants).
- **Out:** a grounded, cited set of **ranked candidate diagnoses**.

## Where it fits

![Where the Rare Disease agent sits — over germline genomics, feeding trials and disease programs](../../assets/infographics/pages/rare-disease-diagnostic-agent-fits.png)
/// caption
It reasons over the germline genomics substrate and hands off to trials and disease programs. Illustrative.
///

It reasons over the germline variants from the
[Genomics Foundations Engine](../engines/genomics-engine.md), hands off to the
[Clinical Trial agent](clinical-trial-intelligence-agent.md), and feeds the disease programs (such as
Tuberous Sclerosis).

## Honest limits

- **Decision support, never diagnosis.** It produces a *candidate* shortlist for a qualified clinician;
  it does not diagnose.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
- **Gene therapy is preclinical.** Where gene-therapy tracking points to a correction, that remains
  preclinical research, not a treatment available today.
