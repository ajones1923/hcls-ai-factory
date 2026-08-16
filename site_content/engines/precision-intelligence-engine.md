## In plain terms

If the [Genomics Foundations Engine](genomics-engine.md) finds **where** a patient's variants are,
the Precision Intelligence Engine explains **what they mean**. It is the factory's interpretation
brain. It takes a raw variant file — millions of rows — and turns it into a short, cited, plain-language
narrative a clinician can actually use: which variants matter, why, and which point to a druggable
target. That interpretation is the shared substrate the eight intelligence agents reason over.

## Why it matters

A variant list is not an answer. The clinically meaningful signal is usually a handful of variants
buried in millions, and the value is in the *context* — is this variant known to be pathogenic, what
does it do to the protein, is there a therapy that targets it. This engine is what turns data into
**decision support**: grounded, sourced, and honest about uncertainty, never an autonomous diagnosis.

*For a patient: the difference between a raw list of millions of variants and a short, plain-language explanation of the few that actually affect their care.*

## How it works

![Inside the Precision Intelligence Engine — annotate, retrieve, reason, hand off](../../assets/infographics/pages/precision-intelligence-engine-how.png)
/// caption
From a variant file to a cited clinical narrative — grounded in retrieved evidence. Illustrative.
///

1. **Annotate** — each variant is tagged against curated knowledge: **ClinVar** clinical significance
   and **AlphaMissense** pathogenicity predictions.
2. **Retrieve** — retrieval-augmented generation (RAG) pulls the relevant evidence from a vector
   database, so the narrative is built on real sources rather than a model's memory.
3. **Reason** — a language model writes a grounded interpretation **with citations**, surfaces
   druggable targets, and is built to refuse rather than fabricate when evidence is thin.
4. **Hand off** — the interpreted result is what the eight specialist agents (pharmacogenomics,
   oncology, rare disease, and the rest) reason over. Chaining engines and agents into one governed
   run is the [workflow composer's](../../factory/platform.md) job, not this engine's.

## What goes in, what comes out

- **In:** a **VCF** of variants (from Engine 1) and a clinical **query**.
- **Out:** a **cited clinical narrative** — the interpreted, contextualized report.

## Where it fits

![Where the Precision Intelligence Engine sits — between genomics and the eight agents](../../assets/infographics/pages/precision-intelligence-engine-fits.png)
/// caption
The interpretation layer: it consumes the genomics substrate and produces what the eight agents
reason over. Illustrative.
///

It sits directly downstream of genomics and upstream of the agents — the interpretation layer the
rest of the clinical reasoning depends on.

## Honest limits

- **Decision support, not diagnosis.** Every output supports a qualified clinician; it never diagnoses
  or prescribes on its own.
- **Grounded, and it says when it can't be.** The interpretation is retrieval-grounded and cited; when
  its vector database or model key are absent it returns an honest degraded response, never invented
  clinical content.
- **VEP is optional and manual.** An Ensembl VEP container ships as an *optional, manually-run*
  annotator — it is **not** in the automated interpretation path.
