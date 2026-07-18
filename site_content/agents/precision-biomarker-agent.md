## In plain terms

The Precision Biomarker Intelligence Agent turns a patient's **multi-omics** and clinical data into a
grounded **risk profile**. It reasons over biological-age clocks and multi-domain risk, genotype-aware,
and explains its answer with citations — decision support for a clinician, not a verdict.

## Why it matters

Biomarkers are only useful in context: a value means one thing in one genotype and another in a
different one, and the signal that matters is usually spread across omics layers. Joining them into a
single, interpretable risk picture is what makes biomarkers actionable.

## How it works

![How the Biomarker agent reasons — multi-omics to biological-age clocks to genotype-aware risk](../../assets/infographics/pages/precision-biomarker-agent-how.png)
/// caption
Multi-omics joined into a genotype-aware risk profile. Research-use, decision support. Illustrative.
///

1. **Gather** — multi-omics and clinical data for the patient.
2. **Compute clocks** — biological-age estimates (**PhenoAge / GrimAge**-style) and **9-domain** risk.
3. **Join, genotype-aware** — a multi-omics join that keeps the patient's genotype in the loop.
4. **Ground the answer** — a cited risk profile; it refuses to fabricate where evidence is thin.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (multi-omics + clinical data).
- **Out:** a grounded, cited biomarker **risk profile**.

## Where it fits

![Where the Biomarker agent sits — the multi-omics join across the factory](../../assets/infographics/pages/precision-biomarker-agent-fits.png)
/// caption
It performs the multi-omics join that other capabilities draw on; part of the cardiometabolic story. Illustrative.
///

It provides the **multi-omics join** other capabilities lean on and is central to the
cardiometabolic / longevity demonstration.

## Honest limits

- **Research-use biomarkers.** The biological-age clocks and multi-domain risk are research-use
  frameworks, not routine clinical diagnostics.
- **Decision support, never diagnosis.** It informs a qualified clinician; it does not diagnose.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
