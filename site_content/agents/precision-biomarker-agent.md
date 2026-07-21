## In plain terms

A **biomarker** is a measurable biological signal — something in your blood, genes, or cells — that
hints at your health or your risk of disease. The Precision Biomarker Intelligence Agent pulls together
many such signals from different layers of a patient's biology — **multi-omics** (looking at DNA, RNA,
proteins, and more at once, rather than just one) — and turns them into one clear, explained **risk
profile** a clinician can act on. Because the same measurement can mean different things depending on a
person's genes, it keeps the patient's genetics in the loop — and it cites its sources rather than
handing down a verdict.

## Why it matters

A single lab value rarely tells the whole story: it can mean one thing in one person's genetic context
and something else in another's, and the signal that matters is usually spread across several layers of
biology at once. Joining them into one interpretable picture is what turns scattered measurements into
something a clinician can actually use.

*For a patient: a fuller, genetics-aware read of their biology — so risk is judged in context, not
from one number in isolation.*

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
