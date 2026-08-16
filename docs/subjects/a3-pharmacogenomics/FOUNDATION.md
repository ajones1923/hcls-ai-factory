# Pharmacogenomics — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Enzymes that metabolise drugs vary between people. If yours works faster or slower than average, a standard dose is the wrong dose.

## What this subject does with that idea

Explains why two patients given the same dose of the same drug can have completely different outcomes.

Genotype-guided dosing has CPIC guidelines behind it. For DPYD and fluorouracil the stakes are toxicity, not efficacy.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

CPIC guidance covers specific gene-drug pairs. Substrate relationships outside those pairs are pharmacology, not guideline-backed dosing advice.

Stating the limit is not a disclaimer bolted on the end — it is how you tell a tool that helps from
a tool that misleads.

## Vocabulary you will meet

| Term | Plain meaning |
|---|---|
| Capability | One named thing the platform can do, registered in `capabilities.json` |
| `live` / `planned` | Whether a capability actually answers today, or is intended |
| Decision support | Output that informs a clinician; it never decides |
| LIVE / REPRESENTATIVE / BURST | Whether a demo ran now, was pre-computed, or ran on remote GPUs |

## Try it

The demo for this subject is **`A3`**:

```bash
.venv/bin/python scripts/run_demo.py A3
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
