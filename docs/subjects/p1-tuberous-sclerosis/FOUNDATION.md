# Tuberous Sclerosis Complex — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

TSC1 and TSC2 act as a brake on mTOR, a growth pathway. Lose the brake and growth runs unchecked — which is why an mTOR inhibitor helps.

## What this subject does with that idea

The whole factory applied to one disease: tuberous sclerosis complex, a genetic condition causing growths across multiple organs.

Five disease-specific agents compose the horizontal engines — variant curation, trajectory, therapeutics, phenotype mapping and TAND surveillance.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Gene-therapy correction of TSC1/TSC2 is preclinical — an open design bench, not a cure today. Pediatric caution at full force. Decision support, never prescribing.

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

The demo for this subject is **`P1`**:

```bash
.venv/bin/python scripts/run_demo.py P1
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
