# Cardiology — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

A risk score is a probability, not a prophecy — and its value is in whether it changes what you do next.

## What this subject does with that idea

Turns cardiac data into a structured risk assessment with the reasoning shown.

Risk only matters if it changes management. Pair it with imaging findings and pharmacogenomics to decide whether a given statin will work.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Risk scores are population statistics applied to an individual. They inform, they do not determine.

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

The demo for this subject is **`E6`**:

```bash
.venv/bin/python scripts/run_demo.py E6
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
