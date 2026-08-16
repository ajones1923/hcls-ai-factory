# CAR-T Intelligence — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

A chimeric antigen receptor is a designed protein that teaches a patient's own T-cells to recognise a target on cancer cells.

## What this subject does with that idea

Reasons about engineered T-cell therapies — which construct, for which patient, with which risks.

CAR-T is powerful and dangerous in the same breath. Construct choice and its safety counterweight belong in the same conversation.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

On-target off-tumour toxicity and cytokine release are real. Any construct suggestion carries its monitoring requirement.

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

The demo for this subject is **`A1`**:

```bash
.venv/bin/python scripts/run_demo.py A1
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
