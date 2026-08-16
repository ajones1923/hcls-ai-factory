# Clinical Imaging — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Structured reporting turns 'there is narrowing' into a category everyone interprets the same way, which is what makes it actionable.

## What this subject does with that idea

Reads medical scans and produces a graded, structured report a radiologist can act on.

Coronary CT to CAD-RADS 2.0 with modifiers, then a named genomic follow-up — imaging and genetics reasoning about the same patient.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

CAD-RADS is a reporting standard, not a diagnosis. The report supports a clinician's read; it does not replace it.

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

The demo for this subject is **`E4`**:

```bash
.venv/bin/python scripts/run_demo.py E4
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
