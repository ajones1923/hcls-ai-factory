# Single-Cell Analysis — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Cluster first, then name. Cells that behave alike group together; matching those groups to known marker genes tells you what they are.

## What this subject does with that idea

Takes a sample of thousands of individual cells and works out how many distinct cell types are present.

Bulk measurements average across cell types and hide the population that matters. Single-cell resolution recovers it.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Cluster annotation is by marker-gene overlap against a canonical panel. It is a strong heuristic, not a ground truth.

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

The demo for this subject is **`E8`**:

```bash
.venv/bin/python scripts/run_demo.py E8
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
