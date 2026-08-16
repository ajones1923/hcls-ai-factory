# Precision Intelligence — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Retrieval-augmented generation: find the relevant evidence first, then reason over it, then cite it. The citation is what makes it checkable.

## What this subject does with that idea

The evidence layer. You ask a question in plain English and it retrieves what is known, with citations.

It is the difference between a variant list and an interpretation. It connects a letter change to published clinical significance.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Retrieval quality is bounded by what has been indexed. An unseeded collection returns nothing — that is a data problem, not a reasoning failure.

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

The demo for this subject is **`E2`**:

```bash
.venv/bin/python scripts/run_demo.py E2
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
