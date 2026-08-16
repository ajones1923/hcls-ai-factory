# Genomic Foundation — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

A variant is a difference from a reference sequence. Most differences mean nothing; a few change a protein's shape, and a very few cause disease.

## What this subject does with that idea

Turns raw sequencing reads into a list of the specific letters where a person's genome differs from the reference.

Every downstream answer — which variant, which drug, which trial — depends on this list being right.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Alignment and variant calling require NVIDIA Parabricks, which is not installed on this box. Results shown today are pre-computed.

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

The demo for this subject is **`E1`**:

```bash
.venv/bin/python scripts/run_demo.py E1
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
