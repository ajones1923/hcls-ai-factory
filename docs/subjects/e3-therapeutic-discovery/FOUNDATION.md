# Therapeutic Discovery — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Most medicines work by locking onto one specific protein. If you know the protein's shape, you can search for molecules that fit it.

## What this subject does with that idea

Designs new candidate molecules shaped to fit a specific protein, and scores how well each might bind.

It attacks the hardest question in drug discovery — which molecule to even try — at the very front of a 10-15 year pipeline.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Molecule generation (MolMIM) and docking (DiffDock) are gated NVIDIA NIMs and are not installed. Candidates shown are pre-computed. This is the flagship demo and the easiest to overclaim.

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

The demo for this subject is **`E3`**:

```bash
.venv/bin/python scripts/run_demo.py E3
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
