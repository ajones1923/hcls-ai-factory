# Structural Biology — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Protein function follows shape. Cryo-electron microscopy and related methods let us see that shape and deposit it publicly in the Protein Data Bank.

## What this subject does with that idea

Works with the three-dimensional shape of proteins — searching, scoring and designing against them.

A drug has to physically fit its target. Knowing the shape, atom by atom, is what makes rational design possible.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

ESMFold prediction and ProteinMPNN design need CUDA, and these services are not currently served from the venv that has it, so structures shown are deposited PDB entries rather than predictions made live. The engine as a whole is registered `planned`: no process binds its aggregate port, and the model-level capabilities are the reachable surface.

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

The demo for this subject is **`E7`**:

```bash
.venv/bin/python scripts/run_demo.py E7
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
