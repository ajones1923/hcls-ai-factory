# Clinical Trial Intelligence — Foundation Learning Guide

For a reader with no background in this area. If you have watched the Foundations film series, this
follows directly from it.

## The one idea to hold

Every approved medicine passed through staged human trials. Matching a patient to an open trial is often the only route to a new option.

## What this subject does with that idea

Finds the trials a patient may be eligible for, and what a candidate molecule would need to enter one.

Trial matching closes the loop: a discovered molecule is only a hypothesis until there is a path to test it in people.

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

## What it cannot do

Matching is informational. Eligibility is determined by the trial investigators, not by this agent.

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

The demo for this subject is **`A6`**:

```bash
.venv/bin/python scripts/run_demo.py A6
```

If it reports `BLOCKED`, the message names exactly what is missing. That is intentional — a demo
that cannot run says so rather than showing you something canned.

## Next

[Advanced Learning Guide](ADVANCED.md) · [Demo Guide](DEMO.md)
