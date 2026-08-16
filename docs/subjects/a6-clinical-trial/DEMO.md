# Clinical Trial Intelligence — Demo Guide

**Demo key:** `A6` · full specification in
[`../../demos/DEMO_CATALOG.md`](../../demos/DEMO_CATALOG.md)

## Run it

```bash
.venv/bin/python scripts/run_demo.py A6
```

Check prerequisites without running anything:

```bash
.venv/bin/python scripts/run_demo.py --check-all
```

## What the audience is seeing

Finds the trials a patient may be eligible for, and what a candidate molecule would need to enter one.

Trial matching closes the loop: a discovered molecule is only a hypothesis until there is a path to test it in people.

## Say this out loud

> **Decision support for a qualified clinician — never autonomous diagnosis or prescribing.** Every output on this page is intended to inform a clinician's judgement, not replace it.

**Matching is informational. Eligibility is determined by the trial investigators, not by this agent.**

The label matters. **LIVE** means it ran, now, on real input, in front of the room. **REPRESENTATIVE**
means a pre-computed or curated result is standing in for a long or gated step — and you say so.
**BURST** means it ran live but on remote GPUs, which is "elastic burst", not "all on one box".

## If it will not run

The runner refuses to execute a demo declared LIVE whose service is unreachable, rather than quietly
degrading to a canned result. The failure message names the missing prerequisite:

- a service not listening → start it (`docker compose -f docker-compose.dgx-spark.yml up -d`)
- a missing Python package → install it
- a gated component → see [`../../build/GATED_SOFTWARE_BUILD_GUIDE.md`](../../build/GATED_SOFTWARE_BUILD_GUIDE.md)

## Transcript

Each run writes `demo/transcripts/A6.txt`. Diff it after a change — a silent regression
shows up there before anyone notices it on stage.
