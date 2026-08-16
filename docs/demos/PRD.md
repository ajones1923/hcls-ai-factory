# PRD — The 17 Demonstrations

**Date:** 2026-08-15 · **Companion:** [`DEMO_CATALOG.md`](DEMO_CATALOG.md) (the specs) ·
[`BUILD_GUIDE.md`](BUILD_GUIDE.md) (how to build and test them)

---

## 1. Problem

The platform has 8,402 passing tests and **nothing a visitor can watch**. There is no runnable demo
for any of the 17 subjects. Fifteen request payloads exist in `demo/requests/`, but no harness runs
them, no expected output is recorded, and nothing states what is live versus pre-computed.

## 2. Goals

| # | Goal | Done when |
|---|---|---|
| DM-1 | One demo per subject | 17 demos exist and each executes start to finish |
| DM-2 | Every demo honestly labelled | each is LIVE / REPRESENTATIVE / BURST, matching what actually ran |
| DM-3 | Reproducible | `scripts/run_demo.py <id>` produces the same output twice |
| DM-4 | Legible to a non-specialist | a clinician follows it without the code |
| DM-5 | Inline honesty | decision-support limit stated the first time clinical output appears |
| DM-6 | Recorded evidence | each demo writes a transcript that can be diffed |

## 3. Non-goals

- Demos that fake a gated capability. A pre-computed MolMIM result is **REPRESENTATIVE**, said out
  loud — never dressed as live generation.
- One demo covering everything. Coverage is a property of the portfolio, not of any single demo.
- Clinical validity claims. Everything remains decision support.

## 4. Requirements

**DR-1 · A single runner.** `scripts/run_demo.py <id>` resolves the subject, checks prerequisites,
executes, prints a transcript, and exits non-zero if a prerequisite is missing. It must **refuse to
run a LIVE demo whose service is unreachable** rather than silently degrading.

**DR-2 · Prerequisite declaration.** Each demo declares what it needs: a service endpoint, a seeded
Milvus collection, an API key, a Python package, a gated component. The runner checks before
starting.

**DR-3 · Label enforcement.** The label is computed, not asserted. If a demo declares LIVE and its
endpoint fails the health probe, the runner fails. This is the same rule that caught
`singlecell-compute` and `structural-biology-engine` registered `live` while unreachable.

**DR-4 · Recorded transcripts.** `demo/transcripts/<id>.txt`, diffable, so a regression is visible.

**DR-5 · Ordering for a room.** The catalogue defines a session arc: D17 (weight) → D4 (proof) →
D8 (something entirely real) → the audience-specific agent → D14 (close the loop).

## 5. Dependencies

| Demo | Blocked by | Unblocks after |
|---|---|---|
| D1 Genomic Foundation | Parabricks (G2) | full FASTQ→VCF live |
| D3 Therapeutic Discovery | MolMIM + DiffDock (G3/G4) | live or burst generation |
| D7 Structural Biology | CUDA torch + ESMFold (G1/G5) | live prediction |
| D8 Single-Cell | `pip install scanpy anndata` — **ungated** | immediately |
| D2, D5, D6, D9–D16 | seeded Milvus + `ANTHROPIC_API_KEY` | immediately |

**Eleven of seventeen are unblocked today** by seeding Milvus and setting one API key. That is the
cheapest path to a demonstrable platform and should happen before any gated work.

## 6. Risks

| Risk | Mitigation |
|---|---|
| A demo claims LIVE and degrades quietly | DR-3 — the runner fails instead |
| Pre-computed results drift from the code | DR-6 transcripts are diffed in CI |
| D3 overclaims because it is the flagship | Catalogue names it explicitly as the one to watch |
| NIM images are x86-only | Relabel BURST and change the narration to "elastic burst" |
| Empty Milvus makes agents look broken | Seeding is a prerequisite, checked by DR-2 |

## 7. Sequence

```
seed Milvus + ANTHROPIC_API_KEY ──► 11 demos LIVE
        │
pip install scanpy anndata ──────► D8 LIVE
        │
gated install (G1-G5) ───────────► D1 D3 D7 LIVE/BURST ──► all 17
```

## 8. Open questions

1. **Is a seeded Milvus corpus available**, or must collections be rebuilt from source data?
2. **Should demos run against containers or host processes?** Unresolved platform-wide.
3. **Recording:** are these demoed live, or captured as video? If captured, they inherit the film
   pipeline's honesty rules.
