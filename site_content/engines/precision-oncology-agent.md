## In plain terms

The Precision Oncology Engine assembles a **molecular tumor board (MTB) packet**: it ranks candidate
therapies and matches clinical trials for a cancer patient, then packages the reasoning for a human
tumor board to review. It is tuned for **pediatrics** — fusion-first, and weighted for the decades a
child survivor will live.

## Why it matters

Childhood cancers are often driven by gene **fusions** rather than the point mutations adult panels
chase, and a pediatric survivor may live sixty years with the consequences of treatment. Ranking
therapies without weighting those **late effects** is the wrong answer. This engine is built around
both facts.

*For a patient: a tumor-board packet tuned to their exact cancer — and, for a child, weighted to protect their decades ahead.*

## How it works

![Inside the Precision Oncology Engine — fusions in, therapy ranking, trial matching, MTB packet out](../../assets/infographics/pages/precision-oncology-agent-how.png)
/// caption
A molecular tumor board packet, fusion-first and late-effects-weighted. Illustrative.
///

1. **Ingest** — tumor variants and **RNA fusions**, with the patient's context.
2. **Rank therapies** — candidate therapies scored with **late-effects weighting** for pediatric
   survivorship.
3. **Match trials** — trial matching against **Pediatric MATCH / COG**.
4. **Assemble** — an MTB packet, built with cross-modal joins to the imaging, biomarker, and trial
   capabilities.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (tumor molecular profile).
- **Out:** an **MTB packet** for a human tumor board.

## Where it fits

![Where the Precision Oncology Engine sits — composing imaging, biomarker, and trial capabilities](../../assets/infographics/pages/precision-oncology-agent-fits.png)
/// caption
It composes imaging, biomarker, and trial capabilities into one tumor-board view. Illustrative.
///

It sits over the horizontal engines and agents, pulling imaging reads, biomarker profiles, and trial
matches into a single packet — the heart of the pediatric molecular-tumor-board demonstration.

## Honest limits

- **A packet for a board, not a prescription.** The output supports a qualified molecular tumor board;
  it is never an autonomous treatment decision.
- **Decision support, always.** Therapy ranking and trial matching are aids to expert judgment, cited
  to their guidelines and sources.
