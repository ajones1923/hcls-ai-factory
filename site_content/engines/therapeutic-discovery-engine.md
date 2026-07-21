## In plain terms

Give this engine a biological **target**, and it designs and ranks candidate small molecules that
might act on it — a ten-stage *generate → dock → score* pipeline. Think of it as an automated,
in-silico medicinal-chemistry bench: it proposes new molecules, predicts how they'd fit the target,
and scores them for drug-likeness, then loops. The molecules it produces are **preclinical design
candidates** — starting points for real lab work, not drugs.

## Why it matters

The slowest, most expensive part of early drug discovery is going from "here is a druggable target"
to "here are a few promising molecules worth making." Compressing that search computationally is how
a single patient's target could, in principle, seed a therapeutic program in an afternoon instead of
over months.

## How it works

![Inside the Therapeutic Discovery Engine — generate, dock, score, reseed](../../assets/infographics/pages/therapeutic-discovery-engine-how.png)
/// caption
Target to ranked candidate molecules across a 10-stage pipeline. Preclinical design, not drugs. Illustrative.
///

1. **Generate** *(runs today)* — new candidate molecules are produced with **RDKit BRICS** (a
   rules-based method that recombines drug-like fragments into new molecules) — the real, live
   generator today.
2. **Dock** *(optional add-on)* — a molecule's fit into the target's binding pocket is predicted; this
   step uses the **DiffDock** model service, which must be running.
3. **Score** *(runs today)* — **RDKit** computes drug-likeness (**QED** — a 0-to-1 "how drug-like is
   it?" score), generates 3-D shapes (conformers), and ranks candidates
   across the ten stages.
4. **Reseed** — the best candidates feed a *generate-score-reseed* loop that iterates toward better
   molecules.

## What goes in, what comes out

- **In:** a **target**.
- **Out:** **ranked candidate molecules** — preclinical design candidates to evaluate.

## Where it fits

![Where the Therapeutic Discovery Engine sits — from druggable target to candidate leads](../../assets/infographics/pages/therapeutic-discovery-engine-fits.png)
/// caption
It consumes druggable targets surfaced upstream and feeds candidate leads to the disease programs. Illustrative.
///

Druggable targets identified by the [Precision Intelligence Engine](precision-intelligence-engine.md)
and the agents flow in; ranked candidate leads flow out toward the disease programs and molecular
tumor boards.

## Honest limits

- **Preclinical design candidates, not drugs.** Generated leads are hypotheses for the lab, nothing
  more.
- **What's real today.** The ten-stage chemistry (RDKit QED, conformers, ranking) is real and live,
  and generation uses the **RDKit BRICS** backend. The **BioNeMo / MolMIM** generation path is **not
  deployed**, and stage-7 docking requires the **DiffDock** model service.
- **Frontier co-folding is separate.** AlphaFold3-class complex co-folding (**Chai-1**) is a distinct,
  `planned` frontier capability — not part of this pipeline yet.
