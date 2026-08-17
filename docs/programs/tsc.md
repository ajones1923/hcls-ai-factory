---
title: Tuberous Sclerosis Complex
description: The flagship — the entire HCLS AI Factory converging on one child in a single afternoon, on one $4,699 box, open-source. Decision support, not diagnosis; gene therapy is preclinical.
---

# Flagship: Tuberous Sclerosis Complex

A child is born with a single change in one gene that can quietly seed growths across five organ
systems — brain, heart, kidneys, skin, lungs. Naming it can take a family *years* of scattered
specialist visits.

This is the flagship: **the entire HCLS AI Factory — genomics, imaging, clinical reasoning, and drug
discovery — converging on one child in a single afternoon, on one $4,699 computer, open for anyone to
run.** TSC is the ideal first patient story: multi-system enough to light up nearly every engine and
agent, with a clean, honest gene-to-drug mechanism. Decision support for a clinician — never a
diagnosis on its own.

<video class="cap-video" controls preload="metadata" playsinline poster="/assets/videos/posters/tsc-usecase.jpg" src="/assets/videos/tsc-usecase.mp4">
  Your browser can't play embedded video — <a href="/assets/videos/tsc-usecase.mp4">download the video</a>.
</video>
/// caption
The flagship story in under a minute.
///

## What TSC is, in plain terms

Tuberous Sclerosis Complex is a genetic condition — roughly **1 in 6,000 births** — in which benign
tumors (**hamartomas**) grow across many organs. Because it touches so many systems and follows a
clean genetic mechanism, it is the ideal first patient story for the factory: enough breadth to
exercise nearly every capability, with a therapy that is real today.

## One condition, the whole factory

![One condition, the whole factory — TSC across brain, heart, kidney, skin, and lungs](../assets/infographics/pages/tsc-multisystem.png)
/// caption
Illustrative.
///

TSC rarely stays in one place. It causes growths and problems in the **brain** (seizures, and effects
on development and behavior), the **heart**, the **kidneys**, the **skin**, and the **lungs** — five
organ systems at once. *(The medical names a specialist would use: subependymal giant cell astrocytoma
and TAND in the brain, cardiac rhabdomyoma, renal angiomyolipoma, skin angiofibromas, and pulmonary
LAM.)* That whole-body reach is exactly why one TSC patient needs the genomics, imaging, cardiology,
pharmacogenomics, rare-disease, and neurology capabilities together — one child exercises almost the
entire platform.

## The gene-to-drug story

![Tuberous Sclerosis Complex — mechanism to medicine](../assets/infographics/tsc-flagship.png)
/// caption
Illustrative.
///

**TSC1 / TSC2 loss → mTORC1 hyperactivation → everolimus (an mTOR inhibitor).** The TSC1 (hamartin)
and TSC2 (tuberin) proteins normally form a complex that restrains **mTORC1**, a master regulator of
cell growth. When a variant knocks that complex out, mTORC1 runs unchecked and hamartomas grow — so an
mTOR inhibitor addresses the mechanism directly. **Everolimus** is real and FDA-approved in TSC (SEGA
2010, renal angiomyolipoma 2012, seizures 2018). The factory helps a clinician reason over a patient's
variants, imaging, and phenotype to *support* that decision. It does not make it.

## The full story — weight, compression, hope

![The arc of the flagship — weight, compression, hope](../assets/infographics/pages/tsc-arc.png)
/// caption
Illustrative.
///

The program is built as three beats: the **weight** of a family facing a multi-system diagnosis; the
**compression** of months of cross-specialty work into a single, governed afternoon; and honest
**hope** — a real, approved therapy today, and an open bench for what comes next.

Watch the whole journey — the diagnostic odyssey, and how the factory changes it:

<video class="cap-video" controls preload="metadata" playsinline poster="/assets/videos/posters/tsc-two-journeys.jpg" src="/assets/videos/tsc-two-journeys.mp4">
  Your browser can't play embedded video — <a href="/assets/videos/tsc-two-journeys.mp4">download the deep-dive film</a>.
</video>
/// caption
Deep dive (~5 min): "One child, two journeys" — the diagnostic odyssey, the hamartin–tuberin / mTORC1 mechanism, mechanism-matched therapy, and designing a new precision molecule, shown with a real mTOR-structure animation. Decision support; molecule design is preclinical (illustrative structure).
///

## How the factory composes here

A disease program is a **vertical** — a thin, deterministic orchestrator that composes the horizontal
[engines](../factory/engines/index.md) and [intelligence agents](../factory/agents/index.md) for one
condition. For a single TSC child, that looks like:

1. **Genomics.** The Genomic Foundation engine calls variants from sequence and looks for the somatic
   **mosaicism** that standard germline pipelines miss — a genuine TSC diagnostic trap.
2. **Interpretation.** [Claude](../factory/frontier-models/anthropic.md), the reasoning layer, weighs
   the variant against ACMG criteria and ClinVar and ties the scattered findings back to one
   mechanism; the **Rare Disease** and **Pharmacogenomics** agents add HPO phenotype matching and
   mTOR-inhibitor dosing context.
3. **Structure & discovery.** For the drug side, the Therapeutic Discovery engine generates and docks
   candidate molecules against the mTOR pathway — **preclinical design, not a therapy**. Co-folding
   the top candidate with **Chai-1** as an independent structural check is *planned*, not live; when
   it stands up it will run by [elastic burst](../run/hardware.md) to a remote GPU.
4. **Surveillance.** The [Tuberous Sclerosis engine](../factory/engines/tuberous-sclerosis-engine.md)'s
   disease-specific agents — variant curator, phenome mapper, trajectory modeler, TAND surveillance,
   therapeutics strategist — turn the results into reviewable, clinician-facing surfaces.

Every step is composable, inspectable, and honestly labeled — and the whole thing runs on one box,
bursting to a remote GPU only for the heaviest models.

## Where to go next

- **Clinicians & geneticists** — see the [demonstrations](../demos/index.md) and the
  [honesty & governance](../honesty/index.md) posture.
- **Builders** — [run the whole factory yourself](../run/index.md) on one box, under Apache-2.0.
- **Researchers** — explore the [engines](../factory/engines/index.md) and
  [intelligence agents](../factory/agents/index.md) the program composes.
- **Families & advocates** — start with the [mission](../about.md) behind the project.

!!! warning "Honesty ledger — carried at full force"
    - **Gene therapy for TSC1/TSC2 is preclinical** — an open design/analysis bench, not a treatment
      available today.
    - **Everolimus is real and approved** — the factory supports, it does not prescribe.
    - **Molecule design is preclinical**, and **Chai-1 is planned** — not yet live.
    - **Decision support, not diagnosis**, and **pediatric / vulnerable-population caution at full
      force.** See [Honesty & Governance](../honesty/index.md).

*Replication roadmap: NF1, NF2, Rett, Williams, and the broader mTORopathies — same pattern, same
horizontal foundation.*
