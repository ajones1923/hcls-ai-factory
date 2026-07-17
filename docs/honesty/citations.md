---
title: Citations & Evidence
description: Every clinical claim the HCLS AI Factory makes traces to a primary or authoritative source. If a claim can't be cited, it doesn't go on screen.
---

# Citations & Evidence

The rule is simple: **if a clinical claim can't be cited, it doesn't go on screen.** This page
collects the primary and authoritative sources behind the factory's clinical reasoning. Citations are
provided so an expert can verify each claim against the source — please confirm the exact reference
edition/version for any downstream use.

## Variant interpretation & genomics

- **ClinVar** — Landrum MJ, et al. *ClinVar: improving access to variant interpretations and
  supporting evidence.* Nucleic Acids Res. 2018. (clinical significance annotation)
- **ACMG/AMP variant classification** — Richards S, et al. *Standards and guidelines for the
  interpretation of sequence variants.* Genet Med. 2015;17(5):405–424.
- **ACMG secondary-findings gene list** — Miller DT, et al. *ACMG SF v3.x list of genes for reporting
  of secondary findings in clinical exome/genome sequencing.* Genet Med. 2023. (the E1 panel refresh)
- **AlphaMissense** — Cheng J, et al. *Accurate proteome-wide missense variant effect prediction with
  AlphaMissense.* Science. 2023;381(6664):eadg7492. (missense pathogenicity prior)
- **Human Phenotype Ontology (HPO)** — Köhler S, et al. *The Human Phenotype Ontology in 2021.*
  Nucleic Acids Res. 2021. (phenotype matching)
- **GIAB / HG002** — Zook JM, et al. Genome in a Bottle benchmark samples & truth sets. (the real
  data the [variant store](../factory/engines/genomics-engine.md) was **verified** against — Ts/Tv ≈ 2.0)

## Single-cell

- **pbmc3k (10x Genomics)** — the 3k-PBMC reference dataset, the canonical scanpy/Seurat tutorial
  data the [single-cell compute engine](../factory/engines/singlecell-compute.md) was **verified**
  against: 2,700 cells → the expected PBMC cell types (CD4 T, B, NK, CD14+ & FCGR3A+ monocytes,
  dendritic, megakaryocytes) with correct marker-gene evidence.

## Pharmacogenomics

- **CPIC** — Relling MV, Klein TE. *CPIC: Clinical Pharmacogenetics Implementation Consortium of the
  Pharmacogenomics Research Network.* Clin Pharmacol Ther. 2011; and the per-gene CPIC guidelines
  (cpicpgx.org). (the dosing/phenotype tables the Pharmacogenomics agent reasons over)

## Oncology

- **CIViC** — Griffith M, et al. *CIViC is a community knowledgebase for expert crowdsourcing the
  clinical interpretation of variants in cancer.* Nat Genet. 2017. (therapy-actionability evidence)

## Cardiology

- **ASCVD / Pooled Cohort Equations** — Goff DC Jr, et al. *2013 ACC/AHA Guideline on the Assessment
  of Cardiovascular Risk.* Circulation. 2014;129(25 Suppl 2):S49–S73.
- **Coronary artery calcium (Agatston score)** — Agatston AS, et al. *Quantification of coronary
  artery calcium using ultrafast computed tomography.* J Am Coll Cardiol. 1990;15(4):827–832.

## Tuberous Sclerosis Complex (flagship)

- **Everolimus in TSC** — Franz DN, et al. *EXIST-1* (SEGA; Lancet 2013); Bissler JJ, et al.
  *EXIST-2* (renal angiomyolipoma; Lancet 2013); French JA, et al. *EXIST-3* (TSC-associated
  seizures; Lancet 2016). Everolimus is **FDA-approved** in TSC (SEGA 2010; renal AML 2012;
  adjunctive treatment of TSC-associated partial-onset seizures 2018).
- **Gene therapy for TSC1/TSC2** — **preclinical.** No approved gene therapy exists; the factory is an
  open design/analysis bench for this direction, not a treatment. See the
  [Honesty & Governance ledger](index.md).

!!! note
    Clinical outputs are **decision support for a qualified clinician**, never autonomous diagnosis or
    prescribing. Sources are cited to let experts check the work — the factory's job is to make the
    evidence legible and traceable, not to replace the clinician who weighs it.
