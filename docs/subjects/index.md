---
title: Subject Guides
description: Four guides and an infographic for each of the 8 engines, 8 intelligence agents and the TSC disease program.
---

# Subject Guides

Four guides plus an infographic for each of the 17 subjects — **85 assets**. Every port, test count
and capability status is generated from the repository by `scripts/subject_facts.py`, so these
pages cannot drift from the code without the generator noticing.

All clinical output described here is **decision support for a qualified clinician — never
autonomous diagnosis or prescribing.**

## Engines

| | Subject | Guides |
|---|---|---|
| E1 | [Genomic Foundation](e1-genomic-foundation/OVERVIEW.md) | [Foundation](e1-genomic-foundation/FOUNDATION.md) · [Advanced](e1-genomic-foundation/ADVANCED.md) · [Demo](e1-genomic-foundation/DEMO.md) |
| E2 | [Precision Intelligence](e2-precision-intelligence/OVERVIEW.md) | [Foundation](e2-precision-intelligence/FOUNDATION.md) · [Advanced](e2-precision-intelligence/ADVANCED.md) · [Demo](e2-precision-intelligence/DEMO.md) |
| E3 | [Therapeutic Discovery](e3-therapeutic-discovery/OVERVIEW.md) | [Foundation](e3-therapeutic-discovery/FOUNDATION.md) · [Advanced](e3-therapeutic-discovery/ADVANCED.md) · [Demo](e3-therapeutic-discovery/DEMO.md) |
| E4 | [Clinical Imaging](e4-clinical-imaging/OVERVIEW.md) | [Foundation](e4-clinical-imaging/FOUNDATION.md) · [Advanced](e4-clinical-imaging/ADVANCED.md) · [Demo](e4-clinical-imaging/DEMO.md) |
| E5 | [Precision Oncology](e5-precision-oncology/OVERVIEW.md) | [Foundation](e5-precision-oncology/FOUNDATION.md) · [Advanced](e5-precision-oncology/ADVANCED.md) · [Demo](e5-precision-oncology/DEMO.md) |
| E6 | [Cardiology](e6-cardiology/OVERVIEW.md) | [Foundation](e6-cardiology/FOUNDATION.md) · [Advanced](e6-cardiology/ADVANCED.md) · [Demo](e6-cardiology/DEMO.md) |
| E7 | [Structural Biology](e7-structural-biology/OVERVIEW.md) | [Foundation](e7-structural-biology/FOUNDATION.md) · [Advanced](e7-structural-biology/ADVANCED.md) · [Demo](e7-structural-biology/DEMO.md) |
| E8 | [Single-Cell Analysis](e8-single-cell/OVERVIEW.md) | [Foundation](e8-single-cell/FOUNDATION.md) · [Advanced](e8-single-cell/ADVANCED.md) · [Demo](e8-single-cell/DEMO.md) |

## Intelligence Agents

| | Subject | Guides |
|---|---|---|
| A1 | [CAR-T Intelligence](a1-cart/OVERVIEW.md) | [Foundation](a1-cart/FOUNDATION.md) · [Advanced](a1-cart/ADVANCED.md) · [Demo](a1-cart/DEMO.md) |
| A2 | [Precision Biomarker](a2-precision-biomarker/OVERVIEW.md) | [Foundation](a2-precision-biomarker/FOUNDATION.md) · [Advanced](a2-precision-biomarker/ADVANCED.md) · [Demo](a2-precision-biomarker/DEMO.md) |
| A3 | [Pharmacogenomics](a3-pharmacogenomics/OVERVIEW.md) | [Foundation](a3-pharmacogenomics/FOUNDATION.md) · [Advanced](a3-pharmacogenomics/ADVANCED.md) · [Demo](a3-pharmacogenomics/DEMO.md) |
| A4 | [Precision Autoimmune](a4-precision-autoimmune/OVERVIEW.md) | [Foundation](a4-precision-autoimmune/FOUNDATION.md) · [Advanced](a4-precision-autoimmune/ADVANCED.md) · [Demo](a4-precision-autoimmune/DEMO.md) |
| A5 | [Neurology Intelligence](a5-neurology/OVERVIEW.md) | [Foundation](a5-neurology/FOUNDATION.md) · [Advanced](a5-neurology/ADVANCED.md) · [Demo](a5-neurology/DEMO.md) |
| A6 | [Clinical Trial Intelligence](a6-clinical-trial/OVERVIEW.md) | [Foundation](a6-clinical-trial/FOUNDATION.md) · [Advanced](a6-clinical-trial/ADVANCED.md) · [Demo](a6-clinical-trial/DEMO.md) |
| A7 | [Rare Disease Diagnostic](a7-rare-disease-diagnostic/OVERVIEW.md) | [Foundation](a7-rare-disease-diagnostic/FOUNDATION.md) · [Advanced](a7-rare-disease-diagnostic/ADVANCED.md) · [Demo](a7-rare-disease-diagnostic/DEMO.md) |
| A8 | [Single-Cell Intelligence](a8-single-cell/OVERVIEW.md) | [Foundation](a8-single-cell/FOUNDATION.md) · [Advanced](a8-single-cell/ADVANCED.md) · [Demo](a8-single-cell/DEMO.md) |

## Disease Program

| | Subject | Guides |
|---|---|---|
| P1 | [Tuberous Sclerosis Complex](p1-tuberous-sclerosis/OVERVIEW.md) | [Foundation](p1-tuberous-sclerosis/FOUNDATION.md) · [Advanced](p1-tuberous-sclerosis/ADVANCED.md) · [Demo](p1-tuberous-sclerosis/DEMO.md) |

## Verify any number on these pages

```bash
.venv/bin/python scripts/subject_facts.py     # the source of every figure
.venv/bin/python scripts/run_all_tests.py     # 17 subjects, 8,402 tests
```
