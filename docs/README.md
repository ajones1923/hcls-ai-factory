# HCLS AI Factory — Documentation

Start with the root [README](../README.md) and [CONTRIBUTING](../CONTRIBUTING.md).
This index covers the platform-wide documentation.

**New here?** The fastest orientation is one subject's four guides — try
[A3 · Pharmacogenomics](subjects/a3-pharmacogenomics/OVERVIEW.md), the most rigorously verified
content in the platform.

---

## Per-subject guides — 8 engines · 8 agents · 1 disease program

Four guides plus an infographic for each of the 17 subjects (85 assets). Every port, test count and
capability status in them is generated from the repo by `scripts/subject_facts.py`, so they cannot
drift from the code without the generator noticing.

| | Engines | | Agents |
|---|---|---|---|
| E1 | [Genomic Foundation](subjects/e1-genomic-foundation/OVERVIEW.md) | A1 | [CAR-T Intelligence](subjects/a1-cart/OVERVIEW.md) |
| E2 | [Precision Intelligence](subjects/e2-precision-intelligence/OVERVIEW.md) | A2 | [Precision Biomarker](subjects/a2-precision-biomarker/OVERVIEW.md) |
| E3 | [Therapeutic Discovery](subjects/e3-therapeutic-discovery/OVERVIEW.md) | A3 | [Pharmacogenomics](subjects/a3-pharmacogenomics/OVERVIEW.md) |
| E4 | [Clinical Imaging](subjects/e4-clinical-imaging/OVERVIEW.md) | A4 | [Precision Autoimmune](subjects/a4-precision-autoimmune/OVERVIEW.md) |
| E5 | [Precision Oncology](subjects/e5-precision-oncology/OVERVIEW.md) | A5 | [Neurology](subjects/a5-neurology/OVERVIEW.md) |
| E6 | [Cardiology](subjects/e6-cardiology/OVERVIEW.md) | A6 | [Clinical Trial](subjects/a6-clinical-trial/OVERVIEW.md) |
| E7 | [Structural Biology](subjects/e7-structural-biology/OVERVIEW.md) | A7 | [Rare Disease Diagnostic](subjects/a7-rare-disease-diagnostic/OVERVIEW.md) |
| E8 | [Single-Cell Analysis](subjects/e8-single-cell/OVERVIEW.md) | A8 | [Single-Cell Intelligence](subjects/a8-single-cell/OVERVIEW.md) |

**P1 · [Tuberous Sclerosis Complex](subjects/p1-tuberous-sclerosis/OVERVIEW.md)** — the disease
program, composed from the horizontal engines and agents. Not a ninth engine.

Each subject has: **Overview** · **Foundation Learning** · **Advanced Learning** · **Demo**.

## Paper

- [**Honesty by Construction**](paper/index.md) — preprint describing the platform and the
  capability registry that makes its claims machine-checkable

## Build, audit and operate

Produced 2026-08-15. Every figure was measured by running the code.

- [**Build documentation index**](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/README.md) — start here
- [Gap Analysis](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/GAP_ANALYSIS.md) — what exists vs what is needed
- [PRD](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PRD.md) · [Build Guide](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/BUILD_GUIDE.md) — cold clone to running platform
- [Port Map](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md) — the UI / UI+1 convention, CI-enforced
- [Best-Practices & Accuracy Audit](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/BEST_PRACTICES_AUDIT.md) — security, design, clinical facts
- [Layout Review](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/LAYOUT_REVIEW.md) · [Scorecard](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/SCORECARD.md) · [Resource Report](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/RESOURCE_REPORT.md)
- [Gated Software PRD](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/GATED_SOFTWARE_PRD.md) · [Build Guide](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/GATED_SOFTWARE_BUILD_GUIDE.md) — also as [.docx](https://github.com/ajones1923/hcls-ai-factory/tree/main/docs/build/docx)

## Demonstrations

Two axes, deliberately distinct — see the reconciliation note in the catalogue.

- [**Portfolio D1–D7**](demos/index.md) — patient-story demonstrations; the coverage contract
- [**Catalogue E1–E8 / A1–A8 / P1**](demos/DEMO_CATALOG.md) — one per subject
- [Demo PRD](demos/PRD.md) · [Demo Build Guide](demos/BUILD_GUIDE.md)

```bash
.venv/bin/python scripts/run_demo.py --check-all   # what is ready
.venv/bin/python scripts/run_demo.py E8            # runs for real today
```

## Architecture & structure

- [STRUCTURE.md](STRUCTURE.md) — repository layout
- [HCLS_AI_FACTORY_PROJECT_BIBLE.md](HCLS_AI_FACTORY_PROJECT_BIBLE.md) — comprehensive reference
- [HCLS_AI_FACTORY_ARCHITECTURE_DIAGRAM_RESEARCH.md](HCLS_AI_FACTORY_ARCHITECTURE_DIAGRAM_RESEARCH.md)
- [HCLS_AI_FACTORY_MINDMAP.md](HCLS_AI_FACTORY_MINDMAP.md) — the platform at a glance

## Deploy & operate

- [DGX Spark Deployment Guide](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md)
- [RUNBOOK.md](RUNBOOK.md)
- [Build Guide](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/BUILD_GUIDE.md) — the current, measured bring-up path

## Honesty & governance

- [Honesty framework](honesty/index.md) · [Capability maturity matrix](honesty/maturity-matrix.md) *(generated at build time by `scripts/site_gen_matrix.py`)* · [Honesty ledger](honesty/ledger.md) · [Decision-support posture](honesty/decision-support.md) · [Citations](honesty/citations.md) · [Governance](honesty/governance.md)

Three rules the documentation itself is held to:

1. A capability marked `live` must answer a health probe. Two were found registered `live` with
   nothing bound to their ports; both are now resolved.
2. A demo labelled **LIVE** must have run, now, on real input. `run_demo.py` fails rather than
   degrading to a canned result.
3. All clinical output is **decision support for a qualified clinician** — never diagnosis or
   prescribing.

## Learning path (platform-wide)

- [Foundations](HCLS_AI_FACTORY_LEARNING_GUIDE_FOUNDATIONS.md) · [Foundations (unified)](HCLS_AI_FACTORY_LEARNING_GUIDE_FOUNDATIONS_UNIFIED.md) · [Advanced](HCLS_AI_FACTORY_LEARNING_GUIDE_ADVANCED.md)

## Reference & background

- [White paper](HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md) · [Branding research](HCLS_AI_FACTORY_UNIFIED_BRANDING_RESEARCH_GUIDE.md) · [v1.3.0 release report](HCLS_AI_FACTORY_v1.3.0_RELEASE_REPORT.md)

> Some research/branding/demo documents predate the current **Eight Engines · Eight Agents** framing
> and the `core/` restructure. The canonical references are the root README, CONTRIBUTING,
> [STRUCTURE.md](STRUCTURE.md), and the [build documentation](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/README.md).

## Verify any claim in these docs

```bash
.venv/bin/python scripts/run_all_tests.py        # 17 subjects · 8,402 passed
.venv/bin/python scripts/validate_registry.py    # port convention + supervisor cross-check
.venv/bin/python scripts/subject_facts.py        # the numbers the guides are generated from
```
