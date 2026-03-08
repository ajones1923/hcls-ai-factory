# Precision Biomarker Intelligence Agent — Demo Guide

**Author:** Adam Jones
**Date:** March 2026

> Step-by-step walkthrough of all 8 Streamlit UI tabs using both sample patients. Use this guide for live demos.

---

## Prerequisites

Before starting the demo, ensure:

1. **Milvus** is running on `localhost:19530` with 14 seeded collections (652 vectors)
2. **Streamlit UI** is running on `http://localhost:8528`
3. **Anthropic API key** is configured (for RAG queries in Tab 5)

Verify everything is ready:

```bash
cd ai_agent_adds/precision_biomarker_agent
python3 scripts/demo_validation.py
# → 65/65 checks passed
```

---

## Demo Flow Overview

| Tab | Time | Key Talking Points |
|---|---|---|
| 1. Biomarker Analysis | 3 min | Critical alerts, discordance detection, multi-lab comparison |
| 2. Biological Age | 2 min | PhenoAge/GrimAge, aging drivers, confidence intervals |
| 3. Disease Trajectories | 3 min | 9-domain risk cards, genotype integration |
| 4. Pharmacogenomics | 3 min | CYP2D6 *1/*4, drug recommendations, CPIC guidelines |
| 5. Evidence Explorer | 2 min | RAG search across 14 collections, citation scoring |
| 6. Export | 2 min | FHIR R4 validation, PDF report, Markdown |
| 7. Genotype Adjustments | 2 min | Age-stratified ranges, genotype-aware thresholds |
| 8. Patient Data | 1 min | Full clinical context, family history, medications |
| **Total** | **~18 min** | |

---

## Patient 1: Male, 45, Ashkenazi Jewish (HCLS-BIO-2026-00001)

### Key Narrative

> A 45-year-old Ashkenazi Jewish male software engineer with ApoE E3/E4 genotype, MTHFR C677T heterozygous, CYP2D6 *1/*4 intermediate metabolizer. Family history of MI (father at 58), Alzheimer's (paternal grandmother at 74), and T2DM (mother at 52). Currently on 7 medications including atorvastatin and lisinopril. Presents with fatigue, joint stiffness, and intermittent brain fog.

### Tab 1: Biomarker Analysis

1. **Load sample patient** — Click "Load Sample Patient" button → selects Male, 45
2. **Click "Run Analysis"**
3. **Critical Values section** — Point out: No critical alerts on this patient (values within normal ranges). Demonstrate what happens with extreme values by mentioning the glucose=450 test.
4. **Discordance Detection** — Show any detected cross-biomarker discordances. Key example: elevated Lp(a) of 85 nmol/L with borderline LDL of 138 — multiplicative cardiovascular risk.
5. **Lab Range Comparison** — Show the three-way Quest vs LabCorp vs Function Health comparison. Point out where "lab normal" differs from "optimal" — e.g., Vitamin D at 38 is Quest-normal but Function Health wants >50.

**Key talking point:** *"Standard lab reports would mark most of these as normal. But when we layer in the ApoE E4 genotype, the LDL of 138 actually needs to be under 100. This is what genomics-informed interpretation means."*

### Tab 2: Biological Age

1. Navigate to Tab 2
2. **PhenoAge** — Show biological age vs chronological age (45). Discuss age acceleration.
3. **Top Aging Drivers** — Show which biomarkers contribute most to aging (positive = aging, negative = protective). Point out albumin and CRP contributions.
4. **Confidence Interval** — Show the 95% CI and explain standard error of ~4.9 years with 7/9 biomarkers.
5. **GrimAge** — Note this returns None for this patient (no plasma markers like GDF-15 or Cystatin C in standard panel). Explain this would require a specialty lab panel.

**Key talking point:** *"PhenoAge uses 9 clinical biomarkers to estimate biological age. This patient's age acceleration tells us whether his cellular aging is faster or slower than expected for 45."*

### Tab 3: Disease Trajectories

1. Navigate to Tab 3
2. **9 risk cards** — Show all 9 disease domains sorted by risk level
3. **Cardiovascular (MODERATE)** — LDL 138 + ApoE E4 + Lp(a) 85 + family history MI at 58. Point out the genotype-adjusted LDL threshold.
4. **Type 2 Diabetes (MODERATE)** — HbA1c 5.6 (pre-diabetes range) + TCF7L2 CT (1 risk allele) + HOMA-IR 1.98. Show genotype-adjusted glucose thresholds.
5. **Cognitive (MODERATE)** — ApoE E3/E4 + family history Alzheimer's at 74. Show modifiable risk factor count.
6. **Nutritional** — Omega-3 Index 5.8 (suboptimal, target >8%), Vitamin D 38 (adequate but not optimal).

**Key talking point:** *"Traditional risk calculators don't account for genotype. TCF7L2 CT shifts the diabetes trajectory left — meaning this patient hits risk thresholds at lower glucose levels than a wild-type individual."*

### Tab 4: Pharmacogenomics

1. Navigate to Tab 4
2. **Enter star alleles:** CYP2D6 = `*1/*4`, CYP2C19 = `*1/*2`, TPMT = `*1/*1`
3. **Enter genotypes:** MTHFR_rs1801133 = `CT`
4. **Click "Map PGx"**
5. **CYP2D6 Intermediate Metabolizer** — Show affected drugs: codeine (reduced efficacy), tramadol (dose adjustment), tamoxifen (reduced activation).
6. **CYP2C19 Intermediate Metabolizer** — Show clopidogrel warning (reduced activation).
7. **Drug-drug interactions** — Show any detected interactions across PGx recommendations.
8. **CPIC guideline versions** — Point out the audit trail with publication dates and PMIDs.

**Key talking point:** *"This patient is a CYP2D6 intermediate metabolizer. Codeine would have reduced efficacy because he can't convert it to morphine efficiently. This is exactly the kind of finding that prevents adverse drug events."*

### Tab 5: Evidence Explorer

1. Navigate to Tab 5
2. **Query:** "ApoE E4 carrier with elevated LDL — cardiovascular risk management"
3. **Show results** — Point out collection sources, relevance scores, and citation levels
4. **Filter by collection** — Demonstrate filtering to just `biomarker_genetic_variants` or `biomarker_clinical_evidence`
5. **Second query:** "MTHFR C677T heterozygous folate metabolism homocysteine"

**Key talking point:** *"The RAG engine searches across 14 specialized collections simultaneously — from genetic variant databases to clinical evidence to drug interaction tables — and synthesizes a grounded answer with citations."*

### Tab 6: Export

1. Navigate to Tab 6
2. **FHIR R4** — Generate FHIR bundle, show validation passes (0 errors)
3. **Show bundle structure** — Patient resource, DiagnosticReport, Observations
4. **PDF** — Download clinical report (requires reportlab)
5. **Markdown** — Show plain-text clinical summary

**Key talking point:** *"The FHIR R4 export produces a structurally validated diagnostic report that can be ingested by any EHR system. Every reference resolves within the bundle — this is interoperability-ready."*

### Tab 7: Genotype Adjustments

1. Navigate to Tab 7
2. **Genotype adjustments** — Show how ApoE E4 modifies LDL target from <130 to <100
3. **Age-stratified ranges** — Show creatinine, TSH, and other ranges for the 40-59 bracket
4. **Comparison** — Point out where standard ranges and age-adjusted ranges differ

**Key talking point:** *"A creatinine of 1.3 in a 25-year-old is concerning. In a 70-year-old, it's expected. Age-stratified ranges prevent false positives and catch age-inappropriate values that standard ranges miss."*

### Tab 8: Patient Data

1. Navigate to Tab 8
2. **Browse patient data** — Show full demographics, medications (7), family history
3. **Family history** — Father MI at 58, paternal grandmother Alzheimer's at 74 (ApoE E4/E4 homozygous), maternal uncle colorectal cancer at 61
4. **Medications** — Show how atorvastatin, lisinopril, and supplements map to biomarker effects

---

## Patient 2: Female, 38, Ashkenazi Jewish (HCLS-BIO-2026-00002)

### Key Narrative

> A 38-year-old Ashkenazi Jewish genetic counselor with active preconception planning (12-18 months). Mother has BRCA1 185delAG confirmed breast cancer at 48. BRCA1 status NOT YET TESTED — URGENT. Ferritin of 28 ng/mL is critically low for preconception. GBA carrier risk 50% (paternal grandmother had Gaucher Disease Type 1).

### Key Demo Points for Patient 2

1. **Tab 1 — Ferritin 28** — Show lab ranges: Quest says "normal" (12-150), but Function Health optimal for preconception is >50. This is the value of multi-lab comparison.

2. **Tab 3 — Iron trajectory** — Ferritin 28 + transferrin saturation 18% → iron depletion stage. For preconception, target ferritin >50-70.

3. **Tab 3 — Nutritional trajectory** — Omega-3 Index 4.9% (target >8% for pregnancy), Vitamin D 32 (target 40-60).

4. **Tab 4 — PGx** — All normal metabolizers (CYP2D6 *1/*1, CYP2C19 *1/*1, TPMT *1/*1). Important for preconception medication safety.

5. **Tab 3 — Diabetes trajectory** — TCF7L2 TT (2 risk alleles) — highest genetic risk category. Father has T2DM. Even though HbA1c is 5.2 and normal now, this patient needs early monitoring.

6. **Tab 8 — OB/GYN context** — BRCA1 NOT YET TESTED, preconception ACTIVE, GBA carrier risk 50%. This drives urgency for genetic testing.

**Key talking point:** *"This is where precision medicine becomes actionable. Her Quest lab report says ferritin 28 is normal. But for a woman actively planning pregnancy, ferritin needs to be above 50. And BRCA1 testing is urgent given her family history — this should happen before any pregnancy planning."*

---

## Troubleshooting

| Issue | Resolution |
|---|---|
| Streamlit not starting | Check port 8528 is free: `lsof -i :8528` |
| "No collections found" | Run `python3 scripts/seed_all.py` |
| FHIR validation errors | Verify export.py has Patient resource in bundle |
| RAG returns empty | Check Milvus is running: `curl localhost:19530/healthz` |
| PGx shows no results | Ensure star_alleles use format `*1/*4` (with asterisks) |
| GrimAge returns None | Expected — requires specialty plasma markers not in standard panels |

---

## Quick Reset

If you need to reset the demo environment:

```bash
# Re-seed all collections (drops and recreates)
python3 scripts/seed_all.py

# Verify
python3 scripts/demo_validation.py
# → 65/65 passed

# Restart UI
# Ctrl+C the streamlit process, then:
streamlit run app/biomarker_ui.py --server.port 8528
```
