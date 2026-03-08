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
| 3. Disease Risk | 3 min | 9-domain risk cards, genotype integration |
| 4. PGx Profile | 3 min | CYP2D6 *1/*4, drug recommendations, CPIC guidelines |
| 5. Evidence Explorer | 2 min | RAG search across 14 collections, citation scoring |
| 6. Reports | 2 min | FHIR R4 validation, PDF report, Markdown |
| 7. Patient 360 | 2 min | Unified cross-agent dashboard |
| 8. Longitudinal | 1 min | Multi-visit biomarker trending |
| **Total** | **~18 min** | |

---

## Patient 1: Male, 45, Ashkenazi Jewish (HCLS-BIO-2026-00001)

### Key Narrative

> A 45-year-old Ashkenazi Jewish male software engineer with ApoE E3/E4 genotype, MTHFR C677T heterozygous, CYP2D6 *1/*4 intermediate metabolizer. Family history of MI (father at 58), Alzheimer's (paternal grandmother at 74), and T2DM (mother at 52). Currently on 7 medications including atorvastatin and lisinopril. Presents with fatigue, joint stiffness, and intermittent brain fog.

### Tab 1: Biomarker Analysis

1. **Load sample patient** — Click **"Load Male Patient (HG002)"** button
2. **Click "Run Full Analysis"**
3. **Critical Values section** — Point out: No critical alerts on this patient (values within normal ranges). Demonstrate what happens with extreme values by mentioning the glucose=450 test.
4. **Discordance Detection** — Show any detected cross-biomarker discordances. Key example: elevated Lp(a) of 85 nmol/L with borderline LDL of 138 — multiplicative cardiovascular risk.
5. **Lab Range Comparison** — Show the three-way Quest vs LabCorp vs Function Health comparison. Point out where "lab normal" differs from "optimal" — e.g., Vitamin D at 38 is Quest-normal but Function Health wants >50.

**Key talking point:** *"Standard lab reports would mark most of these as normal. But when we layer in the ApoE E4 genotype, the LDL of 138 actually needs to be under 100. This is what genomics-informed interpretation means."*

### Tab 2: Biological Age

1. Navigate to Tab 2
2. **PhenoAge** — Show biological age vs chronological age (45). Discuss age acceleration.
3. **Top Aging Drivers** — Show which biomarkers contribute most to aging (positive = aging, negative = protective). Point out albumin and CRP contributions.
4. **Confidence Interval** — Show the 95% CI. This patient has all 9/9 PhenoAge biomarkers available, giving a standard error of ~4.9 years — the tightest possible estimate.
5. **GrimAge** — Note this returns None for this patient (no plasma markers like GDF-15 or Cystatin C in standard panel). Explain this would require a specialty lab panel.

**Key talking point:** *"PhenoAge uses 9 clinical biomarkers to estimate biological age. This patient's age acceleration tells us whether his cellular aging is faster or slower than expected for 45."*

### Tab 3: Disease Risk

1. Navigate to Tab 3
2. **9 risk cards** — Show all 9 disease domains sorted by risk level
3. **Cardiovascular (MODERATE)** — LDL 138 + ApoE E4 + Lp(a) 85 + family history MI at 58. Point out the genotype-adjusted LDL threshold.
4. **Nutritional (MODERATE)** — Omega-3 Index 5.8 (suboptimal, target >8%), Vitamin D 38 (adequate but not optimal).
5. **Cognitive (MODERATE)** — ApoE E3/E4 + family history Alzheimer's at 74. Show modifiable risk factor count.
6. **Type 2 Diabetes (LOW)** — HbA1c 5.6 is just below the pre-diabetes threshold (5.7%), but note the TCF7L2 CT genotype (1 risk allele) means this patient should be monitored more closely than a wild-type individual with the same HbA1c.

**Key talking point:** *"Three domains flag MODERATE risk. The cardiovascular finding is particularly important — traditional LDL cutoffs say 138 is borderline, but ApoE E4 carriers need LDL under 100. Without genotype context, this patient would be told he's fine."*

### Tab 4: PGx Profile

1. Navigate to Tab 4
2. **Enter star alleles:** CYP2D6 = `*1/*4`, CYP2C19 = `*1/*2`, TPMT = `*1/*1`
3. **Enter genotypes:** MTHFR_rs1801133 = `CT`
4. **Click "Map Drug Interactions"**
5. **CYP2D6 Intermediate Metabolizer** — Show affected drugs: codeine (reduced efficacy), tramadol (dose adjustment), tamoxifen (reduced activation).
6. **CYP2C19 Intermediate Metabolizer** — Show clopidogrel warning (reduced activation).
7. **Drug-drug interactions** — Show any detected interactions across PGx recommendations.
8. **CPIC guideline versions** — Point out the audit trail with publication dates and PMIDs.

**Key talking point:** *"This patient is a CYP2D6 intermediate metabolizer. Codeine would have reduced efficacy because he can't convert it to morphine efficiently. This is exactly the kind of finding that prevents adverse drug events."*

### Tab 5: Evidence Explorer

1. Navigate to Tab 5
2. **Query:** "ApoE E4 carrier with elevated LDL — cardiovascular risk management"
3. **Show results** — Point out collection sources, relevance scores, and citation levels
4. **Filter by collection** — Expand "Collection Filters" and demonstrate filtering to just `biomarker_genetic_variants` or `biomarker_clinical_evidence`
5. **Second query:** "MTHFR C677T heterozygous folate metabolism homocysteine"

**Key talking point:** *"The RAG engine searches across 14 specialized collections simultaneously — from genetic variant databases to clinical evidence to drug interaction tables — and synthesizes a grounded answer with citations."*

### Tab 6: Reports

1. Navigate to Tab 6
2. **FHIR R4** — Generate FHIR bundle, show validation passes (0 errors)
3. **Show bundle structure** — Patient resource, DiagnosticReport, Observations
4. **PDF** — Download clinical report (requires reportlab)
5. **Markdown** — Show plain-text clinical summary

**Key talking point:** *"The FHIR R4 export produces a structurally validated diagnostic report that can be ingested by any EHR system. Every reference resolves within the bundle — this is interoperability-ready."*

### Tab 7: Patient 360

1. Navigate to Tab 7
2. **Unified dashboard** — This tab provides a cross-agent intelligence view combining genomics, biomarkers, drug candidates, and clinical evidence from across the HCLS AI Factory platform.
3. **Click "Load Demo Patient 360"** to populate the dashboard with the current patient's data.
4. **Point out cross-pipeline integration** — Show how biomarker results connect to genomic variants and potential drug candidates from other pipeline stages.

**Key talking point:** *"Patient 360 is where all three pipeline stages converge — genomic variants from Stage 1, biomarker intelligence from Stage 2, and drug candidates from Stage 3 — giving clinicians a single unified view."*

### Tab 8: Longitudinal

1. Navigate to Tab 8
2. **Multi-visit tracking** — Show biomarker trends across multiple time points
3. **Trend analysis** — Point out improving, stable, and crisis patterns
4. **Clinical context** — Explain how longitudinal tracking reveals trajectories that single-point-in-time analysis misses

**Key talking point:** *"A single lab draw is a snapshot. Longitudinal tracking reveals the trajectory — is ferritin trending down toward depletion? Is HbA1c creeping up despite medication? These trends drive proactive intervention."*

---

## Patient 2: Female, 38, Ashkenazi Jewish (HCLS-BIO-2026-00002)

### Key Narrative

> A 38-year-old Ashkenazi Jewish genetic counselor with active preconception planning (12-18 months). Mother has BRCA1 185delAG confirmed breast cancer at 48. BRCA1 status NOT YET TESTED — URGENT. Ferritin of 28 ng/mL is critically low for preconception. GBA carrier risk 50% (paternal grandmother had Gaucher Disease Type 1).

### Key Demo Points for Patient 2

1. **Tab 1 — Ferritin 28** — Show lab ranges: Quest says "normal" (12-150), but Function Health optimal for preconception is >50. This is the value of multi-lab comparison.

2. **Tab 3 — Iron trajectory (LOW)** — Ferritin 28 + transferrin saturation 18%. While the engine classifies this as LOW risk, point out that for preconception planning, the clinical target is ferritin >50-70 — making this a priority for optimization.

3. **Tab 3 — Nutritional (MODERATE)** — Omega-3 Index 4.9% (target >8% for pregnancy), Vitamin D 32 (target 40-60). This is the only MODERATE-risk domain for this patient.

4. **Tab 4 — PGx** — All normal metabolizers (CYP2D6 *1/*1, CYP2C19 *1/*1, TPMT *1/*1). Important for preconception medication safety.

5. **Tab 3 — Diabetes trajectory (LOW)** — TCF7L2 TT (2 risk alleles) is the highest genetic risk category, and father has T2DM. HbA1c 5.2 is reassuringly normal now, but the genetic burden warrants monitoring during and after pregnancy when insulin resistance naturally increases.

6. **Tab 7 — Patient 360** — BRCA1 NOT YET TESTED, preconception ACTIVE, GBA carrier risk 50%. This drives urgency for genetic testing before pregnancy planning proceeds.

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
