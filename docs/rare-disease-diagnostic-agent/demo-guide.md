# Rare Disease Diagnostic Agent -- Demo Guide

**Version:** 1.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Demo Overview](#1-demo-overview)
2. [Starting the Agent](#2-starting-the-agent)
3. [Tab 1: Patient Intake](#3-tab-1-patient-intake)
4. [Tab 2: Differential Diagnosis](#4-tab-2-differential-diagnosis)
5. [Tab 3: Variant Review](#5-tab-3-variant-review)
6. [Tab 4: Therapeutic Options](#6-tab-4-therapeutic-options)
7. [Tab 5: Reports](#7-tab-5-reports)
8. [Sample HPO Queries](#8-sample-hpo-queries)
9. [ACMG Classification Demo](#9-acmg-classification-demo)
10. [Gene Therapy Matching Demo](#10-gene-therapy-matching-demo)
11. [API Demo Scenarios](#11-api-demo-scenarios)
12. [Demo Tips](#12-demo-tips)

---

## 1. Demo Overview

The Rare Disease Diagnostic Agent demonstrates AI-powered rare disease diagnosis across 5 interactive tabs. The demo showcases phenotype-driven differential diagnosis, ACMG variant classification, orphan drug matching, gene therapy eligibility assessment, and multi-format report generation.

### Demo Flow (Recommended)

```
1. Patient Intake (HPO terms + clinical notes)
      |
2. Differential Diagnosis (ranked disease candidates)
      |
3. Variant Review (ACMG classification of key variants)
      |
4. Therapeutic Options (orphan drugs + gene therapies + trials)
      |
5. Reports (export diagnostic summary)
```

### Demo Duration

| Version | Duration | Focus |
|---|---|---|
| Quick Demo | 5 minutes | Phenotype query -> differential -> therapy |
| Standard Demo | 15 minutes | All 5 tabs with 2-3 scenarios |
| Deep Dive | 30 minutes | ACMG classification, gene therapy matching, API walkthrough |

---

## 2. Starting the Agent

```bash
# Navigate to agent directory
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/rare_disease_diagnostic_agent

# Option 1: Docker Compose (recommended)
docker compose up -d

# Option 2: Manual start
uvicorn api.main:app --host 0.0.0.0 --port 8134 &
streamlit run app/diagnostic_ui.py --server.port 8544 &
```

### Verify Services

| Service | URL | Expected |
|---|---|---|
| API | http://localhost:8134/health | `{"status": "healthy"}` |
| UI | http://localhost:8544 | Streamlit interface |
| Collections | http://localhost:8134/collections | 14 collections listed |

---

## 3. Tab 1: Patient Intake

### Purpose

Enter patient demographic data, HPO phenotype terms, clinical notes, family history, and optional VCF path for genomic analysis.

### Demo Steps

1. Open http://localhost:8544
2. Navigate to **Patient Intake** tab
3. Enter HPO terms (see Sample Queries below)
4. Add clinical notes describing the presentation
5. Set urgency level and workflow type (or leave auto-detect)
6. Click **Submit for Analysis**

### Demo Scenario: Suspected Marfan Syndrome

```
HPO Terms: HP:0001166, HP:0001519, HP:0004382, HP:0001083
Clinical Notes: "16-year-old male, tall stature with arm span exceeding
  height by 8 cm. Long, thin fingers. Slit lamp exam reveals bilateral
  ectopia lentis. Echocardiogram shows aortic root dilation (42 mm,
  Z-score 3.2) with mitral valve prolapse. Father has similar habitus."
Age: 16 years
Sex: Male
Family History: "Father with similar body habitus, underwent aortic
  root replacement at age 38. Paternal grandmother died suddenly at 42."
```

---

## 4. Tab 2: Differential Diagnosis

### Purpose

View ranked disease candidates matched to the patient's phenotype profile.

### What to Highlight

- **Similarity scores:** Each candidate disease shows a 0-1 score reflecting phenotype match quality
- **Matched phenotypes:** HPO terms the patient shares with the disease
- **Unmatched phenotypes:** Disease features not observed in the patient (important for targeted evaluation)
- **Inheritance pattern:** Autosomal dominant, recessive, X-linked, etc.
- **Causal genes:** Known genes to guide confirmatory testing

### Expected Output (Marfan Scenario)

| Rank | Disease | Score | Matched HPO | Inheritance |
|---|---|---|---|---|
| 1 | Marfan syndrome | 0.92 | HP:0001166, HP:0001519, HP:0004382, HP:0001083 | AD |
| 2 | Loeys-Dietz syndrome | 0.71 | HP:0001166, HP:0004382 | AD |
| 3 | Homocystinuria | 0.58 | HP:0001166, HP:0001519, HP:0001083 | AR |
| 4 | Ehlers-Danlos (vascular) | 0.45 | HP:0001166 | AD |

---

## 5. Tab 3: Variant Review

### Purpose

Classify genetic variants using ACMG/AMP criteria.

### Demo Steps

1. Navigate to **Variant Review** tab
2. Enter variant details:
   - Gene: FBN1
   - Variant type: missense
   - Population frequency: 0.00001
   - De novo: Yes (confirmed)
   - ClinVar: pathogenic
   - Computational prediction: damaging
3. Click **Classify Variant**

### Expected Output

```
Classification: PATHOGENIC
Criteria Met: PS1, PS2, PM2, PP3
Score: Path=11, Benign=0
Evidence Summary: "PS1: Same amino acid change as established pathogenic
  (score=4). PS2: De novo confirmed (score=4). PM2: Absent from controls,
  freq 0.00001 < 0.0001 (score=2). PP3: Computational evidence supports
  deleterious (score=1). Total pathogenic score: 11. Classification: pathogenic."
```

### Additional Variant Scenarios

**Benign variant (BA1):**
- Gene: Any
- Population frequency: 0.08 (8%)
- Result: BENIGN (BA1 standalone)

**VUS:**
- Gene: FBN1
- Variant type: missense
- Population frequency: 0.001
- Computational: "tolerated"
- Result: VUS (conflicting evidence)

---

## 6. Tab 4: Therapeutic Options

### Purpose

Match the patient's diagnosis and genotype to available therapies, gene therapies, and clinical trials.

### Demo Steps

1. Navigate to **Therapeutic Options** tab
2. Enter disease: "Spinal Muscular Atrophy"
3. Enter genotype (optional): "SMN1 deletion, 2 SMN2 copies"
4. Click **Search Therapies**

### Expected Output

| Therapy | Status | Mechanism | Gene |
|---|---|---|---|
| Nusinersen (Spinraza) | FDA Approved | SMN2 splicing modifier | SMN1 |
| Risdiplam (Evrysdi) | FDA Approved | SMN2 splicing modifier | SMN1 |
| Zolgensma (onasemnogene) | FDA Approved | AAV9 gene replacement | SMN1 |

### Gene Therapy Highlight

Point out that Zolgensma is a one-time gene therapy that delivers a functional copy of SMN1 via AAV9 vector. Eligibility typically requires age < 2 years and confirmed biallelic SMN1 mutations.

---

## 7. Tab 5: Reports

### Purpose

Export diagnostic findings as structured reports.

### Demo Steps

1. Navigate to **Reports** tab
2. Select format: Markdown, JSON, or PDF
3. Click **Generate Report**
4. Download the generated file

### Report Contents

- Patient summary
- HPO terms with IC scores
- Ranked differential diagnosis
- Variant classifications (if available)
- Therapeutic options
- Recommended next steps
- Guideline references

---

## 8. Sample HPO Queries

### Scenario 1: Dravet Syndrome (Neurogenetic)

```
HPO Terms: HP:0001250, HP:0001263, HP:0001249, HP:0001252
Clinical Notes: "8-month-old female with first prolonged febrile seizure
  lasting 45 minutes. Previously healthy. Subsequent afebrile seizures at
  10 and 12 months. Developmental plateau noted at 14 months."
Age: 14 months
```

**Expected top candidate:** Dravet Syndrome (SCN1A)

### Scenario 2: Cystic Fibrosis (Metabolic)

```
HPO Terms: HP:0002205, HP:0002110, HP:0001508
Clinical Notes: "3-year-old male with recurrent pneumonias, chronic cough,
  and failure to thrive. Positive newborn screening for IRT. Sweat chloride
  72 mmol/L (elevated). Sputum culture positive for Pseudomonas aeruginosa."
Age: 3 years
```

**Expected top candidate:** Cystic Fibrosis (CFTR)

### Scenario 3: Gaucher Disease (Hematologic/Metabolic)

```
HPO Terms: HP:0001744, HP:0002240, HP:0001882
Clinical Notes: "12-year-old Ashkenazi Jewish female with progressive
  splenomegaly and hepatomegaly. Pancytopenia on CBC. Bone pain in
  bilateral femurs. No neurological symptoms."
Age: 12 years
```

**Expected top candidate:** Gaucher Disease Type 1 (GBA1)

### Scenario 4: Ehlers-Danlos Syndrome (Connective Tissue)

```
HPO Terms: HP:0001382, HP:0000974, HP:0001252
Clinical Notes: "25-year-old female with lifelong joint hypermobility
  (Beighton score 8/9), hyperextensible skin, easy bruising. Multiple
  joint dislocations. Chronic pain syndrome."
Age: 25 years
```

**Expected top candidate:** Ehlers-Danlos Syndrome, Classical (COL5A1)

### Scenario 5: Long QT Syndrome (Cardiac)

```
HPO Terms: HP:0001657, HP:0004756, HP:0001279, HP:0001695
Clinical Notes: "14-year-old male with syncope during swimming. Family
  history: mother with prolonged QT on ECG, maternal uncle died suddenly
  at age 22. Patient ECG shows QTc 520 ms."
Age: 14 years
```

**Expected top candidate:** Long QT Syndrome Type 1 (KCNQ1)

---

## 9. ACMG Classification Demo

### Full ACMG Walkthrough

Demonstrate the variant classification pipeline with a classic pathogenic variant:

```json
{
  "gene": "SCN1A",
  "variant_type": "nonsense",
  "population_frequency": 0.0,
  "is_de_novo": true,
  "de_novo_confirmed": true,
  "computational_prediction": "damaging",
  "in_clinvar": true,
  "clinvar_classification": "pathogenic",
  "phenotype_specific": true
}
```

**Walk through each criterion:**
1. **PVS1 (+8):** Null variant (nonsense) in SCN1A, which is LOF-intolerant
2. **PS1 (+4):** ClinVar reports this variant as pathogenic
3. **PS2 (+4):** De novo, confirmed with parental testing
4. **PM2 (+2):** Absent from population databases (freq = 0)
5. **PP3 (+1):** Computational tools predict damaging
6. **PP4 (+1):** Phenotype (Dravet syndrome) is highly specific for SCN1A

**Total: 20 points -> PATHOGENIC**

---

## 10. Gene Therapy Matching Demo

### Showcase Available Gene Therapies

Demonstrate therapy matching for multiple diseases:

**SMA:**
```
Disease: Spinal Muscular Atrophy
-> Zolgensma (gene replacement), Spinraza (splicing modifier), Evrysdi (splicing modifier)
```

**Sickle Cell Disease:**
```
Disease: Sickle Cell Disease
-> Casgevy (CRISPR editing), Lyfgenia (lentiviral gene addition)
```

**Hemophilia B:**
```
Disease: Hemophilia B
-> Hemgenix (AAV5 gene replacement for Factor IX)
```

**Highlight genotype-specific matching:**
```
Disease: Cystic Fibrosis
Genotype: F508del homozygous
-> Trikafta (eligible: at least one F508del allele)
-> Orkambi (eligible: F508del homozygous specific)
```

---

## 11. API Demo Scenarios

### 11.1 Health Check

```bash
curl http://localhost:8134/health | python -m json.tool
```

### 11.2 Knowledge Version

```bash
curl http://localhost:8134/v1/diagnostic/knowledge-version | python -m json.tool
```

### 11.3 Disease Categories

```bash
curl http://localhost:8134/v1/diagnostic/disease-categories | python -m json.tool
```

### 11.4 Gene Therapies

```bash
curl http://localhost:8134/v1/diagnostic/gene-therapies | python -m json.tool
```

### 11.5 Diagnostic Query

```bash
curl -X POST http://localhost:8134/v1/diagnostic/diagnose \
  -H "Content-Type: application/json" \
  -d '{
    "hpo_terms": ["HP:0001250", "HP:0001263", "HP:0001249"],
    "clinical_notes": "Infant with seizures and developmental delay",
    "age": "8 months",
    "top_k": 5
  }' | python -m json.tool
```

### 11.6 ACMG Classification

```bash
curl -X POST http://localhost:8134/v1/diagnostic/variants/interpret \
  -H "Content-Type: application/json" \
  -d '{
    "variants": [{
      "gene": "FBN1",
      "variant_type": "missense",
      "population_frequency": 0.00001,
      "is_de_novo": true,
      "de_novo_confirmed": true,
      "computational_prediction": "damaging"
    }]
  }' | python -m json.tool
```

---

## 12. Demo Tips

### Talking Points

1. **The diagnostic odyssey:** "Rare disease patients wait an average of 5-7 years for diagnosis. This agent compresses that timeline by integrating phenotype matching, variant classification, and therapy identification into a single platform."

2. **HPO-based matching:** "Each phenotype is weighted by its Information Content -- rare, specific findings like ectopia lentis carry more diagnostic weight than common ones like seizures."

3. **ACMG rigor:** "Variant classification follows the same ACMG/AMP framework used by clinical labs, with 28 criteria systematically evaluated and scored."

4. **Gene therapy revolution:** "12 gene therapies are now approved for rare diseases. The agent matches patients to eligible therapies based on disease, genotype, and access pathway."

5. **Platform integration:** "The agent connects to the broader HCLS AI Factory -- cardiology, pharmacogenomics, and clinical trial agents can be consulted automatically."

### Common Questions

**Q: Is this for clinical use?**
A: This is a research and demonstration platform. Clinical use would require HIPAA compliance, regulatory validation, and integration with clinical workflows. The underlying algorithms (ACMG, HPO matching) follow published clinical standards.

**Q: How current is the knowledge base?**
A: The knowledge base version is 1.0.0 (March 2026). The RAG architecture allows incremental updates by adding records to vector collections without retraining.

**Q: Can it handle novel diseases?**
A: The phenotype-driven approach can identify candidate genes even for uncharacterized diseases by matching symptom patterns to known gene-phenotype associations. The "Undiagnosed Disease" workflow is specifically designed for these cases.

---

*Apache 2.0 License -- HCLS AI Factory*
