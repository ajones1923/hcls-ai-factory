# Learning Guide -- Foundations for the Precision Biomarker Intelligence Agent

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0

---

## Table of Contents

1. [Welcome](#chapter-1-welcome)
2. [What Are Biomarkers?](#chapter-2-what-are-biomarkers)
3. [The Data Challenge](#chapter-3-the-data-challenge)
4. [What Is RAG?](#chapter-4-what-is-rag)
5. [System Overview](#chapter-5-system-overview)
6. [Your First Query](#chapter-6-your-first-query)
7. [Understanding Collections](#chapter-7-understanding-collections)
8. [The Knowledge Graph](#chapter-8-the-knowledge-graph)
9. [Pharmacogenomics Explained](#chapter-9-pharmacogenomics-explained)
10. [Setting Up Locally](#chapter-10-setting-up-locally)
11. [Exploring the API](#chapter-11-exploring-the-api)
12. [Understanding the Codebase](#chapter-12-understanding-the-codebase)
13. [Next Steps](#chapter-13-next-steps)
14. [Glossary](#chapter-14-glossary)

---

## Chapter 1: Welcome

### Who This Guide Is For

This guide is written for developers, data scientists, bioinformaticians, and
anyone curious about how software can transform raw blood work and genetic data
into actionable health intelligence. You do not need a medical degree to follow
along. You do not need to be a machine-learning expert. If you can read Python
and have used a REST API before, you have everything you need.

The Precision Biomarker Intelligence Agent is one of five intelligence agents in
the HCLS AI Factory platform. It sits between the Genomics Pipeline (which turns
raw DNA sequencing into variant calls) and the Drug Discovery Pipeline (which
finds candidate molecules). Its job: take a patient's biomarker values and
genetic variants and produce a comprehensive, evidence-grounded health
intelligence report in seconds.

### What You Will Learn

By the end of this guide you will understand:

- What biomarkers are and why they matter for precision health.
- Why population-average reference ranges are insufficient for individual care.
- How Retrieval-Augmented Generation (RAG) works and why it matters here.
- The architecture of the Biomarker Intelligence Agent -- its 14 vector
  collections, 4 analysis modules, and the plan-analyze-search-synthesize
  pipeline.
- How to run your first query through the Streamlit UI.
- Every collection in the system and what its fields represent.
- The knowledge graph that powers domain reasoning.
- Pharmacogenomics: star alleles, metabolizer phenotypes, and CPIC guidelines.
- How to set up the system locally with Docker.
- How to explore the REST API with curl.
- The complete project structure, file by file.

### Prerequisites

- **Python 3.10+** -- You should be comfortable reading Python code.
- **Docker and Docker Compose** -- For running Milvus and the application stack.
- **Basic command-line skills** -- curl, environment variables, port management.
- **Curiosity** -- The most important prerequisite. Biomarker science is deep,
  and this guide will give you footholds, not final answers.

No prior knowledge of biology, genetics, or clinical laboratory science is
required. Every domain term is explained when first introduced and collected in
the Glossary at the end.

---

## Chapter 2: What Are Biomarkers?

A **biomarker** is any measurable substance in your body that indicates
something about your health. When a doctor orders blood work, every number on
that report is a biomarker. They are the quantitative language of medicine.

### Five Categories of Biomarkers

Biomarkers span a wide range of body systems. This agent organizes them into
five broad categories:

**1. Metabolic Biomarkers**

These track how your body processes energy. They include fasting glucose (how
much sugar is in your blood after an overnight fast), HbA1c (a three-month
average of blood sugar), insulin, and HOMA-IR (a calculated index of insulin
resistance). A fasting glucose of 100 mg/dL is normal. At 126 mg/dL, it
indicates diabetes. The numbers in between -- the "pre-diabetes" zone -- are
where early detection saves lives.

**2. Hematologic Biomarkers**

These describe the cellular composition of your blood. A complete blood count
(CBC) measures white blood cells (WBC), red blood cell width (RDW), mean
corpuscular volume (MCV), lymphocyte percentage, and hemoglobin. These markers
appear in both infection detection and -- perhaps surprisingly -- biological age
estimation. The PhenoAge clock uses WBC, RDW, MCV, and lymphocyte percentage as
inputs.

**3. Hormonal Biomarkers**

Hormones are chemical messengers. Thyroid-stimulating hormone (TSH), free T3,
and free T4 regulate metabolism. Insulin and leptin regulate appetite and energy
storage. Cortisol tracks stress. These biomarkers often have complex feedback
loops -- a high TSH with a low free T4 suggests hypothyroidism, but interpreting
TSH alone without the T4 can be misleading. This is why the agent has a
dedicated discordance detection module.

**4. Inflammatory Biomarkers**

Chronic inflammation underlies many diseases long before symptoms appear.
High-sensitivity C-reactive protein (hs-CRP) is the most widely used marker.
Elevated hs-CRP increases cardiovascular risk. Homocysteine, linked to MTHFR
gene variants, is both an inflammatory and cardiovascular marker. Ferritin, an
iron storage protein, rises with inflammation even when iron stores are normal,
making interpretation genotype-dependent.

**5. Genetic Biomarkers**

These are not measured from blood chemistry but from DNA. A single-nucleotide
polymorphism (SNP) like MTHFR C677T (rs1801133) changes how your body processes
folate. An APOE e4 allele changes your cardiovascular and neurological risk
profile. The HFE C282Y variant predisposes to hereditary hemochromatosis. Genetic
biomarkers do not change over time, but they permanently alter how every other
biomarker should be interpreted.

### Why Biomarkers Matter for Precision Health

Traditional medicine treats symptoms. Precision health uses biomarkers to detect
disease trajectories years before symptoms appear. A patient with an HbA1c of
5.7% is "normal" by standard definitions, but if they carry two risk alleles in
TCF7L2, their diabetes risk is significantly elevated, and the agent adjusts the
threshold downward to 5.5%. That distinction -- invisible without genetic
context -- is the entire point of this system.

---

## Chapter 3: The Data Challenge

### Population Averages vs. Individual Variation

Every reference range you see on a lab report (e.g., "Glucose: 70-100 mg/dL")
comes from a population study. Researchers measured thousands of apparently
healthy people and drew boundaries around the middle 95%. If your value falls
inside those boundaries, you are told you are "normal."

The problem: you are not a population. You are an individual with a unique
genome, ancestry, age, sex, diet, medication history, and environmental
exposure. A ferritin of 300 ng/mL is "normal" by most lab ranges (typically
24-336 for men) but could indicate early hereditary hemochromatosis in someone
carrying HFE C282Y homozygous. Conversely, a ferritin of 30 ng/mL is "normal"
but may indicate functional iron deficiency in a woman with heavy menstruation.

Population averages mask individual risk.

### The Pharmacogenomics Gap

Roughly 99% of prescriptions are written without checking the patient's
pharmacogenomic profile. A patient who is a CYP2D6 poor metabolizer (about 6-10%
of Caucasians) will not convert codeine into morphine and will get no pain
relief. A CYP2D6 ultra-rapid metabolizer given the same dose may convert it too
fast and experience respiratory depression. The drug is identical. The dose is
identical. The outcomes are opposite.

The Clinical Pharmacogenetics Implementation Consortium (CPIC) publishes
evidence-based guidelines mapping gene-drug pairs to dosing recommendations.
These guidelines exist. They are freely available. Yet they are rarely consulted
at the point of care because the data is scattered across multiple databases and
the interpretation requires specialized knowledge. This agent bridges that gap.

### Scattered Reference Data

A clinician trying to give a truly personalized interpretation needs to consult:

- **Standard reference ranges** from the lab report.
- **Genotype-adjusted ranges** from published studies (e.g., MTHFR impact on
  homocysteine, PNPLA3 impact on ALT).
- **Age- and sex-stratified ranges** (TSH thresholds differ for a 25-year-old
  vs. a 75-year-old).
- **CPIC guidelines** for any medications the patient takes.
- **Disease trajectory models** to estimate pre-symptomatic risk.
- **PhenoAge/GrimAge algorithms** to estimate biological aging.
- **Clinical literature** for the latest evidence on any of the above.
- **Nutritional genomics** data (MTHFR and methylfolate, FADS1 and omega-3
  conversion, VDR and vitamin D absorption).
- **Carrier screening panels** (especially for Ashkenazi Jewish populations
  with elevated carrier frequencies for conditions like Tay-Sachs, Gaucher,
  BRCA1/2).

This data lives in dozens of separate databases, PDFs, journal articles, and
clinical guidelines. No human can synthesize it all in a five-minute appointment.
This is the problem the Biomarker Intelligence Agent solves.

---

## Chapter 4: What Is RAG?

RAG stands for **Retrieval-Augmented Generation**. It is the core pattern that
allows this agent to answer questions with grounded, evidence-based responses
rather than hallucinated ones. Let us break it down into its components.

### Step 1: Embeddings

An **embedding** is a list of numbers (a vector) that captures the meaning of a
piece of text. The Biomarker Intelligence Agent uses the
`BAAI/bge-small-en-v1.5` model to convert text into 384-dimensional vectors.

When we say "HbA1c is a measure of three-month average blood glucose" and convert
that sentence into an embedding, we get a vector of 384 floating-point numbers.
When we convert "glycated hemoglobin reflects long-term glycemic control" into
an embedding, we get a different vector -- but one that is mathematically
*close* to the first, because the two sentences mean similar things.

This is the key insight: **semantic similarity becomes geometric proximity**.
Texts that mean similar things end up near each other in vector space.

### Step 2: Vector Similarity Search

Once we have converted all our reference data into embeddings and stored them in
a vector database (Milvus, in this system), we can search by meaning rather than
by keyword.

When a user asks "What does my HbA1c of 5.8% mean?", the system:

1. Converts the question into an embedding.
2. Searches all 14 collections in parallel for the vectors closest to that
   question embedding.
3. Returns the top results ranked by similarity score.

This is fundamentally different from keyword search. A keyword search for
"HbA1c" would miss documents that only mention "glycated hemoglobin" or
"A1c." Vector search finds them all because the embeddings are close.

### Step 3: Retrieval-Augmented Generation

Once the system has retrieved the most relevant evidence from its collections,
it passes that evidence to a large language model (Claude, from Anthropic) as
context. The LLM then generates a response that is grounded in the retrieved
evidence.

This is the "augmented" part. Without the retrieval step, an LLM would answer
from its general training data, which may be outdated, imprecise, or simply
wrong for a specific patient. With retrieval, the LLM has access to curated,
domain-specific evidence including CPIC guidelines, clinical literature,
genotype-specific reference ranges, and disease trajectory models.

The RAG pattern gives us the best of both worlds: the reasoning ability of an
LLM combined with the factual accuracy of a curated knowledge base.

### How the Agent Uses RAG

The Biomarker Intelligence Agent extends basic RAG in several important ways:

- **Multi-collection search**: 14 collections are searched simultaneously using
  parallel threads, not just one.
- **Weighted scoring**: Each collection has a configurable weight (e.g.,
  biomarker_reference at 0.12, clinical_evidence at 0.09) that influences how
  results are ranked.
- **Knowledge graph augmentation**: Before the LLM sees the evidence, the system
  enriches it with structured domain knowledge -- disease domain context, PGx
  gene context, PhenoAge biomarker context.
- **Citation scoring**: Each retrieved result is tagged as high, medium, or low
  relevance based on its similarity score (thresholds: 0.75 for high, 0.60 for
  medium).
- **Disease area filtering**: If the question mentions "diabetes," only the
  disease trajectory and genetic variant collections that have a disease_area
  field are filtered to diabetes-relevant results.

---

## Chapter 5: System Overview

### Architecture

The Precision Biomarker Intelligence Agent follows a four-phase pipeline:

```
Plan --> Analyze --> Search --> Synthesize
```

**Phase 1: Plan** -- The agent analyzes the incoming question using keyword
detection to identify disease areas (diabetes, cardiovascular, liver, thyroid,
iron, nutritional), relevant analysis modules (biological_age,
disease_trajectory, pharmacogenomics, genotype_adjustment), and decomposes
complex questions into sub-queries.

**Phase 2: Analyze** -- If a patient profile is provided, the agent runs four
analysis modules:

- `BiologicalAgeCalculator` -- PhenoAge (Levine 2018) and GrimAge surrogate
  (Lu 2019) estimation from routine blood biomarkers.
- `DiseaseTrajectoryAnalyzer` -- Pre-symptomatic risk assessment across 6+
  disease categories with estimated years to clinical onset.
- `PharmacogenomicMapper` -- Star allele to metabolizer phenotype mapping using
  CPIC Level 1A guidelines for 14 pharmacogenes.
- `GenotypeAdjuster` -- Genotype-adjusted reference ranges for biomarkers
  affected by genetic variants (MTHFR, APOE, PNPLA3, HFE, DIO2, VDR, FADS1).

Additionally, three safety modules run:

- `CriticalValueEngine` -- Checks biomarker values against critical thresholds
  requiring immediate clinical action.
- `DiscordanceDetector` -- Identifies cross-biomarker patterns that are
  clinically inconsistent and may indicate errors or rare conditions.
- `LabRangeInterpreter` -- Compares values against both standard and optimal
  ranges, flagging markers that are "normal" by lab standards but suboptimal by
  evidence-based criteria.

**Phase 3: Search** -- The RAG engine embeds the question using
BGE-small-en-v1.5, searches all 14 Milvus collections in parallel using
ThreadPoolExecutor, applies disease-area and year-range filters, deduplicates
results, and ranks them by weighted similarity score. Up to 30 merged results
are retained after deduplication.

**Phase 4: Synthesize** -- The agent builds an enhanced prompt that combines:

- Retrieved evidence with citation labels and relevance tags.
- Knowledge graph context (disease domains, PGx, PhenoAge, cross-modal links).
- Patient profile data (age, sex, biomarkers, genotypes, star alleles).
- Analysis module results (biological age, disease trajectories, PGx findings,
  adjusted ranges, critical alerts).

This prompt is sent to Claude with a specialized system prompt that instructs the
LLM to cite sources, specify units, provide genotype-specific interpretation,
highlight critical findings, and flag cross-modal triggers.

### The 14 Collections as Specialized Libraries

Think of the 14 Milvus collections as 14 specialized libraries, each containing
a different type of knowledge:

| # | Collection | What It Holds |
|---|-----------|---------------|
| 1 | biomarker_reference | Definitions, units, ranges for every tracked biomarker |
| 2 | biomarker_genetic_variants | SNPs and genes that modify biomarker levels |
| 3 | biomarker_pgx_rules | CPIC drug-gene dosing recommendations |
| 4 | biomarker_disease_trajectories | Disease progression stages with biomarker patterns |
| 5 | biomarker_clinical_evidence | Published clinical literature with PubMed IDs |
| 6 | biomarker_nutrition | Genotype-aware nutritional guidelines |
| 7 | biomarker_drug_interactions | Gene-drug interaction effects and alternatives |
| 8 | biomarker_aging_markers | PhenoAge/GrimAge clock coefficients and interpretation |
| 9 | biomarker_genotype_adjustments | Adjusted reference ranges by genotype |
| 10 | biomarker_monitoring | Condition-specific monitoring protocols |
| 11 | biomarker_critical_values | Threshold values requiring immediate clinical action |
| 12 | biomarker_discordance_rules | Cross-biomarker inconsistency patterns |
| 13 | biomarker_aj_carrier_screening | Ashkenazi Jewish carrier screening panel |
| 14 | genomic_evidence | Shared read-only collection from the Genomics Pipeline |

When a question arrives, all 14 libraries are searched simultaneously. The
results are merged, deduplicated, and ranked. The highest-scoring evidence from
each relevant collection is included in the LLM prompt.

---

## Chapter 6: Your First Query

This chapter walks you through running your first biomarker analysis, from
opening the UI to interpreting the results.

### Step 1: Open the UI

The Streamlit interface runs on port **8528**. Open your browser and navigate to:

```
http://localhost:8528
```

You will see the Biomarker Intelligence Agent interface with 8 tabs:

1. **Biomarker Analysis** -- Full patient analysis pipeline.
2. **Biological Age** -- PhenoAge calculator.
3. **Disease Risk** -- Disease trajectory analysis.
4. **PGx Profile** -- Pharmacogenomic drug interaction mapping.
5. **Evidence Explorer** -- RAG Q&A with collection filtering.
6. **Reports** -- PDF and FHIR R4 export.
7. **Patient 360** -- Unified cross-agent intelligence dashboard.
8. **Longitudinal** -- Biomarker trend tracking across multiple visits.

### Step 2: Load the Demo Patient

Click the **Biomarker Analysis** tab. You will see an option to load a sample
patient. Click it to load demo patient **HCLS-BIO-2026-00001**.

This demo patient comes pre-populated with:

- **Demographics**: Age 45, Male.
- **Biomarkers**: A full panel including albumin, creatinine, glucose, hs-CRP,
  lymphocyte percentage, MCV, RDW, alkaline phosphatase, WBC, HbA1c, ferritin,
  TSH, total cholesterol, LDL, HDL, triglycerides, ALT, AST, GGT, homocysteine,
  and vitamin D.
- **Genotypes**: Key SNPs including MTHFR C677T (rs1801133), APOE (rs429358),
  TCF7L2 (rs7903146), PNPLA3 (rs738409), HFE C282Y (rs1800562), and others.
- **Star Alleles**: CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, and DPYD.

### Step 3: Run the Analysis

Click the **Analyze** button. The system will:

1. Calculate biological age using PhenoAge (and GrimAge surrogate if plasma
   protein markers are available).
2. Run disease trajectory analysis for diabetes, cardiovascular, liver, thyroid,
   iron, and nutritional conditions.
3. Map star alleles to metabolizer phenotypes and drug recommendations.
4. Compute genotype-adjusted reference ranges.
5. Check for critical values and cross-biomarker discordances.
6. Search all 14 collections for evidence relevant to the patient's profile.
7. Synthesize a comprehensive report with citations.

### Step 4: Understanding the Results

The results page shows several sections:

**Biological Age Assessment**: You will see a chronological age (45) alongside a
biological age (e.g., 42.3), an age acceleration score (e.g., -2.7 years,
meaning aging slower than expected), and the top aging drivers -- the biomarkers
contributing most to your biological age calculation, with their direction
(protective or aging).

**Disease Trajectories**: Each disease category (diabetes, cardiovascular, liver,
etc.) shows a risk level (normal, low, moderate, high, critical), relevant
current biomarker values, genetic risk factors, estimated years to onset (if
applicable), and intervention recommendations.

**Pharmacogenomic Profile**: Each pharmacogene shows the star allele result, the
metabolizer phenotype (ultra-rapid, normal, intermediate, poor), and affected
drugs with CPIC-level dosing recommendations.

**Genotype-Adjusted Ranges**: Biomarkers whose reference ranges change based on
genotype are highlighted with both the standard and adjusted ranges. For example,
if the patient has MTHFR CT genotype, the homocysteine upper limit drops from
15 to 12 umol/L.

**Critical Alerts**: Any findings requiring immediate attention are displayed
prominently. These include severe age acceleration, critical disease risk
trajectories, PGx drug safety alerts (e.g., DPYD poor metabolizer with
fluoropyrimidine contraindication), and cross-biomarker discordances.

---

## Chapter 7: Understanding Collections

This chapter describes each of the 14 Milvus vector collections in detail,
including the fields stored in each one. Every collection also contains an
`embedding` field (384-dimensional float vector from BGE-small-en-v1.5) used for
similarity search, which is omitted from the field lists below for brevity.

### 1. biomarker_reference

The foundational collection. Every biomarker tracked by the system has an entry
here defining what it is, how it is measured, what the normal range is, and why
it matters.

| Field | Description |
|-------|-------------|
| name | Display name (e.g., "Albumin", "HbA1c") |
| category | Panel category: CBC, CMP, Lipids, Thyroid, Inflammation, Nutrients |
| unit | Measurement unit (e.g., g/dL, mg/dL, %, mIU/L) |
| standard_range | Population reference range (ref_range_min to ref_range_max) |
| optimal_range | Evidence-based optimal range (often narrower than standard) |
| description | Text chunk describing the biomarker for embedding |
| clinical_significance | Clinical interpretation and significance |

Additional fields include `epigenetic_clock` (PhenoAge or GrimAge coefficient
if applicable) and `genetic_modifiers` (comma-separated gene symbols that
modify this biomarker).

**Example**: Albumin has a standard range of 3.5-5.5 g/dL, but the PhenoAge
algorithm uses it as a protective factor (negative coefficient -0.0336 applied
to g/L) -- higher albumin is associated with younger biological age.

### 2. biomarker_genetic_variants

Stores genetic variants (SNPs) that affect biomarker levels. This is where the
system learns that MTHFR C677T reduces folate metabolism or that APOE e4
increases LDL cholesterol.

| Field | Description |
|-------|-------------|
| gene | Gene symbol (e.g., MTHFR, APOE, PNPLA3) |
| variant | Variant name (e.g., C677T, I148M) |
| rsid | dbSNP identifier (e.g., rs1801133) |
| effect | Effect on biomarker levels (e.g., "20-30% reduction in enzyme activity") |
| clinical_significance | Clinical meaning and risk implications |
| population_frequency | Allele frequency in major populations |

Additional fields include `risk_allele`, `protective_allele`, `mechanism`
(molecular explanation), and `disease_associations`.

### 3. biomarker_pgx_rules

Pharmacogenomic dosing rules from CPIC guidelines. Each record maps a gene +
star allele combination + drug to a specific dosing recommendation.

| Field | Description |
|-------|-------------|
| gene | Pharmacogene (e.g., CYP2D6, CYP2C19, DPYD) |
| drug | Drug name (e.g., codeine, clopidogrel, warfarin) |
| phenotype | Metabolizer classification (ultra_rapid, normal, intermediate, poor) |
| recommendation | CPIC dosing recommendation text |
| cpic_level | Evidence level (1A, 1B, 2A, 2B, 3) |
| source | CPIC guideline URL or PharmGKB entry |

Additional fields include `star_alleles` (the specific allele combination like
`*1/*2`) and `evidence_url`.

### 4. biomarker_disease_trajectories

Models disease progression over time. Each record represents a stage in a
disease trajectory -- from earliest detectable biomarker changes through
pre-clinical diagnosis to clinical disease.

| Field | Description |
|-------|-------------|
| disease | Disease category (diabetes, cardiovascular, liver, thyroid, iron, nutritional) |
| category | Disease classification |
| biomarker_pattern | JSON string of biomarker thresholds defining this stage |
| risk_score_formula | Algorithm used to calculate risk at this stage |
| stage | Stage name (e.g., pre-diabetes, early, advanced) |
| progression_rate | Typical rate of progression between stages |

Additional fields include `years_to_diagnosis`, `intervention_window`, and
`risk_reduction_pct` (potential risk reduction with intervention).

### 5. biomarker_clinical_evidence

Published clinical literature indexed for retrieval. Each record is a chunk from
a journal article, systematic review, or clinical guideline.

| Field | Description |
|-------|-------------|
| title | Publication title |
| text_chunk | Text passage for embedding and retrieval |
| source_type | Type of source (journal article, review, guideline, meta-analysis) |
| year | Publication year (enables temporal filtering) |
| biomarker | Primary biomarker discussed |
| disease | Disease area or specialty |
| evidence_level | Strength of evidence (I, II, III, IV) |

Additional fields include `pmid` (PubMed ID for citation linking) and `finding`
(key finding summary). The RAG engine can filter this collection by year range
(e.g., only evidence from 2020-2026).

### 6. biomarker_nutrition

Genotype-aware nutritional guidelines. These go beyond generic dietary advice
to account for how genetic variants affect nutrient metabolism.

| Field | Description |
|-------|-------------|
| nutrient | Nutrient name (e.g., Folate, Vitamin D, Omega-3) |
| biomarker | Related biomarker (e.g., homocysteine, 25-OH Vitamin D) |
| effect | How the genotype affects nutrient metabolism |
| recommended_intake | Genotype-adjusted intake recommendation |
| deficiency_threshold | Threshold below which deficiency is likely |
| food_sources | Dietary sources for this nutrient |

Additional fields include `genetic_context` (e.g., "MTHFR C677T heterozygous"),
`recommended_form` (e.g., "methylfolate" instead of "folic acid" for MTHFR
carriers), and `dose_range`.

### 7. biomarker_drug_interactions

Gene-drug interactions that affect biomarker levels. Distinct from PGx rules
(which focus on drug dosing), these records describe how a drug alters a
biomarker measurement.

| Field | Description |
|-------|-------------|
| drug | Drug name |
| biomarker | Affected biomarker |
| direction | Direction of change (increase, decrease) |
| magnitude | Size of the effect (mild, moderate, significant) |
| mechanism | How the drug affects the biomarker |
| clinical_significance | Whether the change requires clinical action |

Additional fields include `gene` (gene involved in the interaction),
`interaction_type` (substrate, inhibitor, inducer), `severity`, and
`alternative` (alternative drug recommendation).

### 8. biomarker_aging_markers

Stores the coefficients and interpretation guides for epigenetic aging clock
biomarkers. These are the markers used by the PhenoAge and GrimAge algorithms
to estimate biological age.

| Field | Description |
|-------|-------------|
| biomarker | Marker name (e.g., Albumin, Creatinine, Glucose) |
| age_coefficient | Coefficient weight in the clock algorithm |
| phenoage_weight | Specific PhenoAge coefficient value |
| optimal_range_by_age | How the optimal range shifts with chronological age |

Additional fields include `clock_type` (PhenoAge or GrimAge), `unit`, and
`interpretation` (clinical explanation of the marker's role in aging).

### 9. biomarker_genotype_adjustments

Stores the adjusted reference ranges for biomarkers based on specific genotypes.
When a patient carries a variant that shifts what "normal" means, this collection
provides the corrected range.

| Field | Description |
|-------|-------------|
| biomarker | Biomarker name (e.g., homocysteine, ALT, TSH) |
| genotype | Specific genotype (e.g., CT, TT, CG, GG) |
| ancestry | Population/ancestry context for the adjustment |
| adjustment_factor | Multiplier or offset applied to the standard range |
| reference | Published source for the adjustment |

Additional fields include `gene`, `rs_id`, `genotype_ref` (reference genotype),
`genotype_het` (heterozygous), `genotype_hom` (homozygous alternate),
`adjusted_min`, `adjusted_max`, and `rationale`.

### 10. biomarker_monitoring

Condition-specific monitoring protocols. Once a patient is identified as being
at risk, how often should they be tested, and what should trigger escalation?

| Field | Description |
|-------|-------------|
| biomarker | Biomarker to monitor |
| condition | Medical condition requiring monitoring |
| frequency | How often to test (e.g., "every 3 months") |
| alert_threshold | Value that triggers a clinical alert |
| escalation_protocol | What to do when the alert threshold is crossed |

Additional fields include `biomarker_panel` (comma-separated list of all markers
in the monitoring panel) and `trigger_values` (JSON string of threshold values).

### 11. biomarker_critical_values

Critical value thresholds that require immediate clinical action. These are the
"panic values" -- results so far outside normal that they represent an
immediate health risk.

| Field | Description |
|-------|-------------|
| biomarker | Biomarker name |
| low_critical | Below this value, immediate action is required |
| high_critical | Above this value, immediate action is required |
| unit | Measurement unit |
| action_required | Clinical action to take |

Additional fields include `loinc_code` (standardized lab code), `severity`
(critical/urgent/warning), `escalation_target` (who to notify), and
`cross_checks` (other biomarkers to verify before acting).

### 12. biomarker_discordance_rules

Patterns of cross-biomarker inconsistency that are clinically meaningful. When
two biomarkers that should move together diverge, it may indicate a specific
condition, lab error, or rare disease.

| Field | Description |
|-------|-------------|
| biomarker_a | First biomarker in the pair |
| biomarker_b | Second biomarker in the pair |
| expected_relationship | How these markers normally relate (e.g., "both elevated in inflammation") |
| discordance_threshold | How far apart they must be to trigger a flag |
| clinical_implication | What the discordance might mean clinically |

Additional fields include `name` (rule name), `condition` (triggering
condition), `differential_diagnosis`, `agent_handoff` (which other agent to
consult), and `priority` (high/medium/low).

### 13. biomarker_aj_carrier_screening

Ashkenazi Jewish carrier screening panel data. This collection covers genetic
conditions with elevated carrier frequency in the Ashkenazi Jewish population.

| Field | Description |
|-------|-------------|
| gene | Gene symbol (e.g., BRCA1, BRCA2, GBA, HEXA, FANCC) |
| variant | Specific pathogenic variant(s) |
| condition | Associated disease (e.g., Tay-Sachs, Gaucher, Fanconi Anemia) |
| carrier_frequency | Carrier frequency in AJ population |
| clinical_action | Recommended clinical follow-up |

Additional fields include `inheritance` (autosomal recessive, etc.),
`general_carrier_frequency`, `common_mutations`, `loinc_code`, `method`,
`clinical_significance`, `reproductive_implications`, and `compound_risks`
(e.g., GBA + APOE e4 interaction for Parkinson's risk).

### 14. genomic_evidence

A **read-only** shared collection populated by the Genomics Pipeline. Contains
VCF-derived variant data for the patient. The Biomarker Intelligence Agent
searches this collection but does not write to it.

This collection bridges the gap between raw genomic data (the Genomics Pipeline
output) and the biomarker interpretation layer. When the agent detects a
pharmacogenomic variant or disease-associated SNP in genomic_evidence, it links
that finding to the relevant PGx rule, disease trajectory, or genotype
adjustment.

---

## Chapter 8: The Knowledge Graph

Beyond the 14 searchable collections, the agent maintains a structured knowledge
graph in `src/knowledge.py` (1,326 lines). This is not stored in Milvus; it is
compiled into the application and injected into LLM prompts as additional
context when relevant topics are detected.

### 6 Disease Domains

The knowledge graph defines 7 disease domains (6 primary, plus additional
categories), each containing:

- **Key biomarkers**: The specific markers relevant to this disease.
- **Genetic modifiers**: Genes that modify disease risk or biomarker
  interpretation.
- **Intervention targets**: Actionable levers for risk reduction.
- **Clinical thresholds**: Genotype-specific cutoff values.

The primary domains are:

1. **Diabetes** -- HbA1c, fasting glucose, insulin, HOMA-IR. Genetic modifiers
   include TCF7L2 (the strongest type 2 diabetes risk gene). The knowledge graph
   stores genotype-specific thresholds: for example, TCF7L2 homozygous risk
   allele carriers have the HbA1c concern threshold lowered from 6.0% to 5.5%.

2. **Cardiovascular** -- LDL, HDL, triglycerides, ApoB, Lp(a), homocysteine.
   Genetic modifiers include APOE (lipid metabolism), PCSK9 (LDL receptor
   regulation), and MTHFR (homocysteine clearance). Lp(a) above 50 mg/dL
   is an independent cardiovascular risk factor, genetically determined and
   largely unresponsive to statins.

3. **Liver** -- ALT, AST, GGT, FIB-4 index. Genetic modifier: PNPLA3 I148M
   (rs738409), which dramatically increases NAFLD/NASH risk. The knowledge graph
   adjusts ALT upper limits by genotype: 56 U/L for CC, 45 for CG, 35 for GG.

4. **Thyroid** -- TSH, free T3, free T4. Genetic modifier: DIO2 Thr92Ala
   (rs225014), affecting T4-to-T3 conversion. TSH upper limits are adjusted
   by genotype and age.

5. **Iron** -- Ferritin, transferrin saturation, serum iron. Genetic modifier:
   HFE C282Y/H63D for hereditary hemochromatosis. TMPRSS6 variants affect
   iron-refractory anemia.

6. **Nutritional** -- Folate, B12, vitamin D, omega-3 index. Multiple genetic
   modifiers: MTHFR (folate), FADS1 (omega-3 conversion), VDR (vitamin D
   receptor), BCMO1 (beta-carotene to vitamin A), FUT2 (B12 absorption).

### PhenoAge Coefficients

The knowledge graph stores the complete PhenoAge algorithm context from Levine
et al. 2018 (PMID:29676998). The 9 biomarkers used are albumin, creatinine,
glucose, hs-CRP (log-transformed), lymphocyte percentage, MCV, RDW, alkaline
phosphatase, and WBC. Each has a coefficient that represents its contribution
to biological age estimation.

Key insights encoded in the knowledge graph:

- Albumin is protective (coefficient -0.0336 on g/L) -- higher albumin
  associates with younger biological age.
- RDW has the largest positive coefficient (0.3306 on %) -- elevated RDW is a
  strong aging signal.
- The algorithm also includes chronological age as an input (coefficient 0.0804).
- PhenoAge uses Gompertz mortality modeling to convert a linear predictor into
  a biological age estimate.

### PGx Knowledge

The knowledge graph stores structured data for 14 pharmacogenes:

1. CYP2D6 -- codeine, tramadol, tamoxifen, SSRIs
2. CYP2C19 -- clopidogrel, PPIs, voriconazole
3. CYP2C9 -- warfarin, NSAIDs, phenytoin
4. VKORC1 -- warfarin dose sensitivity
5. SLCO1B1 -- statin myopathy risk (simvastatin)
6. DPYD -- fluoropyrimidine toxicity (5-FU, capecitabine)
7. CYP3A5 -- tacrolimus dosing
8. TPMT -- thiopurine dosing (azathioprine, 6-MP)
9. MTHFR -- folate metabolism (informational, not CPIC Level 1A)

Each gene entry includes key drugs, CPIC guideline versions, and metabolizer
phenotype implications. Additional genes tracked in the pharmacogenomics module
include HLA-B*57:01 (abacavir hypersensitivity), G6PD, HLA-B*58:01, UGT1A1,
and NUDT15.

### Ancestry Adjustments

The knowledge graph includes ancestry-aware context. Allele frequencies and
disease risk differ significantly across populations. For example:

- MTHFR C677T homozygosity: ~10% in Europeans, ~25% in Hispanics.
- CYP2D6 poor metabolizer: ~6-10% in Europeans, ~1% in East Asians.
- HFE C282Y homozygosity: ~0.5% in Northern Europeans, extremely rare in
  Africans and East Asians.
- APOL1 high-risk genotype: ~13% in African Americans, very rare in Europeans.

The agent uses the patient's self-reported ancestry (if provided) to select
appropriate population-specific reference ranges and risk calculations.

### Cross-Modal Links

The knowledge graph defines 8 cross-modal links that trigger communication with
other HCLS AI Factory agents:

- Elevated Lp(a) or iron overload triggers the Imaging Intelligence Agent for
  cardiac calcium scoring or liver MRI.
- PGx drug safety alerts (e.g., DPYD poor metabolizer on fluoropyrimidines)
  trigger the Oncology or CAR-T Intelligence Agent.
- Novel or uncertain variants trigger re-analysis by the Genomics Pipeline.

---

## Chapter 9: Pharmacogenomics Explained

Pharmacogenomics (PGx) is the study of how genetic variation affects drug
response. This chapter explains the concepts you need to understand the agent's
pharmacogenomic module, implemented in `src/pharmacogenomics.py` (1,503 lines).

### Star Alleles

Pharmacogenes like CYP2D6 are described using **star allele** nomenclature. Each
star allele represents a specific combination of genetic variants in that gene:

- **\*1** -- The reference (normal function) allele. Most people carry at least
  one copy.
- **\*2** -- A variant with normal or near-normal function.
- **\*4** -- A common non-functional allele (the most frequent CYP2D6 null
  allele in Europeans).
- **\*5** -- Gene deletion (no protein produced at all).
- **\*17** -- Reduced function allele (more common in Africans).
- **\*xN** -- Gene duplication (extra copies, leading to ultra-rapid
  metabolism).

A patient's result is written as a diplotype: two alleles separated by a slash.
For example, `CYP2D6 *1/*4` means one normal allele and one non-functional
allele.

### The 14 Pharmacogenes

The agent maps the following 14 genes to drug recommendations:

**1. CYP2D6** (Cytochrome P450 2D6) -- Metabolizes ~25% of all drugs. Key drugs:
codeine (converted to morphine), tramadol, tamoxifen (converted to endoxifen),
SSRIs (paroxetine, fluoxetine). Poor metabolizers get no pain relief from
codeine. Ultra-rapid metabolizers may overdose.

**2. CYP2C19** -- Metabolizes clopidogrel (the blood thinner Plavix). Poor
metabolizers cannot activate clopidogrel, leading to treatment failure and
potential stent thrombosis. Alternative: prasugrel or ticagrelor. Also affects
PPIs (omeprazole) and voriconazole.

**3. CYP2C9** -- Metabolizes warfarin, NSAIDs, and phenytoin. Poor metabolizers
accumulate warfarin, increasing bleeding risk. Dose reduction required.

**4. VKORC1** -- Not an enzyme but warfarin's drug target. The -1639G>A variant
(rs9923231) reduces VKORC1 expression, making patients sensitive to warfarin.
Combined with CYP2C9 genotype, it determines 40-50% of warfarin dose
variability.

**5. SLCO1B1** -- Encodes a hepatic transporter. The *5 variant (rs4149056)
reduces statin uptake into the liver, increasing blood levels and myopathy risk.
CPIC recommends avoiding simvastatin >20 mg in carriers.

**6. DPYD** -- Metabolizes fluoropyrimidines (5-FU, capecitabine). This is the
highest-stakes PGx interaction: DPYD poor metabolizers given standard-dose 5-FU
can develop life-threatening toxicity (severe mucositis, neutropenia, death).
Pre-treatment DPYD testing is now mandatory in many institutions.

**7. CYP3A5** -- Affects tacrolimus dosing in transplant patients. Expressers
(carriers of CYP3A5 *1) metabolize tacrolimus faster and require higher doses
to maintain therapeutic levels.

**8. TPMT** -- Metabolizes thiopurines (azathioprine, 6-mercaptopurine) used in
autoimmune diseases and leukemia. Poor metabolizers risk severe
myelosuppression. Dose reduction to 10% of standard is recommended.

**9. MTHFR** -- Encodes methylenetetrahydrofolate reductase. C677T and A1298C
variants reduce enzyme activity, affecting folate metabolism and homocysteine
levels. While not a formal CPIC Level 1A gene, it has significant implications
for methotrexate response and nutritional supplementation.

### Metabolizer Phenotypes

Star alleles are translated into four metabolizer phenotypes:

- **Ultra-Rapid Metabolizer (UM)**: Converts drug too fast. Risk of toxicity
  from active metabolites (codeine) or treatment failure from rapid clearance.
- **Normal Metabolizer (NM)**: Standard drug response. No dose adjustment
  needed.
- **Intermediate Metabolizer (IM)**: Reduced but not absent enzyme activity.
  May need moderate dose adjustment.
- **Poor Metabolizer (PM)**: Little to no enzyme activity. Drug accumulates
  (risk of toxicity) or prodrug is not activated (treatment failure).

### CPIC Levels

The Clinical Pharmacogenetics Implementation Consortium assigns evidence levels:

- **Level 1A**: Strong evidence, prescribing action recommended. Published CPIC
  guideline exists with strong evidence and strong recommendation.
- **Level 1B**: Strong evidence, action recommended. Published guideline with
  strong evidence but moderate recommendation.
- **Level 2A**: Moderate evidence, action suggested. Guideline exists with
  moderate evidence.
- **Level 2B**: Weak evidence, action may be considered.
- **Level 3**: Annotation only. Evidence exists but is insufficient for a
  clinical recommendation.

The agent's pharmacogenomic module implements CPIC Level 1A guidelines for
CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, DPYD, CYP3A5, TPMT, and others.
Guideline versions are tracked in `CPIC_GUIDELINE_VERSIONS` with PMIDs for
traceability.

---

## Chapter 10: Setting Up Locally

### Prerequisites

- **Docker** and **Docker Compose** installed.
- At least **8 GB RAM** available for Docker.
- An **Anthropic API key** for LLM queries (optional -- analysis modules work
  without it, but RAG queries require it).

### Environment Variables

All configuration uses the `BIOMARKER_` prefix. Key variables:

```bash
# Required for LLM queries
export BIOMARKER_ANTHROPIC_API_KEY="sk-ant-..."

# Milvus connection (defaults shown)
export BIOMARKER_MILVUS_HOST="localhost"
export BIOMARKER_MILVUS_PORT=19530

# Service ports (defaults shown)
export BIOMARKER_API_PORT=8529        # FastAPI REST API
export BIOMARKER_STREAMLIT_PORT=8528  # Streamlit UI

# Optional: API key for authentication (empty = no auth)
export BIOMARKER_API_KEY=""
```

The settings are managed by Pydantic `BaseSettings` in `config/settings.py`
(140 lines), which reads from environment variables and `.env` files
automatically.

### Ports

The agent uses three ports:

| Port | Service | Protocol |
|------|---------|----------|
| **8528** | Streamlit UI | HTTP (browser) |
| **8529** | FastAPI REST API | HTTP (JSON) |
| **19530** | Milvus vector database | gRPC |

Milvus also requires etcd (port 2379) and MinIO (ports 9000/9001) for its
internal operations. These are managed by Docker Compose.

### Starting with Docker Compose

```bash
# From the project root
cd ai_agent_adds/precision_biomarker_agent

# Start the full stack (Milvus + API + UI)
docker compose up -d

# Verify services are running
curl http://localhost:8529/health
```

### Seeding the Collections

On first startup, the collections need to be seeded with reference data.
The seeding scripts in `scripts/` populate all 13 writable collections with
biomarker definitions, genetic variants, PGx rules, disease trajectories,
clinical evidence, and more.

```bash
# Seed all collections
python scripts/seed_collections.py

# Verify seeding
curl http://localhost:8529/collections
```

The response will show each collection name and its record count.

### Running Without Docker

For development, you can run the services directly:

```bash
# Install dependencies
pip install -r requirements.txt

# Start the API server
uvicorn api.main:app --host 0.0.0.0 --port 8529 --reload

# In a separate terminal, start the UI
streamlit run app/biomarker_ui.py --server.port 8528
```

This requires Milvus to be running separately (either via Docker or a managed
instance).

---

## Chapter 11: Exploring the API

The FastAPI server at port 8529 provides 19+ endpoints organized into four
groups: status, analysis, reports, and events. Interactive API documentation is
available at `http://localhost:8529/docs` (Swagger UI).

### Status Endpoints

**GET /health** -- Service health with collection and vector counts.

```bash
curl http://localhost:8529/health
```

```json
{
  "status": "healthy",
  "collections": 14,
  "total_vectors": 2847,
  "agent_ready": true
}
```

**GET /collections** -- Collection names and record counts.

```bash
curl http://localhost:8529/collections
```

**GET /knowledge/stats** -- Knowledge graph statistics.

```bash
curl http://localhost:8529/knowledge/stats
```

```json
{
  "disease_domains": 7,
  "total_biomarkers": 42,
  "total_genetic_modifiers": 18,
  "pharmacogenes": 14,
  "pgx_drug_interactions": 27,
  "phenoage_markers": 9,
  "cross_modal_links": 8
}
```

**GET /metrics** -- Prometheus-compatible metrics.

```bash
curl http://localhost:8529/metrics
```

### Analysis Endpoints

**POST /v1/analyze** -- Full patient analysis using all modules.

```bash
curl -X POST http://localhost:8529/v1/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "age": 45,
    "sex": "M",
    "biomarkers": {
      "albumin": 4.2,
      "creatinine": 0.9,
      "glucose": 95,
      "hs_crp": 1.2,
      "lymphocyte_pct": 30,
      "mcv": 88,
      "rdw": 13.5,
      "alkaline_phosphatase": 65,
      "wbc": 6.5,
      "HbA1c": 5.7,
      "ferritin": 180,
      "tsh": 2.5
    },
    "genotypes": {
      "rs1801133": "CT",
      "rs7903146": "CT"
    },
    "star_alleles": {
      "CYP2D6": "*1/*2",
      "CYP2C19": "*1/*1"
    }
  }'
```

**POST /v1/biological-age** -- Biological age calculation only.

```bash
curl -X POST http://localhost:8529/v1/biological-age \
  -H "Content-Type: application/json" \
  -d '{
    "age": 45,
    "biomarkers": {
      "albumin": 4.2,
      "creatinine": 0.9,
      "glucose": 95,
      "hs_crp": 1.2,
      "lymphocyte_pct": 30,
      "mcv": 88,
      "rdw": 13.5,
      "alkaline_phosphatase": 65,
      "wbc": 6.5
    }
  }'
```

**POST /v1/disease-risk** -- Disease trajectory analysis.

```bash
curl -X POST http://localhost:8529/v1/disease-risk \
  -H "Content-Type: application/json" \
  -d '{
    "age": 45,
    "sex": "M",
    "biomarkers": {"HbA1c": 5.8, "fasting_glucose": 105, "triglycerides": 180},
    "genotypes": {"rs7903146": "CT"}
  }'
```

**POST /v1/pgx** -- Pharmacogenomic mapping.

```bash
curl -X POST http://localhost:8529/v1/pgx \
  -H "Content-Type: application/json" \
  -d '{
    "star_alleles": {
      "CYP2D6": "*4/*4",
      "CYP2C19": "*1/*2",
      "DPYD": "*1/*1"
    },
    "genotypes": {
      "rs9923231": "AG"
    }
  }'
```

**POST /v1/query** -- RAG Q&A query with optional patient context.

```bash
curl -X POST http://localhost:8529/v1/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What does an HbA1c of 5.8% mean for someone with TCF7L2 CT genotype?",
    "year_min": 2020
  }'
```

### Report Endpoints

**POST /v1/report/generate** -- Generate a full patient report.

```bash
curl -X POST http://localhost:8529/v1/report/generate \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "age": 45,
    "sex": "M",
    "biomarkers": {"albumin": 4.2, "HbA1c": 5.7}
  }'
```

**GET /v1/report/{report_id}** -- Retrieve a generated report.

**GET /v1/report/{report_id}/pdf** -- Download report as PDF.

**GET /v1/report/{report_id}/fhir** -- Export as FHIR R4 DiagnosticReport.

### Event Endpoints

**POST /v1/events/inbound** -- Receive a cross-modal event from another agent.

```bash
curl -X POST http://localhost:8529/v1/events/inbound \
  -H "Content-Type: application/json" \
  -d '{
    "source_agent": "imaging_intelligence_agent",
    "event_type": "imaging_finding",
    "payload": {"finding": "hepatic steatosis on ultrasound"}
  }'
```

**GET /v1/events/outbound** -- List outbound alerts sent by this agent.

**GET /v1/events/inbound** -- List received inbound events.

### Authentication

If `BIOMARKER_API_KEY` is set, all endpoints (except `/health` and `/metrics`)
require the key in the `X-API-Key` header:

```bash
curl -H "X-API-Key: your-secret-key" http://localhost:8529/v1/query ...
```

---

## Chapter 12: Understanding the Codebase

### Project Structure

```
precision_biomarker_agent/
|-- api/
|   |-- main.py                    (465 lines)  FastAPI app, lifespan, CORS, middleware
|   |-- routes/
|       |-- analysis.py            (495 lines)  /v1/analyze, biological-age, disease-risk, pgx, query
|       |-- reports.py             (296 lines)  /v1/report/generate, PDF, FHIR export
|       |-- events.py             (326 lines)  /v1/events/inbound, outbound
|
|-- app/
|   |-- biomarker_ui.py           (1,863 lines)  Streamlit UI with 8 tabs
|
|-- config/
|   |-- settings.py                (140 lines)  Pydantic BaseSettings, env vars, weights
|
|-- src/
|   |-- agent.py                   (610 lines)  PrecisionBiomarkerAgent: plan-analyze-search-synthesize
|   |-- rag_engine.py              (573 lines)  BiomarkerRAGEngine: multi-collection RAG
|   |-- collections.py            (1,391 lines)  BiomarkerCollectionManager: 14 Milvus schemas
|   |-- knowledge.py              (1,326 lines)  Knowledge graph: domains, PGx, PhenoAge, cross-modal
|   |-- models.py                  (786 lines)  Pydantic models for all collections and analysis
|   |-- pharmacogenomics.py       (1,503 lines)  PharmacogenomicMapper: star alleles to drug recs
|   |-- disease_trajectory.py     (1,421 lines)  DiseaseTrajectoryAnalyzer: 6+ disease categories
|   |-- genotype_adjustment.py    (1,225 lines)  GenotypeAdjuster: genotype-adjusted ranges
|   |-- biological_age.py          (408 lines)  BiologicalAgeCalculator: PhenoAge + GrimAge
|   |-- export.py                 (1,392 lines)  PDF and FHIR R4 export
|   |-- report_generator.py        (993 lines)  12-section clinical report generation
|   |-- critical_values.py         (179 lines)  CriticalValueEngine: panic value detection
|   |-- discordance_detector.py    (299 lines)  DiscordanceDetector: cross-biomarker patterns
|   |-- lab_range_interpreter.py   (221 lines)  Standard vs optimal range comparison
|   |-- audit.py                    (83 lines)  Audit logging for patient data access
|   |-- translation.py             (217 lines)  Multi-language report translation
|
|-- data/                          Reference data, cache, seeding files
|-- scripts/                       Seeding and utility scripts
|-- tests/                         Test suite
|-- docs/                          Documentation (including this guide)
|-- docker-compose.yml             Full stack: Milvus + etcd + MinIO + API + UI
|-- Dockerfile                     Container image definition
|-- requirements.txt               Python dependencies
|-- pyproject.toml                 Project metadata
```

### File-by-File Walkthrough

**config/settings.py** (140 lines) -- The central configuration file. Uses
Pydantic `BaseSettings` with the `BIOMARKER_` env prefix. Defines Milvus
connection parameters, embedding model (`BAAI/bge-small-en-v1.5`, 384
dimensions), LLM settings (Anthropic Claude), 14 collection weights that sum
to ~1.0, search parameters (top_k=5 per collection, score threshold=0.4),
citation thresholds (0.75 high, 0.60 medium), CORS origins, request size
limits, and authentication settings.

**src/models.py** (786 lines) -- Pydantic data models for every entity in the
system. Defines 7 enums (RiskLevel, ClockType, DiseaseCategory,
MetabolizerPhenotype, CPICLevel, Zygosity, plus source types). Contains models
for all 13 collection record types (BiomarkerReference, GeneticVariant, PGxRule,
DiseaseTrajectory, ClinicalEvidence, NutritionGuideline, DrugInteraction,
AgingMarker, GenotypeAdjustment, MonitoringProtocol, CriticalValue,
DiscordanceRule, AJCarrierScreeningEntry). Each model has a `to_embedding_text()`
method for generating the text fed to BGE-small. Also defines PatientProfile
(input), BiologicalAgeResult, DiseaseTrajectoryResult, PGxResult,
GenotypeAdjustmentResult, SearchHit, CrossCollectionResult, AgentQuery,
AgentResponse, AnalysisResult, PatientHistory (longitudinal tracking), and
WearableData (wearable device integration).

**src/collections.py** (1,391 lines) -- The Milvus collection manager. Defines
the schema (field types, dimensions, index parameters) for all 14 collections
using pymilvus. Implements `connect()`, `disconnect()`, `ensure_collections()`
(creates collections if they do not exist), `search_all()` (parallel search
across all collections using ThreadPoolExecutor), `insert()`, `delete()`,
`get_collection_stats()`, and bulk seeding methods. Each collection uses
COSINE similarity with IVF_FLAT index type.

**src/knowledge.py** (1,326 lines) -- The structured knowledge graph. Contains
`BIOMARKER_DOMAINS` (7 disease domains with biomarkers, genetic modifiers,
and intervention targets), `PHENOAGE_KNOWLEDGE` (PhenoAge clock biomarker
descriptions and coefficients), `PGX_KNOWLEDGE` (14 pharmacogenes with key
drugs and CPIC guidance), `CROSS_MODAL_LINKS` (8 links to other agents),
`GENOTYPE_THRESHOLDS` (shared clinical thresholds for cross-module
consistency), `AGE_SEX_REFERENCE_RANGES` (age- and sex-stratified ranges),
and `BIOMARKER_PLAUSIBLE_RANGES` (input validation). Helper functions:
`get_domain_context()`, `get_pgx_context()`, `get_biomarker_context()`,
`get_knowledge_stats()`.

**src/rag_engine.py** (573 lines) -- The multi-collection RAG engine.
`BiomarkerRAGEngine` orchestrates: embedding queries with BGE instruction
prefix, parallel search across all 14 collections, disease-area detection and
filtering, year-range filtering for clinical evidence, result deduplication
and ranking (capped at 30), knowledge graph augmentation, citation formatting
(with PubMed links for clinical evidence), prompt construction with evidence
sections and relevance tags, and LLM generation (both synchronous and
streaming). Also implements `find_related()` for cross-collection entity
linking ("show me everything about MTHFR").

**src/agent.py** (610 lines) -- The `PrecisionBiomarkerAgent` class. Implements
the plan-analyze-search-synthesize pipeline. `search_plan()` analyzes questions
to identify disease areas and relevant modules. `analyze_patient()` runs all
4 analysis modules plus critical value checking, discordance detection, lab
range interpretation, and age-stratified adjustments. `run()` orchestrates the
full pipeline: plan, analyze (if profile provided), search via RAG, evaluate
evidence quality, run sub-queries if evidence is insufficient, build enhanced
prompt with analysis results, and generate LLM answer.

**src/pharmacogenomics.py** (1,503 lines) -- The `PharmacogenomicMapper`. Maps
star alleles to metabolizer phenotypes for all 14 pharmacogenes using CPIC
guidelines. Contains the complete star allele to function mapping for CYP2D6,
CYP2C19, CYP2C9, VKORC1, SLCO1B1, DPYD, CYP3A5, TPMT, and MTHFR. Tracks
CPIC guideline versions with PMIDs. `map_all()` processes a patient's complete
star allele and genotype data and returns phenotypes, affected drugs, and dosing
recommendations.

**src/disease_trajectory.py** (1,421 lines) -- The `DiseaseTrajectoryAnalyzer`.
Analyzes pre-symptomatic disease risk across 6+ categories. Uses biomarker
patterns, genetic risk factors, and demographic data to estimate disease stage,
risk level, years to clinical onset, and intervention recommendations. Includes
detailed trajectory models for type 2 diabetes, cardiovascular disease,
NAFLD/NASH, thyroid disorders, iron overload, and nutritional deficiencies.

**src/biological_age.py** (408 lines) -- The `BiologicalAgeCalculator`.
Implements PhenoAge (Levine et al. 2018) using Gompertz mortality modeling with
9 blood biomarkers. Handles unit conversion (US clinical units to SI for
coefficient application). Also implements GrimAge surrogate estimation from
plasma protein markers (GDF-15, Cystatin C, leptin, PAI-1, TIMP-1, ADM) with
confidence scoring based on marker coverage. Computes 95% confidence intervals
and risk classification with confidence qualifiers.

**src/genotype_adjustment.py** (1,225 lines) -- The `GenotypeAdjuster`.
Calculates genotype-adjusted reference ranges for biomarkers modified by genetic
variants. Covers MTHFR (homocysteine), APOE (lipids), PNPLA3 (ALT), HFE
(ferritin), DIO2 (TSH), VDR (vitamin D), FADS1 (omega-3), TCF7L2 (glucose/
HbA1c), and APOL1 (eGFR). Also applies age- and sex-stratified adjustments
using ranges from NHANES III and Harrison's Principles of Internal Medicine.

**src/export.py** (1,392 lines) -- Export functionality. Generates PDF reports
using ReportLab with clinical formatting (headers, tables, charts). Exports
FHIR R4 DiagnosticReport resources with proper coding (LOINC, SNOMED CT) for
interoperability with electronic health record systems.

**src/report_generator.py** (993 lines) -- Generates comprehensive 12-section
clinical reports combining biological age assessment, disease trajectories,
pharmacogenomic profile, genotype-adjusted ranges, critical alerts, nutritional
recommendations, monitoring protocols, and cross-modal triggers.

**src/critical_values.py** (179 lines) -- `CriticalValueEngine`. Checks
biomarker values against critical (panic) thresholds and generates alerts for
values requiring immediate clinical action.

**src/discordance_detector.py** (299 lines) -- `DiscordanceDetector`. Identifies
cross-biomarker inconsistencies (e.g., elevated ferritin with low iron
saturation, which may suggest inflammation rather than iron overload).

**src/lab_range_interpreter.py** (221 lines) -- `LabRangeInterpreter`. Compares
values against both standard lab ranges and evidence-based optimal ranges,
identifying biomarkers that are "normal" by lab standards but suboptimal for
health optimization.

**src/audit.py** (83 lines) -- Audit logging for patient data access. Records
who accessed what patient data, when, and what operations were performed.

**src/translation.py** (217 lines) -- Multi-language translation support for
generated reports.

**api/main.py** (465 lines) -- FastAPI application. Manages the application
lifespan (startup initialization of Milvus, embedder, LLM, and all analysis
modules; graceful shutdown). Configures CORS, API key authentication, and
request size limiting middleware. Includes status endpoints (/health,
/collections, /knowledge/stats, /metrics). Includes route modules for
analysis, reports, and events.

**app/biomarker_ui.py** (1,863 lines) -- Streamlit interface with 8 tabs. The
most user-facing part of the system. Provides forms for patient data entry,
sample patient loading, result visualization with charts and tables, and
report download.

---

## Chapter 13: Next Steps

Now that you understand the foundations, here are some paths for deeper
exploration:

### For Developers

1. **Read the API docs**: Open `http://localhost:8529/docs` and try every
   endpoint with the interactive Swagger UI.
2. **Trace a query**: Put a breakpoint in `agent.py:run()` and follow a question
   through plan, analyze, search, and synthesize. Watch how evidence flows from
   Milvus through the RAG engine to the LLM.
3. **Add a new collection**: Follow the pattern in `collections.py` to add a new
   domain collection. Create its Pydantic model in `models.py`, add its weight
   in `settings.py`, and register it in `rag_engine.py:COLLECTION_CONFIG`.
4. **Write tests**: The `tests/` directory contains the test suite. The critical
   analysis modules (biological_age, pharmacogenomics, disease_trajectory,
   genotype_adjustment) are pure computation with no external dependencies,
   making them ideal for unit testing.

### For Data Scientists

1. **Experiment with weights**: Adjust the 14 collection weights in
   `config/settings.py` and observe how they change the ranking of evidence in
   RAG responses.
2. **Analyze the PhenoAge algorithm**: Read `src/biological_age.py` alongside
   the Levine et al. 2018 paper (PMID:29676998). Verify the coefficients. Try
   modifying the Gompertz parameters and observe the effect on biological age
   estimates.
3. **Build a validation cohort**: Use the biological age calculator on known
   patient panels and compare against published PhenoAge results.

### For Clinicians

1. **Run the demo patient**: Walk through Chapter 6 with the demo patient and
   verify that the clinical interpretations match your medical knowledge.
2. **Test edge cases**: Enter patients with extreme biomarker values, rare
   genotype combinations, or complex multi-drug regimens. See how the system
   handles them.
3. **Review the PGx recommendations**: Compare the agent's CPIC recommendations
   against the current guidelines at cpicpgx.org.

### For Everyone

1. **Read the knowledge graph**: Open `src/knowledge.py` and read through the
   disease domain definitions. This is the most readable file in the codebase
   and provides a medical education in itself.
2. **Explore cross-modal links**: Understand how this agent connects to the
   other four HCLS AI Factory agents (Imaging, Oncology, CAR-T, Autoimmune).
3. **Contribute**: File issues for incorrect clinical data, missing drug-gene
   pairs, or collection schema improvements.

---

## Chapter 14: Glossary

**Allele** -- One of the possible versions of a gene. Most genes have two
alleles (one from each parent). Variants between alleles can affect protein
function.

**ApoB** -- Apolipoprotein B. A protein component of LDL particles. Considered a
more accurate predictor of cardiovascular risk than LDL-C because it counts
actual atherogenic particles.

**APOE** -- Apolipoprotein E gene. Three common variants (e2, e3, e4) affect
lipid metabolism and Alzheimer's disease risk. E4 carriers have elevated LDL
and increased cardiovascular and neurological risk.

**Biological Age** -- An estimate of physiological age based on biomarker values,
as opposed to chronological (calendar) age. A person who is 50 years old
chronologically may have a biological age of 45 (aging slower) or 55 (aging
faster).

**BGE-small-en-v1.5** -- BAAI General Embedding model. The embedding model used
by this agent to convert text into 384-dimensional vectors for similarity
search.

**Biomarker** -- Any measurable substance or characteristic in the body that
indicates a biological state, risk factor, or disease process.

**CBC** -- Complete Blood Count. A standard blood test measuring WBC, RBC,
hemoglobin, hematocrit, MCV, RDW, platelets, and differential counts.

**CMP** -- Comprehensive Metabolic Panel. A blood test measuring glucose,
electrolytes, kidney function (creatinine, BUN), liver enzymes (ALT, AST,
alkaline phosphatase), and proteins (albumin, total protein).

**CPIC** -- Clinical Pharmacogenetics Implementation Consortium. An organization
that publishes evidence-based guidelines for translating genetic test results
into actionable prescribing decisions.

**CRP (hs-CRP)** -- C-Reactive Protein (high-sensitivity). An inflammatory
marker produced by the liver. Elevated levels indicate systemic inflammation
and increased cardiovascular risk.

**Cross-Modal Trigger** -- A signal sent from one HCLS AI Factory agent to
another when a finding in one domain warrants investigation in another (e.g.,
elevated Lp(a) triggering cardiac imaging).

**CYP2D6** -- Cytochrome P450 2D6. A liver enzyme that metabolizes
approximately 25% of all prescription drugs. Genetic variation causes wide
inter-individual differences in drug metabolism.

**Diplotype** -- The combination of two alleles at a gene locus (one from each
parent), written as allele1/allele2 (e.g., CYP2D6 *1/*4).

**Discordance** -- A pattern where two biomarkers that normally move together
show an unexpected divergence, potentially indicating a specific clinical
condition, lab error, or rare disease.

**Disease Trajectory** -- The projected course of a disease from earliest
detectable biomarker changes through pre-clinical stages to clinical diagnosis,
enabling early intervention.

**DPYD** -- Dihydropyrimidine dehydrogenase gene. Metabolizes fluoropyrimidine
chemotherapy drugs. Poor metabolizers face life-threatening toxicity from
standard doses.

**Embedding** -- A mathematical representation of text as a vector (list of
numbers) in high-dimensional space, where similar meanings produce similar
vectors.

**FADS1** -- Fatty Acid Desaturase 1. Affects conversion of short-chain to
long-chain omega-3 fatty acids. Variants may increase dietary omega-3
requirements.

**FHIR R4** -- Fast Healthcare Interoperability Resources, Release 4. A standard
for exchanging healthcare information electronically. The agent can export
reports in FHIR R4 DiagnosticReport format.

**FIB-4** -- Fibrosis-4 Index. A calculated score using age, AST, ALT, and
platelet count to estimate liver fibrosis severity.

**Genotype** -- The specific combination of alleles at a given genetic locus.
Written as two letters (e.g., CT for heterozygous MTHFR C677T).

**GrimAge** -- An epigenetic clock (Lu et al. 2019) that predicts lifespan and
healthspan using DNA methylation patterns. This agent uses a simplified
surrogate based on plasma protein markers.

**HbA1c** -- Glycated Hemoglobin. Reflects average blood glucose over the
preceding 2-3 months. Used to diagnose and monitor diabetes.

**HFE** -- The hemochromatosis gene. C282Y and H63D variants cause hereditary
hemochromatosis (iron overload disorder).

**HOMA-IR** -- Homeostatic Model Assessment of Insulin Resistance. Calculated
from fasting glucose and insulin to estimate insulin resistance.

**Homocysteine** -- An amino acid in the blood. Elevated levels (often due to
MTHFR variants) are associated with cardiovascular disease, stroke, and
cognitive decline.

**IVF_FLAT** -- An index type used by Milvus for approximate nearest neighbor
search. Provides good recall with moderate memory usage.

**Knowledge Graph** -- A structured representation of domain knowledge
(disease domains, gene-drug relationships, biomarker coefficients) that
supplements the vector-searchable collections.

**LDL** -- Low-Density Lipoprotein cholesterol. The primary target for
cardiovascular risk reduction. Elevated LDL is a major risk factor for
atherosclerosis.

**LOINC** -- Logical Observation Identifiers Names and Codes. A standardized
coding system for laboratory tests and clinical observations.

**Lp(a)** -- Lipoprotein(a). A genetically determined cardiovascular risk factor
largely unresponsive to statins. Values above 50 mg/dL indicate elevated risk.

**MCV** -- Mean Corpuscular Volume. The average size of red blood cells.
Elevated MCV (macrocytosis) can indicate B12 or folate deficiency. Used in
the PhenoAge algorithm.

**Metabolizer Phenotype** -- Classification of enzyme activity based on
genetics: ultra-rapid, normal, intermediate, or poor. Determines how fast a
person processes a drug.

**Milvus** -- An open-source vector database designed for similarity search on
embedding vectors. The storage backend for all 14 collections.

**MTHFR** -- Methylenetetrahydrofolate reductase. C677T and A1298C variants
reduce enzyme activity, affecting folate metabolism and homocysteine levels.

**PGx (Pharmacogenomics)** -- The study of how genetic variation affects
individual responses to drugs. Used to optimize drug selection and dosing.

**PhenoAge** -- An epigenetic aging clock (Levine et al. 2018) that estimates
biological age from 9 routine blood biomarkers using Gompertz mortality
modeling.

**PNPLA3** -- Patatin-like phospholipase domain-containing protein 3. The I148M
variant (rs738409) dramatically increases risk of NAFLD, NASH, and liver
fibrosis.

**Population Frequency** -- The proportion of people in a defined population
who carry a particular allele or genotype.

**RAG** -- Retrieval-Augmented Generation. A pattern where an LLM generates
responses grounded in retrieved evidence rather than relying solely on its
training data.

**RDW** -- Red Cell Distribution Width. Measures variation in red blood cell
size. Elevated RDW is associated with inflammation, nutritional deficiency, and
accelerated aging. Has the largest positive coefficient in PhenoAge.

**Reference Range** -- The range of values considered normal for a biomarker in a
healthy population. Typically defined as the central 95% of a reference
population.

**rsID** -- Reference SNP cluster ID. A standardized identifier for a specific
single-nucleotide polymorphism (e.g., rs1801133 for MTHFR C677T).

**Similarity Score** -- A number (typically 0 to 1) indicating how closely a
search query matches a stored record in vector space. Higher scores mean closer
semantic match.

**SLCO1B1** -- Solute Carrier Organic Anion Transporter. The *5 variant reduces
hepatic statin uptake, increasing myopathy risk with simvastatin.

**SNOMED CT** -- Systematized Nomenclature of Medicine -- Clinical Terms. A
standardized clinical terminology used in electronic health records.

**SNP** -- Single-Nucleotide Polymorphism. A variation at a single position in
DNA. The most common type of genetic variation, with millions identified across
the human genome.

**Star Allele** -- A standardized nomenclature for pharmacogene variants,
written as an asterisk followed by a number (e.g., *1, *4, *17). Each star
allele represents a specific haplotype with defined functional consequences.

**TCF7L2** -- Transcription Factor 7-Like 2. The strongest genetic risk factor
for type 2 diabetes. The rs7903146 T allele increases risk by approximately
40% per copy.

**ThreadPoolExecutor** -- A Python concurrency mechanism used by the agent to
search all 14 collections in parallel, significantly reducing query latency.

**TPMT** -- Thiopurine S-methyltransferase. Metabolizes thiopurine drugs
(azathioprine, 6-MP). Poor metabolizers risk severe myelosuppression.

**TSH** -- Thyroid-Stimulating Hormone. The primary screening test for thyroid
function. Reference ranges vary by age, with higher upper limits acceptable in
elderly patients.

**Vector Database** -- A database optimized for storing and searching
high-dimensional vectors (embeddings). Enables semantic similarity search
rather than keyword matching.

**VDR** -- Vitamin D Receptor. Genetic variants affect vitamin D absorption
and metabolism, potentially requiring higher supplementation doses.

**VKORC1** -- Vitamin K epoxide reductase complex subunit 1. Warfarin's drug
target. Genetic variants determine warfarin dose sensitivity.

**WBC** -- White Blood Cell count. Measures immune system activity. Elevated WBC
contributes positively to PhenoAge (associated with accelerated aging).

---

*This guide is part of the HCLS AI Factory documentation. For questions,
contributions, or corrections, please open an issue in the project repository.*
