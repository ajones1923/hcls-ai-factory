# Precision Oncology Intelligence Agent: Foundations Learning Guide

> **An educational primer for clinicians, bioinformaticians, and software engineers
> who want to understand how this agent works -- no prior oncology or AI
> background required.**

---

## Welcome

This guide explains the ideas behind the Precision Oncology Intelligence Agent,
a software system that helps clinicians make treatment decisions for cancer
patients by combining genomic data, clinical evidence, and artificial
intelligence.

If you have ever wondered how a computer can read a patient's DNA, find the
right clinical trials, rank therapies, and package everything into a report
that a tumor board can act on -- this document walks through every step in
plain language with diagrams and examples.

---

## Who This Guide Is For

This guide serves three distinct audiences. Each chapter is written so that
all three can follow along, but look for the persona tags when a section
speaks directly to one group.

### Persona 1 -- The Oncologist

You treat cancer patients every day. You attend molecular tumor boards. You
want to know what this agent does, what evidence it relies on, and how to
interpret its output. You do not need to write code.

### Persona 2 -- The Bioinformatician

You work with genomic data -- VCF files, variant annotations, gene panels.
You want to understand how the agent ingests and searches that data, and how
to extend it with new collections or annotation sources.

### Persona 3 -- The Software Engineer

You build and deploy applications. You want to understand the architecture --
FastAPI, Streamlit, Milvus, Docker, Pydantic models -- and how to run,
modify, or contribute to the codebase.

---

## What You Will Learn

By the end of this guide you will be able to:

1. Explain what precision oncology is and why it matters.
2. Describe the role of a molecular tumor board (MTB) and how the agent
   supports it.
3. List the major genomic variant types (SNV, indel, CNV, fusion,
   rearrangement) and explain why they matter for treatment.
4. Distinguish five evidence levels (A through E) and map them to the
   AMP/ASCO/CAP classification framework.
5. Name the 11 Milvus vector collections the agent searches and explain
   the purpose of each one.
6. Explain how RAG (Retrieval-Augmented Generation) combines a vector
   database with a large language model.
7. Trace the agent workflow from a question through planning, searching,
   evaluating, and synthesizing an answer.
8. Describe how therapies are ranked with resistance awareness.
9. Explain how clinical trial matching combines deterministic filters and
   semantic search with composite scoring.
10. Identify key oncology biomarkers (TMB, MSI, PD-L1, HRD) and the
    therapies they predict.
11. Run a first query against the agent and read the resulting report.
12. Read an exported MTB packet in Markdown, JSON, PDF, or FHIR R4 format.

---

## Prerequisites

- **Oncologists:** A working web browser. That is all.
- **Bioinformaticians:** Familiarity with VCF files and command-line tools.
- **Software Engineers:** Python 3.10+, Docker, and a text editor.

No special hardware is needed to read this guide. To run the agent itself,
the reference platform is an NVIDIA DGX Spark ($3,999) with CUDA 12.x,
but any machine with Docker and a CPU will work for testing.

---

# Chapter 1: What Is Precision Oncology?

## The Old Approach

For decades, cancer treatment was based primarily on where the tumor was
located. Lung cancer patients got lung cancer drugs. Breast cancer patients
got breast cancer drugs. The same chemotherapy cocktail was given to everyone
with the same tumor site, even though response rates were often below 30%.

Think of it like prescribing the same pair of eyeglasses to every patient in
an optometrist's office -- some will see better, some will see worse, and
some will not benefit at all.

## The New Approach

Precision oncology flips this model. Instead of asking "where is the tumor?"
it asks "what is driving the tumor?" The answer comes from the tumor's DNA.

Here is the core idea:

```
  Traditional oncology         Precision oncology
  =====================        =====================
  Tumor site  --> Drug         Tumor DNA --> Driver mutation --> Matched drug
  (one-size-fits-all)          (tailored to each patient)
```

A patient with a BRAF V600E mutation in melanoma receives the same class of
drug (a BRAF inhibitor like dabrafenib) as a patient with BRAF V600E in
colorectal cancer -- because the molecular driver is the same, even though
the tumor sites are different.

## Why It Works

Cancer is fundamentally a disease of the genome. Tumors grow because specific
genes acquire mutations that tell cells to divide without stopping. If you
can identify the exact mutation driving growth, you can select a drug
designed to block that specific signal.

Key statistics that illustrate the shift:

- A single tumor can harbor tens of thousands of mutations, but typically
  only a handful are "drivers" that fuel growth.
- Over 40 genes now have FDA-approved targeted therapies.
- Some drugs (like larotrectinib for NTRK fusions) are approved regardless
  of tumor site -- these are called "tissue-agnostic" approvals.

## The Data Problem

The promise of precision oncology runs headlong into a data problem. To make
a good treatment decision, a clinician needs to synthesize information from
at least five different domains:

1. **Genomic data** -- the patient's tumor and germline variants.
2. **Clinical evidence** -- what published studies say about those variants.
3. **Treatment guidelines** -- what NCCN, ESMO, and ASCO recommend.
4. **Clinical trials** -- what experimental options are available.
5. **Resistance data** -- which drugs might fail because of co-mutations.

No single database holds all of this. The information is scattered across
CIViC, OncoKB, ClinicalTrials.gov, PubMed, and institutional guidelines.
Reading and cross-referencing these sources manually takes hours per patient.

That is the problem this agent solves.

---

# Chapter 2: The Molecular Tumor Board

## What Is an MTB?

A molecular tumor board (MTB) is a regular meeting where oncologists,
pathologists, geneticists, pharmacists, and other specialists review a
cancer patient's genomic test results and decide on a treatment plan.

Think of it as a panel of experts gathered around a conference table, each
bringing a different piece of the puzzle:

```
  +-----------+    +--------------+    +-------------+
  | Oncologist|    | Pathologist  |    | Geneticist  |
  |           |    |              |    |             |
  | "What drug|    | "The tumor is|    | "The VCF    |
  |  do we    |    |  grade III   |    |  shows EGFR |
  |  give?"   |    |  with high   |    |  L858R and  |
  |           |    |  mitotic     |    |  TP53 R248W"|
  |           |    |  index"      |    |             |
  +-----------+    +--------------+    +-------------+
        |                |                    |
        +----------- MTB Decision -----------+
                         |
                    Treatment Plan
```

## The Bottleneck

An MTB typically reviews 10-20 cases per session. For each case, someone
must prepare a packet that includes:

- A summary of the patient's variants and their clinical significance.
- Relevant published evidence for each actionable variant.
- A ranked list of therapy options.
- Matching clinical trials.
- Known resistance mechanisms for proposed treatments.

Preparing one packet manually can take 30-60 minutes. Multiply that by 15
cases and you have a full day of preparation for a single MTB session.

## How the Agent Helps

The Precision Oncology Intelligence Agent automates packet preparation. It
takes in a patient's genomic profile and produces a structured MTB packet
(the `MTBPacket` data model in the codebase) that includes:

- **Variant table** -- each variant with gene, type, evidence level, and
  associated drugs.
- **Evidence table** -- retrieved literature, guidelines, and outcomes.
- **Therapy ranking** -- prioritized treatments sorted by evidence level,
  with resistance flags.
- **Trial matches** -- clinical trials the patient may qualify for.
- **Open questions** -- uncertainties the board should discuss.
- **Citations** -- links to PubMed, ClinicalTrials.gov, and guidelines.

The packet can be exported to Markdown, JSON, PDF, or FHIR R4 format.

This does not replace the MTB. It gives the board better-prepared input so
clinicians can spend their time on judgment calls rather than data gathering.

### What an MTB Packet Contains

The `MTBPacket` model (`src/models.py`) fields include: `case_id`,
`patient_summary`, `patient_id`, `cancer_type`, `stage`, `variant_table`,
`evidence_table`, `therapy_ranking`, `trial_matches`, `open_questions`,
`generated_at`, and `citations`. The packet is generated automatically by
calling `agent.run()` and typically completes in under 5 seconds.

---

# Chapter 3: Genomic Variants and Actionability

## What Is a Genomic Variant?

Every cell in your body contains DNA -- a long string of chemical "letters"
(A, T, C, G) that encodes instructions for building proteins. A genomic
variant is any place where a patient's DNA differs from the standard
reference sequence.

Most variants are harmless. Some cause disease. A small number drive cancer.

## Types of Variants

The agent recognizes seven variant types, defined in the `VariantType` enum
in `src/models.py`:

```
  Variant Type        What It Means                     Example
  ==================  ================================  ====================
  SNV                 Single letter change              BRAF V600E (T>A)
  Indel               Small insertion or deletion       EGFR exon 19 del
  CNV amplification   Extra copies of a gene            HER2 amplification
  CNV deletion        Lost copies of a gene             CDKN2A homozygous del
  Fusion              Two genes joined together         EML4-ALK fusion
  Rearrangement       Part of a chromosome rearranged   RET rearrangement
  Structural variant  Large-scale chromosomal change    Complex SV
```

### Analogy: The Cookbook

Imagine DNA as a cookbook with 20,000 recipes (genes). Each recipe tells the
cell how to build a specific protein.

- An **SNV** is a single typo in a recipe -- "add 1 cup of sugar" becomes
  "add 1 cup of salt." Sometimes the dish still works; sometimes it is
  ruined.
- An **indel** is a few words inserted or removed -- "bake at 350 degrees
  for 30 minutes" becomes "bake at 350 degrees." Missing information.
- A **CNV amplification** is like photocopying one recipe 50 times. The cell
  makes far too much of that protein.
- A **fusion** is like taping two recipes together so the cell reads them as
  one, producing a hybrid protein that was never meant to exist.

## What Makes a Variant "Actionable"?

A variant is actionable when there is clinical evidence that it predicts
response (or resistance) to a specific therapy. Not every mutation is
actionable. In fact, most are not.

The agent's knowledge base (`src/knowledge.py`) tracks approximately 40
actionable gene targets. Here are ten of the most important:

| Gene   | Key Variant(s)       | Cancer Type(s)     | Matched Therapy        |
|--------|---------------------|--------------------|------------------------|
| EGFR   | L858R, ex19 del     | NSCLC              | Osimertinib            |
| BRAF   | V600E               | Melanoma, CRC      | Dabrafenib+trametinib  |
| ALK    | EML4-ALK fusion     | NSCLC              | Alectinib              |
| KRAS   | G12C                | NSCLC              | Sotorasib, adagrasib   |
| HER2   | Amplification       | Breast, gastric    | Trastuzumab deruxtecan |
| NTRK   | NTRK fusions        | Tissue-agnostic    | Larotrectinib          |
| RET    | Fusions, mutations  | NSCLC, thyroid     | Selpercatinib          |
| ROS1   | CD74-ROS1 fusion    | NSCLC              | Entrectinib            |
| BRCA1/2| Pathogenic variants | Ovarian, breast    | Olaparib (PARP inh.)   |
| PIK3CA | H1047R, E545K       | Breast             | Alpelisib              |

## Variant Classification: AMP/ASCO/CAP

The Association for Molecular Pathology (AMP), the American Society of
Clinical Oncology (ASCO), and the College of American Pathologists (CAP)
jointly published a four-tier system for classifying somatic variants:

```
  Tier I   -- Strong clinical significance
              Variants with FDA-approved therapies or included in
              professional guidelines (NCCN, ESMO).

  Tier II  -- Potential clinical significance
              Variants with clinical evidence from well-powered studies,
              clinical trials, or strong biological rationale.

  Tier III -- Unknown significance
              Variants not observed in existing databases or with
              conflicting evidence.

  Tier IV  -- Benign or likely benign
              Observed at high frequency in the general population;
              no disease association.
```

The agent's five-level evidence system (A through E) maps to a more granular
version of this framework. See Chapter 7 for the full mapping.

---

# Chapter 4: The Data Challenge

## Five Silos, One Patient

To make a treatment decision for a single patient, a clinician ideally
consults at least five major data sources:

```
  +-------------------+     +-------------------+     +-------------------+
  |  CIViC / OncoKB   |     |    PubMed / PMC   |     | ClinicalTrials.gov|
  |                   |     |                   |     |                   |
  | Curated variant-  |     | 35M+ biomedical   |     | 400,000+ trial   |
  | drug associations |     | research articles |     | registrations     |
  +-------------------+     +-------------------+     +-------------------+
           |                         |                         |
           +----------+--------------+--------------+----------+
                      |                             |
              +-------------------+     +-------------------+
              |   NCCN / ESMO /   |     |  Institutional    |
              |   ASCO Guidelines |     |  Tumor Boards     |
              |                   |     |                   |
              | Standard-of-care  |     | Case histories,   |
              | recommendations   |     | outcome records   |
              +-------------------+     +-------------------+
```

### The Problems

1. **Different formats.** CIViC uses a REST API. PubMed returns XML.
   ClinicalTrials.gov has its own JSON schema. Guidelines are often PDFs.
2. **Different vocabularies.** One source calls it "ERBB2"; another calls
   it "HER2." One says "non-small cell lung cancer"; another says "NSCLC."
3. **Constant change.** New papers are published daily. Trial statuses
   change. Guidelines are updated annually. A search done last month may
   already be outdated.
4. **Volume.** There are over 4.1 million records in ClinVar alone. Reading
   even a fraction manually is impossible.
5. **Cross-referencing.** The real insight comes from connecting variant
   data to therapy data to trial data to resistance data. No single source
   makes those connections.

### The Human Cost

A recent study found that preparing a single MTB case requires consulting an
average of 3-5 databases and takes 30-60 minutes of a specialist's time.
With 10-20 cases per session, that is a full day of preparation -- time
that could be spent with patients.

---

# Chapter 5: How RAG Solves the Data Challenge

## What Is RAG?

RAG stands for Retrieval-Augmented Generation. It is a technique that
combines two technologies:

1. **A search engine** that finds relevant documents from a knowledge base.
2. **A large language model (LLM)** that reads those documents and writes a
   coherent, cited answer.

### Analogy: The Research Librarian

Imagine you walk into a medical library and ask the librarian a question:

> "What targeted therapies are available for KRAS G12C-mutant NSCLC, and
> what resistance mechanisms should I watch for?"

The librarian does two things:

1. **Retrieves** relevant books, journals, and guidelines from the shelves
   (the "R" in RAG).
2. **Generates** a summary answer in plain language, citing the sources
   (the "G" in RAG).

RAG works the same way, except the "library" is a vector database (Milvus)
and the "librarian" is an LLM (Claude Sonnet 4.6).

## How Vector Search Works

Traditional search engines match keywords. If you search for "lung cancer
EGFR," they find documents containing those exact words. But what about a
document that says "non-small cell pulmonary adenocarcinoma with epidermal
growth factor receptor mutations"? That document is highly relevant but
shares no keywords with your query.

Vector search solves this by converting text into numbers. Here is how:

### Step 1: Embedding

An embedding model reads a piece of text and converts it into a list of
numbers (a "vector"). The agent uses the BAAI/bge-small-en-v1.5 model,
which produces a vector of 384 numbers for any input text.

```
  "EGFR L858R mutation in NSCLC"
          |
          |  Embedding model (BGE-small-en-v1.5)
          v
  [0.023, -0.112, 0.087, ..., 0.045]   <-- 384 numbers
```

Texts with similar meaning produce vectors that are close together in
384-dimensional space. Texts with different meanings produce vectors that
are far apart.

### Step 2: Indexing

When the agent ingests data, it embeds every document and stores the
resulting vector alongside the original text in Milvus. The agent uses
IVF_FLAT indexing with COSINE similarity -- a method that groups similar
vectors into clusters so searches run fast even over millions of records.

### Step 3: Searching

When a user asks a question, the agent:

1. Embeds the question into a 384-dimensional vector.
2. Searches Milvus for the vectors closest to the question vector.
3. Returns the original documents attached to those vectors.

```
  User question                    Vector database (Milvus)
  ==============                   ========================

  "EGFR resistance                 +------+------+------+
   mechanisms in NSCLC"   ----->   | Lit  | Var  | Res  | ...11 collections
                          embed    +------+------+------+
                          & search       |
                                         v
                                  Top-K most similar documents
                                         |
                                         v
                                  Claude Sonnet 4.6 (LLM)
                                         |
                                         v
                                  Cited, structured answer
```

### Step 4: Weighted Scoring

Not all collections are equally important. The agent assigns a weight to
each collection so that variant data and literature carry more influence
than older case records. The weights (from `config/settings.py`) are:

| Collection         | Weight | Purpose                          |
|--------------------|--------|----------------------------------|
| onco_variants      | 0.18   | Actionable variant evidence      |
| onco_literature    | 0.16   | Published research               |
| onco_therapies     | 0.14   | Drug mechanisms and indications  |
| onco_guidelines    | 0.12   | NCCN/ESMO/ASCO recommendations  |
| onco_trials        | 0.10   | Clinical trial summaries         |
| onco_biomarkers    | 0.08   | Predictive/prognostic biomarkers |
| onco_resistance    | 0.07   | Resistance mechanisms            |
| onco_pathways      | 0.06   | Signaling pathway context        |
| onco_outcomes      | 0.04   | Real-world treatment outcomes    |
| genomic_evidence   | 0.03   | VCF-derived genomic evidence     |
| onco_cases         | 0.02   | De-identified patient snapshots  |

These weights sum to 1.00 and can be tuned via environment variables with
the `ONCO_` prefix.

### Step 5: Synthesis

The top search results (up to 30 evidence items) are passed as context to
Claude Sonnet 4.6, which reads them and produces a structured, cited answer.

---

# Chapter 6: The 11 Knowledge Collections

The agent organizes its knowledge into 11 Milvus vector collections. Each
collection stores a different type of oncology information. All collections
use the same embedding model (BGE-small-en-v1.5, 384 dimensions) and the
same index type (IVF_FLAT with COSINE similarity).

Think of these collections as 11 specialized filing cabinets in the medical
library, each containing a different kind of document.

## Collection 1: onco_variants

**Purpose:** Clinically annotated genomic variants from CIViC and OncoKB.

**What it stores:** For each variant -- the gene name, variant name (e.g.,
"V600E"), variant type (SNV, fusion, etc.), cancer type, evidence level,
associated drugs, CIViC ID, clinical significance, and a text summary.

**Example record:**
```
  Gene:        EGFR
  Variant:     L858R
  Type:        SNV
  Cancer:      NSCLC
  Evidence:    A (FDA-approved)
  Drugs:       osimertinib, erlotinib, gefitinib, afatinib
  Summary:     "EGFR L858R is a sensitizing mutation found in ~40%
                of EGFR-mutant NSCLC cases..."
```

**Who cares:** Oncologists looking up whether a patient's variant has an
approved therapy. Bioinformaticians curating variant databases.

## Collection 2: onco_literature

**Purpose:** Chunks of published research from PubMed, PMC, and preprints.

**What it stores:** Title, text chunk, source type, year, cancer type,
gene, variant, keywords, and journal name.

**Why chunks?** A single paper may be 10 pages long. Embedding and
searching the entire paper at once would dilute the signal. Instead, the
agent splits each paper into smaller chunks (passages) so that a search
can pinpoint the most relevant paragraph.

**Who cares:** Everyone. Literature is the backbone of evidence-based
medicine.

## Collection 3: onco_therapies

**Purpose:** Approved and investigational therapies with mechanism of
action data.

**What it stores:** Drug name, therapy category (targeted, immunotherapy,
chemotherapy, ADC, bispecific, etc.), molecular targets, approved
indications, known resistance mechanisms, evidence level, and a text
summary.

**Example record:**
```
  Drug:        trastuzumab deruxtecan (T-DXd / Enhertu)
  Category:    ADC (antibody-drug conjugate)
  Targets:     HER2
  Indications: HER2+ breast, HER2+ gastric, HER2-mutant NSCLC
  Resistance:  HER2 loss, drug efflux, payload resistance
  Evidence:    A
```

**Who cares:** Oncologists selecting treatments. Pharmacists checking
mechanisms of action.

## Collection 4: onco_guidelines

**Purpose:** Recommendations from NCCN, ASCO, ESMO, and other
guideline-issuing organizations.

**What it stores:** Issuing organization, cancer type, guideline version,
year, key recommendations, text summary, and evidence level.

**Who cares:** Oncologists who need to confirm that a proposed treatment
is guideline-concordant.

## Collection 5: onco_trials

**Purpose:** Clinical trial summaries from ClinicalTrials.gov.

**What it stores:** NCT ID, title, text summary, phase, recruitment
status, sponsor, cancer types, biomarker eligibility criteria, enrollment
count, start year, and outcome summary.

**Example record:**
```
  ID:          NCT04136756
  Title:       "Phase III study of sotorasib vs. docetaxel in KRAS
                G12C-mutant NSCLC"
  Phase:       Phase 3
  Status:      Recruiting
  Biomarkers:  KRAS G12C
  Cancer:      NSCLC
```

**Who cares:** Oncologists matching patients to trials. Trial coordinators
evaluating eligibility.

## Collection 6: onco_biomarkers

**Purpose:** Predictive and prognostic biomarkers with testing context.

**What it stores:** Biomarker name, functional type (predictive,
prognostic, diagnostic, etc.), cancer types, predictive value, testing
method, clinical cutoff, text summary, and evidence level.

**Who cares:** Pathologists ordering tests. Oncologists interpreting
results.

## Collection 7: onco_resistance

**Purpose:** Documented mechanisms of therapeutic resistance.

**What it stores:** Primary therapy, gene, resistance mechanism, bypass
pathway, alternative therapies, and a text summary.

**Example record:**
```
  Therapy:      osimertinib
  Gene:         EGFR
  Mechanism:    C797S mutation in cis with T790M
  Bypass:       MET amplification
  Alternatives: amivantamab + lazertinib, platinum-based chemo
```

**Who cares:** Oncologists planning second-line therapy. Researchers
studying resistance biology.

## Collection 8: onco_pathways

**Purpose:** Oncogenic signaling pathways with therapeutic context.

**What it stores:** Pathway name, key genes, therapeutic targets,
cross-talk relationships, and a text summary.

The agent tracks 13 pathways: MAPK, PI3K/AKT/mTOR, DNA Damage Repair,
Cell Cycle, Apoptosis, WNT, NOTCH, Hedgehog, JAK/STAT, Angiogenesis,
Hippo, NF-kB, and TGF-beta.

**Why pathways matter:** Cancer cells often escape targeted therapy by
activating a parallel signaling pathway ("bypass resistance"). Knowing
which pathways cross-talk helps predict which escape routes a tumor
might take.

**Who cares:** Bioinformaticians building pathway models. Oncologists
understanding why a drug stopped working.

## Collection 9: onco_outcomes

**Purpose:** Real-world and trial-derived treatment outcome records.

**What it stores:** Case ID, therapy given, cancer type, response
category (RECIST criteria: CR, PR, SD, PD, NE), duration in months,
toxicities, baseline biomarkers, and a text summary.

**Who cares:** Oncologists comparing expected outcomes. Researchers
analyzing real-world evidence.

## Collection 10: onco_cases

**Purpose:** De-identified patient case snapshots for similarity search.

**What it stores:** Patient ID, cancer type, stage, variants, biomarkers,
prior therapies, and a text summary.

**How it is used:** When a new patient comes in, the agent can search
this collection to find similar historical cases -- patients with the
same cancer type, similar mutations, and similar biomarker profiles --
and see what treatments were tried and how they responded.

**Who cares:** Oncologists who want to learn from prior cases.

## Collection 11: genomic_evidence (Read-Only)

**Purpose:** VCF-derived genomic evidence from the upstream genomics
pipeline (Stage 1 of the HCLS AI Factory).

**What it stores:** Chromosome, position, reference/alternate alleles,
quality score, gene, consequence, impact, genotype, clinical
significance, rsID, disease associations, AlphaMissense pathogenicity
score and class, and a text summary.

**Why read-only:** This collection is populated by the genomics pipeline
(Parabricks + DeepVariant). The oncology agent reads it but never writes
to it. This separation ensures that the genomics pipeline is the single
source of truth for variant calls.

**Who cares:** Bioinformaticians who manage the variant calling pipeline.
Engineers maintaining data flow between stages.

## All 11 at a Glance

```
  +--------------------+    +--------------------+    +--------------------+
  | 1. onco_variants   |    | 2. onco_literature |    | 3. onco_therapies  |
  | Actionable somatic |    | PubMed / PMC /     |    | Drugs, MOA,        |
  | & germline variants|    | preprint chunks    |    | indications        |
  +--------------------+    +--------------------+    +--------------------+

  +--------------------+    +--------------------+    +--------------------+
  | 4. onco_guidelines |    | 5. onco_trials     |    | 6. onco_biomarkers |
  | NCCN / ESMO / ASCO |    | ClinicalTrials.gov |    | TMB, MSI, PD-L1,  |
  | recommendations    |    | summaries          |    | HRD, assay details |
  +--------------------+    +--------------------+    +--------------------+

  +--------------------+    +--------------------+    +--------------------+
  | 7. onco_resistance |    | 8. onco_pathways   |    | 9. onco_outcomes   |
  | Resistance mechs & |    | MAPK, PI3K, DDR,   |    | RECIST responses,  |
  | bypass strategies  |    | 13 pathways total  |    | duration, toxicity |
  +--------------------+    +--------------------+    +--------------------+

  +--------------------+    +--------------------+
  | 10. onco_cases     |    | 11. genomic_       |
  | De-identified      |    |     evidence       |
  | patient snapshots  |    | (read-only, VCF)   |
  +--------------------+    +--------------------+
```

---

# Chapter 7: Evidence Levels and Variant Classification

## The Five Evidence Levels

The agent uses a five-level evidence system, defined in the `EvidenceLevel`
enum in `src/models.py`. Each level indicates how strong the clinical
evidence is for a particular variant-drug association.

```
  Level   Name                   What It Means
  ======  =====================  ==========================================
  A       FDA-Approved           The variant-drug link has been validated
                                 by the FDA. There is a companion
                                 diagnostic. This is standard of care.

  B       Clinical Evidence      Well-powered clinical studies (phase II/III)
                                 show a significant association. May be
                                 included in guidelines but not yet FDA-
                                 approved for this specific indication.

  C       Case Reports           Evidence comes from case reports, small
                                 series, or early-phase trials. Promising
                                 but not yet confirmed in large studies.

  D       Preclinical            Evidence comes from cell lines, animal
                                 models, or in-vitro experiments. Not yet
                                 tested in patients.

  E       Computational          Evidence comes from computational
                                 prediction (e.g., AlphaMissense, protein
                                 structure modeling). No wet-lab or
                                 clinical validation yet.
```

### Analogy: The Confidence Thermometer

Think of evidence levels as a thermometer measuring confidence:

```
     HOT  |===| A  FDA says "this works"
          |===| B  Large studies confirm it
          |===| C  A few reports suggest it
          |===| D  Works in a lab dish
     COLD |===| E  A computer predicted it
```

Level A evidence means you can prescribe with high confidence. Level E
evidence means the association is speculative and should only be considered
when all other options are exhausted.

## Mapping to AMP/ASCO/CAP Tiers

The agent's evidence levels align with the AMP/ASCO/CAP classification:

```
  Agent Level   AMP/ASCO/CAP Tier   Clinical Action
  ============  ==================  ================================
  A             Tier I-A            FDA-approved therapy or guideline
  B             Tier I-B / II-C     Clinical trial evidence, consensus
  C             Tier II-D           Case reports, small studies
  D             Tier III            Unknown significance, preclinical
  E             Tier III            Computational prediction only
```

Variants classified as Tier IV (benign/likely benign) in the AMP/ASCO/CAP
system are filtered out before they reach the agent -- they are not
clinically actionable.

## VCF Files: Where Variant Data Comes From

Genomic variants reach the agent through VCF (Variant Call Format) files.
A VCF file is a standard text format produced by variant calling tools
like GATK, DeepVariant, Mutect2, and Parabricks.

### What a VCF Line Looks Like

```
  #CHROM  POS       ID           REF  ALT  QUAL   FILTER  INFO
  chr7    55181378  rs121913529  T    A    99     PASS    ANN=A|missense|HIGH|EGFR|...
```

Each line represents one variant. The columns tell you:

- **CHROM/POS** -- where in the genome the variant is (chromosome 7,
  position 55,181,378).
- **REF/ALT** -- what the reference base is (T) and what the patient has
  instead (A).
- **QUAL** -- confidence score for the variant call.
- **FILTER** -- whether the variant passed quality filters.
- **INFO** -- annotations added by tools like SnpEff, VEP, or GENEINFO.

### How the Agent Parses VCF Files

The VCF parser (`src/utils/vcf_parser.py`) extracts gene names and variant
consequences from the INFO field using three annotation formats:

1. **SnpEff (ANN= field)** -- extracts gene name and functional impact
   from the standard SnpEff annotation format.
2. **VEP (CSQ= field)** -- extracts gene symbol from the Ensembl Variant
   Effect Predictor annotation.
3. **GENEINFO** -- extracts gene name from the NCBI GENEINFO field.

The parser handles all common VCF 4.x formats produced by major variant
calling pipelines.

## How Evidence Level Affects Therapy Ranking

When the `TherapyRanker` scores candidate therapies (Chapter 8), evidence
level is the primary sorting criterion. The ranking order in the codebase
(`src/therapy_ranker.py`) is:

```python
  EVIDENCE_LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}
```

A therapy backed by Level A evidence will always rank above a therapy
backed by Level C evidence, regardless of other factors.

---

# Chapter 8: Therapy Ranking and Resistance

## How the Therapy Ranker Works

The `TherapyRanker` class (`src/therapy_ranker.py`) produces a prioritized
list of candidate therapies for a given patient. It uses a six-step process:

```
  Patient Profile
  (cancer type, variants, biomarkers, prior therapies)
        |
        v
  Step 1: Identify variant-driven therapies
        |   Look up each variant in ACTIONABLE_TARGETS
        v
  Step 2: Identify biomarker-driven therapies
        |   Check MSI-H, TMB-H, HRD, PD-L1 status
        v
  Step 3: Rank by evidence level (A > B > C > D > E)
        |
        v
  Step 4: Flag resistance
        |   Check if prior therapies predict resistance
        v
  Step 5: Flag contraindications
        |   Same drug class as a previously failed therapy
        v
  Step 6: Attach supporting evidence
        |   Retrieve from onco_therapies and onco_literature
        v
  Ranked Therapy List with Citations
```

### Example Walkthrough

Consider a patient with:
- Cancer: NSCLC, Stage IV
- Variants: EGFR L858R (SNV), TP53 R248W (SNV)
- Biomarkers: PD-L1 TPS 80%, TMB 8.2 mut/Mb
- Prior therapies: erlotinib (progressed after 11 months)

**Step 1:** EGFR L858R maps to osimertinib, erlotinib, gefitinib, afatinib.

**Step 2:** PD-L1 >= 50% maps to pembrolizumab. TMB 8.2 is below the
TMB-High threshold (typically >= 10 mut/Mb), so no TMB-driven therapy.

**Step 3:** All EGFR therapies are Level A. Pembrolizumab for PD-L1 is
also Level A.

**Step 4:** The patient previously received erlotinib and progressed.
This raises a resistance flag for erlotinib-class drugs. Known resistance
mechanisms include T790M mutation and MET amplification. The ranker
checks for these.

**Step 5:** Erlotinib is flagged as a contraindication (same drug, prior
failure). Gefitinib and afatinib (same generation TKIs) are downranked.

**Step 6:** Osimertinib rises to the top because it is a third-generation
TKI effective against T790M resistance.

**Final ranking:**
```
  Rank  Therapy            Evidence  Notes
  ====  =================  ========  ============================
  1     Osimertinib        A         Overcomes erlotinib resistance
  2     Pembrolizumab      A         PD-L1 >= 50%
  3     Amivantamab +      B         EGFR-directed bispecific +
        lazertinib                   3rd-gen TKI combination
  4     Erlotinib          A         CONTRAINDICATED (prior failure)
```

## Resistance Mechanisms

Resistance is one of the biggest challenges in precision oncology. A drug
that works brilliantly for 12 months can suddenly stop working because the
tumor evolves.

The agent tracks three categories of resistance:

### 1. On-Target Resistance

The drug's target gene acquires a new mutation that prevents the drug from
binding.

```
  Example: EGFR T790M
  - Patient on erlotinib (1st-gen EGFR TKI)
  - Tumor acquires T790M "gatekeeper" mutation
  - Erlotinib can no longer bind EGFR
  - Solution: Switch to osimertinib (3rd-gen, binds despite T790M)
```

### 2. Bypass Resistance

The tumor activates a parallel signaling pathway that circumvents the
blocked target.

```
  Example: MET amplification during EGFR therapy
  - Patient on osimertinib
  - Tumor amplifies MET gene, activating PI3K/AKT pathway
  - EGFR is still inhibited, but the cell survives via MET
  - Solution: Add MET inhibitor (capmatinib, tepotinib)
```

### 3. Lineage Plasticity

The tumor transforms into a different cell type entirely, changing from
(for example) an adenocarcinoma to a small cell carcinoma.

```
  Example: Small cell transformation in EGFR-mutant NSCLC
  - ~5-15% of EGFR-mutant NSCLC patients
  - Tumor loses EGFR dependency entirely
  - Targeted therapy no longer effective
  - Solution: Switch to small cell chemo (platinum + etoposide)
```

The `onco_resistance` collection stores all of these mechanisms so the
agent can proactively warn clinicians about potential resistance.

---

# Chapter 9: Clinical Trial Matching

## Why Trial Matching Matters

Clinical trials offer access to new therapies before they are FDA-approved.
For patients who have exhausted standard options, a trial may be their best
chance. But finding the right trial is hard:

- There are over 400,000 registered trials on ClinicalTrials.gov.
- Eligibility criteria are complex and written in free text.
- Trial statuses change frequently (recruiting, closed, suspended).
- Patients may qualify for a trial they never hear about.

## How the Trial Matcher Works

The `TrialMatcher` class (`src/trial_matcher.py`) uses a hybrid approach
that combines rule-based filtering with semantic search.

```
  Patient Profile
  (cancer type, variants, biomarkers, stage, age)
        |
        v
  Step 1: Deterministic filter
        |   Match cancer type (with fuzzy aliases)
        |   Filter to open statuses (Recruiting, Active, etc.)
        v
  Step 2: Semantic search
        |   Embed patient profile as a query vector
        |   Search onco_trials for similar trial descriptions
        v
  Step 3: Composite scoring
        |   Score = w1*biomarker_match + w2*semantic_sim
        |         + w3*phase_weight + w4*status_weight
        v
  Step 4: Age-based eligibility filtering
        |
        v
  Step 5: Structured explanation for each match
        |
        v
  Ranked Trial List
```

### Composite Scoring Weights

The trial match score is a weighted combination of four factors, configured
in `config/settings.py`:

```
  Factor            Weight   What It Measures
  ================  ======   =================================
  Biomarker match   0.40     Does the trial require this patient's
                             specific biomarkers?
  Semantic sim.     0.25     How similar is the trial description
                             to the patient's profile?
  Phase weight      0.20     Higher phases (III > II > I) get
                             higher weight.
  Status weight     0.15     Actively recruiting trials rank
                             above "not yet recruiting."
```

### Phase Weights

The matcher assigns weights to trial phases to prefer later-phase trials,
which are typically closer to producing actionable results:

```
  Phase         Weight
  ============  ======
  Phase 3       1.0
  Phase 2/3     0.9
  Phase 2       0.8
  Phase 1/2     0.7
  Phase 1       0.6
  Phase 4       0.5
```

### Status Weights

Trials that are actively enrolling patients are more useful than those
that are closed:

```
  Status                      Weight
  ==========================  ======
  Recruiting                  1.0
  Enrolling by invitation     0.8
  Active, not recruiting      0.6
  Not yet recruiting          0.4
```

Trials with statuses like Completed, Terminated, Suspended, or Withdrawn
are excluded from matching entirely (they are not in the OPEN_STATUSES
set).

---

# Chapter 10: Biomarkers in Oncology

## What Is a Biomarker?

A biomarker is a measurable characteristic that tells you something about
a patient's disease or their likely response to a treatment. In oncology,
biomarkers are used to select drugs, predict outcomes, and monitor
treatment response.

### Analogy: The Dashboard Gauges

Think of biomarkers as gauges on a car's dashboard. The speedometer tells
you how fast you are going. The fuel gauge tells you how far you can go.
The temperature gauge warns you of overheating. Similarly:

- **TMB** tells you how many mutations the tumor has (the "mutation load
  gauge").
- **MSI** tells you if the tumor's DNA repair machinery is broken (the
  "repair status gauge").
- **PD-L1** tells you if the tumor is hiding from the immune system (the
  "immune evasion gauge").
- **HRD** tells you if the tumor cannot repair DNA breaks properly (the
  "DNA repair gauge").

## The Eight Biomarker Types

The agent classifies biomarkers into eight functional categories, defined
in the `BiomarkerType` enum in `src/models.py`:

```
  Type                   What It Tells You
  =====================  ===============================================
  Predictive             Predicts response to a specific therapy.
                         Example: EGFR mutation predicts response to
                         osimertinib.

  Prognostic             Predicts disease outcome regardless of therapy.
                         Example: TP53 mutation is associated with worse
                         prognosis across many cancer types.

  Diagnostic             Helps confirm what type of cancer it is.
                         Example: PSA (prostate-specific antigen) for
                         prostate cancer screening.

  Monitoring             Tracks disease status over time.
                         Example: ctDNA (circulating tumor DNA) levels
                         to detect recurrence.

  Resistance             Indicates likely resistance to a therapy.
                         Example: EGFR T790M predicts resistance to
                         erlotinib.

  Pharmacodynamic        Measures whether the drug is hitting its target.
                         Example: phospho-ERK reduction during BRAF
                         inhibitor therapy.

  Screening              Used in population-level screening before
                         diagnosis. Example: BRCA1/2 testing in high-
                         risk individuals.

  Therapeutic Selection  Directly determines which therapy to use.
                         Example: HER2 status determines whether
                         trastuzumab is indicated.
```

## Key Biomarker-Therapy Mappings

The agent maintains a curated set of biomarker-to-therapy mappings. These
are among the most well-established in precision oncology:

```
  Biomarker         Cutoff / Status     Therapy                 Evidence
  ================  ==================  ======================  ========
  MSI-H             Microsatellite      Pembrolizumab           A
                    instability-high    (tissue-agnostic)

  TMB-H             >= 10 mut/Mb        Pembrolizumab           A
                                        (tissue-agnostic)

  HRD               Homologous          PARP inhibitors         A
                    recombination       (olaparib, niraparib,
                    deficiency          rucaparib, talazoparib)

  PD-L1 >= 50%      TPS >= 50%          Pembrolizumab           A
                                        (first-line NSCLC)

  NTRK fusion       Any NTRK1/2/3       Larotrectinib or        A
                    fusion              entrectinib
                                        (tissue-agnostic)

  EGFR mutation     L858R, ex19 del     Osimertinib             A

  BRAF V600E        V600E mutation      Dabrafenib +            A
                                        trametinib

  ALK fusion        EML4-ALK or other   Alectinib               A
                    ALK fusion

  KRAS G12C         G12C mutation       Sotorasib or            A
                                        adagrasib

  RET fusion/mut    Fusions or          Selpercatinib           A
                    activating muts
```

---

# Chapter 11: The Agent Workflow

## The Four-Phase Loop

The `OncoIntelligenceAgent` class (`src/agent.py`) implements a four-phase
workflow: Plan, Search, Evaluate, and Synthesize. If the evidence is
insufficient, the agent retries with broader queries (up to 2 retries).

```
                         User Question
                              |
                              v
                    +-------------------+
                    |   Phase 1: PLAN   |
                    |                   |
                    | - Parse question  |
                    | - Identify genes  |
                    | - Identify cancer |
                    |   types           |
                    | - Choose strategy |
                    | - Decompose into  |
                    |   sub-questions   |
                    +-------------------+
                              |
                              v
                    +-------------------+
                    |  Phase 2: SEARCH  |
                    |                   |
                    | - Embed query     |
                    | - Search 11       |
                    |   collections     |
                    | - Apply weighted  |
                    |   scoring         |
                    | - Collect hits    |
                    +-------------------+
                              |
                              v
                    +-------------------+
                    | Phase 3: EVALUATE |
                    |                   |
                    | - Count hits      |
                    | - Check diversity |
                    |   (collections)   |
                    | - Check quality   |
                    |   (avg score)     |
                    |                   |
                    |  sufficient? --+  |
                    |  partial?    --+  |
                    |  insufficient -+  |
                    +-------------------+
                        |           |
              sufficient|           |insufficient
                        |           |  (retry with
                        v           |   broader queries,
                    +-------------------+  up to 2x)
                    | Phase 4:SYNTHESIZE|<-----+
                    |                   |
                    | - Pass top-30     |
                    |   evidence items  |
                    |   to Claude LLM   |
                    | - Generate cited  |
                    |   answer          |
                    | - Format report   |
                    +-------------------+
                              |
                              v
                    AgentResponse
                    (answer + evidence + report)
```

## Phase 1: Plan

The agent analyzes the user's question to extract structured information:

- **Target genes:** Recognized from a set of 30+ known oncology genes
  (BRAF, EGFR, ALK, KRAS, HER2, NTRK, RET, MET, FGFR, PIK3CA, BRCA1,
  BRCA2, TP53, PTEN, and more).
- **Cancer types:** Recognized from 25+ canonical names and 60+ aliases
  (e.g., "lung cancer" maps to "NSCLC", "colon" maps to "COLORECTAL").
- **Topics:** Identified by keyword matching (resistance, biomarker,
  immunotherapy, clinical trial, combination therapy, etc.).
- **Search strategy:** One of three options:
  - `broad` -- no specific gene or cancer type identified; search
    everything.
  - `targeted` -- both a gene and cancer type are identified; focus the
    search.
  - `comparative` -- the question asks to compare two entities ("vs",
    "compare", "difference between").

Complex questions are decomposed into sub-questions. For example:

> "What are the resistance mechanisms to EGFR TKIs in NSCLC and what
> clinical trials are available?"

Becomes:
1. "Mechanisms of resistance to EGFR inhibitors"
2. "Active clinical trials targeting EGFR in NSCLC"

## Phase 2: Search

The RAG engine (`src/rag_engine.py`) embeds each query and searches all 11
collections in parallel using a thread pool (up to 6 workers). Results are
ranked by weighted cosine similarity.

The search is cross-collection: the same question is sent to all 11
collections, and results from every collection are merged into a single
ranked list.

## Phase 3: Evaluate

The agent assesses whether the retrieved evidence is sufficient:

- **Sufficient:** At least 3 hits from at least 2 different collections.
- **Partial:** Some hits but not enough diversity.
- **Insufficient:** No relevant hits at all.

If the evidence is insufficient or partial and retries remain, the agent
broadens the search. If the strategy was "targeted," it switches to "broad."
It also generates fallback queries (e.g., "{gene} oncology therapeutic
implications").

The minimum similarity score threshold is 0.30 -- hits below this cutoff
are discarded as irrelevant.

## Phase 4: Synthesize

The top evidence items (capped at 30) are passed to Claude Sonnet 4.6
along with a specialized oncology system prompt. The LLM reads the evidence
and generates a structured, cited answer.

The system prompt instructs the LLM to:
1. Cite evidence with PubMed or ClinicalTrials.gov links.
2. Connect genomic variants to therapy options, trials, and resistance.
3. Proactively highlight resistance mechanisms and safety concerns.
4. Reference NCCN, ESMO, or ASCO guidelines.
5. Acknowledge uncertainty and evidence gaps.

The final output is an `AgentResponse` containing the answer text, the
evidence used, the search plan, and a formatted Markdown report.

---

# Chapter 12: Running Your First Query

## Option 1: The Streamlit UI (Recommended for Oncologists)

The agent includes a Streamlit web interface (`app/oncology_ui.py`) that
provides a user-friendly chat-style interaction.

**Step 1:** Start the services (Docker required):

```bash
docker compose up -d
```

**Step 2:** Open your browser to `http://localhost:8526`.

**Step 3:** Type a question in the chat box. For example:

> "What targeted therapies are available for BRAF V600E melanoma, and what
> resistance mechanisms should I monitor?"

**Step 4:** The agent will:
1. Plan the search (identify BRAF, melanoma, resistance topic).
2. Search all 11 collections.
3. Evaluate the evidence.
4. Synthesize a cited answer.

**Step 5:** Review the report. It will include:
- A summary of BRAF V600E in melanoma.
- Ranked therapies (dabrafenib + trametinib at #1).
- Known resistance mechanisms (MEK1/2 mutations, NRAS activation, BRAF
  amplification).
- Relevant clinical trials.
- Citations to published literature.

## Option 2: The FastAPI Endpoint (For Engineers and Integrations)

The agent exposes a REST API on port 8527.

**Query endpoint:**

```bash
curl -X POST http://localhost:8527/api/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the significance of KRAS G12C in NSCLC?",
    "cancer_type": "nsclc",
    "gene": "KRAS"
  }'
```

The response is a JSON object matching the `AgentResponse` model:

```json
{
  "question": "What is the significance of KRAS G12C in NSCLC?",
  "answer": "KRAS G12C is found in approximately 13% of NSCLC cases...",
  "evidence": {
    "query": "...",
    "hits": [ ... ],
    "total_collections_searched": 11,
    "search_time_ms": 342.5
  },
  "knowledge_used": ["KRAS", "NSCLC", "sotorasib", "adagrasib"],
  "timestamp": "2026-03-12T14:22:01Z",
  "report": "# Oncology Intelligence Report\n..."
}
```

## Option 3: The Python API (For Bioinformaticians and Developers)

```python
from src.agent import OncoIntelligenceAgent
from src.rag_engine import OncoRAGEngine
from src.collections import OncoCollectionManager

# Connect to Milvus
manager = OncoCollectionManager(host="localhost", port=19530)
manager.connect()

# Initialize the RAG engine
rag = OncoRAGEngine(collection_manager=manager)

# Create the agent
agent = OncoIntelligenceAgent(rag_engine=rag)

# Run a query
response = agent.run("What ALK inhibitors are effective after crizotinib resistance?")

# Print the report
print(response.report)
```

## Interpreting the Output

Every agent response includes:

| Field              | What It Contains                                     |
|--------------------|------------------------------------------------------|
| `question`         | Your original question.                              |
| `answer`           | The synthesized, cited answer text.                  |
| `evidence`         | The `CrossCollectionResult` with all search hits.    |
| `knowledge_used`   | Genes, drugs, and concepts referenced in the answer. |
| `plan`             | The `SearchPlan` showing strategy and sub-questions. |
| `report`           | A formatted Markdown report ready for review.        |
| `timestamp`        | When the response was generated (UTC).               |

## Exporting Results

The `export` module (`src/export.py`) supports four output formats:

| Format   | Use Case                                  |
|----------|-------------------------------------------|
| Markdown | Quick review, pasting into documents      |
| JSON     | Integration with other systems            |
| PDF      | Formal reports for the medical record     |
| FHIR R4  | Clinical interoperability (HL7 standard)  |

### FHIR R4: Why It Matters

FHIR (Fast Healthcare Interoperability Resources) is the international
standard for exchanging healthcare data electronically. When the agent
exports a report in FHIR R4 format, it produces a structured document
that any FHIR-compliant electronic health record (EHR) system can import.

The export uses standard LOINC codes for genomic observations:

```
  Code       Meaning
  =========  =====================================
  81247-9    Master HL7 genomic report
  48018-6    Gene studied
  69548-6    Genetic variant assessment
```

This means an oncologist can generate a report in the agent and have it
automatically appear in the patient's medical record at their hospital --
no manual data entry required.

### PDF Reports

PDF exports include branded headers (NVIDIA green by default, configurable
via `ONCO_PDF_BRAND_COLOR_R/G/B` environment variables), formatted variant
tables, therapy rankings, and trial match summaries.

## Understanding the Technical Stack

For software engineers new to the project, here is how the major components
fit together:

```
  +-------------------+     +-------------------+     +-------------------+
  |   Streamlit UI    |     |   FastAPI REST     |     |   Python API      |
  |   (port 8526)     |     |   (port 8527)      |     |   (direct import) |
  |                   |     |                   |     |                   |
  | Browser-based     |     | JSON endpoints    |     | Programmatic      |
  | chat interface    |     | for integrations  |     | access from code  |
  +--------+----------+     +--------+----------+     +--------+----------+
           |                         |                         |
           +------------+------------+------------+------------+
                        |
               +--------v----------+
               | OncoIntelligence  |
               | Agent             |
               | (src/agent.py)    |
               +--------+----------+
                        |
               +--------v----------+
               | OncoRAGEngine     |
               | (src/rag_engine)  |
               +--------+----------+
                        |
           +------------+------------+
           |                         |
  +--------v----------+     +--------v----------+
  | BGE-small-en-v1.5 |     | Claude Sonnet 4.6 |
  | (embedding model) |     | (LLM synthesis)   |
  | 384-dim vectors   |     | Anthropic API     |
  +--------+----------+     +-------------------+
           |
  +--------v----------+
  | Milvus             |
  | (port 19530)       |
  | 11 collections     |
  | IVF_FLAT / COSINE  |
  +-------------------+
```

### Pydantic Models

Every data structure in the agent is a Pydantic model (`src/models.py`).
Pydantic provides automatic type validation, JSON serialization, and clear
error messages when data does not match the expected schema.

### Docker Containerization

The entire stack (Milvus, etcd, MinIO, the agent, the UI, monitoring)
runs in Docker containers via `docker-compose.dgx-spark.yml`. Start
everything with `docker compose up -d`; tear down with `docker compose down`.

---

# Appendix A: Glossary

Below are 60 terms used throughout this guide and the codebase, defined
in plain language.

| # | Term | Definition |
|---|------|------------|
| 1 | **Actionable variant** | A genomic mutation for which there is a clinically relevant therapy or clinical trial. |
| 2 | **ADC (antibody-drug conjugate)** | A drug that links an antibody (which finds the tumor cell) to a toxic payload (which kills it). Example: trastuzumab deruxtecan. |
| 3 | **Allele** | One of two copies of a gene at a given position in the genome. |
| 4 | **AlphaMissense** | A machine learning model from DeepMind that predicts whether a missense variant is pathogenic. |
| 5 | **AMP/ASCO/CAP** | The Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists -- jointly published the standard variant classification tiers. |
| 6 | **Biomarker** | A measurable characteristic (e.g., TMB, MSI, PD-L1) used to guide treatment decisions. |
| 7 | **Bispecific antibody** | An antibody engineered to bind two different targets simultaneously. Example: amivantamab (EGFR + MET). |
| 8 | **BRAF** | A gene that encodes a protein in the MAPK signaling pathway. BRAF V600E is the most common oncogenic BRAF mutation. |
| 9 | **Bypass resistance** | When a tumor activates an alternative signaling pathway to circumvent a blocked target. |
| 10 | **CIViC** | Clinical Interpretation of Variants in Cancer -- an open-source database of curated variant-drug associations. |
| 11 | **Claude** | An LLM (large language model) developed by Anthropic. This agent uses Claude Sonnet 4.6 for synthesis. |
| 12 | **Clinical trial** | A research study that tests a new treatment in human patients. Trials progress through phases I-IV. |
| 13 | **ClinVar** | An NCBI database of genomic variants and their clinical significance (over 4.1 million records). |
| 14 | **CNV (copy number variation)** | A change in the number of copies of a gene. Amplification means extra copies; deletion means lost copies. |
| 15 | **Companion diagnostic** | An FDA-approved test required before prescribing a specific drug (e.g., EGFR mutation test before osimertinib). |
| 16 | **Cosine similarity** | A mathematical measure of how similar two vectors are. Values range from -1 (opposite) to 1 (identical). |
| 17 | **ctDNA** | Circulating tumor DNA -- fragments of tumor DNA found in the bloodstream. Used for liquid biopsy. |
| 18 | **Docker** | A platform that packages software into containers so it runs consistently on any machine. |
| 19 | **Driver mutation** | A mutation that actively contributes to cancer growth (as opposed to a passenger mutation). |
| 20 | **Embedding** | A fixed-length numerical representation of text produced by a machine learning model. |
| 21 | **EGFR** | Epidermal Growth Factor Receptor -- a gene commonly mutated in NSCLC. |
| 22 | **ESMO** | European Society for Medical Oncology -- publishes treatment guidelines. |
| 23 | **Evidence level** | A rating (A through E) indicating the strength of clinical evidence for a variant-drug link. |
| 24 | **Exon** | A segment of a gene that encodes protein. Mutations in exons often affect protein function. |
| 25 | **FastAPI** | A modern Python web framework used by the agent's REST API (port 8527). |
| 26 | **FHIR R4** | Fast Healthcare Interoperability Resources, Release 4 -- an HL7 standard for exchanging healthcare data electronically. |
| 27 | **Fusion** | A genomic event where two separate genes become joined, producing a hybrid protein. Example: EML4-ALK. |
| 28 | **Gatekeeper mutation** | A mutation at a specific position in a kinase that blocks drug binding. Example: EGFR T790M. |
| 29 | **Germline variant** | A variant inherited from a parent, present in every cell of the body. Example: BRCA1 mutation. |
| 30 | **HRD (homologous recombination deficiency)** | A defect in DNA repair that makes tumors vulnerable to PARP inhibitors. |
| 31 | **Immunotherapy** | Treatment that helps the immune system recognize and attack cancer cells. Example: pembrolizumab (anti-PD-1). |
| 32 | **Indel** | A small insertion or deletion of nucleotides in the DNA sequence. |
| 33 | **IVF_FLAT** | A vector index type in Milvus that partitions vectors into clusters for fast approximate nearest-neighbor search. |
| 34 | **KRAS** | A gene in the RAS/MAPK pathway. KRAS G12C is targetable with sotorasib and adagrasib. |
| 35 | **LLM (large language model)** | An AI model trained on large amounts of text that can read, summarize, and generate natural language. |
| 36 | **Lineage plasticity** | When a tumor transforms from one cell type to another (e.g., adenocarcinoma to small cell carcinoma). |
| 37 | **MAPK pathway** | A signaling cascade (RAS-RAF-MEK-ERK) that controls cell growth. Frequently mutated in cancer. |
| 38 | **Milvus** | An open-source vector database used to store and search embeddings. |
| 39 | **Missense variant** | An SNV that changes one amino acid to another in the resulting protein. |
| 40 | **MSI (microsatellite instability)** | A condition where the DNA mismatch repair system is defective, leading to hypermutation. MSI-H tumors respond well to immunotherapy. |
| 41 | **MTB (molecular tumor board)** | A meeting where specialists review a patient's genomic data and decide on treatment. |
| 42 | **NCCN** | National Comprehensive Cancer Network -- publishes widely used treatment guidelines. |
| 43 | **NSCLC** | Non-small cell lung cancer -- the most common type of lung cancer (~85% of cases). |
| 44 | **OncoKB** | A precision oncology knowledge base maintained by Memorial Sloan Kettering Cancer Center. |
| 45 | **PARP inhibitor** | A drug that blocks the PARP enzyme, preventing DNA repair in HRD-positive tumors. Examples: olaparib, niraparib. |
| 46 | **Passenger mutation** | A mutation that does not contribute to cancer growth. Present by chance. |
| 47 | **PD-L1** | Programmed Death-Ligand 1 -- a protein that tumors use to hide from the immune system. High PD-L1 expression predicts response to anti-PD-1 immunotherapy. |
| 48 | **Pydantic** | A Python library for data validation using type annotations. The agent uses Pydantic models for all data structures. |
| 49 | **RAG (retrieval-augmented generation)** | A technique that combines document retrieval with LLM text generation to produce grounded, cited answers. |
| 50 | **RECIST** | Response Evaluation Criteria In Solid Tumors -- a standard for measuring tumor response: CR (complete response), PR (partial response), SD (stable disease), PD (progressive disease), NE (not evaluable). |
| 51 | **Resistance mechanism** | A biological change in the tumor that causes a previously effective drug to stop working. |
| 52 | **SNV (single nucleotide variant)** | A change of a single DNA letter. The simplest type of mutation. |
| 53 | **Somatic variant** | A mutation acquired during a person's lifetime, present only in tumor cells. |
| 54 | **Streamlit** | A Python framework for building interactive web applications. The agent's UI runs on port 8526. |
| 55 | **Targeted therapy** | A drug designed to block a specific molecular target (e.g., osimertinib blocks EGFR). |
| 56 | **TMB (tumor mutational burden)** | The number of mutations per megabase of DNA in a tumor. TMB-H (high) tumors may respond to immunotherapy. |
| 57 | **TKI (tyrosine kinase inhibitor)** | A class of targeted therapy that blocks tyrosine kinase enzymes. Examples: osimertinib, alectinib, sotorasib. |
| 58 | **Tissue-agnostic** | An FDA approval that applies regardless of tumor location, based on a molecular biomarker (e.g., NTRK fusion, MSI-H). |
| 59 | **VCF (Variant Call Format)** | A standard file format for storing genomic variant data. The agent parses VCF files annotated by SnpEff, VEP, or GENEINFO. |
| 60 | **Vector** | A list of numbers representing a piece of text in high-dimensional space. Used for similarity search. |

---

# Appendix B: Key Statistics

## Agent Codebase

| Metric                  | Value                                |
|-------------------------|--------------------------------------|
| Python source files     | 66                                   |
| Lines of code           | ~20,490                              |
| Test count              | 556                                  |
| Milvus collections      | 11                                   |
| Actionable gene targets | ~40                                  |
| Supported cancer types  | 26                                   |
| Evidence levels         | 5 (A through E)                      |
| Therapy categories      | 9                                    |
| Biomarker types         | 8                                    |
| Oncogenic pathways      | 13                                   |
| Response categories     | 5 (CR, PR, SD, PD, NE)              |
| Embedding model         | BAAI/bge-small-en-v1.5               |
| Embedding dimensions    | 384                                  |
| LLM                     | Claude Sonnet 4.6                    |
| Vector index type       | IVF_FLAT (COSINE)                    |
| Max evidence to LLM     | 30 items                             |
| Search parallelism      | 6 threads                            |
| Retry attempts          | Up to 2                              |
| Min. similarity score   | 0.30                                 |
| Min. sufficient hits    | 3 (from 2+ collections)              |

## Platform

| Component               | Value                                |
|-------------------------|--------------------------------------|
| Reference hardware      | NVIDIA DGX Spark ($3,999)            |
| GPU/CUDA                | CUDA 12.x                            |
| Embedding batch size    | 32                                   |
| API port                | 8527 (FastAPI)                       |
| UI port                 | 8526 (Streamlit)                     |
| Milvus port             | 19530                                |
| Export formats           | Markdown, JSON, PDF, FHIR R4         |

## Upstream Pipeline (HCLS AI Factory Context)

| Metric                  | Value                                |
|-------------------------|--------------------------------------|
| Demo genome variants    | 11.7 million                         |
| Searchable vectors      | 3.56 million                         |
| ClinVar records         | 4.1 million                          |
| AlphaMissense preds.    | 71 million                           |
| Gene coverage           | 201 genes                            |
| Therapeutic areas       | 13                                   |
| Druggable targets       | 171                                  |
| GPU variant calling     | 120-240 min (vs 24-48 hrs on CPU)    |
| Query response time     | < 5 seconds                          |

## Collection Weight Distribution

```
  onco_variants     |==================  0.18
  onco_literature   |================    0.16
  onco_therapies    |==============      0.14
  onco_guidelines   |============        0.12
  onco_trials       |==========          0.10
  onco_biomarkers   |========            0.08
  onco_resistance   |=======             0.07
  onco_pathways     |======              0.06
  onco_outcomes     |====                0.04
  genomic_evidence  |===                 0.03
  onco_cases        |==                  0.02
                    +----+----+----+----+
                    0   0.05  0.10 0.15 0.20
```

## Enum Quick Reference

### CancerType (26 values)
`nsclc`, `sclc`, `breast`, `colorectal`, `melanoma`, `pancreatic`,
`ovarian`, `prostate`, `renal`, `bladder`, `head_neck`, `hepatocellular`,
`gastric`, `glioblastoma`, `aml`, `cml`, `all`, `cll`, `dlbcl`,
`multiple_myeloma`, `cholangiocarcinoma`, `endometrial`, `thyroid`,
`mesothelioma`, `sarcoma`, `other`

### EvidenceLevel (5 values)
`A` (FDA-approved), `B` (clinical evidence), `C` (case reports),
`D` (preclinical), `E` (computational)

### TherapyCategory (9 values)
`targeted`, `immunotherapy`, `chemotherapy`, `hormonal`, `combination`,
`radiotherapy`, `cell_therapy`, `adc`, `bispecific`

### BiomarkerType (8 values)
`predictive`, `prognostic`, `diagnostic`, `monitoring`, `resistance`,
`pharmacodynamic`, `screening`, `therapeutic_selection`

### PathwayName (13 values)
`mapk`, `pi3k_akt_mtor`, `dna_damage_repair`, `cell_cycle`, `apoptosis`,
`wnt`, `notch`, `hedgehog`, `jak_stat`, `angiogenesis`, `hippo`,
`nf_kb`, `tgf_beta`

### ResponseCategory (5 values -- RECIST criteria)
`complete_response` (CR), `partial_response` (PR), `stable_disease` (SD),
`progressive_disease` (PD), `not_evaluable` (NE)

### VariantType (7 values)
`snv`, `indel`, `cnv_amplification`, `cnv_deletion`, `fusion`,
`rearrangement`, `structural_variant`

### TrialPhase (8 values)
`Early Phase 1`, `Phase 1`, `Phase 1/Phase 2`, `Phase 2`,
`Phase 2/Phase 3`, `Phase 3`, `Phase 4`, `N/A`

### GuidelineOrg (8 values)
`NCCN`, `ESMO`, `ASCO`, `WHO`, `CAP/AMP`, `FDA`, `EMA`, `AACR`

---

> **End of Learning Guide Foundations**
>
> For questions or contributions, see the project repository or contact
> the development team.
