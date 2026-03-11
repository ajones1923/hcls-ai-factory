# Learning Guide -- Foundations

**Precision Autoimmune Intelligence Agent | HCLS AI Factory**

Author: Adam Jones
Date: March 2026
License: Apache 2.0

---

## Welcome

You are reading the foundational learning guide for the Precision Autoimmune Intelligence Agent, an AI-powered clinical decision support system that unifies fragmented autoimmune disease data into a single, searchable intelligence platform. This system is part of the HCLS AI Factory, an end-to-end precision medicine platform that runs on a single NVIDIA DGX Spark.

Autoimmune diseases affect roughly 8% of the global population -- over 600 million people. Patients wait an average of 4.5 years and see 4 or more specialists before receiving a correct diagnosis. This agent exists to shorten that journey by putting all the relevant clinical, laboratory, genetic, and therapeutic evidence in one place, searchable with a single question.

### Who this guide is for

This guide is written for three audiences:

- **Biologists and clinicians** who understand autoimmune diseases and immunology but are new to AI, vector databases, and retrieval-augmented generation.
- **Data scientists and ML engineers** who understand embeddings and large language models but are new to rheumatology, autoantibody panels, HLA typing, and disease activity scoring.
- **Software developers** who want to understand how the system works end to end -- from Streamlit UI to Milvus vector search to Claude LLM synthesis.

You do not need to be an expert in all three areas. The whole point of this guide is to bring you up to speed on whichever parts are new to you.

### What you will learn

By the end of this guide, you will understand:

1. What autoimmune diseases are and why they are so hard to diagnose
2. Why data fragmentation across specialists is the central challenge
3. How Retrieval-Augmented Generation (RAG) works, from first principles
4. How this system searches 14 collections simultaneously using weighted parallel retrieval
5. How to use the UI to ask clinical questions and interpret the answers
6. What each of the 14 data collections contains and why it exists
7. How the knowledge base, five clinical engines, and diagnostic pipeline work together
8. How to set up and run the system locally
9. How to use the REST API
10. How the codebase is organized, file by file

### Prerequisites

- Basic Python knowledge (you can read a Python function and understand what it does)
- A computer with a terminal
- Curiosity about either autoimmune disease or AI-assisted clinical reasoning (or both)

No prior knowledge of rheumatology, vector databases, or large language models is required. We will build every concept from the ground up.

---

## Chapter 1: What Are Autoimmune Diseases?

### The immune system gone wrong

Your immune system is a defense force. It distinguishes "self" from "non-self" and attacks foreign invaders -- bacteria, viruses, fungi. When it works properly, you barely notice it. When it goes wrong, it attacks your own tissues.

An autoimmune disease occurs when the immune system loses tolerance for its own cells. It begins producing antibodies (called **autoantibodies**) that target the body's own proteins, or it activates immune cells that infiltrate and damage healthy organs. The result is chronic inflammation, tissue destruction, and a bewildering variety of symptoms that can affect almost any organ system.

Think of it like a security team that has been given the wrong photo on its watch list. Instead of watching for intruders, it starts detaining employees. The damage is real, the confusion is genuine, and the longer it continues unchecked, the worse the consequences.

### The 13 diseases covered by this agent

This system covers 13 autoimmune diseases spanning multiple organ systems:

| Disease | Primary Target | Key Autoantibodies |
|---------|---------------|-------------------|
| **Rheumatoid Arthritis (RA)** | Joints (synovium) | RF, anti-CCP |
| **Systemic Lupus Erythematosus (SLE)** | Multiple organs (kidneys, skin, joints, brain) | ANA, anti-dsDNA, anti-Smith |
| **Multiple Sclerosis (MS)** | Central nervous system (myelin) | Oligoclonal bands (CSF) |
| **Type 1 Diabetes (T1D)** | Pancreatic beta cells | anti-GAD, anti-IA-2, anti-insulin |
| **Inflammatory Bowel Disease (IBD)** | GI tract (Crohn's/UC) | ANCA (p-ANCA in UC) |
| **Psoriasis** | Skin (keratinocytes) | None specific (T-cell mediated) |
| **Ankylosing Spondylitis (AS)** | Spine, sacroiliac joints | None (HLA-B27 is key) |
| **Sjogren's Syndrome** | Exocrine glands (salivary, lacrimal) | anti-SSA/Ro, anti-SSB/La |
| **Systemic Sclerosis (SSc)** | Skin, lungs, vasculature | anti-Scl-70, anti-centromere |
| **Myasthenia Gravis (MG)** | Neuromuscular junction | AChR, anti-MuSK |
| **Celiac Disease** | Small intestine (villi) | anti-tTG, anti-endomysial |
| **Graves' Disease** | Thyroid (stimulatory) | TSI (thyroid-stimulating immunoglobulin) |
| **Hashimoto's Thyroiditis** | Thyroid (destructive) | anti-TPO, anti-thyroglobulin |

These 13 diseases are not 13 separate islands. They share genetic risk factors (especially HLA alleles), overlap in symptoms, co-occur in the same patients, and respond to some of the same biologic therapies. That interconnection is precisely why a unified intelligence platform is valuable.

### The diagnostic odyssey problem

Autoimmune diseases are notoriously difficult to diagnose:

- **Average time to diagnosis: 4.5 years.** Patients see an average of 4+ specialists before receiving a correct diagnosis. Many are initially told their symptoms are psychosomatic, anxiety-related, or simply unexplained.
- **Symptom overlap is enormous.** Fatigue, joint pain, skin rashes, and abnormal blood work can point to dozens of different conditions. A patient with fatigue, joint pain, and a positive ANA could have lupus, rheumatoid arthritis, Sjogren's syndrome, mixed connective tissue disease, or several other conditions.
- **Tests are imperfect.** ANA (antinuclear antibody) is positive in 95% of lupus patients, but it is also positive in 5-15% of healthy individuals. A positive ANA alone does not mean autoimmune disease.
- **Specialists see their own piece.** The rheumatologist sees the joints. The dermatologist sees the skin. The nephrologist sees the kidneys. The neurologist sees the nervous system. No single specialist sees the whole picture, and clinical records are often fragmented across separate electronic health record systems.

Consider a real pattern that this system was designed to catch: a patient named Sarah Mitchell (one of our nine demo patients) visited her PCP, a dermatologist, and a rheumatologist over 14 months before being diagnosed with lupus nephritis. Each specialist documented important findings, but those findings were scattered across different records. When the system ingests all 35 of Sarah's PDFs and analyzes them together, the pattern becomes obvious within seconds.

### Why AI can help

AI does not replace the rheumatologist. It does what a human physically cannot:

1. **Read everything at once.** The system ingests every clinical document -- lab results, imaging reports, progress notes, referral letters -- and makes them all searchable simultaneously.
2. **Cross-reference across specialists.** It connects a dermatologist's biopsy finding to a nephrologist's lab result to a rheumatologist's clinical assessment.
3. **Apply diagnostic criteria systematically.** The system knows the ACR/EULAR classification criteria for all 13 diseases and can evaluate them against the patient's data.
4. **Detect temporal patterns.** A rising anti-dsDNA titer combined with falling complement levels over 6 months suggests an impending lupus flare. The system tracks these patterns across longitudinal records.
5. **Work in seconds.** What takes a clinician hours of chart review takes the system seconds.

---

## Chapter 2: The Data Challenge

### Where autoimmune knowledge lives

A clinician evaluating a complex autoimmune patient needs information from at least 10 different categories:

1. **Clinical documents** -- Progress notes, consultation letters, discharge summaries from every specialist who has seen the patient
2. **Laboratory results** -- CBC, CMP, inflammatory markers (ESR, CRP), complement levels (C3, C4), urinalysis, and dozens of specialized tests
3. **Autoantibody panels** -- ANA with pattern, anti-dsDNA, anti-Smith, anti-SSA/SSB, RF, anti-CCP, and 18 more antibodies
4. **HLA typing** -- Genetic markers that confer disease susceptibility (HLA-B27 for AS, HLA-DRB1*04:01 for RA, etc.)
5. **Disease classification criteria** -- ACR/EULAR criteria that define when a set of findings qualifies as a specific diagnosis
6. **Disease activity scores** -- DAS28-CRP for RA, SLEDAI-2K for lupus, BASDAI for AS, and 17 more scoring systems
7. **Flare patterns** -- Biomarker trajectories that predict disease exacerbation before symptoms become severe
8. **Treatment databases** -- 22 biologic therapies with mechanism of action, approved indications, pharmacogenomic considerations, and monitoring requirements
9. **Published literature** -- Tens of thousands of papers spanning all 13 diseases
10. **Clinical trials** -- Active and completed trials for new treatments, biomarkers, and diagnostic approaches

### The problem: fragmented clinical records across specialists

Each of these data sources lives in a different system. The rheumatologist's EMR does not talk to the dermatologist's EMR. The geneticist's HLA report is a PDF in one system. The lab results are in another. The imaging is in a PACS system. The published literature is on PubMed. The clinical trials are on ClinicalTrials.gov.

A physician who wants to answer a cross-specialty question -- for example, "Given this patient's ANA pattern, HLA type, lab trajectory, and symptom timeline, what is the most likely diagnosis and what treatment should we consider?" -- would need to:

1. Pull records from every specialist visit
2. Extract and trend all lab values over time
3. Interpret the autoantibody panel in context
4. Check HLA associations
5. Apply classification criteria for multiple candidate diseases
6. Score disease activity
7. Review treatment guidelines and PGx considerations
8. Search for relevant clinical trials

This is slow, error-prone, and in practice, it rarely happens comprehensively. The result is delayed diagnosis, missed connections, and suboptimal treatment decisions.

### The solution: unified intelligence platform

The Precision Autoimmune Intelligence Agent solves this by:

1. **Ingesting** data from all of these sources into a single system
2. **Embedding** every piece of text as a 384-dimensional vector (a list of 384 numbers that captures the meaning of the text)
3. **Storing** these vectors in 14 purpose-built Milvus collections
4. **Searching** all 14 collections simultaneously when you ask a question, with collection-specific weights
5. **Augmenting** the search results with structured knowledge from a hand-curated knowledge base (22 HLA alleles, 24 autoantibodies, 20 disease activity scores, 13 flare patterns, 22 biologic therapies)
6. **Analyzing** the results through five specialized clinical engines (autoantibody interpretation, HLA association, disease activity scoring, flare prediction, biologic recommendation)
7. **Synthesizing** a grounded answer using Claude, with citations back to the original sources

The result: you ask one question, and the system searches all 14 collections in parallel, applies domain-specific knowledge, runs clinical analysis engines, and produces a comprehensive answer with cited evidence -- all in seconds.

---

## Chapter 3: What Is RAG?

### The limitation of language models

Large Language Models (LLMs) like Claude are trained on vast amounts of text. They can write coherently, reason about complex topics, and follow instructions. But they have a fundamental limitation: **they do not have access to your specific data.**

If you ask Claude about lupus, it can give you a general answer based on its training data. But it cannot tell you what is in your Milvus database. It cannot cite a specific patient's lab trend. It cannot look up whether Sarah Mitchell's anti-dsDNA titer has been rising over the past 6 months.

RAG solves this.

### RAG = Retrieval-Augmented Generation

RAG is a three-step pattern:

```
1. RETRIEVAL
   You ask a question. The system finds the most relevant documents
   from your database.

2. AUGMENTATION
   The retrieved documents are added to the prompt as context.
   The LLM now has your specific data in front of it.

3. GENERATION
   The LLM reads the evidence and generates an answer that is
   grounded in your data, with citations.
```

Think of it like this: imagine you are a medical resident preparing for a case presentation. Without RAG, you are answering from memory (which may be outdated or incomplete). With RAG, you have the patient's full chart, the latest lab results, the relevant classification criteria, and the treatment guidelines all spread out on the desk in front of you. You are answering from evidence.

### How retrieval works: embeddings and vector similarity

The retrieval step is the most technically interesting part, so let us break it down.

#### What is an embedding?

An embedding is a way of representing text as a list of numbers (a "vector") such that texts with similar meanings have similar numbers.

This system uses a model called **BAAI/bge-small-en-v1.5**, which converts any piece of text into a vector of **384 numbers** (384 dimensions). For example:

```
"Elevated anti-dsDNA with low C3/C4 suggests active lupus nephritis"
    --> [0.023, -0.156, 0.891, 0.044, ..., -0.312]  (384 numbers)

"Rising double-stranded DNA antibodies with complement consumption in SLE"
    --> [0.019, -0.148, 0.883, 0.051, ..., -0.298]  (384 numbers)
```

These two texts have very similar vectors because they express similar meanings, even though the words are different.

Conversely:

```
"HLA-B27 positive with sacroiliac joint erosion on MRI"
    --> [0.512, 0.078, -0.234, 0.667, ..., 0.112]  (384 numbers)
```

This vector looks very different from the lupus nephritis vectors because the topic is different (ankylosing spondylitis vs. lupus).

#### An analogy: GPS coordinates for meaning

Think of embeddings like GPS coordinates for meaning. Just as GPS coordinates place a physical location in a two-dimensional space (latitude and longitude), embeddings place a piece of text in a 384-dimensional meaning space. Texts about similar topics end up at nearby coordinates. Texts about unrelated topics end up far apart.

You cannot visualize 384 dimensions (nobody can), but the math works the same way as it does in two dimensions. To find texts similar to your question, you measure the "distance" between vectors.

#### How similarity search works

When you type a question into the Precision Autoimmune Intelligence Agent:

1. Your question is embedded into a 384-dimensional vector.
2. That vector is compared to every vector in the database using **cosine similarity** (a measure of how close two vectors are in direction, regardless of length).
3. The vectors with the highest similarity scores are returned as the most relevant results.

Cosine similarity ranges from 0 (completely unrelated) to 1 (identical meaning). In this system:

- Scores >= 0.80 indicate **high relevance**
- Scores >= 0.60 indicate **medium relevance**
- Scores < 0.60 indicate **low relevance**
- Scores below **0.40** are filtered out entirely (the SCORE_THRESHOLD)

The database that stores these vectors and performs fast similarity searches is **Milvus**, a purpose-built vector database. Milvus uses an indexing algorithm called **IVF_FLAT** (Inverted File with Flat quantization) with **nprobe=16** to search vectors in milliseconds rather than scanning them one by one.

### How the full RAG pipeline works in this system

Here is the complete pipeline, from question to answer:

```
User types: "What does Sarah Mitchell's rising anti-dsDNA and low C3 suggest?"
  |
  v
[1] EMBED THE QUESTION
    BGE-small-en-v1.5 converts the question to a 384-dim vector.
    The model prepends a special instruction prefix:
    "Represent this sentence for searching relevant passages: ..."
    This asymmetric prefix improves retrieval quality.
  |
  v
[2] SEARCH ALL 14 COLLECTIONS (in parallel, with weights)
    The query vector is sent to Milvus, which searches all 14
    collections simultaneously using ThreadPoolExecutor:
      - autoimmune_clinical_documents   (weight: 0.18)
      - autoimmune_patient_labs         (weight: 0.14)
      - autoimmune_autoantibody_panels  (weight: 0.12)
      - autoimmune_hla_associations     (weight: 0.08)
      - autoimmune_disease_criteria     (weight: 0.08)
      - autoimmune_disease_activity     (weight: 0.07)
      - autoimmune_flare_patterns       (weight: 0.06)
      - autoimmune_biologic_therapies   (weight: 0.06)
      - autoimmune_pgx_rules            (weight: 0.04)
      - autoimmune_clinical_trials      (weight: 0.05)
      - autoimmune_literature           (weight: 0.05)
      - autoimmune_patient_timelines    (weight: 0.03)
      - autoimmune_cross_disease        (weight: 0.02)
      - genomic_evidence                (weight: 0.02)

    Each collection returns its top-5 most similar results (TOP_K=5).
    Results are weighted by collection importance.
    Maximum evidence items returned: 30 (MAX_EVIDENCE).
  |
  v
[3] KNOWLEDGE BASE AUGMENTATION
    The system detects "anti-dsDNA" and "C3" in the question and adds
    structured knowledge:
      - anti-dsDNA: SLE-specific (specificity >95%), correlates with
        disease activity, especially nephritis
      - C3 consumption: active complement pathway, lupus nephritis marker
      - SLEDAI-2K scoring context
      - Flare prediction biomarker patterns for SLE
  |
  v
[4] MERGE, DEDUPLICATE, AND RANK
    All results are merged, duplicates removed, and ranked by
    weighted score. Citation relevance is assigned:
      - Score >= 0.80: HIGH relevance
      - Score >= 0.60: MEDIUM relevance
      - Score < 0.60: LOW relevance
  |
  v
[5] BUILD THE PROMPT
    The top-ranked evidence is formatted into a structured prompt:
      - Section per collection (Clinical Documents, Labs, Autoantibody
        Panels, etc.)
      - Each evidence item includes collection source and score
      - Knowledge base context is appended
      - The user's original question is stated
  |
  v
[6] LLM SYNTHESIS
    Claude receives the prompt with the system instruction
    (a detailed persona prompt covering autoimmune domain expertise)
    and generates a comprehensive, citation-rich answer.
    Response is streamed token by token to the UI.
  |
  v
[7] DISPLAY
    The Streamlit UI shows:
      - The streaming LLM response with citations
      - An expandable evidence panel with collection badges
      - Citation relevance indicators (HIGH/MEDIUM/LOW)
      - Download buttons (Markdown, FHIR R4, PDF)
```

---

## Chapter 4: System Overview

### Architecture at a high level

The Precision Autoimmune Intelligence Agent has four main layers:

```
+----------------------------------------------------+
|                  USER INTERFACE                      |
|       Streamlit (port 8531) + FastAPI (8532)        |
|  10 tabs: Chat, Patient Analysis, Differential Dx,  |
|  Labs, Timeline, Collections, Knowledge, Export,     |
|  Settings, Help                                      |
+-------------------------+--------------------------+
                          |
                          v
+----------------------------------------------------+
|                   RAG ENGINE                        |
|  AutoimmuneRAGEngine (src/rag_engine.py, 597 lines)|
|  - Weighted parallel multi-collection search        |
|  - Domain-specific reranking                        |
|  - Knowledge base augmentation                      |
|  - Citation relevance scoring                       |
|  - Conversation memory                              |
+-------------------------+--------------------------+
                          |
              +-----------+-----------+
              v                       v
+---------------------+  +--------------------------+
|   MILVUS VECTOR DB  |  |   KNOWLEDGE BASE v2.0.0  |
|   14 Collections    |  |   5 Domains:              |
|   IVF_FLAT / COSINE |  |   - 22 HLA alleles        |
|   384 dimensions    |  |   - 24 autoantibodies      |
|   nprobe = 16       |  |   - 20 activity scores     |
|   (BGE-small)       |  |   - 13 flare patterns      |
|                     |  |   - 22 biologic therapies  |
+---------------------+  +--------------------------+
              |
              v
+----------------------------------------------------+
|           FIVE CLINICAL ENGINES                     |
|  1. Autoantibody Interpretation                     |
|  2. HLA Association Analysis                        |
|  3. Disease Activity Scoring                        |
|  4. Flare Prediction                                |
|  5. Biologic Therapy Recommendation                 |
+-------------------------+--------------------------+
                          |
                          v
+----------------------------------------------------+
|                  CLAUDE LLM                         |
|  Claude (Anthropic API)                             |
|  System prompt: autoimmune domain expert persona    |
|  Streaming token generation                         |
+----------------------------------------------------+
```

### How the pieces connect

1. **The user** types a clinical question in the Streamlit chat interface (port 8531) or sends a POST request to the FastAPI REST API (port 8532).
2. **The RAG engine** embeds the question using BGE-small-en-v1.5, searches all 14 Milvus collections in parallel with collection-specific weights, and retrieves knowledge base context.
3. **Milvus** (port 19530) performs fast cosine-similarity search using IVF_FLAT indexes, returning the most relevant evidence from each collection.
4. **The knowledge base** adds structured facts -- HLA associations with odds ratios, autoantibody specificity data, disease activity scoring thresholds, flare biomarker patterns, and biologic therapy profiles with pharmacogenomic rules.
5. **The five clinical engines** run specialized analyses: interpreting autoantibody panels, evaluating HLA risk, scoring disease activity, predicting flares, and recommending biologic therapies.
6. **Claude** receives the evidence and knowledge context in a carefully constructed prompt, and generates a grounded, citation-rich response.

### The 14 collections as "specialized libraries"

Think of each collection as a specialized section in a medical library:

- **Clinical Documents** is the patient chart room -- all ingested clinical records (progress notes, consults, imaging reports)
- **Patient Labs** is the laboratory archive -- lab results with trending and flag analysis
- **Autoantibody Panels** is the immunology bench -- autoantibody test results and pattern interpretation
- **HLA Associations** is the genetics reading room -- HLA allele to disease risk mappings
- **Disease Criteria** is the classification standards shelf -- ACR/EULAR diagnostic criteria
- **Disease Activity** is the scoring reference desk -- activity score definitions and thresholds
- **Flare Patterns** is the early warning system -- biomarker trajectories that predict flares
- **Biologic Therapies** is the pharmacy formulary -- drug profiles with mechanisms, indications, and PGx
- **PGx Rules** is the pharmacogenomics lab -- genotype-to-dosing rules
- **Clinical Trials** is the research registry -- active and completed autoimmune trials
- **Literature** is the research library -- published papers and reviews
- **Patient Timelines** is the longitudinal record room -- diagnostic journey timelines
- **Cross-Disease** is the overlap syndromes desk -- connections between co-occurring conditions
- **Genomic Evidence** is the genome center -- shared read-only variant data from the genomics pipeline

When you ask a question, the system does not search just one library section -- it searches all fourteen simultaneously and cross-references the findings.

### The five-stage clinical analysis pipeline

Beyond simple question answering, the system runs a structured clinical analysis pipeline:

```
Stage 1: Autoantibody Interpretation
  Input:  Autoantibody panel results (ANA, anti-dsDNA, RF, etc.)
  Output: Disease associations with specificity/sensitivity data

Stage 2: HLA Association Analysis
  Input:  HLA typing results (e.g., HLA-B*27:05, HLA-DRB1*04:01)
  Output: Disease risk with odds ratios and population frequencies

Stage 3: Disease Activity Scoring
  Input:  Lab values and clinical parameters
  Output: Calculated activity scores (DAS28, SLEDAI, BASDAI, etc.)
          with severity classification (remission/low/moderate/high)

Stage 4: Flare Prediction
  Input:  Longitudinal biomarker trends
  Output: Flare risk assessment (low/moderate/high) with
          biomarker pattern analysis

Stage 5: Biologic Therapy Recommendation
  Input:  Diagnosis, activity level, prior treatments, HLA/PGx data
  Output: Ranked therapy options with mechanism, evidence, PGx
          considerations, and monitoring requirements
```

---

## Chapter 5: Your First Query

This chapter walks you through the experience of using the system for the first time.

### Opening the UI

Once the system is running (see Chapter 9 for setup), open your browser and navigate to:

```
http://localhost:8531
```

You will see the Precision Autoimmune Intelligence Agent interface with a clean layout and a sidebar with configuration options.

### The 10 tabs explained

The UI is organized into 10 tabs, each serving a distinct clinical workflow:

| Tab | Purpose |
|-----|---------|
| **Chat** | Ask free-text clinical questions with RAG-powered answers |
| **Patient Analysis** | Run the full five-stage clinical pipeline on a patient |
| **Differential Dx** | Generate differential diagnosis from symptoms and labs |
| **Labs** | Explore and trend laboratory results over time |
| **Timeline** | Visualize a patient's diagnostic journey chronologically |
| **Collections** | View and manage the 14 Milvus collections |
| **Knowledge** | Browse the structured knowledge base (HLA, antibodies, drugs) |
| **Export** | Export analysis results in Markdown, FHIR R4, or PDF format |
| **Settings** | Configure search parameters, API keys, and display options |
| **Help** | Documentation, example queries, and system information |

### Asking a clinical question

Navigate to the **Chat** tab and type a question. For your first query, try:

```
What does a positive ANA with speckled pattern, elevated anti-dsDNA, and
low C3/C4 suggest in a 34-year-old female with joint pain and malar rash?
```

Press Enter. Here is what happens:

1. **Search status** appears, showing "Searching across autoimmune data sources..."
2. The system reports how many results it found and from which collections.
3. An **Evidence Sources** expander appears below the status, showing the raw evidence cards with collection badges.
4. The **LLM response** streams in token by token, with markdown formatting and citations.

### Understanding the response and evidence

The response will contain several types of content:

**Cited evidence:** References to specific clinical documents, lab results, knowledge base entries, and literature findings. Each citation includes the collection it came from and the relevance score.

**Cross-specialty insights:** The response connects evidence from different domains. For example, it might explain how autoantibody results (from the panels collection) relate to disease classification criteria (from the criteria collection), HLA risk factors (from the HLA collection), and treatment options (from the therapies collection).

**Structured clinical reasoning:** For complex cases, the response includes categorized sections (e.g., "Autoantibody Interpretation," "Differential Diagnosis," "Disease Activity Assessment," "Treatment Considerations").

### Reading citation scores

| Score Range | Relevance | What it means |
|-------------|-----------|--------------|
| 0.80 - 1.00 | **HIGH** | Strong semantic match. This evidence is directly relevant to your question. |
| 0.60 - 0.79 | **MEDIUM** | Partial match. The evidence is related but may address a subtopic or adjacent concept. |
| 0.40 - 0.59 | **LOW** | Weak match. The evidence has thematic overlap but may not directly answer your question. |
| Below 0.40 | Filtered out | Not returned. Below the SCORE_THRESHOLD. |

### Using sidebar controls

The sidebar gives you control over how searches are performed:

- **Disease Filter**: Restrict results to a specific disease (e.g., SLE, RA, AS)
- **Collection Toggles**: Enable or disable specific collections, each showing live record counts
- **Search Depth**: Adjust TOP_K (results per collection) and MAX_EVIDENCE (total results)
- **Score Threshold**: Set the minimum similarity score for included results
- **Patient Selector**: Choose a demo patient for analysis

### Downloading results

After each response, download buttons appear offering three formats:

- **Markdown**: A formatted `.md` file with the query, response, evidence, and metadata
- **FHIR R4**: A standards-compliant FHIR R4 Bundle for healthcare system integration
- **PDF**: A formatted PDF report suitable for clinical documentation

---

## Chapter 6: Understanding the 14 Collections

### Overview table

| # | Collection | Weight | Description |
|---|-----------|--------|-------------|
| 1 | `autoimmune_clinical_documents` | 0.18 | Ingested patient records (PDFs) |
| 2 | `autoimmune_patient_labs` | 0.14 | Lab results with flag analysis |
| 3 | `autoimmune_autoantibody_panels` | 0.12 | Autoantibody test result panels |
| 4 | `autoimmune_hla_associations` | 0.08 | HLA allele to disease risk mapping |
| 5 | `autoimmune_disease_criteria` | 0.08 | ACR/EULAR classification criteria |
| 6 | `autoimmune_disease_activity` | 0.07 | Activity scoring reference |
| 7 | `autoimmune_flare_patterns` | 0.06 | Flare prediction biomarker patterns |
| 8 | `autoimmune_biologic_therapies` | 0.06 | Biologic drug database with PGx |
| 9 | `autoimmune_pgx_rules` | 0.04 | Pharmacogenomic dosing rules |
| 10 | `autoimmune_clinical_trials` | 0.05 | Autoimmune clinical trials |
| 11 | `autoimmune_literature` | 0.05 | Published literature |
| 12 | `autoimmune_patient_timelines` | 0.03 | Patient diagnostic timelines |
| 13 | `autoimmune_cross_disease` | 0.02 | Cross-disease / overlap syndromes |
| 14 | `genomic_evidence` | 0.02 | Shared read-only genomic data |

The weights sum to 1.00. They determine how much influence each collection has on the final ranked results. Clinical documents and labs have the highest weights because they are most directly relevant to patient-specific questions.

### Collection details

#### 1. autoimmune_clinical_documents (weight: 0.18)

**What it contains:** Ingested patient clinical records -- progress notes, consultation letters, imaging reports, pathology reports, discharge summaries, referral letters, and medication reconciliations. Each record is chunked into semantically coherent segments (2,500 characters with 200-character overlap) and embedded.

**Why it matters:** This is the patient's actual chart. When you ask about a specific patient, this collection contains the primary clinical documentation. It has the highest weight because patient-specific questions almost always require chart data.

**Example questions it helps answer:**
- "What did Dr. Kim document in Sarah's rheumatology consult?"
- "When was the malar rash first noted in the clinical record?"
- "What were the findings on the renal biopsy?"

**Demo data:** Sarah Mitchell (35 PDFs), Maya Rodriguez (28 PDFs), David Park (26 PDFs), Linda Chen (20 PDFs), Rachel Thompson (22 PDFs), plus Emma Williams, James Cooper, Karen Foster, and Michael Torres.

#### 2. autoimmune_patient_labs (weight: 0.14)

**What it contains:** Structured laboratory results including CBC, CMP, inflammatory markers (ESR, CRP), complement levels (C3, C4), immunoglobulin levels, urinalysis, and specialized tests. Each result includes the value, reference range, flag (normal/abnormal/critical), and date.

**Why it matters:** Lab trends tell the story of disease progression. A single lab value is a snapshot; a series of lab values over time reveals whether a patient is improving, stable, or heading toward a flare. This collection enables temporal queries like "How has the anti-dsDNA titer changed over the past year?"

**Example questions it helps answer:**
- "Show the trend of C3 and C4 levels for this patient"
- "Which labs were flagged as critical in the last 6 months?"
- "Is the CRP trending up or down?"

#### 3. autoimmune_autoantibody_panels (weight: 0.12)

**What it contains:** Autoantibody test results organized as panels. Each panel records the antibodies tested, results (positive/negative/titer), ANA pattern (homogeneous, speckled, nucleolar, centromere), and clinical interpretation context.

**Why it matters:** Autoantibodies are the cornerstone of autoimmune diagnosis. A positive anti-CCP with high RF strongly suggests RA. A positive anti-dsDNA with homogeneous ANA pattern points to SLE. The pattern matters as much as the positivity, and this collection captures both.

**The 24 autoantibodies tracked by this system:**
ANA (with patterns: homogeneous, speckled, nucleolar, centromere), anti-dsDNA, anti-Smith, anti-Scl-70, anti-centromere, anti-SSA/SSB, anti-Ro/La, anti-Jo-1, RF, anti-CCP, AChR, anti-MuSK, anti-tTG, anti-endomysial, TSI, anti-TPO, anti-thyroglobulin, ANCA (c-ANCA/PR3, p-ANCA/MPO).

**Example questions it helps answer:**
- "What does a speckled ANA pattern with positive anti-SSA mean?"
- "Which autoantibodies were positive in Sarah's workup?"
- "What is the specificity of anti-CCP for rheumatoid arthritis?"

#### 4. autoimmune_hla_associations (weight: 0.08)

**What it contains:** HLA allele to disease risk mappings, including the allele name, associated disease, odds ratio, population frequency, and supporting literature (PMID references).

**Why it matters:** HLA alleles are the strongest genetic risk factors for most autoimmune diseases. HLA-B*27:05 confers an odds ratio of 87.4 for ankylosing spondylitis -- one of the strongest known genetic associations in all of medicine. This collection lets the system evaluate genetic risk when HLA typing is available.

**Key HLA associations in this system:**

| Allele | Disease | Odds Ratio |
|--------|---------|-----------|
| HLA-B*27:05 | Ankylosing Spondylitis | 87.4 |
| HLA-C*06:02 | Psoriasis | 10.0 |
| HLA-DQB1*02:01 | Celiac Disease | 7.0 |
| HLA-DRB1*04:01 | Rheumatoid Arthritis | 4.2 |
| HLA-DRB1*03:01 | Type 1 Diabetes | 3.6 |
| HLA-DRB1*15:01 | Multiple Sclerosis | 3.1 |

**Example questions it helps answer:**
- "What does HLA-B27 positivity mean for this patient's back pain?"
- "Which HLA alleles are associated with Sjogren's syndrome?"
- "How strong is the HLA-DRB1*04:01 association with RA?"

#### 5. autoimmune_disease_criteria (weight: 0.08)

**What it contains:** ACR/EULAR classification criteria for the 13 autoimmune diseases. Each record describes the criteria set, the required number of points or features, and the sensitivity/specificity of the criteria.

**Why it matters:** Classification criteria are the formal rules that define when a set of clinical findings qualifies as a specific diagnosis. The 2019 ACR/EULAR SLE criteria, for example, require a positive ANA plus a weighted score of 10 or more from clinical and immunological domains. The system uses these criteria to evaluate whether a patient's data meets the threshold for a given diagnosis.

**Example questions it helps answer:**
- "Does this patient meet the ACR/EULAR criteria for SLE?"
- "What are the 2010 ACR/EULAR criteria for RA?"
- "How are the modified New York criteria for AS applied?"

#### 6. autoimmune_disease_activity (weight: 0.07)

**What it contains:** Disease activity scoring systems with thresholds for remission, low, moderate, and high activity. Each record includes the score name, disease, component variables, calculation method, and threshold values.

**Why it matters:** Activity scores quantify how active a disease is at a given point in time. They guide treatment decisions: a patient in remission may continue current therapy, while a patient with high activity may need escalation to biologics.

**The 20 activity scores in this system:**

| Score | Disease | Thresholds (Remission / Low / Moderate / High) |
|-------|---------|------------------------------------------------|
| DAS28-CRP | RA | <2.6 / 2.6-3.2 / 3.2-5.1 / >5.1 |
| DAS28-ESR | RA | <2.6 / 2.6-3.2 / 3.2-5.1 / >5.1 |
| CDAI | RA | <=2.8 / 2.8-10 / 10-22 / >22 |
| SDAI | RA | <=3.3 / 3.3-11 / 11-26 / >26 |
| SLEDAI-2K | SLE | 0 / 1-5 / 6-10 / >10 |
| BASDAI | AS | <1 / 1-3 / 3-4 / >4 |
| ASDAS | AS | <1.3 / 1.3-2.1 / 2.1-3.5 / >3.5 |
| Mayo Score | UC | 0-1 / 2-4 / 5-7 / 8-12 |
| Harvey-Bradshaw | Crohn's | <5 / 5-7 / 8-16 / >16 |
| PASI | Psoriasis | 0 / <5 / 5-10 / >10 |
| DAPSA | PsA | <=4 / 4-14 / 14-28 / >28 |
| ESSDAI | Sjogren's | 0 / 1-4 / 5-13 / >=14 |
| mRSS | SSc | -- (0-51 scale, higher = worse) |
| EDSS | MS | 0 / 1-3.5 / 4-6.5 / >=7 |
| QMGS | MG | 0-39 scale |
| MG-ADL | MG | 0-24 scale |
| Marsh Score | Celiac | 0-3c scale |
| Burch-Wartofsky | Graves | <25 / 25-44 / >=45 |

**Example questions it helps answer:**
- "What is the DAS28-CRP score for a patient with these parameters?"
- "Is a SLEDAI-2K of 8 considered moderate or high activity?"
- "What BASDAI score threshold triggers biologic therapy consideration?"

#### 7. autoimmune_flare_patterns (weight: 0.06)

**What it contains:** Biomarker trajectory patterns that predict disease flares for each of the 13 diseases. Each pattern describes the key biomarkers to monitor, the direction of change that signals impending flare, the typical lead time (how far in advance the biomarker change appears), and the recommended monitoring frequency.

**Why it matters:** Catching a flare early can prevent irreversible organ damage. If a lupus patient's anti-dsDNA is rising while their C3/C4 is falling, a flare may be 4-8 weeks away -- even if the patient feels fine right now. The system tracks these patterns longitudinally.

**Example questions it helps answer:**
- "What biomarker patterns predict an SLE flare?"
- "How far in advance can we detect an RA flare from lab trends?"
- "What should I monitor in a Crohn's patient to predict relapse?"

#### 8. autoimmune_biologic_therapies (weight: 0.06)

**What it contains:** Profiles for 22 biologic therapies, including the drug name, mechanism of action (TNF inhibitor, IL-6 blocker, B-cell depleter, etc.), approved indications, pharmacogenomic considerations, key clinical trials, monitoring requirements, and common adverse effects.

**Why it matters:** Biologic therapy selection is one of the most consequential clinical decisions in autoimmune disease. The right drug can induce remission; the wrong drug wastes months of time and risks adverse effects. PGx considerations (e.g., TPMT status for azathioprine, HLA-B*58:01 for allopurinol) can prevent serious toxicity.

**Example questions it helps answer:**
- "What biologics are approved for both RA and psoriasis?"
- "What PGx testing should be done before starting a TNF inhibitor?"
- "Compare rituximab and belimumab for lupus nephritis"

#### 9. autoimmune_pgx_rules (weight: 0.04)

**What it contains:** Pharmacogenomic dosing rules linking specific genetic variants to drug dose adjustments, contraindications, or monitoring requirements. Rules follow CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines.

**Why it matters:** Pharmacogenomics can prevent life-threatening adverse drug reactions. A patient who is a TPMT poor metabolizer should not receive standard-dose azathioprine. A patient with HLA-B*58:01 should avoid allopurinol due to the risk of severe cutaneous adverse reactions.

**Example questions it helps answer:**
- "What TPMT testing is needed before azathioprine?"
- "How should methotrexate dosing be adjusted for MTHFR variants?"
- "Which HLA alleles affect drug safety in autoimmune therapy?"

#### 10. autoimmune_clinical_trials (weight: 0.05)

**What it contains:** Active and completed clinical trials for autoimmune therapies, biomarkers, and diagnostic approaches. Records include trial ID, title, phase, status, disease, intervention, sponsor, and enrollment.

**Why it matters:** Clinical trials represent the cutting edge of autoimmune treatment. When a patient has failed standard therapy, finding a relevant trial may be the next step. This collection enables trial discovery by disease, mechanism, or clinical question.

**Example questions it helps answer:**
- "Are there active Phase 3 trials for lupus nephritis?"
- "What JAK inhibitor trials are recruiting for RA?"
- "Which trials are studying CAR-T for autoimmune diseases?"

#### 11. autoimmune_literature (weight: 0.05)

**What it contains:** Published literature including research papers, review articles, case reports, and meta-analyses related to the 13 autoimmune diseases.

**Why it matters:** Published evidence is the foundation of clinical decision-making. When the system cites a treatment recommendation or a diagnostic threshold, the literature collection provides the supporting evidence.

**Example questions it helps answer:**
- "What is the evidence for hydroxychloroquine in lupus?"
- "What do recent meta-analyses say about TNF inhibitor safety in pregnancy?"
- "What are the latest diagnostic criteria updates for Sjogren's?"

#### 12. autoimmune_patient_timelines (weight: 0.03)

**What it contains:** Chronological diagnostic journey records for each patient, capturing when symptoms first appeared, when each specialist was seen, when key tests were ordered, and when diagnoses were made or revised.

**Why it matters:** The diagnostic timeline reveals patterns that individual documents cannot. It shows the 14-month gap between Sarah Mitchell's first ANA test and her lupus diagnosis. It shows the 3 years Maya Rodriguez spent being told her symptoms were anxiety before POTS was diagnosed. These patterns illuminate the diagnostic odyssey problem.

**Example questions it helps answer:**
- "How long did it take from first symptom to diagnosis for this patient?"
- "What was the sequence of specialist referrals?"
- "Were there any dismissed symptoms that turned out to be significant?"

#### 13. autoimmune_cross_disease (weight: 0.02)

**What it contains:** Overlap syndrome data documenting connections between co-occurring autoimmune diseases. Records describe which diseases commonly co-occur, shared genetic risk factors, shared autoantibodies, and implications for treatment.

**Why it matters:** Autoimmune diseases rarely occur in isolation. A patient with Sjogren's syndrome has a higher risk of developing lymphoma. A patient with Hashimoto's thyroiditis has a higher risk of developing celiac disease. Understanding these connections helps clinicians anticipate complications and screen proactively.

**Example questions it helps answer:**
- "What conditions commonly co-occur with Sjogren's?"
- "What is the overlap between lupus and antiphospholipid syndrome?"
- "Should I screen for celiac disease in a Type 1 Diabetes patient?"

#### 14. genomic_evidence (weight: 0.02)

**What it contains:** Genomic variant data from the HCLS AI Factory genomics pipeline. This is a shared, read-only collection that contains millions of annotated variants from whole-genome sequencing.

**Why it matters:** Genomic data provides the deepest layer of patient-specific information. Variants in immune-related genes (HLA, cytokine receptors, complement components) can explain why a patient developed a specific autoimmune disease or why they respond differently to treatment.

**Example questions it helps answer:**
- "Are there pathogenic variants in complement genes for this patient?"
- "What immune-related variants were identified in the WGS data?"
- "Do any identified variants affect drug metabolism?"

### Mapping questions to collections

When you ask a question, the system searches all 14 collections. But certain types of questions draw more heavily from specific collections:

| Question Type | Primary Collections |
|--------------|-------------------|
| "What did the doctor say about...?" | clinical_documents, patient_labs |
| "What does this antibody mean?" | autoantibody_panels, disease_criteria |
| "What is the genetic risk?" | hla_associations, genomic_evidence |
| "How active is the disease?" | disease_activity, patient_labs |
| "Is a flare coming?" | flare_patterns, patient_labs |
| "What treatment should we try?" | biologic_therapies, pgx_rules, clinical_trials |
| "What does the literature say?" | literature, clinical_trials |
| "What is the patient's history?" | patient_timelines, clinical_documents |

---

## Chapter 7: The Knowledge Base

### What it is

The knowledge base (version 2.0.0, last updated 2026-03-10) is a hand-curated, structured repository of autoimmune domain knowledge. Unlike the vector collections (which contain embedded text searchable by semantic similarity), the knowledge base contains precise, structured data: specific HLA alleles with exact odds ratios, specific autoantibodies with known disease associations, specific scoring thresholds with exact cutoff values.

The knowledge base does not replace vector search -- it augments it. When the RAG engine finds relevant evidence in the vector collections, it also checks the knowledge base for structured facts that can enrich the response. This combination of unstructured retrieval (vector search) and structured knowledge produces answers that are both comprehensive and precise.

### The five knowledge domains

#### 1. HLA-Disease Associations (22 alleles)

The system maps 22 HLA alleles to their associated autoimmune diseases, complete with odds ratios and literature references. An odds ratio of 87.4 for HLA-B*27:05 and ankylosing spondylitis means that a person carrying this allele is 87 times more likely to develop AS than a person without it.

Key entries include:
- HLA-B*27:05 -- AS (OR=87.4), reactive arthritis (OR=20.0)
- HLA-DRB1*04:01 -- RA (OR=4.2, shared epitope hypothesis)
- HLA-DRB1*03:01 -- SLE (OR=2.4), Sjogren's (OR=3.1), T1D (OR=3.6), Graves' (OR=2.2), celiac (OR=7.0)
- HLA-DRB1*15:01 -- MS (OR=3.1)
- HLA-C*06:02 -- psoriasis (OR=10.0)
- HLA-DQB1*02:01 -- celiac disease (OR=7.0)

#### 2. Autoantibody-Disease Map (24 antibodies)

Each of the 24 autoantibodies is mapped to its associated diseases, with specificity and sensitivity data where available. For example:
- Anti-dsDNA is >95% specific for SLE and correlates with nephritis activity
- Anti-CCP is >95% specific for RA and can be positive years before clinical onset
- AChR antibodies are ~85% sensitive for generalized myasthenia gravis
- Anti-tTG (IgA) is >95% sensitive and >95% specific for celiac disease

#### 3. Disease Activity Thresholds (20 scores)

Twenty scoring systems are defined with their component variables, calculation methods, and threshold classifications. The system uses these thresholds to interpret calculated scores: a DAS28-CRP of 4.5 is automatically classified as "moderate activity," triggering consideration of treatment escalation.

#### 4. Flare Biomarker Patterns (13 diseases)

Each of the 13 diseases has a defined set of biomarker patterns that predict flares. For SLE, the classic pattern is rising anti-dsDNA + falling C3/C4 + rising ESR with stable CRP. For RA, it is rising CRP + rising ESR + rising RF titer. These patterns enable early warning before clinical symptoms become severe.

#### 5. Biologic Therapies (22 drugs with PGx)

Twenty-two biologic therapies are profiled with mechanism of action, approved indications, pharmacogenomic considerations, key trial evidence, monitoring requirements, and common adverse effects. Drugs include TNF inhibitors (adalimumab, etanercept, infliximab, certolizumab, golimumab), IL-6 blockers (tocilizumab, sarilumab), B-cell depleters (rituximab, belimumab), JAK inhibitors (tofacitinib, baricitinib, upadacitinib), and many more.

### How knowledge augments vector search

When a user asks "What does David Park's HLA-B27 result mean for his back pain?", the system:

1. Searches all 14 vector collections for relevant evidence (finding David's clinical documents, lab results, and imaging reports)
2. Detects "HLA-B27" in the query and retrieves structured knowledge: HLA-B*27:05 is associated with ankylosing spondylitis (OR=87.4) and reactive arthritis (OR=20.0)
3. Adds disease criteria context: modified New York criteria for AS require radiographic sacroiliitis plus at least one clinical criterion
4. Combines unstructured evidence (David's actual records) with structured knowledge (HLA statistics, diagnostic criteria) to produce a precise, well-grounded answer

---

## Chapter 8: The Five Clinical Engines

The Precision Autoimmune Agent is not just a search engine. It contains five specialized clinical analysis engines, implemented in `src/agent.py` (437 lines) and `src/diagnostic_engine.py` (519 lines), that perform structured clinical reasoning.

### 1. Autoantibody Interpretation

**What it does:** Takes an autoantibody panel (list of antibodies tested with results) and interprets each result in the context of autoimmune disease. It maps positive antibodies to associated diseases, evaluates ANA patterns, and identifies antibody combinations that strengthen or narrow the differential diagnosis.

**How it works:** The engine consults the AUTOANTIBODY_DISEASE_MAP from the knowledge base. A positive anti-dsDNA alone suggests SLE. A positive anti-dsDNA plus positive anti-Smith is highly specific for SLE. A positive ANA with centromere pattern suggests limited systemic sclerosis (CREST syndrome). The engine evaluates these combinations and produces a ranked list of disease associations with confidence levels.

**Clinical value:** Autoantibody interpretation requires expertise that many non-rheumatologists lack. A PCP who orders an ANA panel may not know what a speckled pattern at 1:640 titer means. This engine provides that interpretation instantly.

### 2. HLA Association Analysis

**What it does:** Takes HLA typing results and evaluates genetic susceptibility to autoimmune diseases. It calculates risk based on odds ratios from the knowledge base and considers multi-allele interactions.

**How it works:** The engine consults HLA_DISEASE_ASSOCIATIONS. For each HLA allele in the patient's profile, it retrieves all known disease associations with odds ratios. When multiple risk alleles are present (e.g., HLA-DRB1*03:01 + HLA-DRB1*04:01 for RA), it evaluates the compound risk.

**Clinical value:** HLA typing is increasingly available from commercial genetic tests, but interpreting the results in the context of autoimmune risk requires specialized knowledge. This engine bridges that gap.

### 3. Disease Activity Scoring

**What it does:** Calculates disease activity scores from laboratory values and clinical parameters. It classifies the result as remission, low, moderate, or high activity using the thresholds from DISEASE_ACTIVITY_THRESHOLDS.

**How it works:** The engine takes input parameters (e.g., tender joint count, swollen joint count, CRP, patient global assessment for DAS28-CRP) and applies the scoring formula. It then compares the result against defined thresholds and returns both the numeric score and the classification.

**Clinical value:** Activity scores guide treatment decisions. A DAS28-CRP > 5.1 (high activity) in a patient already on methotrexate is a strong signal to escalate to biologic therapy. The engine calculates this automatically and flags the clinical implication.

### 4. Flare Prediction

**What it does:** Analyzes longitudinal biomarker trends to assess flare risk. It compares current biomarker trajectories against the flare patterns defined in FLARE_BIOMARKER_PATTERNS and assigns a risk level (low, moderate, or high).

**How it works:** The engine retrieves a patient's lab results over time, calculates trends (rising, falling, stable), and matches the trend pattern against known flare signatures. For SLE, if anti-dsDNA is rising AND C3 is falling AND C4 is falling, the engine classifies flare risk as high.

**Clinical value:** Flare prediction is one of the most valuable applications of longitudinal data analysis. Catching a flare 4-8 weeks before clinical deterioration allows preemptive treatment adjustment, potentially preventing organ damage.

### 5. Biologic Therapy Recommendation

**What it does:** Recommends biologic therapies based on the patient's diagnosis, disease activity level, prior treatment history, HLA type, and pharmacogenomic data. Recommendations are ranked by evidence strength and include monitoring requirements.

**How it works:** The engine consults BIOLOGIC_THERAPIES and AUTOIMMUNE_PGX_RULES. It filters therapies by approved indication, removes any that are contraindicated based on PGx data, ranks the remainder by evidence strength and disease activity level, and returns a prioritized list with rationale.

**Clinical value:** Biologic therapy selection involves balancing efficacy, safety, PGx considerations, insurance formulary position, and patient preference. The engine provides an evidence-based starting point for that decision, ensuring that PGx considerations are never overlooked.

---

## Chapter 9: Setting Up Locally

### Prerequisites

Before you begin, ensure you have:

- **Python 3.10+** installed
- **Docker** installed and running (for Milvus)
- **Milvus** running on port 19530 (the HCLS AI Factory `start-factory.sh` script starts Milvus automatically; alternatively, use the Docker Compose file)
- **An Anthropic API key** for Claude LLM features (optional -- the system runs in degraded mode without it, providing vector search without LLM synthesis)
- **4+ GB RAM** available (Milvus + embedding model + application)

### Step-by-step setup

**Step 1: Navigate to the agent directory**

```bash
cd ai_agent_adds/precision_autoimmune_agent
```

**Step 2: Create a `.env` file**

```bash
cat > .env << 'EOF'
ANTHROPIC_API_KEY=sk-ant-your-key-here
MILVUS_HOST=localhost
MILVUS_PORT=19530
AUTO_STREAMLIT_PORT=8531
AUTO_API_PORT=8532
EMBEDDING_MODEL=BAAI/bge-small-en-v1.5
LOG_LEVEL=INFO
EOF
```

**Step 3: Run the setup script (creates collections and seeds knowledge)**

```bash
./run.sh --setup
```

This creates all 14 Milvus collections with their schemas and indexes, and seeds the knowledge base entries (HLA associations, autoantibodies, disease criteria, activity thresholds, flare patterns, and biologic therapies).

**Step 4: Load demo patient data**

```bash
./venv/bin/python scripts/generate_demo_patients.py
```

This ingests the demo PDFs for all 9 patients (Sarah Mitchell, Maya Rodriguez, David Park, Linda Chen, Rachel Thompson, Emma Williams, James Cooper, Karen Foster, Michael Torres) into the appropriate collections.

**Step 5: Start the system**

```bash
# Start UI only (port 8531)
./run.sh

# Start API only (port 8532)
./run.sh --api

# Start both UI and API
./run.sh --both
```

**Step 6: Open the UI**

Navigate to `http://localhost:8531` in your browser.

### Using run.sh modes

| Mode | Command | What it starts |
|------|---------|---------------|
| Default | `./run.sh` | Streamlit UI on port 8531 |
| API only | `./run.sh --api` | FastAPI on port 8532 (2 workers) |
| Both | `./run.sh --both` | UI + API with graceful shutdown |
| Setup | `./run.sh --setup` | Collection creation + knowledge seeding |

---

## Chapter 10: Exploring the API

### What the REST API is

The FastAPI REST API (port 8532) provides programmatic access to all of the agent's capabilities. You can query, search, analyze patients, run differential diagnosis, ingest documents, manage collections, and export results -- all without opening the UI.

The API is defined in `api/main.py` (583 lines) and provides 14 endpoints.

### The 14 endpoints

| Method | Path | Purpose |
|--------|------|---------|
| GET | `/` | Root -- returns API info and version |
| GET | `/health` | Health check with service status |
| GET | `/healthz` | Kubernetes-compatible health probe |
| POST | `/query` | RAG-powered clinical query (returns JSON) |
| POST | `/query/stream` | Streaming RAG query (returns SSE stream) |
| POST | `/search` | Vector search only (no LLM synthesis) |
| POST | `/analyze` | Full five-stage clinical analysis pipeline |
| POST | `/differential` | Differential diagnosis from symptoms/labs |
| POST | `/ingest/upload` | Upload and ingest a clinical document (PDF) |
| POST | `/ingest/demo-data` | Load demo patient data into collections |
| GET | `/collections` | List all collections with record counts |
| POST | `/collections/create` | Create collections and seed knowledge |
| POST | `/export` | Export analysis results (Markdown, FHIR, PDF) |
| GET | `/metrics` | Prometheus-compatible metrics |

### Testing with curl

Below are examples you can run from your terminal. These assume the API is running on `localhost:8532`.

#### Check health

```bash
curl -s http://localhost:8532/health | python3 -m json.tool
```

Expected response:

```json
{
    "status": "healthy",
    "services": {
        "milvus": "OK",
        "embedder": "OK",
        "llm": "OK"
    },
    "uptime_seconds": 3421.7,
    "version": "1.0.0"
}
```

#### Ask a clinical question

```bash
curl -s -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -d '{
    "query": "What does a positive anti-dsDNA with low complement suggest?",
    "patient_id": "sarah_mitchell",
    "top_k": 5
  }' | python3 -m json.tool
```

The response includes the LLM-generated answer, evidence items with collection sources and scores, knowledge base context, and timing metadata.

#### Search without LLM synthesis

```bash
curl -s -X POST http://localhost:8532/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": "HLA-B27 ankylosing spondylitis sacroiliac",
    "collections": ["autoimmune_hla_associations", "autoimmune_clinical_documents"],
    "top_k": 3
  }' | python3 -m json.tool
```

This returns raw vector search results without LLM synthesis -- useful for debugging or building custom pipelines.

#### Run a full patient analysis

```bash
curl -s -X POST http://localhost:8532/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sarah_mitchell"
  }' | python3 -m json.tool
```

This runs the complete five-stage pipeline: autoantibody interpretation, HLA analysis, disease activity scoring, flare prediction, and biologic therapy recommendation.

#### Generate a differential diagnosis

```bash
curl -s -X POST http://localhost:8532/differential \
  -H "Content-Type: application/json" \
  -d '{
    "symptoms": ["joint pain", "fatigue", "malar rash", "photosensitivity"],
    "labs": {
        "ANA": "positive, 1:640, homogeneous",
        "anti-dsDNA": "positive, 120 IU/mL",
        "C3": "low, 62 mg/dL",
        "C4": "low, 8 mg/dL"
    }
  }' | python3 -m json.tool
```

#### List collections

```bash
curl -s http://localhost:8532/collections | python3 -m json.tool
```

Expected response (abbreviated):

```json
{
    "collections": [
        {
            "name": "autoimmune_clinical_documents",
            "record_count": 1247,
            "weight": 0.18,
            "status": "loaded"
        },
        {
            "name": "autoimmune_patient_labs",
            "record_count": 892,
            "weight": 0.14,
            "status": "loaded"
        }
    ],
    "total_collections": 14,
    "total_records": 8421
}
```

#### Export results

```bash
curl -s -X POST http://localhost:8532/export \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sarah_mitchell",
    "format": "fhir_r4"
  }' --output sarah_mitchell_export.json
```

Supported formats: `markdown`, `fhir_r4`, `pdf`.

---

## Chapter 11: Understanding the Codebase

### Project structure

```
precision_autoimmune_agent/
|-- api/
|   |-- __init__.py
|   |-- main.py               (583 lines)  FastAPI server, 14 endpoints
|   +-- routes/
|       +-- __init__.py
|-- app/
|   |-- __init__.py
|   +-- autoimmune_ui.py      (1,160 lines) Streamlit UI, 10 tabs
|-- config/
|   |-- __init__.py
|   |-- logging.py             Centralized logging config
|   +-- settings.py            Pydantic settings from .env
|-- src/
|   |-- __init__.py
|   |-- agent.py               (437 lines)  Main orchestrator, 5 clinical engines
|   |-- rag_engine.py          (597 lines)  RAG engine, weighted parallel search
|   |-- collections.py         (562 lines)  14 Milvus collection schemas + manager
|   |-- knowledge.py           (855 lines)  Knowledge base v2.0.0, 5 domains
|   |-- diagnostic_engine.py   (519 lines)  Diagnostic pipeline, differential Dx
|   |-- document_processor.py  (435 lines)  PDF ingestion, chunking, embedding
|   |-- export.py              (389 lines)  Markdown, FHIR R4, PDF export
|   |-- models.py              (238 lines)  Pydantic models and data classes
|   +-- timeline_builder.py    (251 lines)  Patient timeline construction
|-- scripts/
|   |-- setup_collections.py   Collection creation and seeding
|   |-- generate_demo_patients.py  Demo data generation
|   |-- patient_sarah.py       Sarah Mitchell's 35 PDFs
|   |-- patient_maya.py        Maya Rodriguez's 28 PDFs
|   |-- pdf_engine.py          PDF generation engine
|   +-- ...                    (additional patient scripts)
|-- tests/
|   |-- __init__.py
|   |-- test_api.py
|   |-- test_autoimmune.py
|   |-- test_collections.py
|   |-- test_diagnostic_engine.py
|   +-- ...                    (455 tests across 8 files)
|-- demo_data/
|   |-- sarah_mitchell/        (35 PDFs spanning 3.5 years)
|   |-- maya_rodriguez/        (28 PDFs spanning 4 years)
|   |-- david_park/            (26 PDFs spanning 6 years)
|   +-- ...                    (additional patient folders)
|-- docs/
|   |-- API_REFERENCE.md
|   |-- ARCHITECTURE_GUIDE.md
|   |-- DEMO_GUIDE.md
|   |-- PROJECT_BIBLE.md
|   +-- PRECISION_AUTOIMMUNE_AGENT_RESEARCH_PAPER.md
|-- run.sh                     Startup script (4 modes)
|-- pyproject.toml             Project metadata and dependencies
+-- .env                       Configuration (not committed)
```

### File-by-file walkthrough of core files

#### src/agent.py (437 lines) -- The Orchestrator

This is the central coordinator. The `AutoimmuneAgent` class ties together the five clinical engines:

- `analyze_patient()` runs the full five-stage pipeline on a patient profile
- `interpret_autoantibodies()` maps antibody results to disease associations
- `analyze_hla_profile()` evaluates genetic susceptibility from HLA typing
- `score_disease_activity()` calculates activity scores and classifies severity
- `predict_flare()` analyzes biomarker trends against known flare patterns
- `recommend_biologics()` ranks therapy options considering PGx data

The agent imports its knowledge from `knowledge.py` (the five domain dictionaries) and its data models from `models.py`.

#### src/rag_engine.py (597 lines) -- The RAG Engine

This is the retrieval-augmented generation engine. The `AutoimmuneRAGEngine` class handles:

- Embedding queries using BGE-small-en-v1.5
- Parallel multi-collection search with `ThreadPoolExecutor`
- Collection-specific weight application
- Result deduplication and ranking
- Knowledge base context retrieval
- Prompt construction with evidence formatting
- LLM synthesis via Claude API
- Streaming response generation
- Conversation memory (using a deque for history)

Key data classes: `SearchHit` (individual result with collection, score, text, relevance tag) and `CrossCollectionResult` (aggregated results from all collections with timing).

#### src/collections.py (562 lines) -- The Collection Manager

This file defines the schemas for all 14 Milvus collections and provides the `AutoimmuneCollectionManager` class for creating, loading, searching, and managing collections.

Key constants:
- `INDEX_PARAMS`: COSINE metric, IVF_FLAT index, nlist=1024
- `SEARCH_PARAMS`: COSINE metric, nprobe=16
- `_DIM = 384` (BGE-small-en-v1.5 embedding dimension)

Each collection has a unique schema defined with helper functions (`_pk`, `_embedding`, `_varchar`, `_int`, `_float`) and registered via `_register()`.

#### src/knowledge.py (855 lines) -- The Knowledge Base

The largest source file, containing all structured domain knowledge:

- `KNOWLEDGE_VERSION`: Version metadata (v2.0.0), source references, and statistics
- `HLA_DISEASE_ASSOCIATIONS`: 22 alleles mapped to diseases with odds ratios
- `AUTOANTIBODY_DISEASE_MAP`: 24 antibodies mapped to diseases with sensitivity/specificity
- `DISEASE_ACTIVITY_THRESHOLDS`: 20 scoring systems with thresholds
- `FLARE_BIOMARKER_PATTERNS`: 13 disease-specific flare prediction patterns
- `BIOLOGIC_THERAPIES`: 22 drugs with mechanism, indications, PGx, and monitoring

This file is essentially a clinical reference database encoded as Python dictionaries. It is hand-curated from published guidelines and peer-reviewed literature.

#### src/diagnostic_engine.py (519 lines) -- The Diagnostic Engine

The `DiagnosticEngine` class combines the agent's clinical engines with the RAG engine's retrieval capabilities to perform comprehensive diagnostic analysis:

- `run_differential()` generates a differential diagnosis from symptoms and labs
- `evaluate_criteria()` checks patient data against classification criteria
- `comprehensive_analysis()` runs a full workup combining vector search, knowledge base, and clinical engines

It uses the `AutoimmuneAgent` for clinical reasoning and the `AutoimmuneRAGEngine` for evidence retrieval, combining both into structured diagnostic reports.

#### src/document_processor.py (435 lines) -- The Document Processor

The `DocumentProcessor` class handles the ingestion pipeline:

1. **PDF text extraction** -- Extracts text from uploaded PDF files
2. **Semantic chunking** -- Splits text into chunks of up to 2,500 characters with 200-character overlap, respecting sentence boundaries
3. **Entity extraction** -- Identifies medical entities (diseases, drugs, labs, antibodies) in each chunk
4. **Embedding** -- Converts each chunk to a 384-dimensional vector using BGE-small-en-v1.5
5. **Milvus insertion** -- Inserts embedded chunks into the appropriate collection

The chunking parameters (2,500 chars, 200 overlap) are chosen to balance between capturing enough context in each chunk and keeping chunks small enough for precise retrieval.

#### src/export.py (389 lines) -- The Export Engine

Supports three export formats:

- **Markdown** -- Human-readable reports with formatted evidence and citations
- **FHIR R4** -- Standards-compliant FHIR R4 Bundles (DiagnosticReport, Observation, Condition resources) for integration with healthcare systems
- **PDF** -- Formatted PDF reports suitable for clinical documentation

#### src/models.py (238 lines) -- Data Models

Pydantic models and dataclasses that define the data structures used throughout the system:

- `AutoimmunePatientProfile` -- Patient demographics, diagnoses, labs, antibodies, HLA
- `AutoantibodyPanel` -- Autoantibody test results
- `HLAProfile` -- HLA typing results
- `DiseaseActivityScore` / `DiseaseActivityLevel` -- Activity score with classification
- `FlarePredictor` / `FlareRisk` -- Flare risk assessment
- `BiologicTherapy` -- Therapy recommendation
- `AutoimmuneAnalysisResult` -- Complete analysis pipeline output

#### api/main.py (583 lines) -- The FastAPI Server

The API server with lifecycle management, middleware (authentication, timing, request size limiting), and 14 endpoint handlers. Key architectural decisions:

- Global state dictionary (`_state`) holds initialized components (collection manager, embedder, LLM client, RAG engine, agent, document processor, diagnostic engine, timeline builder)
- Lifespan context manager handles startup validation (Milvus connection, embedder loading, LLM client initialization) with graceful degradation
- Optional API key authentication (skip if not configured)
- Request timing header on every response

#### app/autoimmune_ui.py (1,160 lines) -- The Streamlit UI

The largest file in the project. Implements 10 tabs:

1. **Chat** -- Free-text clinical queries with streaming LLM responses
2. **Patient Analysis** -- Five-stage pipeline with visual results
3. **Differential Dx** -- Symptom/lab input with ranked differential
4. **Labs** -- Lab result exploration with trending charts
5. **Timeline** -- Chronological diagnostic journey visualization
6. **Collections** -- Collection management with record counts
7. **Knowledge** -- Interactive knowledge base browser
8. **Export** -- Multi-format export with preview
9. **Settings** -- Configuration management
10. **Help** -- Documentation and example queries

### How the pieces connect: a trace through a query

Let us trace what happens when a user types "What does Sarah's anti-dsDNA trend suggest?" in the Chat tab:

```
1. autoimmune_ui.py receives the input in the Chat tab
   -> Calls the RAG engine's query method

2. rag_engine.py embeds the question
   -> BGE-small-en-v1.5 produces a 384-dim vector
   -> Detects "anti-dsDNA" and "Sarah" in the query

3. rag_engine.py searches all 14 collections in parallel
   -> collections.py executes Milvus searches via ThreadPoolExecutor
   -> Each collection returns up to 5 results (TOP_K=5)
   -> Results below 0.40 score are filtered (SCORE_THRESHOLD)

4. rag_engine.py retrieves knowledge base context
   -> knowledge.py provides anti-dsDNA disease associations,
      SLE flare patterns, and relevant activity score context

5. rag_engine.py merges, deduplicates, and ranks results
   -> Applies collection weights (clinical_documents=0.18, labs=0.14, etc.)
   -> Assigns relevance tags (HIGH >= 0.80, MEDIUM >= 0.60, LOW < 0.60)
   -> Caps total evidence at 30 items (MAX_EVIDENCE)

6. rag_engine.py builds the prompt
   -> Formats evidence by collection section
   -> Appends knowledge base context
   -> Includes the original question

7. rag_engine.py calls Claude API
   -> Streams the response token by token

8. autoimmune_ui.py displays the result
   -> Shows streaming response with citations
   -> Shows expandable evidence panel
   -> Shows download buttons (Markdown, FHIR R4, PDF)
```

---

## Chapter 12: Next Steps

### Where to go from here

Now that you understand the foundations, here are paths forward based on your role:

**For clinicians:**
- Try all 9 demo patients and compare the system's analysis to your clinical judgment
- Use the Differential Dx tab with your own symptom/lab combinations
- Explore the Knowledge tab to browse HLA associations, autoantibodies, and biologic therapies
- Test edge cases: overlap syndromes, atypical presentations, rare antibody patterns

**For data scientists and ML engineers:**
- Explore the embedding space by examining similarity scores across collections
- Experiment with different TOP_K and SCORE_THRESHOLD values
- Analyze the collection weight distribution and how it affects result ranking
- Build custom analysis pipelines using the `/search` API endpoint
- Review the 455 tests across 8 test files for validation coverage

**For software developers:**
- Read the full API Reference (`docs/API_REFERENCE.md`) for endpoint details
- Trace the code from `api/main.py` through `rag_engine.py` to `collections.py`
- Understand the Milvus schema design in `collections.py`
- Explore the document processing pipeline in `document_processor.py`
- Study the export formats (especially FHIR R4) in `export.py`

### The Advanced Learning Guide

The Advanced Learning Guide (a companion to this document) covers:

- Embedding model fine-tuning for autoimmune domain
- Collection weight optimization using relevance feedback
- Custom knowledge base extension
- Multi-agent integration (connecting with the Biomarker Agent for inflammation monitoring and the Imaging Agent for joint assessment)
- Production deployment on DGX Spark
- Performance tuning (Milvus index parameters, batch embedding, caching)

### How to contribute

This project is open source under the Apache 2.0 license. Contributions are welcome:

1. **Data contributions:** Additional curated knowledge entries, clinical trial records, or literature references
2. **Code contributions:** Bug fixes, new clinical engines, improved chunking strategies, or additional export formats
3. **Testing contributions:** Additional test cases, edge case scenarios, or validation against clinical gold standards
4. **Documentation contributions:** Corrections, clarifications, or translations

---

## Glossary

| Term | Definition |
|------|-----------|
| **ACR** | American College of Rheumatology. Professional organization that publishes classification criteria for rheumatic diseases. |
| **AChR** | Acetylcholine Receptor. Antibodies against AChR are diagnostic for myasthenia gravis (~85% sensitivity in generalized MG). |
| **ANA** | Antinuclear Antibody. A screening test for autoimmune disease. Positive in ~95% of SLE patients, but also positive in 5-15% of healthy individuals. The pattern (homogeneous, speckled, nucleolar, centromere) provides diagnostic specificity. |
| **ANCA** | Anti-Neutrophil Cytoplasmic Antibody. Two patterns: c-ANCA (cytoplasmic, targeting PR3) associated with granulomatosis with polyangiitis; p-ANCA (perinuclear, targeting MPO) associated with microscopic polyangiitis and some IBD. |
| **Anti-CCP** | Anti-Cyclic Citrullinated Peptide. Highly specific (>95%) for rheumatoid arthritis. Can be positive years before clinical onset. |
| **Anti-dsDNA** | Anti-double-stranded DNA. Highly specific (>95%) for SLE. Titers correlate with disease activity, especially lupus nephritis. |
| **Anti-Smith** | Antibodies against the Smith nuclear antigen. Highly specific for SLE (~99%) but low sensitivity (~25%). |
| **Anti-SSA/Ro** | Antibodies associated with Sjogren's syndrome, neonatal lupus, and subacute cutaneous lupus. Can cross the placenta and cause congenital heart block. |
| **Anti-SSB/La** | Antibodies typically co-occurring with anti-SSA. More specific for Sjogren's syndrome than anti-SSA alone. |
| **Anti-tTG** | Anti-tissue Transglutaminase. IgA anti-tTG is >95% sensitive and >95% specific for celiac disease. First-line serologic test. |
| **ASDAS** | Ankylosing Spondylitis Disease Activity Score. Uses CRP, back pain, morning stiffness, peripheral pain/swelling, and patient global assessment. |
| **BASDAI** | Bath Ankylosing Spondylitis Disease Activity Index. Six-question patient-reported outcome measure for AS activity. Score 0-10. |
| **BGE-small-en-v1.5** | BAAI General Embedding model (small, English, version 1.5). Produces 384-dimensional embeddings. Used for all text embedding in this system. |
| **Biologic therapy** | A class of drugs derived from living organisms that target specific components of the immune system (e.g., TNF, IL-6, B-cells, T-cells). |
| **CDAI** | Clinical Disease Activity Index. A simplified RA activity score that does not require laboratory values. Remission <= 2.8. |
| **Celiac Disease** | Autoimmune condition triggered by gluten in genetically susceptible individuals (HLA-DQ2/DQ8). Causes villous atrophy of the small intestine. |
| **Complement (C3/C4)** | Proteins in the complement cascade. Low C3 and C4 levels indicate complement consumption, commonly seen in active SLE, especially lupus nephritis. |
| **Cosine similarity** | A mathematical measure of the angle between two vectors, ranging from 0 (orthogonal/unrelated) to 1 (identical direction/meaning). Used for vector similarity search. |
| **CPIC** | Clinical Pharmacogenetics Implementation Consortium. Publishes peer-reviewed guidelines for genotype-to-phenotype translation and drug dosing. |
| **CRP** | C-Reactive Protein. An acute-phase reactant that rises with inflammation. Used in multiple disease activity scores (DAS28-CRP, ASDAS). |
| **DAS28** | Disease Activity Score using 28 joints. The standard composite measure for RA activity. Variants: DAS28-CRP and DAS28-ESR. Remission < 2.6. |
| **EDSS** | Expanded Disability Status Scale. The standard measure of MS disability. Ranges from 0 (normal) to 10 (death due to MS). |
| **Embedding** | A numerical representation of text as a vector (list of numbers) such that semantically similar texts have similar vectors. |
| **ESR** | Erythrocyte Sedimentation Rate. A nonspecific marker of inflammation. Used in DAS28-ESR, ASDAS, and as a general inflammation indicator. |
| **ESSDAI** | EULAR Sjogren's Syndrome Disease Activity Index. Measures systemic disease activity across 12 organ domains. |
| **EULAR** | European Alliance of Associations for Rheumatology. Co-publishes classification criteria and treatment recommendations with ACR. |
| **FHIR R4** | Fast Healthcare Interoperability Resources, Release 4. An international standard for exchanging healthcare information electronically. |
| **Flare** | A period of increased disease activity in a chronic autoimmune condition, often preceded by characteristic biomarker changes. |
| **Graves' Disease** | Autoimmune thyroid disease caused by thyroid-stimulating immunoglobulins (TSI) that activate the TSH receptor, causing hyperthyroidism. |
| **Hashimoto's Thyroiditis** | Autoimmune thyroid disease caused by anti-TPO and anti-thyroglobulin antibodies that destroy thyroid tissue, causing hypothyroidism. |
| **HLA** | Human Leukocyte Antigen. A gene complex encoding cell-surface proteins critical for immune regulation. Specific HLA alleles confer risk for specific autoimmune diseases. |
| **IBD** | Inflammatory Bowel Disease. Encompasses Crohn's disease (can affect any part of the GI tract) and ulcerative colitis (affects the colon). |
| **IVF_FLAT** | Inverted File with Flat quantization. A Milvus indexing algorithm that partitions vectors into clusters for efficient approximate nearest neighbor search. |
| **JAK inhibitor** | Janus Kinase inhibitor. A class of oral small-molecule drugs (tofacitinib, baricitinib, upadacitinib) that block intracellular signaling pathways involved in inflammation. |
| **LLM** | Large Language Model. An AI model trained on vast text data capable of understanding and generating human language. Claude is the LLM used in this system. |
| **Milvus** | An open-source vector database designed for similarity search over large-scale embedding data. Used to store and search all 14 collections. |
| **mRSS** | Modified Rodnan Skin Score. A clinical measure of skin thickness in systemic sclerosis. Ranges from 0 to 51 (17 body areas scored 0-3). |
| **Myasthenia Gravis** | Autoimmune disease affecting the neuromuscular junction, caused by antibodies against acetylcholine receptors (AChR) or MuSK. Causes fluctuating weakness. |
| **nprobe** | A Milvus search parameter that controls how many index partitions (clusters) are searched. Higher nprobe = more thorough search but slower. This system uses nprobe=16. |
| **PASI** | Psoriasis Area and Severity Index. Measures the extent and severity of psoriasis based on erythema, induration, and scaling across four body regions. |
| **PGx** | Pharmacogenomics. The study of how genetic variation affects drug response. Used to guide dosing and prevent adverse reactions (e.g., TPMT testing before azathioprine). |
| **RAG** | Retrieval-Augmented Generation. A pattern where relevant documents are retrieved from a database and provided as context to an LLM, enabling grounded, citation-backed answers. |
| **RF** | Rheumatoid Factor. An autoantibody (IgM against IgG Fc) associated with RA but also present in other autoimmune diseases and 5-10% of healthy elderly. |
| **SDAI** | Simplified Disease Activity Index. Similar to CDAI but adds CRP. Remission <= 3.3. |
| **Shared epitope** | A sequence motif in the HLA-DRB1 third hypervariable region (positions 70-74) that confers susceptibility to RA by binding citrullinated peptides. |
| **SLEDAI-2K** | Systemic Lupus Erythematosus Disease Activity Index 2000. A weighted scoring system for lupus activity. Score 0 = inactive; >10 = high activity. |
| **SLE (Lupus)** | Systemic Lupus Erythematosus. A chronic autoimmune disease that can affect virtually any organ system. Characterized by ANA positivity and multi-organ inflammation. |
| **TOP_K** | The number of results returned per collection in a vector search. This system uses TOP_K=5, meaning each of the 14 collections returns up to 5 results per query. |
| **SCORE_THRESHOLD** | The minimum cosine similarity score for a search result to be included. This system uses 0.40. Results below this threshold are discarded. |
| **TSI** | Thyroid-Stimulating Immunoglobulin. The causative autoantibody in Graves' disease. Mimics TSH and stimulates the thyroid to produce excess hormone. |
| **Vector** | A list of numbers representing a point in multi-dimensional space. In this system, a 384-dimensional vector represents the semantic meaning of a text chunk. |
| **Vector database** | A database optimized for storing and searching high-dimensional vectors using approximate nearest neighbor algorithms. Milvus is the vector database used here. |
| **VCF** | Variant Call Format. A standard file format for storing genomic variant data (SNPs, insertions, deletions). Used by the genomic_evidence collection. |

---

*This guide is part of the HCLS AI Factory documentation. For questions, issues, or contributions, visit the project repository.*
