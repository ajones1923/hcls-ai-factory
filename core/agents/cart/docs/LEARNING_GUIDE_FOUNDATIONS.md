# Learning Guide -- Foundations

**CAR-T Intelligence Agent | HCLS AI Factory**

Author: Adam Jones
Date: March 2026
License: Apache 2.0

---

## Welcome

You are reading the foundational learning guide for the CAR-T Intelligence Agent, an AI-powered research assistant that breaks down data silos across the entire CAR-T cell therapy development lifecycle. This system is part of the HCLS AI Factory, an end-to-end precision medicine platform that runs on a single NVIDIA DGX Spark ($4,699).

### Who this guide is for

This guide is written for three audiences:

- **Biologists and clinicians** who know CAR-T therapy but are new to AI/ML, vector databases, and retrieval-augmented generation.
- **Data scientists and ML engineers** who know embeddings and LLMs but are new to cell therapy, immunology, and clinical development.
- **Software developers** who want to understand how the system works end to end -- from Streamlit UI to Milvus vector search to Claude LLM synthesis.

You do not need to be an expert in all three areas. The whole point of this guide is to bring you up to speed on whichever parts are new to you.

### What you will learn

By the end of this guide, you will understand:

1. What CAR-T cell therapy is and why it matters
2. Why data fragmentation is the central challenge in CAR-T development
3. How Retrieval-Augmented Generation (RAG) works, from first principles
4. How this system searches 11 collections containing 3.5 million vectors simultaneously
5. How to use the UI to ask questions and interpret the answers
6. What each of the 11 data collections contains and why it exists
7. How the knowledge graph, query expansion, and comparative analysis features work
8. How to set up and run the system locally
9. How to use the REST API
10. How the codebase is organized, file by file

### Prerequisites

- Basic Python knowledge (you can read a Python function and understand what it does)
- A computer with a terminal
- Curiosity about either cancer immunotherapy or AI-assisted research (or both)

No prior knowledge of CAR-T biology, vector databases, or large language models is required. We will build every concept from the ground up.

---

## Chapter 1: What Is CAR-T Cell Therapy?

### The core idea

Your immune system already fights disease. White blood cells called T-cells patrol your body, looking for cells that display foreign markers on their surface. When a T-cell finds a match, it kills the target cell. This is how your body fights infections and, to some extent, cancer.

The problem is that cancer cells are sneaky. They look similar to normal cells, and they evolve ways to hide from T-cells. For decades, researchers have tried to help the immune system recognize cancer more effectively. CAR-T therapy is one of the most dramatic breakthroughs in that effort.

**CAR-T = Chimeric Antigen Receptor T-cell therapy.**

The word "chimeric" means "made from parts of different origins" -- like the chimera of Greek mythology, which had the head of a lion, the body of a goat, and the tail of a serpent. A CAR is an engineered protein assembled from parts of different immune molecules, stitched together into something nature never built.

### How it works, step by step

Imagine a patient with leukemia -- cancer of the blood. Here is how CAR-T therapy works:

```
Step 1: COLLECT
  A physician draws blood from the patient and separates out the T-cells
  using a process called leukapheresis. These T-cells are the raw
  material.

Step 2: ENGINEER
  In a laboratory, scientists use a viral vector (typically a modified,
  harmless virus) to insert a new gene into the T-cells. This gene
  encodes the CAR protein. Once the gene is active, the T-cell
  manufactures the CAR protein and displays it on its surface. The CAR
  acts like a synthetic antenna that can recognize a specific marker
  (antigen) on the cancer cell.

Step 3: EXPAND
  The engineered T-cells are grown in a bioreactor for 7 to 14 days,
  multiplying from millions to hundreds of millions or billions of cells.

Step 4: INFUSE
  The patient receives lymphodepletion chemotherapy (to make room for the
  new cells), and then the CAR-T cells are infused back into the
  patient's bloodstream.

Step 5: SEEK AND DESTROY
  The CAR-T cells circulate through the body. When a CAR-T cell
  encounters a cancer cell displaying the target antigen, it locks on
  and kills it. A single CAR-T cell can kill many cancer cells, and the
  cells can persist in the body for months or years.
```

### The anatomy of a CAR

A CAR protein has four main regions, stacked from the outside of the T-cell to the inside:

```
  OUTSIDE THE CELL
  ┌─────────────────────────┐
  │   scFv Binding Domain   │  Recognizes the target antigen (e.g., CD19)
  │   (from an antibody)    │  This is the "targeting" part
  ├─────────────────────────┤
  │   Hinge / Spacer        │  Provides flexibility and reach
  ├─────────────────────────┤
  │   Transmembrane Domain  │  Anchors the CAR in the cell membrane
  ├─────────────────────────┤
  │   Costimulatory Domain  │  Tells the T-cell to persist and multiply
  │   (CD28 or 4-1BB)       │  (like a sustain signal)
  ├─────────────────────────┤
  │   Signaling Domain      │  Tells the T-cell to activate and kill
  │   (CD3-zeta)            │  (like a go signal)
  └─────────────────────────┘
  INSIDE THE CELL
```

Think of it this way: the scFv is the eye that spots the enemy. The costimulatory domain is the stamina that keeps the soldier fighting. The signaling domain is the trigger that fires the weapon.

### The six FDA-approved CAR-T products

As of 2026, the U.S. FDA has approved six CAR-T products for cancer treatment:

| Product | Generic Name | Target | Costimulatory | Manufacturer | Initial Approval |
|---------|-------------|--------|---------------|--------------|-----------------|
| **Kymriah** | tisagenlecleucel | CD19 | 4-1BB | Novartis | Aug 2017 |
| **Yescarta** | axicabtagene ciloleucel | CD19 | CD28 | Kite/Gilead | Oct 2017 |
| **Tecartus** | brexucabtagene autoleucel | CD19 | CD28 | Kite/Gilead | Jul 2020 |
| **Breyanzi** | lisocabtagene maraleucel | CD19 | 4-1BB | BMS | Feb 2021 |
| **Abecma** | idecabtagene vicleucel | BCMA | 4-1BB | BMS | Mar 2021 |
| **Carvykti** | ciltacabtagene autoleucel | BCMA | 4-1BB | J&J/Legend | Feb 2022 |

Notice two patterns:

1. **Four of the six products target CD19**, a protein found on B-cells (a type of white blood cell that becomes cancerous in leukemia and lymphoma). The other two target BCMA, a protein found on plasma cells (which become cancerous in multiple myeloma).

2. **The costimulatory domain is either CD28 or 4-1BB.** This is not a trivial choice -- CD28-based CARs tend to expand rapidly and produce strong short-term responses, while 4-1BB-based CARs tend to persist longer. Choosing between them is one of the central design decisions in CAR engineering.

### Why CAR-T matters

CAR-T therapy can produce complete remissions in patients who have failed every other treatment. In the ELIANA trial, 82% of pediatric patients with relapsed B-ALL achieved complete remission with Kymriah. Some of these patients were told they had weeks to live. Some are still in remission years later.

But CAR-T therapy is not without challenges:

- **Cytokine Release Syndrome (CRS):** When CAR-T cells activate, they can trigger a massive inflammatory response. Symptoms range from fever to life-threatening organ failure. CRS occurs in 50-95% of patients (any grade) and requires careful management with tocilizumab and corticosteroids.

- **Neurotoxicity (ICANS):** Immune Effector Cell-Associated Neurotoxicity Syndrome can cause confusion, aphasia, seizures, and rarely, fatal cerebral edema. It occurs in 20-65% of patients.

- **Relapse:** Cancer can come back. One major mechanism is "antigen escape" -- the cancer cells stop expressing the target antigen (e.g., CD19), so the CAR-T cells can no longer find them.

- **Manufacturing complexity:** Each patient's T-cells must be individually collected, engineered, expanded, and shipped. The process takes 3-6 weeks, costs over $300,000, and sometimes fails.

- **Limited targets:** Most success has been in blood cancers with well-defined surface antigens. Solid tumors remain much harder because of hostile tumor microenvironments and the difficulty of finding truly tumor-specific targets.

These are the challenges that drive ongoing research -- and the challenges that this intelligence system is designed to help address.

---

## Chapter 2: The Data Challenge

### Where CAR-T knowledge lives

A researcher designing a new CAR-T construct needs to consider:

- **Published literature:** Over 100,000 papers on PubMed mention CAR-T or chimeric antigen receptor.
- **Clinical trials:** Thousands of CAR-T trials are registered on ClinicalTrials.gov, spanning dozens of targets and indications.
- **Patent filings:** Novel construct designs, manufacturing methods, and combination approaches are protected by thousands of patents.
- **Safety databases:** The FDA's FAERS (FDA Adverse Event Reporting System) contains post-market safety reports, including rare events like secondary malignancies.
- **Manufacturing records:** Process parameters (transduction efficiency, expansion fold, viability, VCN) are documented in regulatory filings, publications, and internal records.
- **Genomic data:** Patient-level variant data (from whole-genome sequencing) can reveal why certain patients respond differently.
- **Biomarker studies:** Predictive markers (ferritin, CRP, IL-6) and monitoring approaches (MRD, ctDNA) are published across many journals.
- **Regulatory decisions:** FDA approval letters, EMA decisions, label updates, and boxed warnings.
- **Real-world evidence:** Post-market registries like CIBMTR (Center for International Blood and Marrow Transplant Research) track outcomes in the real world, outside of controlled trials.
- **Sequence and structural data:** scFv binding domains, CDR sequences, binding affinities, and humanization strategies.

### The problem: data silos

Each of these data sources lives in a different database, uses a different format, and is searched with a different interface. A researcher who wants to answer a cross-functional question -- for example, "Why does my CAR-T construct fail?" -- would need to:

1. Search PubMed for resistance mechanisms
2. Search ClinicalTrials.gov for trial outcomes
3. Check FDA labels for known toxicities
4. Review manufacturing specifications
5. Look up assay results from pre-clinical testing
6. Check post-market safety databases
7. Consult biomarker literature for predictive signals

This is slow, error-prone, and nearly impossible to do comprehensively. A human researcher might spend days manually searching across these sources and still miss relevant information buried in a FAERS report or a registry analysis.

### The solution: one unified intelligence platform

The CAR-T Intelligence Agent solves this by:

1. **Ingesting** data from all of these sources into a single system
2. **Embedding** every piece of text as a 384-dimensional vector (a list of 384 numbers that captures the meaning of the text)
3. **Storing** these vectors in 11 purpose-built Milvus collections
4. **Searching** all 11 collections simultaneously when you ask a question
5. **Augmenting** the search results with structured knowledge from a hand-curated knowledge graph
6. **Synthesizing** a grounded answer using Claude, with citations back to the original sources

The result: you ask one question, and the system searches 3.5 million vectors across 11 data domains in under 20 milliseconds, then produces a cross-functional answer that cites PubMed papers, clinical trials, safety reports, biomarker data, and manufacturing parameters -- all in one response.

---

## Chapter 3: What Is RAG?

### The limitation of language models

Large Language Models (LLMs) like Claude are trained on vast amounts of text. They can write coherently, reason about complex topics, and follow instructions. But they have a fundamental limitation: **they do not have access to your specific data.**

If you ask Claude about CAR-T therapy, it can give you a general answer based on its training data. But it cannot tell you what is in your Milvus database. It cannot cite a specific PubMed paper. It cannot look up the binding affinity of the FMC63 scFv clone from your curated sequence records.

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

Think of it like this: imagine you are taking an open-book exam. The LLM is the student. RAG is the process of finding the right pages in the textbook (retrieval), putting them on the desk in front of the student (augmentation), and then asking the student to write an answer using those pages (generation).

Without RAG, the student is answering from memory (which may be outdated or vague). With RAG, the student is answering from evidence (which is specific, current, and citable).

### How retrieval works: embeddings and vector similarity

The retrieval step is the most technically interesting part, so let us break it down.

#### What is an embedding?

An embedding is a way of representing text as a list of numbers (a "vector") such that texts with similar meanings have similar numbers.

This system uses a model called **BGE-small-en-v1.5**, which converts any piece of text into a vector of **384 numbers** (384 dimensions). For example:

```
"CD19 CAR-T therapy for B-ALL"
    --> [0.023, -0.156, 0.891, 0.044, ..., -0.312]  (384 numbers)

"CAR-T cells targeting CD19 in acute lymphoblastic leukemia"
    --> [0.019, -0.148, 0.883, 0.051, ..., -0.298]  (384 numbers)
```

These two texts have very similar vectors because they express similar meanings, even though the words are different.

Conversely:

```
"Cryopreservation protocol for cell therapy products"
    --> [0.512, 0.078, -0.234, 0.667, ..., 0.112]  (384 numbers)
```

This vector looks very different from the CD19 vectors because the topic is different.

#### An analogy: GPS coordinates for meaning

Think of embeddings like GPS coordinates for meaning. Just as GPS coordinates place a physical location in a two-dimensional space (latitude and longitude), embeddings place a piece of text in a 384-dimensional meaning space. Texts about similar topics end up at nearby coordinates. Texts about unrelated topics end up far apart.

You cannot visualize 384 dimensions (nobody can), but the math works the same way as it does in two dimensions. To find texts similar to your question, you measure the "distance" between vectors.

#### How similarity search works

When you type a question into the CAR-T Intelligence Agent:

1. Your question is embedded into a 384-dimensional vector.
2. That vector is compared to every vector in the database using **cosine similarity** (a measure of how close two vectors are in direction, regardless of length).
3. The vectors with the highest similarity scores are returned as the most relevant results.

Cosine similarity ranges from 0 (completely unrelated) to 1 (identical meaning). In practice, scores above 0.75 indicate high relevance, and scores between 0.60 and 0.75 indicate moderate relevance.

The database that stores these vectors and performs fast similarity searches is **Milvus**, a purpose-built vector database. Milvus uses an indexing algorithm called **IVF_FLAT** (Inverted File with Flat quantization) to search millions of vectors in milliseconds rather than scanning them one by one.

### How the full RAG pipeline works in this system

Here is the complete pipeline, from question to answer:

```
User types: "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"
  |
  v
[1] EMBED THE QUESTION
    BGE-small-en-v1.5 converts the question to a 384-dim vector.
    The model prepends a special instruction prefix:
    "Represent this sentence for searching relevant passages: ..."
    This asymmetric prefix improves retrieval quality.
  |
  v
[2] SEARCH ALL 11 COLLECTIONS (in parallel)
    The query vector is sent to Milvus, which searches all 11
    collections simultaneously using ThreadPoolExecutor:
      - cart_literature (5,047 vectors)
      - cart_trials (973 vectors)
      - cart_constructs (41 vectors)
      - cart_assays (75 vectors)
      - cart_manufacturing (56 vectors)
      - cart_safety (71 vectors)
      - cart_biomarkers (60 vectors)
      - cart_regulatory (40 vectors)
      - cart_sequences (40 vectors)
      - cart_realworld (54 vectors)
      - genomic_evidence (3,561,170 vectors)

    Each collection returns its top-5 most similar results.
    Total search time: 12-16 ms (cached).
  |
  v
[3] QUERY EXPANSION
    The system detects "CD19" in the question and expands the
    search to include related terms: "B-ALL", "DLBCL", "FMC63",
    "Kymriah", "Yescarta", "tisagenlecleucel", etc.
    Additional searches are run with these expanded terms.
  |
  v
[4] KNOWLEDGE GRAPH AUGMENTATION
    The system detects "CD19" and adds structured knowledge:
      - Protein: B-Lymphocyte Antigen CD19 (UniProt P15391)
      - Approved products: Kymriah, Yescarta, Tecartus, Breyanzi
      - Key trials: ELIANA, ZUMA-1, ZUMA-2, TRANSFORM, TRANSCEND
      - Resistance: CD19 loss, lineage switch, trogocytosis
      - Toxicity profile: CRS 30-90%, ICANS 20-65%
  |
  v
[5] MERGE, DEDUPLICATE, AND RANK
    All results are merged, duplicates removed, and ranked by
    weighted score. Citation relevance is assigned:
      - Score >= 0.75: high relevance
      - Score >= 0.60: medium relevance
      - Score < 0.60: low relevance
  |
  v
[6] BUILD THE PROMPT
    The top-ranked evidence is formatted into a structured prompt:
      - Section per collection (Literature, Trial, Construct, etc.)
      - Each evidence item includes a clickable citation link
      - Knowledge graph context is appended
      - The user's original question is stated
  |
  v
[7] LLM SYNTHESIS
    Claude Sonnet 4.6 receives the prompt with the system instruction
    (a detailed persona prompt covering all 12 expertise domains)
    and generates a comprehensive, citation-rich answer.
    Response is streamed token by token to the UI.
  |
  v
[8] DISPLAY
    The Streamlit UI shows:
      - The streaming LLM response with clickable citations
      - An expandable evidence panel with collection badges
      - Citation relevance indicators (high/medium/low)
      - Download buttons (Markdown, JSON, PDF)
```

---

## Chapter 4: System Overview

### Architecture at a high level

The CAR-T Intelligence Agent has four main layers:

```
┌──────────────────────────────────────────────────┐
│                  USER INTERFACE                    │
│         Streamlit (port 8521) + FastAPI (8522)    │
│  Chat | Knowledge Graph Explorer | Image Analysis │
└───────────────────────┬──────────────────────────┘
                        |
                        v
┌──────────────────────────────────────────────────┐
│                   RAG ENGINE                      │
│  CARTRAGEngine (src/rag_engine.py)               │
│  - Query expansion (12 domain maps)               │
│  - Parallel multi-collection search               │
│  - Comparative analysis mode                      │
│  - Citation relevance scoring                     │
│  - Conversation memory                            │
└───────────────────────┬──────────────────────────┘
                        |
            ┌───────────┴───────────┐
            v                       v
┌─────────────────────┐  ┌────────────────────────┐
│   MILVUS VECTOR DB  │  │   KNOWLEDGE GRAPH      │
│   11 Collections    │  │   3 Domains:           │
│   3,567,622 vectors │  │   - 34 Targets         │
│   IVF_FLAT / COSINE │  │   - 17 Toxicities      │
│   384 dimensions    │  │   - 20 Manufacturing   │
│   (BGE-small)       │  │   - 23 Biomarkers      │
│                     │  │   - 6 Regulatory       │
│                     │  │   - 6 Immunogenicity   │
└─────────────────────┘  └────────────────────────┘
            |
            v
┌──────────────────────────────────────────────────┐
│                  CLAUDE LLM                       │
│  Claude Sonnet 4.6 (Anthropic API)               │
│  System prompt: 12-domain CAR-T expert persona   │
│  Streaming token generation                      │
└──────────────────────────────────────────────────┘
```

### How the pieces connect

1. **The user** types a question in the Streamlit chat interface (or sends a POST to the FastAPI REST API).
2. **The RAG engine** embeds the question using BGE-small-en-v1.5, searches all 11 Milvus collections in parallel, expands the query using 12 domain-specific expansion maps, and retrieves knowledge graph context.
3. **Milvus** performs fast cosine-similarity search using IVF_FLAT indexes, returning the most relevant evidence from each collection.
4. **The knowledge graph** adds structured facts (target antigen profiles, toxicity grading scales, manufacturing parameters, biomarker cutoffs, regulatory timelines) that complement the vector search results.
5. **Claude** receives the evidence and knowledge context in a carefully constructed prompt, and generates a grounded, citation-rich response.

### The 11 collections as "specialized libraries"

Think of each collection as a specialized library shelf:

- **Literature** is the research library -- published papers and patent filings
- **Trials** is the clinical registry -- active and completed clinical trials
- **Constructs** is the engineering workshop -- the exact designs of the 6 FDA-approved products
- **Assays** is the laboratory notebook -- experimental results from key studies
- **Manufacturing** is the factory floor -- process parameters and specifications
- **Safety** is the pharmacovigilance office -- adverse event reports and safety signals
- **Biomarkers** is the diagnostic lab -- predictive and monitoring markers
- **Regulatory** is the FDA filing cabinet -- approval dates, designations, and label changes
- **Sequences** is the structural biology archive -- scFv sequences and binding affinities
- **Real-World** is the outcomes research center -- registry data and post-market outcomes
- **Genomic Evidence** is the genome center -- 3.5 million patient variant records

When you ask a question, the system does not just search one library -- it searches all eleven simultaneously and then cross-references the findings.

---

## Chapter 5: Your First Query

This chapter walks you through the experience of using the system for the first time.

### Opening the UI

Once the system is running (see Chapter 9 for setup), open your browser and navigate to:

```
http://localhost:8521
```

You will see the CAR-T Intelligence Agent interface with:

- A dark theme (NVIDIA black + green)
- A chat input at the bottom of the page
- A sidebar with configuration options, collection statistics, and demo queries

### Asking a question

Type your question in the chat input. For your first query, try one of these:

```
Why do CD19 CAR-T therapies fail in relapsed B-ALL?
```

Press Enter. Here is what happens:

1. **Search status** appears, showing "Searching across CAR-T data sources..."
2. The system reports how many results it found and from which collections.
3. An **Evidence Sources** expander appears below the status, showing the raw evidence cards.
4. The **LLM response** streams in token by token, with markdown formatting and clickable citations.

### Understanding the response

The response will contain several types of content:

**Cited evidence:** Clickable links like `[Literature:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/)` that open the original source in PubMed. Trial citations look like `[Trial:NCT12345678](https://clinicaltrials.gov/study/NCT12345678)`. These are real links to real papers and trials.

**Cross-functional insights:** The response connects evidence from different domains. For example, it might explain how a resistance mechanism (from literature) relates to a specific assay result (from the assays collection) and a clinical trial outcome (from the trials collection).

**Structured analysis:** For complex topics, the response often includes categorized sections (e.g., "Antigen Escape Mechanisms," "Manufacturing Factors," "Patient-Intrinsic Factors").

### Understanding the evidence panel

Click the "Evidence Sources" expander to see the raw evidence. Each evidence card shows:

- **Collection badge** (color-coded): Literature (blue), Trial (purple), Construct (green), Assay (yellow), Manufacturing (orange), Safety (red), Biomarker (teal), Regulatory (indigo), Sequence (pink), RealWorld (brown), Genomic (cyan)
- **Record ID**: The PMID, NCT number, or internal record identifier
- **Similarity score**: How closely the evidence matches your question (0.0-1.0)
- **Relevance tag**: `[high]`, `[medium]`, or `[low]` based on the score
- **Source link**: For literature and trials, a clickable link to PubMed or ClinicalTrials.gov
- **Text snippet**: The first 200 characters of the evidence text

### Reading citation scores

| Score Range | Relevance | What it means |
|-------------|-----------|--------------|
| 0.75 - 1.00 | **High** | Strong semantic match to your question. This evidence is directly relevant. |
| 0.60 - 0.74 | **Medium** | Partial match. The evidence is related but may address a subtopic or adjacent concept. |
| 0.40 - 0.59 | **Low** | Weak match. The evidence has some thematic overlap but may not directly answer your question. |
| Below 0.40 | Filtered out | Not returned. Below the minimum score threshold. |

### Using sidebar controls

The sidebar gives you fine-grained control over searches:

- **Target Antigen Filter**: Restrict results to a specific target (e.g., CD19, BCMA). This adds a Milvus field filter: `target_antigen == "CD19"`.
- **Development Stage**: Focus on a specific phase (Target Identification, CAR Design, Vector Engineering, Testing, Clinical).
- **Date Range**: Filter evidence by publication or trial start year.
- **Collection toggles**: Enable or disable specific collections. Each toggle shows the live record count.
- **Deep Research Mode**: Enables the autonomous agent pipeline, which decomposes complex questions into sub-queries and evaluates evidence quality before generating the response.

### Downloading results

After each response, three download buttons appear:

- **Markdown**: A formatted `.md` file with the query, response, evidence, and metadata
- **JSON**: A structured `.json` file suitable for programmatic processing
- **PDF**: A formatted PDF report (generated with ReportLab)

---

## Chapter 6: Understanding Collections

### The 11 collections at a glance

| # | Collection | Records | Source | Updated |
|---|-----------|---------|--------|---------|
| 1 | `cart_literature` | 5,047 | PubMed abstracts + patent filings | Weekly ingest via NCBI E-utilities |
| 2 | `cart_trials` | 973 | ClinicalTrials.gov API v2 | Weekly ingest |
| 3 | `cart_constructs` | 41 | CAR-T product designs and engineering approaches | Manual curation |
| 4 | `cart_assays` | 75 | Landmark publications (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1) | Manual curation |
| 5 | `cart_manufacturing` | 56 | Published CMC/process data | Manual curation |
| 6 | `cart_safety` | 71 | FAERS, trial safety data, label safety sections | Manual curation + FAERS ingest |
| 7 | `cart_biomarkers` | 60 | Published biomarker studies | Manual curation |
| 8 | `cart_regulatory` | 40 | FDA/EMA approval records, designations, label updates | Manual curation |
| 9 | `cart_sequences` | 40 | scFv/CAR sequence and structural data | Manual curation |
| 10 | `cart_realworld` | 54 | CIBMTR registry, institutional series | Manual curation |
| 11 | `genomic_evidence` | 3,561,170 | VCF variant data from the rag-chat-pipeline | Shared (read-only) |

**Total: 3,567,622 vectors**

### Collection details

#### 1. cart_literature (5,047 records)

**What it contains:** Published research papers from PubMed and patent filings related to CAR-T therapy. Each record includes the title, abstract text, publication year, target antigen, disease indication, CAR-T development stage, and journal name.

**Why it matters:** Published literature is the foundation of scientific knowledge. When you ask "What are the resistance mechanisms to BCMA-targeted CAR-T?", the literature collection contains the peer-reviewed evidence.

**Example questions it helps answer:**
- "What mechanisms drive antigen escape after CD19 CAR-T?"
- "How does tonic signaling lead to T-cell exhaustion?"
- "What is the evidence for dual-targeting CARs?"

**Key fields:** `id` (PMID), `title`, `text_chunk`, `source_type`, `year`, `cart_stage`, `target_antigen`, `disease`, `keywords`, `journal`

#### 2. cart_trials (973 records)

**What it contains:** Clinical trial records from ClinicalTrials.gov, including the trial title, summary, phase, recruitment status, sponsor, target antigen, CAR generation, costimulatory domain, disease indication, enrollment size, start year, and outcome summary.

**Why it matters:** Clinical trials represent the bridge between laboratory research and patient treatment. This collection lets you find active trials for specific targets, compare sponsor strategies, and identify trends in trial design.

**Example questions it helps answer:**
- "How many CAR-T trials are targeting CD22 in Phase 2?"
- "What is the enrollment status of BCMA trials by Janssen?"
- "Which companies are developing allogeneic CAR-T products?"

**Key fields:** `id` (NCT number), `title`, `text_summary`, `phase`, `status`, `sponsor`, `target_antigen`, `car_generation`, `costimulatory`, `disease`, `enrollment`, `start_year`

#### 3. cart_constructs (41 records)

**What it contains:** Detailed design specifications for the six FDA-approved CAR-T products. Each record documents the target antigen, scFv origin (antibody clone), costimulatory domain, signaling domain, hinge/transmembrane region, vector type, FDA status, and known toxicities.

**Why it matters:** Understanding the exact design of approved products is essential for designing new ones. This collection lets you compare how Kymriah (4-1BB costimulation, lentiviral vector) differs from Yescarta (CD28 costimulation, retroviral vector) and what those design choices mean for clinical outcomes.

**Example questions it helps answer:**
- "What scFv does Carvykti use?"
- "Compare the hinge region of Kymriah vs Yescarta"
- "Which approved products use lentiviral vs retroviral vectors?"

**Key fields:** `id`, `name`, `text_summary`, `target_antigen`, `scfv_origin`, `costimulatory_domain`, `signaling_domain`, `generation`, `hinge_tm`, `vector_type`, `fda_status`, `known_toxicities`

#### 4. cart_assays (75 records)

**What it contains:** Laboratory test results from landmark CAR-T publications. Records cover cytotoxicity assays (how well the CAR-T cells kill target cells), cytokine release measurements, proliferation data, persistence studies, and exhaustion marker analysis.

**Why it matters:** Preclinical assay data determines whether a CAR-T construct advances to clinical trials. By collecting standardized assay results, the system can help researchers benchmark their own results against published data.

**Example questions it helps answer:**
- "What is the typical E:T ratio for CD19 cytotoxicity assays?"
- "How does IFN-gamma secretion compare between 4-1BB and CD28 constructs?"
- "What cell lines are used to test BCMA-directed CARs?"

**Key fields:** `id`, `text_summary`, `assay_type`, `construct_id`, `target_antigen`, `cell_line`, `effector_ratio`, `key_metric`, `metric_value`, `outcome`

#### 5. cart_manufacturing (56 records)

**What it contains:** Manufacturing process parameters across the entire CAR-T production workflow: transduction efficiency, expansion protocols, harvest specifications, cryopreservation conditions, release testing criteria, and logistics data.

**Why it matters:** Manufacturing is the most common bottleneck in CAR-T therapy. Understanding typical process parameters, failure modes, and acceptance criteria helps researchers optimize their own manufacturing workflows.

**Example questions it helps answer:**
- "What transduction efficiency is typical for lentiviral vectors?"
- "What are the release testing requirements for CAR expression?"
- "How does point-of-care manufacturing reduce vein-to-vein time?"

**Key fields:** `id`, `text_summary`, `process_step`, `vector_type`, `parameter`, `parameter_value`, `target_spec`, `met_spec`, `batch_id`

#### 6. cart_safety (71 records)

**What it contains:** Adverse event reports and safety data for CAR-T products. Records include the product name, event type (CRS, ICANS, cytopenia, infection, secondary malignancy), severity grade, onset timing, incidence rate, management protocol, outcome, and reporting source.

**Why it matters:** Safety monitoring is critical for both clinical development and post-market surveillance. This collection enables questions about comparative safety profiles, temporal patterns of adverse events, and management protocols.

**Example questions it helps answer:**
- "What is the incidence of grade 3+ CRS with Carvykti?"
- "How does ICANS management differ between CD19 and BCMA products?"
- "What is the timeline for secondary malignancy reports?"

**Key fields:** `id`, `text_summary`, `product`, `event_type`, `severity_grade`, `onset_timing`, `incidence_rate`, `management_protocol`, `outcome`, `reporting_source`, `year`

#### 7. cart_biomarkers (60 records)

**What it contains:** Predictive, prognostic, pharmacodynamic, monitoring, and resistance biomarkers used in CAR-T therapy. Each record includes the biomarker name, type, assay method, clinical cutoff value, predictive value, associated outcome, and evidence level.

**Why it matters:** Biomarkers help predict which patients will respond, which will develop severe toxicity, and when relapse is imminent. This collection consolidates scattered biomarker evidence into a searchable resource.

**Example questions it helps answer:**
- "Which biomarkers predict CRS severity?"
- "What is the clinical cutoff for ferritin in CRS prediction?"
- "How does MRD negativity at day 28 correlate with PFS?"

**Key fields:** `id`, `text_summary`, `biomarker_name`, `biomarker_type`, `assay_method`, `clinical_cutoff`, `predictive_value`, `associated_outcome`, `evidence_level`

#### 8. cart_regulatory (40 records)

**What it contains:** FDA and EMA regulatory milestones for CAR-T products. Records cover BLA filings, breakthrough therapy designations, RMAT designations, initial approvals, supplemental approvals, label updates, REMS requirements, and post-marketing requirements.

**Why it matters:** Understanding the regulatory landscape is essential for clinical development strategy. This collection enables questions about approval timelines, regulatory designations, and precedent-setting decisions.

**Example questions it helps answer:**
- "When was Breyanzi first approved and for what indication?"
- "Which CAR-T products received RMAT designation?"
- "What post-marketing requirements were imposed on Abecma?"

**Key fields:** `id`, `text_summary`, `product`, `regulatory_event`, `date`, `agency`, `indication`, `decision`, `conditions`, `pivotal_trial`

#### 9. cart_sequences (40 records)

**What it contains:** Molecular and structural data for CAR construct binding domains. Records include scFv clone names, target antigens, binding affinity (Kd), variable heavy and light chain information, framework type, species origin (murine, humanized, fully human), immunogenicity risk, and structural notes.

**Why it matters:** The binding domain is the business end of the CAR. Understanding scFv sequences, binding characteristics, and immunogenicity risks is critical for designing new constructs and improving existing ones.

**Example questions it helps answer:**
- "What is the binding affinity of FMC63?"
- "How do murine scFvs compare to humanized versions for immunogenicity?"
- "What nanobody-based CARs are in development?"

**Key fields:** `id`, `text_summary`, `construct_name`, `target_antigen`, `scfv_clone`, `binding_affinity_kd`, `species_origin`, `immunogenicity_risk`, `structural_notes`

#### 10. cart_realworld (54 records)

**What it contains:** Real-world evidence and outcomes data from post-market registries, institutional series, and claims databases. Records include the study type, data source, product, indication, population size, follow-up duration, primary endpoint, outcome value, care setting (academic vs community), and special populations studied.

**Why it matters:** Clinical trial data is collected under controlled conditions. Real-world evidence shows how therapies perform in routine practice, including in patients who might not have qualified for the original trials (elderly patients, patients with comorbidities, community treatment centers).

**Example questions it helps answer:**
- "How do real-world outcomes for Yescarta compare between academic and community centers?"
- "What is the real-world OS for CAR-T in elderly patients (age 65+)?"
- "How does bridging therapy affect outcomes in routine practice?"

**Key fields:** `id`, `text_summary`, `study_type`, `data_source`, `product`, `indication`, `population_size`, `median_followup_months`, `primary_endpoint`, `outcome_value`, `setting`, `special_population`

#### 11. genomic_evidence (3,561,170 records)

**What it contains:** Genomic variant data from patient whole-genome sequencing, shared from the HCLS AI Factory's rag-chat-pipeline. Each record represents a single variant with chromosome, position, reference/alternate alleles, quality score, gene, consequence, impact level, genotype, clinical significance (from ClinVar), and AlphaMissense pathogenicity prediction.

**Why it matters:** Patient genomic data can reveal why certain patients respond differently to CAR-T therapy. For example, variants in immune-related genes might affect T-cell function, and variants in target-antigen genes might affect antigen expression.

**Example questions it helps answer:**
- "What genomic variants in CD19 pathway genes might affect CAR-T response?"
- "Are there pathogenic variants in immune checkpoint genes?"
- "What is the AlphaMissense classification for variants in the BCMA gene?"

**Key fields:** `id`, `chrom`, `pos`, `ref`, `alt`, `qual`, `gene`, `consequence`, `impact`, `genotype`, `text_summary`, `clinical_significance`, `rsid`, `disease_associations`, `am_pathogenicity`, `am_class`

**Note:** This collection is read-only from the CAR-T Intelligence Agent's perspective. It is created and populated by the rag-chat-pipeline (Stage 2 of the HCLS AI Factory).

### Mapping questions to collections

| If you are asking about... | Primary collections | Supporting collections |
|---------------------------|--------------------|-----------------------|
| Mechanism of action / biology | Literature | Assays, Constructs |
| Clinical outcomes / trial design | Trials | Literature, RealWorld |
| Product design / engineering | Constructs | Sequences, Literature |
| Laboratory results / preclinical data | Assays | Literature, Constructs |
| Manufacturing / CMC | Manufacturing | Literature |
| Adverse events / safety signals | Safety | Biomarkers, Literature, Regulatory |
| Predictive markers / monitoring | Biomarkers | Literature, Trials |
| FDA approvals / regulatory | Regulatory | Literature, Safety |
| Binding domains / molecular design | Sequences | Constructs, Literature |
| Post-market outcomes | RealWorld | Safety, Trials |
| Patient genomics | Genomic Evidence | Literature |

The system searches all collections simultaneously, so you do not need to choose -- but understanding which collections are most relevant helps you interpret the results.

---

## Chapter 7: The Knowledge Graph

### What it is

The knowledge graph is a hand-curated database of structured facts about CAR-T cell therapy. Unlike the vector database (which stores free-text evidence as embeddings), the knowledge graph stores explicit, structured relationships.

Think of the difference this way:

- **Vector database:** "There is a document that mentions CD19, B-ALL, CRS, and Kymriah. It is probably relevant to your question."
- **Knowledge graph:** "CD19 is a B-Lymphocyte Antigen (UniProt P15391). It is expressed on B-cell lineage cells. Four approved products target it. Known resistance mechanisms include CD19 loss, lineage switch, and trogocytosis. CRS incidence is 30-90%."

The knowledge graph provides precise, factual context that helps the LLM generate more accurate and specific answers.

### The three knowledge domains

#### 1. Targets (34 antigens)

For each of 34 target antigens (CD19, BCMA, CD22, CD20, CD30, CD33, CD38, CD123, GD2, HER2, GPC3, EGFR, EGFRvIII, Mesothelin, Claudin18.2, MUC1, PSMA, ROR1, GPRC5D, IL13Ra2, DLL3, B7-H3, NKG2D ligands, CD7, CD5, and 8 additional targets), the knowledge graph stores:

- Full protein name and UniProt identifier
- Expression pattern (where the protein is found)
- Disease indications (which cancers it is relevant for)
- Approved products (if any)
- Key clinical trials
- Known resistance mechanisms
- Toxicity profile (CRS rate, ICANS rate, on-target/off-tumor effects)
- Normal tissue expression (safety implications)

#### 2. Toxicities (17 profiles)

For CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, tumor lysis syndrome, GvHD, and on-target/off-tumor toxicity, the knowledge graph stores:

- Full name and mechanism of action
- Grading system (ASTCT 2019 consensus for CRS/ICANS)
- Incidence rates
- Typical onset timing
- Management protocols (first-line and second-line treatments)
- Relevant biomarkers
- Risk factors

#### 3. Manufacturing (20 processes)

For lentiviral transduction, retroviral transduction, T-cell activation, ex vivo expansion, leukapheresis, cryopreservation, release testing, point-of-care manufacturing, lymphodepletion, and vein-to-vein time, the knowledge graph stores:

- Description and typical parameters
- Critical process parameters
- Failure modes
- Release criteria (where applicable)
- Products that use each approach

#### 4. Biomarkers (23 markers)

For ferritin, CRP, IL-6, sIL-2R, CAR-T expansion (Cmax), Tcm%, CD4:CD8 ratio, LDH, PD-1, LAG-3, TIM-3, MRD (flow), ctDNA, sBCMA, and IFN-gamma, the knowledge graph stores:

- Full name and biomarker type (predictive, prognostic, pharmacodynamic, monitoring, resistance)
- Assay method
- Clinical cutoff values
- Predictive value (with effect sizes where available)
- Associated clinical outcome
- Evidence level (validated, emerging, exploratory)
- Key references (PMIDs)

#### 5. Regulatory (6 products)

For each FDA-approved product, the knowledge graph stores:

- Generic name and manufacturer
- Initial FDA approval date and indication
- Pivotal clinical trial
- Regulatory designations (Breakthrough Therapy, RMAT, Priority Review, Orphan Drug)
- Subsequent supplemental approvals (with dates, indications, and supporting trials)
- REMS requirements
- Post-marketing commitments
- EMA approval date

#### 6. Immunogenicity (6 topics)

For murine scFv immunogenicity, humanization strategies, ADA clinical impact, HLA-restricted epitopes, immunogenicity testing paradigms, and allogeneic HLA considerations, the knowledge graph stores:

- Topic description
- Key constructs and risk levels
- ADA incidence rates
- Clinical impact data
- Management and testing approaches
- Computational tools and methods

### How the knowledge graph differs from the vector database

| Feature | Vector Database (Milvus) | Knowledge Graph (Python dicts) |
|---------|-------------------------|-------------------------------|
| Data type | Free-text embeddings | Structured key-value facts |
| Size | 3.5 million vectors | ~70 entities |
| Search method | Cosine similarity | Keyword matching |
| Strength | Finding relevant text you did not know existed | Providing precise facts about known entities |
| Weakness | May return tangentially related text | Only covers pre-curated entities |
| Role in RAG | Retrieval (finding evidence) | Augmentation (adding structured context) |

Both are essential. The vector database finds evidence. The knowledge graph ensures the LLM has the right factual context.

---

## Chapter 8: Query Expansion

### The problem

When you search for "CRS," you want results that also mention:

- "cytokine release syndrome" (the full name)
- "cytokine storm" (a related term)
- "tocilizumab" (the primary treatment)
- "IL-6" (the key cytokine)
- "ferritin" (a predictive biomarker)
- "ASTCT grading" (the grading system)

A pure vector similarity search will catch some of these (because the embeddings capture semantic relationships), but it may miss others -- especially when the related term uses completely different words.

Query expansion solves this by explicitly mapping keywords to their related terms.

### The 12 expansion maps

The system contains 12 hand-curated expansion dictionaries, each covering a specific domain:

| # | Map | Keywords | Total Terms | Domain |
|---|-----|----------|-------------|--------|
| 1 | Target Antigen | 26 | ~260 | CD19, BCMA, CD22, etc. |
| 2 | Disease | 16 | ~200 | B-ALL, DLBCL, myeloma, etc. |
| 3 | Toxicity | 12 | ~165 | CRS, ICANS, HLH, cytopenias, etc. |
| 4 | Manufacturing | 16 | ~200 | transduction, expansion, cryo, etc. |
| 5 | Mechanism | 14 | ~190 | resistance, exhaustion, persistence, etc. |
| 6 | Construct | 16 | ~200 | scFv, hinge, armored CAR, bispecific, etc. |
| 7 | Safety | 8 | ~80 | adverse events, REMS, secondary malignancy, etc. |
| 8 | Biomarker | 12 | ~100 | ferritin, MRD, exhaustion markers, etc. |
| 9 | Regulatory | 8 | ~60 | FDA, BLA, RMAT, label updates, etc. |
| 10 | Sequence | 8 | ~65 | scFv, CDR, binding affinity, nanobody, etc. |
| 11 | Real-World | 10 | ~75 | CIBMTR, registry, community, disparities, etc. |
| 12 | Immunogenicity | 11 | ~100 | ADA, humanization, HLA, HAMA, etc. |

### How expansion works

When you ask a question, the `expand_query()` function:

1. Converts the question to lowercase
2. Scans it for any keyword that appears in any of the 12 maps
3. Collects all related terms for every matched keyword
4. Deduplicates and returns the expansion terms

For example, given the query "Why do patients get CRS after CD19 CAR-T?":

- **"crs"** matches the Toxicity map, expanding to: "CRS", "cytokine release syndrome", "cytokine storm", "tocilizumab", "Actemra", "siltuximab", "IL-6", "ferritin", "CRP", "Lee grading", "ASTCT grading", "vasopressors", "dexamethasone", and more.
- **"cd19"** matches the Target Antigen map, expanding to: "CD19", "B-ALL", "DLBCL", "tisagenlecleucel", "axicabtagene ciloleucel", "Kymriah", "Yescarta", "Breyanzi", "Tecartus", "FMC63", "SJ25C1", and more.

### How expanded terms are used

Expanded terms are used in two ways:

1. **Known target antigens** (like CD19, BCMA) are used as **field filters** on collections that have a `target_antigen` field. This is a precise filter: `target_antigen == "CD19"`.

2. **Non-antigen terms** (like "cytokine release syndrome", "tocilizumab") are **re-embedded** as separate queries and searched across all collections. This is a semantic search that finds additional evidence the original query might have missed.

The results from expanded searches are merged with the original results, deduplicated, and scored at a lower weight (0.7-0.8x) to prioritize the original direct-match results.

---

## Chapter 9: Setting Up Locally

### Prerequisites

You will need:

- **Python 3.10+** (check with `python3 --version`)
- **Milvus 2.4** running on `localhost:19530` (see below for Docker setup)
- **ANTHROPIC_API_KEY** environment variable (required for LLM synthesis; retrieval works without it)
- **Approximately 4 GB of disk space** for the vector database
- **Internet access** for the initial PubMed and ClinicalTrials.gov data ingestion

### Step 1: Clone the repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd hcls-ai-factory/ai_agent_adds/cart_intelligence_agent
```

### Step 2: Install Python dependencies

```bash
pip install -r requirements.txt
```

This installs the core stack:
- `pymilvus` (vector database client)
- `sentence-transformers` (BGE-small-en-v1.5 embedding model)
- `anthropic` (Claude API client)
- `streamlit` (UI framework)
- `fastapi` + `uvicorn` (REST API)
- `pydantic` + `pydantic-settings` (data models and configuration)
- `biopython` + `lxml` (PubMed data parsing)
- `reportlab` (PDF export)
- `pyvis` (knowledge graph visualization)
- `loguru` (structured logging)
- `apscheduler` (automated data refresh)

### Step 3: Start Milvus

If you have Docker installed, the simplest approach is:

```bash
# Start Milvus standalone (creates a local vector database on port 19530)
docker compose up -d milvus-standalone
```

Or if you are using the HCLS AI Factory's full docker-compose:

```bash
cd /home/adam/projects/hcls-ai-factory
docker compose up -d milvus-standalone milvus-etcd milvus-minio
```

Verify Milvus is running:

```bash
curl http://localhost:19530/v1/vector/collections
```

### Step 4: Create collections and seed FDA construct data

```bash
python3 scripts/setup_collections.py --seed-constructs
```

This creates 11 Milvus collections (10 CAR-T-specific + 1 shared genomic) with IVF_FLAT indexes and inserts the 6 FDA-approved CAR-T product records (Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti) into the `cart_constructs` collection.

### Step 5: Ingest PubMed literature (~15 minutes)

```bash
python3 scripts/ingest_pubmed.py --max-results 5000
```

This fetches CAR-T-related abstracts from PubMed via NCBI E-utilities (esearch + efetch), classifies each paper by CAR-T development stage, extracts target antigens, embeds the text with BGE-small-en-v1.5, and stores the results in `cart_literature`.

### Step 6: Ingest clinical trials (~3 minutes)

```bash
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

This fetches CAR-T clinical trials from the ClinicalTrials.gov API v2, extracts phase, status, sponsor, target antigen, CAR generation, and other metadata, embeds the trial summaries, and stores them in `cart_trials`.

### Step 7: Seed domain-specific data

```bash
# Seed 75 curated assay records from landmark papers
python3 scripts/seed_assays.py

# Seed 56 curated manufacturing/CMC records
python3 scripts/seed_manufacturing.py

# Seed safety records
python3 scripts/seed_safety.py

# Seed biomarker records
python3 scripts/seed_biomarkers.py

# Seed regulatory milestone records
python3 scripts/seed_regulatory.py

# Seed scFv/CAR sequence data
python3 scripts/seed_sequences.py

# Seed real-world evidence records
python3 scripts/seed_realworld.py

# Seed knowledge graph to JSON (optional, for inspection)
python3 scripts/seed_knowledge.py
```

### Step 8: Validate the setup

```bash
python3 scripts/validate_e2e.py
```

This runs end-to-end validation: collection statistics, single-collection search, multi-collection `search_all()`, filtered search, and all demo queries.

### Step 9: (Optional) Test the full RAG pipeline

```bash
# Requires ANTHROPIC_API_KEY to be set
python3 scripts/test_rag_pipeline.py
```

This tests the complete pipeline: embed query, search all collections, augment with knowledge graph, and generate a Claude LLM response.

### Step 10: Launch the Streamlit UI

```bash
streamlit run app/cart_ui.py --server.port 8521
```

Open `http://localhost:8521` in your browser. You should see the CAR-T Intelligence Agent interface with live collection statistics in the sidebar.

### Step 11: (Optional) Launch the REST API

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8522 --reload
```

The API documentation is automatically available at `http://localhost:8522/docs` (Swagger UI).

---

## Chapter 10: Exploring the API

### What the REST API is

The REST API wraps the same RAG engine that powers the Streamlit UI, but exposes it as a set of HTTP endpoints. This enables:

- Integration with other applications (dashboards, notebooks, pipelines)
- Programmatic access from any programming language
- Automated testing and monitoring

### The 13 endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with collection count and total vector count |
| GET | `/collections` | All collection names and their record counts |
| POST | `/query` | Full RAG query (retrieve + LLM synthesis) |
| POST | `/search` | Evidence-only retrieval (no LLM, fast) |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics |

### Testing with curl

**Health check:**

```bash
curl http://localhost:8522/health
```

Expected response:

```json
{
  "status": "healthy",
  "collections": 11,
  "total_vectors": 3567622
}
```

**List collections:**

```bash
curl http://localhost:8522/collections
```

Expected response:

```json
{
  "collections": [
    {"name": "cart_literature", "record_count": 5047},
    {"name": "cart_trials", "record_count": 973},
    {"name": "cart_constructs", "record_count": 41},
    ...
  ],
  "total": 11
}
```

**Evidence-only search (fast, no LLM):**

```bash
curl -X POST http://localhost:8522/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "CD19 antigen escape mechanisms",
    "target_antigen": "CD19"
  }'
```

**Full RAG query (with LLM synthesis):**

```bash
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Why do CD19 CAR-T therapies fail in relapsed B-ALL?",
    "target_antigen": "CD19",
    "year_min": 2020
  }'
```

**Cross-collection entity search:**

```bash
curl -X POST http://localhost:8522/find-related \
  -H "Content-Type: application/json" \
  -d '{
    "entity": "Yescarta",
    "top_k": 3
  }'
```

This returns everything the system knows about Yescarta across all 11 collections -- literature, trials, construct design, safety data, regulatory milestones, and more.

**Knowledge graph statistics:**

```bash
curl http://localhost:8522/knowledge/stats
```

Expected response:

```json
{
  "target_antigens": 34,
  "targets_with_approved_products": 2,
  "toxicity_profiles": 17,
  "manufacturing_processes": 20,
  "biomarkers": 23,
  "regulatory_products": 6
}
```

### API documentation

The FastAPI framework automatically generates interactive API documentation. Once the API is running, open:

- **Swagger UI:** `http://localhost:8522/docs`
- **ReDoc:** `http://localhost:8522/redoc`

These pages let you explore endpoints, view request/response schemas, and test API calls directly from the browser.

---

## Chapter 11: Understanding the Codebase

### Project structure

```
cart_intelligence_agent/
├── src/                          # Core source code (the engine)
│   ├── __init__.py
│   ├── models.py                 # Pydantic data models (16 models + enums)
│   ├── collections.py            # 11 Milvus collection schemas + manager
│   ├── knowledge.py              # Knowledge graph (3 domains)
│   ├── query_expansion.py        # 12 expansion maps
│   ├── rag_engine.py             # Multi-collection RAG engine + comparative
│   ├── agent.py                  # Autonomous agent (plan -> search -> synthesize)
│   ├── export.py                 # Markdown, JSON, PDF export
│   ├── metrics.py                # Prometheus metrics
│   ├── scheduler.py              # Automated data refresh
│   ├── ingest/                   # Data ingestion parsers
│   │   ├── base.py               # Base ingest pipeline
│   │   ├── literature_parser.py  # PubMed NCBI E-utilities
│   │   ├── clinical_trials_parser.py  # ClinicalTrials.gov API v2
│   │   ├── construct_parser.py   # CAR construct data
│   │   ├── assay_parser.py       # Assay data
│   │   ├── manufacturing_parser.py # Manufacturing/CMC data
│   │   ├── safety_parser.py      # Pharmacovigilance data
│   │   ├── biomarker_parser.py   # Biomarker data
│   │   ├── regulatory_parser.py  # FDA regulatory data
│   │   ├── sequence_parser.py    # Molecular/structural data
│   │   ├── realworld_parser.py   # Real-world evidence
│   │   ├── faers_parser.py       # FAERS adverse event reports
│   │   ├── dailymed_parser.py    # DailyMed label data
│   │   ├── uniprot_parser.py     # UniProt protein data
│   │   └── cibmtr_parser.py      # CIBMTR registry data
│   └── utils/
│       └── pubmed_client.py      # NCBI E-utilities HTTP client
├── app/
│   └── cart_ui.py                # Streamlit chat interface (v2.0)
├── api/
│   ├── __init__.py
│   ├── main.py                   # FastAPI REST API (9 endpoints)
│   └── routes/                   # Additional route modules (4 endpoints)
├── config/
│   └── settings.py               # Pydantic BaseSettings configuration
├── data/
│   └── reference/
│       ├── assay_seed_data.json
│       └── manufacturing_seed_data.json
├── scripts/                      # CLI scripts for setup and data ingestion
│   ├── setup_collections.py
│   ├── ingest_pubmed.py
│   ├── ingest_clinical_trials.py
│   ├── seed_assays.py
│   ├── seed_manufacturing.py
│   ├── seed_safety.py
│   ├── seed_biomarkers.py
│   ├── seed_regulatory.py
│   ├── seed_sequences.py
│   ├── seed_realworld.py
│   ├── seed_knowledge.py
│   ├── seed_patents.py
│   ├── seed_immunogenicity.py
│   ├── validate_e2e.py
│   └── test_rag_pipeline.py
├── tests/                        # 415 tests
│   ├── conftest.py
│   ├── test_models.py
│   ├── test_rag_engine.py
│   ├── test_agent.py
│   ├── test_knowledge.py
│   ├── test_query_expansion.py
│   └── test_export.py
├── docs/                         # Documentation
├── Docs/
│   └── CART_Intelligence_Agent_Design.md
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
├── .streamlit/
├── LICENSE                       # Apache 2.0
└── README.md
```

### File-by-file walkthrough of the core files

#### config/settings.py (113 lines)

The single source of truth for all configuration. Uses Pydantic `BaseSettings` so every value can be overridden via environment variables (prefixed with `CART_`) or a `.env` file.

Key settings:
- `MILVUS_HOST` / `MILVUS_PORT`: Vector database connection (default: localhost:19530)
- `EMBEDDING_MODEL`: BGE-small-en-v1.5 (384 dimensions)
- `LLM_MODEL`: Claude Sonnet 4.6
- `TOP_K_PER_COLLECTION`: Maximum results per collection per query (default: 5)
- `SCORE_THRESHOLD`: Minimum cosine similarity score (default: 0.4)
- `WEIGHT_*`: Per-collection search weights (Literature=0.20, Trials=0.16, etc.)
- `CITATION_HIGH_THRESHOLD` / `CITATION_MEDIUM_THRESHOLD`: Score thresholds for relevance tagging (0.75 / 0.60)

#### src/models.py (484 lines)

All Pydantic data models used throughout the system. Contains:

- **13 enums:** `CARTStage`, `SourceType`, `TrialPhase`, `TrialStatus`, `CARGeneration`, `AssayType`, `ProcessStep`, `FDAStatus`, `SafetyEventType`, `BiomarkerType`, `EvidenceLevel`, `RegulatoryEvent`, `RWEStudyType`
- **10 collection models** (one per owned collection): `CARTLiterature`, `ClinicalTrial`, `CARConstruct`, `AssayResult`, `ManufacturingRecord`, `SafetyRecord`, `BiomarkerRecord`, `RegulatoryRecord`, `SequenceRecord`, `RealWorldRecord`
- **Search result models:** `SearchHit`, `CrossCollectionResult`, `ComparativeResult`
- **Agent models:** `AgentQuery`, `AgentResponse`

Each collection model has a `to_embedding_text()` method that generates the text string used for BGE-small embedding. This is important -- the quality of the embedding depends on how the text is constructed from the structured fields.

#### src/collections.py (1,004 lines)

Manages the 11 Milvus collections. Contains:

- **Schema definitions** for all 11 collections, defining every field's name, data type, and constraints
- **COLLECTION_SCHEMAS** registry mapping collection names to their schemas
- **CARTCollectionManager** class with methods for:
  - `connect()` / `disconnect()`: Milvus connection lifecycle
  - `create_collection()` / `create_all_collections()`: Schema creation with IVF_FLAT index
  - `get_collection_stats()`: Record counts per collection
  - `insert_batch()`: Bulk data insertion
  - `search()`: Single-collection vector similarity search
  - `search_all()`: Parallel search across ALL collections using `ThreadPoolExecutor`

The `search_all()` method is the workhorse of the system. It launches concurrent searches across all 11 collections and merges the results. This is what enables the system to search 3.5 million vectors in 12-16 milliseconds.

#### src/rag_engine.py (754 lines)

The multi-collection RAG engine. This is the central piece of the system. It orchestrates:

1. **Query embedding** (with BGE asymmetric query prefix)
2. **Parallel collection search** (via the collection manager)
3. **Query expansion** (via the expansion module)
4. **Knowledge graph augmentation** (from all 3 domains: targets, toxicities, manufacturing)
5. **Result merging and ranking** (deduplication, weighted scoring, relevance tagging)
6. **Prompt construction** (formatting evidence, knowledge context, and the question into a structured LLM prompt)
7. **LLM generation** (synchronous and streaming modes)
8. **Comparative analysis** (detecting "vs" queries, parsing entities, running dual retrievals, building comparison prompts)
9. **Cross-collection entity linking** (`find_related()` method)

The system prompt (`CART_SYSTEM_PROMPT`) defines a detailed 12-domain expert persona for Claude. It instructs the LLM to cite evidence using clickable markdown links, think cross-functionally, highlight failure modes, suggest optimizations, and acknowledge uncertainty.

#### src/knowledge.py (2,249 lines)

The knowledge graph. Contains four large Python dictionaries:

- `CART_TARGETS`: 34 target antigen profiles
- `CART_TOXICITIES`: 17 toxicity profiles with grading, management, biomarkers
- `CART_MANUFACTURING`: 20 manufacturing process specifications
- `CART_BIOMARKERS`: 23 biomarker profiles with cutoffs and evidence levels
- `CART_REGULATORY`: 6 FDA-approved product regulatory timelines
- `CART_IMMUNOGENICITY`: 6 immunogenicity topic profiles
- `ENTITY_ALIASES`: Mapping of product names, generic names, and domain terms to canonical entities (used for comparative analysis entity resolution)

Public API functions:
- `get_target_context()`: Formatted text for a target antigen
- `get_toxicity_context()`: Formatted text for a toxicity profile
- `get_manufacturing_context()`: Formatted text for a manufacturing process
- `get_biomarker_context()`: Formatted text for a biomarker
- `get_regulatory_context()`: Formatted text for a product's regulatory history
- `get_immunogenicity_context()`: Formatted text for an immunogenicity topic
- `resolve_comparison_entity()`: Resolves a raw text string to a known entity for comparative analysis
- `get_comparison_context()`: Builds side-by-side knowledge for two entities

#### src/query_expansion.py (1,592 lines)

The 12 domain-specific expansion dictionaries and the expansion functions. The main function is `expand_query()`, which takes a raw user question and returns a deduplicated list of related terms.

Also provides:
- `expand_query_by_category()`: Returns terms grouped by category (useful for weighted per-collection expansion)
- `get_expansion_stats()`: Returns keyword and term counts per map (useful for health checks)

#### src/agent.py (309 lines)

The autonomous CAR-T Intelligence Agent. Implements the plan-search-synthesize pattern:

1. **`search_plan()`**: Analyzes the question to identify target antigens, relevant development stages, and search strategy (broad, targeted, or comparative). Decomposes complex "why did X fail?" questions into sub-queries.
2. **`evaluate_evidence()`**: Assesses evidence quality as "sufficient" (3+ collections, 10+ hits), "partial" (2+ collections, 5+ hits), or "insufficient" (fewer).
3. **`run()`**: Orchestrates the full pipeline: plan, search, evaluate, optionally expand with sub-questions, generate answer.
4. **`generate_report()`**: Creates a formatted markdown report from the results.

#### app/cart_ui.py (1,162 lines)

The Streamlit chat interface. Features:

- **Three tabs:** Chat, Knowledge Graph Explorer, Image Analysis
- **NVIDIA dark theme** with custom CSS (black background, green accents)
- **Sidebar controls:** Deep Research toggle, target antigen filter, development stage filter, date range filter, collection toggles (with live record counts), demo queries
- **Evidence cards** with collection badges, relevance indicators, and clickable source links
- **Comparative analysis mode** with entity-grouped evidence and "VS" divider
- **Streaming LLM response** with real-time token display
- **Knowledge graph visualization** using pyvis (interactive network graph)
- **Image analysis** using Claude Vision for claim extraction and evidence verification
- **Export buttons** for Markdown, JSON, and PDF download
- **Conversation memory** for contextual follow-up queries

#### api/main.py (588 lines)

The FastAPI REST API. Contains:

- **Lifespan management:** Initializes the embedding model, LLM client, collection manager, and RAG engine on startup; disconnects on shutdown.
- **13 endpoints** (detailed in Chapter 10)
- **Pydantic request/response schemas:** `QueryRequest`, `QueryResponse`, `SearchResponse`, `FindRelatedRequest`, `FindRelatedResponse`, etc.
- **CORS middleware:** Restricts cross-origin requests to 3 configured origins (localhost:8080, 8521, 8522)
- **Basic request metrics** exposed at `/metrics` in Prometheus format

### How the pieces connect: a trace through the code

When you type "Compare 4-1BB vs CD28" in the UI:

1. `app/cart_ui.py` receives the input in `st.chat_input()`
2. The UI calls `engine._is_comparative("Compare 4-1BB vs CD28")` which returns `True`
3. The UI calls `engine.retrieve_comparative("Compare 4-1BB vs CD28", ...)`
4. `rag_engine.py` calls `_parse_comparison_entities()` which:
   - Extracts "4-1BB" and "CD28" using regex
   - Calls `knowledge.resolve_comparison_entity("4-1BB")` which returns `{"type": "costimulatory", "canonical": "4-1BB (CD137)"}`
   - Calls `knowledge.resolve_comparison_entity("CD28")` which returns `{"type": "costimulatory", "canonical": "CD28"}`
5. Two separate `retrieve()` calls are made -- one for each entity
6. Each `retrieve()` call:
   - Embeds the query via BGE-small
   - Searches all 11 collections in parallel via `collections.search_all()`
   - Runs query expansion via `query_expansion.expand_query()`
   - Retrieves knowledge graph context via `knowledge.get_target_context()` etc.
7. Results are merged into a `ComparativeResult`
8. The UI calls `engine._build_comparative_prompt()` which creates a structured prompt asking Claude to produce a comparison table, advantages, limitations, and clinical context
9. Claude generates the comparative response, streamed to the UI

---

## Chapter 12: Next Steps

### Where to go from here

Now that you understand the foundations, here are pathways for further learning:

**If you are a biologist / clinician:**
- Try asking the system questions from your own research area
- Explore the knowledge graph tab to see how entities are related
- Use the image analysis feature to verify claims from presentations
- Look at `src/knowledge.py` to see what structured data is available for your targets of interest

**If you are a data scientist / ML engineer:**
- Read `src/rag_engine.py` to understand the retrieval-augmentation-generation pipeline in detail
- Experiment with collection weights in `config/settings.py` to see how they affect result ranking
- Look at the query expansion maps in `src/query_expansion.py` and consider what additional domain maps might help
- Try modifying the system prompt in `rag_engine.py` to adjust the LLM's behavior

**If you are a software developer:**
- Read `src/collections.py` to understand Milvus schema design and index configuration
- Look at `api/main.py` to see how the RAG engine is exposed as a REST API
- Explore the ingest parsers in `src/ingest/` to understand how data flows from PubMed and ClinicalTrials.gov into the vector database
- Run the test suite: `pytest tests/ -v`

### How to contribute

The CAR-T Intelligence Agent is open-source under the Apache 2.0 license. Contributions are welcome in these areas:

- **New data sources**: Additional parsers for FDA FAERS, CIBMTR, SEER, or other databases
- **New collections**: Additional specialized collections for specific data types
- **Knowledge graph expansion**: Adding new target antigens, toxicity profiles, or manufacturing processes
- **Query expansion refinement**: Improving term coverage and precision
- **UI enhancements**: Better visualization, filtering, and export options
- **Performance optimization**: Faster search, better caching, reduced latency
- **Testing**: Additional unit tests, integration tests, and benchmark queries

### Resources for learning more

**CAR-T cell therapy:**
- Maude et al., "Tisagenlecleucel in Children and Young Adults with B-Cell Lymphoblastic Leukemia," NEJM 2018 (ELIANA trial)
- Neelapu et al., "Axicabtagene Ciloleucel CAR T-Cell Therapy in Refractory Large B-Cell Lymphoma," NEJM 2017 (ZUMA-1 trial)
- June and Sadelain, "Chimeric Antigen Receptor Therapy," NEJM 2018 (comprehensive review)
- Lee et al., "ASTCT Consensus Grading for Cytokine Release Syndrome and Neurologic Toxicity," Biology of Blood and Marrow Transplantation 2019

**Vector databases and embeddings:**
- Milvus documentation: https://milvus.io/docs
- BGE embedding models: https://huggingface.co/BAAI/bge-small-en-v1.5
- "Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks" (Lewis et al., 2020) -- the original RAG paper

**This project:**
- HCLS AI Factory GitHub: https://github.com/ajones1923/hcls-ai-factory
- Architecture design: `Docs/CART_Intelligence_Agent_Design.md` (in this repository)
- Full README: `README.md` (in this repository)

---

## Glossary

| Term | Definition |
|------|-----------|
| **4-1BB (CD137)** | A costimulatory receptor. When included as the costimulatory domain in a CAR, it promotes T-cell persistence and memory formation through NF-kB signaling. Used in Kymriah, Breyanzi, Abecma, and Carvykti. |
| **ADA** | Anti-Drug Antibody. An immune response against a therapeutic protein (in this case, the CAR construct). Can reduce CAR-T persistence and impair re-dosing. |
| **Antigen** | A molecule (usually a protein) on the surface of a cell that can be recognized by the immune system. CAR-T therapy targets specific antigens on cancer cells. |
| **Antigen escape** | A resistance mechanism where cancer cells stop expressing the target antigen, making them invisible to CAR-T cells. |
| **Apheresis** | The process of collecting specific blood components (in this case, T-cells) from a patient. Also called leukapheresis. |
| **B-ALL** | B-cell Acute Lymphoblastic Leukemia. A blood cancer that is one of the primary indications for CD19-targeted CAR-T therapy. |
| **BCMA** | B-Cell Maturation Antigen (TNFRSF17). A protein on mature B-cells and plasma cells. The target for Abecma and Carvykti (multiple myeloma treatments). |
| **BGE-small-en-v1.5** | The embedding model used by this system. Produced by BAAI (Beijing Academy of Artificial Intelligence). Converts text to 384-dimensional vectors optimized for semantic similarity search. |
| **BLA** | Biologics License Application. The regulatory submission required for FDA approval of a biological product like a CAR-T therapy. |
| **CAR** | Chimeric Antigen Receptor. An engineered protein that, when expressed on a T-cell, redirects the T-cell to recognize and kill cells displaying a specific antigen. |
| **CAR-T** | Chimeric Antigen Receptor T-cell therapy. A form of immunotherapy where a patient's T-cells are genetically modified to express a CAR, expanded, and infused back into the patient. |
| **CD19** | A protein expressed on B-cell lineage cells. The most common target for CAR-T therapy. Targeted by four of six FDA-approved products. |
| **CD28** | A costimulatory receptor. When included as the costimulatory domain in a CAR, it drives rapid T-cell expansion and strong effector function. Used in Yescarta and Tecartus. |
| **CD3-zeta** | The signaling domain of the T-cell receptor complex. Included in all CARs as the primary activation signal. Contains ITAMs (immunoreceptor tyrosine-based activation motifs) that trigger T-cell killing. |
| **CIBMTR** | Center for International Blood and Marrow Transplant Research. A registry that tracks outcomes of CAR-T therapy and stem cell transplantation in clinical practice. |
| **ClinVar** | A public database of relationships between human genetic variants and clinical conditions, maintained by the NCBI. |
| **Collection** | In Milvus, a collection is a table that stores vectors and their associated metadata. This system has 11 collections. |
| **Cosine similarity** | A measure of similarity between two vectors based on the angle between them. Ranges from 0 (completely different) to 1 (identical direction). Used to find relevant evidence in the vector database. |
| **CRP** | C-Reactive Protein. An inflammatory biomarker. Elevated CRP within 72 hours of CAR-T infusion predicts severe CRS. |
| **CRS** | Cytokine Release Syndrome. An inflammatory response caused by massive cytokine release from activated CAR-T cells. The most common serious side effect of CAR-T therapy. Graded 1-4 by the ASTCT consensus. |
| **ctDNA** | Circulating Tumor DNA. Cell-free DNA released by tumor cells into the bloodstream. Used as a liquid biopsy marker to monitor treatment response. |
| **DLBCL** | Diffuse Large B-Cell Lymphoma. The most common type of non-Hodgkin lymphoma and a major indication for CD19 CAR-T therapy. |
| **DGX Spark** | An NVIDIA desktop computer with a GB10 GPU, 128 GB unified LPDDR5x memory, and 20 ARM cores. The hardware target for this system ($4,699). |
| **Embedding** | A numerical representation of text as a vector (list of numbers) in a high-dimensional space. Similar texts have similar embeddings. |
| **FAERS** | FDA Adverse Event Reporting System. A database of post-market safety reports submitted to the FDA. |
| **FastAPI** | A modern Python web framework for building REST APIs. Used for the CAR-T Intelligence Agent's API server. |
| **FDA** | U.S. Food and Drug Administration. The regulatory agency that approves drugs and biological products in the United States. |
| **FMC63** | A murine monoclonal antibody clone that binds CD19. The scFv from FMC63 is used in the binding domain of Kymriah, Yescarta, Tecartus, and Breyanzi. |
| **GvHD** | Graft-versus-Host Disease. A complication primarily relevant to allogeneic (donor-derived) CAR-T products, where the infused cells attack the patient's normal tissues. |
| **HCLS AI Factory** | Healthcare and Life Sciences AI Factory. The parent platform that includes genomics, RAG, and drug discovery pipelines. |
| **HLH** | Hemophagocytic Lymphohistiocytosis. A severe inflammatory condition that can occur as a complication of CRS, characterized by extremely high ferritin levels and multi-organ dysfunction. |
| **ICANS** | Immune Effector Cell-Associated Neurotoxicity Syndrome. Neurological toxicity associated with CAR-T therapy, caused by blood-brain barrier disruption and CNS inflammation. Graded using the ICE score. |
| **ICE score** | Immune Effector Cell Encephalopathy score. A 10-point assessment tool used to grade ICANS severity. Evaluates orientation, naming, writing, attention, and command following. |
| **IVF_FLAT** | Inverted File with Flat quantization. The Milvus indexing algorithm used by this system. It partitions vectors into clusters (nlist=1024) and searches the nearest clusters (nprobe=16) for fast approximate nearest neighbor search. |
| **Knowledge graph** | A structured database of facts and relationships. In this system, the knowledge graph contains curated data about CAR-T targets, toxicities, manufacturing, biomarkers, regulatory history, and immunogenicity. |
| **Lentiviral vector** | A modified lentivirus (derived from HIV-1) used to deliver the CAR gene into T-cells. Used by Kymriah, Breyanzi, Abecma, and Carvykti. |
| **LLM** | Large Language Model. An AI model trained on vast amounts of text that can generate human-like responses. This system uses Claude Sonnet 4.6 by Anthropic. |
| **Lymphodepletion** | Pre-infusion chemotherapy (typically fludarabine + cyclophosphamide) that depletes the patient's existing immune cells to create space for CAR-T cell expansion. |
| **Milvus** | An open-source vector database designed for similarity search on large-scale embedding datasets. This system uses Milvus 2.4. |
| **MRD** | Minimal Residual Disease. The small number of cancer cells that may remain after treatment. MRD negativity (no detectable cancer cells) after CAR-T therapy is a strong predictor of long-term remission. |
| **Pydantic** | A Python library for data validation using type annotations. Used throughout this system for data models, API schemas, and configuration. |
| **Query expansion** | The process of augmenting a search query with related terms to improve recall. For example, expanding "CRS" to also search for "cytokine release syndrome," "tocilizumab," and "IL-6." |
| **RAG** | Retrieval-Augmented Generation. A technique that combines information retrieval (finding relevant documents) with language model generation (writing an answer based on those documents). |
| **REMS** | Risk Evaluation and Mitigation Strategy. An FDA-required safety program. All CAR-T products have a shared REMS requiring certified treatment centers and trained healthcare providers. |
| **RMAT** | Regenerative Medicine Advanced Therapy. An FDA expedited program designation for cell therapy, gene therapy, and tissue-engineered products. |
| **scFv** | Single-chain Variable Fragment. The antigen-binding portion of a CAR, derived from the variable heavy (VH) and variable light (VL) chains of an antibody, connected by a flexible linker. |
| **Streamlit** | A Python framework for building data applications and dashboards. Used for the CAR-T Intelligence Agent's web UI. |
| **Tocilizumab** | An IL-6 receptor antagonist (brand name: Actemra). The first-line treatment for grade 2+ CRS after CAR-T therapy. |
| **Trogocytosis** | A process where CAR-T cells strip the target antigen from cancer cell surfaces through membrane contact, reducing antigen density and potentially allowing antigen-low cells to escape. |
| **Vector** | In the context of this system, a list of numbers (384 floating-point values) that represents the meaning of a piece of text. |
| **Vector database** | A database optimized for storing and searching high-dimensional vectors. Enables fast similarity search across millions of embeddings. |
| **VCN** | Vector Copy Number. The average number of CAR transgene copies integrated per T-cell. FDA guidelines typically require VCN below 5 copies/cell for safety. |

---

*CAR-T Intelligence Agent | HCLS AI Factory | Apache 2.0 | Adam Jones | March 2026*
