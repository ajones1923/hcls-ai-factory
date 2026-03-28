# HCLS AI Factory -- Learning Guide Foundations

**The Complete Beginner's Guide to Precision Medicine on NVIDIA DGX Spark**

> **Author:** Adam Jones
> **Date:** March 2026
> **License:** Apache 2.0

> This is **Part 1** of a two-part unified learning guide. Part 1 covers the
> platform introduction, architecture, and the three core pipelines. Part 2
> covers the eleven intelligence agents in depth.

---

## Who This Guide Is For

This guide serves three distinct audiences. Every chapter is written so that
all three can follow along. Look for the persona tags when a section speaks
directly to one group.

### Persona 1 -- The Clinician

You care for patients. You may be an oncologist attending a molecular tumor
board, a cardiologist reviewing a pharmacogenomic profile, or a primary care
physician evaluating biomarker results. You want to understand what this
platform does, what evidence it relies on, and how to interpret its output.
You do not need to write code.

### Persona 2 -- The Data Scientist / Bioinformatician

You work with genomic data -- VCF files, variant annotations, gene panels,
embeddings, and machine learning models. You want to understand how the
platform ingests, transforms, and searches data across its three stages and
eleven intelligence agents. You want to know how to extend it with new
collections, annotation sources, or scoring models.

### Persona 3 -- The Software Engineer

You build and deploy applications. You want to understand the architecture --
FastAPI, Streamlit, Milvus, Docker, Nextflow, Pydantic models -- and how to
run, modify, scale, or contribute to the codebase. You care about ports,
containers, health checks, and CI/CD.

---

## What You Will Learn

By the end of this two-part guide you will be able to:

**Part 1 -- Foundations (this document)**

1. Explain what precision medicine is and why it reduces treatment timelines
   from months to hours.
2. Describe the NVIDIA DGX Spark hardware and why GPU acceleration matters
   for genomics, AI, and molecular simulation.
3. Trace a patient's DNA through the three-engine pipeline: Genomic Foundation Engine,
   Precision Intelligence Network, and Therapeutic Discovery Engine.
4. Read and interpret FASTQ, VCF, SMILES, and PDB file formats at a
   conceptual level.
5. Explain how RAG (Retrieval-Augmented Generation) combines a vector
   database with a large language model to produce grounded, cited answers.
6. Describe the shared `genomic_evidence` collection and how all eleven
   agents consume it.
7. Explain how the intelligence agent architecture works: the five-phase
   plan-search-evaluate-synthesize-report loop.

**Part 2 -- Intelligence Agents (separate document)**

8. Describe each of the eleven intelligence agents, their domain-specific
   collections, and their clinical workflows.
9. Compare agent capabilities across oncology, biomarkers, CAR-T, imaging,
   autoimmune, pharmacogenomics, and cardiology.
10. Run a first query against any agent and read the resulting report.

---

## Prerequisites

- **Clinicians:** A web browser. That is all.
- **Data Scientists:** Familiarity with Python, VCF files, and basic
  statistics.
- **Software Engineers:** Python 3.10+, Docker, and a terminal.

No special hardware is needed to read this guide. To run the platform itself,
the reference hardware is an NVIDIA DGX Spark ($4,699), but any machine with
Docker and a CPU will work for testing the intelligence agents.

---

# Chapter 1: The Precision Medicine Revolution

## 1.1 The Old Approach

For most of modern medical history, treatment has been population-based.
Physicians prescribe the same drug at the same dose to every patient with
the same diagnosis. The reasoning is statistical: if a drug works for 30% of
patients with a given condition, it is worth trying -- even though 70% will
not benefit and some will experience side effects for no gain.

In oncology, first-line chemotherapy response rates hover around 30%. In
cardiology, standard heart failure regimens work well for some patients and
poorly for others, with no easy way to predict who will respond. In
autoimmune disease, patients often cycle through three or four medications
over several years before finding one that controls their symptoms.

```
  Traditional Medicine
  ====================

  All patients with      Same drug        30% respond
  Diagnosis X      --->  Same dose  --->  70% don't
                         Same schedule
```

Think of it like prescribing the same pair of eyeglasses to every patient
in an optometrist's office. Some will see better. Some will see worse. Some
will not benefit at all.

## 1.2 The New Approach

Precision medicine flips this model. Instead of asking "what disease does the
patient have?" it asks "what molecular profile does the patient have?" The
answer comes from the patient's DNA.

```
  Precision Medicine
  ==================

  Patient DNA ---> Molecular   ---> Matched       ---> 60-80%
                   Profile          Therapy             response
```

A patient with a BRAF V600E mutation receives a BRAF inhibitor regardless of
whether the tumor is in the skin, colon, or lung -- because the molecular
driver is the same. A patient with a specific CYP2D6 metabolizer status
receives a dose adjustment for their heart medication. A patient with HLA-B
risk alleles avoids a drug that would cause a severe immune reaction.

The key insight: treating the molecule, not just the organ.

## 1.3 The Time Problem

Precision medicine works. The evidence is overwhelming. But there is a
bottleneck: time.

Here is what the traditional precision medicine pipeline looks like:

```
  Traditional Precision Medicine Timeline
  ========================================

  Step 1: DNA Sequencing           1-3 days
  Step 2: Bioinformatics Analysis  1-2 weeks
  Step 3: Variant Interpretation   1-4 weeks (manual literature review)
  Step 4: Clinical Decision        1-2 weeks (tumor board, consult)
  Step 5: Drug Identification      1-6 months (trials, compassionate use)
  -----------------------------------------------
  Total:                           3-12 months
```

For a cancer patient whose tumor is growing, three months is an eternity.
For a heart failure patient decompensating, weeks of delay can mean
hospitalization or death. For a rare disease patient, the diagnostic odyssey
can last years.

The HCLS AI Factory compresses this timeline:

```
  HCLS AI Factory Timeline
  =========================

  Stage 1: Genomic Foundation Engine      2-4 hours (GPU-accelerated)
  Stage 2: Precision Intelligence Network < 5 seconds (RAG query)
  Stage 3: Therapeutic Discovery Engine   8-16 minutes (generative AI)
  -----------------------------------------------
  Total:                           < 5 hours
```

All of this runs on a single $4,699 desktop computer.

**Analogy:** Think of GPS navigation versus paper maps. Both can get you from
New York to Los Angeles. The paper map requires you to plan the route
manually, check for road closures by calling ahead, and recalculate if you
take a wrong turn. GPS does all of that in real time, rerouting around
traffic and construction automatically. The HCLS AI Factory is GPS for
precision medicine -- it does not replace the driver (the clinician), but it
eliminates the hours spent planning the route.

## 1.4 The Three-Stage Pipeline

The platform is organized into three sequential stages, each building on the
output of the previous one. Eleven intelligence agents branch from Stage 2
to provide domain-specific clinical decision support.

```
  +-----------------------------------------------------------------+
  |                    HCLS AI Factory Pipeline                      |
  +-----------------------------------------------------------------+
  |                                                                  |
  |  Patient DNA (FASTQ)                                             |
  |       |                                                          |
  |       v                                                          |
  |  +---------------------------+                                   |
  |  | STAGE 1: Genomic Foundation |                                 |
  |  | BWA-MEM2 + DeepVariant    |                                   |
  |  | Parabricks GPU-accelerated|                                   |
  |  +---------------------------+                                   |
  |       |                                                          |
  |       v                                                          |
  |  VCF (11.7M variants)                                            |
  |       |                                                          |
  |       +---> Annotate (ClinVar, AlphaMissense, VEP)               |
  |       |                                                          |
  |       +---> Embed (BGE-small-en-v1.5, 384 dims)                  |
  |       |                                                          |
  |       v                                                          |
  |  Milvus: genomic_evidence (3.56M vectors)                        |
  |       |                                                          |
  |       +---> [shared read-only access by all 11 agents]           |
  |       |                                                          |
  |       v                                                          |
  |  +-------------------------------+                               |
  |  | STAGE 2: Precision Intel Net |                               |
  |  | RAG + Claude AI + Knowledge   |                               |
  |  +-------------------------------+                               |
  |       |         |       |       |       |       |       |        |
  |       v         v       v       v       v       v       v        |
  |   +------+ +------+ +-----+ +-----+ +-----+ +-----+ +------+   |
  |   |Onco- | |Bio-  | |CAR-T| |Imag-| |Auto-| |Phar-| |Cardi-|   |
  |   |logy  | |marker| |     | |ing  | |immu-| |maco-| |ology |   |
  |   |Agent | |Agent | |Agent| |Agent| |ne   | |gen. | |Agent |   |
  |   |:8503 | |:8502 | |:8504| |:8505| |:8506| |:8507| |:8527 |   |
  |   +------+ +------+ +-----+ +-----+ +-----+ +-----+ +------+   |
  |       |                                                          |
  |       v                                                          |
  |  Therapeutic Targets                                             |
  |       |                                                          |
  |       v                                                          |
  |  +-------------------------------+                               |
  |  | STAGE 3: Therapeutic Discovery|                               |
  |  | MolMIM + DiffDock + RDKit     |                               |
  |  +-------------------------------+                               |
  |       |                                                          |
  |       v                                                          |
  |  Ranked Drug Candidates (PDF)                                    |
  +-----------------------------------------------------------------+
```

Each stage is self-contained. You can run Stage 1 alone to get a VCF. You
can run Stage 2 alone (with a pre-existing VCF) to get clinical
intelligence. You can run Stage 3 alone (with a known target) to generate
drug candidates.

---

# Chapter 2: The Hardware -- NVIDIA DGX Spark

## 2.1 What Is DGX Spark?

The NVIDIA DGX Spark is a desktop workstation designed for AI workloads. It
is the smallest member of the DGX family -- the same product line used in
the world's largest supercomputers. At $4,699, it brings supercomputer-class
AI capabilities to a form factor that fits on a desk.

Key specifications:

```
  NVIDIA DGX Spark
  =================

  Chip:       GB10 Grace Blackwell
  GPU:        Blackwell-generation GPU cores
  CPU:        20 ARM Cortex cores (Grace)
  Memory:     128 GB unified LPDDR5x (shared CPU+GPU)
  Interconnect: NVLink-C2C (chip-to-chip, 900 GB/s)
  Storage:    NVMe SSD
  OS:         Ubuntu (ARM64)
  Price:      $4,699
```

## 2.2 Why GPU Matters for Medicine

A CPU (Central Processing Unit) is designed to execute complex instructions
one at a time, very fast. A GPU (Graphics Processing Unit) is designed to
execute simple instructions thousands at a time, in parallel.

**Analogy:** A CPU is a master chef who can prepare any dish perfectly but
works on one dish at a time. A GPU is a kitchen with 1,000 line cooks, each
trained to do one simple task -- chop, stir, plate -- simultaneously. When
you need to prepare one elaborate dish, the master chef wins. When you need
to prepare 1,000 simple dishes, the line cooks win by a landslide.

In precision medicine, nearly every computational task involves doing the
same operation millions or billions of times:

- **Sequence alignment:** Each of ~1 billion short reads must be compared
  against a 3.1-billion-letter reference genome. Every comparison is
  independent -- perfect for parallel execution.

- **Variant calling (DeepVariant):** A convolutional neural network examines
  each position in the genome. The same neural network architecture is
  applied to millions of positions independently.

- **Vector search:** Computing cosine similarity between a query vector and
  3.56 million stored vectors. Each comparison is independent.

- **Molecular docking:** Scoring how well each of 100+ candidate molecules
  fits into a protein binding pocket. Each molecule is scored independently.

In every case, the GPU's ability to run thousands of operations in parallel
turns hours into minutes.

## 2.3 Why "Unified Memory" Matters

In most computers, the CPU and GPU have separate pools of memory. Data must
be copied from CPU memory to GPU memory before the GPU can use it, and
results must be copied back. This copying is slow and wastes time.

The DGX Spark's GB10 chip uses **unified memory**: the CPU and GPU share
the same 128 GB of LPDDR5x RAM, connected by NVLink-C2C at 900 GB/s. No
copying is needed. The GPU can access the same data the CPU is using, and
vice versa.

For genomics, this means the entire reference genome, all sequencing reads,
and the neural network weights can sit in the same memory pool without
copying overhead.

## 2.4 What Runs Where

| Component              | Runs On  | Why                                    |
|------------------------|----------|----------------------------------------|
| BWA-MEM2 alignment     | GPU      | Billions of read-to-reference comparisons |
| DeepVariant calling    | GPU      | CNN inference on millions of positions |
| Variant annotation     | CPU      | Database lookups (ClinVar, AlphaMissense) |
| Vector embedding       | CPU/GPU  | BGE-small-en-v1.5 model inference      |
| Milvus vector search   | CPU+GPU  | ANN index search over 3.56M vectors    |
| Claude LLM inference   | API call | Runs on Anthropic cloud servers         |
| MolMIM generation      | GPU      | Generative neural network for molecules |
| DiffDock scoring       | GPU      | Diffusion model for binding prediction |
| RDKit drug-likeness    | CPU      | Chemical property calculations          |
| Streamlit / FastAPI     | CPU      | Web server and API handling             |
| Nextflow orchestration  | CPU      | Workflow coordination                  |
| Prometheus / Grafana    | CPU      | Monitoring and dashboards              |

## 2.5 Scaling Up

The same software runs unchanged on larger NVIDIA systems:

| System         | Price          | Use Case                               |
|----------------|----------------|----------------------------------------|
| DGX Spark      | $4,699         | Single researcher, clinic, or lab      |
| DGX B200       | ~$500K-$1M     | Hospital department, multiple patients |
| DGX SuperPOD   | $7M-$60M+     | Health system, population-scale studies |

The Docker containers, Nextflow workflows, and agent configurations are
identical across all three tiers. Only the hardware and resource limits
change.

---

# Chapter 3: Stage 1 -- Genomic Foundation Engine

## 3.1 What Is DNA Sequencing?

DNA is a molecule found in every cell of your body. It is a long string of
chemical "letters" -- just four of them: A (adenine), T (thymine),
C (cytosine), and G (guanine). The complete set of letters in your DNA is
called your **genome**, and it is approximately 3.1 billion letters long,
organized into 23 pairs of chromosomes.

Your genome is your body's instruction manual. It contains about 20,000
genes, each of which is a section of DNA that encodes the instructions for
building a specific protein. Proteins are the molecular machines that do
nearly everything in your body -- from carrying oxygen (hemoglobin) to
fighting infections (antibodies) to contracting muscles (actin and myosin).

**DNA sequencing** is the process of reading those 3.1 billion letters. A
machine called a DNA sequencer (typically made by Illumina) does this by:

1. Breaking the DNA into millions of small fragments.
2. Reading each fragment -- typically 150-250 letters at a time.
3. Using computers to reassemble the fragments into a complete picture.

The result is billions of short "reads" stored in a file called a FASTQ.

**Analogy:** Imagine shredding 1,000 copies of a 3.1-billion-page book into
strips of 250 characters each, then reassembling the book by finding
overlapping fragments. That is essentially what genomic sequencing and
alignment do. The reason you shred 1,000 copies (called "30x coverage") is
so that each position is read about 30 times, which helps catch errors in
any single reading.

## 3.2 The FASTQ Format

The FASTQ file is the starting point of the entire pipeline. It contains
the raw output of the DNA sequencer. Each read consists of exactly four
lines:

```
  Line 1: @READ_ID                    (identifier for this read)
  Line 2: ACGTACGTACGTACGT...          (the DNA sequence, 150-250 letters)
  Line 3: +                            (separator)
  Line 4: IIIIIHHHHHGGGGG...           (quality score for each letter)
```

Here is an actual example:

```
  @HWI-ST1234:100:ABC12AAXX:1:1101:1234:2100
  GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCG
  +
  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

The quality score on line 4 tells you how confident the sequencer is in each
letter. Each character maps to a number (Phred score). A score of `I`
corresponds to Phred 40, which means the probability of error is 1 in
10,000 -- very high confidence.

| Phred Score | Error Probability | Character | Interpretation     |
|-------------|-------------------|-----------|--------------------|
| 10          | 1 in 10           | `+`       | Poor               |
| 20          | 1 in 100          | `5`       | Acceptable         |
| 30          | 1 in 1,000        | `?`       | Good               |
| 40          | 1 in 10,000       | `I`       | Excellent          |

For a typical whole-genome sequencing run at 30x coverage, the FASTQ files
total approximately **200 GB** -- roughly the same as 50 HD movies.

## 3.3 Read Alignment (BWA-MEM2)

Once you have billions of short reads, the next step is figuring out where
each read came from in the genome. This is called **alignment** or
**mapping**.

The reference genome is a standardized "template" of human DNA called
GRCh38 (Genome Reference Consortium Human Build 38). It represents the
consensus sequence assembled from multiple individuals.

**BWA-MEM2** is the alignment tool used by the platform. For each of the
~1 billion reads in a FASTQ file, it:

1. Finds the most likely position in the 3.1-billion-letter reference genome
   where that read belongs.
2. Records how well the read matches (allowing for small differences that
   represent real variants versus sequencing errors).
3. Writes the result to a BAM (Binary Alignment Map) file.

On a traditional CPU, alignment takes 4-8 hours for a whole genome. On the
DGX Spark's GPU using NVIDIA Parabricks, it takes **20-45 minutes** -- a
10x speedup.

```
  Alignment
  =========

  FASTQ reads (billions)       Reference genome (GRCh38)
  ========================     ==========================
  ...ACGTACGTACGT...           ......ACGTACGTACGT.........
  ...TGCATGCATGCA...           Position: chr7:55,249,071
  ...GATCGATCGATC...

       |
       |  BWA-MEM2 (GPU-accelerated)
       v

  BAM file: each read mapped to its position in the reference
```

## 3.4 Variant Calling (DeepVariant)

With all reads aligned, the next step is identifying positions where the
patient's DNA differs from the reference. These differences are called
**variants**.

**Google DeepVariant** is a variant-calling tool that uses a convolutional
neural network (CNN) -- the same type of AI used for image recognition. It
works by converting the aligned reads at each position into a visual
"pileup image" and then using the CNN to classify what it sees:

```
  DeepVariant: How It Works
  ==========================

  Aligned reads at position chr7:55,249,071:

  Reference:  ...T G C A T G C A...
  Read 1:     ...T G C A T G C A...  (matches)
  Read 2:     ...T G C A T G C A...  (matches)
  Read 3:     ...T G C G T G C A...  (A -> G at this position!)
  Read 4:     ...T G C G T G C A...  (A -> G at this position!)
  Read 5:     ...T G C A T G C A...  (matches)

       |
       |  Convert to pileup image
       v
  +-------------------+
  |  CNN classifies:  |
  |  Heterozygous     |
  |  A -> G (SNV)     |
  |  QUAL: 35.2       |
  +-------------------+
```

DeepVariant recognizes several types of variants:

| Variant Type        | What Happened                     | Example              |
|---------------------|-----------------------------------|----------------------|
| SNV                 | Single letter changed             | A -> G at chr7:55M   |
| Insertion           | Extra letters added               | -- -> ACG            |
| Deletion            | Letters removed                   | ACG -> --            |
| Multi-nucleotide    | Multiple adjacent letters changed | AC -> GT             |

DeepVariant achieves **>99% accuracy** on well-characterized benchmark
genomes. On the DGX Spark GPU, variant calling takes **10-35 minutes**.

## 3.5 The VCF File

The output of variant calling is a **VCF** (Variant Call Format) file. This
is the lingua franca of genomics -- the standard format that every
downstream tool understands.

A VCF file has a header (lines starting with `#`) followed by one line per
variant:

```
  #CHROM  POS       ID           REF  ALT  QUAL    FILTER  INFO
  chr9    35065263  rs188935092  G    A    42.1    PASS    DP=45;AF=0.50
  chr7    55249071  rs121913529  T    G    38.7    PASS    DP=38;AF=0.45
  chr17   7674220   rs28934578   G    A    35.2    PASS    DP=52;AF=0.48
```

Key fields explained:

| Field    | Meaning                                                     |
|----------|-------------------------------------------------------------|
| CHROM    | Chromosome (chr1 through chr22, chrX, chrY)                 |
| POS      | Position in the chromosome (base-pair coordinate)           |
| ID       | Known variant identifier (e.g., rs188935092 from dbSNP)     |
| REF      | The reference allele (what the reference genome has)        |
| ALT      | The alternate allele (what the patient has instead)         |
| QUAL     | Quality score -- higher means more confidence               |
| FILTER   | PASS means it met quality thresholds                        |
| INFO     | Additional data: DP=read depth, AF=allele frequency         |

**For clinicians:** You do not need to read VCF files directly. The
intelligence agents parse them automatically and present the results in
plain-language reports. But understanding what a VCF contains helps you
evaluate the quality of the evidence behind a recommendation.

## 3.6 Variant Annotation

A raw VCF file tells you that a variant exists, but not what it means. Is
this variant disease-causing? Is it common in the population? Does it affect
protein function? **Annotation** adds this context.

The platform annotates each variant against three databases:

### ClinVar (4.1 million records)

ClinVar is a public database maintained by the NIH (National Institutes of
Health). Researchers and clinical labs submit their findings about specific
variants, and ClinVar classifies each one:

| Classification       | Meaning                                          |
|----------------------|--------------------------------------------------|
| Pathogenic           | Known to cause disease                           |
| Likely pathogenic    | Strong evidence of causing disease               |
| VUS                  | Variant of Uncertain Significance -- unknown     |
| Likely benign        | Probably harmless                                |
| Benign               | Known to be harmless                             |

### AlphaMissense (71 million predictions)

AlphaMissense is an AI tool from Google DeepMind (the same lab that created
AlphaFold, which won the 2024 Nobel Prize in Chemistry for predicting
protein structures). AlphaMissense predicts whether a missense variant
(a single amino acid change) is likely to cause disease.

Each variant receives a score from 0 to 1:
- Above **0.564**: Likely pathogenic
- Below **0.564**: Likely benign

### VEP -- Variant Effect Predictor

VEP (from Ensembl) classifies the functional impact of each variant:

| Impact Level | Meaning                              | Example            |
|-------------|--------------------------------------|--------------------|
| HIGH        | Probably destroys protein function   | Stop-gain, frameshift |
| MODERATE    | May change protein function          | Missense           |
| LOW         | Unlikely to change protein           | Synonymous         |
| MODIFIER    | Non-coding region, regulatory        | Intronic, UTR      |

### The Annotation Funnel

Starting with 11.7 million raw variants, the platform progressively filters:

```
  11,700,000  Raw variants detected
       |
       v
   3,500,000  High-quality (pass filter)
       |
       v
      35,616  Match ClinVar records
       |
       v
       6,831  Have AlphaMissense scores
       |
       v
      ~2,400  High impact + pathogenic
       |
       v
         847  In druggable genes (171 known targets)
```

## 3.7 By the Numbers

| Metric                          | Value                  |
|---------------------------------|------------------------|
| Raw variants per genome         | 11.7 million           |
| High-quality after filtering    | 3.5 million            |
| Variants in druggable genes     | 847                    |
| Alignment time (GPU)            | 20-45 minutes          |
| Alignment time (CPU)            | 4-8 hours              |
| Variant calling time (GPU)      | 10-35 minutes          |
| Variant calling time (CPU)      | 8-12 hours             |
| Total Stage 1 time (GPU)        | 120-240 minutes        |
| Total Stage 1 time (CPU)        | 24-48 hours            |
| FASTQ file size                 | ~200 GB                |
| Reference genome                | GRCh38                 |
| Coverage depth                  | 30x                    |
| DeepVariant accuracy            | >99%                   |

## 3.8 How This Connects to Stage 2

The VCF file is the handoff point between Stage 1 and Stage 2. But raw
variant records are not searchable by natural language. To make them
searchable, the platform:

1. **Annotates** each high-quality variant with ClinVar, AlphaMissense, and
   VEP data.
2. **Converts** each annotated variant into a text description (e.g.,
   "chr9:35065263 G>A in VCP gene, ClinVar Pathogenic, AlphaMissense 0.87,
   HIGH impact, associated with frontotemporal dementia").
3. **Embeds** each text description into a 384-dimensional vector using the
   BGE-small-en-v1.5 embedding model.
4. **Stores** the vector + metadata in the Milvus `genomic_evidence`
   collection.

The result: **3.56 million annotated variant vectors** that all eleven
intelligence agents can query using natural language.

```
  Stage 1 Output (VCF)
       |
       |  Annotate + Embed
       v
  +---------------------------------------------+
  |  Milvus: genomic_evidence collection         |
  |  3.56M vectors, 384 dimensions each          |
  |  Read-only access by all eleven agents         |
  +---------------------------------------------+
       |         |       |       |    |    |    |
       v         v       v       v    v    v    v
    Oncology  Biomarker CAR-T  Imaging ...  Cardiology
    Agent     Agent    Agent   Agent        Agent
```

---

# Chapter 4: Stage 2 -- Precision Intelligence Network (RAG/Chat)

## 4.1 The Data Challenge

To make a clinical decision about a single patient, a physician ideally
consults at least five major data sources:

```
  +-------------------+     +-------------------+     +-------------------+
  |  ClinVar          |     |  PubMed / PMC     |     | ClinicalTrials.gov|
  |                   |     |                   |     |                   |
  | 4.1M disease-     |     | 35M+ biomedical   |     | 400,000+ trial   |
  | variant records   |     | research articles |     | registrations     |
  +-------------------+     +-------------------+     +-------------------+
           |                         |                         |
           +----------+--------------+--------------+----------+
                      |                             |
              +-------------------+     +-------------------+
              |  AlphaMissense    |     |  Institutional    |
              |                   |     |  Knowledge        |
              | 71M AI-predicted  |     |                   |
              | pathogenicity     |     | Guidelines, case  |
              | scores            |     | histories, SOPs   |
              +-------------------+     +-------------------+
```

### The Problems

1. **Different formats.** ClinVar is a structured database. PubMed returns
   XML. Guidelines are often PDFs. AlphaMissense is a massive TSV file.
   No common query language spans all of them.

2. **Different vocabularies.** One source says "ERBB2." Another says "HER2."
   One says "myocardial infarction." Another says "heart attack." Keyword
   search misses these synonyms.

3. **Constant change.** New papers appear daily. Trial statuses change.
   Guidelines are updated. A search performed last month may be outdated.

4. **Volume.** ClinVar alone has 4.1 million records. AlphaMissense has 71
   million predictions. No human can read even a fraction of this manually.

5. **Cross-referencing.** The real insight comes from connecting variant data
   to protein function to drug targets to clinical evidence. No single
   source makes all these connections.

**Analogy:** Imagine a detective investigating a case where evidence is
spread across five different police departments in five different cities.
Each department uses a different filing system, a different case numbering
scheme, and a different language. The detective must travel to each city,
learn each system, find the relevant files, and piece together the story
manually. That is what a clinician does today when interpreting genomic
results across multiple databases. RAG is the system that unifies all five
departments into a single searchable archive.

## 4.2 What Is RAG?

RAG stands for **Retrieval-Augmented Generation**. It is a technique that
combines two technologies:

1. **A search engine** (Milvus vector database) that finds relevant
   documents from a knowledge base.
2. **A large language model** (Claude from Anthropic) that reads those
   documents and writes a coherent, cited answer.

Here is the step-by-step process:

```
  RAG in Six Steps
  =================

  Step 1: User asks a question
          "What are the therapeutic implications of a VCP R155H
           mutation in frontotemporal dementia?"

  Step 2: Question is embedded into a vector
          [0.023, -0.112, 0.087, ..., 0.045]  (384 numbers)

  Step 3: Vector is compared to 3.56M stored vectors
          Cosine similarity: find the closest matches

  Step 4: Top-K most relevant documents are retrieved
          "ClinVar: VCP R155H, Pathogenic, IBMPFD..."
          "AlphaMissense: VCP R155H, score 0.87..."
          "PubMed: Watts et al. 2004, VCP mutations in FTD..."

  Step 5: Documents + question are sent to Claude (LLM)
          System prompt + retrieved evidence + user question

  Step 6: Claude synthesizes a grounded answer with citations
          "The VCP R155H mutation is classified as Pathogenic
           by ClinVar [1] and scores 0.87 on AlphaMissense [2],
           indicating high pathogenicity. Current therapeutic
           approaches include..."
```

The critical word is **grounded**. The LLM does not make up information.
It can only cite evidence that was actually retrieved from the database.
This dramatically reduces hallucination compared to asking an LLM the
same question without any evidence.

## 4.3 Vector Embeddings Explained

A vector embedding is a list of numbers that captures the *meaning* of a
piece of text. The platform uses the BGE-small-en-v1.5 model, which
produces a vector of 384 numbers for any input text.

The key property: texts with similar meaning produce vectors that are
close together in 384-dimensional space, even if they use completely
different words.

```
  Text                              Vector (simplified to 3D for illustration)
  ================================  ==========================================
  "Heart failure"                   [0.82, 0.15, 0.43]
  "Cardiac insufficiency"           [0.80, 0.17, 0.41]   <-- very close!
  "Lung cancer"                     [0.12, 0.91, 0.33]   <-- far away
  "Non-small cell pulmonary
   adenocarcinoma"                  [0.14, 0.89, 0.35]   <-- close to lung cancer
```

**Analogy:** Think of embeddings like coordinates on a map. Paris and Lyon
are closer to each other than either is to Tokyo. Similarly, "heart failure"
and "cardiac insufficiency" are closer to each other than either is to "lung
cancer." The embedding model learns these relationships by reading millions
of biomedical texts during training.

This is why vector search is more powerful than keyword search. A keyword
search for "heart failure" would miss a document about "cardiac
insufficiency." A vector search finds it because the vectors are close.

## 4.4 Milvus Vector Database

Milvus is the vector database that stores and searches all embeddings in the
platform. It runs as a Docker container with two supporting services (etcd
for metadata and MinIO for object storage).

Key technical details:

| Parameter               | Value                   |
|-------------------------|-------------------------|
| Milvus version          | v2.4.0                  |
| Embedding dimensions    | 384                     |
| Embedding model         | BAAI/bge-small-en-v1.5  |
| Index type              | IVF_FLAT                |
| Similarity metric       | COSINE                  |
| Port                    | 19530                   |
| Memory limit            | 8 GB                    |

### How IVF_FLAT Works

IVF_FLAT stands for Inverted File with Flat storage. It is an approximate
nearest neighbor (ANN) algorithm that speeds up search by organizing vectors
into clusters:

1. **Training:** Milvus clusters all vectors into groups (called "Voronoi
   cells") based on similarity.
2. **Searching:** When a query arrives, Milvus first identifies which
   clusters are most likely to contain relevant vectors, then searches only
   those clusters instead of scanning all 3.56 million vectors.

This is not perfectly accurate (it might miss a relevant vector that ended
up in a different cluster), but it is dramatically faster -- searching
3.56 million vectors in milliseconds instead of seconds.

**Analogy:** Imagine a library with 3.56 million books. Searching every book
would take days. But if the books are organized into sections (cardiology,
oncology, neurology, etc.), you can go directly to the right section and
search only the relevant shelves. You might miss a cardiology book that was
accidentally shelved in neurology, but you save enormous amounts of time.

## 4.5 The Shared genomic_evidence Collection

The `genomic_evidence` collection is the backbone of the entire platform.
It contains **3.56 million variant vectors** generated by Stage 1, and it
is shared read-only by all eleven intelligence agents.

```
  genomic_evidence Collection
  ============================

  Records: 3,560,000 annotated variant vectors
  Access:  Read-only by all eleven agents
  Source:  Stage 1 genomic pipeline output
  Schema:  vector (384 dims) + metadata fields

  Metadata per record:
  - variant_id    (e.g., rs188935092)
  - chromosome    (e.g., chr9)
  - position      (e.g., 35065263)
  - ref_allele    (e.g., G)
  - alt_allele    (e.g., A)
  - gene          (e.g., VCP)
  - clinvar_class (e.g., Pathogenic)
  - alpha_missense_score (e.g., 0.87)
  - vep_impact    (e.g., HIGH)
  - consequence   (e.g., missense_variant)
```

Each agent adds its own **domain-specific collections** on top of the
shared genomic evidence. The oncology agent adds 10 collections covering
therapies, resistance, trials, and guidelines. The cardiology agent adds
12 collections covering risk calculators, medications, and cardiac imaging.
In total, the platform maintains approximately 80+ specialized collections
across all agents.

## 4.6 The Knowledge Graph

Beyond vector search, the platform maintains a structured knowledge graph
that maps relationships between biological entities:

```
  Knowledge Graph Structure
  ==========================

  201 Genes
       |
       v
  Proteins (encoded by genes)
       |
       v
  Pathways (signaling cascades: MAPK, PI3K/AKT, Wnt, etc.)
       |
       v
  Diseases (13 therapeutic areas)
       |
       v
  Drugs (171 druggable targets)
```

The 13 therapeutic areas span:

| # | Therapeutic Area         | Example Genes          |
|---|--------------------------|------------------------|
| 1 | Oncology                 | EGFR, BRAF, ALK, KRAS  |
| 2 | Cardiology               | SCN5A, MYBPC3, KCNQ1   |
| 3 | Neurology                | VCP, APP, PSEN1, MAPT  |
| 4 | Immunology / Autoimmune  | HLA-B, TNF, IL6, JAK2  |
| 5 | Hematology               | JAK2, CALR, MPL, BCR-ABL |
| 6 | Endocrinology            | GCK, HNF1A, INS        |
| 7 | Pulmonology              | CFTR, SERPINA1, SFTPC   |
| 8 | Nephrology               | PKD1, PKD2, COL4A5     |
| 9 | Gastroenterology         | APC, MLH1, MSH2        |
| 10| Dermatology              | FLG, KRT14, COL7A1     |
| 11| Ophthalmology            | RHO, RPE65, PAX6       |
| 12| Musculoskeletal          | DMD, SMN1, COL1A1      |
| 13| Rare Disease             | CFTR, SMN1, GBA, HTT   |

Of the 201 genes in the knowledge graph, **171 (85%)** are known to be
"druggable" -- meaning there is a known mechanism to target them with a
small molecule, antibody, or gene therapy.

## 4.7 Claude AI -- The Synthesis Layer

Claude is a large language model (LLM) from Anthropic that serves as the
reasoning and synthesis engine for all 11 intelligence agents. It is
important to understand exactly what Claude does and does not do in this
system.

### What Claude Does

- **Synthesizes** evidence from multiple retrieved documents into a coherent
  narrative.
- **Cites** specific sources so the clinician can verify every claim.
- **Explains** complex genomic findings in language appropriate for the
  audience (clinician, researcher, or patient).
- **Structures** output into the required format (Markdown report, FHIR R4,
  JSON, PDF).

### What Claude Does NOT Do

- **Does NOT calculate** risk scores. Those are computed by deterministic
  algorithms (e.g., ASCVD Pooled Cohort Equation, CHA2DS2-VASc score).
- **Does NOT make** diagnostic decisions. It presents evidence; the
  clinician decides.
- **Does NOT hallucinate** drug names or dosages -- it can only reference
  drugs that appear in the retrieved evidence.
- **Does NOT store** patient data. All inference happens via API calls to
  Anthropic's servers, and the platform can be configured to use local
  models for environments requiring data sovereignty.

**Analogy:** Claude is like a research assistant who has read every paper in
the library. When you ask a question, the assistant does not guess -- it
goes to the shelves, pulls the relevant papers, reads them, and writes you
a summary with page numbers. If the relevant paper is not in the library,
the assistant says "I could not find evidence for this" rather than making
something up.

## 4.8 By the Numbers

| Metric                               | Value                  |
|--------------------------------------|------------------------|
| Annotated variant vectors            | 3.56 million           |
| ClinVar records available            | 4.1 million            |
| AlphaMissense predictions available  | 71 million             |
| Embedding dimensions                 | 384                    |
| Embedding model                      | BGE-small-en-v1.5      |
| Query response time                  | < 5 seconds            |
| Genes in knowledge graph             | 201                    |
| Druggable targets                    | 171 (85%)              |
| Therapeutic areas                    | 13                     |
| Total agent collections              | 80+                    |
| LLM                                  | Claude (Anthropic)     |
| Vector database                      | Milvus v2.4.0          |

---

# Chapter 5: Stage 3 -- Therapeutic Discovery Engine

## 5.1 From Target to Drug Candidate

Stage 2 identifies therapeutic targets -- specific proteins that are
altered by disease-causing variants and can potentially be blocked or
modulated by a drug. Stage 3 takes those targets and generates novel
molecules that could bind to them.

In traditional pharmaceutical research, this process takes years and costs
hundreds of millions of dollars. The HCLS AI Factory performs the initial
computational steps -- molecular generation, docking simulation, and
drug-likeness scoring -- in **8-16 minutes**.

This does not replace the full drug development pipeline (which includes
laboratory testing, animal studies, and clinical trials). It accelerates
the very first step: identifying promising chemical starting points.

```
  Stage 3: From Target to Candidates
  ====================================

  Protein target (from Stage 2)
       |
       v
  Retrieve 3D structure (PDB/AlphaFold)
       |
       v
  Find seed compound (existing drug/inhibitor)
       |
       v
  Generate 100+ novel analogues (MolMIM)
       |
       v
  Score binding affinity (DiffDock)
       |
       v
  Evaluate drug-likeness (RDKit)
       |
       v
  Ranked candidates with PDF report
```

## 5.2 Protein Structures

Before you can design a drug to fit into a protein, you need to know what
the protein looks like in three dimensions.

Proteins are long chains of amino acids (typically hundreds to thousands of
them) that fold into complex 3D shapes. The shape determines the function.
A protein's **binding site** is a pocket or groove on its surface where
other molecules can attach.

**Analogy:** A protein is like a 3D lock. Drug discovery is the process of
finding -- or making -- the right key. The binding site is the keyhole.
If the key fits precisely, it can turn the lock (activate or block the
protein). If it does not fit, nothing happens.

Protein structures are determined experimentally using:
- **X-ray crystallography** -- shooting X-rays at protein crystals
- **Cryo-EM** (cryo-electron microscopy) -- imaging frozen protein samples

These structures are stored in the **RCSB Protein Data Bank (PDB)**, a
public database with over 200,000 entries. Each structure has a resolution
measured in angstroms (A) -- lower is better:

| Resolution | Quality        | Suitable For       |
|------------|----------------|--------------------|
| < 2.0 A    | Excellent      | Drug design        |
| 2.0-3.0 A  | Good           | Drug design        |
| 3.0-4.0 A  | Moderate       | Docking studies    |
| > 4.0 A    | Low            | General shape only |

For the platform's VCP demo, the structure 5FTK (resolution 2.5 A) is used
because it shows VCP bound to the existing inhibitor CB-5083, revealing the
exact binding pocket that new drugs should target.

## 5.3 Molecular Representation: SMILES

Molecules in the pipeline are represented using **SMILES** (Simplified
Molecular-Input Line-Entry System) -- a text notation that encodes a 3D
molecular structure as a string of characters.

```
  SMILES Examples
  ================

  Water:           O
  Ethanol:         CCO
  Aspirin:         CC(=O)Oc1ccccc1C(=O)O
  Caffeine:        Cn1c(=O)c2c(ncn2C)n(C)c1=O
  CB-5083 (VCP     COc1ccc(cn1)c2cc(NC(=O)c3cccc(c3)C(F)(F)F)ccc2
    inhibitor):
```

SMILES notation follows simple rules: uppercase letters are atoms (C, N, O,
S, F), lowercase letters in parentheses represent aromatic rings, `=` is a
double bond, and numbers indicate ring closures. The advantage of SMILES is
that it lets AI models treat molecules like text -- the same techniques used
for language generation can be applied to molecular generation.

## 5.4 Generative Molecular Design (BioNeMo MolMIM)

**BioNeMo MolMIM** (Molecular Masked Inverse Modeling) is an NVIDIA AI
model that generates novel molecules. It uses the same principle as masked
language modeling -- the technique behind models like BERT and Claude --
but applied to molecular structures instead of words.

Here is how it works:

1. **Input:** A seed molecule (e.g., CB-5083, the existing VCP inhibitor)
   represented as a SMILES string.
2. **Masking:** The model randomly masks (hides) parts of the molecule.
3. **Prediction:** The model predicts what could fill in the masked regions,
   exploring chemical space.
4. **Output:** 100+ novel molecules that are structurally related to the
   seed but chemically distinct.

```
  MolMIM Generation
  ==================

  Seed: COc1ccc(cn1)c2cc(NC(=O)c3cccc(c3)C(F)(F)F)ccc2  (CB-5083)
       |
       |  Mask + Predict (x100)
       v
  Candidate 1: COc1ccc(cn1)c2cc(NC(=O)c3cccc(c3)C(Cl)F)ccc2
  Candidate 2: COc1ccc(cn1)c2cc(NC(=O)c3ccnc(c3)C(F)(F)F)ccc2
  Candidate 3: COc1cnc(cn1)c2cc(NC(=O)c3cccc(c3)CF)ccc2
  ...
  Candidate 100: ...
```

Each candidate is a valid chemical structure that could potentially be
synthesized in a laboratory. The generation process takes **2-5 minutes**
for 100+ candidates.

## 5.5 Molecular Docking (DiffDock)

Once you have 100+ candidate molecules, you need to know which ones
actually fit into the target protein's binding pocket. This is called
**molecular docking**.

**BioNeMo DiffDock** uses a diffusion model (the same class of AI used in
image generators like DALL-E and Stable Diffusion) to predict the 3D
binding pose of each molecule in the protein's pocket. It outputs a
**docking score** in kcal/mol -- a measure of binding energy.

The more negative the score, the stronger the predicted binding:

| Docking Score    | Interpretation         |
|------------------|------------------------|
| -6 to -7 kcal/mol | Weak binding          |
| -7 to -8 kcal/mol | Moderate binding      |
| -8 to -10 kcal/mol | Strong binding       |
| -10 to -12 kcal/mol | Very strong binding |
| < -12 kcal/mol   | Exceptional (verify)   |

```
  DiffDock: Predicting Binding
  =============================

  Protein (VCP)                     Candidate molecule
  +------------------+              +--------+
  |                  |              |        |
  |    +--------+    |              | (drug) |
  |    | Binding|    |    Dock      |        |
  |    | Pocket |<---|-----------   +--------+
  |    |   ?    |    |
  |    +--------+    |
  |                  |
  +------------------+

  Score: -11.4 kcal/mol  (very strong!)
```

Docking takes **5-10 minutes** for 100 candidates, running on the GPU.

## 5.6 Drug-Likeness Scoring

Not every molecule that binds to a protein would make a good drug. A
molecule might bind perfectly but be impossible to manufacture, toxic to
the liver, or unable to survive the digestive system.

**RDKit** is an open-source chemistry toolkit that calculates drug-likeness
properties for each candidate.

### Lipinski's Rule of Five

Christopher Lipinski (a medicinal chemist at Pfizer) observed that most
successful oral drugs share four properties. A molecule that violates more
than one of these rules is unlikely to work as a pill:

| Property                  | Threshold         | Why It Matters              |
|---------------------------|-------------------|-----------------------------|
| Molecular weight          | <= 500 Da         | Larger molecules are harder to absorb |
| LogP (fat solubility)     | <= 5              | Too greasy = poor solubility |
| H-bond donors             | <= 5              | Too many = cannot cross membranes |
| H-bond acceptors          | <= 10             | Too many = cannot cross membranes |

### QED -- Quantitative Estimate of Drug-likeness

QED combines multiple drug-likeness properties into a single score from
0 to 1:

| QED Score | Interpretation          |
|-----------|-------------------------|
| 0.0 - 0.3 | Poor drug-likeness     |
| 0.3 - 0.5 | Below average          |
| 0.5 - 0.67 | Acceptable            |
| 0.67 - 0.85 | Good (drug-like)     |
| 0.85 - 1.0 | Excellent             |

### Additional Metrics

| Metric          | What It Measures                         | Good Range   |
|-----------------|------------------------------------------|--------------|
| TPSA            | Topological polar surface area (membrane crossing) | 20-140 A^2 |
| SA Score        | Synthetic accessibility (ease of manufacture) | 1-4 (easy)   |
| Rotatable bonds | Molecular flexibility                    | <= 10        |

## 5.7 The Full Pipeline

Putting all the pieces together:

```
  Stage 3: Drug Discovery Pipeline
  ==================================

  Input: Protein target + seed compound
         (from Stage 2 therapeutic intelligence)
         |
         v
  +---------------------------+
  | 1. Structure Retrieval    |  Query PDB for target protein structure
  |    (RCSB API)             |  Select best resolution with bound ligand
  +---------------------------+
         |
         v
  +---------------------------+
  | 2. Molecule Generation    |  Generate 100+ novel analogues
  |    (BioNeMo MolMIM)       |  from seed compound SMILES
  |    2-5 min, GPU           |
  +---------------------------+
         |
         v
  +---------------------------+
  | 3. Binding Prediction     |  Predict binding pose and score
  |    (BioNeMo DiffDock)     |  for each candidate in protein pocket
  |    5-10 min, GPU          |
  +---------------------------+
         |
         v
  +---------------------------+
  | 4. Drug-Likeness          |  Calculate Lipinski, QED, TPSA,
  |    (RDKit)                |  SA score, rotatable bonds
  |    < 1 min, CPU           |
  +---------------------------+
         |
         v
  +---------------------------+
  | 5. Composite Scoring      |  Rank candidates by weighted score:
  |    & Ranking              |  30% generation + 40% docking + 30% QED
  +---------------------------+
         |
         v
  +---------------------------+
  | 6. Report Generation      |  PDF report with ranked candidates,
  |    (Python + LaTeX)       |  2D/3D visualizations, and scores
  +---------------------------+
```

## 5.8 By the Numbers

| Metric                             | Value                 |
|------------------------------------|-----------------------|
| Candidates per seed compound       | 100+                  |
| Molecule generation time           | 2-5 minutes           |
| Docking time (100 candidates)      | 5-10 minutes          |
| Drug-likeness calculation          | < 1 minute            |
| Total Stage 3 time                 | 8-16 minutes          |
| Composite scoring weights          | 30/40/30 (gen/dock/QED) |
| Demo target                        | VCP protein (5FTK)    |
| Demo seed compound                 | CB-5083               |
| Best demo docking score            | -11.4 kcal/mol        |
| Best demo QED                      | 0.81                  |
| Improvement over seed (composite)  | 39%                   |

## 5.9 Limitations

Stage 3 is a computational hypothesis generator. It answers the question:
"What molecules *might* work?" It does not answer: "Will this molecule
actually cure the disease?"

Before any computationally designed molecule can become a drug, it must
pass through:

1. **In vitro testing** -- does it work in a test tube?
2. **In vivo testing** -- does it work in animal models?
3. **Phase I clinical trial** -- is it safe in humans?
4. **Phase II clinical trial** -- does it work in patients?
5. **Phase III clinical trial** -- does it work better than existing drugs?
6. **FDA/EMA approval** -- regulatory review.

This process typically takes 10-15 years and costs $1-2 billion. The
HCLS AI Factory accelerates Step 0 -- generating the initial candidates --
from months (traditional medicinal chemistry) to minutes (AI-based
generation).

---

# Chapter 6: The Intelligence Agent Architecture

## 6.1 What Is an Intelligence Agent?

An intelligence agent is an autonomous reasoning system that goes beyond
simple question-answering. Unlike a basic chatbot that generates a single
response to a single question, an intelligence agent:

- **Plans** a multi-step search strategy before executing
- **Searches** across multiple data collections in parallel
- **Evaluates** the quality and completeness of retrieved evidence
- **Synthesizes** findings into a coherent, cited narrative
- **Reports** in multiple formats (Markdown, PDF, FHIR R4, JSON)

Each agent is specialized for a clinical domain (oncology, cardiology,
etc.) but shares the same underlying architecture. All 11 agents in the
HCLS AI Factory follow the same five-phase pattern.

**Analogy:** Think of each agent as a specialist physician who also happens
to be a tireless researcher. When you ask the oncology agent about a KRAS
G12C mutation, it does not just look up "KRAS G12C" in a single textbook.
It formulates a research plan ("I need to check variant databases,
literature, guidelines, trials, and resistance mechanisms"), executes the
plan across 11 different collections, evaluates whether the evidence is
sufficient, writes a comprehensive report with citations, and formats it
for the tumor board. All in under 5 seconds.

## 6.2 The Agent Pattern

All eleven agents follow the same five-phase loop:

```
  The Five-Phase Agent Loop
  ==========================

  Phase 1: PLAN
  +---------------------------+
  | - Parse the user question |
  | - Identify entities       |
  |   (genes, variants, drugs)|
  | - Select search strategy  |
  | - Choose collections      |
  +---------------------------+
           |
           v
  Phase 2: SEARCH
  +---------------------------+
  | - Query domain-specific   |
  |   collections in parallel |
  | - Query shared            |
  |   genomic_evidence        |
  | - Apply collection weights|
  | - Retrieve Top-K results  |
  +---------------------------+
           |
           v
  Phase 3: EVALUATE
  +---------------------------+
  | - Score evidence quality   |
  | - Check completeness      |
  | - Flag gaps or conflicts  |
  | - Re-search if needed     |
  +---------------------------+
           |
           v
  Phase 4: SYNTHESIZE
  +---------------------------+
  | - Send evidence + question|
  |   to Claude LLM           |
  | - Generate structured     |
  |   response with citations |
  | - Apply domain-specific   |
  |   reasoning templates     |
  +---------------------------+
           |
           v
  Phase 5: REPORT
  +---------------------------+
  | - Format output           |
  |   (Markdown, PDF, FHIR,  |
  |    JSON)                  |
  | - Publish events (SSE)    |
  | - Cache for future queries|
  +---------------------------+
```

### Phase 1: Plan

The agent analyzes the incoming question using entity recognition:

- **Gene names:** EGFR, BRAF, VCP, SCN5A, etc.
- **Variant identifiers:** V600E, R155H, rs188935092, etc.
- **Drug names:** osimertinib, sotorasib, entresto, etc.
- **Disease terms:** NSCLC, heart failure, FTD, etc.

Based on the entities found, the agent selects which collections to search
and in what order. A question about drug resistance would prioritize
resistance and therapy collections. A question about clinical trials would
prioritize trial collections.

### Phase 2: Search

The agent executes parallel vector searches across its domain-specific
collections plus the shared `genomic_evidence` collection. Each collection
has a weight that determines how much its results influence the final
answer.

Typical search parameters:
- **Top-K:** 10-30 results per collection
- **Similarity threshold:** Cosine similarity > 0.7
- **Parallel execution:** All collections searched simultaneously

### Phase 3: Evaluate

The agent checks whether the retrieved evidence is sufficient:

- Are all identified entities covered in the results?
- Are there conflicting pieces of evidence that need reconciliation?
- Is the evidence recent enough (publication date)?
- Is there enough evidence from high-quality sources?

If gaps are found, the agent may execute a second search round with
modified parameters (broader terms, different collections, lower
similarity threshold).

### Phase 4: Synthesize

The retrieved evidence is passed to Claude along with a domain-specific
system prompt that instructs the model to:

- Ground all claims in the retrieved evidence
- Cite specific sources by number
- Structure the response according to the domain template
- Flag uncertainties explicitly
- Use terminology appropriate for the requesting persona

### Phase 5: Report

The agent formats the synthesized response into one or more output formats:

| Format    | Use Case                                    |
|-----------|---------------------------------------------|
| Markdown  | Display in Streamlit UI, copy to EHR notes  |
| PDF       | Tumor board packets, clinical reports        |
| FHIR R4   | Integration with EHR systems (Epic, Cerner) |
| JSON      | Programmatic access, downstream pipelines    |

## 6.3 How Agents Share Data

The intelligence agents operate independently but share data through two
mechanisms:

### Mechanism 1: Shared Collections (Read-Only)

The `genomic_evidence` collection (3.56M variant vectors) is populated by
Stage 1 and consumed by all 11 agents. No agent writes to this
collection after the initial population -- it is strictly read-only during
operation.

### Mechanism 2: Cross-Agent Event Publishing (SSE)

Agents can publish events via Server-Sent Events (SSE) to notify other
agents of findings. For example:

- The **imaging agent** detects a suspicious lung nodule on a CT scan and
  publishes an event: `{"type": "finding", "modality": "CT",
  "region": "lung", "finding": "nodule_5mm"}`.
- The **oncology agent** subscribes to imaging events and automatically
  searches for lung-cancer-associated genomic variants in the same
  patient's `genomic_evidence` data.

This cross-modal triggering enables workflows that span multiple clinical
domains without manual intervention.

```
  Cross-Agent Data Flow
  ======================

  +---------------------+
  | genomic_evidence    |  <--- Written by Stage 1 (once)
  | 3.56M variants      |
  +---------------------+
    |    |    |    |    |    |    |
    v    v    v    v    v    v    v
  +----+----+----+----+----+----+----+
  |Onc |Bio |CART|Img |Auto|PGx |Card|+4  <-- 11 agents read
  +----+----+----+----+----+----+----+
    |         ^         |
    |   SSE   |   SSE   |
    +---------+---------+               <--- Cross-agent events
```

## 6.4 The Eleven Intelligence Agents at a Glance

| Agent | Port | Domain Collections | Domain | Key Differentiator |
|-------|------|--------------------|--------|--------------------|
| Precision Oncology | 8503/8103 | 10 | Molecular tumor board decision support | Therapy ranking with resistance awareness, AMP/ASCO/CAP tiering |
| Precision Biomarker | 8502/8102 | 10 | Genotype-aware biomarker interpretation | Biological age estimation (PhenoAge/GrimAge), disease trajectory |
| CAR-T Intelligence | 8504/8104 | 11 | Cellular immunotherapy intelligence | Construct comparison (4-1BB vs CD28), manufacturing protocols |
| Imaging Intelligence | 8505/8105 | 10 | Medical imaging AI (CT, MRI, X-ray) | NVIDIA NIM integration (VISTA-3D, MAISI, VILA-M3), DICOM ingestion |
| Precision Autoimmune | 8506/8106 | 10 | Autoimmune and immune-mediated conditions | Immune pathway analysis, flare prediction |
| Pharmacogenomics | 8507/8107 | 15 | Drug-gene interaction and dosing | 25 pharmacogenes, 100+ drugs, 9 dosing algorithms, 15 HLA associations |
| Cardiology Intelligence | 8527/8126 | 12 | Cardiovascular clinical decision support | 6 risk calculators (ASCVD, HEART, CHA2DS2-VASc), GDMT optimizer |
| Clinical Trial Intelligence | 8538/8128 | 10 | Clinical trial matching and enrollment | Trial eligibility matching, protocol analysis, enrollment optimization |
| Rare Disease Diagnostic | 8134/8544 | 10 | Rare disease differential diagnosis | Gene panel analysis, phenotype-genotype correlation, diagnostic odyssey reduction |
| Neurology Intelligence | 8528/8529 | 10 | Neurological condition assessment | Neurodegeneration pathways, treatment planning, cognitive assessment integration |
| Single-Cell Intelligence | 8540/8130 | 10 | Single-cell transcriptomics analysis | Cell-type identification, expression profiling, spatial transcriptomics |

### Architecture per Agent

Every agent follows the same deployment pattern:

```
  Agent Deployment Pattern
  =========================

  +------------------+     +------------------+
  |  Streamlit UI    |     |  FastAPI Backend  |
  |  (port 850x)     |<--->|  (port 810x)     |
  +------------------+     +------------------+
                                  |
                                  v
                           +-------------+
                           |   Milvus    |
                           | (port 19530)|
                           +-------------+
                                  |
                           +------+------+
                           |             |
                      Domain         Shared
                      Collections    genomic_evidence
                      (10-15 each)   (3.56M vectors)
```

Each agent runs as a Docker container with:
- **Streamlit frontend** for interactive clinical use
- **FastAPI backend** for programmatic API access
- **Health check endpoint** monitored by the platform watchdog
- **Pydantic models** for input validation and output serialization

## 6.5 The Technology Stack

For reference, here is the complete technology stack underlying the
platform:

```
  Technology Stack
  =================

  Layer 1: Hardware
  +-----------------------------------------------+
  | NVIDIA DGX Spark (GB10, 128GB, 20 ARM cores)  |
  +-----------------------------------------------+

  Layer 2: Operating System
  +-----------------------------------------------+
  | Ubuntu (ARM64) + CUDA 12.x + Docker            |
  +-----------------------------------------------+

  Layer 3: Infrastructure
  +-------+--------+--------+---------+------------+
  | Milvus | etcd   | MinIO  | Nextflow| Prometheus |
  | v2.4.0 |        |        | DSL2    | + Grafana  |
  +-------+--------+--------+---------+------------+

  Layer 4: Genomics
  +----------+--------------+----------+
  | Parabricks| DeepVariant  | BWA-MEM2 |
  | 4.6       | (CNN)        |          |
  +----------+--------------+----------+

  Layer 5: AI / LLM
  +---------+------------------+----------+
  | Claude  | BGE-small-en-v1.5| BioNeMo  |
  | (Anthro)| (embeddings)     | NIMs     |
  +---------+------------------+----------+

  Layer 6: Chemistry
  +---------+----------+---------+
  | MolMIM  | DiffDock | RDKit   |
  +---------+----------+---------+

  Layer 7: Applications
  +----------+----------+--------+
  | Streamlit| FastAPI  | Flask  |
  | (UIs)    | (APIs)   | (Hub)  |
  +----------+----------+--------+

  Layer 8: Shared Library
  +-----------------------------------------------+
  | hcls_common (23 modules): config, milvus, LLM,|
  | security, embedding, health, logging, etc.     |
  +-----------------------------------------------+
```

---

## What Comes Next: Part 2

Part 2 of this guide covers each of the eleven intelligence agents in
depth:

- **Chapter 7:** Imaging Intelligence Agent -- NVIDIA NIM models, DICOM
  workflows, cross-modal triggers.
- **Chapter 8:** Precision Oncology Agent -- Molecular tumor boards,
  therapy ranking, trial matching, resistance awareness.
- **Chapter 9:** Precision Biomarker Agent -- Biological age, disease
  trajectory, pharmacogenomic profiling.
- **Chapter 10:** CAR-T Intelligence Agent -- Construct design, manufacturing
  protocols, comparative analysis.
- **Chapter 11:** Precision Autoimmune Agent -- Immune pathways, flare
  prediction, treatment response.
- **Chapter 12:** Pharmacogenomics Agent -- Drug-gene interactions, dosing
  algorithms, HLA associations.
- **Chapter 13:** Cardiology Intelligence Agent -- Risk calculators, GDMT
  optimization, cardiac genomics.
- **Chapter 14:** Clinical Trial Intelligence Agent -- Trial matching,
  eligibility analysis, enrollment optimization.
- **Chapter 15:** Rare Disease Diagnostic Agent -- Differential diagnosis,
  gene panel analysis, diagnostic odyssey reduction.
- **Chapter 16:** Neurology Intelligence Agent -- Neurodegeneration pathways,
  treatment planning, cognitive assessment.
- **Chapter 17:** Single-Cell Intelligence Agent -- Cell-type identification,
  expression profiling, spatial transcriptomics.

Each chapter follows the same structure: domain introduction, clinical
workflow, data collections, agent-specific features, example queries,
and output interpretation.

---

## Review Questions (Part 1)

### Chapter 1: The Precision Medicine Revolution

1. What is the fundamental difference between traditional medicine and
   precision medicine?
2. How long does the traditional precision medicine pipeline take, and how
   does the HCLS AI Factory compress that timeline?
3. What are the three stages of the HCLS AI Factory pipeline?

### Chapter 2: The Hardware

4. What is the DGX Spark's price, and what type of chip does it use?
5. Why is unified memory important for genomics workloads?
6. Name two tasks that run on the GPU and two that run on the CPU.

### Chapter 3: Genomic Foundation Engine

7. What is a FASTQ file, and approximately how large is one for a whole
   genome at 30x coverage?
8. What does BWA-MEM2 do, and how much faster is it on GPU vs CPU?
9. How does DeepVariant use a convolutional neural network to call variants?
10. What are the key fields in a VCF file?
11. Describe the annotation funnel: how many variants survive each
    filtering step?

### Chapter 4: Precision Intelligence Network

12. What does RAG stand for, and what are its two core technologies?
13. Why is vector search better than keyword search for medical queries?
14. What is the `genomic_evidence` collection, and how many vectors does
    it contain?
15. What does Claude do in the RAG pipeline, and what does it NOT do?
16. How many genes are in the knowledge graph, and what percentage are
    druggable?

### Chapter 5: Drug Discovery

17. What is SMILES notation, and why is it useful for AI-based drug design?
18. What is a docking score, and what range indicates strong binding?
19. Name three of Lipinski's Rule of Five criteria.
20. What is the composite scoring formula for ranking drug candidates?

### Chapter 6: Intelligence Agents

21. Name the five phases of the agent reasoning loop.
22. How do agents share genomic data, and how do they communicate events?
23. How many intelligence agents does the platform include, and what
    domains do they cover?

---

# PART 2 -- The Intelligence Agents

> Part 2 covers the eleven intelligence agents in clinical depth. Each chapter
> follows the same pedagogical pattern: concept, analogy, clinical relevance,
> technical detail, and agent integration. Part 2 concludes with a complete
> patient journey, glossary, and quick reference.

---

# Chapter 7: Imaging Intelligence Agent

## 7.1 Medical Imaging for Non-Specialists

Medical imaging lets clinicians see inside the body without surgery. Four main modalities dominate cardiac and general radiology:

**Analogy:** Think of imaging modalities as different types of cameras. X-ray is a flash photo -- it captures bones and dense structures. CT is a 3D scanner -- it takes hundreds of cross-sectional slices. MRI is a long-exposure shot in different lighting -- it reveals soft tissue in extraordinary detail. Nuclear imaging is a heat camera -- it shows metabolic activity, not anatomy.

| Modality | What It Shows | Radiation | Speed | Best For |
|---|---|---|---|---|
| X-ray / CXR | Bones, lung fields, heart silhouette | Low | Seconds | Rapid screening, pneumonia, fractures |
| CT | Cross-sectional anatomy, vascular detail | Moderate | Minutes | Trauma, PE, lung nodules, coronary calcium |
| MRI | Soft tissue contrast, T1/T2 relaxation | None | 30-60 min | Brain, cardiac tissue characterization, joints |
| Nuclear (SPECT/PET) | Metabolic activity, perfusion | Low-Moderate | 30-60 min | Myocardial perfusion, cancer staging, sarcoidosis |

### 7.2 AI in Radiology: NVIDIA NIMs

The Imaging Intelligence Agent integrates three NVIDIA NIM microservices:

- **VISTA-3D** (port 8530): Segments 132 anatomical classes in 3D CT volumes -- from liver and kidneys to individual vertebrae. Think of it as automatic organ labeling.
- **MAISI** (port 8531): Generates synthetic medical images for training and augmentation. Useful when real patient data is limited or restricted.
- **VILA-M3** (port 8532): A vision-language model that can describe what it sees in a medical image, bridging the gap between pixels and clinical reports.

```
  VISTA-3D Segmentation
  ======================

  Input: CT volume (512 x 512 x 300 slices)
         |
         v
  VISTA-3D (132-class model)
         |
         v
  Output: Labeled volume
         - Liver  (class 5)  -> volume: 1,450 mL
         - Spleen (class 6)  -> volume: 210 mL
         - L.Kidney (class 7) -> volume: 165 mL
         - R.Kidney (class 8) -> volume: 158 mL
         - 128 more structures...
```

### 7.3 The 10 Collections

| # | Collection | Records | Purpose |
|---|---|---|---|
| 1 | imaging_literature | 2,678 | PubMed radiology research papers |
| 2 | imaging_trials | 12 | Active imaging clinical trials |
| 3 | imaging_findings | 124 | Seed reference findings (normal/abnormal patterns) |
| 4 | imaging_protocols | -- | Acquisition parameters by modality |
| 5 | imaging_devices | -- | FDA-cleared AI/ML radiology devices |
| 6 | imaging_anatomy | -- | Anatomical reference for segmentation classes |
| 7 | imaging_benchmarks | -- | Model performance (LUNA, TCIA, PhysioNet datasets) |
| 8 | imaging_guidelines | -- | ACR/RSNA appropriateness criteria |
| 9 | imaging_report_templates | -- | Structured reporting templates |
| 10 | imaging_datasets | -- | Public dataset catalog (TCIA, PhysioNet) |
| +1 | genomic_evidence | 3.5M | Shared genomic variants (read-only) |

### 7.4 Clinical Workflows

| Workflow | Target Latency | What It Does |
|---|---|---|
| CT Head Hemorrhage | <90 seconds | Detects intracranial hemorrhage, classifies type, alerts stroke team |
| CXR Rapid Findings | <30 seconds | Screens for pneumonia, cardiomegaly, pleural effusion, pneumothorax |
| CT Lung Nodule | <5 minutes | Measures nodule size, tracks growth, applies Lung-RADS classification |
| MRI Brain MS | <5 minutes | Detects demyelinating lesions, tracks lesion load over time |

### 7.5 Cross-Modal Integration

When the Imaging Agent detects a high-risk lung nodule (Lung-RADS 4A+), it automatically queries the shared genomic_evidence collection for relevant mutations -- EGFR, ALK, ROS1, KRAS -- bridging radiology and molecular pathology.

### 7.6 Federated Learning

The agent supports NVIDIA FLARE (Federated Learning Application Runtime
Environment) for privacy-preserving multi-site model training. Hospitals can
improve AI models collaboratively without sharing patient data -- only model
weight updates are exchanged.

```
  Federated Learning with FLARE
  ==============================

  Hospital A          Hospital B          Hospital C
  +----------+        +----------+        +----------+
  | Local    |        | Local    |        | Local    |
  | Training |        | Training |        | Training |
  | (weights)|        | (weights)|        | (weights)|
  +----+-----+        +----+-----+        +----+-----+
       |                    |                    |
       +--------> Central Aggregator <-----------+
                     |
                     v
               Global Model (shared back)
```

### 7.7 How the Agent Processes a Query

```
User Query: "6mm lung nodule management"
    |
    v
[Query Expansion] → "lung nodule" + "Lung-RADS" + "follow-up" + "screening"
    |
    v
[Parallel Search] → imaging_guidelines (0.20) + imaging_literature (0.18)
                   + imaging_findings (0.15) + imaging_protocols (0.12)
    |
    v
[Evidence Scoring] → Top-K results ranked by cosine similarity
    |
    v
[Cross-Modal Check] → Lung-RADS 4A? → Query genomic_evidence for EGFR/ALK
    |
    v
[Claude Synthesis] → Guideline-grounded response with citations
```

### 7.8 Sample Response

```json
{
  "answer": "A 6mm solid pulmonary nodule in a low-risk patient (Lung-RADS 3) should be followed with a low-dose CT at 6 months...",
  "evidence": [
    {"collection": "imaging_guidelines", "score": 0.89, "text": "ACR Lung-RADS v1.1: Category 3..."},
    {"collection": "imaging_literature", "score": 0.82, "text": "Fleischner Society 2017 guidelines..."}
  ],
  "confidence": 0.87
}
```

### 7.9 Common Questions

**Q: Does the agent interpret actual DICOM images?**
A: Not directly. The RAG engine searches text-based imaging findings, protocols, and guidelines. NVIDIA NIMs (VISTA-3D, VILA-M3) handle pixel-level analysis when deployed. The agent synthesizes clinical knowledge about imaging, not the images themselves.

**Q: How does federated learning work without sharing patient data?**
A: Each hospital trains a local model on its own data. Only the model weight updates (gradients) are sent to a central server for aggregation. No patient images or reports leave the institution.

**Q: Can I add custom imaging protocols?**
A: Yes. Add entries to the imaging_protocols collection via the ingest parser or direct Milvus insertion. The RAG engine will include them in future searches automatically.

### 7.10 Running Your First Query

```bash
# Search for CT lung nodule management guidelines
curl -X POST http://localhost:8524/v1/imaging/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What is the recommended follow-up for a 6mm solid lung nodule?", "top_k": 5}'
```

**Expected output:** A RAG-synthesized response citing ACR Lung-RADS guidelines, recommending 6-12 month follow-up CT depending on risk factors, with references to the NELSON and NLST trials. The response includes evidence passages from `imaging_guidelines` and `imaging_literature` collections with relevance scores.

```bash
# Check available workflows
curl http://localhost:8524/workflows
```

---

## Chapter 8: Precision Oncology Agent

### 8.1 Cancer Is a Disease of the Genome

Every cancer begins with DNA damage. When specific genes acquire mutations, cells lose their normal growth controls and proliferate unchecked. Not all mutations matter -- most are harmless "passengers." The critical few that drive cancer growth are called **driver mutations**, and these are the targets of precision oncology.

**Analogy:** DNA is like a car's instruction manual. A passenger mutation is a typo in the paint color chapter -- harmless. A driver mutation is changing "STOP at red lights" to "GO at red lights" -- dangerous, and the one thing you must fix.

### 8.2 The Molecular Tumor Board

A Molecular Tumor Board (MTB) is a multidisciplinary team that reviews a patient's genomic profile to recommend targeted therapies. The traditional process takes 30-60 minutes per case -- reviewing variant databases, checking clinical trials, cross-referencing guidelines, and debating evidence levels.

The Precision Oncology Agent compresses this to minutes by searching across 11 collections simultaneously, ranking therapies by evidence strength, and generating MTB-ready reports.

### 8.3 Evidence Tiers (AMP/ASCO/CAP)

| Tier | Level | Meaning | Example |
|---|---|---|---|
| IA | Strong | FDA-approved, same tumor type | Vemurafenib for BRAF V600E melanoma |
| IB | Strong | FDA-approved, different tumor type | Pembrolizumab for MSI-H (any solid tumor) |
| IIC | Potential | Clinical trial evidence | ALK inhibitor in ALK+ NSCLC Phase III |
| IID | Potential | Preclinical/case reports | PIK3CA inhibitor in endometrial (early data) |
| III | Unknown | Uncertain significance | VUS in BRCA2 |

### 8.4 Key Actionable Genes

| Gene | Alteration | Matched Therapy | Tumor Types |
|---|---|---|---|
| EGFR | L858R, exon 19 del | Osimertinib | NSCLC |
| BRAF | V600E | Dabrafenib + Trametinib | Melanoma, NSCLC, CRC |
| ALK | Fusion | Alectinib, Lorlatinib | NSCLC |
| KRAS | G12C | Sotorasib, Adagrasib | NSCLC, CRC |
| HER2 | Amplification | Trastuzumab deruxtecan | Breast, gastric |
| NTRK | Fusion | Larotrectinib, Entrectinib | Tissue-agnostic |
| BRCA1/2 | Loss of function | Olaparib, Rucaparib | Ovarian, breast, prostate |
| RET | Fusion | Selpercatinib | NSCLC, thyroid |
| ROS1 | Fusion | Crizotinib, Entrectinib | NSCLC |
| PIK3CA | H1047R, E545K | Alpelisib | Breast |

### 8.5 Resistance Mechanisms

Targeted therapies eventually fail because tumors evolve. The agent tracks
known resistance mechanisms and suggests next-line therapies:

```
  Resistance Example: EGFR
  =========================

  1st line: Osimertinib (targets EGFR L858R)
       |
       | Tumor evolves (6-18 months)
       v
  Resistance: EGFR C797S (binding site mutation)
       |
       v
  2nd line options:
    - Amivantamab + lazertinib (bispecific + TKI)
    - Clinical trial (4th-gen EGFR inhibitor)
    - Combination with MET inhibitor (if MET amplified)
```

### 8.6 The 11 Collections

| #  | Collection                   | Records  | Description                              |
|----|------------------------------|----------|------------------------------------------|
| 1  | oncology_therapies           | 15,200   | Approved targeted/immuno/chemo therapies |
| 2  | oncology_biomarkers          | 8,900    | Actionable genomic biomarkers            |
| 3  | oncology_trials              | 42,000   | Active clinical trials (ClinicalTrials.gov) |
| 4  | oncology_resistance          | 6,400    | Resistance mechanisms and bypass pathways|
| 5  | oncology_guidelines          | 3,800    | NCCN, ESMO, ASCO guidelines             |
| 6  | tumor_profiling              | 22,500   | Tumor mutational burden, MSI, TME data   |
| 7  | oncology_pathways            | 4,100    | Signaling pathway maps (MAPK, PI3K, etc.)|
| 8  | oncology_prognosis           | 7,600    | Stage-specific survival and recurrence   |
| 9  | oncology_combinations        | 5,300    | Combination therapy evidence             |
| 10 | oncology_toxicity            | 3,900    | Treatment-related adverse events         |
| +1 | genomic_evidence (shared)    | 3,560,000| Platform-wide variant vectors            |

### 8.6 How the Agent Processes a Variant

```
Input: BRAF V600E in melanoma
    |
    v
[Entity Detection] → Gene: BRAF | Variant: V600E | Cancer: melanoma
    |
    v
[Evidence Search] → onco_variants (CIViC/OncoKB) → Level IA match
                   → onco_therapies → dabrafenib + trametinib
                   → onco_trials → active Phase III trials
                   → onco_resistance → MAP kinase reactivation
    |
    v
[Therapy Ranking] → Score by: evidence level × biomarker match × trial availability
    |
    v
[MTB Report] → Structured output with tiers, alternatives, resistance warnings
```

### 8.7 Sample Response

```json
{
  "answer": "BRAF V600E is a Level IA actionable target in melanoma...",
  "therapies_ranked": [
    {"rank": 1, "drug": "Dabrafenib + Trametinib", "evidence": "Level IA", "response_rate": "64%"},
    {"rank": 2, "drug": "Vemurafenib + Cobimetinib", "evidence": "Level IA", "response_rate": "68%"},
    {"rank": 3, "drug": "Encorafenib + Binimetinib", "evidence": "Level IA", "response_rate": "63%"}
  ],
  "resistance_mechanisms": ["NRAS Q61K/R secondary mutation", "MEK1/2 amplification"],
  "active_trials": 12,
  "confidence": 0.95
}
```

### 8.8 Common Questions

**Q: What's the difference between Level IA and Level IB evidence?**
A: Level IA means FDA-approved for the same tumor type (e.g., dabrafenib for BRAF V600E melanoma). Level IB means FDA-approved but for a different tumor type (e.g., pembrolizumab for MSI-H in any solid tumor — the approval is tissue-agnostic).

**Q: How does the agent handle variants of uncertain significance (VUS)?**
A: VUS are classified as Level III. The agent surfaces any available functional data, in-silico predictions, and population frequency but explicitly notes the uncertainty. It will not recommend therapy changes based solely on a VUS.

**Q: Can the agent process a full VCF file?**
A: Yes. Upload a VCF through the Streamlit UI or POST to the API. The agent parses all variants, filters for those in known cancer genes, and generates a prioritized analysis of actionable findings.

### 8.9 Running Your First Query

```bash
# Ask about treatment options for a BRAF V600E melanoma
curl -X POST http://localhost:8527/v1/onco/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What are the first-line targeted therapies for BRAF V600E metastatic melanoma?", "top_k": 5}'
```

**Expected output:** A structured MTB-style response citing Level IA evidence for dabrafenib + trametinib (COMBI-d/COMBI-v trials) and vemurafenib + cobimetinib (coBRIM trial), with response rates, PFS data, and resistance mechanism warnings. Evidence passages sourced from `onco_variants`, `onco_therapies`, and `onco_guidelines` collections.

```bash
# Check a specific variant's actionability
curl -X POST http://localhost:8527/v1/onco/query \
  -d '{"question": "Is KRAS G12C actionable in non-small cell lung cancer?"}'
```

---

## Chapter 9: Precision Biomarker Agent

### 9.1 What Are Biomarkers?

A biomarker is any measurable characteristic that indicates a biological state. From simple (blood glucose tells you about diabetes) to complex (a panel of 9 blood tests predicting your biological age).

**Analogy:** Biomarkers are like the dashboard gauges in your car. You don't need to open the hood to know if the engine is overheating -- you check the temperature gauge. Similarly, you don't need a biopsy to know if a patient's liver is stressed -- you check ALT and AST levels.

### 9.2 Biological Age vs Chronological Age

Your birthday tells you how many years you have been alive (chronological age). Your
biology tells a different story. The Precision Biomarker Agent implements the PhenoAge
algorithm (Levine et al., Aging 2018), which calculates biological age from 9 routine
blood biomarkers:

| Biomarker | What It Measures | Direction if Biologically Older |
|-----------|-----------------|-------------------------------|
| Albumin | Liver/nutrition | Lower |
| Creatinine | Kidney function | Higher |
| Glucose | Metabolic health | Higher |
| C-Reactive Protein | Inflammation | Higher |
| Lymphocyte % | Immune function | Lower |
| Mean Cell Volume | Red blood cell size | Higher |
| Red Cell Dist. Width | RBC size variation | Higher |
| Alkaline Phosphatase | Liver/bone health | Higher |
| White Blood Cell Count | Immune activation | Higher |

A 50-year-old with a biological age of 45 has the disease risk profile of a 45-year-old
and is aging well. One with a biological age of 60 has accelerated aging and higher
mortality risk. A 5-year gap (biologically older) is associated with a 1.5x increase
in all-cause mortality risk.

### 9.3 Disease Trajectory Detection

By tracking biomarker trends over time, the agent detects pre-symptomatic disease
across 6 categories. It identifies patients heading toward disease thresholds before
symptoms appear -- enabling preventive intervention.

```
  DISEASE TRAJECTORY DETECTION

  Biomarker Level
       |     Normal Range
       |   +--------------+
       |   |              |  <-- Patient trajectory
       |   |         ....'|......... Alert threshold
       |   |    ....'     |
       |   |...'          |
       |   +--------------+
       |
  -----+------------------------------------------> Time
       Visit 1  Visit 2  Visit 3  Visit 4

  The agent detects the TREND before the value leaves normal range
```

| Disease Category | Key Biomarkers Tracked | Early Signal |
|-----------------|----------------------|-------------|
| Cardiometabolic | HbA1c, fasting glucose, triglycerides | Rising HbA1c trend (5.4 -> 5.6 -> 5.8) |
| Hepatic | ALT, AST, GGT, bilirubin | Rising ALT with normal AST (fatty liver) |
| Renal | eGFR, creatinine, cystatin C | eGFR declining >3 mL/min/year |
| Hematologic | Hemoglobin, MCV, ferritin, B12 | Progressive microcytic anemia pattern |
| Inflammatory | CRP, ESR, IL-6, ferritin | Rising CRP baseline (chronic inflammation) |
| Thyroid | TSH, free T4, free T3 | TSH trending up with normal T4 (subclinical) |

### 9.4 Genotype-Adjusted Reference Ranges

Standard lab reference ranges are derived from population averages, but genetic
variants can shift what is "normal" for an individual. The biomarker agent maintains
a `biomarker_genotype_adjustments` collection that modifies reference ranges.

**Example: PNPLA3 and Liver Function.**

The PNPLA3 I148M variant (rs738409 C>G) is strongly associated with non-alcoholic
fatty liver disease. Carriers have higher baseline ALT levels even without disease:

| PNPLA3 Genotype | Standard ALT Range | Adjusted ALT Range | Clinical Note |
|-----------------|-------------------|-------------------|--------------|
| CC (wild type) | 7-56 U/L | 7-56 U/L | Standard range applies |
| CG (heterozygous) | 7-56 U/L | 7-65 U/L | Slightly elevated baseline expected |
| GG (homozygous I148M) | 7-56 U/L | 7-78 U/L | Significantly elevated baseline; do not over-investigate |

Without genotype adjustment, a GG carrier with an ALT of 70 U/L would be flagged as
abnormal and potentially subjected to unnecessary liver biopsy. With adjustment, the
agent recognizes this as within their genotype-expected range.

### 9.5 The 10 Collections

| # | Collection | Description | Weight |
|---|-----------|-------------|--------|
| 1 | `biomarker_reference` | Reference biomarker definitions and ranges | 0.15 |
| 2 | `biomarker_genetic_variants` | Genetic variants affecting biomarkers | 0.12 |
| 3 | `biomarker_pgx_rules` | Pharmacogenomic dosing rules (CPIC) | 0.10 |
| 4 | `biomarker_disease_trajectories` | Disease progression trajectories | 0.12 |
| 5 | `biomarker_clinical_evidence` | Published clinical evidence | 0.10 |
| 6 | `biomarker_nutrition` | Genotype-aware nutrition guidelines | 0.05 |
| 7 | `biomarker_drug_interactions` | Gene-drug interactions | 0.08 |
| 8 | `biomarker_aging_markers` | Epigenetic aging clock markers | 0.08 |
| 9 | `biomarker_genotype_adjustments` | Genotype-based reference range adjustments | 0.10 |
| 10 | `biomarker_monitoring` | Condition-specific monitoring protocols | 0.10 |
| +1 | `genomic_evidence` | Shared read-only genomic variants | -- |

### 9.6 How Biological Age Is Calculated

```
9 Blood Biomarkers
    |
    v
[PhenoAge Algorithm]
    Albumin (g/dL)         → weight: -0.0336
    Creatinine (mg/dL)     → weight: +0.0095
    Glucose (mg/dL)        → weight: +0.1953
    C-Reactive Protein     → weight: +0.0954
    Lymphocyte %           → weight: -0.0120
    Mean Cell Volume (fL)  → weight: +0.0268
    Red Cell Dist Width    → weight: +0.3306
    Alkaline Phosphatase   → weight: +0.0019
    White Blood Cell Count → weight: +0.0554
    |
    v
[Linear Combination] → Mortality risk score
    |
    v
[Age Conversion] → Biological Age (e.g., 58.3 years for a 55-year-old)
    |
    v
[Delta] → +3.3 years = accelerated aging (inflammation-driven)
```

### 9.7 Sample Response

```json
{
  "biological_age": 58.3,
  "chronological_age": 55,
  "aging_delta": 3.3,
  "aging_status": "accelerated",
  "top_drivers": [
    {"biomarker": "CRP", "value": 2.8, "impact": "primary driver of accelerated aging"},
    {"biomarker": "RDW", "value": 13.5, "impact": "mildly elevated, suggests chronic inflammation"}
  ],
  "disease_trajectories": {
    "cardiovascular": {"risk": "moderate", "trend": "rising"},
    "metabolic": {"risk": "low", "trend": "stable"}
  },
  "recommendations": ["Address CRP elevation (inflammation source workup)", "Recheck in 6 months"]
}
```

### 9.8 Common Questions

**Q: Is biological age clinically validated?**
A: PhenoAge (Levine 2018) has been validated in multiple large cohorts (NHANES, InCHIANTI) as a predictor of all-cause mortality, disease onset, and functional decline independent of chronological age. It uses only routine blood tests — no specialized assays needed.

**Q: Why does genotype affect "normal" lab values?**
A: Genetic variants in metabolic enzymes change baseline biomarker levels. For example, PNPLA3 I148M carriers (common in Hispanic populations) have naturally higher ALT/AST. Without genotype adjustment, these patients may be falsely flagged for liver disease.

**Q: How often should biological age be recalculated?**
A: Every 6-12 months, or after significant lifestyle or treatment changes. The trajectory (improving vs worsening) matters more than any single measurement.

### 9.9 Running Your First Query

```bash
# Calculate biological age from routine labs
curl -X POST http://localhost:8529/v1/biomarker/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "demo-001",
    "age": 55, "sex": "M",
    "biomarkers": {
      "albumin": 4.2, "creatinine": 1.1, "glucose": 105,
      "crp": 2.8, "lymphocyte_pct": 28, "mcv": 88,
      "rdw": 13.5, "alkaline_phosphatase": 72, "wbc": 6.8
    }
  }'
```

**Expected output:** A biological age estimate (e.g., 58.3 years for a 55-year-old — 3.3 years of accelerated aging), aging drivers identified (elevated CRP suggesting chronic inflammation), disease trajectory risk scores across 6 categories, and personalized recommendations for modifiable risk factors.

```bash
# Query biomarker evidence
curl -X POST http://localhost:8529/v1/biomarker/query \
  -d '{"question": "What does elevated sST2 indicate in heart failure prognosis?"}'
```

---

## Chapter 10: CAR-T Intelligence Agent

### 10.1 What Is CAR-T Therapy?

Chimeric Antigen Receptor T-cell therapy takes a patient's own immune cells, engineers them in a lab to recognize cancer, and infuses them back. The "chimeric antigen receptor" is a synthetic protein that gives T-cells a new targeting system.

**Analogy:** Imagine giving soldiers (T-cells) special night-vision goggles (the CAR receptor) that let them see and destroy enemy combatants (cancer cells) that were previously invisible to the immune system.

### 10.2 CAR Protein Structure

```
[scFv] --- [Hinge] --- [TM] --- [Costimulatory] --- [CD3ζ]
  |            |          |           |                  |
  v            v          v           v                  v
Targets    Flexible    Membrane   Sustained         Activation
antigen     spacer     anchor     signaling           signal
```

- **scFv**: Derived from an antibody; recognizes the target antigen (e.g., CD19)
- **Costimulatory domain**: 4-1BB (slower, more persistent T-cells) vs CD28 (faster, stronger initial response)
- **CD3-zeta**: The activation signal that tells the T-cell to kill

### 10.2a Costimulation: 4-1BB vs CD28

The choice of costimulatory domain dramatically affects CAR-T behavior:

| Property            | 4-1BB (CD137)           | CD28                      |
|---------------------|-------------------------|---------------------------|
| Expansion speed     | Slower (days)           | Faster (hours)            |
| Peak T-cell count   | Lower                   | Higher                    |
| Persistence         | Months to years         | Weeks to months           |
| Exhaustion risk     | Lower                   | Higher                    |
| Memory formation    | Central memory (long)   | Effector memory (short)   |
| CRS severity        | Generally milder        | Generally more acute      |
| Best for            | Sustained remission     | Rapid tumor debulking     |
| Example products    | Kymriah, Breyanzi       | Yescarta, Tecartus        |

### 10.2b The 5-Stage Lifecycle

CAR-T therapy is not a pill -- it is a living drug manufactured from a patient's own
cells. Each stage has unique data requirements:

```
  CAR-T THERAPY LIFECYCLE

  Stage 1           Stage 2              Stage 3
  COLLECTION        MANUFACTURING        QUALITY CONTROL
  +---------+       +------------+       +------------+
  | Leuka-  |------>| T-cell     |------>| Potency    |
  | pheresis|       | activation,|       | assays,    |
  | (blood  |       | transduction|      | sterility, |
  | draw)   |       | expansion  |       | identity   |
  +---------+       +------------+       +------------+
                                               |
  Stage 5           Stage 4                    v
  MONITORING        INFUSION             Release
  +---------+       +------------+       criteria met?
  | CRS/    |<------| Lympho-    |<------+
  | ICANS   |       | depleting  |
  | tracking|       | chemo,     |
  | long-   |       | CAR-T      |
  | term    |       | infusion   |
  +---------+       +------------+
```

| Stage | Duration | Key Data Points | Agent Collection |
|-------|----------|-----------------|-----------------|
| Collection | 3-4 hours | CD3+ count, lymphocyte %, viability | cart_biomarkers |
| Manufacturing | 9-14 days | Transduction efficiency, fold expansion, vector copy number | cart_manufacturing |
| Quality Control | 3-5 days | Potency (% lysis), sterility, endotoxin, identity (CD3+CAR+) | cart_assays |
| Infusion | Day 0 | Lymphodepletion regimen, dose, pre-meds | cart_trials |
| Monitoring | Day 0 to years | CRS grade, ICANS grade, B-cell aplasia, response | cart_safety |

### 10.2c CAR Design Trade-offs

The choice of costimulatory domain fundamentally shapes CAR-T behavior:

| Feature | CD28 Costimulation | 4-1BB Costimulation |
|---------|-------------------|-------------------|
| T-cell expansion | Rapid, large peak | Slower, sustained |
| Effector phenotype | More effector memory | More central memory |
| Persistence | Shorter (weeks-months) | Longer (months-years) |
| Exhaustion risk | Higher | Lower |
| CRS severity | Often higher grade | Often lower grade |
| Example product | Yescarta (axi-cel) | Kymriah (tisa-cel) |
| Best for | Aggressive disease needing fast response | Indolent disease needing durability |

### 10.3 The 6 FDA-Approved Products

| Product | Target | Costim | Indication | Year |
|---|---|---|---|---|
| Kymriah | CD19 | 4-1BB | ALL, DLBCL | 2017 |
| Yescarta | CD19 | CD28 | DLBCL, FL | 2017 |
| Tecartus | CD19 | CD28 | MCL | 2020 |
| Breyanzi | CD19 | 4-1BB | DLBCL | 2021 |
| Abecma | BCMA | 4-1BB | Multiple Myeloma | 2021 |
| Carvykti | BCMA | 4-1BB | Multiple Myeloma | 2022 |

### 10.3b The 11 Collections

| # | Collection | Description | Key Fields |
|---|-----------|-------------|------------|
| 1 | `cart_literature` | Published research and patents (5,047 papers) | title, cart_stage, target_antigen |
| 2 | `cart_trials` | ClinicalTrials.gov records (973 records) | phase, target_antigen, car_generation |
| 3 | `cart_constructs` | CAR construct designs and approved products | scfv_origin, costimulatory_domain, vector_type |
| 4 | `cart_assays` | In vitro / in vivo assay results | assay_type, key_metric, metric_value |
| 5 | `cart_manufacturing` | Manufacturing / CMC process records | process_step, parameter, target_spec |
| 6 | `cart_safety` | Pharmacovigilance and post-market safety | event_type, severity_grade, incidence_rate |
| 7 | `cart_biomarkers` | Predictive and pharmacodynamic biomarkers | biomarker_type, clinical_cutoff, evidence_level |
| 8 | `cart_regulatory` | FDA regulatory milestones and approvals | regulatory_event, agency, decision |
| 9 | `cart_sequences` | Molecular/structural data (scFv, binding affinity) | scfv_clone, binding_affinity_kd, framework |
| 10 | `cart_realworld` | Real-world evidence and outcomes | study_type, primary_endpoint, outcome_value |
| 11 | `genomic_evidence` | Read-only genomic variants from Stage 1 | gene, consequence, clinical_significance |

### 10.3c CRS and ICANS Toxicity Grading

| Grade | CRS Symptoms | ICANS Symptoms | Management |
|-------|-------------|----------------|------------|
| 1 | Fever (>38C) | Mild confusion, ICE 7-9 | Supportive care |
| 2 | Fever + hypotension (no vasopressors) | Moderate encephalopathy, ICE 3-6 | Tocilizumab, consider dex |
| 3 | Fever + vasopressors needed | Seizures, ICE 0-2, raised ICP | Tocilizumab + dexamethasone |
| 4 | Life-threatening, ventilator needed | Prolonged seizures, coma, cerebral edema | ICU, high-dose steroids |

### 10.4 Comparative Analysis Mode

The agent uniquely supports "compare X vs Y" queries. When asked "Compare 4-1BB vs CD28 costimulation," it performs dual entity resolution, parallel retrieval across all collections for both entities, and generates a structured side-by-side comparison table -- a capability not found in standard RAG systems.

### 10.5 Toxicity: CRS and ICANS

**Cytokine Release Syndrome (CRS):** When CAR-T cells engage tumor cells
they release a storm of inflammatory cytokines (IL-6, IFN-gamma, TNF-alpha).
This causes fever, hypotension, and in severe cases, organ failure.

| Grade | Symptoms                          | Management                    |
|-------|-----------------------------------|-------------------------------|
| 1     | Fever >= 38 C                     | Supportive care               |
| 2     | Hypotension (no vasopressors)     | Tocilizumab (anti-IL-6)       |
| 3     | Hypotension (vasopressors needed) | Tocilizumab + dexamethasone   |
| 4     | Life-threatening                  | ICU, high-dose steroids       |

**Immune Effector Cell-Associated Neurotoxicity Syndrome (ICANS):** CAR-T
cells can cause neurological symptoms ranging from confusion to seizures.

| Grade | Symptoms                          | Management                    |
|-------|-----------------------------------|-------------------------------|
| 1     | ICE score 7-9 (mild confusion)    | Monitoring, supportive care   |
| 2     | ICE score 3-6 (moderate)          | Dexamethasone 10 mg q6h       |
| 3     | ICE score 0-2, seizure            | High-dose methylprednisolone  |
| 4     | Prolonged seizure, cerebral edema | ICU, anti-epileptics          |

### 10.5b The Agent's Comparative Analysis Pipeline

```
Query: "Compare 4-1BB vs CD28 costimulation"
    |
    v
[Dual Entity Resolution]
    Entity A: 4-1BB (aliases: CD137, TNFRSF9)
    Entity B: CD28 (aliases: T-cell costimulator)
    |
    v
[Parallel Retrieval]
    Entity A search: 24 results across 11 collections
    Entity B search: 22 results across 11 collections
    Total: 46 evidence passages in 365ms
    |
    v
[Structured Comparison]
    | Attribute     | 4-1BB           | CD28              |
    |---------------|-----------------|-------------------|
    | Expansion     | Slower, gradual | Rapid, intense    |
    | Persistence   | Months-years    | Weeks-months      |
    | Exhaustion    | Resistant       | Prone             |
    | Products      | Kymriah,Breyanzi| Yescarta,Tecartus |
    |
    v
[Claude Synthesis] → Narrative with clinical implications
```

### 10.5c Sample Comparative Response

```json
{
  "comparison_type": "costimulatory_domain",
  "entity_a": "4-1BB",
  "entity_b": "CD28",
  "dimensions": {
    "t_cell_expansion": {"4-1BB": "gradual, sustained", "CD28": "rapid, intense peak"},
    "persistence": {"4-1BB": "superior long-term (months-years)", "CD28": "shorter (weeks-months)"},
    "exhaustion_resistance": {"4-1BB": "more resistant via oxidative metabolism", "CD28": "prone via glycolytic shift"},
    "clinical_products": {"4-1BB": "Kymriah, Breyanzi, Abecma, Carvykti", "CD28": "Yescarta, Tecartus"},
    "best_for": {"4-1BB": "indolent/relapsed disease requiring durable response", "CD28": "aggressive/bulky disease requiring rapid cytoreduction"}
  },
  "evidence_count": 46,
  "confidence": 0.91
}
```

### 10.5d Common Questions

**Q: Why can't I just use the "best" costimulatory domain for every patient?**
A: It depends on the clinical scenario. Aggressive bulky disease may benefit from CD28's rapid expansion. Patients needing long-term disease control (e.g., indolent lymphomas) benefit from 4-1BB's persistence. This is why the agent supports comparative analysis — the "best" answer depends on context.

**Q: What is CRS grading and when is it dangerous?**
A: Grade 1 (fever only) is managed with supportive care. Grade 2 (hypotension responding to fluids) may need tocilizumab. Grade 3-4 (vasopressors, ICU) requires tocilizumab ± corticosteroids. Most CRS occurs within 1-14 days post-infusion; earlier onset correlates with higher severity.

**Q: How does the agent track manufacturing data?**
A: The cart_manufacturing collection stores 30+ records covering leukapheresis, transduction efficiency, expansion protocols, cryopreservation, release testing criteria, and lot-specific quality metrics. Queries about manufacturing automatically search this collection with boosted weight.

### 10.6 Running Your First Query

```bash
# Compare costimulatory domains
curl -X POST http://localhost:8521/v1/cart/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Compare 4-1BB vs CD28 costimulatory domains for CAR-T persistence and exhaustion"}'
```

**Expected output:** A structured comparison table showing 4-1BB (slower expansion, longer persistence, less exhaustion, used in Kymriah/Breyanzi/Abecma/Carvykti) vs CD28 (rapid expansion, stronger initial response, faster exhaustion, used in Yescarta/Tecartus). Evidence cited from `cart_literature` (ELIANA, ZUMA-1 trials) and `cart_constructs` collections.

```bash
# Query CAR-T manufacturing
curl -X POST http://localhost:8521/v1/cart/query \
  -d '{"question": "What are the critical quality attributes for CAR-T release testing?"}'
```

---

## Chapter 11: Precision Autoimmune Agent

### 11.1 The Diagnostic Odyssey

The average autoimmune patient waits 4.5 years and sees 4+ doctors before receiving a correct diagnosis. Symptoms overlap across diseases, lab results fluctuate, and no single test is definitive.

**Analogy:** Diagnosing autoimmune disease is like assembling a jigsaw puzzle where the pieces come from different boxes, some pieces are missing, the picture on the box keeps changing, and multiple puzzles may be mixed together.

### 11.2 Autoantibodies as Clues

Autoantibodies are immune proteins that mistakenly attack the body's own tissues. Different patterns suggest different diseases:

| Autoantibody | Associated Disease | Sensitivity |
|---|---|---|
| Anti-dsDNA | Systemic Lupus (SLE) | 70% |
| Anti-CCP | Rheumatoid Arthritis | 95% specific |
| RF (Rheumatoid Factor) | RA (also infections, aging) | 70% |
| Anti-SSA/SSB | Sjogren's Syndrome | 60-70% |
| Anti-Scl-70 | Systemic Sclerosis (diffuse) | 40% |
| Anti-Jo-1 | Dermatomyositis/Polymyositis | 20-30% |
| AChR | Myasthenia Gravis | 85% |

### 11.3 HLA Associations

Human Leukocyte Antigen (HLA) genes determine immune system targeting. Certain HLA alleles dramatically increase disease risk:

| HLA Allele | Disease | Odds Ratio |
|---|---|---|
| HLA-B*27:05 | Ankylosing Spondylitis | 87.4x |
| HLA-DRB1*04:01 | Rheumatoid Arthritis | 4.2x |
| HLA-DQ2/DQ8 | Celiac Disease | 7-10x |
| HLA-B*51:01 | Behcet's Disease | 5.8x |

### 11.4 Disease Activity Scores

| Score | Disease | Remission | Low | Moderate | High |
|---|---|---|---|---|---|
| DAS28-CRP | RA | <2.6 | 2.6-3.2 | 3.2-5.1 | >5.1 |
| SLEDAI-2K | SLE | 0 | 1-5 | 6-10 | >10 |
| CDAI | Crohn's | <150 | 150-219 | 220-450 | >450 |
| BASDAI | AS | <2 | 2-4 | 4-6 | >6 |

### 11.5 The 14 Collections

| # | Collection | Description |
|---|-----------|-------------|
| 1 | `autoimmune_clinical_documents` | Ingested patient records (PDFs) |
| 2 | `autoimmune_patient_labs` | Lab results with flag analysis |
| 3 | `autoimmune_autoantibody_panels` | Autoantibody test result panels |
| 4 | `autoimmune_hla_associations` | HLA allele to disease risk mapping |
| 5 | `autoimmune_disease_criteria` | ACR/EULAR classification criteria |
| 6 | `autoimmune_disease_activity` | Activity scoring reference (DAS28, SLEDAI, etc.) |
| 7 | `autoimmune_flare_patterns` | Flare prediction biomarker patterns |
| 8 | `autoimmune_biologic_therapies` | Biologic drug database with PGx |
| 9 | `autoimmune_pgx_rules` | Pharmacogenomic dosing rules |
| 10 | `autoimmune_clinical_trials` | Autoimmune clinical trials |
| 11 | `autoimmune_literature` | Published literature |
| 12 | `autoimmune_patient_timelines` | Patient diagnostic timelines |
| 13 | `autoimmune_cross_disease` | Cross-disease / overlap syndromes |
| 14 | `genomic_evidence` | Shared read-only genomic variants |

### 11.5a Flare Prediction

Autoimmune flares can cause permanent organ damage if not treated promptly. The
`autoimmune_flare_patterns` collection contains biomarker patterns that predict
flares before symptoms appear:

| Risk Level | Biomarker Pattern | Action |
|-----------|------------------|--------|
| Low | Stable labs, no trend changes | Continue current therapy, routine monitoring |
| Moderate | CRP rising >2x baseline, complement C3/C4 declining | Schedule early follow-up, consider dose adjustment |
| High | Anti-dsDNA doubling, proteinuria appearing, falling C3 | Urgent rheumatology review, consider pulse steroids |
| Critical | Rapid multi-system deterioration, pancytopenia | Emergency admission, IV methylprednisolone |

For SLE, the combination of rising anti-dsDNA + falling complement + increasing
proteinuria has a positive predictive value of ~80% for renal flare within 2-4 weeks.

```
  Flare Prediction: SLE Example
  ===============================

  Months:  -3      -2      -1      Flare
           |       |       |       |
  Anti-dsDNA:  1:80    1:160   1:320   1:640    (rising)
  C3:          95      85      72      58       (falling)
  C4:          22      18      14      9        (falling)
  ESR:         15      22      38      65       (rising)

  Agent alert at Month -1:
  "Anti-dsDNA rising with falling complement. Pattern consistent
   with impending lupus flare. Consider preemptive dose adjustment."
```

### 11.5aa Biologic Therapies

| Drug Class | Mechanism | Key Drugs | Primary Indications |
|-----------|-----------|-----------|-------------------|
| TNF inhibitors | Block TNF-alpha | Adalimumab, etanercept, infliximab | RA, AS, Crohn's, psoriasis |
| Anti-CD20 | Deplete B-cells | Rituximab, ocrelizumab | RA, ANCA vasculitis, MS |
| IL-6R blockers | Block IL-6 signaling | Tocilizumab, sarilumab | RA, GCA, systemic JIA |
| IL-17A inhibitors | Block IL-17A | Secukinumab, ixekizumab | Psoriasis, AS, PsA |
| JAK inhibitors | Block JAK-STAT pathway | Tofacitinib, baricitinib, upadacitinib | RA, PsA, UC, AD |
| Anti-BLyS | Block B-cell survival | Belimumab | SLE |
| CTLA-4 Ig | Block T-cell costimulation | Abatacept | RA |

### 11.5b How the Agent Resolves a Diagnostic Odyssey

```
Patient Timeline: 3 years of symptoms
    |
    v
[Timeline Integration]
    Year 1: Fatigue, joint pain → "anxiety" diagnosis
    Year 2: Malar rash, positive ANA 1:640 → "rosacea" diagnosis
    Year 3: Proteinuria, low C3/C4 → referred to rheumatology
    |
    v
[Pattern Recognition]
    ANA 1:640 homogeneous    → Entry criterion met
    Anti-dsDNA positive      → +6 points (immunology)
    Arthritis (2+ joints)    → +6 points (clinical)
    Malar rash               → +6 points (clinical)
    Low complement C3/C4     → +3 points (immunology)
    Proteinuria              → +4 points (renal)
    Total: 25 points         → SLE diagnosis (threshold: 10)
    |
    v
[Activity Scoring] → SLEDAI-2K = 12 (high activity)
    |
    v
[Therapy Recommendation] → Hydroxychloroquine + mycophenolate (Class III/IV nephritis)
```

### 11.5c Sample Response

```json
{
  "diagnosis": "Systemic Lupus Erythematosus (SLE)",
  "classification_score": 25,
  "threshold": 10,
  "criteria_met": ["ANA entry criterion", "anti-dsDNA", "arthritis", "acute cutaneous lupus", "low complement", "proteinuria"],
  "activity_score": {"sledai_2k": 12, "category": "high_activity"},
  "organ_involvement": ["renal (Class III/IV nephritis suspected)", "musculoskeletal", "dermatologic"],
  "recommendations": [
    "Hydroxychloroquine 200mg BID (background therapy, Class I evidence)",
    "Mycophenolate mofetil for lupus nephritis induction (ALMS trial)",
    "Renal biopsy to confirm nephritis class",
    "Monitor anti-dsDNA titers and complement as activity markers"
  ],
  "confidence": 0.93
}
```

### 11.5d Common Questions

**Q: Why does autoimmune diagnosis take so long?**
A: Three reasons: (1) symptoms overlap across diseases (fatigue, joint pain are common to RA, SLE, Sjogren's, fibromyalgia), (2) autoantibodies can be positive years before clinical disease manifests, and (3) many conditions wax and wane, so symptoms may not be present during a clinic visit. The agent integrates longitudinal data that individual encounters miss.

**Q: What is an overlap syndrome?**
A: When a patient meets criteria for two or more autoimmune diseases simultaneously — e.g., "rhupus" (RA + SLE) or mixed connective tissue disease (MCTD). The agent's cross_disease collection specifically indexes these overlaps.

**Q: How does HLA typing help?**
A: HLA alleles don't diagnose — they stratify risk. HLA-B*27:05 makes ankylosing spondylitis 87x more likely, but ~8% of the general population carries it without disease. The agent uses HLA data as a weighted factor alongside clinical and serologic findings, not as a standalone test.

### 11.6 Running Your First Query

```bash
# Assess a patient with positive ANA and joint symptoms
curl -X POST http://localhost:8532/v1/autoimmune/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Patient has ANA 1:640 homogeneous pattern, anti-dsDNA positive, joint pain, and malar rash. What is the likely diagnosis and recommended workup?"}'
```

**Expected output:** A structured assessment identifying probable Systemic Lupus Erythematosus (SLE) based on 2019 EULAR/ACR classification criteria (ANA entry criterion met, anti-dsDNA +4 points, arthritis +6 points, malar rash +6 points = score well above 10-point threshold). Recommendations include complement levels (C3/C4), CBC, urinalysis for nephritis screening, and SLEDAI-2K activity scoring. Evidence from `autoimmune_disease_criteria` and `autoimmune_autoantibody_panels` collections.

```bash
# Check HLA disease risk
curl -X POST http://localhost:8532/v1/autoimmune/query \
  -d '{"question": "What diseases are associated with HLA-B27 positivity?"}'
```

---

## Chapter 12: Pharmacogenomics Intelligence Agent

### 12.1 Your Genes Affect Your Drugs

Pharmacogenomics (PGx) is the study of how genetic variations affect drug response. The same dose of codeine can be ineffective in one patient (poor CYP2D6 metabolizer -- cannot convert to morphine) and dangerously potent in another (ultra-rapid metabolizer -- produces excessive morphine).

**Analogy:** People process medication like cars process fuel. A sports car (ultra-rapid metabolizer) burns through fuel so fast it runs hot. A diesel truck (normal metabolizer) processes it at a steady rate. A hybrid (intermediate metabolizer) uses less. An electric car (poor metabolizer) can barely use the fuel at all. Same fuel, vastly different outcomes.

### 12.2 Star Alleles and Metabolic Phenotypes

Pharmacogenes are described using **star allele** nomenclature:
- **\*1** = normal function (wild-type)
- **\*2, \*3, etc.** = reduced or no function (variant alleles)
- A patient's two copies (diplotype) determine their phenotype:
  - *1/*1 = Normal metabolizer
  - *1/*2 = Intermediate metabolizer
  - *2/*2 = Poor metabolizer
  - *1/*17 = Ultra-rapid metabolizer (CYP2D6-specific)

### 12.3 The 7 Key Pharmacogenes

| Gene | Enzyme | Drugs Affected | Clinical Impact |
|---|---|---|---|
| CYP2D6 | Debrisoquine hydroxylase | Codeine, tramadol, tamoxifen, SSRIs | Poor metabolizers: no pain relief from codeine |
| CYP2C19 | S-mephenytoin hydroxylase | Clopidogrel, PPIs, voriconazole | Poor metabolizers: clopidogrel ineffective (stent thrombosis) |
| CYP2C9 | Tolbutamide hydroxylase | Warfarin, phenytoin, NSAIDs | Poor metabolizers: warfarin bleeding risk |
| CYP3A5 | Nifedipine oxidase | Tacrolimus, cyclosporine | Expressers need higher tacrolimus doses |
| SLCO1B1 | OATP1B1 transporter | Statins (simvastatin, atorvastatin) | rs4149056 C allele: myopathy risk |
| VKORC1 | Vitamin K epoxide reductase | Warfarin | -1639 G>A: requires lower warfarin dose |
| MTHFR | Methylenetetrahydrofolate reductase | Methotrexate, 5-FU | C677T: toxicity risk, folate supplementation needed |

### 12.3b The 15 Collections

| # | Collection | Description |
|---|-----------|-------------|
| 1 | `pgx_gene_reference` | Pharmacogene star allele definitions and activity scores |
| 2 | `pgx_drug_guidelines` | CPIC/DPWG clinical prescribing guidelines |
| 3 | `pgx_drug_interactions` | Drug-gene interaction records (PharmGKB) |
| 4 | `pgx_hla_hypersensitivity` | HLA-mediated adverse drug reaction screening |
| 5 | `pgx_phenoconversion` | Metabolic phenoconversion via drug-drug interactions |
| 6 | `pgx_dosing_algorithms` | Genotype-guided dosing algorithms and formulas |
| 7 | `pgx_clinical_evidence` | Published PGx clinical evidence and outcomes |
| 8 | `pgx_population_data` | Population-specific allele frequency data |
| 9 | `pgx_clinical_trials` | PGx-related clinical trials |
| 10 | `pgx_fda_labels` | FDA pharmacogenomic labeling information |
| 11 | `pgx_drug_alternatives` | Genotype-guided therapeutic alternatives |
| 12 | `pgx_patient_profiles` | Patient diplotype-phenotype profiles |
| 13 | `pgx_implementation` | Clinical PGx implementation programs |
| 14 | `pgx_education` | PGx educational resources and guidelines |
| 15 | `genomic_evidence` | Shared read-only genomic variants |

### 12.4 Phenoconversion

Phenoconversion occurs when a drug inhibits a metabolic enzyme, effectively changing
a patient's metabolizer phenotype without changing their DNA. This is dangerous and
under-recognized.

**Example: Fluoxetine + Codeine.**

```
  PHENOCONVERSION EXAMPLE

  Genotype:    CYP2D6 *1/*1 (Normal Metabolizer)
  Comedication: Fluoxetine (strong CYP2D6 inhibitor)

  BEFORE fluoxetine:
  Codeine ---[CYP2D6]--> Morphine ----> Pain Relief
                OK                       OK

  AFTER fluoxetine (phenoconversion):
  Codeine ---[CYP2D6]--> Morphine ----> Pain Relief
              BLOCKED                    NONE
              by fluoxetine

  Effective phenotype: Normal -> Poor Metabolizer
  Clinical result: Codeine is INEFFECTIVE for pain
```

The agent's `pgx_phenoconversion` collection catalogs strong, moderate, and weak
inhibitors for each pharmacogene, and flags phenoconversion risk whenever a patient's
medication list includes an enzyme inhibitor alongside a substrate of that enzyme.

### 12.5 HLA Hypersensitivity Screening

Some of the most dangerous adverse drug reactions are mediated by HLA alleles. These
reactions are unpredictable by dose -- they are all-or-nothing immune responses:

| HLA Allele | Drug | Reaction | Severity | Prevalence |
|-----------|------|----------|----------|-----------|
| HLA-B\*57:01 | Abacavir (HIV) | Hypersensitivity syndrome | Fatal if rechallenged | 5-8% European |
| HLA-B\*15:02 | Carbamazepine | Stevens-Johnson Syndrome / TEN | Up to 30% mortality | 8-15% SE Asian |
| HLA-B\*58:01 | Allopurinol | SJS/TEN | Up to 25% mortality | 6-8% SE Asian |
| HLA-A\*31:01 | Carbamazepine | DRESS syndrome | Organ damage | 2-5% European |
| HLA-B\*13:01 | Dapsone | Hypersensitivity syndrome | Multi-organ | 2-6% SE Asian |

**Why screening saves lives.** Before abacavir prescribing was universally preceded by
HLA-B\*57:01 testing, ~5-8% of patients experienced a potentially fatal reaction.
After mandatory screening, the incidence dropped to essentially zero.

### 12.5b How the Agent Translates Genotype to Prescribing

```
Input: CYP2C19 *2/*2 patient prescribed clopidogrel
    |
    v
[Star Allele Interpretation]
    *2 = loss-of-function allele (splicing defect)
    *2/*2 diplotype = Poor Metabolizer
    |
    v
[CPIC Guideline Lookup]
    Drug: clopidogrel
    Gene: CYP2C19
    Phenotype: Poor Metabolizer
    Recommendation: "Use alternative antiplatelet agent"
    Evidence Level: CPIC Level A (strong)
    |
    v
[Alternative Selection]
    ✓ Prasugrel (not CYP2C19 dependent)
    ✓ Ticagrelor (not CYP2C19 dependent)
    ✗ Clopidogrel (AVOID — 3x higher stent thrombosis risk)
    |
    v
[Interaction Check] → No conflicts with current medications
    |
    v
[Report] → "Switch to prasugrel or ticagrelor. Monitor for bleeding."
```

### 12.5c Sample Response

```json
{
  "drug": "clopidogrel",
  "gene": "CYP2C19",
  "diplotype": "*2/*2",
  "phenotype": "Poor Metabolizer",
  "recommendation": "AVOID clopidogrel — use alternative antiplatelet",
  "evidence_level": "CPIC Level A",
  "clinical_impact": "Poor metabolizers have 3x higher risk of stent thrombosis and major adverse cardiac events (TRITON-TIMI 38)",
  "alternatives": [
    {"drug": "prasugrel", "suitability": "preferred", "note": "Not CYP2C19 dependent; contraindicated if age >75, weight <60kg, or prior stroke/TIA"},
    {"drug": "ticagrelor", "suitability": "preferred", "note": "Not CYP2C19 dependent; requires BID dosing; may cause dyspnea"}
  ],
  "monitoring": ["Signs of bleeding", "Platelet function testing if clinical concern"],
  "guideline_source": "CPIC 2022 CYP2C19-Clopidogrel Guideline"
}
```

### 12.5d Common Questions

**Q: Should every patient get PGx testing before starting medications?**
A: Pre-emptive panel testing (testing once, using results for life) is increasingly supported. The cost of a 7-gene panel (~$250-500) is a fraction of an adverse drug reaction hospitalization ($15,000-50,000+). Several health systems (St. Jude, Vanderbilt, Mayo) now implement pre-emptive PGx.

**Q: What is phenoconversion and why does it matter?**
A: A patient genotyped as CYP2D6 normal metabolizer (*1/*1) who takes fluoxetine (a strong CYP2D6 inhibitor) becomes a functional poor metabolizer. If they're also on codeine, it won't be converted to morphine — no pain relief. The agent's phenoconversion calculator detects these drug-drug-gene interactions that static genotyping misses.

**Q: Why aren't all drugs PGx-tested?**
A: Only drugs where genetic variation significantly changes efficacy or toxicity have CPIC guidelines. Many drugs have wide therapeutic windows where genotype doesn't meaningfully change outcomes. The agent focuses on the ~100+ drugs with FDA PGx labeling and CPIC/DPWG guideline support.

### 12.6 Running Your First Query

```bash
# Check drug-gene interaction for clopidogrel
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "clopidogrel", "gene": "CYP2C19", "phenotype": "poor_metabolizer"}'
```

**Expected output:** A CPIC Level A recommendation to avoid clopidogrel in CYP2C19 poor metabolizers due to reduced conversion to active metabolite (increased risk of stent thrombosis and cardiovascular events). Recommends prasugrel or ticagrelor as alternatives. Cites CPIC 2022 guideline, TRITON-TIMI 38 and PLATO trial data. Evidence from `pgx_drug_guidelines` and `pgx_drug_alternatives` collections.

```bash
# Screen for HLA hypersensitivity
curl -X POST http://localhost:8107/v1/pgx/hla-screen \
  -d '{"drug": "abacavir", "hla_allele": "HLA-B*57:01"}'
```

**Expected output:** STOP alert -- HLA-B*57:01 positive patients must NOT receive abacavir. 100% predictive for severe hypersensitivity reaction. Pre-prescription screening is standard of care and FDA black box requirement.

---

## Chapter 13: Cardiology Intelligence Agent

### 13.1 Cardiovascular Disease: The Global #1 Killer

Cardiovascular disease kills 18 million people annually -- more than cancer, infectious diseases, or accidents. Rapid, guideline-concordant cardiac assessment saves lives. But the complexity of modern cardiology (multi-modal imaging, dozens of drugs, evolving guidelines across ACC/AHA/ESC) makes it impossible for any single clinician to keep current.

**Analogy:** The heart is like a house with 4 rooms (chambers), 4 doors (valves), an electrical system (conduction), and plumbing (coronary arteries). Heart failure is when the rooms cannot pump efficiently. Valve disease is a stuck or leaky door. Arrhythmia is a wiring problem. Coronary artery disease is clogged pipes. The Cardiology Agent monitors all four systems simultaneously using genomic, biomarker, and imaging data.

```
  The Heart as a House
  =====================

  Electrical System (SA node, AV node, His-Purkinje)
       |
       v
  +-----+-----+        +-----+-----+
  | Right      |        | Left       |
  | Atrium     |------->| Atrium     |   <- Rooms (4 chambers)
  | (RA)       |  Lung  | (LA)       |
  +-----+------+  circ  +-----+------+
        |                      |
     Tricuspid              Mitral       <- Doors (4 valves)
     Valve                  Valve
        |                      |
  +-----v------+        +-----v------+
  | Right      |        | Left       |
  | Ventricle  |------->| Ventricle  |
  | (RV)       |  Lung  | (LV)       |
  +-----+------+  circ  +-----+------+
        |                      |
     Pulmonic               Aortic
     Valve                  Valve
        |                      |
        v                      v
     Lungs              Body (aorta)     <- Plumbing (arteries)
```

### 13.2 The 6 Risk Calculators

| Calculator | Purpose | Input | Output |
|---|---|---|---|
| ASCVD (PCE) | 10-year heart attack/stroke risk | Age, sex, race, lipids, BP, DM, smoking | 0-100% risk with statin recommendation |
| HEART | Chest pain triage in ED | History, ECG, age, risk factors, troponin | Score 0-10 → low/moderate/high MACE risk |
| CHA2DS2-VASc | Stroke risk in atrial fibrillation | CHF, HTN, age, DM, stroke hx, vascular, sex | Score 0-9 → anticoagulation recommendation |
| HAS-BLED | Bleeding risk on anticoagulation | HTN, renal/liver, stroke, bleeding, age, drugs | Score 0-9 → bleeding risk category |
| MAGGIC | Heart failure mortality | Age, EF, NYHA, BP, BMI, creatinine, comorbidities | 1-year and 3-year mortality % |
| EuroSCORE II | Cardiac surgical mortality | 18 patient/cardiac/operative factors | Predicted operative mortality % |

### 13.3 Heart Failure and GDMT

Guideline-Directed Medical Therapy (GDMT) for heart failure with reduced ejection fraction (HFrEF) has 4 pillars -- each proven to reduce mortality independently:

| Pillar | Drug Class | Key Drug | Mortality Reduction | Key Trial |
|---|---|---|---|---|
| 1 | ARNI/ACEi/ARB | Sacubitril/valsartan | 20% (PARADIGM-HF) | PARADIGM-HF |
| 2 | Beta-blocker | Carvedilol | 35% (COPERNICUS) | COPERNICUS |
| 3 | MRA | Spironolactone | 30% (RALES) | RALES |
| 4 | SGLT2i | Dapagliflozin | 17% (DAPA-HF) | DAPA-HF |

The agent optimizes all 4 pillars simultaneously, checking 14 drugs for contraindications against patient-specific K+, eGFR, BP, and HR, with titration guidance and drug interaction screening.

### 13.4 The 11 Clinical Workflows

CAD Assessment, Heart Failure, Valvular Disease, Arrhythmia, Cardiac MRI, Stress Test, Preventive Risk, Cardio-Oncology, Acute Decompensated HF, Post-MI, and Myocarditis/Pericarditis -- each following the preprocess → execute → postprocess pattern with input validation and cross-modal triggers.

### 13.4b The 13 Collections

| # | Collection | Description | Weight | Est. Records |
|---|-----------|-------------|--------|-------------|
| 1 | `cardio_literature` | Published cardiology research literature | 0.10 | 3,000 |
| 2 | `cardio_trials` | Cardiovascular clinical trials with outcomes | 0.08 | 500 |
| 3 | `cardio_imaging` | Cardiac imaging protocols, findings, criteria | 0.10 | 200 |
| 4 | `cardio_electrophysiology` | ECG/Holter/EP/device electrophysiology | 0.08 | 150 |
| 5 | `cardio_heart_failure` | HF guidelines by type, NYHA, ACC stage | 0.10 | 150 |
| 6 | `cardio_valvular` | Valvular disease severity and interventions | 0.08 | 120 |
| 7 | `cardio_prevention` | CV prevention guidelines and risk factors | 0.10 | 150 |
| 8 | `cardio_interventional` | Interventional/structural procedures | 0.07 | 100 |
| 9 | `cardio_oncology` | Cardio-oncology toxicity monitoring | 0.06 | 100 |
| 10 | `cardio_devices` | AI diagnostic, implantable, wearable devices | 0.04 | 80 |
| 11 | `cardio_guidelines` | ACC/AHA/ESC clinical practice guidelines | 0.10 | 150 |
| 12 | `cardio_hemodynamics` | Invasive/non-invasive hemodynamic parameters | 0.06 | 80 |
| 13 | `genomic_evidence` | Shared read-only genomic variants | -- | 3,560,000 |

Each workflow activates **workflow-specific collection weights** -- for example, the
Heart Failure workflow boosts `cardio_heart_failure` to 0.25 while reducing
`cardio_oncology` to 0.02.

### 13.5 Cross-Modal Genomic Triggers

The agent implements 18 genomic trigger patterns that link cardiac findings to genetic
workup. Key examples:

| Trigger | Criteria | Gene Panel | Urgency |
|---------|----------|-----------|---------|
| Unexplained LVH | Wall thickness >=15mm, no cause | MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1, MYL2, MYL3, GLA, PRKAG2, LAMP2 | High |
| Unexplained DCM | LVEF <45%, dilated, age <60 | TTN, LMNA, RBM20, MYH7, TNNT2, DSP, FLNC, BAG3, SCN5A, PLN | High |
| Long QT | QTc >480ms or >460ms + syncope | KCNQ1, KCNH2, SCN5A, KCNJ2, CALM1-3, TRDN, ANK2 | Critical |
| Brugada pattern | Type 1 Brugada ECG (coved ST) | SCN5A, CACNA1C, CACNB2, SCN1B, SCN2B, SCN3B | Critical |
| CPVT suspected | Exercise-induced polymorphic VT | RYR2, CASQ2, TRDN, CALM1, TECRL | Critical |
| Premature CAD | CAD age <55M/<65F or LDL >=190 | LDLR, PCSK9, APOB, APOE, LDLRAP1 | Moderate |
| Aortic dilation | Root >=4.0cm at age <50 | FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2, MYH11, COL3A1 | High |
| Cardiac amyloid | LVH + diastolic dysfunction + low voltage | TTR | High |
| Arrhythmogenic CM | RV dysfunction + VT + fibro-fatty CMR | PKP2, DSP, DSG2, DSC2, JUP, TMEM43, PLN, FLNC | High |

```
  CROSS-MODAL TRIGGER: IMAGING -> GENOMICS

  CMR Finding: Unexplained LVH (18mm septal thickness)
       |
       v
  Trigger: "unexplained_lvh" activated
       |
       v
  Gene Panel: MYH7, MYBPC3, TNNT2, TNNI3, TPM1...
       |
       v
  Query genomic_evidence (3.5M variants)
       |
       v
  Result: MYBPC3 c.1504C>T (p.Arg502Trp)
          Pathogenic | ClinVar | HCM-associated
       |
       v
  Clinical Action: Confirm HCM diagnosis, family cascade screening,
                   exercise restriction counseling, SCD risk assessment
```

### 13.5aa Cardiac Imaging Integration

The Cardiology Agent integrates measurements from four cardiac imaging modalities:

**Echocardiography (Echo):**

| Parameter | Normal Range | Clinical Significance |
|-----------|-------------|----------------------|
| LVEF | 52-72% (M), 54-74% (F) | <40% = HFrEF, 40-49% = HFmrEF, >=50% = HFpEF |
| LV GLS | -18% to -22% | >-16% suggests subclinical dysfunction |
| LAVI | <34 mL/m2 | >34 = diastolic dysfunction marker |
| TAPSE | >=17mm | <17mm = RV systolic dysfunction |
| E/e' ratio | <14 | >=14 = elevated filling pressures |

**CT Coronary Calcium Scoring:**

| Agatston Score | Risk Category | 10-year Event Rate |
|---------------|--------------|-------------------|
| 0 | Very low risk | <1% |
| 1-99 | Low risk | ~5% |
| 100-399 | Moderate risk | ~10% |
| >=400 | High risk | ~15-20% |

### 13.5b How the GDMT Optimizer Works

```
Input: LVEF 28%, NYHA III, on metoprolol 50mg only
    |
    v
[EF Classification] → HFrEF (LVEF ≤40%)
    |
    v
[4-Pillar Gap Analysis]
    Pillar 1 (ARNI/ACEi/ARB): ✗ NOT STARTED
    Pillar 2 (Beta-blocker):  ▲ SUB-TARGET (50mg, target 200mg)
    Pillar 3 (MRA):           ✗ NOT STARTED
    Pillar 4 (SGLT2i):        ✗ NOT STARTED
    |
    v
[Contraindication Check]
    SBP 108 → Cautious ARNI initiation (start 24/26mg BID)
    K+ 4.1  → MRA safe (spironolactone 25mg)
    eGFR 58 → SGLT2i safe (dapagliflozin 10mg)
    HR 72   → Beta-blocker uptitration safe
    |
    v
[Drug Interaction Screen] → No conflicts detected
    |
    v
[Titration Plan]
    1. Start sacubitril/valsartan 24/26mg BID (check BP in 1 week)
    2. Start spironolactone 25mg daily (check K+ and creatinine in 1 week)
    3. Start dapagliflozin 10mg daily (no titration needed)
    4. Uptitrate metoprolol: 50mg → 100mg → 200mg (every 2 weeks)
```

### 13.5c Sample GDMT Response

```json
{
  "hf_phenotype": "HFrEF",
  "lvef": 28,
  "nyha_class": "III",
  "four_pillars_status": {
    "arni_acei_arb": "NOT STARTED — recommend sacubitril/valsartan 24/26mg BID",
    "beta_blocker": "SUB-TARGET — metoprolol 50mg, target 200mg daily",
    "mra": "NOT STARTED — recommend spironolactone 25mg daily",
    "sglt2i": "NOT STARTED — recommend dapagliflozin 10mg daily"
  },
  "recommendations": [
    "Initiate sacubitril/valsartan 24/26mg BID (SBP 108 — start low, uptitrate in 2 weeks)",
    "Initiate spironolactone 25mg daily — check K+ and creatinine in 1 week",
    "Initiate dapagliflozin 10mg daily — no titration needed (DAPA-HF: 26% reduction in HF events)",
    "Uptitrate metoprolol succinate: 50mg → 100mg → 200mg (every 2 weeks, target HR 60-70)"
  ],
  "monitoring_plan": "Recheck BMP (K+, creatinine) in 1 week after MRA initiation, then monthly x3",
  "device_assessment": "ICD indicated if LVEF remains ≤35% after 90 days of optimized GDMT (SCD-HeFT)"
}
```

### 13.5d Common Questions

**Q: Why does the agent prioritize all 4 pillars simultaneously?**
A: The 2022 AHA/ACC/HFSA guideline recommends initiating all 4 pillars as rapidly as tolerated, rather than the traditional sequential approach. The STRONG-HF trial (2022) showed that rapid uptitration within 2 weeks of discharge reduced HF death/readmission by 34%. Each pillar has independent mortality benefit — delaying any one delays survival improvement.

**Q: What happens if a patient can't tolerate a pillar?**
A: The agent documents intolerance reasons and suggests alternatives within the same class (e.g., eplerenone if spironolactone causes gynecomastia) or marks the pillar as contraindicated with the specific reason. It never leaves a gap without explanation.

**Q: How do the 18 cross-modal genomic triggers work in practice?**
A: When a workflow detects a clinical finding meeting predefined criteria (e.g., unexplained LVH >=15mm on echo), the agent immediately generates a CrossModalTrigger object specifying the gene panel (e.g., HCM panel: MYH7, MYBPC3, TNNT2 -- 11 genes), urgency level, estimated cost ($1,500), and turnaround time (21 days). This appears in the workflow output and can trigger a genomic_evidence collection search for known variants.

### 13.6 Running Your First Query

```bash
# Calculate ASCVD 10-year risk
curl -X POST http://localhost:8126/v1/cardio/risk/ascvd \
  -H "Content-Type: application/json" \
  -d '{
    "age": 55, "sex": "male", "race": "white",
    "total_cholesterol": 240, "hdl_cholesterol": 42,
    "systolic_bp": 145, "bp_treatment": true,
    "diabetes": false, "smoker": false
  }'
```

**Expected output:** A 10-year ASCVD risk percentage (e.g., 14.2% — intermediate risk), with guideline-concordant recommendations: moderate-intensity statin, consider coronary artery calcium score for shared decision-making, optimize blood pressure. Risk category, interpretation, and 2019 ACC/AHA Primary Prevention Guideline citation included.

```bash
# Optimize heart failure GDMT
curl -X POST http://localhost:8126/v1/cardio/gdmt/optimize \
  -H "Content-Type: application/json" \
  -d '{
    "lvef": 28, "nyha_class": "III",
    "current_medications": [{"name": "metoprolol succinate", "dose": "50mg daily"}],
    "patient_data": {"systolic_bp": 108, "heart_rate": 72, "potassium": 4.1, "creatinine": 1.2, "egfr": 58}
  }'
```

**Expected output:** HFrEF phenotype identified. 4-pillar analysis: beta-blocker initiated (uptitrate metoprolol to 200mg target), ARNI not started (recommend sacubitril/valsartan 24/26mg BID), MRA not started (recommend spironolactone 25mg, check K+ in 1 week), SGLT2i not started (recommend dapagliflozin 10mg). Contraindication check: SBP 108 — cautious with ARNI initiation, start low. Drug interaction screening: no conflicts detected.

```bash
# Ask a clinical question
curl -X POST http://localhost:8126/v1/cardio/query \
  -d '{"question": "When should I consider CRT implantation in a heart failure patient?"}'
```

---

## Chapter 14: Putting It All Together

### 14.1 A Complete Patient Journey

1. **Patient DNA arrives** as FASTQ files (200 GB raw sequencing data)
2. **Stage 1 (Genomics, 120 min):** BWA-MEM2 aligns reads → DeepVariant calls variants → 11.7M variants in VCF
3. **Annotation:** ClinVar matches (35,616), AlphaMissense predictions (6,831), VEP consequences → 3.5M searchable vectors in genomic_evidence
4. **Stage 2 (Clinical Intelligence):** Clinician queries "What therapeutic targets exist for this patient?"
5. **Oncology Agent** identifies BRAF V600E (Level IA evidence) → recommends dabrafenib + trametinib
6. **Biomarker Agent** calculates biological age, detects cardiovascular trajectory risk
7. **Pharmacogenomics Agent** identifies CYP2D6 poor metabolizer → flags codeine contraindication
8. **Cardiology Agent** checks cardiac safety of proposed therapies → no LVEF concerns
9. **Stage 3 (Drug Discovery, 5 min):** MolMIM generates 100+ BRAF inhibitor analogues → DiffDock scores binding → Top 10 candidates ranked by QED
10. **Report generated:** PDF with ranked drug candidates ready for medicinal chemistry

Total time: **<5 hours** on a single $4,699 DGX Spark.

### 14.1b Cross-Agent Coordination

The agents do not operate in isolation. They coordinate through three mechanisms:

**1. Shared genomic_evidence collection.** All agents have read-only access to the
same 3.56 million variant records created by Stage 1. This is the common data
substrate that enables cross-modal reasoning.

**2. Server-Sent Events (SSE).** Each agent publishes events when significant findings
are detected. Other agents subscribe to relevant event streams:

```
  EVENT PUBLISHING ARCHITECTURE

  Imaging Agent                     Oncology Agent
  +-----------+    SSE Event:       +-----------+
  | Lung      |--->"lung_nodule    | Variant   |
  | nodule    |    detected,       | matcher   |
  | detected  |    Lung-RADS 4B"   | activates |
  +-----------+         |          +-----------+
                        |
                        v
  PGx Agent         Cardiology Agent
  +-----------+     +-----------+
  | Checks    |     | Pre-chemo |
  | drug-gene |     | cardiac   |
  | before Tx |     | baseline  |
  +-----------+     +-----------+
```

**3. Patient 360 Dashboard.** The biomarker agent's Patient 360 view aggregates
findings from all agents into a unified patient summary.

### 14.1c Collection Inventory

All collections across all 11 intelligence agents:

| Agent | Collections | Total |
|-------|------------|-------|
| Imaging | imaging_literature, imaging_trials, imaging_findings, imaging_protocols, imaging_devices, imaging_anatomy, imaging_benchmarks, imaging_guidelines, imaging_report_templates, imaging_datasets | 10 |
| Oncology | onco_literature, onco_trials, onco_variants, onco_biomarkers, onco_therapies, onco_pathways, onco_guidelines, onco_resistance, onco_outcomes, onco_cases | 10 |
| Biomarker | biomarker_reference, biomarker_genetic_variants, biomarker_pgx_rules, biomarker_disease_trajectories, biomarker_clinical_evidence, biomarker_nutrition, biomarker_drug_interactions, biomarker_aging_markers, biomarker_genotype_adjustments, biomarker_monitoring | 10 |
| CAR-T | cart_literature, cart_trials, cart_constructs, cart_assays, cart_manufacturing, cart_safety, cart_biomarkers, cart_regulatory, cart_sequences, cart_realworld | 10 |
| Autoimmune | autoimmune_clinical_documents through autoimmune_cross_disease | 13 |
| PGx | pgx_gene_reference through pgx_education | 14 |
| Cardiology | cardio_literature through cardio_hemodynamics | 12 |
| Shared | genomic_evidence | 1 |
| **Total** | | **80** |

### 14.2 The Technology Stack

| Component | Technology | Purpose |
|---|---|---|
| GPU Compute | NVIDIA DGX Spark (GB10) | Hardware acceleration for all stages |
| Genomics | Parabricks 4.6 + DeepVariant | GPU-accelerated variant calling |
| Vector Database | Milvus 2.4 | Semantic search across all collections |
| Embeddings | BGE-small-en-v1.5 (384-dim) | Text-to-vector conversion |
| LLM | Claude Sonnet 4.6 (Anthropic) | Natural language synthesis with citations |
| Drug Generation | BioNeMo MolMIM | Generative molecular design |
| Molecular Docking | DiffDock | Binding pose prediction |
| Chemistry | RDKit | Drug-likeness scoring, SMILES parsing |
| Web UI | Streamlit | Interactive clinical interfaces |
| API | FastAPI | REST endpoints for programmatic access |
| Monitoring | Prometheus + Grafana | Observability and dashboards |
| Orchestration | Docker Compose + Nextflow | Multi-service deployment and pipeline management |

### 14.3 Complete Port Map

| Port | Service | Agent/Pipeline |
|---|---|---|
| 8080 | Landing Page / Hub | Platform |
| 5000 | Genomic Foundation Engine | Stage 1 |
| 5001 | RAG/Chat API | Stage 2 |
| 8501 | RAG Chat Interface | Stage 2 |
| 8505 | Drug Discovery | Stage 3 |
| 8524/8525 | Imaging Agent (API/UI) | Intelligence |
| 8526/8527 | Oncology Agent | Intelligence |
| 8528/8529 | Biomarker Agent | Intelligence |
| 8521 | CAR-T Agent | Intelligence |
| 8531/8532 | Autoimmune Agent | Intelligence |
| 8507/8107 | Pharmacogenomics Agent | Intelligence |
| 8126/8536 | Cardiology Agent (API/UI) | Intelligence |
| 19530 | Milvus (shared) | Infrastructure |
| 9099 | Prometheus | Monitoring |
| 3000 | Grafana | Monitoring |

### 14.4 Getting Started

```bash
# Start the full platform
docker compose -f docker-compose.dgx-spark.yml up -d

# Open the landing page
open http://localhost:8080

# Check health of all services
curl http://localhost:8080/api/check-services
```

---

## Glossary

| Term | Definition |
|---|---|
| ASCVD | Atherosclerotic Cardiovascular Disease -- heart attacks and strokes caused by plaque buildup |
| BAM | Binary Alignment Map -- compressed file format for aligned sequencing reads |
| BGE | BAAI General Embedding -- the sentence transformer model used for vector embeddings |
| BioNeMo | NVIDIA's platform for biomolecular AI models |
| CAR-T | Chimeric Antigen Receptor T-cell therapy -- engineered immune cells for cancer |
| CDS | Clinical Decision Support -- systems that aid clinical judgment |
| ClinVar | NCBI database of clinically annotated genomic variants |
| COSINE | Cosine similarity -- metric for comparing vector similarity (0-1 scale) |
| CPIC | Clinical Pharmacogenetics Implementation Consortium |
| CRS | Cytokine Release Syndrome -- immune overactivation from CAR-T therapy |
| CTRCD | Cancer Therapy-Related Cardiac Dysfunction |
| DiffDock | Diffusion-based molecular docking model |
| DGX Spark | NVIDIA's $4,699 GPU workstation for AI |
| DPWG | Dutch Pharmacogenetics Working Group |
| FASTQ | Raw sequencing data format (reads + quality scores) |
| FHIR | Fast Healthcare Interoperability Resources -- standard for health data exchange |
| GDMT | Guideline-Directed Medical Therapy -- evidence-based HF treatment |
| GLS | Global Longitudinal Strain -- cardiac deformation measure |
| GPU | Graphics Processing Unit -- parallel processor for AI workloads |
| GRCh38 | Genome Reference Consortium Human Build 38 -- current reference genome |
| HLA | Human Leukocyte Antigen -- immune system genes |
| ICANS | Immune Effector Cell-Associated Neurotoxicity Syndrome |
| IVF_FLAT | Index type for Milvus vector search (inverted file with flat vectors) |
| LGE | Late Gadolinium Enhancement -- cardiac MRI technique for scar detection |
| LLM | Large Language Model -- AI for natural language understanding |
| LVEF | Left Ventricular Ejection Fraction -- measure of heart pumping strength |
| Milvus | Open-source vector database for similarity search |
| MolMIM | Masked molecular modeling -- generative AI for drug design |
| MTB | Molecular Tumor Board -- multidisciplinary cancer treatment planning |
| NIM | NVIDIA Inference Microservice -- containerized AI model |
| NYHA | New York Heart Association -- heart failure symptom classification |
| PDB | Protein Data Bank -- repository of 3D protein structures |
| PGx | Pharmacogenomics -- study of how genes affect drug response |
| QED | Quantitative Estimate of Drug-likeness (0-1 scale) |
| RAG | Retrieval-Augmented Generation -- combining search with LLM synthesis |
| SMILES | Simplified Molecular Input Line Entry System -- text representation of molecules |
| SNV | Single Nucleotide Variant -- a single DNA base change |
| SSE | Server-Sent Events -- real-time event streaming protocol |
| VCF | Variant Call Format -- standard file format for genomic variants |
| Vector | A list of numbers (384 dimensions) representing the meaning of text |
| VISTA-3D | NVIDIA's 3D medical image segmentation model (132 classes) |
| ANA | Antinuclear Antibody -- screening test for autoimmune diseases; patterns (homogeneous, speckled, nucleolar, centromere) suggest specific conditions |
| Costimulatory Domain | The signaling region in a CAR protein (4-1BB or CD28) that determines T-cell persistence and expansion characteristics |
| Diagnostic Odyssey | The multi-year, multi-doctor journey autoimmune patients endure before receiving a correct diagnosis (average 4.5 years) |
| FLARE (NVIDIA) | Federated Learning Application Runtime Environment -- privacy-preserving multi-site AI model training |
| Lake Louise Criteria | CMR diagnostic criteria for myocarditis requiring 2 of 3: T2 edema, T1/ECV elevation, non-ischemic LGE pattern |
| Lipinski Rule of Five | Drug-likeness filter: MW <500, logP <5, HBD <5, HBA <10 -- molecules violating multiple rules are unlikely to be orally bioavailable |
| MAISI | NVIDIA NIM for synthetic medical image generation |
| Parabricks | NVIDIA's GPU-accelerated genomics toolkit (BWA-MEM2, DeepVariant, samtools) |
| Phenoconversion | When a drug inhibits a metabolic enzyme, changing a patient's effective metabolic phenotype without changing their DNA |
| Prometheus | Open-source monitoring system for metrics collection and alerting |
| scFv | Single-chain variable fragment -- the antigen-binding portion of a CAR protein, derived from an antibody |
| Star Allele | Nomenclature for pharmacogene variants (\*1 = normal function, \*2/\*3 = reduced/absent function) |
| TPSA | Topological Polar Surface Area -- molecular descriptor predicting oral absorption and blood-brain barrier penetration |
| VILA-M3 | NVIDIA vision-language model for medical image interpretation |
| ACEi | Angiotensin-Converting Enzyme inhibitor. Blood pressure drug that blocks angiotensin II production |
| ADC | Apparent Diffusion Coefficient. MRI measure of water molecule movement in tissue |
| ARNI | Angiotensin Receptor-Neprilysin Inhibitor. HF drug class (e.g., sacubitril-valsartan) |
| ARVC | Arrhythmogenic Right Ventricular Cardiomyopathy. Genetic heart muscle disease |
| BASDAI | Bath Ankylosing Spondylitis Disease Activity Index |
| BCMA | B-Cell Maturation Antigen. CAR-T target in multiple myeloma |
| DAS28 | Disease Activity Score using 28 joints. RA activity measurement |
| ECG | Electrocardiogram. Recording of heart electrical activity |
| FLARE | Federated Learning Application Runtime Environment. NVIDIA privacy-preserving multi-site AI training |
| HCM | Hypertrophic Cardiomyopathy. Genetic thickening of heart muscle |
| ICD | Implantable Cardioverter-Defibrillator. Device that shocks fatal arrhythmias |
| Lake Louise | CMR diagnostic criteria for myocarditis (T2 edema, T1/ECV, non-ischemic LGE) |
| Lipinski | Rule of Five. Drug-likeness filter: MW<500, logP<5, HBD<5, HBA<10 |
| MRA | Mineralocorticoid Receptor Antagonist. HF drug class (spironolactone, eplerenone) |
| Parabricks | NVIDIA GPU-accelerated genomics toolkit (BWA-MEM2, DeepVariant, samtools) |
| PCI | Percutaneous Coronary Intervention. Catheter-based coronary stent procedure |
| PCE | Pooled Cohort Equations. ASCVD 10-year risk calculator algorithm |
| Phenoconversion | When a drug inhibits a metabolic enzyme, changing functional phenotype without changing DNA |
| scFv | Single-chain variable fragment. Antibody-derived targeting domain in CAR constructs |
| SGLT2i | Sodium-Glucose Co-Transporter 2 inhibitor. HF/diabetes drug (dapagliflozin, empagliflozin) |
| SLEDAI-2K | SLE Disease Activity Index 2000. Lupus activity measurement |
| TAVR | Transcatheter Aortic Valve Replacement. Minimally invasive valve procedure |
| VUS | Variant of Uncertain Significance. Genomic variant with insufficient evidence |

---

## Quick Reference: Complete Collection Inventory

### Summary by Agent

| Agent             | Owned Collections | + Shared | Total |
|-------------------|-------------------|----------|-------|
| Imaging           | 10                | 1        | 11    |
| Oncology          | 10                | 1        | 11    |
| Biomarker         | 10                | 1        | 11    |
| CAR-T             | 11                | 1        | 12    |
| Autoimmune        | 14                | 1        | 15    |
| Pharmacogenomics  | 15                | 1        | 16    |
| Cardiology        | 13                | 1        | 14    |
| **Total**         | **83**            | **1 shared** | **84 unique** |

### Full Collection List (alphabetical by agent)

```
AUTOIMMUNE (14 owned)
  autoantibody_profiles          autoimmune_biomarkers
  autoimmune_diseases            autoimmune_guidelines
  biologic_therapies             disease_activity_scores
  environmental_triggers         flare_patterns
  hla_associations               immunogenetics
  overlap_syndromes              pediatric_autoimmune
  pregnancy_autoimmune           treatment_sequencing

BIOMARKER (10 owned)
  assay_methods                  biomarker_guidelines
  biomarker_interactions         biomarker_panels
  biomarker_reference            genotype_adjustments
  pharmacodynamic_markers        phenoage_models
  screening_protocols            trajectory_patterns

CARDIOLOGY (13 owned)
  cardiac_biomarkers             cardiac_devices
  cardiac_electrophysiology      cardiac_genes
  cardiac_guidelines             cardiac_imaging
  cardiac_interventions          cardiac_medications
  cardiac_prevention             cardiac_rehabilitation
  cardiac_risk_models            cardiac_trials
  cardio_oncology

CAR-T (11 owned)
  cart_biomarkers                cart_combinations
  cart_constructs                cart_manufacturing
  cart_next_gen                  cart_outcomes
  cart_products                  cart_resistance
  cart_targets                   cart_toxicity
  cart_trials

IMAGING (10 owned)
  ai_model_registry              anatomy_atlas
  contrast_protocols             imaging_biomarkers
  imaging_guidelines             imaging_protocols
  pathology_patterns             radiation_dose
  radiogenomics                  radiology_findings

ONCOLOGY (10 owned)
  oncology_biomarkers            oncology_combinations
  oncology_guidelines            oncology_pathways
  oncology_prognosis             oncology_resistance
  oncology_therapies             oncology_toxicity
  oncology_trials                tumor_profiling

PHARMACOGENOMICS (15 owned)
  pgx_cardiology                 pgx_clinical_annotations
  pgx_cpic_guidelines            pgx_diplotypes
  pgx_dosing                     pgx_dpwg_guidelines
  pgx_gene_drug                  pgx_hla_hypersensitivity
  pgx_implementation             pgx_interactions
  pgx_oncology                   pgx_pathways
  pgx_pediatric                  pgx_population_frequencies
  pgx_psychiatry

SHARED (1)
  genomic_evidence               (3,560,000 vectors, read by all)
```

---

## Review Questions (Part 2)

### Chapter 7: Imaging

24. Name the four imaging modalities and what each is best at detecting.
25. What does VISTA-3D do, and how many anatomical classes does it segment?
26. How does federated learning (NVIDIA FLARE) protect patient privacy?

### Chapter 8: Oncology

27. What is the difference between a driver mutation and a passenger mutation?
28. What does AMP/ASCO/CAP Level 1A evidence mean?
29. Name three actionable genes and their corresponding targeted therapies.

### Chapter 9: Biomarker

30. What is PhenoAge, and how many blood biomarkers does it use?
31. What is a genotype-adjusted reference range, and why does it matter?
32. Name the six disease trajectory categories the agent monitors.

### Chapter 10: CAR-T

33. What are the five domains of a CAR protein, from outside to inside?
34. What is the key trade-off between 4-1BB and CD28 costimulation?
35. What is CRS, and what is the first-line treatment for Grade 2?

### Chapter 11: Autoimmune

36. How long is the average diagnostic odyssey for autoimmune disease?
37. What is the strongest HLA-disease association in the platform, and what
    is its odds ratio?
38. How does the agent predict flares before they occur?

### Chapter 12: Pharmacogenomics

39. What is a star allele, and what does \*1 represent?
40. Explain phenoconversion with a clinical example.
41. Why must HLA-B\*57:01 be tested before prescribing abacavir?

### Chapter 13: Cardiology

42. Name the four pillars of GDMT for heart failure.
43. What is the CHA2DS2-VASc score used for?
44. Name three genomic triggers that activate automatic cardiac assessment.

### Chapter 14: Integration

45. How long does the complete patient journey take, from DNA to drug
    candidates?
46. What two mechanisms do agents use to coordinate with each other?
47. How many total Milvus collections does the platform maintain?

---

*HCLS AI Factory Learning Guide Foundations -- Unified Edition*
*Version 1.0.0 | March 2026 | Adam Jones*
*Apache 2.0 License | NVIDIA DGX Spark Platform*

---

!!! warning "Clinical Decision Support Disclaimer"
    The HCLS AI Factory platform and all intelligence agents described in this document are clinical decision support research tools. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
