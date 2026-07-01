# HCLS AI Factory — Pediatric Oncology Demo Guide (GUI Walkthrough)

> **Version:** 1.0 | **Date:** 2026-03-23
> **Platform:** NVIDIA DGX Spark | **Runtime:** Single node, ~1.1 TB footprint
> **Total Demo Duration:** 75–90 minutes (all 6 demos) or run individually
> **Audience:** Clinicians, researchers, hospital IT, NVIDIA partners, investors

---

## Table of Contents

1. [Pre-Demo Setup & Checklist](#pre-demo-setup--checklist)
2. [Demo 1: End-to-End Precision Medicine Pipeline — Evelyn, 8yo, B-ALL](#demo-1-end-to-end-precision-medicine-pipeline)
3. [Demo 2: Multi-Agent Tumor Board — Evelyn, Day 29, MRD+](#demo-2-multi-agent-tumor-board)
4. [Demo 3: Cardiotoxicity Prevention — Marcus, 6yo, Neuroblastoma](#demo-3-cardiotoxicity-prevention)
5. [Demo 4: Rare Disease + Cancer Predisposition — Aurora, 4yo, Retinoblastoma](#demo-4-rare-disease--cancer-predisposition)
6. [Demo 5: CAR-T Decision — Ethan, 12yo, Relapsed ALL](#demo-5-car-t-decision)
7. [Demo 6: Medulloblastoma + Drug Discovery — Aiden, 10yo, SHH Subtype](#demo-6-medulloblastoma--drug-discovery)
8. [Post-Demo: Closing Statement](#post-demo-closing-statement)
9. [Appendix A: Agent Port Reference](#appendix-a-agent-port-reference)
10. [Appendix B: Troubleshooting](#appendix-b-troubleshooting)
11. [Appendix C: Presenter Quick Reference Card](#appendix-c-presenter-quick-reference-card)

---

## Pre-Demo Setup & Checklist

### 30 Minutes Before Demo

**1. Start all services:**

```bash
cd /home/adam/projects/hcls-ai-factory
docker compose -f docker-compose.dgx-spark.yml up -d
```

**2. Verify infrastructure is healthy:**

```bash
bash health-monitor.sh --check-all
```

**3. Open the Landing Page and confirm agent status:**

- Navigate to **http://localhost:8080**
- Verify: All 11 agent tiles show **green** status indicators
- Verify: Three engine cards (Genomics, RAG/Chat, Drug Discovery) are visible
- If any agent shows red/yellow, restart it individually before proceeding

**4. Prepare your browser:**

| Tab # | URL | Label for Tab |
|-------|-----|---------------|
| 1 | http://localhost:8080 | Landing Page |
| 2 | http://localhost:5000 | Genomics Portal |
| 3 | http://localhost:8501 | RAG Chat |
| 4 | http://localhost:8526 | Oncology |
| 5 | http://localhost:8510 | Drug Discovery |
| 6 | http://localhost:8528 | Biomarker |
| 7 | http://localhost:8521 | CAR-T |
| 8 | http://localhost:8130 | Single-Cell |
| 9 | http://localhost:8128 | Clinical Trial |
| 10 | http://localhost:8507 | PGx |
| 11 | http://localhost:8536 | Cardiology |
| 12 | http://localhost:8529 | Neurology |
| 13 | http://localhost:8531 | Autoimmune |
| 14 | http://localhost:8544 | Rare Disease |
| 15 | http://localhost:8525 | Imaging |

> **PRESENTER NOTE:** Pre-load all tabs before the audience arrives. Each Streamlit agent
> takes 3-5 seconds on first load. You do not want loading spinners during the live demo.

**5. Display setup:**

- Use a single external monitor or projector at 1920x1080 minimum
- Set browser zoom to 90-100% so Streamlit sidebar and main panel are both visible
- Disable browser notifications and OS popups
- Close Slack, email, and any other applications that might produce alerts

**6. Final pre-flight check:**

- [ ] All 15 browser tabs loaded without errors
- [ ] Landing Page shows all agents green
- [ ] RAG Chat loads with variant stats in sidebar
- [ ] Oncology agent shows all 5 tabs
- [ ] Drug Discovery portal is responsive
- [ ] Your speaker notes (this document) are on a second screen or printed

---

## Demo 1: "DNA to Drug" — End-to-End Precision Medicine Pipeline

### Patient: Evelyn Chen, 8 years old, Female, B-cell Acute Lymphoblastic Leukemia

**Duration:** 15-20 minutes
**Story arc:** Follow one child's DNA from sequencing through diagnosis, molecular classification, therapy ranking, trial matching, and drug candidate generation.

**Browser tabs needed:**
- Landing Page (http://localhost:8080)
- Genomics Portal (http://localhost:5000)
- RAG Chat (http://localhost:8501)
- Oncology (http://localhost:8526)
- Drug Discovery (http://localhost:8510)

### Foundations Learning — End-to-End Precision Medicine Pipeline

**Precision medicine** is a paradigm shift from the traditional "one-size-fits-all" approach to treatment. Instead of prescribing the same therapy to every patient with the same diagnosis, precision medicine uses an individual's molecular profile — their DNA, RNA, and protein signatures — to tailor treatment. The goal is to give the right drug, at the right dose, to the right patient, at the right time.

**DNA sequencing** is the process of reading the letters (bases) in a patient's genome. Modern sequencers produce FASTQ files, which contain millions of short DNA fragments called "reads," typically 100-300 base pairs long. Each read includes quality scores indicating the sequencer's confidence in each base call. A single whole-genome sequencing run generates approximately 100 gigabytes of FASTQ data.

**The reference genome (GRCh38)** is the standardized map of the human genome maintained by the Genome Reference Consortium. It contains approximately 3.2 billion base pairs organized across 24 chromosomes (22 autosomes plus X and Y). When we analyze a patient's DNA, we compare their reads against this reference to find differences — the variants that make each person unique and that may drive disease.

**Alignment** is the computational process of taking each short DNA read and determining where it belongs on the reference genome. BWA-MEM2 is the industry-standard alignment algorithm, optimized for speed and accuracy. It produces a BAM file — a binary, sorted, indexed file that records the position and orientation of every aligned read. This step transforms unordered raw reads into a structured, position-aware representation of the patient's genome.

**Variant calling** identifies the positions where the patient's DNA differs from the reference genome. DeepVariant, developed by Google, applies a deep convolutional neural network to this task. It converts read pileups at each genomic position into image-like tensors and classifies them as homozygous reference, heterozygous variant, or homozygous variant. The output is a VCF (Variant Call Format) file.

**A VCF file** is a structured text file that catalogs every variant detected in the patient's genome. Each line records the chromosome, position, reference allele, alternate allele, quality score, and filter status. A typical whole-genome VCF contains 4-5 million variants. Annotations are added to indicate clinical significance, population frequency, and predicted functional impact. This file becomes the foundation for all downstream clinical analysis.

**Somatic versus germline variants** represent two fundamentally different categories of genetic change. Germline variants are inherited — they exist in every cell of the body and were passed down from one or both parents. Somatic variants are acquired during a person's lifetime, typically found only in tumor cells. In cancer genomics, distinguishing the two is critical because somatic variants drive tumor growth while germline variants may indicate hereditary cancer predisposition.

**B-cell Acute Lymphoblastic Leukemia (B-ALL)** is the most common cancer in children, accounting for approximately 25% of all pediatric malignancies. It arises from the malignant proliferation of immature B-lymphocyte precursors (blasts) in the bone marrow, crowding out normal blood cell production. Symptoms include fatigue, bruising, infections, and bone pain. Cure rates now exceed 90% with modern risk-adapted therapy, but treatment intensity varies enormously based on molecular subtype.

**ETV6-RUNX1** is a chromosomal translocation, written as t(12;21), where portions of chromosomes 12 and 21 fuse together. This is the most common genetic alteration in pediatric ALL, found in approximately 25% of cases. The fusion creates an aberrant transcription factor that disrupts normal blood cell development. It is associated with a favorable prognosis — children with ETV6-RUNX1 have cure rates exceeding 95% with standard chemotherapy.

**IKZF1** (Ikaros) is a transcription factor essential for normal B-cell development. Deletions of IKZF1 are found in approximately 15% of pediatric B-ALL cases and serve as an unfavorable prognostic modifier. When IKZF1 deletion co-occurs with an otherwise favorable subtype like ETV6-RUNX1, it creates a tension that requires careful risk stratification. The combination is part of the "IKZF1plus" profile, which may warrant treatment intensification.

**Retrieval-Augmented Generation (RAG)** is an AI architecture that grounds large language model responses in factual evidence. Instead of relying solely on the model's training data (which may be outdated or hallucinated), the system first retrieves relevant documents from a curated knowledge base, then feeds those documents to the LLM as context for generating a response. This approach dramatically reduces hallucination and provides verifiable citations for every claim.

**Milvus** is an open-source vector database purpose-built for similarity search at scale. In this platform, it stores 3.56 million variant embeddings — mathematical representations of genomic variants and their clinical annotations. When a clinician asks a question, the system converts the query into the same vector space, finds the most similar variants in milliseconds, and returns the relevant clinical evidence. Milvus handles the computational geometry of high-dimensional search that would be impractical with traditional databases.

**MolMIM** is NVIDIA's molecular generation model from the BioNeMo framework. It uses a variational autoencoder architecture to generate novel drug-like molecules conditioned on a target protein structure. Given a protein binding pocket, MolMIM can explore chemical space to propose molecules that are predicted to bind effectively, have drug-like properties, and satisfy specified constraints like molecular weight and solubility.

**DiffDock** is a deep learning model for molecular docking — predicting how a small molecule (drug candidate) physically binds to a protein target. Unlike traditional docking software that exhaustively samples binding poses, DiffDock uses a score-based diffusion generative model to directly predict the most likely binding configuration. It generates multiple binding pose candidates and ranks them by predicted confidence, enabling rapid virtual screening of drug candidates.

**Pediatric safety filters** are computational checks specific to treating children, addressing the biological differences between pediatric and adult patients. These include blood-brain barrier (BBB) penetration assessment for brain tumors, cardiac liability screening via hERG channel binding prediction, hepatic immaturity considerations (children's livers metabolize drugs differently than adults'), and growth plate safety evaluation. A drug that is perfectly safe for an adult may be devastating for a growing child.

### Advanced Learning — End-to-End Precision Medicine Pipeline

**NVIDIA Parabricks architecture** provides GPU-accelerated equivalents of standard genomics tools from the GATK toolkit. By parallelizing alignment, sorting, duplicate marking, and variant calling across GPU cores, Parabricks achieves 30x whole-genome sequencing alignment in approximately 45 minutes compared to 8-12 hours on a high-end CPU cluster. The key acceleration comes from porting the Smith-Waterman alignment kernel and the pair-HMM calculations to CUDA, achieving near-linear speedup with GPU count. On the DGX Spark, the full FASTQ-to-VCF pipeline completes in under 2 hours.

**DeepVariant's CNN architecture** treats variant calling as an image classification problem. At each candidate variant site, it constructs a multi-channel pileup image encoding read bases, quality scores, strand information, and mapping quality. This image is fed through an Inception v3-derived convolutional neural network that classifies the site as homozygous reference, heterozygous, or homozygous variant. The image-based approach captures subtle alignment patterns (like strand bias and mapping artifacts) that rule-based callers miss, achieving a reported F1 score above 99.7% on benchmark datasets.

**BGE-small-en-v1.5** is a 384-dimensional dense embedding model from the Beijing Academy of AI (BAAI). Trained on a diverse English corpus using contrastive learning, it maps text passages to a high-dimensional vector space where semantically similar passages cluster together. For biomedical text retrieval, it outperforms older sentence-transformer models by better capturing domain-specific terminology and relationships. Its compact size (33M parameters) enables fast inference while maintaining competitive retrieval accuracy on the MTEB benchmark.

**IVF_FLAT indexing in Milvus** uses an Inverted File structure with flat (uncompressed) quantization to balance recall and search speed for the platform's 3.56 million vectors. The index partitions the vector space into a configurable number of Voronoi cells (typically 1024-4096), then at query time searches only the nearest cells rather than the entire database. The "FLAT" qualifier means vectors within each cell are stored without lossy compression, preserving exact distance calculations. This achieves sub-millisecond query latency while maintaining recall above 99% with an nprobe setting of 32-64 cells.

**COSINE similarity versus L2 (Euclidean) distance** affects how the system measures the relevance of retrieved documents. For text embeddings from models like BGE-small-en-v1.5, cosine similarity is preferred because it normalizes for document length — a short, focused clinical note and a long, detailed review article about the same variant will have high cosine similarity despite different magnitudes. L2 distance, by contrast, is sensitive to vector magnitude, which can cause longer documents to appear more "distant" purely due to their length rather than their semantic content.

**Claude LLM integration** in this platform follows a deterministic-probabilistic split architecture. All quantitative calculations — variant pathogenicity scores, drug interaction probabilities, risk stratifications — are performed by deterministic calculators with fixed, auditable logic. The LLM is used exclusively for natural language synthesis: assembling calculator outputs into coherent clinical narratives, generating explanations, and answering follow-up questions. This architecture ensures that clinical numbers never depend on probabilistic generation, while still providing the fluency and flexibility of a large language model.

**MolMIM's architecture** is a conditional variational autoencoder (CVAE) that operates in molecular latent space. The encoder maps SMILES strings to a continuous latent representation, conditioned on target protein features derived from the binding pocket. During generation, the model samples from this conditioned latent space and decodes back to SMILES strings, producing novel molecules that are biased toward the target's pharmacophore requirements. Latent space interpolation between known active molecules can generate structurally novel analogs with predicted binding affinity.

**DiffDock's diffusion process** is a score-based generative model applied to the SE(3) space of protein-ligand docking. Starting from a random initial placement of the ligand relative to the protein, the model iteratively denoises the configuration through learned score functions that guide translation, rotation, and torsion angle adjustments. Unlike traditional docking methods (AutoDock, Glide) that rely on physics-based scoring functions and exhaustive sampling, DiffDock learns the binding landscape from crystallographic data, achieving higher success rates on blind docking benchmarks while running 10-100x faster.

**Pediatric pharmacokinetic considerations** require fundamentally different modeling than adult drug development. Allometric scaling uses body surface area (BSA) rather than weight for dosing because metabolic rate scales with surface area. CYP enzyme ontogeny is critical: CYP3A4 reaches adult activity levels only by age 1-2, CYP1A2 by age 4-5, and CYP2D6 shows adult-level polymorphic expression by age 5. Renal maturation (GFR reaches adult levels by age 2) affects clearance of renally excreted drugs. The blood-brain barrier in young children is more permeable than in adults, which is both an opportunity (better drug penetration for brain tumors) and a risk (greater CNS toxicity from systemic drugs).

---

### Step 1: Landing Page Overview

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** This is your opening. Set the stage. Speak slowly. Let the audience absorb the platform.

**Actions:**

1. Switch to the **Landing Page** tab (http://localhost:8080)
2. Point to the **three engine cards** at the top:
   - Genomic Foundation Engine
   - RAG/Chat Intelligence Engine
   - Drug Discovery Engine
3. Scroll down to the **agent grid** — show all 11 agent tiles
4. Point out the **green health indicators** on each agent
5. Click on the **pipeline visualization** if available to show the data flow

**Say:**

> "This is the HCLS AI Factory. One platform. One machine. Eight engines and eight specialized AI agents. Everything you see here is running on a single NVIDIA DGX Spark. Today we are going to follow five children with cancer through this platform. Let's start with Evelyn."

**What you will see:** The Landing Page displays pipeline cards for Genomics, RAG/Chat, and Drug Discovery across the top. Below, an agent grid shows all 8 intelligence agents with real-time health status. Green circles indicate all services are operational.

> 📸 **SCREENSHOT PLACEHOLDER:** [Landing Page at http://localhost:8080 showing three engine cards (Genomics, RAG/Chat, Drug Discovery) across the top and the 11-agent grid below with all green health indicators]

> ⚠️ **IF THIS FAILS:** If the Landing Page does not load or agents show red/yellow status, run `bash health-monitor.sh --check-all` in a terminal. If individual agents are down, restart them with `docker compose -f docker-compose.dgx-spark.yml restart <service-name>`. While waiting, you can narrate the architecture from memory and move to Step 2.

> 🎤 **SAY:** "This is your command center. Every green dot represents a specialized AI agent running live on this single machine. In a traditional hospital, coordinating eleven specialists takes weeks of scheduling. Here, they are all available simultaneously, 24/7, for every patient."

---

### Step 2: Genomic Foundation Engine

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** If you have a live FASTQ run available, show it. If not, explain the
> process and skip to Step 3 with the pre-processed VCF. Most demo audiences will not need
> the full genomics walkthrough — the RAG and Oncology agents are the stars.

**Option A — If starting from FASTQ (live genomics run):**

1. Switch to the **Genomics Portal** tab (http://localhost:5000)
2. Show the genomics portal interface
3. Point out the pipeline stages: alignment (BWA-MEM2), variant calling (DeepVariant), annotation

**Say:**

> "Evelyn's tumor DNA and germline DNA have been sequenced. The raw FASTQ files — about 100 gigabytes — are being processed through NVIDIA Parabricks. On a DGX Spark, this takes about 2 hours. On a traditional CPU cluster, it would take 24 to 48 hours. The output is a VCF file — a catalog of every variant in Evelyn's genome."

**Option B — If starting from pre-processed VCF:**

1. Skip the Genomics Portal
2. Briefly acknowledge the genomics stage

**Say:**

> "Evelyn's genome has already been processed. We have her VCF file with 11.7 million variants cataloged and annotated. Let's see what the AI makes of it."

**What you will see:** The Genomics Portal shows pipeline status, processing stages, and output file locations. If a run is active, progress bars indicate current stage completion.

> 📸 **SCREENSHOT PLACEHOLDER:** [Genomics Portal at http://localhost:5000 showing pipeline stages (alignment, variant calling, annotation) with either active progress bars or completed VCF output summary]

> ⚠️ **IF THIS FAILS:** If the Genomics Portal is unresponsive, skip to Option B (pre-processed VCF). Say: "The genomics pipeline has already completed Evelyn's analysis." The Genomics Portal is the least critical tab for this demo — the RAG and Oncology agents are the stars. If you need to verify the VCF exists: `ls -la /home/adam/projects/hcls-ai-factory/data/*.vcf.gz`.

> 🎤 **SAY:** "What you are looking at is the first engine — Genomics. It takes raw DNA sequencing data, 100 gigabytes of it, and compresses it into a structured catalog of every variant in this child's genome. On GPU, this runs in two hours. On traditional hardware, it takes two days. That speed difference can be the difference between a timely diagnosis and a missed window."

---

### Step 3: RAG Chat — Initial Query

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** This is the first "wow" moment. The RAG engine will pull evidence from
> 3.56 million embedded variants in Milvus, combine it with Claude's reasoning, and return
> a clinically detailed response. Let the audience read the response for a moment before
> speaking.

**Actions:**

1. Switch to the **RAG Chat** tab (http://localhost:8501)
2. Note the sidebar: it shows variant database stats (3.56M vectors, 4.1M ClinVar records)
3. If demo question buttons are visible in the interface, you may click one — OR type the query below
4. Click into the **chat input field** at the bottom of the screen
5. Type the following query exactly:

```
8-year-old female with B-cell ALL. ETV6-RUNX1 fusion positive with IKZF1 deletion
detected. What is the molecular risk classification and what are the therapeutic options?
```

6. Press **Enter**
7. Wait for the response to stream in (typically 5-15 seconds)

**What you will see:** The RAG engine retrieves variant evidence from the Milvus vector database, synthesizes it with Claude, and returns a structured clinical response. The response will include:

- **Molecular risk classification:** ETV6-RUNX1 is standard risk, but IKZF1 deletion upgrades to high risk (IKZF1plus phenotype)
- **Therapeutic options:** Intensified chemotherapy per COG protocols, consideration of blinatumomab, MRD-directed therapy
- **Evidence sources:** ClinVar entries, published literature, variant annotations
- **Confidence scoring:** Displayed alongside the response

> **PRESENTER NOTE:** Point out the evidence citations. This is not a hallucination — the
> system is grounding every claim in retrieved documents. That is the difference between
> a chatbot and a clinical decision support system.

**Say:**

> "Notice the evidence sources. The RAG engine searched 3.56 million variant embeddings and pulled specific ClinVar annotations, published guidelines, and molecular classification criteria. This is retrieval-augmented generation — every claim is grounded in evidence."

> 📸 **SCREENSHOT PLACEHOLDER:** [RAG Chat at http://localhost:8501 showing the query about Evelyn's ETV6-RUNX1/IKZF1 profile and the structured clinical response with evidence citations, confidence scores, and molecular risk classification visible]

> ⚠️ **IF THIS FAILS:** If the RAG Chat returns an error or hangs, check Milvus connectivity: `curl http://localhost:19530/v1/vector/collections`. If Milvus is down, the RAG engine may still return a knowledge-base-grounded response using Claude alone — it will just lack the vector-retrieved citations. If the chat is completely unresponsive, restart: `docker compose -f docker-compose.dgx-spark.yml restart rag-chat`. While waiting, narrate what the response would contain and proceed to the Oncology agent.

> 🎤 **SAY:** "This is not a chatbot guessing. Every sentence you see is backed by retrieved evidence from 3.56 million embedded variants. The system identified that Evelyn's ETV6-RUNX1 fusion, normally a favorable marker, becomes high-risk when combined with IKZF1 deletion. It pulled the specific ClinVar entries, the COG guidelines, and the published literature to support that conclusion. This is the difference between AI as a toy and AI as a clinical tool."

---

### Step 4: Precision Oncology Agent — Molecular Classification

> **TIMING: 5 minutes**
> **PRESENTER NOTE:** This is the centerpiece of Demo 1. You are building a complete molecular
> tumor board packet from the GUI. Take your time entering the data — the audience needs to
> see that this is a structured clinical workflow, not just a chat interface.

**Actions:**

1. Switch to the **Oncology Agent** tab (http://localhost:8526)
2. Click the **"Case Workbench"** tab (first of 5 tabs: Case Workbench, Evidence Explorer, Trial Finder, Therapy Ranker, Outcomes Dashboard)

**Fill in the case form:**

| Field | Value | Notes |
|-------|-------|-------|
| **Patient ID** | `PEDS-ALL-001` | Type into the text input |
| **Cancer Type** | Select from the 19-type selectbox | Choose the closest leukemia/ALL option |
| **Stage** | Select from the Stage selectbox | Choose the appropriate stage |
| **TMB** | `2.1` | Type into the number input (low for ALL) |
| **MSI Status** | Select `MSS` (Microsatellite Stable) | From the selectbox |
| **PD-L1 TPS** | `5` | Type into the number input |
| **HRD Score** | `0` | Type into the number input (not applicable for ALL) |
| **Prior Therapies** | Leave empty or select if applicable | Multiselect field |

3. In the **Variant entry** section:
   - Select the **"Manual"** radio button (not VCF)
   - Add the following variants (one row at a time):

| Gene | Variant | Type |
|------|---------|------|
| `ETV6-RUNX1` | `t(12;21)(p13;q22)` | `fusion` |
| `IKZF1` | `exon 4-7 deletion` | `deletion` |
| `CDKN2A` | `homozygous deletion` | `deletion` |

4. Click the **"Create Case"** button
5. Wait for the case to be created (green confirmation message)
6. Click the **"Generate MTB Packet"** button
7. Wait for the molecular tumor board packet to generate (10-20 seconds)

**What you will see:** A complete molecular tumor board (MTB) packet generates with:

- **Risk classification:** ETV6-RUNX1 is typically favorable prognosis, BUT the IKZF1 deletion upgrades this to high risk under the IKZF1plus classification
- **Therapy ranking:** Prioritized treatment options with evidence levels
- **Trial matches:** Automatically matched clinical trials based on molecular profile
- **Actionable variants:** Each variant annotated with clinical significance
- **Molecular summary:** Integrated interpretation of the variant landscape

> **PRESENTER NOTE:** This is the key clinical insight to emphasize: ETV6-RUNX1 alone would
> be standard risk (excellent prognosis). But the combination with IKZF1 deletion changes
> everything. The platform catches this interaction automatically. In a manual review, this
> could be missed.

**Say:**

> "Look at the risk classification. ETV6-RUNX1 fusion alone would put Evelyn in the standard risk category — 90% cure rate. But the IKZF1 deletion changes the picture entirely. This is the IKZF1plus phenotype, which upgrades her to high risk. The platform caught this interaction automatically and is recommending treatment intensification. In a manual chart review, this combination might be missed."

> 📸 **SCREENSHOT PLACEHOLDER:** [Oncology Agent Case Workbench at http://localhost:8526 showing the completed case form for PEDS-ALL-001 with ETV6-RUNX1, IKZF1, and CDKN2A variants entered, and the generated MTB Packet displaying risk classification as high-risk IKZF1plus, therapy ranking, and actionable variants]

> ⚠️ **IF THIS FAILS:** If the case creation or MTB Packet generation fails, check the Oncology agent health: `curl http://localhost:8527/health`. If the agent is running but the packet times out, click "Generate MTB Packet" again — the first attempt sometimes takes longer as the model warms up. If the variant entry fields do not accept the format shown, try entering just the gene name (e.g., "ETV6-RUNX1") without the cytogenetic notation. You can also switch to the Evidence Explorer tab and query the same clinical scenario as a fallback.

> 🎤 **SAY:** "What you are seeing here is the platform identifying that Evelyn's ETV6-RUNX1 fusion — normally a favorable marker — is complicated by an IKZF1 deletion. This changes her risk classification from standard to high, which means she needs more intensive therapy. The system caught this in seconds. In a traditional workflow, a molecular pathologist would spend 30 to 60 minutes assembling this tumor board packet manually, and the critical interaction between these two variants could be overlooked."

---

### Step 5: Oncology — Trial Finder

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** Quick transition. Show the Trial Finder tab to demonstrate that the
> platform doesn't just diagnose — it connects patients to trials.

**Actions:**

1. Stay on the **Oncology Agent** (http://localhost:8526)
2. Click the **"Trial Finder"** tab
3. Select the Cancer Type from the available options
4. Enter Age: `8`
5. Click **"Find Trials"** (or equivalent search button)

**What you will see:** A list of matched clinical trials with:

- Trial identifiers (NCT numbers)
- Eligibility scores (percentage match)
- Phase and status
- Key inclusion/exclusion criteria highlights
- Distance to trial sites (if location data is available)

**Say:**

> "The Trial Finder matched Evelyn to active clinical trials based on her exact molecular profile. Notice the eligibility scores — the platform has already pre-screened inclusion and exclusion criteria."

> 📸 **SCREENSHOT PLACEHOLDER:** [Oncology Agent Trial Finder tab at http://localhost:8526 showing matched clinical trials for pediatric ALL with NCT numbers, eligibility scores, phase/status, and key inclusion criteria highlighted]

> ⚠️ **IF THIS FAILS:** If the Trial Finder returns no results or errors, verify the Oncology agent is responsive: `curl http://localhost:8527/health`. If the trial database is not loaded, the agent may return a smaller set of results or a general recommendation. You can also use the Clinical Trial Agent (http://localhost:8128) as a backup for trial matching — it has a dedicated Patient-Trial Matching tab.

> 🎤 **SAY:** "The platform does not stop at diagnosis. It automatically matches Evelyn's exact molecular profile — ETV6-RUNX1 positive, IKZF1 deleted, 8 years old — against active clinical trials. The eligibility scores tell the oncologist which trials this child qualifies for without manually reviewing hundreds of inclusion and exclusion criteria. That screening process, which normally takes hours, happened in seconds."

---

### Step 6: Drug Discovery Portal

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** This is the third engine. Show the drug discovery interface but don't
> run a full generation during the live demo (it takes time). Use pre-loaded results if
> available, or walk through the interface explaining each component.

**Actions:**

1. Switch to the **Drug Discovery Portal** tab (http://localhost:8510)
2. Show the drug discovery interface layout
3. Walk through the three components:
   - **Target Selection:** Select the molecular target (e.g., ETV6-RUNX1 fusion product)
   - **MolMIM Generation:** Show the molecular generation interface (BioNeMo NIM)
   - **DiffDock Scoring:** Show binding affinity prediction
4. Point out the **pediatric safety filters** if visible (growth plate impact, BBB penetration, hepatotoxicity, hERG channel)

**What you will see:** The Drug Discovery portal displays:

- Target protein structure visualization
- Generated molecular candidates from MolMIM
- DiffDock binding scores
- RDKit property calculations (ADMET, Lipinski)
- Pediatric safety filter results (critical for pediatric oncology)

**Say:**

> "This is the third engine — Drug Discovery. When standard therapies are not enough, the platform can generate novel molecular candidates using BioNeMo's MolMIM, score their binding affinity with DiffDock, and filter for pediatric safety. This is generative chemistry on demand."

> 📸 **SCREENSHOT PLACEHOLDER:** [Drug Discovery Portal at http://localhost:8510 showing the three-component interface: target protein structure visualization, MolMIM molecular generation panel, and DiffDock binding scores with pediatric safety filter results]

> ⚠️ **IF THIS FAILS:** If the Drug Discovery Portal is slow or unresponsive, the BioNeMo NIMs may need additional startup time. Check: `curl http://localhost:8510/health`. If the portal loads but generation is slow, walk through the interface and explain each component without running a live generation. Pre-loaded results (if available) are sufficient for the demo. Say: "The generation takes a few minutes, so let me show you results from a previous run."

> 🎤 **SAY:** "This is where the platform goes beyond what any human team can do in real time. When standard chemotherapy is not enough, the Drug Discovery engine can generate entirely new molecular candidates, predict how they bind to the target protein, and — critically for pediatric patients — filter for safety in growing children. Growth plate toxicity, cardiac liability, liver immaturity — all checked automatically before a candidate is ever proposed."

---

### Step 7: Demo 1 Summary

> **TIMING: 1 minute**

**Actions:**

1. Switch back to the **Landing Page** (http://localhost:8080)
2. Gesture across the three engine cards

**Say:**

> "In 15 minutes, we took Evelyn from raw DNA through molecular classification, therapy ranking, trial matching, and drug candidate generation. Three engines. One platform. One child's life. Let's see what happens when the tumor board convenes."

> 📸 **SCREENSHOT PLACEHOLDER:** [Landing Page at http://localhost:8080 showing the three engine cards with the full pipeline visible — this is the "hero shot" summarizing Demo 1]

> ⚠️ **IF THIS FAILS:** This is a summary step with no live interaction. If the Landing Page has become unresponsive during the demo, simply deliver the closing statement without switching tabs. The verbal summary is the most important element here.

> 🎤 **SAY:** "In 15 minutes, we followed one child's DNA through the entire precision medicine pipeline. From raw sequencing data to molecular classification to therapy ranking to clinical trial matching to novel drug candidate generation. Three engines, eight agents, one machine. Now let us see what happens when Evelyn's treatment hits a complication and the tumor board needs to convene."

---

## Demo 2: "The 30-Second Tumor Board" — Multi-Agent Coordination

### Patient: Evelyn Chen (continued), Day 29 of Induction, MRD Positive (0.1%)

**Duration:** 12-15 minutes
**Story arc:** Evelyn's Day 29 MRD is positive. Five specialist agents convene a virtual tumor board to determine next steps.

**Browser tabs needed:**
- Oncology (http://localhost:8526)
- Biomarker (http://localhost:8528)
- CAR-T (http://localhost:8521)
- Single-Cell (http://localhost:8130)
- Clinical Trial (http://localhost:8128)

### Foundations Learning — Multi-Agent Tumor Board

**A tumor board** is a multidisciplinary conference where cancer specialists gather to review a patient's case and recommend treatment. A typical pediatric oncology tumor board includes a medical oncologist, pathologist, radiologist, geneticist, pharmacist, and often a surgeon. Each brings a different perspective — the pathologist interprets biopsy results, the geneticist reviews molecular profiling, the pharmacist evaluates drug interactions and dosing. The collective expertise produces better decisions than any single specialist.

**Why tumor boards are slow** comes down to scheduling logistics and information fragmentation. Coordinating 5-7 busy specialists often takes 1-2 weeks from referral to meeting. Each specialist must independently review the case beforehand, pulling information from different systems (imaging, pathology, genetics, pharmacy). Critical treatment decisions are delayed while administrative coordination happens. This platform compresses that process by having AI agents perform the specialist reviews simultaneously.

**Minimal Residual Disease (MRD)** is the single most important predictor of relapse in pediatric ALL. After initial chemotherapy (induction), flow cytometry or PCR-based assays measure the fraction of remaining leukemia cells in the bone marrow. Even a tiny fraction of surviving leukemia cells can seed relapse. MRD testing has revolutionized risk stratification by providing a direct, functional readout of how well the leukemia responded to treatment.

**MRD thresholds** define clinical decision boundaries. An MRD level below 0.01% (one leukemia cell per 10,000 normal cells, or 10^-4) is considered MRD negative — these patients have excellent outcomes with standard therapy. An MRD level at or above 0.01% is MRD positive, indicating a higher relapse risk that may warrant treatment intensification. MRD at or above 1% represents very high residual disease burden and is associated with particularly poor prognosis, often triggering consideration of immunotherapy or transplant.

**CD19** is a cell surface protein (receptor) found on virtually all B lymphocytes, from early development through maturity. Because B-ALL cells are malignant B-cell precursors, they almost universally express CD19, making it an ideal therapeutic target. CD19 is the target antigen for blinatumomab (a bispecific antibody), tisagenlecleucel (CAR-T therapy), and several other immunotherapies. Its near-universal expression on B-ALL cells and relative restriction to the B-cell lineage make it the cornerstone of B-ALL immunotherapy.

**CD22** is another cell surface protein expressed on B cells, functioning as an inhibitory receptor that modulates B-cell activation. It serves as a critical backup target when CD19-directed therapies fail due to antigen escape — a phenomenon where leukemia cells lose CD19 expression under therapeutic pressure. Inotuzumab ozogamicin (an antibody-drug conjugate targeting CD22) and CD22-directed CAR-T cells are being evaluated as second-line immunotherapies. Dual CD19/CD22 targeting strategies aim to prevent antigen escape entirely.

**CAR-T therapy (Chimeric Antigen Receptor T-cell therapy)** is a form of immunotherapy where a patient's own T cells are harvested, genetically engineered in a laboratory to express a synthetic receptor targeting a specific cancer antigen, expanded to hundreds of millions of cells, and then infused back into the patient. The CAR construct typically includes an extracellular antigen-binding domain (derived from an antibody), a transmembrane domain, and intracellular signaling domains that activate the T cell upon antigen binding. Tisagenlecleucel (Kymriah) was the first FDA-approved CAR-T for pediatric ALL.

**Cytokine Release Syndrome (CRS)** is the most common serious adverse effect of CAR-T therapy. When engineered T cells encounter and bind their target antigen on tumor cells, they activate and release large quantities of inflammatory cytokines — primarily IL-6, interferon-gamma, and TNF-alpha. This systemic inflammatory response causes fever (often the first sign), hypotension, capillary leak, and in severe cases, organ dysfunction. CRS typically occurs within the first 1-14 days after infusion and correlates with tumor burden — higher disease burden means more targets for CAR-T cells to engage simultaneously.

**The tumor microenvironment (TME)** is the complex ecosystem surrounding and interacting with tumor cells. It includes immune cells (T cells, macrophages, natural killer cells), stromal cells (fibroblasts), blood vessels, extracellular matrix, and a diverse milieu of signaling molecules (cytokines, chemokines, growth factors). The composition of the TME profoundly influences treatment response — a "hot" TME with abundant immune cells may respond well to immunotherapy, while a "cold" or immunosuppressive TME may resist it.

**An "immune desert"** is a tumor microenvironment characterized by very few infiltrating immune cells. This phenotype suggests that the immune system is either not recognizing the tumor as foreign or is being actively excluded by the tumor's defense mechanisms. In the context of CAR-T therapy, an immune desert TME is actually less concerning than an immunosuppressive TME, because CAR-T cells bring their own effector function. However, the TME composition still affects CAR-T persistence and expansion at the tumor site.

**Cross-agent coordination** in this platform works through a RESTful HTTP orchestration layer. The `/integrated-assessment` endpoint on the Oncology agent acts as the coordinator — when triggered, it sends parallel HTTP requests to peer agents (Biomarker, CAR-T, Single-Cell, Clinical Trial), collects their individual assessments, and synthesizes a unified recommendation. Each agent operates independently with its own specialized knowledge base and analytical tools, but the orchestration layer creates the collaborative dynamic of a human tumor board.

### Advanced Learning — Multi-Agent Tumor Board

**IKZF1plus** is a refined prognostic classifier that goes beyond simple IKZF1 deletion. It is defined by IKZF1 deletion occurring simultaneously with at least one additional genetic alteration: deletion of CDKN2A (a cell cycle regulator), deletion of PAX5 (a B-cell transcription factor), or specific cytokine receptor pathway alterations (CRLF2 rearrangement, JAK mutations). The IKZF1plus profile reclassifies patients who would otherwise be standard risk to very high risk, warranting treatment intensification. Studies from the AIEOP-BFM consortium demonstrated that IKZF1plus patients had a 5-year event-free survival of only 53% versus 87% for IKZF1-deleted patients without the additional hits.

**MRD kinetics** provide more prognostic information than any single MRD measurement. The rate of MRD clearance — comparing levels at Day 8 (peripheral blood), Day 15 (bone marrow), and Day 29 (end of induction) — models an exponential decay curve whose parameters predict outcome more accurately than the Day 29 value alone. Rapid early clearance (achieving MRD negativity by Day 15) identifies excellent-prognosis patients who may be candidates for treatment de-escalation. Conversely, slow clearance (still MRD-positive at Day 29) identifies patients needing intensification regardless of the final MRD level achieved.

**CD19 antigen density** measured by flow cytometry as Mean Fluorescence Intensity (MFI) directly impacts the efficacy of CD19-directed therapies. An MFI above 1,000 correlates strongly with CAR-T response, while values below this threshold are associated with significantly reduced efficacy — the CAR construct simply cannot engage enough antigen molecules to trigger robust T-cell activation. Antigen density is not binary (present/absent) but exists on a continuum, and this quantitative assessment is critical for predicting which patients will benefit from CD19-directed therapy versus needing alternative approaches.

**Blinatumomab's mechanism** as a Bispecific T-cell Engager (BiTE) is elegantly simple in concept. The molecule has two antibody-derived binding domains connected by a flexible linker: one arm binds CD3 on any nearby T cell, and the other binds CD19 on the leukemia cell. By physically cross-linking these two cells, blinatumomab creates an artificial immunological synapse — forcing the T cell to activate and kill the leukemia cell regardless of the T cell's native antigen specificity. This bypasses the need for MHC presentation, which leukemia cells often downregulate as an immune evasion strategy.

**Single-cell RNA sequencing in leukemia** provides unprecedented resolution into tumor heterogeneity. By profiling gene expression in thousands of individual cells, scRNA-seq can identify pre-existing resistant subclones before treatment begins. For CD19-directed therapy, this means detecting CD19-dim populations — small subsets of leukemia cells with low CD19 transcript levels that would survive CAR-T therapy and seed relapse. Trajectory analysis can also reveal differentiation states that predict therapy resistance, and clustering algorithms can identify rare populations (down to 0.1% of cells) that would be invisible to bulk sequencing methods.

**TME deconvolution using CIBERSORTx** is a computational method that estimates the fraction of different immune cell types present in a tumor sample from bulk RNA sequencing data. The algorithm uses reference gene expression signatures for 22 immune cell types and applies support vector regression to infer cell type proportions. This provides immune cell composition information without the expense of single-cell sequencing or flow cytometry. CIBERSORTx has been validated against flow cytometry ground truth in multiple tumor types and can reveal whether a patient's TME is inflamed, immune-excluded, or desert — information that guides immunotherapy decisions.

**Cross-agent HTTP orchestration** in this platform uses asynchronous requests with configurable timeout settings (default 30 seconds per agent). When the Oncology agent's `/integrated-assessment` endpoint is called, it dispatches concurrent requests to all participating agents using Python's asyncio/aiohttp stack. Each agent processes its query independently against its specialized vector store and analytical tools, then returns a structured JSON response. If any agent is unavailable or times out, the system implements graceful degradation — the assessment proceeds with available results and flags the missing perspective. The synthesis step uses a weighted aggregation model where agent confidence scores modulate their contribution to the final recommendation.

---

### Step 1: Oncology — Evidence Explorer

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** Set the clinical scenario. MRD positivity at Day 29 is a critical
> decision point in ALL treatment. The audience needs to feel the urgency.

**Actions:**

1. Switch to the **Oncology Agent** tab (http://localhost:8526)
2. Click the **"Evidence Explorer"** tab (second of 5 tabs)
3. Type the following query into the search/question field:

```
Evelyn is day 29 of induction for B-ALL with MRD 0.1%. ETV6-RUNX1 positive with
IKZF1 deletion. What are the treatment intensification options?
```

4. Click **"Ask"** (or press Enter)
5. Wait for the response (5-15 seconds)

**What you will see:** Evidence-based recommendations for MRD-positive ALL including:

- Definition of MRD positivity threshold (0.01% vs 0.1%)
- COG protocol-specific intensification options
- Blinatumomab consolidation data
- Transplant vs. immunotherapy decision framework
- IKZF1plus-specific MRD guidance

**Say:**

> "Evelyn is 29 days into treatment and her MRD — minimal residual disease — is 0.1%. That means there are still leukemia cells hiding. The Oncology agent is pulling COG protocol guidelines, IKZF1plus-specific data, and blinatumomab trial results to recommend intensification options. But we need more information. Let's consult the other specialists."

> 📸 **SCREENSHOT PLACEHOLDER:** [Oncology Agent Evidence Explorer tab at http://localhost:8526 showing the MRD query response with COG protocol intensification options, blinatumomab data, and IKZF1plus-specific MRD guidance]

> ⚠️ **IF THIS FAILS:** If the Evidence Explorer does not return results, check the Oncology agent: `curl http://localhost:8527/health`. If the query times out, try a shorter query: "MRD positive ALL treatment options IKZF1plus". If the Oncology agent is completely down, you can use the RAG Chat (http://localhost:8501) with the same query as a fallback — it draws from the same evidence base.

> 🎤 **SAY:** "Day 29 is decision day for leukemia treatment. Evelyn still has measurable disease — 0.1% MRD. That number sounds tiny, but it means there are still millions of leukemia cells in her bone marrow. The Oncology agent immediately pulls the relevant COG protocols and shows us that for an IKZF1plus patient with persistent MRD, treatment intensification is indicated. But one agent is not enough — we need the full tumor board."

---

### Step 2: Biomarker Agent — CD19/CD22 Assessment

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** This step bridges to the CAR-T discussion. CD19 expression is the
> prerequisite for CAR-T therapy.

**Actions:**

1. Switch to the **Biomarker Agent** tab (http://localhost:8528)
2. Note the 8 tabs available: Biomarker Analysis, Biological Age, Disease Risk, PGx Profile, Evidence Explorer, Reports, Patient 360, Longitudinal
3. Click the **"Biomarker Analysis"** tab (or **"Evidence Explorer"** tab)
4. If a **sample patient quick-load** button is visible in the sidebar, you may use it to pre-populate fields
5. Enter query:

```
CD19 and CD22 expression levels for 8-year-old with B-cell ALL. Is this patient
eligible for CD19-targeted immunotherapy?
```

6. Submit the query

**What you will see:** Biomarker panel results showing:

- CD19 expression status and intensity
- CD22 expression as backup target
- Immunophenotype summary
- Eligibility assessment for immunotherapy (blinatumomab, CAR-T)
- Relevant biomarker thresholds and clinical cut-points

**Say:**

> "Before we can consider CAR-T therapy, we need to confirm the target is there. The Biomarker agent confirms CD19 expression — the handle that CAR-T cells will grab onto. It also checks CD22 as a backup target in case of CD19 escape. Both are positive. Let's talk to the CAR-T specialist."

> 📸 **SCREENSHOT PLACEHOLDER:** [Biomarker Agent at http://localhost:8528 showing CD19/CD22 expression analysis with immunophenotype summary, expression intensity levels, and immunotherapy eligibility assessment]

> ⚠️ **IF THIS FAILS:** If the Biomarker agent is unresponsive, check: `curl http://localhost:8528/health`. If the query returns no results, try the Evidence Explorer tab instead of Biomarker Analysis. You can also state the expected result verbally: "The biomarker panel confirms CD19 positivity with strong expression — the target for CAR-T therapy is present." Then proceed to the CAR-T agent.

> 🎤 **SAY:** "Before we can even discuss CAR-T therapy, we need to confirm the target exists on Evelyn's leukemia cells. CD19 is the surface protein that CAR-T cells are engineered to recognize and attack. The Biomarker agent confirms it is there, with strong expression. It also checks CD22 as a backup — because if the leukemia ever loses CD19 to escape the CAR-T cells, we need a second target ready."

---

### Step 3: CAR-T Agent — Eligibility Assessment

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** The CAR-T agent uses a chat interface with collection filtering and
> Deep Research mode. For the demo, use the standard chat. Save Deep Research mode for
> Demo 5 when we do the full CAR-T workup.

**Actions:**

1. Switch to the **CAR-T Agent** tab (http://localhost:8521)
2. Note the interface: chat input, collection filtering options, Deep Research mode toggle, citation scoring
3. Leave Deep Research mode **off** for now
4. Click into the chat input field
5. Type:

```
12-year-old with relapsed/refractory B-ALL, CD19 positive >95%, failed 2 prior lines.
Is this patient eligible for tisagenlecleucel? What is the CRS risk?
```

> **PRESENTER NOTE:** We are using a slightly older patient profile here (12yo) to match
> the tisagenlecleucel label. Evelyn is 8, but we are exploring the CAR-T pathway as a
> potential future option if she relapses.

6. Press **Enter**
7. Wait for the response with citations

**What you will see:** Detailed CAR-T eligibility assessment including:

- Tisagenlecleucel (Kymriah) eligibility criteria check
- Age requirement (up to 25 years for ALL)
- Prior therapy requirement (2+ lines — met)
- CD19 positivity confirmation
- CRS risk profile (Grade 1-2: 50-60%, Grade 3-4: 20-30%)
- ICANS risk profile
- Citation scores linking to source evidence

**Say:**

> "The CAR-T agent has assessed eligibility against the tisagenlecleucel label. CD19-positive, failed two prior lines, age-appropriate — eligible. But look at the risk profile: 20-30% chance of severe cytokine release syndrome. This is why we need the full team."

> 📸 **SCREENSHOT PLACEHOLDER:** [CAR-T Agent at http://localhost:8521 showing the eligibility assessment response with tisagenlecleucel criteria checklist, CRS risk profile (Grade 1-4 percentages), ICANS risk, and citation scores]

> ⚠️ **IF THIS FAILS:** If the CAR-T agent is unresponsive, check: `curl http://localhost:8521/health`. If the chat hangs, refresh the page and re-enter the query. If the agent is completely down, narrate the eligibility criteria from the presenter notes and skip to the Single-Cell agent. The key point to convey is: eligible but with significant CRS risk.

> 🎤 **SAY:** "The CAR-T agent checks every box on the tisagenlecleucel eligibility criteria — age, prior therapy lines, CD19 status, performance status. The patient qualifies. But look at the risk numbers: 20 to 30 percent chance of severe cytokine release syndrome, where the immune system essentially goes into overdrive. This is exactly why a tumor board exists — no single specialist makes this decision alone."

---

### Step 4: Single-Cell — TME Profiling

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** The TME Profiler shows the tumor microenvironment — critical for
> understanding why some patients respond to immunotherapy and others don't.

**Actions:**

1. Switch to the **Single-Cell Agent** tab (http://localhost:8130)
2. Note the 5 tabs: Dashboard, Evidence Explorer, TME Profiler, Workflow Runner, Reports & Export
3. Click the **"TME Profiler"** tab
4. Select **Tumor Type** from the dropdown (choose the closest leukemia/ALL option)
5. Adjust the **cell proportion sliders** to reflect ALL blast predominance:

| Cell Type | Proportion | Notes |
|-----------|------------|-------|
| Malignant | `0.90` | Dominant — this is active leukemia |
| CD8+ T Cells | `0.02` | Severely depleted |
| CD4+ T Cells | `0.03` | Reduced |
| Tregs | `0.01` | Minimal |
| NK Cells | `0.01` | Minimal |
| Macrophages (M1) | `0.01` | Present |
| Macrophages (M2) | `0.01` | Present |
| B Cells | `0.01` | Suppressed by disease |
| Other | `0.01` | Stromal/other |

6. Click **"Profile TME"** (or equivalent analysis button)

**What you will see:** TME classification and analysis showing:

- TME category (likely "immune desert" or "immune excluded" given the blast predominance)
- Immune score (low)
- Predicted response to checkpoint inhibitors (poor — blast predominance suppresses immune infiltration)
- Predicted response to CAR-T (better — engineered cells bypass microenvironment)
- Therapy predictions based on TME composition

**Say:**

> "The tumor microenvironment tells us why standard immunotherapy might struggle. 90% blasts, barely any T cells left. The immune system is overwhelmed. But CAR-T therapy doesn't rely on the existing immune system — it engineers new soldiers. This TME profile actually supports the CAR-T pathway."

> 📸 **SCREENSHOT PLACEHOLDER:** [Single-Cell Agent TME Profiler at http://localhost:8130 showing cell proportion sliders set to 90% malignant/2% CD8+, TME classification as "immune desert," immune score, and therapy predictions favoring CAR-T over checkpoint inhibitors]

> ⚠️ **IF THIS FAILS:** If the Single-Cell agent does not load or the TME Profiler tab is missing, check: `curl http://localhost:8130/health`. If the sliders do not respond, try refreshing the page. If the agent is completely unavailable, describe the TME concept verbally: "In Evelyn's bone marrow, 90% of the cells are leukemia blasts. The immune system is overwhelmed. This is why CAR-T works where other immunotherapies fail — it brings its own army." Then proceed to Step 5.

> 🎤 **SAY:** "This is the tumor microenvironment — the ecosystem around the cancer. Evelyn's marrow is 90 percent blasts. Her own immune cells are almost gone. Standard immunotherapy like checkpoint inhibitors would fail here because there are no T cells left to unleash. But CAR-T therapy does not need the existing immune system — it engineers brand-new cancer-killing cells. This microenvironment data actually makes the case for CAR-T stronger."

---

### Step 5: Clinical Trial — Matching

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Clinical Trial Agent** tab (http://localhost:8128)
2. Note the 5 tabs: Dashboard, Evidence Explorer, Patient-Trial Matching, Protocol Optimization, Reports
3. Click the **"Patient-Trial Matching"** tab (or use the Evidence Explorer)
4. Enter query:

```
Pediatric ALL clinical trials with MRD-directed therapy for 8-year-old with
ETV6-RUNX1 positive, IKZF1 deletion, MRD 0.1% at Day 29
```

5. Submit the query

**What you will see:** Matched clinical trials including:

- COG AALL1732 or successor protocols for high-risk ALL
- Blinatumomab consolidation trials
- MRD-directed de-escalation/escalation trials
- Phase and enrollment status
- Key eligibility criteria

**Say:**

> "The Clinical Trial agent found COG trials specifically designed for MRD-positive, high-risk ALL. These are the trials Evelyn's oncologist would discuss at the real tumor board."

> 📸 **SCREENSHOT PLACEHOLDER:** [Clinical Trial Agent Patient-Trial Matching tab at http://localhost:8128 showing matched pediatric ALL trials with MRD-directed therapy, COG protocol numbers, phase/enrollment status, and key eligibility criteria]

> ⚠️ **IF THIS FAILS:** If the Clinical Trial agent does not return results, check: `curl http://localhost:8128/health`. Try the Evidence Explorer tab as a fallback with a simpler query: "pediatric ALL MRD positive clinical trials". If the agent is completely down, you can use the Oncology Agent's Trial Finder tab (shown in Demo 1 Step 5) as an alternative source for trial matching.

> 🎤 **SAY:** "The fifth specialist at our virtual tumor board is the Clinical Trial agent. It matched Evelyn to COG protocols specifically designed for MRD-positive, IKZF1plus ALL. These are not generic search results — they are trials where Evelyn meets the molecular and clinical eligibility criteria. Her oncologist would present these exact options at the real tumor board meeting."

---

### Step 6: Demo 2 Summary

> **TIMING: 1 minute**
> **PRESENTER NOTE:** This is the multi-agent payoff. Arrange your browser windows so the
> audience can see multiple agents at once if possible (tile 2-3 windows side by side).

**Actions:**

1. If possible, tile 2-3 agent windows side by side on screen
2. Gesture across all of them

**Say:**

> "This is what a virtual tumor board looks like. Five specialist agents — oncology, biomarker, CAR-T, single-cell, and clinical trial — all working from the same evidence base, all coordinating in minutes instead of weeks. No scheduling conflicts. No waiting for the next multidisciplinary meeting. Every specialist available, every time, for every child."

> 📸 **SCREENSHOT PLACEHOLDER:** [Multiple agent windows tiled side by side (2-3 visible) showing Oncology, CAR-T, and Clinical Trial agents with their respective assessments — the "virtual tumor board" view]

> ⚠️ **IF THIS FAILS:** This is a summary step. If you cannot tile windows, simply switch between 2-3 agent tabs quickly to show the breadth of assessments. The verbal summary is the key deliverable here — the visual tiling is a bonus.

> 🎤 **SAY:** "What you just witnessed took 12 minutes. Five specialist agents — oncology, biomarker, CAR-T, single-cell, and clinical trial — all consulted, all drawing from the same evidence base, all coordinated automatically. In a real hospital, scheduling these five specialists into the same room takes one to two weeks. For a child with persistent leukemia, those weeks matter. This platform makes every specialist available, every time, for every child."

---

## Demo 3: "Protecting the Survivor" — Cardiotoxicity Prevention in Pediatric Chemotherapy

### Patient: Marcus Williams, 6 years old, Male, Stage 4 MYCN-Amplified Neuroblastoma

**Duration:** 15 minutes
**Story arc:** Marcus needs aggressive chemotherapy, but every drug has toxicity. Six agents coordinate to prevent cardiotoxicity, neuropathy, and immune-related adverse events before they happen.

**Browser tabs needed:**
- PGx (http://localhost:8507)
- Biomarker (http://localhost:8528)
- Cardiology (http://localhost:8536)
- Neurology (http://localhost:8529)
- Autoimmune (http://localhost:8531)
- Oncology (http://localhost:8526)

### Foundations Learning — Cardiotoxicity Prevention

**Anthracyclines** are a class of chemotherapy drugs — including doxorubicin, daunorubicin, and epirubicin — that rank among the most effective anticancer agents ever developed. They work by intercalating into DNA and inhibiting topoisomerase II, preventing cancer cells from replicating. Anthracyclines are backbone drugs in protocols for leukemia, lymphoma, sarcoma, and neuroblastoma. However, their effectiveness comes with a serious cost: they are directly toxic to the heart muscle, and this damage is cumulative and largely irreversible.

**How anthracyclines damage the heart** involves multiple molecular mechanisms, but the primary pathway is generation of reactive oxygen species (free radicals). When anthracyclines interact with iron inside heart muscle cells (cardiomyocytes), they catalyze the formation of hydroxyl radicals through the Fenton reaction. These free radicals damage mitochondrial membranes, disrupt calcium handling, and trigger cardiomyocyte apoptosis (programmed cell death). Unlike most other cell types, adult cardiomyocytes have virtually no regenerative capacity — each cell lost is lost permanently.

**Cumulative dose toxicity** means that heart damage from anthracyclines accumulates with each administration. The risk of clinical heart failure remains low below 300 mg/m^2 of doxorubicin-equivalent dosing (approximately 1-2%), but rises sharply at higher cumulative doses — reaching 5% at 400 mg/m^2, 26% at 550 mg/m^2, and 48% at 700 mg/m^2. For pediatric patients, this cumulative risk extends across their entire lifetime: subclinical damage from childhood treatment may not manifest as heart failure until decades later, making prevention during treatment far more effective than treatment after the fact.

**Dexrazoxane** is the only FDA-approved cardioprotective agent for anthracycline chemotherapy. It is an EDTA analog that chelates intracellular iron, preventing the iron-mediated Fenton reaction that generates the most damaging free radicals. Current guidelines from the Children's Oncology Group (COG) recommend dexrazoxane when cumulative anthracycline doses exceed 300 mg/m^2. It is administered as an intravenous infusion 15-30 minutes before each anthracycline dose, at a ratio of 10:1 (dexrazoxane to doxorubicin).

**Left Ventricular Ejection Fraction (LVEF)** is the standard measure of heart pump function — the percentage of blood in the left ventricle that is ejected with each heartbeat. Normal LVEF is 55% or greater. An LVEF below 50% is concerning and may warrant dose modification or cardiology consultation. In pediatric oncology, LVEF is monitored by echocardiography before starting anthracyclines, at regular intervals during treatment, and throughout long-term follow-up. However, LVEF is a relatively insensitive marker — by the time it drops, significant myocardial damage has already occurred.

**TPMT (Thiopurine S-methyltransferase)** is a metabolic enzyme responsible for inactivating 6-mercaptopurine (6-MP), a critical maintenance therapy drug in ALL treatment. Genetic variants in the TPMT gene reduce or eliminate enzyme activity, causing 6-MP to accumulate to toxic levels. Approximately 10% of the population carries one reduced-function allele, and 0.3% carry two — these individuals can develop life-threatening myelosuppression (bone marrow failure) at standard 6-MP doses.

**TPMT phenotypes** guide 6-MP dosing decisions. Normal metabolizers (genotype *1/*1, approximately 89% of population) receive standard doses. Intermediate metabolizers (*1/*3A or *1/*3C, approximately 10%) start at 50% of the standard dose and titrate based on blood counts. Poor metabolizers (*3A/*3A, *3A/*3C, or *3C/*3C, approximately 0.3%) start at just 10% of the standard dose. Failure to test for TPMT before starting 6-MP can result in severe, potentially fatal bone marrow suppression.

**Vincristine neuropathy** is a dose-limiting toxicity of vincristine, a vinca alkaloid used in virtually every pediatric cancer protocol. Vincristine binds to tubulin and disrupts microtubule assembly, which not only kills cancer cells but also damages peripheral nerves that depend on microtubule-based axonal transport. Symptoms include numbness and tingling in hands and feet, foot drop (difficulty lifting the front of the foot), jaw pain, and constipation (from autonomic nerve involvement). Vincristine doses are capped at an absolute maximum of 2 mg per dose regardless of body surface area.

**Pharmacogenomics** is the science of using genetic information to predict drug response and toxicity. It moves beyond trial-and-error prescribing to a model where the patient's genotype informs drug selection and dosing from the first prescription. In pediatric oncology, pharmacogenomic testing for genes like TPMT, NUDT15, CYP3A5, and DPYD can prevent life-threatening toxicities. The Clinical Pharmacogenetics Implementation Consortium (CPIC) publishes evidence-based guidelines translating genotype results into specific dosing recommendations.

**Long-term survivorship** is a triumph and a challenge of modern pediatric oncology. Approximately 85% of children diagnosed with cancer today will survive, creating a growing population of long-term survivors. However, 60-70% of these survivors develop at least one chronic health condition by age 40, and 30-40% develop severe or life-threatening conditions — including heart failure, secondary cancers, endocrine disorders, and neurocognitive deficits. Prevention of treatment-related toxicity during active therapy is far more effective than treating complications decades later.

### Advanced Learning — Cardiotoxicity Prevention

**Anthracycline cardiotoxicity at the molecular level** involves a mechanism distinct from the drug's anti-tumor action. While anthracyclines kill cancer cells primarily through topoisomerase 2-alpha (Top2a) inhibition and DNA intercalation, cardiac damage occurs through topoisomerase 2-beta (Top2b) inhibition in cardiomyocytes. Top2b is the predominant topoisomerase isoform in post-mitotic cardiomyocytes, and its inhibition triggers mitochondrial dysfunction, reactive oxygen species generation, and activation of p53-mediated apoptotic pathways. This mechanistic distinction explains why cardiotoxicity cannot be fully prevented by reducing the anti-tumor dose.

**Dexrazoxane's mechanism** involves conversion from its prodrug form (ICRF-187) to its active, ring-opened form, which is an EDTA analog. Once inside the cell, it chelates free and loosely bound intracellular iron (Fe2+), preventing iron from participating in the Fenton reaction (Fe2+ + H2O2 -> Fe3+ + OH- + OH radical). The hydroxyl radicals generated by this reaction are the most reactive and damaging of all oxygen species. By removing iron from the equation, dexrazoxane breaks the catalytic cycle that converts anthracycline-generated superoxide into the hydroxyl radicals that directly damage mitochondrial membranes and myofilaments.

**CBR3 V244M polymorphism** affects the metabolism of doxorubicin into its cardiotoxic metabolite doxorubicinol. The carbonyl reductase 3 (CBR3) enzyme converts doxorubicin to doxorubicinol via carbonyl reduction, and the V244M variant (rs1056892) significantly increases this conversion rate. Patients homozygous for the Val allele (Val/Val) have approximately 3-fold higher doxorubicinol exposure and correspondingly higher cardiotoxicity risk. This polymorphism, identified through pharmacogenomic studies in childhood cancer survivors, can be tested prospectively to identify patients who may need enhanced cardioprotection or alternative anthracycline regimens.

**RARG rs2229774** is a variant in the retinoic acid receptor gamma gene that was identified through a genome-wide association study by the Canadian Pharmacogenomics Network for Drug Safety (CPNDS). The risk allele is associated with a 5-fold increased risk of anthracycline-induced cardiotoxicity in pediatric ALL patients. The proposed mechanism involves altered RARG-mediated repression of Top2b expression in cardiomyocytes — carriers of the risk allele may have higher Top2b levels, making their heart cells more vulnerable to anthracycline-induced damage. This pharmacogenomic marker is not yet part of standard pre-treatment testing but is being evaluated for clinical implementation.

**TPMT pharmacogenetics in detail** reveals that the *3A allele contains two single nucleotide polymorphisms: G238C (Ala80Pro) in exon 7 and A719G (Tyr240Cys) in exon 10. Both substitutions destabilize the protein through enhanced proteasomal degradation, with the combination reducing enzyme activity to less than 5% of wild-type. The *3C allele carries only A719G and is the most common reduced-function allele in East Asian populations (frequency 2-3%), while *3A is most common in European populations (frequency 3-5%). Understanding the population-specific allele frequencies is essential for interpreting pharmacogenomic results in diverse pediatric populations.

**NUDT15 R139C** is a pharmacogenomic variant in the nudix hydrolase 15 gene that has become essential for thiopurine dosing, particularly in patients of East Asian ancestry. NUDT15 dephosphorylates thiopurine active metabolites (TGTP and TdGTP), and the R139C variant (*3 allele, rs116855232) dramatically reduces this protective activity. Homozygous carriers (*3/*3, C/C genotype) develop severe myelosuppression even at 10% of the standard 6-MP dose. The *3 allele frequency is approximately 10% in East Asian populations, 4% in Hispanic populations, and less than 1% in European populations, making NUDT15 testing particularly important for ensuring equitable pharmacogenomic care.

**CYP3A5 and vincristine metabolism** highlights a pharmacogenomic interaction that affects efficacy rather than toxicity. CYP3A5 is the primary enzyme metabolizing vincristine, and the *1 allele confers CYP3A5 expression (the "expresser" phenotype). CYP3A5*1 carriers metabolize vincristine faster, potentially reducing drug exposure and efficacy — these patients may need dose monitoring to ensure adequate drug levels. The *3/*3 genotype (non-expresser, most common in European populations at ~80%) results in standard vincristine metabolism. CYP3A5*1 frequency varies dramatically by ancestry: approximately 70% in African populations, 30% in East Asian populations, and 15% in European populations.

**Echocardiographic strain imaging** using Global Longitudinal Strain (GLS) is emerging as a more sensitive marker of subclinical cardiotoxicity than traditional LVEF. GLS measures the percentage of myocardial fiber shortening during systole using speckle tracking technology. A GLS decline of more than 15% from baseline (or an absolute value worse than -18%) can detect myocardial damage months before LVEF shows any decline. In pediatric oncology, serial GLS monitoring during anthracycline therapy provides an early warning window during which cardioprotective interventions can be initiated or intensified.

**Children's Oncology Group (COG) Long-Term Follow-Up Guidelines** provide evidence-based recommendations for lifetime cardiac surveillance based on cumulative anthracycline exposure and other risk factors. Patients who received less than 250 mg/m^2 without chest radiation undergo echocardiography every 5 years. Those who received 250 mg/m^2 or more, or any dose combined with chest radiation, require echocardiography every 2 years, with increased frequency if abnormalities are detected. These guidelines reflect the understanding that anthracycline cardiotoxicity is a lifelong risk — new heart failure diagnoses continue to accrue for decades after treatment completion, with a cumulative incidence that does not plateau.

---

### Step 1: PGx — Drug-Gene Screening

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** Start with pharmacogenomics. This is the "safety first" message.
> Before giving any drug, check the patient's genetics for dangerous interactions.

**Actions:**

1. Switch to the **PGx Agent** tab (http://localhost:8507)
2. Note the drug interaction tabs with phenotype interpretation
3. Query about TPMT genotyping:

```
TPMT genotyping for 6-mercaptopurine in a 6-year-old male with neuroblastoma.
Patient genotype: TPMT *1/*3A. What dose adjustment is required?
```

4. Submit the query

**What you will see:** Pharmacogenomic assessment showing:

- **Genotype:** TPMT *1/*3A
- **Phenotype:** Intermediate Metabolizer
- **Clinical recommendation:** 50% dose reduction for 6-mercaptopurine and 6-thioguanine
- **Risk without adjustment:** Severe myelosuppression, potentially fatal pancytopenia
- **CPIC guideline level:** Strong recommendation
- **Alternative agents:** If dose reduction insufficient

> **PRESENTER NOTE:** This is a life-or-death finding. Full-dose 6-mercaptopurine in a TPMT
> intermediate metabolizer can cause fatal bone marrow suppression. The PGx agent catches
> this before the first dose.

**Say:**

> "Before Marcus receives a single dose of chemotherapy, we check his pharmacogenomics. His TPMT genotype is *1/*3A — intermediate metabolizer. If we gave him full-dose 6-mercaptopurine, it could destroy his bone marrow. The PGx agent catches this and recommends a 50% dose reduction. This is a life-saving intervention that happens before treatment even begins."

> 📸 **SCREENSHOT PLACEHOLDER:** [PGx Agent at http://localhost:8507 showing TPMT *1/*3A genotype result, Intermediate Metabolizer phenotype designation, 50% dose reduction recommendation for 6-mercaptopurine, and CPIC guideline level indicator]

> ⚠️ **IF THIS FAILS:** If the PGx agent is unresponsive, check: `curl http://localhost:8507/health`. If the query returns a generic response without TPMT-specific data, try rephrasing: "TPMT *1/*3A intermediate metabolizer 6-MP dosing CPIC guidelines". If the agent is completely down, state the finding verbally — this is a critical safety message that must be delivered regardless: "TPMT intermediate metabolizers need a 50% dose reduction to prevent fatal myelosuppression."

> 🎤 **SAY:** "This is pharmacogenomics saving a life before treatment even starts. Marcus has a TPMT variant that means his body cannot break down 6-mercaptopurine at the normal rate. A full dose would accumulate to toxic levels and destroy his bone marrow. The PGx agent caught this and prescribed a 50% dose reduction — following the same CPIC guidelines that the world's best pharmacists use, but in seconds instead of hours."

---

### Step 2: Biomarker — Neuroblastoma Panel

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Biomarker Agent** tab (http://localhost:8528)
2. Click the **"Biomarker Analysis"** tab or **"Evidence Explorer"** tab
3. Enter query:

```
Neuroblastoma prognostic biomarkers for 6-year-old with MYCN amplified stage 4
disease. Include LDH, ferritin, and INSS/INRGSS classification.
```

4. Submit the query

**What you will see:** Comprehensive neuroblastoma biomarker panel:

- **MYCN status:** Amplified (>10 copies) — high risk
- **LDH:** Elevated — adverse prognostic factor
- **Ferritin:** Elevated — adverse prognostic factor
- **Risk classification:** High risk per COG/INRGSS
- **Segmental chromosomal aberrations:** 1p deletion, 17q gain
- **Expected 5-year survival:** 40-50% with intensive multimodal therapy

**Say:**

> "MYCN amplification is the single worst prognostic factor in neuroblastoma. Marcus is high risk. He needs aggressive, multimodal therapy — but every drug we give him carries toxicity. Let's see what the Cardiology agent says about his heart."

> 📸 **SCREENSHOT PLACEHOLDER:** [Biomarker Agent at http://localhost:8528 showing neuroblastoma prognostic panel with MYCN amplified status, elevated LDH and ferritin, high-risk COG/INRGSS classification, segmental chromosomal aberrations, and 5-year survival estimate]

> ⚠️ **IF THIS FAILS:** If the Biomarker agent does not return neuroblastoma-specific results, check: `curl http://localhost:8528/health`. Try the Evidence Explorer tab with a simpler query: "MYCN amplified neuroblastoma prognosis high risk". If the agent is down, deliver the key message verbally: "MYCN amplification puts Marcus in the high-risk category with a 40-50% five-year survival. He needs aggressive therapy, but every drug has a cost."

> 🎤 **SAY:** "MYCN amplification is the single most important prognostic factor in neuroblastoma. It tells us Marcus's cancer is aggressive, and he needs the most intensive therapy available. But here is the dilemma — the drugs that will save his life can also damage his heart, his nerves, and his immune system. The next four agents are going to help us thread that needle."

---

### Step 3: Cardiology — Cardiotoxicity Assessment

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** The Cardiology agent has 10 tabs. Focus on Cardio-Oncology for the
> demo — it is purpose-built for this exact scenario.

**Actions:**

1. Switch to the **Cardiology Agent** tab (http://localhost:8536)
2. Note the 10 tabs: Dashboard, Clinical Query, Risk Calculator, Heart Failure, CAD Assessment, Arrhythmia, Imaging, Cardio-Oncology, Evidence Explorer, Report Generator
3. Click the **"Cardio-Oncology"** tab
4. Enter query:

```
Anthracycline cardiotoxicity risk assessment for a 6-year-old male planned for
cumulative daunorubicin 300 mg/m2 as part of neuroblastoma induction. No prior
cardiac history. Baseline LVEF 65%.
```

5. Submit the query

**What you will see:** Cardiotoxicity assessment including:

- **Cumulative dose risk:** 300 mg/m2 approaches the pediatric threshold (lifetime max ~450 mg/m2)
- **Dexrazoxane recommendation:** Cardioprotectant recommended at cumulative doses >150 mg/m2
- **Monitoring schedule:** Echocardiogram at baseline, after 150 mg/m2, after 300 mg/m2, then every 100 mg/m2
- **Baseline echo requirements:** LVEF, GLS (global longitudinal strain), troponin-I
- **Long-term surveillance:** Annual cardiac monitoring for life (anthracycline survivors)
- **Risk factors:** Young age increases lifetime risk of late cardiotoxicity

> **PRESENTER NOTE:** You can also click the **"Risk Calculators"** tab to show a structured
> cardiac risk calculator form, if the audience is interested in the tool's depth.

**Say:**

> "Marcus is 6 years old. If we cure his cancer, he has 70+ years of life ahead. But anthracyclines damage the heart, and the damage is cumulative and permanent. The Cardiology agent recommends dexrazoxane cardioprotection starting now, serial echocardiograms, and — critically — lifelong cardiac surveillance. We are protecting his future heart while fighting his current cancer."

> 📸 **SCREENSHOT PLACEHOLDER:** [Cardiology Agent Cardio-Oncology tab at http://localhost:8536 showing anthracycline cardiotoxicity risk assessment with cumulative dose threshold, dexrazoxane recommendation, echocardiogram monitoring schedule, baseline LVEF 65%, and lifelong surveillance plan]

> ⚠️ **IF THIS FAILS:** If the Cardiology agent is unresponsive, check: `curl http://localhost:8536/health`. If the Cardio-Oncology tab is not available, try the Evidence Explorer tab or the Risk Calculator tab with the same clinical parameters. If the agent is completely down, the key numbers to state verbally are: "300 mg/m2 cumulative anthracycline approaching the pediatric threshold, dexrazoxane cardioprotection indicated, echocardiograms at baseline and every 100-150 mg/m2 thereafter."

> 🎤 **SAY:** "Marcus is six years old. If we cure his neuroblastoma — and that is the goal — he has 70-plus years of life ahead of him. Anthracyclines are essential to his treatment, but they permanently damage heart muscle cells that cannot regenerate. The Cardiology agent is prescribing dexrazoxane to protect his heart during treatment, serial echocardiograms to catch damage early, and lifelong cardiac surveillance. We are not just fighting his cancer today — we are protecting his heart for the next seven decades."

---

### Step 4: Neurology — Vincristine Neuropathy

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Neurology Agent** tab (http://localhost:8529)
2. Note the 5 tabs: Dashboard, Evidence Explorer, Clinical Scales (10 scales with dynamic inputs), Workflow Runner (8 workflows), Reports & Export
3. Click the **"Evidence Explorer"** tab
4. Select the **Domain Focus** if available (choose "neuromuscular" or "general")
5. Type the following query:

```
Vincristine neuropathy risk in 6-year-old with neuroblastoma. What monitoring
is recommended? What is the dose cap?
```

6. Click the search/ask button

**What you will see:** Vincristine neuropathy assessment including:

- **Incidence:** 30-40% of pediatric patients develop peripheral neuropathy
- **Dose cap:** 2 mg absolute maximum per dose (regardless of body surface area)
- **Risk factors:** Young age, concurrent azole antifungals (CYP3A4 inhibition increases vincristine levels)
- **Monitoring:** Baseline and serial neurological exams, deep tendon reflexes, grip strength
- **Dose modification criteria:** Grade 2 neuropathy = hold dose; Grade 3 = discontinue
- **Recovery:** Usually reversible within 3-6 months of discontinuation

> **PRESENTER NOTE:** If you want to show more depth, click the **"Clinical Scales"** tab
> and demonstrate one of the 10 neurological assessment scales with dynamic inputs. The
> NIHSS has 11 items, GCS has 3 — each with interactive scoring.

**Say:**

> "Vincristine causes neuropathy in 30-40% of children. The Neurology agent sets the dose cap at 2mg, flags drug interactions that could worsen toxicity, and establishes a monitoring protocol. For a 6-year-old who needs to run, play, and go to school — protecting those nerves is not optional."

> 📸 **SCREENSHOT PLACEHOLDER:** [Neurology Agent Evidence Explorer at http://localhost:8529 showing vincristine neuropathy risk assessment with 30-40% incidence, 2 mg absolute dose cap, risk factors including concurrent azole antifungals, monitoring protocol, and dose modification criteria by neuropathy grade]

> ⚠️ **IF THIS FAILS:** If the Neurology agent is unresponsive, check: `curl http://localhost:8529/health`. If the Evidence Explorer does not return results, try selecting a different Domain Focus or removing the domain filter entirely. If the agent is down, deliver the key facts verbally: "Vincristine neuropathy affects 30-40% of pediatric patients, absolute dose cap is 2 mg, and concurrent azole antifungals can dramatically worsen toxicity by inhibiting vincristine metabolism."

> 🎤 **SAY:** "Vincristine damages peripheral nerves — the nerves that let Marcus feel his fingers, lift his feet, and control his bowels. It happens in 30 to 40 percent of children. The Neurology agent caps the dose at 2 milligrams absolute, flags drug interactions that would make the neuropathy worse, and builds a monitoring protocol. For a 6-year-old who needs to run, play, write, and go to school, protecting those nerves is not optional — it is essential to his quality of life."

---

### Step 5: Autoimmune — Dinutuximab irAE Profiling

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Autoimmune Agent** tab (http://localhost:8531)
2. Note the 10 tabs: Clinical Query, Patient Analysis, Document Ingest, Diagnostic Odyssey, Autoantibody Panel, HLA Analysis, Disease Activity, Flare Prediction, Therapy Advisor, Knowledge Base
3. Click the **"Clinical Query"** tab
4. Enter query:

```
What are the immune-related adverse events of dinutuximab (anti-GD2) immunotherapy
in pediatric neuroblastoma? Include incidence rates and management protocols.
```

5. Submit the query

**What you will see:** Immune-related adverse event (irAE) profiling:

- **Neuropathic pain:** 85% incidence — the most common and distressing side effect
- **Capillary leak syndrome:** 25% incidence — can be life-threatening
- **Hypersensitivity/anaphylaxis:** 10-15% incidence
- **Fever:** Very common
- **Hypotension:** 15-20% requiring vasopressor support
- **Management protocols:** Morphine for pain, IV fluid support, anti-GD2 infusion rate adjustments
- **Pre-medication requirements:** Detailed protocol

**Say:**

> "Dinutuximab is the immunotherapy that has improved neuroblastoma survival by 20%. But 85% of children experience severe neuropathic pain during infusion. The Autoimmune agent profiles every immune-related adverse event, with incidence rates and management protocols, so the care team is prepared before the first infusion."

> 📸 **SCREENSHOT PLACEHOLDER:** [Autoimmune Agent Clinical Query at http://localhost:8531 showing dinutuximab irAE profiling with neuropathic pain (85%), capillary leak syndrome (25%), hypersensitivity (10-15%), management protocols for each, and pre-medication requirements]

> ⚠️ **IF THIS FAILS:** If the Autoimmune agent is unresponsive, check: `curl http://localhost:8531/health`. If the query does not return dinutuximab-specific data, try: "anti-GD2 immunotherapy adverse events neuroblastoma management". If the agent is completely down, deliver the key statistic verbally: "85% of children on dinutuximab experience severe neuropathic pain — the care team must be prepared with morphine protocols and infusion rate adjustments before the first dose."

> 🎤 **SAY:** "Dinutuximab has improved neuroblastoma survival by 20 percent — it is a game-changer. But 85 percent of children experience severe pain during the infusion. Twenty-five percent get capillary leak syndrome. The Autoimmune agent maps every immune-related adverse event with exact incidence rates and management protocols. The care team knows exactly what to expect and how to respond before Marcus ever receives his first dose."

---

### Step 6: Oncology — Integrated Protocol

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** This is the synthesis step. The Oncology agent takes all the safety
> data from the other 5 agents and produces a modified protocol.

**Actions:**

1. Switch to the **Oncology Agent** tab (http://localhost:8526)
2. Click the **"Evidence Explorer"** tab
3. Enter query:

```
Dose-adjusted COG ANBL1232 protocol for 6-year-old with MYCN-amplified
neuroblastoma. TPMT intermediate metabolizer requiring 50% 6-MP dose reduction.
Anthracycline cardiac risk requires dexrazoxane. Vincristine neuropathy monitoring needed.
```

4. Submit the query

**What you will see:** Integrated, modified treatment protocol showing:

- Standard COG ANBL1232 backbone with specific modifications
- PGx-guided 6-MP dose reduction (50%)
- Dexrazoxane co-administration schedule
- Vincristine dose cap and monitoring checkpoints
- Dinutuximab irAE management plan integrated into the timeline
- Overall treatment timeline with safety checkpoints

**Say:**

> "This is the power of multi-agent coordination. The Oncology agent has taken input from PGx, Cardiology, Neurology, Autoimmune, and Biomarker agents and produced an integrated protocol with every safety modification built in. Marcus gets the most aggressive effective therapy — with every toxicity anticipated and mitigated. Six agents. One protocol. Zero preventable harm."

> 📸 **SCREENSHOT PLACEHOLDER:** [Oncology Agent Evidence Explorer at http://localhost:8526 showing the integrated dose-adjusted COG ANBL1232 protocol with PGx-guided 6-MP reduction (50%), dexrazoxane schedule, vincristine dose cap, dinutuximab irAE management plan, and safety checkpoints on the treatment timeline]

> ⚠️ **IF THIS FAILS:** If the Oncology agent does not synthesize all the modifications, check: `curl http://localhost:8527/health`. If the response is incomplete, you can manually summarize the modifications from the previous five steps: "TPMT-guided 50% 6-MP reduction, dexrazoxane cardioprotection, vincristine 2mg cap with neuropathy monitoring, and dinutuximab irAE protocols." The synthesis message is what matters most — the audience has already seen each individual agent deliver its assessment.

> 🎤 **SAY:** "This is the payoff of multi-agent coordination. Six specialized agents — PGx, Biomarker, Cardiology, Neurology, Autoimmune, and Oncology — have each contributed their expertise, and the Oncology agent has woven it all together into one integrated protocol. Marcus gets the most aggressive therapy his cancer demands, with every toxicity anticipated and every dose adjusted for his specific genetics. Six agents. One protocol. Zero preventable harm. That is the standard of care this platform enables."

---

## Demo 4: "One Gene, One Family" — Rare Disease with Cancer Predisposition

### Patient: Aurora Martinez, 4 years old, Female, Hereditary Retinoblastoma

**Duration:** 12 minutes
**Story arc:** Aurora's retinoblastoma diagnosis triggers germline variant classification, lifetime surveillance planning, and sibling testing — all coordinated across four agents.

**Browser tabs needed:**
- Rare Disease (http://localhost:8544)
- Oncology (http://localhost:8526)
- Imaging (http://localhost:8525)
- Clinical Trial (http://localhost:8128)

### Foundations Learning — Rare Disease + Cancer Predisposition

**A rare disease** is defined in the United States as a condition affecting fewer than 200,000 people. While each individual rare disease is uncommon, collectively there are approximately 7,000 known rare diseases affecting an estimated 30 million Americans — roughly 1 in 10 people. Approximately 80% of rare diseases have a genetic basis, and about half of all rare disease patients are children. Many rare diseases involve cancer predisposition, developmental delay, or multi-system organ dysfunction, making early diagnosis critical for intervention and surveillance.

**The diagnostic odyssey** describes the often agonizing journey families experience seeking a diagnosis for a rare disease. On average, it takes 5-7 years and visits to 7 or more specialists before a definitive diagnosis is reached. Approximately 30% of rare disease patients never receive a molecular diagnosis at all. During this odyssey, patients undergo repeated testing, misdiagnoses, and inappropriate treatments. Genomic sequencing combined with AI-powered variant interpretation — as demonstrated in this platform — has the potential to compress this odyssey from years to hours.

**Human Phenotype Ontology (HPO)** is a standardized vocabulary of over 16,000 precisely defined clinical features (phenotypes) organized in a hierarchical structure. Each term has a unique identifier (e.g., HP:0000478 for "Abnormality of the eye"), a definition, and relationships to parent and child terms. HPO allows computers to systematically match a patient's clinical features against known disease gene associations, enabling algorithmic differential diagnosis. The structured nature of HPO eliminates the ambiguity inherent in free-text clinical descriptions.

**ACMG/AMP variant classification** is the gold standard framework for determining whether a genetic variant is disease-causing. Developed jointly by the American College of Medical Genetics and Genomics (ACMG) and the Association for Molecular Pathology (AMP), it defines five categories: Pathogenic, Likely Pathogenic, Variant of Uncertain Significance (VUS), Likely Benign, and Benign. Classification is based on multiple lines of evidence — population frequency, computational predictions, functional studies, segregation in families, and clinical databases — combined using standardized rules.

**ACMG evidence codes** are the individual pieces of evidence used to classify variants. Each code has a prefix indicating strength (PVS = very strong, PS = strong, PM = moderate, PP = supporting) and direction (P = pathogenic, B = benign). Key codes include PVS1 (null variant in a gene where loss of function is a known mechanism of disease), PM2 (absent from population databases like gnomAD), PP3 (computational evidence supports a deleterious effect), and PP4 (patient's phenotype is highly specific for a disease with a single genetic etiology). Combining these codes determines the final classification.

**Retinoblastoma** is a cancer of the retina that almost exclusively affects children under 5 years of age, with most cases diagnosed before age 3. It occurs in approximately 1 in 15,000-20,000 live births, making it the most common intraocular malignancy of childhood. Retinoblastoma presents as leukocoria (white pupillary reflex), strabismus (misaligned eyes), or vision loss. Survival rates exceed 95% in developed countries when diagnosed early, but late diagnosis can require enucleation (eye removal) or be fatal if the tumor spreads beyond the eye.

**The RB1 gene** was the first tumor suppressor gene ever discovered, identified through studies of retinoblastoma families in the 1970s-1980s. Located on chromosome 13q14, RB1 encodes the retinoblastoma protein (pRb), which controls the cell cycle by inhibiting the E2F family of transcription factors. When both copies of RB1 are inactivated (lost or mutated), cells lose a critical brake on proliferation. In hereditary retinoblastoma, one defective copy is inherited (germline mutation), and the disease develops when the remaining normal copy is lost in a retinal cell.

**Knudson's two-hit hypothesis** is a foundational concept in cancer genetics, first proposed by Alfred Knudson in 1971 based on his study of retinoblastoma. The hypothesis states that tumor suppressor genes require inactivation of both alleles (two "hits") for tumor development. In hereditary retinoblastoma, the first hit is inherited (germline) and the second is acquired (somatic). In sporadic retinoblastoma, both hits must occur somatically in the same cell — a much rarer event, explaining why sporadic cases are typically unilateral and present later. This model applies to many other tumor suppressor genes beyond RB1.

**Cascade testing** is the process of systematically testing family members for a known pathogenic genetic variant found in a proband (the first person diagnosed). For a child with hereditary retinoblastoma carrying a germline RB1 mutation, each sibling has a 50% chance of carrying the same mutation. Cascade testing allows early identification of at-risk siblings before they develop tumors, enabling prophylactic screening that can detect retinoblastoma at its earliest, most treatable stage. Siblings who test negative can be spared years of unnecessary screening.

**Trilateral retinoblastoma** refers to the association of bilateral retinoblastoma with an intracranial midline primitive neuroectodermal tumor, most commonly a pinealoblastoma in the pineal gland. This occurs in approximately 5-15% of children with hereditary bilateral retinoblastoma and carries a poor prognosis. The pineal gland contains photoreceptor-like cells that express RB1, making it susceptible to the same two-hit inactivation that causes retinal tumors. Current guidelines recommend regular brain MRI screening for all children with germline RB1 mutations to detect pinealoblastoma early.

### Advanced Learning — Rare Disease + Cancer Predisposition

**ACMG point-based scoring** provides a quantitative framework for variant classification, moving beyond the original categorical rules to a more reproducible numeric system. Under this model, evidence codes are assigned point values: PVS1 = 8 points, PS-level evidence = 4 points each, PM-level = 2 points each, and PP-level = 1 point each (with negative values for benign evidence). The classification thresholds are: Pathogenic requires 10 or more points, Likely Pathogenic requires 6-9, VUS spans -1 to 5, Likely Benign covers -6 to -2, and Benign requires -7 or fewer. This numeric approach reduces inter-laboratory classification discordance and enables computational implementation.

**RB1 penetrance** is not absolute despite being often cited as a classic "complete penetrance" gene. For bilateral retinoblastoma, penetrance is approximately 90% — meaning 10% of germline RB1 mutation carriers never develop bilateral disease. For unilateral retinoblastoma, penetrance is approximately 45%. Low-penetrance variants are particularly challenging: promoter mutations, in-frame deletions, and mosaic mutations may cause milder phenotypes or reduced penetrance that complicates genetic counseling. Mosaic carriers (where the mutation is present in only a fraction of cells) may have a negative blood test despite carrying the variant in their retinal cells.

**Second cancer risk in hereditary retinoblastoma** is substantial and lifelong. Carriers of germline RB1 mutations have dramatically elevated risks of subsequent malignancies: osteosarcoma (approximately 500-fold increased risk), soft tissue sarcoma (approximately 100-fold), and melanoma (approximately 20-fold). The risk is particularly elevated in tissues that received radiation therapy — external beam radiation of the orbit, once standard treatment, significantly increases the risk of sarcomas in the radiation field. This has driven the shift toward chemoreduction, laser therapy, and plaque brachytherapy to avoid external beam radiation whenever possible.

**Prenatal genetic testing** for RB1 mutations includes preimplantation genetic testing for monogenic disorders (PGT-M), which can be offered to parents with known RB1 mutations who are planning future pregnancies. In PGT-M, embryos created through in vitro fertilization are biopsied at the blastocyst stage and tested for the specific familial RB1 mutation. Only embryos that do not carry the mutation are selected for transfer. Prenatal diagnosis through chorionic villus sampling (CVS) or amniocentesis is also available, allowing families to prepare for immediate postnatal screening if the fetus carries the mutation.

**Functional assessment of RB1 VUS (variants of uncertain significance)** is critical because the VUS classification creates clinical uncertainty that impedes decision-making. Several laboratory assays can provide functional evidence: mini-gene splicing assays determine whether intronic or exonic variants disrupt normal mRNA splicing, RB protein stability assays measure whether missense variants produce a stable protein or cause proteasomal degradation, and cell cycle arrest functional readouts test whether the variant protein retains its ability to inhibit E2F-mediated G1/S transition. These functional data can upgrade a VUS to Likely Pathogenic or downgrade it to Likely Benign.

**PM5 criterion implementation** in this platform demonstrates automated ACMG evidence assessment. The PM5 criterion states: "Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before." Our platform implements this by checking the ClinVar database for all known pathogenic missense variants at the same codon position. If the patient's variant affects the same residue but introduces a different amino acid substitution, the system automatically scores +2 PM-level points. This computation, which would take a human analyst 10-15 minutes of database searching, executes in under a second.

---

### Step 1: Rare Disease — Patient Intake & Diagnosis

> **TIMING: 4 minutes**
> **PRESENTER NOTE:** The Rare Disease agent is structured differently from the others.
> It uses HPO (Human Phenotype Ontology) terms for standardized phenotype entry. This is
> the language of rare disease diagnosis.

**Actions:**

1. Switch to the **Rare Disease Agent** tab (http://localhost:8544)
2. Note the 5 tabs: Patient Intake, Diagnostic Dashboard, Variant Review, Therapeutic Options, Report Generator
3. Click the **"Patient Intake"** tab
4. Fill in the intake form:

| Field | Value | Notes |
|-------|-------|-------|
| **HPO Terms** | `HP:0009919 \| Retinoblastoma \| Infantile onset \| Severe` | Type or paste into the HPO field |
| **Genomic Variants** | `RB1 \| c.958C>T p.Arg320Ter \| heterozygous` | Type or paste into the variant field |
| **Age** | `4.0` | Type into the age field |
| **Sex** | `female` | Select from dropdown |

> **PRESENTER NOTE:** HP:0009919 is the HPO term for retinoblastoma. The pipe-delimited
> format may vary — enter the HPO ID and descriptors as the form allows. If the form has
> separate fields for HPO ID, term name, onset, and severity, fill each individually.

5. Click **"Submit for Diagnosis"** (or equivalent submit button)
6. Wait for processing (5-15 seconds)
7. Click the **"Diagnostic Dashboard"** tab

**What you will see:** The Diagnostic Dashboard displays:

- **Differential diagnosis list** with hereditary retinoblastoma at the top (highest probability)
- Confidence scores for each differential
- RB1 gene association highlighted
- Inheritance pattern: Autosomal Dominant (AD)
- Penetrance information: ~90% for germline RB1 mutations
- **Diagnostic odyssey risk:** Low (clear phenotype-genotype correlation)
- Recommendations for confirmatory testing
- Family screening recommendations (siblings, parents)

**Say:**

> "Aurora has a white reflex in her left eye — leukocoria. Her ophthalmologist suspects retinoblastoma. The Rare Disease agent takes her HPO-coded phenotype and germline variant and immediately identifies hereditary retinoblastoma as the top diagnosis. But more importantly, look at the family screening recommendation. This germline RB1 mutation means every sibling and future child needs testing."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Diagnostic Dashboard showing hereditary retinoblastoma as the top differential with confidence score, RB1 gene association, and the family screening recommendation panel]

> ⚠️ **IF THIS FAILS:** If the Rare Disease agent does not return a differential or shows a loading error, refresh the page at http://localhost:8544 and re-submit. If the HPO term format is rejected, try entering just "Retinoblastoma" as free text in the phenotype field. If the agent remains unresponsive, say: "The Rare Disease agent is processing Aurora's case — let me show you a pre-computed result" and describe the expected output verbally while moving to Step 2.

> 🎤 **SAY:** "This is where rare disease meets cancer predisposition. Aurora is four years old, and a routine eye exam revealed a white glow in her pupil. We entered her HPO-coded phenotype — retinoblastoma with infantile onset — along with her RB1 germline variant. In seconds, the platform identified hereditary retinoblastoma, flagged the autosomal dominant inheritance pattern, and immediately recommended cascade testing for her two older siblings. This is not just Aurora's diagnosis — it is a family-wide alert."

---

### Step 2: Rare Disease — ACMG Variant Classification

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** ACMG classification is the gold standard for variant interpretation.
> The platform automates what typically takes a molecular geneticist 30-60 minutes per variant.

**Actions:**

1. Stay on the **Rare Disease Agent** (http://localhost:8544)
2. Click the **"Variant Review"** tab
3. Fill in the classification form:

| Field | Value |
|-------|-------|
| **Gene Symbol** | `RB1` |
| **Variant HGVS** | `c.958C>T` |
| **Zygosity** | `heterozygous` |
| **Inheritance** | `AD` (Autosomal Dominant) |

4. Click **"Classify Variant"** (or equivalent button)

**What you will see:** ACMG classification result:

- **Classification: PATHOGENIC**
- **Criteria met:**
  - **PVS1:** Null variant (nonsense) in a gene where loss of function is a known mechanism of disease
  - **PM2:** Absent or extremely low frequency in population databases (gnomAD)
  - **PP4:** Patient phenotype highly specific for RB1-associated disease
- **ACMG criteria checklist:** Visual checklist showing which of the 23 implemented ACMG criteria are met
- **ClinVar concordance:** Matches existing ClinVar pathogenic classifications
- **Clinical actionability:** High — affects treatment, surveillance, and family planning

**Say:**

> "The Variant Review tab classifies Aurora's RB1 variant using ACMG criteria. c.958C>T creates a premature stop codon — a nonsense mutation that truncates the retinoblastoma protein. PVS1, PM2, PP4 — pathogenic. This classification triggers three things: a treatment plan for Aurora, a surveillance protocol for her lifetime, and genetic testing for her entire family."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Variant Review tab showing the PATHOGENIC classification badge, the ACMG criteria checklist with PVS1/PM2/PP4 highlighted, and the ClinVar concordance section]

> ⚠️ **IF THIS FAILS:** If the ACMG classifier returns an error or times out, verify the variant HGVS is entered as `c.958C>T` (case-sensitive). If the Variant Review tab is blank, click back to Patient Intake and re-submit, then return to Variant Review. As a fallback, say: "The variant classification engine is computing — for a nonsense mutation in RB1, the expected result is Pathogenic with PVS1, PM2, and PP4 criteria met" and proceed to Step 3.

> 🎤 **SAY:** "Now watch what happens when the platform classifies Aurora's RB1 variant. c.958C>T is a nonsense mutation — it creates a premature stop codon that truncates the retinoblastoma protein. The ACMG classifier scores this as Pathogenic with three evidence codes: PVS1 for a null variant in a known loss-of-function gene, PM2 because it is absent from population databases, and PP4 because Aurora's phenotype is a textbook match. This classification, which takes a molecular geneticist 30 to 60 minutes, happened in seconds."

---

### Step 3: Imaging — Trilateral Screening

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** Trilateral retinoblastoma is a rare but devastating complication —
> a pinealoblastoma (brain tumor) that occurs in children with hereditary retinoblastoma.
> The Imaging agent knows this screening protocol.

**Actions:**

1. Switch to the **Imaging Agent** tab (http://localhost:8525)
2. Note the 9 tabs: Evidence Explorer, Workflow Runner (3 pre-loaded cases), Image Gallery, Protocol Advisor, Device & AI Ecosystem, Dose Intelligence, Reports & Export, Patient 360, Benchmarks
3. Note: The sidebar contains a **9-step demo guide** — you may reference it but follow this guide for the demo
4. Click the **"Evidence Explorer"** tab
5. Enter query:

```
Trilateral retinoblastoma screening protocol for 4-year-old with germline RB1
mutation. Brain MRI schedule and pinealoblastoma surveillance.
```

6. Submit the query

**What you will see:** Screening protocol including:

- **Brain MRI every 6 months** until age 5
- **Then annually** until age 10 (some protocols to age 15)
- **Pinealoblastoma incidence:** 5-10% of hereditary retinoblastoma patients
- **Median age of presentation:** 23 months (but can occur later)
- **MRI protocol specifics:** With gadolinium, thin-section through pineal region
- **Screening reduces mortality** by enabling early detection

**Say:**

> "Children with germline RB1 mutations don't just get retinoblastoma in their eyes. Five to ten percent develop pinealoblastoma — a brain tumor in the pineal gland. It is called trilateral retinoblastoma. The Imaging agent automatically prescribes brain MRI every 6 months until age 5. This surveillance protocol saves lives by catching brain tumors early."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Evidence Explorer response showing the trilateral retinoblastoma screening protocol, brain MRI schedule (every 6 months until age 5, then annually), and pinealoblastoma incidence statistics]

> ⚠️ **IF THIS FAILS:** If the Imaging agent returns a generic response without the trilateral protocol, rephrase the query to: "Pinealoblastoma screening schedule for germline RB1 mutation carrier age 4." If the agent is down entirely, say: "The Imaging agent coordinates surveillance protocols — for Aurora, that means brain MRI every 6 months to screen for pinealoblastoma, a brain tumor that occurs in 5-10% of hereditary retinoblastoma patients" and move to Step 4.

> 🎤 **SAY:** "Here is something most people outside pediatric oncology do not know. Children with hereditary retinoblastoma are not just at risk for eye tumors. Five to ten percent develop a brain tumor called pinealoblastoma — in the pineal gland, which has photoreceptor-like cells that express the same RB1 gene. The Imaging agent knows this and automatically prescribes brain MRI surveillance every 6 months. Aurora is four years old. Without this protocol, a pinealoblastoma could go undetected until it is too late."

---

### Step 4: Clinical Trial — Retinoblastoma Trials

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Clinical Trial Agent** tab (http://localhost:8128)
2. Navigate to **"Patient-Trial Matching"** or **"Evidence Explorer"** tab
3. Enter query:

```
Retinoblastoma clinical trials for a 4-year-old female with germline RB1 mutation
(c.958C>T). Include eye-preserving therapies and systemic chemotherapy trials.
```

4. Submit the query

**What you will see:** Matched clinical trials:

- Eye-preserving therapy trials (intra-arterial chemotherapy, intravitreal melphalan)
- Systemic chemoreduction trials
- Novel targeted therapy trials
- Long-term survivor studies

**Say:**

> "The Clinical Trial agent matches Aurora to both treatment trials and long-term surveillance studies. For retinoblastoma, the goal is not just survival — it is preserving vision whenever possible."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Patient-Trial Matching results showing matched retinoblastoma trials, including eye-preserving therapy trials (intra-arterial chemotherapy, intravitreal melphalan) and their eligibility scores]

> ⚠️ **IF THIS FAILS:** If the Clinical Trial agent returns no matches, broaden the query by removing the specific mutation detail and searching for "retinoblastoma clinical trials pediatric." If the agent is unresponsive, say: "The Clinical Trial agent is matching Aurora to active retinoblastoma trials — the priority is eye-preserving therapies like intra-arterial chemotherapy and intravitreal melphalan, which can treat the tumor without removing the eye" and proceed to the summary.

> 🎤 **SAY:** "For a four-year-old with retinoblastoma, survival is not the only goal — we want to save her eyes. The Clinical Trial agent matches Aurora to eye-preserving therapy trials, including intra-arterial chemotherapy delivered directly to the retinal blood supply, and intravitreal melphalan injected into the eye itself. These approaches can treat the tumor without enucleation. The platform also matches her to long-term survivor studies, because Aurora will need monitoring for secondary cancers for the rest of her life."

---

### Step 5: Demo 4 Summary

> **TIMING: 1 minute**

**Say:**

> "Aurora's retinoblastoma diagnosis triggered automatic germline classification, ACMG pathogenic designation, trilateral screening protocol, and family testing recommendations — all coordinated across four agents. Her siblings will be tested. Her brain will be screened. Her variant is classified. One diagnosis, rippling outward to protect an entire family."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Landing Page or the last agent screen visible, suitable for a transition slide showing the full workflow Aurora went through across four agents]

> ⚠️ **IF THIS FAILS:** This is a verbal summary step with no live agent interaction, so there is no failure mode. If previous steps had issues, use this moment to verbally recap what the platform would have shown and re-establish confidence before moving to Demo 5.

> 🎤 **SAY:** "Let me tie Aurora's story together. One four-year-old girl. One germline RB1 mutation. And from that single finding, the platform cascaded outward: ACMG pathogenic classification, trilateral retinoblastoma brain surveillance, eye-preserving clinical trials, and — critically — a recommendation to test her siblings. Her older brother and sister each have a 50% chance of carrying this same mutation. Without cascade testing, one of them could develop retinoblastoma before anyone thinks to look. This platform does not just treat the patient in front of you — it protects the family behind them."

---

## Demo 5: "Last Line of Defense" — CAR-T Therapy Decision and Monitoring

### Patient: Ethan Rodriguez, 12 years old, Male, Relapsed/Refractory B-ALL

**Duration:** 15 minutes
**Story arc:** Ethan has relapsed after two lines of chemotherapy. The team must decide: Is CAR-T the right path? Six agents evaluate eligibility, safety, and alternatives.

**Browser tabs needed:**
- Single-Cell (http://localhost:8130)
- CAR-T (http://localhost:8521)
- Cardiology (http://localhost:8536)
- Autoimmune (http://localhost:8531)
- PGx (http://localhost:8507)
- Clinical Trial (http://localhost:8128)

### Foundations Learning — CAR-T Decision

**How CAR-T works** involves a multi-step process that transforms a patient's own immune cells into precision cancer-killing weapons. First, T cells are collected from the patient's blood through a process called leukapheresis (typically a 3-4 hour blood draw). These T cells are shipped to a specialized manufacturing facility where they are genetically modified using a viral vector to express a chimeric antigen receptor (CAR) — a synthetic protein that combines an antibody's ability to recognize a specific cancer antigen with a T cell's ability to kill. The modified cells are then expanded in culture to hundreds of millions or billions and shipped back to the hospital for infusion.

**Manufacturing timeline** for CAR-T is one of its most significant logistical challenges. From the day of leukapheresis to the day of infusion, the typical turnaround is approximately 22 days (range 17-30 days depending on the manufacturer and product). During this waiting period, the patient's leukemia continues to grow, so "bridging therapy" — usually low-intensity chemotherapy or steroids — is given to control disease without damaging the precious T cells that were just collected. Vein-to-vein time management is a critical aspect of CAR-T treatment planning.

**Cytokine Release Syndrome (CRS)** is the most significant acute toxicity of CAR-T therapy. When the infused CAR-T cells encounter and bind their target antigen on tumor cells, they activate and release massive quantities of inflammatory cytokines — primarily interleukin-6 (IL-6), interferon-gamma (IFN-gamma), and tumor necrosis factor-alpha (TNF-alpha). This triggers a systemic inflammatory cascade that resembles sepsis: fever (often the first sign, typically within 1-7 days post-infusion), hypotension, capillary leak leading to edema, and in severe cases, multi-organ dysfunction. The severity of CRS generally correlates with tumor burden.

**CRS grading** follows a standardized scale that guides clinical management. Grade 1 is defined by fever alone (temperature 38 degrees C or higher) without hemodynamic compromise. Grade 2 involves hypotension that responds to intravenous fluids, or hypoxia requiring low-flow nasal cannula. Grade 3 involves hypotension requiring vasopressor medications (drugs that raise blood pressure), or hypoxia requiring high-flow oxygen or positive pressure ventilation. Grade 4 indicates life-threatening hemodynamic instability requiring multiple vasopressors, or respiratory failure requiring mechanical ventilation. Each grade has specific management algorithms.

**Tocilizumab** is a monoclonal antibody that blocks the interleukin-6 (IL-6) receptor, and it is the first-line treatment for CRS. Because IL-6 is the central cytokine driving CRS pathophysiology, blocking its receptor rapidly resolves symptoms in most patients. Dosing is weight-based: 12 mg/kg for patients weighing less than 30 kg, and 8 mg/kg for patients at 30 kg or above. Tocilizumab can be repeated every 8 hours if needed, up to a maximum of 4 doses. Importantly, tocilizumab does not appear to reduce CAR-T efficacy, making it safe to use without compromising the anti-cancer response.

**ICANS (Immune effector Cell-Associated Neurotoxicity Syndrome)** is the second most common serious adverse effect of CAR-T therapy. It presents with a constellation of neurological symptoms including confusion, word-finding difficulty, aphasia (inability to speak), tremor, and in severe cases seizures or cerebral edema. ICANS is thought to result from cytokine-mediated endothelial activation in the brain and disruption of the blood-brain barrier. Unlike CRS, ICANS does not respond well to tocilizumab (which may actually worsen it by transiently increasing CNS IL-6 levels). The primary treatment is dexamethasone, a corticosteroid that crosses the blood-brain barrier.

**CD19 antigen escape** is the primary mechanism of relapse after CD19-directed CAR-T therapy. Under the selective pressure of CAR-T cells that kill any cell expressing CD19, leukemia cells that have lost or downregulated CD19 expression gain a survival advantage and expand. This "antigen escape" results in a CD19-negative relapse that is resistant to further CD19-directed therapy. Antigen escape occurs in approximately 10-20% of patients who relapse after CD19 CAR-T, and it represents a fundamental limitation of targeting a single antigen.

**A bispecific approach** addresses antigen escape by targeting two different antigens simultaneously — most commonly CD19 and CD22 in B-ALL. By requiring a leukemia cell to lose both antigens to escape, the probability of resistance drops dramatically (the product of two independent escape probabilities rather than a single one). Bispecific CAR-T constructs encode two different antigen-binding domains on a single T cell, and clinical trials of CD19/CD22 dual-targeting CAR-T have shown promising results in reducing antigen-negative relapse rates.

**Lymphodepletion** is the chemotherapy regimen administered in the days immediately before CAR-T infusion, typically consisting of fludarabine (30 mg/m^2 daily for 4 days) and cyclophosphamide (500 mg/m^2 daily for 2 days). Lymphodepletion serves multiple purposes: it reduces the patient's existing lymphocytes to "make room" for the infused CAR-T cells, eliminates regulatory T cells that would suppress CAR-T expansion, and creates a cytokine milieu (elevated IL-7 and IL-15) that promotes CAR-T proliferation. Adequate lymphodepletion is strongly associated with better CAR-T expansion and clinical response.

**B-cell aplasia** is an expected, on-target side effect of CD19 CAR-T therapy. Because CD19 is expressed on all normal B cells (not just leukemia cells), functioning CAR-T cells will eliminate the entire B-cell compartment — healthy and malignant alike. Without B cells, the patient cannot produce new antibodies, creating a state of humoral immunodeficiency. This is managed with regular intravenous immunoglobulin (IVIG) infusions, typically monthly, for as long as the CAR-T cells persist. Paradoxically, B-cell aplasia serves as a useful pharmacodynamic marker: its persistence indicates ongoing CAR-T cell function.

### Advanced Learning — CAR-T Decision

**CAR construct generations** represent the evolution of chimeric antigen receptor engineering. First-generation CARs contained only the CD3-zeta signaling domain, which provided T-cell activation but poor persistence. Second-generation CARs — the current clinical standard — add a costimulatory domain (either 4-1BB or CD28) that provides survival and proliferation signals. Third-generation CARs include two costimulatory domains (e.g., CD28 + 4-1BB) for enhanced signaling, though clinical superiority over second-generation has not been definitively demonstrated. Fourth-generation "armored" CARs are engineered to additionally secrete cytokines (IL-12, IL-15) or express checkpoint-blocking molecules to overcome the immunosuppressive tumor microenvironment.

**4-1BB versus CD28 costimulation** creates CAR-T products with fundamentally different pharmacokinetic and clinical profiles. 4-1BB costimulation (used in tisagenlecleucel/Kymriah) activates NF-kB and promotes oxidative metabolism, resulting in slower initial T-cell expansion, enhanced central memory differentiation, and longer in vivo persistence (months to years). CD28 costimulation (used in axicabtagene ciloleucel/Yescarta) activates PI3K/AKT and promotes glycolytic metabolism, producing faster expansion, higher peak levels, more CRS, but potentially shorter persistence. The choice of costimulatory domain fundamentally shapes the therapeutic window and toxicity profile.

**CAR-T pharmacokinetics** differ from traditional drug pharmacokinetics because the "drug" is a living, replicating cell population. Key PK parameters include Cmax (peak expansion, typically occurring day 7-14 post-infusion), AUC (total exposure, which integrates expansion and persistence), and what some researchers call BIM (a decay rate constant reflecting the rate at which CAR-T cells decline after peak expansion). Higher Cmax and AUC generally correlate with better response rates, while slower BIM (longer persistence) predicts durable remission. These parameters are measured by quantitative PCR for the CAR transgene in peripheral blood.

**Single-cell CAR-T profiling** using scRNA-seq of the infusion product has revealed critical quality attributes that predict clinical outcomes. Exhaustion markers — including PD-1 (PDCD1), LAG-3 (LAG3), TIM-3 (HAVCR2), and TIGIT — when highly expressed on infusion product T cells, strongly correlate with poor expansion and treatment failure. Conversely, a naive/stem cell memory phenotype (TCF7+, LEF1+, CCR7+) in the infusion product predicts robust expansion and durable response. These single-cell biomarkers are being developed as release criteria for next-generation CAR-T manufacturing, potentially enabling quality-based go/no-go decisions before infusion.

**CD22 biology** is more complex than CD19 and creates unique challenges for CAR-T targeting. CD22 is a sialic acid-binding immunoglobulin-like lectin (Siglec) that functions as an inhibitory co-receptor on B cells. Upon antibody binding, CD22 is rapidly internalized from the cell surface — a property that can limit CAR-T engagement time. Furthermore, the efficacy of CD22 CAR-T is dose-dependent on surface antigen density: clinical data from the NIH show that leukemia cells with CD22 site density below approximately 2,000 molecules per cell are resistant to CD22 CAR-T. This density threshold is higher than for CD19 CAR-T, making target validation by flow cytometry essential before pursuing CD22-directed therapy.

**Armored CAR-T for neuroblastoma** represents a frontier approach to solid tumor CAR-T therapy, which has historically been far less successful than hematologic applications. GD2-targeting CARs for neuroblastoma face the hostile solid tumor microenvironment, which suppresses T-cell function through regulatory T cells, myeloid-derived suppressor cells, and inhibitory checkpoint ligands. Armored approaches include constitutive expression of IL-7 and CCL19 (which recruit and sustain endogenous immune cells), the C7R construct (a constitutively active IL-7 receptor that provides survival signaling without exogenous cytokine), and dominant-negative TGF-beta receptors that render CAR-T cells resistant to the immunosuppressive TGF-beta rich in the neuroblastoma TME.

**Real-world ELIANA outcomes** provide the most mature dataset for pediatric CD19 CAR-T therapy. The pivotal ELIANA trial of tisagenlecleucel in pediatric and young adult relapsed/refractory B-ALL reported an overall remission rate of 82% (all MRD-negative), with 12-month event-free survival (EFS) of 50%, 12-month overall survival (OS) of 76%, and 24-month OS of 66%. Long-term follow-up shows that patients who remain in remission beyond 12 months have a durable response plateau, with approximately 50% of initial responders maintaining remission at 3+ years. These data establish CAR-T as a transformative therapy while also highlighting that approximately half of patients eventually relapse, underscoring the need for improved approaches.

---

### Step 1: Single-Cell — CD19 Validation via CAR-T Target Workflow

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** Before starting CAR-T, we need to validate the target. The Single-Cell
> agent has a specific "CAR-T Target Validation" workflow. This is purpose-built for this
> exact clinical question.

**Actions:**

1. Switch to the **Single-Cell Agent** tab (http://localhost:8130)
2. Click the **"Workflow Runner"** tab
3. From the workflow dropdown, select **"CAR-T Target Validation"** (one of 10 available workflows)
4. Fill in the workflow parameters:

| Field | Value |
|-------|-------|
| **Target Gene** | `CD19` |
| **Tumor Type** | Select the closest leukemia/ALL option |

5. Click **"Run Workflow"**
6. Wait for results (10-30 seconds)

**What you will see:** CAR-T target validation report:

- **CD19 expression:** Confirmed positive with expression level
- **On-target/off-tumor assessment:** CD19 expressed on normal B cells (expected B-cell aplasia)
- **Antigen density analysis:** Mean fluorescence intensity (MFI) assessment
- **Escape risk:** CD19-negative subclone analysis
- **Backup targets:** CD22, CD20 expression status
- **Recommendation:** Suitable for CD19-directed CAR-T

**Say:**

> "Before we engineer CAR-T cells, we validate the target. The Single-Cell agent confirms CD19 is expressed on Ethan's leukemia cells, assesses the risk of antigen escape — where cancer cells lose CD19 to evade the CAR-T — and identifies CD22 as a backup target. The target is confirmed. Let's assess eligibility."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the CAR-T Target Validation workflow results showing CD19 expression confirmed positive, antigen density analysis, escape risk assessment, and CD22 backup target status]

> ⚠️ **IF THIS FAILS:** If the "CAR-T Target Validation" workflow is not visible in the dropdown, try selecting a general analysis workflow and entering "CD19 expression validation for B-ALL CAR-T" as the input. If the Single-Cell agent is unresponsive, restart it with `docker restart single-cell-agent` and allow 30 seconds for reload. As a fallback, say: "The Single-Cell agent validates that Ethan's blasts express CD19 at high density — over 95% positive — confirming he is a valid CAR-T candidate" and move to Step 2.

> 🎤 **SAY:** "Before we commit to manufacturing CAR-T cells for Ethan — a process that takes 22 days and costs over $400,000 — we need to confirm the target is actually there. The Single-Cell agent validates that CD19 is expressed on Ethan's leukemia cells at high density. It also checks for CD19-negative subclones that could escape the therapy, and it identifies CD22 as a backup target in case of relapse. You do not start a 22-day manufacturing process without verifying the target first."

---

### Step 2: CAR-T — Full Eligibility Assessment

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** This is the deep dive. Give a detailed, specific query to get the
> most comprehensive response.

**Actions:**

1. Switch to the **CAR-T Agent** tab (http://localhost:8521)
2. In the chat interface, type:

```
12-year-old male with relapsed/refractory B-ALL after 2 prior lines of therapy.
CD19 positive >95%, MFI 8500. Performance status ECOG 1. No active infections.
No prior CAR-T. Evaluate for tisagenlecleucel (Kymriah). What is CRS and ICANS
risk? What is the manufacturing timeline? What is the CD22 backup strategy?
```

3. Press **Enter**
4. Wait for the response with citations

**What you will see:** Comprehensive CAR-T eligibility assessment:

- **Eligibility:** ELIGIBLE for tisagenlecleucel
- **FDA-approved indication:** r/r B-ALL, age up to 25, 2+ prior lines — all met
- **CRS Risk:** Moderate (Grade 1-2: 55%, Grade 3-4: 25%)
  - High tumor burden increases CRS risk
  - Tocilizumab and dexamethasone standing protocols
- **ICANS Risk:** Lower than CRS (Grade 1-2: 20%, Grade 3-4: 10%)
- **Manufacturing timeline:** 22 days (leukapheresis to infusion)
  - Bridging therapy may be needed during manufacturing
- **CD22 backup strategy:** If CD19 relapse occurs post-CAR-T, CD22-directed therapy (inotuzumab or CD22 CAR-T trial)
- **Lymphodepletion:** Fludarabine/cyclophosphamide regimen 5 days pre-infusion
- **Citation scores** linking each recommendation to evidence

**Say:**

> "Ethan is eligible. The CAR-T agent has checked every criterion: age, prior lines, CD19 status, performance status, infection screen. The manufacturing timeline is 22 days — during which Ethan may need bridging chemotherapy to keep the leukemia in check. CRS risk is moderate. And critically, there is a backup plan: if the leukemia comes back CD19-negative, we pivot to CD22. Let's check his heart."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the CAR-T eligibility assessment showing ELIGIBLE status, CRS risk percentage (Grade 1-2: 55%, Grade 3-4: 25%), manufacturing timeline of 22 days, and the CD22 backup strategy section]

> ⚠️ **IF THIS FAILS:** If the CAR-T agent returns an incomplete response or times out, try shortening the query to: "12-year-old male, relapsed B-ALL, CD19+, evaluate for tisagenlecleucel eligibility." If the agent is down, restart with `docker restart car-t-agent`. As a fallback, say: "Ethan meets all FDA criteria for tisagenlecleucel — age under 25, two prior lines failed, CD19 confirmed. The key considerations are CRS risk management and the 22-day manufacturing window" and proceed.

> 🎤 **SAY:** "Ethan is twelve years old and his leukemia has come back twice. Two rounds of chemotherapy have failed. The CAR-T agent evaluates him against every FDA eligibility criterion for tisagenlecleucel — age, prior therapies, target expression, performance status, infection screening — and confirms he qualifies. But eligibility is just the beginning. The agent also maps out the 22-day manufacturing timeline, calculates his CRS risk at roughly 25% for severe grades, and lays out a CD22 backup strategy in case the leukemia escapes by losing CD19. This is a complete decision-support package, not just a yes-or-no answer."

---

### Step 3: Cardiology — Pre-Lymphodepletion Assessment

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Cardiology Agent** tab (http://localhost:8536)
2. Click the **"Evidence Explorer"** tab
3. Enter query:

```
Pre-CAR-T lymphodepletion cardiac assessment for 12-year-old male with prior
anthracycline exposure 250 mg/m2. Current LVEF 58%. Planned fludarabine/
cyclophosphamide lymphodepletion. Is cardiac clearance appropriate?
```

4. Submit the query

**What you will see:** Cardiac clearance assessment:

- **LVEF 58%:** Borderline (normal >55%, but reduced from baseline)
- **Prior anthracycline exposure:** 250 mg/m2 — significant cumulative dose
- **Lymphodepletion cardiac risk:** Cyclophosphamide can cause acute cardiac toxicity at high doses
- **Recommendation:** Proceed with cardiac monitoring; echocardiogram immediately pre-lymphodepletion and 7 days post-infusion
- **CRS cardiac considerations:** CRS can cause myocardial depression; cardiology should be on the CAR-T response team
- **Troponin monitoring:** During CRS episodes

**Say:**

> "Ethan's LVEF is 58% — technically normal, but reduced from what we would expect in a healthy 12-year-old. Prior anthracycline exposure has already taken a toll. The Cardiology agent clears him for CAR-T but with enhanced cardiac monitoring and a standing cardiology consult during any CRS episodes."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Cardiology Evidence Explorer response showing LVEF 58% assessment, prior anthracycline exposure flag at 250 mg/m2, cardiac clearance recommendation with enhanced monitoring, and CRS cardiac considerations]

> ⚠️ **IF THIS FAILS:** If the Cardiology agent returns a generic response without CAR-T-specific cardiac assessment, add "cyclophosphamide cardiotoxicity risk" to the query to trigger a more specific response. If the agent is unresponsive, say: "The Cardiology agent evaluates Ethan's heart before lymphodepletion — his ejection fraction is 58%, borderline after prior anthracycline exposure, cleared with enhanced monitoring" and move to Step 4.

> 🎤 **SAY:** "Before we can give Ethan CAR-T cells, he needs lymphodepletion chemotherapy — fludarabine and cyclophosphamide — to make room in his immune system. But Ethan has already received 250 milligrams per square meter of anthracyclines from his prior treatments. His heart has taken a hit. The Cardiology agent checks his ejection fraction, assesses the cumulative cardiac burden, and clears him — but with a critical condition: enhanced cardiac monitoring throughout the CAR-T process and cardiology on standby during any CRS episodes, because cytokine storm can depress the heart even further."

---

### Step 4: Autoimmune — CRS/irAE Risk Profiling

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Autoimmune Agent** tab (http://localhost:8531)
2. Click the **"Clinical Query"** tab
3. Enter query:

```
What immune-related adverse events should be monitored after CD19 CAR-T therapy in
a 12-year-old? Include CRS management, ICANS grading, and delayed autoimmune
complications including B-cell aplasia and autoimmune cytopenias.
```

4. Submit the query

**What you will see:** Comprehensive irAE profiling:

- **CRS management ladder:**
  - Grade 1: Supportive care, monitoring
  - Grade 2: Tocilizumab 8 mg/kg
  - Grade 3: Tocilizumab + dexamethasone 10 mg IV
  - Grade 4: Tocilizumab + methylprednisolone 1 g IV
- **ICANS grading:** ICE (Immune Effector Cell-Associated Encephalopathy) score
- **Delayed autoimmune complications:**
  - B-cell aplasia: 100% expected (CD19 CAR-T kills normal B cells too)
  - IVIG replacement: Monthly, indefinitely
  - Autoimmune cytopenias: 30-40% incidence (ITP, AIHA)
  - Hypogammaglobulinemia: Universal, requires monitoring
- **Monitoring schedule:** Daily for 14 days, weekly for 30 days, monthly for 1 year

**Say:**

> "CAR-T is not a one-time treatment — it is a new chronic condition. Ethan will have B-cell aplasia for months to years. He will need monthly IVIG infusions. Thirty to forty percent of patients develop autoimmune cytopenias. The Autoimmune agent maps out every immune complication and the monitoring timeline. This is the hidden cost of CAR-T — and the platform ensures nothing is missed."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Autoimmune Clinical Query response showing the CRS management ladder (Grade 1-4), ICANS grading, B-cell aplasia timeline, IVIG replacement schedule, and the 30-40% autoimmune cytopenia risk]

> ⚠️ **IF THIS FAILS:** If the Autoimmune agent returns a response focused on autoimmune disease rather than CAR-T complications, rephrase the query to emphasize "CAR-T toxicity monitoring" and "cytokine release syndrome management protocol." If the agent is down, say: "The Autoimmune agent maps every immune complication of CAR-T — CRS grading, ICANS monitoring, B-cell aplasia management, and delayed autoimmune cytopenias that affect 30-40% of patients" and continue to Step 5.

> 🎤 **SAY:** "Here is what families often do not hear in the initial CAR-T conversation. CAR-T is not a single infusion and you go home. Ethan's engineered T cells will kill every B cell in his body — cancerous and healthy alike. He will have no antibody production for months, possibly years. He will need monthly immunoglobulin infusions. Thirty to forty percent of patients develop autoimmune cytopenias — his own immune system attacking his blood cells. The Autoimmune agent lays out the complete monitoring timeline: daily for 14 days, weekly for a month, monthly for a year. This is the full picture of what CAR-T really means for a 12-year-old's life."

---

### Step 5: PGx — Supportive Care Pharmacogenomics

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **PGx Agent** tab (http://localhost:8507)
2. Enter query:

```
Tocilizumab and corticosteroid pharmacogenomics for CAR-T supportive care in a
12-year-old. Any CYP450 interactions? Does dexamethasone metabolism vary by genotype?
```

3. Submit the query

**What you will see:** PGx assessment for supportive care:

- **Tocilizumab:** No significant CYP450 interactions (monoclonal antibody, not hepatically metabolized via CYP system)
- **Dexamethasone:** Metabolized via CYP3A4; polymorphisms may affect clearance but not clinically actionable at this time
- **Fludarabine:** Minimal CYP involvement
- **Cyclophosphamide:** CYP2B6/CYP2C19 polymorphisms can affect activation
- **Key interactions to monitor:** Concurrent medications, azole antifungals (CYP3A4 inhibition)

**Say:**

> "Good news: tocilizumab — the main CRS rescue drug — has no significant pharmacogenomic interactions. It is a monoclonal antibody, cleared by the immune system, not the liver. But the PGx agent also flags cyclophosphamide metabolism variants and drug interactions to watch during the lymphodepletion phase."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the PGx agent response showing tocilizumab clearance pathway (no CYP450 involvement), dexamethasone CYP3A4 metabolism note, cyclophosphamide CYP2B6/CYP2C19 pharmacogenomic flags, and drug interaction warnings]

> ⚠️ **IF THIS FAILS:** If the PGx agent returns limited results or does not recognize tocilizumab, try a simplified query: "CYP450 drug interactions for cyclophosphamide and dexamethasone in a 12-year-old." If the agent is unresponsive, say: "The PGx agent confirms tocilizumab — our primary CRS rescue drug — has no significant pharmacogenomic interactions because it is a monoclonal antibody cleared by the immune system, not CYP enzymes" and proceed.

> 🎤 **SAY:** "When Ethan develops CRS — and with his tumor burden, there is a significant chance he will — the rescue drug is tocilizumab. The PGx agent checks whether Ethan's pharmacogenomic profile creates any complications with the drugs he will need. The good news: tocilizumab is a monoclonal antibody, so it bypasses the liver's CYP450 system entirely. No gene-drug interactions to worry about for the rescue drug. But the agent also flags cyclophosphamide — part of his lymphodepletion regimen — which does have CYP2B6 and CYP2C19 variants that can affect how fast his body activates the drug."

---

### Step 6: Clinical Trial — CAR-T Trials

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Clinical Trial Agent** tab (http://localhost:8128)
2. Navigate to **"Patient-Trial Matching"** tab
3. Enter query:

```
CAR-T clinical trials for pediatric relapsed/refractory B-ALL. 12-year-old male,
CD19 positive, failed 2 prior lines. Include both approved products and
investigational dual-targeting approaches.
```

4. Submit the query

**What you will see:** Matched trials:

- **ELIANA (NCT02435849):** Tisagenlecleucel pivotal trial (may be closed to enrollment but referenced for data)
- **ZUMA-4 (NCT02625480):** Axicabtagene ciloleucel in pediatric ALL
- **Dual CD19/CD22 CAR-T trials:** Next-generation approaches to prevent antigen escape
- **Armored CAR-T trials:** Enhanced persistence constructs
- **Allogeneic (off-the-shelf) CAR-T trials:** No manufacturing delay

**Say:**

> "The Clinical Trial agent matches Ethan to both approved CAR-T products and next-generation approaches. The dual CD19/CD22 trials are particularly interesting — they attack two targets simultaneously, reducing the chance of antigen escape. And the allogeneic trials eliminate the 22-day manufacturing wait."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Patient-Trial Matching results showing ELIANA (NCT02435849), ZUMA-4 (NCT02625480), dual CD19/CD22 CAR-T trials, and allogeneic off-the-shelf CAR-T trials with eligibility scores]

> ⚠️ **IF THIS FAILS:** If the Clinical Trial agent returns few or no matches, broaden the query to: "CAR-T clinical trials pediatric ALL relapsed refractory." If the agent times out, say: "The Clinical Trial agent matches Ethan to landmark trials like ELIANA and ZUMA-4, as well as next-generation dual-targeting CD19/CD22 approaches that reduce antigen escape risk" and proceed to the summary.

> 🎤 **SAY:** "The Clinical Trial agent matches Ethan to the full landscape of CAR-T options. ELIANA is the landmark tisagenlecleucel trial that proved CAR-T works in children — 82% remission rate. ZUMA-4 is testing an alternative CAR-T product. But the most exciting matches are the dual CD19/CD22 trials. Remember, one of the biggest risks of CAR-T is antigen escape — the leukemia loses CD19 and becomes invisible to the therapy. Dual-targeting trials attack two antigens at once, making escape exponentially harder. And the allogeneic trials — off-the-shelf CAR-T — would eliminate the 22-day manufacturing wait entirely."

---

### Demo 5 Summary

> **TIMING: 1 minute**

**Say:**

> "Ethan's CAR-T decision required six agents working in concert. Target validation. Eligibility assessment. Cardiac clearance. Immune complication profiling. Pharmacogenomic safety. Trial matching. In a traditional workflow, these consultations take weeks of scheduling. Here, they took 15 minutes. When a 12-year-old has relapsed leukemia, weeks matter."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Landing Page or last visible agent screen as a transition point, showing the breadth of the six-agent CAR-T workflow Ethan went through]

> ⚠️ **IF THIS FAILS:** This is a verbal summary step with no live agent interaction. If prior steps encountered issues, use this moment to verbally emphasize the six-agent coordination concept and re-establish the narrative before transitioning to Demo 6.

> 🎤 **SAY:** "Think about what just happened. Six different specialist agents — single-cell biology, CAR-T therapy, cardiology, autoimmune monitoring, pharmacogenomics, and clinical trials — all evaluated the same twelve-year-old boy. In a real hospital, these consultations would require scheduling six different specialists across multiple departments, waiting days or weeks for each opinion. Here, it happened in 15 minutes on one machine. For a child whose leukemia has already relapsed twice, that kind of speed is not a convenience — it is a lifeline."

---

## Demo 6: "When the Standard Drug Can't Be Used" — Novel Drug Discovery for a Growing Child

### Patient: Aiden Park, 10 years old, Male, SHH-Subtype Medulloblastoma

**Duration:** 15-20 minutes
**Story arc:** Aiden's medulloblastoma has a druggable target — but the standard drug would harm a growing child. The platform identifies the contraindication, generates novel alternatives, and matches to trials.

**Browser tabs needed:**
- Oncology (http://localhost:8526)
- Neurology (http://localhost:8529)
- Imaging (http://localhost:8525)
- Drug Discovery (http://localhost:8510)
- Clinical Trial (http://localhost:8128)

### Foundations Learning — Medulloblastoma + Drug Discovery

**Medulloblastoma** is the most common malignant brain tumor in children, accounting for approximately 20% of all pediatric central nervous system tumors. It arises in the cerebellum — the structure at the back of the brain responsible for coordination, balance, and motor control. Medulloblastoma typically presents with headache, vomiting, and unsteady gait due to increased intracranial pressure and cerebellar dysfunction. Treatment involves surgical resection, craniospinal irradiation, and chemotherapy, with overall 5-year survival rates of approximately 70-80%, though outcomes vary dramatically by molecular subtype.

**Molecular subtypes** have revolutionized medulloblastoma classification and treatment. The current WHO classification recognizes four molecular groups: WNT (approximately 10% of cases, with the best prognosis and greater than 95% survival — clinical trials are actively de-escalating therapy for this group), SHH or Sonic Hedgehog (approximately 30%, intermediate prognosis, characterized by hedgehog pathway activation), Group 3 (approximately 25%, worst prognosis, frequently MYC-amplified with high metastatic risk), and Group 4 (approximately 35%, intermediate prognosis, the most common but least understood group). Subtype determines prognosis, treatment intensity, and eligibility for targeted therapies.

**The SHH (Sonic Hedgehog) pathway** is a fundamental signaling cascade that controls cell proliferation and differentiation during embryonic development. In normal cerebellar development, SHH signaling drives the massive proliferation of granule neuron precursors — the most abundant neurons in the brain. This pathway is normally silenced after development is complete. In SHH medulloblastoma, mutations reactivate this pathway (most commonly through loss of the inhibitory receptor PTCH1 or activating mutations in SMO), causing uncontrolled proliferation of cerebellar progenitor cells that form the tumor.

**PTCH1 (Patched 1)** is the central tumor suppressor in the SHH pathway, functioning as the "brake" on pathway activation. Under normal conditions, PTCH1 sits in the cell membrane and actively inhibits SMO (Smoothened), keeping the pathway off. When SHH ligand binds to PTCH1, it relieves this inhibition, allowing SMO to activate downstream signaling. In SHH medulloblastoma, loss-of-function mutations in PTCH1 remove the brake entirely — SMO becomes constitutively active regardless of SHH ligand, driving continuous cell proliferation. PTCH1 mutations are found in approximately 40% of SHH medulloblastomas.

**Vismodegib** is an FDA-approved small molecule inhibitor of SMO, the key signal transducer in the SHH pathway. By binding directly to SMO's transmembrane domain, vismodegib blocks signal transduction even in the absence of PTCH1's inhibitory function. It has shown efficacy in SHH-driven cancers including advanced basal cell carcinoma and adult SHH medulloblastoma. However, its use in pediatric patients is severely limited by a devastating side effect: permanent growth plate fusion.

**Why children cannot take vismodegib** is rooted in developmental biology. The SHH pathway plays an essential role in growth plate chondrocyte proliferation — the process by which children's bones elongate. Growth plates (physes) are cartilaginous regions near the ends of long bones where active SHH signaling drives chondrocyte division. Blocking SMO with vismodegib in a growing child causes premature and irreversible fusion of these growth plates, permanently stunting bone growth. Case reports have documented children losing 5-10 cm of predicted adult height after even brief vismodegib exposure, making it contraindicated in skeletally immature patients.

**Posterior fossa syndrome (PFS)**, also called cerebellar mutism syndrome, is one of the most distressing complications of medulloblastoma surgery. It occurs in 25-30% of children following surgical resection of a posterior fossa (cerebellar) tumor, typically manifesting 1-5 days after surgery. Symptoms include mutism (complete inability to speak), emotional lability (inappropriate laughing or crying), hypotonia, ataxia, and cognitive changes. PFS is usually temporary, with speech recovering over 1-6 months, but many children have lasting cerebellar cognitive deficits. The exact mechanism likely involves damage to the dentato-thalamo-cortical pathway during surgery.

**Craniospinal irradiation (CSI)** is the standard radiation approach for medulloblastoma because the tumor has a propensity to disseminate through the cerebrospinal fluid to the brain and spinal cord. CSI delivers radiation to the entire craniospinal axis (brain and spinal cord) followed by a boost to the posterior fossa tumor bed. While effective at preventing leptomeningeal spread, CSI causes significant long-term toxicities in children: neurocognitive decline, growth hormone deficiency, thyroid dysfunction, hearing loss, and increased risk of secondary cancers. These late effects are most severe in younger children, creating an urgent need for approaches that reduce or eliminate CSI.

**Proton versus photon radiation** represents a critical technological advance for pediatric brain tumors. Conventional photon radiation (X-rays) deposits energy along the entire beam path, irradiating healthy tissue both before and behind the tumor. Proton beams, by contrast, deposit most of their energy at a specific depth (the Bragg peak) with minimal exit dose. For craniospinal irradiation, this physical property translates to 50-60% less radiation exposure to the heart, lungs, gut, and other organs anterior to the spine. In the brain, proton therapy reduces dose to the cochlea (hearing preservation) and hippocampus (memory preservation) by up to 80-95% compared to photons.

**The blood-brain barrier (BBB)** is a selective permeability barrier formed by tight junctions between brain capillary endothelial cells, astrocyte foot processes, and pericytes. It protects the brain from toxins and pathogens but also prevents most drugs from reaching brain tumors. For drug candidates targeting medulloblastoma, BBB penetration is a critical design constraint: molecules generally need to be small (molecular weight below 500 daltons), moderately lipophilic (logP between 1 and 3), have low polar surface area, and minimal hydrogen bond donors. These physicochemical constraints significantly narrow the chemical space available for brain tumor drug design.

### Advanced Learning — Medulloblastoma + Drug Discovery

**The SHH pathway signaling cascade** operates through a precisely choreographed series of molecular interactions. In the absence of Hedgehog ligand, PTCH1 inhibits SMO by preventing its accumulation in the primary cilium. When SHH ligand binds PTCH1, PTCH1 is internalized and degraded, allowing SMO to translocate to the primary cilium and activate. Active SMO disrupts the SUFU (Suppressor of Fused) complex that normally sequesters GLI transcription factors, releasing GLI1 and GLI2 to enter the nucleus. Nuclear GLI activates target genes including MYCN, cyclin D1, BCL2, and GLI1 itself (positive feedback). Understanding this cascade reveals multiple potential drug targets beyond SMO.

**SMO crystal structure (PDB: 5L7D)** reveals that Smoothened is a class F GPCR-like seven-transmembrane receptor with a unique pharmacology. The drug binding pocket resides within the transmembrane domain, distinct from the extracellular cysteine-rich domain (CRD) where cholesterol and oxysterols bind as endogenous modulators. The crystal structure in complex with vismodegib shows the drug occupying the same pocket as cyclopamine (the natural product SMO inhibitor from Veratrum californicum), with key contacts to residues D473, E518, and W281. This structural information is essential for rational design of next-generation inhibitors.

**Vismodegib resistance mechanisms** illustrate the challenges of single-agent targeted therapy. The most common mechanism is point mutations in SMO that disrupt drug binding — SMO D473H is the most frequently reported, mutating a key drug-contact residue that reduces vismodegib binding affinity by over 100-fold. However, resistance can also occur downstream of SMO: loss of SUFU (which normally inhibits GLI in the absence of SMO signaling) makes cells pathway-active regardless of SMO status. GLI2 amplification provides excess transcription factor that overwhelms residual SMO inhibition. These downstream mechanisms cannot be overcome by improved SMO inhibitors alone, highlighting the need for GLI-targeting strategies.

**Novel SMO inhibitor design** for pediatric patients must solve a unique pharmacological puzzle: achieving tumor-killing SHH pathway inhibition in the brain while avoiding growth plate toxicity in the skeleton. Theoretical approaches include partial agonists that reduce (but do not eliminate) pathway activity below the oncogenic threshold while maintaining enough signaling for growth plate function. Pathway-selective inhibitors could exploit tissue-specific differences in SHH signaling components. GLI-targeting compounds (such as GANT-61 or arsenic trioxide) act downstream of SMO entirely, potentially avoiding the growth plate issue since growth plate chondrocytes may be more dependent on SMO-level signaling than GLI-level signaling.

**MolMIM latent space exploration for SHH inhibitors** in this platform demonstrates how AI-driven drug discovery addresses the pediatric SMO inhibitor challenge. The generative model is conditioned on the SMO binding pocket pharmacophore — the 3D arrangement of chemical features (hydrogen bond donors/acceptors, hydrophobic contacts, aromatic stacking) required for binding. Constraints are applied for BBB penetration: molecular weight under 500 daltons, computed logP between 1 and 3, topological polar surface area under 90 square angstroms, and fewer than 3 hydrogen bond donors. Flat aromatic systems (which often have poor brain penetration due to efflux transporter recognition) are penalized in the generation objective. The resulting candidates occupy novel chemical space distinct from vismodegib.

**DiffDock scoring for SMO binding** evaluates generated candidates against the crystallographic binding pocket. The model predicts binding poses and estimates binding free energy, which is compared to the vismodegib reference value of -9.4 kcal/mol. Candidates scoring within 2 kcal/mol of vismodegib are considered promising leads. An important selectivity consideration is distinguishing SMO binding from Frizzled receptor binding — Frizzled receptors share the seven-transmembrane GPCR-like architecture and CRD domain with SMO, and off-target Frizzled inhibition could disrupt WNT signaling with unpredictable consequences. DiffDock can evaluate selectivity by counter-screening candidates against Frizzled receptor structures.

**Pediatric brain tumor radiation dosimetry** illustrates the clinical impact of proton versus photon therapy. For standard-risk medulloblastoma CSI at 23.4 Gy (with a posterior fossa boost to 54 Gy), proton therapy reduces the mean cochlear dose from approximately 33 Gy (photons) to less than 2 Gy (protons), preserving hearing in over 95% of patients compared to less than 50% with photons. Hippocampal dose is reduced by approximately 80%, which is critical because the hippocampus is the primary site of memory consolidation and adult neurogenesis. High-risk medulloblastoma requiring 36 Gy CSI sees even larger absolute dose reductions to organs at risk. These dosimetric advantages translate directly into preserved quality of life.

**Neurocognitive late effects** of craniospinal irradiation represent the most impactful long-term burden for medulloblastoma survivors. Studies consistently show an IQ decline of 2-4 points per year following CSI, with the most severe effects in children treated before age 7 — this population may lose 10-20 IQ points by adulthood. The mechanism involves radiation-induced damage to white matter tracts (particularly oligodendrocyte precursors), hippocampal neurogenesis, and cerebral vasculature. Processing speed and working memory are disproportionately affected. Proton therapy mitigates but does not eliminate these effects, as the craniospinal fields still necessarily irradiate much of the developing brain. Active research into reduced-dose CSI regimens, focal radiation approaches, and neuroprotective agents (memantine, hippocampal-avoidant techniques) aims to further reduce this burden.

---

### Step 1: Oncology — SHH Pathway Analysis

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** This is the most dramatic clinical moment in all 6 demos. Vismodegib
> is FDA-approved for SHH-driven cancers in adults. But in children, it causes permanent
> growth plate closure. The platform must catch this contraindication.

**Actions:**

1. Switch to the **Oncology Agent** tab (http://localhost:8526)
2. Click the **"Evidence Explorer"** tab
3. Type the following query:

```
10-year-old male with SHH-subtype medulloblastoma. PTCH1 nonsense mutation, TP53
wild-type, no MYCN amplification. What targeted therapy options exist? Is vismodegib
appropriate for a child?
```

4. Submit the query
5. Wait for the response (5-15 seconds)

**What you will see:** SHH pathway analysis with a CRITICAL safety flag:

- **SHH pathway:** Activated via PTCH1 loss-of-function mutation
- **Vismodegib (Erivedge):** FDA-approved SMO inhibitor for SHH-driven cancers
- **CRITICAL CONTRAINDICATION:** Vismodegib causes **irreversible premature growth plate fusion** in children
  - Documented in pediatric MATCH trial
  - Causes permanent short stature
  - Growth plates do not reopen after drug discontinuation
- **TP53 wild-type:** Favorable — TP53-mutant SHH medulloblastoma has worse prognosis
- **Alternative approaches:** Reduced-dose craniospinal radiation, novel SHH inhibitors in trials, downstream pathway targeting
- **Molecular risk:** Standard risk (TP53 wt, no MYCN amp)

> **PRESENTER NOTE:** Pause here. Let the audience absorb this. The drug that would work
> best would permanently stunt a 10-year-old's growth. This is the tension that drives
> the Drug Discovery demo in Step 4.

**Say:**

> "Look at this — CRITICAL flag. Vismodegib targets the exact pathway driving Aiden's tumor. In an adult, this would be the obvious choice. But Aiden is 10. His growth plates are open. Vismodegib would permanently fuse them. He would stop growing. Forever. The platform caught this contraindication and is telling us: you need a different approach. That is exactly what the Drug Discovery engine is for."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Evidence Explorer response showing the SHH pathway analysis, the CRITICAL CONTRAINDICATION flag for vismodegib with growth plate fusion warning, PTCH1 loss-of-function identification, and the alternative approaches section]

> ⚠️ **IF THIS FAILS:** If the Oncology agent does not flag the vismodegib contraindication, add to the query: "Specifically address growth plate toxicity of vismodegib in a skeletally immature patient." If the CRITICAL flag does not appear, verbally emphasize: "The platform identifies that vismodegib — the obvious targeted therapy — would cause irreversible growth plate fusion in a 10-year-old. This is the clinical dilemma that drives the rest of this demo." If the agent is down, this point is too important to skip — describe the contraindication in full and proceed to Step 2.

> 🎤 **SAY:** "This is the most important moment in all six demos. Aiden has an SHH-subtype medulloblastoma with a PTCH1 mutation — the Sonic Hedgehog pathway is driving his tumor. There is an FDA-approved drug for exactly this pathway: vismodegib. In an adult, you would prescribe it today. But Aiden is ten. His bones are still growing. The SHH pathway does not just drive his cancer — it drives the growth plates in every long bone in his body. Give him vismodegib, and those growth plates fuse permanently. He stops growing. Case reports show children losing 5 to 10 centimeters of predicted adult height. The platform caught this, flagged it as CRITICAL, and now we need to find another way."

---

### Step 2: Neurology — Posterior Fossa Syndrome

> **TIMING: 3 minutes**
> **PRESENTER NOTE:** Posterior fossa syndrome (cerebellar mutism) occurs in 25-30% of
> children after medulloblastoma surgery. The Neurology agent's Workflow Runner has
> pre-built workflows for brain tumor assessment.

**Actions:**

1. Switch to the **Neurology Agent** tab (http://localhost:8529)
2. Click the **"Workflow Runner"** tab (one of 8 available workflows)
3. Select the **"Brain Tumor Grading"** workflow from the dropdown
4. Fill in the workflow parameters:

| Field | Value |
|-------|-------|
| **Histology** | `medulloblastoma` |
| **Molecular subtype** | `SHH/PTCH1` |
| **Location** | `posterior fossa` |

5. Click **"Run Workflow"**

**What you will see:** Brain tumor grading and complication assessment:

- **Posterior fossa syndrome risk:** 25-30% after surgical resection
  - Cerebellar mutism: inability to speak, lasting days to months
  - Ataxia, emotional lability
  - Risk factors: midline tumor location, brainstem invasion
- **Cranial radiation neurocognitive effects:**
  - IQ decline: 2-4 points per year post-radiation
  - Executive function impairment
  - Memory and processing speed deficits
- **Proton therapy recommendation:** Reduces radiation scatter to hippocampus and temporal lobes
- **Neurocognitive monitoring:** Baseline and annual neuropsychological testing
- **Rehabilitation:** Early speech therapy, occupational therapy, educational support

**Say:**

> "Even if we cure the cancer, the treatment leaves scars. Posterior fossa syndrome affects a quarter of these children — Aiden could wake up from surgery unable to speak. And cranial radiation will cost him 2-4 IQ points per year. The Neurology agent recommends proton therapy to minimize brain damage and a comprehensive neurocognitive monitoring plan. We are not just treating cancer. We are protecting a 10-year-old's ability to think, learn, and grow up."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Brain Tumor Grading workflow results showing posterior fossa syndrome risk at 25-30%, cranial radiation neurocognitive effects (IQ decline 2-4 points/year), proton therapy recommendation, and the neurocognitive monitoring plan]

> ⚠️ **IF THIS FAILS:** If the "Brain Tumor Grading" workflow is not in the dropdown, try selecting a general workflow and entering "medulloblastoma SHH posterior fossa complication assessment." If the Neurology agent is unresponsive, restart with `docker restart neurology-agent`. As a fallback, say: "The Neurology agent maps surgical and radiation complications — 25-30% risk of posterior fossa syndrome after surgery, and 2-4 IQ points lost per year from cranial radiation. Proton therapy is recommended to minimize damage" and continue.

> 🎤 **SAY:** "Even when we cure medulloblastoma, the treatment itself exacts a devastating toll. The Neurology agent runs the Brain Tumor Grading workflow and tells us two things. First: after surgery to remove the tumor, there is a 25-30% chance Aiden wakes up unable to speak — a condition called posterior fossa syndrome or cerebellar mutism. It usually recovers over weeks to months, but some cognitive deficits are permanent. Second: the cranial radiation he needs will cost him 2-4 IQ points every single year. By the time he finishes high school, he could have lost 15-20 IQ points. The Neurology agent recommends proton therapy to reduce the damage and a lifetime neurocognitive monitoring plan. We are fighting for Aiden's future — not just his survival, but his ability to learn, to think, to grow up whole."

---

### Step 3: Imaging — Staging Protocol

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Imaging Agent** tab (http://localhost:8525)
2. Click the **"Evidence Explorer"** tab
3. Enter query:

```
Medulloblastoma staging protocol for 10-year-old male. Brain and spinal MRI
requirements. Post-operative imaging timeline. Leptomeningeal dissemination screening.
```

4. Submit the query

**What you will see:** Complete staging protocol:

- **Pre-operative brain MRI:** With and without gadolinium, thin-section posterior fossa
- **Spinal MRI:** Entire spine, with gadolinium, to detect drop metastases
  - Must be done **before surgery** (post-operative blood products obscure imaging for 2-3 weeks)
- **Post-operative brain MRI:** Within **48 hours** of surgery (residual tumor assessment)
- **CSF cytology:** Lumbar puncture at least 14 days post-surgery (avoid false positives from surgical blood)
- **Chang staging system:** M0-M4 based on dissemination
- **Leptomeningeal screening:** Critical for risk stratification (M+ = high risk regardless of molecular subtype)

**Say:**

> "Staging determines everything. The Imaging agent prescribes the exact protocol: brain and spine MRI before surgery to catch drop metastases, a 48-hour post-op scan to measure residual tumor, and CSF cytology 2 weeks later. If Aiden has leptomeningeal spread, his risk category changes completely — regardless of his favorable molecular profile."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Evidence Explorer response showing the complete staging protocol: pre-operative brain and spinal MRI requirements, 48-hour post-operative scan timeline, CSF cytology at 14 days, and the Chang staging system M0-M4 classification]

> ⚠️ **IF THIS FAILS:** If the Imaging agent returns a generic MRI protocol without medulloblastoma-specific staging, add "Chang staging system and leptomeningeal dissemination screening" to the query. If the agent is unresponsive, say: "The Imaging agent prescribes the staging protocol — pre-operative brain and spine MRI to detect drop metastases, a 48-hour post-op scan, and CSF cytology two weeks after surgery. If there is leptomeningeal spread, Aiden's risk category escalates regardless of his favorable molecular profile" and proceed.

> 🎤 **SAY:** "Before Aiden goes to surgery, the Imaging agent lays out the exact staging protocol. Brain MRI with gadolinium for the primary tumor. Full-spine MRI to look for drop metastases — tumor cells that have seeded down the spinal cord through the cerebrospinal fluid. This spinal MRI must happen before surgery, because post-operative blood products will obscure the images for two to three weeks. Then a 48-hour post-op brain scan to measure residual tumor, and CSF cytology at 14 days. If any of these reveal leptomeningeal dissemination, Aiden's prognosis changes dramatically, and his treatment intensifies."

---

### Step 4: Drug Discovery — Novel SHH Inhibitors

> **TIMING: 5 minutes**
> **PRESENTER NOTE:** This is the climax of the entire demo series. The standard drug is
> contraindicated. The platform generates novel alternatives. This is generative chemistry
> for a real clinical need.

**Actions:**

1. Switch to the **Drug Discovery Portal** tab (http://localhost:8510)
2. Walk through the interface, explaining each component:

**Target Selection:**
- Target: **SMO receptor** (Smoothened — the target of vismodegib in the SHH pathway)
- PDB structure: **5L7D** (crystal structure of human Smoothened)

**MolMIM Generation:**
- Show the BioNeMo MolMIM interface
- Explain: "MolMIM is generating novel small molecules designed to bind the SMO receptor — the same target as vismodegib — but with different pharmacological properties"
- Show the molecular candidates being generated

**DiffDock Scoring:**
- Show DiffDock predicting binding poses and affinity scores
- Explain: "DiffDock scores how well each candidate binds to the SMO receptor. Higher affinity means lower doses needed, which means less toxicity"

**Pediatric Safety Filters:**
- This is the critical differentiator. Show the safety filter results:

| Safety Filter | What It Checks | Why It Matters |
|---------------|----------------|----------------|
| **BBB Penetration** | Can the drug cross the blood-brain barrier? | Medulloblastoma is in the brain — the drug must reach it |
| **Growth Plate Safety** | Does the drug affect chondrocyte proliferation? | The reason vismodegib is contraindicated |
| **hERG Channel** | Risk of cardiac arrhythmia? | QT prolongation is dangerous in children |
| **Hepatotoxicity** | Liver toxicity risk? | Children's livers are still developing |

3. Show the top 3 candidates ranked by composite score (binding affinity + safety profile)

**What you will see:** Drug discovery results showing:

- **3-5 novel molecular candidates** generated by MolMIM
- **DiffDock binding scores** for each candidate against SMO receptor
- **RDKit property calculations:** Molecular weight, LogP, TPSA, Lipinski rule compliance
- **Pediatric safety scores** for each filter category
- **Top candidates** ranked by composite score (efficacy + safety)
- **Structural visualization** of top candidate bound to SMO receptor

> **PRESENTER NOTE:** Emphasize that these are novel molecules — they do not exist in any
> drug database. They were generated by AI specifically for this child's clinical need.
> This is not drug repurposing. This is drug creation.

**Say:**

> "Since vismodegib would harm a growing child, we asked the Drug Discovery engine to create something new. MolMIM generated novel SMO antagonists — molecules that have never existed before. DiffDock scored their binding affinity against the same pocket vismodegib uses. And then the pediatric safety filter screened every candidate for growth plate toxicity, brain penetration, cardiac safety, and liver toxicity. These are not drugs pulled from a database. These are molecules invented by AI for one specific clinical need: a 10-year-old boy whose cancer needs a drug that does not exist yet."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Drug Discovery Portal showing the top 3 generated molecular candidates with their DiffDock binding scores against SMO receptor, RDKit property calculations (MW, LogP, TPSA), pediatric safety filter results (BBB penetration, growth plate safety, hERG, hepatotoxicity), and the structural visualization of the top candidate bound to SMO]

> ⚠️ **IF THIS FAILS:** If MolMIM generation is slow or the Drug Discovery portal is not responding, show the interface layout and explain the workflow conceptually: "The portal takes the SMO crystal structure, generates novel molecules conditioned on the binding pocket, scores them with DiffDock, and filters for pediatric safety." If pre-loaded results are available, use those. If the portal is completely down, restart with `docker restart drug-discovery-portal` and allow 60 seconds. This step is the climax of the demo — do not skip it. Describe the process in detail even if the live demo is unavailable.

> 🎤 **SAY:** "This is why the Drug Discovery engine exists. Vismodegib would work on Aiden's tumor but destroy his growth. So we asked MolMIM — NVIDIA's generative chemistry model — to invent something new. The model takes the crystal structure of the SMO receptor, the same target vismodegib binds, and generates novel small molecules designed to fit that same binding pocket. DiffDock then scores how well each candidate binds. And here is what makes this pediatric-specific: every candidate is screened through four safety filters. Can it cross the blood-brain barrier to reach a brain tumor? Does it affect growth plate chondrocytes? Does it cause cardiac arrhythmia? Is it toxic to a developing liver? These molecules have never existed before. They were created by AI, right now, for one specific child who needs a drug that the pharmaceutical industry has not made yet."

---

### Step 5: Clinical Trial — Brain Tumor Trials

> **TIMING: 2 minutes**

**Actions:**

1. Switch to the **Clinical Trial Agent** tab (http://localhost:8128)
2. Navigate to **"Patient-Trial Matching"** tab
3. Enter query:

```
Pediatric medulloblastoma SHH subtype clinical trials for a 10-year-old male.
PTCH1 mutant, TP53 wild-type. Include SMO inhibitor trials with pediatric safety
data and reduced-intensity radiation trials.
```

4. Submit the query

**What you will see:** Matched clinical trials:

- **SJMB12 (St. Jude):** Risk-adapted therapy for medulloblastoma by molecular subtype
- **Pediatric MATCH:** SHH arm with targeted agents
- **ACNS1422 or similar COG trials:** Reduced craniospinal radiation for WNT/SHH favorable subtypes
- **Novel SHH inhibitor trials:** Newer-generation SMO inhibitors with potentially better safety profiles
- **Proton therapy trials:** Neurocognitive outcome comparisons

**Say:**

> "While the Drug Discovery engine works on novel candidates, the Clinical Trial agent matches Aiden to existing trials. SJMB12 is a risk-adapted protocol that may allow reduced-intensity radiation for favorable SHH medulloblastoma. And there are newer-generation SMO inhibitors in trial that may have better pediatric safety profiles than vismodegib."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Patient-Trial Matching results showing SJMB12 (St. Jude), Pediatric MATCH SHH arm, COG reduced craniospinal radiation trials, novel SHH inhibitor trials, and proton therapy neurocognitive outcome trials with eligibility scores]

> ⚠️ **IF THIS FAILS:** If the Clinical Trial agent returns limited results, try broadening to: "Pediatric medulloblastoma clinical trials, any molecular subtype, age 10." If the agent is unresponsive, say: "The Clinical Trial agent matches Aiden to St. Jude's SJMB12 risk-adapted protocol, Pediatric MATCH for SHH-targeted agents, and reduced-intensity radiation trials — all designed to cure medulloblastoma while preserving as much of a child's development as possible" and proceed to the summary.

> 🎤 **SAY:** "While the Drug Discovery engine works on molecules that do not exist yet, the Clinical Trial agent matches Aiden to trials that do exist right now. SJMB12 out of St. Jude is a risk-adapted protocol — because Aiden's SHH medulloblastoma is TP53 wild-type, he may qualify for reduced-intensity craniospinal radiation, which means less neurocognitive damage. Pediatric MATCH has an SHH arm testing newer targeted agents. And there are proton therapy trials specifically comparing neurocognitive outcomes. The platform gives Aiden both paths: existing trials today, and novel drug candidates for tomorrow."

---

### Step 6: Demo 6 Summary

> **TIMING: 1 minute**
> **PRESENTER NOTE:** This is the emotional peak. Let it land.

**Say:**

> "Aiden's medulloblastoma had a druggable target — but the drug would have stolen his growth. The platform identified the contraindication that could have been missed. It generated novel alternatives that do not exist in any pharmacy on Earth. It matched him to clinical trials. And it planned for the neurocognitive damage that treatment will cause, even after cure. Five agents. One Drug Discovery engine. One child who deserves a future with both survival and growth."

> 📸 **SCREENSHOT PLACEHOLDER:** [Capture the Landing Page with all agent tiles visible, suitable for a closing shot that represents the full five-agent plus Drug Discovery workflow Aiden went through]

> ⚠️ **IF THIS FAILS:** This is a verbal summary step with no live agent interaction. If any earlier steps in Demo 6 had issues, use this moment to verbally recap the full arc — contraindication caught, novel drugs generated, trials matched, neurocognitive plan established — to ensure the audience leaves with the complete story regardless of any technical hiccups.

> 🎤 **SAY:** "Let me tell you why Aiden's demo matters more than any of the others. Every other demo showed the platform doing something faster or more comprehensively than a human team. Aiden's demo showed the platform doing something that could not have been done at all without it. The standard drug for his tumor pathway would have permanently stunted his growth. The platform caught that. Then it generated novel molecular candidates — drugs that have never existed — designed specifically for a child whose bones are still growing. It screened them for brain penetration, growth plate safety, cardiac toxicity, and liver safety. And it matched him to clinical trials while the novel compounds begin the long road from computation to clinic. This is not just decision support. This is drug creation for an unmet need. A ten-year-old boy who needed a medicine that did not exist, and a platform that started making one."

---

## Post-Demo: Closing Statement

> **TIMING: 2 minutes**
> **PRESENTER NOTE:** Return to the Landing Page. Let the full platform be visible on screen.
> Speak slowly. This is the moment the audience will remember.

**Actions:**

1. Switch to the **Landing Page** (http://localhost:8080)
2. Let the full agent grid be visible
3. Pause for 3 seconds before speaking

**Say:**

> "What you just saw is a single, open-source platform — running on one NVIDIA DGX Spark — that took five children with cancer through complete precision medicine workflows.
>
> From DNA to drug candidates. From diagnosis to clinical trial matching. From risk assessment to toxicity prevention.
>
> Every engine. Every agent. Working together.
>
> Evelyn's ALL — classified, risk-stratified, and matched to trials.
> Marcus's neuroblastoma — every toxicity anticipated before the first dose.
> Aurora's retinoblastoma — germline classified, family protected, brain surveilled.
> Ethan's relapsed leukemia — CAR-T eligibility confirmed, complications mapped, backup strategy ready.
> Aiden's medulloblastoma — a contraindicated drug caught, novel alternatives generated, growth protected.
>
> And it is Apache 2.0. Free. For every hospital on Earth.
>
> That is the HCLS AI Factory."

---

## Appendix A: Agent Port Reference

| Agent | Port | URL | Tabs |
|-------|------|-----|------|
| Landing Page | 8080 | http://localhost:8080 | Pipeline cards, agent grid, health |
| Genomics Portal | 5000 | http://localhost:5000 | Pipeline status |
| RAG Chat | 8501 | http://localhost:8501 | Chat interface |
| Drug Discovery | 8510 | http://localhost:8510 | Target, generation, scoring |
| PGx | 8507 | http://localhost:8507 | Drug interaction, phenotype |
| CAR-T | 8521 | http://localhost:8521 | Chat, collection filter, Deep Research |
| Imaging | 8525 | http://localhost:8525 | 9 tabs |
| Oncology | 8526 | http://localhost:8526 | 5 tabs |
| Biomarker | 8528 | http://localhost:8528 | 8 tabs |
| Neurology | 8529 | http://localhost:8529 | 5 tabs |
| Autoimmune | 8531 | http://localhost:8531 | 10 tabs |
| Single-Cell | 8130 | http://localhost:8130 | 5 tabs |
| Clinical Trial | 8128 | http://localhost:8128 | 5 tabs |
| Cardiology | 8536 | http://localhost:8536 | 10 tabs |
| Rare Disease | 8544 | http://localhost:8544 | 5 tabs |

---

## Appendix B: Troubleshooting

### Agent Not Loading (Streamlit Spinner)

**Symptom:** Browser tab shows infinite Streamlit loading spinner.

**Fix:**
```bash
# Check if the container is running
docker ps | grep <agent-name>

# Restart the specific agent
docker compose -f docker-compose.dgx-spark.yml restart <agent-service-name>

# Wait 10-15 seconds, then refresh the browser tab
```

### Landing Page Shows Red Status

**Symptom:** One or more agents show red health indicator on Landing Page.

**Fix:**
```bash
# Run full health check
bash health-monitor.sh --check-all

# The health monitor will attempt auto-restart of failed services
# If an agent remains red after 60 seconds, manually restart it
docker compose -f docker-compose.dgx-spark.yml restart <service>
```

### RAG Chat Returns Empty or Generic Responses

**Symptom:** RAG Chat does not retrieve variant evidence; responses seem generic.

**Fix:**
```bash
# Verify Milvus is running and has data
docker ps | grep milvus

# Check variant count — should show 3.56M vectors
# The RAG Chat sidebar should display variant stats
# If stats show 0 vectors, the embedding database needs to be reloaded
```

### Slow Response Times

**Symptom:** Agent responses take more than 30 seconds.

**Possible causes:**
- GPU memory contention if multiple agents are running inference simultaneously
- Milvus vector search under heavy load
- Claude API rate limiting

**Fix:**
- Space out queries — don't submit to 3 agents simultaneously during live demo
- Pre-warm agents by submitting a test query to each before the demo starts

### Drug Discovery Portal Not Responding

**Symptom:** Port 8510 not accessible.

**Fix:**
```bash
docker compose -f docker-compose.dgx-spark.yml restart drug-discovery
# Drug Discovery depends on BioNeMo NIMs — verify NIM containers are running
docker ps | grep nim
```

---

## Appendix C: Presenter Quick Reference Card

Print this page and keep it at your podium.

### Demo Flow (Total: 75-90 min)

| # | Demo | Patient | Duration | Key Agents |
|---|------|---------|----------|------------|
| 1 | End-to-End Pipeline | Evelyn, 8, B-ALL | 15-20 min | Landing, Genomics, RAG, Oncology, Drug Discovery |
| 2 | Multi-Agent Tumor Board | Evelyn, Day 29, MRD+ | 12-15 min | Oncology, Biomarker, CAR-T, Single-Cell, Clinical Trial |
| 3 | Cardiotoxicity Prevention | Marcus, 6, Neuroblastoma | 15 min | PGx, Biomarker, Cardiology, Neurology, Autoimmune, Oncology |
| 4 | Rare Disease + Predisposition | Aurora, 4, Retinoblastoma | 12 min | Rare Disease, Oncology, Imaging, Clinical Trial |
| 5 | CAR-T Decision | Ethan, 12, Relapsed ALL | 15 min | Single-Cell, CAR-T, Cardiology, Autoimmune, PGx, Clinical Trial |
| 6 | Medulloblastoma + Drug Discovery | Aiden, 10, SHH Medullo | 15-20 min | Oncology, Neurology, Imaging, Drug Discovery, Clinical Trial |

### Key Dramatic Moments

1. **Demo 1, Step 4:** IKZF1 deletion upgrades ETV6-RUNX1 from standard to high risk
2. **Demo 3, Step 1:** TPMT intermediate metabolizer — full dose would be fatal
3. **Demo 4, Step 2:** RB1 pathogenic classification triggers family-wide testing
4. **Demo 5, Step 2:** CAR-T CRS risk profile — 25% severe
5. **Demo 6, Step 1:** Vismodegib CRITICAL contraindication — growth plate closure
6. **Demo 6, Step 4:** Novel molecules generated for a need no existing drug can fill

### Emergency Shortcuts

If an agent fails during the live demo:
- Acknowledge it: "This is a live system — let me switch to another agent while that restarts"
- Move to the next step and return to the failed agent later
- Never apologize excessively — the audience understands live demos

### Timing Signals

- If running long: skip Genomics Portal (Demo 1 Step 2) and drug discovery walkthrough details
- If running short: add Clinical Scales demo in Neurology (10 interactive scales)
- If audience is technical: show Docker health commands and Milvus vector stats
- If audience is clinical: spend more time on clinical rationale, less on infrastructure

---

*Document generated for the HCLS AI Factory Pediatric Oncology Demo Series.*
*Platform version: DGX Spark deployment | Agents: 11 active | Engines: 3*
*License: Apache 2.0*
