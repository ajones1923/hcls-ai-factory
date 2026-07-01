# Genomics to Drug Discovery Pipeline - Complete Phases Overview

## Pipeline Architecture

```
Phase 1-3: Genomics Pipeline (Data Acquisition & Processing)
           ↓
Phase 4:   RAG Chat (Evidence Search & Query)
           ↓
Phase 5:   Target Selection (Hypothesis & Structure)
           ↓
Phase 6:   Molecule Generation (Drug Candidates)
```

---

## Phase 1: Environment Setup & Prerequisites

**Objective**: Prepare the compute environment for GPU-accelerated genomics processing.

```
┌─────────────────────────────────────────────────────────────┐
│                    Phase 1: Prerequisites                    │
├─────────────────────────────────────────────────────────────┤
│  Script: 00-setup-check.sh                                  │
│                                                             │
│  ├── Docker daemon verification                             │
│  ├── NVIDIA Container Runtime check                         │
│  ├── GPU detection and CUDA validation                      │
│  └── Disk space verification (500GB+)                       │
└─────────────────────────────────────────────────────────────┘
```

### Deliverables

| Component | Output |
|-----------|--------|
| Docker | Container runtime operational |
| NVIDIA Runtime | GPU passthrough enabled |
| GPU | CUDA-capable device detected |
| Storage | 500GB+ available space confirmed |

### Commands
```bash
./run.sh check        # Run prerequisites check
./run.sh status       # View system status
```

---

## Phase 2: Authentication & Container Setup

**Objective**: Authenticate with NVIDIA NGC and pull the Parabricks container.

```
┌─────────────────────────────────────────────────────────────┐
│                   Phase 2: Authentication                    │
├─────────────────────────────────────────────────────────────┤
│  Script: 01-ngc-login.sh                                    │
│                                                             │
│  ├── NGC API key authentication                             │
│  ├── Docker login to nvcr.io                                │
│  └── Pull clara-parabricks:4.6.0-1 (~15GB)                  │
└─────────────────────────────────────────────────────────────┘
```

### Deliverables

| Component | Output |
|-----------|--------|
| NGC Auth | Valid API token configured |
| Docker Registry | Authenticated to nvcr.io |
| Parabricks Image | clara-parabricks:4.6.0-1 pulled locally |

### Commands
```bash
./run.sh login        # NGC authentication + container pull
```

---

## Phase 3: Data Acquisition & Genome Processing

**Objective**: Download reference data, acquire sample data, and run GPU-accelerated variant calling.

```
┌─────────────────────────────────────────────────────────────┐
│              Phase 3: Data & Processing                      │
├─────────────────────────────────────────────────────────────┤
│  Scripts: 02-download-data.sh                               │
│           03-setup-reference.sh                             │
│           04-run-chr20-test.sh                              │
│           05-run-full-genome.sh                             │
│                                                             │
│  ├── Download GIAB HG002 FASTQ (~200GB)                     │
│  ├── Setup GRCh38 reference genome                         │
│  ├── GPU-accelerated alignment (BWA-MEM2)                   │
│  ├── Coordinate sorting & duplicate marking                 │
│  ├── DeepVariant CNN variant calling                        │
│  └── VCF generation with quality scores                     │
└─────────────────────────────────────────────────────────────┘
```

### Processing Pipeline

```
FASTQ (Raw Reads)
    ↓
┌───────────────────────────────────┐
│           fq2bam (GPU)            │
│  ├── BWA-MEM2 alignment           │
│  ├── Coordinate sorting           │
│  └── Duplicate marking            │
└───────────────────────────────────┘
    ↓
BAM (Aligned Reads)
    ↓
┌───────────────────────────────────┐
│        DeepVariant (GPU)          │
│  ├── Pileup image generation      │
│  ├── CNN inference                │
│  └── Genotype classification      │
└───────────────────────────────────┘
    ↓
VCF (Variant Calls)
```

### Deliverables

| Component | Output |
|-----------|--------|
| FASTQ Data | HG002_R1.fastq.gz, HG002_R2.fastq.gz (~200GB) |
| Reference | GRCh38.fa with BWA index files |
| BAM | HG002.genome.bam (~100GB aligned reads) |
| VCF | HG002.genome.vcf.gz (~4M variants) |

### Commands
```bash
./run.sh download     # Download GIAB HG002 data
./run.sh reference    # Setup reference genome
./run.sh test         # Run chr20 test (~20 min)
./run.sh full         # Full genome processing (~2-3 hrs)
```

---

## Phase 4: RAG Chat System

**Objective**: Build a retrieval-augmented generation system for querying genomic evidence.

```
┌─────────────────────────────────────────────────────────────┐
│                    Phase 4: RAG Chat                         │
├─────────────────────────────────────────────────────────────┤
│  Components:                                                │
│                                                             │
│  ├── Milvus vector database (Docker)                        │
│  ├── VCF → Evidence object extraction                       │
│  ├── Sentence-transformer embeddings                        │
│  ├── Streamlit chat interface                               │
│  └── LLM backend (vLLM or API)                              │
└─────────────────────────────────────────────────────────────┘
```

### Architecture

```
User Query: "What EGFR variants are in HG002?"
    ↓
┌───────────────────────────────────┐
│        Streamlit Chat UI          │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│     Embedding (bge-small-en)      │
│     Query → 384-dim vector        │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│        Milvus Vector Search       │
│     Top-K similar evidence        │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│           LLM Response            │
│   Context + Query → Answer        │
└───────────────────────────────────┘
    ↓
Evidence-backed response to user
```

### Evidence Object Structure

```python
{
    "variant_key": "chr7:55249071:G:A",
    "chrom": "chr7",
    "pos": 55249071,
    "ref": "G",
    "alt": "A",
    "qual": 45.2,
    "gene": "EGFR",
    "tags": ["giab", "vcf", "variant"],
    "source": "GIAB WGS VCF (GRCh38)",
    "summary_text": "Variant at chr7:55249071 G>A in EGFR gene..."
}
```

### Deliverables

| Component | Output |
|-----------|--------|
| Milvus | Running vector database (port 19530) |
| Evidence Collection | Indexed VCF variants with embeddings |
| Chat Interface | Streamlit app (port 8501) |
| LLM Integration | Query-response pipeline |

### Project Structure
```
~/rag-chat-pipeline/
├── docker-compose.yml          # Milvus standalone
├── ingest_evidence.py          # VCF → embeddings → Milvus
├── chat_app.py                 # Streamlit interface
├── llm_backend.py              # vLLM or API wrapper
└── requirements.txt
```

---

## Phase 5: Target Selection

**Objective**: Form drug target hypothesis using RAG evidence and integrate structural biology data.

```
┌─────────────────────────────────────────────────────────────┐
│                  Phase 5: Target Selection                   │
├─────────────────────────────────────────────────────────────┤
│  Components:                                                │
│                                                             │
│  ├── Hypothesis formation from RAG queries                  │
│  ├── Cryo-EM structure fetch (PDB:7SYE)                     │
│  └── Binding site identification                            │
└─────────────────────────────────────────────────────────────┘
```

### Workflow

```
RAG Chat Query
"What variants affect EGFR in this sample?"
    ↓
┌───────────────────────────────────┐
│       Evidence Retrieval          │
│  Returns: EGFR variants list      │
│  chr7:55249071, chr7:55259515...  │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│      Hypothesis Formation         │
│  "EGFR is a viable drug target    │
│   based on variant evidence"      │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│     Cryo-EM Structure Fetch       │
│  PDB:7SYE (EGFR kinase domain)    │
│  Resolution: 2.8Å                 │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│    Binding Site Identification    │
│  ATP binding pocket coordinates   │
│  Key residues: K745, E762, M793   │
└───────────────────────────────────┘
    ↓
Target profile ready for molecule generation
```

### Deliverables

| Component | Output |
|-----------|--------|
| RAG Query | "EGFR variants in HG002" → evidence list |
| PDB Fetch | 7SYE structure file (~3MB) |
| Binding Site | Coordinates for molecule docking constraints |

### Key Data

```
Target: EGFR (Epidermal Growth Factor Receptor)
PDB ID: 7SYE
Structure: Cryo-EM, 2.8Å resolution
Binding Pocket: ATP-competitive site
Key Residues:
  - K745 (catalytic lysine)
  - E762 (salt bridge)
  - M793 (gatekeeper)
  - T790 (resistance mutation site)
```

---

## Phase 6: Molecule Generation

**Objective**: Generate drug candidate molecules using NVIDIA BioNeMo and MolMIM.

```
┌─────────────────────────────────────────────────────────────┐
│                Phase 6: Molecule Generation                  │
├─────────────────────────────────────────────────────────────┤
│  Components:                                                │
│                                                             │
│  ├── BioNeMo MolMIM model                              │
│  ├── Constraint-based molecule filtering                    │
│  └── SMILES/SDF export for wet lab                          │
└─────────────────────────────────────────────────────────────┘
```

### Architecture

```
Target Profile (from Phase 5)
    ↓
┌───────────────────────────────────┐
│         MolMIM               │
│  Transformer-based generation     │
│  SMILES string output             │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│      Constraint Filtering         │
│  ├── Molecular weight (< 500 Da)  │
│  ├── LogP (drug-likeness)         │
│  ├── Binding affinity score       │
│  └── Toxicity prediction          │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│        Ranked Candidates          │
│  Top molecules by combined score  │
└───────────────────────────────────┘
    ↓
┌───────────────────────────────────┐
│           Export                  │
│  ├── SMILES strings               │
│  ├── SDF 3D structures            │
│  └── CSV with properties          │
└───────────────────────────────────┘
    ↓
Ready for wet lab synthesis/testing
```

### Deliverables

| Component | Output |
|-----------|--------|
| MolMIM | Candidate molecules (SMILES strings) |
| Filtering | Ranked molecules by binding affinity |
| Export | CSV/SDF for wet lab or further simulation |

### Example Output

```
Rank  SMILES                                      MW      LogP   Affinity
────  ──────────────────────────────────────────  ──────  ─────  ────────
1     Cc1cccc(Nc2ncnc3cc(OC)c(OC)cc23)c1          325.4   3.2    -9.8
2     COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OC        347.8   3.5    -9.4
3     Cn1cnc2cc3c(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)   421.9   4.1    -9.1
...
```

---

## Complete Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    GENOMICS TO DRUG DISCOVERY PIPELINE                   │
└─────────────────────────────────────────────────────────────────────────┘

Phase 1: Prerequisites          Phase 2: Authentication
┌─────────────────────┐         ┌─────────────────────┐
│ Docker + GPU Check  │ ──────► │ NGC Login + Pull    │
└─────────────────────┘         └─────────────────────┘
                                          │
                                          ▼
                                Phase 3: Processing
                        ┌─────────────────────────────────┐
                        │  FASTQ → BAM → VCF              │
                        │  (GPU-accelerated Parabricks)   │
                        └─────────────────────────────────┘
                                          │
                                          ▼
                                  Phase 4: RAG Chat
                        ┌─────────────────────────────────┐
                        │  Milvus + Evidence + Streamlit  │
                        │  Query: "EGFR variants?"        │
                        └─────────────────────────────────┘
                                          │
                                          ▼
                              Phase 5: Target Selection
                        ┌─────────────────────────────────┐
                        │  Hypothesis + Cryo-EM (7SYE)    │
                        │  Binding site identification    │
                        └─────────────────────────────────┘
                                          │
                                          ▼
                            Phase 6: Molecule Generation
                        ┌─────────────────────────────────┐
                        │  BioNeMo MolMIM            │
                        │  SMILES → Wet Lab               │
                        └─────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════
 Input: Raw FASTQ reads (~200GB)
 Output: Ranked drug candidate molecules (SMILES/SDF)
═══════════════════════════════════════════════════════════════════════════
```

---

## Technology Stack

| Phase | Primary Technologies |
|-------|---------------------|
| 1-3 | Docker, NVIDIA Parabricks, CUDA, BWA-MEM2, DeepVariant |
| 4 | Milvus, sentence-transformers, Streamlit, vLLM/OpenAI |
| 5 | PDB/RCSB API, PyMOL/ChimeraX, BioPython |
| 6 | NVIDIA BioNeMo, MolMIM, RDKit |

---

## Hardware Requirements

| Phase | GPU | RAM | Storage |
|-------|-----|-----|---------|
| 1-3 | NVIDIA GPU (16GB+ VRAM) | 64GB+ | 500GB+ NVMe |
| 4 | Optional (CPU embeddings OK) | 32GB+ | 50GB |
| 5 | Optional | 16GB+ | 10GB |
| 6 | NVIDIA GPU (24GB+ VRAM recommended) | 64GB+ | 100GB |

---

## Status Checklist

- [x] Phase 1: Prerequisites - Complete
- [x] Phase 2: Authentication - Complete
- [x] Phase 3: Data & Processing - Complete
- [ ] Phase 4: RAG Chat - Pending
- [ ] Phase 5: Target Selection - Pending
- [ ] Phase 6: Molecule Generation - Pending
