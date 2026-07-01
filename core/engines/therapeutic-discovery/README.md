# Drug Discovery Pipeline

[![NVIDIA BioNeMo](https://img.shields.io/badge/NVIDIA-BioNeMo-76B900?style=flat&logo=nvidia)](https://www.nvidia.com/en-us/clara/bionemo/)
[![DGX Spark](https://img.shields.io/badge/NVIDIA-DGX%20Spark-76B900?style=flat&logo=nvidia)](https://www.nvidia.com/en-us/data-center/dgx-spark/)
[![Streamlit](https://img.shields.io/badge/Streamlit-UI-FF4B4B?style=flat&logo=streamlit)](https://streamlit.io/)
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](../LICENSE)

**Stage 3 of the Precision Medicine to Drug Discovery AI Factory**

> Structure-based drug design pipeline that transforms validated therapeutic targets into novel drug candidates using NVIDIA BioNeMo NIM microservices, Cryo-EM structural evidence, and AI-powered molecule generation.

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                    PRECISION MEDICINE TO DRUG DISCOVERY AI FACTORY                   │
├──────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                      │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────────────┐   │
│  │  GENOMICS   │    │  RAG/CHAT   │    │   CRYO-EM   │    │ MOLECULE GENERATION │   │
│  │  PIPELINE   │───▶│  PIPELINE   │───▶│  EVIDENCE   │───▶│     (BioNeMo)       │   │
│  │             │    │             │    │             │    │    (This Repo)      │   │
│  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────────────┘   │
│    FASTQ→VCF         VCF→Target        Target→Structure    Structure→Molecules      │
│    Parabricks        Milvus+Claude     PDB/EMDB            MolMIM+DiffDock          │
│                                                                                      │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Table of Contents

- [Overview](#overview)
- [From Genomics to Drug Candidates](#from-genomics-to-drug-candidates)
- [Key Features](#key-features)
- [Architecture](#architecture)
- [The VCP/FTD Demo Workflow](#the-vcpftd-demo-workflow)
- [Cryo-EM Structure Evidence](#cryo-em-structure-evidence)
- [Molecule Generation with BioNeMo](#molecule-generation-with-bionemo)
- [Drug-Likeness Scoring](#drug-likeness-scoring)
- [PDF Report Generation](#pdf-report-generation)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Directory Structure](#directory-structure)
- [Services](#services)
- [Monitoring Stack](#monitoring-stack)
- [Performance Benchmarks](#performance-benchmarks)
- [Troubleshooting](#troubleshooting)
- [Related Pipelines](#related-pipelines)
- [References](#references)

---

## Overview

This pipeline is the **final stage** of the Precision Medicine to Drug Discovery AI Factory. It takes validated therapeutic targets from Stage 2 (RAG/Chat Pipeline) and transforms them into actionable drug candidates through:

1. **Cryo-EM Structure Retrieval**: Fetch high-resolution protein structures from RCSB PDB
2. **Known Inhibitor Analysis**: Analyze existing drugs as generation seeds
3. **AI-Powered Molecule Generation**: Generate novel candidates using NVIDIA BioNeMo MolMIM
4. **Molecular Docking**: Predict binding poses with DiffDock
5. **Drug-Likeness Scoring**: Rank candidates by Lipinski, QED, and ADMET properties
6. **Professional PDF Reports**: Generate executive-ready reports with visualizations

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STAGE 3: DRUG DISCOVERY PIPELINE                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Target Hypothesis (from RAG/Chat)                                         │
│   "VCP is a druggable target for Frontotemporal Dementia"                  │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │              PHASE 5: STRUCTURAL EVIDENCE                     │        │
│   │                                                               │        │
│   │   ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐ │        │
│   │   │  8OOI    │   │  9DIL    │   │  7K56    │   │  5FTK    │ │        │
│   │   │WT Hexamer│   │ Mutant   │   │ Complex  │   │+CB-5083  │ │        │
│   │   │ 2.9 Å    │   │ 3.2 Å    │   │ 2.5 Å    │   │ 2.3 Å    │ │        │
│   │   └──────────┘   └──────────┘   └──────────┘   └──────────┘ │        │
│   │                                                               │        │
│   │   Binding Site Analysis: D2 ATPase domain, ATP-competitive   │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │              PHASE 6: MOLECULE GENERATION                     │        │
│   │                                                               │        │
│   │   Seed: CB-5083 (Phase I VCP inhibitor)                      │        │
│   │        │                                                      │        │
│   │        ▼                                                      │        │
│   │   ┌─────────────┐   ┌─────────────┐   ┌─────────────┐        │        │
│   │   │   MolMIM    │   │  DiffDock   │   │   RDKit     │        │        │
│   │   │ Generation  │──▶│  Docking    │──▶│  Scoring    │        │        │
│   │   │  (BioNeMo)  │   │  (BioNeMo)  │   │  (QED/Lip)  │        │        │
│   │   └─────────────┘   └─────────────┘   └─────────────┘        │        │
│   │        │                  │                  │                │        │
│   │        ▼                  ▼                  ▼                │        │
│   │   100 Analogues    Binding Poses     Ranked Candidates       │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   Drug Candidate Report (PDF) → Ready for medicinal chemistry              │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## From Genomics to Drug Candidates

### The Complete Journey

What traditionally takes pharmaceutical companies **months to years** can now be explored in **hours**:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     END-TO-END DRUG DISCOVERY TIMELINE                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   TRADITIONAL                          AI FACTORY                           │
│   ────────────                         ──────────                           │
│                                                                             │
│   Sequencing: 2-4 weeks                Parabricks: 120-240 min             │
│   Variant Analysis: 2-4 weeks          RAG/Chat: Interactive               │
│   Target ID: 3-6 months                Clinker: Instant                     │
│   Structure Analysis: 1-2 months       PDB Fetch: < 1 min                  │
│   Lead Discovery: 6-12 months          MolMIM: 2-5 min                     │
│   Lead Optimization: 1-2 years         DiffDock: 5-10 min                  │
│                                                                             │
│   TOTAL: 2-3 years                     TOTAL: < 5 hours                    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### The Demo Narrative

> "Genomics tells us what changed, Cryo-EM shows us how it changed, and generative AI helps us design what could fix it."

This pipeline demonstrates how NVIDIA's accelerated computing transforms drug discovery from an art into a systematic, data-driven science.

---

## Key Features

### Cryo-EM Structure Integration
- **PDB/EMDB Access**: Automatic retrieval of high-resolution protein structures
- **Structure Gallery**: Visual display of available conformations
- **Binding Site Analysis**: ATP-competitive pocket characterization
- **Resolution Tracking**: Quality metrics for each structure (2.3-3.2 Å)

### NVIDIA BioNeMo NIM Microservices
- **MolMIM**: Generative AI for molecule design (masked modeling)
- **DiffDock**: Diffusion-based molecular docking
- **GPU Acceleration**: 10-100x faster than CPU-based methods
- **REST API**: Easy integration via HTTP endpoints

### Known Inhibitor Analysis
- **CB-5083**: Clinical-stage VCP inhibitor as generation seed
- **Structure-Activity Relationships**: Analyze what makes inhibitors work
- **Multi-Seed Support**: Generate from multiple reference compounds

### Drug-Likeness Scoring
- **Lipinski's Rule of Five**: MW, LogP, HBD, HBA compliance
- **QED Score**: Quantitative estimate of drug-likeness (0-1)
- **Synthetic Accessibility**: Ease of chemical synthesis
- **Tanimoto Similarity**: Distance from seed compounds

### Professional PDF Reports
- **Executive-Ready**: Designed for VP-level presentations
- **Cryo-EM Visuals**: Embedded structure images from PDB
- **Molecule Graphics**: 2D structure renderings from PubChem
- **NVIDIA Branding**: Professional color scheme matching DGX Spark

### Real-Time Monitoring
- **Grafana Dashboards**: GPU utilization, memory, power
- **Prometheus Metrics**: Time-series data collection
- **DCGM Exporter**: NVIDIA-specific GPU metrics

---

## Architecture

### System Architecture

```
┌────────────────────────────────────────────────────────────────────────────────┐
│                          DRUG DISCOVERY PIPELINE                                │
├────────────────────────────────────────────────────────────────────────────────┤
│                                                                                │
│   ┌──────────────────────────────────────────────────────────────────────┐    │
│   │                    STREAMLIT UI (Port 8505)                           │    │
│   │                                                                       │    │
│   │   ┌──────────────┐  ┌──────────────┐  ┌──────────────┐              │    │
│   │   │   Target     │  │  Structure   │  │  Molecule    │              │    │
│   │   │  Hypothesis  │  │   Gallery    │  │  Generation  │              │    │
│   │   └──────────────┘  └──────────────┘  └──────────────┘              │    │
│   │                                                                       │    │
│   │   ┌──────────────┐  ┌──────────────┐  ┌──────────────┐              │    │
│   │   │   Docking    │  │   Scoring    │  │    Report    │              │    │
│   │   │   Results    │  │   & Ranking  │  │  Generation  │              │    │
│   │   └──────────────┘  └──────────────┘  └──────────────┘              │    │
│   │                                                                       │    │
│   └───────────────────────────────┬───────────────────────────────────────┘    │
│                                   │                                            │
│   ┌───────────────────────────────▼───────────────────────────────────────┐    │
│   │                        PIPELINE CORE                                   │    │
│   │                                                                       │    │
│   │   ┌────────────────────────────────────────────────────────────────┐ │    │
│   │   │                    src/pipeline.py                              │ │    │
│   │   │                                                                 │ │    │
│   │   │  Target Import → Structure Fetch → Generation → Docking → Score│ │    │
│   │   │                                                                 │ │    │
│   │   └────────────────────────────────────────────────────────────────┘ │    │
│   │                                                                       │    │
│   └───────────┬───────────────────┬───────────────────┬───────────────────┘    │
│               │                   │                   │                        │
│               ▼                   ▼                   ▼                        │
│   ┌───────────────────┐ ┌───────────────────┐ ┌───────────────────┐           │
│   │    RCSB PDB       │ │   BioNeMo NIMs    │ │      RDKit        │           │
│   │                   │ │                   │ │                   │           │
│   │  Structure Data   │ │  MolMIM (8001)    │ │  Cheminformatics  │           │
│   │  8OOI, 5FTK, etc  │ │  DiffDock (8002)  │ │  QED, Lipinski    │           │
│   │                   │ │                   │ │  SMILES parsing   │           │
│   └───────────────────┘ └───────────────────┘ └───────────────────┘           │
│                                                                                │
│   ┌──────────────────────────────────────────────────────────────────────┐    │
│   │                      MONITORING STACK                                 │    │
│   │                                                                       │    │
│   │   Grafana (3000) ←── Prometheus (9099) ←── DCGM Exporter (9400)     │    │
│   │                                                                       │    │
│   └──────────────────────────────────────────────────────────────────────┘    │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              DATA FLOW                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   RAG/Chat Pipeline Output                                                  │
│   targets_for_phase5.json                                                  │
│        │                                                                    │
│        │  {                                                                 │
│        │    "gene": "VCP",                                                 │
│        │    "variant": "rs188935092",                                      │
│        │    "disease": "Frontotemporal Dementia",                          │
│        │    "druggability": "HIGH"                                         │
│        │  }                                                                 │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 1. TARGET IMPORT                                              │        │
│   │    • Load hypothesis from RAG/Chat                           │        │
│   │    • Validate gene symbol and variant                        │        │
│   │    • Check druggability annotation                           │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 2. STRUCTURE RETRIEVAL                                        │        │
│   │    • Query RCSB PDB for gene (VCP → p97)                     │        │
│   │    • Fetch structure metadata (resolution, method)           │        │
│   │    • Download structure images for visualization             │        │
│   │    • Identify inhibitor-bound structures (5FTK + CB-5083)   │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 3. MOLECULE GENERATION (MolMIM)                               │        │
│   │    • Extract seed SMILES from known inhibitor                │        │
│   │    • Call BioNeMo MolMIM API                                 │        │
│   │    • Generate 100 structural analogues                       │        │
│   │    • Filter by similarity threshold (0.5-0.8)               │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 4. MOLECULAR DOCKING (DiffDock)                               │        │
│   │    • Load protein structure (5FTK)                           │        │
│   │    • Dock each candidate molecule                            │        │
│   │    • Predict binding pose and affinity                       │        │
│   │    • Score: -8 to -12 kcal/mol (good binders)               │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   ┌──────────────────────────────────────────────────────────────┐        │
│   │ 5. SCORING & RANKING                                          │        │
│   │    • Calculate QED score (drug-likeness)                     │        │
│   │    • Check Lipinski Rule of Five                             │        │
│   │    • Compute Tanimoto similarity to seed                     │        │
│   │    • Rank by composite score                                  │        │
│   └──────────────────────────────────────────────────────────────┘        │
│        │                                                                    │
│        ▼                                                                    │
│   Output: Ranked Drug Candidates                                           │
│   VCP_Drug_Candidate_Report.pdf                                            │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## The VCP/FTD Demo Workflow

### Complete Patient-to-Molecule Journey

This demo shows the full AI Factory workflow using a real genomic variant:

#### Step 1: Genomic Discovery (Stage 1)
```
Patient: HG002 (GIAB reference genome)
Variant: rs188935092 (chr9:35065263 G>A)
Gene: VCP (Valosin-Containing Protein)
Impact: Missense mutation
AlphaMissense: 0.89 (Likely Pathogenic)
```

#### Step 2: Target Validation (Stage 2)
```
Query: "What variants are associated with frontotemporal dementia?"
Result: VCP identified as druggable target
Evidence: 13 VCP variants in genome
Diseases: FTD, ALS, IBMPFD (Inclusion Body Myopathy)
Druggability: HIGH (ATP-competitive site validated)
```

#### Step 3: Structure-Based Design (This Pipeline)
```
Structures: 4 VCP Cryo-EM structures from PDB
Seed: CB-5083 (clinical VCP inhibitor)
Generated: 100 novel analogues
Top Candidates: 4 with QED > 0.35, Lipinski compliant
Best Docking: -10.95 kcal/mol
```

### Why VCP?

VCP (also known as **p97**) is a AAA+ ATPase that plays critical roles in:
- **Protein quality control**: Extracts misfolded proteins for degradation
- **Autophagy**: Clears damaged organelles
- **DNA repair**: Removes proteins from chromatin

Mutations in VCP cause:
- **Frontotemporal Dementia (FTD)**: Progressive brain degeneration
- **ALS**: Motor neuron disease
- **IBMPFD**: Muscle, bone, and brain disorder

VCP is an **ideal drug target** because:
- Well-characterized ATP-binding pocket
- Multiple high-resolution structures available
- Clinical-stage inhibitors exist (CB-5083)
- Clear disease mechanism (loss of protein homeostasis)

---

## Cryo-EM Structure Evidence

### Available VCP Structures

| PDB ID | Description | Method | Resolution | Key Feature |
|--------|-------------|--------|------------|-------------|
| **8OOI** | VCP/p97 wild-type hexamer | Cryo-EM | 2.9 Å | Native conformation |
| **9DIL** | VCP with disease mutation | Cryo-EM | 3.2 Å | Mutant structure |
| **7K56** | VCP-cofactor complex | Cryo-EM | 2.5 Å | Binding mechanism |
| **5FTK** | VCP + CB-5083 inhibitor | X-ray | 2.3 Å | Drug binding site |

### Structure Analysis

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        VCP STRUCTURE GALLERY                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐          │
│   │      8OOI       │   │      5FTK       │   │      7K56       │          │
│   │                 │   │                 │   │                 │          │
│   │   [Hexamer]     │   │   [+CB-5083]    │   │   [Complex]     │          │
│   │                 │   │                 │   │                 │          │
│   │   Wild-type     │   │   Inhibitor     │   │   Cofactor      │          │
│   │   2.9 Å         │   │   2.3 Å         │   │   2.5 Å         │          │
│   └─────────────────┘   └─────────────────┘   └─────────────────┘          │
│                                                                             │
│   Binding Site: D2 ATPase Domain                                           │
│   Mode: ATP-competitive                                                     │
│   Key Residues: ALA464, GLY479, ASP320, GLY215                            │
│   Pocket Volume: ~450 Å³                                                   │
│   Druggability Score: 0.92                                                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Image Caching

The pipeline automatically downloads and caches structure images:

```python
from src.cryoem_evidence import CryoEMEvidence

evidence = CryoEMEvidence()
structures = evidence.get_structures("VCP")

# Images cached in data/structures/image_cache/
# - 8ooi_structure.jpeg
# - 5ftk_structure.jpeg
# - 7k56_structure.jpeg
```

---

## Molecule Generation with BioNeMo

### MolMIM: Masked Modeling for Molecules

NVIDIA BioNeMo's MolMIM uses transformer-based masked modeling to generate novel molecules:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           MolMIM GENERATION                                  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Input: CB-5083 SMILES                                                    │
│   CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5  │
│                                                                             │
│   Process:                                                                  │
│   1. Tokenize SMILES into substructures                                    │
│   2. Randomly mask 15-25% of tokens                                        │
│   3. Model predicts replacements                                           │
│   4. Generate diverse completions                                          │
│   5. Filter by validity and novelty                                        │
│                                                                             │
│   Output: 100 Novel Analogues                                              │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐ │
│   │ VCP-AI-001: CC(C)C1=C(C=C(C=C1)NC2=NC3=...  Sim: 0.98  QED: 0.387 │ │
│   │ VCP-AI-002: CC(N)c1ccc(Nc2ncc3c(ccn3C)n2)... Sim: 0.85  QED: 0.365 │ │
│   │ VCP-AI-003: Cc1ccc(NC2=NC3=C(C=N2)N(C=C3)... Sim: 0.72  QED: 0.454 │ │
│   │ VCP-AI-004: CC(C)c1ccc(NC2=NC3=C(C=N2)N(... Sim: 0.75  QED: 0.387 │ │
│   └─────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### DiffDock: Diffusion-Based Docking

DiffDock predicts how generated molecules bind to the target:

```python
from src.molecule_generator import DiffDockClient

client = DiffDockClient(endpoint="http://localhost:8002")
poses = client.dock(
    protein_pdb="data/structures/5ftk.pdb",
    ligand_smiles="CC(C)C1=C(C=C(C=C1)..."
)

# Returns:
# - Binding pose (PDB coordinates)
# - Predicted affinity (kcal/mol)
# - Confidence score
```

### API Endpoints

| Service | Port | Endpoint | Description |
|---------|------|----------|-------------|
| **MolMIM** | 8001 | POST /generate | Generate molecules from seed |
| **DiffDock** | 8002 | POST /dock | Predict binding poses |

### Example API Call

```bash
# Generate molecules with MolMIM
curl -X POST http://localhost:8001/generate \
  -H "Content-Type: application/json" \
  -d '{
    "seed_smiles": "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
    "num_samples": 50,
    "temperature": 0.8,
    "similarity_threshold": 0.7
  }'
```

---

## Drug-Likeness Scoring

### Lipinski's Rule of Five

Predicts oral bioavailability based on molecular properties:

| Rule | Threshold | VCP-AI-001 | Status |
|------|-----------|------------|--------|
| Molecular Weight | ≤ 500 Da | 484.6 Da | ✓ PASS |
| LogP | ≤ 5 | 4.92 | ✓ PASS |
| H-Bond Donors | ≤ 5 | 2 | ✓ PASS |
| H-Bond Acceptors | ≤ 10 | 6 | ✓ PASS |

### QED Score

Quantitative Estimate of Drug-likeness (0-1 scale):

```python
from rdkit.Chem.QED import qed

# Calculate QED
qed_score = qed(mol)

# Interpretation:
# > 0.67: Drug-like
# 0.49-0.67: Moderately drug-like
# < 0.49: Less drug-like
```

### Composite Scoring

Candidates are ranked by a weighted combination:

```python
# Normalize docking score: lower is better, map to [0, 1]
dock_normalized = max(0, min(1, (10 + dock_score) / 20))

# Composite score (weights sum to 1.0)
score = (
    0.3 * generation_score +
    0.4 * dock_normalized +
    0.3 * qed_score
)
```

### Top Candidates

| Rank | ID | Docking (kcal/mol) | QED | MW (Da) | LogP | Score |
|------|----|--------------------|-----|---------|------|-------|
| 1 | VCP-AI-001 | -8.62 | 0.387 | 484.6 | 4.92 | 0.444 |
| 2 | VCP-AI-002 | -8.26 | 0.365 | 485.6 | 3.82 | 0.399 |
| 3 | VCP-AI-003 | -9.86 | 0.454 | 456.5 | 4.10 | 0.364 |
| 4 | VCP-AI-004 | -10.95 | 0.387 | 484.6 | 4.92 | 0.356 |

---

## PDF Report Generation

### Professional Executive Reports

The pipeline generates stunning PDF reports suitable for VP-level presentations:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    VCP DRUG CANDIDATE REPORT                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   PRECISION MEDICINE TO DRUG                                                │
│   DISCOVERY                                                                  │
│   AI Factory Pipeline Report                                                │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐ │
│   │ Target: VCP | Patient: HG002 | Generated: January 14, 2026         │ │
│   └─────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│   ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐                     │
│   │GENOMICS │─▶│RAG/CHAT │─▶│STRUCTURE│─▶│MOLECULES│                     │
│   │Phase 1-3│  │Phase 4  │  │Phase 5  │  │Phase 6  │                     │
│   └─────────┘  └─────────┘  └─────────┘  └─────────┘                     │
│                                                                             │
│   1. GENOMIC VARIANT DETECTION                                             │
│      VCP missense variant rs188935092                                      │
│      AlphaMissense: 0.89 (LIKELY PATHOGENIC)                              │
│                                                                             │
│   2. RAG/CHAT TARGET HYPOTHESIS                                            │
│      VCP/p97 confirmed as high-priority target                            │
│                                                                             │
│   3. STRUCTURAL EVIDENCE                                                    │
│      ┌───────────────┐  ┌───────────────┐                                 │
│      │ [8OOI Image]  │  │ [5FTK Image]  │                                 │
│      │ 2.9 Å Cryo-EM │  │ 2.3 Å X-ray   │                                 │
│      └───────────────┘  └───────────────┘                                 │
│                                                                             │
│   4. GENERATED DRUG CANDIDATES                                             │
│      ┌─────────────────────────────────────────────────────────────────┐ │
│      │ Rank │ ID        │ Docking  │ QED   │ Score │                   │ │
│      │──────│───────────│──────────│───────│───────│                   │ │
│      │ #1   │ VCP-AI-001│ -8.62    │ 0.387 │ 0.444 │                   │ │
│      │ #2   │ VCP-AI-002│ -8.26    │ 0.365 │ 0.399 │                   │ │
│      └─────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│   5. EXECUTIVE SUMMARY                                                      │
│      • Pathogenic variant identified                                       │
│      • Validated drug target                                               │
│      • 4 novel candidates generated                                        │
│                                                                             │
│   ─────────────────────────────────────────────────────────────────────── │
│   HCLS AI Factory | Powered by NVIDIA Accelerated Computing               │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Report Features

- **Cryo-EM Structure Images**: Downloaded from RCSB PDB
- **Molecule Graphics**: 2D structures from PubChem
- **KeepTogether Tables**: No awkward page breaks
- **NVIDIA Color Scheme**: #76B900 green accents
- **Professional Typography**: Clean, readable fonts

### Generate Report

```python
from generate_vcp_report_enhanced import VCPReportGeneratorEnhanced

generator = VCPReportGeneratorEnhanced(
    output_path="outputs/VCP_Drug_Candidate_Report.pdf"
)
generator.generate()
```

---

## Quick Start

### Prerequisites

- **Python 3.10+**
- **NVIDIA GPU** with 16GB+ VRAM (for NIM services)
- **Docker** (for BioNeMo NIMs and monitoring)
- **NGC API Key** (for pulling NIM containers)
- Completed RAG/Chat Pipeline with target hypothesis

### Installation

```bash
# Clone the repository
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd drug-discovery-pipeline

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Configure environment
cp .env.example .env
nano .env  # Add your NGC_API_KEY
```

### Start the Pipeline

```bash
# Option 1: Streamlit Discovery UI (Recommended)
streamlit run app/discovery_ui.py --server.port 8505

# Option 2: Generate PDF Report
python generate_vcp_report_enhanced.py
```

Access the UI at: **http://localhost:8505**

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd drug-discovery-pipeline
```

### Step 2: Create Virtual Environment

```bash
python -m venv venv
source venv/bin/activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- `streamlit`: Web UI framework
- `rdkit`: Cheminformatics toolkit
- `reportlab`: PDF generation
- `requests`: HTTP client for NIMs
- `pydantic`: Data validation and models

### Step 4: Configure Environment

```bash
cp .env.example .env
```

Edit `.env`:

```bash
# NGC API Key (for BioNeMo NIMs)
NGC_API_KEY=your_key_here

# NIM Endpoints
MOLMIM_URL=http://localhost:8001
DIFFDOCK_URL=http://localhost:8002

# Mock Mode (set to true if NIMs unavailable)
NIM_ALLOW_MOCK_FALLBACK=true
```

### Step 5: Start BioNeMo NIMs (Optional)

```bash
docker-compose up -d molmim diffdock
```

### Step 6: Start Monitoring (Optional)

```bash
cd monitoring
docker-compose up -d
```

---

## Usage

### Streamlit Discovery UI

The main interface for interactive drug discovery:

```bash
streamlit run app/discovery_ui.py --server.port 8505
```

**Features:**
- Step-by-step guided workflow
- Interactive 3D structure viewer
- Real-time molecule generation
- Export to PDF report

### PDF Report Generation

Generate executive-ready reports:

```bash
python generate_vcp_report_enhanced.py
```

Output: `outputs/VCP_Drug_Candidate_Report.pdf`

### Python API

```python
from src.pipeline import DrugDiscoveryPipeline
from src.models import PipelineConfig

config = PipelineConfig(
    target_gene="VCP",
    reference_compound_smiles="CC(C)C1=C(C=C(C=C1)...",
    num_molecules=50,
)
pipeline = DrugDiscoveryPipeline(config=config, use_mock=True)

# Run complete workflow
run = pipeline.run_pipeline()

# Access results
print(run.ranked_candidates)
print(run.status)
```

---

## Configuration

### Environment Variables

```bash
# API Keys
NGC_API_KEY=nvapi-...

# NIM Endpoints
MOLMIM_URL=http://localhost:8001
DIFFDOCK_URL=http://localhost:8002

# Generation Parameters
NUM_MOLECULES=100
SIMILARITY_THRESHOLD=0.7
TEMPERATURE=0.8

# Scoring Weights (must sum to 1.0)
DOCKING_WEIGHT=0.4
GENERATION_WEIGHT=0.3
QED_WEIGHT=0.3

# Mock Mode
NIM_ALLOW_MOCK_FALLBACK=true
```

### Molecule Generation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_MOLECULES` | 100 | Candidates to generate |
| `TEMPERATURE` | 0.8 | Generation diversity (0-1) |
| `SIMILARITY_THRESHOLD` | 0.7 | Min Tanimoto to seed |
| `MAX_MW` | 550 | Maximum molecular weight |
| `MAX_LOGP` | 5.5 | Maximum LogP |

### Docking Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_POSES` | 10 | Poses per molecule |
| `EXHAUSTIVENESS` | 8 | Search thoroughness |
| `BINDING_THRESHOLD` | -6.0 | Min binding affinity |

---

## Directory Structure

```
drug-discovery-pipeline/
├── app/
│   └── discovery_ui.py              # Main Streamlit UI (Port 8505)
├── src/
│   ├── pipeline.py                  # Main pipeline orchestrator (includes scoring)
│   ├── target_import.py             # Import from RAG/Chat
│   ├── cryoem_evidence.py           # Cryo-EM structure handling
│   ├── structure_viewer.py          # 3D visualization
│   ├── molecule_generator.py        # BioNeMo MolMIM integration
│   ├── nim_clients.py               # NIM API clients
│   ├── cli.py                       # Typer CLI interface
│   └── models.py                    # Data models (Pydantic)
├── generate_vcp_report_enhanced.py  # PDF report generator
├── data/
│   ├── structures/                  # Structure metadata
│   │   ├── vcp_structures.json      # VCP PDB entries
│   │   └── image_cache/             # Cached structure images
│   ├── targets/                     # Imported target hypotheses
│   └── molecules/                   # Reference molecules
├── outputs/
│   ├── VCP_Drug_Candidate_Report.pdf  # Generated reports
│   ├── candidates/                    # Ranked candidates
│   └── docking/                       # DiffDock results
├── monitoring/
│   ├── docker-compose.yml           # Prometheus + Grafana
│   ├── grafana/
│   │   └── dashboards/              # GPU monitoring dashboards
│   └── prometheus/
│       └── prometheus.yml           # Scrape configuration
├── docker-compose.yml               # BioNeMo NIM services
├── requirements.txt                 # Python dependencies
├── .env.example                     # Environment template
└── README.md                        # This documentation
```

---

## Services

| Service | Port | Description | URL |
|---------|------|-------------|-----|
| **Discovery UI** | 8505 | Main Streamlit interface | http://localhost:8505 |
| **MolMIM NIM** | 8001 | Molecule generation | http://localhost:8001 |
| **DiffDock NIM** | 8002 | Molecular docking | http://localhost:8002 |
| **Grafana** | 3000 | GPU monitoring dashboards | http://localhost:3000 |
| **Prometheus** | 9099 | Metrics collection | http://localhost:9099 |
| **DCGM Exporter** | 9400 | NVIDIA GPU metrics | http://localhost:9400 |
| **Node Exporter** | 9100 | System metrics | http://localhost:9100 |

---

## Monitoring Stack

### Start Monitoring

```bash
cd monitoring
docker-compose up -d
```

### Grafana Dashboard

Access GPU monitoring:
- **URL**: http://localhost:3000
- **Username**: admin
- **Password**: dgxspark

### Available Metrics

| Metric | Description |
|--------|-------------|
| `DCGM_FI_DEV_GPU_UTIL` | GPU compute utilization |
| `DCGM_FI_DEV_MEM_COPY_UTIL` | Memory bandwidth utilization |
| `DCGM_FI_DEV_FB_USED` | GPU memory used |
| `DCGM_FI_DEV_POWER_USAGE` | Power consumption (W) |
| `DCGM_FI_DEV_SM_CLOCK` | SM clock speed |

---

## Performance Benchmarks

### Expected Timings (DGX Spark GB10)

| Step | Time | Notes |
|------|------|-------|
| Target Import | < 1 sec | Load from JSON |
| Structure Fetch | 5-10 sec | PDB API + image download |
| MolMIM Generation | 2-5 min | 100 molecules |
| DiffDock Docking | 5-10 min | 10 poses each |
| Scoring | < 30 sec | RDKit calculations |
| PDF Report | 10-30 sec | With image embedding |
| **Total** | **8-16 min** | End-to-end |

### GPU Utilization

| Phase | GPU Util | Memory | Power |
|-------|----------|--------|-------|
| MolMIM | 70-85% | 8-12 GB | 200-300W |
| DiffDock | 80-95% | 12-16 GB | 250-350W |
| Scoring | 5-10% | 2 GB | 50-100W |

### Comparison to Traditional Methods

| Method | Time | Hardware |
|--------|------|----------|
| GPU Pipeline (This) | 8-16 min | DGX Spark |
| CPU Pipeline | 2-4 hours | 32-core Xeon |
| Manual Discovery | 6-12 months | Chemistry lab |

---

## Troubleshooting

### BioNeMo NIM Not Responding

```bash
# Check container status
docker ps | grep bionemo

# View logs
docker logs dd-molmim

# Restart NIMs
docker-compose restart molmim diffdock
```

### GPU Out of Memory

```bash
# Reduce batch size
python src/cli.py --batch-size 10

# Check GPU memory
nvidia-smi

# Clear GPU memory
nvidia-smi --gpu-reset
```

### PDF Report Generation Failed

```bash
# Check ReportLab installation
pip install --upgrade reportlab

# Verify image cache
ls -la data/structures/image_cache/

# Clear cache and regenerate
rm -rf data/structures/image_cache/*
python generate_vcp_report_enhanced.py
```

### Structure Images Not Loading

```bash
# Check RCSB PDB connectivity
curl https://cdn.rcsb.org/images/structures/8o/8ooi/8ooi_assembly-1.jpeg

# Clear and refetch images
rm -rf data/structures/image_cache/*
python -c "from src.cryoem_evidence import CryoEMEvidence; e = CryoEMEvidence(); e.get_structures('VCP')"
```

### Monitoring Not Showing Data

```bash
# Check Prometheus targets
curl http://localhost:9099/api/v1/targets

# Verify DCGM exporter
curl http://localhost:9400/metrics | grep DCGM

# Restart monitoring stack
cd monitoring && docker-compose down && docker-compose up -d
```

---

## Related Pipelines

| Stage | Pipeline | Description |
|-------|----------|-------------|
| **1** | [Genomics Pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/genomics-pipeline) | FASTQ → VCF with Parabricks |
| **2** | [RAG/Chat Pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/rag-chat-pipeline) | VCF → Target Hypothesis |
| **3** | **Drug Discovery Pipeline** (This repo) | Target → Molecule Candidates |

### Integration Flow

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│    GENOMICS     │     │    RAG/CHAT     │     │ DRUG DISCOVERY  │
│    PIPELINE     │────▶│    PIPELINE     │────▶│    PIPELINE     │
├─────────────────┤     ├─────────────────┤     ├─────────────────┤
│ Input:          │     │ Input:          │     │ Input:          │
│  FASTQ files    │     │  VCF file       │     │  Target JSON    │
│                 │     │                 │     │                 │
│ Process:        │     │ Process:        │     │ Process:        │
│  Parabricks     │     │  Milvus + Claude│     │  BioNeMo NIMs   │
│  DeepVariant    │     │  Clinker        │     │  RDKit          │
│                 │     │                 │     │                 │
│ Output:         │     │ Output:         │     │ Output:         │
│  HG002.vcf.gz   │     │  VCP target     │     │  Candidates PDF │
│  (11.7M vars)   │     │  hypothesis     │     │  (100 molecules)│
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

---

## References

### NVIDIA Technologies

- [NVIDIA BioNeMo](https://www.nvidia.com/en-us/clara/bionemo/)
- [NVIDIA DGX Spark](https://www.nvidia.com/en-us/data-center/dgx-spark/)
- [NVIDIA Clara](https://www.nvidia.com/en-us/clara/)
- [NGC Catalog](https://catalog.ngc.nvidia.com/)

### Structural Biology

- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [EMDB (Electron Microscopy Data Bank)](https://www.ebi.ac.uk/emdb/)
- [PDB File Format](https://www.wwpdb.org/documentation/file-format)

### Cheminformatics

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [SMILES Notation](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
- [QED Score Paper](https://www.nature.com/articles/nchem.1243)
- [Lipinski's Rule of Five](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five)

### VCP/p97 Biology

- [VCP in Protein Quality Control](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3856598/)
- [VCP-Associated Diseases](https://www.ncbi.nlm.nih.gov/books/NBK1476/)
- [CB-5083 Clinical Trial](https://clinicaltrials.gov/ct2/show/NCT02243917)

---

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](../LICENSE) file for details.

---

## Acknowledgments

- **NVIDIA** for BioNeMo NIMs, DGX Spark, and accelerated computing platform
- **RCSB PDB** for open structural biology data
- **RDKit** for open-source cheminformatics toolkit
- **Google DeepMind** for AlphaMissense pathogenicity predictions
- **Streamlit** for the interactive web framework
- **Anthropic** for Claude AI reasoning capabilities

---

**Status**: Production Ready | Optimized for NVIDIA DGX Spark | FTD → Drug Discovery Demo
