---
title: Therapeutic Discovery Engine
description: AI-driven generative drug design, molecular docking, and drug-likeness scoring
---

# Therapeutic Discovery Engine

**Stage 3 of the HCLS AI Factory Pipeline**

The Therapeutic Discovery Engine transforms validated protein targets from the [Precision Intelligence Network](precision-intelligence.md) into ranked drug candidates using NVIDIA BioNeMo generative chemistry, DiffDock molecular docking, and RDKit drug-likeness scoring. It closes the loop from patient DNA to novel therapeutic candidates.

## What It Does

Protein target → BioNeMo MolMIM generative design → DiffDock molecular docking → RDKit ADMET scoring → Ranked drug candidates with binding predictions

## Pipeline Flow

```
Validated Target (from Precision Intelligence Network)
    |
    v
[Target Preparation]
Structure retrieval (RCSB PDB) + binding site identification
    |
    v
[BioNeMo MolMIM Generative Design]
Generate 10-100 novel molecular candidates per target
    |
    v
[DiffDock Molecular Docking]
Predict binding poses and affinities for each candidate
    |
    v
[RDKit Drug-Likeness Scoring]
Lipinski Rule of 5, QED, SA Score, ADMET properties
    |
    v
[Composite Ranking]
Binding affinity + drug-likeness + synthetic accessibility
    |
    v
Ranked Drug Candidates (PDF report with structures)
```

## Key Numbers

| Metric | Value |
|--------|-------|
| Molecule generation | 10-100 candidates in seconds |
| Docking predictions | Per-candidate binding pose + affinity |
| Drug-likeness filters | Lipinski Ro5, QED, SA Score |
| Structure sources | RCSB PDB (Cryo-EM, X-ray) |
| Output formats | PDF, JSON, Markdown |
| End-to-end time | Minutes per target |

## How It Connects

The Therapeutic Discovery Engine receives validated targets from the Precision Intelligence Network:

- **Precision Oncology Agent** identifies actionable mutations (BRAF V600E, EGFR T790M) with therapy gaps
- **CAR-T Intelligence Agent** flags targets for novel CAR construct design
- **Rare Disease Agent** identifies protein targets for rare disease therapeutics
- **Cross-modal triggers** automatically route high-confidence targets from Stage 2 to Stage 3

## Technology

- **NVIDIA BioNeMo MolMIM** -- Generative molecular design using masked inverse modeling
- **NVIDIA BioNeMo DiffDock** -- Diffusion-based molecular docking for binding prediction
- **RDKit** -- Cheminformatics library for drug-likeness scoring and molecular property calculation
- **RCSB PDB** -- Protein structure database (Cryo-EM, X-ray crystallography)
- **Streamlit UI** -- Interactive drug discovery interface (port 8505)

## Drug-Likeness Scoring

Each generated candidate is evaluated against multiple pharmacological criteria:

| Criterion | Description | Threshold |
|-----------|-------------|-----------|
| Lipinski Rule of 5 | Oral bioavailability prediction | MW <500, LogP <5, HBD <5, HBA <10 |
| QED (Quantitative Estimate of Drug-likeness) | Composite desirability score | 0-1 (higher is better) |
| SA Score (Synthetic Accessibility) | Ease of chemical synthesis | 1-10 (lower is easier) |
| Binding Affinity | Predicted target binding strength | kcal/mol (more negative is stronger) |

## Demo: VCP and Frontotemporal Dementia

The platform includes a complete demonstration targeting VCP (Valosin-containing protein) for Frontotemporal Dementia:

1. **Target identification** -- 13 VCP variants including rs188935092 identified in Stage 1
2. **Evidence synthesis** -- Variants connected to FTD through knowledge graph in Stage 2
3. **Structure retrieval** -- Cryo-EM structures (8OOI, 9DIL, 7K56, 5FTK) from RCSB PDB
4. **Molecule generation** -- Novel VCP inhibitor candidates via BioNeMo MolMIM
5. **Binding prediction** -- Candidates docked to ATP-binding pocket via DiffDock
6. **Ranking** -- Candidates scored by binding affinity, QED, and synthetic accessibility
7. **Report** -- PDF with ranked candidates, structures, and pharmacological properties

## Getting Started

- [Deployment Guide](../HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md) -- Full installation and configuration
- [Architecture](../architecture.md) -- System architecture overview
- [Drug Discovery Pipeline](../drug-discovery-pipeline/README.md) -- Stage 3 pipeline details
- [Performance Benchmarks](../PERFORMANCE.md) -- DGX Spark benchmark data

---

!!! warning "Clinical Decision Support Disclaimer"
    The Therapeutic Discovery Engine is a research tool for computational drug design. Generated molecules are computational candidates only and have not undergone preclinical or clinical validation. All results require expert review and experimental validation. Apache 2.0 License.
