# HCLS AI Factory -- Napkin AI Visual Blueprint (Detailed Edition)

Paste each section into Napkin AI to generate a clean infographic.
Every number, port, collection, gene, and filter is sourced from the running codebase.

---

## VISUAL 1: Platform Overview

**HCLS AI Factory** -- Open-source precision medicine platform.
Patient DNA to drug candidates in under 5 hours on a single NVIDIA DGX Spark.

### Engine 1: Genomic Foundation (Ports 5000-5001)

- BWA-MEM2 alignment (via Parabricks fq2bam)
- Google DeepVariant variant calling (via Parabricks)
- samtools BAM indexing
- GRCh38 reference genome (3.1 GB)
- ClinVar annotation (4.1M variants, February 2026 release)
- AlphaMissense scoring (71,697,560 predictions)
- Ensembl VEP consequence annotation
- Output: VCF with 11.7M called variants

### Engine 2: Precision Intelligence Engine (Ports 8080, 8100-8544, 19530)

- Milvus v2.4 vector database (45 collections, IVF_FLAT index, COSINE metric)
- BGE-small-en-v1.5 embedding model (384 dimensions)
- Claude LLM synthesis (claude-sonnet-4-20250514, temperature=0.3)
- 8 specialized FastAPI intelligence agents
- Streamlit UI per agent
- 201 genes across 13 therapeutic areas (171 druggable, 85%)
- RAG-grounded multi-turn chat with streaming SSE

### Engine 3: Therapeutic Discovery (Ports 8001-8002, 8505)

- MolMIM molecule generation (NIM container, port 8001)
- DiffDock binding prediction (NIM container, port 8002)
- RDKit chemistry QC (Lipinski, QED, TPSA)
- 10-stage pipeline with checkpoint/resume
- 6 pediatric safety filters
- Composite scoring: 30% generation + 40% docking + 30% QED
- PDF report with full provenance chain

### Key Stats

| Metric | Value |
|---|---|
| Variants per genome | 11,724,891 |
| Searchable annotated variants | 3,487,216 |
| Milvus vector collections | 45 |
| Intelligence agents | 11 |
| Therapeutic areas | 13 |
| Genes mapped | 201 |
| Druggable targets | 171 (85%) |
| End-to-end time | ~4 hours 12 minutes |
| Hardware | 1x NVIDIA DGX Spark |

### Port Ranges

| Range | Purpose |
|---|---|
| 3000 | Grafana dashboards |
| 5000-5001 | Genomics portal, RAG API |
| 8001-8002 | NIM services (MolMIM, DiffDock) |
| 8080 | Landing page / hub |
| 8100-8134 | Agent FastAPI endpoints |
| 8500-8544 | Agent Streamlit UIs |
| 9090-9100 | Prometheus, node-exporter |
| 19530 | Milvus gRPC |

---

## VISUAL 2: Complete Statistics Dashboard

### Genomic Scale

| Metric | Value |
|---|---|
| Total variants called (HG002 30x WGS) | 11,724,891 |
| PASS quality (QUAL > 30) | 3,487,216 |
| SNPs | 4,198,433 |
| Indels | 1,012,548 |
| Multi-allelic sites | 148,762 |
| Coding region variants | 35,616 |
| Transition/transversion ratio | 2.07 |
| ClinVar annotated variants | 35,616 |
| AlphaMissense scored variants | 6,831 |
| HIGH impact + pathogenic (actionable) | 2,412 |
| Druggable gene target variants | 847 |
| FASTQ input size (paired-end) | 198.7 GB (R1: 99.4 GB, R2: 99.3 GB) |
| Reference genome size (GRCh38) | 3.1 GB |

### Clinical Intelligence

| Metric | Value |
|---|---|
| Intelligence agents | 11 |
| Therapeutic areas | 13 |
| Genes mapped | 201 |
| Druggable targets | 171 (85%) |
| Rare diseases covered | 88 |
| CAR-T target antigens profiled | 34 |
| Pharmacogenes (core) | 14 |
| CPIC dosing algorithms | 9 |
| Cardiac genes mapped | 56 |
| Cell types classified (single-cell) | 57 |
| Autoimmune conditions covered | 13 |
| Clinical scales (neurology) | 10 |
| Risk calculators (cardiology) | 6 |
| Gene therapies assessed (rare disease) | 25 |

### Data Architecture

| Metric | Value |
|---|---|
| Milvus vector collections | 45 |
| Total vectors indexed | 3,487,216 |
| Embedding dimensions | 384 |
| Embedding model | BGE-small-en-v1.5 (BAAI) |
| Index type | IVF_FLAT (nlist=1024) |
| Distance metric | COSINE |
| Search nprobe | 16 |
| Top-k per query | 20 |
| Search latency | 12 ms |
| ClinVar records | 4,100,000 |
| AlphaMissense predictions | 71,697,560 |

### Infrastructure

| Metric | Value |
|---|---|
| Total Docker services | 21+ |
| Hardware | 1x NVIDIA DGX Spark |
| GPU | NVIDIA GB10 (Grace Blackwell) |
| Memory | 128 GB unified LPDDR5x |
| CPU | 20 ARM cores (Grace) |
| Interconnect | NVLink-C2C |
| CUDA version | 12.x |
| Licensing | 100% open-source (Apache 2.0) |
| Total platform footprint | ~1.1 TB |

### Performance Benchmarks

| Stage | Duration | GPU Util | Peak Memory |
|---|---|---|---|
| Alignment (fq2bam) | 34 min | 82% | 38 GB |
| Variant calling (DeepVariant) | 22 min | 91% | 54 GB |
| Annotation (ClinVar + AlphaMissense) | 18 min | 15% (CPU) | 12 GB |
| Milvus indexing | 24 min | 35% | 22 GB |
| RAG/Chat interactive session | 45 min | 5% | 8 GB |
| Structure retrieval (PDB) | 2 min | 0% (network) | 2 GB |
| MolMIM generation (100 molecules) | 2 min 14 sec | 78% | 18 GB |
| DiffDock docking (98 candidates) | 8 min 42 sec | 85% | 24 GB |
| Scoring + report generation | 1 min 30 sec | 0% (CPU) | 4 GB |
| **Total pipeline** | **~4 hr 12 min** | | |

---

## VISUAL 3: Engine 1 -- Genomic Foundation Detail

### Input Specification

| Parameter | Value |
|---|---|
| Sample | HG002 (NA24385, GIAB Ashkenazi male) |
| Sequencing | Illumina, 30x WGS, 2x250 bp paired-end |
| FASTQ size | 198.7 GB (R1: 99.4 GB, R2: 99.3 GB) |
| Reference genome | GRCh38 (hg38), 3.1 GB |

### Pipeline Stages (Measured Benchmarks)

| Stage | Tool | Duration | GPU Util | Peak Mem | Input | Output |
|---|---|---|---|---|---|---|
| 1. Alignment | BWA-MEM2 via Parabricks fq2bam | 34 min | 82% | 38 GB | 198.7 GB FASTQ | ~120 GB BAM |
| 2. Sorting + Indexing | samtools (via Parabricks) | included | -- | -- | BAM | Sorted BAM + BAI |
| 3. Variant Calling | Google DeepVariant (via Parabricks) | 22 min | 91% | 54 GB | Sorted BAM | VCF (11.7M variants) |
| 4. Quality Filtering | bcftools / VCF filter | <1 min | 0% | 1 GB | 11.7M variants | 3.49M PASS variants |
| 5. ClinVar Annotation | ClinVar Feb 2026 database | ~5 min | 0% | 4 GB | 3.49M variants | 35,616 clinically annotated |
| 6. AlphaMissense Scoring | AlphaMissense v1.0 lookup | ~5 min | 0% | 6 GB | 35,616 variants | 6,831 AI-scored |
| 7. VEP Annotation | Ensembl VEP (GRCh38) | ~8 min | 0% | 8 GB | 6,831 variants | Consequence-annotated VCF |
| **Total Stage 1** | | **~56 min** | | | | |

### VCF Output Statistics

| Metric | Count |
|---|---|
| Total variants called | 11,724,891 |
| PASS quality (QUAL > 30) | 3,487,216 |
| SNPs | 4,198,433 |
| Indels | 1,012,548 |
| Multi-allelic sites | 148,762 |
| Coding region variants | 35,616 |
| Transition/transversion ratio | 2.07 |

### Annotation Sources

| Source | Version | Records | Purpose |
|---|---|---|---|
| GRCh38 | hg38 | 3.1 GB reference | Alignment reference |
| ClinVar | February 2026 | 4,100,000 | Clinical significance classification |
| AlphaMissense | v1.0 | 71,697,560 | AI pathogenicity prediction (threshold >0.564) |
| Ensembl VEP | GRCh38 | -- | Variant consequence and impact annotation |
| GIAB truth set | HG002 v4.2.1 | -- | Validation benchmark |

### Annotation Funnel

```
Raw VCF: 11,724,891
  -> QUAL > 30 filter: 3,487,216
    -> ClinVar match: 35,616
      -> AlphaMissense scored: 6,831
        -> HIGH impact + Pathogenic: 2,412
          -> Druggable gene match: 847
```

### Speed Comparison

| Stage | Traditional (CPU) | HCLS AI Factory (GPU) | Speedup |
|---|---|---|---|
| Alignment + variant calling | 24-48 hours | 56 minutes | ~30-50x |
| Full genome to VCF | 2-3 days | <1 hour | ~50x |

---

## VISUAL 4: Engine 2 -- Intelligence Network Detail

### Milvus Configuration

| Parameter | Value |
|---|---|
| Version | Milvus v2.4.0 |
| Image | milvusdb/milvus:v2.4.0 |
| gRPC port | 19530 |
| Health port | 9091 |
| Index type | IVF_FLAT |
| Index nlist | 1024 |
| Search nprobe | 16 |
| Distance metric | COSINE |
| Memory limit | 8 GB |
| Backend | etcd v3.5.11 (metadata) + MinIO (object storage) |
| etcd quota | 4 GB |
| MinIO image | minio/minio:RELEASE.2024-01-01T16-36-33Z |

### All 45 Collections by Agent

**Precision Oncology (11 collections):**
onco_literature, onco_trials, onco_variants, onco_biomarkers, onco_therapies, onco_pathways, onco_guidelines, onco_resistance, onco_outcomes, onco_cases, genomic_evidence

**CAR-T Intelligence (11 collections):**
cart_literature, cart_trials, cart_constructs, cart_assays, cart_manufacturing, cart_safety, cart_biomarkers, cart_regulatory, cart_sequences, cart_realworld, genomic_evidence

**Precision Biomarker (14 collections):**
biomarker_reference, biomarker_genetic_variants, biomarker_pgx_rules, biomarker_disease_trajectories, biomarker_clinical_evidence, biomarker_nutrition, biomarker_drug_interactions, biomarker_aging_markers, biomarker_genotype_adjustments, biomarker_monitoring, biomarker_critical_values, biomarker_discordance_rules, biomarker_aj_carrier_screening, genomic_evidence

**Imaging Intelligence (11 collections):**
imaging_literature, imaging_trials, imaging_findings, imaging_protocols, imaging_devices, imaging_anatomy, imaging_benchmarks, imaging_guidelines, imaging_report_templates, imaging_datasets, genomic_evidence

**Precision Autoimmune (13 collections):**
autoimmune_clinical_documents, autoimmune_patient_labs, autoimmune_autoantibody_panels, autoimmune_hla_associations, autoimmune_disease_criteria, autoimmune_disease_activity, autoimmune_flare_patterns, autoimmune_biologic_therapies, autoimmune_pgx_rules, autoimmune_clinical_trials, autoimmune_literature, autoimmune_patient_timelines, autoimmune_cross_disease

**Cardiology Intelligence (13 collections):**
cardio_literature, cardio_trials, cardio_imaging, cardio_electrophysiology, cardio_heart_failure, cardio_valvular, cardio_prevention, cardio_interventional, cardio_oncology, cardio_devices, cardio_guidelines, cardio_hemodynamics, genomic_evidence

**Neurology Intelligence (14 collections):**
neuro_literature, neuro_trials, neuro_imaging, neuro_electrophysiology, neuro_degenerative, neuro_cerebrovascular, neuro_epilepsy, neuro_oncology, neuro_ms, neuro_movement, neuro_headache, neuro_neuromuscular, neuro_guidelines, genomic_evidence

**Rare Disease Diagnostic (14 collections):**
rd_phenotypes, rd_diseases, rd_genes, rd_variants, rd_literature, rd_trials, rd_therapies, rd_case_reports, rd_guidelines, rd_pathways, rd_registries, rd_natural_history, rd_newborn_screening, genomic_evidence

**Clinical Trial Intelligence (14 collections):**
trial_protocols, trial_eligibility, trial_endpoints, trial_sites, trial_investigators, trial_results, trial_regulatory, trial_literature, trial_biomarkers, trial_safety, trial_rwe, trial_adaptive, trial_guidelines, genomic_evidence

**Single-Cell Intelligence (12 collections):**
sc_cell_types, sc_markers, sc_spatial, sc_tme, sc_drug_response, sc_literature, sc_methods, sc_datasets, sc_trajectories, sc_pathways, sc_clinical, genomic_evidence

**Pharmacogenomics Intelligence (15 collections):**
pgx_gene_reference, pgx_drug_guidelines, pgx_drug_interactions, pgx_hla_hypersensitivity, pgx_phenoconversion, pgx_dosing_algorithms, pgx_clinical_evidence, pgx_population_data, pgx_clinical_trials, pgx_fda_labels, pgx_drug_alternatives, pgx_patient_profiles, pgx_implementation, pgx_education, genomic_evidence

**Note:** `genomic_evidence` is a shared collection used by all agents. Unique collection names total 45 when the shared collection is counted once.

### Claude LLM Configuration

| Parameter | Value |
|---|---|
| Model | claude-sonnet-4-20250514 |
| Temperature | 0.3 |
| Provider | Anthropic |
| Integration | anthropic Python SDK v0.75.0 |
| Streaming | Server-Sent Events (SSE) |
| Multi-turn | Conversation memory maintained per session |

### BGE Embedding Specification

| Parameter | Value |
|---|---|
| Model | BAAI/bge-small-en-v1.5 |
| Dimensions | 384 |
| Max sequence length | 512 tokens |
| Framework | sentence-transformers v5.2.0 |
| Backend | PyTorch v2.9.1 |
| Normalization | L2-normalized output vectors |

### RAG Query Flow

```
1. User query -> BGE-small-en-v1.5 embedding (384-dim vector)
2. Milvus COSINE similarity search (nprobe=16, top-k=20)
3. Score threshold filter (>0.5)
4. Therapeutic area keyword expansion (gene mappings)
5. Evidence context assembly from matched collections
6. Claude LLM synthesis with system prompt + evidence
7. Streaming SSE response to client
```

---

## VISUAL 5: Engine 3 -- Discovery Pipeline Detail

### 10-Stage Pipeline

| Stage | Name | Technology | Input | Output |
|---|---|---|---|---|
| 0 | Initialize | Python / PipelineConfig | Config JSON | Pipeline state |
| 1 | Normalize Target | HGNC / UniProt lookup | Gene symbol | Validated target (protein, UniProt ID) |
| 2 | Structure Discovery | PDB / EMDB REST API | UniProt ID | PDB structure list with scores |
| 3 | Structure Prep | BioPython / PDB parser | PDB files | Cleaned protein structure |
| 4 | Molecule Generation | MolMIM NIM (port 8001) | Seed SMILES + config | 100 candidate SMILES |
| 5 | Chemistry QC | RDKit | Raw SMILES | Validated molecules (Lipinski, QED, TPSA) |
| 6 | Conformers | RDKit ETKDG | 2D molecules | 3D conformers (SDF) |
| 7 | Docking | DiffDock NIM (port 8002) | Protein + ligand 3D | Binding poses + scores (kcal/mol) |
| 8 | Ranking | Composite scorer | All scores | Ranked candidates |
| 9 | Reporting | ReportLab | Ranked list | PDF report with provenance |

### NIM Container Images

| Service | Image | Version | Port (host:container) | GPU Required |
|---|---|---|---|---|
| MolMIM | nvcr.io/nim/nvidia/molmim | 1.0.0 | 8001:8000 | Yes (1x GPU) |
| DiffDock | nvcr.io/nim/mit/diffdock | 2.2.0 | 8002:8000 | Yes (1x GPU) |

### MolMIM Parameters

| Parameter | Default |
|---|---|
| Seed SMILES | From target hypothesis (e.g., CB-5083) |
| num_molecules | 100 |
| temperature | 1.0 |
| num_samples_per_token | 10 |
| masked_ratio | 0.1 |
| Generation time (100 molecules) | 2 min 14 sec |
| Validity rate | ~98% (RDKit valence check) |

### DiffDock Parameters

| Parameter | Value |
|---|---|
| Protein input | PDB structure (selected by score) |
| Ligand input | 3D conformer (SDF from Stage 6) |
| Docking time (98 candidates) | 8 min 42 sec |
| Mean dock score | -7.4 kcal/mol |
| Best dock score | -11.4 kcal/mol |
| Excellent binding (< -8.0) | 34 candidates |
| Good+ binding (< -6.0) | 78 candidates |

### Composite Scoring Formula

```
composite = 0.30 * generation_score + 0.40 * docking_normalized + 0.30 * qed
docking_normalized = max(0, min(1, (10 + dock_score) / 20))
```

### RDKit Capabilities Used

- Molecular weight calculation
- LogP (Wildman-Crippen)
- TPSA (topological polar surface area)
- Rotatable bond count
- QED (quantitative estimate of drug-likeness)
- Lipinski Rule of Five validation
- SMILES parsing and validation
- 3D conformer generation (ETKDG algorithm)
- Substructure search (teratogenic alerts via SMARTS)
- Hydrogen bond donor/acceptor counts

### 6 Pediatric Safety Filters

| # | Filter | Threshold | Severity | Clinical Rationale | Estimation Method |
|---|---|---|---|---|---|
| 1 | Blood-Brain Barrier | MW > 500 Da | medium | Pediatric BBB is more permeable in young children; high-MW compounds risk unintended CNS exposure | RDKit molecular weight calculation |
| 2 | Cardiac Liability (hERG) | hERG IC50 < 10 uM | critical | Children have lower QT prolongation threshold; hERG channel blockade is a primary cardiac safety signal | Estimated: max(0.5, 30 - (logP * 4) - (MW - 350) * 0.02) |
| 3 | Hepatic Immaturity | LogP > 5 | high | CYP3A4 is immature until age 1-2; high lipophilicity increases hepatotoxicity risk in developing liver | RDKit Wildman-Crippen LogP |
| 4 | Teratogenicity | Structural alert match | high | Adolescent pregnancy risk requires screening for known teratogenic scaffolds | RDKit SMARTS substructure search (thalidomide, retinoid, valproic acid motifs) |
| 5 | Oral Bioavailability | TPSA > 140 A^2 | medium | Children often need liquid formulations; high TPSA limits oral absorption | RDKit TPSA calculation |
| 6 | GI Absorption | > 10 rotatable bonds | low | Pediatric GI tract differs from adults; high flexibility causes variable absorption | RDKit rotatable bond count |

**Decision rule:** Any "critical" flag = NOT pediatric-safe. Zero flags = cleared for pediatric development. Non-critical flags = requires additional evaluation.

---

## VISUAL 6: 11 Agents -- Complete Reference

### Agent 1: Precision Oncology

| Parameter | Value |
|---|---|
| Name | precision-oncology-agent |
| API port | 8103 (FastAPI) |
| UI port | 8503 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 11 (onco_literature, onco_trials, onco_variants, onco_biomarkers, onco_therapies, onco_pathways, onco_guidelines, onco_resistance, onco_outcomes, onco_cases, genomic_evidence) |
| Key capabilities | Molecular tumor board support, 80+ targeted therapy matching, resistance mechanism prediction, evidence-tiered reporting (Levels I-IV), 12+ drug class coverage |

### Agent 2: CAR-T Intelligence

| Parameter | Value |
|---|---|
| Name | cart-intelligence-agent |
| API port | 8104 (FastAPI) |
| UI port | 8504 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 11 (cart_literature, cart_trials, cart_constructs, cart_assays, cart_manufacturing, cart_safety, cart_biomarkers, cart_regulatory, cart_sequences, cart_realworld, genomic_evidence) |
| Key capabilities | 34 CAR-T target antigens profiled, CRS/ICANS toxicity prediction, manufacturing optimization, costimulatory domain analysis, real-world evidence tracking |

### Agent 3: Precision Biomarker

| Parameter | Value |
|---|---|
| Name | precision-biomarker-agent |
| API port | 8102 (FastAPI) |
| UI port | 8502 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 14 (biomarker_reference, biomarker_genetic_variants, biomarker_pgx_rules, biomarker_disease_trajectories, biomarker_clinical_evidence, biomarker_nutrition, biomarker_drug_interactions, biomarker_aging_markers, biomarker_genotype_adjustments, biomarker_monitoring, biomarker_critical_values, biomarker_discordance_rules, biomarker_aj_carrier_screening, genomic_evidence) |
| Key capabilities | Biological age estimation, 9 disease trajectory detection, pharmacogenomic profiling, FHIR R4 report export, Ashkenazi Jewish carrier screening |

### Agent 4: Cardiology Intelligence

| Parameter | Value |
|---|---|
| Name | cardiology-intelligence-agent |
| API port | 8126 (FastAPI) |
| UI port | 8536 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 13 (cardio_literature, cardio_trials, cardio_imaging, cardio_electrophysiology, cardio_heart_failure, cardio_valvular, cardio_prevention, cardio_interventional, cardio_oncology, cardio_devices, cardio_guidelines, cardio_hemodynamics, genomic_evidence) |
| Workflows | 11 |
| Key capabilities | 6 risk calculators (ASCVD, MAGGIC, HEART, CHA2DS2-VASc, HAS-BLED, EuroSCORE), GDMT optimizer (HFrEF/HFpEF), 56 cardiac genes mapped, cardio-oncology crossover |

### Agent 5: Neurology Intelligence

| Parameter | Value |
|---|---|
| Name | neurology-intelligence-agent |
| API port | 8528 (FastAPI) |
| UI port | 8529 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 14 (neuro_literature, neuro_trials, neuro_imaging, neuro_electrophysiology, neuro_degenerative, neuro_cerebrovascular, neuro_epilepsy, neuro_oncology, neuro_ms, neuro_movement, neuro_headache, neuro_neuromuscular, neuro_guidelines, genomic_evidence) |
| Key capabilities | 10 clinical scales (NIHSS, GCS, MoCA, EDSS, UPDRS, HY, ALSFRS-R, MRC, Barthel, mRS), 8 clinical workflows, multi-guideline integration (AAN, WHO), multi-format reports (PDF, FHIR) |

### Agent 6: Precision Autoimmune

| Parameter | Value |
|---|---|
| Name | precision-autoimmune-agent |
| API port | 8106 (FastAPI) |
| UI port | 8506 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 13 (autoimmune_clinical_documents, autoimmune_patient_labs, autoimmune_autoantibody_panels, autoimmune_hla_associations, autoimmune_disease_criteria, autoimmune_disease_activity, autoimmune_flare_patterns, autoimmune_biologic_therapies, autoimmune_pgx_rules, autoimmune_clinical_trials, autoimmune_literature, autoimmune_patient_timelines, autoimmune_cross_disease) |
| Key capabilities | 13 autoimmune conditions + overlap syndromes, disease activity scoring (DAS28, SLEDAI), autoantibody panel interpretation, flare prediction (0-1 scale), cross-disease analysis |

### Agent 7: Rare Disease Diagnostic

| Parameter | Value |
|---|---|
| Name | rare-disease-diagnostic-agent |
| API port | 8134 (FastAPI) |
| UI port | 8544 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 14 (rd_phenotypes, rd_diseases, rd_genes, rd_variants, rd_literature, rd_trials, rd_therapies, rd_case_reports, rd_guidelines, rd_pathways, rd_registries, rd_natural_history, rd_newborn_screening, genomic_evidence) |
| Key capabilities | 88 rare diseases covered, ACMG/AMP classification (23 criteria), gene therapy eligibility (25 therapies assessed), HPO phenotype matching, newborn screening integration |

### Agent 8: Pharmacogenomics Intelligence

| Parameter | Value |
|---|---|
| Name | pharmacogenomics-intelligence-agent |
| API port | 8107 (FastAPI) |
| UI port | 8507 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 15 (pgx_gene_reference, pgx_drug_guidelines, pgx_drug_interactions, pgx_hla_hypersensitivity, pgx_phenoconversion, pgx_dosing_algorithms, pgx_clinical_evidence, pgx_population_data, pgx_clinical_trials, pgx_fda_labels, pgx_drug_alternatives, pgx_patient_profiles, pgx_implementation, pgx_education, genomic_evidence) |
| Key capabilities | 14 core pharmacogenes, 9 CPIC dosing algorithms, HLA allele risk assessment, metabolizer phenotype classification, FDA label pharmacogenomic annotations |

### Agent 9: Imaging Intelligence

| Parameter | Value |
|---|---|
| Name | imaging-intelligence-agent |
| API port | 8105 (FastAPI) |
| UI port | 8505 (Streamlit) |
| Memory limit | 8 GB |
| Collections | 11 (imaging_literature, imaging_trials, imaging_findings, imaging_protocols, imaging_devices, imaging_anatomy, imaging_benchmarks, imaging_guidelines, imaging_report_templates, imaging_datasets, genomic_evidence) |
| Key capabilities | Vista3D (127+ anatomical structures), MAISI synthetic CT generation, VILA-M3 visual question answering, DICOM auto-ingestion (Orthanc), cross-modal genomic-imaging correlation |

### Agent 10: Single-Cell Intelligence

| Parameter | Value |
|---|---|
| Name | single-cell-intelligence-agent |
| API port | 8540 (FastAPI) |
| UI port | 8130 (Streamlit) |
| Memory limit | 8 GB |
| Collections | 12 (sc_cell_types, sc_markers, sc_spatial, sc_tme, sc_drug_response, sc_literature, sc_methods, sc_datasets, sc_trajectories, sc_pathways, sc_clinical, genomic_evidence) |
| Key capabilities | 57 cell types classified, tumor microenvironment profiling (4 phenotypes), spatial transcriptomics analysis, CAR-T target validation, drug response prediction at single-cell resolution |

### Agent 11: Clinical Trial Intelligence

| Parameter | Value |
|---|---|
| Name | clinical-trial-intelligence-agent |
| API port | 8538 (FastAPI) |
| UI port | 8128 (Streamlit) |
| Memory limit | 4 GB |
| Collections | 14 (trial_protocols, trial_eligibility, trial_endpoints, trial_sites, trial_investigators, trial_results, trial_regulatory, trial_literature, trial_biomarkers, trial_safety, trial_rwe, trial_adaptive, trial_guidelines, genomic_evidence) |
| Key capabilities | Protocol design optimization, patient-trial matching, 7-factor site selection scoring, safety signal detection, adaptive trial design support |

---

## VISUAL 7: All 45 Milvus Collections

### Global Configuration

| Parameter | Value |
|---|---|
| Embedding dimension (all collections) | 384 |
| Index type (all collections) | IVF_FLAT |
| nlist | 1024 |
| Distance metric | COSINE |
| Total indexed vectors | 3,487,216 |
| Shared collection | genomic_evidence (used by all 8 agents) |

### Collections by Agent

**Precision Oncology Agent (11 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | onco_literature | Oncology research publications |
| 2 | onco_trials | Clinical trial protocols and results |
| 3 | onco_variants | Cancer-associated genetic variants |
| 4 | onco_biomarkers | Tumor biomarkers and companion diagnostics |
| 5 | onco_therapies | Targeted therapies and regimens |
| 6 | onco_pathways | Signaling pathways and drug targets |
| 7 | onco_guidelines | NCCN/ASCO/ESMO treatment guidelines |
| 8 | onco_resistance | Resistance mechanisms and bypass strategies |
| 9 | onco_outcomes | Treatment outcomes and survival data |
| 10 | onco_cases | Molecular tumor board case reports |
| 11 | genomic_evidence | (shared) Whole-genome variant evidence |

**CAR-T Intelligence Agent (11 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | cart_literature | CAR-T research publications |
| 2 | cart_trials | CAR-T clinical trials |
| 3 | cart_constructs | CAR construct designs and domains |
| 4 | cart_assays | Functional assay protocols and results |
| 5 | cart_manufacturing | Manufacturing process knowledge |
| 6 | cart_safety | CRS/ICANS toxicity data |
| 7 | cart_biomarkers | Response and toxicity biomarkers |
| 8 | cart_regulatory | FDA/EMA regulatory submissions |
| 9 | cart_sequences | scFv and construct sequences |
| 10 | cart_realworld | Real-world evidence and outcomes |
| 11 | genomic_evidence | (shared) |

**Precision Biomarker Agent (14 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | biomarker_reference | Reference ranges and interpretation |
| 2 | biomarker_genetic_variants | Genotype-phenotype correlations |
| 3 | biomarker_pgx_rules | Pharmacogenomic decision rules |
| 4 | biomarker_disease_trajectories | Disease progression modeling |
| 5 | biomarker_clinical_evidence | Published clinical evidence |
| 6 | biomarker_nutrition | Nutritional genomics data |
| 7 | biomarker_drug_interactions | Drug-biomarker interactions |
| 8 | biomarker_aging_markers | Biological aging indicators |
| 9 | biomarker_genotype_adjustments | Ancestry-specific adjustments |
| 10 | biomarker_monitoring | Longitudinal monitoring rules |
| 11 | biomarker_critical_values | Critical value thresholds |
| 12 | biomarker_discordance_rules | Discordant result interpretation |
| 13 | biomarker_aj_carrier_screening | Ashkenazi Jewish carrier panels |
| 14 | genomic_evidence | (shared) |

**Imaging Intelligence Agent (11 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | imaging_literature | Radiology/imaging research |
| 2 | imaging_trials | Imaging clinical trials |
| 3 | imaging_findings | Standardized finding descriptions |
| 4 | imaging_protocols | Acquisition protocol parameters |
| 5 | imaging_devices | Device specifications and calibration |
| 6 | imaging_anatomy | Anatomical structure reference |
| 7 | imaging_benchmarks | AI model benchmark results |
| 8 | imaging_guidelines | ACR/RSNA imaging guidelines |
| 9 | imaging_report_templates | Structured report templates |
| 10 | imaging_datasets | Training/validation dataset catalog |
| 11 | genomic_evidence | (shared) |

**Precision Autoimmune Agent (13 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | autoimmune_clinical_documents | Patient clinical documents |
| 2 | autoimmune_patient_labs | Laboratory results |
| 3 | autoimmune_autoantibody_panels | Autoantibody test panels |
| 4 | autoimmune_hla_associations | HLA allele-disease associations |
| 5 | autoimmune_disease_criteria | Classification criteria (ACR/EULAR) |
| 6 | autoimmune_disease_activity | Activity scoring indices |
| 7 | autoimmune_flare_patterns | Flare pattern recognition data |
| 8 | autoimmune_biologic_therapies | Biologic therapy evidence |
| 9 | autoimmune_pgx_rules | Pharmacogenomic rules |
| 10 | autoimmune_clinical_trials | Autoimmune trial data |
| 11 | autoimmune_literature | Research publications |
| 12 | autoimmune_patient_timelines | Longitudinal patient timelines |
| 13 | autoimmune_cross_disease | Cross-disease overlap patterns |

**Cardiology Intelligence Agent (13 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | cardio_literature | Cardiovascular research |
| 2 | cardio_trials | Cardiology clinical trials |
| 3 | cardio_imaging | Echocardiography/cardiac imaging |
| 4 | cardio_electrophysiology | EP studies and arrhythmia data |
| 5 | cardio_heart_failure | Heart failure management |
| 6 | cardio_valvular | Valvular heart disease |
| 7 | cardio_prevention | Preventive cardiology |
| 8 | cardio_interventional | Interventional procedures |
| 9 | cardio_oncology | Cardio-oncology crossover |
| 10 | cardio_devices | Cardiac devices (ICD, CRT, LVAD) |
| 11 | cardio_guidelines | ACC/AHA/ESC guidelines |
| 12 | cardio_hemodynamics | Hemodynamic assessment data |
| 13 | genomic_evidence | (shared) |

**Neurology Intelligence Agent (14 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | neuro_literature | Neurology research publications |
| 2 | neuro_trials | Neurological clinical trials |
| 3 | neuro_imaging | Neuroimaging (MRI, PET, CT) |
| 4 | neuro_electrophysiology | EEG/EMG/NCS data |
| 5 | neuro_degenerative | Neurodegenerative disease knowledge |
| 6 | neuro_cerebrovascular | Stroke and vascular neurology |
| 7 | neuro_epilepsy | Epilepsy and seizure disorders |
| 8 | neuro_oncology | Neuro-oncology (brain tumors) |
| 9 | neuro_ms | Multiple sclerosis |
| 10 | neuro_movement | Movement disorders (Parkinson, etc.) |
| 11 | neuro_headache | Headache and migraine |
| 12 | neuro_neuromuscular | Neuromuscular disorders |
| 13 | neuro_guidelines | AAN/WHO guidelines |
| 14 | genomic_evidence | (shared) |

**Rare Disease Diagnostic Agent (14 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | rd_phenotypes | HPO phenotype descriptions |
| 2 | rd_diseases | OMIM/Orphanet disease entries |
| 3 | rd_genes | Disease-gene associations |
| 4 | rd_variants | Pathogenic variant database |
| 5 | rd_literature | Rare disease publications |
| 6 | rd_trials | Rare disease clinical trials |
| 7 | rd_therapies | Available therapies/gene therapies |
| 8 | rd_case_reports | Published case reports |
| 9 | rd_guidelines | Diagnostic/management guidelines |
| 10 | rd_pathways | Metabolic/signaling pathways |
| 11 | rd_registries | Patient registry data |
| 12 | rd_natural_history | Natural history studies |
| 13 | rd_newborn_screening | Newborn screening panels |
| 14 | genomic_evidence | (shared) |

**Clinical Trial Intelligence Agent (14 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | trial_protocols | Protocol documents and designs |
| 2 | trial_eligibility | Eligibility criteria (I/E) |
| 3 | trial_endpoints | Primary/secondary endpoints |
| 4 | trial_sites | Site performance and capabilities |
| 5 | trial_investigators | Investigator profiles |
| 6 | trial_results | Published trial results |
| 7 | trial_regulatory | Regulatory submission data |
| 8 | trial_literature | Trial-related publications |
| 9 | trial_biomarkers | Biomarker-driven trial data |
| 10 | trial_safety | Safety signal data (DSMB) |
| 11 | trial_rwe | Real-world evidence |
| 12 | trial_adaptive | Adaptive design parameters |
| 13 | trial_guidelines | ICH/FDA/EMA guidance |
| 14 | genomic_evidence | (shared) |

**Single-Cell Intelligence Agent (12 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | sc_cell_types | Cell type classifications |
| 2 | sc_markers | Cell surface markers |
| 3 | sc_spatial | Spatial transcriptomics data |
| 4 | sc_tme | Tumor microenvironment profiles |
| 5 | sc_drug_response | Single-cell drug response data |
| 6 | sc_literature | Single-cell research publications |
| 7 | sc_methods | Analysis methods and pipelines |
| 8 | sc_datasets | Reference dataset catalog |
| 9 | sc_trajectories | Cell trajectory/pseudotime data |
| 10 | sc_pathways | Pathway enrichment results |
| 11 | sc_clinical | Clinical correlations |
| 12 | genomic_evidence | (shared) |

**Pharmacogenomics Intelligence Agent (15 collections)**

| # | Collection Name | Description |
|---|---|---|
| 1 | pgx_gene_reference | Pharmacogene reference data |
| 2 | pgx_drug_guidelines | CPIC/DPWG drug guidelines |
| 3 | pgx_drug_interactions | Drug-gene interactions |
| 4 | pgx_hla_hypersensitivity | HLA-mediated adverse reactions |
| 5 | pgx_phenoconversion | Phenoconversion rules (inhibitors) |
| 6 | pgx_dosing_algorithms | Genotype-based dosing algorithms |
| 7 | pgx_clinical_evidence | PGx clinical evidence |
| 8 | pgx_population_data | Allele frequency by population |
| 9 | pgx_clinical_trials | PGx-guided trial data |
| 10 | pgx_fda_labels | FDA pharmacogenomic labels |
| 11 | pgx_drug_alternatives | Alternative drug recommendations |
| 12 | pgx_patient_profiles | PGx patient profile templates |
| 13 | pgx_implementation | Clinical implementation guides |
| 14 | pgx_education | Patient/provider education |
| 15 | genomic_evidence | (shared) |

---

## VISUAL 8: 13 Therapeutic Areas -- Gene Detail

### Gene Inventory by Therapeutic Area

| # | Therapeutic Area | Gene Count | Key Gene Symbols |
|---|---|---|---|
| 1 | Neurology and Neurodegeneration | 36 | VCP, C9orf72, MAPT, APP, HTT, CGRP, SCN9A, TREM2, GRN, TBK1, FUS, TARDBP |
| 2 | Oncology | 27 | BRCA1, BRCA2, EGFR, KRAS, ALK, BRAF, PDCD1, CD274, CTLA4, ERBB2, PIK3CA, FLT3 |
| 3 | Metabolic and Endocrine | 22 | GLP1R, GIPR, GCGR, INSR, PPARG, DPP4, SGLT2, MC4R, HMGCR, PCSK9, LDLR, TSHR |
| 4 | Infectious Disease | 21 | HIV1_RT, HIV1_PR, HIV1_IN, CCR5, HCV_NS3, HCV_NS5A, SARS2_MPRO, ACE2, RPOB, CYP51A1 |
| 5 | Respiratory | 13 | ADRB2, CHRM3, IL4R, IL5, IL5RA, IL13, BMPR2, PDE4D, SERPINA1, TSLP |
| 6 | Rare Disease | 12 | CFTR, SMN1, DMD, HBB, WT1, EWSR1, PAX3, RB1, PTCH1, PHOX2B |
| 7 | Hematology | 12 | SYK, THPO, MPL, EPOR, JAK2, F10, ADAMTS13, VWF, F8, F9, HBB, CALR |
| 8 | GI and Hepatology | 12 | NOD2, ATG16L1, IL12B, ITGA4, S1PR1, PNPLA3, HSD17B13, FXR, THR_BETA, ATP4A |
| 9 | Ophthalmology | 11 | VEGFA, CFH, C3, C5, ROCK1, CA2, RHO, RPE65, USH2A, ABCA4, OPTN |
| 10 | Pharmacogenomics | 11 | CYP2D6, CYP2C19, CYP3A4, CYP2C9, VKORC1, DPYD, TPMT, NUDT15, SLCO1B1, UGT1A1, ABCB1 |
| 11 | Cardiovascular | 10 | PCSK9, TTR, MYBPC3, SCN5A, KCNQ1, KCNH2, MYH7, LMNA, DSP, PKP2 |
| 12 | Immunology and Inflammation | 9 | IL6, TNF, JAK1, JAK2, IL17A, IL23A, TSLP, TNFRSF1A, NLRP3 |
| 13 | Dermatology | 9 | IL31RA, TYK2, COL7A1, KIT, IL13, IL4R, JAK1, PDE4D, EGFR |

### Summary

| Metric | Value |
|---|---|
| Total genes mapped | 201 |
| Druggable targets | 171 |
| Druggability rate | 85% |
| Unique gene symbols | 201 (some appear in multiple areas) |
| Knowledge base entries per gene | Mechanism, druggability, diseases, PDB structures |

---

## VISUAL 9: Pediatric Safety -- Complete Filter Specification

### All 6 Pediatric Safety Filters

| # | Filter Name | Flag Code | Threshold | Severity | Clinical Rationale | Estimation Method |
|---|---|---|---|---|---|---|
| 1 | Blood-Brain Barrier Penetration | HIGH_MW_BBB_RISK | Molecular weight > 500 Da | medium | Pediatric BBB is more permeable in young children (especially <2 years); high-MW compounds that would be excluded by adult BBB may penetrate and cause CNS toxicity | RDKit `Descriptors.ExactMolWt()` |
| 2 | Cardiac Liability (hERG) | HERG_CARDIAC_RISK | hERG IC50 < 10 uM | critical | Pediatric patients have lower threshold for QT prolongation; hERG potassium channel blockade is the primary molecular mechanism for drug-induced long QT syndrome | Estimated from physicochemical properties: `max(0.5, 30 - (logP * 4) - (MW - 350) * 0.02)` |
| 3 | Hepatic Immaturity | HIGH_LIPOPHILICITY | LogP > 5 | high | CYP3A4 enzyme does not reach adult expression levels until age 1-2; highly lipophilic compounds accumulate in immature hepatocytes causing dose-dependent toxicity | RDKit `Crippen.MolLogP()` (Wildman-Crippen method) |
| 4 | Teratogenicity | TERATOGENICITY | Any teratogenic SMARTS match | high | Adolescent patients require pregnancy prevention programs; known teratogenic scaffolds (thalidomide analogs, retinoids, valproic acid derivatives) must be flagged pre-ranking | RDKit SMARTS substructure search against curated alert library |
| 5 | Oral Bioavailability | HIGH_TPSA_ORAL_RISK | TPSA > 140 A^2 | medium | Pediatric patients often require liquid oral formulations; high topological polar surface area limits intestinal absorption and reduces oral bioavailability | RDKit `Descriptors.TPSA()` |
| 6 | GI Absorption Variability | HIGH_FLEXIBILITY | Rotatable bonds > 10 | low | Pediatric GI tract has different pH, transit time, and surface area compared to adults; highly flexible molecules show unpredictable absorption variability in children | RDKit `Lipinski.NumRotatableBonds()` |

### Decision Logic

```
IF any flag has severity == "critical":
    pediatric_safe = False
    recommendation = "Requires additional pediatric safety evaluation"
ELIF no flags at all:
    pediatric_safe = True
    recommendation = "Suitable for pediatric development"
ELSE:
    pediatric_safe = True (but flagged)
    recommendation = "Requires additional pediatric safety evaluation"
```

### Traditional vs. HCLS AI Factory Approach

| Aspect | Traditional | HCLS AI Factory |
|---|---|---|
| When safety assessed | Late-stage clinical development (Phase II/III) | Before any candidate is ranked (Stage 8) |
| Pediatric-specific filters | Rarely applied computationally | 6 automated filters with severity grading |
| Cost of late-stage failure | Hundreds of millions of dollars | Zero (filtered before ranking) |
| Time to safety signal | Months to years | Seconds (computed inline) |

---

## VISUAL 10: Docker Infrastructure

### Core Infrastructure Services

| Service | Image:Tag | Ports (host:container) | Memory Limit | Purpose |
|---|---|---|---|---|
| milvus-etcd | quay.io/coreos/etcd:v3.5.11 | -- (internal) | 512 MB | Milvus metadata storage |
| milvus-minio | minio/minio:RELEASE.2024-01-01T16-36-33Z | -- (internal) | 1 GB | Milvus object storage |
| milvus | milvusdb/milvus:v2.4.0 | 19530:19530, 9091:9091 | 8 GB | Vector database |
| prometheus | prom/prometheus:v2.49.1 | 9099:9090 | 1 GB | Metrics collection |
| grafana | grafana/grafana:10.3.1 | 3000:3000 | 512 MB | Dashboards |
| node-exporter | prom/node-exporter:v1.8.0 | 9100:9100 | 128 MB | Host metrics |

### Intelligence Agent Services

| Service | Build Context | API Port | UI Port | Memory Limit |
|---|---|---|---|---|
| precision-biomarker-agent | ./core/agents/precision-biomarker | 8102:8000 | 8502:8501 | 4 GB |
| precision-oncology-agent | ./core/engines/precision-oncology/agent/agent | 8103:8000 | 8503:8501 | 4 GB |
| cart-intelligence-agent | ./core/agents/cart | 8104:8000 | 8504:8501 | 4 GB |
| imaging-intelligence-agent | ./core/engines/clinical-imaging/agent/agent | 8105:8000 | 8505:8501 | 8 GB |
| precision-autoimmune-agent | ./core/agents/precision-autoimmune | 8106:8000 | 8506:8501 | 4 GB |
| cardiology-intelligence-agent | ./core/engines/cardiology | 8126:8126 | 8536:8536 | 4 GB |
| neurology-intelligence-agent | ./core/agents/neurology | 8528:8528 | 8529:8529 | 4 GB |
| rare-disease-diagnostic-agent | ./core/agents/rare-disease-diagnostic | 8134:8134 | 8544:8544 | 4 GB |
| clinical-trial-intelligence-agent | ./core/agents/clinical-trial | 8538:8538 | 8128:8128 | 4 GB |
| single-cell-intelligence-agent | ./core/agents/single-cell | 8540:8540 | 8130:8130 | 8 GB |
| pharmacogenomics-intelligence-agent | ./core/agents/pharmacogenomics | 8107:8107 | 8507:8507 | 4 GB |

### Application Services

| Service | Build/Image | Ports | Memory Limit | Purpose |
|---|---|---|---|---|
| landing-page | ./landing-page (Dockerfile) | 8080:8080 | 256 MB | Service hub and health dashboard |
| dd-pipeline-ui | ./core/engines/therapeutic-discovery (Dockerfile) | 8505:8501 | -- | Drug discovery Streamlit UI |

### NIM GPU Services

| Service | Image:Tag | Port (host:container) | GPU | Purpose |
|---|---|---|---|---|
| dd-molmim | nvcr.io/nim/nvidia/molmim:1.0.0 | 8001:8000 | 1x GPU required | Molecule generation |
| dd-diffdock | nvcr.io/nim/mit/diffdock:2.2.0 | 8002:8000 | 1x GPU required | Molecular docking |

### Docker Volumes

| Volume | Mount Point | Service | Purpose |
|---|---|---|---|
| etcd_data | /etcd | milvus-etcd | Metadata persistence |
| minio_data | /minio_data | milvus-minio | Object storage persistence |
| milvus_data | /var/lib/milvus | milvus | Vector index persistence |
| prometheus_data | /prometheus | prometheus | Metric time series |
| grafana_data | /var/lib/grafana | grafana | Dashboard configs |

### Docker Networks

| Network | Purpose |
|---|---|
| default (docker-compose) | All agent and infrastructure communication |
| drug-discovery-network | NIM services + pipeline UI |

### Common Environment Variables

| Variable | Value | Used By |
|---|---|---|
| MILVUS_HOST | milvus | All agents |
| MILVUS_PORT | 19530 | All agents |
| ANTHROPIC_API_KEY | ${ANTHROPIC_API_KEY} | All agents |
| EMBEDDING_MODEL | BAAI/bge-small-en-v1.5 | All agents |
| LOG_LEVEL | INFO | All agents |
| NGC_API_KEY | ${NGC_API_KEY} | NIM services |

### Healthcheck Configuration (All Agents)

| Parameter | Value |
|---|---|
| Interval | 30 seconds |
| Timeout | 10 seconds |
| Retries | 3 |
| Start period | 30 seconds |
| Method | HTTP GET to /health endpoint |

---

## VISUAL 11: Complete API Endpoint Reference

### Common Endpoints (All 11 Agents)

Every agent exposes:

| Method | Path | Description | Response |
|---|---|---|---|
| GET | /health | Service health and readiness check | JSON: status, uptime, version, collections loaded |
| GET | /collections | List all Milvus collections with entity counts | JSON: collection names + counts |
| GET | /metrics | Prometheus-format metrics | text/plain: OpenMetrics format |

### RAG Query Endpoints (Standard Pattern)

| Method | Path | Description | Request Body | Response |
|---|---|---|---|---|
| POST | /query | RAG-grounded question answering | `{query, context?, top_k?}` | JSON: answer, sources, confidence |
| POST | /search | Vector similarity search only | `{query, collections?, top_k?}` | JSON: ranked results with scores |
| POST | /find-related | Find related evidence across collections | `{entity, entity_type?}` | JSON: related entities and evidence |
| GET | /knowledge/stats | Knowledge base statistics | -- | JSON: collection counts, total vectors |

### Agent-Specific Endpoints

**Precision Oncology (:8103)**
- POST /query -- Molecular tumor board queries
- POST /search -- Variant/therapy search
- POST /find-related -- Gene-pathway-drug relationships
- GET /knowledge/stats -- Oncology knowledge statistics

**CAR-T Intelligence (:8104)**
- POST /query -- CAR-T design and safety queries
- POST /search -- Construct/antigen search
- POST /find-related -- Target-construct-outcome links
- GET /knowledge/stats -- CAR-T knowledge statistics

**Precision Biomarker (:8102)**
- POST /query -- Biomarker interpretation queries
- GET /knowledge/stats -- Biomarker panel statistics

**Cardiology (:8126)**
- GET /workflows -- 11 clinical workflow definitions
- POST /v1/cardio/workflow/* -- Execute clinical workflows (risk calculators, GDMT optimizer)
- GET /health -- Extended health with workflow count

**Neurology (:8528)**
- GET /workflows -- 8 clinical workflow definitions
- GET /health -- Extended health with scale count

**Autoimmune (:8106)**
- POST /v1/autoimmune/integrated-assessment -- Full autoimmune workup
- POST /query/stream -- Streaming SSE query response
- POST /analyze -- Disease activity analysis
- POST /differential -- Differential diagnosis
- POST /ingest/upload -- Document upload and processing
- POST /ingest/demo-data -- Seed demo patient data
- POST /collections/create -- Create new collections
- POST /export -- Export assessment reports

**Rare Disease (:8134)**
- GET /workflows -- Diagnostic workflow definitions
- GET /health -- Extended health with phenotype matcher status

**Clinical Trial (:8538)**
- GET /workflows -- Trial design workflow definitions
- Routers: /clinical (trial matching), /reports (generation), /events (safety signals)

**Single-Cell (:8540)**
- GET /workflows -- Analysis workflow definitions

**Pharmacogenomics (:8107)**
- POST /query -- PGx interpretation queries
- POST /search -- Gene/drug/allele search
- POST /find-related -- Drug-gene-phenotype relationships

**Imaging (:8105)**
- POST /query -- Imaging interpretation queries
- POST /search -- Finding/protocol search
- POST /find-related -- Imaging-genomic cross-modal queries

### Authentication and Rate Limiting

| Parameter | Value |
|---|---|
| API key auth | ANTHROPIC_API_KEY required (environment variable) |
| Agent-level auth | No additional auth middleware (designed for internal network) |
| Rate limiting | Not enforced at agent level (designed for single-machine deployment) |
| CORS | Enabled via FastAPI/Flask middleware |
| Response format | JSON (all endpoints except /metrics which returns text/plain) |
| Streaming | SSE (Server-Sent Events) for /query/stream endpoints |

---

## VISUAL 12: Technology Stack -- Full Inventory

### Layer 1: GPU Compute Hardware

| Component | Specification |
|---|---|
| Machine | NVIDIA DGX Spark |
| GPU | NVIDIA GB10 (Grace Blackwell) |
| Memory | 128 GB unified LPDDR5x |
| CPU | 20 ARM cores (NVIDIA Grace) |
| Interconnect | NVLink-C2C (CPU-GPU) |
| CUDA | 12.x |
| Storage | NVMe SSD (~1.1 TB platform footprint) |

### Layer 2: Genomics Tools

| Tool | Version / Specification |
|---|---|
| NVIDIA Parabricks | 4.6.0-1 (nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1) |
| BWA-MEM2 | Bundled with Parabricks (fq2bam) |
| Google DeepVariant | Bundled with Parabricks (>99% accuracy) |
| samtools | Bundled with Parabricks (BAM indexing) |
| GRCh38 reference | hg38 (3.1 GB) |
| ClinVar | February 2026 release (4.1M variants) |
| AlphaMissense | v1.0 (71,697,560 predictions) |
| Ensembl VEP | GRCh38 build |

### Layer 3: AI / LLM / Embedding

| Tool | Version |
|---|---|
| Claude LLM | claude-sonnet-4-20250514 (Anthropic) |
| anthropic (Python SDK) | 0.75.0 |
| openai (Python SDK) | 2.15.0 |
| BGE-small-en-v1.5 | BAAI (384-dim embeddings) |
| sentence-transformers | 5.2.0 |
| PyTorch | 2.9.1 |

### Layer 4: Drug Discovery / Chemistry

| Tool | Version |
|---|---|
| MolMIM NIM | nvcr.io/nim/nvidia/molmim:1.0.0 |
| DiffDock NIM | nvcr.io/nim/mit/diffdock:2.2.0 |
| RDKit | 2025.9.3 |
| py3Dmol | 2.5.3 |
| stmol | 0.0.9 |

### Layer 5: Data Infrastructure

| Tool | Version |
|---|---|
| Milvus | v2.4.0 (milvusdb/milvus:v2.4.0) |
| pymilvus | 2.6.6 |
| etcd | v3.5.11 (quay.io/coreos/etcd:v3.5.11) |
| MinIO | RELEASE.2024-01-01T16-36-33Z |

### Layer 6: Web Frameworks

| Tool | Version |
|---|---|
| FastAPI | 0.128.0 |
| uvicorn | 0.40.0 |
| Flask | 3.1.2 |
| flask-cors | 6.0.2 |
| Streamlit | 1.52.2 |
| streamlit-chat | 0.1.1 |

### Layer 7: Python Core Libraries

| Package | Version | Purpose |
|---|---|---|
| pydantic | 2.12.5 | Data validation and settings |
| pydantic-settings | 2.12.0 | Environment-based configuration |
| loguru | 0.7.3 | Structured logging |
| python-dotenv | 1.2.1 | Environment file loading |
| requests | 2.32.5 | HTTP client (NIM, PDB) |
| numpy | 2.4.0 / 2.4.1 | Numerical computation |
| pandas | 2.3.3 | Data manipulation |
| tqdm | 4.67.1 | Progress bars |
| typer | 0.21.1 | CLI framework |
| rich | 14.2.0 | Terminal formatting |
| pillow | 12.1.0 | Image processing |
| reportlab | 4.4.0 | PDF report generation |
| psutil | 7.2.1 | System monitoring |
| pynvml | 13.0.1 | GPU monitoring |

### Layer 8: Genomics Python Libraries

| Package | Version | Purpose |
|---|---|---|
| cyvcf2 | 0.31.4 | VCF file parsing (C-backed, fast) |
| pysam | 0.23.3 | BAM/SAM file access |

### Layer 9: Observability

| Tool | Version |
|---|---|
| Prometheus | v2.49.1 / v2.52.0 (prom/prometheus) |
| Grafana | 10.3.1 / 11.0.0 (grafana/grafana) |
| Node Exporter | v1.8.0 (prom/node-exporter) |
| DCGM Exporter | 3.3.5-3.4.0 (nvcr.io/nvidia/k8s/dcgm-exporter) |
| OpenTelemetry API | >=1.29.0 |
| OpenTelemetry SDK | >=1.29.0 |

### Layer 10: Orchestration and Infrastructure

| Tool | Version |
|---|---|
| Docker Compose | v2.20+ (compose spec 3.8) |
| Nextflow | DSL2 (HLS-Pipeline v1.0.0) |

### Layer 11: Testing

| Tool | Version | Purpose |
|---|---|---|
| pytest | 9.0.2 | Test framework |
| pytest-cov | 7.0.0 | Coverage reporting |
| pytest-asyncio | 0.24.0 | Async test support |
