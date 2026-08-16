# HCLS AI Factory

## Patient DNA to Drug Candidates in <5 Hours

## Engine 1: Genomic Foundation

### Input: FASTQ (R1 + R2, ~200 GB Paired-End, 30x WGS)
#### Benchmark Genome: GIAB HG002 (Ashkenazi Son)
#### Sequencing: Illumina Short-Read, Paired-End

### Parabricks 4.6 (GPU-Accelerated)
#### Stage 1: QC and Preprocessing
##### Adapter Trimming + Quality Filtering
##### GPU-Accelerated FastQC-Equivalent Metrics
#### Stage 2: fq2bam -- BWA-MEM2 Alignment
##### Reference: GRCh38/hg38 (with Alt Contigs, ~3.2 GB)
##### In-Stream Duplicate Marking
##### Output: Sorted, Indexed BAM
##### Performance: 45-90 min GPU vs 8-16 hrs CPU (8-12x Speedup)
#### Stage 3: DeepVariant -- Variant Calling (v1.6+)
##### Deep Learning SNP + Indel Calling
##### Trained on Genome in a Bottle Truth Sets
##### Outputs Per-Sample gVCF
##### Performance: 30-60 min GPU vs 6-12 hrs CPU (6-10x Speedup)
#### Stage 4: Structural Variant Detection
##### Deletions, Duplications, Inversions, Translocations
##### CNV Analysis with Read-Depth Normalization
#### Stage 5: samtools/htslib -- BAM Indexing + QC
##### BAM/VCF File Manipulation
##### bcftools VCF Processing and Filtering
#### Stage 6: Reference Genome
##### GRCh38 Human Genome Assembly
##### Alt Contigs Included

### Output: Annotated VCF (11.7M Variant Records)

### Annotation
#### ClinVar (4.1M Clinical Records)
#### AlphaMissense (71M Missense Pathogenicity Predictions)
#### gnomAD v4 (Population Allele Frequencies)
#### CADD, REVEL, SpliceAI Functional Impact Scores
#### Gene-Disease Association Mapping
#### Result: 3.56M Annotated Searchable Variants

### Performance
#### Total Pipeline: 2-4 hrs GPU vs 24-48 hrs CPU (6-12x Speedup)
#### Alignment: 45-90 min GPU vs 8-16 hrs CPU
#### Variant Calling: 30-60 min GPU vs 6-12 hrs CPU
#### Annotation: 15-30 min GPU vs 2-4 hrs CPU
#### Working Storage per Sample: ~200 GB Temporary, ~50 GB Final
#### Annotation Databases: ~50 GB (ClinVar + AlphaMissense + gnomAD)

---

## Engine 2: Precision Intelligence Engine

### Milvus Vector Database (v2.4.0, Standalone Mode)

#### Infrastructure
##### etcd v3.5.11 (:2379, 512 MB Memory, 4 GB Backend Quota)
##### MinIO RELEASE.2024-01-01 (:9000, Console :9001, 1 GB Memory)
##### Milvus gRPC (:19530, Health :9091, 8 GB Memory)

#### Index Configuration
##### Index Type: IVF_FLAT
##### Similarity Metric: COSINE
##### Embedding Dimension: 384
##### nlist: 128

#### 45 Collections (~4,948,660 Total Vectors)

##### Shared Collection (Read-Only, All Agents)
###### genomic_evidence (3,560,000 vectors) -- Annotated Variants from VCF Pipeline

##### Biomarker Agent Collections (13 Own)
###### biomarker_reference
###### biomarker_genetic_variants
###### biomarker_pgx_rules
###### biomarker_disease_trajectories
###### biomarker_clinical_evidence
###### biomarker_nutrition
###### biomarker_drug_interactions
###### biomarker_aging_markers
###### biomarker_genotype_adjustments
###### biomarker_monitoring
###### biomarker_critical_values
###### biomarker_discordance_rules
###### biomarker_aj_carrier_screening
###### Subtotal: ~7,650 Vectors

##### Oncology Agent Collections (10 Own)
###### onco_literature
###### onco_trials
###### onco_variants
###### onco_biomarkers
###### onco_therapies
###### onco_pathways
###### onco_guidelines
###### onco_resistance
###### onco_outcomes
###### onco_cases
###### Subtotal: ~11,900 Vectors

##### CAR-T Agent Collections (10 Own)
###### cart_literature
###### cart_trials
###### cart_constructs
###### cart_assays
###### cart_manufacturing
###### cart_safety
###### cart_biomarkers
###### cart_regulatory
###### cart_sequences
###### cart_realworld
###### Subtotal: ~5,750 Vectors

##### Imaging Agent Collections (10 Own)
###### imaging_literature
###### imaging_trials
###### imaging_findings
###### imaging_protocols
###### imaging_devices
###### imaging_anatomy
###### imaging_benchmarks
###### imaging_guidelines
###### imaging_report_templates
###### imaging_datasets
###### Subtotal: ~6,350 Vectors

##### Autoimmune Agent Collections (13 Own)
###### autoimmune_clinical_documents
###### autoimmune_patient_labs
###### autoimmune_autoantibody_panels
###### autoimmune_hla_associations
###### autoimmune_disease_criteria
###### autoimmune_disease_activity
###### autoimmune_flare_patterns
###### autoimmune_biologic_therapies
###### autoimmune_pgx_rules
###### autoimmune_clinical_trials
###### autoimmune_literature
###### autoimmune_patient_timelines
###### autoimmune_cross_disease
###### Subtotal: ~8,750 Vectors

##### Pharmacogenomics Agent Collections (14 Own)
###### pgx_gene_reference
###### pgx_drug_guidelines
###### pgx_drug_interactions
###### pgx_hla_hypersensitivity
###### pgx_phenoconversion
###### pgx_dosing_algorithms
###### pgx_clinical_evidence
###### pgx_population_data
###### pgx_clinical_trials
###### pgx_fda_labels
###### pgx_drug_alternatives
###### pgx_patient_profiles
###### pgx_implementation
###### pgx_education
###### Subtotal: ~7,350 Vectors

##### Cardiology Agent Collections (12 Own)
###### cardio_literature
###### cardio_trials
###### cardio_imaging
###### cardio_electrophysiology
###### cardio_heart_failure
###### cardio_valvular
###### cardio_prevention
###### cardio_interventional
###### cardio_oncology
###### cardio_devices
###### cardio_guidelines
###### cardio_hemodynamics
###### Subtotal: ~4,830 Vectors

##### Clinical Trial Agent Collections (13 Own)
###### trial_protocols
###### trial_eligibility
###### trial_endpoints
###### trial_sites
###### trial_investigators
###### trial_results
###### trial_regulatory
###### trial_literature
###### trial_biomarkers
###### trial_safety
###### trial_rwe
###### trial_adaptive
###### trial_guidelines
###### Subtotal: ~151,500 Vectors

##### Rare Disease Agent Collections (13 Own)
###### rd_phenotypes
###### rd_diseases
###### rd_genes
###### rd_variants
###### rd_literature
###### rd_trials
###### rd_therapies
###### rd_case_reports
###### rd_guidelines
###### rd_pathways
###### rd_registries
###### rd_natural_history
###### rd_newborn_screening
###### Subtotal: ~624,580 Vectors

##### Neurology Agent Collections (13 Own)
###### neuro_literature
###### neuro_trials
###### neuro_imaging
###### neuro_electrophysiology
###### neuro_degenerative
###### neuro_cerebrovascular
###### neuro_epilepsy
###### neuro_oncology
###### neuro_ms
###### neuro_movement
###### neuro_headache
###### neuro_neuromuscular
###### neuro_guidelines
###### Subtotal: ~355,000 Vectors

##### Single-Cell Agent Collections (11 Own)
###### sc_cell_types
###### sc_markers
###### sc_spatial
###### sc_tme
###### sc_drug_response
###### sc_literature
###### sc_methods
###### sc_datasets
###### sc_trajectories
###### sc_pathways
###### sc_clinical
###### Subtotal: ~205,000 Vectors

### Embedding Model: BGE-small-en-v1.5 (BAAI)
#### Dimension: 384
#### Batch Size: 32
#### Instruction Prefix: "Represent this sentence for searching relevant passages: "

### Claude LLM (Anthropic)
#### Model: claude-sonnet-4-20250514
#### Max Tokens: 4,096
#### Temperature: 0.2
#### Streaming: SSE (Server-Sent Events)
#### Max Retries: 3
#### Multi-Turn Conversation Memory: 3 Exchanges

### RAG Search Configuration
#### Top-K Per Collection: 5
#### Score Threshold: 0.40
#### Citation High Threshold: 0.75
#### Citation Medium Threshold: 0.60

### Knowledge Graph
#### 201 Genes | 171 Druggable (85%)
#### 13 Therapeutic Areas
#### Mapping: Variant -> Gene -> Protein -> Pathway -> Disease -> Drug
#### Cross-References: HGNC Gene Symbols, OMIM IDs, NCT IDs, HPO Terms, ClinVar Accessions, PubMed IDs

### 11 Intelligence Agents

#### Precision Biomarker Agent
##### API: :8529 (FastAPI/Uvicorn)
##### UI: :8528 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 14 (13 Own + genomic_evidence)
##### Endpoints: /health, /collections, /knowledge/stats, /metrics
##### Capabilities
###### Biological Age Estimation
###### Disease Trajectory Detection
###### Pharmacogenomic Profiling
###### FHIR R4 Report Export
###### Ashkenazi Jewish Carrier Screening
###### Critical Value Detection
###### Discordance Rules Analysis
###### Cross-Agent: Oncology, CAR-T, PGx, Cardiology, Trial

#### Precision Oncology Agent
##### API: :8527 (FastAPI/Uvicorn)
##### UI: :8526 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 11 (10 Own + genomic_evidence)
##### Endpoints: /health, /collections, /query, /search, /find-related, /knowledge/stats, /metrics
##### Capabilities
###### Molecular Tumor Board Support
###### 80+ Targeted Therapy Matching
###### Resistance Mechanism Prediction
###### Evidence-Tiered Reporting (Level I-IV)
###### Cross-Modal Enrichment (Genomic + Imaging)
###### CIViC Integration (civicdb.org/api)
###### Cross-Agent: CAR-T, Biomarker, Trial, Cardiology, Neurology, PGx, Imaging, Single-Cell

#### CAR-T Intelligence Agent
##### API: :8522 (FastAPI/Uvicorn)
##### UI: :8521 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 11 (10 Own + genomic_evidence)
##### Endpoints: /health, /collections, /query, /search, /find-related, /knowledge/stats, /metrics
##### Capabilities
###### 34 Target Antigens Profiled
###### CRS/ICANS Toxicity Prediction
###### Manufacturing Optimization
###### Costimulatory Domain Analysis (4-1BB vs CD28)
###### Regulatory Intelligence
###### Real-World Evidence Integration
###### Cross-Agent: Oncology, Biomarker, Single-Cell, Cardiology, Trial

#### Imaging Intelligence Agent
##### API: :8524 (FastAPI/Uvicorn)
##### UI: :8525 (Streamlit)
##### Docker Memory: 8 GB
##### Collections: 11 (10 Own + genomic_evidence)
##### Endpoints: /health, /collections, /query, /search, /find-related, /knowledge/stats, /metrics
##### NIM Integration
###### Vista3D (:8530) -- 127 Anatomical Structures
###### MAISI (:8531) -- Synthetic CT Generation
###### VILA-M3 (:8532) -- Visual Question Answering
##### DICOM: Orthanc Server (:8042), OHIF Viewer (:8526)
##### Capabilities
###### DICOM Auto-Ingestion (5s Poll Interval)
###### 3D Volume Preview (MP4, 8 FPS, Max 200 Frames)
###### Cross-Agent: Oncology, Trial, Cardiology, Neurology

#### Precision Autoimmune Agent
##### API: :8532 (FastAPI/Uvicorn)
##### UI: :8531 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 14 (13 Own + genomic_evidence)
##### Endpoints: /health, /healthz, /v1/autoimmune/integrated-assessment, /query, /query/stream, /search, /analyze, /differential, /ingest/upload, /ingest/demo-data, /collections, /collections/create, /export, /metrics
##### Capabilities
###### 13 Conditions + Overlap Syndromes
###### Disease Activity Scoring (DAS28, SLEDAI, BVAS, BASDAI)
###### Autoantibody Panel Interpretation
###### Flare Prediction (0-1 Scale: Imminent >0.8, High >0.6, Moderate >0.4)
###### HLA Association Mapping
###### Biologic Therapy Selection
###### PDF Upload + Document Processing (Max 50 MB, 2500-char Chunks, 200 Overlap)
###### Cross-Agent: Oncology, Cardiology, Neurology, PGx, Biomarker, Trial

#### Pharmacogenomics Agent
##### API: :8107 (FastAPI/Uvicorn)
##### UI: :8507 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 15 (14 Own + genomic_evidence)
##### Endpoints: /health, /collections, /query, /search, /find-related, /knowledge/stats, /metrics
##### Capabilities
###### 14 Core Pharmacogenes (CYP2D6, CYP2C19, CYP3A4, CYP2C9, VKORC1, DPYD, TPMT, UGT1A1, SLCO1B1, NUDT15, HNF1A)
###### 9 CPIC Dosing Algorithms
###### HLA Allele Hypersensitivity Risk Assessment
###### Metabolizer Phenotype Classification
###### Drug Interaction Detection
###### Drug Alternative Recommendations
###### FDA Label Integration
###### Phenoconversion Analysis
###### Cross-Agent: Oncology, Cardiology, Neurology, Trial

#### Cardiology Intelligence Agent
##### API: :8126 (FastAPI/Uvicorn)
##### UI: :8536 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 13 (12 Own + genomic_evidence)
##### Endpoints: /health, /collections, /workflows, /metrics
##### Capabilities
###### 6 Risk Calculators (ASCVD, MAGGIC, HEART, CHA2DS2-VASc, HAS-BLED, TIMI)
###### 11 Clinical Workflows
###### GDMT Optimizer (HFrEF/HFpEF)
###### 56 Cardiac Genes Mapped
###### Electrophysiology Analysis
###### Hemodynamic Assessment
###### Cardio-Oncology Support
###### Cross-Agent: Oncology, Trial, Biomarker, Neurology, Imaging

#### Neurology Intelligence Agent
##### API: :8536 (FastAPI/Uvicorn)
##### UI: :8535 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 14 (13 Own + genomic_evidence)
##### Endpoints: /health, /collections, /workflows, /metrics
##### Capabilities
###### 10 Clinical Scales (NIHSS, GCS, MoCA, UPDRS, EDSS, ALSFRS-R, mRS, Hunt-Hess, BIMS, CAM)
###### 8 Clinical Workflows
###### Multi-Guideline Integration (AAN, WHO, EAN)
###### Multi-Format Reports (PDF, FHIR)
###### Sub-Specialty Coverage: Degenerative, Cerebrovascular, Epilepsy, Oncology, MS, Movement, Headache, Neuromuscular
###### Cross-Agent: Genomics, Imaging, Cardiology, Biomarker, Trial, Rare Disease

#### Rare Disease Diagnostic Agent
##### API: :8134 (FastAPI/Uvicorn)
##### UI: :8544 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 14 (13 Own + genomic_evidence)
##### Endpoints: /health, /collections, /workflows, /metrics
##### Capabilities
###### 88 Diseases Covered
###### ACMG/AMP Variant Classification (23 Criteria)
###### Gene Therapy Eligibility Assessment (25 Therapies)
###### HPO Phenotype Matching
###### Orphanet Integration
###### Newborn Screening Panel
###### Natural History Data
###### Patient Registry Linkage
###### Cross-Agent: Genomics, PGx, Cardiology, Biomarker, Trial, Imaging

#### Clinical Trial Intelligence Agent
##### API: :8538 (FastAPI/Uvicorn)
##### UI: :8128 (Streamlit)
##### Docker Memory: 4 GB
##### Collections: 14 (13 Own + genomic_evidence)
##### Endpoints: /health, /collections, /workflows, /metrics
##### Capabilities
###### Protocol Design Optimization
###### Patient-Trial Matching
###### Site Selection (7-Factor Score)
###### Safety Signal Detection
###### Adaptive Trial Design
###### Real-World Evidence Integration
###### Regulatory Intelligence
###### Biomarker-Driven Enrichment
###### ClinicalTrials.gov API v2 Integration
###### Cross-Agent: Oncology, PGx, Cardiology, Biomarker, Rare Disease, Neurology, Single-Cell, Imaging

#### Single-Cell Intelligence Agent
##### API: :8540 (FastAPI/Uvicorn)
##### UI: :8130 (Streamlit)
##### Docker Memory: 8 GB
##### GPU Memory Limit: 120 GB
##### Collections: 12 (11 Own + genomic_evidence)
##### Endpoints: /health, /collections, /workflows, /metrics
##### Capabilities
###### 57 Cell Types Classified
###### TME Profiling (4 Immune Phenotypes)
###### Spatial Transcriptomics
###### CAR-T Target Validation
###### Drug Response Prediction
###### Cell Trajectory Analysis
###### CellxGene API Integration
###### Cross-Agent: Genomics, Biomarker, Oncology, Trial

---

## Engine 3: Therapeutic Discovery

### 10-Stage Pipeline
#### Stage 0: Initialize -- Load Config, Check NIM Services
#### Stage 1: Normalize Target -- HGNC/UniProt Validation, Gene Symbol Enrichment
#### Stage 2: Structure Discovery -- PDB/EMDB Fetch, RCSB Auto-Load, Cryo-EM Evidence
#### Stage 3: Structure Prep -- Remove Water/Ions, Add Hydrogens, Identify Binding Site
#### Stage 4: Molecule Generation -- MolMIM (Seed SMILES -> Novel Candidates)
#### Stage 5: Chemistry QC -- Lipinski Rule of Five, QED Score, Drug-Likeness Filters
#### Stage 6: 3D Conformer Generation -- RDKit EmbedMolecule + MMFF Optimization
#### Stage 7: Molecular Docking -- DiffDock Protein-Ligand Binding Prediction
#### Stage 8: Composite Ranking -- Weighted Score (Docking 0.4, Generation 0.3, QED 0.3) + Pediatric Safety
#### Stage 9: PDF Report -- JSON + Provenance, Candidate Summary, Stage Timings

### MolMIM (BioNeMo NIM)
#### Container: nvcr.io/nim/nvidia/molmim:1.0.0
#### Cloud API: https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate
#### Local NIM Port: :8001
#### Fallback: RDKit Atom-Swap Generation (MockMolMIMClient)
#### Default: 50 Molecules per Run, Diversity 0.3

### DiffDock (BioNeMo NIM)
#### Container: nvcr.io/nim/mit/diffdock:2.2.0
#### Cloud API: https://health.api.nvidia.com/v1/biology/mit/diffdock
#### Local NIM Port: :8002
#### Poses per Molecule: 10
#### Output: Docking Score, Confidence, RMSD, Hydrogen Bonds, Contact Residues

### RDKit (Chemistry Engine)
#### Lipinski Rule of Five (MW <= 550, LogP <= 5.0, HBD <= 5, HBA <= 10, Max 1 Violation)
#### QED Score (Quantitative Estimate of Drug-Likeness, 0-1)
#### Synthetic Accessibility Score
#### MMFF Conformer Optimization (AllChem.MMFFOptimizeMolecule)
#### Molecular Descriptors: MW, LogP, TPSA, Rotatable Bonds, HBD, HBA
#### SMARTS Structural Alerts (Glutarimide/Thalidomide Class, Retinoid Scaffolds)

### 6 Pediatric Safety Filters
#### BBB Penetration -- MW > 500 Flags CNS Penetration Risk (Severity: Medium)
#### Cardiac Liability (hERG) -- IC50 < 10 uM Flags QT Prolongation (Severity: Critical)
#### Hepatic Immaturity -- LogP > 5 Flags Hepatotoxicity in Immature Liver (Severity: High)
#### Teratogenicity Screen -- SMARTS Match for Glutarimide/Retinoid Scaffolds (Severity: High)
#### Oral Bioavailability -- TPSA > 140 Flags Pediatric Formulation Challenges (Severity: Medium)
#### GI Absorption -- Rotatable Bonds > 10 Flags Variable Pediatric Absorption (Severity: Low)

### Pipeline Configuration Defaults
#### Max MW: 550
#### Max LogP: 5.0
#### Max Lipinski Violations: 1
#### Top N Candidates: 10
#### Checkpoint Resume: Supported (Per-Stage)

### Demo Target: VCP (Frontotemporal Dementia)
#### Protein: Valosin-Containing Protein (p97), UniProt P55072
#### Mechanism: AAA+ ATPase Inhibition
#### Reference Compound: CB-5083 (Phase I)
#### PDB Structures: 5FTK, 8OOI, 9DIL, 7K56

---

## 13 Therapeutic Areas (201 Genes)

### Neurology / Neurodegeneration (36 Genes)
#### VCP, GRN, C9orf72, MAPT, TBK1, FUS, TARDBP, SOD1
#### APP, PSEN1, PSEN2, LRRK2, SNCA, HTT, PARK7, PINK1
#### PRKN, OPTN, SQSTM1, ATXN1, ATXN2, ATXN3, MS4A6A, TREM2
#### APOE, MS_CD20, IL2RA, CGRP, CALCA, TRPV1, SCN9A, SCN10A
#### GRIN2B, GABRA1, SV2A, GBA (shared with Rare Disease)

### Oncology (27 Genes)
#### BRCA1, BRCA2, EGFR, ALK, BRAF, KRAS, ERBB2, PIK3CA
#### BCR-ABL1, RET, MET, ROS1, NTRK1, IDH1, IDH2, FLT3
#### BTK, BCL2, PDCD1, CD274, CTLA4, TP53, PTEN, MTOR
#### CDK4, CDK6, KIT (shared with Dermatology)

### Metabolic / Endocrine (22 Genes)
#### GLP1R, GIPR, GCGR, INSR, PPARG, PPARA, DPP4, SGLT2
#### KCNJ11, ABCC8, GCK, HNF1A, HNF4A, TSHR, MC4R, LEPR
#### POMC, PCSK1, NPC1L1, ANGPTL3, CETP, HMGCR

### Infectious Disease (21 Genes)
#### HIV1_RT, HIV1_PR, HIV1_IN, CCR5, CD4
#### HCV_NS3, HCV_NS5A, HCV_NS5B, HBV_RT
#### SARS2_MPRO, SARS2_RDRP, SARS2_SPIKE, ACE2, TMPRSS2
#### DHFR, DHPS, GYRA, RPOB, INHA, CYP51A1, FKS1

### Respiratory / Pulmonary (13 Genes)
#### ADRB2, CHRM3, PDE4D, IL5, IL5RA, IL4R, IL13
#### FCER1A, SERPINA1, BMPR2, EDNRA, PDE5A, GUCY1A1

### Rare Disease (12 Genes)
#### CFTR, SMN1, DMD, HBB, F8, F9, GAA, GBA
#### IDUA, PAH, SCN1A, MECP2

### Hematology (12 Genes)
#### SYK, THPO, MPL, EPOR, HIF2A, CALR
#### CXCR4, ADAMTS13, VWF, SERPINC1, F10, F2

### GI / Hepatology (12 Genes)
#### ATP4A, KCNQ1_GI, NOD2, ATG16L1, IL12B, ITGA4
#### ITGB7, S1PR1, PNPLA3, HSD17B13, FXR, THR_BETA

### Ophthalmology (11 Genes)
#### VEGFA, CFH, C3, C5, RHO, RPE65
#### BEST1, ABCA4, USH2A, ROCK1, CA2

### Pharmacogenomics (11 Genes)
#### CYP2C19, CYP2D6, CYP3A4, CYP2C9, VKORC1
#### DPYD, TPMT, UGT1A1, SLCO1B1, NUDT15, HNF1A (shared with Metabolic)

### Cardiovascular (10 Genes)
#### LDLR, PCSK9, APOB, TTR, MYBPC3
#### MYH7, SCN5A, KCNH2, KCNQ1, LMNA

### Immunology / Inflammation (9 Genes)
#### IL6, TNF, JAK1, JAK2, IL17A
#### IL23A, TSLP, TNFRSF1A, IL4R (shared with Respiratory)

### Dermatology (9 Genes)
#### IL31RA, IL36RN, OSM, TYK2, KIT
#### COL7A1, TGM1, SPINK5, FLG

---

## Data Architecture

### ClinVar (4.1M Clinical Variant Records)
### AlphaMissense (71M Missense Pathogenicity Predictions)
### gnomAD v4 (Population Allele Frequencies)
### GIAB HG002 Benchmark Genome (Ashkenazi Son)
### PDB/EMDB Protein Structures
### ClinicalTrials.gov API v2
### CIViC (Clinical Interpretations of Variants in Cancer)
### Orphanet (Rare Disease Ontology)
### HPO (Human Phenotype Ontology)
### CellxGene (Single-Cell Expression Atlases)
### CPIC (Clinical Pharmacogenetics Implementation Consortium)

---

## Infrastructure

### Compute Platform: NVIDIA DGX Spark
#### GPU: GB10 (Grace Blackwell)
#### CPU: 20 ARM Grace Cores
#### Memory: 128 GB Unified LPDDR5x
#### Interconnect: NVLink-C2C
#### CUDA: 12.x

### Docker Compose Orchestration (docker-compose.dgx-spark.yml)

#### Infrastructure Containers
##### milvus-etcd: quay.io/coreos/etcd:v3.5.11 (512 MB)
##### milvus-minio: minio/minio:RELEASE.2024-01-01T16-36-33Z (1 GB)
##### milvus: milvusdb/milvus:v2.4.0 (:19530, :9091, 8 GB)

#### Agent Containers (11)
##### precision-biomarker-agent (:8529 API, :8528 UI, 4 GB)
##### precision-oncology-agent (:8527 API, :8526 UI, 4 GB)
##### cart-intelligence-agent (:8522 API, :8521 UI, 4 GB)
##### imaging-intelligence-agent (:8524 API, :8525 UI, 8 GB)
##### precision-autoimmune-agent (:8532 API, :8531 UI, 4 GB)
##### pharmacogenomics-intelligence-agent (:8107 API, :8507 UI, 4 GB)
##### cardiology-intelligence-agent (:8126 API, :8536 UI, 4 GB)
##### neurology-intelligence-agent (:8536 API, :8535 UI, 4 GB)
##### rare-disease-diagnostic-agent (:8134 API, :8544 UI, 4 GB)
##### clinical-trial-intelligence-agent (:8538 API, :8128 UI, 4 GB)
##### single-cell-intelligence-agent (:8540 API, :8130 UI, 8 GB)

#### Application Containers
##### landing-page (:8080, Flask Hub with Real-Time Health Monitoring)

#### Observability Containers
##### prometheus: prom/prometheus:v2.49.1 (:9099, 1 GB)
##### grafana: grafana/grafana:10.3.1 (:3000, 512 MB)

#### Named Volumes
##### etcd_data, minio_data, milvus_data, prometheus_data, grafana_data

### Additional Pipeline Services (Non-Docker, Process-Based)
#### Genomics Portal (:5000, Python/Flask)
#### RAG/Chat API (:5001, Python/FastAPI)
#### RAG Chat UI (:8501, Streamlit)
#### Drug Discovery UI (:8505, Streamlit)
#### Discovery Portal (:8510, Streamlit)

### Monitoring
#### Prometheus v2.49.1 (:9099) -- Metrics Collection
#### Grafana v10.3.1 (:3000) -- Dashboards (hcls-overview.json)
#### DCGM GPU Exporter (:9400) -- GPU Metrics
#### Node Exporter (:9100) -- Host Metrics

### Nextflow DSL2 Workflows (hcls-orchestrator/)
#### Genomics Pipeline Orchestration
#### MolMIM Container: nvcr.io/nim/nvidia/molmim:1.0.0
#### DiffDock Container: nvcr.io/nim/mit/diffdock:2.2.0

### Security
#### HTTPS via Caddy + TLS Termination
#### API Key Authentication (Per-Agent, ENV-Configurable)
#### XSS/Injection Prevention
#### CORS Origin Whitelisting (Per-Agent)
#### Max Request Size: 10 MB (50 MB for Autoimmune PDF Upload)
#### FDA 21 CFR Part 11 Compliance Design

### Health Monitor (health-monitor.sh, 19 KB)
#### 24 Services Monitored (11 Agents + UIs + Infrastructure)
#### Health Check: HTTP Endpoint or TCP Port Probe
#### HTTP Timeout: 3s, TCP Timeout: 2s
#### Auto-Recovery: max_attempts 30, 2s Between Attempts (60s Window)
#### Cron Schedule: */5 * * * * (Every 5 Minutes)
#### Watchdog Mode: 60s Continuous Loop
#### Log Rotation: 10 MB Max
#### Commands: status [--json], fix, watch, restart <svc|all>, stop, log, install, uninstall

### Startup/Shutdown
#### start-factory.sh (19 KB) -- Master Startup for All Services
#### stop-factory.sh -- Graceful Shutdown by Port

### Total Footprint: ~1.1 TB (Models + Data + Annotation DBs)
### Open Source on DGX Spark
