---
search:
  exclude: true
---

# HCLS AI Factory - Precision Medicine to Drug Discovery

## Patient Data Source
### DNA Sample
- Blood/Saliva/Tissue collection
- Patient consent and tracking

### Illumina Sequencer
- Short-Read Sequencing
- 2x250bp Paired-End reads
- High-throughput processing

### FASTQ Output Files
#### FASTQ R1 (Forward Reads)
- ~100GB file size
- ~400M Read Pairs
- Quality encoded

#### FASTQ R2 (Reverse Reads)
- ~100GB file size
- Paired with R1
- Phred Quality Scores

## HCLS Landing Page (Port 8080)
### Flask Server Core
#### server.py
- Flask + CORS enabled
- Auto-Start Dependencies
- Service orchestration

#### Health Monitor
- 10 Services tracked
- 30s Refresh Interval
- Real-time status

#### REST API
- /api/check-services
- /api/check-service/<id>
- JSON responses

### Web Interface
#### Pipeline Cards
- 5 Pipeline Interfaces
- Click-to-Launch navigation
- Visual status indicators

#### Monitor Cards
- 4 Infrastructure Services
- Real-Time Status display
- Health indicators

#### Platform Statistics
- 36,000 Lines of Code
- Animated Counters
- Key metrics display

### Auto-Start Services
#### Genomics Portal
- Port 5000
- Automatic startup
- Health monitoring

#### RAG/Chat API
- Port 5001
- Automatic startup
- Dependency check

## Stage 1: Genomics Pipeline (FASTQ to VCF, 120-240 Minutes)
### Input Layer
#### Input FASTQ Files
- HG002_R1.fastq.gz
- HG002_R2.fastq.gz
- ~200GB Total size

#### Reference Genome
- GRCh38.fa (3.1GB)
- BWA Index Files
- FASTA Index (.fai)

### Web Portal (Port 5000)
#### Flask Server
- app/server.py
- Real-Time Monitoring
- SSE streaming support

#### Web Interface
- templates/index.html
- Click-to-Run Steps
- Progress visualization

#### Frontend Logic
- static/js/app.js
- SSE Streaming
- Real-time updates

#### Styling
- static/css/style.css
- Bootstrap 5
- Responsive design

### Docker Container
#### clara-parabricks:4.6.0-1
- NVIDIA Container Runtime
- GPU Passthrough
- Optimized for genomics

### Step 1: Alignment (20-45 min)
#### pbrun fq2bam
- BWA-MEM2 (GPU accelerated)
- Coordinate Sorting
- PCR Duplicate Marking
- 10-50x faster than CPU

### Step 2: Indexing (2-5 min)
#### samtools index
- BAM Index (.bai) creation
- samtools flagstat
- Alignment QC Statistics

### Step 3: Variant Calling (10-35 min)
#### pbrun deepvariant
- Google DeepVariant
- CNN-Based Caller
- GPU Accelerated
- >99% Accuracy

### Output Layer
#### BAM File
- HG002.genome.bam
- ~100GB Aligned Reads
- Sorted + Deduplicated

#### VCF File
- HG002.genome.vcf.gz
- ~11.7M Variants
- SNPs + Indels

## Stage 2: RAG/Chat Pipeline (VCF to Target Hypothesis, Interactive)
### Annotation Layer
#### ClinVar Database
- 4.1M Clinical Variants
- Pathogenicity Status
- Disease Associations
- Clinical significance

#### AlphaMissense
- 71M AI Predictions
- Pathogenicity Scores
- 0.0-1.0 Range
- DeepMind model

#### VEP (Variant Effect Predictor)
- Functional Consequences
- Gene Impact assessment
- Protein Changes

#### Annotated Variants
- 35,616 ClinVar Matches
- 6,831 AlphaMissense hits
- Combined Evidence

### Embedding Layer
#### BGE-small-en-v1.5
- BAAI Embedding Model
- 384 Dimensions
- Semantic Encoding
- Fast inference

#### Variant Embeddings
- 3.5M Vectors
- Semantic Representation
- Query-Ready index

### Vector Database (Port 19530)
#### Milvus
- Vector Similarity Search
- Millisecond Latency
- Hybrid Filtering
- Scalable architecture

#### genomic_variants Collection
- Collection Schema
- Metadata + Vectors
- IVF_FLAT Index

### Knowledge Layer
#### Clinker Knowledge Base
- 201 Target Genes
- 150+ Diseases
- 13 Therapeutic Areas
- Expert curated

#### Gene Coverage
- Oncology: 45 genes
- Neurology: 38 genes
- Rare Disease: 52 genes
- Cardiovascular: 28 genes

#### Druggability Assessment
- 171 Druggable Targets
- 85% Druggability Rate
- Known Inhibitors mapped

### API Portal (Port 5001)
#### Flask API Server
- portal/app/server.py
- REST Endpoints
- JSON responses

#### /api/search
- Semantic Search
- Metadata Filtering
- Fast retrieval

#### /api/query
- Natural Language input
- RAG Pipeline execution
- Evidence synthesis

### Chat Interface (Port 8501)
#### Streamlit UI
- app/chat_ui.py
- Interactive Chat
- Session management

#### Natural Language Input
- Example: "What pathogenic variants are associated with VCP?"
- Free-form queries
- Context-aware

#### AI Response
- Grounded in Evidence
- Citations Included
- Structured output

### LLM Layer
#### Claude (Anthropic)
- claude-sonnet-4-20250514
- RAG Grounding
- Evidence Synthesis
- Advanced reasoning

#### System Prompt
- Genomics Expert Role
- Citation Requirements
- Structured Output format

### Output Layer
#### Target Hypothesis
- Example: "VCP is a druggable target for FTD"
- Evidence-backed
- Actionable insights

#### Supporting Evidence
- Variant Details
- Clinical Significance
- Literature References

## Stage 3: Drug Discovery Pipeline (Target to Drug Candidates, Minutes)
### Phase 5: Structure Evidence
#### RCSB PDB API
- Protein Data Bank
- Real-Time Fetch
- Structure retrieval

#### Cryo-EM Structures
- 8OOI: WT Hexamer 2.9A
- 9DIL: Mutant 3.2A
- 7K56: Complex 2.5A
- 5FTK: +CB-5083 2.3A

#### Binding Site Analysis
- D2 ATPase Domain
- ATP-Competitive Pocket
- Key Residues Mapped

#### Structure Cache
- PDB Files (.pdb)
- Structure Images (.jpeg)
- Local Storage

### Seed Molecule
#### CB-5083
- Known VCP Inhibitor
- Phase I Clinical
- Reference Structure

#### SMILES Encoding
- Molecular String
- Generation Seed
- Chemical representation

### BioNeMo NIM Microservices
#### MolMIM
- Molecule Generation
- Masked Language Model
- Novel Analogs creation
- NVIDIA BioNeMo

#### DiffDock
- Molecular Docking
- Diffusion-Based
- Binding Pose Prediction
- GPU accelerated

#### 3D Conformers
- RDKit Generation
- Energy Minimization
- SDF Output format

### Drug-Likeness Scoring
#### Lipinski's Rule of 5
- MW <= 500 Da
- LogP <= 5
- HBD <= 5
- HBA <= 10

#### QED Score
- Quantitative Estimate of Drug-likeness
- 0.0-1.0 Scale
- Multi-parameter optimization

#### ADMET Properties
- Absorption
- Distribution
- Metabolism
- Excretion
- Toxicity

#### Candidate Ranking
- Binding Affinity
- Drug-likeness score
- Synthetic Feasibility

### Main UI (Port 8505)
#### Streamlit Interface
- app/discovery_ui.py
- Structure Visualization
- Interactive controls

#### 3Dmol.js Viewer
- Interactive 3D rendering
- Binding Site View
- Protein-ligand display

#### Generation Controls
- Similarity Threshold
- Number of Candidates
- Scoring Weights

### Discovery Portal (Port 8510)
#### Management Dashboard
- portal/app.py
- Pipeline Orchestration
- Overview display

#### Target Management
- Active Targets List
- Progress Tracking
- Status monitoring

#### Generation History
- Previous Runs
- Result Comparison
- Export options

### Report Generation
#### ReportLab
- PDF Generation
- Professional Layout
- Custom branding

#### Drug Candidate Report
- VCP_Drug_Candidate_Report.pdf
- Executive Summary
- Ranked Candidates
- Structure Images
- Scoring Details

## Monitoring and Infrastructure
### Grafana (Port 3000)
#### Dashboard
- nvidia-dgx-spark
- GPU Monitoring
- Real-time metrics

#### Dashboard Panels
- GPU Utilization
- CPU Utilization
- GPU Temperature
- GPU Power Usage
- Memory Bandwidth
- NVMe Throughput

### Prometheus (Port 9099)
#### Metrics Server
- Metrics Collection
- Time Series DB
- 15s Scrape Interval

#### Scrape Targets
- Node Exporter
- DCGM Exporter
- Application Metrics

### Metric Exporters
#### Node Exporter (Port 9100)
- System Metrics
- CPU monitoring
- RAM monitoring
- Disk monitoring
- Network monitoring

#### DCGM Exporter (Port 9400)
- GPU Metrics
- NVIDIA Data Center GPU Manager
- Temperature and power

## NVIDIA DGX Spark Infrastructure
### GPU Resources
#### NVIDIA GB10 GPU
- 128GB HBM3 Memory
- Blackwell Architecture
- CUDA 12.x support

#### Tensor Cores
- DeepVariant Inference
- AI Acceleration
- Matrix operations

#### CUDA Cores
- BWA-MEM2 Alignment
- General Compute
- Parallel processing

### System Resources
#### CPU
- ARM64 Cores
- ARM Architecture
- Parallel Processing

#### System RAM
- 128GB unified LPDDR5x
- BAM Processing
- Large Dataset handling

#### NVMe Storage
- 2TB+ Capacity
- High IOPS
- Fast I/O throughput

### Container Runtime
#### Docker
- Version 24.0+
- Container Orchestration
- Image management

#### NVIDIA Container Runtime
- GPU Passthrough
- CUDA Support
- Device mapping

## HCLS AI Factory Orchestration
### Startup Scripts
#### start-services.sh
- Master Startup script
- All Services launch
- Dependency ordering

#### --status Flag
- Service Health Check
- Port verification
- Status display

#### --stop Flag
- Graceful Shutdown
- Resource cleanup
- Process termination

### Documentation
#### README.md
- 650+ Lines
- Complete Guide
- Quick start instructions

#### Product Documentation
- 3,200+ Lines
- Technical Reference
- API documentation

#### Executive Summary
- Business Overview
- Key Metrics
- Value proposition

### Configuration
#### Environment Variables
- ANTHROPIC_API_KEY
- NGC_API_KEY
- Service Ports

#### Port Assignments
- 8080: Landing Page
- 5000: Genomics Portal
- 5001: RAG API
- 8501: Chat Interface
- 8505: Drug Discovery UI
- 8510: Discovery Portal
- 19530: Milvus
- 3000: Grafana
- 9099: Prometheus

## Data Flow Summary
### Pipeline Connections
- Patient DNA -> Genomics Pipeline
- FASTQ files -> Parabricks processing
- VCF output -> RAG/Chat Pipeline
- Target Hypothesis -> Drug Discovery
- Drug Candidates -> PDF Report

### Key Metrics
- Lines of Code: 36,000
- Target Genes: 201
- Variant Embeddings: 3.5M
- ClinVar Variants: 4.1M
- AlphaMissense Predictions: 71M
- End-to-End Time: ~5 hours
- Therapeutic Areas: 13
- Druggable Targets: 171 (85%)
