# Genomics Pipeline: Detailed Technical Overview

## Executive Summary

The Genomics Pipeline is a GPU-accelerated bioinformatics workflow designed to transform raw DNA sequencing data (FASTQ format) into clinically actionable variant calls (VCF format). Built on NVIDIA Parabricks, it leverages GPU parallelization to achieve 10-50x speedups compared to traditional CPU-based pipelines, completing whole-genome analysis in under 2 hours rather than the typical 24-48 hours.

---

## Pipeline Architecture Overview

### High-Level Data Flow

The pipeline follows a linear workflow through four major phases: **Setup**, **Data Acquisition**, **Processing**, and **Output Generation**. Each phase builds upon the previous, with checkpoint validation ensuring data integrity before proceeding to computationally expensive steps.

**Phase 1 - Environment Setup** establishes the computational foundation by validating hardware prerequisites (NVIDIA GPU, Docker runtime, sufficient storage), authenticating with NVIDIA's NGC container registry, and pulling the Parabricks container image. This phase is critical because GPU-accelerated bioinformatics requires specific driver versions, container runtime configurations, and authentication tokens that must be verified before any data processing begins.

**Phase 2 - Data Acquisition** handles the retrieval and preparation of two essential datasets: the raw sequencing reads (FASTQ files) and the reference genome (GRCh38). The FASTQ files represent the "query" - billions of short DNA fragments read from a biological sample. The reference genome provides the "template" - the canonical human genome sequence against which the query fragments will be aligned. This phase employs parallel downloading via aria2c with automatic resume capability, crucial for the ~200GB dataset transfers involved.

**Phase 3 - Core Processing** executes the computationally intensive bioinformatics algorithms. This is where GPU acceleration provides the most significant performance gains. The processing sub-pipeline includes read alignment (mapping each DNA fragment to its corresponding location in the reference genome), coordinate sorting (organizing aligned reads by genomic position), duplicate marking (identifying PCR artifacts), and variant calling (detecting positions where the sample differs from the reference).

**Phase 4 - Output Generation** produces the final deliverables: a BAM file containing all aligned reads with quality metrics, and a VCF file containing the detected genetic variants. These outputs are indexed for efficient random access and validated against expected quality thresholds.

---

## Detailed Component Analysis

### 1. Prerequisites and Environment Validation

The pipeline begins with comprehensive system validation through the `00-setup-check.sh` script. This component performs five critical checks:

**Docker Validation** confirms that Docker Engine is installed and the daemon is running. The pipeline is entirely containerized, meaning all bioinformatics tools execute within Docker containers rather than requiring local installation. This ensures reproducibility across different host systems and simplifies dependency management.

**NVIDIA Container Runtime Validation** tests that the nvidia-container-toolkit is properly configured by attempting to run a test container with GPU access. This runtime allows Docker containers to access host GPU resources, which is essential for Parabricks acceleration. The test specifically verifies that `--gpus all` flag successfully exposes GPU devices inside containers.

**GPU Detection** runs nvidia-smi to confirm GPU availability, driver version compatibility, and available memory. Parabricks requires specific minimum driver versions and GPU architectures (Volta or newer). The script reports GPU model, driver version, and memory capacity to help users understand their performance expectations.

**Disk Space Assessment** calculates available storage in the project directory. Whole-genome analysis requires approximately 500GB: ~200GB for input FASTQ files, ~100GB for the reference genome and indexes, ~100GB for output BAM files, and ~100GB temporary working space. Insufficient storage is a common failure mode that this check prevents.

**Network Connectivity** (implicit) - While not explicitly tested, subsequent steps require network access to NGC and data repositories. The setup check's Docker pull operations implicitly verify network availability.

### 2. NGC Authentication and Container Management

NVIDIA GPU Cloud (NGC) serves as the distribution platform for Parabricks containers. The `01-ngc-login.sh` script manages authentication using Docker's credential system:

**Authentication Flow**: Users provide their NGC API key (obtained from ngc.nvidia.com), which the script passes to `docker login nvcr.io`. The special username `$oauthtoken` indicates OAuth-based authentication rather than traditional username/password. Successful authentication stores credentials in `~/.docker/config.json` for subsequent container pulls.

**Container Pull**: Upon successful authentication, Docker can pull the Parabricks image (`nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1`). This ~15GB image contains all required bioinformatics tools pre-configured for GPU execution: BWA-MEM2 for alignment, samtools for BAM manipulation, DeepVariant for variant calling, and supporting utilities.

**License Considerations**: Parabricks is commercial software with usage tracking. The NGC authentication enables license validation and usage metering. Academic and evaluation licenses are available through NVIDIA's developer programs.

### 3. Data Acquisition: FASTQ Download and Preparation

The `02-download-data.sh` script orchestrates retrieval of the Genome in a Bottle (GIAB) HG002 reference sample, a gold-standard benchmark dataset maintained by NIST for validating genomics pipelines.

**Dataset Characteristics**: HG002 (NA24385) is an Ashkenazi Jewish male whose genome has been extensively characterized through multiple sequencing technologies. The Illumina 2x250bp dataset provides high-coverage whole-genome sequencing data - approximately 800 million paired-end reads totaling ~200GB compressed.

**Download Architecture**: The script employs aria2c, a high-performance download utility supporting parallel connections and automatic resume. Configuration uses 8 connections per file (`-x 8`), 8 parallel file downloads (`-j 4`), and segment-based downloading (`-s 8`) to maximize throughput. Downloaded files are verified via checksums embedded in the index file.

**Lane Merging**: Illumina sequencers produce data across multiple "lanes" - parallel sequencing channels that increase throughput. The raw download contains separate files for each lane (e.g., `L001`, `L002`, etc.). The script concatenates these into unified R1 (forward reads) and R2 (reverse reads) files required by downstream tools. This merging preserves read pairing - the nth read in R1.fastq.gz corresponds to the nth read in R2.fastq.gz.

**Validation**: Post-merge validation counts reads in both files to confirm matching counts. Mismatched counts indicate corruption or incomplete merging, which would cause alignment failures. The script reports total read counts (typically ~400 million paired reads for 30x coverage).

### 4. Reference Genome Setup

The `03-setup-reference.sh` script prepares the GRCh38 human reference genome and associated index files required for read alignment.

**Reference Genome**: GRCh38 (Genome Reference Consortium Human Build 38) is the current standard human reference assembly. The FASTA file contains ~3.1 billion nucleotides across 24 chromosomes plus alternate contigs and decoy sequences. This ~3GB file serves as the coordinate system for all downstream analysis.

**Index Generation**: Read alignment algorithms require pre-computed index structures for efficient searching:

- **BWA Index** (`.bwt`, `.pac`, `.sa`, `.amb`, `.ann`): The Burrows-Wheeler Transform index enables BWA-MEM2 to rapidly locate candidate alignment positions. Index construction is computationally expensive (~1 hour on CPU) but the pipeline uses pre-built indexes from the Parabricks sample bundle.

- **FASTA Index** (`.fai`): A simple tabular index mapping chromosome names to file offsets, enabling random access to any genomic region without scanning the entire file.

- **Sequence Dictionary** (`.dict`): SAM/BAM header information listing all reference sequences with their lengths and MD5 checksums. Required for GATK-compatible tools.

**Source**: The pipeline downloads NVIDIA's Parabricks sample bundle (~11GB), which contains the reference genome with pre-built indexes optimized for GPU alignment. This avoids the hours-long index construction step.

### 5. Core Processing: Alignment Pipeline (fq2bam)

The `fq2bam` command is Parabricks' flagship tool, combining multiple alignment steps into a single GPU-accelerated operation:

**Read Alignment (BWA-MEM2)**: Each read pair is aligned to the reference genome using the BWA-MEM2 algorithm. This involves:
1. Seeding: Finding exact matches between read subsequences and the reference index
2. Chaining: Connecting nearby seeds into candidate alignment regions
3. Smith-Waterman: Computing optimal local alignments with gap penalties
4. Scoring: Assigning mapping quality scores based on alignment uniqueness

GPU acceleration parallelizes these operations across thousands of CUDA cores. A single GPU can process ~100 million reads per hour, compared to ~5 million on a high-end CPU.

**Coordinate Sorting**: Aligned reads are sorted by genomic position (chromosome, then coordinate). This ordering is essential for downstream analysis - variant callers examine reads covering each position sequentially. GPU-accelerated sorting uses parallel merge-sort algorithms optimized for GPU memory hierarchies.

**Duplicate Marking**: PCR amplification during library preparation can create identical copies of the same original DNA fragment. These duplicates inflate coverage metrics and can bias variant calling. The duplicate marking algorithm identifies reads with identical alignment coordinates and flags all but one as duplicates. GPU implementation uses hash-based duplicate detection parallelized across read groups.

**Output**: A single coordinate-sorted, duplicate-marked BAM file containing all aligned reads with quality metrics, alignment scores, and duplicate flags.

### 6. Core Processing: Variant Calling (DeepVariant)

DeepVariant represents a paradigm shift in variant calling, using deep learning rather than statistical models:

**Pileup Image Generation**: For each candidate variant position, DeepVariant constructs a "pileup image" - a visual representation encoding:
- Reference sequence context (surrounding nucleotides)
- Aligned read bases at the position
- Base quality scores (color intensity)
- Mapping quality scores
- Strand orientation
- Read pairing information

**Convolutional Neural Network**: The pileup images are processed by an Inception-v3 CNN trained on millions of validated variant calls. The network outputs probabilities for three genotype states: homozygous reference (0/0), heterozygous (0/1), or homozygous alternate (1/1).

**GPU Acceleration**: CNN inference is highly parallelizable on GPUs. DeepVariant processes thousands of candidate positions simultaneously, with GPU tensor cores accelerating the matrix operations underlying neural network computation.

**Variant Types**: DeepVariant detects:
- SNVs (Single Nucleotide Variants): Single base substitutions
- Small Indels: Insertions and deletions up to ~50bp
- The output VCF includes quality scores (QUAL), genotype likelihoods (PL), and read depth (DP) for each variant

### 7. Web Portal Architecture

The web portal provides a graphical interface for pipeline management, built on a Flask backend with Server-Sent Events for real-time updates:

**Backend (Flask)**: The `server.py` application exposes REST endpoints for:
- `/api/status`: System and pipeline state monitoring
- `/api/run/<step>`: Initiating pipeline steps
- `/api/stop`: Process termination
- `/api/config`: Configuration management
- `/api/stream`: Real-time log streaming via SSE

**Process Management**: Pipeline scripts execute as subprocesses with stdout/stderr captured and streamed to connected clients. The backend maintains global state tracking current step, status, and accumulated logs.

**System Monitoring**: Real-time metrics collection via:
- `psutil`: CPU utilization, memory usage, disk space
- `pynvml`: GPU utilization, memory, temperature, power draw
- Docker API: Container status

**Frontend**: Bootstrap 5 UI with JavaScript handling SSE connections, dynamic status updates, and console output rendering. The interface provides click-to-run workflow steps, configuration editing, and log browsing.

---

## Data Flow Summary

```
Input Data                    Processing                      Output
─────────────────────────────────────────────────────────────────────────
FASTQ R1 (100GB) ─┐
                  ├──► BWA-MEM2 ──► Sort ──► MarkDup ──► BAM (100GB)
FASTQ R2 (100GB) ─┘         │                              │
                            │                              │
GRCh38.fa (3GB) ────────────┘                              │
                                                           │
                                                           ▼
                            DeepVariant ◄──────────────────┘
                                │
                                ▼
                            VCF (1-2GB)
```

---

## Performance Characteristics

| Metric | Chr20 Test | Full Genome |
|--------|------------|-------------|
| Input Size | ~10GB subset | ~200GB |
| Alignment Time | 2-10 min | 20-75 min |
| Variant Calling | 1-6 min | 10-35 min |
| Total Runtime | 5-20 min | 30-110 min |
| GPU Memory | 8-12 GB | 12-16 GB |
| Peak Throughput | 10 GB/hr | 45-60 GB/hr |

---

## Quality Assurance

The pipeline incorporates multiple validation checkpoints:

1. **Input Validation**: Read count verification after FASTQ merging
2. **Alignment QC**: samtools flagstat reports mapping rate, duplicate rate, properly paired percentage
3. **File Verification**: Size and existence checks before proceeding to expensive steps
4. **Retry Logic**: Automatic retries with exponential backoff for transient failures
5. **Output Validation**: VCF variant counts and index verification

---

## Integration Points

**Upstream**: Accepts standard Illumina paired-end FASTQ files from any sequencing provider

**Downstream**: Produces standard VCF 4.2 compatible with:
- Pharmacogenomics analysis (PharmCAT, Stargazer)
- Clinical annotation (ClinVar, ACMG classification)
- Population analysis (gnomAD comparison)
- Research workflows (GWAS, burden testing)
