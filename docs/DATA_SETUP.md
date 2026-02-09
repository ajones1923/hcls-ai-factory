---
search:
  exclude: true
---

# Stage 0: Data Acquisition

**The foundation step before any pipeline can run.**

This guide covers downloading and verifying all data required to run the HCLS AI Factory end-to-end (~500 GB). The `setup-data.sh` script automates this process, but this document provides context, troubleshooting, and manual procedures for when things go wrong.

Stage 0 is a one-time setup step that acquires all reference data, sequencing data, and annotation databases needed by the three processing stages. Once complete, the pipeline can run repeatedly without re-downloading.

---

## Quick Start

```bash
# Download everything (run once, ~500 GB total)
./setup-data.sh --all

# Or download stage by stage
./setup-data.sh --stage2    # ClinVar + AlphaMissense (~2 GB, 5 min)
./setup-data.sh --stage1    # HG002 FASTQ + reference (~300 GB, 2-6 hours)
./setup-data.sh --stage3    # PDB structure cache (optional)

# Check what's downloaded
./setup-data.sh --status

# Re-verify all checksums
./setup-data.sh --verify
```

---

## What Gets Downloaded

| Stage | Component | Size | Source | Checksum |
|-------|-----------|------|--------|----------|
| 1 | HG002 FASTQ files (68 files) | ~200 GB | NCBI GIAB FTP | MD5 per file |
| 1 | GRCh38 reference genome + BWA index | ~11 GB | S3 (Parabricks bundle) | -- |
| 1 | Merged FASTQ (R1 + R2) | ~100 GB | Generated locally | -- |
| 2 | ClinVar variant_summary.txt.gz | ~394 MB | NCBI FTP | gzip integrity |
| 2 | ClinVar VCF + tabix index | ~85 MB | NCBI FTP | gzip integrity |
| 2 | AlphaMissense_hg38.tsv.gz | ~614 MB | Google Cloud Storage | gzip integrity |
| 3 | PDB structures | ~3 MB/gene | RCSB PDB (auto-fetched) | -- |
| **Total** | | **~500 GB** | | |

### Time Estimates

| Component | 100 Mbps | 1 Gbps | 10 Gbps |
|-----------|----------|--------|---------|
| FASTQ download | 6 hours | 45 min | 5 min |
| Reference genome | 20 min | 2 min | 15 sec |
| FASTQ merge | 30-60 min (CPU-bound) | 30-60 min | 30-60 min |
| Stage 2 databases | 15 min | 2 min | 15 sec |

---

## Stage-by-Stage Reference

### Stage 1: HG002 FASTQ Files

**What**: Whole-genome sequencing data for GIAB HG002 (Ashkenazi male reference standard). 30x coverage, Illumina 2x250bp paired-end.

**Source**: NCBI Genome in a Bottle (GIAB) FTP
- Index URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina_2x250bps_06012016_updated.HG002`
- Data: 34 R1 files + 34 R2 files across 2 sequencing lanes (L001, L002)

**Checksum**: MD5 checksums provided in the GIAB index TSV file (columns 2 and 4).

**Download tool**: aria2c with parallel connections. Falls back to wget on retry.

**Output paths**:
- Raw chunks: `genomics-pipeline/data/input/giab_hg002/reads/*.fastq.gz`
- Merged: `genomics-pipeline/data/input/giab_hg002/HG002_R1.fastq.gz` and `HG002_R2.fastq.gz`
- Pipeline input: `genomics-pipeline/data/input/HG002_R1.fastq.gz` and `HG002_R2.fastq.gz`

### Stage 1: Reference Genome (GRCh38)

**What**: Human reference genome build 38 with pre-built BWA-MEM2 index files for NVIDIA Parabricks alignment.

**Source**: NVIDIA Parabricks sample bundle on S3
- URL: `https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz`
- Size: ~11 GB compressed

**Contents after extraction**:
- `GRCh38.fa` — Reference FASTA (3.1 GB)
- `GRCh38.fa.bwt` — BWT suffix array (3.0 GB)
- `GRCh38.fa.pac` — Packed sequence (768 MB)
- `GRCh38.fa.sa` — Suffix array (1.5 GB)
- `GRCh38.fa.amb` — Ambiguous bases
- `GRCh38.fa.ann` — Annotations
- `GRCh38.fa.fai` — FASTA index (created via samtools)
- `GRCh38.dict` — Sequence dictionary (created via samtools)

**Output path**: `genomics-pipeline/data/ref/`

### Stage 1: FASTQ Merge

**What**: Combines 68 individual FASTQ chunks (from 2 sequencing lanes) into 2 merged files.

**Process**: Decompresses all chunks with `zcat`, re-compresses with `pigz` (parallel gzip using all CPU cores).

**Time**: 30-60 minutes depending on CPU speed and I/O.

**Requires**: `pigz` for parallel compression (falls back to `gzip` if unavailable, but significantly slower).

### Stage 2: ClinVar

**What**: NCBI Clinical Variants database. Contains clinical significance classifications for 4.1 million human variants.

**Files**:
- `clinvar_variant_summary.txt.gz` (~394 MB) — Tab-delimited format used by the ClinVar annotator
- `clinvar.vcf.gz` (~85 MB) — VCF format for secondary annotation paths
- `clinvar.vcf.gz.tbi` — Tabix index for VCF

**Source**: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`

**Output path**: `rag-chat-pipeline/data/annotations/`

### Stage 2: AlphaMissense

**What**: DeepMind AlphaMissense predictions. AI-predicted pathogenicity scores for 71 million missense variants.

**Source**: `https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz`

**Size**: ~614 MB compressed. Requires ~8-10 GB RAM when loaded into memory during ingestion.

**License**: CC-BY-4.0 (Google DeepMind)

**Output path**: `rag-chat-pipeline/data/annotations/AlphaMissense_hg38.tsv.gz`

### Stage 3: PDB Structures

**What**: Protein structure files from the RCSB Protein Data Bank. Used for molecular docking in Stage 3.

**Behavior**: Auto-fetched from RCSB during pipeline execution. The `CryoEMEvidenceManager` automatically downloads and caches structures for any gene with PDB IDs.

**VCP demo structures** (optionally pre-fetched by setup-data.sh):
- 5FTK — VCP + CB-5083 inhibitor (2.9 MB)
- 7K56 — VCP + cofactor complex (2.9 MB)
- 8OOI — VCP wild-type hexamer (2.9 MB)
- 9DIL — VCP disease mutant (950 KB)

**Output path**: `drug-discovery-pipeline/data/structures/pdb_cache/`

---

## Troubleshooting

### FASTQ downloads fail checksum verification

This is the most common issue. GIAB FASTQ files are large (3-6 GB each) and the NCBI FTP server can be unreliable.

**Symptoms**: `MD5 mismatch` errors during download, or files that download but fail verification.

**Solutions**:

1. **Re-run the script** — It's idempotent. Only failed files are re-downloaded:
   ```bash
   ./setup-data.sh --stage1
   ```

2. **Use fewer parallel connections** — Reduces load on FTP server:
   ```bash
   ./setup-data.sh --stage1 --connections 4
   ```

3. **Wait and retry** — NCBI rate-limits aggressive downloaders. Wait 15-30 minutes:
   ```bash
   sleep 1800 && ./setup-data.sh --stage1
   ```

4. **Try off-peak hours** — NCBI servers are less loaded during US Eastern 2-6 AM.

5. **Check NCBI status** — Verify the server isn't down: https://www.ncbi.nlm.nih.gov/Status/

6. **Verify existing downloads** — See which files are good and which need re-download:
   ```bash
   ./setup-data.sh --verify
   ```

### NCBI FTP server is slow or unreachable

**Symptoms**: Downloads stall, timeouts, "connection refused" errors.

**Solutions**:

1. **Check DNS resolution**:
   ```bash
   nslookup ftp-trace.ncbi.nlm.nih.gov
   ```

2. **Check if HTTPS works** (some networks block FTP):
   ```bash
   curl -I https://ftp-trace.ncbi.nlm.nih.gov/
   ```

3. **Check for firewall/proxy issues**:
   ```bash
   wget --spider https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/
   ```

4. **Use a VPN** — If your institution throttles FTP traffic.

5. **Check if your institution has a GIAB mirror** — Some academic centers maintain local copies.

### Disk space issues

**Symptoms**: "No space left on device" or pre-flight check failure.

**Required space by stage**:
- Stage 1: ~350 GB (200 GB chunks + 11 GB reference + 100 GB merged + workspace)
- Stage 2: ~2 GB
- Stage 3: ~50 MB

**Solutions**:

1. **Download stages independently** — Stage 2 only needs 2 GB:
   ```bash
   ./setup-data.sh --stage2
   ```

2. **Symlink data directories to a larger volume**:
   ```bash
   ln -s /mnt/large-drive/genomics-data genomics-pipeline/data
   ```

3. **Skip the merge** to save 100 GB (merged files can be created later):
   ```bash
   ./setup-data.sh --stage1 --skip-merge
   ```

4. **Free space**:
   ```bash
   sudo apt-get clean
   docker system prune
   ```

### Download interrupted — how to resume

The script is fully **idempotent**. Simply re-run the same command:

```bash
./setup-data.sh --stage1
```

**How it works**:
1. A state file (`.data-setup-state`) tracks which files have been verified.
2. On re-run, verified files are skipped instantly (no re-download, no re-checksum).
3. Only missing or corrupt files are downloaded.
4. The `--verify` flag forces a full re-verification (ignores state file).

**To reset and start fresh**: Delete the state file:
```bash
rm .data-setup-state
```

### DGX Spark (ARM64) specific notes

The HCLS AI Factory is designed to run on the NVIDIA DGX Spark with its ARM64 (aarch64) processor.

- All required tools (`aria2c`, `wget`, `pigz`, `md5sum`) are available in ARM64 Ubuntu repos
- Same commands work on both x86_64 and aarch64
- Docker images used for samtools indexing support multi-architecture
- The DGX Spark's 128 GB unified memory is sufficient for all pipeline stages

**Install tools on DGX Spark**:
```bash
sudo apt-get update
sudo apt-get install -y aria2 pigz
```

### VAST AI OS deployment notes

For VAST R&D porting to VAST AI OS:

- **DataStore**: Point data directories to VAST DataStore paths instead of local filesystem
- **DataEngine**: FASTQ arrival can trigger automatic pipeline execution via event triggers
- **DataBase**: Milvus can be replaced with VAST DataBase (unified SQL + vector)
- **Storage**: VAST's parallel I/O significantly accelerates FASTQ download and merge operations
- Data paths in `genomics-pipeline/config/pipeline.env` and `rag-chat-pipeline/.env` should be updated to VAST mount points

---

## Data Sources and Licensing

| Dataset | Provider | License | URL |
|---------|----------|---------|-----|
| HG002 WGS | NIST (GIAB) | Public Domain | https://www.nist.gov/programs-projects/genome-bottle |
| GRCh38 Reference | NCBI | Public Domain | https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/ |
| ClinVar | NCBI | Public Domain | https://www.ncbi.nlm.nih.gov/clinvar/ |
| AlphaMissense | Google DeepMind | CC-BY-4.0 | https://github.com/google-deepmind/alphamissense |
| Parabricks Bundle | NVIDIA | Free to download | https://docs.nvidia.com/clara/parabricks/ |
| PDB Structures | RCSB | CC0 1.0 | https://www.rcsb.org/ |

---

## Prerequisites

### Required tools

| Tool | Purpose | Install |
|------|---------|---------|
| `aria2c` | Parallel FASTQ downloads | `sudo apt-get install -y aria2` |
| `wget` | Reference genome + Stage 2 downloads | Usually pre-installed |
| `pigz` | Parallel FASTQ compression during merge | `sudo apt-get install -y pigz` |
| `md5sum` | Checksum verification | Pre-installed on Linux |
| `Docker` | Samtools indexing, Milvus, BioNeMo NIMs | https://docs.docker.com/get-docker/ |

### API keys

| Key | Required For | Get It |
|-----|-------------|--------|
| NGC API Key | Parabricks container, BioNeMo NIMs | https://ngc.nvidia.com/ |
| Anthropic API Key | Claude AI in RAG pipeline | https://console.anthropic.com/ |
| HuggingFace Token | Local LLM (optional, if using Llama) | https://huggingface.co/settings/tokens |
