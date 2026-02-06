# Annotation Data

This directory holds the annotation databases used by the RAG pipeline to enrich genomic variants with clinical and pathogenicity information. These files are too large to include in the repository and must be downloaded separately.

## Required Files

### ClinVar VCF

- **File:** `clinvar.vcf.gz` + `clinvar.vcf.gz.tbi`
- **Size:** ~394 MB (compressed)
- **Source:** [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)
- **Description:** Clinical variant interpretations with pathogenicity classifications

```bash
# Download ClinVar VCF (GRCh38)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

### AlphaMissense Predictions

- **File:** `AlphaMissense_hg38.tsv.gz`
- **Size:** ~614 MB (compressed)
- **Source:** [Google DeepMind - AlphaMissense](https://zenodo.org/records/8208688)
- **Description:** Proteome-wide missense variant pathogenicity predictions (71M variants)

```bash
# Download AlphaMissense predictions (requires Zenodo access)
wget https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz
```

## Setup

After downloading, place both files in this directory:

```
rag-chat-pipeline/data/annotations/
├── clinvar.vcf.gz
├── clinvar.vcf.gz.tbi
├── AlphaMissense_hg38.tsv.gz
└── README.md
```

The annotation pipeline will automatically detect and use these files during VCF ingestion.
