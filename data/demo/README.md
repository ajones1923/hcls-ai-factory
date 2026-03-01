# Demo Variant Dataset

## What is this?

`demo_variants.vcf` is a **synthetic** VCF file containing 90 curated variants across 13 therapeutic areas. It allows users to exercise the full Stage 2 (RAG/Chat) pipeline without downloading the 200 GB HG002 reference dataset.

**This file is for demonstration purposes only. Do not use it for clinical, diagnostic, or treatment decisions.**

## How to use with Stage 2

From the `rag-chat-pipeline/` directory, ingest the demo VCF into Milvus:

```bash
# Start Milvus first
docker-compose up -d milvus

# Ingest demo variants (no external annotation files needed -- annotations are embedded in the VCF)
python scripts/ingest_vcf.py \
  --vcf ../data/demo/demo_variants.vcf \
  --no-annotation \
  --drop-existing

# Launch the chat UI
streamlit run app/chat_ui.py
```

Because CSQ annotations and AlphaMissense scores are already embedded in the INFO field, the `--no-annotation` flag skips the ClinVar/AlphaMissense database downloads. The pipeline's fallback parser will read genotype, quality, and rsID directly from the VCF columns, then the RAG system can surface these variants in chat.

## Variant distribution

| Therapeutic Area     | Variants | Key Genes                              |
|----------------------|---------:|----------------------------------------|
| Neurology            |       14 | VCP (5), APP, PSEN1, MAPT, HTT, TARDBP |
| Oncology             |       13 | BRCA1, BRCA2, TP53, EGFR, BRAF, KRAS, PIK3CA |
| Rare Disease         |        8 | CFTR, SMN1, DMD, FMR1                  |
| Metabolic            |        7 | GCK, HNF1A, PCSK9, LDLR               |
| Cardiovascular       |        7 | MYBPC3, SCN5A, KCNQ1                   |
| Pharmacogenomics     |        7 | CYP2D6, CYP2C19, DPYD                  |
| Hematology           |        6 | HBB, JAK2, CALR                        |
| Immunology           |        6 | HLA-B, CTLA4, IL6                      |
| GI / Hepatology      |        5 | APC, MLH1                              |
| Dermatology          |        5 | FLG, MC1R                              |
| Respiratory          |        4 | SERPINA1, DNAI1                         |
| Ophthalmology        |        4 | RHO, RPE65                             |
| Infectious Disease   |        4 | CCR5, IFITM3, TLR7                     |
| **Total**            |   **90** | **40+ genes**                          |

### Variant type breakdown

- SNPs: 67 (74%)
- Indels: 16 (18%)
- Multi-allelic: 7 (8%)

### Impact level breakdown

- HIGH: 32 (36%)
- MODERATE: 43 (48%)
- LOW: 10 (11%)
- MODIFIER: 5 (6%)

### Clinical significance breakdown

- Pathogenic: 46
- Likely pathogenic: 12
- Uncertain significance (VUS): 13
- Likely benign: 9
- Benign: 10

### Notable variants included

- **rs188935092** (VCP p.Arg412His) -- the project demo target for Frontotemporal Dementia
- **rs334** (HBB p.Glu7Val) -- Sickle cell disease
- **rs80357906** (BRCA1 p.Glu1756Ter) -- Hereditary breast/ovarian cancer
- **rs28934578** (TP53 p.Arg248Trp) -- Li-Fraumeni syndrome hotspot
- **rs113488022** (BRAF p.Val600Glu) -- Melanoma driver
- **rs75961395** (CFTR p.Phe508del) -- Cystic Fibrosis
- **rs3918290** (DPYD IVS14+1G>A) -- 5-FU toxicity pharmacogenomic marker
- **rs63750847** (APP p.Ala717Val) -- Early-onset Alzheimer disease

## Format details

- VCF version 4.2, GRCh38/hg38 reference
- VEP-style CSQ annotations in the INFO field
- AlphaMissense pathogenicity scores (AM= tag) on 52 missense variants
- Clinical significance (CLNSIG= tag) on all variants
- Single sample with GT:DP:GQ:AD format fields

## License

Apache 2.0, consistent with the parent HCLS AI Factory project.
