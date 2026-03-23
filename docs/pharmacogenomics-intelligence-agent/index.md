# Pharmacogenomics Intelligence Agent

> **Part of the [Precision Intelligence Network](../engines/precision-intelligence.md)** — one of 11 specialized agents sharing a common molecular foundation within the HCLS AI Factory.

The **Pharmacogenomics Intelligence Agent** is a domain-specialized retrieval-augmented generation (RAG) system that translates patient genetic data into actionable drug prescribing recommendations. It searches 15 Milvus vector collections containing pharmacogene references, CPIC/DPWG clinical guidelines, drug-gene interactions, HLA hypersensitivity associations, phenoconversion models, validated dosing algorithms, clinical evidence, and population allele frequency data.

## Key Capabilities

| Capability | Detail |
|-----------|--------|
| **Pharmacogenes** | 25 genes including CYP2D6, CYP2C19, CYP2C9, DPYD, TPMT, NUDT15, SLCO1B1 |
| **Drugs Covered** | 100+ across 12 therapeutic categories |
| **Dosing Algorithms** | 9 genotype-guided (IWPC warfarin, tacrolimus, fluoropyrimidine, thiopurine, clopidogrel, simvastatin, SSRI, phenytoin, TCA) |
| **HLA Screening** | 15 HLA-drug hypersensitivity associations |
| **Phenoconversion** | 60 CYP inhibitors/inducers across CYP2D6, CYP2C19, CYP2C9, CYP3A4/5 |
| **Vector Collections** | 15 Milvus collections with BGE-small-en-v1.5 embeddings (384-dim) |
| **Clinical Workflows** | 8 workflows across 10 Streamlit UI tabs |
| **Test Suite** | 1,001 tests across 16 test files |

## Architecture

```
STREAMLIT UI (8507)  ──►  FastAPI REST API (8107)  ──►  PGxRAGEngine
     │                          │                          │
     └── 10 Tabs                └── 16+ Endpoints          └── 15 Collections
         Dashboard                  /health                    pgx_gene_reference
         Drug Check                 /query                     pgx_drug_guidelines
         Medication Review          /v1/pgx/drug-check         pgx_drug_interactions
         Warfarin Dosing            /v1/pgx/warfarin-dose      pgx_hla_hypersensitivity
         Chemo Safety               /v1/pgx/hla-screen         pgx_phenoconversion
         HLA Screening              /v1/pgx/phenoconversion    pgx_dosing_algorithms
         Report Generator           ...                        ...
         Evidence Explorer
         Phenoconversion Modeler
         Population Analytics
```

## Clinical Pipelines

- **Star Allele Calling** -- VCF variants to pharmacogene star allele nomenclature
- **Phenotype Translation** -- Diplotypes to CPIC-standardized metabolizer phenotypes via activity scores
- **Drug-Gene Matching** -- CPIC guideline lookup with alert severity classification
- **Phenoconversion Detection** -- CYP inhibitor/inducer phenotype adjustment
- **HLA Screening** -- Pre-prescription HLA-drug hypersensitivity screening
- **Genotype-Guided Dosing** -- 9 validated dosing algorithms with patient-specific calculations

## Documentation

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } **Project Bible**

    ---

    Complete system reference covering architecture, file inventory, collections, clinical pipelines, and configuration.

    [:octicons-arrow-right-24: Project Bible](project-bible.md)

-   :material-file-document:{ .lg .middle } **White Paper**

    ---

    Technical white paper with clinical validation, performance benchmarks, and comparison with existing PGx solutions.

    [:octicons-arrow-right-24: White Paper](white-paper.md)

-   :material-play-circle:{ .lg .middle } **Demo Guide**

    ---

    Step-by-step walkthrough of all 10 UI tabs with sample inputs, expected outputs, and talking points.

    [:octicons-arrow-right-24: Demo Guide](demo-guide.md)

-   :material-rocket-launch:{ .lg .middle } **Deployment Guide**

    ---

    Docker Compose deployment, manual setup, Milvus tuning, security checklist, and monitoring configuration.

    [:octicons-arrow-right-24: Deployment Guide](deployment-guide.md)

-   :material-school:{ .lg .middle } **Learning Guide -- Foundations**

    ---

    Introduction to pharmacogenomics: CYP450 enzymes, star alleles, diplotypes, activity scores, CPIC guidelines.

    [:octicons-arrow-right-24: Learning Guide -- Foundations](learning-guide-foundations.md)

-   :material-school-outline:{ .lg .middle } **Learning Guide -- Advanced**

    ---

    Advanced topics: phenoconversion modeling, multi-gene interactions, IWPC algorithm, HLA pharmacovigilance, RAG architecture.

    [:octicons-arrow-right-24: Learning Guide -- Advanced](learning-guide-advanced.md)

</div>

## Tech Stack

| Component | Technology |
|-----------|-----------|
| LLM | Claude Sonnet 4.6 (Anthropic) |
| Vector DB | Milvus 2.4 with IVF_FLAT / COSINE |
| Embeddings | BGE-small-en-v1.5 (384-dim) |
| API | FastAPI |
| UI | Streamlit with NVIDIA dark theme |
| Compute | NVIDIA DGX Spark |
| Export | Markdown, JSON, PDF, FHIR R4 |

---

!!! warning "Clinical Decision Support Disclaimer"
    This agent is a clinical decision support research tool. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
