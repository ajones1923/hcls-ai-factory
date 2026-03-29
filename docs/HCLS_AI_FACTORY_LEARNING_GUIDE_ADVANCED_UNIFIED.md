# HCLS AI Factory — Learning Guide: Unified Advanced

Author: Adam Jones | License: Apache 2.0 | Date: March 2026

---

## Part 1: Advanced Genomic Analysis

### Custom Variant Filtering Strategies

The default genomics pipeline (Parabricks 4.6 + DeepVariant) produces a comprehensive VCF containing every called variant. Real clinical workflows demand surgical filtering. The strategies below go beyond `bcftools view -f PASS`.

**Tiered Clinical Filtering**

A three-tier approach mirrors clinical reporting standards:

| Tier | Criteria | Typical Count |
|------|----------|---------------|
| Tier 1 — Actionable | QUAL >= 100, DP >= 30, GQ >= 40, ClinVar Pathogenic/Likely Pathogenic | ~50-200 |
| Tier 2 — Potentially Significant | QUAL >= 50, DP >= 20, GQ >= 30, ClinVar VUS + AlphaMissense >= 0.8 | ~500-2,000 |
| Tier 3 — Observed | PASS filter, DP >= 10, population AF < 0.01 | ~20,000-50,000 |

```bash
# Tier 1: Actionable variants only
bcftools view -f PASS input.vcf.gz \
  | bcftools filter -i 'QUAL>=100 && INFO/DP>=30 && FORMAT/GQ>=40' \
  | bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG \
  | bcftools filter -i 'INFO/CLNSIG~"Pathogenic"'
```

**Gene Panel Filtering**

When the clinical question targets a specific condition, restrict analysis to curated gene panels:

```bash
# Extract regions from a BED file (e.g., hereditary cancer panel — 84 genes)
bcftools view -R cancer_panel.bed -f PASS sample.vcf.gz > panel_variants.vcf
```

The HCLS AI Factory ships with pre-built BED files for 13 therapeutic areas aligned to the intelligence agents.

**Population Frequency Filtering**

Filter against gnomAD allele frequencies to surface rare variants:

```bash
# Retain variants with AF < 0.001 or absent from gnomAD
bcftools filter -i 'INFO/gnomAD_AF<0.001 || INFO/gnomAD_AF="."' input.vcf.gz
```

### Multi-Sample and Trio Analysis Patterns

**Trio Analysis (Proband + Parents)**

Parabricks 4.6 supports joint calling. For de novo variant detection:

```bash
# Joint calling with Parabricks DeepVariant
pbrun deepvariant \
  --ref reference.fasta \
  --in-bam proband.bam,mother.bam,father.bam \
  --out-variants trio_joint.vcf

# De novo filtering: variant present in proband, absent in both parents
bcftools view -s proband trio_joint.vcf \
  | bcftools filter -i 'GT="0/1" || GT="1/1"' \
  | bcftools isec -C -w1 - mother_only.vcf father_only.vcf
```

**Multi-Sample Cohort Analysis**

For cohort studies (e.g., tumor-normal pairs or family groups), use merged VCFs with sample-specific genotype fields. The RAG pipeline can ingest multi-sample VCFs and maintain per-sample lineage in the `genomic_evidence` Milvus collection.

### Structural Variant Detection

The default pipeline focuses on SNVs and small indels. Structural variants (SVs) require dedicated tools:

| Tool | SV Types | Integration Point |
|------|----------|-------------------|
| Parabricks pbsv | DEL, INS, DUP, INV, BND | Direct GPU-accelerated calling |
| Manta | DEL, INS, DUP, INV, BND | Post-alignment, CPU |
| CNVkit | Copy number variants | BAM coverage analysis |

SVs are annotated against ClinVar SV records and can be ingested into the RAG pipeline using the `genomic_evidence` collection with `variant_type` metadata filtering.

### ClinVar + AlphaMissense Annotation Deep Dive

The platform cross-references two major annotation databases:

**ClinVar (4.1M records)**

- Variant-level clinical significance: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign
- Review status (0-4 stars) indicates evidence strength
- Condition associations with OMIM and MedGen identifiers
- Submission history tracks reclassification over time

**AlphaMissense (71M predictions)**

- Machine learning pathogenicity scores (0.0 to 1.0) for all possible single amino acid substitutions
- Thresholds: >= 0.564 likely pathogenic, <= 0.340 likely benign
- Covers missense variants that ClinVar has not yet classified
- Particularly valuable for VUS reclassification support

**Cross-Reference Strategy**

The RAG pipeline performs a layered annotation:

1. Check ClinVar first — if Pathogenic/Likely Pathogenic with >= 2 stars, use as primary evidence
2. For VUS or missing ClinVar entries, check AlphaMissense score
3. Concordance between ClinVar and AlphaMissense strengthens confidence
4. Discordance flags the variant for manual review

### VCF Quality Metrics Interpretation

| Metric | Field | Good | Marginal | Poor |
|--------|-------|------|----------|------|
| Variant Quality | QUAL | >= 100 | 30-100 | < 30 |
| Read Depth | DP | >= 30x | 15-30x | < 15x |
| Genotype Quality | GQ | >= 40 | 20-40 | < 20 |
| Allele Depth | AD | Ref/Alt balanced for het | Skewed > 4:1 | Single-strand only |
| Mapping Quality | MQ | >= 50 | 30-50 | < 30 |
| Strand Bias | FS | < 60 | 60-200 | > 200 |

**Interpreting AD (Allele Depth) for heterozygous calls:**

A heterozygous SNV should show roughly 50/50 reference/alternate reads. A call showing `AD=95,5` at `DP=100` is suspicious — possibly a sequencing artifact or somatic mosaic variant. The expected range for germline heterozygous calls is 30-70% alternate allele fraction.

### Parabricks 4.6 Performance Tuning on DGX Spark

The NVIDIA DGX Spark (Grace Blackwell GB10, 128GB unified memory) runs the full genomics pipeline. Key tuning parameters:

| Parameter | Default | Tuned | Impact |
|-----------|---------|-------|--------|
| `--num-gpus` | 1 | 1 (GB10 is single-GPU) | N/A |
| `--bwa-options` | default | `-K 100000000 -Y` | Improves reproducibility |
| `--run-partition-size` | 25M | 50M | Fewer partitions, less overhead |
| `--memory-limit` | auto | 100G | Reserve 28GB for system + Milvus |

**Benchmark on DGX Spark:**

| Stage | Time (30x WGS) |
|-------|-----------------|
| BWA-MEM2 alignment | ~45 min |
| DeepVariant calling | ~60 min |
| Annotation + VCF post-processing | ~15 min |
| **Total** | **~120 min** |

This compares to 24-48 hours on a 32-core CPU-only server.

---

## Part 2: Advanced RAG Architecture

### Milvus 2.4 Collection Design Patterns

Every intelligence agent in the HCLS AI Factory follows a consistent Milvus collection design. Understanding these patterns is essential for extending or customizing the platform.

**Index Configuration**

All collections use the same core index parameters:

```python
index_params = {
    "index_type": "IVF_FLAT",
    "metric_type": "COSINE",
    "params": {"nlist": 128}
}

search_params = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16}
}
```

- **IVF_FLAT**: Inverted file index with flat (exact) distance computation within clusters. Chosen over HNSW for its lower memory footprint and deterministic results.
- **COSINE**: Cosine similarity metric. Normalized embeddings make this equivalent to inner product but more interpretable (0.0 to 1.0 scale).
- **nlist=128**: Number of Voronoi cells. With collections ranging from hundreds to millions of vectors, 128 provides a good balance between recall and speed.
- **nprobe=16**: Number of cells to search at query time (12.5% of nlist). Increasing to 32 improves recall by ~2% but doubles search latency.

**Schema Pattern**

Every collection follows this field layout:

| Field | Type | Purpose |
|-------|------|---------|
| `id` | INT64 (auto-increment) | Primary key |
| `embedding` | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| `text` | VARCHAR(65535) | Source text chunk |
| `source` | VARCHAR(512) | Provenance (PubMed ID, URL, file path) |
| `metadata` | VARCHAR(65535) | JSON-encoded domain metadata |
| `created_at` | VARCHAR(64) | ISO 8601 timestamp |
| Domain-specific fields | Various | Filterable attributes (e.g., `cancer_type`, `target_antigen`, `gene_symbol`) |

**Shared `genomic_evidence` Collection**

All 11 agents mount the `genomic_evidence` collection as read-only. This collection holds 3.56M vectors derived from the Stage 2 RAG pipeline (ClinVar + AlphaMissense annotations). Agents query it to enrich domain-specific analysis with genomic context without duplicating the data.

**Collection Sizing Across Agents**

| Agent | Owned Collections | Approx. Vectors | Shared |
|-------|-------------------|-----------------|--------|
| Precision Oncology | 10 | ~10,000+ | +genomic_evidence |
| CAR-T Intelligence | 10 | 6,266 | +genomic_evidence |
| Precision Biomarker | 10 | ~320 | +genomic_evidence |
| Cardiology | 12 | ~5,000+ | +genomic_evidence |
| Neurology | 12 | ~5,000+ | +genomic_evidence |
| Autoimmune | 13 | ~5,000+ | +genomic_evidence |
| Rare Disease | 13 | ~5,000+ | +genomic_evidence |
| Pharmacogenomics | 14 | ~5,000+ | +genomic_evidence |
| Imaging | 10 | 2,814+ | +genomic_evidence |
| Single-Cell | 12 | ~5,000+ | +genomic_evidence |
| Clinical Trial | 13 | ~5,000+ | +genomic_evidence |

### BGE-small-en-v1.5 Embedding Model

**Why BGE-small-en-v1.5?**

| Property | Value |
|----------|-------|
| Dimensions | 384 |
| Model size | 33M parameters (~130 MB) |
| Max sequence length | 512 tokens |
| Embedding speed (DGX Spark) | ~1,200 texts/sec |
| MTEB retrieval score | 51.68 |

The model was chosen for three reasons:

1. **Efficiency** — 384 dimensions keeps memory manageable at scale (3.56M vectors x 384 dims x 4 bytes = ~5.2 GB for genomic_evidence alone)
2. **Quality** — Top-tier retrieval performance among small models; sufficient for domain-specific retrieval when combined with query expansion
3. **On-device inference** — Runs entirely on the DGX Spark without requiring a separate embedding service

**Limitations**

- 512-token context window means long clinical documents must be chunked (default: 2,500 characters with 200-character overlap)
- General-purpose English model — not fine-tuned for biomedical text. Query expansion compensates for vocabulary gaps.
- Does not understand structured data (tables, lab results) natively. The ingest pipeline converts structured data to natural language descriptions before embedding.

**Asymmetric Query Prefix**

BGE models use an instruction prefix for queries (but not for documents):

```python
# At query time
query_embedding = model.encode("Represent this sentence for searching relevant passages: " + query)

# At ingest time (no prefix)
doc_embedding = model.encode(document_text)
```

This asymmetric encoding improves retrieval accuracy by ~3-5% on domain benchmarks.

### Multi-Collection Parallel Search Strategy

The core architectural pattern across all 11 agents is parallel multi-collection search. A single user query fans out to all collections simultaneously.

**Search Flow**

```
User Query
    |
    v
[Query Expansion] -- add domain synonyms and related terms
    |
    v
[BGE Embedding] -- 384-dim vector
    |
    v
[Parallel Search] -- fan out to N collections simultaneously
    |   |   |   |   ... |
    v   v   v   v       v
 Coll1 Coll2 Coll3 Coll4 ... CollN
    |   |   |   |   ... |
    v   v   v   v       v
[Merge + Deduplicate] -- by content hash
    |
    v
[Score Normalization] -- weight by collection relevance
    |
    v
[Top-K Selection] -- default: 5 per collection, 30 total max
    |
    v
[Claude Synthesis] -- grounded response with citations
```

**Weighted Collection Scoring**

Each agent assigns weights to its collections. The final relevance score combines cosine similarity with collection weight:

```
final_score = cosine_similarity * collection_weight
```

The Autoimmune Agent provides a good example with weights ranging from 0.18 (clinical documents) down to 0.02 (genomic evidence and cross-disease). This ensures that highly relevant clinical notes outrank tangentially related genomic data.

**Performance on DGX Spark**

| Metric | Value |
|--------|-------|
| Single-collection search (top-5) | 2-5 ms |
| 11-collection parallel search | 12-16 ms |
| Dual retrieval (comparative mode) | ~365 ms |
| Full RAG with Claude synthesis | ~24 sec |
| End-to-end with streaming | First token ~3 sec |

Total query-to-response time stays under 5 seconds for retrieval; the majority of latency comes from LLM synthesis.

### Query Expansion with Domain Keyword Maps

Raw user queries often miss domain-specific terminology. Each agent maintains expansion maps that inject synonyms and related terms.

**Example: CAR-T Agent Expansion**

The CAR-T agent maintains 12 expansion maps covering 169 seed keywords that expand to 1,496 related terms:

| Map Category | Seed Keywords | Expanded Terms |
|-------------|---------------|----------------|
| Target Antigen | CD19, BCMA, CD22 | B-cell maturation antigen, TNFRSF17, ... |
| Toxicity | CRS, ICANS | cytokine release syndrome, neurotoxicity, tocilizumab, ... |
| Manufacturing | transduction, expansion | lentiviral vector, T-cell activation, ... |

**Expansion Logic**

```python
def expand_query(query: str) -> str:
    expanded_terms = []
    for keyword, synonyms in expansion_maps.items():
        if keyword.lower() in query.lower():
            expanded_terms.extend(synonyms[:5])  # Top 5 synonyms
    return f"{query} {' '.join(expanded_terms)}"
```

The expanded query is embedded alongside the original, and both embedding vectors are used in the parallel search. This typically improves recall by 15-25% on domain-specific queries.

### Claude LLM Prompt Engineering for Clinical RAG Synthesis

All agents use a structured prompt template for Claude synthesis:

**System Prompt Structure**

```
You are a {domain} clinical decision support AI.
You must:
1. Ground every claim in the retrieved evidence
2. Cite sources using [Collection:ID] format
3. Flag clinical alerts (contraindications, critical values)
4. Acknowledge uncertainty when evidence is limited
5. Never fabricate references or statistics

Patient context: {patient_context}
Conversation history: {last_N_turns}
```

**Evidence Injection Format**

Retrieved evidence is injected between XML-style tags:

```
<evidence>
[Literature:PMID_12345678] (score: 0.87) EGFR mutations in non-small cell lung
cancer respond to erlotinib with a median PFS of 10.4 months...

[Trials:NCT04012345] (score: 0.82) Phase III trial of osimertinib vs.
comparator in EGFR T790M-positive NSCLC...
</evidence>
```

**Key Prompt Engineering Patterns**

1. **Citation enforcement** — The system prompt requires `[Collection:Source]` format. Claude rarely hallucinates citations when the format is explicitly specified and evidence is provided.
2. **Confidence calibration** — Agents instruct Claude to state "Based on limited evidence..." or "No relevant evidence found..." rather than generating speculative answers.
3. **Streaming** — All agents support SSE streaming for real-time token delivery. The first token typically arrives within 3 seconds.

### Cross-Modal Triggers Between Agents

Agents communicate through an event-driven architecture. When one agent detects a finding that is relevant to another domain, it publishes a cross-modal event.

**Example Trigger Chains**

| Source Agent | Trigger Condition | Target Agent | Action |
|-------------|-------------------|--------------|--------|
| Imaging | Lung-RADS 4A+ nodule | Oncology | Query EGFR/ALK/ROS1 variants in genomic_evidence |
| Oncology | BRCA1/2 pathogenic variant | Clinical Trial | Search for PARP inhibitor trials |
| Rare Disease | HPO phenotype match score > 0.8 | Pharmacogenomics | Check PGx implications for candidate therapies |
| Autoimmune | Flare risk > 0.8 (imminent) | Biomarker | Pull longitudinal CRP/ESR/complement trends |
| Single-Cell | Antigen escape detected | CAR-T | Flag resistance mechanism and alternative targets |
| Cardiology | EF < 35% detected | Oncology | Check cardiotoxicity risk for proposed chemotherapy |

**Event Format**

```json
{
  "event_type": "cross_modal_trigger",
  "source_agent": "imaging_intelligence",
  "target_agent": "precision_oncology",
  "trigger": "lung_rads_4a_detected",
  "payload": {
    "finding": "Spiculated nodule 18mm RUL",
    "query": "EGFR ALK ROS1 KRAS lung adenocarcinoma variants"
  }
}
```

### Performance Optimization

**Achieving Sub-5-Second Search Latency**

1. **Pre-load collections** — All collections are loaded into memory at agent startup. Milvus `load()` is called during the health check lifecycle.
2. **Connection pooling** — A single `MilvusClient` connection is reused across requests. Connection setup (~200ms) is amortized.
3. **Embedding cache** — Frequently repeated queries (e.g., gene names, drug names) are cached with an LRU cache (default 1,024 entries).
4. **Batch embedding** — During ingest, texts are embedded in batches of 64 to maximize GPU throughput.
5. **Collection partitioning** — Large collections (genomic_evidence at 3.56M vectors) use partition keys (e.g., chromosome) for targeted search.

**Memory Budget on DGX Spark (128GB unified)**

| Component | Memory |
|-----------|--------|
| Milvus (all collections loaded) | ~25 GB |
| BGE-small-en-v1.5 model | ~0.5 GB |
| Python agent processes (11 agents) | ~8 GB |
| System + Docker overhead | ~10 GB |
| Parabricks (when running) | ~60 GB |
| **Available headroom** | **~24 GB** |

---

## Part 3: Deep Dive — All 11 Intelligence Agents

Each agent follows the same foundational architecture: FastAPI REST server, Streamlit UI, multi-collection Milvus RAG engine, BGE-small-en-v1.5 embeddings (384-dim), and Claude LLM synthesis. What differs is the domain knowledge, clinical workflows, and cross-agent integration points.

---

### 3.1 Precision Oncology Agent (:8527)

**Overview and Clinical Purpose**

The Precision Oncology Agent is a closed-loop clinical decision support system designed for molecular tumor board (MTB) workflows. Given a patient's somatic and germline VCF, it identifies actionable variants, matches them against oncology knowledge bases, ranks candidate therapies, surfaces relevant clinical trials, and assembles a structured MTB-ready report.

**Milvus Collections (11)**

| Collection | Description | Source |
|-----------|-------------|--------|
| `onco_literature` | PubMed/PMC chunks by cancer type | PubMed E-utils |
| `onco_trials` | Clinical trial summaries with biomarker criteria | ClinicalTrials.gov |
| `onco_variants` | Actionable somatic/germline variants | CIViC, OncoKB |
| `onco_biomarkers` | Predictive/prognostic biomarkers and assays | CIViC, literature |
| `onco_therapies` | Approved and investigational therapies with MOA | OncoKB, FDA labels |
| `onco_pathways` | Signaling pathways, cross-talk, druggable nodes | KEGG, Reactome |
| `onco_guidelines` | NCCN/ASCO/ESMO guideline recommendations | Guideline PDFs |
| `onco_resistance` | Resistance mechanisms and bypass strategies | CIViC, literature |
| `onco_outcomes` | Real-world treatment outcome records | De-identified RWD |
| `onco_cases` | De-identified patient case snapshots | Synthetic/RWD |
| `genomic_evidence` | VCF-derived evidence (shared, read-only) | Stage 1 pipeline |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8527/health

# Free-text oncology RAG query
curl -X POST http://localhost:8527/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What targeted therapies are available for EGFR L858R in NSCLC?"}'

# Create a patient case
curl -X POST http://localhost:8527/cases \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "P001", "cancer_type": "NSCLC", "variants": ["EGFR L858R"]}'

# Run full MTB analysis
curl -X POST http://localhost:8527/cases/P001/analyze

# Generate MTB report
curl -X POST http://localhost:8527/reports/mtb \
  -H "Content-Type: application/json" \
  -d '{"case_id": "P001"}'

# Match molecular profile to clinical trials
curl -X POST http://localhost:8527/trials/match \
  -H "Content-Type: application/json" \
  -d '{"variants": ["EGFR L858R"], "cancer_type": "NSCLC"}'
```

**Clinical Workflow**

1. Clinician uploads somatic/germline VCF or enters variant list
2. Agent annotates variants against CIViC and OncoKB with evidence-level tiering (Level 1-4)
3. TherapyRanker scores approved and investigational therapies using variant profile, biomarkers, resistance mechanisms, and guideline concordance
4. TrialMatcher identifies eligible clinical trials based on molecular profile and cancer type
5. MTB report is generated with tiered evidence citations and FHIR R4 export

**Cross-Agent Integration**

- Receives genomic evidence from Stage 1 pipeline via shared `genomic_evidence` collection
- Triggers Clinical Trial Agent when novel actionable variants are identified
- Receives cardiotoxicity alerts from Cardiology Agent for proposed chemotherapy regimens
- Connects to Drug Discovery pipeline (Stage 3) for candidate molecule docking

---

### 3.2 CAR-T Intelligence Agent (:8522)

**Overview and Clinical Purpose**

The CAR-T Intelligence Agent breaks down data silos across the 5 stages of CAR-T cell therapy development: target selection, construct design, manufacturing, clinical evaluation, and post-market surveillance. It features a unique comparative analysis mode that auto-detects "X vs Y" queries and produces structured side-by-side comparisons.

**Milvus Collections (11)**

| Collection | Records | Source |
|-----------|---------|--------|
| `cart_literature` | 5,047 | PubMed abstracts via NCBI E-utilities |
| `cart_trials` | 973 | ClinicalTrials.gov API v2 |
| `cart_constructs` | 6 | 6 FDA-approved CAR-T products |
| `cart_assay_results` | 45 | Curated from landmark papers (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1) |
| `cart_manufacturing` | 30 | CMC/process data (transduction, expansion, release, cryo, logistics) |
| `cart_safety` | — | Pharmacovigilance, CRS/ICANS profiles |
| `cart_biomarkers` | — | CRS prediction, exhaustion monitoring |
| `cart_regulatory` | — | Approval timelines, post-marketing requirements |
| `cart_sequences` | — | Molecular binding, scFv sequences |
| `cart_rwe` | — | Registry outcomes, real-world data |
| `genomic_evidence` | 3.56M | Shared from Stage 2 RAG pipeline (read-only) |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8522/health

# Cross-functional RAG query
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"}'

# Comparative analysis (auto-detected)
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"}'

# Collection statistics
curl http://localhost:8522/collections/stats
```

**Clinical Workflow**

1. User submits a query about CAR-T development (any stage)
2. Comparative detection engine checks for "vs/versus/compare" patterns
3. If comparative: parse entities, resolve against knowledge graph (25 antigens, 6 products, 39+ aliases), run dual retrieval
4. If standard: query expansion (12 maps, 169 keywords to 1,496 terms), parallel search across 11 collections
5. Claude synthesizes grounded response with clickable PubMed and ClinicalTrials.gov citations

**Cross-Agent Integration**

- Receives antigen escape alerts from Single-Cell Agent
- Queries Oncology Agent for tumor mutational burden context
- Shares CRS/ICANS toxicity data with Biomarker Agent for monitoring protocols

---

### 3.3 Precision Biomarker Agent (:8529)

**Overview and Clinical Purpose**

The Precision Biomarker Agent interprets patient biomarker panels with genotype awareness. It estimates biological age using PhenoAge/GrimAge algorithms, detects disease trajectories across 6 categories, provides pharmacogenomic profiling for 7 key pharmacogenes, and generates genotype-adjusted reference ranges.

**Milvus Collections (11)**

| Collection | Description | Records |
|-----------|-------------|---------|
| `biomarker_reference` | Reference biomarker definitions and ranges | ~60 |
| `biomarker_genetic_variants` | Genetic variants affecting biomarker levels | ~30 |
| `biomarker_pgx_rules` | CPIC pharmacogenomic dosing rules | ~50 |
| `biomarker_disease_trajectories` | Disease progression stage definitions | ~30 |
| `biomarker_clinical_evidence` | Published clinical evidence | ~40 |
| `biomarker_nutrition` | Genotype-aware nutrition guidelines | ~25 |
| `biomarker_drug_interactions` | Gene-drug interactions | ~35 |
| `biomarker_aging_markers` | Epigenetic aging clock markers | ~15 |
| `biomarker_genotype_adjustments` | Genotype-based reference range adjustments | ~20 |
| `biomarker_monitoring` | Condition-specific monitoring protocols | ~15 |
| `genomic_evidence` | Shared genomic variants (read-only) | 3.56M |

**Key API Endpoints**

```bash
# Full patient analysis (biomarkers + genotypes)
curl -X POST http://localhost:8529/analyze \
  -H "Content-Type: application/json" \
  -d '{"biomarkers": {"CRP": 2.1, "HbA1c": 6.8, "LDL": 145}, "genotypes": {"CYP2D6": "*1/*4"}}'

# Biological age calculation
curl -X POST http://localhost:8529/biological-age \
  -H "Content-Type: application/json" \
  -d '{"age": 55, "albumin": 4.2, "creatinine": 0.9, "glucose": 105, "crp": 1.8}'

# Pharmacogenomic profile
curl -X POST http://localhost:8529/pgx \
  -H "Content-Type: application/json" \
  -d '{"diplotypes": {"CYP2D6": "*1/*4", "CYP2C19": "*1/*2"}}'

# Disease risk trajectory
curl -X POST http://localhost:8529/disease-risk \
  -H "Content-Type: application/json" \
  -d '{"biomarkers": {"HbA1c": 6.2, "fasting_glucose": 115}}'
```

**Clinical Workflow**

1. Patient biomarker panel and genotype data are submitted
2. Biological Age Engine computes PhenoAge and GrimAge estimates from 9 routine blood biomarkers
3. Disease Trajectory Analyzer screens 6 disease categories with genotype-stratified thresholds
4. Pharmacogenomic Mapper interprets star alleles for CYP2D6, CYP2C19, CYP2C9, CYP3A5, SLCO1B1, VKORC1, MTHFR, plus HLA-B*57:01
5. Genotype Adjustment Engine modifies standard reference ranges (e.g., PNPLA3, TCF7L2, APOE variants)
6. 12-section clinical report generated with PDF and FHIR R4 export

**Cross-Agent Integration**

- Receives inflammation monitoring requests from Autoimmune Agent during flare events
- Provides biological age context to Cardiology Agent for ASCVD risk refinement
- Shares pharmacogenomic profiles with Pharmacogenomics Agent for detailed dosing guidance

---

### 3.4 Cardiology Intelligence Agent (:8126)

**Overview and Clinical Purpose**

The Cardiology Intelligence Agent synthesizes cardiac imaging, electrophysiology, hemodynamics, heart failure management, valvular disease, preventive cardiology, interventional data, and cardio-oncology surveillance into ACC/AHA/ESC guideline-aligned recommendations. It includes 6 validated risk calculators and 11 clinical workflows.

**Milvus Collections (13)**

| Collection | Description |
|-----------|-------------|
| `cardio_literature` | Published cardiovascular research, reviews, meta-analyses |
| `cardio_trials` | Cardiovascular clinical trials and landmark results |
| `cardio_imaging` | Cardiac imaging protocols, findings, measurements |
| `cardio_electrophysiology` | ECG interpretation, arrhythmia classification, EP data |
| `cardio_heart_failure` | HF classification, GDMT protocols, management algorithms |
| `cardio_valvular` | Valvular heart disease assessment and intervention criteria |
| `cardio_prevention` | Preventive cardiology: risk stratification, lipid management |
| `cardio_interventional` | Interventional procedures, techniques, outcomes |
| `cardio_oncology` | Cardio-oncology surveillance, cardiotoxicity detection |
| `cardio_devices` | FDA-cleared cardiovascular AI devices, implantables |
| `cardio_guidelines` | ACC/AHA/ESC/HRS clinical practice guidelines |
| `cardio_hemodynamics` | Catheterization data, pressure tracings, derived calculations |
| `genomic_evidence` | Shared genomic evidence (read-only, 3.56M variants) |

**6 Risk Calculators**

| Calculator | Clinical Use |
|-----------|-------------|
| ASCVD (PCE) | 10-year atherosclerotic cardiovascular disease risk |
| HEART Score | Chest pain risk stratification in ED |
| CHA2DS2-VASc | Atrial fibrillation stroke risk |
| HAS-BLED | Anticoagulation bleeding risk |
| MAGGIC | Heart failure mortality risk |
| EuroSCORE II | Cardiac surgical mortality risk |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8126/health

# RAG query
curl -X POST http://localhost:8126/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Optimal GDMT titration for HFrEF with EF 25%"}'

# ASCVD risk calculation
curl -X POST http://localhost:8126/risk/ascvd \
  -H "Content-Type: application/json" \
  -d '{"age": 55, "sex": "male", "total_cholesterol": 220, "hdl": 45, "systolic_bp": 140, "diabetes": false, "smoker": false}'

# CHA2DS2-VASc score
curl -X POST http://localhost:8126/risk/chadsvasc \
  -H "Content-Type: application/json" \
  -d '{"age": 72, "sex": "female", "chf": true, "hypertension": true, "diabetes": false, "stroke_history": false, "vascular_disease": true}'
```

**Clinical Workflow (11 workflows)**

Workflows span coronary artery disease assessment, heart failure classification and GDMT optimization, valvular disease quantification, arrhythmia detection, cardiac MRI tissue characterization, stress test interpretation, preventive risk stratification, cardio-oncology surveillance, acute decompensated HF, post-MI management, and myocarditis/pericarditis evaluation.

**Cross-Agent Integration**

- Sends cardiotoxicity alerts to Oncology Agent when chemotherapy candidates show cardiac risk
- Receives biomarker trends (BNP, troponin) from Biomarker Agent
- Queries genomic_evidence for cardiomyopathy-associated variants (TTN, LMNA, MYH7)

---

### 3.5 Neurology Intelligence Agent (:8528)

**Overview and Clinical Purpose**

The Neurology Intelligence Agent supports 8 clinical workflows covering acute stroke triage, dementia evaluation, epilepsy classification, brain tumor grading, MS monitoring, Parkinson's assessment, headache classification, and neuromuscular evaluation. It integrates guidelines from AAN, AHA/ASA, ILAE, ICHD-3, WHO CNS 2021, McDonald 2017, and MDS criteria.

**10 Clinical Scale Calculators**

| Scale | Purpose |
|-------|---------|
| NIHSS | National Institutes of Health Stroke Scale |
| GCS | Glasgow Coma Scale |
| MoCA | Montreal Cognitive Assessment |
| MDS-UPDRS Part III | Parkinson's motor examination |
| EDSS | Expanded Disability Status Scale (MS) |
| mRS | Modified Rankin Scale (functional outcome) |
| HIT-6 | Headache Impact Test |
| ALSFRS-R | ALS Functional Rating Scale — Revised |
| ASPECTS | Alberta Stroke Programme Early CT Score |
| Hoehn-Yahr | Parkinson's disease staging |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8528/health

# RAG query
curl -X POST http://localhost:8528/v1/neuro/query \
  -H "Content-Type: application/json" \
  -d '{"question": "McDonald 2017 criteria for MS diagnosis with one clinical attack"}'

# Clinical scale calculation
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{"scale": "NIHSS", "items": {"consciousness": 0, "gaze": 1, "visual_fields": 0, "facial_palsy": 2, "motor_arm_left": 3}}'

# Acute stroke triage
curl -X POST http://localhost:8528/v1/neuro/stroke/triage \
  -H "Content-Type: application/json" \
  -d '{"onset_time": "2026-03-29T08:30:00", "nihss_score": 14, "ct_findings": "no hemorrhage"}'

# Dementia evaluation
curl -X POST http://localhost:8528/v1/neuro/dementia/evaluate \
  -H "Content-Type: application/json" \
  -d '{"moca_score": 18, "age": 72, "symptoms": ["memory loss", "word finding difficulty"]}'
```

**Clinical Workflow**

1. Clinician selects a workflow (stroke, dementia, epilepsy, tumor, MS, Parkinson's, headache, neuromuscular)
2. Structured input forms collect workflow-specific parameters
3. Clinical scale calculators provide standardized scoring
4. RAG engine retrieves guideline-concordant evidence across 12 collections
5. Multi-format report generated (Markdown, JSON, PDF, FHIR R4 DiagnosticReport)

**Cross-Agent Integration**

- Queries Imaging Agent for MRI findings (lesion counts, brain volumetrics)
- Receives pharmacogenomic context from PGx Agent for antiepileptic drug selection
- Shares stroke triage events with Cardiology Agent for cardiac workup (AF screening)

---

### 3.6 Precision Autoimmune Agent (:8532)

**Overview and Clinical Purpose**

The Precision Autoimmune Agent addresses the diagnostic odyssey that autoimmune patients face — often 3-6 years across multiple specialists before diagnosis. It interprets autoantibody panels, HLA typing, biomarker trends, and genomic data across 13 autoimmune conditions. It features disease activity scoring (DAS28-CRP, SLEDAI-2K, CDAI, BASDAI), flare prediction, biologic therapy recommendations with PGx context, and diagnostic odyssey analysis.

**13 Supported Conditions**

Rheumatoid arthritis, SLE, multiple sclerosis, type 1 diabetes, inflammatory bowel disease, psoriasis/psoriatic arthritis, ankylosing spondylitis, Sjogren's syndrome, systemic sclerosis, myasthenia gravis, celiac disease, Graves' disease, Hashimoto's thyroiditis. Also detects the POTS/hEDS/MCAS triad and overlap syndromes.

**Milvus Collections (14)**

| Collection | Weight | Description |
|-----------|--------|-------------|
| `autoimmune_clinical_documents` | 0.18 | Progress notes, discharge summaries |
| `autoimmune_patient_labs` | 0.14 | Lab results with flags and reference ranges |
| `autoimmune_autoantibody_panels` | 0.12 | ANA, anti-dsDNA, anti-CCP, RF, and 14+ types |
| `autoimmune_hla_associations` | 0.08 | HLA allele-disease risk (50+ alleles with odds ratios) |
| `autoimmune_disease_criteria` | 0.08 | ACR/EULAR classification criteria |
| `autoimmune_disease_activity` | 0.07 | DAS28-CRP, SLEDAI-2K, CDAI, BASDAI scoring |
| `autoimmune_flare_patterns` | 0.06 | Biomarker patterns preceding flares |
| `autoimmune_biologic_therapies` | 0.06 | TNF inhibitors, anti-CD20, IL-6R, JAK inhibitors |
| `autoimmune_pgx_rules` | 0.04 | Pharmacogenomic rules (CYP2C19, FCGR3A, HLA) |
| `autoimmune_clinical_trials` | 0.05 | Active and recent clinical trials |
| `autoimmune_literature` | 0.05 | Published literature with year filtering |
| `autoimmune_patient_timelines` | 0.03 | Longitudinal patient event timelines |
| `autoimmune_cross_disease` | 0.02 | Overlap syndromes, cross-disease mechanisms |
| `genomic_evidence` | 0.02 | Shared genomic evidence (read-only) |

**Key API Endpoints**

```bash
# Full patient analysis
curl -X POST http://localhost:8532/analyze \
  -H "Content-Type: application/json" \
  -d '{"antibodies": {"ANA": "positive", "anti_dsDNA": 85}, "hla": ["B*27:05"], "biomarkers": {"CRP": 18.5, "ESR": 42}}'

# Disease activity scoring
curl -X POST http://localhost:8532/disease-activity \
  -H "Content-Type: application/json" \
  -d '{"disease": "RA", "score_type": "DAS28-CRP", "tender_joints": 6, "swollen_joints": 4, "crp": 2.1, "patient_global": 55}'

# Flare prediction
curl -X POST http://localhost:8532/flare-risk \
  -H "Content-Type: application/json" \
  -d '{"disease": "SLE", "biomarkers": {"CRP": 12, "C3": 65, "C4": 8, "anti_dsDNA": 120}}'

# Ingest demo patient data
curl -X POST http://localhost:8532/ingest/demo-data
```

**Clinical Workflow**

1. Clinical documents (progress notes, lab reports, imaging) are ingested and chunked
2. Autoantibody panel interpretation maps positives to disease associations with sensitivity/specificity
3. HLA analysis evaluates typing against 50+ allele-disease associations (e.g., HLA-B*27:05 and AS at OR=87.4)
4. Disease activity scoring calculates standardized scores with remission/low/moderate/high/very high levels
5. Flare prediction analyzes CRP, ESR, IL-6, complement C3/C4, calprotectin patterns on a 0-1 risk scale
6. Biologic therapy recommendations filtered by diagnosis and PGx context
7. Diagnostic odyssey analysis surfaces patterns missed across fragmented multi-specialist records

**Cross-Agent Integration**

- Receives longitudinal biomarker trends from Biomarker Agent
- Requests imaging assessment from Imaging Agent for joint/organ evaluation
- Publishes diagnosis events to other agents via event bus

---

### 3.7 Rare Disease Diagnostic Agent (:8134)

**Overview and Clinical Purpose**

The Rare Disease Diagnostic Agent provides differential diagnosis, ACMG/AMP variant interpretation (implementing a subset of 28 ACMG criteria), HPO-based phenotype matching with information content (IC) weighted similarity scoring, therapeutic option search (orphan drugs, gene therapy, enzyme replacement), and clinical trial eligibility assessment. It covers 13 disease categories with 88 diseases in its knowledge base.

**Milvus Collections (14)**

| Collection | Description |
|-----------|-------------|
| `rd_phenotypes` | HPO-coded phenotype descriptions |
| `rd_diseases` | Disease definitions (OMIM, Orphanet) |
| `rd_genes` | Gene-disease associations |
| `rd_variants` | Pathogenic variant records |
| `rd_literature` | Rare disease literature |
| `rd_trials` | Rare disease clinical trials |
| `rd_therapies` | Orphan drugs, gene therapy, ERT |
| `rd_case_reports` | Published case reports |
| `rd_guidelines` | Clinical management guidelines |
| `rd_pathways` | Molecular pathways |
| `rd_registries` | Patient registries |
| `rd_natural_history` | Natural history studies |
| `rd_newborn_screening` | Newborn screening panels |
| `genomic_evidence` | Shared genomic evidence (read-only) |

**Key API Endpoints**

```bash
# Differential diagnosis from HPO terms
curl -X POST http://localhost:8134/v1/diagnostic/diagnose \
  -H "Content-Type: application/json" \
  -d '{"hpo_terms": ["HP:0001250", "HP:0001263", "HP:0002079"], "age_of_onset": "infantile"}'

# ACMG variant interpretation
curl -X POST http://localhost:8134/v1/diagnostic/variants/interpret \
  -H "Content-Type: application/json" \
  -d '{"variant": "NM_000492.4(CFTR):c.1521_1523del", "gene": "CFTR", "consequence": "in_frame_deletion"}'

# HPO-to-disease phenotype matching
curl -X POST http://localhost:8134/v1/diagnostic/phenotype/match \
  -H "Content-Type: application/json" \
  -d '{"hpo_terms": ["HP:0001250", "HP:0001263"], "max_results": 10}'

# Therapeutic option search
curl -X POST http://localhost:8134/v1/diagnostic/therapy/search \
  -H "Content-Type: application/json" \
  -d '{"disease": "Gaucher disease type 1", "therapy_types": ["ERT", "SRT"]}'

# Clinical trial eligibility
curl -X POST http://localhost:8134/v1/diagnostic/trial/match \
  -H "Content-Type: application/json" \
  -d '{"disease": "Duchenne muscular dystrophy", "age": 8, "gene": "DMD"}'
```

**Clinical Workflow**

1. Patient phenotypes entered as HPO terms via the 5-tab Streamlit UI
2. HPO-to-disease matching uses IC-weighted similarity against OMIM and Orphanet
3. Differential diagnosis ranked by phenotype match score with confidence levels
4. Candidate variants classified using ACMG/AMP criteria (PVS1, PM1-PM6, PP1-PP5, BA1, BS1-BS4, BP1-BP7)
5. Therapeutic options searched across orphan drugs, gene therapy, and enzyme replacement
6. Report exported as Markdown, JSON, PDF, FHIR R4, or GA4GH Phenopacket v2

**Cross-Agent Integration**

- Queries Pharmacogenomics Agent for drug metabolism context in rare disease therapeutics
- Shares phenotype-genotype correlations with Biomarker Agent
- Triggers Clinical Trial Agent for rare disease trial matching

---

### 3.8 Pharmacogenomics Intelligence Agent (:8107)

**Overview and Clinical Purpose**

The Pharmacogenomics (PGx) Intelligence Agent translates patient genotype data into actionable prescribing guidance. It interprets star alleles for 14 pharmacogenes, applies CPIC and DPWG guidelines, screens for HLA-mediated adverse drug reactions, models phenoconversion from drug-drug interactions, and provides genotype-guided dosing algorithms.

**Milvus Collections (15)**

| Collection | Description |
|-----------|-------------|
| `pgx_gene_reference` | Pharmacogene star allele definitions and activity scores |
| `pgx_drug_guidelines` | CPIC/DPWG clinical prescribing guidelines |
| `pgx_drug_interactions` | Drug-gene interaction records (PharmGKB) |
| `pgx_hla_hypersensitivity` | HLA-mediated adverse drug reaction screening |
| `pgx_phenoconversion` | Metabolic phenoconversion via DDI |
| `pgx_dosing_algorithms` | Genotype-guided dosing algorithms and formulas |
| `pgx_clinical_evidence` | Published PGx clinical evidence and outcomes |
| `pgx_population_data` | Population-specific allele frequency data |
| `pgx_clinical_trials` | PGx-related clinical trials |
| `pgx_fda_labels` | FDA pharmacogenomic labeling information |
| `pgx_drug_alternatives` | Genotype-guided therapeutic alternatives |
| `pgx_patient_profiles` | Patient diplotype-phenotype profiles |
| `pgx_implementation` | Clinical PGx implementation programs |
| `pgx_education` | PGx educational resources and guidelines |
| `genomic_evidence` | Shared genomic evidence (read-only) |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8107/health

# Full PGx profile from diplotypes
curl -X POST http://localhost:8107/profile \
  -H "Content-Type: application/json" \
  -d '{"diplotypes": {"CYP2D6": "*1/*4", "CYP2C19": "*1/*2", "CYP2C9": "*1/*3", "SLCO1B1": "*1a/*5"}}'

# Drug-specific dosing guidance
curl -X POST http://localhost:8107/dosing \
  -H "Content-Type: application/json" \
  -d '{"drug": "codeine", "diplotype": {"CYP2D6": "*4/*4"}, "indication": "pain"}'

# HLA screening
curl -X POST http://localhost:8107/hla-screen \
  -H "Content-Type: application/json" \
  -d '{"hla_alleles": ["B*57:01", "B*58:01"], "proposed_drugs": ["abacavir", "allopurinol"]}'

# RAG query
curl -X POST http://localhost:8107/query \
  -H "Content-Type: application/json" \
  -d '{"question": "CPIC guidelines for CYP2D6 poor metabolizers prescribed tramadol"}'
```

**Clinical Workflow**

1. Patient diplotype data submitted (from VCF star allele calling or manual entry)
2. Star alleles mapped to metabolizer phenotypes (ultra-rapid, normal, intermediate, poor)
3. CPIC/DPWG guidelines matched for each gene-drug pair
4. Phenoconversion check: current medications that may alter metabolizer status via enzyme inhibition/induction
5. HLA screening for hypersensitivity risk (HLA-B*57:01 and abacavir, HLA-B*58:01 and allopurinol, etc.)
6. Dosing algorithms applied (e.g., IWPC warfarin algorithm)
7. Alternative drug recommendations for poor/ultra-rapid metabolizers

**Cross-Agent Integration**

- Provides PGx context to all agents that recommend drug therapies
- Receives dosing queries from Oncology Agent for chemotherapy metabolism
- Shares metabolizer status with Biomarker Agent for genotype-adjusted reference ranges

---

### 3.9 Imaging Intelligence Agent (:8524)

**Overview and Clinical Purpose**

The Imaging Intelligence Agent provides clinical decision support for radiology through a RAG knowledge system, 4 NVIDIA NIM microservices (VISTA-3D, MAISI, VILA-M3, Llama-3), and 4 reference clinical workflows. It supports Orthanc DICOM auto-ingestion, cross-modal genomics enrichment, NVIDIA FLARE federated learning, and FHIR R4 interoperability.

**NVIDIA NIM Integration**

| NIM Service | Port | Capability |
|------------|------|------------|
| VISTA-3D | 8530 | 3D medical image segmentation (132 anatomical classes) |
| MAISI | 8531 | Synthetic CT volume generation with paired segmentation masks |
| VILA-M3 | 8532 | Vision-language model for radiology image understanding |
| Llama-3 8B | 8520 | Clinical reasoning and report generation |

**Milvus Collections (11)**

| Collection | Content |
|-----------|---------|
| `imaging_literature` | 2,678 PubMed research papers |
| `imaging_trials` | 12 clinical trials |
| `imaging_findings` | Imaging finding templates and patterns |
| `imaging_protocols` | Acquisition protocols and parameters |
| `imaging_devices` | FDA-cleared AI/ML medical devices |
| `imaging_anatomy` | Anatomical structure references |
| `imaging_benchmarks` | Model performance benchmarks |
| `imaging_guidelines` | ACR, RSNA, NCCN guidelines |
| `imaging_report_templates` | Structured radiology report templates |
| `imaging_datasets` | Public imaging datasets (TCIA, PhysioNet) |
| `genomic_evidence` | Shared from Stage 2 (read-only, 3.56M) |

**4 Reference Workflows**

| Workflow | Modality | Target Latency |
|---------|----------|---------------|
| CT Head Hemorrhage Triage | CT | < 90 sec |
| CT Chest Lung Nodule Tracking | CT | < 5 min |
| CXR Rapid Findings | X-ray | < 30 sec |
| MRI Brain MS Lesion Tracking | MRI | < 5 min |

**Key API Endpoints**

```bash
# Meta-agent question answering
curl -X POST http://localhost:8524/api/ask \
  -H "Content-Type: application/json" \
  -d '{"question": "What is the sensitivity of AI-assisted CXR for pneumothorax detection?"}'

# Run VISTA-3D segmentation
curl -X POST http://localhost:8524/nim/vista3d/segment \
  -H "Content-Type: application/json" \
  -d '{"image_path": "/data/ct_scan.nii.gz", "target_classes": ["liver", "spleen", "kidney"]}'

# Execute a workflow
curl -X POST http://localhost:8524/workflow/ct_head_hemorrhage/run \
  -H "Content-Type: application/json" \
  -d '{"study_path": "/data/ct_head/", "patient_id": "P001"}'

# DICOM webhook (auto-triggered by Orthanc)
curl -X POST http://localhost:8524/events/dicom-webhook \
  -H "Content-Type: application/json" \
  -d '{"study_uid": "1.2.3.4.5", "modality": "CT", "body_part": "HEAD"}'
```

**Cross-Agent Integration**

- When Lung-RADS 4A+ nodules are detected, automatically queries `genomic_evidence` for EGFR/ALK/ROS1/KRAS variants
- Receives MRI requests from Neurology Agent for MS lesion tracking
- Provides cardiac imaging context to Cardiology Agent

---

### 3.10 Single-Cell Intelligence Agent (:8540)

**Overview and Clinical Purpose**

The Single-Cell Intelligence Agent provides cell-type annotation, tumor microenvironment (TME) profiling, drug response prediction, subclonal architecture analysis, spatial transcriptomics niche mapping, trajectory inference, ligand-receptor interaction analysis, biomarker discovery, CAR-T target validation, and treatment monitoring. Its knowledge base covers 57 cell types, 30 drugs, 75 markers, and 4 spatial platforms (Visium, MERFISH, Xenium, CosMx).

**Milvus Collections (13)**

| Collection | Description |
|-----------|-------------|
| `sc_cell_types` | Cell type reference profiles and markers |
| `sc_markers` | Marker gene signatures per cell type |
| `sc_literature` | Single-cell genomics literature |
| `sc_trials` | Clinical trials involving single-cell analysis |
| `sc_drugs` | Drug response signatures (GDSC/DepMap) |
| `sc_pathways` | Pathway activity signatures |
| `sc_spatial` | Spatial transcriptomics references |
| `sc_trajectories` | Differentiation trajectory templates |
| `sc_interactions` | Ligand-receptor pair databases (CellPhoneDB/NicheNet) |
| `sc_tme` | TME classification profiles |
| `sc_clinical` | Clinical correlation data |
| `sc_methods` | Computational method references |
| `genomic_evidence` | Shared genomic evidence (read-only) |

**Key API Endpoints**

```bash
# Cell type annotation
curl -X POST http://localhost:8540/v1/sc/annotate \
  -H "Content-Type: application/json" \
  -d '{"markers": {"CD3D": 5.2, "CD8A": 4.1, "GZMB": 3.8}, "tissue": "tumor"}'

# TME profiling
curl -X POST http://localhost:8540/v1/sc/tme-profile \
  -H "Content-Type: application/json" \
  -d '{"cell_type_proportions": {"CD8_T_cell": 0.15, "Macrophage": 0.25, "Treg": 0.08, "Fibroblast": 0.30}}'

# Drug response prediction
curl -X POST http://localhost:8540/v1/sc/drug-response \
  -H "Content-Type: application/json" \
  -d '{"cell_type": "CD8_T_cell", "drug": "pembrolizumab", "expression_profile": {"PD1": 4.2, "PDL1": 3.1}}'

# CAR-T target validation
curl -X POST http://localhost:8540/v1/sc/cart-validate \
  -H "Content-Type: application/json" \
  -d '{"target": "CD19", "tumor_type": "B-ALL", "expression_data": {"tumor": 0.95, "normal_b_cells": 0.88, "other_tissues": 0.02}}'

# Spatial niche mapping
curl -X POST http://localhost:8540/v1/sc/spatial-niche \
  -H "Content-Type: application/json" \
  -d '{"platform": "Visium", "cell_types": {"spot_1": "T_cell", "spot_2": "Macrophage", "spot_3": "Tumor"}}'
```

**Clinical Workflow**

1. Single-cell expression data submitted (gene expression matrix or marker panel)
2. Multi-strategy cell type annotation: reference-based, marker-based, and LLM-augmented
3. TME classification: hot, cold, excluded, or immunosuppressive microenvironment
4. Drug response prediction using GDSC/DepMap cell-type-resolved signatures
5. Subclonal architecture analysis with CNV-based detection and antigen escape risk
6. Spatial niche identification for spatial transcriptomics data
7. Reports exported in Markdown, JSON, PDF, or FHIR R4

**Cross-Agent Integration**

- Sends antigen escape alerts to CAR-T Agent when off-tumor expression detected
- Provides TME context to Oncology Agent for immunotherapy selection
- Shares cell-type-resolved drug sensitivity data with Pharmacogenomics Agent

---

### 3.11 Clinical Trial Intelligence Agent (:8538)

**Overview and Clinical Purpose**

The Clinical Trial Intelligence Agent provides AI-driven support across the full clinical trial lifecycle: protocol optimization, patient-trial matching, site selection, eligibility optimization, adaptive design evaluation, safety signal detection, regulatory document generation, competitive intelligence, diversity assessment, and decentralized trial planning. It supports 10 specialized workflows plus general RAG query capability.

**Milvus Collections (14)**

| Collection | Description |
|-----------|-------------|
| `trial_protocols` | Protocol designs, amendments, SOAs |
| `trial_eligibility` | Inclusion/exclusion criteria |
| `trial_endpoints` | Primary, secondary, exploratory endpoints |
| `trial_sites` | Site feasibility, enrollment history |
| `trial_investigators` | Investigator profiles, experience |
| `trial_results` | Published trial results, CSRs |
| `trial_regulatory` | FDA/EMA guidance, IND/NDA data |
| `trial_literature` | PubMed clinical trial publications |
| `trial_biomarkers` | Biomarker-driven trial designs |
| `trial_safety` | Safety data, DSMB reports, AE databases |
| `trial_rwe` | Real-world evidence, registries |
| `trial_adaptive` | Adaptive design references |
| `trial_guidelines` | ICH, FDA, EMA guidelines |
| `genomic_evidence` | Shared genomic evidence (cross-agent) |

**10 Specialized Workflows**

| Workflow | Description |
|---------|-------------|
| Protocol Design | Complexity scoring, SOA review, endpoint optimization |
| Patient Matching | AI eligibility screening with genomic/biomarker integration |
| Site Selection | Feasibility scoring, enrollment forecasting, diversity metrics |
| Eligibility Optimization | Population impact modeling, competitor benchmarking |
| Adaptive Design | Bayesian interim analysis, dose-response, futility assessment |
| Safety Signal | Disproportionality analysis (PRR/ROR), causality assessment |
| Regulatory Docs | IND, CSR, briefing doc generation (FDA/EMA/PMDA) |
| Competitive Intel | Landscape analysis, enrollment race tracking |
| Diversity Assessment | FDA diversity guidance compliance, gap analysis |
| Decentralized Planning | DCT component feasibility, hybrid model design |

**Key API Endpoints**

```bash
# Health check
curl http://localhost:8538/health

# RAG query
curl -X POST http://localhost:8538/v1/trial/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What adaptive designs are used in Phase II oncology trials?"}'

# Patient-trial matching
curl -X POST http://localhost:8538/v1/trial/workflow/patient_matching/run \
  -H "Content-Type: application/json" \
  -d '{"patient_profile": {"age": 55, "diagnosis": "NSCLC", "biomarkers": ["EGFR L858R"], "prior_therapies": ["carboplatin"]}}'

# Protocol design review
curl -X POST http://localhost:8538/v1/trial/workflow/protocol_design/run \
  -H "Content-Type: application/json" \
  -d '{"protocol": {"phase": "II", "indication": "NSCLC", "primary_endpoint": "ORR", "sample_size": 120}}'

# Safety signal detection
curl -X POST http://localhost:8538/v1/trial/workflow/safety_signal/run \
  -H "Content-Type: application/json" \
  -d '{"drug": "experimental_agent_X", "adverse_events": [{"event": "hepatotoxicity", "count": 8, "total_exposed": 200}]}'

# List available workflows
curl http://localhost:8538/workflows
```

**Clinical Workflow**

1. User selects a workflow or submits a free-form query
2. Workflow engine routes to the appropriate analysis module
3. RAG engine retrieves evidence from 14 collections with domain-specific query expansion
4. For patient matching: genomic and biomarker data cross-referenced against trial eligibility criteria
5. For protocol design: complexity scoring against benchmark protocols, SOA review
6. For safety signals: disproportionality analysis (PRR/ROR) with historical comparator data
7. Regulatory document templates generated with auto-populated sections

**Cross-Agent Integration**

- Receives trial matching requests from Oncology Agent when actionable variants are found
- Receives rare disease trial queries from Rare Disease Agent
- Accesses genomic_evidence for biomarker-driven trial design support
- Integrates with Biomarker Agent for eligibility biomarker interpretation

---

## Part 4: Advanced Drug Discovery

### The 10-Stage Therapeutic Discovery Pipeline

The Drug Discovery Pipeline transforms a clinically actionable variant into ranked therapeutic candidates. Each stage feeds the next, producing an auditable chain of evidence from gene to molecule.

| Stage | Name | Technology | Typical Time |
|-------|------|------------|-------------|
| 1 | Variant Reception | VCF parser + ClinVar lookup | <1 s |
| 2 | Target Identification | UniProt mapping, PDB structure fetch | 2-5 s |
| 3 | Binding Site Analysis | fpocket cavity detection, druggability scoring | 5-10 s |
| 4 | Seed Compound Retrieval | ChEMBL/DrugBank similarity search | 3-8 s |
| 5 | Molecular Generation | BioNeMo MolMIM conditional generation | 15-30 s |
| 6 | Chemical Feasibility | RDKit Lipinski, QED, SA score | 1-2 s |
| 7 | Pediatric Safety Filtering | 6-filter pediatric safety gate | <1 s |
| 8 | Binding Prediction | BioNeMo DiffDock pose sampling | 30-90 s |
| 9 | Composite Scoring | Weighted multi-objective ranking | <1 s |
| 10 | Report Generation | Claude LLM narrative + evidence assembly | 5-10 s |

**End-to-end time per variant:** 60-150 seconds depending on the number of generated molecules and DiffDock sampling depth.

### MolMIM Molecular Generation

MolMIM (Molecular Masked Inverse Modeling) is a BioNeMo NIM that generates novel molecules conditioned on a seed SMILES string. The HCLS AI Factory sends generation requests to the local MolMIM endpoint.

**Seed SMILES Selection**

The pipeline selects seed compounds from the top ChEMBL/DrugBank hits in Stage 4. Selection criteria:

- Known activity against the target gene product (IC50 or Ki reported)
- Molecular weight between 200-500 Da (optimal generation range)
- Existing SAR data available for the scaffold class

**Scaffold Constraints**

MolMIM supports scaffold-constrained generation, where a substructure is preserved while peripheral groups are modified:

```python
generation_params = {
    "algorithm": "CMA-ES",
    "num_molecules": 30,
    "property_name": "QED",
    "min_similarity": 0.4,
    "particles": 30,
    "iterations": 10,
    "smi": seed_smiles  # e.g., "CC(=O)Oc1ccccc1C(=O)O"
}
response = requests.post(f"{MOLMIM_URL}/generate", json=generation_params)
```

The `min_similarity` parameter (Tanimoto on Morgan fingerprints) controls how far generated molecules can drift from the seed. Lower values (0.3) explore more chemical space; higher values (0.7) stay closer to known actives.

**Diversity Control**

To avoid redundant candidates, the pipeline applies MaxMin diversity picking after generation. From 30 raw candidates, the top 10 most structurally diverse molecules are retained for downstream evaluation.

### DiffDock Binding Prediction

DiffDock is a diffusion-based molecular docking model that predicts binding poses without requiring a pre-defined search box. The HCLS AI Factory uses the BioNeMo NIM deployment.

**Inputs and Outputs**

- **Input:** Protein structure (PDB), ligand (SDF/MOL2), number of poses
- **Output:** Ranked poses with confidence scores, per-pose RMSD estimates, contact residue lists

**Confidence Score Interpretation**

| Confidence Range | Interpretation | Action |
|-----------------|----------------|--------|
| > 0.8 | High confidence binding | Strong candidate |
| 0.5 - 0.8 | Moderate confidence | Proceed with caution |
| < 0.5 | Low confidence | Likely non-binder, discard |

**Binding Energy Conversion**

DiffDock confidence scores are converted to approximate binding energies using a calibrated regression:

```
docking_score_kcal = -2.0 + (confidence * -12.0)
```

A confidence of 0.78 yields approximately -11.4 kcal/mol, consistent with the VCP demo results.

**Contact Residue Analysis**

For each predicted pose, the pipeline extracts residues within 4.0 Angstroms of the ligand. These contact residues are cross-referenced against:

- Known active site residues from UniProt annotations
- Mutation sites from the patient's VCF
- Conserved residues across orthologs

### RDKit Chemical Property Analysis

Every generated molecule passes through a battery of RDKit property calculations before scoring.

**Lipinski's Rule of Five**

Molecules violating more than one Lipinski rule are flagged but not automatically discarded (many approved oncology drugs violate Lipinski):

| Property | Threshold | Calculation |
|----------|-----------|-------------|
| Molecular Weight | <= 500 Da | `Descriptors.MolWt(mol)` |
| LogP | <= 5.0 | `Descriptors.MolLogP(mol)` |
| H-Bond Donors | <= 5 | `Descriptors.NumHDonors(mol)` |
| H-Bond Acceptors | <= 10 | `Descriptors.NumHAcceptors(mol)` |

**Quantitative Estimate of Drug-likeness (QED)**

QED combines eight molecular descriptors into a single 0-1 score. The HCLS AI Factory uses the weighted QED variant:

```python
from rdkit.Chem.QED import qed
score = qed(mol)  # 0.0 (least drug-like) to 1.0 (most drug-like)
```

Typical thresholds: QED >= 0.5 passes initial screen; QED >= 0.7 is considered favorable.

**Synthetic Accessibility (SA) Score**

The SA score (1-10 scale, lower is easier to synthesize) helps prioritize molecules that can realistically be made:

```python
from rdkit.Chem import RDConfig
from rdkit.Contrib.SA_Score import sascorer
sa = sascorer.calculateScore(mol)  # 1.0 (trivial) to 10.0 (very hard)
```

Molecules with SA > 6.0 are flagged as synthetically challenging.

**MMFF Conformer Generation**

For 3D analysis and docking preparation, the pipeline generates low-energy conformers:

```python
from rdkit.Chem import AllChem
AllChem.EmbedMultipleConfs(mol, numConfs=50, pruneRmsThresh=0.5)
AllChem.MMFFOptimizeMoleculeConfs(mol)
```

The lowest-energy conformer is selected as the docking input for DiffDock.

### Pediatric Safety Filters

The HCLS AI Factory applies six pediatric-specific safety filters. A molecule must pass all six to be considered a viable pediatric candidate. These filters go beyond standard adult drug safety screens to address the unique vulnerabilities of developing organ systems.

**Filter 1: Blood-Brain Barrier (BBB) Penetration Risk**

- **Threshold:** MW > 500 Da triggers concern
- **Rationale:** Large molecules that cross the BBB pose neurotoxicity risk in developing brains
- **Implementation:** `Descriptors.MolWt(mol) > 500`
- **Action:** Flag, do not auto-reject (some CNS-targeted therapies require BBB penetration)

**Filter 2: Cardiac Safety (hERG Liability)**

- **Threshold:** Predicted hERG IC50 < 10 uM
- **Rationale:** QT prolongation risk is higher in pediatric patients with smaller cardiac mass
- **Implementation:** Structure-based hERG prediction using Morgan fingerprint similarity to known hERG blockers
- **Action:** Hard reject if predicted IC50 < 1 uM; flag if 1-10 uM

**Filter 3: Hepatotoxicity Risk**

- **Threshold:** LogP > 5.0
- **Rationale:** Highly lipophilic compounds accumulate in hepatic tissue; pediatric livers have immature metabolic capacity
- **Implementation:** `Descriptors.MolLogP(mol) > 5.0`
- **Action:** Hard reject

**Filter 4: Teratogenicity Screen**

- **Threshold:** SMARTS pattern match against known teratogenic substructures
- **Rationale:** Adolescent patients of reproductive age require teratogenicity screening
- **Implementation:** Library of 47 SMARTS patterns derived from FDA pregnancy category X drugs
- **Action:** Hard reject on match with known high-risk pharmacophores

**Filter 5: Oral Bioavailability**

- **Threshold:** TPSA > 140 Angstroms squared
- **Rationale:** Pediatric formulation strongly prefers oral dosing; high TPSA reduces absorption
- **Implementation:** `Descriptors.TPSA(mol) > 140`
- **Action:** Flag (injectable formulations remain possible but less preferred)

**Filter 6: Gastrointestinal Tolerability**

- **Threshold:** Rotatable bonds > 10
- **Rationale:** Flexible molecules often cause GI distress; pediatric patients are more sensitive
- **Implementation:** `Descriptors.NumRotatableBonds(mol) > 10`
- **Action:** Flag

### Composite Scoring Methodology

Candidates that survive safety filtering receive a composite score combining three weighted components:

```
composite_score = (docking_weight * normalized_docking) +
                  (generation_weight * normalized_generation) +
                  (qed_weight * normalized_qed)
```

**Default Weights:**

| Component | Weight | Source | Normalization |
|-----------|--------|--------|---------------|
| Docking Score | 0.4 | DiffDock confidence | Min-max across candidate set |
| Generation Score | 0.3 | MolMIM similarity to seed | Already 0-1 (Tanimoto) |
| QED Score | 0.3 | RDKit QED | Already 0-1 |

Weights are configurable per therapeutic area. Oncology workflows may increase docking weight to 0.5; rare disease workflows may increase generation weight to favor novelty.

### VCP Demo Results

The VCP (Valosin-Containing Protein) variant demonstration showcases the full pipeline:

- **Input Variant:** VCP p.R155H (rs121909334), associated with IBMPFD
- **Target Structure:** PDB 5KIU, ATPase domain
- **Best Candidate Docking Score:** -11.4 kcal/mol (DiffDock confidence 0.78)
- **Best Candidate QED:** 0.81
- **Candidates Generated:** 30 via MolMIM, 10 after diversity picking
- **Candidates Passing Safety:** 7 of 10 passed all 6 pediatric filters
- **Top Composite Score:** 0.84
- **Pipeline Runtime:** 127 seconds end-to-end

---

## Part 5: Cross-Agent Coordination

### Shared Genomic Evidence Architecture

All 11 intelligence agents share access to the `genomic_evidence` Milvus collection. This collection contains 35,678 curated variant-evidence vectors derived from ClinVar, AlphaMissense, and the patient's annotated VCF.

**Collection Schema:**

| Field | Type | Description |
|-------|------|-------------|
| id | INT64 (primary) | Auto-generated unique ID |
| text | VARCHAR(4096) | Evidence text (variant description + clinical significance) |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| source | VARCHAR(256) | Origin database (ClinVar, AlphaMissense, VCF) |
| gene | VARCHAR(64) | Gene symbol |
| variant_id | VARCHAR(128) | rs number or HGVS notation |
| clinical_significance | VARCHAR(128) | ClinVar classification |
| therapeutic_area | VARCHAR(128) | Mapped therapeutic area(s) |

**Access Pattern:** Read-only. Agents query `genomic_evidence` using cosine similarity search with domain-specific query expansion. Write access is restricted to the RAG ingestion pipeline (Stage 2).

**Query Expansion by Agent Type:**

Each agent type prepends domain-specific context to the raw user query before embedding:

- **Oncology Agent:** Appends "cancer somatic mutation driver gene tumor"
- **Cardiology Agent:** Appends "cardiac arrhythmia cardiomyopathy channelopathy"
- **Neurology Agent:** Appends "neurodegeneration epilepsy demyelination"
- **Pharmacogenomics Agent:** Appends "drug metabolism CYP450 pharmacokinetic"

This expansion biases the similarity search toward domain-relevant evidence without requiring separate collections.

### Cross-Modal Trigger Examples

Agents communicate through structured trigger events. When one agent's analysis produces a finding that is relevant to another agent's domain, it emits a trigger that the orchestrator routes to the appropriate downstream agent.

**Imaging Finding to Genomic Variant Query**

Scenario: The Imaging Intelligence Agent detects a suspicious mass pattern in a brain MRI.

```json
{
  "trigger_type": "imaging_to_genomic",
  "source_agent": "imaging_intelligence",
  "target_agent": "precision_oncology",
  "payload": {
    "finding": "ring-enhancing lesion, temporal lobe",
    "modality": "MRI_T1_contrast",
    "urgency": "high",
    "suggested_query": "IDH1 IDH2 ATRX TP53 glioma variants"
  }
}
```

The Oncology Agent receives this trigger and queries `genomic_evidence` for variants in IDH1, IDH2, ATRX, and TP53. Results are correlated with the imaging finding to produce an integrated assessment.

**Pharmacogenomic Alert to Dosing Recommendation**

Scenario: The Pharmacogenomics Agent identifies a CYP2D6 poor metabolizer genotype.

```json
{
  "trigger_type": "pgx_to_dosing",
  "source_agent": "pharmacogenomics_intelligence",
  "target_agent": "precision_oncology",
  "payload": {
    "gene": "CYP2D6",
    "phenotype": "poor_metabolizer",
    "star_alleles": ["*4", "*4"],
    "affected_drugs": ["tamoxifen", "codeine", "ondansetron"],
    "recommendation": "avoid_or_reduce_dose"
  }
}
```

The receiving agent adjusts therapeutic recommendations, avoiding CYP2D6-metabolized prodrugs and suggesting alternatives.

**Rare Disease Variant to Family Cascade**

Scenario: The Rare Disease Diagnostic Agent identifies a pathogenic variant in a recessive disorder gene.

```json
{
  "trigger_type": "rare_to_cascade",
  "source_agent": "rare_disease_diagnostic",
  "target_agent": "precision_biomarker",
  "payload": {
    "variant": "CFTR p.F508del",
    "zygosity": "heterozygous",
    "inheritance": "autosomal_recessive",
    "cascade_recommendation": "test_parents_and_siblings",
    "carrier_frequency": "1_in_25_European"
  }
}
```

The Biomarker Agent prepares a carrier screening panel and family counseling summary.

### Multi-Agent Query Orchestration Patterns

**Pattern 1: Fan-Out / Fan-In**

A single patient query is broadcast to multiple agents simultaneously. Each agent returns domain-specific findings, and the orchestrator merges results:

```
User Query --> Orchestrator --> [Oncology, Cardiology, PGx, Biomarker]
                                      |         |        |       |
                                      v         v        v       v
                                 Findings   Findings  Findings  Findings
                                      \        |        |       /
                                       --> Merged Report <---
```

**Pattern 2: Sequential Pipeline**

Results flow through agents in a defined order, each enriching the analysis:

```
Variant --> Biomarker Agent --> Oncology Agent --> Drug Discovery --> Clinical Trial Agent
         (classify variant)  (assess cancer risk)  (find candidates)  (match trials)
```

**Pattern 3: Conditional Routing**

The orchestrator examines variant annotations and routes to the most relevant agent:

- Variants in cardiac genes (SCN5A, KCNQ1, MYH7) route to the Cardiology Agent
- Variants in neurological genes (APP, PSEN1, HTT) route to the Neurology Agent
- Variants with oncogenic annotations route to the Oncology Agent
- Unclassified rare variants route to the Rare Disease Agent

### Building Custom Cross-Agent Workflows

To define a custom workflow, create a JSON workflow definition:

```json
{
  "workflow_name": "pediatric_cns_tumor_workup",
  "description": "Integrated CNS tumor assessment for pediatric patients",
  "steps": [
    {
      "agent": "imaging_intelligence",
      "action": "analyze_neuroimaging",
      "timeout_seconds": 30
    },
    {
      "agent": "precision_oncology",
      "action": "variant_assessment",
      "depends_on": ["imaging_intelligence"],
      "timeout_seconds": 20
    },
    {
      "agent": "pharmacogenomics_intelligence",
      "action": "drug_interaction_check",
      "parallel_with": ["precision_oncology"],
      "timeout_seconds": 15
    },
    {
      "agent": "clinical_trial_intelligence",
      "action": "match_trials",
      "depends_on": ["precision_oncology"],
      "timeout_seconds": 20
    }
  ],
  "output": "merged_report"
}
```

Register the workflow via the orchestrator API:

```bash
curl -X POST http://localhost:8080/api/workflows \
  -H "Content-Type: application/json" \
  -d @pediatric_cns_tumor_workup.json
```

---

## Part 6: Deployment & Operations

### Docker Compose Architecture

The HCLS AI Factory runs as 21 interconnected services defined in `docker-compose.dgx-spark.yml`. Services are organized into four tiers:

**Tier 1: Infrastructure (always-on)**

| Service | Image | Port | Purpose |
|---------|-------|------|---------|
| etcd | quay.io/coreos/etcd:v3.5.x | 2379 | Milvus metadata store |
| minio | minio/minio:latest | 9000, 9001 | Milvus object storage |
| milvus-standalone | milvusdb/milvus:v2.4.x | 19530, 9091 | Vector database |

**Tier 2: Core Services**

| Service | Port | Purpose |
|---------|------|---------|
| landing-page | 8080 | Flask hub with health dashboard |
| rag-api | 8888 | RAG query engine |
| genomics-portal | 8085 | Genomics web interface |
| drug-discovery-api | 8890 | Drug discovery pipeline API |

**Tier 3: Intelligence Agents**

| Service | Port | Purpose |
|---------|------|---------|
| precision-oncology-agent | 8501 | Oncology analysis |
| cart-intelligence-agent | 8502 | CAR-T therapy design |
| imaging-intelligence-agent | 8503 | Medical imaging analysis |
| precision-autoimmune-agent | 8504 | Autoimmune assessment |
| precision-biomarker-agent | 8505 | Biomarker interpretation |
| cardiology-intelligence-agent | 8506 | Cardiac variant analysis |
| pharmacogenomics-intelligence-agent | 8507 | PGx and drug interaction |
| neurology-intelligence-agent | 8508 | Neurological assessment |
| rare-disease-diagnostic-agent | 8509 | Rare disease diagnosis |
| single-cell-intelligence-agent | 8510 | Single-cell analysis |
| clinical-trial-intelligence-agent | 8511 | Clinical trial matching |

**Tier 4: Monitoring**

| Service | Port | Purpose |
|---------|------|---------|
| prometheus | 9090 | Metrics collection |
| grafana | 3000 | Dashboards and alerting |

### Health Monitoring

The `health-monitor.sh` script (19 KB) provides comprehensive health monitoring across 22 service endpoints. It runs via cron every 5 minutes.

**Monitored Checks per Service:**

1. TCP port reachability (timeout: 5 seconds)
2. HTTP health endpoint response (expected: 200 OK)
3. Process memory consumption (threshold: 90% of allocated)
4. Response latency (threshold: 10 seconds)

**Auto-Restart Behavior:**

When a service fails 3 consecutive health checks, the monitor executes:

```bash
# For Docker services
docker restart <service_name>

# For native processes
kill -TERM <pid> && sleep 2 && <start_command>
```

**Watchdog Mode:**

The health monitor itself is supervised by a systemd watchdog. If `health-monitor.sh` fails to write a heartbeat file within 10 minutes, systemd restarts it.

**Log Locations:**

- Health check results: `logs/cron-health.log`
- Service-specific logs: `logs/<service-name>.log`
- Auto-restart events: `logs/health-restarts.log`

### Milvus Administration

The platform manages 139 Milvus collections across all agents and core services.

**Collection Categories:**

| Category | Count | Example |
|----------|-------|---------|
| Genomic evidence | 1 | `genomic_evidence` (35,678 vectors) |
| Agent knowledge bases | 55 | `oncology_guidelines`, `cart_protocols` |
| ClinVar annotations | 12 | `clinvar_pathogenic`, `clinvar_vus` |
| Drug/compound data | 8 | `drug_interactions`, `chembl_targets` |
| Clinical trials | 14 | `trial_eligibility`, `trial_endpoints` |
| Literature | 22 | `pubmed_oncology`, `pubmed_cardiology` |
| Patient data (demo) | 5 | `demo_patient_variants` |
| Imaging features | 10 | `imaging_features_mri`, `imaging_features_ct` |
| Other | 12 | `system_prompts`, `workflow_templates` |

**Backup Strategy:**

```bash
# Full backup of all collections (run weekly)
python3 -c "
from pymilvus import connections, utility
connections.connect('default', host='localhost', port='19530')
collections = utility.list_collections()
for coll in collections:
    utility.do_bulk_insert(coll, files=[], backup=True)
"

# Incremental backup via Milvus backup tool
milvus-backup create --name weekly_$(date +%Y%m%d)
```

**Flush Rate Tuning:**

By default, Milvus flushes segments every 1 second. For bulk ingestion (e.g., loading 3.56M annotated variants), increase the flush interval:

```python
from pymilvus import Collection
collection = Collection("genomic_evidence")
collection.set_properties({"collection.flush.interval": "30"})  # 30 seconds
```

Reset to default after ingestion completes.

### Agent Startup Sequence

Each intelligence agent follows a deterministic startup sequence:

1. **Configuration Load** (< 1 s): Read environment variables and config files
2. **Embedding Model Load** (5-15 s): Load BGE-small-en-v1.5 into GPU memory (~130 MB)
3. **Milvus Connection** (1-3 s): Connect to Milvus standalone, verify collection availability
4. **Collection Warm-up** (2-5 s): Load frequently-queried collections into memory
5. **Health Endpoint Registration** (< 1 s): Register with the health monitor
6. **API Ready** (< 1 s): Begin accepting requests

**Total cold start time:** 10-25 seconds per agent. When running all 11 agents, stagger starts by 3 seconds each to avoid GPU memory contention during embedding model loading.

### Performance Tuning

**GPU Memory Management**

The DGX Spark provides 128 GB unified memory. Recommended allocation:

| Component | Memory | Notes |
|-----------|--------|-------|
| Milvus indexes | 16-24 GB | Depends on loaded collections |
| BGE-small embedding (x11) | ~1.4 GB | 130 MB per agent instance |
| MolMIM NIM | 4-8 GB | Active during drug discovery |
| DiffDock NIM | 4-8 GB | Active during drug discovery |
| Parabricks | 32-64 GB | Active during genomics pipeline |
| System / OS | 8-16 GB | Linux kernel, Docker overhead |

**Concurrent Agent Tuning**

For systems with limited memory, reduce the number of concurrent agents:

```yaml
# In docker-compose.dgx-spark.yml, add resource limits
services:
  precision-oncology-agent:
    deploy:
      resources:
        limits:
          memory: 4G
        reservations:
          memory: 2G
```

**Milvus Index Tuning**

The default IVF_FLAT index uses `nlist=1024` and `nprobe=16`. Adjust for your workload:

| Scenario | nlist | nprobe | Trade-off |
|----------|-------|--------|-----------|
| Low latency (< 10 ms) | 1024 | 10 | Slightly lower recall |
| Balanced (default) | 1024 | 16 | Good recall, ~15 ms |
| High recall | 2048 | 32 | Higher latency (~30 ms) |
| Maximum recall | 4096 | 64 | ~50 ms, near-exhaustive |

```python
search_params = {"metric_type": "COSINE", "params": {"nprobe": 32}}
results = collection.search(
    data=[query_embedding],
    anns_field="embedding",
    param=search_params,
    limit=20
)
```

### HTTPS via Caddy

For production deployments, use Caddy as a reverse proxy with automatic TLS:

```
# Caddyfile
hcls-ai-factory.example.com {
    handle /api/rag/* {
        reverse_proxy localhost:8888
    }
    handle /api/genomics/* {
        reverse_proxy localhost:8085
    }
    handle /api/drugs/* {
        reverse_proxy localhost:8890
    }
    handle /agents/* {
        reverse_proxy localhost:8501-8511
    }
    handle {
        reverse_proxy localhost:8080
    }
    header {
        X-Content-Type-Options nosniff
        X-Frame-Options DENY
        Referrer-Policy strict-origin-when-cross-origin
    }
}
```

### Security

**Portal Authentication**

The landing page (`localhost:8080`) supports configurable authentication:

```python
# In landing-page/app.py
PORTAL_AUTH_ENABLED = os.getenv("PORTAL_AUTH_ENABLED", "false").lower() == "true"
PORTAL_USERNAME = os.getenv("PORTAL_USERNAME", "admin")
PORTAL_PASSWORD_HASH = os.getenv("PORTAL_PASSWORD_HASH")  # bcrypt hash
```

**XSS Protection**

All agent Streamlit interfaces and the Flask landing page implement output sanitization:

- HTML entities escaped in all rendered text
- Content-Security-Policy headers set to `default-src 'self'`
- User inputs validated and sanitized before Milvus query construction

**Pickle Hardening**

The platform avoids `pickle.loads()` on untrusted data. Model artifacts use SafeTensors format. Where pickle is unavoidable (legacy sklearn models), a restricted unpickler is used:

```python
import pickle
import io

ALLOWED_CLASSES = {"numpy.ndarray", "sklearn.ensemble._forest.RandomForestClassifier"}

class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        full_name = f"{module}.{name}"
        if full_name not in ALLOWED_CLASSES:
            raise pickle.UnpicklingError(f"Blocked: {full_name}")
        return super().find_class(module, name)

def safe_load(data: bytes):
    return RestrictedUnpickler(io.BytesIO(data)).load()
```

---

## Part 7: Extending the Platform

### Adding a New Intelligence Agent

Follow these steps to add a custom intelligence agent (e.g., a "Reproductive Health Agent"):

**Step 1: Create Milvus Collections**

Define the knowledge collections your agent needs:

```python
from pymilvus import Collection, CollectionSchema, FieldSchema, DataType, connections

connections.connect("default", host="localhost", port="19530")

fields = [
    FieldSchema(name="id", dtype=DataType.INT64, is_primary=True, auto_id=True),
    FieldSchema(name="text", dtype=DataType.VARCHAR, max_length=4096),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=384),
    FieldSchema(name="source", dtype=DataType.VARCHAR, max_length=256),
    FieldSchema(name="category", dtype=DataType.VARCHAR, max_length=128),
]

schema = CollectionSchema(fields, description="Reproductive health knowledge base")
collection = Collection("reproductive_health_knowledge", schema)

# Create index
index_params = {"index_type": "IVF_FLAT", "metric_type": "COSINE", "params": {"nlist": 1024}}
collection.create_index("embedding", index_params)
```

**Step 2: Ingest Domain Knowledge**

Prepare knowledge sources and embed them:

```python
from sentence_transformers import SentenceTransformer

model = SentenceTransformer("BAAI/bge-small-en-v1.5")

documents = [
    {"text": "Preeclampsia affects 2-8% of pregnancies...", "source": "ACOG", "category": "guidelines"},
    # ... more documents
]

for doc in documents:
    embedding = model.encode(doc["text"]).tolist()
    collection.insert([[doc["text"]], [embedding], [doc["source"]], [doc["category"]]])

collection.flush()
```

**Step 3: Create the Agent API**

Use the standard agent template (Streamlit + FastAPI):

```
reproductive-health-agent/
    app.py              # Streamlit UI
    api.py              # FastAPI health + query endpoints
    knowledge/          # Ingestion scripts and raw data
    requirements.txt
    Dockerfile
```

The `api.py` must expose:

- `GET /health` -- returns `{"status": "healthy", "collections": [...], "uptime": ...}`
- `POST /query` -- accepts `{"query": "...", "patient_context": {...}}` and returns findings

**Step 4: Register in the UI**

Add the agent to the landing page configuration:

```python
# In landing-page/config.py
AGENTS = {
    # ... existing agents ...
    "reproductive_health": {
        "name": "Reproductive Health Agent",
        "port": 8512,
        "icon": "baby",
        "description": "Maternal-fetal medicine and reproductive genomics",
        "collections": ["reproductive_health_knowledge", "genomic_evidence"]
    }
}
```

**Step 5: Add to Docker Compose**

```yaml
reproductive-health-agent:
  build: ./reproductive-health-agent
  ports:
    - "8512:8512"
  environment:
    - MILVUS_HOST=milvus-standalone
    - MILVUS_PORT=19530
    - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
  depends_on:
    - milvus-standalone
  restart: unless-stopped
```

**Step 6: Write Tests**

At minimum, every agent needs:

- Unit tests for query processing and response formatting
- Integration test confirming Milvus connectivity and collection access
- End-to-end test with a sample query and expected output structure

```python
def test_reproductive_health_query():
    response = client.post("/query", json={"query": "BRCA1 and pregnancy risk"})
    assert response.status_code == 200
    assert "findings" in response.json()
    assert len(response.json()["findings"]) > 0
```

### Creating Custom Milvus Collections with Dynamic Fields

Milvus supports dynamic fields for flexible schema evolution:

```python
schema = CollectionSchema(fields, enable_dynamic_field=True)
collection = Collection("flexible_collection", schema)

# Insert with extra fields not in the schema
collection.insert([{
    "text": "Example document",
    "embedding": [0.1] * 384,
    "custom_score": 0.95,       # dynamic field
    "reviewer": "Dr. Smith"     # dynamic field
}])
```

Dynamic fields are stored as JSON and can be filtered in search expressions:

```python
results = collection.search(
    data=[query_embedding],
    anns_field="embedding",
    param=search_params,
    limit=10,
    expr='custom_score > 0.8 and reviewer == "Dr. Smith"'
)
```

### Integrating New Knowledge Sources

**PubMed Integration**

Use the NCBI E-utilities API to fetch and ingest PubMed abstracts:

```python
import requests

def fetch_pubmed_abstracts(query: str, max_results: int = 100):
    # Search for PMIDs
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_resp = requests.get(search_url, params={
        "db": "pubmed", "term": query, "retmax": max_results, "retmode": "json"
    })
    pmids = search_resp.json()["esearchresult"]["idlist"]

    # Fetch abstracts
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_resp = requests.get(fetch_url, params={
        "db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml"
    })
    return parse_abstracts(fetch_resp.text)
```

**ClinicalTrials.gov Integration**

The ClinicalTrials.gov API v2 provides structured trial data:

```python
def fetch_trials(condition: str, status: str = "RECRUITING"):
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": condition,
        "filter.overallStatus": status,
        "pageSize": 50,
        "fields": "NCTId,BriefTitle,EligibilityCriteria,InterventionName"
    }
    response = requests.get(url, params=params)
    return response.json()["studies"]
```

Ingest trial data into agent-specific collections for real-time trial matching.

### Custom Clinical Workflows

Beyond the built-in workflows, you can define condition-specific analysis pipelines:

```python
from hcls_common.workflow import WorkflowEngine

workflow = WorkflowEngine()

@workflow.step("variant_triage")
def triage_variants(context):
    """Filter and prioritize variants by clinical significance."""
    variants = context["variants"]
    return [v for v in variants if v["clinical_significance"] in ("Pathogenic", "Likely_pathogenic")]

@workflow.step("literature_search")
def search_literature(context):
    """Query PubMed and internal knowledge bases for supporting evidence."""
    genes = set(v["gene"] for v in context["triaged_variants"])
    evidence = {}
    for gene in genes:
        evidence[gene] = query_milvus(f"{gene} clinical significance treatment")
    return evidence

@workflow.step("generate_report")
def generate_report(context):
    """Produce a clinical summary using Claude."""
    prompt = build_clinical_prompt(context["triaged_variants"], context["evidence"])
    return call_claude(prompt)

# Execute
result = workflow.run(
    steps=["variant_triage", "literature_search", "generate_report"],
    initial_context={"variants": patient_variants}
)
```

### Contributing to the Open-Source Project

The HCLS AI Factory is released under the Apache 2.0 License. Contributions are welcome.

**Getting Started:**

1. Fork the repository at `https://github.com/ajones1923/hcls-ai-factory`
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Follow the existing code style (Black formatting, type hints, docstrings)
4. Add tests for any new functionality
5. Run the test suite: `pytest tests/ -v`
6. Submit a pull request with a clear description of changes

**Areas Seeking Contributions:**

- New intelligence agents for underserved therapeutic areas
- Improved embedding models for biomedical text
- Additional safety filters for specific patient populations
- Performance optimizations for large-scale variant sets
- Documentation and tutorials
- Internationalization of clinical content

**Code Review Process:**

All pull requests require:

- Passing CI checks (linting, tests, type checking)
- At least one maintainer review
- No reduction in test coverage
- Updated documentation if user-facing behavior changes

---

!!! warning "Clinical Decision Support Disclaimer"
    The HCLS AI Factory platform and its components are clinical decision support research tools. They are not FDA-cleared and are not intended as standalone diagnostic devices. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
