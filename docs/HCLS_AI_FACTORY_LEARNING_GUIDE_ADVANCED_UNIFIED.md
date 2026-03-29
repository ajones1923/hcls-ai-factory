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

The Precision Oncology Agent is a closed-loop clinical decision support system designed for molecular tumor board (MTB) workflows. Given a patient's somatic and germline VCF, it identifies actionable variants, matches them against oncology knowledge bases, ranks candidate therapies, surfaces relevant clinical trials, and assembles a structured MTB-ready report. The entire pipeline runs on a single NVIDIA DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_AGENT_PORT` | `8527` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size; increase for concurrent MTB sessions |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model ID; upgrade to `claude-opus-4-20250514` for complex cases |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | Sentence-transformer for 384-dim embeddings |
| `TOP_K` | `15` | Max chunks retrieved per collection search |
| `EVIDENCE_THRESHOLD` | `0.65` | Minimum cosine similarity for evidence inclusion |
| `CROSS_AGENT_TIMEOUT` | `30` | Seconds to wait for cross-agent responses |

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with component status |
| GET | `/collections` | List all loaded collections with record counts |
| GET | `/knowledge/stats` | Knowledge base statistics (total vectors, last ingest) |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/query` | Free-text oncology RAG query with evidence synthesis |
| POST | `/search` | Direct vector similarity search across collections |
| POST | `/find-related` | Find related concepts given an entity (gene, drug, pathway) |
| POST | `/v1/onco/integrated-assessment` | Full multi-collection integrated assessment |
| POST | `/api/ask` | Meta-agent question answering with chain-of-thought |
| POST | `/api/cases` | Create a new patient case with cancer type and variants |
| GET | `/api/cases/{case_id}` | Retrieve case details and analysis status |
| POST | `/api/cases/{case_id}/mtb` | Run full molecular tumor board analysis for a case |
| GET | `/api/cases/{case_id}/variants` | List annotated variants for a case |
| POST | `/api/trials/match` | Match molecular profile to clinical trials |
| POST | `/api/trials/match-case/{case_id}` | Match a stored case to clinical trials |
| POST | `/api/therapies/rank` | Rank therapies by evidence level and biomarker fit |
| POST | `/api/reports/generate` | Generate MTB report (Markdown, JSON, PDF, FHIR R4) |
| GET | `/api/reports/{case_id}/{fmt}` | Download generated report in specified format |
| GET | `/api/events` | List cross-agent events (SSE stream) |
| GET | `/api/events/{event_id}` | Retrieve a specific event by ID |

**Detailed Clinical Workflow Walkthrough**

1. **Variant submission** -- Clinician uploads somatic/germline VCF via Streamlit UI or POSTs variant list to `/api/cases`. The VCF parser extracts gene, protein change, allele frequency, and quality metrics.
2. **Entity detection** -- The agent identifies gene names (EGFR, BRAF, ALK), variant notations (L858R, V600E), and cancer types from free-text or structured input.
3. **Evidence-level annotation** -- Variants are annotated against CIViC and OncoKB with AMP/ASCO/CAP evidence-level tiering (Level IA through Level III). Each variant receives a clinical significance classification.
4. **Multi-collection retrieval** -- Parallel search across all 11 collections: therapies, biomarkers, trials, resistance mechanisms, pathways, guidelines, outcomes, and cases. Each result carries a cosine similarity score.
5. **TherapyRanker scoring** -- Approved and investigational therapies are scored using a composite of evidence level (40%), biomarker match (25%), resistance profile (15%), guideline concordance (10%), and trial availability (10%).
6. **TrialMatcher execution** -- Molecular profile is matched against eligibility criteria in `onco_trials`. Results include NCT IDs, phase, enrollment status, and distance to nearest trial sites.
7. **Claude LLM synthesis** -- All retrieved evidence is assembled into a structured prompt. Claude generates a narrative summary with inline citations, therapy recommendations, and resistance warnings.
8. **MTB report generation** -- Structured report exported as Markdown, JSON, PDF, or FHIR R4 DiagnosticReport with tiered evidence, therapy rankings, trial matches, and resistance alerts.

**Example Query and Response**

```bash
curl -X POST http://localhost:8527/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What targeted therapies are available for EGFR L858R in NSCLC?"}'
```

```json
{
  "answer": "EGFR L858R is a Level IA actionable target in NSCLC. Osimertinib (3rd-gen TKI) is the preferred first-line therapy per NCCN 2026 guidelines, with a median PFS of 18.9 months (FLAURA trial). Second-line options include amivantamab + lazertinib for C797S-mediated resistance...",
  "evidence": [
    {"collection": "onco_therapies", "score": 0.94, "text": "Osimertinib: 3rd-generation EGFR TKI, FDA-approved for EGFR-mutant NSCLC..."},
    {"collection": "onco_guidelines", "score": 0.91, "text": "NCCN NSCLC v3.2026: EGFR L858R — Preferred first-line: osimertinib..."},
    {"collection": "onco_resistance", "score": 0.87, "text": "C797S mutation confers resistance to osimertinib in 10-15% of cases..."},
    {"collection": "onco_trials", "score": 0.83, "text": "NCT04487080: Phase III, osimertinib + savolitinib for MET-amplified progression..."}
  ],
  "therapies_ranked": [
    {"rank": 1, "drug": "Osimertinib", "evidence": "Level IA", "score": 0.96},
    {"rank": 2, "drug": "Amivantamab + Lazertinib", "evidence": "Level IB", "score": 0.82},
    {"rank": 3, "drug": "Erlotinib + Ramucirumab", "evidence": "Level IB", "score": 0.74}
  ],
  "active_trials": 8,
  "confidence": 0.93,
  "collections_searched": 11,
  "processing_time_ms": 3420
}
```

**Cross-Agent Integration**

- **Receives** genomic evidence from Stage 1 pipeline via shared `genomic_evidence` collection (3.56M variants)
- **Triggers** Clinical Trial Agent (:8538) when novel actionable variants are identified, passing variant list and cancer type for trial matching
- **Receives** cardiotoxicity alerts from Cardiology Agent (:8126) for proposed chemotherapy regimens (anthracyclines, trastuzumab) with pre-treatment ejection fraction data
- **Sends** antigen expression queries to Single-Cell Agent (:8540) for immunotherapy candidate validation (TME profiling, PD-L1 expression)
- **Connects** to Drug Discovery pipeline (Stage 3) for candidate molecule docking when no approved therapies match the variant profile

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Single collection search (top-15) | 80-150 ms |
| Full 11-collection parallel search | 200-400 ms |
| Embedding generation (BGE-small) | 15-25 ms |
| Claude LLM synthesis | 2,000-3,500 ms |
| Full MTB analysis (end-to-end) | 3,500-5,000 ms |
| Report generation (PDF) | 1,500-2,500 ms |

**Common Pitfalls and Troubleshooting**

- **"No actionable variants found"** -- Check that the VCF contains PASS-filtered variants with gene annotations. Raw VCF without functional annotation will not match CIViC/OncoKB entries. Run `bcftools view -f PASS` first.
- **Therapy ranking returns empty** -- Verify that `onco_therapies` collection is loaded (`GET /collections`). If the collection shows 0 records, re-run the ingest script.
- **Slow LLM responses (>10s)** -- Reduce `TOP_K` from 15 to 8 to shrink the context window. Alternatively, switch to a faster Claude model via `LLM_MODEL`.
- **Cross-agent timeout** -- If the Cardiology or Trial agent is down, the oncology agent logs a warning but continues without cross-agent enrichment. Check `CROSS_AGENT_TIMEOUT` and verify target agents are healthy.
- **Duplicate case IDs** -- Case creation is idempotent on `patient_id` + `cancer_type`. Submitting the same combination returns the existing case rather than creating a duplicate.

---

### 3.2 CAR-T Intelligence Agent (:8522)

**Overview and Clinical Purpose**

The CAR-T Intelligence Agent breaks down data silos across the 5 stages of CAR-T cell therapy development: target selection, construct design, manufacturing, clinical evaluation, and post-market surveillance. It features a unique comparative analysis mode that auto-detects "X vs Y" queries and produces structured side-by-side comparisons. The agent covers all 6 FDA-approved CAR-T products, 25 target antigens, and 39+ product aliases, running on a single DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `CART_AGENT_PORT` | `8522` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model; use `claude-opus-4-20250514` for nuanced comparative analyses |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | 384-dim sentence-transformer |
| `TOP_K` | `15` | Max chunks per collection search |
| `COMPARATIVE_THRESHOLD` | `0.70` | Minimum confidence to trigger comparative mode |
| `QUERY_EXPANSION_MAPS` | `12` | Number of domain-specific expansion maps (169 keywords to 1,496 terms) |

**Milvus Collections (11)**

| Collection | Records | Source |
|-----------|---------|--------|
| `cart_literature` | 5,047 | PubMed abstracts via NCBI E-utilities |
| `cart_trials` | 973 | ClinicalTrials.gov API v2 |
| `cart_constructs` | 6 | 6 FDA-approved CAR-T products |
| `cart_assay_results` | 45 | Curated from landmark papers (ELIANA, ZUMA-1, KarMMa, CARTITUDE-1) |
| `cart_manufacturing` | 30 | CMC/process data (transduction, expansion, release, cryo, logistics) |
| `cart_safety` | -- | Pharmacovigilance, CRS/ICANS profiles |
| `cart_biomarkers` | -- | CRS prediction, exhaustion monitoring |
| `cart_regulatory` | -- | Approval timelines, post-marketing requirements |
| `cart_sequences` | -- | Molecular binding, scFv sequences |
| `cart_rwe` | -- | Registry outcomes, real-world data |
| `genomic_evidence` | 3.56M | Shared from Stage 2 RAG pipeline (read-only) |

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus and engine status |
| GET | `/collections` | List loaded collections with record counts |
| GET | `/knowledge/stats` | Knowledge base statistics and last ingest time |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/query` | Cross-functional RAG query (standard or comparative, auto-detected) |
| POST | `/search` | Direct vector similarity search with collection filtering |
| POST | `/find-related` | Find related entities (antigens, products, mechanisms) |
| POST | `/v1/cart/integrated-assessment` | Full multi-collection integrated CAR-T assessment |
| POST | `/api/ask` | Meta-agent question answering with reasoning chain |
| GET | `/api/reports/{patient_id}` | Retrieve generated report for a patient |
| GET | `/api/reports/{patient_id}/{fmt}` | Download report in specified format (md, json, pdf, fhir) |
| GET | `/api/events` | List cross-agent events |
| GET | `/api/events/{event_id}` | Retrieve a specific cross-agent event |

**Detailed Clinical Workflow Walkthrough**

1. **Query submission** -- User submits a question about any stage of CAR-T development via Streamlit UI or POST to `/query`. Input can be free-text or structured with explicit parameters.
2. **Comparative detection** -- The engine scans for "vs", "versus", "compare", "difference between" patterns. If detected, it enters comparative mode; otherwise, standard RAG mode.
3. **Entity resolution (comparative)** -- Both entities are parsed and resolved against the knowledge graph containing 25 antigens, 6 FDA-approved products, and 39+ aliases (e.g., "tisa-cel" resolves to "tisagenlecleucel / Kymriah").
4. **Dual retrieval (comparative)** -- Separate searches are run for each entity across all 11 collections. Results are aligned by category (efficacy, safety, manufacturing, cost) for side-by-side comparison.
5. **Query expansion (standard)** -- 12 domain-specific expansion maps transform 169 seed keywords into 1,496 search terms. For example, "CRS management" expands to include "cytokine release syndrome", "tocilizumab", "IL-6 blockade", "grading criteria".
6. **Parallel multi-collection search** -- All 11 collections are searched concurrently. Results are scored by cosine similarity and weighted by collection relevance to the query domain.
7. **Claude LLM synthesis** -- Retrieved evidence is assembled with comparative formatting when applicable. Responses include clickable PubMed IDs (PMID) and ClinicalTrials.gov NCT numbers.
8. **Report generation** -- Exportable as Markdown, JSON, PDF, or FHIR R4 DiagnosticReport.

**Example Query and Response**

```bash
curl -X POST http://localhost:8522/query \
  -H "Content-Type: application/json" \
  -d '{"question": "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"}'
```

```json
{
  "answer": "4-1BB (CD137) and CD28 costimulatory domains represent the two main signaling architectures in FDA-approved CAR-T products for DLBCL...",
  "comparison": {
    "entity_a": "4-1BB (CD137)",
    "entity_b": "CD28",
    "dimensions": [
      {"category": "Persistence", "entity_a": "Superior T-cell persistence (>6 months detectable)", "entity_b": "Rapid expansion, shorter persistence (~3 months)"},
      {"category": "CRS Rate", "entity_a": "Lower grade 3+ CRS (2-22%)", "entity_b": "Higher grade 3+ CRS (13-15%)"},
      {"category": "FDA Products", "entity_a": "Tisagenlecleucel (Kymriah), Lisocabtagene (Breyanzi)", "entity_b": "Axicabtagene (Yescarta), Brexucabtagene (Tecartus)"},
      {"category": "ORR (DLBCL)", "entity_a": "52-73%", "entity_b": "72-83%"}
    ]
  },
  "evidence": [
    {"collection": "cart_constructs", "score": 0.93, "text": "4-1BB signaling promotes memory T-cell differentiation..."},
    {"collection": "cart_literature", "score": 0.89, "text": "ZUMA-1: axi-cel (CD28) achieved 83% ORR in r/r DLBCL (PMID: 28982775)..."},
    {"collection": "cart_safety", "score": 0.86, "text": "CRS incidence comparison across costimulatory domains..."}
  ],
  "mode": "comparative",
  "confidence": 0.91,
  "processing_time_ms": 3180
}
```

**Cross-Agent Integration**

- **Receives** antigen escape alerts from Single-Cell Agent (:8540) when off-tumor expression is detected for a target antigen (e.g., CD19 loss in B-ALL relapse)
- **Queries** Oncology Agent (:8527) for tumor mutational burden and microsatellite instability context that may influence CAR-T eligibility
- **Shares** CRS/ICANS toxicity profiles with Biomarker Agent (:8529) for real-time monitoring protocols (IL-6, ferritin, CRP thresholds)
- **Triggers** Clinical Trial Agent (:8538) when a query involves investigational CAR-T constructs not yet FDA-approved

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Standard single-collection search | 80-150 ms |
| Full 11-collection parallel search | 200-400 ms |
| Comparative dual-entity retrieval | 350-600 ms |
| Query expansion (12 maps) | 5-10 ms |
| Claude LLM synthesis (standard) | 2,000-3,000 ms |
| Claude LLM synthesis (comparative) | 2,500-4,000 ms |
| End-to-end query (standard) | 2,500-3,500 ms |
| End-to-end query (comparative) | 3,000-5,000 ms |

**Common Pitfalls and Troubleshooting**

- **Comparative mode not triggering** -- The detection engine requires explicit comparison language ("vs", "versus", "compare", "difference between"). Rephrase ambiguous queries to include one of these trigger words.
- **Entity resolution failure** -- If a product name is misspelled or uses a non-standard alias, the knowledge graph may not resolve it. Use canonical names (e.g., "axicabtagene ciloleucel" not "axi-cel") or add custom aliases via the configuration.
- **Empty `cart_constructs` collection** -- This is a small curated collection (6 records). If it returns 0, re-run `scripts/seed_constructs.py` to reload the 6 FDA-approved product profiles.
- **Slow comparative queries** -- Comparative mode runs two full retrieval passes. If latency exceeds 5s, reduce `TOP_K` from 15 to 10 or limit the collection set.
- **Missing PubMed citations** -- The `cart_literature` collection requires periodic refresh via the NCBI E-utilities ingest script. Stale data will miss recent publications.

---

### 3.3 Precision Biomarker Agent (:8529)

**Overview and Clinical Purpose**

The Precision Biomarker Agent interprets patient biomarker panels with genotype awareness. It estimates biological age using PhenoAge/GrimAge algorithms, detects disease trajectories across 6 categories, provides pharmacogenomic profiling for 7 key pharmacogenes, and generates genotype-adjusted reference ranges. All processing runs on a single DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `BIOMARKER_AGENT_PORT` | `8529` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model ID |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | 384-dim sentence-transformer |
| `TOP_K` | `15` | Max chunks per collection search |
| `AGING_ALGORITHM` | `phenoage` | Biological age algorithm: `phenoage` or `grimage` |
| `DISEASE_CATEGORIES` | `6` | Number of disease trajectory categories screened |
| `PGX_GENES` | `7` | Number of pharmacogenes profiled (CYP2D6, CYP2C19, CYP2C9, CYP3A5, SLCO1B1, VKORC1, MTHFR) |

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus and engine status |
| GET | `/collections` | List loaded collections with record counts |
| GET | `/knowledge/stats` | Knowledge base statistics |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/v1/biomarker/integrated-assessment` | Full multi-collection integrated biomarker assessment |
| POST | `/v1/analyze` | Full patient analysis (biomarkers + genotypes combined) |
| POST | `/v1/biological-age` | Biological age calculation (PhenoAge/GrimAge) |
| POST | `/v1/disease-risk` | Disease risk trajectory analysis across 6 categories |
| POST | `/v1/pgx` | Pharmacogenomic profile from star allele diplotypes |
| POST | `/v1/query` | Free-text RAG query across all biomarker collections |
| POST | `/v1/query/stream` | Streaming RAG query (Server-Sent Events) |
| POST | `/v1/report/generate` | Generate clinical report (Markdown, JSON, PDF, FHIR R4) |
| GET | `/v1/report/{report_id}/pdf` | Download generated PDF report |
| POST | `/v1/report/fhir` | Export analysis as FHIR R4 DiagnosticReport |
| POST | `/v1/events/cross-modal` | Receive cross-modal events from other agents |
| POST | `/v1/events/biomarker-alert` | Receive biomarker alert events (threshold breaches) |
| GET | `/v1/events/cross-modal` | List received cross-modal events |
| GET | `/v1/events/biomarker-alert` | List received biomarker alerts |

**Detailed Clinical Workflow Walkthrough**

1. **Data submission** -- Patient biomarker panel (CRP, HbA1c, LDL, albumin, creatinine, glucose, etc.) and genotype data (star alleles or VCF-derived diplotypes) are submitted via Streamlit UI or POST to `/v1/analyze`.
2. **Biological Age Engine** -- Computes PhenoAge and GrimAge estimates from 9 routine blood biomarkers (albumin, creatinine, glucose, CRP, lymphocyte %, mean cell volume, red cell distribution width, alkaline phosphatase, white blood cell count). Output: chronological age, biological age, acceleration (positive = aging faster).
3. **Disease Trajectory Analyzer** -- Screens 6 disease categories (cardiovascular, metabolic/diabetes, hepatic, renal, inflammatory, hematologic) using genotype-stratified thresholds. Each trajectory returns a stage (normal, borderline, early, established, advanced) with progression rate.
4. **Pharmacogenomic Mapper** -- Interprets star alleles for CYP2D6, CYP2C19, CYP2C9, CYP3A5, SLCO1B1, VKORC1, MTHFR, plus HLA-B*57:01 screening. Maps diplotypes to metabolizer phenotypes (ultra-rapid, normal, intermediate, poor) with CPIC-level dosing recommendations.
5. **Genotype Adjustment Engine** -- Modifies standard reference ranges based on genomic context. For example: PNPLA3 I148M carriers have lower ALT thresholds for NAFLD screening; TCF7L2 rs7903146 T/T carriers have elevated fasting glucose baselines; APOE e4/e4 carriers have higher LDL reference ceilings.
6. **Evidence retrieval** -- Parallel search across all 11 collections to ground recommendations in published literature, CPIC guidelines, and clinical evidence.
7. **Report generation** -- 12-section clinical report covering biological age summary, disease trajectories (6 categories), pharmacogenomic profile, genotype-adjusted ranges, nutrition recommendations, drug interaction alerts, and monitoring schedule. Exported as PDF and FHIR R4.

**Example Query and Response**

```bash
curl -X POST http://localhost:8529/v1/analyze \
  -H "Content-Type: application/json" \
  -d '{"biomarkers": {"CRP": 2.1, "HbA1c": 6.8, "LDL": 145, "albumin": 4.2, "creatinine": 0.9, "glucose": 118}, "genotypes": {"CYP2D6": "*1/*4", "CYP2C19": "*1/*2"}}'
```

```json
{
  "biological_age": {
    "chronological_age": 55,
    "phenoage": 58.3,
    "grimage": 57.1,
    "acceleration": "+3.3 years",
    "percentile": 72,
    "interpretation": "Moderately accelerated aging. CRP and glucose are primary contributors."
  },
  "disease_trajectories": [
    {"category": "Metabolic", "stage": "Borderline", "risk_score": 0.62, "key_markers": ["HbA1c: 6.8 (pre-diabetic)", "glucose: 118"]},
    {"category": "Cardiovascular", "stage": "Early", "risk_score": 0.48, "key_markers": ["LDL: 145", "CRP: 2.1"]},
    {"category": "Inflammatory", "stage": "Borderline", "risk_score": 0.35, "key_markers": ["CRP: 2.1"]}
  ],
  "pharmacogenomics": {
    "CYP2D6": {"diplotype": "*1/*4", "phenotype": "Intermediate Metabolizer", "drugs_affected": ["codeine", "tramadol", "tamoxifen"]},
    "CYP2C19": {"diplotype": "*1/*2", "phenotype": "Intermediate Metabolizer", "drugs_affected": ["clopidogrel", "omeprazole", "escitalopram"]}
  },
  "genotype_adjustments": [
    {"marker": "CYP2D6", "adjustment": "Codeine: avoid or reduce dose by 50%; consider morphine alternative"}
  ],
  "confidence": 0.88,
  "processing_time_ms": 3850
}
```

**Cross-Agent Integration**

- **Receives** inflammation monitoring requests from Autoimmune Agent (:8532) during flare events, tracking CRP, ESR, IL-6, and calprotectin trends
- **Provides** biological age context to Cardiology Agent (:8126) for ASCVD risk refinement -- accelerated biological age increases the 10-year risk estimate
- **Shares** pharmacogenomic profiles with Pharmacogenomics Agent (:8107) for detailed CPIC/DPWG dosing guidance
- **Sends** biomarker alerts to Oncology Agent (:8527) when tumor markers (CEA, CA-125, PSA) exceed surveillance thresholds
- **Receives** genotype-adjusted reference range requests from Rare Disease Agent (:8134) for rare metabolic conditions

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Single collection search (top-15) | 80-150 ms |
| Full 11-collection parallel search | 200-400 ms |
| Biological age calculation | 50-100 ms |
| Disease trajectory analysis (6 categories) | 100-200 ms |
| PGx diplotype interpretation | 30-60 ms |
| Claude LLM synthesis | 2,000-3,500 ms |
| Full analysis (end-to-end) | 3,000-5,000 ms |
| PDF report generation | 1,500-2,500 ms |

**Common Pitfalls and Troubleshooting**

- **Biological age returns null** -- The PhenoAge algorithm requires at least 7 of the 9 input biomarkers. If fewer are provided, the engine falls back to a partial estimate with a warning flag. Ensure albumin, creatinine, and glucose are included at minimum.
- **PGx phenotype "Indeterminate"** -- This occurs when the star allele combination is not in the CPIC lookup table. Verify that diplotype notation follows PharmVar conventions (e.g., `*1/*4` not `*1/4`).
- **Genotype adjustments not applied** -- The adjustment engine requires both the biomarker value and the relevant genotype. Submitting biomarkers without genotypes produces standard (non-adjusted) reference ranges.
- **Disease trajectory showing "Unknown" stage** -- Some trajectories require specific marker combinations. For example, the hepatic trajectory needs ALT, AST, and albumin together. A single liver enzyme without context is insufficient.
- **Slow report generation** -- PDF rendering adds 1-2s overhead. For faster iteration, request Markdown format first, then generate PDF only for the final version.

---

### 3.4 Cardiology Intelligence Agent (:8126)

**Overview and Clinical Purpose**

The Cardiology Intelligence Agent synthesizes cardiac imaging, electrophysiology, hemodynamics, heart failure management, valvular disease, preventive cardiology, interventional data, and cardio-oncology surveillance into ACC/AHA/ESC guideline-aligned recommendations. It includes 6 validated risk calculators (ASCVD, HEART Score, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II), a GDMT optimizer, and 11 clinical workflows. All processing runs on a single DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `CARDIO_AGENT_PORT` | `8126` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model ID |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | 384-dim sentence-transformer |
| `TOP_K` | `15` | Max chunks per collection search |
| `GDMT_TITRATION_STEPS` | `4` | GDMT medication classes optimized simultaneously (ARNI/ACEi, beta-blocker, MRA, SGLT2i) |
| `RISK_CALCULATOR_COUNT` | `6` | Number of validated risk calculators |

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus, RAG engine, GDMT optimizer, and calculator status |
| GET | `/collections` | List loaded collections with record counts |
| GET | `/workflows` | List all 11 clinical workflow definitions |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/v1/cardio/integrated-assessment` | Full multi-collection integrated cardiac assessment |
| POST | `/v1/cardio/query` | Free-text RAG query across all cardiology collections |
| POST | `/v1/cardio/search` | Direct vector similarity search with collection filtering |
| POST | `/v1/cardio/find-related` | Find related entities (drugs, conditions, biomarkers) |
| POST | `/v1/cardio/risk/ascvd` | ASCVD 10-year risk (Pooled Cohort Equations) |
| POST | `/v1/cardio/risk/heart-score` | HEART Score for ED chest pain risk stratification |
| POST | `/v1/cardio/risk/cha2ds2-vasc` | CHA2DS2-VASc stroke risk for atrial fibrillation |
| POST | `/v1/cardio/risk/has-bled` | HAS-BLED anticoagulation bleeding risk |
| POST | `/v1/cardio/risk/maggic` | MAGGIC heart failure mortality risk |
| POST | `/v1/cardio/risk/euroscore` | EuroSCORE II cardiac surgical mortality |
| POST | `/v1/cardio/gdmt/optimize` | GDMT titration optimization for heart failure |
| POST | `/v1/cardio/workflow/cad` | Coronary artery disease assessment workflow |
| POST | `/v1/cardio/workflow/heart-failure` | Heart failure classification and management |
| POST | `/v1/cardio/workflow/valvular` | Valvular heart disease workflow |
| POST | `/v1/cardio/workflow/arrhythmia` | Arrhythmia and EP assessment |
| POST | `/v1/cardio/workflow/cardiac-mri` | Cardiac MRI tissue characterization |
| POST | `/v1/cardio/workflow/stress-test` | Stress testing protocol |
| POST | `/v1/cardio/workflow/prevention` | Cardiovascular prevention and risk stratification |
| POST | `/v1/cardio/workflow/cardio-oncology` | Cardio-oncology surveillance |
| POST | `/v1/cardio/cross-modal/evaluate` | Cross-modal genomic-cardiac evaluation |
| GET | `/v1/cardio/guidelines` | List available clinical guideline references |
| GET | `/v1/cardio/conditions` | List supported cardiovascular conditions |
| GET | `/v1/cardio/biomarkers` | List tracked cardiac biomarkers |
| GET | `/v1/cardio/drugs` | List cardiovascular drugs in knowledge base |
| GET | `/v1/cardio/genes` | List cardiomyopathy-associated genes |
| GET | `/v1/cardio/knowledge-version` | Knowledge base version and last update timestamp |
| POST | `/v1/reports/generate` | Generate clinical report (Markdown, JSON, PDF, FHIR R4) |
| GET | `/v1/reports/formats` | List available report export formats |
| GET | `/v1/events/stream` | SSE stream for cross-agent events |
| GET | `/v1/events/health` | Event system health status |

**Detailed Clinical Workflow Walkthrough (11 workflows)**

1. **Coronary Artery Disease (CAD)** -- Calcium score interpretation, CAD-RADS classification, plaque characterization, revascularization planning (PCI vs CABG), and ischemia-guided management.
2. **Heart Failure Management** -- HFrEF/HFpEF/HFmrEF classification by ejection fraction, GDMT optimization (4-pillar: ARNI/ACEi, beta-blocker, MRA, SGLT2i), device evaluation (ICD/CRT), and hemodynamic profiling.
3. **Valvular Heart Disease** -- Severity grading (mild/moderate/severe) by echocardiographic criteria, intervention timing per ACC/AHA guidelines, and prosthetic valve follow-up.
4. **Arrhythmia and EP** -- AF management with CHA2DS2-VASc and HAS-BLED scoring, ablation candidacy assessment, and device programming optimization.
5. **Cardiac MRI** -- LGE pattern interpretation, T1/T2 mapping analysis, strain assessment, and tissue characterization for cardiomyopathy diagnosis.
6. **Stress Testing** -- Exercise and pharmacologic stress interpretation, Duke treadmill score, perfusion defect classification, and ischemia quantification.
7. **Cardiovascular Prevention** -- ASCVD risk stratification, lipid management (statin/PCSK9i/ezetimibe selection), lifestyle Rx, and risk enhancer assessment.
8. **Cardio-Oncology** -- Chemotherapy cardiotoxicity surveillance, GLS tracking (global longitudinal strain), troponin/BNP monitoring, and CTRCD detection.
9. **Acute Decompensated HF** -- Hemodynamic profiling (warm-wet/cold-wet/cold-dry/warm-dry), IV diuretic dosing, inotrope selection, and MCS escalation criteria.
10. **Post-MI Secondary Prevention** -- Reperfusion assessment, DAPT strategy duration, beta-blocker/statin/ACEi titration, cardiac rehab referral, and ICD timing.
11. **Myocarditis/Pericarditis** -- Lake Louise CMR criteria evaluation, biopsy indications, NSAIDs/colchicine dosing, and activity restriction guidelines.

**Example Query and Response**

```bash
curl -X POST http://localhost:8126/v1/cardio/risk/heart-score \
  -H "Content-Type: application/json" \
  -d '{"age": 62, "history": "moderately suspicious", "ecg": "nonspecific ST changes", "troponin": "normal", "risk_factors": 3}'
```

```json
{
  "calculator": "HEART Score",
  "score": 6,
  "risk_category": "Moderate",
  "interpretation": "HEART Score 6/10 (Moderate risk). 30-day MACE rate: 12-17%. Recommend observation, serial troponins, and non-invasive testing within 72 hours.",
  "components": {
    "history": {"score": 1, "detail": "Moderately suspicious for ACS"},
    "ecg": {"score": 1, "detail": "Nonspecific ST-segment changes"},
    "age": {"score": 2, "detail": "Age 62 (>= 45 and < 65)"},
    "risk_factors": {"score": 2, "detail": "3+ risk factors (>= 3)"},
    "troponin": {"score": 0, "detail": "Normal troponin"}
  },
  "guideline_reference": "ACC/AHA 2021 Chest Pain Guideline, Section 5.3",
  "confidence": 0.95,
  "processing_time_ms": 280
}
```

**Cross-Agent Integration**

- **Sends** cardiotoxicity alerts to Oncology Agent (:8527) when chemotherapy candidates show cardiac risk (pre-treatment EF, GLS baseline, anthracycline cumulative dose tracking)
- **Receives** biomarker trends (BNP, troponin, CRP) from Biomarker Agent (:8529) for longitudinal heart failure monitoring
- **Queries** `genomic_evidence` for cardiomyopathy-associated variants (TTN, LMNA, MYH7, MYBPC3, SCN5A) to integrate genetic risk into clinical assessment
- **Receives** stroke triage events from Neurology Agent (:8528) to trigger cardiac workup (AF screening, echocardiography) after cryptogenic stroke
- **Provides** cardiotoxicity risk context to CAR-T Agent (:8522) for patients receiving lymphodepleting chemotherapy

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Risk calculator (any of 6) | 50-150 ms |
| GDMT optimization | 200-500 ms |
| Single collection search (top-15) | 80-150 ms |
| Full 13-collection parallel search | 250-450 ms |
| Claude LLM synthesis | 2,000-3,500 ms |
| Full integrated assessment | 3,500-5,500 ms |
| Report generation (PDF) | 1,500-2,500 ms |

**Common Pitfalls and Troubleshooting**

- **ASCVD risk returns "out of range"** -- The Pooled Cohort Equations are validated for ages 40-79. Inputs outside this range receive a warning and extrapolated estimate.
- **HEART Score confusion with Framingham** -- The HEART Score is for acute chest pain triage in the ED (History, ECG, Age, Risk factors, Troponin), not for long-term cardiovascular risk. Use ASCVD (PCE) for 10-year risk prediction.
- **GDMT optimizer returns "contraindicated"** -- The optimizer checks for absolute contraindications (hyperkalemia for MRA, bradycardia for beta-blockers, angioedema history for ACEi). Review the contraindication reason in the response before overriding.
- **Missing workflow endpoints** -- Workflows for acute decompensated HF, post-MI, and myocarditis/pericarditis use the generic `/v1/cardio/query` endpoint with structured parameters rather than dedicated workflow routes.
- **Cross-modal evaluation empty** -- The `/v1/cardio/cross-modal/evaluate` endpoint requires `genomic_evidence` collection to be loaded. Verify via `GET /collections` that it shows 3.56M records.

---

### 3.5 Neurology Intelligence Agent (:8528)

**Overview and Clinical Purpose**

The Neurology Intelligence Agent supports 8 specialized clinical workflows plus a general neurology RAG query mode (9 total), covering acute stroke triage, dementia evaluation, epilepsy classification, brain tumor grading, MS monitoring, Parkinson's assessment, headache classification, and neuromuscular evaluation. It integrates guidelines from AAN, AHA/ASA, ILAE, ICHD-3, WHO CNS 2021, McDonald 2017, and MDS criteria, with 10 validated clinical scale calculators. All processing runs on a single DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `NEURO_AGENT_PORT` | `8528` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model ID |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | 384-dim sentence-transformer |
| `TOP_K` | `15` | Max chunks per collection search |
| `STROKE_WINDOW_HOURS` | `4.5` | tPA eligibility window in hours from symptom onset |
| `THROMBECTOMY_WINDOW_HOURS` | `24` | Mechanical thrombectomy eligibility window (DAWN/DEFUSE-3) |

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
| ALSFRS-R | ALS Functional Rating Scale -- Revised |
| ASPECTS | Alberta Stroke Programme Early CT Score |
| Hoehn-Yahr | Parkinson's disease staging |

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus, RAG engine, and workflow engine status |
| GET | `/collections` | List loaded collections with record counts |
| GET | `/workflows` | List all 9 workflow definitions |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/v1/neuro/integrated-assessment` | Full multi-collection integrated neurological assessment |
| POST | `/v1/neuro/query` | Free-text RAG query across all neurology collections |
| POST | `/v1/neuro/search` | Direct vector similarity search with collection filtering |
| POST | `/v1/neuro/scale/calculate` | Clinical scale calculation (any of 10 scales) |
| POST | `/v1/neuro/stroke/triage` | Acute stroke triage with NIHSS, ASPECTS, thrombolysis/thrombectomy eligibility |
| POST | `/v1/neuro/dementia/evaluate` | Dementia evaluation with MoCA, ATN staging, differential diagnosis |
| POST | `/v1/neuro/epilepsy/classify` | Epilepsy classification (ILAE 2017), syndrome identification, EEG-MRI concordance |
| POST | `/v1/neuro/tumor/grade` | Brain tumor grading (WHO CNS 2021) with molecular markers (IDH, MGMT, 1p/19q) |
| POST | `/v1/neuro/ms/assess` | MS disease monitoring with EDSS, NEDA-3/4, DMT escalation evaluation |
| POST | `/v1/neuro/parkinsons/assess` | Parkinson's assessment with MDS-UPDRS III, Hoehn-Yahr, motor subtype classification |
| POST | `/v1/neuro/headache/classify` | Headache classification (ICHD-3), HIT-6 scoring, red flag screening |
| POST | `/v1/neuro/neuromuscular/evaluate` | Neuromuscular evaluation with ALSFRS-R, EMG/NCS pattern analysis |
| POST | `/v1/neuro/workflow/{workflow_type}` | Generic workflow endpoint for any of 9 workflow types |
| GET | `/v1/neuro/domains` | List supported neurological domains |
| GET | `/v1/neuro/scales` | List available clinical scales with input requirements |
| GET | `/v1/neuro/guidelines` | List integrated clinical guidelines |
| GET | `/v1/neuro/knowledge-version` | Knowledge base version and last update timestamp |
| POST | `/v1/reports/generate` | Generate clinical report (Markdown, JSON, PDF, FHIR R4) |
| GET | `/v1/reports/formats` | List available report export formats |
| GET | `/v1/events/stream` | SSE stream for cross-agent events |
| GET | `/v1/events/health` | Event system health status |

**Detailed Clinical Workflow Walkthrough**

1. **Acute Stroke Triage** -- Clinician enters onset time, NIHSS score, and CT findings. The agent calculates time-to-treatment, evaluates tPA eligibility (<4.5 hours per AHA/ASA guidelines), assesses thrombectomy candidacy (<24 hours, DAWN/DEFUSE-3 criteria), computes ASPECTS score for anterior circulation strokes, and localizes vascular territory.
2. **Dementia Evaluation** -- MoCA score, age, and symptom list are submitted. The agent performs ATN biomarker staging (Amyloid/Tau/Neurodegeneration), generates a differential diagnosis ranked by probability (AD, FTD, LBD, VaD, NPH, CJD), recommends further workup (PET, CSF, MRI volumetrics), and suggests treatment options.
3. **Epilepsy Classification** -- Seizure semiology and EEG/MRI findings are entered. The agent classifies seizures per ILAE 2017 (focal, generalized, unknown onset), identifies epilepsy syndromes, performs EEG-MRI concordance analysis for surgical candidacy, and recommends ASMs based on seizure type with PGx context.
4. **Brain Tumor Grading** -- Histology and molecular markers (IDH mutation, MGMT methylation, 1p/19q codeletion, ATRX, H3K27M) are submitted. The agent applies WHO CNS 2021 classification, determines integrated diagnosis (e.g., "Astrocytoma, IDH-mutant, grade 3"), and maps to NCCN treatment guidelines.
5. **MS Disease Monitoring** -- EDSS score, relapse history, and MRI lesion data are entered. The agent assesses NEDA-3/4 status (No Evidence of Disease Activity), evaluates DMT escalation criteria, calculates relapse risk stratification, and tracks lesion burden over time.
6. **Parkinson's Assessment** -- MDS-UPDRS Part III motor items are scored. The agent computes total score, determines Hoehn-Yahr stage, classifies motor subtype (tremor-dominant, PIGD, indeterminate), and optimizes medication regimen (levodopa, dopamine agonists, MAO-B inhibitors).
7. **Headache Classification** -- Attack characteristics (duration, location, quality, associated symptoms) are entered. The agent applies ICHD-3 criteria, computes HIT-6 disability score, screens for red flags (thunderclap, positional, progressive), and recommends preventive therapy if indicated.
8. **Neuromuscular Evaluation** -- Weakness pattern, reflexes, and EMG/NCS findings are submitted. The agent computes ALSFRS-R functional score, analyzes EMG/NCS patterns (neuropathic vs myopathic vs NMJ), guides genetic testing, and recommends specialty referral.
9. **General Neurology Query** -- Free-form RAG-powered Q&A across all neurology knowledge collections, not bound to a specific workflow.

**Example Query and Response**

```bash
curl -X POST http://localhost:8528/v1/neuro/stroke/triage \
  -H "Content-Type: application/json" \
  -d '{"onset_time": "2026-03-29T08:30:00", "current_time": "2026-03-29T10:45:00", "nihss_score": 14, "ct_findings": "no hemorrhage", "aspects_score": 8}'
```

```json
{
  "triage_result": {
    "time_from_onset_minutes": 135,
    "tpa_eligible": true,
    "tpa_window_remaining_minutes": 135,
    "thrombectomy_eligible": true,
    "thrombectomy_window_remaining_minutes": 1305,
    "nihss_severity": "Moderate-Severe",
    "aspects_interpretation": "Score 8/10 -- favorable for intervention (>= 6 threshold)",
    "vascular_territory": "MCA (likely M1 or M2 occlusion based on NIHSS pattern)",
    "recommendation": "URGENT: Administer IV tPA immediately (within 4.5-hour window). Activate neurointerventional team for mechanical thrombectomy evaluation. Obtain CTA head and neck to confirm large vessel occlusion."
  },
  "guideline_references": [
    "AHA/ASA 2019 Acute Ischemic Stroke Guidelines, Section 3.2 (tPA)",
    "DAWN Trial: thrombectomy 6-24 hours with clinical-core mismatch",
    "DEFUSE-3 Trial: thrombectomy 6-16 hours with perfusion mismatch"
  ],
  "evidence": [
    {"collection": "neuro_stroke", "score": 0.93, "text": "IV alteplase within 4.5 hours of symptom onset..."},
    {"collection": "neuro_guidelines", "score": 0.91, "text": "AHA/ASA Class I recommendation for mechanical thrombectomy..."}
  ],
  "confidence": 0.96,
  "processing_time_ms": 1850
}
```

**Cross-Agent Integration**

- **Queries** Imaging Agent (:8524) for MRI findings (MS lesion counts, brain volumetrics, DWI/ADC maps for stroke)
- **Receives** pharmacogenomic context from PGx Agent (:8107) for antiepileptic drug selection (HLA-B*15:02 and carbamazepine, CYP2C9 and phenytoin)
- **Shares** stroke triage events with Cardiology Agent (:8126) for cardiac workup -- AF screening, echocardiography, and Holter monitoring after cryptogenic stroke
- **Queries** `genomic_evidence` for neurogenetic variants (APP, PSEN1/2, GBA, LRRK2, PARK2, SMN1) to integrate genetic risk into clinical assessment
- **Receives** tumor molecular markers from Oncology Agent (:8527) for integrated brain tumor classification

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Clinical scale calculation (any of 10) | 20-80 ms |
| Single collection search (top-15) | 80-150 ms |
| Full multi-collection parallel search | 250-450 ms |
| Stroke triage (end-to-end, no LLM) | 200-500 ms |
| Claude LLM synthesis | 2,000-3,500 ms |
| Full workflow assessment (end-to-end) | 3,000-5,500 ms |
| Report generation (PDF) | 1,500-2,500 ms |

**Common Pitfalls and Troubleshooting**

- **Stroke triage returns "window expired"** -- Verify that `onset_time` and `current_time` are in ISO 8601 format with timezone. If `current_time` is omitted, the server uses UTC now, which may differ from local time.
- **NIHSS score mismatch** -- The agent accepts either a pre-computed total score or individual item scores. If both are provided, individual items take precedence and the total is recalculated. Ensure all 15 items are included for accurate scoring.
- **MoCA score interpretation varies by education** -- The agent applies the standard +1 point adjustment for education <= 12 years. Provide the raw score without manual adjustment; the agent handles the correction.
- **MS NEDA-3 assessment incomplete** -- NEDA-3 requires three inputs: relapse count (0), new/enlarging T2 lesions (0), and disability progression (stable EDSS). Missing any component results in a partial assessment.
- **Workflow type not found** -- Valid workflow types are: `acute_stroke_triage`, `dementia_evaluation`, `epilepsy_focus`, `brain_tumor_grading`, `ms_monitoring`, `parkinsons_assessment`, `headache_classification`, `neuromuscular_evaluation`, `general`. Use the exact ID from `GET /workflows`.

---

### 3.6 Precision Autoimmune Agent (:8532)

**Overview and Clinical Purpose**

The Precision Autoimmune Agent addresses the diagnostic odyssey that autoimmune patients face -- often 3-6 years across multiple specialists before diagnosis. It interprets autoantibody panels, HLA typing, biomarker trends, and genomic data across 13 autoimmune conditions. It features disease activity scoring (DAS28-CRP, SLEDAI-2K, CDAI, BASDAI), flare prediction, biologic therapy recommendations with PGx context, and diagnostic odyssey analysis. Also detects the POTS/hEDS/MCAS triad and overlap syndromes. All processing runs on a single DGX Spark ($4,699, March 2026).

**Advanced Configuration Options**

| Variable | Default | Description |
|----------|---------|-------------|
| `AUTOIMMUNE_AGENT_PORT` | `8532` | FastAPI listen port |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MILVUS_POOL_SIZE` | `10` | Connection pool size |
| `ANTHROPIC_API_KEY` | (required) | Claude API key for LLM synthesis |
| `LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model ID; upgrade to `claude-opus-4-20250514` for complex diagnostic odyssey analyses |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | 384-dim sentence-transformer |
| `TOP_K` | `15` | Max chunks per collection search |
| `FLARE_RISK_THRESHOLD` | `0.65` | Risk score threshold for flare alerts (0-1 scale) |
| `HLA_ODDS_RATIO_MIN` | `2.0` | Minimum odds ratio for HLA-disease associations to flag |
| `COLLECTION_WEIGHTS` | (see table) | Weighted relevance for each collection in composite scoring |

**13 Supported Conditions**

Rheumatoid arthritis, SLE, multiple sclerosis, type 1 diabetes, inflammatory bowel disease, psoriasis/psoriatic arthritis, ankylosing spondylitis, Sjogren's syndrome, systemic sclerosis, myasthenia gravis, celiac disease, Graves' disease, Hashimoto's thyroiditis.

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus and engine status |
| GET | `/healthz` | Kubernetes-style liveness probe |
| GET | `/collections` | List loaded collections with record counts |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/v1/autoimmune/integrated-assessment` | Full multi-collection integrated autoimmune assessment |
| POST | `/query` | Free-text RAG query across all autoimmune collections |
| POST | `/query/stream` | Streaming RAG query (Server-Sent Events) |
| POST | `/search` | Direct vector similarity search with collection filtering |
| POST | `/analyze` | Full patient analysis (antibodies + HLA + biomarkers combined) |
| POST | `/differential` | Differential diagnosis from symptom/antibody profile |
| POST | `/ingest/upload` | Upload and ingest clinical documents (PDF, DOCX, text) |
| POST | `/ingest/demo-data` | Load demo patient dataset for testing |
| POST | `/collections/create` | Create custom collections for institution-specific data |
| POST | `/export` | Export analysis results (Markdown, JSON, PDF, FHIR R4) |

**Detailed Clinical Workflow Walkthrough**

1. **Document ingestion** -- Clinical documents (progress notes, lab reports, imaging reports) are uploaded via `/ingest/upload` or Streamlit UI. Documents are chunked, embedded with BGE-small-en-v1.5, and stored in `autoimmune_clinical_documents`.
2. **Autoantibody panel interpretation** -- Antibody results (ANA, anti-dsDNA, anti-CCP, RF, anti-SSA/SSB, anti-Scl-70, anti-Jo-1, ANCA, anti-centromere, anti-Smith, anti-RNP, anti-TPO, TSI, anti-AChR, anti-tTG) are mapped to disease associations with published sensitivity and specificity. For example: anti-CCP positive with RF positive has 95% specificity for RA.
3. **HLA typing analysis** -- HLA alleles are evaluated against 50+ allele-disease associations with odds ratios. Key associations: HLA-B*27:05 and ankylosing spondylitis (OR=87.4), HLA-DRB1*04:01 and RA (OR=4.2), HLA-DQ2/DQ8 and celiac disease (OR=7.0).
4. **Disease activity scoring** -- Standardized scores are calculated based on disease:
   - RA: DAS28-CRP (remission <2.6, low 2.6-3.2, moderate 3.2-5.1, high >5.1)
   - SLE: SLEDAI-2K (no activity 0, mild 1-5, moderate 6-10, high 11-19, very high >=20)
   - IBD: CDAI for Crohn's (remission <150, mild 150-219, moderate 220-450, severe >450)
   - AS: BASDAI (low <4, high >=4)
5. **Flare prediction** -- Analyzes longitudinal biomarker patterns (CRP, ESR, IL-6, complement C3/C4, calprotectin, ferritin) against known flare precursors. Returns a 0-1 risk score with contributing factors. Scores above 0.65 trigger an alert to the Biomarker Agent.
6. **Biologic therapy recommendations** -- Filtered by confirmed diagnosis, disease activity level, prior therapy failures, and PGx context (FCGR3A V/F polymorphism for rituximab response, TPMT for azathioprine dosing). Recommendations include TNF inhibitors, anti-CD20, IL-6R blockers, JAK inhibitors, and IL-17/23 inhibitors.
7. **Diagnostic odyssey analysis** -- Surfaces patterns missed across fragmented multi-specialist records. Identifies symptom constellations that span years and multiple providers, flags potential overlap syndromes (e.g., RA + Sjogren's, SLE + APS), and detects the POTS/hEDS/MCAS triad.

**Example Query and Response**

```bash
curl -X POST http://localhost:8532/analyze \
  -H "Content-Type: application/json" \
  -d '{"antibodies": {"ANA": "positive", "anti_dsDNA": 85, "anti_Smith": "positive"}, "hla": ["DRB1*15:01"], "biomarkers": {"CRP": 18.5, "ESR": 42, "C3": 62, "C4": 7, "anti_dsDNA_titer": 120}}'
```

```json
{
  "primary_diagnosis": {
    "disease": "Systemic Lupus Erythematosus (SLE)",
    "confidence": 0.92,
    "acr_eular_score": 18,
    "classification_criteria": "ACR/EULAR 2019 SLE Classification (threshold >= 10)"
  },
  "disease_activity": {
    "score_type": "SLEDAI-2K",
    "score": 12,
    "level": "High Activity",
    "active_domains": ["immunologic (anti-dsDNA, low complement)", "serositis risk"]
  },
  "flare_risk": {
    "score": 0.78,
    "risk_level": "High",
    "contributing_factors": ["Rising anti-dsDNA (120 IU/mL)", "Low C4 (7 mg/dL)", "Low C3 (62 mg/dL)"],
    "recommendation": "Consider preemptive dose adjustment; monitor weekly for 4 weeks"
  },
  "hla_associations": [
    {"allele": "DRB1*15:01", "disease": "SLE", "odds_ratio": 2.1, "significance": "Moderate risk association"}
  ],
  "therapy_recommendations": [
    {"rank": 1, "drug": "Belimumab (anti-BLyS)", "rationale": "Add-on for active SLE with high anti-dsDNA and low complement"},
    {"rank": 2, "drug": "Mycophenolate mofetil", "rationale": "Steroid-sparing immunosuppression for moderate-severe SLE"},
    {"rank": 3, "drug": "Anifrolumab (anti-IFNAR1)", "rationale": "Type I IFN pathway blockade for refractory skin/joint SLE"}
  ],
  "evidence": [
    {"collection": "autoimmune_disease_criteria", "score": 0.94, "text": "ACR/EULAR 2019: Anti-dsDNA + low complement = 10 points..."},
    {"collection": "autoimmune_biologic_therapies", "score": 0.89, "text": "BLISS-76: Belimumab + SOC improved SRI-4 response vs placebo..."}
  ],
  "processing_time_ms": 4120
}
```

**Cross-Agent Integration**

- **Receives** longitudinal biomarker trends from Biomarker Agent (:8529) for disease monitoring -- CRP, ESR, complement, calprotectin tracked over time with trend analysis
- **Requests** imaging assessment from Imaging Agent (:8524) for joint erosion scoring (MRI/ultrasound), organ evaluation (renal, pulmonary), and skin lesion characterization
- **Publishes** diagnosis events to other agents via SSE event bus -- new autoimmune diagnoses trigger downstream assessments in Cardiology (pericarditis risk), Neurology (CNS lupus, MS overlap), and PGx (immunosuppressant metabolism)
- **Queries** `genomic_evidence` for HLA alleles and autoimmune-associated variants (CTLA4, PTPN22, IL2RA, STAT4) to integrate genetic risk
- **Sends** flare alerts to Biomarker Agent (:8529) when biomarker patterns indicate impending flare, triggering intensified monitoring protocols

**Performance Characteristics**

| Operation | Typical Latency |
|-----------|----------------|
| Single collection search (top-15) | 80-150 ms |
| Full 14-collection weighted search | 300-500 ms |
| Disease activity calculation | 30-80 ms |
| Flare risk prediction | 50-120 ms |
| Autoantibody panel interpretation | 40-100 ms |
| HLA association lookup | 20-50 ms |
| Claude LLM synthesis | 2,500-4,000 ms |
| Full analysis (end-to-end) | 3,500-5,500 ms |
| Document ingestion (per page) | 500-1,000 ms |

**Common Pitfalls and Troubleshooting**

- **ANA "positive" without titer** -- The agent accepts qualitative ("positive"/"negative") or quantitative (titer, e.g., "1:640") ANA values. Quantitative values enable more precise disease association scoring. Include titer and pattern (homogeneous, speckled, nucleolar, centromere) when available.
- **Disease activity score returns "insufficient data"** -- Each score requires specific inputs. DAS28-CRP needs tender joint count (28 joints), swollen joint count (28 joints), CRP (mg/L), and patient global assessment (0-100 VAS). Missing any component prevents calculation.
- **HLA alleles not recognized** -- Use standard nomenclature (e.g., "B*27:05" not "B27" or "HLA-B27"). The agent accepts both 2-field (B*27:05) and 1-field (B*27) resolution but 2-field provides more accurate odds ratios.
- **Flare prediction unreliable with single timepoint** -- The prediction model works best with longitudinal data (3+ timepoints). A single snapshot can estimate current risk but cannot detect trajectory-based patterns (rising anti-dsDNA, falling complement).
- **Overlap syndromes missed** -- The diagnostic odyssey analysis requires documents from multiple visits and specialties. Single-visit data limits the agent's ability to detect cross-specialist patterns. Upload historical records for comprehensive analysis.
- **Demo data ingestion fails** -- Ensure the `/data/demo` directory exists and contains the sample patient files. Run `curl http://localhost:8532/health` first to verify Milvus connectivity before ingesting.

---

### 3.7 Rare Disease Diagnostic Agent (:8134)

**Overview and Clinical Purpose**

The Rare Disease Diagnostic Agent provides differential diagnosis, ACMG/AMP variant interpretation (implementing 23 ACMG criteria), HPO-based phenotype matching with information content (IC) weighted similarity scoring, therapeutic option search (orphan drugs, gene therapy, enzyme replacement), and clinical trial eligibility assessment. It covers 13 disease categories with 88 diseases in its knowledge base and tracks 25 approved or investigational gene therapies.

**Advanced Configuration Options**

All environment variables use the `RD_` prefix. Key settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `RD_MILVUS_HOST` | `localhost` | Milvus vector database host |
| `RD_MILVUS_PORT` | `19530` | Milvus port |
| `RD_API_PORT` | `8134` | FastAPI REST server port |
| `RD_STREAMLIT_PORT` | `8544` | 5-tab Streamlit UI port |
| `RD_LLM_MODEL` | `claude-sonnet-4-6` | Claude model for synthesis |
| `RD_ANTHROPIC_API_KEY` | (none) | Required for LLM; search-only mode without it |
| `RD_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity for evidence retrieval |
| `RD_TOP_K_PHENOTYPES` | `50` | Max results from `rd_phenotypes` per query |
| `RD_TOP_K_VARIANTS` | `100` | Max results from `rd_variants` per query |
| `RD_INGEST_SCHEDULE_HOURS` | `24` | Auto-ingest cycle (set `RD_INGEST_ENABLED=true`) |
| `RD_MAX_CONVERSATION_CONTEXT` | `3` | Prior exchanges injected for multi-turn |
| `RD_CITATION_HIGH_THRESHOLD` | `0.75` | Score above which citations are marked high-confidence |
| `RD_API_KEY` | (empty) | Set to require X-API-Key header on all endpoints |
| `RD_CROSS_AGENT_TIMEOUT` | `30` | Seconds before cross-agent HTTP calls time out |

Collection weight tuning: 14 `RD_WEIGHT_*` variables (e.g., `RD_WEIGHT_PHENOTYPES=0.12`) control the relative influence of each collection during multi-collection fusion. Weights must sum to approximately 1.0 (tolerance 0.05).

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with component readiness and vector counts |
| GET | `/collections` | Collection names and record counts |
| GET | `/workflows` | 11 available diagnostic workflow definitions |
| GET | `/metrics` | Prometheus-compatible metrics export |
| POST | `/v1/diagnostic/query` | RAG Q&A across all rare disease collections |
| POST | `/v1/diagnostic/search` | Multi-collection evidence search (no LLM) |
| POST | `/v1/diagnostic/diagnose` | Differential diagnosis from HPO terms |
| POST | `/v1/diagnostic/variants/interpret` | ACMG/AMP variant classification |
| POST | `/v1/diagnostic/phenotype/match` | HPO-to-disease phenotype matching |
| POST | `/v1/diagnostic/therapy/search` | Therapeutic option search (ERT, SRT, gene therapy) |
| POST | `/v1/diagnostic/trial/match` | Clinical trial eligibility assessment |
| POST | `/v1/diagnostic/workflow/{type}` | Generic workflow dispatch (11 types) |
| POST | `/v1/diagnostic/integrated-assessment` | Cross-agent multi-agent assessment |
| GET | `/v1/diagnostic/disease-categories` | Reference catalog of 13 disease categories |
| GET | `/v1/diagnostic/gene-therapies` | 25 approved/investigational gene therapies |
| GET | `/v1/diagnostic/acmg-criteria` | 23 ACMG criteria reference |
| GET | `/v1/diagnostic/hpo-categories` | HPO top-level ontology terms |
| GET | `/v1/diagnostic/knowledge-version` | Knowledge base version metadata |
| POST | `/v1/reports/generate` | Multi-format report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | Server-sent event stream |

**Detailed Clinical Workflow Walkthrough**

1. **Phenotype entry** -- Clinician enters patient phenotypes as HPO terms (e.g., HP:0001250 Seizure, HP:0001263 Global developmental delay) via the 5-tab Streamlit UI or the `/v1/diagnostic/diagnose` API. Age of onset category (antenatal, neonatal, infantile, childhood, juvenile, adult) is provided as context.

2. **HPO-to-disease matching** -- The agent embeds the HPO term descriptions using BGE-small-en-v1.5, searches the `rd_phenotypes` collection (top_k=50), and computes IC-weighted Resnik similarity against OMIM and Orphanet disease phenotype profiles. Each disease receives a composite phenotype match score.

3. **Differential diagnosis ranking** -- Candidate diseases are ranked by phenotype match score. The LLM synthesizes evidence from `rd_diseases`, `rd_literature`, and `rd_case_reports` to assign confidence levels (high/moderate/low) and highlight discriminating features.

4. **Variant interpretation** -- If VCF data is available, candidate variants in genes associated with top-ranked diseases are classified using 23 ACMG/AMP criteria: PVS1, PM1-PM6, PP1-PP5, BA1, BS1-BS4, BP1-BP7. ClinVar, gnomAD allele frequency, AlphaMissense scores, and in-silico predictors are integrated.

5. **Therapeutic search** -- For confirmed or suspected diagnoses, the agent queries `rd_therapies` for orphan drugs (FDA/EMA approved), gene therapies (25 tracked), enzyme replacement therapies, substrate reduction therapies, and investigational compounds.

6. **Trial eligibility** -- The agent cross-references patient age, gene, diagnosis, and variant status against `rd_trials` and triggers the Clinical Trial Agent for expanded matching.

7. **Report export** -- Final report exported as Markdown, JSON, PDF, FHIR R4 DiagnosticReport, or GA4GH Phenopacket v2.

**Example Query and Response**

```bash
curl -X POST http://localhost:8134/v1/diagnostic/diagnose \
  -H "Content-Type: application/json" \
  -d '{"hpo_terms": ["HP:0001250", "HP:0001263", "HP:0002079"], "age_of_onset": "infantile"}'
```

```json
{
  "status": "completed",
  "differential_diagnosis": [
    {
      "rank": 1,
      "disease": "Dravet syndrome",
      "omim": "607208",
      "orphanet": "ORPHA:33069",
      "phenotype_match_score": 0.87,
      "confidence": "high",
      "key_gene": "SCN1A",
      "discriminating_features": ["Febrile seizures progressing to afebrile", "Temperature sensitivity"]
    },
    {
      "rank": 2,
      "disease": "CLN2 disease (late-infantile neuronal ceroid lipofuscinosis)",
      "omim": "204500",
      "orphanet": "ORPHA:228346",
      "phenotype_match_score": 0.74,
      "confidence": "moderate",
      "key_gene": "TPP1",
      "discriminating_features": ["Progressive vision loss", "Cerebellar atrophy on MRI"]
    }
  ],
  "hpo_terms_matched": 3,
  "diseases_evaluated": 88,
  "evidence_sources": 42,
  "processing_time_ms": 2340.5
}
```

**Cross-Agent Integration**

| Direction | Partner Agent | Trigger/Scenario |
|-----------|--------------|------------------|
| Outbound | Pharmacogenomics (:8107) | Queries drug metabolism context when rare disease therapeutics have PGx implications (e.g., ERT dosing) |
| Outbound | Clinical Trial (:8538) | Triggers expanded trial matching for rare disease patients via `/v1/trial/workflow/patient_matching/run` |
| Outbound | Imaging (:8524) | Requests imaging assessment for rare diseases with organ involvement (e.g., lysosomal storage diseases) |
| Outbound | Cardiology (:8126) | Queries cardiac phenotype context for cardiomyopathy-associated rare diseases (e.g., Fabry, Pompe) |
| Inbound | Biomarker (:8529) | Receives phenotype-genotype correlation requests |
| Bidirectional | All agents | Publishes diagnosis events via SSE event bus; reads `genomic_evidence` (shared) |

**Performance Characteristics**

| Operation | Typical Latency | Notes |
|-----------|----------------|-------|
| Phenotype match (3 HPO terms) | 800-1,200 ms | Vector search across `rd_phenotypes` (top_k=50) |
| Differential diagnosis (full) | 2,000-3,500 ms | Multi-collection search + LLM synthesis |
| ACMG variant interpretation | 1,500-2,500 ms | ClinVar + gnomAD + AlphaMissense lookup + LLM classification |
| Therapeutic search | 600-1,000 ms | `rd_therapies` search only |
| Report generation (PDF) | 3,000-5,000 ms | LLM narrative + template rendering |

**Common Pitfalls and Troubleshooting**

- **HPO term formatting**: Terms must use the `HP:NNNNNNN` format (7 digits, zero-padded). Using free-text symptom descriptions instead of HPO codes bypasses the IC-weighted matching and falls back to text similarity, which is less precise.
- **ACMG criteria incomplete**: The agent implements 23 of 28 ACMG criteria. Functional assay evidence (PS3) and segregation data (PP1 with quantitative scoring) require manual review.
- **Low phenotype match scores**: If all differential diagnoses score below 0.5, consider whether the phenotype is incomplete. Adding more HPO terms (especially specific sub-terms rather than general parent terms) significantly improves matching.
- **Gene therapy list currency**: The 25 gene therapies are current as of March 2026. New FDA/EMA approvals require a data refresh via `scripts/seed_knowledge.py`.
- **Cross-agent timeout**: If the PGx or Trial agent is down, the integrated assessment degrades gracefully (returns partial results). Check `RD_CROSS_AGENT_TIMEOUT` if timeouts are frequent.

---

### 3.8 Pharmacogenomics Intelligence Agent (:8107)

**Overview and Clinical Purpose**

The Pharmacogenomics (PGx) Intelligence Agent translates patient genotype data into actionable prescribing guidance. It interprets star alleles for 14 pharmacogenes, applies 9 CPIC guideline algorithms, screens for 12 HLA-mediated adverse drug reaction associations, models phenoconversion from drug-drug interactions, and provides genotype-guided dosing algorithms including the IWPC warfarin dosing calculator.

**Advanced Configuration Options**

All environment variables use the `PGX_` prefix. Key settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `PGX_MILVUS_HOST` | `localhost` | Milvus vector database host |
| `PGX_MILVUS_PORT` | `19530` | Milvus port |
| `PGX_API_PORT` | `8107` | FastAPI REST server port |
| `PGX_STREAMLIT_PORT` | `8507` | Streamlit UI port |
| `PGX_LLM_MODEL` | `claude-sonnet-4-6` | Claude model for synthesis |
| `PGX_ANTHROPIC_API_KEY` | (none) | Required for LLM features |
| `PGX_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity for evidence retrieval |
| `PGX_TOP_K_PER_COLLECTION` | `5` | Default results per collection |
| `PGX_INGEST_SCHEDULE_HOURS` | `168` | Weekly auto-ingest cycle (7 days) |
| `PGX_API_KEY` | (empty) | Set to require X-API-Key header |
| `PGX_CROSS_AGENT_TIMEOUT` | `30` | Cross-agent HTTP timeout in seconds |
| `PGX_RAG_PIPELINE_ROOT` | `/app/rag-chat-pipeline` | Path to shared RAG pipeline |

Collection weight tuning: 15 `PGX_WEIGHT_*` variables control multi-collection fusion. `PGX_WEIGHT_DRUG_GUIDELINES` defaults to 0.14 (highest weight) because CPIC/DPWG guideline matching is the most clinically relevant signal.

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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with Milvus and LLM status |
| GET | `/collections` | Collection names and record counts |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/query` | RAG Q&A across all PGx collections |
| POST | `/search` | Multi-collection evidence search (no LLM) |
| POST | `/profile` | Full PGx profile from diplotype map |
| POST | `/dosing` | Drug-specific dosing guidance |
| POST | `/hla-screen` | HLA hypersensitivity screening |
| POST | `/v1/pgx/drug-check` | Single-drug PGx check with CPIC alerts |
| POST | `/v1/pgx/medication-review` | Polypharmacy review with DDI detection |
| POST | `/v1/pgx/dosing/warfarin` | IWPC warfarin dosing algorithm |
| POST | `/v1/pgx/hla-screen` | HLA adverse reaction screening |
| POST | `/v1/pgx/phenoconversion` | Phenoconversion analysis from concomitant drugs |
| POST | `/v1/pgx/integrated-assessment` | Cross-agent multi-agent assessment |
| GET | `/v1/pgx/genes` | List all 14 pharmacogenes with metadata |
| GET | `/v1/pgx/drugs` | List all tracked drugs by therapeutic category |
| POST | `/v1/reports/generate` | Multi-format report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | Server-sent event stream |

**Detailed Clinical Workflow Walkthrough**

1. **Diplotype submission** -- Clinician submits patient diplotypes (e.g., `CYP2D6: *1/*4`, `CYP2C19: *1/*2`) from VCF-based star allele calling (Stargazer, Cyrius) or laboratory report. The `/profile` or `/v1/pgx/drug-check` endpoint receives the data.

2. **Metabolizer phenotype assignment** -- Each diplotype is mapped to an activity score and metabolizer phenotype: ultra-rapid metabolizer (UM), normal metabolizer (NM), intermediate metabolizer (IM), or poor metabolizer (PM). The mapping follows CPIC standardized terms.

3. **CPIC/DPWG guideline matching** -- For each gene-drug pair in the patient's medication list, the agent queries `pgx_drug_guidelines` for applicable CPIC Level A/B recommendations. 9 CPIC algorithms are implemented covering CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, SLCO1B1, and UGT1A1.

4. **Phenoconversion check** -- Current medications are screened against CYP inhibitor/inducer databases. A strong CYP2D6 inhibitor (e.g., fluoxetine, paroxetine) in a genotypic NM converts the effective phenotype to PM. The `/v1/pgx/phenoconversion` endpoint returns the precipitant drug, inhibitor strength, converted phenotype, and affected substrate drugs.

5. **HLA screening** -- 12 HLA-drug associations are checked. Key associations include HLA-B*57:01/abacavir (mandatory screening), HLA-B*58:01/allopurinol, HLA-B*15:02/carbamazepine, and HLA-A*31:01/carbamazepine. Population-specific prevalence data is returned.

6. **Dosing algorithm** -- For warfarin, the IWPC algorithm (Klein et al., NEJM 2009) calculates predicted weekly dose incorporating VKORC1, CYP2C9, age, height, weight, race, amiodarone use, and enzyme inducer status. Safety bounds clamp output to 3-105 mg/week.

7. **Alternative recommendations** -- For PM/UM patients, the agent retrieves genotype-appropriate alternative drugs from `pgx_drug_alternatives` and generates a structured report.

**Example Query and Response**

```bash
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "codeine", "gene": "CYP2D6", "phenotype": "poor_metabolizer", "diplotype": "*4/*4"}'
```

```json
{
  "drug": "codeine",
  "gene": "CYP2D6",
  "phenotype": "poor_metabolizer",
  "alerts": [
    {
      "alert_level": "critical",
      "gene": "CYP2D6",
      "drug": "codeine",
      "phenotype": "poor_metabolizer",
      "recommendation": "CPIC guideline exists for CYP2D6-codeine. Patient phenotype: poor_metabolizer. Consider alternative therapy or significant dose reduction.",
      "guideline_body": "CPIC",
      "cpic_level": "A",
      "alternative_drugs": ["morphine", "non-opioid analgesics", "acetaminophen"]
    }
  ],
  "knowledge_context": "CYP2D6 encodes cytochrome P450 2D6, which metabolizes codeine to morphine via O-demethylation. Poor metabolizers (*4/*4) produce insufficient morphine for analgesic effect.",
  "processing_time_ms": 45.2
}
```

**Cross-Agent Integration**

| Direction | Partner Agent | Trigger/Scenario |
|-----------|--------------|------------------|
| Outbound | Oncology (:8527) | Queries planned chemotherapy context for metabolism assessment (e.g., tamoxifen and CYP2D6) |
| Outbound | Cardiology (:8126) | Queries cardiac drug interaction context (e.g., clopidogrel and CYP2C19) |
| Outbound | Neurology (:8528) | Queries neurotoxic drug interaction context (e.g., carbamazepine and HLA-B*15:02) |
| Outbound | Clinical Trial (:8538) | Queries PGx-guided clinical trials |
| Inbound | All prescribing agents | Any agent recommending drug therapy can query PGx for metabolism context |
| Inbound | Biomarker (:8529) | Receives requests for genotype-adjusted biomarker reference ranges |

**Performance Characteristics**

| Operation | Typical Latency | Notes |
|-----------|----------------|-------|
| Single drug check | 30-80 ms | Knowledge graph lookup only, no vector search |
| Polypharmacy review (5 drugs) | 100-250 ms | Iterates CYP inhibitor/inducer databases |
| Warfarin IWPC dosing | 10-30 ms | Pure mathematical calculation |
| HLA screening | 20-60 ms | Dictionary lookup against 12 associations |
| Phenoconversion analysis | 40-120 ms | CYP inhibitor/inducer scan per drug |
| Full RAG query | 1,500-3,000 ms | 15-collection vector search + LLM synthesis |

**Common Pitfalls and Troubleshooting**

- **Star allele formatting**: Diplotypes must use the `*N/*N` format (e.g., `*1/*4`). Suballeles like `*1a` and `*5` are supported for SLCO1B1. Submitting raw VCF positions instead of called star alleles will not produce useful results.
- **Phenoconversion missed**: The phenoconversion engine only detects interactions from the concomitant drug list explicitly provided. It does not infer medications from other clinical data. Ensure the full medication list is submitted.
- **Population-specific allele frequencies**: HLA prevalence data varies significantly by ancestry. The agent reports population-specific frequencies (e.g., HLA-B*57:01 prevalence: 6-8% European, <1% East Asian) but does not automatically select the relevant population.
- **IWPC warfarin clamping**: Doses below 3 mg/week or above 105 mg/week are clamped to safety bounds. If the algorithm produces extreme values, verify input genotypes and clinical parameters.
- **Missing CPIC guideline**: Not all gene-drug pairs have CPIC Level A guidelines. The agent returns an `info` alert listing which drugs do have guidelines for the queried gene.

---

### 3.9 Imaging Intelligence Agent (:8524)

**Overview and Clinical Purpose**

The Imaging Intelligence Agent provides clinical decision support for radiology through a RAG knowledge system, 4 NVIDIA NIM microservices (VISTA-3D with 127 anatomical structures, MAISI, VILA-M3, Llama-3), and 4 reference clinical workflows. It supports Orthanc DICOM auto-ingestion, cross-modal genomics enrichment, NVIDIA FLARE federated learning, and FHIR R4 interoperability.

**Advanced Configuration Options**

All environment variables use the `IMAGING_` prefix. Key settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_MILVUS_HOST` | `localhost` | Milvus vector database host |
| `IMAGING_MILVUS_PORT` | `19530` | Milvus port |
| `IMAGING_API_PORT` | `8524` | FastAPI REST server port |
| `IMAGING_STREAMLIT_PORT` | `8525` | 9-tab Streamlit UI port |
| `IMAGING_LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model for synthesis |
| `IMAGING_NIM_MODE` | `local` | NIM mode: `local`, `cloud`, or `mock` |
| `IMAGING_NIM_LLM_URL` | `http://localhost:8520/v1` | Local Llama-3 NIM endpoint |
| `IMAGING_NIM_VISTA3D_URL` | `http://localhost:8530` | VISTA-3D segmentation endpoint |
| `IMAGING_NIM_MAISI_URL` | `http://localhost:8531` | MAISI synthetic CT endpoint |
| `IMAGING_NIM_VILAM3_URL` | `http://localhost:8532` | VILA-M3 vision-language endpoint |
| `IMAGING_NIM_ALLOW_MOCK_FALLBACK` | `true` | Fall back to mock responses if NIM unavailable |
| `IMAGING_NVIDIA_API_KEY` | (none) | Required for cloud NIM mode |
| `IMAGING_ORTHANC_URL` | `http://localhost:8042` | Orthanc DICOM server |
| `IMAGING_DICOM_AUTO_INGEST` | `false` | Enable automatic DICOM webhook processing |
| `IMAGING_DICOM_WATCH_INTERVAL` | `5` | Seconds between Orthanc `/changes` polls |
| `IMAGING_CROSS_MODAL_ENABLED` | `false` | Enable imaging-to-genomics cross-modal triggers |
| `IMAGING_PREVIEW_DEFAULT_FPS` | `8` | Frames per second for NIfTI preview videos |

**NVIDIA NIM Integration**

| NIM Service | Port | Capability |
|------------|------|------------|
| VISTA-3D | 8530 | 3D medical image segmentation (127 anatomical structures) |
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

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with collection stats and NIM status |
| GET | `/` | Service info with links to docs and health |
| GET | `/collections` | Collection names, counts, and labels |
| GET | `/metrics` | Prometheus metrics (Counter, Histogram) |
| GET | `/knowledge/stats` | Imaging domain knowledge graph statistics |
| POST | `/query` | Full RAG query: multi-collection retrieval + LLM synthesis |
| POST | `/search` | Evidence-only search (no LLM synthesis) |
| POST | `/find-related` | Cross-collection entity linking |
| POST | `/api/ask` | Meta-agent question answering |
| POST | `/nim/vista3d/segment` | VISTA-3D 3D segmentation (127 structures) |
| POST | `/nim/maisi/generate` | MAISI synthetic CT generation |
| POST | `/nim/vilam3/analyze` | VILA-M3 image understanding |
| POST | `/nim/llm/generate` | Llama-3 text generation |
| GET | `/nim/status` | NIM service availability check |
| POST | `/workflow/{name}/run` | Execute a clinical workflow |
| GET | `/workflows` | List available workflows |
| POST | `/events/dicom-webhook` | DICOM auto-ingest webhook (Orthanc) |
| GET | `/events/stream` | SSE event stream |
| POST | `/reports/generate` | Multi-format report generation |
| GET | `/reports/formats` | Supported export formats |
| GET | `/preview/{study_id}` | NIfTI slice preview (MP4/GIF) |
| GET | `/demo-cases` | Pre-loaded demo imaging cases |
| POST | `/protocol/optimize` | Protocol optimization advisor |
| POST | `/dose/estimate` | Radiation dose intelligence |

**Detailed Clinical Workflow Walkthrough (CT Head Hemorrhage Triage)**

1. **DICOM arrival** -- A CT head study arrives at Orthanc. If `DICOM_AUTO_INGEST` is enabled, the `/events/dicom-webhook` endpoint is triggered with the study UID, modality, and body part.

2. **Study retrieval** -- The agent fetches the DICOM series from Orthanc, converts to NIfTI format, and generates a slice preview video (8 FPS MP4).

3. **VISTA-3D segmentation** -- The NIfTI volume is sent to the VISTA-3D NIM at port 8530. The model segments 127 anatomical structures, identifying brain parenchyma, ventricles, midline structures, and extra-axial spaces.

4. **Hemorrhage detection** -- The workflow applies density thresholds and morphological analysis on the segmented volume to detect hyperdense regions consistent with acute hemorrhage. Classification includes epidural, subdural, subarachnoid, intraparenchymal, and intraventricular subtypes.

5. **RAG evidence enrichment** -- The agent searches `imaging_literature`, `imaging_guidelines`, and `imaging_findings` for evidence relevant to the detected pattern. Genomic evidence from `genomic_evidence` is queried if cross-modal is enabled (e.g., coagulation factor variants).

6. **Report generation** -- A structured radiology report is generated using the Llama-3 NIM or Claude, populated with findings, measurements, classification, and relevant literature citations. FHIR R4 DiagnosticReport format is available.

7. **Triage alert** -- If hemorrhage is detected, a high-priority SSE event is published to the event stream for downstream consumers.

**Example Query and Response**

```bash
curl -X POST http://localhost:8524/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What is the sensitivity of AI-assisted CXR for pneumothorax detection?", "top_k": 5}'
```

```json
{
  "question": "What is the sensitivity of AI-assisted CXR for pneumothorax detection?",
  "answer": "AI-assisted chest X-ray systems demonstrate sensitivity of 85-95% for pneumothorax detection, with specificity of 90-98%. FDA-cleared devices include Annalise.ai CXR (sensitivity 94.2%), Qure.ai qXR (sensitivity 91.8%), and Zebra Medical Vision (sensitivity 89.5%). Performance is highest for moderate-to-large pneumothoraces (>2 cm) and decreases for small or loculated collections. Current ACR guidelines recommend AI as a triage aid rather than a standalone diagnostic tool.",
  "evidence_count": 18,
  "collections_searched": 8,
  "search_time_ms": 1245.3,
  "nim_services_used": ["llm"]
}
```

**Cross-Agent Integration**

| Direction | Partner Agent | Trigger/Scenario |
|-----------|--------------|------------------|
| Outbound | Oncology (:8527) | Lung-RADS 4A+ nodules trigger genomic_evidence query for EGFR/ALK/ROS1/KRAS variants |
| Outbound | Cardiology (:8126) | Cardiac CT findings shared for integrated cardiovascular assessment |
| Inbound | Neurology (:8528) | Receives MRI brain requests for MS lesion tracking workflow |
| Inbound | Rare Disease (:8134) | Receives imaging assessment requests for organ-involved rare diseases |
| Bidirectional | All agents | Publishes imaging events via SSE; reads `genomic_evidence` (shared) |

**Performance Characteristics**

| Operation | Typical Latency | Notes |
|-----------|----------------|-------|
| RAG query (full) | 1,500-3,000 ms | 11-collection search + LLM synthesis |
| Evidence search (no LLM) | 400-800 ms | Vector search only |
| VISTA-3D segmentation | 30-90 sec | GPU-dependent; 127-class volumetric segmentation |
| MAISI synthetic generation | 60-180 sec | Full CT volume synthesis |
| VILA-M3 image analysis | 5-15 sec | Single image understanding |
| CXR Rapid Findings workflow | < 30 sec | End-to-end triage pipeline |
| NIfTI preview generation | 2-5 sec | MP4/GIF slice animation |

**Common Pitfalls and Troubleshooting**

- **NIM service unavailable**: If NIM containers are not running, set `IMAGING_NIM_ALLOW_MOCK_FALLBACK=true` to get structured mock responses for development. Check `GET /nim/status` for per-service availability.
- **DICOM webhook not triggering**: Verify `IMAGING_DICOM_AUTO_INGEST=true` and that Orthanc is configured with a Lua script pointing its `OnStableStudy` callback to `http://imaging-agent:8524/events/dicom-webhook`.
- **Cross-modal disabled by default**: The `IMAGING_CROSS_MODAL_ENABLED` flag defaults to `false`. Enable it to allow imaging findings to automatically query genomic evidence.
- **Large DICOM series**: Series with >500 slices may exceed the 10MB request size limit. Increase `IMAGING_MAX_REQUEST_SIZE_MB` or process via the Orthanc webhook path (which streams data).
- **Cloud vs local NIM**: When `IMAGING_NIM_MODE=cloud`, requests go to `integrate.api.nvidia.com` and require a valid `IMAGING_NVIDIA_API_KEY`. Cloud mode uses different model identifiers (e.g., `meta/llama-3.1-8b-instruct`).

---

### 3.10 Single-Cell Intelligence Agent (:8540)

**Overview and Clinical Purpose**

The Single-Cell Intelligence Agent provides cell-type annotation, tumor microenvironment (TME) profiling, drug response prediction, subclonal architecture analysis, spatial transcriptomics niche mapping, trajectory inference, ligand-receptor interaction analysis, biomarker discovery, CAR-T target validation, and treatment monitoring. Its knowledge base covers 57 cell types, 30 drugs, 75 markers, 4 TME phenotypes (hot, cold, excluded, immunosuppressive), and 4 spatial platforms (Visium, MERFISH, Xenium, CosMx). It supports 10 specialized workflows plus general RAG query.

**Advanced Configuration Options**

All environment variables use the `SC_` prefix. Key settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `SC_MILVUS_HOST` | `localhost` | Milvus vector database host |
| `SC_MILVUS_PORT` | `19530` | Milvus port |
| `SC_API_PORT` | `8540` | FastAPI REST server port |
| `SC_STREAMLIT_PORT` | `8130` | Streamlit UI port |
| `SC_LLM_MODEL` | `claude-sonnet-4-6` | Claude model for synthesis |
| `SC_ANTHROPIC_API_KEY` | (none) | Required for LLM features |
| `SC_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity |
| `SC_TOP_K_CELL_TYPES` | `50` | Max results from `sc_cell_types` per query |
| `SC_TOP_K_MARKERS` | `40` | Max results from `sc_markers` per query |
| `SC_TOP_K_SPATIAL` | `30` | Max results from `sc_spatial` per query |
| `SC_INGEST_SCHEDULE_HOURS` | `24` | Auto-ingest cycle |
| `SC_GPU_MEMORY_LIMIT_GB` | `120` | GPU memory limit for RAPIDS-based analysis |
| `SC_CELLXGENE_API_URL` | (none) | Optional CellxGene Census API endpoint |
| `SC_CROSS_AGENT_TIMEOUT` | `30` | Cross-agent HTTP timeout in seconds |

Collection weight tuning: 12 `SC_WEIGHT_*` variables control multi-collection fusion. `SC_WEIGHT_CELL_TYPES` defaults to 0.14 (highest) followed by `SC_WEIGHT_MARKERS` at 0.12, reflecting the primacy of cell identity in single-cell analysis.

**Milvus Collections (13)**

| Collection | Description |
|-----------|-------------|
| `sc_cell_types` | Cell type reference profiles and markers (57 types) |
| `sc_markers` | Marker gene signatures per cell type (75 markers) |
| `sc_literature` | Single-cell genomics literature |
| `sc_trials` | Clinical trials involving single-cell analysis |
| `sc_drugs` | Drug response signatures (GDSC/DepMap, 30 drugs) |
| `sc_pathways` | Pathway activity signatures |
| `sc_spatial` | Spatial transcriptomics references (4 platforms) |
| `sc_trajectories` | Differentiation trajectory templates |
| `sc_interactions` | Ligand-receptor pair databases (CellPhoneDB/NicheNet) |
| `sc_tme` | TME classification profiles (4 phenotypes) |
| `sc_clinical` | Clinical correlation data |
| `sc_methods` | Computational method references |
| `genomic_evidence` | Shared genomic evidence (read-only) |

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with component readiness |
| GET | `/collections` | Collection names and record counts |
| GET | `/workflows` | 10 available workflow definitions |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/v1/sc/query` | RAG Q&A across all single-cell collections |
| POST | `/v1/sc/search` | Multi-collection evidence search (no LLM) |
| POST | `/v1/sc/annotate` | Cell type annotation from marker expression |
| POST | `/v1/sc/tme-profile` | TME classification from cell proportions |
| POST | `/v1/sc/drug-response` | Drug response prediction per cell type |
| POST | `/v1/sc/cart-validate` | CAR-T target validation (on/off-tumor expression) |
| POST | `/v1/sc/spatial-niche` | Spatial niche mapping (Visium/MERFISH/Xenium/CosMx) |
| POST | `/v1/sc/trajectory` | Trajectory inference from pseudotime data |
| POST | `/v1/sc/interactions` | Ligand-receptor interaction analysis |
| POST | `/v1/sc/biomarker-discovery` | Biomarker discovery from differential expression |
| POST | `/v1/sc/subclonal` | Subclonal architecture and CNV analysis |
| POST | `/v1/sc/treatment-monitor` | Treatment response monitoring over time |
| POST | `/v1/sc/workflow/{type}` | Generic workflow dispatch (10 types) |
| POST | `/v1/sc/integrated-assessment` | Cross-agent multi-agent assessment |
| GET | `/v1/sc/cell-types` | Reference catalog of 57 cell types |
| GET | `/v1/sc/markers` | Reference catalog of 75 markers |
| GET | `/v1/sc/tme-classes` | TME classification definitions (4 phenotypes) |
| GET | `/v1/sc/spatial-platforms` | Supported spatial platforms |
| POST | `/v1/reports/generate` | Multi-format report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | Server-sent event stream |

**Detailed Clinical Workflow Walkthrough (TME-Guided Immunotherapy Selection)**

1. **Expression data submission** -- Clinician submits a gene expression marker panel or cell type proportions from a tumor biopsy via the `/v1/sc/tme-profile` endpoint. Input can be raw marker values (e.g., `CD3D: 5.2, CD8A: 4.1, GZMB: 3.8`) or pre-computed cell type proportions.

2. **Cell type annotation** -- If raw markers are provided, the agent performs multi-strategy annotation: reference-based matching against the 57-cell-type atlas, marker-based scoring using the 75-marker database, and LLM-augmented classification for ambiguous profiles. Confidence levels (high/medium/low) are assigned.

3. **TME classification** -- Cell type proportions are analyzed against the 4 TME phenotype profiles:
   - **Hot** (immune-inflamed): high CD8+ T cells (>15%), high PD-L1 expression, active IFN-gamma signaling
   - **Cold** (immune-desert): low immune infiltrate (<5%), absent T cell signatures
   - **Excluded** (immune-excluded): immune cells at tumor periphery, high stromal content
   - **Immunosuppressive**: high Treg (>8%), high M2 macrophage, elevated TGF-beta/IL-10

4. **Drug response prediction** -- Based on the TME phenotype and cell-type-resolved expression, the agent queries `sc_drugs` for GDSC/DepMap sensitivity signatures across 30 drugs. Checkpoint inhibitors (pembrolizumab, nivolumab, atezolizumab) are ranked for hot TMEs; combination strategies are suggested for cold/excluded TMEs.

5. **Subclonal architecture** -- If CNV data is available, the agent performs subclonal detection to identify resistant subpopulations and assess antigen escape risk. Clones with loss of heterozygosity at HLA loci are flagged.

6. **Spatial context** -- For spatial transcriptomics data (Visium, MERFISH, Xenium, CosMx), the agent maps immune cell niches relative to tumor boundaries, quantifying immune exclusion zones and tertiary lymphoid structures.

7. **Report generation** -- Results exported as Markdown, JSON, PDF, or FHIR R4.

**Example Query and Response**

```bash
curl -X POST http://localhost:8540/v1/sc/tme-profile \
  -H "Content-Type: application/json" \
  -d '{"cell_type_proportions": {"CD8_T_cell": 0.15, "Macrophage": 0.25, "Treg": 0.08, "Fibroblast": 0.30, "Tumor": 0.18, "NK_cell": 0.04}}'
```

```json
{
  "status": "completed",
  "tme_classification": "immunosuppressive",
  "confidence": "high",
  "rationale": "Elevated Treg proportion (8%) combined with high macrophage content (25%) suggests immunosuppressive microenvironment. CD8+ T cell infiltration (15%) is present but likely functionally exhausted based on the Treg:CD8 ratio.",
  "immunotherapy_recommendations": [
    {
      "drug": "ipilimumab + nivolumab",
      "rationale": "Combination anti-CTLA-4/PD-1 to overcome Treg-mediated suppression",
      "evidence_level": "Level 1A (CheckMate 067)"
    },
    {
      "drug": "anti-CCR4 (mogamulizumab)",
      "rationale": "Treg depletion strategy for immunosuppressive TME",
      "evidence_level": "Level 2B (Phase II data)"
    }
  ],
  "escape_risk": "moderate",
  "collections_searched": 8,
  "processing_time_ms": 1850.3
}
```

**Cross-Agent Integration**

| Direction | Partner Agent | Trigger/Scenario |
|-----------|--------------|------------------|
| Outbound | CAR-T (:8522) | Sends antigen escape alerts when off-tumor expression detected in CAR-T target validation |
| Outbound | Oncology (:8527) | Provides TME context for immunotherapy selection and combination strategy |
| Outbound | Biomarker (:8529) | Shares cell-type-resolved biomarker panels for panel enrichment |
| Outbound | Drug Discovery | Queries compound integration for cell-type-selective drug candidates |
| Inbound | Imaging (:8524) | Receives spatial-imaging correlation requests |
| Inbound | Pharmacogenomics (:8107) | Receives requests for cell-type-resolved drug sensitivity data |
| Bidirectional | All agents | Reads `genomic_evidence` (shared); publishes events via SSE |

**Performance Characteristics**

| Operation | Typical Latency | Notes |
|-----------|----------------|-------|
| Cell type annotation (marker panel) | 500-1,200 ms | Multi-strategy: reference + marker + LLM |
| TME profiling | 800-1,500 ms | Proportion analysis + evidence retrieval |
| Drug response prediction | 600-1,200 ms | GDSC/DepMap signature matching |
| CAR-T target validation | 400-800 ms | Expression ratio analysis |
| Spatial niche mapping | 1,000-2,500 ms | Platform-specific spatial analysis |
| Full RAG query | 1,500-3,500 ms | 13-collection search + LLM synthesis |
| Subclonal architecture | 2,000-4,000 ms | CNV detection + clone assignment |

**Common Pitfalls and Troubleshooting**

- **Cell type proportion normalization**: Proportions submitted to `/v1/sc/tme-profile` should sum to approximately 1.0. The agent normalizes internally, but significantly unbalanced inputs (e.g., sum > 2.0) may indicate double-counting.
- **Marker gene naming**: Gene symbols must follow HGNC nomenclature (e.g., `CD3D` not `CD3`). The agent attempts alias resolution via `ENTITY_ALIASES` but non-standard names may fail silently.
- **TME classification edge cases**: Tumors with mixed hot/cold regions (spatially heterogeneous) may be classified inconsistently depending on the biopsy site. Spatial transcriptomics data is recommended for heterogeneous tumors.
- **GDSC/DepMap coverage**: Drug response predictions are limited to the 30 drugs in the sensitivity database. Novel agents or combination regimens fall back to LLM-based inference from literature.
- **GPU memory for RAPIDS**: If `SC_GPU_MEMORY_LIMIT_GB` is set too low for large expression matrices, RAPIDS operations will fall back to CPU, significantly increasing latency.

---

### 3.11 Clinical Trial Intelligence Agent (:8538)

**Overview and Clinical Purpose**

The Clinical Trial Intelligence Agent provides AI-driven support across the full clinical trial lifecycle: protocol optimization, patient-trial matching, site selection with 7-factor scoring, eligibility optimization, adaptive design evaluation, safety signal detection using PRR/ROR disproportionality analysis, regulatory document generation, competitive intelligence, diversity assessment, and decentralized trial planning. It supports 10 specialized workflows plus general RAG query capability, and its knowledge base includes 40 landmark trial references.

**Advanced Configuration Options**

All environment variables use the `TRIAL_` prefix. Key settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `TRIAL_MILVUS_HOST` | `localhost` | Milvus vector database host |
| `TRIAL_MILVUS_PORT` | `19530` | Milvus port |
| `TRIAL_API_PORT` | `8538` | FastAPI REST server port |
| `TRIAL_STREAMLIT_PORT` | `8128` | Streamlit UI port |
| `TRIAL_LLM_MODEL` | `claude-sonnet-4-6` | Claude model for synthesis |
| `TRIAL_ANTHROPIC_API_KEY` | (none) | Required for LLM features |
| `TRIAL_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity |
| `TRIAL_TOP_K_PER_COLLECTION` | `5` | Default results per collection |
| `TRIAL_INGEST_SCHEDULE_HOURS` | `24` | Auto-ingest cycle |
| `TRIAL_API_KEY` | (empty) | Set to require X-API-Key header |
| `TRIAL_CROSS_AGENT_TIMEOUT` | `30` | Cross-agent HTTP timeout in seconds |
| `TRIAL_CLINICALTRIALS_API_KEY` | (none) | Optional ClinicalTrials.gov API key |

The agent connects to 8 partner agents for cross-agent integration, with dedicated URL variables: `TRIAL_ONCOLOGY_AGENT_URL`, `TRIAL_PGX_AGENT_URL`, `TRIAL_CARDIOLOGY_AGENT_URL`, `TRIAL_BIOMARKER_AGENT_URL`, `TRIAL_RARE_DISEASE_AGENT_URL`, `TRIAL_NEUROLOGY_AGENT_URL`, `TRIAL_SINGLE_CELL_AGENT_URL`, `TRIAL_IMAGING_AGENT_URL`.

**Milvus Collections (14)**

| Collection | Description |
|-----------|-------------|
| `trial_protocols` | Protocol designs, amendments, SOAs |
| `trial_eligibility` | Inclusion/exclusion criteria |
| `trial_endpoints` | Primary, secondary, exploratory endpoints |
| `trial_sites` | Site feasibility, enrollment history |
| `trial_investigators` | Investigator profiles, experience |
| `trial_results` | Published trial results, CSRs (40 landmark trials) |
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
| Site Selection | 7-factor feasibility scoring, enrollment forecasting, diversity metrics |
| Eligibility Optimization | Population impact modeling, competitor benchmarking |
| Adaptive Design | Bayesian interim analysis, dose-response, futility assessment |
| Safety Signal | Disproportionality analysis (PRR/ROR), causality assessment |
| Regulatory Docs | IND, CSR, briefing doc generation (FDA/EMA/PMDA) |
| Competitive Intel | Landscape analysis, enrollment race tracking |
| Diversity Assessment | FDA diversity guidance compliance, gap analysis |
| Decentralized Planning | DCT component feasibility, hybrid model design |

**Complete API Endpoint Table**

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with component readiness |
| GET | `/collections` | Collection names and record counts |
| GET | `/workflows` | 10 workflow definitions + general query |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/v1/trial/query` | RAG Q&A across all trial collections |
| POST | `/v1/trial/search` | Multi-collection evidence search (no LLM) |
| POST | `/v1/trial/workflow/protocol_design/run` | Protocol design complexity scoring |
| POST | `/v1/trial/workflow/patient_matching/run` | Patient-trial matching with genomics |
| POST | `/v1/trial/workflow/site_selection/run` | 7-factor site feasibility scoring |
| POST | `/v1/trial/workflow/eligibility_optimization/run` | Population impact modeling |
| POST | `/v1/trial/workflow/adaptive_design/run` | Bayesian adaptive design evaluation |
| POST | `/v1/trial/workflow/safety_signal/run` | PRR/ROR disproportionality analysis |
| POST | `/v1/trial/workflow/regulatory_docs/run` | Regulatory document generation |
| POST | `/v1/trial/workflow/competitive_intel/run` | Competitive landscape analysis |
| POST | `/v1/trial/workflow/diversity_assessment/run` | FDA diversity guidance compliance |
| POST | `/v1/trial/workflow/decentralized_planning/run` | DCT feasibility assessment |
| POST | `/v1/trial/integrated-assessment` | Cross-agent multi-agent assessment |
| GET | `/v1/trial/therapeutic-areas` | Reference catalog of therapeutic areas |
| GET | `/v1/trial/phases` | Trial phase definitions |
| GET | `/v1/trial/landmark-trials` | 40 landmark trial references |
| GET | `/v1/trial/regulatory-agencies` | Supported regulatory agencies (FDA/EMA/PMDA) |
| POST | `/v1/reports/generate` | Multi-format report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | Server-sent event stream |

**Detailed Clinical Workflow Walkthrough (Safety Signal Detection)**

1. **Adverse event submission** -- Safety team submits adverse event data via `/v1/trial/workflow/safety_signal/run`. Input includes the drug name, a list of adverse events with counts (exposed and unexposed), and optional background rate data.

2. **Disproportionality analysis** -- The agent computes two standard pharmacovigilance metrics:
   - **PRR (Proportional Reporting Ratio)**: ratio of the proportion of a specific AE for the drug vs. all other drugs in the database. PRR > 2 with chi-squared > 4 and N >= 3 triggers a signal.
   - **ROR (Reporting Odds Ratio)**: odds ratio of the AE in exposed vs. unexposed populations. The 95% confidence interval lower bound > 1 indicates a statistically significant signal.

3. **Historical comparator retrieval** -- The agent searches `trial_safety` for historical AE rates from similar drug classes and indications. Published DSMB reports and FDA safety communications are retrieved from `trial_regulatory`.

4. **Causality assessment** -- Evidence from `trial_literature` and `trial_results` is used to assess temporal relationship, dose-response, biological plausibility, and consistency across studies. The LLM synthesizes a causality narrative referencing the WHO-UMC and Naranjo scales.

5. **Benchmark against landmark trials** -- The detected signal is compared against AE profiles from relevant landmark trials in `trial_results` (40 stored). For example, hepatotoxicity in a checkpoint inhibitor trial would be benchmarked against CheckMate 067 and KEYNOTE 024 safety data.

6. **Report generation** -- A structured safety signal report is generated with PRR/ROR calculations, confidence intervals, historical comparator data, causality assessment, and recommended actions (continue monitoring, dose modification, or study hold).

**Example Query and Response**

```bash
curl -X POST http://localhost:8538/v1/trial/workflow/safety_signal/run \
  -H "Content-Type: application/json" \
  -d '{
    "drug": "experimental_agent_X",
    "adverse_events": [
      {"event": "hepatotoxicity", "count": 8, "total_exposed": 200},
      {"event": "rash", "count": 15, "total_exposed": 200},
      {"event": "fatigue", "count": 45, "total_exposed": 200}
    ],
    "comparator_class": "checkpoint_inhibitor"
  }'
```

```json
{
  "status": "completed",
  "workflow": "safety_signal",
  "signals_detected": [
    {
      "event": "hepatotoxicity",
      "prr": 3.2,
      "prr_ci_lower": 1.8,
      "ror": 3.5,
      "ror_ci_lower": 1.6,
      "signal_strength": "moderate",
      "historical_rate_class": "2-5% for checkpoint inhibitors",
      "observed_rate": "4.0%",
      "assessment": "Signal detected. Observed hepatotoxicity rate (4.0%) is within the expected range for checkpoint inhibitors (2-5%) but PRR exceeds threshold. Recommend enhanced liver function monitoring at 2-week intervals.",
      "causality_score": "possible (Naranjo: 4)"
    },
    {
      "event": "rash",
      "prr": 1.4,
      "signal_strength": "none",
      "assessment": "No signal. Rash rate (7.5%) is below the expected range for checkpoint inhibitors (10-20%)."
    }
  ],
  "landmark_comparators": ["CheckMate 067", "KEYNOTE 024", "IMpower 150"],
  "evidence_sources": 28,
  "processing_time_ms": 3120.7
}
```

**7-Factor Site Selection Scoring**

The site selection workflow evaluates investigative sites using a composite score across 7 factors:

| Factor | Weight | Data Sources |
|--------|--------|-------------|
| Disease prevalence in catchment area | 0.20 | `trial_sites`, census data |
| Prior enrollment performance | 0.20 | `trial_investigators`, historical enrollment |
| Investigator experience (publications, prior trials) | 0.15 | `trial_investigators`, `trial_literature` |
| Regulatory readiness (IRB turnaround, startup time) | 0.15 | `trial_sites` |
| Diversity index (demographic representation) | 0.10 | `trial_sites`, FDA diversity guidance |
| Competitor trial burden (active competing studies) | 0.10 | `trial_protocols`, competitive intelligence |
| Infrastructure score (lab, imaging, pharmacy) | 0.10 | `trial_sites` |

**Cross-Agent Integration**

| Direction | Partner Agent | Trigger/Scenario |
|-----------|--------------|------------------|
| Inbound | Oncology (:8527) | Receives trial matching requests when actionable variants found |
| Inbound | Rare Disease (:8134) | Receives rare disease trial eligibility queries |
| Inbound | Biomarker (:8529) | Receives biomarker-driven trial design queries |
| Inbound | PGx (:8107) | Receives queries for PGx-guided clinical trials |
| Outbound | Oncology (:8527) | Queries tumor profiling for patient eligibility context |
| Outbound | PGx (:8107) | Queries metabolizer status for trial-specific dosing arms |
| Outbound | Cardiology (:8126) | Queries cardiac safety context for cardiovascular endpoint trials |
| Outbound | Biomarker (:8529) | Queries biomarker interpretation for eligibility criteria |
| Bidirectional | All agents | Reads `genomic_evidence` (shared); publishes events via SSE |

**Performance Characteristics**

| Operation | Typical Latency | Notes |
|-----------|----------------|-------|
| RAG query | 1,500-3,000 ms | 14-collection search + LLM synthesis |
| Patient-trial matching | 2,000-4,000 ms | Multi-collection search + eligibility scoring |
| Protocol design review | 2,500-4,500 ms | Complexity scoring + benchmark comparison |
| Safety signal (PRR/ROR) | 1,500-3,500 ms | Statistical calculation + evidence retrieval |
| Site selection (7-factor) | 3,000-5,000 ms | Multi-factor scoring + enrollment forecasting |
| Regulatory document generation | 5,000-10,000 ms | Template population + LLM narrative |
| Competitive intelligence | 2,000-4,000 ms | Landscape analysis across `trial_protocols` |

**Common Pitfalls and Troubleshooting**

- **Patient matching requires complete profiles**: The patient-trial matching workflow performs best with age, diagnosis, biomarker status, prior therapy history, and ECOG performance status. Partial profiles produce fewer matches with lower confidence.
- **PRR/ROR minimum counts**: Disproportionality analysis requires at least 3 events (N >= 3) to generate a statistically meaningful signal. Submitting events with count < 3 will return "insufficient data" rather than a false signal.
- **Landmark trial coverage**: The 40 landmark trials are curated references. New landmark trials (e.g., recently published Phase III results) require manual addition via `scripts/seed_knowledge.py`.
- **ClinicalTrials.gov rate limiting**: If `TRIAL_CLINICALTRIALS_API_KEY` is not set, the ClinicalTrials.gov v2 API rate-limits to 3 requests/second. Set the key for production workloads.
- **Regulatory document templates**: Generated IND/CSR documents are templates requiring clinical review. They auto-populate statistical sections, safety summaries, and study schemas but do not replace regulatory affairs expertise.
- **Cross-agent integration breadth**: This agent connects to 8 partner agents (more than any other agent). If multiple partner agents are unavailable, the integrated assessment still returns partial results from available agents.

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
