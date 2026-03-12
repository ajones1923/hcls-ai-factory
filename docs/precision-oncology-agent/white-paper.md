# Multi-Collection RAG Architecture for Molecular Tumor Board Intelligence: Evidence-Based Therapy Ranking, Trial Matching, and Resistance-Aware Recommendations

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**License:** Apache 2.0

---

## Abstract

Molecular tumor boards (MTBs) require synthesis of heterogeneous
evidence spanning somatic variant databases, clinical guidelines,
trial registries, pharmacogenomic literature, and resistance
mechanisms -- a task that exceeds the cognitive bandwidth of any
single clinician and resists consolidation by conventional
informatics pipelines. We present the Precision Oncology
Intelligence Agent, a retrieval-augmented generation (RAG) system
that federates eleven Milvus vector collections under a weighted,
multi-strategy search planner to deliver therapy-ranked,
resistance-aware, trial-matched decision support for 26 cancer
types and over 40 actionable gene targets. The system ingests
patient VCF files annotated by SnpEff, VEP, or GENE/GENEINFO
pipelines; constructs per-case knowledge graphs linking variants
to biomarkers, therapies, resistance pathways, and open trials;
ranks candidate therapies using an AMP/ASCO/CAP-aligned evidence
tier (levels A through E); matches patients to clinical trials
via a composite score blending biomarker relevance (0.40), semantic
similarity (0.25), trial phase (0.20), and recruitment status
(0.15); and exports MTB packets as Markdown, JSON, PDF (via
ReportLab), and FHIR R4 DiagnosticReport bundles. All embeddings
use 384-dimensional BAAI/bge-small-en-v1.5 vectors indexed with
IVF_FLAT and COSINE similarity. Synthesis is performed by Claude
Sonnet 4.6. The complete system comprises 66 Python files totaling
approximately 20,490 lines of code, validated by 556 tests across
10 test modules, and runs on a single NVIDIA DGX Spark at a
hardware cost of $3,999.

---

## 1. Introduction

### 1.1 The Molecular Tumor Board Challenge

Precision oncology has shifted the standard of care from
histology-driven protocols to biomarker-driven therapy selection.
Molecular tumor boards -- multidisciplinary panels that review
next-generation sequencing (NGS) results and recommend targeted
treatments -- are now established at most academic cancer centers.
Yet the informatics burden on these panels is acute and growing.

A single tumor genome may harbor hundreds of somatic variants.
For each variant, the MTB must consult:

- **Variant databases** (CIViC, OncoKB) cataloging the clinical
  significance of somatic and germline alterations.
- **Clinical guidelines** (NCCN, ESMO, ASCO) encoding
  category-of-evidence recommendations for biomarker-therapy
  pairings.
- **Trial registries** (ClinicalTrials.gov) listing thousands of
  open studies with complex eligibility criteria.
- **Primary literature** (PubMed, PMC) where new evidence appears
  daily.
- **Resistance literature** documenting acquired and intrinsic
  mechanisms that render first-line therapies ineffective.

This information is fragmented across incompatible formats, updated
on different cadences, and scaled beyond manual review. A typical
MTB case review lasts 10-15 minutes, during which oncologists,
pathologists, geneticists, and pharmacists must reach consensus on
therapy, trials, and follow-up -- often with incomplete evidence.

### 1.2 Limitations of Existing Approaches

Current decision support tools fall into two categories. Rule-based
annotation engines (e.g., OpenCRAVAT, Ensembl VEP) attach known
annotations to variants but do not synthesize evidence across
sources or adapt to clinical context. Commercial platforms (e.g.,
Foundation Medicine, Tempus) provide integrated reports but operate
as closed systems with per-patient licensing costs that limit
accessibility. Neither approach provides real-time, multi-source
evidence synthesis with resistance awareness, trial matching, and
exportable MTB packets in an open-source, locally deployable
package.

### 1.3 Contribution

This paper describes the Precision Oncology Intelligence Agent, a
RAG-powered MTB decision support system that:

1. Federates 11 vector collections covering the full evidence
   landscape.
2. Implements weighted, strategy-adaptive search across all
   collections.
3. Ranks therapies with resistance-aware, evidence-tiered scoring.
4. Matches patients to clinical trials using hybrid
   deterministic-semantic composite scores.
5. Generates exportable MTB packets in four formats including
   FHIR R4.
6. Runs on a single NVIDIA DGX Spark ($3,999), making
   institutional-grade decision support accessible to community
   oncology practices.
7. Is fully open-source under the Apache 2.0 license.

---

## 2. System Architecture

### 2.1 Deployment Context

The agent operates within the HCLS AI Factory, a three-stage
precision medicine pipeline transforming raw patient DNA (FASTQ)
into drug candidates in under five hours:

1. **Genomics Pipeline** -- Parabricks/DeepVariant/BWA-MEM2,
   producing annotated VCF files.
2. **RAG/Chat Pipeline** -- Milvus vector database with Claude AI
   for variant interpretation.
3. **Drug Discovery Pipeline** -- BioNeMo MolMIM/DiffDock/RDKit
   for molecular docking and lead optimization.

The oncology agent sits at the intersection of stages 1 and 2,
consuming VCF output from the genomics pipeline and producing
structured therapy recommendations that feed the drug discovery
stage.

### 2.2 Multi-Collection Vector Store

The system organizes knowledge into 11 Milvus collections, each
serving a distinct evidentiary role. All collections use
384-dimensional embeddings from BAAI/bge-small-en-v1.5, indexed
with IVF_FLAT, searched by COSINE similarity:

| # | Collection | Weight | Source Domain |
|---|---|---|---|
| 1 | `onco_variants` | 0.18 | CIViC/OncoKB actionable variants |
| 2 | `onco_literature` | 0.16 | PubMed/PMC literature chunks |
| 3 | `onco_therapies` | 0.14 | FDA-approved therapies + mechanisms |
| 4 | `onco_guidelines` | 0.12 | NCCN/ESMO/ASCO clinical guidelines |
| 5 | `onco_trials` | 0.10 | ClinicalTrials.gov records |
| 6 | `onco_biomarkers` | 0.08 | Predictive/prognostic biomarkers |
| 7 | `onco_resistance` | 0.07 | Resistance mechanisms/bypass paths |
| 8 | `onco_pathways` | 0.06 | Oncogenic signaling (13 pathways) |
| 9 | `onco_outcomes` | 0.04 | Treatment outcomes (RECIST scale) |
| 10 | `onco_cases` | 0.02 | Patient case snapshots |
| 11 | `genomic_evidence` | 0.03 | Read-only shared from Stage 1 |

The weight distribution reflects clinical priority: actionable
variant evidence and peer-reviewed literature receive the highest
weights, while historical case snapshots serve as supplementary
context. The `genomic_evidence` collection is read-only and shared
with other agents in the HCLS AI Factory ecosystem.

### 2.3 Embedding Strategy

BAAI/bge-small-en-v1.5 (33M parameters, 384-dimensional output)
was selected for biomedical retrieval competence, deployment
efficiency on constrained hardware, and compatibility with IVF_FLAT
indexing for sub-second search latency across millions of vectors.

### 2.4 Agent Workflow

The agent follows a plan-search-evaluate-synthesize loop:

```
plan -> search -> evaluate -> (retry if insufficient) -> synthesize
```

**SearchPlan.** Given a patient query or VCF-derived variant set,
the planner identifies relevant genes, cancer types, and topics,
then selects a search strategy:

- **Broad**: queries all 11 collections for exploratory or complex
  cases.
- **Targeted**: focuses on the highest-weighted collections for
  well-characterized variants.
- **Comparative**: searches therapies, resistance, and outcomes
  collections for head-to-head therapy comparison.

**Evidence Evaluation.** Retrieved evidence is classified as:

- **Sufficient**: 3 or more hits from 2 or more distinct
  collections.
- **Partial**: some relevant hits below the sufficiency threshold.
- **Insufficient**: fewer than 3 total hits.

If evidence is insufficient, the agent retries with broadened
queries or relaxed similarity thresholds, up to MAX_RETRIES = 2.

### 2.5 Multi-Module Architecture

The agent comprises five core modules:

1. **RAG Engine** -- collection access, query embedding, weighted
   retrieval, and evidence assembly.
2. **Therapy Ranker** -- evidence tiering, resistance penalties,
   biomarker overrides, and ranked output.
3. **Trial Matcher** -- composite scoring, biomarker matching,
   semantic search, and trial ranking.
4. **Case Manager** -- VCF ingestion, patient profile construction,
   and MTB packet orchestration.
5. **Knowledge Graph** -- structured relationships between genes,
   variants, therapies, resistance mechanisms, pathways, and
   biomarkers.

### 2.6 Cross-Agent Integration

The agent participates in the HCLS AI Factory event system:

- **ONCOLOGY_CASE_CREATED** -- emitted when a new patient case is
  registered, triggering downstream drug discovery docking.
- **THERAPY_RANKED** -- emitted when therapy ranking completes,
  notifying the dashboard and subscribed agents.

---

## 3. Knowledge Graph

### 3.1 Actionable Gene Targets

The knowledge graph encodes over 40 actionable gene targets curated
from CIViC, OncoKB, and NCCN guidelines. Each entry includes:

- Gene symbol and HGNC identifier.
- Associated cancer types (from the 26 supported types).
- Known actionable variants (e.g., EGFR L858R, BRAF V600E).
- Linked therapies with evidence levels.
- Known resistance mechanisms and bypass pathways.

### 3.2 Therapy Map

The therapy map links each FDA-approved targeted therapy to its
molecular targets, approved indications, and mechanism of action.
Entries are sourced from FDA label data and supplemented by
guideline-level recommendations from NCCN and ESMO.

### 3.3 Resistance Map

The resistance map catalogs known mechanisms of acquired and
intrinsic resistance for each gene-therapy pair:

- **On-target resistance**: secondary mutations in the drug target
  (e.g., EGFR T790M conferring resistance to first-generation EGFR
  TKIs).
- **Bypass pathway activation**: upregulation of alternative
  signaling (e.g., MET amplification bypassing EGFR inhibition).
- **Downstream reactivation**: mutations in downstream effectors
  (e.g., KRAS mutations reactivating MAPK signaling under BRAF
  inhibition).

### 3.4 Pathway Map

The pathway map models 13 oncogenic signaling pathways: MAPK/ERK,
PI3K/AKT/mTOR, WNT/beta-catenin, Notch, Hedgehog, JAK/STAT, DNA
damage repair (HRD), cell cycle (CDK4/6-RB), apoptosis (BCL-2),
angiogenesis (VEGF/VEGFR), chromatin remodeling, immune checkpoint
(PD-1/PD-L1/CTLA-4), and receptor tyrosine kinase
(EGFR/ALK/ROS1/RET/NTRK). Each node links to member genes,
oncogenic alterations, and pathway-targeted therapies.

### 3.5 Biomarker Panels

Structured biomarker panels encode deterministic therapy-selection
rules aligned with NCCN and FDA companion diagnostic approvals:

| Biomarker | Therapies |
|---|---|
| MSI-H / TMB-H | pembrolizumab, nivolumab |
| HRD+ | olaparib, rucaparib, niraparib, talazoparib |
| PD-L1 TPS >= 50% | pembrolizumab |
| NTRK fusion | larotrectinib, entrectinib |
| EGFR+ | osimertinib, erlotinib, gefitinib |
| BRAF V600E | vemurafenib, dabrafenib, encorafenib |
| ALK+ | crizotinib, alectinib |
| KRAS G12C | sotorasib, adagrasib |
| RET fusion | selpercatinib, pralsetinib |

These rules serve as hard constraints in the therapy ranking
engine, ensuring guideline-concordant recommendations are never
suppressed by lower-confidence RAG-derived evidence.

---

## 4. Therapy Ranking Engine

### 4.1 Evidence-Level Tiering

Therapies are scored using a tier system aligned with the
AMP/ASCO/CAP joint consensus guidelines for somatic variant
classification:

| Tier | Level | Description |
|---|---|---|
| I | A | FDA-approved therapy, companion diagnostic approved |
| I | B | Well-powered studies, consensus guideline support |
| II | C | FDA-approved for different tumor type, case series |
| II | D | Preclinical data, case reports, biological rationale |
| III | E | Investigational, conflicting evidence |

Each therapy-variant pair receives a tier assignment based on the
strongest evidence retrieved across collections. The tier drives
the base score for therapy ranking.

### 4.2 Resistance Awareness

The ranking engine applies resistance penalties at three levels:

1. **Known resistance mutation detected**: if the patient's VCF
   contains a variant cataloged as a resistance mechanism for a
   candidate therapy, that therapy is flagged and penalized.
2. **Pathway-level bypass risk**: if alterations in bypass pathway
   genes are detected, therapies targeting the primary pathway
   receive a risk annotation.
3. **Historical resistance prevalence**: population-level resistance
   rates from the outcomes collection inform prior probabilities of
   treatment failure.

### 4.3 Biomarker-Driven Recommendations

Biomarker rules from Section 3.5 are applied as deterministic
overrides. When a patient's molecular profile matches a biomarker
panel entry, associated therapies are promoted to the top of the
ranked list regardless of RAG-derived scores, ensuring
standard-of-care recommendations are never ranked below
investigational options.

### 4.4 Ranking Output

The final ranked therapy list includes, for each candidate:

- Drug name and mechanism of action.
- Evidence level (A-E) and supporting sources.
- Resistance flags and bypass pathway warnings.
- Biomarker match status.
- Relevant clinical trials (cross-referenced from the trial
  matcher).

---

## 5. Clinical Trial Matching

### 5.1 Hybrid Scoring

The trial matcher blends deterministic biomarker matching with
semantic similarity, addressing a limitation of pure semantic
search: eligibility criteria contain structured fields requiring
exact matching alongside free-text benefiting from semantic
understanding. The composite score is computed as:

The composite trial match score is computed as:

```
score = (0.40 * biomarker_score) + (0.25 * semantic_score)
      + (0.20 * phase_score)     + (0.15 * status_score)
```

| Component | Weight | Description |
|---|---|---|
| Biomarker match | 0.40 | Exact match to trial molecular eligibility |
| Semantic score | 0.25 | COSINE similarity of profile to trial text |
| Phase score | 0.20 | Phase 3 > Phase 2 > Phase 1 |
| Status score | 0.15 | Recruiting > Not yet recruiting > Active |

### 5.2 Matching Components

The **biomarker component** performs deterministic matching between
patient variants/biomarkers and trial molecular eligibility:
specific gene mutations, gene-level alterations, biomarker status
(MSI-H, TMB-H, PD-L1), and fusion events. The **semantic
component** embeds the patient's clinical summary and computes
COSINE similarity against trial descriptions in `onco_trials`,
capturing nuanced signals like "progression on prior immunotherapy"
that resist rule-based extraction.

---

## 6. Case Management

### 6.1 VCF Parsing

The case manager ingests annotated VCF files supporting three
annotation formats:

- **SnpEff**: parses ANN fields for gene, variant effect, impact,
  and HGVS notation.
- **VEP**: parses CSQ fields for gene, consequence, SIFT/PolyPhen
  predictions, and clinical significance.
- **GENE/GENEINFO**: parses simpler formats providing gene symbol
  and basic variant information.

The parser extracts chromosome, position, ref/alt alleles, quality
scores, genotype, and all annotation fields. Variants are filtered
by impact severity, actionable target membership, and pathogenicity
classifications.

### 6.2 Patient Profile and MTB Packet Generation

From the parsed VCF and supplementary clinical data, the case
manager constructs a structured profile (detected variants,
biomarker panel results, cancer type, prior treatments, family
history) and orchestrates the MTB packet workflow: parse VCF,
retrieve evidence via RAG engine, run therapy ranker, execute trial
matcher, assemble structured packet, and export in requested
format(s).

---

## 7. Export and Interoperability

**Markdown** -- default human-readable format with patient summary,
ranked therapies, matched trials, resistance warnings, and PubMed
references. **JSON** -- machine-readable structured output for EHR
integration and downstream analytics. **PDF** -- publication-quality
documents via ReportLab with formatted tables and institutional
branding, for official medical records. **FHIR R4** --
standards-compliant DiagnosticReport bundle containing Observation
(variants/biomarkers), MedicationRequest (therapies), and
ResearchStudy (trials) resources, enabling interoperability with
any FHIR-capable EHR per the HL7 Genomics Implementation Guide.

---

## 8. Implementation

### 8.1 Codebase Metrics

The system is implemented in Python across 66 source files totaling
approximately 20,490 lines of code:

| Module | Responsibility |
|---|---|
| RAG Engine | Collection management, query routing, weighted retrieval |
| Therapy Ranker | Evidence tiering, resistance penalties, biomarker overrides |
| Trial Matcher | Composite scoring, biomarker matching, trial ranking |
| Case Manager | VCF parsing, profile construction, MTB orchestration |
| Knowledge Graph | Gene targets, therapy/resistance/pathway maps |
| Export Module | Markdown, JSON, PDF, FHIR R4 generation |
| Event Integration | Cross-agent event emission and consumption |
| API Layer | FastAPI endpoints for external access |

### 8.2 Test Coverage

556 tests across 10 test files cover VCF parsing, RAG engine
integration, therapy ranking, trial matching scores, FHIR R4
validation, export correctness, knowledge graph integrity, and
end-to-end workflows. All 556 tests pass.

### 8.3 Dependencies and Hardware

Key dependencies: Milvus (vector DB), BAAI/bge-small-en-v1.5
(embeddings), Claude Sonnet 4.6 (LLM synthesis), ReportLab (PDF),
FastAPI (API layer), RDKit (optional). The complete system runs on
a single NVIDIA DGX Spark ($3,999) with sufficient GPU memory for
concurrent embedding and vector search across all 11 collections.

---

## 9. Clinical Validation

Four representative patient profiles demonstrate system capability
across common precision oncology scenarios, exercising the full
pipeline from VCF ingestion through MTB packet generation.

### 9.1 Case 1: NSCLC with EGFR L858R

**Profile:** 62-year-old female, stage IIIB NSCLC adenocarcinoma.
EGFR L858R (exon 21) missense mutation. PD-L1 TPS 30%, TMB-low.

**Output:** Tier I-A therapies: osimertinib (preferred), erlotinib,
gefitinib -- all FDA-approved for EGFR-mutant NSCLC. Resistance
flag: T790M not detected but noted as primary resistance mechanism
to first-generation TKIs; osimertinib ranked highest due to T790M
activity. Trial matches: 3 Phase 3 trials for EGFR-mutant NSCLC
including anti-angiogenic combinations. Pathway context: EGFR
signaling through MAPK/ERK and PI3K/AKT/mTOR annotated; MET
amplification bypass risk documented.

### 9.2 Case 2: Breast Cancer with BRCA1 Pathogenic Variant

**Profile:** 45-year-old female, stage II triple-negative breast
cancer (TNBC). BRCA1 c.68_69delAG (185delAG) frameshift. HRD+.

**Output:** Tier I-A: olaparib, talazoparib -- FDA-approved PARP
inhibitors for gBRCA/HER2-negative breast cancer. Tier I-B:
rucaparib, niraparib (NCCN Category 2A). Biomarker override: HRD+
triggers PARP inhibitor promotion. Trial matches: 2 Phase 2 PARP +
checkpoint inhibitor trials for TNBC. Resistance context: BRCA1
reversion mutations documented as primary PARP inhibitor resistance
mechanism.

### 9.3 Case 3: Colorectal Cancer with KRAS G12C

**Profile:** 58-year-old male, stage IV colorectal adenocarcinoma
with liver metastases. KRAS G12C. MSS, TMB-low.

**Output:** Tier I-A: sotorasib (FDA-approved for KRAS G12C NSCLC,
emerging CRC data). Tier II-C: adagrasib (cross-indication
evidence). Resistance flag: KRAS and MET amplification as emerging
resistance mechanisms; RAS-MAPK reactivation risk. Trial matches:
4 open trials including sotorasib + panitumumab combinations.
Negative biomarker: MSS contraindicates pembrolizumab; agent
correctly excludes checkpoint immunotherapy.

### 9.4 Case 4: Melanoma with BRAF V600E

**Profile:** 51-year-old male, stage IIIC cutaneous melanoma. BRAF
V600E. PD-L1 TPS 60%, TMB-high.

**Output:** Tier I-A (targeted): dabrafenib + trametinib,
encorafenib + binimetinib, vemurafenib + cobimetinib. Tier I-A
(immunotherapy): pembrolizumab, nivolumab -- supported by PD-L1
>= 50% and TMB-H rules. Treatment sequencing note: agent flags
targeted-first vs. immunotherapy-first debate with supporting
literature from `onco_literature`. Resistance context: MAPK
reactivation (MEK1/2, NRAS mutations) documented; MEK inhibitor
combination partially addresses this. Trial matches: 3 trials
including BRAF + MEK + anti-PD-1 triplet combinations.

---

## 10. Discussion

### 10.1 Strengths

The architecture offers five key advantages: (1) **evidence
completeness** through 11 federated collections spanning the full
decision landscape; (2) **weighted prioritization** ensuring
clinically critical evidence takes precedence; (3) **resistance
awareness** addressing a critical gap in existing tools; (4)
**accessibility** via single DGX Spark deployment at $3,999; and
(5) **interoperability** through FHIR R4 export for EHR
integration.

### 10.2 Limitations

1. **Evidence currency.** The system depends on collection freshness;
   rapidly evolving fields may outpace the update cycle.
2. **LLM synthesis fidelity.** While Claude Sonnet 4.6 demonstrates
   strong biomedical reasoning, outputs may occasionally introduce
   imprecision. All results carry a decision-support disclaimer.
3. **Validation scope.** The four demo cases exercise common
   scenarios; rare tumors, multi-driver cases, and pediatric
   oncology require additional validation.
4. **Prospective validation.** Current validation is retrospective
   and synthetic. Prospective studies comparing agent recommendations
   to MTB consensus are needed to establish clinical utility.
5. **Regulatory status.** The system is a research tool, not
   FDA-cleared or CE-marked, and is not intended for independent
   clinical use without physician oversight.

### 10.3 Future Work

Planned extensions: automated collection refresh via scheduled
CIViC/ClinicalTrials.gov/PubMed ingestion; multi-omics integration
(RNA-seq, CNV, methylation); prospective multi-site validation;
expanded coverage for rare tumors and pediatric malignancies;
treatment sequencing optimization using real-world outcomes; and
pharmacogenomic integration (DPYD, UGT1A1) for toxicity flagging.

---

## 11. Conclusion

The Precision Oncology Intelligence Agent demonstrates that a
multi-collection RAG architecture can provide comprehensive,
resistance-aware, trial-matched molecular tumor board decision
support on commodity hardware. By federating 11 vector collections
under a weighted search planner, applying AMP/ASCO/CAP-aligned
evidence tiering, and exporting results in interoperable formats
including FHIR R4, the system addresses the critical information
synthesis challenge facing molecular tumor boards.

With 556 passing tests, support for 26 cancer types and 40+
actionable gene targets, and deployment on a single NVIDIA DGX
Spark at $3,999, the system is positioned as an accessible,
open-source alternative to commercial precision oncology platforms.
The Apache 2.0 license ensures that community oncology practices,
resource-limited institutions, and academic researchers can deploy
and extend the system without licensing barriers.

Precision oncology is fundamentally an information problem. The
volume, heterogeneity, and velocity of molecular evidence exceed
human cognitive capacity. RAG-powered decision support that combines
structured biomedical knowledge with large language model synthesis
offers a path toward democratizing the molecular tumor board --
ensuring that every patient, regardless of institutional resources,
benefits from comprehensive, evidence-based therapy selection.

---

## 12. References

1. Li, M. M., et al. "Standards and guidelines for the
   interpretation and reporting of sequence variants in cancer."
   *Journal of Molecular Diagnostics* 19.1 (2017): 4-23.

2. Griffith, M., et al. "CIViC is a community knowledgebase for
   expert crowdsourcing the clinical interpretation of variants
   in cancer." *Nature Genetics* 49.2 (2017): 170-174.

3. Chakravarty, D., et al. "OncoKB: a precision oncology knowledge
   base." *JCO Precision Oncology* 1 (2017): 1-16.

4. Xiao, S., et al. "C-Pack: Packaged resources to advance general
   Chinese embedding." *arXiv:2309.07597* (2023).

5. NCCN Clinical Practice Guidelines in Oncology. National
   Comprehensive Cancer Network. https://www.nccn.org/guidelines

6. ClinicalTrials.gov. U.S. National Library of Medicine.
   https://clinicaltrials.gov

7. Benson, A. B., et al. "Colon cancer, version 2.2021." *JNCCN*
   19.3 (2021): 329-359.

8. Planchard, D., et al. "Metastatic non-small cell lung cancer:
   ESMO Clinical Practice Guidelines." *Annals of Oncology* 29
   (2018): iv192-iv237.

9. Skoulidis, F., et al. "Sotorasib for lung cancers with KRAS
   p.G12C mutation." *NEJM* 384.25 (2021): 2371-2381.

10. Eisenhauer, E. A., et al. "New response evaluation criteria in
    solid tumours: revised RECIST guideline (version 1.1)."
    *Eur J Cancer* 45.2 (2009): 228-247.

11. HL7 FHIR Genomics Implementation Guide. Health Level Seven
    International. https://www.hl7.org/fhir/genomics.html

12. Wang, J., et al. "Milvus: A purpose-built vector data management
    system." *ACM SIGMOD* (2021).

13. Long, G. V., et al. "Dabrafenib plus trametinib in BRAF
    V600-mutant metastatic melanoma." *Lancet Oncol* 16.13 (2015):
    1515-1527.

14. Robson, M., et al. "Olaparib for metastatic breast cancer in
    patients with a germline BRCA mutation." *NEJM* 377.6 (2017):
    523-533.

15. Drilon, A., et al. "Efficacy of larotrectinib in TRK
    fusion-positive cancers." *NEJM* 378.8 (2018): 731-739.

---

*This work is part of the HCLS AI Factory, an open-source precision
medicine platform. Source code available under the Apache 2.0
license.*
