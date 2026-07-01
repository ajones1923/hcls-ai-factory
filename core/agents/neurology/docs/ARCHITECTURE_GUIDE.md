# Neurology Intelligence Agent -- Architecture Guide

Part of the HCLS AI Factory Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine) that compose the end-to-end precision medicine platform on NVIDIA DGX Spark.

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## 1. System Architecture

The Neurology Intelligence Agent follows a layered architecture with clear separation between presentation, intelligence, and data layers.

```
+================================================================+
|                     PRESENTATION LAYER                          |
|  +---------------------+   +-----------------------------+     |
|  | Streamlit Chat UI   |   | FastAPI REST API             |     |
|  | Port 8529           |   | Port 8528                    |     |
|  | Interactive Q&A     |   | Versioned endpoints (v1)     |     |
|  | Scale calculators   |   | CORS, auth, rate limiting    |     |
|  +----------+----------+   +-------------+---------------+     |
+================================================================+
               |                            |
+================================================================+
|                    INTELLIGENCE LAYER                           |
|  +----------------+  +----------------+  +------------------+  |
|  | Query Expander |  | Workflow Engine |  | Clinical Scale   |  |
|  | 251+ aliases   |  | 8+1 workflows  |  | Calculators      |  |
|  | 16 synonym     |  | domain-specific|  | 10 validated     |  |
|  | maps           |  | weight boost   |  | instruments      |  |
|  +-------+--------+  +-------+--------+  +--------+---------+  |
|          |                   |                     |            |
|  +-------v-------------------v---------------------v---------+  |
|  |              Neurology RAG Engine                         |  |
|  | ThreadPoolExecutor parallel search across 14 collections  |  |
|  | Workflow-specific collection weight boosting              |  |
|  | Citation scoring (high/medium/low)                        |  |
|  | Conversation memory (24h TTL)                             |  |
|  +---------------------------+-------------------------------+  |
|                              |                                  |
|  +---------------------------v-------------------------------+  |
|  |              Claude LLM (Anthropic)                       |  |
|  | Evidence synthesis with clinical system prompt            |  |
|  | Guideline-grounded recommendations (AAN/AHA/ILAE/MDS)    |  |
|  +-----------------------------------------------------------+  |
+================================================================+
               |
+================================================================+
|                        DATA LAYER                              |
|  +-----------------------------------------------------------+  |
|  |              Milvus 2.4 Vector Database                   |  |
|  |  14 collections | BGE-small 384-dim | IVF_FLAT/COSINE    |  |
|  |  855K estimated records | etcd + MinIO backend            |  |
|  +-----------------------------------------------------------+  |
+================================================================+
```

---

## 2. Clinical Scale Calculator Architecture

Each of the 10 clinical scale calculators follows a consistent design pattern:

```
Input (Dict/scalar)
        |
        v
  +------------------+
  | Input Validation  |  Clamp values to valid ranges
  +--------+---------+
           |
           v
  +--------+---------+
  | Score Computation |  Sum items, apply weights/rules
  +--------+---------+
           |
           v
  +--------+---------+
  | Threshold Eval    |  Compare against clinical cutoffs
  +--------+---------+
           |
           v
  +--------+---------+
  | Interpretation    |  Category, severity, free-text
  +--------+---------+
           |
           v
  +--------+---------+
  | Recommendations   |  Evidence-based action items
  +--------+---------+
           |
           v
  ScaleResult(
      scale_type,
      score,
      max_score,
      interpretation,
      severity_category,
      thresholds,
      recommendations
  )
```

### Scale Calculator Specifications

| Calculator | Input Type | Score Range | Key Decision Points |
|---|---|---|---|
| `NIHSSCalculator` | Dict[str, int] (15 items) | 0-42 | >=1 tPA, >=6 CTA/thrombectomy, >25 ICU |
| `GCSCalculator` | eye, verbal, motor ints | 3-15 | <=8 intubate, <=12 CT required |
| `MoCACalculator` | Dict[str, int] (7 domains) + education | 0-30 | <26 MCI, <18 AD biomarkers, <10 severe |
| `UPDRSCalculator` | Dict[str, int] (33 items) | 0-132 | >=40 optimize meds, >=59 DBS eval |
| `EDSSCalculator` | Dict[str, int] (7 FS) + step | 0-10 | >=4 rehab, >=6 walking aid, >=7 wheelchair |
| `mRSCalculator` | Single int | 0-6 | <=2 good outcome, >=4 poor outcome |
| `HIT6Calculator` | List[str] (6 responses) | 36-78 | >=56 preventive Rx, >=60 CGRP |
| `ALSFRSCalculator` | Dict[str, int] (12 items) | 0-48 | <30 multidisciplinary, >1.0/mo rapid |
| `ASPECTSCalculator` | List[str] (affected regions) | 0-10 | >=6 thrombectomy, <6 large core |
| `HoehnYahrCalculator` | Single float | 1-5 | >=3 DBS eval, >=4 multidisciplinary |

---

## 3. Stroke Triage Pipeline

The acute stroke triage pipeline is the most time-critical workflow. It integrates NIHSS, ASPECTS, and treatment eligibility into a unified decision pathway.

```
Stroke Alert
     |
     v
+----+----+
| NIHSS   |  Calculate stroke severity (0-42)
| Calc    |  Determine tPA consideration threshold
+---------+
     |
     +---> NIHSS >= 6? ---> Order CTA Head/Neck for LVO
     |
     v
+----+----+
| ASPECTS |  Score CT for ischemic changes (10 MCA regions)
| Calc    |  Determine infarct core size
+---------+
     |
     +---> ASPECTS >= 6? ---> Favorable for intervention
     |
     v
+---------+
| Time    |  Determine treatment window
| Window  |
+---------+
     |
     +---> 0-4.5h? -------> tPA eligibility check
     |                       (NINDS/ECASS III criteria)
     +---> 6-16h? --------> DEFUSE-3 criteria
     |                       (core <70mL, mismatch >=1.8)
     +---> 6-24h? --------> DAWN criteria
     |                       (clinical-imaging mismatch)
     v
+---------+
| RAG     |  Search neuro_cerebrovascular (0.25),
| Search  |  neuro_imaging (0.18), neuro_guidelines (0.12)
+---------+
     |
     v
+---------+
| LLM     |  Synthesize: severity, eligibility,
| Synth   |  contraindications, secondary prevention
+---------+
     |
     v
WorkflowResult with:
  - stroke_type (ischemic/hemorrhagic/TIA/SAH)
  - nihss_score + interpretation
  - aspects_score + affected regions
  - tpa_eligible (bool + reasoning)
  - thrombectomy_eligible (bool + criteria met)
  - toast_classification
  - recommendations (prioritized)
  - guideline_references (AHA/ASA)
```

---

## 4. ATN Biomarker Staging

The dementia evaluation workflow implements the NIA-AA ATN (Amyloid-Tau-Neurodegeneration) framework for Alzheimer's disease biological staging.

```
Patient Data
     |
     v
+----+-----------+
| A (Amyloid)    |  CSF Abeta42/40 ratio, amyloid PET
| Binary: +/-    |  Threshold: Abeta42 < 500 pg/mL = A+
+----------------+
     |
     v
+----+-----------+
| T (Tau)        |  CSF p-tau 181 or 217, tau PET
| Binary: +/-    |  Threshold: p-tau181 > 27 pg/mL = T+
+----------------+
     |
     v
+----+-----------+
| N (Neuro-      |  CSF total tau, NfL, FDG-PET,
|   degeneration)|  hippocampal volume on MRI
| Binary: +/-    |  Threshold: t-tau > 400 pg/mL = N+
+----------------+
     |
     v
+----+-----------+
| ATN Stage      |  8 possible combinations:
| Classification |  A-T-N-, A+T-N-, A+T+N-, A+T+N+,
|                |  A+T-N+, A-T+N-, A-T-N+, A-T+N+
+----------------+
     |
     v
Clinical Interpretation:
  A-T-N-  -->  Normal AD biomarkers
  A+T-N-  -->  Alzheimer's pathologic change (preclinical)
  A+T+N-  -->  Alzheimer's disease (early biological)
  A+T+N+  -->  Alzheimer's disease (full biological)
  A+T-N+  -->  Alzheimer's with neurodegeneration (non-tau)
  A-T+N-  -->  Non-AD pathologic change (SNAP or primary tauopathy)
  A-T-N+  -->  Non-AD neurodegeneration
  A-T+N+  -->  Non-AD pathologic change with neurodegeneration
     |
     v
Treatment Eligibility:
  A+T+/-N+/- with MoCA 18-25 --> Anti-amyloid therapy candidate
  (lecanemab/donanemab if amyloid confirmed on PET or CSF)
```

---

## 5. Multi-Collection Search Architecture

```
Query
  |
  v
+-------------------+
| Query Expansion   |  1. Entity alias resolution (251+ mappings)
| Module            |  2. Synonym expansion (16 domain maps)
|                   |  3. Workflow term injection
+--------+----------+
         |
         v
+--------+----------+
| Workflow Router   |  Select workflow type (9 options)
|                   |  Load collection weight profile
+--------+----------+
         |
         v
+--------+----------+
| BGE-small Encoder |  Encode expanded query to 384-dim vector
+--------+----------+
         |
         v
+--------+---------------------------+
| ThreadPoolExecutor                 |
| (max_workers = 14)                 |
|                                    |
|  neuro_literature    (weight=0.08) |
|  neuro_trials        (weight=0.06) |
|  neuro_imaging       (weight=0.09) |
|  neuro_electrophys.  (weight=0.07) |
|  neuro_degenerative  (weight=0.09) |     Boosted weights per workflow:
|  neuro_cerebrovasc.  (weight=0.09) |     e.g., acute_stroke boosts
|  neuro_epilepsy      (weight=0.08) |     neuro_cerebrovascular to 0.25
|  neuro_oncology      (weight=0.06) |
|  neuro_ms            (weight=0.07) |
|  neuro_movement      (weight=0.07) |
|  neuro_headache      (weight=0.06) |
|  neuro_neuromuscular (weight=0.06) |
|  neuro_guidelines    (weight=0.07) |
|  genomic_evidence    (weight=0.05) |
+--------+---------------------------+
         |
         v  (parallel results, each scored)
+--------+----------+
| Score Aggregation  |  weighted_score = cosine_sim * weight
| & Ranking          |  Filter by SCORE_THRESHOLD (0.4)
|                    |  Sort by weighted_score descending
+--------+----------+
         |
         v
+--------+----------+
| Citation Scoring   |  >= 0.75 = high relevance
|                    |  >= 0.60 = medium relevance
|                    |  < 0.60 = low relevance
+--------+----------+
         |
         v
Top-K results with metadata, citations, evidence levels
```

---

## 6. Workflow Engine Design

```python
class BaseNeuroWorkflow(ABC):
    """Abstract base for all neurology clinical workflows."""

    workflow_type: NeuroWorkflowType
    domain: Optional[NeuroDomain] = None

    def run(self, inputs: dict) -> WorkflowResult:
        """Template method: preprocess -> execute -> postprocess."""
        preprocessed = self.preprocess(inputs)
        result = self.execute(preprocessed)
        return self.postprocess(result)

    @abstractmethod
    def preprocess(self, inputs: dict) -> dict: ...

    @abstractmethod
    def execute(self, data: dict) -> WorkflowResult: ...

    @abstractmethod
    def postprocess(self, result: WorkflowResult) -> WorkflowResult: ...
```

Each workflow implements domain-specific logic:
- **Preprocess:** Validate inputs, compute relevant scales, extract structured data
- **Execute:** Run RAG search with workflow-specific weights, call LLM for synthesis
- **Postprocess:** Add cross-modal triggers, format recommendations, attach guidelines

---

## 7. Data Flow Diagrams

### Query Processing Pipeline

```
User Query (REST/Streamlit)
  |
  +--> Input validation (Pydantic NeuroQuery model)
  |
  +--> Workflow inference (if not specified)
  |      |
  |      +--> Score each workflow's terms against query
  |      +--> Select highest-scoring workflow
  |
  +--> Query expansion (aliases + synonyms + workflow terms)
  |
  +--> RAG search (parallel across 14 collections)
  |
  +--> Evidence ranking (weighted cosine similarity)
  |
  +--> Clinical scale computation (if applicable)
  |
  +--> LLM synthesis (Claude with neuro system prompt)
  |
  +--> Response assembly (NeuroResponse with citations)
  |
  +--> Return to user
```

### Data Ingestion Pipeline

```
Source Data (PubMed, ClinicalTrials.gov, Guidelines)
  |
  +--> Domain-specific parser (pubmed/imaging/eeg)
  |
  +--> Text extraction and normalization
  |
  +--> BGE-small embedding (384-dim)
  |
  +--> Milvus upsert to target collection
  |
  +--> Index rebuild (if needed)
```

---

## 8. Security Architecture

```
External Request
  |
  +--> Rate Limiter (100 req/60s per IP)
  |
  +--> Request Size Check (<= 10 MB)
  |
  +--> API Key Validation (X-API-Key header)
  |      Skip: /health, /healthz, /metrics
  |
  +--> CORS Validation (allowlist from settings)
  |
  +--> Pydantic Input Validation
  |
  +--> Business Logic
  |
  +--> Response (with error handling)
```

---

## 9. Conversation Memory Architecture

```
User Message (with session_id)
  |
  +--> Load conversation from disk
  |      File: data/cache/conversations/{session_id}.json
  |      TTL: 24 hours
  |
  +--> Append to history (max 3 turns)
  |
  +--> Include history in LLM context
  |
  +--> Save updated conversation to disk
  |
  +--> Return response with session_id
```

---

## 10. Cross-Agent Communication

### 10.1 Platform Context: 11 Intelligence Agents

The Neurology Intelligence Agent is one of 11 specialized intelligence agents within the HCLS AI Factory Precision Intelligence Network:

| # | Agent | UI Port | API Port | Domain |
|---|-------|---------|----------|--------|
| 1 | Precision Biomarker | 8528 | 8529 | Biomarker interpretation and trends |
| 2 | Precision Oncology | 8525 | 8526 | Cancer genomics and treatment |
| 3 | CAR-T Intelligence | 8521 | 8522 | CAR-T cell therapy |
| 4 | Medical Imaging | 8523 | 8524 | Multi-modal imaging analysis |
| 5 | Precision Autoimmune | 8531 | 8532 | Autoimmune disease analysis |
| 6 | Cardiology Intelligence | 8536 | 8126 | Cardiovascular decision support |
| 7 | Clinical Trial Intelligence | 8538 | 8537 | Trial matching and eligibility |
| 8 | Neurology Intelligence | 8529 | 8528 | Neurological decision support |
| 9 | Rare Disease Intelligence | 8135 | 8134 | Rare disease diagnosis |
| 10 | Pharmacogenomics (PGx) | 8541 | 8540 | Drug-gene interactions |
| 11 | Pediatric Oncology | 8543 | 8542 | Pediatric cancer intelligence |

### 10.2 Cross-Agent Calls

The Neurology Intelligence Agent calls 5 sibling agents for integrated assessments:

```
Neurology Agent
  |
  +--> Imaging Agent (8524)
  |      Neuroimaging correlation: MRI lesion characterization, DWI/FLAIR mismatch
  |      Requests: brain MRI findings, perfusion data, tumor segmentation
  |
  +--> Cardiology Agent (8126)
  |      Cardioembolic stroke evaluation: AF detection, PFO assessment, echo findings
  |      Requests: cardiac rhythm data, echocardiography results
  |
  +--> Biomarker Agent (8529)
  |      CSF biomarker trending: amyloid, tau, NfL trajectories
  |      Requests: longitudinal biomarker data, threshold alerts
  |
  +--> Clinical Trial Agent (8537)
  |      Matches neurology patients to active trials (anti-amyloid, gene therapy, DMTs)
  |      Requests: trial eligibility assessment
  |
  +--> Rare Disease Agent (8134)
         Rare neurogenetic disorder evaluation: inherited ataxias, leukodystrophies
         Requests: rare disease differential, genetic variant interpretation
```

### 10.3 Integrated Assessment Endpoint

```
POST /v1/neuro/integrated-assessment

Request:
  patient_id: str
  include_agents: List[str]  # ["imaging", "cardiology", "biomarker", "trial", "rare_disease"]
  clinical_context: dict

Response:
  neurology_assessment: NeuroResponse
  cross_agent_findings: List[CrossAgentResult]
  pediatric_neurotoxicity: Optional[PediatricNeurotoxicityAssessment]
  integrated_recommendations: List[str]
```

### 10.4 Legacy Cross-Agent Flow

```
Neurology Agent
  |
  +--> Cross-modal trigger detected
  |      (e.g., stroke patient with genetic variant)
  |
  +--> HTTP request to target agent
  |      Timeout: 30 seconds
  |
  +--> Merge external evidence into response
  |
  +--> SSE event broadcast to subscribers
```

---

## 11. Deployment Topology

### Standalone Mode

```
+----- Docker Host -----+
|                        |
|  neuro-network         |
|  +-----------+         |
|  | etcd      |         |
|  | MinIO     |         |
|  | Milvus    |-59530   |
|  +-----------+         |
|  | neuro-api |-8528    |
|  | neuro-ui  |-8529    |
|  | neuro-    |         |
|  |   setup   |         |
|  +-----------+         |
+------------------------+
```

### Integrated Mode (DGX Spark)

```
+----- DGX Spark --------+
|                         |
|  Shared Milvus (19530)  |
|  Shared etcd + MinIO    |
|                         |
|  neuro-api    (8528)    |
|  neuro-ui     (8529)    |
|  + 8 other agents       |
|  + landing page (8080)  |
|  + Prometheus + Grafana  |
+-------------------------+
```

---

*Neurology Intelligence Agent -- Architecture Guide v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
