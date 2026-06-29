# Clinical Trial Intelligence Agent -- Demo Guide

**Version:** 2.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Pre-Demo Checklist](#1-pre-demo-checklist)
2. [Starting the System](#2-starting-the-system)
3. [UI Walkthrough: 5 Tabs](#3-ui-walkthrough-5-tabs)
4. [Demo Scenario 1: Oncology Protocol Design](#4-demo-scenario-1-oncology-protocol-design)
5. [Demo Scenario 2: Patient-Trial Matching](#5-demo-scenario-2-patient-trial-matching)
6. [Demo Scenario 3: Eligibility Optimization](#6-demo-scenario-3-eligibility-optimization)
7. [Demo Scenario 4: Competitive Intelligence](#7-demo-scenario-4-competitive-intelligence)
8. [Demo Scenario 5: Adaptive Design Selection](#8-demo-scenario-5-adaptive-design-selection)
9. [API Demo Queries](#9-api-demo-queries)
10. [Talking Points](#10-talking-points)
11. [FAQ](#11-faq)
12. [Recovery Procedures](#12-recovery-procedures)

---

## 1. Pre-Demo Checklist

Run through this checklist before every demo:

| # | Check | Command | Expected |
|---|---|---|---|
| 1 | Tests pass | `python -m pytest tests/ -q` | 769 passed in 0.47s |
| 2 | API starts | `curl localhost:8538/health` | `{"status": "ok"}` |
| 3 | UI loads | Open `http://localhost:8128` | NVIDIA-themed 5-tab UI |
| 4 | Collections exist | `curl localhost:8538/collections` | 14 collections listed |
| 5 | Knowledge loaded | `curl localhost:8538/v1/trial/knowledge-version` | version 2.0.0 |
| 6 | Workflows available | `curl localhost:8538/workflows` | 10+ workflow types |
| 7 | Milvus running | Check port 19530 | Connected status in health |
| 8 | Browser configured | Close unrelated tabs | Clean browser window |
| 9 | Terminal ready | Two terminals open | API + backup commands |
| 10 | API key set | Check .env | ANTHROPIC_API_KEY present |

---

## 2. Starting the System

### Quick Start (Integrated)

```bash
cd /home/adam/projects/hcls-ai-factory
./start-factory.sh
```

### Quick Start (Standalone)

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/clinical_trial_intelligence_agent

# Terminal 1: API
uvicorn api.main:app --host 0.0.0.0 --port 8538

# Terminal 2: UI
streamlit run app/trial_ui.py --server.port 8128
```

### Verify Both Services

```bash
curl -s localhost:8538/health | python -m json.tool
# Open http://localhost:8128 in browser
```

---

## 3. UI Walkthrough: 5 Tabs

### Tab 1: Trial Intelligence

**Purpose:** RAG-powered Q&A across all 14 collections with workflow-aware routing.

**What to show:**
- Type a clinical trial question in the text input
- System detects the appropriate workflow automatically
- Results show: answer, citations with collection source and relevance scores, guideline references, confidence score
- Point out the NVIDIA dark theme and responsive layout

**Sample queries to type:**
- "What are the standard Phase 3 endpoints for NSCLC immunotherapy trials?"
- "Compare adaptive design options for a rare disease Phase 2/3 study"
- "What regulatory pathways are available for breakthrough therapy designation?"

### Tab 2: Patient Matching

**Purpose:** Match a patient profile against clinical trial eligibility criteria.

**What to show:**
- Enter patient demographics (age, sex, diagnosis)
- Add biomarkers (e.g., PD-L1 TPS 80%, EGFR L858R)
- Add genomic variants and current medications
- Submit to see per-criterion match scores
- Highlight the overall match score and confidence level
- Show nearby trial sites if geographic location is provided

### Tab 3: Protocol Optimizer

**Purpose:** Protocol complexity scoring and design recommendations.

**What to show:**
- Enter trial parameters: indication, phase, procedure count, visit count, endpoint count, eligibility criteria count
- Submit to see complexity score and percentile rank
- Point out comparisons to Tufts CSDD industry benchmarks
- Show endpoint recommendations by indication
- Demonstrate historical success rate lookup

### Tab 4: Competitive Landscape

**Purpose:** Competitive threat assessment and landscape visualization.

**What to show:**
- Enter an indication and mechanism of action
- See competitor trial listing with threat scores
- Explain the 4-factor threat model (phase, enrollment, sponsor, differentiation)
- Highlight threat classification (critical/high/moderate/low/minimal)
- Show enrollment progress comparison

### Tab 5: Dashboard

**Purpose:** System health and operational metrics.

**What to show:**
- Collection health indicators (14 collections with record counts)
- Query volume metrics
- Workflow execution breakdown
- Cross-agent integration status
- Knowledge base version and last update

---

## 4. Demo Scenario 1: Oncology Protocol Design

**Story:** "A pharma company wants to design a Phase 3 trial for a novel PD-L1/VEGF bispecific antibody in first-line NSCLC. Let's use the agent to generate an evidence-based protocol blueprint."

### Step 1: Ask the Question

Go to **Tab 1 (Trial Intelligence)** and type:

> "Design a Phase 3 protocol for a PD-L1/VEGF bispecific antibody in first-line metastatic NSCLC with PD-L1 TPS >= 50%"

### Expected Output

The system should return:

- **Workflow detected:** Protocol Design
- **Recommended primary endpoint:** Progression-Free Survival (PFS) or Overall Survival (OS)
- **Recommended secondary endpoints:** ORR, DoR, OS (if PFS is primary), quality of life
- **Comparator recommendation:** Pembrolizumab monotherapy (based on KEYNOTE-024)
- **Sample size estimate:** 300-600 patients based on expected HR and event rate
- **Eligibility framework:** PD-L1 TPS >= 50%, ECOG 0-1, no prior systemic therapy
- **Adaptive design suggestion:** Group sequential with interim analysis for OS
- **Regulatory pathway:** Potential for Breakthrough Therapy Designation; RTOR eligible
- **Landmark trial references:** KEYNOTE-024 (pembrolizumab), IMpower110 (atezolizumab)
- **Historical success rate:** ~36% for Phase 3 oncology

### Talking Points

- "Notice how the system referenced KEYNOTE-024 as the benchmark trial -- that's the current standard of care for this indication"
- "The PFS and OS dual primary endpoint strategy is exactly what FDA expects for first-line NSCLC"
- "The 36% historical success rate is calibrated from BIO/QLS Advisors data across thousands of oncology trials"

---

## 5. Demo Scenario 2: Patient-Trial Matching

**Story:** "A 58-year-old male with metastatic NSCLC, PD-L1 TPS 80%, EGFR wild-type, no prior systemic therapy, and ECOG 1 wants to find matching clinical trials."

### Step 1: Enter Patient Profile

Go to **Tab 2 (Patient Matching)** and enter:

| Field | Value |
|---|---|
| Age | 58 |
| Sex | Male |
| Diagnosis | Metastatic non-small cell lung cancer |
| Biomarkers | PD-L1 TPS 80%, EGFR wild-type, ALK negative |
| Medications | None (treatment-naive) |
| Genomic variants | KRAS G12C |
| Comorbidities | Hypertension (controlled) |
| Location | Boston, MA |

### Expected Output

- Multiple trial matches with per-criterion scoring
- High match score for first-line IO trials (PD-L1 >= 50%)
- KRAS G12C targeted therapy trials flagged
- Controlled hypertension should pass cardiac exclusion criteria
- Nearby Boston-area trial sites listed

### Talking Points

- "Each eligibility criterion gets its own confidence score -- the system doesn't just say yes or no"
- "Notice the KRAS G12C variant triggered an additional match for targeted therapy trials"
- "The cross-agent trigger sent this to the oncology agent for molecular match confirmation"

---

## 6. Demo Scenario 3: Eligibility Optimization

**Story:** "A trial protocol has 28 eligibility criteria and is struggling with enrollment. The medical team wants to know which criteria could be broadened without compromising safety."

### Step 1: Submit Criteria for Analysis

Go to **Tab 1 (Trial Intelligence)** and type:

> "Analyze these eligibility criteria for enrollment impact: ECOG 0-1, no prior immunotherapy, no CNS metastases, hemoglobin >= 10 g/dL, no autoimmune disease, creatinine clearance >= 30 mL/min, no prior organ transplant, ejection fraction >= 50%, pregnancy or lactation excluded"

### Expected Output

The system should return criteria ranked by population impact:

| Criterion | Population Impact | Recommendation |
|---|---|---|
| Pregnancy or lactation | 50% | RETAIN with monitoring |
| ECOG 0-1 | 25% | REVIEW: consider ECOG 0-2 |
| No prior immunotherapy | 20% | BROADEN: weak justification for some lines |
| No CNS metastases | 15% | REVIEW: stable CNS may be includable |
| Ejection fraction >= 50% | 8% | RETAIN: cardiac safety justified |
| Creatinine clearance >= 30 | 15% | REVIEW: could lower threshold |
| No autoimmune disease | 8% | REVIEW: mild autoimmune may be safe |

### Talking Points

- "The system identified that 'no prior immunotherapy' excludes 20% of the population but has weak scientific justification for this particular trial"
- "Each criterion's population impact is calibrated against real-world data -- not guesswork"
- "RETAIN recommendations have strong scientific backing; BROADEN recommendations have high impact with weak justification"

---

## 7. Demo Scenario 4: Competitive Intelligence

**Story:** "The sponsor wants to understand the competitive landscape for GLP-1 receptor agonists in obesity before starting their Phase 3 program."

### Step 1: Query the Landscape

Go to **Tab 4 (Competitive Landscape)** or type in Tab 1:

> "What is the competitive landscape for GLP-1 agonists in obesity Phase 3 trials? Include semaglutide, tirzepatide, and survodutide."

### Expected Output

- Competitor profiles for SURMOUNT-1 (tirzepatide), SELECT (semaglutide), and others
- Threat scores with classification:
  - Tirzepatide/Zepbound: **Critical** (0.85+) -- Phase 3 completed, 22.5% weight loss
  - Semaglutide/Wegovy: **Critical** (0.90+) -- Approved, MACE benefit proven
  - Survodutide: **High** (0.65) -- Phase 3 ongoing
- Enrollment progress comparison
- Differentiation analysis (mechanism similarity)

### Talking Points

- "The threat model considers four factors: phase, enrollment, sponsor resources, and mechanism similarity"
- "Semaglutide scores as critical because it's already approved AND has cardiovascular outcome data"
- "This analysis helps sponsors decide whether to differentiate on efficacy, safety, convenience, or a new endpoint"

---

## 8. Demo Scenario 5: Adaptive Design Selection

**Story:** "A biotech company is developing a novel gene therapy for a rare disease (hemophilia A) and wants to know which adaptive design is best for their Phase 2/3 program."

### Step 1: Ask About Adaptive Designs

Type in **Tab 1**:

> "What adaptive trial design should I use for a Phase 2/3 seamless gene therapy trial in hemophilia A? Consider the small patient population and uncertain dose-response."

### Expected Output

- **Recommended design:** Seamless Phase 2/3 with sample size re-estimation
- **Alternative designs:** Response-adaptive randomization, biomarker-adaptive enrichment
- **Rationale:** Small population (rare disease), uncertain effect size, long follow-up needed
- **Regulatory guidance:** FDA Guidance on Adaptive Designs (2019), RMAT designation eligible
- **Precedent trials:** SPRINT_SMA (gene therapy adaptive design), Zolgensma experience
- **Statistical considerations:** Bayesian framework, external historical controls, LTFU requirements
- **Success rate estimate:** ~55% for Phase 3 rare disease (higher than average)

### Talking Points

- "The system recommended seamless Phase 2/3 because it eliminates the inter-trial gap -- critical for rare diseases where every patient matters"
- "FDA RMAT designation makes this eligible for accelerated timelines"
- "The 55% Phase 3 success rate for rare diseases is significantly higher than the overall average of ~50%, reflecting smaller trials and often dramatic clinical effects"

---

## 9. API Demo Queries

For technical audiences, demonstrate the REST API directly:

### 9.1 RAG Query

```bash
curl -X POST http://localhost:8538/v1/trial/query \
    -H "Content-Type: application/json" \
    -d '{
        "question": "What biomarker strategies are used in NSCLC immunotherapy trials?",
        "workflow_type": "protocol_design",
        "top_k": 5
    }'
```

### 9.2 Patient Matching

```bash
curl -X POST http://localhost:8538/v1/trial/match \
    -H "Content-Type: application/json" \
    -d '{
        "patient": {
            "age": 65,
            "sex": "female",
            "diagnosis": "HER2-low metastatic breast cancer",
            "biomarkers": ["HER2 IHC 1+", "HR-positive", "BRCA wild-type"]
        }
    }'
```

### 9.3 Safety Signal Detection

```bash
curl -X POST http://localhost:8538/v1/trial/safety/signal \
    -H "Content-Type: application/json" \
    -d '{
        "events": [
            {"event_type": "hepatotoxicity", "severity": "grade_3", "frequency": 0.08},
            {"event_type": "neutropenia", "severity": "grade_4", "frequency": 0.12}
        ],
        "trial_id": "NCT00000001"
    }'
```

### 9.4 Knowledge Version

```bash
curl -s http://localhost:8538/v1/trial/knowledge-version | python -m json.tool
```

### 9.5 Therapeutic Areas

```bash
curl -s http://localhost:8538/v1/trial/therapeutic-areas | python -m json.tool
```

---

## 10. Talking Points

### For Executive Audiences

- "This agent reduces the time to generate an evidence-based protocol blueprint from weeks to seconds"
- "40 landmark trials from KEYNOTE-024 to RECOVERY to CLARITY-AD serve as real-world design precedents"
- "The system covers 13 therapeutic areas, from oncology to gene therapy, with calibrated success rates"
- "Running on DGX Spark, this integrates with genomics, drug discovery, and four peer AI agents"

### For Clinical Operations Audiences

- "Patient-trial matching evaluates every eligibility criterion with confidence scoring"
- "Eligibility optimization identifies criteria that exclude patients without scientific justification"
- "Site selection considers enrollment history, diversity index, and screen failure rates"
- "Enrollment predictions use prevalence, competition, and phase-specific factors"

### For Regulatory Affairs Audiences

- "The system references ICH E6(R3), E9(R1), and E8(R1) guidelines"
- "9 regulatory agencies covered with approval pathways and expedited programs"
- "Adaptive design recommendations include specific FDA guidance citations"
- "Every recommendation includes evidence level classification (A1 through E)"

### For Data Science Audiences

- "14 Milvus vector collections with IVF_FLAT indexing and COSINE similarity"
- "384-dimensional BGE-small-en-v1.5 embeddings for clinical trial text"
- "Workflow-specific collection weight boosting for domain-relevant retrieval"
- "Calibrated confidence model: evidence base (0.3) + raw confidence (0.3) + documents (0.2) + agreement (0.2)"
- "769 tests in 0.47 seconds -- 100% pass rate"

---

## 11. FAQ

**Q: Does this replace clinical trial design teams?**
A: No. This is a decision support tool that accelerates evidence retrieval and provides calibrated recommendations. All outputs should be reviewed by qualified clinical, regulatory, and statistical professionals.

**Q: How current is the trial data?**
A: The knowledge base contains 40 curated landmark trials and comprehensive reference data. Real-time ClinicalTrials.gov integration is available through the ingest pipeline on a configurable schedule (default 24 hours).

**Q: What happens if Milvus is down?**
A: The system degrades gracefully. All 10 workflows, 5 decision support engines, and the knowledge base continue to function. Only vector search results are unavailable.

**Q: What happens without an LLM API key?**
A: The system operates in search-only mode. Vector searches, workflow execution, and decision support engines all function normally. Natural language synthesis is unavailable.

**Q: Can this be used for regulatory submissions?**
A: The system generates draft documents and evidence summaries. All regulatory submissions require expert review, validation, and formatting per agency-specific requirements.

**Q: What therapeutic areas are supported?**
A: All 13 major therapeutic areas: oncology, cardiovascular, neuroscience, immunology, infectious disease, rare diseases, metabolic/endocrinology, respiratory, hematology, gastroenterology, dermatology, ophthalmology, and gene/cell therapy.

---

## 12. Recovery Procedures

### If the API Crashes

```bash
# Restart API
uvicorn api.main:app --host 0.0.0.0 --port 8538
```

### If the UI Freezes

```bash
# Restart Streamlit (it auto-reconnects)
streamlit run app/trial_ui.py --server.port 8128
```

### If Milvus Is Unreachable

```bash
# Check Milvus status
docker ps | grep milvus

# Restart if needed
docker restart milvus-standalone

# Verify
curl localhost:8538/health
```

### If Searches Return Empty

```bash
# Re-seed knowledge base
python scripts/seed_knowledge.py

# Re-create collections if needed
python scripts/setup_collections.py
python scripts/seed_knowledge.py
```

---

*Clinical Trial Intelligence Agent v2.0.0 -- Demo Guide -- March 2026*
