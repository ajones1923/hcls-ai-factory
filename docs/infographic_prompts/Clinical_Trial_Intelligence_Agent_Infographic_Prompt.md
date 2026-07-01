# Nano Banana Pro — Clinical Trial Intelligence Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Clinical Trial Intelligence Agent — the HCLS AI Factory's most operationally comprehensive agent, covering the complete clinical trial lifecycle: protocol design optimization, automated patient-trial matching, site selection scoring, enrollment prediction, eligibility optimization, adaptive design evaluation, safety signal detection, regulatory document generation, competitive landscape analysis, biomarker strategy, RWE analysis, diversity assessment, and decentralized trial planning. It features 13 therapeutic areas, 40 landmark trials, 9 regulatory agencies, 9 endpoint types, 9 adaptive designs, 9 biomarker strategies, 9 DCT components, 6 decision support engines, 19 workflow types, 14 Milvus collections, and a 5-tab Streamlit UI. It connects to 8 peer agents (the most of any non-orchestrator agent) and appears in 4 of 6 demos — more than any other non-Oncology agent. Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as all other HCLS AI Factory agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a clinical operations command center poster — comprehensive enough for a VP of Clinical Development, detailed enough for a biostatistician designing an adaptive trial, relevant enough for a pediatric oncologist matching patients to COG protocols.

**Typography:**
- Title: Large, bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Section headers: Bold, navy, green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body: Dark gray (#333333), 8-10pt, SHORT PHRASES only
- Metric badges: White text on colored pills
- All text MUST be legible at 1920x1080

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Red: #DC2626 — safety signals, failed criteria, high risk
- Amber: #F5A623 — moderate match scores, conditional eligibility
- Emerald: #059669 — high match scores, met criteria, enrolled
- Blue: #2196F3 — trial/protocol collections, regulatory
- Purple: #7B2D8E — adaptive designs, biomarker strategy
- Orange: #FF9800 — competitive intelligence, RWE
- Pink: #EC4899 — patient matching, diversity
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — clipboard, patient, globe, safety shield, graph, regulatory stamp motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for safety signals
- Color-coded match scores (emerald=high, amber=moderate, red=low)
- Trial phase badges (I, II, III, IV color-coded)
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Clinical Trial Intelligence Agent" (green #76B900, white text)
- "19 Workflows | 8 Peer Agents" (amber #F5A623, white text)

**Center:**
- Line 1 (large, bold, navy): **Clinical Trial Intelligence Agent**
- Line 2 (gray): "Complete Clinical Trial Lifecycle Intelligence on NVIDIA DGX Spark"
- Line 3 (teal): "Part of the Precision Intelligence Engine — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB | 13 Therapeutic Areas | 40 Landmark Trials | 19 Workflows | 6 Decision Engines | 14 Collections"

**Right side — Capability badge:**
```
TRIAL LIFECYCLE COVERAGE
━━━━━━━━━━━━━━━━━━━━━━━━━
Design → Match → Select → Enroll
Monitor → Analyze → Report → Submit
19 Workflow Types
6 Decision Support Engines
13 Therapeutic Areas
Pediatric COG Protocols
In 4 of 6 Demos
8 Cross-Agent Connections
```

---

### ━━━ LEFT COLUMN: TRIAL KNOWLEDGE BASE ━━━

**Stacked knowledge domain cards:**

1. **13 Therapeutic Areas** [globe icon]
   Oncology, cardiology, neurology, immunology,
   rare disease, infectious disease, metabolic,
   respiratory, dermatology, ophthalmology,
   hematology, psychiatry, musculoskeletal
   [13 areas] (blue badge)

2. **40 Landmark Trials** [star icon]
   ELIANA (CAR-T ALL), KEYNOTE-024 (pembrolizumab),
   DAPA-HF (dapagliflozin), PARADIGM-HF (sacubitril),
   CheckMate-067 (nivo+ipi), IMpower150 (atezo+bev),
   DESTINY-Breast03, SOLO1, AALL0932 (COG)...
   [40 trials] (purple badge)

3. **9 Regulatory Agencies** [stamp icon]
   FDA, EMA, PMDA, NMPA, TGA, Health Canada,
   MHRA, Swissmedic, ANVISA
   [9 agencies] (blue badge)

4. **9 Endpoint Types** [target icon]
   OS, PFS, ORR, DOR, EFS, MRD, PRO, biomarker, composite
   [9 types] (teal badge)

5. **9 Adaptive Designs** [branch icon]
   Bayesian, basket, umbrella, platform,
   dose-finding, group-sequential, response-adaptive,
   biomarker-adaptive, seamless phase II/III
   [9 designs] (purple badge)

6. **9 Biomarker Strategies** [DNA+target icon]
   Enrichment, stratification, prognostic,
   predictive, pharmacodynamic, surrogate,
   companion diagnostic, complementary, monitoring
   [9 strategies] (emerald badge)

7. **9 DCT Components** [home icon]
   E-consent, telemedicine, home nursing,
   wearables, direct-to-patient shipping,
   local labs, digital endpoints, remote monitoring,
   hybrid visit schedules
   [9 DCT] (teal badge)

8. **Safety Signal Metrics** [shield icon]
   Disproportionality analysis (PRR, ROR, BCPNN),
   Bayesian signal detection, time-to-onset,
   dose-response, organ class analysis
   [6 metrics] (red badge)

Arrows rightward: "14 Collections | BGE-small-en-v1.5 | 384-dim"

---

### ━━━ CENTER-TOP: MILVUS COLLECTIONS — 14 ━━━

```
Milvus 2.4 | 14 Collections | IVF_FLAT / COSINE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**Collection badges (2 rows):**

Row 1:
[trial_protocols] [trial_eligibility] [trial_endpoints] [trial_sites] [trial_investigators] [trial_results] [trial_regulatory]
  (blue)            (pink)               (teal)            (emerald)     (blue)                 (emerald)      (blue)

Row 2:
[trial_literature] [trial_biomarkers] [trial_safety] [trial_rwe] [trial_adaptive] [trial_guidelines] [genomic_evidence]
  (blue)             (purple)            (red)          (orange)    (purple)          (blue)             (emerald)

---

### ━━━ CENTER: 6 DECISION SUPPORT ENGINES ━━━

**Section header (navy, green underline):** "6 Decision Support Engines — Quantitative Trial Intelligence"

**Visual: 6 engine cards in a 3x2 grid:**

**Row 1:**

```
[calibrate]                  [complexity]                 [enrollment]
Confidence Calibrator        Protocol Complexity           Enrollment Predictor
━━━━━━━━━━━━━━━━━━━          ━━━━━━━━━━━━━━━━━━━━          ━━━━━━━━━━━━━━━━━━━
Match score calibration      Procedure count, visit        Recruitment rate
Evidence-level weighting     frequency, endpoint count,    estimation based on
Bayesian confidence          criteria complexity, lab       site capacity, disease
interval estimation          burden → composite score       prevalence, competition

Multi-source evidence        [Low/Med/High score]          [months to target N]
aggregation
```

**Row 2:**

```
[eligibility]                [competitive]                [historical]
Eligibility Analyzer         Competitive Threat            Historical Success
━━━━━━━━━━━━━━━━━━           ━━━━━━━━━━━━━━━━━━           ━━━━━━━━━━━━━━━━━━
Criterion-by-criterion       Active trials in same         Phase success
patient screening            indication + mechanism        probability based on
Met/Not Met/Partial          Enrollment competition        therapeutic area,
for each criterion           Differentiation scoring       mechanism, biomarker,
                                                           endpoint, population
[traffic light per           [threat: low/med/high]        [P(success) estimate]
 criterion]
```

---

### ━━━ CENTER: TRIAL LIFECYCLE VISUAL (the unique differentiator) ━━━

**Horizontal lifecycle flow showing all 19 workflow capabilities mapped to trial phases:**

```
CLINICAL TRIAL LIFECYCLE — 19 Workflows Across All Phases
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DESIGN                    CONDUCT                    ANALYZE                  SUBMIT
Phase 0/I                 Phase I-III                Phase II-IV              Regulatory
━━━━━━━━                  ━━━━━━━━                   ━━━━━━━━                 ━━━━━━━━
Protocol Design           Patient Matching           Safety Signal            Regulatory Docs
Eligibility Optimization  Site Selection             RWE Analysis             Regulatory Strategy
Endpoint Strategy         Recruitment Optimization   Competitive Intel
Adaptive Design           Safety Monitoring          Biomarker Strategy
Biomarker Strategy        Diversity Assessment
DCT Planning              Eligibility Analysis

[19 workflows total — covering every stage from concept to submission]
```

**Phase success rate callout (data-driven):**
```
Historical Phase Success Rates (industry average):
Phase I → II: 52%    Phase II → III: 29%    Phase III → Filing: 58%    Filing → Approval: 90%
Overall: ~8% success rate from Phase I to Approval

This agent improves every transition through data-driven optimization.
```

---

### ━━━ CENTER-BOTTOM: PEDIATRIC TRIALS + CROSS-AGENT ━━━

**Pediatric Trial Intelligence box (emerald border):**
```
PEDIATRIC CLINICAL TRIAL INTELLIGENCE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ClinicalTrials.gov Age Filters:
• Child: birth-17 | Adolescent: 13-17 | Adult: 18-64

Children's Oncology Group (COG) Protocols:
• AALL0932: Standard-risk B-ALL
• AALL1231: High-risk B-ALL with MRD-directed therapy
• ANBL1232: High-risk neuroblastoma with immunotherapy
• ARET0321: Retinoblastoma

Precision Programs:
• Pediatric MATCH (NCI-COG): Molecularly targeted for relapsed solid tumors
• INFORM: International pediatric precision oncology registry

Pediatric Trial Design:
• Rolling-6 design (not traditional 3+3 for Phase I)
• EFS preferred endpoint over OS
• Parental consent + child assent (age ≥7)
• Long-term follow-up required (COG LTFU)
```

**Cross-Agent Integration box (teal border — 8 agents, most of any non-orchestrator):**
```
Cross-Agent (/v1/trial/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)       Molecular matches, tumor profiling context
→ PGx Agent (:8107)            Drug-gene interaction screening for trial drugs
→ Cardiology Agent (:8126)     Cardiac eligibility criteria, monitoring requirements
→ Biomarker Agent (:8529)      Biomarker-driven eligibility enrichment
→ Rare Disease Agent (:8134)   Hereditary predisposition trial matching
→ Neurology Agent (:8528)      Neuro-oncology trial requirements
→ Single-Cell Agent (:8540)    Biomarker validation for targeted trials
→ Imaging Agent (:8524)        Trial-specific imaging requirements

8 peer agents — most cross-agent connections of any non-orchestrator agent
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

1. **Patient-Trial Matches** [match icon]
   Ranked trial list with confidence scores
   Criterion-by-criterion eligibility (green/amber/red)
   NCT ID + phase + status + distance

2. **Protocol Optimization** [blueprint icon]
   Design recommendations
   Complexity score reduction
   Adaptive design suggestions

3. **Site Recommendations** [map icon]
   Ranked sites by capability + enrollment
   Geographic optimization
   Diversity assessment

4. **Safety Signal Report** [alert icon]
   Disproportionality analysis
   Signal detection alerts
   Organ class mapping

5. **Competitive Landscape** [radar icon]
   Active competing trials
   Differentiation scoring
   Enrollment threat assessment

6. **Integrated Assessment** [network icon]
   /v1/trial/integrated-assessment
   8-agent coordinated output

7. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4

**Streamlit UI (:8128) — 5 Tabs:**
```
1. Dashboard (health, metrics, therapeutic areas)
2. Trial Explorer (RAG Q&A, collection search)
3. Patient Matcher
   Patient profile input
   Molecular data, demographics, prior therapy
   → Ranked trial matches with eligibility detail
4. Protocol Analyzer
   Protocol upload/input
   Complexity scoring
   Optimization recommendations
5. Competitive Dashboard
   Therapeutic area landscape
   Enrollment competition
   Differentiation analysis
```

**FastAPI (:8538):**
```
/protocol/optimize | /match | /match/batch
/site/recommend | /eligibility/optimize
/adaptive/evaluate | /safety/signal
/regulatory/generate | /competitive/landscape
/diversity/assess | /dct/plan
/v1/trial/integrated-assessment
```

---

### ━━━ BOTTOM STRIP: DEMO PARTICIPATION (4 of 6 demos) ━━━

**Header (navy, white text):** "Clinical Trial Intelligence Agent — In 4 of 6 Demos"

```
Demo 2 ✓ "The 30-Second Tumor Board" — Evelyn, MRD+ ALL
  → COG AALL1731 matched (0.97 score)
  → Pediatric MATCH arm identified
  → AMELIA dual CD19/CD22 CAR-T trial

Demo 4 ✓ "One Gene, One Family" — Aurora, 4yo, Retinoblastoma
  → COG ARET0321 matched
  → NCT05765045, NCT04587544 identified
  → Age-appropriate trial filtering

Demo 5 ✓ "Last Line of Defense" — Ethan, 12yo, CAR-T
  → ELIANA (tisagenlecleucel) matched
  → ZUMA-4 (axicabtagene) matched
  → Dual CD19/CD22 (NCT03448393) matched

Demo 6 ✓ "When Standard Drug Can't Be Used" — Aiden, 10yo, Medulloblastoma
  → SJMB12 (St. Jude reduced-dose CSI) matched
  → Pediatric MATCH SHH arm
  → Arsenic trioxide GLI1 trial
```

**Right — 4 of 6 = most appearances of any non-Oncology agent**
```
Trial Matching Across All Patient Stories:
[Evelyn] 3 trials | [Aurora] 3 trials | [Ethan] 3 trials | [Aiden] 3 trials
12 total trial matches across 4 demos
All with real NCT IDs from ClinicalTrials.gov
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

Green (#76B900), white text:
```
DGX Spark (GB10) $4,699 | Milvus 2.4, 14 collections | BGE-small 384-dim | Claude Sonnet 4 | Ports 8128, 8538 | 769 tests
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
COG: AALL0932, AALL1231, ANBL1232, ARET0321 | Pediatric MATCH (NCI-COG) | INFORM Registry
ClinicalTrials.gov API v2 | ICH E6(R2) GCP | FDA 21 CFR Part 11 | CDISC SDTM/ADaM
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Complete trial lifecycle** — Design through submission, 19 workflows covering every operational need. The lifecycle visual must show this breadth.
2. **6 decision support engines** — Confidence Calibrator, Protocol Complexity, Enrollment Predictor, Eligibility Analyzer, Competitive Threat, Historical Success. Quantitative, not subjective.
3. **13 therapeutic areas, 40 landmark trials** — not just oncology. Cardiology, neurology, rare disease, immunology, and more. Real trial IDs.
4. **14 Milvus collections** — protocols, eligibility, endpoints, sites, investigators, results, regulatory, literature, biomarkers, safety, RWE, adaptive, guidelines, genomic.
5. **Pediatric trial intelligence** — COG protocols (AALL0932, AALL1231, ANBL1232, ARET0321), Pediatric MATCH, rolling-6 design, EFS endpoints, consent/assent requirements.
6. **8 cross-agent connections** — the most of any non-orchestrator agent. Oncology, PGx, Cardiology, Biomarker, Rare Disease, Neurology, Single-Cell, Imaging.
7. **4 of 6 demos** — more appearances than any non-Oncology agent. 12 total trial matches with real NCT IDs.
8. **Patient matching with criterion-level detail** — not just "match/no match" but green/amber/red per eligibility criterion.
9. **769 tests** — second-most tested agent in the platform.
10. **Dense enough for a clinical operations war room, precise enough for an FDA pre-IND meeting, and accessible enough to match a child with the trial that could save their life.**
