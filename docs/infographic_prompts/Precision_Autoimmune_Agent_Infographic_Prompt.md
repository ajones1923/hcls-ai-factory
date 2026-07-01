# Nano Banana Pro — Precision Autoimmune Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Precision Autoimmune Agent — the HCLS AI Factory's immunology intelligence engine featuring 50+ HLA-disease associations, 24 autoantibody-disease mappings, 22 biologic therapies, 20 disease activity scoring thresholds, flare prediction biomarker patterns, a diagnostic engine with ACR/EULAR classification criteria, a diagnostic odyssey timeline builder, clinical document ingestion, and a 10-tab Streamlit UI — the most tabs of any agent in the platform. It connects to 6 peer agents and plays critical roles in Demo 3 ("Protecting the Survivor" — dinutuximab irAEs) and Demo 5 ("Last Line of Defense" — CRS/irAE monitoring post-CAR-T). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as all other HCLS AI Factory agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a rheumatology/immunology department reference poster — comprehensive enough for a clinical immunologist, relevant enough for an immuno-oncologist managing checkpoint inhibitor toxicity, clean enough for a hospital executive.

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
- Red: #DC2626 — critical flare alerts, severe disease activity
- Amber: #F5A623 — moderate disease activity, watch indicators
- Emerald: #059669 — remission, low disease activity, treatment response
- Purple: #7B2D8E — HLA associations, autoantibody patterns
- Blue: #2196F3 — biologic therapies, trial collections
- Orange: #FF9800 — immunotherapy irAE indicators
- Pink: #EC4899 — flare prediction, diagnostic odyssey
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — immune cell, antibody, shield, DNA, flame motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for flare alerts
- Color-coded disease activity (emerald/amber/red for remission/moderate/severe)
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Precision Autoimmune Agent" (green #76B900, white text)
- "10-Tab UI | 6 Peer Agents" (teal #1AAFCC, white text)

**Center:**
- Line 1 (large, bold, navy): **Precision Autoimmune Agent**
- Line 2 (gray): "Immunology Intelligence & Immune-Related Adverse Event Profiling"
- Line 3 (teal): "Part of the Precision Intelligence Network — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB | 50+ HLA Associations | 24 Autoantibodies | 22 Biologics | 20 Activity Scores"

**Right side — Clinical capability badge:**
```
IMMUNOLOGY ENGINE
━━━━━━━━━━━━━━━━━━━━━
50+ HLA-Disease Associations
24 Autoantibody-Disease Maps
22 Biologic Therapies
20 Disease Activity Scores
Flare Prediction Biomarkers
Diagnostic Odyssey Timeline
ACR/EULAR Classification
Demo 3 + Demo 5
```

---

### ━━━ LEFT COLUMN: IMMUNOLOGY KNOWLEDGE DOMAINS ━━━

**Stacked input cards:**

1. **HLA-Disease Associations** [DNA-shield icon]
   50+ HLA allele-disease pairs
   HLA-B*27:05 → Ankylosing Spondylitis
   HLA-DRB1*04:01 → Rheumatoid Arthritis
   HLA-DQB1*02:01 → Celiac Disease
   HLA-C*06:02 → Psoriasis
   HLA-B*51:01 → Behcet's Disease
   [50+ associations] (purple badge)

2. **Autoantibody Database** [antibody icon]
   24 autoantibody-disease mappings
   ANA, dsDNA, anti-CCP, ANCA
   Anti-Jo-1, anti-Scl-70, anti-centromere
   Pattern recognition algorithms
   [24 autoantibodies] (purple badge)

3. **Biologic Therapies** [syringe icon]
   22 biologic agents profiled
   TNF inhibitors, IL-6, IL-17, IL-23
   B-cell depletion (rituximab)
   JAK inhibitors, CTLA-4, PD-1/PD-L1
   [22 biologics] (blue badge)

4. **Disease Activity Scores** [gauge icon]
   20 validated scoring thresholds
   DAS28, CDAI, SDAI (RA)
   SLEDAI, BILAG (SLE)
   BASDAI (SpA), PASI (Psoriasis)
   [20 scores] (teal badge)

5. **Flare Biomarker Patterns** [flame icon]
   Predictive biomarker signatures
   CRP, ESR, complement, cytokine panels
   Pre-flare trajectory detection
   [13 flare patterns] (amber badge)

6. **Clinical Documents** [document icon]
   PDF clinical note ingestion
   Diagnostic odyssey extraction
   Temporal event timeline
   [document_processor.py]

7. **Immunotherapy irAEs** [warning-shield icon] [NEW]
   Checkpoint inhibitor toxicity
   CAR-T CRS/ICANS monitoring
   Dinutuximab adverse events
   [pediatric irAE profiling]

Arrows rightward: "BGE-small-en-v1.5 | 384-dim | Milvus Collections"

---

### ━━━ CENTER-TOP: DIAGNOSTIC ENGINE (unique feature) ━━━

**Section header (navy, green underline):** "Diagnostic Engine — ACR/EULAR Classification & Differential Diagnosis"

**Three processing boxes in a row:**

```
[puzzle]                        [timeline]                      [overlap]
ACR/EULAR Classification        Diagnostic Odyssey              Cross-Disease Overlap
━━━━━━━━━━━━━━━━━━━━━━          ━━━━━━━━━━━━━━━━━               ━━━━━━━━━━━━━━━━━━━━
Classification criteria          Timeline construction          Lupus-Sjogren overlap
evaluation for 13 conditions    from clinical documents        RA-SLE overlap detection
Probabilistic scoring            Temporal event extraction      Shared autoantibody patterns
Criteria met / not met           Years-to-diagnosis tracking    Undifferentiated CTD
                                 Diagnostic delay quantification
```

**Below — HLA-Autoantibody Integration visual:**
```
HLA Typing → Autoantibody Panel → Disease Activity Score → Flare Risk → Therapy Selection

Example flow:
HLA-DRB1*04:01 + Anti-CCP+ → RA confirmed (ACR/EULAR 2010)
→ DAS28-CRP 5.1 (High activity) → Flare risk: elevated CRP trend
→ Therapy: Methotrexate → Adalimumab → JAK inhibitor escalation
```

---

### ━━━ CENTER: AUTOANTIBODY PATTERN RECOGNITION + DISEASE ACTIVITY ━━━

**Two side-by-side sections:**

#### Left half: Autoantibody Pattern Map

```
Autoantibody Pattern Recognition
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Antibody          → Disease(s)                    Specificity
ANA (homogeneous) → SLE, drug-induced lupus       Moderate
Anti-dsDNA        → SLE (active nephritis)        High
Anti-CCP          → Rheumatoid Arthritis           Very High
c-ANCA/PR3        → GPA (Wegener's)               High
p-ANCA/MPO        → MPA, EGPA                     Moderate
Anti-Jo-1         → Antisynthetase syndrome        High
Anti-Scl-70       → Diffuse systemic sclerosis     High
Anti-centromere   → Limited scleroderma (CREST)    High
Anti-SSA/Ro       → Sjogren's, neonatal lupus      Moderate
Anti-SSB/La       → Sjogren's                      High
Anti-phospholipid → APS                            Moderate
Anti-TTG/EMA      → Celiac disease                 Very High
```

#### Right half: Disease Activity Dashboard

```
Disease Activity Scoring
━━━━━━━━━━━━━━━━━━━━━━━━

Disease     Score      Remission   Low        Moderate     High
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RA          DAS28      <2.6        2.6-3.2    3.2-5.1      >5.1
RA          CDAI       ≤2.8       2.9-10     10.1-22      >22
SLE         SLEDAI     0           1-5        6-10         >10
SLE         BILAG      All D/E     1C         1B           1A
SpA         BASDAI     <2          2-4        4-6          >6
Psoriasis   PASI       <3          3-10       10-20        >20
Vasculitis  BVAS       0           1-9        10-19        ≥20
IBD         Harvey-BI  <5          5-7        8-16         >16

[Emerald]    [Emerald]  [Amber]     [Red]
```
Color-coded: remission/low = emerald, moderate = amber, high = red

---

### ━━━ CENTER-BOTTOM: PEDIATRIC irAE + CROSS-AGENT ━━━

**Immunotherapy irAE Profiling box (orange #FF9800 border — signals immunotherapy context):**
```
IMMUNOTHERAPY ADVERSE EVENT PROFILING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Checkpoint Inhibitor irAEs (Pediatric Oncology):
• Autoimmune colitis: 10-20% | Diarrhea grade 3-4
• Hepatitis: 5-10% | Monitor ALT/AST
• Pneumonitis: 3-5% | Dyspnea, ground-glass CT
• Myocarditis: 1-2% | 25-50% mortality if grade 3-4
• Autoimmune encephalitis: <1% | Altered mental status

Dinutuximab (Anti-GD2) irAEs — Demo 3 Marcus:
• Neuropathic pain: 85% | Morphine infusion protocol
• Capillary leak syndrome: 25% | Daily weights, albumin
• Hypersensitivity: 15% | Premedication protocol

CAR-T CRS Monitoring — Demo 5 Ethan:
• CRS: 77% any grade, 47% grade 3-4
• Autoimmune cytopenias: 30-40% post-CAR-T
• B-cell aplasia: 100% (lifelong IVIG)
```

**Cross-Agent Integration box (teal border):**
```
Cross-Agent Integration (/v1/autoimmune/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)     irAE risk from cancer immunotherapy
→ Cardiology Agent (:8126)   Myocarditis risk from checkpoint inhibitors
→ Neurology Agent (:8528)    Autoimmune encephalitis from immunotherapy
→ PGx Agent (:8107)          Drug-gene interactions for biologics (TPMT, NUDT15)
→ Biomarker Agent (:8529)    Autoantibody-disease activity correlation
→ Clinical Trial Agent (:8538) Autoimmune/immunotherapy trial matching

6 peer agents — most cross-agent connections of any specialty agent
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

**Output cards:**

1. **Diagnostic Report** [document icon]
   ACR/EULAR classification result
   Differential diagnosis ranked
   Autoantibody pattern interpretation

2. **Disease Activity Dashboard** [gauge icon]
   Current activity score + trend
   Color-coded severity
   Flare risk prediction

3. **Diagnostic Odyssey Timeline** [timeline icon]
   Chronological event visualization
   Years-to-diagnosis tracking
   Missed diagnosis identification

4. **Therapy Recommendation** [pill icon]
   Biologic therapy selection
   Escalation pathway
   PGx-guided dosing

5. **Integrated Assessment** [network icon]
   /v1/autoimmune/integrated-assessment
   6-agent coordinated output

6. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4

**Streamlit UI (prominent — 10 tabs, most of any agent):**
```
Streamlit UI (:8531) — 10 Tabs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 1. Clinical Query (RAG Q&A)
 2. Patient Analysis (full workup)
 3. Document Ingest (PDF upload)
 4. Diagnostic Odyssey (timeline)
 5. Autoantibody Panel (interpretation)
 6. HLA Analysis (allele-disease lookup)
 7. Disease Activity (scoring dashboard)
 8. Flare Prediction (biomarker trends)
 9. Therapy Advisor (biologic selection)
10. Knowledge Base (collection explorer)

10 TABS — most of any agent
```

**FastAPI REST:**
```
Port :8532
/query | /search | /analyze | /differential
/ingest/upload | /ingest/demo-data | /export
/v1/autoimmune/integrated-assessment
```

---

### ━━━ BOTTOM STRIP: DEMO PARTICIPATION ━━━

**Header (navy, white text):** "Precision Autoimmune Agent in HCLS AI Factory"

**Left — Demo participation:**
```
Demo 3 ✓ "Protecting the Survivor" — Marcus, 6yo, Neuroblastoma
  → Dinutuximab (anti-GD2) irAE profiling
  → Neuropathic pain 85%, capillary leak 25%, hypersensitivity 15%
  → Cross-agent flags to Cardiology + Neurology

Demo 5 ✓ "Last Line of Defense" — Ethan, 12yo, CAR-T
  → CRS cytokine cascade monitoring (IL-6, IFN-gamma, TNF-alpha)
  → Autoimmune cytopenias: 30-40% post-CAR-T
  → B-cell aplasia management (lifelong IVIG replacement)
```

**Right — Unique capabilities visual:**
```
What Makes This Agent Unique:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[timeline] Diagnostic Odyssey — Tracks years of misdiagnosis
[upload]   Document Ingest — Reads clinical PDFs directly
[flame]    Flare Prediction — Pre-flare biomarker trajectories
[overlap]  Cross-Disease — Detects overlapping autoimmune conditions
[shield]   irAE Profiling — Checkpoint + CAR-T + anti-GD2 toxicity
[10 tabs]  Most Comprehensive UI — 10 specialized tabs
[6 agents] Most Connected Specialty — 6 cross-agent integrations
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

Green (#76B900), white text:
```
DGX Spark (GB10) $4,699 | Milvus 2.4 | BGE-small 384-dim | Claude Sonnet 4 | Ports 8531, 8532 | 200+ tests
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
HLA: PMID:28622507 | ACR/EULAR: 2010/2019 Classification Criteria | irAEs: ASCO/NCCN Immunotherapy Guidelines
Dinutuximab: COG ANBL1232 | CAR-T CRS: Lee et al. Blood 2019 ASTCT Consensus Grading
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **50+ HLA-disease associations** — the deepest immunogenetics database of any agent. HLA-B27, DRB1, DQB1 pairs with specific diseases.
2. **24 autoantibody-disease mappings** — from ANA to anti-TTG, each with disease specificity and clinical interpretation.
3. **22 biologic therapies** — TNF inhibitors through JAK inhibitors, with escalation pathways.
4. **20 disease activity scoring thresholds** — DAS28, SLEDAI, BASDAI, PASI — all with color-coded severity (remission/low/moderate/high).
5. **Diagnostic Engine** — ACR/EULAR classification criteria evaluation, differential diagnosis scoring, cross-disease overlap detection. This agent diagnoses, not just informs.
6. **Diagnostic Odyssey Timeline** — unique capability that tracks years of misdiagnosis from ingested clinical documents. No other agent has this.
7. **Flare Prediction** — pre-flare biomarker trajectory detection. Identifies flares before they happen.
8. **Immunotherapy irAE profiling** — checkpoint inhibitors, dinutuximab, CAR-T CRS. Bridges autoimmunity and oncology.
9. **10-tab Streamlit UI + 6 peer agents** — most tabs and most cross-agent connections of any specialty agent.
10. **Dense enough for a rheumatology clinic wall, relevant enough for an immuno-oncology tumor board, and protective enough to monitor a child's immune system through cancer treatment.**
