# Nano Banana Pro — Precision Oncology Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Precision Oncology Agent — the central orchestrator of the HCLS AI Factory's Precision Intelligence Network. It is the agent that coordinates all other agents for cancer care, with pediatric oncology (11 cancer types including B-ALL, neuroblastoma, medulloblastoma) as the primary use case. It connects to 8 peer agents and serves as the hub for Demo 1 ("DNA to Drug") and Demo 2 ("The 30-Second Tumor Board"). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9). Same visual language as the CAR-T Intelligence Agent and Imaging Intelligence Agent infographics in this series. White canvas, structured horizontal bands, data sources on left, processing in center, outputs on right, cross-agent strip at bottom.

**Canvas:** White background (#FFFFFF). Dense but organized. Every section carries information. The diagram should feel like a molecular tumor board reference poster — comprehensive enough for a genomic oncologist, clean enough for a healthcare executive.

**Typography:**
- Title: Large, bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Section headers: Bold, navy, thin green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body: Dark gray (#333333), 8-10pt, SHORT PHRASES only
- Metric badges: White text on colored pills
- All text MUST be legible — no garbled text

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar, metric pills
- Navy: #1B2333 — title, headers, dark bars
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Red: #DC2626 — critical/resistance findings, high-risk indicators
- Amber: #F5A623 — warning indicators, moderate-risk badges
- Purple: #7B2D8E — clinical trials, guidelines collections
- Blue: #2196F3 — variant/literature collections
- Orange: #FF9800 — therapy/pathway collections
- Emerald: #059669 — favorable findings, safe indicators
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome)
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for critical alerts
- Color-coded collection badges
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Precision Oncology Agent" (green #76B900, white text)
- "Cross-Agent Orchestrator" (teal #1AAFCC, white text)

**Center:**
- Line 1 (large, bold, navy): **Precision Oncology Agent**
- Line 2 (medium, gray): "Molecular Tumor Board Intelligence on NVIDIA DGX Spark"
- Line 3 (teal): "Cross-Agent Orchestrator — Hub of the Precision Intelligence Network"
- Line 4 (small, gray): "GB10 Superchip | 128 GB Unified Memory | 28 Cancer Types | 11 Pediatric Cancers | 40+ Actionable Targets"

**Right side — Orchestrator badge (prominent, unique to this agent):**
```
ORCHESTRATOR HUB
━━━━━━━━━━━━━━━━━━━
Coordinates 8 peer agents
/v1/onco/integrated-assessment
Demo 1: "DNA to Drug"
Demo 2: "30-Second Tumor Board"
```
Teal (#1AAFCC) border, white background. This badge signals that this agent is different from the others — it's the coordinator.

---

### ━━━ LEFT COLUMN: DATA SOURCES & INGESTION (stacked vertically) ━━━

**8 data source cards with thin-line icons, count badges, and ingestion pipeline labels:**

1. **CIViC Database** [evidence icon]
   Curated variant evidence
   Clinical interpretations
   [civic_parser.py]

2. **ClinicalTrials.gov** [trial icon]
   Active oncology trials
   Phase I-IV, pediatric filters
   [clinical_trials_parser.py]

3. **OncoKB** [star icon]
   MSK precision oncology KB
   Levels of evidence (1-4, R1-R2)
   [oncokb_parser.py]

4. **PubMed / Literature** [journal icon]
   Oncology literature corpus
   Pediatric cancer publications
   [literature_parser.py]

5. **NCCN / Clinical Guidelines** [guidelines icon]
   Treatment guidelines
   Risk stratification protocols
   [guideline_parser.py]

6. **Pathway Databases** [network icon]
   KEGG, Reactome, WikiPathways
   Signaling cascades
   [pathway_parser.py]

7. **Resistance Mechanisms** [shield-break icon]
   Known resistance mutations
   Bypass pathway activations
   [resistance_parser.py]

8. **Real-World Outcomes** [chart icon]
   Treatment outcomes data
   Survival curves, response rates
   [outcome_parser.py]

**Below data sources — Ingestion Pipeline label:**
```
9 Specialized Parsers → BGE-small-en-v1.5 → 384-dim Embeddings → Milvus
```

---

### ━━━ CENTER-TOP: MILVUS VECTOR DATABASE — 11 COLLECTIONS ━━━

**Rounded container:**

```
Milvus 2.4 | Vector Database — 11 Collections | IVF_FLAT / COSINE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**Collection badges in colored rectangles (each a small card):**

Row 1:
[onco_literature] [onco_trials] [onco_variants] [onco_biomarkers] [onco_therapies]
  (blue)            (purple)      (blue)           (green)           (orange)

Row 2:
[onco_pathways] [onco_guidelines] [onco_resistance] [onco_outcomes] [onco_cases] [genomic_evidence]
  (orange)        (purple)          (red)              (green)         (amber)       (emerald — shared)

**Below — Cancer Type Coverage badges (28 types in a wrapped row):**
```
28 Cancer Types:
B-ALL | T-ALL | AML | CML | NSCLC | SCLC | Breast | Colorectal | Pancreatic | Melanoma |
Glioblastoma | Medulloblastoma | Neuroblastoma | Rhabdomyosarcoma | Wilms | Ewing |
Retinoblastoma | Ovarian | Prostate | Bladder | HCC | RCC | HNSCC | Gastric | Esophageal |
Cholangiocarcinoma | Thyroid | Other
```
Pediatric types highlighted with emerald (#059669) text or underline.

---

### ━━━ CENTER: PROCESSING ARCHITECTURE ━━━

**Three-layer processing stack:**

#### Layer 1: Query & Evidence Retrieval

```
User Query                    Evidence Retrieval              Knowledge Synthesis
"8yo female, B-ALL,    →     Parallel search across    →    Claude Sonnet 4
ETV6-RUNX1+, IKZF1          11 collections                  Streaming RAG
deletion. Risk?"              Workflow-weighted boost         Evidence-cited responses
                              Pediatric collection 2x        Deterministic-probabilistic split
```

#### Layer 2: Molecular Tumor Board Intelligence (the core)

**5 MTB workflow boxes in a row:**

```
Molecular           Therapy            Resistance         Trial              Risk
Classification  →   Ranking        →   Prediction     →   Matching       →   Stratification
                    
Variant → Gene      Evidence-level     Known escape       NCT matching       NCI/COG/INRG
→ Protein →         ranking            mutations          Biomarker-driven   criteria
Pathway → Disease   NCCN/ESMO          Bypass pathways    Phase filtering    Pediatric-specific
                    guideline-linked   Combination        Age-appropriate    IKZF1plus
                                       strategies                            MRD-directed
```

**Central callout box (large, navy border, prominent):**
```
40+ Actionable Targets
━━━━━━━━━━━━━━━━━━━━━━
Adult: BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK, RET, MET, FGFR, PIK3CA,
       IDH1/2, BRCA1/2, TP53, PTEN, CDKN2A, STK11, ESR1, MSI-H, POLE...
Pediatric [NEW]: ETV6-RUNX1, BCR-ABL1, NOTCH1, MYCN, PTCH1, PAX3-FOXO1,
                 EWSR1-FLI1, WT1, CTNNB1, RB1, KMT2A
```
Pediatric targets highlighted in emerald (#059669).

#### Layer 3: Pediatric Dosing Intelligence [NEW]

```
PEDIATRIC_DOSING Knowledge Base
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
BSA-scaled dosing (not fixed adult)
Anthracycline limit: 450 mg/m² lifetime
Vincristine cap: 2 mg absolute
MTX intrathecal: age-adjusted (8/10/12 mg)
Dexrazoxane at >300 mg/m²
```
Emerald (#059669) border, light green background.

---

### ━━━ CENTER-RIGHT: CROSS-AGENT ORCHESTRATION (the unique feature) ━━━

**This is what makes the Oncology Agent different from all others. It's the HUB.**

**Large, prominent box with teal (#1AAFCC) border and subtle glow effect:**

```
Cross-Agent Orchestrator — 8 Peer Agents
/v1/onco/integrated-assessment
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**8 agent connection cards radiating outward (or stacked), each with an icon, name, port, and purpose:**

```
[T-cell]    CAR-T Agent (:8522)
            Immunotherapy eligibility, CRS risk, CD19/CD22 targets

[test tube] Biomarker Agent (:8529)
            CD19 expression, MRD status, prognostic panel

[heart]     Cardiology Agent (:8126)
            Anthracycline cardiotoxicity, dexrazoxane, baseline LVEF

[brain]     Neurology Agent (:8528)
            MTX leukoencephalopathy, vincristine neuropathy

[pill]      PGx Agent (:8107)
            TPMT/NUDT15 dosing, CYP3A5 vincristine metabolism

[scan]      Imaging Agent (:8524)
            Staging requirements, response criteria, RANO/RECIST

[cluster]   Single-Cell Agent (:8540)
            TME profiling, blast immunophenotyping, escape clones

[clipboard] Clinical Trial Agent (:8538)
            COG protocol matching, Pediatric MATCH, age filters
```

**Integration flow label:**
```
query_*_agent() → collect responses → integrate_cross_agent_results() → unified assessment
Timeout: 30s per agent | Graceful degradation | Evidence-cited
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

**Output cards stacked vertically:**

1. **MTB Packet** [document icon]
   Complete molecular tumor board report
   Variant classification + therapy ranking
   Trial matches + risk stratification

2. **Integrated Assessment** [network icon] [NEW]
   /v1/onco/integrated-assessment
   Multi-agent unified recommendation
   8-agent coordinated output

3. **Therapy Ranking** [ranked list icon]
   Evidence-level ranked
   NCCN/ESMO guideline-linked
   Resistance-aware

4. **Trial Matches** [clipboard icon]
   NCT-linked trials
   COG pediatric protocols
   Age-appropriate filtering

5. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4
   Phenopacket (rare disease)

6. **Streamlit UI** [screen icon]
   Port :8526
   5 Tabs: Case Workbench | Evidence Explorer
   Trial Finder | Therapy Ranker | Outcomes Dashboard
   Demo patient loader
   28 cancer types in selectbox

7. **FastAPI REST** [api icon]
   Port :8527
   /api/ask | /api/cases | /api/trials
   /v1/onco/integrated-assessment

---

### ━━━ BOTTOM STRIP: DEMO WORKFLOWS ━━━

**Header (navy background, white text):** "Precision Oncology Agent in 6 Pediatric Demos"

**6 demo cards in a row, showing which demos the Oncology Agent participates in:**

```
Demo 1 ✓           Demo 2 ✓           Demo 3 ✓           Demo 4 ✓           Demo 5              Demo 6 ✓
"DNA to Drug"       "30-Second         "Protecting        "One Gene,         "Last Line          "When Standard
                    Tumor Board"        the Survivor"      One Family"        of Defense"         Drug Can't Be Used"
Evelyn, 8, ALL      Evelyn, MRD+       Marcus, 6, NB      Aurora, 4, RB      Ethan, 12          Aiden, 10, MB
                                                                              (CAR-T focused)
[ORCHESTRATOR]      [ORCHESTRATOR]     Dose adjustment    RB management      —                   SHH analysis
                                                                                                  Vismodegib flag
```

Green checkmarks (✓) on the 5 demos where Oncology participates. The most of any agent.

---

### ━━━ INFRASTRUCTURE BAR ━━━

**Green (#76B900) background, white text:**

```
DGX Spark (GB10)    Milvus 2.4         BGE-small-en-v1.5   Claude Sonnet 4    Docker Compose
$4,699              11 collections      384-dim vectors      Streaming RAG      Ports 8526, 8527
128 GB unified      IVF_FLAT/COSINE    Async embedding      + Pediatric prompt 9 ingest parsers
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
Cross-Agent Orchestrator: Coordinates CAR-T, Biomarker, Cardiology, Neurology, PGx, Imaging, Single-Cell, Clinical Trial
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **This is the HUB agent** — it orchestrates 8 peer agents, more than any other. The cross-agent box must be visually prominent.
2. **28 cancer types** — both adult and pediatric, with 11 pediatric types highlighted. Most comprehensive oncology coverage in the platform.
3. **40+ actionable targets** — from BRAF/EGFR (adult) to ETV6-RUNX1/MYCN (pediatric). All visible.
4. **11 Milvus collections** — literature, trials, variants, biomarkers, therapies, pathways, guidelines, resistance, outcomes, cases, genomic_evidence.
5. **9 specialized parsers** — CIViC, OncoKB, ClinicalTrials.gov, PubMed, NCCN, pathways, resistance, outcomes, literature.
6. **Pediatric dosing is built in** — BSA-scaled dosing, vincristine 2mg cap, MTX intrathecal, anthracycline limits.
7. **5-tab Streamlit UI** — Case Workbench with variant entry, Trial Finder, Therapy Ranker, Evidence Explorer, Outcomes Dashboard.
8. **Participates in 5 of 6 demos** — the most of any agent. Orchestrates Demo 1 and Demo 2.
9. **It runs on one $4,699 DGX Spark** — entire molecular tumor board capability on a desktop.
10. **Dense enough for a genomic oncology laboratory wall, clear enough to present at ASCO.**
