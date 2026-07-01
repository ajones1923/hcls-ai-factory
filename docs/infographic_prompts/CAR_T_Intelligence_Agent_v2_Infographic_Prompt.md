# Nano Banana Pro — CAR-T Intelligence Agent on NVIDIA DGX Spark v2.0

## IMPORTANT: Read this entire prompt before generating. This is an UPDATE to the CAR-T Intelligence Agent infographic (v1.0). The original showed the CAR-T agent with 5 data sources, Milvus 2.4 with 11 collections, query expansion + comparative mode + parallel search, and RAG response outputs. v2.0 preserves that exact visual style but updates for: pediatric CAR-T as primary use case (tisagenlecleucel for pediatric ALL), cross-agent integration with 5 peer agents, PEDIATRIC_CART and PEDIATRIC_CRS_ICANS_MANAGEMENT knowledge bases, the role in Demo 2 ("30-Second Tumor Board") and Demo 5 ("Last Line of Defense"), and production-grade endpoints. ONE canvas, landscape 16:9, same dense technical poster as v1.0.

---

## REFERENCE: What v1.0 looked like (preserve this style EXACTLY)

The v1.0 had these characteristics:
- White canvas, landscape 16:9
- "CAR-T Intelligence Agent" badge top-left (green)
- Title top-center: "CAR-T Intelligence Agent" with "Cross-Functional CAR-T Intelligence on NVIDIA DGX Spark"
- "5 CAR-T Development Stages in Every Query" callout
- LEFT COLUMN: 5 data source cards stacked vertically (PubMed, ClinicalTrials.gov, FDA Products, Landmark Papers, Manufacturing)
- Each data source with a count badge (5,047 abstracts, 973 trials, etc.)
- CENTER-TOP: Milvus 2.4 box showing 11 collections with color-coded category badges (Literature, Trials, Constructs, Assays, Manufacturing, Safety, Biomarkers, Regulatory, Sequences, Real-World)
- CENTER: Three processing boxes in a row (Query Expansion → Comparative Mode → Parallel Search)
- "User Query" entry point on left of processing row
- Claude Sonnet 4 box below processing row (Anthropic API, Streaming RAG)
- RIGHT COLUMN: Output boxes (RAG Response, Comparative Analysis, Knowledge Graph, Image/Claude Vision, Streamlit UI, FastAPI REST, Report Generation, Export)
- BOTTOM STRIP: "HCLS AI Factory — Cross-Pipeline Integration" with 6 connected boxes (Genomics, RAG/Chat, Drug Discovery, Imaging Agent, Biomarker Agent, Open Source)
- VERY BOTTOM: Infrastructure bar (DGX Spark, Milvus 2.4, BGE-small, Claude Sonnet, Docker Compose)
- Dense, organized, every section labeled with metric counts

## WHAT CHANGES FROM v1.0 TO v2.0:

1. **Badge** → "v2.0 Production" (was just "CAR-T Intelligence Agent")
2. **Subtitle** → Adds "Part of the Precision Intelligence Network" and "Pediatric ALL Primary Use Case"
3. **Data sources** → EXPANDED with pediatric CAR-T literature and ELIANA trial data
4. **Collections** → 11 collections PLUS new pediatric-specific knowledge bases
5. **NEW: Pediatric CAR-T Knowledge section** — PEDIATRIC_CART, PEDIATRIC_CRS_ICANS_MANAGEMENT
6. **NEW: Cross-Agent Integration box** — connections to Biomarker, Oncology, Single-Cell, Cardiology, Clinical Trial
7. **Processing row** → Updated with pediatric eligibility assessment and CRS risk modeling
8. **Outputs** → NEW: /integrated-assessment endpoint, pediatric CRS/ICANS management protocols
9. **Cross-Pipeline strip** → Updated to show all 11 agents, not just 6 connections
10. **Demo references** → Demo 2 ("30-Second Tumor Board") and Demo 5 ("Last Line of Defense")

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background. Dense technical poster. IDENTICAL visual language to v1.0.

**Typography (matching v1.0):**
- Title: Large, bold, navy (#1B2333)
- Subtitle: Gray (#666666), smaller
- Collection badges: Small colored rectangles with white text
- Body labels: Dark gray (#333333), 8-10pt, SHORT PHRASES
- Metric counts: Bold, inside colored pills
- All text MUST be legible — no garbled small text

**Color Palette (matching v1.0):**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — data flow lines, sub-headers
- Purple: #7B2D8E — literature/research collections
- Blue: #2196F3 — trials collections
- Orange: #F5A623 — constructs/manufacturing
- Green: #4CAF50 — assays/biomarkers
- Red: #DC2626 — safety collections, critical alerts
- Emerald: #059669 — verified/cleared indicators
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements (matching v1.0):**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px)
- Color-coded collection badges (small rectangles)
- Directional arrows for data flow
- Metric count badges throughout
- Dual-path arrows where applicable

**CRITICAL TEXT RULES:**
- Use SHORT PHRASES (3-5 words) not sentences
- Metric badges for numbers, not inline text
- Agent descriptions as bullet points not paragraphs

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Rounded badge: "CAR-T Intelligence Agent" (green #76B900) + "v2.0 Production" (navy) [UPDATED]

**Center:**
- Title (large, bold, navy): **CAR-T Intelligence Agent**
- Subtitle line 1 (gray): "Cross-Functional CAR-T Intelligence on NVIDIA DGX Spark"
- Subtitle line 2 (teal): "Part of the Precision Intelligence Network — HCLS AI Factory" [NEW]
- Subtitle line 3 (gray, small): "GB10 Superchip | 128 GB Unified Memory | Pediatric ALL Primary Use Case" [UPDATED]

**Right side — Stage badges (matching v1.0 style, updated):**
```
7 CAR-T Decision Stages [UPDATED from 5]
● Target Validation    ● Eligibility Assessment
● CRS/ICANS Risk      ● Manufacturing Planning
● Cardiac Clearance   ● Trial Matching
● Bridging Therapy
```

**Far right — Data flow key (matching v1.0):**
```
● Literature (5,047)  ● Trials (973)  ● Constructs (6)
● Assays (45)  ● Manufacturing (30)  ● Safety
● Pediatric CAR-T [NEW]  ● CRS/ICANS Mgmt [NEW]
→ Data Flow  - → Cross-Agent  ●→ Comparative
```

---

### ━━━ LEFT COLUMN: DATA SOURCES (matching v1.0, expanded) ━━━

**Stacked input cards with count badges. Matching v1.0 visual style exactly:**

1. **PubMed** [journal icon]
   5,047 abstracts via NCBI E-utilities
   + Pediatric CAR-T literature [NEW badge]
   [5,047]

2. **ClinicalTrials.gov** [trial icon]
   973 trials via API v2
   + ELIANA, ZUMA-4, Dual CD19/CD22 [NEW badge]
   [973]

3. **FDA Products** [shield icon]
   6 approved CAR-T therapies
   Tisagenlecleucel (Kymriah) — Pediatric ALL [HIGHLIGHTED]
   [6 approved]

4. **Landmark Papers** [star icon]
   45 seminal studies
   + Maude et al. NEJM 2018 (ELIANA) [NEW badge]
   + Fry et al. Nature Medicine 2018 (CD22) [NEW badge]
   [45 papers]

5. **Manufacturing** [factory icon]
   30 CMC process records
   + Pediatric apheresis protocols [NEW badge]
   + 22-day manufacturing timeline [NEW badge]
   [30 records]

**NEW DATA SOURCE:**
6. **Pediatric Knowledge Bases** [child icon] [NEW — emerald border]
   PEDIATRIC_CART: CD19/CD22/CD30/GD2 targets
   PEDIATRIC_CRS_ICANS_MANAGEMENT: Weight-based dosing
   ELIANA outcomes: 82% CR, 12-mo EFS 50%, OS 76%
   [5 pediatric KBs]

Arrows flow rightward into the Milvus / processing section. Label: "Embedding via BGE-small-en-v1.5 | 384-dim vectors"

---

### ━━━ CENTER-TOP: MILVUS VECTOR DATABASE (matching v1.0, expanded) ━━━

**Rounded container matching v1.0 style:**

```
Milvus 2.4 | Vector Database — 11+ Collections
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**Collection badges in colored rectangles (matching v1.0 color-coding):**

Row 1 (original):
[Literature: 5,047] [Trials: 973] [Constructs: 6] [Assays: 45] [Manufacturing: 30]

Row 2 (original):
[Safety: 43] [Biomarkers: 43] [Regulatory: 25] [Sequences: 27] [Real-World: 8]

Row 3 [NEW — teal/emerald badges]:
[Pediatric CAR-T Targets] [Pediatric CRS/ICANS] [Pediatric Trial Outcomes]

**Below collections:**
```
IVF_FLAT | COSINE | nlist=1024
Total indexed: ~6,200+ documents across all collections
```

---

### ━━━ CENTER: PROCESSING ROW (matching v1.0 layout, updated content) ━━━

**Left entry point:**
```
User Query
Example: "Is this 12-year-old with
r/r B-ALL eligible for tisagenlecleucel?
CD19+ 97.2%, failed 2 prior lines,
LVEF 58%"
```

**Three processing boxes connected by arrows (matching v1.0):**

**Box 1: Query Expansion** (matching v1.0)
```
Query Expansion
169 synonyms → 1.4% entity recall boost
+ Pediatric aliases [NEW]
  "peds ALL" → B-ALL
  "Kymriah" → tisagenlecleucel
  "tisa-cel" → tisagenlecleucel
[~7 ms]
```

**YES ← Dual Path →**

**Box 2: Comparative Mode** [UPDATED]
```
Comparative Mode + Eligibility Assessment [UPDATED]
Dual-entity/side-by-side analysis
+ Pediatric Eligibility Checklist [NEW]:
  ☑ Age ≤25  ☑ r/r B-ALL  ☑ CD19+
  ☑ Adequate organ function  ☑ No active GVHD
  ☑ Prior lines ≥2
+ CRS Risk Modeling [NEW]:
  Tumor burden → CRS grade prediction
  Weight-based tocilizumab dosing
[~215 ms]
```

**Box 3: Parallel Search** (matching v1.0, updated)
```
Parallel Search
11+ collections searched simultaneously
+ Pediatric-weighted collection boost [NEW]:
  PEDIATRIC_CART: 2.5x weight
  PEDIATRIC_CRS_ICANS: 2.0x weight
  cart_trials: 1.5x weight (ELIANA filter)
```

---

### ━━━ CENTER-BOTTOM: CLAUDE LLM + CROSS-AGENT (below processing row) ━━━

**Claude box (matching v1.0):**
```
Claude Sonnet 4
Anthropic API | Streaming RAG + citations
Pediatric CAR-T system prompt [NEW]
<24 ms TTFT | >300 tok/sec
```

**NEW: Cross-Agent Integration box (teal border, prominent):** [NEW in v2]
```
Cross-Agent Integration (/v1/cart/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Biomarker Agent (:8529)    CD19/CD22 expression validation
→ Oncology Agent (:8527)     Tumor profile, disease context
→ Single-Cell Agent (:8540)  TME profiling, target density
→ Cardiology Agent (:8126)   Pre-lymphodepletion cardiac clearance
→ Clinical Trial Agent (:8538) ELIANA, ZUMA-4, dual CD19/CD22 trials

Timeout: 30s | Graceful degradation | integrate_cross_agent_results()
```

---

### ━━━ RIGHT COLUMN: OUTPUTS (matching v1.0, expanded) ━━━

**Matching v1.0 output boxes, with additions:**

1. **RAG Response** [chat icon]
   Streaming with citations
   NVIDIA dark theme
   + Pediatric-specific language [NEW]

2. **Comparative Analysis** [table icon]
   Side-by-side tables
   + Tisagenlecleucel vs Axi-cel vs CD22 [NEW]
   + Pediatric vs Adult outcomes [NEW]

3. **Eligibility Report** [checklist icon] [NEW]
   7-point eligibility checklist
   CRS risk grade prediction
   Tocilizumab dosing calculation
   CD22 backup strategy

4. **Knowledge Graph** [network icon]
   Interactive network visualization

5. **Report Generation** [document icon]
   Markdown | JSON | PDF | FHIR
   + Pediatric CRS monitoring protocol [NEW]

6. **Streamlit UI** [screen icon]
   Chat interface (:8521)
   Deep Research mode
   Collection filtering
   Citation relevance scoring

7. **FastAPI REST** [api icon]
   /health | /query | /search
   /v1/cart/integrated-assessment [NEW]
   Port 8522

---

### ━━━ BOTTOM STRIP: CROSS-PIPELINE INTEGRATION (updated from v1.0) ━━━

**Header (navy background, white text):** "HCLS AI Factory — Cross-Agent Integration (11 Agents)"

**Updated from v1.0's 6 boxes to show specific CAR-T integration paths:**

1. **CAR-T → Genomics** — Variant context for target antigen genes, germline screening
2. **CAR-T → Oncology** — Tumor molecular profile, treatment history, disease context
3. **CAR-T → Biomarker** — CD19/CD22 expression quantification, MRD monitoring
4. **CAR-T → Single-Cell** — TME profiling, blast immunophenotyping, escape clone detection
5. **CAR-T → Cardiology** — Pre-lymphodepletion cardiac clearance, CRS cardiac monitoring
6. **CAR-T → Clinical Trial** — ELIANA, ZUMA-4, dual CD19/CD22, Pediatric MATCH
7. **CAR-T → Drug Discovery** — Bridging therapy candidates for 22-day manufacturing window

Dashed teal arrows connecting each box.

**Demo reference badges:**
[Demo 2: "30-Second Tumor Board"] [Demo 5: "Last Line of Defense" — Ethan, 12yo]

---

### ━━━ INFRASTRUCTURE BAR (matching v1.0) ━━━

**Green (#76B900) background, white text:**

```
DGX Spark          Milvus 2.4        BGE-small          Claude Sonnet 4     Docker Compose
GB10, $4,699       IVF_FLAT/COSINE   384-dim vectors     Streaming RAG       6 services
                   11+ collections    Async embedding     + Pediatric prompt  Ports 8521, 8522
```

---

### ━━━ BOTTOM FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org | Part of HCLS AI Factory v2.0 | March 2026
Pediatric CAR-T: Tisagenlecleucel 82% CR (ELIANA) | CD22 73% CR (Fry 2018) | GD2 Neuroblastoma (Investigational)
```

---

## WHAT THIS v2.0 MUST COMMUNICATE vs v1.0

| v1.0 | v2.0 |
|------|------|
| Generic CAR-T agent | Pediatric ALL as primary use case |
| 5 data sources | 6 data sources + pediatric knowledge bases |
| 11 collections | 11+ collections with pediatric-specific additions |
| No peer agent connections | 5 cross-agent integrations with specific ports |
| 5 CAR-T stages | 7 CAR-T decision stages (added cardiac, bridging) |
| Generic query example | "12-year-old with r/r B-ALL, CD19+ 97.2%" |
| No eligibility checklist | Pediatric eligibility 7-point checklist with CRS risk |
| 6 cross-pipeline boxes | 7 cross-agent boxes with specific integration paths |
| No demo references | Demo 2 + Demo 5 badges |
| Generic comparative mode | + Tisagenlecleucel vs CD22 vs dual targeting |

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Same architecture as v1.0** — recognizably the same poster, naturally evolved
2. **Pediatric ALL is the primary use case** — Ethan's case, tisagenlecleucel, ELIANA data
3. **Pediatric knowledge bases are new** — PEDIATRIC_CART, PEDIATRIC_CRS_ICANS prominently shown
4. **5 peer agents connected** — Biomarker, Oncology, Single-Cell, Cardiology, Clinical Trial
5. **7 decision stages** — expanded from 5 to include cardiac clearance and bridging therapy
6. **ELIANA outcomes are specific** — 82% CR, 12-mo EFS 50%, OS 76% visible as badges
7. **CRS risk modeling is built in** — weight-based tocilizumab, grade prediction
8. **CD22 backup strategy** — antigen escape addressed
9. **Still runs on one DGX Spark** — $4,699 machine
10. **Dense enough for a CAR-T manufacturing facility wall, clear enough for a pediatric oncologist**
