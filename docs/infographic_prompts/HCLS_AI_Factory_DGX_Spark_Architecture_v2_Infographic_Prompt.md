# Nano Banana Pro — HCLS AI Factory on NVIDIA DGX Spark Architectural Infographic v2.0

## IMPORTANT: Read this entire prompt before generating. This is an UPDATE to an existing infographic (v1.0). The original showed 3 pipelines with 5 agents. This v2.0 reflects the current platform: 3 RENAMED engines, 11 intelligence agents in 4 domain groups, pediatric oncology as primary use case, cross-agent integration, and pediatric safety filters in drug discovery. The visual style, layout density, and architectural poster aesthetic of v1.0 must be PRESERVED — this should look like a natural evolution, not a redesign. ONE canvas, landscape 16:9, reference architecture poster density.

---

## REFERENCE: What v1.0 looked like (preserve this style)

The original v1.0 had these visual characteristics that MUST be maintained:
- White canvas with structured horizontal bands
- Data sources on the LEFT edge (FASTQ, GRCh38, ClinVar, AlphaMissense, RCSB PDB, Cryo-EM)
- Three pipeline bands stacked vertically in the CENTER (Genomics → RAG → Drug Discovery)
- Output consumers on the RIGHT edge (Researcher, Clinician, Med Chemist)
- Milvus Vector Database as a shared layer between RAG and Drug Discovery
- NVIDIA infrastructure bar at the BOTTOM
- Performance metrics badge in bottom-left
- Color-coded pipeline bands (green for genomics, blue/teal for RAG, purple for drug discovery)
- Small process-flow boxes connected by arrows within each pipeline
- Dense but organized — every section carries information
- Warm, professional aesthetic with subtle background tints per pipeline band

## WHAT CHANGES FROM v1.0 TO v2.0:

1. **Pipeline names** → Engine names: "GPU-Accelerated Genomics" → "Genomic Foundation Engine" | "RAG/Chat Evidence Retrieval" → "Precision Intelligence Network" | "Generative Drug Discovery" → "Therapeutic Discovery Engine"
2. **5 agents** → **11 agents in 4 domain groups** shown explicitly in the Precision Intelligence Network band
3. **Demo use case** → "Pediatric Oncology" (was "VCP & Frontotemporal Dementia")
4. **New: Cross-agent integration** lines between agents
5. **New: Pediatric Safety Filter** added to drug discovery pipeline
6. **New: 6 Demo Workflows** referenced
7. **Updated metrics** in performance badge
8. **New: Agent output consumers** (Clinical Trial Matching, Cardiotoxicity Prevention, etc.)

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual language as v1.0.

**Typography:**
- Title: Large, bold, sans-serif (Inter/Helvetica), navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Pipeline/Engine names: Bold, colored headers matching their band tint
- Process boxes: Bold label (navy, 9pt) + detail text (gray #666666, 8pt)
- Metric badges: White text on green (#76B900) or teal (#1AAFCC) rounded pills
- All text MUST be legible — use short phrases, not sentences

**Color Palette (same as v1.0):**
- NVIDIA Green: #76B900 — Genomic Foundation Engine band, infrastructure bar, metric badges
- Teal/Blue: #1AAFCC — Precision Intelligence Network band, agent connections
- Purple: #7B2D8E — Therapeutic Discovery Engine band
- Navy: #1B2333 — Title, headers, text
- Light Gray: #F5F5F5 — Card backgrounds
- White: #FFFFFF — Canvas
- Amber: #F5A623 — Clinical alert indicators
- Red: #DC2626 — Critical finding markers
- Emerald: #059669 — Safety/verified indicators
- Gray: #666666 — Secondary text

**Visual Elements (matching v1.0):**
- Rounded-corner rectangles (8px radius) for all process boxes
- Thin-line icons (16-24px, monochrome) for data sources and outputs
- Directional arrows: solid colored for primary flow, dashed for cross-connections
- Pipeline bands with subtle background color tints
- Metric badge pills throughout
- NVIDIA logo in title and infrastructure bar

---

## CANVAS STRUCTURE (Top to Bottom)

### ━━━ BAND 1: TITLE BAR ━━━

**Left side — Demo Use Case Badge (updated):**
```
Demo Use Case Badge
━━━━━━━━━━━━━━━━━━━━━━━
Pediatric Oncology
5 Patient Stories | 6 Demo Workflows
Evelyn (ALL) • Marcus (Neuroblastoma)
Aurora (Retinoblastoma) • Ethan (CAR-T)
Aiden (Medulloblastoma)
```
Light green (#E8F5E9) background, thin green border, small text.

**Center:**
- Line 1 (large, bold, navy): **Healthcare & Life Sciences AI Factory**
- Line 2 (medium, teal): **Precision Medicine Platform on NVIDIA DGX Spark**
- Line 3 (small, gray): Patient DNA → Drug Candidates in <5 Hours | Apache 2.0 Open Source
- NVIDIA logo (right of title, matching v1.0 placement)

**Right side — Pipeline Key (matching v1.0 style):**
```
━━━ Genomic Foundation Engine
━━━ Precision Intelligence Network
━━━ Therapeutic Discovery Engine
→   Data Flow
- → Cross-Agent Coordination
```
Color-coded lines matching each engine's band color.

---

### ━━━ LEFT COLUMN: DATA SOURCES (spanning vertically alongside Bands 2-4) ━━━

**Matching v1.0 layout exactly — stacked input cards on the left edge:**

1. **FASTQ** [DNA icon]
   Raw Sequencing
   Tumor + Germline

2. **GRCh38** [genome icon]
   Reference Genome
   Parabricks

3. **Known Sites** [database icon]
   dbSNP, Mills
   Variant Recalibration

4. **ClinVar** [medical icon]
   4.1M Records
   Variant Classification

5. **AlphaMissense** [prediction icon]
   71M Predictions
   Pathogenicity Scoring

6. **RCSB PDB** [structure icon]
   Protein Structures
   Docking Targets

7. **Cryo-EM** [microscope icon]
   Structure Data
   Binding Pockets

Arrows flow rightward from each into the appropriate engine band.

---

### ━━━ BAND 2: GENOMIC FOUNDATION ENGINE (green tint band) ━━━

**Band header (left side, vertical green bar):** "Genomic Foundation Engine"
**Background tint:** Very light green (#E8F5E9)

**Header box centered at top of band:**
```
[DNA icon] HCLS AI Factory Pipeline
Three Integrated AI Engines
```

**Process flow (left to right, small rounded boxes connected by arrows):**

```
FASTQ → fq2bam/BWA-MEM2 → Sort → Mark Duplicates → BQSR → DeepVariant → VCF Output
         Parabricks                                  Parabricks  CNN Calling
```

**Metric badges along bottom of band:**
- [10-50x Speedup] [>99% Accuracy] [2-4 Hours] [11.7M Variants]

**Below process flow — Annotation sub-row:**
```
Variant Annotation          Embedding Generation
ClinVar/AlphaMissense →    BGE-small-en-v1.5 →    Milvus Ingest
VEP                        384 dimensions          3.56M Vectors
```

---

### ━━━ BAND 3: PRECISION INTELLIGENCE NETWORK (teal tint band — LARGEST BAND) ━━━

**Band header:** "Precision Intelligence Network"
**Background tint:** Very light teal (#E0F7FA)

**This band is EXPANDED from v1.0 to show all 11 agents.** It has two sub-layers:

#### Sub-layer 3A: RAG Architecture (top of band)

```
Retrieval Flow:
User Query → Embedding → Milvus Vector Search → Metadata Filter → Hybrid Search
                         [3.56M Variants]

Clinker Knowledge Layer:
Variant → Gene → Protein → Pathway → Disease → Drug
201 Genes | 13 Therapeutic Areas | 171 Druggable Targets
```

Claude LLM box (right side of sub-layer):
```
Claude Sonnet 4
Evidence Synthesis
Deterministic-Probabilistic Split
LLM explains, never computes
```

#### Sub-layer 3B: 11 Intelligence Agents (bottom of band — the new addition)

**Section header:** "11 Intelligence Agents — Cross-Agent Coordination"

**Four domain columns, matching the One-Pager layout:**

**Column 1: ONCOLOGY & IMMUNOTHERAPY** (green left accent)
```
Precision Oncology Agent
• Tumor profiling • Therapy matching
• 11 pediatric cancers
[cross-agent hub]

CAR-T Intelligence Agent
• Eligibility • CRS/ICANS
• Manufacturing optimization
[ELIANA: 82% CR]

Precision Biomarker Agent
• Biomarker discovery
• MRD detection • Biological age
```

**Column 2: SPECIALTY MEDICINE** (teal left accent)
```
Cardiology Intelligence Agent
• 6 calculators • 45 conditions
• 56 genes • Cardiotoxicity prevention
[dexrazoxane protocols]

Neurology Intelligence Agent
• 10 scales • 8 workflows
• Stroke to neuromuscular
[pediatric neurotoxicity]

Precision Autoimmune Agent
• 13 conditions • Flare prediction
• irAE profiling
```

**Column 3: DIAGNOSTICS & GENOMICS** (purple left accent)
```
Rare Disease Diagnostic Agent
• 88 diseases • 23 ACMG criteria
• Gene therapy eligibility
[cascade testing]

Pharmacogenomics Agent
• Drug-gene interactions
• TPMT/NUDT15 • HLA risk
[pediatric dosing]

Imaging Intelligence Agent
• Vista3D • MAISI • ViLaM3
• NIM-powered analysis

Single-Cell Intelligence Agent
• 57 cell types • TME profiling
• CAR-T target validation
```

**Column 4: CLINICAL OPERATIONS** (amber left accent)
```
Clinical Trial Intelligence Agent
• Patient-trial matching
• COG protocols • Pediatric MATCH
• Safety signal detection
```

**Dashed teal lines** connecting agents across columns (subtle):
- Oncology ↔ Cardiology (cardiotoxicity)
- Oncology ↔ CAR-T
- CAR-T ↔ Single-Cell (target validation)
- Rare Disease ↔ Imaging (screening)
- All agents → "/integrated-assessment" label

---

### ━━━ BAND 4: THERAPEUTIC DISCOVERY ENGINE (purple tint band) ━━━

**Band header:** "Therapeutic Discovery Engine"
**Background tint:** Very light purple (#F3E8FF)

**Process flow (left to right):**

```
Target Import → Structure Retrieval → Seed Molecule → MolMIM → Pediatric → DiffDock → Composite → Results
PDB / EMBL     Filtering by           CB-5083        BioNeMo   Safety       BioNeMo    Scoring     & Export
                resolution            SMILES                    Filter       NIM
                                                                [NEW in v2]
```

**NEW box (highlighted with emerald #059669 border) — Pediatric Safety Filter:**
```
Pediatric Safety Filter [NEW]
• BBB penetration (MW <500)
• hERG cardiac liability
• Hepatic immaturity
• Growth plate / teratogenicity
• 500 generated → ~89 pass → Top 3 ranked
[Pediatric Safe ✓]
```

**Metric badges:**
- [2 Minutes] [50-100 Candidates] [Complete Lineage] [Pediatric Filtered]

---

### ━━━ SHARED LAYER: MILVUS VECTOR DATABASE (between Bands 3-4, matching v1.0) ━━━

```
Milvus Vector Database
PyMilvus Client | IVF_FLAT Index | 384-dim Vectors
3.56M Variants Indexed | Hybrid Search
140+ Collections across 11 Agents
```

---

### ━━━ RIGHT COLUMN: OUTPUT CONSUMERS (spanning vertically alongside Bands 2-4) ━━━

**Updated from v1.0 — expanded to reflect new capabilities:**

1. **Researcher** [microscope icon]
   Output & Results
   Evidence-Cited Reports

2. **Clinician** [doctor icon]
   Clinical Decision Support
   Risk Calculators & Scales

3. **Tumor Board** [group icon] **[NEW]**
   30-Second Virtual MTB
   Multi-Agent Assessment

4. **Med Chemist** [molecule icon]
   Target Hypotheses
   Ranked Candidates

5. **Pediatric Oncologist** [child icon] **[NEW]**
   Cardiotoxicity Prevention
   PGx-Guided Dosing

6. **Reports** [document icon]
   Markdown • JSON • PDF
   FHIR R4 • Phenopacket

7. **Clinical Trials** [clipboard icon] **[NEW]**
   COG Protocol Matching
   Pediatric MATCH

Output arrows flow from engine bands rightward to these consumers. Labeled:
- SDF/Molecules (from Drug Discovery)
- FHIR/PDF (from Intelligence Network)
- VCF (from Genomics)

---

### ━━━ BOTTOM: PERFORMANCE METRICS BADGE (bottom-left, matching v1.0) ━━━

```
Performance Metrics Badge
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
<5 Hours, Not Months
Patient DNA → Drug Candidates
3 Engines | 11 Agents | 55+ Cross-Agent Functions
250K+ Lines | 5,000+ Tests | 24,043 Files
6 Pediatric Oncology Demos
Apache 2.0 — Free for Every Hospital
```

Green (#76B900) border, white background, small but dense text.

---

### ━━━ BOTTOM: NVIDIA INFRASTRUCTURE BAR (full width, matching v1.0) ━━━

**Green (#76B900) background, white text, NVIDIA logo on left.**

```
NVIDIA ▸ Built for NVIDIA DGX Spark

GPU Compute (GB10) | Grace (72 ARM cores) | 128 GB Unified Memory | NVMe Storage | ConnectX-7 | Parabricks | BioNeMo NIMs | MONAI | Docker | FastAPI | Streamlit UI
```

---

### ━━━ VERY BOTTOM: FOOTER LINE ━━━

```
Open Source v2.0 — Created by Adam Jones | Apache 2.0 | hcls-ai-factory.org | March 2026 | GitHub Available
```
Small, gray, centered.

---

## WHAT THIS v2.0 MUST COMMUNICATE vs v1.0

| v1.0 Said | v2.0 Says |
|-----------|-----------|
| 3 pipelines | 3 engines (renamed) |
| 5 agents (not shown individually) | 11 agents in 4 domain groups (all visible) |
| VCP/FTD demo use case | Pediatric oncology (5 patients, 6 demos) |
| Researcher + Clinician + Med Chemist outputs | + Tumor Board + Pediatric Oncologist + Clinical Trials |
| No cross-agent integration shown | Dashed teal lines between agents |
| No pediatric safety | Explicit Pediatric Safety Filter in Drug Discovery |
| "35,600 Clinker Matches / 3.5M Variants" | "3.56M Variants / 140+ Collections / 11 Agents" |
| "MVP" badge | "v2.0" badge with full production metrics |

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Same elegant architecture as v1.0** — recognizably the same platform, evolved
2. **Three engines renamed** — Genomic Foundation, Precision Intelligence, Therapeutic Discovery
3. **11 agents are now VISIBLE** — not hidden behind a generic "Knowledge Layer" box
4. **Pediatric oncology is the mission** — use case badge, patient names, demo count
5. **Cross-agent coordination exists** — dashed lines show agents talking to each other
6. **Pediatric safety filter is explicit** — new box in Drug Discovery with BBB, hERG, hepatic, growth plate
7. **More output consumers** — Tumor Board, Pediatric Oncologist, Clinical Trials added
8. **The numbers have grown** — 11 agents, 140+ collections, 55+ cross-agent functions, 250K+ lines
9. **Open source, free, Apache 2.0** — prominent in metrics badge and footer
10. **This runs on one DGX Spark** — same desktop machine, 10x more capability

The overall impression: v1.0 was a proof of concept. v2.0 is a production platform that could change pediatric oncology worldwide. Same architecture. Same machine. 10x the capability. Free for everyone.
