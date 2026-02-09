---
search:
  exclude: true
---

# Napkin AI Pro — HCLS AI Factory MVP on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the MVP proof build running on a single NVIDIA DGX Spark. This is what runs on the desk during the proof build — a complete genomics-to-drug-discovery pipeline.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of a polished enterprise technical white paper (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of an NVIDIA solutions architecture reference diagram.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a reference architecture poster a solutions architect pins to their wall during a proof build.

**Typography:**
- Title: Large, bold, sans-serif (Inter, Helvetica, or similar), deep navy (#1B2333)
- Subtitle: Smaller, medium gray (#666666), directly below title
- Section headers: Bold, dark navy (#1B2333) on white, with a thin NVIDIA green (#76B900) left-border accent or underline
- Sub-headers: Bold, teal (#1AAFCC)
- Body text: Small (8-10pt equivalent), clean sans-serif, dark gray (#333333)
- Metric callouts: Bold, slightly larger than body, inside small rounded green (#76B900) or teal (#1AAFCC) pill badges with white text

**Color Palette (exact):**
- NVIDIA Green: #76B900 — primary accent for all NVIDIA components, pipeline headers, infrastructure bar, metric badges
- Deep Navy: #1B2333 — title text, dark section bars, footer
- Teal: #1AAFCC — secondary accent for data flow lines, sub-headers, connection lines
- Light Gray: #F5F5F5 — card backgrounds, pipeline row backgrounds
- White: #FFFFFF — canvas, text on dark backgrounds, card interiors
- Orange: #F5A623 — drug candidate output badges, seed compound indicators
- Red: #DC2626 — pathogenic variant indicators, critical findings
- Purple: #7B2D8E — annotation database badges (ClinVar, AlphaMissense, VEP)
- Medium Gray: #666666 — metadata text, secondary labels
- Emerald Green: #059669 — knowledge base badges, drug-like status indicators

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons in a clean enterprise white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-modal triggers, bold NVIDIA green (#76B900) for primary pipeline flow
- Color-coded pipeline stages with distinct light background tints
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer
- Three distinct pipeline stage bands with clear left-to-right flow

---

## CANVAS STRUCTURE (Top to Bottom, 7 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "MVP Proof Build" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "HCLS AI Factory"
- **Subtitle line 1 (medium, gray #666666, centered):** "From Patient DNA to Novel Drug Candidates — MVP on NVIDIA DGX Spark"
- **Subtitle line 2 (smaller, gray #666666, centered):** "GB10 Grace Blackwell Superchip | 128 GB Unified Memory | < 5 Hours End-to-End"
- **Date / Author line (smallest, gray #999999, centered):** "February 2026 | Apache 2.0 | Author: Adam Jones"

**Right side — Key/Legend box** (small, top-right corner, thin gray border, white background):
```
Key
——————————————————
● Stage 1 — Genomics Pipeline (green)
● Stage 2 — RAG/Chat Pipeline (teal)
● Stage 3 — Drug Discovery Pipeline (orange)
→ Primary Data Flow (solid)
- → Cross-Modal Trigger (dashed)
◆ NVIDIA Component (green badge)
```
Use small colored dots: NVIDIA green for Stage 1, teal for Stage 2, orange for Stage 3.

---

### ━━━ BAND 2: PATIENT DATA SOURCE (left column, spanning vertically alongside Bands 3-5) ━━━

**Position:** Left edge of canvas, vertically stacked column of input cards. Each card has a thin-line icon, bold label, and 1-2 lines of detail. White background (#FFFFFF) with thin teal (#1AAFCC) left-border accent.

**Cards (top to bottom):**

1. **Patient DNA Sample** [DNA helix icon]
   Clinical specimen
   Extracted genomic DNA

2. **Illumina Sequencer** [sequencer icon]
   2×250 bp paired-end
   30× WGS coverage

3. **FASTQ Files** [file icon]
   ~200 GB per sample
   Raw sequencing reads

4. **GRCh38 Reference** [genome icon]
   3.1 GB reference genome
   Pre-indexed for BWA-MEM2

**Arrows:** Vertical downward flow between cards. From "FASTQ Files" → thick green arrow right into Band 3 (Stage 1). From "GRCh38 Reference" → thin teal arrow right into BWA-MEM2.

---

### ━━━ BAND 3: STAGE 1 — GENOMICS PIPELINE (first horizontal stage band) ━━━

**Background:** Very light green tint (#F0F9E0) with a green (#76B900) left-border accent bar.

**Header (left-aligned in the band):**
- Section label: "STAGE 1: GENOMICS PIPELINE" in bold navy (#1B2333)
- Subtitle: "NVIDIA Parabricks 4.6 — GPU-Accelerated" in teal (#1AAFCC)
- Time badge: rounded green pill with white text: "120-240 min"

**Process Cards (left to right flow within band):**

1. **NVIDIA Parabricks 4.6** [GPU chip icon, green badge]
   nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
   GPU-accelerated bioinformatics
   ◆ NVIDIA Component

2. **BWA-MEM2 Alignment** [alignment icon]
   fq2bam process
   Metric badges: "20-45 min" (green pill), "70-90% GPU" (teal pill)
   Output: Sorted BAM + BAI index

3. **Google DeepVariant** [neural network icon]
   CNN-based variant calling
   Metric badges: "10-35 min" (green pill), ">99% accuracy" (teal pill), "80-95% GPU" (teal pill)

4. **VCF Output** [file output icon]
   Metric badge: "~11.7M variants" (green pill)
   Sub-metrics: "3.5M QUAL>30 | 4.2M SNPs | 1.0M indels | 35K coding"

5. **Genomics Portal** [dashboard icon, small card]
   Flask :5000
   Pipeline progress, BAM stats

**Arrows:** Thick green directional arrows left-to-right between process cards. VCF Output has a thick green arrow down into Band 4 (Stage 2).

---

### ━━━ BAND 4: STAGE 2 — RAG/CHAT PIPELINE (second horizontal stage band) ━━━

**Background:** Very light teal tint (#E8F8FC) with a teal (#1AAFCC) left-border accent bar.

**Header:**
- Section label: "STAGE 2: RAG/CHAT PIPELINE" in bold navy (#1B2333)
- Subtitle: "Milvus + Claude — RAG-Grounded Target Identification" in teal (#1AAFCC)
- Time badge: "Interactive" (teal pill)

**This band is wider and has three sub-rows:**

**Sub-row 1: Variant Annotation (left-to-right)**

1. **ClinVar** [database icon, purple badge]
   NCBI clinical variants
   Metric badge: "4.1M variants" (purple pill)
   ~35,616 patient matches

2. **AlphaMissense** [AI brain icon, purple badge]
   DeepMind pathogenicity predictions
   Metric badge: "71M predictions" (purple pill)
   Thresholds: >0.564 pathogenic | 0.34-0.564 ambiguous | <0.34 benign

3. **Ensembl VEP** [function icon, purple badge]
   Functional consequences
   Impact: HIGH | MODERATE | LOW | MODIFIER

**Sub-row 2: Vector Database (left-to-right)**

4. **BGE-small-en-v1.5** [embedding icon]
   384-dim text embeddings
   Variant text summaries

5. **Milvus 2.4** [vector database icon, teal badge]
   Port 19530 | IVF_FLAT | nlist=1024
   Metric badge: "3.5M embeddings" (teal pill)
   17-field schema | COSINE metric

6. **Knowledge Base** [book icon, emerald badge]
   Metric badges: "201 genes" (emerald pill), "13 areas" (emerald pill), "171 druggable" (emerald pill)
   "85% druggability rate"

**Sub-row 3: LLM Reasoning (left-to-right)**

7. **Anthropic Claude** [LLM icon, navy badge]
   claude-sonnet-4-20250514 | temp=0.3
   RAG-grounded clinical reasoning
   10 therapeutic area query expansion maps

8. **Streamlit Chat** [chat icon, light gray card]
   Port 8501
   Interactive natural language interface

9. **Target Hypothesis** [target icon, orange badge]
   Gene + UniProt ID + evidence chain
   Confidence assessment + therapeutic area
   → Feeds into Stage 3

**Arrows:** Sub-row 1 feeds into sub-row 2 (annotation → embedding/indexing). Sub-row 2 feeds into sub-row 3 (embeddings → Claude RAG). Knowledge Base also feeds into Claude. Target Hypothesis has thick green arrow down into Band 5 (Stage 3).

---

### ━━━ BAND 5: STAGE 3 — DRUG DISCOVERY PIPELINE (third horizontal stage band) ━━━

**Background:** Very light orange tint (#FFF8F0) with an orange (#F5A623) left-border accent bar.

**Header:**
- Section label: "STAGE 3: DRUG DISCOVERY PIPELINE" in bold navy (#1B2333)
- Subtitle: "BioNeMo NIM + RDKit — AI-Driven Drug Discovery" in teal (#1AAFCC)
- Time badge: "8-16 min" (orange pill)

**Process Cards (left to right, 10-stage pipeline condensed):**

1. **RCSB PDB** [protein structure icon, purple badge]
   Cryo-EM / X-ray structures
   VCP demo: 8OOI, 9DIL, 7K56, 5FTK

2. **Seed Compound** [molecule icon, orange badge]
   CB-5083 (VCP inhibitor)
   Phase I clinical candidate
   SMILES input

3. **MolMIM NIM** [generative AI icon, green badge]
   Port 8001 | POST /v1/generate
   Molecule generation
   Metric badge: "100 candidates" (green pill)
   ◆ NVIDIA BioNeMo

4. **RDKit QC** [chemistry icon, teal badge]
   Lipinski: MW≤500, LogP≤5, HBD≤5, HBA≤10
   QED > 0.67 (drug-like)
   TPSA < 140 Å²

5. **DiffDock NIM** [docking icon, green badge]
   Port 8002 | POST /v1/dock
   Molecular docking
   Metric badge: "Diffusion-based" (green pill)
   ◆ NVIDIA BioNeMo

6. **Composite Ranking** [ranking icon, teal badge]
   30% generation + 40% docking + 30% QED
   dock_norm = max(0, min(1, (10+dock)/20))

7. **100 Drug Candidates** [trophy/results icon, orange badge]
   Ranked by composite score
   Metric badges: "87 Lipinski PASS" (emerald pill), "72 QED>0.67" (emerald pill)
   Top: -11.4 kcal/mol | 0.89 composite

8. **PDF Report** [report icon, light gray card]
   ReportLab
   Full provenance

9. **Discovery UI** [dashboard icon, light gray card]
   Streamlit :8505
   Interactive results explorer

**Arrows:** Green directional arrows left-to-right. RCSB PDB and Seed Compound both feed into MolMIM. MolMIM → RDKit QC → DiffDock → Composite Ranking → 100 Drug Candidates → PDF Report. Candidates also arrow to Discovery UI.

---

### ━━━ BAND 6: ORCHESTRATION + MONITORING (bottom infrastructure band) ━━━

**Background:** Light gray (#F5F5F5) with a thin green (#76B900) top-border accent.

**Left section — Nextflow DSL2:**
- Card: "HLS-Pipeline v1.0.0" in bold navy
- Sub-cards: genomics.nf | rag_chat.nf | drug_discovery.nf | reporting.nf
- Modes badge: "full | target | drug | demo | genomics_only"
- Dashed teal lines connecting up to each stage band

**Center section — Landing Page:**
- Card: "HCLS AI Factory Landing Page" in bold navy
- "Flask :8080 | 10-Service Health Monitor"
- Small grid of 10 green dots (all healthy)

**Right section — Monitoring Stack:**
- Card: "Observability" in bold navy
- Sub-cards in a row:
  - "Grafana :3000" (green badge) — Dashboards
  - "Prometheus :9099" (teal badge) — Metrics TSDB
  - "Node Exporter :9100" — CPU/RAM/Disk
  - "DCGM Exporter :9400" — GPU Metrics
- Arrows: Node + DCGM → Prometheus → Grafana

---

### ━━━ BAND 7: HARDWARE FOOTER (bottom of canvas) ━━━

**Background:** Dark navy (#1B2333) full-width bar.

**Left:** NVIDIA logo mark + "NVIDIA DGX Spark" in white

**Center (metric badges, green pills on navy background):**
- "GB10 GPU" | "128 GB LPDDR5x" | "ARM64 Cores" | "NVMe Storage" | "NVLink-C2C"

**Right:**
- "$3,999" in large bold white text
- "Desktop-Class AI Factory" in smaller gray (#999999) text below

---

### ━━━ CROSS-MODAL OVERLAY ━━━

**Position:** Right edge of canvas, spanning vertically alongside Bands 3-5, mirror of Band 2.

**Title:** "Cross-Modal Integration" in bold navy with teal left-border

**Cards (top to bottom):**

1. **Imaging Intelligence Agent** [medical imaging icon, navy badge]
   CT / MRI / X-Ray AI
   MONAI Deploy on DGX Spark

2. **Lung-RADS 4B+ Trigger** [alert icon, red badge]
   FHIR ServiceRequest
   → Genomics analysis
   Dashed red arrow left into Stage 1

3. **NVIDIA FLARE** [federated icon, green badge]
   Federated Learning
   Phase 3: Multi-site deployment
   Data stays local

**Arrows:** Dashed teal arrows from Imaging Agent → Stage 1 (trigger genomics). Dashed teal from Stage 3 → Imaging Agent (combined report).

---

### ━━━ DEPLOYMENT ROADMAP (bottom-right corner, small inset) ━━━

**Position:** Small card in the bottom-right corner, above the footer bar.

**Title:** "Deployment Roadmap" in bold navy, small font

**Three-row mini-table:**
| Phase | Hardware | Scale |
|---|---|---|
| 1 — Proof Build | DGX Spark | Docker Compose, $3,999 |
| 2 — Departmental | DGX B200 | Kubernetes, $500K-$1M |
| 3 — Enterprise | DGX SuperPOD | FLARE Federated, $7M+ |

Green left-border on Phase 1 (current), teal on Phase 2, navy on Phase 3.

---

## VCP/FTD DEMO CALLOUT (floating annotation card)

**Position:** Centered overlay or side callout, connected by dashed lines to relevant components.

**Card design:** White background, thin green (#76B900) border, slight drop shadow.

**Title:** "VCP/FTD Demo — Live Walkthrough" in bold teal

**Content (compact, dense):**
```
Target: VCP (p97) — UniProt P55072
Disease: Frontotemporal Dementia, ALS, IBMPFD
Variant: rs188935092 (chr9:35065263 G>A)
ClinVar: Pathogenic | AlphaMissense: 0.87
Seed: CB-5083 (Phase I VCP inhibitor)
Structures: 8OOI | 9DIL | 7K56 | 5FTK
Binding Site: D2 ATPase (~450 Å³, druggability 0.92)
Result: 100 candidates | Top composite: 0.89
        +39% improvement over CB-5083 seed
```

Dashed teal lines from this callout to: VCF Output (Stage 1), Claude (Stage 2), MolMIM (Stage 3), and 100 Drug Candidates (Stage 3).

---

## FINAL NOTES FOR THE AI GENERATOR

1. **This is ONE infographic on ONE canvas.** Not a slide deck. Not separate pages. Everything on one dense poster.
2. **Left-to-right primary flow** within each stage band. **Top-to-bottom** flow between stages.
3. **Every metric badge mentioned above must appear.** These are the key numbers audiences remember.
4. **Color consistency is critical.** Green = NVIDIA/GPU. Teal = data flow/secondary. Orange = drug discovery outputs. Purple = annotation databases. Red = pathogenic/critical. Navy = structural/titles. Emerald = knowledge base/healthy.
5. **The three stage bands are the visual centerpiece.** They should be the largest elements, clearly showing the pipeline flow.
6. **The VCP/FTD demo callout** should be prominent — it's what the audience sees live.
7. **All port numbers must be visible** on their respective components.
8. **The $3,999 price point** in the footer should be visually prominent — it's a key selling point.
9. **No proprietary tools.** Everything shown is open or NVIDIA-licensed.
10. **The cross-modal section** on the right shows the connection to the Imaging Intelligence Agent — this is the broader HCLS AI Factory ecosystem.
