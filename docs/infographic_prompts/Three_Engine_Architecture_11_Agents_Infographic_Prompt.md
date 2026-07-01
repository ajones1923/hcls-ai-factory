# Nano Banana Pro — HCLS AI Factory: Three-Engine Architecture & 11 Intelligence Agents

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck, NOT multiple pages. Every element described below appears on ONE canvas. This diagram is the definitive visual representation of the HCLS AI Factory platform — showing the three GPU-accelerated engines and all 11 intelligence agents that together transform patient DNA into drug candidates in under 5 hours. It should be the single image that makes someone understand the entire platform architecture, scope, and clinical capability in 30 seconds of looking at it.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a reference architecture poster a solutions architect pins to their wall — comprehensive enough to be the single source of truth for the platform, clean enough to present to executives.

**Typography:**
- Title: Large, bold, sans-serif (Inter, Helvetica, or similar), deep navy (#1B2333)
- Subtitle: Smaller, medium gray (#666666), directly below title
- Section headers: Bold, dark navy (#1B2333) on white, with a thin NVIDIA green (#76B900) left-border accent or underline
- Sub-headers (engine names): Bold, teal (#1AAFCC), slightly larger than body
- Agent names: Bold, dark navy (#1B2333), 9-10pt equivalent
- Agent descriptions: Regular, dark gray (#333333), 8pt equivalent
- Domain group labels: Bold, teal (#1AAFCC), uppercase, small, letter-spaced
- Body text: Small (8-10pt equivalent), clean sans-serif, dark gray (#333333)
- Metric callouts: Bold, slightly larger than body, inside small rounded green (#76B900) or teal (#1AAFCC) pill badges with white text
- Footer text: Small, gray (#999999)

**Color Palette (exact — use these hex values):**
- NVIDIA Green: #76B900 — primary accent for engine headers, pipeline arrows, metric badges, infrastructure bar
- Deep Navy: #1B2333 — title text, agent names, section bars, footer background
- Teal: #1AAFCC — engine sub-headers, domain group labels, cross-agent connection lines, secondary badges
- Light Gray: #F5F5F5 — card backgrounds, engine section backgrounds
- White: #FFFFFF — canvas, text on dark backgrounds, card interiors
- Amber/Orange: #F5A623 — warning/clinical alert indicators
- Red: #DC2626 — critical finding indicators
- Purple: #7B2D8E — rare disease / genetic domain accent
- Medium Gray: #666666 — metadata text, secondary labels, descriptions
- Emerald Green: #059669 — positive/verified indicators

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box and agent card
- Thin-line icons (16x16 to 24x24) next to engines and agents — simple, monochrome line icons matching the the platform white paper style (not emoji, not filled icons)
- Directional arrows: bold NVIDIA green (#76B900) for primary engine-to-engine data flow, dashed teal (#1AAFCC) for cross-agent connections, solid medium gray (#999999) for internal data paths
- Color-coded domain groups with distinct subtle left-border accents on agent cards
- Metric badges: small rounded pills (16px height, auto-width) with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer
- "Apache 2.0" badge prominently visible in title bar or footer

---

## CANVAS STRUCTURE (Top to Bottom, 6 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "v1.0.0" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "HCLS AI Factory"
- **Subtitle line 1 (medium, teal #1AAFCC, centered):** "Patient DNA → Drug Candidates in <5 Hours"
- **Subtitle line 2 (smaller, gray #666666, centered):** "Three-Engine Architecture | 11 Intelligence Agents | Open Source (Apache 2.0)"
- **Subtitle line 3 (smallest, gray #999999, centered):** "Designed with Pediatric Oncology as the Primary Use Case — Built for NVIDIA DGX Spark"

**Right side — Key/Legend box** (small, top-right corner, thin gray border, white background):
```
Key
━━━━━━━━━━━━━━━━━━━━━
→   Primary Data Flow (green)
- → Cross-Agent Coordination (teal dashed)
█   Engine Block
▪   Agent Card
●   Metric Badge
```

---

### ━━━ BAND 2: THREE-ENGINE PIPELINE (the architectural backbone — ~25% of canvas height) ━━━

**Position:** Full-width horizontal band below title. This is the core architectural story: three engines connected left-to-right showing the complete DNA-to-drug pipeline. Each engine is a large rounded-corner container (8px radius) with a light gray (#F5F5F5) background, thin NVIDIA green (#76B900) top-border accent (3px), and structured internal content.

**Layout:** Three engine blocks of equal width arranged horizontally, separated by bold green directional arrows (→) showing data flow left-to-right.

#### ENGINE BLOCK 1 (left) — Genomic Foundation Engine

**Header bar (inside the block, top):**
- Icon: DNA helix (thin-line, 24x24, navy)
- Text: "Genomic Foundation Engine" (bold, teal #1AAFCC, 14pt equivalent)
- Small badge right-aligned: "STAGE 1" (navy pill, white text)

**Content — 4 process cards arranged horizontally within the block:**

Card 1:
- Label: "FASTQ Input" (bold, navy)
- Detail: "Raw sequencing reads" (gray, 8pt)
- Detail: "Tumor + Germline DNA" (gray, 8pt)
- Icon: file/document (thin-line)

→ (internal gray arrow)

Card 2:
- Label: "BWA-MEM2 Alignment" (bold, navy)
- Detail: "GPU-accelerated via Parabricks" (gray, 8pt)
- Metric badge: "45 min" (green pill) — below: "vs 12 hrs CPU" (gray strikethrough text)
- Icon: align/layers (thin-line)

→ (internal gray arrow)

Card 3:
- Label: "DeepVariant Calling" (bold, navy)
- Detail: "CNN-based variant classification" (gray, 8pt)
- Metric badge: ">99.7% accuracy" (green pill)
- Icon: neural network (thin-line)

→ (internal gray arrow)

Card 4:
- Label: "Annotated VCF" (bold, navy)
- Detail: "ClinVar + AlphaMissense" (gray, 8pt)
- Metric badge: "11.7M variants" (green pill)
- Icon: document with checkmark (thin-line)

**Footer line (bottom of block, small gray text):**
"Powered by NVIDIA Parabricks on DGX Spark | GRCh38 Reference Genome"

---

#### BOLD GREEN ARROW (→) connecting Engine 1 to Engine 2

Large, prominent, NVIDIA green (#76B900) arrow with label above:
**"VCF + Annotations"** (small, green text)

---

#### ENGINE BLOCK 2 (center, slightly larger than others — this is the intelligence hub) — Precision Intelligence Network

**Header bar:**
- Icon: Brain/network (thin-line, 24x24, navy)
- Text: "Precision Intelligence Network" (bold, teal #1AAFCC, 14pt equivalent)
- Small badge right-aligned: "STAGE 2" (navy pill, white text)

**Content — 3 core component cards arranged horizontally:**

Card 1:
- Label: "Milvus Vector Database" (bold, navy)
- Detail: "3.56M annotated variant embeddings" (gray, 8pt)
- Detail: "BGE-small-en-v1.5 | 384 dimensions" (gray, 8pt)
- Detail: "IVF_FLAT / COSINE indexing" (gray, 8pt)
- Metric badge: "3.56M vectors" (teal pill)
- Icon: database cylinder (thin-line)

Card 2:
- Label: "Claude LLM Synthesis" (bold, navy)
- Detail: "Evidence-grounded reasoning" (gray, 8pt)
- Detail: "Deterministic-probabilistic split" (gray, 8pt)
- Detail: "LLM explains, never computes scores" (gray, 8pt)
- Metric badge: "RAG Architecture" (teal pill)
- Icon: brain/AI (thin-line)

Card 3:
- Label: "11 Intelligence Agents" (bold, navy)
- Detail: "Cross-agent coordination via HTTP" (gray, 8pt)
- Detail: "/integrated-assessment endpoints" (gray, 8pt)
- Detail: "55+ cross-agent functions" (gray, 8pt)
- Metric badge: "11 agents" (green pill)
- Icon: grid of nodes (thin-line)

**Below the 3 cards — a thin dashed teal line points DOWN to Band 3** with label:
**"Agent Detail Below ↓"** (small, teal)

**Footer line:**
"Retrieval-Augmented Generation | Evidence-Cited Responses | All Clinical Calculations Deterministic"

---

#### BOLD GREEN ARROW (→) connecting Engine 2 to Engine 3

Large, prominent, NVIDIA green (#76B900) arrow with label above:
**"Actionable Targets"** (small, green text)

---

#### ENGINE BLOCK 3 (right) — Therapeutic Discovery Engine

**Header bar:**
- Icon: Molecule/flask (thin-line, 24x24, navy)
- Text: "Therapeutic Discovery Engine" (bold, teal #1AAFCC, 14pt equivalent)
- Small badge right-aligned: "STAGE 3" (navy pill, white text)

**Content — 4 process cards arranged horizontally:**

Card 1:
- Label: "MolMIM Generation" (bold, navy)
- Detail: "NVIDIA BioNeMo" (gray, 8pt)
- Detail: "500 candidates per target" (gray, 8pt)
- Metric badge: "500 molecules" (green pill)
- Icon: molecular generation (thin-line)

→ (internal gray arrow)

Card 2:
- Label: "DiffDock Scoring" (bold, navy)
- Detail: "Protein-ligand docking" (gray, 8pt)
- Detail: "Binding affinity prediction" (gray, 8pt)
- Metric badge: "kcal/mol ranked" (green pill)
- Icon: protein-ligand (thin-line)

→ (internal gray arrow)

Card 3:
- Label: "RDKit ADMET" (bold, navy)
- Detail: "Drug-likeness scoring" (gray, 8pt)
- Detail: "Chemical analysis" (gray, 8pt)
- Icon: chemistry (thin-line)

→ (internal gray arrow)

Card 4:
- Label: "Pediatric Safety Filter" (bold, navy)
- Detail: "BBB penetration" (gray, 8pt)
- Detail: "hERG cardiac liability" (gray, 8pt)
- Detail: "Hepatic immaturity" (gray, 8pt)
- Detail: "Growth plate / teratogenicity" (gray, 8pt)
- Metric badge: "Pediatric Safe" (emerald #059669 pill)
- Icon: shield with checkmark (thin-line)

**Footer line:**
"NVIDIA BioNeMo + RDKit | Novel Drug Candidates with Pediatric Safety Assessment"

---

### ━━━ BAND 3: 11 INTELLIGENCE AGENTS (the detail layer — ~40% of canvas height) ━━━

**Position:** Full-width band below the three engines. This is the agent constellation. It connects visually to Engine Block 2 (Precision Intelligence Network) via the dashed teal line.

**Section header bar (full width, navy #1B2333 background, white text):**
"11 Intelligence Agents — Precision Intelligence Network"

**Layout:** Four domain groups arranged horizontally. Each group has a teal (#1AAFCC) uppercase domain label at top, then agent cards stacked vertically below.

**Agent card design (consistent for all 11):**
- Rounded-corner rectangle (8px radius)
- Light gray (#F5F5F5) background
- Thin left-border accent in the domain's color
- Agent name: bold, navy (#1B2333), 9pt
- Description: regular, dark gray (#333333), 8pt
- 1-2 metric badges where applicable (green or teal pills)
- Thin-line icon (16x16) left of agent name

---

#### DOMAIN GROUP 1 (leftmost column) — Left border color: #76B900 (NVIDIA Green)

**Domain label (teal, uppercase, letter-spaced):** "ONCOLOGY & IMMUNOTHERAPY"

**Agent Card 1.1:**
- Icon: target/crosshair (thin-line)
- Name: "Precision Oncology Agent"
- Description: "Molecular tumor profiling, targeted therapy matching, resistance mechanism prediction, clinical trial correlation."
- Badge: "11 pediatric cancer types" (green pill)

**Agent Card 1.2:**
- Icon: T-cell (thin-line)
- Name: "CAR-T Intelligence Agent"
- Description: "CAR-T candidate evaluation, manufacturing protocol optimization, CRS/ICANS toxicity prediction, response biomarker tracking."
- Badge: "ELIANA: 82% CR" (teal pill)

**Agent Card 1.3:**
- Icon: test tube (thin-line)
- Name: "Precision Biomarker Agent"
- Description: "Biomarker discovery and analysis, biological age modeling, disease trajectory prediction, pharmacogenomic mapping."
- Badge: "MRD detection" (teal pill)

---

#### DOMAIN GROUP 2 (center-left column) — Left border color: #1AAFCC (Teal)

**Domain label:** "SPECIALTY MEDICINE"

**Agent Card 2.1:**
- Icon: heart with pulse (thin-line)
- Name: "Cardiology Intelligence Agent"
- Description: "6 validated risk calculators (ASCVD, MAGGIC, EuroSCORE II, CHA₂DS₂-VASc, HAS-BLED, HEART), 11 clinical workflows, 45 conditions, 56 cardiac genes, 51 guideline references."
- Badges: "6 calculators" (green pill) + "45 conditions" (teal pill)

**Agent Card 2.2:**
- Icon: brain (thin-line)
- Name: "Neurology Intelligence Agent"
- Description: "10 validated clinical scales (NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn & Yahr), 8 workflows spanning stroke, dementia, epilepsy, MS, Parkinson's, brain tumors, headache, and neuromuscular disease."
- Badges: "10 scales" (green pill) + "8 workflows" (teal pill)

**Agent Card 2.3:**
- Icon: shield/immune (thin-line)
- Name: "Precision Autoimmune Agent"
- Description: "Autoimmune disease classification across 13 conditions, treatment pathway analysis, flare prediction, immunological profiling."
- Badge: "13 conditions" (teal pill)

---

#### DOMAIN GROUP 3 (center-right column) — Left border color: #7B2D8E (Purple)

**Domain label:** "DIAGNOSTICS & GENOMICS"

**Agent Card 3.1:**
- Icon: DNA magnifying glass (thin-line)
- Name: "Rare Disease Diagnostic Agent"
- Description: "Phenotype-driven differential diagnosis across 88 rare diseases, ACMG/AMP variant classification with 23 of the 28 criteria, gene therapy eligibility assessment, diagnostic algorithm recommendation."
- Badges: "88 diseases" (purple #7B2D8E pill, white text) + "23 ACMG criteria" (teal pill)

**Agent Card 3.2:**
- Icon: pill/DNA (thin-line)
- Name: "Pharmacogenomics Intelligence Agent"
- Description: "Drug-gene interaction analysis, HLA allele risk assessment, dosing optimization, adverse reaction prediction."
- Badge: "TPMT/NUDT15/CYP" (teal pill)

**Agent Card 3.3:**
- Icon: scan/x-ray (thin-line)
- Name: "Imaging Intelligence Agent"
- Description: "NVIDIA NIM-powered medical image analysis using Vista3D, MAISI, and ViLaM3 for segmentation, synthetic generation, and visual question answering."
- Badge: "3 NVIDIA NIMs" (green pill)

**Agent Card 3.4:**
- Icon: cell cluster (thin-line)
- Name: "Single-Cell Intelligence Agent"
- Description: "Cell type annotation across 57 cell types, tumor microenvironment profiling, spatial transcriptomics, ligand-receptor interaction mapping, subclonal architecture analysis, CAR-T target validation."
- Badges: "57 cell types" (teal pill) + "TME profiling" (green pill)

---

#### DOMAIN GROUP 4 (rightmost column) — Left border color: #F5A623 (Amber)

**Domain label:** "CLINICAL OPERATIONS"

**Agent Card 4.1:**
- Icon: clipboard/checklist (thin-line)
- Name: "Clinical Trial Intelligence Agent"
- Description: "Protocol design optimization, automated patient-trial matching, site selection scoring, enrollment prediction, safety signal detection, competitive landscape analysis."
- Badges: "COG trials" (amber #F5A623 pill, white text) + "Pediatric MATCH" (teal pill)

---

#### CROSS-AGENT CONNECTION LINES (overlaid on Band 3)

Thin dashed teal (#1AAFCC) lines connecting agents that coordinate with each other. These should be subtle — not overwhelming — but show the web of integration:

- Oncology ↔ CAR-T (within Oncology & Immunotherapy group)
- Oncology ↔ Cardiology (cross-group: cardiotoxicity)
- Oncology ↔ Neurology (cross-group: neurotoxicity)
- Oncology ↔ PGx (cross-group: drug-gene interactions)
- Oncology ↔ Clinical Trial (cross-group: trial matching)
- CAR-T ↔ Single-Cell (cross-group: target validation)
- CAR-T ↔ Biomarker (cross-group: CD19/CD22 expression)
- Rare Disease ↔ Imaging (cross-group: trilateral screening)

Label one representative line: "Cross-Agent /integrated-assessment" (small teal text)

---

### ━━━ BAND 4: METRICS BAR (narrow, high-impact) ━━━

**Position:** Full-width horizontal strip. NVIDIA green (#76B900) background. White bold text. This is the credibility bar — the numbers that make the platform undeniable.

**Seven metrics evenly spaced across the bar, centered:**

| <5 hrs | 250K+ | 3.56M | 11 | 60+ | 5,000+ | 24,043 |
|---|---|---|---|---|---|---|
| DNA to Drug | Lines of Code | Searchable Variants | AI Agents | Clinical Workflows | Automated Tests | Files, Zero Failures |

Each metric rendered as:
- Large number/value (white, bold, 18pt equivalent)
- Small label below (white, regular, 8pt)
- Thin white vertical divider lines between each metric

---

### ━━━ BAND 5: TECHNOLOGY STACK BAR (narrow, informational) ━━━

**Position:** Full-width horizontal strip below metrics. Light gray (#F5F5F5) background. Small text showing the complete technology stack.

**Two rows of technology labels, organized by category:**

Row 1 (left to right):
- **Compute:** "NVIDIA DGX Spark | GB10 Grace Blackwell | 128 GB Unified Memory"
- **Genomics:** "Parabricks 4.6 (DeepVariant + BWA-MEM2)"
- **AI/LLM:** "Claude (Anthropic) | BioNeMo NIMs"

Row 2 (left to right):
- **Database:** "Milvus (vectors) | ClinVar | AlphaMissense"
- **Chemistry:** "RDKit | MolMIM | DiffDock"
- **Infrastructure:** "Docker | FastAPI | Streamlit | Prometheus | Grafana"

Each category label in teal (#1AAFCC, bold), values in dark gray (#333333, regular), separated by pipes.

---

### ━━━ BAND 6: FOOTER BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. Navy (#1B2333) background. White and gray text.

**Left:** NVIDIA logo mark (small, white) + "Built for NVIDIA DGX Spark"

**Center:** "Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org"

**Right:** "March 2026 | Pediatric Oncology as Primary Use Case"

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges scattered throughout (small rounded green pills with white text):**
- "45 min alignment" on BWA-MEM2 card
- ">99.7% accuracy" on DeepVariant card
- "11.7M variants" on VCF output card
- "3.56M vectors" on Milvus card
- "11 agents" on agent component card
- "500 molecules" on MolMIM card
- "Pediatric Safe" on safety filter card
- Domain-specific badges on agent cards (as specified above)

**Deterministic-probabilistic annotation** (small dashed box, bottom-right of Engine Block 2):
- "Clinical Safety Architecture"
- "Calculators compute scores (deterministic)"
- "LLM synthesizes explanations (probabilistic)"
- "No clinical value ever generated by LLM"

**Apache 2.0 annotation** (small badge, prominent, near title):
- Rounded green (#76B900) badge with white text: "Apache 2.0 — Free for Every Hospital on Earth"

**Data flow arrow style guide:**
- **Bold green (#76B900) arrows:** Engine-to-engine primary data flow (FASTQ → VCF → Targets → Candidates). These are the largest, most prominent arrows on the canvas.
- **Solid gray (#999999) arrows:** Internal process flow within engine blocks.
- **Dashed teal (#1AAFCC) lines:** Cross-agent connections in Band 3 and the link from Engine 2 down to the agent detail.

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Three engines, one pipeline** — Data flows left to right: DNA goes in, drug candidates come out. The three engines are visually distinct, equally prominent, and connected by bold arrows.
2. **11 agents organized by clinical domain** — Four domain groups (Oncology & Immunotherapy, Specialty Medicine, Diagnostics & Genomics, Clinical Operations) with every agent visible and described.
3. **The agents are the intelligence layer** — They live inside the Precision Intelligence Network (Engine 2) and are the bridge between genomic data and therapeutic discovery.
4. **Cross-agent coordination is real** — Dashed teal lines show agents communicating across domains. This isn't 11 isolated tools — it's a coordinated network.
5. **Every number is specific** — 11.7M variants, 3.56M vectors, 88 rare diseases, 57 cell types, 10 clinical scales, 6 risk calculators, 23 ACMG criteria. These aren't round numbers. They're verified counts.
6. **Pediatric safety is built in** — The Therapeutic Discovery Engine explicitly shows pediatric safety filters (BBB, cardiac, hepatic, growth plate). This isn't an afterthought.
7. **The technology stack is identified** — NVIDIA Parabricks, BioNeMo, DGX Spark, Milvus, Claude, RDKit — everything is named and placed.
8. **It's open source** — Apache 2.0 badge is prominent. "Free for Every Hospital on Earth" is visible.
9. **One person built this** — "Created by Adam Jones" in the footer. The scope-to-individual ratio is extraordinary.
10. **It runs on one machine** — "NVIDIA DGX Spark" is mentioned in the title, the engine blocks, and the infrastructure bar. This is desktop-scale, not data-center-scale.

The overall impression should be: a complete, technically precise, architecturally elegant precision medicine platform — dense enough to be a reference poster, clean enough to present to Jensen Huang, and clinically grounded enough to hang on the wall at St. Jude Children's Research Hospital.
