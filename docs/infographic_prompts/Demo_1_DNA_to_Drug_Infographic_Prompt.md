# Napkin AI Pro — DNA to Drug: End-to-End Precision Medicine Pipeline

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the complete HCLS AI Factory pipeline: Patient DNA to Drug Candidates in under 5 hours, running on a single NVIDIA DGX Spark. Three engines — Genomic Foundation, Precision Intelligence, and Therapeutic Discovery — connected end-to-end for a single pediatric oncology patient.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

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
- Amber/Orange: #F5A623 — clinical finding badges, elevated indicators
- Red: #DC2626 — critical finding indicators, bold alert arrows, IKZF1plus reclassification
- Emerald Green: #059669 — safety-pass indicators, normal status
- Medium Gray: #666666 — metadata text, secondary labels

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons matching the the platform white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-agent communication, bold red (#DC2626) for critical finding paths
- Color-coded pipeline rows with distinct light background tints
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer

---

## CANVAS STRUCTURE (Top to Bottom, 6 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "Demo 1 of 6" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "DNA to Drug — End-to-End Precision Medicine Pipeline"
- **Subtitle line 1 (medium, gray #666666, centered):** "Patient DNA to Drug Candidates in Under 5 Hours"
- **Subtitle line 2 (smaller, gray #666666, centered):** "NVIDIA DGX Spark | GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "March 2026 | HCLS AI Factory"

**Right side — Patient Card** (prominent, white background, thin red #DC2626 left-border accent, rounded corners):
```
PATIENT
——————————————————
Name: Evelyn (de-identified)
Age: 8 years old
Dx: B-cell Acute Lymphoblastic Leukemia (B-ALL)
Sample: Whole genome sequencing (FASTQ)
Goal: Identify actionable variants → drug candidates
```
Small thin-line icon of a child silhouette next to the name. Red (#DC2626) small pill badge: "Pediatric Oncology"

---

### ━━━ BAND 2: GENOMIC FOUNDATION ENGINE (full width, largest band) ━━━

**Position:** Full-width horizontal band below the title bar. Light green background tint (#E8F5E9). This is the first of three engines.

**Left label (vertical text, NVIDIA green #76B900 background, white bold text):** "ENGINE 1 — Genomic Foundation"

**Contents (left to right, connected by solid gray arrows):**

**Card 1: FASTQ Input** [DNA helix icon]
- Thin green-bordered (#76B900) rounded rectangle, white interior
- Header: "FASTQ Input" in bold navy
- Body: Whole genome sequencing | Paired-end reads | ~120 GB raw data
- Small teal pill badge: "Evelyn's DNA"

**Arrow (solid gray #999999, rightward):** labeled "Raw Reads"

**Card 2: BWA-MEM2 Alignment** [align icon]
- Thin green-bordered rounded rectangle
- Header: "BWA-MEM2 Alignment" in bold navy
- Body: GPU-accelerated read mapping | Parabricks 4.6 | Reference genome GRCh38
- Metric badge (green pill): "< 45 min alignment"

**Arrow (solid gray, rightward):** labeled "Aligned BAM"

**Card 3: DeepVariant Calling** [magnifying glass icon]
- Thin green-bordered rounded rectangle
- Header: "DeepVariant Variant Calling" in bold navy
- Body: Deep learning variant caller | SNVs + Indels | GPU-accelerated inference
- Metric badge (green pill): "> 99.7% accuracy"

**Arrow (solid gray, rightward):** labeled "Called Variants"

**Card 4: VCF Output** [document icon]
- Thin green-bordered rounded rectangle
- Header: "VCF Output" in bold navy
- Body: Variant Call Format | Annotated with ClinVar (4.1M records) | AlphaMissense (71M predictions)
- Metric badge (green pill): "11.7M variants"

**Arrow (solid gray, rightward):** labeled "Annotated Variants"

**Card 5: Milvus Vector Store** [database icon]
- Thin teal-bordered (#1AAFCC) rounded rectangle
- Header: "Milvus Vector Database" in bold teal
- Body: BGE-small-en-v1.5 embeddings | 384-dim vectors | Indexed for semantic search
- Metric badge (teal pill): "3.56M searchable vectors"

**Below the pipeline row, a thin dashed line annotation:** "Parabricks 4.6 — GPU-accelerated genomics on DGX Spark"

---

### ━━━ BAND 3: PRECISION INTELLIGENCE NETWORK (full width) ━━━

**Position:** Full-width horizontal band below Band 2. Light blue background tint (#E3F2FD).

**Left label (vertical text, navy #1B2333 background, white bold text):** "ENGINE 2 — Precision Intelligence"

**Contents (left to right, connected by dashed teal arrows):**

**Card 1: Milvus Retrieval** [search icon]
- Thin teal-bordered rounded rectangle
- Header: "Vector Retrieval" in bold teal
- Body: Semantic similarity search | ClinVar + AlphaMissense annotations | Evidence ranking
- Metric badge (teal pill): "3.56M vectors queried"

**Dashed arrow (teal #1AAFCC, rightward):** labeled "Top-K Evidence"

**Card 2: Precision Oncology Agent** [brain-circuit icon]
- Larger card, prominent, thin red (#DC2626) left-border accent
- Header: "Precision Oncology Agent" in bold navy
- Body: Multi-source reasoning | Variant pathogenicity classification | Therapeutic actionability scoring | Guideline-aligned interpretation
- Sub-badge row (amber pills #F5A623):
  - "ETV6-RUNX1 fusion detected"
  - "IKZF1 deletion confirmed"
- Red pill badge (#DC2626): "IKZF1plus reclassification: High Risk"

**Dashed arrow (teal, rightward):** labeled "Structured Findings"

**Card 3: RAG Evidence Synthesis** [document-stack icon]
- Thin teal-bordered rounded rectangle
- Header: "RAG Evidence Synthesis" in bold teal
- Body: Retrieval-Augmented Generation | Claude LLM reasoning | ClinVar + literature grounding | Confidence scoring
- Small annotation: "Evidence-grounded, not hallucinated"

**Bold red arrow (#DC2626) from "IKZF1plus" badge downward to Band 4:** labeled "Critical: Risk Reclassification Drives Therapy"

**Annotation box (dashed border, bottom-right of band):**
- "Key Finding: ETV6-RUNX1 is typically favorable prognosis"
- "But: IKZF1 deletion + ETV6-RUNX1 = IKZF1plus"
- "Result: Reclassified from Standard Risk to High Risk"
- "Impact: Changes treatment protocol entirely"

---

### ━━━ BAND 4: THERAPEUTIC DISCOVERY ENGINE (full width) ━━━

**Position:** Full-width horizontal band below Band 3. Light amber background tint (#FFF8E1).

**Left label (vertical text, NVIDIA green #76B900 background, white bold text):** "ENGINE 3 — Therapeutic Discovery"

**Contents (left to right, connected by solid gray arrows):**

**Card 1: MolMIM Generation** [molecule icon]
- Thin green-bordered rounded rectangle
- Header: "MolMIM Molecule Generation" in bold navy
- Body: BioNeMo foundation model | De novo molecular generation | Target: IKZF1-pathway proteins | Scaffold-constrained optimization
- Metric badge (green pill): "500 candidates generated"

**Arrow (solid gray, rightward):** labeled "Candidate Library"

**Card 2: DiffDock Scoring** [lock-and-key icon]
- Thin green-bordered rounded rectangle
- Header: "DiffDock Binding Prediction" in bold navy
- Body: Diffusion-based molecular docking | Binding affinity estimation | Pose generation + ranking
- Metric badge (green pill): "-9.1 kcal/mol best score"

**Arrow (solid gray, rightward):** labeled "Ranked by Affinity"

**Card 3: Pediatric Safety Filter** [shield icon]
- Thin emerald-bordered (#059669) rounded rectangle, light green (#D1FAE5) background tint
- Header: "Pediatric Safety Filter" in bold emerald
- Body: Molecular weight constraint | hERG cardiac liability screen | Developmental toxicity flags | Blood-brain barrier permeability | CYP450 interaction check
- Metric badges (emerald pills):
  - "MW < 500 Da"
  - "hERG > 10 uM"
- Small amber badge: "Pediatric-specific thresholds"

**Arrow (solid gray, rightward):** labeled "Safety-Cleared"

**Card 4: Top 3 Candidates** [trophy icon]
- Larger card, prominent, thin NVIDIA green (#76B900) border, light green background
- Header: "Top 3 Ranked Candidates" in bold navy
- Three sub-cards stacked vertically:
  1. **Candidate 1:** Binding: -9.1 kcal/mol | MW: 423 | hERG: 14.2 uM | RDKit: PASS
  2. **Candidate 2:** Binding: -8.7 kcal/mol | MW: 387 | hERG: 11.8 uM | RDKit: PASS
  3. **Candidate 3:** Binding: -8.3 kcal/mol | MW: 461 | hERG: 12.5 uM | RDKit: PASS
- Annotation below: "All candidates pass Lipinski's Rule of Five"

---

### ━━━ BAND 5: OUTPUT BAR (below engines) ━━━

**Position:** Full-width horizontal strip. Navy (#1B2333) background. White text. Centered.

**Center — large, bold, white text:**
"Patient DNA --> Drug Candidates in < 5 Hours"

**Below, three timing segments in a horizontal row (green pill badges on navy):**

| Genomic Foundation | Precision Intelligence | Therapeutic Discovery |
|---|---|---|
| ~2.5 hours | ~15 minutes | ~2 hours |
| FASTQ to annotated VCF | Variant interpretation + risk | Molecule generation + docking |

**Right side:** Small white text: "Traditional timeline: 6-18 months"

**Comparison arrow:** A bold visual comparison — a thin white timeline bar showing "< 5 Hours" (green) vs a much longer bar showing "6-18 Months" (red, faded), making the acceleration viscerally clear.

---

### ━━━ BAND 6: NVIDIA DGX SPARK INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | Parabricks 4.6 | BioNeMo NIMs | AI Software | Open Source |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | GPU-accelerated genomics | MolMIM molecule generation | Milvus vector database | Apache 2.0 licensed |
| Desktop AI supercomputer | BWA-MEM2 + DeepVariant | DiffDock molecular docking | BGE-small-en-v1.5 embeddings | ~1.1 TB total footprint |
| 128 GB unified LPDDR5x | < 45 min whole genome | RDKit cheminformatics | Claude LLM reasoning | Fully reproducible |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "< 45 min" on BWA-MEM2 alignment
- "> 99.7% accuracy" on DeepVariant
- "11.7M variants" on VCF output
- "3.56M vectors" on Milvus
- "500 candidates" on MolMIM
- "-9.1 kcal/mol" on DiffDock
- "MW < 500" on Pediatric Safety Filter
- "hERG > 10 uM" on Pediatric Safety Filter
- "< 5 hours end-to-end" on output bar

**Critical finding annotation** (small dashed box, connecting Band 3 to Band 4):
- "IKZF1plus Reclassification"
- ETV6-RUNX1 (favorable) + IKZF1 deletion = IKZF1plus (high risk)
- This finding changes the entire treatment approach
- Discovered by Precision Oncology Agent in < 15 minutes

**Data flow arrows style guide:**
- **Solid arrows:** Primary data flow (FASTQ through pipeline). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-agent communication and evidence retrieval. Teal (#1AAFCC) dashes.
- **Bold arrows:** Critical finding escalation paths. Red (#DC2626).

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Three engines in sequence** — Genomic Foundation, Precision Intelligence, Therapeutic Discovery — form a complete pipeline from patient DNA to drug candidates
2. **A real patient scenario** — Evelyn, 8 years old, B-ALL — grounds every element in clinical reality
3. **The IKZF1plus discovery** is the critical pivot — what appears to be a favorable-prognosis leukemia is reclassified to high risk, changing the entire treatment plan
4. **Quantitative metrics at every stage** — alignment time, accuracy, vector counts, binding scores, safety thresholds — make the pipeline credible and auditable
5. **Under 5 hours total** vs 6-18 months traditional — the acceleration is the headline
6. **Pediatric safety is not an afterthought** — a dedicated filter with age-appropriate thresholds gates every candidate
7. **Everything runs on one NVIDIA DGX Spark** — a desktop supercomputer, not a data center
8. **Open source and reproducible** — Apache 2.0 licensed, ~1.1 TB footprint, fully documented
9. **NVIDIA infrastructure anchors everything** — Parabricks, BioNeMo, DGX Spark hardware at the foundation
10. **The flow is left-to-right, top-to-bottom** — DNA enters at top-left, drug candidates emerge at bottom-right, with clear visual progression through all three engines

The overall impression should be: a complete, technically precise, production-quality precision medicine pipeline that transforms a child's DNA into ranked drug candidates in under 5 hours — dense enough to be a reference poster, clean enough to present to NVIDIA executives.
