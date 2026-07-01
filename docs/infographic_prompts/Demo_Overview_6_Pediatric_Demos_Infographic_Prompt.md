# Nano Banana Pro — HCLS AI Factory: 6 Pediatric Oncology Demos Overview

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical overview infographic — NOT a slide deck, NOT six separate diagrams. Every element described below appears on ONE canvas. This diagram is a visual summary of all six pediatric oncology demonstration workflows running on the HCLS AI Factory platform — showing how 3 engines and 11 AI agents coordinate across 5 children with cancer. It should be the single image someone sees and immediately understands the scope, depth, and clinical impact of the platform.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional clinical technology infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished clinical research poster — clean, structured, authoritative, emotionally grounded. Match the aesthetic of the HCLS AI Factory white paper series (bold sans-serif headers, structured card-based layouts with subtle borders, metric badges) combined with the clinical warmth of a pediatric oncology conference poster where every data point represents a child.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. Six demo cards arranged as the dominant visual element, each telling a complete patient story. The diagram should feel like a reference poster that a pediatric oncologist, a bioinformatician, and an NVIDIA executive could all read and understand.

**Typography:**
- Title: Large, bold, sans-serif (Inter, Helvetica, or similar), deep navy (#1B2333)
- Subtitle: Smaller, medium gray (#666666), directly below title — emotional tagline
- Section headers: Bold, dark navy (#1B2333) on white, with a thin NVIDIA green (#76B900) left-border accent or underline
- Sub-headers: Bold, teal (#1AAFCC)
- Demo titles (short names): Bold, large, navy (#1B2333) — these are the emotional hooks
- Demo subtitles (long names): Medium, teal (#1AAFCC), directly below short names
- Patient lines: Bold, teal (#1AAFCC), slightly smaller — name, age, diagnosis
- Body text: Small (8-10pt equivalent), clean sans-serif, dark gray (#333333)
- Agent flow text: Small, NVIDIA green (#76B900), monospace or distinct from body
- Metric callouts: Bold, slightly larger than body, inside small rounded green (#76B900) or teal (#1AAFCC) pill badges with white text
- Footer text: Small, gray (#999999)

**Color Palette (exact — use these hex values):**
- NVIDIA Green: #76B900 — primary accent for metric badges, agent flow lines, verification bar, infrastructure
- Deep Navy: #1B2333 — title text, demo short names, dark section bars
- Teal: #1AAFCC — sub-headers, patient lines, demo long names, agent connection lines
- Light Gray: #F5F5F5 — demo card backgrounds
- White: #FFFFFF — canvas, text on dark backgrounds
- Amber/Orange: #F5A623 — warning indicators (vismodegib contraindication, CRS risk)
- Red: #DC2626 — critical finding indicators (IKZF1plus reclassification, failed prior lines)
- Purple: #7B2D8E — rare disease / genetic cascade indicators
- Medium Gray: #666666 — metadata text, secondary labels
- Emerald Green: #059669 — positive outcomes, GO checkmarks, safety clearances

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every demo card and component box
- Thin-line icons (16x16 to 24x24) next to each demo — simple, monochrome line icons (not emoji): DNA helix for Demo 1, group/team for Demo 2, shield for Demo 3, family tree for Demo 4, T-cell for Demo 5, molecule for Demo 6
- Directional arrows: solid medium gray (#999999) for agent flow sequences within each demo
- Color-coded demo cards with distinct subtle background tints per domain
- Metric badges: small rounded pills with white text on green or teal background
- Patient avatar silhouettes: small child silhouette icons (age-appropriate) in each demo card header
- NVIDIA logo mark in the title bar and infrastructure footer
- "Apache 2.0" badge prominently visible

---

## CANVAS STRUCTURE (Top to Bottom, 5 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below: "Pediatric Oncology" in teal (#1AAFCC) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "6 Pediatric Oncology Demos"
- **Tagline (medium, teal #1AAFCC, centered, italic):** "Built for the hardest cases in medicine. Works for all of them."
- **Description (smaller, gray #666666, centered):** "Each demo tells the story of a child with cancer moving through the full HCLS AI Factory pipeline. All 3 engines and all 8 agents are exercised across the 6 workflows."

**Right side — Platform Summary box** (small, top-right corner, thin navy border, white background):
```
Platform
━━━━━━━━━━━━━━━━━━━
3 GPU-Accelerated Engines
11 Intelligence Agents
250,000+ Lines of Code
5,000+ Automated Tests
Apache 2.0 Open Source
NVIDIA DGX Spark
```
Use small NVIDIA green dots (●) before each line.

---

### ━━━ BAND 2: THREE-ENGINE STRIP (narrow horizontal band below title) ━━━

**Position:** Full-width horizontal strip. Three equal columns showing the three engines as a visual pipeline. Light gray (#F5F5F5) background with NVIDIA green (#76B900) top and bottom border lines.

**Three columns, connected by bold green arrows (→):**

**Column 1 — Genomic Foundation Engine**
- Icon: DNA helix (thin-line, navy)
- Header: "Genomic Foundation Engine" (bold, navy)
- Line 1: "FASTQ → VCF" in teal
- Line 2: "Parabricks • DeepVariant • BWA-MEM2" (small, gray)
- Metric badge: "2-4 hrs" (green pill) vs "24-48 hrs CPU" (gray pill, strikethrough)

**Bold green arrow →**

**Column 2 — Precision Intelligence Engine**
- Icon: Network/brain (thin-line, navy)
- Header: "Precision Intelligence Engine" (bold, navy)
- Line 1: "11 AI Agents" in teal
- Line 2: "Milvus • Claude LLM • 3.56M Vectors" (small, gray)
- Metric badge: "Cross-Agent Coordination" (teal pill)

**Bold green arrow →**

**Column 3 — Therapeutic Discovery Engine**
- Icon: Molecule (thin-line, navy)
- Header: "Therapeutic Discovery Engine" (bold, navy)
- Line 1: "Targets → Drug Candidates" in teal
- Line 2: "MolMIM • DiffDock • RDKit" (small, gray)
- Metric badge: "Pediatric Safety Filters" (green pill)

---

### ━━━ BAND 3: SIX DEMO CARDS (the main visual — largest section, ~55% of canvas) ━━━

**Position:** Center of canvas. Six demo cards arranged in a 3×2 grid (3 columns, 2 rows). Each card is a rounded-corner rectangle (8px radius) with a light background tint, thin left-border accent in the demo's theme color, and structured internal layout.

**Card dimensions:** Equal width, roughly 30% of canvas width each. Equal height. Generous internal padding. Each card is self-contained — a reader should understand the demo from the card alone.

**Grid layout:**

| Row 1 | Demo 1 (top-left) | Demo 2 (top-center) | Demo 3 (top-right) |
| Row 2 | Demo 4 (bottom-left) | Demo 5 (bottom-center) | Demo 6 (bottom-right) |

---

#### DEMO CARD 1 (top-left) — Background tint: #E8F5E9 (light green), left border: #76B900

**Card header:**
- Icon: DNA helix (thin-line, 24x24, navy)
- Short name (bold, navy, 14pt equivalent): "DNA to Drug"
- Long name (teal, 10pt): "End-to-End Precision Medicine Pipeline"

**Patient line (teal, bold, 9pt):**
"Evelyn, Age 8 | Newly Diagnosed B-Cell ALL"

**Narrative (dark gray #333333, 8pt, 4-5 lines):**
"Evelyn's tumor and germline DNA enter the Genomic Foundation Engine as FASTQ files. Parabricks aligns against GRCh38 in 45 minutes. DeepVariant identifies an ETV6-RUNX1 fusion (favorable) with an IKZF1 deletion (unfavorable). The Oncology Agent reclassifies to high risk (IKZF1plus). The Therapeutic Discovery Engine generates novel CREBBP inhibitors with pediatric safety filters."

**Key findings (3 small pill badges in a row):**
- "ETV6-RUNX1" (green pill)
- "IKZF1plus → High Risk" (red pill)
- "<5 hrs DNA to Drug" (green pill)

**Agent flow (small, green #76B900, bottom of card):**
"Genomic Foundation → Oncology → RAG → Therapeutic Discovery"

---

#### DEMO CARD 2 (top-center) — Background tint: #E3F2FD (light blue), left border: #1AAFCC

**Card header:**
- Icon: Group/team (thin-line, 24x24, navy)
- Short name: "The 30-Second Tumor Board"
- Long name: "Multi-Agent Coordination"

**Patient line:** "Evelyn, Day 29 | MRD Positive (0.1%)"

**Narrative:**
"Five agents convene a virtual tumor board in under 30 seconds. Biomarker quantifies CD19 (MFI 45,200) and MRD kinetics. CAR-T evaluates tisagenlecleucel eligibility. Single-Cell identifies a 21% CD19-dim escape subclone. Clinical Trial matches 3 COG protocols. What takes 1-2 weeks to schedule happens in a single query."

**Key findings:**
- "MRD 0.1%" (amber pill)
- "CD19-dim 21% Escape" (red pill)
- "<30 sec 5 Agents" (teal pill)

**Agent flow:**
"Oncology → Biomarker → CAR-T → Single-Cell → Clinical Trial → Therapeutic Discovery"

---

#### DEMO CARD 3 (top-right) — Background tint: #FFF8E1 (light amber), left border: #F5A623

**Card header:**
- Icon: Shield (thin-line, 24x24, navy)
- Short name: "Protecting the Survivor"
- Long name: "Cardiotoxicity Prevention"

**Patient line:** "Marcus, Age 6 | High-Risk Neuroblastoma (MYCN Amplified)"

**Narrative:**
"Before the first dose, 7 agents coordinate. PGx: TPMT *1/*3A → 50% 6-MP reduction. Biomarker: MYCN amplified, elevated LDH. Cardiology: dexrazoxane at 300 mg/m². Neurology: vincristine 2mg cap. Autoimmune: dinutuximab irAEs (pain 85%, capillary leak 25%). Therapeutic Discovery: novel ALK inhibitors with cardiac safety. 85% of survivors develop chronic conditions by 40."

**Key findings:**
- "TPMT 50% Dose" (green pill)
- "Dexrazoxane Required" (amber pill)
- "7 Agents Protecting" (green pill)

**Agent flow:**
"PGx → Biomarker → Cardiology → Neurology → Autoimmune → Oncology → Therapeutic Discovery"

---

#### DEMO CARD 4 (bottom-left) — Background tint: #F3E8FF (light purple), left border: #7B2D8E

**Card header:**
- Icon: Family tree (thin-line, 24x24, navy)
- Short name: "One Gene, One Family"
- Long name: "Rare Disease + Cancer Predisposition"

**Patient line:** "Aurora, Age 4 | Bilateral Retinoblastoma"

**Narrative:**
"Rare Disease classifies RB1 c.958C>T as Pathogenic (ACMG: PVS1+PM2+PP4). Her 2-year-old sibling has 50% chance — cascade testing critical. Oncology: globe-sparing therapy, avoid radiation (3-4× second cancer risk). Imaging: trilateral screening — brain MRI every 6 months until age 5. Therapeutic Discovery: CDK4/6 inhibitors for ocular delivery. One diagnosis changes one family forever."

**Key findings:**
- "RB1 Pathogenic" (purple pill)
- "Sibling 50% Risk" (red pill)
- "Lifetime Surveillance" (teal pill)

**Agent flow:**
"Rare Disease → Oncology → Imaging → Clinical Trial → Therapeutic Discovery"

---

#### DEMO CARD 5 (bottom-center) — Background tint: #FFEBEE (light red), left border: #DC2626

**Card header:**
- Icon: T-cell / engineered cell (thin-line, 24x24, navy)
- Short name: "Last Line of Defense"
- Long name: "CAR-T Therapy Decision"

**Patient line:** "Ethan, Age 12 | Relapsed/Refractory B-ALL (Failed 2 Lines)"

**Narrative:**
"Seven agents evaluate CAR-T candidacy. Single-Cell: CD19 97.2% (MFI 8,500), immune desert TME. CAR-T: tisagenlecleucel eligible, CRS moderate, CD22 backup. Cardiology: LVEF 58%, cleared for conditioning. Autoimmune: CRS monitoring, 30-40% cytopenia risk. PGx: all supportive meds cleared. Clinical Trial: ELIANA, ZUMA-4, dual CD19/CD22. Therapeutic Discovery: bridging candidates for 22-day manufacturing."

**Key findings:**
- "CD19+ 97.2%" (green pill)
- "CRS Moderate" (amber pill)
- "7 Gates Cleared" (emerald pill)

**Agent flow:**
"Single-Cell → CAR-T → Cardiology → Autoimmune → PGx → Clinical Trial → Therapeutic Discovery"

---

#### DEMO CARD 6 (bottom-right) — Background tint: #E0F2F1 (light teal), left border: #059669

**Card header:**
- Icon: Molecule / flask (thin-line, 24x24, navy)
- Short name: "When the Standard Drug Can't Be Used"
- Long name: "Novel Drug Discovery for a Growing Child"

**Patient line:** "Aiden, Age 10 | SHH Medulloblastoma (PTCH1 Mutation)"

**Narrative:**
"Oncology flags CRITICAL: vismodegib causes permanent growth plate fusion in children. Neurology: posterior fossa syndrome 25-30%, proton therapy recommended. Imaging: T3aM0, spinal surveillance. Therapeutic Discovery generates 500 SMO antagonists via MolMIM (PDB: 5L7D), pediatric safety filter eliminates 411, DiffDock scores 89, top 3 achieve IC50 8-18 nM with growth plate safety. Clinical Trial: SJMB12, Pediatric MATCH."

**Key findings:**
- "Vismodegib BLOCKED" (red pill)
- "500 → 89 → 3 Candidates" (green pill)
- "Growth Plate Safe" (emerald pill)

**Agent flow:**
"Oncology → Neurology → Imaging → Therapeutic Discovery → Clinical Trial"

---

### ━━━ BAND 4: VERIFICATION & METRICS BAR (below demo cards) ━━━

**Position:** Full-width horizontal strip. NVIDIA green (#76B900) background. White bold text. This is the credibility bar.

**Single row of metrics, evenly spaced, centered:**

| All 3 Engines | All 11 Agents | All 6 Demos | 25/25 Claims Verified | 24,043 Files | Zero Failures |
|---|---|---|---|---|---|

Each metric rendered as:
- Large number/text (white, bold, 16pt equivalent)
- Small label below (white, 8pt)
- Thin white vertical divider lines between each metric

---

### ━━━ BAND 5: FOOTER BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. Navy (#1B2333) background. White and gray text. NVIDIA logo mark on left side.

**Left:** NVIDIA logo mark (small) + "Built for NVIDIA DGX Spark"

**Center:** "Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org"

**Right:** "March 2026 | HCLS AI Factory v1.0.0"

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Demo card visual consistency:**
- All 6 cards must be exactly the same dimensions
- Internal layout must be identical across cards (icon position, name position, patient line position, narrative area, badges row, agent flow line)
- The only differences between cards are: background tint color, left border color, icon, and content
- This visual uniformity reinforces that these are 6 instances of one platform, not 6 different products

**Agent flow arrows within each card:**
- Small, thin, green (#76B900) text showing the agent sequence
- Use → characters between agent names
- Positioned consistently at the bottom of each card

**Metric badge pill style:**
- Rounded pill shape (16px height, auto-width based on text)
- White text on colored background
- Green (#76B900) for positive findings
- Teal (#1AAFCC) for informational
- Amber (#F5A623) for warnings
- Red (#DC2626) for critical findings
- Purple (#7B2D8E) for genetic/rare disease
- Emerald (#059669) for safety clearances

**Connection between Band 2 (engines) and Band 3 (demos):**
- Thin dashed teal (#1AAFCC) lines connecting each engine column to the demo cards that use it
- Since all demos use all 3 engines, these lines fan out from each engine to all 6 cards
- This visually reinforces the "all engines, all demos" message

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Six children with cancer** — each demo has a name, an age, and a diagnosis. These are human stories, not technical abstractions.
2. **Three engines, one pipeline** — raw DNA enters on the left, drug candidates exit on the right, intelligence happens in the middle.
3. **Eleven agents, fully coordinated** — every demo shows multiple agents working together, not in isolation.
4. **Every demo uses all three engines** — the Therapeutic Discovery Engine appears in every workflow, generating novel drug candidates.
5. **Clinical accuracy is verified** — 25/25 claims verified, 24,043 files with zero failures. This is not a prototype.
6. **Pediatric-specific** — growth plate safety, BSA-scaled dosing, TPMT genotyping, dexrazoxane cardioprotection. Every detail matters for a growing child.
7. **Speed** — DNA to drug in <5 hours, tumor board in <30 seconds. The traditional timeline is measured in weeks or months.
8. **Open source** — Apache 2.0, free for every hospital on Earth. This is prominently displayed, not hidden.
9. **Built on NVIDIA** — DGX Spark, Parabricks, BioNeMo, NIMs. The infrastructure is real and identified.
10. **One person built this** — "Created by Adam Jones" in the footer. The scope-to-individual ratio is part of the story.

The overall impression should be: a comprehensive, clinically rigorous, emotionally grounded overview of six precision medicine workflows for children with cancer — dense enough to be a conference poster, clean enough to be a one-page executive handout, and human enough to make someone stop scrolling.
