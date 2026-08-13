# Napkin AI Pro — The 30-Second Tumor Board

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram shows five intelligence agents orchestrated by the Precision Oncology Agent, producing a comprehensive tumor board consultation in under 30 seconds — a process that traditionally takes 1-2 weeks of multidisciplinary coordination. All running on a single NVIDIA DGX Spark.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a reference architecture poster a solutions architect pins to their wall.

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
- Red: #DC2626 — critical finding indicators, bold alert arrows
- Emerald Green: #059669 — positive status indicators, GO signals
- Medium Gray: #666666 — metadata text, secondary labels

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons matching the the platform white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-agent communication, bold red (#DC2626) for critical finding paths
- Color-coded agent cards with distinct left-border accents
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer

---

## CANVAS STRUCTURE (Top to Bottom, 5 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "Demo 2 of 6" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "The 30-Second Tumor Board"
- **Subtitle line 1 (medium, gray #666666, centered):** "Five Intelligence Agents, One Orchestrated Consultation"
- **Subtitle line 2 (smaller, gray #666666, centered):** "NVIDIA DGX Spark | GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "March 2026 | HCLS AI Factory"

**Right side — Patient Card** (prominent, white background, thin red #DC2626 left-border accent, rounded corners):
```
PATIENT STATUS
——————————————————
Name: Evelyn (de-identified)
Age: 8 years old
Dx: B-ALL (ETV6-RUNX1 / IKZF1plus)
Day 29 MRD: 0.1% (positive)
Status: Post-induction, decision point
Question: What does the tumor board recommend?
```
Small thin-line icon of a clock next to a red pill badge: "Decision Needed Now"

---

### ━━━ BAND 2: CENTRAL HUB — PRECISION ONCOLOGY AGENT (center of canvas, dominant element) ━━━

**Position:** Center of the canvas, visually prominent. This is the orchestrator that coordinates all five spoke agents.

**Large central card** (double-bordered, navy #1B2333 outer border, teal #1AAFCC inner border, white interior, subtle shadow):
- **Header (bold, large, navy):** "Precision Oncology Agent"
- **Sub-header (teal):** "Orchestrator — Tumor Board Coordinator"
- **Body content in structured rows:**

  **Input Row:**
  - Patient genome: 11.7M variant records (VCF)
  - Clinical status: Day 29 MRD 0.1%
  - Key findings: ETV6-RUNX1 fusion, IKZF1 deletion
  - Risk classification: IKZF1plus (High Risk)

  **Orchestration Logic:**
  - Dispatches queries to 5 specialized agents simultaneously
  - Aggregates findings with confidence scores
  - Resolves conflicts between agent recommendations
  - Produces unified tumor board report

  **Output Row:**
  - Comprehensive tumor board consultation
  - Ranked treatment recommendations
  - Evidence citations for every finding
  - Minority opinions flagged

- Metric badge (green pill): "< 30 seconds total"
- Metric badge (teal pill): "5 agents queried in parallel"

**Five dashed teal arrows (#1AAFCC) radiate outward** from the central hub — one to each spoke agent in Band 3. Each arrow is labeled with the query sent.

---

### ━━━ BAND 3: FIVE SPOKE AGENTS (arranged in a semicircle or horizontal row around the hub) ━━━

**Position:** Five agent cards arranged in a horizontal row below (or surrounding) the central Precision Oncology Agent hub. Each card has a distinct left-border color accent, and each connects back to the hub via dashed teal arrows (outbound query, inbound findings).

**Agent Card 1: Biomarker Discovery Agent** [DNA helix icon]
- Left-border accent: NVIDIA green (#76B900)
- Thin green-bordered rounded rectangle, light green tint (#E8F5E9)
- Header: "Biomarker Discovery Agent" in bold navy
- Query from hub (label on inbound arrow): "Identify actionable biomarkers"
- **Findings (stacked vertically):**
  - Amber pill: "ETV6-RUNX1 fusion — favorable baseline"
  - Red pill: "IKZF1 deletion — modifier detected"
  - Red pill: "IKZF1plus composite — HIGH RISK reclassification"
  - Emerald pill: "No TP53 mutation — good prognostic sign"
- Return arrow label: "4 biomarker findings"

**Agent Card 2: CAR-T Cell Therapy Agent** [T-cell icon]
- Left-border accent: Teal (#1AAFCC)
- Thin teal-bordered rounded rectangle, light blue tint (#E3F2FD)
- Header: "CAR-T Cell Therapy Agent" in bold navy
- Query from hub: "Assess CAR-T eligibility if needed"
- **Findings:**
  - Emerald pill: "CD19 expression: POSITIVE (98%)"
  - Amber pill: "CD22 co-expression: 67% — backup target"
  - Teal pill: "Tisagenlecleucel: eligible if relapse"
  - Gray pill: "Current status: Not indicated yet (first-line)"
- Return arrow label: "CAR-T readiness assessed"

**Agent Card 3: Single-Cell Intelligence Agent** [cell cluster icon]
- Left-border accent: Amber (#F5A623)
- Thin amber-bordered rounded rectangle, light amber tint (#FFF8E1)
- Header: "Single-Cell Intelligence Agent" in bold navy
- Query from hub: "Characterize tumor heterogeneity"
- **Findings:**
  - Amber pill: "3 distinct blast subpopulations identified"
  - Red pill: "Minor clone (8%) with drug-resistant signature"
  - Amber pill: "MRD population: phenotypically distinct from bulk"
  - Teal pill: "Recommended: monitor resistant clone evolution"
- Return arrow label: "Clonal architecture mapped"

**Agent Card 4: Clinical Trial Intelligence Agent** [clipboard-search icon]
- Left-border accent: Emerald (#059669)
- Thin emerald-bordered rounded rectangle, light emerald tint (#D1FAE5)
- Header: "Clinical Trial Intelligence Agent" in bold navy
- Query from hub: "Find matching trials for IKZF1plus B-ALL"
- **Findings:**
  - Emerald pill: "AALL1732 — IKZF1plus-specific arm: MATCH"
  - Emerald pill: "COG AALL0434 — intensified therapy: ELIGIBLE"
  - Amber pill: "NCT04697706 — blinatumomab consolidation: OPEN"
  - Teal pill: "3 trials matched, 2 actively enrolling"
- Return arrow label: "3 trial matches"

**Agent Card 5: Therapeutic Discovery Agent** [molecule icon]
- Left-border accent: NVIDIA green (#76B900)
- Thin green-bordered rounded rectangle, light green tint (#E8F5E9)
- Header: "Therapeutic Discovery Agent" in bold navy
- Query from hub: "Generate candidates targeting IKZF1 pathway"
- **Findings:**
  - Green pill: "MolMIM: 500 candidates generated"
  - Green pill: "Top candidate: -9.1 kcal/mol binding"
  - Emerald pill: "3 candidates pass pediatric safety filter"
  - Amber pill: "Novel IKZF1-pathway modulators identified"
- Return arrow label: "3 drug candidates ranked"

---

### ━━━ BAND 4: COMPARISON BAR (below spoke agents) ━━━

**Position:** Full-width horizontal strip. Split into two halves with a bold vertical divider.

**Left half — "Traditional Tumor Board" (red-tinted background #FEE2E2):**
- Header (bold, red #DC2626): "Traditional Multidisciplinary Tumor Board"
- Content rows with thin-line icons:
  - [calendar icon] Schedule coordination: 3-5 days
  - [people icon] 6-12 specialists must align calendars
  - [document icon] Manual literature review: 2-4 hours per specialist
  - [clock icon] Meeting duration: 30-60 minutes
  - [hourglass icon] Total time to recommendations: **1-2 weeks**
- Large red pill badge: "1-2 WEEKS"

**Right half — "HCLS AI Factory Tumor Board" (green-tinted background #D1FAE5):**
- Header (bold, emerald #059669): "HCLS AI Factory — AI Tumor Board"
- Content rows with thin-line icons:
  - [lightning icon] All 5 agents queried in parallel
  - [brain icon] 3.56M vectors searched simultaneously
  - [document icon] Evidence-grounded, citation-backed findings
  - [shield icon] Confidence scores on every recommendation
  - [clock icon] Total time to recommendations: **< 30 seconds**
- Large emerald pill badge: "< 30 SECONDS"

**Center divider:** Bold vertical line with a large "vs" in a circle. Below: "2,880x faster" in bold navy text.

---

### ━━━ BAND 5: NVIDIA DGX SPARK INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | Agent Framework | Vector Intelligence | LLM Reasoning | Open Source |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | 5 specialized intelligence agents | Milvus: 3.56M vectors | Claude via Anthropic API | Apache 2.0 licensed |
| Desktop AI supercomputer | Parallel orchestration | BGE-small-en-v1.5 embeddings | RAG evidence grounding | ~1.1 TB total footprint |
| 128 GB unified LPDDR5x | FastAPI microservices | ClinVar + AlphaMissense | Confidence scoring | Fully reproducible |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "< 30 seconds" on central hub
- "5 agents in parallel" on orchestrator
- "3.56M vectors" on vector search references
- "11.7M variant records" on patient genome input
- "2,880x faster" on comparison bar
- "98% CD19+" on CAR-T agent
- "3 trial matches" on Clinical Trial agent
- "500 candidates" on Therapeutic Discovery agent

**Orchestration annotation** (small dashed box, below central hub):
- "Parallel Dispatch Pattern"
- All 5 agents receive queries simultaneously
- Results aggregated with confidence weighting
- Conflicting findings flagged for clinician review
- Complete audit trail preserved

**Data flow arrows style guide:**
- **Solid arrows:** Primary data flow (patient data into orchestrator). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-agent communication (hub to spoke and back). Teal (#1AAFCC) dashes.
- **Bold arrows:** Critical finding escalation (IKZF1plus reclassification). Red (#DC2626).

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **One orchestrator, five specialist agents** — the Precision Oncology Agent coordinates Biomarker, CAR-T, Single-Cell, Clinical Trial, and Therapeutic Discovery agents in parallel
2. **Hub-and-spoke architecture** — the central hub dispatches, aggregates, and resolves findings from all five agents
3. **A real patient at a real decision point** — Evelyn, Day 29 MRD-positive, needs a tumor board recommendation NOW
4. **Each agent delivers specific, actionable findings** — not vague summaries, but concrete biomarkers, trial matches, drug candidates, and clonal architecture
5. **Under 30 seconds** vs 1-2 weeks traditional — the comparison bar makes this viscerally clear
6. **Every finding is evidence-grounded** — RAG retrieval, confidence scores, citations
7. **The IKZF1plus reclassification** threads through multiple agents as the critical finding
8. **CAR-T readiness is pre-assessed** — even before it's needed, the system has already evaluated eligibility
9. **Clinical trials are matched in real time** — not weeks of manual searching
10. **Everything runs on one NVIDIA DGX Spark** — anchored by the infrastructure bar at the bottom

The overall impression should be: a complete, technically precise tumor board consultation delivered by five coordinated AI agents in under 30 seconds — replacing weeks of manual coordination with a parallel, evidence-grounded intelligence network on a single desktop supercomputer.
