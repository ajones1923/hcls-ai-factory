# Napkin AI Pro — Last Line of Defense: CAR-T Eligibility Decision Tree

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram depicts a vertical decision tree with seven sequential gates — each gate is an intelligence agent that must return GO before the patient can proceed to CAR-T cell therapy. This is a portrait-oriented (vertical) infographic, reflecting the top-to-bottom decision cascade. The patient is a 12-year-old with relapsed/refractory B-ALL who has exhausted all other options. All running on a single NVIDIA DGX Spark.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in portrait orientation (9:16 aspect ratio, vertical). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a clinical decision algorithm — a step-by-step clearance process that must complete top-to-bottom before therapy can begin.

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
- Amber/Orange: #F5A623 — caution indicators, conditional-pass findings
- Red: #DC2626 — critical findings, urgency indicators, NO-GO signals
- Emerald Green: #059669 — GO signals, cleared status, checkmarks
- Medium Gray: #666666 — metadata text, secondary labels

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons matching the the platform white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary decision flow (downward), dashed teal (#1AAFCC) for cross-agent communication, bold red (#DC2626) for urgency paths
- GO/NO-GO badges: large, prominent, emerald (#059669) checkmarks for GO, red (#DC2626) X marks for NO-GO
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer

---

## CANVAS STRUCTURE (Top to Bottom, 10 horizontal bands — portrait orientation)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "Demo 5 of 6" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "Last Line of Defense"
- **Subtitle line 1 (medium, gray #666666, centered):** "Seven-Gate CAR-T Eligibility Decision Tree"
- **Subtitle line 2 (smaller, gray #666666, centered):** "NVIDIA DGX Spark | GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "March 2026 | HCLS AI Factory"

**Right side — Patient Card** (prominent, white background, thin red #DC2626 left-border accent, rounded corners, urgency indicator):
```
PATIENT — URGENT
——————————————————
Name: Ethan (de-identified)
Age: 12 years old
Dx: Relapsed/Refractory B-ALL
Prior Lines: 3 (induction, re-induction, salvage)
Status: All standard therapies exhausted
MRD: 4.7% (rising)
Question: Is Ethan eligible for CAR-T therapy?
```
Red pill badge (pulsing/bold): "LAST LINE OF DEFENSE"
Amber pill badge: "Time-critical decision"

---

### ━━━ BAND 2: GATE 1 — SINGLE-CELL INTELLIGENCE AGENT ━━━

**Position:** Full-width horizontal band. Light amber background tint (#FFF8E1). Numbered gate badge on left.

**Left column — Gate Badge:**
- Large circle with "1" inside, navy background, white text
- Label: "GATE 1" in bold navy
- Sub-label: "Tumor Characterization" in teal

**Center — Agent Card:** [cell cluster icon]
- Header: "Single-Cell Intelligence Agent" in bold navy
- Body: Characterizes blast population at single-cell resolution | Identifies clonal heterogeneity | Assesses antigen expression for CAR-T targeting

**Right column — Findings + GO/NO-GO:**
- Amber pill: "3 blast subpopulations identified"
- Red pill: "Minor clone (12%) with CD19-dim phenotype"
- Emerald pill: "Dominant population (88%): CD19-bright"
- Teal pill: "Clonal diversity: moderate — monitor for antigen escape"
- **Large GO badge:** Emerald (#059669) rounded rectangle with white checkmark and bold text "GO"
  - Annotation: "Dominant clone is CD19-bright — targetable"

**Solid gray arrow pointing downward to Gate 2**

---

### ━━━ BAND 3: GATE 2 — CAR-T CELL THERAPY AGENT ━━━

**Position:** Full-width horizontal band. Light blue background tint (#E3F2FD).

**Left column — Gate Badge:**
- Large circle with "2" inside
- Label: "GATE 2" in bold navy
- Sub-label: "CAR-T Eligibility" in teal

**Center — Agent Card:** [T-cell icon]
- Header: "CAR-T Cell Therapy Agent" in bold navy
- Body: Evaluates target antigen expression | Assesses manufacturing feasibility | Selects CAR-T product | Checks prior CAR-T history

**Right column — Findings + GO/NO-GO:**
- Emerald pill: "CD19 expression: 92% of blasts — ABOVE threshold"
- Emerald pill: "No prior CAR-T therapy — first-time eligible"
- Emerald pill: "Tisagenlecleucel (Kymriah): FDA-approved for r/r B-ALL age <25"
- Teal pill: "Absolute lymphocyte count: adequate for apheresis"
- Amber pill: "CD22 backup: 71% expression if CD19 escape"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "CD19+ dominant, Kymriah eligible, no prior CAR-T"

**Solid gray arrow pointing downward to Gate 3**

---

### ━━━ BAND 4: GATE 3 — CARDIOLOGY INTELLIGENCE AGENT ━━━

**Position:** Full-width horizontal band. Light red background tint (#FEE2E2).

**Left column — Gate Badge:**
- Large circle with "3" inside
- Label: "GATE 3" in bold navy
- Sub-label: "Cardiac Clearance" in teal

**Center — Agent Card:** [heart icon]
- Header: "Cardiology Intelligence Agent" in bold navy
- Body: Screens for pre-existing cardiac conditions | Cytokine Release Syndrome (CRS) cardiac risk | Anthracycline cumulative dose assessment | Echocardiographic baseline

**Right column — Findings + GO/NO-GO:**
- Amber pill: "Cumulative doxorubicin: 300 mg/m2 — moderate exposure"
- Emerald pill: "Echocardiogram: EF 58% — above CRS threshold (>50%)"
- Emerald pill: "No arrhythmias on baseline ECG"
- Amber pill: "CRS cardiac risk: MODERATE — requires ICU monitoring"
- Teal pill: "Troponin and BNP baselines: within normal limits"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "Cardiac function sufficient — ICU monitoring recommended during CRS"

**Solid gray arrow pointing downward to Gate 4**

---

### ━━━ BAND 5: GATE 4 — AUTOIMMUNE INTELLIGENCE AGENT ━━━

**Position:** Full-width horizontal band. Light amber background tint (#FFF8E1).

**Left column — Gate Badge:**
- Large circle with "4" inside
- Label: "GATE 4" in bold navy
- Sub-label: "Autoimmune Clearance" in teal

**Center — Agent Card:** [antibody icon]
- Header: "Autoimmune Intelligence Agent" in bold navy
- Body: Screens for pre-existing autoimmune conditions | HLA-mediated risk profiling | B-cell aplasia tolerance assessment | Immunoglobulin replacement planning

**Right column — Findings + GO/NO-GO:**
- Emerald pill: "No pre-existing autoimmune conditions"
- Emerald pill: "HLA profile: no high-risk autoimmune alleles"
- Amber pill: "Expected: prolonged B-cell aplasia post-CAR-T"
- Teal pill: "IVIG replacement protocol: pre-authorized"
- Emerald pill: "Infection screening: CMV/EBV/HBV — all negative or controlled"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "No autoimmune contraindications — IVIG plan in place"

**Solid gray arrow pointing downward to Gate 5**

---

### ━━━ BAND 6: GATE 5 — PHARMACOGENOMICS (PGx) AGENT ━━━

**Position:** Full-width horizontal band. Light purple background tint (#F3E8FF).

**Left column — Gate Badge:**
- Large circle with "5" inside
- Label: "GATE 5" in bold navy
- Sub-label: "Pharmacogenomics" in teal

**Center — Agent Card:** [pill-DNA icon]
- Header: "Pharmacogenomics (PGx) Agent" in bold navy
- Body: Screens variants affecting CAR-T supportive medications | CRS management drug metabolism | Tocilizumab/corticosteroid response prediction

**Right column — Findings + GO/NO-GO:**
- Emerald pill: "CYP3A4 *1/*1 — Normal metabolizer (tocilizumab OK)"
- Emerald pill: "CYP2D6 *1/*2 — Normal metabolizer"
- Amber pill: "IL6R variant detected — may affect tocilizumab binding"
- Teal pill: "RECOMMEND: standard tocilizumab dosing with close monitoring"
- Emerald pill: "No contraindications to supportive care regimen"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "Supportive medications cleared — monitor tocilizumab response"

**Solid gray arrow pointing downward to Gate 6**

---

### ━━━ BAND 7: GATE 6 — CLINICAL TRIAL INTELLIGENCE AGENT ━━━

**Position:** Full-width horizontal band. Light emerald background tint (#D1FAE5).

**Left column — Gate Badge:**
- Large circle with "6" inside
- Label: "GATE 6" in bold navy
- Sub-label: "Trial Matching" in teal

**Center — Agent Card:** [clipboard-search icon]
- Header: "Clinical Trial Intelligence Agent" in bold navy
- Body: Matches patient to CAR-T clinical trials | Checks if commercial product or trial is optimal | Identifies post-CAR-T maintenance trials

**Right column — Findings + GO/NO-GO:**
- Emerald pill: "Tisagenlecleucel: FDA-approved — commercial product available"
- Amber pill: "NCT04195789 — CD19/CD22 bispecific CAR-T: ELIGIBLE (Phase II)"
- Amber pill: "NCT05480449 — Post-CAR-T blinatumomab maintenance: ELIGIBLE"
- Teal pill: "Recommendation: Proceed with commercial Kymriah + enroll in maintenance trial"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "Commercial product available — maintenance trial identified"

**Solid gray arrow pointing downward to Gate 7**

---

### ━━━ BAND 8: GATE 7 — THERAPEUTIC DISCOVERY AGENT ━━━

**Position:** Full-width horizontal band. Light green background tint (#E8F5E9).

**Left column — Gate Badge:**
- Large circle with "7" inside
- Label: "GATE 7" in bold navy
- Sub-label: "Backup Therapeutics" in teal

**Center — Agent Card:** [molecule icon]
- Header: "Therapeutic Discovery Agent" in bold navy
- Body: Generates backup compounds if CAR-T fails | Targets resistant clone (CD19-dim 12%) | Pre-computes escape-scenario therapeutics

**Right column — Findings + GO/NO-GO:**
- Green pill: "MolMIM: 500 candidates targeting CD19-dim pathway"
- Green pill: "Top candidate: -8.6 kcal/mol, MW 398"
- Emerald pill: "3 candidates pass pediatric safety filter"
- Amber pill: "Backup strategy: ready if antigen escape occurs"
- **Large GO badge:** Emerald (#059669) "GO"
  - Annotation: "Escape-scenario therapeutics pre-computed"

**Bold solid green arrow pointing downward to Band 9**

---

### ━━━ BAND 9: DECISION OUTPUT (prominent, high-impact band) ━━━

**Position:** Full-width horizontal band. Emerald green (#059669) background. White text. This is the culmination of the entire decision tree.

**Center — large, bold, white text:**
"ALL 7 GATES: GO"

**Below, in slightly smaller white bold text:**
"PROCEED WITH TISAGENLECLEUCEL (KYMRIAH)"

**Below, three output cards in a horizontal row (white background, emerald border):**

| Treatment Plan | Safety Monitoring | Backup Strategy |
|---|---|---|
| Tisagenlecleucel (commercial) | ICU monitoring during CRS window | CD19-dim escape: 3 candidates ready |
| Lymphodepleting conditioning: Flu/Cy | Tocilizumab on standby (monitor IL6R) | CD22-directed therapy: backup target |
| Target infusion: Day 0 | IVIG replacement post-B-cell aplasia | Maintenance trial: NCT05480449 |
| Apheresis → Manufacturing → Infusion | Cardiac monitoring: EF checks Day 7, 14, 28 | Bispecific CAR-T trial: eligible |

**Below cards, a single-line bold annotation:**
"7 intelligence agents, 7 gates cleared, 1 decision — in under 5 minutes"

---

### ━━━ BAND 10: NVIDIA DGX SPARK INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | Agent Framework | Vector Intelligence | Drug Discovery | Open Source |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | 7 sequential decision gates | Milvus: 3.56M vectors | MolMIM escape-scenario compounds | Apache 2.0 licensed |
| Desktop AI supercomputer | GO/NO-GO at every gate | ClinVar + AlphaMissense | DiffDock + RDKit scoring | ~1.1 TB total footprint |
| 128 GB unified LPDDR5x | Complete audit trail | BGE-small-en-v1.5 embeddings | Pediatric safety filter | Fully reproducible |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "7/7 gates cleared" on the decision output
- "< 5 minutes total" on the decision output
- "92% CD19+" on CAR-T agent
- "EF 58%" on Cardiology agent
- "500 backup candidates" on Therapeutic Discovery agent
- "3 trials matched" on Clinical Trial agent
- "47 pharmacogenes screened" on PGx agent

**Urgency annotation** (small dashed box, top-right area):
- "TIME-CRITICAL DECISION"
- MRD rising: 4.7%
- 3 prior therapy lines exhausted
- CAR-T manufacturing takes 3-4 weeks
- Every day of delay = disease progression
- This decision tree runs in < 5 minutes vs days of manual review

**Gate failure annotation** (small dashed box, left side of canvas):
- "If any gate returns NO-GO:"
- Decision tree halts at that gate
- Alternative pathways presented
- Clinical team notified immediately
- All findings preserved in audit trail
- NO-GO does not mean no treatment — it means different treatment

**Data flow arrows style guide:**
- **Solid arrows:** Primary decision flow (top to bottom through gates). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-agent data sharing between gates. Teal (#1AAFCC) dashes.
- **Bold arrows:** Critical urgency indicators and final decision output. Red (#DC2626) for urgency, Emerald (#059669) for GO.

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Seven sequential gates** — each an intelligence agent that must clear before CAR-T can proceed
2. **Portrait/vertical orientation** — the decision flows top-to-bottom like a clinical algorithm, reinforcing the sequential nature
3. **A real patient in crisis** — Ethan, 12, has exhausted all options; CAR-T is his last line of defense
4. **Every gate shows GO** — the emerald checkmarks cascading downward create a powerful visual of clearance
5. **Each gate contributes specific, critical assessments** — tumor characterization, target expression, cardiac function, autoimmune risk, drug metabolism, trial matching, backup therapeutics
6. **The escape scenario is pre-computed** — Gate 7 (Therapeutic Discovery) already has backup compounds ready if CD19 escape occurs
7. **Under 5 minutes** for a decision that traditionally takes days of multidisciplinary review
8. **The urgency is palpable** — rising MRD, exhausted therapies, manufacturing delay ahead
9. **ICU monitoring and supportive care are pre-planned** — this is not just eligibility but a complete treatment plan
10. **Everything runs on one NVIDIA DGX Spark** — anchored by the infrastructure bar at the bottom

The overall impression should be: a rigorous, gate-by-gate eligibility assessment for the most consequential therapy decision in pediatric oncology — seven intelligence agents sequentially clearing a child for CAR-T cell therapy, with backup strategies pre-computed, all on a single desktop supercomputer.
