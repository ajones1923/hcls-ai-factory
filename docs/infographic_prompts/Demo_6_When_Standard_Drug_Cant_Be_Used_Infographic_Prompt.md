# Napkin AI Pro — When the Standard Drug Can't Be Used: Pediatric Drug Contraindication and De Novo Molecule Discovery

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram is split into two halves — THE PROBLEM (left, red-tinted) showing why the standard drug vismodegib cannot be used in a growing child, and THE SOLUTION (right, green-tinted) showing how MolMIM generates novel alternatives that avoid the contraindication. The split-screen format makes the before/after contrast viscerally clear. All running on a single NVIDIA DGX Spark.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF) with a clear vertical split — left half has a subtle red tint (#FEF2F2), right half has a subtle green tint (#F0FDF4). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries.

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
- Amber/Orange: #F5A623 — caution indicators, conditional findings
- Red: #DC2626 — THE PROBLEM side accent, contraindication indicators, danger signals
- Emerald Green: #059669 — THE SOLUTION side accent, safety-pass indicators, success signals
- Medium Gray: #666666 — metadata text, secondary labels

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons matching the the platform white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-agent communication, bold red (#DC2626) for contraindication paths
- Split-screen divider: bold vertical line (navy #1B2333) down the center of the canvas
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer

---

## CANVAS STRUCTURE (Top to Bottom, with left/right split)

### ━━━ BAND 1: TITLE BAR (top of canvas, full width) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "Demo 6 of 6" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "When the Standard Drug Can't Be Used"
- **Subtitle line 1 (medium, gray #666666, centered):** "Pediatric Drug Contraindication and De Novo Molecule Discovery"
- **Subtitle line 2 (smaller, gray #666666, centered):** "NVIDIA DGX Spark | GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "March 2026 | HCLS AI Factory"

**Right side — Patient Card** (prominent, white background, thin red #DC2626 left-border accent, rounded corners):
```
PATIENT
——————————————————
Name: Aiden (de-identified)
Age: 10 years old
Dx: SHH-subtype Medulloblastoma
Standard Drug: Vismodegib (Hedgehog pathway inhibitor)
Problem: Vismodegib causes irreversible growth plate fusion
Status: Growth plates still open — drug is contraindicated
```
Red pill badge: "Drug Contraindicated"
Amber pill badge: "Growth Plates Open — Active Growth"

---

### ━━━ LEFT HALF: THE PROBLEM (red-tinted, #FEF2F2 background) ━━━

**Header bar spanning left half (red #DC2626 background, white bold text):**
"THE PROBLEM — Why Vismodegib Cannot Be Used"

#### ── Section A: The Drug ──

**Card: Vismodegib (Erivedge)** [pill icon]
- Large card, red (#DC2626) left-border accent, white interior
- Header: "Vismodegib (Erivedge)" in bold navy
- Sub-header: "FDA-Approved Hedgehog Pathway Inhibitor" in red
- Body:
  - Standard of care for SHH medulloblastoma in adults
  - Inhibits Smoothened (SMO) receptor
  - Blocks Hedgehog signaling pathway
  - Effective: tumor response rate ~35-40%
- Emerald pill: "Effective against SHH pathway"
- Red pill: "BUT: Devastating pediatric side effect"

#### ── Section B: The Contraindication ──

**Card: Growth Plate Fusion** [bone icon with X mark]
- Large card, red (#DC2626) border, light red tint (#FEE2E2)
- Header: "Irreversible Growth Plate Fusion" in bold red
- Body:
  - Hedgehog pathway is essential for bone growth
  - Vismodegib permanently fuses growth plates in children
  - Growth arrest: child stops growing permanently
  - Documented in clinical trials and case reports
  - Effect is IRREVERSIBLE — cannot be undone
- Bold annotation: "A 10-year-old given vismodegib may never grow past 4'8\""

**Visual: Growth plate diagram** (simplified, thin-line illustration):
- Two bone cross-sections side by side
- Left: "Normal" — growth plate open, labeled "Active growth zone" in emerald
- Right: "After Vismodegib" — growth plate fused shut, labeled "PERMANENTLY FUSED" in red
- Red X overlay on the right diagram

#### ── Section C: The Clinical Dilemma ──

**Card: The Dilemma** [scales/balance icon]
- Red-bordered rounded rectangle
- Header: "The Clinical Dilemma" in bold navy
- Body:
  - Aiden has SHH medulloblastoma — the Hedgehog pathway drives his tumor
  - The only FDA-approved targeted therapy (vismodegib) will stunt his growth forever
  - Alternative: cytotoxic chemotherapy + radiation (significant toxicity, neurocognitive effects)
  - What if we could find a molecule that blocks the tumor but spares the growth plates?
- Red pill: "Standard drug: CONTRAINDICATED"
- Red pill: "Alternative (chemo/radiation): HIGH TOXICITY"
- Amber pill: "Unmet clinical need: pediatric-safe Hedgehog inhibitor"

**Bold red arrow pointing rightward across the divider to THE SOLUTION side, labeled:**
"Can AI find an alternative?"

---

### ━━━ CENTER DIVIDER ━━━

**Bold vertical line (navy #1B2333, 3px width) running from below the title bar to above the infrastructure bar.**

**Center medallion (circle on the divider line):** Navy background, white text: "vs"

---

### ━━━ RIGHT HALF: THE SOLUTION (green-tinted, #F0FDF4 background) ━━━

**Header bar spanning right half (emerald #059669 background, white bold text):**
"THE SOLUTION — De Novo Molecule Discovery with MolMIM"

#### ── Section D: Molecule Generation ──

**Card: MolMIM Generation** [molecule-plus icon]
- Large card, NVIDIA green (#76B900) left-border accent, white interior
- Header: "MolMIM — De Novo Molecule Generation" in bold navy
- Sub-header: "BioNeMo Foundation Model" in teal
- Body:
  - Target: Smoothened (SMO) receptor — same as vismodegib
  - Constraint: Must NOT cross-react with growth plate Hedgehog signaling
  - Method: Scaffold-constrained generation with selectivity filters
  - Input: Vismodegib structure as starting scaffold
  - Output: Structurally diverse analogs with modified binding profile
- Metric badge (green pill): "500 candidates generated"
- Metric badge (teal pill): "< 45 minutes generation time"

**Solid gray arrow downward, labeled "500 Raw Candidates"**

#### ── Section E: Multi-Stage Safety Filtering ──

**Funnel visualization** (wide at top, narrow at bottom, showing progressive filtering):

**Stage 1: Drug-Likeness Filter** [funnel icon]
- Green-bordered card
- Body: Lipinski's Rule of Five | Molecular weight | LogP | H-bond donors/acceptors
- Input: 500 candidates
- Output: 312 pass
- Metric badge: "312 / 500 pass"

**Arrow downward**

**Stage 2: Selectivity Filter (THE KEY FILTER)** [target icon]
- Emerald-bordered card, slightly larger, prominent
- Header: "Growth Plate Selectivity Filter" in bold emerald
- Body: Screens for selectivity between tumor SMO and growth plate SMO | Differential binding affinity modeling | Must show >10x selectivity for tumor context
- This is the filter that separates this pipeline from standard drug discovery
- Input: 312 candidates
- Output: 142 pass selectivity
- Metric badge (emerald pill): "142 / 312 pass selectivity"
- Bold annotation: "THIS IS THE CRITICAL FILTER — spares growth plates"

**Arrow downward**

**Stage 3: hERG Cardiac Liability Screen** [heart icon]
- Green-bordered card
- Body: hERG potassium channel inhibition | Must be > 10 uM (low cardiac risk) | Pediatric cardiac safety threshold
- Input: 142 candidates
- Output: 104 pass
- Metric badge: "104 / 142 pass hERG"

**Arrow downward**

**Stage 4: Pediatric Toxicity Panel** [shield icon]
- Emerald-bordered card
- Body: CYP450 interactions | Developmental toxicity | Blood-brain barrier permeability (required — brain tumor) | Neurotoxicity screen
- Input: 104 candidates
- Output: 89 pass
- Metric badge (emerald pill): "89 / 104 pass pediatric panel"

**Arrow downward**

#### ── Section F: DiffDock Binding + Ranking ──

**Card: DiffDock Scoring** [lock-and-key icon]
- NVIDIA green-bordered card
- Header: "DiffDock Binding Prediction" in bold navy
- Body: Diffusion-based molecular docking | Binding affinity scoring against SMO | Pose generation and ranking
- Input: 89 safety-cleared candidates
- Metric badge: "89 candidates docked"

**Arrow downward, labeled "Ranked by Binding Affinity + Selectivity"**

#### ── Section G: Top 3 Candidates ──

**Card: Top 3 Candidates** [trophy icon]
- Larger card, prominent, double-bordered (NVIDIA green outer, emerald inner), light green background (#D1FAE5)
- Header: "Top 3 Pediatric-Safe SHH Inhibitors" in bold navy

**Three sub-cards arranged vertically:**

**Candidate 1 (best):**
- Green pill: "Binding: -9.3 kcal/mol"
- Green pill: "IC50: 8 nM against SHH"
- Emerald pill: "Growth plate selectivity: 47x"
- Emerald pill: "hERG: 18.2 uM — SAFE"
- Emerald pill: "MW: 438 | BBB: Permeable"

**Candidate 2:**
- Green pill: "Binding: -8.9 kcal/mol"
- Green pill: "IC50: 12 nM against SHH"
- Emerald pill: "Growth plate selectivity: 32x"
- Emerald pill: "hERG: 15.7 uM — SAFE"
- Emerald pill: "MW: 411 | BBB: Permeable"

**Candidate 3:**
- Green pill: "Binding: -8.5 kcal/mol"
- Green pill: "IC50: 18 nM against SHH"
- Emerald pill: "Growth plate selectivity: 28x"
- Emerald pill: "hERG: 22.1 uM — SAFE"
- Emerald pill: "MW: 465 | BBB: Permeable"

**Summary annotation below top 3:**
- "All 3 candidates: IC50 8-18 nM | Growth-plate selective | BBB-permeable | Pediatric-safe"
- "Next step: In vitro validation in SHH medulloblastoma cell lines"

---

### ━━━ BAND BOTTOM: COMPARISON SUMMARY (full width, below the split) ━━━

**Position:** Full-width horizontal strip below the split-screen section. Navy (#1B2333) background. White text.

**Three-column comparison:**

| Vismodegib (Standard) | Chemo + Radiation (Alternative) | MolMIM Candidates (AI-Generated) |
|---|---|---|
| Effective but CONTRAINDICATED | Effective but HIGH TOXICITY | Effective AND Pediatric-SAFE |
| Growth plate fusion: PERMANENT | Neurocognitive damage | Growth plate selectivity: 28-47x |
| Not usable in children | Long-term morbidity | IC50: 8-18 nM (potent) |
| Single mechanism | Broad cytotoxicity | BBB-permeable (brain tumor) |

**Center annotation (bold, white):** "AI-generated molecules solve a problem that no existing drug can address for growing children."

---

### ━━━ INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | BioNeMo NIMs | Cheminformatics | Safety Filtering | Open Source |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | MolMIM molecule generation | RDKit property calculation | Pediatric-specific thresholds | Apache 2.0 licensed |
| Desktop AI supercomputer | DiffDock molecular docking | Lipinski / Veber filters | Growth plate selectivity screen | ~1.1 TB total footprint |
| 128 GB unified LPDDR5x | Scaffold-constrained generation | BBB permeability prediction | hERG cardiac liability | Fully reproducible |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "500 candidates" on MolMIM generation
- "89 pass all safety filters" on the funnel output
- "IC50 8-18 nM" on top 3 candidates
- "28-47x growth plate selectivity" on selectivity filter
- "hERG > 10 uM" on all candidates
- "BBB-permeable" on all candidates (required for brain tumor)
- "< 2 hours total pipeline" on the solution side

**The critical insight annotation** (small dashed box, center of canvas near the divider):
- "WHY THIS MATTERS"
- Vismodegib works but permanently stunts growth in children
- No existing drug solves this problem
- AI-generated molecules can be designed with selectivity constraints
- The Hedgehog pathway has different roles in tumors vs growth plates
- Selective inhibition is theoretically possible — MolMIM finds the molecules

**Filtering funnel annotation** (small dashed box, alongside the funnel):
- "Progressive Safety Filtering"
- 500 → 312 → 142 → 104 → 89 candidates
- Each stage eliminates a different risk category
- The selectivity filter (Stage 2) is unique to this use case
- 89 candidates in the final pool = high hit rate (17.8%)

**Data flow arrows style guide:**
- **Solid arrows:** Primary data flow (generation through filtering to ranking). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-reference between problem and solution sides. Teal (#1AAFCC) dashes.
- **Bold arrows:** Contraindication warning paths (left side) and the critical "Can AI find an alternative?" bridge arrow. Red (#DC2626).

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Split-screen format** — THE PROBLEM (red, left) vs THE SOLUTION (green, right) makes the contrast immediate
2. **A real clinical dilemma** — the standard drug works but will permanently stunt a 10-year-old's growth
3. **Growth plate fusion is irreversible** — the bone diagram makes this viscerally clear
4. **MolMIM generates 500 novel candidates** as alternatives to the contraindicated drug
5. **The filtering funnel is rigorous** — 500 down to 89, with the growth plate selectivity filter being the critical innovation
6. **Top 3 candidates are potent AND safe** — IC50 8-18 nM, growth plate selectivity 28-47x, all BBB-permeable
7. **BBB permeability is required** — this is a brain tumor, so the drug must cross the blood-brain barrier
8. **The three-column comparison at the bottom** makes it clear that AI-generated molecules solve what no existing drug can
9. **Under 2 hours** for the entire generation-to-ranking pipeline
10. **Everything runs on one NVIDIA DGX Spark** — anchored by the infrastructure bar at the bottom

The overall impression should be: a clear, compelling visual argument that AI-driven molecule generation can solve pediatric drug contraindications that have no existing solution — transforming an impossible clinical dilemma into a ranked list of safe, potent candidates on a single desktop supercomputer.
