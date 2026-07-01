# Napkin AI Pro — One Gene, One Family: Hereditary Cancer and Lifetime Surveillance

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram traces a single RB1 mutation through a family — from a child's bilateral retinoblastoma diagnosis, through cascade genetic testing of parents and siblings, to five intelligence agents generating personalized surveillance plans, and finally a lifetime monitoring timeline from age 0 to 50. All running on a single NVIDIA DGX Spark.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of the platform's white paper series (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a genetic counseling reference poster — a single sheet that captures the entire hereditary cancer story for one family.

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
- Amber/Orange: #F5A623 — carrier status indicators, watch-level findings
- Red: #DC2626 — affected status indicators, critical findings, mutation-positive markers
- Emerald Green: #059669 — mutation-negative indicators, cleared status
- Medium Gray: #666666 — metadata text, secondary labels

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons matching the the platform white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-agent communication, bold red (#DC2626) for critical mutation inheritance paths
- Family tree connections: solid lines for parent-child, dashed for sibling
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer

---

## CANVAS STRUCTURE (Top to Bottom, 5 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "HCLS AI Factory" — plus a second smaller badge below it: "Demo 4 of 6" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "One Gene, One Family"
- **Subtitle line 1 (medium, gray #666666, centered):** "RB1 Hereditary Retinoblastoma — From Diagnosis to Lifetime Surveillance"
- **Subtitle line 2 (smaller, gray #666666, centered):** "NVIDIA DGX Spark | GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "March 2026 | HCLS AI Factory"

**Right side — Patient Card** (prominent, white background, thin red #DC2626 left-border accent, rounded corners):
```
INDEX PATIENT
——————————————————
Name: Aurora (de-identified)
Age: 4 years old
Dx: Bilateral retinoblastoma
Gene: RB1 c.958C>T (p.Arg320Ter)
Classification: Pathogenic (ACMG Class 5)
Inheritance: Autosomal dominant, 50% risk to offspring
```
Small thin-line icon of a family silhouette. Red pill badge: "Hereditary Cancer Syndrome"

---

### ━━━ BAND 2: FAMILY PEDIGREE — RB1 MUTATION CASCADE (full width) ━━━

**Position:** Full-width horizontal band below the title bar. Light gray background (#F5F5F5). This band shows the family tree with genetic test results cascading through each member.

**Header bar (navy #1B2333 background, white text):** "Family Pedigree — RB1 Cascade Genetic Testing"

**Family tree layout (centered, classic pedigree format):**

**Generation 1 (top row) — Parents:**

**Father Card** (left, rounded rectangle, thin border):
- [male silhouette icon]
- Header: "Father" in bold navy
- Body: Age 34 | Healthy | No eye history
- Genetic test result: Large emerald (#059669) bordered box
- Emerald pill: "RB1: WILD-TYPE (Negative)"
- Emerald pill: "No carrier status"
- Annotation: "De novo mutation confirmed in Aurora"

**Mother Card** (right, rounded rectangle, thin border):
- [female silhouette icon]
- Header: "Mother" in bold navy
- Body: Age 32 | Healthy | No eye history
- Genetic test result: Large emerald (#059669) bordered box
- Emerald pill: "RB1: WILD-TYPE (Negative)"
- Emerald pill: "No carrier status"

**Vertical lines descend from parents to children, meeting at a horizontal bar**

**Generation 2 (bottom row) — Children:**

**Aurora Card** (left child, larger, prominent, red #DC2626 left-border accent):
- [female child silhouette icon]
- Header: "Aurora — Index Patient" in bold navy
- Body: Age 4 | Bilateral retinoblastoma
- Genetic test result: Large red (#DC2626) bordered box
- Red pill: "RB1 c.958C>T (p.Arg320Ter)"
- Red pill: "PATHOGENIC — ACMG Class 5"
- Red pill: "De novo germline mutation"
- Metric badge (red pill): "50% transmission risk to future offspring"

**Sibling Card** (right child, rounded rectangle):
- [male child silhouette icon]
- Header: "Sibling — Leo" in bold navy
- Body: Age 7 | Healthy
- Genetic test result: Large emerald (#059669) bordered box
- Emerald pill: "RB1: WILD-TYPE (Negative)"
- Emerald pill: "Standard population risk"
- Annotation: "No enhanced surveillance needed"

**Cascade arrows** (bold red #DC2626 from Aurora, dashed teal #1AAFCC to family members):
- Bold red downward arrow from Aurora labeled: "De novo pathogenic variant identified"
- Dashed teal arrows to Father and Mother labeled: "Cascade testing triggered"
- Dashed teal arrow to Sibling labeled: "Targeted testing: RB1 c.958C>T"

**Annotation box (dashed border, bottom of pedigree):**
- "De novo mutation: Neither parent carries the variant"
- "Aurora's children will have 50% risk — prenatal testing recommended"
- "Sibling cleared — no enhanced surveillance required"
- "Mosaicism in parents cannot be fully excluded — recommend periodic exam"

---

### ━━━ BAND 3: FIVE INTELLIGENCE AGENTS (full width, arranged as flowing cards) ━━━

**Position:** Full-width horizontal band below the pedigree. Light blue background tint (#E3F2FD).

**Header bar (navy #1B2333 background, white text):** "Intelligence Agent Network — Personalized for Aurora"

**Five agent cards arranged left-to-right, each receiving a dashed teal arrow from Aurora's pedigree card above. Each card flows output downward to Band 4.**

**Agent Card 1: Rare Disease ACMG Classification Agent** [document-check icon]
- Left-border accent: NVIDIA green (#76B900)
- Light green tint (#E8F5E9)
- Header: "Rare Disease ACMG Agent" in bold navy
- Sub-header: "Variant Classification" in teal
- **Findings:**
  - Red pill: "RB1 c.958C>T — Pathogenic (Class 5)"
  - Body: Nonsense mutation | Premature stop codon | Loss of function
  - Amber pill: "ClinVar: 47 submissions, all Pathogenic"
  - Green pill: "AlphaMissense: 0.97 pathogenicity score"
  - Body: ACMG criteria met: PVS1 + PS3 + PM2 + PP5
- Metric badge (green pill): "5 ACMG criteria satisfied"

**Agent Card 2: Precision Oncology Agent** [radar icon]
- Left-border accent: Red (#DC2626)
- Light red tint (#FEE2E2)
- Header: "Precision Oncology Agent" in bold navy
- Sub-header: "Cancer Risk Profiling" in teal
- **Findings:**
  - Red pill: "Bilateral retinoblastoma confirmed — germline RB1"
  - Red pill: "Lifetime second malignancy risk: osteosarcoma (highest)"
  - Amber pill: "Soft tissue sarcoma risk: elevated"
  - Amber pill: "Melanoma risk: moderately elevated"
  - Amber pill: "Pinealoblastoma (trilateral Rb): MRI surveillance"
- Metric badge (red pill): "5 second-cancer risks identified"

**Agent Card 3: Imaging Intelligence Agent** [eye icon]
- Left-border accent: Teal (#1AAFCC)
- Light blue tint (#E3F2FD)
- Header: "Imaging Intelligence Agent" in bold navy
- Sub-header: "Ocular + Systemic Surveillance" in teal
- **Findings:**
  - Red pill: "EUA (exam under anesthesia) every 3-4 weeks until age 3"
  - Amber pill: "Brain MRI every 6 months until age 5 (trilateral Rb screen)"
  - Amber pill: "Annual whole-body MRI from age 8"
  - Emerald pill: "Current imaging: bilateral tumors responding to chemoreduction"
- Metric badge (teal pill): "3 imaging protocols scheduled"

**Agent Card 4: Clinical Trial Intelligence Agent** [clipboard-search icon]
- Left-border accent: Emerald (#059669)
- Light emerald tint (#D1FAE5)
- Header: "Clinical Trial Intelligence Agent" in bold navy
- Sub-header: "Trial Matching" in teal
- **Findings:**
  - Emerald pill: "NCT04300517 — Intra-arterial chemotherapy: ELIGIBLE"
  - Emerald pill: "NCT03284684 — Rb with germline RB1: MATCH"
  - Amber pill: "NCT05199584 — Second malignancy prevention: ELIGIBLE at age 10"
  - Teal pill: "3 trials matched across treatment + prevention"
- Metric badge (emerald pill): "3 clinical trials matched"

**Agent Card 5: Therapeutic Discovery Agent** [molecule icon]
- Left-border accent: NVIDIA green (#76B900)
- Light green tint (#E8F5E9)
- Header: "Therapeutic Discovery Agent" in bold navy
- Sub-header: "Novel Compound Screening" in teal
- **Findings:**
  - Green pill: "MolMIM: 500 candidates targeting RB1-pathway"
  - Green pill: "MDM2 inhibitors: 14 candidates generated"
  - Green pill: "Top candidate: -8.8 kcal/mol, MW 412"
  - Emerald pill: "3 pass pediatric safety filter (hERG >10 uM)"
- Metric badge (green pill): "500 candidates screened"

---

### ━━━ BAND 4: LIFETIME SURVEILLANCE TIMELINE (full width) ━━━

**Position:** Full-width horizontal band below the agent cards. White background with a prominent horizontal timeline.

**Header bar (navy #1B2333 background, white text):** "Lifetime Surveillance Timeline — Aurora, RB1 Germline Carrier"

**Horizontal timeline spanning the full width of the band, marked with age milestones from 0 to 50:**

**Timeline bar:** Thick horizontal line, gradient from red (#DC2626) on left (intensive early surveillance) fading to amber (#F5A623) in middle to emerald (#059669) on right (maintained but less intensive).

**Age markers and surveillance events (positioned above and below the timeline):**

**Age 0-3 (red zone):**
- Above: "EUA every 3-4 weeks" [eye icon]
- Below: "Brain MRI every 6 months" [brain icon]
- Badge: "Highest intensity"

**Age 3-5 (red-amber zone):**
- Above: "EUA every 2-3 months" [eye icon]
- Below: "Brain MRI every 6 months → annual" [brain icon]
- Additional: "Chemoreduction monitoring"

**Age 5-10 (amber zone):**
- Above: "EUA every 4-6 months" [eye icon]
- Below: "Annual brain MRI" [brain icon]
- Additional: "Bone health screening begins" [bone icon]

**Age 10-18 (amber zone):**
- Above: "Annual ophthalmology" [eye icon]
- Below: "Annual whole-body MRI" [body scan icon]
- Additional: "Osteosarcoma peak risk period" [warning icon]
- Red annotation: "Peak osteosarcoma risk: ages 10-20"

**Age 18-30 (amber-emerald zone):**
- Above: "Annual ophthalmology" [eye icon]
- Below: "Annual whole-body MRI" [body scan icon]
- Additional: "Genetic counseling for family planning" [DNA icon]
- Amber annotation: "50% transmission risk — prenatal testing available"

**Age 30-50 (emerald zone):**
- Above: "Annual comprehensive screening" [clipboard icon]
- Below: "Soft tissue sarcoma surveillance" [body scan icon]
- Additional: "Melanoma skin checks annually" [skin icon]
- Emerald annotation: "Ongoing but lower intensity"

**Below the timeline, a summary row:**
- "Powered by: 5 intelligence agents continuously updating protocols as new evidence emerges"
- "Every screening recommendation is evidence-grounded with literature citations"

---

### ━━━ BAND 5: NVIDIA DGX SPARK INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | Genomic Analysis | Agent Framework | Vector Intelligence | Open Source |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | Parabricks 4.6 + DeepVariant | 5 specialized intelligence agents | Milvus: 3.56M vectors | Apache 2.0 licensed |
| Desktop AI supercomputer | ACMG variant classification | Cascade testing orchestration | ClinVar + AlphaMissense | ~1.1 TB total footprint |
| 128 GB unified LPDDR5x | Family-aware analysis | Lifetime surveillance generation | BGE-small-en-v1.5 embeddings | Fully reproducible |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "1 gene, 4 lives assessed" on the pedigree
- "ACMG Class 5 — Pathogenic" on the variant
- "50% transmission risk" on Aurora's card
- "5 second-cancer risks" on Oncology agent
- "3 imaging protocols" on Imaging agent
- "3 clinical trials" on Clinical Trial agent
- "500 candidates" on Therapeutic Discovery agent
- "Age 0-50 surveillance" on the timeline

**Genetic inheritance annotation** (small dashed box, right of pedigree):
- "Autosomal Dominant Inheritance"
- One copy of mutant RB1 inherited (or de novo)
- Second hit (somatic) → tumor formation (Knudson's two-hit hypothesis)
- 50% risk to each offspring of a carrier
- Penetrance: ~90% for retinoblastoma

**Data flow arrows style guide:**
- **Solid arrows:** Primary data flow (genome sequencing → variant calling → agent analysis). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-agent communication and cascade testing triggers. Teal (#1AAFCC) dashes.
- **Bold arrows:** Critical mutation inheritance paths and high-risk alerts. Red (#DC2626).

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **One gene affects an entire family** — RB1 mutation in Aurora cascades to testing her parents and sibling
2. **De novo mutation** — neither parent carries the variant, but Aurora's future children have 50% risk
3. **The pedigree is central** — a classic genetic counseling family tree with color-coded test results
4. **Five intelligence agents** personalize surveillance for Aurora's specific variant and cancer risks
5. **Osteosarcoma is the biggest second-cancer risk** — prominently flagged by the Oncology agent
6. **Lifetime surveillance from age 0 to 50** — the timeline makes the duration and intensity of follow-up viscerally clear
7. **Early years are the most intensive** — the red-to-green gradient shows surveillance intensity decreasing over time
8. **Clinical trials are matched** — not just for current treatment, but for future second-malignancy prevention
9. **Therapeutic Discovery is generating novel compounds** — even for a rare disease like hereditary retinoblastoma
10. **Everything runs on one NVIDIA DGX Spark** — from genome sequencing through family cascade testing to lifetime surveillance planning

The overall impression should be: a complete hereditary cancer story told through one family — from a child's diagnosis through cascade genetic testing of every family member to a lifetime surveillance plan generated by five coordinated intelligence agents, all on a single desktop supercomputer.
