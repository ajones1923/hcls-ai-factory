# Nano Banana Pro — From Patient DNA to New Medicines v2.0 (Updated for 11 Agents & Pediatric Oncology)

## IMPORTANT: Read this entire prompt before generating. This is an UPDATE to an existing storytelling infographic (v1.0). The original was a narrative-driven, general-audience-friendly diagram showing the journey from patient DNA to drug candidates in 3 stages. v2.0 preserves that storytelling approach but updates for: 3 renamed engines, 11 intelligence agents (was generic), pediatric oncology as primary use case (was generic genomics), cross-agent coordination, and pediatric safety filters. ONE canvas, landscape 16:9. This is NOT a technical architecture poster — it's a STORY told visually for a general audience, clinician, or executive.

---

## REFERENCE: What v1.0 looked like (preserve this style)

The original v1.0 had these characteristics that MUST be maintained:
- Dark navy (#1B2333) top banner with title in white
- Left-to-right narrative flow across 3 stages
- "It Starts with DNA" opening section on the far left with DNA illustration
- Stage 1: "Reading the DNA" — Genomic Foundation Engine
- Stage 2: "Finding What Matters" — funnel visualization narrowing 11.7M variants down to actionable targets
- Stage 3: "Designing New Medicines" — molecule generation to ranked candidates
- Bottom strip with comparison boxes ("AI-Designed vs Original Drug", "A Supercomputer on Your Desk", "Part of a Larger Vision")
- Green footer bar with key metrics
- Warm, approachable, educational tone — NOT intimidating technical jargon
- Small illustrative icons throughout (DNA helix, funnel, molecules, computer)
- The feeling of a magazine infographic or museum exhibit panel

## WHAT CHANGES FROM v1.0 TO v2.0:

1. **Subtitle** → Updated to mention pediatric oncology
2. **"It Starts with DNA"** → Now about a child with cancer (Evelyn, 8, B-ALL)
3. **Stage 1** → "Genomic Foundation Engine" (was "Reading the DNA / GPU-Accelerated Genomics")
4. **Stage 2** → "Precision Intelligence Network" with 11 agents visible (was generic funnel)
5. **Stage 3** → "Therapeutic Discovery Engine" with pediatric safety filter (was generic drug discovery)
6. **Bottom strip** → Updated: "11 Agents Working Together" replaces "AI vs Original Drug" | "6 Pediatric Demos" added | "Part of a Larger Vision" updated with 11 agents
7. **Footer metrics** → Updated numbers
8. **VCP/FTD demo reference** → Replaced with pediatric oncology

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. STORYTELLING infographic for general audiences. Warm, educational, approachable.

**Canvas:** White/light background for the 3 stage areas. Dark navy (#1B2333) top banner and section headers. The diagram should feel like it belongs in a children's hospital lobby display or a TED talk slide — informative but not intimidating.

**CRITICAL TEXT RENDERING RULES:**
- All text must be LEGIBLE at 1920x1080
- Use SHORT PHRASES and LARGE NUMBERS as visual anchors
- Prefer icons + labels over paragraphs
- The narrative should be followable in 60 seconds of scanning left to right

**Typography:**
- Title: Large, bold, white on navy, Inter/Helvetica (24pt+)
- Subtitle: Teal (#1AAFCC) or white, 14pt
- Stage headers: Bold, navy or white on colored background, 16pt
- Body labels: Dark gray (#333333), 10-11pt — SHORT PHRASES ONLY
- Big numbers: Very large (28pt+), NVIDIA green (#76B900) — these are the visual anchors
- Metric badges: White text on green pills
- Footer: White text on green bar

**Color Palette:**
- Navy: #1B2333 — top banner, headers, dark backgrounds
- NVIDIA Green: #76B900 — metric numbers, footer bar, accent badges
- Teal: #1AAFCC — subtitle, agent connections, Stage 2 accents
- White: #FFFFFF — text on dark backgrounds, light areas
- Light Gray: #F5F5F5 — stage area backgrounds
- Purple: #7B2D8E — Stage 3 / drug discovery accents
- Amber: #F5A623 — alert indicators
- Emerald: #059669 — safety/verified badges
- Warm Gray: #666666 — secondary text

---

## CANVAS STRUCTURE

### ━━━ TOP BANNER (full width, navy #1B2333 background) ━━━

**Left corner:** Small green badge: "Open Source • Apache 2.0" + smaller: "Powered by NVIDIA DGX Spark"

**Center (large, commanding, white text):**
- Line 1 (28pt, bold, white): **From Patient DNA to New Medicines**
- Line 2 (14pt, white): How AI Transforms a Child's Genetic Data into Drug Candidates in Under 5 Hours
- Line 3 (11pt, teal #1AAFCC): The HCLS AI Factory on NVIDIA DGX Spark — 3 Engines, 11 Intelligence Agents

**Right corner:** Small badge: "Pediatric Oncology Primary Use Case" in teal on white + "v2.0 • March 2026"

---

### ━━━ MAIN CONTENT AREA (3 stages left to right, ~65% of canvas) ━━━

#### FAR LEFT: "It Starts with DNA" (intro panel)

**Visual:** Stylized DNA double helix illustration, warm and approachable (not clinical). Small child silhouette icon.

**Text (stacked, large readable text):**
```
It Starts with DNA

Meet Evelyn. She's 8 years old.
She has leukemia.

Her DNA contains 3.1 billion letters.
When those letters change, it can
cause disease — or reveal how to
treat it.

[DNA helix icon]
3.1 Billion Letters
A, T, C, G
23 Chromosome Pairs
200 GB of Data
```

**Arrow labeled "Step 1: Read the DNA →"**

---

#### STAGE 1: GENOMIC FOUNDATION ENGINE (left-center panel)

**Header bar (green #76B900 background, white text):**
"Stage 1: Reading the DNA"
"Genomic Foundation Engine — NVIDIA Parabricks"

**Visual: Process flow with icons and BIG NUMBERS:**

```
[Sequencer icon]         [GPU chip icon]         [Document icon]
DNA Sequencer    →    GPU Alignment     →    Variant Calls
FASTQ files           BWA-MEM2               VCF Output
                      Parabricks

         [45 MINUTES]              [>99.7%]
         vs 12 hours CPU           accuracy
```

**Funnel/narrowing visual below:**
```
3.1 Billion bases
    ↓
11.7 Million variants found
    ↓
Annotated with ClinVar (4.1M records)
& AlphaMissense (71M predictions)
```

**Big number callout (NVIDIA green, very large):**
**11.7M** Variants Identified

**Arrow labeled "Step 2: Find What Matters →"**

---

#### STAGE 2: PRECISION INTELLIGENCE NETWORK (center panel — LARGEST)

**Header bar (teal #1AAFCC background, white text):**
"Stage 2: Finding What Matters"
"Precision Intelligence Network — 11 AI Agents"

**Visual: Updated funnel narrowing variants to actionable targets, but NOW showing the 11 agents as the intelligence layer.**

**Top of funnel:**
```
11.7 Million variants enter
         ↓
```

**Middle: The Intelligence Layer (the new addition)**

Visual: A rounded box showing 11 agent icons arranged in a grid or arc, each with a short label:

```
┌─────────────────────────────────────────────────────┐
│              11 Intelligence Agents                  │
│                                                      │
│  [target] Oncology    [heart] Cardiology   [DNA] Rare Disease    │
│  [cell] CAR-T        [brain] Neurology    [pill] Pharmacogenomics│
│  [tube] Biomarker    [shield] Autoimmune  [scan] Imaging         │
│  [cluster] Single-Cell  [clipboard] Clinical Trial               │
│                                                      │
│  Milvus: 3.56M Vectors | Claude LLM | RAG Evidence  │
│  Cross-Agent Coordination • /integrated-assessment   │
└─────────────────────────────────────────────────────┘
```

**Below agents — funnel narrows:**
```
What are the most important variants?
    ↓
Clinker Knowledge Layer
Variant → Gene → Protein → Pathway → Disease → Drug
201 Genes | 13 Therapeutic Areas
    ↓
For Evelyn:
• ETV6-RUNX1 fusion (favorable)
• IKZF1 deletion (unfavorable)
• IKZF1plus → High Risk reclassification
• TPMT *1/*3A → 50% dose reduction needed
• CD19+ → CAR-T eligible
    ↓
```

**Big number callout (teal, very large):**
From **11.7M** → **5 Actionable Findings**

**Arrow labeled "Step 3: Design New Medicines →"**

---

#### STAGE 3: THERAPEUTIC DISCOVERY ENGINE (right panel)

**Header bar (purple #7B2D8E background, white text):**
"Stage 3: Designing New Medicines"
"Therapeutic Discovery Engine — From Targets to Candidates in Minutes"

**Visual: Molecule generation flow with NEW pediatric safety filter:**

```
[target icon]              [molecule icon]          [shield icon]           [ranked list icon]
Target                →    MolMIM              →    Pediatric          →    Ranked
Selection                  Generation               Safety Filter           Candidates

Evelyn's CREBBP            500 candidates           BBB penetration
HAT domain                 generated by             hERG cardiac
                           BioNeMo                  Hepatic immaturity
                                                    Growth plate safety

                           [DiffDock]               500 → 89 pass
                           Binding scored            → Top 3 ranked
```

**Big number callout (purple, very large):**
**100** Ranked Drug Candidates (with **Pediatric Safety** badge in emerald)

**Output icons (right edge):**
```
[document] Reports (PDF, FHIR, Phenopacket)
[molecule] SDF Molecules
[list] Ranked Candidates with Safety Scores
```

---

### ━━━ BOTTOM STRIP: THREE INFO BOXES (matching v1.0 layout) ━━━

Three equal-width boxes side by side below the 3 stages. Light gray backgrounds, thin borders.

#### Box 1 (left): "11 Agents Working Together" [NEW — replaces "AI vs Original Drug"]

```
11 Agents Working Together

For Evelyn, the platform coordinates:
• Oncology identifies ETV6-RUNX1 + IKZF1plus
• Cardiology prevents anthracycline toxicity
• PGx adjusts 6-MP dose (TPMT)
• CAR-T evaluates CD19 eligibility
• Clinical Trial matches COG protocols

[icon] Traditional tumor board: 1-2 weeks
[icon] HCLS AI Factory: 30 seconds
```

#### Box 2 (center): "A Supercomputer on Your Desk" [UPDATED from v1.0]

```
A Supercomputer on Your Desk

[DGX Spark image/icon]

GB10 Grace Blackwell Superchip
ARM CPU + Blackwell GPU
128 GB Unified Memory
NVMe Storage

Everything you just saw runs on this
single $4,699 machine.

3 Engines • 11 Agents • 6 Demos
All on one desktop.
```

#### Box 3 (right): "Part of a Larger Vision" [UPDATED]

```
Part of a Larger Vision

[icon] Genomic Foundation Engine
       Parabricks • DeepVariant

[icon] Precision Intelligence Network
       11 Agents • Milvus • Claude

[icon] Therapeutic Discovery Engine
       MolMIM • DiffDock • Pediatric Safety

6 Pediatric Demos:
"DNA to Drug" • "30-Second Tumor Board"
"Protecting the Survivor" • "One Gene, One Family"
"Last Line of Defense" • "When the Standard Drug Can't Be Used"

Scaling: DGX Spark → DGX B200 → DGX SuperPOD
```

---

### ━━━ FOOTER BAR (full width, NVIDIA green #76B900 background) ━━━

**White bold text, evenly spaced metrics:**

```
HCLS AI Factory    |    11.7M Variants    |    201 Genes    |    11 Agents    |    100 Drug Candidates    |    Under 5 Hours    |    Apache 2.0
```

**Below (very small, gray on white):**
"Created by Adam Jones | hcls-ai-factory.org | Designed with Pediatric Oncology as Primary Use Case | March 2026"

---

## WHAT THIS v2.0 MUST COMMUNICATE vs v1.0

| v1.0 | v2.0 |
|------|------|
| Generic genomics use case | Evelyn, 8 years old, leukemia — a real patient story |
| "GPU-Accelerated Genomics" | "Genomic Foundation Engine" |
| Generic funnel in Stage 2 | 11 named agents visible in the funnel |
| "Generative Drug Discovery" | "Therapeutic Discovery Engine" with Pediatric Safety Filter |
| "AI vs Original Drug" bottom box | "11 Agents Working Together" with tumor board comparison |
| Generic "Part of a Larger Vision" | Named 6 demos + scaling path |
| No pediatric mention | Pediatric oncology throughout — title, patient, findings, safety |
| VCP/FTD demo reference | Evelyn's specific findings (ETV6-RUNX1, IKZF1, TPMT, CD19) |

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **A child's DNA goes in on the left, drug candidates come out on the right** — the story is linear, simple, and human
2. **Three stages are clear and labeled** — Reading (green), Finding (teal), Designing (purple)
3. **11 agents are VISIBLE in Stage 2** — not hidden, each one named with an icon
4. **The numbers tell the story** — 3.1B bases → 11.7M variants → 5 actionable findings → 100 candidates → Top 3 with pediatric safety
5. **Evelyn is real** — her name, her age, her specific findings appear throughout
6. **Pediatric safety is explicit** — BBB, hERG, hepatic, growth plate shown in Stage 3
7. **It runs on one $4,699 machine** — DGX Spark box in bottom strip
8. **It's faster** — "Under 5 Hours" vs traditional months
9. **It's free** — Apache 2.0 in title banner and footer
10. **Someone who knows nothing about genomics can follow the story from left to right in 60 seconds**

The overall impression: a child has cancer, her DNA is sequenced, AI finds what matters and designs new medicines to help her — all on one machine, in under 5 hours, and it's free for every hospital on Earth. That story, told beautifully on one canvas.
