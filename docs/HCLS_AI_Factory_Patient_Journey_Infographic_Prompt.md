# Napkin AI Pro — HCLS AI Factory: The Patient Journey
# From DNA to Drug Candidates — A Visual Walkthrough

## IMPORTANT: Read this entire prompt before generating. This describes a single, visually rich infographic designed for NON-TECHNICAL audiences — healthcare executives, hospital administrators, conference attendees, investors, and policymakers. Unlike a technical architecture diagram, this infographic TELLS A STORY. It guides the viewer through a patient's journey from a DNA sample to 100 novel drug candidates, using simple language, visual metaphors, and clear cause-and-effect. Every element appears on ONE canvas. The demo case is VCP/Frontotemporal Dementia (FTD). The entire pipeline runs on a single NVIDIA DGX Spark desktop workstation ($3,999).

---

## OVERALL LAYOUT AND STYLE

Create a visually stunning, story-driven infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a premium healthcare innovation poster — the kind displayed at HIMSS, GTC, or NVIDIA keynote events. Clean, modern, approachable. It should make a non-scientist say "I understand this" while making a scientist say "This is accurate."

**Design Philosophy:** Guide the viewer's eye from left to right (or top to bottom) through a clear narrative journey. Each stage should feel like turning a page in a story. Use visual metaphors that connect unfamiliar scientific concepts to everyday understanding. The infographic should be self-explanatory — someone should be able to understand the entire pipeline without any verbal explanation.

**Canvas:** Clean white background (#FFFFFF) with generous white space between sections. Not dense like a technical diagram — open, breathable, and inviting. Each stage should have enough visual weight to stand alone but clearly flow into the next.

**Typography:**
- Main title: Large, bold, modern sans-serif (Inter, Montserrat, or similar), deep navy (#1B2333)
- Stage headers: Bold, slightly smaller than title, in stage-specific colors
- Explanatory text: Clean, readable sans-serif, dark gray (#333333), written at a high school reading level — no jargon without immediate explanation
- Callout numbers: Extra-bold, oversized, gradient-colored — these are the "wow" numbers audiences remember
- Annotations/labels: Small but clear, medium gray (#666666)

**Color Palette (exact):**
- NVIDIA Green: #76B900 — primary accent for GPU/NVIDIA elements, success indicators, CTA elements
- Deep Navy: #1B2333 — titles, dark sections, authority and trust
- Teal: #00B4D8 — data flow, AI elements, digital/tech feeling
- Warm Orange: #F5A623 — drug discovery elements, results, achievements
- Soft Purple: #8B5CF6 — annotation/database elements, intelligence, knowledge
- Alert Red: #DC2626 — disease/pathogenic indicators (used sparingly for emphasis)
- Emerald: #059669 — positive outcomes, drug-likeness, healthy indicators
- Light backgrounds: Very subtle tints of each stage color for section backgrounds
- White: #FFFFFF — canvas, card interiors, breathing space

**Visual Elements:**
- Large, friendly icons (not technical line icons — more illustrative, slightly rounded)
- Smooth curved arrows and flowing connection lines (not rigid straight arrows)
- Rounded-corner cards (12-16px radius) with subtle drop shadows for depth
- Gradient overlays on key sections for visual richness
- "Before vs. After" comparison callouts to emphasize advantages
- A clear visual "path" or "journey line" connecting all stages

---

## THE STORY ARC

The infographic tells this story in five acts:

**Act 1:** "Every person's DNA is unique" — The patient starting point
**Act 2:** "A supercomputer on your desk reads the DNA" — GPU-powered genomics
**Act 3:** "AI finds the needle in the haystack" — Variant annotation and target discovery
**Act 4:** "AI designs new medicines" — Drug candidate generation
**Act 5:** "Results that used to take months — delivered in hours" — Impact and outcomes

---

## CANVAS STRUCTURE

### ━━━ SECTION 1: TITLE AND CONTEXT BAR (top of canvas) ━━━

**Layout:** Full-width banner across the top. Dark navy (#1B2333) background with a subtle diagonal gradient to slightly lighter navy. Clean and authoritative.

**Left side:**
- Small NVIDIA green (#76B900) rounded badge: "Open Source — Apache 2.0"
- Below it: "Powered by NVIDIA DGX Spark" in white text

**Center (dominant element):**
- **Main title (large, bold, white):** "From Patient DNA to New Medicines"
- **Subtitle (medium, light gray #CCCCCC):** "How AI Transforms Genetic Data into Drug Candidates in Under 5 Hours"
- **Tagline (smaller, NVIDIA green #76B900):** "The HCLS AI Factory on NVIDIA DGX Spark"

**Right side:**
- Author: "Adam Jones | February 2026" in small light gray text
- Small badge: "VCP / Frontotemporal Dementia Demo" in teal (#00B4D8) with white text

---

### ━━━ SECTION 2: THE PATIENT STARTING POINT (left side of canvas) ━━━

**Layout:** A visually warm, human-centered opening section. This is where the story begins — with a person.

**Visual:** A friendly, abstract illustration of a patient (silhouette or simple figure — NOT a photo). Next to the figure, a stylized DNA double helix icon, large and recognizable.

**Header:** "It Starts with DNA" in bold navy (#1B2333)

**Explanatory text (simple language, conversational):**
> "Every cell in your body contains DNA — a molecule with 3.1 billion letters that makes you unique. When something goes wrong in those letters, it can cause disease. The first step to finding a cure is reading the DNA."

**Key fact cards (small rounded cards arranged near the figure):**

Card 1 — DNA icon:
- "3.1 Billion Letters"
- "Your complete genetic instruction manual"

Card 2 — Microscope icon:
- "A, T, C, G"
- "The 4-letter alphabet of life"

Card 3 — Chromosome icon:
- "23 Pairs of Chromosomes"
- "46 packages of DNA in every cell"

Card 4 — Document stack icon:
- "~200 GB of Data"
- "One person's DNA fills 50 HD movies' worth of storage"

**Transition element:** A large, flowing curved arrow (NVIDIA green #76B900) sweeps from the patient figure rightward into Section 3, labeled: "Step 1: Read the DNA"

---

### ━━━ SECTION 3: STAGE 1 — READING THE DNA (first major pipeline section) ━━━

**Background:** Very light green tint (#F0F9E8) with a subtle NVIDIA green left-border accent.

**Header:** "Stage 1: Reading the DNA" in bold NVIDIA green (#76B900)
**Subheader:** "GPU-Accelerated Genomics — What Used to Take Days Now Takes Minutes" in dark gray

**Visual storytelling approach — Three illustrated steps flowing left to right:**

**Step 1.1 — The Sequencing Machine**
- Visual: Stylized illustration of a DNA sequencer (a boxy lab instrument)
- Label: "DNA Sequencer"
- Explanation: "A machine chops the DNA into millions of tiny pieces and reads each one. Each piece is about 250 letters long."
- Output callout: "Billions of short reads → ~200 GB of raw data (FASTQ files)"

**Step 1.2 — Putting the Puzzle Together (Alignment)**
- Visual: Puzzle-piece metaphor — scattered pieces flowing into a completed puzzle
- Label: "BWA-MEM2 Alignment on GPU"
- Explanation: "Like assembling a 3.1-billion-piece jigsaw puzzle — the computer matches each short read to the right position in a reference human genome."
- Speed comparison callout (split card, green vs gray):
  - Left (gray, crossed out): "CPU: 12-24 hours"
  - Right (green, checkmark): "GPU: 20-45 minutes"
  - Bottom: "10-20x faster with NVIDIA Parabricks"

**Step 1.3 — Finding the Differences (Variant Calling)**
- Visual: Magnifying glass over a DNA strand, highlighting a single changed letter (e.g., showing "G" where "A" should be, highlighted in red)
- Label: "Google DeepVariant — AI Variant Calling"
- Explanation: "An AI trained to spot real differences in DNA — the same kind of AI that recognizes faces in photos, but applied to genetics. Over 99% accurate."
- Speed comparison callout:
  - Left (gray, crossed out): "CPU: 8-12 hours"
  - Right (green, checkmark): "GPU: 10-35 minutes"

**Stage 1 output — Large highlighted result card:**
- Big number: "11.7 Million" in bold gradient green text
- Label: "Variants Discovered"
- Sub-label: "Differences found between this patient's DNA and the reference genome"
- Smaller detail: "3.5 million pass quality filters"
- File icon: "Output: VCF File (Variant Call Format)"

**Stage 1 timing badge:** Green rounded pill — "Total: 1-4 Hours (vs. 1-2 Days on CPU)"

**Transition arrow:** Flowing curved arrow from the VCF output downward or rightward into Section 4, labeled: "Step 2: Find What Matters"

---

### ━━━ SECTION 4: STAGE 2 — FINDING WHAT MATTERS (second major pipeline section) ━━━

**Background:** Very light teal tint (#E8F7FC) with a teal (#00B4D8) left-border accent.

**Header:** "Stage 2: Finding What Matters" in bold teal (#00B4D8)
**Subheader:** "AI-Powered Analysis — Finding the Needle in an 11.7-Million-Piece Haystack" in dark gray

**Visual storytelling approach — A funnel metaphor is the centerpiece:**

**The Annotation Funnel (CENTRAL VISUAL ELEMENT of this section)**

Design a large, visually striking funnel graphic. Wide at the top, narrow at the bottom. Each level of the funnel is color-coded and labeled with counts. The funnel should be the most prominent element in this section — it's the "wow" visual that shows how AI narrows millions of variants to actionable targets.

Funnel levels (top to bottom, widest to narrowest):

Level 1 (widest, light gray):
- "11.7 Million Variants" — "Everything found in Stage 1"

Level 2 (slightly narrower, light green):
- "3.5 Million" — "Pass quality filters (high confidence)"

Level 3 (narrower, light purple):
- "35,616" — "Match known clinical variants (ClinVar database)"

Level 4 (narrower still, medium purple):
- "6,831" — "Have AI pathogenicity predictions (AlphaMissense)"

Level 5 (narrow, light orange):
- "~2,400" — "High-impact, disease-causing variants"

Level 6 (narrowest, bold emerald):
- "847" — "In genes we can target with drugs"

**Three knowledge sources feeding INTO the funnel (arranged beside or above the funnel):**

Database Card 1 — ClinVar (purple accent):
- Icon: Hospital/medical database icon
- "ClinVar — NIH Clinical Database"
- "4.1 million variants studied by scientists worldwide"
- "Classifies variants: Pathogenic → Likely Pathogenic → Uncertain → Likely Benign → Benign"
- Simple analogy: "Like a medical encyclopedia for DNA changes"

Database Card 2 — AlphaMissense (purple accent):
- Icon: AI brain/neural network icon
- "AlphaMissense — DeepMind AI"
- "Predicts whether a variant causes disease"
- "Scores 0 to 1 — above 0.564 = likely disease-causing"
- Simple analogy: "From the same team that won the Nobel Prize for protein structure prediction"

Database Card 3 — VEP (purple accent):
- Icon: Function/impact icon
- "VEP — Variant Effect Predictor"
- "Tells you what the variant does to the protein"
- "Impact: HIGH | MODERATE | LOW"
- Simple analogy: "Like a damage assessment report for each DNA change"

**Vector Database explanation (small card below funnel):**
- Icon: Database/search icon in teal
- "Smart Search: Milvus Vector Database"
- "Each variant is converted to a mathematical fingerprint (384 numbers) so AI can search 3.5 million variants instantly using plain English questions"
- Simple analogy: "Like Google Search, but for your DNA"

**AI Reasoning section (below or beside the funnel):**
- Visual: Chat interface mockup showing a question and answer
- Label: "Claude AI — Grounded in Real Evidence"
- Explanation: "A medical AI assistant reads all the evidence and identifies the best drug targets. It can only cite variants that actually exist in the patient's data — no guessing."
- Example query shown in a chat bubble: "What are the most promising drug targets for neurodegenerative disease?"
- Example response summary: "VCP gene identified — Pathogenic variant linked to Frontotemporal Dementia"

**Knowledge Base callout (emerald accent):**
- Three emerald pill badges: "201 Genes" | "13 Therapeutic Areas" | "171 Druggable Targets (85%)"
- Explanation: "The AI draws from a curated knowledge base spanning neurology, oncology, cardiovascular, and 10 more disease categories"

**VCP Discovery callout (PROMINENT — this is the demo "aha" moment):**
- Highlighted card with teal border and subtle glow/shadow
- Header: "Target Found: VCP Gene"
- Icon: Bullseye/target icon
- Details in simple language:
  - "Variant: A single letter change at position 35,065,263 on chromosome 9 (G changed to A)"
  - "ClinVar says: Pathogenic (disease-causing)"
  - "AlphaMissense score: 0.87 out of 1.0 (well above the 0.564 danger threshold)"
  - "Disease: Frontotemporal Dementia (FTD) — a brain disease affecting personality, behavior, and language"
  - "Good news: This protein is 'druggable' — scientists know how to design medicines that target it"

**Transition arrow:** Flowing curved arrow from the VCP target card into Section 5, labeled: "Step 3: Design New Medicines"

---

### ━━━ SECTION 5: STAGE 3 — DESIGNING NEW MEDICINES (third major pipeline section) ━━━

**Background:** Very light orange tint (#FFF8F0) with a warm orange (#F5A623) left-border accent.

**Header:** "Stage 3: Designing New Medicines" in bold warm orange (#F5A623)
**Subheader:** "AI-Driven Drug Discovery — From Target to 100 Candidate Molecules in Minutes" in dark gray

**Visual storytelling approach — Five illustrated steps flowing as a clear pathway:**

**Step 3.1 — Finding the Protein Shape**
- Visual: 3D protein structure illustration (ribbon diagram style, simplified and colorful — not too technical)
- Label: "Protein Structure from PDB"
- Explanation: "Before designing a drug, you need to know what the target protein looks like in 3D. Scientists have mapped the VCP protein's shape using techniques like X-ray crystallography and Cryo-EM (electron microscopy of frozen samples)."
- Detail card: "Best structure: 5FTK — shows VCP with an existing drug (CB-5083) already attached, revealing exactly where new drugs should bind"
- Simple analogy: "Like getting the blueprint of a lock before designing a key"

**Step 3.2 — The Starting Point (Seed Compound)**
- Visual: Simple 2D molecular structure icon (hexagonal ring shapes, friendly style)
- Label: "Seed Compound: CB-5083"
- Explanation: "CB-5083 is an existing drug that was tested in clinical trials. It binds to VCP but isn't perfect. The AI uses it as a starting point to design something better."
- Simple analogy: "Like giving an architect an existing house plan and asking them to design 100 improved versions"

**Step 3.3 — AI Generates New Molecules**
- Visual: One molecule on the left branching into many molecules on the right (tree/branching visual). Stars or sparkle effects to indicate AI generation.
- Label: "BioNeMo MolMIM — NVIDIA AI"
- Explanation: "This AI generates 100 brand-new molecule designs. It's the same type of AI that powers text generators, but instead of writing sentences, it writes molecular structures."
- Callout badge: "100 Novel Candidates Generated"
- Timing: "Under 5 minutes"

**Step 3.4 — Testing If They Fit (Virtual Docking)**
- Visual: A key fitting into a lock — molecule shape fitting into a protein pocket
- Label: "BioNeMo DiffDock — NVIDIA AI"
- Explanation: "Another AI predicts whether each new molecule will actually stick to the VCP protein. It uses the same kind of AI that powers image generators like DALL-E, but for 3D molecular shapes."
- Key metric card:
  - "Docking Score — measures binding strength"
  - "Below -8.0 kcal/mol = excellent"
  - Big highlight number: "-11.4 kcal/mol" in bold green
  - "Best candidate binds 41% stronger than the original drug"

**Step 3.5 — Would It Make a Good Pill? (Drug-Likeness)**
- Visual: Pill/capsule icon with a checklist overlay
- Label: "RDKit — Drug-Likeness Screening"
- Explanation: "Not every molecule that sticks to a protein would make a good medicine. RDKit checks practical questions: Is it the right size? Can it be absorbed? Can it cross cell membranes? Can it be taken as a pill?"
- Three check items with green checkmarks:
  - "Lipinski's Rule of Five — Predicts if it can be a pill (87 of 98 pass)"
  - "QED Score — Overall drug-likeness: 0.81 out of 1.0 (above 0.67 = drug-like)"
  - "TPSA — Can cross cell membranes: Yes"

**Stage 3 output — Large highlighted result card (the "money shot"):**
- Visual: Trophy or medal icon, warm orange and gold gradient border
- Big number: "100" in bold gradient orange text
- Label: "Ranked Drug Candidates"
- Sub-label: "Novel VCP inhibitor molecules, ranked by a composite score"
- Scoring formula shown visually (three colored bars adding up):
  - Green bar (30%): "Generation Confidence — How confident the AI is in the molecule design"
  - Teal bar (40%): "Binding Strength — How well it sticks to the protein"
  - Orange bar (30%): "Drug-Likeness — How practical it is as a medicine"
- Final output: "PDF Report with Full Provenance — Every tool version, every database, every threshold documented"

**Stage 3 timing badge:** Orange rounded pill — "Total: 8-16 Minutes"

---

### ━━━ SECTION 6: RESULTS COMPARISON — THE "WOW" MOMENT (prominent section) ━━━

**Layout:** This section should be visually striking — it's the climax of the story. Dark navy (#1B2333) background for contrast and drama.

**Header (white text on navy):** "The AI-Designed Molecule vs. the Original Drug"

**Visual:** Side-by-side comparison card — the most impactful visual in the entire infographic.

**Left card — "Original: CB-5083" (gray border, muted styling):**
- Molecular structure icon
- "Existing drug from clinical trials"
- Binding: -8.1 kcal/mol
- Drug-Likeness (QED): 0.62
- Composite Score: 0.64

**Right card — "AI-Designed: Top Candidate" (green/gold border, glowing styling, slightly larger):**
- Molecular structure icon with sparkle/star effect
- "Novel molecule designed by AI in minutes"
- Binding: -11.4 kcal/mol
- Drug-Likeness (QED): 0.81
- Composite Score: 0.89

**Between the two cards — Improvement arrows with large percentages:**
- Arrow pointing up: "+41% Stronger Binding"
- Arrow pointing up: "+31% More Drug-Like"
- Arrow pointing up: "+39% Better Overall Score"

**Important caveat (small text at bottom of this section, honest and transparent):**
- "These are computational predictions — promising starting points for laboratory testing, not finished medicines. Real drug development requires years of laboratory and clinical validation."

---

### ━━━ SECTION 7: TIMELINE COMPARISON — TRADITIONAL vs. HCLS AI FACTORY ━━━

**Layout:** Wide horizontal section with a dramatic before-and-after comparison. Light gray (#F5F5F5) background.

**Header:** "What Used to Take Months — Now Takes Hours"

**Visual design:** Two parallel horizontal timelines, one above the other.

**Top timeline — "Traditional Approach" (gray, muted, long):**
- A long bar stretching most of the width, with labeled segments:
  - "Sequence Alignment: 12-24 hours" (gray segment)
  - "Variant Calling: 8-12 hours" (gray segment)
  - "Annotation: Days (manual)" (gray segment)
  - "Target Identification: Weeks (literature review)" (gray segment)
  - "Drug Candidate Design: Months (medicinal chemistry)" (gray segment)
- Total label: "WEEKS TO MONTHS" in bold gray
- Dollar sign indicator: "Requires expensive compute clusters and large teams"

**Bottom timeline — "HCLS AI Factory" (green, vibrant, short):**
- A dramatically shorter bar, with labeled segments:
  - "GPU Alignment: 20-45 min" (green segment)
  - "AI Variant Calling: 10-35 min" (green segment)
  - "Automated Annotation: Minutes" (teal segment)
  - "AI Target Discovery: Minutes" (teal segment)
  - "AI Drug Design: 8-16 min" (orange segment)
- Total label: "UNDER 5 HOURS" in bold NVIDIA green
- Badge: "On a single $3,999 desktop workstation"

**The visual contrast between the two timeline lengths IS the message.** The traditional bar should be 5-10x longer than the HCLS AI Factory bar.

---

### ━━━ SECTION 8: THE HARDWARE — APPROACHABLE AND IMPRESSIVE ━━━

**Layout:** Clean section with a visual of the DGX Spark. White background.

**Header:** "A Supercomputer on Your Desk" in bold navy

**Visual:** Stylized illustration of the DGX Spark (a compact desktop form factor — looks like a premium PC, not a server rack). Should look approachable and impressive simultaneously.

**Key specs arranged around the illustration (friendly badge style):**

Badge 1 (NVIDIA green): "GB10 GPU — Grace Blackwell Superchip"
- "The chip that makes everything fast"

Badge 2 (teal): "128 GB Unified Memory"
- "CPU and GPU share memory — no data copying bottleneck"

Badge 3 (navy): "144 ARM64 CPU Cores"
- "For tasks that don't need the GPU"

Badge 4 (emerald): "NVMe Storage"
- "Fast enough for 200 GB of genomic data"

**Price callout (PROMINENT — this is a key selling point):**
- Very large, bold: "$3,999"
- Below it: "Everything you just saw runs on this"
- Subtle comparison: "Traditional genomics compute: $50,000-$500,000+"

**Scaling roadmap (three small cards in a row, showing growth path):**

Card 1 (green border, "YOU ARE HERE" badge):
- "DGX Spark — $3,999"
- "Proof of concept, research labs, clinics"
- "What you just saw in this infographic"

Card 2 (teal border):
- "DGX B200 — ~$500K-$1M"
- "Hospital departments"
- "Multiple patients simultaneously"

Card 3 (navy border):
- "DGX SuperPOD — $7M-$60M+"
- "Large health systems"
- "Thousands of patients, federated learning"

---

### ━━━ SECTION 9: THE BIGGER PICTURE (bottom section before footer) ━━━

**Layout:** Medium section showing the broader ecosystem. Light teal tint background.

**Header:** "Part of a Larger Vision" in bold navy

**Three capability cards in a row:**

Card 1 — Genomics (green accent, checkmark badge indicating "available now"):
- Icon: DNA helix
- "Genomic Intelligence"
- "From DNA to drug candidates"
- "What this infographic shows"

Card 2 — Imaging (teal accent, "coming soon" badge):
- Icon: CT scan/brain scan image
- "Imaging Intelligence"
- "CT scans, MRI, X-rays analyzed by AI"
- "MONAI Deploy on DGX Spark"

Card 3 — Cross-Modal (orange accent, "coming soon" badge):
- Icon: Connected nodes/integration
- "Cross-Modal Triggers"
- "A suspicious finding on a CT scan automatically triggers genomic analysis"
- "Imaging + Genomics working together"

**Federated learning callout (small card below):**
- Icon: Network/connected hospitals
- "NVIDIA FLARE — Federated Learning"
- "Multiple hospitals train AI models together without sharing patient data — privacy preserved"

---

### ━━━ SECTION 10: FOOTER (bottom of canvas) ━━━

**Layout:** Full-width dark navy (#1B2333) bar.

**Left:** "HCLS AI Factory" in white, "Open Source — Apache 2.0" in NVIDIA green below

**Center:** Key badges in a row:
- "11.7M Variants" (green pill)
- "201 Genes" (teal pill)
- "100 Drug Candidates" (orange pill)
- "< 5 Hours" (emerald pill)
- "$3,999" (white pill on navy)

**Right:**
- "Author: Adam Jones | February 2026"
- "hcls-ai-factory.org" in NVIDIA green

---

## VCP / FRONTOTEMPORAL DEMENTIA — STORY THREAD

Throughout the infographic, the VCP/FTD case should be woven as a consistent narrative thread. Use a subtle visual indicator (like a small teal dot or "VCP Demo" mini-badge) next to every element that directly relates to the demo case, so viewers can follow the specific example through the entire pipeline:

- Section 2: "In this demo, we analyze a patient's whole genome"
- Section 3: "11.7 million variants found in this patient's DNA"
- Section 4: "AI discovers a disease-causing change in the VCP gene → linked to Frontotemporal Dementia"
- Section 5: "AI designs 100 new molecules targeting VCP"
- Section 6: "The best AI-designed molecule outperforms the existing drug by 39%"

This thread gives the infographic a concrete, real example rather than abstract descriptions.

---

## VISUAL FLOW AND NAVIGATION

The viewer's eye should follow a clear path through the infographic:

1. **Title** (top) — "What is this about?"
2. **Patient/DNA** (left) — "Where does it start?"
3. **Stage 1** (reading DNA) — "How do we read it?"
4. **Stage 2** (finding targets) — "What do we find?"
5. **Stage 3** (designing drugs) — "What do we do about it?"
6. **Results comparison** (climax) — "How good is it?"
7. **Timeline comparison** — "How fast is it?"
8. **Hardware** — "What does it run on?"
9. **Bigger picture** — "What's next?"

Use large, flowing curved arrows between major sections. The arrows should feel organic and inviting, not rigid. Color them with a gradient that transitions between stage colors (green → teal → orange) to reinforce the journey progression.

---

## LANGUAGE GUIDELINES

Every piece of text in this infographic should follow these rules:

1. **No unexplained jargon.** If a technical term is used, immediately explain it in parentheses or in a subtitle. Example: "Variant calling (finding the differences in your DNA)"
2. **Use analogies.** Connect unfamiliar concepts to everyday experiences:
   - DNA sequencing → "Reading a 3.1-billion-letter book by looking at random pages"
   - Alignment → "Assembling a 3.1-billion-piece jigsaw puzzle"
   - Variant calling → "Proofreading for differences"
   - Vector database → "Google Search for your DNA"
   - Molecule generation → "An architect designing 100 variations of a house"
   - Docking → "Testing if a key fits a lock"
   - Drug-likeness → "Would this actually work as a pill?"
3. **Lead with outcomes, not process.** Say "finds disease-causing changes with 99% accuracy" before explaining the CNN architecture.
4. **Use specific numbers.** Audiences remember "11.7 million variants" better than "millions of variants."
5. **Acknowledge limitations honestly.** Include the caveat that these are computational predictions, not finished drugs.

---

## KEY NUMBERS THAT MUST APPEAR (these are the metrics audiences remember)

These numbers should be visually prominent — large, bold, color-coded:

| Number | Context | Where It Appears |
|---|---|---|
| 3.1 billion | Letters in the human genome | Section 2 |
| ~200 GB | Data per patient | Section 2, Section 3 |
| 11.7 million | Variants called | Section 3 output, Section 4 funnel top |
| 3.5 million | High-quality variants | Section 4 funnel |
| >99% | DeepVariant accuracy | Section 3 |
| 35,616 | ClinVar matches | Section 4 funnel |
| 847 | Variants in druggable genes | Section 4 funnel bottom |
| 201 | Genes in knowledge base | Section 4 |
| 171 | Druggable targets (85%) | Section 4 |
| 13 | Therapeutic areas covered | Section 4 |
| 0.87 | AlphaMissense score for VCP variant | Section 4 VCP callout |
| 100 | Drug candidates generated | Section 5 output |
| -11.4 kcal/mol | Best docking score | Section 5, Section 6 |
| 0.89 | Best composite score | Section 5, Section 6 |
| +39% | Improvement over seed compound | Section 6 |
| +41% | Stronger binding | Section 6 |
| +31% | More drug-like | Section 6 |
| < 5 hours | End-to-end runtime | Section 7, Footer |
| $3,999 | DGX Spark price | Section 8, Footer |

---

## WHAT THIS INFOGRAPHIC IS NOT

1. **NOT a technical architecture diagram.** No port numbers, no container images, no API endpoints, no schema details. Those belong in the separate technical infographic.
2. **NOT a slide deck.** Everything is on one canvas, one poster.
3. **NOT a marketing brochure.** It's an honest, accurate, visually rich explanation of real technology with real results. Include the caveat about computational predictions vs. finished drugs.
4. **NOT overwhelming.** Despite the amount of content, the visual hierarchy should make it easy to scan. Someone should be able to get the key message in 10 seconds ("DNA → AI → Drug Candidates in 5 hours on a $3,999 computer") and then dive deeper into any section that interests them.

---

## FINAL NOTES FOR THE AI GENERATOR

1. **This is ONE infographic on ONE canvas.** Landscape 16:9. Not a slide deck. Not separate pages.
2. **The primary audience is non-technical.** Healthcare executives, hospital administrators, conference attendees at HIMSS (25,000+ attendees), GTC, and similar events. They need to understand the VALUE and the JOURNEY, not the technical implementation.
3. **The VCP/FTD demo is the concrete example** woven throughout — it gives the audience something specific to follow.
4. **The three "wow" moments** viewers should walk away with:
   - "11.7 million variants narrowed to 847 druggable targets" (the funnel)
   - "AI-designed molecule is 39% better than the existing drug" (the comparison)
   - "All of this runs on a $3,999 desktop in under 5 hours" (the hardware)
5. **Color consistency matters.** Green = GPU/NVIDIA, Teal = AI/data, Orange = drug discovery, Purple = databases/knowledge, Navy = structure/authority, Emerald = positive outcomes.
6. **The annotation funnel in Section 4 is the visual centerpiece.** It's the single most powerful visual for showing how AI finds signal in noise.
7. **The timeline comparison in Section 7 is the second most powerful visual.** The dramatic length difference between traditional and HCLS AI Factory timelines tells the story instantly.
8. **The $3,999 price point should surprise people.** Make it visually prominent in Section 8 — it's a key differentiator.
9. **White space is your friend.** Unlike the technical architecture infographic, this one should breathe. Let each section have room. The audience is scanning, not studying.
10. **Every section should be understandable in isolation.** Someone's eye might land on any section first — each one should make sense on its own while contributing to the overall flow.
