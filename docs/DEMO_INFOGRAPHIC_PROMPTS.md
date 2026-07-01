# HCLS AI Factory — Demo Infographic Prompts for Napkin Pro

## Design Guidelines (Apply to All 6)

**Brand Colors:**
- Navy: #1B2333 (backgrounds, headers)
- Teal: #1AAFCC (accents, connections, highlights)
- NVIDIA Green: #76B900 (success indicators, key metrics)
- White: #FFFFFF (text on dark backgrounds)
- Light Gray: #F5F5F5 (card backgrounds)
- Red: #E74C3C (alerts, critical flags only)

**Style:** Clean, modern, medical-professional. Dark navy backgrounds with teal/green accent elements. Rounded cards for each agent/engine. Flow arrows showing data movement. Icons for each component. Patient avatar at the center or start of each flow.

---

## Infographic 1: "DNA to Drug" — End-to-End Precision Medicine Pipeline

**Prompt:**

Create a detailed horizontal pipeline infographic titled "DNA to Drug — End-to-End Precision Medicine Pipeline" for the HCLS AI Factory platform. Dark navy (#1B2333) background throughout.

At the top left, show a patient card: "Evelyn, Age 8, B-Cell ALL" with a small child silhouette icon in white.

The pipeline flows left to right through 3 major engine blocks, each as a large rounded rectangle with a subtle glow border:

**Block 1 — Genomic Foundation Engine (teal #1AAFCC border)**
Show 4 sequential steps inside with small icons:
- FASTQ icon (DNA helix) → "Raw Sequencing Data"
- GPU chip icon → "BWA-MEM2 Alignment (45 min GPU vs 12 hrs CPU)"
- Neural network icon → "DeepVariant Variant Calling"
- Document icon → "Annotated VCF: 11.7M variants"
Label at bottom: "Powered by NVIDIA Parabricks on DGX Spark"

**Connecting arrow (green #76B900, animated-style dotted line)**

**Block 2 — Precision Intelligence Network (teal border)**
Show in center a Milvus database cylinder icon labeled "3.56M Vectors" with the Precision Oncology Agent card connecting to it:
- Agent card: "Precision Oncology Agent" with a tumor icon
- Inside the card show 3 findings as small pills/badges:
  - "ETV6-RUNX1 Fusion Detected (Favorable)"
  - "IKZF1 Deletion (Unfavorable Modifier)"
  - "Risk: Standard → High (IKZF1plus)"
- Below: "Claude LLM Evidence Synthesis" with RAG retrieval icon
Label at bottom: "11 AI Agents | BGE-small-en-v1.5 Embeddings | 384 Dimensions"

**Connecting arrow (green)**

**Block 3 — Therapeutic Discovery Engine (teal border)**
Show 3 sequential steps:
- Molecule generation icon → "MolMIM: 500 candidates generated"
- Protein-ligand docking icon → "DiffDock: Binding affinity scored"
- Shield/safety icon (green) → "Pediatric Safety Filter: BBB, hERG, hepatic"
- Final output: 3 ranked drug candidate cards showing "Candidate 1: -9.1 kcal/mol ✓ Pediatric Safe"
Label at bottom: "NVIDIA BioNeMo + RDKit"

**Bottom banner (green #76B900 background):**
"Patient DNA → Drug Candidates in <5 Hours | All Three Engines | Open Source (Apache 2.0)"

**Top right corner:** HCLS AI Factory logo, "Built for NVIDIA DGX Spark"

---

## Infographic 2: "The 30-Second Tumor Board" — Multi-Agent Coordination

**Prompt:**

Create a radial/hub-and-spoke infographic titled "The 30-Second Tumor Board — Multi-Agent Coordination" for the HCLS AI Factory platform. Dark navy (#1B2333) background.

**Center hub:** Large circle with patient card: "Evelyn, Age 8, Day 29 Induction, MRD Positive (0.1%)" — use a heartbeat/medical icon. The hub is labeled "Precision Oncology Agent — Orchestrator" with a teal (#1AAFCC) glowing border.

**5 spoke agents radiating outward**, each as a rounded card connected by teal lines with directional arrows flowing both ways:

**Spoke 1 (top):** "Precision Biomarker Agent" — Biomarker icon
- Key finding badge: "CD19+ (MFI 45,200) | CD22+ | MRD 0.1%"

**Spoke 2 (top right):** "CAR-T Intelligence Agent" — T-cell icon
- Key finding badge: "Tisagenlecleucel Eligible | CRS Risk: Moderate (Grade 1-2)"

**Spoke 3 (right):** "Single-Cell Intelligence Agent" — Cell cluster icon
- Key finding badge: "97% CD19+ Blasts | 21% CD19-dim Escape Risk | Immune Desert TME"

**Spoke 4 (bottom right):** "Clinical Trial Intelligence Agent" — Clipboard/trial icon
- Key finding badge: "3 Matched Trials: AALL1731, Pediatric MATCH, AMELIA"

**Spoke 5 (bottom):** "Therapeutic Discovery Engine" — Molecule icon
- Key finding badge: "3 Novel CREBBP Inhibitors | PDB: 5BN4 | Pediatric Safe"

**Bottom section:** A timeline bar showing:
- "Traditional Tumor Board: 1-2 weeks to schedule, hours to discuss"
- vs (with dramatic visual contrast)
- "HCLS AI Factory: <30 seconds, 5 agents, all evidence cited"

**Bottom banner (green):** "5 Specialists. 30 Seconds. Every Recommendation Evidence-Based."

---

## Infographic 3: "Protecting the Survivor" — Cardiotoxicity Prevention

**Prompt:**

Create a vertical shield/protection-themed infographic titled "Protecting the Survivor — Cardiotoxicity Prevention in Pediatric Chemotherapy" for the HCLS AI Factory. Dark navy (#1B2333) background.

**Top:** Patient card: "Marcus, Age 6, High-Risk Neuroblastoma (MYCN Amplified, Stage 4)"

**Visual concept:** A large shield shape in the center, divided into 7 horizontal sections, each representing an agent that protects Marcus. The shield has a green (#76B900) border symbolizing protection.

**Shield Section 1 — Pharmacogenomics Agent:**
- Icon: DNA/pill
- "TPMT *1/*3A → Intermediate Metabolizer → 50% 6-MP Dose Reduction"
- "CYP3A5: Vincristine metabolism screening"

**Shield Section 2 — Precision Biomarker Agent:**
- Icon: Test tube
- "MYCN Amplified | LDH 1,450 U/L | Ferritin 280 ng/mL"
- "COG/INRG Risk: HIGH"

**Shield Section 3 — Cardiology Intelligence Agent:**
- Icon: Heart with shield
- "Cumulative Daunorubicin 300 mg/m² → Dexrazoxane Mandatory"
- "Baseline Echo: LVEF 62%, FS 34% — Normal"
- "Monitoring: Echo every 3 months during therapy"

**Shield Section 4 — Neurology Intelligence Agent:**
- Icon: Brain
- "Vincristine Neuropathy Risk: 30-40%"
- "Absolute Dose Cap: 2 mg"
- "Monitor: Foot drop, constipation, jaw pain"

**Shield Section 5 — Precision Autoimmune Agent:**
- Icon: Immune cell
- "Dinutuximab (Anti-GD2) irAEs:"
- "Neuropathic Pain 85% | Capillary Leak 25% | Hypersensitivity 15%"

**Shield Section 6 — Precision Oncology Agent:**
- Icon: Target
- "Dose-Adjusted COG ANBL1232 Protocol"
- "Integrates all agent findings into modified treatment plan"

**Shield Section 7 — Therapeutic Discovery Engine:**
- Icon: Molecule
- "Novel ALK Inhibitors with Cardiac Safety Profile"
- "PDB: 2XP2 | hERG IC50 >30 μM"

**Bottom section — Statistics bar:**
- "85% of childhood cancer survivors have at least one chronic health condition by age 40"
- "This platform prevents toxicity BEFORE treatment starts"

**Bottom banner (green):** "7 Agents. 1 Child. Every Toxicity Anticipated Before the First Dose."

---

## Infographic 4: "One Gene, One Family" — Rare Disease with Cancer Predisposition

**Prompt:**

Create a family tree / cascade infographic titled "One Gene, One Family — Rare Disease with Cancer Predisposition" for the HCLS AI Factory. Dark navy (#1B2333) background.

**Top center:** Patient card: "Aurora, Age 4, Bilateral Retinoblastoma"

**Visual concept:** A simplified family tree structure showing genetic cascade, with the platform's analysis flowing downward through 5 agent blocks.

**Family tree section (top third):**
- Aurora (center, highlighted in teal): "RB1 c.958C>T p.Arg320Ter — PATHOGENIC"
- Mother (left): "Testing recommended" with question mark icon
- Father (right): "Testing recommended" with question mark icon
- Sibling — 2-year-old brother (below parents): "50% chance of carrying RB1 mutation — CASCADE TESTING CRITICAL" highlighted in red (#E74C3C)
- Future pregnancies note: "PGT-M (preimplantation genetic testing) available"

**Agent flow (bottom two-thirds) — 5 connected cards flowing downward:**

**Card 1 — Rare Disease Diagnostic Agent:**
- "ACMG/AMP Classification: PATHOGENIC"
- "Evidence: PVS1 (null variant) + PM2 (absent from population) + PP4 (phenotype match)"
- "Score: 12/28 criteria met"

**Card 2 — Precision Oncology Agent:**
- "Retinoblastoma Management"
- "Globe-sparing therapy preferred (left eye)"
- "AVOID external beam radiation — increases second cancer risk 3-4x in RB1 carriers"

**Card 3 — Imaging Intelligence Agent:**
- "Trilateral Retinoblastoma Screening"
- "Brain MRI every 6 months until age 5"
- "Pinealoblastoma risk: 5-15% in bilateral RB"

**Card 4 — Clinical Trial Intelligence Agent:**
- "3 Matched Trials"
- "COG ARET0321 | NCT05765045 | NCT04587544"

**Card 5 — Therapeutic Discovery Engine:**
- "Novel CDK4/6 Inhibitors for Ocular Delivery"
- "Optimized for intravitreal injection | MW <400"

**Bottom section — Lifetime surveillance timeline:**
Horizontal timeline from age 0 to 50 showing:
- Age 0-5: Brain MRI q6 months (pinealoblastoma)
- Age 5-20: Annual bone scan (osteosarcoma risk 500x)
- Age 10-30: Soft tissue monitoring (sarcoma risk 100x)
- Age 20-50: Melanoma screening (20x risk)
- Label: "One gene mutation → lifetime of coordinated surveillance"

**Bottom banner (green):** "1 Diagnosis. 1 Family Changed. Surveillance That Lasts a Lifetime."

---

## Infographic 5: "Last Line of Defense" — CAR-T Therapy Decision

**Prompt:**

Create a decision-pathway infographic titled "Last Line of Defense — CAR-T Therapy Decision and Monitoring" for the HCLS AI Factory. Dark navy (#1B2333) background.

**Top:** Patient card with urgency indicator (red border): "Ethan, Age 12, Relapsed/Refractory B-ALL — Failed 2 Prior Lines"

**Visual concept:** A vertical decision tree/pathway flowing downward through 7 decision gates, each represented by an agent. Each gate shows a "GO/NO-GO" assessment. Green checkmarks for cleared gates.

**Gate 1 — Single-Cell Intelligence Agent:**
- Icon: Microscope
- Question: "Is the CAR-T target expressed?"
- Answer: "CD19+ on 97.2% of blasts | MFI 8,500 | ✓ TARGET VALIDATED"
- Secondary: "TME: Immune desert (2% CD8+ T cells)"
- Status: GREEN CHECKMARK

**Gate 2 — CAR-T Intelligence Agent:**
- Icon: Engineered T-cell
- Question: "Is the patient eligible?"
- Answer: "Tisagenlecleucel: ELIGIBLE"
- Details: "ELIANA criteria met | CRS risk: Moderate (Grade 1-2 predicted) | Manufacturing: 22 days"
- "Backup: CD22 CAR-T if CD19 escape"
- Status: GREEN CHECKMARK

**Gate 3 — Cardiology Intelligence Agent:**
- Icon: Heart + ECG
- Question: "Is the heart ready?"
- Answer: "LVEF 58% (normal) | Prior anthracycline 250 mg/m² | Troponin normal"
- "Flu/Cy conditioning: Cardiac risk LOW"
- Status: GREEN CHECKMARK

**Gate 4 — Precision Autoimmune Agent:**
- Icon: Immune shield
- Question: "What immune risks exist?"
- Answer: "CRS monitoring plan: IL-6, IFN-gamma, TNF-alpha"
- "Autoimmune cytopenias: 30-40% expected"
- "B-cell aplasia: 100% (lifelong IVIG replacement)"
- Status: YELLOW WARNING (manageable)

**Gate 5 — Pharmacogenomics Agent:**
- Icon: Pill + DNA
- Question: "Are supportive meds safe?"
- Answer: "Tocilizumab: No CYP interactions ✓"
- "Corticosteroids: CYP3A4 normal ✓"
- "All supportive medications cleared"
- Status: GREEN CHECKMARK

**Gate 6 — Clinical Trial Intelligence Agent:**
- Icon: Clipboard
- Question: "Are there better options?"
- Answer: "3 trials matched: ELIANA, ZUMA-4, Dual CD19/CD22 (NCT03448393)"
- Status: GREEN CHECKMARK

**Gate 7 — Therapeutic Discovery Engine:**
- Icon: Molecule
- Question: "Bridging therapy during 22-day manufacturing?"
- Answer: "3 CD22 modulator candidates | T-cell sparing | Short half-life"
- Status: GREEN CHECKMARK

**Decision output (bottom, large green box):**
"PROCEED WITH TISAGENLECLEUCEL"
"All 7 gates cleared | Manufacturing initiated | Bridging therapy active"

**Bottom banner (green):** "7 Checks. 7 Agents. One Decision That Could Save Ethan's Life."

---

## Infographic 6: "When the Standard Drug Can't Be Used" — Novel Drug Discovery

**Prompt:**

Create a problem-solution split infographic titled "When the Standard Drug Can't Be Used — Novel Drug Discovery for a Growing Child" for the HCLS AI Factory. Dark navy (#1B2333) background.

**Top:** Patient card: "Aiden, Age 10, SHH-Subtype Medulloblastoma (PTCH1 Mutation)"

**Visual concept:** Split into two halves — LEFT side shows the PROBLEM (red-tinted), RIGHT side shows the SOLUTION (green-tinted). A large "X" or barrier symbol divides them.

**LEFT SIDE — The Problem (subtle red #E74C3C tint):**

**Oncology Agent finding:**
- "SHH Pathway: PTCH1 → SMO → GLI → Tumor Growth"
- "Standard targeted therapy: Vismodegib (SMO inhibitor)"

**CRITICAL WARNING (large, impossible to miss):**
- Red border box with warning icon
- "VISMODEGIB CONTRAINDICATED IN GROWING CHILDREN"
- "Blocks Hedgehog signaling in growth plates"
- "Causes PERMANENT growth plate fusion"
- "Aiden's bones would stop growing"

**Additional context cards:**
- Neurology Agent: "Posterior fossa syndrome risk: 25-30% post-surgery"
- Neurology Agent: "Craniospinal radiation: IQ decline 2-4 points/year"
- Imaging Agent: "MRI staging: T3aM0, no spinal drop metastases"

**CENTER DIVIDER:** Large arrow labeled "The Platform Finds Another Way"

**RIGHT SIDE — The Solution (subtle green #76B900 tint):**

**Therapeutic Discovery Engine (large, prominent):**

Step 1 — "MolMIM generates 500 novel SMO antagonists"
- Icon: Molecular generation
- "Conditioned on SMO crystal structure (PDB: 5L7D)"

Step 2 — "Pediatric Safety Filter eliminates 411 candidates"
- Icon: Shield/filter
- Checkboxes: "✓ MW <500 (BBB penetration) | ✓ hERG IC50 >10 μM | ✓ No growth plate toxicity | ✓ Hepatic safety for immature liver"

Step 3 — "DiffDock scores remaining 89 candidates"
- Icon: Protein-ligand docking
- "Binding affinity range: -7.2 to -9.8 kcal/mol"

Step 4 — "Top 3 candidates ranked"
- 3 molecule cards showing:
  - "SMO-PED-001: IC50 8 nM | -9.4 kcal/mol | ✓ Growth plate safe"
  - "SMO-PED-002: IC50 12 nM | -8.9 kcal/mol | ✓ Growth plate safe"  
  - "SMO-PED-003: IC50 18 nM | -8.7 kcal/mol | ✓ Growth plate safe"

**Clinical Trial Agent card (bottom right):**
- "3 Matched Trials: SJMB12 (St. Jude), Pediatric MATCH SHH Arm, Arsenic Trioxide GLI1"

**Bottom banner (green, spanning full width):**
"The standard drug would have harmed a growing child. The platform generated safe alternatives in minutes."

**Final line:** "3 Engines. 5 Agents. 500 Candidates. 3 Safe for a Child. Open Source."

---

## General Notes for All Infographics

- Use clean, minimal icons (line-art style, white or teal on navy)
- Patient names and ages should be prominent — these are stories about children
- Agent names should be clearly labeled with their function
- All clinical data shown is factually verified
- Include "HCLS AI Factory | Apache 2.0 | Built for NVIDIA DGX Spark" on every infographic
- Resolution: Print-quality (300 DPI minimum)
- Aspect ratio: 16:9 for presentation slides, or tall portrait for social media posts
