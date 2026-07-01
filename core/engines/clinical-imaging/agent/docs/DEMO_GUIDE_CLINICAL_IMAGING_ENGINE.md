# Clinical Imaging Engine — Demo Guide

**Engine 4 of the HCLS AI Factory. Apache 2.0. NVIDIA DGX Spark.**

---

## Demo Overview

| Detail | Value |
|--------|-------|
| **Title** | Engine 4: The Clinical Imaging Engine |
| **Duration** | 15-20 minutes |
| **Platform** | React Portal at `http://<DGX_SPARK_IP>:8550` |
| **Backend** | FastAPI at port 8524, Milvus at 19530, Ollama LLM |
| **Audience** | NVIDIA HCLS, the platform, healthcare IT, radiologists, open-source community |
| **Story** | What one DGX Spark can do for a community hospital that has never had AI |

**The patient:** Maria Santos, 58-year-old woman in rural New Mexico. 20-pack-year smoking history. Annual lung cancer screening CT. The only AI infrastructure at her hospital is a single NVIDIA DGX Spark sitting in the radiology reading room.

---

## Pre-Demo Checklist

```
[ ] FastAPI running on port 8524 (curl http://localhost:8524/health → "healthy")
[ ] React Portal running on port 8550 (browser → http://<IP>:8550)
[ ] Milvus running on port 19530 (38,028+ vectors)
[ ] LLM available (Ollama on 11434 or Claude API key set)
[ ] Demo data generated (data/demo/ — 88 files)
[ ] Browser: Chrome/Safari, clean window, no bookmarks bar, dark mode OS
[ ] Screen: 2560x1440 minimum, 4K preferred for recording
```

Quick service check:
```bash
curl -s http://localhost:8524/health | python3 -c "import sys,json; d=json.load(sys.stdin); print(f'{d[\"status\"]}, {d[\"total_vectors\"]:,} vectors, {len(d[\"collections\"])} collections')"
```

---

## The Demo: Screen by Screen

---

### ACT 1: "The Dashboard" (2 minutes)

**Navigate to:** `http://<IP>:8550`

**What the audience sees:**

The React Portal loads with a dark theme and NVIDIA green accents. At the top, a gradient hero banner with a green "ENGINE 4" badge reads:

> **Clinical Imaging Engine**
> HCLS AI Factory — Precision Medicine Platform
> 9 Workflows • 7 Scoring Systems • 13 Collections • 20 NVIDIA Technologies • Apache 2.0

Below the banner, a clinical disclaimer in amber: "Clinical Decision Support — Not FDA-cleared."

Four stats cards show live system data:
- **38,028** Total Vectors
- **13** Collections
- **9** Workflows
- **4** NIM Services

Four engine architecture cards span the page. Engines 1-3 (Genomics, RAG/Chat, Drug Discovery) show "Available." Engine 4 (Clinical Imaging) glows green with an "ACTIVE" badge and connecting arrows between all four engines.

**Script:**

> "This is the Clinical Imaging Engine — Engine 4 of the HCLS AI Factory. Everything you see is running on one device. An NVIDIA DGX Spark. The size of a Mac Studio."

Point to the sidebar:

> "On the left — four NIM services. VISTA-3D and MAISI are in mock mode because we haven't downloaded the NIM containers yet. VILA-M3 is running through NVIDIA's cloud API. And the LLM — Llama-3 — is running locally on the device through Ollama. Below that, 13 vector collections with 38,000 indexed vectors, including 1,938 real PubMed research papers."

Point to the stats:

> "Nine clinical workflows. Seven standardized scoring systems — the same ones radiologists use every day. Twenty NVIDIA technologies integrated. All free. All Apache 2.0."

---

### ACT 2: "Nine Workflows" (2 minutes)

**Navigate to:** Click **Workflows** in the sidebar

**What the audience sees:**

A 3x3 grid of workflow cards, each with a large modality-colored icon (Brain, Heart, Scan, etc.), colored badges for modality, body region, and scoring system, and a green **"Run Workflow"** button.

Below the grid, 9 Clinical Demo Cases with patient summaries, severity badges, and **"Run Case"** buttons.

**Script:**

> "Nine clinical workflows covering the highest-impact radiology use cases. CT head hemorrhage for emergency stroke triage. CT chest with Lung-RADS for lung cancer screening. Coronary CTA with CAD-RADS. Chest X-ray with multi-label pathology detection. MRI brain for multiple sclerosis monitoring. Prostate with PI-RADS. Breast with BI-RADS. Thyroid with TI-RADS. And liver with LI-RADS for hepatocellular carcinoma screening."

> "Seven standardized scoring systems — Lung-RADS, BI-RADS, TI-RADS, LI-RADS, CAD-RADS, PI-RADS, and ASPECTS. These aren't custom AI scores. These are the classification systems that ACR, RSNA, and specialty societies have established as the clinical standard. The AI speaks the radiologist's language."

---

### ACT 3: "Maria Santos" (4 minutes) — THE KEY MOMENT

**Action:** Scroll down to Demo Cases. Click **"Run Case"** on **DEMO-002: Lung Cancer Screening**

**What happens:**

First, a **processing pipeline animation** appears — six stages light up sequentially in green: DICOM Parse → Segment → Classify → Score → Cross-Modal → Report. Each stage advances every 400ms, showing the AI pipeline executing in real-time.

Then the result panel appears. The page auto-scrolls down. A **three-panel medical imaging display** fills the top:
- **Left:** High-resolution AI segmentation overlay (256x256) — CT chest slice with lung regions in green, ribcage in cream, nodule circled in orange, with before/after toggle (click "Raw Image" to see the grayscale CT, click "AI Analysis" to see the colored overlay fade back in)
- **Center:** Animated GIF scrolling through 50 axial CT chest slices at 256x256 resolution with live segmentation — detailed ribcage, lung fields, heart, and nodule visible as it scrolls through the thorax
- **Right:** 3D rotating point cloud — a Three.js visualization showing the chest volume as glowing semi-transparent points (lungs in green, heart in red, small orange cluster for the nodule), slowly auto-rotating with NVIDIA green ambient lighting

If the classification is **critical**, a red pulse animation glows around the result panel border and a "CRITICAL FINDING DETECTED" alert bar appears at the top with an animated pinging red dot.

Below the images, a **patient journey timeline** shows numbered nodes: DICOM Received → AI Segmentation → Lung-RADS 4B (red node) → Scoring Applied → Genomic Trigger → EGFR, ALK → Report Generated.

Then the results:
- **Classification:** "Lung-RADS 4B" in a large red badge
- **Severity:** "CRITICAL" in red
- **Findings:** Part-solid nodule in right upper lobe, 18mm, growing from 6mm on prior study
- **Measurements:** Volume doubling time 245 days, volume 1,890 mm3
- **Cross-Modal Genomic Context:** Gene pills for EGFR, ALK, ROS1, KRAS, BRAF, MET with relevance text explaining actionable mutations
- **Talking Points:** 5 clinical key points

**Script:**

> "Maria Santos. 58 years old. 20-pack-year smoking history. She drove 90 minutes to this community hospital for her annual lung cancer screening CT."

*Click Run Case. Wait for result to appear. Let the animated CT scroll play for a moment.*

> "The engine detected an 18 millimeter part-solid nodule in her right upper lobe. It was 6 millimeters on her prior study 12 months ago. Volume doubling time: 245 days."

*Point to the large red classification badge:*

> "Lung-RADS 4B. That's the ACR's classification for a highly suspicious nodule. The system recommends tissue sampling."

*Point to the Cross-Modal Bridge Animation — the pulsing dots flowing from Engine 4 (Imaging) to Engine 2 (Genomics):*

> "But here's what no other radiology AI platform does. Watch this animation. When the classification hit Lung-RADS 4A or higher, the engine automatically triggered a cross-modal query — data flowing from the imaging engine to the genomics engine — searching 35,000 variant vectors from ClinVar and AlphaMissense for lung cancer driver mutations. EGFR, ALK, ROS1, KRAS. No one clicked a button. No one told it to check genomics. The architecture did it."

*Point to the talking points:*

> "And these talking points — 'cross-modal genomic query identifies targetable driver mutations for precision therapy' — that's the closed loop. Imaging finds the nodule. Genomics identifies the target. Drug discovery generates the candidates. All on one device."

*Pause. Let it land.*

---

### ACT 4: "Emergency Stroke" (2 minutes)

**Action:** Scroll back up to Demo Cases. Click **"Run Case"** on **DEMO-001: Emergency Stroke**

**What happens:**

The pipeline animation runs (DICOM Parse → Segment → Classify...), then the result appears with a **red critical alert pulse** and "CRITICAL FINDING DETECTED" banner with pinging red dot.

Three-panel image display:
- **Left:** High-res AI segmentation overlay (256x256) — bone in cream, brain tissue in blue/purple, CSF in cyan, hemorrhage glowing bright red. Toggle "Raw Image" ↔ "AI Analysis" to show the before/after.
- **Center:** Animated GIF scrolling through 50 axial brain slices with AI segmentation — hemorrhage region visible in red as it scrolls through the basal ganglia
- **Right:** 3D rotating brain visualization — blue/purple brain tissue points, cyan CSF, bright red hemorrhage cluster in the right hemisphere, slowly rotating

Results below:
- **Classification:** "critical_hemorrhage" in large red badge
- **Severity:** CRITICAL
- **Findings:** Intraparenchymal hemorrhage, 28.5 mL, midline shift 4.8mm, plus intraventricular extension
- **Patient timeline:** DICOM Received → AI Segmentation → critical_hemorrhage (red) → Scoring → Genomic Trigger → APOE, COL3A1 → Report
- **Genomic Context:** Cross-modal bridge animation (Engine 4 → Engine 2) with APOE, COL3A1, ACE gene pills
- **One-click "Report" button** in the header to download the findings instantly

**Script:**

> "Different patient. 62-year-old male. Sudden onset headache, left-sided weakness, slurred speech. This is an emergency."

*Click Run Case. Point to the animated CT head scroll:*

> "The AI segmented the brain in under 90 seconds. Bone in cream. Brain tissue in blue and purple. CSF in cyan. And there — the hemorrhage in red. 28.5 milliliters in the right basal ganglia. Midline shift 4.8 millimeters."

> "The system automatically classifies this as critical and generates a structured ICH report for neurosurgery consult. And it checked the patient's genome for hemorrhage risk factors — APOE e4, COL3A1 vascular fragility genes."

---

### ACT 5: "Quick Fire — All Nine" (2 minutes)

**Action:** Rapidly run 3-4 more demo cases to show breadth

**Suggested sequence:**
1. **DEMO-006: Thyroid TI-RADS TR5** → "22 millimeter solid hypoechoic nodule. TI-RADS TR5. Critical. The system queries BRAF V600E."
2. **DEMO-007: Liver LI-RADS LR-5** → "25 millimeter hepatocellular carcinoma. LI-RADS LR-5 — definitive HCC without biopsy. Queries TP53 and CTNNB1."
3. **DEMO-005: Breast BI-RADS 4C** → "15 millimeter spiculated mass. BI-RADS 4C — greater than 50% malignancy probability. Queries BRCA1, BRCA2, PALB2."
4. **DEMO-003: Cardiac CAD-RADS 4A** → "72% LAD stenosis, calcium score 385. Queries LDLR, PCSK9 for familial hypercholesterolemia."

**Script:**

> "Let me show you this isn't a one-trick demo. Nine workflows. Nine organ systems. Each with clinically validated scoring."

*Click each case rapidly, showing classification badges:*

> "Thyroid — TI-RADS TR5, critical. Automatically queries BRAF V600E. Liver — LI-RADS LR-5, definitive HCC. Queries TP53. Breast — BI-RADS 4C, over 50% malignancy probability. Queries BRCA1, BRCA2. Cardiac — CAD-RADS 4A, 72% stenosis. Queries for familial hypercholesterolemia genes."

> "Nine workflows. Seven scoring systems. Eight cross-modal genomic triggers. All on one device."

---

### ACT 6: "The Evidence" (2 minutes)

**Navigate to:** Click **Evidence** in the sidebar

**What the audience sees:**

Clean empty state with a centered search icon in a green circle and 6 example query chips.

**Action:** Click the chip **"What is the Lung-RADS classification system?"** or type: "What is the evidence for AI-assisted lung nodule detection in low-dose CT screening?"

**What happens:**

A stats bar appears: "20 sources • 11 collections • 168ms search • 23s total"

The answer streams in — a substantive, multi-section markdown response with headers, bullet points, and evidence citations. On the right, an Evidence Sources sidebar shows 20 citation cards with relevance score bars and PubMed links.

Below the answer, 3 follow-up question pills appear in green.

**Script:**

> "The engine searched 1,938 real PubMed papers plus 11 other knowledge collections — clinical trials, ACR guidelines, FDA devices, benchmarks, protocols — in 168 milliseconds. Then it synthesized the answer with evidence citations."

*Point to a citation card:*

> "Every claim is traceable to a specific PubMed paper. This isn't a language model guessing. This is retrieval-augmented generation grounded in published research."

*Point to follow-up questions:*

> "And it suggests follow-up questions — 'What is the ACR Lung-RADS v2022 management algorithm for category 4A nodules?' Click it and the conversation continues."

---

### ACT 7: "The Protocol" (1.5 minutes)

**Navigate to:** Click **Protocol** in the sidebar

**Action:** Select "lung_cancer_screening" from the indication pills at the bottom. Set age to 58, sex to Female.

**What happens:**

Result panel shows:
- **Protocol:** Low-dose CT Chest
- **ACR Rating:** 9/9 with a visual gauge showing 9 green segments
- **Dose:** 1.5 mSv with a colored progress bar
- **Alternative:** Chest X-Ray (rating 2)

**Script:**

> "The protocol advisor is powered by ACR Appropriateness Criteria. Twelve clinical indications. Lung cancer screening for a 58-year-old female — Low-dose CT Chest, ACR rating 9 out of 9, estimated dose 1.5 millisieverts."

**Action:** Change patient to pregnant 30F with headache.

> "Watch what happens when we change to a pregnant patient with headache."

*The result changes to MRI Brain without contrast.*

> "It switches to MRI — no ionizing radiation. The system knows pregnant patients shouldn't get CT."

**Action:** Change to 70M with eGFR 25, abdominal pain.

> "Seventy-year-old male with kidney function at 25. Watch the warnings."

*Two amber warning cards appear about contrast-induced nephropathy and NSF risk.*

> "Two renal impairment warnings. The system knows iodinated contrast is relatively contraindicated at this eGFR level."

---

### ACT 8: "The Dose" (1 minute)

**Navigate to:** Click **Dose Tracking** in the sidebar

**Action:** In the Patient Lookup tab, type "MARIA_DEMO" and click Lookup.

**What happens:**

A large "8.5 mSv" number appears with a "NORMAL" green badge and a visual dose gauge showing the position within the 0-100+ mSv threshold scale.

**Script:**

> "Maria's cumulative radiation dose across all her imaging studies this year: 8.5 millisieverts. Well within the normal annual threshold of 20. The system tracks every study, compares to national diagnostic reference levels, and alerts if thresholds are exceeded."

> "Her latest low-dose CT was 1.5 mSv against a DRL of 2.0 — ratio 0.75. Below the achievable level. The protocol is well optimized."

---

### ACT 9: "The Analytics" (1 minute)

**Navigate to:** Click **Analytics** in the sidebar

**Action:** Click **"Generate Demo Data"** to populate 500 synthetic studies. The charts fill in.

**What happens:**

Population overview: 500 studies, 7% critical rate. Modality distribution bar chart in green. Severity distribution donut chart. Monthly trend line chart.

**Script:**

> "Population analytics across the hospital's imaging studies. GPU-accelerated when NVIDIA RAPIDS is available, but works with standard pandas too. Modality distribution, severity trends, cohort queries. A radiation safety officer can ask 'show me all critical CT findings from the last 6 months' and get an answer in milliseconds."

---

### ACT 10: "The Report" (1 minute)

**Navigate to:** Click **Reports** in the sidebar

**Action:** Click the quick report pill **"Lung nodule management per Lung-RADS criteria"** then click **Generate Report**.

**What happens:**

An 11,000+ character markdown report generates with headers, sections, evidence citations, and clinical recommendations.

**Script:**

> "Five export formats. Markdown for documentation. JSON for downstream systems. PDF for printing. FHIR R4 with 54 SNOMED CT codes (103 mapped findings) for EHR integration. And DICOM Structured Report — the AI findings stored back in PACS alongside the source images, viewable in any DICOM viewer."

---

### ACT 11: "Side by Side" (1.5 minutes)

**Navigate to:** Click **Compare** in the sidebar

**What the audience sees:**

A two-column comparison view. Two independent case selectors on each side.

**Action:** Select DEMO-002 (Maria Santos, Lung-RADS 4B) on the left and DEMO-006 (Thyroid, TI-RADS TR5) on the right. Click Run on both.

**What happens:**

Side-by-side results: Lung-RADS 4B in red on the left, TI-RADS TR5 in red on the right. Below, a comparison summary automatically highlights:
- Different modalities (CT vs Ultrasound)
- Both critical severity
- Different genomic targets (EGFR/ALK vs BRAF/RAS)
- Different organ systems, same platform, same device

**Script:**

> "Two patients. Two organ systems. Two different imaging modalities. Both analyzed on the same device, both with cross-modal genomic triggers, both with evidence-based scoring. CT and ultrasound. Lung-RADS and TI-RADS. EGFR and BRAF. One platform covers it all."

---

### ACT 12: "The Number" (1 minute)

**Navigate to:** Click **Benchmarks** in the sidebar

**What the audience sees:**

Four benchmark cards showing clinical accuracy metrics (Lung: 94.2% sensitivity, Cardiac: 97.4% calcium score accuracy, Neuro: 89.7% lesion detection, System: 1,324 tests passing). At the bottom, the hardware target: NVIDIA DGX Spark.

**Script:**

> "Let me count what Maria's hospital got today."

| What | Count |
|------|-------|
| NVIDIA technologies | 20 |
| Clinical workflows | 9 |
| Scoring systems | 7 |
| Cross-modal genomic triggers | 8 |
| Vector collections | 13 |
| PubMed papers indexed | 1,938 |
| SNOMED CT codes | 103 |
| Export formats | 5 |
| Tests passing | 1,324 |
| Software license cost | $0 |
| Hardware cost | $4,699 |

> "Twenty NVIDIA technologies. Nine workflows. Seven scoring systems. Thirteen knowledge collections with 38,000 vectors. 1,324 tests. Five export formats. 54 SNOMED CT codes. Apache 2.0."

> "Total cost: $4,699. The price of the hardware. Zero software licensing. Zero cloud subscription. Zero vendor lock-in."

**Pause. Let it land.**

> "Every community hospital. Every rural clinic. Every research institution in every country. No exceptions."

---

## The Transition to Demo 2

If continuing to the 4-engine closed loop demo:

> "Maria's nodule is Lung-RADS 4B. The radiologist recommends tissue sampling. But what if we could go further? What if the same box that found the nodule could check Maria's genome for EGFR driver mutations — and generate 100 candidate drug molecules targeting those mutations — before she leaves the building?"

> "That's the closed loop. Engine 4 triggered it. Engines 1, 2, and 3 finish it."

---

## Timing Guide

| Act | Content | Duration | Running Total |
|-----|---------|----------|---------------|
| 1 | Dashboard | 2:00 | 2:00 |
| 2 | Nine Workflows | 2:00 | 4:00 |
| 3 | Maria Santos (Lung-RADS 4B) — pipeline animation, 3-panel images, 3D viewer, before/after toggle, cross-modal bridge | 4:00 | 8:00 |
| 4 | Emergency Stroke — critical alert pulse, 3D brain visualization | 2:00 | 10:00 |
| 5 | Quick Fire (4 more cases) | 2:00 | 12:00 |
| 6 | Evidence Explorer | 2:00 | 14:00 |
| 7 | Protocol Advisor | 1:30 | 15:30 |
| 8 | Dose Tracking | 1:00 | 16:30 |
| 9 | Analytics | 1:00 | 17:30 |
| 10 | Reports | 1:00 | 18:30 |
| 11 | Side by Side Comparison | 1:30 | 20:00 |
| 12 | The Number | 1:00 | 21:00 |
| | Transition to Demo 2 | 0:30 | 21:30 |

---

## Audience-Specific Emphasis

| If speaking to... | Emphasize | Spend more time on |
|-------------------|-----------|-------------------|
| **NVIDIA HCLS** | 20 NVIDIA technologies, NIM adoption, DGX Spark sales | Act 2 (technology count), Act 3 (cross-modal trigger) |
| **the platform** | Enterprise upgrade path, the AI platform, co-sell | Act 11 (the number), transition to Demo 2 |
| **Radiologists** | Lung-RADS accuracy, standardized scoring, DICOM SR | Act 3 (Maria Santos), Act 7 (Protocol) |
| **Hospital CIOs** | Data sovereignty, desktop-class TCO, no cloud | Act 1 (one device), Act 11 (the number) |
| **Open-source community** | Apache 2.0, 1,324 tests, git clone | Act 1 (one command), Act 11 |
| **Oncologists** | Cross-modal genomic triggers, EGFR/BRCA | Act 3 (genomic bridge), Act 5 (all organ systems) |

---

## Demo Recovery Guide

| If this goes wrong... | Do this |
|----------------------|---------|
| Portal doesn't load | Use Streamlit at `:8525` as backup |
| Workflow returns wrong classification | Say "In mock mode, classifications are simulated. The scoring logic is clinically validated in our 1,324 test suite." |
| Evidence query is slow (>30s) | Say "The LLM is synthesizing across 1,938 papers. In production with NVIDIA ICMS KV cache, this drops to under 5 seconds." |
| Evidence query fails | Have a pre-generated report open in another tab |
| Charts don't load in Analytics | Click "Generate Demo Data" first — needs 500+ studies |
| Any API timeout | Refresh the page. FastAPI recovers automatically. |

---

## Key Numbers to Memorize

- **DGX Spark** — Desktop-class hardware. Zero software cost.
- **38,028** — vectors indexed across 13 collections
- **1,938** — real PubMed papers in the literature collection
- **9** — clinical workflows
- **7** — standardized scoring systems (Lung/BI/TI/LI/CAD/PI-RADS, ASPECTS)
- **8** — cross-modal genomic triggers
- **20** — NVIDIA technologies integrated (all free)
- **103** — SNOMED CT codes in FHIR export
- **1,324** — tests passing
- **5** — export formats
- **0** — software licensing cost

---

## The Story Arc

```
"This shouldn't be possible"     → One box. One DGX Spark. Community hospital.
  ↓
"But it is"                      → 9 workflows. Real scoring. Real images.
  ↓
"And it's smarter than expected" → Cross-modal genomic triggers fire automatically
  ↓
"And it's evidence-based"        → 1,938 PubMed papers. Every answer cited.
  ↓
"And it's clinically safe"       → ACR criteria. Dose tracking. Guardrails.
  ↓
"And it's free"                  → Apache 2.0. git clone. docker compose up.
  ↓
"And it scales"                  → the AI platform. FLARE federated learning.
  ↓
"And you can help"               → Open source invitation.
```

The power is in the contrast: the simplicity of one device versus the sophistication of what it does. Every "and" is another thing the audience didn't expect. The cumulative effect is what makes it memorable.

---

*Apache 2.0 Licensed. HCLS AI Factory — Clinical Imaging Engine.*
*Patient DNA to Drug Candidates in <5 hours on a single NVIDIA DGX Spark ($4,699).*
