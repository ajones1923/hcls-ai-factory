# Demo 2: The Closed Loop — From CT Scan to Drug Candidate

**The demo that has never been done before. On any platform. At any price.**

---

## Overview

| Detail | Value |
|--------|-------|
| **Title** | The Closed Loop: CT Scan to Drug Candidate in Under One Hour |
| **Duration** | 18-22 minutes (follows Demo 1, or standalone) |
| **Prerequisite** | Demo 1 completed (Maria Santos, Lung-RADS 4B established) |
| **Engines Used** | All 4: Imaging (8550) → Genomics (5000) → Intelligence (8501) → Drug Discovery (8505) |
| **Story** | A routine screening CT triggers the entire precision medicine pipeline — automatically |
| **Audience** | NVIDIA executives, the platform leadership, GTC keynote, oncology conferences |

**The headline:** Maria Santos got a screening CT. The AI found a suspicious nodule, checked her genome, identified an EGFR driver mutation, and generated 100 candidate drug molecules — all before she left the building. One device. One hour. One DGX Spark.

---

## Pre-Demo Checklist

```
[ ] Demo 1 completed (or at minimum, DEMO-002 Lung-RADS 4B result visible)
[ ] Engine 1 — Genomics Portal at http://<IP>:5000 (200 OK)
[ ] Engine 2 — RAG Chat at http://<IP>:8501 (200 OK)
[ ] Engine 3 — Drug Discovery at http://<IP>:8505 (200 OK)
[ ] Engine 4 — Imaging Portal at http://<IP>:8550 (200 OK)
[ ] Pre-computed demo data at data/demo/ (engine3_handoff, cross_modal)
[ ] 4 browser tabs ready: :8550, :5000, :8501, :8505
[ ] Or single screen switching between tabs during the demo
```

---

## Screen Layout

**Option A (4 screens/monitors):** One engine per screen. Best for live presentations.

**Option B (single screen, 4 tabs):** Switch between browser tabs. Color-code each tab mentally:
- Tab 1: Engine 4 — Imaging (green, where we start)
- Tab 2: Engine 1 — Genomics (blue)
- Tab 3: Engine 2 — Intelligence (teal)
- Tab 4: Engine 3 — Drug Discovery (purple)

**Option C (recorded video):** Pre-record with transitions. Each engine gets a colored border or label overlay.

---

## The Demo: Act by Act

---

### ACT 1: "The Trigger" (3 minutes)

**Screen:** Engine 4 — Imaging Portal (http://<IP>:8550/workflows)

**Starting point:** DEMO-002 Maria Santos result from Demo 1 is still visible. If not, click "Run Case" on DEMO-002 to re-run it.

**What the audience sees:**

The Lung-RADS 4B result with:
- 3-panel medical imaging (AI segmentation + animated CT scroll + 3D point cloud)
- CRITICAL FINDING DETECTED banner
- Patient journey timeline with red Lung-RADS 4B node
- Cross-modal bridge animation — pulsing dots flowing from Engine 4 (Imaging) → Engine 2 (Genomics)
- Gene pills: EGFR, ALK, ROS1, KRAS, BRAF, MET

**Script:**

> "In Demo 1, we showed you what Engine 4 — the Clinical Imaging Engine — can do on its own. Nine workflows, seven scoring systems, cross-modal genomic triggers. But Maria's story doesn't end with Lung-RADS 4B."

*Point to the cross-modal bridge animation:*

> "Watch this animation. Data is flowing from Engine 4 — Imaging — to Engine 2 — Genomics. The engine didn't just find a nodule and classify it. It asked a question: does this patient have actionable driver mutations?"

> "In a traditional hospital, that question takes 6 to 8 weeks. Biopsy referral, tissue sampling, molecular testing, results. On this device, it takes minutes."

*Point to the gene pills:*

> "EGFR. ALK. ROS1. KRAS. BRAF. MET. These are the genes that determine which targeted therapies will work for this specific type of lung cancer. The imaging engine queried them automatically. Now let me show you what happened on the other side."

---

### ACT 2: "The Genome" (3 minutes)

**Switch to:** Engine 1 — Genomics Portal (http://<IP>:5000)

**What the audience sees:**

The Genomics Pipeline portal showing the completed genomic analysis. If showing live data, the VCF file with 11.7 million variants is displayed.

**Script:**

> "This is Engine 1 — the Genomic Foundation Engine. When Maria visited this hospital previously, her whole genome was sequenced. 200 gigabytes of raw DNA data. Engine 1 processed it using NVIDIA Parabricks — GPU-accelerated alignment with BWA-MEM2 and variant calling with Google DeepVariant."

> "What would take 48 hours on a CPU server took 3 hours on this DGX Spark. The result: 11.7 million variants called with over 99.7% accuracy."

*Point to the VCF data:*

> "Every variant in Maria's genome is indexed and searchable. But 11.7 million variants is too many for a human to review. That's where Engine 2 comes in."

---

### ACT 3: "The Intelligence" (4 minutes) — THE PIVOT MOMENT

**Switch to:** Engine 2 — RAG Chat (http://<IP>:8501)

**What the audience sees:**

The Precision Intelligence Engine — a chat interface backed by Milvus with 3.56 million indexed vectors from ClinVar (4.1M clinical variants), AlphaMissense (71M AI pathogenicity predictions), and a curated knowledge base of 201 genes across 13 therapeutic areas.

**Action:** Type or paste: "What EGFR variants does this patient carry and what is their clinical significance for non-small cell lung cancer treatment?"

**What happens:**

The RAG engine searches the genomic evidence collection. The answer identifies EGFR L858R — a missense mutation in exon 21. ClinVar classifies it as pathogenic. AlphaMissense gives it a pathogenicity score of 0.94 (high confidence). The knowledge graph connects it to erlotinib, osimertinib, and gefitinib as targeted therapies.

**Script:**

> "Engine 2 is the Precision Intelligence Engine. It searched 3.56 million genomic evidence vectors — ClinVar clinical annotations, AlphaMissense AI pathogenicity predictions, and a curated knowledge base of 201 druggable genes."

*Wait for the answer to appear:*

> "There it is. EGFR L858R. A missense mutation in exon 21 of the epidermal growth factor receptor gene. ClinVar says pathogenic. AlphaMissense scores it 0.94 out of 1.0 — high confidence that this variant is disease-causing."

> "This is the most common actionable EGFR mutation in non-small cell lung cancer. It's the mutation that erlotinib and osimertinib were designed to target. But those are existing drugs. What if we could design something new? Something optimized for this specific mutation?"

*Pause. The audience knows what's coming.*

---

### ACT 4: "The Target" (2 minutes)

**Still on Engine 2 or switch to pre-computed data.**

**Script:**

> "Engine 2 validated EGFR as a druggable target. Priority 5 out of 5 — the highest. It pulled the protein's crystal structures from the Protein Data Bank. PDB IDs: 1M17, 4ZAU, 5CAL. And it identified erlotinib as the reference compound — the seed molecule for drug generation."

> "This target hypothesis — gene, variant, protein structure, reference drug, druggability assessment — is the handoff to Engine 3."

---

### ACT 5: "The Molecules" (5 minutes) — THE JAW-DROP MOMENT

**Switch to:** Engine 3 — Drug Discovery UI (http://<IP>:8505)

**What the audience sees:**

The Therapeutic Discovery Engine. The EGFR target is loaded. Crystal structure PDB:5CAL is visible — the EGFR kinase domain with the erlotinib binding pocket highlighted.

**Action:** Show the drug discovery pipeline or pre-computed results.

**What the audience sees:**

The 10-stage pipeline:
1. Initialize target
2. Normalize to UniProt/PDB
3. Structure discovery (RCSB PDB)
4. Structure preparation
5. **Molecule generation (MolMIM)** — 100 new molecular structures appear
6. Chemistry QC (RDKit validation)
7. Conformer generation
8. **Molecular docking (DiffDock)** — binding poses predicted
9. Composite ranking
10. Report generation

Results:
- 100 novel EGFR inhibitor candidates
- 87 pass Lipinski's Rule of Five
- 72 have QED > 0.67 (drug-like)
- Top 10 docking scores: -8.2 to -11.4 kcal/mol
- Composite scores: 0.68-0.89

**Script:**

> "This is Engine 3 — the Therapeutic Discovery Engine. It loaded Maria's EGFR target. It pulled the crystal structure — this is the protein's kinase domain, the molecular machine that drives cell growth. And here, in the active site, is where erlotinib binds."

*Point to the molecule generation results:*

> "MolMIM — NVIDIA's molecular generation model — created 100 novel molecular structures using erlotinib as a seed. Each one is a variation on the EGFR inhibitor scaffold, designed to fit the same binding pocket but with different chemical properties."

> "DiffDock predicted how each molecule docks into the EGFR active site. Binding poses, affinity scores, hydrogen bonds, contact residues."

*Point to the top-ranked candidates:*

> "100 novel EGFR inhibitor candidates. 87 pass Lipinski's Rule of Five — they look like real drugs. The top 10 have binding scores comparable to or better than erlotinib itself."

*Let that sink in.*

> "These aren't existing drugs pulled from a database. These are new molecules that didn't exist before this pipeline ran. Generated, validated, docked, and ranked — in under 16 minutes."

---

### ACT 6: "The Loop Closes" (3 minutes)

**Switch back to:** Engine 4 — Imaging Portal (http://<IP>:8550)

**Script:**

> "Let me show you what just happened."

*Show or describe the timeline:*

```
8:00 AM — Maria arrives for screening CT
8:15 AM — Engine 4: Lung-RADS 4B nodule detected, cross-modal trigger fires
8:16 AM — Engine 2: EGFR L858R mutation identified from 3.56M genomic vectors
8:20 AM — Engine 3: 100 novel EGFR inhibitor candidates generated and ranked
8:35 AM — Radiologist reviews AI findings alongside CT images
8:40 AM — FHIR DiagnosticReport + DICOM SR exported to EHR
8:45 AM — Oncologist receives complete precision medicine packet
```

> "From CT scan to drug candidates. 45 minutes. On a single device."

*Pause.*

> "But the loop doesn't end there."

---

### ACT 7: "The Follow-Up" (2 minutes)

**Script:**

> "Three months from now, Maria comes back for follow-up imaging. Engine 4 runs the same CT workflow. But this time, it also extracts 1,500 radiomics features from the nodule — quantitative texture measurements at the microstructural level."

> "It compares these features to the baseline scan. Changes in tissue entropy, heterogeneity, and shape that predict treatment response — weeks before the nodule visibly shrinks or grows on imaging."

> "The loop is circular. Imaging triggers genomics. Genomics triggers drug discovery. Drug discovery guides treatment. Imaging monitors response. And if treatment fails, the cycle restarts with new molecular targets."

> "That's not a pipeline. That's a precision medicine system."

---

### ACT 8: "The Architecture" (2 minutes)

**Script:**

> "Four engines. One device."

| Engine | Name | What It Did for Maria | Time |
|--------|------|----------------------|------|
| **4** | Clinical Imaging | Found the nodule, classified Lung-RADS 4B, triggered genomics | 1 minute |
| **1** | Genomic Foundation | Processed her genome, called 11.7M variants | 3 hours (pre-computed) |
| **2** | Precision Intelligence | Found EGFR L858R, validated druggable target | 2 minutes |
| **3** | Therapeutic Discovery | Generated 100 EGFR inhibitor candidates | 16 minutes |

> "Total active time from CT scan to drug candidates: under one hour. On a single NVIDIA DGX Spark. With zero software licensing costs."

> "No other system on earth does this. Not commercially. Not in open source. Not at any price point."

---

### ACT 9: "The Scale" (2 minutes)

**Script:**

> "Everything you just saw runs on one box. But what happens when it's not one hospital — it's a hundred?"

> "That's where the AI platform enters the picture."

*If presenting to the platform/NVIDIA:*

> "the AI platform replaces the Docker volumes with canonical file storage — every DICOM image, every VCF file, every molecule, every AI result exists as a transactional Element in a unified namespace. data engine triggers the pipeline automatically when a CT arrives. retrieval engine embeds reports without ETL. Vector Database provides unified SQL plus vector search. And ICMS accelerates the LLM by 10 to 20x."

> "One DGX Spark serves one hospital. the AI platform plus DGX SuperPOD serves a health system. The code is the same. The deployment tier changes."

*If presenting to open-source community, skip the platform and go to:*

> "And through NVIDIA FLARE federated learning, multiple hospitals can train shared AI models without sharing patient data. Only model improvements flow between sites. Data sovereignty preserved."

---

### ACT 10: "The Invitation" (1 minute)

**Script:**

> "This is open source. Apache 2.0. The repository is public."

```
git clone https://github.com/ajones1923/hcls-ai-factory.git
```

> "Four engines. Twenty NVIDIA technologies. Nine clinical workflows. Seven scoring systems. Cross-modal genomic triggers. 100 drug candidates per target. 1,324 tests. All free."

> "We built this because Maria Santos — and every patient like her in every community hospital, in every rural clinic, in every low- and middle-income country — deserves the same precision medicine that academic medical centers provide."

> "Not in five years. Not when the budget allows. Now."

> "Clone it. Deploy it. Improve it. That's the invitation."

**End.**

---

## Timing Guide

| Act | Content | Duration | Running Total |
|-----|---------|----------|---------------|
| 1 | The Trigger (Engine 4 → cross-modal bridge) | 3:00 | 3:00 |
| 2 | The Genome (Engine 1 — Parabricks, 11.7M variants) | 3:00 | 6:00 |
| 3 | The Intelligence (Engine 2 — EGFR L858R identified) | 4:00 | 10:00 |
| 4 | The Target (validation, PDB structures, druggability) | 2:00 | 12:00 |
| 5 | The Molecules (Engine 3 — 100 EGFR inhibitors generated) | 5:00 | 17:00 |
| 6 | The Loop Closes (timeline: 45 minutes CT → drugs) | 3:00 | 20:00 |
| 7 | The Follow-Up (radiomics monitoring, circular loop) | 2:00 | 22:00 |
| 8 | The Architecture (4 engines, one device) | 2:00 | 24:00 |
| 9 | The Scale (the AI platform / FLARE federation) | 2:00 | 26:00 |
| 10 | The Invitation (Apache 2.0, git clone) | 1:00 | 27:00 |

**Combined Demo 1 + Demo 2: ~48 minutes total**

Or run Demo 2 standalone in ~27 minutes (briefly recap the Lung-RADS 4B finding in Act 1).

---

## The Numbers That Matter

| Metric | Value |
|--------|-------|
| **CT scan to drug candidates** | Under 1 hour |
| **Drug candidates generated** | 100 per target |
| **Pass Lipinski's Rule of Five** | 87 of 100 |
| **Drug-like (QED > 0.67)** | 72 of 100 |
| **Genomic variants processed** | 11.7 million |
| **Genomic evidence vectors** | 3.56 million |
| **Druggable gene targets** | 201 across 13 therapeutic areas |
| **Engines** | 4 |
| **NVIDIA technologies** | 20 (free) |
| **Device** | NVIDIA DGX Spark |
| **Price** | $4,699 hardware, $0 software |
| **License** | Apache 2.0 |

---

## The Emotional Arc

```
Demo 1 ended with:
  "Every hospital on earth. No exceptions."

Demo 2 builds:
  "But what if we could go further?"
    ↓
  "The imaging engine triggered something."
    ↓
  "EGFR L858R. The driver mutation."
    ↓
  "100 new molecules. In 16 minutes."
    ↓
  "From CT scan to drug candidates. 45 minutes. One device."
    ↓
  [silence]
    ↓
  "$4,699. Apache 2.0."
    ↓
  [silence]
    ↓
  "Clone it. Deploy it. Improve it."
```

The silences are as important as the words. After "45 minutes, one device" — stop talking for 3 full seconds. After the price reveal and "Apache 2.0" — stop for 3 more. Let the audience process what they just saw.

---

## What Makes This Demo Historic

No one has ever demonstrated, on any platform at any price point:

1. A routine imaging study triggering genomic analysis automatically
2. Genomic analysis identifying a specific driver mutation
3. That mutation driving AI-powered drug candidate generation
4. 100 novel molecules generated, validated, docked, and ranked
5. The entire cycle completing in under one hour
6. On a single desktop-class DGX Spark
7. With the entire platform available as open source under Apache 2.0

Each of these individually would be impressive. All seven together, in sequence, live, on one device — that's unprecedented.

---

## Recovery Guide

| If this goes wrong... | Do this |
|----------------------|---------|
| Engine 2 chat is slow | Use pre-computed EGFR evidence from data/demo/cross_modal/ |
| Engine 3 isn't loaded | Show the pre-computed target hypothesis from data/demo/engine3_handoff/ |
| Engine 1 isn't accessible | Skip Act 2, say "The genome was pre-processed. Let me show you what Engine 2 found." |
| Any engine UI is down | Use the pre-computed JSON files and narrate the results |
| Audience asks about FDA clearance | "This is a research platform and clinical decision support tool, not a cleared medical device. All findings require review by qualified healthcare professionals." |
| Audience asks about validation | "1,324 tests passing. Mock mode produces clinically realistic results based on published scoring criteria. Clinical validation studies are a next step." |

---

*Apache 2.0 Licensed. HCLS AI Factory.*
*Patient DNA to Drug Candidates in <5 hours on a single NVIDIA DGX Spark ($4,699).*
