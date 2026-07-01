# Nano Banana Pro — Neurology Intelligence Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Neurology Intelligence Agent — the most scale-rich agent in the HCLS AI Factory, featuring 10 validated clinical scales (NIHSS, GCS, MoCA, UPDRS-III, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr), 8 clinical workflows spanning stroke to neuromuscular disease, 14 Milvus collections across 10 neurological domains, 42 neuro drugs, 38 neuro genes, and a 5-tab Streamlit UI. It connects to 5 peer agents and plays critical roles in Demo 3 ("Protecting the Survivor" — vincristine neuropathy) and Demo 6 ("When the Standard Drug Can't Be Used" — posterior fossa syndrome, proton therapy). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as all other HCLS AI Factory agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a neurology department reference poster — comprehensive enough for a stroke neurologist, detailed enough for an epileptologist, clean enough for a hospital executive.

**Typography:**
- Title: Large, bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Section headers: Bold, navy, green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body: Dark gray (#333333), 8-10pt, SHORT PHRASES only
- Metric badges: White text on colored pills
- Scale score outputs: Color-coded green/amber/red for severity
- All text MUST be legible at 1920x1080

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Red: #DC2626 — severe/critical scale scores, acute stroke alerts
- Amber: #F5A623 — moderate severity, watch indicators
- Emerald: #059669 — mild/normal, favorable outcomes
- Blue: #2196F3 — cerebrovascular, imaging collections
- Purple: #7B2D8E — degenerative, oncology, genomic collections
- Orange: #FF9800 — epilepsy, movement disorder collections
- Pink: #EC4899 — headache, MS collections
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — brain, neuron, spine, EEG, MRI motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for acute alerts
- Color-coded severity scales (green to amber to red) on each calculator output
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Neurology Intelligence Agent" (green #76B900, white text)
- "10 Validated Clinical Scales" (teal #1AAFCC, white text)

**Center:**
- Line 1 (large, bold, navy): **Neurology Intelligence Agent**
- Line 2 (gray): "Comprehensive Neurological Intelligence on NVIDIA DGX Spark"
- Line 3 (teal): "Part of the Precision Intelligence Engine — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB | 10 Scales | 8 Workflows | 10 Domains | 42 Drugs | 38 Genes | 14 Collections"

**Right side — Clinical capability badge:**
```
NEUROLOGICAL DEPTH
━━━━━━━━━━━━━━━━━━━━━
10 Validated Clinical Scales
8 Specialty Workflows
10 Neurological Domains
Pediatric Neuro-Oncology
Deterministic scale scoring
Demo 3: "Protecting the Survivor"
Demo 6: "When Standard Drug Can't Be Used"
```

---

### ━━━ LEFT COLUMN: 10 NEUROLOGICAL DOMAINS ━━━

**Stacked input cards representing the 10 domains that feed 14 Milvus collections:**

1. **Cerebrovascular** [brain+vessel icon]
   Stroke (ischemic, hemorrhagic), TIA, SAH
   [neuro_cerebrovascular]

2. **Neurodegenerative** [neuron-decay icon]
   Alzheimer's, Parkinson's, ALS, FTD, Huntington's, MSA, DLB
   15 diseases | [neuro_degenerative]

3. **Epilepsy** [lightning-brain icon]
   12 syndromes, EEG patterns, AED selection
   [neuro_epilepsy]

4. **Neuro-Oncology** [brain+tumor icon]
   Glioblastoma, medulloblastoma, DMG (H3K27M)
   [neuro_oncology]

5. **Multiple Sclerosis** [myelin icon]
   RRMS, SPMS, PPMS, 3 DMT tiers, NEDA-3
   [neuro_ms]

6. **Movement Disorders** [tremor icon]
   Parkinson's, essential tremor, dystonia, chorea
   [neuro_movement]

7. **Headache** [head-pain icon]
   8 ICHD-3 classifications, migraine, cluster, TACs
   [neuro_headache]

8. **Neuromuscular** [muscle-nerve icon]
   MG, GBS, CIDP, SMA, motor neuron disease
   [neuro_neuromuscular]

9. **Electrophysiology** [EEG icon]
   35+ EEG patterns, EMG/NCS protocols
   [neuro_electrophysiology]

10. **Neuroimaging** [MRI icon]
    55 imaging protocols, brain/spine MRI
    [neuro_imaging]

Additional: [neuro_literature] [neuro_trials] [neuro_guidelines] [genomic_evidence]

Arrows rightward: "14 Collections | BGE-small-en-v1.5 | 384-dim"

---

### ━━━ CENTER-TOP: 10 CLINICAL SCALES (the headline feature) ━━━

**Section header (navy, green underline):** "10 Validated Clinical Scales — Deterministic, Published Criteria"

**Visual: 10 scale cards in a 5x2 grid. Each shows name, range, severity thresholds (green/amber/red), clinical use.**

**Row 1:**

```
NIHSS                GCS                  MoCA                 UPDRS-III            EDSS
Stroke Severity      Consciousness        Cognition            Parkinson's Motor    MS Disability
Range: 0-42          Range: 3-15          Range: 0-30          Range: 0-132         Range: 0-10
11 items             3 items (E+V+M)      7 domains + edu      13 motor items       8 functional sys

[0-4 Minor]          [13-15 Mild]         [26-30 Normal]       [0-32 Mild]          [0-3.5 Mild]
[5-15 Moderate]      [9-12 Moderate]      [18-25 MCI]          [33-58 Moderate]     [4.0-6.5 Mod]
[16-20 Mod-Severe]   [3-8 Severe]         [10-17 Moderate]     [59+ Severe]         [7.0+ Severe]
[21-42 Severe]                            [<10 Severe]
```

**Row 2:**

```
mRS                  HIT-6                ALSFRS-R             ASPECTS              Hoehn-Yahr
Rankin Outcome       Headache Impact      ALS Function         Stroke CT            PD Staging
Range: 0-6           Range: 36-78         Range: 0-48          Range: 0-10          Range: 0-5
Single score         6 items x 5pt        12 items             10 MCA regions       Stages 0-5

[0-1 Good]           [36-49 Little]       [36-48 Normal]       [8-10 Small core]    [0-2 Early]
[2-3 Moderate]       [50-55 Some]         [24-35 Moderate]     [4-7 Moderate]       [2.5-3 Mid]
[4-5 Severe]         [56-59 Substantial]  [0-23 Severe]        [0-3 Large core]     [4-5 Late]
[6 Dead]             [60-78 Severe]
```

**Callout (emerald border):**
```
ALL 10 SCALES DETERMINISTIC — Published scoring, Python calculators, LLM never computes
```

---

### ━━━ CENTER: 8 WORKFLOWS + KNOWLEDGE ━━━

**Left half: 8 Workflows**

```
1. Acute Stroke Triage      → NIHSS + onset + ASPECTS + LVO → tPA/thrombectomy
2. Dementia Evaluation       → MoCA + age + edu + symptoms → AD/FTD/DLB/vascular
3. Epilepsy Classification   → Seizure + onset + EEG + MRI → Syndrome + AED
4. Brain Tumor Grading       → Histology + molecular + KPS → Grade + treatment
5. MS Assessment             → EDSS + relapses + lesions + JCV → NEDA-3 + DMT tier
6. Parkinson's Assessment    → UPDRS + H-Y + duration + meds → Stage + optimization
7. Headache Classification   → Description + duration + aura → ICHD-3 + treatment
8. Neuromuscular Evaluation  → Weakness + EMG + CK + antibodies → Diagnosis
```

**Right half: Knowledge Stats**
```
[42] Neuro Drugs    [38] Neuro Genes    [15] Degenerative Diseases
[12] Epilepsy Syndromes    [6] Stroke Protocols    [8] Headache Types
[3] MS DMT Tiers    [35+] EEG Patterns    [55] Imaging Protocols
```

---

### ━━━ CENTER-BOTTOM: PEDIATRIC NEURO-ONCOLOGY + CROSS-AGENT ━━━

**Pediatric box (emerald border):**
```
PEDIATRIC NEURO-ONCOLOGY
• MTX Leukoencephalopathy: 3-10% of pediatric ALL, MTHFR risk, periventricular T2 on MRI
• Vincristine Neuropathy: 30-40%, 2mg cap, foot drop/constipation, recovery ~3 months
• Posterior Fossa Syndrome: 25-30% post-resection, cerebellar mutism, 1-6 month recovery
• Cranial Radiation: IQ decline 2-4 pts/year, proton reduces cochlear 95%, hippocampal 80%
• Pediatric Tumors: Medulloblastoma (WNT/SHH/Grp3/Grp4), DMG (H3K27M, <10% survival)
```

**Cross-Agent box (teal border):**
```
Cross-Agent (/v1/neuro/integrated-assessment)
→ Imaging (:8524)      Brain/spine MRI baseline
→ Cardiology (:8126)   Shared toxicity monitoring
→ Biomarker (:8529)    Neurotoxicity genotypes (MTHFR, CYP3A5)
→ Trial (:8538)        Brain tumor trials (SJMB12, Pediatric MATCH)
→ Rare Disease (:8134) Neurogenetic conditions
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

1. **Scale Assessment** — 10-scale color-coded severity report
2. **Workflow Results** — 8 workflow outputs with evidence
3. **Integrated Assessment** — 5-agent coordinated output
4. **Reports** — Markdown | JSON | PDF | FHIR R4

**Streamlit UI (:8529) — 5 Tabs:**
```
1. Dashboard (health, metrics, domains)
2. Evidence Explorer (domain focus selector: 10 options)
3. Clinical Scales (10 dynamic calculators, unique input per scale)
4. Workflow Runner (8 workflows, dynamic inputs)
5. Reports & Export (multi-format)
```

**FastAPI (:8528):**
```
/scale/calculate | /stroke/triage | /dementia/evaluate
/epilepsy/classify | /tumor/grade | /ms/assess
/parkinsons/assess | /headache/classify | /neuromuscular/evaluate
/v1/neuro/integrated-assessment
```

---

### ━━━ BOTTOM: DEMO + SCALE MAPPING ━━━

**Left — Demos:**
```
Demo 3 ✓ "Protecting the Survivor" — Vincristine neuropathy 30-40%, 2mg cap
Demo 6 ✓ "When Standard Drug Can't Be Used" — Posterior fossa, proton therapy, IQ decline
```

**Right — Scale-to-Decision map:**
```
NIHSS → tPA/thrombectomy    GCS → Intubation/ICP    MoCA → Dementia type
UPDRS → Med adjustment      EDSS → DMT escalation   mRS → Rehab level
HIT-6 → Preventive therapy  ALSFRS → Vent/PEG       ASPECTS → Thrombectomy
Hoehn-Yahr → Care planning
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

Green (#76B900), white text:
```
DGX Spark (GB10) $4,699 | Milvus 2.4, 14 collections | BGE-small 384-dim | Claude Sonnet 4 | Ports 8529, 8528 | 208 tests
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 | hcls-ai-factory.org | v2.0 March 2026
Scales: NIHSS (Brott 1989) | GCS (Teasdale 1974) | MoCA (Nasreddine 2005) | UPDRS (Goetz 2008) | EDSS (Kurtzke 1983)
mRS (van Swieten 1988) | HIT-6 (Kosinski 2003) | ALSFRS-R (Cedarbaum 1999) | ASPECTS (Barber 2000) | Hoehn-Yahr (1967)
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **10 clinical scales** — most of any agent. Each with range, severity thresholds (green/amber/red), and clinical decision link
2. **8 specialty workflows** — acute stroke to neuromuscular, each with specific inputs and outputs
3. **10 neurological domains** — broadest domain coverage of any specialty agent
4. **14 Milvus collections** — one per domain plus literature, trials, guidelines, genomic
5. **42 drugs, 38 genes, 15 degenerative diseases, 12 epilepsy syndromes** — verified counts
6. **Pediatric neuro-oncology** — MTX leukoencephalopathy, vincristine neuropathy, posterior fossa syndrome, proton therapy
7. **5 peer agents** — Imaging, Cardiology, Biomarker, Trial, Rare Disease
8. **Severity scoring is visual** — green/amber/red on every scale
9. **Demo 3 + Demo 6** — vincristine neuropathy for Marcus, posterior fossa and proton therapy for Aiden
10. **Dense enough for a neurology ward, protective enough to guard a child's developing brain**
