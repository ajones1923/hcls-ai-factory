# Neurology Intelligence Agent -- White Paper

## Bridging the Neurological Data Fragmentation Gap with RAG-Driven Clinical Intelligence

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones
**Affiliation:** HCLS AI Factory

---

## Abstract

Neurological disorders affect approximately 3 billion people worldwide, making them the leading cause of disability globally (Lancet Neurology, 2024). Despite this burden, neurological clinical decision-making remains fragmented across subspecialty silos -- cerebrovascular, neurodegenerative, epilepsy, movement disorders, multiple sclerosis, headache, neuromuscular, and neuro-oncology -- each with its own terminology, classification systems, and evidence base. This paper presents the Neurology Intelligence Agent, a Retrieval-Augmented Generation (RAG) system that unifies 14 domain-specific vector collections, 10 validated clinical scale calculators, and 8 evidence-based clinical workflows into a single platform that delivers guideline-grounded recommendations in under five seconds. Deployed on NVIDIA DGX Spark as part of the HCLS AI Factory, the agent demonstrates that AI-driven evidence synthesis can address the fragmentation problem at the point of care.

---

## 1. The Problem: Neurological Data Fragmentation

### 1.1 Scale of the Challenge

The World Health Organization's 2024 Global Burden of Disease study estimates that neurological conditions affect 3.4 billion people -- nearly 43% of the global population. Stroke alone kills 6.6 million annually. Alzheimer's disease affects 55 million people worldwide with projections reaching 139 million by 2050. Epilepsy affects 50 million. Parkinson's disease affects 8.5 million. Multiple sclerosis affects 2.8 million. Migraine, the most prevalent neurological disorder, affects over 1 billion.

### 1.2 Fragmented Evidence Landscape

Neurological evidence is distributed across:

- **10+ subspecialty domains**, each with distinct classification systems (ILAE for epilepsy, ICHD-3 for headache, WHO 2021 for CNS tumors, McDonald 2017 for MS)
- **8+ guideline-issuing bodies** (AAN, AHA/ASA, ILAE, IHS, MDS, NCCN, NIA-AA, ACNS), each publishing independently
- **Multiple data modalities**: neuroimaging (MRI, CT, PET), electrophysiology (EEG, EMG, NCS), cerebrospinal fluid biomarkers, genomic variants, and clinical scales
- **Rapidly evolving therapeutics**: anti-amyloid antibodies (lecanemab, donanemab), CGRP-targeted therapies, gene therapies (tofersen for SOD1-ALS, nusinersen for SMA), and novel surgical approaches (LITT, focused ultrasound)

A neurologist evaluating a patient with new-onset seizures must simultaneously consider the ILAE 2017 classification, EEG patterns, MRI findings, genetic testing results, drug interaction profiles, and surgical candidacy criteria. This information exists across separate databases, journals, and guidelines.

### 1.3 Time-Critical Decisions

Acute stroke management epitomizes the urgency: IV tPA must be administered within 4.5 hours (NINDS, ECASS III), mechanical thrombectomy within 6-24 hours (DAWN, DEFUSE-3), and every minute of delay costs 1.9 million neurons. The clinician must compute NIHSS, evaluate ASPECTS, check contraindications, and determine eligibility criteria -- all while coordinating imaging, interventional neurology, and ICU resources.

### 1.4 The Information Overload Problem

Neurologists face approximately 15,000 new publications annually in their field. No individual can maintain current knowledge across all 10 domains, 43+ drug mechanisms, 38+ neurological genes, 12 epilepsy syndromes, 15 neurodegenerative diseases, 6 stroke protocols, and 8 headache classifications while seeing patients.

---

## 2. The Solution: RAG-Driven Clinical Intelligence

### 2.1 Architectural Approach

The Neurology Intelligence Agent addresses fragmentation through a three-layer architecture:

1. **Domain-Specific Vector Collections (14):** Rather than a single monolithic database, neurological knowledge is organized into 14 specialized Milvus collections. Each collection has domain-specific schemas capturing the unique metadata fields relevant to its subspecialty (e.g., `time_window` and `scoring_scales` for cerebrovascular disease; `eeg_pattern` and `first_line_aed` for epilepsy).

2. **Workflow-Specific Weight Boosting (9 profiles):** When a clinical question is classified into a specific workflow (e.g., acute stroke), the search weight for the most relevant collection (neuro_cerebrovascular: 0.25) is boosted relative to less relevant collections (neuro_headache: 0.03). This ensures domain-relevant evidence surfaces preferentially.

3. **Validated Clinical Scale Calculators (10):** Integrated calculators eliminate ad-hoc scoring. The NIHSS calculator accepts 15 item scores and returns not just a total but an interpretation, severity category, clinical thresholds, and prioritized recommendations including tPA dosing guidance.

### 2.2 Query Expansion for Medical Terminology

Neurological terminology is exceptionally complex, with extensive abbreviation systems, brand-to-generic drug name mappings, and syndrome eponyms. The agent's query expansion system addresses this through:

- **251+ entity aliases** mapping abbreviations to canonical forms (e.g., "bvFTD" to "behavioral variant frontotemporal dementia")
- **16 synonym maps** with domain-specific term clusters (e.g., "stroke" maps to "cerebrovascular accident", "CVA", "brain attack", "cerebral infarction")
- **Workflow term injection** adding relevant search terms based on detected or specified clinical context

### 2.3 Multi-Collection Parallel Search

The RAG engine searches all 14 collections simultaneously using Python's ThreadPoolExecutor. Results are scored by weighted cosine similarity, filtered by a configurable threshold (default 0.4), and ranked for presentation. This approach reduces search latency from 14 sequential queries (~700ms) to a single parallel execution (~200-350ms).

### 2.4 Evidence-Grounded LLM Synthesis

Retrieved evidence is synthesized by Claude (Anthropic) using a neurology-specific system prompt that instructs the model to cite guidelines, classify evidence levels, and flag clinical urgency. The LLM operates as a final synthesis layer, never as the primary source of clinical knowledge.

---

## 3. Clinical Capabilities

### 3.1 Acute Stroke Decision Support

The acute stroke workflow computes NIHSS and ASPECTS scores, evaluates tPA eligibility (NINDS/ECASS III criteria), assesses thrombectomy candidacy (DAWN 6-24h, DEFUSE-3 6-16h), classifies stroke etiology (TOAST), and generates time-critical recommendations with door-to-needle and door-to-groin benchmarks.

### 3.2 Memory Clinic / Alzheimer's Assessment

The dementia workflow implements the NIA-AA ATN biomarker framework, scoring Amyloid (A), Tau (T), and Neurodegeneration (N) binary markers to stage patients across 8 possible combinations. It integrates MoCA cognitive screening, identifies anti-amyloid therapy candidates (lecanemab/donanemab), and generates differential diagnosis across 9 dementia subtypes.

### 3.3 Drug-Resistant Epilepsy Evaluation

The epilepsy workflow applies ILAE 2017 seizure classification, identifies epilepsy syndromes from 12 defined patterns, evaluates EEG-MRI concordance, assesses surgical candidacy (anterior temporal lobectomy, laser ablation, VNS, RNS), and flags contraindicated anti-seizure medications for specific syndromes (e.g., sodium channel blockers in Dravet).

### 3.4 Brain Tumor Molecular Grading

The neuro-oncology workflow applies WHO 2021 CNS tumor classification with integrated molecular markers (IDH, MGMT methylation, 1p/19q codeletion, H3K27M, TERT, ATRX, BRAF V600E, EGFR amplification). Treatment protocols are matched to molecular profile and tumor grade.

### 3.5 MS Disease Monitoring

The MS workflow tracks NEDA-3 (No Evidence of Disease Activity) status: no relapses, no new/enlarging T2 lesions, no EDSS progression. It evaluates DMT escalation based on three-tiered efficacy classification (platform, moderate, high), assesses JCV/PML risk for natalizumab, and monitors serum NfL as a biomarker of disease activity.

---

## 4. Technical Implementation

### 4.1 Embedding Strategy

All collections use BAAI/bge-small-en-v1.5 (384 dimensions), selected for its strong performance on medical text benchmarks while fitting within the memory constraints of DGX Spark deployment. IVF_FLAT indexing with nlist=128 provides sub-100ms search latency at the estimated 855K record scale.

### 4.2 Collection Schema Design

Each collection is designed with domain-specific metadata fields. For example, `neuro_cerebrovascular` includes `time_window` (critical for treatment eligibility), `treatment_acute` and `treatment_secondary` (distinct management phases), and `scoring_scales` (linking to relevant calculators). This enables metadata-based filtering in addition to semantic search.

### 4.3 Validated Scale Implementations

All 10 clinical scale calculators are implemented as classmethod-based calculators that accept structured input, validate ranges, compute scores, and return `ScaleResult` objects with interpretation, thresholds, and recommendations. Each calculator is backed by 3-5 unit tests covering boundary conditions and clinical decision points.

### 4.4 Test Coverage

The agent includes 209 automated tests across 12 modules, covering:
- All 18 enums and 12 Pydantic models (55 model tests)
- All 10 scale calculators at clinical thresholds (35 scale tests)
- Knowledge base integrity (30 knowledge tests verifying drug/gene/disease counts)
- Configuration validation (18 settings tests)
- Integration testing (16 end-to-end workflow tests)

---

## 5. Results and Impact

### 5.1 Coverage Metrics

| Metric | Value |
|---|---|
| Disease domains covered | 10 |
| Clinical conditions modeled | 58 |
| Drugs with brand/generic mapping | 43 |
| Genes with disease associations | 38 |
| Clinical scale calculators | 10 |
| Epilepsy syndromes with ASM guidance | 12 |
| Stroke protocols with eligibility criteria | 6 |
| Neurodegenerative diseases with staging | 15 |

### 5.2 Performance

| Metric | Value |
|---|---|
| Query-to-response latency | 2-4 seconds |
| Parallel search latency (14 collections) | 200-350 ms |
| Scale calculation latency | < 5 ms |
| Knowledge base coverage | 855K estimated records |

### 5.3 Clinical Utility

The agent has been validated against published clinical cases demonstrating:
- Correct NIHSS scoring and tPA eligibility determination for acute stroke scenarios
- Accurate ATN staging and anti-amyloid therapy eligibility for memory clinic cases
- Appropriate identification of contraindicated ASMs for specific epilepsy syndromes
- Correct molecular classification and treatment matching for CNS tumors

---

## 6. Three-Engine Architecture and Cross-Agent Coordination

### 6.1 HCLS AI Factory Three-Engine Architecture

The Neurology Intelligence Agent operates within the HCLS AI Factory's three-engine architecture:

1. **Genomic Foundation Engine**: Parabricks/DeepVariant/BWA-MEM2 for FASTQ-to-VCF processing. Provides genomic variant data including neurogenetic variants (SOD1 for ALS, HTT for Huntington's, FMR1 for Fragile X) to the neurology agent via the shared `genomic_evidence` collection.

2. **Precision Intelligence Network**: 11 specialized RAG agents including the Neurology Intelligence Agent. Cross-agent communication enables integrated assessments spanning neuroimaging, cardiac evaluation, biomarker trending, clinical trial matching, and rare disease diagnosis.

3. **Therapeutic Discovery Engine**: BioNeMo/MolMIM/DiffDock/RDKit for drug candidate generation. Neurological therapeutic targets identified by this agent (e.g., novel anti-amyloid epitopes, LRRK2 kinase inhibitor targets) feed into the drug discovery pipeline.

### 6.2 Cross-Agent Coordination

The neurology agent calls 5 sibling agents for integrated patient assessments:

- **Imaging Agent**: MRI lesion characterization, DWI/FLAIR mismatch quantification, tumor segmentation
- **Cardiology Agent**: Cardioembolic stroke evaluation (AF detection, PFO assessment, echocardiography)
- **Biomarker Agent**: CSF biomarker trending (amyloid, tau, NfL trajectories)
- **Clinical Trial Agent**: Trial matching for anti-amyloid therapies, gene therapies, novel DMTs
- **Rare Disease Agent**: Rare neurogenetic disorder evaluation (inherited ataxias, leukodystrophies)

### 6.3 Pediatric Neuro-Oncology Applications

The neurology agent provides specialized support for neurotoxicity in pediatric oncology patients:

- **Methotrexate leukoencephalopathy**: Occurs in 3-10% of pediatric ALL patients receiving intrathecal/high-dose MTX. The agent identifies white matter changes on MRI and correlates with MTX exposure timelines.
- **Vincristine neuropathy**: Affects 30-40% of children on vincristine-containing regimens. The agent tracks peripheral neuropathy severity using age-adapted scales and recommends dose modifications.
- **Posterior reversible encephalopathy syndrome (PRES)**: Associated with chemotherapy, calcineurin inhibitors, and hypertension in pediatric transplant patients. The agent flags characteristic imaging patterns and risk factors.
- **L-asparaginase cerebral thrombosis**: Monitors for coagulopathy-related CNS events in pediatric ALL, including sagittal sinus thrombosis and cortical vein thrombosis.
- **Cranial radiation effects**: Tracks cognitive decline, leukoencephalopathy, and secondary tumor risk in children who received cranial radiation, with long-term neurocognitive monitoring protocols.
- **CNS tumors in children**: Specialized support for medulloblastoma (molecular subgrouping: WNT, SHH, Group 3, Group 4) and diffuse midline glioma (H3K27M-altered), the most common pediatric brain tumors.

---

## 7. Limitations and Future Work

### 7.1 Current Limitations

- Knowledge base relies on seed data from guidelines and landmark trials; institutional clinical data integration is required for production deployment
- No direct EHR integration; clinical data must be manually entered
- Single-node Milvus deployment limits horizontal scalability
- LLM synthesis depends on external API availability (Anthropic Claude)

### 7.2 Future Directions

- **EHR Integration:** FHIR R4 DiagnosticReport output is already supported; bidirectional EHR connectivity is planned for v1.1
- **Real-Time Monitoring:** Continuous EEG pattern classification and stroke alert integration
- **Multi-Institutional Knowledge Federation:** Federated learning across hospital systems while preserving patient privacy
- **Longitudinal Patient Tracking:** Disease progression modeling with automated scale re-assessment
- **Imaging AI Integration:** DWI/FLAIR mismatch quantification, volumetric hippocampal analysis, tumor segmentation

---

## 8. Conclusion

The Neurology Intelligence Agent demonstrates that domain-specific RAG architecture -- with carefully designed collection schemas, workflow-aware search weight boosting, validated clinical scale calculators, and comprehensive query expansion -- can address the fundamental challenge of neurological data fragmentation. By unifying evidence from 10 subspecialty domains, 8 guideline bodies, and multiple data modalities into a single platform that responds in under five seconds, the agent brings the promise of AI-driven clinical decision support to the 3 billion people affected by neurological disease.

---

## References

1. GBD 2021 Nervous System Disorders Collaborators. Global, regional, and national burden of disorders affecting the nervous system, 1990-2021. Lancet Neurology 2024;23:344-381.
2. Powers WJ, et al. Guidelines for the Early Management of Patients With Acute Ischemic Stroke. Stroke 2019;50:e344-e418.
3. Nogueira RG, et al. Thrombectomy 6 to 24 Hours after Stroke with a Mismatch between Deficit and Infarct. NEJM 2018;378:11-21. (DAWN)
4. Albers GW, et al. Thrombectomy for Stroke at 6 to 16 Hours with Selection by Perfusion Imaging. NEJM 2018;378:708-718. (DEFUSE-3)
5. Jack CR Jr, et al. NIA-AA Research Framework: Toward a biological definition of Alzheimer's disease. Alzheimers Dement 2018;14:535-562.
6. van Dyck CH, et al. Lecanemab in Early Alzheimer's Disease. NEJM 2023;388:9-21. (CLARITY AD)
7. Fisher RS, et al. Operational classification of seizure types by the International League Against Epilepsy. Epilepsia 2017;58:522-530.
8. Thompson AJ, et al. Diagnosis of multiple sclerosis: 2017 revisions of the McDonald criteria. Lancet Neurology 2018;17:162-173.
9. Louis DN, et al. The 2021 WHO Classification of Tumors of the Central Nervous System. Neuro Oncol 2021;23:1231-1251.
10. Headache Classification Committee of the IHS. The International Classification of Headache Disorders, 3rd edition. Cephalalgia 2018;38:1-211.

---

*Neurology Intelligence Agent -- White Paper v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
