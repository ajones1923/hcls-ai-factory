# The TSC Intelligence Engine

### From First Principles to a Deployable Disease Engine within the HCLS AI Factory — A Build and Institutional Blueprint for Cincinnati Children's

**Author:** Adam M. Jones — Architect, HCLS AI Factory; Product Marketing Manager, HPC, the platform
**Platform:** HCLS AI Factory — **Engine 7: The TSC Intelligence Engine** · *Seven Engines. Eight Intelligence Agents. One Platform.*
**License:** Apache 2.0 · **Status:** Master research volume (supersedes the prior TSC Explainer, Winslow Initiative, Builder's Guide, and FAQ) · 2026

> **Scope honesty, up front.** Everything described as *running today* runs on a single NVIDIA DGX Spark (with RunPod for extra GPU as needed) against a **synthetic** TSC cohort, watermarked throughout. Integration with Epic, the biobank LIMS, and imaging-AI is described architecturally but is **institutional Phase-1 work, not built**. Every agent is decision support with a human review gate; nothing here is FDA-cleared or validated on real patients. This posture is deliberate — it is the credibility the work depends on.

---

## Implementation Status — v0.1 Built & Verified (June 2026)

Since this document was written, the engine it specifies has been **built and verified** on the NVIDIA DGX Spark. It is no longer only a design.

- **Implemented** — the full repository (`ai_agent_adds/tuberous_clerosis_engine/`, ~63 Python files, ~2,870 LOC): the event-sourced core (13 event types, append-only log + materialized projections), the deterministic orchestrator (dependency-ordered dispatch, conservative failure handling, demand-driven surface assembly), all five agents, the synthetic-cohort pipeline, the RAG retriever, three Streamlit clinician surfaces, and an evaluation harness.
- **Tested** — 41 automated tests pass in under one second; the three-act demonstration runs end-to-end against a deterministic, version-controlled 50-patient synthetic cohort generated on disk.
- **Deterministic and classical components are fully real** — the ACMG-AMP combinatorial classifier (Richards 2015), the Gaussian-process trajectory forecasts, the Marshall-Hagedorn discourse-marker detection, the ITSC surveillance-gap analyzer, and the ClinicalTrials.gov eligibility matcher.
- **Frontier-model reasoning verified live** — with an Anthropic API key configured, the TSC-Variant Curator produced a genuine `claude-opus-4-8` ACMG per-criterion narrative for the featured mosaic patient (correctly explaining the PVS1 strength downgrade for the low-VAF mosaic context), while the deterministic validator remained authoritative for the classification — exactly as this document specifies.
- **Verified demonstration outcomes** — Patient A's 8.3%-VAF TSC2 mosaic variant is recovered from the VCF and classified Likely Pathogenic; Patient B's SEGA forecast crosses the discussion threshold in the 12–18-month window and the academic TAND cluster surfaces as briefing material (never an alert); Patient C's six-section, source-attributed therapeutics brief assembles.

**What remains (the honest edges are unchanged).** All data is synthetic and watermarked. Epic/Clarity, biobank LIMS, and imaging-AI integration remain institutional Phase-1 work (not built). The faithful BAMSurgeon→Parabricks blinded-calling substrate is a documented GPU upgrade. Every output is decision support behind a clinician review gate; nothing is FDA-cleared. Real-data validation, a clinician's review of the synthetic notes, and delivery at Cincinnati Children's are the steps that remain — and they are the right ones to do next.

### Update — expanded build (June 2026)

Work since this status carried the engine well beyond the v0.1 walking skeleton; the file/test counts above are superseded. The repository is now **~86 application Python files (~5,700 LOC) with a 92-test suite (offline, deterministic)**, and the following are implemented and verified on the Spark:

- **Ontology-grounded phenome** — every HPO term validated against the official Human Phenotype Ontology release (19,389 terms), labels corrected, obsolete/alt IDs remapped, non-ontology codes dropped (this caught two miscoded terms in the cohort itself).
- **A measured evaluation harness** (`/eval`) reporting, against the cohort's known ground truth, variant-classification accuracy 1.00 (construct validity), a **+12-point diagnostic detection uplift** from 100% recovery of the six sub-threshold mosaics, zero truncating-variants-called-benign, and 100% provenance completeness — each carrying a "construct validity, not prospective accuracy" caveat.
- **Read-level genomics** — strand-resolved VCFs plus a strand-biased artifact the curator must reject, so "zero false-positive pathogenic" is a real discrimination.
- **Multi-quantity trajectories** — SEGA/AML plus renal function (eGFR) and seizure frequency, with a Bayesian growth slope, a survival-style crossing probability, a graded crossing, and a surveillance-cadence recommendation.
- **Span-grounded clinical NLP** — character-offset HPO/TAND spans with polarity/temporality, a present/absent discordance log, and weighted, recency-aware, negation-filtered TAND scoring.
- **Governance & operations** — a per-tier token/cost ledger (`/system/usage`), an AI-labeled DRAFT report with a molecular-geneticist sign-off gate, incremental note dispatch, and a per-clinician sliding-window alert counter.
- **Production substrate** — statsmodels mixed-effects, the LangGraph runtime (byte-equivalent), and PostgreSQL/Redis/MinIO/Milvus behind env flags, with runs-here fallbacks so the engine runs on a bare Spark.
- **An Omniverse / OpenUSD digital-twin layer (Surface d)** — the Spark authors four scene kinds from the real projections (lesion-trajectory twin with the uncertainty envelope = the prediction intervals; mosaic "powers-of-ten" cell field; whole-child organ atlas; 50-patient population array with the seven recoveries ringed in gold), with MDL materials for film-quality RTX. Authoring is CPU-side and dependency-free; rendering is a RunPod RTX step. The companion *TSC Anatomical Digital Twin* research paper and PRD specify this layer.

The honest edges are unchanged: synthetic data throughout; the faithful BAMSurgeon→Parabricks substrate, Epic/Clarity/LIMS/PACS/IAM/SIEM, frontier-model note/imaging *generation*, and patient-specific (segmented) anatomy remain GPU- or institution-bound; the RTX render runs off-box; nothing is FDA-cleared. Prospective validation on real Cincinnati specimens is the next step.

---

## Table of Contents

**Front Matter** — Abstract, Executive Summary & How to Read This Paper

**Part I — Foundations: TSC and the Opportunity**
1. Introduction & Thesis
2. The Biology, From First Principles
3. The Disease: TSC Across Organ Systems
4. The Clinical & Scientific Opportunity

**Part II — The Platform and Engine 7**
5. The HCLS AI Factory Platform
6. The TSC Intelligence Engine as Engine 7

**Part III — Engine Architecture**
7. Engine Architecture Overview
8. Agent 1 — TSC-Variant Curator
9. Agent 2 — TSC-Phenome Mapper
10. Agent 3 — TSC-Trajectory Modeler
11. Agent 4 — TAND Surveillance Agent (and the Marshall-Hagedorn Methodology)
12. Agent 5 — TSC-Therapeutics Strategist
13. The Orchestrator and Clinician-Facing Surfaces

**Part IV — Building on the DGX Spark + RunPod**
14. Building It on the DGX Spark (with RunPod for Extra GPU)
15. The Synthetic Cohort
16. Model Tiering, the TSC RAG Corpus, Audit & Evaluation
17. Build Sequencing & Engineering Discipline

**Part V — The Cincinnati Children's Case**
18. Cincinnati Children's: The Beachhead, the Winslow Pavilion, and the Touchpoints
19. The Demonstration and the Engagement Path
20. Governance, Ethics, Regulatory & Privacy
21. The Winslow Intelligence Initiative: A 24-Month Institutional Path

**Part VI — Beyond TSC**
22. Beyond TSC: The Replication Thesis
23. Risks, Limitations & Open Questions
24. Conclusion & Call to Action

**Appendices** — A. Glossary · B. Agent I/O Schemas · C. Synthetic Cohort Specification · D. References · E. Mapping to the Prior TSC Papers

---




## Abstract, Executive Summary & How to Read This Paper

### Abstract

Tuberous Sclerosis Complex (TSC) is a multisystem genetic disorder (~1 in 6,000 births) in which dysregulated mTOR signaling drives hamartomas across brain, kidney, lung, heart, and skin. Its care is hard not because the underlying biology is mysterious but because the relevant evidence is fragmented across modalities and time: molecular genetics misses the ~10-15% mosaic "no-mutation-identified" (NMI) cohort that standard blood testing cannot see; structured phenotype data is scattered through years of encounters; longitudinal trajectories (SEGA growth, AML size, seizure burden) are reconstructed by hand at each visit; and TSC-Associated Neuropsychiatric Disorders (TAND), present in ~90% of patients, go 30-50% unaddressed (TOSCA). We describe the TSC Intelligence Engine, Engine 7 of the open-source HCLS AI Factory, a coordinated system of five disease-specific agents and a deterministic orchestrator that turn this dispersed evidence into reviewable clinician-facing surfaces. The agents recover low-VAF somatic mosaic variants and synthesize ACMG-AMP classifications (Variant Curator); extract time-anchored HPO phenotypes (Phenome Mapper); forecast organ trajectories with classical statistics (Trajectory Modeler); surface under-recognized TAND signals via a published diagnostic-uncertainty discourse taxonomy (TAND Surveillance Agent); and assemble source-attributed therapeutic options briefs (Therapeutics Strategist). Every output carries full provenance and a human review gate; nothing is autonomous or diagnostic. The Engine runs today on a single NVIDIA DGX Spark with elastic RunPod GPUs, demonstrated against a deterministic 50-patient synthetic cohort with explicit evaluation targets. We argue this is a credible, replicable intelligence layer for institutional biobanks and TSC programs, developed against the Cincinnati Children's Hospital Medical Center (CCHMC) and Winslow Research Pavilion opportunity, and designed to generalize across rare disease. The build is the argument.

### Executive Summary

#### The problem, stated plainly

TSC1/TSC2 form a complex that inhibits mTOR, the cellular growth switch. A loss-of-function mutation releases the brake, and hamartomas form across organs: cortical tubers and epilepsy (~85% of patients, two-thirds refractory), subependymal giant cell astrocytomas (SEGA) that threaten the foramen of Monro, renal angiomyolipomas (AML, ~80%, bleeding risk above ~4 cm), pulmonary LAM in adult women, and the pervasive but quietly under-treated TAND phenotype. The 2021 International TSC Consensus (ITSC) surveillance guidelines define a demanding multi-organ, multi-year monitoring cadence. Meeting that standard requires a clinician to integrate genetics, imaging, labs, notes, medications, and the literature at every visit, for every organ, across the patient's life. The evidence exists; the integration does not scale.

#### What the Engine is

The TSC Intelligence Engine is Engine 7 of the HCLS AI Factory, joining Genomic Foundation, Precision Intelligence, Therapeutic Discovery, Clinical Imaging, Precision Oncology, and Cardiology Intelligence under the platform headline "Seven Engines. Eight Intelligence Agents. One Platform." TSC is an engine rather than an agent because it is five coordinated agents plus a deterministic orchestrator, operating cross-modally and coordinating shared state:

| # | Agent | What it does | Model posture |
|---|-------|--------------|---------------|
| 1 | TSC-Variant Curator | BAM to VCF (Parabricks/GATK, mosaic-aware ≥5% VAF), deterministic annotation, evidence aggregation, ACMG-AMP synthesis; recovers the NMI mosaic cohort | Deterministic calling; Sonnet aggregation; Opus classification, validated vs rules |
| 2 | TSC-Phenome Mapper | Time-anchored HPO phenotypes from Clarity-shaped structured data and notes; discordance log; ITSC surveillance-gap report | Sonnet per-note; Haiku normalization; Opus conflict resolution |
| 3 | TSC-Trajectory Modeler | Forecasts SEGA/AML growth, seizure burden, renal function at 6/12/18 months with prediction intervals and threshold alerts | Classical statistics, not LLM-driven; Haiku prose only |
| 4 | TAND Surveillance Agent | Surfaces under-recognized TAND signals across six clusters using the Marshall-Hagedorn diagnostic-uncertainty discourse taxonomy | Sonnet discourse analysis; deterministic scoring; Opus briefing |
| 5 | TSC-Therapeutics Strategist | Six-section options brief integrating all prior agents plus meds, adverse events, PubMed/PMC RAG, ClinicalTrials.gov, FDA | Opus-class, non-negotiable |

A deterministic LangGraph orchestrator (not an LLM) routes 13 event types over an append-only, event-sourced PostgreSQL state with Redis ephemerals, enforcing dependency-ordered enrollment, incremental updates, demand-driven surface assembly, and conservative failure handling. Outputs render to three clinician surfaces: a one-screen pre-visit briefing, a four-quadrant in-visit dashboard, and a disciplined async alert surface. Surfaces are standalone watermarked-synthetic web apps, not Epic.

#### How it runs, and what is and is not built

The Engine runs first on Adam Jones's NVIDIA DGX Spark (GB10 Grace Blackwell, ~1,000 TOPS, 128 GB unified memory) with elastic RunPod GPUs for Parabricks variant calling and parallel cohort generation. It reuses the HCLS AI Factory v1.3.0 substrate (LangGraph, Milvus RAG with a TSC corpus partition, tiered Claude models with a local Llama 3.1 70B fallback, PostgreSQL/Redis/MinIO, audit and provenance). The net-new pieces are the multi-agent orchestrator and shared event-sourced state, the synthetic cohort pipeline, and the five agents themselves.

What runs now is a synthetic-data demonstration. A deterministic, version-controlled 50-patient cohort (30 TSC2 germline, 12 TSC1 germline, 7 mosaic including 1 NMI) is built through a four-layer pipeline: Synthea with TSC modules, BAMSurgeon genomic substrate yielding realistic VCFs, frontier-model clinical notes, and frontier-model imaging reports. The demonstration is evaluated against synthetic ground truth, not validated clinically: the Variant Curator must recover all seven mosaic variants at ≥5% VAF with correct ACMG class and no false-positive Pathogenic calls in under five minutes per case; the Phenome Mapper must hit ≥90% recall and ≥85% precision; the Trajectory Modeler must forecast Patient B's SEGA threshold crossing in the 12-18 month window without false alarms. Epic Clarity/Caboodle, biobank LIMS, and imaging-AI integration are described architecturally and are explicitly not built; they are institutional Phase-1 work. The system is not FDA-cleared, its SaMD posture is undetermined, and every output is decision-support behind a human gate. The platform is Apache 2.0 and non-commercial in intent.

#### The CCHMC and Winslow opportunity

The work is developed against a concrete institutional opportunity at Cincinnati Children's. The Winslow Research Pavilion is the physical envelope; five CCHMC areas are sources that feed the Engine rather than parts of it: the Discover Together Biobank (banked tuber/AML/SEGA tissue, the molecular substrate for mosaic recovery), Biomedical Informatics under Dr. Philip Hagedorn (the TAND methodology, output surfacing, and sponsorship), the TSC clinical and research program (patient concentration), and the Epic Clarity/Caboodle and LIMS data plumbing (Phase-1, not in the demo). The TAND Surveillance Agent is a direct extension of Hagedorn's own published clinical-NLP research on documented diagnostic uncertainty, not an external graft. The governing mantra: a biobank without an intelligence layer is a freezer full of tubes.

#### Why it replicates

The architecture separates the wiring (agents, orchestrator, surfaces, provenance) from the box labels (a specific institution's biobank, informatics group, and clinical program). Replication is "swap the box labels, keep the wiring" to TGen, City of Hope, and other centers, and the same pattern extends to NF1/NF2, Rett, Williams, and other mTORopathies and rare diseases where dispersed multimodal evidence is the core obstacle.

### How to Read This Paper

This paper is written for two audiences at once and asks both to read with a dual lens: the skeptical clinical and informatics reader who needs to trust the outputs, and the engineer who will build them. Where it matters, each part gives both how it is built and runs on Spark plus RunPod, and what it means for a center like CCHMC.

- §1 **Introduction & Thesis** states the clinical problem, the integration gap, and the core claim that follows this front matter.
- The agent chapters specify each of the five agents in build-ready detail: inputs, model tiers, deterministic versus LLM boundaries, outputs, and provenance.
- The orchestration and state chapter covers the event-sourced design and the 13 event types.
- The surfaces chapter defines the pre-visit briefing, in-visit dashboard, and alert surface.
- The synthetic cohort and evaluation chapters define the 50-patient build, the three featured patients, and the explicit eval targets that make claims falsifiable.
- The institutional chapters develop the CCHMC/Winslow opportunity and the replication path.

A reader who wants the thesis can stop after §1. A reader who wants to build can read the agent, orchestration, and cohort chapters. A reader weighing the institutional case can go straight to the CCHMC and replication material. Throughout, the discipline is the same: every claim is sourced, every output is reviewable, and the boundary between what runs today and what is Phase-1 institutional work is marked wherever it is relevant. We turn first to the problem and the thesis.


---

# Part I — Foundations: TSC and the Opportunity


## 1. Introduction & Thesis

Tuberous Sclerosis Complex (TSC) is a disease that medicine understands well and manages poorly. The molecular mechanism is textbook clean: a loss-of-function mutation in *TSC1* or *TSC2* releases the brake on the mTOR growth pathway, and the result is hamartomatous growth across the brain, kidneys, heart, lungs, and skin. We can sequence the gene. We have a targeted therapy that works — mTOR inhibitors shrink subependymal giant cell astrocytomas (SEGA), renal angiomyolipomas (AML), and pulmonary lymphangioleiomyomatosis (LAM), and everolimus became, via the EXIST-3 trial, the first targeted therapy approved for a genetically defined epilepsy. The International TSC Consensus Group has published detailed surveillance guidelines since 2012, updated in 2021. By the standards of rare disease, TSC is an embarrassment of riches.

And yet the care of an individual TSC patient leaks at every seam. The leaks are not failures of knowledge. They are failures of *integration* — of getting the right piece of an enormous, multi-organ, multi-decade clinical picture in front of the right clinician at the right moment. This paper describes an engine built to close those leaks, why we built it for TSC first, why we built it now, and exactly what it does and does not yet do.

### 1.1 The gaps

Four gaps in TSC care motivate this work. Each is well documented in the literature, and each is the kind of problem that an integration layer — not a new drug, not a new scanner — is positioned to address.

**The no-mutation-identified (NMI) cohort.** Roughly 10–15% of patients who meet clinical diagnostic criteria for TSC return a negative result on standard blood-based genetic testing. The dominant explanation is somatic mosaicism: the pathogenic variant is present at low variant allele fraction (VAF), sometimes only in affected tissue, and is simply below the sensitivity of conventional germline assays run on peripheral blood (Tyburczy et al. 2015; Giannikou et al. 2016; Lim et al. 2017). For these families a negative test is not reassurance — it is a dead end that forecloses genetic counseling, cascade testing, and trial eligibility that often requires a molecular diagnosis. Recovering these variants requires deep, mosaic-aware sequencing of the right tissue, and an interpretation pipeline tuned to take a 6% VAF call seriously rather than discard it as noise.

**TSC-Associated Neuropsychiatric Disorders (TAND).** TAND affects roughly 90% of patients across behavioral, psychiatric, intellectual, academic, neuropsychological, and psychosocial domains. The TOSCA natural-history registry found that 30–50% of TAND features go unidentified and unaddressed in routine care. The reason is structural: TAND manifestations surface as offhand mentions in clinical notes — a parent's worry about school, a deferred referral, a behavior described and then not followed up — distributed across years of documentation and across specialties that each see one facet. The TAND consensus framework and the TAND-L lifetime checklist exist precisely because this signal is otherwise lost. No single visit reveals it; only the longitudinal record does.

**Surveillance burden.** A child with TSC accrues a calendar of recurring obligations: serial brain MRI to watch for SEGA growth near the foramen of Monro, renal imaging to track AML size against the ~4 cm bleeding-risk threshold, EEG, ophthalmologic exam, dermatologic and cardiac follow-up, neuropsychological assessment. The ITSC guidelines specify the cadence; real clinics, with real no-shows and competing acute issues, drift from it. Surveillance gaps are common, consequential, and invisible until something crosses a threshold that earlier imaging would have flagged.

**Trial and therapeutic access.** The therapeutic landscape is moving — next-generation selective mTORC1 inhibitors are in development, and ClinicalTrials.gov carries an evolving set of TSC studies — but matching a specific patient's genotype, organ involvement, prior therapy response, and adverse-event history against eligibility criteria and current evidence is a synthesis task that exceeds what a busy specialist can do reliably for every patient at every visit.

None of these gaps is mysterious. Each is a known, published, quantified shortfall. What unites them is that each is an *information-routing* problem layered on top of solid biology and solid guidelines. That is the wedge this engine drives into.

### 1.2 The argument for an integrated engine

The instinct in clinical AI is to build a point tool: a variant classifier, a note summarizer, a growth predictor. We have built point tools too, and they appear here. But the central claim of this paper is that the gaps above are not independent, and a federation of disconnected point tools reproduces the very fragmentation it is meant to cure.

The connections are concrete. A recovered mosaic variant (gap one) changes which surveillance schedule applies and which trials a patient qualifies for (gaps three and four). A longitudinal phenotype profile assembled to surface TAND (gap two) is the same substrate a trajectory model needs to forecast SEGA growth (gap three) and the same evidence base a therapeutics synthesis draws on (gap four). Build these as five strangers and you spend your integration budget gluing them together at the surface, badly. Build them on a shared, event-sourced state with a deterministic orchestrator routing work between them, and integration is the architecture rather than an afterthought.

This is why the system is an **engine** and not an agent. Within the HCLS AI Factory platform an *agent* is a single coordinated capability; an *engine* is a set of agents plus a deterministic orchestrator that is cross-modal and coordinating. The TSC Intelligence Engine comprises five agents — TSC-Variant Curator, TSC-Phenome Mapper, TSC-Trajectory Modeler, TAND Surveillance Agent, and TSC-Therapeutics Strategist — plus a deterministic TSC-Orchestrator and three clinician surfaces. It spans genomic, phenotypic, temporal, and free-text modalities. It coordinates. It is, formally, **Engine 7 of the HCLS AI Factory**, joining the six that precede it: (1) Genomic Foundation, (2) Precision Intelligence, (3) Therapeutic Discovery, (4) Clinical Imaging, (5) Precision Oncology, and (6) Cardiology Intelligence. The platform headline becomes *Seven Engines. Eight Intelligence Agents. One Platform.*

A design principle runs through all five agents and is worth stating before the biology: **the model tier follows the task, and statistics are not outsourced to a language model.** Deterministic work stays deterministic. Variant calling is Parabricks / BWA-MEM / GATK HaplotypeCaller plus Mutect2, mosaic-aware down to ~5% VAF; annotation is snpEff/VEP; these are reproducible pipelines, not prompts. Forecasting in the Trajectory Modeler is *classical* statistics — mixed-effects models, Gaussian process regression, survival analysis, Bayesian updating — emphatically not an LLM dressed up as a forecaster; a language model touches it only to write the prose summary. Where judgment over unstructured evidence is genuinely required, we route by stakes across the Claude tier: Haiku for high-volume normalization (ICD-10 and lab codes to HPO terms), Sonnet for per-note discourse and evidence extraction, Opus for the highest-stakes synthesis — ACMG-AMP classification, rare conflict resolution, and the entire Therapeutics Strategist. Every output carries provenance and a human gate. The engine augments clinicians; it never diagnoses, never recommends autonomously, never replaces a sign-off.

### 1.3 Why TSC, and why now

**Why TSC as the beachhead.** TSC is unusually well suited to be the first disease an integration engine proves itself on, for reasons that are not sentimental but architectural. The genotype-to-phenotype link is mechanistically clear, so a recovered variant has unambiguous clinical meaning. The disease is multi-organ and longitudinal, so it exercises every modality the engine spans — it is a hard enough problem to be a real test, not a toy. The surveillance guidelines are explicit and consensus-backed, giving us ground truth to forecast against. The published gaps (NMI, TAND under-recognition, surveillance drift) are *measured*, so success is measurable. And the literature on each gap — Tyburczy 2015, Giannikou 2016, Lim 2017 on mosaicism; TOSCA on TAND; the TAND consensus framework; ITSC 2021; EXIST-3 — is mature enough to anchor every agent in citable evidence rather than assertion. A disease that is simultaneously tractable, multi-modal, guideline-rich, and gap-quantified is the ideal proving ground for an engine whose whole thesis is integration.

**Why now.** Three things converge. First, frontier language models have crossed the threshold where careful, source-attributed extraction from clinical text and stakes-tiered synthesis are reliable enough to put in front of skeptical clinicians as *draft* material under a human gate — not before. Second, the hardware to run this without a data center now fits on a desk: the engine runs first on an NVIDIA DGX Spark (GB10 Grace Blackwell, ~1,000 TOPS, 128 GB unified LPDDR5x, 4 TB NVMe), with additional GPU capacity drawn from RunPod on demand for GPU-accelerated Parabricks variant calling, parallel synthetic-cohort generation, and heavier local-LLM inference. An engine that a single hospital — or a single architect — can stand up changes the economics of who gets to build precision-medicine infrastructure. Third, the institutional substrate is ready, and that brings us to Cincinnati.

### 1.4 The dual purpose

This work is built with two audiences and two purposes held simultaneously, and the paper does not pretend otherwise.

The first purpose is to **demonstrate real TSC advancements** — to show, on running software, that mosaic variants the standard blood test misses can be recovered and correctly classified, that under-recognized TAND signal can be surfaced as briefing material, that SEGA growth can be forecast across a threshold-crossing window, and that therapeutic options can be synthesized with full attribution. These are engineering claims, and they are demonstrated on a synthetic cohort against synthetic ground truth, with explicit evaluation targets stated later in this paper. This is the build, and the build is the argument.

The second purpose is to **make the institutional case for Cincinnati Children's Hospital Medical Center (CCHMC)** as the place where this engine moves from synthetic demonstration to real-data Phase-1 work. CCHMC is a top pediatric academic medical center — roughly 16,000 staff, more than 600 faculty, over $300M a year in research — with a concentrated TSC clinical and research program. It has, in the Winslow Research Pavilion (opened 2024, ~45,000 sq ft in Avondale), a physical infrastructure envelope housing the Discover Together Biobank and the institution's first centralized biospecimen freezer archive. And it has, in Dr. Philip A. Hagedorn, MD MBI — Chief Health Information Officer and head of the Division of Biomedical Informatics — a published body of clinical-NLP research on detecting diagnostic uncertainty in documentation (Marshall, Nickels, Brady, Hagedorn 2023; Nickels et al. 2024; Ipsaro et al. 2021; Orenstein et al. 2021). The TAND Surveillance Agent's discourse-marker taxonomy is a direct extension of that methodology — an extension of his own published research, not an external graft. Dr. Hagedorn has taken interest in the work and offered to engage his team and the faculty TSC lead.

The architecture of that institutional case is specific and worth stating up front, because it governs the honest-scope posture below. The Winslow Pavilion is the *envelope*; five CCHMC areas are *sources that feed the engine, not parts of it* — the Discover Together Biobank, Biomedical Informatics under Dr. Hagedorn, the TSC clinical and research program, and the Epic Clarity/Caboodle and biobank LIMS systems (the CCHMC chapter details each). The mantra is deliberate: *a biobank without an intelligence layer is a freezer full of tubes.* The engine is the kitchen, not the meal. And because the wiring is institution-agnostic, replication elsewhere — to a TGen or a City of Hope — is a matter of swapping the box labels and keeping the wiring.

### 1.5 The honest-scope posture

Credibility is the entire game with this audience, so the boundary between what is built and what is described must be drawn sharply and kept that way throughout the paper.

What runs **now**, on the Spark with RunPod overflow, is the **synthetic-data engine**: a deterministic, version-controlled 50-patient cohort and the five agents plus orchestrator operating on it end to end. The Epic Clarity / Caboodle integration, the biobank LIMS connection, the real-tissue sequencing, and the imaging-AI pipelines are described *architecturally* and are explicitly **not built**. They are institutional Phase-1 work, contingent on the CCHMC engagement and the IRB approval that real-data work requires; the synthetic demonstration requires no IRB. Wherever this paper describes the institutional data flow, it marks that boundary.

Three further constraints are stated once here and honored everywhere. The engine is **not FDA-cleared**; its Software-as-a-Medical-Device posture is undetermined and is itself institutional work. Every agent output is **decision-support with a human gate** — the Variant Curator produces a draft molecular-genetics report for a board-certified molecular geneticist to sign off, never an autonomous result; the TAND agent surfaces patterns as pre-visit briefing material, never interruptive alerts and never a diagnosis. And the platform is **Apache 2.0, open, and non-commercial in intent**; the clinician-facing surfaces are standalone, persistently watermarked-synthetic web applications, not Epic modules.

These constraints are not apologies. They are the design. An engine that is honest about its boundary is the only kind worth putting in front of clinicians who will, rightly, distrust anything that overclaims. With the thesis, the designation, the dual purpose, and the scope now fixed, the next section turns to the biology the engine is built around — the mTOR mechanism and the specific clinical facts each agent is engineered to act on.


## 2. The Biology, From First Principles

Everything the TSC Intelligence Engine does downstream, recovering a low-frequency variant from tuber tissue, classifying it against ACMG-AMP criteria, forecasting the growth of a subependymal giant cell astrocytoma, reasoning about why an mTOR inhibitor helped one patient and not another, rests on a small number of biological facts that are worth stating precisely before we build anything on top of them. This section lays them out from the bottom up: the cell, its DNA, the proteins that DNA encodes, the ways that DNA can break, the difference between a change you inherit and a change you acquire, the special case of mosaicism, and finally the one signaling pathway and the two genes that turn all of this into Tuberous Sclerosis Complex. Where an analogy helps, we use one; where the analogy would mislead the variant discussion later, we drop it and stay literal.

### 2.1 The cell as the unit of the disease

A human body is built from roughly thirty trillion cells. Each one is a self-contained chemical factory: it takes in nutrients, decides whether conditions are good enough to grow and divide, builds the molecular machinery it needs, and disposes of what it no longer wants. For most of this paper the relevant cells are the ones lining the surface of the brain, the filtering units of the kidney, the smooth muscle of the airways, and the cells of the skin and heart, because those are the tissues where TSC declares itself. But the logic is the same everywhere: a cell that mistakenly believes conditions are always favorable for growth will grow when it should not, and the visible consequence of that mistake is a tumor.

TSC is, at its root, a disease of one such mistake. The cell loses a brake on its own growth decision and behaves as though the green light is permanently on. The tumors of TSC, cortical tubers in the brain, angiomyolipomas (AMLs) in the kidney, lymphangioleiomyomatosis (LAM) in the lung, rhabdomyomas in the heart, angiofibromas in the skin, are called hamartomas: disorganized overgrowths of cells that are native to the organ but present in the wrong amount and the wrong arrangement. They are not, in the usual case, aggressive cancers that spread to distant sites. They are local failures of the growth brake. Understanding the brake is therefore the whole game, and the brake is built from DNA, transcribed into protein, embedded in a pathway. We take those in order.

### 2.2 DNA: the instruction set

Inside almost every cell is a nucleus, and inside the nucleus is the cell's copy of the genome: about three billion base pairs of DNA, organized into 46 chromosomes (23 inherited from each parent). DNA is a long chain written in a four-letter alphabet, the bases A, C, G, and T. The order of those letters is information, and a stretch of DNA that codes for a single protein is a gene. Humans have roughly twenty thousand protein-coding genes spread across that three-billion-letter sequence; the coding portions are a surprisingly small fraction of the whole, interrupted by long non-coding stretches and regulatory regions that decide when and where each gene is read.

Two properties of this arrangement matter for TSC. First, we are diploid: for almost every gene we carry two copies, one on the maternal chromosome and one on the paternal. That redundancy is the reason a single broken copy is often survivable, and it is central to how TSC tumors actually form, as we will see. Second, the sequence is read in three-letter words called codons, each of which specifies one amino acid, the building block of protein. A protein is a chain of amino acids that folds into a specific three-dimensional shape, and the shape is the function. Change the sequence of letters and you can change the amino acid chain; change the chain and you can change or destroy the shape; destroy the shape and you destroy the function. This is the chain of causation that the TSC-Variant Curator agent (Section 6) is ultimately reasoning about when it asks whether a given DNA change matters.

### 2.3 From gene to protein, and the ways the message can break

It helps to fix the central pipeline before describing how it fails. DNA is transcribed into messenger RNA (mRNA), a working copy of one gene. The mRNA is then translated by the ribosome, which walks along it three letters at a time, reads each codon, and adds the corresponding amino acid to a growing protein chain. Three codons are special "stop" signals that tell the ribosome the protein is finished. The protein then folds and goes to work.

A mutation is any change to the DNA sequence. Several distinct kinds matter for TSC, and the engine treats them differently because they have different consequences:

- **Missense.** A single base change that swaps one amino acid for another. The protein is still full length, but one of its building blocks is wrong. The effect ranges from harmless to severe depending on which amino acid, in which part of the protein. These are the hard ones to classify, because the protein is structurally intact and you must reason about whether the substitution actually breaks function.
- **Nonsense.** A single base change that turns an amino-acid codon into a premature "stop." The ribosome quits early and produces a truncated, usually nonfunctional protein. In the synthetic cohort, Patient B carries exactly this: TSC2 c.3037C>T, p.Arg1013Ter, a C-to-T change that converts the codon for arginine at position 1013 into a stop. The protein is cut off roughly halfway through and cannot function.
- **Frameshift.** A small insertion or deletion of bases that is not a multiple of three. Because the ribosome reads in three-letter words, adding or removing one or two letters shifts the reading frame, and every codon downstream is now read in the wrong grouping, producing garbage and, usually, a premature stop soon after. Patient A's variant is a TSC2 frameshift, and the engine treats frameshifts and nonsense changes as presumptively loss-of-function (the basis of the ACMG PVS1 criterion invoked in Section 6).
- **Splice-site.** Genes contain coding segments (exons) interleaved with non-coding segments (introns) that must be cut out of the mRNA before translation. The cut points are marked by specific sequences. A mutation at one of those marks can cause the cell to retain an intron or skip an exon, scrambling the protein.
- **Large deletions and contiguous-gene events.** Sometimes whole exons, whole genes, or stretches spanning multiple genes are deleted. TSC2 sits immediately adjacent to PKD1 on chromosome 16, and a deletion large enough to remove both produces the TSC2/PKD1 contiguous-gene syndrome with severe early polycystic kidney disease. This is why copy-number, not just single-base resolution, has to be on the table.

The practical point for the engine is that "the patient has a TSC2 mutation" is not one thing. The class of mutation drives the interpretation, and the interpretation drives everything downstream. A frameshift and a missense change in the same gene can sit at opposite ends of the pathogenicity scale.

### 2.4 Germline versus somatic: when the change happened

Two cells in your body can carry different DNA, and the reason is timing. You begin as a single fertilized egg. Every cell in the adult body descends from that egg by repeated division, and every division copies the genome. Copying is extraordinarily accurate but not perfect, so changes accumulate.

A **germline** mutation is present in that very first cell, the egg or sperm that formed you, and is therefore copied into every cell of the body, including your own eggs or sperm. It is the same in your blood, your skin, your kidney, your brain. Because it is in your gametes, you can pass it to a child. About one-third of TSC cases are inherited germline variants of this kind. The other two-thirds are *de novo* germline events: a brand-new mutation that arose in a parental gamete or in the earliest divisions of the embryo, present throughout the child's body but absent in both parents. This is why most children with TSC have unaffected parents, and it is the single most important fact about TSC inheritance for genetic counseling.

A **somatic** mutation happens later, in a single body cell at some point after that first division, during development or across a lifetime. It is copied only into that cell's descendants, so it is present in some tissues or some fraction of cells and absent in others. A somatic mutation cannot be passed to a child unless it happens to occur in the germ-cell lineage. Standard genetic testing on a blood sample sees the genome of the blood; a somatic change confined to brain or kidney tissue is, by construction, invisible to it.

This germline-versus-somatic distinction is not a technicality. It determines which tissue you must sequence to find the answer, and it sets up both the second-hit mechanism that builds individual TSC tumors and the mosaicism that defines the most diagnostically difficult patients.

### 2.5 The two-hit mechanism: why one broken copy is not enough

Here the diploid genome and the somatic clock come together. TSC1 and TSC2 are *tumor suppressor* genes: their normal job is to restrain growth (Section 2.7). Tumor suppressors generally follow Knudson's two-hit logic. Because you carry two copies of the gene, losing one copy is often survivable; the remaining copy can still produce enough functional protein to hold the brake. A person with germline TSC carries one broken copy in every cell, the "first hit," present body-wide. That is enough to predispose, but a given cell only loses growth control entirely when the *second*, still-working copy is also knocked out by a somatic event in that particular cell. The descendants of that doubly-hit cell form a tumor.

This explains a pattern that would otherwise be puzzling: a person can carry a single TSC mutation in every cell of the body yet develop discrete tumors in specific locations rather than uniform overgrowth everywhere. Each tumor marks a cell that happened to suffer a second hit. It also explains why the burden of disease is so variable between patients and even between organs in the same patient: the germline hit is fixed, but the somatic second hits are stochastic.

### 2.6 Mosaicism and the no-mutation-identified cohort

Now consider what happens when even the *first* hit is somatic. If the founding mutation arises not in the gamete but in an early embryonic cell, only the descendants of that cell carry it. The result is a mosaic individual: a patchwork of cells, some carrying the TSC variant and some not, with the proportion depending on how early the mutation occurred and which lineages that early cell gave rise to. The fraction of sequencing reads that carry the variant at a given position, the **variant allele fraction (VAF)**, is the readout of this patchwork. A clean germline heterozygous variant sits near 50% VAF, because one of two copies carries it in essentially every cell. A mosaic variant can sit anywhere below that, sometimes far below: a variant present in only a tenth of the cells in a sample appears at roughly 5% VAF.

This is the biological origin of the **no-mutation-identified (NMI)** cohort, the roughly 10 to 15% of patients who meet the clinical diagnostic criteria for TSC but whose standard blood test comes back negative (Tyburczy 2015; Giannikou 2016; Lim 2017). In most of these patients the mutation is real but mosaic, present at low frequency, often concentrated in the affected tissue and scarce or absent in blood. Two things conspire to hide it. First, if the mosaicism is tissue-biased, the blood may simply not contain the variant at appreciable frequency, so you are sequencing the wrong sample. Second, even when the variant is present, it sits at a low VAF that conventional germline-calling workflows, tuned to expect the clean ~50%/100% signatures of inherited variants, are prone to dismiss as sequencing noise.

Patient A in the synthetic cohort is constructed to be exactly this case: a 4-year-old who meets clinical criteria, whose blood testing is negative, and who carries a TSC2 frameshift at 8.3% VAF detectable only in resected tuber tissue. Recovering that variant, distinguishing a real 8.3% signal from artifact, classifying it correctly, and recommending an orthogonal confirmation (droplet digital PCR, ddPCR), is the live demonstration of Act One. The biology says the answer is there in the tissue; the engineering problem is being mosaic-aware enough to see it and disciplined enough not to hallucinate variants that are not there. Both halves of that problem are biological at root, which is why we have spent this long on VAF.

### 2.7 The mTOR pathway: the cell's growth decision

We now have DNA, protein, the kinds of breakage, and the timing of breakage. The remaining question is what TSC1 and TSC2 actually *do*, and the answer is that they sit on top of one of the most important control circuits in cell biology: the mTOR pathway.

mTOR, the mechanistic target of rapamycin, is a protein kinase, an enzyme that switches other proteins on or off by attaching phosphate groups to them. It functions as the cell's master growth controller. When mTOR is active, the cell ramps up the energy-expensive business of growing and dividing: it builds proteins, manufactures lipids for new membranes, and suppresses autophagy (the recycling program a cell runs when starving). When mTOR is quiet, the cell holds back, conserves resources, and waits. mTOR integrates signals about whether growth is a good idea right now: are nutrients available, is there energy, are growth-factor signals arriving from outside, is there enough oxygen. Only when the answers are favorable does it give the order to grow.

The useful mental model is an accelerator pedal. Active mTOR is the foot on the accelerator: grow now. The pathway above it is full of sensors deciding how hard to press. And, crucially, there is a brake.

### 2.8 TSC1 and TSC2 as the brake

The brake is built from the proteins that TSC1 and TSC2 encode, hamartin (TSC1) and tuberin (TSC2). Together with a third partner (TBC1D7) they assemble into the **TSC protein complex**, and the complex's job is to keep mTOR (specifically the mTOR complex 1, mTORC1) switched off until the upstream sensors agree that growth is warranted. Mechanistically, the complex acts on a small switch protein called Rheb: tuberin is the business end that turns Rheb off, and turning Rheb off turns mTORC1 off. So the TSC complex is the foot on the brake, and it presses the brake by default, releasing only when growth signals override it.

Hamartin and tuberin work as an obligate pair. Hamartin stabilizes tuberin and holds the complex together; tuberin carries the enzymatic activity that actually shuts down Rheb. Disable either one and the complex falls apart and stops functioning. This is why mutations in *either* TSC1 *or* TSC2 produce the same disease: they break the same brake from different ends.

Now the full picture assembles. A loss-of-function mutation in TSC1 or TSC2 means a cell can no longer build a working brake. With the brake released, Rheb stays on, mTORC1 stays on, and the cell grows and divides as though the accelerator were stuck down, regardless of whether conditions actually warrant it. The visible result is a hamartoma: native cells of the organ, growing where and when they should not. The two-hit logic of Section 2.5 sets the location: the germline-mutant background predisposes, and a cell that loses its second working copy releases the brake completely and seeds a tumor. The mosaic logic of Section 2.6 sets the diagnostic difficulty: when the founding hit is somatic and tissue-biased, the evidence of the broken brake lives in the affected tissue at low frequency, not in the blood.

Two consequences follow directly, and both anchor later sections. First, because TSC2 is the larger and more complex gene and tends to produce more severe disease, the cohort is weighted toward TSC2 (Section 9 details the 60/24/10/4/2 split across germline TSC2, germline TSC1, mosaic TSC2, mosaic TSC1, and NMI). Second, and more consequentially for therapeutics, the brake-and-accelerator picture is *directly actionable*. If the disease is a stuck accelerator caused by a missing brake, you can intervene downstream by pressing on the accelerator itself. That is precisely what mTOR inhibitors such as everolimus and sirolimus do: they are chemical descendants of rapamycin (the "TOR" in mTOR is named for it), they inhibit mTORC1 directly, and so they substitute pharmacologically for the brake that the mutation removed. This is why a single mechanistic insight ties together SEGA, renal AML, and pulmonary LAM under one drug class, and why understanding the pathway is not optional background but the load-bearing logic of the TSC-Therapeutics Strategist (Section 6).

### 2.9 What this buys us downstream

To summarize the chain we will rely on. Cells decide whether to grow; the mTOR pathway is the decision circuit; TSC1 and TSC2 build the brake on that circuit; loss-of-function mutations in either gene release the brake and produce hamartomas. Whether a mutation is germline or somatic decides which tissue carries the answer and whether a parent can pass it on; the two-hit mechanism decides where individual tumors form; mosaicism, the somatic-first-hit case, produces the low-VAF, blood-negative NMI patients who are the hardest to diagnose and the clearest place for the engine to add value. The class of mutation, missense versus nonsense versus frameshift versus splice versus large deletion, determines how confidently we can call it pathogenic.

None of this yet describes what the disease *does* to a patient across a lifetime: the epilepsy, the SEGAs and their hydrocephalus risk, the renal bleeding threshold, the pulmonary involvement, the under-recognized neuropsychiatric burden, or the surveillance cadence that tries to stay ahead of all of it. That clinical picture, and the consensus guidelines that structure it, is the subject of Section 3.


## 3. The Disease: TSC Across Organ Systems

Tuberous Sclerosis Complex (TSC) is the disease the engine is built around, and almost every design decision downstream traces back to a clinical fact in this section. Before describing what the engine does, it is worth being precise about what the disease does, because the two are tightly coupled: the multi-organ spread, the silent-until-dangerous growth dynamics, the under-recognized neuropsychiatric burden, and the cohort of patients whose mutation hides below the limit of routine testing are not incidental complications. They are the structure of the problem. A reader who understands TSC will recognize each of the five agents introduced later as a direct response to a specific, well-documented gap in how the disease is currently surveilled and managed.

TSC is a single-gene disorder with a multi-system reach that is unusually broad even among the genetic syndromes. It affects roughly **1 in 6,000 births**, which places it among the more common Mendelian disorders a pediatric academic medical center will see, and well within the volume that justifies a dedicated clinical and research program of the kind Cincinnati Children's runs. Inheritance is split: roughly **two-thirds of cases are de novo** (a new mutation in the affected child, with unaffected parents) and roughly **one-third are inherited** from an affected parent, typically in an autosomal-dominant pattern with highly variable expressivity. That variability is the first clinical reality worth holding onto. Two patients carrying the same pathogenic variant, even within the same family, can present with radically different organ involvement and severity. The disease is genetically deterministic at the level of the molecular lesion and stubbornly probabilistic at the level of the patient.

### 3.1 The molecular lesion and why it touches everything

As established in the biology chapter, a loss-of-function mutation in **TSC1** (hamartin) or **TSC2** (tuberin) releases the brake their complex holds on **mTOR**, the cell's central growth switch; the result is uncontrolled growth that manifests as **hamartomas** — benign but space-occupying, often progressive tumors — in essentially any organ the disease reaches. That single mechanism explains the multi-organ phenotype.

Two features of this mechanism shape the rest of the paper. First, because the lesion is a loss of a tumor-suppressor function, TSC follows a **two-hit** logic at the tissue level: a germline first hit plus a somatic second hit in a given cell lineage produces the focal lesions, which is why a patient with a single germline variant develops discrete tubers and tumors rather than uniform tissue change. Second, and clinically decisive, the same pathway that the mutation deranges is **druggable**. mTOR inhibitors (everolimus, sirolimus) act directly on the switch the broken brake fails to control. TSC is therefore one of the relatively rare genetic diseases in which the molecular diagnosis is not merely a label but a direct pointer to mechanism-targeted therapy. That fact is what makes the Therapeutics Strategist (§later) a coherent capability rather than a generic literature-search tool.

### 3.2 Brain

The neurological manifestations dominate the clinical picture for most families, both in severity and in the surveillance attention they demand.

**Cortical tubers** are the lesions the disease is named for: focal malformations of cortical development, present from before birth, that serve as epileptogenic foci. Their number and location correlate loosely with neurological severity but not deterministically, which is part of why imaging-based surveillance is so important and so individualized.

**Epilepsy** is the most common neurological feature, affecting roughly **85% of patients**. The presentation is age-dependent and consequential. In infancy, the characteristic seizure type is **infantile spasms**, a catastrophic epilepsy syndrome where early recognition and treatment materially change developmental outcome. Across the lifespan, roughly **two-thirds of TSC epilepsy is drug-refractory**, meaning seizures persist despite adequate trials of two or more appropriate anti-seizure medications. Refractory epilepsy is the single largest driver of the neurodevelopmental burden, and it is also where the therapeutic landscape has recently shifted (the EXIST-3 trial, discussed below).

**Subependymal nodules** line the ventricles and are usually stable, but a subset grow into **subependymal giant cell astrocytomas (SEGAs)**. SEGAs are the lesion that turns a chronic disease into an acute one. Located at or near the **foramen of Monro**, a growing SEGA can obstruct cerebrospinal-fluid flow and produce **hydrocephalus**, raised intracranial pressure, and, if unrecognized, neurological emergency. The clinical management of SEGA is fundamentally a problem of **growth-rate surveillance**: a lesion measured at 0.8 cm is unremarkable; the same lesion at 1.3 cm and trending toward the foramen is a decision point about whether to start or escalate an mTOR inhibitor versus proceed to neurosurgery. SEGAs typically enlarge on the order of a few millimeters per year, slowly enough that any single MRI looks reassuring and only the **trajectory across serial imaging** carries the signal. This is precisely the gap the Trajectory Modeler is built for, and our featured Patient B — a 12-year-old with a SEGA progressing 0.8 → 1.1 → 1.3 cm at the foramen of Monro — is constructed to exercise it.

### 3.3 Kidney

Renal disease is the second major organ system and, in adults, a leading cause of TSC-related mortality.

**Angiomyolipomas (AMLs)** — benign tumors of vascular, smooth-muscle, and fat tissue — occur in roughly **80% of patients**. They are usually asymptomatic when small, but they carry a size-dependent risk of spontaneous hemorrhage that becomes clinically significant above roughly **4 cm**, where aneurysmal vessels within the lesion can rupture, sometimes catastrophically. As with SEGA, the management question is one of growth and threshold-crossing: an AML tracked across serial renal imaging that is approaching or exceeding 4 cm changes the conversation toward mTOR-inhibitor therapy or selective embolization. Patient C in our synthetic cohort (18-year-old, AML around 4 cm) is positioned at exactly this decision boundary.

**Renal cysts** are also common and usually benign. The clinically important subtlety is the **TSC2/PKD1 contiguous gene deletion syndrome**: TSC2 and PKD1 (the major polycystic-kidney-disease gene) are physically adjacent on chromosome 16, and in roughly **1–2% of patients** a deletion removes both. These patients present with severe, early-onset polycystic kidney disease in addition to the standard TSC phenotype, and their renal prognosis is materially worse. This is a small but important fraction precisely because it is a genotype that changes surveillance and counseling, and it is the kind of finding the Variant Curator must represent faithfully rather than collapse into a generic "TSC2 pathogenic" call.

### 3.4 Lung

**Lymphangioleiomyomatosis (LAM)** is a progressive cystic lung disease that occurs almost exclusively in **adult women** with TSC. It is frequently asymptomatic for years and then presents with dyspnea or pneumothorax. LAM matters for surveillance design because it is a manifestation with a defined at-risk subpopulation (post-pubertal females) and a defined screening cadence in the consensus guidelines, and because it is one of the manifestations for which mTOR inhibitors are an established, evidence-supported therapy. It is a clean example of the more general TSC pattern: a manifestation that is invisible until it is dangerous, with a known at-risk group and a known intervention, where the failure mode is not absence of knowledge but absence of timely, individualized attention.

### 3.5 Heart

**Cardiac rhabdomyomas** are benign muscle tumors that are often the earliest detectable sign of TSC — frequently found on prenatal or neonatal echocardiography, sometimes prompting the initial diagnosis. Their clinical course is reassuring in most cases: they **typically regress** over the first years of life. Their importance is diagnostic and temporal rather than chronic-management: a rhabdomyoma on a fetal echo is often the first thread in the diagnostic story, and it anchors the timeline of a patient's TSC history. For the Phenome Mapper, a regressing rhabdomyoma is a good test of time-anchored phenotype extraction, where the same finding is salient at birth and unremarkable at age five.

### 3.6 Skin

Cutaneous manifestations are nearly universal and are heavily weighted in the clinical diagnostic criteria because they are visible and specific. They include **hypomelanotic macules** (ash-leaf spots, often the earliest visible sign, best seen under Wood's lamp), **facial angiofibromas**, **shagreen patches** (connective-tissue nevi, usually lumbosacral), and **ungual fibromas** (periungual growths that tend to appear later, in adolescence and adulthood). Dermatological findings rarely threaten life, but they are central to clinical diagnosis and they are richly documented in clinical notes, which makes them valuable, well-attested phenotype anchors for the Phenome Mapper to extract and time-stamp.

### 3.7 TAND: the burden that gets missed

The manifestations above are the ones medicine surveils well. The **TSC-Associated Neuropsychiatric Disorders (TAND)** are the ones it does not, and they are the disease's largest unaddressed burden. TAND is an umbrella term, formalized by international consensus, for the full range of behavioral, psychiatric, intellectual, academic, neuropsychological, and psychosocial difficulties associated with TSC. The consensus framework organizes TAND into **six clusters**:

1. **Behavioral** — overactivity, aggression, self-injury, sleep disturbance, mood difficulties.
2. **Psychiatric** — formally diagnosable disorders, prominently autism spectrum disorder, ADHD, anxiety, and depression.
3. **Intellectual** — the global cognitive-ability profile, ranging from unaffected to profound intellectual disability.
4. **Academic** — specific difficulties with reading, writing, spelling, and mathematics.
5. **Neuropsychological** — domain-specific deficits in attention, memory, executive function, and language.
6. **Psychosocial** — self-esteem, family stress, relationship and quality-of-life impacts that extend to caregivers and siblings.

Two numbers define the problem. TAND features affect roughly **90% of patients** across the lifespan — TAND is not a complication of a subset but a near-universal dimension of the disease. And yet, per the **TOSCA** international registry, an estimated **30–50% of TAND features go unidentified or unaddressed** in routine care. The gap is not a knowledge gap in the literature; the natural history and the recommendation to screen are both well established. It is a gap in **recognition and follow-through at the point of care**. TAND features surface in clinic visits as soft signals — a parent's offhand remark about sleep, a teacher's note relayed third-hand, a clinician's hedged observation that a behavior "might be worth watching" — and these signals routinely fail to crystallize into assessment, referral, or a documented plan.

The consensus response to this is the **TAND Checklist (TAND-L)**, a lifetime screening instrument designed to be completed periodically across a patient's life to systematically surface features that would otherwise be missed. The TAND-L is the right clinical instrument; the problem is that completing it consistently, across a busy multi-specialty TSC clinic with competing acute concerns, is exactly the kind of structured-attention task that erodes under real-world clinic pressure.

This is the precise seam the TAND Surveillance Agent (introduced later) is designed to sit in, and it is worth being explicit about *why* the seam is tractable. The unaddressed TAND signal is, very often, **already documented** — it exists in the longitudinal clinical notes as hedged, deferred, or third-party-attributed language that never escalated to a coded diagnosis or a referral. Detecting that pattern is a clinical-NLP problem about **diagnostic uncertainty in documentation**, which is exactly the research program Dr. Philip Hagedorn's group at Cincinnati Children's has published on (Marshall, Nickels, Brady, Hagedorn 2023; Nickels, Marshall, Edgerton, Brady, Hagedorn, Lee 2024). The agent applies that discourse-marker methodology to the six TAND clusters. The disease fact that makes this honest rather than speculative is that **the signal is there to be found** — the failure is downstream of documentation, in recognition and action, not upstream in observation. We return to this in detail later; here the point is only that TAND, the largest missed burden in TSC, is missed in a way that leaves a documentary trace.

### 3.8 The NMI cohort and the diagnostic odyssey

TSC can be diagnosed clinically, on the basis of major and minor features, without a genetic confirmation. But molecular confirmation matters: it sharpens prognosis, enables cascade testing of family members, qualifies patients for genotype-specific surveillance and trials, and resolves diagnostic ambiguity in atypical presentations. The problem is that an estimated **10–15% of patients who meet clinical diagnostic criteria test negative on standard genetic testing** — the **no-mutation-identified (NMI)** cohort.

In the large majority of these cases, the explanation is **somatic mosaicism**. The pathogenic variant arose after fertilization and is present in only a fraction of the patient's cells — and critically, often at a **low variant allele fraction (VAF)** in the tissues that routine testing samples. Standard clinical genetic testing is performed on **blood**, and a mosaic variant that is abundant in affected tissue (a resected tuber, an AML, a SEGA specimen) may be present in blood at a fraction below the effective limit of detection of standard pipelines, which are tuned for germline variants near 50% VAF. The variant is real; the assay simply cannot see it where it looks. This mechanism is well documented in the somatic-mosaicism literature (Tyburczy 2015; Giannikou 2016; Lim 2017), which established that deep sequencing of **affected tissue** recovers pathogenic variants in a substantial fraction of NMI patients, frequently at VAFs in the single-digit-percent range.

For the family, an NMI result is the start of a **diagnostic odyssey**. The clinical picture says TSC; the genetic test says nothing found. That discordance undermines confidence in the diagnosis, complicates eligibility for genotype-anchored care, blocks informative family testing, and leaves a patient and family in an evidentiary limbo that can persist for years. The resolution — deep tissue sequencing with a mosaic-aware variant-calling pipeline capable of confidently distinguishing a true low-VAF variant from sequencing artifact — is technically achievable but is not routinely deployed, in part because the analytic burden of validating a single-digit-VAF call is real.

This is the gap the **TSC-Variant Curator** is built to close, and it is the centerpiece of the demonstration's opening act. Our synthetic cohort deliberately includes **seven mosaic patients** out of fifty — five TSC2 mosaic, two TSC1 mosaic — plus one true NMI patient, and the featured **Patient A** is a 4-year-old whose pathogenic TSC2 frameshift sits at an **8.3% VAF** in tuber tissue: invisible to a blood test, recoverable from tissue, and classifiable to **Likely Pathogenic** (PVS1 + PM2 + PP4) with a recommendation for orthogonal ddPCR validation. The disease fact here is exact and unforgiving: if the pipeline's effective sensitivity floor is above the VAF where these variants actually live, the entire capability is theater. That is why the engineering target is stated as recovery of every mosaic variant at **VAF ≥ 5%** with **no false-positive Pathogenic calls**, and why the molecular-geneticist sign-off gate is non-negotiable — a low-VAF call carries enough consequence that it must be a draft for review, never an autonomous output.

### 3.9 Surveillance burden and the family's lived reality

The picture that emerges from the organ-by-organ tour is of a disease that is **multi-system, lifelong, and dominated by slow-moving lesions that are dangerous only at thresholds and only in trajectory**. The clinical response, codified in the **2021 International TSC Consensus surveillance guidelines (ITSC)**, is a schedule of periodic, organ-specific monitoring: serial brain MRI for SEGA, renal imaging for AML and cysts, pulmonary screening for at-risk women, neurodevelopmental and TAND assessment, dermatological and ophthalmological review, and seizure management — each on its own recommended cadence, often under a different subspecialist.

For a single patient and family, this translates into a relentless logistical and emotional load: multiple specialists, multiple imaging studies on staggered schedules, the cognitive burden of integrating findings that no single clinician sees in full, and the chronic background anxiety of waiting to learn whether this year's MRI shows the SEGA has grown. The information that matters most is rarely contained in any single encounter. It lives in the **comparison across encounters** — the growth rate, the new soft TAND signal that echoes one from two visits ago, the surveillance interval that quietly lapsed. And that cross-encounter, longitudinal integration is exactly what the human system does worst under load, and what is most often dropped: the **ITSC surveillance-gap** where a recommended study simply did not get ordered, the TAND signal that never closed the loop, the trajectory that was reassuring on each individual scan and alarming only in aggregate.

These are not failures of knowledge or of clinician diligence. They are failures of **attention at scale** — the predictable consequence of asking busy humans to hold a multi-year, multi-organ, multi-specialist picture in their heads across fragmented encounters. That reframing is the bridge to the rest of this paper. Each of the gaps named in this section — the low-VAF variant below the assay floor, the SEGA dangerous only in trajectory, the near-universal but half-missed TAND burden, the lapsed surveillance interval — is a place where a disciplined, provenance-bearing, decision-support layer can recover signal that already exists in the data without ever displacing the clinician's judgment. The next section turns from the disease to that opportunity: what specifically becomes possible when these gaps are treated not as inevitable but as addressable, and why TSC at an institution like Cincinnati Children's is the right place to prove it.


## 4. The Clinical & Scientific Opportunity

The case for the TSC Intelligence Engine does not rest on the novelty of any single technique. It rests on a set of gaps in Tuberous Sclerosis Complex (TSC) care that are already documented in the literature, already named in consensus guidelines, and already costing patients diagnoses and timely interventions. Each gap is a place where information that exists somewhere in the record, the tissue, or the trajectory fails to reach the clinician at the moment a decision is made. The Engine targets five of these gaps. None of them is hypothetical. The contribution is not "AI can do TSC"; it is that the right model tier, applied to the right substrate, with a human gate and full provenance, can surface signal that is currently lost to fragmentation, under-recognition, and timing.

This section frames the five gaps as documented unmet needs and ties each to a plausible mechanism of impact. It deliberately stops short of claiming clinical validation: the demonstration described later in this paper runs against a synthetic cohort with known ground truth, which establishes that the pipeline behaves correctly, not that it improves outcomes in patients. The opportunity is real and well-attested; the evidence that the Engine closes it is the institutional Phase-1 work, not the demo.

### 4.1 Mosaic recovery: the no-mutation-identified cohort

As the disease chapter establishes, roughly 10–15% of patients who meet the clinical diagnostic criteria for TSC have negative conventional genetic testing — the no-mutation-identified (NMI) cohort — most often because of somatic mosaicism: a *TSC1* or *TSC2* variant present in only a fraction of cells, often below the detection floor of standard blood-based panels designed around the ~50% variant allele fraction (VAF) expected of a germline heterozygous variant. A blood draw samples the wrong tissue at the wrong sensitivity, and the causal variant is simply not seen.

The mosaic literature is now specific enough to design against. Tyburczy et al. (2015) demonstrated that deep, targeted next-generation sequencing recovers low-level mosaic and intronic *TSC1/TSC2* variants in a substantial fraction of previously NMI patients, with mutant allele fractions reaching into the low single-digit percentages. Giannikou et al. (2016) characterized the somatic mutational landscape of TSC-associated tumors, showing that lesional tissue — cortical tubers, angiomyolipomas (AML), subependymal giant cell astrocytomas (SEGA) — frequently carries the causal hit at far higher VAF than blood. Lim et al. (2017) extended this to brain mosaicism, recovering pathogenic variants from resected tuber tissue in patients whose blood was uninformative. The through-line is consistent: the variant is recoverable if you sequence the right tissue deeply and call variants in a way that does not silently discard low-VAF reads as noise.

This is precisely where standard germline-tuned variant-calling fails. A pipeline configured for ~50% VAF germline heterozygotes will treat a 6–8% VAF mosaic allele as sequencing artifact and filter it out. The unmet need is therefore not "more sequencing" in the abstract; it is mosaic-aware calling — the somatic callers, low-VAF threshold, and evidence aggregation detailed in the Variant Curator chapter, with orthogonal confirmation by droplet digital PCR (ddPCR).

The mechanism of impact is concrete and bounded. For an NMI patient, recovering a pathogenic mosaic variant converts a clinical-only diagnosis into a molecularly confirmed one. That conversion changes what happens next: it ends the diagnostic odyssey for that family, it enables informed reproductive and cascade-testing conversations (a mosaic variant in proband tissue carries very different recurrence implications than an inherited germline variant), and it can anchor surveillance and therapy decisions in the genotype. A molecular geneticist still owns the call — the Variant Curator produces an AI-labeled draft interpretation with full provenance, never an autonomous classification — but the draft puts a recoverable variant in front of a human who would otherwise never have seen it. Featured Patient A in the synthetic cohort makes this tangible: a 4-year-old with an 8.3% VAF *TSC2* frameshift in tuber tissue, classified Likely Pathogenic under ACMG-AMP (PVS1 + PM2 + PP4), with a ddPCR confirmation recommendation attached. The recovery is the entire point of Act One.

### 4.2 TAND surfacing: the most common, least addressed manifestation

TSC-Associated Neuropsychiatric Disorders (TAND) are the manifestation patients and families most often report as the hardest to live with, and they are also the most systematically missed. As detailed in the disease chapter, TAND spans six clusters, affects approximately 90% of patients across the lifespan, and goes unrecognized or unaddressed in routine care on the order of 30–50% of the time — the gap the TOSCA natural-history registry quantified and the TAND framework and Checklist (TAND-L) were created to close.

The reason TAND is under-recognized is structural, not negligent. TSC clinics are organized around the manifestations with imaging and laboratory surrogates — SEGA growth, AML size, renal function, seizure frequency. TAND signal, by contrast, lives in clinical narrative: a parent's offhand remark about school, a note that a behavioral concern was raised but deferred, an attention problem mentioned and then never formalized into an assessment. These signals are diffuse, longitudinal, and rarely coded. No single note flags TAND; the pattern emerges only across many notes, and only if someone reads them as a series with TAND in mind. In a 20-minute visit anchored on an MRI result, that reading rarely happens.

This is the gap the TAND Surveillance Agent addresses, and it is the gap where this project has its strongest scientific footing. The agent applies the Marshall–Hagedorn diagnostic-uncertainty discourse-marker taxonomy — hedging, deferral, third-party attribution, conditional framing, and follow-up-without-formalization — developed in clinical-NLP research at Cincinnati Children's (Marshall, Nickels, Brady, Hagedorn 2023, *Journal of Hospital Medicine*; Nickels, Marshall, Edgerton, Brady, Hagedorn, Lee 2024, *Applied Linguistics*; Ipsaro, Patel, Marshall, Hagedorn 2021, *Hospital Pediatrics*). Those markers are exactly the linguistic fingerprints of an under-recognized concern: a clinician who writes "mother reports some difficulty at school, will monitor" has left a deferral-and-third-party-attribution trace that, repeated across visits, is a TAND signal. Detecting that pattern is not a graft of generic sentiment analysis onto TSC; it is a direct extension of published methodology into a disease where the consensus literature already says the under-recognition exists.

The mechanism of impact is recognition, not diagnosis. The agent surfaces candidate TAND patterns as pre-visit briefing material — never an interruptive alert, never a diagnosis — so that a clinician walks into the room already aware that the academic and behavioral clusters have accumulated signal worth a direct conversation and, where indicated, a formal TAND-L review or referral. Recognition is the precondition for everything downstream; the TOSCA gap is fundamentally a recognition gap, and recognition is what an attentive longitudinal read of the narrative provides. Patient B in the synthetic cohort carries scattered, under-recognized TAND signals embedded across his notes for exactly this reason: the demo's eval target is that the agent detects the embedded signals without raising spurious flags.

### 4.3 Trajectory forecasting: making surveillance anticipatory

TSC surveillance is currently reactive by design. The international consensus guidelines (Northrup et al., 2021; the ITSC 2021 surveillance recommendations) specify what to monitor and at what baseline cadence — brain MRI for SEGA, renal imaging for AML, and so on — but the cadence is fixed by guideline rather than tuned to the individual's measured trajectory. A patient whose SEGA is growing steadily and a patient whose SEGA is stable receive, by default, the same scan interval. The clinically consequential events — a SEGA approaching the foramen of Monro and the threat of obstructive hydrocephalus, an AML crossing the ~4 cm bleeding-risk threshold, a measurable decline in renal function — are events the guidelines tell us to watch for but do not help us anticipate from a given patient's own longitudinal measurements.

The unmet need is anticipation: given a patient's serial measurements, when is a clinically meaningful threshold likely to be crossed, and is the current surveillance cadence adequate to catch it in time? This is a classical statistical estimation problem, and the Engine treats it as one. The Trajectory Modeler is deliberately not LLM-driven; it uses mixed-effects models, Gaussian process regression, survival analysis, and Bayesian updating to forecast SEGA and AML growth, seizure burden, and renal function at 6/12/18 months with explicit 50% and 90% prediction intervals, then issues threshold-crossing alerts and surveillance-cadence recommendations. A language model writing a plausible-sounding number here would be a liability; a calibrated interval with an honest uncertainty band is the correct tool, and the use of a large model is confined to prose summary and the interpretation of genuinely unusual trajectories.

The mechanism of impact is better-timed intervention. SEGA, in particular, is a manifestation where timing is everything: a lesion caught while it is small and growing is a candidate for elective mTOR-inhibitor therapy or planned resection; a lesion that obstructs cerebrospinal fluid flow before the next scheduled scan is a neurosurgical emergency. A forecast that Patient B's SEGA — measured at 0.8 → 1.1 → 1.3 cm at the foramen of Monro — is likely to cross an actionable threshold within a 12–18-month window, delivered with a calibrated interval and a recommendation to shorten the imaging interval, converts a fixed-cadence schedule into an individualized one. The demo eval target is exactly that: forecast the Patient B crossing within the window, with no false alarms. Anticipation that cannot be trusted is worse than none, which is why calibration and false-alarm suppression are first-class requirements rather than afterthoughts.

### 4.4 Trial matching and emerging therapeutics: keeping pace with a moving field

TSC therapeutics have moved from supportive care to targeted intervention within a clinical generation, and the pace is the problem. mTOR inhibitors (everolimus, sirolimus) are now standard for SEGA, AML, and lymphangioleiomyomatosis (LAM). The EXIST-3 trial established everolimus as adjunctive therapy for TSC-associated refractory seizures — notable as one of the first targeted therapies approved in a genetically defined epilepsy. Next-generation selective mTORC1 inhibitors are in development, and the trial landscape on ClinicalTrials.gov, the emerging evidence in PubMed/PMC, and FDA actions all shift on a timescale faster than any individual clinician can track across a full TSC panel of organ systems and comorbidities.

The unmet need is synthesis at the point of decision. A TSC patient on a partial everolimus response with a dose-limiting toxicity, an AML near the intervention threshold, and refractory seizures presents a genuinely multi-dimensional therapeutic question — and the relevant trials, the relevant emerging evidence, and the relevant safety signals are scattered across registries and literature that no clinic has time to canonically re-search at every visit. Patient C in the synthetic cohort embodies this: an 18-year-old with a *TSC1* variant, a partial everolimus response complicated by mucositis requiring dose reduction, an AML around 4 cm, and refractory focal seizures.

The Engine addresses this with the Therapeutics Strategist, the one agent where the model tier is non-negotiably Opus-class because the reasoning integrates all four prior agents plus medications and adherence, adverse events, a PubMed/PMC retrieval-augmented corpus, a ClinicalTrials.gov snapshot, and FDA actions. Its output is a six-section structured options brief — Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions — in which every claim is source-attributed and the framing is explicitly decision-support, not recommendation. The mechanism of impact is twofold: surfacing trial-eligibility matches the clinic would otherwise have to find by hand, and assembling the current evidence picture so the human decision is made against the field as it stands rather than as it stood at the last literature review. The discipline that makes this credible is the attribution requirement and the appropriately hedged language; the demo eval targets correct trial matches, appropriate hedging, and full source attribution under a three-minute latency budget.

### 4.5 Surveillance adherence: the gap between guideline and practice

Beneath the four clinical gaps sits a fifth that conditions all of them: the distance between what the ITSC 2021 guidelines recommend and what actually happens in a busy clinic across a multi-system disease. The guidelines are comprehensive precisely because TSC touches the brain, kidneys, lungs, heart, skin, and neuropsychiatric domains, and that comprehensiveness is hard to execute patient-by-patient when the recommended surveillance is distributed across specialties, intervals, and modalities. Surveillance items get missed not through inattention but through the sheer combinatorial load of tracking, for each patient, which recommended assessment is due, which is overdue, and which manifestation has changed enough to warrant deviation from the default cadence.

The unmet need here is bookkeeping made visible: a per-patient view of the recommended surveillance set against what has actually been done, with the gaps surfaced. This is the connective tissue the Phenome Mapper provides as a byproduct of building the longitudinal HPO-coded phenotype profile — its surveillance-gap report (mapped against the ITSC recommendations) tells a clinician, before the visit, which guideline-recommended assessments have lapsed. The mechanism of impact is the most prosaic of the five and arguably the highest-yield: a missed renal ultrasound or a lapsed MRI interval is a preventable failure, and surfacing it costs nothing but attention. Crucially, this gap is also where the strict "do not overclaim" posture bites hardest — the surveillance-gap report in the demonstration is computed over synthetic structured data; against real Epic Clarity/Caboodle and biobank LIMS feeds it is institutional Phase-1 work, not something the synthetic demo builds.

### 4.6 Why these five, and why now

The five gaps are not independent; they compound. Mosaic recovery anchors surveillance and therapy in genotype. The phenotype profile that drives surveillance-gap detection is the same substrate the trajectory and TAND analyses build on. Trajectory forecasts inform the therapeutic question the Strategist synthesizes. Each agent is individually defensible against the literature, and the coordination across them — handled by a deterministic orchestrator rather than a model — is what makes this an Engine rather than five disconnected tools. The "now" is supplied by three convergences: a mosaic literature specific enough to design calling pipelines against; a published discourse-marker methodology that turns TAND under-recognition into a tractable NLP problem; and a hardware substrate — a single DGX Spark with elastic RunPod GPUs — on which the whole pipeline can be built and run by one person.

What follows is how that pipeline is actually constructed. Section 5 turns from the opportunity to the platform: the reused HCLS AI Factory v1.3.0 substrate, the net-new multi-agent orchestrator and event-sourced state, and the synthetic-cohort pipeline that lets every claim in this section be tested against known ground truth before any of it touches a patient.


---

# Part II — The Platform and Engine 7


## 5. The HCLS AI Factory Platform

The TSC Intelligence Engine does not stand alone. It is the seventh engine of the HCLS AI Factory, an open-source (Apache 2.0) precision-medicine platform whose substrate it reuses wholesale: the same orchestration library, the same retrieval stack, the same tiered-model routing, the same audit and provenance discipline, the same data services, the same hardware target. This section describes the platform on its own terms — what it is, how it is layered, and what runs where — so that §6 can do the narrower work of showing why five coordinated TSC agents plus a deterministic orchestrator constitute an engine rather than an agent.

The discipline that matters here is honesty about the boundary between what is built and what is described. The HCLS AI Factory v1.3.0 substrate is real, exercised, and version-controlled. Where this platform reaches institutional scale — the platform storage, Epic Clarity plumbing, hospital-grade deployment — those are described as the envelope the engine is designed to drop into, not as things the synthetic-data demo exercises. The build is the argument, and the argument only holds if the line between the two is drawn plainly.

### 5.1 The hierarchy: Platform → Engine → Agent → Tool

The platform organizes capability into four nested tiers. The vocabulary is load-bearing, because the distinction between an engine and an agent is exactly what §6 turns on.

| Tier | Definition | Example |
|------|------------|---------|
| **Platform** | The whole: shared substrate, conventions, hardware target, license. One thing. | HCLS AI Factory v1.3.0 |
| **Engine** | A coordinated assembly — multiple agents plus a deterministic orchestrator — that is cross-modal and produces an integrated work product across a domain. | TSC Intelligence Engine (Engine 7) |
| **Agent** | A single bounded reasoning unit: one input contract, one output contract, one job. May call multiple tools and multiple model tiers, but it answers one question. | TAND Surveillance Agent |
| **Tool** | A deterministic capability an agent invokes: a caller, an annotator, a retriever, a statistical routine, a database write. No autonomy of its own. | Parabricks variant calling; snpEff/VEP annotation; a Milvus similarity query |

The two facts that separate an engine from an agent are coordination and cross-modality. An agent has one input contract and one output contract; it does one job, however sophisticated. An engine coordinates several such agents through a deterministic state layer, draws on more than one data modality (genomic, structured clinical, free-text, longitudinal), and assembles their outputs into a single coherent work product. The TSC Intelligence Engine has five agents and an orchestrator and spans BAM-derived variants, HPO-coded phenotypes, longitudinal trajectories, free-text TAND signals, and a literature-grounded therapeutics brief. That is an engine. A pharmacogenomics agent that maps a star-allele genotype to a dosing recommendation is an agent. The platform contains both kinds of thing, and names them differently on purpose.

### 5.2 Seven Engines. Eight Intelligence Agents. One Platform.

The platform headline, as of v1.3.0, is **Seven Engines. Eight Intelligence Agents. One Platform.** The engines are the coordinated, cross-modal assemblies; the agents are the standalone bounded units. The TSC Intelligence Engine is the seventh engine — the first to be specified disease-first and the first built around a named institutional methodology (Marshall-Hagedorn, §6).

The seven engines:

1. **Genomic Foundation** — primary genomic processing: alignment and variant calling from sequence to annotated variants, the deterministic floor the rest of the platform stands on.
2. **Precision Intelligence** — variant interpretation against curated knowledge (ClinVar, gnomAD, AlphaMissense), evidence aggregation, and ACMG-style classification synthesis.
3. **Therapeutic Discovery** — generative molecule design and structure-based docking for target-to-candidate workflows.
4. **Clinical Imaging** — image-derived measurement and longitudinal change detection across radiology modalities.
5. **Precision Oncology** — tumor-context interpretation, biomarker integration, and oncology decision-support assembly.
6. **Cardiology Intelligence** — cardiac phenotype and risk integration across structured and signal data.
7. **TSC Intelligence Engine** — the subject of this paper; five coordinated agents plus a deterministic orchestrator for Tuberous Sclerosis Complex.

The eight intelligence agents — bounded, single-contract units that run standalone:

1. **CAR-T** — cell-therapy candidacy and construct context.
2. **Precision Biomarker** — biomarker extraction and panel-level interpretation.
3. **Pharmacogenomics** — star-allele genotype to dosing and interaction guidance.
4. **Precision Autoimmune** — autoimmune phenotype and serology integration.
5. **Neurology** — neurological phenotype interpretation.
6. **Clinical Trial** — trial eligibility matching against structured criteria.
7. **Rare Disease Diagnostic** — phenotype-driven differential narrowing for undiagnosed disease.
8. **Single-Cell** — single-cell expression interpretation.

Several TSC-engine capabilities rhyme with platform agents — the Variant Curator shares lineage with Precision Intelligence's classification logic; the Phenome Mapper's HPO work resembles the Rare Disease Diagnostic agent; the Therapeutics Strategist's trial matching parallels the Clinical Trial agent. That is the point of a shared substrate. The TSC engine does not reinvent these capabilities; it specializes, coordinates, and cross-references them under one deterministic orchestrator, which is what makes the assembly an engine rather than a basket of agents.

### 5.3 Model tiering: route to the cheapest tier that is correct

The platform does not use one model for everything. It routes each unit of work to the least expensive tier that is demonstrably correct for that work, reserving the most capable (and most expensive) model for the steps where judgment under uncertainty genuinely matters. Four tiers are in play.

| Tier | Model | Where it is used |
|------|-------|------------------|
| **Fast** | Claude Haiku | High-volume normalization and prose summary: ICD-10/lab-to-HPO mapping (Phenome Mapper), trajectory prose summaries (Trajectory Modeler). |
| **Analysis** | Claude Sonnet | Per-note extraction and discourse analysis: per-note HPO extraction, TAND discourse-marker analysis across the six clusters, unusual-trajectory interpretation, variant evidence aggregation. |
| **High** | Claude Opus | Judgment under uncertainty where the answer must be defensible: ACMG-AMP classification synthesis (validated against combinatorial rules), rare phenotype conflict resolution, TAND briefing summaries, and the entire six-section Therapeutics Strategist brief (Opus-class, non-negotiable). |
| **Local fallback** | Llama 3.1 70B Instruct (Ollama) | On-Spark inference when API access is unavailable or when a step must run fully local; correctness-checked against the API tiers, not a silent substitute. |

Two design rules govern tier assignment. First, statistics is not an LLM job: the Trajectory Modeler's forecasting — mixed-effects models, Gaussian process regression, survival analysis, Bayesian updating — runs as classical statistics, with Haiku confined to turning the numbers into a sentence and Sonnet reserved for flagging a trajectory that does not fit the model. Second, the highest tier is gated to the steps that carry clinical or interpretive risk, because an Opus call costs an order of magnitude more than a Haiku call and the platform processes a cohort, not a single case. The cost ledger is per-tier and per-call, so the economics of a cohort run are auditable rather than assumed.

### 5.4 Retrieval: Milvus, biomedical embeddings, corpus partitions

Every LLM step that makes a literature- or knowledge-grounded claim retrieves its evidence rather than recalling it from weights. The platform's retrieval-augmented generation (RAG) stack is built on **Milvus 2.4** as the vector store, with a two-embedding strategy: **BAAI/bge-large-en-v1.5** for general biomedical text and a **BiomedBERT-derived clinical embedding** for clinical-language passages, so that note-shaped text and literature-shaped text are each embedded by a model tuned for it.

The corpus is partitioned, not pooled. A TSC corpus partition holds the literature anchors the engine reasons over — the ITSC 2021 consensus surveillance guidelines, the somatic-mosaicism literature (Tyburczy 2015, Giannikou 2016, Lim 2017), the EXIST-3 trial and the broader mTOR-inhibitor evidence base, the TOSCA registry, and the TAND consensus framework — alongside a ClinicalTrials.gov snapshot and an FDA-actions snapshot. Partitioning matters for two reasons: it keeps the Therapeutics Strategist's retrievals scoped to relevant evidence rather than the whole platform corpus, and it lets each retrieved passage carry a stable URI back to its source. That URI is what makes the provenance record (§5.6) checkable. When the Therapeutics Strategist asserts that everolimus is indicated for SEGA and AML, the assertion is attached to the retrieved passage that supports it, by URI, in the brief itself.

### 5.5 Deterministic orchestration: LangGraph, not an LLM

The platform's orchestration substrate is **LangGraph**, and the single most important architectural decision in the TSC engine is that the orchestrator is deterministic. The TSC-Orchestrator is an event router, not a reasoning agent: it has no LLM in its decision path. It consumes 13 event types and applies fixed, inspectable rules — dependency-ordered enrollment (the Phenome Mapper runs first, because the other agents build on its longitudinal HPO profile), incremental-update minimization (a new note re-runs only the agents that depend on it), demand-driven surface assembly (a surface is composed when a clinician asks for it), and conservative failure handling (a failed agent leaves its surface region showing "pending" or a staleness marker, never a silently missing output).

This is a deliberate inversion of the common "agentic" pattern in which a planner LLM decides what to do next. A planner LLM is non-deterministic, hard to audit, and capable of surprising you in a clinical setting — exactly the properties you do not want in the component that decides which agent sees which patient's data. Putting deterministic LangGraph in that seat means the control flow is a finite, testable state machine; the LLMs do the reasoning inside agents, where their outputs are bounded by an input and output contract and gated by a human.

State is event-sourced. PostgreSQL holds an append-only event log plus materialized current-state projections; Redis holds ephemeral working state with a TTL; a YAML demo configuration parameterizes the run. Event-sourcing means the engine's entire history is replayable: any surface a clinician saw can be reconstructed from the event log, which is both an audit property and a debugging property. The orchestrator is net-new for TSC — it is one of the three things the engine adds to the substrate, alongside the synthetic-cohort pipeline and the five agents themselves — but it is built on the platform's existing LangGraph dependency and follows the platform's event-sourcing conventions.

### 5.6 Audit and provenance: every output is a record

The platform treats provenance as a first-class output, not a log line. Every agent output carries a provenance record:

```json
{
  "agent": "tsc-therapeutics-strategist",
  "model_id": "claude-opus-4-8",
  "prompt_template_version": "ths-2.3.0",
  "rag_sources": [
    {"uri": "pmc://EXIST-3", "score": 0.88},
    {"uri": "ctgov://NCT0XXXXXXX", "score": 0.81}
  ],
  "input_hash": "sha256:…",
  "latency_ms": 41280,
  "event_id": "evt_…"
}
```

The record names the model and version, the prompt-template version, every retrieved RAG source with its URI, a hash of the input, and the latency. It is append-only and queryable. This is the property that earns trust with a skeptical informatics audience: an output is not "the AI said so," it is a reconstructible artifact tied to a specific model, a specific prompt, specific evidence, and a specific input. The synthetic provenance demo (Act One, §[demo]) shows this live — the molecular-genetics draft for Patient A opens its full audit trail, model and prompt version and retrieved sources visible, alongside the draft itself.

Provenance is also where the human gate is encoded. Every agent output is decision-support with a human in the loop: the Variant Curator emits an AI-labeled draft molecular-genetics report for a board-certified molecular geneticist to sign off, never an autonomous classification; the TAND Surveillance Agent surfaces patterns as pre-visit briefing material, never a diagnosis and never an interruptive alert. The provenance record is what the human reviews against. Nothing on the platform is autonomous, and the record is how that claim is kept honest.

### 5.7 Data substrate: PostgreSQL, Redis, MinIO — and the platform at institutional scale

The platform's data services are deliberately ordinary, because ordinary is auditable and portable. **PostgreSQL 16** holds the event log and structured state. **Redis 7** holds ephemeral working state. **MinIO** provides S3-compatible object storage for artifacts — generated VCFs, draft reports, surface snapshots. The TSC engine adds no exotic infrastructure; it uses the same three services every other engine uses.

At institutional scale the object tier changes. In a hospital deployment the natural substrate for the biospecimen-linked data the engine is designed to consume — the banked tuber, AML, and SEGA tissue in a biobank; the imaging archives; the longitudinal Clarity extracts — is enterprise storage built for that volume and that access pattern, which is where **the platform** storage sits in the architecture. This is described as the envelope, not as something the demo exercises: the synthetic-data demo runs on local NVMe with MinIO, and the institutional storage tier is Phase-1 work. The same line applies to the data plumbing itself — Epic Clarity/Caboodle and biobank LIMS integration are described architecturally as how the Phenome Mapper, TAND, and Trajectory agents would draw real longitudinal data, and are explicitly not built in the synthetic demo. The engine is designed to drop into that envelope; it does not pretend to already live there.

### 5.8 Hardware: the DGX Spark, with RunPod for burst

The platform's hardware target — and the machine this engine is built and demonstrated on first — is a single **NVIDIA DGX Spark**: a GB10 Grace Blackwell system delivering roughly 1,000 TOPS, with **128 GB of unified LPDDR5x memory**, a **4 TB NVMe** drive, running DGX OS. The unified-memory architecture is what makes a one-box build tractable: the Grace CPU and Blackwell GPU share a single address space, so the genomic, statistical, and inference workloads do not pay a host-to-device copy tax to move data between stages.

The Spark is sufficient for the steady-state engine — orchestration, retrieval, the tiered Claude calls (which are API-bound, not local-GPU-bound), the classical-statistics trajectory work, and Llama 3.1 70B local fallback inference. It is not sufficient for the bursty, embarrassingly parallel jobs in the build, and the platform does not pretend otherwise: GPU-accelerated Parabricks variant calling, parallel synthetic-cohort generation, and heavier local-LLM inference burst to **additional GPUs via RunPod** as needed. The synthetic-cohort regeneration — a deterministic, version-controlled ~12-hour run over 50 patients — is the canonical example of work that fans out to RunPod GPUs and then collapses back to the Spark for the steady-state demo. "Runs on a Spark with RunPod for burst" is a precise claim, not a marketing one: the steady state is genuinely one box, and the spikes are genuinely offloaded.

### 5.9 What the substrate gives the engine

Taken together, the platform hands the TSC engine almost everything it needs before a line of TSC-specific code is written: LangGraph for deterministic orchestration, Milvus with biomedical embeddings for grounded retrieval, three Claude tiers plus a local fallback for tier-routed reasoning, PostgreSQL/Redis/MinIO for event-sourced state and artifacts, the provenance and human-gate conventions, and a one-box hardware target with a burst path. The engine is scaffolded from the existing precision_oncology_agent template — Streamlit UI, FastAPI, Milvus, Pydantic, Claude — so even the application skeleton is inherited rather than invented.

What is net-new for TSC is small and specific: the multi-agent orchestrator with its shared event-sourced state, the synthetic-cohort pipeline, and the five TSC agents. Everything else is platform. That economy is the substance of the claim that TSC is the seventh engine of a real platform and not a bespoke one-off — and it is the claim the next section takes up directly, arguing why those five coordinated agents and one deterministic orchestrator add up to an engine.


## 6. The TSC Intelligence Engine as Engine 7

The preceding sections established the disease, the documentation problem, the clinical-NLP methodology we extend, and the institutional setting. This section makes the architectural claim that organizes everything that follows: what we are building is an *engine*, and it takes its place as Engine 7 of the HCLS AI Factory. That is a specific structural claim, not a marketing label, and it carries obligations. An engine is more than a clever model behind a web form: it is a coordinated system of agents with a deterministic control plane, it spans more than one data modality, and it composes capabilities that already exist in the platform rather than reimplementing them. The TSC Intelligence Engine meets that bar. The rest of this section says why, draws the boundary of what it does and does not contain, shows how it stands on the three foundational engines, and lays out the replication thesis that makes the effort worth more than a single disease.

### 6.1 Where it sits: "Seven Engines. Eight Intelligence Agents. One Platform."

The HCLS AI Factory is organized as a small number of *engines* — durable, cross-modal, coordinating systems — and a larger number of *intelligence agents* — focused single-purpose analyzers that run inside or alongside the engines. The current platform headline is **Seven Engines. Eight Intelligence Agents. One Platform.** The six pre-existing engines are (1) Genomic Foundation, (2) Precision Intelligence, (3) Therapeutic Discovery, (4) Clinical Imaging, (5) Precision Oncology, and (6) Cardiology Intelligence. The eight agents are CAR-T, Precision Biomarker, Pharmacogenomics, Precision Autoimmune, Neurology, Clinical Trial, Rare Disease Diagnostic, and Single-Cell.

The TSC Intelligence Engine is **Engine 7**. It is the first engine in the platform organized around a single disease rather than a modality (genomics, imaging) or a therapeutic posture (oncology, cardiology). That choice is deliberate, and it is the whole thesis of the TSC cluster: Tuberous Sclerosis Complex is a multi-system disorder whose management *is* a coordination problem across modalities and over time. An engine, not an agent, is the right shape for it.

### 6.2 Why it is an engine, not an agent

Three properties distinguish an engine from an agent in this platform. The TSC work has all three.

**It is five coordinated agents, not one.** An agent does one thing. The TSC Intelligence Engine is a federation of five purpose-built agents, each owning a distinct analytic surface, each with its own model tier and its own evidence contract:

| Agent | Owns | Primary modality | Model posture |
| --- | --- | --- | --- |
| 1. TSC-Variant Curator | BAM → VCF → annotation → ACMG-AMP synthesis; mosaic recovery | Genomic | Deterministic calling/annotation; Sonnet aggregation; Opus classification |
| 2. TSC-Phenome Mapper | Time-anchored HPO phenotype extraction from structured data + notes | Clinical text / structured | Sonnet per-note; Haiku normalization; Opus conflict resolution |
| 3. TSC-Trajectory Modeler | SEGA/AML growth, seizure burden, renal forecasting with prediction intervals | Longitudinal quantitative | Classical statistics (not LLM); Haiku/Sonnet for prose only |
| 4. TAND Surveillance Agent | Surfacing under-recognized neuropsychiatric signals via discourse markers | Clinical text | Sonnet per-note; deterministic scoring; Opus briefing |
| 5. TSC-Therapeutics Strategist | Six-section options brief integrating all prior agents + literature + trials | Cross-modal synthesis | Opus-class (non-negotiable) |

No single agent here is doing the others' job, and no single agent could be reasonably stretched to cover the span. The Variant Curator reasons over reads and strand bias; the Phenome Mapper reasons over note discourse and ICD-10 mappings; the Trajectory Modeler runs mixed-effects models and Gaussian process regression and deliberately does *not* let a language model touch the forecast. These are different disciplines with different correctness criteria. Bundling them under one prompt would be both a worse system and a dishonest one.

**It is cross-modal.** The engine ingests and reasons over genomic reads, structured clinical data shaped like Epic Clarity/Caboodle, free-text clinical notes, longitudinal quantitative measurements (lesion dimensions, seizure counts, renal function), and the published literature via retrieval. A modality-bound analyzer cannot answer the question the engine exists to answer — *what is the state of this TSC patient, and what should the upcoming visit attend to?* — because that question is inherently cross-modal. The mosaic variant recovered from tuber tissue (Agent 1) only becomes clinically meaningful next to the HPO phenotype timeline (Agent 2), the SEGA growth forecast (Agent 3), the under-recognized TAND signal (Agent 4), and the therapeutic option space (Agent 5). Modality crossing is not a feature we added; it is the reason the system is shaped the way it is.

**It coordinates.** The five agents do not fire independently into a void. A deterministic control plane — the **TSC-Orchestrator** — sequences them, manages their shared state, and assembles their outputs into clinician surfaces on demand. The orchestrator is a LangGraph event router, **not an LLM**. It handles thirteen event types over an event-sourced state model: an append-only PostgreSQL event log with materialized current-state projections, Redis for ephemeral TTL-scoped state, and a YAML demo configuration. Its coordination patterns are explicit and conservative: dependency-ordered enrollment (the Phenome Mapper runs first because other agents build on its HPO foundation), incremental-update minimization (re-run only what new evidence touches), demand-driven surface assembly (build a surface when a clinician asks for it, not speculatively), and conservative failure handling (a failed agent yields a surface that shows "pending" or a staleness marker — never a silently missing output). The orchestrator's full design is the subject of §7; what matters here is that *coordination is a first-class, deterministic component*, and that is the third property that makes this an engine.

A single LLM agent has none of these three properties at once. The TSC Intelligence Engine has all three. That is the structural justification for the Engine 7 designation, and it is also a credibility position: we are not relabeling a chatbot.

### 6.3 Scope and boundaries: the engine draws from CCHMC, it does not contain CCHMC

The most important boundary in this paper is the one between the engine and its institutional context. Cincinnati Children's Hospital Medical Center (CCHMC) is the setting, the beachhead, and the source of the engine's substrate — but the engine does not *contain* CCHMC's areas. It draws from them.

The physical envelope is the Winslow Research Pavilion (~45,000 sq ft, opened 2024 in Avondale), which houses the Discover Together Biobank, the institution's first centralized biospecimen freezer archive, the Gamble Vaccine Research Center, and roughly 70 research staff. The Pavilion is the envelope; it is not the engine. Within and around it, **five CCHMC areas are sources that feed the engine, not parts of it**:

1. **Discover Together Biobank** — banked tuber, AML, and SEGA tissue is the molecular substrate the TSC-Variant Curator needs to recover low-VAF somatic mosaic variants that blood testing misses. The biobank's three stated missions (find biomarkers, identify disease-causing DNA changes, understand resilience) map cleanly onto the engine's analytic surfaces.
2. **Biomedical Informatics / Hagedorn** — supplies the Marshall-Hagedorn diagnostic-uncertainty discourse methodology that the TAND Surveillance Agent extends, plus output-surfacing discipline and institutional sponsorship.
3. **TSC clinical and research program** — the concentration of TSC patients that gives all five agents something real to reason about in a Phase-1 deployment.
4. **Epic Clarity/Caboodle + biobank LIMS** — the data plumbing that would feed the Phenome Mapper, TAND, and Trajectory agents in production.

The data flow is one-directional and clean: **Source areas → TSC Engine (agents + orchestrator) → clinician surfaces.** The engine is the intelligence layer that sits on top of the sources and turns them into something a clinician can act on at the point of care. Our recurring mantras name this exactly: *a biobank without an intelligence layer is a freezer full of tubes*, and the engine is *the kitchen, not the meal*.

This boundary carries a hard honesty obligation that runs through the entire paper. **The synthetic-data demonstration is what runs now** on Adam's NVIDIA DGX Spark with extra RunPod GPUs as needed. The Epic Clarity/Caboodle integration, the biobank LIMS connections, and the imaging-AI pipelines are described *architecturally* in this paper, but they are **explicitly not built** — they are institutional Phase-1 work that requires real-data governance, IRB approval, and CCHMC engineering that the synthetic demo neither has nor claims. Where source area 4 (data plumbing) appears in this paper, it is a Phase-1 design, not a running system. We mark this wherever it is relevant because the credibility of the whole effort with a skeptical informatics audience depends on never blurring the line between "we built this" and "this is how it would connect."

### 6.4 How the engine composes the foundational engines

Engine 7 is not built on bare metal. It composes capabilities that already exist as foundational engines and reused HCLS AI Factory v1.3.0 substrate (LangGraph orchestration, Milvus RAG with BAAI/bge-large-en-v1.5 plus BiomedBERT-derived clinical embeddings over a TSC corpus partition, tiered Claude models with a local Llama 3.1 70B Instruct fallback via Ollama, PostgreSQL + Redis + MinIO, and the audit/provenance layer). The composition over the three foundational engines is direct:

- **Genomic Foundation → TSC-Variant Curator.** The platform's genomics engine provides the deterministic calling and annotation backbone the Variant Curator stands on: Parabricks/BWA-MEM alignment, GATK HaplotypeCaller and Mutect2 calling tuned mosaic-aware to recover variants at VAF ≥ 5%, and snpEff/VEP annotation. The TSC-specific work — evidence aggregation, ACMG-AMP classification synthesis against combinatorial rules, the mosaic flag, the ddPCR validation recommendation, the AI-labeled draft molecular-genetics report — is the net-new layer. The calling pipeline itself is composed, not rebuilt.
- **Precision Intelligence → TSC-Phenome Mapper and the RAG-dependent agents.** The platform's Precision Intelligence engine (formerly named the "Precision Intelligence Network" — a name we no longer use) provides the retrieval and phenotyping machinery: the Milvus vector store, the clinical embedding stack, the HPO/SNOMED-CT normalization tooling, and the literature corpus that the TAND and Therapeutics agents query. The Phenome Mapper, the TAND Surveillance Agent's PubMed/PMC grounding, and the Therapeutics Strategist's ClinicalTrials.gov and FDA-action retrieval all sit on this engine.
- **Therapeutic Discovery → TSC-Therapeutics Strategist.** The platform's therapeutic-discovery engine supplies the structured-options framing and the mechanism-aware reasoning the Strategist uses to assemble its six-section brief (Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions). The Strategist does not recommend; it composes a decision-support brief with every claim source-attributed, framed against the TSC therapeutic landscape (mTOR inhibitors such as everolimus and sirolimus, EXIST-3-era targeted epilepsy therapy, next-generation selective mTORC1 inhibitors in development, plus ketogenic diet, surgery, and VNS).

What is genuinely net-new for TSC, and therefore what this paper is actually about, is a short list: the **multi-agent orchestrator and shared event-sourced state**, the **synthetic cohort pipeline**, and the **five TSC agents** themselves. Everything underneath is composed. That is what it means to be an engine in a platform rather than a standalone application: you inherit the foundations and you build the disease-specific coordination on top. The Trajectory Modeler is the one place where the engine reaches *outside* the LLM-centric platform substrate entirely, into classical statistics (mixed-effects models, Gaussian process regression, survival analysis, Bayesian updating), precisely because forecasting SEGA and AML growth and seizure burden is a quantitative problem that a language model should not be allowed to hallucinate its way through. Composing the right tool — including a non-LLM tool — is itself an engineering position.

### 6.5 The replication thesis: swap the box labels, keep the wiring

The final reason this is an engine, and not a one-off, is that the architecture is disease-shaped but not disease-locked. TSC was chosen as the first target for reasons developed elsewhere in this paper: a tractable genetics (TSC1/TSC2 → mTOR), a well-characterized somatic-mosaicism problem with a clear no-mutation-identified cohort, an established surveillance literature (ITSC 2021 consensus), a documented neuropsychiatric under-recognition problem (TAND), and a methodology-owning sponsor at CCHMC. But the *shape* of the engine — variant curation with mosaic awareness, longitudinal phenotype mapping, quantitative trajectory forecasting, discourse-driven surveillance of an under-recognized comorbidity, and an integrating therapeutics brief, all coordinated by a deterministic event router — is a general shape for any multi-system genetic disorder managed longitudinally across modalities.

The replication mantra is **"swap the box labels, keep the wiring."** The wiring is the orchestrator, the event-sourced state model, the surface-assembly patterns, the provenance discipline, and the agent contracts. The box labels are the disease-specific knowledge each agent carries: which genes, which annotation rules, which HPO terms, which growth model, which discourse taxonomy, which therapeutic landscape. To retarget the engine, you replace the labels and keep the wiring.

The natural replication candidates are other dominantly inherited, multi-system, longitudinally-managed genetic syndromes with under-recognized neurobehavioral components and active surveillance regimes — neurofibromatosis (NF1/NF2), Rett syndrome, Williams syndrome, Cornelia de Lange syndrome, and the other mTORopathies, where the Variant Curator's mosaic-aware calling and the Therapeutics Strategist's mTOR-pathway reasoning transfer with the least relabeling. The replication chapter (§22) develops these second-wave targets in full.

The same logic applies institutionally. The CCHMC source map — biobank tissue, a biomedical-informatics methodology owner, a disease clinical program, and Epic/LIMS plumbing under a research-pavilion envelope — is itself a reusable template. Replicating to another center (for example, TGen or City of Hope) is the same "swap the box labels, keep the wiring" move at the institutional layer: different freezer, different sponsor, different EHR instance, same engine sitting on top turning banked tubes into point-of-care intelligence.

This is why the investment is worth more than one disease. Building Engine 7 well means building a replicable pattern for an entire class of rare multi-system disorders, each of which has the same fundamental problem — a fragmented, longitudinal, cross-modal management burden that no single clinician can hold in working memory across the cadence of a surveillance schedule.

With the engine defined — its agents, its boundary against the institution, its composition of the foundational engines, and its replication thesis — the next section turns to the control plane that makes the coordination real: the architecture of the TSC-Orchestrator and the shared event-sourced state on which the five agents and three clinician surfaces depend.


---

# Part III — Engine Architecture


## 7. Engine Architecture Overview

The preceding sections established why Tuberous Sclerosis Complex (TSC) is the right disease to build for and why Cincinnati Children's (CCHMC) is the right place to build it; this section assembles the answer to *how* into a single view. The TSC Intelligence Engine is Engine 7 of the HCLS AI Factory, and it earns the "engine" designation rather than the "agent" designation for three structural reasons: it is **five coordinated agents** rather than one, it is **cross-modal** (genomic, phenotypic, longitudinal-quantitative, free-text, and literature inputs reasoned over together), and it is **coordinating** — a deterministic orchestrator sequences the agents and assembles their outputs into clinician-facing surfaces. An agent answers one question. An engine runs a clinic-shaped workflow.

What follows is the whole machine in one frame: the components, the data flow that connects CCHMC's institutional sources to the clinician's screen, the cross-modal design, the four-tier model hierarchy as it lands on TSC specifically, and the single principle that governs every build decision — put reasoning where reasoning adds value and predictability where predictability is the requirement. Sections 8 through 13 then open each component up. This section is the map; those are the territory.

### 7.1 The three layers

The engine has exactly three layers, and the discipline of the architecture comes from refusing to let them blur.

1. **Sources** — CCHMC institutional data and tissue that *feed* the engine. These are not parts of the engine. They are upstream.
2. **The engine** — five TSC agents plus one deterministic orchestrator operating over shared, event-sourced state.
3. **Surfaces** — three clinician-facing applications that *consume* engine output: a pre-visit briefing, an in-visit dashboard, and an asynchronous alert surface.

Keeping sources outside the engine boundary is not a diagram convenience. It is the claim that makes replication honest. The five CCHMC source areas — described in §4 and revisited below — are the things a peer institution such as TGen or City of Hope would swap out; the engine and surfaces are the things they would keep. "Swap the box labels, keep the wiring." A biobank without an intelligence layer is a freezer full of tubes; the engine is the kitchen, not the meal.

### 7.2 Architecture diagram

```
                         CCHMC SOURCE AREAS  (feed the engine — not part of it)
                         EXPLICITLY NOT BUILT in the synthetic demo (Phase-1 institutional work)
  ┌──────────────────────────────────────────────────────────────────────────────────┐
  │  (1) Discover Together     (2) Biomedical          (3) TSC clinical &   (4) Epic    │
  │      Biobank                   Informatics /            research            Clarity/ │
  │      banked tuber/AML/         Hagedorn (TAND           program             Caboodle │
  │      SEGA tissue               methodology +            (patient            + biobank│
  │                                surfacing +              concentration)      LIMS     │
  │                                sponsorship)                                          │
  └───────┬────────────────────────────┬────────────────────┬──────────────────┬───────┘
          │ tissue/BAM                  │ method             │ cohort           │ data plumbing
          ▼                             ▼                    ▼                  ▼
  ════════════════════════════════ THE TSC INTELLIGENCE ENGINE ════════════════════════════
  │                                                                                        │
  │   TSC-ORCHESTRATOR  (deterministic LangGraph router — NOT an LLM)                       │
  │   13 event types · event-sourced state (PostgreSQL append-only log + projections)      │
  │   Redis ephemeral (TTL) · YAML demo config · conservative failure handling             │
  │   ─────────────────────────────────────────────────────────────────────────────────  │
  │                                                                                        │
  │   ┌──────────────┐  enrolled first   ┌──────────────┐   ┌──────────────┐               │
  │   │ Agent 2      │ ────────────────▶ │ Agent 1      │   │ Agent 3      │               │
  │   │ Phenome      │  (foundation)     │ Variant      │   │ Trajectory   │               │
  │   │ Mapper       │                   │ Curator      │   │ Modeler      │               │
  │   │ HPO + spans  │                   │ BAM→VCF→ACMG │   │ classical    │               │
  │   └──────┬───────┘                   └──────┬───────┘   │ statistics   │               │
  │          │                                  │           └──────┬───────┘               │
  │          │            ┌──────────────┐      │                  │                       │
  │          └──────────▶ │ Agent 4      │ ◀────┘                  │                       │
  │                       │ TAND         │                         │                       │
  │                       │ Surveillance │                         │                       │
  │                       │ Marshall-    │                         │                       │
  │                       │ Hagedorn     │                         │                       │
  │                       └──────┬───────┘                         │                       │
  │                              │       ┌──────────────┐          │                       │
  │                              └─────▶ │ Agent 5      │ ◀────────┘                       │
  │            all four prior +          │ Therapeutics │                                   │
  │            meds/AEs + RAG +          │ Strategist   │  6-section options brief          │
  │            trials + FDA       ─────▶ │ (Opus-class) │                                   │
  │                                      └──────┬───────┘                                   │
  ══════════════════════════════════════════════│══════════════════════════════════════════
                                                 │ demand-driven surface assembly
                       ┌─────────────────────────┼─────────────────────────┐
                       ▼                         ▼                         ▼
              ┌─────────────────┐      ┌─────────────────┐      ┌─────────────────┐
              │ Pre-visit       │      │ In-visit        │      │ Async alert     │
              │ briefing        │      │ dashboard       │      │ surface         │
              │ (1-screen,      │      │ (4-quadrant,    │      │ (4 categories,  │
              │  mobile)        │      │  prog. disclose)│      │  strict disc.)  │
              └─────────────────┘      └─────────────────┘      └─────────────────┘
              Standalone web apps · persistently watermarked synthetic · not Epic
```

The arrows inside the engine are not data buses; they are dependency relationships that the orchestrator enforces. The Phenome Mapper runs first because every other agent reasons over the time-anchored phenotype profile it produces. The Therapeutics Strategist runs last because it integrates all four upstream agents plus medication, adverse-event, literature, trial, and regulatory inputs. The orchestrator is the only thing that decides what runs when.

### 7.3 The five agents in one paragraph each

These are the deep-dive subjects of §8 through §12; here is just enough to read the diagram.

- **Agent 1 — TSC-Variant Curator** takes a patient BAM through deterministic variant calling (Parabricks / BWA-MEM / GATK HaplotypeCaller + Mutect2, mosaic-aware down to VAF ≥ 5%), deterministic snpEff/VEP annotation, Sonnet-tier evidence aggregation, and an Opus-tier ACMG-AMP classification synthesis validated against combinatorial rules. Its differentiating job is recovering low-VAF somatic mosaic variants in tissue that standard blood testing misses — the ~10–15% no-mutation-identified (NMI) cohort. Output is a draft molecular-genetics report for a board-certified molecular geneticist to sign; never autonomous.
- **Agent 2 — TSC-Phenome Mapper** extracts time-anchored, HPO-coded phenotypes from synthetic Epic Clarity-shaped structured data and clinical notes (Sonnet per-note; Haiku for ICD-10 / lab → HPO normalization; Opus only for rare conflict resolution). It is the foundation the other agents build on, and it emits a longitudinal HPO profile with evidence spans, a discordance log, and an ITSC surveillance-gap report.
- **Agent 3 — TSC-Trajectory Modeler** is classical statistics, not an LLM: mixed-effects models, Gaussian process regression, survival analysis, and Bayesian updating forecast SEGA/AML growth, seizure burden, and renal function at 6/12/18 months with 50% and 90% prediction intervals. It raises threshold-crossing alerts and surveillance-cadence recommendations. An LLM only writes the prose summary (Haiku) or interprets an unusual trajectory (Sonnet).
- **Agent 4 — TAND Surveillance Agent** surfaces under-recognized TSC-Associated Neuropsychiatric Disorders from longitudinal notes by applying the Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy (hedging, deferral, third-party attribution, conditional, follow-up-without-formalization) across the six TAND clusters. Sonnet does per-note discourse analysis; a deterministic layer scores and aggregates; Opus writes the briefing. It surfaces patterns as pre-visit briefing material, never as interruptive alerts, and never diagnoses.
- **Agent 5 — TSC-Therapeutics Strategist** is Opus-class, non-negotiably. It integrates all four prior agents plus medications/adherence, adverse events, PubMed/PMC RAG, a ClinicalTrials.gov snapshot, and FDA actions into a six-section structured options brief — Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions — with every claim source-attributed and the whole thing framed as decision-support, not recommendation.

### 7.4 The orchestrator is deterministic on purpose

The TSC-Orchestrator is a LangGraph event router, and it is emphatically **not** an LLM. It recognizes 13 event types, holds an event-sourced state (a PostgreSQL append-only event log plus materialized current-state projections), uses Redis for ephemeral TTL state, and reads a YAML demo configuration. Four behavior patterns define it:

1. **Dependency-ordered enrollment** — Phenome Mapper first, then the agents that consume its output.
2. **Incremental-update minimization** — a new note triggers only the agents whose inputs changed, not a full cohort recompute.
3. **Demand-driven surface assembly** — surfaces are built when requested from current-state projections, not pushed eagerly.
4. **Conservative failure handling** — a failed agent makes its surface region show "pending" or a staleness marker; it never produces a silent missing output.

A clinician must be able to trust that the dashboard reflects the data, that two identical inputs produce two identical run orders, and that nothing was quietly dropped. Routing is a correctness property, not a reasoning problem, so the router does not reason. This is the architecture's first concrete expression of the right-tool principle.

### 7.5 Cross-modal by design

TSC is a cross-modal disease and the engine is built to match. A single TSC patient generates genomic data (the tuber or AML tissue BAM), structured EHR data (encounters, problems, labs, medications), free-text clinical notes across neurology, nephrology, genetics, dermatology, and developmental specialties, longitudinal quantitative measurements (serial SEGA dimensions, seizure counts, renal function), and a moving body of external evidence (literature, trials, FDA actions). No single modality tells the story. The mosaic variant in Patient A's tissue (8.3% VAF, TSC2 frameshift) is meaningless without the phenotype that motivated tissue sequencing; the SEGA growth curve for Patient B is meaningless without the variant and the surveillance cadence; the under-recognized TAND signals scattered through Patient B's notes are invisible to any agent that only reads structured data.

The cross-modal integration happens through the shared event-sourced state, not through agents calling each other directly. Each agent writes typed outputs to the event log; downstream agents read the materialized projections. This is why the Therapeutics Strategist can integrate four upstream agents without being coupled to any of them, and why a peer institution can replace a source without rewiring the engine.

### 7.6 The four-tier hierarchy applied to TSC

The HCLS AI Factory v1.3.0 substrate gives the engine a four-tier compute-and-reasoning hierarchy. The TSC Engine reuses it deliberately, and the assignment of work to tiers is the second concrete expression of the right-tool principle.

| Tier | What it is | TSC usage |
|------|-----------|-----------|
| **Deterministic / classical** | Pipelines and statistical models with reproducible output | Variant calling (Parabricks/GATK), snpEff/VEP annotation, the entire Trajectory Modeler (mixed-effects, GP regression, survival, Bayesian), the orchestrator, ACMG combinatorial-rule validation, all scoring/aggregation layers |
| **Haiku** | Fast, cheap LLM for normalization and short prose | ICD-10 / lab → HPO normalization (Phenome Mapper); prose summaries of Trajectory forecasts |
| **Sonnet** | Mid-tier reasoning for per-document analysis | Per-note phenotype extraction (Phenome Mapper); per-note discourse analysis (TAND); evidence aggregation (Variant Curator); unusual-trajectory interpretation |
| **Opus** | Highest-stakes synthesis, validated against rules | ACMG-AMP classification synthesis (Variant Curator); rare conflict resolution (Phenome Mapper); TAND briefing summary; the entire Therapeutics Strategist (non-negotiable) |

Two reads of this table matter. First, the Opus tier is reserved for synthesis where a clinician's decision is downstream and where the model's judgment is genuinely additive — and even there it is gated (ACMG against combinatorial rules, Therapeutics against source attribution and human review). Second, the most numerically consequential agent, the Trajectory Modeler, uses no LLM for its forecasts at all. A SEGA growth projection with a 90% prediction interval is a statistical claim; it must be reproducible, calibratable, and defensible to a skeptic, and a language model is the wrong instrument for it. Haiku writes the sentence that describes the curve; it does not compute the curve.

### 7.7 The right tool for the job

Stated plainly, the governing principle is: **use an LLM where probabilistic reasoning over messy language or evidence adds value; use deterministic or classical methods where predictability, reproducibility, or numerical rigor is the requirement.** The architecture applies it at every level:

- **Routing** is deterministic (correctness, not reasoning).
- **Variant calling and annotation** are deterministic (reproducibility, and these are mature, validated tools).
- **Trajectory forecasting** is classical statistics (calibrated uncertainty, defensible intervals).
- **Phenotype and discourse extraction from notes** are LLM tasks (the input is unstructured language; the value is in reading it well).
- **ACMG synthesis and therapeutic strategy** are LLM tasks *with hard gates* (the value is in integrating heterogeneous evidence; the gate is human review and source attribution).

This is also the credibility argument. The audience is a skeptical clinical-informatics readership that has watched LLMs hallucinate citations and miscount. Putting an LLM where it can fabricate a number would forfeit trust on first contact. Putting it where it reads a note for hedging language — precisely the task the Marshall-Hagedorn methodology was published to address — is using it for what it is good at. The architecture is, in this sense, the argument: every tier assignment is a statement about what we are and are not willing to let a probabilistic model decide.

### 7.8 What is built, and what is not

One boundary must be stated here because it governs how the rest of the paper should be read. Everything inside the engine box and the surfaces — the five agents, the orchestrator, the shared state, the three clinician applications — runs **now**, on the synthetic 50-patient cohort, on Adam's NVIDIA DGX Spark (GB10 Grace Blackwell, 128 GB unified memory) with extra GPUs from RunPod for Parabricks variant calling and parallel cohort generation. The CCHMC source areas — biobank tissue ingest, Epic Clarity/Caboodle and LIMS data plumbing, and the imaging-AI pipelines — are described architecturally but are **explicitly not built**. They are institutional Phase-1 work, gated on IRB and on real-data agreements that the synthetic demo does not require. Wherever a later section describes a source-area integration, it is describing the envelope, not the running code. The synthetic cohort is persistently watermarked, the demo is not FDA-cleared, and the SaMD posture is undetermined institutional work.

With the whole machine in view — three layers, five agents, one deterministic orchestrator, three surfaces, a cross-modal data flow from CCHMC sources to the clinician's screen, and a four-tier hierarchy governed by the right-tool principle — the remaining sections can go deep on each component. Section 8 begins with the foundation the rest stands on: the TSC-Variant Curator and the mosaic-recovery problem that opens the demo.


## 8. Agent 1 — TSC-Variant Curator

The molecular diagnosis of Tuberous Sclerosis Complex is, for most patients, a solved problem. A blood draw, a targeted gene panel covering *TSC1* and *TSC2*, and standard germline variant interpretation will identify a pathogenic variant in roughly 85–90% of clinically diagnosed individuals. The TSC-Variant Curator is not built for those patients. It is built for the remainder: the 10–15% who carry a confident clinical diagnosis of TSC under the 2021 International Tuberous Sclerosis Complex (ITSC) consensus criteria and yet return a negative blood test. This is the **no-mutation-identified (NMI)** cohort, and the literature is now unambiguous about why it exists. In the large majority of these cases the causal variant is present, but it is a low-level somatic mosaic event — confined to a fraction of cells, often the affected tissue itself, at a variant allele fraction (VAF) well below what germline-tuned pipelines are designed to detect (Tyburczy et al. 2015; Giannikou et al. 2016; Lim et al. 2017). A germline caller expecting heterozygous variants near 50% VAF, with hard quality filters tuned to suppress sequencing artifact, will treat a true 6–8% mosaic variant as noise and discard it.

The clinical cost of that discard is real. A molecular diagnosis in TSC is not a formality. It enables predictive testing of relatives, informs reproductive counseling (including for the ~2/3 of cases that are de novo, where recurrence risk hinges on the possibility of parental gonadal mosaicism), establishes the genotype context for surveillance intensity, and — increasingly — anchors eligibility for genotype-aware trials and mTOR-directed therapy. An NMI patient is a patient whose family carries unquantified uncertainty and whose record lacks the molecular substrate the rest of this engine reasons over.

The TSC-Variant Curator addresses this by doing two things that a conventional panel pipeline does not. First, it is **mosaic-aware by construction**: it calls down to VAF ≥ 5% and treats low-VAF events as candidates to be adjudicated rather than artifacts to be filtered. Second, it operates on **tissue** — banked tuber, angiomyolipoma (AML), and subependymal giant cell astrocytoma (SEGA) specimens drawn from the institutional biobank substrate — where a mosaic variant absent from blood is enriched and recoverable. In the CCHMC architecture, the Discover Together Biobank is the source that makes this agent possible: it holds the tissue, and the Variant Curator is part of the intelligence layer that turns banked tubes into an interpretable molecular finding. A biobank without that layer is, as the program's framing puts it, a freezer full of tubes.

Two framing constraints govern everything below. The agent's output is a **draft for review by a board-certified molecular geneticist** — it is never autonomous and never a released clinical result. And in the build that runs today, it operates on the **synthetic cohort** described in §6: BAM files derived from NA12878 with TSC1/TSC2 variants inserted by BAMSurgeon at controlled mosaic (4–12%) and germline (~50%) VAFs. The pipeline architecture is identical to what would run on real banked tissue; the data is synthetic and persistently watermarked, and direct biobank/LIMS ingestion is institutional Phase-1 work, explicitly not built here.

### 8.1 The four-step pipeline

The Variant Curator is deliberately layered so that the deterministic, validatable steps carry the genomic burden and the language-model steps carry only the synthesis burden — and even then under explicit validation. The pipeline takes a BAM as input and emits a structured, fully provenanced draft molecular-genetics report.

```
  BAM (tissue or blood, NA12878-derived in the demo)
    │
    ▼
 [1] Variant calling          ── deterministic ──  Parabricks / BWA-MEM
     germline + somatic                            GATK HaplotypeCaller + Mutect2
     mosaic-aware (VAF ≥ 5%)                        tumor-only / low-VAF mode
    │
    ▼
 [2] Annotation               ── deterministic ──  snpEff / VEP
     functional consequence                        gnomAD v4, ClinVar, LOVD-TSC,
     population + clinical DB                       dbSNP cross-reference
    │
    ▼
 [3] Evidence aggregation      ── Claude Sonnet ──  per-variant evidence dossier
     structure, not invent                         (one dossier per candidate)
    │
    ▼
 [4] ACMG-AMP synthesis        ── Claude Opus  ──   classification reasoning
     validated vs. rule engine                      + combinatorial rule check
    │
    ▼
 Draft molecular-genetics report  →  molecular geneticist sign-off gate
```

#### Step 1 — Deterministic variant calling

Calling runs on GPU-accelerated **Parabricks** where a RunPod GPU is attached for throughput, falling back to the equivalent open-source stack (**BWA-MEM** alignment, **GATK HaplotypeCaller** for germline and **Mutect2** in tumor-only / low-VAF configuration for somatic and mosaic events). The germline and somatic call sets are unioned rather than treated separately, because the entire point of this agent is that a single sample may carry a near-50% germline event in one patient and a 6% mosaic event in another, and the pipeline must not pre-commit to one model of the truth.

The critical departure from a standard panel pipeline is the threshold and filter posture. We call to **VAF ≥ 5%** and we do not apply the aggressive low-VAF rejection that germline pipelines use to keep their false-positive rate down. That choice imports false positives by design; steps 2 and 3 exist to discharge them. Each candidate carries forward its raw evidence — VAF, supporting and total read depth, strand bias, base-quality and mapping-quality distributions, and proximity to known artifact-prone contexts — because the mosaic-versus-artifact decision is made downstream on exactly those features.

This step is deterministic and version-pinned. Given the same BAM and the same tool versions, it produces the same VCF. That property is non-negotiable for a clinical-adjacent artifact and is what makes the synthetic-cohort eval meaningful: ground truth is known because BAMSurgeon inserted it.

#### Step 2 — Deterministic annotation

The unioned VCF is annotated with **snpEff/VEP** for functional consequence (frameshift, nonsense, canonical splice, missense, in-frame indel, synonymous), then cross-referenced against **gnomAD v4** for population frequency, **ClinVar** and **LOVD-TSC** for prior clinical assertions, and **dbSNP** for known-variant status. This produces, per candidate, the raw material that ACMG-AMP criteria are evaluated against: the predicted effect on the *TSC1*/*TSC2* protein, whether the position falls in a documented hotspot or recurrently reported residue, and whether the variant has ever been seen in a population database (and at what frequency). This step too is deterministic and version-pinned; the database snapshots are recorded in provenance so a classification can be reproduced against the exact reference state it was made under.

#### Step 3 — Evidence aggregation (Claude Sonnet)

Step 3 is the first language-model step, and its job is narrow on purpose: **structure, do not invent**. For each surviving candidate, Sonnet assembles an evidence dossier from the deterministic outputs of steps 1 and 2. It does not call variants, it does not assign pathogenicity, and it is prompt-constrained to cite only fields present in its input. The dossier it produces is a normalized, machine-readable object:

```json
{
  "variant_id": "TSC2:c.4537C>T",
  "gene": "TSC2",
  "hgvs_p": "p.Arg1513Ter",
  "consequence": "stop_gained",
  "vaf": 0.061,
  "depth_total": 412,
  "depth_alt": 25,
  "strand_balance": 0.48,
  "mosaic_assessment": {
    "is_candidate_mosaic": true,
    "vaf_band": "5-12%",
    "artifact_flags": [],
    "strand_bias_pass": true
  },
  "population": { "gnomad_v4_af": 0.0, "in_dbsnp": false },
  "clinical_db": { "clinvar": "not_reported", "lovd_tsc": "not_reported" },
  "evidence_spans": ["depth_alt/depth_total", "snpEff:stop_gained", "gnomad_v4"]
}
```

Sonnet is the right tier here. The task is high-volume (every candidate in every case), pattern-regular, and bounded by its inputs; it rewards a fast, capable model that follows a fixed schema, and it does not require the deeper reasoning reserved for the final classification. The `mosaic_assessment` block is where the agent makes its explicit, auditable call on the artifact-versus-real question, surfacing strand balance and artifact flags so a human can interrogate the decision rather than trust it.

#### Step 4 — ACMG-AMP synthesis (Claude Opus), validated against a rule engine

The final step assigns a clinical classification under the ACMG-AMP framework (Richards et al. 2015) — Pathogenic, Likely Pathogenic, Variant of Uncertain Significance (VUS), Likely Benign, or Benign — together with the explicit list of criteria that fire and the reasoning that justifies each. This is the only **Opus**-tier step in the agent, because it is the step where genuine interpretive judgment is required: weighing the strength of a null-variant call against population and segregation evidence, recognizing when a mosaic context modifies the confidence of a phenotype-match criterion, and producing reasoning a molecular geneticist will actually read and either accept or contest.

Opus does not get the last word unchecked. Its proposed criteria set is **validated against a deterministic combinatorial rule engine** that encodes the ACMG-AMP combining rules (e.g., one Very Strong plus one Moderate yields Pathogenic; the published combinations that produce Likely Pathogenic, VUS, and the benign tiers). If Opus proposes a set of firing criteria, the rule engine independently computes what classification that set *must* yield and compares it to Opus's stated classification. A mismatch is a hard error: the case is flagged, the discrepancy is recorded, and the draft is escalated rather than released. This is the guardrail that lets a language model participate in a clinical-adjacent classification at all — the model proposes and reasons, but the arithmetic of criteria combination is enforced by code, not by the model's say-so.

### 8.2 ACMG-AMP criteria as applied to TSC null variants

The criteria the agent reasons over are the standard ACMG-AMP set, applied with the gene-specific particulars of *TSC1*/*TSC2*. The most consequential for this cohort are:

- **PVS1 (null variant in a gene where loss of function is a known mechanism).** *TSC1* and *TSC2* are classic loss-of-function genes — haploinsufficiency releases the brake on mTOR. A nonsense, canonical-splice, or frameshift variant that is predicted to trigger nonsense-mediated decay or remove functional protein meets PVS1, modulo the usual caveats about last-exon escape and rescue. PVS1 is Very Strong and, paired with a single Moderate criterion, is sufficient for Likely Pathogenic.
- **PM2 (absent or extremely low frequency in population databases).** A variant absent from gnomAD v4 satisfies PM2 (applied at Supporting-to-Moderate strength per current ClinGen guidance). De novo dominant TSC variants are expected to be private; absence from population data is the norm, not the exception.
- **PP4 (phenotype highly specific for a single-gene disorder).** TSC under the ITSC clinical criteria is a strongly specific phenotype. A patient meeting clinical diagnostic criteria, carrying a candidate variant in *TSC1*/*TSC2*, supports PP4. This criterion is where the engine's other agents matter: the Phenome Mapper's HPO profile (§9) is the structured phenotype evidence that justifies PP4, and the provenance ties the criterion to specific phenotype evidence spans rather than an unexplained assertion.
- **PS3 / PM1 / PM5 / PP3** are evaluated where applicable (functional studies, mutational hotspot, novel change at a known pathogenic residue, computational concordance) but are secondary for the null-variant cases that dominate the NMI cohort.

The mosaic context does not change which criteria are available, but it does change the burden on step 3's `mosaic_assessment`: a low-VAF call that fails strand-balance or carries artifact flags is not allowed to anchor a Pathogenic classification, regardless of how clean its consequence annotation looks. The agent is explicitly tuned to avoid false-positive Pathogenic calls, because in this setting a false positive is a wrong molecular diagnosis delivered to a family.

### 8.3 Inputs, outputs, and the human gate

**Inputs.** A BAM (tissue-derived in the intended deployment, NA12878-derived with inserted variants in the demo), the pinned reference build, and the annotation database snapshots. Optionally, the patient's HPO profile from the Phenome Mapper to support PP4.

**Outputs.** A single structured draft comprising: (1) a ClinVar-spec interpretation per reportable variant (HGVS at cDNA and protein level, classification, criteria with reasoning); (2) a **mosaic flag** carrying VAF, supporting/total reads, strand balance, and artifact assessment for any low-VAF event; (3) an explicit **ddPCR validation recommendation** for mosaic candidates, because orthogonal confirmation of a low-VAF event is standard practice before it informs care; (4) an **AI-labeled draft molecular-genetics report** in clinician-readable prose; and (5) full provenance.

**The human gate.** Nothing the Variant Curator produces is a released result. Every draft carries an explicit AI-generated label and routes to a board-certified molecular geneticist for sign-off. The agent's role is to surface the candidate a germline pipeline missed, assemble the evidence, and propose a defensible classification — compressing the analyst's time, not replacing the analyst's judgment or signature.

### 8.4 Provenance

Every Variant Curator output carries the engine-wide provenance record: model id and version for each LLM step (e.g., the Sonnet and Opus model strings), prompt-template version, the retrieved RAG/database sources with URIs and snapshot dates (gnomAD v4, ClinVar, LOVD-TSC), the input BAM hash, tool versions for the deterministic steps (BWA-MEM/Parabricks, GATK, snpEff/VEP), the rule-engine validation result, and per-step latency. The record is append-only and queryable. For a clinical-adjacent artifact this is not optional instrumentation — it is what allows a geneticist to reconstruct exactly how a classification was reached and against which reference state, and it is what makes the eval below reproducible.

### 8.5 Evaluation targets (demo, against synthetic ground truth)

These targets measure the agent against the BAMSurgeon-inserted truth of the synthetic cohort. They are demo correctness checks, **not clinical validation**:

| Target | Threshold |
|---|---|
| Mosaic variant recovery | Recover all **7** mosaic variants at VAF ≥ 5% (5 TSC2-mosaic, 2 TSC1-mosaic) |
| Classification accuracy | Correct ACMG-AMP class for every reportable variant |
| False-positive Pathogenic | **Zero** false-positive Pathogenic calls |
| Latency | **< 5 min** per case |

The seven-mosaic recovery target is the single most important number for this agent. The cohort was built with exactly seven mosaic carriers precisely so that this agent has a defined, falsifiable bar: if it misses one, that is a visible failure against known truth, not a soft judgment call.

### 8.6 Worked example — Patient A (Act One of the demo)

Patient A is a 4-year-old girl meeting ITSC clinical criteria for TSC who is **NMI on blood testing** — the canonical case this agent exists for. The input is a BAM derived from banked **tuber tissue** (synthetic: BAMSurgeon-inserted into an NA12878-derived BAM at a controlled low VAF).

1. **Calling.** The mosaic-aware caller recovers a *TSC2* frameshift variant at **8.3% VAF** — well below germline expectation, above the 5% floor, with adequate supporting reads and balanced strand support. A standard germline panel on blood, where the variant is absent or far rarer, returns the NMI result that brought the patient here.
2. **Annotation.** snpEff/VEP classifies the variant as a frameshift (predicted loss of function); gnomAD v4 shows it absent from the population; ClinVar and LOVD-TSC have no prior assertion.
3. **Aggregation (Sonnet).** The evidence dossier records VAF 8.3%, a clean `mosaic_assessment` (no artifact flags, strand balance within tolerance, VAF band 5–12%), frameshift consequence, gnomAD AF = 0, and the supporting evidence spans.
4. **Synthesis (Opus, rule-validated).** The criteria that fire are **PVS1** (frameshift null in a known loss-of-function gene), **PM2** (absent from gnomAD v4), and **PP4** (TSC clinical phenotype highly specific, supported by the HPO profile). The combinatorial rule engine confirms PVS1 + PM2 + PP4 yields **Likely Pathogenic**, matching Opus's stated classification — no discrepancy, so the draft proceeds. The output attaches a **mosaic flag** and a **ddPCR validation recommendation** for orthogonal confirmation of the 8.3% event.

The result is a draft molecular-genetics report, explicitly AI-labeled, that converts an NMI patient into a patient with a Likely Pathogenic mosaic *TSC2* finding and a clear confirmation path — pending the molecular geneticist's sign-off. In the live demo (Act One), the same case runs in front of the audience with the **full audit trail open**: the recovered call, the evidence dossier, the firing criteria, the rule-engine check, and the provenance record. The argument is not "trust the model"; it is "watch the deterministic floor catch the variant, watch the model reason within an enforced rule set, and inspect every step."

### 8.7 Known limitations

The honest boundaries of this agent are part of its credibility:

- **VAF floor.** The 5% threshold is a deliberate balance between sensitivity and the artifact burden the downstream steps can responsibly discharge. True mosaic variants below 5% — which the literature documents do occur — are out of scope for this configuration. Lowering the floor would require deeper sequencing, duplex/UMI-based error suppression, and a correspondingly heavier artifact-adjudication burden; that is acknowledged future work, not a current claim.
- **Synthetic data.** Everything here runs on BAMSurgeon-inserted variants in NA12878-derived BAMs. The pipeline wiring is real and the eval against known truth is meaningful, but the demo does not establish performance on real banked tissue, with its true noise structure, contamination, and sequencing-chemistry artifacts. Direct biobank/LIMS ingestion is Phase-1 institutional work.
- **ddPCR is recommended, not performed.** The agent recommends orthogonal validation of mosaic calls; it does not perform it. Confirmation remains a wet-lab step under the geneticist's direction.
- **Not autonomous, not cleared.** Every output is a draft for sign-off. The agent is not FDA-cleared, its SaMD posture is undetermined and is institutional work, and no output is a released clinical result.
- **Scope of evidence.** PP4 is only as strong as the phenotype evidence behind it; PS3 functional evidence is rarely available for private variants. The agent is built to under-call rather than over-call, accepting more VUS outcomes in exchange for zero false-positive Pathogenic results.

With a defensible molecular finding in hand — including, for the first time, the NMI patients the engine is designed to recover — the question becomes what the rest of the clinical record says about each patient over time. That longitudinal phenotype, structured and HPO-coded, is the foundation the next four agents build on, and it is the subject of §9: the TSC-Phenome Mapper.


## 9. Agent 2 — TSC-Phenome Mapper

The Variant Curator (§8) answers *what is broken in the genome*. The Phenome Mapper answers *what the body has actually done about it, over time, as documented by the people who watched.* These are different questions with different evidence bases, and TSC is a disease where the second one is chronically under-captured. The TAND consensus work and the TOSCA registry both make the same observation from different angles: a large fraction of what a TSC patient experiences never makes it into a structured, queryable, time-anchored form. It lives in free-text notes, scattered across years, across subspecialty clinics, in the hedged and conditional language clinicians use when something is noticed but not yet acted on. The Phenome Mapper exists to turn that documentary record into a clean longitudinal phenotype that the other agents can stand on.

This is why it enrolls first. The TSC-Orchestrator (§14) imposes a dependency-ordered enrollment in which the Phenome Mapper runs before the Trajectory Modeler, the TAND Surveillance Agent, and the Therapeutics Strategist for any given patient. The Trajectory Modeler needs time-anchored measurements to fit a curve. The TAND agent needs the longitudinal note corpus and the HPO scaffold to know which clusters are already represented and which are conspicuously absent. The Therapeutics Strategist needs the full phenotype to reason about indication and contraindication. None of them can do their work against raw notes; they consume the Phenome Mapper's structured output. It is the foundation, and a weak foundation propagates quietly into every downstream forecast and brief. That places a particular burden on this agent: its errors are the most expensive in the system, so its design leans hard on conservatism and on an auditable evidence trail.

### 9.1 What it does

The Phenome Mapper ingests two classes of input for each synthetic patient and emits one structured longitudinal artifact plus two diagnostic byproducts.

**Inputs.** First, structured data shaped like an Epic Clarity / Caboodle extract: problem-list entries with ICD-10 codes, encounter records with dates and departments, laboratory results with LOINC-style codes and values, medication orders, and procedure/imaging order records. In the synthetic cohort this is generated by Synthea with the TSC modules (§7) into a relational, Clarity-shaped form; it is not a live Epic feed, and that distinction is load-bearing — the real Clarity/Caboodle plumbing is institutional Phase-1 work, explicitly not built in this demo. Second, the clinical note corpus: the ~600–1,000 frontier-model-generated, watermarked-synthetic notes spanning neurology, nephrology, genetics, dermatology, ophthalmology, and general pediatrics encounters, each carrying a date, an author role, and an encounter type.

**Output.** The primary artifact is a **longitudinal HPO profile**: a time-ordered set of phenotype assertions, each coded to a Human Phenotype Ontology (HPO) term, each anchored to one or more dated source events, and each carrying an evidence span — the exact text or structured field that supports it — plus a provenance record. A phenotype is not a single boolean; it is a series of observations over time. "Subependymal giant cell astrocytoma" (HPO term for SEGA-class lesion) asserted at age 11 from a radiology note, re-asserted at 11.5 and 12 from follow-up MRIs, is three time-anchored observations of one phenotype, and the downstream Trajectory Modeler will want all three.

The two byproducts are a **discordance log** and an **ITSC surveillance-gap report**, described in §9.3 and §9.4.

### 9.2 The extraction pipeline and model tiering

The Phenome Mapper is the agent where the tiered-model strategy is most visible, because the work decomposes cleanly into tasks of very different difficulty and the wrong model on the wrong task is either expensive or unreliable.

**Per-note extraction — Claude Sonnet.** Each clinical note is processed independently by Sonnet under a versioned prompt template that asks for a structured list of asserted, denied, and historical phenotypes, each with the supporting text span and an assertion polarity (present / absent / resolved / family-history). Sonnet is the right tier here: per-note clinical-language extraction with negation and temporality is harder than mechanical normalization but does not, in the common case, require the heaviest reasoning. The note's own date supplies the time anchor; relative temporal phrases ("first noted two years ago") are captured as text but are *not* trusted as precise anchors — they are flagged for the discordance pass rather than silently resolved into a date.

**Normalization — Claude Haiku.** Structured-data terms (ICD-10 problem codes, LOINC-coded labs) and the raw phenotype strings emitted by the Sonnet pass are normalized to HPO terms by Haiku against a controlled crosswalk. This is high-volume, low-ambiguity mapping work — "tuberous sclerosis" → its HPO/SNOMED anchors, an elevated renal-imaging finding → the appropriate angiomyolipoma term — and Haiku is both fast enough to run it across the whole cohort and accurate enough when the crosswalk constrains the target vocabulary. Where a term has no confident HPO mapping it is parked in the discordance log rather than force-fit.

**Conflict resolution — Claude Opus.** A small minority of cases involve genuine semantic conflict: two notes from the same week asserting present and absent for what may or may not be the same finding; a structured problem-list code that contradicts the narrative; a phenotype whose HPO placement is ambiguous because the clinical description spans two terms. These go to Opus, which is reserved for exactly this rare, high-stakes reasoning and is not run on the common path. Opus produces a resolution *with a rationale*, and that rationale is written into the discordance log so a reviewer can see not just the answer but the argument.

The aggregation and time-ordering layer that assembles per-note and per-record outputs into the longitudinal profile is deterministic code, not a model. Deduplication (the same phenotype seen in three consecutive visits collapses to one phenotype with three observations), temporal sorting, and evidence-span attachment are mechanical and must be reproducible; putting them in a model would sacrifice provenance for no benefit.

```json
{
  "patient_id": "B",
  "hpo_term": "HP:0009734",
  "hpo_label": "Subependymal giant cell astrocytoma (class)",
  "first_observed": "2024-03-12",
  "polarity": "present",
  "observations": [
    {
      "date": "2024-03-12",
      "source_event": "MRI_BRAIN_RPT_00417",
      "source_type": "imaging_report",
      "evidence_span": "subependymal nodule near the foramen of Monro measuring 0.8 cm",
      "extracted_by": "sonnet:phenome-extract:v3",
      "normalized_by": "haiku:hpo-norm:v2"
    },
    {
      "date": "2024-09-30",
      "source_event": "MRI_BRAIN_RPT_00461",
      "source_type": "imaging_report",
      "evidence_span": "interval growth to 1.1 cm",
      "extracted_by": "sonnet:phenome-extract:v3",
      "normalized_by": "haiku:hpo-norm:v2"
    }
  ],
  "provenance": {
    "input_hash": "sha256:…",
    "assembled_by": "phenome-aggregator:v4",
    "review_status": "draft"
  }
}
```

Every observation carries the model id and prompt-template version that produced it, the source event it came from, and the literal evidence span — the append-only, queryable provenance the platform requires of every agent output (§5). A reviewer who distrusts a single phenotype can trace it to one sentence in one dated note.

### 9.3 The discordance log

Phenotyping from real-world documentation is full of disagreement, and the honest thing to do with disagreement is to surface it rather than launder it into false confidence. The discordance log is the Phenome Mapper's record of every case where the evidence did not cleanly agree:

- **Polarity conflicts**: present in one source, absent in another, for the same or overlapping HPO term within a time window.
- **Structured-vs-narrative conflicts**: a problem-list ICD-10 code with no narrative support, or a narrative finding never coded.
- **Temporal ambiguity**: relative time phrases that could not be anchored to a date with confidence.
- **Normalization parks**: terms with no confident HPO mapping.
- **Opus resolutions**: every conflict escalated to Opus, with its rationale.

The discordance log is not an error report to be cleared; it is a clinical-reasoning artifact. Some discordances are real history — a finding that was present, treated, and resolved — and the log preserves that nuance instead of overwriting it. It is also the agent's primary failure-transparency mechanism: when the Phenome Mapper is uncertain, the uncertainty is visible and attributable, which is the posture this audience expects and which the conservative-failure design of the orchestrator depends on.

### 9.4 The ITSC surveillance-gap report

The longitudinal profile makes a second analysis cheap, and it is one of the most directly useful outputs the engine produces. The 2021 International TSC Consensus (ITSC) surveillance guidelines specify recommended monitoring cadences by organ system and age: brain MRI intervals for SEGA surveillance, renal imaging for angiomyolipoma, EEG and seizure review, ophthalmologic exam, dermatologic exam, blood-pressure monitoring, TAND screening at key developmental points, and so on. The Phenome Mapper already knows, from the structured encounter and order records, *when each surveillance category was last performed* for a patient. A deterministic rules layer compares the observed cadence against the ITSC-recommended cadence for that patient's age and known phenotype and emits a gap report: which surveillance is current, which is approaching due, and which is overdue, with the date of the last qualifying event for each.

This is intentionally a deterministic comparison against published guideline cadences, not a model judgment — the guideline is the authority, and the value is in reliably noticing what a busy multi-subspecialty care team can lose track of across years and clinics. The gap report feeds the pre-visit briefing surface (§13) as concrete, low-controversy, high-utility content: a TSC patient overdue for renal imaging is a fact a clinician wants surfaced before the visit, framed as a watchlist item, never as an alert and never as a diagnosis.

### 9.5 Patient B walkthrough

Patient B is the 12-year-old male carrying the TSC2 nonsense variant c.3037C>T (p.Arg1013Ter), and he is the cohort's longitudinal anchor. His record spans roughly four years of encounters across neurology, nephrology, genetics, and general pediatrics. The Phenome Mapper assembles, from his notes and structured data, a longitudinal HPO profile that includes: the molecularly confirmed TSC diagnosis; a SEGA-class subependymal lesion near the foramen of Monro observed at 0.8, then 1.1, then 1.3 cm across serial MRIs; bilateral renal angiomyolipomas with a dominant lesion near 2.8 cm; well-controlled focal seizures on a stable regimen; and cutaneous findings consistent with TSC. Each is time-anchored to its source imaging reports and clinic notes, with evidence spans attached.

Two things matter about this profile downstream. First, the three dated SEGA observations are exactly the input the Trajectory Modeler (§10) needs to fit a growth curve and forecast a threshold crossing in the 12–18-month window — the Phenome Mapper does no forecasting itself, but it produces the time series that makes forecasting possible. Second, scattered through Patient B's notes are the *under-recognized TAND signals* — hedged, conditional, deferred mentions of behavioral and academic concerns that were noted but never formalized. The Phenome Mapper's job here is narrow and important: it records the longitudinal note corpus and the HPO scaffold, and it notes where TAND-relevant structured findings exist, but it deliberately does *not* attempt the discourse-marker analysis that surfaces the soft signals. That is the TAND Surveillance Agent's distinct methodology (§12, the Marshall–Hagedorn taxonomy), and conflating the two would blur a clean separation of concerns. The Phenome Mapper hands TAND a well-organized substrate; TAND brings the specialized reading.

### 9.6 Evaluation

Against the synthetic cohort's ground truth — and this is demo-grade evaluation against a known-answer synthetic record, not clinical validation — the Phenome Mapper targets **phenotype recall ≥ 90% and precision ≥ 85%**, computed over the per-patient HPO assertions, with the full 50-patient cohort processed in **under one hour** end to end. Recall is weighted toward the surveillance-relevant phenotypes that drive downstream agents (SEGA, AML, seizure burden, renal function); a missed cutaneous finding and a missed SEGA observation are not equally costly, and the eval reflects that. Precision is enforced hard because false phenotypes propagate: a spurious phenotype can trigger a spurious trajectory forecast or a spurious surveillance-gap entry. Time-anchoring is scored separately — an observation attached to the wrong date is counted as an error even if the phenotype is correct, because the Trajectory Modeler is sensitive to it.

### 9.7 Limitations

The honest limitations are several and we state them plainly. The agent extracts from synthetic notes written by a frontier model from published templates; real clinical documentation is messier, more abbreviated, more idiosyncratic per author, and the recall/precision figures here do not transfer to real text without re-evaluation against a real (IRB-governed) corpus. The structured input is Clarity-*shaped*, not a real Clarity/Caboodle extract, and the mapping from a production EHR's quirks to this clean schema is institutional Phase-1 work that this demo does not undertake. HPO normalization is only as good as the crosswalk, and rare or compound clinical descriptions can be parked rather than coded. The surveillance-gap report assumes the structured order/encounter record faithfully reflects what was actually done — true by construction in the synthetic cohort, an assumption to be validated against real data. And the entire output is a *draft* phenotype: like every agent in the engine, the Phenome Mapper augments, it does not adjudicate, and its profile is reviewable, sourced, and never the final clinical word.

With a time-anchored phenotype in hand — most importantly the serial SEGA, AML, and seizure measurements for Patient B — the engine can move from describing the present to forecasting the future. That is the work of Agent 3.


## 10. Agent 3 — TSC-Trajectory Modeler

The two preceding agents establish *what is true now*: the Variant Curator resolves the molecular substrate, and the Phenome Mapper assembles the time-anchored HPO profile. The Trajectory Modeler answers a different question, the one a TSC clinician actually carries into a surveillance visit: *where is this patient heading, and how fast?* A subependymal giant cell astrocytoma (SEGA) measured at 1.1 cm today is not the clinical concern; the concern is whether it will reach the size and location at which it threatens the foramen of Monro and precipitates hydrocephalus, and whether that crossing falls inside or outside the next surveillance interval. An angiomyolipoma (AML) at 3.6 cm is watched; one trending toward the ~4 cm bleeding-risk threshold changes the cadence of imaging and the conversation about embolization or mTOR inhibition.

This is a forecasting problem with quantified uncertainty, longitudinal repeated measures, censored time-to-event structure, and prior knowledge drawn from a literature with decades of natural-history data. It is, deliberately and explicitly, **not** a problem we hand to a large language model.

### 10.1 Why classical statistics, not an LLM

The architectural decision here is the most opinionated in the entire engine, so it deserves a direct defense. The Trajectory Modeler is built from mixed-effects regression, Gaussian process (GP) regression, parametric and semi-parametric survival models, and Bayesian sequential updating. An LLM appears only at the periphery, and never touches a numeric forecast. Four reasons.

**1. Calibrated uncertainty is the deliverable, not a garnish.** A surveillance recommendation lives or dies on the prediction interval, not the point estimate. "SEGA reaches 1.5 cm in roughly 14 months (50% PI 11–17 mo, 90% PI 8–22 mo)" is an actionable statement; a confidently phrased single number is a liability. Mixed-effects and GP models produce intervals with defined coverage semantics that we can audit against held-out trajectories. An LLM asked to "predict the growth" produces fluent numbers whose stated confidence has no relationship to any frequentist or Bayesian coverage guarantee. For an audience that includes Dr. Hagedorn's informatics group, a forecast we cannot characterize the calibration of is not a forecast.

**2. Determinism and provenance.** Run the same patient series twice and a classical model returns the identical forecast, byte for byte, with the same coefficients, the same kernel hyperparameters, and the same likelihood. That is a hard requirement for anything that influences a surveillance schedule, and it is something even a temperature-0 LLM does not reliably guarantee across model versions. Every Trajectory Modeler output carries the model family, the fitted parameter vector, the data points consumed (with input hash), and the software/library versions — the same provenance discipline applied across the engine, but here it actually reconstructs the math.

**3. Data efficiency at n = a-handful-of-scans.** A real TSC patient might have three to six SEGA measurements over several years. Classical hierarchical models are designed for exactly this regime: the population (fixed-effects) prior, learned across the cohort, *regularizes* the individual (random-effects) fit so a patient with two data points still receives a sensible, shrinkage-stabilized trajectory. An LLM has no principled mechanism to borrow strength across patients in a quantitatively correct way.

**4. The literature is a prior, and we can encode it directly.** SEGA growth, AML growth rates, and renal function decline have published natural-history ranges. We fold these in as informative priors and as cohort-level fixed effects rather than hoping a language model has internalized them correctly and will reproduce them without drift.

The honest counter-argument: classical models impose functional-form assumptions (linearity on a transformed scale, a chosen GP kernel, a proportional-hazards assumption) that can be wrong. We address that with model checking, with the GP's flexibility where linearity is implausible, and — this is where the LLM legitimately earns a seat — with a Sonnet-tier review step that flags trajectories the model fits *well numerically but implausibly clinically* (Section 10.5). The division of labor is the point: math owns the numbers; language owns the prose and the "does this even make sense" sniff test.

### 10.2 The four forecasting targets and their models

The agent forecasts four clinically governed quantities, each at **6, 12, and 18 months**, each reported with **50% and 90% prediction intervals**. Model choice follows the data-generating structure of each quantity.

| Target | Quantity | Model family | Why |
|---|---|---|---|
| SEGA growth | Max diameter (mm) over time | GP regression per lesion + mixed-effects population prior | Nonlinear, possibly accelerating; few points; needs smooth interpolation with widening tails |
| AML growth | Max diameter (mm), per kidney/lesion | Linear/log-linear mixed-effects | Growth often slow and approximately log-linear; population shrinkage stabilizes sparse series |
| Seizure burden | Monthly countable-seizure rate | Negative-binomial mixed-effects (count, overdispersed) | Counts with overdispersion; random patient intercept/slope |
| Renal function | eGFR (mL/min/1.73 m²) trajectory | Linear mixed-effects + survival model for threshold | Slope estimation plus time-to-CKD-stage crossing |

**SEGA — Gaussian process regression.** Each SEGA lesion is modeled as a latent growth function with a GP prior, using a composite kernel (a smooth long-scale component plus a small white-noise term for measurement error, with measurement error informed by inter-reader variability in the synthetic imaging reports). The population mixed-effects fit supplies the mean function and length-scale prior so a lesion with three measurements does not produce an absurdly confident or absurdly wide forecast. The GP is the right tool because SEGA growth is not reliably linear — it can accelerate — and the GP's posterior naturally widens as we extrapolate from 6 to 18 months, which is exactly the honesty we want in the 90% band.

**AML — log-linear mixed-effects.** `log(diameter) ~ time` with random intercept and slope per lesion, nested within patient. AMLs tend to grow slowly and approximately exponentially over the windows of interest; the log scale linearizes this and keeps prediction intervals non-negative on the diameter scale.

**Seizure burden — negative-binomial mixed-effects.** Monthly seizure counts (where the Phenome Mapper has extracted them with evidence spans) are overdispersed counts, so Poisson is inadequate; a negative-binomial GLMM with a random patient intercept and a time slope captures both the rate and its trend, and propagates the dispersion into the predictive interval.

**Renal function — linear mixed-effects + survival.** eGFR trajectory is modeled as a linear mixed-effects slope; a parallel time-to-event model (parametric Weibull, with a Cox cross-check) estimates time to crossing a clinically meaningful eGFR threshold, handling right-censoring correctly.

### 10.3 Threshold crossing and surveillance cadence

Forecasts are not the product; *decisions* are. The agent maps each forecast onto the clinically governed thresholds drawn from the ITSC 2021 consensus surveillance framework and standard TSC practice, and computes the **probability and expected timing of crossing**:

- **SEGA**: growth velocity and proximity to the foramen of Monro; sustained growth in a previously stable lesion; absolute size in context. The headline output is the probability the lesion crosses a pre-specified attention threshold *within each surveillance window*.
- **AML**: approach to the ~4 cm bleeding-risk threshold.
- **Renal**: time-to eGFR stage crossing.
- **Seizure**: a sustained worsening trend above a per-patient baseline.

For each, the agent emits `P(cross before 6/12/18 mo)` as a calibrated probability, plus a **surveillance-cadence recommendation**: when the model places a meaningful probability mass of a threshold crossing *before* the next guideline-default imaging interval, it recommends tightening the cadence (e.g., "consider 6-month rather than 12-month brain MRI") and shows the math behind it. Crucially, the recommendation is *relative to* the ITSC default, never a replacement for it — the agent surfaces the gap between guideline cadence and forecast risk, and a clinician decides. This is decision-support with a human gate, identical in posture to every other agent in the engine.

A representative output object:

```json
{
  "patient_id": "B",
  "target": "SEGA_right_foramen_of_Monro",
  "model": {"family": "gaussian_process", "kernel": "RBF+White",
            "length_scale": 9.4, "noise_var": 0.06,
            "lib": "scikit-learn 1.4.2"},
  "observed": [{"t_mo": 0, "mm": 8}, {"t_mo": 12, "mm": 11}, {"t_mo": 22, "mm": 13}],
  "forecast": [
    {"horizon_mo": 6,  "point_mm": 14.0, "pi50": [13.3, 14.8], "pi90": [12.4, 15.7]},
    {"horizon_mo": 12, "point_mm": 15.1, "pi50": [13.9, 16.4], "pi90": [12.2, 18.0]},
    {"horizon_mo": 18, "point_mm": 16.2, "pi50": [14.3, 18.2], "pi90": [11.9, 20.6]}
  ],
  "threshold": {"mm": 15.0, "label": "attention_threshold",
                "p_cross_by": {"6_mo": 0.07, "12_mo": 0.52, "18_mo": 0.78}},
  "cadence_recommendation": "Forecast places ~0.52 probability of crossing the 15 mm "
    "attention threshold before the next default 12-month MRI. Consider 6-month interval.",
  "provenance": {"input_hash": "…", "fit_timestamp": "…", "data_points": 3}
}
```

### 10.4 Patient B SEGA — the worked example and the eval anchor

Patient B (12-year-old male, *TSC2* c.3037C>T p.Arg1013Ter) is the synthetic cohort's trajectory exemplar and the eval target for this agent. His SEGA series is constructed at **0.8 → 1.1 → 1.3 cm** over the longitudinal window at the foramen of Monro — a clinically realistic ~2–4 mm/yr trajectory in a location where continued growth matters. He also carries bilateral AMLs (~2.8 cm) and well-controlled focal seizures, giving the agent live multi-target work.

The eval criterion is concrete and falsifiable against synthetic ground truth: **the agent must forecast Patient B's SEGA crossing its attention threshold within the 12–18-month window, and must do so without raising false alarms** on his stable AMLs or his controlled seizures. The GP fit on the three SEGA points should place the bulk of the threshold-crossing probability in that window (as in the object above), recommend tightening MRI cadence, and simultaneously report his AMLs as low crossing probability against the 4 cm threshold and his seizure trend as flat. That combination — a true positive on the lesion that is moving and clean negatives on the two domains that are not — is exactly the discrimination a surveillance clinician needs, and it is precisely what an uncalibrated narrative forecaster would blur. The Patient B forecast is what populates the **trajectory-forecast quadrant** of the in-visit dashboard (Section 14) and feeds threshold-crossing items into the pre-visit briefing's watchlist.

### 10.5 Where the LLM is allowed in — and where it is not

Two narrow, peripheral roles, both downstream of frozen numbers:

- **Haiku — prose summary only.** Once the models have produced forecasts and intervals, Haiku renders a one- or two-sentence plain-language gloss for the briefing surface ("SEGA is growing slowly but on a track that may cross the attention threshold before the next scheduled MRI"). It receives the computed numbers as structured input and is constrained, by prompt template and by a post-generation numeric-consistency check, to restate them without inventing or altering any value. If the prose and the JSON disagree, the prose is rejected.
- **Sonnet — unusual-trajectory interpretation.** When a fit triggers a heuristic flag — an outlier residual, an implausible velocity, a contradiction between two targets (e.g., worsening renal slope with no corresponding clinical signal in the HPO profile), or a poor model-fit diagnostic — Sonnet produces an interpretive note for the clinician describing *why the pattern is unusual and what it might warrant*. It never edits the numeric forecast; it annotates it.

No numeric forecast, no prediction interval, no threshold-crossing probability, and no cadence recommendation is ever produced, adjusted, or "smoothed" by an LLM. That boundary is the agent's defining design commitment and is enforced in code, not by convention: the LLM calls take the model outputs as immutable inputs.

### 10.6 Evaluation and limitations

Evaluation is against the synthetic cohort's known generative trajectories — this is **demo validation, not clinical validation**, and we say so plainly. The targets:

- **Calibration**: empirical coverage of the 50% and 90% prediction intervals across the cohort's forecastable targets should approximate nominal coverage (we report the coverage table, not just the point hits).
- **Patient B SEGA**: threshold crossing forecast in the 12–18-month window; AMLs and seizures reported as non-alarming (no false alarms), as above.
- **Determinism**: bit-identical forecasts on repeated runs of the same input series.

The limitations are real and we lead with them. First, the cohort is synthetic, so the natural-history shapes are only as good as the published ranges we generated them from; calibration on synthetic data does not transfer to a calibration claim on real patients. Second, sparse series (two or three measurements) yield wide, prior-dominated intervals — correctly, but a clinician must read those wide bands as "we don't yet know," not as precision. Third, the proportional-hazards and functional-form assumptions can fail in individuals; the Sonnet flag layer mitigates but does not eliminate this. Fourth, the engine starts at BAM and at extracted measurements — it does **not** ingest DICOM, so all imaging-derived sizes come from the synthetic radiology reports the cohort pipeline produced, and a production deployment's measurement noise would depend on real radiology workflows and image-analysis pipelines that are explicit institutional Phase-1 work, not built here.

What the agent does deliver, today, on the Spark, is a deterministic, provenance-complete, uncertainty-calibrated forecasting layer that turns the Phenome Mapper's timeline into forward-looking, threshold-aware surveillance guidance — the quantitative spine beneath the clinician surfaces. With trajectories established, the engine turns from numbers back to language: the next agent reads the longitudinal notes not for what they measure, but for what they quietly hesitate to say.


## 11. Agent 4 — TAND Surveillance Agent (and the Marshall-Hagedorn Methodology)

The Trajectory Modeler answers "what is measured, where is it heading?" This agent answers a harder question: "what is in the chart that no one has named yet?" TSC-Associated Neuropsychiatric Disorders (TAND) are the most prevalent manifestation of Tuberous Sclerosis Complex and the most reliably under-managed. This is the section where the engine's design and Cincinnati Children's own published informatics research are the same idea, so it warrants a deeper treatment than the other agents.

### 11.1 Why TAND is a documentation problem before it is a diagnostic one

TAND spans the behavioral, psychiatric, intellectual, academic, neuropsychological, and psychosocial difficulties seen in TSC, and affects roughly 90% of patients across the lifespan. As established in the disease chapter, 30 to 50% of TAND features are never formally assessed in routine care, and ITSC 2021 responded by recommending at-minimum-annual TAND screening and lifetime use of the TAND-L checklist.

The failure mode is not clinician ignorance. A neurologist managing a child's refractory epilepsy, a nephrologist tracking a 4 cm angiomyolipoma, and a developmental pediatrician each see a slice of a child whose attention, mood, sleep, and social function drift across years. The signal is real and it is in the record. It is just distributed across many notes, many authors, and many encounters, and it is rarely escalated to a formal diagnosis or referral in any single visit. Critically, the way clinicians write about something they have noticed but not yet acted on is linguistically distinctive. They hedge. They defer. They attribute the concern to a parent or teacher. They make it conditional on a future visit that may or may not happen.

That observation — that diagnostic uncertainty and unaddressed concern leave a recognizable *documentation* fingerprint, and that this fingerprint is computationally detectable — is the load-bearing premise of this agent. It is also precisely the premise of Dr. Philip Hagedorn's published clinical-NLP research at CCHMC. The TAND Surveillance Agent is an extension of his own work, not an external graft, and that is the reason it leads the CCHMC alignment case rather than following it.

### 11.2 The six TAND clusters as the analysis scaffold

The agent organizes everything around the six consensus TAND clusters. They are the unit of analysis, the unit of scoring, and the unit of presentation in the briefing.

| Cluster | Examples the agent watches for in notes |
|---|---|
| Behavioral | Aggression, self-injury, hyperactivity, sleep disturbance, mood lability |
| Psychiatric | Anxiety, depressive symptoms, autism spectrum features, ADHD |
| Intellectual | Global developmental delay, intellectual disability, regression |
| Academic | Reading/writing/math difficulty, IEP/504 mentions, school accommodation gaps |
| Neuropsychological | Attention, memory, executive-function, visuospatial deficits |
| Psychosocial | Self-esteem, family stress, peer relationships, caregiver burden |

The agent does not score TAND severity, administer the TAND-L, or assign a diagnosis. It surfaces *where, in the longitudinal record, a cluster shows accumulating but unaddressed signal* — so a clinician can decide whether a formal assessment is warranted. The distinction is enforced in the prompt contracts and again in the surfacing layer (§11.6).

### 11.3 Per-note discourse analysis (Sonnet)

The input is the same time-anchored corpus the Phenome Mapper consumed: synthetic clinical notes shaped from published TSC encounter templates, persistently watermarked, with evidence spans preserved. The TAND agent processes notes individually with Claude Sonnet, one note per call, for two reasons. First, per-note isolation keeps provenance clean: every signal the agent emits points to a single note, a character span, and the model/prompt version that produced it. Second, the relevant judgment is local — it is about how a particular sentence is phrased — so there is no accuracy benefit to cross-note context at this stage, and considerable provenance cost.

For each note, Sonnet performs two passes folded into one structured output:

1. **Cluster-relevant content extraction.** Which of the six clusters does this note touch, and on what spans? A sentence about a child "having a hard time settling at bedtime and during transitions at school" touches Behavioral and Academic.
2. **Discourse-marker classification.** For each touched span, does the language carry a diagnostic-uncertainty marker — and which one? This is where the Marshall-Hagedorn taxonomy lives (§11.4).

The output is strictly typed JSON validated against a Pydantic schema before it is allowed into the event log:

```json
{
  "note_id": "synthnote_000412",
  "patient_id": "B",
  "note_date": "2025-09-18",
  "author_role": "developmental_pediatrics",
  "signals": [
    {
      "cluster": "behavioral",
      "span": {"start": 1142, "end": 1268},
      "text": "Mom reports he's been more irritable and harder to settle; we'll keep an eye on it and revisit at the next visit.",
      "uncertainty_markers": ["third_party_attribution", "follow_up_without_formalization"],
      "marker_confidence": 0.82,
      "addressed": false
    }
  ],
  "provenance": {
    "model": "claude-sonnet-4-5-20250929",
    "prompt_template": "tand_per_note_v3",
    "input_hash": "sha256:…",
    "latency_ms": 2140
  }
}
```

The `addressed` boolean is the hinge. A note can mention an attention problem *and* document a referral to neuropsychology in the same paragraph — that signal is addressed and should not accumulate as a surveillance gap. Sonnet is instructed to mark a signal `addressed: true` only when the same note contains an explicit action (referral placed, assessment ordered, medication started, formal diagnosis recorded). Everything else stays `false` and flows to aggregation. Haiku is not used here; the discourse judgment is too subtle for the fast tier. Opus is reserved for the briefing synthesis (§11.5), not per-note work, to keep the cohort pass affordable.

### 11.4 Deep dive: the Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy

The intellectual core of this agent is a taxonomy of *how clinicians linguistically encode uncertainty and unresolved concern in documentation*. It is drawn directly from a body of published research led from CCHMC's Division of Biomedical Informatics.

The relevant literature anchors:

- **Marshall, Nickels, Brady, Hagedorn (2023), *Journal of Hospital Medicine*** — characterized how diagnostic uncertainty is expressed in inpatient pediatric documentation, establishing that uncertainty is communicated through identifiable linguistic patterns rather than explicit "I am uncertain" statements.
- **Nickels, Marshall, Edgerton, Brady, Hagedorn, Lee (2024), *Applied Linguistics*** — a deeper linguistic treatment of the discourse markers themselves, formalizing the categories of hedging and uncertainty expression in clinical text.
- **Ipsaro, Patel, Marshall, Hagedorn (2021), *Hospital Pediatrics*** — clinical context for diagnostic uncertainty and its documentation in pediatric inpatient care.
- **Orenstein et al., including Hagedorn (2021), *JAMIA*** — on alert burden and the cost of interruptive notifications, which directly informs this agent's anti-alert posture (§11.6).

From this work the agent operationalizes five discourse-marker categories. Each is a label Sonnet must assign by class, with a span, so the briefing can later explain *why* a signal was surfaced in the clinician's own vocabulary.

| Marker | What it sounds like in a note | Why it indicates an unaddressed TAND signal |
|---|---|---|
| **Hedging** | "possible," "appears," "some concern for," "borderline," "mild ?attention difficulty" | The clinician noticed something but stopped short of asserting it — classic pre-diagnostic state. |
| **Deferral** | "will reassess," "monitor for now," "too early to say," "let's see how things go" | The decision is explicitly pushed to a future encounter that may not be the one that closes the loop. |
| **Third-party attribution** | "Mom reports," "teacher notes," "per school," "family concerned about" | The concern enters the record secondhand and is at high risk of never being independently formalized. |
| **Conditional** | "if this continues, we would consider," "should the behaviors worsen, refer to…" | The action is real but gated on a condition no one is systematically tracking. |
| **Follow-up without formalization** | "discussed at length, will revisit," "keep an eye on," "address at next visit" | The intent to act exists but no concrete order, referral, or diagnosis was placed. |

This taxonomy is the right instrument for TAND specifically because TAND under-recognition is overwhelmingly a *follow-up-without-formalization* and *third-party-attribution* problem. Parents and teachers raise developmental and behavioral concerns; those concerns are dutifully recorded; and they are then deferred, visit after visit, by clinicians whose primary cognitive load is the seizure or the renal mass in front of them. The taxonomy gives the agent a principled way to distinguish "this was mentioned and handled" from "this has been mentioned three times across eighteen months and never acted on" — exactly the longitudinal pattern the aggregation layer is built to find.

This is also why the agent is the strongest CCHMC fit in the engine. The Variant Curator borrows the institution's banked tissue; the Phenome Mapper borrows its Clarity data. The TAND agent borrows its *method*. The discourse-marker categories are not invented for this paper — they are the published output of the same division and the same investigator the engine is being built alongside. Hagedorn has taken interest in the work and has offered to engage his team and the faculty TSC lead, and the natural reading is that he would be evaluating an implementation of his own research applied to a disease where its premise is unusually load-bearing.

A scoping note in the anti-overclaim spirit of this paper: the agent applies the *taxonomy and its categories* as published. It does not reproduce any proprietary CCHMC annotated corpus or trained model, which would be institutional Phase-1 work under appropriate data governance. In the synthetic demo the discourse markers are detected zero-shot by Sonnet under a prompt that encodes the five categories with definitions and few-shot examples authored from the public literature. Calibration against a real annotated corpus is explicitly future work, not a claim of the demo.

### 11.5 Deterministic aggregation, scoring, and the Opus briefing

The per-note signals are LLM-derived; what happens to them next is not. The aggregation and scoring layer is ordinary deterministic code operating over the event-sourced state in PostgreSQL. This separation is deliberate and matches the pattern used throughout the engine: a language model makes the local linguistic judgment, and auditable code makes the longitudinal decision about what rises to a clinician's attention.

For each patient and each of the six clusters the aggregator computes a cluster surveillance state from the unaddressed signals:

```python
@dataclass
class ClusterState:
    cluster: str
    unaddressed_signal_count: int        # signals with addressed == False
    distinct_encounters: int             # number of separate visits with a signal
    span_months: float                   # first to most recent unaddressed signal
    distinct_authors: int                # corroboration across clinicians
    dominant_markers: list[str]          # most frequent uncertainty markers
    last_formalization: date | None      # last referral/dx/assessment, if any
    surveillance_gap_score: float        # deterministic composite, 0–1
```

The `surveillance_gap_score` is a transparent weighted function — longer span, more distinct encounters, more distinct authors, and longer time since the last formalization all increase it, while recent formalization sharply discounts it. The exact weights live in version-controlled YAML, not in a model, so the score is reproducible and the demo's behavior is explainable line by line. The design intent: a single hedged mention should never surface; a concern corroborated across three authors over a year with no formalization should.

Only after the deterministic layer has decided *which* cluster states cross the surfacing threshold does Opus enter. Opus does one job: turn a small set of qualifying cluster states into pre-visit **briefing prose** — a few sentences per cluster that name the pattern, quote the clinician's own language with dates and authors, and stop. Opus is constrained to write from the structured `ClusterState` and the underlying signal spans only; it may not introduce clinical claims, may not suggest a diagnosis, and must attribute every assertion to a dated note. The output is framed as material a clinician reads in the ninety seconds before walking into the room.

For Patient B (the 12-year-old with a slow-growing SEGA and bilateral AMLs), the synthetic record carries scattered, under-recognized TAND signals deliberately seeded across encounters — mild attention concerns raised by a teacher, intermittent sleep and irritability reports from a parent, none ever formalized. The agent's success criterion on Patient B is that the briefing surfaces these as a coherent academic/behavioral surveillance gap, in the clinicians' own words, without inventing severity or diagnosis. Patient B is the demo's proof that the method finds the quiet thing in the chart.

### 11.6 Surfaced as briefing, never as an alert

The single most important product decision in this agent is what it is *not*: it is not an alert. The output reaches the clinician only through the pre-visit briefing surface (the one-screen, mobile-readable header / what's-new / 0–3 action items / watchlist / links layout) and, where appropriate, the in-visit dashboard's TAND quadrant. It never fires an interruptive notification.

This is a direct consequence of the Orenstein/Hagedorn (2021) JAMIA work on alert burden. The fastest way to make a TAND-surveillance tool clinically worthless is to turn it into another interruptive popup that clinicians learn to dismiss. So the discipline is strict:

- TAND signals appear as **watchlist** items and, at most, as one of the briefing's 0–3 action items when a cluster's gap score is high and persistent.
- The async alert surface's strict-discipline rule applies (recalibrate if the system is producing more than roughly three alerts per clinician per week); TAND rarely warrants an async alert at all and defaults to the briefing.
- The language is always "this pattern is present in the record; consider whether a formal TAND assessment is warranted," never "patient has TAND" and never "refer now."

The clinician remains the gate. The agent's contribution is to compress eighteen months of distributed, hedged, secondhand mentions into something a busy specialist can absorb before a visit and act on if they judge it warranted. That is decision-support, and it is the only posture this agent is permitted to take.

### 11.7 Evaluation against synthetic ground truth, and limitations

Evaluation uses the synthetic cohort's known seeded signals as ground truth — this is build validation against a controlled cohort, not clinical validation. The targets:

- **Detection.** The agent must surface the TAND signals deliberately embedded in Patient B's record and across the cohort's seeded cases, mapped to the correct cluster.
- **No spurious flags.** The agent must not manufacture TAND surveillance gaps where the seeded record contains only addressed or absent signals. Given the briefing-not-alert posture, a false surveillance gap is the costliest error, and the deterministic threshold is tuned conservatively against it.
- **Provenance integrity.** Every surfaced signal must trace to a note, a span, a marker class, and a model/prompt version, queryable from the append-only event log.

Limitations stated plainly, in keeping with the rest of this paper:

1. **The discourse markers are detected zero-shot, not corpus-calibrated.** Sonnet applies the published taxonomy under a prompt; it has not been tuned against a labeled CCHMC corpus. Real calibration, with inter-annotator agreement against the literature's own annotation scheme, is institutional Phase-1 work and would run against real Clarity notes under IRB. The demo demonstrates the *mechanism*, not a validated performance figure on real documentation.
2. **Synthetic notes are easier than real notes.** They are coherent, templated, and free of the messiness — copy-forward, conflicting authors, OCR'd scans — that real charts carry. Detection numbers on the synthetic cohort should be read as an upper-bound feasibility signal, not an estimate of real-world recall.
3. **The taxonomy is necessary but not sufficient.** A clinician can fail to address a concern they never wrote down at all; this agent can only see what is in the record. It reduces the documentation-driven slice of TAND under-recognition, which the literature suggests is large, but it does not claim to close the gap.
4. **It is not a diagnostic instrument.** It does not administer the TAND-L, does not produce severity scores, and does not replace neuropsychological assessment. It tells a clinician where to look.

With the longitudinal phenotype mapped, the trajectory of measurable disease forecast, and the quiet neuropsychiatric signal surfaced from the clinician's own words, the engine has assembled a complete picture of where a TSC patient is and where they are heading. Section 12 turns that picture into a structured set of therapeutic options — the work of the TSC-Therapeutics Strategist.


## 12. Agent 5 — TSC-Therapeutics Strategist

Where Agent 4 surfaces patterns a clinician might otherwise miss, Agent 5 answers a harder question: given everything now known about this patient, what are the defensible therapeutic options, and what does the literature actually say about each? The TSC-Therapeutics Strategist is the integrative endpoint of the engine. It does not produce new measurements. It consumes the outputs of the four prior agents, layers in medications, adherence, and adverse-event history, retrieves the relevant external evidence, and synthesizes a structured options brief that a treating physician can read, interrogate, and override. It is the only agent in the cluster that touches the question of what to *do*, and for that reason it is the agent held to the highest standard of restraint.

### 12.1 Why Opus is non-negotiable here

Within the engine's per-step tiering policy (§16) — deterministic code wherever the task is rule-shaped, Opus reserved for the few tasks where the cost of a subtly wrong synthesis is high — the Therapeutics Strategist is the clearest case for Opus. It is the one place in the system where the model must hold five heterogeneous evidence streams in working context simultaneously — a variant interpretation, a longitudinal phenotype, a set of quantitative trajectory forecasts with prediction intervals, a TAND briefing, and a medication/adverse-event history — and reason over their *interactions* (an mTOR inhibitor that helps the SEGA but whose mucositis intersects an existing adherence problem; a trial whose inclusion criteria conflict with a renal trajectory crossing a threshold). That is integrative clinical reasoning under uncertainty, and it is exactly the regime where weaker models produce fluent text that is locally plausible and globally wrong.

So Agent 5 runs on Claude Opus, non-negotiably, with the local Llama 3.1 70B Instruct (Ollama) path available only as a degraded fallback that is explicitly labeled as such in the provenance record and never used for a delivered brief without a banner. This is the most expensive call in the engine per invocation, and it is invoked the least often — once per patient per synthesis event, demand-driven, not on every data change. The orchestrator (§13) only triggers a Therapeutics synthesis when an upstream output materially changes (a new variant classification, a trajectory threshold-crossing alert, a new adverse event), which keeps the Opus spend bounded.

### 12.2 Inputs: the five streams plus three

The Strategist's context is assembled deterministically by the orchestrator before the Opus call, never gathered by the model itself. The assembled context for a patient is a structured object, every field carrying its own provenance:

| Stream | Source | What it contributes |
|---|---|---|
| Variant interpretation | Agent 1 (Variant Curator) | Gene (TSC1/TSC2), variant, ACMG-AMP class, mosaic flag, contiguous-deletion flag (PKD1) |
| Longitudinal phenotype | Agent 2 (Phenome Mapper) | Time-anchored HPO profile, organ involvement, surveillance-gap report |
| Trajectory forecasts | Agent 3 (Trajectory Modeler) | SEGA/AML growth and seizure-burden forecasts at 6/12/18 mo with 50%/90% intervals, threshold alerts |
| TAND briefing | Agent 4 (TAND Surveillance) | Surfaced TAND-cluster signals as briefing material (never as diagnoses) |
| Medications & AEs | Synthetic structured data | Current/prior agents (everolimus, sirolimus, AEDs), dose history, adherence flags, documented adverse events |
| PubMed/PMC RAG | Milvus TSC corpus partition | Retrieved literature passages with URIs (mTOR-inhibitor efficacy, EXIST-3, next-gen mTORC1, surveillance) |
| ClinicalTrials.gov snapshot | Version-pinned local snapshot | Open TSC-relevant trials with structured eligibility criteria |
| FDA actions | Version-pinned local snapshot | Approvals, label changes, safety communications relevant to TSC therapeutics |

Two points of discipline. First, the ClinicalTrials.gov and FDA "snapshots" are exactly that — a dated, version-pinned local copy retrieved as part of the synthetic-demo build, not a live query. The brief states the snapshot date on its face, and the demo makes no claim of real-time trial currency. Second, the RAG retrieval runs against the TSC corpus partition in Milvus (BAAI/bge-large-en-v1.5 plus the BiomedBERT-derived clinical embeddings described in §3), so retrieved passages are TSC-scoped and carry resolvable URIs back to the source abstract or full-text section.

### 12.3 Output: the six-section options brief

The Strategist produces a single structured artifact with a fixed six-section schema. The fixed structure is itself a credibility device: a clinician learns the shape once and can navigate any patient's brief, and a fixed schema makes it far harder for the model to bury a weak claim in unstructured prose. Every section is constrained to source-attributed statements; a sentence that asserts a clinical fact without an attached source reference is a schema violation caught by the post-generation validator, not a stylistic preference.

1. **Current Therapy.** A factual restatement of what the patient is on, at what dose, with what documented response and tolerability. No recommendation — this section grounds the rest.
2. **Optimization.** Options within the current regimen: dose titration toward therapeutic everolimus troughs, AE management (e.g., mucositis prophylaxis), adherence support. Each framed against the trajectory forecasts and the literature.
3. **Combination.** Where evidence supports adding a modality (e.g., mTOR inhibitor alongside an AED for refractory seizures, with reference to EXIST-3), the option is laid out with its supporting and contradicting evidence.
4. **Trial Matching.** Eligibility-criterion-level matching against the ClinicalTrials.gov snapshot (detailed in §12.4).
5. **Emerging Evidence.** Literature and pipeline signals not yet standard of care — e.g., next-generation selective mTORC1 inhibitors in development — flagged explicitly as not-yet-established.
6. **Open Questions.** What the engine cannot resolve: data the brief lacks, conflicting evidence, decisions that genuinely require the treating physician's judgment. This section is mandatory and must be non-empty; a brief with no open questions is a red flag, not a triumph.

Each section header carries a one-line confidence/sourcing summary, and the whole brief is wrapped in the standard provenance envelope (model id and version, prompt-template version, full list of retrieved RAG URIs, input hash over the assembled context, and end-to-end latency).

```json
{
  "patient_id": "C",
  "synthesis_event_id": "evt_7f3a...",
  "generated_at": "2026-07-15T14:22:09Z",
  "snapshot_dates": { "clinicaltrials_gov": "2026-06-30", "fda": "2026-06-30" },
  "provenance": {
    "model": "claude-opus-4-x",
    "prompt_template": "therapeutics_brief.v3",
    "input_hash": "sha256:9b1c...",
    "rag_sources": ["pmc:PMC1234567#results", "ct:NCT0XXXXXXX"],
    "latency_ms": 168400,
    "fallback_used": false
  },
  "sections": {
    "current_therapy": [ { "text": "...", "sources": ["..."] } ],
    "optimization": [ ... ],
    "combination": [ ... ],
    "trial_matching": [ ... ],
    "emerging_evidence": [ ... ],
    "open_questions": [ ... ]
  }
}
```

### 12.4 Trial matching at the eligibility-criterion level

Naive trial matching — "this patient has TSC, here are TSC trials" — is worse than useless to a clinician because it generates a list they then have to manually adjudicate. The Strategist matches at the level of the individual eligibility criterion. Each open trial in the snapshot carries its inclusion and exclusion criteria as discrete, structured items. For each criterion, the engine emits one of three statuses, and crucially, the *reason*:

- **eligible** — the patient's data satisfies the criterion, with the satisfying evidence cited (e.g., "TSC2 pathogenic variant confirmed — Agent 1, ACMG class").
- **ineligible** — the patient's data violates the criterion, with the violating evidence cited (e.g., "prior everolimus exposure" against an mTOR-inhibitor-naive inclusion criterion).
- **requires-clarification** — the engine cannot determine status from available data (e.g., a criterion requiring a specific echocardiographic parameter the synthetic record does not contain). This is treated as a first-class outcome, never silently coerced to eligible or ineligible.

The mechanics are deliberately split between deterministic code and the model. Structured criteria that map cleanly to structured patient data (age band, gene, prior-therapy flags, lab thresholds) are adjudicated by deterministic rules — these are auditable and never hallucinated. Free-text criteria that require interpretation against the longitudinal record ("clinically significant TAND," "stable disease for ≥3 months") are passed to Opus, which must return both a status and a quoted span of patient evidence supporting it. A trial is only ever surfaced as a candidate in the brief if no criterion is *ineligible*; trials with one or more `requires-clarification` criteria are surfaced with those criteria flagged as the specific things the clinician must verify. The brief never says "the patient qualifies." It says "the patient appears to meet criteria 1–6; criteria 7 and 8 require clarification of X and Y," which is the form a clinician can act on.

### 12.5 Source attribution and the decision-support frame

Two non-negotiable framing constraints govern every word of the output, and both are enforced mechanically rather than left to model temperament.

First, **every clinical claim is source-attributed.** A statement about mTOR-inhibitor efficacy carries the PubMed/PMC URI of the passage that supports it; a statement about a trial carries its NCT identifier; a statement about FDA status carries the snapshot reference. Claims that derive from the patient's own data carry the upstream agent and field as their source. The post-generation validator parses the brief, confirms that every assertion in sections 2–5 has at least one resolvable source attached, and rejects the brief back to a regeneration with a stricter prompt if any does not. This is the same provenance discipline applied across the engine, but it bites hardest here because this is where unsupported claims would be most dangerous.

Second, **the brief is decision-support, never a recommendation.** The language is constrained to options and evidence — "option," "supported by," "would require," "the literature reports," "not yet established." It does not say "you should," "the recommended therapy is," or "start." The Open Questions section exists precisely to make the residual judgment visible rather than papered over. The treating physician is the decision-maker; the engine's job is to make sure the relevant evidence and the relevant gaps are in front of them, organized and cited. Like every output in the engine, the brief is persistently watermarked as synthetic-data, AI-generated, decision-support material, and surfaces (§13) present it as such.

### 12.6 Worked example: Patient C

Patient C is the synthetic 18-year-old female with a germline TSC1 variant, a partial everolimus response complicated by mucositis that forced a dose reduction, a roughly 4 cm renal AML, and refractory focal seizures. She is the cohort's deliberate "hard therapeutics" case, and her brief exercises every section.

- **Current Therapy.** Everolimus at a reduced dose following grade-2 mucositis; partial response noted in the AML trajectory; refractory focal seizures on her current AED regimen. Drawn from medications/AE history and Agents 2–3, fully attributed.
- **Optimization.** Two literature-supported options: (a) mucositis prophylaxis (dexamethasone mouthwash) to permit re-titration toward a therapeutic everolimus trough, since the dose reduction is the plausible driver of the partial AML response; (b) trough-guided titration with surveillance cadence informed by Agent 3's AML forecast. The ~4 cm AML sits near the bleeding-risk threshold, which the brief flags as raising the stakes of the optimization question — drawn directly from the Trajectory Modeler, with its prediction interval shown.
- **Combination.** The refractory seizures plus an existing mTOR inhibitor make adjunctive mTOR-pathway evidence relevant; the section lays out the EXIST-3 evidence for everolimus as adjunctive therapy in TSC-associated seizures, with the citation, and notes the patient is already on an mTOR inhibitor so the question is optimization of the existing agent rather than addition.
- **Trial Matching.** Criterion-level matching against the snapshot. A next-generation selective mTORC1 inhibitor trial is surfaced; the brief shows she is *eligible* on gene, age, and TSC diagnosis criteria, but flags *requires-clarification* on a washout criterion (her current everolimus exposure may or may not satisfy the trial's washout definition) — exactly the kind of item the engine refuses to guess.
- **Emerging Evidence.** Next-generation selective mTORC1 inhibitors are described as in development and not yet standard of care, explicitly labeled emerging.
- **Open Questions.** Whether mucositis-limited dosing can be overcome before considering a trial; the bleeding-risk decision on the ~4 cm AML, which is a treating-physician call the engine does not make; and the washout ambiguity from the trial-matching section.

The brief gives the clinician a navigable map of Patient C's therapeutic landscape with every claim cited and every gap named. It does not tell them what to do.

### 12.7 Evaluation against synthetic ground truth

As with every agent in this engine, the targets below are demo-quality checks against the synthetic cohort's ground truth, not clinical validation. The cohort was authored with known therapeutic situations, so the Strategist's output can be scored against intended answers.

- **Trial matching correctness.** For the patients with authored trial-eligibility ground truth (notably Patient C), the engine assigns the correct eligible/ineligible/requires-clarification status to each criterion, with no criterion that should be *ineligible* marked eligible (the safety-critical direction).
- **Appropriate hedging.** Emerging and not-yet-established options are flagged as such; the brief contains no overclaimed efficacy. Reviewed by sampling against authored expectations.
- **Full attribution.** 100% of clinical claims in sections 2–5 carry a resolvable source. Enforced by the validator and confirmed in eval as a hard gate.
- **Latency.** Synthesis completes in under 3 minutes per patient on the build target (Opus via API; the heavier upstream variant-calling and cohort-generation work is what consumes the RunPod GPU budget, not this call).

### 12.8 Limitations

The honest boundaries on this agent are sharp and worth stating plainly. The Strategist reasons only over what the engine knows: it has no pharmacy or claims feed, no patient-reported outcomes, no neuropsychological test scores, and no real imaging — all explicitly absent from the synthetic cohort. Its trial and FDA views are pinned snapshots, not live data, and the brief says so on its face; real-time currency is institutional Phase-1 work, not part of the demo. The RAG corpus is a curated TSC partition, so genuinely novel evidence outside that partition will not surface. And the integrative reasoning, while Opus-class, is precisely the place where a fluent-but-wrong synthesis is most plausible — which is why the human gate is absolute here. The brief is a draft for a treating physician, every claim cited, every gap named, and nothing in it is autonomous. None of these constraints are bugs to be apologized for; in this audience they are the reason the brief is worth reading at all.

With all five agents now specified, what remains is the connective tissue that makes them behave as one engine rather than five tools — the deterministic orchestrator that sequences enrollment, minimizes recomputation, and decides when a Therapeutics synthesis like Patient C's is even warranted, and the three clinician surfaces on which all of this is finally presented. That is §13.


## 13. The Orchestrator and Clinician-Facing Surfaces

The five agents described in the preceding sections each produce structured output in isolation. The Engine is not five agents; it is five agents plus the connective tissue that decides when each one runs, what it consumes, what depends on it, and how its output is composed into something a clinician actually looks at during a fifteen-minute encounter. That connective tissue is two distinct layers: a deterministic orchestrator that governs execution and state, and three clinician-facing surfaces that govern presentation. This section specifies both. The design principle that runs through all of it is conservatism — the orchestrator never guesses, never silently drops output, and never lets a model decide control flow; the surfaces never interrupt, never overstate, and never present a forecast without its provenance. With Hagedorn's audience, the failure mode that ends the conversation is a system that looks confident when it should look uncertain. The orchestrator and surfaces are where that discipline is enforced mechanically rather than promised rhetorically.

### 13.1 Why the orchestrator is deterministic

A reasonable first instinct, given that four of the five agents are LLM-backed, is to make the coordinator an LLM too — an "agent of agents" that reads the patient state, reasons about what to do next, and dispatches the sub-agents. We reject that design explicitly. The orchestrator is a deterministic LangGraph event router. It contains no model call. Given the same event log it produces the same execution, every time, and that property is the point.

There are three reasons. First, **auditability**. Every section of this paper that touches the demo audience returns to the same theme: the build is the argument, and the argument only holds if a skeptical reviewer can replay it. A model-driven controller would make "why did the Engine run the Trajectory Modeler before the TAND agent on this patient?" a question with a probabilistic, prompt-dependent, sometimes-irreproducible answer. A deterministic router makes it a question answerable by reading 200 lines of Python and the event log. Second, **failure semantics**. When an LLM agent times out or returns malformed output, the system must do something sane and predictable. A deterministic router can encode "if the Variant Curator fails, mark the variant quadrant pending and proceed" as an explicit, tested branch. A model asked to handle the same failure might improvise, hallucinate a partial result, or stall. Third, **cost and latency**. The control plane fires on every event; putting a model in that loop would add tokens, dollars, and seconds to bookkeeping that is pure routing. The expensive cognition belongs in the agents, where it is gated, provenance-tracked, and reviewed. The plumbing should be cheap and boring.

This is the same instinct that put classical statistics rather than an LLM inside the Trajectory Modeler. Where a problem is well-specified — routing, forecasting from longitudinal measurements — we use the deterministic tool that is correct and inspectable, and reserve model cognition for the genuinely ambiguous synthesis tasks (ACMG classification, discourse analysis, therapeutic options framing). Determinism is not a limitation we tolerate; it is a feature we chose.

### 13.2 Event model

The orchestrator is event-sourced. The authoritative record of what has happened to a patient is an append-only log of events in PostgreSQL; the "current state" a surface reads is a materialized projection derived from that log. Nothing in the system mutates patient state in place. To correct a value you append a corrective event; the history remains. This buys two things at once: a complete audit trail for free (the log *is* the audit trail), and trivial replay (re-run the projection over the log to reconstruct any past state, or to test a code change against recorded history).

Thirteen event types span the lifecycle from enrollment through surface assembly. They fall into four families: lifecycle, data-ingest, agent, and surface.

| # | Event type | Family | Emitted by | Triggers / effect |
|---|-----------|--------|-----------|-------------------|
| 1 | `PATIENT_ENROLLED` | lifecycle | enrollment driver | Registers patient; kicks off dependency-ordered agent scheduling (Phenome Mapper first) |
| 2 | `PATIENT_DEACTIVATED` | lifecycle | enrollment driver | Halts further scheduling; state retained, surfaces frozen |
| 3 | `GENOMIC_DATA_AVAILABLE` | data-ingest | cohort loader | New/updated BAM or VCF registered for patient; eligible input for Variant Curator |
| 4 | `CLINICAL_DATA_UPDATED` | data-ingest | cohort loader | New notes / Clarity-shaped structured rows; eligible input for Phenome Mapper, TAND |
| 5 | `IMAGING_REPORT_AVAILABLE` | data-ingest | cohort loader | New imaging report (synthetic); eligible input for Phenome Mapper, Trajectory Modeler |
| 6 | `AGENT_RUN_REQUESTED` | agent | orchestrator | Router determines an agent's inputs are ready and its output is stale; enqueues a run |
| 7 | `AGENT_RUN_STARTED` | agent | agent worker | Marks run in-flight; records model id/version, prompt template version, input hash, start time |
| 8 | `AGENT_OUTPUT_PRODUCED` | agent | agent worker | Carries the structured output + full provenance; updates the patient-state projection |
| 9 | `AGENT_RUN_FAILED` | agent | agent worker | Records error, latency, retry count; triggers conservative failure handling for that agent's surface region |
| 10 | `THRESHOLD_CROSSED` | agent | Trajectory Modeler | Forecast crosses a clinical threshold (e.g., SEGA approaching surgical consideration); candidate for alert surface |
| 11 | `BRIEFING_REQUESTED` | surface | clinician action / scheduler | Demand-driven assembly of the pre-visit briefing for an upcoming encounter |
| 12 | `DASHBOARD_OPENED` | surface | clinician action | Demand-driven assembly / refresh of the in-visit 4-quadrant dashboard |
| 13 | `ALERT_RAISED` | surface | orchestrator | A `THRESHOLD_CROSSED` or agent finding passes alert-discipline gating and is posted to the async alert surface |

A subtlety worth naming: `THRESHOLD_CROSSED` and `ALERT_RAISED` are deliberately separate. An agent *detecting* a threshold crossing is a fact about the data and is always logged. An alert being *raised to a clinician* is a decision about attention, made by the orchestrator under the discipline rule of §13.6. Decoupling them means the system can record everything it noticed while showing the clinician only what is worth an interruption — and an auditor can later see the difference between the two.

### 13.3 State and storage

Three stores, each with a distinct job, and no overlap in responsibility:

- **PostgreSQL — durable truth.** Two logical pieces. (1) The append-only `event_log` table: one row per event, with `event_id`, `patient_id`, `event_type`, monotonically-increasing `sequence`, `payload` (JSONB), `provenance` (JSONB), and `created_at`. Append-only is enforced at the application layer and reinforced by withholding UPDATE/DELETE grants on the table. (2) Materialized current-state projections — per-patient, per-agent latest-output tables (e.g., `current_variant_interpretation`, `current_hpo_profile`, `current_trajectory_forecast`) rebuilt by folding the log. Projections are derived and disposable; the log is canonical.
- **Redis — ephemeral coordination.** In-flight run locks (so a patient's Phenome Mapper does not run twice concurrently), surface-assembly caches, and short-lived computed artifacts, all under explicit TTLs. Nothing in Redis is a source of truth; flushing it loses nothing that cannot be rebuilt from PostgreSQL. This separation keeps the durable store clean of transient coordination noise.
- **YAML demo config — declared scenario.** The 50-patient cohort, the featured Patients A/B/C, the demo clock, alert thresholds, and which surfaces are pre-warmed are all declared in version-controlled YAML rather than baked into code. The demo is reproducible because its configuration is a file in the repo, not state accumulated by clicking around. Re-pointing the Engine at a different cohort is a config edit, not a code change — the same property that makes "swap the box labels, keep the wiring" replication to another institution tractable.

### 13.4 Orchestration patterns

Four patterns govern how the router turns events into agent runs and surfaces. Each is a small, testable rule, not an emergent behavior.

**Dependency-ordered enrollment.** When `PATIENT_ENROLLED` fires, agents are not all dispatched at once. The Phenome Mapper runs first, because its longitudinal HPO profile is the foundation the Trajectory Modeler and TAND agent build on, and because the Therapeutics Strategist integrates all four prior agents. The dependency graph the router enforces:

```
PATIENT_ENROLLED
        │
        ▼
  Phenome Mapper ──────────────┐
        │                      │
        ├──► Trajectory Modeler │   (Variant Curator runs in
        │            │          │    parallel, gated only on
        ├──► TAND Agent          │    GENOMIC_DATA_AVAILABLE)
        │            │          ▼
        └────────────┴──► Therapeutics Strategist
```

The Variant Curator is the one agent not downstream of the Phenome Mapper; it depends only on genomic data availability and runs in parallel. The Therapeutics Strategist is the join point — the router will not request its run until the four upstream outputs it integrates are present (or explicitly marked unavailable). This ordering is encoded once, in the router's dependency table, and exercised by every enrollment in the cohort.

**Incremental-update minimization.** Re-running every agent on every data change would be wasteful and, for the Opus-class Therapeutics Strategist, expensive. The router consults the input hash recorded on each agent's last `AGENT_OUTPUT_PRODUCED`. A new `CLINICAL_DATA_UPDATED` event recomputes the hash of the Phenome Mapper's inputs; if it is unchanged (e.g., an administrative note with no extractable phenotype), no run is requested. If the Phenome Mapper output changes, the router propagates: downstream agents whose inputs now differ are re-requested; those unaffected are not. The unit of work is "what actually changed," not "everything, on a timer."

**Demand-driven surface assembly.** Surfaces are assembled when requested, not maintained continuously. `BRIEFING_REQUESTED` and `DASHBOARD_OPENED` pull the current-state projections and compose the surface on demand, caching the result briefly in Redis. This means a surface always reflects the latest available agent output without a background process churning to keep stale screens warm, and it means surface assembly logic is exercised exactly when a clinician looks — the only moment it matters.

**Conservative failure handling.** When `AGENT_RUN_FAILED` fires, the affected region of every surface that would have shown that agent's output renders an explicit state — `pending`, or a staleness indicator with the timestamp of the last good output — never a silent omission. A clinician must never be shown a four-quadrant dashboard where the trajectory quadrant is blank because the modeler crashed and the reader is left to assume "nothing to report." The distinction between "the Engine has nothing to say here" and "the Engine could not run here" is preserved on screen. The router retries failed runs a bounded number of times with backoff; persistent failure surfaces as visible staleness, logged with full error provenance for replay.

### 13.5 The three surfaces

The Engine ships three standalone web applications, not Epic embeddings. They read the current-state projections and are persistently watermarked synthetic — every screen carries the watermark so no screenshot can be mistaken for a real patient record. Each surface has a different relationship to the clinician's attention, and that relationship dictates its structure.

#### Pre-visit briefing — earn ten seconds

The briefing is built for the moment before a clinician walks into the room: scanned on a phone or in a hallway, in well under a minute. It is one screen, mobile-readable, with a fixed five-part structure:

1. **Header** — patient identifier (synthetic), age, TSC genotype in one line (e.g., *TSC2 c.3037C>T p.Arg1013Ter, germline*), upcoming visit context.
2. **What's new** — the delta since the last encounter, in plain language. The point of the briefing is that the clinician should not re-read the whole chart; the Engine tells them what changed.
3. **Action items (0–3, hard cap)** — concrete, optional next steps. Zero is a valid and common answer. The cap is a forcing function: if everything is an action item, nothing is.
4. **Watchlist** — trajectories and findings being monitored that do not yet warrant action (e.g., a SEGA growing within expected range), with their forecast context.
5. **Links** — deep links into the in-visit dashboard for any item the clinician wants to open fully.

The demo configures the briefing for Patient B: the SEGA series (0.8 → 1.1 → 1.3 cm at the foramen of Monro) appears on the watchlist with the Trajectory Modeler's forecast, and the scattered under-recognized TAND signals appear as briefing material, never as a diagnosis. This is the surface that makes the "briefing material, not alert" posture concrete — the TAND agent's findings reach the clinician here, pre-visit, framed for consideration, and not as an interruptive pop-up during care.

#### In-visit four-quadrant dashboard — the whole picture, drillable

The dashboard is the surface for the encounter itself, opened on `DASHBOARD_OPENED`. It lays the Engine's full output across four quadrants on one view, with progressive disclosure so the top level is scannable and any element drills to its evidence.

| Quadrant | Source agent | Top-level content | Drill-down |
|----------|-------------|-------------------|-----------|
| Upper-left | Variant Curator | Variant interpretation: genotype, ACMG-AMP class, mosaic flag if present | ClinVar-spec detail, VAF/reads/strand/artifact evidence, ddPCR validation rec, AI-labeled draft report, full provenance |
| Upper-right | Phenome Mapper | Longitudinal HPO timeline, organ-system roll-up | Per-phenotype evidence spans linked to source notes, discordance log, ITSC surveillance-gap report |
| Lower-left | Trajectory Modeler | SEGA/AML growth, seizure burden, renal function forecasts at 6/12/18 mo with 50%/90% prediction intervals | Model family used, threshold-crossing detail, surveillance-cadence recommendation |
| Lower-right | TAND + Therapeutics | TAND briefing summary across the six clusters; six-section therapeutic options brief | Marshall-Hagedorn discourse evidence per TAND signal; per-claim source attribution for every therapeutic option |

Every quadrant supports source navigation: a clinician can move from a displayed conclusion to the underlying evidence and the provenance record (model id/version, prompt template version, retrieved RAG sources with URIs, input hash, latency) in one or two clicks. The dashboard never shows a conclusion the reader cannot trace. The demo drives this quadrant set with Patient B for the variant/HPO/trajectory/TAND picture and with Patient C for the lower-right therapeutics brief — Patient C being the 18-year-old with a partial everolimus response, mucositis-driven dose reduction, an ~4 cm AML, and refractory focal seizures, where the six-section brief (Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions) does its most interesting work.

#### Async alert surface — interruption is a budget

The third surface is for findings that arrive between encounters and may warrant attention before the next scheduled visit. It is deliberately the most disciplined of the three, because an interruptive system that cries wolf is worse than no system. `ALERT_RAISED` events post here, in four categories:

| Category | Example trigger | Typical source |
|----------|----------------|----------------|
| Threshold approach | Forecast indicates SEGA approaching a size where surgical/dosing reconsideration is discussed | Trajectory Modeler |
| Surveillance gap | An ITSC-recommended assessment is overdue given the patient's profile | Phenome Mapper |
| New finding requiring review | A draft variant interpretation ready for molecular-geneticist sign-off | Variant Curator |
| Emerging evidence | A newly indexed trial or FDA action materially relevant to the patient's therapy | Therapeutics Strategist |

#### 13.6 The alert-discipline rule

The single most important number governing the alert surface is the rate ceiling: **no more than roughly three alerts per clinician per week.** This is not a soft preference; it is a calibration target the orchestrator enforces. If alert volume for a clinician trends above that ceiling, the gating thresholds are recalibrated upward — the bar to raise an alert rises — rather than the clinician being asked to absorb more interruptions. The lineage here is direct: this is Orenstein/.../Hagedorn 2021 (JAMIA) on alert burden carried into the design. A surveillance system that desensitizes its clinicians has failed regardless of how good its detection is. The decoupling of `THRESHOLD_CROSSED` (always logged) from `ALERT_RAISED` (gated) exists precisely so the Engine can be sensitive in what it notices and conservative in what it surfaces. Everything noticed remains available in the dashboard and briefing on the clinician's own schedule; only the genuinely time-sensitive minority earns an interruption.

Two further disciplines apply across all three surfaces. First, the human gate is structural: the Variant Curator's output reaches the clinician as a *draft for molecular-geneticist sign-off*, the TAND agent's output as *briefing material*, the Therapeutics Strategist's as a *decision-support options brief* — never as a recommendation the surface presents as settled. The surfaces are built to make the augment-not-replace posture visible in their wording, not merely true in their architecture. Second, the synthetic watermark is persistent and non-removable on every view, and the Epic/Clarity/LIMS integration that would feed these surfaces from real institutional data is, as elsewhere in this paper, described architecturally but explicitly not built in the demo — it is institutional Phase-1 work. What runs now on the Spark reads the synthetic cohort's projections and nothing else.

Taken together, the orchestrator and surfaces are where the Engine's stated values stop being adjectives and become mechanisms: determinism as replayable Python, conservatism as visible `pending` states, augment-not-replace as draft-for-review wording, and respect for clinician attention as a hard rate ceiling tied to published alert-burden research. The remaining question is the substrate all of this executes on — the single DGX Spark that runs the demo and the RunPod GPUs it reaches for when variant calling and cohort generation exceed one box. That is the subject of §14.


---

# Part IV — Building on the DGX Spark + RunPod


## 14. Building It on the DGX Spark (with RunPod for Extra GPU)

The preceding sections described what the five agents and the deterministic orchestrator *do*. This one answers the load-bearing question any skeptical informatics reader asks next: where does this actually run, what does it cost, and what happens when you move from a 50-patient synthetic demo to a real institutional cohort? The honest answer is that the entire demo — every agent, the orchestrator, the RAG layer, the three clinician surfaces — runs first on a single workstation that costs $4,699, and bursts to rented GPU only for a small set of clearly-bounded, embarrassingly-parallel jobs. The build is the argument. If it needed a datacenter to demonstrate, the demonstration would be making a different claim than the one we want to make.

### 14.1 The box on the desk

The primary compute substrate is an NVIDIA DGX Spark: a GB10 Grace Blackwell superchip, roughly 1,000 TOPS of AI compute, 128 GB of unified LPDDR5x memory shared coherently between the Grace ARM CPU and the Blackwell GPU over NVLink-C2C, a 4 TB NVMe local disk, and DGX OS. The unified-memory architecture matters more than the raw TOPS number for this workload. A 70B-parameter local model quantized to 4-bit needs roughly 40-45 GB of weights resident plus KV cache; on a discrete-GPU laptop or a 24 GB card that is a non-starter, but on the Spark the model lives in the same 128 GB pool the rest of the stack draws from, with no PCIe copy across a host/device boundary. That single design choice is what lets a sub-$5k device hold a 70B model *and* a vector database *and* a Postgres instance *and* the agent processes at the same time.

What the Spark runs locally, continuously, with no external dependency other than the Claude API:

| Component | Role | Resident footprint (approx.) |
|---|---|---|
| LangGraph orchestrator + 5 agent processes | TSC-Orchestrator event router + the five agents (FastAPI) | 2-4 GB |
| Milvus 2.4 | TSC RAG corpus partition (PubMed/PMC abstracts, ITSC 2021, TAND consensus, ClinicalTrials.gov snapshot) embedded with `BAAI/bge-large-en-v1.5` + BiomedBERT-derived clinical embeddings | 6-10 GB depending on corpus size |
| PostgreSQL 16 | Append-only event log + materialized current-state projections | 1-2 GB working set |
| Redis 7 | Ephemeral per-event state (TTL'd) | <1 GB |
| MinIO | BAMs, VCFs, synthetic notes/imaging reports, provenance blobs, watermarked artifacts | disk-bound (4 TB NVMe) |
| Ollama: Llama 3.1 70B Instruct (4-bit) | Local-LLM fallback for Haiku/Sonnet-tier work when offline or cost-capped | 40-45 GB when loaded |
| Streamlit | The three clinician surfaces (briefing / dashboard / alert) | <1 GB |

This is the HCLS AI Factory v1.3.0 substrate, reused as-is. Nothing in the list above is net-new for TSC; the net-new code is the multi-agent orchestrator and its event-sourced shared state, the synthetic-cohort pipeline, and the five agents themselves. The point of reuse is not thrift for its own sake — it is that the operational characteristics (how Milvus is partitioned, how provenance is written, how the model tiers are routed) are already understood and already exercised by the six other engines on the platform, so the TSC work inherits a known-good baseline rather than inventing infrastructure under demo pressure.

Two clarifications keep this section honest. First, the model tiers used in the *demo* are predominantly the Claude API (Haiku, Sonnet, Opus), not the local Llama. Llama 3.1 70B is the fallback path — for offline operation, for cost-capping a long batch run, or for the Phase-1 institutional posture where a site may not want PHI-adjacent text leaving the building. The demo runs on Claude because the eval targets in §13 were set against Claude-tier reasoning, and the ACMG-AMP synthesis (Opus, validated against combinatorial rules) and the six-section Therapeutics brief (Opus, non-negotiable) are where the reasoning quality has to be highest. Second, the Spark is a development and demonstration device, not the deployment target. The same stack — LangGraph, Milvus, Postgres, Redis, MinIO, tiered models — is exactly what runs in a the platform-backed institutional deployment, just with the storage layer and the vector index sitting on shared infrastructure instead of the Spark's local NVMe. "Runs on the box on my desk" and "runs in the Winslow Pavilion's data path" are the same software with a different storage substrate underneath; that is the deliberate design, and it is why the demo is a credible preview of the institutional build rather than a throwaway prototype.

### 14.2 When the Spark is not enough, and exactly why

The Spark handles the *steady-state* of the engine without difficulty: serving the surfaces, running an agent over one patient's record in response to an orchestrator event, answering RAG queries, holding the cohort. The places it strains are three specific, bounded, batch-shaped jobs. For these we burst to rented GPU on RunPod — on-demand A100 80GB or H100 instances, spun up for the duration of the job and torn down after. This is not "the Spark can't do AI"; it is "three jobs are wall-clock-bound and embarrassingly parallel, and renting eight GPUs for two hours is cheaper and faster than waiting on one."

**1. GPU-accelerated variant calling (Parabricks).** Agent 1, the TSC-Variant Curator, starts from BAMs (synthetic BAMs for the demo, produced by BAMSurgeon — see §15) and runs deterministic calling: BWA-MEM alignment refinement, GATK HaplotypeCaller for germline, Mutect2 for low-VAF somatic mosaic recovery down to ≥5% VAF. NVIDIA Parabricks accelerates this end-to-end, but Parabricks expects datacenter-class GPUs and benefits enormously from running multiple samples in parallel. On the Spark, a single mosaic-aware calling run on one patient is tolerable for the demo's live Act-One moment (Patient A, the 8.3% VAF TSC2 frameshift). But generating the *initial* called VCFs for all 50 cohort patients, or for a real 200-500 patient institutional cohort, is exactly the kind of job you fan out: one GPU per sample (or per few samples), all at once, on RunPod. The output is deterministic and version-controlled, so this is a one-time-per-regeneration cost, not a per-demo cost.

**2. Parallel synthetic-cohort generation.** The four-layer cohort pipeline (§15) has two GPU-heavy layers: BAMSurgeon spike-ins producing the genomic substrate, and the frontier-model generation of ~600-1,000 clinical notes plus the imaging reports. The Synthea skeleton (layer 1) is CPU-bound and runs locally. The genomic substrate (layer 2) and the note/report generation (layers 3-4) parallelize cleanly across patients. The full cohort regeneration is budgeted at ~12 hours; the realistic way to hit that number is to fan the per-patient work across rented GPUs rather than serialize 50 patients on one device. Again: deterministic, version-controlled, generated once. After regeneration the cohort sits in MinIO on the Spark and no further RunPod time is needed to *use* it.

**3. Heavier local-LLM inference at batch scale.** Running Llama 3.1 70B for a single agent call on the Spark is fine. Running it across an entire cohort — for example, re-deriving the Phenome Mapper's HPO extraction over every note for all patients using the local model instead of the Sonnet API path — is throughput-bound on one GPU. When the institutional posture requires the local model (PHI not leaving the building) *and* the workload is a full-cohort batch, RunPod gives you the throughput. For the synthetic demo this path is rarely exercised, because the demo uses the Claude API; it is documented here because it is the bridge to the Phase-1 institutional build where the local-model batch path becomes real.

What does *not* need RunPod, stated explicitly so as not to overclaim in the other direction: serving any of the three surfaces; running any single agent over any single patient in response to a live orchestrator event; the entire 3-act demo, which operates on the already-generated cohort; ACMG synthesis and the Therapeutics brief, which are Claude API calls with no local-GPU dependency at all. The live demo touches RunPod for *nothing*. RunPod is a build-time and regeneration-time accelerator, not a runtime dependency of the demonstration. That distinction is the one to hold onto: a skeptical reviewer watching Act One should understand that the mosaic-recovery they are seeing live ran on the Spark in front of them, and that RunPod's only role was to produce the synthetic BAM corpus once, weeks earlier.

### 14.3 The local stack mirrors a the platform-backed deployment

The reason to belabor the storage and data-path detail is replication. As established from the CCHMC touchpoints (§18), the Winslow Research Pavilion is the *envelope*, and five CCHMC areas are *sources* that feed the engine — the Discover Together Biobank's banked tuber/AML/SEGA tissue, Biomedical Informatics' TAND methodology, the TSC clinical and research program, and the Epic Clarity/Caboodle plus biobank LIMS data plumbing. None of those institutional data paths are built in the synthetic demo; they are explicitly Phase-1 institutional work. But the *software* that consumes them is identical between the Spark demo and the institutional deployment. On the Spark, Milvus indexes a synthetic TSC corpus on local NVMe; in the Pavilion, the same Milvus configuration indexes the same corpus on the platform-backed shared storage, and the agents do not know or care which. On the Spark, MinIO holds watermarked synthetic BAMs; in the institutional build, the object store fronts real biobank-derived sequence under IRB. "Swap the box labels, keep the wiring." The Spark proves the wiring.

This is also the cleanest way to frame the move to a second site — TGen, City of Hope — without rebuilding. The engine's interface to the world is a small set of source contracts (a BAM source for the Variant Curator, a structured-EHR-plus-notes source for the Phenome Mapper and TAND agents, a meds/AE/literature source for the Therapeutics Strategist). Point those contracts at a new institution's sources and the agents and orchestrator are unchanged. A biobank without an intelligence layer is a freezer full of tubes; the intelligence layer is portable precisely because it was built to run on a $4,699 box first.

### 14.4 A cost model

The cost story has two regimes: the per-patient inference cost (dominated by Claude API usage, with a local-Llama floor) and the one-time generation/regeneration cost (dominated by RunPod GPU-hours). They scale differently and should be reasoned about separately.

**Per-patient inference cost.** A full pass of the engine over one patient touches the tiers roughly as follows: Haiku for high-volume normalization (ICD-10/lab→HPO mapping in the Phenome Mapper, prose summaries in the Trajectory Modeler), Sonnet for the per-note discourse and phenotype extraction work (Phenome Mapper and TAND, the highest token-volume tier because it reads every note), and Opus for the two reasoning-critical outputs (ACMG-AMP synthesis in the Variant Curator, the six-section Therapeutics brief). The Trajectory Modeler is classical statistics, not LLM inference, so its compute cost is negligible (CPU mixed-effects and Gaussian-process fits). A representative per-patient full pass, for a patient with on the order of 12-20 longitudinal notes:

| Tier | Primary use | Indicative tokens/patient | Indicative cost/patient |
|---|---|---|---|
| Haiku | HPO/ICD-10/lab normalization, prose summaries | ~150-300K | ~$0.05-0.10 |
| Sonnet | Per-note phenotype + TAND discourse extraction | ~400-700K | ~$1.50-2.80 |
| Opus | ACMG-AMP synthesis + 6-section Therapeutics brief | ~60-120K | ~$1.20-2.40 |
| **Total (Claude API path)** | | | **~$3-5 / patient / full pass** |

These are order-of-magnitude figures for planning, not a billed invoice; actual tokens depend on note length and how much RAG context the Therapeutics Strategist pulls. The number to internalize is single-digit dollars per patient for a *complete* multi-agent pass, and that a full pass is not the common case — the orchestrator's incremental-update discipline (§13) means that after the first pass, most events re-run only the affected agent over the changed slice, not the whole engine. Steady-state cost per patient-month is well below the first-pass number.

The local-Llama path changes the shape entirely: marginal API cost goes to zero, replaced by the amortized cost of GPU-hours (Spark electricity at runtime, or RunPod GPU-hours for batch). For a site running everything locally for PHI reasons, the per-patient *marginal* cost is effectively the power draw of the Spark, with the trade being lower reasoning quality on the Opus-tier tasks and slower throughput. The realistic institutional posture is hybrid: local model for the high-volume per-note extraction that touches the most PHI, Claude API (or a self-hosted frontier model) for the lower-volume, higher-stakes synthesis.

**One-time generation cost.** Regenerating the 50-patient synthetic cohort is the RunPod-heavy event. Budgeting conservatively — say 8 GPU-hours for the BAMSurgeon-plus-Parabricks genomic substrate fanned across patients, and a few more for the frontier-model note and imaging-report generation — a full regeneration lands in the low tens of dollars of RunPod time, on the order of a single working session, consistent with the ~12-hour wall-clock budget when parallelized. Because the cohort is deterministic and version-controlled, this cost is paid only when the cohort definition changes, not per demo and not per development day.

**Linear scaling, 50 → 1,000 → 10,000.** The per-patient inference cost is, to first order, linear in patients, because each patient is processed independently. Fifty patients at ~$3-5/full-pass is ~$150-250 of API cost for a complete cohort sweep; 1,000 patients is ~$3,000-5,000; 10,000 is ~$30,000-50,000 — for a *full* pass of all five agents over the entire cohort. Three things bend that curve downward in practice. First, incremental updates: real cohorts do not get full-passed daily; the orchestrator re-runs only changed slices, so the recurring cost tracks new notes and new events, not cohort size. Second, the hybrid local path moves the high-volume Sonnet-tier extraction onto amortized local GPU, which is where most of the per-patient token cost lives. Third, the generation cost does *not* scale linearly with patients in the institutional build, because real institutions provide real data — there is no BAMSurgeon spike-in or frontier-model note synthesis once the source contracts point at the biobank and Clarity. The synthetic-generation cost is a demo artifact that disappears at the real-data boundary.

The defensible summary for Hagedorn's audience: the entire demonstrable system runs on a single $4,699 device with no runtime cloud GPU dependency; rented GPU is used only to manufacture the synthetic cohort and could be replaced by real institutional data; per-patient reasoning cost is single-digit dollars and scales linearly, with a local-model path that drives marginal API cost toward zero at the expense of throughput and top-tier reasoning quality. The constraint of building on a small box is not a limitation we are apologizing for; it is the feature that makes the cost and replication story honest.

The 50-patient synthetic cohort that all of this generates and consumes is the subject of the next section, which specifies its composition, its four-layer pipeline, the three featured patients, and — equally important — what it deliberately does not contain.


## 15. The Synthetic Cohort

Everything described so far — five agents, a deterministic orchestrator, three clinician surfaces — needs patients to run against. It also needs patients we can show to Hagedorn's team, write about openly, version-control, and regenerate from scratch without a single approval signature. Real Cincinnati Children's data cannot do any of those things on the timeline of an eight-week build. So the demonstration runs on a cohort of fifty synthetic patients that we generate, own, and publish. This section explains why that is the right call rather than a compromise, exactly how the cohort is built, what it contains, and — equally important for credibility — what it does not.

The cohort is itself a deliverable. There is no public, multi-modal TSC patient corpus that pairs mosaic-aware genomics with longitudinal phenotyping, imaging trajectories, and clinical notes. The synthetic cohort we describe here is released under Apache 2.0 alongside the engine, with full generation provenance, so that anyone can reproduce the demo or build against it. That contribution is real even before a line of the institutional Phase-1 work exists.

### 15.1 Why synthetic

Three constraints push in the same direction.

**Practical.** The demo target is early Q3 2026, built first on a single DGX Spark with RunPod GPUs spun up as needed. A real-data path would require an executed data-use agreement, an IRB protocol, a CCHMC honest-broker de-identification pass, and a secure compute enclave — none of which fits an eight-week MVP and none of which is needed to prove the engineering. The build is the argument; the build needs data on day one.

**Ethical.** TSC is a rare disease (~1 in 6,000 births). A cohort that is rare-disease-specific, longitudinal, and multi-modal is, by construction, re-identifiable even after standard de-identification — small denominators, distinctive imaging trajectories, and idiosyncratic note language defeat safe-harbor stripping. Demonstrating a system on real children with a rare disease, before the system has institutional governance, is the wrong default. Synthetic data carries zero re-identification risk because there is no person behind the record.

**Legal.** A synthetic corpus has no PHI, needs no IRB for the demo, and can be committed to a public repository and screen-shared in a 30-minute pitch without a single agreement in place. It lets us hand the cohort to a skeptical informatics audience and invite them to inspect every record.

The trade is honest and explicit: synthetic data proves the *engineering* — that the agents run, that the orchestrator routes correctly, that ACMG synthesis matches combinatorial rules, that the TAND discourse markers fire on planted signals — but it is **not clinical validation**. Eval targets in §16 are measured against synthetic ground truth, not patient outcomes. We say this everywhere it could be misread.

### 15.2 Cohort composition

Fifty patients, with a genetic distribution chosen to mirror the real TSC population and, deliberately, to over-represent the cases the engine exists to handle. Roughly two-thirds of real TSC is *de novo* and one-third inherited; ~10–15% of clinically diagnosed patients are no-mutation-identified (NMI) on blood testing, usually because of somatic mosaicism. We seat enough mosaic and NMI patients to exercise the Variant Curator's low-VAF recovery path rather than leaving it as a single token case.

| Genotype class | Count | Share | What it exercises |
|---|---:|---:|---|
| TSC2 germline | 30 | 60% | Baseline variant calling + ACMG; most severe-end phenotypes |
| TSC1 germline | 12 | 24% | Generally milder phenotype; gene-aware interpretation |
| TSC2 mosaic | 5 | 10% | Low-VAF somatic recovery in tissue; mosaic flagging |
| TSC1 mosaic | 2 | 4% | Low-VAF recovery, second gene |
| NMI (mosaic, deep tissue) | 1 | 2% | The headline case: blood-negative, tuber-positive |
| **Total** | **50** | **100%** | — |

That is **seven mosaic patients** in total (5 + 2), plus the NMI patient whose only detectable variant lives in tuber tissue at low VAF — eight cases whose genomics the standard blood-only workflow would mishandle. The Variant Curator eval target is to recover all seven planted mosaic variants at VAF ≥ 5% with the correct ACMG class and no false-positive Pathogenic calls.

Phenotype severity is sampled to span the clinical range rather than cluster at the mean: SEGA present in a minority near the foramen of Monro, renal AMLs of varying size including some past the ~4 cm bleeding-risk threshold, epilepsy across the spectrum from well-controlled focal to refractory, and TAND features present in ~90% of the cohort consistent with the TOSCA registry — with a deliberate fraction of those features left *under-recognized in the notes* so the TAND Surveillance Agent has something real to surface.

### 15.3 The four-layer generation pipeline

The cohort is built in four layers, each producing a distinct modality, each version-controlled, and each carrying provenance into the next. The whole pipeline is deterministic given a fixed seed and a pinned set of model versions; a full regeneration takes roughly 12 hours, most of it BAMSurgeon and Parabricks-equivalent variant calling on RunPod GPUs.

```
Layer 1  Synthea + TSC modules        -> demographic + clinical skeleton (FHIR R4 + Clarity-shaped relational)
Layer 2  BAMSurgeon -> Parabricks-equiv -> realistic VCFs (germline + mosaic substrate)
Layer 3  Frontier-model clinical notes -> ~600-1000 longitudinal notes (watermarked)
Layer 4  Frontier-model imaging reports -> brain MRI / renal US / echo / ophtho (longitudinal series)
```

#### Layer 1 — Synthea skeleton

Synthea (MIT-licensed) generates the demographic and clinical backbone: birth, encounters, conditions, medications, lab orders, and an encounter timeline per patient. We author TSC-specific Synthea modules that lay down disease-appropriate care patterns — the surveillance cadence implied by the ITSC 2021 consensus guidelines (serial brain MRI for SEGA, renal imaging for AML, ophthalmology, cardiology in infancy), mTOR-inhibitor prescribing where indicated, and seizure-management encounters. Synthea emits both **FHIR R4** bundles and a **Clarity-shaped relational** export. The relational shape matters: it is the structure the Phenome Mapper expects to read in the synthetic demo, and it is the structure that a real Epic Clarity/Caboodle feed would later present — so the demo plumbing and the eventual institutional plumbing share a contract. To be unambiguous: we are building against *Clarity-shaped synthetic tables*, not against Clarity. The real Clarity/Caboodle and biobank-LIMS integration is institutional Phase-1 work and is explicitly not built here.

#### Layer 2 — BAMSurgeon genomic substrate

This is where the cohort earns the right to demonstrate mosaic recovery. We start from publicly available NA12878-derived alignments (we start at **BAM**, not raw FASTQ — see §15.6) and use **BAMSurgeon** to spike in TSC1 and TSC2 variants drawn from ClinVar and LOVD-TSC at controlled variant allele fractions:

- **Germline** variants at VAF ≈ 50% (heterozygous) for the 42 germline patients.
- **Mosaic** variants at **VAF 4–12%** for the seven mosaic patients, with the NMI patient's pathogenic frameshift placed near the low end so that blood-equivalent coverage would plausibly miss it while tissue-depth coverage recovers it.

The edited BAMs are then run through the Variant Curator's own deterministic calling path (Parabricks/BWA-MEM/GATK HaplotypeCaller + Mutect2, mosaic-aware down to VAF ≥ 5%) to produce realistic per-patient **VCFs** with authentic read support, strand balance, and the artifacts a curator would actually weigh. Crucially, the cohort generator records the planted variant — gene, transcript coordinate, intended VAF — as **ground truth** in a held-out manifest, so §16's eval can score recovery and classification against a known answer. We plant variants at biologically defensible VAFs and let the calling pipeline rediscover them; we do not hand the caller the answer.

#### Layer 3 — Frontier-model clinical notes

Roughly 600–1000 longitudinal clinical notes across the cohort — progress notes, neurology and nephrology follow-ups, genetics counseling notes, neuropsychology summaries — generated by a frontier model from published, de-identified TSC note templates and the patient's own Synthea timeline and genotype. The notes are the substrate for the Phenome Mapper (HPO extraction with evidence spans) and the TAND Surveillance Agent (discourse-marker analysis), so two things are engineered into them deliberately:

1. **Time-anchored clinical signal** consistent with each patient's genotype and imaging — seizures, developmental milestones, renal findings — that the Phenome Mapper must recover at recall ≥ 90% / precision ≥ 85%.
2. **Under-recognized TAND signals** planted as the diagnostic-uncertainty discourse markers from the Marshall-Hagedorn taxonomy — hedging, deferral, third-party attribution, conditional language, follow-up-without-formalization — scattered across notes so the TAND agent has genuine, non-obvious patterns to surface and the eval can check for both detection and the absence of spurious flags.

Every generated note is reviewed by clinician sampling for plausibility before it enters the cohort, and every note is **persistently watermarked as synthetic** (see §15.5). The notes are generated from templates and timelines, not free-associated, which keeps them clinically coherent and reproducible.

#### Layer 4 — Frontier-model imaging reports

The cohort contains **imaging reports**, not images. A frontier model generates radiology and ophthalmology report text — brain MRI (tubers, subependymal nodules, SEGA), renal ultrasound (AML size and count, cysts), echocardiography (rhabdomyomas in the youngest patients), and ophthalmology — anchored to each patient's genotype and severity sampling. The clinically load-bearing feature is **longitudinal series**: for SEGA-bearing patients we generate a sequence of MRI reports across visits with realistic interval growth (~2–4 mm/yr), because the Trajectory Modeler must forecast threshold crossing from that series and the eval requires Patient B's SEGA to cross within a 12–18-month window. The cohort contains no DICOM, no pixel data, and no real images (§15.6); imaging *AI* over real pixels is a separately described, not-built, institutional pipeline.

### 15.4 The three featured patients

Fifty patients give statistical texture; three carry the demo narrative. Each is constructed to make one capability undeniable on screen.

#### Patient A — the mosaic recovery (Act One)

4-year-old female, clinically diagnosed TSC, **blood testing negative** — a member of the ~10–15% NMI cohort. Tuber tissue (modeled as sourced from the Discover Together Biobank in the institutional framing) carries a **TSC2 frameshift variant at 8.3% VAF**. The Variant Curator's mosaic-aware path recovers it, the evidence aggregator assembles the case, and the ACMG-AMP synthesis returns **Likely Pathogenic on PVS1 + PM2 + PP4**, with a ddPCR orthogonal-validation recommendation and a mosaic flag carrying VAF, read support, strand balance, and artifact assessment. This is the only patient whose causative variant is invisible to standard blood workup, which is exactly why she opens the demo: Act One runs her live and then opens the audit trail. The output is a draft molecular-genetics report for board-certified geneticist sign-off, never an autonomous call.

#### Patient B — the longitudinal multi-agent picture (Act Two)

12-year-old male, **TSC2 c.3037C>T (p.Arg1013Ter)** — a clean nonsense variant, so the genomics is not the story here; the *integration over time* is. He carries a **SEGA growing 0.8 → 1.1 → 1.3 cm** across serial MRIs at the foramen of Monro (the hydrocephalus-risk location), **bilateral renal AMLs (~2.8 cm)**, well-controlled focal seizures, and **scattered under-recognized TAND signals** seeded into his notes. He is the four-quadrant in-visit dashboard made concrete: Variant Curator interpretation, Phenome Mapper HPO timeline, Trajectory Modeler forecast (his SEGA is engineered to cross the intervention-discussion threshold inside the 12–18-month window the eval checks, with calibrated 50%/90% prediction intervals and no false alarms on the rest of the cohort), and TAND briefing material. He shows what it means to have five agents looking at one child over years.

#### Patient C — the therapeutics strategist (Act Two)

18-year-old female, **TSC1**, with a real therapeutic management problem: a **partial response to everolimus** complicated by **mucositis that forced a dose reduction**, an **AML approaching ~4 cm** (nearing the bleeding-risk threshold), and **refractory focal seizures**. She is the case that exercises the Opus-class Therapeutics Strategist's six-section options brief — Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions — every claim source-attributed to PubMed/PMC, a ClinicalTrials.gov snapshot, or FDA actions, framed as decision-support with appropriate hedging, not a recommendation. Her complexity (efficacy/toxicity trade-off, threshold-adjacent AML, refractory epilepsy) is what makes a structured options brief worth reading.

The three together cover the engine's span: recovery of the invisible (A), longitudinal multi-agent integration (B), and therapeutic decision-support under genuine ambiguity (C).

### 15.5 Fidelity, plausibility, sufficiency, and watermarking

Four properties define whether the cohort is fit for purpose, and we hold ourselves to each explicitly.

**Structural fidelity.** The data lands in the shapes the engine and a real institution would use: FHIR R4 and Clarity-shaped relational tables (Layer 1), genuinely called VCFs with real read evidence rather than hand-asserted variant lists (Layer 2), HPO and SNOMED-CT-codeable note content (Layer 3). The Phenome Mapper does not get an easier parsing job on synthetic data than it would on real data — that is the point of conforming to the same contract.

**Clinical plausibility.** Genotype, phenotype, imaging trajectory, and treatment history are internally consistent per patient (a TSC2 nonsense patient with a foramen-of-Monro SEGA gets the surveillance MRIs and the neurology follow-ups that implies), severity is sampled across the real range rather than to a convenient mean, and clinician sampling reviews the notes and imaging reports before they enter the cohort. Plausibility is the credibility gate with Hagedorn's audience; an implausible note is a tell, and a tell is fatal in this room.

**Demonstration sufficiency.** The cohort is sized and shaped to exercise every agent and every eval target — seven mosaic variants for the Variant Curator, ≥90%/≥85% phenotype extraction for the Phenome Mapper, a forecastable SEGA trajectory for the Trajectory Modeler, planted-but-not-obvious TAND signals for the Surveillance Agent, and a genuinely ambiguous therapeutic case for the Strategist. Fifty is enough to make the demo and the eval meaningful while staying inside a ~12-hour regeneration budget.

**Honest watermarking.** Every record is persistently labeled synthetic. Notes and imaging reports carry an inline synthetic watermark; the surfaces render against an unmistakable synthetic-data banner; AI-generated content is labeled as such (the Variant Curator's report is "AI-labeled draft"). The surfaces are standalone web apps, **not** Epic, by design — there is no path by which a synthetic record could be mistaken for or written back to a real chart. We would rather over-label than ever let a viewer forget what they are looking at.

### 15.6 Reproducibility and what the cohort does not contain

**Reproducibility.** The cohort is generated **once**, **deterministically** (fixed seeds, pinned model versions and prompt-template versions), and **version-controlled** end to end — the Synthea modules, the BAMSurgeon variant manifest with intended VAFs, the note and imaging-report generation prompts, the ground-truth answer key, and the resulting artifacts. A full regeneration is ~12 hours on the Spark + RunPod. Anyone with the repository can reproduce the exact cohort and therefore reproduce the eval. This is what makes the demo an experiment rather than a screenshot.

**What it deliberately does not contain.** Stating the boundaries plainly is part of not overclaiming:

- **No real imaging or DICOM** — imaging is report *text* only; no pixel data, no real scans.
- **No raw FASTQ** — the genomic pipeline starts at BAM; we do not synthesize or distribute raw reads.
- **No neuropsychological test scores** — TAND signal comes from note language, not from instrument data.
- **No pedigree** beyond what a note mentions — no structured family-history graph.
- **No pharmacy or claims data** — medications and adherence come from the Synthea/notes layer, not a pharmacy feed.
- **No patient-reported-outcome (PRO) scores.**

And, at the system level: there is **no Epic/Clarity/Caboodle integration, no biobank LIMS integration, and no imaging-AI pipeline over real pixels** in this build. Those are described architecturally elsewhere in the paper and are institutional Phase-1 work. The synthetic-data demo is what runs now on the Spark and RunPod; everything bank-, chart-, and image-connected is future, governed work.

With the patients in hand, §16 turns to the machinery that reasons over them: the tiered model assignments, the TSC RAG corpus, the audit and provenance layer, and the eval harness that scores every agent against the ground truth this cohort was built to carry.


## 16. Model Tiering, the TSC RAG Corpus, Audit & Evaluation

The preceding sections described what each agent does and how the orchestrator wires them together. This section covers four cross-cutting concerns that determine whether the engine is trustworthy rather than merely functional: which model runs each step and why, what the agents read when they reason, what record every invocation leaves behind, and how we know the outputs are right against a ground truth we control. These are the parts a skeptical informatics reader inspects first. A demo that produces a plausible briefing is easy. A demo whose every claim can be traced to a model version, a prompt template, and a retrieved source, and whose accuracy is measured against a held-out answer key, is the thing worth building.

### 16.1 Model tiering policy

The engine does not use one model. It routes each step to the cheapest tier that can do it reliably, and pushes the work that is irreversible or clinically consequential up to the most capable tier. The principle is blunt: deterministic work never touches an LLM, normalization and extraction run on fast models, synthesis and conflict resolution run on Opus, and the steps that produce a clinician-facing claim are Opus by policy and not negotiable. Where a frontier model is genuinely required and a network is unavailable, the local Llama 3.1 70B Instruct served via Ollama is the documented fallback, with the substitution recorded in the audit log (see §16.3) because a fallback is a material change to the reasoning substrate, not an implementation detail.

The tiers, reused from the HCLS AI Factory v1.3.0 substrate, are Claude Haiku (fast, cheap normalization), Claude Sonnet (per-note extraction, evidence aggregation, unusual-case interpretation), Claude Opus (synthesis, conflict resolution, clinician-facing briefs), and a non-LLM tier for everything deterministic. The table below is the controlling policy. Each agent's own section gives the clinical rationale; here we state the routing rule and the reason it sits where it does.

| Agent | Step | Tier | Rationale |
|---|---|---|---|
| Variant Curator | BAM→VCF calling (BWA-MEM/GATK HaplotypeCaller + Mutect2, mosaic-aware) | Deterministic (Parabricks) | Variant calling is reproducible computation; an LLM here would be malpractice. |
| Variant Curator | snpEff/VEP annotation | Deterministic | Coordinate-based annotation; fixed tool versions, hashed. |
| Variant Curator | Evidence aggregation across ClinVar/gnomAD/LOVD-TSC | Sonnet | Reconciles structured evidence into a coherent narrative; not the final call. |
| Variant Curator | ACMG-AMP classification synthesis | Opus, validated vs combinatorial rules | The classification is the consequential output; Opus reasons over criteria, then a deterministic rule engine checks that the criteria combine to the asserted class. Disagreement is flagged, not hidden. |
| Phenome Mapper | ICD-10 / LOINC lab → HPO normalization | Haiku | High-volume, low-ambiguity mapping against fixed crosswalks. |
| Phenome Mapper | Per-note phenotype extraction with evidence spans | Sonnet | Discourse-level reading; the bulk of the cohort runs here for cost reasons. |
| Phenome Mapper | Rare cross-note conflict resolution | Opus | Reserved for the small fraction of genuinely discordant cases. |
| Trajectory Modeler | Mixed-effects / Gaussian process / survival / Bayesian forecasting | Non-LLM (classical statistics) | Forecasts with calibrated 50%/90% intervals must come from statistical models, not a language model's prose. |
| Trajectory Modeler | Plain-language summary of the forecast | Haiku | Renders numbers the statistical layer already produced; adds no inference. |
| Trajectory Modeler | Interpretation of unusual trajectories | Sonnet | Flagging, not forecasting; the numbers remain the statistical model's. |
| TAND Surveillance | Per-note discourse analysis (Marshall-Hagedorn markers across 6 TAND clusters) | Sonnet | Pattern detection over hedging/deferral/attribution language; the methodological core. |
| TAND Surveillance | Scoring / aggregation across notes | Deterministic | Counts and thresholds are reproducible; the model finds signals, the rules score them. |
| TAND Surveillance | Pre-visit briefing summary | Opus | Clinician-facing prose; must hedge correctly and never imply diagnosis. |
| Therapeutics Strategist | Six-section options brief (with PubMed/PMC RAG, ClinicalTrials.gov, FDA) | Opus, non-negotiable | The highest-stakes output in the engine; every claim source-attributed, framed as decision-support. |
| Orchestrator | 13 event-type routing, state projection | Deterministic (LangGraph) | The router is plumbing; an LLM in the control path would make the system non-reproducible. |

Two routing decisions deserve emphasis because they are where most systems quietly cut corners. First, the Trajectory Modeler is classical statistics end to end; the only LLM touches are a Haiku prose wrapper over numbers the statistical layer already computed and a Sonnet interpreter for outliers. We do not let a language model forecast a SEGA growth curve. A prediction interval that cannot be derived from a fitted model is not a prediction interval. Second, ACMG-AMP classification is Opus reasoning followed by a deterministic check: Opus assigns criteria (PVS1, PM2, PP4, and so on) and proposes a class, and a combinatorial rule engine independently verifies that those criteria, under the published ACMG-AMP combining rules, actually yield the asserted class. When Opus and the rule engine disagree, the report carries the disagreement to the molecular geneticist rather than silently picking one. The model proposes; the rules audit; the human signs.

### 16.2 The TSC RAG corpus

Three agents read literature: the Variant Curator (for variant-level evidence framing), the TAND Surveillance Agent (for the methodological grounding of its markers), and above all the Therapeutics Strategist, whose entire output is a synthesis of evidence that must be attributable. They read from a single curated corpus held as a dedicated partition in the existing Milvus instance, kept separate from the platform's other engines so that a TSC query never retrieves an oncology passage and so the partition can be re-embedded and versioned independently.

The corpus is deliberately small and curated rather than a scrape of the literature. Precision matters more than recall for this audience: a Therapeutics brief that cites the ITSC consensus and EXIST-3 correctly is worth more than one that surfaces forty marginally relevant abstracts. The seed corpus is organized into named source classes, each with an ingestion provenance record (source, retrieval date, license posture, document hash):

| Source class | Contents | Primary consumer |
|---|---|---|
| Surveillance guidelines | ITSC 2021 consensus surveillance and management recommendations; TAND consensus framework and the TAND-L lifetime checklist | Phenome Mapper (surveillance-gap framing), TAND Agent, Therapeutics |
| Somatic mosaicism literature | Tyburczy 2015, Giannikou 2016, Lim 2017, and the deep-sequencing methods they describe | Variant Curator |
| Targeted-therapy trials & registries | EXIST-3 (everolimus in TSC-associated refractory epilepsy) and the broader EXIST program; TOSCA registry findings on TAND under-recognition | Therapeutics, TAND Agent |
| Diagnostic-uncertainty NLP methodology | Marshall/Nickels/Brady/Hagedorn 2023 (J Hospital Medicine); Nickels et al. 2024 (Applied Linguistics); Ipsaro et al. 2021 (Hospital Pediatrics); Orenstein et al. 2021 (JAMIA, alert burden) | TAND Agent |
| Trial landscape | ClinicalTrials.gov snapshot filtered to TSC, mTOR inhibitors, and TSC-associated epilepsy (point-in-time, dated) | Therapeutics |
| Regulatory actions | FDA labeling and action snapshots for everolimus/sirolimus and TSC indications (point-in-time, dated) | Therapeutics |

Two facts about this corpus are non-negotiable for credibility. First, the ClinicalTrials.gov and FDA sources are explicitly snapshots, frozen at a stated retrieval date and stamped as such in every passage's metadata. A Therapeutics brief never implies it is reading a live registry. It reads a dated snapshot, and the date is shown to the clinician. Second, the inclusion of the Hagedorn group's own publications is not decoration. The TAND Surveillance Agent's discourse-marker taxonomy is an extension of that published work, and the agent retrieves and cites those papers so that the methodological lineage is visible in the output, not asserted in a slide.

#### Embeddings and retrieval

The corpus is embedded with the substrate's two-model scheme: BAAI/bge-large-en-v1.5 for general biomedical prose and a BiomedBERT-derived clinical embedding for the clinically dense passages (guideline criteria, trial endpoints). Documents are chunked at section granularity with overlap, and each chunk carries its source class, document hash, page or section anchor, and retrieval date as Milvus payload fields. Retrieval is filtered by source class at query time. The Therapeutics Strategist, for example, retrieves trial and regulatory snapshots with their dates attached so that the citation it emits carries the same date the clinician sees. Every retrieved chunk's URI and hash flow into the audit record (§16.3); a citation in a brief is therefore not free text but a pointer back to a specific, hashed passage in a specific corpus version.

The corpus is version-controlled as a manifest (source list, retrieval dates, hashes, embedding model versions) checked into the repository alongside the synthetic cohort. Re-embedding produces a new corpus version id, and that id is recorded on every invocation that used it. This is the same discipline applied to the cohort: the inputs to reasoning are versioned artifacts, not ambient state.

### 16.3 Audit and provenance

This is the part Hagedorn's team checks first, and it is the right thing to check first. The Orenstein et al. 2021 JAMIA work on alert burden and the broader diagnostic-uncertainty research that grounds the TAND agent share a premise: a clinical informatics tool earns trust by being inspectable, and loses it the moment a clinician cannot answer "why did it say that?" So every agent invocation in the engine writes an append-only provenance record before its output is allowed to surface. There is no path by which a claim reaches a clinician surface without a corresponding audit row. This is enforced structurally in the event-sourced state layer (the PostgreSQL append-only event log described in the orchestrator section), not by convention.

The per-invocation record captures the full reasoning context. A representative schema:

```json
{
  "invocation_id": "uuid",
  "event_id": "uuid",
  "patient_id": "synthetic-A",
  "agent": "therapeutics_strategist",
  "step": "options_brief_v1",
  "timestamp": "2026-07-08T14:22:07Z",
  "model": {
    "id": "claude-opus-4-8",
    "tier": "high",
    "fallback_used": false,
    "fallback_model": null
  },
  "prompt_template": {
    "name": "therapeutics_six_section",
    "version": "v3"
  },
  "inputs": {
    "input_hash": "sha256:…",
    "upstream_invocations": ["uuid-variant", "uuid-phenome", "uuid-trajectory", "uuid-tand"],
    "cohort_version": "cohort-2026.07.01",
    "rag_corpus_version": "tsc-corpus-2026.07.05"
  },
  "retrieved_sources": [
    {"uri": "pmc:EXIST-3", "source_class": "trials", "chunk_hash": "sha256:…", "retrieval_date": "2026-07-05", "score": 0.83},
    {"uri": "ctgov:NCT…", "source_class": "trial_landscape", "snapshot_date": "2026-07-01", "chunk_hash": "sha256:…", "score": 0.79}
  ],
  "output_hash": "sha256:…",
  "latency_ms": 41280,
  "validation": {
    "acmg_rule_check": null,
    "human_gate": "pending_molecular_geneticist"
  },
  "watermark": "SYNTHETIC"
}
```

The fields are not arbitrary. Each answers a question a reviewer will actually ask:

- **model.id / tier / fallback_used.** Which model produced this, at what tier, and was the local Llama fallback substituted? A fallback is flagged because it changes the reasoning substrate; a reviewer must be able to discount or re-run any output that did not run on the intended tier.
- **prompt_template.name / version.** The exact instruction the model received. Prompt templates are versioned in the repository, so a record from July can be replayed against the prompt that actually generated it, not whatever the prompt has since become.
- **input_hash / upstream_invocations / cohort_version / rag_corpus_version.** The complete input provenance. Because the Therapeutics brief integrates all four upstream agents, its record names the specific upstream invocations it consumed, so a reviewer can walk the chain from a therapeutic claim back to the variant call and the HPO profile that informed it. The cohort and corpus version ids make the whole thing reproducible.
- **retrieved_sources.** Every RAG chunk with its URI, source class, hash, retrieval or snapshot date, and similarity score. This is what turns a citation in the brief into a verifiable pointer. A clinician who doubts a trial-matching claim clicks through to the dated snapshot the model actually read.
- **output_hash / latency_ms.** Tamper-evidence and the performance figure that the eval harness checks against the per-agent latency targets.
- **validation.** For the Variant Curator, the result of the deterministic ACMG rule check against Opus's proposed class. For every clinician-facing output, the human-gate state, which begins as pending and is never auto-cleared. The record exists to show that a human signed, or that no human has yet.

The record is append-only and queryable. "Show me every Pathogenic call this week and the criteria and reviewer for each" is a SQL query, not a forensic exercise. That property, more than any single model choice, is what makes the engine defensible to the audience it is built for.

### 16.4 The evaluation harness

The synthetic cohort exists so that we can grade the engine against an answer key we wrote. Every patient was generated with known ground truth: the inserted TSC1/TSC2 variants and their VAFs, the seeded HPO phenotypes, the parameterized growth trajectories, and the deliberately embedded TAND signals. The harness runs each agent over the cohort and scores its output against that ground truth, per agent, with the targets fixed in advance.

| Agent | Eval target (vs synthetic ground truth) | What a pass means |
|---|---|---|
| Variant Curator | Recover all 7 mosaic variants at VAF ≥ 5%; correct ACMG class for each; zero false-positive Pathogenic calls; < 5 min/case | The mosaic-recovery claim (the engine's headline) holds, and it does not manufacture pathogenicity. |
| Phenome Mapper | Recall ≥ 90%, precision ≥ 85% on seeded HPO terms; full cohort < 1 hr | The phenotype foundation other agents build on is accurate and tractable at cohort scale. |
| Trajectory Modeler | Forecast Patient B's SEGA crossing the action threshold within the 12–18 month window; no false alarms across stable patients | The forecasting layer's intervals are calibrated enough to be actionable without crying wolf. |
| TAND Surveillance | Detect the embedded TAND signals; zero spurious flags on patients with none | The Marshall-Hagedorn extension surfaces real signal without inventing it. |
| Therapeutics Strategist | Correct trial matches; appropriate hedging; full source attribution; < 3 min | The six-section brief is accurate, properly cautious, and every claim is traceable. |

The mosaic targets are the most important and the most exposed. The seven mosaic patients (five TSC2 mosaic, two TSC1 mosaic, plus the one NMI case whose tuber-tissue variant sits at 8.3% VAF) were built with BAMSurgeon at VAFs between 4% and 12% precisely so that the ≥ 5% recovery target is a real test and not a foregone conclusion. Recovering them while making zero false-positive Pathogenic calls across the germline majority is the claim the whole demo rests on, and it is the first number we report.

#### What it proves and what it does not

This is a synthetic-data evaluation, and the harness is honest about its ceiling. What it proves: the pipeline is correct against a known answer key, the model tiering produces the right classifications and forecasts, the retrieval and attribution machinery works end to end, and the latency budgets hold on the Spark-plus-RunPod build. It demonstrates that the engineering does what it claims, on data whose truth we defined.

What it does not prove, and we say so plainly: this is not clinical validation. The synthetic notes were written from published templates and sampled by a clinician for plausibility, but they are not real charts; the variants were inserted, not sequenced from patients; the trajectories were parameterized, not observed. Recall and precision against seeded HPO terms are not the same as recall against the messy reality of a CCHMC clinic note, and a SEGA forecast that lands on a synthetic curve says nothing yet about a real child. Establishing clinical performance requires the institutional Phase-1 work described elsewhere in this paper: real data behind the Epic Clarity and biobank LIMS integrations that are explicitly not built in this demo, under IRB, with prospective evaluation. The synthetic harness is the necessary evidence that the engine is built correctly. It is not, and is never presented as, evidence that it is clinically validated.

With the model routing, the corpus, the audit substrate, and the evaluation targets fixed, the remaining question is the order in which to build it. The next section sequences the eight-week build so that the highest-risk claim, mosaic recovery on Patient A, is provable as early as possible.


## 17. Build Sequencing & Engineering Discipline

The preceding sections described what the TSC Intelligence Engine does and why each design choice holds. This section describes how it gets built in eight weeks on a single DGX Spark with RunPod GPUs on demand, and the engineering discipline that keeps the demo honest. The principle throughout is simple: the build is the argument. A working system that does five things and admits the boundary of the sixth is more persuasive to a skeptical informatics audience than a slide deck that promises six. What follows is the summary arc; the companion Product Requirements Document carries the task-level detail (acceptance tests, schemas, file manifests, and per-agent prompt-template versions).

### The eight-week arc

The schedule is sequenced by data dependency, not by agent number. The synthetic cohort is the substrate everything else consumes, so it leads. Among the agents, the TSC-Phenome Mapper is built early because the Trajectory Modeler, TAND Surveillance Agent, and Therapeutics Strategist all read its longitudinal HPO profile. The TSC-Variant Curator is independent of the phenome work and runs as a parallel track from Week 2, gated only by the genomic substrate (BAMSurgeon-edited BAMs) being ready. The TSC-Therapeutics Strategist is last because it integrates the outputs of all four prior agents.

| Week | Primary deliverable | Featured patient / eval anchor |
|------|---------------------|-------------------------------|
| W1 | Cohort foundation: Synthea + TSC modules wired; BAMSurgeon pipeline stood up; orchestrator event schema + PostgreSQL append-only event log defined | Synthetic harness reproducible; deterministic seed locked |
| W2 | Full 50-patient cohort generated (genomic + structured + draft notes); Variant Curator core (BWA-MEM/GATK HaplotypeCaller+Mutect2, mosaic-aware low-VAF path) | All 7 mosaic variants present in BAMs at VAF 4–12% |
| W3 | Variant Curator ACMG-AMP synthesis (Opus, validated against combinatorial rules); Phenome Mapper extraction (Sonnet per-note, Haiku ICD-10/lab→HPO) | Patient A NMI 8.3% VAF → Likely Pathogenic (PVS1+PM2+PP4) |
| W4 | Phenome Mapper hardening (recall/precision targets); Trajectory Modeler (mixed-effects, GP regression, survival, Bayesian updating) | Patient B SEGA 0.8→1.1→1.3 cm forecast to cross threshold in 12–18 mo |
| W5 | TAND Surveillance Agent (Marshall-Hagedorn discourse taxonomy, deterministic scoring layer); Therapeutics Strategist scaffolding + RAG corpus | Patient B embedded TAND signals detected, no spurious flags |
| W6 | Therapeutics Strategist six-section brief complete (PubMed/PMC RAG + ClinicalTrials.gov + FDA); pre-visit briefing + in-visit dashboard surfaces | Patient C TSC1, partial everolimus response, AML ~4 cm, trial matching |
| W7 | Surface completion (async alert surface, progressive disclosure, source navigation); orchestrator integration (13 event types, dependency-ordered enrollment) | End-to-end: enroll → all five agents → three surfaces |
| W8 | Cohort regeneration from clean seed (~12 hr); dry runs; clinician spot-review of synthetic notes; provenance audit; delivery | Full 3-act demo rehearsed; eval targets confirmed |

The eval anchors in the right column are the contract for each week. They are measured against the cohort's own synthetic ground truth, which is recorded at generation time, not inferred after the fact. This is demo-scale validation, not clinical validation, and the paper says so plainly wherever the claim appears.

### Parallel work streams

Eight weeks is not enough for a strictly serial build, so four streams run concurrently after Week 1, with explicit hand-off points where one stream's output becomes another's input.

- **Cohort stream (W1–W2, regen W8).** Owns the 4-layer synthetic pipeline: Synthea+TSC modules → BAMSurgeon genomic substrate → Parabricks-equivalent calling → frontier-model clinical notes and imaging reports. Produces the ground-truth ledger. Once the cohort is frozen after W2, this stream goes quiet until the W8 clean regeneration, which exists to prove the cohort is deterministically reproducible from seed (a reviewer can rebuild it).
- **Agent stream (W2–W6).** Two sub-tracks. Track A is the Variant Curator, dependent only on the BAMs, so it can run in parallel with phenome work from Week 2. Track B is the dependency chain: Phenome Mapper (W3–W4) → Trajectory Modeler (W4) and TAND Agent (W5) in parallel → Therapeutics Strategist (W5–W6), which closes the chain by integrating all four.
- **Orchestrator stream (W1, integration W7).** The deterministic LangGraph event router and the event-sourced state schema are defined in Week 1 so the agents emit and consume events from the start, against stubbed peers where needed. Full integration — dependency-ordered enrollment, incremental-update minimization, demand-driven surface assembly, conservative failure handling — lands in Week 7.
- **Surface stream (W6–W7).** The three clinician surfaces (pre-visit briefing, in-visit 4-quadrant dashboard, async alert surface) are built late because they are demand-driven assemblies over agent outputs; building them before the agents produce real structured outputs would mean building against fiction. They are standalone web apps, persistently watermarked synthetic, not Epic-embedded.

The hand-off discipline matters more than the parallelism. Each stream publishes a versioned output contract (the event schema, the HPO profile shape, the agent output envelope) so a downstream stream can develop against the contract before the upstream implementation is finished. When the contract changes, it changes in one place and the event log absorbs it.

### De-scope priority order

The schedule will slip somewhere; planning for it is part of the discipline. The de-scope order is decided now, in calm, rather than under Week 7 pressure, and it is ordered to protect the two things the demo cannot lose: the live mosaic recovery (Act One) and a coherent in-visit picture (Act Two).

1. **Cohort 50 → 30 patients.** The first lever. Thirty patients still contains all 7 mosaic cases and all three featured patients (A, B, C); the cohort exists to demonstrate behavior, not to power a statistic. Halving generation and review time buys back the most schedule for the least demonstrative loss.
2. **Simplify the async alert surface.** The pre-visit briefing and in-visit dashboard carry the demo; the alert surface is the most defensible to thin (fewer categories, coarser recalibration logic) without weakening the core story.
3. **Simplify the TAND cluster set.** Reduce from the full six TAND clusters (behavioral, psychiatric, intellectual, academic, neuropsychological, psychosocial) to the two or three that carry Patient B's embedded signals. The Marshall-Hagedorn discourse method is still demonstrated; the breadth is trimmed.
4. **Last resort: Therapeutics six → four sections.** The Therapeutics Strategist's six-section brief (Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions) drops Combination and Emerging Evidence if nothing else has recovered enough time. This is last because the Strategist is the integrative payoff of the whole engine, and trial matching plus source-attributed hedging is the part Hagedorn's audience will scrutinize.

What is never on the de-scope list: provenance, the human-review gates, the synthetic watermark, and the Variant Curator's mosaic recovery path. Those are load-bearing for credibility, not features.

### Engineering-discipline themes

Four themes run through every stream and every week. They are the difference between a demo that survives expert scrutiny and one that does not.

**Honest scope.** The synthetic-data demo is what runs now on Spark plus RunPod. The Epic Clarity/Caboodle and biobank LIMS plumbing, and the imaging-AI pipelines, are described architecturally and explicitly not built — they are CCHMC Phase-1 institutional work, gated on IRB and real data. Every artifact that touches that boundary marks it. The five CCHMC source areas (Discover Together Biobank, Biomedical Informatics, the TSC clinical/research program, and the data plumbing) feed the engine; they are sources, not components, and the demo simulates their outputs rather than claiming integration. A biobank without an intelligence layer is a freezer full of tubes; the demo shows the intelligence layer over simulated freezer contents, and says exactly that.

**Watermark everything.** Every synthetic note, imaging report, dashboard, and exported brief is persistently and visibly watermarked as synthetic. This is not decoration. With a clinical audience, an unwatermarked synthetic note is a credibility hazard the moment it is screenshotted out of context. The watermark is applied at generation for notes and at render for surfaces, so it cannot be lost in transit.

**Human-review gates.** No agent output is autonomous or terminal. The Variant Curator emits an AI-labeled draft molecular-genetics report for sign-off by a board-certified molecular geneticist, with the mosaic flag and ddPCR validation recommendation attached. The TAND Agent surfaces patterns as pre-visit briefing material — never an interruptive alert, never a diagnosis. Every surface carries the decision-support framing in the interface itself, not just in documentation. The gate is in the product, not in a disclaimer.

**Conservative alerts and conservative failure.** Alert discipline is a measured target: if the async surface produces more than roughly three alerts per clinician per week, the thresholds recalibrate. Alert fatigue is a documented failure mode in this exact institution's prior informatics work (Orenstein et al. 2021, JAMIA), so under-alerting is the safer error. The same conservatism governs failure: when an agent fails, the orchestrator's affected surface shows "pending" or a staleness marker — never a silent missing output that a clinician might read as a true negative. Provenance backs all of it: every output records model id and version, prompt-template version, retrieved RAG sources with URIs, input hash, and latency, in an append-only, queryable log.

Taken together, these themes are why the de-scope order protects what it protects. Provenance, gates, watermarks, and mosaic recovery are not the parts of the system that look impressive in a thirty-minute demo; they are the parts that let the demo claim anything at all. With the build sequence and its discipline established, the next section turns from how the engine is built to where it lands first: the Cincinnati Children's beachhead.


---

# Part V — The Cincinnati Children's Case


## 18. Cincinnati Children's: The Beachhead, the Winslow Pavilion, and the Touchpoints

Everything to this point in the paper has been about the engine in the abstract: five agents, a deterministic orchestrator, three clinician surfaces, a synthetic cohort that lets all of it run today on a single workstation. None of that is institution-specific. The architecture would slot into any pediatric academic medical center with a TSC program and a biobank. But an engine that runs anywhere runs nowhere in particular, and "runs nowhere in particular" is not a deployment. This section names the particular place: Cincinnati Children's Hospital Medical Center (CCHMC), and the specific physical, data, and human infrastructure there that the engine connects to. The case for CCHMC as the first deployment — the beachhead — rests on a coincidence of three things that rarely co-locate: a clinical patient concentration deep enough to be worth modeling, an institutional research apparatus capable of operating real biospecimens and real records, and a leadership figure whose own published methodology is one of the engine's five agents. CCHMC has all three.

### 18.1 Why CCHMC is the right beachhead

The first criterion for a beachhead is that the disease is actually present in volume. TSC affects roughly 1 in 6,000 births, which makes it rare in the population but not rare in a referral center that has organized itself around it. CCHMC runs a multidisciplinary TSC clinical and research program: neurology, nephrology, genetics, neurosurgery, dermatology, and neuropsychology coordinated around a single patient population, with the longitudinal follow-up that TSC demands. SEGA growth is measured in millimeters per year; AML bleeding risk turns on whether a lesion has crossed roughly 4 cm; TAND features accumulate across a lifetime and go unaddressed in 30 to 50 percent of patients. A platform whose entire premise is longitudinal integration — surveillance cadence, trajectory forecasting, time-anchored phenotyping — needs patients who have been followed for years across multiple subspecialties. A center without that depth would force the engine to model trajectories from two data points. CCHMC's program supplies the depth.

The second criterion is research capacity: not just patients, but the institutional machinery to turn patients into governed, analyzable substrate. CCHMC is among the largest pediatric academic medical centers in the United States, with roughly 16,000 staff, more than 600 faculty, and over $300M per year in research expenditure. That scale matters in a specific, unglamorous way. The institutional Phase-1 work this engine eventually requires — IRB protocols, a data-use governance path, an Epic Clarity/Caboodle extract, a biobank LIMS query, a SaMD determination — is the kind of work that only exists at institutions that already do it routinely. A small center might have the patients and the enthusiasm and still be unable to stand up a governed real-data pipeline in any reasonable timeframe. CCHMC has done versions of this many times.

The third criterion is leadership alignment, and here CCHMC is not merely adequate but unusually well-suited, because the methodology underneath one of the five agents was authored there. That is the subject of §18.5. Taken together, the three criteria converge: enough patients to model, enough institution to govern real data, and a sponsor whose own research the engine extends rather than imports. The honest version of the beachhead argument is that we did not choose TSC and then go looking for a hospital. We chose the hospital where the disease, the data apparatus, and the methodology already coincided, and let that coincidence define the disease.

### 18.2 The Winslow Research Pavilion: the infrastructure envelope

The physical center of gravity for this work is the Winslow Research Pavilion, a roughly 45,000-square-foot research building that opened in 2024 in the Avondale neighborhood adjacent to the main CCHMC campus. The Pavilion is best understood not as a single lab but as an *envelope*: a building that houses several distinct research functions under one institutional roof and one set of operating norms. It contains the Gamble Vaccine Research Center, the Discover Together Biobank, the institution's first centralized biospecimen freezer archive, and roughly 70 research staff. The CCHMC Research Foundation, directed by Dr. Tina Cheng, is the institutional home for the research enterprise the Pavilion concentrates.

The distinction between "envelope" and "engine" is load-bearing for the rest of this section, so it is worth stating precisely. The Winslow Pavilion is the envelope: a physical and institutional container. The TSC Intelligence Engine — the five agents and the orchestrator — is not inside the Pavilion in any architectural sense. The Pavilion supplies *sources*: tissue in freezers, a biobank with a mission, an informatics division, a clinical program, data plumbing. The engine consumes those sources and produces clinician surfaces. The relationship is supplier-to-consumer, not part-to-whole. This matters because it is the thing that makes the engine portable. If the Pavilion were the engine, the engine could not be replicated; because the Pavilion is an envelope of sources, the engine can be lifted to TGen or City of Hope by swapping which sources fill which slots. We return to that "swap the box labels, keep the wiring" property at the end of the section.

The most consequential single asset inside the envelope is the freezer archive — specifically, the banked tissue that becomes the molecular substrate for the TSC-Variant Curator. Which brings us to the biobank.

### 18.3 The Discover Together Biobank: tissue as substrate, and the freezer-without-a-layer problem

The Discover Together Biobank is CCHMC's institutional biospecimen resource, housed in the Pavilion and built around three stated missions:

1. **Find biomarkers** — discover measurable signals that predict, stratify, or track disease.
2. **Identify disease-causing DNA changes** — connect banked specimens to the genetic variants that explain phenotype.
3. **Understand resilience** — learn why some individuals with a given genetic burden do better than others.

Each of those missions maps cleanly onto something the engine does. Biomarker discovery is, in part, what the Trajectory Modeler operationalizes when it learns which early measurements predict later threshold crossings. Resilience is the kind of question the longitudinal HPO profiles from the Phenome Mapper and the trajectory forecasts are built to support. But the mission that matters most for the demonstration — the one that closes the loop on the engine's most distinctive claim — is the second: identifying disease-causing DNA changes. That is the TSC-Variant Curator's entire job, and the biobank is where its hardest cases come from.

Those hardest cases are the no-mutation-identified (NMI) patients developed earlier in the paper: the 10-to-15 percent who meet clinical TSC criteria yet return negative on standard blood-based testing because the causative TSC1/TSC2 variant sits at low VAF in affected tissue, below what a blood draw and a standard pipeline will call. The variant is real; it simply is not in the blood at a detectable fraction.

The recovery of those variants requires two things the biobank uniquely supplies: deep sequencing, and the *right tissue*. The biobank banks resected tuber, AML, and SEGA tissue — exactly the affected lesions where a mosaic variant is enriched. The Variant Curator is designed to be mosaic-aware: low-VAF-sensitive calling (Mutect2 alongside HaplotypeCaller, calling down to VAF ≥ 5%), a mosaic flag carrying VAF, supporting reads, strand balance, and artifact assessment, a ddPCR validation recommendation, and an ACMG-AMP classification synthesis that treats the low-fraction variant as the candidate explanation rather than discarding it as noise. None of that capability does anything without tissue to run it on. The biobank is the tissue.

This is the place for the mantra that organizes the whole CCHMC argument: **a biobank without an intelligence layer is a freezer full of tubes.** The Discover Together Biobank has done the genuinely hard, expensive, slow work — consent, collection, processing, cold-chain, indexing, governance — of accumulating banked tuber, AML, and SEGA specimens from a TSC population. What it does not have, and what no freezer has, is a layer that turns a resected tuber into a draft molecular-genetics interpretation: a deterministic calling-and-annotation pipeline, an evidence-aggregation step, an ACMG synthesis, a provenance trail, a mosaic flag, a validation recommendation, a draft report ready for a board-certified molecular geneticist to sign. The freezer is the kitchen, not the meal. The engine is the cooking. The biobank's second mission — identify disease-causing DNA changes — is precisely the mission an intelligence layer makes executable at scale, and it is the mission whose value is most visibly demonstrated by recovering a mosaic variant that blood testing missed. That demonstration is Act One of §19, run live on synthetic Patient A (the 4-year-old NMI case carrying an 8.3% VAF TSC2 frameshift in tuber tissue, classified Likely Pathogenic on PVS1 + PM2 + PP4 with a ddPCR recommendation). In production, the same path runs on a real banked tuber. In the demo, it does not touch the biobank at all — a boundary we make explicit in §18.6.

### 18.4 The touchpoints architecture: five sources feeding one engine

CCHMC is not a single point of contact for this work. It is five, and they play different roles. The architecture that emerged from the institutional touchpoints treats these five areas as *sources* that feed the engine — upstream suppliers — rather than as components of the engine itself. The data flow is one-directional and simple to state: **source areas → TSC Engine (five agents + deterministic orchestrator) → clinician surfaces.** The sources supply substrate, methodology, patients, and plumbing; the engine transforms; the surfaces deliver. Keeping the sources outside the engine is what keeps the engine portable.

The five source areas, what each supplies, and which agents it feeds:

| # | CCHMC source area | What it supplies | Connects to / role |
|---|---|---|---|
| 1 | **Discover Together Biobank** | Banked tuber / AML / SEGA tissue; the molecular substrate for mosaic recovery | → **TSC-Variant Curator**. The tissue is the input to low-VAF-aware calling and ACMG synthesis; the biobank's second mission (identify disease-causing DNA changes) is what the Curator operationalizes. |
| 2 | **Division of Biomedical Informatics (Dr. Hagedorn)** | The diagnostic-uncertainty discourse-marker methodology; output-surfacing discipline; executive sponsorship | → **TAND Surveillance Agent** (methodology home) and the **clinician surfaces** (how output is surfaced without becoming alert fatigue). Also the institutional sponsor. |
| 3 | **TSC clinical & research program** | Patient concentration; longitudinal multi-subspecialty follow-up; the clinical questions worth answering | → **all five agents**. This is the population the engine models and the clinicians the surfaces serve. |
| 4 | **Epic Clarity / Caboodle + biobank LIMS** | Structured EHR data (encounters, labs, meds, imaging-report text) and specimen metadata | → **Phenome Mapper, Trajectory Modeler, TAND Surveillance Agent** (the data plumbing). **Explicitly NOT built in the synthetic demo** — institutional Phase-1 work. |
| 5 | **Winslow Pavilion (envelope)** | The physical + institutional container for sources 1, 2, and the freezer archive; ~70 research staff; Research Foundation governance (Dr. Tina Cheng) | The envelope, not a feed. Houses the sources; provides the operating norms, governance, and physical infrastructure under which a future real-data deployment runs. |

A few properties of this table are worth drawing out, because they are the difference between an architecture diagram and an honest one.

First, **source 3 (the clinical program) feeds all five agents**, which is just to say that the patients are the point. Every other source exists to serve the clinicians who care for those patients. The Variant Curator answers a geneticist's question; the Phenome Mapper, Trajectory Modeler, TAND agent, and Therapeutics Strategist answer the questions a TSC clinic asks across a patient's life. The surfaces — pre-visit briefing, in-visit dashboard, async alert surface — are built around that clinic's actual workflow, not around the data.

Second, **source 4 is the one that is not built.** Epic Clarity/Caboodle and the biobank LIMS are the real-world plumbing that would feed the Phenome Mapper, Trajectory Modeler, and TAND agent in production. In the synthetic demo, that plumbing is replaced wholesale by the synthetic cohort: Synthea-generated, Clarity-*shaped* structured data and frontier-model clinical notes, deterministically generated and persistently watermarked. The shape is right so that the integration is a swap rather than a rewrite, but the connection to the live institutional systems is described architecturally and not implemented. We are explicit about this rather than hopeful about it.

Third, **source 5 is the envelope, not a feed**, and we keep it in the table precisely so that no one reads "Winslow Pavilion" as a sixth agent or a data source in its own right. The Pavilion is where sources 1 and 2 physically and institutionally live. It supplies governance and operating norms, not data.

The replication property falls directly out of this structure. Because the five sources are external suppliers with well-defined interfaces — tissue in, methodology in, patients in, EHR/LIMS extracts in, governance around — a second institution can be wired up by re-pointing each slot. **Swap the box labels, keep the wiring.** A TGen or a City of Hope brings its own biobank, its own informatics group, its own clinical program, its own Epic/LIMS, its own research-building envelope; the engine, the orchestrator, the agents, the surfaces, and the synthetic-cohort scaffolding are unchanged. The engine is the constant; the sources are the variables. CCHMC is the first set of values, not a special case in the code.

### 18.5 Dr. Hagedorn and Biomedical Informatics: methodology and sponsorship in one place

The second source area deserves its own treatment, because it is the reason one of the five agents exists in the form it does, and the reason the whole project has an institutional sponsor rather than a cold pitch.

Dr. Philip A. Hagedorn, MD, MBI, is the Chief Health Informatics Officer at CCHMC and leads the Division of Biomedical Informatics. His published research program is, in part, about a problem that sounds narrow and turns out to be central to TSC surveillance: detecting *diagnostic uncertainty* in clinical documentation — the linguistic ways clinicians signal that something is not yet settled. The relevant body of work includes Marshall, Nickels, Brady, and Hagedorn (2023, *Journal of Hospital Medicine*); Nickels, Marshall, Edgerton, Brady, Hagedorn, and Lee (2024, *Applied Linguistics*); Ipsaro, Patel, Marshall, and Hagedorn (2021, *Hospital Pediatrics*); and Orenstein and colleagues including Hagedorn (2021, *JAMIA*) on alert burden. The first three establish a discourse-marker approach to uncertainty in notes; the fourth is, in effect, a warning about what happens when you surface findings to clinicians badly.

The TAND Surveillance Agent is a direct extension of this methodology, not an external graft onto it. TSC-Associated Neuropsychiatric Disorders affect roughly 90% of patients and are documented and addressed inconsistently — the 30-to-50-percent under-recognition figure is the headline number, but the mechanism underneath it is often linguistic. A clinician notices a behavioral or academic concern and records it in a way that signals uncertainty rather than action: hedging, deferral, third-party attribution ("mother reports"), conditional framing, or a follow-up plan that never gets formalized. The TAND agent applies what the paper calls the **Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy** — those same marker classes — across longitudinal notes, per the six TAND clusters (behavioral, psychiatric, intellectual, academic, neuropsychological, psychosocial), using Sonnet for per-note discourse analysis, a deterministic scoring and aggregation layer, and Opus only for the briefing summary. It surfaces these under-recognized signals as pre-visit *briefing material*, never as interruptive alerts and never as a diagnosis. The Orenstein-et-al alert-burden work is the reason the surfacing discipline is built in from the start rather than discovered after the fact.

That lineage is the whole point of putting source 2 in the table. We are not asking Hagedorn's division to evaluate a method it has never seen; we are showing it its own method, applied to a disease where the under-documentation problem is acute, with the surfacing discipline its own alert-burden research demands. Hagedorn has taken interest in the work and has offered to engage his team and the faculty TSC lead. That is the sponsorship: methodology and executive support living in the same source area. It is also a credibility constraint — an audience that wrote the papers will not be impressed by overclaiming about them, which is part of why the do-not-overclaim posture (§18.6) is non-negotiable.

### 18.6 The demo-versus-institutional boundary

The most important thing this section can do for a skeptical reader is draw a hard line between what runs now and what is institutional Phase-1 work, and refuse to blur it. The credibility of everything else depends on this boundary being unambiguous.

**What runs now**, on Adam's NVIDIA DGX Spark with extra GPUs via RunPod as needed, is the synthetic-data demonstration: the 50-patient synthetic cohort (version-controlled, deterministic, regenerable in roughly 12 hours), the five agents, the deterministic orchestrator, and the three clinician surfaces. It runs on Synthea-generated Clarity-shaped data, BAMSurgeon-derived genomic substrate (variants inserted into NA12878-derived BAMs at germline ~50% and mosaic 4–12% VAF), frontier-model clinical notes and imaging reports — all persistently watermarked synthetic. No real patient, no real specimen, no live institutional system is touched.

**What is described architecturally but explicitly not built**:

- **Source 4 — Epic Clarity/Caboodle and biobank LIMS integration.** The synthetic cohort is *shaped* like Clarity output so that the eventual integration is a swap, but there is no live extract, no connection to CCHMC's EHR, and no LIMS query. This is institutional Phase-1 work.
- **The biobank tissue path.** The Variant Curator runs on synthetic BAMs in the demo. Running it on a real banked tuber from the Discover Together Biobank is production, and requires the governance and consent path that comes with real specimens.
- **Imaging-AI pipelines.** Imaging *reports* in the demo are frontier-model-generated text. The engine consumes no real DICOM and runs no imaging model. A real imaging-AI capability is separate institutional work.

The governance posture follows from the boundary. The synthetic demo does not require IRB approval, because it touches no human-subjects data; a real-data Phase-1 deployment does require IRB approval and a data-use governance path through the Research Foundation. The platform is not FDA-cleared, and its Software-as-a-Medical-Device posture is undetermined — that determination is itself institutional work, not a claim the demo makes. Every agent output is decision-support behind a human gate: the Variant Curator produces a draft for a board-certified molecular geneticist to sign, never an autonomous diagnosis; the TAND agent surfaces patterns, never diagnoses. Every output carries full provenance — model id and version, prompt-template version, retrieved RAG sources with URIs, input hash, latency — in an append-only, queryable log. The platform is Apache 2.0 and open, with non-commercial intent.

Stating the boundary this plainly is not a hedge; it is the argument. An audience that includes the people who published the alert-burden literature will trust a system that is precise about what it has not built far more than one that gestures at a finished institutional integration. The build is the argument, and the build is the synthetic demo — which is exactly what §19 now walks through, on Patients A, B, and C, in three acts.


## 19. The Demonstration and the Engagement Path

A precision-medicine engine for a rare disease cannot earn the trust of a skeptical informatics division through a slide deck. The argument has to run. Everything described in the preceding sections — the mosaic-aware Variant Curator, the HPO-coded Phenome Mapper, the classical-statistics Trajectory Modeler, the Marshall-Hagedorn TAND Surveillance Agent, the Opus-class Therapeutics Strategist, and the deterministic orchestrator that wires them together — exists so that a single live demonstration can be honest about what it shows. This section describes that demonstration: the three-act arc, the runbook discipline that keeps it credible under failure, the precise boundary between what the demo proves and what it cannot, and the engagement path that follows from Dr. Hagedorn's offer to convene his team and the Cincinnati Children's TSC faculty lead.

The demonstration is not a pitch. It is a discovery instrument run in front of the people who would have to live with the consequences of the build, and the synthetic cohort is what makes that posture possible.

### 19.1 The Three-Act Arc

The demonstration is scoped to 30 minutes of presentation plus 15 minutes of discussion, run on Adam's NVIDIA DGX Spark (GB10 Grace Blackwell, 128 GB unified LPDDR5x) with RunPod GPUs attached for any Parabricks-equivalent variant calling that exceeds the local envelope. Three patients carry the narrative: Patient A (the live mosaic recovery), Patient B (the longitudinal dashboard), and Patient C (the therapeutics brief). All three are version-controlled members of the 50-patient synthetic cohort, and every screen they appear on carries a persistent synthetic-data watermark.

#### Opening — Problem Framing (≈5 minutes)

The opening does no software. It states the three gaps that motivate the engine, each tied to a literature anchor and to a number a clinician will recognize:

1. The no-mutation-identified (NMI) cohort: roughly 10–15% of patients meet clinical diagnostic criteria for TSC yet test negative on standard blood-based genetic testing, most often because the causative variant is a low-VAF somatic mosaic confined to lesional tissue (Tyburczy 2015; Giannikou 2016; Lim 2017). Blood misses it; deep tissue sequencing finds it.
2. The TAND gap: TSC-Associated Neuropsychiatric Disorders affect roughly 90% of patients, yet the TOSCA registry data indicate 30–50% of TAND features go unaddressed in routine care. The signal is in the notes; it is rarely surfaced.
3. The surveillance and trajectory gap: ITSC 2021 consensus surveillance is cadence-based, not patient-trajectory-based, and a slowly enlarging SEGA at the foramen of Monro is exactly the kind of finding that benefits from forecasting rather than fixed-interval imaging.

The opening closes with one honest framing sentence: what runs in the next 25 minutes runs entirely on synthetic data, on hardware that costs less than a mid-range workstation, and is decision-support with a human gate at every output. That sentence is load-bearing. It sets the credibility contract for everything that follows.

#### Act One — Live Mosaic Recovery on Patient A (≈8 minutes)

Patient A is a 4-year-old girl in the NMI cohort: clinically diagnosed, blood-test negative, with an 8.3% VAF TSC2 frameshift variant present only in resected tuber tissue. This is the act that is run live, end to end, because it is the single most defensible claim the engine makes and the one a molecular geneticist in the room can check.

The operator submits Patient A's tuber-tissue BAM to the TSC-Variant Curator. The deterministic stages run first and visibly: Parabricks/BWA-MEM alignment context, GATK HaplotypeCaller plus Mutect2 in mosaic-aware mode with the low-VAF threshold set at ≥5%, then snpEff/VEP annotation. The pipeline recovers the 8.3% VAF frameshift, flags it as mosaic with the supporting evidence (VAF, supporting reads, strand balance, artifact assessment), and recommends ddPCR orthogonal validation. The Sonnet evidence-aggregation step and the Opus ACMG-AMP synthesis then produce the classification — Likely Pathogenic on PVS1 + PM2 + PP4 — as an AI-labeled draft molecular-genetics report explicitly marked for board-certified molecular-geneticist sign-off.

The act does not end at the answer. It ends at the audit trail. The operator opens the provenance record for the call: model id and version, prompt-template version, the input BAM hash, retrieved RAG sources with URIs, and per-stage latency, all from the append-only event log. The point being made to the informatics audience is not "the AI got it right." The point is "you can see exactly how it got there, and a human still signs it." The eval target this act demonstrates against synthetic ground truth — recover the mosaic variant at VAF ≥5%, correct ACMG class, no false-positive Pathogenic call, under 5 minutes per case — is stated out loud so the audience knows what success looked like before they saw it.

#### Act Two — Patient B Dashboard, Patient C Therapeutics, Pre-Visit Briefing (≈12 minutes)

Act Two shifts from the single dramatic recovery to the longitudinal, multi-agent picture that is the engine's actual day-to-day value.

**Patient B** is a 12-year-old boy with a TSC2 c.3037C>T (p.Arg1013Ter) variant, a SEGA at the foramen of Monro measured at 0.8 → 1.1 → 1.3 cm across the synthetic imaging series, bilateral angiomyolipomas (≈2.8 cm), well-controlled focal seizures, and scattered under-recognized TAND signals embedded in his notes. He drives the in-visit four-quadrant dashboard: variant interpretation, the HPO timeline from the Phenome Mapper, the Trajectory Modeler's SEGA/AML/seizure/renal forecasts at 6/12/18 months with 50% and 90% prediction intervals, and the combined TAND-plus-therapeutics quadrant. The Trajectory Modeler is classical statistics — mixed-effects models, Gaussian-process regression, survival analysis, Bayesian updating — not an LLM, and the demonstration says so explicitly, because a forecast a clinician can interrogate beats a forecast a clinician must take on faith. The target moment is the model placing Patient B's SEGA threshold crossing inside the 12–18 month window without raising a false alarm, with progressive disclosure letting the clinician drill from the forecast to the source measurements.

The TAND quadrant for Patient B is where the Marshall-Hagedorn lineage shows. The Surveillance Agent has run Sonnet per-note discourse analysis across the six TAND clusters, applying the diagnostic-uncertainty discourse-marker taxonomy — hedging, deferral, third-party attribution, conditional language, follow-up-without-formalization — that comes directly from Hagedorn's published clinical-NLP work (Marshall/Nickels/Brady/Hagedorn 2023; Nickels/Marshall/Edgerton/Brady/Hagedorn/Lee 2024). The surfaced patterns appear as pre-visit briefing material, never as interruptive alerts and never as a diagnosis. This is the moment to say plainly, to the person whose method it is: the TAND agent is an extension of your own published research, not an external graft.

**Patient C** is an 18-year-old woman with a TSC1 variant, a partial everolimus response complicated by mucositis that forced a dose reduction, an AML near 4 cm, and refractory focal seizures. She drives the Therapeutics Strategist. The Opus-class agent integrates the four prior agents' outputs plus medications, adherence, and adverse events with PubMed/PMC RAG, a ClinicalTrials.gov snapshot, and FDA actions, producing the six-section structured options brief: Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions. Every claim is source-attributed. The brief is framed as decision-support, not recommendation, and the demonstration highlights the agent's hedging where the evidence is thin — appropriate uncertainty is a feature being shown, not an embarrassment being hidden. The eval target: correct trial matches against ground truth, appropriate hedging, full attribution, under 3 minutes.

Act Two closes on the **pre-visit briefing surface** — the one-screen, mobile-readable assembly the orchestrator produces on demand: header, what's-new, zero-to-three action items, a watchlist, and source links. This is the surface a busy TSC clinic would actually touch first, and showing it last in the act makes the point that all five agents exist to feed a single calm screen rather than five competing dashboards.

#### Act Three — Infrastructure, Cost, Scaling, License (≈5 minutes)

The final act answers the question the informatics division is already asking. It states the hardware: a single DGX Spark, RunPod GPUs rented for variant-calling bursts and synthetic-cohort regeneration. It states the cost honestly — frontier-model API calls plus modest rented GPU time, not a seven-figure platform. It states the reused substrate (HCLS AI Factory v1.3.0: LangGraph orchestration, Milvus RAG, tiered Claude models with a local Llama 3.1 70B fallback, PostgreSQL/Redis/MinIO, audit/provenance) versus the net-new TSC work (the multi-agent orchestrator and shared event-sourced state, the synthetic cohort pipeline, the five agents).

Then it states the scaling architecture in the institution's own terms. The engine is the kitchen, not the meal; a biobank without an intelligence layer is a freezer full of tubes. The five Cincinnati Children's areas — the Discover Together Biobank's banked tuber/AML/SEGA tissue, Biomedical Informatics and the TAND methodology, the TSC clinical and research program's patient concentration, and the Epic Clarity/Caboodle plus LIMS data plumbing — are **sources that feed the engine, not parts of it**. Replication to a TGen or a City of Hope is "swap the box labels, keep the wiring." The act ends on the license: Apache 2.0, open, non-commercial intent. The closing sentence returns to the contract from the opening: this is the synthetic-data build that runs today; the Epic/Clarity/LIMS and imaging-AI integrations are described architecturally and are explicitly institutional Phase-1 work, not built.

### 19.2 Runbook Discipline

A live demo of a multi-agent system is a managed risk, and the runbook treats it that way.

**Dry runs.** Week 8 of the build is reserved for cohort regeneration, full dry runs, clinician review of the synthetic notes, and delivery. The dry runs are end-to-end timed walkthroughs against the frozen, version-controlled cohort, so the latency numbers stated in each act (under 5 minutes for Patient A's variant call, under 3 minutes for Patient C's brief) are observed numbers from rehearsal, not aspirations.

**Failure-mode handling.** The orchestrator's conservative failure posture is what makes a live run safe: if an agent fails, the surface it feeds shows "pending" or a staleness marker — never a silently missing output. In a demonstration this matters twice over. If the live Act One variant call hits a transient RunPod or API problem, the runbook falls back to the pre-computed, provenance-complete result for Patient A captured during the final Week-8 dry run, and the operator says so explicitly rather than pretending. Acts Two and Three run against pre-materialized current-state projections from the event log, so a live failure in one agent degrades a single quadrant to "pending" rather than collapsing the dashboard. The discipline is that the demo never hides a failure; surfacing the failure honestly is itself part of what the engine is supposed to do in production.

**Watermarking.** Every screen, every exported brief, and every PDF generated during the demo carries the persistent synthetic-data watermark. The cohort contains no real imaging or DICOM, no raw FASTQ (it starts at BAM), no real neuropsych scores, no pedigree beyond the notes, and no pharmacy or claims data. The synthetic clinical notes were generated from published templates and passed through clinician-sampled review. Stating these exclusions during Act Three is not a disclaimer to rush past; with this audience it is a credibility asset, because it tells them precisely what the engine has and has not been exposed to.

### 19.3 What the Demonstration Proves — and What It Cannot

The demonstration is evaluated against synthetic ground truth. That is a real, falsifiable test of the engineering, and it is not clinical validation. The distinction must be stated plainly in the room.

**What it proves.** The engine can recover all seven mosaic variants in the cohort at VAF ≥5% with correct ACMG classification and no false-positive Pathogenic calls; the Phenome Mapper can hit recall ≥90% and precision ≥85% against the HPO ground truth with a cohort pass under an hour; the Trajectory Modeler can place a known threshold crossing in the correct forecast window without false alarms; the TAND agent can detect the embedded discourse signals without spurious flags; the Therapeutics Strategist can match the correct trials with appropriate hedging and full attribution. It proves the orchestration is deterministic and dependency-ordered, that state is event-sourced and queryable, that provenance is attached to every output, and that the whole thing runs on a single Spark plus rented GPUs.

**What it cannot prove.** It cannot prove clinical accuracy on real patients — synthetic ground truth is a generator's truth, not a patient's. It cannot prove generalization to messy production data: the synthetic notes are cleaner and the synthetic genomes more controlled than real Epic Clarity exports and real tissue BAMs. It says nothing about regulatory status — the engine is not FDA-cleared, the SaMD posture is undetermined and is institutional work, and IRB review (unnecessary for the synthetic demo) is required before any real-data Phase 1. And it does not prove that the Epic/Clarity/LIMS plumbing or the imaging-AI pipelines work, because those are not built. The demonstration's honesty about this second list is what gives the first list its weight.

### 19.4 The Engagement Path

Dr. Hagedorn has taken an interest in the work and offered to convene his team and the Cincinnati Children's TSC faculty lead. That offer sets the posture for everything after the demo: this is discovery, not a sale.

**Discovery, not pitch.** The convening's purpose is to let the people who own the data, the methodology, and the patients tell us where the synthetic model is wrong. The TAND agent's lineage in Hagedorn's own publications makes him a critic with standing, not a prospect; the right outcome of the meeting is a list of ways the discourse-marker taxonomy needs to be adapted to real CCHMC documentation, not a signature.

**Synthetic-data de-risking.** Because the entire demonstration runs on watermarked synthetic data on hardware Adam already owns, the conversation can happen with zero institutional data exposure, no IRB dependency, and no procurement commitment. There is nothing to protect in the demo cohort, which is exactly why it can be shown freely. The institution incurs risk only if and when it chooses to move to real data.

**A scoped next step.** The natural, narrow next step after the convening is a single-area Phase-1 scoping engagement rather than a platform deployment — most plausibly the Discover Together Biobank's banked tuber/AML/SEGA tissue feeding a real-data Variant Curator pilot under IRB, since that is the area where the synthetic demo's strongest claim (mosaic recovery) meets the institution's strongest asset (the freezer archive). That keeps the first real-data step small, governed, and aligned to a question the institution already cares about: how much of its own NMI cohort is hiding a recoverable mosaic variant.

Which leads directly to the governance the engine must carry into any real-data step — the audit, oversight, and human-gate model that the next section makes explicit.


## 20. Governance, Ethics, Regulatory & Privacy

The preceding section described how the engine behaves; this one describes the rules it must live under, and — more importantly — names the places where those rules are not yet satisfied. With an audience like Hagedorn's, the candor is the credibility. A clinical informatics team that has spent years detecting diagnostic uncertainty in documentation will detect it instantly in a vendor pitch. So this section is written the way the rest of the paper is written: the build is the argument, and the seams are stated out loud rather than papered over. Several of the items below are explicitly *not* solved in the synthetic-data demonstration and are flagged as institutional Phase-1 work. That distinction — what runs now on a DGX Spark with synthetic patients versus what requires a hospital, an IRB, and a CLIA-certified laboratory — runs through every subsection.

### 20.1 What the demonstration is, and the privacy posture that follows

The artifact that runs today is a 50-patient synthetic cohort generated once, version-controlled, deterministic, and persistently watermarked as synthetic on every surface. No real patient touches it. Every clinical note, imaging report, and variant call descends from Synthea skeletons, BAMSurgeon-edited NA12878-derived BAMs, and frontier-model generation from published templates. The featured patients — A (4yo, NMI mosaic recovery), B (12yo, SEGA trajectory), C (18yo, everolimus optimization) — are fictional composites, not de-identified real people. This matters legally: synthetic data that was never derived from a real individual is not protected health information (PHI), and the synthetic demo therefore sits outside HIPAA, outside the Common Rule, and outside the need for an IRB protocol. That is the correct posture for an MVP, and it is also the *only* honest claim available, because nothing in the current build has the data plumbing to ingest a real record.

The moment the engine touches a real Cincinnati Children's (CCHMC) patient, every line below changes. Epic Clarity/Caboodle extraction, biobank LIMS linkage, and any imaging-AI pipeline are described architecturally in this paper and are **explicitly not built**. They are Phase-1 institutional work, and they are the trigger that converts almost every "not applicable to the demo" answer here into a hard requirement.

### 20.2 HIPAA, the data plane, and where it bites

For the synthetic demo HIPAA does not apply, and we do not claim a HIPAA-compliant architecture — claiming compliance for data that is not PHI is a category error and reviewers notice it. What we *can* state is that the substrate was designed so that the real-data transition is an integration task rather than a re-architecture: the event-sourced PostgreSQL log, the Redis TTL layer, MinIO object storage, and the provenance record on every agent output are the same primitives a covered-entity deployment would harden.

When real data arrives, the obligations are concrete and CCHMC-specific:

- **Minimum necessary.** The Phenome Mapper consumes Clarity-shaped structured fields and clinical notes; the Variant Curator consumes BAMs. Each agent should receive only the data classes it needs, enforced at the extraction layer, not filtered after the fact inside the engine.
- **De-identification vs. limited data set.** A research deployment under an IRB will likely operate on a limited data set with dates and ages retained (the Trajectory Modeler is useless without real dates) under a data use agreement, not on a Safe Harbor de-identified set. This is a decision for CCHMC's privacy office, not for us to assert.
- **Audit and access.** The append-only event log already records who/what wrote each state transition; mapping that to authenticated CCHMC identities and to the institution's access-review cadence is Phase-1.
- **Boundary of the standalone apps.** The three clinician surfaces are standalone web apps, not Epic-embedded. In a real deployment that means a second system holding clinical data, which the institution's security review must treat as in-scope. We do not pretend the standalone posture makes the privacy question smaller; it makes it a distinct attack surface.

### 20.3 IRB and the research/operations line

There is no IRB protocol for the synthetic demo, and there should not be — there are no human subjects. Phase-1, on real CCHMC patients, requires IRB review, and the protocol design has a genuine subtlety the team should anticipate: parts of this engine could be framed as quality improvement / operations (surfacing surveillance gaps already mandated by the ITSC 2021 consensus guidelines) and parts are unambiguously research (modeling SEGA/AML trajectories, validating mosaic recovery against tissue ground truth). The honest reading is that the engine as a whole is research and should go through full IRB review rather than seeking a QI carve-out for the convenient pieces. The Discover Together Biobank's existing consent framework and the institution's biorepository governance are the natural anchors for the tissue side; Hagedorn's Division of Biomedical Informatics is the natural sponsor for the documentation-NLP side.

### 20.4 Pediatric consent: assent, parental permission, and the transition no one schedules

TSC is overwhelmingly a pediatric disease at presentation, and a children's hospital is the deployment site, so the consent model is not the standard adult model and cannot be treated as a footnote.

- **Parental permission + assent.** Enrollment of a minor requires parental (or guardian) permission, plus the child's assent where developmentally appropriate. TAND means a meaningful fraction of this population has intellectual disability (TAND affects ~90% of patients, with 30–50% of features unrecognized per TOSCA), which complicates assent capacity and argues for documented, age- and capacity-appropriate assent procedures rather than a checkbox.
- **The minor-to-adult transition.** This is the seam clinicians will press on. A patient enrolled at 4 (Patient A) reaches the age of majority while the longitudinal record — and the trajectory model trained on it — continues. At 18 the patient must be re-consented as an adult in their own right; permission granted by a parent does not silently convert. The engine's event-sourced design helps here (consent state is itself an event with effective dates), but the *governance* obligation is a re-consent workflow, and we flag that the demo does not implement it because the demo has no real subjects to re-consent.
- **Family-member consent and incidental genetic findings.** Roughly one-third of TSC is inherited, and the Variant Curator works on tissue/germline-equivalent variants. A pathogenic TSC1/TSC2 finding has implications for parents and siblings. Any real deployment needs a defined policy — coordinated with CCHMC genetic counseling — on whether and how relatives are addressed, and the consent must contemplate familial implications up front. This is squarely Phase-1 and is not modeled in the synthetic cohort.

### 20.5 Return of results

Return-of-results policy is where the augment-not-replace principle becomes operational. The Variant Curator produces an AI-labeled *draft* molecular-genetics report; it never returns a result to a patient and never finalizes one. Under any real deployment, results reach patients only through the existing clinical pathway: a board-certified molecular geneticist signs the report, and genetic counseling handles disclosure, including the familial and reproductive implications and the minor's own future interest in their genetic information. The same gate applies to the Therapeutics Strategist's options brief — it is decision-support for a clinician, never a recommendation delivered to a patient.

A specific sub-question deserves naming: the engine recovers low-VAF somatic mosaic variants in tissue (the ~10–15% NMI cohort). A mosaic finding in tuber tissue carries a different counseling message than a germline finding — recurrence risk, transmissibility, and the limits of what a tissue-restricted variant predicts are genuinely different conversations. The engine flags mosaicism (VAF, reads, strand, artifact assessment) and recommends ddPCR validation precisely so the geneticist and counselor have what they need; it does not attempt to script the disclosure.

### 20.6 Algorithmic bias and subgroup monitoring

Three of the five agents have model-dependent components (Variant Curator's evidence aggregation and ACMG synthesis; Phenome Mapper's note extraction; TAND's discourse analysis; the Trajectory Modeler is deliberately classical statistics and the Therapeutics Strategist is RAG-grounded). Each model-dependent component is a potential vector for inequitable performance.

The honest limitation up front: a 50-patient synthetic cohort cannot characterize subgroup performance. It can demonstrate function; it cannot establish fairness. Several bias surfaces are specific and worth stating:

- **Documentation bias propagates.** The Phenome Mapper and TAND agents read clinical notes. If TAND features are systematically under-documented for some groups — and the 30–50% miss rate is not evenly distributed across language, socioeconomic status, or caregiver advocacy — an agent trained to surface what is *written* will inherit and can amplify that gap. The Marshall-Hagedorn discourse-marker taxonomy detects diagnostic uncertainty in the text; it cannot detect a concern that was never recorded.
- **Reference-database ancestry skew.** Variant interpretation leans on gnomAD v4, ClinVar, and LOVD-TSC. Population allele-frequency evidence (e.g., PM2) is less reliable for under-represented ancestries, which can bias ACMG classification toward or away from "pathogenic" by ancestry. This is a known, field-wide problem; we do not solve it, we flag it and route it to the human geneticist.
- **The plan, not the claim.** Phase-1 must include a pre-specified subgroup monitoring protocol (by age band, sex, ancestry where ethically collectable, and primary language for the NLP agents), with performance reported per subgroup against ground truth, and a stopping/escalation rule if a subgroup falls below target. We commit to the protocol; we do not claim fairness the demo cannot support.

### 20.7 SaMD and FDA posture

The engine is **not FDA-cleared**, and its Software-as-a-Medical-Device (SaMD) classification is **undetermined**. That is a deliberate, accurate statement, not an evasion. Whether any component meets the device definition turns on intended use and claims, and those are exactly what an institution decides as it moves from synthetic demonstration to clinical use. Determining SaMD posture, and any predicate or De Novo pathway analysis, is institutional Phase-1 regulatory work that would be done with CCHMC's regulatory affairs function and is out of scope for this paper to assert.

What we *can* say about how the design lowers regulatory risk rather than raising it:

- Every model-dependent agent output is framed as decision-support with a mandatory human gate. The Variant Curator drafts; the geneticist signs. The Therapeutics Strategist presents options with attribution; the clinician decides. Human-in-the-loop decision-support that does not drive automatic action is a materially different regulatory object than autonomous diagnosis.
- The Trajectory Modeler is the component most likely to look device-like (it forecasts SEGA/AML growth with prediction intervals and emits threshold-crossing alerts). We flag it as the highest-attention item for a future SaMD analysis and note that its outputs are surveillance-cadence recommendations to a clinician, not orders.

We make no clinical-validation claim of any kind. The eval targets in this paper are measured against *synthetic ground truth* and are statements about software behavior, not clinical performance.

### 20.8 Two boundary questions reviewers will raise

**Variant interpretation vs. CLIA.** The Variant Curator performs variant calling and ACMG-AMP classification synthesis on research-grade synthetic data. It is not a clinical laboratory and produces nothing CLIA-reportable. In a real deployment, clinical molecular results must originate from a CLIA-certified, CAP-accredited laboratory; the engine's role is to *draft* and *organize evidence for* a result that the certified lab and its molecular geneticist own and sign. The ddPCR validation recommendation it emits is precisely the hand-off to that orthogonal, validated assay. We state plainly that the engine does not and must not substitute for CLIA-regulated testing.

**TAND surveillance vs. psychiatric screening.** The TAND Surveillance Agent is the most ethically sensitive component because it touches mental-health and neuropsychiatric domains in children. The design constraints are governance choices, not incidental UX:

- It surfaces *patterns already present in the documentation* as pre-visit **briefing material**, never interruptive alerts (Orenstein/.../Hagedorn 2021 on alert burden is a direct motivation), and **never a diagnosis or a screening score**.
- It is not a validated psychiatric screening instrument and does not impersonate one. It does not generate risk scores for self-harm or any acute psychiatric outcome; if such signals appeared in real notes, the appropriate response is the existing clinical pathway, not an algorithmic flag.
- Its discipline includes recalibration if it produces more than roughly three alerts per clinician per week, the explicit guard against the alert-fatigue failure mode its own intellectual lineage documented.

This boundary is why the TAND agent is positioned as a direct extension of Hagedorn's published diagnostic-uncertainty work rather than as a novel psychiatric tool: it detects *how clinicians write about uncertainty*, which is a documentation phenomenon, not a patient diagnosis.

### 20.9 Indemnification, liability, and the open-source gap

This is the seam most pitches hide. The platform is Apache 2.0, open-source, with non-commercial intent. The Apache 2.0 license disclaims warranties and limits liability — appropriately for source code, and bluntly inadequate as the basis for clinical deployment. There is **no indemnification** flowing from this project to a deploying institution, no clinical liability coverage, and no vendor standing behind a result. An institution that puts any version of this in front of real patient care assumes the corresponding clinical and operational liability itself, through its own malpractice, governance, and risk structures, with its own clinicians as the accountable decision-makers at every human gate.

We state this rather than soften it because the human-gate architecture is, in part, a response to it: liability tracks the clinician who signs the report and makes the decision, which is exactly where TSC care already places it. The engine does not insert itself into that chain of accountability; it organizes evidence for the human who is already in it. Any move toward a posture with stronger guarantees — a commercial entity, indemnification, a SaMD submission — is an institutional decision well beyond this paper, and we do not gesture at it as though it were planned.

### 20.10 Single-maintainer sustainability

The final honest seam is organizational rather than legal. This engine, and the HCLS AI Factory substrate it reuses (v1.3.0: LangGraph, Milvus, tiered Claude models with a local Llama 3.1 70B fallback, PostgreSQL/Redis/MinIO), is built and maintained by one person on one DGX Spark with RunPod burst capacity. That is a feature for a fast, credible MVP and a risk for anything an institution would depend on. We name it directly:

- **Bus factor.** A single maintainer is a single point of failure for security patching, dependency upkeep, and model-deprecation migrations (the model tiers are versioned in provenance precisely so a migration is auditable, but migration still requires a maintainer).
- **The mitigation is the partnership, not a promise.** The path out of single-maintainer risk is the Winslow Initiative described next — engaging Hagedorn's team and the faculty TSC lead so that the work becomes institutionally owned rather than personally heroic. Hagedorn has offered to engage his team and the TSC lead; the right reading of that offer is as the beginning of a sustainability answer, not a validation badge.
- **Reproducibility as insurance.** The deterministic, version-controlled synthetic cohort (~12-hour regeneration), the event-sourced state, and the provenance-on-every-output design mean the system's behavior is reconstructable by someone other than its author. That is the most a single maintainer can responsibly offer, and we offer it.

Taken together, these seams are not reasons the engine should not be built; several are reasons it should be built *this* way — synthetic-first, human-gated, provenance-complete, and open. They are also the reason the next section exists. Closing the gaps named here — IRB, real-data consent including the pediatric transition, CLIA hand-off, subgroup monitoring, SaMD determination, and institutional ownership that retires the bus-factor risk — is not engineering work that fits in the eight-week build. It is the institutional program that the Winslow Initiative lays out.


## 21. The Winslow Intelligence Initiative: A 24-Month Institutional Path

The synthetic demo described in the preceding sections is the argument, not the destination. It runs today on a DGX Spark with RunPod overflow, on a 50-patient cohort that contains no real patient, no real imaging, and no real molecular result. That constraint is deliberate and it is honest: it lets us show the wiring without touching protected data or making clinical claims we cannot yet support. But the demo's value is precisely that it makes the next question concrete. If the engine recovers all seven synthetic mosaic variants, surfaces TAND signals a clinician would otherwise have to assemble by hand, and forecasts Patient B's SEGA crossing a surveillance threshold, then the institution's real question is no longer "does this work in principle" but "what would it take to point it at our own patients, our own banked tissue, and our own Epic instance, safely."

This section answers that question as a 24-month phased plan we call the Winslow Intelligence Initiative, named for the Winslow Research Pavilion that houses the Discover Together Biobank, the Gamble Vaccine Research Center, and the institution's first centralized biospecimen freezer archive. The Pavilion is the infrastructure envelope; the engine is the intelligence layer that envelope has so far lacked. A biobank without an intelligence layer is a freezer full of tubes. The plan below is how those tubes acquire a query interface — incrementally, with decision gates, and with the honesty that the synthetic demo bought us.

Each phase ends with an explicit go/no-go decision gate. Phases do not auto-advance. Budget figures are indicative ranges intended for planning conversations, not quotes; they assume the open-source HCLS AI Factory v1.3.0 substrate (Apache 2.0) is reused rather than rebuilt, that compute is a mix of the existing Spark and institutional or RunPod GPU, and that the dominant cost is people, not licenses.

### Data classes and the honest broker

Before any phase touches real data, one structural commitment governs everything: an honest-broker boundary between identified clinical/genomic source data and the engine's working environment. The engine never sees direct identifiers. CCHMC's Biomedical Informatics group, under Dr. Hagedorn, is the natural home for the honest-broker function, which fits its existing role in clinical-NLP and data-de-identification work.

We distinguish four data classes, because they have different governance, different latency, and different IRB posture:

| Class | Examples | Source area | Where it lives | Gate |
| --- | --- | --- | --- | --- |
| A — Reference/public | ClinVar, gnomAD v4, LOVD-TSC, HPO, PubMed/PMC, ClinicalTrials.gov, FDA actions | Public | Milvus TSC partition, local mirror | None |
| B — Synthetic | The 50-patient cohort, watermarked | Generated | Engine workspace | None (IRB-exempt) |
| C — De-identified real | Banked-tissue VCFs, de-identified notes, longitudinal labs/imaging reports via honest broker | Biobank LIMS, Epic Clarity/Caboodle | Limited-access enclave | IRB; honest broker |
| D — Identified real | Linkage to a living patient, return of a draft molecular report into the chart | Epic Genomics, molecular-genetics service | Never in engine; only the human-gated output crosses back | IRB; molecular-geneticist sign-off; institutional SaMD review |

The engine operates on Class A, B, and C. Class D is only ever the human-reviewed output crossing back into the clinical world through a board-certified molecular geneticist (for variant interpretation) or a treating clinician (for surfaces) — never an autonomous write. This boundary is what makes the surfaces "standalone web apps, not Epic" a feature rather than a limitation in early phases: the engine can be useful against Class C data long before it is cleared to write anything into a chart.

### Phase 0 — Foundation and discovery (months 0–3)

Phase 0 builds nothing for production. It establishes the conditions under which building is responsible.

**Sponsorship and governance.** The Initiative needs a named executive sponsor and a faculty TSC clinical lead in addition to Dr. Hagedorn's informatics sponsorship. A small steering committee — informatics, the TSC clinical/research program, the Discover Together Biobank, molecular genetics, and a research-IRB liaison — meets monthly and owns the decision gates. The committee's first job is to convert the demo's implicit architecture into an institutional product requirements document (PRD): which agents deploy in which order, against which data classes, on which surfaces, with which success criteria.

**Baseline.** We cannot claim improvement without a baseline. Phase 0 measures, retrospectively and from existing data, two things: (1) the size and current disposition of the institution's no-mutation-identified (NMI) TSC cohort — clinically diagnosed patients with negative blood-based genetic testing, the candidates for mosaic recovery from banked tissue; and (2) a documentation baseline for TAND recognition, i.e., how often the six TAND clusters are addressed in existing notes, which is the denominator against which the TOSCA-reported 30–50% miss rate will be checked locally.

**IRB and data governance.** A protocol covering Class C de-identified secondary use of biobank tissue and Clarity/Caboodle extracts is drafted and submitted. The synthetic demo required no IRB; everything past it does.

**Decision gate G0.** Proceed to Phase 1 only if: sponsor and faculty TSC lead are named; the steering committee is seated; an NMI cohort of at least ~30 historical cases with retrievable banked tissue is confirmed; and the IRB protocol is submitted with a credible path to approval. Indicative budget: **$0.3–0.6M**, almost entirely a fraction of existing staff time plus IRB and governance overhead.

### Phase 1 — Variant Curator and the foundation layer (months 3–9)

Phase 1 is deliberately the narrowest, most defensible, most measurable thing the engine does: somatic-mosaic variant recovery from banked tissue, against the cohort whose standard blood testing already came back negative. This is where the engine has the clearest right to exist, because the comparator is "no answer."

**What gets built and integrated.** The TSC-Variant Curator moves from the synthetic BAMSurgeon substrate to real, de-identified banked-tissue sequencing. The deterministic core — Parabricks/BWA-MEM alignment, GATK HaplotypeCaller plus Mutect2 mosaic-aware calling at VAF ≥5%, snpEff/VEP annotation — is unchanged from the demo; the work is plumbing it to the biobank LIMS for specimen provenance and to Epic Genomics for the existing negative blood result, so the AI-labeled draft molecular-genetics report carries full provenance and the mosaic flag (VAF, supporting reads, strand bias, artifact assessment) and a ddPCR validation recommendation. The Sonnet evidence-aggregation and Opus ACMG-AMP classification-synthesis layers are validated against the institution's existing combinatorial-rule pipeline before any draft reaches a geneticist.

This is the phase where the Epic and LIMS integrations the demo explicitly did not build actually get built — for one agent, against one data class, behind the honest broker.

**The re-analysis.** The substantive deliverable is a blinded re-analysis of at least 30 historical NMI cases for which banked tuber, AML, or SEGA tissue exists. Every candidate variant the engine surfaces is reviewed by a board-certified molecular geneticist; nothing is autonomous. The clinical question is simple and consequential: in how many genuinely NMI patients does deep tissue sequencing with mosaic-aware calling and AI-synthesized ACMG classification recover a likely-pathogenic or pathogenic variant the blood test missed?

**First paper.** Phase 1 produces the first peer-reviewed output: the mosaic-recovery study, methods and results, with the engine described honestly as decision support whose every output a geneticist signed. Literature anchors — Tyburczy 2015, Giannikou 2016, Lim 2017 — establish that the biology is real; the contribution is an operationalized, provenance-complete, reproducible pipeline against a defined NMI cohort.

**Decision gate G1.** Proceed only if: the geneticist-adjudicated recovery rate is non-trivial and the false-positive Pathogenic rate is essentially zero (the synthetic eval target of no false-positive Pathogenic calls must hold against real tissue); turnaround per case is clinically reasonable; and the molecular-genetics service judges the draft reports useful rather than additional burden. A null or weak result is itself publishable and is a legitimate reason to stop or rescope — see the negative-result commitment below. Indicative budget: **$0.8–1.5M** (sequencing of banked specimens dominates the non-personnel cost; ddPCR confirmation; engineering for LIMS/Epic-Genomics integration; geneticist review time).

### Phase 2 — Phenome, Trajectory, TAND, and SMART-on-FHIR surfaces (months 9–18)

Phase 2 adds the longitudinal-intelligence agents and, for the first time, puts a surface in front of clinicians during real care.

**Agents.** The TSC-Phenome Mapper extracts time-anchored HPO-coded phenotypes from de-identified Clarity/Caboodle structured data and notes (Sonnet per-note; Haiku for ICD-10/lab-to-HPO normalization; Opus for rare conflict resolution), producing the longitudinal HPO profile, discordance log, and ITSC surveillance-gap report that the other agents build on. The TSC-Trajectory Modeler — classical statistics, not LLM-driven: mixed-effects models, Gaussian-process regression, survival analysis, Bayesian updating — forecasts SEGA/AML growth and seizure burden at 6/12/18 months with 50%/90% prediction intervals and threshold-crossing alerts. The TAND Surveillance Agent applies the Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy (hedging, deferral, third-party attribution, conditional language, follow-up-without-formalization) over the six TAND clusters, with a deterministic scoring layer and an Opus briefing summary. The TAND agent is, by design, a direct extension of Dr. Hagedorn's own published clinical-NLP work (Marshall/Nickels/Brady/Hagedorn 2023; Nickels et al. 2024; Ipsaro et al. 2021) rather than an external graft — which is both why it is methodologically credible and why it belongs at CCHMC first.

**Surfaces and the SMART-on-FHIR step.** In Phase 2 the three clinician surfaces — pre-visit briefing, in-visit four-quadrant dashboard, async alert surface — graduate from standalone watermarked web apps to SMART-on-FHIR applications that can launch in the Epic context and read FHIR resources, while still writing nothing back autonomously. TAND output remains pre-visit briefing material, never an interruptive alert, and never a diagnosis.

**The alert-burden guardrail.** This is the phase's most important non-feature. Orenstein et al. 2021 (JAMIA) — again, Hagedorn-lineage work — documents how alert burden degrades clinical attention, and the engine's async alert surface is governed by a hard discipline: across all four alert categories, if a clinician receives more than roughly three alerts per week, the thresholds recalibrate. This guardrail is measured and reported as a primary safety/usability outcome of Phase 2, not an afterthought. A surface that clinicians silence is worse than no surface.

**Validation against the local baseline.** The Phenome Mapper is held to the synthetic eval targets translated to real data (recall ≥90%, precision ≥85%, sampled against chart review); the Trajectory Modeler's prediction intervals are checked for calibration against observed growth; the TAND agent's surfaced signals are reviewed against the Phase-0 documentation baseline to estimate how many genuinely under-recognized features it brings forward, and at what false-positive cost.

**Decision gate G2.** Proceed only if: surface adoption is real (clinicians open the briefing before TSC visits); the alert-burden guardrail holds without the engine becoming useless; Phenome/Trajectory metrics meet the translated targets; and TAND signals are judged clinically meaningful and non-noisy by the TSC clinical team. Indicative budget: **$1.5–3.0M** (the largest phase: SMART-on-FHIR engineering and Epic app review, ongoing API/inference cost across three model tiers, statistical validation, and clinician evaluation time).

### Phase 3 — Therapeutics Strategist, research engagement, and first external collaboration (months 18–24)

Phase 3 closes the loop with the highest-stakes agent and turns the Initiative outward.

**The Therapeutics Strategist.** This Opus-class agent integrates all four prior agents plus medications and adherence, adverse events, the PubMed/PMC RAG corpus, a ClinicalTrials.gov snapshot, and FDA actions into the six-section structured options brief — Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions. Every claim is source-attributed; the brief is framed throughout as decision support, not recommendation. In a real deployment its most concrete near-term value is trial matching: connecting patients to mTOR-inhibitor and next-generation selective mTORC1 protocols, and to targeted-epilepsy trials in the lineage of EXIST-3, for which the institution's TSC concentration makes it an unusually strong site. The Strategist writes nothing; it assembles, attributes, and hedges, and a clinician decides.

**Research engagement.** With the foundation, longitudinal, and therapeutic layers live against Class C data, the engine becomes a research instrument. The Discover Together Biobank's three stated missions — find biomarkers, identify disease-causing DNA changes, and understand resilience — map directly onto Phase-1 mosaic recovery, Phase-2 phenotype/trajectory modeling, and the cross-modal questions the full engine can now ask. Phase 3 funds the first investigator-initiated studies that use the engine as infrastructure rather than studying the engine itself.

**First external collaboration.** Phase 3 also proves replication by exporting the architecture, not the data, to one external partner — a TGen or a City of Hope. The replication principle is "swap the box labels, keep the wiring": the five CCHMC source areas (biobank, informatics/Hagedorn methodology, TSC clinical program, Epic/LIMS plumbing) are institutional roles, and a partner site fills the same roles with its own people and systems. No patient data crosses institutions; the Apache-2.0 codebase and the architecture do. Section 22 develops replication in full; here it is simply the Phase-3 gate's outward proof.

**Decision gate G3.** Proceed to sustained operation if: the Therapeutics briefs are judged useful and appropriately hedged by treating clinicians; at least one investigator-initiated study is funded and running on the engine; and one external site has stood up the architecture against its own data. Indicative budget: **$1.2–2.5M** (Therapeutics RAG and Opus inference at volume, research-study support, and the engineering and documentation cost of the first external transfer).

### Indicative 24-month envelope

| Phase | Months | Indicative range | Dominant cost |
| --- | --- | --- | --- |
| 0 — Foundation/discovery | 0–3 | $0.3–0.6M | Staff time, governance, IRB |
| 1 — Variant Curator | 3–9 | $0.8–1.5M | Banked-tissue sequencing, ddPCR, integration |
| 2 — Phenome/Trajectory/TAND + surfaces | 9–18 | $1.5–3.0M | SMART-on-FHIR engineering, inference, validation |
| 3 — Therapeutics + research + external | 18–24 | $1.2–2.5M | Opus inference at volume, research, transfer |
| **Total** | **0–24** | **~$3.8–7.6M** | **People, not licenses** |

These are planning ranges, not quotes. The wide bands reflect honest uncertainty about sequencing volume, the scope of Epic app review at the institution, and inference cost at production volume. The substrate is open-source and reused; the spend is overwhelmingly people and the sequencing of real specimens.

### The negative-result publication commitment

One commitment runs through every gate and is non-negotiable: we publish the negative results. If Phase 1's mosaic-recovery rate in genuinely NMI patients turns out to be lower than the somatic-mosaicism literature would predict, that is a finding the field needs. If the TAND agent's surfaced signals prove too noisy to be clinically useful, or if the alert-burden guardrail forces thresholds so high that the async surface adds nothing, those are findings too. If a decision gate stops the Initiative, the reason for stopping is written up.

This is the only posture consistent with the engineering-direct, anti-hype voice the whole project depends on, and the one most likely to earn the trust of a skeptical informatics audience. AI-in-medicine has an abundance of positive demos and a shortage of honest accounts of what did not work. With Dr. Hagedorn's audience, credibility is the entire game, and a published null result is worth more to that credibility than a quietly abandoned pilot.

A final scoping note for honesty: none of this is FDA-cleared, and the software-as-a-medical-device posture for any agent that influences care is undetermined and is itself institutional work, surfacing first in the Phase-1 molecular-genetics review and again at the Phase-2 surface-deployment gate. The synthetic demo required no IRB; every phase here does, and every phase preserves the human gate. The plan is ambitious in scope and conservative in claim, by design.

This is what one institution does with the engine. What follows is how a second, and a tenth, does the same without rebuilding it — the replication path that turns a CCHMC capability into a portable architecture.


---

# Part VI — Beyond TSC


## 22. Beyond TSC: The Replication Thesis

Everything in the preceding sections was built for Tuberous Sclerosis Complex, but very little of it is *about* TSC. That is the point. The TSC Intelligence Engine is the first instance of a pattern, and the pattern, not the instance, is what carries strategic weight. The case for building Engine 7 the way we built it rests on a claim that has to survive scrutiny: the second disease engine costs roughly a third to a half of what the first cost, the third costs less than the second, and the marginal cost continues to fall until what remains is mostly the irreducible biology and the clinician relationships. If that claim holds, then a single institution can stand up a portfolio of rare-disease intelligence engines on the same substrate, and the substrate itself becomes the durable asset. To borrow the mantra that has run through this paper: the institutional pattern is the prize.

This section makes the cost-decay argument concrete, names the second-wave disease candidates and the rationale for each, and describes how the same engine adapts across institutions without a rebuild. It does not relitigate the risks of building on synthetic data or the governance posture; §23 takes those up directly.

### What gets built once

The honest accounting starts by separating what is genuinely TSC-specific from what is shared infrastructure that any disease engine on this platform inherits for free. Six layers were paid for during the eight-week TSC build, and each is reusable without modification or with only configuration-level change.

**Hardware and runtime.** The DGX Spark (GB10 Grace Blackwell, ~1,000 TOPS, 128 GB unified LPDDR5x, 4 TB NVMe) plus the RunPod burst pattern for GPU-accelerated Parabricks variant calling and parallel synthetic-cohort generation is a fixed cost that the first engine amortizes entirely. A second disease engine runs on the same box. There is no incremental hardware spend until concurrent production load forces it, and that is an institutional Phase-1 problem, not a per-disease one.

**Data and synthetic-cohort tooling.** The four-layer synthetic pipeline — Synthea for the demographic and clinical skeleton, BAMSurgeon for the genomic substrate over NA12878-derived BAMs, frontier-model clinical-note and imaging-report generation, all watermarked and version-controlled — is disease-agnostic in its plumbing. Adapting it to a new disease means writing new Synthea disease modules, choosing the variants BAMSurgeon inserts, and curating new note and report templates from the relevant literature. The deterministic ~12-hour regeneration harness, the version-control discipline, the watermarking, and the provenance hooks all transfer unchanged.

**Governance and provenance.** Every agent output already carries model id and version, prompt-template version, retrieved RAG sources with URIs, input hash, and latency, written to an append-only, queryable log. The human-gate discipline (draft-for-review, never autonomous; surfaces patterns, never diagnoses) and the do-not-overclaim posture are encoded in how the surfaces are built, not in TSC content. A new engine inherits the entire compliance story — the part that is hardest to retrofit and most expensive to get wrong — on day one.

**Agent abstractions and orchestration.** This is the largest and least obvious reuse. The five TSC agents are instances of four recurring archetypes:

| Archetype | TSC instance | Reusable across diseases |
|---|---|---|
| Deterministic-pipeline + LLM-synthesis curator | TSC-Variant Curator | Any monogenic/oligogenic disease with a variant-calling + ACMG-AMP path |
| Longitudinal phenotype extractor | TSC-Phenome Mapper | Any disease with HPO-codable, time-anchored clinical features |
| Classical-statistics trajectory modeler | TSC-Trajectory Modeler | Any disease with measurable progressive lesions or quantitative endpoints |
| Discourse-marker surveillance agent | TAND Surveillance Agent | Any under-recognized comorbidity cluster documented in notes |
| Opus-class options-brief synthesizer | TSC-Therapeutics Strategist | Any disease with an evolving therapeutic and trial landscape |

The TSC-Orchestrator — the deterministic LangGraph event router with its 13 event types, event-sourced PostgreSQL state, Redis ephemeral layer, and conservative failure handling — is not rewritten per disease. It is reconfigured: which agents enroll, in what dependency order (the phenotype extractor first remains a near-universal invariant), which events fire. The orchestration *logic* is paid for once.

**RAG substrate.** The Milvus deployment, the BAAI/bge-large-en-v1.5 plus BiomedBERT-derived clinical embedding stack, and the partitioned-corpus design carry over directly. A new disease engine adds a new corpus partition (its consensus guidelines, its mosaicism and natural-history literature, its trial registry slice) rather than standing up new retrieval infrastructure.

**Credibility.** This one resists quantification but is the most valuable. The first engine had to earn the trust of a skeptical clinical-informatics audience from zero. A second engine delivered on the same substrate, with the same provenance and the same honesty about what is and is not built, starts from a track record. Within an institution, credibility compounds.

### The cost-decay claim, stated carefully

Putting numbers on this requires discipline, because overclaiming here would undercut the whole paper. The eight-week TSC build was dominated by work that does not recur: standing up the orchestrator and shared event-sourced state, building the synthetic-cohort pipeline from scratch, and establishing the governance and provenance scaffolding. Call that, roughly, 50–65% of the effort. The genuinely TSC-specific remainder — the disease modules, the ACMG synthesis tuned for TSC1/TSC2, the Marshall-Hagedorn TAND taxonomy, the literature corpus, the eval ground truth — is the part a second disease must redo.

That yields the central estimate: a second disease engine on this substrate is a **~1/3 to ~1/2** effort relative to the first, with the variance driven mostly by how closely the new disease's agent set maps onto the four existing archetypes. A disease that needs all five archetypes and no new ones lands near the bottom of that range. A disease that needs a genuinely new agent type — say, a quantitative imaging-segmentation agent we have not yet built — lands near the top, or above it, because that agent is itself a one-time investment that then becomes reusable. The decay is real but not free, and the first disease that demands a new archetype pays for it on behalf of all the diseases that follow.

### Second-wave candidates

Candidate selection is not arbitrary. The strongest second-wave diseases share three properties with TSC: a tractable genetic or quantitative-phenotype substrate that exercises the existing archetypes, a documented surveillance or under-recognition gap that the engine is uniquely suited to close, and a plausible patient concentration at CCHMC or a partner institution. Five candidates meet the bar, in rough priority order.

1. **Neurofibromatosis type 1 and type 2 (NF1/NF2).** The closest analog to TSC and the obvious second engine. NF1 is a mTOR/RAS-pathway monogenic disorder with a high de-novo rate, somatic mosaicism that standard blood testing misses, multi-organ hamartomatous and tumor burden requiring longitudinal surveillance, and a well-documented neuropsychiatric and learning-disability comorbidity load that is chronically under-addressed. All four archetypes map almost without modification; the TAND surveillance pattern transfers to NF1's cognitive and behavioral profile with the same discourse-marker methodology. NF1 is the lowest-risk, highest-confidence replication target and the one most likely to land near the floor of the cost-decay range.

2. **Rett syndrome.** Predominantly *MECP2*-driven, with a defined natural-history trajectory (regression, motor and respiratory decline) that the classical-statistics Trajectory Modeler is built to forecast, and a rich longitudinal phenotype that the Phenome Mapper handles. The therapeutic landscape is now active enough — including the first approved disease-specific agent — to give the Therapeutics Strategist real material. Rett exercises four archetypes cleanly and adds a natural-history forecasting use case that strengthens the platform's trajectory story.

3. **Williams syndrome.** A contiguous-gene deletion disorder with a distinctive and well-characterized phenotype spanning cardiovascular (supravalvular aortic stenosis), cognitive, and behavioral domains. Williams tests the engine's range: its substrate is a structural deletion rather than point mutation, which stresses the Variant Curator's archetype usefully, and its behavioral profile is a natural fit for the surveillance agent. It is the candidate most likely to surface a needed new capability, which makes it valuable precisely for stretching the platform.

4. **Cornelia de Lange syndrome.** Genetically heterogeneous (multiple cohesin-pathway genes), with significant mosaicism that, as with TSC, blood testing under-detects — directly exercising the mosaic-aware low-VAF recovery that is the Variant Curator's signature capability. The multi-system phenotype and the under-recognized neurobehavioral burden round out a strong archetype fit. Cornelia de Lange is the candidate that most reinforces the "deep-tissue mosaic recovery" thesis beyond TSC.

5. **The broader mTORopathies.** Beyond TSC, the family of mTOR-pathway disorders — focal cortical dysplasia, hemimegalencephaly, *PIK3CA*-related and *DEPDC5*-related epilepsies, the broader spectrum of mosaic mTOR activation in epileptic brain — shares both the molecular biology and the somatic-mosaicism diagnostic challenge. These are not a single disease but a coherent cluster, and they let the engine reuse not just its archetypes but a substantial fraction of its actual TSC content: the mTOR-pathway corpus, the mosaic-recovery tuning, the everolimus/sirolimus therapeutic knowledge. The mTORopathies are where the cost-decay curve is steepest, because the *biology* itself is largely shared, not merely the infrastructure.

The pattern across all five is deliberate: each one either confirms an existing strength (NF1, Cornelia de Lange, the mTORopathies) or productively stretches the platform into a new capability that then becomes reusable (Rett's natural-history forecasting, Williams's structural substrate). A portfolio assembled in roughly this order spends its early replication budget on near-certain wins and its later budget on capability expansion.

### Cross-institution adaptation: swap the box labels, keep the wiring

The replication thesis has a second axis. The same engine that replicates across diseases within CCHMC also replicates across institutions for the same disease, and the architecture was drawn from the start to make that move cheap.

Recall the CCHMC source-area architecture from earlier sections: the Winslow Research Pavilion is the physical envelope, and five institutional areas are *sources* that feed the engine rather than parts of it — the Discover Together Biobank's banked tuber, AML, and SEGA tissue feeding the Variant Curator; Biomedical Informatics and the Marshall-Hagedorn methodology feeding the TAND agent; the TSC clinical and research program supplying patient concentration; the Epic Clarity/Caboodle and biobank LIMS plumbing feeding the Phenome Mapper, TAND, and Trajectory agents. The engine sits downstream of these sources and upstream of the clinician surfaces. Crucially, the engine does not care *which* institution's sources it is wired to, only that the source *types* are present.

That is what "swap the box labels, keep the wiring" means in practice. Adapting the TSC engine to TGen or City of Hope is a matter of repointing each source connection: their biospecimen archive in place of the Discover Together Biobank, their informatics group's documentation conventions in place of CCHMC's, their EHR's clinical-data export in place of Clarity/Caboodle, their TSC or rare-disease program supplying the patient cohort. The agents, the orchestrator, the RAG substrate, the surfaces, and the governance scaffolding are untouched. What changes is configuration and connector code at the source boundary — by design the thinnest, most replaceable layer in the system.

This is also the honest place to restate a hard rule from earlier sections. The Epic/Clarity, LIMS, and biobank integrations are described here architecturally; they are *not built* in the synthetic-data demo and are institutional Phase-1 work wherever the engine lands. The cross-institution claim is that the *integration surface is small and well-defined*, not that the integration is done. A biobank without an intelligence layer is a freezer full of tubes; the engine is the kitchen, not the meal — and the kitchen has been designed so that a different institution can bring its own ingredients.

### Why the pattern is the prize

Put the two axes together. A single substrate, paid for once, supports a portfolio of disease engines (NF1/NF2, Rett, Williams, Cornelia de Lange, the mTORopathies) each at a fraction of the first engine's cost, and each of those engines can in turn be re-sourced to a partner institution at the cost of connector code. The compounding is multiplicative: every new disease makes the next disease cheaper, and every new institutional deployment makes the next deployment cheaper, while credibility accrues across both axes.

TSC was the right first instance not because it is the largest market — it is a 1-in-6,000 rare disease — but because it exercises every archetype the platform needs, sits inside CCHMC's genuine clinical and biobank strengths, and extends a CCHMC investigator's own published methodology rather than grafting something foreign onto the institution. It is, in the fullest sense, a proof case. The argument of this paper is that the proof, once delivered, is general. What CCHMC would be standing up is not a TSC tool but a rare-disease intelligence platform that happens to start with TSC.

With the upside stated plainly, the obligation is to be equally plain about what could go wrong. §23 turns to the risks — technical, clinical, institutional, and ethical — that this replication thesis has to survive.


## 23. Risks, Limitations & Open Questions

A system whose only job is to earn the trust of skeptical clinicians and informaticians cannot survive a dishonest risk section. The credibility of the entire engagement rests on naming the things that can go wrong before anyone in the room does. This section is therefore deliberately uncomfortable. It separates two categories that are too often conflated: the risks of *building and demonstrating* the engine on the eight-week timeline, and the risks of *operationalizing* it inside Cincinnati Children's Hospital Medical Center (CCHMC) over the institutional Phase-1 horizon. It closes by stating the demonstration's inherent limitations plainly, so that no claim in the preceding twenty-two sections can be mistaken for more than it is.

The governing posture throughout: the synthetic-data demonstration is what runs now on the DGX Spark with RunPod overflow. Everything touching real Epic Clarity/Caboodle data, the Discover Together Biobank LIMS, real imaging, and any FDA or Software-as-a-Medical-Device (SaMD) posture is described architecturally and is explicitly *not built*. Treating that boundary as fuzzy is itself the largest risk to credibility, and we draw it hard.

### 23.1 Build risks (the eight-week MVP)

These are the risks to delivering a working, defensible demo for early Q3 2026.

**Schedule compression.** Eight weeks for a five-agent engine plus a deterministic orchestrator plus a 50-patient synthetic cohort plus three clinician surfaces is aggressive. The single largest schedule dependency is the cohort itself (W1–W2): every downstream agent eval depends on having version-controlled ground truth. If the four-layer cohort pipeline (Synthea → BAMSurgeon → frontier-model notes → frontier-model imaging reports) slips, the slip cascades. *Mitigations:* the cohort is generated once and version-controlled, so it is built before any agent is wired against it; the explicit de-scope order is pre-committed (cohort 50→30; simplify the async alert surface; simplify the TAND cluster set; and only as a last resort cut the Therapeutics brief from six sections to four). The de-scope order protects the two load-bearing demo moments — Act One mosaic recovery on Patient A and Act Two's Patient B dashboard — which never get cut.

**API cost.** The Therapeutics Strategist is Opus-class and non-negotiable; the Variant Curator's ACMG synthesis and the Phenome Mapper's rare-conflict resolution also reach for Opus. Across a 50-patient cohort with ~600–1000 generated notes, repeated eval runs and cohort regenerations (~12 hr each) can accumulate real spend. *Mitigations:* tiered model routing is enforced (Haiku for ICD-10/lab→HPO normalization and prose summaries, Sonnet for per-note extraction, Opus reserved for synthesis and conflict resolution); Milvus RAG retrieval is cached per corpus partition; the local Llama 3.1 70B Instruct fallback via Ollama absorbs development iteration that does not need frontier quality. Cost is a tracked line item, not a surprise.

**Synthetic realism.** This is the build risk most likely to be challenged in the room, and it deserves the most candor. BAMSurgeon variants spiked into NA12878-derived BAMs at mosaic VAF 4–12% and germline ~50% produce *plausible* reads, but spike-in artifacts (strand bias, edge effects near insertion sites, unrealistic error profiles) can either flatter the Variant Curator (variants too easy to call) or mislead it (artifacts the mosaic-aware caller learns to over-trust). Frontier-model clinical notes can drift toward template uniformity, under-representing the messy hedging and omission that the TAND Surveillance Agent specifically exists to detect — a circularity risk, because the same class of model that writes the notes also reads them. *Mitigations:* clinician-sampled review of generated notes before they enter the cohort; notes built from *published* TSC documentation templates rather than free generation; persistent watermarking so synthetic provenance is never lost; and an eval design that scores the Variant Curator against the *known* spiked VAFs and the TAND agent against *deliberately embedded* signals, with explicit attention to false positives (no FP Pathogenic calls; no spurious TAND flags). We state openly that strong eval numbers on this cohort are evidence of correct *engineering*, not of *clinical* accuracy.

**Agent quality and the human gate.** Each agent can fail in a domain-specific way: the Variant Curator can mis-synthesize an ACMG-AMP classification relative to the combinatorial rules it is validated against; the Phenome Mapper can miss a phenotype (recall target ≥90%) or hallucinate one (precision target ≥85%); the Trajectory Modeler — classical statistics, not LLM — can produce overconfident prediction intervals if the synthetic longitudinal series is too clean; the TAND agent can over-flag; the Therapeutics Strategist can under-hedge or drop a source attribution. *Mitigations:* every agent output is draft-for-review with a named human gate (the Variant Curator's molecular-genetics report is a draft for a board-certified molecular geneticist's sign-off, never autonomous); the Trajectory Modeler's 50%/90% prediction intervals are reported as intervals, not point estimates; full provenance (model id/version, prompt-template version, retrieved RAG URIs, input hash, latency) accompanies every output so any error is traceable rather than mysterious.

**Technical drift.** Frontier model versions change under us; prompt templates evolve; the Milvus TSC corpus partition can go stale against ClinicalTrials.gov and FDA snapshots. A demo that worked in W6 can degrade by W8. *Mitigations:* model ids are pinned and recorded in provenance; prompt templates are versioned; the cohort and corpus snapshot are frozen for the demo build; W8 is reserved for a full regen, dry runs, and clinician review precisely so that drift is caught before delivery, not during it.

**Demo logistics.** A 30-minute live demonstration with Act One running mosaic recovery live on Patient A carries the ordinary risks of any live system: network dependence for API calls, latency spikes, a cold cache, a RunPod instance that does not spin up. *Mitigations:* the eval latency targets (Variant Curator <5 min/case; Therapeutics <3 min) leave headroom; a pre-warmed, recorded fallback of each act exists; the open audit trail shown in Act One is a static, queryable artifact that does not depend on a live call to render.

### 23.2 Institutional and long-term risks

These risks live beyond the demo, in the Phase-1 institutional work, and they are not ours alone to mitigate — most require CCHMC partnership. We name them honestly rather than minimize them.

**Adoption.** The hardest problem is not technical. Clinicians do not lack for screens, and alert fatigue is real and well-documented in this very institution's informatics literature (Orenstein et al. 2021, JAMIA). An engine that adds another interruptive surface will be ignored or resented. *Mitigations baked into the design:* the TAND output is pre-visit *briefing material*, never an interruptive alert and never a diagnosis; the async alert surface enforces a strict discipline and recalibrates if it exceeds ~3 alerts per clinician per week; surfaces are demand-driven and progressively disclosed. Whether this is sufficient is an open question that only real clinical use can answer.

**Governance.** Real-data Phase-1 work requires IRB approval (the synthetic demo does not). It also requires answering who owns agent outputs, how the append-only event log is retained and audited, how model updates are governed against drift, and how a draft molecular-genetics report's provenance survives in the medical record. These are institutional decisions, not engineering defaults.

**Sustainability.** The platform is Apache 2.0 and non-commercial in intent. That is a strength for trust and a risk for longevity: a system maintained by a single architect on a single Spark is not an institutional dependency anyone should take lightly. The replication thesis — "swap the box labels, keep the wiring" to TGen or City of Hope — only matters if the wiring is maintained. *Mitigation posture:* the reuse of the existing HCLS AI Factory v1.3.0 substrate (LangGraph, Milvus, PostgreSQL/Redis/MinIO, tiered models) means TSC is net-new orchestration and agents on a maintained base, not a bespoke stack; the path to institutional ownership is itself Phase-1 scope.

**FDA / regulatory shift.** The SaMD posture is undetermined and is explicitly institutional work, not a claim we make. The regulatory landscape for clinical-decision-support and for LLM-in-the-loop tools is moving; a surface that is decision-support today could be re-scoped as a regulated device under a future interpretation, particularly anything touching the Variant Curator's molecular-genetics draft. *Mitigation posture:* the human gate on every output and the decision-support (not recommendation) framing of the Therapeutics brief are deliberate regulatory-conservative choices, but they do not pre-empt a regulatory determination, and we do not pretend they do.

### 23.3 The demonstration's inherent limitations

These are not failures to be fixed; they are the *scope* of the demonstration, stated so that the eval results in §22 are read correctly.

- **It is entirely synthetic.** Fifty patients, deterministically generated, persistently watermarked. No real patient data of any kind is present. Eval numbers measure engineering correctness against known ground truth, not clinical concordance.
- **There is no real imaging.** The cohort starts at BAM, not FASTQ; it contains frontier-model imaging *reports*, not DICOM. The longitudinal SEGA series that drives Patient B's trajectory forecast is a generated narrative, not pixels. No imaging-AI pipeline is built.
- **There is no Epic, Clarity, Caboodle, or LIMS integration.** The Phenome Mapper, TAND, and Trajectory agents consume Clarity-*shaped* synthetic relational data. The actual data plumbing into CCHMC systems is described architecturally and is Phase-1 institutional work, not demonstrated.
- **The mosaic-recovery result is not clinically validated.** Recovering all seven mosaic variants at VAF ≥5% on this cohort, with correct ACMG class and no false-positive Pathogenic calls, demonstrates that the mosaic-aware pipeline does what it is designed to do on data where we planted the answer. It does not establish performance on real tuber, AML, or SEGA tissue from the Biobank — that concordance study is the central Phase-1 scientific question.
- **The CCHMC source areas are sources, not built integrations.** The Discover Together Biobank, Biomedical Informatics, the TSC clinical program, and the Clarity/LIMS plumbing feed the engine in the architecture; in the demo they are represented by synthetic stand-ins.

### 23.4 Open questions

Three questions remain genuinely open, and we would rather pose them than paper over them. First, does mosaic-aware recovery hold on real biobanked tissue — the question that justifies the whole Variant Curator? Second, does the Marshall-Hagedorn discourse-marker taxonomy, extended here to the six TAND clusters, transfer from the inpatient diagnostic-uncertainty corpus on which it was developed to longitudinal outpatient TSC notes without losing precision? Third, can briefing-style surfacing actually change surveillance behavior, or does it join the long list of well-intentioned tools that clinicians route around? None of these is answerable by the demonstration. All three are answerable by the Phase-1 partnership the demonstration is built to earn — which is the subject of the conclusion that follows.


## 24. Conclusion & Call to Action

The preceding twenty-three sections described, in deliberate detail, a system that does not yet exist as deployed clinical software. That is the point. This paper is a build specification written as an argument, and the argument is that the build is the only honest form the argument can take. So before the invitation, a synthesis of what has actually been claimed, what has been demonstrated, and what is being asked.

### What the TSC Intelligence Engine is

The TSC Intelligence Engine is Engine 7 of the open-source HCLS AI Factory: five coordinated agents and a deterministic orchestrator that turn the genomic and clinical data already accumulating around a Tuberous Sclerosis Complex patient into something a clinician can use in the eight minutes before a visit. It is an engine rather than a single agent because it spans modalities and coordinates: the TSC-Variant Curator recovers low-VAF somatic mosaic variants from tissue that standard blood testing misses; the TSC-Phenome Mapper builds the longitudinal HPO-coded substrate the others stand on; the TSC-Trajectory Modeler forecasts SEGA and AML growth and seizure burden with classical statistics and honest prediction intervals; the TAND Surveillance Agent surfaces the neuropsychiatric features that TOSCA tells us are missed or unaddressed in 30 to 50 percent of patients; and the TSC-Therapeutics Strategist assembles a source-attributed, six-section options brief. The TSC-Orchestrator is a LangGraph event router with an append-only event log, not a language model — the connective tissue is deterministic by design, because a clinician auditing why a surface said what it said deserves a trace, not a confabulation.

None of these agents diagnoses. Each produces decision-support with a human gate. The Variant Curator output is a draft for a board-certified molecular geneticist to sign; the TAND agent surfaces patterns as pre-visit briefing material, never an interruptive alert and never a label. Every output carries provenance: model id and version, prompt template version, retrieved RAG sources with URIs, input hash, latency. This is not decoration. With the audience this paper is written for, the provenance trail is the difference between a tool and a liability.

### What it demonstrates now, and what it explicitly does not

What runs now, on a single NVIDIA DGX Spark with RunPod GPUs attached for the heavy variant-calling and cohort-generation passes, is a synthetic-data demonstration. Fifty version-controlled synthetic patients — thirty TSC2 germline, twelve TSC1 germline, seven mosaic, one NMI — built from Synthea clinical skeletons, BAMSurgeon-seeded genomic substrate, and frontier-model notes and imaging reports, every artifact persistently watermarked as synthetic. Against that known ground truth the demo targets are concrete and falsifiable: recover all seven mosaic variants at VAF greater than or equal to five percent with correct ACMG classification and no false-positive Pathogenic calls in under five minutes per case; Phenome Mapper recall of at least ninety percent at precision at least eighty-five; a Trajectory forecast that puts Patient B's SEGA crossing its threshold in the 12-to-18-month window without false alarms; a Therapeutics brief in under three minutes with appropriate hedging and full attribution. The three-act demo — mosaic recovery live on Patient A with the audit trail open, the four-quadrant dashboard and pre-visit briefing on Patients B and C, then the infrastructure and cost story — is thirty minutes plus discussion, targeted for early Q3 2026.

What does not run now must be said as plainly as what does. The Epic Clarity and Caboodle integration, the biobank LIMS plumbing, and the imaging-AI pipelines are described architecturally and are explicitly not built. They are institutional Phase-1 work. The engine is not FDA-cleared; its software-as-a-medical-device posture is undetermined and is itself institutional work. No IRB is required for a synthetic demonstration; real-data Phase-1 will require one. These are not caveats buried in a footnote. They are the boundary of the claim, and keeping that boundary visible is what makes the rest credible.

### What the Cincinnati Children's engagement could become

The reason this engine is scoped to TSC and aimed at Cincinnati Children's is not opportunism. It is that the pieces already fit. The Discover Together Biobank holds banked tuber, AML, and SEGA tissue — the molecular substrate the Variant Curator needs to do the one thing a blood draw cannot, which is recover the mosaic variant in the roughly ten to fifteen percent NMI cohort. A biobank without an intelligence layer is a freezer full of tubes; the engine is the layer. The Division of Biomedical Informatics under Dr. Philip Hagedorn has already published the diagnostic-uncertainty discourse methodology (Marshall et al. 2023; Nickels et al. 2024) that the TAND Surveillance Agent applies; the agent is an extension of that work, not an external graft onto it. The TSC clinical and research program concentrates the patients. The Winslow Research Pavilion is the envelope around all of it. Five institutional areas feed the engine; the engine feeds three clinician surfaces; the wiring does not change when the box labels do.

That last point is the one that scales. The same architecture — swap TSC1/TSC2 for NF1, NF2, the broader mTORopathies, Rett, Williams — replicates wherever a biobank, an informatics division, and a disease program coexist. The Phase-1 institutional plan turns the synthetic demonstration into a real-data deployment behind the firewall, on real Clarity exports and real banked tissue, under IRB, with the molecular-genetics and TAND outputs validated against actual clinician judgment rather than synthetic ground truth. None of that is promised here as accomplished. It is the work the demonstration is meant to earn the right to start.

### Why the cumulative effect is the point

A single recovered mosaic variant changes one family's answer from "we don't know why" to a Likely Pathogenic call with a ddPCR confirmation path. A single surfaced TAND signal moves one psychosocial concern from the margin of a note into a pre-visit briefing where it can be addressed. Neither is revolutionary in isolation. The argument is in the accumulation: across the patients in one TSC program over one year, across the diseases that share the wiring, across the institutions that share the pattern, a small consistent recovery of signal that is currently lost compounds. The engine does not try to be smarter than the clinician. It tries to make sure the clinician is never working from less than the data already contains. That is a modest claim per patient and a large one in aggregate, and it is the only kind of claim a synthetic demonstration can responsibly set up.

### The invitation

Dr. Hagedorn has offered to convene his Biomedical Informatics team and the faculty TSC lead. This is the concrete next step, and it is specific: a working session at Cincinnati Children's around the live demonstration on the Spark, run end to end on the fifty-patient synthetic cohort, with the audit trail open for inspection. The agenda is three questions. First, do the Variant Curator and TAND outputs meet the standard of the people who would have to sign or act on them — and where do they not. Second, which of the five institutional source areas — biobank tissue, Clarity and Caboodle plumbing, the clinical program, the informatics methodology — is the right first real-data integration, and what does its Phase-1 scope and IRB pathway look like. Third, what would constitute a genuine validation, against real clinician judgment, that this engine recovers signal that is otherwise lost.

The code is Apache 2.0 and open. The substrate is reused from HCLS AI Factory v1.3.0; the net-new work is the orchestrator, the synthetic cohort, and the five agents. The demonstration is built to be inspected, not believed. The appendices that follow give the schemas, prompts, eval harnesses, and cohort specifications in the detail required to reproduce and to interrogate every claim above. The invitation is to do exactly that.


---

## Appendices

These appendices collect the reference material the body points to but does not interrupt itself to enumerate: a shared vocabulary, the concrete data contracts between agents, the exact synthetic-cohort specification, the literature that anchors the clinical and methodological claims, and a crosswalk to the three prior TSC documents and the FAQ. Nothing here is new design; it is the design made checkable. Where a schema or a parameter appears below, it is the same one the build targets — the appendices exist so that a skeptical reader can verify the engine against itself.

### Appendix A: Glossary

Acronyms are defined once here and used without expansion elsewhere in the paper. Terms are grouped by domain.

#### A.1 Biology and genetics

| Term | Definition |
|---|---|
| **TSC1 / TSC2** | The two genes whose loss-of-function causes Tuberous Sclerosis Complex. Their protein products (hamartin / tuberin, with TBC1D7) form a complex that inhibits mTOR. A pathogenic variant in either releases the brake on cell growth. |
| **mTOR** | Mechanistic target of rapamycin; the growth-signaling kinase the TSC complex normally restrains. The therapeutic rationale for everolimus and sirolimus. |
| **mTORC1** | The rapamycin-sensitive mTOR complex; target of next-generation selective inhibitors now in development. |
| **Hamartoma** | A benign, disorganized overgrowth of cells native to the affected organ — the unifying lesion type in TSC. |
| **Mosaicism** | The presence of two or more genetically distinct cell populations arising from a single zygote. In TSC, a variant present only in a fraction of cells, often undetectable in blood but recoverable in affected tissue. |
| **VAF** | Variant allele fraction; the proportion of sequencing reads supporting the variant allele. Germline heterozygous variants sit near 50%; mosaic variants present at low VAF (the engine targets recovery at VAF ≥ 5%). |
| **NMI** | No mutation identified; the ~10–15% of clinically diagnosed TSC patients with negative conventional (blood) genetic testing, usually explained by somatic mosaicism. |
| **De novo** | A variant arising newly in the patient, not inherited (~2/3 of TSC cases). |
| **ddPCR** | Droplet digital PCR; an orthogonal, high-sensitivity assay used to confirm low-VAF mosaic calls. The Variant Curator recommends it; it does not perform it. |
| **ACMG-AMP** | The American College of Medical Genetics and Genomics / Association for Molecular Pathology variant classification framework (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign) built from weighted evidence criteria (e.g., PVS1, PM2, PP4). |
| **ClinVar / gnomAD / LOVD / dbSNP / HPO / SNOMED-CT** | Reference resources: clinical variant interpretations (ClinVar); population allele frequencies (gnomAD v4); locus-specific TSC variant database (LOVD-TSC); variant identifiers (dbSNP); phenotype ontology (Human Phenotype Ontology); clinical terminology (SNOMED-CT). |

#### A.2 Disease and clinical

| Term | Definition |
|---|---|
| **SEGA** | Subependymal giant cell astrocytoma; a periventricular tumor (often at the foramen of Monro) that can obstruct CSF flow and cause hydrocephalus. Growth trajectory is a primary surveillance concern. |
| **AML** | Angiomyolipoma; benign renal lesions present in ~80% of patients, with bleeding risk rising above ~4 cm. |
| **LAM** | Lymphangioleiomyomatosis; cystic lung disease seen in adult women with TSC. |
| **Cortical tuber** | Focal cortical malformation; the brain hamartoma most associated with epilepsy. |
| **Rhabdomyoma** | Cardiac hamartoma, frequently present at birth and usually regressing. |
| **TAND** | TSC-Associated Neuropsychiatric Disorders; the behavioral, psychiatric, intellectual, academic, neuropsychological, and psychosocial manifestations affecting ~90% of patients, with 30–50% of features unaddressed (TOSCA). |
| **TAND-L** | The TAND Lifetime Checklist; the consensus instrument for structured TAND screening. |
| **PKD1** | Polycystic kidney disease gene; a TSC2/PKD1 contiguous gene deletion (~1–2% of cases) produces severe early polycystic kidney disease. |
| **ITSC** | International Tuberous Sclerosis Complex consensus group; author of the 2021 surveillance guidelines the engine encodes. |
| **TOSCA** | TuberOus SClerosis registry to increase disease Awareness; the international natural-history registry. |
| **EXIST-3** | The trial establishing adjunctive everolimus for TSC-associated refractory seizures; cited as the first targeted therapy in a genetically defined epilepsy. |

#### A.3 AI, infrastructure, and methodology

| Term | Definition |
|---|---|
| **The Engine** | The TSC Intelligence Engine; Engine 7 of the HCLS AI Factory. Five agents plus one deterministic orchestrator plus three clinician surfaces. |
| **Agent** | A bounded analytic unit (Variant Curator, Phenome Mapper, Trajectory Modeler, TAND Surveillance, Therapeutics Strategist) with defined inputs, outputs, and model tiers. |
| **Orchestrator** | The deterministic LangGraph event router (not an LLM) that sequences agents and assembles surfaces over event-sourced state. |
| **Surface** | A clinician-facing product: pre-visit briefing, in-visit dashboard, async alert surface. Standalone web apps, never Epic, persistently watermarked synthetic. |
| **Model tiers** | Claude **Haiku** (fast normalization/summary), **Sonnet** (per-note analysis, evidence aggregation), **Opus** (high-stakes synthesis: ACMG classification, therapeutics brief), with local **Llama 3.1 70B Instruct** via Ollama as fallback. |
| **RAG** | Retrieval-augmented generation; here a Milvus vector store (BAAI/bge-large-en-v1.5 plus BiomedBERT-derived clinical embeddings) over a TSC literature corpus partition. |
| **Provenance** | The record attached to every agent output: model id/version, prompt template version, retrieved sources with URIs, input hash, latency. Append-only and queryable. |
| **Event sourcing** | The state model: a PostgreSQL append-only event log with materialized current-state projections; Redis for ephemeral TTL state. |
| **Marshall-Hagedorn taxonomy** | The diagnostic-uncertainty discourse-marker framework (hedging, deferral, third-party attribution, conditional, follow-up-without-formalization) the TAND agent applies; see Appendix D. |
| **DGX Spark** | The first build target: NVIDIA GB10 Grace Blackwell, ~1,000 TOPS, 128 GB unified LPDDR5x, 4 TB NVMe, DGX OS. **RunPod** supplies burst GPUs (Parabricks variant calling, cohort generation, heavier local inference). |
| **Synthea / BAMSurgeon / snpEff / VEP** | Cohort tooling: clinical-record simulator (Synthea); spike-in variant editor for BAMs (BAMSurgeon); variant annotators (snpEff, VEP). |

#### A.4 Institutional

| Term | Definition |
|---|---|
| **CCHMC** | Cincinnati Children's Hospital Medical Center; the beachhead institution. |
| **Winslow Research Pavilion** | The ~45,000 sq ft Avondale research building (opened 2024) that houses the Discover Together Biobank and the Gamble Vaccine Research Center — the institution's physical infrastructure envelope. |
| **Discover Together Biobank** | CCHMC's biospecimen archive; the source of banked tuber/AML/SEGA tissue that is the molecular substrate for mosaic recovery. |
| **Epic Clarity / Caboodle** | The reporting databases behind the Epic EHR; the data plumbing for Phenome Mapper, TAND, and Trajectory in institutional Phase-1 work (explicitly **not built** in the synthetic demo). |
| **CHIO** | Chief Health Informatics Officer; Dr. Philip A. Hagedorn's role at CCHMC. |

### Appendix B: Agent input/output JSON schemas

These are the data contracts between agents and the orchestrator. Fields shown are illustrative of the build, not exhaustive; every output object carries a shared `provenance` block (omitted from each schema after the first for brevity). All outputs are decision-support drafts for human review.

#### B.1 Shared provenance block

```json
{
  "provenance": {
    "agent": "tsc-variant-curator",
    "agent_version": "0.4.0",
    "model_calls": [
      {"role": "acmg_synthesis", "model": "claude-opus-4-8", "prompt_template": "acmg_synth_v3"}
    ],
    "rag_sources": [
      {"uri": "pmc:PMC4537774", "title": "Tyburczy 2015", "score": 0.88}
    ],
    "input_hash": "sha256:9f2c...e1",
    "latency_ms": 41280,
    "generated_at": "2026-07-14T16:22:09Z",
    "human_review_status": "pending"
  }
}
```

#### B.2 TSC-Variant Curator

```json
// input
{
  "patient_id": "SYNTH-A-001",
  "specimen": {"type": "tuber_tissue", "source": "biobank_synthetic"},
  "bam_uri": "minio://cohort/SYNTH-A-001/tuber.bam",
  "reference": "GRCh38",
  "mosaic_aware": true,
  "min_vaf": 0.05
}
// output
{
  "patient_id": "SYNTH-A-001",
  "variants": [
    {
      "gene": "TSC2",
      "hgvs_c": "c.4639del",
      "hgvs_p": "p.(Ser1547fs)",
      "consequence": "frameshift",
      "vaf": 0.083,
      "supporting_reads": 47,
      "total_depth": 566,
      "strand_balance": 0.52,
      "mosaic_flag": true,
      "artifact_score": 0.04,
      "acmg_classification": "Likely Pathogenic",
      "acmg_criteria": ["PVS1", "PM2", "PP4"],
      "classification_basis": "frameshift in TSC2 + absent gnomAD + phenotype match",
      "ddpcr_recommended": true,
      "clinvar_spec_interpretation": "..."
    }
  ],
  "draft_report_uri": "minio://reports/SYNTH-A-001/molgen_draft.pdf",
  "ai_labeled": true,
  "review_required_by": "board_certified_molecular_geneticist"
}
```

#### B.3 TSC-Phenome Mapper

```json
// output (per patient)
{
  "patient_id": "SYNTH-B-002",
  "hpo_profile": [
    {
      "hpo_id": "HP:0002133",
      "label": "Status epilepticus",
      "onset_date": "2019-03-02",
      "evidence_spans": [
        {"note_id": "N-1182", "char_start": 410, "char_end": 478, "source_model": "claude-sonnet-4-5"}
      ],
      "normalization_path": "icd10:G40.901 -> HP:0002133",
      "confidence": "high"
    }
  ],
  "discordance_log": [
    {"hpo_id": "HP:0009721", "note_id": "N-1207", "issue": "conflicting laterality", "resolved_by": "claude-opus-4-8"}
  ],
  "surveillance_gaps": [
    {"itsc_recommendation": "renal MRI q1-3y", "last_observed": "2022-08", "gap_status": "overdue"}
  ]
}
```

#### B.4 TSC-Trajectory Modeler

Classical statistics only; no LLM in the numeric path. Haiku produces the prose `summary`; Sonnet is invoked only when `flag == "unusual_trajectory"`.

```json
// output
{
  "patient_id": "SYNTH-B-002",
  "forecasts": [
    {
      "metric": "sega_diameter_cm",
      "method": "gaussian_process_regression",
      "observed": [{"date": "2024-09", "value": 0.8}, {"date": "2025-09", "value": 1.1}, {"date": "2026-03", "value": 1.3}],
      "horizons": [
        {"months": 6,  "point": 1.42, "pi50": [1.36, 1.49], "pi90": [1.22, 1.63]},
        {"months": 12, "point": 1.58, "pi50": [1.47, 1.70], "pi90": [1.25, 1.92]},
        {"months": 18, "point": 1.74, "pi50": [1.58, 1.91], "pi90": [1.28, 2.21]}
      ],
      "threshold_crossing": {"threshold_cm": 1.5, "window": "12-18mo", "probability": 0.71}
    }
  ],
  "surveillance_recommendation": "shorten brain MRI cadence to q6mo",
  "flag": "threshold_watch"
}
```

#### B.5 TAND Surveillance Agent

Output is briefing material, never an alert and never a diagnosis. The numeric `cluster_scores` come from the deterministic aggregation layer over Sonnet-extracted discourse markers; Opus writes only the `briefing_summary`.

```json
// output
{
  "patient_id": "SYNTH-B-002",
  "tand_clusters": [
    {
      "cluster": "psychiatric",
      "signal_strength": "emerging",
      "markers": [
        {"type": "hedging", "note_id": "N-1340", "span": "parents wonder if anxiety is increasing"},
        {"type": "follow_up_without_formalization", "note_id": "N-1377", "span": "will revisit at next visit"}
      ],
      "score": 0.61
    }
  ],
  "tand_l_coverage": {"clusters_documented": 3, "clusters_total": 6},
  "briefing_summary": "Under-recognized signals in 2 clusters; consider TAND-L at next visit.",
  "presentation": "pre_visit_briefing",
  "is_diagnosis": false
}
```

#### B.6 TSC-Therapeutics Strategist

Opus-class, non-negotiable. The six sections are fixed; every claim carries a source attribution and the brief is framed as decision-support.

```json
// output
{
  "patient_id": "SYNTH-C-003",
  "options_brief": {
    "current_therapy": {"text": "...", "sources": ["pmc:..."]},
    "optimization": {"text": "...", "sources": ["fda:..."]},
    "combination": {"text": "...", "sources": ["pmc:..."]},
    "trial_matching": [
      {"nct_id": "NCT0XXXXXXX", "title": "...", "match_basis": "TSC1 + refractory focal seizures", "eligibility_flags": ["age_ok", "prior_mTORi_review_needed"]}
    ],
    "emerging_evidence": {"text": "...", "sources": ["pmc:..."]},
    "open_questions": ["mucositis recurrence risk on dose re-escalation?"]
  },
  "framing": "decision_support_not_recommendation",
  "all_claims_attributed": true
}
```

### Appendix C: Synthetic cohort specification

The cohort is generated once, version-controlled, deterministic (fixed seeds), and regenerable in roughly 12 hours. It is the ground truth against which the demo evaluations are scored; it is not clinical data and is persistently watermarked synthetic.

#### C.1 Composition (n = 50)

| Stratum | Count | Share | Notes |
|---|---|---|---|
| TSC2 germline | 30 | 60% | Pathogenic/likely-pathogenic at germline VAF (~50%) |
| TSC1 germline | 12 | 24% | |
| TSC2 mosaic | 5 | 10% | Low-VAF spike-in (4–12%) |
| TSC1 mosaic | 2 | 4% | Low-VAF spike-in (4–12%) |
| NMI | 1 | 2% | Blood-negative; recoverable only in tissue |
| **Mosaic total** | **7** | **14%** | The mosaic-recovery evaluation target |

#### C.2 Four-layer generation pipeline

1. **Clinical skeleton** — Synthea (MIT) with custom TSC modules emits demographics, encounters, problems, medications, and labs as FHIR R4 and Clarity-shaped relational tables.
2. **Genomic substrate** — BAMSurgeon inserts TSC1/TSC2 variants into NA12878-derived BAMs (germline ~50% VAF; mosaic 4–12% VAF), followed by Parabricks-equivalent calling (BWA-MEM/GATK HaplotypeCaller + Mutect2) to produce realistic VCFs. The cohort starts at BAM; there is no raw FASTQ.
3. **Clinical notes** — ~600–1000 frontier-model notes generated from published templates, clinician-sampled for review, watermarked.
4. **Imaging reports** — frontier-model brain MRI / renal ultrasound / echo / ophthalmology reports, including longitudinal SEGA series at ~2–4 mm/yr. These are **text reports only** — no DICOM, no pixels.

The cohort explicitly **does not contain**: real imaging or DICOM, raw FASTQ, neuropsychological test scores, pedigree data beyond what notes imply, pharmacy/claims data, or patient-reported-outcome scores.

#### C.3 Featured patients

| | Patient A | Patient B | Patient C |
|---|---|---|---|
| Profile | 4yo F | 12yo M | 18yo F |
| Genetics | NMI; TSC2 frameshift at 8.3% VAF in tuber tissue | TSC2 c.3037C>T p.Arg1013Ter | TSC1 |
| Result | Likely Pathogenic (PVS1+PM2+PP4) + ddPCR rec | SEGA 0.8→1.1→1.3 cm at foramen of Monro; bilateral AML 2.8 cm; well-controlled focal seizures; scattered under-recognized TAND | partial everolimus response, mucositis-driven dose reduction; AML ~4 cm; refractory focal seizures |
| Demo role | Act One: live mosaic recovery + open audit trail | Act Two: 4-quadrant dashboard, SEGA trajectory, TAND briefing | Act Two: therapeutics brief |

### Appendix D: Selected references

Cited by author/year inline throughout the paper; resolved here. This is a working anchor list, not an exhaustive bibliography.

**Surveillance and natural history**
- International Tuberous Sclerosis Complex Consensus Group (Krueger, Northrup, et al.), 2021 — updated diagnostic criteria and surveillance/management recommendations (ITSC 2021).
- Kingswood et al., TOSCA — the international TSC natural-history registry; source for the 30–50% TAND under-recognition figure.

**Somatic mosaicism and variant detection**
- Tyburczy et al., 2015 — mosaic and intronic TSC mutations in the "no mutation identified" cohort; the empirical basis for tissue-level deep sequencing.
- Giannikou et al., 2016 — whole-exome characterization of TSC2 mosaicism.
- Lim et al., 2017 — somatic mosaicism and low-VAF detection in TSC.

**Variant interpretation**
- Richards et al., 2015 — ACMG-AMP standards and guidelines for the interpretation of sequence variants (PVS1/PM2/PP4 framework).

**Therapeutics**
- French et al. (EXIST-3), 2016 — adjunctive everolimus for TSC-associated treatment-resistant focal seizures; first targeted therapy in a genetically defined epilepsy.

**TAND framework**
- de Vries et al. — the TAND consensus framework and the TAND-L lifetime checklist.

**Marshall-Hagedorn diagnostic-uncertainty methodology (the direct lineage of the TAND agent)**
- Marshall, Nickels, Brady, Hagedorn, 2023 — detecting diagnostic uncertainty in clinical documentation (Journal of Hospital Medicine).
- Nickels, Marshall, Edgerton, Brady, Hagedorn, Lee, 2024 — linguistic markers of diagnostic uncertainty (Applied Linguistics).
- Ipsaro, Patel, Marshall, Hagedorn, 2021 — uncertainty in inpatient documentation (Hospital Pediatrics).
- Orenstein, …, Hagedorn, 2021 — alert burden and signal discipline (JAMIA); the empirical caution behind the engine's "recalibrate if >~3 alerts/clinician/week" rule.

The TAND Surveillance Agent is an extension of this published body of work, not an external graft onto it.

### Appendix E: Mapping to the prior TSC papers

This paper is the consolidated technical reference. It sits alongside three earlier documents and an FAQ, each written for a different reader and at a different altitude. The table shows where a reader of each prior document finds the corresponding depth here, and flags where this paper supersedes earlier framing.

| Prior document | Audience / purpose | Where it maps in this paper | Notes on reconciliation |
|---|---|---|---|
| **Explainer** (the disease-and-opportunity overview) | Clinical and lay-informatics readers; why TSC, why now | Disease facts (body) + Appendix A glossary; the agent overview sections | The Explainer's narrative of the NMI cohort and TAND under-recognition is made operational here as the Variant Curator and TAND agent specs. |
| **Winslow Initiative** (Paper II — institutional plan) | CCHMC leadership; 24-month, four-phase build with decision gates | Architecture sections + the "sources feed the engine" model; Appendix C (Phase-1 vs. demo boundary) | The Winslow plan's Phase-1 Epic Clarity/Caboodle and LIMS integration is the work explicitly **not built** in this synthetic demo. This paper marks that boundary wherever it appears. |
| **Builder's Guide** (Paper III — the 8-week MVP) | Engineers building on DGX Spark + RunPod | The agent schemas (Appendix B), cohort spec (Appendix C), eval targets and the W1–W8 plan (body) | This paper is the canonical home of the schemas and eval numbers the Builder's Guide referenced narratively. |
| **FAQ** | Mixed skeptical readers; the do-not-overclaim answers | The do-not-overclaim posture throughout; provenance block (Appendix B.1); the human-gate language in every agent output | The FAQ's credibility commitments are encoded structurally here — `human_review_status`, `is_diagnosis: false`, `framing: decision_support_not_recommendation`. |

Two framing reconciliations apply across all four. First, the platform designation: earlier TSC documents were written in the "three engines / eleven agents" generation of the HCLS AI Factory and named a "Precision Intelligence Network." This paper is authoritative for the current framing — **Seven Engines, Eight Intelligence Agents, One Platform**, with the TSC Intelligence Engine as Engine 7 and no entity named "Precision Intelligence Network." Second, the institutional posture: any prior promotional or validation language (vendor RFP figures, "most complete pipeline" claims, named-institution endorsements) is, in this paper, stated at most as "presented at" or "in discussion with," or omitted. Where a prior document asserted more, this paper is the correction of record.

That is the full reference apparatus for the engine: the words, the contracts, the cohort, the literature, and the lineage to what came before. The argument of the paper is the build, and these appendices are what make it checkable.
