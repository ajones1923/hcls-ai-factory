# The Pediatric Patient Potential of the TSC Intelligence Engine

### What an Open, Give-It-Away Disease Engine Could Reach — A Sized, Sourced, and Honest Estimate

**Author:** Adam M. Jones — Architect, HCLS AI Factory; Product Marketing Manager, HPC, the platform
**Companion to:** *The TSC Intelligence Engine* (master research volume) and its build PRD · HCLS AI Factory, **Engine 7**
**License:** Apache 2.0 · 2026

> **What this document is, and is not.** This is an estimate of *reach* — how many children sit inside the documented clinical gaps the TSC Intelligence Engine is built to close — not a prediction of outcomes. Today the Engine runs on synthetic data; real-patient benefit requires the institutional Phase-1 validation described in the master paper, and every agent output is decision support behind a clinician's review. Nothing here is FDA-cleared, and no specific outcome is promised. The numbers below are order-of-magnitude estimates derived from published incidence and prevalence figures and the gap percentages cited in the master paper. They are meant to size an opportunity honestly, with the arithmetic shown.

---

## Implementation Status (June 2026): Built & Verified

The engine whose reach this document sizes is no longer only a plan. A working v0.1 is **implemented and running** on the NVIDIA DGX Spark — all five agents, the deterministic orchestrator, the synthetic 50-patient cohort, and three clinician surfaces — with 41 automated tests passing and the three-act demonstration running end-to-end. The deterministic components (the ACMG-AMP classifier, the Gaussian-process trajectory forecasts, the Marshall-Hagedorn discourse detection) are fully real, and the frontier-model reasoning has been verified live: the Variant Curator produced a genuine `claude-opus-4-8` ACMG narrative for the featured mosaic patient. The honest framing below is unchanged — all data is synthetic, real-patient benefit requires the institutional Phase-1 validation, every output is clinician-reviewed, and nothing is FDA-cleared. In short: it runs today, which is what makes the reach below a credible *could* rather than a hope.

---

## 1. Why "reach," and why pediatric

The case for this work has never been a single headline number. It is that a specific, well-documented set of gaps in Tuberous Sclerosis Complex (TSC) care affects a large, identifiable population of children — and that an open, replicable engine can address those gaps at a marginal cost low enough that the limiting factor stops being economics and becomes only the decision to deploy.

TSC is, for practical purposes, a pediatric disease in its highest-stakes window. It is diagnosed prenatally, in infancy, or in childhood. The events that most shape a life — infantile spasms in the first year, subependymal giant cell astrocytoma (SEGA) growth through childhood and adolescence, and the developmental and neuropsychiatric trajectory of TSC-Associated Neuropsychiatric Disorders (TAND) — are pediatric events. So when we size who this could reach, the great majority are children.

This document sizes that reach in three layers: the direct TSC pediatric population; the children inside each specific gap the five agents target; and the multiplier that comes from giving the work away under Apache 2.0 and replicating the pattern to other pediatric rare diseases.

## 2. Layer one — the direct TSC pediatric population

TSC occurs in approximately **1 in 6,000 live births**. The commonly cited prevalence figures, from the TSC Alliance, are roughly **50,000 people in the United States and on the order of 1 million worldwide** (some estimates range higher, toward 2 million). At ~3.6 million US births per year, that is roughly **600 new TSC diagnoses per year in the United States** alone.

Because diagnosis and the most consequential interventions concentrate in childhood, the majority of the actively surveilled, intervention-relevant population at any moment is pediatric. This is the base population the Engine serves wherever it is deployed.

## 3. Layer two — the children inside the gaps

The honest core of "potential" is not "lives saved." It is the count of children who sit inside the documented gaps each agent is purpose-built to close. The figures below apply the published gap percentages to the US (~50,000) and worldwide (~1,000,000) TSC populations.

| Gap (the agent that addresses it) | Prevalence within TSC | Estimated pediatric-weighted reach (US / worldwide) |
|---|---|---|
| **No-mutation-identified / somatic mosaicism** (TSC-Variant Curator) | ~10-15% have no genetic diagnosis on standard blood testing | **~5,000-7,500 / ~100,000-150,000.** A confirmed molecular diagnosis unlocks surveillance authorization, cascade family testing, trial eligibility, and prognostic clarity. Recoverable in roughly half of cases where disease-affected tissue is available. |
| **Under-recognized TAND** (TAND Surveillance Agent) | ~90% have TAND; **30-50% of features go unaddressed** (TOSCA) | **~13,000-22,000 / ~270,000-450,000** with neuropsychiatric features going unsurfaced. Families consistently rank TAND the single greatest day-to-day burden of the disease. This is the largest-reach gap. |
| **Surveillance gaps and lesion trajectory** (Phenome Mapper + Trajectory Modeler) | ~85% epilepsy; SEGA in ~10-20%; AML in ~80% by adulthood; lifelong multi-organ monitoring | **Touches the whole population.** Better-timed intervention can mean catching a SEGA approaching the foramen of Monro before a hydrocephalus emergency, not after; flagging an overdue surveillance scan before a lesion is missed. |
| **Fragmented trial access** (Therapeutics Strategist) | active but dispersed TSC trial landscape | Raises appropriate trial enrollment among the epilepsy, SEGA, and AML subsets — children who could benefit but are never matched. |

Bounded to TSC alone, the Engine addresses gaps that touch **tens of thousands of children in the United States and hundreds of thousands worldwide**, with TAND surveillance and surveillance-gap detection reaching nearly the entire population.

## 4. Layer three — the multiplier of giving it away

This is the layer that changes the number's character, and it is the part that matters most about the Apache 2.0 license.

### Open-source removes the gatekeeping

A commercial product reaches the patients of paying customers. An open engine reaches **the patients of every center that chooses to run it** — Cincinnati Children's, TGen, City of Hope, a TSC program at any academic medical center, and, because it runs on a single desktop-class NVIDIA DGX Spark with elastic burst compute, potentially centers in lower-resourced health systems that could never license enterprise clinical AI. The addressable population is not a customer base; it is the TSC population that any adopting center serves, at no licensing cost and with no vendor lock-in.

### Replication turns one disease into a class

The master paper's architecture deliberately separates the wiring — the agents, the deterministic orchestrator, the clinician surfaces, the provenance layer — from the "box labels" of a specific disease and institution. Replication is, in the paper's phrase, "swap the box labels, keep the wiring." The named second-wave targets are also overwhelmingly pediatric, and several are **more common than TSC**:

| Disease | Approximate incidence | Note |
|---|---|---|
| Neurofibromatosis type 1 (NF1) | ~1 in 3,000 births | ~2-3× as common as TSC; on the order of **2 million worldwide**. Shares the monogenic, multi-organ, surveillance-heavy profile. |
| Rett syndrome | ~1 in 10,000-15,000 female births | MECP2; strong genotype-phenotype and neurodevelopmental parallel to TAND. |
| Williams syndrome | ~1 in 7,500-10,000 births | Copy-number interpretation; distinctive cardiovascular/behavioral/cognitive profile. |
| Cornelia de Lange syndrome | ~1 in 10,000-30,000 births | Multi-system developmental disorder; substantial care-coordination need. |
| mTORopathy focal epilepsies (FCD II, hemimegalencephaly, DEPDC5/NPRL2/NPRL3) | large epilepsy subpopulations | The TSC genomic and trajectory pipeline extends with modest adaptation. |

Each additional disease reuses the hardware, the data-layer integrations, the governance machinery, the agent abstractions, and the institutional credibility built once — at an estimated one-third to one-half the cost of the first. Summed across the replication class, the engine pattern addresses a pediatric rare-disease population in the **several millions worldwide**, all sharing the same core obstacle the Engine attacks: *the evidence exists, but it is dispersed across modalities and time and is not integrated at the point of care.*

## 5. From reach to outcome — what it would take, and what it would look like

Reach is the population in the gaps. Outcome is the subset for whom closing a gap changes a life. Converting one to the other is exactly what the institutional engagement in the master paper is designed to do: real-data validation, Epic and biobank integration, governance, and prospective measurement. The mechanisms of impact, stated as the paper states them — plausibly, not as promises — are concrete:

- **Earlier and more complete diagnosis.** Recovering a mosaic diagnosis for a child in the no-mutation-identified cohort ends a diagnostic odyssey and unlocks surveillance, family testing, and trials.
- **Recognized rather than missed TAND.** Surfacing a pattern of academic and anxiety signals that no single visit formalized turns an unaddressed burden into an addressable one.
- **Better-timed intervention.** A forecast that a SEGA will approach the intervention threshold in a defined window converts a reactive emergency into a planned discussion.
- **Appropriate trial enrollment and surveillance adherence.** Matching a child to an eligible trial, or flagging an overdue scan, gives more children access to what already exists.

For any one child, the impact is one of those four things. Because the Engine is free and replicable, that single-child impact can compound across hundreds of thousands of children with TSC and millions across the rare-disease class — anywhere a center chooses to run it.

## 6. The honest bottom line

The truest statement of the pediatric potential is two sentences. **For one child: an earlier diagnosis, a TAND feature finally surfaced, a SEGA caught in time, or a trial finally matched.** **Across the field: because it is open and replicable, that single-child benefit can compound across a TSC population of ~1 million — mostly children — and a broader pediatric rare-disease class numbering in the millions, limited not by who can pay but only by who chooses to deploy.**

That compounding, accumulated across many children, many diseases, and many years, is the point. As the master paper closes: the impact is not in any single intervention; it is in the steady accumulation of better-informed care across many patients across many diseases across many years.

---

## Appendix — Methodology, sources, and caveats

**Population figures.** TSC incidence (~1 in 6,000 live births) is the standard figure used in the International TSC Consensus literature. US prevalence (~50,000) and worldwide prevalence (~1 million, with higher estimates toward 2 million) are the commonly cited TSC Alliance estimates. US annual birth count (~3.6 million) is used to derive ~600 new US diagnoses per year.

**Gap percentages.** The no-mutation-identified share (~10-15%), TAND prevalence (~90%) and under-recognition (~30-50% of features, per TOSCA registry findings), epilepsy prevalence (~85%), and AML prevalence (~80% by adulthood) are the figures cited and sourced in the master research volume (§3-§4), drawn from the 2021 ITSC consensus guidelines, the TOSCA natural-history registry, and the somatic-mosaicism literature (Tyburczy 2015; Giannikou 2016; Lim 2017).

**Reach arithmetic.** Gap-level reach is the gap percentage applied to the US (~50,000) and worldwide (~1,000,000) populations. The TAND under-recognition reach interprets the "30-50% of features unaddressed" figure as a proxy for the share of the ~90% TAND population with under-recognized features; this is an approximation, not a patient-level epidemiological measurement. The "recoverable in roughly half" note for mosaic recovery reflects that tissue-based detection requires disease-affected tissue (e.g., from epilepsy surgery) to be available.

**Replication incidences.** NF1 (~1 in 3,000), NF2 (~1 in 25,000-33,000), Rett (~1 in 10,000-15,000 female births), Williams (~1 in 7,500-10,000), and Cornelia de Lange (~1 in 10,000-30,000) are standard published incidence ranges; the "several million worldwide" cumulative figure for the replication class is an order-of-magnitude sum, not a precise total.

**Caveats.** These estimates size *potential reach*, not delivered clinical benefit. Realizing benefit for real patients depends on the institutional Phase-1 validation, real-data integration, governance, and prospective measurement described in the master paper. Every Engine output is decision support reviewed by a qualified clinician; the system is not FDA-cleared, and its Software-as-a-Medical-Device posture is undetermined. Figures are rounded and should be read as ranges, not point estimates.
