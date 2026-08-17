---
title: "The HCLS AI Factory: An Open-Source Precision-Medicine Platform on a Single Workstation"
description: Eight computational engines, eight retrieval-augmented intelligence agents and one disease program — raw sequencing reads to ranked drug candidates on one workstation, with every capability claim machine-checkable.
---

# The HCLS AI Factory

### An Open-Source Precision-Medicine Platform on a Single Workstation

**Preprint · 2026** · Apache-2.0 · [github.com/ajones1923/hcls-ai-factory](https://github.com/ajones1923/hcls-ai-factory)

---

## Abstract

Precision medicine is fragmented across specialist tools that rarely compose: a variant caller here,
an annotation service there, a docking package elsewhere, each with its own data model, and none
aware of the patient. Institutions that want the whole path — sequencing reads to a therapeutic
hypothesis — assemble it themselves, at a scale of cost and integration effort that excludes most of
them.

We present the **HCLS AI Factory**, an open-source platform that composes that path on a single
~$4,700 workstation. It comprises **eight computational engines**, **eight retrieval-augmented
intelligence agents**, and a **disease program** that composes them vertically for one condition.
Engines compute; agents interpret; disease programs orchestrate. The separation is deliberate and is
what makes each layer independently verifiable.

The platform performs GPU-accelerated germline variant calling, variant annotation against curated
clinical databases with retrieval-augmented interpretation, a ten-stage therapeutic-discovery
pipeline producing ranked candidate molecules, DICOM imaging analysis with cross-modal reasoning
against genomics, molecular tumour-board packet generation, cardiovascular risk workflows,
structural-biology prediction and design, and single-cell transcriptomic analysis. Eight agents
provide cited clinical reasoning across pharmacogenomics, CAR-T, biomarkers, autoimmunity,
neurology, rare disease, trial matching and single-cell interpretation.

Every capability is declared in a machine-readable registry that is the input to continuous
integration, not a description of it. Documentation, including the public capability-maturity
matrix, is generated from that registry. We report **8,402 passing tests across all 17 subjects**,
44 pharmacogenomic gene–drug pairs verified against consortium guidance, and a single-cell workflow
reproducible from a clean clone on commodity CPU.

All clinical output is decision support for a qualified clinician.

---

## 1. Introduction

The distance between a patient's genome and a therapeutic hypothesis is not conceptually difficult.
It is a sequence of well-understood steps: align reads, call variants, annotate against clinical
knowledge, identify a druggable target, obtain or predict its structure, generate candidate
molecules, score them.

In practice each step lives in a different tool, maintained by a different community, with a
different data model and a different notion of what a "patient" is. The integration burden falls on
the institution, and it is substantial enough that the complete path is assembled only where there
is dedicated engineering capacity. The result is that precision medicine is, in practice, available
in proportion to institutional wealth.

This platform is an argument that the whole path fits on one machine — a 20-core ARM workstation
with an integrated GPU and 128 GB of unified memory, at roughly $4,700 — and that the code to do it
can be given away.

### 1.1 Design commitments

**Composition over monoliths.** Capabilities are independently addressable services with declared
inputs and outputs, not modules of a single application. A workflow composer wires them.

**Computation separated from interpretation.** Engines produce numbers; agents reason over them.
A clustering result can be reproduced without a language model; an interpretation can be traced back
to the numbers it rests on.

**Verticals compose horizontals.** A disease program is not a ninth engine. It is an orchestration
of the eight horizontal engines and eight agents for one condition.

**Claims are machine-checkable.** What the platform says about itself is generated from the same
manifest that continuous integration validates.

---

## 2. System architecture

### 2.1 The eight engines

| # | Engine | Function |
|---|---|---|
| E1 | **Genomic Foundation** | GPU germline variant calling — accelerated alignment (BWA-MEM2 via `fq2bam`) and DeepVariant, FASTQ → BAM → VCF; a queryable variant store with Ts/Tv quality gating, VAF and somatic-mosaicism handling; ACMG secondary-findings joins; single-box GWAS |
| E2 | **Precision Intelligence** | Variant annotation against curated clinical databases and missense-pathogenicity predictions, with retrieval-augmented interpretation that identifies and contextualises druggable targets |
| E3 | **Therapeutic Discovery** | Target → ranked candidates across a ten-stage pipeline: generation, conformer enumeration, docking, and cheminformatic scoring |
| E4 | **Clinical Imaging** | DICOM analysis with FHIR R4 export and cross-modality reasoning against genomics |
| E5 | **Precision Oncology** | Molecular tumour-board packet generation, therapy ranking and trial matching; fusion-first paediatric logic with late-effects-weighted ranking |
| E6 | **Cardiology** | Eleven clinical workflows and six risk calculators spanning prevention, intervention and rhythm; coronary-calcium and ten-year risk logic following professional-society guidance |
| E7 | **Structural Biology** | Protein structure prediction and design — folding, embedding-based search with Smith-Waterman re-ranking, sequence design, developability and immunogenicity |
| E8 | **Single-Cell** | scRNA-seq: QC → normalisation → highly-variable-gene selection → PCA → neighbourhood graph → Leiden clustering → marker differential expression, with marker-based cell-type annotation |

### 2.2 The eight intelligence agents

Each agent performs retrieval-augmented reasoning over a curated domain corpus and returns cited
output: **pharmacogenomics** (with drug-metabolism filtering), **CAR-T** (including construct and
sequence stores), **precision biomarker**, **precision autoimmune**, **neurology** (spanning ten
domains), **rare-disease diagnostics**, **clinical trial matching**, and **single-cell
interpretation**.

The last is the clearest illustration of the engine/agent split: E8 computes cell-type clusters; the
agent interprets them — tumour-microenvironment profiling, drug-response hypotheses,
ligand–receptor analysis. The two are separate services precisely so the computation can be checked
without the interpretation.

### 2.3 Composition and orchestration

A capability registry gives every capability an identifier, a type, a status, a network endpoint and
declared input/output shapes. A workflow composer uses those shapes to wire capabilities into
pipelines; a tool surface exposes them to agents. Because ports and shapes are declared rather than
discovered, a composition that cannot work is rejected before it runs.

### 2.4 Hardware and elastic burst

The reference machine is a single workstation. Components that are architecturally incompatible or
computationally excessive run on remote accelerators over a private mesh. We describe this as
**elastic burst** wherever it appears, because "runs on one box" is a claim that must survive
inspection: patient data remains local, and only derived, non-identifying work bursts.

---

## 3. The workflow

### 3.1 Reads to variants

E1 performs accelerated alignment and variant calling, producing a QC'd variant substrate. Quality
is gated on transition/transversion ratio; the store handles variant allele fraction and mosaicism
rather than assuming germline heterozygosity. ACMG secondary-findings logic runs as a join against
the called set.

Reference verification uses a publicly consented benchmark genome. No patient data is used.

### 3.2 Variants to interpretation

E2 annotates against curated clinical-significance databases and missense-pathogenicity predictions,
then applies retrieval-augmented interpretation over that annotated set. The output is a cited
interpretation naming candidate druggable targets — the citation being the part a clinician can
check.

### 3.3 Target to structure

E7 supplies the structural context a designed molecule must fit: folding for sequences without a
deposited structure, embedding-based similarity search with Smith-Waterman re-ranking for
homologues, and sequence design against a fixed backbone. Developability and immunogenicity scoring
apply to biologic candidates.

### 3.4 Structure to candidates

E3 runs a ten-stage pipeline from target to ranked candidates. Molecule generation, conformer
enumeration, property calculation and ranking are performed with an open cheminformatics toolkit;
docking places candidates in the target site. Generated leads are **preclinical design candidates,
not drugs** — the pipeline compresses the earliest and most combinatorial step of discovery,
choosing what to synthesise, not the trials that follow.

### 3.5 Cross-modal joins

The workflow is not a single line. E4 reads imaging and hands findings to genomics for cross-modal
reasoning; E6 contextualises cardiovascular risk and hands to pharmacogenomics for whether a given
drug will work in that patient; E5 assembles tumour-board packets drawing on imaging, biomarker and
trial capabilities simultaneously. E8's clusters feed oncology's microenvironment analysis.

Composition, not a pipeline, is the architecture.

---

## 4. The disease program

A **disease program** is a vertical that composes the horizontal engines and agents for one
condition. The reference implementation targets tuberous sclerosis complex — a genetic condition in
which loss of function in either of two genes releases a brake on a growth-regulating pathway,
producing lesions across brain, kidney, heart, skin and lung.

It is a good test of composition precisely because it is multi-system: a condition confined to one
organ would not exercise the platform.

The program orchestrates five disease-specific sub-agents — variant curation, trajectory modelling,
therapeutics strategy, phenotype mapping, and neuropsychiatric surveillance — each of which calls
horizontal capabilities rather than reimplementing them. The variant curator uses E1's store and
E2's annotation; the therapeutics strategist uses E3 and the pharmacogenomics agent; trial matching
uses the trial agent.

Three constraints are stated in the program itself, not in a footnote: output is decision support
and never autonomous diagnosis or prescribing; paediatric caution applies at full force; and
gene-therapy correction of the underlying genes is **preclinical** — an open design and analysis
bench, not a treatment available today. An approved pathway-inhibiting therapy exists and is
handled as such.

---

## 5. Verification

Clinical platforms are difficult to evaluate from outside because documentation and implementation
are produced by the same people and nothing structurally prevents divergence. We treat that as an
architectural problem.

### 5.1 The registry as source of truth

The capability manifest — **42 capabilities** at the time of writing — is the input to continuous
integration. CI fails if an engine or agent directory has no registered capability, if two
capabilities claim the same port, if the port convention is violated, or if the process supervisor
disagrees with the manifest. The public capability-maturity matrix is generated from the same file.

### 5.2 Status follows evidence

A status of `live` asserts that a service answers on a stated port. Enforcing that against reality
found **two capabilities registered `live` with nothing bound to their advertised ports**. Neither
was a mock; both were real, tested code that had never been given a deployment path. One was
containerised; the other was demoted, because it comprised five independent services with no
aggregator on the advertised port.

The instructive part is that no single artefact was wrong. The code passed its tests. The
documentation was accurate. Only the *conjunction* of registry, network state and documentation was
false — which is why a structural check finds it and review does not.

### 5.3 Demonstration labelling

Each demonstration is labelled **LIVE** (running now on real input), **REPRESENTATIVE** (a
pre-computed result standing in for a long or unavailable step) or **BURST** (running live on remote
accelerators). The runner **refuses to execute a LIVE demonstration whose service is unreachable**
rather than substituting a cached result.

### 5.4 Generated documentation

Per-subject documentation — 68 guides across the 17 subjects — is generated from repository facts:
ports, test counts, registered capabilities and module inventories. A guide cannot state a port the
registry does not allocate.

---

## 6. Results

### 6.1 Test health

```
17 subjects · 8,402 passed · 0 failed · 0 errors
```

From a clean clone, without local reference data, the same command reports **8,194 passed with 213
skips** — tests that assert on unpublished reference data skip rather than fail, so a fresh clone
reports honestly rather than appearing broken.

### 6.2 Verified clinical content

- **44 pharmacogenomic gene–drug pairs across 13 genes**, each checked against consortium guidance;
  all genuine. Variant identifiers verified individually.
- **Zero gene→cytoband errors** across every such assertion in the codebase.
- One corpus error found and corrected: an assertion attributing a tacrolimus dosing guideline to
  the wrong cytochrome gene. We note the failure mode — a plausible, adjacent, single-sentence error
  inside a corpus that agents retrieve from — as one that automated checks catch poorly and domain
  review catches immediately.

### 6.3 Reproducible computation

A single-cell workflow runs end-to-end from a clean clone with no accelerator and no credentialed
software, resolving a public 2,700-cell peripheral-blood dataset into **9 clusters and 7 annotated
cell types** consistent with the canonical marker panel. The dataset is fetched by the analysis
library itself, so the result reproduces on any clone.

### 6.4 Resource envelope

A full 17-subject test sweep on the reference workstation peaks at **19 GB of 119 GB** of unified
memory and 34% of 20 cores. The machine is not the constraint at this stage of the workload.

---

## 7. Limitations

- **Accelerated inference is enabled in one environment, not across the platform.** A CUDA build
  of the tensor runtime is installed and reports the integrated GPU correctly, but it backs a
  single service environment; the platform interpreter and the remaining services carry a
  CPU-only build, so most model paths do not reach the GPU today.
- **Accelerated secondary analysis requires a licensed container.** Genomics demonstrations
  currently serve pre-computed results and are labelled accordingly.
- **One generation backend of two is deployed.** The open cheminformatics generation path is live;
  an accelerated inference-microservice path is not, and stage-7 docking depends on a credentialed
  service.
- **Imaging modality coverage is partial.** Chest-radiograph analysis is verified on real DICOM
  pixels; the CT reasoning path is planned, and the cardiac CT workflow is explicitly representative.
- **The structural-biology engine is aggregated only at the model level.** Five services are
  reachable individually; the engine-level endpoint is not yet served.
- **Nothing is continuously deployed.** All 17 subjects are declared and buildable; the reference
  machine does not run them persistently.
- **Authentication ships disabled.** Every service can enforce an API key and fails closed once
  configured; the default preserves an open trusted-network posture.
- **Test depth is uneven**, from 4 to 1,966 tests per subject. Pass rate alone is uninformative.
- **All output is decision support.** No component performs diagnosis or prescribing.

---

## 8. Discussion

Two claims in this work are separable.

The first is that the precision-medicine path fits on one affordable machine. That is an engineering
claim, and the evidence is the platform: eight engines, eight agents, a disease program, and a test
suite a reader can run.

The second is that a platform of this kind should be able to state what it does *not* do, in a form
that is checked rather than asserted. The mechanisms are inexpensive — the registry is a JSON file,
the port guard is roughly forty lines, the demonstration runner's central rule is one conditional.
What they cost is the willingness to publish a number you did not choose.

We suggest the useful question to ask of any clinical AI platform is not what it claims, but what
would have to be true for the claim to be false, and whether anything checks that.

---

## 9. Availability

Source, documentation and the complete audit trail are public under Apache-2.0 at
[github.com/ajones1923/hcls-ai-factory](https://github.com/ajones1923/hcls-ai-factory).

```bash
git clone https://github.com/ajones1923/hcls-ai-factory && cd hcls-ai-factory
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -e lib/hcls_common
.venv/bin/python scripts/run_all_tests.py        # 17 subjects
.venv/bin/python scripts/validate_registry.py    # registry, ports, supervisor
.venv/bin/python scripts/run_demo.py E8          # a live demonstration
```

Reference genomic material used for verification is publicly consented. No patient data was used at
any point in development or evaluation.
