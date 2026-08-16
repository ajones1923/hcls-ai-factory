---
title: "Honesty by Construction: A Verifiable Open-Source Precision-Medicine Platform on a Single Workstation"
description: A preprint describing the HCLS AI Factory — eight engines, eight intelligence agents and one disease program on a single workstation — and the capability registry that makes its claims machine-checkable.
---

# Honesty by Construction

### A Verifiable Open-Source Precision-Medicine Platform on a Single Workstation

**Preprint · 2026** · Apache-2.0 · [github.com/ajones1923/hcls-ai-factory](https://github.com/ajones1923/hcls-ai-factory)

---

## Abstract

Clinical AI platforms routinely describe capabilities that cannot be exercised by a reader. The
resulting credibility gap is not primarily technical — it is evidentiary. We present the **HCLS AI
Factory**, an open-source precision-medicine platform comprising eight computational engines, eight
retrieval-augmented intelligence agents and one disease program, deployable on a single ~$4,700
workstation.

Our contribution is less the pipeline than the **mechanism that constrains what the pipeline may
claim**. A machine-readable capability registry is the single source of truth for what exists, what
is served and on which port. Continuous integration fails if a documented capability contradicts the
registry, if two services claim the same port, or if the process supervisor disagrees with the
manifest. Public documentation — including the capability maturity matrix — is *generated from* that
registry rather than written alongside it. Demonstrations are labelled **LIVE**, **REPRESENTATIVE**
or **BURST**, and the demonstration runner refuses to execute a demonstration labelled LIVE whose
service is unreachable rather than returning a cached result.

We report measurements from applying these mechanisms to a codebase of ~385,000 lines: **8,402
passing tests across all 17 subjects**, and the detection of two capabilities registered as `live`
while no process was bound to their advertised ports. We describe what the platform cannot currently
do, including that its GPU was measured at **0% utilisation across a 203-minute working session**
because the installed tensor runtime is CPU-only.

We argue that verifiable negative claims — an explicit, enforced account of what a system does *not*
do — are a practical requirement for clinical AI, and that they are cheap to implement once
capability state is machine-readable.

---

## 1. Introduction

A clinician evaluating a computational platform faces a specific epistemic problem: the
documentation is written by the same people who wrote the code, and nothing structurally prevents
the two from diverging. Screenshots persist after features regress. Roadmap items acquire the
present tense. A capability marked "available" may be served by a mock.

This is usually treated as a documentation-hygiene problem. We treat it as an architectural one.

The platform described here performs a recognisable precision-medicine workflow: sequencing reads to
variants, variants to interpretation, interpretation to candidate molecules, with specialist agents
reasoning over curated corpora. That workflow is not novel. What we believe is worth reporting is a
set of enforcement mechanisms that make the platform's public claims **checkable by a third party
without trusting the authors**, and the empirical result of turning those mechanisms on a
substantial codebase.

### Contributions

1. A **capability registry** as the single source of truth for capability existence, status and
   network location, enforced in CI.
2. **Documentation generated from that registry**, so published claims cannot silently drift.
3. A **three-valued demonstration label** (LIVE / REPRESENTATIVE / BURST) with a runner that fails
   rather than degrades.
4. An **empirical account** of applying these to ~385,000 lines, including defects found and the
   platform's current limits.

---

## 2. Architecture

### 2.1 Composition

| Layer | Count | Role |
|---|---:|---|
| Engines | 8 | Deterministic computation — genomics, interpretation, therapeutic discovery, imaging, oncology, cardiology, structural biology, single-cell |
| Intelligence agents | 8 | Retrieval-augmented reasoning over curated corpora, returning cited output |
| Disease programs | 1 | A vertical composing the horizontal engines and agents for one condition |

The engine/agent split is deliberate and load-bearing. Engines compute; agents interpret. The
single-cell pair illustrates it: the engine clusters an expression matrix and annotates clusters by
marker-gene overlap; the agent explains what those clusters imply. Keeping computation and
interpretation in separate services means each can be checked independently — a numerical result can
be reproduced without an LLM, and an interpretation can be traced to the numbers it rests on.

A disease program is **not** a ninth engine. It is a vertical that composes existing horizontal
capabilities.

### 2.2 The capability registry

Every capability is declared in a single JSON manifest with an identifier, type, status
(`live` / `planned`), a network endpoint and a description. At the time of writing the registry holds
**42 capabilities**.

The manifest is not documentation. It is the input to:

- **coverage validation** — every engine and agent directory must map to a registered capability;
- **a port drift-guard** — no two capabilities may claim the same port;
- **a convention guard** — the registry advertises the UI port and the API is UI + 1;
- **a supervisor cross-check** — the process supervisor is parsed and compared against the manifest;
- **site generation** — the public capability maturity matrix is emitted from the manifest.

All of these run in continuous integration. A contradiction fails the build.

### 2.3 Deployment target

The reference target is a single workstation: a 20-core ARM CPU, an integrated GPU, and 128 GB of
unified memory, at approximately $4,700. The claim being tested is that a clinically meaningful
precision-medicine stack fits within one such machine's reach — not that every component runs
locally. Components that are architecturally incompatible or computationally excessive run on remote
accelerators over a private mesh, and are described as **elastic burst** wherever they appear.

---

## 3. Honesty by construction

### 3.1 Status must follow evidence

A capability's status is a factual claim about the world: `live` asserts that a real service answers
on a real port. We enforce that a `live` capability answer a health probe, and treat a failure as a
build defect rather than a documentation nit.

This is not hypothetical. Applying the check found **two capabilities registered `live` with nothing
bound to their advertised ports**. Neither was a mock; both were real, tested code that had never
been given a deployment path. One was containerised; the other was demoted to `planned`, because it
comprised five independent services with no aggregator on the port the registry advertised.

We note the asymmetry that makes this class of error persistent: the code passed its tests, the
documentation described it accurately, and only the *conjunction* of registry, network state and
documentation was false. No single artefact was wrong.

### 3.2 Three-valued demonstration labels

A demonstration is labelled exactly one of:

| Label | Assertion |
|---|---|
| **LIVE** | Ran now, on real input, on the machine in front of the audience |
| **REPRESENTATIVE** | A pre-computed or curated result stands in for a long or unavailable step |
| **BURST** | Ran live, on remote accelerators over a private mesh |

The label is enforced by the runner, which **refuses to execute a LIVE demonstration whose service
is unreachable** rather than silently substituting a cached result. A demonstration that cannot run
reports what is missing.

Of 17 demonstrations, **6 are currently LIVE**, 9 are REPRESENTATIVE pending gated components, and
2 are mixed. We report that distribution rather than the distribution we would prefer.

### 3.3 Generated documentation

Per-subject documentation — 68 guides across the 17 subjects — is generated from repository facts:
ports, test counts, registered capabilities, module inventories and declared dependencies. A guide
cannot state a port the registry does not allocate, because it does not know one.

During generation we found a class of error this prevents: single-service portals have no separate
API port, and naively deriving "API = UI + 1" would have published one portal's API as another's UI.

---

## 4. Evaluation

### 4.1 Test health

A unified harness executes all 17 subject suites:

```
17 subjects · 8,402 passed · 0 failed · 0 errors
```

Two properties of the harness are worth reporting because both produced **false readings** before
being identified:

1. Eleven subjects ship a module whose filename shadows a Python standard-library module. Placing
   those directories on the import path terminates the interpreter *before test collection*. A naive
   harness reports the platform as entirely broken.
2. Two scientific packages register test-framework plugins that import application configuration at
   start-up, aborting collection for three subjects. Disabling plugin autoloading wholesale fixes
   those three and breaks 38 asynchronous tests elsewhere; only disabling the two specific plugins
   is correct.

We report these because a uniform result across every subject is a signal that the *harness* is
wrong, not that the platform is uniformly broken — an inference we initially failed to make.

### 4.2 Reproducibility on a clean clone

The published quickstart was executed on a fresh clone before publication. On such a clone the
harness reports:

```
17 subjects · 8,194 passed · 0 failed · 0 errors · 213 skipped
```

The 213 skips are correct and deliberate: reference data is not published, and tests that assert on
it skip rather than fail. Earlier those tests failed hard, producing 44 failures on any clone — a
result that reads as a broken project while the project is sound. We consider "fails cleanly for the
right reason" a reportable property.

One demonstration runs end-to-end on a clean clone with no accelerator and no credentialed
software: a single-cell workflow over a public 2,700-cell dataset resolving to 9 clusters and 7
annotated cell types.

### 4.3 Clinical content

Machine-checkable clinical assertions were verified against primary sources:

- **44 gene–drug pairs** across 13 genes checked against consortium guidance; all genuine.
- **Zero gene→cytoband errors** across every such assertion in the codebase.
- Variant identifiers verified individually.

One error was found and corrected: a retrieval corpus asserted that consortium guidance references
*CYP3A4* for tacrolimus dosing, when the actionable guideline is *CYP3A5*-based. We note the failure
mode — a plausible, adjacent, single-sentence error inside a corpus that agents retrieve from — as
one that automated checks catch poorly and that domain review catches immediately.

### 4.4 Resource characterisation

Sampled at 30-second intervals across a 203-minute working session:

| Metric | Value |
|---|---|
| CPU | mean 4.4%, peak 34.2% of 20 cores |
| Memory | peak 19 GB of 119 GB |
| **GPU utilisation** | **0%** |
| GPU power | 12.3–16.0 W (idle) |

**The accelerator was never engaged.** The installed tensor runtime is a CPU-only build, so every
model path is degraded or non-functional. A full 8,402-test sweep consumed 16% of available memory.
We report this because it is the most consequential limitation of the current deployment and is
invisible from the test results.

---

## 5. Limitations

We state these as the platform's own honesty framework requires.

1. **No accelerated inference.** The GPU is unused; a credentialed tensor runtime is required.
2. **Accelerated secondary analysis is unavailable.** Alignment and variant calling depend on a
   licensed container. Genomics demonstrations serve pre-computed results and are labelled
   REPRESENTATIVE.
3. **Generative molecular design and docking are unavailable.** Both depend on credentialed
   inference microservices; the flagship discovery demonstration is REPRESENTATIVE.
4. **Nothing is continuously deployed.** All 17 subjects are declared; the reference machine runs
   none of them persistently.
5. **Authentication ships disabled.** Every service can enforce an API key and fails closed once
   configured, but the default preserves an open trusted-network posture.
6. **Test depth is uneven** — from 4 to 1,966 tests per subject. Pass rate is uniform and therefore
   uninformative.
7. **All output is decision support.** No component performs diagnosis or prescribing. Gene-therapy
   approaches referenced for the disease program are preclinical.

---

## 6. Discussion

The mechanisms described are neither difficult nor expensive. A capability registry is a JSON file.
The port drift-guard is roughly forty lines. The demonstration runner's central rule — fail rather
than degrade — is one conditional.

What they cost is the willingness to publish a number you did not choose. Our demonstration
distribution is 6 LIVE of 17; our GPU utilisation is 0%; our clean-clone test count is lower than
our developer-machine count. Each of these is a worse-looking figure than the alternative, and each
is the one a reader can verify.

We suggest the relevant question for a clinical AI platform is not *what does it claim* but *what
would have to be true for the claim to be false, and does anything check that*. For the platform
described here, the answers are enumerable: a health probe, a port comparison, a supervisor parse, a
generated matrix, and a demonstration runner that refuses to pretend.

---

## 7. Availability

Source, documentation and the complete audit trail — including the defects reported here — are
public under Apache-2.0 at
[github.com/ajones1923/hcls-ai-factory](https://github.com/ajones1923/hcls-ai-factory).

Every figure in this paper is reproducible from a clean clone:

```bash
git clone https://github.com/ajones1923/hcls-ai-factory && cd hcls-ai-factory
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -e lib/hcls_common
.venv/bin/python scripts/run_all_tests.py        # test health
.venv/bin/python scripts/validate_registry.py    # registry, ports, supervisor
.venv/bin/python scripts/run_demo.py E8          # a live demonstration
```

Reference genomic material used for verification is publicly consented; no patient data was used at
any point.
