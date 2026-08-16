# HCLS AI Factory — Scorecard (DGX Spark)

**Date:** 2026-08-15 · **Scope:** the platform **as it exists on the DGX Spark**. GitHub and
hcls-ai-factory.org are **not scored** — neither was touched in this pass, and scoring them from
memory would be guessing.

---

## Overall: **6.5 / 10**

An unusually strong codebase sitting on a platform layer that does not yet enforce what the project
claims. The science is the best part; the controls are the weakest.

| Dimension | Score | Basis |
|---|---:|---|
| Clinical / scientific accuracy | **9.5** | 44 CPIC pairs all genuine; 0 gene-locus errors; correct variant IDs |
| Test health | **9.0** | 8,402 passing, 0 failing across 17/17 |
| Code quality | 8.0 | ruff-clean on real-bug rules; 12 stdlib-shadowing modules |
| Documentation (internal) | 7.5 | 17/17 have README + docs/; 5 lack the decision-support frame |
| Registry & honesty | 7.5 | now 0 live-but-unreachable, CI-gated — but only after 2 were found |
| Deployability | 7.0 | 17/17 declared; **nothing actually running**; images unbuilt |
| Layout consistency | 6.5 | 12/17 conformant; 3 legitimate patterns, 1 real defect |
| Security | **4.0** | **1 of 12** endpoints authenticated; governance gate adopted by 1 of 12 |
| Hardware utilisation | **2.0** | GPU at **0%** the entire session; CPU-only torch on a GB10 |

---

## Where the marks come from

### 9.5 — Clinical and scientific accuracy

The strongest dimension and the one that matters most for the audience. Every checkable claim held:

- **44 pharmacogenomic pairs across 13 genes, all genuine CPIC** — no spurious pair found.
- **Zero gene→locus errors** (APOE 19q13.32, CYP2C19 10q23.33, LDLR 19p13.2, LMNA 1q22,
  MYH7 14q11.2, TTN 2q31.2).
- Variant identifiers correct: VKORC1 **rs9923231**, SLCO1B1 **rs4149056**.

Deduction: CYP3A4→midazolam/cyclosporine are presented alongside CPIC pairs. They are correct
pharmacology but not CPIC guideline pairs and should be labelled as substrate relationships.

### 4.0 — Security

**1 of 12 FastAPI entrypoints enforces authentication.** The other eleven — every intelligence
agent plus the cardiology, imaging and oncology engines — expose clinical decision-support endpoints
with no auth dependency.

Worse, the platform *has* a governance module and **1 of 12 services uses it**. And
`install_governance` emits an `X-HCLS-Governed` response header while performing no gating: its own
`/governance` payload describes input-validation and output-honesty as things a handler must call.
A header that asserts governance on an ungoverned request is the wrong failure mode for this project.

Credit where due: **zero real secrets** in tracked code, `.gitignore` correct, all 27 Dockerfiles
pin base tags, and the two unhardened service images were fixed this session.

### 2.0 — Hardware utilisation

```
torch 2.10.0+cpu     torch.cuda.is_available() = False
GPU across 200 minutes: 0% utilisation, 12.3–16.0 W (idle floor)
```

The DGX Spark's entire premise is the GB10 and 128 GB of unified memory. Peak memory used during a
full platform test sweep was **19 GB of 119 (16%)**, and the GPU was never engaged. Gated behind
NGC credentials — see `GATED_SOFTWARE_PRD.md`.

### 7.0 — Deployability

Moved from 11/17 to **17/17 declared** this session, with the port map resolved and gated. But
declaring is not running: **no engine or agent container exists on the box**, and per-service venvs
exist for only 4 of ~20 supervised services. The platform is cold.

---

## Recommendations, ranked by value

| # | Action | Lifts | Effort |
|---|---|---|---|
| 1 | **Authenticate the 11 open clinical endpoints** | Security 4.0 → 8 | medium |
| 2 | **Stop emitting `X-HCLS-Governed` unless the request was gated** | Honesty | low |
| 3 | **Install CUDA PyTorch** (gated) | Hardware 2.0 → 8 | gated |
| 4 | Adopt `install_governance` + `require_valid_input` across all 12 | Security, honesty | medium |
| 5 | Bring the platform up and keep it up | Deployability 7 → 9 | medium |
| 6 | Add the decision-support frame to genomic-foundation and precision-intelligence | Docs, honesty | low |
| 7 | Choose one runtime model (compose) | Layout, ops | decision |
| 8 | Rename the 12 stdlib-shadowing modules | Code quality | medium |
| 9 | Agree a test-depth floor (range is 4 → 1,966) | Confidence | medium |
| 10 | Label CYP3A4 pairs as substrate, not CPIC | Accuracy 9.5 → 10 | trivial |

**Items 1 and 2 are the highest-value changes in this document.** The platform's credibility rests
on honesty-by-construction, and it currently ships a header claiming governance it does not perform,
in front of unauthenticated clinical endpoints.

---

## What changed during this session

| | Before | After |
|---|---|---|
| Merge gate | could not execute | passes |
| Subject tests | unknown | 8,402 / 0 fail |
| Compose coverage | 11/17 | 17/17 |
| Port conflicts | 2 agents on the same 2 ports | one convention, CI-gated |
| `live` but unreachable | 2 | 0 |
| Unhardened service images | 2 | 0 |

## Not scored, and why

**GitHub** and **hcls-ai-factory.org** were not reviewed, updated or deployed in this pass. Scoring
them would require reading the published state, which I have not done. They remain open items
(asks 6–11).
