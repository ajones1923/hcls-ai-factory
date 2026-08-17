# HCLS AI Factory — Scorecard (DGX Spark)

**Date:** 2026-08-15 · **Scope:** the platform **as it exists on the DGX Spark**. GitHub and
hcls-ai-factory.org are **not scored** — neither was touched in this pass, and scoring them from
memory would be guessing.

---

## Overall: **7.6 / 10**  *(6.5 before the 2026-08-16 security work; 7.3 before the hardware
re-measurement below)*

An unusually strong codebase sitting on a platform layer that does not yet enforce what the project
claims. The science is the best part; the controls are the weakest.

| Dimension | Score | Basis |
|---|---:|---|
| Clinical / scientific accuracy | **9.8** | 44 CPIC pairs all genuine; 0 gene-locus errors; correct variant IDs |
| Test health | **9.0** | 8,402 passing, 0 failing across 17/17 |
| Code quality | 8.0 | ruff-clean on real-bug rules; 12 stdlib-shadowing modules |
| Documentation (internal) | 7.5 | 17/17 have README + docs/; 5 lack the decision-support frame |
| Registry & honesty | 7.5 | now 0 live-but-unreachable, CI-gated — but only after 2 were found |
| Deployability | 7.0 | 17/17 declared; **nothing actually running**; images unbuilt |
| Layout consistency | 6.5 | 12/17 conformant; 3 legitimate patterns, 1 real defect |
| Security | **7.5** *(was 4.0)* | 12/12 endpoints can enforce auth, fail-closed when configured; header no longer overclaims |
| Hardware utilisation | **5.0** | CUDA works, but only one service venv has it; GPU idle all session |

---

## Where the marks come from

### 9.8 — Clinical and scientific accuracy

The strongest dimension and the one that matters most for the audience. Every checkable claim held:

- **44 pharmacogenomic pairs across 13 genes, all genuine CPIC** — no spurious pair found.
- **Zero gene→locus errors** (APOE 19q13.32, CYP2C19 10q23.33, LDLR 19p13.2, LMNA 1q22,
  MYH7 14q11.2, TTN 2q31.2).
- Variant identifiers correct: VKORC1 **rs9923231**, SLCO1B1 **rs4149056**.

One real error found and fixed: the knowledge base asserted *"CPIC guidelines reference CYP3A4 for
tacrolimus dosing."* CPIC's tacrolimus guideline is **CYP3A5**-based. (The earlier claim that CYP3A4
inhibitor tables overstated CPIC backing was retracted — those tables were correctly structured.)

### 7.5 — Security *(was 4.0)*

**Now 12 of 12 entrypoints can enforce authentication**, via one shared module rather than twelve
copies, fail-closed once `HCLS_API_KEY` is set and verified live (401 without a key, 401 with a
wrong key, routing reached with the right one). The `X-HCLS-Governed` header now reports only gates
that actually ran.

Not yet 10: the gate is **off by default** (deliberate — it preserves the trusted-network demo
posture), and the handlers still do not call `require_valid_input()` / `honesty_flags()`, so
input-validation and output-honesty remain opt-in per route.

The original finding, for the record:

**1 of 12 FastAPI entrypoints enforced authentication.** The other eleven — every intelligence
agent plus the cardiology, imaging and oncology engines — expose clinical decision-support endpoints
with no auth dependency.

Worse, the platform *has* a governance module and **1 of 12 services uses it**. And
`install_governance` emits an `X-HCLS-Governed` response header while performing no gating: its own
`/governance` payload describes input-validation and output-honesty as things a handler must call.
A header that asserts governance on an ungoverned request is the wrong failure mode for this project.

Credit where due: **zero real secrets** in tracked code, `.gitignore` correct, all 27 Dockerfiles
pin base tags, and the two unhardened service images were fixed this session.

### 5.0 — Hardware utilisation *(was scored 2.0 on a bad measurement)*

```
platform .venv                  torch 2.10.0+cpu    cuda.is_available() = False
therapeutic-discovery/venv      torch 2.12.1+cu130  cuda.is_available() = True   (GB10, sm_121)
GPU across 200 minutes: 0% utilisation, 12.3–16.0 W (idle floor)
```

The DGX Spark's entire premise is the GB10 and 128 GB of unified memory. Peak memory used during a
full platform test sweep was **19 GB of 119 (16%)**, and the GPU was never engaged during it.

> **Correction (2026-08-16).** This section originally read *"CPU-only torch on a GB10"* and scored
> **2.0**, and the recommendations below listed "install CUDA PyTorch" as a gated blocker. That was
> wrong: it measured the **platform venv only**. `core/engines/therapeutic-discovery/venv` carries
> `torch 2.12.1+cu130` and reports the GB10 at `sm_121`, so CUDA is present and working on this box.
> The real gap is narrower and different — the CUDA build is in **one service venv**, not the
> platform interpreter or the other services, so most code paths still cannot reach the GPU. The
> lesson is the one this project keeps relearning: a measurement taken through a single interpreter
> is a claim about that interpreter, not about the machine.

**Also observed:** with ~95 GiB of the 119 GiB unified memory held by the host page cache, CUDA
allocation fails outright (`cudaMemGetInfo` → out of memory). On unified memory the page cache and
the GPU compete for the same pages, so a cache drop is part of the bring-up, not an afterthought.

### 7.0 — Deployability

Moved from 11/17 to **17/17 declared** this session, with the port map resolved and gated. But
declaring is not running: **no engine or agent container exists on the box**, and per-service venvs
exist for only 4 of ~20 supervised services. The platform is cold.

---

## Recommendations, ranked by value

| # | Action | Lifts | Effort |
|---|---|---|---|
| ~~1~~ | ~~Authenticate the 11 open clinical endpoints~~ | **DONE** — 4.0 → 7.5 | — |
| ~~2~~ | ~~Stop emitting `X-HCLS-Governed` unless gated~~ | **DONE** | — |
| 1 | **Set `HCLS_API_KEY` in production** — the gate ships off | Security 7.5 → 9 | trivial |
| 3 | **Put the working CUDA build behind the platform venv and the other services** — it already exists in one service venv | Hardware 5.0 → 8 | low |
| 4 | Adopt `install_governance` + `require_valid_input` across all 12 | Security, honesty | medium |
| 5 | Bring the platform up and keep it up | Deployability 7 → 9 | medium |
| 6 | Add the decision-support frame to genomic-foundation and precision-intelligence | Docs, honesty | low |
| 7 | Choose one runtime model (compose) | Layout, ops | decision |
| 8 | Rename the 12 stdlib-shadowing modules | Code quality | medium |
| 9 | Agree a test-depth floor (range is 4 → 1,966) | Confidence | medium |


**The two struck-through items were the highest-value changes in this document** — the platform's
credibility rests on honesty-by-construction, and it was shipping a header claiming governance it
did not perform, in front of unauthenticated clinical endpoints. Both are now fixed.

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
