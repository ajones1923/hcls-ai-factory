# Best-Practices & Factual-Accuracy Audit

**Scope:** 8 Engines · 8 Intelligence Agents · 1 Disease Program (TSC).
**Date:** 2026-08-15 · **Method:** every figure below was produced by running a check against the
tree, not by reading it. Where a first result was wrong, the correction is stated.

---

## Summary

| Dimension | Result |
|---|---|
| Secrets in tracked code | **clean** — 0 real secrets |
| Container hardening (service images) | 16/18 → **18/18** after fixes |
| API authentication | 1 of 12 → **12 of 12** *(fixed 2026-08-16)* |
| Platform governance gate adoption | 1 of 12 → **12 of 12** *(fixed 2026-08-16)* |
| Gene → locus claims | **0 mismatches** |
| Pharmacogenomics (CPIC) | **44 pairs, all genuine** |
| Decision-support disclaimer in docs | **12 of 17** subjects |
| Test health | 17/17 green, 8,402 passing |

The clinical *content* is markedly stronger than the platform *controls*. The science checks out; the
service layer does not enforce the guarantees the project claims.

---

## 1. Security

### 1.1 Secrets — clean

A scan for hardcoded credentials across `core/`, `lib/`, `scripts/`, `hcls-orchestrator/` returned
**8 hits, all benign**: docstring examples (`api_key="your_ncbi_key"`, `anthropic_api_key="sk-ant-..."`)
and `api_key="not-needed"` placeholders for **local NIM / vLLM endpoints that genuinely take no
auth**. No live credential is committed. `.gitignore` covers `.env`, `.env.*` and `*.env`.

### 1.2 Authentication — FIXED 2026-08-16

**Was: 1 of 12 FastAPI entrypoints enforced authentication** (the TSC program). The other eleven —
every intelligence agent, plus the cardiology, imaging and oncology engines — exposed their
endpoints with no auth dependency at all.

**Now: 12 of 12.** `lib/hcls_common/api_auth.py` generalises the pattern TSC already used, so there
is one implementation rather than twelve. Posture is deliberately unchanged by default:

    HCLS_API_KEY unset  ->  open (the trusted-network / synthetic-demo posture the platform ships)
    HCLS_API_KEY set    ->  FAIL CLOSED on every route except /health, /docs, /openapi.json

Per-service override: `HCLS_API_KEY_<SLUG>`. Key comparison is `hmac.compare_digest` — a plain `==`
leaks key length and prefix through timing.

**Proven on a live service**, not just in unit tests (cart agent, `HCLS_API_KEY` set):

| Request | Result |
|---|---|
| `/query` no key | **401** |
| `/query` wrong key | **401** |
| `/query` correct key | 405 — passed the gate and reached routing |
| `/health` | 200 — probes must not fail closed |

These are clinical decision-support endpoints. They accept patient context and return variant
interpretations, therapy suggestions and trial matches.

> **Correction made during this audit.** I first counted 2/12 as "authenticated" by matching
> `api_gate`. That is wrong: `lib/hcls_common/api_gate.py::install_governance` is **not**
> authentication. It attaches a request ID, timing, and an `X-HCLS-Governed` response header, plus a
> `/governance` info endpoint. Its own payload describes the gates as things a handler must *call* —
> `require_valid_input()` and `honesty_flags()` — not middleware enforcement.

That produced a second, sharper issue: **a response header asserting `X-HCLS-Governed` on a request
that was not gated** — instrumentation dressed as assurance.

**FIXED 2026-08-16.** The middleware now emits `X-HCLS-Service` (true — it does add the service
name, a request id and timing) and emits `X-HCLS-Governed` **only when a gate actually executed**,
listing which ones. `require_valid_input()` and `honesty_flags()` record themselves into a
per-request set that the middleware reads.

A subtlety worth keeping: the first implementation rebound a `ContextVar` inside the handler, and
the header never appeared — Starlette runs sync endpoints in a threadpool with a copied context, so
the middleware never saw the change. It now mutates one shared set. **A test asserts the header
actually appears**, which is what caught it.

### 1.3 Governance gate adoption — 1 of 12

| Control | Adoption |
|---|---|
| `install_governance` / `create_governed_app` | 1/12 (cart) |
| `require_valid_input()` | 1/12 (cart) |
| `honesty_flags()` / `assert_publishable()` | 1/12 (cart) |

`core/agents/cart` is the reference implementation and **nothing else adopted it**. This confirms
the outstanding "gate-middleware" remediation item with numbers.

### 1.4 Container hardening — fixed

27 Dockerfiles. They split into two populations, and treating them alike gives a false reading:

- **18 long-running service images** — must run unprivileged with a healthcheck.
- **9 MONAI Application Packages** (`maps/*`) — batch jobs executed by `monai-deploy` with
  `--gpus all` over mounted `/input` `/output`. Root is the MONAI base-image convention and a
  HEALTHCHECK on a batch container is meaningless. **Not findings.**

Of the 18 service images, two failed:

| Image | Problem | Action |
|---|---|---|
| `core/disease-programs/tuberous-sclerosis` | ran as **root**, no HEALTHCHECK | non-root `tscuser` + healthcheck added |
| `core/engines/structural-biology` | no HEALTHCHECK | added *(this image was written earlier in this same session — my own omission)* |

All 27 pin an explicit base-image tag; none use `:latest`. TSC's `COPY . .` is broad but its
`.dockerignore` excludes `.env`, `venv/`, `.venv/`, `data/` and `.git/`.

---

## 2. Design

### 2.1 Standard-library shadowing — 12 files

`src/collections.py` in 11 subjects, plus `precision-autoimmune/config/logging.py`. Imports are
package-qualified (`from src.collections import …`) so the code works today, but any configuration
that places `src/` on `sys.path` — a routine Docker and pytest pattern — kills the interpreter
before collection:

```
ImportError: cannot import name 'namedtuple' from partially initialized module 'collections'
```

This is not hypothetical: it produced **two false audit results** during this session before being
identified. Rename to `vector_collections.py` (PRD R9).

### 2.2 Test depth is very uneven

4 tests (single-cell engine) to 1,966 (cardiology). Pass rate is uniform and therefore uninformative;
coverage is the real signal. Subjects that emit clinical output need an agreed floor.

### 2.3 Two half-maintained deployment models

`docker-compose.dgx-spark.yml` (now 17/17) versus `health-monitor.sh` under cron. The compose header
calls itself "the declarative/portable target" and says the live system is host-process based — but
the per-service `venv/` directories that path needs **do not exist**, so it cannot start anything.
Unresolved (PRD open decision 2).

---

## 3. Factual accuracy

This is the strongest part of the platform.

### 3.1 Gene → locus — 0 mismatches

Every gene/cytoband assertion found in code was checked against HGNC/NCBI: APOE 19q13.32,
CYP2C19 10q23.33, LDLR 19p13.2, LMNA 1q22, MYH7 14q11.2, TTN 2q31.2. **All correct.**

### 3.2 Pharmacogenomics — 44 pairs, all genuine CPIC

13 genes with correctly paired drugs:

| Gene | Drugs asserted |
|---|---|
| CYP2D6 | codeine, tramadol, tamoxifen, ondansetron, atomoxetine, nortriptyline, amitriptyline, paroxetine, fluvoxamine |
| CYP2C19 | clopidogrel, voriconazole, omeprazole, escitalopram, citalopram, sertraline, amitriptyline |
| DPYD | fluorouracil, capecitabine, 5-FU, tegafur |
| TPMT / NUDT15 | azathioprine, mercaptopurine, thioguanine |
| SLCO1B1 | simvastatin, atorvastatin, rosuvastatin |
| VKORC1 / CYP2C9 | warfarin (+ phenytoin, celecoxib) |
| HLA-B | abacavir, carbamazepine, allopurinol, phenytoin |
| UGT1A1 | irinotecan, atazanavir |
| G6PD | rasburicase, primaquine |
| IFNL3 | peginterferon |
| CFTR | ivacaftor |

**No spurious pair was found.** Variant identifiers spot-checked and correct: VKORC1 **rs9923231**
(−1639G>A), SLCO1B1 **rs4149056** (c.521T>C, \*5). CYP3A5→tacrolimus is genuine CPIC; CYP3A4→
midazolam/cyclosporine are pharmacologically correct substrate relationships rather than CPIC
guideline pairs — accurate, but should not be presented as guideline-backed.

### 3.3 Decision-support framing — 12 of 17

Five subjects **never** state the limit in any documentation:

```
genomic-foundation · precision-intelligence · therapeutic-discovery
structural-biology · single-cell
```

All five are **engines** (compute layer) rather than agents, which is a partial defence — but two of
them are not merely infrastructure. `genomic-foundation` emits variant calls and ACMG secondary
findings; `precision-intelligence` performs variant interpretation over ClinVar/AlphaMissense. Both
reach clinical interpretation and should carry the frame.

Coverage where it matters most is good: TSC 70 of 145 documents, clinical-imaging 21 of 48,
precision-oncology 19 of 21.

---

## 4. Prioritised recommendations

| # | Recommendation | Severity |
|---|---|---|
| ~~A1~~ | ~~Enforce authentication on the 11 unauthenticated endpoints~~ | **DONE** |
| ~~A2~~ | ~~Stop emitting `X-HCLS-Governed` unless the request was gated~~ | **DONE** |
| ~~A3~~ | ~~Adopt `install_governance` across all 12 entrypoints~~ | **DONE** (11/12; TSC keeps its own) |
| A3b | Call `require_valid_input()` / `honesty_flags()` in the handlers themselves | high |
| A4 | Add the decision-support frame to genomic-foundation and precision-intelligence | high |
| A5 | Rename the 12 stdlib-shadowing modules | medium |
| A6 | Agree a test-depth floor for clinical-output subjects | medium |
| A7 | Choose one deployment model | medium |
| A8 | Label CYP3A4 pairs as substrate relationships, not CPIC | low |

A1–A3 are one change in shape: the platform already has the gate; it is adopted by one service.

## 5. What was fixed during this audit

- TSC container: root → unprivileged `tscuser`, HEALTHCHECK added.
- structural-biology container: HEALTHCHECK added.

## 6. Errors made during this audit

1. Counted `api_gate` as authentication. It is observability, not enforcement — corrected in §1.2,
   and the correction produced a better finding than the original claim.
2. Initially flagged all 27 Dockerfiles against service criteria. Nine are MONAI batch packages
   where root and no-healthcheck are correct; separating them cut 10 false findings to 2 real ones.
