# Layout Review — 8 Engines · 8 Agents · 1 Disease Program

**Date:** 2026-08-15 · **Method:** structural inventory of all 17 subject directories, plus the
platform's deployment and registry surfaces.

---

## 1. Structural conformance

Twelve of seventeen subjects follow the same eleven-element template; five do not. The template is
`api/ app/ src/ tests/ config/ docs/ scripts/ data/ README.md requirements.txt Dockerfile`.

| Element | Adoption |
|---|---|
| `src/`, `tests/`, `docs/`, `README.md` | **17/17** |
| `data/`, `requirements.txt`, `Dockerfile` | 16/17 |
| `scripts/` | 15/17 |
| `app/`, `config/` | 14/17 |
| `api/` | **12/17** |

**The five non-conformant subjects are exactly the five that were also missing from
`docker-compose.dgx-spark.yml`:**

| Subject | Missing |
|---|---|
| genomic-foundation | `api/`, `app/`; Dockerfile + requirements live under `web-portal/` |
| precision-intelligence | `api/` |
| therapeutic-discovery | `api/`, `config/` |
| structural-biology | `api/`, `app/`, `config/`, `scripts/`, `data/` |
| single-cell (engine) | `api/`, `app/`, `config/`, `scripts/` |

That correlation is the finding: **conformance and deployability moved together.** The subjects
built to the template were declared, containerised and supervised; the five that diverged were
invisible to the launcher until this session.

> I checked whether the five shared an origin date, expecting a "pre-template cohort". They do not —
> all but genomic-foundation were first committed on the same day as the conformant subjects. The
> divergence is architectural, not chronological, and the theory was dropped.

## 2. Architectural shape — three genuine patterns

Not every divergence is a defect. Three distinct shapes exist and only one is wrong:

**A · Service subjects (12).** One API + one UI, template-conformant. Correct as-is.

**B · Portal engines (2).** `genomic-foundation` and `precision-intelligence` are Flask portals on
5000/5001 rather than API+UI pairs. They have no `api/` because they do not need one. Their layout
is *different*, not *deficient* — though genomic-foundation nesting its Dockerfile and
`requirements.txt` under `web-portal/` while its `src/`, `tests/` and `config/` sit at the top level
is genuinely confusing and worth flattening.

**C · Multi-service engines (2).** `structural-biology` is **five** FastAPI services (esmfold 8570,
protein-search 8571, developability 8576, immunogenicity 8577, proteinmpnn 8578) sharing one
codebase. `single-cell` is one compute service on 8573. Neither fits the API+UI template.

Pattern C exposed a real defect: the registry advertised `structural-biology-engine` at **:8581 with
no process binding it**. There is no aggregator. Now marked `planned`.

## 3. Deployment surface — resolved this session

| | Before | After |
|---|---|---|
| Declared in compose | 11/17 | **17/17** |
| Port convention | none; 8 of 13 registry/supervisor disagreements | **UI / UI+1, CI-gated** |
| Capabilities `live` but unreachable | 2 | **0** |

Added: `genomic-foundation-engine`, `precision-intelligence-engine`, `therapeutic-discovery-engine`,
`singlecell-compute`, `tuberous-sclerosis-program`, `structural-biology-search`, plus Dockerfiles for
single-cell and structural-biology.

## 4. Runtime model — still two, and it should be one

| | Compose | `health-monitor.sh` (cron, 5-min) |
|---|---|---|
| Coverage | 17/17 subjects | ~20 services |
| Interpreter | image-internal | `./venv/bin/python` |
| Venvs present | n/a | **4 of ~20** — precision-intelligence, therapeutic-discovery, tuberous-sclerosis, genomic-foundation/web-portal |
| Agents | all declared | **none have a venv** |

`health-monitor.sh::rebuild_venv()` can create the missing venvs on demand, so the host path is
recoverable rather than broken. But maintaining both models is why the port map drifted in the first
place: two sources of truth, neither authoritative.

**Recommendation: standardise on compose**, and reduce `health-monitor.sh` to container health
checks. This is PRD open decision 2 and needs Adam.

## 5. Recommendations

| # | Recommendation | Effort |
|---|---|---|
| L1 | Choose one runtime model (compose) | decision |
| L2 | Flatten `genomic-foundation` — move `web-portal/Dockerfile` + `requirements.txt` to the subject root | low |
| L3 | Build an aggregator on :8581 **or** keep `structural-biology-engine` `planned` | medium |
| L4 | Document patterns A/B/C so divergence is deliberate, not accidental | low |
| L5 | Add `api/` to therapeutic-discovery and single-cell if they are to expose programmatic surfaces | medium |
| L6 | Rename the 12 stdlib-shadowing modules (`src/collections.py`) | medium |

L1 is the one that matters. The rest are tidiness; L1 is why things drift.

## 6. What is already optimal

- `src/` + `tests/` + `docs/` + `README.md` at **17/17** — no subject is undocumented or untested.
- Base-image tags pinned in **27/27** Dockerfiles; none use `:latest`.
- Registry coverage complete: every engine and agent directory maps to a registered capability, and
  `validate_registry.py` fails if a new directory appears without one.
- 8,402 tests green across all 17.
