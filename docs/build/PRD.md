# PRD — Bringing the HCLS AI Factory to Full Operation

**Status:** draft for review · **Date:** 2026-08-15 · **Author:** generated during batch-01
**Companion documents:** [`GAP_ANALYSIS.md`](GAP_ANALYSIS.md) (evidence) ·
[`BUILD_GUIDE.md`](BUILD_GUIDE.md) (how) · `GATED_SOFTWARE_PRD.md` (what needs credentials)

---

## 1. Problem statement

The HCLS AI Factory has **8,402 passing tests across 385,000 lines** of engine, agent and disease-
program code — and **zero of it is running**. The work required is not feature development. It is
provisioning, deployment, and closing a small number of correctness defects that block operation.

A user who clones this repo today cannot start it: there is no environment template, half the
subjects are absent from the compose file, the shared library will not import, and two agents claim
the same two ports.

## 2. Goals

| # | Goal | Measure of done |
|---|---|---|
| P0-1 | Any clone can run the merge gate | `ruff` + `pytest` + `validate_registry.py` pass from a fresh checkout |
| P0-2 | Any clone can start the platform | `docker compose up -d` brings up all 17 subjects healthy |
| P0-3 | No capability is registered `live` that cannot be reached | health probe on every `live` endpoint returns 200 |
| P0-4 | The GPU is usable | `torch.cuda.is_available()` is `True`; a real workload shows >0% GPU |
| P1-1 | Ports are unambiguous and gated | one documented convention; CI fails on conflict |
| P1-2 | Every subject has a runnable demo | 17 demos, each executable start to finish |
| P2-1 | Test depth is proportional to risk | no clinical-output subject below an agreed floor |

## 3. Non-goals

- Rewriting working engines. Their tests pass; leave them.
- Making gated software un-gated. Parabricks, BioNeMo NIMs, Chai-1 and CUDA wheels need NGC
  credentials; this PRD *documents* their acquisition, it does not work around it.
- Claiming clinical validity. Everything remains decision support, never diagnosis.

## 4. Requirements

### P0 — the platform cannot start without these

**R1 · Reproducible Python environment.** *(done)*
`lib/hcls_common` was unimportable: a stale editable install with an empty `MAPPING`, unfixable
because PEP 668 blocks pip into system Python and no venv existed. A repo-level `.venv` is now the
supported interpreter. `start-factory.sh` still expects per-service `venv/` directories that do not
exist — either create them or move fully to containers (**decision required**).

**R2 · Environment template.** *(done)* `.env.example` added. Compose requires four variables and
failed to parse without them; there was no way to discover this but to read the YAML.

**R3 · Every subject declared in compose.** *(done for 5 of 6)*
Added `genomic-foundation-engine`, `precision-intelligence-engine`, `therapeutic-discovery-engine`,
`singlecell-compute`, `tuberous-sclerosis-program`. Compose now validates and declares 16 of 17.
**`structural-biology` still has no Dockerfile of its own** (the two found earlier belong to vendored
RFdiffusion) — R4.

**R4 · Containerise structural-biology.** Needs a Dockerfile + compose entry on its registry port
8579. Its own suite passes (34 tests) once `vendor_rfdiffusion` is excluded; the vendored tree needs
gated GPU packages and must not be built into the image.

**R5 · Resolve the neurology / precision-biomarker port collision.** *(done — Adam's decision,
2026-08-15)* Convention adopted: **registry endpoint = UI, API = UI + 1**. Applying it exposed two
further collisions (adjacent UI ports) and two clashes with running processes, so three services
were re-seated: imaging 8525→8523, neurology 8529→8535, structural-biology 8579→8581. 8534 is
returned to the NV-Segment-CT NIM. Full allocation: [`PORT_MAP.md`](PORT_MAP.md).

**R6 · CUDA-enabled PyTorch for GB10 / aarch64.** `torch 2.10.0+cpu`, `cuda_is_available() = False`
on a machine whose entire premise is the GB10. Every ML path is degraded or dead. Gated — see the
gated-software PRD.

### P1 — correctness and trustworthiness

**R7 · Make `singlecell-compute` genuinely live.** *(done)* It was registered `live` on `:8573` with
no Dockerfile, no requirements and no compose entry, so nothing could bind the port. The code was
always real (scanpy QC → normalise → HVG → PCA → neighbours → Leiden → marker DE, verified on
PBMC 3k). Added requirements, Dockerfile and compose entry.

**R8 · Gate the port map.** *(done)* `validate_registry.py` now enforces UI/UI+1 **and** parses
`health-monitor.sh`, failing on any supervised port the registry does not allocate. Both failure
modes were negative-tested rather than assumed.

**R9 · Remove stdlib shadowing.** 12 files named for stdlib modules, 11 of them `src/collections.py`.
Imports are package-qualified (`from src.collections import …`) so it works today, but any config
that puts `src/` on `sys.path` kills the interpreter before collection — this happened three times
during the audit. Rename to `vector_collections.py`.

**R10 · Health probes for every live capability.** A registry status of `live` should be verified by
a probe in CI, not asserted in JSON.

### P2 — depth

**R11 · Level up test depth.** Range is 4 tests (single-cell engine) to 1,966 (cardiology). Pass rate
is not the signal; coverage is. Agree a floor for any subject that emits clinical output.

**R12 · 17 demos.** One per engine, agent and the TSC program — each runnable end to end, each
labelled LIVE / REPRESENTATIVE / BURST.

**R13 · 85 documentation assets.** Overview, Foundation, Advanced and Demo guide per subject.

## 5. Sequencing

```
R1 R2 R3 R5 R7 R8 ──► R4 ──► deploy 17 ──► R9 R10 ──► R12 ──► R13
                                   │
                       R6 (gated, parallel — unblocks GPU demos)
```

R5 and R8 are **closed**. R4 (containerise structural-biology) is the last thing between here and a
full bring-up. R6 gates every GPU demo.

## 6. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| CUDA wheels unavailable for GB10/sm_121 | No GPU acceleration at all | Verify NGC channel early; fall back to NIM containers that ship their own runtime |
| Port rewiring breaks the running cardiac demo | Regression on approved work | Nothing is running now — do it before bring-up |
| `health-monitor.sh` fights manual starts | Services flap every 5 min | Disable the cron entry during bring-up |
| Registry drifts from reality again | Honesty gate erodes | R8 + R10 make it CI-enforced |

## 7. Open decisions for Adam

1. ~~R5 — port convention.~~ **Decided 2026-08-15:** registry = UI, API = UI + 1. Applied.
2. **Containers or host processes?** The compose header says agents run as native processes under
   `health-monitor.sh`; the compose file is described as "the declarative/portable target." Two
   deployment models are half-maintained. Pick one.
3. **Test-depth floor** for clinical-output subjects (R11).
