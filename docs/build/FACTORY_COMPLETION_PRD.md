# PRD — Finishing the HCLS AI Factory

**Date:** 2026-08-17 · **Owner:** Adam Jones · **Companion:** [FACTORY_COMPLETION_WORKBOOK.md](FACTORY_COMPLETION_WORKBOOK.md)
**Status of this document:** repo-only. It is not published to hcls-ai-factory.org.

---

## 1. The problem in one paragraph

The HCLS AI Factory is **built but not running**. Forty-two capabilities are registered, 8,402 tests
pass across all seventeen subjects, the merge gate is green, and the public site is accurate. But on
the reference machine **two of roughly twenty services are up**, there is **no `.env` file at all**,
the vector database on the expected port belongs to a different stack, and **7 of 17 demonstrations
can run**. The distance between "the code is real" and "a visitor can watch it work" is the whole of
this document.

Nothing here is a rewrite. The gap is bring-up, credentials, gated software, and seeding — in that
order.

---

## 2. Measured starting state (2026-08-17)

Every figure below was produced by running a check, not by reading a doc.

| Dimension | State |
|---|---|
| Registry | 42 capabilities — **25 live · 9 verified · 7 planned · 1 gated** |
| Tests | 8,402 passing, 0 failing, 17/17 subjects |
| Services responding | **2** (`:5000` genomics, `:5001` precision-intelligence) |
| `.env` | **absent entirely**; 0 of 4 required values set |
| Milvus `:19530` | occupied by `imaging-milvus-standalone`, **12 imaging collections only** |
| Demos with prerequisites met | **7 of 17** |
| Per-service venvs | 3 of ~20 supervised services |
| HCLS engine/agent container images | **0** built (149 unrelated images on the box) |
| Compose services declared | 26 |
| Python packages missing from platform venv | 13 of 25 checked (`pysam`, `cyvcf2`, `rdkit`, `monai`, `SimpleITK`, `e3nn`, `openai`, `langgraph`, `minio`, `redis`, `pynvml`, `plotly`, `fhir.resources`) |
| System binaries missing | **8 of 13** (`samtools`, `bcftools`, `tabix`, `bgzip`, `nextflow`, `gatk`, `ffmpeg`, `pbrun`) |
| CUDA | works — `torch 2.12.1+cu130`, GB10 at `sm_121` — but **only in `core/engines/therapeutic-discovery/venv`** |
| GPU allocation | **fails**: 99 GiB of 119 GiB unified memory held by the host page cache |
| Reference data | ClinVar + AlphaMissense present (1.5 GB `hcls-ai-factory-core-data/`); **GRCh38 absent**; HG002 FASTQ are **0-byte placeholders** |
| Security | `HCLS_API_KEY` unset — all 12 entrypoints open |
| Governance gates | middleware 11/12 · input gate **1/12** · honesty gate **1/12** |
| CI | green overall, but the subject-suite job is **failing** under `continue-on-error: true` |

---

## 3. What "finished" means

Five acceptance criteria. Each is a command whose output settles it — no judgement calls.

| # | Done when | Command that proves it |
|---|---|---|
| **A1** | Every registered `live` capability answers a health probe | `scripts/validate_registry.py --probe` (to be added, R11) |
| **A2** | All 17 demonstrations run start to finish and write a transcript | `run_demo.py --check-all` reports **17/17**; `run_demo.py --all` exits 0 |
| **A3** | The platform survives a reboot unattended | reboot, wait 5 min, `--check-all` still **17/17** |
| **A4** | A stranger reproduces the quickstart on a clean clone | `run_all_tests.py` → 17 subjects, 0 failed; `run_demo.py E8` → PASS |
| **A5** | Nothing on the site claims more than the box delivers | `mkdocs build --strict` green **and** every `live` badge backed by A1 |

**The honesty rule stands throughout:** a capability is promoted only after it answers a probe; a
demo is labelled LIVE only after it has run on real input. If a phase cannot reach its criterion,
the registry status changes — not the criterion.

---

## 4. Requirements

Priority: **P0** blocks everything · **P1** blocks a phase · **P2** quality.

### 4.1 Foundation (no accounts needed)

| # | Requirement | Pri | Done when |
|---|---|---|---|
| R1 | `.env` exists from `.env.example` with all four required values | **P0** | `docker compose config` resolves with no warnings |
| R2 | `HCLS_API_KEY` set; auth verified fail-closed | **P0** | 401 without key, 401 wrong key, 200/405 with key |
| R3 | 13 missing Python packages installed into the platform venv | P1 | `import` check passes for all 25 |
| R4 | 6 free system binaries installed (`samtools` `bcftools` `tabix`/`bgzip` `nextflow` `gatk` `ffmpeg`) | P1 | each on `PATH`, `--version` responds |
| R5 | Page-cache drop added to bring-up so CUDA can allocate | **P0** | `torch.cuda.mem_get_info()` returns without error |
| R6 | A dedicated Milvus for HCLS that does not collide with the imaging stack | **P0** | HCLS connects and lists its own collections |

### 4.2 Bring-up

| # | Requirement | Pri | Done when |
|---|---|---|---|
| R7 | One runtime model chosen — compose **or** host processes, not both | **P0** | decision recorded in this PRD §7 |
| R8 | All 26 compose services build and start | P1 | `docker compose ps` shows all healthy |
| R9 | Agent corpora seeded for all eight agents | **P0** | each agent's collections exist and return hits |
| R10 | Services restart automatically after reboot | P1 | A3 passes |
| R11 | `validate_registry.py --probe` added — fails CI on a `live` capability that does not answer | P1 | flipping a service off turns CI red |

### 4.3 Gated software

| # | Requirement | Pri | Done when |
|---|---|---|---|
| R12 | NGC account + `docker login nvcr.io` succeeds | **P0** for this phase | *Login Succeeded* |
| R13 | CUDA PyTorch behind the platform venv and every GPU service | P1 | `cuda.is_available()` true from the platform interpreter |
| R14 | Parabricks pulled; FASTQ→VCF runs on real HG002 | P1 | VCF written, Ts/Tv ≈ 2.0 |
| R15 | MolMIM + DiffDock NIMs deployed **or** the arch answer recorded and demos relabelled BURST | P1 | either containers healthy, or §7 decision + relabel |
| R16 | Real HG002 FASTQ downloaded (~200 GB) and GRCh38 obtained | P1 | files non-zero, index built |

### 4.4 Demonstrations

| # | Requirement | Pri | Done when |
|---|---|---|---|
| R17 | 17/17 demos have prerequisites met | **P0** | `--check-all` reports 17/17 |
| R18 | Every demo writes a diffable transcript | P1 | `demo/transcripts/*.txt` for all 17 |
| R19 | Each demo's label matches what actually ran | **P0** | LIVE only where the service answered |
| R20 | Demo catalogue counts regenerated from the runner, not hand-written | P2 | catalogue figures come from `--check-all` |

### 4.5 Quality debt

| # | Requirement | Pri | Done when |
|---|---|---|---|
| R21 | Subject-suite CI job made blocking and green | P1 | `continue-on-error` removed, run green |
| R22 | `require_valid_input()` / `honesty_flags()` called by handlers, not just available | P1 | coverage 12/12 on the governance page |
| R23 | 12 stdlib-shadowing modules renamed (`src/collections.py` → `vector_collections.py`) | P2 | bare `pytest` works without the harness |
| R24 | structural-biology aggregator on `:8581`, or the engine stays `planned` | P2 | probe answers, or status unchanged |
| R25 | Test-depth floor agreed for clinical-output subjects | P2 | recorded; range is currently 4 → 1,966 |

---

## 5. Phases

Each phase has an exit gate. **Do not start the next phase until the gate passes** — the reason is
in §8, and it is the difference between debugging one thing and debugging four.

```
Phase 0  Foundation        R1-R6        no accounts        ~1 day
   │     gate: .env resolves, API key fails closed, CUDA allocates, HCLS Milvus up
   ▼
Phase 1  Bring-up          R7-R11       free               ~2-4 days
   │     gate: all services healthy, corpora seeded, 13+ demos ready, survives reboot
   ▼
Phase 2  Free accounts     B-tier       HuggingFace, NCBI, Anthropic   ~1 day
   │     gate: BGE embeds, Claude answers, agents return cited output
   ▼
Phase 3  Gated software    R12-R16      NGC entitlement    unknown — the long pole
   │     gate: Parabricks produces a real VCF; NIM arch question answered
   ▼
Phase 4  Demonstrations    R17-R20      free               ~3-5 days
   │     gate: 17/17 run, transcripts recorded, labels honest
   ▼
Phase 5  Quality debt      R21-R25      free               ongoing
         gate: CI blocking and green; governance coverage stated as a number
```

**Phases 0, 1, 2, 4 and 5 need no gated software.** That is roughly 80% of the remaining work and it
is all unblocked today. Phase 3 is the only part waiting on someone else.

---

## 6. Non-goals

Stating these prevents scope creep into work that does not move A1–A5.

- **Not** a multi-node or Kubernetes deployment. One box plus elastic burst is the thesis.
- **Not** clinical validation, regulatory submission, or a medical-device claim.
- **Not** new capabilities. Forty-two is enough; finish them.
- **Not** raising test counts for their own sake. Depth where output is clinical (R25), not everywhere.
- **Not** re-architecting the engine/agent split, the registry, or the port convention.
- **Not** publishing the internal build docs to the public site. They stay in the repo.

---

## 7. Open decisions — you must answer these

These are yours, not mine. Each blocks a requirement, and guessing wrong costs days.

| # | Decision | Blocks | Why it matters |
|---|---|---|---|
| **D1** | **Compose or host processes?** The compose file calls itself the portable target; `health-monitor.sh` runs the live system. Both are half-maintained. | R7, R8, R10 | Two runtime models means every service is configured twice and drifts. Pick one; delete or clearly demote the other. |
| **D2** | **Are MolMIM and DiffDock x86-only?** | R15 | If yes, they can never run locally, the flagship drug-discovery demo is permanently **BURST**, and every script must say "elastic burst". Determine before building demos around them. |
| **D3** | **Which Milvus does HCLS use?** `:19530` is the imaging stack's. Share it (namespaced collections) or stand up a second on another port. | R6, R9 | Sharing risks one stack's reset wiping the other's corpora. |
| **D4** | **AlphaMissense is CC BY-NC-SA — non-commercial.** It is the most-referenced external artefact in the codebase (455 refs) and the platform is Apache-2.0 with "commercial use welcome". | R16, site copy | A licence conflict at the centre of the value proposition. Either scope the claim, isolate the dependency, or accept and document. |
| **D5** | **Is 200 GB of HG002 FASTQ worth the disk?** The placeholders are 0 bytes today. | R14, R16 | A downsampled chromosome subset may prove the pipeline at a fraction of the cost. |
| **D6** | **Does the video build tree stay local-only?** `tmp/` is gitignored — every script and cached take lives on one disk. | — | Losing it loses the ability to edit your primary assets. |

---

## 8. Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| NGC entitlement never arrives | medium | Phase 3 blocked | Phases 0–2, 4–5 are independent. Relabel affected demos REPRESENTATIVE and ship. |
| MolMIM/DiffDock are x86-only | **likely** | flagship demo cannot run locally | Answer D2 early; write the burst path into the narration now, not later. |
| Page cache re-fills and CUDA fails again | high | intermittent GPU failure | R5 makes the drop part of bring-up, not a manual step. |
| Shared Milvus wiped by the other stack | medium | all agent corpora lost | D3; back up collections before any reset. |
| Bringing up four things at once | high | cannot tell which failed | The phase gates exist for exactly this. |
| Registry drifts from reality once services run | medium | the honesty claim breaks | R11 makes CI enforce it. |
| AlphaMissense licence surfaces publicly first | low | credibility | D4 now, on your terms. |

---

## 9. What this buys

When A1–A5 pass, the site's central claim — *fork it and run it yourself* — becomes literally true
for a stranger, the 17 demonstrations become watchable rather than described, and the maturity
matrix becomes a live readout rather than a periodic audit. That is the difference between a project
that is impressive to read and one that is impressive to *use*, which is the gap the current site
cannot close on its own.

---

*Work the [Workbook](FACTORY_COMPLETION_WORKBOOK.md) in order. It carries the exact commands, the
verification after each step, and the traps that have already cost time on this machine.*
