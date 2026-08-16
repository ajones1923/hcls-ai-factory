# Gated Software — PRD

**Purpose.** Everything in the HCLS AI Factory that **cannot be installed without credentials,
a licence, or an accepted model licence** — what it is, what breaks without it, and what "done"
looks like. This is the document to work from during the final upgrade pass.

**Date:** 2026-08-15 · **Target:** NVIDIA DGX Spark — GB10, 128 GB unified LPDDR5x, 20 Grace cores,
aarch64 · **Companion:** `GATED_SOFTWARE_BUILD_GUIDE.md` (step-by-step install)

---

## 1. Why this document exists

The platform's tests pass (**8,402 green across 17 subjects**) because they exercise logic, not
models. The moment a demo needs a real structure prediction, a real docking pose, a real
secondary-analysis run or a real generated molecule, it hits software that is **not on this box and
cannot be fetched anonymously**.

The single most consequential fact:

```
torch 2.10.0+cpu        torch.cuda.is_available() = False
```

A CPU-only PyTorch on a machine whose entire premise is the GB10. Measured GPU utilisation across a
3-hour working session: **0%, 12.3–15.1 W — idle throughout.**

## 2. Inventory

| # | Component | Gate | Blocks | Priority |
|---|---|---|---|---|
| **G1** | **CUDA PyTorch** (aarch64 / sm_121) | NVIDIA channel | ESM-2, ESMFold, ProteinMPNN, LoRA, every local model | **P0** |
| **G2** | **NVIDIA Parabricks** | NGC, licensed container | Secondary analysis — BWA-MEM2 + DeepVariant (**94 code refs**) | **P0** |
| **G3** | **BioNeMo MolMIM NIM** | NGC NIM registry | Molecule generation (**131 refs**) | **P0** |
| **G4** | **BioNeMo DiffDock NIM** | NGC NIM registry | Docking poses (**121 refs**) | **P0** |
| G5 | ESMFold weights | HuggingFace licence acceptance | Structure prediction `:8570` | P1 |
| G6 | Chai-1 co-folding | Chai Discovery gated weights | Co-folding; registry status `planned` | P1 |
| G7 | RFdiffusion + DGL | GitHub + gated weights; **needs G1** | De novo backbone design (vendored) | P2 |
| G8 | MHCflurry models | PyPI + model download | Immunogenicity `:8577` (`planned`) | P2 |

### Not gated — simply missing

These need no credentials and should be installed first, because they unblock real work immediately:

```
scanpy  anndata  pysam  cyvcf2  rdkit  nextflow
```

`nextflow` is the orchestrator named throughout the architecture and is **not installed**.

## 3. What each gate actually blocks

**G1 — CUDA PyTorch.** The keystone. Without it: ESM-2 embeddings run on Grace cores (slow but
functional), ESMFold is impractical, ProteinMPNN and RFdiffusion are non-functional, and no demo can
honestly claim GPU acceleration. It also blocks G7, which needs a CUDA build to compile DGL.

**G2 — Parabricks.** The "Patient DNA → drug candidates in under five hours" claim rests on
GPU-accelerated alignment and variant calling. Without it the genomics engine can serve
pre-computed results only, which must be labelled **REPRESENTATIVE**, never LIVE.

**G3/G4 — MolMIM and DiffDock.** The therapeutic-discovery engine's two headline capabilities. With
252 combined code references, these are not peripheral. Without them the drug-discovery demo has no
honest content.

**G5 — ESMFold.** Structure prediction at `:8570`. Note the engine capability
`structural-biology-engine:8581` is now `planned` because **no process binds 8581** — the engine is
five separate services with no aggregator.

**G6 — Chai-1.** Already `planned` in the registry. Correctly labelled; do not promote until
weights are in place.

## 4. Requirements

| # | Requirement | Done when |
|---|---|---|
| GR-1 | NGC account with the correct entitlements | `docker login nvcr.io` succeeds; Parabricks + BioNeMo images pull |
| GR-2 | CUDA PyTorch installed | `torch.cuda.is_available()` is `True` **and** a real workload shows >0% GPU in `nvidia-smi` |
| GR-3 | Parabricks runs a real sample | FASTQ → BAM → VCF completes on HG002 with recorded wall-clock |
| GR-4 | MolMIM + DiffDock serve | `/health` 200 on both; one real generation and one real pose recorded |
| GR-5 | ESMFold serves | `:8570` returns a structure for a known sequence |
| GR-6 | Registry reflects reality | every `live` capability answers a health probe; anything unmet stays `planned` |
| GR-7 | Demos labelled honestly | each of the 17 marked LIVE / REPRESENTATIVE / BURST, matching what actually ran |

**GR-6 is the one that protects the project.** Two capabilities were found registered `live` while
unreachable (`singlecell-compute`, `structural-biology-engine`). Installing gated software must not
become an excuse to mark things `live` ahead of evidence.

## 5. Verification — mandatory before any status is promoted

```bash
# 1. the GPU is genuinely in use, not merely present
nvidia-smi --query-gpu=utilization.gpu,power.draw --format=csv -l 1
# run a real workload in parallel; utilisation must leave 0% and power must rise above ~15 W idle

# 2. every live capability answers
.venv/bin/python - <<'PY'
import json, urllib.request
for c in json.load(open("lib/hcls_common/capabilities.json"))["capabilities"]:
    if c.get("status")!="live" or not c.get("endpoint"): continue
    try: code=urllib.request.urlopen(f"http://{c['endpoint']}/health", timeout=3).status
    except Exception as e: code=type(e).__name__
    print(f"{c['id']:42s} {code}")
PY

# 3. the platform still passes its own gate
.venv/bin/python scripts/run_all_tests.py
.venv/bin/python scripts/validate_registry.py
```

## 6. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| No CUDA wheel for GB10 / sm_121 aarch64 | G1 fails; G7 unreachable | Use NGC PyTorch **containers**, which ship their own matched runtime, rather than pip wheels |
| Parabricks licence does not cover this box | G2 blocked | Confirm entitlement before planning demos around it |
| NIM images are x86-only | G3/G4 blocked locally | Elastic burst to remote GPUs over the private mesh — and say "elastic burst", never "all on one box" |
| Gated install slips, demos ship anyway | Honesty-gate breach | GR-6/GR-7 — status follows evidence, not intent |

## 7. Sequence

```
ungated installs (scanpy, rdkit, pysam, cyvcf2, nextflow)   ← do this first, no credentials needed
        │
GR-1 NGC entitlement ──► G1 CUDA torch ──► G5 ESMFold ──► G7 RFdiffusion/DGL
        │                      │
        └──► G2 Parabricks     └──► G3 MolMIM · G4 DiffDock ──► 17 demos
```

## 8. Open questions for Adam

1. **Does the NGC account carry Parabricks and BioNeMo entitlements**, or is a separate licence
   needed?
2. **Local or burst for the NIMs?** If MolMIM/DiffDock images are x86-only, the demos must say
   "elastic burst" — this changes the narration of the drug-discovery story.
3. **Chai-1 weights** — has access been granted? It is `planned` in the registry and should stay
   there until it has.
