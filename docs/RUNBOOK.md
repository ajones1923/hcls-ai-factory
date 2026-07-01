# HCLS AI Factory — Failure-Mode Runbook (A8)

Reproducible-environment failure modes and their fixes, captured from real builds. Supports
21 CFR Part 11 reproducibility: a known failure has a documented, repeatable remedy.

## Hardware / CUDA (the Blackwell GB10 wall)
| Symptom | Cause | Fix |
|---|---|---|
| Model pins `torch 1.9 / CUDA 11.1` (e.g. RFdiffusion) won't install/run | CUDA 11.1 maxes at Ampere `sm_86`; GB10 is Blackwell `sm_121` (needs CUDA 13) | Run on RunPod (x86 + standard CUDA), or use a modern-torch equivalent. See the RunPod plan. |
| `pip install parasail` → "error: building wheel failed" | No aarch64 wheel; source build needs the parasail C lib | Substitute Biopython `PairwiseAligner` (pure, ARM-clean) for exact SW. |
| `rapids-singlecell` installs but won't import | No CUDA-13 / Blackwell / aarch64 RAPIDS wheels | Keep scanpy (CPU); run RAPIDS on RunPod x86. |
| GPU model needs the GB10 but is on CPU | service started from the CPU venv | Start from the cu130 venv (`core/engines/therapeutic-discovery/venv`, torch 2.12+cu130); confirm `torch.cuda.is_available()`. |

## Service lifecycle
| Symptom | Cause | Fix |
|---|---|---|
| FastAPI returns 422 "field required" on a valid body | `from __future__ import annotations` makes the request-model annotation a string FastAPI can't resolve when the model is defined inside `create_app()` | Remove `from __future__ import annotations` from the `*_service.py` module. |
| Killing a service exits the shell (exit 144) | `pkill -f 'mod:_app_factory'` matches the bash subshell too | Kill by **port**: `ss -tlnp \| grep ":PORT " \| grep -oP 'pid=\K[0-9]+' \| xargs kill`. |
| Service won't start, no error | uvicorn factory path wrong | `python -m uvicorn module:_app_factory --factory --host 127.0.0.1 --port N`; check `*.log`. |

## Python env / dependencies
| Symptom | Cause | Fix |
|---|---|---|
| `ImportError: cannot import name 'DisjunctiveConstraint'` | transformers 5.x removed it (SAFE/older libs expect 4.x) | Isolate the library in its own venv, or use a graceful fallback (e.g. BRICS for generation). Do not downgrade transformers in the shared venv. |
| `ModuleNotFoundError: pydantic_settings` when importing `hcls_common.mlops` from a non-lib venv | the package `__init__` eager-imports the full stack | Load the stdlib-only module directly via `importlib.util.spec_from_file_location` (bypass the package init). |
| `mhcflurry` "Missing downloadable file" | weights not fetched / wrong invocation | Use the console script `mhcflurry-downloads fetch models_class1_presentation`, not `python -m mhcflurry.downloads`. |
| Geneformer install "finishes" but won't import | build/subprocess error; needs git-lfs + pinned deps | Dedicated setup (git-lfs, pinned torch/transformers), or run on RunPod x86. |

## Correctness
| Symptom | Cause | Fix |
|---|---|---|
| Ts/Tv ≈ 0.9 (should be ~2.0) | computed over gVCF including RefCall sites | Filter to `filter='PASS'` biallelic SNVs before Ts/Tv. |
| Local SW shows 100% identity for divergent sequences | local alignment covers only the conserved core | Use the SW **score** (not local %identity) as the discriminator; guard the empty-alignment case. |

## Executing external/third-party code
- The auto-mode classifier blocks running agent-cloned external repos (e.g. ProteinMPNN, RFdiffusion).
  This is intended — get explicit approval before executing third-party model code, then add a Bash
  permission rule if it should be standing.
