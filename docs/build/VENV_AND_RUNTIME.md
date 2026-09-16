# Service venvs and the runtime model

**Status:** operational note, repo-only. Not published to hcls-ai-factory.org.
**Written:** 2026-09-15, after taking the box from 9/24 to 24/24 healthy services.

`venv/` is gitignored, so none of what follows is visible in a diff. It is recorded here because
every one of these facts cost time to rediscover.

## The shape

`health-monitor.sh` starts each service with a **per-service `./venv/bin/…`**. There are three
different kinds of venv behind that one convention, and they are not interchangeable:

| Kind | Which services | What it is |
|---|---|---|
| **Symlink to the platform `.venv`** | the 8 agents, cardiology, precision-oncology/agent, rare-disease | `venv -> /home/adam/projects/hcls-ai-factory/.venv` |
| **Own venv, CPU** | precision-intelligence | real venv, its own interpreter |
| **Own venv, CUDA** ⚠️ | therapeutic-discovery | **torch 2.12.1+cu130, `cuda.is_available() == True`** |

⚠️ **Do not replace `core/engines/therapeutic-discovery/venv` with a symlink to `.venv`.** It is
the only interpreter on this box with working CUDA. The platform `.venv` carries torch
**2.10.0+cpu**. Measuring GPU capability through the wrong interpreter is how the scorecard's
"hardware utilisation 5.0" finding happened.

## Why the agents share one venv

13 of 21 services had no `venv/` at all, which is what kept the fleet down and the supervisor
crash-looping. Building 10 separate ones would have installed torch ten times (~30 GB) for dep
sets that differ trivially: across the 10 agent-class services there are 33 distinct packages and
exactly **one** true `==` conflict (tqdm 4.66.2 vs 4.66.4). The platform `.venv` already satisfied
every runtime import except `plotly`.

It is also the interpreter `scripts/run_all_tests.py` has been running all 17 suites against — 8,402
passing tests are the evidence that one shared environment works for this code.

So: install `plotly`, symlink, done. One place to patch, no duplication.

## Two traps already paid for

**1. Rename-severed console scripts.** `precision-intelligence` and `therapeutic-discovery` were
once `rag-chat-pipeline/` and `drug-discovery-pipeline/`. A venv's console scripts hard-code an
absolute interpreter path in their shebang, so after the rename all 141 of them pointed at
directories that no longer contain a Python:

```
$ head -1 core/engines/therapeutic-discovery/venv/bin/streamlit
#!/home/adam/projects/hcls-ai-factory/drug-discovery-pipeline/venv/bin/python3   # gone
```

The failure surfaces as `cannot execute: required file not found` — which reads like a missing
file, not a broken shebang. Repaired by rewriting line 1 of every script in those two `venv/bin/`
directories. **If either engine is moved again, repair the shebangs or the UIs die silently.**
(The empty `rag-chat-pipeline/` and `drug-discovery-pipeline/` directories still exist at the repo
root, root-owned, holding only stale config. They are not the engines.)

**2. A module is not a console script.** `.venv` could `import streamlit` fine while
`.venv/bin/streamlit` did not exist, so every Streamlit UI failed with
`No such file or directory` despite the package being installed. `pip install --force-reinstall
--no-deps streamlit` restores the entry point. Checking `importlib.util.find_spec` proves the
import works; it says nothing about the launcher.

## Rebuilding from scratch

```bash
.venv/bin/pip install plotly streamlit
for d in core/agents/*/ core/engines/cardiology core/engines/precision-oncology/agent; do
  ln -sfn "$PWD/.venv" "$d/venv"
done
./health-monitor.sh fix && ./health-monitor.sh status   # expect 24/24
```

Related: `docs/build/FACTORY_COMPLETION_WORKBOOK.md` (Phase 1), `docs/build/PORT_MAP.md`.

## Services the supervisor does not know about

`health-monitor.sh`'s SERVICES table covers 24 processes. Several registered capabilities are
**not in it** and must be started by hand, which is why they sat "registered `live`, nothing
listening" — the code was fine, nothing was ever asked to run it. All use a `create_app()`
factory (no `__main__` block), so they need uvicorn's `--factory`:

```bash
V=$PWD/.venv/bin/python
( cd core/engines/genomic-foundation/src   && $V -m uvicorn --factory variant_store_service:create_app --port 8575 & )
( cd core/engines/single-cell/src          && $V -m uvicorn --factory single_cell_service:create_app   --port 8573 & )
( cd core/engines/structural-biology/src   && $V -m uvicorn --factory proteinmpnn_service:create_app   --port 8578 & )
( cd core/engines/therapeutic-discovery/small-molecule/src \
                                           && $V -m uvicorn --factory molecule_gen_service:create_app --port 8574 & )
# the flagship disease program has its own venv and a normal app object:
( cd core/disease-programs/tuberous-sclerosis && ./venv/bin/python -m uvicorn api.main:app --port 8560 & )
```

Still unserved, and still registered `live`:

| Capability | Port | Why |
|---|---|---|
| `chemprop-admet` ⚠️ verified | 8572 | `chemprop` not installed |
| `esm2-search` ⚠️ verified | 8571 | `esm` not installed |
| `esmfold-model` ⚠️ verified | 8570 | model weights + CUDA |
| `molmim-nim` | 8001 | gated NIM (NGC) |
| `diffdock-nim` | 8002 | gated NIM (NGC) |

⚠️ **Three of these carry `maturity: verified`**, which the published honesty ledger defines as
"live **and additionally proven against real, recorded input**". The site's maturity matrix is
generated from the registry, so those badges are live on hcls-ai-factory.org while nothing serves
them. Either start the service or move the status — `scripts/validate_registry.py --probe` will
tell you which, and fails while the claim and the machine disagree.

**Adding any of these to `health-monitor.sh`** means adding a SERVICES row *and* respecting the
UI/UI+1 port convention, which `validate_registry.py` cross-checks. See `PORT_MAP.md`.
