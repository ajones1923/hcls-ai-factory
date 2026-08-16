# Build Guide — from a cold clone to a running HCLS AI Factory

Companion to [`PRD.md`](PRD.md) and [`GAP_ANALYSIS.md`](GAP_ANALYSIS.md).
Target: **NVIDIA DGX Spark** — GB10, 128 GB unified LPDDR5x, 20 Grace cores, aarch64.

Every command below was run on the box on 2026-08-15. Where a step is blocked, it says so and why.

---

## 0. Prerequisites

```bash
uname -m            # aarch64
nproc               # 20
free -g             # ~119 GB total
nvidia-smi --query-gpu=name,power.draw --format=csv   # NVIDIA GB10
docker --version && docker compose version
```

---

## 1. Python environment  *(P0-1)*

The repo has **no** committed virtualenv, and system Python is PEP 668 "externally managed", so
`pip install -e` into it fails. A repo-level venv is the supported interpreter.

```bash
cd /home/adam/projects/hcls-ai-factory
python3 -m venv --system-site-packages .venv     # reuse the heavy aarch64 wheels already present
.venv/bin/pip install -e lib/hcls_common
.venv/bin/python -c "import hcls_common; print(hcls_common.__file__)"
```

> **If this fails with `No module named 'hcls_common'` while a `.pth` exists in
> `~/.local/lib/python3.12/site-packages/`:** a stale editable install is shadowing it. Check
> `__editable___hcls_common_0_1_0_finder.py` — if its `MAPPING` dict is `{}` the install is corrupt.
> Delete both `__editable__*hcls_common*` files and redo the step above. This is exactly what
> blocked the merge gate before this guide existed.

Extra runtime dependencies not in `pyproject.toml` but required by two engines:

```bash
.venv/bin/pip install duckdb statsmodels     # genomic-foundation: variant store + GWAS
.venv/bin/pip install biopython peft         # structural-biology: sequence IO + LoRA
```

## 2. Verify the merge gate  *(P0-1)*

```bash
.venv/bin/python -m ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && ../../.venv/bin/python -m pytest -q )
.venv/bin/python scripts/validate_registry.py
```

Expected:

```
All checks passed!
372 passed
registry: 42 capabilities, 9 typed 'engine'
OK — manifest valid and every engine/agent directory is registered.
```

> `9 typed 'engine'` is correct, not a bug. TSC is `type: engine` + `tags: [disease-program]`, and
> `scripts/site_gen_matrix.py::_is_program()` excludes it from the eight. Do not "fix" it — retyping
> drops TSC out of every `types=("engine","agent")` query, including the port drift-guard.

## 3. Whole-platform test sweep

```bash
.venv/bin/python scripts/run_all_tests.py
```

Expected: **17 subjects · 8,402 passed · 0 failed · 0 errors**.

Two traps this script exists to handle — do not hand-roll a replacement without them:

- **`src/collections.py` shadows the standard library** in 11 subjects. Putting their `src/` on
  `PYTHONPATH` kills the interpreter *before collection*:
  `cannot import name 'namedtuple' from partially initialized module 'collections'`.
- **`structural-biology/vendor_rfdiffusion/`** is vendored third-party code needing gated GPU
  packages (`dgl`, `rfdiffusion`). Excluded — it is not ours.

> If **every** subject reports the same result, suspect the harness, not the platform. During the
> audit a missing `pytest-timeout` made all 17 look broken, and a bad `PYTHONPATH` made all 8 agents
> look empty. cart alone has 415 tests.

## 4. Environment file  *(P0-2)*

```bash
cp .env.example .env && ${EDITOR:-nano} .env
```

Four variables are required or compose will not even parse:
`ANTHROPIC_API_KEY`, `HCLS_MINIO_USER`, `HCLS_MINIO_PASSWORD`, `GRAFANA_PASSWORD`.
Do not use `minioadmin/minioadmin`. `.env` is gitignored and must stay so.

Validate:

```bash
docker compose -f docker-compose.dgx-spark.yml config --quiet && echo VALID
```

## 5. Port convention  *(resolved 2026-08-15)*

**The registry advertises the UI port; the API is UI + 1.** Adopted by Adam and now enforced.

Neurology and precision-biomarker previously each claimed **both** 8528 and 8529 and could not run
together. Applying the convention exposed two further adjacency collisions, so three services were
re-seated — imaging 8525→8523, neurology 8529→8535, structural-biology 8579→8581 — and 8534 was
returned to the NV-Segment-CT NIM.

Full allocation: [`PORT_MAP.md`](PORT_MAP.md). `validate_registry.py` now fails the build on either
a convention violation or a `health-monitor.sh` port the registry does not allocate; both failure
modes were negative-tested.

```bash
.venv/bin/python scripts/validate_registry.py   # must print OK before bring-up
```

## 6. Bring up the platform  *(P0-2)*

`health-monitor.sh` runs from cron every 5 minutes and will restart services underneath you.
**Disable it first:**

```bash
crontab -l > /tmp/crontab.bak
crontab -l | grep -v health-monitor.sh | crontab -
```

Then:

```bash
docker compose -f docker-compose.dgx-spark.yml build          # first run: slow, many aarch64 wheels
docker compose -f docker-compose.dgx-spark.yml up -d milvus-etcd milvus-minio milvus
docker compose -f docker-compose.dgx-spark.yml up -d
docker compose -f docker-compose.dgx-spark.yml ps
```

Restore supervision when you are done:

```bash
crontab /tmp/crontab.bak
```

## 7. Verify every `live` capability answers  *(P0-3)*

```bash
.venv/bin/python - <<'PY'
import json, urllib.request
caps = json.load(open("lib/hcls_common/capabilities.json"))["capabilities"]
for c in caps:
    if c.get("status") != "live" or not c.get("endpoint"): continue
    url = f"http://{c['endpoint']}/health"
    try:
        code = urllib.request.urlopen(url, timeout=3).status
    except Exception as e:
        code = type(e).__name__
    print(f"{c['id']:42s} {c['endpoint']:22s} {code}")
PY
```

A capability that is `live` in the registry and unreachable here is an honesty-gate failure — fix
the service or change the status. This is how `singlecell-compute` was caught.

## 8. Deployment model — pick one  *(open decision)*

Two half-maintained paths exist:

| | Compose | `health-monitor.sh` |
|---|---|---|
| What | 17 subjects + infra | ~20 native host processes, cron every 5 min |
| Interpreter | image-internal | `./venv/bin/python` per service |
| Status | covers **17/17** | re-seated onto the convention; now gated |
| Blocker | none once §5 is settled | venvs exist for only **4 of ~20** services (agents have none) |

The compose header calls itself "the declarative/portable target" and says the live system is
host-process based. In practice the host-process path can start only the 4 services that have a
`venv/` (precision-intelligence, therapeutic-discovery, tuberous-sclerosis,
genomic-foundation/web-portal); the 8 agents and cardiology/imaging/oncology have none.
`health-monitor.sh::rebuild_venv()` can create them on demand, so the path is recoverable.
**Recommendation: standardise on compose** and reduce `health-monitor.sh` to container health
checks.

## 9. GPU enablement  *(P0-4, gated)*

```bash
.venv/bin/python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
# 2.10.0+cpu False        <-- current state
```

A CPU-only build on a GB10. See `GATED_SOFTWARE_PRD.md` for the NGC-credentialed install of CUDA
PyTorch for aarch64/sm_121, Parabricks, BioNeMo NIMs and Chai-1. Until then every ML path — ESM-2,
ESMFold, LoRA fine-tuning, local NIMs — is CPU-bound or non-functional.

## 10. What still cannot be done on this box

| Item | Why | Where |
|---|---|---|
| CUDA PyTorch | gated NGC wheels for GB10 | gated PRD |
| Parabricks (94 refs) | licensed container | gated PRD |
| BioNeMo NIMs — MolMIM (131), DiffDock (121) | NGC registry auth | gated PRD |
| Chai-1 co-folding | gated weights | gated PRD |
| RFdiffusion / dgl | GPU-only, vendored | gated PRD |
| `nextflow` | not installed | `curl -s https://get.nextflow.io \| bash` |
