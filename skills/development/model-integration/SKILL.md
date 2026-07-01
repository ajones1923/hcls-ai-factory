---
name: model-integration
description: >-
  How to add a model, NIM, or frontier model to the HCLS AI Factory cleanly — choose the serving mode,
  register a typed capability with an honest status, verify a canonical input→output before flipping to
  live, and front the service with the governance gates. Use when standing up a new model endpoint,
  wiring a NIM, bursting a heavy/ARM-incompatible model to remote GPU, or clearing a model's license.
---

# Model Integration — add a model as a governed, registered capability

A model only exists to the factory once it is a **registered capability** with a truthful `serving`
mode and `status`. Local-first on the GB10; heavy or ARM-incompatible models burst to remote GPU over
a private Tailscale mesh and register exactly like a local one — transparent to callers. Nothing is
advertised that isn't real.

## 1. Choose the serving mode
The registry `serving` enum is the contract for *where and how* a model runs:

- **`native`** — a host FastAPI/Streamlit process on the GB10 (the default for our GPU model
  services): ESMFold `:8570 /fold`, ESM-2 search `:8571 /search`, ChemProp ADMET `:8572 /admet`,
  single-cell compute `:8573 /analyze`, molecule generator `:8574 /generate`, ProteinMPNN `:8578
  /design`.
- **`container`** — a docker service (Parabricks genomics, monitoring).
- **`local_nim`** — an NVIDIA NIM container on this box (MolMIM `:8001`, DiffDock `:8002`).
- **`cloud_nim`** — a hosted accelerated endpoint (gated/partnership frontier models, e.g. Chai-2).
- **`none`** — a library / in-process capability, no endpoint (LoRA fine-tune, ACMG SF rules, GWAS).
- **`mock`** — simulated; **never** advertised as real and may not back a `live` capability.

**Local-first, then burst.** Try the GB10 first (cu130 torch, bf16/half-precision). If the model has
no aarch64 + cu130 wheel or is too heavy for 128 GB unified memory, host it on remote GPU (RunPod
A100/L40S/H100) over the private Tailscale mesh and **still register it as a native capability** with a
remote `endpoint` (a Tailscale name/env var, never a hardcoded IP) and honest `cost_class`/`gpu`. The
frontier models that burst this way — **Chai-1, Chai-2, RFdiffusion, Evo 2** — are registry entries
like any other; the composer must not special-case them. Say "elastic burst," never imply everything
runs locally.

## 2. License gate — clear it before you build
Prefer **Apache-2.0 / MIT** weights so the platform stays open and forkable — the license is a
first-class gate, not an afterthought. Chai-1 is **Apache-2.0 code + weights, commercial-cleared**, so
it integrates freely. Chai-2 is **gated / partnership access (not open weights)** — it registers as
`cloud_nim`, `status: planned`, and its integration is explicitly contingent on securing access. If a
model's license forbids commercial or redistribution use, flag it and prefer an open alternative
(ESM-2, ProteinMPNN, MHCflurry, Chai-1 are the open reference set).

## 3. Register a typed capability
Add an entry to `lib/hcls_common/capabilities.json` with a full, honest contract:
```json
{
  "id": "my-model", "type": "model", "name": "My Model",
  "domain": "proteins", "serving": "native", "gpu": true, "cost_class": "medium",
  "status": "planned",
  "endpoint": "localhost:8579", "invoke_path": "/predict",
  "inputs":  [{ "name": "sequence", "shape": "scalar", "required": true }],
  "outputs": [{ "name": "structure", "shape": "structure" }]
}
```
- **Typed `Port`s.** Every input/output declares a `ValueShape`:
  `scalar`/`list`/`list_of_objects`/`map`/`file`/`structure`. Add `required`/`enum`/`minimum`/
  `maximum`/`default` where useful — these drive both composer wiring and the input-validation gate.
- Set `endpoint` + `invoke_path`, `serving`, `gpu`, `cost_class`, and an honest `status`.
- **A not-yet-live model is `planned`** — register it as a stub and flip to `live` only when the
  service actually answers on the box (the **Chai-1 pattern**: registered, build-in-progress, `status:
  planned`, `endpoint`/`invoke_path` already declared, flipped once `chai_lab` serves `/fold`).
- Map the model's directory in `scripts/validate_registry.py` (`COVERAGE`) and run it until it prints
  `OK`.

## 4. Verify before flipping to `live`
Feed a **canonical input** and record the **expected output** in the entry's description before
flipping `planned → live` — the way every live model already carries its proof: ChemProp ADMET
`caffeine BBB=0.98`, ESM-2 `ubiquitin cosine=1.0`, MHCflurry `NLVPMVATV vs HLA-A*02:01 → 0.97/17nM`,
ProteinMPNN `5L33 → ~45% seq-recovery`. The honesty rule is enforced at registration: **a `live`
capability may never be `mock`-served** (`capability_registry.py` rejects `status == LIVE and serving
== MOCK`). Any fallback (a NIM that can't start, a dropped Tailscale tunnel) returns a **clearly
labeled** mock/degraded result behind a circuit breaker — never a silent fake, and never reported as
`live`.

## 5. Wire the service behind the governance gates
Front the model's routes with the shared governance surface so validation + honesty are inherited:
```python
from hcls_common.api_gate import create_governed_app, require_valid_input
app = create_governed_app("my-model", capability_id="my-model")   # in api/main.py

@app.post("/predict")
def predict(payload: dict):
    payload = require_valid_input("my-model", payload)            # 422 on bad input
    ...
```
- **Every model service exposes `/health`** — compose healthchecks and the MCP `health` verb depend
  on it; an unhealthy service must report unhealthy, not answer with stale/mock data.
- **One model, one service, one port**; keep the port map stable (`:857x` protein/chem/cell, `:800x`
  NIMs, `:19530` vector database) and wire it into `docker-compose.dgx-spark.yml`.
- **GPU is a shared, finite resource.** Long inferences run as tracked MLOps jobs
  (`hcls_common.get_mlops_store()`), not inline in a blocking request; stage/unload big models to
  avoid OOM on the 128 GB unified memory.
- Gate heavy stacks (torch, bionemo, monai) behind **optional extras**, never a hard requirement of a
  light service.

## 6. Drive it through the platform
Once registered + live, the model is automatically callable from the MCP tool-surface
(`python -m hcls_common.mcp_server` → `list / describe / health / invoke / plan`) and wireable by the
workflow composer via input/output `ValueShape` — no per-model glue. That is the payoff for a truthful
typed contract.

## Do / Don't
**Do:** choose the honest `serving` mode; clear the license (prefer Apache-2.0/MIT); try the GB10
before remote GPU; register a full typed contract with `planned` until it answers; verify a canonical
input→recorded output before `live`; front routes with `create_governed_app`; expose `/health`; run
heavy inference as tracked jobs; label every fallback as mock.

**Don't:** serve a `live` capability from a mock; advertise a `planned` model as running; hardcode a
RunPod IP (use Tailscale names/env); imply a burst model runs locally; skip `require_valid_input`;
run folding/docking inline in a request; integrate a non-commercial-licensed model as if it were open;
special-case a remote model in the composer.

## Integration checklist
- [ ] **Serving** — mode chosen (`native`/`container`/`local_nim`/`cloud_nim`/`none`); GB10 tried before burst.
- [ ] **License** — cleared; Apache-2.0/MIT preferred; gated models noted contingent.
- [ ] **Register** — `capabilities.json` typed entry + `invoke_path`/`gpu`/`cost_class`/`status`; `COVERAGE` mapped; `validate_registry.py` = OK.
- [ ] **Verify** — canonical input → recorded expected output; `planned → live` only after it answers.
- [ ] **Govern** — `create_governed_app` + `require_valid_input`; `/health` exposed; port wired into compose.
- [ ] **Honesty** — fallbacks labeled mock behind a circuit breaker; live ≠ mock.

## Related
- `08-inference-serving` — Pillar 08: serving modes, the honesty rule, stakes-routed Claude tiers.
- `01-compute-dgx-spark-remote-gpus` — GB10 local-first + elastic burst to remote GPU over Tailscale.
- `09-ai-orchestration` — the registry contract, MCP surface, and compose-by-shape wiring.
- `capability-delivery-playbook` — the end-to-end ship checklist this specializes.
