---
name: 08-inference-serving
description: >-
  Best-practice standards for Pillar 08 (Inference Serving) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing how models are served — local GB10 services, Claude via API, NVIDIA
  NIMs, and on-demand remote GPU. Triggers: standing up a model endpoint, declaring `serving`/`status` in
  the registry, stakes-routing a Claude tier, adding a `/health` route, or wiring a mock fallback.
---

# Pillar 08 — Inference Serving

How every model becomes a callable, governed capability: served locally on the GB10, via the Claude
API in stakes-routed tiers, as NVIDIA NIMs, or on-demand on remote GPU — each declared honestly in the
capability registry so nothing is advertised that isn't real.

## In the HCLS AI Factory
- **Local on the GB10 (native, GPU)** — real FastAPI services on Blackwell (cu130 torch, half/bf16):
  ESMFold structure prediction :8570 (`/fold`), ESM-2 sequence search :8571 (`/search`), ChemProp
  ADMET :8572 (`/admet`, 104 endpoints), single-cell compute :8573 (`/analyze`, scanpy), molecule
  generator :8574 (`/generate`, BRICS), developability + design :8576 (`/score /optimize`), MHCflurry
  pMHC-I :8577 (`/scan`), ProteinMPNN inverse folding :8578 (`/design`).
- **Claude via API — three stakes-routed tiers.** Opus for high-stakes / send-ready clinical text and
  verification, Sonnet for analysis, Haiku for fast/cheap. Routed in `llm_client.py`; send-ready text
  must pass `assert_publishable(...)`/the verify gate.
- **NVIDIA NIMs** — MolMIM generator :8001 (`serving: local_nim`), DiffDock docking :8002
  (`local_nim`), GenMol (`planned`). Local NIM container on the box, or a cloud NIM endpoint.
- **On-demand / remote** — Chai-1 (Apache-2.0 code+weights, complex co-folding, `planned`, build in
  progress): try the GB10 first, else host on RunPod (A100/L40S/H100) over the private Tailscale mesh
  and register as a native capability, transparent to callers. Chai-2 is access-gated (`cloud_nim`,
  roadmap).
- **The `serving` field** (registry enum): `native` (host FastAPI/Streamlit) · `container` (docker) ·
  `local_nim` (NIM on this box) · `cloud_nim` (hosted accelerated endpoint) · `mock` (simulated, never
  advertised as real) · `none` (library/in-process). `status` is `live | planned | mock`.

## The honesty rule (enforced, non-negotiable)
A `live` capability may **never** be `mock`-served — `capability_registry.py` rejects the combination
at registration (`status == LIVE and serving == MOCK` → error). Any mock fallback is always labeled as
mock in the response and the registry (see the MolMIM/DiffDock/imaging entries: "mock fallback
labeled"). Never advertise a simulated result as real; a graceful fallback is fine only when it says so.

## Best-practice standards
- **Declare `serving` + `status` truthfully in the registry** for every model/service, with a typed
  `endpoint`, `invoke_path`, and `inputs`/`outputs`. If it isn't running, it's `planned`; if only a
  mock exists, it's `mock` — and it may not claim `live`.
- **Stakes-route the model tiers.** High-stakes/send-ready → Opus, analysis → Sonnet, fast/bulk →
  Haiku. Don't burn Opus on bulk classification or serve send-ready clinical text from Haiku.
- **Graceful, labeled fallback.** When a NIM/remote model is unavailable, return a clearly labeled
  mock/degraded result (or fail) — never a silent fake. Wrap remote calls in the circuit breaker.
- **Every service exposes `/health`.** Compose healthchecks and the MCP `health` verb depend on it;
  an unhealthy service must report unhealthy, not answer with stale/mock data.
- **One model, one service, one port.** Keep the port map stable (:857x protein/chem/cell, :800x NIMs,
  :19530 vector db) and non-overlapping; wire it into `docker-compose.dgx-spark.yml`.
- **GPU is a shared, finite resource.** Batch, queue, and bound concurrency; long inferences run as
  tracked MLOps jobs, not inline in a blocking request. Load/unload big models deliberately.
- **Verify before `live`.** Flip `planned → live` only after a real input produces a real, recorded
  result (the registry entries carry their verification, e.g. caffeine BBB=0.98, ubiquitin cosine=1.0).
- **Prefer permissive licenses.** Apache-2.0/MIT model weights (ESM-2, Chai-1, MHCflurry, ProteinMPNN)
  keep the platform open; the license gate enforces this.

## Do / Don't
**Do:** register `serving`+`status`+`endpoint` accurately; stakes-route Claude tiers; label every
fallback as mock; expose `/health`; run heavy inference as tracked jobs; try the GB10 before RunPod;
flip to `live` only after real verification.
**Don't:** serve a `live` capability from a mock (the registry blocks it); return an unlabeled
simulated result; send high-stakes clinical text from Haiku; run folding/docking inline in a request;
hardcode a RunPod IP (use Tailscale names/env); claim a model is running when the port is dead.

## Wiring it in
- Register the service (reference `esmfold-model` / `molmim-nim`):
  ```json
  { "id": "my-model", "type": "model", "endpoint": "localhost:8579",
    "invoke_path": "/predict", "serving": "native", "status": "live", "gpu": true }
  ```
  Map its directory in `scripts/validate_registry.py`; `python scripts/validate_registry.py` → `OK`.
- Health + invoke are uniform:
  ```bash
  curl -s localhost:8570/health
  curl -s -X POST localhost:8572/admet -d '{"smiles":"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"}'
  ```
- Drive any registered model from an assistant via the MCP surface
  (`python -m hcls_common.mcp_server` → `list / describe / health / invoke / plan`), or let the
  workflow composer wire it by input/output shape.
- Claude tiers route through `llm_client.py`; gate send-ready text with `assert_publishable(text,
  llm=...)` / the verify gate before it leaves the box.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **ARM/GB10 wheels.** Not everything has an aarch64 + cu130 build; some models resist the GB10 and
  belong on RunPod. Don't assume a PyPI wheel exists — verify on the box, and gate heavy stacks behind
  optional extras.
- **Unified 128 GB is shared** across GPU inference, Milvus, and DuckDB scans. Concurrently loading
  ESMFold + DiffDock + a large collection can OOM — stage loads, unload idle models, bound concurrency.
- **Remote (RunPod over Tailscale) adds latency + failure modes** — always behind the circuit breaker,
  with a labeled fallback and a health check; a dropped tunnel must not masquerade as a result.
- **NIMs need the right runtime.** local_nim containers assume a compatible GPU/driver; if the NIM
  can't start, the labeled mock path engages — never let that mock report as `live`.
- **Dead port, live claim.** The commonest honesty bug: registry says `live` but the port is down.
  Health-check in CI/compose and keep the registry synced with reality.

## Related
- Pillars: 07-message-bus-and-async, 06-vector-databases-and-embeddings, 05-structured-databases
- build-housekeeping-standards
