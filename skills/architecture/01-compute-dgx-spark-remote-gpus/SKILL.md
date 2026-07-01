---
name: 01-compute-dgx-spark-remote-gpus
description: >-
  Best-practice standards for Pillar 01 (Compute, DGX Spark & Remote GPUs) of the HCLS AI Factory.
  Use when designing, building, operating, or reviewing where a model or job actually runs — the
  local GB10 box versus a remote burst GPU. Concrete triggers: adding a GPU capability, an ARM64+CUDA
  install failure, picking `serving`/`gpu` registry fields, or wiring a RunPod-over-Tailscale endpoint.
---

# Pillar 01 — Compute, DGX Spark & Remote GPUs

The factory's compute substrate: one NVIDIA DGX Spark (GB10 Blackwell GPU, 128 GB unified LPDDR5x
shared CPU+GPU, 20 ARM Grace cores, NVLink-C2C, ~$4,699) that runs the whole DNA→drug pipeline in
hours, with heavy or ARM-incompatible models bursting to remote x86 GPUs (RunPod) over a private mesh.

## In the HCLS AI Factory
- **Local, on the GB10:** the default. Modern-PyTorch + bfloat16 workloads run in-box on Blackwell —
  ESMFold (`facebook/esmfold_v1`, `core/engines/structural-biology/src/esmfold_service.py`, :8570,
  VCP-scale in ~1 s), scanpy single-cell compute (:8573), the RAG/agent stack, DuckDB, Milvus.
- **Unified memory** means the 128 GB is shared between CPU and GPU — there is no separate VRAM budget;
  one oversized model or a leaky container can starve everything else on the box.
- **Remote burst (RunPod over Tailscale):** x86 + standard CUDA for anything the ARM/Blackwell box
  can't serve — RFdiffusion (`vendor_rfdiffusion/`), Evo 2, Chai-1 (if it resists local), RAPIDS
  (`rapids-singlecell`), Geneformer, heavy folding. Documented in `docs/RUNBOOK.md` and the engine READMEs.
- **The registry is the contract:** every compute-bearing capability declares `gpu: true|false`, a
  `serving` mode (`native` / `container` / `local_nim` / `cloud_nim` / `mock` / `none`, from
  `lib/hcls_common/capability_registry.py`), and a `cost_class` — so where it runs is discoverable and
  a remote endpoint is transparent to callers.

## Best-practice standards
- **Local-first, burst-remote deliberately.** Try the GB10 first (modern torch + bf16 usually works);
  reach for RunPod only when an install/import genuinely fails or the job won't fit in unified memory.
- **Declare compute honestly in the registry.** Set `gpu`, `serving`, and `cost_class` truthfully; a
  `cloud_nim`/remote endpoint registers as a native capability so the composer and MCP call it unchanged.
- **Gate heavy/GPU stacks behind optional extras**, never a hard dependency of a light service (as
  `hcls_common` does) — a laptop clone must still `pip install` the platform lib.
- **Budget unified memory explicitly.** Treat 128 GB as one shared pool (CPU + GPU + page cache);
  size batch/precision to leave headroom for the co-resident agents, Milvus, and DuckDB.
- **Prefer bf16 on Blackwell**; don't assume fp16/older CUDA kernels — pin CUDA-13-class toolchains.
- **Keep remote endpoints stateless and idempotent** — the burst box can vanish; the caller retries.
- **Never hardcode a RunPod host/IP.** Endpoints come from env/registry; reach them over the Tailscale
  mesh (private DNS), not a public address.

## Do / Don't
**Do:** run ESMFold/scanpy/RAG locally; burst RFdiffusion/Evo 2/RAPIDS to RunPod x86; keep the same
client interface for local vs. remote (GPU acceleration is a drop-in swap behind the client); size to
the shared 128 GB; register `serving`/`gpu`/`cost_class` accurately.
**Don't:** assume x86 CUDA wheels install on ARM64; fight a CUDA 11.1 / `sm_86`-pinned model onto
Blackwell (`sm_121`); leave a mock serving a `live` capability; hardcode a burst-host address; run a
heavy folding/RAPIDS job in-box "just to see" and OOM the whole machine.

## Wiring it in
```jsonc
// lib/hcls_common/capabilities.json — a remote-burst capability, transparent to callers
{ "id": "rfdiffusion-backbone", "type": "model", "gpu": true,
  "serving": "cloud_nim",           // runs on RunPod x86; endpoint over Tailscale
  "cost_class": "high", "status": "planned",
  "endpoint": "http://<tailscale-host>:PORT/generate" }  // from env, never a literal LAN/public IP
```
- ARM/Blackwell failure playbook lives in `docs/RUNBOOK.md` (CUDA 11.1 → Blackwell `sm_121`, RAPIDS
  aarch64 gap, Geneformer build). Consult it before declaring a model "local".
- Local GPU service pattern: `uvicorn esmfold_service:_app_factory --factory --port 8570` (structural-biology).

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **CUDA 11.1 / torch 1.9 pins max out at Ampere `sm_86`; GB10 is Blackwell `sm_121` (needs CUDA 13)** —
  such models must go to RunPod x86 or a modern-torch equivalent.
- **`rapids-singlecell` installs but won't import** — no CUDA-13/Blackwell/aarch64 RAPIDS wheels; keep
  scanpy (CPU) local, run RAPIDS remote.
- **`x86_64` CUDA wheels silently absent on ARM64** — `pip` resolves to sdists that fail to build; gate.
- **Unified memory has no VRAM safety valve** — a GPU OOM can take CPU-side services down with it.
- **Grace has 20 cores, not a big x86 socket** — CPU-bound preprocessing is comparatively slow; profile.

## Related
- Pillars: 04-containers-and-orchestration-runtime, 08-inference-serving, 12-cost-and-economics
- build-housekeeping-standards
