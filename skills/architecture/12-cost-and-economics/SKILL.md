---
name: 12-cost-and-economics
description: >-
  Best-practice standards for Pillar 12 (Cost & Economics) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing where compute and LLM spend go and how it is tracked. Concrete
  triggers: tagging a capability's cost_class/gpu, deciding local-vs-RunPod for a heavy model, choosing
  a Claude tier for a call, or wiring a token/cost budget or metric.
---

# Pillar 12 — Cost & Economics

The factory's core economic claim is one **~$4,699 DGX Spark** running patient-DNA→drug-candidate
end-to-end, versus recurring cloud spend. Holding that claim means every capability declares its cost
posture, remote GPU is a deliberate exception, and LLM calls are routed by stakes and metered.

## In the HCLS AI Factory
- **One box is the baseline.** GB10 Blackwell GPU, 128 GB unified LPDDR5x, 20 ARM Grace cores — a
  fixed one-time cost the whole platform runs on. The economic argument is fixed-cost local compute vs
  metered cloud.
- **Every capability is cost-tagged.** `lib/hcls_common/capability_registry.py` carries `cost_class`
  (`low`/`medium`/`high`) and a `gpu` boolean on each `Capability`; `capabilities.json` uses all three
  classes today (~5 low / ~7 medium / ~3 high) with ~12 GPU-flagged entries. These tags let the MCP
  surface and composer reason about cost, and `registry.find(gpu=...)` filters accordingly.
- **RunPod is burst-only.** Heavy or ARM-incompatible models (RFdiffusion, Evo 2, Chai-1, …) run on
  RunPod over a private Tailscale mesh and register as native capabilities with a remote `endpoint` and
  honest `cost_class="high"` / `gpu=true` — transparent to the rest of the factory, but never the default.
- **LLM spend is tiered and metered.** Claude is used in three stakes-routed tiers — **Opus** for
  high-stakes send-ready clinical text + verification, **Sonnet** for analysis, **Haiku** for fast/cheap
  work. `lib/hcls_common/llm_client.py` carries a thread-safe daily token budget (`_TokenBudget`, `0 =
  unlimited`) that refuses calls once exhausted, a prompt cache, and Prometheus counters
  (`hcls_llm_input_tokens_total`, `hcls_llm_output_tokens_total`, api-calls, latency).

## Best-practice standards
- **Tag `cost_class` + `gpu` on every capability.** An untagged capability is invisible to cost
  reasoning; be honest — a GPU folding service is `high`/`gpu:true`, a stdlib calculator is `low`.
- **Prefer local; try the GB10 first.** Blackwell + bfloat16 runs many "big" models (ESMFold-scale, and
  potentially Chai-1) locally. Attempt local before reaching for remote GPU.
- **Burst remotely on purpose, and only for cause** — genuinely ARM-incompatible or too large for
  128 GB unified memory. Register the burst target like any native capability; don't hide the cost.
- **Route LLM calls by stakes, not habit.** Haiku for extraction/classification/routing, Sonnet for
  analysis, Opus only for send-ready clinical output and adversarial verification. Don't pay Opus
  rates for a classification.
- **Meter and budget tokens.** Configure `_TokenBudget` with a daily limit in production; scrape the
  token counters into Grafana so per-tier spend is visible (ties to Pillar 10).
- **Cache and validate before spend.** Use the prompt cache; run the input-validation gate first so a
  malformed request never burns a GPU job or an Opus call.
- **Right-size batch runs.** On one box, GPU jobs serialize — batch and schedule heavy work rather than
  fanning out concurrent GPU capabilities that thrash unified memory.

## Do / Don't
**Do:** set `cost_class` + `gpu` truthfully in the registry; try the GB10 before RunPod; register burst
targets as native capabilities; route by stakes (Haiku→Sonnet→Opus); set a daily token budget and scrape
the counters; cache and validate before spending.
**Don't:** default heavy models to RunPod out of habit; leave capabilities untagged; use Opus for
low-stakes work; run token spend unmetered; launch concurrent GPU jobs that exceed 128 GB; hide remote
cost behind a "local-looking" capability.

## Wiring it in
Registry (declare cost posture):
```json
{ "id": "rfdiffusion-design", "type": "model", "serving": "cloud_nim",
  "gpu": true, "cost_class": "high", "status": "planned",
  "endpoint": "runpod-host:8080", "invoke_path": "/design" }
```
```python
from hcls_common.capability_registry import get_registry
reg = get_registry()
[c.id for c in reg.all() if c.cost_class == "high"]   # audit high-cost surface
reg.find(gpu=True)                                    # everything that needs a GPU

# meter LLM spend
from hcls_common.llm_client import _budget           # daily token budget
_budget.configure(daily_limit=5_000_000)              # 0 = unlimited
```
Route by stakes at the call site (Opus = high / Sonnet = analysis / Haiku = fast), and scrape
`hcls_llm_*_tokens_total` into the Grafana dashboards from Pillar 10.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **RunPod is the recurring-cost leak.** The whole $4,699 story survives only if burst is the exception;
  a capability that quietly defaults to RunPod turns a fixed-cost box into a metered bill.
- **Unified memory is the real GPU budget.** 128 GB is shared CPU+GPU; two `high`/`gpu` capabilities
  running at once can OOM where each fits alone — cost and reliability are the same constraint here.
- **ARM wheels can force a remote burst you didn't price in.** A model with no Grace/ARM build silently
  becomes a `high`-cost remote job; tag it honestly rather than pretending it runs local.
- **Unbudgeted tokens are unbounded spend.** `_TokenBudget` defaults to `0 = unlimited`; a runaway
  agent loop on Opus is the one line item that can dwarf the hardware — set a limit and alert on it.
- **Opus-by-default is the silent overspend.** The stakes router exists precisely so most calls land on
  Haiku/Sonnet; a service pinned to the high tier for everything erodes the economic case.

## Related
- Pillars: 09-ai-orchestration, 10-observability, 11-security-and-secrets
- build-housekeeping-standards
