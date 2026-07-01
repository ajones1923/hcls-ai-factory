---
name: 10-observability
description: >-
  Best-practice standards for Pillar 10 (Observability) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing how the platform is monitored, traced, and made reproducible.
  Concrete triggers: adding /health or /metrics to a service, adding a Prometheus scrape or Grafana
  panel, instrumenting a span, or emitting a run-lineage + reproducibility manifest for a pipeline.
---

# Pillar 10 — Observability

Observability is knowing what the factory is doing on one box in real time — health, saturation,
latency, GPU/memory pressure — and being able to reconstruct exactly what any past run did.
On a single-box, clinical-adjacent platform, that spans live metrics, distributed traces, and
per-run reproducibility manifests.

## In the HCLS AI Factory
- **Metrics stack** (`docker-compose.dgx-spark.yml` + `monitoring/`): Prometheus at host **:9099**
  (container 9090, config `monitoring/prometheus.yml`) scrapes each agent's `/metrics`, plus
  `node-exporter` (**:9100**) and Milvus (**:9091**); Grafana at **:3000** auto-loads
  `monitoring/grafana/dashboards/` (provisioned datasource + dashboards). GPU telemetry comes from a
  DCGM exporter (**:9400**) as a scrape target.
- **Every service exposes `/health` and `/metrics`.** The governance middleware
  (`hcls_common.api_gate.install_governance`) also stamps every response with `X-Request-ID`,
  `X-HCLS-Governed`, and `X-HCLS-Duration-ms` for request tracing.
- **Distributed tracing** — `lib/hcls_common/tracing.py` (OpenTelemetry): `init_tracing(service)`,
  the `@traced` decorator (sync + async), and W3C `inject_trace_context()` /
  `extract_trace_context()` for propagation across HTTP and Nextflow stages. It **degrades to a
  no-op** if OTel isn't installed, so downstream code never guards imports.
- **LLM telemetry** — `lib/hcls_common/llm_client.py` exports Prometheus counters
  (`hcls_llm_input_tokens_total`, `hcls_llm_output_tokens_total`, api-calls, latency) and a daily
  token budget (see Pillar 12).
- **Run lineage** — `lib/hcls_common/mlops.py` (SQLite, default `data/mlops.db`) records every run's
  params/metrics/status through the `submitted→started→running→complete/failed` state machine, with
  `search_runs`/`best_run` and a staged model registry.
- **Reproducibility manifests** — `lib/hcls_common/reproducibility.py` snapshots hardware
  (`detect_gpu` via `nvidia-smi`), software versions, model weights, knowledge-base versions, inputs,
  and timing per run — framed for 21 CFR Part 11 / scientific reproducibility.
- **Host watchdog** — `health-monitor.sh` at the repo root.

## Best-practice standards
- **Two endpoints, always.** Every engine/agent/service exposes `/health` (liveness) and `/metrics`
  (Prometheus). No `/metrics` ⇒ it is invisible on the Spark; wire a scrape job in
  `monitoring/prometheus.yml` in the same change.
- **Trace the hot paths.** Decorate cross-service and long-running functions with `@traced(...)`; call
  `init_tracing(service_name)` once at startup; propagate context on every outbound call
  (`inject_trace_context()` → headers / Nextflow metadata).
- **Every model run gets a manifest.** Capture a `ReproducibilityManifest` (hardware + software +
  weights + reference-data versions) and an MLOps run record — never leave a clinical-adjacent result
  un-reconstructable.
- **Dashboard GPU and unified memory first.** The GB10's 128 GB is shared CPU+GPU; panel VRAM/RAM
  saturation, GPU utilization/temperature (DCGM), and per-service latency together.
- **Alert on saturation, not just liveness.** A `/health` 200 while unified memory is at 95% is the
  failure that actually hurts on one box — alert on the pressure signal.
- **Advance run status through the state machine.** Use `set_status()`; a run stuck in `running` for
  hours is your earliest signal of a wedged GPU job.
- **Honest, structured logs.** Log the `X-Request-ID` so a trace, its metrics, and its logs join up.

## Do / Don't
**Do:** expose `/health` + `/metrics` on every service; register a Prometheus scrape job; `@traced` the
hot paths; emit a reproducibility manifest + MLOps run per pipeline; dashboard GPU/unified-memory
saturation; alert on pressure.
**Don't:** treat a 200 on `/health` as "healthy" while memory/GPU is saturated; hardcode an OTLP
endpoint or LAN IP (use `OTEL_EXPORTER_OTLP_ENDPOINT`); let tracing become a hard dependency (keep the
no-op fallback); skip the manifest because "it's just a demo run."

## Wiring it in
```python
# startup
from hcls_common.tracing import init_tracing, traced, inject_trace_context
init_tracing("my-engine")                                  # OTLP via env, else console/no-op

@traced(name="fold", attributes={"component": "esmfold"})
async def fold(seq: str): ...

# per run: lineage + reproducibility
from hcls_common.mlops import get_mlops_store
from hcls_common.reproducibility import ReproducibilityManifest
store = get_mlops_store()                                  # data/mlops.db
rid = store.start_run("fold-vcp", capability="my-engine-fold", params={"seq_len": 806})
m = ReproducibilityManifest(run_id=rid); m.detect_gpu(); m.detect_python_packages()
m.save(f"outputs/{rid}/manifest.json")
store.set_status(rid, "complete")
```
Add the scrape target to `monitoring/prometheus.yml` (`targets: ['my-engine:8000']`, `metrics_path:
/metrics`); Grafana picks up any dashboard JSON dropped in `monitoring/grafana/dashboards/`.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **Unified memory is the shared cliff.** CPU and GPU draw from the same 128 GB; a node-exporter RAM
  panel and a DCGM VRAM panel can both look fine individually while their sum OOMs — dashboard the sum.
- **OTel exporter wheels are ARM-sensitive.** If `opentelemetry-exporter-otlp-proto-grpc` isn't
  installed on Grace, `init_tracing` falls back to the console exporter — expected, not a failure;
  don't paper over it with a fake endpoint.
- **`mlops.db` and dashboards are local state.** They live in the 711 GB working tree, are gitignored,
  and are **not** backed by a durable named volume by default — they vanish on a fresh clone. Treat run
  history as machine-local unless you back it up.
- **Host `localhost:9100` vs container targets.** node-exporter is scraped at `localhost:9100` while
  agents are scraped by container name — get the network context right or the panel silently reads empty.
- **RunPod-burst spans cross the Tailscale mesh.** Propagate traceparent to remote GPU jobs too, or the
  most expensive spans in the factory become invisible.

## Related
- Pillars: 09-ai-orchestration, 12-cost-and-economics, 11-security-and-secrets
- build-housekeeping-standards
