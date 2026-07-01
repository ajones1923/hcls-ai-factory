---
name: 04-containers-and-orchestration-runtime
description: >-
  Best-practice standards for Pillar 04 (Containers & the Orchestration Runtime) of the HCLS AI Factory.
  Use when designing, building, operating, or reviewing how services are packaged and run — Docker
  Compose services and the Nextflow pipeline. Concrete triggers: adding a compose service, writing a
  Dockerfile, pinning an image, setting healthchecks/limits, or wiring a stage into the orchestrator.
---

# Pillar 04 — Containers & the Orchestration Runtime

How the factory is packaged and driven: a canonical Docker Compose that brings up the whole platform on
one box with one command, per-component composes for focused work, and a Nextflow DSL2 orchestrator that
runs the cross-stage DNA→drug pipeline and the trigger fabric.

## In the HCLS AI Factory
- **`docker-compose.dgx-spark.yml` is the canonical, exemplary compose.** It embodies the standards:
  required creds via `${VAR:?}` (`${HCLS_MINIO_USER:?}`, `${HCLS_MINIO_PASSWORD:?}`, `${GRAFANA_PASSWORD:?}`),
  pinned image tags (`milvusdb/milvus:v2.4.0`, `etcd:v3.5.11`, `minio:RELEASE.2024-01-01T...`,
  `prom/prometheus:v2.49.1`, `grafana:10.3.1`), a healthcheck on every service, `deploy.resources.limits.memory`
  per service, `restart: unless-stopped`, and no `privileged`. Shared config is DRY'd via YAML anchors
  (`x-common-env`, `x-agent-healthcheck`).
- **Per-component composes** exist for focused bring-up (`core/engines/*/docker-compose.yml`,
  `core/agents/*/docker-compose.yml`, and a `docker-compose.runpod.yml` for the burst path).
- **Nextflow DSL2 orchestrator (`hcls-orchestrator/`)** runs the cross-stage pipeline; the trigger fabric
  in `lib/hcls_common/` (`bidirectional_triggers.py`, `event_bus.py`, `circuit_breaker.py`) fans stage
  completions out to downstream engines/agents.
- **Deployment reality (documented in the compose header):** the *containerized* compose is the
  declarative/portable target; the *live* system runs the intelligence agents as native host
  processes (Python/Streamlit) supervised by `health-monitor.sh` via cron. Treat the compose as truth
  for ports, limits, and dependencies even where the host-process form is what's running.

## Best-practice standards
- **One-command bring-up:** `docker compose -f docker-compose.dgx-spark.yml up -d` must stand up the
  platform; new services join that file so nothing is left out of the golden path.
- **Pin every image to an exact tag — never `:latest`.** Reproducible builds on ARM64 depend on it.
- **Healthcheck every service** (reuse the `x-agent-healthcheck` anchor); gate dependents on
  `condition: service_healthy`, as the agents do on `milvus`.
- **Set `deploy.resources.limits.memory` on every service** — the box's 128 GB is shared; an unbounded
  container can starve the rest. Give memory-heavy services (imaging, single-cell) more, deliberately.
- **Creds via required env substitution (`${VAR:?}`), never defaulted in the file** — the compose fails
  loudly if a secret is missing rather than shipping a weak default.
- **`restart: unless-stopped`, no `privileged`, least capability.** Persist state to named volumes.
- **DRY shared config with anchors** (`x-common-env`) instead of copy-pasting env blocks.
- **Route multi-stage runs through Nextflow**, not shell glue — the orchestrator owns stage ordering,
  retries, and the trigger fabric.

## Do / Don't
**Do:** pin tags; healthcheck + memory-limit every service; `${VAR:?}` for creds; `depends_on:
condition: service_healthy`; wire new services into `docker-compose.dgx-spark.yml`; use YAML anchors;
run cross-stage work under Nextflow.
**Don't:** use `:latest` or unpinned digests; ship a service with no healthcheck or no memory limit;
default a password inside the compose (`${VAR:-changeme}` for a secret); run `privileged`; add a stage
as ad-hoc shell instead of an orchestrator process; assume the compose form is what's live without
checking the host-process note.

## Wiring it in
```yaml
# docker-compose.dgx-spark.yml — the shape every new service copies
x-agent-healthcheck: &agent-healthcheck
  interval: 30s
  timeout: 10s
  retries: 3
  start_period: 30s

your-new-engine:
  build: { context: ./core/engines/your-engine, dockerfile: Dockerfile }
  image: your-engine:0.1.0            # pinned, never :latest
  restart: unless-stopped
  ports: ["8577:8577"]                # claim a port in the documented map
  environment:
    <<: *common-env
    SOME_SECRET: ${SOME_SECRET:?set SOME_SECRET in .env}
  depends_on: { milvus: { condition: service_healthy } }
  deploy: { resources: { limits: { memory: 4G } } }
  healthcheck:
    test: ["CMD","python3","-c","import urllib.request;urllib.request.urlopen('http://localhost:8577/health')"]
    <<: *agent-healthcheck
```
- Bring-up: `docker compose -f docker-compose.dgx-spark.yml up -d`. Cross-stage pipeline:
  `hcls-orchestrator/` (Nextflow DSL2). Host-process supervision: `health-monitor.sh` (cron).

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **Multi-arch images:** base images and tags must have ARM64 (`linux/arm64`) variants or the build
  fails on the GB10 — verify before pinning.
- **No memory limit = shared-pool starvation.** Unlike a multi-node cluster, one greedy container hurts
  every other service; limits are load-bearing here, not cosmetic.
- **Compose ≠ live.** The agents run as native host processes under cron in the current deployment; if
  you "fix" something only in the compose, the running system may not change — mirror both.
- **`condition: service_healthy` needs a real healthcheck on the dependency** — Milvus gates every agent;
  a missing/incorrect Milvus healthcheck stalls the whole bring-up.
- **Anonymous volumes lose state on recreate** — name every stateful volume.

## Related
- Pillars: 01-compute-dgx-spark-remote-gpus, 03-networking-and-ingress, 13-reliability-and-operations
- build-housekeeping-standards
