---
name: 13-reliability-and-operations
description: >-
  Best-practice standards for Pillar 13 (Reliability & Operations) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing service uptime, health, backup/restore, or failure handling on the single
  DGX Spark box. Concrete triggers: adding a healthcheck or restart policy, wiring a circuit breaker around an
  external call, editing start/stop/health-monitor scripts, writing a RUNBOOK entry, or planning backups.
---

# Pillar 13 — Reliability & Operations

Keeping the factory up, honest, and recoverable when things fail — on one box with no cluster to fail over
to. Reliability here means: everything is health-checked, external calls are circuit-broken, degradation is
graceful *and labeled*, and every failure mode has a documented, tested remedy.

## In the HCLS AI Factory
- **Compose-level health + restart:** `docker-compose.dgx-spark.yml` defines a shared
  `x-agent-healthcheck` anchor (`interval: 30s`, `timeout: 10s`, `retries: 3`, `start_period: 30s`) applied
  to services, and `restart: unless-stopped` on all containers. Every service exposes `/health` (or
  `/healthz`).
- **Process lifecycle:** `start-factory.sh` (venv validation + auto-rebuild, port-wait, post-start health
  verification, `HEALTHY/TOTAL` summary), `stop-factory.sh` (kills by port), and `health-monitor.sh`
  (`status` / `fix` auto-restart / `watch` watchdog / `install` as a 5-min cron).
- **Circuit breakers:** `lib/hcls_common/circuit_breaker.py` — thread-safe `CircuitBreaker` with
  CLOSED/OPEN/HALF_OPEN states, usable as decorator, sync, or async context manager. Wrap Milvus, NIM, and
  LLM API calls.
- **Backups:** `core/engines/precision-intelligence/scripts/milvus_backup.py` for the vector database; a
  daily backup cron for state that matters (vector collections, DuckDB variant store, MLOps SQLite).
- **Runbook:** `docs/RUNBOOK.md` — real, reproducible failure modes (GB10/CUDA wall, aarch64 wheel gaps,
  422 annotation bug, Ts/Tv QC) each mapped to a fix. Supports the 21 CFR Part 11 reproducibility claim.
- **Graceful degradation:** `NIM_ALLOW_MOCK_FALLBACK` lets drug discovery fall back to mock data — always
  labeled, never silent. A `live` capability is never mock-served (the registry rejects it).

## Best-practice standards
- Healthcheck **every** service via the shared anchor; a service with no `/health` endpoint is not done.
- Set `restart: unless-stopped` on every container; never leave a service with no restart policy.
- Circuit-break every call that leaves the process (Milvus, NIMs, LLM/Claude API, RunPod) — reuse
  `CircuitBreaker`, don't hand-roll retry loops.
- Degrade gracefully and **label the degradation** in the response payload/UI; a fallback that looks real is
  a correctness bug, not resilience.
- Back up stateful volumes on a schedule **and test the restore** — an untested backup is a guess.
- Every failure you hit twice earns a `docs/RUNBOOK.md` row: symptom → cause → fix.
- Prefer idempotent start/stop: `start-factory.sh` no-ops on already-running ports; make new services do the
  same.
- Emit logs to a known path (`/tmp/hcls-factory/`, `logs/`) and cap growth (`health-monitor.sh` rotates at
  10 MB) — an unbounded log fills the shared 128 GB and takes the box down.

## Do / Don't
**Do:** wire `CircuitBreaker` around external I/O; add a healthcheck + `restart: unless-stopped` with every
new service; run `./health-monitor.sh fix` after any crash; document new failure modes in the RUNBOOK; test
that a backup actually restores.
**Don't:** ship a silent mock fallback; serve a `live` capability from mocks; hardcode `/home/<user>` or a
LAN IP in an ops script (derive from `$(dirname "${BASH_SOURCE[0]}")`); assume the box has failover — it has
none; let logs grow unbounded.

## Wiring it in
```python
from hcls_common.circuit_breaker import CircuitBreaker
cb = CircuitBreaker("milvus", failure_threshold=5, reset_timeout=30)

@cb
def search(query):        # rejects fast while OPEN, probes in HALF_OPEN
    return client.search(query)
```
```yaml
# docker-compose.dgx-spark.yml — inherit the shared anchor
myservice:
  restart: unless-stopped
  healthcheck:
    test: ["CMD", "curl", "-f", "http://localhost:PORT/health"]
    <<: *agent-healthcheck
```
```bash
./health-monitor.sh status          # one-shot health (add --json for machine use)
./health-monitor.sh fix             # check + auto-restart failed services
./health-monitor.sh install         # register the 5-minute watchdog cron
```

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **No UPS on the single box — state this honestly.** A power event is an unclean shutdown of *everything*;
  it is the strongest argument for tested backups and idempotent restart, and it is the known gap to name in
  any reliability write-up.
- One box = one blast radius: a runaway service starving the 128 GB unified LPDDR5x can take down unrelated
  services. Cap resources and watch memory, not just liveness.
- GPU services must start from the correct cu130 venv (see RUNBOOK) — a healthcheck can pass while the model
  silently ran on CPU; verify `torch.cuda.is_available()`, not just the port.
- `pkill -f` on a factory pattern can match the invoking shell (exit 144) — kill by **port**, as
  `stop-factory.sh` does.
- Milvus is TCP-only (`:19530`) — probe the port, not an HTTP path.

## Related
- Pillars: 04-containers-and-orchestration, 10-observability, 12-cost-and-economics
- build-housekeeping-standards
