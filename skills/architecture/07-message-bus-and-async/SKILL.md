---
name: 07-message-bus-and-async
description: >-
  Best-practice standards for Pillar 07 (Message Bus & Async) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing event-driven coordination, cross-stage handoffs, async endpoints, or
  batch/scheduled work. Triggers: emitting or handling a pipeline event, wiring a bidirectional trigger,
  adding an async FastAPI route, a per-agent scheduler, or batch job status/search.
---

# Pillar 07 — Message Bus & Async

The coordination fabric that decouples pipeline stages and agents: a file-audited event bus with
bidirectional cross-stage triggers, async FastAPI endpoints, and per-agent schedulers. On one box this
is filesystem-backed and in-process by design — no external broker.

## In the HCLS AI Factory
- **Event bus** (`event-bus`, `serving: none`, `event_bus.py`). File-based JSON-lines persistence for
  an audit trail (no Postgres/Kafka), a `heapq` priority queue, a process-wide singleton
  (`get_event_bus()`), Prometheus counters for events emitted/processed, sync publish for immediate
  handlers + async publish via a background thread, and **event replay from the audit log for crash
  recovery**. Typed `EventType` (`vcf_ready`, `annotation_complete`, `drug_candidates_ready`,
  `cart_analysis_complete`, `pgx_result_ready`, `drug_design_request`, …) with `EventPriority` and
  `EventStatus` (`completed`, …).
- **Bidirectional triggers** (`bidirectional_triggers.py`). Cross-stage rules that auto-hand-off
  between stages and agents: VCF ready → annotation + Milvus embedding; pathogenic variant → target
  hypothesis → MolMIM generation; RAG ↔ CAR-T evidence exchange; top candidates → report generation;
  and dynamic variant-driven target ID (replacing the hardcoded gene lookup).
- **Circuit breaker** (`circuit_breaker.py`) guards handlers that call flaky/remote services (e.g.
  RunPod-hosted models over Tailscale).
- **Async endpoints & batch status.** Services expose async FastAPI routes; long jobs are tracked as
  MLOps runs whose status walks the state machine `submitted → started → running → complete | failed`
  (`RUN_STATES` in `mlops.py`), making batch status queryable/searchable. Per-agent schedulers drive
  recurring work.

## Best-practice standards
- **Idempotent handlers.** Every subscriber must tolerate the same event twice (dedupe on `event_id`
  / a natural key) — replay-on-recovery and at-least-once delivery mean redelivery is normal, not
  exceptional.
- **At-least-once + dedupe, not exactly-once.** Don't design for exactly-once; guarantee the effect is
  idempotent and record what you've already processed.
- **Decouple stages through events, not direct calls.** A stage emits a typed event and moves on;
  downstream stages subscribe. Producers never block on, or import, their consumers.
- **Everything is audited.** Emit through the bus so the JSON-lines audit trail captures it; don't
  fire side effects that leave no event record. The audit log is the recovery source of truth.
- **Type your events.** Use the `EventType` enum and a typed payload; no stringly-typed ad-hoc topics.
  Add new types to the enum, don't invent free-form strings.
- **Set priority deliberately.** Use `EventPriority` so clinically urgent handoffs (e.g. a pathogenic
  finding) outrank bulk work in the queue.
- **Async I/O, offload CPU/GPU.** FastAPI routes are `async` for I/O; heavy compute (folding, docking,
  embedding) is dispatched as a tracked batch job, not run inline in the request.
- **Handlers fail loud but isolated.** A failing subscriber sets `EventStatus` failed and trips the
  circuit breaker for remote calls — it must not crash the bus or block sibling handlers.
- **Track every long job as a run.** Submit → status transitions → complete/failed via the MLOps store
  so status/search works and nothing is a silent background thread.

## Do / Don't
**Do:** publish typed events through `get_event_bus()`; make handlers idempotent and dedupe on
`event_id`; register bidirectional triggers for cross-stage handoffs; wrap remote/RunPod calls in the
circuit breaker; track long jobs as MLOps runs; keep async routes non-blocking.
**Don't:** couple stages with direct synchronous calls; assume exactly-once or single delivery; invent
untyped topic strings; run folding/docking inline in a request handler; swallow handler errors without
setting status; stand up Kafka/RabbitMQ/Redis-streams — the box uses the file-audited bus.

## Wiring it in
- Publish and subscribe through the singleton:
  ```python
  from hcls_common.event_bus import get_event_bus, EventType, EventPriority
  bus = get_event_bus()
  bus.publish(EventType.VCF_READY, payload={"vcf": path, "sample": sid},
              priority=EventPriority.HIGH)

  @bus.subscribe(EventType.DRUG_CANDIDATES_READY)
  def on_candidates(evt):
      if already_processed(evt.event_id):        # idempotent / dedupe
          return
      generate_report(evt.payload)
  ```
- Cross-stage automation lives in `bidirectional_triggers.py` — add a rule there rather than hand-wiring
  producer→consumer calls; inspect flow with `event_monitor.py`.
- Batch/long jobs walk the state machine via the MLOps store (`RUN_STATES`), so `/batch` status and
  search stay truthful.
- The bus/triggers/circuit-breaker are `serving: none` platform capabilities already in the registry —
  reuse them; don't re-register a parallel coordinator.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **In-process singleton, not multi-node.** The bus is one process's coordinator; it does not fan out
  across machines. RunPod-hosted models integrate as capabilities called over Tailscale, guarded by
  the circuit breaker — not as remote bus peers.
- **The audit log is on local disk.** It's the recovery source of truth — put it on a backed volume and
  rotate it; losing it loses replay. Don't let it grow unbounded.
- **Background threads share the 20 Grace cores and 128 GB pool.** Too many eager async handlers doing
  CPU/GPU work will contend with inference; offload to tracked jobs and bound concurrency.
- **Trigger cycles.** Bidirectional rules can loop (A triggers B triggers A) — dedupe on event id and
  add cycle/termination guards, mirroring the composer's cycle detection.
- **Non-idempotent replay double-acts.** Crash recovery replays the log; a handler that isn't
  idempotent will re-generate reports/molecules — always guard on a processed-key.

## Related
- Pillars: 08-inference-serving, 05-structured-databases, 06-vector-databases-and-embeddings
- build-housekeeping-standards
