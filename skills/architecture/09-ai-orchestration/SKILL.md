---
name: 09-ai-orchestration
description: >-
  Best-practice standards for Pillar 09 (AI Orchestration) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing how engines, agents, models, and services are discovered, wired, and
  driven as one platform. Concrete triggers: adding a capability to the registry, exposing a tool over
  MCP, composing a multi-step pipeline from a goal, or wiring governance into an execution path.
---

# Pillar 09 — AI Orchestration

Orchestration is how the eight engines, eight agents, models, and platform services stop being
loose microservices and become **one drivable factory**: a typed registry every capability
declares itself to, an MCP tool-surface any assistant calls, and a composer that turns a goal
into a validated, governed pipeline.

## In the HCLS AI Factory
- **Capability registry** — `lib/hcls_common/capability_registry.py` + the JSON manifest
  `lib/hcls_common/capabilities.json` are the single source of truth. Every capability declares an
  `id`, `type` (`engine`/`agent`/`model`/`nim`/`stage`/`service`), typed `inputs`/`outputs` (each a
  `Port` with a `ValueShape`: `scalar`/`list`/`list_of_objects`/`map`/`file`/`structure`), an
  `endpoint` + `invoke_path`, `serving`, `gpu`, `cost_class`, and an honest `status`
  (`live`/`planned`/`mock`). `get_registry()` is the process-wide singleton.
- **MCP tool-surface** — `lib/hcls_common/mcp_server.py`. `FactoryTools` (transport-agnostic, fully
  testable with injected network adapters) plus `build_server()` (FastMCP, stdio or streamable-HTTP)
  expose six tools that drive the whole factory from Claude / Cursor / any MCP client:
  `list_capabilities`, `describe_capability`, `check_factory_health`, `invoke_capability`,
  `plan_pipeline`, `compose_workflow`.
- **Workflow composer** — `lib/hcls_common/workflow_composer.py`. `WorkflowComposer.compose(goal)`
  drafts a DAG (deterministic shape-graph backward-chaining, or an injectable LLM planner), then
  `repair()`s it (drops unknown/non-live nodes, strips shape-incompatible edges with a changelog),
  `validate()`s it (required inputs, port-shape matching, cycle detection via `_topo`), and `run()`s
  it in topological order through the tool-surface, with `root_cause()` giving a data-vs-system verdict
  on failure.
- **Governance wrapping** — execution flows through `execute_capability`/`invoke_capability`, which
  inherit the input-validation gate (`validate_inputs`), and services front their routes with
  `hcls_common.api_gate` (`create_governed_app` / `install_governance`). Orchestration is always governed.
- **Drift-guard** — `scripts/validate_registry.py` (in the CI merge gate) fails the build if any
  `core/engines/*` or `core/agents/*` directory has no registered capability.

## Best-practice standards
- **Register every capability.** New engine/agent/model/service ⇒ a `capabilities.json` entry with an
  accurate shape contract, before you call the work done. If it isn't in the registry, it is invisible
  to the MCP surface, the composer, and health checks.
- **Compose by shapes, never hardcoded endpoints.** Wire producers→consumers by matching `ValueShape`
  (`can_connect`, `producers_of`, `consumers_of`). A pipeline that names host:port literals defeats the
  registry and breaks the moment a service moves.
- **Keep `status` honest.** A `live` capability may never be `mock`-served — the registry rejects it
  (`_validate`). The MCP surface refuses to `invoke` anything that isn't `live` and labels mocks;
  `planned` stubs (e.g. `chai1-structure`) stay `planned` until the service actually answers.
- **Declare the full contract.** Give every `Port` a `shape`, `required`, and where useful `enum` /
  `minimum` / `maximum` / `default` — these drive both composer wiring and the input-validation gate.
- **Govern the execution seam, don't reinvent it.** Route all invocation through
  `FactoryTools.execute_capability` (the single shared executor) so the tool-surface and the composer
  never diverge and both inherit input validation.
- **Let the composer self-repair, then read the checklist.** Treat `compose()`'s `repair_log` +
  `checklist` as the definition of "runnable"; never `run()` a pipeline with `error`-severity issues.
- **Prefer the deterministic planner in tests/CI.** The composer works offline with no model spend;
  keep an injectable LLM planner optional, never a hard dependency.

## Do / Don't
**Do:** register a capability with a typed contract; drive services via MCP tools; wire by `ValueShape`;
keep `status`/`serving` truthful; run `python scripts/validate_registry.py` until it prints `OK`.
**Don't:** hardcode endpoints in orchestration logic; advertise a `planned`/`mock` capability as `live`;
bypass the executor with ad-hoc `httpx` calls; ship an engine/agent directory with no registered capability;
invent a shape not in the `ValueShape` enum.

## Wiring it in
Register (add to `lib/hcls_common/capabilities.json`):
```json
{
  "id": "my-engine-fold", "type": "model", "name": "My Fold Service",
  "domain": "proteins", "serving": "container", "gpu": true, "cost_class": "high",
  "status": "live", "endpoint": "my-engine:8570", "invoke_path": "/fold",
  "inputs":  [{"name": "sequence", "shape": "scalar", "required": true}],
  "outputs": [{"name": "structure", "shape": "structure"}]
}
```
Then map the directory in `scripts/validate_registry.py` (`COVERAGE`) and run it. Drive / compose:
```python
from hcls_common.mcp_server import FactoryTools
from hcls_common.workflow_composer import WorkflowComposer
ft = FactoryTools()
ft.check_factory_health()                        # ping every live endpoint
pipe, meta = WorkflowComposer(ft.registry, tools=ft).compose("fold VCP and score developability")
if meta["checklist"]["runnable"]:
    WorkflowComposer(ft.registry, tools=ft).run(pipe)
```
The orchestrator (`hcls-orchestrator/`, Nextflow DSL2) is the batch counterpart; W3C traceparent
propagation for inter-stage calls comes from `tracing.inject_trace_context()` (see Pillar 10).

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **One box, real ports.** `check_factory_health` opens real TCP sockets; on the Spark a capability
  can be "registered but down" simply because its container isn't started — that's a `system` verdict,
  not a bug in the composer.
- **RunPod-burst capabilities are still native registry entries.** A heavy/ARM-incompatible model
  bursting to a remote GPU registers with a remote `endpoint` and honest `cost_class`/`gpu`; the
  composer must not special-case it — transparency is the whole point.
- **`planned` is not a lie, `live`-on-a-mock is.** Flipping `planned → live` is a two-line change
  (`status` + `endpoint`), but only after the service actually answers on the box.
- **Cycles are silent until `_topo`.** Backward-chaining can propose a cycle if two capabilities share
  a shape both ways; always run `validate()`/`checklist()` before `run()`.

## Related
- Pillars: 10-observability, 11-security-and-secrets, 12-cost-and-economics
- build-housekeeping-standards
