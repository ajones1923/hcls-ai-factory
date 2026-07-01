# hcls_common — the HCLS AI Factory platform layer

The shared library that ties the eight engines and eight intelligence agents together into
**one platform**: a typed capability registry, an assistant tool-surface (MCP), an AI workflow
composer, single-box MLOps, and the governance gates. Everything here is stdlib-first with heavy
dependencies (Milvus, Anthropic, embeddings, RDKit, OpenTelemetry) gated behind optional extras.

## Install

```bash
pip install -e "lib/hcls_common"                 # core (pydantic, numpy, requests)
pip install -e "lib/hcls_common[milvus,anthropic]"   # add a vector DB + LLM client
pip install -e "lib/hcls_common[all]"            # everything
pip install -e "lib/hcls_common[dev]"            # tests + lint + types
```

After install, `import hcls_common` works from anywhere — no `sys.path` hacks.

## What's inside

| Area | Module(s) | What it does |
|---|---|---|
| **Capability Registry** | `capability_registry.py` | One typed source of truth for every engine/agent/model/service, with JSON-Schema-style port contracts. The honesty rule is enforced at registration: a `live` capability can never be `mock`-served. |
| **Assistant Tool-Surface (MCP)** | `mcp_server.py` | Drive the whole factory from Claude / Cursor / any MCP client: list / describe / health / invoke / plan / compose. |
| **AI Workflow Composer** | `workflow_composer.py` | Natural-language goal → a validated, executable, governed pipeline (shape-based wiring, cycle detection, self-repair, structured root-cause). |
| **Governance gates** | `verify_gate.py`, `governance.py`, `license_gate.py` | Clinical-overclaim honesty gate, dual-model second-opinion verification, and a permissive-license/dependency-pin gate. |
| **Single-box MLOps** | `mlops.py` | Experiment tracking, model registry, run lineage, stage transitions (SQLite — no warehouse). |
| **Reproducibility** | `reproducibility.py` | 21 CFR Part 11-style run manifests (GPU + package capture). |
| **Shared infra** | `config.py`, embedding/LLM/Milvus clients, `event_monitor.py` | Configuration, vector-DB client, LLM wrappers, event bus. |

## Run the tests

```bash
cd lib/hcls_common && pip install -e ".[dev]" && pytest
```

The pure-logic test suite runs in ~1s with no GPU/Milvus/torch required.

## Extending the platform

New engines and agents should register their capabilities here and route their API through the
governance gates rather than re-implementing them. See the repository `CONTRIBUTING.md` for the
"Add an Engine / Register a Capability" guide.
