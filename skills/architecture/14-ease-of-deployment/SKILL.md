---
name: 14-ease-of-deployment
description: >-
  Best-practice standards for Pillar 14 (Ease of Deployment) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing how quickly an open-source user can understand and stand up the factory.
  Concrete triggers: editing the compose file or Quickstart, adding a service that needs setup, writing a
  README or `.env.example`, shipping demo/fixture data, or touching the one-command bring-up path.
---

# Pillar 14 — Ease of Deployment

Because this is an Apache-2.0 open-source project, deployment *is* adoption: a newcomer must be able to
understand and run the factory as fast as possible. The bar is one command to a running stack and a real
first result with no hidden manual steps.

## In the HCLS AI Factory
- **One-command bring-up:** `docker compose -f docker-compose.dgx-spark.yml up -d` starts the whole stack;
  services report healthy within ~2–5 minutes (build cache dependent).
- **Front-door docs:** `QUICKSTART.md` (prerequisites → clone → key → up → logs → first demo curl) and the
  `README.md` Quickstart block. Full docs at https://hcls-ai-factory.org.
- **Minimal prerequisites:** Python 3.11+, Docker + Compose v2; NVIDIA GPU/CUDA 12.x optional (required only
  for the genomics engine). Everything else pulls via containers.
- **Config templates, not secrets:** a `.env.example` ships per component (15+ across engines/agents/
  programs); the only hard requirement is `ANTHROPIC_API_KEY`, with `GRAFANA_PASSWORD` optional.
- **Sensible defaults:** composes use `${VAR:-default}` substitution; `start-factory.sh` degrades to local
  Ollama when no `ANTHROPIC_API_KEY` is set, so the stack runs without paid keys.
- **Demo data ships in-repo:** small synthetic/public fixtures (`demo/`, `data/`, a tiny demo VCF) give a
  fast first result — no 1 TB download to see something work.
- **Drive-it-from-an-assistant one-liner:** `python -m hcls_common.mcp_server` exposes the whole factory as
  MCP tools (`list`/`describe`/`health`/`invoke`/`plan`/`compose_workflow`) to Claude, Cursor, any client.

## Best-practice standards
- Preserve the single-command up: any new service registers in `docker-compose.dgx-spark.yml` with a port,
  healthcheck, and defaults — never a separate manual start step.
- Every service ships a `.env.example` with placeholders and a `README.md` stating purpose, port, and how to
  run — match the existing engines/agents, all of which have one.
- Default to runnable: use `${VAR:-default}`, gate heavy/GPU stacks behind optional extras, and provide a
  no-GPU / no-paid-key path so the demo works on a laptop-class checkout.
- Optimize time-to-first-result: ship a tiny fixture and a copy-paste `curl` (the Quickstart cardiology
  `assess-risk` call) that returns something real in seconds.
- No hidden steps: if a human has to do it, it is in `QUICKSTART.md` or scripted — not tribal knowledge.
- Keep prerequisites minimal and explicit; if a step needs the GPU or a large download, say so and make it
  skippable for the demo path.
- Keep the two front doors (`README.md`, `QUICKSTART.md`) in sync with the actual ports and commands.

## Do / Don't
**Do:** add a `.env.example` + README + compose entry with every service; provide a working default and a
no-key/no-GPU fallback; ship a small fixture; give a one-line invocation that produces a real result.
**Don't:** require an undocumented manual step; commit real secrets (ship `.env.example` placeholders only);
make the demo depend on a multi-hundred-GB download; assume a GPU for the light path; let README/Quickstart
drift from the compose file.

## Wiring it in
```bash
# the entire deploy story a newcomer needs
git clone <repo> && cd hcls-ai-factory
export ANTHROPIC_API_KEY="sk-ant-..."          # or omit → local Ollama fallback
docker compose -f docker-compose.dgx-spark.yml up -d
docker compose -f docker-compose.dgx-spark.yml logs -f --tail=50

# first real result (seconds, from QUICKSTART.md)
curl -s http://localhost:8126/health | python3 -m json.tool

# drive the whole factory from an assistant
python -m hcls_common.mcp_server
```
New service checklist: compose entry (port + healthcheck + `${VAR:-default}`) → `.env.example` →
`README.md` (purpose/port/run) → registry entry → a Quickstart line if it is demo-facing.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- The full data footprint is ~1.1 TB (models + reference genomes + vectors) — never make it a prerequisite
  for *understanding* or a first demo; keep a fixture path that is orders of magnitude smaller.
- ARM aarch64 breaks pip installs that lack wheels (parasail, RAPIDS, older-CUDA models) — prefer
  containerized/pre-built images and pure-Python fallbacks so `up -d` doesn't force a source build on the
  user's box (see `docs/RUNBOOK.md`).
- GPU is optional but easy to accidentally make mandatory — keep the light services CPU-runnable so the
  Quickstart works before anyone provisions CUDA.
- Don't bake machine-specific paths or a LAN IP into compose/scripts; derive from the repo root so a fresh
  clone runs unchanged.

## Related
- Pillars: 15-github-structure-and-presentation, 04-containers-and-orchestration, 16-github-netlify-site
- build-housekeeping-standards
