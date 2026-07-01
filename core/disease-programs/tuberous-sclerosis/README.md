# TSC Intelligence Engine — Engine 7 (HCLS AI Factory)

Five coordinated agents plus a deterministic orchestrator that turn dispersed TSC
evidence (genomics, notes, imaging, trajectories) into reviewable clinician-facing
surfaces. Apache 2.0. Runs first on a single NVIDIA DGX Spark, with RunPod GPUs added
only for the cohort-build spike.

> **Scope.** The synthetic-data demonstration is what runs now. Epic/Clarity, biobank
> LIMS, and imaging-AI integration are institutional Phase-1 work, not built. Every
> agent output is decision support behind a clinician review gate. Not FDA-cleared.

See `docs/TSC_INTELLIGENCE_ENGINE_RESEARCH_PAPER.pdf` (master volume) and
`docs/TSC_INTELLIGENCE_ENGINE_PRD.pdf` (build spec). This package implements the PRD.

## Status

| Layer | State |
|---|---|
| Event-sourced core (13 event types, append-only log + projections) | **built, tested** |
| Deterministic orchestrator (dependency-ordered dispatch, conservative failure, surfaces) | **built, tested** |
| Uniform agent contract + provenance + tiered model router (offline-capable) | **built, tested** |
| 5 agents | **contracts wired** (Trajectory has a real classical fit); real reasoning = W2-W6 |
| Synthetic cohort pipeline | W1-W2 |
| RAG corpus, Streamlit surfaces, eval harness | W5-W8 |

## Run the walking skeleton (no external services needed)

```bash
cd ai_agent_adds/tuberous_clerosis_engine
python -m pytest tests/unit -q          # core + orchestrator tests
python scripts/dry_run_demo.py          # enroll featured patients A/B/C; assemble all 3 surfaces
```

With `TSC_ANTHROPIC_API_KEY` set, the LLM-tier agent steps call real Claude models;
without it, they return deterministic watermarked stubs so the engine runs end-to-end.

## Layout

`src/orchestrator` (events, state, graph, policies) · `src/agents/<agent>` (runner +
prompts + schema, uniform contract) · `src/cohort` (synthetic generation) · `src/rag`
· `src/utils` (provenance, model router) · `api` (FastAPI) · `app` (3 Streamlit
surfaces) · `config` (settings, model_policy.yaml, demo.yaml) · `tests`.
