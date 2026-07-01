# HCLS AI Factory — working instructions

This file is loaded into every Claude Code session in this repo. It carries the north star and
routes to the skills; the skills themselves (in `skills/`, auto-discovered from `.claude/skills/`)
carry the detail and load on demand.

## North star (always)

> **No one should die of a disease we could have understood in time.**

Every decision — what to build, how to frame it, what to promise — serves this mission. **When a
technical choice and the mission conflict, the mission wins.** Full founding vision and principles:
the `hcls-core-vision-mission` skill.

## Standing guardrails (apply to everything)

- **Honest by construction.** Label preclinical / roadmap / research-use plainly; a `live`
  capability is never mock-served; cite clinical claims to a primary source. Overstatement betrays
  the mission more than any missing feature.
- **Decision support, not diagnosis.** All clinical output supports a qualified clinician — never
  autonomous diagnosis or prescribing.
- **"Elastic burst," not "all on one box."** Heavy / ARM-incompatible models burst to RunPod over
  a private mesh — say so; don't imply everything runs locally.
- **Real, never mocked.** Verify a capability with a real input + recorded output before calling it
  done.
- **Neutral repo.** No proprietary vendor or alternate-edition branding anywhere in tracked files
  (the pre-commit guard enforces it). Data / weights / secrets stay local — only code + docs publish.

## Which skill to consult (invoke with `/name`, or just describe the task)

Framework: **why → what → who → how → standards.**

- **Why** — `hcls-core-vision-mission` (the north star).
- **What** — `demo-foundation-alignment`: anchor every capability to the D1–D7 demo portfolio + keep
  the coverage matrix green.
- **How** — `build-housekeeping-standards`: the build gate (structure, governed app, tests,
  registry, neutrality, docs).
- **Ship a capability** — `capability-delivery-playbook` · new disease program —
  `disease-program-authoring` · add a model/NIM/frontier model — `model-integration`.
- **Honesty & compliance** — `clinical-claim-honesty` · `regulatory-compliance-posture`.
- **Clinical correctness** — `clinical-genomics-standards` · `pharmacogenomics-cpic` · `oncology-mtb`.
- **Communication** — `demo-script-authoring` · `release-and-site-publishing` · audiences:
  `broad-general-persona` + the five `persona-*` drill-downs.
- **Technical layers** — `01`–`16` (one skill per Core Pillar, under `skills/architecture/`).

Full index: `skills/README.md`. After adding or editing a skill, run `./scripts/sync-skills.sh`.

## Before finishing any change

Run the gate (also the CI merge gate — no `|| true`):

```bash
ruff check --select E9,F82,F811,F706,F707 core lib scripts   # real-bug lint
( cd lib/hcls_common && pytest -q )                          # platform tests
python scripts/validate_registry.py                          # registry coverage
```

Commit small, coherent units; branch off `main` for non-trivial work; confirm before anything
outward-facing (push, release, external calls). Details: `build-housekeeping-standards`.

## Repo shape

`core/{engines,agents,disease-programs}` (the 8 engines + 8 agents + verticals) · `lib/hcls_common`
(the One Platform: capability registry, MCP tool-surface, workflow composer, MLOps, governance
gates) · `hcls-orchestrator` (Nextflow) · `skills/` (this framework). See `docs/STRUCTURE.md`.
