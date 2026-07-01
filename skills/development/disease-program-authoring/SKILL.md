---
name: disease-program-authoring
description: >-
  How to add a new disease-program vertical under core/disease-programs/<name>/ — a disease-specific
  engine plus a cluster of disease-specific agents and a deterministic orchestrator that COMPOSE the
  eight horizontal engines and eight agents for one condition. Use when scoping or building a new
  program (NF1, NF2, Rett, Williams, the broader mTORopathies) with Tuberous Sclerosis as the working
  template. A disease program is a vertical, never a ninth engine.
---

# Disease-Program Authoring — build a vertical, not a ninth engine

Disease programs are the **vertical** solutions the horizontal engines + agents power: each takes one
condition and composes the platform's live capabilities into a coherent, clinician-facing program.
The canonical rule (`docs/STRUCTURE.md`): a disease program **composes** the eight engines and eight
agents — it does not add a ninth engine. Engine 7 is structural-biology, Engine 8 is single-cell; a
program's own engine is a *disease* engine (an orchestrator over the horizontals), scoped to the
condition. Tuberous Sclerosis Complex (`core/disease-programs/tuberous-sclerosis/`) is the first
program and the template every new one follows.

## Is this a strong program? (selection criteria)
Pick conditions that light up the whole factory *legitimately* — not by forcing capabilities that
don't fit. A strong candidate has all three:

- **Multi-system reach.** The disease touches many organs, so composing many engines/agents is
  clinically real, not contrived. TSC spans brain, kidney, heart, skin, and lungs — one patient
  legitimately recruits genomics, imaging, oncology-adjacent surveillance, trajectory modeling, and
  several agents at once.
- **A clean gene → mechanism → drug story.** A tractable molecular arc makes the program provable
  end-to-end. TSC is ideal: `TSC1/TSC2` loss-of-function → mTORC1 hyperactivation → mTOR inhibition
  (everolimus/sirolimus) as an on-label, mechanism-matched therapy.
- **A clinician beachhead.** A real specialist audience/site that would use it (for TSC: a pediatric
  center of excellence). The program should open a doorway into the whole factory for that audience.

Roadmap targets that fit the pattern: **NF1, NF2, Rett, Williams,** and the broader **mTORopathies**
(which share TSC's mTOR-pathway spine and can reuse much of its scaffolding).

## Layout — mirror the TSC template
Each program is a self-contained folder under `core/disease-programs/<name>/` (kebab-case) so it can
be deployed or replicated on its own:

```
core/disease-programs/<name>/
├── engine/ (or src/orchestrator/)  # deterministic orchestrator: dependency-ordered dispatch,
│                                   #   conservative failure, event-sourced state, review surfaces
├── agents/ (src/agents/<agent>/)   # disease-specific agents, uniform contract (runner+prompts+schema)
├── config/                         # which core engines/agents it composes; settings, model policy,
│                                   #   disease scoring configs (e.g. TSC's tand_scoring.yaml)
├── data/                           # disease reference data, synthetic cohort, ontology
├── api/  app/  tests/  docs/
├── README.md  requirements.txt  Dockerfile  .env.example
```

The TSC engine orchestrates **five** coordinated sub-agents (variant curator, trajectory modeler,
therapeutics strategist, phenome mapper, TAND surveillance) over an event-sourced core and a tiered,
offline-capable model router. Match that shape: a deterministic orchestrator + a small cluster of
disease-specific agents, each with a uniform contract and provenance.

## Compose, don't duplicate
The core discipline: the program's agents **call** the horizontal engines/agents through the registry
— they do not reimplement variant calling, folding, docking, annotation, or RAG. Declare the
capabilities the program composes in `config/` (its dependency set), and wire calls through the
registry / MCP tool-surface / workflow composer by `ValueShape`, never hardcoded endpoints. If you
find yourself rebuilding a horizontal capability inside a program, stop and compose the existing one
(add or extend it upstream if it's missing). → `09-ai-orchestration`.

## Register the program's capabilities
The program's own engine (and any distinct served surface) must appear in
`lib/hcls_common/capabilities.json` with `type: "engine"` and honest I/O — see the
`tuberous-sclerosis-engine` entry (`domain: "tsc"`, tags `disease-program`, orchestrates its 5
sub-agents, input `patient_id` → output `assessment`). Give it a demo home and a coverage-matrix row.
The disease-specific *sub-agents* live inside the program and are driven by its orchestrator; the
composed horizontal capabilities are already registered — don't re-register them.

## Give it a flagship demo (weight → compression → hope)
The unit of impact is the disease program itself: **one patient, the whole factory, one governed
afternoon.** Build the program toward a three-beat arc — **weight** (the burden / the stakes for this
patient), **compression** (the factory converges many modalities in minutes, not months), **hope** (a
mechanism-matched, structure-validated, trial-aware path forward). TSC is the flagship demo (D3) for
exactly this reason. → `demo-foundation-alignment`.

## Label honesty out loud
Disease programs carry strong claims, so label plainly:
- **Gene therapy / gene correction for the disease is typically preclinical** — the factory is the
  open design/analysis bench, not a cure today (mark it `preclinical`).
- Synthetic-cohort demonstrations are what run now; EHR/biobank/imaging-AI integrations are
  institutional phase work, not built — say so.
- Every agent output is **decision support behind a clinician review gate; not FDA-cleared; not
  autonomous diagnosis or prescribing.** Cite clinical claims.
- Reconcile any historical framing in the folder (e.g. TSC's docs predate the 6/8 model and say
  "Engine 7") toward "disease-program vertical" incrementally.

## Do / Don't
**Do:** compose the horizontal engines/agents through the registry; keep the program self-contained
and replicable; mirror the TSC orchestrator + sub-agent shape; declare composed capabilities in
`config/`; register the program engine and give it a weight→compression→hope flagship demo; label
preclinical/roadmap items honestly.

**Don't:** create a ninth engine; duplicate a horizontal capability inside a program; hardcode
endpoints instead of wiring by `ValueShape`; ship a program with no registry entry or demo home; imply
gene therapy is available today; present agent output as autonomous diagnosis; leave "Engine 7"-era
framing unreconciled.

## Authoring checklist
- [ ] **Selection** — multi-system reach + clean gene→mechanism→drug story + clinician beachhead.
- [ ] **Layout** — `core/disease-programs/<name>/` mirrors the TSC template (engine, agents, config, data).
- [ ] **Compose** — `config/` declares composed core capabilities; calls wired via the registry.
- [ ] **Register** — program engine in `capabilities.json` + `COVERAGE`; `validate_registry.py` = OK.
- [ ] **Demo** — flagship weight→compression→hope arc; coverage matrix updated.
- [ ] **Honesty** — preclinical/roadmap labeled; clinician review gate; not FDA-cleared; claims cited.
- [ ] **Gate** — build to housekeeping standards; lint + platform tests + registry validation pass.

## Related
- `capability-delivery-playbook` — the end-to-end ship checklist this specializes.
- `demo-foundation-alignment` — the flagship disease-program demo arc (D3 = TSC).
- `09-ai-orchestration` — compose-by-shape over the registry, don't duplicate.
- `build-housekeeping-standards` — the self-contained governed layout every program follows.
