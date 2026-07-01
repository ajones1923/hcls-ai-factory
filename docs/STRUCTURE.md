# HCLS AI Factory — Structure & Migration

## Current layout (live, all services running)
Engines and the shared library live at top level as independent, cleanly-named components:
`genomics-pipeline/ rag-chat-pipeline/ drug-discovery-pipeline/ large-molecule-pipeline/
single-cell-pipeline/ small-molecule-pipeline/ hls-orchestrator/ lib/hcls_common/`,
the agents under `ai_agent_adds/`, and infra under `deploy/ monitoring/ docs/ scripts/`.

Phase 0 + Phase 1 cleanup is complete: secrets quarantined, 249 GB backups moved off-tree,
the nested site-clones relocated to siblings, working-file cruft archived to `~/hcls-archive/`,
and a top-level `README.md` added.

## Target layout (deliberate migration)
```
engines/   genomic-foundation/ precision-intelligence/ therapeutic-discovery/
           clinical-imaging/ precision-oncology/ cardiology/ tuberous-sclerosis/
agents/    cart/ precision-biomarker/ pharmacogenomics/ precision-autoimmune/
           neurology/ rare-disease-diagnostic/ single-cell/ clinical-trial/
libs/hcls_common/   deploy/   infra/   docs/   examples/
```

## Why this is a migration, not an in-place rename
Nine capability services bind to current paths and the shared library is imported via
relative `sys.path` from them. The safe cutover, run deliberately (not while live):

1. `pip install -e libs/hcls_common` (remove the `sys.path` hacks) — the hard prerequisite.
2. `git mv` each component into `engines/` / `agents/` on a branch.
3. Update each service's start path + the registry endpoints (ports are unchanged).
4. Restart services one at a time; verify `/health` + the MCP `check_factory_health`.
5. `mkdocs build` + Netlify preview before publishing.

Full rationale + the per-component map: `HCLS_STRUCTURE_OPTIMIZATION_PLAN.md` (planning workspace).
