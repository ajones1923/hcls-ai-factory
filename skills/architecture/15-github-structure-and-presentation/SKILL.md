---
name: 15-github-structure-and-presentation
description: >-
  Best-practice standards for Pillar 15 (GitHub Structure, Visibility & Presentability) of the HCLS AI
  Factory. Use when designing, building, operating, or reviewing repository layout, front-door docs, CI,
  releases, and first-impression presentability. Concrete triggers: restructuring `core/`/`lib/`, editing the
  README/CONTRIBUTING/SECURITY/CITATION, adding a CI job, cutting a release tag, or checking a newcomer's
  first browse.
---

# Pillar 15 — GitHub Structure, Visibility & Presentability

The repository is the product's storefront: a clean, neutral, CI-gated monorepo that a newcomer can browse in
minutes and trust. This pillar is the housekeeping foundation *made visible* — good structure that a stranger
can see and a real CI that a stranger can believe.

## In the HCLS AI Factory
- **Clean layout:** `core/engines/`, `core/agents/`, `core/disease-programs/`, shared platform in
  `lib/hcls_common/`, orchestrator in `hcls-orchestrator/`. Directories kebab-case, Python modules
  snake_case; each engine/agent self-contained (`src/ api/ app/ config/ tests/ README.md`). Layout is
  documented in `docs/STRUCTURE.md`.
- **Front-door docs:** `README.md` (the *Eight Engines · Eight Intelligence Agents · One Platform* hero +
  port map + repo map), `CONTRIBUTING.md`, `LICENSE` (Apache-2.0), `SECURITY.md`, `CITATION.cff` (v2.0.0),
  `CODE_OF_CONDUCT.md`, and `.github/` issue + PR templates.
- **Real CI, no theater:** `.github/workflows/ci.yml` runs three jobs for real (never `|| true`) — real-bug
  lint (`ruff --select E9,F82,F811,F706,F707`), the platform library test matrix (`hcls_common`, Python
  3.11/3.12), and capability-registry validation (`scripts/validate_registry.py`). These are the merge gate.
- **Tagged releases:** semver, currently v2.0.0 (also in `CITATION.cff`).
- **Lean repo:** tracked content is small (data, weights, and venvs are gitignored) — the 500 GB+ working
  tree stays local; only code + docs publish. Largest tracked file stays well under 1 MB.
- **Local guard:** `scripts/install-hooks.sh` installs `scripts/pre-commit-hook.sh` — a neutrality/secret/
  large-file guard that blocks non-neutral branding, real key shapes, and oversized files at commit time.

## Best-practice standards
- Keep the `core/{engines,agents,disease-programs}` + `lib/` shape; add capabilities *within* it — don't
  invent top-level directories or a "ninth engine" (verticals compose the horizontals under
  `disease-programs/`).
- Every engine/agent has its own `README.md` (purpose, port, how to run); the root README carries the 8/8
  hero framing and stays accurate.
- CI stays real and green: no `|| true`, no skipped gate. If a job can't run reliably yet, say so in a
  comment (as the ci.yml note does for per-component jobs) rather than faking a pass.
- Every engine/agent directory maps to a registered capability — CI fails otherwise; this is what keeps the
  structure honest and discoverable.
- Cut tagged, semver releases and bump `CITATION.cff` in the same change.
- Keep the repo lean: `data/`, weights, `outputs/`, `results/`, venvs, `*.bak` gitignored; never `git add -f`
  artifacts.
- Neutral by construction: no proprietary/vendor or alternate-edition branding in tracked files; install the
  pre-commit guard on every clone.
- Docs match code: move/rename/restructure updates `docs/STRUCTURE.md`, the README repo map, and any `cd`
  paths in the **same** change — no dead links, no stale paths.

## Do / Don't
**Do:** run `./scripts/install-hooks.sh` on clone; add a README + registry entry with every new
component; keep CI real; browse the repo as a stranger before shipping structural changes; update docs in
lockstep with layout.
**Don't:** add `|| true` to a CI step; commit data/weights/venvs; introduce blocked non-neutral tokens (the
guard rejects them); leave a directory with no README or no registered capability; let `docs/STRUCTURE.md`
lag the tree.

## Wiring it in
```bash
# the exact local mirror of the merge gate (run before every commit)
ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && pip install -e ".[dev]" -q && pytest -q )
python scripts/validate_registry.py            # every engine/agent dir → a capability

./scripts/install-hooks.sh                     # one-time: neutrality/secret/large-file guard
```
Front-door files to keep current: `README.md`, `CONTRIBUTING.md`, `SECURITY.md`, `CITATION.cff`,
`CODE_OF_CONDUCT.md`, `docs/STRUCTURE.md`, `.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- The working tree is 500 GB+ of local data — the danger is accidentally publishing it. Keep the gitignore +
  large-file hook strict; the tracked repo must stay lean and cloneable.
- Two engines historically nest a deployable app under `agent/` (`clinical-imaging/agent/`,
  `precision-oncology/agent/`) — document these exceptions in `docs/STRUCTURE.md` rather than "fixing" them
  blindly; they carry real history.
- `single-cell` appears in both `engines/` and `agents/` **by design** (compute vs. reasoning) — disambiguate
  by registry ID, don't dedupe it thinking it's a mistake.
- This is the neutral repo; edition-specific material belongs in the sibling repo and never merges here — the
  pre-commit guard enforces the token blocklist.

## Related
- Pillars: 14-ease-of-deployment, 16-github-netlify-site, 11-security-and-secrets
- build-housekeeping-standards
