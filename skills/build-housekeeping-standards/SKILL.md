---
name: build-housekeeping-standards
description: >-
  Best-practices and housekeeping standards for building the next generation of the HCLS AI
  Factory. Use whenever adding or modifying an engine, agent, capability, service, dependency,
  test, or doc in /home/adam/projects/hcls-ai-factory — so new work stays aligned with the
  consolidated, neutral, CI-gated foundation (registry drift-guard, governance middleware,
  dependency baseline, neutrality + secret guards, and docs that match the code).
---

# Build Housekeeping Standards — HCLS AI Factory

The factory is one neutral, Apache-2.0 monorepo: **Eight Engines · Eight Intelligence Agents ·
One Platform**, plus disease-program verticals. Engines and agents are horizontal capabilities
under `core/`; the shared platform is `lib/hcls_common/`. This skill is the standing checklist
for keeping every new build aligned with the foundation. When in doubt, the canonical sources are
`README.md`, `CONTRIBUTING.md`, and `docs/STRUCTURE.md`.

## 0. The non-negotiable gate (run before every commit)

A change is not done until all of these pass locally — they are also the CI merge gate
(`.github/workflows/ci.yml`), which runs for real (never `|| true`):

```bash
# real-bug lint (syntax, undefined names, redefinition) — not style noise
ruff check --select E9,F82,F811,F706,F707 core lib scripts
# platform library must stay green
( cd lib/hcls_common && pip install -e ".[dev]" -q && pytest -q )
# registry must cover every engine/agent directory
python scripts/validate_registry.py
```

Install the local guard once per clone so oversized files, secrets, and non-neutral material
are blocked at commit time: `./scripts/install-hooks.sh`. Bypass only with a real reason
(`git commit --no-verify`), and say why.

## 1. Golden path — add an Engine

Engines live in `core/engines/<engine-name>/` (kebab-case). Use `cardiology` as the template.

1. Scaffold the self-contained layout: `src/ api/ tests/ config/ README.md requirements.txt Dockerfile`.
2. Build a **pre-governed** app so the gates are inherited by construction:
   ```python
   from hcls_common.api_gate import create_governed_app
   app = create_governed_app("your-engine", capability_id="your-engine-id")
   ```
3. Register a capability (§3) — CI fails if the new directory has no registered capability.
4. Add a `README.md` (purpose, port, how to run) — all 8 engines have one; keep it that way.
5. Add tests; wire the service into `docker-compose.dgx-spark.yml` with a port + healthcheck.

## 2. Golden path — add an Agent

Agents live in `core/agents/<agent-name>/` (kebab-case). The 8 agents are near-identical — match them.

- `config/settings.py` (Pydantic), `src/agent.py`, `src/knowledge.py`, `src/rag_engine.py`.
- `api/main.py` + `api/routes/`, `app/<name>_ui.py`, `tests/`, `.env.example`, `README.md`.
- Wire governance (§4) and register a capability (§3).

Disease-specific verticals go under `core/disease-programs/<name>/` — they **compose** the
horizontal engines/agents, they are not a ninth engine.

## 3. Register a capability (the drift-guard)

Every engine/agent must appear in the registry so it is discoverable, composable, and governed.

1. Add an entry to `lib/hcls_common/capabilities.json`: unique `id`, a `type`
   (`engine`/`agent`/`model`/`nim`/`stage`/`service`), typed `inputs`/`outputs`, an `endpoint`,
   and `status` (`live`/`planned`). **A `live` capability may never be `mock`-served** — the
   registry rejects it. Never advertise a mock as real; label graceful fallbacks honestly.
2. Map the directory in `scripts/validate_registry.py` (`COVERAGE`) to the capability id(s).
3. Run `python scripts/validate_registry.py` → must print `OK`.

## 4. Governance — wire the gates

The input-validation and output-honesty gates live in `hcls_common.api_gate`. Use them so
governance is real, not optional (`cart` is the reference):

```python
from hcls_common.api_gate import install_governance, require_valid_input, honesty_flags
install_governance(app, service="cart", capability_id="cart-intelligence-agent")   # in api/main.py
payload = require_valid_input("cart-intelligence-agent", payload)   # 422 on bad input
flags = honesty_flags(answer_text)                                  # deterministic overclaim scan
```

For send-ready clinical text use `assert_publishable(text, llm=...)`. Reproducibility manifests
and MLOps run-tracking live in `hcls_common` too — use them, don't reinvent them.

## 5. Dependencies

- Target **Python 3.11+**.
- New work follows the shared baseline in root [`constraints.txt`]. Install with it applied:
  `pip install -r requirements.txt -c ../../constraints.txt`.
- Prefer compatible-release ranges (`>=X,<X+1`) in `requirements.txt`; keep exact pins in the
  constraints/lockfile, not scattered across components.
- Gate heavy/GPU stacks (torch, monai, scanpy, bionemo) behind **optional extras**, as
  `hcls_common` does — never a hard requirement of a light service.
- Prefer building on `hcls_common` over copy-pasting shared logic.

## 6. Neutrality, secrets, data & path hygiene

- **Neutral repo.** No proprietary vendor or alternate-edition branding anywhere in tracked
  files. The exact blocked terms are enforced by the neutrality guard in
  `scripts/pre-commit-hook.sh` — run it to see the current list before writing product prose.
  The `-vast` sibling repo is where edition-specific material belongs; it never merges here.
- **No secrets in git.** Keys/tokens go in env vars or local `.env` (gitignored); ship only
  `.env.example` templates with placeholders. The hook blocks real key shapes.
- **No machine-specific paths.** Never hardcode `/home/<user>/...` or a LAN IP in shipping code —
  derive from `Path(__file__)` / the repo root, or an env var like `$HCLS_ROOT`.
- **No default creds in composes.** Use `${VAR:-default}` env substitution (see the top-level
  `docker-compose.dgx-spark.yml`).
- **No data/weights/artifacts in git.** The 711 GB working tree stays local — only code + docs
  publish. Generated outputs (`outputs/`, `results/`, `test_outputs/`, logs, `*.bak`) are
  gitignored; never `git add -f` them. Largest tracked file should stay well under 1 MB.

## 7. Docs stay true

- Keep docs matching the code. When you move/rename/restructure, update the docs in the **same**
  change — no stale `cd` commands, no dead links, no path that no longer exists.
- Hold the **Eight Engines · Eight Intelligence Agents · One Platform** framing. "Three engines"
  is only correct for the DNA→drug *core pipeline* (stages 1–3), not the platform total.
- Every engine/agent has a `README.md`; `docs/` has an index.
- Honesty register: don't overclaim. Label anything not real; a "live" capability is never
  mock-served. Cite sources for clinical calculators (PMID/DOI/guideline).

## 8. Verify, then commit

- **Real, never mocked.** New capabilities must actually run; verify with a real input and record
  the result before calling it done.
- Test before commit; add/extend tests with every feature or fix; don't let the platform suite
  regress.
- Commit in small, coherent, well-messaged units. Branch off `main` for anything non-trivial.
  End commit messages with the project's `Co-Authored-By` trailer.
- Anything outward-facing (push, release, external calls) — confirm the manifest first, then act.

## 9. Quick reference

```bash
./scripts/install-hooks.sh                                   # one-time: local guard
ruff check --select E9,F82,F811,F706,F707 core lib scripts   # lint
( cd lib/hcls_common && pytest -q )                          # platform tests
python scripts/validate_registry.py                          # registry coverage
pip install -e "lib/hcls_common[dev]" -c constraints.txt     # shared lib + baseline
```

## 10. Alignment with the 16 Core Pillars

This skill operationalizes the housekeeping side of `skills/AI Factory Core Pillars.txt` —
especially **11 Security & Secrets**, **13 Reliability & Operations**, **14 Ease of Deployment**,
and **15 GitHub Structure / Visibility**. When a build touches a pillar, keep both this skill and
the pillar's intent satisfied.
