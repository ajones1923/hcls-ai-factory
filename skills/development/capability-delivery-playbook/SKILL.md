---
name: capability-delivery-playbook
description: >-
  The one end-to-end checklist for shipping any new capability into the HCLS AI Factory — engine,
  agent, model, NIM, stage, or platform service. Use it whenever you start, review, or finish a
  build: it ties the north-star mission, the demo portfolio, the housekeeping standards, the
  capability registry, and the CI gate into a single ordered flow so nothing ships half-wired,
  undemoed, unregistered, or overclaimed.
---

# Capability Delivery Playbook — ship anything, end to end

This is the orchestrating checklist that ties the whole framework together: **vision → demo home →
build → register → verify → honesty → docs → gate → commit.** Each step has its own deeper skill;
this playbook is the order you run them in and the invariants that must all hold before a capability
is "done." Skipping a step is how the platform drifts, overclaims, or grows a capability no demo can
show.

## The flow (run in order — each step gates the next)

### 1. Mission check — does it serve the north star?
Before writing code, confirm the capability advances the mission: *no one should die of a disease we
could have understood in time*. Favor integration over silos, access/affordability over complexity,
honesty over polish. If a technical choice and the mission conflict, the mission wins. If it doesn't
widen reach, deepen the genome→molecule story, or make a patient's afternoon more real, reconsider.
→ `hcls-core-vision-mission`.

### 2. Demo home — which D1–D7 does it serve?
Nothing ships that cannot be shown, honestly, inside the demo portfolio. Pick the capability's home
demo(s) from the D1–D7 matrix (or justify a new demo with a weight→compression→hope arc). Then
**update the coverage matrix** so the "every capability has a demo home" invariant holds — treat a
coverage gap like a failing test. If it maps to no demo, it's premature or a demo is missing; resolve
that first. → `demo-foundation-alignment`.

### 3. Build to the housekeeping standards
Build it neutral, self-contained, and governed by construction:
- **kebab-case layout** under `core/engines/<name>/`, `core/agents/<name>/`, or
  `core/disease-programs/<name>/` (`src/ api/ tests/ config/ README.md requirements.txt Dockerfile`).
- **Pre-governed app** so the input-validation + output-honesty gates are inherited:
  ```python
  from hcls_common.api_gate import create_governed_app
  app = create_governed_app("your-service", capability_id="your-capability-id")
  ```
- **Python 3.11+**, dependencies against the root `constraints.txt`, heavy/GPU stacks behind
  optional extras. Build on `hcls_common`, don't copy shared logic.
- **Tests** with every feature; expose `/health` on any service.
→ `build-housekeeping-standards`.

### 4. Register in the capability registry (the drift-guard)
If it isn't registered it is invisible to the MCP tool-surface, the workflow composer, and health
checks — and CI fails the build. Add an entry to `lib/hcls_common/capabilities.json`:
- unique `id`, a `type` (`engine`/`agent`/`model`/`nim`/`stage`/`service`), a `domain`;
- **typed** `inputs`/`outputs` — each a `Port` with a `ValueShape`
  (`scalar`/`list`/`list_of_objects`/`map`/`file`/`structure`), plus `required`/`enum`/`minimum`/
  `maximum`/`default` where useful (these drive both composer wiring and the input gate);
- `endpoint` + `invoke_path`, `serving`, `gpu`, `cost_class`, and an honest `status`.

Then map the directory in `scripts/validate_registry.py` (`COVERAGE`) to the capability id(s) — this
is the CI drift-guard that fails if any `core/engines/*` or `core/agents/*` directory has no
registered capability. → `09-ai-orchestration`, `model-integration`.

### 5. Verify — real, never mocked
Prove the capability actually runs before you call it live: feed a **canonical input** and record the
**expected output** in the entry's description, the way the registry already does
(`caffeine BBB=0.98`, `ubiquitin cosine=1.0`, `E31K instability 36.0→28.2`). This recorded
verification is what lets a reviewer trust the `live` label. The honesty rule is enforced at
registration: **a `live` capability may never be `mock`-served** (`capability_registry.py` rejects
`status == LIVE and serving == MOCK`). Flip `planned → live` only after the service answers on the
box with a real result. → `08-inference-serving`.

### 6. Honesty label — and the ledger
Assign an honest `status`: **live / planned / preclinical / research-use / roadmap / gated.** If it
is anything other than plainly live, add it to the honesty ledger in `demo-foundation-alignment` and
make sure any demo script says so out loud. Label graceful fallbacks as mock in the response. Clinical
output is **decision support for a qualified clinician — never autonomous diagnosis or prescribing**;
cite clinical claims (PMID/DOI/guideline).

### 7. Docs — keep them true to the code
Add/refresh the component `README.md` (purpose, port, how to run) and update any doc that names the
capability in the **same** change — no stale `cd`, dead links, or paths that no longer exist. Hold the
**Eight Engines · Eight Intelligence Agents · One Platform** framing ("three engines" is only the
DNA→drug core pipeline, not the platform total). → `build-housekeeping-standards §7`.

### 8. Pass the gate, then commit
Run the CI merge gate locally until all three are green, then commit in small, coherent units:
```bash
ruff check --select E9,F82,F811,F706,F707 core lib scripts   # real-bug lint
( cd lib/hcls_common && pytest -q )                          # platform tests
python scripts/validate_registry.py                          # registry + coverage → OK
```
Install the local guard once (`./scripts/install-hooks.sh`) so oversized files, secrets, and
non-neutral material are blocked at commit. End commit messages with the project's `Co-Authored-By`
trailer. Anything outward-facing (push, release) — confirm first, then act.

## Do / Don't
**Do:** start from the mission and a demo need; give every capability a demo home and a coverage-matrix
row; build governed-by-construction; register a full typed contract; verify with a canonical
input→recorded output before flipping to `live`; label anything not plainly live; keep docs in the same
change; make all three gate checks pass before committing.

**Don't:** build a capability with no demo home or no registry entry; advertise a `planned`/`mock`
capability as `live`; serve a `live` capability from a mock; ship without a recorded verification;
leave the coverage matrix or docs stale; present any output as autonomous diagnosis; commit past a red
gate without a stated reason.

## One-glance checklist
- [ ] **Mission** — advances the north star; integration/access/honesty upheld.
- [ ] **Demo** — home demo picked (D1–D7 or new); coverage matrix updated.
- [ ] **Build** — kebab-case, `create_governed_app`, tests, `/health`, Py 3.11 + `constraints.txt`, neutral.
- [ ] **Register** — `capabilities.json` entry (typed I/O, endpoint, serving, gpu, cost_class, status) + `COVERAGE` map.
- [ ] **Verify** — canonical input → recorded expected output; live ≠ mock.
- [ ] **Honesty** — status labeled; ledger updated if not plainly live; clinical claims cited.
- [ ] **Docs** — README + affected docs updated in the same change.
- [ ] **Gate** — lint + platform tests + `validate_registry.py` = OK; commit.

## Related
- `hcls-core-vision-mission` — step 1, why we build.
- `demo-foundation-alignment` — step 2, the D1–D7 coverage contract.
- `build-housekeeping-standards` — steps 3–8, how we build.
- `model-integration` / `disease-program-authoring` — the two specialized delivery paths.
