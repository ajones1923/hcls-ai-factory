# Repository & Site Review

**Date:** 2026-08-16 · Covers batch-01 asks **7** (GitHub layout and best practices) and **10**
(deep site content review). Every check was run against the live repo and the published site.

---

## Part 1 — GitHub repository

### Community health: complete

| File | |
|---|---|
| README, LICENSE (Apache-2.0), CONTRIBUTING, CODE_OF_CONDUCT, SECURITY, CITATION.cff | ✅ all present |
| `.github/` — issue templates (bug, feature), PR template, `dependabot.yml` | ✅ |

Public · issues enabled · Apache-2.0 · default branch `main`. Nothing missing.

### CI: well-built, and it proved itself

`.github/workflows/ci.yml` runs four jobs — lint (real-bug ruff rules), the platform library on a
**3.11/3.12 matrix**, registry validation, and a **strict docs build**.

That last job earned its keep during this work. The first push (`f40281b`) shipped a `.gitignore`
bug that excluded `docs/build/` while `mkdocs.yml` referenced it:

```
f40281b   CI  failure     <- caught it
c599dd4   CI  success     <- after the fix
```

**Process note worth keeping:** the failure was visible in CI within a minute. Ten minutes were
spent instead polling the live site for a page that was never going to appear. *Check CI before
checking the deploy.*

### Gap closed: the subject suites were not gated

CI ran the **382** platform-library tests and none of the **8,402** subject tests. A change could
break all 17 engines/agents and still merge green.

Added a `subjects` job invoking `scripts/run_all_tests.py`. It uses the harness rather than a bare
`pytest` because two traps make a naive invocation wrong — `src/collections.py` shadows the stdlib
in 11 subjects, and `zarr`/`fast-array-utils` register pytest plugins that abort collection for
three. Both are documented inline in the workflow.

> **First run failed and the job is currently ADVISORY** (`continue-on-error: true`). The suites are
> green on the DGX Spark (aarch64); the runner is x86 with a different dependency set. It reports the
> signal without blocking merges — an always-red required check trains people to ignore CI. Remove
> `continue-on-error` once it is green.
>
> The same push also **broke the platform-library job**, which had been green: the new security tests
> import `fastapi` and `httpx` via `TestClient`, and those were not in the package's `[dev]` extra —
> so they passed locally (dev machines have FastAPI) and failed in CI. Fixed by declaring them, and
> verified in a clean venv built exactly as CI builds it.

### Recommendations

| # | Item | Effort |
|---|---|---|
| R1 | **Set the repo homepage** to `https://hcls-ai-factory.org` — currently unset | trivial |
| R2 | **Add topics** — none set. Hurts discoverability for an Apache-2.0 project (`precision-medicine`, `genomics`, `bioinformatics`, `nvidia-dgx`, `healthcare-ai`) | trivial |
| R3 | Enable branch protection on `main` requiring the CI checks | low |
| R4 | Triage the 15 open issues | low |
| R5 | Repo is **537 MB**, dominated by video blobs. LFS was evaluated and **rejected** — it trades repo size for per-build bandwidth on a git-linked Netlify site. The discipline is to commit only videos that changed | ongoing |

---

## Part 2 — Site content review

Checks were run against the registry and the live pages, not by reading prose.

### Accuracy: clean

| Check | Result |
|---|---|
| Pages naming a port owned by a different capability | **0** |
| Capability status claims contradicting the registry | **0** |
| Stale roster counts on published pages | **0** |
| Gene → locus claims | 0 mismatches (verified in the code audit) |
| Pharmacogenomics pairs | 44, all genuine CPIC |

The one apparent status mismatch was a false positive: `docs/demos/PRD.md` contains the sentence
*"registered `live` while unreachable"* — accurate prose describing the historical finding, not a
current claim.

### Mission wording: correct

**0 published pages** use the internal *"die of a disease we could have understood in time"*
phrasing. **4 pages** use the public softened line: *"No one should wait years for a disease we
could have understood in a day."* That distinction was a deliberate decision and it has held.

### Decision-support framing: strong where it matters

Live pages, counted on the rendered HTML:

| Page | Mentions |
|---|---|
| `/factory/agents/pharmacogenomics-intelligence-agent/` | 12 |
| `/factory/engines/genomics-engine/` | 9 |
| `/programs/tsc/` | 8 |

### Fixed

- **`$4,000` → `$5,000` (DGX Spark, $4,699)** in the branding research guide. The hardware price
  rose to $4,699, so "deploy for under $4,000" was a factual error about the platform's own cost.
  54 other references already said $4,699.

### Not published — and correctly so

Stale agent counts ("6 agents", "7 intelligence agents") survive in
`docs/infographic_prompts/` and root-level `HCLS_AI_FACTORY_*.md`. Both are excluded from the build
via `exclude_docs`, so **0 of these reach the site**. They are working documents. Worth a cleanup
pass eventually; not a publication risk.

### Recommendations

| # | Item | Effort |
|---|---|---|
| S1 | Refresh `docs/infographic_prompts/` to the 8·8·1 roster, or archive it | low |
| S2 | Add the decision-support frame to genomic-foundation and precision-intelligence docs — 5 of 17 subjects still never state it | low |
| S3 | Label CYP3A4→midazolam/cyclosporine as substrate relationships, not CPIC guideline pairs | trivial |
| S4 | The site now carries 22 videos in Adam's clone; re-verify after any future narration change | ongoing |

---

## Combined state

| | Before batch-01 | Now |
|---|---|---|
| CI coverage | 382 tests | **8,784** (382 lib + 8,402 subjects) |
| API authentication | 1 of 12 | **12 of 12** |
| `X-HCLS-Governed` honesty | claimed on every request | only when a gate ran |
| Site accuracy defects | unknown | 0 ports, 0 statuses, 1 price (fixed) |
| Repo health files | complete | complete |
