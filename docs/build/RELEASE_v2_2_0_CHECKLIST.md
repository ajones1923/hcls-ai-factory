# Cutting v2.2.0 — what is ready, what must be re-measured, what is yours to decide

**Status:** repo-only working note. **Written:** 2026-11-03, against `fix/ten-out-of-ten` (PR #158).

The release machinery follows `skills/communication/release-and-site-publishing`: semver, an
**annotated** tag, a GitHub Release with human notes, CHANGELOG and `CITATION.cff` bumped in the
same change, green CI first. This file is the part that skill cannot know: whether *this* release
is honest yet.

## Why v2.2.0 and not v2.1.1 or v3.0.0

**MINOR.** The roster is unchanged — the same eight engines, eight agents and TSC program, 42
capabilities, 32 live and 10 planned. Nothing breaks: no API removed, no renamed field, no
changed call signature. But `retrieval_coverage`, `collections_contributing`, `api_endpoint`,
`gate_text()`, `gate_token_stream()`, `corpus_coverage` and two new operator scripts are all
additive capability, which is more than a patch.

One caveat worth stating in the notes rather than hiding: agents that previously answered
confidently from an empty corpus will now say they retrieved nothing. **No API contract changes,
but observable behaviour does** — for the better, and visibly.

## Ready now

- [x] **Lint** — `ruff check --select E9,F82,F811,F706,F707 core lib scripts` clean.
- [x] **Platform tests** — 650 pass, 9 skipped (`lib/hcls_common`).
- [x] **Subject suites** — 8,687 pass / 0 fail / 211 skip across 18 suites.
- [x] **Registry** — `validate_registry.py` OK: 42 capabilities, every engine/agent registered,
      and the new reverse guard confirms every registry-allocated port is supervised.
- [x] **Test-depth ratchet** — passes, no clinical subject regressed. (It had never run before
      this branch; see CHANGELOG.)
- [x] **CI green on a clean runner** — all six jobs, x86_64, Python 3.11 and 3.12. PR #158.
- [x] **CHANGELOG.md** written with the v2.2.0 entry under `[Unreleased]`.
- [x] **`CITATION.cff`** now carries `version` and `date-released` at all — it carried neither
      before, so the citation widget had no version to resolve.

## Must happen before the tag

1. **Merge PR #158 to `main`.** The tag must come off `main`, and CI's merge gate only truly
   runs there and on PRs into it.
2. **Rename the CHANGELOG section** from `[Unreleased] — proposed v2.2.0` to
   `## v2.2.0 — <date>`.
3. **Bump `CITATION.cff`** to `version: 2.2.0` and the release date, in that same commit. The
   file now says why, inline.
4. **Re-measure anything the notes claim.** See below — this is the part that matters.
5. `git tag -a v2.2.0 -m "..."` → push → publish the GitHub Release.

## Claims from the v2.1.0 notes that must NOT be copied forward

The v2.1.0 release table is the template, and three of its five rows would be false today:

| v2.1.0 claimed | Status for v2.2.0 |
|---|---|
| Capabilities: 42 — 32 live, 10 planned | **True**, verified from the registry. |
| Tests: 8,191 | **Now 8,687**, and for the first time verified by CI rather than one aarch64 box. |
| Services: **32/32 healthy** | **Do not restate.** Measured 2026-11-03: 28 listeners up, 11 down — nine agent UIs and the Discovery Portal among them. The fix is in this release (`health-monitor.sh`), so the number is only true *after* deploying it. Re-measure post-deploy or omit the row. |
| Demos: **17/17 running on real input** | **Unverified for this code.** Not re-run on this branch. |
| Clinical accuracy: **26/27 cases correct, every answer grounded in retrieved evidence** | **Do not restate either half.** `scripts/run_clinical_eval.py` has not been run against this code; and "every answer grounded" was the claim this release exists to make true — restating it as already-achieved would be the exact overclaim the honesty gate was built to catch. |

If you want those rows, they are three commands: deploy then re-probe for services,
`scripts/run_demo.py` across the roster for demos, `scripts/run_clinical_eval.py` for accuracy.
If you do not want to run them, **omit the rows** rather than carry last release's numbers.

## Decisions that are yours, and whether they block

| Item | Blocks v2.2.0? |
|---|---|
| **Merge #158** | **Yes.** Nothing else can proceed. |
| **`CLINICAL_REVIEW` §8 PhenoAge CRP** (2.40 years), **§12 MAGGIC** (creatinine ÷100 vs ÷88.4), and ten smaller entries | **No**, if the notes say thirteen entries are open and point at the file. That is honest for a research-use platform and better than a silent guess. §3 CYP2D6 — the one that read worst, because the repo contradicted itself — was **settled against the primary source** and is no longer on this list. |
| **Cardiology's twelve absent collections** | **No.** A data/deployment gap, declared in the Known Limitations and in `docs/build/EMPTY_COLLECTIONS.md`. The calculators are unaffected. Fixing it is one `setup_collections.py` run plus seeding, and would be a strong v2.2.1. |
| **The 5 patient-level collections** (synthetic cohort vs nothing vs IRB) | **No.** "Nothing, and the agent says so" is a defensible shipped state — and as of this release the agents actually do say so. |
| **Deploying to the box** | **No** — the tag describes the repo, not the running instance. But the live tree still runs 2026-09-17 code, so until it is deployed the released code is not the demonstrated code. |

## The honest one-line answer

**The code is release-ready; the release *notes* are the remaining work**, because the v2.1.0
table's service, demo and clinical-accuracy rows cannot be carried forward without re-measuring,
and re-measuring the service row means deploying first. Merge, drop or re-measure those three
rows, and v2.2.0 can be tagged the same day.
