# Cutting v2.2.0 — what is ready, what must be re-measured, what is yours to decide

**Status:** repo-only working note. **Written:** 2026-11-03 against `fix/ten-out-of-ten` (PR #158).
**Revised:** 2026-11-04 against `integrate/v2.2.0-remediation` (PR #160), after an overnight
audit of every layer. PR #158 is merged; **#159 is subsumed by #160** (its NUDT15 commit is an
ancestor) and should be closed rather than merged.

> **The verdict changed overnight.** This file used to end with "the code is release-ready; the
> release *notes* are the remaining work." That is withdrawn. A ten-auditor read found roughly
> 300 defects, about twenty-five of them in a class where a clinician acting on the output could
> be harmed and nothing downstream would look wrong — none of it visible to a green test suite.
> Eight branches fixed them and are merged into #160. **Do not tag until #160 lands.**

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
- [x] **Platform tests** — 673 pass, 9 skipped (`lib/hcls_common`), measured on #160.
- [x] **Subject suites** — **10,225 pass / 0 fail / 0 err / 223 skip** across 18 suites, measured
      on the *merged* #160 rather than on any single branch. +1,298 over `main`; the
      orchestrator subject went 15 → 56.
- [x] **Registry** — `validate_registry.py` OK: 42 capabilities, every engine/agent registered,
      and the new reverse guard confirms every registry-allocated port is supervised.
- [x] **Test-depth ratchet** — passes, no clinical subject regressed. (It had never run before
      this branch; see CHANGELOG.)
- [x] **CI green on a clean runner** — all six jobs, x86_64, Python 3.11 and 3.12, on PR #160
      (`b8f44e5`). The runner counts **10,225 passed / 0 failed / 223 skipped — the same as a
      developer box**, which is the only evidence that a pass reflects the code rather than
      what happened to be installed. Getting there took fixing 37 failures that appeared only
      on the runner: nine gate tests that required a live Milvus, and twenty-eight that needed
      RDKit, which the remediation made a hard requirement by refusing to fake chemistry QC.
- [x] **CHANGELOG.md** written with the v2.2.0 entry under `[Unreleased]`.
- [x] **`CITATION.cff`** now carries `version` and `date-released` at all — it carried neither
      before, so the citation widget had no version to resolve.

## Must happen before the tag

1. **Merge PR #160 to `main`.** The tag must come off `main`, and CI's merge gate only truly
   runs there and on PRs into it.
2. ~~**Rename the CHANGELOG section**~~ — done. `## v2.2.0 — 2026-11-04`, with the overnight
   audit's fixes written up in four groups of their own.
3. ~~**Bump `CITATION.cff`**~~ — done: `version: 2.2.0`, `date-released: "2026-11-04"`.
   **If the tag slips past 2026-11-04, change both the CHANGELOG heading and this date**, or the
   citation will claim a release date that never happened.
4. **Re-measure anything the notes claim.** See below — this is the part that matters.
5. `git tag -a v2.2.0 -m "..."` → push → publish the GitHub Release.

## Claims from the v2.1.0 notes that must NOT be copied forward

The v2.1.0 release table is the template, and three of its five rows would be false today:

| v2.1.0 claimed | Status for v2.2.0 |
|---|---|
| Capabilities: 42 — 32 live, 10 planned | **True**, verified from the registry. |
| Tests: 8,191 | **Now 10,212**, and for the first time verified by CI rather than one aarch64 box. |
| Services: **32/32 healthy** | **Do not restate.** Measured 2026-11-03: 28 listeners up, 11 down — nine agent UIs and the Discovery Portal among them. The fix is in this release (`health-monitor.sh`), so the number is only true *after* deploying it. Re-measure post-deploy or omit the row. |
| Demos: **17/17 running on real input** | **Unverified for this code.** Not re-run on this branch. |
| Clinical accuracy: **26/27 cases correct, every answer grounded in retrieved evidence** | **Do not restate either half.** `scripts/run_clinical_eval.py` has not been run against this code; and "every answer grounded" was the claim this release exists to make true — restating it as already-achieved would be the exact overclaim the honesty gate was built to catch. |

If you want those rows, they are three commands: deploy then re-probe for services,
`scripts/run_demo.py` across the roster for demos, `scripts/run_clinical_eval.py` for accuracy.
If you do not want to run them, **omit the rows** rather than carry last release's numbers.

## Decisions that are yours, and whether they block

| Item | Blocks v2.2.0? |
|---|---|
| **Merge #160** | **Yes.** Nothing else can proceed. It carries #158's successor work plus all eight remediation branches. |
| **`CLINICAL_REVIEW` §8 PhenoAge CRP** (2.40 years), **§12 MAGGIC** (creatinine ÷100 vs ÷88.4), and the rest of the now-32 sections | **No**, if the notes say how many are open and point at the file. That is honest for a research-use platform and better than a silent guess. §3 CYP2D6 — the one that read worst, because the repo contradicted itself — was **settled against the primary source** and is no longer on this list. |
| **Cardiology's twelve absent collections** | **No.** A data/deployment gap, declared in the Known Limitations and in `docs/build/EMPTY_COLLECTIONS.md`. The calculators are unaffected. Fixing it is one `setup_collections.py` run plus seeding, and would be a strong v2.2.1. |
| **The 5 patient-level collections** (synthetic cohort vs nothing vs IRB) | **No.** "Nothing, and the agent says so" is a defensible shipped state — and as of this release the agents actually do say so. |
| **Deploying to the box** | **No** — the tag describes the repo, not the running instance. But the live tree still runs 2026-09-17 code, so until it is deployed the released code is not the demonstrated code. |

## Owner-only items — all three resolved 2026-11-04

Decided rather than deferred, with the reasoning in each commit:

| Item | Resolution |
|---|---|
| `therapeutic-discovery-engine` registered `status: live` | **Status kept.** The service is supervised on :8505, the 10-stage RDKit chemistry is real, and generation runs the real RDKit BRICS backend; MolMIM and DiffDock are separate capabilities already marked `planned`. Downgrading would have understated a capability that does real work. What was wrong was narrower: `tags` advertised `molmim` and `diffdock` beside `live`, and the caveat lived only in prose. Tags corrected, `maturity: preclinical` added. |
| Fabricated GPU telemetry | **Removed.** Five of six fields were a theoretical maximum scaled by utilization and jittered with `hash(str(time.time()))`. Now NVML readings only, `None` (rendering `--`) for anything unmeasurable, and the one surviving estimate declares itself in its own payload. Six regression tests, all failing against `main`. |
| HTTP routes bypassing the fixed deterministic workflows | **Split decision.** ACMG **wired**: `/v1/diagnostic/variants/interpret` claimed to classify and reported the model's prose label, so the corrected Table 5 classifier now governs it — the model proposes criteria, Table 5 decides the class. Trial eligibility **not wired, deliberately**: `/v1/trials/match` does not claim to screen, returns `criteria_evaluated: false`, and names the library class; adding a screening endpoint is a feature, not a correction. Recorded in the CHANGELOG's Known limitations so no release note can claim it. |

## Superseded notes

### Earlier framing of the same three items


Three things no agent can close, because each is a statement about intent rather than a defect:

| Item | Blocks v2.2.0? |
|---|---|
| **`therapeutic-discovery-engine` is registered `status: live`** | **Your call.** `status` is a hand-written string that overrides every honest signal the code emits. If it is not live, `mock` or `planned` is the honest value. |
| **Fabricated GPU telemetry** — `hash(time.time())` jitter rendered as instrument readings, `genomic-foundation/web-portal/app/server.py` | **Yes, in my view.** Invented numbers presented as measurements are the one thing this release is explicitly about not doing. |
| **HTTP routes that bypass the fixed deterministic workflows** — `app.state.workflow_engine` is an LLM passthrough | **Yes, for any note that claims the fix.** The corrected ACMG classifier and tri-state trial eligibility are real and tested, but **no HTTP caller reaches them**. Claiming them in the notes would be false; shipping them unclaimed is fine. |

## Documents that would close the most open work

1. **Klein 2009 Supplementary Appendix 1** — restores IWPC warfarin dosing (§21). Highest value by a distance: it is the difference between a disabled calculator and a working one.
2. **ACR Lung-RADS v2022 Table 1** — §6f, part-solid and non-solid nodules.
3. **Pocock 2013 Table 2** — §12, the MAGGIC SBP bands and the HF-duration direction.
4. **ClinGen Dosage Sensitivity / MYH7 VCEP** — §30, restores PVS1, which currently fires for no caller.
5. **Revised Ghent nosology** — §31.

## The honest one-line answer

**The code is not yet what the notes would say about it, and that gap is now measured rather
than guessed.** Eight branches of remediation are merged into #160 with 10,212 tests passing,
every clinical change citing a primary source and every unsourced change refused into
`CLINICAL_REVIEW.md`. Land #160, decide the three owner-only items above, and drop or
re-measure the v2.1.0 service, demo and clinical-accuracy rows — then v2.2.0 is honest.

**What changed overnight is not the amount of work left. It is that the work left is written
down.**
