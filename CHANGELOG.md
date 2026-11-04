# Changelog

All notable changes to the HCLS AI Factory, newest first. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Releases before v2.2.0 predate this file; their notes live on the
[GitHub Releases page](https://github.com/ajones1923/hcls-ai-factory/releases) and are summarised
here from the tags rather than reconstructed commit by commit.

This platform is **research-use decision support for a qualified clinician — never autonomous
diagnosis**. Entries that change a clinical value cite the primary source that justifies it.

## [Unreleased] — proposed v2.2.0

A correctness and honesty release. No new engines, agents or disease programs; **the eight
engines, eight agents and the TSC program are the same roster as v2.1.0**. What changed is how
much of what they say can be trusted, and two clinical values that were wrong.

Everything below is on `fix/ten-out-of-ten` (PR #158) and gated on a green CI run — the first
independent verification this work has had, since `ci.yml` triggers only on push-to-main and
PR-to-main.

### Fixed — clinical

- **EuroSCORE II scored poor ventricular function with the moderate coefficient.** Four
  coefficients disagreed with Nashef 2012 Table 6, two of them holding copy-pasted values from
  the renal block. For a 72-year-old with EF 21–30, NYHA III, urgent, two procedures, predicted
  operative mortality moved from **2.93% to 4.83%**. All 2,046 cardiology tests passed before the
  fix and after it, because monotonic ordering was preserved — which is why range assertions
  could not see a coefficient wrong by a factor of 2.6.
  (Nashef et al., *Eur J Cardiothorac Surg* 2012;41(4):734–745, Table 6.)
- **HEART score treated a missing troponin as a normal one**, scoring 0 and moving patients
  toward the low-risk discharge band. It now raises.
- **PhenoAge alkaline phosphatase coefficient** 0.0019 → **0.00188**, the published value
  (Levine 2018 and its correction).
- **tPA was not cleared for an unknown stroke onset time**; an absent `onset_hours` read as 0.0.
- **`/variants/interpret` returned Pathogenic for almost any input.**
- **TSC PP4 was granted for having any phenotype at all**, and two citations were wrong.
- **Therapy ranking borrowed evidence from other variants and other tumours.**
- **A total coronary occlusion graded CAD-RADS 4B**; high-risk plaque was called on one feature.
- **A cardiac workflow that failed reported "completed".**

### Fixed — honesty

- **The output-honesty gate could only see one answer key, and only if the prose was long
  enough.** It now walks the payload (bounded: depth 6, 2,000 nodes) across 20 answer keys, and
  covers all three surfaces: the FastAPI response middleware (JSON), `gate_text()` in-process,
  and `gate_token_stream()` for SSE.
- **`X-HCLS-Governed` was read from a contextvar set in a different middleware context**, so the
  automatic gate never announced itself. `X-HCLS-Honesty`, documented since 2026-09-15, was never
  set at all.
- **Five agents asked the model for citations when nothing had been retrieved** and reported 0.5
  confidence on answers with no evidence behind them.
- **Canned imaging findings were served with `is_mock=False`.**
- **The pipeline fabricated success:** removed the Nextflow processes that invented chemistry
  under real tool names, the orchestrator's fabricated stages 4 and 7, the two portal buttons
  that simulated a run, and a `_generate_bionemo()` that returned a hand-written list without
  calling BioNeMo.
- **A DiffDock confidence was ranked as a kcal/mol affinity.**

### Fixed — security and ops

- **The API auth gate failed open.**
- **An uploaded filename chose where the imaging service wrote.**
- **Session ids were used directly as filesystem paths.**
- **MCP `invoke_capability` posted agent calls at Streamlit**: the registry advertises the UI port
  and the API is UI+1. Capabilities now carry an explicit `api_endpoint`.
- **Nine agent UIs were advertised `status: live` and supervised by nothing** — added to
  `health-monitor.sh`, each verified on a scratch port first.
- **The Discovery Portal could never start**, and bring-up re-exposed the UIs.
- **The test-depth ratchet crashed on every invocation** with `KeyError: 'orchestrator'` — a
  merge gate that had never run, invisible because nothing triggered CI.

### Added — telling the truth about retrieval

- **`hcls_common.corpus_coverage`**: how much of the retrieval budget a workflow declared
  actually answered. A search of an empty collection succeeds and returns zero rows, which is
  indistinguishable at every layer from a search that found nothing relevant.
- **`retrieval_coverage` and `collections_contributing` on responses** from clinical-trial,
  cardiology, neurology, rare-disease-diagnostic, single-cell and cart, with the disclosure
  appended to the answer itself once coverage is degraded. `collections_searched` counts the
  question asked; it reported 14 whether the corpus held 328 rows or none.
- **`scripts/check_corpus_weights.py`**: the same arithmetic against the live instance.
  Measured 2026-11-03 — **56 of 70 workflows retrieve from under three-quarters of the evidence
  they declare**, worst being cardiology at 2%, whose twelve `cardio_*` collections turned out
  never to have been created.
- **`scripts/ingest_literature.py`**: PubMed into any collection, schema-projected.
- **`docs/build/EMPTY_COLLECTIONS.md`**: all 40 empty collections classified by what can fill
  them — 17 with tooling already here, 13 needing a fetcher, 4 needing per-instrument curation,
  and **5 that no public source can fill because they are patient-level**.
- **125+ clinical reference cases** — a described patient in, a published number out, checkable
  by hand against the citation — for CHA₂DS₂-VASc (Lip 2010), HAS-BLED (Pisters 2010), HEART
  (Six 2008 / Backus 2013), ASCVD PCE (Goff 2013), EuroSCORE II (Nashef 2012), ACR TI-RADS 2017,
  Lung-RADS v2022, CAD-RADS 2.0 (Cury 2022), ACMG/AMP 2015 Table 5 (Richards), PhenoAge
  (Levine 2018).
- **`scripts/validate_registry.py` reverse drift guard**: a registry-allocated port that nothing
  supervises now fails validation.
- **This CHANGELOG.**

### Changed

- **Test count: 8,191 → 8,687** across the 18 suites, plus 650 platform tests. Green on a clean
  x86_64 runner across Python 3.11 and 3.12.
- **PhenoAge now declares its CRP unit convention** in every result as `crp_unit_convention`.
  Levine 2018, its published correction and Liu 2018 all print the equation with **no units for
  any biomarker**; mg/L versus mg/dL is a constant **2.40 years** on every result. The value is
  unchanged and the assumption is now stated rather than implied — see `CLINICAL_REVIEW.md` §8.
- **Docs corrected** where they overstated: the PDB record, three overstated capability claims,
  and the TAND prevalence figures now marked unverified.

### Known limitations

Stated here because a release is a promise the repo can keep.

- **The corpus is thin: 4,764 vectors across 114 collections, 40 of them empty.** The grounding
  fixes in this release make the affected agents say so out loud, which will look like a
  regression in a demo and is not one.
- **Cardiology's twelve collections do not exist in Milvus.** Its risk calculators are pure
  functions and unaffected; its retrieved output has always been empty.
- **Nine entries in `CLINICAL_REVIEW.md` are open** and await a reviewer's ruling, including
  CYP2D6 activity score 1.0 (where the repo contradicts itself) and MAGGIC's creatinine
  conversion. They are recorded in the repo rather than silently resolved.
- **No clinical evaluation figure is claimed in this release.** v2.1.0's "26/27 cases correct"
  is not restated, because `scripts/run_clinical_eval.py` has not been re-run against this code.

## v2.1.0 — 2026-09-17

Open Project release: a signed, dated anchor over the work at the point of public release.
42 capabilities in a typed registry (32 live, 10 planned), 8,191 tests.
[Release notes](https://github.com/ajones1923/hcls-ai-factory/releases/tag/v2.1.0).

## v2.0.0 — 2026-07-01

Consolidated single-repo release.
[Release notes](https://github.com/ajones1923/hcls-ai-factory/releases/tag/v2.0.0).

## v1.3.0 — 2026-03-30

See `docs/HCLS_AI_FACTORY_v1.3.0_RELEASE_REPORT.md` and the
[release notes](https://github.com/ajones1923/hcls-ai-factory/releases/tag/v1.3.0).

## v1.2.0 — 2026-03-01

Intelligence Agents, Quickstart and publications.
[Release notes](https://github.com/ajones1923/hcls-ai-factory/releases/tag/v1.2.0).

## v1.1.0 and earlier

See the [Releases page](https://github.com/ajones1923/hcls-ai-factory/releases).
