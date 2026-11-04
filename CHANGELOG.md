# Changelog

All notable changes to the HCLS AI Factory, newest first. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

Releases before v2.2.0 predate this file; their notes live on the
[GitHub Releases page](https://github.com/ajones1923/hcls-ai-factory/releases) and are summarised
here from the tags rather than reconstructed commit by commit.

This platform is **research-use decision support for a qualified clinician — never autonomous
diagnosis**. Entries that change a clinical value cite the primary source that justifies it.

## v2.2.0 — 2026-11-04

A correctness and honesty release. No new engines, agents or disease programs; **the eight
engines, eight agents and the TSC program are the same roster as v2.1.0**. What changed is how
much of what they say can be trusted, and two clinical values that were wrong.

This release has two halves. The first landed as PR #158. The second came from an audit of every
layer of the platform on 2026-11-04 that found roughly 300 further defects — about twenty-five of
them in a class where a clinician acting on the output could hurt a patient and nothing
downstream would look wrong — and landed as PR #160, eight remediation branches merged together.

**The finding behind all the other findings: every one of those defects was invisible to a green
test suite.** 8,914 tests passed on a clean CI runner, across two Python versions, with the
ratchet and registry guards green, while all of it was true.

### Fixed — clinical

- **EuroSCORE II scored poor ventricular function with the moderate coefficient.** Four
  coefficients disagreed with Nashef 2012 Table 6, two of them holding copy-pasted values from
  the renal block. For a 72-year-old with EF 21–30, NYHA III, urgent, two procedures, predicted
  operative mortality moved from **2.93% to 4.83%**. All 2,046 cardiology tests passed before the
  fix and after it, because monotonic ordering was preserved — which is why range assertions
  could not see a coefficient wrong by a factor of 2.6.
  (Nashef et al., *Eur J Cardiothorac Surg* 2012;41(4):734–745, Table 6.)
- **CYP2D6 activity score 1.0 was classified Normal Metabolizer; the consensus calls it
  Intermediate** — and the opposite boundary was wrong too: AS 2.25 was classified Ultrarapid,
  where 2.25 is the inclusive top of Normal. An AS 1.0 patient (`*1/*4`, `*1/*5`, `*41/*41`)
  therefore received standard tricyclic dosing where CPIC asks for a 25% reduction, and an
  AS 2.25 patient (`*1x2/*41`) received the ultrarapid recommendation to change drug class. Both
  corrected in `pgx_pipeline.py` and in `dosing.py`, which carried its own copy of the same
  thresholds with the same two errors. The repo had contradicted itself: `knowledge.py` already
  stated NM as 1.25–2.25, and a test pinned the wrong value "per CPIC 2023".
  (Caudle et al., *Clin Transl Sci* 2020;13(1):116–124, Table 3, PMID 31647186 — "reclassify AS
  of 1 from NM to IM".)
- **NUDT15 `*2` was scored as a decreased-function allele where CPIC calls it no-function**, so
  `*2/*2` — a CPIC **poor metabolizer** — reported "NUDT15 Intermediate". That patient is given
  30–70% of a standard thiopurine dose where CPIC says avoid the drug or give 10%; the failure
  mode is life-threatening myelosuppression. `dosing.py` held the correct value all along — two
  copies of one clinical table, and the one the pipeline reports from was wrong. Also corrected:
  `*4`, `*5` and `*6` carried invented function assignments where CPIC says "unclear function",
  and both the translator and the dose calculator now return **Indeterminate** and name the
  allele rather than scoring it. (Relling et al., *Clin Pharmacol Ther* 2019;105(5):1095–1105,
  Table 1, PMID 30447069.)
- **An unrecognised allele was scored as fully functional** — `.get(allele, 1.0)`, commented
  "default to normal if unknown allele" — in the phenotype translator and in the thiopurine dose
  calculator for both TPMT and NUDT15. For genes whose unrecognised alleles are mostly
  loss-of-function, that is the least safe default available. Both now return Indeterminate.
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
- **`test_allele_tables_do_not_drift.py`**: allele activity values live in three copies in the
  pharmacogenomics agent with nothing reconciling them, and that has now produced two clinical
  defects. This fails the moment two copies disagree about an allele they both define. The
  duplication itself is recorded in `CLINICAL_REVIEW.md` §3c as still open — the fix is one
  table, imported, and belongs in its own change.
- **`scripts/validate_registry.py` reverse drift guard**: a registry-allocated port that nothing
  supervises now fails validation.
- **This CHANGELOG.**

### Fixed — clinical (2026-11-04 audit)

- **IWPC warfarin dosing contained a race-specific overdose, and is now DISABLED.** Every
  coefficient in the block is a real IWPC number attached to the wrong variable — the column is
  shifted. A Black patient picked up **+1.1816**, the enzyme-inducer coefficient, where Klein
  2009 applies **−0.2760**: a predicted dose inflated roughly **1.5–1.7×**, for that patient
  group only, in a narrow-therapeutic-index anticoagulant, in the direction that causes a bleed.
  It was returned under the name "IWPC Pharmacogenetic Warfarin Dosing Algorithm". Eighteen tests
  covered the function and every one asserted only a direction; three pinned the defect as the
  specification, including `test_african_american_increases_dose`. The function now raises, the
  UI explains why and tells the clinician to dose by protocol with INR monitoring, and
  `CLINICAL_REVIEW.md` §21 carries the published mapping so the restore is mechanical. Not
  corrected in place, because Klein 2009's Supplementary Appendix 1 could not be retrieved and
  this repo does not retype a clinical table from memory.
- **The recurring bug class: an absent or unrecognised input scored as the most benign value.**
  Found independently in cardiology, imaging, neurology, autoimmune, pharmacogenomics and TSC.
  `/gdmt/optimize` passed the wrong dict, so every lab, vital and contraindication was invisible
  — an MRA at K⁺ 6.1 and eGFR 18, an ARNI after angioedema. `/stroke/triage` held
  `tpa_eligible: true` at SBP 228 and ordered alteplase *ahead of* the CT that excludes
  hemorrhage. MoCA inflated the total by years of education, so a true 12/30 reported as 28/30
  "Normal cognition". TI-RADS scored an unrecognised feature label as 0, so a 13-point TR5
  nodule reported "TR1, benign, no FNA". CHA₂DS₂-VASc and HAS-BLED defaulted every unsupplied
  comorbidity to absent, so a true score of 5 reported as 1. The EuroSCORE II API dropped five
  terms and scored absent renal/LVEF/NYHA as best: 3.21% against 30.93%. A *negative* autoimmune
  lab scored its criterion as MET. Most of these needed no guideline at all — they needed the
  code to stop inventing.
- **ACMG/AMP classification was an additive point system, not Richards 2015 Table 5.** Table 5 is
  combinatorial over strength buckets, and the footnote that contradictory pathogenic and benign
  evidence yields VUS was absent, as was the BA1 handling that kept an *HFE* C282Y homozygote
  from being called Benign. Rebuilt against Table 5, and **now reached by an HTTP caller**: the
  model proposes criteria from retrieved evidence, Table 5 decides the class.
  (Richards et al., *Genet Med* 2015;17(5):405–424, Table 5.)
- **Trial eligibility is tri-state**, so a criterion no rule understands is UNKNOWN rather than
  met. An empty patient profile previously scored the same as a matched one; it now scores
  strictly below.
- **Pharmacogenomics allele tables rebuilt from CPIC's maintained machine-readable tables**
  rather than from the guideline papers, with a refusal path for alleles CPIC cannot score. The
  lesson is recorded at the top of `CLINICAL_REVIEW.md`: prefer the maintained table over the
  paper, cite it with a retrieval date, and check the whole table rather than a sample.

### Fixed — honesty and fabrication (2026-11-04 audit)

- **Three honesty-gate surfaces could not fail.** The verification gate stamped
  `status: "verified"` during an LLM outage, because every exception was swallowed into an empty
  verdict set and zero refuted claims read as success. `honesty_check` used `rx.search`, so only
  the first regex match was examined and a third-party sentence laundered a later
  self-referential block claim. Unadjudicated claims counted as supported.
- **Fabricated GPU telemetry rendered as instrument readings.** Five of the six fields in the
  genomics portal's detailed GPU metrics were a theoretical GB10 maximum scaled by utilization
  and then jittered with `hash(str(time.time()))`. Replaced with NVML readings only; a field this
  machine cannot measure is `None` and renders `--`, never a plausible-looking number. The one
  remaining estimate ships `throughput_is_modeled: true` and says so in its own payload, and the
  chart is no longer titled "Real-time Processing Throughput".
- **`therapeutic-discovery-engine` advertised two backends it does not run.** `tags` listed
  `molmim` and `diffdock` beside `status: live`; both NIMs are separately registered `planned`.
  Status stays `live` — the 10-stage RDKit chemistry and BRICS generation are real and the
  service is supervised — and `maturity: preclinical` now carries in a field what had been a
  parenthetical in prose.

### Fixed — security (2026-11-04 audit)

- **Two unauthenticated routes on services with no authorization**: one wrote to a path the
  deployment machinery later reads, one returned an engine's environment file. Both now behind
  the key gate.
- **Two independent key gates per service disagreed.** One read `HCLS_API_KEY`, the other a
  per-service `API_KEY`, with different open-path sets, so which endpoints were credentialed
  depended on which setting the operator happened to use — and `/governance` reported on only
  one of them. Now one resolved key, compared in constant time, with both gates leaving the same
  paths open.
- **An unbounded read of an entire upload into memory** in the autoimmune service.

### Fixed — process and CI (2026-11-04 audit)

- **The README test-count badge guard had never run once.** Its regex required a `+` the badge
  did not contain, and the no-match branch returned PASS. Four counts had drifted apart behind
  it. A badge that exists but does not parse is now a failure.
- **`run_all_tests.py` ignored pytest's exit code**, so zero collection, a conftest crash and a
  timeout each reported a clean subject and exited 0.
- **18 tests were failing in `precision-intelligence/portal/tests/`**, excluded by `testpaths`,
  so nothing ran them.
- **The suite passed on a developer box and failed on a clean runner.** Nine gate tests required
  a live Milvus only because `with TestClient(app)` runs app startup, and twenty-eight needed
  RDKit, which this release makes a hard requirement by refusing to certify unexamined molecules
  as having passed drug-likeness QC. CI and a developer box now report the same number, which is
  the only evidence that a pass reflects the code rather than what happened to be installed.
- **Five wrong ports in the architecture diagram embedded at `README.md:22`** — two of them the
  live API port of a different agent. All sixteen labels now resolve to a supervised entry in
  `health-monitor.sh`.
- **The README overstated the API key gate.** `DEFAULT_OPEN_PATHS` exempts eight paths including
  `/docs`, `/openapi.json`, `/metrics` and `/governance`; the open set is now named explicitly.

### Changed

- **Test count: 8,191 → 10,212** across the 18 suites, plus 673 platform tests. Green on a clean
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
- **`CLINICAL_REVIEW.md` holds 32 sections**, most of them open and awaiting a reviewer's
  ruling. The largest are MAGGIC's creatinine conversion and HF-duration direction (§12), the
  PhenoAge CRP convention (§8), Lung-RADS part-solid and non-solid nodules (§6f), and the
  disabled warfarin algorithm (§21). Settled in this release: EuroSCORE II (§1–2), CYP2D6 (§3),
  DPYD (§4) and the pharmacogenomics allele tables (§4b); §6e records the instruments checked
  that needed no change. They are recorded in the repo rather than silently resolved.
- **Warfarin dosing returns no dose.** `DosingCalculator.warfarin_dose` raises rather than
  returning a mis-transcribed one. This is a deliberate capability regression, and the single
  highest-value thing that would reverse it is Klein 2009's Supplementary Appendix 1.
- **PVS1 currently fires for no caller.** Correcting the ACMG implementation meant deleting a
  hand-estimated pLI list that had been standing in for a per-gene loss-of-function mechanism
  table. That is an accuracy fix that is also a capability regression, recorded as §30; the
  ClinGen Dosage Sensitivity map would close it.
- **Deterministic eligibility screening is a library API, not an HTTP endpoint.**
  `PatientMatchingWorkflow` in `src/clinical_workflows.py` evaluates criteria one by one with
  the tri-state logic above. `/v1/trials/match` does retrieval by text similarity and does not
  screen: it returns `criteria_evaluated: false` on every response and says so in its own
  docstring. No route was added for this release, because adding one is a feature rather than a
  correction, and the endpoint that exists no longer misdescribes itself.
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
