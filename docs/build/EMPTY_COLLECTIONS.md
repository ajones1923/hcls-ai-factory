# The 40 empty collections, and what can actually fill each one

**Status:** operational note, repo-only. **Measured:** 2026-11-03 against the live instance —
114 collections, **40 empty**, 4,764 vectors total.

`CORPUS_SEEDING.md` is the how. This is the what: a line per empty collection naming the source
that can fill it, whether a fetcher for that source already exists, and — for five of them —
why no public source can fill them at all.

The distinction matters because "40 empty" reads like one job. It is four jobs of very different
sizes, and one of them is not a data job but a decision.

## Why this is worth more than the vector count

`check_corpus.py` reports clinical-trial at 328 vectors across 14 collections: thin, working.
Ten of those collections are empty, and `scripts/check_corpus_weights.py` shows what that costs.
Coverage below is weight, computed the way each engine computes it — base weight × workflow
boost, renormalized — against live row counts:

| Subject | Workflows below 75% coverage | Worst |
|---|---:|---|
| cardiology | 12 of 12 | valvular_disease, **2%** |
| clinical-trial | 13 of 13 | recruitment_optimization, 24% |
| rare-disease-diagnostic | 13 of 13 | prenatal_diagnosis, 37% |
| single-cell | 11 of 11 | drug_response, 46% |
| cart | 2 of 5 | car_design, 38% |
| neurology | 5 of 10 | neuromuscular_eval, 70% |
| pharmacogenomics | 0 of 6 | — all at 100% |

**56 of 70 workflows retrieve from under three-quarters of the evidence they declare**, and
nothing anywhere in the retrieval path fails: a search of an empty collection succeeds and
returns zero rows, which is indistinguishable from a search that found nothing relevant.

## 0. cardiology: fourteen collections that were never created

**Not empty — absent.** The engine declares 12 workflows and 12 `cardio_*` collections
(`cardio_literature`, `_trials`, `_imaging`, `_electrophysiology`, `_heart_failure`, `_valvular`,
`_prevention`, `_interventional`, `_oncology`, `_devices`, `_guidelines`, `_hemodynamics`);
**none of them exists in Milvus.** No `cardio_*` collection has ever been created, so every
cardiology RAG query has only ever been able to retrieve from `genomic_evidence`, the shared
collection it borrows — 2% of its declared weight, and the only reason its coverage is not zero.
`check_corpus.py` could not show any of this, because that script can only iterate collections
that exist, and so printed no cardiology row at all.

This is the one item here that needs no sourcing decision — the corpus design is already
written in `core/engines/cardiology/src/vector_collections.py`. It needs
`scripts/setup_collections.py` run, then seeding.

Worth noting what this does *not* affect: the risk calculators (CHA₂DS₂-VASc, HAS-BLED, HEART,
ASCVD, EuroSCORE II, MAGGIC) are pure functions over their inputs and do not touch Milvus. The
engine's calculated output is sound; its *retrieved* output has always been empty.

## 1. Fillable now, with a fetcher that already exists (13)

`scripts/ingest_literature.py` (PubMed, any collection, schema-projected) and the
ClinicalTrials.gov v2 API are both working and verified on 2026-11-03.

| Collection | Source | Query / endpoint |
|---|---|---|
| `rd_case_reports` | PubMed | `"<disease>"[tiab] AND case reports[pt]` |
| `rd_natural_history` | PubMed | `natural history[tiab] AND <disease>` |
| `rd_registries` | Orphanet + PubMed | registry listings; `patient registry[tiab]` |
| `trial_results` | PubMed + CTG | `results[pt]`; CTG `hasResults=true` |
| `trial_rwe` | PubMed | `real-world evidence[tiab] OR observational` |
| `trial_adaptive` | PubMed | `adaptive design[tiab] AND clinical trial[pt]` |
| `trial_biomarkers` | PubMed | `biomarker enrichment[tiab] AND trial design` |
| `sc_methods` | PubMed | `single-cell[tiab] AND (method OR tool OR pipeline)` |
| `autoimmune_cross_disease` | PubMed | `overlap syndrome[tiab]`, shared-pathway terms |
| `rd_trials` | ClinicalTrials.gov | `query.cond=<disease>` |
| `cart_trials` | ClinicalTrials.gov | `query.intr=CAR-T`, `query.term=chimeric antigen` |
| `autoimmune_clinical_trials` | ClinicalTrials.gov | `query.cond=<autoimmune disease>` |
| `neuro_trials` | ClinicalTrials.gov | `query.cond=<neuro condition>` |

Four more are reachable from the same CTG records, which carry eligibility, outcomes, locations
and sponsor/investigator fields in the protocol section — these are **projections of one fetch**,
not four separate jobs: `trial_eligibility`, `trial_endpoints`, `trial_sites`,
`trial_investigators`. (`trial_sites` and `trial_eligibility` are the two heaviest empty
collections in the whole platform by retrieval weight: 0.25 each in their workflows.)

## 2. Fillable from a public source, needs a fetcher written (13)

| Collection | Source | License note |
|---|---|---|
| `rd_variants` | ClinVar (NCBI FTP / E-utilities) | public domain |
| `rd_genes` | NCBI Gene, HGNC, GenCC | public / CC0 |
| `rd_pathways` | Reactome | CC0 |
| `sc_pathways` | Reactome, MSigDB | Reactome CC0; MSigDB registration |
| `rd_newborn_screening` | RUSP (HHS) | US government, public domain |
| `trial_guidelines` | FDA guidance + ICH | **public domain** — unlike the clinical societies |
| `trial_safety` | CTG results + FAERS | public |
| `autoimmune_pgx_rules` | CPIC + PharmGKB | CPIC free; PharmGKB CC-BY-SA, attribution required |
| `cart_constructs` | FDA labels via DailyMed | public |
| `sc_datasets` | CELLxGENE Discover, HCA | public APIs, metadata only |
| `sc_drug_response` | GDSC, DepMap + PubMed | GDSC/DepMap free for academic use |
| `sc_spatial`, `sc_trajectories` | PubMed + method repos | abstracts only |
| `sc_clinical` | PubMed | abstracts only |

**Abstracts stay local.** PubMed abstracts are retrieved through E-utilities under NCBI's terms
and remain publisher-copyrighted. That is already this repo's rule — data and weights stay on the
box, only code and docs publish (`CLAUDE.md`) — so ingesting them locally is consistent, and
redistributing them in the public repo would not be.

Two entries need a licensing decision before anyone writes the fetcher: **OMIM** requires a
registered API key (free for academic use, licensed for other use), and **KEGG** restricts bulk
download. Neither is required — Orphanet, HPO and Reactome cover the same ground on open terms.

## 3. Curation, not fetching: the society criteria (4)

| Collection | What it wants |
|---|---|
| `autoimmune_disease_criteria` | ACR/EULAR classification criteria |
| `neuro_guidelines` | AAN / AHA / ILAE recommendations |
| `rd_guidelines` | ACMG / ACOG clinical genetics guidance |
| `neuro_electrophysiology` | EEG/EMG/NCS patterns and normal values |

These are not a fetch job. The published criteria *facts* — items, thresholds, score cut-points —
are citable and belong here with their primary source, exactly as the clinical reference-case
tests already cite CHA₂DS₂-VASc to Lip 2010. The *documents* are copyrighted and must not be
bulk-ingested.

So each record is a curated entry carrying its citation, which under this repo's rules
(`CLAUDE.md`: cite clinical claims to a primary source) is **Adam's call per instrument**, the
same standing as the open `CLINICAL_REVIEW` entries. Note that `neuro_guidelines` is the heaviest
empty collection in seven of the nine neurology workflows.

## 4. No public source can fill these — they need a labeled synthetic cohort (5)

| Collection | Fields that make it patient-level |
|---|---|
| `autoimmune_clinical_documents` | `patient_id`, `visit_date`, `provider`, `source_file` |
| `autoimmune_patient_labs` | `patient_id`, `collection_date`, `flag` |
| `autoimmune_patient_timelines` | `patient_id`, `event_date`, `specialty` |
| `imaging_reports` | `patient_id`, `study_date`, `findings_text`, `impression_text` |
| `imaging_radiomics` | `patient_id`, `study_date`, feature vectors per region |

There is no public corpus of individual patient records, and there should not be. These five can
hold one of exactly three things:

1. **Nothing** — and the agents that search them say so, which is what the grounding fixes on
   `fix/ten-out-of-ten` now make them do.
2. **A labeled synthetic cohort**, which the TSC program already does for 50 patients. Honest,
   demonstrable, and never to be presented as observed data.
3. **Real records under an IRB and a DUA**, which is a Cincinnati Children's conversation, not a
   seeding task.

**This is a decision, not a backlog item.** The autoimmune "diagnostic odyssey" workflow is built
on `autoimmune_patient_timelines`; the imaging agent's report retrieval is built on
`imaging_reports`. Until one of those three is chosen, those features retrieve nothing — and
option 1 is a perfectly defensible state for a research-use platform, as long as the UI says so.

## One thing to know before reading any weight table in this repo

Five subjects — clinical-trial, neurology, rare-disease-diagnostic, single-cell and cardiology —
carry **two** collection weight tables each:

| Table | Where | Read at retrieval time? |
|---|---|---|
| `WORKFLOW_COLLECTION_WEIGHTS` | `src/vector_collections.py` | **No** |
| `WORKFLOW_COLLECTION_BOOST` | `src/agent.py` | Yes |

The first is the one that looks authoritative: named for weights, fourteen entries summing to
1.00, exposed through a `get_search_weights()` accessor, with tests asserting the sum. Outside
those tests, `get_search_weights()` has no caller. What reaches Milvus is
`_get_boosted_weights`, which multiplies each collection's base weight from `COLLECTION_CONFIG`
by its `WORKFLOW_COLLECTION_BOOST` multiplier and renormalizes.

The two usually agree to within a few points, so this is not a defect — but it is a trap for
anyone reasoning about retrieval from the file that reads like the specification.
`scripts/check_corpus_weights.py` follows the live path, and a test in the clinical-trial suite
pins its arithmetic to `_get_boosted_weights` so the reimplementation cannot drift.

## The honest summary

- **1 subject** (cardiology, 12 collections) needs a script run that already exists.
- **17 collections** can be filled with tooling in the repo today.
- **13** need a fetcher written against an open source.
- **4** need per-instrument curation against primary sources — Adam's ruling.
- **5** cannot be filled with real data at all, and need a decision about synthetic or nothing.

Nothing in this file is a code defect. The code now *reports* the gap —
`scripts/check_corpus_weights.py` for the operator, `retrieval_coverage` on the response for the
reader, and `hcls_common.corpus_coverage` for anything else that needs to ask.

**Which subjects report it.** Six of the seven that declare a weight map now carry
`retrieval_coverage` and `collections_contributing` on their response, and append the disclosure
to the answer itself once coverage is degraded:

| Subject | Reports coverage | Note |
|---|---|---|
| clinical-trial | yes | |
| cardiology | yes | the 2% case |
| neurology | yes | |
| rare-disease-diagnostic | yes | |
| single-cell | yes | |
| cart | yes | verdict rides on `CrossCollectionResult`, prose rendered in `query()` |
| pharmacogenomics | not wired | selects collections from a list, no weights; measured at 100% coverage on all six workflows, so there is nothing to disclose |
| precision-autoimmune, precision-biomarker, clinical-imaging | n/a | no weight table; flat collection set |

The three with no weight table are the remaining gap in this mechanism. Two of them
(precision-biomarker at 0 empty collections, clinical-imaging at 2 of 12) have little to
disclose; precision-autoimmune has seven empty of thirteen, **three of which are the
patient-level collections in §4 that no public source can fill** — so its coverage number would
be a standing report of a decision not yet made rather than a data gap.
