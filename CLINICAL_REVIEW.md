# Clinical and scientific values awaiting review

Values in this repository that a review on 2026-11-03 could not settle against a primary
source. They are listed rather than changed, because a plausible-looking correction is the
same failure as a plausible-looking error. Each entry gives the current value, what the check
found, and the source to settle it with.

Nothing here is a clinical endorsement of the surrounding code. This platform is
decision-support for research use, not a medical device, and no value below has been reviewed
by a clinician in the course of this work.

**How to use this file:** settle an entry against its source, change the code with that source
cited in a comment next to the value, and delete the entry. If a value is confirmed correct,
delete the entry and add the citation anyway — the next reader will ask the same question.

**Prefer the maintained table over the guideline paper.** §3d is a correction that was correctly
sourced to a 2019 guideline paper and still wrong, because the issuing body had since revised the
table the paper printed. Where a body publishes a machine-readable table it maintains — CPIC
serves allele function, diplotype-to-phenotype and drug recommendations as JSON at
`api.cpicpgx.org` — cite **that**, with the retrieval date, and prefer checking a whole table
exhaustively to spot-checking rows. "Checked against its source" has a shelf life.

---

## 1–2. EuroSCORE II — SETTLED 2026-11-03, four coefficients were wrong

`core/engines/cardiology/src/risk_calculators.py`

The suspicion in this entry was correct and understated. Verified against the primary source
(Nashef SAM et al., "EuroSCORE II", Eur J Cardiothorac Surg 2012;41(4):734-745, **Table 6**,
the final logistic model) and cross-checked against two independent calculator
implementations. Four coefficients were wrong, not one:

| Field | Was | Table 6 |
|---|---|---|
| `creat_cc_lt_50` | 0.6521653 | **0.8592256** |
| `creat_cc_50_85` | 0.3150652 | **0.303553** |
| `lvef_31_50` | 0.1084443 | **0.3150652** |
| `lvef_21_30` | 0.3150652 | **0.8084096** |

The LV-moderate coefficient (0.3150652) had been copied into both the LV-poor and the renal
moderate slots. **Poor ventricular function was therefore scored with the moderate
coefficient — 0.3151 instead of 0.8084.** For a 72-year-old with poor LV, NYHA III, urgent,
two procedures, predicted operative mortality was 2.93% where it should be 4.83%. EuroSCORE II
decides whether to operate, so this understated surgical risk for exactly the patients in whom
the number is decisive.

The age term (§2) was also wrong, as suspected: it returned 0 at age ≤60 and
`coeff × (age − 60)` above. The published coding is a count that is **1 at age ≤60 and rises by
one per year after**, so every patient's linear predictor was short by one unit of 0.0285181.
Small, systematic, and in the direction of understating risk.

All 2,046 cardiology tests passed both before and after. Monotonic ordering was preserved, so
every output stayed plausible — which is why a suite of range assertions could not see it.
Pinned now in `tests/test_euroscore_ii_reference.py` (45 cases, every Table 6 coefficient).

⚠️ **Still open — three fourth-decimal differences.** The primary table was read once for
these and the repo may well be right; the effect of each is negligible (fifth significant
figure), so they were NOT changed on a single reading:

| Field | Repo | Single reading of Table 6 |
|---|---|---|
| `poor_mobility` | 0.2407818 | 0.2407181 (digits transposed either way) |
| `critical_preoperative_state` | 1.0862550 | 1.086517 |
| `urgency_salvage` | 1.3626972 | 1.362947 |

- **Settle with:** a second reading of Nashef 2012 Table 6, or euroscore.org's documentation.
- ⚠️ **Do not substitute the paper's UNIVARIABLE coefficients** (e.g. CC ≤50 at 1.688774,
  CC 50–85 at 0.6674962). Those belong to the derivation analysis, not the final model.

## 3. CYP2D6 activity score 1.0 — SETTLED 2026-11-03, and a second boundary was wrong too

**Resolved against the primary source**, which was retrieved rather than inferred this time:
the consensus paper itself, whose Discussion says in as many words **"reclassify AS of 1 from NM
to IM"**, and whose Table 3 gives PM = 0, **IM = 0 < AS < 1.25**, NM = 1.25 ≤ AS ≤ 2.25,
UM > 2.25.

> Caudle KE, Sangkuhl K, Whirl-Carrillo M, Swen JJ, Haidar CE, Klein TE, Gammal RS, Relling MV,
> Scott SA, Hertz DL, Guchelaar HJ, Gaedigk A. "Standardizing CYP2D6 Genotype to Phenotype
> Translation: Consensus Recommendations from the Clinical Pharmacogenetics Implementation
> Consortium and Dutch Pharmacogenetics Working Group." *Clin Transl Sci* 2020;13(1):116–124,
> Table 3. doi:10.1111/cts.12692, PMID 31647186.
> https://pmc.ncbi.nlm.nih.gov/articles/PMC6951851/

Reading the table rather than a summary of it surfaced a **second** error at the other end:

| Activity score | Was | Consensus | Who that is |
|---|---|---|---|
| 1.00 | Normal | **Intermediate** | `*1/*4`, `*1/*5`, `*41/*41` — common diplotypes |
| 2.25 | Ultrarapid | **Normal** | `*1x2/*41` — 2.25 is the inclusive *top* of Normal |

Clinically: an AS 1.0 patient was getting standard tricyclic dosing where CPIC asks for a 25%
reduction, and an AS 2.25 patient was getting the ultrarapid recommendation to change drug
class. Both ends corrected in `pgx_pipeline.py` (`PHENOTYPE_THRESHOLDS`) and `dosing.py` (the
TCA path carried its own copy of the same thresholds, with the same two errors), each with the
citation at the change.

The contradiction §3 was originally about is also closed: `knowledge.py` already stated
NM = 1.25–2.25, so the knowledge base and the code now agree.
`tests/test_cyp2d6_consensus_reference.py` holds **both** sides to the table, including the test
that previously asserted the wrong value "per CPIC 2023".

### 3a. `*41` allele activity value — OPEN, and probably harmless

While reading the consensus table: this repo scores `*41` at **0.5**, where CPIC's
allele-functionality table assigns **0.25**. `*10` is already 0.25 here, consistent with CPIC.

It appears to change no phenotype in practice, because the diplotypes it affects land in the
same bucket either way (`*41/*41` → AS 1.0 here vs 0.5 there; Intermediate in both), but it will
shift AS values the agent *reports*, and a `*1x2/*41` patient sits exactly on the NM/UM boundary
where 0.25 vs 0.5 decides 2.25 against 2.5.

- **Settle with:** the CPIC CYP2D6 Allele Functionality Table (the companion file to the paper
  above, not the genotype-to-phenotype table) — I did not retrieve it, so I did not touch the
  allele scores.
- **Scope if it changes:** `ACTIVITY_SCORE_TABLES["CYP2D6"]["allele_scores"]` in `knowledge.py`,
  `PhenotypeTranslator.ACTIVITY_SCORES` in `pgx_pipeline.py`, `_CYP2D6_ACTIVITY` in `dosing.py`
  — three copies of the same table, which is its own small problem.

### 3c. NUDT15 — SETTLED 2026-11-03, then PARTLY SUPERSEDED 2026-11-04 — see §3d

Found by acting on §3a's observation that the allele tables exist in three copies: the two
copies **disagreed**, and the one a clinician reads was the wrong one.

`PhenotypeTranslator.ACTIVITY_SCORES` scored NUDT15 `*2` at **0.5**. CPIC Table 1 lists `*2` as
a **no-function** allele. So `*2/*2` — two no-function alleles, a CPIC **poor metabolizer** —
summed to 1.0 and reported **"NUDT15 Intermediate"**. `dosing.py` had `*2` at 0.0 and got the
same patient right.

**Why it matters more than most:** a NUDT15 poor metabolizer reported as intermediate is given
30–70% of a standard thiopurine dose where CPIC says avoid azathioprine/6-MP or give 10% with
intensive monitoring. The failure mode is life-threatening myelosuppression, and NUDT15 PM is
not rare in East Asian and Hispanic populations.

Three further errors surfaced from reading Table 1 rather than summarising it:

| | Was | CPIC Table 1 |
|---|---|---|
| `*2` (translator) | 0.5, decreased | **no function** |
| `*4`, `*5` (dosing) | 0.0, no function | **unclear function** |
| `*6` (both) | 1.0, normal | **unclear function** |
| any unrecognised allele | 1.0, normal | — |

CPIC assigns **Indeterminate** to a diplotype carrying an unclear-function allele. An activity
score cannot express that, so both the translator and the dose calculator now refuse to produce
one, name the allele they could not call, and say "do not infer normal function from the absence
of a call". The same change removed `.get(allele, 1.0)` — *default to normal if unknown allele* —
which was the least safe default available for a gene whose unrecognised alleles are mostly
loss-of-function. TPMT got the same treatment in the same calculator: it names five alleles
where CPIC's table names dozens.

> Relling MV, Schwab M, Whirl-Carrillo M, et al. "Clinical Pharmacogenetics Implementation
> Consortium Guideline for Thiopurine Dosing Based on *TPMT* and *NUDT15* Genotypes: 2018
> Update." *Clin Pharmacol Ther* 2019;105(5):1095–1105, Table 1. doi:10.1002/cpt.1304,
> PMID 30447069. https://pmc.ncbi.nlm.nih.gov/articles/PMC6576267/

**The duplication itself is still open.** Allele activity values live in three places
(`pgx_pipeline.py`, `dosing.py`, `knowledge.py`) with nothing reconciling them, and it has now
produced two clinical defects — this one and §3's CYP2D6 thresholds, where *both* copies carried
the same two wrong boundaries. `tests/test_allele_tables_do_not_drift.py` fails the moment two
copies disagree about an allele they both define, with an allowlist of the asymmetries that are
real. That is a guard, not a fix; the fix is one table, imported, and it belongs in its own
change rather than riding along with a clinical correction.

### 3d. NUDT15 again — SETTLED 2026-11-04, §3c was built on a superseded table

§3c is right about `*2`, and the `.get(allele, 1.0)` removal it made is right. What it got wrong
is the *source*, and the error it introduced is the same harm §3c was written to prevent,
reached through a different door.

§3c sourced `*4`–`*9` to **Relling 2019 Table 1**, which in 2018 designated them "unclear
function", and so declared them unscoreable. **CPIC has since assigned them.** From CPIC's
maintained allele-functionality table:

| Allele | §3c treated as | CPIC current |
|---|---|---|
| `*4`, `*6`, `*7`, `*9` | uncertain → Indeterminate | **No function** |
| `*5`, `*8` | uncertain → Indeterminate | **Decreased function** |

and CPIC's diplotype table follows without ambiguity: `*1/*4` and `*1/*6` are **Intermediate**
Metabolizers, `*2/*4`, `*3/*6` and `*6/*6` are **Poor** Metabolizers, and `*1/*5` and `*1/*8`
are **Normal**. CPIC reserves Indeterminate for `*11`, `*12`, `*15`, `*16` and `*20`–`*23` —
none of which were in §3c's set, and all of which the "not in the allele table" branch was
already catching correctly.

**The harm.** A 9-year-old of Han Chinese ancestry starting 6-mercaptopurine maintenance for ALL,
NUDT15 `*3/*6` — a CPIC Poor Metabolizer — was told "No genotype-guided dose can be given …
NUDT15 `*6` has uncertain function in CPIC's guideline … **dose on clinical grounds** with the
standard CBC monitoring". For a thiopurine in a child that is a full weight-based dose, and what
follows is profound, prolonged myelosuppression. `*1/*4` and `*1/*6`, Intermediate Metabolizers
who should start reduced, got the same answer. The caveat text made it worse rather than better:
it was *factually wrong* as of CPIC's current table, and it read as a confident justification for
doing nothing.

**The general lesson, and it applies to every table in this file.** A guideline paper and the
table the issuing body maintains and republishes are different artefacts, and **for allele
function the maintained table governs**. A value "checked against its source" has a shelf life.
Where a maintained machine-readable table exists, cite that, with the retrieval date — which is
what §4b now does for every pharmacogenomics gene.

> CPIC NUDT15 allele-functionality table:
> https://api.cpicpgx.org/v1/allele?genesymbol=eq.NUDT15 (retrieved 2026-11-03)
> CPIC NUDT15 diplotype-to-phenotype table (276 rows, all reproduced):
> https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.NUDT15 (retrieved 2026-11-03)

**Also corrected here:** §3c's dose text. Its Poor Metabolizer branch said "AVOID thiopurines if
possible. If essential, reduce dose to 10% of standard (e.g., azathioprine 0.5 mg/kg/day thrice
weekly)". CPIC puts the 10-fold reduction and the thrice-weekly schedule on **mercaptopurine and
thioguanine**; its **azathioprine** row for a poor metabolizer at either locus reads, in full,
"Consider alternative nonthiopurine immunosuppressant therapy." And `min(tpmt, nudt15)` cannot
express a two-gene recommendation: CPIC gives a single intermediate locus 30–80% of standard and
a double heterozygote 20–50%, where the code gave both "30–70%, start at lower end" with a
structured field of 50 — a headline range running 20 points above CPIC's ceiling for the double
heterozygote, with titrate-upward advice attached.

> CPIC azathioprine recommendation table:
> https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm%3A1256 (retrieved 2026-11-03)
> CPIC mercaptopurine recommendation table:
> https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm%3A103 (retrieved 2026-11-03)

## 4. DPYD activity score 0.5 — SETTLED 2026-11-04, `pgx_pipeline.py` was right

`dosing.py` labelled DPYD AS 0.5 "DPD Intermediate Activity (Partial Deficiency)" with a 50%
reduction; `pgx_pipeline.py` called the same score Deficient; `knowledge.py` gave a third answer
again, "25–50% dose with intensive monitoring". Three copies, three answers.

CPIC settles it. Its fluorouracil and capecitabine recommendation rows for DPYD activity score
0.5, verbatim:

> "**Avoid** use of 5-fluorouracil or 5-fluorouracil prodrug-based regimens. In the event, based
> on clinical advice, alternative agents are not considered a suitable therapeutic option,
> 5-fluorouracil should be administered at a strongly reduced dose with early therapeutic drug
> monitoring."

CPIC calls AS 0.5 a **Poor Metabolizer** — its diplotype table assigns Poor Metabolizer to every
no-function + decreased-function pair — so `pgx_pipeline.py` ("DPD Deficient") was right and
`dosing.py` was wrong. All three copies now carry *avoid*.

The same retrieval settled the neighbouring row: CPIC gives activity score **1.5 the same 50%
starting-dose reduction as 1.0**, where `dosing.py` had 25% and `knowledge.py` had 25%. AS 1.5
is DPYD `*1`/c.2846A>T, the commonest actionable DPYD genotype in European-ancestry patients.

And DPYD HapB3 is present in the caller's allele set; CPIC lists it as the single allele
"c.1129-5923C>G, c.1236G>A (HapB3) | Decreased function | 0.5", which is how all three copies
now spell it. `knowledge.py` had listed c.1236G>A separately as "*5", so a patient genotyped for
both HapB3 markers scored two decreased-function alleles instead of one.

> CPIC fluorouracil recommendation table:
> https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm%3A4492 (retrieved 2026-11-03)
> CPIC capecitabine recommendation table:
> https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm%3A194000 (retrieved 2026-11-03)
> CPIC DPYD allele-functionality table:
> https://api.cpicpgx.org/v1/allele?genesymbol=eq.DPYD (retrieved 2026-11-03)

## 4b. Pharmacogenomics allele tables — SETTLED 2026-11-04 against CPIC's maintained tables

Every allele table in the pharmacogenomics agent was re-derived from
`api.cpicpgx.org/v1/allele`, and the phenotype thresholds were checked against
`api.cpicpgx.org/v1/diplotype` exhaustively rather than by sampling. The values and thresholds
now in `dosing.py` and `pgx_pipeline.py` reproduce CPIC's own phenotype for **all** 276 NUDT15,
1,225 TPMT, 3,570 DPYD, 2,850 CYP2C9, 351 CYP2C19, 36 UGT1A1, 21 CYP3A5 and 741 SLCO1B1
diplotypes built from alleles CPIC assigns a function to.

What was wrong, beyond §3c and §4 above:

| Gene | Allele | Was | CPIC |
|---|---|---|---|
| DPYD | `*7` (c.295_298delTCAT) | 0.5 decreased | **No function**, 0.0 |
| DPYD | `*10` (c.2983G>T) | 0.5 decreased | **No function**, 0.0 |
| DPYD | `*9A` (c.85T>C) | 0.5 decreased | **Normal function**, 1.0 |
| DPYD | `*4` (c.1601G>A, `knowledge.py`) | 0.5 decreased | **Normal function**, 1.0 |
| CYP2C9 | `*5` | 0.0 no function | **Decreased function**, 0.5 |
| CYP2C19 | `*9`, `*10` | 0.0 no function | **Decreased function**, 0.5 |
| UGT1A1 | `*37` | 0.0 no function | **Decreased function**, 0.5 |
| UGT1A1 | `*80` | 0.5 | **Unknown function** → Indeterminate |

Two threshold findings came out of the same retrieval and are worth recording, because neither
is guessable from a guideline PDF:

1. **TPMT and NUDT15 do not share a rule.** CPIC's diplotype tables make a no-function +
   decreased-function pair a NUDT15 **Poor** Metabolizer but a TPMT **Intermediate** one, and a
   decreased + normal pair **Normal** for both. One shared rule was applied to both genes.
2. **An increased-function allele does not rescue a reduced-function one.** CPIC calls UGT1A1
   `*28/*36` (Decreased + Increased) an **Intermediate** Metabolizer and CYP2C19 `*10/*17` a
   (Likely) **Intermediate** Metabolizer, where the summed activity score 2.0 reads Normal. An
   additive activity score cannot express that on its own, so the increased allele now
   contributes a normal 1.0 when paired with a reduced-function one.

### 4c. Open pharmacogenomics questions, not settled here

1. **SLCO1B1 `*17`.** `dosing.py` assigned phenotypes to `*1/*17` (Decreased Function) and
   `*5/*17` (Poor Function). `*17` appears in **neither** CPIC's SLCO1B1 allele-functionality
   table (which runs `*1`–`*58` but skips `*3`, `*7`, `*12`, `*13`, `*17`, `*18`, `*21`, `*22`,
   `*35`, `*38`, `*46`) nor its 741-row diplotype table. Both entries were removed rather than
   corrected — removing a fabricated value needs no source — so a `*17` genotype now returns
   Indeterminate. **Settle with:** PharmVar's SLCO1B1 allele definitions, then CPIC's SLCO1B1
   table, to confirm whether CPIC retired or renamed the haplotype. If `*17` is a real
   decreased-function haplotype, patients carrying it currently get "No recommendation" where a
   dose ceiling would be safer.
2. **HLA-A\*31:01 and carbamazepine severity.** `hla_screener.py` emits CONTRAINDICATED for this
   pair (FATAL and SEVERE both map to CONTRAINDICATED) where the knowledge-base entry's own text
   says "HIGH RISK". CPIC's carbamazepine rows in the API carry no phenotype keys, so the
   recommendation text could not be attributed to the allele. **Settle with:** Phillips EJ et
   al., "CPIC Guideline for HLA Genotype and Use of Carbamazepine and Oxcarbazepine: 2017
   Update", *Clin Pharmacol Ther* 2018;103:574–581, Table 2.
3. **CYP2D6 `*41`** — still §3a above, re-confirmed as 0.5 here against CPIC's 0.25 and
   deliberately not changed, because this agent's CYP2D6 table was otherwise exact and the
   change belongs with the §3a retrieval.
4. **Allele coverage.** The tables are now correct for every allele they name, and an allele
   they do not name returns Indeterminate rather than a normal-metabolizer guess. But coverage is
   still partial for the largest genes — CYP2D6 scores 21 of CPIC's >100 — so a patient with an
   uncommon CYP2D6 allele gets "No recommendation" rather than the phenotype CPIC would assign.
   **Settle with:** loading the tables from CPIC's data at build time rather than transcribing
   them, which is the refactor §3c's closing note describes.

## 5. Tuberous sclerosis — surveillance and treatment thresholds

`core/disease-programs/tuberous-sclerosis/`

| Value | Where | Question |
|---|---|---|
| AML treatment threshold 4.0 cm | `trajectory_config.yaml` | Consensus guidance for asymptomatic *growing* AML is first-line mTOR inhibitor at **>3 cm**; `trials.py` already uses ≥3 cm, so the two disagree |
| SEGA 1.8 cm "discussion threshold" | trajectory config | No guideline basis was found for a size threshold; management keys on growth and symptoms |
| Renal surveillance 12 months | `trajectory_config.yaml` | `itsc.py` uses 24 months, the corpus text says 1–3 years |
| Unconditional annual ophthalmology | `itsc.py` | Check against the surveillance table |

- **Settle with:** Northrup et al., *Updated International TSC Diagnostic Criteria and
  Surveillance and Management Recommendations*, Pediatr Neurol 2021
  (doi:10.1016/j.pediatrneurol.2021.07.011 — title confirmed via Crossref), and the 2012
  management recommendations it builds on.

## 5b. Tuberous sclerosis — five more uncited constants, found by the claim-integrity review

Opened 2026-11-04 by the per-subject review of `tuberous-sclerosis` (check 4). §5 above covers
the SEGA 1.8 cm and AML 4.0 cm thresholds and the surveillance cadence. These five are in the
same tree, carry the same clinical weight, and were not recorded anywhere.

| Value | Where | Question |
|---|---|---|
| SEGA growth slope **0.026 cm/month** (~3 mm/yr) | `trajectory_modeler/population.py` | The constant the whole trajectory prior rests on. Its name said `_PUBLISHED_FALLBACK` and no publication was named. Renamed `_ASSUMED_FALLBACK` pending a source |
| eGFR threshold **60.0** mL/min/1.73m² | `config/trajectory_config.yaml` | Matches KDIGO CKD stage G3a, but nothing in the file says so. Confirm and cite |
| Seizure frequency threshold **15.0 / month** | `config/trajectory_config.yaml` | No basis found. This fires a crossing alert; what guideline, if any, puts the line here? |
| Surveillance cadence multipliers **0.5 / 0.75 / 1.0**, floor **3 months** | `config/trajectory_config.yaml` | Shortens a real surveillance interval by up to half on a model-predicted crossing. Is a halved interval defensible, and is a 3-month floor the right clamp? |
| Mosaic VAF ceiling **0.30** | `variant_curator/runner.py` | Mosaicism call boundary. `src/eval/harness.py` separately uses `STANDARD_VAF_LIMIT = 0.10`; confirm the two are answering different questions |

**None of these were changed.** A clinical value moves only with a primary source cited in a
code comment, and there is none for any of them yet.

- **Settle with:** Northrup et al. 2021 (as §5) for surveillance and thresholds; KDIGO 2024
  CKD guideline for the eGFR boundary; for the SEGA growth rate, the TSC Alliance Natural
  History Database publications, or the EXIST-1/EXIST-2 placebo arms, which are the obvious
  published sources for untreated SEGA growth.

**Related, and not a clinical question:** the trajectory prior reported itself as an "NHD
mixed-effects" fit with `population_informed: true`, while being fitted to 24 simulated
patients generated from the 0.026 constant above. That was a false claim rather than a wrong
value, so it was fixed immediately rather than recorded here — see
`tests/unit/test_trajectory_prior_is_honest.py`. The constant itself stays open, above.

## 6. Tuberous sclerosis — ACMG criteria application

`src/variant_curator/annotation.py`: PVS1 is downgraded for low VAF, and PS1 is used for "same
variant known pathogenic". PS1 is defined as the same *amino acid change* as an established
pathogenic variant with different nucleotide change, so this reads as a misuse. The check
suggested it changes no classification in the current cohort, which is worth confirming.

- **Settle with:** Richards et al., ACMG/AMP 2015, Table 3 (criteria definitions). The
  combining rules in `acmg.py` were checked against Table 5 and matched.

## 6b. Tuberous sclerosis — PP4 now applies the diagnostic criteria, with one gap

`src/agents/variant_curator/annotation.py`

PP4 used to be granted for any recorded phenotype at all, so a lone nonspecific seizure
earned a supporting pathogenic criterion. It now requires the recorded features to meet the
TSC clinical diagnostic criteria — 2 major, or 1 major plus ≥2 minor (Northrup 2021).

What could not be settled: three of the major features are defined with counts — ≥3
hypomelanotic macules, ≥3 angiofibromas, ≥2 angiomyolipomas — and the cohort's phenotype
records carry an HPO term with no count, so presence is all the code can read. A patient
with one hypomelanotic macule is currently counted as having the major feature.

- **Settle with:** the Northrup 2021 criteria table, then either add counts to the phenotype
  records (`src/cohort/clinical.py`) or mark those three features as count-unverified in the
  criterion's rationale.
- **Also:** the `HPO` dict in `src/cohort/clinical.py` is commented "Canonical TSC
  major-feature codes" but includes seizure, infantile spasms, intellectual disability and
  autism, none of which are TSC diagnostic criteria. The comment, not the codes, is wrong.

## 6c. Coronary CTA — the HRP modifier threshold

`core/engines/clinical-imaging/agent/src/workflows/ct_coronary_angiography.py`

The HRP modifier used to fire on a **single** high-risk plaque feature. It is now gated at two
(`HRP_MIN_FEATURES = 2`), which is what two independent secondary sources state: HRP replaced
the term "vulnerable plaque" and kept its definition of **at least two** of low-attenuation
plaque (<30 HU), positive remodeling, spotty calcification and the napkin-ring sign.

Two things were **not** settled:

- The primary modifier table (Cury et al., CAD-RADS 2.0, J Cardiovasc Comput Tomogr 2022 —
  already cited by this engine's test suite) was not read for the change.
- Both secondary sources qualify the rule as ≥2 features **in the same plaque**. This workflow
  cannot evaluate that: `high_risk_plaque_features` is a case-level list with no plaque
  identity, so two features found in two different plaques currently count as HRP.

The demo case carries three features, so its headline (`CAD-RADS 4A/P3/HRP`) is unchanged and
the 2026-08-13 accuracy pass still holds.

- **Settle with:** the CAD-RADS 2.0 modifier table. If plaque identity matters, the feature
  list needs to become per-plaque before the rule can be applied properly.
- **Secondary sources checked:**
  https://radiologyassistant.nl/cardiovascular/cad-rads/coronary-artery-disease-reporting-and-data-system
  and https://radiopaedia.org/articles/cad-rads

## 6d. Radiology scoring — two edges the published charts do not define

Both of these were checked against their charts on 2026-11-03 and the point tables, size bands
and management thresholds all matched. These are the two edges the charts leave open, now
written into the code rather than left as silent choices.

**ACR TI-RADS, a total of exactly 1 point.** The chart names 0, 2, 3, 4–6 and 7+. One point is
reachable (mixed cystic/solid composition with no other suspicious feature) and is not in the
list. `points_to_tirads` treats it as TR2. Management is identical either way — neither 1 nor 2
points warrants FNA or follow-up at any size — so nothing turns on it clinically, but it is a
choice and should be a documented one.

**ACR Lung-RADS v2022, baseline only.** `classify_solid_nodule_lung_rads` takes a diameter and
nothing else, and it is now applied to solid nodules only (§6f). Lung-RADS also keys on exam context: at follow-up, a NEW 4–6 mm solid nodule is
category 3 rather than 2, a new 6–8 mm nodule is 4A rather than 3, and a growing nodule moves
up a category. So on a follow-up LDCT this can under-call by one category. The fix is to give
the workflow a prior exam to compare against, not to change the thresholds.

- **Settle with:** Tessler FN et al., J Am Coll Radiol 2017;14(5):587-595 (TI-RADS); the ACR
  Lung-RADS v2022 assessment categories.

## 6e. Risk scores verified against their publications (no action needed)

Recorded so the next reader knows what has been checked rather than re-deriving it. Reference
cases live in `core/engines/cardiology/tests/test_score_reference_cases.py`,
`core/engines/clinical-imaging/agent/tests/test_radiology_reference_cases.py` and
`core/disease-programs/tuberous-sclerosis/tests/unit/test_acmg_table5_reference.py`.

| Instrument | Source | Result |
|---|---|---|
| CHA2DS2-VASc | Lip, Chest 2010 | point table and both age bands correct; one defect found and fixed (every comorbidity defaulted to absent, so the score came from age and sex alone) |
| HAS-BLED | Pisters, Chest 2010 | point table correct, including renal/liver and drugs/alcohol scoring separately and E being age >65; same defaulting defect found and fixed |
| HEART | Six 2008 / Backus 2013 | bands correct; one defect found and fixed (missing troponin scored as normal) |
| ACR TI-RADS | Tessler 2017 | all five feature tables, level mapping and FNA size thresholds correct |
| ACR Lung-RADS v2022 | ACR | solid-nodule size bands correct (baseline; see §6d). **Part-solid and non-solid tables are NOT implemented and must not be guessed — see §6f** |
| ACMG/AMP combining rules | Richards 2015, Table 5 | every row correct, including conflicting-evidence → VUS |
| ASCVD Pooled Cohort Equations | Goff 2013, Appendix 7 + guideline vignette | all four cohorts' coefficients correct; vignette reproduced within 0.09 pp |
| EuroSCORE II | Nashef 2012, Table 6 | **4 coefficients wrong, corrected** (see §1–2); 3 fourth-decimal values still open |

**ASCVD note for anyone comparing against the ACC/AHA calculator:** the guideline's vignette
(55-year-old, TC 213, HDL 50, SBP 120 untreated, non-smoker, non-diabetic) gives 5.3% for white
men and this platform displays **5.4%**. The individual coefficient sum is 60.6995; the
guideline prints 60.69, which is that truncated, and truncating through to the end reproduces
5.3%. The full-precision path is the more defensible one, so this is a presentational
difference of one tenth of a point, not drift. The other three cohorts display exactly the
published 6.1% / 2.1% / 3.0%.

Still unverified: **MAGGIC** and **PhenoAge** (§8). Both are regression/point models needing
their published tables and a worked example — a case invented from the repo's own coefficients
would only pin current behaviour.

## 6f. ACR Lung-RADS v2022 — OPEN: no part-solid or non-solid table is implemented

`core/engines/clinical-imaging/agent/src/workflows/ct_chest_lung_nodule.py`

**What was wrong.** `classify_solid_nodule_lung_rads` — a function whose name, docstring and
cited source all say *solid* — was being handed every nodule the workflow detected, with only
the choice of *which* diameter varying by density. Lung-RADS v2022 grades part-solid nodules on
criteria that differ from the solid table and non-solid (ground-glass) nodules on a different
scale again, so the solid table was wrong in both directions. Measured on the release
candidate:

| Nodule | Reported | Consequence |
|---|---|---|
| part-solid, 25 mm total, **9.3 mm solid component** | **4A**, "3-month LDCT or PET/CT" | the classic invasive-adenocarcinoma appearance, and the exact nodule shipped as the demo case; three months of growth is stage migration |
| part-solid, 18 mm, 6.5 mm solid component | **3**, 6-month LDCT | under-called |
| pure ground-glass, 20 mm | **4B, CRITICAL**, "tissue sampling recommended" | a biopsy of an indolent lesion that belongs on annual screening |

**What was changed (2026-11-03).** Nodule density is now explicit and required
(`resolve_nodule_type`), and `classify_nodule_lung_rads` grades a density only if this module
holds its published table. It holds **solid only**. A part-solid nodule, a ground-glass nodule,
an unstated density and an unrecognised density are each **refused** — `LungRADSGradingError`,
no category, no screening interval — and the workflow reports the study
`WorkflowStatus.FAILED`, severity `UNKNOWN`, with the ungraded nodule carried as an explicit
indeterminate finding. The demo payload deliberately still contains the part-solid nodule, so
the gap is visible in the demo rather than hidden by it.

**The open question — the thresholds themselves. Do not guess them.** No part-solid or
non-solid band is written into this code, because the issuing body's own document could not be
retrieved: the ACR asset
(`https://www.acr.org/-/media/ACR/Files/RADS/Lung-RADS/Lung-RADS-2022.pdf`) 302-redirects to a
Sitecore tenant page and serves no document. Reconstructing the bands from secondary sources
would be the same failure as the defect, one step further from the source.

- **Settle with:** Christensen J et al. (ACR Lung-RADS Committee), "ACR Lung-RADS v2022:
  Assessment Categories and Management Recommendations", *Chest* 2024;165(3):738-753 (also
  *J Am Coll Radiol* 2024;21(3):473-488), **Table 1** — the issuing body's own table
  reproduced in the literature.
- **Then:** add `classify_part_solid_lung_rads` and `classify_non_solid_lung_rads` with those
  thresholds cited at the values, extend `_IMPLEMENTED_DENSITY_TABLES`, and add one reference
  case per density to
  `core/engines/clinical-imaging/agent/tests/test_radiology_reference_cases.py`. Until then
  the refusal is the correct behaviour, not a placeholder.

## 7. CAR-T knowledge base — approval dates and REMS

`core/agents/cart/src/knowledge.py` dates Yescarta's ZUMA-7 indication 2021-03-05 and ZUMA-5
2024-03-15. A check of FDA's approval history found April 2021 actions and an approval letter
dated 2022-04-01, and no action on 2021-03-05. Separately, REMS requirements for approved CAR-T
products were removed on 2025-06-26, so any REMS text in this corpus is likely stale, and
Aucatzyl is absent from the product list.

- **Settle with:** FDA's Yescarta approval history
  (https://www.fda.gov/BiologicsBloodVaccines/CellularGeneTherapyProducts/ApprovedProducts/ucm581222.htm)
  and the current prescribing information.

## 8. PhenoAge — the CRP unit convention is worth 2.40 years, and no source settles it

`core/agents/precision-biomarker/src/biological_age.py`

**This is the entry to rule on.** The module applies `ln(CRP in mg/L)`. The alternative reading
is `ln(CRP in mg/dL)`, which is mg/L ÷ 10. Because the term is `0.0954 × ln(CRP)`, the choice
shifts the linear predictor by a constant `0.0954 × ln(10) = 0.2197` — measured as **exactly
2.40 years on the reported biological age, at every CRP level**. A healthy 50-year-old reads
46.2 under the current convention and 43.8 under the other.

Checked 2026-11-03 and **not settled**:

- Levine et al. 2018 (PMID:29676998), its published **correction**, and Liu et al. 2018 all
  print the equation with **no units for any biomarker**. I read all three.
- `dayoonkwon/BioAge`, which this module cites as the source of its coefficients, takes a
  pre-computed `lncrp` and states no unit either.
- **For mg/dL:** several independent public calculators explicitly divide mg/L by 10 before the
  log; NHANES, the derivation cohort, reports CRP in mg/dL.
- **For mg/L (current):** the other three converted terms here are SI (g/L, μmol/L, mmol/L),
  and mg/L is the SI unit for CRP, so a consistent-SI reading supports the current behaviour.

It was **not changed**, because changing a clinical value on an unsettled source is the same
error as leaving a wrong one, only harder to spot later. Instead the assumption is now stated in
the module docstring, carried in every result as `crp_unit_convention`, and pinned by tests that
fail if the 2.40-year figure drifts.

- **Settle with:** the Levine 2018 supplementary material, or the NHANES III lab documentation
  for the CRP variable the model was fit on. One line in `UNIT_CONVERSIONS` switches it.
- **Fixed in passing:** the alkaline phosphatase coefficient was 0.0019; the published
  correction prints **0.00188**. Changed (effect is in the fourth decimal of the predictor).
- **Still open from the original entry:** `report_generator.py` advises avoiding synthetic folic
  acid on an MTHFR variant, which runs against ACMG and CDC positions on MTHFR testing; and
  "Function Health optimal" ranges appear beside clinical reference ranges.

## 9. Structure and chemistry facts

These are database lookups rather than clinical judgments, and the ones below were corrected
against the primary record on 2026-11-03. Listed here so the next reader can see what was
checked:

- **CB-5083 SMILES** — the repo's string parses to C28H32N6O2, which is not CB-5083. PubChem
  CID 73051434 gives C24H23N5O2, MW 413.5,
  `CC1=CC2=C(C=CC=C2N1C3=NC4=C(COCC4)C(=N3)NCC5=CC=CC=C5)C(=O)N`
  (1-[4-(benzylamino)-7,8-dihydro-5H-pyrano[4,3-d]pyrimidin-2-yl]-2-methylindole-4-carboxamide).
- **PDB 5FTK** — described in the repo as X-ray, 2.3 Å, "VCP + CB-5083 inhibitor". The RCSB
  record is *Cryo-EM structure of human p97 bound to ADP*, electron microscopy, 2.4 Å. Related
  entries, all cryo-EM: 8OOI 2.61 Å, 9DIL 3.3 Å, 7K56 3.9 Å.

Any remaining per-entry resolution/method/ligand claims in
`core/agents/precision-biomarker/app/protein_viewer.py` and the gene→PDB mappings in
`core/engines/precision-intelligence/src/knowledge.py` should be checked the same way, one
entry at a time, against https://data.rcsb.org/rest/v1/core/entry/<ID>.

## 10. TSC corpus — prevalence figures in the TAND entries

`core/disease-programs/tuberous-sclerosis/src/rag/corpus.py`

Both TAND entries had DOIs that resolved to unrelated papers; those are corrected (the TOSCA
TAND findings paper and de Vries 2015's TAND Checklist paper, both verified via Crossref).
What was *not* verified is the text each entry carries:

- "TAND features are present in ~90% of TSC patients … 30-50% of features go unrecognized in
  routine care."
- "The TAND consensus framework defines six clusters … and the TAND-L lifetime checklist."

The cluster list and the checklist are consistent with the de Vries work as far as this check
could tell, but the percentages were not traced to a figure or table.

- **Settle with:** the two papers now cited in that file. If a figure does not support a
  number, change the text rather than the citation.

## 11. Demo VCF — HGVS internal consistency

`data/demo/demo_variants.vcf` pairs `c.463C>T` with `p.Arg155His` (a C>T at that position gives
Arg→Cys; R155H needs c.464G>A) and `c.1234G>A` with `p.Arg412His` (a G at the codon's first
base codes Gly, not Arg). A third entry, `c.476G>A` / R159H, is self-consistent. The file is
labelled synthetic and not for clinical use, so this is a correctness issue in demo data rather
than a patient-facing one.

- **Settle with:** the transcript sequence for the gene in question, then correct either the
  cDNA or the protein notation.

## 12. MAGGIC — creatinine bands converted with the wrong factor, and an extra EF point

`core/engines/cardiology/src/risk_calculators.py`

Checked against the published integer score (Pocock SJ et al., Eur Heart J
2013;34(19):1404-1413) via an independent calculator implementation; the primary Table was not
retrieved, so **nothing here was changed**. Two findings, the first of which is plain arithmetic
and does not depend on the band values at all.

**1. Creatinine thresholds use ÷100 where the unit conversion is ÷88.4.** The published bands
are in **μmol/L**: <90, 90–109, 110–129, 130–149, 150–169, 170–209, 210–249, ≥250 scoring
0,1,2,3,4,5,6,8. This code works in mg/dL with thresholds 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.5 —
i.e. the μmol/L numbers divided by **100**, not by the 88.4 μmol/L per mg/dL conversion. The
6-point edge is also 190 where the published band starts at 210. Net effect, computed:

| Creatinine | = μmol/L | repo pts | published pts |
|---|---|---|---|
| 0.9 mg/dL | 79.6 | 1 | **0** |
| 1.0 | 88.4 | 1 | **0** |
| 1.2 | 106.1 | 2 | **1** |
| 1.5 | 132.6 | 4 | **3** |
| 1.9 | 168.0 | 6 | **4** |
| 2.8 | 247.5 | 8 | **6** |

So renal function is over-scored by 1–2 points across most of the clinical range, including for
a perfectly normal creatinine of 1.0 mg/dL.

**2. An extra point for EF 40–44%.** Published EF points are <20→7, 20–24→6, 25–29→5, 30–34→3,
35–39→2, **≥40→0**. This code adds a `<45 → 1` band, so every HFmrEF patient with EF 40–44%
carries one point the published score does not give.

Combined, a typical EF 42% patient with creatinine 1.0 mg/dL scores **2 points above** the
published MAGGIC, and the score maps directly onto 1- and 3-year mortality estimates.

**3. Not verified at all:** the systolic BP bands, and whether SBP points are stratified by EF
in the published score the way age points are. This code stratifies age by EF but not SBP.

- **Settle with:** Pocock 2013 Table 2/Table 4, or the calculator at heartfailurerisk.org.

**4. The HF-duration term reaches the calculator now, but its direction is still inverted —
SETTLED by the source, fix NOT applied here.** `api/routes/cardio_clinical.py` used to send
`extra={"hf_duration_18m": ...}` while `calculate_maggic` reads `extra["hf_duration_months"]`,
so the duration fell back to 0 on every API call and the "HF diagnosed >18 months ago" checkbox
had no effect on the score. The route now sends `hf_duration_months`, which is what makes the
remaining defect visible: `risk_calculators.py:1458-1461` reads

```python
# HF first diagnosed within 18 months = +2 (higher risk for new diagnosis)
hf_new_pts = 2 if not hf_duration_gte_18 else 0
```

The published model's variable is **"HF duration >18 months"** — longer-standing heart failure
carries the HIGHER risk — so the two points belong to the chronic patient, not the newly
diagnosed one. The magnitude (2) is right; only the direction is wrong.

- **Primary source:** Pocock SJ et al., *Eur Heart J* 2013;34(19):1404-1413, **Table 2**
  (multivariable model, all 39,372 patients): "HF duration >18 months — rate ratio 1.188
  (1.139, 1.240), log rate ratio 0.173, Z 7.96, P<0.0001". The integer rule is in the same
  paper's Methods (p.1405): "the score increases by an integer amount for each risk factor
  level above the lowest risk ... making log rate ratio 0.1 equivalent to 1 point", so
  0.173 → 2 points, awarded for duration >18 months.
- **Why it was not changed here:** `src/risk_calculators.py` is owned by a parallel change in
  this release and was deliberately not touched. Until the direction is corrected there, a
  patient with HF for six years scores 2 points BELOW the published MAGGIC (the mirror image of
  the pre-fix behaviour, which charged every API-scored patient the +2 "new HF" point). 20% is
  the threshold at which this function emits "consider referral to advanced heart failure /
  transplant centre", so the sign matters at the point of referral. The regression test added
  for the route (`tests/test_cardio_safety_regressions.py::TestA2MAGGICDurationReachesTheCalculator`)
  asserts only that the input reaches the calculator and changes the score, so it is
  direction-agnostic and will stay green either way — **it does not prove the direction**.

---

## 13. CRT Class 2a QRS lower bound — converged on 120 ms, full text not retrievable

`core/engines/cardiology/src/hf_thresholds.py`

`GDMTOptimizer` used 120 ms and `HeartFailureWorkflow` used 130 ms for the QRS duration at
which CRT becomes Class 2a in LBBB, both citing the 2022 AHA/ACC/HFSA HF guideline. A patient
with LVEF 32%, QRS 125 ms and LBBB therefore got a Class 2a CRT recommendation from
`/gdmt/optimize` and nothing at all from `/workflow/heart-failure`. The two surfaces now share
one constant, `CRT_QRS_CLASS_IIA_MS = 120`.

- **Why 120:** both surfaces cite the 2022 AHA/ACC/HFSA guideline, and the only retrievable
  reproductions of that guideline's own recommendation text put the Class 2a LBBB band at
  120–149 ms ("For patients who have LVEF <=35%, sinus rhythm, LBBB with a QRS duration of 120
  to 149 ms, and NYHA class II, III, or ambulatory IV symptoms on GDMT, CRT can be useful...",
  COR 2a, LOE B-NR), corroborated independently by the 2023 HRS physiologic-pacing guideline's
  2a (B-R) recommendation for the same band. 130 ms is the **2021 ESC** figure (McDonagh TA et
  al., *Eur Heart J* 2021;42:3599-3726), not a value of the guideline the code cites.
- **What is unsettled:** the publisher's full text could not be retrieved from this
  environment — ahajournals.org and jacc.org both return HTTP 403 — so this rests on
  reproductions of the recommendation rather than on the article.
- **Settle with:** Heidenreich PA et al., 2022 AHA/ACC/HFSA Guideline for the Management of
  Heart Failure, *Circulation* 2022;145:e895-e1032, device-therapy recommendation table
  (section 7.4.1), recommendations 8 and 13.
- **A second question in the same table, NOT acted on:** the same reproductions restrict the
  **non-LBBB 120–149 ms** band to Class 2b *and to NYHA III or ambulatory IV*.
  `_assess_device_therapy` gives Class IIb for non-LBBB QRS >=120 at NYHA II as well, which
  would over-recommend a device. Settle with the same table before changing it.

---

## 14. Cardiology defaults and gaps found during the November 2026 fix pass, left as they are

`core/engines/cardiology/`

Behaviours noticed while fixing the audited defects, each outside the scope of that fix list and
none of them changed:

1. **`_parse_nyha("unknown")` and `_classify_nyha("unknown")` both return NYHA II**
   (`src/gdmt_optimizer.py`, `src/clinical_workflows.py`), and NYHA II is the gate for ICD
   eligibility and for ivabradine, so an unparseable functional class silently becomes
   device-eligible. Two tests assert this value
   (`tests/test_gdmt_optimizer.py:156`, `tests/test_clinical_workflows.py:331`), so changing it
   is a deliberate behaviour change, not a bug fix. The ICD path no longer asserts eligibility
   on absent waiting periods, which blunts the consequence but does not remove it.
2. **Every comorbidity flag across this engine defaults to `False`** (CHA2DS2-VASc, HAS-BLED,
   `bp_treatment` in the preventive workflow, the EuroSCORE II booleans), so "not stated" and
   "absent" are the same value. This is the audit's finding A-9 and is not confined to
   cardiology.
3. **SGLT2i initiation consults no renal function.** `_recommend_pillar_initiation` now refuses
   to start a beta-blocker, an ARNI or an MRA without the relevant vital or lab, but the SGLT2i
   branch reads nothing at all, so it does not invent a value — it simply never checks one.
   Dapagliflozin and empagliflozin both carry eGFR thresholds in their labels. Settle with the
   DAPA-HF and EMPEROR-Reduced labels before adding a guard.
4. **A beta-blocker already at target is not re-checked against the current vitals.**
   `_check_contraindications` is consulted only on the titration path, so a patient at
   carvedilol 25 mg BID with HR 48 and SBP 82 gets no "hold" advice from `/gdmt/optimize`.
5. **`_assess_current_gdmt` resolves a bare `"metoprolol"` to `metoprolol_succinate`** via
   substring matching (`src/gdmt_optimizer.py`), which would score an immediate-release
   tartrate patient as being on a pillar drug that is not indicated in heart failure.
6. **`SeverityLevel` has no "unknown" member.** Workflows whose inputs were not supplied now
   report `INFORMATIONAL` with the indeterminacy stated in the findings, because the alternative
   was to map an unmeasured patient onto a disease severity. A dedicated level would be better.

## 21. Warfarin IWPC — DISABLED 2026-11-04, the coefficient column is shifted and the error is race-specific

**This is the most harmful defect the pre-release audit found, and the function no longer returns
a dose.** `DosingCalculator.warfarin_dose` now raises; the Streamlit page shows the reason and
tells the clinician to dose by protocol with INR monitoring.

Every number in the block *is* an IWPC coefficient, attached to the wrong variable:

| Coefficient | Published variable (Klein 2009) | What the code multiplied it by |
|---|---|---|
| −0.2188 | CYP2C9 genotype **unknown** | `cyp4f2_1_3` |
| −0.2760 | **Black or African American** | `cyp4f2_3_3` |
| +1.1816 | **enzyme inducer status** | `african_american` |
| +1.2799 | *not an IWPC coefficient* | `enzyme_inducer` |
| −0.2029 | *not the published amiodarone value* (−0.5503) | `amiodarone` |
| +0.2107 | *not an IWPC coefficient* — IWPC has no smoking term | `smoker` |

IWPC also has **no CYP4F2 term at all**, and the block is missing all three "genotype unknown"
terms (VKORC1 unknown −0.4854, CYP2C9 unknown −0.2188, missing/mixed race −0.1032).

**Why it had to be switched off rather than left with a note.** A Black patient picks up
**+1.1816** on the square-root-of-weekly-dose scale where the algorithm applies **−0.2760**.
That is ~1.5–1.7× the predicted dose, for one patient group only, in a narrow-therapeutic-index
anticoagulant, in the direction that causes a bleed. Nothing downstream would look wrong: the
number is plausible, the units are right, and the result is returned under the name
"IWPC Pharmacogenetic Warfarin Dosing Algorithm".

**Why it is not simply corrected here.** This repo changes a clinical value only against a
primary source in hand. The mapping above is the structural argument — those numbers are
IWPC's, in the wrong slots, which is enough to know the code is wrong without knowing the whole
table. It is *not* enough to retype the table from memory. Supplementary Appendix 1 of Klein
2009 and the warfarindosing.org workbook both carry it; neither could be retrieved that night.

- **Settle with:** Klein TE et al., "Estimation of the Warfarin Dose with Clinical and
  Pharmacogenetic Data", *N Engl J Med* 2009;360:753–764, **Supplementary Appendix 1 §1** (the
  CPIC warfarin guideline supplement reproduces it verbatim).
- **Then:** restore all coefficients and add the three unknown-genotype terms, replace
  `tests/test_dosing.py::TestWarfarinDose` with **reference cases** (a described patient in, the
  published dose out), and delete the invented CYP4F2 and smoking parameters from the signature
  or implement them from a cited source.
- **Also fix:** the agent's own research paper and its readiness report each state a *different*
  wrong version of this algorithm, with the coefficient column shifted by one row — the same
  transcription error, propagated into the documentation.

### What the old tests were doing

Eighteen tests covered this function and **every one asserted only a direction** — "this input
moves the dose down". A coefficient attached to the wrong variable cannot fail a direction test.
Three of them pinned the defect as the specification, by name:

```python
def test_african_american_increases_dose(self):
    assert aa["predicted_weekly_dose"] > other["predicted_weekly_dose"]
```

plus `test_smoker_increases_dose` and `test_cyp4f2_effect`, which assert two terms the model
does not contain. This is the EuroSCORE II blind spot (§1–2) exactly: range and direction
assertions cannot see a wrong coefficient, which is the argument for reference cases.

---

## 22. `precision-oncology/agent/src/knowledge.py` — 51 per-gene evidence levels and indication lists, declared unverified

`ACTIONABLE_TARGETS` (2,053-line module, 51 gene entries). Each entry carries an
`evidence_level` and a `cancer_types` list, and `therapy_ranker.py` used the latter as the
gate deciding whether a therapy is "indicated" for the patient's tumour. **No part of this
table had ever been checked against a label or a guideline.**

### What a spot-check of the twelve highest-consequence entries found (2026-11-03)

Checked against the drugs' own FDA labels, retrieved from DailyMed (section 1 INDICATIONS AND
USAGE). The citation for each is on the matching `VERIFIED_DRUG_INDICATIONS` entry in
`therapy_ranker.py`.

| Gene | `cancer_types` (level) | Verdict against the labels |
|---|---|---|
| **EGFR** | NSCLC, head and neck, colorectal, glioblastoma (A) | **Over-broad.** All six listed drugs (osimertinib, erlotinib, gefitinib, afatinib, dacomitinib, amivantamab) are NSCLC-only for an EGFR driver. Head and neck, colorectal and glioblastoma are unsupported; the colorectal entry conflates anti-EGFR *antibodies* given on RAS-wild-type status with TKIs given on an activating mutation. |
| **ROS1** | NSCLC, cholangiocarcinoma, glioblastoma (A) | **Over-broad.** Every ROS1 label is NSCLC-only. Worse, `lorlatinib` is in the drug list and has **no ROS1 indication at all** (LORBRENA is ALK-positive NSCLC). |
| **MET** | NSCLC, renal cell, hepatocellular, gastric (A) | **Over-broad.** capmatinib is MET-exon-14 NSCLC only; `crizotinib` is in the list but its label carries ALK and ROS1 indications, not MET. |
| **PIK3CA** | breast, endometrial, cervical, head and neck, colorectal (A) | **Over-broad.** PIQRAY is PIK3CA-mutant HR+/HER2− breast cancer with fulvestrant; the other four are unsupported. |
| **KRAS** | NSCLC, colorectal, pancreatic (A) | **Over-broad.** sotorasib covers NSCLC and mCRC (with panitumumab); adagrasib's retrieved label is NSCLC only; pancreatic is unsupported by either. |
| **ALK** | NSCLC, anaplastic large cell lymphoma, neuroblastoma (A) | **Partly wrong.** NSCLC correct for all five drugs; ALCL correct for crizotinib only; **neuroblastoma is unsupported**, and crizotinib's third indication — inflammatory myofibroblastic tumour — is missing. |
| **BRAF** | melanoma, NSCLC, colorectal, thyroid, hairy cell leukemia (A) | **Nearly right as a union, wrong per drug.** Melanoma/NSCLC/colorectal/thyroid are each covered by *some* drug; **hairy cell leukemia is not in any of the three labels** (ZELBORAF's second indication is Erdheim-Chester disease). TAFINLAR explicitly carries a Limitation of Use **excluding colorectal cancer**, which the gene-level list cannot express. |
| **HER2** | breast, gastric, NSCLC, colorectal, endometrial (A) | **Correct as a union** (ENHERTU covers breast, gastric, NSCLC and HER2 IHC3+ solid tumours; TUKYSA covers breast and RAS-wild-type CRC) — but no single drug covers the list. |
| **NTRK** | tissue-agnostic + sites (A) | **Correct.** larotrectinib, entrectinib and repotrectinib all carry site-agnostic NTRK-fusion indications. |
| **BRCA1 / BRCA2** | breast, ovarian, prostate, pancreatic (A) | **Correct as a union** (LYNPARZA covers all four) — wrong per drug: RUBRACA has no breast or pancreatic indication, ZEJULA is ovarian only, TALZENNA is breast and prostate. |
| **IDH1** | glioma, AML, cholangiocarcinoma (A) | **Correct** (TIBSOVO: AML, MDS, cholangiocarcinoma; VORANIGO: grade-2 astrocytoma/oligodendroglioma). MDS is missing from the list. |
| **RET** | NSCLC, medullary thyroid, papillary thyroid (A) | **Correct, and narrower than the label** — RETEVMO also has a site-agnostic RET-fusion indication. Erring narrow is the safe direction. |

**The systematic pattern:** `cancer_types` is a *gene-relevance* list, not an indication list,
and it is the wider of the two. Seven of the twelve were over-broad; of the five that were
defensible as a union, three were wrong drug-by-drug. So the defect is in the mechanism, not in
individual rows, and that is what was fixed: the gate now decides per **(drug, driver)** pair
against `VERIFIED_DRUG_INDICATIONS`, which carries the label citation and the Limitations of
Use. A pair with no retrieved label reports `indication_unverified=True` and a note naming the
document needed — the old behaviour, now declared.

### What is still unverified, and what would settle it

- **The remaining 39 `ACTIONABLE_TARGETS` entries' `cancer_types` lists and all 51
  `evidence_level` values.** Nothing in this work verified an *evidence level*: the levels were
  left as found and only *lowered* where a label contradicted the indication.
- **`THERAPY_MAP` (30 drugs) `approved_indications` free text**, `RESISTANCE_MAP`,
  `PATHWAY_MAP`, `BIOMARKER_PANELS` (20 panels, with their own `evidence_level` and
  `cancer_types`), `ENTITY_ALIASES` and `PEDIATRIC_DOSING` — none read in this pass.
  `PEDIATRIC_DOSING` is the most consequential of those: it is dosing.
- **Settle with:** each drug's FDA label via Drugs@FDA (`accessdata.fda.gov/scripts/cder/daf/`)
  or DailyMed, section 1; the NCCN guideline for the tumour; and OncoKB's own published level
  definitions (`oncokb.org/levels`) for what each letter is supposed to mean — the A–E ladder in
  this repo is not defined anywhere in it.
- **Then:** move each verified pair into `VERIFIED_DRUG_INDICATIONS` with its citation. The
  `indication_unverified` flag makes the remaining surface countable; it should reach zero.

### Open questions this pass deliberately did not answer

1. **Platinum chemotherapy in HRD/BRCA** is emitted at level B on "BRCA/HRD-positive tumors
   show enhanced platinum sensitivity". No platinum label selects patients on BRCA or HRD
   status, so there is no label to gate it with. Left at level B, flagged
   `indication_unverified`. *Settle with:* the NCCN ovarian/breast guidelines, or the trial that
   the level B is meant to rest on.
2. **Atezolizumab for TMB-H** was emitted at level B on "evidence supports activity"
   (IMvigor211, POPLAR). The TECENTRIQ label does not mention tumour mutational burden at all,
   so the gate now treats it as not indicated and drops it two tiers. *Open:* what tier a
   literature signal with no labelled indication deserves. The two-tier drop is a ranking
   policy chosen here for consistency, not a value read from a source.
3. **dostarlimab's MSI-H vs dMMR scope.** JEMPERLI §1.3 covers **dMMR** recurrent or advanced
   solid tumours site-agnostically, while its MSI-H indication (§1.2) is locally advanced rectal
   cancer. The engine treats MSI-H and dMMR as one biomarker, so dostarlimab is gated as
   site-agnostic. *Settle with:* the JEMPERLI label §1.1–1.3 plus whether the engine should
   separate MSI-H from dMMR at all.
4. **The `BIOMARKER_PANELS` loop in `_identify_biomarker_therapies` is dead code.** No panel
   carries the `marker`, `threshold` or `drugs` keys it reads, so it has never emitted a
   recommendation. It has been gated anyway so it cannot reproduce finding A-8 when those keys
   appear, but someone should decide whether it is meant to work at all.

## 23. `precision-intelligence/src/knowledge.py` — 201-gene knowledge graph, first read 2026-11-03

2,684 lines, never read by any previous audit. Two structures:

- **`KNOWLEDGE_CONNECTIONS`** — 201 genes across 13 therapeutic areas, each with `protein`,
  `function`, `pathway`, `diseases`, `drugs`, `drug_status`, `pdb_ids`, `druggable`. It is
  injected into the LLM's context by `rag_engine._build_prompt` via
  `format_knowledge_for_prompt`, so every field is read to the clinician through the model.
- **`GENE_REFERENCE_DATA`** — 73 genes with `uniprot_id`, `reference_smiles`, `reference_drug`.
  These seed the drug-discovery stage: a wrong SMILES or a wrong structure silently generates
  chemistry against the wrong target.

**Good news, verified so the next reader need not re-check it.** Unlike the oncology table,
this module asserts **no evidence levels and no indications** — it is descriptive. And its
identifiers check out:

- **73 / 73 `uniprot_id` values are correct** (UniProtKB REST, 2026-11-03). 64 resolve to the
  human gene they are filed under; the other nine are pathogen targets correctly mapped to
  their organism's accession — `GYRA`→P0AES4 (*E. coli* DNA gyrase A), `INHA`→P9WGR1 (*M.
  tuberculosis*), `HIV1_PR`/`HIV1_RT`/`HIV1_IN`→P03367/P03366/Q76353 (HIV-1 Gag-Pol),
  `HCV_NS3`/`HCV_NS5B`→P26664/Q9WMX2, `SARS2_MPRO`/`SARS2_RDRP`→P0DTD1 (both cleaved from
  pp1ab, so one accession for two targets is right).
- **23 / 23 `pdb_ids` checked resolve to the correct protein** (RCSB PDB REST, polymer-entity
  UniProt cross-reference, for VCP, C9orf72, MAPT, EGFR, KRAS, BRCA1, ALK, CFTR, PCSK9, TTR,
  CYP2D6, TPMT). The demo structures are right: 5FTK, 7K56 and 8OOI all map to **P55072**, VCP.

**What is NOT verified, and what would settle it:**

- **`drug_status`** is a regulatory claim in free text ("FDA approved", "Clinical development",
  "Phase I") for 201 genes. Not one was checked. A stale "FDA approved" is the same class of
  error as the oncology indication lists. *Settle with:* Drugs@FDA / DailyMed per drug, and
  ClinicalTrials.gov for the phase claims.
- **`drugs`** — the per-gene drug lists, including compound codes with phases (e.g.
  "CB-5083 (Phase I)"). CB-5083's development status in particular should be re-checked before
  any demo narration rests on it.
- **`diseases`, `pathway`, `function`, `protein`** — descriptive prose read into the LLM
  context. Plausible on inspection, unchecked. *Settle with:* UniProt entry function sections
  and OMIM / MONDO for the disease associations.
- **`reference_smiles`** (73 structures) — not parsed or compared with the named
  `reference_drug` in this pass. *Settle with:* PubChem / ChEMBL canonical SMILES for each
  `reference_drug`; an RDKit round-trip plus an InChIKey comparison would make this a test
  rather than a review item.
- `src/target_hypothesis.py`, `src/llm_client.py` and `src/milvus_client.py` remain unread.

## 30. Rare disease — PVS1 needs a per-gene loss-of-function MECHANISM table

`core/agents/rare-disease-diagnostic/src/decision_support.py`,
`ACMGVariantClassifier.classify` (PVS1 branch).

**Settled on 2026-11-03 for the generic defect; the replacement table is OPEN.**

Richards 2015 Table 3 (PMC4544753, retrieved) conditions PVS1 on *"a gene where loss of
function (LOF) is a known mechanism of disease"*. The code used a hand-written
`_LOF_INTOLERANT_GENES` set commented *"pLI > 0.9 in gnomAD"*. Gene-level LOF **constraint**
is a different concept from disease **mechanism**, and the list contained genes whose
mechanism is explicitly not LOF: *FGFR3* (achondroplasia is gain-of-function; *FGFR3* LOF
causes CATSHL, a different disorder), *HTT* (CAG expansion; *HTT* LOF is not Huntington
disease), *BRAF*/*RAF1*/*PTPN11* (activating RASopathies), *MYH7* (dominant-negative
missense).

The list is deleted. PVS1 now requires the caller to assert `lof_is_mechanism` for the gene and
otherwise **refuses to apply the criterion**, recording "whether loss of function is a known
disease mechanism for <gene> was not asserted" in the evidence summary. That is correct and
conservative, but it means PVS1 is never applied by any current caller, which is a capability
regression even though it is an accuracy improvement.

**What is open:** the curated per-gene mechanism table that would let PVS1 fire again.

**Source to settle it with:** the ClinGen Dosage Sensitivity curations
(https://search.clinicalgenome.org/kb/gene-dosage) give a haploinsufficiency score per gene,
which is the closest curated expression of "LOF is a known mechanism"; a score of 3
("sufficient evidence for haploinsufficiency") is the defensible threshold. Where a ClinGen
Variant Curation Expert Panel specification exists for the gene it governs instead. Neither was
retrieved in this work.

**Specifically for *MYH7*, also open:** whether PVS1 may be applied at all. Settle against
Kelly MA et al., *"Adaptation and validation of the ACMG/AMP variant classification framework
for MYH7-associated inherited cardiomyopathies"*, Genet Med 2018;20:351–359 (PMC5563454) — the
ClinGen Hypertrophic/Dilated Cardiomyopathy VCEP specification. Not retrieved. The generic
Table 3 defect above stands without it.

---

## 31. Rare disease — the revised Ghent nosology family-history criteria are NOT implemented

`core/agents/rare-disease-diagnostic/src/clinical_workflows.py`,
`ConnectiveTissueWorkflow.execute`.

Line 1348 was a bare expression statement — `ghent.get("family_history", False)` — so the value
was computed and bound to nothing, and **every family-history criterion of the nosology was
unimplemented while the code reported "Ghent criteria NOT fully met"**. A 16-year-old with
bilateral ectopia lentis whose father has genetically confirmed Marfan syndrome, aortic root not
yet dilated, got "Ghent criteria NOT fully met (aortic Z=0, …)" and no aortic surveillance,
beta-blocker/ARB, ophthalmology referral or activity restriction. The `aortic_z` default of `0`
additionally asserted a **normal aortic root for a patient whose aorta was never measured**.

Both missing-value defects are fixed (an unmeasured Z-score is `None` and reported as not
evaluable; `family_history` is assigned and, when positive, the workflow states that the
family-history arm is NOT machine-evaluated, recommends evaluating it by hand before excluding
Marfan syndrome, and raises a cardiology trigger).

**What is open:** implementing the criteria.

**Source to settle it with:** Loeys BL, Dietz HC, Braverman AC, et al., *"The revised Ghent
nosology for the Marfan syndrome"*, J Med Genet 2010;47(7):476–485, PMID 20591885 — the
diagnostic-criteria table. **jmg.bmj.com returns HTTP 403 to automated fetch and the paper was
not retrieved.** The only reachable version was an *adapted* table in a secondary review (Pepe
G et al., Appl Clin Genet 2016;9:55–65, PMC4869846, Table 1), and a diagnostic criterion is not
transcribed from a secondary summary — that is the §3c mechanism of this very file.

Three criteria are needed: ectopia lentis + FH; systemic score ≥7 + FH; aortic dilatation + FH
at an age-stratified Z threshold. Note that the secondary table also indicates the age
stratification applies **only** to the family-history criterion, which would mean the existing
no-FH branch's flat Z ≥ 2 is correct — but that too needs the primary source before anything is
changed on its basis.

Also open from the same block: whether criterion "EL + FBN1" requires an *FBN1* variant
**previously associated with aortic disease** (the code requires only `fbn1_mutation`).

**Removed, needing no source:** the citation `"ACMG Practice Resource: Marfan syndrome and
related disorders (2021)"` carried no author, journal or DOI, could not be identified in the
ACMG practice-resource index, and nothing in the workflow derived from it. Removing a
fabrication is not a clinical change.

---

## 32. Rare disease — PM2 strength, and the status of PP5/BP6

`core/agents/rare-disease-diagnostic/src/decision_support.py`.

PM2 is applied at Moderate, which is Richards 2015 Table 3 as published. The ClinGen Sequence
Variant Interpretation Working Group has recommended downgrading it to **PM2_Supporting**, and
has separately considered retiring PP5 and BP6. This branch changed the PS1 branch to **PP5 at
Supporting** (correct per Table 3 as retrieved) and did not touch PM2's strength.

**OPEN.** Settle against the ClinGen SVI recommendations index,
https://clinicalgenome.org/working-groups/sequence-variant-interpretation/, for (a) the
PM2_Supporting downgrade and (b) the current status of PP5/BP6. Not retrieved. Needed before
any change to PM2's weight or any removal of PP5/BP6 — and note that under the new
combinatorial classifier a strength change alters the label, so this is not cosmetic.

Still **not implemented at all**, and therefore structurally unreachable rule sets in Table 5:
BS3, BS4, PP2, BP2, BP5. With only BS1 and BS2 available, Table 5's ">=2 Strong benign ->
Benign" rule can be reached only one way.

---

## 33. Rare disease — criteria still granted on evidence weaker than Table 3 defines

`core/agents/rare-disease-diagnostic/src/decision_support.py` and
`src/clinical_workflows.py`. **OPEN — not addressed on this branch (audit §24).**

These four were out of scope here and remain as the audit found them. They matter more now,
not less, because the combining rule is correct: a wrongly-granted criterion now lands on the
right rule set and changes the label cleanly.

- **PP4** (`clinical_workflows.py`, passed into the classifier as `phenotype_specific`) is
  granted when the gene appears in **any** disease in the 21-entry `_HPO_DISEASE_DB` sharing
  **one** HPO term with the patient. Table 3 requires a phenotype *"highly specific for a
  disease with a single genetic etiology"*. A single shared HP:0001250 (Seizures) — present in
  MSUD, PKU, Dravet, galactosemia, OTC deficiency and MMA in that same table — satisfies the
  code, and the table's own Long QT entry lists three genes, so "single genetic etiology" fails
  by construction.
- **PM1**: `"TP53": [(125, 300)]` treats the entire 175-residue DNA-binding domain as a
  mutational hot spot. Table 3 requires a hot spot *"without benign variation"*; the TP53
  DNA-binding domain contains benign variation.
- **BP7**: `if synonymous and not splice_impact` — no conservation test. Table 3 BP7 ends
  *"AND the nucleotide is not highly conserved"*. The clause is missing from the implementation
  **and** from both documentation tables (`src/knowledge.py` ACMG_CRITERIA and the
  `/v1/diagnostic/acmg-criteria` route).
- **BS1**: `if freq > 0.01` — one blanket 1% threshold for every disorder. Table 3 BS1 is
  *"Allele frequency is greater than expected for disorder"* and the agent's own documentation
  table says *"(using disease-specific thresholds)"*. A *GBA1* founder variant at 2% in
  Ashkenazi Jews (carrier frequency roughly 1 in 15 for Gaucher disease) earns BS1. Under the
  new classifier BS1 alone is correctly a VUS rather than Likely benign, so the harm is reduced
  — but BS1 plus any one Supporting benign criterion is still Likely benign.

Sources: Richards 2015 Table 3 (PMC4544753) for all four definitions; the disease-specific
threshold requirement is reinforced by the gene-specific modification mechanism in Ghosh et al.
2018 (PMC6188666). For per-disease BS1 thresholds, the ClinGen VCEP specifications are the
practical source.

---

## 34. Rare disease — HPO information content is hand-estimated, not derived

`core/agents/rare-disease-diagnostic/src/decision_support.py`, `HPOToGeneMatcher`.

Fixed on this branch: an HPO term with no frequency on file used to default to `0.001`, giving
any unknown string — including `"banana"` and six terms from the class's own gene map — an IC of
9.97 against 2.74 for HP:0001250 (Seizures). The matcher's highest-weighted evidence was
systematically its least-known evidence. `_compute_ic` now returns `None`, contributes 0.0 to
similarity, and names the term in `unscored_terms`.

**Still open:** `_TERM_FREQUENCIES` is 43 hand-written *estimates* (the comment says so), not
annotation frequencies, and `_bma_similarity` is an IC-weighted set **intersection** — exact
membership, no ancestor, descendant or MICA logic — while the class docstring describes
best-match-average semantic similarity. Scores are also divided by the gene's annotation count,
so sparsely annotated genes outrank well annotated ones.

Settle by computing IC from real annotation frequency over the loaded ontology
(https://ontology.jax.org/api/hp/terms/ and the HPO annotation files), or rename the method to
what it computes. Until an ontology is loaded, the ranking reflects which terms the author
added to a 43-entry dict.

---

## 35. Clinical trial — a criterion no rule parses is now UNKNOWN, which is a capability gap

`core/agents/clinical-trial/src/clinical_workflows.py`,
`PatientMatchingWorkflow._evaluate_criterion`.

Fixed: verdicts are tri-state and an unparsed criterion is `None`, never met. This is correct —
ClinicalTrials.gov publishes `eligibilityCriteria` as free-text `markup` with no
machine-readable per-criterion representation (API v2 `GET /studies/metadata`, retrieved
2026-11-03), so an unparsed criterion has no verdict.

**The consequence to be aware of before release:** on a realistic oncology trial most criteria
are now "not evaluated", because only six dimensions are implemented (age, sex, diagnosis,
biomarker, comorbidity, medication). ECOG performance status, laboratory values (platelets,
haemoglobin, creatinine clearance, hepatic function), prior therapies, measurable disease per
RECIST and washout periods are all unparsed. A match score is **withheld entirely** when more
than half the criteria are unevaluable. The workflow is now honest but much less able to score
a trial, and that threshold (0.5) is an engineering choice with no clinical source behind it.

Three further gaps, all explicitly reported in the evidence string rather than guessed:

- The structured CTG eligibility gates — `minimumAge`/`maximumAge` (`NormalizedTime`), `sex`
  (`Sex`), `stdAges` (`StandardAge[]`) — are ingested by
  `src/ingest/clinicaltrials_parser.py` and **never read by the matcher**, which still regexes
  free text. Audit **F20** also notes `sex` defaults to the permissive `"All"` when absent.
  Both are unfixed here.
- Diagnosis matching requires whole-token containment rather than coded concepts. A partial
  overlap is reported `None`, so a legitimate synonym ("NSCLC" vs "non-small cell lung cancer")
  is not matched. Settle by normalising against MeSH/ICD-10/SNOMED or at minimum the CTG
  `conditionsModule.conditions` array the parser already stores.
- `genomic_variants` is collected on `PatientProfile` and never read; the biomarker rule
  compares the gene symbol, not the variant, so "EGFR exon 19 deletion" and "EGFR T790M" are
  indistinguishable to it.

---

## 36. Clinical trial — exclusion severity qualifiers are handled by keyword, not by concept

`core/agents/clinical-trial/src/clinical_workflows.py`,
`PatientMatchingWorkflow._evaluate_criterion`.

Fixed: comorbidity matching now requires every content word of the patient's comorbidity to
appear in the criterion as a whole word, which drops the two measured false positives
("controlled hypothyroidism" tripping "severe **un**controlled autoimmune thyroid disease", and
"seasonal allergic rhinitis" tripping "history of **allergic** reaction to platinum agents").

A severity qualifier in the criterion (`severe`, `uncontrolled`, `active`, `clinically
significant`, …) is now handled **asymmetrically**: a patient record that contradicts the
qualifier does not trip the criterion, while a record merely silent on it returns `None` and
asks for clinician review. That asymmetry is a deliberate safety choice — clearing a patient
whose severity is unrecorded is the unsafe direction — but it is a judgement, not a sourced
rule, and the qualifier list is hand-written.

**Open:** code comorbidities and medications to SNOMED CT / RxNorm concepts and compare concepts
rather than words. Until then "rifampin" and "rifampicin" are different strings to this matcher
(the test `test_medication_substring_does_not_trip_an_exclusion` pins that as the *current*
behaviour, and it is a false negative on an exclusion — the unsafe direction).

---

## 37. Clinical trial — remaining honesty and table defects not addressed on this branch

All **OPEN**, all confirmed by the audit, all out of scope for `fix/rare-trial`.

- **F12** — four incompatible phase vocabularies and no CTG mapping.
  `ProtocolDesignWorkflow.preprocess` still silently warns a CTG `"Phase 1"` into `phase_ii`
  because its `valid_phases` set lacks `early_phase_1` and `not_applicable`. The *scoring* path
  no longer has a permissive default (an unmapped phase now refuses to report a probability of
  success), but the vocabularies still disagree. The CTG `Phase` enum is exactly `NA`,
  `EARLY_PHASE1`, `PHASE1`, `PHASE2`, `PHASE3`, `PHASE4`, and a combined study is the
  two-element array `["PHASE1","PHASE2"]` (`GET /studies/enums`, retrieved 2026-11-03). There
  is no combined enum value.
- **F19** — `/v1/trial/workflow/{type}` is wired to `api.main._WorkflowEngine`, an LLM
  passthrough, so `PatientMatchingWorkflow` and the other deterministic workflows are
  **unreachable over HTTP**. The tri-state fix in this branch therefore does not yet reach any
  HTTP caller. `api/main.py:209` also calls `.get()` on `TrialSearchResult` dataclasses and
  swallows the `AttributeError`, so `context` is always `""` and `evidence_used` always `False`
  while `status` reports `"completed"`.
- **F11** — `get_search_weights`, `_DEFAULT_SEARCH_WEIGHTS` and `WORKFLOW_COLLECTION_WEIGHTS`
  are dead; `PRODUCTION_READINESS_REPORT.md` documents the dead table as live; six live
  workflows have no entry in `WORKFLOW_COLLECTION_BOOST` and run unboosted; and `search()`
  recomputes weights from the base table, so the workflow boost never reaches ranking at all.
- **F16** — `_RESTRICTIVE_PATTERNS` matches `"age"` as a bare substring, so "St**age** IIIB"
  and "m**anage** oral medication" are reported as narrow age ranges at 18% population impact
  each. (The *matcher*'s age rule is fixed to whole-word matching on this branch; the
  eligibility-optimisation pattern list is not.)
- **F18** — `src/cross_modal.py` reports `"N agent(s) consulted with no safety concerns"` when
  agents were unreachable. An absence of data reported as an absence of risk.
- **F21/F22/F25** — two disagreeing protocol-complexity tables with fabricated percentiles
  (`percentile_rank = complexity_score * 85 + 10`, so a protocol of zero complexity is reported
  at the 10th percentile, against no distribution at all); `src/decision_support.py` has **no
  production importer** and its six classes are dead while the live paths reimplement three of
  their tables with different numbers; and response models advertise fields that are always
  empty, including `/dct/plan` returning `{"status": "evaluated", "overall": "feasible"}`
  unconditionally.
- **F24** — site scoring gives a site with no data a real, rankable score from six clinical and
  operational defaults, and the `diversity_index` default of 0.3 sits exactly on the `< 0.3`
  low-diversity test, so an unknown-diversity site reads as acceptable.
- **F26** — secondary endpoints are stored but never embedded, so a semantic search for a
  secondary endpoint cannot retrieve the trial; `outcomesModule.otherOutcomes` is dropped; and
  `filterable_fields` is advertised in docstrings with no reader anywhere.

### Guideline citation strings still unverified

The ICH and FDA titles attached to the remaining workflows were not retrieved. On this branch
the two workflows I touched cite only what they use (`PatientMatchingWorkflow` cites nothing,
`ProtocolDesignWorkflow` and `DiversityAssessmentWorkflow` cite their one retrieved source), but
the others still carry fixed lists. Two look wrong on their face:

- *"FDA Guidance: Broadening Eligibility Criteria (2019)"* — I am not aware of an FDA guidance
  by that title. FDA's are the five **"Cancer Clinical Trial Eligibility Criteria"** guidances
  (final July 2020) and **"Enhancing the Diversity of Clinical Trial Populations — Eligibility
  Criteria, Enrollment Practices, and Trial Designs"** (final November 2020). Settle against
  https://www.fda.gov/regulatory-information/search-fda-guidance-documents.
- *"ICH E6(R3) GCP — Eligibility Assessment"*, *"— Investigator Qualifications"*, *"—
  Technology in Clinical Trials"* — whether sections by these names exist. Settle against the
  ICH E6(R3) Step 4 document (adopted 6 January 2025) table of contents at
  `database.ich.org`. Also `/v1/trial/guidelines` lists **ICH E20** as available; E20 (adaptive
  designs) was at Step 2 draft, not Step 4, as of my knowledge cutoff — settle against
  https://www.ich.org/page/efficacy-guidelines.

### Tables with no retrieved source

- **Tufts CSDD protocol-complexity benchmarks** (`src/decision_support.py`) — medians/p75/p90
  for procedure, visit, endpoint, criteria and amendment counts, attributed to "Industry
  benchmarks (Tufts CSDD)" naming no issue. Settle against the Tufts CSDD *Impact Report* issue
  on protocol design complexity; not retrievable without a subscription. The code must name the
  volume, number and year.
- **CYP450 allele frequencies by ancestry** (`src/clinical_workflows.py` `_CYP_ANCESTRY`) —
  settle against the CPIC allele-frequency table supplements (`cpicpgx.org/guidelines/`) or
  PharmGKB biogeographical-group tables, after mapping the agent's ancestry labels
  (`european/african/east_asian/south_asian/hispanic`) onto CPIC's biogeographical groups.
- **Population-impact values** (`_POPULATION_IMPACT`, e.g. `ecog 0-1` 0.25, `hemoglobin` 0.12,
  `creatinine clearance` 0.15) — the live and dead copies agree on overlapping keys, which
  means only that they were copied from each other. Provenance unknown.
- **FDA approval dates** in `src/ingest/regulatory_parser.py` `REGULATORY_MILESTONES` — all
  individually plausible, none verified. Settle against Drugs@FDA.
- **BIO rare-disease phase transition rates** — `rare_disease` was removed from
  `_HISTORICAL_SUCCESS` because Figure 2 does not break the area out. BIO analyses rare disease
  in a separate section of the same report, which was not read. Settle there before
  reintroducing a rare-disease row.

### Also noted: 4 wrong NCT IDs in the OTHER LANDMARK_TRIALS copy

`src/knowledge.py:LANDMARK_TRIALS` is the better copy and is what `seed_knowledge.py` should
eventually read, but it has its own errors, confirmed by the audit and **not corrected here**
because nothing imports it: `"ADVANCE (Aducanumab)"` → NCT02477800 is **ENGAGE** (221AD301, and
there was no aducanumab Phase 3 named ADVANCE); `"MAGELLAN (Zolbetuximab)"` → NCT03504397 is
**SPOTLIGHT**; `"CREST (Semaglutide in Obesity)"` → NCT03693430 is **STEP 5** (registry
enrollment 304 against a claimed 1,961, which is STEP 1's number); `"PANORAMIC"` → NCT04998396
is a 26-participant AAV9 gene-therapy study against a claimed 26,411 and Phase 3. `SUNFISH` is
claimed Phase 2/3 where the registry says `PHASE2`. Several enrollment figures across the table
are the randomized rather than the enrolled count.

The CI gate added in `core/agents/clinical-trial/tests/test_landmark_trial_ids.py` covers only
the parser copy. Pointing it at both tables is a ~15 s check for 76 IDs at the unkeyed 3 req/s
limit.

---

## 38. Both agents — FHIR, LOINC and Phenopacket conformance

`core/agents/rare-disease-diagnostic/api/routes/reports.py`.

**Settled on this branch:** `DiagnosticReport.status`. The resource was emitted `"final"`, whose
normative definition is *"The report is complete and verified by an authorized person"* — false
on its face for an automated tool. It is now `partial` (*"data in the report may be incomplete
or unverified"*) when there is a conclusion and `registered` (*"the existence of the report is
registered, but there is nothing yet available"*) when there is not, and an empty `conclusion`
is omitted rather than emitted as `""`. Source: HL7 FHIR R4 `DiagnosticReportStatus` value set,
https://hl7.org/fhir/R4/valueset-diagnostic-report-status.html, retrieved 2026-11-03.
`preliminary` was deliberately **not** chosen: its definition, *"verified early results are
available"*, also asserts verification.

**Still OPEN:**

- Whether **LOINC 81247-9** ("Master HL7 genetic variant reporting panel") is the correct panel
  code for this content. Settle against https://loinc.org/81247-9/. Not retrieved.
- The **GA4GH Phenopacket v2** required-field list, in particular whether `MetaData.resources`
  may be empty. This branch emptied it rather than assert an HPO ontology version
  (`"2024-04-26"`) for an agent that loads no ontology — but an empty `resources` array may
  itself be non-conformant, in which case the right fix is to load the ontology, not to omit the
  declaration. Settle against https://phenopacket-schema.readthedocs.io/. Not retrieved.
- `interpretations` is always `[]`, so a Phenopacket from a diagnostic agent carries no
  interpretation.
- The agent advertises `v1.0.0` in `app/diagnostic_ui.py`, `api/routes/reports.py` and
  `src/export.py` in a v2.2.0 release.

---

## 39. Rare disease — remaining defects in scope-adjacent code, not addressed

All **OPEN**, confirmed by the audit, out of scope for `fix/rare-trial`.

- **§12** `/diagnose` adds a gene-match bonus of `+5` to the *term-overlap* counter and then
  reports that counter as `phenotype_overlap` and as `confidence`: one phenotype term and one
  gene symbol returned *"Gaucher Disease, confidence 0.95, phenotype_overlap 5.0, evidence
  ['Matched on 5 term(s)']"* with **zero** phenotype terms matched. The same route computes
  `diagnostic_odyssey_risk` from the number of lines the clinician typed.
- **§15** `/phenotype/match` claims "semantic similarity", splits an HPO ID into `["hp",
  "0001250"]` and substring-searches those tokens against disease free-text vocabulary. Seizures
  plus global developmental delay returns `matches: []` with `total_diseases_screened: 48`,
  where the correct answer from the agent's own table was seven candidates. The IC-weighted
  matcher in `src/decision_support.py` is wired to no route.
- **§18** every Pathogenic/Likely-pathogenic call publishes a `critical_alert` onto the
  cross-agent bus with no minimum-evidence gate and no `research_use_only` flag. The corrected
  classifier fires far fewer of these (a 20-case probe fired 12 pre-fix, including from two-key
  inputs), but the gate is still absent.
- **§19** the dysmorphology matcher scores **negated** features as present (`"no cleft palate"`,
  `"absent hypertelorism"` → Treacher Collins as the top match, driving a *TCOF1* panel order),
  and a missing growth percentile is scored as the 50th, so "not measured" and "measured,
  normal" are indistinguishable.
- **§20** the Streamlit Variant Review tab raises `TypeError` on every classification —
  `f"{None:.3f}"` from `result.get('pathogenicity_score', 0)` where the key is present with
  value `None` — dying before the ACMG criteria block, and it never renders the route's honesty
  `notes`. So the branch's central honesty messages are invisible on the only surface a
  clinician uses, and a refusal looks like a crash.
- **§21** `api/main.py` returns `{"status": "completed", "note": "LLM unavailable; returning
  placeholder result."}` — `status` is the machine-readable field — and when Milvus is down
  *every* workflow takes that branch. The real `src/clinical_workflows.WorkflowEngine` (ten
  workflows) is referenced by nothing outside tests and has no `execute` method, so the ACMG
  fixes in this branch do not yet reach an HTTP caller.
- **§22** confidence is computed from input *length* across the agent:
  `min(0.95, 0.3 + 0.1 * len(evidence_items))` on `/query` gives a floor of **0.30 with zero
  retrieved passages**, and `_clamp(0.3 + 0.1 * len(findings))` in several workflows computes
  confidence from the number of log lines printed. One is circular:
  `_clamp(0.3 + 0.1 * pathogenic_count)` — the more variants the classifier calls Pathogenic,
  the more confident it says it is in having called them. Fixed only in the two trial workflows
  touched here.
- **§26** `src/knowledge.py` maps **HP:0010972** to *"Abnormality of the genitourinary system"*;
  the HPO term service returns HP:0010972 = *"Anemia of inadequate production"*. The correct ID
  is HP:0000119, which is already in the same dict, so the dict carries one label twice and one
  haematology term mislabelled. Served over HTTP via `/v1/diagnostic/hpo-categories`.
  (Retrieved: https://ontology.jax.org/api/hp/terms/HP%3A0010972 and `…HP%3A0000119`.)
  Separately `clinical_workflows.py` gives Gaucher `OMIM:232300` and galactosemia
  `OMIM:230800` while `agent.py` gives Gaucher `OMIM:230800` and Pompe `OMIM:232300` — the two
  files contradict each other and these IDs flow into the FHIR and Phenopacket outputs.
  omim.org returns HTTP 403 to automated fetch, so settle against
  https://www.omim.org/entry/230800, `/232300` and `/230400` by hand.
- **§25 (partly open)** disease-name matching now requires normalised equality, which fixes the
  wrong-drug and wrong-prognosis results. The **drug indication wordings themselves** are
  unverified: settle Camzyos (mavacamten) and mexiletine against DailyMed before any fix text
  restates an indication.
