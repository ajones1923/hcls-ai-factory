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
nothing else. Lung-RADS also keys on exam context: at follow-up, a NEW 4–6 mm solid nodule is
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
| CHA2DS2-VASc | Lip, Chest 2010 | all components and both age bands correct |
| HAS-BLED | Pisters, Chest 2010 | correct, including renal/liver and drugs/alcohol scoring separately and E being age >65 |
| HEART | Six 2008 / Backus 2013 | bands correct; one defect found and fixed (missing troponin scored as normal) |
| ACR TI-RADS | Tessler 2017 | all five feature tables, level mapping and FNA size thresholds correct |
| ACR Lung-RADS v2022 | ACR | solid-nodule size bands correct (baseline; see §6d) |
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
