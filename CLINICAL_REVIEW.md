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

### 3c. NUDT15 — SETTLED 2026-11-03, a poor metabolizer was reported as intermediate

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

## 4. DPYD activity score 0.5 — phenotype and dose

`core/agents/pharmacogenomics/src/dosing.py` labels DPYD AS 0.5 Intermediate with a 50%
reduction; `pgx_pipeline.py` calls the same score Deficient. The two disagree, so at least one
is wrong regardless of which the guideline supports.

- **Settle with:** the CPIC fluoropyrimidines/DPYD guideline dosing table (the guideline page
  moved; find the current URL from cpicpgx.org/guidelines).
- **Also:** DPYD HapB3 appears absent from the caller's allele set.

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
