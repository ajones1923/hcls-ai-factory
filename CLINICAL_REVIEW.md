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
