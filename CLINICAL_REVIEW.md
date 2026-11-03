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

## 1. EuroSCORE II — LVEF coefficients look mis-slotted

`core/engines/cardiology/src/risk_calculators.py` (coefficient table)

```
"lvef_21_30": 0.3150652
"lvef_31_50": 0.1084443
"lvef_le_20":  0.9346919
```

`lvef_21_30` is byte-identical to `creat_cc_50_85` in the same table, which is the signature of
a copy-paste. A secondary calculator source tabulates the LV-function betas as moderate
(EF 31–50) **0.3151**, poor (EF 21–30) **0.8084**, very poor (≤20) **0.9347** — i.e. the repo's
`21_30` value appears to be the *moderate* coefficient, and 0.1084443 matches nothing found.

Every other EuroSCORE II beta in the table matched that source exactly, which makes this trio
the outlier. Monotonic ordering is preserved, so outputs stay plausible and a reference case is
the only way to catch it.

- **Settle with:** Nashef et al., *EuroSCORE II*, Eur J Cardiothorac Surg 2012, coefficient
  table; or the official documentation at euroscore.org.
- **Secondary source checked:** https://metricgate.com/docs/euroscore-ii-cardiac-surgery/
- **Also fix when settled:** `tests/test_risk_calculators.py` range-asserts this calculator
  (`0 ≤ v ≤ 100`) rather than checking a published reference case, which is why the suite is
  green. Add a reference case for ASCVD, MAGGIC and EuroSCORE II at the same time.

## 2. EuroSCORE II — age term coding

`risk_calculators.py` uses `coeff × (age − 60)`, contributing 0 below age 60. The published
coding, as far as the check could establish, is 1 for age ≤ 60 and +1 per year after, i.e.
`× max(1, age − 59)`. Effect is small but systematic.

- **Settle with:** the same Nashef 2012 table.

## 3. CYP2D6 activity score 1.0 — phenotype boundary

`core/agents/pharmacogenomics/src/pgx_pipeline.py` and `src/dosing.py` place activity score
1.0 in **Normal Metabolizer**.

The CPIC/DPWG consensus (Caudle et al. 2020) places AS 1.0 in **Intermediate Metabolizer**:
0 = poor, >0 to 1.0 = intermediate, 1.25–2.25 = normal, >2.25 = ultrarapid. The repo
contradicts itself — `pgx knowledge.py` states NM is 1.25–2.25, and the biomarker agent
classifies `*1/*4` (AS 1.0) as IM. A test pins the current value "per CPIC 2023".

This one is close to settled against secondary sources but was left alone because it changes a
dosing phenotype and the CPIC table itself was not retrieved.

- **Settle with:** the CPIC CYP2D6 genotype-to-phenotype table at cpicpgx.org.
- **Secondary sources checked:** Frontiers Pharmacol 2021 on the consensus
  (https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2021.650750/epub);
  a lab-report study showing AS 1 reclassified NM → IM
  (https://um.edu.mt/library/oar/bitstream/123456789/60804/1/FIP%20Poster%20Submission.pdf).
- **When settled:** fix `pgx_pipeline.py`, `dosing.py`, the pinning test, and the comment that
  states "IM=0.25-0.75, NM≥1.0".

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

## 7. CAR-T knowledge base — approval dates and REMS

`core/agents/cart/src/knowledge.py` dates Yescarta's ZUMA-7 indication 2021-03-05 and ZUMA-5
2024-03-15. A check of FDA's approval history found April 2021 actions and an approval letter
dated 2022-04-01, and no action on 2021-03-05. Separately, REMS requirements for approved CAR-T
products were removed on 2025-06-26, so any REMS text in this corpus is likely stale, and
Aucatzyl is absent from the product list.

- **Settle with:** FDA's Yescarta approval history
  (https://www.fda.gov/BiologicsBloodVaccines/CellularGeneTherapyProducts/ApprovedProducts/ucm581222.htm)
  and the current prescribing information.

## 8. Biomarker agent — PhenoAge inputs and advice text

- CRP is applied as ln(mg/L); the original PhenoAge formulation may use mg/dL. A unit error here
  shifts every result.
- `src/report_generator.py` advises avoiding synthetic folic acid on an MTHFR variant, which
  runs against ACMG and CDC positions on MTHFR testing and supplementation.
- "Function Health optimal" ranges appear as a reference range alongside clinical ones.

- **Settle with:** Levine et al. 2018 (PhenoAge) for units; ACMG's MTHFR statement for the
  advice text.

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

## 10. Demo VCF — HGVS internal consistency

`data/demo/demo_variants.vcf` pairs `c.463C>T` with `p.Arg155His` (a C>T at that position gives
Arg→Cys; R155H needs c.464G>A) and `c.1234G>A` with `p.Arg412His` (a G at the codon's first
base codes Gly, not Arg). A third entry, `c.476G>A` / R159H, is self-consistent. The file is
labelled synthetic and not for clinical use, so this is a correctness issue in demo data rather
than a patient-facing one.

- **Settle with:** the transcript sequence for the gene in question, then correct either the
  cDNA or the protein notation.
