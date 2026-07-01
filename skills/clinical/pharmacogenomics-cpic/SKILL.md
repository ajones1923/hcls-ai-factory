---
name: pharmacogenomics-cpic
description: >-
  The pharmacogenomics standards that keep the Pharmacogenomics agent (8507) clinically correct and
  honestly framed. Consult whenever building, reviewing, or demoing genotype-to-drug logic —
  star-allele calling, phenotype translation, CPIC dosing guidance, or a pre-prescription safety
  interlock (statins, thiopurines, carbamazepine, everolimus). Every output is decision support for
  a qualified prescriber — it flags and recommends, it never prescribes.
---

# Pharmacogenomics Standards (CPIC) — genotype-guided safety, not prescribing

The Pharmacogenomics agent (`localhost:8507`) turns a genotype into a **CPIC-backed dosing/selection
recommendation** and, most importantly, runs **safety interlocks** — a genotype check that fires
*before* a drug is considered, so a foreseeable adverse event is caught up front. Every result is
**decision support for a qualified prescriber; the factory never prescribes.**

## Standards — the genotype → recommendation chain

The canonical pipeline is **diplotype (star alleles) → predicted phenotype → CPIC recommendation**:

1. **Star-allele / diplotype calling** from the variant call set, using PharmVar allele definitions.
2. **Phenotype translation** to the standardized metabolizer/function terms (e.g. Poor / Intermediate
   / Normal / Rapid / Ultrarapid Metabolizer; or decreased/normal/increased function) per the
   CPIC/CPIC-PharmGKB term-standardization consensus.
3. **CPIC guideline recommendation** — dose adjustment, alternative agent, or standard dosing — with
   the CPIC strength of recommendation (Strong / Moderate / Optional) carried through.

Always cite the specific CPIC guideline (and its PharmGKB/PMID) for the gene–drug pair. CPIC
guidelines are the authoritative, peer-reviewed, freely available standard (cpicpgx.org); do not
substitute an unsourced heuristic.

## Key gene–drug pairs the agent reasons over

- **SLCO1B1 rs4149056 (c.521T>C, *5) → simvastatin.** Decreased-function transporter → higher statin
  systemic exposure → **myopathy/SAMS risk**. This is the highest-evidence statin PGx pair in the
  **CPIC 2022 statin guideline (Cooper-DeHoff et al., *Clin Pharmacol Ther* 2022)** — the flagship
  D2 example. Recommendation: lower dose or alternative statin per genotype and the desired intensity.
- **ABCG2 rs2231142 (Q141K) → rosuvastatin.** Poor-function → higher rosuvastatin exposure; CPIC 2022
  advises **rosuvastatin ≤ 20 mg/day or an alternative** in poor-function patients.
- **CYP2C9 → statins (e.g. fluvastatin) and beyond** (warfarin, NSAIDs, phenytoin). Reduced-function
  alleles raise exposure; combine with SLCO1B1 for the statin decision per CPIC.
- **CYP3A4/CYP3A5 + P-glycoprotein (ABCB1) → everolimus.** mTOR-inhibitor exposure (central to the
  TSC program) is CYP3A/P-gp driven; strong CYP3A inhibitors/inducers and transporter status shift
  levels — everolimus is **therapeutic-drug-monitoring (trough-targeted) dosed**, and PGx/interaction
  context informs the starting point and the interaction watch-list, not a fixed dose.
- **TPMT / NUDT15 → thiopurines (azathioprine, 6-MP) — a pre-treatment interlock.** CPIC
  (Relling et al., updated) recommends genotype **before starting** a thiopurine immunosuppressant;
  poor-metabolizer status → **severe, potentially fatal myelosuppression** → substantially reduced
  dose or alternative. NUDT15 is essential for correct assessment in Asian and Hispanic populations.
- **HLA-B\*15:02 → carbamazepine (and oxcarbazepine) — a hard interlock.** CPIC/FDA: screen **before**
  carbamazepine in at-risk (esp. Southeast Asian) ancestry; a positive result contraindicates the
  drug due to **Stevens-Johnson syndrome / toxic epidermal necrolysis (SCAR)**. Also note
  **HLA-B\*57:01 → abacavir** and **HLA-A\*31:01 → carbamazepine** as companion HLA interlocks.
- **Levodopa-response context (Parkinson's, D7).** Reason over reported pharmacogenetic response
  modifiers as **research-grade context**, clearly labeled — CPIC-level actionability here is limited,
  so frame it as background, not a dosing directive.

## Safety interlocks — the pattern that matters most

An interlock is a genotype gate evaluated **before a drug enters the plan**, not a note added after.
Design PGx logic so that ordering carbamazepine surfaces the HLA-B\*15:02 check, a thiopurine surfaces
TPMT/NUDT15, and a statin surfaces SLCO1B1/ABCG2 — automatically. Absence of a result is reported as
**"not assessed," never as "clear."** The interlock's job is to make a foreseeable, genotype-linked
adverse event impossible to miss; the prescriber still makes the call.

## Do / Don't

**Do:** run diplotype → phenotype → CPIC recommendation and cite the specific CPIC guideline + its
recommendation strength; fire interlocks *before* the drug is considered; carry ancestry-specific
alleles (NUDT15, HLA-B\*15:02) so at-risk populations aren't missed; frame everolimus as TDM-dosed;
report missing genotypes as "not assessed"; keep output as decision support.

**Don't:** emit a dose without the CPIC source; treat AlphaMissense/in-silico calls as a substitute
for validated star-allele definitions; imply the agent prescribes or replaces TDM; report "no
variant found" as "safe to dose"; present levodopa-response or other low-actionability PGx as a
firm directive; drop the NUDT15/HLA population-specific alleles.

## Interlock checklist (per drug considered)

1. Is there a genotype-linked CPIC pair for this drug? If yes, is the genotype present or "not
   assessed"?
2. Diplotype → phenotype done with PharmVar/CPIC-standard terms?
3. CPIC recommendation applied **with** its recommendation strength and citation?
4. HLA/SCAR and myelosuppression interlocks (HLA-B\*15:02, TPMT/NUDT15) evaluated pre-prescription?
5. Interaction/transporter context (CYP3A/P-gp for everolimus) surfaced with the TDM caveat?
6. Output labeled decision support for the prescriber, not an order?

## Related
- `clinical-claim-honesty` — the honesty register every recommendation and "not assessed" label passes.
- `clinical-genomics-standards` — where the star-allele-defining variants come from (E1 calling/QC).
- `oncology-mtb` — therapy ranking that must also respect PGx safety interlocks.
- `demo-foundation-alignment` — D2 (SLCO1B1 statin PGx) demo home and the honesty ledger.
