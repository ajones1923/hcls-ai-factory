---
name: clinical-genomics-standards
description: >-
  The variant-interpretation standards that keep the Genomic Foundation Engine (E1) and the Rare
  Disease agent clinically correct and honestly framed. Consult whenever building, reviewing, or
  demoing anything that calls a variant pathogenic/benign, annotates from ClinVar/AlphaMissense,
  reports secondary findings, does trio/de-novo or mosaicism analysis, or scores polygenic risk —
  so classifications map to ACMG/AMP, cite their sources, and stay decision support for a
  clinician, never autonomous diagnosis.
---

# Clinical Genomics Standards — variant interpretation, honestly framed

Germline (and mosaic) variant interpretation in the factory follows the community-standard rule
sets so that every "pathogenic" or "benign" call on screen is defensible, cited, and framed as
**decision support for a qualified clinician — never a diagnosis**. This skill is the standing
reference for E1 (Genomic Foundation Engine, `localhost:5000` → VCF; variant store `:8575`) and the
Rare Disease Diagnostic agent (`:8544`).

## Standards — how a variant is classified

- **ACMG/AMP 2015 (Richards et al., *Genet Med* 2015; PMID 25741868) is the classification
  backbone.** Every reported variant lands in a 5-tier system: **Pathogenic / Likely Pathogenic /
  Variant of Uncertain Significance (VUS) / Likely Benign / Benign**, built from weighted evidence
  criteria (PVS1, PS1–4, PM1–6, PP1–5 for pathogenic; BA1, BS1–4, BP1–7 for benign). Prefer the
  quantitative **Bayesian/points refinement (Tavtigian et al., *Genet Med* 2018; PMID 30192042)**
  and current **ClinGen** gene/disease-specific rule specifications where they exist.
- **ClinVar is an annotation input, not a verdict.** Report the aggregate clinical significance
  **with its review status (gold stars)** — a 0-star single submitter is not a 3-star expert-panel
  or practice-guideline assertion. Never launder a low-confidence ClinVar entry into a factory
  "Pathogenic" call; carry the star level through to the output.
- **AlphaMissense (Cheng et al., *Science* 2023) supplies missense pathogenicity priors** (benign /
  ambiguous / pathogenic likelihood, 0–1). It is **computational evidence (ACMG PP3/BP4-strength
  at most)** — one input among many, never sufficient alone to classify. Label it as in-silico.
- **VUS is an honest answer, not a failure.** A VUS is reported as uncertain and **must not** be
  reported or acted on as pathogenic; the 2015 guideline is explicit that a VUS should not drive
  clinical action. State what evidence is missing and what would resolve it (segregation, functional
  data, phenotype fit).

## In the factory — analyses E1 and the Rare Disease agent run

- **VAF analysis to catch somatic mosaicism.** Standard germline pipelines assume a ~0.5 (het) or
  ~1.0 (hom) variant-allele fraction and can **miss low-VAF mosaic variants** — a genuine diagnostic
  trap in TSC (TSC1/TSC2 mosaicism is well documented and is a leading cause of "no mutation
  identified" results on standard testing). Inspect VAF explicitly, flag sub-heterozygous calls
  (e.g. VAF well below ~0.35), and do not let a low-VAF true positive be filtered as noise. This is a
  named differentiator vs. a vanilla germline caller.
- **Ts/Tv as a QC gate.** Whole-genome human PASS SNVs should show a transition/transversion ratio
  of **~2.0–2.1** (higher, ~3.0+, in exonic/coding regions). The variant store computes it
  (verified Ts/Tv = 2.40 on real HG002 PASS calls); a ratio far from expectation signals a
  false-positive-heavy call set — treat it as a red flag before any interpretation.
- **ACMG SF v3.2 secondary-findings panel.** The `acmg-secondary-findings` stage flags reportable
  (likely-)pathogenic variants in the **ACMG SF v3.2 medically-actionable panel (82 genes, ~90
  conditions; Miller et al., *Genet Med* 2023)**. Only report established (likely-)pathogenic
  significance in these genes; honor patient opt-out; **validate the gene list against the official
  publication before any clinical use** (panel membership changes version to version).
- **Trio / de-novo vs inherited.** Where parental samples exist, use trio analysis to establish
  inheritance — de-novo status is itself ACMG evidence (PS2/PM6) and reshapes recurrence-risk
  counseling. Confirm apparent de-novo calls against mosaicism and coverage before asserting them.
- **HPO phenotype matching.** Encode the patient phenotype in **Human Phenotype Ontology** terms and
  match against gene–phenotype associations to prioritize candidate variants — this is how a VUS in
  a phenotype-consistent gene gets prioritized for follow-up (not upgraded to pathogenic).
- **Single-box GWAS / polygenic risk.** The `gwas-association` stage runs per-variant
  logistic-regression association (statsmodels, MAF filter, covariates) and PRS against a reference
  cohort at single-box scale. Frame PRS as a **population-relative risk estimate**, ancestry-
  dependent and **not diagnostic**; state the reference cohort and its ancestry limits.
- **Evo 2 (on-demand, RunPod) for VUS context.** For a VUS with no ClinVar/functional evidence, the
  Evo 2 genomic foundation model can give **zero-shot variant-effect context**. This is
  **research-use only, elastic burst** (not local) — it informs prioritization, it does **not**
  classify a variant. Say "elastic burst," never imply it runs on the box.

## Do / Don't

**Do:** classify to ACMG/AMP 5-tier with weighted criteria; carry ClinVar review-status (stars) and
label AlphaMissense as in-silico PP3/BP4-strength evidence; inspect VAF for mosaicism; gate on
Ts/Tv; validate the ACMG SF panel against the source publication; report a VUS as uncertain and
name the missing evidence; cite the guideline/PMID with every classification.

**Don't:** call a variant pathogenic on one line of evidence (a single AlphaMissense score, a 0-star
ClinVar entry) or act on a VUS as if pathogenic; assume 0.5/1.0 VAF and silently drop mosaic calls;
present PRS, Evo 2, or any annotation as a diagnosis; skip Ts/Tv QC; ship an ACMG SF panel without
verifying its version; imply Evo 2 runs locally.

## Interpretation checklist

1. Call set passes QC (Ts/Tv ~2.0–2.1 PASS SNVs; coverage adequate)?
2. VAF inspected — mosaic / low-VAF variants preserved and flagged, not filtered as noise?
3. Each reported variant classified ACMG/AMP 5-tier with the specific criteria that fired?
4. ClinVar significance reported **with review status**; AlphaMissense labeled in-silico?
5. Inheritance (trio / de-novo) established where samples allow; phenotype encoded in HPO?
6. Secondary findings limited to ACMG SF v3.2, panel verified, opt-out honored?
7. VUS reported as uncertain with the resolving evidence named; Evo 2 flagged research-use/burst?
8. Every classification cites its guideline/source; output framed as decision support?

## Related
- `clinical-claim-honesty` — the honesty register every classification and label passes.
- `pharmacogenomics-cpic` — the other genotype→clinical-standard skill (star-alleles → CPIC).
- `oncology-mtb` — germline findings feed the molecular tumor board (Li-Fraumeni, etc.).
- `demo-foundation-alignment` — D3 (TSC, mosaicism) and D4 (germline) demo homes and honesty ledger.
