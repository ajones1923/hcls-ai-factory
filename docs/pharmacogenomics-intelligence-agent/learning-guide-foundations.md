# Pharmacogenomics Intelligence Agent -- Learning Guide: Foundations

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [What Is Pharmacogenomics?](#1-what-is-pharmacogenomics)
2. [Why Pharmacogenomics Matters](#2-why-pharmacogenomics-matters)
3. [The CYP450 Enzyme System](#3-the-cyp450-enzyme-system)
4. [Star Allele Nomenclature](#4-star-allele-nomenclature)
5. [Diplotypes](#5-diplotypes)
6. [Activity Scores and Phenotypes](#6-activity-scores-and-phenotypes)
7. [CPIC Guidelines](#7-cpic-guidelines)
8. [HLA Pharmacogenomics](#8-hla-pharmacogenomics)
9. [Drug Metabolism Basics](#9-drug-metabolism-basics)
10. [Clinical Scenarios](#10-clinical-scenarios)
11. [Why PGx Testing Matters](#11-why-pgx-testing-matters)
12. [How the PGx Intelligence Agent Uses These Concepts](#12-how-the-pgx-intelligence-agent-uses-these-concepts)
13. [Glossary](#13-glossary)

---

## 1. What Is Pharmacogenomics?

Pharmacogenomics (PGx) is the study of how an individual's genetic makeup affects their response to medications. The term combines "pharmacology" (the study of drugs) and "genomics" (the study of genes and their functions).

Every person inherits two copies of each gene -- one from each parent. Variations in these genes can change how the body processes drugs, leading to differences in:

- **Drug efficacy**: Whether a medication works as intended
- **Drug toxicity**: Whether a medication causes harmful side effects
- **Optimal dosing**: How much of a medication is needed for the desired effect

The goal of pharmacogenomics is to move from a "one size fits all" approach to prescribing toward **precision medicine**, where the right drug is given at the right dose to the right patient based on their genetic profile.

### Key Statistics

- Approximately **95% of people** carry at least one actionable pharmacogenomic variant.
- Adverse drug reactions cause over **100,000 deaths per year** in the United States alone.
- Pharmacogenomic variation contributes to **30-60% of treatment failures** for common drug classes.
- The global economic burden of adverse drug reactions exceeds **$528 billion per year**.
- ADRs account for approximately **7% of all hospital admissions**.
- The FDA lists pharmacogenomic biomarkers in the labeling of **350+ drugs**.
- CPIC has published evidence-based guidelines for **100+ gene-drug pairs**.

### The Three Pillars of Pharmacogenomics

Pharmacogenomics is built on three interconnected pillars:

```
PILLAR 1: PHARMACOKINETICS (PK)
"What the body does to the drug"
  - How fast the drug is metabolized (CYP enzymes)
  - How the drug is transported across membranes (SLCO1B1, ABCB1)
  - How the drug is conjugated for excretion (UGT1A1, TPMT, DPYD)

PILLAR 2: PHARMACODYNAMICS (PD)
"What the drug does to the body"
  - Drug target sensitivity (VKORC1 for warfarin)
  - Receptor binding affinity (OPRM1 for opioids)
  - Downstream signaling (COMT for catecholamine metabolism)

PILLAR 3: IMMUNOGENETICS
"How the immune system reacts to the drug"
  - HLA-mediated hypersensitivity (HLA-B*15:02 + carbamazepine)
  - Immune-mediated drug reactions (SJS/TEN, DRESS)
  - Idiosyncratic hepatotoxicity (HLA-B*57:01 + flucloxacillin)
```

Most current CPIC guidelines focus on Pillar 1 (pharmacokinetic PGx), but all three pillars are represented in the PGx Intelligence Agent's knowledge base.

### A Brief History of Pharmacogenomics

The roots of pharmacogenomics extend back to the 1950s, when clinicians first observed that some patients experienced unexpected reactions to standard drug doses:

- **1956**: Primaquine sensitivity linked to glucose-6-phosphate dehydrogenase (G6PD) deficiency, one of the earliest pharmacogenetic associations.
- **1957**: Succinylcholine-induced prolonged apnea linked to butyrylcholinesterase variants.
- **1977**: Debrisoquine/sparteine hydroxylation polymorphism discovered (later mapped to CYP2D6).
- **1988**: N-acetyltransferase (NAT2) "slow acetylator" phenotype linked to isoniazid toxicity.
- **2000**: Human Genome Project completed, enabling genome-wide pharmacogenomic studies.
- **2004**: FDA approves AmpliChip CYP450, the first pharmacogenomic test.
- **2005**: CPIC established to create evidence-based PGx guidelines.
- **2010**: PharmVar consortium formed to standardize star allele nomenclature.
- **2023**: PREPARE trial (Lancet) demonstrates 30% ADR reduction with pre-emptive PGx testing.

---

## 2. Why Pharmacogenomics Matters

### Real-World Examples

**Codeine and CYP2D6**: Codeine is a prodrug that must be converted to morphine by the CYP2D6 enzyme to produce pain relief. Patients who are CYP2D6 poor metabolizers get no benefit from codeine (they cannot convert it to morphine). Conversely, ultra-rapid metabolizers convert codeine too quickly, producing dangerously high morphine levels that can cause respiratory depression -- deaths have been reported in children. In 2013, the FDA added a contraindication for codeine use in children after tonsillectomy following multiple pediatric deaths linked to CYP2D6 ultra-rapid metabolism.

**Clopidogrel and CYP2C19**: Clopidogrel (Plavix) is an antiplatelet prodrug activated by CYP2C19. Poor metabolizers cannot activate clopidogrel effectively, leaving them at increased risk of stent thrombosis and cardiovascular events after coronary stent placement. The FDA label includes a boxed warning about CYP2C19 poor metabolizer status. Studies show that CYP2C19 poor metabolizers have a 3.4-fold increased risk of stent thrombosis compared to normal metabolizers (Mega et al., NEJM 2009).

**Carbamazepine and HLA-B*15:02**: Patients carrying the HLA-B*15:02 allele (prevalent in Southeast Asian populations) have a dramatically increased risk of Stevens-Johnson syndrome (SJS) and toxic epidermal necrolysis (TEN) -- life-threatening skin reactions -- when treated with carbamazepine. The FDA requires HLA-B*15:02 testing before prescribing carbamazepine. The relative risk of SJS/TEN in HLA-B*15:02 carriers is estimated at 100-fold compared to non-carriers.

**Abacavir and HLA-B*57:01**: Abacavir, an HIV antiretroviral, causes hypersensitivity syndrome in approximately 5-8% of patients carrying HLA-B*57:01. Mandatory HLA-B*57:01 testing before abacavir prescription has virtually eliminated this reaction. This is widely considered the most successful example of clinical pharmacogenomics implementation -- the PREDICT-1 trial showed that HLA-B*57:01 screening has a 100% negative predictive value for abacavir hypersensitivity.

---

## 3. The CYP450 Enzyme System

The cytochrome P450 (CYP450) enzymes are a superfamily of proteins primarily found in the liver. They are responsible for metabolizing (breaking down) the majority of medications. Understanding the CYP450 system is fundamental to pharmacogenomics.

### Phase I Metabolism: The CYP450 Oxidative System

CYP450 enzymes catalyze **Phase I oxidative reactions**, which modify drug molecules through:

- **Hydroxylation**: Adding an -OH group (most common CYP reaction)
- **Dealkylation**: Removing alkyl groups (N-dealkylation, O-dealkylation)
- **Oxidation**: Adding oxygen to nitrogen or sulfur atoms
- **Epoxidation**: Creating epoxide intermediates
- **Dehalogenation**: Removing halogen atoms

These reactions typically make drug molecules more polar (water-soluble), facilitating their excretion by the kidneys.

### The Seven Key CYP450 Enzymes in Pharmacogenomics

| Enzyme | Chromosome | Drugs Metabolized | % of All Drugs | PGx Relevance | Key Variants |
|--------|-----------|------------------|----------------|---------------|-------------|
| CYP2D6 | 22q13.2 | Opioids, antidepressants, tamoxifen, beta-blockers | ~25% | Most polymorphic CYP; >150 alleles; gene deletions/duplications | *3, *4, *5, *6, *10, *17, *41 |
| CYP2C19 | 10q23.33 | Clopidogrel, PPIs, antidepressants, voriconazole | ~10% | Loss-of-function alleles affect antiplatelet therapy | *2, *3, *17 |
| CYP2C9 | 10q23.33 | Warfarin, phenytoin, NSAIDs, sulfonylureas | ~15% | Critical for warfarin dosing | *2, *3, *5, *6, *8, *11 |
| CYP3A4 | 7q22.1 | Immunosuppressants, statins, calcium channel blockers | ~30% | Most abundant hepatic CYP; few high-impact variants | *22 |
| CYP3A5 | 7q22.1 | Tacrolimus | Subset of CYP3A4 | Expresser vs non-expresser affects tacrolimus dose | *3, *6, *7 |
| CYP1A2 | 15q24.1 | Caffeine, theophylline, clozapine | ~5% | Inducible by smoking, cruciferous vegetables | *1F, *1C |
| CYP2B6 | 19q13.2 | Efavirenz, methadone, cyclophosphamide | ~5% | Important in HIV treatment | *6, *18 |

### How CYP Enzymes Work

```
                    CYP450 Enzyme
                         |
                         v
Drug (parent) ---[oxidation/reduction]--> Metabolite
                                              |
                                              v
                                    More polar (water-soluble)
                                              |
                                              v
                                    Excretion (kidneys, bile)
```

### CYP Enzyme Location and Expression

CYP450 enzymes are primarily expressed in the liver (hepatocytes), but are also found in:

- **Small intestine**: CYP3A4 in enterocytes contributes to first-pass metabolism of orally administered drugs
- **Kidneys**: CYP enzymes contribute to renal drug metabolism
- **Lungs**: CYP1A1, CYP2B6 expressed in bronchial epithelium
- **Brain**: CYP2D6 expression in the brain may contribute to local neurotransmitter metabolism

The abundance of each CYP enzyme in the liver varies:

```
CYP3A4/5  ||||||||||||||||||||||||||||||||  30% of total hepatic CYP
CYP2C9    ||||||||||||||||||                 20%
CYP2D6    ||||                                5% (but metabolizes 25% of drugs)
CYP2C19   ||||                                5%
CYP1A2    ||||||||||||                        13%
CYP2E1    |||||||                              7%
CYP2B6    ||||                                 5%
Other CYP ||||||||||||||                       15%
```

CYP2D6 is notably the enzyme that metabolizes the most drugs relative to its expression level -- it represents only ~5% of hepatic CYP content but metabolizes approximately 25% of all drugs.

### Active Drug vs. Prodrug: A Critical Distinction

```
ACTIVE DRUG (e.g., omeprazole):
  Active parent compound --> CYP2C19 --> Inactive metabolite --> Excretion
  PM impact: Drug ACCUMULATES (increased effect/toxicity)

PRODRUG (e.g., codeine):
  Inactive parent compound --> CYP2D6 --> Active metabolite (morphine) --> Therapeutic effect
  PM impact: Drug NOT ACTIVATED (treatment failure)
```

For most drugs, metabolism by CYP enzymes **inactivates** the drug. However, for **prodrugs** (codeine, clopidogrel, tamoxifen), CYP metabolism **activates** the drug. This distinction is critical for understanding PGx implications:

- **Active drugs + Poor Metabolizer**: Drug accumulates in the body (increased toxicity risk, prolonged effect)
- **Active drugs + Ultra-Rapid Metabolizer**: Drug cleared too quickly (sub-therapeutic levels, treatment failure)
- **Prodrugs + Poor Metabolizer**: Drug cannot be activated (no therapeutic effect)
- **Prodrugs + Ultra-Rapid Metabolizer**: Drug activated too quickly (excessive active metabolite, toxicity)

### Phase II Metabolism: Conjugation Enzymes

After Phase I modification, many drugs undergo **Phase II conjugation reactions** that attach water-soluble groups for excretion. Key Phase II enzymes in PGx:

| Enzyme | Reaction Type | PGx Drugs | Clinical Impact |
|--------|--------------|----------|----------------|
| UGT1A1 | Glucuronidation | Irinotecan | *28 homozygotes: severe neutropenia |
| NAT2 | Acetylation | Isoniazid, hydralazine, procainamide | Slow acetylators: increased toxicity |
| TPMT | Methylation | Azathioprine, mercaptopurine, thioguanine | Deficiency: life-threatening myelosuppression |
| DPYD | Pyrimidine catabolism | Fluorouracil, capecitabine | Deficiency: severe/fatal fluoropyrimidine toxicity |

### Drug Transporters

Transporter proteins move drugs across cell membranes, affecting absorption, distribution, and excretion:

| Transporter | Gene | PGx Impact |
|-----------|------|-----------|
| OATP1B1 | SLCO1B1 | Statin uptake into liver; rs4149056 (Val174Ala) = myopathy risk with simvastatin |
| P-glycoprotein | ABCB1 | Efflux pump; affects absorption of many drugs at the gut, blood-brain barrier, and kidneys |

---

## 4. Star Allele Nomenclature

Pharmacogenes use a unique naming system called **star allele nomenclature** to describe genetic variants. Star alleles are designated with an asterisk (*) followed by a number.

### How Star Alleles Work

- **\*1**: The "reference" or "wild-type" allele, representing normal enzyme function. By convention, *1 is the allele present in the human reference genome.
- **\*2, \*3, \*4, ...**: Variant alleles, each defined by a specific set of DNA changes (SNPs, insertions, deletions, or structural rearrangements).

### Star Allele Assignment

Each star allele is defined by one or more **defining variants** (specific DNA changes). A patient's genotype at these positions determines which star allele(s) they carry:

```
Step 1: Sequence the pharmacogene region
Step 2: Identify all variants (SNPs, indels) compared to reference genome
Step 3: Match the variant combination to a known star allele definition
Step 4: Assign star alleles (one per chromosome = two per patient)
```

### Examples: CYP2D6 Star Alleles

| Star Allele | Function | Defining Variant(s) | rsID(s) | Mechanism | Clinical Impact |
|-------------|----------|---------------------|---------|-----------|----------------|
| *1 | Normal function | Reference sequence | N/A | Fully functional enzyme | Standard metabolism |
| *2 | Normal function | Multiple SNPs | rs16947 | Amino acid changes, functional | Standard metabolism |
| *3 | No function | 2549del (frameshift) | rs35742686 | Truncated, non-functional protein | Cannot metabolize CYP2D6 substrates |
| *4 | No function | 1846G>A (splice site) | rs3892097 | Aberrant mRNA splicing | Most common no-function allele in Europeans (~20-25%) |
| *5 | No function | Whole gene deletion | N/A | Complete absence of CYP2D6 | No enzyme produced |
| *6 | No function | 1707delT (frameshift) | rs5030655 | Truncated protein | Rare no-function allele |
| *9 | Decreased function | 2615-2617delAAG | rs5030656 | In-frame deletion of Lys281 | Reduced substrate binding |
| *10 | Decreased function | 100C>T (Pro34Ser) + 4180G>C (Ser486Thr) | rs1065852, rs1135840 | Unstable protein, reduced activity | Most common decreased allele in East Asians (~35-45%) |
| *17 | Decreased function | 1023C>T (Thr107Ile) + 2850C>T (Arg296Cys) | rs28371706, rs16947 | Altered substrate specificity | Most common decreased allele in Africans (~20-30%) |
| *29 | Decreased function | Multiple variants | Multiple | Reduced activity | Common in African populations |
| *41 | Decreased function | 2988G>A (splice variant) | rs28371725 | Reduced mRNA expression | Moderate decrease in activity |
| *1xN | Increased function | Gene duplication (2-13 copies) | N/A | Multiple gene copies | Ultra-rapid metabolism |
| *2xN | Increased function | Gene duplication | N/A | Multiple functional copies | Ultra-rapid metabolism |

### Examples: CYP2C19 Star Alleles

| Star Allele | Function | Defining Variant | rsID | Clinical Impact |
|-------------|----------|-----------------|------|----------------|
| *1 | Normal function | Reference | N/A | Standard metabolism |
| *2 | No function | 681G>A (splice site) | rs4244285 | Most common no-function allele globally |
| *3 | No function | 636G>A (premature stop) | rs4986893 | Common in East Asian populations |
| *17 | Increased function | -806C>T (promoter) | rs12248560 | Increased CYP2C19 expression; rapid metabolism |

### Examples: CYP2C9 Star Alleles

| Star Allele | Function | Defining Variant | rsID | Clinical Impact |
|-------------|----------|-----------------|------|----------------|
| *1 | Normal function | Reference | N/A | Standard warfarin metabolism |
| *2 | Decreased function | Arg144Cys | rs1799853 | ~30% reduced activity; lower warfarin dose |
| *3 | Decreased function | Ile359Leu | rs1057910 | ~80% reduced activity; significantly lower warfarin dose |
| *5 | No function | Asp360Glu | rs28371686 | No-function; very low warfarin dose |
| *6 | No function | Frame deletion | rs9332131 | No-function allele |
| *8 | Decreased function | Arg150His | rs7900194 | Common in African populations |
| *11 | Decreased function | Arg335Trp | rs28371685 | Moderately decreased activity |

### Star Alleles Are Gene-Specific

Star allele numbering is independent across genes. CYP2D6\*4 and CYP2C19\*4 are completely different variants with different clinical implications. Star allele definitions are curated by the **Pharmacogene Variation Consortium (PharmVar)** at pharmvar.org.

PharmVar maintains the authoritative catalog of star allele definitions. As of 2026, PharmVar catalogs:
- **CYP2D6**: >150 star alleles (most polymorphic CYP enzyme)
- **CYP2C19**: >35 star alleles
- **CYP2C9**: >60 star alleles
- **CYP2B6**: >40 star alleles

### Reading a Star Allele Report

When a patient receives PGx testing results, the report typically contains:

```
+--------------------------------------------------+
|  PHARMACOGENOMIC TEST RESULT                      |
|                                                   |
|  Gene:      CYP2D6                               |
|  Genotype:  *1/*4                                |
|  Phenotype: Intermediate Metabolizer             |
|                                                   |
|  Allele 1:  *1 (Normal function)                 |
|  Allele 2:  *4 (No function)                     |
|             Defining variant: rs3892097 (G>A)    |
|             Splice site defect                    |
|                                                   |
|  Activity Score: 1.0 (1.0 + 0.0)                |
|                                                   |
|  Clinical Significance:                          |
|  Reduced CYP2D6 metabolism. Consider dose        |
|  adjustment for CYP2D6 substrates. AVOID         |
|  codeine (prodrug - insufficient activation).    |
+--------------------------------------------------+
```

### The Difference Between rsIDs and Star Alleles

rsIDs (Reference SNP IDs, from dbSNP) refer to individual DNA positions. Star alleles refer to combinations of DNA positions that define a named variant. The relationship is:

```
rsIDs (individual DNA positions):
  rs3892097  = single nucleotide change at CYP2D6 position 1846
  rs16947    = single nucleotide change at CYP2D6 position 2850

Star alleles (combinations of rsIDs):
  CYP2D6*4  = defined by rs3892097 (G>A) + several other variants
  CYP2D6*2  = defined by rs16947 (G>A) + several other variants

One star allele may require multiple rsIDs to define it.
One rsID may appear in multiple star allele definitions.
```

---

## 5. Diplotypes

Because humans inherit two copies of each gene (one from each parent), a patient's pharmacogenomic genotype is described as a **diplotype** -- the combination of their two star alleles.

### Diplotype Notation

Diplotypes are written as two star alleles separated by a forward slash:

```
Gene     Diplotype     Meaning
CYP2D6   *1/*4        One normal allele + one no-function allele
CYP2C19  *1/*2        One normal allele + one no-function allele
CYP2C9   *1/*1        Two normal alleles (reference/reference)
CYP2D6   *4/*4        Two no-function alleles (homozygous no-function)
CYP2D6   *1/*1xN      One normal + one duplicated normal (ultra-rapid)
CYP2D6   *10/*10      Two decreased alleles (common in East Asians)
CYP2D6   *4/*41       No-function + decreased-function
CYP2C19  *2/*3        Two different no-function alleles (compound het)
```

### Common Diplotype Patterns

```
Homozygous reference:     *1/*1    (two normal copies)
Heterozygous carrier:     *1/*4    (one normal, one variant)
Homozygous variant:       *4/*4    (two variant copies)
Compound heterozygous:    *4/*5    (two different variant copies)
Duplicated allele:        *1/*2xN  (normal + duplicated variant)
```

### From VCF to Diplotype

The Pharmacogenomics Intelligence Agent's `StarAlleleCaller` class performs this translation:

```
VCF File (variant calls from genomic sequencing)
     |
     v
Identify variants in pharmacogene regions (25 genes)
     |
     v
Match variant combinations to PharmVar star allele definitions
     |
     v
Assemble diplotype (two alleles per gene)
     |
     v
CYP2D6: rs3892097 G/A --> heterozygous for *4 defining variant
        No other defining variants detected
        Diplotype: *1/*4

CYP2C19: rs4244285 G/A --> heterozygous for *2 defining variant
         rs12248560 C/C --> reference (no *17)
         Diplotype: *1/*2
```

---

## 6. Activity Scores and Phenotypes

### Activity Scores

Each star allele is assigned a numerical **activity score** reflecting its enzyme function:

| Function Status | Activity Score | CYP2D6 Examples | CYP2C19 Examples |
|----------------|---------------|-----------------|------------------|
| No function | 0.0 | *3, *4, *5, *6 | *2, *3 |
| Decreased function | 0.5 | *9, *10, *17, *29, *41 | N/A |
| Normal function | 1.0 | *1, *2 | *1 |
| Increased function | 1.5 | N/A | *17 |
| Increased function (duplication) | 2.0+ | *1xN, *2xN | N/A |

### Diplotype Activity Score Calculation

The diplotype activity score is the **sum** of both allele scores:

```
EXAMPLE 1: CYP2D6 *1/*4
  Allele 1: *1  = 1.0 (normal function)
  Allele 2: *4  = 0.0 (no function)
  Diplotype AS  = 1.0 + 0.0 = 1.0

EXAMPLE 2: CYP2D6 *4/*4
  Allele 1: *4  = 0.0 (no function)
  Allele 2: *4  = 0.0 (no function)
  Diplotype AS  = 0.0 + 0.0 = 0.0

EXAMPLE 3: CYP2D6 *1/*41
  Allele 1: *1  = 1.0 (normal function)
  Allele 2: *41 = 0.5 (decreased function)
  Diplotype AS  = 1.0 + 0.5 = 1.5

EXAMPLE 4: CYP2D6 *1/*1xN (gene duplication)
  Allele 1: *1  = 1.0 (normal function)
  Allele 2: *1xN = 2.0+ (duplicated normal)
  Diplotype AS  = 1.0 + 2.0 = 3.0

EXAMPLE 5: CYP2D6 *10/*10
  Allele 1: *10 = 0.5 (decreased function)
  Allele 2: *10 = 0.5 (decreased function)
  Diplotype AS  = 0.5 + 0.5 = 1.0

EXAMPLE 6: CYP2D6 *4/*41
  Allele 1: *4  = 0.0 (no function)
  Allele 2: *41 = 0.5 (decreased function)
  Diplotype AS  = 0.0 + 0.5 = 0.5
```

### Metabolizer Phenotypes

The activity score maps to a standardized **metabolizer phenotype** using CPIC-defined thresholds:

| Phenotype | Abbreviation | CYP2D6 Activity Score | Example Diplotypes | Clinical Meaning |
|-----------|-------------|----------------------|-------------------|-----------------|
| Ultra-Rapid Metabolizer | UM | > 2.25 | *1/*2xN, *2/*2xN | Metabolizes drugs too fast; prodrugs = toxicity risk, regular drugs = sub-therapeutic |
| Rapid Metabolizer | RM | 2.0 - 2.25 | *1/*1xN (with 2 copies) | Slightly faster than normal; limited clinical significance for most drugs |
| Normal Metabolizer | NM | 1.25 - 2.0 | *1/*1, *1/*2, *1/*41 | Standard drug response expected |
| Intermediate Metabolizer | IM | 0.25 - 1.0 | *1/*4, *4/*41, *10/*10 | Reduced metabolism; may need dose adjustment |
| Poor Metabolizer | PM | 0.0 | *4/*4, *4/*5, *3/*4 | Cannot metabolize via this enzyme; prodrugs fail, regular drugs accumulate |

Note: The exact activity score boundaries for CYP2D6 phenotype assignment follow the 2020 CPIC standardization by Caudle et al. (J Clin Pharmacol 2020). Boundaries differ slightly between genes.

### Understanding the Phenotype Spectrum

```
PHENOTYPE SPECTRUM (CYP2D6)

Ultra-Rapid (UM)     Rapid (RM)     Normal (NM)     Intermediate (IM)     Poor (PM)
|<--- Too Fast --->|<-- Faster -->|<-- Standard -->|<--- Slower --->|<--- No Activity --->|
|                   |              |                |                |                      |
| Prodrug toxicity  | Monitor      | Expected       | Dose adjust    | Prodrug failure      |
| Active drug:      |              | response       | Active drug:   | Active drug:         |
|  sub-therapeutic  |              |                |  accumulation  |  severe accumulation  |
|                   |              |                |                |                      |
|  AS > 2.25        | AS 2.0-2.25  | AS 1.25-2.0    | AS 0.25-1.0    | AS = 0.0             |
```

### CYP2C19 Activity Score System

CYP2C19 uses a slightly different activity score system due to the presence of *17 (an increased-function allele):

| Allele | Activity Score | Function |
|--------|---------------|----------|
| *1 | 1.0 | Normal |
| *2 | 0.0 | No function |
| *3 | 0.0 | No function |
| *17 | 1.5 | Increased function |

| Diplotype | AS Calculation | Total AS | Phenotype |
|-----------|---------------|----------|-----------|
| *1/*1 | 1.0 + 1.0 | 2.0 | NM |
| *1/*17 | 1.0 + 1.5 | 2.5 | RM |
| *17/*17 | 1.5 + 1.5 | 3.0 | UM |
| *1/*2 | 1.0 + 0.0 | 1.0 | IM |
| *2/*17 | 0.0 + 1.5 | 1.5 | IM |
| *2/*2 | 0.0 + 0.0 | 0.0 | PM |
| *2/*3 | 0.0 + 0.0 | 0.0 | PM |

### Complete CYP2D6 Diplotype-to-Phenotype Walkthrough

| Diplotype | AS Calculation | Total AS | Phenotype |
|-----------|---------------|----------|-----------|
| *1/*1 | 1.0 + 1.0 | 2.0 | NM |
| *1/*2 | 1.0 + 1.0 | 2.0 | NM |
| *1/*4 | 1.0 + 0.0 | 1.0 | IM |
| *1/*5 | 1.0 + 0.0 | 1.0 | IM |
| *1/*9 | 1.0 + 0.5 | 1.5 | NM |
| *1/*10 | 1.0 + 0.5 | 1.5 | NM |
| *1/*17 | 1.0 + 0.5 | 1.5 | NM |
| *1/*41 | 1.0 + 0.5 | 1.5 | NM |
| *2/*2 | 1.0 + 1.0 | 2.0 | NM |
| *4/*4 | 0.0 + 0.0 | 0.0 | PM |
| *4/*5 | 0.0 + 0.0 | 0.0 | PM |
| *4/*10 | 0.0 + 0.5 | 0.5 | IM |
| *4/*41 | 0.0 + 0.5 | 0.5 | IM |
| *10/*10 | 0.5 + 0.5 | 1.0 | IM |
| *10/*17 | 0.5 + 0.5 | 1.0 | IM |
| *10/*41 | 0.5 + 0.5 | 1.0 | IM |
| *1/*1xN (3 copies) | 1.0 + 3.0 | 4.0 | UM |
| *1/*2xN (3 copies) | 1.0 + 3.0 | 4.0 | UM |

### Population Distribution (CYP2D6)

| Phenotype | European | East Asian | African | South Asian |
|-----------|---------|-----------|---------|------------|
| UM | 1-2% | <1% | 2-5% | 1-2% |
| NM | 70-80% | 50-60% | 60-70% | 65-75% |
| IM | 10-15% | 30-40% | 20-30% | 15-25% |
| PM | 5-10% | 1-2% | 2-5% | 1-3% |

---

## 7. CPIC Guidelines

The **Clinical Pharmacogenetics Implementation Consortium (CPIC)** develops evidence-based, peer-reviewed guidelines for translating genetic test results into prescribing decisions.

### How CPIC Guidelines Work

Each CPIC guideline covers a specific **gene-drug pair** and provides:

1. **Genotype-to-phenotype translation**: How to convert star allele diplotypes to metabolizer phenotypes
2. **Prescribing recommendations**: What to do for each phenotype (standard dose, adjusted dose, alternative drug, avoid)
3. **Evidence level classification**: How strong the evidence is
4. **Population considerations**: Allele frequency data and ancestry-specific nuances

### CPIC Evidence Levels

| Level | Meaning | Strength of Evidence | Clinical Action | Example Gene-Drug Pairs |
|-------|---------|---------------------|----------------|------------------------|
| A | Strong evidence, actionable | Based on replicated studies with consistent results; meta-analyses and clinical trials support | Change prescribing recommended; strong gene-drug association | CYP2D6/codeine, CYP2C19/clopidogrel, HLA-B*57:01/abacavir |
| A/B | Strong to moderate evidence | Strong evidence for some phenotypes, moderate for others | Change prescribing may be recommended | CYP2D6/tramadol, CYP2C19/escitalopram |
| B | Moderate evidence | Based on significant findings in well-designed studies, but lacking replication or having some inconsistency | Consider alternative if available | CYP2D6/fluvoxamine, CYP2C19/sertraline |
| C | Weak evidence | Based on case reports, pharmacokinetic studies, or limited clinical data | Monitor, no prescribing change required | Some newer gene-drug associations |
| D | Insufficient evidence | Data are conflicting, of poor quality, or absent | No clinical action needed; informational only | Emerging or poorly studied pairs |

### Example: CPIC Guideline for CYP2D6 and Codeine (Level A)

| CYP2D6 Phenotype | Activity Score | CPIC Recommendation | Rationale |
|-------------------|---------------|---------------------|-----------|
| Ultra-Rapid Metabolizer | > 2.25 | **Avoid codeine**; use non-CYP2D6 analgesic | Excessive morphine formation; respiratory depression and death risk |
| Rapid Metabolizer | 2.0-2.25 | Use codeine with caution; monitor for toxicity | Slightly increased morphine formation |
| Normal Metabolizer | 1.25-2.0 | Use codeine per standard dosing | Standard metabolism expected |
| Intermediate Metabolizer | 0.25-1.0 | Use codeine with caution; monitor for reduced efficacy | Reduced morphine formation; may have inadequate pain relief |
| Poor Metabolizer | 0.0 | **Avoid codeine**; use non-CYP2D6 analgesic | Cannot form morphine; no analgesic effect from codeine |

### Example: CPIC Guideline for CYP2C19 and Clopidogrel (Level A)

| CYP2C19 Phenotype | Recommendation | Rationale |
|-------------------|---------------|-----------|
| Ultra-Rapid/Rapid Metabolizer | Standard clopidogrel dosing | Enhanced activation; adequate antiplatelet effect |
| Normal Metabolizer | Standard clopidogrel dosing | Normal activation |
| Intermediate Metabolizer | Consider alternative (prasugrel, ticagrelor) if high cardiovascular risk | Reduced activation; moderately increased thrombosis risk |
| Poor Metabolizer | **Use alternative**: prasugrel or ticagrelor | Severely reduced activation; high thrombosis risk; FDA Boxed Warning |

### Example: CPIC Guideline for CYP2C9/VKORC1 and Warfarin (Level A)

| CYP2C9 + VKORC1 Status | Dose Adjustment |
|------------------------|----------------|
| CYP2C9 *1/*1 + VKORC1 GG | Standard dose (~35 mg/week) |
| CYP2C9 *1/*2 + VKORC1 AG | Reduce dose 20-30% |
| CYP2C9 *1/*3 + VKORC1 AG | Reduce dose 30-40% |
| CYP2C9 *2/*3 + VKORC1 AA | Reduce dose 60-70% |
| CYP2C9 *3/*3 + VKORC1 AA | Reduce dose 70-80% (very low dose) |

### Other Guideline Bodies

| Organization | Abbreviation | Region | Focus |
|-------------|-------------|--------|-------|
| Clinical Pharmacogenetics Implementation Consortium | CPIC | USA/Global | Evidence-based prescribing guidelines |
| Dutch Pharmacogenetics Working Group | DPWG | Netherlands/Europe | Therapeutic dose recommendations |
| Food and Drug Administration | FDA | USA | Drug labeling (pharmacogenomic biomarkers) |
| Canadian Pharmacogenomics Network for Drug Safety | CPNDS | Canada | ADR prevention guidelines |

---

## 8. HLA Pharmacogenomics

The **Human Leukocyte Antigen (HLA)** system is part of the immune system. Certain HLA alleles can trigger severe immune-mediated drug hypersensitivity reactions.

### The MHC and HLA System

HLA genes are located on chromosome 6p21 within the **Major Histocompatibility Complex (MHC)**. They encode cell-surface proteins that present peptide fragments to T cells, enabling the immune system to distinguish "self" from "foreign."

```
Chromosome 6p21: The MHC Region
+----------------------------------------------------+
|  Class I          |  Class III      |  Class II      |
|  HLA-A            |  Complement     |  HLA-DR        |
|  HLA-B            |  TNF genes      |  HLA-DQ        |
|  HLA-C            |  Heat shock     |  HLA-DP        |
+----------------------------------------------------+
        |                                      |
        v                                      v
  Present peptides to               Present peptides to
  CD8+ T cells (cytotoxic)         CD4+ T cells (helper)
```

In drug hypersensitivity, certain drugs (or their metabolites) bind to specific HLA proteins and create a "neo-antigen" that is recognized as foreign by T cells. This triggers an immune response that can range from mild skin rash to fatal toxic epidermal necrolysis.

### Mechanism of HLA-Mediated Drug Reactions

Three models explain how drugs interact with HLA molecules:

1. **Hapten model**: Drug metabolite binds covalently to self-peptides, altering their presentation
2. **Pharmacological interaction (p-i) model**: Drug binds directly to the HLA molecule non-covalently, mimicking a foreign antigen
3. **Altered peptide repertoire model**: Drug binds inside the HLA peptide-binding groove, changing which self-peptides are presented

### The 12 Key HLA-Drug Associations

| HLA Allele | Drug | Reaction | Severity | Mortality | Testing Status |
|-----------|------|----------|----------|-----------|---------------|
| HLA-B*57:01 | Abacavir | Hypersensitivity syndrome | Severe | 1-5% | FDA required |
| HLA-B*15:02 | Carbamazepine | SJS/TEN | Fatal/Severe | 10-30% | FDA recommended |
| HLA-B*15:02 | Phenytoin | SJS/TEN | Fatal/Severe | 10-30% | FDA recommended |
| HLA-B*15:02 | Oxcarbazepine | SJS/TEN | Severe | 5-15% | FDA recommended |
| HLA-A*31:01 | Carbamazepine | DRESS/SJS | Severe | 5-10% | CPIC recommended |
| HLA-B*58:01 | Allopurinol | SJS/TEN | Fatal/Severe | 10-25% | CPIC recommended |
| HLA-B*15:11 | Carbamazepine | SJS/TEN | Severe | 5-15% | CPIC Level B |
| HLA-A*02:01 | Carbamazepine | SJS/TEN | Severe | 5-15% | Emerging evidence |
| HLA-B*57:01 | Flucloxacillin | DILI | Moderate | 1-2% | Emerging evidence |
| HLA-DRB1*07:01 | Lapatinib | Hepatotoxicity | Moderate | <1% | Emerging evidence |
| HLA-A*33:03 | Ticlopidine | Hepatotoxicity | Moderate | <1% | Emerging evidence |
| HLA-B*13:01 | Dapsone | Hypersensitivity | Severe | 2-5% | Emerging evidence |

### Understanding the Reactions

- **SJS (Stevens-Johnson Syndrome)**: Severe skin reaction with blistering and skin detachment affecting < 10% of body surface area. Mortality 1-5%. Onset typically 1-3 weeks after drug initiation.
- **TEN (Toxic Epidermal Necrolysis)**: Extreme form of SJS with > 30% skin detachment. Mortality 25-35%. Requires ICU-level care.
- **SJS/TEN overlap**: Skin detachment of 10-30% of body surface area. Mortality 10-15%.
- **DRESS (Drug Reaction with Eosinophilia and Systemic Symptoms)**: Multi-organ hypersensitivity with skin rash, fever, lymphadenopathy, and internal organ involvement (liver, kidney, heart). Mortality 5-10%. Onset typically 2-8 weeks after drug initiation.
- **DILI (Drug-Induced Liver Injury)**: Immune-mediated hepatotoxicity.

### Population Prevalence

HLA allele frequencies vary dramatically across populations, which has direct implications for pre-treatment testing recommendations:

| HLA Allele | European | Southeast Asian | African | East Asian | South Asian |
|-----------|---------|-----------------|---------|-----------|------------|
| HLA-B*15:02 | < 1% | 2-15% | < 1% | 2-6% | 2-4% |
| HLA-B*57:01 | 6-8% | 1-2% | < 1% | 1-2% | 2-5% |
| HLA-B*58:01 | 1-2% | 6-8% | 3-6% | 4-8% | 4-6% |
| HLA-A*31:01 | 2-5% | 1-3% | 1-2% | 2-4% | 2-4% |

### The Cost of Missing an HLA Screen

The clinical and economic consequences of missing an HLA-mediated drug reaction are severe:

| Outcome | Cost | Impact |
|---------|------|--------|
| SJS requiring hospitalization | $20,000-100,000+ | ICU care, skin grafting, 2-4 week stay |
| TEN requiring burn unit care | $100,000-500,000+ | ICU, ventilatory support, skin grafting, 30-50 day stay |
| DRESS with organ involvement | $30,000-80,000 | Liver/kidney monitoring, immunosuppression |
| Fatal SJS/TEN | $200,000+ | Loss of life; litigation |
| Abacavir hypersensitivity (hospital) | $10,000-30,000 | Drug discontinuation, supportive care |

The cost of a single HLA test ($50-150) is negligible compared to the cost of treating a severe immune-mediated drug reaction. Pre-prescription HLA screening is one of the most cost-effective PGx interventions available.

### How HLA Reactions Differ from Dose-Dependent ADRs

```
DOSE-DEPENDENT ADR (Type A):                    HLA-MEDIATED ADR (Type B):
  Predictable from drug pharmacology              Unpredictable from pharmacology
  Occurs at higher doses                           Occurs at any dose (even first)
  Affects any patient at sufficient dose           Only genetically susceptible patients
  Reversible with dose reduction                   May be irreversible / fatal
  Example: Warfarin bleeding at high INR           Example: Carbamazepine SJS/TEN

  Frequency: Common (>1%)                          Frequency: Rare (<1% overall)
                                                   but HIGH in HLA carriers (5-30%)
```

The key insight is that Type B (HLA-mediated) reactions are **rare in the general population** but **common in genetically susceptible individuals**. HLA screening identifies those susceptible individuals before drug exposure.

---

## 9. Drug Metabolism Basics

### Phase I Metabolism (CYP450)

Phase I reactions modify the drug molecule through oxidation, reduction, or hydrolysis. The CYP450 enzymes are the primary Phase I metabolizers.

```
PHASE I METABOLISM PATHWAY

Drug absorbed from GI tract
     |
     v
Portal vein --> LIVER
     |
     v
Hepatocyte (liver cell)
     |
     +-- Endoplasmic reticulum
     |     |
     |     v
     |   CYP450 enzymes (embedded in ER membrane)
     |     |
     |     v
     |   Oxidized metabolite
     |
     v
Phase II conjugation (or direct excretion if sufficiently polar)
```

### Phase II Metabolism (Conjugation)

Phase II reactions attach a water-soluble group to the drug (or its Phase I metabolite), making it easier to excrete:

```
Phase I metabolite (or parent drug)
     |
     +-- Glucuronidation (UGT enzymes) --> Glucuronide conjugate
     |
     +-- Acetylation (NAT enzymes) ------> Acetylated conjugate
     |
     +-- Methylation (TPMT, COMT) -------> Methylated conjugate
     |
     +-- Sulfation (SULT enzymes) --------> Sulfate conjugate
     |
     +-- Glutathione conjugation (GST) ---> GSH conjugate
     |
     v
Water-soluble conjugate --> Kidneys --> Urine excretion
                        --> Bile --> Fecal excretion
```

### First-Pass Metabolism

When a drug is taken orally, it passes through the liver before reaching the systemic circulation. This "first-pass effect" can significantly reduce the amount of active drug that reaches its target:

```
Oral drug --> GI absorption --> Portal vein --> LIVER --> Systemic circulation
                                                 |
                                                 v
                                           First-pass metabolism
                                           (CYP450 + Phase II)
                                                 |
                                                 v
                                           Some drug inactivated
                                           before reaching target
```

Drugs with high first-pass metabolism (e.g., codeine, propranolol) are particularly sensitive to CYP enzyme variation because a larger fraction of the dose is metabolized in the liver before reaching the bloodstream.

### Drug-Drug Interactions Affecting Metabolism

Beyond genetic variation, **drug-drug interactions** can alter CYP enzyme activity. These interactions are classified as inhibition or induction:

**Enzyme Inhibition:**

```
Normal metabolism:
Drug A --> CYP2D6 --> Metabolite A (clearance rate: 100%)

With CYP2D6 inhibitor (Drug B) added:
Drug A --> CYP2D6 (BLOCKED by Drug B) --> Drug A accumulates
Metabolite A production: REDUCED
Drug A plasma levels: INCREASED (toxicity risk)
```

Types of inhibition:
- **Competitive (reversible)**: Inhibitor competes with substrate for enzyme binding site. Effect proportional to inhibitor concentration. Resolves when inhibitor is discontinued (within 3-5 half-lives).
- **Mechanism-based (irreversible)**: Inhibitor permanently inactivates the enzyme. Effect persists until new enzyme is synthesized (days to weeks). Examples: paroxetine (CYP2D6), clarithromycin (CYP3A4).

**Enzyme Induction:**

```
Normal metabolism:
Drug A --> CYP3A4 --> Metabolite A (clearance rate: 100%)

With CYP3A4 inducer (Drug B) added:
Gene expression: CYP3A4 mRNA INCREASED (via PXR/CAR nuclear receptors)
Protein synthesis: More CYP3A4 enzyme produced
Drug A --> CYP3A4 (MORE ENZYME) --> Metabolite A (clearance rate: 200-500%)
Drug A plasma levels: DECREASED (sub-therapeutic)
```

Induction takes 7-14 days to reach full effect and 7-14 days to resolve after the inducer is stopped (time required for enzyme synthesis and degradation).

### The Therapeutic Index and PGx Significance

The **therapeutic index** (TI) is the ratio between the toxic dose and the therapeutic dose. Drugs with a narrow therapeutic index are most sensitive to pharmacogenomic variation:

```
WIDE THERAPEUTIC INDEX (e.g., amoxicillin):
  |----therapeutic range----|
  |                         |
  |  [SAFE ZONE is large]   |
  |                         |
  Minimum effective dose    Toxic dose
  (large safety margin; PGx variation less clinically significant)

NARROW THERAPEUTIC INDEX (e.g., warfarin):
  |----therapeutic range----|
  |                         |
  |  [SAFE ZONE is small]   |
  |                         |
  Minimum effective dose    Toxic dose
  (small safety margin; PGx variation HIGHLY clinically significant)
```

Drugs with narrow therapeutic indices where PGx is most important:

| Drug | Therapeutic Index | PGx Genes | Clinical Risk |
|------|------------------|-----------|---------------|
| Warfarin | Very narrow (INR 2-3) | CYP2C9, VKORC1 | Bleeding vs thrombosis |
| Tacrolimus | Narrow (trough 5-15 ng/mL) | CYP3A5 | Rejection vs nephrotoxicity |
| Phenytoin | Narrow (10-20 mcg/mL) | CYP2C9, HLA-B*15:02 | Toxicity vs seizure breakthrough |
| 5-Fluorouracil | Narrow | DPYD | Toxicity vs treatment failure |
| Thiopurines | Narrow | TPMT, NUDT15 | Myelosuppression vs inadequate immunosuppression |
| Codeine (as morphine prodrug) | Narrow | CYP2D6 | Respiratory depression vs no analgesia |

### Key Terminology for Drug Metabolism

```
SUBSTRATE:   A drug that is metabolized BY a specific enzyme
             Example: Codeine is a CYP2D6 SUBSTRATE

INHIBITOR:   A drug that REDUCES the activity of a specific enzyme
             Example: Fluoxetine is a CYP2D6 INHIBITOR

INDUCER:     A drug that INCREASES the expression of a specific enzyme
             Example: Rifampin is a CYP3A4 INDUCER

PRODRUG:     A drug that requires metabolism to become active
             Example: Codeine is a PRODRUG (requires CYP2D6 to form morphine)

METABOLITE:  The product of drug metabolism
             Example: Morphine is the active METABOLITE of codeine
```

---

## 10. Clinical Scenarios

### Scenario 1: Warfarin Initiation in a 68-Year-Old Patient

**Patient:** 68-year-old Caucasian female, 65 kg, 160 cm, newly diagnosed atrial fibrillation. No amiodarone, no enzyme inducers.

**PGx results:** CYP2C9 *1/*3, VKORC1 A/G

**Analysis:**
- CYP2C9 *3 is a decreased-function allele (~80% reduced CYP2C9 activity)
- VKORC1 A/G indicates intermediate warfarin sensitivity
- IWPC algorithm predicts ~25 mg/week (vs population average ~35 mg/week)
- Starting at the population average dose would put this patient at high risk for over-anticoagulation and bleeding

**PGx-guided action:** Start warfarin at approximately 3.5 mg/day (25 mg/week) instead of the standard 5 mg/day. Monitor INR closely during initiation. Expect stable dose to be 25-30% below population average.

**Without PGx testing:** Patient started at standard 5 mg/day. INR overshoots therapeutic range within 3-5 days. Risk of bleeding event during dose titration.

### Scenario 2: Codeine Prescribed to a CYP2D6 Ultra-Rapid Metabolizer

**Patient:** 28-year-old mother, post-cesarean delivery, prescribed codeine for pain. Breastfeeding.

**PGx results:** CYP2D6 *1/*2xN (ultra-rapid metabolizer, activity score > 2.25)

**Analysis:**
- Ultra-rapid metabolizer converts codeine to morphine at an accelerated rate
- Higher-than-expected plasma morphine levels in mother
- Morphine enters breast milk
- Infant exposed to supra-therapeutic morphine levels

**PGx-guided action:** Avoid codeine entirely. Use acetaminophen, ibuprofen, or if opioid needed, use morphine at standard doses (does not require CYP2D6 activation) with monitoring.

**Without PGx testing:** In 2006, a 13-day-old breastfed infant died from morphine toxicity after the mother was prescribed codeine post-delivery. The mother was later found to be a CYP2D6 ultra-rapid metabolizer.

### Scenario 3: Abacavir Initiation in an HIV Patient

**Patient:** 35-year-old male, newly diagnosed HIV, treatment-naive. Clinician plans to start abacavir-containing regimen (Epzicom).

**PGx results:** HLA-B*57:01 POSITIVE

**Analysis:**
- HLA-B*57:01 carriers have 5-8% risk of abacavir hypersensitivity syndrome
- Reaction includes fever, rash, GI symptoms, respiratory symptoms
- Can be fatal on rechallenge

**PGx-guided action:** Do NOT prescribe abacavir. Use alternative NRTI backbone (e.g., tenofovir alafenamide/emtricitabine).

**Without PGx testing:** Patient develops fever and rash on day 10. If not recognized as abacavir hypersensitivity and drug is continued or re-challenged, can be fatal.

### Scenario 4: 5-Fluorouracil Toxicity in a DPYD-Deficient Patient

**Patient:** 58-year-old male, stage III colon cancer, scheduled for FOLFOX chemotherapy (fluorouracil + oxaliplatin + leucovorin).

**PGx results:** DPYD *1/*2A (heterozygous carrier of the IVS14+1G>A splice site variant)

**Analysis:**
- DPYD *2A is a no-function allele (splice site abolished, no DPD enzyme produced from that allele)
- Patient has ~50% of normal DPD activity (Gene Activity Score = 1.0)
- 5-FU catabolism is reduced, leading to accumulation of cytotoxic metabolites
- Standard-dose 5-FU carries a 73% risk of grade 3+ toxicity in this patient

**PGx-guided action:** Reduce 5-fluorouracil dose by 50%. Monitor for toxicity with complete blood counts. Dose can be escalated if tolerated, guided by therapeutic drug monitoring (TDM) of 5-FU plasma levels.

**Without PGx testing:** Patient receives standard-dose 5-FU. Develops severe mucositis (unable to eat or drink), grade 4 neutropenia (absolute neutrophil count < 500), and hand-foot syndrome requiring hospitalization. Chemotherapy delayed by 4 weeks. Mortality risk: 10-20%.

### Scenario 5: Simvastatin Myopathy Risk with SLCO1B1 Variant

**Patient:** 62-year-old female, hyperlipidemia, started on simvastatin 80 mg daily.

**PGx results:** SLCO1B1 rs4149056 TC (heterozygous, *1/*5 decreased function)

**Analysis:**
- SLCO1B1 encodes the OATP1B1 transporter that moves statins into liver cells
- Decreased function variant (rs4149056 C allele) reduces hepatic statin uptake
- Simvastatin plasma levels increase 2-3 fold in heterozygous carriers
- Risk of statin-induced myopathy increases from ~1% to ~5% at high doses

**PGx-guided action:** Avoid simvastatin 80 mg dose. Use lower dose simvastatin (20-40 mg) or switch to alternative statin with lower SLCO1B1 sensitivity (rosuvastatin, pravastatin, or fluvastatin).

**Without PGx testing:** Patient develops progressive muscle weakness and elevated creatine kinase (CK) after 3 months. Risk of rhabdomyolysis with renal failure if not detected early.

### Scenario 6: CYP2C19 Poor Metabolizer Prescribed Clopidogrel After Cardiac Stent

**Patient:** 55-year-old male, acute myocardial infarction, drug-eluting stent placed. Standard discharge prescription: aspirin 81 mg daily + clopidogrel 75 mg daily for 12 months.

**PGx results:** CYP2C19 *2/*2 (poor metabolizer, activity score = 0.0)

**Analysis:**
- Clopidogrel is a prodrug that requires CYP2C19 activation to produce its active metabolite
- CYP2C19 poor metabolizers produce essentially no active metabolite from clopidogrel
- FDA Boxed Warning states: "effectiveness of Plavix is dependent on activation by CYP2C19"
- Poor metabolizers have a 3.4-fold increased risk of cardiovascular death, MI, or stroke after stent placement (Mega et al., NEJM 2009)
- Stent thrombosis risk increases from ~1% to ~3-4% in poor metabolizers

**PGx-guided action:** Do NOT use clopidogrel. Switch to alternative P2Y12 inhibitor that does not require CYP2C19 activation:
- **Prasugrel 10 mg daily** (if no history of stroke/TIA, age <75, weight >60 kg)
- **Ticagrelor 90 mg twice daily** (if prasugrel contraindicated)

**Without PGx testing:** Patient receives clopidogrel but achieves inadequate platelet inhibition. Risk of stent thrombosis -- a catastrophic event with 20-40% mortality rate.

### Scenario 7: Psychiatric Patient with CYP2D6 Ultra-Rapid Metabolizer Status

**Patient:** 32-year-old female, major depressive disorder, inadequate response to two prior antidepressants (sertraline and fluoxetine). Clinician plans to start nortriptyline (a tricyclic antidepressant).

**PGx results:** CYP2D6 *1/*2xN (ultra-rapid metabolizer, activity score = 3.0)

**Analysis:**
- Nortriptyline is metabolized primarily by CYP2D6
- Ultra-rapid metabolizers have accelerated CYP2D6 metabolism
- Standard nortriptyline doses will produce sub-therapeutic plasma levels
- CPIC guideline recommends: avoid TCAs in ultra-rapid metabolizers, or if TCA necessary, increase dose by 25% and use therapeutic drug monitoring (TDM)
- Prior treatment failures with sertraline and fluoxetine may have been partly due to rapid CYP2D6 metabolism (both are CYP2D6 substrates)

**PGx-guided action:** Consider an antidepressant that is NOT primarily metabolized by CYP2D6:
- **Escitalopram** (primarily CYP2C19, minimal CYP2D6 involvement)
- **Bupropion** (primarily CYP2B6, negligible CYP2D6 involvement)
- If TCA required, use nortriptyline with TDM; target trough 70-170 ng/mL; may need 150-200% of standard dose

**Without PGx testing:** Patient is labeled as "treatment-resistant depression" after failing a third antidepressant. The underlying pharmacokinetic reason (ultra-rapid metabolism) goes unrecognized, leading to unnecessary medication trials, prolonged suffering, and potentially inappropriate escalation to more invasive treatments (ECT, ketamine).

### Scenario 8: Allopurinol-Induced Severe Cutaneous Adverse Reaction

**Patient:** 52-year-old Southeast Asian male (Thai descent), recurrent gout, clinician plans to start allopurinol for urate-lowering therapy.

**PGx results:** HLA-B*58:01 POSITIVE

**Analysis:**
- HLA-B*58:01 is associated with allopurinol-induced SJS/TEN and DRESS syndrome
- Prevalence of HLA-B*58:01: 6-8% in Southeast Asian populations, 3-4% in African American populations, <2% in European populations
- SJS/TEN mortality rate: 10-30%
- ACR (American College of Rheumatology) conditionally recommends HLA-B*58:01 testing before allopurinol initiation in Southeast Asian and African American patients

**PGx-guided action:** Do NOT prescribe allopurinol. Use alternative xanthine oxidase inhibitor:
- **Febuxostat** (not associated with HLA-B*58:01-related SJS/TEN)
- Alternative: uricosuric agents (probenecid, lesinurad) if appropriate

**Without PGx testing:** Patient starts allopurinol 100 mg daily. After 2-6 weeks, develops fever, diffuse erythematous rash progressing to epidermal detachment (SJS/TEN). Hospitalization in burn unit required. Mortality risk: 10-30%.

---

## 10a. Understanding PGx Test Reports

### Anatomy of a PGx Lab Report

A clinical pharmacogenomic test report typically includes the following elements:

```
============================================================================
              PHARMACOGENOMIC TEST REPORT
============================================================================
Patient: [Name]                    DOB: [Date]          MRN: [Number]
Ordering Clinician: [Name]         Collection Date: [Date]
Report Date: [Date]                Lab: [CLIA# XXX]

GENE        GENOTYPE      PHENOTYPE              ACTIVITY SCORE
--------    ----------    --------------------   --------------
CYP2D6      *1/*4         Intermediate Metabolizer    1.0
CYP2C19     *1/*1         Normal Metabolizer          2.0
CYP2C9      *1/*2         Intermediate Metabolizer    1.5
CYP3A5      *3/*3         Poor Metabolizer            0.0
VKORC1      -1639 A/A     Low Warfarin Dose           --
DPYD        *1/*1         Normal Metabolizer          2.0
TPMT        *1/*1         Normal Metabolizer          --
NUDT15      *1/*1         Normal Metabolizer          --
SLCO1B1     *1/*1         Normal Function             --
HLA-B       Negative for *57:01, *58:01, *15:02       --
HLA-A       Negative for *31:01                       --

DRUG-SPECIFIC RECOMMENDATIONS:
-------------------------------
1. CODEINE: Intermediate metabolizer (CYP2D6). Use label-
   recommended dose. Monitor for reduced efficacy. Consider
   alternative analgesic if inadequate response.

2. CLOPIDOGREL: Normal metabolizer (CYP2C19). Standard dosing
   expected to be effective.

3. WARFARIN: Intermediate metabolizer (CYP2C9) + Low dose
   (VKORC1). Consider reduced initial dose. IWPC-predicted
   dose: ~28 mg/week.
============================================================================
```

### How to Read the Report

1. **Genotype column**: Shows the two alleles (diplotype) for each gene. *1 is typically the reference/normal allele. Other numbers indicate specific variants.

2. **Phenotype column**: The functional interpretation. This is the most actionable column for clinicians -- it tells you whether the patient metabolizes drugs faster or slower than expected.

3. **Activity score column**: A numeric representation of enzyme function. Higher scores generally mean more enzyme activity. Not all genes use the activity score system (e.g., VKORC1, HLA genes).

4. **Drug-specific recommendations**: The most clinically useful section. Translates genotype/phenotype into prescribing guidance for specific drugs.

### Common Pitfalls in Report Interpretation

| Pitfall | Explanation | Example |
|---------|-------------|---------|
| Ignoring phenoconversion | Report shows genotype-based phenotype, but concurrent medications may shift actual phenotype | CYP2D6 NM on report, but patient takes fluoxetine (strong CYP2D6 inhibitor) -- functionally a PM |
| Applying wrong drug to gene | Assuming all drugs in a class are metabolized by the same enzyme | Omeprazole (CYP2C19) vs pantoprazole (less CYP2C19 dependent) |
| One-time vs lifetime result | Not recognizing that germline PGx results are lifelong | Repeating a PGx test that was done 5 years ago -- results don't change |
| Missing gene-drug pairs | Assuming "no alert" means "no interaction" | Report may not cover every drug the patient takes |
| Population context | Not considering that reference allele (*1) is defined by European genomes | A patient from a non-European background may carry "missing" alleles classified as *1 |

---

## 11. Why PGx Testing Matters

### Current State of PGx Implementation

- The FDA lists pharmacogenomic biomarkers in the labeling of **350+ drugs**.
- CPIC has published guidelines for **100+ gene-drug pairs**.
- Yet fewer than **5% of patients** receive PGx testing before starting a high-risk medication.

### Barriers to Adoption

| Barrier | Description |
|---------|------------|
| Knowledge gap | Many clinicians lack training in PGx interpretation |
| Turnaround time | Traditional PGx tests take days to weeks |
| Fragmented data | Results scattered across multiple databases |
| Workflow integration | No seamless connection to EHR clinical decision support |
| Cost perception | Insurance coverage varies; pre-emptive panels cost $200-500 |
| Guideline awareness | Clinicians may not know CPIC guidelines exist for their drugs |
| IT infrastructure | EHR systems may lack PGx CDS alert capabilities |

### The Case for Pre-emptive Testing

Pre-emptive PGx testing (testing before a drug is needed) is more efficient than reactive testing (testing after a prescription is written):

- A **single multi-gene panel** covers most actionable pharmacogenes (typically 12-25 genes)
- Results are **available at the point of care** when a drug is prescribed (no turnaround time delay)
- **Cost-effective**: The average patient encounters 3-5 PGx-guided prescribing decisions over their lifetime
- **One-time test**: Germline pharmacogenomic variants do not change over a patient's lifetime (unlike tumor genomics)
- Panel cost ($200-500) is offset by preventing a single ADR-related hospitalization ($2,000-20,000)

### The PREPARE Trial Evidence

The most definitive evidence for pre-emptive PGx testing comes from the PREPARE study (Swen et al., Lancet 2023):

- 6,944 patients across 7 European medical centers
- 12-gene pre-emptive panel
- **30% reduction in ADRs** in the PGx-guided arm (OR 0.70; 95% CI 0.54-0.91)
- Results available at point of care via electronic CDS alerts
- No increase in healthcare costs

### Types of PGx Tests Available

| Test Type | Genes | Turnaround | Cost | Best For |
|-----------|-------|-----------|------|---------|
| Single-gene test | 1 gene | 1-3 days | $50-150 | Reactive testing (e.g., HLA-B*57:01 before abacavir) |
| Focused panel | 5-10 genes | 3-5 days | $100-250 | Specialty-specific (e.g., psychiatry panel) |
| Comprehensive panel | 12-25 genes | 5-10 days | $200-500 | Pre-emptive testing (recommended) |
| Whole exome/genome | All genes | 10-30 days | $500-2,000 | Research; includes non-PGx findings |

### PGx Testing Technologies

| Technology | Principle | Star Allele Coverage | Structural Variants |
|-----------|----------|---------------------|-------------------|
| SNP array (microarray) | Hybridization to known variant probes | Common alleles only | No |
| Targeted sequencing | PCR amplification + sequencing of gene regions | Comprehensive for targeted genes | Limited |
| Whole genome sequencing | Complete genome sequencing | Complete | Yes (with specialized tools) |
| Long-read sequencing | PacBio/ONT long reads through gene regions | Complete | Excellent (CYP2D6 hybrids, CNV) |

### Who Should Get PGx Testing?

Current evidence supports PGx testing in the following scenarios:

| Scenario | Priority | Genes to Test |
|----------|---------|--------------|
| Before starting warfarin | High | CYP2C9, VKORC1 |
| Before starting clopidogrel (post-stent) | High | CYP2C19 |
| Before starting abacavir | Mandatory (FDA) | HLA-B*57:01 |
| Before starting carbamazepine | High (FDA recommended) | HLA-B*15:02, HLA-A*31:01 |
| Before starting fluoropyrimidine chemo | High (EMA mandated) | DPYD |
| Before starting thiopurines | High | TPMT, NUDT15 |
| Before starting high-dose simvastatin | Moderate | SLCO1B1 |
| Pre-emptive panel (any patient) | Recommended | All 25 pharmacogenes |
| Psychiatric medication selection | Moderate | CYP2D6, CYP2C19 |
| Chronic pain management | Moderate | CYP2D6, OPRM1 |

---

## 12. How the PGx Intelligence Agent Uses These Concepts

The Pharmacogenomics Intelligence Agent operationalizes every concept in this guide:

| Concept | Agent Component | How It Works |
|---------|----------------|-------------|
| Star allele nomenclature | `StarAlleleCaller` (pgx_pipeline.py) | Resolves VCF variants to star alleles for 25 genes |
| Diplotype-to-phenotype | `PhenotypeTranslator` (pgx_pipeline.py) | Activity score summation per CPIC standardized terms |
| CPIC guidelines | `pgx_drug_guidelines` collection | 240 seed records with guideline recommendations |
| HLA screening | `HLAScreener` (hla_screener.py) | 12 HLA-drug associations with severity and alternatives |
| Drug metabolism | `DrugGeneMatcher` (pgx_pipeline.py) | Cross-references phenotype profiles against medication lists |
| Phenoconversion | `PhenoconversionDetector` (phenoconversion.py) | 30+ CYP inhibitors/inducers with phenotype shift modeling |
| Dosing algorithms | `DosingCalculator` (dosing.py) | 4 validated algorithms (IWPC warfarin, tacrolimus, fluoropyrimidine, thiopurine) |
| Population genetics | `pgx_population_data` collection | Allele frequencies across all major populations |
| Knowledge graph | `knowledge.py` (2,512 lines) | 25 pharmacogenes, 100+ drugs, 12 HLA associations |
| Evidence retrieval | `PGxRAGEngine` (rag_engine.py) | 15-collection parallel search with weighted ranking |
| Clinical reports | `export.py` | Markdown, JSON, PDF, FHIR R4 DiagnosticReport Bundle |

---

## 13. Glossary

| Term | Definition |
|------|-----------|
| **Activity Score (AS)** | Numerical value assigned to a star allele reflecting enzyme function (0.0 = no function, 0.5 = decreased, 1.0 = normal, 2.0+ = increased/duplicated). Diplotype AS is the sum of both allele scores. |
| **ADR** | Adverse Drug Reaction -- any unintended, harmful response to a medication at normal doses used for prophylaxis, diagnosis, or treatment |
| **Allele** | One of two or more versions of a gene; each person inherits one allele from each parent |
| **Allele Frequency** | The proportion of a specific allele in a population (expressed as a percentage or decimal) |
| **Biomarker** | A measurable biological characteristic used to assess health or disease status; in PGx, typically a genetic variant |
| **CDS** | Clinical Decision Support -- electronic alerts integrated into EHR systems to guide prescribing |
| **CPIC** | Clinical Pharmacogenetics Implementation Consortium -- develops evidence-based PGx prescribing guidelines |
| **CYP450** | Cytochrome P450 -- superfamily of liver enzymes responsible for metabolizing ~75% of all drugs |
| **Diplotype** | The combination of two alleles inherited for a gene (e.g., CYP2D6 *1/*4); one from each parent |
| **DPD** | Dihydropyrimidine Dehydrogenase -- enzyme encoded by DPYD gene; catabolizes fluoropyrimidine chemotherapy |
| **DPWG** | Dutch Pharmacogenetics Working Group -- European PGx guideline body |
| **DRESS** | Drug Reaction with Eosinophilia and Systemic Symptoms -- severe multi-organ hypersensitivity reaction |
| **First-Pass Metabolism** | Hepatic metabolism of a drug before it reaches systemic circulation (after oral administration) |
| **Genotype** | An individual's genetic makeup at a specific gene or variant position |
| **HLA** | Human Leukocyte Antigen -- immune system genes in the MHC region; variants associated with drug hypersensitivity |
| **IM** | Intermediate Metabolizer -- reduced enzyme activity; may need dose adjustment |
| **MHC** | Major Histocompatibility Complex -- chromosomal region (6p21) containing HLA genes |
| **NM** | Normal Metabolizer -- standard enzyme activity; typical drug response expected |
| **Pharmacogenomics (PGx)** | The study of how inherited genetic variation affects an individual's response to medications |
| **PharmGKB** | Pharmacogenomics Knowledge Base (Stanford University) -- curated database of PGx knowledge |
| **PharmVar** | Pharmacogene Variation Consortium -- maintains the definitive catalog of star allele definitions |
| **Phenoconversion** | Change in effective metabolizer phenotype due to concomitant drug-drug interaction (CYP inhibition/induction) |
| **Phenotype** | The observable effect of a genotype (e.g., metabolizer status: PM, IM, NM, UM) |
| **PM** | Poor Metabolizer -- minimal or absent enzyme activity; highest risk for drug-gene interactions |
| **Prodrug** | A medication that must be metabolized (activated) by a CYP enzyme to become pharmacologically active (e.g., codeine -> morphine via CYP2D6) |
| **RM** | Rapid Metabolizer -- slightly higher than normal enzyme activity |
| **SJS** | Stevens-Johnson Syndrome -- severe immune-mediated skin reaction with blistering (<10% BSA detachment) |
| **Star Allele** | Named variant of a pharmacogene (e.g., CYP2D6*4); defined by specific DNA changes; curated by PharmVar |
| **TEN** | Toxic Epidermal Necrolysis -- life-threatening skin detachment (>30% BSA); mortality 25-35% |
| **UM** | Ultra-Rapid Metabolizer -- excessive enzyme activity, often due to gene duplication |
| **VCF** | Variant Call Format -- standard bioinformatics file format for recording genomic variants identified by sequencing |
| **VKORC1** | Vitamin K Epoxide Reductase Complex Subunit 1 -- target of warfarin; genetic variants affect warfarin dose requirements |
| **CDS** | Clinical Decision Support -- electronic alerts integrated into EHR systems to guide prescribing decisions |
| **Copy Number Variation (CNV)** | Presence of a different number of copies of a gene than the normal two; CYP2D6 can have 0-13+ copies |
| **DILI** | Drug-Induced Liver Injury -- hepatotoxicity caused by a medication |
| **EHR** | Electronic Health Record -- digital version of a patient's paper chart |
| **Endoxifen** | Active metabolite of tamoxifen, formed via CYP2D6 metabolism; plasma levels correlate with tamoxifen efficacy |
| **FHIR** | Fast Healthcare Interoperability Resources -- a standard for exchanging healthcare information electronically |
| **G6PD** | Glucose-6-Phosphate Dehydrogenase -- enzyme whose deficiency causes hemolytic anemia with certain drugs |
| **Hapten** | A small molecule that binds to a protein to form an antigen, triggering an immune response |
| **Heterozygous** | Having two different alleles at a genetic locus (e.g., CYP2D6 *1/*4) |
| **Homozygous** | Having two identical alleles at a genetic locus (e.g., CYP2D6 *4/*4) |
| **IGNITE** | Implementing Genomics in Practice -- network of institutions implementing PGx programs |
| **INR** | International Normalized Ratio -- standardized measure of blood clotting time; target range for warfarin is typically 2.0-3.0 |
| **IWPC** | International Warfarin Pharmacogenetics Consortium -- developed the PGx-guided warfarin dosing algorithm |
| **LOINC** | Logical Observation Identifiers Names and Codes -- universal standard for identifying medical laboratory observations |
| **Myelosuppression** | Reduced bone marrow activity leading to low blood cell counts (neutropenia, anemia, thrombocytopenia) |
| **NUDT15** | Nudix Hydrolase 15 -- enzyme that dephosphorylates thioguanine nucleotides; deficiency increases thiopurine toxicity |
| **Ontogeny** | The developmental changes in enzyme expression from birth through adulthood; CYP enzymes mature at different rates |
| **OATP1B1** | Organic Anion Transporting Polypeptide 1B1 -- hepatic uptake transporter encoded by SLCO1B1; affects statin disposition |
| **Polygenic** | Involving multiple genes; polygenic PGx considers variants across many genes simultaneously |
| **RAG** | Retrieval-Augmented Generation -- AI architecture that retrieves relevant documents before generating a response |
| **rsID** | Reference SNP ID -- unique identifier for a single nucleotide polymorphism in the dbSNP database (e.g., rs3892097) |
| **SNP** | Single Nucleotide Polymorphism -- a variation at a single DNA base position |
| **Substrate** | A drug that is metabolized by a specific enzyme (e.g., codeine is a CYP2D6 substrate) |
| **TPMT** | Thiopurine S-Methyltransferase -- enzyme that inactivates thiopurine drugs; deficiency causes severe myelosuppression |
| **Therapeutic Index** | The ratio between the toxic dose and the therapeutic dose; narrow therapeutic index drugs are most sensitive to PGx variation |
| **WGS** | Whole Genome Sequencing -- sequencing of the entire genome; provides comprehensive pharmacogenomic data |

---

## 14. Visual Learning Aids

### The PGx Decision Process

```
PATIENT PRESENTS FOR NEW MEDICATION

Step 1: Is PGx testing available for this drug?
        Check CPIC guidelines (100+ gene-drug pairs)
     |
     v
Step 2: Does the patient have PGx results on file?
     |
     +-- YES --> Proceed to Step 4
     |
     +-- NO  --> Step 3: Should PGx testing be ordered?
                  Consider: drug risk, patient population,
                  available turnaround time
                  |
                  +-- Order test --> Wait for results --> Step 4
                  |
                  +-- No test --> Prescribe with standard monitoring
     |
     v
Step 4: Interpret PGx results
     |
     +-- Check metabolizer phenotype (star alleles -> AS -> phenotype)
     +-- Check HLA status (if applicable)
     +-- Check phenoconversion (current medication list)
     +-- Calculate adjusted dose (if dosing algorithm available)
     |
     v
Step 5: Apply CPIC guideline recommendation
     |
     +-- Standard dose (NM, no interactions)
     +-- Adjusted dose (IM, specific algorithms)
     +-- Alternative drug (PM or UM for prodrugs)
     +-- CONTRAINDICATED (HLA positive, PM + contraindicated drug)
     |
     v
Step 6: Document and monitor
     +-- Record PGx results in EHR
     +-- Set up appropriate monitoring
     +-- Educate patient about their PGx profile
```

### The Drug Metabolism Pathway Summary

```
ORAL DRUG ADMINISTRATION

[1] ABSORPTION (GI tract)
    |
    +-- Intestinal CYP3A4 (first-pass)
    +-- P-glycoprotein (ABCB1) efflux
    |
    v
[2] PORTAL CIRCULATION --> LIVER
    |
    +-- Phase I: CYP450 oxidation
    |   CYP2D6, CYP2C19, CYP2C9, CYP3A4, CYP1A2, CYP2B6
    |   Result: More polar metabolite
    |
    +-- Phase II: Conjugation
    |   UGT1A1, TPMT, NAT2, DPYD
    |   Result: Water-soluble conjugate
    |
    v
[3] SYSTEMIC CIRCULATION
    |
    +-- DISTRIBUTION to target tissue
    |   Affected by: protein binding, transporters (SLCO1B1)
    |
    +-- PHARMACODYNAMIC EFFECT
    |   Affected by: receptor variants (VKORC1, OPRM1)
    |
    v
[4] ELIMINATION
    |
    +-- Renal excretion (kidneys)
    +-- Biliary excretion (bile --> feces)
    +-- Enterohepatic recirculation (microbiome)

PHARMACOGENOMIC VARIATION CAN AFFECT EVERY STEP
```

### The 25 Pharmacogenes at a Glance

```
PHASE I ENZYMES (CYP450):
  CYP2D6  *** (most polymorphic; 25% of drugs; opioids, antidepressants)
  CYP2C19 *** (clopidogrel, PPIs, antidepressants)
  CYP2C9  *** (warfarin, phenytoin, NSAIDs)
  CYP3A4  **  (most abundant; statins, immunosuppressants)
  CYP3A5  **  (tacrolimus dosing)
  CYP2B6  *   (efavirenz, methadone)
  CYP1A2  *   (caffeine, clozapine, theophylline)
  CYP4F2  *   (vitamin K metabolism; warfarin fine-tuning)

PHASE II ENZYMES:
  UGT1A1  **  (irinotecan toxicity)
  NAT2    *   (isoniazid acetylation)
  TPMT    *** (thiopurine dosing)
  DPYD    *** (fluoropyrimidine toxicity)

TRANSPORTERS:
  SLCO1B1 **  (statin myopathy)
  ABCB1   *   (P-glycoprotein; drug absorption/distribution)

HLA GENES:
  HLA-A   *** (carbamazepine DRESS, ticlopidine hepatotoxicity)
  HLA-B   *** (abacavir, carbamazepine, allopurinol SJS/TEN)
  HLA-DRB1 *  (lapatinib hepatotoxicity)

DRUG TARGETS & OTHER:
  VKORC1  *** (warfarin target; dose sensitivity)
  NUDT15  **  (thiopurine toxicity in East Asians)
  G6PD    **  (hemolytic anemia with certain drugs)
  IFNL3   *   (hepatitis C treatment response)
  RYR1    *   (malignant hyperthermia risk)
  CACNA1S *   (malignant hyperthermia risk)
  CYP2A6  *   (nicotine metabolism)
  COMT    *   (catecholamine/pain metabolism)

*** = High clinical impact, well-established CPIC guidelines
**  = Moderate impact, emerging or specific guidelines
*   = Lower impact or emerging evidence
```

---

## 15. Frequently Asked Questions (FAQ)

### Patient-Facing Questions

**Q: Do I need to repeat PGx testing?**
A: No. Germline pharmacogenomic variants are inherited and do not change over a patient's lifetime. A PGx test taken at any age provides results valid for the rest of the patient's life. This is fundamentally different from tumor genomic testing, which reflects acquired somatic mutations that evolve over time.

**Q: Will my PGx results change if I gain or lose weight?**
A: No. PGx genotype and metabolizer phenotype are determined by inherited DNA. However, weight does affect drug dosing calculations (e.g., the IWPC warfarin algorithm includes body surface area as a variable). A patient's PGx-guided dose recommendation may change with significant weight change, but the underlying genotype does not.

**Q: Can PGx testing tell me if a drug will work for me?**
A: PGx testing predicts how your body processes a drug (pharmacokinetics) and, in some cases, how your body responds to a drug (pharmacodynamics). It can identify drugs that are likely to be ineffective (e.g., codeine in a CYP2D6 poor metabolizer) or toxic (e.g., 5-FU in a DPYD-deficient patient). However, it does not guarantee that a drug will work -- drug response depends on many factors beyond genetics, including disease severity, adherence, diet, and concurrent medications.

**Q: Is PGx testing covered by insurance?**
A: Coverage varies by payer and indication. FDA-required tests (HLA-B*57:01 before abacavir) are generally covered. Pre-emptive panels are increasingly covered, especially when ordered for high-risk prescribing situations (post-stent antiplatelet selection, psychiatric medication selection). Medicare and many commercial payers now cover multi-gene PGx panels when ordered by a qualified prescriber with documented clinical indication.

**Q: Should I get a direct-to-consumer (DTC) PGx test?**
A: DTC tests (e.g., 23andMe) provide limited PGx information. They typically test only selected SNPs for 2-5 genes and do not provide the comprehensive star allele coverage, structural variant detection, or clinical interpretation that CLIA-certified clinical PGx panels offer. DTC results should not be used for clinical decision-making without confirmation by a clinical-grade test.

### Clinician-Facing Questions

**Q: Which patients should I prioritize for PGx testing?**
A: Highest priority: patients about to start a high-risk PGx-actionable drug (fluoropyrimidines, thiopurines, clopidogrel post-stent, abacavir, carbamazepine). Moderate priority: patients with a history of ADRs or treatment failure, polypharmacy patients (5+ medications), and patients starting psychiatric medications. Low-impact: patients on medications with low PGx sensitivity.

**Q: How do I integrate PGx results into my prescribing workflow?**
A: The recommended workflow is: (1) Order a pre-emptive multi-gene PGx panel, (2) Store results in the EHR PGx section (not just a PDF in the chart), (3) Enable CDS alerts that fire when a PGx-actionable drug is prescribed, (4) Consult CPIC guidelines (or the PGx Intelligence Agent) for specific dosing/drug recommendations, (5) Document PGx-guided decisions in the progress note.

**Q: What if my patient has a genotype not covered by CPIC guidelines?**
A: CPIC guidelines assign functional status to well-characterized alleles. For novel or rare alleles (classified as "uncertain function" or not in PharmVar), clinical judgment is required. The PGx Intelligence Agent flags uncertain-function alleles and provides context from the literature, but a definitive recommendation may not be possible. Consider consultation with a clinical pharmacogenomics specialist.

**Q: How long do I need to wait after stopping a CYP inhibitor before the patient's phenotype reverts to their genetic baseline?**
A: It depends on the type of inhibitor. Competitive (reversible) inhibitors: phenotype reverts within 3-5 half-lives of the inhibitor (typically 1-3 days). Mechanism-based (irreversible) inhibitors: recovery requires synthesis of new enzyme, typically 2-4 weeks. The PGx Intelligence Agent's phenoconversion module accounts for both types.

---

## 16. Resources for Further Learning

### Online Resources

| Resource | URL | Description |
|----------|-----|-------------|
| CPIC Guidelines | cpicpgx.org | Peer-reviewed, evidence-based PGx prescribing guidelines |
| PharmGKB | pharmgkb.org | Curated PGx knowledge base (clinical annotations, pathways, drug labels) |
| PharmVar | pharmvar.org | Definitive catalog of star allele definitions for all pharmacogenes |
| FDA PGx Biomarkers | fda.gov/drugs/table-pharmacogenomic-biomarkers-drug-labeling | Table of 350+ drugs with PGx biomarkers in labeling |
| DPWG Guidelines | knmp.nl | Dutch PGx guidelines (complementary to CPIC) |
| PharmCAT | pharmcat.org | Open-source star allele caller and CPIC guideline annotator |

### Textbooks and Key References

1. **Pharmacogenomics: Challenges and Opportunities in Therapeutic Implementation** (Relling & Evans, 2015) -- foundational review of PGx implementation science
2. **CPIC Guideline Series** (2011-present) -- the definitive clinical PGx guidelines; 25+ gene-drug pair publications
3. **Standardizing CYP2D6 Genotype to Phenotype Translation** (Caudle et al., 2020) -- the activity score consensus paper
4. **The PREPARE Study** (Swen et al., Lancet 2023) -- landmark pre-emptive PGx trial demonstrating 30% ADR reduction
5. **Addressing Phenoconversion: The Achilles' Heel of Personalized Medicine** (Shah & Smith, 2015) -- comprehensive review of phenoconversion challenges
