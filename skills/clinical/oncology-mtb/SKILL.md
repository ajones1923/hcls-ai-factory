---
name: oncology-mtb
description: >-
  The molecular-tumor-board standards that keep the Precision Oncology Engine (E5, 8526) clinically
  correct and honestly framed — with pediatric oncology built in as a first-class, DIFFERENT case,
  not an afterthought. Consult whenever building, reviewing, or demoing MTB packet generation,
  therapy ranking, trial matching, or cross-modal joins (esp. the D4 pediatric MTB). Every packet
  is decision support for a molecular tumor board — never an autonomous treatment decision.
---

# Molecular Tumor Board Standards — pediatric-aware, decision support only

The Precision Oncology Engine (`localhost:8526`) assembles the **MTB packet**: variant/fusion
findings, ranked therapy options, matched trials, and cross-modal context (imaging, single-cell TME,
germline). It is **decision support prepared for a qualified molecular tumor board — never an
autonomous treatment or prescribing decision.** The load-bearing rule: **pediatric oncology is
biologically different from adult and must be built in, not bolted on.**

## Adult MTB backbone (the baseline)

- Rank therapies against the driver alteration with a **level-of-evidence tier** (align to
  **AMP/ASCO/CAP tiers I–IV; Li et al., *J Mol Diagn* 2017; PMID 27993330**, and OncoKB levels).
  Carry the tier on every recommendation — an FDA on-label match and a preclinical rationale are not
  the same claim.
- Distinguish **on-label / off-label / trial-only / preclinical** for each option, and surface
  **resistance biology** (see below), not just first-line matches.
- Match trials from ClinicalTrials.gov with real eligibility logic (biomarker, line of therapy, age).

## Pediatric MTB — the differences that must be encoded

1. **Fusion-first, not mutation-first.** Pediatric tumors are frequently driven by **gene fusions
   and structural events** on a quiet mutational background, so a mutation-centric adult pipeline
   will miss the driver. Encode the canonical drivers:
   - **ALK amplification/mutation & MYCN amplification** — high-risk neuroblastoma
   - **EWSR1-FLI1** (and EWSR1 variants) — Ewing sarcoma
   - **PAX3-FOXO1 / PAX7-FOXO1** — alveolar rhabdomyosarcoma
   - **KMT2A (MLL) rearrangements** — infant / high-risk ALL
   - **NTRK1/2/3 fusions** — tissue-agnostic → **larotrectinib / entrectinib**
   - **H3 K27M** — diffuse midline glioma (DMG)
   - **SMARCB1 / INI1 loss** — rhabdoid tumors (ATRT, MRT)
2. **Germline is first-class.** In pediatric oncology a large share of cases carry a pathogenic
   germline predisposition variant that **changes treatment and triggers family screening**:
   **TP53 (Li-Fraumeni** — avoid/limit radiation and radiomimetics), **RB1, DICER1, NF1**, and
   others. The engine must join **E1 germline** findings (see `clinical-genomics-standards`) into the
   packet — a somatic-only MTB is incomplete for a child.
3. **Late-effects-weighted ranking.** For a patient with decades of survivorship ahead, ranking
   weighs **long-term toxicity**, not just response: **anthracycline cardiotoxicity, neurocognitive
   impact (cranial radiation), growth and endocrine effects, fertility preservation, and
   secondary-malignancy risk** (acute against a lifetime). Two equally active options are **not**
   equally ranked if one mortgages the child's future.
4. **Pediatric-appropriate trials.** Match to the pediatric trial ecosystem — **NCI-COG Pediatric
   MATCH** and **Children's Oncology Group (COG)** protocols — not adult basket trials the child is
   ineligible for. Age/weight/developmental eligibility is part of the match.

## Resistance scenarios to encode

- **Ph+ ALL with BCR-ABL1 T315I** → **ponatinib** (the T315I-active TKI when imatinib/dasatinib fail).
- **ALK-driven neuroblastoma, crizotinib resistance** → **lorlatinib** (next-generation ALK
  inhibitor with CNS activity).

The packet should show the resistance path, not present a first-line agent as if resistance won't
emerge.

## Honest gap — state it out loud

**Pediatric MTBs lean on RNA fusion detection** (**Arriba / STAR-Fusion over RNA-seq**) to catch the
fusion drivers above; robust RNA-fusion calling is a **recommended near-term registry addition**, not
yet a fully wired live capability. Flag it plainly so a **St. Jude / Cincinnati Children's** demo of
D4 is airtight rather than overclaimed — this is exactly the item carried in the
`demo-foundation-alignment` honesty ledger. Do not present fusion coverage as complete until that
capability is registered, verified, and live.

## Do / Don't

**Do:** carry AMP/ASCO/CAP + OncoKB evidence tiers on every therapy; switch to **fusion-first** logic
and encode the pediatric drivers for any pediatric case; join germline predisposition (TP53/RB1/
DICER1/NF1) and trigger family-screening notes; weight ranking by **late effects**; match COG /
Pediatric MATCH trials; show resistance paths (T315I→ponatinib, ALK→lorlatinib); state the RNA-fusion
gap; keep the packet as MTB decision support.

**Don't:** run a mutation-only pipeline on a child and miss the fusion driver; drop germline findings
from the packet; rank purely on response while ignoring cardiotoxicity/neurocognition/fertility/
secondary-malignancy risk; match a child to an adult-only trial; imply RNA-fusion detection is fully
live when it's a roadmap addition; present any ranking as an autonomous treatment decision.

## MTB packet checklist

1. Somatic **and** germline findings included; predisposition genes flagged with family-screening note?
2. Pediatric case → fusion-first logic run; canonical drivers screened (ALK/MYCN, EWSR1-FLI1,
   PAX3/7-FOXO1, KMT2A, NTRK, H3K27M, SMARCB1)?
3. Every therapy carries an AMP/ASCO/CAP + OncoKB evidence tier and on-/off-label/trial/preclinical
   status?
4. Ranking weights late effects for pediatric survivorship, not response alone?
5. Trials matched to the correct ecosystem (COG / Pediatric MATCH for children) with real eligibility?
6. Resistance scenarios surfaced where relevant?
7. RNA-fusion coverage stated honestly (roadmap vs live)? Packet framed as MTB decision support?

## Related
- `clinical-claim-honesty` — the honesty register the RNA-fusion gap and evidence tiers pass.
- `clinical-genomics-standards` — germline classification (Li-Fraumeni etc.) feeding the packet.
- `pharmacogenomics-cpic` — chemo/targeted-agent safety interlocks the ranking must respect.
- `demo-foundation-alignment` — D4 (pediatric MTB) demo home and the RNA-fusion honesty-ledger item.
