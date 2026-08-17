# The 17 Demonstrations

One per engine, per agent, and for the TSC disease program. Each is written to the same shape: a
person at stake, the factory compressing their odyssey, and an honest statement of what the room is
actually seeing.

**Every demo carries exactly one label.** This is not decoration — it is the project's honesty gate:

| Label | Meaning |
|---|---|
| **LIVE** | Running now, in front of the audience, on real input, on this box |
| **REPRESENTATIVE** | A pre-computed or curated result standing in for a long or gated step — said out loud |
| **BURST** | Running live, but on remote GPUs over the private mesh — say "elastic burst", never "all on one box" |

**Current distribution:** **13 LIVE · 3 REPRESENTATIVE · 1 MIXED.** The three representative ones
(E1, E3, E7) wait on gated software — Parabricks, the BioNeMo NIMs, and a CUDA-served ESMFold —
and flip to LIVE when it lands.

**A LIVE label is a claim about the demo, not a promise the box is up.** Every LIVE demo needs the
platform running, Milvus seeded, and `ANTHROPIC_API_KEY` set. Ask the runner rather than trusting
this page:

```bash
.venv/bin/python scripts/run_demo.py --check-all
```

It refuses to run a demo labelled LIVE whose service is unreachable rather than returning a canned
result — so the honest count on a cold box is lower than the count above, and it says so.

Every demo states **decision support, not diagnosis** the first time clinical output appears.

---

## Engines

### E1 · Genomic Foundation — "The variant that was always there"
**Label: REPRESENTATIVE** (LIVE after Parabricks · G2)
**Audience:** clinicians, general · **Runtime:** 6 min

*Weight.* A child has had three years of tests and no answer. Their genome was sequenced eighteen
months ago; the answer was in the file the whole time.

*Compression.* Load HG002 — a **publicly consented** reference sample, never a patient — into the
variant store. Query the region. Apply ACMG secondary-findings logic. Show the variant surface: SNVs,
indels, CNVs, structural variants, and the transition/transversion ratio as a quality signal.

*What actually runs:* DuckDB variant store, ACMG SF module, GWAS association (156 tests green).
*What does not:* alignment and variant calling. Parabricks is not installed, so FASTQ→BAM→VCF is
**pre-computed**. Say so. The "under five hours" claim belongs to Parabricks, not to this demo.

*Hope.* The same store answers the next child in seconds.

---

### E2 · Precision Intelligence — "Ask the evidence layer a question"
**Label: LIVE** (needs `ANTHROPIC_API_KEY` + seeded Milvus)
**Audience:** clinicians · **Runtime:** 5 min

Ask in plain English: *"What variants are associated with frontotemporal dementia?"* The engine
retrieves across ClinVar and AlphaMissense vectors, reasons over them, and returns cited evidence —
including **VCP on chromosome 9**, the thread E3 and E7 pick up.

*Impressive because* the citation is checkable. Click through to the source record.
*Honesty:* retrieval quality depends on what is indexed; an empty collection returns nothing, and
the demo should show that failure mode rather than hide it.

---

### E3 · Therapeutic Discovery — "From one protein to a hundred candidates"
**Label: REPRESENTATIVE** (LIVE/BURST after MolMIM + DiffDock · G3/G4)
**Audience:** general — this is the flagship · **Runtime:** 8 min

Seed compound CB-5083 against **p97/VCP**. Generate novel analogues, dock them, score, rank.

*What actually runs today:* RDKit property calculation, the ranking and filtering logic (124 tests).
*What does not:* MolMIM generation and DiffDock posing — **252 combined code references, neither
installed**. The molecules and poses shown are pre-computed. If the NIMs are x86-only, this becomes
**BURST** and the narration must say "elastic burst", not "one box".

*Do not run this as LIVE until §4 of the gated build guide is complete.* It is the most-watched
demo and the easiest to overclaim.

---

### E4 · Clinical Imaging — "The scan that had already answered"
**Label: LIVE** · **Audience:** cardiologists, radiologists · **Runtime:** 7 min

Coronary CT → measured stenosis → **CAD-RADS 2.0** category with modifiers (P, HRP) → cross-modal
hand-off to genomics.

*This is the strongest demo in the catalogue* and the only one already proven end-to-end on this box:
API on :8524, 1,365 tests green, accuracy verified 2026-08-13. It reads a scan, produces a graded
report, then names the genomic follow-up — 9p21, LPA — and hands to **Engine 2 · Precision
Intelligence**.

*Honesty inline:* decision support for a qualified clinician; CAD-RADS is a reporting standard, not
a diagnosis.

---

### E5 · Precision Oncology — "The molecular tumour board, in one afternoon"
**Label: LIVE** (needs seeded Milvus) · **Audience:** oncologists · **Runtime:** 7 min

A pediatric solid-tumour case (`demo/pediatric_oncology_case.json`). Variant → actionable target →
therapy options → evidence tier. 556 tests green.

*Impressive because* it shows the **tier**, not just the answer: what is FDA-approved, what is
off-label, what is trial-only. Clinicians trust the gradation more than the recommendation.

---

### E6 · Cardiology — "Risk that changes management"
**Label: LIVE** (needs seeded Milvus) · **Audience:** cardiologists · **Runtime:** 5 min

`demo/cardiology_risk.json` → structured risk with the reasoning shown. 1,966 tests green — the
deepest suite in the platform.

*Pair with E4.* Imaging finds the lesion; cardiology contextualises the risk; pharmacogenomics (A3)
decides whether the statin will work. Three engines, one patient — that sequence is the point.

---

### E7 · Structural Biology — "The shape you have to fit"
**Label: REPRESENTATIVE** (LIVE after ESMFold · G5)
**Audience:** protein scientists · **Runtime:** 6 min

**PDB 5FTK** — cryo-EM p97/VCP with **ADP bound** (a nucleotide, not a drug). Show the druggable
pockets a molecule must fit.

*What actually runs:* protein search and re-ranking, developability scoring (34 tests).
*What does not:* ESMFold prediction and ProteinMPNN design — both need CUDA, which the environment
serving them does not have.
Structures shown are **deposited PDB entries**, not predictions made in the room.

> The engine capability `structural-biology-engine:8581` is `planned` — nothing binds that port.
> Demo the **model-level** services (8570/8571/8578), not the engine endpoint.

---

### E8 · Single-Cell — "Nine populations from one sample"
**Label: LIVE** · **Audience:** researchers · **Runtime:** 5 min

Real scanpy on the bundled **PBMC 3k**: QC → normalise → HVG → PCA → neighbours → Leiden → marker DE,
then annotate clusters by marker overlap. Verified output: **2,700 cells → 9 clusters →** CD4 T,
CD8 T, B, NK, CD14+ and FCGR3A+ monocytes, dendritic, megakaryocytes.

*Impressive because it is small and completely real* — no gated model, no GPU, runs on the Grace
cores in minutes. After E3's honesty caveats, this one lands harder precisely because nothing is
withheld.

*Prerequisite:* `pip install scanpy anndata` (ungated, step 0 of the gated build guide).

---

## Intelligence Agents

Each agent reasons over a curated corpus and returns cited output. All eight are **LIVE** once
Milvus is seeded and `ANTHROPIC_API_KEY` is set; none needs gated software. Payloads already exist
in `demo/requests/`.

### A1 · CAR-T — "Why this construct, for this patient"
**LIVE** · 6 min · `cart_query.json`
Eleven collections — constructs, assays, manufacturing, safety, biomarkers, regulatory. Show
construct selection *with its safety counterweight*: on-target off-tumour risk and CRS/ICANS
monitoring. 415 tests. **The reference implementation** — the only service using the governance gate.

### A2 · Precision Biomarker — "The marker that changes the decision"
**LIVE** · 5 min · `biomarker_query.json` · 709 tests
Predictive vs prognostic, stated explicitly. Clinicians conflate them; the agent does not.

### A3 · Pharmacogenomics — "Two patients, same dose, different outcome"
**LIVE** · 6 min · `pharmacogenomics_query.json` · 1,001 tests
**The most rigorous content in the platform** — 44 CPIC pairs across 13 genes, all verified genuine.
Show CYP2C19 → clopidogrel, or DPYD → fluorouracil where the stakes are toxicity, not efficacy.
*Honesty:* CYP3A4 substrate relationships are pharmacology, not CPIC guidance — label them so.

### A4 · Precision Autoimmune — "Before the third flare"
**LIVE** · 5 min · `autoimmune_query.json` · 455 tests

### A5 · Neurology — "NIHSS, and what comes next"
**LIVE** · 5 min · `neurology_nihss.json` · 208 tests
Now on UI **8535** / API **8536** after the port re-seat.

### A6 · Clinical Trial — "The trial that was open all along"
**LIVE** · 6 min · `trial_match.json` · 769 tests
*Close the loop:* hand E3's discovered molecule to this agent and ask what it would need to enter
trials. That hand-off is the single most effective closing move in the catalogue.

### A7 · Rare Disease — "Ending the odyssey"
**LIVE** · 6 min · `rare_disease_diagnose.json` · 206 tests
HPO phenotype terms → ranked differential with evidence. The demo closest to the mission.

### A8 · Single-Cell Intelligence — "The engine computes, the agent interprets"
**LIVE** · 5 min · `single_cell_query.json` · 185 tests
**Run immediately after E8.** The engine produced nine clusters; this agent explains what they mean
clinically. The clearest illustration of the engine/agent split in the whole platform.

---

## Disease Program

### P1 · Tuberous Sclerosis Complex — "The whole factory, one child"
**Label: MIXED — LIVE agents, REPRESENTATIVE discovery** · **Audience:** everyone · **Runtime:** 12 min

The flagship. An infant with seizures and cardiac rhabdomyomas. Five disease-specific agents
composing the horizontal engines: variant curator → trajectory modeller → therapeutics strategist →
phenome mapper → TAND surveillance.

*Weight.* The need for this program came from paediatric clinicians. TSC1/TSC2 act as a **brake** on mTOR; lose the
brake and growth runs unchecked.

*Compression.* Genome → variant curated → mTOR pathway → everolimus considered with its real
tolerability profile → trial matched → TAND surveillance scheduled.

*Honesty, inline and non-negotiable:*
- decision support for a qualified clinician, never diagnosis or prescribing;
- **gene-therapy correction of TSC1/TSC2 is preclinical** — the open design bench, not a cure today;
- pediatric caution at full force.

*Hope.* Open-source, one box, and the next child's afternoon is faster because of this one.

*Public wording:* **"No one should wait years for a disease we could have understood in a day."**
Not the internal phrasing.

---

## Honest summary

| Label | Count | Which |
|---|---:|---|
| **LIVE** | 6 | E2, E4, E5, E6, E8, + the 8 agents once Milvus is seeded |
| **REPRESENTATIVE** | 9 | E1, E3, E7 and the agents until Milvus is seeded |
| **MIXED** | 2 | P1, and E3 if the NIMs burst |

The honest count improves in one step: execute `GATED_SOFTWARE_BUILD_GUIDE.md`, and E1, E3 and E7
become LIVE or BURST.

**The strongest three to show today:** **E4** (imaging — proven, end-to-end, on this box),
**E8** (single-cell — small and entirely real), **A3** (pharmacogenomics — 44 verified CPIC pairs).

**The one to be most careful with:** **E3**. It is the flagship, the most watched, and currently the
least live.

---

## Relationship to the D1–D7 portfolio

**These are a different axis, not a replacement.** `docs/demos/index.md` defines the established
**D1–D7 portfolio**: seven *patient-story* demonstrations (D3 is the TSC flagship, D1 is secondary
genomics → novel molecule, D2 is imaging → cardiac → pharmacogenomics, and so on). That portfolio is
the coverage contract the honesty framework and the public site refer to.

This catalogue is keyed **per subject** — one demo for each of the 8 engines (E1–E8), 8 agents
(A1–A8) and the TSC program (P1) — because the brief asked for one per capability.

| Axis | Keys | Organising principle | Use it for |
|---|---|---|---|
| **Portfolio** (`docs/demos/index.md`) | D1–D7 | one patient, many capabilities | proving the factory to a room |
| **Catalogue** (this document) | E1–E8, A1–A8, P1 | one capability, shown properly | proving a single subject works |

They compose: portfolio **D3** (TSC flagship) is the same program as catalogue **P1**, shown from the
patient's side rather than the capability's. Portfolio **D1** draws on catalogue **E1 → E2 → E3**.

**Do not renumber D1–D7.** It is referenced by `docs/honesty/maturity-matrix.md`, the site, and the
demo-foundation alignment. This catalogue was re-keyed to E/A/P precisely to avoid colliding with it.
