---
name: engines-agents-disease-programs
description: >-
  The canonical capability picture of the HCLS AI Factory — the authoritative, current roster of the
  8 engines, 8 intelligence agents, and the 1 disease program (TSC), with ports and what each does.
  Consult whenever naming or counting engines/agents, or creating ANY artifact (doc, figure, demo,
  video, README, deck) that must reflect the canonical framing — so nothing reverts to a stale count
  (e.g. "six engines") or invents a capability. A north star for WHAT the factory is; pairs with
  hcls-core-vision-mission (WHY) and demo-foundation-alignment (how it's proven).
---

# Eight Engines · Eight Intelligence Agents · One Platform — the canon

The authoritative roster. Patient DNA → therapeutic candidates, in hours, on a single NVIDIA DGX
Spark (~$4,699), open-source under Apache-2.0. "One box, **elastic burst**" — heavy / ARM-
incompatible models run on RunPod over a private mesh and register as native capabilities. **All
clinical outputs are decision support for a qualified clinician — never autonomous diagnosis or
prescribing.** Full detail: the "8 Engines - 8 Agents Focus" paper.

**The model — Platform → Engine → Agent → Tool:** *engines* are the 8 horizontal compute
capabilities; *agents* are the 8 clinical reasoning layers over them; *disease programs* are the
verticals that compose them for one condition; *One Platform* is the glue (registry, MCP, composer,
MLOps, governance).

## The 8 Engines
| # | Engine | Port(s) | What it does |
|---|---|---|---|
| 1 | Genomic Foundation | 5000 (+ store 8575) | GPU variant calling + annotation (Parabricks · DeepVariant + HaplotypeCaller · ClinVar/AlphaMissense · DuckDB Ts/Tv · ACMG SF · VAF/mosaicism · GWAS) |
| 2 | Precision Intelligence | 5001 | RAG clinical interpretation over the variant foundation; coordinates the 8 agents |
| 3 | Therapeutic Discovery | 8505 (+8574,8572) | Small-molecule design (MolMIM/BRICS → DiffDock → ADMET → generate-score-reseed; Chai-1 co-fold) |
| 4 | Clinical Imaging | 8525 | DICOM analysis (VISTA-3D · VILA-M3 · FHIR R4) + cross-modal reasoning with genomics |
| 5 | Precision Oncology | 8526 | MTB packets, therapy ranking, trial matching; fusion-first pediatric |
| 6 | Cardiology | 8536 | 11 workflows + 6 risk calculators (prevention → intervention → rhythm; CAC + ASCVD) |
| 7 | Large-Molecule / Structural Biology ⭐ | 8570–8578 | ESMFold · ESM-2 + Smith-Waterman · ProteinMPNN · developability · MHCflurry · ESM-2 LoRA · Chai-1 |
| 8 | Single-Cell Compute ⭐ | 8573 | scanpy compute → cell-type clusters + TME map (shared service) |

## The 8 Intelligence Agents
| Agent | Domain | Port | Key capabilities |
|---|---|---|---|
| CAR-T Intelligence | Cell therapy | 8521 | cross-collection evidence, comparative analysis, deep research |
| Precision Biomarker | Biomarkers | 8528 | PhenoAge/GrimAge, 9-domain risk, genotype-aware; multi-omics join |
| Pharmacogenomics | Drug–gene | 8507 | star-allele calling, CPIC, 9 dosing algorithms (safety interlock) |
| Precision Autoimmune | Autoimmune | 8531 | autoantibody interpretation, HLA analysis, flare prediction |
| Neurology Intelligence | Neurology | 8529 | stroke triage, dementia eval, EDSS; Parkinson's S+N+G staging |
| Clinical Trial Intelligence | Clinical trials | 8128 | trial optimization, adaptive design, biomarker strategy |
| Rare Disease Diagnostic | Rare disease | 8544 | HPO matching, ACMG classification, trio, gene-therapy tracking |
| Single-Cell Intelligence | Single-cell | 8130 | cell-type annotation, TME profiling, drug-response prediction |

**Single-cell — engine vs. agent (not a duplicate):** the *engine* (Engine 8, `singlecell-compute`,
:8573) is deterministic scanpy *compute*; the *agent* (`single-cell-intelligence-agent`, :8130) is
the RAG *reasoning* layer. The engine computes; the agent interprets. **Presentation:** in lists
name them **"Single-Cell Compute Engine"** and **"Single-Cell Intelligence Agent"** (the Compute vs.
Intelligence contrast reads as two roles, not a duplicate); in diagrams **stack them as one
capability** — a reasoning cap on a compute base — never two identical side-by-side nodes.

## The 1 Disease Program — Tuberous Sclerosis Complex (flagship)
Disease programs are the verticals the horizontal foundation powers; **TSC is the first — the
pediatric clinical beachhead.** A dedicated TSC engine + **five disease-specific agents** (variant
curator · trajectory modeler · therapeutics strategist · phenome mapper · TAND surveillance),
composing the horizontal engines/agents. Chosen because it is multi-system (brain/kidney/heart/
skin/lungs → lights up the whole factory) and a clean gene→drug story: **TSC1/TSC2 loss →
mTORC1 hyperactivation → everolimus** (FDA-approved SEGA 2010, renal AML 2012, seizures 2018). The
flagship follows the **weight → compression → hope** arc over one infant, one governed afternoon.
On disk: a self-contained, portable vertical under `core/disease-programs/`. *Honesty:* TSC1/TSC2
gene therapy is **preclinical**; decision support; pediatric caution at full force. **Replication
roadmap:** NF1, NF2, Rett, Williams, the broader mTORopathies — same pattern, same foundation.

## One Platform
Capability Registry (typed source of truth; a `live` capability is never `mock`-served) · MCP
tool-surface · AI Workflow Composer (NL → validated governed pipeline) · single-box MLOps · the
governance gates · shared infra (vector database, event bus, monitoring, multi-omics join) ·
on-demand remote GPU (RunPod / Tailscale).

## Frontier models
**Chai-1** (registered `planned`, build in progress) — AlphaFold3-class complex co-folding,
Apache-2.0 code+weights, commercial-cleared; flips to `live` when stood up. **Chai-2** (planned,
gated/partnership) — zero-shot de novo binder design for the CAR-T / autoimmune binder family.

## Standing honesty ledger (carry everywhere)
Gene therapy (TSC1/TSC2 & general) = preclinical · MAISI synthetic imaging = research/QA, never
diagnostic · single-cell atlas similarity + cell embeddings = roadmap · Chai-2 = gated · αSyn-SAA /
NSD-ISS / SynNeurGe = research/trial-use · pediatric RNA fusion detection = recommended addition ·
"elastic burst" not "all on one box" · **all clinical output = decision support, not diagnosis.**

## Do / Don't
**Do:** use this exact roster and these counts (8 engines, 8 agents, 1 disease program) in every
doc/figure/demo/video; name engines/agents as here; carry the honesty labels; verify against the
capability registry (`lib/hcls_common/capabilities.json`) when in doubt.
**Don't:** reproduce a stale "six engines" (or any other) count; invent a capability, port, or
model not in the roster/registry; imply a preclinical/roadmap/gated item ships; treat single-cell's
engine and agent as one thing.

## Related
- `hcls-core-vision-mission` — the mission (WHY) this capability set serves.
- `demo-foundation-alignment` — the D1–D7 demos that prove these capabilities.
- `clinical-claim-honesty` — the honesty ledger every use of this roster carries.
- `figurelabs-illustrations` · `demo-script-authoring` — artifacts that must match this canon.
