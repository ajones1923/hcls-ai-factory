# Cardiac Workup — imaging to genomics walkthrough (technical cut)

**Demo:** DEMO-003, Cardiac Workup — Coronary Artery Disease Assessment (Engine 4, Clinical Imaging)
**Run time:** ~12 minutes · **Room:** clinicians, builders, skeptics
**Surfaces:** portal `:3002` → Workflows · API `:8524`
**Visual version:** published artifact (screens embedded), see `tmp/ct-23.png` … `tmp/ct-29.png`

Built to the three-beat arc — weight → compression → hope. Honesty flags are **inline, spoken in
the beat they apply to**, never saved for a closing disclaimer.

---

## Pre-flight

1. **Milvus first.** `docker ps` must show `imaging-milvus-standalone` healthy. It restarts itself;
   the API and portal do not.
2. **Start API + portal**, then confirm the log reads `Milvus connected` — not merely that
   `/health` returns `healthy`. With Milvus down the API still reports healthy and reads every
   collection as 0.
3. **Warm the index.** First search after a restart takes ~15 s (collections load into memory);
   every one after is ~100 ms.
   ```bash
   curl -s -X POST http://127.0.0.1:8524/search -H 'Content-Type: application/json' \
     -d '{"question":"warmup","top_k":1}' > /dev/null
   ```
4. **Hard-refresh the portal.** Panels are regenerated in place; a cached tab can show a figure
   that disagrees with the data.

---

## Beat 1 — WEIGHT (no screen yet)

> A 55-year-old man has chest pain that only turns up when he walks uphill. Three weeks of it. His
> father had a heart attack at fifty. His LDL is 185. His stress test came back equivocal — which
> is the worst answer available, because it means nobody knows anything yet.
>
> He is already on a statin. Atorvastatin, 20 milligrams. Someone was paying attention.
>
> So the question is not *does this man have a blockage.* We are going to find one. The question is
> whether he is about to repeat his father's fifties — and whether the drug he is already taking is
> the right drug at the right dose. Those are different questions, and only one of them is answered
> by a picture.

---

## Beat 2 — COMPRESSION (seven screens)

### 01 · Dashboard — the factory tells you what it is (ct-23, 60 s)

**Say:** One engine of eight, and everything I'm about to show you runs on the box under this
desk. Before running anything, note what it volunteers: 440 vectors across 13 collections, nine workflows, four model services — two labelled
*mock* in the product, where a customer would see it.

**Do:** Point at Total Vectors 440 and the collection rail · the NIM list (VISTA-3D *mock*, MAISI
*mock*, VILA-M3 *cloud*) · the "Not FDA-cleared" banner.

- `DECISION SUPPORT` — say it the first time any clinical surface appears.
- `LIVE` — counts read from Milvus at page load. If someone spots the three zeros in the rail: radiomics and reports fill at runtime, `genomic_evidence` is owned by another engine and read-only here. Ten collections hold the 440.
- `CAREFUL` — say "this workflow runs on one box," not "everything runs on one box." Heavy or ARM-incompatible models elsewhere in the factory burst to remote GPUs.

### 02 · Workflows — nine workflows, seven scoring systems (ct-24, 45 s)

**Say:** Seven of the nine report in the scoring system a radiologist actually uses — Lung-RADS,
BI-RADS, PI-RADS, LI-RADS, TI-RADS, Marshall CT, CAD-RADS. A system that invents its own severity
scale is one no clinician can sign.

**Do:** Sweep the scoring chips · scroll to Clinical Demo Cases · Run Case on the cardiac workup.

- `REPRESENTATIVE` — curated demo cases, not live patients.

### 03 · Panels A & B — a measurement, not an assertion (ct-25, 2 min)

**Say:** And there it is. Note the badge: it says *cached result*, not a latency — don't pretend a
precomputed case is a live inference. What matters is the number underneath. 77.9%. Not "severe" —
a number with a decimal, and beneath it the two lengths it came from: lumen 0.60 mm against a reference calibre of 2.70 mm. Geometric measurement off a segmented
vessel surface by shrinking-ball medial axis. Recomputable, identical every time.

**Do:** Read panel A's subtitle aloud (it names its own source) · point at panel B's caliper key.

- `LIVE` — 77.9% / 95.1% computed from the mesh.
- `NOT RUN` — coronary inference is not wired; SegResNet and VISTA-3D are labelled *not run*. We
  measure an expert's segmentation, we don't produce one.
- `REPRESENTATIVE` — **vessel names are scripted.** The geometry cannot resolve branch identity.
  Only the magnitude is measured. **Volunteer this before a cardiologist asks.**

### 04 · Panels C & D — the number in red isn't the one that kills him (ct-26, 2.5 min)

**Say:** CAD-RADS 4A/P3/HRP — three separate statements. *4A* is the narrowing (70–99%): the
plumbing, and why he hurts walking uphill. *P3* is plaque burden, 301–999 Agatston — severe. *HRP*
is high-risk plaque, and it's the one to lose sleep over.

HRP is not a vibe: CAD-RADS 2.0 requires **at least two** features in the same plaque, from low
attenuation, positive remodeling, spotty calcification, napkin-ring sign. The findings list
low-attenuation plaque and positive remodeling — two — **separately** from the percentage. Because
the lesion that kills is very often not the tight one: a soft, inflamed plaque under a thin cap
tears, clots in seconds, and can close an artery that was only modestly narrowed the day before. Stenosis tells you why he's symptomatic; plaque features tell you whether he's about to
have an event.

**Do:** Point at the badge, then the matching strip on every panel (one report string, reused) ·
land on the high-risk-plaque finding · panel D is the story in one frame.

- `LIVE` — grade derived deterministically from the measurement, same function the API uses.
- `REPRESENTATIVE` — Agatston 385 labelled *(repr.)*; a calcium score needs Hounsfield units a
  surface mesh cannot carry. Plaque features scripted.
- `REPRESENTATIVE` — panel D is generated art. Its numbers are pixels, so the regeneration script
  compares them against the analysis and fails loudly on drift.

### 05 · The gap, and the hand-off (ct-27, 3 min)

**Say:** The finding that matters most is not an imaging finding. He is on atorvastatin 20 mg —
*moderate* intensity, 30–49% LDL reduction. Established CAD at 4A calls for *high* intensity
(atorvastatin 40–80, rosuvastatin 20–40; ≥50%). Someone did the right thing; nobody escalated.

The system did not prescribe. It surfaced a gap, cited 2018 AHA/ACC, and wrote "decision support
only — not a prescribing instruction" into the finding itself.

Then the hand-off: Engine 4 is imaging, handing off to Engine 2, Precision Intelligence — the
evidence layer that owns the genomic knowledge base. Premature disease + family history +
LDL 185 is the familial-hypercholesterolaemia pattern, so the imaging result nominates the genes:
LDLR, PCSK9, APOB, LPA, 9p21. **The picture changed which genes were worth sequencing.**

**Do:** Read the therapy-gap finding aloud including its closing clause · walk the measurements
table (3 vessels diseased, **1** obstructive ≥70% — which is why this is 4A, not 4B) · point at the
cross-modal trigger band (**Engine 4 · Imaging → Engine 2 · Precision Intelligence**) · click **LDLR**.

**Why it matters (the question every cardiologist asks).** A statin does not clear the blockage; it
changes what the blockage does to you. Two jobs: lower LDL so less new plaque forms — the CTT meta-analysis of 26 trials puts that at a 22% relative reduction in major vascular events
per 1 mmol/L (~39 mg/dL) lowered — and stabilise existing plaque, thickening the cap and shrinking the lipid core (SATURN,
ASTEROID by IVUS; PARADIGM by serial CCTA). And blockages are rarely "removed": a stent props the
artery open, a bypass reroutes around it; the diseased wall and every unimaged plaque remain.
Revascularisation is a local fix, the drug is systemic — which is why COURAGE and ISCHEMIA found
stenting improved symptoms in stable disease without beating good medical therapy on hard outcomes.

- `DECISION SUPPORT` — repeat it here, with a drug and a dose on screen.
- `LIVE` — retrieval behind the gene queries is real; ~100 ms warm.

### 06 · LDLR, chromosome view — the gene the image recruited (ct-28, 1.5 min)

**Say:** LDLR, the receptor that clears LDL from the bloodstream. 19p13.2, ~11.09 Mb, marked on the
real ideogram at the real position. The card carries its provenance on its face — the transcript
it's numbered against, GRCh38 coordinates, the OMIM entry, and a ClinVar link so the room can check
us rather than trust us.

**Do:** Point at the marker position, then `*606945` and the phenotype entry · read the amber line
at the bottom aloud.

- `REPRESENTATIVE` — the card says it: *illustrative allele for this demo case — no sequencing was
  performed in this workflow.*

### 07 · LDLR, clinical details — a whole family in scope (ct-29, 2 min)

**Say:** If this is FH, everything reframes: lifetime LDL burden rather than midlife drift, which
means high-intensity statin from the start, likely ezetimibe added, PCSK9 inhibitor or inclisiran in
the conversation. The card gives that escalation in guideline order, because the guideline puts
ezetimibe first and payers generally expect that step.

Then the sentence that changes the arithmetic: FH affects roughly 1 in 250. Each first-degree
relative — his children, his siblings — has a 50% chance of carrying it. He came in with chest
pain; he may leave having put his kids into screening.

**Do:** Read the therapeutic implication, ending on "decision support for a clinician, not a
prescription" · point at Evidence & Provenance and the Verify in ClinVar link · land on 1-in-250
and 50%.

**Why it's a genomics question, not a pharmacy one.** Two genes, two jobs. LDLR is *cause* — it
explains the LDL and the early onset, and puts the family in scope. SLCO1B1 is *tolerability*: it
encodes the hepatic transporter OATP1B1, and c.521T>C carriers get higher systemic exposure and
materially higher muscle-toxicity risk — strongest with simvastatin, weaker with rosuvastatin or
pravastatin, with a CPIC guideline for exactly this. Pick the wrong statin for that genotype and
the patient gets myalgia, stops the drug, and loses the protection entirely.

- `DECISION SUPPORT` — genetic testing here is a recommendation to order, not a result we report.

---

## Beat 3 — HOPE

> What happened in the last ten minutes was not "we found a blockage and gave a pill."
>
> An imaging finding changed which genes were worth sequencing. The genotype changes which drug and
> which dose. And the diagnosis put an entire family into screening — people who are not patients
> yet and, if this works, never will be.
>
> That crossing is the whole point. It is also the one thing no single modality can do alone.
>
> Everything you just watched runs on one DGX Spark — 128 GB unified memory, about $4,699. The
> code is Apache-2.0:
> clone it tonight and read every line that produced every number you just watched, and the 1,365
> tests that hold it to the guideline.
>
> And the next patient gets the same afternoon.

---

## Honesty ledger — keep this answer in one breath

| On screen | Status | The honest sentence |
|---|---|---|
| 77.9% / 95.1% stenosis | **LIVE** | Computed from a segmented vessel surface, shrinking-ball medial axis. Reproducible. |
| CAD-RADS 4A/P3/HRP | **LIVE** | Derived deterministically from that measurement, same function the API calls. |
| Therapy-gap finding | **LIVE** | Rule over measurement + medication list, citing 2018 AHA/ACC. |
| Retrieval / 440 vectors | **LIVE** | Milvus, 13 collections (10 populated), ~100 ms warm. |
| Vessel identities (LAD/LCx/RCA) | REPR. | Scripted. Geometry cannot resolve branch identity. |
| Agatston 385 | REPR. | Needs Hounsfield units a surface mesh cannot carry. Labelled "(repr.)". |
| Plaque types & HRP features | REPR. | Scripted for this case. |
| LDLR exemplar variant | REPR. | No sequencing was performed in this workflow. |
| Panel D pathway figure | REPR. | Illustration — not this patient's anatomy, not to scale. |
| SegResNet, VISTA-3D | **NOT RUN** | Coronary inference not yet wired; labelled "not run" in the product. |

---

## What the room will ask

**"How did you get LAD out of an unlabelled mesh?"** We didn't. Magnitude measured, vessel name
scripted — the analysis file says so in its own caveats. Say it before they ask; it converts the
skeptic faster than anything else in the demo.

**"Does the statin clear the blockage?"** No, and that's the point. It lowers LDL so less new
plaque forms, and stabilises what's there — thicker cap, smaller lipid core. Less likely to rupture
even when it doesn't shrink much. *(CTT; SATURN/ASTEROID by IVUS; PARADIGM by serial CCTA.)*

**"If he's getting a stent, why does the drug matter?"** A stent is a local fix for a local flow
problem. The diseased wall and every unimaged plaque remain — and those are statistically more
likely to rupture. The statin is systemic. *(COURAGE, ISCHEMIA.)*

**"You never sequenced anything."** Correct, and the card says so. Imaging changed the pre-test
probability enough to make the genetic question worth asking. The output is a recommendation to
test plus what a result would mean. Wiring a real variant call in is the next step, not a claim.

**"Why is the calcium score representative?"** Agatston needs Hounsfield units; a surface mesh
carries geometry, not radiodensity. Same reason the percentile is "above the 90th" rather than
exact — a true MESA percentile also needs race/ethnicity, which this case doesn't specify.

**"Which statin, and why is that genomic?"** SLCO1B1 c.521T>C (rs4149056) → OATP1B1 reduced
function → higher exposure → higher myopathy risk, strongest with simvastatin. CPIC publishes
dosing guidance. The failure it prevents is mundane: muscle pain, patient stops the drug, all
protection lost. *(CPIC statins/SLCO1B1-ABCG2-CYP2C9 2022; Link et al., NEJM 2008.)*

**"Does it all really run on one box?"** This workflow does. Across the wider factory, heavy or
ARM-incompatible models burst to remote GPUs over a private mesh — say "elastic burst" and name it
rather than implying everything is local.

---

Clinical content adapted from the "Statin to Variant" note. Decision support for a qualified
clinician — not diagnosis, not prescribing, not FDA-cleared.
