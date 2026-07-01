# What the TSC Intelligence Engine Gives Cincinnati Children's

### A head start in pediatric rare-disease intelligence — built with Cincinnati Children's in mind, and given freely.

**From:** Adam M. Jones — Architect, HCLS AI Factory; Product Marketing Manager, HPC, the platform
**Prepared for:** Dr. Philip A. Hagedorn (CHIO) and the Cincinnati Children's TSC team
**Platform:** HCLS AI Factory, **Engine 7** — the TSC Intelligence Engine · **License:** Apache 2.0 (free, open, non-commercial) · 2026

---

> **It is already running.** As of June 2026 a working v0.1 is built and verified on my NVIDIA DGX Spark — all five agents, the orchestrator, a synthetic 50-patient cohort, and three clinician surfaces, with 41 automated tests passing and the three-act demonstration running end-to-end. The frontier-model reasoning is verified live (the Variant Curator produced a genuine Claude-Opus ACMG narrative for the featured mosaic patient). All data is synthetic and watermarked; every output is decision support for clinician review; nothing is FDA-cleared. This is not a slide deck — it is something your team can watch run.

## The offer, plainly

For more than a decade I have been building an open-source precision-medicine platform, the HCLS AI Factory. After our conversation about Tuberous Sclerosis Complex, the Discover Together Biobank, and the Winslow Research Pavilion, I built a TSC-specific engine within it. It is Apache 2.0 — free, open, and not something I will ever commercialize. I built it with Cincinnati Children's in mind, and I would like to give it to you in the way that does the most good: as a head start for your TSC program.

This page is the short version. A full research volume and an engineering build specification sit behind it, but the case is simple enough to state here.

## What it is

The TSC Intelligence Engine is five coordinated AI agents plus a deterministic orchestrator that take the evidence that is already in a TSC patient's record — genomics, clinical notes, imaging, longitudinal measurements — and turn it into reviewable, clinician-facing surfaces: a one-screen pre-visit briefing, an in-visit dashboard, and a disciplined set of alerts. The agents recover low-frequency mosaic variants that standard blood testing misses, map longitudinal HPO phenotypes and flag surveillance gaps, forecast SEGA and AML trajectories, surface under-recognized TAND signals, and assemble source-attributed therapeutic options briefs. Every output carries full provenance and is decision support behind a clinician's review — nothing is autonomous or diagnostic. It runs on a single desktop-class NVIDIA system.

## What it gives Cincinnati Children's

- **Reference-deployment leadership.** Cincinnati Children's would be the first institution to operationalize a disease-specific intelligence engine in pediatric rare disease. When other centers adopt the pattern, they replicate *from you*. This is leadership, not adoption.
- **A productive Winslow Pavilion and biobank.** You built the building and the freezer archive; this is the computational layer that turns banked tissue into discovery. *We built the biobank and the intelligence to interrogate it* is a stronger grant, donor, and U.S. News story than either investment alone.
- **An extension of your own program, and publications.** The TAND Surveillance Agent is a direct downstream application of Dr. Hagedorn's published work on diagnostic uncertainty in clinical documentation — an extension of an existing Cincinnati Children's research program, not an outside system imposed on it. The first somatic-mosaic-recovery analysis and the TAND-surveillance work are publications waiting to happen, with co-authorship belonging here.
- **Zero acquisition friction.** Apache 2.0, no licensing fees, no vendor lock-in, and small enough to run on a desktop appliance. Cincinnati Children's can evaluate, adopt, extend, and *own* it without a procurement gauntlet — a speed almost nothing in enterprise clinical AI allows.

## What it draws on here — and what it does not touch

The Engine integrates with, rather than contains, Cincinnati Children's. The Winslow Pavilion is the physical envelope; five areas feed the Engine as sources: the Discover Together Biobank (banked tuber, AML, and SEGA tissue — the molecular substrate that lets the Variant Curator recover diagnoses blood testing misses), the Division of Biomedical Informatics under Dr. Hagedorn (the methodology and the disciplined surfacing of output as briefings, not alert fatigue), the TSC clinical program (the patient concentration that makes the work meaningful), and Epic and the biobank LIMS (the data plumbing). A biobank without an intelligence layer is a freezer full of tubes; this is the layer.

## What is real today, and the honest line

Today the Engine runs as a complete, working demonstration on a **synthetic**, watermarked patient cohort — no protected health information, and nothing that requires an IRB to see. Real-data validation, Epic and biobank-LIMS integration, and institutional governance are the next phase, and I have scoped them candidly rather than glossed them. The Engine is not FDA-cleared, its outputs are reviewed by qualified clinicians, and it is designed to make clinicians more effective, never to replace them. I would rather show you something honest about its edges than oversell it — that honesty is the point, and it is what makes the rest credible.

## The ask

Only a conversation. I would value the chance to show the Engine to you and to the colleague leading your TSC work — a short demonstration and an open discussion about what is most clinically useful, what to prioritize, and what you would want done differently. There is no commitment and no cost. I would rather shape this with your input than work in isolation and hope it lands.

If it helps even a handful of children get an earlier diagnosis or better-monitored care, that is where my heart has been since I started this work. I would be honored to give it a home at Cincinnati Children's.

— **Adam**
