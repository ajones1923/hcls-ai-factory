## In plain terms

The Genomics Foundations Engine is the factory's **front door**. A DNA sequencer never hands you a
person's genome as a tidy readout — it hands you hundreds of millions of short, overlapping fragments
of text called *reads*. This engine turns that raw pile into a clean, trustworthy list of exactly
where a patient's DNA differs from the standard reference human genome. Those differences — called
**variants** — are the raw material that every other engine and agent in the factory reasons about.

It does the job the way a modern clinical genomics lab does, only faster: it runs **GPU-accelerated**
tools — NVIDIA Parabricks for read alignment and Google's **DeepVariant** for variant calling — so a
whole genome goes from raw reads to a finished variant file in **hours instead of the day or two** a
CPU pipeline would take.

## Why it matters

Everything downstream — interpreting a variant's meaning, matching a patient to a therapy, designing
a molecule, running a disease program — is only ever as good as this first step. A missed variant here
is a missed diagnosis later; a false one is a false alarm. The engine's job is therefore not just
*speed* but **correctness you can check**: it produces a result whose quality can be measured and
reproduced, not merely asserted.

## How it works

![Inside the Genomics Foundations Engine — align, call, quality-check, hand off](../../assets/infographics/pages/genomics-engine-how.png)
/// caption
Germline variant calling: FASTQ → BAM → VCF, in hours not days. Illustrative.
///

1. **Align** — the millions of raw reads are mapped onto the GRCh38 reference genome (Parabricks
   `fq2bam`, built on BWA-MEM2), so every fragment is placed where it belongs.
2. **Call variants** — **DeepVariant**, a convolutional neural network, reads the stacked-up
   ("pileup") evidence at each position and decides where the patient genuinely differs from the
   reference, producing germline SNVs and small insertions/deletions.
3. **Quality-check** — the run is sanity-checked before anyone trusts it. A headline metric is the
   **transition/transversion (Ts/Tv) ratio**, which lands near **2.0** for a healthy whole-genome
   call set — a quick, standard signal that calling worked.
4. **Hand off** — the QC-passed variant file (VCF) is handed to the
   [Precision Intelligence Engine](precision-intelligence-engine.md), which annotates and interprets
   it. This engine finds *where* the variants are; the next engine explains *what they mean*.

## What goes in, what comes out

- **In:** paired-end sequencing reads (**FASTQ**) and the **GRCh38** reference genome.
- **Out:** a **VCF** of called variants and the **BAM** of aligned reads — the substrate the rest of
  the factory builds on.

## Where it fits

![The front door of the pipeline — genomics feeds interpretation, agents, and disease programs](../../assets/infographics/pages/genomics-engine-fits.png)
/// caption
The genomics substrate flows to interpretation, then to the agents and disease programs. Illustrative.
///

The engine sits at the head of the pipeline. It does **not** annotate or interpret variants itself —
that is deliberately the [Precision Intelligence Engine](precision-intelligence-engine.md)'s job
(ClinVar / AlphaMissense annotation on `:5001`). Keeping *calling* and *interpretation* separate is
what lets each stay honest and independently verifiable.

## Honest limits

- **Germline, not tumor.** This engine calls inherited (germline) variants. Somatic/tumor analysis
  lives in the [Precision Oncology Engine](precision-oncology-agent.md).
- **Script-invoked today.** It currently runs as a container pipeline (Parabricks in Docker via
  `run.sh` / `run_caller.sh`), **not yet a request/response API** — that API is on the roadmap.
- **Elastic burst.** Parabricks is an x86-only CUDA container, so on an ARM DGX Spark this step
  **bursts to a remote GPU** over a private mesh; only non-identifying work leaves the box.
- **Proven on a public benchmark.** The companion **Variant Store** capability is marked `verified`
  against the public **HG002** reference sample (Ts/Tv ≈ 2.0) — real-data evidence, not a claim.
- **Decision support, not diagnosis.** The output supports a qualified clinician's judgment; it never
  diagnoses on its own.
