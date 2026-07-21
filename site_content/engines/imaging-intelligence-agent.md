## In plain terms

The Clinical Imaging Engine reads medical images and reasons about them **alongside** the patient's
genomics. It turns a scan — delivered as standard **DICOM** — into structured findings a clinician can
act on, and it can cross-reference those findings against what the genome says. Image and gene
together often say more than either alone.

## Why it matters

A radiology read and a genomic finding usually live in separate systems and separate heads. Bringing
them into one reasoning step — a calcium score next to a lipid genotype, a tumor's imaging next to its
fusion — is where cross-modal precision medicine actually happens.

*For a patient: their scan and their genetics read together, so nothing that matters falls through the gap between two systems.*

## How it works

![Inside the Clinical Imaging Engine — ingest DICOM, analyze pixels, reason cross-modally, export FHIR](../../assets/infographics/pages/imaging-intelligence-agent-how.png)
/// caption
DICOM in, structured findings out — reasoned against genomics. Illustrative.
///

1. **Ingest** — a standard **DICOM** imaging study comes in (referenced through the patient context).
2. **Analyze** — the live, verified modality today is **chest X-ray analysis via a DenseNet-121 model
   on real DICOM pixels**.
3. **Reason** — that live read is joined **cross-modally** with the patient's genomics (the deeper
   image-and-report reasoning, VILA-M3 on CT, is on the planned path — see *Honest limits*).
4. **Export** — findings are written out as **FHIR R4** so they slot into clinical systems.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (which references the DICOM study).
- **Out:** a structured **image read**.

## Where it fits

![Where the Clinical Imaging Engine sits — feeding oncology, cardiology, and cross-modal reasoning](../../assets/infographics/pages/imaging-intelligence-agent-fits.png)
/// caption
Imaging joins the tumor board, the calcium-score-to-statin path, and genomic reasoning. Illustrative.
///

Its reads feed the [Precision Oncology Engine](precision-oncology-agent.md) (imaging joins the tumor
board) and the [Cardiology Engine](cardiology-intelligence-agent.md) (a coronary-calcium input — a
CT-derived value, on the planned CT path, not produced by the live chest-X-ray model), and they
reason cross-modally with genomics throughout.

## Honest limits

- **What's live vs. planned.** The **chest X-ray / DenseNet-121** path is the live-verified modality.
  The **VISTA-3D segmentation and VILA-M3 CT** path is **planned / placeholder — not yet live**.
- **Synthetic imaging is not diagnostic.** **MAISI** synthetic image generation is for
  research / augmentation / QA only — **never a diagnostic source**.
- **Decision support, not diagnosis.** The read supports a qualified clinician; it does not diagnose
  on its own.
