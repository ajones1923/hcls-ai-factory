## In plain terms

The Cardiology Engine is a cardiovascular **risk-and-prevention** engine: **11 clinical workflows** and
**6 risk calculators** spanning prevention, intervention, and rhythm. Its signature move is precision
*prevention* — combine a coronary calcium score with the patient's genetics, then hand off to
pharmacogenomics to pick the statin that is actually safe for them.

## Why it matters

The highest-leverage cardiology happens before the heart attack. A calcium score plus a ten-year risk
estimate tells you *whether* to treat; the patient's pharmacogenomics tells you *which* drug won't harm
them. Doing both, together, is precision prevention rather than one-size-fits-all.

*For a patient: prevention matched to their heart and their genes — the right risk assessment, then a drug that is safe for them specifically.*

## How it works

![Inside the Cardiology Engine — imaging and genotype in, risk scored, genetically safe drug out](../../assets/infographics/pages/cardiology-intelligence-agent-how.png)
/// caption
Coronary calcium and 10-year ASCVD risk, handed to pharmacogenomics for the safe drug. Illustrative.
///

1. **Ingest** — imaging, genotype, and labs.
2. **Score** — coronary artery calcium (**CAC**) and **10-year ASCVD** risk, using the **ACC/AHA**
   prevention logic.
3. **Recommend** — across 11 workflows and 6 calculators, from prevention through intervention to
   rhythm.
4. **Hand off** — to the [Pharmacogenomics Intelligence Agent](../agents/pharmacogenomics-intelligence-agent.md)
   for the genetically safe drug and dose.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (imaging, genotype, labs).
- **Out:** **risk scores** and prevention recommendations.

## Where it fits

![Where the Cardiology Engine sits — joining imaging and pharmacogenomics](../../assets/infographics/pages/cardiology-intelligence-agent-fits.png)
/// caption
Calcium score to genetically safe statin: imaging in, pharmacogenomics out. Illustrative.
///

It joins the [Clinical Imaging Engine](imaging-intelligence-agent.md) (for the calcium score) and the
Pharmacogenomics agent (for the safe drug) — the spine of the cross-modality prevention demonstration.

## Honest limits

- **Guideline-cited.** The calculators implement published **ACC/AHA** logic and cite it.
- **Decision support, not prescribing.** Risk scores and recommendations support a clinician's
  decision; the engine never prescribes.
