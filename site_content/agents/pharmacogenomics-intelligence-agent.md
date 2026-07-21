## In plain terms

The Pharmacogenomics Intelligence Agent answers one deceptively hard question: **given this patient's
genes, which drug and what dose is safe?** It reads a patient's pharmacogenes, works out how fast they
metabolize a drug, and grounds a dosing recommendation in the published **CPIC** guidelines — as
decision support for the prescriber, never an autonomous prescription.

## Why it matters

The same standard dose can be ineffective in one person and toxic in another, purely because of
inherited differences in drug metabolism. Adverse drug reactions are a major cause of preventable
harm. Getting the drug-gene match right *before* the first dose is one of the most concrete wins in
precision medicine.

*For a patient: the right drug at the right dose from the first prescription — matched to how their body actually processes it.*

## How it works

![How the Pharmacogenomics agent reasons — genotype to star allele to phenotype to CPIC dose](../../assets/infographics/pages/pharmacogenomics-intelligence-agent-how.png)
/// caption
Genotype → metabolizer phenotype → guideline dose. Decision support, never autonomous prescribing. Illustrative.
///

1. **Call star alleles** — the patient's drug-metabolism gene (e.g. *CYP2D6*, *TPMT*) is resolved into
   its named versions — *star alleles*, the standard shorthand (like \*2 or \*3) for the variants that
   change how fast that gene's enzyme works.
2. **Assign phenotype** — those alleles map to a **metabolizer phenotype** (poor → ultra-rapid).
3. **Apply the guideline** — the phenotype is matched to **CPIC** dosing guidance, with drug-metabolism
   filtering and a safety interlock.
4. **Ground the answer** — retrieval-augmented generation returns a cited recommendation, refusing to
   fabricate when the evidence isn't there.

## What goes in, what comes out

- **In:** a clinical **query** and the **patient context** (genotype).
- **Out:** a grounded, cited dosing **answer** for the clinician.

## Where it fits

![Where the Pharmacogenomics agent sits — the safety interlock before any prescription](../../assets/infographics/pages/pharmacogenomics-intelligence-agent-fits.png)
/// caption
The right-drug-right-dose safety checkpoint the other engines hand off to. Illustrative.
///

It is the **safety interlock**: the [Cardiology Engine](../engines/cardiology-intelligence-agent.md)
and oncology hand off to it for the genetically safe drug, and it can flag any prescribing decision
against the patient's metabolism — the prescriber still decides.

## Honest limits

- **Decision support, never prescribing.** It informs a qualified prescriber; it does not prescribe.
- **Grounded, and honest when it can't be.** Like every intelligence agent it is a retrieval-augmented
  service: it needs a populated vector database and an LLM API key at runtime, and returns an honest
  degraded response (HTTP 503) rather than inventing clinical content when they're absent.
- **Guideline-cited.** Recommendations trace to **CPIC**; where guidance is absent, it says so.
