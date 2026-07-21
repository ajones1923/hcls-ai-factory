## In plain terms

Proteins are the machines of biology — and the targets and the *modality* of large-molecule medicines
like antibodies and cell therapies. The Large-Molecule / Structural Biology Engine is the factory's
protein workbench: it predicts a protein's 3-D **structure** from its sequence, **searches** for
related proteins, **designs** new sequences, and **assesses** whether a candidate is manufacturable
and safe to the immune system.

## Why it matters

Reading a protein's structure and reasoning about how to change it used to take a lab weeks. Doing it
in seconds, on real models, is what makes computational biologics and binder design possible — and it
is the substrate the therapeutic-discovery and cell-therapy work builds on.

*For a patient: the protein groundwork behind the next generation of biologics and cell therapies designed for their disease.*

## How it works

![Inside the Structural Biology Engine — fold, search, design, assess](../../assets/infographics/pages/structural-biology-engine-how.png)
/// caption
Sequence → structure → design → developability, on real GPU-served models. Illustrative.
///

1. **Fold** — **ESMFold** predicts a single protein's 3-D structure from its sequence, with a
   per-residue confidence score (pLDDT). *(verified on real sequences.)*
2. **Search** — **ESM-2** embeddings plus Smith-Waterman alignment find related proteins by sequence
   and by learned similarity. *(verified.)*
3. **Design** — **ProteinMPNN** designs new amino-acid sequences that should fold to a desired
   backbone. *(verified.)*
4. **Assess** — developability and MHC-immunogenicity checks flag whether a candidate is manufacturable
   and unlikely to provoke an immune response.

## What goes in, what comes out

- **In:** a **protein sequence** (or a backbone to design onto).
- **Out:** a predicted **structure**, similarity **search** hits, **designed sequences**, and
  developability / immunogenicity **scores**.

## Where it fits

![Where the Structural Biology Engine sits — feeding therapeutic discovery and cell-therapy design](../../assets/infographics/pages/structural-biology-engine-fits.png)
/// caption
The protein substrate under structure-based drug design and biologic/binder work. Illustrative.
///

Its structures and designs feed the [Therapeutic Discovery Engine](therapeutic-discovery-engine.md)
(structure-based design) and the [CAR-T Intelligence Agent](../agents/cart-intelligence-agent.md)
(binder evidence), and support the oncology, autoimmune, and rare-disease work.

## Honest limits

- **What's verified vs. planned.** **ESMFold**, **ESM-2 search**, and **ProteinMPNN** are `verified`
  on real data. MHC **immunogenicity** (MHCflurry), **ESM-2 fine-tune**, and **protein developability**
  are `planned` — their dependencies aren't yet installed in the runtime.
- **Frontier co-folding is separate.** AlphaFold3-class **complex** co-folding (**Chai-1**) — which
  ESMFold's single-sequence folding can't do — is `planned` and **bursts to a remote GPU**;
  de-novo binder design (**Chai-2**) is `gated`.
- **Research-use / decision support.** Structures and designs are research and design outputs for a
  qualified scientist — not a therapy, and not a clinical diagnosis.
