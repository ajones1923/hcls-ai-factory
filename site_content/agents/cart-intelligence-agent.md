## In plain terms

The CAR-T Intelligence Agent is a research companion for **cell-therapy design**. It gathers and
compares the evidence behind a CAR-T construct — the antigen it targets, the binder (scFv/antibody)
sequences that recognize it, and what the literature says has worked — and hands a clinician or
scientist a grounded evidence brief. It reasons; it does not design a therapy on its own.

## Why it matters

CAR-T can be curative, but designing and choosing one is an evidence-heavy, high-stakes exercise
spread across scattered literature and sequence databases. Pulling that evidence into one comparative,
cited view is where an intelligence layer earns its place.

## How it works

![How the CAR-T agent reasons — target antigen to comparative evidence to design brief](../../assets/infographics/pages/cart-intelligence-agent-how.png)
/// caption
Cross-collection evidence and comparative analysis for CAR-T design. Decision support. Illustrative.
///

1. **Frame the target** — the target antigen and patient context.
2. **Retrieve across collections** — cross-collection evidence, including stored **scFv/antibody**
   sequences.
3. **Compare** — comparative analysis and deeper research across candidate constructs.
4. **Ground the answer** — a cited evidence brief; it refuses to fabricate where evidence is thin.

## What goes in, what comes out

- **In:** a **query** and the **patient context**.
- **Out:** a grounded, cited CAR-T design evidence **brief**.

## Where it fits

![Where the CAR-T agent sits — reasoning with structural biology and single-cell](../../assets/infographics/pages/cart-intelligence-agent-fits.png)
/// caption
It draws on structural biology and single-cell context; anchor of the CAR-T demonstration. Illustrative.
///

It reasons alongside the structural-biology and single-cell capabilities and anchors the
oncology / CAR-T design demonstration.

## Honest limits

- **Decision support, never diagnosis.** It supports a qualified clinician or scientist; it does not
  make treatment decisions.
- **Grounded, and honest when it can't be.** As a retrieval-augmented service it needs a populated
  vector database and an LLM API key at runtime, returning an honest degraded response (HTTP 503)
  rather than inventing content when they're absent.
- **Design frontier is separate.** De-novo binder design (Chai-2) is a distinct, gated frontier
  capability — this agent reasons over evidence, it does not generate binders.
