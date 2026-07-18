## In plain terms

The Single-Cell Intelligence Agent is the **reasoning layer** over the
[Single-Cell Analysis Engine](../engines/singlecell-compute.md). The engine computes cell-type
clusters; this agent **interprets** them — annotating cell types, profiling the tumor
microenvironment, and reasoning about drug response — grounded in the literature. The engine computes;
the agent interprets. They are two roles, not a duplicate.

## Why it matters

A clustered single-cell dataset is a starting point, not an answer. The clinically interesting
questions — *what* are these cells, is this tumor immunologically hot or cold, which subclone is
driving it, what might it respond to — are interpretation problems. Separating the deterministic
compute from the grounded reasoning is what keeps both honest.

## How it works

![How the Single-Cell agent reasons — from computed clusters to interpretation](../../assets/infographics/pages/single-cell-intelligence-agent-how.png)
/// caption
Interpreting the engine's clusters: annotation, TME profiling, drug-response hypotheses. Illustrative.
///

1. **Take the clusters** — the cell-type clusters computed by the engine, plus patient context.
2. **Annotate** — cell-type annotation grounded in markers and literature.
3. **Profile the TME** — tumor-microenvironment profiling (hot / cold / excluded), subclonal
   architecture, and ligand-receptor analysis.
4. **Reason about response** — drug-response reasoning, returned as a cited answer that refuses to
   fabricate — and whose predictions are hypotheses to test.

## What goes in, what comes out

- **In:** a **query** and the **patient context** (the engine's computed clusters).
- **Out:** a grounded, cited single-cell **interpretation**.

## Where it fits

![Where the Single-Cell agent sits — the reasoning cap on the compute engine](../../assets/infographics/pages/single-cell-intelligence-agent-fits.png)
/// caption
The reasoning cap on the compute base — the engine computes, the agent interprets. Illustrative.
///

It sits directly on top of the Single-Cell Analysis Engine — the engine produces the clusters, the
agent reasons over them — and feeds oncology's tumor-microenvironment view.

## Honest limits

- **Not a duplicate of the engine.** The engine computes clusters (verified on PBMC 3k); this agent
  interprets them. Two distinct roles.
- **Drug-response predictions are hypotheses.** They are leads to test, never treatment decisions.
- **Decision support, never diagnosis.** As a retrieval-augmented service it needs a populated vector
  database and an LLM API key at runtime, returning an honest degraded response (HTTP 503) rather than
  inventing content when they're absent.
