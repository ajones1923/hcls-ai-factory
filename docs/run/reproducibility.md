---
title: Reproducibility & Datasets
description: Pinned, verifiable inputs — the reproducible build and the version-pinned public-datasets manifest behind every capability.
---

# Reproducibility & Datasets

A result you can't reproduce isn't evidence. The factory is built so its inputs are pinned and its
public datasets are versioned — so a `verified` badge means the same thing tomorrow as today.

## Reproducible by construction

- **Pinned build** — the docs site and the platform pin their generators; CI runs the real gate
  (lint · platform tests · registry validation) with no `|| true`, and the strict docs build fails
  before publishing anything stale.
- **Deterministic lineage** — governed runs carry a manifest (inputs → capability chain → serving →
  composed honesty tier → a deterministic lineage hash), so a result traces to exactly what made it.
  See [Governance & Lineage](../honesty/governance.md).

## Public-datasets manifest (version-pinned)

Every public dataset used across the engines and agents is pinned to a release, so cross-omics joins
stay reproducible over time. Data itself is never committed — only its provenance.

| Dataset | Version / release | Used by |
|---|---|---|
| GRCh38 primary assembly | Ensembl release-112 | foundation / genomics |
| HGNC complete set · HPO · Cell Ontology · Reactome | current releases | BioKey resolver, phenotype/pathway |
| ClinVar (GRCh38) | weekly VCF | E1 annotation · **ACMG secondary-findings (verified)** |
| AlphaMissense hg38 | v1 (2023) | missense pathogenicity |
| ChEMBL | chembl_37 (SQLite) | therapeutic discovery reference |
| CPIC (diplotype/recommendation/allele tables) | API snapshot | pharmacogenomics |
| Orphanet (Orphadata product 1/4/6) | current `en_` | rare disease |
| CIViC (evidence/variant/gene) | nightly | precision oncology |
| 10x pbmc3k | scanpy example | **single-cell compute (verified)** |

*The full manifest and a resumable downloader live in the neutral repo's central public-datasets
store; only code and docs publish — data, weights, and secrets stay local.*

!!! note
    Version pins are the difference between "it worked once" and "it reproduces." Where a capability is
    marked [`verified`](../honesty/ledger.md), the exact dataset + version it was proven against is
    recorded on [Citations & Evidence](../honesty/citations.md).
