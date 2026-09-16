# Data licences

**The code in this repository is Apache-2.0. The data it reads is not.**

This file exists because those two facts are easy to conflate, and the conflation runs in the
dangerous direction: a reader sees "Apache-2.0, commercial use welcome" on the project and assumes
it covers everything the platform touches. It does not. The platform is a set of pipelines that
*read* public biomedical datasets, each carrying its own terms, several of which are
**non-commercial**.

## The governing fact

**No third-party dataset is redistributed by this repository.** Verified, and you can verify it:

```bash
git ls-files | grep -iE 'alphamissense|clinvar|oncokb|cosmic'   # returns nothing
```

Data, model weights and credentials stay local and are never committed. Everything below is
something *you* download from its source, under that source's terms, onto your own machine. This
repository ships the code that reads it and the instructions for obtaining it.

That is why there is no licence conflict *in the repository*. The obligation this file discharges
is the other one: telling you what you are about to pull, before you pull it.

## The one that most often surprises people

> **AlphaMissense is CC BY-NC-SA 4.0 — non-commercial and share-alike.**
> It is the most-referenced external artefact in this codebase. If your deployment is commercial,
> you need permission from the rights holder, or you run without it.

The platform degrades cleanly without it, by construction rather than by luck: if the file is
absent, `core/engines/precision-intelligence/scripts/ingest_vcf.py` checks `.exists()`, logs
*"AlphaMissense annotation will be skipped"*, and continues on ClinVar significance alone. The
score is `Optional` everywhere it is carried. It is a missing column, not a broken pipeline.

Two further datasets are gated in practice rather than merely attributed:

| Dataset | Practical restriction |
|---|---|
| **OncoKB** | Requires a **licensed API key**; commercial use requires a licence. The client refuses to run without one rather than silently degrading — see `core/engines/precision-oncology/agent/src/ingest/oncokb_parser.py`. |
| **COSMIC** | Free for academic use; **commercial use requires a paid licence.** |
| **OMIM** | Bulk download requires registration and a licence agreement. |

## Every external dataset this platform reads

The list is complete as to *which* datasets are referenced. The **Terms** column links to the
authority — check it there rather than trusting a table in a repository, including this one.

| Dataset | Used for | Terms |
|---|---|---|
| ClinVar | Variant clinical significance | <https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/> |
| AlphaMissense | Missense pathogenicity prediction | **CC BY-NC-SA 4.0** — <https://zenodo.org/records/8360242> |
| dbSNP | Variant identifiers | <https://www.ncbi.nlm.nih.gov/snp/docs/> |
| gnomAD | Population allele frequencies | <https://gnomad.broadinstitute.org/policies> |
| OncoKB | Oncology variant actionability | **Licence required** — <https://www.oncokb.org/apiAccess> |
| COSMIC | Somatic mutation catalogue | **Commercial licence required** — <https://cancer.sanger.ac.uk/cosmic/licensing> |
| CIViC | Clinical interpretation of variants | <https://civicdb.org/about> |
| PharmGKB | Pharmacogenomic annotations | <https://www.pharmgkb.org/page/dataUsagePolicy> |
| CPIC | Pharmacogenomic dosing guidelines | <https://cpicpgx.org/> |
| OMIM | Mendelian disease catalogue | **Registration + licence** — <https://www.omim.org/help/agreement> |
| Orphanet | Rare disease nomenclature | <https://www.orphadata.com/> |
| ClinicalTrials.gov | Trial registrations and eligibility | <https://clinicaltrials.gov/about-site/terms-conditions> |
| PubMed / PMC | Literature abstracts and open-access text | <https://www.ncbi.nlm.nih.gov/home/about/policies/> · full text varies **per publisher** |
| GRCh38 reference | Alignment reference genome | <https://www.ncbi.nlm.nih.gov/grc> |
| HG002 / GIAB | Benchmark sample FASTQ | <https://www.nist.gov/programs-projects/genome-bottle> |
| Protein Data Bank | Structures for docking and visualisation | <https://www.rcsb.org/pages/usage-policy> |
| MedMNIST | Imaging demo data | CC BY 4.0 — <https://medmnist.com/> |

Model weights and containers carry separate terms and are inventoried in
[`docs/build/ACQUISITION_MANIFEST.md`](docs/build/ACQUISITION_MANIFEST.md) (NGC-gated items,
HuggingFace licence acceptance, Parabricks). The Imaging engine keeps its own per-technology
inventory in `core/engines/clinical-imaging/agent/LICENSES.md`.

## If you are deploying commercially

1. Decide about **AlphaMissense** first — it is the widest dependency and the clearest restriction.
2. **OncoKB** and **COSMIC** need licences; neither is required for the platform to run.
3. Re-read the **Terms** column above for anything else you intend to ingest. Terms change.

*This is a pointer to each dataset's own terms, maintained in good faith. It is not legal advice
and it is not a substitute for reading them.*
