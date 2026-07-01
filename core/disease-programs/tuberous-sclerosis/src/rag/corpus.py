"""
Seed TSC literature corpus (PRD §2.5.4; master paper §16). A small, source-attributed
set of the literature anchors, partitioned by source. The full PubMed/PMC, ClinicalTrials.gov,
and FDA ingestion at scale is the W5 upgrade; every chunk carries the source_uri that
flows into the Therapeutics Strategist's attribution.
"""
from __future__ import annotations

SEED_CORPUS = [
    {"id": "itsc2021-surveillance", "partition": "itsc_guidelines", "pub_year": 2021,
     "section": "surveillance",
     "source_uri": "https://doi.org/10.1016/j.pediatrneurol.2021.07.011",
     "text": "2021 International TSC Consensus: surveillance recommends brain MRI every 1-3 years "
             "for SEGA monitoring, renal imaging every 1-3 years for angiomyolipoma, and annual "
             "screening for TSC-Associated Neuropsychiatric Disorders (TAND)."},
    {"id": "exist3-everolimus", "partition": "pubmed_pmc", "pub_year": 2016,
     "section": "therapeutics",
     "source_uri": "https://doi.org/10.1016/S0140-6736(16)31419-2",
     "text": "EXIST-3 demonstrated adjunctive everolimus significantly reduced seizure frequency in "
             "treatment-resistant focal seizures associated with TSC, the first targeted therapy in a "
             "genetically defined epilepsy; mTOR inhibitors also reduce SEGA and renal AML volume."},
    {"id": "tosca-tand", "partition": "pubmed_pmc", "pub_year": 2018,
     "section": "tand",
     "source_uri": "https://doi.org/10.1016/j.ejpn.2018.05.006",
     "text": "TOSCA natural-history registry: TAND features are present in ~90% of TSC patients but are "
             "formally evaluated and addressed in only a minority; 30-50% of features go unrecognized "
             "in routine care."},
    {"id": "mosaicism-tyburczy2015", "partition": "pubmed_pmc", "pub_year": 2015,
     "section": "genetics",
     "source_uri": "https://doi.org/10.1371/journal.pgen.1005637",
     "text": "Deep sequencing recovers low-level mosaic TSC1/TSC2 variants in patients meeting clinical "
             "criteria with negative conventional testing, explaining a large share of the no-mutation-"
             "identified cohort; affected tissue and high read depth are required."},
    {"id": "mtorc1-next-gen", "partition": "fda_actions", "pub_year": 2024,
     "section": "emerging",
     "source_uri": "https://clinicaltrials.gov/",
     "text": "Next-generation selective mTORC1 inhibitors are in clinical development for TSC, aiming to "
             "improve the therapeutic index relative to rapalogs; trial eligibility commonly requires a "
             "confirmed TSC1/TSC2 molecular diagnosis and a stable anti-seizure regimen."},
    {"id": "tand-consensus", "partition": "itsc_guidelines", "pub_year": 2018,
     "section": "tand",
     "source_uri": "https://doi.org/10.1016/j.ejpn.2018.05.014",
     "text": "The TAND consensus framework defines six clusters (behavioral, psychiatric, intellectual, "
             "academic, neuropsychological, psychosocial) and the TAND-L lifetime checklist for "
             "systematic screening."},
]
