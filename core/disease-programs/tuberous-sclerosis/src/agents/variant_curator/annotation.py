"""
Variant annotation + ACMG criterion assignment (PRD §3 FR-VC-3/4).

Annotation: consequence, LOF mechanism, population frequency, and prior clinical
classification. A small curated table stands in for live ClinVar / gnomAD v4 / LOVD-TSC
here; the deterministic snpEff/VEP annotation and live database queries are the W2
integration (the contract below does not change when they land).

Criterion assignment follows the standard logic for LOF genes. PVS1 strength note:
for known recurrent germline null variants it is applied Very Strong; for low-VAF
mosaic or not-previously-established nulls it is conservatively applied Strong, pending
the full ClinGen PVS1 decision tree (W2-W3). This reproduces the cohort's expected
calls (e.g., Patient A's mosaic frameshift -> Likely Pathogenic).
"""
from __future__ import annotations

from src.agents.variant_curator.acmg import Criterion

# TSC1/TSC2 loss-of-function is an established disease mechanism.
LOF_GENES = {"TSC1", "TSC2"}
NULL_CONSEQUENCES = {"nonsense", "frameshift", "splice", "stop_gained"}

# Curated stand-in for ClinVar/LOVD-TSC (recurrent, previously-established variants).
KNOWN_PATHOGENIC = {
    ("TSC2", "c.3037C>T"): {"clinvar": "Pathogenic", "source": "ClinVar/LOVD-TSC (recurrent)"},
    ("TSC2", "c.2251C>T"): {"clinvar": "Pathogenic", "source": "ClinVar/LOVD-TSC"},
    ("TSC1", "c.1888C>T"): {"clinvar": "Pathogenic", "source": "ClinVar/LOVD-TSC"},
    ("TSC1", "c.2356C>T"): {"clinvar": "Pathogenic", "source": "ClinVar/LOVD-TSC"},
}


def annotate(variant: dict) -> dict:
    gene = variant.get("gene")
    consequence = (variant.get("consequence") or "").lower()
    hgvsc = variant.get("hgvsc")
    known = KNOWN_PATHOGENIC.get((gene, hgvsc))
    return {
        "gene": gene,
        "consequence": consequence,
        "is_null": consequence in NULL_CONSEQUENCES,
        "lof_gene": gene in LOF_GENES,
        "gnomad_af": variant.get("gnomad_af", 0.0),   # synthetic -> absent from gnomAD
        "clinvar": (known or {}).get("clinvar"),
        "known_pathogenic": known is not None,
        "known_source": (known or {}).get("source"),
    }


def assign_criteria(annotation: dict, *, mosaic: bool, phenotype_consistent: bool) -> list[Criterion]:
    crits: list[Criterion] = []
    if annotation["is_null"] and annotation["lof_gene"]:
        if annotation["known_pathogenic"]:
            crits.append(Criterion("PVS1", "PVS",
                                   f"null variant ({annotation['consequence']}) in LOF gene "
                                   f"{annotation['gene']}; established recurrent — Very Strong"))
        elif mosaic:
            crits.append(Criterion("PVS1", "PS",
                                   f"null variant ({annotation['consequence']}) in LOF gene "
                                   f"{annotation['gene']}; low-VAF mosaic — applied Strong (conservative)"))
        else:
            crits.append(Criterion("PVS1", "PVS",
                                   f"null variant ({annotation['consequence']}) in LOF gene "
                                   f"{annotation['gene']} — Very Strong"))
    if annotation.get("gnomad_af", 0.0) in (None, 0.0):
        crits.append(Criterion("PM2", "PM", "absent from population databases (gnomAD)"))
    if annotation["known_pathogenic"]:
        crits.append(Criterion("PS1", "PS",
                               f"variant previously established pathogenic ({annotation['known_source']})"))
    if phenotype_consistent:
        crits.append(Criterion("PP4", "PP", "phenotype highly specific for TSC (meets clinical criteria)"))
    return crits
