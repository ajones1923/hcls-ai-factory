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


# ── PP4: "phenotype highly specific for the disease" ────────────────────────
# ACMG/AMP 2015 (Richards et al., Genet Med 17:405, Table 3) defines PP4 as a patient
# phenotype or family history highly specific for a disease with a single genetic
# aetiology. The caller used to grant it for `bool(cohort["phenotypes"])` -- i.e. for any
# recorded phenotype at all -- so a patient whose only entry was "Seizure" earned a
# supporting pathogenic criterion on a finding that is not specific to anything.
#
# Specificity here is the published TSC clinical diagnostic criterion the reason string
# already invokes: a definite clinical diagnosis requires two major features, or one major
# plus two or more minor features (Northrup et al., Updated International Tuberous Sclerosis
# Complex Diagnostic Criteria and Surveillance and Management Recommendations, Pediatr Neurol
# 2021;123:50-66, doi:10.1016/j.pediatrneurol.2021.07.011).
#
# Seizures, infantile spasms, intellectual disability and autism are deliberately absent
# below. They are common in TSC and they are not diagnostic criteria -- TAND and epilepsy
# are consequences of the disease, not features specific to it.
#
# LIMITATION, recorded rather than papered over: three of the major features are defined
# with counts (>=3 hypomelanotic macules, >=3 angiofibromas, >=2 angiomyolipomas) and the
# cohort's phenotype records carry an HPO term with no count, so presence is all that can be
# read. See CLINICAL_REVIEW.md.
TSC_MAJOR_HPO = {
    "HP:0009717": "cortical tubers / radial migration lines",
    "HP:0009718": "subependymal giant-cell astrocytoma",
    "HP:0009719": "hypomelanotic macules (>=3; count not recorded)",
    "HP:0009720": "facial angiofibromas (>=3; count not recorded)",
    "HP:0006772": "renal angiomyolipoma (>=2; count not recorded)",
    "HP:0009721": "shagreen patch",
    "HP:0001629": "cardiac rhabdomyoma",
    "HP:0030004": "lymphangioleiomyomatosis",
    "HP:0009722": "ungual fibromas (>=2)",
    "HP:0009737": "multiple retinal hamartomas",
}

TSC_MINOR_HPO = {
    "HP:0007449": "confetti skin lesions",
    "HP:0006297": "dental enamel pits (>3)",
    "HP:0000167": "intraoral fibromas (>=2)",
    "HP:0007968": "retinal achromic patch",
    "HP:0000107": "multiple renal cysts",
}


def meets_tsc_clinical_criteria(phenotypes) -> tuple[bool, str]:
    """Whether the recorded phenotypes meet the TSC clinical diagnostic criteria.

    Returns `(met, basis)`; `basis` names the features counted, so the criterion's reason
    string states what the call rests on instead of asserting specificity.
    """
    major, minor = [], []
    for p in phenotypes or []:
        hpo_id = (p or {}).get("hpo_id") if isinstance(p, dict) else None
        if hpo_id in TSC_MAJOR_HPO:
            if TSC_MAJOR_HPO[hpo_id] not in major:
                major.append(TSC_MAJOR_HPO[hpo_id])
        elif hpo_id in TSC_MINOR_HPO:
            if TSC_MINOR_HPO[hpo_id] not in minor:
                minor.append(TSC_MINOR_HPO[hpo_id])
    met = len(major) >= 2 or (len(major) >= 1 and len(minor) >= 2)
    if met:
        basis = f"major: {', '.join(major)}" + (f"; minor: {', '.join(minor)}" if minor else "")
    else:
        basis = (f"{len(major)} major, {len(minor)} minor recorded — below the 2-major / "
                 f"1-major-plus-2-minor threshold")
    return met, basis


def assign_criteria(annotation: dict, *, mosaic: bool, phenotype_consistent: bool,
                    phenotype_basis: str = "") -> list[Criterion]:
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
        reason = "phenotype meets the TSC clinical diagnostic criteria (Northrup 2021)"
        crits.append(Criterion("PP4", "PP", f"{reason} — {phenotype_basis}" if phenotype_basis
                               else reason))
    return crits
