"""
ACMG Secondary Findings (SF) categorization — E6.

A rules layer over annotated variants: flag variants in genes on the ACMG SF panel whose
clinical significance is reportable (Pathogenic / Likely pathogenic). This is a deterministic
join — gene membership x significance — not a model. It consumes the platform's existing
annotation output (gene + clinical_significance), so it adds a clinically meaningful capability
the reference designs lacked, with no new infrastructure.

PROVENANCE / VALIDATION: ACMG_SF_PANEL below is a curated representation of the ACMG SF v3.2
gene list (Miller et al., Genet Med 2023). The canonical list is updated periodically (v3.3 now
supersedes v3.2). Before any clinical use, validate/refresh the panel and its
phenotype/inheritance metadata against the current official ACMG publication. The categorization
logic is version-independent; only the data table changes.
"""
from __future__ import annotations

from typing import Any

PANEL_VERSION = "ACMG SF v3.2 (curated; validate against the official publication before clinical use)"

# Curated ACMG SF v3.2 panel, grouped by condition area. Reportable per ACMG recommendations.
ACMG_SF_PANEL: dict[str, str] = {
    # Hereditary cancer
    "BRCA1": "Hereditary breast and ovarian cancer", "BRCA2": "Hereditary breast and ovarian cancer",
    "PALB2": "Hereditary breast cancer", "TP53": "Li-Fraumeni syndrome",
    "STK11": "Peutz-Jeghers syndrome", "MLH1": "Lynch syndrome", "MSH2": "Lynch syndrome",
    "MSH6": "Lynch syndrome", "PMS2": "Lynch syndrome", "APC": "Familial adenomatous polyposis",
    "MUTYH": "MUTYH-associated polyposis", "BMPR1A": "Juvenile polyposis", "SMAD4": "Juvenile polyposis",
    "PTEN": "PTEN hamartoma tumor syndrome", "RB1": "Retinoblastoma", "MEN1": "Multiple endocrine neoplasia type 1",
    "RET": "Multiple endocrine neoplasia type 2", "VHL": "Von Hippel-Lindau syndrome",
    "NF2": "Neurofibromatosis type 2", "TSC1": "Tuberous sclerosis complex", "TSC2": "Tuberous sclerosis complex",
    "WT1": "WT1-related Wilms tumor", "SDHD": "Hereditary paraganglioma-pheochromocytoma",
    "SDHAF2": "Hereditary paraganglioma-pheochromocytoma", "SDHC": "Hereditary paraganglioma-pheochromocytoma",
    "SDHB": "Hereditary paraganglioma-pheochromocytoma", "MAX": "Hereditary paraganglioma-pheochromocytoma",
    "TMEM127": "Hereditary paraganglioma-pheochromocytoma", "CDH1": "Hereditary diffuse gastric cancer",
    # Cardiovascular
    "LDLR": "Familial hypercholesterolemia", "APOB": "Familial hypercholesterolemia",
    "PCSK9": "Familial hypercholesterolemia", "MYH7": "Cardiomyopathy", "MYBPC3": "Cardiomyopathy",
    "TNNT2": "Cardiomyopathy", "TNNI3": "Cardiomyopathy", "TPM1": "Cardiomyopathy", "MYL3": "Cardiomyopathy",
    "ACTC1": "Cardiomyopathy", "PRKAG2": "Cardiomyopathy", "MYL2": "Cardiomyopathy", "LMNA": "Cardiomyopathy",
    "BAG3": "Dilated cardiomyopathy", "DES": "Dilated cardiomyopathy", "RBM20": "Dilated cardiomyopathy",
    "TNNC1": "Cardiomyopathy", "FLNC": "Cardiomyopathy", "PKP2": "Arrhythmogenic right ventricular cardiomyopathy",
    "DSP": "Arrhythmogenic right ventricular cardiomyopathy", "DSC2": "Arrhythmogenic right ventricular cardiomyopathy",
    "DSG2": "Arrhythmogenic right ventricular cardiomyopathy", "TMEM43": "Arrhythmogenic right ventricular cardiomyopathy",
    "RYR2": "Catecholaminergic polymorphic VT", "CASQ2": "Catecholaminergic polymorphic VT",
    "TRDN": "Catecholaminergic polymorphic VT", "CALM1": "Catecholaminergic polymorphic VT",
    "CALM2": "Catecholaminergic polymorphic VT", "CALM3": "Catecholaminergic polymorphic VT",
    "KCNQ1": "Long QT syndrome", "KCNH2": "Long QT syndrome", "SCN5A": "Long QT / Brugada syndrome",
    "FBN1": "Marfan / thoracic aortic aneurysm", "TGFBR1": "Loeys-Dietz syndrome", "TGFBR2": "Loeys-Dietz syndrome",
    "SMAD3": "Loeys-Dietz syndrome", "ACTA2": "Familial thoracic aortic aneurysm", "MYH11": "Familial thoracic aortic aneurysm",
    "LOX": "Familial thoracic aortic aneurysm", "COL3A1": "Ehlers-Danlos syndrome, vascular type",
    # Metabolic / other (treatable)
    "GLA": "Fabry disease", "GAA": "Pompe disease", "ATP7B": "Wilson disease",
    "OTC": "Ornithine transcarbamylase deficiency", "BTD": "Biotinidase deficiency",
    "RYR1": "Malignant hyperthermia susceptibility", "CACNA1S": "Malignant hyperthermia susceptibility",
    "HFE": "Hereditary hemochromatosis", "TTR": "Hereditary transthyretin amyloidosis",
    "HNF1A": "Maturity-onset diabetes of the young", "RPE65": "RPE65-related retinopathy",
    "ACVRL1": "Hereditary hemorrhagic telangiectasia", "ENG": "Hereditary hemorrhagic telangiectasia",
}

_REPORTABLE_SIGNIFICANCE = {"pathogenic", "likely_pathogenic", "likely pathogenic",
                            "pathogenic/likely_pathogenic", "pathogenic/likely pathogenic"}


def is_on_panel(gene: str) -> bool:
    return bool(gene) and gene.upper() in ACMG_SF_PANEL


def is_reportable(gene: str, significance: str | None) -> bool:
    """Reportable secondary finding: gene on the ACMG SF panel AND a (likely) pathogenic call."""
    if not is_on_panel(gene) or not significance:
        return False
    return significance.strip().lower() in _REPORTABLE_SIGNIFICANCE


def secondary_findings(annotated_variants: list[dict[str, Any]],
                       gene_key: str = "gene", sig_key: str = "clinical_significance") -> list[dict[str, Any]]:
    """Filter annotated variants to reportable ACMG SF secondary findings.

    Each input variant is a dict with at least a gene symbol and a clinical-significance field
    (as produced by the annotation layer). Returns the reportable subset, each enriched with the
    panel condition.
    """
    out = []
    for v in annotated_variants:
        gene, sig = v.get(gene_key, ""), v.get(sig_key)
        if is_reportable(gene, sig):
            out.append({**v, "acmg_sf_condition": ACMG_SF_PANEL[gene.upper()],
                        "acmg_sf_panel": PANEL_VERSION})
    return out


def panel_summary() -> dict[str, Any]:
    conditions: dict[str, int] = {}
    for cond in ACMG_SF_PANEL.values():
        conditions[cond] = conditions.get(cond, 0) + 1
    return {"version": PANEL_VERSION, "n_genes": len(ACMG_SF_PANEL), "n_conditions": len(conditions)}
