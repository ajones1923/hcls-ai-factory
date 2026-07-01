"""
Synthetic cohort roster (PRD §3 FR-CG-1; master paper §15).

Deterministic 50-patient composition matching real TSC epidemiology, with the three
featured patients (A/B/C) given exact specs. Same seed -> same roster. This module
defines WHO is in the cohort; genomic/clinical artifacts are generated downstream.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import yaml

from config.settings import settings

# Representative TSC variants (cDNA, consequence) with approximate GRCh38 coordinates.
# TSC2: chr16p13.3 (~2.05-2.09 Mb); TSC1: chr9q34.13 (~132.89-132.95 Mb).
TSC2_VARIANTS = [
    {"cdna": "c.3037C>T", "protein": "p.Arg1013Ter", "kind": "nonsense", "chrom": "chr16", "pos": 2074000, "ref": "C", "alt": "T"},
    {"cdna": "c.5024C>T", "protein": "p.Pro1675Leu", "kind": "missense", "chrom": "chr16", "pos": 2085000, "ref": "C", "alt": "T"},
    {"cdna": "c.2251C>T", "protein": "p.Arg751Ter", "kind": "nonsense", "chrom": "chr16", "pos": 2069000, "ref": "C", "alt": "T"},
    {"cdna": "c.4255del", "protein": "p.frameshift", "kind": "frameshift", "chrom": "chr16", "pos": 2080100, "ref": "AG", "alt": "A"},
]
TSC1_VARIANTS = [
    {"cdna": "c.1888C>T", "protein": "p.Arg630Ter", "kind": "nonsense", "chrom": "chr9", "pos": 132905000, "ref": "C", "alt": "T"},
    {"cdna": "c.2356C>T", "protein": "p.Arg786Ter", "kind": "nonsense", "chrom": "chr9", "pos": 132910000, "ref": "C", "alt": "T"},
    {"cdna": "c.733del", "protein": "p.frameshift", "kind": "frameshift", "chrom": "chr9", "pos": 132898050, "ref": "TC", "alt": "T"},
]


@dataclass
class PatientSpec:
    patient_id: str
    gene: str | None              # "TSC1" | "TSC2" | None (NMI)
    zygosity: str                 # "germline" | "mosaic" | "none"
    variant: dict | None          # one of the *_VARIANTS dicts
    vaf: float                    # alt allele fraction in the sequenced tissue
    tissue_available: bool        # mosaic recovery needs affected tissue
    sex: str
    age: int
    severity: str
    featured: str | None = None   # "A" | "B" | "C" | None
    expected_acmg: str = "Pathogenic"
    extras: dict = field(default_factory=dict)   # SEGA series, TAND signals, etc.


def _demo() -> dict:
    with open(settings.DEMO_CONFIG_PATH) as f:
        return yaml.safe_load(f)


def build_roster(seed: int | None = None) -> list[PatientSpec]:
    cfg = _demo()
    comp = cfg["cohort"]["composition"]
    vlo, vhi = cfg["cohort"]["mosaic_vaf_range"]
    rng = random.Random(seed if seed is not None else cfg["cohort"]["seed"])

    plan: list[tuple[str, str]] = (
        [("TSC2", "germline")] * comp["tsc2_germline"]
        + [("TSC1", "germline")] * comp["tsc1_germline"]
        + [("TSC2", "mosaic")] * comp["tsc2_mosaic"]
        + [("TSC1", "mosaic")] * comp["tsc1_mosaic"]
        + [("__NMI__", "none")] * comp["no_mutation_identified"]
    )
    specs: list[PatientSpec] = []
    for i, (gene, zyg) in enumerate(plan, start=1):
        pid = f"TSC-{i:04d}"
        if gene == "__NMI__":
            specs.append(PatientSpec(pid, None, "none", None, 0.0, False,
                                     rng.choice("MF"), rng.randint(1, 25), "moderate",
                                     expected_acmg="No variant identified"))
            continue
        variants = TSC2_VARIANTS if gene == "TSC2" else TSC1_VARIANTS
        variant = rng.choice(variants)
        # expected_acmg is the GROUND-TRUTH ACMG call a clinical geneticist would make for this
        # variant + zygosity (consumed only by the evaluation harness). It follows the actual
        # ACMG/AMP rules, not optimism: a bare missense with only PM2 is a VUS; an established
        # recurrent null is Pathogenic even at low VAF; a NOVEL null seen only as low-VAF mosaic
        # gets the conservative PVS1->Strong -> Likely Pathogenic. (In this synthetic set the
        # established recurrent nulls are the nonsense variants; the frameshifts are novel.)
        if zyg == "mosaic":
            vaf = round(rng.uniform(vlo, vhi), 3)
            if variant["kind"] == "missense":
                acmg = "Variant of Uncertain Significance"
            elif variant["kind"] == "nonsense":          # established recurrent null
                acmg = "Pathogenic"
            else:                                         # novel frameshift, low-VAF mosaic
                acmg = "Likely Pathogenic (PVS1 Strong, PM2, PP4)"
        else:
            vaf = 0.5
            acmg = "Pathogenic" if variant["kind"] in ("nonsense", "frameshift") else "Variant of Uncertain Significance"
        specs.append(PatientSpec(pid, gene, zyg, variant, vaf, zyg == "mosaic",
                                 rng.choice("MF"), rng.randint(1, 25),
                                 rng.choice(["moderate", "moderate-severe", "severe"]),
                                 expected_acmg=acmg))

    _apply_featured(specs, cfg, rng)
    return specs


def _apply_featured(specs: list[PatientSpec], cfg: dict, rng: random.Random) -> None:
    """Override three slots with the exact featured-patient specs (master paper §15)."""
    by_profile = {
        "A": next(s for s in specs if s.gene == "TSC2" and s.zygosity == "mosaic"),
        "B": next(s for s in specs if s.gene == "TSC2" and s.zygosity == "germline"),
        "C": next(s for s in specs if s.gene == "TSC1" and s.zygosity == "germline"),
    }
    fp = cfg["featured_patients"]
    # A — 4yo F, NMI on blood, 8.3% VAF TSC2 frameshift in resected tuber tissue
    a = by_profile["A"]
    a.featured, a.age, a.sex, a.vaf, a.tissue_available = "A", fp["A"]["age"], "F", fp["A"]["tissue_vaf"], True
    a.variant = TSC2_VARIANTS[3]  # frameshift
    a.expected_acmg = fp["A"]["expected_acmg"]
    a.extras = {"tand_clusters": {}}   # the quiet-TAND specificity control (no spurious flags)
    # B — 12yo M, TSC2 c.3037C>T p.Arg1013Ter, SEGA series, scattered TAND
    b = by_profile["B"]
    b.featured, b.age, b.sex = "B", fp["B"]["age"], "M"
    b.variant = TSC2_VARIANTS[0]        # c.3037C>T nonsense, established recurrent
    b.expected_acmg = "Pathogenic"      # recompute ground truth after the variant override
    b.extras = {
        "sega_series_cm": fp["B"]["sega_series_cm"], "sega_months": [-24, -12, -6],
        "sega_location": "foramen of Monro", "aml_max_cm": fp["B"]["aml_max_cm"],
        "tand_clusters": {"academic": 3, "psychiatric": 1},
        # well-controlled focal epilepsy — a quantity that stays comfortably below threshold
        "seizure_freq_series": [2, 1, 1], "seizure_months": [-24, -12, -6],
    }
    # C — 18yo F, TSC1, partial everolimus response, AML ~4cm, refractory seizures
    c = by_profile["C"]
    c.featured, c.age, c.sex = "C", fp["C"]["age"], "F"
    c.variant = TSC1_VARIANTS[0]        # c.1888C>T nonsense — definite TSC1 diagnosis (matches narrative)
    c.expected_acmg = "Pathogenic"
    c.extras = {"aml_max_cm": fp["C"]["aml_max_cm"], "therapy": fp["C"]["therapy"],
                "refractory_epilepsy": True, "seizures": "refractory focal",
                "aml_series_cm": [3.2, 3.6, 4.0], "aml_months": [-24, -12, -6],
                # renal decline alongside the growing AML, and refractory rising seizures
                "egfr_series": [95, 84, 76], "egfr_months": [-24, -12, -6],
                "seizure_freq_series": [8, 11, 13], "seizure_months": [-24, -12, -6],
                "trial_eligible": True}
