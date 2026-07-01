"""HLA allele screening module for drug hypersensitivity.

Screens patients for HLA alleles associated with severe adverse drug
reactions (Stevens-Johnson syndrome, toxic epidermal necrolysis, drug
reaction with eosinophilia and systemic symptoms, drug-induced liver
injury, and abacavir hypersensitivity syndrome).

Pre-emptive HLA screening before drug prescription can prevent
life-threatening hypersensitivity reactions.

References:
    - CPIC guideline for HLA genotype and carbamazepine (2013, updated 2018)
    - CPIC guideline for HLA-B*57:01 and abacavir (2012, updated 2014)
    - PharmGKB HLA clinical annotations
    - FDA Table of Pharmacogenomic Biomarkers in Drug Labeling

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════════════
# Enums and data classes
# ═════════════════════════════════════════════════════════════════════════════


class HLAStatus(str, Enum):
    """Risk status for an HLA-drug combination."""
    SAFE = "SAFE"
    CONTRAINDICATED = "CONTRAINDICATED"
    HIGH_RISK = "HIGH_RISK"
    UNKNOWN = "UNKNOWN"


class ReactionSeverity(str, Enum):
    """Severity of the potential hypersensitivity reaction."""
    FATAL = "fatal"
    SEVERE = "severe"
    MODERATE = "moderate"
    MILD = "mild"


@dataclass
class HLAScreenResult:
    """Result of screening a drug against HLA typing."""
    drug: str
    hla_allele: str
    status: HLAStatus
    reaction_type: str
    severity: ReactionSeverity
    recommendation: str
    alternatives: List[str] = field(default_factory=list)
    evidence_level: str = ""   # e.g. "CPIC Level A", "FDA Boxed Warning"
    population_risk: str = ""  # e.g. "High in Southeast Asian"


# ═════════════════════════════════════════════════════════════════════════════
# HLA-Drug Association Knowledge Base
# ═════════════════════════════════════════════════════════════════════════════

# Each entry: hla_allele, reaction_type, severity, recommendation,
#             alternatives, evidence_level, population_risk
HLA_DRUG_ASSOCIATIONS: Dict[str, List[Dict[str, Any]]] = {
    # ── Abacavir (HIV) ──
    "abacavir": [
        {
            "hla_allele": "HLA-B*57:01",
            "reaction_type": "Abacavir Hypersensitivity Syndrome (AHS)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "CONTRAINDICATED. Do NOT prescribe abacavir to HLA-B*57:01 "
                "positive patients. AHS occurs in ~48% of carriers and can "
                "be fatal on rechallenge."
            ),
            "alternatives": [
                "tenofovir disoproxil fumarate (TDF)",
                "tenofovir alafenamide (TAF)",
                "emtricitabine",
            ],
            "evidence_level": "CPIC Level A; FDA Boxed Warning",
            "population_risk": "~6-8% European, ~2-3% African American, ~1% East Asian",
        },
    ],
    # ── Carbamazepine (antiepileptic) ──
    "carbamazepine": [
        {
            "hla_allele": "HLA-B*15:02",
            "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
            "severity": ReactionSeverity.FATAL,
            "recommendation": (
                "CONTRAINDICATED in HLA-B*15:02 carriers. Risk of fatal "
                "SJS/TEN. FDA Boxed Warning requires testing in patients "
                "of Southeast Asian ancestry before prescribing."
            ),
            "alternatives": [
                "levetiracetam",
                "valproic acid",
                "lamotrigine (if HLA-B*15:02 negative and HLA-A*31:01 negative)",
                "lacosamide",
            ],
            "evidence_level": "CPIC Level A; FDA Boxed Warning",
            "population_risk": "~8% Southeast Asian, ~6% South Asian, <1% European",
        },
        {
            "hla_allele": "HLA-A*31:01",
            "reaction_type": "Drug Reaction with Eosinophilia and Systemic Symptoms (DRESS) / SJS/TEN / MPE",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. Avoid carbamazepine in HLA-A*31:01 carriers. "
                "5-26% risk of cutaneous adverse reactions including DRESS "
                "and SJS/TEN."
            ),
            "alternatives": [
                "levetiracetam",
                "valproic acid",
                "lacosamide",
            ],
            "evidence_level": "CPIC Level A",
            "population_risk": "~3-5% European, ~5-10% Japanese, ~2-3% African",
        },
    ],
    # ── Oxcarbazepine (antiepileptic) ──
    "oxcarbazepine": [
        {
            "hla_allele": "HLA-B*15:02",
            "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
            "severity": ReactionSeverity.FATAL,
            "recommendation": (
                "CONTRAINDICATED in HLA-B*15:02 carriers. Cross-reactivity "
                "with carbamazepine. FDA warns against use without testing."
            ),
            "alternatives": [
                "levetiracetam",
                "valproic acid",
                "lacosamide",
            ],
            "evidence_level": "CPIC Level A; FDA Warning",
            "population_risk": "~8% Southeast Asian",
        },
    ],
    # ── Phenytoin (antiepileptic) ──
    "phenytoin": [
        {
            "hla_allele": "HLA-B*15:02",
            "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. Consider avoiding phenytoin in HLA-B*15:02 "
                "carriers. Lower absolute risk than carbamazepine but "
                "still significantly elevated."
            ),
            "alternatives": [
                "levetiracetam",
                "valproic acid",
                "lacosamide",
            ],
            "evidence_level": "CPIC Level A",
            "population_risk": "~8% Southeast Asian",
        },
    ],
    # ── Allopurinol (gout) ──
    "allopurinol": [
        {
            "hla_allele": "HLA-B*58:01",
            "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN) / DRESS",
            "severity": ReactionSeverity.FATAL,
            "recommendation": (
                "CONTRAINDICATED in HLA-B*58:01 carriers. Use alternative "
                "urate-lowering therapy. ACR conditionally recommends "
                "HLA-B*58:01 testing prior to allopurinol, especially in "
                "Southeast Asian and African American patients."
            ),
            "alternatives": [
                "febuxostat",
                "probenecid",
                "pegloticase (refractory)",
            ],
            "evidence_level": "CPIC Level A; ACR Conditional Recommendation",
            "population_risk": "~6-8% Southeast Asian, ~6% Korean, ~3-4% African American, ~1-2% European",
        },
    ],
    # ── Flucloxacillin (antibiotic) ──
    "flucloxacillin": [
        {
            "hla_allele": "HLA-B*57:01",
            "reaction_type": "Drug-Induced Liver Injury (DILI)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. HLA-B*57:01 carriers have ~80-fold increased "
                "risk of flucloxacillin DILI. Consider alternative "
                "anti-staphylococcal agent."
            ),
            "alternatives": [
                "dicloxacillin",
                "cefazolin",
                "clindamycin",
                "trimethoprim-sulfamethoxazole",
            ],
            "evidence_level": "Strong association; prescreening not yet standard",
            "population_risk": "~6-8% European",
        },
    ],
    # ── Lamotrigine (antiepileptic) ──
    "lamotrigine": [
        {
            "hla_allele": "HLA-B*15:02",
            "reaction_type": "Stevens-Johnson Syndrome (SJS)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. Elevated SJS risk in HLA-B*15:02 carriers, "
                "though lower than carbamazepine. Use slow titration if "
                "prescribed; consider alternative."
            ),
            "alternatives": [
                "levetiracetam",
                "valproic acid",
                "lacosamide",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "~8% Southeast Asian",
        },
    ],
    # ── Dapsone (dermatology, leprosy) ──
    "dapsone": [
        {
            "hla_allele": "HLA-B*13:01",
            "reaction_type": "Dapsone Hypersensitivity Syndrome (DHS)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "CONTRAINDICATED in HLA-B*13:01 carriers. DHS includes "
                "fever, rash, hepatitis, lymphadenopathy. Test before "
                "prescribing in at-risk populations."
            ),
            "alternatives": [
                "clofazimine",
                "minocycline",
            ],
            "evidence_level": "Strong association; recommended testing in endemic regions",
            "population_risk": "~2-10% Southeast Asian, ~1-2% South Asian",
        },
    ],
    # ── Ticlopidine (antiplatelet — rarely used) ──
    "ticlopidine": [
        {
            "hla_allele": "HLA-A*33:03",
            "reaction_type": "Drug-Induced Liver Injury (DILI)",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK of hepatotoxicity in HLA-A*33:03 carriers. "
                "Ticlopidine rarely used; prefer clopidogrel/prasugrel/"
                "ticagrelor."
            ),
            "alternatives": [
                "clopidogrel",
                "prasugrel",
                "ticagrelor",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "~10-15% Japanese, ~3-5% East Asian",
        },
    ],
    # ── Nevirapine (HIV) ──
    "nevirapine": [
        {
            "hla_allele": "HLA-B*35:05",
            "reaction_type": "Skin rash / SJS",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. Avoid nevirapine in HLA-B*35:05 carriers. "
                "Use alternative NNRTI (efavirenz, rilpivirine, doravirine) "
                "or integrase inhibitor-based regimen."
            ),
            "alternatives": [
                "efavirenz",
                "rilpivirine",
                "doravirine",
                "dolutegravir (INSTI-based regimen)",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "~5-15% Southeast Asian",
        },
        {
            "hla_allele": "HLA-DRB1*01:01",
            "reaction_type": "Nevirapine-associated hepatotoxicity",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK of hepatotoxicity. Particularly in patients with "
                "higher CD4 counts. Consider alternative ARV regimen."
            ),
            "alternatives": [
                "efavirenz",
                "dolutegravir",
                "rilpivirine",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "Variable by population",
        },
    ],
    # ── Sulfasalazine (DMARD) ──
    "sulfasalazine": [
        {
            "hla_allele": "HLA-B*13:01",
            "reaction_type": "DRESS/SJS",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. Avoid sulfasalazine in HLA-B*13:01 carriers. "
                "Use alternative DMARD. DRESS with hepatic involvement "
                "can be fatal."
            ),
            "alternatives": [
                "methotrexate",
                "leflunomide",
                "hydroxychloroquine",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "~2-10% Southeast Asian, ~1-2% South Asian",
        },
    ],
    # ── Methazolamide (carbonic anhydrase inhibitor) ──
    "methazolamide": [
        {
            "hla_allele": "HLA-B*59:01",
            "reaction_type": "SJS/TEN",
            "severity": ReactionSeverity.FATAL,
            "recommendation": (
                "CONTRAINDICATED in HLA-B*59:01 carriers. Life-threatening "
                "SJS/TEN reported. Test in patients of East Asian ancestry "
                "before prescribing."
            ),
            "alternatives": [
                "acetazolamide (with caution)",
                "dorzolamide (topical)",
            ],
            "evidence_level": "Strong association; region-specific testing recommended",
            "population_risk": "~3% Japanese, ~1% Korean, rare in European",
        },
    ],
    # ── Clozapine (atypical antipsychotic) ──
    "clozapine": [
        {
            "hla_allele": "HLA-DQB1*05:02",
            "reaction_type": "Agranulocytosis",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. HLA-DQB1*05:02 carriers have significantly "
                "elevated risk of clozapine-induced agranulocytosis. "
                "Enhanced monitoring required; consider alternative "
                "atypical antipsychotic if available."
            ),
            "alternatives": [
                "olanzapine",
                "quetiapine",
                "cariprazine",
            ],
            "evidence_level": "Moderate evidence; ANC monitoring still required regardless",
            "population_risk": "~10-15% Ashkenazi Jewish, ~5-8% European, ~2-3% East Asian",
        },
    ],
    # ── Trimethoprim-sulfamethoxazole (antibiotic) ──
    "trimethoprim-sulfamethoxazole": [
        {
            "hla_allele": "HLA-B*38:02",
            "reaction_type": "SJS/TEN",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK of severe cutaneous adverse reactions. "
                "Consider alternative antibiotic."
            ),
            "alternatives": [
                "nitrofurantoin (for UTI)",
                "amoxicillin-clavulanate",
                "fluoroquinolone (with caution)",
            ],
            "evidence_level": "Emerging evidence",
            "population_risk": "~1-3% European, variable by population",
        },
    ],
    # ── Minocycline (tetracycline antibiotic) ──
    "minocycline": [
        {
            "hla_allele": "HLA-B*35:02",
            "reaction_type": "DRESS",
            "severity": ReactionSeverity.SEVERE,
            "recommendation": (
                "HIGH RISK. DRESS syndrome with hepatic involvement "
                "reported in HLA-B*35:02 carriers. Consider alternative "
                "tetracycline or different antibiotic class."
            ),
            "alternatives": [
                "doxycycline",
                "amoxicillin",
                "azithromycin",
            ],
            "evidence_level": "Moderate evidence",
            "population_risk": "~3-5% East Asian, ~1-2% European",
        },
    ],
}


# ═════════════════════════════════════════════════════════════════════════════
# HLA Screener
# ═════════════════════════════════════════════════════════════════════════════


class HLAScreener:
    """Screens for HLA alleles associated with drug hypersensitivity.

    Supports two workflows:
    1. Pre-prescription screening: given a drug, check if the patient
       carries any risk alleles.
    2. Panel screening: given HLA typing results, identify all
       contraindicated/high-risk drugs.
    """

    def __init__(
        self,
        associations: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    ) -> None:
        """Initialize with HLA-drug associations.

        Parameters
        ----------
        associations : dict, optional
            Custom HLA-drug associations. Defaults to the built-in
            HLA_DRUG_ASSOCIATIONS knowledge base.
        """
        self._associations = associations or HLA_DRUG_ASSOCIATIONS

    def screen_drug(
        self,
        drug: str,
        hla_typing: Dict[str, List[str]],
    ) -> List[HLAScreenResult]:
        """Screen a single drug against the patient's HLA typing.

        Parameters
        ----------
        drug : str
            Drug name (case-insensitive).
        hla_typing : dict
            HLA locus → list of alleles.
            Example: {"HLA-B": ["*57:01", "*44:02"], "HLA-A": ["*31:01", "*02:01"]}

        Returns
        -------
        list of HLAScreenResult
            One result per matching HLA-drug association.  Empty list if
            no risk alleles detected (drug is safe to prescribe).
        """
        drug_lower = drug.lower().strip()
        assocs = self._associations.get(drug_lower, [])

        if not assocs:
            return [HLAScreenResult(
                drug=drug,
                hla_allele="N/A",
                status=HLAStatus.UNKNOWN,
                reaction_type="No pharmacogenomic HLA data available",
                severity=ReactionSeverity.MILD,
                recommendation=(
                    f"No HLA-drug association data for {drug}. "
                    f"Prescribe per standard guidelines."
                ),
            )]

        # Flatten patient alleles for matching
        patient_alleles = self._flatten_hla_typing(hla_typing)
        results: List[HLAScreenResult] = []

        for assoc in assocs:
            risk_allele = assoc["hla_allele"]
            is_carrier = self._check_allele_match(risk_allele, patient_alleles)

            if is_carrier:
                status = self._determine_status(assoc["severity"])
                results.append(HLAScreenResult(
                    drug=drug,
                    hla_allele=risk_allele,
                    status=status,
                    reaction_type=assoc["reaction_type"],
                    severity=assoc["severity"],
                    recommendation=assoc["recommendation"],
                    alternatives=assoc.get("alternatives", []),
                    evidence_level=assoc.get("evidence_level", ""),
                    population_risk=assoc.get("population_risk", ""),
                ))

        # If no risk alleles found, drug is safe
        if not results:
            results.append(HLAScreenResult(
                drug=drug,
                hla_allele="None detected",
                status=HLAStatus.SAFE,
                reaction_type="No risk alleles detected",
                severity=ReactionSeverity.MILD,
                recommendation=(
                    f"{drug} is safe to prescribe based on HLA screening. "
                    f"No carrier status for known risk alleles."
                ),
            ))

        return results

    def screen_all_drugs(
        self,
        hla_typing: Dict[str, List[str]],
    ) -> List[HLAScreenResult]:
        """Screen all drugs in the knowledge base against patient HLA typing.

        Useful for generating a comprehensive pre-emptive pharmacogenomic
        HLA report.

        Parameters
        ----------
        hla_typing : dict
            HLA locus → list of alleles.

        Returns
        -------
        list of HLAScreenResult
            All CONTRAINDICATED and HIGH_RISK results, sorted by severity.
        """
        severity_order = {
            ReactionSeverity.FATAL: 0,
            ReactionSeverity.SEVERE: 1,
            ReactionSeverity.MODERATE: 2,
            ReactionSeverity.MILD: 3,
        }

        all_results: List[HLAScreenResult] = []

        for drug in self._associations:
            results = self.screen_drug(drug, hla_typing)
            # Only include actionable results (not SAFE / UNKNOWN)
            for r in results:
                if r.status in (HLAStatus.CONTRAINDICATED, HLAStatus.HIGH_RISK):
                    all_results.append(r)

        all_results.sort(key=lambda r: severity_order.get(r.severity, 99))

        logger.info(
            "HLA panel screen: %d drugs checked, %d contraindicated, %d high-risk",
            len(self._associations),
            sum(1 for r in all_results if r.status == HLAStatus.CONTRAINDICATED),
            sum(1 for r in all_results if r.status == HLAStatus.HIGH_RISK),
        )
        return all_results

    def get_hla_from_vcf(
        self,
        vcf_path: str | Path,
    ) -> Dict[str, List[str]]:
        """Extract HLA types from a VCF file.

        Supports VCF files produced by HLA typing tools such as:
        - OptiType (HLA class I)
        - HLA-HD
        - xHLA
        - HLA*LA

        The method looks for HLA-related records in the VCF by searching
        for variants on chr6 in the MHC region (6p21.3, roughly
        28,477,797-33,448,354 on GRCh38) and for HLA allele annotations
        in the INFO field.

        Parameters
        ----------
        vcf_path : str or Path
            Path to a VCF file containing HLA typing results.

        Returns
        -------
        dict
            Locus → list of alleles.
            Example: {"HLA-A": ["*02:01", "*31:01"], "HLA-B": ["*57:01", "*44:02"]}
        """
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")

        hla_typing: Dict[str, List[str]] = {}
        # MHC region on GRCh38 (extended)
        mhc_start = 28_400_000
        mhc_end = 33_500_000

        # Regex patterns for HLA allele notation in INFO fields
        hla_pattern = re.compile(
            r"(HLA-[A-Z][A-Z0-9]*)\*(\d{2,3}:\d{2,3}(?::\d{2,3})?(?::\d{2,3})?)"
        )

        open_fn = self._open_vcf(vcf_path)

        with open_fn(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    # Check header comments for HLA tool output
                    match = hla_pattern.search(line)
                    if match:
                        locus = match.group(1)
                        allele = f"*{match.group(2)}"
                        hla_typing.setdefault(locus, [])
                        if allele not in hla_typing[locus]:
                            hla_typing[locus].append(allele)
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0].replace("chr", "")
                if chrom != "6":
                    continue

                try:
                    pos = int(fields[1])
                except ValueError:
                    continue

                # Only process MHC region
                if not (mhc_start <= pos <= mhc_end):
                    continue

                # Search INFO field for HLA allele annotations
                info = fields[7] if len(fields) > 7 else ""
                for match in hla_pattern.finditer(info):
                    locus = match.group(1)
                    allele = f"*{match.group(2)}"
                    hla_typing.setdefault(locus, [])
                    if allele not in hla_typing[locus]:
                        hla_typing[locus].append(allele)

                # Check ID field (some tools put HLA alleles here)
                id_field = fields[2]
                for match in hla_pattern.finditer(id_field):
                    locus = match.group(1)
                    allele = f"*{match.group(2)}"
                    hla_typing.setdefault(locus, [])
                    if allele not in hla_typing[locus]:
                        hla_typing[locus].append(allele)

        if not hla_typing:
            logger.warning(
                "No HLA alleles extracted from %s. Ensure VCF contains "
                "HLA typing annotations (e.g., from OptiType or HLA-HD).",
                vcf_path.name,
            )

        logger.info(
            "Extracted HLA typing from %s: %s",
            vcf_path.name,
            {k: v for k, v in hla_typing.items()},
        )
        return hla_typing

    # ── Private helpers ─────────────────────────────────────────────────

    @staticmethod
    def _flatten_hla_typing(hla_typing: Dict[str, List[str]]) -> Set[str]:
        """Flatten HLA typing dict into a set of fully qualified allele names.

        Input:  {"HLA-B": ["*57:01", "*44:02"]}
        Output: {"HLA-B*57:01", "HLA-B*44:02"}
        """
        alleles: Set[str] = set()
        for locus, allele_list in hla_typing.items():
            for allele in allele_list:
                # Normalize: ensure allele starts with *
                if not allele.startswith("*"):
                    allele = f"*{allele}"
                # Build full name: HLA-B*57:01
                full = f"{locus}{allele}"
                alleles.add(full)
        return alleles

    @staticmethod
    def _check_allele_match(risk_allele: str, patient_alleles: Set[str]) -> bool:
        """Check if the patient carries a specific risk allele.

        Supports forward resolution matching only: a patient allele at
        higher resolution (e.g., HLA-B*15:02:01) matches a two-field risk
        allele (HLA-B*15:02).  Reverse matching (one-field patient allele
        matching a two-field risk) is NOT performed to avoid false positives
        (e.g., HLA-B*15:01 should not match HLA-B*15:02 risk).
        """
        if risk_allele in patient_alleles:
            return True

        # Forward match: patient allele at higher resolution matches risk allele
        # e.g., patient HLA-B*15:02:01 matches risk HLA-B*15:02
        for patient_allele in patient_alleles:
            if patient_allele.startswith(risk_allele):
                return True

        return False

    @staticmethod
    def _determine_status(severity: ReactionSeverity) -> HLAStatus:
        """Map reaction severity to HLA status."""
        if severity in (ReactionSeverity.FATAL,):
            return HLAStatus.CONTRAINDICATED
        if severity in (ReactionSeverity.SEVERE,):
            return HLAStatus.CONTRAINDICATED
        if severity == ReactionSeverity.MODERATE:
            return HLAStatus.HIGH_RISK
        return HLAStatus.HIGH_RISK  # default to cautious

    @staticmethod
    def _open_vcf(path: Path):
        """Return the appropriate open function for plain or gzipped VCF."""
        if path.suffix == ".gz" or str(path).endswith(".vcf.gz"):
            import gzip
            return gzip.open
        return open
