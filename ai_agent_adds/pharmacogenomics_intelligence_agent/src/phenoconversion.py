"""Phenoconversion detection module for Pharmacogenomics Intelligence Agent.

Phenoconversion occurs when a concomitant medication inhibits or induces a
CYP enzyme, causing the patient's *effective* metabolizer phenotype to differ
from their *genetic* phenotype.  For example, a CYP2D6 Normal Metabolizer
taking fluoxetine (a strong CYP2D6 inhibitor) effectively becomes a Poor
Metabolizer for all CYP2D6 substrates.

This module detects these drug-drug-gene interactions and adjusts the
effective phenotype accordingly.

References:
    - Shah & Smith, Genes (Basel) 2015; 6(4): 1000-1007
    - Storelli et al., Clin Pharmacol Ther 2018; 104(3): 559-568
    - CPIC phenoconversion guidance

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════════════
# CYP Inhibitor / Inducer Knowledge Base
# ═════════════════════════════════════════════════════════════════════════════

class InhibitorStrength(str, Enum):
    """FDA classification of CYP inhibition potency."""
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"


class InducerStrength(str, Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"


# Drug → CYP enzyme → inhibitor strength
# Sources: FDA Drug Development and Drug Interactions Table,
#          Flockhart cytochrome P450 table (Indiana University)
CYP_INHIBITORS: Dict[str, Dict[str, InhibitorStrength]] = {
    # ── CYP2D6 Inhibitors ──
    "fluoxetine":    {"CYP2D6": InhibitorStrength.STRONG},
    "paroxetine":    {"CYP2D6": InhibitorStrength.STRONG},
    "bupropion":     {"CYP2D6": InhibitorStrength.STRONG},
    "quinidine":     {"CYP2D6": InhibitorStrength.STRONG},
    "cinacalcet":    {"CYP2D6": InhibitorStrength.STRONG},
    "duloxetine":    {"CYP2D6": InhibitorStrength.MODERATE},
    "sertraline":    {"CYP2D6": InhibitorStrength.MODERATE},
    "terbinafine":   {"CYP2D6": InhibitorStrength.STRONG},
    "diphenhydramine": {"CYP2D6": InhibitorStrength.MODERATE},
    "celecoxib":     {"CYP2D6": InhibitorStrength.WEAK},
    "cimetidine":    {"CYP2D6": InhibitorStrength.WEAK, "CYP2C19": InhibitorStrength.WEAK},
    "methadone":     {"CYP2D6": InhibitorStrength.WEAK},
    # ── CYP2C19 Inhibitors ──
    "fluvoxamine":   {"CYP2C19": InhibitorStrength.STRONG, "CYP1A2": InhibitorStrength.STRONG},
    "fluconazole":   {"CYP2C19": InhibitorStrength.MODERATE, "CYP2C9": InhibitorStrength.MODERATE},
    "omeprazole":    {"CYP2C19": InhibitorStrength.MODERATE},
    "esomeprazole":  {"CYP2C19": InhibitorStrength.MODERATE},
    "voriconazole":  {"CYP2C19": InhibitorStrength.STRONG},
    "ticlopidine":   {"CYP2C19": InhibitorStrength.STRONG, "CYP2B6": InhibitorStrength.STRONG},
    # ── CYP2C9 Inhibitors ──
    "amiodarone":    {"CYP2C9": InhibitorStrength.MODERATE},
    "miconazole":    {"CYP2C9": InhibitorStrength.STRONG},
    "metronidazole": {"CYP2C9": InhibitorStrength.WEAK},
    "trimethoprim":  {"CYP2C9": InhibitorStrength.WEAK},
    "sulfamethoxazole": {"CYP2C9": InhibitorStrength.WEAK},
    # ── CYP3A4/5 Inhibitors ──
    "ketoconazole":    {"CYP3A5": InhibitorStrength.STRONG, "CYP3A4": InhibitorStrength.STRONG},
    "itraconazole":    {"CYP3A5": InhibitorStrength.STRONG, "CYP3A4": InhibitorStrength.STRONG},
    "clarithromycin":  {"CYP3A5": InhibitorStrength.STRONG, "CYP3A4": InhibitorStrength.STRONG},
    "erythromycin":    {"CYP3A5": InhibitorStrength.MODERATE, "CYP3A4": InhibitorStrength.MODERATE},
    "diltiazem":       {"CYP3A5": InhibitorStrength.MODERATE, "CYP3A4": InhibitorStrength.MODERATE},
    "verapamil":       {"CYP3A5": InhibitorStrength.MODERATE, "CYP3A4": InhibitorStrength.MODERATE},
    "grapefruit_juice": {"CYP3A5": InhibitorStrength.MODERATE, "CYP3A4": InhibitorStrength.MODERATE},
    "grapefruit juice": {"CYP3A5": InhibitorStrength.MODERATE, "CYP3A4": InhibitorStrength.MODERATE},
    # ── CYP2B6 Inhibitors ──
    "clopidogrel":   {"CYP2B6": InhibitorStrength.MODERATE},
    "ritonavir":     {"CYP2B6": InhibitorStrength.STRONG, "CYP3A4": InhibitorStrength.STRONG, "CYP3A5": InhibitorStrength.STRONG},
    # ── CYP1A2 Inhibitors ──
    "ciprofloxacin": {"CYP1A2": InhibitorStrength.STRONG},
    "enoxacin":      {"CYP1A2": InhibitorStrength.STRONG},
    "mexiletine":    {"CYP1A2": InhibitorStrength.MODERATE},
}

# Drug → CYP enzyme → inducer strength
CYP_INDUCERS: Dict[str, Dict[str, InducerStrength]] = {
    "rifampin":       {"CYP2C19": InducerStrength.STRONG, "CYP3A5": InducerStrength.STRONG,
                       "CYP2C9": InducerStrength.STRONG, "CYP3A4": InducerStrength.STRONG,
                       "CYP2B6": InducerStrength.STRONG},
    "carbamazepine":  {"CYP2C19": InducerStrength.MODERATE, "CYP3A5": InducerStrength.STRONG,
                       "CYP3A4": InducerStrength.STRONG, "CYP2B6": InducerStrength.MODERATE},
    "phenytoin":      {"CYP2C19": InducerStrength.MODERATE, "CYP3A5": InducerStrength.STRONG,
                       "CYP3A4": InducerStrength.STRONG, "CYP2C9": InducerStrength.MODERATE,
                       "CYP2B6": InducerStrength.MODERATE},
    "phenobarbital":  {"CYP2C19": InducerStrength.MODERATE, "CYP3A5": InducerStrength.STRONG,
                       "CYP3A4": InducerStrength.STRONG, "CYP2B6": InducerStrength.MODERATE},
    "st_johns_wort":  {"CYP3A5": InducerStrength.STRONG, "CYP3A4": InducerStrength.STRONG,
                       "CYP2C19": InducerStrength.MODERATE},
    "efavirenz":      {"CYP3A5": InducerStrength.MODERATE, "CYP3A4": InducerStrength.MODERATE,
                       "CYP2C19": InducerStrength.MODERATE, "CYP2B6": InducerStrength.MODERATE},
    "modafinil":      {"CYP3A5": InducerStrength.WEAK, "CYP3A4": InducerStrength.WEAK},
    "dexamethasone":  {"CYP3A5": InducerStrength.WEAK, "CYP3A4": InducerStrength.WEAK},
    # ── CYP1A2 Inducers ──
    "tobacco_smoking": {"CYP1A2": InducerStrength.STRONG},
    "charbroiled_food": {"CYP1A2": InducerStrength.MODERATE},
}

# Substrates: gene → list of drugs metabolized by that enzyme
# Used to identify which co-prescribed drugs are affected by the phenoconversion
CYP_SUBSTRATES: Dict[str, List[str]] = {
    "CYP2D6": [
        "codeine", "tramadol", "oxycodone", "hydrocodone", "tamoxifen",
        "amitriptyline", "nortriptyline", "desipramine", "venlafaxine",
        "paroxetine", "fluoxetine", "aripiprazole",
        "haloperidol", "risperidone", "atomoxetine", "metoprolol",
        "carvedilol", "propafenone", "flecainide", "dextromethorphan",
        "ondansetron",
    ],
    "CYP2C19": [
        "clopidogrel", "omeprazole", "esomeprazole", "lansoprazole",
        "pantoprazole", "voriconazole", "diazepam", "citalopram",
        "escitalopram", "amitriptyline", "clomipramine", "sertraline",
        "phenytoin", "brivaracetam",
    ],
    "CYP2C9": [
        "warfarin", "phenytoin", "losartan", "irbesartan", "celecoxib",
        "flurbiprofen", "ibuprofen", "glipizide", "tolbutamide",
        "siponimod",
    ],
    "CYP3A5": [
        "tacrolimus", "sirolimus", "cyclosporine", "midazolam",
        "triazolam", "alprazolam",
    ],
    "CYP3A4": [
        "tacrolimus", "sirolimus", "cyclosporine", "midazolam",
        "triazolam", "alprazolam", "atorvastatin", "lovastatin",
        "simvastatin", "amlodipine", "felodipine", "nifedipine",
        "apixaban", "rivaroxaban",
    ],
    "CYP1A2": [
        "clozapine", "olanzapine", "theophylline", "caffeine",
        "melatonin", "duloxetine", "fluvoxamine", "ropinirole",
        "ramelteon", "tizanidine", "alosetron", "pirfenidone",
        "tasimelteon",
    ],
    "CYP2B6": [
        "efavirenz", "nevirapine", "methadone", "bupropion",
        "cyclophosphamide", "ifosfamide", "ketamine", "propofol",
        "artemisinin", "sertraline", "prasugrel",
    ],
}

# ── Phenotype shift rules ───────────────────────────────────────────────
# Maps (genetic_phenotype, inhibitor_strength) → effective_phenotype
# Based on CPIC phenoconversion tables.
INHIBITION_SHIFT: Dict[str, Dict[str, str]] = {
    "Ultrarapid Metabolizer": {
        InhibitorStrength.STRONG:   "Normal Metabolizer",
        InhibitorStrength.MODERATE: "Normal Metabolizer",
        InhibitorStrength.WEAK:     "Ultrarapid Metabolizer",
    },
    "Normal Metabolizer": {
        InhibitorStrength.STRONG:   "Poor Metabolizer",
        InhibitorStrength.MODERATE: "Intermediate Metabolizer",
        InhibitorStrength.WEAK:     "Normal Metabolizer",
    },
    "Intermediate Metabolizer": {
        InhibitorStrength.STRONG:   "Poor Metabolizer",
        InhibitorStrength.MODERATE: "Poor Metabolizer",
        InhibitorStrength.WEAK:     "Intermediate Metabolizer",
    },
    "Poor Metabolizer": {
        InhibitorStrength.STRONG:   "Poor Metabolizer",
        InhibitorStrength.MODERATE: "Poor Metabolizer",
        InhibitorStrength.WEAK:     "Poor Metabolizer",
    },
    # Enzyme-specific phenotypes — map to generic for shifting
    "CYP3A5 Expresser": {
        InhibitorStrength.STRONG:   "CYP3A5 Non-expresser",
        InhibitorStrength.MODERATE: "CYP3A5 Intermediate Expresser",
        InhibitorStrength.WEAK:     "CYP3A5 Expresser",
    },
    "CYP3A5 Intermediate Expresser": {
        InhibitorStrength.STRONG:   "CYP3A5 Non-expresser",
        InhibitorStrength.MODERATE: "CYP3A5 Non-expresser",
        InhibitorStrength.WEAK:     "CYP3A5 Intermediate Expresser",
    },
    "CYP3A5 Non-expresser": {
        InhibitorStrength.STRONG:   "CYP3A5 Non-expresser",
        InhibitorStrength.MODERATE: "CYP3A5 Non-expresser",
        InhibitorStrength.WEAK:     "CYP3A5 Non-expresser",
    },
}

# Maps (genetic_phenotype, inducer_strength) → effective_phenotype
INDUCTION_SHIFT: Dict[str, Dict[str, str]] = {
    "Poor Metabolizer": {
        InducerStrength.STRONG:   "Intermediate Metabolizer",
        InducerStrength.MODERATE: "Poor Metabolizer",
        InducerStrength.WEAK:     "Poor Metabolizer",
    },
    "Intermediate Metabolizer": {
        InducerStrength.STRONG:   "Normal Metabolizer",
        InducerStrength.MODERATE: "Normal Metabolizer",
        InducerStrength.WEAK:     "Intermediate Metabolizer",
    },
    "Normal Metabolizer": {
        InducerStrength.STRONG:   "Ultrarapid Metabolizer",
        InducerStrength.MODERATE: "Normal Metabolizer",
        InducerStrength.WEAK:     "Normal Metabolizer",
    },
    "Ultrarapid Metabolizer": {
        InducerStrength.STRONG:   "Ultrarapid Metabolizer",
        InducerStrength.MODERATE: "Ultrarapid Metabolizer",
        InducerStrength.WEAK:     "Ultrarapid Metabolizer",
    },
    "CYP3A5 Non-expresser": {
        InducerStrength.STRONG:   "CYP3A5 Intermediate Expresser",
        InducerStrength.MODERATE: "CYP3A5 Non-expresser",
        InducerStrength.WEAK:     "CYP3A5 Non-expresser",
    },
    "CYP3A5 Intermediate Expresser": {
        InducerStrength.STRONG:   "CYP3A5 Expresser",
        InducerStrength.MODERATE: "CYP3A5 Intermediate Expresser",
        InducerStrength.WEAK:     "CYP3A5 Intermediate Expresser",
    },
    "CYP3A5 Expresser": {
        InducerStrength.STRONG:   "CYP3A5 Expresser",
        InducerStrength.MODERATE: "CYP3A5 Expresser",
        InducerStrength.WEAK:     "CYP3A5 Expresser",
    },
}


# ═════════════════════════════════════════════════════════════════════════════
# Data class
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class PhenoconversionAlert:
    """Alert for a detected phenoconversion event."""
    gene: str
    genetic_phenotype: str
    effective_phenotype: str
    precipitant_drug: str
    inhibitor_strength: str  # InhibitorStrength or InducerStrength value
    affected_substrates: List[str] = field(default_factory=list)
    clinical_significance: str = ""


# ═════════════════════════════════════════════════════════════════════════════
# Detector
# ═════════════════════════════════════════════════════════════════════════════


class PhenoconversionDetector:
    """Detects when concomitant medications alter effective metabolizer status.

    A patient who is genetically a CYP2D6 Normal Metabolizer but is taking
    fluoxetine (strong CYP2D6 inhibitor) should be treated as a CYP2D6 Poor
    Metabolizer for all CYP2D6 substrate dosing decisions.

    This detector scans the medication list for known CYP inhibitors and
    inducers, then adjusts the genetic phenotype accordingly.
    """

    def detect(
        self,
        genetic_phenotype: Dict[str, Dict[str, Any]],
        medication_list: List[str],
    ) -> List[PhenoconversionAlert]:
        """Detect all phenoconversion events for a patient.

        Parameters
        ----------
        genetic_phenotype : dict
            gene → {phenotype, diplotype, ...} from PhenotypeTranslator.
        medication_list : list of str
            Currently prescribed medications (lowercase).

        Returns
        -------
        list of PhenoconversionAlert
        """
        alerts: List[PhenoconversionAlert] = []
        med_lower = [m.lower().strip() for m in medication_list]

        for med in med_lower:
            # Check inhibitors
            if med in CYP_INHIBITORS:
                for gene, strength in CYP_INHIBITORS[med].items():
                    gene_data = genetic_phenotype.get(gene)
                    if not gene_data:
                        continue

                    gen_pheno = gene_data.get("phenotype", "Normal Metabolizer")
                    eff_pheno = self._apply_inhibition(gen_pheno, strength)

                    if eff_pheno != gen_pheno:
                        affected = self._find_affected_substrates(gene, med_lower, med)
                        significance = self._assess_significance(
                            gene, gen_pheno, eff_pheno, strength, affected,
                        )
                        alerts.append(PhenoconversionAlert(
                            gene=gene,
                            genetic_phenotype=gen_pheno,
                            effective_phenotype=eff_pheno,
                            precipitant_drug=med,
                            inhibitor_strength=strength.value,
                            affected_substrates=affected,
                            clinical_significance=significance,
                        ))

            # Check inducers
            if med in CYP_INDUCERS:
                for gene, strength in CYP_INDUCERS[med].items():
                    gene_data = genetic_phenotype.get(gene)
                    if not gene_data:
                        continue

                    gen_pheno = gene_data.get("phenotype", "Normal Metabolizer")
                    eff_pheno = self._apply_induction(gen_pheno, strength)

                    if eff_pheno != gen_pheno:
                        affected = self._find_affected_substrates(gene, med_lower, med)
                        significance = self._assess_significance(
                            gene, gen_pheno, eff_pheno, strength.value, affected,
                        )
                        alerts.append(PhenoconversionAlert(
                            gene=gene,
                            genetic_phenotype=gen_pheno,
                            effective_phenotype=eff_pheno,
                            precipitant_drug=med,
                            inhibitor_strength=strength.value,
                            affected_substrates=affected,
                            clinical_significance=significance,
                        ))

        logger.info(
            "Phenoconversion scan: %d medications, %d conversions detected",
            len(medication_list), len(alerts),
        )
        return alerts

    def get_effective_phenotype(
        self,
        gene: str,
        genetic_phenotype: str,
        medication_list: List[str],
    ) -> str:
        """Return the effective phenotype for a single gene after phenoconversion.

        Applies the strongest inhibitor/inducer effect found in the
        medication list.  Inhibition takes precedence over induction when
        both are present (conservative approach per CPIC guidance).

        Parameters
        ----------
        gene : str
            Pharmacogene (e.g. "CYP2D6").
        genetic_phenotype : str
            Genetic phenotype (e.g. "Normal Metabolizer").
        medication_list : list of str
            Currently prescribed medications.

        Returns
        -------
        str
            Effective phenotype after considering DDIs.
        """
        med_lower = [m.lower().strip() for m in medication_list]
        strongest_inhibitor: Optional[InhibitorStrength] = None
        strongest_inducer: Optional[InducerStrength] = None

        strength_order_inh = {
            InhibitorStrength.WEAK: 0,
            InhibitorStrength.MODERATE: 1,
            InhibitorStrength.STRONG: 2,
        }
        strength_order_ind = {
            InducerStrength.WEAK: 0,
            InducerStrength.MODERATE: 1,
            InducerStrength.STRONG: 2,
        }

        for med in med_lower:
            if med in CYP_INHIBITORS:
                inh = CYP_INHIBITORS[med].get(gene)
                if inh:
                    if strongest_inhibitor is None or (
                        strength_order_inh[inh] > strength_order_inh[strongest_inhibitor]
                    ):
                        strongest_inhibitor = inh

            if med in CYP_INDUCERS:
                ind = CYP_INDUCERS[med].get(gene)
                if ind:
                    if strongest_inducer is None or (
                        strength_order_ind[ind] > strength_order_ind[strongest_inducer]
                    ):
                        strongest_inducer = ind

        # Inhibition takes precedence (conservative)
        if strongest_inhibitor is not None:
            return self._apply_inhibition(genetic_phenotype, strongest_inhibitor)
        if strongest_inducer is not None:
            return self._apply_induction(genetic_phenotype, strongest_inducer)

        return genetic_phenotype

    def get_all_effective_phenotypes(
        self,
        genetic_profile: Dict[str, Dict[str, Any]],
        medication_list: List[str],
    ) -> Dict[str, Dict[str, Any]]:
        """Compute effective phenotypes for all genes in the profile.

        Returns a copy of the genetic profile with effective_phenotype and
        phenoconverted (bool) fields added.

        Parameters
        ----------
        genetic_profile : dict
            gene → {phenotype, diplotype, ...} from PhenotypeTranslator.
        medication_list : list of str

        Returns
        -------
        dict
            gene → {phenotype, effective_phenotype, phenoconverted, ...}
        """
        result: Dict[str, Dict[str, Any]] = {}

        for gene, data in genetic_profile.items():
            gen_pheno = data.get("phenotype", "Normal Metabolizer")
            eff_pheno = self.get_effective_phenotype(gene, gen_pheno, medication_list)
            entry = dict(data)
            entry["effective_phenotype"] = eff_pheno
            entry["phenoconverted"] = eff_pheno != gen_pheno
            result[gene] = entry

        converted = sum(1 for v in result.values() if v.get("phenoconverted"))
        logger.info(
            "Effective phenotypes computed for %d genes (%d phenoconverted)",
            len(result), converted,
        )
        return result

    # ── Private helpers ─────────────────────────────────────────────────

    @staticmethod
    def _apply_inhibition(phenotype: str, strength: InhibitorStrength) -> str:
        """Shift phenotype due to enzyme inhibition."""
        shift_map = INHIBITION_SHIFT.get(phenotype)
        if not shift_map:
            return phenotype
        return shift_map.get(strength, phenotype)

    @staticmethod
    def _apply_induction(phenotype: str, strength: InducerStrength) -> str:
        """Shift phenotype due to enzyme induction."""
        shift_map = INDUCTION_SHIFT.get(phenotype)
        if not shift_map:
            return phenotype
        return shift_map.get(strength, phenotype)

    @staticmethod
    def _find_affected_substrates(
        gene: str, medication_list: List[str], precipitant: str,
    ) -> List[str]:
        """Identify co-prescribed substrates affected by phenoconversion."""
        substrates = CYP_SUBSTRATES.get(gene, [])
        affected = [
            m for m in medication_list
            if m in substrates and m != precipitant
        ]
        return affected

    @staticmethod
    def _assess_significance(
        gene: str,
        genetic_phenotype: str,
        effective_phenotype: str,
        strength: Any,
        affected_substrates: List[str],
    ) -> str:
        """Generate a clinical significance statement."""
        strength_str = strength.value if hasattr(strength, "value") else str(strength)

        if not affected_substrates:
            return (
                f"Phenoconversion detected: {gene} shifted from "
                f"{genetic_phenotype} to {effective_phenotype} due to "
                f"{strength_str} inhibition/induction. No co-prescribed "
                f"{gene} substrates identified — monitor if new substrates added."
            )

        substrate_str = ", ".join(affected_substrates)
        return (
            f"CLINICAL ACTION REQUIRED: {gene} phenoconversion from "
            f"{genetic_phenotype} to {effective_phenotype} ({strength_str} "
            f"effect). Co-prescribed {gene} substrates affected: "
            f"{substrate_str}. Review dosing for these medications based on "
            f"effective phenotype ({effective_phenotype})."
        )
