"""Multi-lab reference range interpreter.

Compares patient biomarker values against Quest Diagnostics, LabCorp,
and Function Health optimal ranges to provide nuanced interpretation
that goes beyond standard lab reporting.

Author: Adam Jones
Date: March 2026
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "reference"


class RangeComparison:
    """Comparison of a biomarker value against multiple lab ranges."""
    def __init__(self, biomarker: str, value: float, unit: str,
                 quest_status: str, labcorp_status: str, optimal_status: str,
                 quest_range: Optional[Dict] = None,
                 labcorp_range: Optional[Dict] = None,
                 optimal_range: Optional[Dict] = None):
        self.biomarker = biomarker
        self.value = value
        self.unit = unit
        self.quest_status = quest_status      # "normal", "low", "high"
        self.labcorp_status = labcorp_status
        self.optimal_status = optimal_status
        self.quest_range = quest_range
        self.labcorp_range = labcorp_range
        self.optimal_range = optimal_range

    @property
    def has_discrepancy(self) -> bool:
        """True if the value is normal by standard labs but not by optimal ranges."""
        return (self.quest_status == "normal" or self.labcorp_status == "normal") and self.optimal_status != "normal"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "biomarker": self.biomarker,
            "value": self.value,
            "unit": self.unit,
            "quest_status": self.quest_status,
            "labcorp_status": self.labcorp_status,
            "optimal_status": self.optimal_status,
            "has_discrepancy": self.has_discrepancy,
        }

    def to_interpretation(self) -> str:
        """Generate human-readable interpretation."""
        if not self.has_discrepancy:
            if self.optimal_status == "normal":
                return f"{self.biomarker}: {self.value} {self.unit} — within optimal range."
            else:
                return (
                    f"{self.biomarker}: {self.value} {self.unit} — "
                    f"flagged as {self.quest_status} by Quest, {self.labcorp_status} by LabCorp."
                )

        # Discrepancy: normal by standard labs, suboptimal by Function Health
        opt = self.optimal_range or {}
        opt_min = opt.get("min", "?")
        opt_max = opt.get("max", "?")
        return (
            f"{self.biomarker}: {self.value} {self.unit} — "
            f"within standard reference range but {self.optimal_status} by Function Health optimal "
            f"(optimal: {opt_min}-{opt_max} {self.unit}). Consider optimization."
        )


class LabRangeInterpreter:
    """Interprets biomarker values against Quest, LabCorp, and Function Health ranges.

    Loads multi-lab ranges from biomarker_lab_ranges.json and provides
    three-way comparison for each biomarker value.
    """

    # Mapping aliases for flexible biomarker name matching
    ALIASES = {
        "glucose": "Glucose Fasting", "fasting_glucose": "Glucose Fasting", "glucose_fasting": "Glucose Fasting",
        "hba1c": "HbA1c", "a1c": "HbA1c",
        "insulin": "Insulin Fasting", "fasting_insulin": "Insulin Fasting", "insulin_fasting": "Insulin Fasting",
        "homa_ir": "HOMA-IR", "homa-ir": "HOMA-IR", "homair": "HOMA-IR",
        "ldl": "LDL Cholesterol", "ldl_c": "LDL Cholesterol", "ldl_cholesterol": "LDL Cholesterol",
        "hdl": "HDL Cholesterol", "hdl_cholesterol": "HDL Cholesterol",
        "triglycerides": "Triglycerides", "tg": "Triglycerides",
        "total_cholesterol": "Total Cholesterol",
        "apob": "ApoB", "apo_b": "ApoB",
        "hscrp": "hs-CRP", "hs_crp": "hs-CRP", "hs-crp": "hs-CRP", "crp": "hs-CRP",
        "homocysteine": "Homocysteine",
        "ferritin": "Ferritin",
        "tsat": "Transferrin Saturation", "transferrin_saturation": "Transferrin Saturation",
        "vitamin_d": "Vitamin D 25-OH", "vit_d": "Vitamin D 25-OH",
        "b12": "Vitamin B12", "vitamin_b12": "Vitamin B12",
        "folate": "Folate",
        "tsh": "TSH",
        "free_t4": "Free T4", "ft4": "Free T4",
        "free_t3": "Free T3", "ft3": "Free T3",
        "creatinine": "Creatinine",
        "egfr": "eGFR", "gfr": "eGFR",
        "cystatin_c": "Cystatin C",
        "uric_acid": "Uric Acid",
        "alt": "ALT", "ast": "AST", "ggt": "GGT",
        "albumin": "Albumin",
        "omega3": "Omega-3 Index", "omega_3": "Omega-3 Index", "omega_3_index": "Omega-3 Index",
        "magnesium": "Magnesium", "zinc": "Zinc", "selenium": "Selenium",
        "sodium": "Sodium", "potassium": "Potassium",
    }

    def __init__(self, data: Optional[Dict] = None):
        if data is not None:
            self._data = data
        else:
            self._data = self._load_data()

    def _load_data(self) -> Dict:
        filepath = REFERENCE_DIR / "biomarker_lab_ranges.json"
        try:
            with open(filepath) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.warning(f"Failed to load {filepath}: {e}; module will return empty results")
            return {"labs": {}}

    def _resolve(self, name: str) -> Optional[str]:
        """Resolve biomarker name to canonical form."""
        canonical = self.ALIASES.get(name.lower())
        if canonical:
            return canonical
        # Direct match check
        labs = self._data.get("labs", {})
        for lab_data in labs.values():
            if name in lab_data.get("ranges", {}):
                return name
        return None

    @staticmethod
    def _check_range(value: float, range_info: Dict) -> str:
        """Check a value against a range dict with min/max."""
        low = range_info.get("min")
        high = range_info.get("max")
        if low is not None and value < low:
            return "low"
        if high is not None and value > high:
            return "high"
        return "normal"

    def interpret(self, biomarkers: Dict[str, float],
                  sex: str = "Male") -> List[RangeComparison]:
        """Compare biomarker values against all three lab ranges.

        Args:
            biomarkers: Dict of biomarker_name -> measured value.
            sex: Patient sex for sex-specific ranges (e.g., HDL, Hemoglobin).

        Returns:
            List of RangeComparison objects.
        """
        labs = self._data.get("labs", {})
        quest = labs.get("quest_diagnostics", {}).get("ranges", {})
        labcorp = labs.get("labcorp", {}).get("ranges", {})
        optimal = labs.get("function_health_optimal", {}).get("ranges", {})

        results = []

        for name, value in biomarkers.items():
            canonical = self._resolve(name)
            if not canonical:
                continue

            # Handle sex-specific biomarkers: try "{name} ({sex})" first, fall back to "{name}"
            sex_specific = f"{canonical} ({sex})" if sex else None
            q_range = (quest.get(sex_specific) if sex_specific else None) or quest.get(canonical)
            l_range = (labcorp.get(sex_specific) if sex_specific else None) or labcorp.get(canonical)
            o_range = (optimal.get(sex_specific) if sex_specific else None) or optimal.get(canonical)

            if not any([q_range, l_range, o_range]):
                continue

            unit = (q_range or l_range or o_range or {}).get("unit", "")

            q_status = self._check_range(value, q_range) if q_range else "unknown"
            l_status = self._check_range(value, l_range) if l_range else "unknown"
            o_status = self._check_range(value, o_range) if o_range else "unknown"

            results.append(RangeComparison(
                biomarker=canonical,
                value=value,
                unit=unit,
                quest_status=q_status,
                labcorp_status=l_status,
                optimal_status=o_status,
                quest_range=q_range,
                labcorp_range=l_range,
                optimal_range=o_range,
            ))

        return results

    def get_discrepancies(self, biomarkers: Dict[str, float],
                          sex: str = "Male") -> List[RangeComparison]:
        """Return only biomarkers where standard and optimal ranges disagree."""
        return [r for r in self.interpret(biomarkers, sex) if r.has_discrepancy]

    def format_report_section(self, biomarkers: Dict[str, float],
                               sex: str = "Male") -> str:
        """Generate a markdown report section for lab range comparison."""
        discrepancies = self.get_discrepancies(biomarkers, sex)
        if not discrepancies:
            return ""

        lines = ["### Lab Range Optimization Opportunities\n"]
        lines.append("The following biomarkers are within standard lab reference ranges but fall outside Function Health optimal ranges:\n")
        for d in discrepancies:
            lines.append(f"- {d.to_interpretation()}")
        return "\n".join(lines)
