"""Critical biomarker value detection engine.

Loads thresholds from biomarker_critical_values.json and evaluates patient
biomarker values against critical/urgent/warning thresholds in real-time.

Author: Adam Jones
Date: March 2026
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# Path to reference data
REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "reference"


class CriticalValueAlert:
    """A triggered critical value alert."""
    def __init__(self, biomarker: str, value: float, threshold: float,
                 direction: str, severity: str, escalation_target: str,
                 clinical_action: str, cross_checks: List[str],
                 loinc_code: str = ""):
        self.biomarker = biomarker
        self.value = value
        self.threshold = threshold
        self.direction = direction  # "high" or "low"
        self.severity = severity
        self.escalation_target = escalation_target
        self.clinical_action = clinical_action
        self.cross_checks = cross_checks
        self.loinc_code = loinc_code

    def to_dict(self) -> Dict[str, Any]:
        return {
            "biomarker": self.biomarker,
            "value": self.value,
            "threshold": self.threshold,
            "direction": self.direction,
            "severity": self.severity,
            "escalation_target": self.escalation_target,
            "clinical_action": self.clinical_action,
            "cross_checks": self.cross_checks,
            "loinc_code": self.loinc_code,
        }

    def to_alert_string(self) -> str:
        dir_label = "above" if self.direction == "high" else "below"
        return (
            f"CRITICAL VALUE [{self.severity.upper()}]: {self.biomarker} = {self.value} "
            f"({dir_label} threshold {self.threshold}). "
            f"Escalate to {self.escalation_target}. "
            f"Action: {self.clinical_action}"
        )


class CriticalValueEngine:
    """Evaluates patient biomarker values against critical thresholds.

    Loads the 21 critical value rules from biomarker_critical_values.json
    and checks incoming biomarker values against them, generating alerts
    for any values that exceed critical/urgent/warning thresholds.
    """

    # Mapping of biomarker names in critical_values.json to common input names
    BIOMARKER_ALIASES = {
        "Platelet Count": ["platelet_count", "platelets", "plt"],
        "Glucose": ["glucose", "glucose_fasting", "fasting_glucose", "glucose fasting"],
        "Potassium": ["potassium", "k"],
        "INR": ["inr", "pt_inr"],
        "Sodium": ["sodium", "na"],
        "Hemoglobin": ["hemoglobin", "hgb", "hb"],
        "Calcium (Total)": ["calcium", "ca", "calcium_total", "total_calcium"],
        "Troponin I (High-Sensitivity)": ["troponin", "troponin_i", "tni", "hs_troponin"],
        "WBC Count": ["wbc", "wbc_count", "white_blood_cells"],
        "Creatinine": ["creatinine", "cr"],
        "Total Bilirubin": ["total_bilirubin", "bilirubin"],
        "eGFR (CKD-EPI)": ["egfr", "gfr", "egfr_ckd_epi"],
        "Free T4": ["free_t4", "ft4"],
        "TSH": ["tsh"],
        "Lactate": ["lactate", "lactic_acid"],
    }

    def __init__(self, rules: Optional[List[Dict]] = None):
        if rules is not None:
            self._rules = rules
        else:
            self._rules = self._load_rules()
        # Build reverse alias map: lowercase alias -> canonical biomarker name
        self._alias_map: Dict[str, str] = {}
        for canonical, aliases in self.BIOMARKER_ALIASES.items():
            self._alias_map[canonical.lower()] = canonical
            for alias in aliases:
                self._alias_map[alias.lower()] = canonical

    def _load_rules(self) -> List[Dict]:
        filepath = REFERENCE_DIR / "biomarker_critical_values.json"
        try:
            with open(filepath) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.warning(f"Failed to load {filepath}: {e}; module will return empty results")
            return []

    def _resolve_biomarker(self, name: str) -> Optional[str]:
        """Resolve a biomarker name to its canonical form."""
        return self._alias_map.get(name.lower())

    def check(self, biomarkers: Dict[str, float]) -> List[CriticalValueAlert]:
        """Check patient biomarkers against all critical value thresholds.

        Args:
            biomarkers: Dict of biomarker_name -> measured value.

        Returns:
            List of CriticalValueAlert for any triggered thresholds.
        """
        alerts = []

        for rule in self._rules:
            rule_biomarker = rule.get("biomarker", "")
            critical_low = rule.get("critical_low")
            critical_high = rule.get("critical_high")

            # Find matching patient biomarker
            patient_value = None
            for name, value in biomarkers.items():
                canonical = self._resolve_biomarker(name)
                if canonical and canonical == rule_biomarker:
                    patient_value = value
                    break

            if patient_value is None:
                continue

            cross_checks = rule.get("cross_checks", [])
            if isinstance(cross_checks, str):
                cross_checks = [c.strip() for c in cross_checks.split(",") if c.strip()]

            # Check high threshold
            if critical_high is not None and patient_value > critical_high:
                alerts.append(CriticalValueAlert(
                    biomarker=rule_biomarker,
                    value=patient_value,
                    threshold=critical_high,
                    direction="high",
                    severity=rule.get("severity", "warning"),
                    escalation_target=rule.get("escalation_target", ""),
                    clinical_action=rule.get("clinical_action", ""),
                    cross_checks=cross_checks,
                    loinc_code=rule.get("loinc_code", ""),
                ))

            # Check low threshold
            if critical_low is not None and patient_value < critical_low:
                alerts.append(CriticalValueAlert(
                    biomarker=rule_biomarker,
                    value=patient_value,
                    threshold=critical_low,
                    direction="low",
                    severity=rule.get("severity", "warning"),
                    escalation_target=rule.get("escalation_target", ""),
                    clinical_action=rule.get("clinical_action", ""),
                    cross_checks=cross_checks,
                    loinc_code=rule.get("loinc_code", ""),
                ))

        # Sort by severity: critical first, then urgent, then warning
        severity_order = {"critical": 0, "urgent": 1, "warning": 2}
        alerts.sort(key=lambda a: severity_order.get(a.severity, 3))

        return alerts

    @property
    def rule_count(self) -> int:
        return len(self._rules)
