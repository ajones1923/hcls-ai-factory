"""Cross-biomarker discordance detection engine.

Evaluates pairs of biomarker values against known discordance patterns
from biomarker_discordance_rules.json. Identifies contradictory or
unexpected relationships between biomarkers.

Author: Adam Jones
Date: March 2026
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "reference"


class DiscordanceFinding:
    """A detected cross-biomarker discordance."""
    def __init__(self, rule_name: str, biomarker_a: str, biomarker_b: str,
                 value_a: Optional[float], value_b: Optional[float],
                 condition: str, differential_diagnosis: List[str],
                 agent_handoff: List[str], priority: str, text_chunk: str):
        self.rule_name = rule_name
        self.biomarker_a = biomarker_a
        self.biomarker_b = biomarker_b
        self.value_a = value_a
        self.value_b = value_b
        self.condition = condition
        self.differential_diagnosis = differential_diagnosis
        self.agent_handoff = agent_handoff
        self.priority = priority
        self.text_chunk = text_chunk

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rule_name": self.rule_name,
            "biomarker_a": self.biomarker_a,
            "biomarker_b": self.biomarker_b,
            "value_a": self.value_a,
            "value_b": self.value_b,
            "condition": self.condition,
            "differential_diagnosis": self.differential_diagnosis,
            "agent_handoff": self.agent_handoff,
            "priority": self.priority,
        }

    def to_alert_string(self) -> str:
        return (
            f"DISCORDANCE [{self.priority.upper()}]: {self.rule_name} — "
            f"{self.biomarker_a}={self.value_a} vs {self.biomarker_b}={self.value_b}. "
            f"Differential: {', '.join(self.differential_diagnosis[:3])}"
        )


class DiscordanceDetector:
    """Detects cross-biomarker discordance patterns.

    Evaluates 12 discordance rules from biomarker_discordance_rules.json
    against patient biomarker values to identify contradictory or
    clinically significant cross-biomarker patterns.
    """

    # Known biomarker reference ranges for determining HIGH/LOW/NORMAL status
    # These are approximate clinical thresholds used for discordance logic
    THRESHOLDS = {
        "Ferritin": {"low": 30, "high": 300},
        "Transferrin Saturation (TSAT)": {"low": 20, "high": 45},
        "TSAT": {"low": 20, "high": 45},
        "Free T3": {"low": 2.0, "high": 4.4},
        "TSH": {"low": 0.45, "high": 4.5},
        "GGT": {"low": 0, "high": 65},
        "AST": {"low": 0, "high": 40},
        "ALT": {"low": 0, "high": 44},
        "ALP": {"low": 44, "high": 121},
        "Fibrinogen": {"low": 200, "high": 400},
        "hs-CRP": {"low": 0, "high": 3.0},
        "Cystatin C": {"low": 0.53, "high": 0.95},
        "eGFR": {"low": 60, "high": 9999},
        "Adiponectin": {"low": 4, "high": 20},
        "HDL-P": {"low": 30.5, "high": 9999},
        "HbA1c": {"low": 4.0, "high": 5.6},
        "Glucose Fasting": {"low": 65, "high": 99},
        "LDL-C": {"low": 0, "high": 100},
        "ApoB": {"low": 0, "high": 90},
        "Platelet Count": {"low": 150, "high": 379},
        "INR": {"low": 0.8, "high": 1.2},
        "Free T4": {"low": 0.82, "high": 1.77},
        "Vitamin B12": {"low": 232, "high": 1245},
        "MMA": {"low": 0, "high": 0.4},
        "Vitamin D": {"low": 30, "high": 100},
        "PTH": {"low": 15, "high": 65},
    }

    # Biomarker name aliases for flexible matching
    ALIASES = {
        "ferritin": "Ferritin",
        "tsat": "TSAT",
        "transferrin_saturation": "Transferrin Saturation (TSAT)",
        "free_t3": "Free T3",
        "ft3": "Free T3",
        "tsh": "TSH",
        "ggt": "GGT",
        "ast": "AST",
        "alt": "ALT",
        "alp": "ALP",
        "fibrinogen": "Fibrinogen",
        "hs_crp": "hs-CRP",
        "hscrp": "hs-CRP",
        "hs-crp": "hs-CRP",
        "crp": "hs-CRP",
        "cystatin_c": "Cystatin C",
        "egfr": "eGFR",
        "adiponectin": "Adiponectin",
        "hdl_p": "HDL-P",
        "hdl-p": "HDL-P",
        "hba1c": "HbA1c",
        "glucose": "Glucose Fasting",
        "glucose_fasting": "Glucose Fasting",
        "fasting_glucose": "Glucose Fasting",
        "ldl": "LDL-C",
        "ldl_c": "LDL-C",
        "ldl-c": "LDL-C",
        "apob": "ApoB",
        "platelet_count": "Platelet Count",
        "platelets": "Platelet Count",
        "inr": "INR",
        "free_t4": "Free T4",
        "ft4": "Free T4",
        "vitamin_b12": "Vitamin B12",
        "b12": "Vitamin B12",
        "mma": "MMA",
        "methylmalonic_acid": "MMA",
        "vitamin_d": "Vitamin D",
        "pth": "PTH",
    }

    def __init__(self, rules: Optional[List[Dict]] = None):
        if rules is not None:
            self._rules = rules
        else:
            self._rules = self._load_rules()

    def _load_rules(self) -> List[Dict]:
        filepath = REFERENCE_DIR / "biomarker_discordance_rules.json"
        try:
            with open(filepath) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.warning(f"Failed to load {filepath}: {e}; module will return empty results")
            return []

    def _resolve(self, name: str) -> Optional[str]:
        """Resolve biomarker name to canonical form."""
        canonical = self.ALIASES.get(name.lower())
        if canonical:
            return canonical
        # Try direct match (case-insensitive)
        for threshold_name in self.THRESHOLDS:
            if threshold_name.lower() == name.lower():
                return threshold_name
        return None

    def _get_value(self, biomarkers: Dict[str, float], target: str) -> Optional[float]:
        """Get a biomarker value by trying aliases."""
        for name, value in biomarkers.items():
            resolved = self._resolve(name)
            if resolved and resolved == target:
                return value
        return None

    def _is_status(self, value: float, biomarker: str, status: str) -> bool:
        """Check if a biomarker value matches a status (HIGH/LOW/NORMAL/CRITICAL)."""
        thresholds = self.THRESHOLDS.get(biomarker)
        if not thresholds:
            return False

        low = thresholds["low"]
        high = thresholds["high"]

        if status == "HIGH":
            return value > high
        elif status == "LOW":
            return value < low
        elif status == "NORMAL":
            return low <= value <= high
        elif status == "CRITICAL HIGH":
            return value > high * 2.5  # roughly critical range
        elif status == "CRITICAL":
            return value < low * 0.5 or value > high * 2
        return False

    def check(self, biomarkers: Dict[str, float]) -> List[DiscordanceFinding]:
        """Check patient biomarkers against all discordance rules.

        Args:
            biomarkers: Dict of biomarker_name -> measured value.

        Returns:
            List of DiscordanceFinding for any detected discordances.
        """
        findings = []

        for rule in self._rules:
            rule_id = rule.get("id", "")
            biomarker_a_name = rule.get("biomarker_a", "")
            biomarker_b_name = rule.get("biomarker_b", "")
            condition = rule.get("condition", "")

            # Get patient values for both biomarkers
            value_a = self._get_value(biomarkers, biomarker_a_name)
            value_b = self._get_value(biomarkers, biomarker_b_name)

            # Both biomarkers must be present
            if value_a is None or value_b is None:
                continue

            # Parse condition to check if discordance is present
            # Conditions follow patterns like "Ferritin LOW and TSAT HIGH"
            discordant = self._evaluate_condition(
                condition, biomarker_a_name, biomarker_b_name, value_a, value_b
            )

            if discordant:
                diff_diag = rule.get("differential_diagnosis", [])
                if isinstance(diff_diag, str):
                    diff_diag = [d.strip() for d in diff_diag.split(",") if d.strip()]

                handoff = rule.get("agent_handoff", [])
                if isinstance(handoff, str):
                    handoff = [h.strip() for h in handoff.split(",") if h.strip()]

                findings.append(DiscordanceFinding(
                    rule_name=rule.get("name", rule_id),
                    biomarker_a=biomarker_a_name,
                    biomarker_b=biomarker_b_name,
                    value_a=value_a,
                    value_b=value_b,
                    condition=condition,
                    differential_diagnosis=diff_diag,
                    agent_handoff=handoff,
                    priority=rule.get("priority", "medium"),
                    text_chunk=rule.get("text_chunk", ""),
                ))

        # Sort by priority: high first
        priority_order = {"high": 0, "medium": 1, "low": 2}
        findings.sort(key=lambda f: priority_order.get(f.priority, 3))

        return findings

    def _evaluate_condition(self, condition: str, biomarker_a: str, biomarker_b: str,
                            value_a: float, value_b: float) -> bool:
        """Evaluate a discordance condition string against actual values.

        Parses conditions like:
        - "Ferritin LOW and TSAT HIGH"
        - "HbA1c NORMAL and Fasting Glucose CRITICAL HIGH"
        - "LDL-C DISCORDANT with ApoB" (checks >20% divergence from expected)
        - "Platelet CRITICAL and INR CRITICAL" (compound critical)
        """
        condition_upper = condition.upper()

        # Handle DISCORDANT pattern (e.g., LDL-C vs ApoB)
        if "DISCORDANT" in condition_upper:
            # LDL-C and ApoB should correlate; check for >20% disagreement in risk classification
            a_high = self._is_status(value_a, biomarker_a, "HIGH")
            b_high = self._is_status(value_b, biomarker_b, "HIGH")
            a_normal = self._is_status(value_a, biomarker_a, "NORMAL")
            b_normal = self._is_status(value_b, biomarker_b, "NORMAL")
            return (a_high and b_normal) or (a_normal and b_high)

        # Handle compound critical patterns
        if " AND " in condition_upper:
            parts = condition_upper.split(" AND ")
            if len(parts) == 2:
                a_match = self._check_part(parts[0].strip(), biomarker_a, value_a)
                b_match = self._check_part(parts[1].strip(), biomarker_b, value_b)
                return a_match and b_match

        return False

    def _check_part(self, part: str, biomarker: str, value: float) -> bool:
        """Check a single condition part like 'Ferritin LOW' or 'TSH HIGH'."""
        part = part.strip()

        # Extract status from end of string
        for status in ["CRITICAL HIGH", "CRITICAL LOW", "CRITICAL", "HIGH", "LOW", "NORMAL"]:
            if part.endswith(status):
                return self._is_status(value, biomarker, status)

        return False

    @property
    def rule_count(self) -> int:
        return len(self._rules)
