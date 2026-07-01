"""Tests for CriticalValueEngine and CriticalValueAlert.

Validates critical value threshold checking, alias resolution,
severity sorting, alert formatting, and graceful degradation when
the reference data file is missing.

Author: Adam Jones
Date: March 2026
"""

import json
from unittest.mock import patch

import pytest

from src.critical_values import CriticalValueAlert, CriticalValueEngine

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def sample_rules():
    """Return a minimal set of critical value rules for testing."""
    return [
        {
            "id": "crit_glucose_high",
            "biomarker": "Glucose",
            "loinc_code": "2345-7",
            "critical_low": None,
            "critical_high": 250.0,
            "severity": "urgent",
            "escalation_target": "endocrinology_oncall",
            "clinical_action": "Check serum ketones and arterial blood gas.",
            "cross_checks": ["HbA1c", "Serum Ketones", "Bicarbonate"],
        },
        {
            "id": "crit_glucose_low",
            "biomarker": "Glucose",
            "loinc_code": "2345-7",
            "critical_low": 50.0,
            "critical_high": None,
            "severity": "critical",
            "escalation_target": "rapid_response_team",
            "clinical_action": "Administer IV dextrose (D50) immediately.",
            "cross_checks": ["Insulin", "C-Peptide"],
        },
        {
            "id": "crit_potassium_high",
            "biomarker": "Potassium",
            "loinc_code": "2823-3",
            "critical_low": None,
            "critical_high": 6.0,
            "severity": "critical",
            "escalation_target": "rapid_response_team",
            "clinical_action": "Obtain stat 12-lead ECG. Administer IV calcium gluconate.",
            "cross_checks": ["Creatinine", "eGFR"],
        },
        {
            "id": "crit_platelet_low",
            "biomarker": "Platelet Count",
            "loinc_code": "777-3",
            "critical_low": None,
            "critical_high": 100.0,
            "severity": "warning",
            "escalation_target": "ordering_provider",
            "clinical_action": "Review medication list for offending agents.",
            "cross_checks": ["INR", "Fibrinogen"],
        },
    ]


@pytest.fixture
def engine(sample_rules):
    """Return a CriticalValueEngine initialized with sample rules."""
    return CriticalValueEngine(rules=sample_rules)


# =====================================================================
# INITIALIZATION
# =====================================================================


class TestInitialization:
    """Tests for CriticalValueEngine construction and rule loading."""

    def test_init_with_provided_rules(self, sample_rules):
        """Engine uses rules passed directly to the constructor."""
        engine = CriticalValueEngine(rules=sample_rules)
        assert engine.rule_count == len(sample_rules)

    def test_init_loads_rules_from_json(self, sample_rules, tmp_path):
        """Engine loads rules from the JSON file when none are provided."""
        json_file = tmp_path / "biomarker_critical_values.json"
        json_file.write_text(json.dumps(sample_rules))

        with patch("src.critical_values.REFERENCE_DIR", tmp_path):
            engine = CriticalValueEngine()
            assert engine.rule_count == len(sample_rules)

    def test_init_builds_alias_map(self, engine):
        """Engine builds a reverse alias map from BIOMARKER_ALIASES."""
        # The alias map should resolve lowercase aliases to canonical names
        assert engine._resolve_biomarker("glucose") == "Glucose"
        assert engine._resolve_biomarker("potassium") == "Potassium"
        assert engine._resolve_biomarker("plt") == "Platelet Count"


# =====================================================================
# CHECK — HIGH THRESHOLD
# =====================================================================


class TestCheckHighThreshold:
    """Tests for values exceeding the critical_high threshold."""

    def test_value_above_critical_high_returns_alert(self, engine):
        """A glucose value above 250 should trigger an urgent alert."""
        alerts = engine.check({"glucose": 300.0})
        assert len(alerts) == 1
        alert = alerts[0]
        assert alert.biomarker == "Glucose"
        assert alert.value == 300.0
        assert alert.threshold == 250.0
        assert alert.direction == "high"
        assert alert.severity == "urgent"

    def test_value_at_threshold_no_alert(self, engine):
        """A value exactly at the threshold should not trigger an alert."""
        alerts = engine.check({"glucose": 250.0})
        # 250.0 is not > 250.0, so no alert
        glucose_alerts = [a for a in alerts if a.biomarker == "Glucose" and a.direction == "high"]
        assert len(glucose_alerts) == 0

    def test_potassium_above_critical_high(self, engine):
        """Potassium above 6.0 should trigger a critical alert."""
        alerts = engine.check({"potassium": 7.2})
        assert len(alerts) == 1
        assert alerts[0].severity == "critical"
        assert alerts[0].direction == "high"


# =====================================================================
# CHECK — LOW THRESHOLD
# =====================================================================


class TestCheckLowThreshold:
    """Tests for values below the critical_low threshold."""

    def test_value_below_critical_low_returns_alert(self, engine):
        """A glucose value below 50 should trigger a critical alert."""
        alerts = engine.check({"glucose": 35.0})
        # Should match the critical_low rule (severity=critical)
        low_alerts = [a for a in alerts if a.direction == "low"]
        assert len(low_alerts) == 1
        alert = low_alerts[0]
        assert alert.biomarker == "Glucose"
        assert alert.value == 35.0
        assert alert.threshold == 50.0
        assert alert.severity == "critical"

    def test_value_at_low_threshold_no_alert(self, engine):
        """A value exactly at the low threshold should not trigger an alert."""
        alerts = engine.check({"glucose": 50.0})
        low_alerts = [a for a in alerts if a.direction == "low"]
        assert len(low_alerts) == 0


# =====================================================================
# CHECK — NORMAL VALUES
# =====================================================================


class TestCheckNormalValues:
    """Tests for values within normal range (no alerts)."""

    def test_normal_glucose_no_alerts(self, engine):
        """A glucose of 95 should not trigger any alerts."""
        alerts = engine.check({"glucose": 95.0})
        assert len(alerts) == 0

    def test_normal_potassium_no_alerts(self, engine):
        """A potassium of 4.0 should not trigger any alerts."""
        alerts = engine.check({"potassium": 4.0})
        assert len(alerts) == 0

    def test_empty_biomarkers_no_alerts(self, engine):
        """An empty biomarker dict should produce no alerts."""
        alerts = engine.check({})
        assert len(alerts) == 0

    def test_unknown_biomarker_ignored(self, engine):
        """A biomarker not in the rules should produce no alerts."""
        alerts = engine.check({"xyz_unknown": 9999.0})
        assert len(alerts) == 0


# =====================================================================
# ALIAS RESOLUTION
# =====================================================================


class TestAliasResolution:
    """Tests for biomarker name alias resolution."""

    def test_glucose_alias(self, engine):
        """'glucose' resolves to 'Glucose Fasting'."""
        assert engine._resolve_biomarker("glucose") == "Glucose"

    def test_fasting_glucose_alias(self, engine):
        """'fasting_glucose' resolves to 'Glucose Fasting'."""
        assert engine._resolve_biomarker("fasting_glucose") == "Glucose"

    def test_case_insensitive(self, engine):
        """Alias resolution is case-insensitive."""
        assert engine._resolve_biomarker("GLUCOSE") == "Glucose"
        assert engine._resolve_biomarker("Potassium") == "Potassium"

    def test_canonical_name_resolves(self, engine):
        """Canonical name itself resolves correctly."""
        assert engine._resolve_biomarker("Glucose") == "Glucose"

    def test_unknown_alias_returns_none(self, engine):
        """An unrecognized name returns None."""
        assert engine._resolve_biomarker("xyz_unknown") is None

    def test_alerts_from_aliased_input(self, engine):
        """Aliased biomarker names correctly trigger threshold checks."""
        alerts = engine.check({"fasting_glucose": 300.0})
        assert len(alerts) == 1
        assert alerts[0].biomarker == "Glucose"

    def test_plt_alias(self, engine):
        """'plt' resolves to 'Platelet Count'."""
        assert engine._resolve_biomarker("plt") == "Platelet Count"


# =====================================================================
# SEVERITY SORTING
# =====================================================================


class TestSeveritySorting:
    """Tests for alert severity ordering."""

    def test_alerts_sorted_critical_first(self, engine):
        """Alerts are returned sorted: critical > urgent > warning."""
        # Trigger multiple alerts across severity levels:
        # potassium 7.2 -> critical (>6.0), glucose 300 -> urgent (>250), platelets 110 -> warning (>100)
        alerts = engine.check({
            "potassium": 7.2,
            "glucose": 300.0,
            "platelets": 110.0,
        })
        assert len(alerts) == 3
        severities = [a.severity for a in alerts]
        assert severities == ["critical", "urgent", "warning"]

    def test_multiple_same_severity(self):
        """Multiple alerts of the same severity are grouped together."""
        rules = [
            {
                "biomarker": "Potassium",
                "critical_high": 6.0,
                "severity": "critical",
                "escalation_target": "team_a",
                "clinical_action": "Action A",
                "cross_checks": [],
            },
            {
                "biomarker": "Sodium",
                "critical_high": 160.0,
                "severity": "critical",
                "escalation_target": "team_b",
                "clinical_action": "Action B",
                "cross_checks": [],
            },
        ]
        engine = CriticalValueEngine(rules=rules)
        alerts = engine.check({"potassium": 7.0, "sodium": 165.0})
        assert len(alerts) == 2
        assert all(a.severity == "critical" for a in alerts)


# =====================================================================
# GRACEFUL HANDLING — MISSING JSON FILE
# =====================================================================


class TestMissingJsonFile:
    """Tests for graceful behavior when the data file is absent."""

    def test_missing_file_returns_empty_rules(self, tmp_path):
        """Engine loads with zero rules when JSON file does not exist."""
        with patch("src.critical_values.REFERENCE_DIR", tmp_path):
            engine = CriticalValueEngine()
            assert engine.rule_count == 0

    def test_missing_file_check_returns_no_alerts(self, tmp_path):
        """check() returns an empty list when no rules are loaded."""
        with patch("src.critical_values.REFERENCE_DIR", tmp_path):
            engine = CriticalValueEngine()
            alerts = engine.check({"glucose": 999.0})
            assert alerts == []


# =====================================================================
# ALERT FORMATTING
# =====================================================================


class TestAlertFormatting:
    """Tests for CriticalValueAlert output methods."""

    @pytest.fixture
    def sample_alert(self):
        """Return a sample CriticalValueAlert."""
        return CriticalValueAlert(
            biomarker="Glucose",
            value=300.0,
            threshold=250.0,
            direction="high",
            severity="urgent",
            escalation_target="endocrinology_oncall",
            clinical_action="Check serum ketones and arterial blood gas.",
            cross_checks=["HbA1c", "Serum Ketones"],
            loinc_code="2345-7",
        )

    def test_to_alert_string_format(self, sample_alert):
        """to_alert_string() includes severity, biomarker, value, threshold, escalation."""
        result = sample_alert.to_alert_string()
        assert "CRITICAL VALUE [URGENT]" in result
        assert "Glucose" in result
        assert "300.0" in result
        assert "above threshold 250.0" in result
        assert "Escalate to endocrinology_oncall" in result
        assert "Action: Check serum ketones" in result

    def test_to_alert_string_low_direction(self):
        """to_alert_string() uses 'below' for low-direction alerts."""
        alert = CriticalValueAlert(
            biomarker="Glucose",
            value=35.0,
            threshold=50.0,
            direction="low",
            severity="critical",
            escalation_target="rapid_response_team",
            clinical_action="Administer IV dextrose.",
            cross_checks=[],
        )
        result = alert.to_alert_string()
        assert "below threshold 50.0" in result

    def test_to_dict_returns_expected_keys(self, sample_alert):
        """to_dict() returns all expected keys."""
        d = sample_alert.to_dict()
        expected_keys = {
            "biomarker", "value", "threshold", "direction",
            "severity", "escalation_target", "clinical_action",
            "cross_checks", "loinc_code",
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_values(self, sample_alert):
        """to_dict() returns correct values."""
        d = sample_alert.to_dict()
        assert d["biomarker"] == "Glucose"
        assert d["value"] == 300.0
        assert d["threshold"] == 250.0
        assert d["direction"] == "high"
        assert d["severity"] == "urgent"
        assert d["loinc_code"] == "2345-7"
        assert d["cross_checks"] == ["HbA1c", "Serum Ketones"]

    def test_to_dict_default_loinc_code(self):
        """to_dict() returns empty string for loinc_code when not provided."""
        alert = CriticalValueAlert(
            biomarker="Test",
            value=1.0,
            threshold=2.0,
            direction="low",
            severity="warning",
            escalation_target="test",
            clinical_action="test",
            cross_checks=[],
        )
        assert alert.to_dict()["loinc_code"] == ""
