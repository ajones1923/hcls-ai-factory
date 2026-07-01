"""Tests for DiscordanceDetector and DiscordanceFinding.

Validates cross-biomarker discordance detection, alias resolution,
priority sorting, alert formatting, and graceful degradation when
the reference data file is missing.

Author: Adam Jones
Date: March 2026
"""

import json
from unittest.mock import patch

import pytest

from src.discordance_detector import DiscordanceDetector, DiscordanceFinding

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def sample_rules():
    """Return a minimal set of discordance rules for testing."""
    return [
        {
            "id": "disc_ferritin_low_tsat_high",
            "name": "Contradictory Iron Studies",
            "biomarker_a": "Ferritin",
            "biomarker_b": "Transferrin Saturation (TSAT)",
            "condition": "Ferritin LOW and TSAT HIGH",
            "expected_relationship": "Ferritin and TSAT typically move in the same direction.",
            "differential_diagnosis": [
                "Specimen mislabeling or pre-analytical error",
                "Acute hepatocellular necrosis",
                "Recent parenteral iron infusion",
            ],
            "agent_handoff": ["lab_quality_agent", "iron_metabolism_agent"],
            "priority": "high",
            "text_chunk": "Ferritin and TSAT normally trend concordantly.",
        },
        {
            "id": "disc_ft3_high_tsh_normal",
            "name": "Thyroid Hormone-TSH Discordance",
            "biomarker_a": "Free T3",
            "biomarker_b": "TSH",
            "condition": "Free T3 HIGH and TSH NORMAL",
            "expected_relationship": "Elevated Free T3 should suppress TSH.",
            "differential_diagnosis": [
                "TSH-secreting pituitary adenoma",
                "Thyroid hormone resistance syndrome",
                "Assay interference",
            ],
            "agent_handoff": ["endocrine_agent", "pituitary_agent"],
            "priority": "high",
            "text_chunk": "Elevated Free T3 with a non-suppressed TSH.",
        },
        {
            "id": "disc_fibrinogen_high_crp_normal",
            "name": "Acute Phase Reactant Discordance",
            "biomarker_a": "Fibrinogen",
            "biomarker_b": "hs-CRP",
            "condition": "Fibrinogen HIGH and hs-CRP NORMAL",
            "expected_relationship": "Both are acute phase reactants.",
            "differential_diagnosis": [
                "Estrogen effect",
                "Nephrotic syndrome",
                "Smoking",
            ],
            "agent_handoff": ["coagulation_agent", "inflammation_agent"],
            "priority": "medium",
            "text_chunk": "Fibrinogen and CRP normally co-elevate.",
        },
        {
            "id": "disc_ldl_apob",
            "name": "LDL-ApoB Discordance",
            "biomarker_a": "LDL-C",
            "biomarker_b": "ApoB",
            "condition": "LDL-C DISCORDANT with ApoB",
            "expected_relationship": "LDL-C and ApoB should correlate.",
            "differential_diagnosis": [
                "Small dense LDL predominance",
                "Residual cholesterol-rich remnants",
            ],
            "agent_handoff": ["lipid_agent"],
            "priority": "low",
            "text_chunk": "LDL-C and ApoB discordance suggests particle mismatch.",
        },
    ]


@pytest.fixture
def detector(sample_rules):
    """Return a DiscordanceDetector initialized with sample rules."""
    return DiscordanceDetector(rules=sample_rules)


# =====================================================================
# INITIALIZATION
# =====================================================================


class TestInitialization:
    """Tests for DiscordanceDetector construction and rule loading."""

    def test_init_with_provided_rules(self, sample_rules):
        """Detector uses rules passed directly to the constructor."""
        detector = DiscordanceDetector(rules=sample_rules)
        assert detector.rule_count == len(sample_rules)

    def test_init_loads_rules_from_json(self, sample_rules, tmp_path):
        """Detector loads rules from the JSON file when none are provided."""
        json_file = tmp_path / "biomarker_discordance_rules.json"
        json_file.write_text(json.dumps(sample_rules))

        with patch("src.discordance_detector.REFERENCE_DIR", tmp_path):
            detector = DiscordanceDetector()
            assert detector.rule_count == len(sample_rules)


# =====================================================================
# CHECK — KNOWN DISCORDANCE PATTERNS
# =====================================================================


class TestCheckDiscordance:
    """Tests for detecting known discordance patterns."""

    def test_ferritin_low_tsat_high_detected(self, detector):
        """Low ferritin + high TSAT triggers the contradictory iron studies rule."""
        # Ferritin LOW threshold = 30, TSAT HIGH threshold = 45
        # Use 'transferrin_saturation' alias which resolves to "Transferrin Saturation (TSAT)"
        findings = detector.check({"ferritin": 15.0, "transferrin_saturation": 55.0})
        assert len(findings) >= 1
        iron_findings = [f for f in findings if f.rule_name == "Contradictory Iron Studies"]
        assert len(iron_findings) == 1
        finding = iron_findings[0]
        assert finding.biomarker_a == "Ferritin"
        assert finding.biomarker_b == "Transferrin Saturation (TSAT)"
        assert finding.value_a == 15.0
        assert finding.value_b == 55.0
        assert finding.priority == "high"

    def test_ft3_high_tsh_normal_detected(self, detector):
        """High Free T3 + normal TSH triggers thyroid discordance rule."""
        # Free T3 HIGH threshold = 4.4, TSH NORMAL = 0.45-4.5
        findings = detector.check({"free_t3": 5.5, "tsh": 2.0})
        thyroid_findings = [f for f in findings if f.rule_name == "Thyroid Hormone-TSH Discordance"]
        assert len(thyroid_findings) == 1

    def test_fibrinogen_high_crp_normal_detected(self, detector):
        """High fibrinogen + normal CRP triggers acute phase discordance."""
        # Fibrinogen HIGH threshold = 400, hs-CRP NORMAL = 0-3.0
        findings = detector.check({"fibrinogen": 500.0, "hs_crp": 1.5})
        phase_findings = [f for f in findings if f.rule_name == "Acute Phase Reactant Discordance"]
        assert len(phase_findings) == 1
        assert phase_findings[0].priority == "medium"

    def test_ldl_apob_discordance_detected(self, detector):
        """LDL-C normal + ApoB high triggers LDL-ApoB discordance."""
        # LDL-C NORMAL = 0-100, ApoB HIGH = >90
        findings = detector.check({"ldl_c": 85.0, "apob": 120.0})
        ldl_findings = [f for f in findings if f.rule_name == "LDL-ApoB Discordance"]
        assert len(ldl_findings) == 1


# =====================================================================
# CHECK — CONCORDANT / NORMAL VALUES
# =====================================================================


class TestCheckConcordant:
    """Tests for concordant or normal values producing no findings."""

    def test_concordant_iron_studies_no_findings(self, detector):
        """Both ferritin and TSAT normal produces no iron discordance."""
        findings = detector.check({"ferritin": 150.0, "tsat": 30.0})
        iron_findings = [f for f in findings if f.rule_name == "Contradictory Iron Studies"]
        assert len(iron_findings) == 0

    def test_concordant_thyroid_no_findings(self, detector):
        """Normal Free T3 + normal TSH produces no thyroid discordance."""
        findings = detector.check({"free_t3": 3.0, "tsh": 2.0})
        thyroid_findings = [f for f in findings if f.rule_name == "Thyroid Hormone-TSH Discordance"]
        assert len(thyroid_findings) == 0

    def test_missing_one_biomarker_no_findings(self, detector):
        """When only one biomarker of a pair is present, no findings."""
        findings = detector.check({"ferritin": 15.0})
        iron_findings = [f for f in findings if f.rule_name == "Contradictory Iron Studies"]
        assert len(iron_findings) == 0

    def test_empty_biomarkers_no_findings(self, detector):
        """Empty biomarker dict produces no findings."""
        findings = detector.check({})
        assert findings == []


# =====================================================================
# ALIAS RESOLUTION
# =====================================================================


class TestAliasResolution:
    """Tests for biomarker name alias resolution."""

    def test_resolve_ferritin_alias(self, detector):
        """'ferritin' resolves to 'Ferritin'."""
        assert detector._resolve("ferritin") == "Ferritin"

    def test_resolve_tsat_alias(self, detector):
        """'tsat' resolves to 'TSAT'."""
        assert detector._resolve("tsat") == "TSAT"

    def test_resolve_transferrin_saturation_alias(self, detector):
        """'transferrin_saturation' resolves to 'Transferrin Saturation (TSAT)'."""
        assert detector._resolve("transferrin_saturation") == "Transferrin Saturation (TSAT)"

    def test_resolve_hs_crp_alias(self, detector):
        """'hs_crp' resolves to 'hs-CRP'."""
        assert detector._resolve("hs_crp") == "hs-CRP"

    def test_resolve_case_insensitive_direct_match(self, detector):
        """Direct match is case-insensitive against THRESHOLDS keys."""
        assert detector._resolve("Ferritin") == "Ferritin"
        assert detector._resolve("FERRITIN") == "Ferritin"

    def test_resolve_unknown_returns_none(self, detector):
        """An unrecognized name returns None."""
        assert detector._resolve("xyz_unknown") is None

    def test_aliased_input_triggers_detection(self, detector):
        """Aliased inputs correctly trigger discordance detection."""
        # Use 'ft3' alias instead of 'Free T3'
        findings = detector.check({"ft3": 5.5, "tsh": 2.0})
        thyroid_findings = [f for f in findings if f.rule_name == "Thyroid Hormone-TSH Discordance"]
        assert len(thyroid_findings) == 1


# =====================================================================
# PRIORITY SORTING
# =====================================================================


class TestPrioritySorting:
    """Tests for findings priority ordering."""

    def test_findings_sorted_high_first(self, detector):
        """Findings are returned sorted: high > medium > low."""
        # Trigger multiple findings across priorities:
        # Iron (high), Fibrinogen-CRP (medium), LDL-ApoB (low)
        findings = detector.check({
            "ferritin": 15.0,
            "transferrin_saturation": 55.0,
            "fibrinogen": 500.0,
            "hs_crp": 1.5,
            "ldl_c": 85.0,
            "apob": 120.0,
        })
        assert len(findings) >= 3
        priorities = [f.priority for f in findings]
        # Verify sorted order: all 'high' before 'medium' before 'low'
        high_indices = [i for i, p in enumerate(priorities) if p == "high"]
        medium_indices = [i for i, p in enumerate(priorities) if p == "medium"]
        low_indices = [i for i, p in enumerate(priorities) if p == "low"]
        if high_indices and medium_indices:
            assert max(high_indices) < min(medium_indices)
        if medium_indices and low_indices:
            assert max(medium_indices) < min(low_indices)


# =====================================================================
# GRACEFUL HANDLING — MISSING JSON FILE
# =====================================================================


class TestMissingJsonFile:
    """Tests for graceful behavior when the data file is absent."""

    def test_missing_file_returns_empty_rules(self, tmp_path):
        """Detector loads with zero rules when JSON file does not exist."""
        with patch("src.discordance_detector.REFERENCE_DIR", tmp_path):
            detector = DiscordanceDetector()
            assert detector.rule_count == 0

    def test_missing_file_check_returns_no_findings(self, tmp_path):
        """check() returns an empty list when no rules are loaded."""
        with patch("src.discordance_detector.REFERENCE_DIR", tmp_path):
            detector = DiscordanceDetector()
            findings = detector.check({"ferritin": 5.0, "tsat": 90.0})
            assert findings == []


# =====================================================================
# FINDING FORMATTING
# =====================================================================


class TestFindingFormatting:
    """Tests for DiscordanceFinding output methods."""

    @pytest.fixture
    def sample_finding(self):
        """Return a sample DiscordanceFinding."""
        return DiscordanceFinding(
            rule_name="Contradictory Iron Studies",
            biomarker_a="Ferritin",
            biomarker_b="Transferrin Saturation (TSAT)",
            value_a=15.0,
            value_b=55.0,
            condition="Ferritin LOW and TSAT HIGH",
            differential_diagnosis=[
                "Specimen mislabeling",
                "Acute hepatocellular necrosis",
                "Recent parenteral iron infusion",
            ],
            agent_handoff=["lab_quality_agent", "iron_metabolism_agent"],
            priority="high",
            text_chunk="Ferritin and TSAT normally trend concordantly.",
        )

    def test_to_alert_string_format(self, sample_finding):
        """to_alert_string() includes priority, rule name, values, and differentials."""
        result = sample_finding.to_alert_string()
        assert "DISCORDANCE [HIGH]" in result
        assert "Contradictory Iron Studies" in result
        assert "Ferritin=15.0" in result
        assert "Transferrin Saturation (TSAT)=55.0" in result
        assert "Specimen mislabeling" in result

    def test_to_alert_string_limits_differentials_to_three(self):
        """to_alert_string() shows at most 3 differential diagnoses."""
        finding = DiscordanceFinding(
            rule_name="Test",
            biomarker_a="A",
            biomarker_b="B",
            value_a=1.0,
            value_b=2.0,
            condition="test",
            differential_diagnosis=["D1", "D2", "D3", "D4", "D5"],
            agent_handoff=[],
            priority="medium",
            text_chunk="",
        )
        result = finding.to_alert_string()
        assert "D1" in result
        assert "D2" in result
        assert "D3" in result
        # D4 and D5 should not appear (sliced to [:3])
        assert "D4" not in result
        assert "D5" not in result

    def test_to_dict_returns_expected_keys(self, sample_finding):
        """to_dict() returns all expected keys."""
        d = sample_finding.to_dict()
        expected_keys = {
            "rule_name", "biomarker_a", "biomarker_b",
            "value_a", "value_b", "condition",
            "differential_diagnosis", "agent_handoff", "priority",
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_values(self, sample_finding):
        """to_dict() returns correct values."""
        d = sample_finding.to_dict()
        assert d["rule_name"] == "Contradictory Iron Studies"
        assert d["biomarker_a"] == "Ferritin"
        assert d["value_a"] == 15.0
        assert d["value_b"] == 55.0
        assert d["priority"] == "high"
        assert len(d["differential_diagnosis"]) == 3
        assert len(d["agent_handoff"]) == 2

    def test_to_dict_excludes_text_chunk(self, sample_finding):
        """to_dict() does not include the text_chunk field."""
        d = sample_finding.to_dict()
        assert "text_chunk" not in d
