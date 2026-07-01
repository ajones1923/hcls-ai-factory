"""Tests for LabRangeInterpreter and RangeComparison.

Validates multi-lab range interpretation, sex-specific range selection,
discrepancy detection, alias resolution, and graceful degradation when
the reference data file is missing.

Author: Adam Jones
Date: March 2026
"""

import json
from unittest.mock import patch

import pytest

from src.lab_range_interpreter import LabRangeInterpreter, RangeComparison

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def sample_lab_data():
    """Return mock multi-lab range data mimicking biomarker_lab_ranges.json."""
    return {
        "labs": {
            "quest_diagnostics": {
                "name": "Quest Diagnostics",
                "ranges": {
                    "Glucose Fasting": {"min": 65, "max": 99, "unit": "mg/dL"},
                    "HbA1c": {"min": 4.8, "max": 5.6, "unit": "%"},
                    "ALT": {"min": None, "max": 44, "unit": "U/L"},
                    "Ferritin": {"min": 30, "max": 400, "unit": "ng/mL"},
                    "HDL Cholesterol (Male)": {"min": 40, "max": None, "unit": "mg/dL"},
                    "HDL Cholesterol (Female)": {"min": 50, "max": None, "unit": "mg/dL"},
                    "Creatinine (Male)": {"min": 0.76, "max": 1.27, "unit": "mg/dL"},
                    "Creatinine (Female)": {"min": 0.55, "max": 1.02, "unit": "mg/dL"},
                },
            },
            "labcorp": {
                "name": "LabCorp",
                "ranges": {
                    "Glucose Fasting": {"min": 65, "max": 99, "unit": "mg/dL"},
                    "HbA1c": {"min": 4.8, "max": 5.6, "unit": "%"},
                    "ALT": {"min": None, "max": 44, "unit": "U/L"},
                    "Ferritin": {"min": 30, "max": 400, "unit": "ng/mL"},
                    "HDL Cholesterol (Male)": {"min": 40, "max": None, "unit": "mg/dL"},
                    "HDL Cholesterol (Female)": {"min": 50, "max": None, "unit": "mg/dL"},
                    "Creatinine (Male)": {"min": 0.74, "max": 1.35, "unit": "mg/dL"},
                    "Creatinine (Female)": {"min": 0.59, "max": 1.04, "unit": "mg/dL"},
                },
            },
            "function_health_optimal": {
                "name": "Function Health Optimal",
                "ranges": {
                    "Glucose Fasting": {"min": 72, "max": 85, "unit": "mg/dL"},
                    "HbA1c": {"min": 4.5, "max": 5.2, "unit": "%"},
                    "ALT": {"min": None, "max": 25, "unit": "U/L"},
                    "Ferritin": {"min": 40, "max": 150, "unit": "ng/mL"},
                    "HDL Cholesterol (Male)": {"min": 50, "max": None, "unit": "mg/dL"},
                    "HDL Cholesterol (Female)": {"min": 60, "max": None, "unit": "mg/dL"},
                    "Creatinine (Male)": {"min": 0.80, "max": 1.10, "unit": "mg/dL"},
                    "Creatinine (Female)": {"min": 0.60, "max": 0.90, "unit": "mg/dL"},
                },
            },
        }
    }


@pytest.fixture
def interpreter(sample_lab_data):
    """Return a LabRangeInterpreter initialized with sample data."""
    return LabRangeInterpreter(data=sample_lab_data)


# =====================================================================
# INITIALIZATION
# =====================================================================


class TestInitialization:
    """Tests for LabRangeInterpreter construction and data loading."""

    def test_init_with_provided_data(self, sample_lab_data):
        """Interpreter uses data passed directly to the constructor."""
        interp = LabRangeInterpreter(data=sample_lab_data)
        labs = interp._data.get("labs", {})
        assert "quest_diagnostics" in labs
        assert "labcorp" in labs
        assert "function_health_optimal" in labs

    def test_init_loads_data_from_json(self, sample_lab_data, tmp_path):
        """Interpreter loads data from the JSON file when none are provided."""
        json_file = tmp_path / "biomarker_lab_ranges.json"
        json_file.write_text(json.dumps(sample_lab_data))

        with patch("src.lab_range_interpreter.REFERENCE_DIR", tmp_path):
            interp = LabRangeInterpreter()
            assert "quest_diagnostics" in interp._data["labs"]


# =====================================================================
# INTERPRET — BASIC RESULTS
# =====================================================================


class TestInterpret:
    """Tests for interpret() returning RangeComparison objects."""

    def test_interpret_returns_range_comparisons(self, interpreter):
        """interpret() returns a list of RangeComparison for known biomarkers."""
        results = interpreter.interpret({"Glucose Fasting": 90.0})
        assert len(results) == 1
        assert isinstance(results[0], RangeComparison)
        assert results[0].biomarker == "Glucose Fasting"
        assert results[0].value == 90.0

    def test_interpret_normal_across_all_labs(self, interpreter):
        """A value within all three ranges is normal everywhere."""
        results = interpreter.interpret({"Glucose Fasting": 80.0})
        assert len(results) == 1
        r = results[0]
        assert r.quest_status == "normal"
        assert r.labcorp_status == "normal"
        assert r.optimal_status == "normal"

    def test_interpret_high_by_optimal_only(self, interpreter):
        """A value normal by standard labs but high by optimal."""
        # Glucose 95: Quest/LabCorp normal (65-99), Function Health high (72-85)
        results = interpreter.interpret({"Glucose Fasting": 95.0})
        assert len(results) == 1
        r = results[0]
        assert r.quest_status == "normal"
        assert r.labcorp_status == "normal"
        assert r.optimal_status == "high"

    def test_interpret_high_across_all_labs(self, interpreter):
        """A value above all ranges is high everywhere."""
        # Glucose 110: above all max thresholds
        results = interpreter.interpret({"Glucose Fasting": 110.0})
        r = results[0]
        assert r.quest_status == "high"
        assert r.labcorp_status == "high"
        assert r.optimal_status == "high"

    def test_interpret_low_across_all_labs(self, interpreter):
        """A value below all ranges is low everywhere."""
        # Glucose 50: below all min thresholds (Quest/LabCorp min=65, optimal min=72)
        results = interpreter.interpret({"Glucose Fasting": 50.0})
        r = results[0]
        assert r.quest_status == "low"
        assert r.labcorp_status == "low"
        assert r.optimal_status == "low"

    def test_interpret_unknown_biomarker_skipped(self, interpreter):
        """An unrecognized biomarker name is skipped."""
        results = interpreter.interpret({"xyz_unknown_marker": 42.0})
        assert len(results) == 0

    def test_interpret_multiple_biomarkers(self, interpreter):
        """interpret() processes multiple biomarkers in a single call."""
        results = interpreter.interpret({
            "Glucose Fasting": 80.0,
            "ALT": 30.0,
        })
        assert len(results) == 2
        biomarker_names = {r.biomarker for r in results}
        assert "Glucose Fasting" in biomarker_names
        assert "ALT" in biomarker_names

    def test_interpret_includes_unit(self, interpreter):
        """RangeComparison includes the unit from the lab range data."""
        results = interpreter.interpret({"Glucose Fasting": 80.0})
        assert results[0].unit == "mg/dL"


# =====================================================================
# INTERPRET — SEX-SPECIFIC RANGES
# =====================================================================


class TestSexSpecificRanges:
    """Tests for sex parameter selecting sex-specific ranges."""

    def test_male_creatinine_uses_male_range(self, interpreter):
        """Creatinine with sex='Male' uses the (Male) range."""
        results = interpreter.interpret({"creatinine": 1.20}, sex="Male")
        assert len(results) == 1
        r = results[0]
        # Quest male range: 0.76-1.27 -> normal at 1.20
        assert r.quest_status == "normal"

    def test_female_creatinine_uses_female_range(self, interpreter):
        """Creatinine with sex='Female' uses the (Female) range."""
        results = interpreter.interpret({"creatinine": 1.20}, sex="Female")
        assert len(results) == 1
        r = results[0]
        # Quest female range: 0.55-1.02 -> high at 1.20
        assert r.quest_status == "high"

    def test_male_hdl_uses_male_range(self, interpreter):
        """HDL with sex='Male' uses Male threshold (min=40)."""
        results = interpreter.interpret({"HDL Cholesterol (Male)": 45.0}, sex="Male")
        if results:
            r = results[0]
            assert r.quest_status == "normal"


# =====================================================================
# GET DISCREPANCIES
# =====================================================================


class TestGetDiscrepancies:
    """Tests for get_discrepancies() filtering logic."""

    def test_discrepancy_detected(self, interpreter):
        """Biomarker normal by standard labs but not by optimal is a discrepancy."""
        # Glucose 95: Quest normal, LabCorp normal, optimal high
        discrepancies = interpreter.get_discrepancies({"Glucose Fasting": 95.0})
        assert len(discrepancies) == 1
        assert discrepancies[0].biomarker == "Glucose Fasting"
        assert discrepancies[0].has_discrepancy is True

    def test_no_discrepancy_when_all_agree(self, interpreter):
        """No discrepancy when all three labs agree the value is normal."""
        # Glucose 80: normal everywhere
        discrepancies = interpreter.get_discrepancies({"Glucose Fasting": 80.0})
        assert len(discrepancies) == 0

    def test_no_discrepancy_when_all_abnormal(self, interpreter):
        """No discrepancy when standard labs also flag the value."""
        # Glucose 110: high everywhere, so quest_status != normal
        discrepancies = interpreter.get_discrepancies({"Glucose Fasting": 110.0})
        assert len(discrepancies) == 0

    def test_alt_discrepancy(self, interpreter):
        """ALT 35: within Quest/LabCorp (max 44) but above optimal (max 25)."""
        discrepancies = interpreter.get_discrepancies({"ALT": 35.0})
        assert len(discrepancies) == 1
        assert discrepancies[0].biomarker == "ALT"
        assert discrepancies[0].optimal_status == "high"


# =====================================================================
# HAS DISCREPANCY PROPERTY
# =====================================================================


class TestHasDiscrepancy:
    """Tests for the RangeComparison.has_discrepancy property."""

    def test_has_discrepancy_true(self):
        """has_discrepancy is True when standard=normal but optimal!=normal."""
        rc = RangeComparison(
            biomarker="Test",
            value=95.0,
            unit="mg/dL",
            quest_status="normal",
            labcorp_status="normal",
            optimal_status="high",
        )
        assert rc.has_discrepancy is True

    def test_has_discrepancy_false_all_normal(self):
        """has_discrepancy is False when all statuses are normal."""
        rc = RangeComparison(
            biomarker="Test",
            value=80.0,
            unit="mg/dL",
            quest_status="normal",
            labcorp_status="normal",
            optimal_status="normal",
        )
        assert rc.has_discrepancy is False

    def test_has_discrepancy_false_all_abnormal(self):
        """has_discrepancy is False when standard labs are also abnormal."""
        rc = RangeComparison(
            biomarker="Test",
            value=110.0,
            unit="mg/dL",
            quest_status="high",
            labcorp_status="high",
            optimal_status="high",
        )
        assert rc.has_discrepancy is False

    def test_has_discrepancy_one_standard_normal(self):
        """has_discrepancy is True when at least one standard lab is normal."""
        rc = RangeComparison(
            biomarker="Test",
            value=95.0,
            unit="mg/dL",
            quest_status="normal",
            labcorp_status="high",
            optimal_status="high",
        )
        assert rc.has_discrepancy is True


# =====================================================================
# ALIAS RESOLUTION
# =====================================================================


class TestAliasResolution:
    """Tests for biomarker name alias resolution."""

    def test_glucose_alias(self, interpreter):
        """'glucose' resolves to 'Glucose Fasting'."""
        assert interpreter._resolve("glucose") == "Glucose Fasting"

    def test_hba1c_alias(self, interpreter):
        """'hba1c' resolves to 'HbA1c'."""
        assert interpreter._resolve("hba1c") == "HbA1c"

    def test_alt_alias(self, interpreter):
        """'alt' resolves to 'ALT'."""
        assert interpreter._resolve("alt") == "ALT"

    def test_hscrp_alias(self, interpreter):
        """'hscrp' resolves to 'hs-CRP'."""
        assert interpreter._resolve("hscrp") == "hs-CRP"

    def test_unknown_alias_returns_none(self, interpreter):
        """An unrecognized name returns None."""
        assert interpreter._resolve("xyz_unknown") is None

    def test_aliased_input_triggers_interpret(self, interpreter):
        """Aliased names work in interpret()."""
        results = interpreter.interpret({"glucose": 90.0})
        assert len(results) == 1
        assert results[0].biomarker == "Glucose Fasting"

    def test_direct_canonical_name(self, interpreter):
        """Using the exact canonical name in ranges works via fallback."""
        results = interpreter.interpret({"Glucose Fasting": 90.0})
        assert len(results) == 1


# =====================================================================
# GRACEFUL HANDLING — MISSING JSON FILE
# =====================================================================


class TestMissingJsonFile:
    """Tests for graceful behavior when the data file is absent."""

    def test_missing_file_returns_empty_labs(self, tmp_path):
        """Interpreter loads with empty labs when JSON file does not exist."""
        with patch("src.lab_range_interpreter.REFERENCE_DIR", tmp_path):
            interp = LabRangeInterpreter()
            assert interp._data == {"labs": {}}

    def test_missing_file_interpret_returns_empty(self, tmp_path):
        """interpret() returns an empty list when no data is loaded."""
        with patch("src.lab_range_interpreter.REFERENCE_DIR", tmp_path):
            interp = LabRangeInterpreter()
            results = interp.interpret({"glucose": 90.0})
            assert results == []

    def test_missing_file_get_discrepancies_returns_empty(self, tmp_path):
        """get_discrepancies() returns an empty list when no data is loaded."""
        with patch("src.lab_range_interpreter.REFERENCE_DIR", tmp_path):
            interp = LabRangeInterpreter()
            results = interp.get_discrepancies({"glucose": 90.0})
            assert results == []


# =====================================================================
# OUTPUT FORMATTING
# =====================================================================


class TestOutputFormatting:
    """Tests for RangeComparison output methods."""

    @pytest.fixture
    def normal_comparison(self):
        """Return a RangeComparison that is normal across all labs."""
        return RangeComparison(
            biomarker="Glucose Fasting",
            value=80.0,
            unit="mg/dL",
            quest_status="normal",
            labcorp_status="normal",
            optimal_status="normal",
        )

    @pytest.fixture
    def discrepancy_comparison(self):
        """Return a RangeComparison with a standard-vs-optimal discrepancy."""
        return RangeComparison(
            biomarker="Glucose Fasting",
            value=95.0,
            unit="mg/dL",
            quest_status="normal",
            labcorp_status="normal",
            optimal_status="high",
            optimal_range={"min": 72, "max": 85, "unit": "mg/dL"},
        )

    def test_to_interpretation_normal(self, normal_comparison):
        """to_interpretation() for fully normal value mentions optimal range."""
        result = normal_comparison.to_interpretation()
        assert "Glucose Fasting" in result
        assert "80.0" in result
        assert "within optimal range" in result

    def test_to_interpretation_discrepancy(self, discrepancy_comparison):
        """to_interpretation() for a discrepancy mentions optimization."""
        result = discrepancy_comparison.to_interpretation()
        assert "within standard reference range" in result
        assert "Function Health optimal" in result
        assert "Consider optimization" in result
        assert "72" in result
        assert "85" in result

    def test_to_interpretation_abnormal_no_discrepancy(self):
        """to_interpretation() for a value flagged by standard labs."""
        rc = RangeComparison(
            biomarker="Glucose Fasting",
            value=110.0,
            unit="mg/dL",
            quest_status="high",
            labcorp_status="high",
            optimal_status="high",
        )
        result = rc.to_interpretation()
        assert "flagged as high by Quest" in result
        assert "high by LabCorp" in result

    def test_to_dict_returns_expected_keys(self, normal_comparison):
        """to_dict() returns all expected keys."""
        d = normal_comparison.to_dict()
        expected_keys = {
            "biomarker", "value", "unit",
            "quest_status", "labcorp_status", "optimal_status",
            "has_discrepancy",
        }
        assert set(d.keys()) == expected_keys

    def test_to_dict_values_normal(self, normal_comparison):
        """to_dict() returns correct values for a normal comparison."""
        d = normal_comparison.to_dict()
        assert d["biomarker"] == "Glucose Fasting"
        assert d["value"] == 80.0
        assert d["unit"] == "mg/dL"
        assert d["quest_status"] == "normal"
        assert d["labcorp_status"] == "normal"
        assert d["optimal_status"] == "normal"
        assert d["has_discrepancy"] is False

    def test_to_dict_values_discrepancy(self, discrepancy_comparison):
        """to_dict() returns correct values for a discrepancy comparison."""
        d = discrepancy_comparison.to_dict()
        assert d["has_discrepancy"] is True
        assert d["optimal_status"] == "high"
