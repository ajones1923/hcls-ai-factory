"""Tests for the dose intelligence engine.

Validates dose recording, cumulative dose calculation, DRL comparison,
alert thresholds (including pediatric), optimization suggestions, and
population-level analytics.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.dose_intelligence import DoseIntelligenceEngine, DoseRecord


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def engine():
    return DoseIntelligenceEngine()


@pytest.fixture
def sample_record():
    return DoseRecord(
        patient_id="PAT-00001",
        study_date="2026-03-15",
        modality="CT",
        protocol="CT Head without contrast",
        body_region="head",
        effective_dose_msv=2.0,
        dlp_mgy_cm=900.0,
        ctdi_vol_mgy=45.0,
    )


@pytest.fixture
def high_dose_record():
    return DoseRecord(
        patient_id="PAT-00002",
        study_date="2026-03-15",
        modality="CT",
        protocol="CT Abdomen multiphase",
        body_region="abdomen",
        effective_dose_msv=28.0,
        dlp_mgy_cm=1800.0,
        ctdi_vol_mgy=30.0,
    )


@pytest.fixture
def pediatric_record():
    return DoseRecord(
        patient_id="PAT-PEDS-001",
        study_date="2026-03-15",
        modality="CT",
        protocol="CT Head without contrast",
        body_region="head",
        effective_dose_msv=1.0,
        pediatric=True,
    )


def _make_record(patient_id, protocol, dose, date="2026-03-15", modality="CT",
                 region="head", pediatric=False):
    return DoseRecord(
        patient_id=patient_id,
        study_date=date,
        modality=modality,
        protocol=protocol,
        body_region=region,
        effective_dose_msv=dose,
        pediatric=pediatric,
    )


# ═══════════════════════════════════════════════════════════════════════
# DOSE RECORDING TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestDoseRecording:
    """Test basic dose recording functionality."""

    def test_record_dose(self, engine, sample_record):
        engine.record_dose(sample_record)
        assert engine.registry_size == 1

    def test_record_multiple_doses(self, engine, sample_record, high_dose_record):
        engine.record_dose(sample_record)
        engine.record_dose(high_dose_record)
        assert engine.registry_size == 2

    def test_record_preserves_data(self, engine, sample_record):
        engine.record_dose(sample_record)
        cumulative = engine.get_cumulative_dose("PAT-00001")
        assert cumulative.total_effective_dose_msv == 2.0


# ═══════════════════════════════════════════════════════════════════════
# CUMULATIVE DOSE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestCumulativeDose:
    """Test cumulative dose calculation."""

    def test_cumulative_dose_single(self, engine, sample_record):
        engine.record_dose(sample_record)
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 2.0
        assert result.study_count == 1

    def test_cumulative_dose_multiple(self, engine):
        for i in range(5):
            engine.record_dose(_make_record("PAT-00001", "CT Head without contrast", 2.0,
                                            date=f"2026-0{i+1}-15"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 10.0
        assert result.study_count == 5

    def test_cumulative_by_modality(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0, modality="CT"))
        engine.record_dose(_make_record("PAT-00001", "Chest XR", 0.02, modality="XR", region="chest"))
        engine.record_dose(_make_record("PAT-00001", "CT Chest", 7.0, modality="CT", region="chest"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert "CT" in result.by_modality
        assert "XR" in result.by_modality
        assert result.by_modality["CT"] == pytest.approx(9.0, abs=0.1)
        assert result.by_modality["XR"] == pytest.approx(0.02, abs=0.01)

    def test_cumulative_by_body_region(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0, region="head"))
        engine.record_dose(_make_record("PAT-00001", "CT Chest", 7.0, region="chest"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert "head" in result.by_body_region
        assert "chest" in result.by_body_region

    def test_cumulative_date_range(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0, date="2026-01-15"))
        engine.record_dose(_make_record("PAT-00001", "CT Chest", 7.0, date="2026-03-20"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.date_range["first"] == "2026-01-15"
        assert result.date_range["last"] == "2026-03-20"

    def test_cumulative_empty_patient(self, engine):
        result = engine.get_cumulative_dose("NONEXISTENT")
        assert result.total_effective_dose_msv == 0.0
        assert result.study_count == 0
        assert result.alert_level == "normal"

    def test_cumulative_filters_by_patient(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0))
        engine.record_dose(_make_record("PAT-00002", "CT Chest", 7.0))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 2.0
        assert result.study_count == 1


# ═══════════════════════════════════════════════════════════════════════
# ALERT LEVEL TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestAlertLevels:
    """Test cumulative dose alert classification."""

    def test_alert_normal(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 5.0))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.alert_level == "normal"

    def test_alert_elevated(self, engine):
        # 20-50 mSv = elevated
        for i in range(5):
            engine.record_dose(_make_record("PAT-00001", "CT A/P", 7.0,
                                            date=f"2026-0{i+1}-15"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 35.0
        assert result.alert_level == "elevated"

    def test_alert_high(self, engine):
        # 50-100 mSv = high
        for i in range(8):
            engine.record_dose(_make_record("PAT-00001", "CT A/P", 10.0,
                                            date=f"2026-0{min(i+1,9)}-{15+i}"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 80.0
        assert result.alert_level == "high"

    def test_alert_critical(self, engine):
        # > 100 mSv = critical
        for i in range(8):
            engine.record_dose(_make_record("PAT-00001", "CT multiphase", 15.0,
                                            date=f"2026-0{min(i+1,9)}-{15+i}"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.total_effective_dose_msv == 120.0
        assert result.alert_level == "critical"
        assert result.alert_message is not None
        assert "CRITICAL" in result.alert_message

    def test_pediatric_thresholds(self, engine):
        # Pediatric thresholds are halved: 10 mSv = elevated (adult: 20 mSv)
        for i in range(3):
            engine.record_dose(_make_record("PAT-PEDS", "CT Head", 5.0,
                                            date=f"2026-0{i+1}-15", pediatric=True))
        result = engine.get_cumulative_dose("PAT-PEDS")
        # 15 mSv: adult=normal, pediatric=elevated (threshold 10 mSv)
        assert result.total_effective_dose_msv == 15.0
        assert result.alert_level == "elevated"
        assert "Pediatric" in result.alert_message

    def test_alert_message_present_when_elevated(self, engine):
        for i in range(4):
            engine.record_dose(_make_record("PAT-00001", "CT A/P", 8.0,
                                            date=f"2026-0{i+1}-15"))
        result = engine.get_cumulative_dose("PAT-00001")
        assert result.alert_message is not None
        assert len(result.alert_message) > 0


# ═══════════════════════════════════════════════════════════════════════
# DRL COMPARISON TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestDRLComparison:
    """Test dose vs DRL comparison."""

    def test_drl_comparison_below_achievable(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 1.5)
        result = engine.compare_to_drl(record)
        assert result.status == "below_achievable"
        assert result.ratio < 1.0

    def test_drl_comparison_below_drl(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 2.2)
        result = engine.compare_to_drl(record)
        assert result.status == "below_drl"

    def test_drl_comparison_above_drl(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 3.0)
        result = engine.compare_to_drl(record)
        assert result.status == "above_drl"
        assert result.ratio > 1.0

    def test_drl_comparison_significantly_above(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 5.0)
        result = engine.compare_to_drl(record)
        assert result.status == "significantly_above"
        assert result.ratio > 1.5

    def test_drl_no_match(self, engine):
        record = _make_record("PAT-00001", "MRI Brain", 0.0, modality="MRI")
        result = engine.compare_to_drl(record)
        assert result.status == "no_drl_available"

    def test_drl_values_populated(self, engine):
        record = _make_record("PAT-00001", "CTA Coronary", 4.0)
        result = engine.compare_to_drl(record)
        assert result.drl_msv == 7.0
        assert result.achievable_dose_msv == 3.5


# ═══════════════════════════════════════════════════════════════════════
# OPTIMIZATION SUGGESTIONS TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestOptimizationSuggestions:
    """Test dose reduction recommendations."""

    def test_optimization_suggestions_above_drl(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 4.0)
        result = engine.compare_to_drl(record)
        assert len(result.optimization_suggestions) > 0

    def test_suggestions_significantly_above(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 5.0)
        result = engine.compare_to_drl(record)
        suggestions = " ".join(result.optimization_suggestions).lower()
        assert "kv" in suggestions or "iterative" in suggestions or "reduce" in suggestions

    def test_suggestions_well_optimized(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 1.5)
        result = engine.compare_to_drl(record)
        suggestions = " ".join(result.optimization_suggestions).lower()
        assert "optimized" in suggestions or "achievable" in suggestions

    def test_get_optimization_suggestions_method(self, engine):
        record = _make_record("PAT-00001", "CT Head without contrast", 4.0)
        suggestions = engine.get_optimization_suggestions(record)
        assert isinstance(suggestions, list)
        assert len(suggestions) > 0

    def test_pediatric_suggestions(self, engine):
        record = _make_record("PAT-PEDS", "CT Head without contrast", 2.0, pediatric=True)
        result = engine.compare_to_drl(record)
        pediatric_suggestions = [s for s in result.optimization_suggestions if "PEDIATRIC" in s]
        assert len(pediatric_suggestions) > 0


# ═══════════════════════════════════════════════════════════════════════
# POPULATION SUMMARY TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestPopulationSummary:
    """Test population-level dose analytics."""

    def test_population_summary_empty(self, engine):
        summary = engine.get_population_dose_summary()
        assert summary["total_records"] == 0
        assert summary["unique_patients"] == 0

    def test_population_summary_with_data(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0))
        engine.record_dose(_make_record("PAT-00002", "CT Chest", 7.0, region="chest"))
        engine.record_dose(_make_record("PAT-00001", "CT A/P", 10.0, region="abdomen"))
        summary = engine.get_population_dose_summary()
        assert summary["total_records"] == 3
        assert summary["unique_patients"] == 2
        assert summary["mean_dose_msv"] > 0
        assert "CT" in summary["by_modality"]

    def test_population_summary_by_protocol(self, engine):
        for _ in range(3):
            engine.record_dose(_make_record("PAT-00001", "CT Head", 2.0))
        for _ in range(2):
            engine.record_dose(_make_record("PAT-00002", "CT Chest", 7.0, region="chest"))
        summary = engine.get_population_dose_summary()
        assert "CT Head" in summary["by_protocol"]
        assert summary["by_protocol"]["CT Head"]["count"] == 3

    def test_population_alert_distribution(self, engine):
        engine.record_dose(_make_record("PAT-00001", "CT Head", 5.0))
        summary = engine.get_population_dose_summary()
        assert "normal" in summary["alert_distribution"]


# ═══════════════════════════════════════════════════════════════════════
# DEMO DATA GENERATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestDemoData:
    """Test synthetic data generation."""

    def test_generate_demo_data(self, engine):
        engine.generate_demo_data(n_records=50)
        assert engine.registry_size == 50

    def test_generate_demo_data_default(self, engine):
        engine.generate_demo_data()
        assert engine.registry_size == 200

    def test_demo_data_has_multiple_patients(self, engine):
        engine.generate_demo_data(n_records=100)
        summary = engine.get_population_dose_summary()
        assert summary["unique_patients"] >= 10

    def test_demo_data_has_multiple_protocols(self, engine):
        engine.generate_demo_data(n_records=100)
        summary = engine.get_population_dose_summary()
        assert len(summary["by_protocol"]) >= 3

    def test_demo_data_has_multiple_modalities(self, engine):
        engine.generate_demo_data(n_records=100)
        summary = engine.get_population_dose_summary()
        assert len(summary["by_modality"]) >= 2

    def test_demo_data_reproducible(self, engine):
        engine.generate_demo_data(n_records=50)
        summary1 = engine.get_population_dose_summary()

        engine2 = DoseIntelligenceEngine()
        engine2.generate_demo_data(n_records=50)
        summary2 = engine2.get_population_dose_summary()

        assert summary1["mean_dose_msv"] == summary2["mean_dose_msv"]
