"""Tests for ImagingAnalyticsEngine.

Validates population analytics, cohort queries, temporal trends,
finding correlations, severity cross-tabulation, demo data generation,
and RAPIDS/pandas fallback behavior.

All tests work WITHOUT RAPIDS — they use the pandas fallback path.

Author: Adam Jones
Date: April 2026
"""

import random
import sys
from pathlib import Path

import pytest

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analytics_engine import (
    CohortCriteria,
    CohortResult,
    ImagingAnalyticsEngine,
    PopulationSummary,
    TrendResult,
    _BODY_REGIONS,
    _FINDING_TYPES,
    _MODALITIES,
    _SEVERITIES,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


def _make_study(
    study_id: str = "STUDY-000001",
    modality: str = "ct",
    body_region: str = "chest",
    severity: str = "routine",
    finding_type: str = "nodule",
    study_date: str = "2025-06-15",
    processing_time_ms: float = 250.0,
    patient_age: int = 55,
    patient_id: str = "PAT-00001",
    finding_size_mm: float = 12.0,
    **kwargs,
):
    """Create a single study dict with sensible defaults."""
    study = {
        "study_id": study_id,
        "patient_id": patient_id,
        "modality": modality,
        "body_region": body_region,
        "severity": severity,
        "finding_type": finding_type,
        "study_date": study_date,
        "processing_time_ms": processing_time_ms,
        "patient_age": patient_age,
        "finding_size_mm": finding_size_mm,
    }
    study.update(kwargs)
    return study


def _make_batch(n: int, seed: int = 42):
    """Generate n studies with deterministic random data."""
    rng = random.Random(seed)
    studies = []
    for i in range(n):
        studies.append(
            _make_study(
                study_id=f"STUDY-{i:06d}",
                patient_id=f"PAT-{rng.randint(1, n // 3 + 1):05d}",
                modality=rng.choice(_MODALITIES),
                body_region=rng.choice(_BODY_REGIONS),
                severity=rng.choice(_SEVERITIES),
                finding_type=rng.choice(_FINDING_TYPES),
                study_date=f"2025-{rng.randint(1, 12):02d}-{rng.randint(1, 28):02d}",
                processing_time_ms=round(rng.uniform(50, 5000), 1),
                patient_age=rng.randint(18, 95),
                finding_size_mm=round(rng.uniform(1.0, 80.0), 1),
            )
        )
    return studies


@pytest.fixture
def engine():
    """Fresh analytics engine for each test."""
    return ImagingAnalyticsEngine()


@pytest.fixture
def populated_engine():
    """Engine with 100 deterministic studies."""
    eng = ImagingAnalyticsEngine()
    eng.register_studies(_make_batch(100))
    return eng


# ═══════════════════════════════════════════════════════════════════════
# REGISTRATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestRegisterStudy:
    """Test single study registration."""

    def test_register_study(self, engine):
        assert engine.study_count == 0
        engine.register_study(_make_study())
        assert engine.study_count == 1

    def test_register_preserves_data(self, engine):
        study = _make_study(modality="mri", severity="critical")
        engine.register_study(study)
        assert engine._study_registry[0]["modality"] == "mri"
        assert engine._study_registry[0]["severity"] == "critical"


class TestRegisterBatch:
    """Test batch study registration."""

    def test_register_batch(self, engine):
        studies = _make_batch(100)
        engine.register_studies(studies)
        assert engine.study_count == 100

    def test_register_batch_additive(self, engine):
        engine.register_studies(_make_batch(50, seed=1))
        engine.register_studies(_make_batch(50, seed=2))
        assert engine.study_count == 100


# ═══════════════════════════════════════════════════════════════════════
# POPULATION SUMMARY TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestPopulationSummary:
    """Test population_summary() aggregations."""

    def test_population_summary(self, populated_engine):
        summary = populated_engine.population_summary()
        assert isinstance(summary, PopulationSummary)
        assert summary.total_studies == 100
        assert sum(summary.modality_distribution.values()) == 100
        assert sum(summary.body_region_distribution.values()) == 100
        assert sum(summary.severity_distribution.values()) == 100
        assert summary.mean_processing_time_ms > 0

    def test_population_summary_distributions_valid(self, populated_engine):
        summary = populated_engine.population_summary()
        # All modalities in distribution should be valid
        for mod in summary.modality_distribution:
            assert mod in _MODALITIES
        for region in summary.body_region_distribution:
            assert region in _BODY_REGIONS
        for sev in summary.severity_distribution:
            assert sev in _SEVERITIES

    def test_population_summary_finding_prevalence(self, populated_engine):
        summary = populated_engine.population_summary()
        # Prevalence should sum to ~100%
        total_prev = sum(summary.finding_prevalence.values())
        assert 99.0 <= total_prev <= 101.0

    def test_population_summary_date_range(self, populated_engine):
        summary = populated_engine.population_summary()
        assert summary.date_range is not None
        assert "start" in summary.date_range
        assert "end" in summary.date_range
        assert summary.date_range["start"] <= summary.date_range["end"]


class TestPopulationSummaryEmpty:
    """Test population_summary() with zero studies."""

    def test_population_summary_empty(self, engine):
        summary = engine.population_summary()
        assert summary.total_studies == 0
        assert summary.modality_distribution == {}
        assert summary.body_region_distribution == {}
        assert summary.severity_distribution == {}
        assert summary.finding_prevalence == {}
        assert summary.mean_processing_time_ms == 0.0
        assert summary.date_range is None
        assert summary.studies_with_critical_findings == 0
        assert summary.critical_finding_rate == 0.0


# ═══════════════════════════════════════════════════════════════════════
# COHORT QUERY TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestCohortQueryByModality:
    """Test cohort filtering by modality."""

    def test_cohort_query_by_modality(self, populated_engine):
        criteria = CohortCriteria(modality="ct")
        result = populated_engine.cohort_query(criteria)
        assert isinstance(result, CohortResult)
        assert result.total_studies == 100
        # All matched studies should be CT
        for study in result.studies:
            assert study["modality"] == "ct"
        assert result.matching_studies == len(result.studies)


class TestCohortQueryBySeverity:
    """Test cohort filtering by severity."""

    def test_cohort_query_by_severity(self, populated_engine):
        criteria = CohortCriteria(severity="critical")
        result = populated_engine.cohort_query(criteria)
        for study in result.studies:
            assert study["severity"] == "critical"
        assert result.matching_studies == len(result.studies)


class TestCohortQueryCombined:
    """Test cohort filtering with multiple criteria."""

    def test_cohort_query_combined(self):
        engine = ImagingAnalyticsEngine()
        # Register known studies
        engine.register_study(
            _make_study(modality="ct", severity="critical", study_date="2025-03-15")
        )
        engine.register_study(
            _make_study(
                study_id="STUDY-000002",
                modality="ct",
                severity="routine",
                study_date="2025-03-20",
            )
        )
        engine.register_study(
            _make_study(
                study_id="STUDY-000003",
                modality="mri",
                severity="critical",
                study_date="2025-03-10",
            )
        )

        criteria = CohortCriteria(
            modality="ct",
            severity="critical",
            start_date="2025-03-01",
            end_date="2025-03-31",
        )
        result = engine.cohort_query(criteria)
        assert result.matching_studies == 1
        assert result.studies[0]["modality"] == "ct"
        assert result.studies[0]["severity"] == "critical"


class TestCohortQueryNoMatches:
    """Test cohort query that matches nothing."""

    def test_cohort_query_no_matches(self, populated_engine):
        criteria = CohortCriteria(modality="fluoroscopy")
        result = populated_engine.cohort_query(criteria)
        # fluoroscopy is not in _MODALITIES used by _make_batch
        assert result.matching_studies == 0
        assert result.studies == []
        assert result.match_rate == 0.0


# ═══════════════════════════════════════════════════════════════════════
# TEMPORAL TREND TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestTemporalTrendsMonthly:
    """Test monthly temporal trend aggregation."""

    def test_temporal_trends_monthly(self, populated_engine):
        result = populated_engine.temporal_trends("total_volume", "monthly")
        assert isinstance(result, TrendResult)
        assert result.granularity == "monthly"
        assert result.metric == "total_volume"
        # Should have multiple months across the year
        assert len(result.data_points) > 1
        # Each point should have valid period format YYYY-MM
        for dp in result.data_points:
            assert len(dp.period.split("-")) == 2
            assert dp.count > 0


class TestTemporalTrendsDirection:
    """Test trend direction detection."""

    def test_increasing_trend(self):
        engine = ImagingAnalyticsEngine()
        # Create increasing critical findings over time
        for month in range(1, 7):
            for i in range(month * 5):  # 5, 10, 15, 20, 25, 30
                engine.register_study(
                    _make_study(
                        study_id=f"INC-{month}-{i}",
                        severity="critical",
                        study_date=f"2025-{month:02d}-15",
                    )
                )
        result = engine.temporal_trends("critical_findings", "monthly")
        assert result.trend_direction == "increasing"

    def test_decreasing_trend(self):
        engine = ImagingAnalyticsEngine()
        # Create decreasing critical findings over time
        for month in range(1, 7):
            for i in range((7 - month) * 5):  # 30, 25, 20, 15, 10, 5
                engine.register_study(
                    _make_study(
                        study_id=f"DEC-{month}-{i}",
                        severity="critical",
                        study_date=f"2025-{month:02d}-15",
                    )
                )
        result = engine.temporal_trends("critical_findings", "monthly")
        assert result.trend_direction == "decreasing"

    def test_stable_trend(self):
        engine = ImagingAnalyticsEngine()
        # Create stable findings: same count each month
        for month in range(1, 7):
            for i in range(10):
                engine.register_study(
                    _make_study(
                        study_id=f"STABLE-{month}-{i}",
                        severity="critical",
                        study_date=f"2025-{month:02d}-15",
                    )
                )
        result = engine.temporal_trends("critical_findings", "monthly")
        assert result.trend_direction == "stable"

    def test_empty_trend(self, engine):
        result = engine.temporal_trends("critical_findings")
        assert result.trend_direction == "stable"
        assert result.data_points == []


# ═══════════════════════════════════════════════════════════════════════
# FINDING CORRELATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestFindingCorrelation:
    """Test co-occurrence rate calculation."""

    def test_finding_correlation(self):
        engine = ImagingAnalyticsEngine()
        # Patient A has both nodule and effusion
        engine.register_study(
            _make_study(patient_id="PAT-A", finding_type="nodule")
        )
        engine.register_study(
            _make_study(
                study_id="S2", patient_id="PAT-A", finding_type="effusion"
            )
        )
        # Patient B has only nodule
        engine.register_study(
            _make_study(study_id="S3", patient_id="PAT-B", finding_type="nodule")
        )
        # Patient C has only effusion
        engine.register_study(
            _make_study(study_id="S4", patient_id="PAT-C", finding_type="effusion")
        )

        result = engine.finding_correlation("nodule", "effusion")
        assert result["finding_a"] == "nodule"
        assert result["finding_b"] == "effusion"
        assert result["finding_a_count"] == 2  # PAT-A and PAT-B
        assert result["finding_b_count"] == 2  # PAT-A and PAT-C
        assert result["co_occurrence_count"] == 1  # PAT-A only
        assert result["co_occurrence_rate"] == 0.5  # 1/2

    def test_finding_correlation_empty(self, engine):
        result = engine.finding_correlation("nodule", "effusion")
        assert result["co_occurrence_count"] == 0
        assert result["co_occurrence_rate"] == 0.0
        assert result["total_studies"] == 0


# ═══════════════════════════════════════════════════════════════════════
# SEVERITY BY MODALITY TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestSeverityByModality:
    """Test severity cross-tabulation by modality."""

    def test_severity_by_modality(self, populated_engine):
        result = populated_engine.severity_by_modality()
        assert isinstance(result, dict)
        # Should have entries for modalities present in the data
        assert len(result) > 0
        for mod, sev_counts in result.items():
            assert mod in _MODALITIES
            assert isinstance(sev_counts, dict)
            for sev in sev_counts:
                assert sev in _SEVERITIES

    def test_severity_by_modality_empty(self, engine):
        result = engine.severity_by_modality()
        assert result == {}

    def test_severity_by_modality_totals(self, populated_engine):
        result = populated_engine.severity_by_modality()
        # Sum of all counts should equal total studies
        total = sum(
            count for sev_counts in result.values() for count in sev_counts.values()
        )
        assert total == 100


# ═══════════════════════════════════════════════════════════════════════
# DEMO DATA GENERATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestGenerateDemoData:
    """Test synthetic data generation."""

    def test_generate_demo_data(self, engine):
        count = engine.generate_demo_data(n_studies=200)
        assert count == 200
        assert engine.study_count == 200

    def test_generate_demo_data_valid_fields(self, engine):
        engine.generate_demo_data(n_studies=100)
        for study in engine._study_registry:
            assert study["modality"] in _MODALITIES
            assert study["body_region"] in _BODY_REGIONS
            assert study["severity"] in _SEVERITIES
            assert study["finding_type"] in _FINDING_TYPES
            assert "study_date" in study
            assert "processing_time_ms" in study
            assert "patient_age" in study
            assert 18 <= study["patient_age"] <= 95

    def test_generate_demo_data_deterministic(self, engine):
        engine.generate_demo_data(n_studies=50, seed=123)
        first_run = [s["modality"] for s in engine._study_registry]

        engine2 = ImagingAnalyticsEngine()
        engine2.generate_demo_data(n_studies=50, seed=123)
        second_run = [s["modality"] for s in engine2._study_registry]

        assert first_run == second_run

    def test_generate_demo_data_additive(self, engine):
        engine.generate_demo_data(n_studies=50)
        engine.generate_demo_data(n_studies=50)
        assert engine.study_count == 100


# ═══════════════════════════════════════════════════════════════════════
# RAPIDS FALLBACK TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestRapidsFallback:
    """Test that engine works when cudf is not installed."""

    def test_rapids_fallback(self, engine):
        # In test environment, cudf is almost certainly not available
        assert engine.backend_name == "pandas"
        assert engine._rapids_available is False
        # All operations should still work
        engine.generate_demo_data(n_studies=50)
        summary = engine.population_summary()
        assert summary.total_studies == 50


# ═══════════════════════════════════════════════════════════════════════
# CRITICAL FINDING RATE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestCriticalFindingRate:
    """Test critical finding rate computation."""

    def test_critical_finding_rate(self):
        engine = ImagingAnalyticsEngine()
        # 3 critical out of 10 total
        for i in range(7):
            engine.register_study(
                _make_study(study_id=f"NORM-{i}", severity="routine")
            )
        for i in range(3):
            engine.register_study(
                _make_study(study_id=f"CRIT-{i}", severity="critical")
            )
        summary = engine.population_summary()
        assert summary.studies_with_critical_findings == 3
        assert summary.critical_finding_rate == 0.3

    def test_critical_finding_rate_zero(self):
        engine = ImagingAnalyticsEngine()
        for i in range(5):
            engine.register_study(
                _make_study(study_id=f"NORM-{i}", severity="normal")
            )
        summary = engine.population_summary()
        assert summary.studies_with_critical_findings == 0
        assert summary.critical_finding_rate == 0.0

    def test_critical_finding_rate_all_critical(self):
        engine = ImagingAnalyticsEngine()
        for i in range(5):
            engine.register_study(
                _make_study(study_id=f"CRIT-{i}", severity="critical")
            )
        summary = engine.population_summary()
        assert summary.studies_with_critical_findings == 5
        assert summary.critical_finding_rate == 1.0


# ═══════════════════════════════════════════════════════════════════════
# AI CONCORDANCE TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestAIConcordance:
    """Test AI vs radiologist agreement computation."""

    def test_concordance_agreement(self):
        engine = ImagingAnalyticsEngine()
        engine.register_study(
            _make_study(
                ai_severity="critical",
                radiologist_severity="critical",
            )
        )
        engine.register_study(
            _make_study(
                study_id="S2",
                ai_severity="urgent",
                radiologist_severity="routine",
            )
        )
        result = engine.ai_concordance()
        assert result["total_paired"] == 2
        assert result["agreement_count"] == 1
        assert result["agreement_rate"] == 0.5
        assert len(result["disagreements"]) == 1

    def test_concordance_empty(self, engine):
        result = engine.ai_concordance()
        assert result["total_paired"] == 0
        assert result["agreement_rate"] == 0.0
