"""Tests for longitudinal tracking, wearable models, and translation."""

import pytest

from src.models import BiomarkerPanel, PatientHistory, WearableData
from src.translation import get_supported_languages, translate_report_headers, translate_term


class TestPatientHistory:
    """Tests for longitudinal patient tracking."""

    def test_empty_history(self):
        history = PatientHistory(patient_id="TEST-001")
        assert history.panel_count == 0
        assert history.date_range == (None, None)
        assert history.age_acceleration_trend() is None

    def test_single_panel(self):
        history = PatientHistory(
            patient_id="TEST-001",
            panels=[
                BiomarkerPanel(
                    date="2025-06-15",
                    biomarkers={"albumin": 4.0, "creatinine": 0.9},
                    biological_age=45.2,
                ),
            ],
        )
        assert history.panel_count == 1
        assert history.date_range == ("2025-06-15", "2025-06-15")
        assert history.age_acceleration_trend() is None  # Need >= 2

    def test_multiple_panels_worsening(self):
        history = PatientHistory(
            patient_id="TEST-001",
            panels=[
                BiomarkerPanel(date="2025-01-15", biomarkers={"albumin": 4.2}, biological_age=44.0),
                BiomarkerPanel(date="2025-06-15", biomarkers={"albumin": 3.8}, biological_age=46.0),
                BiomarkerPanel(date="2025-12-15", biomarkers={"albumin": 3.5}, biological_age=48.5),
            ],
        )
        assert history.panel_count == 3
        assert history.date_range == ("2025-01-15", "2025-12-15")
        assert history.age_acceleration_trend() == "worsening"

    def test_multiple_panels_improving(self):
        history = PatientHistory(
            patient_id="TEST-001",
            panels=[
                BiomarkerPanel(date="2025-01-15", biomarkers={"albumin": 3.5}, biological_age=50.0),
                BiomarkerPanel(date="2025-06-15", biomarkers={"albumin": 4.0}, biological_age=48.0),
                BiomarkerPanel(date="2025-12-15", biomarkers={"albumin": 4.3}, biological_age=46.5),
            ],
        )
        assert history.age_acceleration_trend() == "improving"

    def test_stable_trend(self):
        history = PatientHistory(
            patient_id="TEST-001",
            panels=[
                BiomarkerPanel(date="2025-01-15", biomarkers={"albumin": 4.0}, biological_age=45.0),
                BiomarkerPanel(date="2025-06-15", biomarkers={"albumin": 4.0}, biological_age=45.2),
            ],
        )
        assert history.age_acceleration_trend() == "stable"

    def test_trajectory_extraction(self):
        history = PatientHistory(
            patient_id="TEST-001",
            panels=[
                BiomarkerPanel(date="2025-01-15", biomarkers={}, biological_age=44.0),
                BiomarkerPanel(date="2025-06-15", biomarkers={}, biological_age=None),
                BiomarkerPanel(date="2025-12-15", biomarkers={}, biological_age=46.0),
            ],
        )
        trajectory = history.biological_age_trajectory()
        assert len(trajectory) == 2  # Skips None
        assert trajectory[0]["date"] == "2025-01-15"
        assert trajectory[1]["biological_age"] == 46.0


class TestWearableData:
    """Tests for wearable device data models."""

    def test_valid_wearable(self):
        data = WearableData(
            device_type="Apple Watch Series 9",
            measurement_date="2025-06-15",
            resting_heart_rate=62,
            heart_rate_variability=45.0,
            spo2_average=97.5,
            sleep_duration_hours=7.5,
            deep_sleep_pct=22.0,
            steps=8500,
            vo2_max=42.0,
        )
        assert data.resting_heart_rate == 62
        assert data.vo2_max == 42.0

    def test_partial_wearable(self):
        """Wearable data with only some fields populated."""
        data = WearableData(
            measurement_date="2025-06-15",
            resting_heart_rate=70,
            steps=5000,
        )
        assert data.device_type is None
        assert data.sleep_duration_hours is None

    def test_wearable_validation(self):
        """Test field range validation."""
        with pytest.raises(Exception):
            WearableData(measurement_date="2025-06-15", resting_heart_rate=500)
        with pytest.raises(Exception):
            WearableData(measurement_date="2025-06-15", spo2_average=110)


class TestTranslation:
    """Tests for multi-language translation."""

    def test_supported_languages(self):
        langs = get_supported_languages()
        assert "en" in langs
        assert "es" in langs
        assert "zh" in langs
        assert len(langs) >= 7

    def test_translate_term_spanish(self):
        result = translate_term("Biological Age", "es")
        assert result == "Edad Biológica"

    def test_translate_term_chinese(self):
        result = translate_term("Mortality Risk", "zh")
        assert result == "死亡风险"

    def test_translate_term_english_passthrough(self):
        result = translate_term("Biological Age", "en")
        assert result == "Biological Age"

    def test_translate_term_unknown_term(self):
        result = translate_term("Some Unknown Term", "es")
        assert result == "Some Unknown Term"

    def test_translate_term_unknown_language(self):
        result = translate_term("Biological Age", "xx")
        assert result == "Biological Age"

    def test_translate_report_headers(self):
        report = "# Biological Age Assessment\n\nRisk: HIGH\n\nMortality Risk: LOW"
        translated = translate_report_headers(report, "es")
        assert "Evaluación de Edad Biológica" in translated
        assert "ALTO" in translated
        assert "BAJO" in translated

    def test_translate_report_english_passthrough(self):
        report = "# Biological Age Assessment"
        assert translate_report_headers(report, "en") == report

    def test_all_risk_levels_translated(self):
        for risk in ["LOW", "MODERATE", "HIGH", "CRITICAL"]:
            for lang in ["es", "zh", "hi", "fr", "ar", "pt"]:
                translated = translate_term(risk, lang)
                assert translated != risk, f"{risk} not translated to {lang}"
