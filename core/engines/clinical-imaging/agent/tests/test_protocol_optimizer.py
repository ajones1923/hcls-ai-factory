"""Tests for the evidence-based protocol optimizer.

Validates ACR Appropriateness Criteria matching, patient-specific
adjustments (pregnancy, renal impairment, contrast allergy, pediatric),
dose estimation, and alternative protocol ranking.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.protocol_optimizer import PatientFactors, ProtocolOptimizer, ProtocolRecommendation


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def optimizer():
    return ProtocolOptimizer()


@pytest.fixture
def default_patient():
    return PatientFactors(age=50, weight_kg=70.0, sex="M")


@pytest.fixture
def pregnant_patient():
    return PatientFactors(age=30, weight_kg=65.0, sex="F", pregnant=True)


@pytest.fixture
def renal_impaired_patient():
    return PatientFactors(age=72, weight_kg=80.0, sex="M", renal_function_egfr=25)


@pytest.fixture
def contrast_allergy_patient():
    return PatientFactors(
        age=55, weight_kg=75.0, sex="F",
        contrast_allergy=True, contrast_allergy_severity="severe",
    )


@pytest.fixture
def pediatric_patient():
    return PatientFactors(age=8, weight_kg=25.0, sex="M", pediatric=True)


@pytest.fixture
def obese_patient():
    return PatientFactors(age=60, weight_kg=130.0, sex="M", bmi=42.0)


# ═══════════════════════════════════════════════════════════════════════
# BASIC RECOMMENDATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestBasicRecommendations:
    """Test core recommendation logic without patient adjustments."""

    def test_recommend_headache(self, optimizer):
        rec = optimizer.recommend("headache")
        assert isinstance(rec, ProtocolRecommendation)
        assert rec.recommended_modality == "MRI"
        assert "Brain" in rec.recommended_protocol
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_lung_screening(self, optimizer):
        rec = optimizer.recommend("lung cancer screening")
        assert rec.recommended_modality == "CT"
        assert "Low-dose" in rec.recommended_protocol
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_stroke_acute(self, optimizer):
        rec = optimizer.recommend("stroke")
        assert rec.recommended_modality == "CT"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_abdominal_pain(self, optimizer):
        rec = optimizer.recommend("abdominal pain")
        assert rec.recommended_modality == "CT"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_breast_screening(self, optimizer):
        rec = optimizer.recommend("breast screening")
        assert rec.recommended_modality == "MAMMO"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_pulmonary_embolism(self, optimizer):
        rec = optimizer.recommend("pulmonary embolism")
        assert rec.recommended_modality == "CT"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_liver_lesion(self, optimizer):
        rec = optimizer.recommend("liver lesion")
        assert rec.recommended_modality == "MRI"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_thyroid_nodule(self, optimizer):
        rec = optimizer.recommend("thyroid nodule")
        assert rec.recommended_modality == "US"
        assert rec.acr_appropriateness_rating == 9

    def test_recommend_prostate_cancer(self, optimizer):
        rec = optimizer.recommend("prostate cancer")
        assert rec.recommended_modality == "MRI"
        assert "PI-RADS" in rec.recommended_protocol
        assert rec.acr_appropriateness_rating == 9


# ═══════════════════════════════════════════════════════════════════════
# DOSE ESTIMATION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestDoseEstimation:
    """Test that CT protocols include dose estimates."""

    def test_dose_estimate_included_ct(self, optimizer):
        rec = optimizer.recommend("stroke")
        assert rec.dose_estimate_msv is not None
        assert rec.dose_estimate_msv > 0

    def test_dose_estimate_ldct(self, optimizer):
        rec = optimizer.recommend("lung cancer screening")
        assert rec.dose_estimate_msv is not None
        assert rec.dose_estimate_msv <= 2.0  # Low-dose CT should be < 2 mSv

    def test_dose_estimate_none_for_mri(self, optimizer):
        rec = optimizer.recommend("headache")
        # MRI has no ionizing radiation — dose estimate may be None
        assert rec.dose_estimate_msv is None or rec.dose_estimate_msv == 0.0

    def test_dose_estimate_xray_low(self, optimizer):
        rec = optimizer.recommend("breast screening")
        # Mammography dose should be very low
        assert rec.dose_estimate_msv is not None
        assert rec.dose_estimate_msv < 1.0


# ═══════════════════════════════════════════════════════════════════════
# ALTERNATIVES TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestAlternatives:
    """Test that alternatives are provided with ratings."""

    def test_alternatives_provided(self, optimizer):
        rec = optimizer.recommend("headache")
        assert len(rec.alternatives) >= 1
        for alt in rec.alternatives:
            assert "modality" in alt
            assert "protocol" in alt
            assert "rating" in alt
            assert 1 <= alt["rating"] <= 9

    def test_alternatives_ranked_by_rating(self, optimizer):
        rec = optimizer.recommend("abdominal pain")
        assert len(rec.alternatives) >= 2
        # All alternatives should have rating <= top recommendation
        for alt in rec.alternatives:
            assert alt["rating"] <= rec.acr_appropriateness_rating

    def test_chest_pain_has_multiple_alternatives(self, optimizer):
        rec = optimizer.recommend("chest pain")
        assert len(rec.alternatives) >= 2


# ═══════════════════════════════════════════════════════════════════════
# PATIENT ADJUSTMENT TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestPregnancyAdjustments:
    """Test protocol adjustments for pregnant patients."""

    def test_recommend_with_pregnancy_ct_warning(self, optimizer, pregnant_patient):
        rec = optimizer.recommend("stroke", patient=pregnant_patient)
        assert any("PREGNANCY" in w for w in rec.warnings)

    def test_recommend_with_pregnancy_contrast_warning(self, optimizer, pregnant_patient):
        rec = optimizer.recommend("abdominal pain", patient=pregnant_patient)
        # CT A/P with contrast should warn about both radiation and contrast
        pregnancy_warnings = [w for w in rec.warnings if "PREGNANCY" in w]
        assert len(pregnancy_warnings) >= 1

    def test_pregnancy_no_warning_for_us(self, optimizer, pregnant_patient):
        rec = optimizer.recommend("thyroid nodule", patient=pregnant_patient)
        # Ultrasound is safe in pregnancy — no pregnancy warning
        pregnancy_warnings = [w for w in rec.warnings if "PREGNANCY" in w]
        assert len(pregnancy_warnings) == 0


class TestRenalImpairment:
    """Test protocol adjustments for renal impairment."""

    def test_recommend_with_renal_impairment(self, optimizer, renal_impaired_patient):
        rec = optimizer.recommend("abdominal pain", patient=renal_impaired_patient)
        assert any("RENAL" in w.upper() for w in rec.warnings)

    def test_egfr_below_30_warns_contrast(self, optimizer, renal_impaired_patient):
        rec = optimizer.recommend("liver lesion", patient=renal_impaired_patient)
        # eGFR < 30 should warn about both iodinated and gadolinium contrast
        renal_warnings = [w for w in rec.warnings if "RENAL" in w.upper()]
        assert len(renal_warnings) >= 1

    def test_egfr_moderate_caution(self, optimizer):
        patient = PatientFactors(age=65, renal_function_egfr=45)
        rec = optimizer.recommend("abdominal pain", patient=patient)
        assert any("RENAL CAUTION" in w for w in rec.warnings)


class TestContrastAllergy:
    """Test protocol adjustments for contrast allergy."""

    def test_recommend_with_contrast_allergy(self, optimizer, contrast_allergy_patient):
        rec = optimizer.recommend("abdominal pain", patient=contrast_allergy_patient)
        assert any("ALLERGY" in w.upper() for w in rec.warnings)

    def test_severe_allergy_contraindication(self, optimizer, contrast_allergy_patient):
        rec = optimizer.recommend("liver lesion", patient=contrast_allergy_patient)
        warnings_text = " ".join(rec.warnings).upper()
        assert "CONTRAINDICATED" in warnings_text or "SEVERE" in warnings_text

    def test_mild_allergy_premedication(self, optimizer):
        patient = PatientFactors(
            contrast_allergy=True, contrast_allergy_severity="mild",
        )
        rec = optimizer.recommend("abdominal pain", patient=patient)
        assert any("PREMEDICATION" in w.upper() or "ALLERGY" in w.upper() for w in rec.warnings)


class TestPediatricAdjustments:
    """Test protocol adjustments for pediatric patients."""

    def test_recommend_pediatric(self, optimizer, pediatric_patient):
        rec = optimizer.recommend("headache", patient=pediatric_patient)
        assert any("PEDIATRIC" in w for w in rec.warnings)

    def test_pediatric_dose_reduction(self, optimizer, pediatric_patient):
        adult_rec = optimizer.recommend("stroke")
        peds_rec = optimizer.recommend("stroke", patient=pediatric_patient)
        if adult_rec.dose_estimate_msv and peds_rec.dose_estimate_msv:
            assert peds_rec.dose_estimate_msv < adult_rec.dose_estimate_msv

    def test_pediatric_kv_reduction(self, optimizer, pediatric_patient):
        rec = optimizer.recommend("stroke", patient=pediatric_patient)
        if "kv" in rec.parameters:
            assert rec.parameters["kv"] <= 100

    def test_pediatric_protocol_label(self, optimizer, pediatric_patient):
        rec = optimizer.recommend("stroke", patient=pediatric_patient)
        assert "Pediatric" in rec.recommended_protocol or "PEDIATRIC" in " ".join(rec.warnings)

    def test_young_child_lower_kv(self, optimizer):
        young_child = PatientFactors(age=3, weight_kg=15.0, pediatric=True)
        rec = optimizer.recommend("abdominal pain", patient=young_child)
        if "kv" in rec.parameters:
            assert rec.parameters["kv"] <= 80


class TestBodyHabitus:
    """Test adjustments for large body habitus."""

    def test_obese_patient_kv_increase(self, optimizer, obese_patient):
        rec = optimizer.recommend("abdominal pain", patient=obese_patient)
        assert any("BMI" in w or "habitus" in w for w in rec.warnings)

    def test_heavy_patient_warning(self, optimizer):
        patient = PatientFactors(age=50, weight_kg=125.0)
        rec = optimizer.recommend("stroke", patient=patient)
        assert any("weight" in w.lower() or "120" in w for w in rec.warnings)


# ═══════════════════════════════════════════════════════════════════════
# INDICATION MATCHING TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestIndicationMatching:
    """Test fuzzy indication matching."""

    def test_recommend_unknown_indication(self, optimizer):
        rec = optimizer.recommend("xyzzy_unknown_condition_12345")
        assert isinstance(rec, ProtocolRecommendation)
        # Should return a default/best-effort recommendation
        assert rec.indication == "xyzzy_unknown_condition_12345"

    def test_synonym_matching_migraine(self, optimizer):
        rec = optimizer.recommend("migraine")
        assert rec.recommended_modality == "MRI"

    def test_synonym_matching_pe(self, optimizer):
        rec = optimizer.recommend("PE")
        assert rec.recommended_modality == "CT"

    def test_synonym_matching_chest_pain(self, optimizer):
        rec = optimizer.recommend("chest pain")
        assert rec.recommended_modality == "CT"

    def test_exact_key_match(self, optimizer):
        rec = optimizer.recommend("headache")
        assert rec.acr_appropriateness_rating == 9

    def test_available_indications(self, optimizer):
        indications = optimizer.get_available_indications()
        assert len(indications) >= 10
        assert "headache" in indications
        assert "stroke_acute" in indications
        assert "lung_cancer_screening" in indications


# ═══════════════════════════════════════════════════════════════════════
# RATIONALE AND CONTRAST TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestRationaleAndContrast:
    """Test rationale generation and contrast information."""

    def test_rationale_present(self, optimizer):
        rec = optimizer.recommend("headache")
        assert len(rec.rationale) > 0
        assert "ACR" in rec.rationale

    def test_contrast_info_for_contrast_protocol(self, optimizer):
        rec = optimizer.recommend("abdominal pain")
        assert rec.contrast is not None
        assert "agent" in rec.contrast

    def test_no_contrast_for_noncontrast_protocol(self, optimizer):
        rec = optimizer.recommend("headache")
        # MRI Brain without contrast — no contrast
        assert rec.contrast is None

    def test_contrast_agent_type(self, optimizer):
        rec = optimizer.recommend("liver lesion")
        assert rec.contrast is not None
        assert "gadolinium" in rec.contrast.get("agent", "").lower()
