"""Tests for all enums and Pydantic models in src/models.py.

Covers:
  - Enum member counts and values
  - NeuroWorkflowType members
  - SeverityLevel members
  - EvidenceLevel members
  - ClinicalScaleType members
  - StrokeType, DementiaSubtype, SeizureType, etc.
  - Pydantic model validation (NeuroQuery, WorkflowResult, ScaleResult, etc.)

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    ATNStage,
    ClinicalScaleType,
    DementiaSubtype,
    EvidenceLevel,
    HeadacheType,
    NMJPattern,
    NeuroWorkflowType,
    ScaleResult,
    SeizureType,
    SeverityLevel,
    StrokeType,
    WorkflowResult,
    NeuroQuery,
    NeuroResponse,
    NeuroSearchResult,
)


# ===================================================================
# ENUM TESTS
# ===================================================================


class TestNeuroWorkflowType:
    """Tests for NeuroWorkflowType enum."""

    def test_member_count(self):
        """NeuroWorkflowType must have exactly 9 members."""
        assert len(NeuroWorkflowType) == 9

    def test_acute_stroke_value(self):
        assert NeuroWorkflowType.ACUTE_STROKE.value == "acute_stroke"

    def test_dementia_evaluation_value(self):
        assert NeuroWorkflowType.DEMENTIA_EVALUATION.value == "dementia_evaluation"

    def test_epilepsy_focus_value(self):
        assert NeuroWorkflowType.EPILEPSY_FOCUS.value == "epilepsy_focus"

    def test_brain_tumor_value(self):
        assert NeuroWorkflowType.BRAIN_TUMOR.value == "brain_tumor"

    def test_ms_monitoring_value(self):
        assert NeuroWorkflowType.MS_MONITORING.value == "ms_monitoring"

    def test_parkinsons_assessment_value(self):
        assert NeuroWorkflowType.PARKINSONS_ASSESSMENT.value == "parkinsons_assessment"

    def test_headache_classification_value(self):
        assert NeuroWorkflowType.HEADACHE_CLASSIFICATION.value == "headache_classification"

    def test_neuromuscular_evaluation_value(self):
        assert NeuroWorkflowType.NEUROMUSCULAR_EVALUATION.value == "neuromuscular_evaluation"

    def test_general_value(self):
        assert NeuroWorkflowType.GENERAL.value == "general"

    def test_all_members_are_strings(self):
        for member in NeuroWorkflowType:
            assert isinstance(member.value, str)

    def test_from_value_roundtrip(self):
        for member in NeuroWorkflowType:
            assert NeuroWorkflowType(member.value) == member

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            NeuroWorkflowType("nonexistent_workflow")

    def test_all_values_unique(self):
        values = [m.value for m in NeuroWorkflowType]
        assert len(values) == len(set(values))

    def test_str_enum_behavior(self):
        assert NeuroWorkflowType.GENERAL == "general"


class TestSeverityLevel:
    """Tests for SeverityLevel enum."""

    def test_critical_value(self):
        assert SeverityLevel.CRITICAL.value == "critical"

    def test_informational_value(self):
        assert SeverityLevel.INFORMATIONAL.value == "informational"

    def test_has_moderate(self):
        assert hasattr(SeverityLevel, "MODERATE")

    def test_has_low(self):
        assert hasattr(SeverityLevel, "LOW")

    def test_member_count(self):
        assert len(SeverityLevel) == 5

    def test_all_values_unique(self):
        values = [m.value for m in SeverityLevel]
        assert len(values) == len(set(values))


class TestEvidenceLevel:
    """Tests for EvidenceLevel enum."""

    def test_class_i_value(self):
        assert EvidenceLevel.CLASS_I.value == "class_I"

    def test_class_iv_value(self):
        assert EvidenceLevel.CLASS_IV.value == "class_IV"

    def test_member_count(self):
        assert len(EvidenceLevel) == 4


class TestClinicalScaleType:
    """Tests for ClinicalScaleType enum."""

    def test_nihss_value(self):
        assert ClinicalScaleType.NIHSS.value == "nihss"

    def test_gcs_value(self):
        assert ClinicalScaleType.GCS.value == "gcs"

    def test_moca_value(self):
        assert ClinicalScaleType.MOCA.value == "moca"

    def test_edss_value(self):
        assert ClinicalScaleType.EDSS.value == "edss"

    def test_member_count(self):
        assert len(ClinicalScaleType) == 10

    def test_all_values_unique(self):
        values = [m.value for m in ClinicalScaleType]
        assert len(values) == len(set(values))


class TestStrokeType:
    """Tests for StrokeType enum."""

    def test_has_ischemic(self):
        assert hasattr(StrokeType, "ISCHEMIC")

    def test_has_hemorrhagic(self):
        assert hasattr(StrokeType, "HEMORRHAGIC")

    def test_has_tia(self):
        assert hasattr(StrokeType, "TIA")

    def test_has_sah(self):
        assert hasattr(StrokeType, "SAH")


class TestDementiaSubtype:
    """Tests for DementiaSubtype enum."""

    def test_has_alzheimers(self):
        assert hasattr(DementiaSubtype, "ALZHEIMERS")

    def test_has_lewy_body(self):
        assert hasattr(DementiaSubtype, "LEWY_BODY")

    def test_has_frontotemporal(self):
        assert hasattr(DementiaSubtype, "FRONTOTEMPORAL")


class TestSeizureType:
    """Tests for SeizureType enum."""

    def test_has_focal_aware(self):
        assert hasattr(SeizureType, "FOCAL_AWARE")

    def test_has_generalized_tonic_clonic(self):
        assert hasattr(SeizureType, "GENERALIZED_TONIC_CLONIC")


class TestATNStage:
    """Tests for ATN staging enum."""

    def test_has_all_negative(self):
        assert hasattr(ATNStage, "A_NEG_T_NEG_N_NEG")

    def test_has_all_positive(self):
        assert hasattr(ATNStage, "A_POS_T_POS_N_POS")


class TestHeadacheType:
    """Tests for HeadacheType enum."""

    def test_has_migraine_without_aura(self):
        assert hasattr(HeadacheType, "MIGRAINE_WITHOUT_AURA")

    def test_has_cluster(self):
        assert hasattr(HeadacheType, "CLUSTER")


class TestNMJPattern:
    """Tests for NMJPattern enum."""

    def test_has_demyelinating(self):
        assert hasattr(NMJPattern, "DEMYELINATING")

    def test_has_myopathic(self):
        assert hasattr(NMJPattern, "MYOPATHIC")


# ===================================================================
# PYDANTIC MODEL TESTS
# ===================================================================


class TestScaleResult:
    """Tests for ScaleResult Pydantic model."""

    def test_valid_creation(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=8.0,
            max_score=42.0,
            interpretation="Moderate stroke",
            severity_category="moderate",
        )
        assert result.score == 8.0
        assert result.max_score == 42.0

    def test_with_recommendations(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.GCS,
            score=7.0,
            max_score=15.0,
            interpretation="Severe",
            severity_category="severe",
            recommendations=["Consider intubation", "ICU admission"],
        )
        assert len(result.recommendations) == 2


class TestWorkflowResult:
    """Tests for WorkflowResult Pydantic model."""

    def test_valid_creation(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["Large vessel occlusion"],
            recommendations=["Proceed to thrombectomy"],
        )
        assert result.workflow_type == NeuroWorkflowType.ACUTE_STROKE
        assert len(result.findings) == 1

    def test_default_severity(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.GENERAL,
        )
        assert result.severity == SeverityLevel.INFORMATIONAL

    def test_empty_lists_by_default(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.GENERAL,
        )
        assert result.findings == []
        assert result.recommendations == []
        assert result.guideline_references == []
        assert result.cross_modal_triggers == []


class TestNeuroQuery:
    """Tests for NeuroQuery Pydantic model."""

    def test_valid_creation(self):
        query = NeuroQuery(
            query="Evaluate acute left MCA stroke",
        )
        assert query.query == "Evaluate acute left MCA stroke"

    def test_with_workflow_type(self):
        query = NeuroQuery(
            query="Assess patient for Parkinson's disease",
            workflow=NeuroWorkflowType.PARKINSONS_ASSESSMENT,
        )
        assert query.workflow == NeuroWorkflowType.PARKINSONS_ASSESSMENT


class TestNeuroSearchResult:
    """Tests for NeuroSearchResult Pydantic model."""

    def test_valid_creation(self):
        result = NeuroSearchResult(
            collection="neuro_cerebrovascular",
            content="Large vessel occlusion stroke management",
            score=0.85,
        )
        assert result.score == 0.85
        assert result.collection == "neuro_cerebrovascular"


class TestNeuroResponse:
    """Tests for NeuroResponse Pydantic model."""

    def test_valid_creation(self):
        wf_result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["LVO detected"],
        )
        response = NeuroResponse(
            query="Evaluate stroke",
            workflow=NeuroWorkflowType.ACUTE_STROKE,
            result=wf_result,
        )
        assert response.workflow == NeuroWorkflowType.ACUTE_STROKE

    def test_has_expected_fields(self):
        wf_result = WorkflowResult(
            workflow_type=NeuroWorkflowType.GENERAL,
        )
        response = NeuroResponse(
            query="test query",
            workflow=NeuroWorkflowType.GENERAL,
            result=wf_result,
        )
        assert response.citations == []
        assert response.warnings == []
