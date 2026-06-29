"""Tests for Imaging Intelligence Agent reference workflows.

Validates WORKFLOW_REGISTRY, all four workflows in mock mode,
workflow results structure (findings, measurements, severity),
and get_workflow_info() metadata.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.models import (
    FindingCategory,
    FindingSeverity,
    LungRADS,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows import (
    WORKFLOW_REGISTRY,
    BaseImagingWorkflow,
    CTChestLungNoduleWorkflow,
    CTHeadHemorrhageWorkflow,
    CXRRapidFindingsWorkflow,
    MRIBrainMSLesionWorkflow,
)
from src.workflows.ct_chest_lung_nodule import (
    classify_solid_nodule_lung_rads,
    lung_rads_to_severity,
    lung_rads_recommendation,
)
from src.workflows.mri_brain_ms_lesion import (
    classify_ms_activity,
    ms_activity_to_severity,
    ms_activity_recommendation,
)


# ===================================================================
# WORKFLOW_REGISTRY
# ===================================================================


class TestWorkflowRegistry:
    """Tests for the WORKFLOW_REGISTRY."""

    def test_has_4_entries(self):
        assert len(WORKFLOW_REGISTRY) == 4

    def test_ct_head_hemorrhage_registered(self):
        assert "ct_head_hemorrhage" in WORKFLOW_REGISTRY

    def test_ct_chest_lung_nodule_registered(self):
        assert "ct_chest_lung_nodule" in WORKFLOW_REGISTRY

    def test_cxr_rapid_findings_registered(self):
        assert "cxr_rapid_findings" in WORKFLOW_REGISTRY

    def test_mri_brain_ms_lesion_registered(self):
        assert "mri_brain_ms_lesion" in WORKFLOW_REGISTRY

    def test_all_values_are_workflow_subclasses(self):
        for name, cls in WORKFLOW_REGISTRY.items():
            assert issubclass(cls, BaseImagingWorkflow), (
                f"{name} is not a BaseImagingWorkflow subclass"
            )


# ===================================================================
# CTHeadHemorrhageWorkflow
# ===================================================================


class TestCTHeadHemorrhageWorkflow:
    """Tests for CT head hemorrhage mock workflow."""

    def test_init(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "ct_head_hemorrhage"
        assert wf.MODALITY == "ct"
        assert wf.BODY_REGION == "head"

    def test_mock_run_returns_completed(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        assert isinstance(result, WorkflowResult)
        assert result.status == WorkflowStatus.COMPLETED

    def test_mock_run_detects_hemorrhage(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        assert len(result.findings) >= 1
        assert result.findings[0]["category"] == FindingCategory.HEMORRHAGE.value

    def test_mock_run_has_measurements(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        assert "volume_ml" in result.measurements
        assert "midline_shift_mm" in result.measurements
        assert "max_thickness_mm" in result.measurements

    def test_mock_run_severity_is_urgent(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        # Mock data: 12.5 mL > 5 mL threshold -> URGENT
        assert result.severity == FindingSeverity.URGENT

    def test_mock_run_is_mock_flag(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        assert result.is_mock is True

    def test_mock_run_has_inference_time(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.run()
        assert result.inference_time_ms > 0

    def test_get_workflow_info(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        info = wf.get_workflow_info()
        assert info["name"] == "ct_head_hemorrhage"
        assert info["modality"] == "ct"
        assert info["body_region"] == "head"
        assert info["mock_mode"] is True

    def test_postprocess_no_hemorrhage(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.postprocess({"hemorrhage_detected": False})
        assert result.severity == FindingSeverity.NORMAL
        assert result.classification == "no_hemorrhage"

    def test_postprocess_critical_volume(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.postprocess({
            "hemorrhage_detected": True,
            "hemorrhage_type": "subdural",
            "location": "left convexity",
            "volume_ml": 35.0,
            "midline_shift_mm": 2.0,
            "max_thickness_mm": 8.0,
        })
        assert result.severity == FindingSeverity.CRITICAL

    def test_postprocess_critical_shift(self):
        wf = CTHeadHemorrhageWorkflow(mock_mode=True)
        result = wf.postprocess({
            "hemorrhage_detected": True,
            "hemorrhage_type": "epidural",
            "location": "right temporal",
            "volume_ml": 10.0,
            "midline_shift_mm": 7.0,
            "max_thickness_mm": 5.0,
        })
        assert result.severity == FindingSeverity.CRITICAL


# ===================================================================
# CTChestLungNoduleWorkflow
# ===================================================================


class TestCTChestLungNoduleWorkflow:
    """Tests for CT chest lung nodule mock workflow."""

    def test_init(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "ct_chest_lung_nodule"
        assert wf.MODALITY == "ct"
        assert wf.BODY_REGION == "chest"

    def test_mock_run_returns_completed(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.run()
        assert result.status == WorkflowStatus.COMPLETED

    def test_mock_run_detects_nodules(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.run()
        assert len(result.findings) >= 1
        assert result.findings[0]["category"] == FindingCategory.NODULE.value

    def test_mock_run_has_lung_rads(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.run()
        for finding in result.findings:
            assert "lung_rads" in finding

    def test_mock_run_has_measurements(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.run()
        assert "nodule_count" in result.measurements

    def test_mock_run_classification_includes_lung_rads(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.run()
        assert "Lung-RADS" in result.classification

    def test_get_workflow_info(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        info = wf.get_workflow_info()
        assert info["name"] == "ct_chest_lung_nodule"
        assert "RetinaNet" in str(info["models_used"])


class TestLungRADSClassification:
    """Tests for Lung-RADS helper functions."""

    def test_small_nodule_cat_2(self):
        assert classify_solid_nodule_lung_rads(4.0) == LungRADS.CAT_2

    def test_medium_nodule_cat_3(self):
        assert classify_solid_nodule_lung_rads(7.0) == LungRADS.CAT_3

    def test_suspicious_nodule_cat_4a(self):
        assert classify_solid_nodule_lung_rads(12.0) == LungRADS.CAT_4A

    def test_large_nodule_cat_4b(self):
        assert classify_solid_nodule_lung_rads(18.0) == LungRADS.CAT_4B

    def test_lung_rads_to_severity_cat_1(self):
        assert lung_rads_to_severity(LungRADS.CAT_1) == FindingSeverity.NORMAL

    def test_lung_rads_to_severity_cat_4b(self):
        assert lung_rads_to_severity(LungRADS.CAT_4B) == FindingSeverity.CRITICAL

    def test_lung_rads_recommendation_cat_2(self):
        rec = lung_rads_recommendation(LungRADS.CAT_2)
        assert "annual" in rec.lower()

    def test_lung_rads_recommendation_cat_4b(self):
        rec = lung_rads_recommendation(LungRADS.CAT_4B)
        assert "tissue sampling" in rec.lower() or "pet" in rec.lower()


# ===================================================================
# CXRRapidFindingsWorkflow
# ===================================================================


class TestCXRRapidFindingsWorkflow:
    """Tests for CXR rapid findings mock workflow."""

    def test_init(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "cxr_rapid_findings"
        assert wf.MODALITY == "cxr"
        assert wf.BODY_REGION == "chest"

    def test_mock_run_returns_completed(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()
        assert result.status == WorkflowStatus.COMPLETED

    def test_mock_run_detects_consolidation(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()
        # Mock data has consolidation: 0.87 (above 0.60 threshold)
        finding_categories = [f["category"] for f in result.findings]
        assert FindingCategory.CONSOLIDATION.value in finding_categories

    def test_mock_run_detects_effusion(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()
        # Mock data has pleural_effusion: 0.72 (above 0.55 threshold)
        finding_categories = [f["category"] for f in result.findings]
        assert FindingCategory.EFFUSION.value in finding_categories

    def test_mock_run_has_confidence_measurements(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()
        assert "pneumothorax_confidence" in result.measurements
        assert "consolidation_confidence" in result.measurements

    def test_mock_run_severity_is_at_least_urgent(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.run()
        # Consolidation has URGENT severity
        severity_order = [
            FindingSeverity.NORMAL, FindingSeverity.ROUTINE,
            FindingSeverity.SIGNIFICANT, FindingSeverity.URGENT,
            FindingSeverity.CRITICAL,
        ]
        assert severity_order.index(result.severity) >= severity_order.index(FindingSeverity.URGENT)

    def test_get_workflow_info(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        info = wf.get_workflow_info()
        assert info["name"] == "cxr_rapid_findings"
        assert info["target_latency_sec"] == 30.0


# ===================================================================
# MRIBrainMSLesionWorkflow
# ===================================================================


class TestMRIBrainMSLesionWorkflow:
    """Tests for MRI brain MS lesion mock workflow."""

    def test_init(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "mri_brain_ms_lesion"
        assert wf.MODALITY == "mri"
        assert wf.BODY_REGION == "brain"

    def test_mock_run_returns_completed(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        assert result.status == WorkflowStatus.COMPLETED

    def test_mock_run_detects_lesions(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        assert len(result.findings) >= 1
        assert result.findings[0]["category"] == FindingCategory.LESION.value

    def test_mock_run_has_lesion_measurements(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        assert "total_lesion_count" in result.measurements
        assert "total_lesion_volume_ml" in result.measurements
        assert "new_lesion_count" in result.measurements

    def test_mock_run_classification_includes_ms(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        assert "ms_" in result.classification

    def test_mock_run_severity_is_significant(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        # Mock data: 2 new + 1 enlarging = 3 active lesions -> active -> SIGNIFICANT
        assert result.severity == FindingSeverity.SIGNIFICANT

    def test_mock_run_has_per_lesion_findings(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        result = wf.run()
        # Should have per-lesion findings for new/enhancing lesions
        new_enhancing = [
            f for f in result.findings
            if "NEW" in f.get("description", "") or "ENHANCING" in f.get("description", "")
        ]
        assert len(new_enhancing) >= 1

    def test_get_workflow_info(self):
        wf = MRIBrainMSLesionWorkflow(mock_mode=True)
        info = wf.get_workflow_info()
        assert info["name"] == "mri_brain_ms_lesion"
        assert info["modality"] == "mri"


class TestMSActivityClassification:
    """Tests for MS activity helper functions."""

    def test_stable_with_0_new(self):
        assert classify_ms_activity(0) == "stable"

    def test_active_with_2_new(self):
        assert classify_ms_activity(2) == "active"

    def test_highly_active_with_5_new(self):
        assert classify_ms_activity(5) == "highly_active"

    def test_activity_to_severity_stable(self):
        assert ms_activity_to_severity("stable") == FindingSeverity.ROUTINE

    def test_activity_to_severity_active(self):
        assert ms_activity_to_severity("active") == FindingSeverity.SIGNIFICANT

    def test_activity_to_severity_highly_active(self):
        assert ms_activity_to_severity("highly_active") == FindingSeverity.URGENT

    def test_recommendation_stable(self):
        rec = ms_activity_recommendation("stable")
        assert "12 months" in rec or "routine" in rec.lower()

    def test_recommendation_highly_active(self):
        rec = ms_activity_recommendation("highly_active")
        assert "escalation" in rec.lower() or "higher-efficacy" in rec.lower()
