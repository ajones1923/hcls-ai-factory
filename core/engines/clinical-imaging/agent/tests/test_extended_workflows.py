"""Tests for extended Imaging Intelligence Agent workflows.

Validates Breast BI-RADS, Thyroid TI-RADS, and Liver LI-RADS workflows
in mock mode, including classification logic, severity mapping,
cross-modal triggers, and recommendation text.

Author: Adam Jones
Date: April 2026
"""

import pytest

from src.models import (
    BiRADS,
    FindingCategory,
    FindingSeverity,
    LIRADS,
    TIRADS,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows import (
    WORKFLOW_REGISTRY,
    BaseImagingWorkflow,
    BreastBIRADSWorkflow,
    ThyroidTIRADSWorkflow,
    LiverLIRADSWorkflow,
)
from src.workflows.breast_birads import (
    birads_to_severity,
    birads_recommendation,
    classify_most_suspicious_finding,
)
from src.workflows.thyroid_tirads import (
    calculate_tirads_points,
    points_to_tirads,
    tirads_to_severity,
    tirads_fna_recommendation,
)
from src.workflows.liver_lirads import (
    classify_lirads,
    lirads_to_severity,
    lirads_recommendation,
)


# ===================================================================
# Breast BI-RADS Workflow
# ===================================================================


class TestBreastBIRADSWorkflow:
    """Tests for breast BI-RADS mock workflow."""

    def test_birads_init(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "breast_birads"
        assert wf.MODALITY == "mammography"
        assert wf.BODY_REGION == "breast"
        assert wf.mock_mode is True
        info = wf.get_workflow_info()
        assert info["name"] == "breast_birads"
        assert info["target_latency_sec"] == 300.0

    def test_birads_mock_run_completed(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        result = wf.run()
        assert isinstance(result, WorkflowResult)
        assert result.status == WorkflowStatus.COMPLETED
        assert result.is_mock is True
        assert result.inference_time_ms > 0

    def test_birads_mock_has_findings(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        result = wf.run()
        assert len(result.findings) >= 2
        # Should have mass and calcification findings
        categories = [f["category"] for f in result.findings]
        assert FindingCategory.MASS.value in categories
        assert FindingCategory.CALCIFICATION.value in categories

    def test_birads_classification(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # Mock has BI-RADS 4 findings -> overall BI-RADS 4
        assert "BI-RADS" in result.classification
        assert "4" in result.classification

    def test_birads_severity_mapping_normal(self):
        assert birads_to_severity(BiRADS.CAT_1) == FindingSeverity.NORMAL
        assert birads_to_severity(BiRADS.CAT_2) == FindingSeverity.NORMAL

    def test_birads_severity_mapping_significant(self):
        assert birads_to_severity(BiRADS.CAT_3) == FindingSeverity.SIGNIFICANT

    def test_birads_severity_mapping_urgent(self):
        assert birads_to_severity(BiRADS.CAT_4) == FindingSeverity.URGENT

    def test_birads_severity_mapping_critical(self):
        assert birads_to_severity(BiRADS.CAT_5) == FindingSeverity.CRITICAL
        assert birads_to_severity(BiRADS.CAT_6) == FindingSeverity.CRITICAL

    def test_birads_cross_modal_trigger(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # BI-RADS 4 should trigger cross-modal
        trigger_findings = [
            f for f in result.findings
            if f.get("cross_modal_trigger")
        ]
        assert len(trigger_findings) >= 1
        assert "BRCA1" in trigger_findings[0]["genomic_queries"]
        assert "BRCA2" in trigger_findings[0]["genomic_queries"]

    def test_birads_recommendation_4(self):
        rec = birads_recommendation(BiRADS.CAT_4)
        assert "tissue diagnosis" in rec.lower()

    def test_birads_recommendation_3(self):
        rec = birads_recommendation(BiRADS.CAT_3)
        assert "short-interval follow-up" in rec.lower() or "6 months" in rec.lower()

    def test_birads_recommendation_5(self):
        rec = birads_recommendation(BiRADS.CAT_5)
        assert "malignancy" in rec.lower()

    def test_birads_no_findings_returns_normal(self):
        wf = BreastBIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"findings": []})
        assert result.severity == FindingSeverity.NORMAL
        assert "BI-RADS 1" in result.classification

    def test_birads_classify_most_suspicious(self):
        findings = [
            {"birads_score": 2},
            {"birads_score": 4},
            {"birads_score": 3},
        ]
        assert classify_most_suspicious_finding(findings) == BiRADS.CAT_4

    def test_birads_classify_single_benign(self):
        findings = [{"birads_score": 2}]
        assert classify_most_suspicious_finding(findings) == BiRADS.CAT_2


# ===================================================================
# Thyroid TI-RADS Workflow
# ===================================================================


class TestThyroidTIRADSWorkflow:
    """Tests for thyroid TI-RADS mock workflow."""

    def test_tirads_init(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "thyroid_tirads"
        assert wf.MODALITY == "ultrasound"
        assert wf.BODY_REGION == "neck"
        assert wf.mock_mode is True
        info = wf.get_workflow_info()
        assert info["name"] == "thyroid_tirads"
        assert info["target_latency_sec"] == 180.0

    def test_tirads_mock_run_completed(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.run()
        assert isinstance(result, WorkflowResult)
        assert result.status == WorkflowStatus.COMPLETED
        assert result.is_mock is True
        assert result.inference_time_ms > 0

    def test_tirads_point_calculation_tr1(self):
        # All benign features: 0 points -> TR1
        nodule = {
            "composition": "cystic",
            "echogenicity": "anechoic",
            "shape": "wider_than_tall",
            "margin": "smooth",
            "echogenic_foci": "none",
        }
        points = calculate_tirads_points(nodule)
        assert points == 0
        assert points_to_tirads(points) == TIRADS.TR1

    def test_tirads_point_calculation_tr3(self):
        # solid(2) + isoechoic(1) = 3 points -> TR3
        nodule = {
            "composition": "solid",
            "echogenicity": "isoechoic",
            "shape": "wider_than_tall",
            "margin": "smooth",
            "echogenic_foci": "none",
        }
        points = calculate_tirads_points(nodule)
        assert points == 3
        assert points_to_tirads(points) == TIRADS.TR3

    def test_tirads_point_calculation_tr4(self):
        # solid(2) + hypoechoic(2) + smooth(0) = 4 points -> TR4
        # Also test with 5 points: solid(2) + hypoechoic(2) + macrocalc(1) = 5
        nodule = {
            "composition": "solid",
            "echogenicity": "hypoechoic",
            "shape": "wider_than_tall",
            "margin": "smooth",
            "echogenic_foci": "macrocalcifications",
        }
        points = calculate_tirads_points(nodule)
        assert points == 5
        assert points_to_tirads(points) == TIRADS.TR4

    def test_tirads_point_calculation_tr5(self):
        # solid(2) + hypoechoic(2) + taller(3) + irregular(2) + punctate(3) = 12 -> TR5
        nodule = {
            "composition": "solid",
            "echogenicity": "hypoechoic",
            "shape": "taller_than_wide",
            "margin": "irregular",
            "echogenic_foci": "punctate_echogenic_foci",
        }
        points = calculate_tirads_points(nodule)
        assert points >= 7
        assert points_to_tirads(points) == TIRADS.TR5

    def test_tirads_fna_recommendation_tr5_large(self):
        rec = tirads_fna_recommendation(TIRADS.TR5, 15.0)
        assert "fna" in rec.lower()

    def test_tirads_fna_recommendation_tr5_small(self):
        rec = tirads_fna_recommendation(TIRADS.TR5, 7.0)
        assert "follow-up" in rec.lower()

    def test_tirads_fna_recommendation_tr4_large(self):
        rec = tirads_fna_recommendation(TIRADS.TR4, 20.0)
        assert "fna" in rec.lower()

    def test_tirads_fna_recommendation_tr3_large(self):
        rec = tirads_fna_recommendation(TIRADS.TR3, 30.0)
        assert "fna" in rec.lower()

    def test_tirads_fna_recommendation_tr1(self):
        rec = tirads_fna_recommendation(TIRADS.TR1, 50.0)
        assert "no fna" in rec.lower()

    def test_tirads_cross_modal_trigger(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # Mock has a TR5 nodule -> should trigger cross-modal
        trigger_findings = [
            f for f in result.findings
            if f.get("cross_modal_trigger")
        ]
        assert len(trigger_findings) >= 1
        assert "BRAF V600E" in trigger_findings[0]["genomic_queries"]

    def test_tirads_multiple_nodules(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # Mock data has 2 nodules
        assert len(result.findings) == 2
        assert result.measurements["nodule_count"] == 2.0
        # Each finding should have tirads_points
        for finding in result.findings:
            assert "tirads_points" in finding
            assert "tirads" in finding

    def test_tirads_severity_mapping(self):
        assert tirads_to_severity(TIRADS.TR1) == FindingSeverity.NORMAL
        assert tirads_to_severity(TIRADS.TR2) == FindingSeverity.NORMAL
        assert tirads_to_severity(TIRADS.TR3) == FindingSeverity.SIGNIFICANT
        assert tirads_to_severity(TIRADS.TR4) == FindingSeverity.URGENT
        assert tirads_to_severity(TIRADS.TR5) == FindingSeverity.CRITICAL

    def test_tirads_no_nodules_returns_normal(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [], "thyroid_volume_ml": 15.0})
        assert result.severity == FindingSeverity.NORMAL
        assert "TR1" in result.classification

    def test_tirads_points_to_tr2(self):
        # mixed(1) + isoechoic(1) = 2 -> TR2
        assert points_to_tirads(2) == TIRADS.TR2

    def test_tirads_mock_classification_includes_tr(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.run()
        assert "TI-RADS" in result.classification


# ===================================================================
# Liver LI-RADS Workflow
# ===================================================================


class TestLiverLIRADSWorkflow:
    """Tests for liver LI-RADS mock workflow."""

    def test_lirads_init(self):
        wf = LiverLIRADSWorkflow(mock_mode=True)
        assert wf.WORKFLOW_NAME == "liver_lirads"
        assert wf.MODALITY == "ct"
        assert wf.BODY_REGION == "abdomen"
        assert wf.mock_mode is True
        info = wf.get_workflow_info()
        assert info["name"] == "liver_lirads"
        assert info["target_latency_sec"] == 300.0

    def test_lirads_mock_run_completed(self):
        wf = LiverLIRADSWorkflow(mock_mode=True)
        result = wf.run()
        assert isinstance(result, WorkflowResult)
        assert result.status == WorkflowStatus.COMPLETED
        assert result.is_mock is True
        assert result.inference_time_ms > 0

    def test_lirads_classification(self):
        wf = LiverLIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # Mock has a 25mm observation with APHE + washout + capsule -> LR-5
        assert "LI-RADS" in result.classification
        assert "LR-5" in result.classification

    def test_lirads_lr5_criteria(self):
        """APHE + >= 20mm + washout -> LR-5."""
        obs = {
            "aphe": True,
            "size_mm": 25.0,
            "washout": True,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_5

    def test_lirads_lr5_with_capsule(self):
        """APHE + >= 20mm + capsule -> LR-5."""
        obs = {
            "aphe": True,
            "size_mm": 22.0,
            "washout": False,
            "capsule": True,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_5

    def test_lirads_lr4_aphe_large_no_features(self):
        """APHE + >= 20mm + 0 additional features -> LR-4."""
        obs = {
            "aphe": True,
            "size_mm": 20.0,
            "washout": False,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_4

    def test_lirads_lr4_aphe_small_with_feature(self):
        """APHE + < 10mm + 1 feature -> LR-4."""
        obs = {
            "aphe": True,
            "size_mm": 8.0,
            "washout": True,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_4

    def test_lirads_lr3_aphe_small_no_features(self):
        """APHE + < 10mm + 0 features -> LR-3."""
        obs = {
            "aphe": True,
            "size_mm": 8.0,
            "washout": False,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_3

    def test_lirads_lr3_no_aphe(self):
        """No APHE + < 20mm -> LR-3."""
        obs = {
            "aphe": False,
            "size_mm": 15.0,
            "washout": True,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_3

    def test_lirads_cross_modal_trigger(self):
        wf = LiverLIRADSWorkflow(mock_mode=True)
        result = wf.run()
        # LR-5 should trigger cross-modal
        trigger_findings = [
            f for f in result.findings
            if f.get("cross_modal_trigger")
        ]
        assert len(trigger_findings) >= 1
        assert "TP53" in trigger_findings[0]["genomic_queries"]
        assert "CTNNB1" in trigger_findings[0]["genomic_queries"]

    def test_lirads_severity_mapping(self):
        assert lirads_to_severity(LIRADS.LR_1) == FindingSeverity.NORMAL
        assert lirads_to_severity(LIRADS.LR_2) == FindingSeverity.ROUTINE
        assert lirads_to_severity(LIRADS.LR_3) == FindingSeverity.SIGNIFICANT
        assert lirads_to_severity(LIRADS.LR_4) == FindingSeverity.URGENT
        assert lirads_to_severity(LIRADS.LR_5) == FindingSeverity.CRITICAL
        assert lirads_to_severity(LIRADS.LR_M) == FindingSeverity.CRITICAL
        assert lirads_to_severity(LIRADS.LR_TIV) == FindingSeverity.CRITICAL

    def test_lirads_tiv_detection(self):
        """Tumor in vein -> LR-TIV regardless of other features."""
        obs = {
            "aphe": True,
            "size_mm": 15.0,
            "washout": True,
            "capsule": True,
            "threshold_growth": False,
            "tumor_in_vein": True,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_TIV

    def test_lirads_lrm_targetoid(self):
        """Targetoid appearance -> LR-M."""
        obs = {
            "aphe": False,
            "size_mm": 30.0,
            "washout": False,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": True,
        }
        assert classify_lirads(obs) == LIRADS.LR_M

    def test_lirads_no_observations_returns_normal(self):
        wf = LiverLIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"observations": [], "liver_volume_ml": 1500.0})
        assert result.severity == FindingSeverity.NORMAL
        assert "LR-1" in result.classification

    def test_lirads_recommendation_lr5(self):
        rec = lirads_recommendation(LIRADS.LR_5)
        assert "hcc" in rec.lower() or "multidisciplinary" in rec.lower()

    def test_lirads_recommendation_lr_tiv(self):
        rec = lirads_recommendation(LIRADS.LR_TIV)
        assert "tumor" in rec.lower() or "thrombosis" in rec.lower()

    def test_lirads_lr5_medium_two_features(self):
        """APHE + 10-19mm + 2 features -> LR-5."""
        obs = {
            "aphe": True,
            "size_mm": 15.0,
            "washout": True,
            "capsule": True,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_5

    def test_lirads_lr4_medium_one_feature(self):
        """APHE + 10-19mm + 1 feature -> LR-4."""
        obs = {
            "aphe": True,
            "size_mm": 15.0,
            "washout": True,
            "capsule": False,
            "threshold_growth": False,
            "tumor_in_vein": False,
            "targetoid_appearance": False,
        }
        assert classify_lirads(obs) == LIRADS.LR_4


# ===================================================================
# Registry tests
# ===================================================================


class TestExtendedRegistry:
    """Tests for the updated WORKFLOW_REGISTRY."""

    def test_registry_has_9_entries(self):
        assert len(WORKFLOW_REGISTRY) == 9

    def test_all_new_workflows_registered(self):
        assert "breast_birads" in WORKFLOW_REGISTRY
        assert "thyroid_tirads" in WORKFLOW_REGISTRY
        assert "liver_lirads" in WORKFLOW_REGISTRY

    def test_new_workflows_are_subclasses(self):
        for name in ("breast_birads", "thyroid_tirads", "liver_lirads"):
            assert issubclass(WORKFLOW_REGISTRY[name], BaseImagingWorkflow), (
                f"{name} is not a BaseImagingWorkflow subclass"
            )

    def test_breast_birads_class_in_registry(self):
        assert WORKFLOW_REGISTRY["breast_birads"] is BreastBIRADSWorkflow

    def test_thyroid_tirads_class_in_registry(self):
        assert WORKFLOW_REGISTRY["thyroid_tirads"] is ThyroidTIRADSWorkflow

    def test_liver_lirads_class_in_registry(self):
        assert WORKFLOW_REGISTRY["liver_lirads"] is LiverLIRADSWorkflow
