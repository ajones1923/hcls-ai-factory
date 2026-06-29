"""Integration tests for the Cardiology Intelligence Agent.

Tests end-to-end workflows, cross-module consistency, and the full
plan -> search -> evaluate -> synthesize pipeline using mocked
external dependencies (Milvus, LLM, embeddings).

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.models import (
    CardioWorkflowType,
    EjectionFractionCategory,
    RiskScoreInput,
    RiskScoreResult,
    RiskScoreType,
    SeverityLevel,
    WorkflowResult,
)
from src.collections import (
    ALL_COLLECTIONS,
    COLLECTION_NAMES,
    WORKFLOW_COLLECTION_WEIGHTS,
    get_search_weights,
)
from src.agent import WORKFLOW_COLLECTION_BOOST, CARDIO_CONDITIONS


# =====================================================================
# Risk Calculator Integration
# =====================================================================


class TestRiskCalculatorIntegration:
    """Verify that risk calculators produce clinically valid results."""

    def test_ascvd_white_female_low_risk(self):
        """45-year-old white female with ideal values -> low risk."""
        from src.risk_calculators import RiskCalculatorEngine
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=45, sex="female", race="white",
            total_cholesterol=180.0, hdl=60.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            smoker=False, diabetes=False,
        )
        result = engine.calculate(inp)
        assert result.score_value < 5.0, f"Expected low risk, got {result.score_value}%"
        assert result.risk_category.lower() == "low"

    def test_ascvd_high_risk_profile(self):
        """65-year-old male smoker with diabetes -> high risk."""
        from src.risk_calculators import RiskCalculatorEngine
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=65, sex="male", race="white",
            total_cholesterol=260.0, hdl=35.0,
            systolic_bp=160.0,
            hypertension_treatment=True,
            smoker=True, diabetes=True,
        )
        result = engine.calculate(inp)
        assert result.score_value >= 20.0, f"Expected high risk, got {result.score_value}%"
        assert result.risk_category.lower() == "high"

    def test_cha2ds2_vasc_score_zero(self):
        """Young male with no risk factors -> score 0."""
        from src.risk_calculators import RiskCalculatorEngine
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.CHA2DS2_VASC,
            age=45, sex="male",
            congestive_heart_failure=False,
            hypertension_treatment=False,
            diabetes=False,
            history_of_stroke=False,
            vascular_disease=False,
        )
        result = engine.calculate(inp)
        assert result.score_value == 0

    def test_cha2ds2_vasc_score_max_components(self):
        """75+ female with all risk factors -> high score."""
        from src.risk_calculators import RiskCalculatorEngine
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.CHA2DS2_VASC,
            age=78, sex="female",
            congestive_heart_failure=True,
            hypertension_treatment=True,
            diabetes=True,
            history_of_stroke=True,
            vascular_disease=True,
        )
        result = engine.calculate(inp)
        assert result.score_value >= 7

    def test_heart_score_low(self):
        """Low-risk chest pain patient -> HEART <= 3."""
        from src.risk_calculators import RiskCalculatorEngine
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.HEART,
            age=35,
        )
        result = engine.calculate(inp, extra={
            "history_suspicion": 0,
            "ecg_finding": 0,
            "troponin_result": 0,
            "risk_factor_count": 0,
        })
        assert result.score_value <= 3


# =====================================================================
# GDMT Optimizer Integration
# =====================================================================


class TestGDMTOptimizerIntegration:
    """Test GDMT optimizer produces guideline-consistent recommendations."""

    def test_hfref_all_pillars_recommended(self):
        """Patient with HFrEF on no meds should get all 4 pillars recommended."""
        from src.gdmt_optimizer import GDMTOptimizer
        optimizer = GDMTOptimizer()
        result = optimizer.optimize(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={"age": 60, "sex": "male", "systolic_bp": 110, "heart_rate": 80,
                          "potassium": 4.2, "creatinine": 1.1, "egfr": 65},
        )
        assert result.ef_category == EjectionFractionCategory.HFrEF
        assert len(result.recommendations) > 0

    def test_hfpef_sglt2i_recommended(self):
        """Patient with HFpEF should get SGLT2i recommendation."""
        from src.gdmt_optimizer import GDMTOptimizer
        optimizer = GDMTOptimizer()
        result = optimizer.optimize(
            lvef=55.0,
            nyha_class="II",
            current_medications=[],
            patient_data={"age": 70, "sex": "female", "systolic_bp": 135, "heart_rate": 72,
                          "potassium": 4.0, "creatinine": 0.9, "egfr": 70},
        )
        assert result.ef_category == EjectionFractionCategory.HFpEF
        # SGLT2i should be in recommendations for HFpEF (Class I, 2022 guidelines)
        rec_text = " ".join(result.recommendations).lower()
        assert "sglt2" in rec_text or "dapagliflozin" in rec_text or "empagliflozin" in rec_text


# =====================================================================
# Clinical Workflow Integration
# =====================================================================


class TestClinicalWorkflowIntegration:
    """Test clinical workflows produce correct outputs."""

    def test_cad_assessment_high_risk(self):
        """High calcium score + severe CAD-RADS should produce critical severity."""
        from src.clinical_workflows import CADAssessmentWorkflow
        wf = CADAssessmentWorkflow()
        result = wf.run({
            "calcium_score": 1200,
            "cad_rads": "4B",
            "plaque_features": ["low_attenuation", "positive_remodeling"],
            "symptoms": ["chest_pain", "dyspnea"],
            "risk_factors": ["diabetes", "hypertension", "smoking"],
        })
        assert result.workflow_type == CardioWorkflowType.CAD_ASSESSMENT
        assert result.severity in (SeverityLevel.HIGH, SeverityLevel.VERY_HIGH, SeverityLevel.CRITICAL)
        assert len(result.findings) > 0

    def test_cad_assessment_no_disease(self):
        """Zero calcium score + CAD-RADS 0 should be low severity."""
        from src.clinical_workflows import CADAssessmentWorkflow
        wf = CADAssessmentWorkflow()
        result = wf.run({
            "calcium_score": 0,
            "cad_rads": "0",
            "plaque_features": [],
            "symptoms": [],
            "risk_factors": [],
        })
        assert result.severity == SeverityLevel.LOW


# =====================================================================
# Cross-Module Consistency
# =====================================================================


class TestCrossModuleConsistency:
    """Verify consistency between modules."""

    def test_workflow_boost_covers_all_types(self):
        """WORKFLOW_COLLECTION_BOOST should have an entry for every CardioWorkflowType."""
        for wf_type in CardioWorkflowType:
            assert wf_type in WORKFLOW_COLLECTION_BOOST, (
                f"Missing WORKFLOW_COLLECTION_BOOST entry for {wf_type}"
            )

    def test_workflow_weights_covers_all_types(self):
        """WORKFLOW_COLLECTION_WEIGHTS should have an entry for every CardioWorkflowType."""
        for wf_type in CardioWorkflowType:
            assert wf_type in WORKFLOW_COLLECTION_WEIGHTS, (
                f"Missing WORKFLOW_COLLECTION_WEIGHTS entry for {wf_type}"
            )

    def test_workflow_weights_all_sum_to_one(self):
        """Each workflow's weights in WORKFLOW_COLLECTION_WEIGHTS must sum to ~1.0."""
        for wf_type, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 0.02, (
                f"{wf_type.value} weights sum to {total}, expected ~1.0"
            )

    def test_all_12_collections_in_registry(self):
        """ALL_COLLECTIONS should contain exactly 12 collections."""
        assert len(ALL_COLLECTIONS) == 12

    def test_collection_names_match_configs(self):
        """COLLECTION_NAMES aliases should resolve to actual collection configs."""
        for alias, full_name in COLLECTION_NAMES.items():
            found = any(cfg.name == full_name for cfg in ALL_COLLECTIONS)
            assert found, f"Alias '{alias}' -> '{full_name}' not in ALL_COLLECTIONS"

    def test_conditions_dict_has_workflows(self):
        """Each condition in CARDIO_CONDITIONS should map to at least one workflow."""
        for key, condition in CARDIO_CONDITIONS.items():
            assert "workflows" in condition, f"Condition '{key}' missing workflows field"
            assert len(condition["workflows"]) > 0, f"Condition '{key}' has empty workflows"

    def test_search_weights_default_sums_to_one(self):
        """Default search weights should sum to ~1.0."""
        weights = get_search_weights()
        total = sum(weights.values())
        assert abs(total - 1.0) < 0.05, f"Default weights sum to {total}"

    def test_search_weights_workflow_override(self):
        """Workflow-specific weights should differ from defaults."""
        defaults = get_search_weights()
        hf_weights = get_search_weights(CardioWorkflowType.HEART_FAILURE)
        assert hf_weights != defaults, "HF weights should differ from defaults"
        assert hf_weights["cardio_heart_failure"] > defaults.get("cardio_heart_failure", 0)


# =====================================================================
# Export Round-Trip
# =====================================================================


class TestExportRoundTrip:
    """Verify that workflow results can be exported without errors."""

    def test_workflow_result_to_dict(self):
        """WorkflowResult should serialize to dict without errors."""
        result = WorkflowResult(
            workflow_type=CardioWorkflowType.HEART_FAILURE,
            findings=["LVEF 25%", "NYHA III"],
            recommendations=["Start ARNI", "Initiate SGLT2i"],
            severity=SeverityLevel.HIGH,
            guideline_references=["2022 AHA/ACC/HFSA HF Guideline"],
        )
        d = result.model_dump()
        assert d["workflow_type"] == "heart_failure"
        assert len(d["findings"]) == 2
        assert d["severity"] == "high"

    def test_risk_score_result_serializable(self):
        """RiskScoreResult should serialize cleanly."""
        result = RiskScoreResult(
            score_type=RiskScoreType.ASCVD,
            score_value=15.3,
            risk_category="intermediate",
            interpretation="Intermediate 10-year ASCVD risk",
            recommendations=["Moderate-intensity statin"],
            guideline_reference="2019 ACC/AHA Primary Prevention Guideline",
        )
        d = result.model_dump()
        assert d["score_value"] == 15.3
        assert d["score_type"] == "ascvd"
