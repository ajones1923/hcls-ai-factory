"""Tests for clinical_workflows.py -- Cardiology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Comprehensive tests for all 8 clinical workflow classes, the WorkflowEngine
dispatcher, and helper utilities. Covers calcium score stratification,
CAD-RADS grading, EF classification, NYHA classification, aortic stenosis
severity, arrhythmia detection, cardiac MRI interpretation, Duke Treadmill
Score, preventive risk statin framework, and cardio-oncology GLS/LVEF
monitoring.
"""

import pytest

from src.models import (
    CardioWorkflowType,
    EjectionFractionCategory,
    HeartFailureClass,
    HeartFailureStage,
    SeverityLevel,
    AnticoagulationRecommendation,
    WorkflowResult,
)
from src.clinical_workflows import (
    BaseCardioWorkflow,
    CADAssessmentWorkflow,
    HeartFailureWorkflow,
    ValvularDiseaseWorkflow,
    ArrhythmiaWorkflow,
    CardiacMRIWorkflow,
    StressTestWorkflow,
    PreventiveRiskWorkflow,
    CardioOncologyWorkflow,
    WorkflowEngine,
    _max_severity,
    _trigger_string,
)


# =====================================================================
# WORKFLOW CLASS EXISTENCE & WORKFLOW TYPE
# =====================================================================


class TestWorkflowExistence:
    """Verify all 8 workflow classes exist and have correct workflow_type."""

    def test_cad_assessment_workflow_exists(self):
        wf = CADAssessmentWorkflow()
        assert wf.workflow_type == CardioWorkflowType.CAD_ASSESSMENT

    def test_heart_failure_workflow_exists(self):
        wf = HeartFailureWorkflow()
        assert wf.workflow_type == CardioWorkflowType.HEART_FAILURE

    def test_valvular_disease_workflow_exists(self):
        wf = ValvularDiseaseWorkflow()
        assert wf.workflow_type == CardioWorkflowType.VALVULAR_DISEASE

    def test_arrhythmia_workflow_exists(self):
        wf = ArrhythmiaWorkflow()
        assert wf.workflow_type == CardioWorkflowType.ARRHYTHMIA

    def test_cardiac_mri_workflow_exists(self):
        wf = CardiacMRIWorkflow()
        assert wf.workflow_type == CardioWorkflowType.CARDIAC_MRI

    def test_stress_test_workflow_exists(self):
        wf = StressTestWorkflow()
        assert wf.workflow_type == CardioWorkflowType.STRESS_TEST

    def test_preventive_risk_workflow_exists(self):
        wf = PreventiveRiskWorkflow()
        assert wf.workflow_type == CardioWorkflowType.PREVENTIVE_RISK

    def test_cardio_oncology_workflow_exists(self):
        wf = CardioOncologyWorkflow()
        assert wf.workflow_type == CardioWorkflowType.CARDIO_ONCOLOGY

    def test_all_workflows_are_base_subclass(self):
        for cls in [
            CADAssessmentWorkflow, HeartFailureWorkflow, ValvularDiseaseWorkflow,
            ArrhythmiaWorkflow, CardiacMRIWorkflow, StressTestWorkflow,
            PreventiveRiskWorkflow, CardioOncologyWorkflow,
        ]:
            assert issubclass(cls, BaseCardioWorkflow)

    def test_workflow_types_are_unique(self):
        types = [
            CADAssessmentWorkflow.workflow_type,
            HeartFailureWorkflow.workflow_type,
            ValvularDiseaseWorkflow.workflow_type,
            ArrhythmiaWorkflow.workflow_type,
            CardiacMRIWorkflow.workflow_type,
            StressTestWorkflow.workflow_type,
            PreventiveRiskWorkflow.workflow_type,
            CardioOncologyWorkflow.workflow_type,
        ]
        assert len(types) == len(set(types))


# =====================================================================
# HELPER TESTS
# =====================================================================


class TestHelpers:
    """Test module-level helper functions."""

    def test_max_severity_single(self):
        assert _max_severity(SeverityLevel.LOW) == SeverityLevel.LOW

    def test_max_severity_returns_highest(self):
        result = _max_severity(SeverityLevel.LOW, SeverityLevel.HIGH, SeverityLevel.MODERATE)
        assert result == SeverityLevel.HIGH

    def test_max_severity_with_critical(self):
        result = _max_severity(SeverityLevel.INFORMATIONAL, SeverityLevel.CRITICAL)
        assert result == SeverityLevel.CRITICAL

    def test_trigger_string_format(self):
        result = _trigger_string("test_type", ["GENE1", "GENE2"], "test reason")
        assert "[test_type]" in result
        assert "GENE1" in result
        assert "GENE2" in result
        assert "test reason" in result

    def test_trigger_string_truncates_long_gene_list(self):
        genes = [f"GENE{i}" for i in range(12)]
        result = _trigger_string("type", genes, "reason")
        assert "+4 more" in result


# =====================================================================
# CAD ASSESSMENT WORKFLOW
# =====================================================================


class TestCADAssessmentWorkflow:
    """Test CADAssessmentWorkflow: calcium scoring, CAD-RADS, plaque features."""

    @pytest.fixture
    def wf(self):
        return CADAssessmentWorkflow()

    # -- Calcium score risk stratification --

    def test_calcium_score_zero_very_low_risk(self, wf):
        result = wf.run({"calcium_score": 0})
        assert any("very low" in f.lower() for f in result.findings)

    def test_calcium_score_1_low_risk(self, wf):
        result = wf.run({"calcium_score": 1})
        assert any("low risk" in f.lower() for f in result.findings)

    def test_calcium_score_50_low_risk(self, wf):
        result = wf.run({"calcium_score": 50})
        assert any("low risk" in f.lower() for f in result.findings)

    def test_calcium_score_99_low_risk(self, wf):
        result = wf.run({"calcium_score": 99})
        assert any("low risk" in f.lower() for f in result.findings)

    def test_calcium_score_100_moderate_risk(self, wf):
        result = wf.run({"calcium_score": 100})
        assert any("moderate" in f.lower() for f in result.findings)

    def test_calcium_score_200_moderate_risk(self, wf):
        result = wf.run({"calcium_score": 200})
        assert any("moderate" in f.lower() for f in result.findings)

    def test_calcium_score_399_high_risk(self, wf):
        """CAC 300-999 is high risk per updated tiering."""
        result = wf.run({"calcium_score": 399})
        assert any("high risk" in f.lower() for f in result.findings)

    def test_calcium_score_400_high_risk(self, wf):
        result = wf.run({"calcium_score": 400})
        assert any("high risk" in f.lower() for f in result.findings)

    def test_calcium_score_1000_very_high_risk(self, wf):
        """CAC >=1000 is very high risk per updated tiering."""
        result = wf.run({"calcium_score": 1000})
        assert any("very high" in f.lower() for f in result.findings)

    def test_high_calcium_triggers_statin_rec(self, wf):
        result = wf.run({"calcium_score": 500})
        assert any("statin" in r.lower() for r in result.recommendations)

    def test_high_calcium_triggers_genomic(self, wf):
        result = wf.run({"calcium_score": 500})
        assert len(result.cross_modal_triggers) > 0

    # -- CAD-RADS grading --

    def test_cadrads_0_low_severity(self, wf):
        result = wf.run({"cad_rads": "0", "calcium_score": 0})
        assert any("CAD-RADS 0" in f for f in result.findings)

    def test_cadrads_3_moderate(self, wf):
        result = wf.run({"cad_rads": "3", "calcium_score": 0})
        assert any("50-69%" in f for f in result.findings)

    def test_cadrads_4a_severe(self, wf):
        result = wf.run({"cad_rads": "4A", "calcium_score": 0})
        assert any("70-99%" in f for f in result.findings)

    def test_cadrads_4b_critical(self, wf):
        result = wf.run({"cad_rads": "4B", "calcium_score": 0})
        assert any("left main" in f.lower() for f in result.findings)

    def test_cadrads_5_total_occlusion(self, wf):
        result = wf.run({"cad_rads": "5", "calcium_score": 0})
        assert any("100%" in f for f in result.findings)

    def test_cadrads_normalizes_case(self, wf):
        result = wf.run({"cad_rads": "4a", "calcium_score": 0})
        assert any("CAD-RADS 4A" in f for f in result.findings)

    # -- Plaque features --

    def test_plaque_low_attenuation(self, wf):
        result = wf.run({"plaque_features": ["low_attenuation"], "calcium_score": 0})
        assert any("low-attenuation" in f.lower() for f in result.findings)

    def test_plaque_napkin_ring(self, wf):
        result = wf.run({"plaque_features": ["napkin_ring_sign"], "calcium_score": 0})
        assert any("napkin-ring" in f.lower() for f in result.findings)

    def test_plaque_positive_remodeling(self, wf):
        result = wf.run({"plaque_features": ["positive_remodeling"], "calcium_score": 0})
        assert any("remodelling" in f.lower() or "remodeling" in f.lower() for f in result.findings)

    def test_plaque_spotty_calcification(self, wf):
        result = wf.run({"plaque_features": ["spotty_calcification"], "calcium_score": 0})
        assert any("spotty" in f.lower() for f in result.findings)

    def test_multiple_plaque_features(self, wf):
        features = ["low_attenuation", "napkin_ring_sign", "positive_remodeling"]
        result = wf.run({"plaque_features": features, "calcium_score": 0})
        assert any("3" in f for f in result.findings if "plaque" in f.lower())

    def test_no_plaque_features(self, wf):
        result = wf.run({"plaque_features": [], "calcium_score": 0})
        assert not any("high-risk plaque" in f.lower() for f in result.findings)

    # -- Result structure --

    def test_returns_workflow_result(self, wf):
        result = wf.run({})
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == CardioWorkflowType.CAD_ASSESSMENT

    def test_has_guideline_references(self, wf):
        result = wf.run({})
        assert len(result.guideline_references) > 0

    def test_severity_is_set(self, wf):
        result = wf.run({"calcium_score": 0, "cad_rads": "0"})
        assert result.severity in list(SeverityLevel)


# =====================================================================
# HEART FAILURE WORKFLOW
# =====================================================================


class TestHeartFailureWorkflow:
    """Test HeartFailureWorkflow: EF classification, NYHA, GDMT assessment."""

    @pytest.fixture
    def wf(self):
        return HeartFailureWorkflow()

    # -- EF classification --

    def test_ef_hfref_at_40(self, wf):
        assert wf._classify_ef(40) == EjectionFractionCategory.HFrEF

    def test_ef_hfref_at_30(self, wf):
        assert wf._classify_ef(30) == EjectionFractionCategory.HFrEF

    def test_ef_hfref_at_10(self, wf):
        assert wf._classify_ef(10) == EjectionFractionCategory.HFrEF

    def test_ef_hfmref_at_41(self, wf):
        assert wf._classify_ef(41) == EjectionFractionCategory.HFmrEF

    def test_ef_hfmref_at_49(self, wf):
        assert wf._classify_ef(49) == EjectionFractionCategory.HFmrEF

    def test_ef_hfmref_at_45(self, wf):
        assert wf._classify_ef(45) == EjectionFractionCategory.HFmrEF

    def test_ef_hfpef_at_50(self, wf):
        assert wf._classify_ef(50) == EjectionFractionCategory.HFpEF

    def test_ef_hfpef_at_55(self, wf):
        assert wf._classify_ef(55) == EjectionFractionCategory.HFpEF

    def test_ef_hfpef_at_65(self, wf):
        assert wf._classify_ef(65) == EjectionFractionCategory.HFpEF

    # -- NYHA classification --

    def test_nyha_i(self, wf):
        assert wf._classify_nyha("I") == HeartFailureClass.NYHA_I

    def test_nyha_ii(self, wf):
        assert wf._classify_nyha("II") == HeartFailureClass.NYHA_II

    def test_nyha_iii(self, wf):
        assert wf._classify_nyha("III") == HeartFailureClass.NYHA_III

    def test_nyha_iv(self, wf):
        assert wf._classify_nyha("IV") == HeartFailureClass.NYHA_IV

    def test_nyha_numeric_1(self, wf):
        assert wf._classify_nyha("1") == HeartFailureClass.NYHA_I

    def test_nyha_numeric_2(self, wf):
        assert wf._classify_nyha("2") == HeartFailureClass.NYHA_II

    def test_nyha_numeric_3(self, wf):
        assert wf._classify_nyha("3") == HeartFailureClass.NYHA_III

    def test_nyha_numeric_4(self, wf):
        assert wf._classify_nyha("4") == HeartFailureClass.NYHA_IV

    def test_nyha_default_for_unknown(self, wf):
        assert wf._classify_nyha("unknown") == HeartFailureClass.NYHA_II

    # -- HF stage classification --

    def test_stage_d_high_bnp(self, wf):
        stage = wf._classify_stage(
            EjectionFractionCategory.HFrEF, HeartFailureClass.NYHA_IV,
            bnp=700, nt_probnp=0,
        )
        assert stage == HeartFailureStage.STAGE_D

    def test_stage_d_high_nt_probnp(self, wf):
        stage = wf._classify_stage(
            EjectionFractionCategory.HFrEF, HeartFailureClass.NYHA_III,
            bnp=0, nt_probnp=6000,
        )
        assert stage == HeartFailureStage.STAGE_D

    def test_stage_c_symptomatic(self, wf):
        stage = wf._classify_stage(
            EjectionFractionCategory.HFrEF, HeartFailureClass.NYHA_III,
            bnp=50, nt_probnp=200,
        )
        assert stage == HeartFailureStage.STAGE_C

    def test_stage_b_structural(self, wf):
        stage = wf._classify_stage(
            EjectionFractionCategory.HFrEF, HeartFailureClass.NYHA_I,
            bnp=50, nt_probnp=0,
        )
        assert stage == HeartFailureStage.STAGE_B

    def test_stage_a_no_structural(self, wf):
        stage = wf._classify_stage(
            EjectionFractionCategory.HFpEF, HeartFailureClass.NYHA_I,
            bnp=50, nt_probnp=100,
        )
        assert stage == HeartFailureStage.STAGE_A

    # -- Run workflow --

    def test_hfref_run(self, wf):
        result = wf.run({"lvef": 30, "nyha_class": "III"})
        assert result.workflow_type == CardioWorkflowType.HEART_FAILURE
        assert any("HFrEF" in f for f in result.findings)

    def test_gdmt_gap_detection(self, wf):
        result = wf.run({"lvef": 30, "nyha_class": "II", "current_meds": []})
        assert any("GDMT gaps" in r for r in result.recommendations)

    def test_all_gdmt_present(self, wf):
        meds = ["carvedilol", "sacubitril/valsartan", "spironolactone", "dapagliflozin"]
        result = wf.run({"lvef": 30, "nyha_class": "II", "current_meds": meds})
        assert any("four gdmt pillars" in f.lower() for f in result.findings)

    def test_device_eligibility_icd(self, wf):
        result = wf.run({"lvef": 30, "nyha_class": "II"})
        assert any("ICD" in r for r in result.recommendations)

    def test_crt_recommendation_lbbb(self, wf):
        result = wf.run({
            "lvef": 25, "nyha_class": "II",
            "qrs_duration": 160, "qrs_morphology": "lbbb",
        })
        assert any("CRT" in r for r in result.recommendations)

    def test_young_age_triggers_genetic(self, wf):
        result = wf.run({"lvef": 30, "age": 35})
        assert len(result.cross_modal_triggers) > 0

    def test_hfpef_sglt2i_rec(self, wf):
        result = wf.run({"lvef": 60, "nyha_class": "II"})
        assert any("SGLT2" in r for r in result.recommendations)

    def test_has_guideline_references(self, wf):
        result = wf.run({})
        assert len(result.guideline_references) >= 2


# =====================================================================
# VALVULAR DISEASE WORKFLOW
# =====================================================================


class TestValvularDiseaseWorkflow:
    """Test ValvularDiseaseWorkflow: AS severity, MR grading, intervention criteria."""

    @pytest.fixture
    def wf(self):
        return ValvularDiseaseWorkflow()

    # -- AS severity --

    def test_as_mild(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 2.5, "mean_gradient": 15, "valve_area": 1.8},
        })
        assert any("mild" in f.lower() for f in result.findings)

    def test_as_moderate(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 3.5, "mean_gradient": 30, "valve_area": 1.2},
        })
        assert any("moderate" in f.lower() for f in result.findings if "severity" in f.lower())

    def test_as_severe(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
        })
        assert any("severe" in f.lower() for f in result.findings if "severity" in f.lower())

    def test_as_severe_symptomatic_avr(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
            "symptoms": ["dyspnea", "chest pain"],
        })
        assert any("AVR" in r for r in result.recommendations)

    def test_as_severe_low_ef_avr(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
            "lvef": 40,
        })
        assert any("AVR" in r and "LVEF" in r for r in result.recommendations)

    def test_sts_low_risk(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
            "symptoms": ["dyspnea"], "sts_score": 1.5,
        })
        assert any("low risk" in r.lower() for r in result.recommendations)

    def test_sts_high_risk_tavr(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
            "symptoms": ["dyspnea"], "sts_score": 10.0,
        })
        assert any("TAVR" in r for r in result.recommendations)

    def test_dimensionless_index(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45,
                             "valve_area": 0.8, "dimensionless_index": 0.2},
        })
        assert any("DI" in f for f in result.findings)

    # -- MR severity --

    def test_mr_mild(self, wf):
        result = wf.run({
            "valve": "mitral", "pathology": "regurgitation",
            "measurements": {"regurgitant_volume": 20, "ero": 0.1, "vena_contracta": 0.2},
        })
        assert any("mild" in f.lower() for f in result.findings)

    def test_mr_severe(self, wf):
        result = wf.run({
            "valve": "mitral", "pathology": "regurgitation",
            "measurements": {"regurgitant_volume": 70, "ero": 0.5, "vena_contracta": 0.8},
        })
        assert any("severe" in f.lower() for f in result.findings if "severity" in f.lower())

    def test_mr_severe_symptomatic_surgery(self, wf):
        result = wf.run({
            "valve": "mitral", "pathology": "regurgitation",
            "measurements": {"regurgitant_volume": 70, "ero": 0.5, "vena_contracta": 0.8},
            "symptoms": ["dyspnea"], "lvef": 55,
        })
        assert any("surgery" in r.lower() for r in result.recommendations)

    # -- Generic valve --

    def test_generic_valve_defers(self, wf):
        result = wf.run({
            "valve": "tricuspid", "pathology": "regurgitation",
            "measurements": {"regurgitant_volume": 40},
        })
        assert any("specialist" in f.lower() for f in result.findings)

    # -- Severity mapping --

    def test_mild_maps_to_low_severity(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 2.5, "mean_gradient": 15, "valve_area": 1.8},
        })
        assert result.severity == SeverityLevel.LOW

    def test_severe_maps_to_high_severity(self, wf):
        result = wf.run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.5, "mean_gradient": 45, "valve_area": 0.8},
        })
        assert result.severity == SeverityLevel.HIGH


# =====================================================================
# ARRHYTHMIA WORKFLOW
# =====================================================================


class TestArrhythmiaWorkflow:
    """Test ArrhythmiaWorkflow: AF detection, anticoag, QTc, bradycardia."""

    @pytest.fixture
    def wf(self):
        return ArrhythmiaWorkflow()

    # -- AF detection --

    def test_af_detection_atrial_fibrillation(self, wf):
        result = wf.run({"rhythm": "atrial fibrillation", "age": 70, "sex": "male"})
        assert any("Atrial fibrillation" in f for f in result.findings)

    def test_af_detection_afib(self, wf):
        result = wf.run({"rhythm": "afib"})
        assert any("fibrillation" in f.lower() for f in result.findings)

    def test_af_detection_af(self, wf):
        result = wf.run({"rhythm": "af"})
        assert any("fibrillation" in f.lower() for f in result.findings)

    # -- CHA2DS2-VASc --

    def test_cha2ds2vasc_zero(self, wf):
        score, _ = wf._cha2ds2_vasc(50, "male", {})
        assert score == 0

    def test_cha2ds2vasc_chf(self, wf):
        score, _ = wf._cha2ds2_vasc(50, "male", {"chf": True})
        assert score == 1

    def test_cha2ds2vasc_age_75(self, wf):
        score, breakdown = wf._cha2ds2_vasc(75, "male", {})
        assert score == 2

    def test_cha2ds2vasc_age_65(self, wf):
        score, _ = wf._cha2ds2_vasc(65, "male", {})
        assert score == 1

    def test_cha2ds2vasc_female(self, wf):
        score, _ = wf._cha2ds2_vasc(50, "female", {})
        assert score == 1

    def test_cha2ds2vasc_stroke(self, wf):
        score, _ = wf._cha2ds2_vasc(50, "male", {"stroke_tia": True})
        assert score == 2

    def test_cha2ds2vasc_max(self, wf):
        score, _ = wf._cha2ds2_vasc(80, "female", {
            "chf": True, "hypertension": True, "diabetes": True,
            "stroke_tia": True, "vascular_disease": True,
        })
        assert score == 9

    # -- Anticoag recommendation --

    def test_anticoag_no_anticoag_male_0(self, wf):
        rec, _ = wf._anticoag_recommendation(0, "male")
        assert rec == AnticoagulationRecommendation.NO_ANTICOAG

    def test_anticoag_no_anticoag_female_1(self, wf):
        rec, _ = wf._anticoag_recommendation(1, "female")
        assert rec == AnticoagulationRecommendation.NO_ANTICOAG

    def test_anticoag_consider_male_1(self, wf):
        rec, _ = wf._anticoag_recommendation(1, "male")
        assert rec == AnticoagulationRecommendation.CONSIDER_ANTICOAG

    def test_anticoag_recommended_score_2(self, wf):
        rec, _ = wf._anticoag_recommendation(2, "male")
        assert rec == AnticoagulationRecommendation.ANTICOAG_RECOMMENDED

    def test_anticoag_recommended_high_score(self, wf):
        rec, text = wf._anticoag_recommendation(5, "male")
        assert rec == AnticoagulationRecommendation.ANTICOAG_RECOMMENDED
        assert "DOAC" in text

    # -- QTc assessment --

    def test_qtc_normal_male(self, wf):
        result = wf.run({"qtc": 440, "sex": "male"})
        assert any("normal" in f.lower() for f in result.findings if "QTc" in f)

    def test_qtc_prolonged(self, wf):
        result = wf.run({"qtc": 480, "sex": "male"})
        assert any("borderline" in f.lower() for f in result.findings if "QTc" in f)

    def test_qtc_significantly_prolonged(self, wf):
        result = wf.run({"qtc": 520, "sex": "male"})
        assert any("significantly prolonged" in f.lower() for f in result.findings if "QTc" in f)

    def test_qtc_high_triggers_genetic(self, wf):
        result = wf.run({"qtc": 520})
        assert len(result.cross_modal_triggers) > 0

    # -- Bradycardia --

    def test_bradycardia_detected(self, wf):
        result = wf.run({"heart_rate": 40})
        assert any("Bradycardia" in f for f in result.findings)

    def test_high_grade_block_pacemaker(self, wf):
        result = wf.run({
            "heart_rate": 40,
            "findings": ["complete heart block"],
        })
        assert any("pacemaker" in r.lower() for r in result.recommendations)

    # -- Wide complex tachycardia --

    def test_wide_complex_tachycardia(self, wf):
        result = wf.run({"qrs_duration": 140, "heart_rate": 150})
        assert any("wide-complex" in f.lower() for f in result.findings)

    # -- Severity --

    def test_severity_af_moderate(self, wf):
        result = wf.run({"rhythm": "atrial fibrillation"})
        assert result.severity == SeverityLevel.MODERATE

    def test_severity_high_qtc(self, wf):
        result = wf.run({"qtc": 520})
        assert result.severity == SeverityLevel.HIGH


# =====================================================================
# CARDIAC MRI WORKFLOW
# =====================================================================


class TestCardiacMRIWorkflow:
    """Test CardiacMRIWorkflow: LGE patterns, T1/T2 mapping, ECV."""

    @pytest.fixture
    def wf(self):
        return CardiacMRIWorkflow()

    # -- LGE patterns --

    def test_lge_subendocardial(self, wf):
        result = wf.run({"lge_pattern": "subendocardial", "lge_extent_percent": 10})
        assert any("ischaemic" in f.lower() for f in result.findings)

    def test_lge_mid_wall(self, wf):
        result = wf.run({"lge_pattern": "mid_wall", "lge_extent_percent": 10})
        assert any("non-ischaemic" in f.lower() for f in result.findings)
        assert len(result.cross_modal_triggers) > 0

    def test_lge_epicardial(self, wf):
        result = wf.run({"lge_pattern": "epicardial", "lge_extent_percent": 5})
        assert any("myocarditis" in f.lower() for f in result.findings)

    def test_lge_rv_insertion(self, wf):
        result = wf.run({"lge_pattern": "rv_insertion", "lge_extent_percent": 5})
        assert any("hcm" in f.lower() or "hypertrophic" in f.lower() for f in result.findings)

    def test_lge_diffuse_amyloid(self, wf):
        result = wf.run({"lge_pattern": "diffuse", "lge_extent_percent": 20})
        assert any("amyloid" in f.lower() for f in result.findings)

    def test_lge_none(self, wf):
        result = wf.run({"lge_pattern": "none"})
        assert any("No late gadolinium" in f for f in result.findings)

    def test_extensive_lge_icd_rec(self, wf):
        result = wf.run({"lge_pattern": "mid_wall", "lge_extent_percent": 20, "lvef": 30})
        assert any("ICD" in r for r in result.recommendations)

    # -- T1 mapping --

    def test_t1_elevated(self, wf):
        result = wf.run({"t1_native": 1150})
        assert any("elevated" in f.lower() for f in result.findings if "T1" in f)

    def test_t1_low(self, wf):
        result = wf.run({"t1_native": 850})
        assert any("low" in f.lower() for f in result.findings if "T1" in f)

    def test_t1_normal(self, wf):
        result = wf.run({"t1_native": 1000})
        assert any("normal" in f.lower() for f in result.findings if "T1" in f)

    # -- T2 mapping --

    def test_t2_elevated(self, wf):
        result = wf.run({"t2_value": 60})
        assert any("elevated" in f.lower() for f in result.findings if "T2" in f)

    def test_t2_normal(self, wf):
        result = wf.run({"t2_value": 45})
        assert any("normal" in f.lower() for f in result.findings if "T2" in f)

    # -- ECV --

    def test_ecv_elevated(self, wf):
        result = wf.run({"ecv_percent": 35})
        assert any("elevated" in f.lower() for f in result.findings if "ECV" in f)

    def test_ecv_markedly_elevated(self, wf):
        result = wf.run({"ecv_percent": 45})
        assert any("markedly" in f.lower() for f in result.findings if "ECV" in f)

    def test_ecv_normal(self, wf):
        result = wf.run({"ecv_percent": 27})
        assert any("normal" in f.lower() for f in result.findings if "ECV" in f)

    # -- Severity --

    def test_severity_high_low_ef(self, wf):
        result = wf.run({"lvef": 30, "lge_pattern": "none"})
        assert result.severity == SeverityLevel.HIGH

    def test_severity_informational_normal(self, wf):
        result = wf.run({"lvef": 60, "lge_pattern": "none", "lge_extent_percent": 0})
        assert result.severity == SeverityLevel.INFORMATIONAL


# =====================================================================
# STRESS TEST WORKFLOW
# =====================================================================


class TestStressTestWorkflow:
    """Test StressTestWorkflow: DTS calculation, perfusion, HR recovery."""

    @pytest.fixture
    def wf(self):
        return StressTestWorkflow()

    # -- Duke Treadmill Score --

    def test_dts_low_risk(self, wf):
        # DTS = 12 - 5*0 - 4*0 = 12 (low risk)
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 12.0,
            "max_st_deviation_mm": 0.0, "angina_index": 0,
        })
        assert any("Low risk" in f for f in result.findings if "Duke" in f)

    def test_dts_moderate_risk(self, wf):
        # DTS = 9 - 5*2 - 4*0 = -1 (moderate)
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 9.0,
            "max_st_deviation_mm": 2.0, "angina_index": 0,
        })
        assert any("Moderate risk" in f for f in result.findings if "Duke" in f)

    def test_dts_high_risk(self, wf):
        # DTS = 3 - 5*3 - 4*2 = 3 - 15 - 8 = -20 (high risk)
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 3.0,
            "max_st_deviation_mm": 3.0, "angina_index": 2,
        })
        assert any("High risk" in f for f in result.findings if "Duke" in f)

    def test_dts_formula_correct(self, wf):
        # DTS = exercise_time - 5 * ST_dev - 4 * angina
        # DTS = 10 - 5*1 - 4*1 = 10 - 5 - 4 = 1
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 10.0,
            "max_st_deviation_mm": 1.0, "angina_index": 1,
        })
        assert any("1.0" in f for f in result.findings if "Duke" in f)

    def test_dts_boundary_low(self, wf):
        # DTS = 5 => Low risk boundary
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 5.0,
            "max_st_deviation_mm": 0.0, "angina_index": 0,
        })
        assert any("Low risk" in f for f in result.findings if "Duke" in f)

    def test_dts_boundary_moderate(self, wf):
        # DTS = -10 => Moderate risk boundary
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 0.0,
            "max_st_deviation_mm": 2.0, "angina_index": 0,
        })
        assert any("Moderate" in f or "High" in f for f in result.findings if "Duke" in f)

    # -- Perfusion --

    def test_reversible_defect(self, wf):
        result = wf.run({"perfusion_defect": "reversible"})
        assert any("reversible" in f.lower() for f in result.findings if "perfusion" in f.lower())

    def test_fixed_defect(self, wf):
        result = wf.run({"perfusion_defect": "fixed"})
        assert any("fixed" in f.lower() and "scar" in f.lower() for f in result.findings)

    def test_normal_perfusion(self, wf):
        result = wf.run({"perfusion_defect": "none"})
        assert any("normal" in f.lower() for f in result.findings if "perfusion" in f.lower())

    # -- HR recovery --

    def test_abnormal_hr_recovery(self, wf):
        result = wf.run({"hr_recovery_1min": 8})
        assert any("abnormal" in f.lower() for f in result.findings if "recovery" in f.lower())

    def test_normal_hr_recovery(self, wf):
        result = wf.run({"hr_recovery_1min": 20})
        assert any("normal" in f.lower() for f in result.findings if "recovery" in f.lower())

    # -- BP response --

    def test_hypotensive_response(self, wf):
        result = wf.run({"bp_response": "hypotensive"})
        assert any("hypotensive" in f.lower() for f in result.findings)

    # -- Severity --

    def test_severity_high_st(self, wf):
        result = wf.run({
            "test_type": "exercise", "exercise_time_min": 5,
            "max_st_deviation_mm": 2.5, "angina_index": 0,
        })
        assert result.severity in (SeverityLevel.HIGH, SeverityLevel.CRITICAL)


# =====================================================================
# PREVENTIVE RISK WORKFLOW
# =====================================================================


class TestPreventiveRiskWorkflow:
    """Test PreventiveRiskWorkflow: statin decision framework, ASCVD risk."""

    @pytest.fixture
    def wf(self):
        return PreventiveRiskWorkflow()

    # -- Statin framework --

    def test_ldl_190_statin(self, wf):
        result = wf.run({"ldl": 200})
        assert any("high-intensity statin" in r.lower() for r in result.recommendations)

    def test_ldl_190_triggers_fh(self, wf):
        result = wf.run({"ldl": 200})
        assert len(result.cross_modal_triggers) > 0

    def test_diabetes_moderate_statin(self, wf):
        result = wf.run({"diabetes": True, "age": 55, "ldl": 130})
        assert any("statin" in r.lower() for r in result.recommendations)

    def test_low_risk_lifestyle(self, wf):
        result = wf.run({
            "age": 45, "sex": "male", "race": "white",
            "total_cholesterol": 180, "hdl": 55, "ldl": 100,
            "systolic_bp": 120, "diabetes": False, "smoker": False,
        })
        assert any("lifestyle" in r.lower() for r in result.recommendations)

    # -- Risk enhancers --

    def test_risk_enhancer_lpa(self, wf):
        result = wf.run({"lpa": 60})
        assert any("Lp(a)" in f for f in result.findings)

    def test_risk_enhancer_hscrp(self, wf):
        result = wf.run({"hscrp": 3.0})
        assert any("CRP" in f for f in result.findings)

    def test_risk_enhancer_low_abi(self, wf):
        result = wf.run({"abi": 0.8})
        assert any("ABI" in f for f in result.findings)

    # -- ASCVD outside range --

    def test_age_outside_range(self, wf):
        result = wf.run({"age": 30})
        assert any("outside" in f.lower() for f in result.findings)

    # -- Severity --

    def test_severity_high_ldl_190(self, wf):
        result = wf.run({"ldl": 200})
        assert result.severity == SeverityLevel.HIGH


# =====================================================================
# CARDIO-ONCOLOGY WORKFLOW
# =====================================================================


class TestCardioOncologyWorkflow:
    """Test CardioOncologyWorkflow: GLS, LVEF, cardiotoxicity detection."""

    @pytest.fixture
    def wf(self):
        return CardioOncologyWorkflow()

    # -- GLS decline --

    def test_gls_subclinical_over_15(self, wf):
        # baseline -20, current -16 => 20% decline (>15%)
        result = wf.run({
            "baseline_gls": -20.0, "current_gls": -16.0,
        })
        assert any("SUBCLINICAL" in f for f in result.findings)

    def test_gls_borderline(self, wf):
        # baseline -20, current -18 => 10% decline (8-15%)
        result = wf.run({
            "baseline_gls": -20.0, "current_gls": -18.0,
        })
        assert any("borderline" in f.lower() for f in result.findings)

    def test_gls_acceptable(self, wf):
        result = wf.run({
            "baseline_gls": -20.0, "current_gls": -19.5,
        })
        assert any("acceptable" in f.lower() for f in result.findings)

    # -- LVEF drop --

    def test_lvef_ctrcd_confirmed(self, wf):
        result = wf.run({
            "baseline_lvef": 60, "current_lvef": 45,
        })
        assert any("CTRCD" in f for f in result.findings)

    def test_lvef_significant_decline_above_50(self, wf):
        result = wf.run({
            "baseline_lvef": 65, "current_lvef": 52,
        })
        assert any("significant" in f.lower() for f in result.findings if "LVEF" in f)

    def test_lvef_acceptable(self, wf):
        result = wf.run({
            "baseline_lvef": 60, "current_lvef": 58,
        })
        assert any("acceptable" in f.lower() for f in result.findings if "LVEF" in f)

    # -- Agent protocols --

    def test_anthracycline_protocol(self, wf):
        result = wf.run({"chemotherapy_agent": "anthracycline"})
        assert any("anthracycline" in f.lower() for f in result.findings)

    def test_trastuzumab_protocol(self, wf):
        result = wf.run({"chemotherapy_agent": "trastuzumab"})
        assert any("trastuzumab" in f.lower() for f in result.findings)

    def test_ici_protocol(self, wf):
        result = wf.run({"chemotherapy_agent": "immune_checkpoint_inhibitor"})
        assert any("immune" in f.lower() or "checkpoint" in f.lower() for f in result.findings)

    def test_unknown_agent_defaults(self, wf):
        result = wf.run({"chemotherapy_agent": "unknown_drug"})
        assert any("not in database" in f.lower() for f in result.findings)

    # -- Cumulative dose --

    def test_cumulative_dose_exceeds(self, wf):
        result = wf.run({
            "chemotherapy_agent": "anthracycline",
            "cumulative_dose": 450,
        })
        assert any("exceeds" in f.lower() for f in result.findings)

    # -- Biomarkers --

    def test_elevated_troponin(self, wf):
        result = wf.run({"troponin": 0.05})
        assert any("elevated" in f.lower() for f in result.findings if "Troponin" in f)

    def test_elevated_nt_probnp(self, wf):
        result = wf.run({"nt_probnp": 500})
        assert any("elevated" in f.lower() for f in result.findings if "NT-proBNP" in f)

    # -- Severity --

    def test_severity_critical_ctrcd(self, wf):
        result = wf.run({
            "baseline_lvef": 60, "current_lvef": 40,
            "chemotherapy_agent": "anthracycline", "cumulative_dose": 450,
        })
        assert result.severity == SeverityLevel.CRITICAL


# =====================================================================
# WORKFLOW ENGINE
# =====================================================================


class TestWorkflowEngine:
    """Test WorkflowEngine dispatcher and workflow detection."""

    @pytest.fixture
    def engine(self):
        return WorkflowEngine()

    def test_engine_creation(self):
        engine = WorkflowEngine()
        assert engine is not None

    def test_run_workflow_cad(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.CAD_ASSESSMENT,
            {"calcium_score": 100},
        )
        assert result.workflow_type == CardioWorkflowType.CAD_ASSESSMENT

    def test_run_workflow_hf(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.HEART_FAILURE,
            {"lvef": 30},
        )
        assert result.workflow_type == CardioWorkflowType.HEART_FAILURE

    def test_run_workflow_valvular(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.VALVULAR_DISEASE,
            {"valve": "aortic", "pathology": "stenosis"},
        )
        assert result.workflow_type == CardioWorkflowType.VALVULAR_DISEASE

    def test_run_workflow_arrhythmia(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.ARRHYTHMIA,
            {"rhythm": "sinus"},
        )
        assert result.workflow_type == CardioWorkflowType.ARRHYTHMIA

    def test_run_workflow_mri(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.CARDIAC_MRI,
            {"lge_pattern": "none"},
        )
        assert result.workflow_type == CardioWorkflowType.CARDIAC_MRI

    def test_run_workflow_stress(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.STRESS_TEST,
            {"test_type": "exercise"},
        )
        assert result.workflow_type == CardioWorkflowType.STRESS_TEST

    def test_run_workflow_preventive(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.PREVENTIVE_RISK,
            {"age": 55},
        )
        assert result.workflow_type == CardioWorkflowType.PREVENTIVE_RISK

    def test_run_workflow_cardio_oncology(self, engine):
        result = engine.run_workflow(
            CardioWorkflowType.CARDIO_ONCOLOGY,
            {"chemotherapy_agent": "anthracycline"},
        )
        assert result.workflow_type == CardioWorkflowType.CARDIO_ONCOLOGY

    # -- detect_workflow --

    def test_detect_workflow_cad(self, engine):
        wf = engine.detect_workflow("patient with coronary artery disease")
        assert wf == CardioWorkflowType.CAD_ASSESSMENT

    def test_detect_workflow_hf(self, engine):
        wf = engine.detect_workflow("heart failure with reduced EF")
        assert wf == CardioWorkflowType.HEART_FAILURE

    def test_detect_workflow_valve(self, engine):
        wf = engine.detect_workflow("severe aortic stenosis evaluation")
        assert wf == CardioWorkflowType.VALVULAR_DISEASE

    def test_detect_workflow_arrhythmia(self, engine):
        wf = engine.detect_workflow("atrial fibrillation management")
        assert wf == CardioWorkflowType.ARRHYTHMIA

    def test_detect_workflow_none(self, engine):
        wf = engine.detect_workflow("random unrelated text")
        assert wf is None


# =====================================================================
# SAFETY-CRITICAL REGRESSION TESTS
# =====================================================================


class TestSafetyCriticalClinicalFixes:
    """Regression tests for clinically significant workflow fixes.

    Each test documents the clinical rationale and source guideline.
    """

    # -- HR recovery boundary: <=12 bpm is abnormal (Cole et al. NEJM 1999) --

    def test_hr_recovery_12_is_abnormal(self):
        """HR recovery of exactly 12 bpm at 1 min is abnormal (<=12 threshold).

        Source: Cole CR et al. NEJM 1999;341:1351-7.
        """
        wf = StressTestWorkflow()
        result = wf.run({"hr_recovery_1min": 12})
        assert any("abnormal" in f.lower() for f in result.findings if "recovery" in f.lower())

    def test_hr_recovery_13_is_normal(self):
        """HR recovery of 13 bpm at 1 min is normal (>12 threshold)."""
        wf = StressTestWorkflow()
        result = wf.run({"hr_recovery_1min": 13})
        assert any("normal" in f.lower() for f in result.findings if "recovery" in f.lower())

    # -- QTc: male upper limit 440ms, female 460ms --

    def test_qtc_male_440_normal(self):
        """Male QTc 440ms should be normal (upper limit is 440ms)."""
        wf = ArrhythmiaWorkflow()
        result = wf.run({"qtc": 440, "sex": "male"})
        prolonged = [f for f in result.findings if "prolong" in f.lower()]
        assert len(prolonged) == 0

    def test_qtc_male_441_prolonged(self):
        """Male QTc 441ms should be prolonged (>440ms threshold)."""
        wf = ArrhythmiaWorkflow()
        result = wf.run({"qtc": 441, "sex": "male"})
        assert any("prolong" in f.lower() for f in result.findings)

    # -- CAC scoring tiers --

    def test_cac_299_is_moderate(self):
        """CAC 299 should be moderate risk (100-299 range)."""
        wf = CADAssessmentWorkflow()
        result = wf.run({"calcium_score": 299})
        assert any("moderate" in f.lower() for f in result.findings)

    def test_cac_300_is_high(self):
        """CAC 300 should be high risk (300-999 range)."""
        wf = CADAssessmentWorkflow()
        result = wf.run({"calcium_score": 300})
        assert any("high risk" in f.lower() for f in result.findings)

    def test_cac_999_is_high(self):
        """CAC 999 should be high risk (300-999 range)."""
        wf = CADAssessmentWorkflow()
        result = wf.run({"calcium_score": 999})
        assert any("high risk" in f.lower() for f in result.findings)

    def test_cac_1000_is_very_high(self):
        """CAC >=1000 should be very high risk."""
        wf = CADAssessmentWorkflow()
        result = wf.run({"calcium_score": 1000})
        assert any("very high" in f.lower() for f in result.findings)

    # -- AF detection (flexible matching) --

    def test_detect_af_in_query(self):
        """'atrial fibrillation' should be detected as arrhythmia workflow."""
        engine = WorkflowEngine()
        wf = engine.detect_workflow("patient with atrial fibrillation")
        assert wf == CardioWorkflowType.ARRHYTHMIA

    def test_detect_vt_in_query(self):
        """'ventricular tachycardia' should be detected as arrhythmia workflow."""
        engine = WorkflowEngine()
        wf = engine.detect_workflow("ventricular tachycardia management")
        assert wf == CardioWorkflowType.ARRHYTHMIA

    def test_detect_brugada_in_query(self):
        """'brugada' should be detected as arrhythmia workflow."""
        engine = WorkflowEngine()
        wf = engine.detect_workflow("brugada syndrome evaluation")
        assert wf == CardioWorkflowType.ARRHYTHMIA

    def test_detect_wpw_in_query(self):
        """'wpw' should be detected as arrhythmia workflow."""
        engine = WorkflowEngine()
        wf = engine.detect_workflow("wpw pre-excitation pattern")
        assert wf == CardioWorkflowType.ARRHYTHMIA
