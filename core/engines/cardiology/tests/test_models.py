"""Unit tests for Cardiology Intelligence Agent models.

Tests all enums (values, members, string conversion) and Pydantic models
(creation, defaults, validation, field types) defined in src/models.py.

Author: Adam Jones
Date: March 2026
"""

import pytest
from pydantic import ValidationError

from src.models import (
    AnticoagulationRecommendation,
    CADRADSScore,
    CardioQuery,
    CardioResponse,
    CardioSearchResult,
    CardiotoxicityAssessment,
    CardiotoxicityRisk,
    CardioWorkflowType,
    CrossModalTrigger,
    ECGInterpretation,
    EjectionFractionCategory,
    EvidenceLevel,
    GDMTMedication,
    GDMTPillar,
    GDMTRecommendation,
    GDMTStatus,
    GuidelineClass,
    HeartFailureClass,
    HeartFailureStage,
    ImagingModality,
    ImagingResult,
    LGEPattern,
    RiskScoreInput,
    RiskScoreResult,
    RiskScoreType,
    SearchPlan,
    SeverityLevel,
    ValveAssessment,
    ValveSeverity,
    WorkflowResult,
)


# ===================================================================
# ENUM TESTS: CardioWorkflowType (8 members)
# ===================================================================


class TestCardioWorkflowType:
    """Tests for CardioWorkflowType enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("CAD_ASSESSMENT", "cad_assessment"),
            ("HEART_FAILURE", "heart_failure"),
            ("VALVULAR_DISEASE", "valvular_disease"),
            ("ARRHYTHMIA", "arrhythmia"),
            ("CARDIAC_MRI", "cardiac_mri"),
            ("STRESS_TEST", "stress_test"),
            ("PREVENTIVE_RISK", "preventive_risk"),
            ("CARDIO_ONCOLOGY", "cardio_oncology"),
            ("ACUTE_DECOMPENSATED_HF", "acute_decompensated_hf"),
            ("POST_MI", "post_mi"),
            ("MYOCARDITIS_PERICARDITIS", "myocarditis_pericarditis"),
        ],
    )
    def test_member_value(self, member, value):
        assert CardioWorkflowType[member].value == value

    def test_member_count(self):
        assert len(CardioWorkflowType) == 12

    @pytest.mark.parametrize(
        "member",
        [
            "CAD_ASSESSMENT", "HEART_FAILURE", "VALVULAR_DISEASE",
            "ARRHYTHMIA", "CARDIAC_MRI", "STRESS_TEST",
            "PREVENTIVE_RISK", "CARDIO_ONCOLOGY",
            "ACUTE_DECOMPENSATED_HF", "POST_MI", "MYOCARDITIS_PERICARDITIS",
        ],
    )
    def test_member_exists(self, member):
        assert hasattr(CardioWorkflowType, member)

    def test_string_conversion(self):
        assert str(CardioWorkflowType.HEART_FAILURE) == "CardioWorkflowType.HEART_FAILURE"

    def test_value_access(self):
        wf = CardioWorkflowType("heart_failure")
        assert wf == CardioWorkflowType.HEART_FAILURE

    def test_is_str_subclass(self):
        assert isinstance(CardioWorkflowType.CAD_ASSESSMENT, str)


# ===================================================================
# ENUM TESTS: RiskScoreType (6 members)
# ===================================================================


class TestRiskScoreType:
    """Tests for RiskScoreType enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("ASCVD", "ascvd"),
            ("HEART", "heart"),
            ("CHA2DS2_VASC", "cha2ds2_vasc"),
            ("HAS_BLED", "has_bled"),
            ("MAGGIC", "maggic"),
            ("EUROSCORE_II", "euroscore_ii"),
        ],
    )
    def test_member_value(self, member, value):
        assert RiskScoreType[member].value == value

    def test_member_count(self):
        assert len(RiskScoreType) == 6

    def test_value_lookup(self):
        assert RiskScoreType("ascvd") == RiskScoreType.ASCVD

    def test_is_str_subclass(self):
        assert isinstance(RiskScoreType.HEART, str)


# ===================================================================
# ENUM TESTS: SeverityLevel (5 members)
# ===================================================================


class TestSeverityLevel:
    """Tests for SeverityLevel enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("CRITICAL", "critical"),
            ("HIGH", "high"),
            ("MODERATE", "moderate"),
            ("LOW", "low"),
            ("INFORMATIONAL", "informational"),
        ],
    )
    def test_member_value(self, member, value):
        assert SeverityLevel[member].value == value

    def test_member_count(self):
        assert len(SeverityLevel) == 8


# ===================================================================
# ENUM TESTS: ImagingModality (5 members)
# ===================================================================


class TestImagingModality:
    """Tests for ImagingModality enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("ECHOCARDIOGRAPHY", "echocardiography"),
            ("CARDIAC_CT", "cardiac_ct"),
            ("CARDIAC_MRI", "cardiac_mri"),
            ("NUCLEAR_CARDIOLOGY", "nuclear_cardiology"),
            ("CORONARY_ANGIOGRAPHY", "coronary_angiography"),
        ],
    )
    def test_member_value(self, member, value):
        assert ImagingModality[member].value == value

    def test_member_count(self):
        assert len(ImagingModality) == 5


# ===================================================================
# ENUM TESTS: HeartFailureClass (4 members)
# ===================================================================


class TestHeartFailureClass:
    """Tests for HeartFailureClass enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("NYHA_I", "nyha_i"),
            ("NYHA_II", "nyha_ii"),
            ("NYHA_III", "nyha_iii"),
            ("NYHA_IV", "nyha_iv"),
        ],
    )
    def test_member_value(self, member, value):
        assert HeartFailureClass[member].value == value

    def test_member_count(self):
        assert len(HeartFailureClass) == 4


# ===================================================================
# ENUM TESTS: HeartFailureStage (4 members)
# ===================================================================


class TestHeartFailureStage:
    """Tests for HeartFailureStage enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("STAGE_A", "stage_a"),
            ("STAGE_B", "stage_b"),
            ("STAGE_C", "stage_c"),
            ("STAGE_D", "stage_d"),
        ],
    )
    def test_member_value(self, member, value):
        assert HeartFailureStage[member].value == value

    def test_member_count(self):
        assert len(HeartFailureStage) == 4


# ===================================================================
# ENUM TESTS: EjectionFractionCategory (4 members)
# ===================================================================


class TestEjectionFractionCategory:
    """Tests for EjectionFractionCategory enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("HFrEF", "hfref"),
            ("HFmrEF", "hfmref"),
            ("HFpEF", "hfpef"),
            ("HFimpEF", "hfimpef"),
        ],
    )
    def test_member_value(self, member, value):
        assert EjectionFractionCategory[member].value == value

    def test_member_count(self):
        assert len(EjectionFractionCategory) == 4


# ===================================================================
# ENUM TESTS: ValveSeverity (4 members)
# ===================================================================


class TestValveSeverity:
    """Tests for ValveSeverity enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("MILD", "mild"),
            ("MODERATE", "moderate"),
            ("SEVERE", "severe"),
            ("CRITICAL", "critical"),
        ],
    )
    def test_member_value(self, member, value):
        assert ValveSeverity[member].value == value

    def test_member_count(self):
        assert len(ValveSeverity) == 4


# ===================================================================
# ENUM TESTS: CADRADSScore (7 members)
# ===================================================================


class TestCADRADSScore:
    """Tests for CADRADSScore enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("CAD_RADS_0", "cad_rads_0"),
            ("CAD_RADS_1", "cad_rads_1"),
            ("CAD_RADS_2", "cad_rads_2"),
            ("CAD_RADS_3", "cad_rads_3"),
            ("CAD_RADS_4A", "cad_rads_4a"),
            ("CAD_RADS_4B", "cad_rads_4b"),
            ("CAD_RADS_5", "cad_rads_5"),
        ],
    )
    def test_member_value(self, member, value):
        assert CADRADSScore[member].value == value

    def test_member_count(self):
        assert len(CADRADSScore) == 7


# ===================================================================
# ENUM TESTS: AnticoagulationRecommendation (4 members)
# ===================================================================


class TestAnticoagulationRecommendation:
    """Tests for AnticoagulationRecommendation enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("NO_ANTICOAG", "no_anticoag"),
            ("CONSIDER_ANTICOAG", "consider_anticoag"),
            ("ANTICOAG_RECOMMENDED", "anticoag_recommended"),
            ("ANTICOAG_CONTRAINDICATED", "anticoag_contraindicated"),
        ],
    )
    def test_member_value(self, member, value):
        assert AnticoagulationRecommendation[member].value == value

    def test_member_count(self):
        assert len(AnticoagulationRecommendation) == 4


# ===================================================================
# ENUM TESTS: GDMTPillar (4 members)
# ===================================================================


class TestGDMTPillar:
    """Tests for GDMTPillar enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("BETA_BLOCKER", "beta_blocker"),
            ("ARNI_ACEI_ARB", "arni_acei_arb"),
            ("MRA", "mra"),
            ("SGLT2I", "sglt2i"),
        ],
    )
    def test_member_value(self, member, value):
        assert GDMTPillar[member].value == value

    def test_member_count(self):
        assert len(GDMTPillar) == 4


# ===================================================================
# ENUM TESTS: GDMTStatus (6 members)
# ===================================================================


class TestGDMTStatus:
    """Tests for GDMTStatus enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("NOT_STARTED", "not_started"),
            ("INITIATED", "initiated"),
            ("UPTITRATING", "uptitrating"),
            ("AT_TARGET", "at_target"),
            ("CONTRAINDICATED", "contraindicated"),
            ("INTOLERANT", "intolerant"),
        ],
    )
    def test_member_value(self, member, value):
        assert GDMTStatus[member].value == value

    def test_member_count(self):
        assert len(GDMTStatus) == 6


# ===================================================================
# ENUM TESTS: CardiotoxicityRisk (4 members)
# ===================================================================


class TestCardiotoxicityRisk:
    """Tests for CardiotoxicityRisk enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("LOW", "low"),
            ("MODERATE", "moderate"),
            ("HIGH", "high"),
            ("VERY_HIGH", "very_high"),
        ],
    )
    def test_member_value(self, member, value):
        assert CardiotoxicityRisk[member].value == value

    def test_member_count(self):
        assert len(CardiotoxicityRisk) == 4


# ===================================================================
# ENUM TESTS: LGEPattern (7 members)
# ===================================================================


class TestLGEPattern:
    """Tests for LGEPattern enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("SUBENDOCARDIAL", "subendocardial"),
            ("MID_WALL", "mid_wall"),
            ("EPICARDIAL", "epicardial"),
            ("TRANSMURAL", "transmural"),
            ("RV_INSERTION", "rv_insertion"),
            ("PATCHY", "patchy"),
            ("NONE", "none"),
        ],
    )
    def test_member_value(self, member, value):
        assert LGEPattern[member].value == value

    def test_member_count(self):
        assert len(LGEPattern) == 7


# ===================================================================
# ENUM TESTS: GuidelineClass (5 members)
# ===================================================================


class TestGuidelineClass:
    """Tests for GuidelineClass enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("CLASS_I", "class_i"),
            ("CLASS_IIA", "class_iia"),
            ("CLASS_IIB", "class_iib"),
            ("CLASS_III_NO_BENEFIT", "class_iii_no_benefit"),
            ("CLASS_III_HARM", "class_iii_harm"),
        ],
    )
    def test_member_value(self, member, value):
        assert GuidelineClass[member].value == value

    def test_member_count(self):
        assert len(GuidelineClass) == 5


# ===================================================================
# ENUM TESTS: EvidenceLevel (5 members)
# ===================================================================


class TestEvidenceLevel:
    """Tests for EvidenceLevel enum."""

    @pytest.mark.parametrize(
        "member, value",
        [
            ("LEVEL_A", "level_a"),
            ("LEVEL_B_R", "level_b_r"),
            ("LEVEL_B_NR", "level_b_nr"),
            ("LEVEL_C_LD", "level_c_ld"),
            ("LEVEL_C_EO", "level_c_eo"),
        ],
    )
    def test_member_value(self, member, value):
        assert EvidenceLevel[member].value == value

    def test_member_count(self):
        assert len(EvidenceLevel) == 5


# ===================================================================
# PYDANTIC MODEL TESTS: CardioQuery
# ===================================================================


class TestCardioQuery:
    """Tests for CardioQuery Pydantic model."""

    def test_basic_creation(self):
        q = CardioQuery(question="What is HFrEF?")
        assert q.question == "What is HFrEF?"
        assert q.workflow_type is None
        assert q.patient_context is None

    def test_with_workflow_type(self):
        q = CardioQuery(
            question="How to manage HFrEF?",
            workflow_type=CardioWorkflowType.HEART_FAILURE,
        )
        assert q.workflow_type == CardioWorkflowType.HEART_FAILURE

    def test_with_patient_context(self):
        ctx = {"age": 65, "sex": "male", "lvef": 30}
        q = CardioQuery(question="GDMT optimization", patient_context=ctx)
        assert q.patient_context["lvef"] == 30

    def test_empty_question_raises(self):
        with pytest.raises(ValidationError):
            CardioQuery(question="")

    def test_question_required(self):
        with pytest.raises(ValidationError):
            CardioQuery()


# ===================================================================
# PYDANTIC MODEL TESTS: CardioSearchResult
# ===================================================================


class TestCardioSearchResult:
    """Tests for CardioSearchResult Pydantic model."""

    def test_basic_creation(self):
        r = CardioSearchResult(
            collection="cardio_literature",
            content="Study on HFrEF treatment",
            score=0.95,
        )
        assert r.collection == "cardio_literature"
        assert r.score == 0.95
        assert r.metadata == {}

    def test_with_metadata(self):
        r = CardioSearchResult(
            collection="cardio_trials",
            content="DAPA-HF",
            score=0.8,
            metadata={"nct_id": "NCT03036124"},
        )
        assert r.metadata["nct_id"] == "NCT03036124"

    def test_score_ge_zero_validation(self):
        with pytest.raises(ValidationError):
            CardioSearchResult(
                collection="x", content="y", score=-0.1,
            )

    def test_score_zero_valid(self):
        r = CardioSearchResult(collection="x", content="y", score=0.0)
        assert r.score == 0.0


# ===================================================================
# PYDANTIC MODEL TESTS: RiskScoreInput
# ===================================================================


class TestRiskScoreInput:
    """Tests for RiskScoreInput Pydantic model."""

    def test_minimal_creation(self):
        inp = RiskScoreInput(score_type=RiskScoreType.ASCVD)
        assert inp.score_type == RiskScoreType.ASCVD
        assert inp.age is None

    def test_full_ascvd_input(self):
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=55, sex="male", race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        assert inp.age == 55
        assert inp.sex == "male"
        assert inp.total_cholesterol == 213.0

    def test_sex_validation_lowercase(self):
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            sex="Male",
        )
        assert inp.sex == "male"

    def test_sex_validation_invalid(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(
                score_type=RiskScoreType.ASCVD,
                sex="unknown",
            )

    def test_age_range_validation_low(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.ASCVD, age=17)

    def test_age_range_validation_high(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.ASCVD, age=121)

    def test_age_boundary_18(self):
        inp = RiskScoreInput(score_type=RiskScoreType.ASCVD, age=18)
        assert inp.age == 18

    def test_age_boundary_120(self):
        inp = RiskScoreInput(score_type=RiskScoreType.ASCVD, age=120)
        assert inp.age == 120

    def test_systolic_bp_range(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.ASCVD, systolic_bp=59)

    def test_cholesterol_range(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.ASCVD, total_cholesterol=49)

    def test_lvef_range_low(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.MAGGIC, lvef=4.9)

    def test_lvef_range_high(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.MAGGIC, lvef=90.1)

    def test_lvef_valid_boundaries(self):
        inp = RiskScoreInput(score_type=RiskScoreType.MAGGIC, lvef=5.0)
        assert inp.lvef == 5.0
        inp2 = RiskScoreInput(score_type=RiskScoreType.MAGGIC, lvef=90.0)
        assert inp2.lvef == 90.0

    def test_creatinine_range(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.MAGGIC, creatinine=0.09)

    def test_hdl_range(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.ASCVD, hdl=9.0)

    def test_bmi_range(self):
        with pytest.raises(ValidationError):
            RiskScoreInput(score_type=RiskScoreType.MAGGIC, bmi=9.0)

    def test_nyha_class_field(self):
        inp = RiskScoreInput(
            score_type=RiskScoreType.MAGGIC,
            nyha_class=HeartFailureClass.NYHA_III,
        )
        assert inp.nyha_class == HeartFailureClass.NYHA_III

    def test_defaults_are_none(self):
        inp = RiskScoreInput(score_type=RiskScoreType.ASCVD)
        assert inp.age is None
        assert inp.sex is None
        assert inp.race is None
        assert inp.diabetes is None
        assert inp.smoker is None
        assert inp.lvef is None
        assert inp.history_of_stroke is None


# ===================================================================
# PYDANTIC MODEL TESTS: RiskScoreResult
# ===================================================================


class TestRiskScoreResult:
    """Tests for RiskScoreResult Pydantic model."""

    def test_basic_creation(self):
        r = RiskScoreResult(
            score_type=RiskScoreType.ASCVD,
            score_value=7.5,
            risk_category="Intermediate",
            interpretation="10-year risk is 7.5%",
        )
        assert r.score_value == 7.5
        assert r.risk_category == "Intermediate"
        assert r.recommendations == []
        assert r.guideline_reference == ""

    def test_with_recommendations(self):
        r = RiskScoreResult(
            score_type=RiskScoreType.HEART,
            score_value=3.0,
            risk_category="Low",
            interpretation="Low risk",
            recommendations=["Consider discharge", "Follow up"],
        )
        assert len(r.recommendations) == 2

    def test_with_guideline_reference(self):
        r = RiskScoreResult(
            score_type=RiskScoreType.CHA2DS2_VASC,
            score_value=4.0,
            risk_category="High",
            interpretation="High stroke risk",
            guideline_reference="Lip GYH, Chest 2010",
        )
        assert "Lip" in r.guideline_reference


# ===================================================================
# PYDANTIC MODEL TESTS: GDMTMedication
# ===================================================================


class TestGDMTMedication:
    """Tests for GDMTMedication Pydantic model."""

    def test_basic_creation(self):
        m = GDMTMedication(
            pillar=GDMTPillar.BETA_BLOCKER,
            drug_name="carvedilol",
            target_dose="25 mg BID",
        )
        assert m.pillar == GDMTPillar.BETA_BLOCKER
        assert m.drug_name == "carvedilol"
        assert m.status == GDMTStatus.NOT_STARTED
        assert m.current_dose is None
        assert m.contraindications == []

    def test_with_all_fields(self):
        m = GDMTMedication(
            pillar=GDMTPillar.ARNI_ACEI_ARB,
            drug_name="sacubitril/valsartan",
            current_dose="49/51 mg BID",
            target_dose="97/103 mg BID",
            status=GDMTStatus.UPTITRATING,
            contraindications=["angioedema history"],
        )
        assert m.status == GDMTStatus.UPTITRATING
        assert len(m.contraindications) == 1


# ===================================================================
# PYDANTIC MODEL TESTS: GDMTRecommendation
# ===================================================================


class TestGDMTRecommendation:
    """Tests for GDMTRecommendation Pydantic model."""

    def test_basic_creation(self):
        r = GDMTRecommendation(
            ef_category=EjectionFractionCategory.HFrEF,
        )
        assert r.ef_category == EjectionFractionCategory.HFrEF
        assert r.current_meds == []
        assert r.recommendations == []
        assert r.next_steps == []
        assert r.guideline_references == []

    def test_with_medications(self):
        med = GDMTMedication(
            pillar=GDMTPillar.SGLT2I,
            drug_name="dapagliflozin",
            target_dose="10 mg daily",
            status=GDMTStatus.AT_TARGET,
        )
        r = GDMTRecommendation(
            ef_category=EjectionFractionCategory.HFrEF,
            current_meds=[med],
            recommendations=["Continue SGLT2i"],
        )
        assert len(r.current_meds) == 1
        assert r.current_meds[0].drug_name == "dapagliflozin"


# ===================================================================
# PYDANTIC MODEL TESTS: ValveAssessment
# ===================================================================


class TestValveAssessment:
    """Tests for ValveAssessment Pydantic model."""

    def test_basic_creation(self):
        v = ValveAssessment(
            valve="aortic",
            pathology="aortic stenosis",
            severity=ValveSeverity.SEVERE,
        )
        assert v.valve == "aortic"
        assert v.severity == ValveSeverity.SEVERE
        assert v.intervention_criteria_met is False
        assert v.measurements == {}
        assert v.recommendation == ""

    def test_with_measurements(self):
        v = ValveAssessment(
            valve="mitral",
            pathology="mitral regurgitation",
            severity=ValveSeverity.MODERATE,
            measurements={"EROA": 0.25, "RVol": 40},
            intervention_criteria_met=False,
        )
        assert v.measurements["EROA"] == 0.25


# ===================================================================
# PYDANTIC MODEL TESTS: ECGInterpretation
# ===================================================================


class TestECGInterpretation:
    """Tests for ECGInterpretation Pydantic model."""

    def test_basic_creation(self):
        ecg = ECGInterpretation(
            rhythm="Normal sinus rhythm",
            rate=72,
        )
        assert ecg.rhythm == "Normal sinus rhythm"
        assert ecg.rate == 72
        assert ecg.intervals == {}
        assert ecg.axis == ""
        assert ecg.findings == []
        assert ecg.urgency == SeverityLevel.INFORMATIONAL

    def test_with_findings(self):
        ecg = ECGInterpretation(
            rhythm="Atrial fibrillation",
            rate=110,
            intervals={"QTc": 480},
            findings=["Rapid AF", "ST depression lateral leads"],
            urgency=SeverityLevel.HIGH,
        )
        assert len(ecg.findings) == 2
        assert ecg.urgency == SeverityLevel.HIGH

    def test_rate_validation_low(self):
        with pytest.raises(ValidationError):
            ECGInterpretation(rhythm="test", rate=-1)

    def test_rate_validation_high(self):
        with pytest.raises(ValidationError):
            ECGInterpretation(rhythm="test", rate=401)


# ===================================================================
# PYDANTIC MODEL TESTS: ImagingResult
# ===================================================================


class TestImagingResult:
    """Tests for ImagingResult Pydantic model."""

    def test_basic_creation(self):
        img = ImagingResult(modality=ImagingModality.ECHOCARDIOGRAPHY)
        assert img.modality == ImagingModality.ECHOCARDIOGRAPHY
        assert img.findings == []
        assert img.measurements == {}
        assert img.impression == ""
        assert img.cross_modal_triggers == []

    def test_with_measurements(self):
        img = ImagingResult(
            modality=ImagingModality.CARDIAC_MRI,
            findings=["Mid-wall LGE", "Reduced LVEF"],
            measurements={"LVEF": 35, "native_T1": 1100},
            impression="Dilated cardiomyopathy with fibrosis",
            cross_modal_triggers=["Genetic testing for DCM panel"],
        )
        assert len(img.findings) == 2
        assert img.measurements["native_T1"] == 1100


# ===================================================================
# PYDANTIC MODEL TESTS: CardiotoxicityAssessment
# ===================================================================


class TestCardiotoxicityAssessment:
    """Tests for CardiotoxicityAssessment Pydantic model."""

    def test_basic_creation(self):
        ca = CardiotoxicityAssessment(
            agent="doxorubicin",
            risk_level=CardiotoxicityRisk.HIGH,
        )
        assert ca.agent == "doxorubicin"
        assert ca.risk_level == CardiotoxicityRisk.HIGH
        assert ca.baseline_lvef is None
        assert ca.current_lvef is None
        assert ca.recommendation == ""

    def test_with_lvef_tracking(self):
        ca = CardiotoxicityAssessment(
            agent="trastuzumab",
            risk_level=CardiotoxicityRisk.MODERATE,
            baseline_lvef=62.0,
            current_lvef=50.0,
            gls_decline_percent=18.0,
        )
        assert ca.baseline_lvef == 62.0
        assert ca.current_lvef == 50.0

    def test_lvef_validation(self):
        with pytest.raises(ValidationError):
            CardiotoxicityAssessment(
                agent="doxorubicin",
                risk_level=CardiotoxicityRisk.HIGH,
                baseline_lvef=4.0,
            )


# ===================================================================
# PYDANTIC MODEL TESTS: CrossModalTrigger
# ===================================================================


class TestCrossModalTrigger:
    """Tests for CrossModalTrigger Pydantic model."""

    def test_basic_creation(self):
        t = CrossModalTrigger(
            trigger_source="cardiac MRI",
            finding="Unexplained LVH >=15mm",
        )
        assert t.trigger_source == "cardiac MRI"
        assert t.gene_panel == []
        assert t.conditions == []
        assert t.rationale == ""

    def test_with_gene_panel(self):
        t = CrossModalTrigger(
            trigger_source="echocardiogram",
            finding="Asymmetric septal hypertrophy",
            gene_panel=["MYH7", "MYBPC3", "TNNT2"],
            conditions=["HCM"],
            rationale="Suspect sarcomeric HCM",
        )
        assert "MYH7" in t.gene_panel
        assert len(t.gene_panel) == 3


# ===================================================================
# PYDANTIC MODEL TESTS: WorkflowResult
# ===================================================================


class TestWorkflowResult:
    """Tests for WorkflowResult Pydantic model."""

    def test_basic_creation(self):
        w = WorkflowResult(
            workflow_type=CardioWorkflowType.HEART_FAILURE,
        )
        assert w.workflow_type == CardioWorkflowType.HEART_FAILURE
        assert w.findings == []
        assert w.risk_scores == []
        assert w.recommendations == []
        assert w.severity == SeverityLevel.INFORMATIONAL

    def test_with_risk_scores(self):
        score = RiskScoreResult(
            score_type=RiskScoreType.MAGGIC,
            score_value=25.0,
            risk_category="High",
            interpretation="High 1-year mortality risk",
        )
        w = WorkflowResult(
            workflow_type=CardioWorkflowType.HEART_FAILURE,
            risk_scores=[score],
            severity=SeverityLevel.HIGH,
        )
        assert len(w.risk_scores) == 1
        assert w.severity == SeverityLevel.HIGH


# ===================================================================
# PYDANTIC MODEL TESTS: CardioResponse
# ===================================================================


class TestCardioResponse:
    """Tests for CardioResponse Pydantic model."""

    def test_basic_creation(self):
        r = CardioResponse(answer="HFrEF requires 4-pillar GDMT.")
        assert r.answer == "HFrEF requires 4-pillar GDMT."
        assert r.citations == []
        assert r.risk_scores == []
        assert r.workflow_results == []
        assert r.cross_modal_triggers == []
        assert r.confidence == 0.0

    def test_confidence_range_valid(self):
        r = CardioResponse(answer="test", confidence=0.95)
        assert r.confidence == 0.95

    def test_confidence_too_high(self):
        with pytest.raises(ValidationError):
            CardioResponse(answer="test", confidence=1.1)

    def test_confidence_too_low(self):
        with pytest.raises(ValidationError):
            CardioResponse(answer="test", confidence=-0.1)

    def test_confidence_boundaries(self):
        r0 = CardioResponse(answer="test", confidence=0.0)
        r1 = CardioResponse(answer="test", confidence=1.0)
        assert r0.confidence == 0.0
        assert r1.confidence == 1.0


# ===================================================================
# DATACLASS TESTS: SearchPlan
# ===================================================================


class TestSearchPlan:
    """Tests for SearchPlan dataclass."""

    def test_default_creation(self):
        sp = SearchPlan()
        assert sp.question == ""
        assert sp.conditions == []
        assert sp.drugs == []
        assert sp.imaging_modalities == []
        assert sp.relevant_workflows == []
        assert sp.search_strategy == "broad"
        assert sp.sub_questions == []
        assert sp.identified_topics == []

    def test_creation_with_values(self):
        sp = SearchPlan(
            question="How to manage severe AS?",
            conditions=["Aortic Stenosis"],
            drugs=[],
            imaging_modalities=["echocardiography"],
            relevant_workflows=[CardioWorkflowType.VALVULAR_DISEASE],
            search_strategy="targeted",
        )
        assert sp.question == "How to manage severe AS?"
        assert sp.search_strategy == "targeted"
        assert len(sp.relevant_workflows) == 1

    def test_mutable_defaults_are_independent(self):
        sp1 = SearchPlan()
        sp2 = SearchPlan()
        sp1.conditions.append("HCM")
        assert sp2.conditions == []

    def test_search_strategy_values(self):
        for strategy in ("broad", "targeted", "comparative", "clinical"):
            sp = SearchPlan(search_strategy=strategy)
            assert sp.search_strategy == strategy
