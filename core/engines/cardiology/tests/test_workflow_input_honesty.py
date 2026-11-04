"""Regression tests for the cardiology clinical workflows' safety defects.

Found by the pre-release audit (identifiers B-3, B-4, B-5, B-11, B-12, B-13,
B-15 and B-16).  The common thread is an absent input reported as a
measurement: ten workflow runs with empty or partial input produced ZERO
`[INPUT WARNING]` lines while printing fabricated T1/T2/ECV values, vitals,
calcium scores and culprit vessels as findings.
"""
import pytest

from src.clinical_workflows import (
    AcuteDecompensatedHFWorkflow,
    CADAssessmentWorkflow,
    CardiacMRIWorkflow,
    CardioOncologyWorkflow,
    HeartFailureWorkflow,
    MyocarditisPericarditisWorkflow,
    PostMIWorkflow,
    PreventiveRiskWorkflow,
    StressTestWorkflow,
    ValvularDiseaseWorkflow,
)
from src.gdmt_optimizer import GDMTOptimizer
from src.hf_thresholds import CRT_QRS_CLASS_IIA_MS
from src.models import RiskScoreType, SeverityLevel, ValveSeverity


# =====================================================================
# B-3 — workflows invented normal measurements and printed them as findings
# =====================================================================

class TestB3AbsentInputsAreNotReportedAsMeasurements:
    """Every `preprocess` filled absent inputs with plausible in-range values
    via bare `setdefault`, and the validators that ran afterwards therefore saw
    a valid number and appended no warning.  `execute` then printed the
    placeholder as a finding.
    """

    def test_cad_absent_calcium_score_is_not_a_calcium_score_of_zero(self):
        result = CADAssessmentWorkflow().run({})
        text = " ".join(result.findings)
        assert "Calcium score 0: very low risk (zero plaque burden)" not in text
        assert "CAD-RADS 0: 0% — no plaque or stenosis" not in text
        assert "INDETERMINATE" in text
        assert not any("No further testing" in r for r in result.recommendations)

    def test_cmr_without_parametric_mapping_reports_no_t1_t2_or_ecv_value(self):
        """Mapping is not universally performed, so a CMR report without it is
        routine.  This one came back asserting three normal parametric values, a
        reader could use to rule out amyloid or active myocarditis on numbers
        the scanner never produced."""
        result = CardiacMRIWorkflow().run({"lvef": 32.0, "lge_pattern": "mid_wall"})
        text = " ".join(result.findings)
        assert "Native T1: 1000 ms — within normal range" not in text
        assert "T2: 45 ms — within normal range" not in text
        assert "ECV: 27.0% — within normal range" not in text
        assert "LGE extent: 0.0% of LV mass" not in text
        assert "Native T1 mapping: NOT PERFORMED" in text
        assert "T2 mapping: NOT PERFORMED" in text
        assert "Extracellular volume (ECV): NOT PERFORMED" in text

    def test_adhf_without_vitals_produces_no_furosemide_dose_order(self):
        result = AcuteDecompensatedHFWorkflow().run({})
        findings = " ".join(result.findings)
        recs = " ".join(result.recommendations)
        assert "SBP 110.0 mmHg, HR 90 bpm, LVEF 35.0%" not in findings
        assert "Warm-Wet (most common ADHF profile" not in findings
        assert "furosemide 40-80 mg IV" not in recs
        assert "Target net negative 1-2 L/day" not in recs
        assert "Hemodynamic profile: INDETERMINATE" in findings

    def test_stress_test_without_an_exercise_time_reports_no_duke_score(self):
        result = StressTestWorkflow().run({})
        text = " ".join(result.findings)
        assert "Exercise time: 0.0 min" not in text
        assert "Poor exercise capacity (<5 METs)" not in text
        assert "Duke Treadmill Score: 0.0" not in text
        assert "Duke Treadmill Score: INDETERMINATE" in text

    def test_heart_failure_without_an_ef_assigns_no_phenotype(self):
        result = HeartFailureWorkflow().run({})
        text = " ".join(result.findings)
        assert "LVEF 55.0% — classified as hfpef" not in text
        assert "NYHA class nyha_i" not in text
        assert "ACC/AHA stage stage_a" not in text
        assert "LVEF: INDETERMINATE" in text

    def test_post_mi_without_an_angiogram_invents_no_culprit_vessel(self):
        """Because the culprit defaulted to the LAD, an empty request also
        produced "ACEi/ARB within 24 hours (anterior MI or LVEF <40%)"."""
        result = PostMIWorkflow().run({})
        findings = " ".join(result.findings)
        recs = " ".join(result.recommendations)
        assert "NSTEMI — culprit vessel: LAD" not in findings
        assert "Time from symptom onset: 6.0 hours" not in findings
        assert "LVEF: 50.0%" not in findings
        assert "ACEi/ARB within 24 hours (anterior MI or LVEF <40%)" not in recs

    def test_cardio_oncology_without_an_echo_reports_no_gls_or_lvef(self):
        result = CardioOncologyWorkflow().run(
            {"chemotherapy_agent": "anthracycline"}
        )
        text = " ".join(result.findings)
        assert "Baseline GLS: -20.0%, Current GLS: -19.0%" not in text
        assert "Baseline LVEF: 60.0%, Current LVEF: 58.0%" not in text
        assert "LVEF remains within acceptable range" not in text
        assert "Overall cardiotoxicity risk: low" not in text
        assert "Overall cardiotoxicity risk: INDETERMINATE" in text

    def test_myocarditis_without_inputs_invents_no_presentation_or_ef(self):
        result = MyocarditisPericarditisWorkflow().run({})
        text = " ".join(result.findings)
        assert "Symptom onset: 7 days ago" not in text
        assert "LVEF: 55.0%" not in text

    @pytest.mark.parametrize(
        "workflow",
        [
            CADAssessmentWorkflow, HeartFailureWorkflow, CardiacMRIWorkflow,
            StressTestWorkflow, CardioOncologyWorkflow,
            AcuteDecompensatedHFWorkflow, PostMIWorkflow,
            MyocarditisPericarditisWorkflow, ValvularDiseaseWorkflow,
            PreventiveRiskWorkflow,
        ],
    )
    def test_an_absent_input_raises_an_input_warning(self, workflow):
        """`BaseCardioWorkflow.run` injects `[INPUT WARNING]` findings from
        `_validation_warnings`, but `setdefault` ran first so nothing was ever
        appended for an absent key: a run over all ten workflows with empty or
        partial input produced zero warnings."""
        result = workflow().run({})
        assert any(
            f.startswith("[INPUT WARNING]") and "not supplied" in f
            for f in result.findings
        ), "absent inputs must be announced, not silently substituted"


# =====================================================================
# B-4 — unreported measurements voted for "mild"
# =====================================================================

class TestB4ValveGradingDoesNotCountUnreportedMeasurements:
    """`peak_velocity=0.0`, `mean_gradient=0.0` and `valve_area=999.0` all sit
    inside the MILD windows (`mean_grad (0,19)`, `ava (1.5,999)`), so two absent
    measurements voted for MILD and outvoted one real severe measurement.  A
    symptomatic patient with Vmax 4.2 m/s -- symptomatic severe aortic stenosis,
    roughly 50% two-year mortality untreated -- was graded "mild", severity low,
    with no AVR recommendation of any kind, and a mean gradient of 0.0 mmHg and
    a valve area of 999 cm2 were printed as findings.
    """

    def test_symptomatic_severe_as_on_peak_velocity_alone_gets_avr(self):
        result = ValvularDiseaseWorkflow().run({
            "valve": "aortic", "pathology": "stenosis",
            "measurements": {"peak_velocity": 4.2},
            "symptoms": ["exertional dyspnoea", "angina"],
        })
        findings = " ".join(result.findings)
        recs = " ".join(result.recommendations)
        assert "Aortic stenosis severity: mild" not in findings
        assert "Aortic stenosis severity: severe" in findings
        assert "Mean gradient: 0.0 mmHg" not in findings
        assert "Aortic valve area: 999.0 cm2" not in findings
        assert "AVR indicated (Class I, Level A)" in recs
        assert result.severity == SeverityLevel.HIGH

    def test_severe_mr_on_ero_alone_is_not_graded_mild(self):
        severity, findings = ValvularDiseaseWorkflow()._grade_mitral_regurgitation(
            {"ero": 0.45}
        )
        assert severity == ValveSeverity.SEVERE
        assert not any("Regurgitant volume: 0.0 mL" in f for f in findings)
        assert any("Not reported" in f for f in findings)

    def test_no_measurement_at_all_is_indeterminate_not_mild(self):
        severity, findings = ValvularDiseaseWorkflow()._grade_aortic_stenosis({})
        assert severity is None
        assert any("INDETERMINATE" in f for f in findings)

    def test_three_concordant_measurements_still_grade_as_before(self):
        grade = ValvularDiseaseWorkflow()._grade_aortic_stenosis
        assert grade(
            {"peak_velocity": 2.5, "mean_gradient": 15, "valve_area": 1.8}
        )[0] == ValveSeverity.MILD
        assert grade(
            {"peak_velocity": 3.5, "mean_gradient": 30, "valve_area": 1.2}
        )[0] == ValveSeverity.MODERATE
        assert grade(
            {"peak_velocity": 4.5, "mean_gradient": 50, "valve_area": 0.7}
        )[0] == ValveSeverity.SEVERE

    def test_severe_mr_with_no_ef_does_not_trigger_lv_dysfunction_surgery(self):
        """The LVEF default of exactly 60.0 landed on the `lvef <= 60` trigger,
        so an asymptomatic severe MR with no EF on record was told surgery was
        indicated "for LV dysfunction"."""
        result = ValvularDiseaseWorkflow().run({
            "valve": "mitral", "pathology": "regurgitation",
            "measurements": {"ero": 0.45, "regurgitant_volume": 70,
                             "vena_contracta": 0.8},
        })
        recs = " ".join(result.recommendations)
        assert "Severe MR with LV dysfunction (LVEF <=60%) — surgery indicated" not in recs
        assert "measure LVEF" in recs


# =====================================================================
# B-5 — the preventive workflow invented an ASCVD percentage
# =====================================================================

class TestB5PreventiveWorkflowUsesTheValidatedPCE:
    """`_estimate_ascvd_risk` was a hand-written point model ("calibrated to PCE
    outputs") whose own docstring admitted it was "not a replacement for the full
    validated calculator" -- but its output was published as the real thing,
    typed `RiskScoreType.ASCVD` and cited to the 2019 ACC/AHA Primary Prevention
    Guideline.  It reported 10.7% for a 70-year-old diabetic woman whose risk is
    33.7% by this repo's own validated Pooled Cohort Equations, and flipped the
    statin band in four of six reference cases, every one toward
    under-treatment.  Its race term was a flat `base_risk *= 1.5`, which the
    Pooled Cohort Equations do not contain in any cohort.
    """

    def test_the_hand_written_estimator_is_gone(self):
        assert not hasattr(PreventiveRiskWorkflow, "_estimate_ascvd_risk")

    def test_high_risk_diabetic_woman_gets_a_high_intensity_statin(self):
        result = PreventiveRiskWorkflow().run({
            "age": 70, "sex": "female", "race": "white",
            "total_cholesterol": 200, "hdl": 50, "ldl": 120,
            "systolic_bp": 155, "bp_treatment": True,
            "diabetes": True, "smoker": False,
        })
        risk = next(f for f in result.findings if "10-year ASCVD risk" in f)
        assert "10.7%" not in risk
        assert "High risk category" in risk
        assert any(
            "high-intensity statin" in r.lower() for r in result.recommendations
        )

    def test_risk_matches_the_validated_calculator_exactly(self):
        result = PreventiveRiskWorkflow().run({
            "age": 60, "sex": "male", "race": "white",
            "total_cholesterol": 250, "hdl": 35, "ldl": 150,
            "systolic_bp": 145, "bp_treatment": False,
            "diabetes": False, "smoker": True,
        })
        score = next(
            s for s in result.risk_scores if s.score_type == RiskScoreType.ASCVD
        )

        from src.models import RiskScoreInput
        from src.risk_calculators import calculate_ascvd

        reference = calculate_ascvd(RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=60, sex="male", race="white", total_cholesterol=250, hdl=35,
            systolic_bp=145, hypertension_treatment=False, diabetes=False,
            smoker=True,
        ))
        assert score.score_value == reference.score_value
        assert score.score_value > 20.0

    def test_no_risk_number_is_published_when_the_inputs_are_absent(self):
        result = PreventiveRiskWorkflow().run({})
        assert result.risk_scores == []
        risk = next(f for f in result.findings if "10-year ASCVD risk" in f)
        assert "INDETERMINATE" in risk
        assert not any("Low ASCVD risk" in r for r in result.recommendations)

    def test_age_outside_the_validated_range_is_still_reported(self):
        result = PreventiveRiskWorkflow().run({"age": 30})
        assert any("outside the validated range" in f for f in result.findings)
        assert result.risk_scores == []


# =====================================================================
# B-11 — this workflow disagreed with the optimizer about the CRT Class 2a
# QRS threshold and the SGLT2i class
# =====================================================================

def test_b11_crt_qrs_125_lbbb_is_recommended_by_both_surfaces():
    """`qrs_dur >= 130` here against `>= 120` in GDMTOptimizer: a patient with
    LVEF 32%, QRS 125 ms and LBBB got a Class 2a CRT recommendation from
    /gdmt/optimize and nothing at all from /workflow/heart-failure."""
    workflow = HeartFailureWorkflow().run({
        "lvef": 32.0, "nyha_class": "II",
        "qrs_duration": 125, "qrs_morphology": "lbbb",
    })
    optimizer = GDMTOptimizer().optimize(
        lvef=32.0, nyha_class="II", current_medications=[],
        patient_data={
            "sbp": 120.0, "hr": 70, "potassium": 4.0, "egfr": 60.0,
            "qrs_duration": 125, "qrs_morphology": "LBBB", "rhythm": "sinus",
        },
    )
    assert any("CRT" in r for r in workflow.recommendations)
    assert "CRT" in " ".join(optimizer.recommendations)
    assert CRT_QRS_CLASS_IIA_MS == 120


def test_b11_workflow_sglt2i_class_outside_hfref_is_2a():
    for lvef in (45.0, 60.0):
        result = HeartFailureWorkflow().run({"lvef": lvef, "nyha_class": "II"})
        sglt2i = [r for r in result.recommendations if "SGLT2" in r]
        assert sglt2i
        assert all("Class 2a" in r for r in sglt2i)


# =====================================================================
# B-12 / B-13 / B-15 / B-16 — false claims and dead code
# =====================================================================

def test_b12_symptom_onset_is_not_labelled_door_to_balloon():
    """`time_from_onset_hours` is documented as hours from symptom onset but was
    printed as door-to-balloon / door-to-needle time -- different quality
    metrics, measured from first hospital contact.  A STEMI patient who
    presented 5 hours after onset with a 42-minute door-to-balloon time was
    reported as "Door-to-balloon time 5.0h — exceeds 90-minute target", which
    would corrupt any D2B audit built on this output."""
    result = PostMIWorkflow().run({
        "mi_type": "STEMI", "reperfusion_strategy": "PCI",
        "time_from_onset_hours": 5.0, "door_to_balloon_min": 42,
        "lvef": 45.0, "culprit_vessel": "RCA",
    })
    text = " ".join(result.findings)
    assert "Door-to-balloon time 5.0h" not in text
    assert "Door-to-balloon time 42 min — within 90-minute target" in text
    assert "Symptom-onset-to-presentation interval: 5.0h" in text


def test_b12_absent_door_to_balloon_is_not_asserted_as_a_failure():
    """The default of 6.0 hours meant an empty request asserted a D2B
    failure."""
    result = PostMIWorkflow().run({
        "mi_type": "STEMI", "reperfusion_strategy": "PCI", "lvef": 45.0,
    })
    text = " ".join(result.findings)
    assert "exceeds 90-minute target" not in text
    assert "Door-to-balloon time: not reported" in text


def test_b13_pericarditis_pain_and_rub_can_be_met():
    """Criteria 1 (pleuritic chest pain) and 2 (pericardial rub) were searched
    for in `ecg_findings`, documented as the ECG findings list, so neither could
    ever be met: a textbook patient with pain, a rub, diffuse ST elevation with
    PR depression and a small effusion scored 2/4 -- exactly the diagnostic
    threshold -- with no note that two criteria had never been assessed."""
    result = MyocarditisPericarditisWorkflow().run({
        "presentation": "pericarditis",
        "symptoms": ["pleuritic_chest_pain"],
        "exam_findings": ["pericardial_rub"],
        "ecg_findings": ["diffuse_st_elevation", "pr_depression"],
        "cmr_findings": {"pericardial_effusion": True},
    })
    criteria = next(f for f in result.findings if "diagnostic criteria" in f)
    assert "4/4 met" in criteria
    assert "pleuritic chest pain" in criteria
    assert "pericardial rub" in criteria


def test_b13_unassessed_criteria_are_named_not_counted_as_absent():
    result = MyocarditisPericarditisWorkflow().run({
        "presentation": "pericarditis",
        "ecg_findings": ["diffuse_st_elevation", "pr_depression"],
        "cmr_findings": {"pericardial_effusion": True},
    })
    criteria = next(f for f in result.findings if "diagnostic criteria" in f)
    assert "2/4 met" in criteria
    assert "NOT ASSESSED" in criteria
    assert "pleuritic chest pain" in criteria
    assert "pericardial rub" in criteria


def test_b15_post_mi_beta_blocker_check_is_not_discarded():
    """The `any(...)` result was computed and thrown away by a bare expression
    with no assignment, so "Beta-blocker within 24 hours" was appended
    unconditionally, including for a patient already taking carvedilol."""
    on_bb = PostMIWorkflow().run({
        "mi_type": "NSTEMI", "lvef": 45.0, "culprit_vessel": "RCA",
        "current_medications": ["carvedilol 12.5 mg BID", "aspirin"],
    })
    assert not any(
        "Beta-blocker within 24 hours" in r for r in on_bb.recommendations
    )
    assert any("Already on a beta-blocker" in f for f in on_bb.findings)

    not_on_bb = PostMIWorkflow().run({
        "mi_type": "NSTEMI", "lvef": 45.0, "culprit_vessel": "RCA",
        "current_medications": ["aspirin"],
    })
    assert any(
        "Beta-blocker within 24 hours" in r for r in not_on_bb.recommendations
    )


def test_b16_stress_hrr_prose_matches_its_boundary():
    """`if hr_rec <= 12` printed "(normal >=12 bpm)", calling 12 both normal and
    abnormal in one sentence.  The boundary is right (Cole et al. define
    abnormal as a drop of <=12 bpm); the prose was wrong."""
    result = StressTestWorkflow().run({
        "test_type": "exercise", "exercise_time_min": 9.0,
        "max_st_deviation_mm": 0.0, "angina_index": 0,
        "peak_hr": 150, "percent_mphr": 90.0, "bp_response": "normal",
        "perfusion_defect": "none", "hr_recovery_1min": 12,
    })
    hrr = next(f for f in result.findings if "heart rate recovery" in f)
    assert "Abnormal heart rate recovery" in hrr
    assert "normal >12 bpm" in hrr
    assert "normal >=12 bpm" not in hrr
