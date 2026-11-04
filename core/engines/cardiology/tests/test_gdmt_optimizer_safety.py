"""Regression tests for the GDMT optimizer's clinical-safety defects.

Found by the pre-release audit (identifiers B-1, B-2, B-8, B-9, B-10, B-11 and
B-14), every one of them present with 2,099 tests passing.  Each test names the
wrong output it asserts is gone, and asserts on the DECISION -- the drug order,
the dose, the device class -- rather than on vocabulary appearing somewhere in
the output, which is how these reached a release candidate.
"""
from src.gdmt_optimizer import GDMTOptimizer, _dose_from_med
from src.hf_thresholds import (
    CRT_QRS_CLASS_I_MS,
    CRT_QRS_CLASS_IIA_MS,
    HFPEF_LVEF_CUTOFF,
    HFREF_LVEF_CUTOFF,
    SGLT2I_CLASS_HFMREF,
    SGLT2I_CLASS_HFPEF,
)
from src.models import (
    EjectionFractionCategory,
    GDMTPillar,
    GDMTStatus,
    GuidelineClass,
)

import pytest


# =====================================================================
# B-1 (optimizer side) — an absent lab or vital is not a normal one
# =====================================================================

class TestB1AbsentLabsAndVitalsAreNotTreatedAsNormal:
    """With nothing but an EF and a functional class, the optimizer used to
    issue unconditional "INITIATE BETA-BLOCKER", "INITIATE ARNI" and
    "INITIATE MRA" orders, because every safety check was written
    `if value is not None and value < threshold` -- so absence read as benign.
    """

    @pytest.fixture
    def recs(self):
        result = GDMTOptimizer().optimize(
            lvef=25.0, nyha_class="III", current_medications=[], patient_data={},
        )
        return " ".join(result.recommendations)

    def test_no_beta_blocker_order_without_hr_or_bp(self, recs):
        assert "INITIATE BETA-BLOCKER: Start carvedilol" not in recs
        assert "BETA-BLOCKER: INDETERMINATE" in recs
        assert "resting heart rate" in recs

    def test_no_arni_order_without_bp_or_egfr(self, recs):
        assert "INITIATE ARNI: Start sacubitril/valsartan" not in recs
        assert "ARNI/ACEi/ARB: INDETERMINATE" in recs

    def test_no_mra_order_without_potassium_or_egfr(self, recs):
        assert "INITIATE MRA: Start spironolactone" not in recs
        assert "MRA: INDETERMINATE" in recs
        assert "serum potassium" in recs

    def test_mra_is_still_deferred_when_the_potassium_is_high(self):
        result = GDMTOptimizer().optimize(
            lvef=25.0, nyha_class="III", current_medications=[],
            patient_data={"sbp": 115.0, "hr": 70, "potassium": 6.1, "egfr": 18.0},
        )
        recs = " ".join(result.recommendations)
        assert "INITIATE MRA: Start spironolactone" not in recs
        assert "K+ is 6.1 mEq/L" in recs


# =====================================================================
# B-2 — dose_mg vs dose: a patient at target was told to titrate DOWN
# =====================================================================

class TestB2DoseKeyMismatch:
    """`_assess_current_gdmt` read `med.get("dose")` while the documented
    request schema spells it `{"name": ..., "dose_mg": 25, "frequency": "BID"}`,
    so every API-driven medication read as sub-target: a patient already at the
    carvedilol target of 25 mg BID got
    `TITRATE carvedilol: None -> 25 mg BID ... Next step: 6.25 mg BID`, the
    second rung of the ladder.
    """

    def test_dose_mg_is_read_as_a_dose(self):
        assert _dose_from_med(
            {"name": "carvedilol", "dose_mg": 25, "frequency": "BID"}
        ) == "25 mg BID"
        assert _dose_from_med({"name": "carvedilol", "dose": "25 mg BID"}) == "25 mg BID"
        assert _dose_from_med({"name": "carvedilol"}) == ""

    def test_patient_at_target_is_not_told_to_titrate_down(self):
        result = GDMTOptimizer().optimize(
            lvef=30.0, nyha_class="III",
            current_medications=[
                {"name": "carvedilol", "dose_mg": 25, "frequency": "BID"},
            ],
            patient_data={"sbp": 120.0, "hr": 65, "potassium": 4.2, "egfr": 60.0},
        )
        text = " ".join(result.recommendations + result.next_steps)
        assert "Next step: 6.25 mg BID" not in text
        assert "TITRATE carvedilol" not in text

    def test_both_dose_spellings_give_the_same_status(self):
        optimizer = GDMTOptimizer()
        by_dose_mg = optimizer._assess_current_gdmt(
            [{"name": "carvedilol", "dose_mg": 25, "frequency": "BID"}]
        )
        by_dose = optimizer._assess_current_gdmt(
            [{"name": "carvedilol", "dose": "25 mg BID"}]
        )
        assert by_dose_mg[GDMTPillar.BETA_BLOCKER].status == GDMTStatus.AT_TARGET
        assert by_dose[GDMTPillar.BETA_BLOCKER].status == GDMTStatus.AT_TARGET

    def test_unrecorded_dose_is_not_placed_at_the_bottom_of_the_ladder(self):
        optimizer = GDMTOptimizer()
        meds = optimizer._assess_current_gdmt([{"name": "carvedilol"}])
        plan = optimizer._generate_titration_plan(
            meds[GDMTPillar.BETA_BLOCKER], {"sbp": 120.0, "hr": 65},
        )
        assert plan["dose_not_recorded"] is True
        assert plan["remaining_steps"] == []
        assert "cannot place" in plan["plan"]


# =====================================================================
# B-8 — ARNI starting dose chosen from systolic BP
# =====================================================================

class TestB8ARNIStartingDoseFollowsTheLabel:
    """`starting = "49/51 mg BID" if (sbp and sbp >= 110) else "24/26 mg BID"`
    used a criterion the ENTRESTO label does not contain in section 2 (the
    systolic cut-off appears only in the PARADIGM-HF eligibility criteria in
    section 14.1), and doubled the label dose for RAASi-naive patients and for
    eGFR <30 -- the two groups the label halves it for (sections 2.6 and 2.7).

    Primary source: ENTRESTO FDA-approved labeling, DailyMed setid
    000dc81d-ab91-450c-8eae-8eb74e72296f, sections 2.2, 2.6, 2.7, 2.8.
    """

    def _recs(self, **patient):
        base = {"sbp": 115.0, "hr": 70, "potassium": 4.0, "egfr": 60.0}
        base.update(patient)
        result = GDMTOptimizer().optimize(
            lvef=25.0, nyha_class="III", current_medications=[], patient_data=base,
        )
        return " ".join(result.recommendations)

    def test_raasi_naive_patient_starts_at_half_dose(self):
        recs = self._recs(sbp=115.0, egfr=60.0)
        assert "sacubitril/valsartan 24/26 mg BID" in recs
        assert "sacubitril/valsartan 49/51 mg BID" not in recs
        assert "not currently taking an ACE inhibitor or ARB (label 2.6)" in recs

    def test_egfr_below_30_starts_at_half_dose(self):
        recs = self._recs(sbp=115.0, egfr=22.0)
        assert "sacubitril/valsartan 24/26 mg BID" in recs
        assert "sacubitril/valsartan 49/51 mg BID" not in recs
        assert "label 2.7" in recs

    def test_systolic_bp_alone_does_not_double_the_dose(self):
        """SBP is not a dosing criterion in the label, so SBP 140 and SBP 95
        must give the same starting dose."""
        assert "24/26 mg BID" in self._recs(sbp=140.0, egfr=60.0)
        assert "24/26 mg BID" in self._recs(sbp=95.0, egfr=60.0)

    def test_patient_on_a_full_dose_arb_starts_at_the_usual_dose(self):
        """A patient already on a full-dose ARB meets no label halving
        criterion, so the switch to ARNI starts at the usual 49/51 mg BID.  The
        switch message used to carry no dose at all."""
        result = GDMTOptimizer().optimize(
            lvef=25.0, nyha_class="III",
            current_medications=[{"name": "valsartan", "dose": "160 mg BID"}],
            patient_data={"sbp": 118.0, "hr": 70, "potassium": 4.0, "egfr": 60.0},
        )
        recs = " ".join(result.recommendations)
        assert "Start at 49/51 mg BID" in recs
        assert "usual starting dose per ENTRESTO label 2.2" in recs

    def test_patient_on_a_low_dose_arb_starts_at_half_dose(self):
        result = GDMTOptimizer().optimize(
            lvef=25.0, nyha_class="III",
            current_medications=[{"name": "valsartan", "dose": "40 mg BID"}],
            patient_data={"sbp": 118.0, "hr": 70, "potassium": 4.0, "egfr": 60.0},
        )
        recs = " ".join(result.recommendations)
        assert "Start at 24/26 mg BID" in recs
        assert "previously taking a low dose" in recs


# =====================================================================
# B-9 / B-10 — device eligibility was assumed
# =====================================================================

class TestB9ICDEligibilityIsNotAssumed:
    """`icd_eligible = True` was the starting state, the waiting periods only
    created a barrier `if ... is not None`, and `expected_survival_years`
    defaulted to 2, so the >1-year requirement could never fail on absent data.
    A patient 10 days after an anterior MI with EF 30% was told they met
    Class I primary-prevention ICD criteria -- a device, and its lifetime of
    infection and lead risk, in a patient whose EF commonly recovers in that
    window, which is why the waiting period exists.
    """

    def _recs(self, **patient):
        base = {"sbp": 115.0, "hr": 70, "potassium": 4.0, "egfr": 60.0}
        base.update(patient)
        result = GDMTOptimizer().optimize(
            lvef=30.0, nyha_class="III", current_medications=[], patient_data=base,
        )
        return " ".join(result.recommendations)

    def test_nothing_recorded_means_eligibility_is_not_confirmed(self):
        recs = self._recs()
        assert "Patient meets primary prevention ICD criteria" not in recs
        assert "DEVICE THERAPY - ICD: INDETERMINATE" in recs
        assert "days since most recent MI not recorded" in recs
        assert "days on optimal GDMT not recorded" in recs
        assert "expected survival not recorded" in recs

    def test_expected_survival_no_longer_defaults_to_two_years(self):
        recs = self._recs(days_post_mi=120, days_on_gdmt=120)
        assert "Patient meets primary prevention ICD criteria" not in recs
        assert "expected survival not recorded" in recs

    def test_all_criteria_documented_gives_the_class_i_statement(self):
        recs = self._recs(
            days_post_mi=120, days_on_gdmt=120, expected_survival_years=5,
        )
        assert "DEVICE THERAPY - ICD (Class I, LOE A)" in recs
        assert "Patient meets primary prevention ICD criteria" in recs

    def test_ten_days_post_mi_is_a_barrier(self):
        recs = self._recs(
            days_post_mi=10, days_on_gdmt=120, expected_survival_years=5,
        )
        assert "Patient meets primary prevention ICD criteria" not in recs
        assert "Only 10 days post-MI" in recs


class TestB10CRTReadsTheRhythm:
    """`_assess_device_therapy` gated CRT on LVEF, NYHA and QRS only -- `rhythm`
    was read in the same file for ivabradine but never for CRT -- so a patient
    in permanent AF with LVEF 30%, QRS 160 ms and LBBB was told
    "DEVICE THERAPY - CRT (Class I, LOE A): Patient meets CRT criteria".  That
    class comes from CARE-HF, a sinus-rhythm recommendation, and the module's
    own criteria table reads "sinus rhythm preferred (AF with rate control may
    benefit)".
    """

    def _recs(self, **patient):
        base = {
            "sbp": 115.0, "hr": 70, "potassium": 4.0, "egfr": 60.0,
            "qrs_duration": 160, "qrs_morphology": "LBBB",
        }
        base.update(patient)
        result = GDMTOptimizer().optimize(
            lvef=30.0, nyha_class="III", current_medications=[], patient_data=base,
        )
        return " ".join(result.recommendations)

    def test_permanent_af_does_not_get_a_class_i_crt_statement(self):
        recs = self._recs(rhythm="permanent atrial fibrillation")
        assert "DEVICE THERAPY - CRT (Class I, LOE A)" not in recs
        assert "Patient meets CRT criteria" not in recs
        assert "DEVICE THERAPY - CRT: INDETERMINATE" in recs
        assert "not sinus rhythm" in recs

    def test_unrecorded_rhythm_does_not_get_a_class_i_crt_statement(self):
        recs = self._recs()
        assert "Patient meets CRT criteria" not in recs
        assert "the rhythm was not recorded" in recs

    def test_sinus_rhythm_still_gets_the_class_i_statement(self):
        recs = self._recs(rhythm="sinus")
        assert "DEVICE THERAPY - CRT (Class I, LOE A)" in recs
        assert "Patient meets CRT criteria" in recs

    def test_crt_d_is_not_claimed_when_icd_eligibility_is_unconfirmed(self):
        recs = self._recs(rhythm="sinus")
        assert "DEVICE THERAPY - CRT-D" not in recs


# =====================================================================
# B-11 — the optimizer's copies of the shared thresholds
# =====================================================================

class TestB11OptimizerUsesTheSharedThresholds:
    """The HFrEF LVEF boundary, the CRT Class 2a QRS threshold and the SGLT2i
    class outside HFrEF were restated in four places with three disagreements,
    so the same patient got different advice from different endpoints.  They now
    come from src/hf_thresholds.py; see that module for the sources.
    """

    def test_ef_boundaries_come_from_the_shared_module(self):
        assert GDMTOptimizer._HFREF_EF_CUTOFF == HFREF_LVEF_CUTOFF == 40.0
        assert GDMTOptimizer._HFPEF_EF_CUTOFF == HFPEF_LVEF_CUTOFF == 50.0
        assert GDMTOptimizer()._classify_ef(40.0, {}) == (
            EjectionFractionCategory.HFrEF
        )

    def test_lvef_49_5_is_hfmref_not_hfpef(self):
        """Comparing against the integer 49 put a 49.5% ventricle into HFpEF,
        whose banner then read "HFpEF MANAGEMENT (LVEF >=50%)"."""
        result = GDMTOptimizer().optimize(
            lvef=49.5, nyha_class="II", current_medications=[],
            patient_data={"sbp": 120.0, "hr": 70, "potassium": 4.0, "egfr": 60.0},
        )
        assert result.ef_category == EjectionFractionCategory.HFmrEF
        assert "HFpEF MANAGEMENT (LVEF >=50%)" not in " ".join(result.recommendations)

    def test_crt_qrs_thresholds_come_from_the_shared_module(self):
        assert GDMTOptimizer._CRT_QRS_CLASS_IIA == CRT_QRS_CLASS_IIA_MS
        assert GDMTOptimizer._CRT_QRS_CLASS_I == CRT_QRS_CLASS_I_MS

    def test_sglt2i_outside_hfref_is_class_2a(self):
        assert SGLT2I_CLASS_HFMREF == GuidelineClass.CLASS_IIA
        assert SGLT2I_CLASS_HFPEF == GuidelineClass.CLASS_IIA
        hfpef = GDMTOptimizer().optimize(
            lvef=60.0, nyha_class="II", current_medications=[],
            patient_data={"sbp": 120.0, "hr": 70, "potassium": 4.0, "egfr": 60.0},
        )
        text = " ".join(hfpef.recommendations)
        assert "SGLT2i is Class 2a" in text
        assert "SGLT2i is Class I" not in text


# =====================================================================
# B-14 — MRA prose contradicted the threshold the code tests
# =====================================================================

def test_b14_mra_potassium_prose_matches_the_code():
    """The guard is `k > 5.0` and the message said `(>=5.0)`; at K+ exactly 5.0
    the optimizer initiates the MRA, so the prose contradicted the behaviour."""
    result = GDMTOptimizer().optimize(
        lvef=25.0, nyha_class="III", current_medications=[],
        patient_data={"sbp": 115.0, "hr": 70, "potassium": 5.4, "egfr": 60.0},
    )
    text = " ".join(result.recommendations)
    assert "(>5.0)" in text
    assert "(>=5.0)" not in text

    allowed = GDMTOptimizer().optimize(
        lvef=25.0, nyha_class="III", current_medications=[],
        patient_data={"sbp": 115.0, "hr": 70, "potassium": 5.0, "egfr": 60.0},
    )
    assert "INITIATE MRA" in " ".join(allowed.recommendations)
