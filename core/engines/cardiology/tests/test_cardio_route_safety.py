"""Regression tests for the cardiology API routes' clinical-safety defects.

Found by the pre-release audit (identifiers B-1, B-6, B-7, B-11, A-2, A-3 and
A-4).  These are the faults at the HTTP boundary: a nested request field passed
whole, a pillar reported at target without reading the dose, a swallowed engine
failure, and risk-calculator key names the calculators do not read.
"""
from unittest.mock import patch

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from api.routes.cardio_clinical import router
from src.gdmt_optimizer import GDMTOptimizer
from src.hf_thresholds import HFREF_LVEF_CUTOFF
from src.risk_calculators import RiskCalculatorEngine


class _Lock:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _client(optimizer=None):
    app = FastAPI()
    app.include_router(router)
    app.state.metrics = {
        "gdmt_requests_total": 0, "risk_calc_requests_total": 0,
        "requests_total": 0, "workflow_requests_total": 0,
    }
    app.state.metrics_lock = _Lock()
    app.state.engine = None
    app.state.gdmt_optimizer = optimizer
    app.state.risk_calculators = {"engine": RiskCalculatorEngine()}
    return TestClient(app)


def _rec_text(body: dict) -> str:
    return " ".join(
        r["recommendation"] if isinstance(r, dict) else str(r)
        for r in body["recommendations"]
    )


# =====================================================================
# B-1 — /gdmt/optimize passed the wrong dict, so every lab, vital and
# contraindication the client sent was invisible to the optimizer
# =====================================================================

B1_PAYLOAD = {
    "lvef": 28.0,
    "nyha_class": "III",
    "current_medications": [{"name": "carvedilol", "dose_mg": 25, "frequency": "BID"}],
    "patient_data": {
        "potassium": 6.1,
        "egfr": 18,
        "sbp": 82,
        "hr": 48,
        "contraindications": ["history of angioedema with ACEi"],
        "qrs_duration": 168,
        "qrs_morphology": "LBBB",
        "days_post_mi": 5,
        "days_on_gdmt": 3,
        "rhythm": "sinus",
    },
}


class TestB1NestedPatientDataReachesTheOptimizer:
    """`patient_data` is a NESTED field on `GDMTRequest`, but the route passed
    the whole `request.model_dump()`, so every key the optimizer reads --
    potassium, egfr, sbp, hr, contraindications, qrs_duration, qrs_morphology,
    days_post_mi, days_on_gdmt, rhythm -- sat one level down and resolved to
    None.  Every guard in the optimizer is `if value is not None and ...`, so
    absent read as benign, and one response carried: an MRA order at K+ 6.1 with
    eGFR 18 (the combination that causes fatal hyperkalaemia), an ARNI order
    after ACEi angioedema (a boxed contraindication) at SBP 82, a beta-blocker
    titration at HR 48, a Class I ICD assertion 5 days post-MI, and no CRT
    assessment at all for QRS 168 ms with LBBB.
    """

    @pytest.fixture
    def body(self):
        response = _client(GDMTOptimizer()).post(
            "/v1/cardio/gdmt/optimize", json=B1_PAYLOAD,
        )
        assert response.status_code == 200
        return response.json()

    def test_optimizer_actually_ran(self, body):
        assert body["status"] == "optimizer"
        assert body["error"] is None

    def test_no_mra_order_at_potassium_6_1(self, body):
        text = _rec_text(body)
        assert "INITIATE MRA: Start spironolactone" not in text
        assert "MRA: K+ is 6.1 mEq/L" in text

    def test_no_arni_order_after_acei_angioedema(self, body):
        text = _rec_text(body)
        assert "INITIATE ARNI: Start sacubitril/valsartan" not in text
        assert "History of angioedema" in text

    def test_no_beta_blocker_titration_down_to_6_25_at_hr_48(self, body):
        text = _rec_text(body) + body["summary"]
        assert "Next step: 6.25 mg BID" not in text

    def test_icd_not_asserted_five_days_post_mi(self, body):
        text = _rec_text(body)
        assert "Patient meets primary prevention ICD criteria" not in text
        assert "Only 5 days post-MI" in text
        assert "Only 3 days on optimal GDMT" in text

    def test_crt_is_assessed_for_qrs_168_lbbb(self, body):
        text = _rec_text(body)
        assert "DEVICE THERAPY - CRT (Class I, LOE A)" in text


def test_b1_optimizer_output_is_not_silently_dropped_by_the_response_model():
    """`GDMTResponse.recommendations` is `List[GDMTRecommendation]` while the
    optimizer returns sentences, so handing the strings over raised a pydantic
    ValidationError that the route's bare `except Exception` swallowed: EVERY
    /gdmt/optimize response came from the simplified built-in fallback and
    nothing said so.
    """
    body = _client(GDMTOptimizer()).post(
        "/v1/cardio/gdmt/optimize", json=B1_PAYLOAD,
    ).json()
    assert body["status"] == "optimizer"
    assert "HFrEF MANAGEMENT" in _rec_text(body)
    # The four-pillar tile came back empty on this path.
    assert body["four_pillars_status"]


# =====================================================================
# B-6 — "AT TARGET" without examining the dose
# =====================================================================

class TestB6PillarStatusReadsTheDose:
    """`at_target = True  # Simplified; would check actual dose` reported
    "Beta-blocker: AT TARGET" for a patient on carvedilol 3.125 mg BID, an
    eighth of the 25 mg BID mortality-reducing dose, and emitted no up-titration
    recommendation -- hiding under-titration, the default failure mode of
    heart-failure care, from the tile meant to catch it.
    """

    def test_one_eighth_of_target_is_not_at_target(self):
        body = _client().post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
            "current_medications": [
                {"name": "carvedilol", "dose_mg": 3.125, "frequency": "BID"},
            ],
        }).json()
        assert body["four_pillars_status"]["Beta-blocker"] == "SUBOPTIMAL"
        assert "UP-TITRATE Beta-blocker" in _rec_text(body)

    def test_target_dose_is_reported_at_target(self):
        body = _client().post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
            "current_medications": [
                {"name": "carvedilol", "dose_mg": 25, "frequency": "BID"},
            ],
        }).json()
        assert body["four_pillars_status"]["Beta-blocker"] == "AT TARGET"

    def test_unrecorded_dose_is_neither_at_target_nor_suboptimal(self):
        body = _client().post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
            "current_medications": [{"name": "carvedilol"}],
        }).json()
        assert body["four_pillars_status"]["Beta-blocker"] == "DOSE NOT RECORDED"
        assert "RECORD Beta-blocker DOSE" in _rec_text(body)


# =====================================================================
# B-7 — a swallowed optimizer failure returned a complete-looking plan
# =====================================================================

class TestB7OptimizerFailureIsReported:
    """`except Exception as exc: logger.warning(...)` fell through to the
    built-in logic and returned HTTP 200 with no status and no error field, so
    the caller could not tell the guideline engine had been bypassed and that
    what came back was the route's simplified fallback.
    """

    def test_optimizer_exception_is_reported_as_degraded(self):
        class _Boom:
            def optimize(self, **_kwargs):
                raise RuntimeError("milvus down")

        response = _client(_Boom()).post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
        })
        body = response.json()
        assert response.status_code == 200
        assert body["status"] == "degraded"
        assert body["status"] != "optimizer"
        assert "milvus down" in body["error"]
        assert "DEGRADED" in body["summary"]

    def test_absent_optimizer_says_builtin(self):
        body = _client(None).post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
        }).json()
        assert body["status"] == "builtin"
        assert body["error"] is None

    def test_successful_optimizer_says_so(self):
        body = _client(GDMTOptimizer()).post("/v1/cardio/gdmt/optimize", json={
            "lvef": 25.0, "nyha_class": "II",
        }).json()
        assert body["status"] == "optimizer"
        assert body["error"] is None


# =====================================================================
# B-11 — the route disagreed with the other three surfaces
# =====================================================================

def test_b11_lvef_exactly_40_is_hfref_on_the_route():
    """`if request.lvef < 40` in the route against `<= 40` in GDMTOptimizer,
    HeartFailureWorkflow and risk_calculators: an LVEF of exactly 40.0% dropped
    from four Class I "INITIATE" recommendations to "Consider"."""
    body = _client().post("/v1/cardio/gdmt/optimize", json={
        "lvef": 40.0, "nyha_class": "III",
    }).json()
    assert body["hf_phenotype"] == "HFrEF"
    assert "Consider" not in str(list(body["four_pillars_status"].values()))
    assert HFREF_LVEF_CUTOFF == 40.0


def test_b11_sglt2i_class_in_hfmref_and_hfpef_is_2a_on_the_route():
    for lvef in (45.0, 60.0):
        body = _client().post("/v1/cardio/gdmt/optimize", json={
            "lvef": lvef, "nyha_class": "II",
        }).json()
        assert "Class 2a" in _rec_text(body)


# =====================================================================
# A-2 / A-3 / A-4 — risk-calculator route wiring
# =====================================================================

class TestA2MAGGICDurationReachesTheCalculator:
    """The route sent `extra={"hf_duration_18m": ...}` while `calculate_maggic`
    reads `extra["hf_duration_months"]`, so the duration fell back to 0 on every
    API call: the "HF diagnosed >18 months ago" checkbox the form collects had
    no effect on the score at all, and nothing said the input was discarded.

    These tests are deliberately direction-agnostic.  The sign of the duration
    term is a separate defect in `src/risk_calculators.py` (audit finding A-1),
    recorded in CLINICAL_REVIEW.md section 12 with the Pocock 2013 Table 2
    citation; it is NOT asserted here.
    """

    BASE = {
        "age": 70, "sex": "male", "lvef": 30.0, "nyha_class": 3,
        "systolic_bp": 120.0, "bmi": 26.0, "creatinine": 1.0,
        "beta_blocker": True, "acei_arb": True,
    }

    def test_the_duration_flag_changes_the_score(self):
        client = _client()
        chronic = client.post(
            "/v1/cardio/risk/maggic", json={**self.BASE, "hf_duration_18m": True},
        ).json()["score"]
        new_diagnosis = client.post(
            "/v1/cardio/risk/maggic", json={**self.BASE, "hf_duration_18m": False},
        ).json()["score"]
        assert chronic != new_diagnosis, (
            "the HF-duration input must reach the calculator; it was silently "
            "discarded by a key-name mismatch"
        )

    def test_the_route_sends_the_key_the_calculator_reads(self):
        client = _client()
        with patch.object(
            RiskCalculatorEngine, "calculate", autospec=True,
        ) as calculate:
            calculate.return_value = type(
                "_R", (), {
                    "score_value": 20.0, "risk_category": "Intermediate",
                    "interpretation": "", "recommendations": [],
                    "guideline_reference": "",
                },
            )()
            client.post(
                "/v1/cardio/risk/maggic", json={**self.BASE, "hf_duration_18m": True},
            )
        extra = calculate.call_args.kwargs["extra"]
        assert "hf_duration_months" in extra
        assert extra["hf_duration_months"] >= 18
        assert "hf_duration_18m" not in extra


class TestA3A4EuroSCOREIIRouteWiring:
    """Six wiring faults in one route: four key names the calculator does not
    read (renal, critical preoperative state, thoracic aorta, NYHA),
    `weight_of_procedure` never sent at all, and creatinine clearance passed
    through `RiskScoreInput.creatinine` (bounded 0.1-20.0 mg/dL), which raised a
    pydantic ValidationError that the bare `except Exception` turned into HTTP
    500 for every patient with a CrCl above 20 mL/min -- including the form's own
    default of 65.  It worked only in the dialysis range, where it then scored
    the renal term as if kidney function were normal.

    The 78-year-old referred for valve surgery: 3.21% ("Moderate operative
    risk") as wired, against 30.93% ("Very high ... goals-of-care discussion")
    under the key names the calculator reads -- a tenfold understatement of the
    number a surgeon quotes during consent.  With age and sex only, 1.03%
    ("Standard perioperative management appropriate") against her real 42.35%.
    """

    COMPLETE = {
        "age": 78, "sex": "female", "creatinine_clearance": 18.0,
        "lvef": 25.0, "nyha_class": 4, "critical_preop_state": True,
        "pulmonary_hypertension": "severe", "thoracic_aorta": True,
        "procedure_weight": "two_procedures", "urgency": "elective",
    }

    def test_a_normal_creatinine_clearance_no_longer_500s(self):
        response = _client().post("/v1/cardio/risk/euroscore", json={
            "age": 70, "sex": "male", "creatinine_clearance": 65.0,
            "lvef": 55.0, "nyha_class": 2, "procedure_weight": "isolated_cabg",
        })
        assert response.status_code == 200
        assert response.json()["score"] > 0

    def test_every_model_term_reaches_the_calculator(self):
        body = _client().post(
            "/v1/cardio/risk/euroscore", json=self.COMPLETE,
        ).json()
        assert body["score"] > 20.0, (
            "this patient scored 3.21% while five model terms were dropped"
        )
        assert body["risk_category"] in ("High", "Very High")
        interpretation = body["interpretation"]
        for term in ("renal", "critical", "NYHA"):
            assert term in interpretation

    def test_dropping_critical_preop_state_changes_the_score(self):
        client = _client()
        with_it = client.post(
            "/v1/cardio/risk/euroscore", json=self.COMPLETE,
        ).json()["score"]
        without = client.post("/v1/cardio/risk/euroscore", json={
            **self.COMPLETE, "critical_preop_state": False,
        }).json()["score"]
        assert with_it > without

    def test_dropping_the_procedure_weight_changes_the_score(self):
        client = _client()
        two = client.post(
            "/v1/cardio/risk/euroscore", json=self.COMPLETE,
        ).json()["score"]
        isolated = client.post("/v1/cardio/risk/euroscore", json={
            **self.COMPLETE, "procedure_weight": "isolated_cabg",
        }).json()["score"]
        assert two > isolated

    def test_thoracic_aorta_surgery_changes_the_score(self):
        client = _client()
        with_it = client.post(
            "/v1/cardio/risk/euroscore", json=self.COMPLETE,
        ).json()["score"]
        without = client.post("/v1/cardio/risk/euroscore", json={
            **self.COMPLETE, "thoracic_aorta": False,
        }).json()["score"]
        assert with_it > without

    def test_age_and_sex_alone_is_refused_not_scored_as_healthy(self):
        response = _client().post(
            "/v1/cardio/risk/euroscore", json={"age": 78, "sex": "female"},
        )
        assert response.status_code == 422
        detail = response.json()["detail"]
        assert "indeterminate" in detail.lower()
        for field in ("creatinine_clearance", "lvef", "nyha_class",
                      "procedure_weight"):
            assert field in detail

    @pytest.mark.parametrize(
        "omit", ["lvef", "nyha_class", "procedure_weight", "creatinine_clearance"],
    )
    def test_each_patient_factor_is_required(self, omit):
        payload = {k: v for k, v in self.COMPLETE.items() if k != omit}
        response = _client().post("/v1/cardio/risk/euroscore", json=payload)
        assert response.status_code == 422
        assert omit.split("_")[0] in response.json()["detail"]

    def test_dialysis_patients_need_no_creatinine_clearance(self):
        payload = {k: v for k, v in self.COMPLETE.items()
                   if k != "creatinine_clearance"}
        payload["on_dialysis"] = True
        response = _client().post("/v1/cardio/risk/euroscore", json=payload)
        assert response.status_code == 200
