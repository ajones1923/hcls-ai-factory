"""Tests for the Cardiology Intelligence Agent API routes.

Validates route handler function existence, request model validation
(field ranges, patterns), router configuration (prefix, tags), and
reference catalogue endpoints.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import unittest

from pydantic import ValidationError

from api.routes.cardio_clinical import (
    ASCVDRequest,
    CHA2DS2VAScRequest,
    EuroSCORERequest,
    GDMTRequest,
    HASBLEDRequest,
    HEARTRequest,
    MAGGICRequest,
    QueryRequest,
    RiskResult,
    SearchRequest,
    router,
    # Route handler functions
    cad_assessment,
    calculate_ascvd,
    calculate_cha2ds2vasc,
    calculate_euroscore,
    calculate_hasbled,
    calculate_heart,
    calculate_maggic,
    cardiac_mri_workflow,
    cardio_oncology_workflow,
    cardio_query,
    cardio_search,
    heart_failure_workflow,
    arrhythmia_workflow,
    list_biomarkers,
    list_conditions,
    list_drug_classes,
    list_genes,
    list_guidelines,
    optimize_gdmt,
    prevention_workflow,
    stress_test_workflow,
    valvular_workflow,
)


# =====================================================================
# Route handler function existence
# =====================================================================


class TestRouteHandlersExist(unittest.TestCase):
    """Verify all route handler functions are defined and callable."""

    def test_cardio_query_exists(self):
        self.assertTrue(callable(cardio_query))

    def test_cardio_search_exists(self):
        self.assertTrue(callable(cardio_search))

    def test_calculate_ascvd_exists(self):
        self.assertTrue(callable(calculate_ascvd))

    def test_calculate_heart_exists(self):
        self.assertTrue(callable(calculate_heart))

    def test_calculate_cha2ds2vasc_exists(self):
        self.assertTrue(callable(calculate_cha2ds2vasc))

    def test_calculate_hasbled_exists(self):
        self.assertTrue(callable(calculate_hasbled))

    def test_calculate_maggic_exists(self):
        self.assertTrue(callable(calculate_maggic))

    def test_calculate_euroscore_exists(self):
        self.assertTrue(callable(calculate_euroscore))

    def test_optimize_gdmt_exists(self):
        self.assertTrue(callable(optimize_gdmt))

    def test_cad_assessment_exists(self):
        self.assertTrue(callable(cad_assessment))

    def test_heart_failure_workflow_exists(self):
        self.assertTrue(callable(heart_failure_workflow))

    def test_valvular_workflow_exists(self):
        self.assertTrue(callable(valvular_workflow))

    def test_arrhythmia_workflow_exists(self):
        self.assertTrue(callable(arrhythmia_workflow))

    def test_cardiac_mri_workflow_exists(self):
        self.assertTrue(callable(cardiac_mri_workflow))

    def test_stress_test_workflow_exists(self):
        self.assertTrue(callable(stress_test_workflow))

    def test_prevention_workflow_exists(self):
        self.assertTrue(callable(prevention_workflow))

    def test_cardio_oncology_workflow_exists(self):
        self.assertTrue(callable(cardio_oncology_workflow))

    def test_list_conditions_exists(self):
        self.assertTrue(callable(list_conditions))

    def test_list_biomarkers_exists(self):
        self.assertTrue(callable(list_biomarkers))

    def test_list_drug_classes_exists(self):
        self.assertTrue(callable(list_drug_classes))

    def test_list_genes_exists(self):
        self.assertTrue(callable(list_genes))

    def test_list_guidelines_exists(self):
        self.assertTrue(callable(list_guidelines))


# =====================================================================
# Router configuration
# =====================================================================


class TestRouterConfig(unittest.TestCase):
    """Tests for the APIRouter prefix and tags."""

    def test_router_prefix(self):
        self.assertEqual(router.prefix, "/v1/cardio")

    def test_router_has_tags(self):
        self.assertIn("cardiology", router.tags)


# =====================================================================
# ASCVDRequest validation
# =====================================================================


class TestASCVDRequestValidation(unittest.TestCase):
    """Tests for ASCVDRequest Pydantic model validation."""

    def _valid_params(self):
        return {
            "age": 55,
            "sex": "male",
            "race": "white",
            "total_cholesterol": 200.0,
            "hdl_cholesterol": 50.0,
            "systolic_bp": 130.0,
            "bp_treatment": False,
            "diabetes": False,
            "smoker": False,
        }

    def test_valid_request(self):
        req = ASCVDRequest(**self._valid_params())
        self.assertEqual(req.age, 55)

    def test_age_min_40(self):
        params = self._valid_params()
        params["age"] = 39
        with self.assertRaises(ValidationError):
            ASCVDRequest(**params)

    def test_age_max_79(self):
        params = self._valid_params()
        params["age"] = 80
        with self.assertRaises(ValidationError):
            ASCVDRequest(**params)

    def test_age_at_lower_bound(self):
        params = self._valid_params()
        params["age"] = 40
        req = ASCVDRequest(**params)
        self.assertEqual(req.age, 40)

    def test_age_at_upper_bound(self):
        params = self._valid_params()
        params["age"] = 79
        req = ASCVDRequest(**params)
        self.assertEqual(req.age, 79)

    def test_sex_must_be_male_or_female(self):
        params = self._valid_params()
        params["sex"] = "other"
        with self.assertRaises(ValidationError):
            ASCVDRequest(**params)

    def test_race_valid_values(self):
        for race in ["white", "african_american", "other"]:
            params = self._valid_params()
            params["race"] = race
            req = ASCVDRequest(**params)
            self.assertEqual(req.race, race)

    def test_total_cholesterol_min(self):
        params = self._valid_params()
        params["total_cholesterol"] = 50.0
        with self.assertRaises(ValidationError):
            ASCVDRequest(**params)

    def test_total_cholesterol_max(self):
        params = self._valid_params()
        params["total_cholesterol"] = 500.0
        with self.assertRaises(ValidationError):
            ASCVDRequest(**params)


# =====================================================================
# HEARTRequest validation
# =====================================================================


class TestHEARTRequestValidation(unittest.TestCase):
    """Tests for HEARTRequest Pydantic model validation."""

    def _valid_params(self):
        return {
            "history_score": 1,
            "ecg_score": 1,
            "age": 55,
            "risk_factors": 2,
            "troponin_score": 0,
        }

    def test_valid_request(self):
        req = HEARTRequest(**self._valid_params())
        self.assertEqual(req.history_score, 1)

    def test_history_score_min_0(self):
        params = self._valid_params()
        params["history_score"] = -1
        with self.assertRaises(ValidationError):
            HEARTRequest(**params)

    def test_history_score_max_2(self):
        params = self._valid_params()
        params["history_score"] = 3
        with self.assertRaises(ValidationError):
            HEARTRequest(**params)

    def test_ecg_score_range(self):
        for val in [0, 1, 2]:
            params = self._valid_params()
            params["ecg_score"] = val
            req = HEARTRequest(**params)
            self.assertEqual(req.ecg_score, val)

    def test_ecg_score_out_of_range(self):
        params = self._valid_params()
        params["ecg_score"] = 3
        with self.assertRaises(ValidationError):
            HEARTRequest(**params)

    def test_troponin_score_range(self):
        for val in [0, 1, 2]:
            params = self._valid_params()
            params["troponin_score"] = val
            req = HEARTRequest(**params)
            self.assertEqual(req.troponin_score, val)

    def test_troponin_score_out_of_range(self):
        params = self._valid_params()
        params["troponin_score"] = 3
        with self.assertRaises(ValidationError):
            HEARTRequest(**params)


# =====================================================================
# CHA2DS2VAScRequest validation
# =====================================================================


class TestCHA2DS2VAScRequestValidation(unittest.TestCase):
    """Tests for CHA2DS2VAScRequest Pydantic model validation."""

    def _valid_params(self):
        return {
            "chf": False,
            "hypertension": True,
            "age": 70,
            "diabetes": False,
            "stroke_tia": False,
            "vascular_disease": False,
            "sex": "male",
        }

    def test_valid_request(self):
        req = CHA2DS2VAScRequest(**self._valid_params())
        self.assertEqual(req.age, 70)

    def test_sex_must_be_valid(self):
        params = self._valid_params()
        params["sex"] = "nonbinary"
        with self.assertRaises(ValidationError):
            CHA2DS2VAScRequest(**params)

    def test_boolean_fields_default_false(self):
        req = CHA2DS2VAScRequest(age=65, sex="male")
        self.assertFalse(req.chf)
        self.assertFalse(req.hypertension)
        self.assertFalse(req.diabetes)
        self.assertFalse(req.stroke_tia)
        self.assertFalse(req.vascular_disease)


# =====================================================================
# HASBLEDRequest validation
# =====================================================================


class TestHASBLEDRequestValidation(unittest.TestCase):
    """Tests for HASBLEDRequest Pydantic model validation."""

    def test_valid_request(self):
        req = HASBLEDRequest()
        self.assertFalse(req.hypertension_uncontrolled)

    def test_drugs_alcohol_min_0(self):
        req = HASBLEDRequest(drugs_alcohol=0)
        self.assertEqual(req.drugs_alcohol, 0)

    def test_drugs_alcohol_max_2(self):
        req = HASBLEDRequest(drugs_alcohol=2)
        self.assertEqual(req.drugs_alcohol, 2)

    def test_drugs_alcohol_over_max(self):
        with self.assertRaises(ValidationError):
            HASBLEDRequest(drugs_alcohol=3)

    def test_all_booleans_default_false(self):
        req = HASBLEDRequest()
        self.assertFalse(req.renal_disease)
        self.assertFalse(req.liver_disease)
        self.assertFalse(req.stroke_history)
        self.assertFalse(req.bleeding_history)
        self.assertFalse(req.labile_inr)
        self.assertFalse(req.age_over_65)


# =====================================================================
# QueryRequest validation
# =====================================================================


class TestQueryRequestValidation(unittest.TestCase):
    """Tests for QueryRequest Pydantic model validation."""

    def test_valid_request(self):
        req = QueryRequest(question="What is heart failure?")
        self.assertEqual(req.question, "What is heart failure?")

    def test_question_min_length(self):
        with self.assertRaises(ValidationError):
            QueryRequest(question="Hi")

    def test_default_top_k(self):
        req = QueryRequest(question="Test question here")
        self.assertEqual(req.top_k, 5)

    def test_default_include_guidelines(self):
        req = QueryRequest(question="Test question here")
        self.assertTrue(req.include_guidelines)

    def test_workflow_type_optional(self):
        req = QueryRequest(question="Test question here")
        self.assertIsNone(req.workflow_type)

    def test_patient_context_optional(self):
        req = QueryRequest(question="Test question here")
        self.assertIsNone(req.patient_context)


# =====================================================================
# SearchRequest validation
# =====================================================================


class TestSearchRequestValidation(unittest.TestCase):
    """Tests for SearchRequest Pydantic model validation."""

    def test_valid_request(self):
        req = SearchRequest(question="SGLT2 inhibitors")
        self.assertEqual(req.question, "SGLT2 inhibitors")

    def test_question_min_length(self):
        with self.assertRaises(ValidationError):
            SearchRequest(question="ab")

    def test_default_top_k(self):
        req = SearchRequest(question="Test query text")
        self.assertEqual(req.top_k, 5)

    def test_default_threshold(self):
        req = SearchRequest(question="Test query text")
        self.assertEqual(req.threshold, 0.0)

    def test_collections_optional(self):
        req = SearchRequest(question="Test query text")
        self.assertIsNone(req.collections)


# =====================================================================
# RiskResult model
# =====================================================================


class TestRiskResultModel(unittest.TestCase):
    """Tests for RiskResult response model."""

    def test_create_risk_result(self):
        result = RiskResult(
            calculator="ASCVD",
            score=7.5,
            risk_category="borderline",
            interpretation="Borderline risk",
        )
        self.assertEqual(result.calculator, "ASCVD")
        self.assertEqual(result.score, 7.5)

    def test_default_recommendations(self):
        result = RiskResult(
            calculator="test", score=0.0, risk_category="low", interpretation="test"
        )
        self.assertEqual(result.recommendations, [])

    def test_default_details(self):
        result = RiskResult(
            calculator="test", score=0.0, risk_category="low", interpretation="test"
        )
        self.assertEqual(result.details, {})


# =====================================================================
# GDMTRequest validation
# =====================================================================


class TestGDMTRequestValidation(unittest.TestCase):
    """Tests for GDMTRequest Pydantic model validation."""

    def test_valid_request(self):
        req = GDMTRequest(lvef=30.0, nyha_class="II")
        self.assertEqual(req.lvef, 30.0)

    def test_lvef_min(self):
        with self.assertRaises(ValidationError):
            GDMTRequest(lvef=3.0, nyha_class="II")

    def test_lvef_max(self):
        with self.assertRaises(ValidationError):
            GDMTRequest(lvef=85.0, nyha_class="II")

    def test_nyha_class_valid(self):
        for cls in ["I", "II", "III", "IV"]:
            req = GDMTRequest(lvef=35.0, nyha_class=cls)
            self.assertEqual(req.nyha_class, cls)

    def test_nyha_class_invalid(self):
        with self.assertRaises(ValidationError):
            GDMTRequest(lvef=35.0, nyha_class="V")


# =====================================================================
# MAGGICRequest / EuroSCORERequest
# =====================================================================


class TestMAGGICRequestValidation(unittest.TestCase):
    """Tests for MAGGICRequest Pydantic model validation."""

    def test_valid_request(self):
        req = MAGGICRequest(
            age=65, sex="male", lvef=30.0, nyha_class=2,
            systolic_bp=120.0, bmi=25.0, creatinine=1.2,
        )
        self.assertEqual(req.age, 65)

    def test_nyha_class_range(self):
        with self.assertRaises(ValidationError):
            MAGGICRequest(
                age=65, sex="male", lvef=30.0, nyha_class=5,
                systolic_bp=120.0, bmi=25.0, creatinine=1.2,
            )


class TestEuroSCORERequestValidation(unittest.TestCase):
    """Tests for EuroSCORERequest Pydantic model validation."""

    def test_valid_request(self):
        req = EuroSCORERequest(
            age=70, sex="male", creatinine_clearance=80.0,
        )
        self.assertEqual(req.age, 70)

    def test_pulmonary_hypertension_valid(self):
        for ph in ["none", "moderate", "severe"]:
            req = EuroSCORERequest(
                age=70, sex="male", creatinine_clearance=80.0,
                pulmonary_hypertension=ph,
            )
            self.assertEqual(req.pulmonary_hypertension, ph)

    def test_urgency_valid(self):
        for urg in ["elective", "urgent", "emergency", "salvage"]:
            req = EuroSCORERequest(
                age=70, sex="male", creatinine_clearance=80.0,
                urgency=urg,
            )
            self.assertEqual(req.urgency, urg)


if __name__ == "__main__":
    unittest.main()
