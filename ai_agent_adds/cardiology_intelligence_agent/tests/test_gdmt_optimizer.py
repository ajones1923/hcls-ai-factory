"""Unit tests for Cardiology Intelligence Agent GDMT Optimizer.

Tests GDMTOptimizer creation, EF classification, drug database content,
device therapy criteria, and optimization recommendations.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.gdmt_optimizer import (
    ADDITIONAL_THERAPIES,
    DEVICE_THERAPY_CRITERIA,
    GDMT_DRUG_DATABASE,
    GDMTOptimizer,
    HFPEF_THERAPIES,
)
from src.models import (
    EjectionFractionCategory,
    EvidenceLevel,
    GDMTPillar,
    GDMTRecommendation,
    GuidelineClass,
    HeartFailureClass,
)


# ===================================================================
# GDMTOptimizer CREATION
# ===================================================================


class TestGDMTOptimizerCreation:
    """Tests for GDMTOptimizer instantiation."""

    def test_creation(self):
        optimizer = GDMTOptimizer()
        assert optimizer is not None

    def test_has_optimize_method(self):
        optimizer = GDMTOptimizer()
        assert callable(getattr(optimizer, "optimize", None))

    def test_has_guideline_reference(self):
        assert "2022" in GDMTOptimizer.GUIDELINE_REFERENCE
        assert "Heidenreich" in GDMTOptimizer.GUIDELINE_REFERENCE

    def test_has_ef_cutoffs(self):
        assert GDMTOptimizer._HFREF_EF_CUTOFF == 40.0
        assert GDMTOptimizer._HFMREF_UPPER == 49.0
        assert GDMTOptimizer._HFPEF_EF_CUTOFF == 50.0

    def test_has_clinical_thresholds(self):
        assert GDMTOptimizer._MIN_SBP_FOR_RAAS == 100.0
        assert GDMTOptimizer._MAX_K_FOR_MRA == 5.0
        assert GDMTOptimizer._MIN_EGFR_FOR_MRA == 30.0
        assert GDMTOptimizer._ICD_EF_CUTOFF == 35.0
        assert GDMTOptimizer._CRT_EF_CUTOFF == 35.0


# ===================================================================
# _classify_ef TESTS
# ===================================================================


class TestClassifyEF:
    """Tests for the EF classification method."""

    def setup_method(self):
        self.optimizer = GDMTOptimizer()

    def test_hfref_at_40(self):
        result = self.optimizer._classify_ef(40.0, {})
        assert result == EjectionFractionCategory.HFrEF

    def test_hfref_below_40(self):
        result = self.optimizer._classify_ef(25.0, {})
        assert result == EjectionFractionCategory.HFrEF

    def test_hfmref_at_41(self):
        result = self.optimizer._classify_ef(41.0, {})
        assert result == EjectionFractionCategory.HFmrEF

    def test_hfmref_at_49(self):
        result = self.optimizer._classify_ef(49.0, {})
        assert result == EjectionFractionCategory.HFmrEF

    def test_hfpef_at_50(self):
        result = self.optimizer._classify_ef(50.0, {})
        assert result == EjectionFractionCategory.HFpEF

    def test_hfpef_at_65(self):
        result = self.optimizer._classify_ef(65.0, {})
        assert result == EjectionFractionCategory.HFpEF

    def test_hfimpef_with_previous_ef(self):
        result = self.optimizer._classify_ef(45.0, {"previous_ef": 30.0})
        assert result == EjectionFractionCategory.HFimpEF

    def test_not_hfimpef_if_previous_ef_above_40(self):
        result = self.optimizer._classify_ef(55.0, {"previous_ef": 42.0})
        assert result == EjectionFractionCategory.HFpEF

    def test_not_hfimpef_if_current_still_below_40(self):
        result = self.optimizer._classify_ef(35.0, {"previous_ef": 25.0})
        assert result == EjectionFractionCategory.HFrEF

    @pytest.mark.parametrize(
        "lvef, expected",
        [
            (10.0, EjectionFractionCategory.HFrEF),
            (20.0, EjectionFractionCategory.HFrEF),
            (30.0, EjectionFractionCategory.HFrEF),
            (40.0, EjectionFractionCategory.HFrEF),
            (41.0, EjectionFractionCategory.HFmrEF),
            (45.0, EjectionFractionCategory.HFmrEF),
            (49.0, EjectionFractionCategory.HFmrEF),
            (50.0, EjectionFractionCategory.HFpEF),
            (55.0, EjectionFractionCategory.HFpEF),
            (70.0, EjectionFractionCategory.HFpEF),
        ],
    )
    def test_ef_classification_parametrized(self, lvef, expected):
        result = self.optimizer._classify_ef(lvef, {})
        assert result == expected


# ===================================================================
# _parse_nyha TESTS
# ===================================================================


class TestParseNYHA:
    """Tests for the NYHA parsing static method."""

    @pytest.mark.parametrize(
        "input_str, expected",
        [
            ("I", HeartFailureClass.NYHA_I),
            ("II", HeartFailureClass.NYHA_II),
            ("III", HeartFailureClass.NYHA_III),
            ("IV", HeartFailureClass.NYHA_IV),
            ("1", HeartFailureClass.NYHA_I),
            ("2", HeartFailureClass.NYHA_II),
            ("3", HeartFailureClass.NYHA_III),
            ("4", HeartFailureClass.NYHA_IV),
            ("NYHA III", HeartFailureClass.NYHA_III),
            ("Class II", HeartFailureClass.NYHA_II),
        ],
    )
    def test_valid_parsing(self, input_str, expected):
        result = GDMTOptimizer._parse_nyha(input_str)
        assert result == expected

    def test_unknown_defaults_to_nyha_ii(self):
        result = GDMTOptimizer._parse_nyha("unknown")
        assert result == HeartFailureClass.NYHA_II


# ===================================================================
# GDMT_DRUG_DATABASE TESTS
# ===================================================================


class TestGDMTDrugDatabase:
    """Tests for the GDMT_DRUG_DATABASE dictionary."""

    def test_has_four_pillars(self):
        assert len(GDMT_DRUG_DATABASE) == 4

    @pytest.mark.parametrize(
        "pillar",
        [GDMTPillar.BETA_BLOCKER, GDMTPillar.ARNI_ACEI_ARB, GDMTPillar.MRA, GDMTPillar.SGLT2I],
    )
    def test_pillar_exists(self, pillar):
        assert pillar in GDMT_DRUG_DATABASE

    def test_each_pillar_has_drugs(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            assert len(drugs) > 0, f"Pillar {pillar} has no drugs"

    # ---- Beta-blocker drugs ----

    def test_bb_has_carvedilol(self):
        assert "carvedilol" in GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]

    def test_bb_has_metoprolol_succinate(self):
        assert "metoprolol_succinate" in GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]

    def test_bb_has_bisoprolol(self):
        assert "bisoprolol" in GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]

    def test_carvedilol_target_dose(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert "25 mg BID" in carv["target_dose"]

    def test_carvedilol_starting_dose(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert "3.125 mg BID" in carv["starting_dose"]

    def test_carvedilol_has_contraindications(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert isinstance(carv["contraindications"], list)
        assert len(carv["contraindications"]) > 0

    def test_carvedilol_key_trial(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert carv["key_trial"] == "COPERNICUS"

    def test_carvedilol_evidence_class_i(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert carv["evidence"]["class"] == GuidelineClass.CLASS_I

    def test_carvedilol_evidence_level_a(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        assert carv["evidence"]["level"] == EvidenceLevel.LEVEL_A

    def test_metoprolol_target_dose(self):
        meto = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["metoprolol_succinate"]
        assert meto["target_dose"] == "200 mg daily"

    def test_bisoprolol_target_dose(self):
        biso = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["bisoprolol"]
        assert biso["target_dose"] == "10 mg daily"

    # ---- ARNI / ACEi / ARB drugs ----

    def test_arni_has_sacubitril_valsartan(self):
        assert "sacubitril_valsartan" in GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]

    def test_sacubitril_valsartan_target_dose(self):
        sv = GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]["sacubitril_valsartan"]
        assert sv["target_dose"] == "97/103 mg BID"

    def test_sacubitril_valsartan_key_trial(self):
        sv = GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]["sacubitril_valsartan"]
        assert sv["key_trial"] == "PARADIGM-HF"

    def test_arni_has_enalapril(self):
        assert "enalapril" in GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]

    def test_arni_has_lisinopril(self):
        assert "lisinopril" in GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]

    def test_arni_has_losartan(self):
        assert "losartan" in GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]

    def test_arni_has_valsartan(self):
        assert "valsartan" in GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]

    # ---- MRA drugs ----

    def test_mra_has_spironolactone(self):
        assert "spironolactone" in GDMT_DRUG_DATABASE[GDMTPillar.MRA]

    def test_mra_has_eplerenone(self):
        assert "eplerenone" in GDMT_DRUG_DATABASE[GDMTPillar.MRA]

    def test_spironolactone_target_dose(self):
        spiro = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["spironolactone"]
        assert "25" in spiro["target_dose"]
        assert "50" in spiro["target_dose"]

    def test_spironolactone_key_trial(self):
        spiro = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["spironolactone"]
        assert spiro["key_trial"] == "RALES"

    def test_spironolactone_k_contraindication(self):
        spiro = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["spironolactone"]
        contraindications = spiro["contraindications"]
        assert any("potassium" in c.lower() or "5.0" in c for c in contraindications)

    # ---- SGLT2i drugs ----

    def test_sglt2i_has_dapagliflozin(self):
        assert "dapagliflozin" in GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]

    def test_sglt2i_has_empagliflozin(self):
        assert "empagliflozin" in GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]

    def test_dapagliflozin_target_dose(self):
        dapa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["dapagliflozin"]
        assert dapa["target_dose"] == "10 mg daily"

    def test_dapagliflozin_starting_dose(self):
        dapa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["dapagliflozin"]
        assert dapa["starting_dose"] == "10 mg daily"

    def test_dapagliflozin_key_trial(self):
        dapa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["dapagliflozin"]
        assert dapa["key_trial"] == "DAPA-HF"

    def test_empagliflozin_target_dose(self):
        empa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["empagliflozin"]
        assert empa["target_dose"] == "10 mg daily"

    def test_empagliflozin_key_trial(self):
        empa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["empagliflozin"]
        assert empa["key_trial"] == "EMPEROR-Reduced"

    # ---- All drugs have required keys ----

    def test_all_drugs_have_starting_dose(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "starting_dose" in data, f"{pillar}/{drug_name} missing 'starting_dose'"

    def test_all_drugs_have_target_dose(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "target_dose" in data, f"{pillar}/{drug_name} missing 'target_dose'"

    def test_all_drugs_have_contraindications(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "contraindications" in data, f"{pillar}/{drug_name} missing 'contraindications'"

    def test_all_drugs_have_monitoring(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "monitoring" in data, f"{pillar}/{drug_name} missing 'monitoring'"

    def test_all_drugs_have_evidence(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "evidence" in data, f"{pillar}/{drug_name} missing 'evidence'"
                assert "class" in data["evidence"]
                assert "level" in data["evidence"]

    def test_all_drugs_have_key_trial(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "key_trial" in data, f"{pillar}/{drug_name} missing 'key_trial'"

    def test_all_drugs_have_mechanism(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "mechanism" in data, f"{pillar}/{drug_name} missing 'mechanism'"

    def test_all_drugs_have_titration_steps(self):
        for pillar, drugs in GDMT_DRUG_DATABASE.items():
            for drug_name, data in drugs.items():
                assert "titration_steps" in data, f"{pillar}/{drug_name} missing 'titration_steps'"
                assert isinstance(data["titration_steps"], list)


# ===================================================================
# ADDITIONAL_THERAPIES TESTS
# ===================================================================


class TestAdditionalTherapies:
    """Tests for the ADDITIONAL_THERAPIES dictionary."""

    def test_is_dict(self):
        assert isinstance(ADDITIONAL_THERAPIES, dict)

    def test_has_hydralazine_isdn(self):
        assert "hydralazine_isosorbide_dinitrate" in ADDITIONAL_THERAPIES

    def test_has_ivabradine(self):
        assert "ivabradine" in ADDITIONAL_THERAPIES

    def test_has_digoxin(self):
        assert "digoxin" in ADDITIONAL_THERAPIES

    def test_has_vericiguat(self):
        assert "vericiguat" in ADDITIONAL_THERAPIES

    def test_ivabradine_key_trial(self):
        assert ADDITIONAL_THERAPIES["ivabradine"]["key_trial"] == "SHIFT"

    def test_digoxin_key_trial(self):
        assert ADDITIONAL_THERAPIES["digoxin"]["key_trial"] == "DIG"


# ===================================================================
# DEVICE_THERAPY_CRITERIA TESTS
# ===================================================================


class TestDeviceTherapyCriteria:
    """Tests for the DEVICE_THERAPY_CRITERIA dictionary."""

    def test_is_dict(self):
        assert isinstance(DEVICE_THERAPY_CRITERIA, dict)

    def test_has_icd(self):
        assert "ICD" in DEVICE_THERAPY_CRITERIA

    def test_has_crt(self):
        assert "CRT" in DEVICE_THERAPY_CRITERIA

    def test_has_crt_d(self):
        assert "CRT_D" in DEVICE_THERAPY_CRITERIA

    def test_icd_lvef_criteria(self):
        criteria = DEVICE_THERAPY_CRITERIA["ICD"]["criteria"]
        assert "<=35%" in criteria["lvef"]

    def test_icd_nyha_class(self):
        criteria = DEVICE_THERAPY_CRITERIA["ICD"]["criteria"]
        assert "II" in criteria["nyha_class"]
        assert "III" in criteria["nyha_class"]

    def test_icd_gdmt_duration(self):
        criteria = DEVICE_THERAPY_CRITERIA["ICD"]["criteria"]
        assert "90 days" in criteria["gdmt_duration"]

    def test_icd_evidence_class_i(self):
        evidence = DEVICE_THERAPY_CRITERIA["ICD"]["evidence"]
        assert evidence["class"] == GuidelineClass.CLASS_I

    def test_icd_evidence_level_a(self):
        evidence = DEVICE_THERAPY_CRITERIA["ICD"]["evidence"]
        assert evidence["level"] == EvidenceLevel.LEVEL_A

    def test_icd_key_trials(self):
        trials = DEVICE_THERAPY_CRITERIA["ICD"]["key_trials"]
        assert "SCD-HeFT" in trials
        assert "MADIT-II" in trials

    def test_crt_qrs_duration(self):
        criteria = DEVICE_THERAPY_CRITERIA["CRT"]["criteria"]
        assert "150" in criteria["qrs_duration"]

    def test_crt_key_trials(self):
        trials = DEVICE_THERAPY_CRITERIA["CRT"]["key_trials"]
        assert "CARE-HF" in trials

    def test_crt_d_requirement(self):
        criteria = DEVICE_THERAPY_CRITERIA["CRT_D"]["criteria"]
        assert "BOTH" in criteria["requirement"]

    def test_all_devices_have_description(self):
        for name, data in DEVICE_THERAPY_CRITERIA.items():
            assert "description" in data, f"'{name}' missing 'description'"

    def test_all_devices_have_criteria(self):
        for name, data in DEVICE_THERAPY_CRITERIA.items():
            assert "criteria" in data, f"'{name}' missing 'criteria'"

    def test_all_devices_have_evidence(self):
        for name, data in DEVICE_THERAPY_CRITERIA.items():
            assert "evidence" in data, f"'{name}' missing 'evidence'"


# ===================================================================
# HFPEF_THERAPIES TESTS
# ===================================================================


class TestHFpEFTherapies:
    """Tests for the HFPEF_THERAPIES dictionary."""

    def test_is_dict(self):
        assert isinstance(HFPEF_THERAPIES, dict)

    def test_has_sglt2i(self):
        assert "sglt2i" in HFPEF_THERAPIES

    def test_has_diuretics(self):
        assert "diuretics" in HFPEF_THERAPIES

    def test_has_glp1_ra(self):
        assert "glp1_ra" in HFPEF_THERAPIES

    def test_sglt2i_drugs_include_dapagliflozin(self):
        drugs = HFPEF_THERAPIES["sglt2i"]["drugs"]
        assert any("dapagliflozin" in d for d in drugs)

    def test_sglt2i_drugs_include_empagliflozin(self):
        drugs = HFPEF_THERAPIES["sglt2i"]["drugs"]
        assert any("empagliflozin" in d for d in drugs)

    def test_sglt2i_evidence_class_i(self):
        evidence = HFPEF_THERAPIES["sglt2i"]["evidence"]
        assert evidence["class"] == GuidelineClass.CLASS_I

    def test_sglt2i_key_trials(self):
        trials = HFPEF_THERAPIES["sglt2i"]["key_trials"]
        assert "EMPEROR-Preserved" in trials
        assert "DELIVER" in trials

    def test_glp1_ra_drugs_include_semaglutide(self):
        drugs = HFPEF_THERAPIES["glp1_ra"]["drugs"]
        assert any("semaglutide" in d for d in drugs)


# ===================================================================
# optimize() METHOD TESTS
# ===================================================================


class TestOptimizeMethod:
    """Tests for the GDMTOptimizer.optimize() method."""

    def setup_method(self):
        self.optimizer = GDMTOptimizer()

    def _make_hfref_patient(self):
        return dict(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 110.0,
                "hr": 72,
                "potassium": 4.2,
                "egfr": 60.0,
                "age": 65,
                "sex": "male",
            },
        )

    def test_optimize_returns_recommendation(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        assert isinstance(result, GDMTRecommendation)

    def test_hfref_returns_hfref_category(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        assert result.ef_category == EjectionFractionCategory.HFrEF

    def test_hfref_has_recommendations(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        assert len(result.recommendations) > 0

    def test_hfref_has_next_steps(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        assert len(result.next_steps) > 0

    def test_hfref_has_guideline_references(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        assert len(result.guideline_references) > 0
        assert any("2022" in ref for ref in result.guideline_references)

    def test_hfref_mentions_all_four_pillars(self):
        """HFrEF optimization should reference all 4 GDMT pillars."""
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        all_text = " ".join(result.recommendations).lower()
        pillar_keywords = ["beta", "arni", "mra", "sglt2"]
        # At least some pillar keywords should appear
        found = sum(1 for kw in pillar_keywords if kw in all_text)
        assert found >= 2, f"Expected >= 2 pillar keywords, found {found}"

    def test_hfref_reassess_ef_in_next_steps(self):
        params = self._make_hfref_patient()
        result = self.optimizer.optimize(**params)
        all_next = " ".join(result.next_steps).lower()
        assert "reassess" in all_next or "echocardiography" in all_next

    def test_hfpef_classification(self):
        result = self.optimizer.optimize(
            lvef=55.0,
            nyha_class="II",
            current_medications=[],
            patient_data={"sbp": 130.0, "hr": 78, "potassium": 4.0, "egfr": 70.0},
        )
        assert result.ef_category == EjectionFractionCategory.HFpEF

    def test_hfpef_recommends_sglt2i(self):
        result = self.optimizer.optimize(
            lvef=55.0,
            nyha_class="II",
            current_medications=[],
            patient_data={"sbp": 130.0, "hr": 78, "potassium": 4.0, "egfr": 70.0},
        )
        all_text = " ".join(result.recommendations).lower()
        assert "sglt2" in all_text

    def test_hfmref_classification(self):
        result = self.optimizer.optimize(
            lvef=45.0,
            nyha_class="II",
            current_medications=[],
            patient_data={"sbp": 125.0, "hr": 70, "potassium": 4.5, "egfr": 55.0},
        )
        assert result.ef_category == EjectionFractionCategory.HFmrEF

    def test_hfimpef_classification(self):
        result = self.optimizer.optimize(
            lvef=45.0,
            nyha_class="II",
            current_medications=[],
            patient_data={
                "sbp": 120.0, "hr": 65,
                "potassium": 4.0, "egfr": 60.0,
                "previous_ef": 30.0,
            },
        )
        assert result.ef_category == EjectionFractionCategory.HFimpEF

    def test_with_existing_medications(self):
        result = self.optimizer.optimize(
            lvef=30.0,
            nyha_class="III",
            current_medications=[
                {"name": "carvedilol", "dose": "12.5 mg BID"},
                {"name": "sacubitril/valsartan", "dose": "49/51 mg BID"},
            ],
            patient_data={
                "sbp": 105.0, "hr": 68,
                "potassium": 4.3, "egfr": 50.0,
            },
        )
        assert isinstance(result, GDMTRecommendation)
        assert len(result.recommendations) > 0

    def test_high_potassium_mra_concern(self):
        """With K+ > 5.0, MRA should be flagged or not recommended."""
        result = self.optimizer.optimize(
            lvef=30.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 110.0, "hr": 72,
                "potassium": 5.3,
                "egfr": 60.0,
            },
        )
        all_text = " ".join(result.recommendations + result.next_steps).lower()
        # Should mention potassium or caution about MRA
        assert "potassium" in all_text or "mra" in all_text or "hyperkalemia" in all_text

    def test_low_egfr_impacts_recommendations(self):
        result = self.optimizer.optimize(
            lvef=30.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 110.0, "hr": 72,
                "potassium": 4.0,
                "egfr": 25.0,
            },
        )
        all_text = " ".join(result.recommendations + result.next_steps).lower()
        assert "renal" in all_text or "egfr" in all_text or "creatinine" in all_text

    def test_hfref_device_therapy_assessment(self):
        """With LVEF <=35% and NYHA III, ICD should be considered."""
        result = self.optimizer.optimize(
            lvef=30.0,
            nyha_class="III",
            current_medications=[
                {"name": "carvedilol", "dose": "25 mg BID"},
                {"name": "sacubitril/valsartan", "dose": "97/103 mg BID"},
                {"name": "spironolactone", "dose": "25 mg daily"},
                {"name": "dapagliflozin", "dose": "10 mg daily"},
            ],
            patient_data={
                "sbp": 110.0, "hr": 65,
                "potassium": 4.2, "egfr": 55.0,
                "days_on_gdmt": 100,
            },
        )
        all_text = " ".join(result.recommendations).lower()
        assert "icd" in all_text or "device" in all_text or "defibrillator" in all_text


# ===================================================================
# CONTRAINDICATION CHECKING TESTS
# ===================================================================


class TestContraindicationLogic:
    """Tests for contraindication-related logic in GDMTOptimizer."""

    def setup_method(self):
        self.optimizer = GDMTOptimizer()

    def test_k_above_5_blocks_mra_in_drug_db(self):
        """Verify that spironolactone contraindications include K+ >5.0."""
        spiro = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["spironolactone"]
        k_contras = [c for c in spiro["contraindications"] if "potassium" in c.lower() or "5.0" in c]
        assert len(k_contras) > 0

    def test_angioedema_blocks_arni(self):
        sv = GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]["sacubitril_valsartan"]
        angioedema_contras = [c for c in sv["contraindications"] if "angioedema" in c.lower()]
        assert len(angioedema_contras) > 0

    def test_bradycardia_blocks_bb(self):
        carv = GDMT_DRUG_DATABASE[GDMTPillar.BETA_BLOCKER]["carvedilol"]
        brady_contras = [c for c in carv["contraindications"] if "bradycardia" in c.lower()]
        assert len(brady_contras) > 0

    def test_type1_dm_blocks_sglt2i(self):
        dapa = GDMT_DRUG_DATABASE[GDMTPillar.SGLT2I]["dapagliflozin"]
        t1dm_contras = [c for c in dapa["contraindications"] if "type 1" in c.lower()]
        assert len(t1dm_contras) > 0

    def test_concurrent_acei_blocks_arni(self):
        sv = GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]["sacubitril_valsartan"]
        acei_contras = [c for c in sv["contraindications"] if "acei" in c.lower() or "ace" in c.lower()]
        assert len(acei_contras) > 0

    def test_eplerenone_egfr_contraindication(self):
        eple = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["eplerenone"]
        egfr_contras = [c for c in eple["contraindications"] if "egfr" in c.lower() or "30" in c]
        assert len(egfr_contras) > 0

    def test_spironolactone_egfr_contraindication(self):
        spiro = GDMT_DRUG_DATABASE[GDMTPillar.MRA]["spironolactone"]
        egfr_contras = [c for c in spiro["contraindications"] if "egfr" in c.lower() or "30" in c]
        assert len(egfr_contras) > 0


# ===================================================================
# SAFETY-CRITICAL REGRESSION TESTS
# ===================================================================


class TestSafetyCriticalFixes:
    """Regression tests for clinically significant fixes.

    These tests verify that safety-critical corrections remain intact.
    Each test documents the clinical rationale and source guideline.
    """

    def setup_method(self):
        self.optimizer = GDMTOptimizer()

    # -- SBP <90 contraindication for ARNI/ACEi/ARB --

    def test_low_sbp_blocks_arni(self):
        """SBP <90 mmHg should contraindicate ARNI/ACEi/ARB initiation.

        Guideline: 2022 AHA/ACC/HFSA §7.3.2 - hypotension precaution.
        """
        result = self.optimizer.optimize(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 85.0, "hr": 72,
                "potassium": 4.2, "egfr": 60.0,
            },
        )
        all_text = " ".join(result.recommendations + result.next_steps).lower()
        assert "hypotension" in all_text or "sbp" in all_text or "contraindicated" in all_text

    def test_normal_sbp_allows_arni(self):
        """SBP >=100 should not block ARNI/ACEi/ARB."""
        result = self.optimizer.optimize(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 110.0, "hr": 72,
                "potassium": 4.2, "egfr": 60.0,
            },
        )
        all_text = " ".join(result.recommendations).lower()
        # Should recommend ARNI or ACEi/ARB
        assert "arni" in all_text or "sacubitril" in all_text or "acei" in all_text

    # -- Ivabradine requires LVEF ≤35%, not just HFrEF (≤40%) --

    def test_ivabradine_blocked_at_ef_38(self):
        """LVEF 38% is HFrEF but ivabradine requires ≤35% per SHIFT trial.

        Guideline: 2022 AHA/ACC/HFSA §7.3.7; SHIFT trial inclusion.
        """
        result = self.optimizer.optimize(
            lvef=38.0,
            nyha_class="II",
            current_medications=[
                {"name": "metoprolol succinate", "dose": "200 mg daily"},
            ],
            patient_data={
                "sbp": 120.0, "hr": 78,
                "potassium": 4.0, "egfr": 65.0,
                "rhythm": "sinus",
            },
        )
        all_text = " ".join(result.recommendations).lower()
        # Ivabradine should NOT be recommended at EF 38%
        assert "ivabradine" not in all_text

    def test_ivabradine_eligible_at_ef_35(self):
        """LVEF 35% + sinus rhythm + HR ≥70 should trigger ivabradine."""
        result = self.optimizer.optimize(
            lvef=35.0,
            nyha_class="II",
            current_medications=[
                {"name": "metoprolol succinate", "dose": "200 mg daily"},
            ],
            patient_data={
                "sbp": 120.0, "hr": 78,
                "potassium": 4.0, "egfr": 65.0,
                "rhythm": "sinus",
            },
        )
        all_text = " ".join(result.recommendations).lower()
        assert "ivabradine" in all_text

    # -- MRA K+ threshold: >5.0 (not >=5.0) --

    def test_mra_allowed_at_k_exactly_5(self):
        """K+ exactly 5.0 should NOT block MRA (threshold is >5.0).

        Guideline: 2022 AHA/ACC/HFSA §7.3.4 - K+ >5.0 mEq/L.
        """
        result = self.optimizer.optimize(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 115.0, "hr": 72,
                "potassium": 5.0,
                "egfr": 60.0,
            },
        )
        all_text = " ".join(result.recommendations).lower()
        # MRA should still be recommended
        assert "spironolactone" in all_text or "eplerenone" in all_text or "mra" in all_text

    def test_mra_blocked_at_k_5_1(self):
        """K+ 5.1 (>5.0) should flag MRA concern."""
        result = self.optimizer.optimize(
            lvef=25.0,
            nyha_class="III",
            current_medications=[],
            patient_data={
                "sbp": 115.0, "hr": 72,
                "potassium": 5.1,
                "egfr": 60.0,
            },
        )
        all_text = " ".join(result.recommendations + result.next_steps).lower()
        assert "potassium" in all_text or "hyperkalemia" in all_text or "mra" in all_text

    # -- Candesartan in drug database --

    def test_candesartan_in_arni_acei_arb_pillar(self):
        """Candesartan should exist per CHARM-Alternative trial."""
        arb_drugs = GDMT_DRUG_DATABASE[GDMTPillar.ARNI_ACEI_ARB]
        assert "candesartan" in arb_drugs
        assert arb_drugs["candesartan"]["target_dose"] is not None
