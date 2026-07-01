"""Tests for Biomarker Intelligence Agent disease trajectory detection engine.

Validates each disease analyzer with normal/elevated biomarkers,
genotype-stratified thresholds, and missing data handling.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.disease_trajectory import DiseaseTrajectoryAnalyzer


@pytest.fixture
def analyzer():
    """Return a fresh DiseaseTrajectoryAnalyzer instance."""
    return DiseaseTrajectoryAnalyzer()


# =====================================================================
# TYPE 2 DIABETES
# =====================================================================


class TestType2Diabetes:
    """Tests for type 2 diabetes trajectory analysis."""

    def test_normal_biomarkers_low_risk(self, analyzer):
        """Normal HbA1c and glucose should yield LOW risk."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.2, "fasting_glucose": 88},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "normal"

    def test_prediabetic_hba1c(self, analyzer):
        """HbA1c 5.7-6.4 should yield pre_diabetic stage."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.9},
            genotypes={},
        )
        assert result["risk_level"] == "HIGH"
        assert result["stage"] == "pre_diabetic"

    def test_diabetic_hba1c(self, analyzer):
        """HbA1c >= 6.5 should yield diabetic stage."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 7.2},
            genotypes={},
        )
        assert result["risk_level"] == "CRITICAL"
        assert result["stage"] == "diabetic"

    def test_tcf7l2_tt_tighter_threshold(self, analyzer):
        """TCF7L2 TT genotype should use tighter HbA1c threshold (5.5)."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.6},
            genotypes={"TCF7L2_rs7903146": "TT"},
        )
        assert result["risk_level"] == "MODERATE"
        assert result["stage"] == "early_metabolic_shift"
        assert result["genotype_adjusted_thresholds"]["hba1c"] == 5.5

    def test_tcf7l2_cc_standard_threshold(self, analyzer):
        """TCF7L2 CC genotype should use standard HbA1c threshold (6.0)."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.6},
            genotypes={"TCF7L2_rs7903146": "CC"},
        )
        assert result["risk_level"] == "LOW"
        assert result["genotype_adjusted_thresholds"]["hba1c"] == 6.0

    def test_homa_ir_calculated_from_insulin_glucose(self, analyzer):
        """HOMA-IR is calculated when fasting_insulin and fasting_glucose are provided."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"fasting_insulin": 20.0, "fasting_glucose": 110.0},
            genotypes={},
        )
        expected_homa = round((20.0 * 110.0) / 405.0, 2)
        assert result["current_markers"]["homa_ir"] == expected_homa

    def test_diabetic_fasting_glucose(self, analyzer):
        """Fasting glucose >= 126 should yield diabetic stage."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"fasting_glucose": 130},
            genotypes={},
        )
        assert result["stage"] == "diabetic"
        assert result["risk_level"] == "CRITICAL"

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={},
            genotypes={"TCF7L2_rs7903146": "TT"},
        )
        assert result["disease"] == "type2_diabetes"
        assert result["risk_level"] == "LOW"

    def test_genetic_risk_factors_detected(self, analyzer):
        """Genetic risk factors are listed when present."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.2},
            genotypes={
                "TCF7L2_rs7903146": "TT",
                "PPARG_rs1801282": "CC",
                "SLC30A8_rs13266634": "CT",
            },
        )
        genes = [f["gene"] for f in result["genetic_risk_factors"]]
        assert "TCF7L2_rs7903146" in genes


# =====================================================================
# CARDIOVASCULAR
# =====================================================================


class TestCardiovascular:
    """Tests for cardiovascular disease trajectory analysis."""

    def test_normal_lipids_low_risk(self, analyzer):
        """Normal LDL and Lp(a) should yield LOW risk."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 100, "lpa": 20, "hs_crp": 0.5},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "optimal"

    def test_elevated_lpa(self, analyzer):
        """Lp(a) > 50 should yield elevated risk."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"lpa": 75},
            genotypes={},
        )
        assert result["risk_level"] == "MODERATE"
        assert "elevated" in result["stage"]

    def test_very_high_ldl(self, analyzer):
        """LDL >= 190 should yield high risk."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 200},
            genotypes={},
        )
        assert result["risk_level"] == "HIGH"

    def test_apoe_e4e4_tighter_ldl_threshold(self, analyzer):
        """APOE E4/E4 should use LDL threshold of 100 mg/dL."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 110},
            genotypes={"APOE": "E4/E4"},
        )
        assert result["ldl_threshold"] == 100
        assert result["risk_level"] != "LOW"

    def test_apoe_e3e3_standard_ldl_threshold(self, analyzer):
        """APOE E3/E3 should use standard LDL threshold of 130 mg/dL."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 110},
            genotypes={"APOE": "E3/E3"},
        )
        assert result["ldl_threshold"] == 130
        assert result["risk_level"] == "LOW"

    def test_elevated_crp_inflammatory_risk(self, analyzer):
        """hs-CRP > 3.0 should flag high inflammatory risk."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"hs_crp": 5.0},
            genotypes={},
        )
        findings_text = " ".join(result["findings"])
        assert "inflam" in findings_text.lower()

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_cardiovascular(biomarkers={}, genotypes={})
        assert result["disease"] == "cardiovascular"

    def test_apoe_e2e2_type_iii_warning(self, analyzer):
        """APOE E2/E2 with high triglycerides flags type III dyslipidemia."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"triglycerides": 250},
            genotypes={"APOE": "E2/E2"},
        )
        findings_text = " ".join(result["findings"])
        assert "type III" in findings_text


# =====================================================================
# LIVER
# =====================================================================


class TestLiver:
    """Tests for liver disease (NAFLD/fibrosis) trajectory analysis."""

    def test_normal_liver_low_risk(self, analyzer):
        """Normal ALT and AST should yield LOW risk."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 20, "ast": 18},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "normal"

    def test_pnpla3_gg_lower_alt_threshold(self, analyzer):
        """PNPLA3 GG genotype should lower ALT upper limit to 35."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 40},
            genotypes={"PNPLA3_rs738409": "GG"},
        )
        assert result["alt_upper_limit"] == 35
        assert result["stage"] == "steatosis_risk"

    def test_pnpla3_cc_standard_alt_threshold(self, analyzer):
        """PNPLA3 CC genotype should use standard ALT upper limit of 56."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 40},
            genotypes={"PNPLA3_rs738409": "CC"},
        )
        assert result["alt_upper_limit"] == 56
        assert result["stage"] == "normal"

    def test_fib4_calculation(self, analyzer):
        """FIB-4 score is calculated when age, AST, ALT, platelets provided."""
        result = analyzer.analyze_liver(
            biomarkers={"ast": 40, "alt": 30, "platelets": 200},
            genotypes={},
            age=55.0,
        )
        assert result["fib4_score"] is not None

    def test_high_fib4_advanced_fibrosis(self, analyzer):
        """FIB-4 > 2.67 should indicate advanced fibrosis."""
        result = analyzer.analyze_liver(
            biomarkers={"ast": 80, "alt": 20, "platelets": 100},
            genotypes={},
            age=65.0,
        )
        # FIB-4 = (65 * 80) / (100 * sqrt(20)) = 5200 / 447 ~ 11.6
        assert result["stage"] == "advanced_fibrosis"
        assert result["risk_level"] == "HIGH"

    def test_elevated_ggt_flagged(self, analyzer):
        """GGT > 60 should be flagged in findings."""
        result = analyzer.analyze_liver(
            biomarkers={"ggt": 85},
            genotypes={},
        )
        findings_text = " ".join(result["findings"])
        assert "GGT" in findings_text

    def test_hsd17b13_protective(self, analyzer):
        """HSD17B13 TA variant should be noted as protective."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 25},
            genotypes={"HSD17B13_rs72613567": "T/TA"},
        )
        findings_text = " ".join(result["findings"])
        assert "protective" in findings_text.lower()

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_liver(biomarkers={}, genotypes={})
        assert result["disease"] == "liver"


# =====================================================================
# THYROID
# =====================================================================


class TestThyroid:
    """Tests for thyroid dysfunction trajectory analysis."""

    def test_normal_thyroid_euthyroid(self, analyzer):
        """Normal TSH and free T3 should yield euthyroid status."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 2.0, "free_t3": 3.0, "free_t4": 1.2},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "euthyroid"

    def test_overt_hypothyroidism(self, analyzer):
        """TSH > 10 should indicate overt hypothyroidism."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 15.0},
            genotypes={},
        )
        assert result["stage"] == "overt_dysfunction"
        assert result["risk_level"] == "HIGH"

    def test_dio2_aa_tighter_thresholds(self, analyzer):
        """DIO2 AA genotype should use tighter TSH upper limit (3.0)."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 3.2},
            genotypes={"DIO2_rs225014": "AA"},
        )
        assert result["tsh_upper_limit"] == 3.0
        assert result["stage"] == "subclinical"

    def test_dio2_gg_standard_thresholds(self, analyzer):
        """DIO2 GG genotype should use standard TSH upper limit (4.0)."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 3.2},
            genotypes={"DIO2_rs225014": "GG"},
        )
        assert result["tsh_upper_limit"] == 4.0
        assert result["stage"] == "euthyroid"

    def test_low_free_t3_with_dio2_aa(self, analyzer):
        """Low free T3 with DIO2 AA should flag impaired T4->T3 conversion."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 2.5, "free_t3": 2.2},
            genotypes={"DIO2_rs225014": "AA"},
        )
        assert result["stage"] == "subclinical"
        findings_text = " ".join(result["findings"])
        assert "DIO2" in findings_text or "T3" in findings_text

    def test_missing_free_t3_dio2_aa_warning(self, analyzer):
        """DIO2 AA without free T3 measurement should flag for monitoring."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 2.5},
            genotypes={"DIO2_rs225014": "AA"},
        )
        findings_text = " ".join(result["findings"])
        assert "T3" in findings_text

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_thyroid(biomarkers={}, genotypes={})
        assert result["disease"] == "thyroid"


# =====================================================================
# IRON
# =====================================================================


class TestIron:
    """Tests for iron metabolism disorder trajectory analysis."""

    def test_normal_iron_low_risk(self, analyzer):
        """Normal ferritin and transferrin saturation should yield LOW risk."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 100, "transferrin_saturation": 25},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "normal"

    def test_hfe_c282y_homozygous_high_risk(self, analyzer):
        """HFE C282Y homozygous (AA) with elevated ferritin should be HIGH risk."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 250},
            genotypes={"HFE_rs1800562": "AA"},
        )
        assert result["risk_level"] == "HIGH"
        assert result["stage"] == "iron_overload"

    def test_hfe_cc_standard_threshold(self, analyzer):
        """HFE CC genotype should use standard ferritin upper limit."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 300},
            genotypes={"HFE_rs1800562": "GG"},
            sex="male",
        )
        assert result["ferritin_upper_limit"] == 400
        assert result["stage"] == "normal"

    def test_sex_stratified_ferritin_thresholds(self, analyzer):
        """Female ferritin upper limit should be lower than male."""
        result_male = analyzer.analyze_iron(
            biomarkers={"ferritin": 160},
            genotypes={},
            sex="male",
        )
        result_female = analyzer.analyze_iron(
            biomarkers={"ferritin": 160},
            genotypes={},
            sex="female",
        )
        assert result_male["ferritin_upper_limit"] > result_female["ferritin_upper_limit"]

    def test_elevated_transferrin_saturation(self, analyzer):
        """Transferrin saturation > 45% should flag early accumulation."""
        result = analyzer.analyze_iron(
            biomarkers={"transferrin_saturation": 50},
            genotypes={},
        )
        assert result["stage"] == "early_accumulation"

    def test_very_high_tsat_iron_overload(self, analyzer):
        """Transferrin saturation > 60% should flag iron overload."""
        result = analyzer.analyze_iron(
            biomarkers={"transferrin_saturation": 65},
            genotypes={},
        )
        assert result["stage"] == "iron_overload"
        assert result["risk_level"] == "HIGH"

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_iron(biomarkers={}, genotypes={})
        assert result["disease"] == "iron"


# =====================================================================
# NUTRITIONAL
# =====================================================================


class TestNutritional:
    """Tests for nutritional deficiency trajectory analysis."""

    def test_optimal_nutrients_low_risk(self, analyzer):
        """Optimal nutrient levels should yield LOW risk."""
        result = analyzer.analyze_nutritional(
            biomarkers={"vitamin_d": 45, "vitamin_b12": 600, "omega3_index": 9.0},
            genotypes={},
        )
        assert result["risk_level"] == "LOW"
        assert result["stage"] == "optimal"

    def test_vitamin_d_deficiency(self, analyzer):
        """Vitamin D < 20 should flag deficiency."""
        result = analyzer.analyze_nutritional(
            biomarkers={"vitamin_d": 15},
            genotypes={},
        )
        assert "vitamin_d" in result["deficiencies"]

    def test_mthfr_tt_folate_guidance(self, analyzer):
        """MTHFR TT genotype should recommend L-methylfolate."""
        result = analyzer.analyze_nutritional(
            biomarkers={"folate": 8.0},
            genotypes={"MTHFR_rs1801133": "TT"},
        )
        findings_text = " ".join(result["findings"])
        assert "MTHFR" in findings_text
        recs_text = " ".join(result["recommendations"])
        assert "methylfolate" in recs_text.lower()

    def test_fads1_cc_omega3_target(self, analyzer):
        """FADS1 CC genotype should set omega-3 target to 5.0%."""
        result = analyzer.analyze_nutritional(
            biomarkers={"omega3_index": 4.5},
            genotypes={"FADS1_rs174546": "CC"},
        )
        assert result["omega3_target"] == 5.0
        assert "omega3" in result["deficiencies"]

    def test_multiple_deficiencies_high_risk(self, analyzer):
        """3+ deficiencies should yield HIGH risk."""
        result = analyzer.analyze_nutritional(
            biomarkers={
                "vitamin_d": 12,
                "vitamin_b12": 150,
                "omega3_index": 3.0,
                "magnesium": 1.5,
            },
            genotypes={},
        )
        assert result["risk_level"] == "HIGH"
        assert result["stage"] == "deficient"

    def test_missing_biomarkers_no_crash(self, analyzer):
        """Empty biomarkers should not crash."""
        result = analyzer.analyze_nutritional(biomarkers={}, genotypes={})
        assert result["disease"] == "nutritional"


# =====================================================================
# ANALYZE ALL
# =====================================================================


class TestAnalyzeAll:
    """Tests for the combined analyze_all() method."""

    def test_returns_nine_results(self, analyzer):
        """analyze_all() returns exactly 9 disease trajectory results."""
        results = analyzer.analyze_all(
            biomarkers={"hba1c": 5.4, "ldl_c": 120},
            genotypes={},
        )
        assert len(results) == 9

    def test_sorted_by_risk_level(self, analyzer):
        """Results should be sorted by risk severity (CRITICAL first)."""
        results = analyzer.analyze_all(
            biomarkers={"hba1c": 7.0, "ldl_c": 120, "alt": 20},
            genotypes={},
        )
        risk_order = {"CRITICAL": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3}
        risk_values = [risk_order.get(r["risk_level"], 4) for r in results]
        assert risk_values == sorted(risk_values)

    def test_all_diseases_represented(self, analyzer):
        """All 9 disease categories should be in results."""
        results = analyzer.analyze_all(biomarkers={}, genotypes={})
        diseases = {r["disease"] for r in results}
        expected = {"type2_diabetes", "cardiovascular", "liver", "thyroid", "iron", "nutritional", "kidney", "bone_health", "cognitive"}
        assert diseases == expected
