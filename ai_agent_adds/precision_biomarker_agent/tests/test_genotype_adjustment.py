"""Tests for Biomarker Intelligence Agent genotype-based reference range adjustments.

Validates that genotype-adjusted thresholds are applied correctly for
PNPLA3, TCF7L2, APOE, DIO2, HFE, and FADS1, and that missing genotypes
produce no adjustment.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.disease_trajectory import DiseaseTrajectoryAnalyzer
from src.models import GenotypeAdjustment, GenotypeAdjustmentResult


@pytest.fixture
def analyzer():
    """Return a fresh DiseaseTrajectoryAnalyzer instance."""
    return DiseaseTrajectoryAnalyzer()


# =====================================================================
# PNPLA3 ALT ADJUSTMENT
# =====================================================================


class TestPNPLA3Adjustment:
    """Tests for PNPLA3-adjusted ALT reference ranges."""

    def test_pnpla3_cc_standard_alt_range(self, analyzer):
        """PNPLA3 CC should use standard ALT upper limit of 56 U/L."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 50},
            genotypes={"PNPLA3_rs738409": "CC"},
        )
        assert result["alt_upper_limit"] == 56
        assert result["stage"] == "normal"

    def test_pnpla3_cg_intermediate_alt_range(self, analyzer):
        """PNPLA3 CG should use intermediate ALT upper limit of 45 U/L."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 48},
            genotypes={"PNPLA3_rs738409": "CG"},
        )
        assert result["alt_upper_limit"] == 45
        assert result["stage"] == "steatosis_risk"

    def test_pnpla3_gg_adjusted_alt_upper_35(self, analyzer):
        """PNPLA3 GG should adjust ALT upper limit to 35 U/L."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 38},
            genotypes={"PNPLA3_rs738409": "GG"},
        )
        assert result["alt_upper_limit"] == 35
        assert result["stage"] == "steatosis_risk"

    def test_pnpla3_gg_below_threshold_normal(self, analyzer):
        """PNPLA3 GG with ALT < 35 should still be normal."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 30},
            genotypes={"PNPLA3_rs738409": "GG"},
        )
        assert result["alt_upper_limit"] == 35
        assert result["stage"] == "normal"


# =====================================================================
# TCF7L2 GLUCOSE ADJUSTMENT
# =====================================================================


class TestTCF7L2Adjustment:
    """Tests for TCF7L2-adjusted glucose/HbA1c thresholds."""

    def test_tcf7l2_cc_standard_glucose_threshold(self, analyzer):
        """TCF7L2 CC should use standard fasting glucose threshold of 100."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"fasting_glucose": 98},
            genotypes={"TCF7L2_rs7903146": "CC"},
        )
        assert result["genotype_adjusted_thresholds"]["fasting_glucose"] == 100
        assert result["risk_level"] == "LOW"

    def test_tcf7l2_tt_tighter_glucose_threshold(self, analyzer):
        """TCF7L2 TT should use tighter fasting glucose threshold of 90."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"fasting_glucose": 92},
            genotypes={"TCF7L2_rs7903146": "TT"},
        )
        assert result["genotype_adjusted_thresholds"]["fasting_glucose"] == 90
        # 92 exceeds 90, so findings should note it
        assert len(result["findings"]) > 0

    def test_tcf7l2_ct_intermediate_threshold(self, analyzer):
        """TCF7L2 CT should use intermediate fasting glucose threshold of 95."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"fasting_glucose": 96},
            genotypes={"TCF7L2_rs7903146": "CT"},
        )
        assert result["genotype_adjusted_thresholds"]["fasting_glucose"] == 95

    def test_tcf7l2_tt_hba1c_threshold(self, analyzer):
        """TCF7L2 TT should use tighter HbA1c threshold of 5.5%."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.6},
            genotypes={"TCF7L2_rs7903146": "TT"},
        )
        assert result["genotype_adjusted_thresholds"]["hba1c"] == 5.5
        assert result["risk_level"] == "MODERATE"

    def test_tcf7l2_cc_hba1c_threshold(self, analyzer):
        """TCF7L2 CC should use standard HbA1c threshold of 6.0%."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.6},
            genotypes={"TCF7L2_rs7903146": "CC"},
        )
        assert result["genotype_adjusted_thresholds"]["hba1c"] == 6.0
        assert result["risk_level"] == "LOW"


# =====================================================================
# APOE LDL ADJUSTMENT
# =====================================================================


class TestAPOEAdjustment:
    """Tests for APOE-adjusted LDL targets."""

    def test_apoe_e4e4_aggressive_ldl_target(self, analyzer):
        """APOE E4/E4 should use aggressive LDL threshold of 100 mg/dL."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 105},
            genotypes={"APOE": "E4/E4"},
        )
        assert result["ldl_threshold"] == 100
        assert result["risk_level"] != "LOW"

    def test_apoe_e3e4_intermediate_ldl_target(self, analyzer):
        """APOE E3/E4 should use LDL threshold of 115 mg/dL."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 120},
            genotypes={"APOE": "E3/E4"},
        )
        assert result["ldl_threshold"] == 115

    def test_apoe_e3e3_standard_ldl_target(self, analyzer):
        """APOE E3/E3 should use standard LDL threshold of 130 mg/dL."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 120},
            genotypes={"APOE": "E3/E3"},
        )
        assert result["ldl_threshold"] == 130
        assert result["risk_level"] == "LOW"

    def test_apoe_e2e2_special_case(self, analyzer):
        """APOE E2/E2 should use standard threshold but flag type III risk."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 100, "triglycerides": 250},
            genotypes={"APOE": "E2/E2"},
        )
        assert result["ldl_threshold"] == 130
        findings_text = " ".join(result["findings"])
        assert "type III" in findings_text


# =====================================================================
# DIO2 THYROID ADJUSTMENT
# =====================================================================


class TestDIO2Adjustment:
    """Tests for DIO2-adjusted thyroid reference ranges."""

    def test_dio2_gg_standard_tsh(self, analyzer):
        """DIO2 GG should use standard TSH upper limit of 4.0 mIU/L."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 3.5},
            genotypes={"DIO2_rs225014": "GG"},
        )
        assert result["tsh_upper_limit"] == 4.0
        assert result["stage"] == "euthyroid"

    def test_dio2_ga_intermediate_tsh(self, analyzer):
        """DIO2 GA should use TSH upper limit of 3.5 mIU/L."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 3.6},
            genotypes={"DIO2_rs225014": "GA"},
        )
        assert result["tsh_upper_limit"] == 3.5
        assert result["stage"] == "subclinical"

    def test_dio2_aa_tight_tsh_and_ft3(self, analyzer):
        """DIO2 AA should use TSH upper limit of 3.0 and free T3 lower of 2.8."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 2.5, "free_t3": 2.5},
            genotypes={"DIO2_rs225014": "AA"},
        )
        assert result["tsh_upper_limit"] == 3.0
        assert result["ft3_lower_limit"] == 2.8
        # free_t3 2.5 < 2.8, so should flag subclinical
        assert result["stage"] == "subclinical"


# =====================================================================
# HFE FERRITIN ADJUSTMENT
# =====================================================================


class TestHFEAdjustment:
    """Tests for HFE-adjusted ferritin reference ranges."""

    def test_hfe_gg_standard_ferritin_male(self, analyzer):
        """HFE GG male should use standard ferritin upper limit of 400."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 350},
            genotypes={"HFE_rs1800562": "GG"},
            sex="male",
        )
        assert result["ferritin_upper_limit"] == 400
        assert result["stage"] == "normal"

    def test_hfe_ga_reduced_ferritin_male(self, analyzer):
        """HFE GA male should use reduced ferritin upper limit of 300."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 320},
            genotypes={"HFE_rs1800562": "GA"},
            sex="male",
        )
        assert result["ferritin_upper_limit"] == 300
        assert result["stage"] == "early_accumulation"

    def test_hfe_aa_strict_ferritin_female(self, analyzer):
        """HFE AA female should use strict ferritin upper limit of 100."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 120},
            genotypes={"HFE_rs1800562": "AA"},
            sex="female",
        )
        assert result["ferritin_upper_limit"] == 100


# =====================================================================
# MISSING GENOTYPES
# =====================================================================


class TestMissingGenotypes:
    """Tests that missing genotypes produce no adjustment (standard ranges)."""

    def test_no_genotype_standard_alt(self, analyzer):
        """Missing PNPLA3 genotype should use standard ALT upper limit."""
        result = analyzer.analyze_liver(
            biomarkers={"alt": 50},
            genotypes={},
        )
        assert result["alt_upper_limit"] == 56

    def test_no_genotype_standard_hba1c(self, analyzer):
        """Missing TCF7L2 genotype should use standard HbA1c threshold."""
        result = analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 5.6},
            genotypes={},
        )
        assert result["genotype_adjusted_thresholds"]["hba1c"] == 6.0

    def test_no_genotype_standard_ldl(self, analyzer):
        """Missing APOE genotype should use standard LDL threshold of 130."""
        result = analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 120},
            genotypes={},
        )
        assert result["ldl_threshold"] == 130

    def test_no_genotype_standard_tsh(self, analyzer):
        """Missing DIO2 genotype should use standard TSH upper limit."""
        result = analyzer.analyze_thyroid(
            biomarkers={"tsh": 3.5},
            genotypes={},
        )
        assert result["tsh_upper_limit"] == 4.0

    def test_no_genotype_standard_ferritin(self, analyzer):
        """Missing HFE genotype should use standard ferritin upper limit."""
        result = analyzer.analyze_iron(
            biomarkers={"ferritin": 350},
            genotypes={},
            sex="male",
        )
        assert result["ferritin_upper_limit"] == 400


# =====================================================================
# GENOTYPE ADJUSTMENT MODEL
# =====================================================================


class TestGenotypeAdjustmentModel:
    """Tests for the GenotypeAdjustment Pydantic model."""

    def test_create_genotype_adjustment(self):
        """GenotypeAdjustment model can be created with required fields."""
        adj = GenotypeAdjustment(
            id="adj-pnpla3-alt",
            biomarker="ALT",
            gene="PNPLA3",
            rs_id="rs738409",
            genotype_ref="CC",
            genotype_het="CG",
            genotype_hom="GG",
            adjusted_min=0.0,
            adjusted_max=35.0,
            text_chunk="PNPLA3 GG lowers ALT upper reference to 35 U/L.",
            rationale="PNPLA3 I148M causes hepatic lipid accumulation.",
        )
        assert adj.biomarker == "ALT"
        assert adj.adjusted_max == 35.0

    def test_genotype_adjustment_result_model(self):
        """GenotypeAdjustmentResult can be created and serialized."""
        result = GenotypeAdjustmentResult(
            biomarker="ALT",
            standard_range="0-56 U/L",
            adjusted_range="0-35 U/L",
            genotype="GG",
            gene="PNPLA3",
            rationale="PNPLA3 GG carrier: 2-3x NAFLD risk",
        )
        assert result.biomarker == "ALT"
        assert result.adjusted_range == "0-35 U/L"
        data = result.model_dump()
        assert data["gene"] == "PNPLA3"
