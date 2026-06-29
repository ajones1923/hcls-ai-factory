"""Tests for Biomarker Intelligence Agent pharmacogenomic mapping.

Validates PGx phenotype determination, drug-gene interaction mapping,
star allele interpretation, and critical alert generation.

All tests delegate to the real PharmacogenomicMapper from
src.pharmacogenomics so that drug recommendations stay in sync
with the production module rather than a duplicated local copy.

Author: Adam Jones
Date: March 2026
"""


from src.models import MetabolizerPhenotype, PGxResult
from src.pharmacogenomics import PharmacogenomicMapper

# =====================================================================
# Shared mapper instance used by all test classes
# =====================================================================

_mapper = PharmacogenomicMapper()


# =====================================================================
# CYP2D6 TESTS
# =====================================================================


class TestCYP2D6:
    """Tests for CYP2D6 pharmacogenomic mapping."""

    def test_poor_metabolizer_codeine_avoid(self):
        """CYP2D6 *4/*4 + codeine should recommend AVOID."""
        result = _mapper.map_gene("CYP2D6", star_alleles="*4/*4")
        assert result["phenotype"] == "Poor Metabolizer"
        codeine = next(d for d in result["affected_drugs"] if d["drug"] == "codeine")
        assert "AVOID" in codeine["recommendation"]
        assert codeine["action"] == "AVOID"

    def test_intermediate_metabolizer(self):
        """CYP2D6 *1/*4 should yield Intermediate Metabolizer."""
        result = _mapper.map_gene("CYP2D6", star_alleles="*1/*4")
        assert result["phenotype"] == "Intermediate Metabolizer"

    def test_normal_metabolizer_standard_dosing(self):
        """CYP2D6 *1/*1 should yield STANDARD_DOSING for all drugs."""
        result = _mapper.map_gene("CYP2D6", star_alleles="*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        for d in result["affected_drugs"]:
            assert d["action"] == "STANDARD_DOSING"

    def test_ultra_rapid_metabolizer_codeine_avoid(self):
        """CYP2D6 *1/*1xN should yield Ultra-rapid and AVOID codeine."""
        result = _mapper.map_gene("CYP2D6", star_alleles="*1/*1xN")
        assert result["phenotype"] == "Ultra-rapid Metabolizer"
        codeine = next(d for d in result["affected_drugs"] if d["drug"] == "codeine")
        assert "AVOID" in codeine["recommendation"]
        assert codeine["action"] == "AVOID"


# =====================================================================
# CYP2C19 TESTS
# =====================================================================


class TestCYP2C19:
    """Tests for CYP2C19 pharmacogenomic mapping."""

    def test_poor_metabolizer_clopidogrel_avoid(self):
        """CYP2C19 *2/*2 + clopidogrel should recommend AVOID."""
        result = _mapper.map_gene("CYP2C19", star_alleles="*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"
        clopidogrel = next(d for d in result["affected_drugs"] if d["drug"] == "clopidogrel")
        assert "AVOID" in clopidogrel["recommendation"]
        assert "prasugrel" in clopidogrel["recommendation"]

    def test_poor_metabolizer_omeprazole_dose_reduction(self):
        """CYP2C19 *2/*2 + omeprazole should recommend dose reduction."""
        result = _mapper.map_gene("CYP2C19", star_alleles="*2/*2")
        omeprazole = next(d for d in result["affected_drugs"] if d["drug"] == "omeprazole")
        assert "reduction" in omeprazole["recommendation"].lower()
        assert omeprazole["action"] == "DOSE_REDUCTION"

    def test_normal_metabolizer_standard(self):
        """CYP2C19 *1/*1 should yield standard dosing."""
        result = _mapper.map_gene("CYP2C19", star_alleles="*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        for d in result["affected_drugs"]:
            assert d["action"] == "STANDARD_DOSING"


# =====================================================================
# SLCO1B1 TESTS
# =====================================================================


class TestSLCO1B1:
    """Tests for SLCO1B1 pharmacogenomic mapping."""

    def test_poor_function_simvastatin_avoid(self):
        """SLCO1B1 CC genotype + simvastatin should recommend AVOID."""
        result = _mapper.map_gene("SLCO1B1", genotype="CC")
        assert result["phenotype"] == "Poor Function"
        simva = next(d for d in result["affected_drugs"] if d["drug"] == "simvastatin")
        assert "AVOID" in simva["recommendation"]
        assert simva["action"] == "AVOID"

    def test_normal_function_standard(self):
        """SLCO1B1 TT genotype should yield Normal Function."""
        result = _mapper.map_gene("SLCO1B1", genotype="TT")
        assert result["phenotype"] == "Normal Function"


# =====================================================================
# MTHFR TESTS
# =====================================================================


class TestMTHFR:
    """Tests for MTHFR pharmacogenomic mapping."""

    def test_mthfr_tt_reduced_activity(self):
        """MTHFR TT genotype should map to Reduced Activity."""
        result = _mapper.map_gene("MTHFR", genotype="TT")
        assert result["phenotype"] == "Reduced Activity"

    def test_mthfr_reduced_activity_methotrexate(self):
        """MTHFR Reduced Activity + methotrexate should recommend dose reduction."""
        result = _mapper.map_gene("MTHFR", genotype="TT")
        mtx = next(d for d in result["affected_drugs"] if d["drug"] == "methotrexate")
        assert "reduc" in mtx["recommendation"].lower()
        assert mtx["action"] in ("DOSE_REDUCTION", "AVOID")

    def test_mthfr_normal_standard(self):
        """MTHFR CC genotype should map to Normal Activity."""
        result = _mapper.map_gene("MTHFR", genotype="CC")
        assert result["phenotype"] == "Normal Activity"


# =====================================================================
# HLA-B*57:01 TESTS
# =====================================================================


class TestHLAB5701:
    """Tests for HLA-B*57:01 abacavir contraindication."""

    def test_hla_positive_abacavir_contraindicated(self):
        """HLA-B*57:01 positive + abacavir should be CONTRAINDICATED."""
        result = _mapper.map_gene("HLA-B*57:01", genotype="positive")
        assert result["phenotype"] == "Positive"
        abacavir = next(d for d in result["affected_drugs"] if d["drug"] == "abacavir")
        assert "CONTRAINDICATED" in abacavir["recommendation"].upper() or abacavir["action"] == "CONTRAINDICATED"

    def test_hla_negative_abacavir_safe(self):
        """HLA-B*57:01 negative should not trigger abacavir critical alert."""
        result = _mapper.map_gene("HLA-B*57:01", genotype="negative")
        assert result["phenotype"] == "Negative"
        assert len(result["critical_alerts"]) == 0


# =====================================================================
# TPMT TESTS
# =====================================================================


class TestTPMT:
    """Tests for TPMT pharmacogenomic mapping."""

    def test_tpmt_poor_azathioprine_avoid(self):
        """TPMT *3A/*3A + azathioprine should recommend AVOID."""
        result = _mapper.map_gene("TPMT", star_alleles="*3A/*3A")
        assert result["phenotype"] == "Poor Metabolizer"
        aza = next(d for d in result["affected_drugs"] if d["drug"] == "azathioprine")
        assert "AVOID" in aza["recommendation"]
        assert aza["action"] == "AVOID"

    def test_tpmt_normal_standard(self):
        """TPMT *1/*1 should yield standard dosing."""
        result = _mapper.map_gene("TPMT", star_alleles="*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        for d in result["affected_drugs"]:
            assert d["action"] == "STANDARD_DOSING"


# =====================================================================
# MULTIPLE GENES SIMULTANEOUSLY
# =====================================================================


class TestMultipleGenes:
    """Tests for multiple PGx genes evaluated simultaneously via map_all()."""

    def test_multiple_genes_all_normal(self):
        """All normal metabolizers should yield STANDARD_DOSING across all genes."""
        result = _mapper.map_all(
            star_alleles={
                "CYP2D6": "*1/*1",
                "CYP2C19": "*1/*1",
                "TPMT": "*1/*1",
            },
        )
        for gene_result in result["gene_results"]:
            assert "Normal" in gene_result["phenotype"]
            for d in gene_result["affected_drugs"]:
                assert d["action"] == "STANDARD_DOSING"

    def test_multiple_poor_metabolizers_multiple_alerts(self):
        """Multiple poor metabolizer results should generate multiple AVOID actions."""
        result = _mapper.map_all(
            star_alleles={
                "CYP2D6": "*4/*4",
                "CYP2C19": "*2/*2",
            },
        )
        avoid_count = sum(
            1 for d in result["drugs_to_avoid"]
        )
        # CYP2D6 PM: codeine, tramadol, tamoxifen AVOID; CYP2C19 PM: clopidogrel, voriconazole AVOID
        assert avoid_count >= 3

    def test_mixed_metabolizer_phenotypes(self):
        """Different genes can have different phenotypes simultaneously."""
        result = _mapper.map_all(
            star_alleles={
                "CYP2D6": "*4/*4",   # Poor Metabolizer
                "CYP2C19": "*1/*1",  # Normal Metabolizer
                "TPMT": "*1/*3A",    # Intermediate Metabolizer
            },
        )
        phenotypes = {r["gene"]: r["phenotype"] for r in result["gene_results"]}
        assert phenotypes["CYP2D6"] == "Poor Metabolizer"
        assert phenotypes["CYP2C19"] == "Normal Metabolizer"
        assert phenotypes["TPMT"] == "Intermediate Metabolizer"


# =====================================================================
# PGxResult MODEL
# =====================================================================


class TestPGxResultModel:
    """Tests for the PGxResult Pydantic model."""

    def test_create_pgx_result(self):
        """PGxResult can be created with required fields."""
        result = PGxResult(
            gene="CYP2D6",
            star_alleles="*4/*4",
            phenotype=MetabolizerPhenotype.POOR,
            drugs_affected=[
                {"drug": "Codeine", "recommendation": "AVOID", "cpic_level": "1A"},
            ],
        )
        assert result.gene == "CYP2D6"
        assert result.phenotype == MetabolizerPhenotype.POOR
        assert len(result.drugs_affected) == 1

    def test_default_phenotype_is_normal(self):
        """PGxResult defaults to NORMAL phenotype."""
        result = PGxResult(gene="CYP2D6", star_alleles="*1/*1")
        assert result.phenotype == MetabolizerPhenotype.NORMAL
        assert result.drugs_affected == []


# =====================================================================
# DRUG-DRUG INTERACTION SCREENING
# =====================================================================


class TestDrugInteractions:
    """Tests for cross-recommendation drug-drug interaction detection."""

    def test_no_interactions_when_no_overlapping_drugs(self):
        """Non-overlapping drug sets should produce zero interaction warnings."""
        # CYP2D6 normal (codeine, tramadol, etc.) and TPMT normal (azathioprine, etc.)
        # All standard dosing, no interacting pairs in the table.
        result = _mapper.map_all(
            star_alleles={
                "TPMT": "*1/*1",
            },
            genotypes={
                "SLCO1B1_rs4149056": "TT",
            },
        )
        interactions = result["drug_interactions"]
        assert interactions == []

    def test_detects_clopidogrel_omeprazole_interaction(self):
        """Clopidogrel (CYP2C19) + omeprazole (CYP2C19) should flag interaction."""
        # CYP2C19 *1/*1 covers both clopidogrel and omeprazole
        result = _mapper.map_all(
            star_alleles={
                "CYP2C19": "*1/*1",
            },
        )
        interactions = result["drug_interactions"]
        pair_found = any(
            (i["drug_a"] == "clopidogrel" and i["drug_b"] == "omeprazole")
            or (i["drug_a"] == "omeprazole" and i["drug_b"] == "clopidogrel")
            for i in interactions
        )
        assert pair_found, (
            "Expected clopidogrel-omeprazole interaction warning but got: "
            f"{interactions}"
        )
        interaction = next(
            i for i in interactions
            if {i["drug_a"], i["drug_b"]} == {"clopidogrel", "omeprazole"}
        )
        assert interaction["severity"] == "high"
        assert "CYP2C19" in interaction["genes_involved"]

    def test_detects_fluoxetine_codeine_interaction(self):
        """Fluoxetine (CYP2D6 inhibitor) + codeine (CYP2D6 substrate) should flag."""
        # We need both drugs in the results. CYP2D6 covers codeine and tramadol;
        # fluoxetine is not directly in PGx gene configs, so we build a synthetic
        # map_all result to exercise check_drug_interactions directly.
        synthetic_result = {
            "gene_results": [
                {
                    "gene": "CYP2D6",
                    "phenotype": "Normal Metabolizer",
                    "affected_drugs": [
                        {"drug": "codeine", "action": "STANDARD_DOSING", "alert_level": "INFO"},
                        {"drug": "fluoxetine", "action": "STANDARD_DOSING", "alert_level": "INFO"},
                    ],
                    "critical_alerts": [],
                },
            ],
        }
        interactions = _mapper.check_drug_interactions(synthetic_result)
        pair_found = any(
            {i["drug_a"], i["drug_b"]} == {"fluoxetine", "codeine"}
            for i in interactions
        )
        assert pair_found, (
            "Expected fluoxetine-codeine interaction warning but got: "
            f"{interactions}"
        )
        interaction = next(
            i for i in interactions
            if {i["drug_a"], i["drug_b"]} == {"fluoxetine", "codeine"}
        )
        assert interaction["severity"] == "high"
        assert "CYP2D6" in interaction["mechanism"]


# =====================================================================
# GUIDELINE VERSION TRACKING
# =====================================================================


class TestGuidelineVersions:
    """Test CPIC guideline version tracking."""

    def test_all_genes_have_versions(self):
        from src.pharmacogenomics import CPIC_GUIDELINE_VERSIONS, PharmacogenomicMapper
        PharmacogenomicMapper()
        # Every gene in PGX_GENE_CONFIGS should have a version entry
        from src.pharmacogenomics import PGX_GENE_CONFIGS
        for gene in PGX_GENE_CONFIGS:
            assert gene in CPIC_GUIDELINE_VERSIONS, f"Missing version for {gene}"

    def test_get_guideline_versions_returns_dict(self):
        from src.pharmacogenomics import PharmacogenomicMapper
        mapper = PharmacogenomicMapper()
        versions = mapper.get_guideline_versions()
        assert isinstance(versions, dict)
        assert "CYP2D6" in versions
        assert "level" in versions["CYP2D6"]

    def test_map_all_includes_guideline_versions(self):
        from src.pharmacogenomics import PharmacogenomicMapper
        mapper = PharmacogenomicMapper()
        result = mapper.map_all(star_alleles={"CYP2D6": "*1/*1"})
        assert "guideline_versions" in result
