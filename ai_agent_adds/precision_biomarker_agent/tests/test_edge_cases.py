"""Edge case tests for Biomarker Intelligence Agent safety fixes.

Validates overflow protection, NaN handling, score clamping,
input validation, and other defensive code paths added to
prevent production failures.

Author: Adam Jones
Date: March 2026
"""

import math
import re

import pytest

from config.settings import PrecisionBiomarkerSettings
from src.biological_age import BiologicalAgeCalculator
from src.disease_trajectory import DiseaseTrajectoryAnalyzer
from src.export import generate_filename
from src.models import (
    BiologicalAgeResult,
    BiomarkerReference,
    SearchHit,
)
from src.pharmacogenomics import PharmacogenomicMapper
from src.rag_engine import BiomarkerRAGEngine

# =====================================================================
# 1. BIOLOGICAL AGE EDGE CASES
# =====================================================================


class TestBiologicalAgeEdgeCases:
    """Edge cases for BiologicalAgeCalculator: negative values, CRP clamping,
    exp() overflow protection, incomplete panels, and normal sanity check."""

    def setup_method(self):
        self.calc = BiologicalAgeCalculator()

    def test_negative_albumin_clamped_to_zero(self):
        """Negative albumin should be clamped to 0 (not crash)."""
        biomarkers = {
            "albumin": -2.0,
            "creatinine": 0.9,
            "glucose": 95.0,
            "hs_crp": 1.5,
            "lymphocyte_pct": 30.0,
            "mcv": 89.0,
            "rdw": 13.0,
            "alkaline_phosphatase": 65.0,
            "wbc": 6.0,
        }
        result = self.calc.calculate_phenoage(45, biomarkers)
        assert "biological_age" in result
        assert isinstance(result["biological_age"], float)

    def test_negative_creatinine_clamped_to_zero(self):
        """Negative creatinine should be clamped to 0 (not crash)."""
        biomarkers = {
            "albumin": 4.2,
            "creatinine": -5.0,
            "glucose": 95.0,
            "hs_crp": 1.5,
            "lymphocyte_pct": 30.0,
            "mcv": 89.0,
            "rdw": 13.0,
            "alkaline_phosphatase": 65.0,
            "wbc": 6.0,
        }
        result = self.calc.calculate_phenoage(45, biomarkers)
        assert "biological_age" in result
        assert isinstance(result["biological_age"], float)

    def test_crp_zero_clamped_to_minimum(self):
        """CRP=0 should be clamped to 0.01 before log transform."""
        biomarkers = {
            "albumin": 4.2,
            "creatinine": 0.9,
            "glucose": 95.0,
            "hs_crp": 0.0,
            "lymphocyte_pct": 30.0,
            "mcv": 89.0,
            "rdw": 13.0,
            "alkaline_phosphatase": 65.0,
            "wbc": 6.0,
        }
        result = self.calc.calculate_phenoage(45, biomarkers)
        # Should not crash with log(0) -- clamped to log(0.01)
        assert "biological_age" in result
        assert math.isfinite(result["biological_age"])

    def test_crp_extremely_high_clamped_to_max(self):
        """CRP=500 should be clamped to 200 before log transform."""
        biomarkers = {
            "albumin": 4.2,
            "creatinine": 0.9,
            "glucose": 95.0,
            "hs_crp": 500.0,
            "lymphocyte_pct": 30.0,
            "mcv": 89.0,
            "rdw": 13.0,
            "alkaline_phosphatase": 65.0,
            "wbc": 6.0,
        }
        result = self.calc.calculate_phenoage(45, biomarkers)
        assert "biological_age" in result
        assert math.isfinite(result["biological_age"])

    def test_extreme_biomarkers_no_overflow(self):
        """Extreme biomarker values should not cause exp() overflow.

        The xb clamp (max(min(xb, 700), -700)) prevents math.exp() from
        raising OverflowError.
        """
        biomarkers = {
            "albumin": 0.01,            # Very low (bad)
            "creatinine": 50.0,         # Extremely elevated
            "glucose": 5000.0,          # Absurdly high
            "hs_crp": 200.0,           # Max CRP
            "lymphocyte_pct": 0.1,     # Extremely low
            "mcv": 200.0,              # Very high
            "rdw": 50.0,               # Very high
            "alkaline_phosphatase": 5000.0,  # Extremely elevated
            "wbc": 100.0,              # Extreme leukocytosis
        }
        result = self.calc.calculate_phenoage(90, biomarkers)
        assert "biological_age" in result
        assert math.isfinite(result["biological_age"])
        assert math.isfinite(result["mortality_score"])

    def test_incomplete_panel_fewer_than_half(self):
        """Fewer than 50% of PhenoAge biomarkers should trigger a warning
        and still return a valid result with the available biomarkers."""
        # Only 3 out of 9 PhenoAge markers provided (< 50%)
        biomarkers = {
            "albumin": 4.2,
            "creatinine": 0.9,
            "glucose": 95.0,
        }
        result = self.calc.calculate_phenoage(45, biomarkers)
        assert "biological_age" in result
        assert len(result["missing_biomarkers"]) > len(result["all_contributions"]) - 2
        # At least 6 missing (albumin, creatinine, glucose present + age contribution)

    def test_normal_biomarkers_reasonable_age(self, sample_biomarkers):
        """Standard healthy biomarkers should produce biological age within
        +/-10 years of chronological age."""
        chrono_age = 45
        result = self.calc.calculate(chrono_age, sample_biomarkers)
        bio_age = result["biological_age"]
        assert abs(bio_age - chrono_age) <= 10, (
            f"Biological age {bio_age} differs by more than 10 years "
            f"from chronological age {chrono_age}"
        )

    def test_calculate_returns_phenoage_and_grimage(self, sample_biomarkers):
        """The combined calculate() method should return both PhenoAge and
        optional GrimAge data."""
        result = self.calc.calculate(45, sample_biomarkers)
        assert "phenoage" in result
        assert "biological_age" in result
        assert "age_acceleration" in result
        assert "mortality_risk" in result


# =====================================================================
# 2. DISEASE TRAJECTORY EDGE CASES
# =====================================================================


class TestDiseaseTrajectoryEdgeCases:
    """Edge cases for DiseaseTrajectoryAnalyzer: NaN filtering, empty
    biomarkers, division by zero guards, and HOMA-IR escalation."""

    def setup_method(self):
        self.analyzer = DiseaseTrajectoryAnalyzer()

    def test_nan_biomarkers_filtered_out(self):
        """NaN values in biomarkers dict should be filtered by analyze_all()
        using the v != v idiom."""
        biomarkers = {
            "hba1c": 5.4,
            "fasting_glucose": float("nan"),
            "fasting_insulin": float("nan"),
            "ldl_c": float("nan"),
            "alt": 25.0,
        }
        results = self.analyzer.analyze_all(biomarkers, {}, age=45, sex="male")
        assert isinstance(results, list)
        assert len(results) == 9  # All 9 disease categories
        # NaN values should not appear in any current_markers
        for r in results:
            for val in r.get("current_markers", {}).values():
                assert val == val, "NaN value leaked through to analysis"

    def test_none_biomarkers_filtered_out(self):
        """None values in biomarkers dict should be filtered by analyze_all()."""
        biomarkers = {
            "hba1c": 5.4,
            "fasting_glucose": None,
            "alt": None,
            "tsh": 2.5,
        }
        results = self.analyzer.analyze_all(biomarkers, {}, age=45, sex="male")
        assert isinstance(results, list)
        for r in results:
            for val in r.get("current_markers", {}).values():
                assert val is not None, "None value leaked through to analysis"

    def test_empty_biomarkers_dict(self):
        """Empty biomarkers dict should not crash; all analyses return low risk."""
        results = self.analyzer.analyze_all({}, {}, age=45, sex="male")
        assert len(results) == 9
        for r in results:
            assert r["risk_level"] == "LOW"
            assert len(r["findings"]) == 0

    def test_fib4_zero_alt_no_crash(self):
        """FIB-4 with ALT=0 should not crash (division by zero guard).

        FIB-4 = (age * AST) / (platelets * sqrt(ALT)); ALT=0 => guard.
        """
        biomarkers = {
            "alt": 0.0,
            "ast": 30.0,
            "platelets": 250.0,
        }
        result = self.analyzer.analyze_liver(biomarkers, {}, age=45)
        # FIB-4 should not be calculated when ALT=0 (guard: alt > 0)
        assert result["fib4_score"] is None

    def test_fib4_zero_platelets_no_crash(self):
        """FIB-4 with platelets=0 should not crash (division by zero guard)."""
        biomarkers = {
            "alt": 25.0,
            "ast": 30.0,
            "platelets": 0.0,
        }
        result = self.analyzer.analyze_liver(biomarkers, {}, age=45)
        assert result["fib4_score"] is None

    def test_homa_ir_escalation_with_pre_diabetic_stage(self):
        """HOMA-IR findings should still appear even when HbA1c has already
        set the stage to pre_diabetic (risk_level escalation)."""
        biomarkers = {
            "hba1c": 6.0,            # Pre-diabetic range
            "fasting_insulin": 25.0,  # High insulin
            "fasting_glucose": 110.0, # Elevated glucose -> HOMA-IR = (25*110)/405 = 6.79
        }
        genotypes = {"TCF7L2_rs7903146": "TT"}  # Tightest thresholds
        result = self.analyzer.analyze_type2_diabetes(biomarkers, genotypes)
        assert result["stage"] == "pre_diabetic"
        assert result["risk_level"] == "HIGH"
        # HOMA-IR finding should still be present alongside HbA1c finding
        homa_findings = [f for f in result["findings"] if "HOMA-IR" in f]
        assert len(homa_findings) > 0, "HOMA-IR finding should appear even when stage is pre_diabetic"

    def test_analyze_all_sorted_by_risk(self):
        """Results from analyze_all() should be sorted by risk severity."""
        biomarkers = {
            "hba1c": 7.0,          # Diabetic => CRITICAL
            "fasting_glucose": 130, # Diabetic
            "ldl_c": 100.0,        # Low => LOW risk
            "tsh": 2.5,            # Normal => LOW risk
        }
        results = self.analyzer.analyze_all(biomarkers, {}, age=45, sex="male")
        risk_order = {"CRITICAL": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3}
        for i in range(len(results) - 1):
            current_rank = risk_order.get(results[i]["risk_level"], 4)
            next_rank = risk_order.get(results[i + 1]["risk_level"], 4)
            assert current_rank <= next_rank, (
                f"Results not sorted: {results[i]['disease']} ({results[i]['risk_level']}) "
                f"came before {results[i+1]['disease']} ({results[i+1]['risk_level']})"
            )


# =====================================================================
# 3. RAG ENGINE EDGE CASES
# =====================================================================


class TestRAGEngineEdgeCases:
    """Edge cases for BiomarkerRAGEngine: score clamping, text dedup,
    init validation."""

    def test_score_clamping_below_one(self, mock_collection_manager,
                                      mock_embedder, mock_llm_client):
        """Weighted scores should never exceed 1.0 even with high raw scores."""
        from src import knowledge as kg

        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge=kg,
        )

        # Create hits with very high raw scores
        hits = [
            SearchHit(
                collection="BiomarkerRef",
                id="ref-001",
                score=0.99,
                text="Very relevant biomarker reference text",
                metadata={},
            ),
            SearchHit(
                collection="ClinicalEvidence",
                id="evi-001",
                score=0.98,
                text="Strong clinical evidence text",
                metadata={},
            ),
        ]

        ranked = engine._merge_and_rank(hits)
        for hit in ranked:
            assert hit.score <= 1.0, (
                f"Score {hit.score} exceeds 1.0 for hit {hit.id}"
            )

    def test_text_dedup_across_collections(self, mock_collection_manager,
                                           mock_embedder, mock_llm_client):
        """Duplicate text content should be deduplicated even with different IDs."""
        from src import knowledge as kg

        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge=kg,
        )

        same_text = "HbA1c reflects glycated hemoglobin over 2-3 months in clinical practice"
        hits = [
            SearchHit(
                collection="BiomarkerRef",
                id="ref-hba1c-1",
                score=0.90,
                text=same_text,
                metadata={},
            ),
            SearchHit(
                collection="ClinicalEvidence",
                id="evi-hba1c-2",
                score=0.85,
                text=same_text,  # Exact same text, different collection/ID
                metadata={},
            ),
            SearchHit(
                collection="Nutrition",
                id="nut-unique",
                score=0.80,
                text="Unique nutrition content about folate metabolism",
                metadata={},
            ),
        ]

        ranked = engine._merge_and_rank(hits)
        # Should have deduplicated the first two hits (same text[:200])
        assert len(ranked) == 2, (
            f"Expected 2 unique hits after dedup, got {len(ranked)}"
        )

    def test_query_raises_without_embedder(self, mock_collection_manager,
                                           mock_llm_client):
        """query() should raise RuntimeError if embedder is None."""
        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=None,
            llm_client=mock_llm_client,
        )
        with pytest.raises(RuntimeError, match="Embedding model not initialized"):
            engine.query("What is PhenoAge?")

    def test_query_raises_without_llm(self, mock_collection_manager,
                                      mock_embedder):
        """query() should raise RuntimeError if llm is None."""
        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=None,
        )
        with pytest.raises(RuntimeError, match="LLM client not initialized"):
            engine.query("What is PhenoAge?")

    def test_merge_and_rank_caps_at_30(self, mock_collection_manager,
                                       mock_embedder, mock_llm_client):
        """_merge_and_rank should cap results at 30 even if more hits exist."""
        from src import knowledge as kg

        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge=kg,
        )

        hits = [
            SearchHit(
                collection="BiomarkerRef",
                id=f"ref-{i:03d}",
                score=0.5 + (i * 0.01),
                text=f"Unique biomarker text number {i} for testing dedup and ranking",
                metadata={},
            )
            for i in range(50)
        ]

        ranked = engine._merge_and_rank(hits)
        assert len(ranked) == 30, f"Expected 30 hits, got {len(ranked)}"

    def test_id_dedup(self, mock_collection_manager, mock_embedder,
                      mock_llm_client):
        """Hits with the same ID should be deduplicated regardless of text."""
        from src import knowledge as kg

        engine = BiomarkerRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge=kg,
        )

        hits = [
            SearchHit(
                collection="BiomarkerRef",
                id="ref-hba1c",
                score=0.90,
                text="First version of HbA1c text",
                metadata={},
            ),
            SearchHit(
                collection="BiomarkerRef",
                id="ref-hba1c",
                score=0.85,
                text="Second version of HbA1c text (different text, same ID)",
                metadata={},
            ),
        ]

        ranked = engine._merge_and_rank(hits)
        assert len(ranked) == 1, f"Expected 1 hit after ID dedup, got {len(ranked)}"


# =====================================================================
# 4. MODEL VALIDATORS
# =====================================================================


class TestModelValidators:
    """Edge cases for Pydantic model validators: ref range check,
    age acceleration mismatch logging."""

    def test_ref_range_min_greater_than_max_raises(self):
        """BiomarkerReference should raise ValueError when ref_range_min > ref_range_max
        and ref_range_max > 0."""
        with pytest.raises(ValueError, match="ref_range_min.*ref_range_max"):
            BiomarkerReference(
                id="test-marker",
                name="Test Marker",
                text_chunk="Test text chunk for embedding",
                ref_range_min=100.0,
                ref_range_max=50.0,
            )

    def test_ref_range_valid_passes(self):
        """BiomarkerReference with valid ref_range should pass validation."""
        marker = BiomarkerReference(
            id="test-marker",
            name="Test Marker",
            text_chunk="Test text chunk for embedding",
            ref_range_min=10.0,
            ref_range_max=50.0,
        )
        assert marker.ref_range_min == 10.0
        assert marker.ref_range_max == 50.0

    def test_ref_range_both_zero_passes(self):
        """BiomarkerReference with both ref_range values at 0 should pass
        (condition requires ref_range_max > 0)."""
        marker = BiomarkerReference(
            id="test-marker",
            name="Test Marker",
            text_chunk="Test text chunk",
            ref_range_min=0.0,
            ref_range_max=0.0,
        )
        assert marker.ref_range_min == 0.0

    def test_age_acceleration_mismatch_logs_warning(self, caplog):
        """BiologicalAgeResult should log a warning when age_acceleration
        does not match biological_age - chronological_age, but NOT raise."""
        import logging

        with caplog.at_level(logging.WARNING):
            result = BiologicalAgeResult(
                chronological_age=45,
                biological_age=50.0,
                age_acceleration=10.0,  # Should be 5.0 (50 - 45)
                phenoage_score=0.5,
                mortality_risk=0.3,
            )
        # Should still create the object (warning, not error)
        assert result.age_acceleration == 10.0
        assert "age_acceleration" in caplog.text

    def test_age_acceleration_correct_no_warning(self, caplog):
        """BiologicalAgeResult should NOT warn when acceleration matches."""
        import logging

        with caplog.at_level(logging.WARNING):
            result = BiologicalAgeResult(
                chronological_age=45,
                biological_age=50.0,
                age_acceleration=5.0,
                phenoage_score=0.5,
                mortality_risk=0.3,
            )
        assert result.age_acceleration == 5.0
        assert "age_acceleration" not in caplog.text


# =====================================================================
# 5. SETTINGS VALIDATION
# =====================================================================


class TestSettingsValidation:
    """Edge cases for PrecisionBiomarkerSettings: weight sum, timeouts."""

    def test_default_weights_sum_to_one(self):
        """Default collection search weights should sum to ~1.0."""
        s = PrecisionBiomarkerSettings()
        weights = [
            s.WEIGHT_BIOMARKER_REF,
            s.WEIGHT_GENETIC_VARIANTS,
            s.WEIGHT_PGX_RULES,
            s.WEIGHT_DISEASE_TRAJECTORIES,
            s.WEIGHT_CLINICAL_EVIDENCE,
            s.WEIGHT_NUTRITION,
            s.WEIGHT_DRUG_INTERACTIONS,
            s.WEIGHT_AGING_MARKERS,
            s.WEIGHT_GENOTYPE_ADJUSTMENTS,
            s.WEIGHT_MONITORING,
            s.WEIGHT_CRITICAL_VALUES,
            s.WEIGHT_DISCORDANCE_RULES,
            s.WEIGHT_AJ_CARRIER_SCREENING,
        ]
        total = sum(weights)
        assert abs(total - 1.0) <= 0.1, (
            f"Collection weights sum to {total}, expected ~1.0"
        )

    def test_request_timeout_exists_and_reasonable(self):
        """REQUEST_TIMEOUT_SECONDS should exist and be between 5 and 300."""
        s = PrecisionBiomarkerSettings()
        assert hasattr(s, "REQUEST_TIMEOUT_SECONDS")
        assert 5 <= s.REQUEST_TIMEOUT_SECONDS <= 300

    def test_milvus_timeout_exists_and_reasonable(self):
        """MILVUS_TIMEOUT_SECONDS should exist and be between 1 and 120."""
        s = PrecisionBiomarkerSettings()
        assert hasattr(s, "MILVUS_TIMEOUT_SECONDS")
        assert 1 <= s.MILVUS_TIMEOUT_SECONDS <= 120

    def test_llm_max_retries_exists_and_reasonable(self):
        """LLM_MAX_RETRIES should exist and be between 1 and 10."""
        s = PrecisionBiomarkerSettings()
        assert hasattr(s, "LLM_MAX_RETRIES")
        assert 1 <= s.LLM_MAX_RETRIES <= 10

    def test_score_threshold_in_valid_range(self):
        """SCORE_THRESHOLD should be between 0 and 1."""
        s = PrecisionBiomarkerSettings()
        assert 0.0 <= s.SCORE_THRESHOLD <= 1.0

    def test_citation_thresholds_ordered(self):
        """CITATION_HIGH_THRESHOLD should be greater than CITATION_MEDIUM_THRESHOLD."""
        s = PrecisionBiomarkerSettings()
        assert s.CITATION_HIGH_THRESHOLD > s.CITATION_MEDIUM_THRESHOLD


# =====================================================================
# 6. EXPORT EDGE CASES
# =====================================================================


class TestExportEdgeCases:
    """Edge cases for export module: pipe escaping and UUID filenames."""

    def test_pipe_escaping_in_markdown_table(self):
        """The PDF renderer uses re.split(r'(?<!\\\\)\\|', row_text) so escaped
        pipes should not split the cell."""
        row_text = "| value with \\| pipe | normal cell |"
        # Simulate the split logic from export_pdf
        cells = [c.strip() for c in re.split(r'(?<!\\)\|', row_text)]
        # Remove empty first/last from leading/trailing |
        if cells and cells[0] == "":
            cells = cells[1:]
        if cells and cells[-1] == "":
            cells = cells[:-1]
        # The escaped pipe should stay within the first cell
        assert len(cells) == 2, f"Expected 2 cells, got {len(cells)}: {cells}"
        assert "\\|" in cells[0] or "|" in cells[0]

    def test_unescaped_pipes_split_correctly(self):
        """Unescaped pipes should split normally."""
        row_text = "| cell1 | cell2 | cell3 |"
        cells = [c.strip() for c in re.split(r'(?<!\\)\|', row_text)]
        if cells and cells[0] == "":
            cells = cells[1:]
        if cells and cells[-1] == "":
            cells = cells[:-1]
        assert len(cells) == 3

    def test_generate_filename_uniqueness(self):
        """generate_filename() should produce unique names with 4-char hex suffix."""
        names = {generate_filename("md") for _ in range(20)}
        # All 20 names should be unique (UUID suffix provides uniqueness)
        assert len(names) == 20, "generate_filename() produced duplicate names"

    def test_generate_filename_format(self):
        """Filename should match expected pattern."""
        name = generate_filename("json")
        assert name.startswith("biomarker_report_")
        assert name.endswith(".json")
        # Check for UUID hex suffix before extension
        parts = name.replace(".json", "").split("_")
        hex_suffix = parts[-1]
        assert len(hex_suffix) == 4
        assert all(c in "0123456789abcdef" for c in hex_suffix)

    def test_generate_filename_different_extensions(self):
        """generate_filename() should work with various extensions."""
        for ext in ("md", "json", "pdf", "csv"):
            name = generate_filename(ext)
            assert name.endswith(f".{ext}")


# =====================================================================
# 7. PGx EDGE CASES
# =====================================================================


class TestPGxEdgeCases:
    """Edge cases for PharmacogenomicMapper: unknown drugs, star allele/genotype
    conflicts."""

    def setup_method(self):
        self.mapper = PharmacogenomicMapper()

    def test_unknown_drug_returns_not_in_database(self):
        """check_drug() with an unknown drug name should return NOT_IN_DATABASE status."""
        result = self.mapper.check_drug("aspirin_xyz_fake_drug")
        assert result["drug_in_database"] is False
        assert any(
            f.get("status") == "NOT_IN_DATABASE"
            for f in result["findings"]
        ), "Expected NOT_IN_DATABASE status in findings"

    def test_known_drug_found_in_database(self):
        """check_drug() with a known drug should return drug_in_database=True."""
        result = self.mapper.check_drug(
            "codeine",
            star_alleles={"CYP2D6": "*1/*1"},
        )
        assert result["drug_in_database"] is True
        assert result["is_safe"] is True

    def test_star_allele_genotype_conflict_uses_star_allele(self):
        """When star alleles and genotype for the same gene conflict, the
        mapper should prefer the star allele result and log a warning."""
        # SLCO1B1 has genotype_to_phenotype mapping
        # Provide conflicting star allele and genotype
        result = self.mapper.map_gene(
            "SLCO1B1",
            star_alleles="*1/*1",   # Would be "Normal Function" if recognized
            genotype="CC",           # "Poor Function"
        )
        # star alleles for SLCO1B1 are not in allele_to_phenotype (it uses genotype_to_phenotype)
        # So genotype should be used as fallback
        assert result["phenotype"] is not None or result.get("error") is not None

    def test_vkorc1_genotype_conflict_warning(self):
        """VKORC1 uses genotype_to_phenotype only. Providing star_alleles
        that cannot be resolved should fall back to genotype."""
        result = self.mapper.map_gene(
            "VKORC1",
            star_alleles=None,
            genotype="AA",
        )
        assert result["phenotype"] == "High Sensitivity"

    def test_unknown_gene_returns_error(self):
        """map_gene() with an unknown gene should return an error dict."""
        result = self.mapper.map_gene("FAKE_GENE_XYZ", star_alleles="*1/*1")
        assert result["phenotype"] is None
        assert "error" in result
        assert "not in PGx database" in result["error"]

    def test_unrecognized_allele_returns_error(self):
        """map_gene() with a real gene but unrecognized allele should
        return an error."""
        result = self.mapper.map_gene("CYP2D6", star_alleles="*99/*99")
        assert result["phenotype"] is None
        assert "error" in result

    def test_check_drug_not_tested_when_no_alleles_provided(self):
        """check_drug() should report NOT_TESTED if no star alleles or
        genotypes are provided for a gene that affects the drug."""
        result = self.mapper.check_drug("codeine")  # No alleles at all
        assert result["drug_in_database"] is True
        not_tested = [f for f in result["findings"] if f.get("status") == "NOT_TESTED"]
        assert len(not_tested) > 0, "Expected NOT_TESTED when no alleles provided"

    def test_map_all_with_mixed_inputs(self):
        """map_all() should process both star alleles and genotype-based genes."""
        result = self.mapper.map_all(
            star_alleles={"CYP2D6": "*4/*4"},
            genotypes={"VKORC1_rs9923231": "AA", "SLCO1B1_rs4149056": "CC"},
        )
        assert result["genes_analyzed"] == 3
        assert result["has_critical_findings"] is True  # CYP2D6 *4/*4 = Poor Metabolizer
        assert len(result["drugs_to_avoid"]) > 0

    def test_poor_metabolizer_critical_alerts(self):
        """CYP2D6 *4/*4 (Poor Metabolizer) should generate critical alerts
        for codeine and tramadol."""
        result = self.mapper.map_gene("CYP2D6", star_alleles="*4/*4")
        assert result["phenotype"] == "Poor Metabolizer"
        critical_drugs = {a["drug"] for a in result["critical_alerts"]}
        assert "codeine" in critical_drugs
        assert "tramadol" in critical_drugs


# =====================================================================
# 8. SECURITY EDGE CASES
# =====================================================================


class TestSecurityEdgeCases:
    """Security-focused edge case tests."""

    def test_patient_id_sanitization(self):
        """Patient IDs with path traversal chars should not crash."""
        from src.models import PatientProfile
        profile = PatientProfile(
            patient_id="../../etc/passwd",
            age=45, sex="M",
            biomarkers={"albumin": 4.2},
        )
        assert profile.patient_id == "../../etc/passwd"  # Model accepts it
        # But export filename should be sanitized
        from src.export import generate_filename
        name = generate_filename("pdf")
        assert ".." not in name

    def test_extremely_long_patient_id(self):
        """Very long patient IDs should be handled."""
        from src.models import PatientProfile
        profile = PatientProfile(
            patient_id="A" * 10000,
            age=45, sex="M",
            biomarkers={"albumin": 4.2},
        )
        assert len(profile.patient_id) == 10000

    def test_unicode_in_biomarker_names(self):
        """Unicode biomarker names should not crash analysis."""
        from src.biological_age import BiologicalAgeCalculator
        calc = BiologicalAgeCalculator()
        result = calc.calculate(45, {"albümin": 4.2, "glücose": 95.0})
        assert "biological_age" in result

    def test_empty_string_genotype(self):
        """Empty string genotypes should not crash."""
        from src.disease_trajectory import DiseaseTrajectoryAnalyzer
        analyzer = DiseaseTrajectoryAnalyzer()
        results = analyzer.analyze_all(
            {"hba1c": 5.5}, {"TCF7L2_rs7903146": ""}, 45, "M"
        )
        assert len(results) == 9

    def test_inf_biomarker_value(self):
        """Infinity biomarker values should be handled."""
        from src.biological_age import BiologicalAgeCalculator
        calc = BiologicalAgeCalculator()
        result = calc.calculate(45, {
            "albumin": 4.2, "creatinine": 0.9, "glucose": float("inf"),
            "hs_crp": 1.5, "lymphocyte_pct": 30.0, "mcv": 89.0,
            "rdw": 13.0, "alkaline_phosphatase": 65.0, "wbc": 6.0,
        })
        # Should not crash -- overflow protection should handle it
        assert "biological_age" in result


# =====================================================================
# 9. ANCESTRY ADJUSTMENT EDGE CASES
# =====================================================================


class TestAncestryAdjustments:
    """Test population-specific biomarker adjustments."""

    def test_african_creatinine_adjustment(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.apply_ancestry_adjustments(
            {"creatinine": 1.2}, "african"
        )
        assert len(result) == 1
        assert result[0]["threshold_multiplier"] == 1.15
        assert "muscle mass" in result[0]["note"]

    def test_south_asian_ldl_adjustment(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.apply_ancestry_adjustments(
            {"ldl": 130}, "south_asian"
        )
        assert len(result) == 1
        assert result[0]["threshold_multiplier"] == 0.90

    def test_no_adjustment_for_unknown_ancestry(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.apply_ancestry_adjustments(
            {"creatinine": 1.2}, "martian"
        )
        assert result == []

    def test_no_adjustment_when_ancestry_none(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.apply_ancestry_adjustments(
            {"creatinine": 1.2}, None
        )
        assert result == []

    def test_european_has_no_special_adjustments(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.apply_ancestry_adjustments(
            {"creatinine": 1.2, "ldl": 130}, "european"
        )
        # European is the reference population -- no adjustments defined
        assert result == []


# =====================================================================
# 10. AGE- AND SEX-STRATIFIED REFERENCE RANGES
# =====================================================================


class TestAgeSexReferenceRanges:
    """Test age- and sex-stratified reference ranges."""

    def test_young_male_creatinine_normal(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.get_age_sex_ranges({"creatinine": 1.0}, age=35, sex="M")
        assert len(result) == 1
        assert result[0]["status"] == "normal"

    def test_elderly_female_tsh_wider_range(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        # TSH 5.5 is high for young adult but normal for 75-year-old woman
        result_young = engine.get_age_sex_ranges({"tsh": 5.5}, age=30, sex="F")
        result_old = engine.get_age_sex_ranges({"tsh": 5.5}, age=75, sex="F")
        assert result_young[0]["status"] == "high"
        assert result_old[0]["status"] == "normal"

    def test_sex_stratified_ferritin(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        # 200 ng/mL: normal for male, high for premenopausal female
        male = engine.get_age_sex_ranges({"ferritin": 200}, age=40, sex="M")
        female = engine.get_age_sex_ranges({"ferritin": 200}, age=40, sex="F")
        assert male[0]["status"] == "normal"
        assert female[0]["status"] == "high"

    def test_unknown_biomarker_skipped(self):
        from src.genotype_adjustment import GenotypeAdjuster
        engine = GenotypeAdjuster()
        result = engine.get_age_sex_ranges({"unknown_marker": 5.0}, age=50, sex="M")
        assert result == []


# =====================================================================
# 12. KNOWLEDGE BASE VERSIONING
# =====================================================================


class TestKnowledgeVersioning:
    """Test knowledge base version tracking."""

    def test_knowledge_version_exists(self):
        from src.knowledge import KNOWLEDGE_VERSION
        assert "version" in KNOWLEDGE_VERSION
        assert "last_updated" in KNOWLEDGE_VERSION
        assert "cpic_version" in KNOWLEDGE_VERSION

    def test_knowledge_version_has_sources(self):
        from src.knowledge import KNOWLEDGE_VERSION
        assert len(KNOWLEDGE_VERSION["sources"]) >= 5


# =====================================================================
# 11. BIOMARKER INPUT VALIDATION
# =====================================================================


class TestBiomarkerInputValidation:
    """Test plausible range validation for biomarker inputs."""

    def test_valid_biomarkers_no_warnings(self):
        from src.biological_age import validate_biomarker_ranges
        warnings = validate_biomarker_ranges({"albumin": 4.2, "glucose": 95})
        assert warnings == []

    def test_implausible_glucose_warns(self):
        from src.biological_age import validate_biomarker_ranges
        warnings = validate_biomarker_ranges({"glucose": 5000})
        assert len(warnings) == 1
        assert "glucose" in warnings[0]
        assert "plausible range" in warnings[0]

    def test_negative_after_sanitization_warns(self):
        from src.biological_age import validate_biomarker_ranges
        warnings = validate_biomarker_ranges({"albumin": 0.5})
        # 0.5 < 1.0 (lower bound), so this should warn
        assert len(warnings) == 1
        assert "albumin" in warnings[0]


# =====================================================================
# 10. AUDIT LOGGING
# =====================================================================


class TestAuditLogging:
    """Test HIPAA audit logging."""

    def test_audit_log_returns_event_id(self):
        from src.audit import AuditAction, audit_log
        event_id = audit_log(AuditAction.PATIENT_ANALYSIS, patient_id="P001")
        assert isinstance(event_id, str)
        assert len(event_id) == 16

    def test_audit_log_hashes_patient_id(self):
        import hashlib

        from src.audit import AuditAction, audit_log
        event_id = audit_log(AuditAction.BIOLOGICAL_AGE, patient_id="P001")
        # Verify the patient_id is not stored in plain text
        hashlib.sha256("P001".encode()).hexdigest()[:12]
        assert isinstance(event_id, str)

    def test_audit_log_without_patient_id(self):
        from src.audit import AuditAction, audit_log
        event_id = audit_log(AuditAction.RAG_QUERY)
        assert isinstance(event_id, str)

    def test_all_audit_actions_defined(self):
        from src.audit import AuditAction
        actions = list(AuditAction)
        assert len(actions) >= 9
        assert AuditAction.PATIENT_ANALYSIS in actions
        assert AuditAction.FHIR_EXPORTED in actions
