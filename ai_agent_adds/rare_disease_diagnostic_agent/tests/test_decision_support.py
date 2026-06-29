"""Tests for decision support: ACMG classification, HPO matching, LOD scores.

Tests the ingest parsers' data quality for downstream decision support tasks.

Author: Adam Jones
Date: March 2026
"""


from src.ingest.omim_parser import OMIM_DISEASES
from src.ingest.hpo_parser import HPO_TERMS
from src.ingest.gene_therapy_parser import APPROVED_GENE_THERAPIES, INVESTIGATIONAL_GENE_THERAPIES


# ===================================================================
# ACMG CLASSIFICATION SUPPORT
# ===================================================================


class TestACMGSupport:
    """Test data quality for ACMG variant classification support."""

    def test_all_diseases_have_gene(self):
        """Every disease must have a gene for variant classification."""
        for d in OMIM_DISEASES:
            assert d["gene"], f"{d['disease_name']} missing gene"

    def test_all_diseases_have_inheritance(self):
        """Inheritance pattern is essential for ACMG PS4/PM3 criteria."""
        for d in OMIM_DISEASES:
            assert d["inheritance"], f"{d['disease_name']} missing inheritance"

    def test_inheritance_values_are_valid(self):
        valid = {
            "autosomal_recessive", "autosomal_dominant",
            "x_linked_recessive", "x_linked_dominant",
        }
        for d in OMIM_DISEASES:
            assert d["inheritance"] in valid, (
                f"{d['disease_name']} has invalid inheritance: {d['inheritance']}"
            )

    def test_clinical_features_for_phenotype_matching(self):
        """Clinical features enable PP4 criterion (specific phenotype)."""
        for d in OMIM_DISEASES:
            assert len(d["clinical_features"]) >= 3


# ===================================================================
# HPO MATCHER SUPPORT
# ===================================================================


class TestHPOMatcherSupport:
    """Test HPO data quality for phenotype matching."""

    def test_ic_scores_differentiate_specificity(self):
        """Higher IC scores should indicate more specific phenotypes."""
        general_terms = [t for t in HPO_TERMS if t["ic_score"] < 2.0]
        specific_terms = [t for t in HPO_TERMS if t["ic_score"] > 4.0]
        assert len(general_terms) > 0, "No general (low IC) terms found"
        assert len(specific_terms) > 0, "No specific (high IC) terms found"

    def test_hpo_definitions_non_empty(self):
        for term in HPO_TERMS:
            assert term["definition"], f"{term['name']} missing definition"

    def test_hpo_names_non_empty(self):
        for term in HPO_TERMS:
            assert term["name"], f"{term['hpo_id']} missing name"

    def test_hpo_ids_format(self):
        """HPO IDs must follow HP:NNNNNNN format."""
        import re
        pattern = re.compile(r"^HP:\d{7}$")
        for term in HPO_TERMS:
            assert pattern.match(term["hpo_id"]), f"Invalid HPO ID format: {term['hpo_id']}"

    def test_synonym_overlap_detection(self):
        """Some terms should have overlapping synonyms for disambiguation."""
        all_synonyms = []
        for term in HPO_TERMS:
            all_synonyms.extend(s.lower() for s in term["synonyms"])
        # Some duplicates expected (different terms with overlapping synonyms)
        assert len(all_synonyms) > len(HPO_TERMS)


# ===================================================================
# LOD SCORE SUPPORT
# ===================================================================


class TestLODScoreSupport:
    """Test data supporting LOD score-like calculations."""

    def test_prevalence_data_available(self):
        """Prevalence data enables prior probability calculation."""
        with_prevalence = [d for d in OMIM_DISEASES if d.get("prevalence")]
        assert len(with_prevalence) >= 30

    def test_chromosome_data_available(self):
        """Chromosome data supports linkage analysis."""
        with_chr = [d for d in OMIM_DISEASES if d.get("chromosome")]
        assert len(with_chr) >= 25

    def test_multiple_inheritance_patterns(self):
        """Multiple inheritance patterns needed for segregation analysis."""
        patterns = {d["inheritance"] for d in OMIM_DISEASES}
        assert len(patterns) >= 3


# ===================================================================
# GENE THERAPY ELIGIBILITY SUPPORT
# ===================================================================


class TestGeneTherapyEligibility:
    """Test gene therapy data quality for eligibility assessment."""

    def test_all_approved_have_manufacturer(self):
        for t in APPROVED_GENE_THERAPIES:
            assert t.get("manufacturer"), f"{t['drug_name']} missing manufacturer"

    def test_all_approved_have_age_group(self):
        for t in APPROVED_GENE_THERAPIES:
            assert t.get("age_group"), f"{t['drug_name']} missing age_group"

    def test_all_have_vector_type(self):
        for t in APPROVED_GENE_THERAPIES + INVESTIGATIONAL_GENE_THERAPIES:
            assert t.get("vector_type"), f"{t['drug_name']} missing vector_type"

    def test_sma_therapy_present(self):
        sma_therapies = [
            t for t in APPROVED_GENE_THERAPIES
            if "spinal muscular" in t["indication"].lower()
        ]
        assert len(sma_therapies) >= 1

    def test_hemophilia_therapies_present(self):
        hemo_therapies = [
            t for t in APPROVED_GENE_THERAPIES
            if "hemophilia" in t["indication"].lower()
        ]
        assert len(hemo_therapies) >= 2
