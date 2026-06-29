"""Tests for query expansion and entity detection.

Tests the ingest parsers' ability to provide aliases, maps, and
entity detection for rare disease terminology.

Author: Adam Jones
Date: March 2026
"""


from src.ingest.omim_parser import OMIM_DISEASES, get_omim_genes
from src.ingest.hpo_parser import HPO_TERMS
from src.ingest.orphanet_parser import ORPHANET_DISEASES
from src.ingest.gene_therapy_parser import ALL_GENE_THERAPIES


# ===================================================================
# GENE ALIAS TESTS
# ===================================================================


class TestGeneAliases:
    """Test gene aliases and mappings from OMIM/Orphanet data."""

    def test_omim_gene_list(self):
        genes = get_omim_genes()
        assert "CFTR" in genes
        assert "DMD" in genes
        assert "HTT" in genes

    def test_gene_to_disease_mapping(self):
        """Each OMIM gene should map to at least one disease."""
        gene_diseases = {}
        for d in OMIM_DISEASES:
            gene = d["gene"]
            if gene not in gene_diseases:
                gene_diseases[gene] = []
            gene_diseases[gene].append(d["disease_name"])
        assert len(gene_diseases) > 20

    def test_orphanet_gene_overlap(self):
        """Some genes should appear in both OMIM and Orphanet."""
        omim_genes = {d["gene"] for d in OMIM_DISEASES}
        orphanet_genes = {d["gene"] for d in ORPHANET_DISEASES if d["gene"] != "Unknown"}
        overlap = omim_genes & orphanet_genes
        assert len(overlap) >= 10, f"Only {len(overlap)} overlapping genes"

    def test_therapy_gene_targets(self):
        """Gene therapy targets should overlap with disease genes."""
        therapy_genes = set()
        for t in ALL_GENE_THERAPIES:
            gene = t["target_gene"].split(" ")[0].strip("()")
            therapy_genes.add(gene)
        disease_genes = {d["gene"] for d in OMIM_DISEASES}
        overlap = therapy_genes & disease_genes
        assert len(overlap) >= 3


# ===================================================================
# HPO SYNONYM TESTS
# ===================================================================


class TestHPOSynonyms:
    """Test HPO synonyms and phenotype term expansion."""

    def test_seizures_has_synonyms(self):
        seizure_terms = [t for t in HPO_TERMS if t["name"] == "Seizures"]
        assert len(seizure_terms) >= 1
        assert "Epilepsy" in seizure_terms[0]["synonyms"]

    def test_hypotonia_has_synonyms(self):
        terms = [t for t in HPO_TERMS if t["name"] == "Hypotonia"]
        assert len(terms) >= 1
        synonyms = terms[0]["synonyms"]
        assert any("floppy" in s.lower() for s in synonyms)

    def test_all_terms_have_at_least_one_synonym(self):
        for term in HPO_TERMS:
            assert len(term["synonyms"]) >= 1, f"{term['name']} has no synonyms"


# ===================================================================
# ENTITY DETECTION TESTS
# ===================================================================


class TestEntityDetection:
    """Test entity detection capabilities from seed data."""

    def test_detect_omim_disease_by_name(self):
        disease_names = {d["disease_name"].lower() for d in OMIM_DISEASES}
        assert "cystic fibrosis" in disease_names

    def test_detect_hpo_term_by_name(self):
        term_names = {t["name"].lower() for t in HPO_TERMS}
        assert "seizures" in term_names

    def test_detect_orphanet_disease_by_name(self):
        disease_names = {d["name"].lower() for d in ORPHANET_DISEASES}
        assert "cystic fibrosis" in disease_names

    def test_detect_gene_by_symbol(self):
        omim_genes = {d["gene"] for d in OMIM_DISEASES}
        assert "CFTR" in omim_genes
        assert "DMD" in omim_genes

    def test_detect_inheritance_pattern(self):
        patterns = {d["inheritance"] for d in OMIM_DISEASES}
        assert "autosomal_recessive" in patterns
        assert "autosomal_dominant" in patterns
        assert "x_linked_recessive" in patterns

    def test_ic_score_range(self):
        """IC scores should be positive and reasonable."""
        for term in HPO_TERMS:
            assert 0 < term["ic_score"] < 20, f"{term['name']} IC={term['ic_score']}"
