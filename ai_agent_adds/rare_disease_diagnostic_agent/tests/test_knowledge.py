"""Tests for seed knowledge data in the ingest parsers.

Covers:
  - OMIM disease counts and gene data
  - HPO term counts
  - Orphanet disease counts
  - Gene therapy data (approved + investigational)

Author: Adam Jones
Date: March 2026
"""


from src.ingest.omim_parser import (
    OMIM_DISEASES,
    get_omim_disease_count,
    get_omim_genes,
)
from src.ingest.hpo_parser import (
    HPO_TERMS,
    get_hpo_term_count,
    get_hpo_ids,
)
from src.ingest.orphanet_parser import (
    ORPHANET_DISEASES,
    get_orphanet_disease_count,
    get_orphanet_genes,
)
from src.ingest.gene_therapy_parser import (
    APPROVED_GENE_THERAPIES,
    INVESTIGATIONAL_GENE_THERAPIES,
    ALL_GENE_THERAPIES,
    get_approved_therapy_count,
    get_investigational_therapy_count,
    get_all_therapy_count,
)


# ===================================================================
# OMIM KNOWLEDGE TESTS
# ===================================================================


class TestOMIMKnowledge:
    """Tests for OMIM seed data."""

    def test_omim_disease_count_at_least_30(self):
        assert len(OMIM_DISEASES) >= 30

    def test_get_omim_disease_count(self):
        assert get_omim_disease_count() == len(OMIM_DISEASES)

    def test_all_diseases_have_required_fields(self):
        required = {"omim_id", "disease_name", "gene", "inheritance", "clinical_features"}
        for disease in OMIM_DISEASES:
            for field in required:
                assert field in disease, f"Missing {field} in {disease.get('disease_name', '?')}"

    def test_all_omim_ids_unique(self):
        ids = [d["omim_id"] for d in OMIM_DISEASES]
        assert len(ids) == len(set(ids))

    def test_all_diseases_have_clinical_features(self):
        for disease in OMIM_DISEASES:
            assert len(disease["clinical_features"]) >= 3, (
                f"{disease['disease_name']} has fewer than 3 clinical features"
            )

    def test_get_omim_genes_returns_sorted_unique(self):
        genes = get_omim_genes()
        assert len(genes) > 0
        assert genes == sorted(genes)
        assert len(genes) == len(set(genes))

    def test_known_disease_present(self):
        names = [d["disease_name"] for d in OMIM_DISEASES]
        assert "Cystic Fibrosis" in names
        assert "Duchenne Muscular Dystrophy" in names
        assert "Huntington Disease" in names

    def test_inheritance_patterns(self):
        patterns = {d["inheritance"] for d in OMIM_DISEASES}
        assert "autosomal_recessive" in patterns
        assert "autosomal_dominant" in patterns
        assert "x_linked_recessive" in patterns


# ===================================================================
# HPO KNOWLEDGE TESTS
# ===================================================================


class TestHPOKnowledge:
    """Tests for HPO seed data."""

    def test_hpo_term_count_at_least_50(self):
        assert len(HPO_TERMS) >= 50

    def test_get_hpo_term_count(self):
        assert get_hpo_term_count() == len(HPO_TERMS)

    def test_all_terms_have_required_fields(self):
        required = {"hpo_id", "name", "definition", "synonyms", "ic_score"}
        for term in HPO_TERMS:
            for field in required:
                assert field in term, f"Missing {field} in {term.get('name', '?')}"

    def test_all_hpo_ids_start_with_hp(self):
        for term in HPO_TERMS:
            assert term["hpo_id"].startswith("HP:"), f"Bad HPO ID: {term['hpo_id']}"

    def test_all_ic_scores_positive(self):
        for term in HPO_TERMS:
            assert term["ic_score"] > 0, f"{term['name']} has non-positive IC score"

    def test_get_hpo_ids(self):
        ids = get_hpo_ids()
        assert len(ids) == len(HPO_TERMS)
        assert all(i.startswith("HP:") for i in ids)

    def test_known_term_present(self):
        names = [t["name"] for t in HPO_TERMS]
        assert "Seizures" in names
        assert "Hypotonia" in names
        assert "Intellectual disability" in names

    def test_all_terms_have_synonyms(self):
        for term in HPO_TERMS:
            assert isinstance(term["synonyms"], list)
            assert len(term["synonyms"]) >= 1, f"{term['name']} has no synonyms"


# ===================================================================
# ORPHANET KNOWLEDGE TESTS
# ===================================================================


class TestOrphanetKnowledge:
    """Tests for Orphanet seed data."""

    def test_orphanet_disease_count_at_least_30(self):
        assert len(ORPHANET_DISEASES) >= 30

    def test_get_orphanet_disease_count(self):
        assert get_orphanet_disease_count() == len(ORPHANET_DISEASES)

    def test_all_diseases_have_required_fields(self):
        required = {"orpha_code", "name", "prevalence", "inheritance", "age_onset", "gene"}
        for disease in ORPHANET_DISEASES:
            for field in required:
                assert field in disease, f"Missing {field} in {disease.get('name', '?')}"

    def test_all_orpha_codes_start_with_orpha(self):
        for d in ORPHANET_DISEASES:
            assert d["orpha_code"].startswith("ORPHA:"), f"Bad code: {d['orpha_code']}"

    def test_get_orphanet_genes(self):
        genes = get_orphanet_genes()
        assert len(genes) > 0
        assert "Unknown" not in genes  # Filtered out


# ===================================================================
# GENE THERAPY KNOWLEDGE TESTS
# ===================================================================


class TestGeneTherapyKnowledge:
    """Tests for gene therapy seed data."""

    def test_approved_count_is_6(self):
        assert get_approved_therapy_count() == 6

    def test_investigational_count_is_19(self):
        assert get_investigational_therapy_count() == 19

    def test_total_count(self):
        assert get_all_therapy_count() == 25
        assert get_all_therapy_count() == get_approved_therapy_count() + get_investigational_therapy_count()

    def test_all_therapies_combined(self):
        assert len(ALL_GENE_THERAPIES) == len(APPROVED_GENE_THERAPIES) + len(INVESTIGATIONAL_GENE_THERAPIES)

    def test_all_therapies_have_required_fields(self):
        required = {"drug_name", "target_gene", "indication", "approval_status", "mechanism", "vector_type"}
        for therapy in ALL_GENE_THERAPIES:
            for field in required:
                assert field in therapy, f"Missing {field} in {therapy.get('drug_name', '?')}"

    def test_approved_therapies_have_route(self):
        for therapy in APPROVED_GENE_THERAPIES:
            assert "route" in therapy and therapy["route"], (
                f"{therapy['drug_name']} missing route"
            )

    def test_known_therapies_present(self):
        names = [t["drug_name"] for t in ALL_GENE_THERAPIES]
        drug_names_str = " ".join(names)
        assert "Zolgensma" in drug_names_str
        assert "Luxturna" in drug_names_str
        assert "Casgevy" in drug_names_str

    def test_vector_types_variety(self):
        vectors = {t["vector_type"] for t in ALL_GENE_THERAPIES}
        assert "AAV9" in vectors
        assert "Lentiviral" in vectors
