"""Tests for src/collections.py — PGxCollectionManager and schema definitions.

Tests 15 collection schemas, EMBEDDING_DIM, COLLECTION_SCHEMAS registry,
and PGxCollectionManager class structure.

Author: Adam Jones
Date: March 2026
"""

import pytest
from unittest.mock import MagicMock, patch
from src.collections import (
    EMBEDDING_DIM,
    COLLECTION_SCHEMAS,
    GENE_REFERENCE_SCHEMA,
    DRUG_GUIDELINES_SCHEMA,
    DRUG_INTERACTIONS_SCHEMA,
    HLA_HYPERSENSITIVITY_SCHEMA,
    PHENOCONVERSION_SCHEMA,
    DOSING_ALGORITHMS_SCHEMA,
    CLINICAL_EVIDENCE_SCHEMA,
    POPULATION_DATA_SCHEMA,
    CLINICAL_TRIALS_SCHEMA,
    FDA_LABELS_SCHEMA,
    DRUG_ALTERNATIVES_SCHEMA,
    PATIENT_PROFILES_SCHEMA,
    IMPLEMENTATION_SCHEMA,
    EDUCATION_SCHEMA,
    GENOMIC_EVIDENCE_SCHEMA,
    PGxCollectionManager,
)


# ═══════════════════════════════════════════════════════════════════════
# EMBEDDING_DIM
# ═══════════════════════════════════════════════════════════════════════

class TestEmbeddingDim:
    def test_dimension_is_384(self):
        assert EMBEDDING_DIM == 384


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION_SCHEMAS REGISTRY
# ═══════════════════════════════════════════════════════════════════════

class TestCollectionSchemas:
    def test_has_15_collections(self):
        assert len(COLLECTION_SCHEMAS) == 15

    @pytest.mark.parametrize("name", [
        "pgx_gene_reference", "pgx_drug_guidelines", "pgx_drug_interactions",
        "pgx_hla_hypersensitivity", "pgx_phenoconversion", "pgx_dosing_algorithms",
        "pgx_clinical_evidence", "pgx_population_data", "pgx_clinical_trials",
        "pgx_fda_labels", "pgx_drug_alternatives", "pgx_patient_profiles",
        "pgx_implementation", "pgx_education", "genomic_evidence",
    ])
    def test_collection_present(self, name):
        assert name in COLLECTION_SCHEMAS

    def test_all_schemas_have_id_field(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            field_names = [f.name for f in schema.fields]
            assert "id" in field_names, f"{name} missing 'id' field"

    def test_all_schemas_have_embedding_field(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            field_names = [f.name for f in schema.fields]
            assert "embedding" in field_names, f"{name} missing 'embedding' field"

    def test_embedding_dim_384_in_all(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            emb_field = next(f for f in schema.fields if f.name == "embedding")
            assert emb_field.dim == 384, f"{name} embedding dim != 384"

    def test_id_is_primary_key(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            id_field = next(f for f in schema.fields if f.name == "id")
            assert id_field.is_primary is True, f"{name} 'id' not primary"


# ═══════════════════════════════════════════════════════════════════════
# INDIVIDUAL SCHEMA TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestGeneReferenceSchema:
    def test_has_gene_field(self):
        names = [f.name for f in GENE_REFERENCE_SCHEMA.fields]
        assert "gene" in names

    def test_has_star_allele_field(self):
        names = [f.name for f in GENE_REFERENCE_SCHEMA.fields]
        assert "star_allele" in names

    def test_has_activity_score(self):
        names = [f.name for f in GENE_REFERENCE_SCHEMA.fields]
        assert "activity_score" in names

    def test_has_population_frequencies(self):
        names = [f.name for f in GENE_REFERENCE_SCHEMA.fields]
        assert "allele_frequency_global" in names
        assert "allele_frequency_european" in names
        assert "allele_frequency_african" in names
        assert "allele_frequency_east_asian" in names


class TestDrugGuidelinesSchema:
    def test_has_gene_and_drug(self):
        names = [f.name for f in DRUG_GUIDELINES_SCHEMA.fields]
        assert "gene" in names
        assert "drug" in names

    def test_has_text_chunk(self):
        names = [f.name for f in DRUG_GUIDELINES_SCHEMA.fields]
        assert "text_chunk" in names


class TestHLASchema:
    def test_has_hla_allele(self):
        names = [f.name for f in HLA_HYPERSENSITIVITY_SCHEMA.fields]
        assert "hla_allele" in names

    def test_has_drug(self):
        names = [f.name for f in HLA_HYPERSENSITIVITY_SCHEMA.fields]
        assert "drug" in names


class TestGenomicEvidenceSchema:
    def test_has_chrom_pos(self):
        names = [f.name for f in GENOMIC_EVIDENCE_SCHEMA.fields]
        assert "chrom" in names
        assert "pos" in names

    def test_has_rsid(self):
        names = [f.name for f in GENOMIC_EVIDENCE_SCHEMA.fields]
        assert "rsid" in names


# ═══════════════════════════════════════════════════════════════════════
# PGxCollectionManager
# ═══════════════════════════════════════════════════════════════════════

class TestPGxCollectionManager:
    def test_init_defaults(self):
        mgr = PGxCollectionManager()
        assert mgr.host == "localhost"
        assert mgr.port == 19530
        assert mgr.embedding_dim == 384

    def test_init_custom(self):
        mgr = PGxCollectionManager(host="myhost", port=12345, embedding_dim=768)
        assert mgr.host == "myhost"
        assert mgr.port == 12345
        assert mgr.embedding_dim == 768

    def test_collections_dict_empty_on_init(self):
        mgr = PGxCollectionManager()
        assert mgr._collections == {}

    def test_index_params(self):
        assert PGxCollectionManager.INDEX_PARAMS["metric_type"] == "COSINE"
        assert PGxCollectionManager.INDEX_PARAMS["index_type"] == "IVF_FLAT"
        assert PGxCollectionManager.INDEX_PARAMS["params"]["nlist"] == 1024

    def test_search_params(self):
        assert PGxCollectionManager.SEARCH_PARAMS["metric_type"] == "COSINE"
        assert PGxCollectionManager.SEARCH_PARAMS["params"]["nprobe"] == 16

    @patch("src.collections.connections")
    def test_connect(self, mock_conn):
        mgr = PGxCollectionManager()
        mgr.connect()
        mock_conn.connect.assert_called_once()

    @patch("src.collections.connections")
    def test_disconnect(self, mock_conn):
        mgr = PGxCollectionManager()
        mgr._collections = {"test": MagicMock()}
        mgr.disconnect()
        mock_conn.disconnect.assert_called_once_with("default")
        assert mgr._collections == {}

    @patch("src.collections.utility")
    @patch("src.collections.Collection")
    def test_create_collection_new(self, mock_coll_cls, mock_utility):
        mock_utility.has_collection.return_value = False
        mgr = PGxCollectionManager()
        mock_coll_instance = MagicMock()
        mock_coll_cls.return_value = mock_coll_instance

        result = mgr.create_collection("pgx_gene_reference", GENE_REFERENCE_SCHEMA)
        assert result == mock_coll_instance
        assert "pgx_gene_reference" in mgr._collections

    @patch("src.collections.utility")
    @patch("src.collections.Collection")
    def test_get_output_fields_excludes_embedding(self, mock_coll, mock_util):
        mgr = PGxCollectionManager()
        fields = mgr._get_output_fields("pgx_gene_reference")
        assert "embedding" not in fields
        assert "id" in fields
        assert "gene" in fields

    def test_get_output_fields_unknown_raises(self):
        mgr = PGxCollectionManager()
        with pytest.raises(ValueError, match="Unknown collection"):
            mgr._get_output_fields("nonexistent_collection")
