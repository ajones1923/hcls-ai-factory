"""Tests for collection configurations in scripts/setup_collections.py.

Covers:
  - Collection count (14 expected)
  - Required fields present in each schema
  - Weight sums approximately equal to 1.0
  - Schema completeness

Author: Adam Jones
Date: March 2026
"""

import pytest

from scripts.setup_collections import (
    COLLECTION_SCHEMAS,
    EMBEDDING_DIM,
    INDEX_TYPE,
    METRIC_TYPE,
    get_collection_names,
    get_collection_schema,
)


class TestCollectionCount:
    """Test that all 14 collections are defined."""

    def test_total_count(self):
        assert len(COLLECTION_SCHEMAS) == 14

    def test_get_collection_names(self):
        names = get_collection_names()
        assert len(names) == 14
        assert isinstance(names, list)


class TestCollectionNames:
    """Test expected collection names are present."""

    EXPECTED_NAMES = [
        "trial_protocols",
        "trial_eligibility",
        "trial_endpoints",
        "trial_sites",
        "trial_investigators",
        "trial_results",
        "trial_regulatory",
        "trial_literature",
        "trial_biomarkers",
        "trial_safety",
        "trial_rwe",
        "trial_adaptive",
        "trial_guidelines",
        "genomic_evidence",
    ]

    @pytest.mark.parametrize("name", EXPECTED_NAMES)
    def test_collection_exists(self, name):
        assert name in COLLECTION_SCHEMAS, f"Collection '{name}' not found"


class TestCollectionSchemas:
    """Test schema completeness for each collection."""

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_description(self, name):
        schema = COLLECTION_SCHEMAS[name]
        assert "description" in schema
        assert len(schema["description"]) > 0

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_fields(self, name):
        schema = COLLECTION_SCHEMAS[name]
        assert "fields" in schema
        assert len(schema["fields"]) >= 3  # At least id, embedding, text

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_id_field(self, name):
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f["name"] for f in schema["fields"]]
        assert "id" in field_names

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_embedding_field(self, name):
        schema = COLLECTION_SCHEMAS[name]
        embedding_fields = [
            f for f in schema["fields"]
            if f["name"] == "embedding"
        ]
        assert len(embedding_fields) == 1
        assert embedding_fields[0]["dtype"] == "FLOAT_VECTOR"
        assert embedding_fields[0]["dim"] == EMBEDDING_DIM

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_text_field(self, name):
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f["name"] for f in schema["fields"]]
        assert "text" in field_names

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_source_field(self, name):
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f["name"] for f in schema["fields"]]
        assert "source" in field_names

    @pytest.mark.parametrize("name", list(COLLECTION_SCHEMAS.keys()))
    def test_has_search_weight(self, name):
        schema = COLLECTION_SCHEMAS[name]
        assert "search_weight" in schema
        assert 0.0 <= schema["search_weight"] <= 1.0


class TestCollectionWeights:
    """Test that collection weights sum to approximately 1.0."""

    def test_weights_sum(self):
        total = sum(s["search_weight"] for s in COLLECTION_SCHEMAS.values())
        assert abs(total - 1.0) <= 0.05, (
            f"Collection weights sum to {total:.4f}, expected ~1.0"
        )

    def test_all_weights_positive(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            assert schema["search_weight"] > 0, (
                f"Collection '{name}' has zero or negative weight"
            )


class TestConstants:
    """Test embedding and index constants."""

    def test_embedding_dim(self):
        assert EMBEDDING_DIM == 384

    def test_index_type(self):
        assert INDEX_TYPE == "IVF_FLAT"

    def test_metric_type(self):
        assert METRIC_TYPE == "COSINE"


class TestGetCollectionSchema:
    """Test the get_collection_schema helper."""

    def test_existing_collection(self):
        schema = get_collection_schema("trial_protocols")
        assert schema is not None
        assert "fields" in schema

    def test_nonexistent_collection(self):
        schema = get_collection_schema("nonexistent")
        assert schema == {}
