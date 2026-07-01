"""
Tests for the OncoCollectionManager and collection schemas.
============================================================
Validates that COLLECTION_SCHEMAS, COLLECTION_MODELS, and field
definitions are correctly structured for all 11 Milvus collections.

All Milvus client calls are mocked -- no live connection required.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.collections import (
    COLLECTION_MODELS,
    COLLECTION_SCHEMAS,
    EMBEDDING_DIM,
    GENOMIC_EVIDENCE_FIELDS,
    ONCO_BIOMARKERS_FIELDS,
    ONCO_CASES_FIELDS,
    ONCO_GUIDELINES_FIELDS,
    ONCO_LITERATURE_FIELDS,
    ONCO_OUTCOMES_FIELDS,
    ONCO_PATHWAYS_FIELDS,
    ONCO_RESISTANCE_FIELDS,
    ONCO_THERAPIES_FIELDS,
    ONCO_TRIALS_FIELDS,
    ONCO_VARIANTS_FIELDS,
    OncoCollectionManager,
)
from src.models import (
    CaseSnapshot,
    OncologyBiomarker,
    OncologyGuideline,
    OncologyLiterature,
    OncologyPathway,
    OncologyTherapy,
    OncologyTrial,
    OncologyVariant,
    OutcomeRecord,
    ResistanceMechanism,
)


# ═══════════════════════════════════════════════════════════════════════════
# Schema Registry Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCollectionSchemas:
    """Verify COLLECTION_SCHEMAS has all 11 collections defined."""

    EXPECTED_COLLECTIONS = [
        "onco_literature",
        "onco_trials",
        "onco_variants",
        "onco_biomarkers",
        "onco_therapies",
        "onco_pathways",
        "onco_guidelines",
        "onco_resistance",
        "onco_outcomes",
        "onco_cases",
        "genomic_evidence",
    ]

    def test_schema_count(self):
        """COLLECTION_SCHEMAS should have exactly 11 entries."""
        assert len(COLLECTION_SCHEMAS) == 11

    def test_all_collections_present(self):
        """Every expected collection should be in COLLECTION_SCHEMAS."""
        for name in self.EXPECTED_COLLECTIONS:
            assert name in COLLECTION_SCHEMAS, f"Missing collection: {name}"

    def test_no_extra_collections(self):
        """COLLECTION_SCHEMAS should not contain unexpected collections."""
        for name in COLLECTION_SCHEMAS:
            assert name in self.EXPECTED_COLLECTIONS, f"Unexpected collection: {name}"


class TestCollectionModels:
    """Verify COLLECTION_MODELS maps collections to correct Pydantic models."""

    EXPECTED_MAPPINGS = {
        "onco_literature": OncologyLiterature,
        "onco_trials": OncologyTrial,
        "onco_variants": OncologyVariant,
        "onco_biomarkers": OncologyBiomarker,
        "onco_therapies": OncologyTherapy,
        "onco_pathways": OncologyPathway,
        "onco_guidelines": OncologyGuideline,
        "onco_resistance": ResistanceMechanism,
        "onco_outcomes": OutcomeRecord,
        "onco_cases": CaseSnapshot,
        "genomic_evidence": None,
    }

    def test_model_count(self):
        """COLLECTION_MODELS should have 11 entries."""
        assert len(COLLECTION_MODELS) == 11

    def test_all_mappings_correct(self):
        """Each collection should map to the correct Pydantic model."""
        for name, expected_model in self.EXPECTED_MAPPINGS.items():
            assert COLLECTION_MODELS[name] is expected_model, (
                f"{name} should map to {expected_model}"
            )

    def test_genomic_evidence_is_none(self):
        """genomic_evidence should map to None (read-only, populated by Stage 1)."""
        assert COLLECTION_MODELS["genomic_evidence"] is None


# ═══════════════════════════════════════════════════════════════════════════
# Field Schema Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestFieldSchemas:
    """Verify required fields are present in every collection schema."""

    ALL_FIELD_LISTS = {
        "onco_literature": ONCO_LITERATURE_FIELDS,
        "onco_trials": ONCO_TRIALS_FIELDS,
        "onco_variants": ONCO_VARIANTS_FIELDS,
        "onco_biomarkers": ONCO_BIOMARKERS_FIELDS,
        "onco_therapies": ONCO_THERAPIES_FIELDS,
        "onco_pathways": ONCO_PATHWAYS_FIELDS,
        "onco_guidelines": ONCO_GUIDELINES_FIELDS,
        "onco_resistance": ONCO_RESISTANCE_FIELDS,
        "onco_outcomes": ONCO_OUTCOMES_FIELDS,
        "onco_cases": ONCO_CASES_FIELDS,
        "genomic_evidence": GENOMIC_EVIDENCE_FIELDS,
    }

    @pytest.mark.parametrize("collection_name", ALL_FIELD_LISTS.keys())
    def test_has_id_field(self, collection_name):
        """Every collection must have an 'id' primary key field."""
        fields = self.ALL_FIELD_LISTS[collection_name]
        field_names = [f.name for f in fields]
        assert "id" in field_names, f"{collection_name} missing 'id' field"

    @pytest.mark.parametrize("collection_name", ALL_FIELD_LISTS.keys())
    def test_has_embedding_field(self, collection_name):
        """Every collection must have an 'embedding' vector field."""
        fields = self.ALL_FIELD_LISTS[collection_name]
        field_names = [f.name for f in fields]
        assert "embedding" in field_names, f"{collection_name} missing 'embedding' field"

    @pytest.mark.parametrize("collection_name", [
        "onco_literature", "onco_variants", "onco_biomarkers",
        "onco_therapies", "onco_pathways", "onco_guidelines",
        "onco_resistance", "onco_outcomes", "onco_cases",
        "genomic_evidence",
    ])
    def test_has_text_or_summary_field(self, collection_name):
        """Most collections should have a text_chunk or text_summary field."""
        fields = self.ALL_FIELD_LISTS[collection_name]
        field_names = [f.name for f in fields]
        has_text = any(
            fn in field_names for fn in ("text_chunk", "text_summary", "text")
        )
        assert has_text, f"{collection_name} missing text content field"


class TestEmbeddingDimension:
    """Verify embedding dimension is 384 for all collections."""

    def test_constant_is_384(self):
        """EMBEDDING_DIM constant should be 384."""
        assert EMBEDDING_DIM == 384

    @pytest.mark.parametrize("collection_name", COLLECTION_SCHEMAS.keys())
    def test_embedding_dim_in_schema(self, collection_name):
        """Each collection's embedding field should have dim=384."""
        schema = COLLECTION_SCHEMAS[collection_name]
        embedding_fields = [f for f in schema.fields if f.name == "embedding"]
        assert len(embedding_fields) == 1, f"{collection_name} should have exactly one embedding field"
        assert embedding_fields[0].dim == 384, (
            f"{collection_name} embedding dim should be 384, got {embedding_fields[0].dim}"
        )


# ═══════════════════════════════════════════════════════════════════════════
# Settings Consistency Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCollectionNamesMatchSettings:
    """Verify collection names match the defaults in OncoSettings."""

    def test_settings_collection_names(self):
        """OncoSettings default collection names should match COLLECTION_SCHEMAS keys."""
        try:
            from config.settings import OncoSettings
            s = OncoSettings()
            expected_mapping = {
                s.COLLECTION_LITERATURE: "onco_literature",
                s.COLLECTION_TRIALS: "onco_trials",
                s.COLLECTION_VARIANTS: "onco_variants",
                s.COLLECTION_BIOMARKERS: "onco_biomarkers",
                s.COLLECTION_THERAPIES: "onco_therapies",
                s.COLLECTION_PATHWAYS: "onco_pathways",
                s.COLLECTION_GUIDELINES: "onco_guidelines",
                s.COLLECTION_RESISTANCE: "onco_resistance",
                s.COLLECTION_OUTCOMES: "onco_outcomes",
                s.COLLECTION_CASES: "onco_cases",
                s.COLLECTION_GENOMIC: "genomic_evidence",
            }
            for setting_val, expected in expected_mapping.items():
                assert setting_val == expected, (
                    f"Settings value '{setting_val}' does not match expected '{expected}'"
                )
        except ImportError:
            pytest.skip("config.settings not importable in test environment")


# ═══════════════════════════════════════════════════════════════════════════
# OncoCollectionManager Unit Tests (Milvus mocked)
# ═══════════════════════════════════════════════════════════════════════════


class TestOncoCollectionManagerInit:
    """Test OncoCollectionManager initialization."""

    def test_default_init(self):
        manager = OncoCollectionManager()
        assert manager.host == "localhost"
        assert manager.port == 19530
        assert manager.alias == "default"

    def test_custom_init(self):
        manager = OncoCollectionManager(host="milvus-server", port=19531, alias="test")
        assert manager.host == "milvus-server"
        assert manager.port == 19531
        assert manager.alias == "test"

    def test_internal_collections_empty(self):
        manager = OncoCollectionManager()
        assert manager._collections == {}


class TestOncoCollectionManagerCreateCollection:
    """Test create_collection with mocked Milvus."""

    @patch("src.collections.utility")
    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_create_unknown_collection_raises(self, mock_conn, mock_col, mock_util):
        """Creating an unknown collection should raise ValueError."""
        manager = OncoCollectionManager()
        with pytest.raises(ValueError, match="Unknown collection"):
            manager.create_collection("nonexistent_collection")

    @patch("src.collections.utility")
    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_create_new_collection(self, mock_conn, mock_col_cls, mock_util):
        """Creating a valid new collection should succeed."""
        mock_util.has_collection.return_value = False
        mock_col_instance = MagicMock()
        mock_col_cls.return_value = mock_col_instance

        manager = OncoCollectionManager()
        col = manager.create_collection("onco_literature")

        assert col is mock_col_instance
        mock_col_instance.create_index.assert_called_once()
        assert "onco_literature" in manager._collections
