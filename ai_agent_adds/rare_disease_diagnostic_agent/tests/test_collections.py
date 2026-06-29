"""Tests for Milvus collection schemas in scripts/setup_collections.py.

Covers:
  - 14 collections defined
  - Weight sums approximately 1.0
  - Schema field counts
  - Collection names

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path


# Ensure scripts directory is importable
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from scripts.setup_collections import (
    COLLECTION_SCHEMAS,
    EMBEDDING_DIM,
    get_collection_names,
    get_collection_schema,
)


class TestCollectionSchemas:
    """Tests for collection schema definitions."""

    def test_collection_count(self):
        """Must define exactly 14 collections."""
        assert len(COLLECTION_SCHEMAS) == 14

    def test_expected_collection_names(self):
        expected = {
            "rd_phenotypes", "rd_diseases", "rd_genes", "rd_variants",
            "rd_literature", "rd_trials", "rd_therapies", "rd_case_reports",
            "rd_guidelines", "rd_pathways", "rd_registries", "rd_natural_history",
            "rd_newborn_screening", "genomic_evidence",
        }
        assert set(COLLECTION_SCHEMAS.keys()) == expected

    def test_weight_sum_approximately_one(self):
        """Search weights should sum to approximately 1.0."""
        total = sum(cfg["search_weight"] for cfg in COLLECTION_SCHEMAS.values())
        assert abs(total - 1.0) < 0.05, f"Weight sum is {total}, expected ~1.0"

    def test_all_weights_positive(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            assert cfg["search_weight"] > 0, f"{name} has non-positive weight"

    def test_all_have_embedding_field(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            field_names = [f["name"] for f in cfg["fields"]]
            assert "embedding" in field_names, f"{name} missing embedding field"

    def test_embedding_dim(self):
        assert EMBEDDING_DIM == 384

    def test_all_have_id_field(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            field_names = [f["name"] for f in cfg["fields"]]
            assert "id" in field_names, f"{name} missing id field"

    def test_all_have_source_field(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            field_names = [f["name"] for f in cfg["fields"]]
            assert "source" in field_names, f"{name} missing source field"

    def test_all_have_description(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            assert cfg["description"], f"{name} missing description"

    def test_all_have_estimated_records(self):
        for name, cfg in COLLECTION_SCHEMAS.items():
            assert cfg["estimated_records"] >= 0, f"{name} has invalid estimated_records"

    def test_get_collection_names_returns_list(self):
        names = get_collection_names()
        assert isinstance(names, list)
        assert len(names) == 14

    def test_get_collection_schema_valid(self):
        schema = get_collection_schema("rd_phenotypes")
        assert "description" in schema
        assert "fields" in schema

    def test_get_collection_schema_invalid(self):
        schema = get_collection_schema("nonexistent")
        assert schema == {}

    def test_rd_phenotypes_has_hpo_id(self):
        fields = [f["name"] for f in COLLECTION_SCHEMAS["rd_phenotypes"]["fields"]]
        assert "hpo_id" in fields

    def test_rd_diseases_has_disease_name(self):
        fields = [f["name"] for f in COLLECTION_SCHEMAS["rd_diseases"]["fields"]]
        assert "disease_name" in fields
