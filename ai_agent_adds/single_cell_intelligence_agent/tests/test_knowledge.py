"""Tests for domain knowledge base in src/knowledge.py.

Covers cell type atlas, TME profiles, spatial platforms, immune signatures,
marker genes, foundation models, and knowledge version metadata.

Author: Adam Jones
Date: March 2026
"""


from src.knowledge import (
    CELL_TYPE_ATLAS,
    FOUNDATION_MODELS,
    GPU_BENCHMARKS,
    IMMUNE_SIGNATURES,
    KNOWLEDGE_VERSION,
    SPATIAL_PLATFORMS,
    TME_PROFILES,
)


class TestKnowledgeVersion:
    """Tests for KNOWLEDGE_VERSION metadata."""

    def test_version_exists(self):
        assert KNOWLEDGE_VERSION["version"]

    def test_has_counts(self):
        counts = KNOWLEDGE_VERSION["counts"]
        assert counts["cell_types"] >= 30
        assert counts["tme_profiles"] >= 4

    def test_has_sources(self):
        assert len(KNOWLEDGE_VERSION["sources"]) > 0


class TestCellTypeAtlas:
    """Tests for CELL_TYPE_ATLAS data."""

    def test_cell_type_count(self):
        assert len(CELL_TYPE_ATLAS) >= 30

    def test_expected_cell_types(self):
        expected = ["T_cell", "CD8_T", "B_cell", "Macrophage", "Fibroblast"]
        for ct in expected:
            assert ct in CELL_TYPE_ATLAS, f"Missing cell type: {ct}"

    def test_cell_types_have_markers(self):
        for name, ct in CELL_TYPE_ATLAS.items():
            assert ct["markers"], f"Cell type {name} missing markers"

    def test_cell_types_have_descriptions(self):
        for name, ct in CELL_TYPE_ATLAS.items():
            assert ct["description"], f"Cell type {name} missing description"

    def test_cell_types_have_tissues(self):
        for name, ct in CELL_TYPE_ATLAS.items():
            assert ct["tissues"], f"Cell type {name} missing tissues"

    def test_cell_ontology_ids(self):
        for name, ct in CELL_TYPE_ATLAS.items():
            assert ct.get("cell_ontology_id", "").startswith("CL:"), (
                f"Cell type {name} missing valid Cell Ontology ID"
            )


class TestTMEProfiles:
    """Tests for TME_PROFILES data."""

    def test_tme_profile_count(self):
        assert len(TME_PROFILES) >= 4

    def test_expected_profiles(self):
        expected = ["hot", "cold", "excluded"]
        for p in expected:
            assert p in TME_PROFILES, f"Missing TME profile: {p}"

    def test_profiles_have_signature_genes(self):
        for name, profile in TME_PROFILES.items():
            assert profile.get("signature_genes"), (
                f"TME profile {name} missing signature genes"
            )

    def test_profiles_have_descriptions(self):
        for name, profile in TME_PROFILES.items():
            assert profile.get("description"), (
                f"TME profile {name} missing description"
            )


class TestSpatialPlatforms:
    """Tests for SPATIAL_PLATFORMS data."""

    def test_platform_count(self):
        assert len(SPATIAL_PLATFORMS) >= 4

    def test_expected_platforms(self):
        expected = ["Visium", "MERFISH"]
        for p in expected:
            assert p in SPATIAL_PLATFORMS, f"Missing spatial platform: {p}"


class TestImmuneSignatures:
    """Tests for IMMUNE_SIGNATURES data."""

    def test_signature_count(self):
        assert len(IMMUNE_SIGNATURES) >= 4

    def test_signatures_have_genes(self):
        for name, sig in IMMUNE_SIGNATURES.items():
            assert sig.get("genes") or sig.get("signature_genes"), (
                f"Immune signature {name} missing genes"
            )


class TestFoundationModels:
    """Tests for FOUNDATION_MODELS data."""

    def test_model_count(self):
        assert len(FOUNDATION_MODELS) >= 3

    def test_models_have_info(self):
        for name, model in FOUNDATION_MODELS.items():
            # Models may use "description" or "full_name" as primary descriptor
            has_info = model.get("description") or model.get("full_name")
            assert has_info, (
                f"Foundation model {name} missing description or full_name"
            )


class TestGPUBenchmarks:
    """Tests for GPU_BENCHMARKS data."""

    def test_benchmark_count(self):
        assert len(GPU_BENCHMARKS) >= 3
