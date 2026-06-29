"""Tests for Imaging Intelligence Agent domain knowledge graph.

Validates IMAGING_PATHOLOGIES, IMAGING_MODALITIES, IMAGING_ANATOMY dicts,
and all helper functions (get_*_context, get_nim_recommendation,
resolve_comparison_entity, get_comparison_context, get_knowledge_stats).

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.knowledge import (
    IMAGING_ANATOMY,
    IMAGING_MODALITIES,
    IMAGING_PATHOLOGIES,
    get_anatomy_context,
    get_comparison_context,
    get_knowledge_stats,
    get_modality_context,
    get_nim_recommendation,
    get_pathology_context,
    resolve_comparison_entity,
)


# ===================================================================
# IMAGING_PATHOLOGIES
# ===================================================================


class TestImagingPathologies:
    """Tests for the IMAGING_PATHOLOGIES knowledge dictionary."""

    def test_intracranial_hemorrhage_exists(self):
        assert "intracranial_hemorrhage" in IMAGING_PATHOLOGIES

    def test_lung_nodule_exists(self):
        assert "lung_nodule" in IMAGING_PATHOLOGIES

    def test_pneumonia_exists(self):
        assert "pneumonia" in IMAGING_PATHOLOGIES

    def test_pulmonary_embolism_exists(self):
        assert "pulmonary_embolism" in IMAGING_PATHOLOGIES

    def test_stroke_ischemic_exists(self):
        assert "stroke_ischemic" in IMAGING_PATHOLOGIES

    def test_brain_tumor_exists(self):
        assert "brain_tumor" in IMAGING_PATHOLOGIES

    def test_ms_lesion_exists(self):
        assert "ms_lesion" in IMAGING_PATHOLOGIES

    def test_pneumothorax_exists(self):
        assert "pneumothorax" in IMAGING_PATHOLOGIES

    def test_pleural_effusion_exists(self):
        assert "pleural_effusion" in IMAGING_PATHOLOGIES

    def test_aortic_dissection_exists(self):
        assert "aortic_dissection" in IMAGING_PATHOLOGIES

    def test_fracture_exists(self):
        assert "fracture" in IMAGING_PATHOLOGIES

    def test_has_at_least_10_pathologies(self):
        assert len(IMAGING_PATHOLOGIES) >= 10

    def test_each_pathology_has_required_keys(self):
        required = {"icd10", "display_name", "modalities", "body_region", "subtypes"}
        for key, data in IMAGING_PATHOLOGIES.items():
            missing = required - set(data.keys())
            assert not missing, f"{key} missing keys: {missing}"

    def test_intracranial_hemorrhage_has_ct_modality(self):
        assert "ct" in IMAGING_PATHOLOGIES["intracranial_hemorrhage"]["modalities"]

    def test_lung_nodule_body_region_is_chest(self):
        assert IMAGING_PATHOLOGIES["lung_nodule"]["body_region"] == "chest"


# ===================================================================
# IMAGING_MODALITIES
# ===================================================================


class TestImagingModalities:
    """Tests for the IMAGING_MODALITIES knowledge dictionary."""

    def test_ct_exists(self):
        assert "ct" in IMAGING_MODALITIES

    def test_mri_exists(self):
        assert "mri" in IMAGING_MODALITIES

    def test_xray_exists(self):
        assert "xray" in IMAGING_MODALITIES

    def test_cxr_exists(self):
        assert "cxr" in IMAGING_MODALITIES

    def test_pet_ct_exists(self):
        assert "pet_ct" in IMAGING_MODALITIES

    def test_ultrasound_exists(self):
        assert "ultrasound" in IMAGING_MODALITIES

    def test_mammography_exists(self):
        assert "mammography" in IMAGING_MODALITIES

    def test_fluoroscopy_exists(self):
        assert "fluoroscopy" in IMAGING_MODALITIES

    def test_has_8_modalities(self):
        assert len(IMAGING_MODALITIES) == 8

    def test_each_modality_has_required_keys(self):
        required = {"full_name", "dicom_modality_code", "physics", "strengths", "limitations"}
        for key, data in IMAGING_MODALITIES.items():
            missing = required - set(data.keys())
            assert not missing, f"{key} missing keys: {missing}"

    def test_ct_full_name(self):
        assert IMAGING_MODALITIES["ct"]["full_name"] == "Computed Tomography"

    def test_mri_is_non_ionizing(self):
        assert IMAGING_MODALITIES["mri"]["typical_dose_msv"] is None


# ===================================================================
# IMAGING_ANATOMY
# ===================================================================


class TestImagingAnatomy:
    """Tests for the IMAGING_ANATOMY knowledge dictionary."""

    def test_brain_exists(self):
        assert "brain" in IMAGING_ANATOMY

    def test_lungs_exists(self):
        assert "lungs" in IMAGING_ANATOMY

    def test_heart_exists(self):
        assert "heart" in IMAGING_ANATOMY

    def test_liver_exists(self):
        assert "liver" in IMAGING_ANATOMY

    def test_spine_exists(self):
        assert "spine" in IMAGING_ANATOMY

    def test_has_at_least_10_anatomy_entries(self):
        assert len(IMAGING_ANATOMY) >= 10

    def test_each_anatomy_has_required_keys(self):
        required = {"display_name", "structures", "vista3d_labels", "preferred_modality"}
        for key, data in IMAGING_ANATOMY.items():
            missing = required - set(data.keys())
            assert not missing, f"{key} missing keys: {missing}"

    def test_brain_preferred_modality_is_mri(self):
        assert IMAGING_ANATOMY["brain"]["preferred_modality"] == "mri"

    def test_lungs_preferred_modality_is_ct(self):
        assert IMAGING_ANATOMY["lungs"]["preferred_modality"] == "ct"


# ===================================================================
# HELPER FUNCTIONS
# ===================================================================


class TestGetPathologyContext:
    """Tests for get_pathology_context()."""

    def test_returns_context_for_known_pathology(self):
        ctx = get_pathology_context("intracranial_hemorrhage")
        assert len(ctx) > 0
        assert "Intracranial Hemorrhage" in ctx

    def test_returns_icd10_in_context(self):
        ctx = get_pathology_context("lung_nodule")
        assert "R91.1" in ctx

    def test_case_insensitive(self):
        ctx = get_pathology_context("INTRACRANIAL_HEMORRHAGE")
        assert len(ctx) > 0

    def test_space_to_underscore(self):
        ctx = get_pathology_context("lung nodule")
        assert "Pulmonary Nodule" in ctx

    def test_returns_empty_for_unknown(self):
        ctx = get_pathology_context("nonexistent_disease_xyz")
        assert ctx == ""

    def test_partial_match(self):
        ctx = get_pathology_context("hemorrhage")
        assert len(ctx) > 0


class TestGetModalityContext:
    """Tests for get_modality_context()."""

    def test_returns_context_for_ct(self):
        ctx = get_modality_context("ct")
        assert "Computed Tomography" in ctx
        assert len(ctx) > 0

    def test_returns_context_for_mri(self):
        ctx = get_modality_context("mri")
        assert "Magnetic Resonance" in ctx

    def test_returns_empty_for_unknown(self):
        ctx = get_modality_context("xyz_nonexistent")
        assert ctx == ""

    def test_includes_physics(self):
        ctx = get_modality_context("ct")
        assert "Physics:" in ctx

    def test_includes_strengths(self):
        ctx = get_modality_context("mri")
        assert "Strengths:" in ctx


class TestGetAnatomyContext:
    """Tests for get_anatomy_context()."""

    def test_returns_context_for_brain(self):
        ctx = get_anatomy_context("brain")
        assert "Brain" in ctx
        assert len(ctx) > 0

    def test_returns_context_for_lungs(self):
        ctx = get_anatomy_context("lungs")
        assert "Lungs" in ctx

    def test_returns_empty_for_unknown(self):
        ctx = get_anatomy_context("xyz_nonexistent_anatomy")
        assert ctx == ""

    def test_includes_vista3d_labels(self):
        ctx = get_anatomy_context("brain")
        assert "VISTA-3D Labels:" in ctx

    def test_includes_preferred_modality(self):
        ctx = get_anatomy_context("heart")
        assert "Preferred Modality:" in ctx


class TestGetNimRecommendation:
    """Tests for get_nim_recommendation()."""

    def test_returns_workflow_for_hemorrhage(self):
        result = get_nim_recommendation("intracranial_hemorrhage")
        assert result == "ct_head_hemorrhage"

    def test_returns_workflow_for_lung_nodule(self):
        result = get_nim_recommendation("lung_nodule")
        assert result == "ct_chest_lung_nodule"

    def test_returns_workflow_for_pneumonia(self):
        result = get_nim_recommendation("pneumonia")
        assert result == "cxr_rapid_findings"

    def test_returns_none_for_no_workflow(self):
        result = get_nim_recommendation("pulmonary_embolism")
        assert result is None

    def test_returns_none_for_unknown(self):
        result = get_nim_recommendation("nonexistent_xyz")
        assert result is None


class TestResolveComparisonEntity:
    """Tests for resolve_comparison_entity()."""

    def test_resolves_ct_as_modality(self):
        result = resolve_comparison_entity("ct")
        assert result is not None
        assert result["type"] == "modality"
        assert result["canonical"] == "ct"

    def test_resolves_brain_as_anatomy(self):
        result = resolve_comparison_entity("brain")
        assert result is not None
        assert result["type"] == "anatomy"

    def test_resolves_hemorrhage_as_pathology(self):
        result = resolve_comparison_entity("intracranial_hemorrhage")
        assert result is not None
        assert result["type"] == "pathology"

    def test_returns_none_for_unknown(self):
        result = resolve_comparison_entity("xyzzy_nonexistent")
        assert result is None


class TestGetComparisonContext:
    """Tests for get_comparison_context()."""

    def test_returns_context_for_two_modalities(self):
        entity_a = {"type": "modality", "canonical": "ct"}
        entity_b = {"type": "modality", "canonical": "mri"}
        ctx = get_comparison_context(entity_a, entity_b)
        assert "Computed Tomography" in ctx
        assert "Magnetic Resonance" in ctx

    def test_returns_context_for_pathology_and_modality(self):
        entity_a = {"type": "pathology", "canonical": "lung_nodule"}
        entity_b = {"type": "modality", "canonical": "ct"}
        ctx = get_comparison_context(entity_a, entity_b)
        assert len(ctx) > 0


class TestGetKnowledgeStats:
    """Tests for get_knowledge_stats()."""

    def test_returns_dict_with_expected_keys(self):
        stats = get_knowledge_stats()
        assert "pathologies" in stats
        assert "modalities" in stats
        assert "anatomy_regions" in stats

    def test_pathology_count_at_least_10(self):
        stats = get_knowledge_stats()
        assert stats["pathologies"] >= 10

    def test_modality_count_is_8(self):
        stats = get_knowledge_stats()
        assert stats["modalities"] == 8

    def test_anatomy_count_at_least_10(self):
        stats = get_knowledge_stats()
        assert stats["anatomy_regions"] >= 10

    def test_pathologies_with_nim_workflow(self):
        stats = get_knowledge_stats()
        assert stats["pathologies_with_nim_workflow"] >= 3

    def test_modalities_ionizing_count(self):
        stats = get_knowledge_stats()
        assert stats["modalities_ionizing"] >= 4  # ct, xray, cxr, pet_ct, mammography, fluoroscopy
