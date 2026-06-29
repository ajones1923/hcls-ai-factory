"""Tests for Imaging Intelligence Agent query expansion module.

Validates all 12 expansion maps and the expand_query() function
for synonym resolution, case insensitivity, and multi-map coverage.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.query_expansion import (
    AI_TASK_EXPANSION,
    ALL_EXPANSION_MAPS,
    BODY_REGION_EXPANSION,
    CONTRAST_EXPANSION,
    DATASET_EXPANSION,
    DEVICE_EXPANSION,
    FINDING_EXPANSION,
    GUIDELINE_EXPANSION,
    MEASUREMENT_EXPANSION,
    MODALITY_EXPANSION,
    MODEL_ARCHITECTURE_EXPANSION,
    PATHOLOGY_EXPANSION,
    SEVERITY_EXPANSION,
    expand_query,
)


# ===================================================================
# EXPANSION MAP STRUCTURE
# ===================================================================


class TestExpansionMapStructure:
    """Verify each expansion map has entries and correct types."""

    def test_modality_expansion_has_entries(self):
        assert len(MODALITY_EXPANSION) >= 8
        assert "ct" in MODALITY_EXPANSION
        assert "mri" in MODALITY_EXPANSION

    def test_body_region_expansion_has_entries(self):
        assert len(BODY_REGION_EXPANSION) >= 8
        assert "head" in BODY_REGION_EXPANSION
        assert "chest" in BODY_REGION_EXPANSION

    def test_pathology_expansion_has_entries(self):
        assert len(PATHOLOGY_EXPANSION) >= 10
        assert "hemorrhage" in PATHOLOGY_EXPANSION
        assert "nodule" in PATHOLOGY_EXPANSION

    def test_ai_task_expansion_has_entries(self):
        assert len(AI_TASK_EXPANSION) >= 7
        assert "segmentation" in AI_TASK_EXPANSION
        assert "detection" in AI_TASK_EXPANSION

    def test_severity_expansion_has_entries(self):
        assert len(SEVERITY_EXPANSION) >= 3
        assert "critical" in SEVERITY_EXPANSION

    def test_finding_expansion_has_entries(self):
        assert len(FINDING_EXPANSION) >= 5
        assert "consolidation" in FINDING_EXPANSION

    def test_guideline_expansion_has_entries(self):
        assert len(GUIDELINE_EXPANSION) >= 5
        assert "lung_rads" in GUIDELINE_EXPANSION
        assert "bi_rads" in GUIDELINE_EXPANSION

    def test_device_expansion_has_entries(self):
        assert len(DEVICE_EXPANSION) >= 3
        assert "fda" in DEVICE_EXPANSION

    def test_dataset_expansion_has_entries(self):
        assert len(DATASET_EXPANSION) >= 5
        assert "rsna" in DATASET_EXPANSION
        assert "tcia" in DATASET_EXPANSION

    def test_model_architecture_expansion_has_entries(self):
        assert len(MODEL_ARCHITECTURE_EXPANSION) >= 8
        assert "vista3d" in MODEL_ARCHITECTURE_EXPANSION
        assert "nnunet" in MODEL_ARCHITECTURE_EXPANSION

    def test_measurement_expansion_has_entries(self):
        assert len(MEASUREMENT_EXPANSION) >= 4
        assert "volume" in MEASUREMENT_EXPANSION
        assert "hounsfield" in MEASUREMENT_EXPANSION

    def test_contrast_expansion_has_entries(self):
        assert len(CONTRAST_EXPANSION) >= 3
        assert "gadolinium" in CONTRAST_EXPANSION

    def test_all_expansion_maps_has_12_entries(self):
        assert len(ALL_EXPANSION_MAPS) == 12

    def test_all_expansion_maps_values_are_sets(self):
        for expansion_map in ALL_EXPANSION_MAPS:
            for key, synonyms in expansion_map.items():
                assert isinstance(synonyms, set), (
                    f"Expansion map value for '{key}' is {type(synonyms)}, not set"
                )

    def test_modality_ct_synonyms(self):
        ct_synonyms = MODALITY_EXPANSION["ct"]
        assert "computed tomography" in ct_synonyms
        assert "cat scan" in ct_synonyms

    def test_modality_mri_synonyms(self):
        mri_synonyms = MODALITY_EXPANSION["mri"]
        assert "magnetic resonance imaging" in mri_synonyms

    def test_pathology_hemorrhage_synonyms(self):
        hem_synonyms = PATHOLOGY_EXPANSION["hemorrhage"]
        assert "bleeding" in hem_synonyms
        assert "hematoma" in hem_synonyms


# ===================================================================
# expand_query() FUNCTION
# ===================================================================


class TestExpandQuery:
    """Tests for the expand_query() function."""

    def test_ct_scan_expands_to_ct_synonyms(self):
        result = expand_query("CT scan of the head")
        assert isinstance(result, set)
        assert len(result) > 0
        # Should include CT modality expansions
        assert "computed tomography" in result or "cat scan" in result

    def test_mri_expands_to_mri_synonyms(self):
        result = expand_query("brain MRI findings")
        assert "magnetic resonance imaging" in result or "mr imaging" in result

    def test_brain_query_expands_body_region(self):
        result = expand_query("brain hemorrhage on CT")
        # 'brain' is a key in BODY_REGION_EXPANSION
        assert "cerebral" in result or "intracranial" in result

    def test_lung_nodule_expands(self):
        result = expand_query("lung nodule detection")
        # 'nodule' is a key in PATHOLOGY_EXPANSION
        assert "pulmonary nodule" in result or "lung nodule" in result

    def test_chest_query_expands(self):
        result = expand_query("chest x-ray findings")
        assert "thorax" in result or "thoracic" in result or "lung" in result

    def test_segmentation_task_expands(self):
        result = expand_query("organ segmentation model")
        assert "segment" in result or "delineation" in result or "contour" in result

    def test_hemorrhage_pathology_expands(self):
        result = expand_query("hemorrhage detection")
        assert "bleeding" in result or "hematoma" in result

    def test_no_match_returns_empty_set(self):
        result = expand_query("xyzzy foobar qwerty")
        assert isinstance(result, set)
        assert len(result) == 0

    def test_case_insensitivity_uppercase(self):
        result_lower = expand_query("ct scan")
        result_upper = expand_query("CT SCAN")
        # Both should match 'ct' key
        assert len(result_lower) > 0
        assert len(result_upper) > 0

    def test_case_insensitivity_mixed(self):
        result = expand_query("Brain MRI with contrast")
        assert len(result) > 0

    def test_multiple_maps_matched(self):
        result = expand_query("CT hemorrhage detection segmentation")
        # Should match: ct (modality), hemorrhage (pathology),
        # detection (ai_task), segmentation (ai_task)
        assert len(result) > 5

    def test_guideline_expansion(self):
        result = expand_query("lung rads classification")
        assert "lung-rads" in result or "acr lung-rads" in result

    def test_dataset_expansion(self):
        result = expand_query("RSNA hemorrhage dataset")
        assert "rsna dataset" in result or "rsna challenge" in result

    def test_model_architecture_expansion(self):
        result = expand_query("vista3d segmentation benchmark")
        assert "vista-3d" in result or "nvidia vista" in result

    def test_contrast_expansion(self):
        result = expand_query("gadolinium enhanced MRI")
        assert "gbca" in result or "mri contrast" in result

    def test_measurement_expansion(self):
        result = expand_query("volume measurement of tumor")
        assert "volumetric" in result or "cc" in result

    def test_returns_set_type(self):
        result = expand_query("any query")
        assert isinstance(result, set)
