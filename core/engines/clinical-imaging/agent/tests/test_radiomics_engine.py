"""Tests for RadiomicsEngine and RadiomicsIngestPipeline.

Validates mock feature extraction, multi-region extraction, longitudinal
comparison with significance flags, GPU fallback behavior, configurable
feature classes, and text summary generation for embedding.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.radiomics_engine import (
    ALL_FEATURE_CLASSES,
    DEFAULT_FEATURE_CLASSES,
    LongitudinalComparison,
    LongitudinalDelta,
    RadiomicsEngine,
    RadiomicsFeatureSet,
    _GPU_AVAILABLE,
    _MOCK_FEATURE_MAP,
)
from src.ingest.radiomics_parser import RadiomicsIngestPipeline, RadiomicsRecord


# ===================================================================
# RadiomicsEngine — Mock Feature Extraction
# ===================================================================


class TestExtractFeaturesMock:
    """Test extract_features() in mock mode."""

    def test_returns_dict(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        assert isinstance(features, dict)

    def test_contains_firstorder_features(self):
        engine = RadiomicsEngine(
            feature_classes=["firstorder"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        assert any(k.startswith("firstorder_") for k in features)
        assert "firstorder_Mean" in features
        assert "firstorder_Entropy" in features

    def test_contains_shape_features(self):
        engine = RadiomicsEngine(
            feature_classes=["shape"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        assert any(k.startswith("shape_") for k in features)
        assert "shape_VoxelVolume" in features
        assert "shape_Sphericity" in features

    def test_contains_glcm_features(self):
        engine = RadiomicsEngine(
            feature_classes=["glcm"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        assert any(k.startswith("glcm_") for k in features)
        assert "glcm_Contrast" in features
        assert "glcm_Correlation" in features

    def test_default_classes_combined(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        # Default is firstorder + shape + glcm
        assert any(k.startswith("firstorder_") for k in features)
        assert any(k.startswith("shape_") for k in features)
        assert any(k.startswith("glcm_") for k in features)

    def test_all_values_are_float(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        for key, value in features.items():
            assert isinstance(value, float), f"{key} is {type(value)}, expected float"

    def test_values_within_expected_ranges(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        for key, value in features.items():
            for cls_features in _MOCK_FEATURE_MAP.values():
                if key in cls_features:
                    lo, hi = cls_features[key]
                    assert lo <= value <= hi, (
                        f"{key}={value} outside range [{lo}, {hi}]"
                    )

    def test_deterministic_for_same_label(self):
        engine = RadiomicsEngine(mock=True)
        f1 = engine.extract_features("a.nii.gz", "m.nii.gz", label=1)
        f2 = engine.extract_features("b.nii.gz", "n.nii.gz", label=1)
        assert f1 == f2

    def test_different_for_different_labels(self):
        engine = RadiomicsEngine(mock=True)
        f1 = engine.extract_features("a.nii.gz", "m.nii.gz", label=1)
        f2 = engine.extract_features("a.nii.gz", "m.nii.gz", label=2)
        assert f1 != f2

    def test_feature_count(self):
        engine = RadiomicsEngine(
            feature_classes=["firstorder", "shape", "glcm"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        # Should have features from all three classes
        expected_count = (
            len(_MOCK_FEATURE_MAP["firstorder"])
            + len(_MOCK_FEATURE_MAP["shape"])
            + len(_MOCK_FEATURE_MAP["glcm"])
        )
        assert len(features) == expected_count


# ===================================================================
# RadiomicsEngine — Extract All Regions (Mock)
# ===================================================================


class TestExtractAllRegionsMock:
    """Test extract_all_regions() in mock mode."""

    def test_returns_dict_of_dicts(self):
        engine = RadiomicsEngine(mock=True)
        regions = engine.extract_all_regions("dummy.nii.gz", "mask.nii.gz")
        assert isinstance(regions, dict)
        for label, features in regions.items():
            assert isinstance(label, int)
            assert isinstance(features, dict)

    def test_returns_multiple_regions(self):
        engine = RadiomicsEngine(mock=True)
        regions = engine.extract_all_regions("dummy.nii.gz", "mask.nii.gz")
        assert len(regions) >= 2

    def test_each_region_has_features(self):
        engine = RadiomicsEngine(mock=True)
        regions = engine.extract_all_regions("dummy.nii.gz", "mask.nii.gz")
        for label, features in regions.items():
            assert len(features) > 0
            assert all(isinstance(v, float) for v in features.values())

    def test_regions_have_different_values(self):
        engine = RadiomicsEngine(mock=True)
        regions = engine.extract_all_regions("dummy.nii.gz", "mask.nii.gz")
        labels = list(regions.keys())
        if len(labels) >= 2:
            assert regions[labels[0]] != regions[labels[1]]


# ===================================================================
# RadiomicsEngine — Longitudinal Comparison
# ===================================================================


class TestCompareLongitudinal:
    """Test compare_longitudinal() delta calculation and significance flags."""

    def test_returns_longitudinal_comparison(self):
        engine = RadiomicsEngine(mock=True)
        f1 = engine.extract_features("a.nii.gz", "m.nii.gz", label=1)
        f2 = engine.extract_features("a.nii.gz", "m.nii.gz", label=2)
        result = engine.compare_longitudinal(f1, f2)
        assert isinstance(result, LongitudinalComparison)

    def test_delta_calculation(self):
        engine = RadiomicsEngine(mock=True)
        f1 = {"firstorder_Mean": 50.0, "firstorder_Entropy": 4.0}
        f2 = {"firstorder_Mean": 60.0, "firstorder_Entropy": 4.5}
        result = engine.compare_longitudinal(f1, f2)
        assert result.total_features_compared == 2

        mean_delta = next(d for d in result.deltas if d.feature_name == "firstorder_Mean")
        assert mean_delta.absolute_delta == 10.0
        assert mean_delta.value_t1 == 50.0
        assert mean_delta.value_t2 == 60.0
        assert mean_delta.direction == "increased"

    def test_significance_flags(self):
        engine = RadiomicsEngine(mock=True)
        # firstorder_Mean range is (30.0, 120.0), so SD ~ 22.5, 2*SD ~ 45
        f1 = {"firstorder_Mean": 40.0}
        f2 = {"firstorder_Mean": 100.0}  # delta = 60 > 45
        result = engine.compare_longitudinal(f1, f2)
        delta = result.deltas[0]
        assert delta.is_significant is True
        assert "firstorder_Mean" in result.significant_changes

    def test_non_significant_change(self):
        engine = RadiomicsEngine(mock=True)
        # firstorder_Mean range is (30.0, 120.0), SD ~ 22.5, 2*SD ~ 45
        f1 = {"firstorder_Mean": 50.0}
        f2 = {"firstorder_Mean": 55.0}  # delta = 5 < 45
        result = engine.compare_longitudinal(f1, f2)
        delta = result.deltas[0]
        assert delta.is_significant is False

    def test_stable_direction(self):
        engine = RadiomicsEngine(mock=True)
        f1 = {"firstorder_Mean": 50.0}
        f2 = {"firstorder_Mean": 50.0}
        result = engine.compare_longitudinal(f1, f2)
        delta = result.deltas[0]
        assert delta.direction == "stable"
        assert delta.absolute_delta == 0.0

    def test_decreased_direction(self):
        engine = RadiomicsEngine(mock=True)
        f1 = {"firstorder_Mean": 80.0}
        f2 = {"firstorder_Mean": 50.0}
        result = engine.compare_longitudinal(f1, f2)
        delta = result.deltas[0]
        assert delta.direction == "decreased"

    def test_percent_change(self):
        engine = RadiomicsEngine(mock=True)
        f1 = {"firstorder_Mean": 100.0}
        f2 = {"firstorder_Mean": 120.0}
        result = engine.compare_longitudinal(f1, f2)
        delta = result.deltas[0]
        assert delta.percent_change == 20.0

    def test_no_common_features(self):
        engine = RadiomicsEngine(mock=True)
        f1 = {"firstorder_Mean": 50.0}
        f2 = {"shape_Sphericity": 0.8}
        result = engine.compare_longitudinal(f1, f2)
        assert result.total_features_compared == 0

    def test_summary_generated(self):
        engine = RadiomicsEngine(mock=True)
        f1 = engine.extract_features("a.nii.gz", "m.nii.gz", label=1)
        f2 = engine.extract_features("a.nii.gz", "m.nii.gz", label=2)
        result = engine.compare_longitudinal(f1, f2)
        assert len(result.summary) > 0
        assert "Compared" in result.summary


# ===================================================================
# RadiomicsEngine — GPU Fallback
# ===================================================================


class TestGPUFallback:
    """Test GPU fallback behavior when pyradiomics_cuda is unavailable."""

    def test_gpu_not_available_uses_cpu(self):
        # The standard test environment will not have pyradiomics_cuda
        engine = RadiomicsEngine(gpu_enabled=True, mock=True)
        # _GPU_AVAILABLE is module-level; engine.gpu_enabled reflects
        # whether GPU was both requested AND available
        if not _GPU_AVAILABLE:
            assert engine.is_gpu_accelerated is False
            assert engine.backend_name == "mock"

    def test_gpu_disabled_explicitly(self):
        engine = RadiomicsEngine(gpu_enabled=False, mock=True)
        assert engine.is_gpu_accelerated is False

    def test_mock_mode_backend_name(self):
        engine = RadiomicsEngine(mock=True)
        assert engine.backend_name == "mock"

    def test_status_contains_gpu_info(self):
        engine = RadiomicsEngine(mock=True)
        status = engine.get_status()
        assert "gpu_available" in status
        assert "gpu_enabled" in status
        assert "mock_mode" in status
        assert "backend" in status
        assert status["mock_mode"] is True


# ===================================================================
# RadiomicsEngine — Configurable Feature Classes
# ===================================================================


class TestConfigurableFeatureClasses:
    """Test that only requested feature classes are extracted."""

    def test_single_class(self):
        engine = RadiomicsEngine(feature_classes=["shape"], mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        # All features should be shape_*
        for key in features:
            assert key.startswith("shape_"), f"Unexpected feature: {key}"

    def test_two_classes(self):
        engine = RadiomicsEngine(
            feature_classes=["firstorder", "glcm"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        for key in features:
            assert key.startswith("firstorder_") or key.startswith("glcm_"), (
                f"Unexpected feature: {key}"
            )

    def test_all_seven_classes(self):
        engine = RadiomicsEngine(
            feature_classes=ALL_FEATURE_CLASSES, mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        # Should have features from every class
        prefixes_found = set()
        for key in features:
            prefix = key.split("_", 1)[0]
            prefixes_found.add(prefix)
        for cls in ALL_FEATURE_CLASSES:
            assert cls in prefixes_found, f"Missing features for class: {cls}"

    def test_invalid_class_ignored(self):
        engine = RadiomicsEngine(
            feature_classes=["firstorder", "invalid_class"], mock=True
        )
        assert "firstorder" in engine.feature_classes
        assert "invalid_class" not in engine.feature_classes

    def test_no_valid_classes_uses_defaults(self):
        engine = RadiomicsEngine(
            feature_classes=["bogus", "nonsense"], mock=True
        )
        assert engine.feature_classes == list(DEFAULT_FEATURE_CLASSES)

    def test_excludes_unrequested_classes(self):
        engine = RadiomicsEngine(feature_classes=["ngtdm"], mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        for key in features:
            assert key.startswith("ngtdm_"), f"Unexpected feature: {key}"
        # Should NOT have firstorder, shape, etc.
        assert not any(k.startswith("firstorder_") for k in features)
        assert not any(k.startswith("shape_") for k in features)


# ===================================================================
# RadiomicsEngine — Feature Summary Generation
# ===================================================================


class TestFeatureSummaryGeneration:
    """Test generate_feature_summary() for embedding suitability."""

    def test_returns_string(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(features)
        assert isinstance(summary, str)

    def test_includes_region_name(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(
            features, region_name="tumor_core"
        )
        assert "tumor_core" in summary

    def test_includes_modality(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(
            features, modality="ct"
        )
        assert "ct" in summary.lower() or "Modality: ct" in summary

    def test_includes_firstorder_stats(self):
        engine = RadiomicsEngine(
            feature_classes=["firstorder"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(features)
        assert "mean intensity" in summary.lower() or "First-order" in summary

    def test_includes_shape_info(self):
        engine = RadiomicsEngine(
            feature_classes=["shape"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(features)
        assert "volume" in summary.lower() or "Shape" in summary

    def test_includes_glcm_info(self):
        engine = RadiomicsEngine(
            feature_classes=["glcm"], mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(features)
        assert "GLCM" in summary or "texture" in summary.lower()

    def test_total_feature_count_in_summary(self):
        engine = RadiomicsEngine(mock=True)
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(features)
        assert "Total radiomic features" in summary

    def test_max_length_respected(self):
        engine = RadiomicsEngine(
            feature_classes=ALL_FEATURE_CLASSES, mock=True
        )
        features = engine.extract_features("dummy.nii.gz", "mask.nii.gz")
        summary = RadiomicsEngine.generate_feature_summary(
            features, region_name="tumor", modality="mri"
        )
        assert len(summary) <= 5000

    def test_empty_features(self):
        summary = RadiomicsEngine.generate_feature_summary({})
        assert isinstance(summary, str)
        assert "Total radiomic features: 0" in summary


# ===================================================================
# RadiomicsRecord Model
# ===================================================================


class TestRadiomicsRecord:
    """Test the RadiomicsRecord Pydantic model."""

    def test_valid_record(self):
        record = RadiomicsRecord(
            id="rad-001",
            patient_id="PT-001",
            study_date="2026-04-13",
            region_label="tumor_core",
            modality="ct",
            feature_class="firstorder",
            feature_count=19,
            mean_intensity=75.3,
            volume_mm3=12500.0,
            sphericity=0.82,
            entropy=4.5,
            feature_summary="Region: tumor_core. First-order statistics extracted.",
        )
        assert record.id == "rad-001"
        assert record.feature_count == 19

    def test_to_embedding_text(self):
        record = RadiomicsRecord(
            id="rad-001",
            patient_id="PT-001",
            region_label="tumor_core",
            modality="ct",
            feature_class="firstorder",
            feature_summary="Region: tumor_core. Mean intensity 75.3.",
        )
        text = record.to_embedding_text()
        assert "tumor_core" in text
        assert "ct" in text.lower() or "Modality: ct" in text

    def test_to_embedding_text_empty(self):
        record = RadiomicsRecord(id="rad-empty")
        text = record.to_embedding_text()
        assert isinstance(text, str)
        assert len(text) > 0


# ===================================================================
# RadiomicsIngestPipeline
# ===================================================================


class TestRadiomicsIngestPipeline:
    """Test the RadiomicsIngestPipeline fetch/parse cycle."""

    def test_fetch_returns_dict(self, mock_embedder, mock_collection_manager):
        engine = RadiomicsEngine(mock=True)
        pipeline = RadiomicsIngestPipeline(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            radiomics_engine=engine,
        )
        raw = pipeline.fetch(
            image_path="dummy.nii.gz",
            mask_path="mask.nii.gz",
            patient_id="PT-001",
            modality="ct",
        )
        assert isinstance(raw, dict)
        assert "regions" in raw
        assert "patient_id" in raw
        assert raw["patient_id"] == "PT-001"

    def test_parse_returns_records(self, mock_embedder, mock_collection_manager):
        engine = RadiomicsEngine(mock=True)
        pipeline = RadiomicsIngestPipeline(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            radiomics_engine=engine,
        )
        raw = pipeline.fetch(
            image_path="dummy.nii.gz",
            mask_path="mask.nii.gz",
            patient_id="PT-001",
            modality="ct",
        )
        records = pipeline.parse(raw)
        assert isinstance(records, list)
        assert len(records) > 0
        for record in records:
            assert isinstance(record, RadiomicsRecord)
            assert record.patient_id == "PT-001"
            assert record.modality == "ct"

    def test_records_have_unique_ids(self, mock_embedder, mock_collection_manager):
        engine = RadiomicsEngine(mock=True)
        pipeline = RadiomicsIngestPipeline(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            radiomics_engine=engine,
        )
        raw = pipeline.fetch(
            image_path="dummy.nii.gz",
            mask_path="mask.nii.gz",
            patient_id="PT-001",
        )
        records = pipeline.parse(raw)
        ids = [r.id for r in records]
        assert len(ids) == len(set(ids)), "Record IDs are not unique"

    def test_collection_name(self, mock_embedder, mock_collection_manager):
        pipeline = RadiomicsIngestPipeline(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
        )
        assert pipeline.COLLECTION_NAME == "imaging_radiomics"
