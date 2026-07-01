"""Tests for hcls_common.reproducibility — manifest creation and serialization."""

import json
import tempfile
from pathlib import Path

import pytest

from hcls_common.reproducibility import (
    HardwareSpec,
    ReproducibilityManifest,
    SoftwareVersion,
)


class TestHardwareSpec:
    def test_defaults_populated(self):
        spec = HardwareSpec()
        assert spec.platform  # non-empty string
        assert spec.python_version  # non-empty string
        assert spec.cpu_count >= 0
        assert spec.architecture  # non-empty string

    def test_gpu_fields_default_empty(self):
        spec = HardwareSpec()
        assert spec.gpu_name == ""
        assert spec.gpu_memory_gb == 0.0
        assert spec.cuda_version == ""


class TestSoftwareVersion:
    def test_creation(self):
        sv = SoftwareVersion(component="numpy", version="1.24.0", source="pip")
        assert sv.component == "numpy"
        assert sv.version == "1.24.0"
        assert sv.source == "pip"

    def test_default_source(self):
        sv = SoftwareVersion(component="pymilvus", version="2.4.0")
        assert sv.source == ""


class TestReproducibilityManifest:
    def test_creation(self):
        manifest = ReproducibilityManifest(run_id="test-run-001")
        assert manifest.run_id == "test-run-001"
        assert manifest.pipeline_name == "HCLS AI Factory"
        assert manifest.created_at  # ISO timestamp

    def test_add_software(self):
        manifest = ReproducibilityManifest(run_id="test-run-002")
        manifest.add_software("numpy", "1.24.0", source="pip")
        manifest.add_software("pymilvus", "2.4.0")
        assert len(manifest.software) == 2
        assert manifest.software[0].component == "numpy"

    def test_to_dict(self):
        manifest = ReproducibilityManifest(
            run_id="test-run-003",
            patient_id="HG002",
        )
        d = manifest.to_dict()
        assert d["run_id"] == "test-run-003"
        assert d["patient_id"] == "HG002"
        assert "hardware" in d
        assert "software" in d

    def test_save_and_load(self):
        manifest = ReproducibilityManifest(
            run_id="test-roundtrip",
            patient_id="HG002",
            llm_model="claude-sonnet-4-20250514",
        )
        manifest.add_software("numpy", "1.24.0", source="pip")

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "manifest.json"
            manifest.save(path)

            assert path.exists()
            with open(path) as f:
                data = json.load(f)
            assert data["run_id"] == "test-roundtrip"

            loaded = ReproducibilityManifest.from_file(path)
            assert loaded.run_id == "test-roundtrip"
            assert loaded.patient_id == "HG002"
            assert loaded.llm_model == "claude-sonnet-4-20250514"
            assert len(loaded.software) == 1

    def test_detect_python_packages(self):
        manifest = ReproducibilityManifest(run_id="test-detect")
        manifest.detect_python_packages()
        # numpy should be found since it is a dependency
        components = [s.component for s in manifest.software]
        assert "numpy" in components

    def test_save_creates_parent_dirs(self):
        manifest = ReproducibilityManifest(run_id="test-mkdir")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "subdir" / "nested" / "manifest.json"
            manifest.save(path)
            assert path.exists()
