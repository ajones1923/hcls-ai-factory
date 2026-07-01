"""
Reproducibility Manifest for the HCLS AI Factory.

After each pipeline run, generates a JSON manifest documenting:
- Software versions (all agents, pipelines, libraries)
- Model weights and embeddings used
- Reference data versions
- Hardware specifications
- Timestamps and run identifiers

Essential for clinical validation (FDA 21 CFR Part 11) and
scientific reproducibility.
"""

from __future__ import annotations

import json
import os
import platform
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class HardwareSpec:
    """Hardware specification snapshot."""
    platform: str = field(default_factory=lambda: platform.platform())
    python_version: str = field(default_factory=lambda: sys.version)
    cpu_count: int = field(default_factory=lambda: os.cpu_count() or 0)
    architecture: str = field(default_factory=lambda: platform.machine())

    # GPU info populated at runtime
    gpu_name: str = ""
    gpu_memory_gb: float = 0.0
    cuda_version: str = ""


@dataclass
class SoftwareVersion:
    """A single software component version."""
    component: str
    version: str
    source: str = ""  # e.g., "pip", "docker", "nim"


@dataclass
class ReproducibilityManifest:
    """Complete reproducibility manifest for an HCLS AI Factory run.

    This manifest captures everything needed to reproduce a pipeline run,
    supporting both scientific reproducibility and regulatory compliance.
    """
    run_id: str
    pipeline_name: str = "HCLS AI Factory"
    pipeline_version: str = "1.0.0"
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())

    # Hardware
    hardware: HardwareSpec = field(default_factory=HardwareSpec)

    # Software versions
    software: list[SoftwareVersion] = field(default_factory=list)

    # Knowledge base versions (populated from KnowledgeManifest)
    knowledge_versions: dict[str, str] = field(default_factory=dict)

    # Model weights
    embedding_model: str = "BAAI/bge-small-en-v1.5"
    embedding_dim: int = 384
    llm_model: str = ""
    nim_services: dict[str, str] = field(default_factory=dict)

    # Input data
    input_files: list[str] = field(default_factory=list)
    patient_id: str = ""

    # Output data
    output_dir: str = ""
    output_files: list[str] = field(default_factory=list)

    # Timing
    started_at: str = ""
    completed_at: str = ""
    duration_seconds: float = 0.0

    # Stages executed
    stages_executed: list[str] = field(default_factory=list)
    agents_invoked: list[str] = field(default_factory=list)

    def add_software(self, component: str, version: str, source: str = "") -> None:
        self.software.append(SoftwareVersion(component=component, version=version, source=source))

    def detect_gpu(self) -> None:
        """Attempt to detect GPU hardware."""
        try:
            import subprocess
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split(",")
                self.hardware.gpu_name = parts[0].strip()
                self.hardware.gpu_memory_gb = float(parts[1].strip()) / 1024 if len(parts) > 1 else 0.0
        except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
            pass

    def detect_python_packages(self) -> None:
        """Auto-detect installed Python packages relevant to the pipeline."""
        key_packages = [
            "pydantic", "fastapi", "streamlit", "loguru", "anthropic",
            "sentence-transformers", "pymilvus", "rdkit", "reportlab",
            "numpy", "scipy", "pandas",
        ]
        try:
            from importlib.metadata import version as pkg_version, PackageNotFoundError
            for pkg in key_packages:
                try:
                    v = pkg_version(pkg)
                    self.add_software(pkg, v, source="pip")
                except PackageNotFoundError:
                    pass
        except ImportError:
            pass

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: str | Path) -> Path:
        """Save manifest to JSON file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
        return path

    @classmethod
    def from_file(cls, path: str | Path) -> "ReproducibilityManifest":
        """Load manifest from JSON file."""
        with open(path) as f:
            data = json.load(f)
        manifest = cls(run_id=data["run_id"])
        for k, v in data.items():
            if k == "hardware":
                manifest.hardware = HardwareSpec(**v)
            elif k == "software":
                manifest.software = [SoftwareVersion(**s) for s in v]
            elif hasattr(manifest, k):
                setattr(manifest, k, v)
        return manifest
