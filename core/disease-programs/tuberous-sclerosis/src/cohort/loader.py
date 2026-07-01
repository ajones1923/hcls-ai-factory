"""
Cohort loader (PRD §2.6 cohort_loaded). Reads the built cohort into the `cohort`
dict the agents consume, so the orchestrator runs against version-controlled artifacts
instead of inline demo data.
"""
from __future__ import annotations

import json
from pathlib import Path

from config.settings import settings


def _root(cohort_dir: str | Path | None) -> Path:
    return Path(cohort_dir) if cohort_dir else settings.COHORT_DIR


def load_manifest(cohort_dir: str | Path | None = None) -> dict:
    return json.loads((_root(cohort_dir) / "manifest.json").read_text())


def load_patient(patient_id: str, cohort_dir: str | Path | None = None) -> dict:
    return json.loads((_root(cohort_dir) / patient_id / "patient.json").read_text())


def featured_map(cohort_dir: str | Path | None = None) -> dict[str, str]:
    return load_manifest(cohort_dir)["featured"]
