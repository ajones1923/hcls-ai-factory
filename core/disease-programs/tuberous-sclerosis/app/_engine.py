"""Shared engine bootstrap for the three Streamlit surfaces (PRD §2.4)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # engine root on path

from config.settings import settings                       # noqa: E402
from src.cohort.build import build_cohort                  # noqa: E402
from src.cohort.loader import featured_map, load_manifest, load_patient  # noqa: E402
from src.orchestrator.graph import Orchestrator            # noqa: E402
from src.orchestrator.state import SqliteEventStore        # noqa: E402


def get_engine():
    """Build the cohort if needed, enroll all patients, return (orchestrator, manifest)."""
    if not (settings.COHORT_DIR / "manifest.json").exists():
        build_cohort()
    manifest = load_manifest()
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    for p in manifest["patients"]:
        orch.enroll(p["patient_id"], load_patient(p["patient_id"]))
    return orch, manifest


def featured() -> dict[str, str]:
    return featured_map()
