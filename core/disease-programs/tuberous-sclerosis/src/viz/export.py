"""
Export façade (FR-VZ-10/15). Builds the SceneSpec from a projection and writes a `.usda`
scene. Uses the dependency-free ASCII writer (always available on the Spark); a pxr-based
writer can be added for USDZ/material polish without changing this interface.
"""
from __future__ import annotations

from pathlib import Path

from config.settings import settings
from src.viz.build_spec import (
    build_atlas_scene, build_lesion_scene, build_mosaic_scene, build_population_scene,
)
from src.viz.usd_ascii import (
    write_atlas_usda, write_mosaic_usda, write_population_usda, write_usda,
)

_OUT = Path(settings.DATA_DIR) / "usd"
_WRITERS = {"lesion_trajectory": write_usda, "mosaic": write_mosaic_usda, "atlas": write_atlas_usda}


def scene_for(patient_id: str, projection: dict, scene_kind: str = "lesion_trajectory"):
    if scene_kind == "lesion_trajectory":
        return build_lesion_scene(patient_id, projection)
    if scene_kind == "mosaic":
        return build_mosaic_scene(patient_id, projection)
    if scene_kind == "atlas":
        return build_atlas_scene(patient_id, projection)
    raise ValueError(f"scene_kind '{scene_kind}' not implemented in Phase V0")


def author_scene(patient_id: str, projection: dict, scene_kind: str = "lesion_trajectory",
                 out_dir: Path | None = None) -> dict:
    """Author and write a per-patient `.usda`. Returns {path, usda, spec}."""
    spec = scene_for(patient_id, projection, scene_kind)
    usda = _WRITERS[scene_kind](spec)
    out_dir = Path(out_dir) if out_dir else _OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{patient_id}_{scene_kind}.usda"
    path.write_text(usda, encoding="utf-8")
    return {"patient_id": patient_id, "scene_kind": scene_kind, "path": str(path),
            "usda": usda, "spec": spec.to_dict()}


def author_population(patients: list, out_dir: Path | None = None) -> dict:
    """Author the cohort population scene. `patients` = [(pid, featured, projection)]."""
    spec = build_population_scene(patients)
    usda = write_population_usda(spec)
    out_dir = Path(out_dir) if out_dir else _OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "cohort_population.usda"
    path.write_text(usda, encoding="utf-8")
    return {"scene_kind": "population", "path": str(path), "usda": usda, "spec": spec.to_dict()}
