"""
Digital-twin / visualization endpoint (PRD §3 FR-VZ-15; Surface d). Authors an OpenUSD
scene from a patient's projection and returns the artifact path + the SceneSpec. The Spark
authors; render the `.usda` with RTX in NVIDIA Omniverse (RunPod / RTX workstation).
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request

from src.viz import author_population, author_scene

router = APIRouter(prefix="/viz", tags=["viz"])
_SCENES = {"lesion": "lesion_trajectory", "lesion_trajectory": "lesion_trajectory",
           "mosaic": "mosaic", "atlas": "atlas"}


def _strip(res: dict, inline: bool) -> dict:
    return res if inline else {k: v for k, v in res.items() if k != "usda"}


@router.get("/population")
def population(request: Request, inline: bool = False):
    app = request.app
    eng, manifest, fmap = app.state.engine, app.state.manifest, app.state.featured
    rev = {v: k for k, v in (fmap or {}).items()}
    patients = [(p["patient_id"], rev.get(p["patient_id"]), eng.store.projection(p["patient_id"]))
                for p in manifest["patients"]]
    return _strip(author_population(patients), inline)


@router.get("/{scene}/{patient_id}")
def viz(scene: str, patient_id: str, request: Request, inline: bool = False):
    if scene not in _SCENES:
        raise HTTPException(404, f"unknown scene '{scene}' (Phase V0: {sorted(set(_SCENES))} + population)")
    proj = request.app.state.engine.store.projection(patient_id)
    return _strip(author_scene(patient_id, proj, _SCENES[scene]), inline)
