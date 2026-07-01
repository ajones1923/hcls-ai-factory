"""
Stylized anatomical anchors (FR-VZ-6). A lesion is placed at its canonical landmark in a
low-detail reference geometry and scaled to the engine's measured size. This is an ATLAS,
not patient imaging — the same scene schema accepts segmented patient anatomy in Phase V3
by swapping this layer. Coordinates are in scene units (cm at 1:1 scale).
"""
from __future__ import annotations

# anchor position + the radius of the stylized context organ drawn around it
ANCHORS = {
    "foramen of Monro": {"pos": (0.0, 0.0, 0.0), "context": "brain", "context_radius": 6.0},
    "kidney": {"pos": (0.0, -2.0, 0.0), "context": "kidney", "context_radius": 5.0},
    "cortex": {"pos": (0.0, 3.5, 0.0), "context": "brain", "context_radius": 6.0},
    "clinical": {"pos": (0.0, 0.0, 0.0), "context": "body", "context_radius": 8.0},
}
_DEFAULT = {"pos": (0.0, 0.0, 0.0), "context": "body", "context_radius": 6.0}


def anchor_for(location: str) -> dict:
    return ANCHORS.get(location, _DEFAULT)


# ── whole-child organ atlas (Scene 3) ───────────────────────────────────────
# canonical TSC-affected organ systems: body-local position + display colour
ORGANS = {
    "brain": {"pos": (0.0, 4.2, 0.0), "rgb": (0.55, 0.45, 0.95)},
    "heart": {"pos": (-0.6, 1.7, 0.5), "rgb": (0.90, 0.30, 0.35)},
    "lung": {"pos": (0.7, 2.3, 0.4), "rgb": (0.45, 0.70, 0.85)},
    "kidney": {"pos": (0.9, -0.4, -0.3), "rgb": (0.85, 0.65, 0.30)},
    "skin": {"pos": (0.0, 1.0, 1.1), "rgb": (0.95, 0.80, 0.55)},
}

# HPO term -> affected organ system (drives organ illumination from the phenome profile)
HPO_ORGAN = {
    "HP:0009717": "brain",    # cortical tubers
    "HP:0009718": "brain",    # SEGA
    "HP:0001250": "brain",    # seizure
    "HP:0012469": "brain",    # infantile spasms
    "HP:0001249": "brain",    # intellectual disability
    "HP:0000717": "brain",    # autism
    "HP:0006772": "kidney",   # renal angiomyolipoma
    "HP:0009719": "skin",     # hypomelanotic macule
    "HP:0009720": "skin",     # adenoma sebaceum
    "HP:0009721": "skin",     # shagreen patch
    "HP:0030679": "skin",     # ash-leaf spot
    "HP:0030680": "heart",    # cardiovascular (rhabdomyoma, if present)
    "HP:0002088": "lung",     # abnormal lung (LAM, if present)
}
