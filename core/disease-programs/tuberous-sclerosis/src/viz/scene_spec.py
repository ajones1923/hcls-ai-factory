"""
SceneSpec — the engine-agnostic intermediate representation (FR-VZ-1).

A scene captured as plain data: lesions with anatomical anchors and per-time-code keyframes
(radius, threshold-state colour, and the 50%/90% uncertainty-envelope radii), plus the
threshold surface, provenance, and the SYNTHETIC watermark. No USD or rendering concept
lives here — the exporters (ASCII or pxr) consume this. Every numeric traces to a
projection field (the no-new-claim invariant).
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Keyframe:
    time_code: float          # = month (negative = past, positive = forecast)
    month: float
    radius: float             # lesion sphere radius (cm at 1:1 scale)
    state: str                # below | approaching | at_or_above
    color: tuple              # rgb for the state
    env50_radius: float       # outer radius of the 50% prediction-interval shell
    env90_radius: float       # outer radius of the 90% prediction-interval shell
    forecast: bool            # observed vs forecast horizon


@dataclass
class LesionScene:
    label: str
    location: str
    unit: str
    threshold: float
    anchor: tuple
    context: str
    context_radius: float
    threshold_radius: float
    crossing_grade: str
    keyframes: list[Keyframe] = field(default_factory=list)
    provenance: dict = field(default_factory=dict)


@dataclass
class Cell:
    pos: tuple
    variant: bool             # True = carries the recovered mosaic variant


@dataclass
class MosaicScene:
    """The mosaic 'powers-of-ten' scene (FR-VZ-17): a cellular field in which exactly the
    recovered VAF fraction of cells carries the variant — '8.3% VAF' made countable."""
    patient_id: str
    watermark: str
    gene: str
    variant: str
    vaf: float
    classification: str
    acmg_rule: str
    recovered: bool
    criteria: list            # [{code, bucket, rationale}]
    n_cells: int
    n_variant: int
    cells: list               # list[Cell]

    def to_dict(self) -> dict:
        return {"patient_id": self.patient_id, "scene_kind": "mosaic", "gene": self.gene,
                "variant": self.variant, "vaf": self.vaf, "classification": self.classification,
                "recovered": self.recovered, "n_cells": self.n_cells, "n_variant": self.n_variant,
                "observed_fraction": round(self.n_variant / self.n_cells, 4) if self.n_cells else 0.0,
                "criteria": self.criteria}


@dataclass
class OrganMark:
    organ: str
    pos: tuple
    rgb: tuple
    involved: bool            # lit if the patient's phenome implies involvement


@dataclass
class BodyFigure:
    """One whole-child figure (Scene 3) — body coloured by ACMG class, organs lit by the
    phenome profile, a gold halo if a mosaic variant was recovered."""
    patient_id: str
    featured: str | None
    pos: tuple                # placement (origin for the atlas; grid cell for population)
    body_rgb: tuple
    classification: str
    recovered: bool
    organs: list              # list[OrganMark] (all canonical organs, with involved flag)
    overdue: list             # surveillance types overdue
    burden: int               # count of involved organ systems


@dataclass
class AtlasScene:
    patient_id: str
    watermark: str
    figure: BodyFigure

    def to_dict(self) -> dict:
        f = self.figure
        return {"patient_id": self.patient_id, "scene_kind": "atlas",
                "classification": f.classification, "recovered": f.recovered,
                "organs_involved": [o.organ for o in f.organs if o.involved],
                "overdue": f.overdue, "burden": f.burden}


@dataclass
class PopulationScene:
    watermark: str
    n_patients: int
    n_recovered: int
    cols: int
    spacing: float
    figures: list             # list[BodyFigure]
    distributions: dict

    def to_dict(self) -> dict:
        return {"scene_kind": "population", "n_patients": self.n_patients,
                "n_recovered": self.n_recovered, "distributions": self.distributions,
                "figures": [{"patient_id": f.patient_id, "featured": f.featured,
                             "classification": f.classification, "recovered": f.recovered,
                             "burden": f.burden,
                             "organs_involved": [o.organ for o in f.organs if o.involved]}
                            for f in self.figures]}


@dataclass
class SceneSpec:
    patient_id: str
    scene_kind: str
    watermark: str
    scale: float
    start_tc: float
    end_tc: float
    time_codes_per_second: float
    lesions: list[LesionScene] = field(default_factory=list)
    provenance: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "patient_id": self.patient_id, "scene_kind": self.scene_kind,
            "watermark": self.watermark, "scale": self.scale,
            "start_tc": self.start_tc, "end_tc": self.end_tc,
            "time_codes_per_second": self.time_codes_per_second,
            "lesions": [{
                "label": le.label, "location": le.location, "unit": le.unit,
                "threshold": le.threshold, "crossing_grade": le.crossing_grade,
                "threshold_radius": le.threshold_radius,
                "keyframes": [k.__dict__ for k in le.keyframes],
            } for le in self.lesions],
        }
