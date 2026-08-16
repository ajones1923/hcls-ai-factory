"""CAD-RADS 2.0 grading, vessel accounting and therapy-gap detection.

Every test here pins a defect that was live in the demo. They are written against the guideline
rather than against the previous behaviour, so if the thresholds drift back the tests fail rather
than the screen quietly showing a wrong grade.

Reference: Cury et al., "CAD-RADS 2.0", J Cardiovasc Comput Tomogr 2022.
Statin intensity bands: 2018 AHA/ACC/multisociety cholesterol guideline, Table 3.
"""
import json
import pathlib

import pytest

from src.models import CADRADS
from src.workflows.ct_coronary_angiography import (
    CTCoronaryAngiographyWorkflow,
    cadrads_report,
    classify_stenosis_cadrads,
    parent_vessel,
    plaque_burden_modifier,
    statin_intensity,
)

DEMO_CASES = (pathlib.Path(__file__).resolve().parents[1]
              / "data" / "reference" / "demo_cases.json")



# data/ is deliberately never published (.gitignore; CLAUDE.md: "Data / weights / secrets stay
# local — only code + docs publish"). These tests read that local-only reference data, so on a
# fresh clone it is absent by design. Skip rather than fail — a clean clone should report a green
# run with skips, not errors that read as a broken project.
pytestmark = pytest.mark.skipif(
    not DEMO_CASES.is_file(),
    reason="local-only reference data under data/ is not published (see .gitignore)",
)

# ---------------------------------------------------------------------------------------------
# Stenosis category
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("pct,expected", [
    (0, CADRADS.CAD_0), (15, CADRADS.CAD_1), (40, CADRADS.CAD_2),
    (60, CADRADS.CAD_3), (77.9, CADRADS.CAD_4A), (100, CADRADS.CAD_5),
])
def test_stenosis_bands(pct, expected):
    assert classify_stenosis_cadrads(pct) is expected


def test_three_moderate_lesions_are_not_4b():
    """The regression this exists for.

    4B previously fired on three vessels at >= 50%. CAD-RADS 2.0 defines three-vessel disease at
    the OBSTRUCTIVE threshold (>= 70%), so three moderate lesions are CAD-RADS 3 -- consider
    functional assessment -- not an urgent referral for invasive angiography.
    """
    assert classify_stenosis_cadrads(60, 0, num_vessels_obstructive=0) is CADRADS.CAD_3


def test_three_obstructive_vessels_are_4b():
    assert classify_stenosis_cadrads(75, 0, num_vessels_obstructive=3) is CADRADS.CAD_4B


def test_left_main_boundary_is_inclusive():
    """Exactly 50% left main is significant; a strict > 50% test let the boundary case slip."""
    assert classify_stenosis_cadrads(50, left_main_stenosis_pct=50) is CADRADS.CAD_4B


# ---------------------------------------------------------------------------------------------
# Vessel accounting
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("segment,artery", [
    ("LAD", "LAD"), ("LAD-mid", "LAD"), ("D1", "LAD"),
    ("LCx", "LCx"), ("OM", "LCx"), ("RCA", "RCA"), ("PDA", "RCA"),
    ("Left Main", "Left Main"),
])
def test_segments_fold_onto_their_parent_artery(segment, artery):
    assert parent_vessel(segment) == artery


def test_two_lad_segments_count_as_one_vessel():
    """Counting rows made 4 segments look like 4 vessels, and could have fired the 3-vessel rule
    off two segments of the same LAD."""
    segs = ["LAD", "LAD-mid", "LCx", "RCA"]
    assert len({parent_vessel(s) for s in segs}) == 3


# ---------------------------------------------------------------------------------------------
# Modifiers
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("agatston,band", [
    (0, "P0"), (1, "P1"), (100, "P1"), (101, "P2"),
    (300, "P2"), (301, "P3"), (385, "P3"), (999, "P3"), (1000, "P4"),
])
def test_plaque_burden_bands(agatston, band):
    assert plaque_burden_modifier(agatston) == band


def test_report_carries_modifiers():
    assert cadrads_report(CADRADS.CAD_4A, 385, True) == "CAD-RADS 4A/P3/HRP"
    assert cadrads_report(CADRADS.CAD_4A, 385, False) == "CAD-RADS 4A/P3"


# ---------------------------------------------------------------------------------------------
# Therapy intensity
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("med,intensity", [
    ("atorvastatin_20mg", "moderate"), ("atorvastatin_40mg", "high"),
    ("atorvastatin_80mg", "high"), ("rosuvastatin_5mg", "moderate"),
    ("rosuvastatin_20mg", "high"), ("simvastatin_10mg", "low"),
])
def test_statin_intensity_bands(med, intensity):
    assert statin_intensity(med)["intensity"] == intensity


def test_non_statins_are_ignored():
    assert statin_intensity("aspirin_81mg") is None
    assert statin_intensity("") is None


# ---------------------------------------------------------------------------------------------
# End to end, on the shipped demo case
# ---------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def demo_result():
    case = next(c for c in json.loads(DEMO_CASES.read_text())
                if c["case_id"] == "DEMO-003")
    return CTCoronaryAngiographyWorkflow(
        mock_mode=True, mock_overrides=case["workflow_overrides"]).run()


def test_demo_case_reports_full_cadrads(demo_result):
    assert demo_result.classification == "CAD-RADS 4A/P3/HRP"


def test_demo_case_counts_three_arteries_not_four_rows(demo_result):
    assert demo_result.measurements["num_vessels_diseased"] == 3.0
    assert demo_result.measurements["num_vessels_obstructive_gte70"] == 1.0


def test_measured_stenosis_survives_the_demo_case_override(demo_result):
    """demo_cases.json merges with dict.update(), which used to discard the measurement."""
    assert demo_result.measurements["max_stenosis_pct"] == pytest.approx(77.9)
    assert demo_result.measurements["lad_stenosis_pct"] == pytest.approx(77.9)


def test_therapy_gap_is_surfaced(demo_result):
    gaps = [f for f in demo_result.findings if "Therapy gap" in f["description"]]
    assert len(gaps) == 1
    assert "moderate-intensity" in gaps[0]["description"]
    assert "not a prescribing instruction" in gaps[0]["description"]


def test_no_stenosis_without_plaque(demo_result):
    """A narrowed vessel with 'plaque type: none' is self-contradictory -- stenosis IS plaque."""
    for f in demo_result.findings:
        d = f["description"]
        if "% stenosis" in d and "plaque type:" in d:
            pct = float(d.split("): ")[1].split("%")[0])
            if pct > 0:
                assert "plaque type: none" not in d, d


def test_ejection_fraction_is_labelled_an_estimate(demo_result):
    """Standard coronary CTA is single-phase; LV function needs a functional acquisition."""
    assert "ejection_fraction_pct_estimate" in demo_result.measurements
    assert "ejection_fraction_pct" not in demo_result.measurements
