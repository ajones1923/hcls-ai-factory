"""Reference cases for the published radiology scoring systems: findings in, category out.

Same purpose as the cardiology reference cases. The existing imaging tests exercise these
classifiers; they largely do not check an output against the number a guideline would give, so
a transposed point value or an off-by-one size threshold would pass. Each case below is derived
from the cited chart by adding up the feature points or reading the size band, so every
expected value is checkable by hand against the citation.

The boundaries are covered deliberately. A nodule of exactly 15.0 mm, or a TI-RADS total of
exactly 7, is where an inclusive/exclusive slip changes management -- FNA versus follow-up,
3-month CT versus tissue sampling -- and no range assertion can see it.
"""
import pytest

from src.models import LungRADS, TIRADS
from src.workflows.ct_chest_lung_nodule import (
    classify_solid_nodule_lung_rads,
    lung_rads_recommendation,
)
from src.workflows.thyroid_tirads import (
    calculate_tirads_points,
    points_to_tirads,
    tirads_fna_recommendation,
)


# =====================================================================
# ACR TI-RADS 2017 — Tessler FN et al., J Am Coll Radiol 2017;14(5):587-595
#   Composition: cystic/spongiform 0 · mixed cystic-solid 1 · solid 2
#   Echogenicity: anechoic 0 · hyper/isoechoic 1 · hypoechoic 2 · very hypoechoic 3
#   Shape: wider-than-tall 0 · taller-than-wide 3
#   Margin: smooth/ill-defined 0 · lobulated/irregular 2 · extra-thyroidal extension 3
#   Echogenic foci: none/large comet-tail 0 · macrocalcifications 1 ·
#                   peripheral calcifications 2 · punctate echogenic foci 3
#   Totals: 0 -> TR1 · 2 -> TR2 · 3 -> TR3 · 4-6 -> TR4 · 7+ -> TR5
# =====================================================================
def nodule(composition="solid", echogenicity="isoechoic", shape="wider_than_tall",
           margin="smooth", echogenic_foci="none"):
    return {"composition": composition, "echogenicity": echogenicity, "shape": shape,
            "margin": margin, "echogenic_foci": echogenic_foci}


class TestTiradsPointsFromThePublishedChart:
    def test_a_simple_cyst_scores_zero(self):
        """Cystic 0 + anechoic 0 + wider-than-tall 0 + smooth 0 + none 0."""
        pts = calculate_tirads_points(nodule(composition="cystic", echogenicity="anechoic"))
        assert pts == 0
        assert points_to_tirads(pts) is TIRADS.TR1

    def test_a_spongiform_nodule_scores_zero(self):
        pts = calculate_tirads_points(nodule(composition="spongiform", echogenicity="anechoic"))
        assert points_to_tirads(pts) is TIRADS.TR1

    def test_the_classic_tr5_nodule_scores_seven(self):
        """Solid 2 + very hypoechoic 3 + punctate echogenic foci... is 8; the published
        worked example of a 7 is solid 2 + very hypoechoic 3 + taller-than-wide... = 8 too.
        A clean 7: solid 2 + hypoechoic 2 + taller-than-wide 3."""
        pts = calculate_tirads_points(
            nodule(composition="solid", echogenicity="hypoechoic", shape="taller_than_wide"))
        assert pts == 7
        assert points_to_tirads(pts) is TIRADS.TR5

    def test_the_maximum_is_fourteen(self):
        pts = calculate_tirads_points(nodule(
            composition="solid", echogenicity="very_hypoechoic", shape="taller_than_wide",
            margin="extra_thyroidal_extension", echogenic_foci="punctate_echogenic_foci"))
        assert pts == 2 + 3 + 3 + 3 + 3 == 14
        assert points_to_tirads(pts) is TIRADS.TR5

    @pytest.mark.parametrize("feature,value,expected", [
        ("composition", "mixed_cystic_solid", 1),
        ("composition", "solid", 2),
        ("echogenicity", "hyperechoic", 1),
        ("echogenicity", "hypoechoic", 2),
        ("echogenicity", "very_hypoechoic", 3),
        ("shape", "taller_than_wide", 3),
        ("margin", "lobulated", 2),
        ("margin", "irregular", 2),
        ("margin", "extra_thyroidal_extension", 3),
        ("echogenic_foci", "macrocalcifications", 1),
        ("echogenic_foci", "peripheral_calcifications", 2),
        ("echogenic_foci", "punctate_echogenic_foci", 3),
    ])
    def test_each_feature_contributes_its_published_points(self, feature, value, expected):
        base = nodule(composition="cystic", echogenicity="anechoic")
        base[feature] = value
        assert calculate_tirads_points(base) == expected

    @pytest.mark.parametrize("feature,value", [
        ("margin", "ill_defined"),
        ("echogenic_foci", "large_comet_tail"),
        ("shape", "wider_than_tall"),
    ])
    def test_the_zero_point_features_really_are_zero(self, feature, value):
        """Ill-defined margin is 0, not 2 -- it is the one that reads as if it should score."""
        base = nodule(composition="cystic", echogenicity="anechoic")
        base[feature] = value
        assert calculate_tirads_points(base) == 0


class TestTiradsCategoryBoundaries:
    @pytest.mark.parametrize("points,expected", [
        (0, TIRADS.TR1), (1, TIRADS.TR2), (2, TIRADS.TR2), (3, TIRADS.TR3),
        (4, TIRADS.TR4), (6, TIRADS.TR4), (7, TIRADS.TR5), (14, TIRADS.TR5),
    ])
    def test_point_totals_map_to_the_published_levels(self, points, expected):
        assert points_to_tirads(points) is expected


class TestTiradsFnaThresholds:
    """ACR TI-RADS size thresholds: TR3 FNA >=2.5cm / follow >=1.5cm;
    TR4 FNA >=1.5 / follow >=1.0; TR5 FNA >=1.0 / follow >=0.5."""

    @pytest.mark.parametrize("tr,size_mm,wants_fna", [
        (TIRADS.TR3, 24.9, False), (TIRADS.TR3, 25.0, True),
        (TIRADS.TR4, 14.9, False), (TIRADS.TR4, 15.0, True),
        (TIRADS.TR5, 9.9, False), (TIRADS.TR5, 10.0, True),
    ])
    def test_the_fna_boundary_is_inclusive(self, tr, size_mm, wants_fna):
        rec = tirads_fna_recommendation(tr, size_mm)
        assert ("FNA recommended" in rec) is wants_fna

    @pytest.mark.parametrize("tr,size_mm", [
        (TIRADS.TR3, 15.0), (TIRADS.TR4, 10.0), (TIRADS.TR5, 5.0),
    ])
    def test_the_follow_up_boundary_is_inclusive(self, tr, size_mm):
        rec = tirads_fna_recommendation(tr, size_mm)
        assert "Follow-up ultrasound" in rec

    @pytest.mark.parametrize("tr", [TIRADS.TR1, TIRADS.TR2])
    def test_benign_levels_need_nothing_at_any_size(self, tr):
        rec = tirads_fna_recommendation(tr, 40.0)
        assert "No FNA indicated" in rec


# =====================================================================
# ACR Lung-RADS v2022 — solid nodules, BASELINE screen
#   Category 2 < 6mm · 3 is 6 to <8mm · 4A is 8 to <15mm · 4B >= 15mm
# =====================================================================
class TestLungRadsSizeBands:
    @pytest.mark.parametrize("mm,expected", [
        (3.0, LungRADS.CAT_2),
        (5.9, LungRADS.CAT_2),
        (6.0, LungRADS.CAT_3),     # lower edge of 3
        (7.9, LungRADS.CAT_3),
        (8.0, LungRADS.CAT_4A),    # lower edge of 4A
        (14.9, LungRADS.CAT_4A),
        (15.0, LungRADS.CAT_4B),   # lower edge of 4B
        (30.0, LungRADS.CAT_4B),
    ])
    def test_each_band_boundary(self, mm, expected):
        assert classify_solid_nodule_lung_rads(mm) is expected

    def test_the_management_follows_the_category(self):
        assert "annual" in lung_rads_recommendation(LungRADS.CAT_2).lower()
        assert "6-month" in lung_rads_recommendation(LungRADS.CAT_3)
        assert "3-month" in lung_rads_recommendation(LungRADS.CAT_4A)
        assert "tissue sampling" in lung_rads_recommendation(LungRADS.CAT_4B).lower()


class TestLungRadsStatesItsLimitation:
    def test_the_docstring_says_these_are_baseline_thresholds(self):
        """Size is the only input, and at follow-up a new nodule moves up a category. A
        reader has to be able to find that out without re-deriving it."""
        doc = classify_solid_nodule_lung_rads.__doc__ or ""
        assert "BASELINE" in doc
        assert "follow-up" in doc.lower()
