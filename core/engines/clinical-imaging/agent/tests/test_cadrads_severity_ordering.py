"""Two CAD-RADS 2.0 grading defects, both of which moved a severity.

1. The 4B test ran BEFORE the total-occlusion test. A 100% occlusion in a patient who also
   had left-main disease or three-vessel obstructive disease could never reach category 5 —
   it returned 4B. `cadrads_to_severity` maps 4B to URGENT and 5 to CRITICAL, so the more
   extensive disease silently downgraded the severity of the occlusion.

2. The HRP modifier fired on a single high-risk feature. HRP requires at least two. The
   primary document (Cury et al., CAD-RADS 2.0, JACC Cardiovasc Imaging 2022) was not
   retrieved for this change; two independent secondary sources agree on the >=2 rule, and
   the entry is logged in CLINICAL_REVIEW.md for confirmation.
"""
import pytest

from src.models import FindingSeverity
from src.workflows.ct_coronary_angiography import (
    CADRADS, HRP_MIN_FEATURES, cadrads_report, cadrads_to_severity,
    classify_stenosis_cadrads, high_risk_plaque_applies,
)


class TestTotalOcclusionOutranksExtent:
    def test_occlusion_with_left_main_disease_is_category_5(self):
        assert classify_stenosis_cadrads(
            max_stenosis_pct=100, left_main_stenosis_pct=60) is CADRADS.CAD_5

    def test_occlusion_with_three_vessel_disease_is_category_5(self):
        assert classify_stenosis_cadrads(
            max_stenosis_pct=100, num_vessels_obstructive=3) is CADRADS.CAD_5

    def test_the_severity_is_critical_not_urgent(self):
        """This is what the ordering cost: CRITICAL demoted to URGENT."""
        cat = classify_stenosis_cadrads(max_stenosis_pct=100, left_main_stenosis_pct=60)
        assert cadrads_to_severity(cat) is FindingSeverity.CRITICAL

    def test_4b_still_applies_without_an_occlusion(self):
        assert classify_stenosis_cadrads(
            max_stenosis_pct=85, left_main_stenosis_pct=55) is CADRADS.CAD_4B
        assert classify_stenosis_cadrads(
            max_stenosis_pct=85, num_vessels_obstructive=3) is CADRADS.CAD_4B

    def test_the_4b_threshold_is_still_70_percent_not_50(self):
        """Three moderate 50-69% lesions are CAD-RADS 3, not 4B."""
        assert classify_stenosis_cadrads(
            max_stenosis_pct=60, num_vessels_obstructive=0) is CADRADS.CAD_3

    def test_the_ordinary_ladder_is_unchanged(self):
        for pct, expected in [(0, CADRADS.CAD_0), (10, CADRADS.CAD_1), (30, CADRADS.CAD_2),
                              (60, CADRADS.CAD_3), (85, CADRADS.CAD_4A), (100, CADRADS.CAD_5)]:
            assert classify_stenosis_cadrads(max_stenosis_pct=pct) is expected

    def test_the_extent_is_not_lost_from_the_report(self):
        """5 outranks 4B, but left-main extent is the clinically consequential part."""
        report = cadrads_report(CADRADS.CAD_5, 385, True,
                                extent_note="also meets CAD-RADS 4B extent: left main >=50%")
        assert report.startswith("CAD-RADS 5/")
        assert "left main >=50%" in report


class TestHrpNeedsTwoFeatures:
    def test_one_feature_is_not_hrp(self):
        assert high_risk_plaque_applies(["napkin-ring sign"]) is False

    def test_two_features_are_hrp(self):
        assert high_risk_plaque_applies(["napkin-ring sign", "positive remodeling"]) is True

    @pytest.mark.parametrize("features", [None, []])
    def test_no_features_is_not_hrp(self, features):
        assert high_risk_plaque_applies(features) is False

    def test_the_threshold_is_two(self):
        assert HRP_MIN_FEATURES == 2

    def test_the_modifier_is_absent_from_the_report_below_the_threshold(self):
        one = cadrads_report(CADRADS.CAD_4A, 385, high_risk_plaque_applies(["napkin-ring sign"]))
        assert "HRP" not in one

    def test_the_modifier_is_present_at_or_above_the_threshold(self):
        two = cadrads_report(CADRADS.CAD_4A, 385,
                             high_risk_plaque_applies(["napkin-ring sign", "low-attenuation plaque"]))
        assert two.endswith("/HRP")


class TestTheRepresentativeCaseStillReportsHrp:
    def test_three_features_in_the_demo_case_still_fire(self):
        """The demo case has three features, so the headline it produces does not change."""
        from src.workflows.ct_coronary_angiography import CTCoronaryAngiographyWorkflow
        rep = CTCoronaryAngiographyWorkflow._representative_inference(
            CTCoronaryAngiographyWorkflow.__new__(CTCoronaryAngiographyWorkflow))
        assert high_risk_plaque_applies(rep["high_risk_plaque_features"]) is True
