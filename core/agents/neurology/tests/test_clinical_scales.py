"""Tests for clinical scale calculators.

Tests all 10 scale calculators at their boundary values:
  NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn & Yahr

Author: Adam Jones
Date: March 2026

CORRECTED (2026-11-03). This file claimed to test "all 10 scale calculators with boundary
values" and tested none of them. Every test constructed a `ScaleResult` pydantic model with a
hand-written `score` and asserted `result.score == <that same number>` -- `test_nihss_severe_
boundary` asserted `21 == 21`. They were field-assignment tests: they passed identically with
all ten calculators deleted, which is why they caught none of N1, N4, N5, N8, N9, N17 or N20.

Two further problems were fixed with them:

* Line 15 imported `calculate_scale` from `src.clinical_scales`, and no such name existed
  anywhere in the repo, so `_CS_AVAILABLE` was always False -- and `_CS_AVAILABLE` was never
  referenced by any `skipif`, so the flag was dead and nothing was reported as skipped. The
  name now exists (it is the single shared entry point the route-local `_calculate_*`
  functions delegate to), so the import is real and the flag is gone.
* `test_edss_max_score` "verified" that EDSS 10.0 means death by typing the string "Death due
  to MS" into a model it then read back, while the calculator was labelling 9.5 as death
  (N8). It now asserts against the calculator, at both 9.5 and 10.0.

Each test below keeps the boundary value the original named and asserts the calculator's
output at it. The behaviour-specific reference cases for the defects live in
`test_scale_data_regressions.py`.
"""
import pytest

from src.clinical_scales import (
    ALSFRSCalculator,
    ASPECTSCalculator,
    EDSSCalculator,
    GCSCalculator,
    HIT6Calculator,
    HoehnYahrCalculator,
    MoCACalculator,
    NIHSSCalculator,
    ScaleDataError,
    UPDRSCalculator,
    calculate_scale,
    mRSCalculator,
)
from src.knowledge import CLINICAL_SCALES
from src.models import ClinicalScaleType


def _all_items(item_maxima, value=0):
    """A COMPLETE item set, so the calculator returns a score rather than 'incomplete'."""
    return {item: value for item in item_maxima}


def _nihss_totalling(total):
    scores = _all_items(NIHSSCalculator.ITEMS)
    remaining = total
    for item, max_val in NIHSSCalculator.ITEMS.items():
        take = min(max_val, remaining)
        scores[item] = take
        remaining -= take
    assert remaining == 0, f"NIHSS cannot reach {total}"
    return scores


class TestNIHSSScale:
    """Tests for NIHSS scale boundaries."""

    def test_nihss_exists_in_knowledge(self):
        assert "nihss" in CLINICAL_SCALES

    def test_nihss_min_score_valid(self):
        result = NIHSSCalculator.calculate(_all_items(NIHSSCalculator.ITEMS))
        assert result.score == 0
        assert result.max_score == 42
        assert result.scale_type == ClinicalScaleType.NIHSS
        assert result.severity_category == "No stroke symptoms"

    def test_nihss_max_score_valid(self):
        result = NIHSSCalculator.calculate(dict(NIHSSCalculator.ITEMS))
        assert result.score == 42
        assert result.severity_category == "Very severe stroke"

    def test_nihss_moderate_boundary(self):
        result = NIHSSCalculator.calculate(_nihss_totalling(5))
        assert result.score == 5
        assert result.severity_category == "Moderate stroke"

    def test_nihss_severe_boundary(self):
        result = NIHSSCalculator.calculate(_nihss_totalling(21))
        assert result.score == 21
        assert result.severity_category == "Severe stroke"


class TestGCSScale:
    """Tests for GCS scale boundaries."""

    def test_gcs_exists_in_knowledge(self):
        assert "gcs" in CLINICAL_SCALES

    def test_gcs_min_score(self):
        result = GCSCalculator.calculate(eye=1, verbal=1, motor=1)
        assert result.score == 3
        assert result.severity_category == "Severe"

    def test_gcs_max_score(self):
        result = GCSCalculator.calculate(eye=4, verbal=5, motor=6)
        assert result.score == 15
        assert result.max_score == 15
        assert result.severity_category == "Mild"

    def test_gcs_moderate_boundary(self):
        result = GCSCalculator.calculate(eye=2, verbal=3, motor=4)
        assert result.score == 9
        assert result.severity_category == "Moderate"

    def test_a_missing_component_is_refused_not_defaulted_to_1(self):
        """The route twin defaulted each absent component to 1 -- the worst value for it."""
        with pytest.raises(ScaleDataError):
            calculate_scale("gcs", {"eye": 4, "verbal": 5})


class TestMoCAScale:
    """Tests for MoCA scale boundaries."""

    def test_moca_exists_in_knowledge(self):
        assert "moca" in CLINICAL_SCALES

    def test_moca_normal_score(self):
        domains = dict(MoCACalculator.DOMAINS)
        domains["delayed_recall"] = 3  # 30 - 2 = 28
        result = MoCACalculator.calculate(domains, education_years=16)
        assert result.score == 28
        assert result.severity_category == "Normal"

    def test_moca_mci_boundary(self):
        domains = dict(MoCACalculator.DOMAINS)
        domains["delayed_recall"] = 0  # 30 - 5 = 25
        result = MoCACalculator.calculate(domains, education_years=16)
        assert result.score == 25
        assert result.severity_category == "Mild cognitive impairment"

    def test_moca_min_score(self):
        result = MoCACalculator.calculate(
            _all_items(MoCACalculator.DOMAINS), education_years=16
        )
        assert result.score == 0
        assert result.max_score == 30
        assert result.severity_category == "Severe cognitive impairment"


class TestUPDRSScale:
    """Tests for UPDRS Part III scale boundaries."""

    def test_updrs_exists_in_knowledge(self):
        assert "updrs" in CLINICAL_SCALES

    def test_updrs_min_score(self):
        result = UPDRSCalculator.calculate(_all_items(UPDRSCalculator.ITEMS))
        assert result.score == 0
        assert result.severity_category == "Minimal"

    def test_updrs_max_score(self):
        result = UPDRSCalculator.calculate(_all_items(UPDRSCalculator.ITEMS, 4))
        assert result.score == 132
        assert result.max_score == 132
        assert result.severity_category == "Very severe"


class TestEDSSScale:
    """Tests for EDSS scale boundaries."""

    def test_edss_exists_in_knowledge(self):
        assert "edss" in CLINICAL_SCALES

    def test_edss_min_score(self):
        result = EDSSCalculator.calculate({}, 0.0)
        assert result.score == 0.0
        assert result.severity_category == "Minimal disability"

    def test_edss_max_score(self):
        """EDSS 10.0, and only 10.0, is death. The band used to start at 9.5 (N8)."""
        result = EDSSCalculator.calculate({}, 10.0)
        assert result.score == 10.0
        assert result.max_score == 10.0
        assert result.severity_category == "Death due to MS"

    def test_edss_9_5_is_not_death(self):
        result = EDSSCalculator.calculate({}, 9.5)
        assert "death" not in result.severity_category.lower()

    def test_edss_wheelchair_boundary(self):
        result = EDSSCalculator.calculate({}, 7.0)
        assert result.score == 7.0
        assert result.severity_category == "Wheelchair dependent"


class TestMRSScale:
    """Tests for Modified Rankin Scale boundaries."""

    def test_mrs_exists_in_knowledge(self):
        assert "mrs" in CLINICAL_SCALES

    def test_mrs_min_score(self):
        result = mRSCalculator.calculate(0)
        assert result.score == 0
        assert result.severity_category == "No symptoms"

    def test_mrs_max_score(self):
        result = mRSCalculator.calculate(6)
        assert result.score == 6
        assert result.max_score == 6
        assert result.severity_category == "Dead"


class TestHIT6Scale:
    """Tests for HIT-6 scale boundaries."""

    def test_hit6_exists_in_knowledge(self):
        assert "hit6" in CLINICAL_SCALES

    def test_hit6_min_score(self):
        result = HIT6Calculator.calculate(["never"] * 6)
        assert result.score == 36
        assert result.severity_category == "Little or no impact"

    def test_hit6_severe_boundary(self):
        # never=6, always=13: four "always" + two "rarely" = 52 + 16 = 68; tune to 60.
        result = HIT6Calculator.calculate(
            ["always", "always", "always", "always", "sometimes", "never"]
        )
        assert result.score == 68
        assert result.max_score == 78
        assert result.severity_category == "Severe impact"

    def test_hit6_substantial_impact_boundary(self):
        result = HIT6Calculator.calculate(
            ["sometimes", "sometimes", "sometimes", "sometimes", "rarely", "rarely"]
        )
        assert result.score == 56
        assert result.severity_category == "Substantial impact"


class TestALSFRSScale:
    """Tests for ALSFRS-R scale boundaries."""

    def test_alsfrs_exists_in_knowledge(self):
        assert "alsfrs" in CLINICAL_SCALES

    def test_alsfrs_max_score(self):
        result = ALSFRSCalculator.calculate({item: 4 for item in ALSFRSCalculator.ITEMS})
        assert result.score == 48
        assert result.max_score == 48
        assert result.severity_category == "Mild functional impairment"

    def test_alsfrs_min_score(self):
        result = ALSFRSCalculator.calculate({item: 0 for item in ALSFRSCalculator.ITEMS})
        assert result.score == 0
        assert result.severity_category == "Severe functional impairment"


class TestASPECTSScale:
    """Tests for ASPECTS scale boundaries."""

    def test_aspects_exists_in_knowledge(self):
        assert "aspects" in CLINICAL_SCALES

    def test_aspects_max_score(self):
        result = ASPECTSCalculator.calculate([])
        assert result.score == 10
        assert result.max_score == 10
        assert result.severity_category == "Favorable for intervention"

    def test_aspects_thrombectomy_boundary(self):
        result = ASPECTSCalculator.calculate(["C", "L", "IC", "I"])
        assert result.score == 6
        assert result.severity_category == "Intermediate"

    def test_aspects_min_score(self):
        result = ASPECTSCalculator.calculate(list(ASPECTSCalculator.REGIONS))
        assert result.score == 0
        assert result.severity_category == "Unfavorable for intervention"


class TestHoehnYahrScale:
    """Tests for Hoehn & Yahr scale boundaries."""

    def test_hoehn_yahr_exists_in_knowledge(self):
        assert "hoehn_yahr" in CLINICAL_SCALES

    def test_hoehn_yahr_min_score(self):
        result = HoehnYahrCalculator.calculate(1.0)
        assert result.score == 1.0
        assert "Stage 1" in result.severity_category

    def test_hoehn_yahr_max_score(self):
        result = HoehnYahrCalculator.calculate(5.0)
        assert result.score == 5.0
        assert result.max_score == 5.0
        assert "Stage 5" in result.severity_category


class TestSharedDispatch:
    """`calculate_scale` is the single implementation both call paths resolve to."""

    def test_the_import_that_was_dead_now_resolves(self):
        assert callable(calculate_scale)

    def test_every_knowledge_base_scale_has_a_calculator(self):
        from src.clinical_scales import SUPPORTED_SCALES

        for name in CLINICAL_SCALES:
            assert name in SUPPORTED_SCALES, f"{name} is documented but not calculable"

    def test_an_unknown_scale_name_raises(self):
        with pytest.raises(KeyError):
            calculate_scale("not_a_scale", {})
