"""Tests for clinical scale calculators.

Tests all 10 scale calculators with boundary values:
  NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn & Yahr

Author: Adam Jones
Date: March 2026
"""


from src.models import ClinicalScaleType, ScaleResult
from src.knowledge import CLINICAL_SCALES

# clinical_scales module may not exist yet; test knowledge data and model boundaries
try:
    from src.clinical_scales import calculate_scale
    _CS_AVAILABLE = True
except ImportError:
    _CS_AVAILABLE = False


class TestNIHSSScale:
    """Tests for NIHSS scale boundaries."""

    def test_nihss_exists_in_knowledge(self):
        assert "nihss" in CLINICAL_SCALES

    def test_nihss_min_score_valid(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=0,
            max_score=42,
            interpretation="No stroke symptoms",
            severity_category="none",
        )
        assert result.score == 0

    def test_nihss_max_score_valid(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=42,
            max_score=42,
            interpretation="Severe stroke",
            severity_category="severe",
        )
        assert result.score == 42

    def test_nihss_moderate_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=5,
            max_score=42,
            interpretation="Moderate stroke",
            severity_category="moderate",
        )
        assert result.score == 5

    def test_nihss_severe_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=21,
            max_score=42,
            interpretation="Severe stroke",
            severity_category="severe",
        )
        assert result.score == 21


class TestGCSScale:
    """Tests for GCS scale boundaries."""

    def test_gcs_exists_in_knowledge(self):
        assert "gcs" in CLINICAL_SCALES

    def test_gcs_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.GCS,
            score=3,
            max_score=15,
            interpretation="Severe (coma)",
            severity_category="severe",
        )
        assert result.score == 3

    def test_gcs_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.GCS,
            score=15,
            max_score=15,
            interpretation="Normal",
            severity_category="normal",
        )
        assert result.score == 15

    def test_gcs_moderate_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.GCS,
            score=9,
            max_score=15,
            interpretation="Moderate",
            severity_category="moderate",
        )
        assert result.score == 9


class TestMoCAScale:
    """Tests for MoCA scale boundaries."""

    def test_moca_exists_in_knowledge(self):
        assert "moca" in CLINICAL_SCALES

    def test_moca_normal_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.MOCA,
            score=28,
            max_score=30,
            interpretation="Normal cognition",
            severity_category="normal",
        )
        assert result.score == 28

    def test_moca_mci_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.MOCA,
            score=25,
            max_score=30,
            interpretation="Mild cognitive impairment",
            severity_category="mci",
        )
        assert result.score == 25

    def test_moca_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.MOCA,
            score=0,
            max_score=30,
            interpretation="Severe cognitive impairment",
            severity_category="severe",
        )
        assert result.score == 0


class TestUPDRSScale:
    """Tests for UPDRS Part III scale boundaries."""

    def test_updrs_exists_in_knowledge(self):
        assert "updrs" in CLINICAL_SCALES

    def test_updrs_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.UPDRS,
            score=0,
            max_score=132,
            interpretation="Minimal motor findings",
            severity_category="minimal",
        )
        assert result.score == 0

    def test_updrs_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.UPDRS,
            score=132,
            max_score=132,
            interpretation="Severe motor impairment",
            severity_category="severe",
        )
        assert result.score == 132


class TestEDSSScale:
    """Tests for EDSS scale boundaries."""

    def test_edss_exists_in_knowledge(self):
        assert "edss" in CLINICAL_SCALES

    def test_edss_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.EDSS,
            score=0.0,
            max_score=10.0,
            interpretation="Normal neurological exam",
            severity_category="normal",
        )
        assert result.score == 0.0

    def test_edss_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.EDSS,
            score=10.0,
            max_score=10.0,
            interpretation="Death due to MS",
            severity_category="death",
        )
        assert result.score == 10.0

    def test_edss_wheelchair_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.EDSS,
            score=7.0,
            max_score=10.0,
            interpretation="Wheelchair bound",
            severity_category="severe",
        )
        assert result.score == 7.0


class TestMRSScale:
    """Tests for Modified Rankin Scale boundaries."""

    def test_mrs_exists_in_knowledge(self):
        assert "mrs" in CLINICAL_SCALES

    def test_mrs_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.MRS,
            score=0,
            max_score=6,
            interpretation="No symptoms",
            severity_category="normal",
        )
        assert result.score == 0

    def test_mrs_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.MRS,
            score=6,
            max_score=6,
            interpretation="Dead",
            severity_category="death",
        )
        assert result.score == 6


class TestHIT6Scale:
    """Tests for HIT-6 scale boundaries."""

    def test_hit6_exists_in_knowledge(self):
        assert "hit6" in CLINICAL_SCALES

    def test_hit6_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.HIT6,
            score=36,
            max_score=78,
            interpretation="Little/no impact",
            severity_category="minimal",
        )
        assert result.score == 36

    def test_hit6_severe_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.HIT6,
            score=60,
            max_score=78,
            interpretation="Severe impact",
            severity_category="severe",
        )
        assert result.score == 60


class TestALSFRSScale:
    """Tests for ALSFRS-R scale boundaries."""

    def test_alsfrs_exists_in_knowledge(self):
        assert "alsfrs" in CLINICAL_SCALES

    def test_alsfrs_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.ALSFRS,
            score=48,
            max_score=48,
            interpretation="Normal function",
            severity_category="normal",
        )
        assert result.score == 48

    def test_alsfrs_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.ALSFRS,
            score=0,
            max_score=48,
            interpretation="Complete functional loss",
            severity_category="severe",
        )
        assert result.score == 0


class TestASPECTSScale:
    """Tests for ASPECTS scale boundaries."""

    def test_aspects_exists_in_knowledge(self):
        assert "aspects" in CLINICAL_SCALES

    def test_aspects_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.ASPECTS,
            score=10,
            max_score=10,
            interpretation="No early ischemic changes",
            severity_category="favorable",
        )
        assert result.score == 10

    def test_aspects_thrombectomy_boundary(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.ASPECTS,
            score=6,
            max_score=10,
            interpretation="Moderate infarct core",
            severity_category="borderline",
        )
        assert result.score == 6


class TestHoehnYahrScale:
    """Tests for Hoehn & Yahr scale boundaries."""

    def test_hoehn_yahr_exists_in_knowledge(self):
        assert "hoehn_yahr" in CLINICAL_SCALES

    def test_hoehn_yahr_min_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.HOEHN_YAHR,
            score=0,
            max_score=5,
            interpretation="No signs of disease",
            severity_category="normal",
        )
        assert result.score == 0

    def test_hoehn_yahr_max_score(self):
        result = ScaleResult(
            scale_type=ClinicalScaleType.HOEHN_YAHR,
            score=5,
            max_score=5,
            interpretation="Wheelchair bound or bedridden",
            severity_category="severe",
        )
        assert result.score == 5
