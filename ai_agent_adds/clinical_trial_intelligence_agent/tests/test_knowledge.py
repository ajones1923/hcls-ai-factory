"""Tests for knowledge data completeness.

Covers:
  - Landmark trials count and required fields
  - Regulatory milestones count and required fields
  - Therapeutic area coverage
  - Phase coverage

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.ingest.clinicaltrials_parser import LANDMARK_TRIALS
from src.ingest.regulatory_parser import REGULATORY_MILESTONES


class TestLandmarkTrials:
    """Test landmark trial knowledge data."""

    def test_count(self):
        assert len(LANDMARK_TRIALS) >= 20

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_nct_id(self, trial):
        assert "nct_id" in trial
        assert trial["nct_id"].startswith("NCT")

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_title(self, trial):
        assert "title" in trial
        assert len(trial["title"]) > 10

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_phase(self, trial):
        assert "phase" in trial
        assert "Phase" in trial["phase"]

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_status(self, trial):
        assert "status" in trial

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_conditions(self, trial):
        assert "conditions" in trial
        assert len(trial["conditions"]) >= 1

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_interventions(self, trial):
        assert "interventions" in trial
        assert len(trial["interventions"]) >= 1

    @pytest.mark.parametrize("trial", LANDMARK_TRIALS)
    def test_has_significance(self, trial):
        assert "significance" in trial
        assert len(trial["significance"]) > 10

    def test_unique_nct_ids(self):
        ids = [t["nct_id"] for t in LANDMARK_TRIALS]
        assert len(ids) == len(set(ids)), "Duplicate NCT IDs found"


class TestTherapeuticAreaCoverage:
    """Test that landmark trials cover multiple therapeutic areas."""

    def test_multiple_areas(self):
        areas = {t.get("therapeutic_area", "") for t in LANDMARK_TRIALS}
        assert len(areas) >= 3, f"Only {len(areas)} therapeutic areas covered"

    def test_oncology_present(self):
        areas = {t.get("therapeutic_area", "") for t in LANDMARK_TRIALS}
        assert "oncology" in areas

    def test_cardiology_present(self):
        areas = {t.get("therapeutic_area", "") for t in LANDMARK_TRIALS}
        assert "cardiology" in areas


class TestPhaseCoverage:
    """Test that landmark trials cover multiple phases."""

    def test_phase_iii_present(self):
        phases = {t.get("phase", "") for t in LANDMARK_TRIALS}
        assert "Phase III" in phases

    def test_multiple_phases(self):
        phases = {t.get("phase", "") for t in LANDMARK_TRIALS}
        assert len(phases) >= 2


class TestRegulatoryMilestones:
    """Test regulatory milestone knowledge data."""

    def test_count(self):
        assert len(REGULATORY_MILESTONES) >= 10

    @pytest.mark.parametrize("milestone", REGULATORY_MILESTONES)
    def test_has_drug(self, milestone):
        assert "drug" in milestone
        assert len(milestone["drug"]) > 0

    @pytest.mark.parametrize("milestone", REGULATORY_MILESTONES)
    def test_has_agency(self, milestone):
        assert "agency" in milestone
        assert milestone["agency"] in ("FDA", "EMA", "PMDA", "MHRA", "TGA", "Health Canada")

    @pytest.mark.parametrize("milestone", REGULATORY_MILESTONES)
    def test_has_decision(self, milestone):
        assert "decision" in milestone
        assert len(milestone["decision"]) > 0

    @pytest.mark.parametrize("milestone", REGULATORY_MILESTONES)
    def test_has_date(self, milestone):
        assert "date" in milestone
        # Date should be in YYYY-MM-DD format
        parts = milestone["date"].split("-")
        assert len(parts) == 3

    @pytest.mark.parametrize("milestone", REGULATORY_MILESTONES)
    def test_has_significance(self, milestone):
        assert "significance" in milestone
        assert len(milestone["significance"]) > 0
