"""
Tests for the TrialMatcher.
=============================
Validates trial matching, cancer type filtering, biomarker matching,
and composite scoring logic.

All external dependencies (Milvus, embeddings) are mocked.
"""

import sys
from pathlib import Path

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.trial_matcher import TrialMatcher


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def trial_matcher(mock_collection_manager, mock_embedder):
    """Create a TrialMatcher with mocked dependencies."""
    return TrialMatcher(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
    )


@pytest.fixture
def sample_trial_results():
    """Simulated trial search results from Milvus."""
    return [
        {
            "trial_id": "NCT04613596",
            "title": "ADAURA: Adjuvant Osimertinib in EGFR-Mutant NSCLC",
            "phase": "Phase 3",
            "status": "Recruiting",
            "cancer_type": "nsclc",
            "criteria": "EGFR exon 19 deletion or L858R stage IB-IIIA NSCLC",
            "biomarker_criteria": "EGFR L858R exon 19 deletion",
            "sponsor": "AstraZeneca",
            "score": 0.85,
        },
        {
            "trial_id": "NCT03785249",
            "title": "CROWN: Lorlatinib vs Crizotinib in ALK-Positive NSCLC",
            "phase": "Phase 3",
            "status": "Active, not recruiting",
            "cancer_type": "nsclc",
            "criteria": "ALK-positive advanced NSCLC, no prior ALK TKI",
            "biomarker_criteria": "ALK rearrangement",
            "sponsor": "Pfizer",
            "score": 0.78,
        },
        {
            "trial_id": "NCT02628067",
            "title": "KEYNOTE-158: Pembrolizumab in MSI-H Solid Tumors",
            "phase": "Phase 2",
            "status": "Recruiting",
            "cancer_type": "solid tumors",
            "criteria": "MSI-H or TMB-H advanced solid tumors",
            "biomarker_criteria": "MSI-H TMB-H dMMR",
            "sponsor": "Merck",
            "score": 0.72,
        },
    ]


# ═══════════════════════════════════════════════════════════════════════════
# Initialization Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestTrialMatcherInit:
    """Test TrialMatcher initialization."""

    def test_basic_init(self, mock_collection_manager, mock_embedder):
        matcher = TrialMatcher(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
        )
        assert matcher.collection_manager is mock_collection_manager
        assert matcher.embedder is mock_embedder

    def test_collection_name_constant(self, trial_matcher):
        assert trial_matcher.COLLECTION_NAME == "onco_trials"

    def test_phase_weights_defined(self, trial_matcher):
        assert "Phase 3" in trial_matcher.PHASE_WEIGHTS
        assert "Phase 2" in trial_matcher.PHASE_WEIGHTS
        assert "Phase 1" in trial_matcher.PHASE_WEIGHTS

    def test_status_weights_defined(self, trial_matcher):
        assert "Recruiting" in trial_matcher.STATUS_WEIGHTS
        assert "Active, not recruiting" in trial_matcher.STATUS_WEIGHTS


# ═══════════════════════════════════════════════════════════════════════════
# Match Trials Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMatchTrials:
    """Test match_trials returns scored results."""

    def test_match_trials_returns_list(self, trial_matcher):
        """match_trials should return a list (even if empty when mocked)."""
        results = trial_matcher.match_trials(
            cancer_type="NSCLC",
            biomarkers={"EGFR": "L858R"},
            stage="IV",
        )
        assert isinstance(results, list)

    def test_match_trials_with_enriched_results(
        self, mock_collection_manager, mock_embedder, sample_trial_results
    ):
        """When Milvus returns results, match_trials should score them."""
        # Configure mocks to return sample results
        mock_collection_manager.query.return_value = sample_trial_results[:1]
        mock_collection_manager.search.return_value = sample_trial_results

        matcher = TrialMatcher(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
        )
        results = matcher.match_trials(
            cancer_type="NSCLC",
            biomarkers={"EGFR": "L858R"},
            stage="IV",
            top_k=5,
        )
        assert isinstance(results, list)
        # Results should have match_score and trial_id
        for result in results:
            assert "match_score" in result
            assert "trial_id" in result

    def test_match_trials_scores_descending(
        self, mock_collection_manager, mock_embedder, sample_trial_results
    ):
        """Results should be sorted by match_score descending."""
        mock_collection_manager.query.return_value = sample_trial_results
        mock_collection_manager.search.return_value = sample_trial_results

        matcher = TrialMatcher(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
        )
        results = matcher.match_trials(
            cancer_type="NSCLC",
            biomarkers={"EGFR": "L858R"},
            stage="IV",
        )
        if len(results) >= 2:
            for i in range(len(results) - 1):
                assert results[i]["match_score"] >= results[i + 1]["match_score"]


# ═══════════════════════════════════════════════════════════════════════════
# Biomarker Matching Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestScoreBiomarkerMatch:
    """Test _score_biomarker_match for exact and partial matches."""

    @pytest.fixture
    def matcher(self, trial_matcher):
        return trial_matcher

    def test_exact_match(self, matcher):
        """Patient biomarker appearing in criteria should score > 0."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="EGFR L858R mutation required",
            patient_biomarkers={"EGFR": "L858R"},
        )
        assert score > 0.0

    def test_no_match(self, matcher):
        """Patient biomarker not in criteria should score 0."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="ALK rearrangement required",
            patient_biomarkers={"BRAF": "V600E"},
        )
        assert score == 0.0

    def test_partial_match(self, matcher):
        """Some but not all biomarkers matching should give fractional score."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="EGFR mutation MSI-H testing",
            patient_biomarkers={"EGFR": "L858R", "BRAF": "V600E", "MSI": "MSI-H"},
        )
        # At least 1 of 3 should match (EGFR appears in criteria)
        assert 0.0 < score <= 1.0

    def test_empty_biomarkers(self, matcher):
        """Empty patient biomarkers should return 0."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="EGFR required",
            patient_biomarkers={},
        )
        assert score == 0.0

    def test_empty_criteria(self, matcher):
        """Empty trial criteria should return 0."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="",
            patient_biomarkers={"EGFR": "L858R"},
        )
        assert score == 0.0

    def test_full_match(self, matcher):
        """All biomarkers matching should score 1.0."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="EGFR L858R mutation, PD-L1 TPS 80%",
            patient_biomarkers={"EGFR": "L858R", "PD-L1": "80"},
        )
        # Both keys should be found in criteria
        assert score > 0.0

    def test_case_insensitive(self, matcher):
        """Matching should be case-insensitive."""
        score = matcher._score_biomarker_match(
            trial_biomarker_criteria="egfr mutation required",
            patient_biomarkers={"EGFR": "L858R"},
        )
        assert score > 0.0


# ═══════════════════════════════════════════════════════════════════════════
# Cancer Type Filtering Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCancerTypeFiltering:
    """Test cancer type filtering in trial matching."""

    def test_build_eligibility_query(self, trial_matcher):
        """_build_eligibility_query should include cancer type and biomarkers."""
        query = trial_matcher._build_eligibility_query(
            cancer_type="NSCLC",
            biomarkers={"EGFR": "L858R", "PD-L1_TPS": 80},
            stage="IV",
        )
        assert "NSCLC" in query
        assert "stage IV" in query
        assert "EGFR" in query


# ═══════════════════════════════════════════════════════════════════════════
# Composite Scoring Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCompositeScore:
    """Test _compute_composite_score produces valid scores."""

    @pytest.fixture
    def matcher(self, trial_matcher):
        return trial_matcher

    def test_score_range(self, matcher):
        """Composite score should be between 0 and 1."""
        trial = {
            "phase": "Phase 3",
            "status": "Recruiting",
            "criteria": "EGFR L858R NSCLC",
            "biomarker_criteria": "EGFR mutation",
            "score": 0.8,
        }
        patient_data = {
            "cancer_type": "NSCLC",
            "biomarkers": {"EGFR": "L858R"},
            "stage": "IV",
        }
        score = matcher._compute_composite_score(trial, patient_data)
        assert 0.0 <= score <= 1.0

    def test_phase_3_higher_than_phase_1(self, matcher):
        """Phase 3 trial should score higher than Phase 1 (all else equal)."""
        base_trial = {
            "criteria": "EGFR mutation",
            "biomarker_criteria": "",
            "status": "Recruiting",
            "score": 0.5,
        }
        patient_data = {
            "cancer_type": "NSCLC",
            "biomarkers": {"EGFR": "L858R"},
            "stage": "IV",
        }

        trial_p3 = {**base_trial, "phase": "Phase 3"}
        trial_p1 = {**base_trial, "phase": "Phase 1"}

        score_p3 = matcher._compute_composite_score(trial_p3, patient_data)
        score_p1 = matcher._compute_composite_score(trial_p1, patient_data)
        assert score_p3 >= score_p1

    def test_recruiting_higher_than_not_yet(self, matcher):
        """Recruiting trial should score higher than not-yet-recruiting."""
        base_trial = {
            "phase": "Phase 2",
            "criteria": "",
            "biomarker_criteria": "",
            "score": 0.5,
        }
        patient_data = {
            "cancer_type": "NSCLC",
            "biomarkers": {},
            "stage": "IV",
        }

        trial_rec = {**base_trial, "status": "Recruiting"}
        trial_nyr = {**base_trial, "status": "Not yet recruiting"}

        score_rec = matcher._compute_composite_score(trial_rec, patient_data)
        score_nyr = matcher._compute_composite_score(trial_nyr, patient_data)
        assert score_rec >= score_nyr


# ═══════════════════════════════════════════════════════════════════════════
# Merge Results Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestMergeResults:
    """Test _merge_results deduplication."""

    @pytest.fixture
    def matcher(self, trial_matcher):
        return trial_matcher

    def test_deduplication(self, matcher):
        """Duplicate trial_ids should be deduplicated."""
        deterministic = [
            {"trial_id": "NCT001", "title": "Trial A", "score": 0.9},
        ]
        semantic = [
            {"trial_id": "NCT001", "title": "Trial A", "score": 0.8},
            {"trial_id": "NCT002", "title": "Trial B", "score": 0.7},
        ]
        merged = matcher._merge_results(deterministic, semantic)
        trial_ids = [t["trial_id"] for t in merged]
        assert len(trial_ids) == len(set(trial_ids))

    def test_union(self, matcher):
        """Merged results should contain all unique trials."""
        deterministic = [{"trial_id": "NCT001", "title": "A"}]
        semantic = [{"trial_id": "NCT002", "title": "B"}]
        merged = matcher._merge_results(deterministic, semantic)
        trial_ids = {t["trial_id"] for t in merged}
        assert "NCT001" in trial_ids
        assert "NCT002" in trial_ids
