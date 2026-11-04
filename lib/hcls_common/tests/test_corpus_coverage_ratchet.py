"""The corpus-coverage ratchet must fail for the reasons a ratchet exists.

v2.2.0 made the agents report retrieval coverage honestly; 56 of 70 workflows retrieved from
under three-quarters of the evidence they declared. Seeding fixes that once. A ratchet is what
stops it drifting back, and the test-depth ratchet in this repo is the cautionary example: it
crashed on every invocation with `KeyError: 'orchestrator'` for long enough that nobody noticed
it had never run, and its badge-guard sibling could not fail at all because a regex never
matched and the no-match branch returned PASS.

So these tests are mostly about the ways a gate can quietly stop gating.
"""
import pytest

from hcls_common.corpus_coverage import (
    RATCHET_TOLERANCE,
    CoverageRegression,
    compare_coverage,
)


class TestItCatchesWhatItIsFor:
    def test_a_collection_going_quiet_is_a_regression(self):
        regressions, vanished = compare_coverage(
            {"cardiology/heart_failure": 0.82}, {"cardiology/heart_failure": 0.41}
        )
        assert len(regressions) == 1
        assert regressions[0].subject == "cardiology"
        assert regressions[0].workflow == "heart_failure"
        assert regressions[0].dropped == pytest.approx(0.41)
        assert not vanished

    def test_worst_regression_is_reported_first(self):
        regressions, _ = compare_coverage(
            {"a/small": 0.90, "a/large": 0.90, "a/mid": 0.90},
            {"a/small": 0.85, "a/large": 0.10, "a/mid": 0.50},
        )
        assert [r.workflow for r in regressions] == ["large", "mid", "small"]

    def test_a_workflow_that_stopped_being_measured_is_not_a_pass(self):
        """A vanished workflow and a healthy one are indistinguishable to a naive diff.

        This is the test-depth ratchet's exact failure: a subject that stops reporting looks
        like a subject that never regressed.
        """
        regressions, vanished = compare_coverage(
            {"neurology/stroke": 0.77, "neurology/epilepsy": 0.80},
            {"neurology/epilepsy": 0.80},
        )
        assert vanished == ["neurology/stroke"]
        assert not regressions, "a vanished workflow is not a coverage regression, it is worse"


class TestItDoesNotFireOnNoise:
    def test_improvement_is_never_a_regression(self):
        regressions, _ = compare_coverage({"a/x": 0.50}, {"a/x": 0.95})
        assert not regressions

    def test_a_move_inside_the_tolerance_is_ignored(self):
        """Renormalized weights move slightly when a collection is created or a boost changes."""
        regressions, _ = compare_coverage({"a/x": 0.800}, {"a/x": 0.795})
        assert not regressions

    def test_a_move_just_past_the_tolerance_is_caught(self):
        regressions, _ = compare_coverage(
            {"a/x": 0.80}, {"a/x": 0.80 - (RATCHET_TOLERANCE + 0.001)}
        )
        assert len(regressions) == 1

    def test_a_new_workflow_is_not_a_regression(self):
        """A workflow absent from the baseline has nothing to regress against."""
        regressions, vanished = compare_coverage({"a/x": 0.80}, {"a/x": 0.80, "a/new": 0.05})
        assert not regressions and not vanished


class TestKeyHandling:
    def test_tuple_and_string_keys_are_equivalent(self):
        """JSON round-trips force string keys; callers should not have to stringify twice."""
        as_tuples, _ = compare_coverage({("a", "x"): 0.80}, {("a", "x"): 0.20})
        as_strings, _ = compare_coverage({"a/x": 0.80}, {"a/x": 0.20})
        assert [str(r) for r in as_tuples] == [str(r) for r in as_strings]

    def test_a_workflow_name_containing_a_slash_still_parses_its_subject(self):
        regressions, _ = compare_coverage({("a", "x/y"): 0.80}, {("a", "x/y"): 0.10})
        assert regressions[0].subject == "a"
        assert regressions[0].workflow == "x/y"

    def test_empty_baseline_is_vacuously_clean(self):
        """A first run has no baseline. It must not be read as everything regressing."""
        regressions, vanished = compare_coverage({}, {"a/x": 0.10})
        assert not regressions and not vanished


class TestRenderingSaysSomethingUseful:
    def test_str_names_the_workflow_and_both_numbers(self):
        text = str(CoverageRegression("cardiology", "valvular_disease", 0.75, 0.02))
        assert "cardiology/valvular_disease" in text
        assert "75%" in text and "2%" in text
