"""Retrieval coverage: the arithmetic, and the case it was written for.

The measured case is clinical-trial's PATIENT_MATCHING workflow on 2026-11-03 — ten of fourteen
collections empty, carrying 73% of the declared weight, with `collections_searched` reporting 14.
It is pinned here as a reference case so the number cannot drift unnoticed.
"""
import pytest

from hcls_common.corpus_coverage import DEGRADED_BELOW, coverage

# clinical-trial src/vector_collections.py, WORKFLOW_COLLECTION_WEIGHTS[PATIENT_MATCHING] -- the
# declared budget. The engine multiplies a second table (WORKFLOW_COLLECTION_BOOST) against base
# weights instead and lands at 28.8% for the same corpus; the agent's own test follows that path.
# Either way ~73% of the weight points at nothing, which is the finding.
PATIENT_MATCHING = {
    "trial_eligibility": 0.25,
    "trial_protocols": 0.15,
    "trial_sites": 0.12,
    "trial_biomarkers": 0.10,
    "genomic_evidence": 0.08,
    "trial_endpoints": 0.06,
    "trial_safety": 0.05,
    "trial_results": 0.05,
    "trial_guidelines": 0.04,
    "trial_literature": 0.03,
    "trial_investigators": 0.03,
    "trial_rwe": 0.02,
    "trial_regulatory": 0.01,
    "trial_adaptive": 0.01,
}

# Measured against the live instance the same day: the four that held anything.
HELD_SOMETHING = {
    "trial_protocols": 156,
    "trial_literature": 150,
    "trial_regulatory": 22,
    "genomic_evidence": 60,
}


class TestTheMeasuredCase:
    def test_patient_matching_covered_27_percent_of_its_declared_weight(self):
        v = coverage(PATIENT_MATCHING, HELD_SOMETHING)
        assert v.coverage == pytest.approx(0.27, abs=0.005)

    def test_ten_of_fourteen_sources_were_silent(self):
        v = coverage(PATIENT_MATCHING, HELD_SOMETHING)
        assert len(v.declared) == 14
        assert len(v.contributing) == 4
        assert len(v.silent) == 10

    def test_it_is_degraded(self):
        assert coverage(PATIENT_MATCHING, HELD_SOMETHING).is_degraded

    def test_the_heaviest_silent_source_is_named_first(self):
        """A truncated disclosure must name eligibility, not `trial_adaptive` at 0.01."""
        v = coverage(PATIENT_MATCHING, HELD_SOMETHING)
        assert v.silent[0] == "trial_eligibility"
        assert v.silent[1] == "trial_sites"

    def test_the_disclosure_states_the_share_and_names_the_worst(self):
        d = coverage(PATIENT_MATCHING, HELD_SOMETHING).disclosure()
        assert "10 of 14" in d
        assert "73%" in d          # the weight that returned nothing
        assert "27%" in d          # what the answer actually rests on
        assert "trial_eligibility" in d

    def test_the_disclosure_truncates_but_says_how_many_it_dropped(self):
        d = coverage(PATIENT_MATCHING, HELD_SOMETHING).disclosure(max_named=4)
        assert "and 6 more" in d


class TestCompleteCoverageSaysNothing:
    def test_a_fully_populated_corpus_produces_no_disclosure(self):
        v = coverage(PATIENT_MATCHING, {k: 5 for k in PATIENT_MATCHING})
        assert v.coverage == pytest.approx(1.0)
        assert v.is_complete
        assert not v.is_degraded
        assert v.disclosure() == ""

    def test_one_trivial_source_missing_is_not_degraded_but_is_disclosed(self):
        """0.01 of 1.00 gone: worth recording, not worth a warning."""
        counts = {k: 5 for k in PATIENT_MATCHING if k != "trial_adaptive"}
        v = coverage(PATIENT_MATCHING, counts)
        assert not v.is_degraded
        assert not v.is_complete
        assert "trial_adaptive" in v.disclosure()

    def test_losing_only_eligibility_crosses_the_threshold(self):
        """0.25 is a quarter of the budget in one collection."""
        counts = {k: 5 for k in PATIENT_MATCHING if k != "trial_eligibility"}
        v = coverage(PATIENT_MATCHING, counts)
        assert v.coverage == pytest.approx(0.75, abs=0.005)
        assert v.coverage >= DEGRADED_BELOW  # exactly at the floor is not degraded


class TestNothingRetrievedAtAll:
    def test_an_entirely_empty_corpus_reports_zero_not_an_error(self):
        v = coverage(PATIENT_MATCHING, {})
        assert v.coverage == 0.0
        assert v.is_degraded
        assert len(v.silent) == 14

    def test_an_empty_weight_map_does_not_divide_by_zero(self):
        v = coverage({}, {})
        assert v.coverage == 0.0
        assert v.is_complete          # nothing was declared, so nothing went missing
        assert v.disclosure() == ""


class TestTheThreeResultShapesCallersActuallyHave:
    def test_objects_with_a_collection_attribute(self):
        class R:
            def __init__(self, c):
                self.collection = c

        v = coverage({"a": 0.5, "b": 0.5}, [R("a"), R("a"), R("a")])
        assert v.contributing == ("a",)
        assert v.coverage == pytest.approx(0.5)

    def test_dicts_carrying_a_collection_key(self):
        v = coverage({"a": 0.5, "b": 0.5}, [{"collection": "b", "score": 0.9}])
        assert v.contributing == ("b",)

    def test_name_count_pairs(self):
        v = coverage({"a": 0.5, "b": 0.5}, [("a", 2), ("b", 0)])
        assert v.contributing == ("a",)
        assert v.silent == ("b",)

    def test_a_result_naming_an_undeclared_collection_does_not_inflate_coverage(self):
        """Retrieval reaching outside the weight map is not coverage of the map."""
        v = coverage({"a": 1.0}, [{"collection": "somewhere_else"}])
        assert v.coverage == 0.0
        assert v.silent == ("a",)


class TestWeightSemantics:
    def test_a_zero_weight_collection_is_not_declared(self):
        """A weight map saying 0.0 is saying "not for this workflow"."""
        v = coverage({"a": 1.0, "b": 0.0}, {"a": 3})
        assert v.declared == ("a",)
        assert v.is_complete
        assert v.coverage == pytest.approx(1.0)

    def test_coverage_is_weighted_not_counted(self):
        """Three empty trivial sources beat one empty heavy one on count, and lose on weight."""
        heavy_gone = coverage({"big": 0.7, "s1": 0.1, "s2": 0.1, "s3": 0.1},
                              {"s1": 1, "s2": 1, "s3": 1})
        light_gone = coverage({"big": 0.7, "s1": 0.1, "s2": 0.1, "s3": 0.1}, {"big": 1})
        assert len(heavy_gone.silent) < len(light_gone.silent)
        assert heavy_gone.coverage < light_gone.coverage

    def test_weights_need_not_sum_to_one(self):
        """Coverage is a share of what was declared, whatever that totals."""
        v = coverage({"a": 2.0, "b": 2.0}, {"a": 1})
        assert v.weight_declared == pytest.approx(4.0)
        assert v.coverage == pytest.approx(0.5)


class TestTheSerializedForm:
    def test_as_dict_carries_the_numbers_a_response_needs(self):
        d = coverage(PATIENT_MATCHING, HELD_SOMETHING).as_dict()
        assert d["degraded"] is True
        assert d["collections_declared"] == 14
        assert d["collections_contributing"] == 4
        assert d["coverage"] == pytest.approx(0.27, abs=0.005)
        assert d["silent_collections"][0] == "trial_eligibility"

    def test_as_dict_is_json_serializable(self):
        import json
        json.dumps(coverage(PATIENT_MATCHING, HELD_SOMETHING).as_dict())


class TestRenderingFromTheSerializedForm:
    """`disclosure_from_dict` exists because cart keeps the dict and renders a layer later."""

    def test_it_produces_the_same_paragraph_as_the_verdict(self):
        from hcls_common.corpus_coverage import disclosure_from_dict
        v = coverage(PATIENT_MATCHING, HELD_SOMETHING)
        assert disclosure_from_dict(v.as_dict()) == v.disclosure()

    def test_a_complete_verdict_renders_nothing(self):
        from hcls_common.corpus_coverage import disclosure_from_dict
        v = coverage(PATIENT_MATCHING, {k: 1 for k in PATIENT_MATCHING})
        assert disclosure_from_dict(v.as_dict()) == ""

    def test_an_empty_dict_does_not_raise(self):
        """A caller that never ran retrieval has no verdict, not a broken one."""
        from hcls_common.corpus_coverage import disclosure_from_dict
        assert disclosure_from_dict({}) == ""

    def test_it_truncates_the_same_way(self):
        from hcls_common.corpus_coverage import disclosure_from_dict
        v = coverage(PATIENT_MATCHING, HELD_SOMETHING)
        assert "and 6 more" in disclosure_from_dict(v.as_dict(), max_named=4)
