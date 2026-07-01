"""A5 dual-model second opinion + A6 status state-machine / run search."""
import pytest
from hcls_common.governance import second_opinion, dual_predict
from hcls_common.mlops import MLOpsStore, RUN_STATES


class TestSecondOpinion:
    def test_agreement_is_high_confidence(self):
        r = second_opinion("toxic", "toxic", label="tox")
        assert r["agree"] and r["confidence"] == "high" and r["note"] is None

    def test_disagreement_is_low_confidence(self):
        r = second_opinion("toxic", "safe", label="tox")
        assert not r["agree"] and r["confidence"] == "low" and "disagree" in r["note"]

    def test_numeric_tolerance(self):
        assert second_opinion(0.81, 0.84, numeric_tol=0.05)["agree"]
        assert not second_opinion(0.20, 0.84, numeric_tol=0.05)["agree"]

    def test_dual_predict_extract(self):
        a = lambda x: {"verdict": "caution"}
        b = lambda x: {"verdict": "clean"}
        r = dual_predict(a, b, "CCO", extract=lambda d: d["verdict"], label="admet")
        assert not r["agree"] and r["raw_a"]["verdict"] == "caution"


class TestStatusMachine:
    def test_status_transitions_and_end_time(self):
        s = MLOpsStore(":memory:")
        rid = s.start_run("job", capability="esmfold-model")
        for st in ("submitted", "started", "running"):
            s.set_status(rid, st)
            assert s.get_run(rid)["status"] == st and s.get_run(rid)["ended_at"] is None
        s.set_status(rid, "complete")
        assert s.get_run(rid)["status"] == "complete" and s.get_run(rid)["ended_at"] is not None

    def test_invalid_status_rejected(self):
        s = MLOpsStore(":memory:")
        rid = s.start_run("j")
        with pytest.raises(ValueError, match="invalid status"):
            s.set_status(rid, "done")

    def test_search_runs(self):
        s = MLOpsStore(":memory:")
        a = s.start_run("fold-vcp", capability="esmfold-model"); s.set_status(a, "complete")
        b = s.start_run("fold-ras", capability="esmfold-model"); s.set_status(b, "failed")
        s.start_run("admet-x", capability="chemprop-admet")
        assert len(s.search_runs(name_like="fold")) == 2
        assert len(s.search_runs(capability="esmfold-model", status="failed")) == 1
        assert len(s.search_runs(status="complete")) == 1
