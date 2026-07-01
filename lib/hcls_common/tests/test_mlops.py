"""Tests for single-box MLOps (A3) — in-memory SQLite, no external service."""
import pytest

from hcls_common.mlops import MLOpsStore


def store():
    return MLOpsStore(":memory:")


class TestRuns:
    def test_run_lifecycle_round_trip(self):
        s = store()
        rid = s.start_run("esmfold-trpcage", capability="esmfold-model",
                          params={"sequence_len": 20}, manifest={"model": "esmfold_v1"})
        s.log_metric(rid, "mean_plddt", 0.83)
        s.log_param(rid, "device", "cpu")
        s.end_run(rid)
        run = s.get_run(rid)
        assert run["status"] == "completed" and run["ended_at"] is not None
        assert run["params"]["sequence_len"] == 20 and run["params"]["device"] == "cpu"
        assert run["metrics"]["mean_plddt"] == 0.83
        assert run["manifest"]["model"] == "esmfold_v1"

    def test_list_runs_filtered_by_capability(self):
        s = store()
        s.end_run(s.start_run("a", capability="esmfold-model"))
        s.end_run(s.start_run("b", capability="chemprop-admet"))
        s.end_run(s.start_run("c", capability="esmfold-model"))
        assert len(s.list_runs(capability="esmfold-model")) == 2
        assert len(s.list_runs()) == 3

    def test_best_run_maximize_and_minimize(self):
        s = store()
        for v in (0.5, 0.9, 0.7):
            rid = s.start_run("r", capability="esmfold-model")
            s.log_metric(rid, "mean_plddt", v); s.end_run(rid)
        assert store_best(s, maximize=True) == 0.9
        assert store_best(s, maximize=False) == 0.5


def store_best(s, maximize):
    return s.best_run("mean_plddt", capability="esmfold-model", maximize=maximize)["metrics"]["mean_plddt"]


class TestModelRegistry:
    def test_register_and_get_version(self):
        s = store()
        s.register_model("esmfold", "1.0", source="hf:facebook/esmfold_v1",
                         metrics={"mean_plddt": 0.83})
        v = s.get_model_version("esmfold", "1.0")
        assert v["source"].endswith("esmfold_v1") and v["metrics"]["mean_plddt"] == 0.83
        assert v["stage"] == "none"

    def test_production_transition_archives_incumbent(self):
        s = store()
        s.register_model("admet", "1.0", stage="production")
        s.register_model("admet", "2.0")
        s.transition_stage("admet", "2.0", "production")
        prod = s.get_production("admet")
        assert prod["version"] == "2.0"
        assert s.get_model_version("admet", "1.0")["stage"] == "archived"  # incumbent archived

    def test_invalid_stage_rejected(self):
        s = store()
        with pytest.raises(ValueError, match="invalid stage"):
            s.register_model("x", "1.0", stage="prod")          # typo, not a valid stage
