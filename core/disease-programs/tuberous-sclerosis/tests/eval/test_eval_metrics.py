"""
Validation-harness metrics (src/eval) against the default-seed synthetic cohort. These
assert MEASURED performance, the difference between "looks right" and "is right" — and
guard the headline claims (mosaic recovery, detection uplift, no benign call on a null,
full provenance) from regression. Construct validity only; clinical validation is Phase-1.
"""
from fastapi.testclient import TestClient

from api.main import app
from app._engine import featured, get_engine
from src.eval import evaluate


def _report():
    orch, manifest = get_engine()
    return evaluate(orch, manifest, featured())


def test_classification_is_faithful_and_safe():
    r = _report()
    c = r["variant_classification"]
    assert c["accuracy"] >= 0.95                       # faithful ACMG implementation
    assert c["false_benign_on_null_variant"] == 0      # a truncating variant is never benign


def test_detection_uplift_recovers_subthreshold_mosaics():
    dy = _report()["diagnostic_yield"]
    # standard germline-focused calling misses low-VAF mosaics; the engine recovers them
    assert dy["detection"]["engine_rate"] > dy["detection"]["standard_rate"]
    assert dy["detection"]["uplift_points"] >= 0.10
    assert dy["subthreshold_mosaic_recovery_sensitivity"] == 1.0
    assert dy["detection"]["false_detections_in_true_nmi"] == []   # no phantom variant in true-NMI


def test_molecular_diagnosis_uplift_all_true_positive():
    mdx = _report()["diagnostic_yield"]["molecular_diagnosis"]
    assert mdx["engine_rate"] > mdx["standard_rate"]
    assert mdx["newly_diagnosed_all_true_positive"] is True
    assert len(mdx["newly_diagnosed"]) >= 1


def test_forecast_calibration_probe_covers():
    mc = _report()["forecast_calibration"]
    assert mc["n_held_out_points"] >= 2
    assert mc["empirical_coverage_pi90"] == 1.0        # every held-out scan lands in the 90% band


def test_provenance_is_complete():
    p = _report()["provenance"]
    assert p["completeness"] == 1.0 and p["model_backed_outputs"] > 0


def test_eval_endpoint():
    with TestClient(app) as c:
        e = c.get("/eval").json()
    assert e["headline"]["mosaic_recovery_sensitivity"] == 1.0
    assert "construct validity" in e["cohort"].lower()
    assert e["caveat"]
