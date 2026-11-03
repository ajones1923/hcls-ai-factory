"""A cardiac workflow that did not run must not report that it completed.

`_workflow_stub` backs all eight /workflow/* routes. It returned
`{"status": "completed", "recommendations": [], "evidence": []}` on both fallback paths — an
unknown workflow id, and a WorkflowEngine exception that was logged at warning level and then
discarded. A caller checking `status` saw a completed cardiac assessment containing nothing.

`/cross-modal/evaluate` had the same shape via `except Exception: pass`, reporting
`"status": "evaluated"` when the integration had raised.
"""
from unittest.mock import patch

import pytest
from fastapi import FastAPI, Request
from fastapi.testclient import TestClient

from api.routes.cardio_clinical import _workflow_stub, router


class _Lock:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture
def req():
    """A Request whose app.state carries just what _workflow_stub touches."""
    app = FastAPI()
    app.state.metrics = {"workflow_requests_total": 0}
    app.state.metrics_lock = _Lock()
    app.state.engine = None
    scope = {"type": "http", "app": app, "headers": []}
    return Request(scope)


def test_engine_failure_is_reported_as_failed(req):
    with patch("api.routes.cardio_clinical.WorkflowEngine") as engine:
        engine.return_value.run_workflow.side_effect = RuntimeError("milvus down")
        out = _workflow_stub("cad", "CAD Assessment", {"age": 61}, req)
    assert out["status"] == "failed"
    assert out["status"] != "completed"
    assert "milvus down" in out["error"]
    assert out["recommendations"] == []


def test_unknown_workflow_id_is_not_completed(req):
    out = _workflow_stub("not_a_workflow", "Nonsense", {}, req)
    assert out["status"] == "unknown_workflow"
    assert out["error"] is None


def test_successful_workflow_passes_through(req):
    class _Result:
        def model_dump(self):
            return {"status": "completed", "recommendations": ["statin"]}

    with patch("api.routes.cardio_clinical.WorkflowEngine") as engine:
        engine.return_value.run_workflow.return_value = _Result()
        out = _workflow_stub("cad", "CAD Assessment", {"age": 61}, req)
    assert out["status"] == "completed"
    assert out["recommendations"] == ["statin"]


class TestCrossModal:
    @pytest.fixture
    def client(self):
        app = FastAPI()
        app.include_router(router)
        app.state.metrics = {"workflow_requests_total": 0, "requests_total": 0}
        app.state.metrics_lock = _Lock()
        app.state.engine = None
        return TestClient(app)

    def test_failed_integration_is_not_reported_as_evaluated(self, client):
        class _Engine:
            def cross_modal_evaluate(self, _findings):
                raise RuntimeError("no cardio collections")

        client.app.state.engine = _Engine()
        data = client.post("/v1/cardio/cross-modal/evaluate", json={"echo": {}, "ct": {}}).json()
        assert data["status"] == "failed"
        assert "no cardio collections" in data["error"]

    def test_absent_engine_says_listed_only(self, client):
        data = client.post("/v1/cardio/cross-modal/evaluate", json={"echo": {}}).json()
        assert data["status"] == "listed_only"
        assert data["error"] is None
