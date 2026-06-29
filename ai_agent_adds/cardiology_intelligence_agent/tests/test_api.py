"""FastAPI TestClient tests for Cardiology Intelligence Agent API.

Tests actual HTTP request/response cycles with mocked backends.
"""

import sys
import threading
from contextlib import asynccontextmanager
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ── Helpers ──────────────────────────────────────────────────────────────

@asynccontextmanager
async def _noop_lifespan(app):
    """No-op lifespan that skips Milvus/LLM/model initialization."""
    yield


def _make_client():
    """Build a TestClient with the real app but a no-op lifespan."""
    from fastapi.testclient import TestClient
    from api.main import app

    # Replace lifespan so TestClient startup does NOT try to connect
    # to Milvus, load sentence-transformers, or call Anthropic.
    app.router.lifespan_context = _noop_lifespan

    # Inject safe defaults into app.state
    app.state.engine = None
    app.state.manager = None
    app.state.risk_calculators = {}
    app.state.gdmt_optimizer = None
    app.state.metrics = {
        "requests_total": 0,
        "query_requests_total": 0,
        "search_requests_total": 0,
        "risk_calc_requests_total": 0,
        "workflow_requests_total": 0,
        "gdmt_requests_total": 0,
        "report_requests_total": 0,
        "errors_total": 0,
    }
    app.state.metrics_lock = threading.Lock()

    return TestClient(app, raise_server_exceptions=False)


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture
def client():
    """Create a TestClient with mocked lifespan (no Milvus/LLM required)."""
    from api.main import _rate_limit_store
    _rate_limit_store.clear()

    with _make_client() as c:
        yield c


@pytest.fixture
def client_with_risk_engine():
    """TestClient with the real RiskCalculatorEngine loaded into app.state."""
    from api.main import app, _rate_limit_store
    from src.risk_calculators import RiskCalculatorEngine

    _rate_limit_store.clear()

    tc = _make_client()
    calc_engine = RiskCalculatorEngine()
    app.state.risk_calculators = {
        "engine": calc_engine,
        "ascvd": calc_engine,
        "heart": calc_engine,
        "cha2ds2_vasc": calc_engine,
        "has_bled": calc_engine,
        "maggic": calc_engine,
        "euroscore": calc_engine,
    }

    with tc as c:
        yield c


@pytest.fixture
def client_with_api_key():
    """TestClient where API_KEY is set so auth middleware is active."""
    from api.main import _rate_limit_store
    from config.settings import settings

    _rate_limit_store.clear()

    original_key = settings.API_KEY
    settings.API_KEY = "test-secret-key-12345"
    try:
        with _make_client() as c:
            yield c
    finally:
        settings.API_KEY = original_key


# =====================================================================
# 1. GET /health -- degraded status when no Milvus
# =====================================================================

class TestHealth:
    def test_health_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_health_degraded_without_milvus(self, client):
        data = client.get("/health").json()
        assert data["status"] == "degraded"

    def test_health_has_all_component_fields(self, client):
        data = client.get("/health").json()
        assert "agent" in data
        assert "version" in data
        components = data["components"]
        assert "milvus" in components
        assert "rag_engine" in components
        assert "gdmt_optimizer" in components
        assert "risk_calculators" in components

    def test_health_components_show_unavailable(self, client):
        components = client.get("/health").json()["components"]
        assert components["milvus"] == "unavailable"
        assert components["rag_engine"] == "unavailable"


# =====================================================================
# 2. GET /collections -- 503 when Milvus is unavailable
# =====================================================================

class TestCollections:
    def test_collections_503_without_milvus(self, client):
        resp = client.get("/collections")
        assert resp.status_code == 503


# =====================================================================
# 3. GET /workflows -- returns all workflow definitions
# =====================================================================

class TestWorkflows:
    def test_workflows_returns_200(self, client):
        resp = client.get("/workflows")
        assert resp.status_code == 200

    def test_workflows_lists_all(self, client):
        data = client.get("/workflows").json()
        workflows = data["workflows"]
        assert len(workflows) == 11
        ids = {w["id"] for w in workflows}
        expected = {
            "cad", "heart_failure", "valvular", "arrhythmia",
            "cardiac_mri", "stress_test", "prevention", "cardio_oncology",
            "acute_decompensated_hf", "post_mi", "myocarditis_pericarditis",
        }
        assert ids == expected

    def test_workflows_have_required_fields(self, client):
        workflows = client.get("/workflows").json()["workflows"]
        for w in workflows:
            assert "id" in w
            assert "name" in w
            assert "description" in w
            assert "risk_calculators" in w


# =====================================================================
# 4. GET /metrics -- Prometheus-compatible text
# =====================================================================

class TestMetrics:
    def test_metrics_returns_200(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200

    def test_metrics_text_format(self, client):
        text = client.get("/metrics").text
        # May come from prometheus_client (src.metrics) or fallback counters
        assert "# TYPE" in text
        # Either Prometheus-native cardio metrics or the simple fallback
        assert "cardio_" in text or "cardio_agent_" in text


# =====================================================================
# 5-8. Risk calculator endpoints (no engine loaded)
# =====================================================================

class TestRiskCalculatorsNoEngine:
    """Risk endpoints when RiskCalculatorEngine is NOT loaded."""

    def test_ascvd_works_with_lazy_init(self, client):
        """ASCVD lazy-initializes RiskCalculatorEngine if not pre-loaded."""
        resp = client.post("/v1/cardio/risk/ascvd", json={
            "age": 55, "sex": "male", "race": "white",
            "total_cholesterol": 213, "hdl_cholesterol": 50,
            "systolic_bp": 120, "bp_treatment": False,
            "diabetes": False, "smoker": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["calculator"] == "ASCVD_Pooled_Cohort_Equations"
        assert data["score"] >= 0

    def test_heart_score_works_inline(self, client):
        """HEART score uses inline scoring -- works without engine."""
        resp = client.post("/v1/cardio/risk/heart-score", json={
            "history_score": 1, "ecg_score": 1,
            "age": 55, "risk_factors": 2, "troponin_score": 0,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["calculator"] == "HEART_Score"
        assert "risk_category" in data
        assert "recommendations" in data
        # total = 1 + 1 + 1(age 45-64) + 1(rf 1-2) + 0 = 4 => moderate
        assert data["score"] == 4.0
        assert data["risk_category"] == "moderate"

    def test_cha2ds2_vasc_works_inline(self, client):
        """CHA2DS2-VASc uses inline scoring -- works without engine."""
        resp = client.post("/v1/cardio/risk/cha2ds2-vasc", json={
            "chf": True, "hypertension": True,
            "age": 76, "diabetes": False,
            "stroke_tia": False, "vascular_disease": False,
            "sex": "male",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["calculator"] == "CHA2DS2-VASc"
        # CHF=1, HTN=1, age>=75=2 => 4
        assert data["score"] == 4.0
        assert data["risk_category"] == "high"

    def test_has_bled_works_inline(self, client):
        """HAS-BLED uses inline scoring -- works without engine."""
        resp = client.post("/v1/cardio/risk/has-bled", json={
            "hypertension_uncontrolled": True,
            "renal_disease": False, "liver_disease": False,
            "stroke_history": False, "bleeding_history": True,
            "labile_inr": False, "age_over_65": True,
            "drugs_alcohol": 0,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["calculator"] == "HAS-BLED"
        # HTN=1, bleeding=1, age=1 => 3 => high
        assert data["score"] == 3.0
        assert data["risk_category"] == "high"


# =====================================================================
# 9. ASCVD with real RiskCalculatorEngine loaded
# =====================================================================

class TestRiskCalculatorsWithEngine:
    def test_ascvd_with_engine(self, client_with_risk_engine):
        resp = client_with_risk_engine.post("/v1/cardio/risk/ascvd", json={
            "age": 55, "sex": "male", "race": "white",
            "total_cholesterol": 213, "hdl_cholesterol": 50,
            "systolic_bp": 120, "bp_treatment": False,
            "diabetes": False, "smoker": False,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["calculator"] == "ASCVD_Pooled_Cohort_Equations"
        assert isinstance(data["score"], float)
        assert data["score"] >= 0.0
        assert "risk_category" in data
        assert "interpretation" in data
        assert "recommendations" in data
        assert isinstance(data["recommendations"], list)

    def test_ascvd_high_risk_patient(self, client_with_risk_engine):
        """High-risk patient: older, diabetic, smoker, high BP."""
        resp = client_with_risk_engine.post("/v1/cardio/risk/ascvd", json={
            "age": 70, "sex": "male", "race": "african_american",
            "total_cholesterol": 280, "hdl_cholesterol": 35,
            "systolic_bp": 170, "bp_treatment": True,
            "diabetes": True, "smoker": True,
        })
        assert resp.status_code == 200
        data = resp.json()
        # Should be high risk
        assert data["score"] > 10.0


# =====================================================================
# 10-12. Authentication middleware tests
# =====================================================================

class TestAuth:
    def test_api_key_enables_auth(self, client_with_api_key):
        """Requests without a valid key should be rejected."""
        resp = client_with_api_key.get("/workflows")
        assert resp.status_code == 401

    def test_health_exempt_from_auth(self, client_with_api_key):
        """/health is in _AUTH_SKIP_PATHS and should always be accessible."""
        resp = client_with_api_key.get("/health")
        assert resp.status_code == 200

    def test_valid_key_passes_auth(self, client_with_api_key):
        """Providing the correct X-API-Key header should succeed."""
        resp = client_with_api_key.get(
            "/workflows",
            headers={"X-API-Key": "test-secret-key-12345"},
        )
        assert resp.status_code == 200

    def test_invalid_key_returns_401(self, client_with_api_key):
        resp = client_with_api_key.get(
            "/workflows",
            headers={"X-API-Key": "wrong-key"},
        )
        assert resp.status_code == 401
        assert "Invalid" in resp.json().get("detail", "")

    def test_metrics_exempt_from_auth(self, client_with_api_key):
        """/metrics is in _AUTH_SKIP_PATHS."""
        resp = client_with_api_key.get("/metrics")
        assert resp.status_code == 200


# =====================================================================
# 13. Rate limiting -- > 100 rapid requests triggers 429
# =====================================================================

class TestRateLimiting:
    def test_rate_limit_triggers_429(self, client):
        """Exceeding 100 requests in the rate-limit window returns 429."""
        # Send 101 requests to a lightweight endpoint
        statuses = []
        for _ in range(101):
            resp = client.get("/workflows")
            statuses.append(resp.status_code)

        # At least the 101st should be 429
        assert 429 in statuses, "Expected at least one 429 after 100+ requests"
        # First requests should have succeeded
        assert statuses[0] == 200


# =====================================================================
# 14-16. Reference endpoints
# =====================================================================

class TestReferenceEndpoints:
    def test_guidelines_returns_list(self, client):
        resp = client.get("/v1/cardio/guidelines")
        assert resp.status_code == 200
        data = resp.json()
        assert "guidelines" in data
        assert isinstance(data["guidelines"], list)
        assert len(data["guidelines"]) > 0
        # Check structure of a guideline entry
        g = data["guidelines"][0]
        assert "id" in g
        assert "title" in g

    def test_guidelines_filter_by_condition(self, client):
        resp = client.get("/v1/cardio/guidelines", params={"condition": "heart_failure"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["filter"] == "heart_failure"
        for g in data["guidelines"]:
            assert "heart_failure" in g["conditions"]

    def test_conditions_returns_list(self, client):
        resp = client.get("/v1/cardio/conditions")
        assert resp.status_code == 200
        data = resp.json()
        assert "conditions" in data
        assert isinstance(data["conditions"], list)
        assert data["total"] == 12
        # Check structure
        c = data["conditions"][0]
        assert "id" in c
        assert "name" in c
        assert "subtypes" in c

    def test_knowledge_version_returns_info(self, client):
        resp = client.get("/v1/cardio/knowledge-version")
        assert resp.status_code == 200
        data = resp.json()
        assert "version" in data
        assert "last_updated" in data
        assert "sources" in data


# =====================================================================
# Additional endpoint coverage
# =====================================================================

class TestAdditionalEndpoints:
    def test_biomarkers_returns_list(self, client):
        resp = client.get("/v1/cardio/biomarkers")
        assert resp.status_code == 200
        data = resp.json()
        assert "biomarkers" in data
        assert len(data["biomarkers"]) > 0

    def test_drugs_returns_list(self, client):
        resp = client.get("/v1/cardio/drugs")
        assert resp.status_code == 200
        data = resp.json()
        assert "drug_classes" in data
        assert len(data["drug_classes"]) > 0

    def test_genes_returns_list(self, client):
        resp = client.get("/v1/cardio/genes")
        assert resp.status_code == 200
        data = resp.json()
        assert "genes" in data
        assert len(data["genes"]) > 0

    def test_heart_score_low_risk(self, client):
        """HEART score = 0+0+0+0+0 = 0 => low."""
        resp = client.post("/v1/cardio/risk/heart-score", json={
            "history_score": 0, "ecg_score": 0,
            "age": 30, "risk_factors": 0, "troponin_score": 0,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["score"] == 0.0
        assert data["risk_category"] == "low"

    def test_heart_score_high_risk(self, client):
        """HEART score = 2+2+2+2+2 = 10 => high."""
        resp = client.post("/v1/cardio/risk/heart-score", json={
            "history_score": 2, "ecg_score": 2,
            "age": 70, "risk_factors": 5, "troponin_score": 2,
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["score"] == 10.0
        assert data["risk_category"] == "high"

    def test_cha2ds2_vasc_low_for_female_only(self, client):
        """Female sex alone = score 1, but categorized as low."""
        resp = client.post("/v1/cardio/risk/cha2ds2-vasc", json={
            "chf": False, "hypertension": False,
            "age": 40, "diabetes": False,
            "stroke_tia": False, "vascular_disease": False,
            "sex": "female",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["score"] == 1.0
        assert data["risk_category"] == "low"

    def test_invalid_ascvd_age_validation(self, client):
        """Age outside 40-79 should be rejected by Pydantic."""
        resp = client.post("/v1/cardio/risk/ascvd", json={
            "age": 25, "sex": "male", "race": "white",
            "total_cholesterol": 200, "hdl_cholesterol": 50,
            "systolic_bp": 120,
        })
        assert resp.status_code == 422  # Validation error
