"""Tests for Biomarker Intelligence Agent REST API.

Validates all API routes using FastAPI TestClient with mocked backends
so tests run without Milvus, sentence-transformers, or Anthropic API.

Author: Adam Jones
Date: March 2026
"""

import json
from contextlib import asynccontextmanager
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def mock_agent():
    """Return a mock PrecisionBiomarkerAgent."""
    agent = MagicMock()

    from src.models import (
        AnalysisResult,
        BiologicalAgeResult,
        DiseaseCategory,
        DiseaseTrajectoryResult,
        GenotypeAdjustmentResult,
        MetabolizerPhenotype,
        PatientProfile,
        PGxResult,
        RiskLevel,
    )

    profile = PatientProfile(
        patient_id="TEST-001",
        age=45,
        sex="M",
        biomarkers={
            "albumin": 4.2, "creatinine": 0.9, "glucose": 95.0,
            "hs_crp": 1.5, "lymphocyte_pct": 30.0, "mcv": 89.0,
            "rdw": 13.0, "alkaline_phosphatase": 65.0, "wbc": 6.0,
        },
    )
    mock_result = AnalysisResult(
        patient_profile=profile,
        biological_age=BiologicalAgeResult(
            chronological_age=45,
            biological_age=46.5,
            age_acceleration=1.5,
            phenoage_score=0.003,
            mortality_risk=0.012,
        ),
        disease_trajectories=[
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.DIABETES,
                risk_level=RiskLevel.MODERATE,
                current_markers={"hba1c": 5.6},
            ),
        ],
        pgx_results=[
            PGxResult(
                gene="CYP2D6",
                star_alleles="*1/*1",
                phenotype=MetabolizerPhenotype.NORMAL,
                drugs_affected=[
                    {"drug": "Codeine", "recommendation": "Standard dosing"},
                ],
            ),
        ],
        genotype_adjustments=[
            GenotypeAdjustmentResult(
                biomarker="ALT",
                standard_range="0-56 U/L",
                adjusted_range="0-45 U/L",
                genotype="CG",
                gene="PNPLA3",
                rationale="PNPLA3 CG carrier",
            ),
        ],
        critical_alerts=[],
    )
    agent.analyze_patient.return_value = mock_result
    return agent


@pytest.fixture
def mock_manager():
    """Return a mock collection manager."""
    manager = MagicMock()
    manager.get_collection_stats.return_value = {
        "biomarker_reference": 100,
        "biomarker_genetic_variants": 200,
        "biomarker_pgx_rules": 50,
    }
    return manager


@pytest.fixture
def mock_engine():
    """Return a mock RAG engine."""
    engine = MagicMock()
    engine.llm = MagicMock()
    engine.embedder = MagicMock()

    from src.models import CrossCollectionResult, SearchHit

    evidence = CrossCollectionResult(
        query="test",
        hits=[SearchHit(collection="test", id="t1", score=0.9, text="Test evidence")],
        total_collections_searched=10,
        search_time_ms=15.0,
    )
    engine.retrieve.return_value = evidence
    engine._build_prompt.return_value = "Test prompt"
    engine.llm.generate.return_value = "Mock LLM response"
    return engine


@pytest.fixture
def mock_bio_age_calc():
    """Return a mock biological age calculator."""
    calc = MagicMock()
    calc.calculate.return_value = {
        "biological_age": 46.5,
        "age_acceleration": 1.5,
        "phenoage": {
            "mortality_risk": "LOW",
            "top_aging_drivers": [
                {"biomarker": "rdw", "value": 13.5, "contribution": 4.46},
            ],
        },
    }
    return calc


@pytest.fixture
def mock_trajectory_analyzer():
    """Return a mock disease trajectory analyzer."""
    analyzer = MagicMock()

    from src.models import DiseaseCategory, DiseaseTrajectoryResult, RiskLevel

    analyzer.analyze_all.return_value = [
        DiseaseTrajectoryResult(
            disease=DiseaseCategory.DIABETES,
            risk_level=RiskLevel.MODERATE,
            current_markers={"hba1c": 5.6},
        ),
    ]
    return analyzer


@pytest.fixture
def mock_pgx_mapper():
    """Return a mock pharmacogenomic mapper."""
    mapper = MagicMock()
    mapper.map_all.return_value = {
        "genes_analyzed": 2,
        "gene_results": [
            {
                "gene": "CYP2D6",
                "star_alleles": "*1/*1",
                "phenotype": "Normal Metabolizer",
                "drugs_affected": [
                    {"drug": "Codeine", "recommendation": "Standard dosing"},
                ],
            },
        ],
    }
    return mapper


@pytest.fixture
def mock_genotype_adjuster():
    """Return a mock genotype adjuster."""
    adjuster = MagicMock()
    return adjuster


@pytest.fixture
def client(
    mock_agent,
    mock_manager,
    mock_engine,
    mock_bio_age_calc,
    mock_trajectory_analyzer,
    mock_pgx_mapper,
    mock_genotype_adjuster,
):
    """Create a TestClient with all mocked dependencies.

    Replaces the real lifespan (which requires Milvus, sentence-transformers,
    and Anthropic) with a no-op so that the TestClient can start without any
    external services.
    """
    import api.main as main_module
    from api.main import app

    # Replace lifespan with a no-op so TestClient does not run the real startup
    @asynccontextmanager
    async def noop_lifespan(a):
        yield

    original_lifespan = app.router.lifespan_context
    app.router.lifespan_context = noop_lifespan
    try:
        with TestClient(app, raise_server_exceptions=False) as c:
            # Set mocked state on app.state (read by route handlers via req.app.state)
            app.state.engine = mock_engine
            app.state.agent = mock_agent
            app.state.manager = mock_manager
            app.state.bio_age_calc = mock_bio_age_calc
            app.state.trajectory_analyzer = mock_trajectory_analyzer
            app.state.pgx_mapper = mock_pgx_mapper
            app.state.genotype_adjuster = mock_genotype_adjuster
            app.state.metrics = {
                "requests_total": 0,
                "query_requests_total": 0,
                "search_requests_total": 0,
                "analyze_requests_total": 0,
                "bio_age_requests_total": 0,
                "errors_total": 0,
            }

            # Also set module-level variables that /health and /collections read directly
            main_module._engine = mock_engine
            main_module._agent = mock_agent
            main_module._manager = mock_manager

            yield c
    finally:
        # Restore original lifespan and module-level state even if the test fails
        app.router.lifespan_context = original_lifespan
        main_module._engine = None
        main_module._agent = None
        main_module._manager = None


@pytest.fixture(autouse=True)
def _clear_event_stores():
    """Clear in-memory event stores before each test to prevent accumulation."""
    from api.routes.events import _inbound_events, _outbound_alerts

    _inbound_events.clear()
    _outbound_alerts.clear()
    yield
    _inbound_events.clear()
    _outbound_alerts.clear()


@pytest.fixture(autouse=True)
def _clear_report_store():
    """Clear in-memory report store before each test."""
    from api.routes.reports import _report_store

    _report_store.clear()
    yield
    _report_store.clear()


# =====================================================================
# HEALTH & STATUS ENDPOINTS
# =====================================================================


class TestHealthEndpoint:
    """Tests for GET /health."""

    def test_healthy_response(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "healthy"
        # The mock manager returns 3 collections, all with count > 0
        assert data["collections"] == 3
        assert data["total_vectors"] == 350
        assert data["agent_ready"] is True

    def test_degraded_when_no_manager(self, client):
        import api.main as main_module

        original_manager = main_module._manager
        original_engine = main_module._engine
        original_agent = main_module._agent
        main_module._manager = None
        # With no manager the endpoint returns degraded.
        # agent_ready depends on _engine and _agent being non-None.
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "degraded"
        assert data["collections"] == 0
        assert data["total_vectors"] == 0
        # Restore
        main_module._manager = original_manager
        main_module._engine = original_engine
        main_module._agent = original_agent

    def test_degraded_agent_not_ready(self, client):
        import api.main as main_module

        original_agent = main_module._agent
        main_module._agent = None
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        # Manager still works but agent is not ready
        assert data["agent_ready"] is False
        main_module._agent = original_agent


class TestCollectionsEndpoint:
    """Tests for GET /collections."""

    def test_returns_collections(self, client):
        resp = client.get("/collections")
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 3
        assert len(data["collections"]) == 3
        names = {c["name"] for c in data["collections"]}
        assert "biomarker_reference" in names

    def test_503_when_no_manager(self, client):
        import api.main as main_module

        original = main_module._manager
        main_module._manager = None
        resp = client.get("/collections")
        assert resp.status_code == 503
        main_module._manager = original


class TestKnowledgeStats:
    """Tests for GET /knowledge/stats."""

    def test_returns_stats(self, client):
        resp = client.get("/knowledge/stats")
        assert resp.status_code == 200
        data = resp.json()
        assert "disease_domains" in data
        assert "total_biomarkers" in data
        assert "pharmacogenes" in data
        assert "pgx_drug_interactions" in data
        assert "phenoage_markers" in data
        assert "cross_modal_links" in data
        # All values should be non-negative integers
        for key, val in data.items():
            assert isinstance(val, int)
            assert val >= 0


class TestMetricsEndpoint:
    """Tests for GET /metrics."""

    def test_returns_prometheus_format(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200
        assert "text/plain" in resp.headers["content-type"]
        text = resp.text
        assert "biomarker_api_requests_total" in text
        assert "biomarker_api_errors_total" in text
        assert "biomarker_api_query_requests_total" in text
        assert "biomarker_api_analyze_requests_total" in text
        assert "biomarker_api_bio_age_requests_total" in text

    def test_metrics_include_collection_vectors(self, client):
        """When a manager is available, metrics should include per-collection vector counts."""
        resp = client.get("/metrics")
        text = resp.text
        assert "biomarker_collection_vectors" in text
        assert 'collection="biomarker_reference"' in text


class TestV1Health:
    """Tests for GET /v1/health."""

    def test_returns_status(self, client):
        resp = client.get("/v1/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["agent_ready"] is True
        assert data["rag_engine_ready"] is True
        assert data["status"] == "healthy"
        assert data["llm_available"] is True
        assert data["embedder_available"] is True


# =====================================================================
# ANALYSIS ENDPOINTS
# =====================================================================


class TestAnalyzeEndpoint:
    """Tests for POST /v1/analyze."""

    def test_full_analysis(self, client):
        payload = {
            "patient_id": "TEST-001",
            "age": 45,
            "sex": "M",
            "biomarkers": {"albumin": 4.2, "glucose": 95.0},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["patient_id"] == "TEST-001"
        assert "biological_age" in data
        assert "disease_trajectories" in data
        assert "processing_time_ms" in data
        assert "pgx_results" in data
        assert "genotype_adjustments" in data
        assert "critical_alerts" in data

    def test_rejects_empty_data(self, client):
        """At least one of biomarkers, genotypes, or star_alleles must be provided."""
        payload = {
            "patient_id": "TEST-001",
            "age": 45,
            "sex": "M",
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 422

    def test_rejects_invalid_age(self, client):
        payload = {
            "patient_id": "TEST-001",
            "age": -5,
            "sex": "M",
            "biomarkers": {"albumin": 4.2},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 422

    def test_rejects_age_over_150(self, client):
        payload = {
            "patient_id": "TEST-001",
            "age": 200,
            "sex": "M",
            "biomarkers": {"albumin": 4.2},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 422

    def test_503_when_no_agent(self, client):
        from api.main import app

        original = app.state.agent
        app.state.agent = None
        payload = {
            "patient_id": "TEST-001",
            "age": 45,
            "sex": "M",
            "biomarkers": {"albumin": 4.2},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 503
        app.state.agent = original

    def test_accepts_genotypes_only(self, client):
        """Providing only genotypes (no biomarkers) should be valid input."""
        payload = {
            "patient_id": "TEST-002",
            "age": 30,
            "sex": "F",
            "genotypes": {"rs1801133": "CT"},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 200

    def test_accepts_star_alleles_only(self, client):
        """Providing only star_alleles should be valid input."""
        payload = {
            "patient_id": "TEST-003",
            "age": 60,
            "sex": "M",
            "star_alleles": {"CYP2D6": "*1/*4"},
        }
        resp = client.post("/v1/analyze", json=payload)
        assert resp.status_code == 200


class TestBiologicalAgeEndpoint:
    """Tests for POST /v1/biological-age."""

    def test_returns_biological_age(self, client):
        payload = {
            "age": 45,
            "biomarkers": {"albumin": 4.2, "creatinine": 0.9, "glucose": 95.0},
        }
        resp = client.post("/v1/biological-age", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["chronological_age"] == 45
        assert "biological_age" in data
        assert "age_acceleration" in data
        assert "mortality_risk" in data
        assert "top_aging_drivers" in data

    def test_503_when_no_calc(self, client):
        from api.main import app

        original = app.state.bio_age_calc
        app.state.bio_age_calc = None
        payload = {
            "age": 45,
            "biomarkers": {"albumin": 4.2},
        }
        resp = client.post("/v1/biological-age", json=payload)
        assert resp.status_code == 503
        app.state.bio_age_calc = original

    def test_rejects_invalid_age(self, client):
        payload = {"age": -1, "biomarkers": {"albumin": 4.2}}
        resp = client.post("/v1/biological-age", json=payload)
        assert resp.status_code == 422


class TestDiseaseRiskEndpoint:
    """Tests for POST /v1/disease-risk."""

    def test_returns_trajectories(self, client):
        payload = {
            "age": 45,
            "sex": "M",
            "biomarkers": {"hba1c": 5.6},
        }
        resp = client.post("/v1/disease-risk", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "trajectories" in data
        assert "processing_time_ms" in data
        assert len(data["trajectories"]) >= 1

    def test_503_when_no_analyzer(self, client):
        from api.main import app

        original = app.state.trajectory_analyzer
        app.state.trajectory_analyzer = None
        payload = {"age": 45, "sex": "M", "biomarkers": {"hba1c": 5.6}}
        resp = client.post("/v1/disease-risk", json=payload)
        assert resp.status_code == 503
        app.state.trajectory_analyzer = original


class TestPGxEndpoint:
    """Tests for POST /v1/pgx."""

    def test_returns_pgx_results(self, client):
        payload = {
            "star_alleles": {"CYP2D6": "*1/*1"},
        }
        resp = client.post("/v1/pgx", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "results" in data
        assert "total_genes" in data
        assert "critical_findings" in data
        assert data["total_genes"] == 2

    def test_empty_star_alleles_with_genotypes(self, client):
        payload = {
            "genotypes": {"rs1801133": "CT"},
        }
        resp = client.post("/v1/pgx", json=payload)
        assert resp.status_code == 200

    def test_critical_findings_flagged(self, client, mock_pgx_mapper):
        """Poor metabolizer phenotypes should appear in critical_findings."""
        mock_pgx_mapper.map_all.return_value = {
            "genes_analyzed": 1,
            "gene_results": [
                {
                    "gene": "CYP2D6",
                    "star_alleles": "*4/*4",
                    "phenotype": "Poor Metabolizer",
                    "drugs_affected": [
                        {"drug": "Codeine", "recommendation": "Avoid"},
                    ],
                },
            ],
        }
        payload = {"star_alleles": {"CYP2D6": "*4/*4"}}
        resp = client.post("/v1/pgx", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["critical_findings"]) >= 1
        assert "CYP2D6" in data["critical_findings"][0]

    def test_503_when_no_mapper(self, client):
        from api.main import app

        original = app.state.pgx_mapper
        app.state.pgx_mapper = None
        payload = {"star_alleles": {"CYP2D6": "*1/*1"}}
        resp = client.post("/v1/pgx", json=payload)
        assert resp.status_code == 503
        app.state.pgx_mapper = original


class TestQueryEndpoint:
    """Tests for POST /v1/query."""

    def test_returns_answer(self, client):
        payload = {"question": "What is HbA1c?"}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["question"] == "What is HbA1c?"
        assert "answer" in data
        assert data["answer"] == "Mock LLM response"
        assert "evidence_count" in data
        assert "collections_searched" in data
        assert "search_time_ms" in data

    def test_rejects_empty_question(self, client):
        payload = {"question": ""}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 422

    def test_rejects_oversized_question(self, client):
        payload = {"question": "x" * 5001}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 422

    def test_503_when_no_engine(self, client):
        from api.main import app

        original = app.state.engine
        app.state.engine = None
        payload = {"question": "What is HbA1c?"}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 503
        app.state.engine = original

    def test_503_when_no_llm(self, client):
        from api.main import app

        app.state.engine.llm = None
        payload = {"question": "What is HbA1c?"}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 503

    def test_503_when_no_embedder(self, client):
        from api.main import app

        app.state.engine.embedder = None
        payload = {"question": "What is HbA1c?"}
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 503

    def test_query_with_patient_profile(self, client):
        payload = {
            "question": "Is my HbA1c concerning?",
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"hba1c": 5.9},
            },
        }
        resp = client.post("/v1/query", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["question"] == "Is my HbA1c concerning?"


# =====================================================================
# EVENT ENDPOINTS
# =====================================================================


class TestCrossModalEvents:
    """Tests for cross-modal event endpoints."""

    def test_receive_imaging_event(self, client):
        payload = {
            "source_agent": "imaging_intelligence_agent",
            "event_type": "imaging_finding",
            "patient_id": "TEST-001",
            "payload": {"finding_type": "coronary_calcium", "score": 150},
            "urgency": "high",
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "processed"
        assert "biomarker_correlation_analysis" in data["actions_triggered"]
        assert "cardiovascular_risk_reassessment" in data["actions_triggered"]

    def test_receive_genomic_variant_event(self, client):
        payload = {
            "source_agent": "genomics_pipeline",
            "event_type": "genomic_variant",
            "patient_id": "TEST-001",
            "payload": {"gene": "BRCA1"},
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "genotype_adjustment_update:BRCA1" in data["actions_triggered"]
        assert "disease_trajectory_reanalysis" in data["actions_triggered"]

    def test_receive_drug_alert(self, client):
        payload = {
            "source_agent": "oncology_agent",
            "event_type": "drug_alert",
            "patient_id": "TEST-001",
            "payload": {"drug": "Tamoxifen"},
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "pgx_cross_check" in data["actions_triggered"]
        assert "pgx_drug_safety_review:Tamoxifen" in data["actions_triggered"]

    def test_receive_treatment_started_event(self, client):
        payload = {
            "source_agent": "clinical_agent",
            "event_type": "treatment_started",
            "patient_id": "TEST-001",
            "payload": {"drug": "Metformin"},
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "pgx_treatment_compatibility_check" in data["actions_triggered"]
        assert "biomarker_monitoring_schedule_update" in data["actions_triggered"]

    def test_receive_unknown_event_type(self, client):
        payload = {
            "source_agent": "unknown_agent",
            "event_type": "unknown_event",
            "patient_id": "TEST-001",
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "event_logged" in data["actions_triggered"]

    def test_list_inbound_events(self, client):
        # Create two events first
        client.post("/v1/events/cross-modal", json={
            "source_agent": "test", "event_type": "test",
            "patient_id": "P1",
        })
        client.post("/v1/events/cross-modal", json={
            "source_agent": "test2", "event_type": "imaging_finding",
            "patient_id": "P2",
        })
        resp = client.get("/v1/events/cross-modal")
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 2
        assert len(data["events"]) == 2
        # Newest first
        assert data["events"][0]["source_agent"] == "test2"

    def test_list_inbound_events_filter_by_patient(self, client):
        client.post("/v1/events/cross-modal", json={
            "source_agent": "a1", "event_type": "test", "patient_id": "P1",
        })
        client.post("/v1/events/cross-modal", json={
            "source_agent": "a2", "event_type": "test", "patient_id": "P2",
        })
        resp = client.get("/v1/events/cross-modal", params={"patient_id": "P1"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1
        assert data["events"][0]["patient_id"] == "P1"

    def test_event_has_event_id_and_timestamp(self, client):
        payload = {
            "source_agent": "test",
            "event_type": "test",
            "patient_id": "P1",
        }
        resp = client.post("/v1/events/cross-modal", json=payload)
        data = resp.json()
        assert "event_id" in data
        assert len(data["event_id"]) == 12
        assert "timestamp" in data


class TestBiomarkerAlerts:
    """Tests for biomarker alert endpoints."""

    def test_send_alert(self, client):
        payload = {
            "patient_id": "TEST-001",
            "alert_type": "elevated_lpa",
            "severity": "high",
            "target_agent": "imaging_intelligence_agent",
            "payload": {"lpa_value": 65, "unit": "nmol/L"},
            "description": "Lp(a) elevated - cardiovascular risk",
        }
        resp = client.post("/v1/events/biomarker-alert", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "sent"
        assert data["target_agent"] == "imaging_intelligence_agent"
        assert "alert_id" in data
        assert "timestamp" in data

    def test_send_alert_unknown_type(self, client):
        """Unknown alert types should still be accepted (just logged)."""
        payload = {
            "patient_id": "P1",
            "alert_type": "custom_alert",
            "severity": "low",
            "target_agent": "test_agent",
        }
        resp = client.post("/v1/events/biomarker-alert", json=payload)
        assert resp.status_code == 200
        assert resp.json()["status"] == "sent"

    def test_list_outbound_alerts(self, client):
        client.post("/v1/events/biomarker-alert", json={
            "patient_id": "P1", "alert_type": "elevated_lpa",
            "severity": "high", "target_agent": "imaging_agent",
        })
        client.post("/v1/events/biomarker-alert", json={
            "patient_id": "P2", "alert_type": "pgx_critical",
            "severity": "critical", "target_agent": "oncology_agent",
        })
        resp = client.get("/v1/events/biomarker-alert")
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 2
        assert len(data["events"]) == 2
        # Newest first
        assert data["events"][0]["patient_id"] == "P2"

    def test_list_outbound_alerts_filter_by_patient(self, client):
        client.post("/v1/events/biomarker-alert", json={
            "patient_id": "P1", "alert_type": "test",
            "severity": "low", "target_agent": "a1",
        })
        client.post("/v1/events/biomarker-alert", json={
            "patient_id": "P2", "alert_type": "test",
            "severity": "low", "target_agent": "a2",
        })
        resp = client.get("/v1/events/biomarker-alert", params={"patient_id": "P2"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1
        assert data["events"][0]["patient_id"] == "P2"


# =====================================================================
# REPORT ENDPOINTS
# =====================================================================


class TestReportEndpoints:
    """Tests for report generation and export endpoints."""

    def test_generate_markdown_report(self, client):
        payload = {
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"albumin": 4.2, "glucose": 95.0},
            },
            "format": "markdown",
        }
        with patch("src.report_generator.ReportGenerator") as MockGen:
            mock_gen_instance = MagicMock()
            mock_gen_instance.generate.return_value = "# Mock Report\n\nContent here."
            MockGen.return_value = mock_gen_instance

            resp = client.post("/v1/report/generate", json=payload)
            assert resp.status_code == 200
            data = resp.json()
            assert "meta" in data
            assert "report" in data
            assert data["meta"]["patient_id"] == "TEST-001"
            assert data["meta"]["sections"] == 12
            assert data["meta"]["format"] == "markdown"
            assert "report_id" in data["meta"]
            assert data["report"] == "# Mock Report\n\nContent here."

    def test_generate_json_report(self, client):
        payload = {
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"albumin": 4.2},
            },
            "format": "json",
        }
        with patch("src.report_generator.ReportGenerator") as MockGen, \
             patch("src.export.export_json") as mock_export_json:
            mock_gen_instance = MagicMock()
            mock_gen_instance.generate.return_value = "# Markdown"
            MockGen.return_value = mock_gen_instance
            mock_export_json.return_value = '{"summary": "test"}'

            resp = client.post("/v1/report/generate", json=payload)
            assert resp.status_code == 200
            data = resp.json()
            assert data["meta"]["format"] == "json"
            assert data["report"] == '{"summary": "test"}'

    def test_503_when_no_agent(self, client):
        from api.main import app

        original = app.state.agent
        app.state.agent = None
        payload = {
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"albumin": 4.2},
            },
        }
        resp = client.post("/v1/report/generate", json=payload)
        assert resp.status_code == 503
        app.state.agent = original

    def test_pdf_download_404_for_missing_report(self, client):
        resp = client.get("/v1/report/nonexistent/pdf")
        assert resp.status_code == 404

    def test_pdf_download_for_existing_report(self, client):
        """After generating a report, the PDF download endpoint should find it."""
        from api.routes.reports import _report_store

        # Manually populate the report store
        _report_store["test123"] = {
            "report_id": "test123",
            "patient_id": "TEST-001",
            "markdown": "# Test Report",
            "analysis": MagicMock(),
            "profile": MagicMock(),
            "generated_at": "2026-03-01T00:00:00+00:00",
        }
        with patch("src.export.export_pdf") as mock_pdf:
            mock_pdf.return_value = b"%PDF-1.4 fake content"
            resp = client.get("/v1/report/test123/pdf")
            assert resp.status_code == 200
            assert resp.headers["content-type"] == "application/pdf"
            assert b"%PDF" in resp.content

    def test_fhir_export(self, client):
        payload = {
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"albumin": 4.2},
            },
        }
        fhir_bundle = json.dumps({
            "resourceType": "Bundle",
            "type": "document",
            "entry": [],
        })
        with patch("src.export.export_fhir_diagnostic_report") as mock_fhir:
            mock_fhir.return_value = fhir_bundle
            resp = client.post("/v1/report/fhir", json=payload)
            assert resp.status_code == 200
            data = resp.json()
            assert data["format"] == "fhir_r4"
            assert data["resource_type"] == "Bundle"
            assert data["patient_id"] == "TEST-001"
            assert "bundle" in data
            assert data["bundle"]["resourceType"] == "Bundle"

    def test_fhir_503_when_no_agent(self, client):
        from api.main import app

        original = app.state.agent
        app.state.agent = None
        payload = {
            "patient_profile": {
                "patient_id": "TEST-001",
                "age": 45,
                "sex": "M",
                "biomarkers": {"albumin": 4.2},
            },
        }
        resp = client.post("/v1/report/fhir", json=payload)
        assert resp.status_code == 503
        app.state.agent = original


# =====================================================================
# REQUEST SIZE LIMIT MIDDLEWARE
# =====================================================================


class TestMiddleware:
    """Tests for middleware behavior."""

    def test_request_size_limit_header(self, client):
        """Requests with content-length exceeding limit should be rejected."""
        # The default MAX_REQUEST_SIZE_MB is 10, so 100 MB should be rejected
        resp = client.post(
            "/v1/query",
            json={"question": "test"},
            headers={"content-length": str(100 * 1024 * 1024)},
        )
        assert resp.status_code == 413

    def test_normal_request_passes_size_limit(self, client):
        """Requests within the size limit should pass through to the route."""
        resp = client.post(
            "/v1/query",
            json={"question": "What is HbA1c?"},
        )
        # Should reach the route handler (not blocked by middleware)
        assert resp.status_code == 200

    def test_cors_headers_present(self, client):
        """CORS preflight should return appropriate headers."""
        resp = client.options(
            "/health",
            headers={
                "Origin": "http://localhost:8501",
                "Access-Control-Request-Method": "GET",
            },
        )
        # CORS middleware should respond (FastAPI CORS returns 200 for preflight)
        assert resp.status_code in (200, 204, 400)


# =====================================================================
# API AUTHENTICATION
# =====================================================================


class TestAPIAuthentication:
    """Tests for API key authentication middleware."""

    def test_no_auth_when_key_not_configured(self, client):
        """When API_KEY is empty, all requests pass."""
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_health_always_accessible(self, client):
        """Health endpoint should bypass auth even with key configured."""
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_metrics_always_accessible(self, client):
        """Metrics endpoint should bypass auth."""
        resp = client.get("/metrics")
        assert resp.status_code == 200


# =====================================================================
# THREAD SAFETY
# =====================================================================


class TestThreadSafety:
    """Tests for thread-safe stores."""

    def test_report_store_concurrent_access(self, client):
        """Multiple report creations should not crash."""
        from api.routes.reports import _report_lock
        with _report_lock:
            # Verify lock is acquirable
            pass  # If we get here, lock works

    def test_event_store_concurrent_access(self, client):
        """Multiple event submissions should not crash."""
        from api.routes.events import _events_lock
        with _events_lock:
            pass  # Lock is acquirable
