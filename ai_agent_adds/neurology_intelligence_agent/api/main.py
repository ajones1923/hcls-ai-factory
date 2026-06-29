"""Neurology Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine, clinical scale calculators,
neurology-specific workflows (stroke triage, dementia evaluation,
epilepsy classification, tumor grading, MS assessment, Parkinson's
assessment, headache classification, neuromuscular evaluation), and
reference catalogues as a production-ready REST API.

Endpoints:
    GET  /health           -- Service health with collection and vector counts
    GET  /collections      -- Collection names and record counts
    GET  /workflows        -- Available neurology workflows
    GET  /metrics          -- Prometheus-compatible metrics (placeholder)

    Versioned routes (via api/routes/):
    POST /v1/neuro/query               -- RAG Q&A query
    POST /v1/neuro/search              -- Multi-collection search
    POST /v1/neuro/scale/calculate     -- Clinical scale calculator
    POST /v1/neuro/stroke/triage       -- Acute stroke triage
    POST /v1/neuro/dementia/evaluate   -- Dementia evaluation
    POST /v1/neuro/epilepsy/classify   -- Epilepsy classification
    POST /v1/neuro/tumor/grade         -- Brain tumor grading
    POST /v1/neuro/ms/assess           -- MS assessment
    POST /v1/neuro/parkinsons/assess   -- Parkinson's assessment
    POST /v1/neuro/headache/classify   -- Headache classification
    POST /v1/neuro/neuromuscular/evaluate -- Neuromuscular evaluation
    POST /v1/neuro/workflow/{workflow_type} -- Generic workflow dispatch
    GET  /v1/neuro/domains             -- Neurology domain catalogue
    GET  /v1/neuro/scales              -- Available clinical scales
    GET  /v1/neuro/guidelines          -- Guideline reference
    GET  /v1/neuro/knowledge-version   -- Version metadata
    POST /v1/reports/generate          -- Report generation
    GET  /v1/reports/formats           -- Supported export formats
    GET  /v1/events/stream             -- SSE event stream

Port: 8528 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8528 --reload

Author: Adam Jones
Date: March 2026
"""

import os
import sys
import time
import threading
from collections import defaultdict
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Dict, List

from fastapi import FastAPI, HTTPException, Request
from loguru import logger
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, PlainTextResponse

# =====================================================================
# Path setup -- ensure project root is importable
# =====================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Load API key from environment variables
_api_key = (
    os.environ.get("ANTHROPIC_API_KEY")
    or os.environ.get("NEURO_ANTHROPIC_API_KEY")
)
if _api_key:
    os.environ["ANTHROPIC_API_KEY"] = _api_key

from config.settings import settings

# System prompt for LLM fallback
_NEURO_SYSTEM_PROMPT = (
    "You are a neurology intelligence system. "
    "Provide evidence-based recommendations for neurological diagnosis, "
    "treatment, and monitoring citing AAN, AHA/ASA, ILAE, IHS, and NCCN "
    "guidelines and published neurological research."
)

# Route modules
from api.routes.neuro_clinical import router as clinical_router
from api.routes.reports import router as reports_router
from api.routes.events import router as events_router

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine = None          # NeuroRAGEngine
_manager = None         # Collection manager
_workflow_engine = None  # Workflow engine

# Simple request counters for /metrics
_metrics: Dict[str, int] = {
    "requests_total": 0,
    "query_requests_total": 0,
    "search_requests_total": 0,
    "scale_requests_total": 0,
    "workflow_requests_total": 0,
    "report_requests_total": 0,
    "errors_total": 0,
}
_metrics_lock = threading.Lock()


# =====================================================================
# Lightweight Milvus collection manager
# =====================================================================

class _CollectionManager:
    """Thin wrapper around pymilvus for collection management."""

    def __init__(self, host: str = "localhost", port: int = 19530):
        self.host = host
        self.port = port
        self._connections = None

    def connect(self):
        """Connect to Milvus. Degrades gracefully if pymilvus is absent."""
        try:
            from pymilvus import connections
            self._connections = connections
            connections.connect(alias="default", host=self.host, port=str(self.port))
        except Exception as exc:
            logger.warning(f"_CollectionManager.connect failed: {exc}")
            self._connections = None

    def disconnect(self):
        """Disconnect from Milvus if connected."""
        try:
            if self._connections is not None:
                self._connections.disconnect(alias="default")
        except Exception as exc:
            logger.debug("Milvus disconnect ignored: %s", exc)

    def list_collections(self) -> List[str]:
        """Return collection names from Milvus."""
        try:
            from pymilvus import utility
            return utility.list_collections()
        except Exception:
            return []

    def get_stats(self) -> Dict[str, int]:
        """Return dict with collection_count and total_vectors."""
        try:
            from pymilvus import Collection, utility
            names = utility.list_collections()
            total = 0
            for name in names:
                try:
                    col = Collection(name)
                    total += col.num_entities
                except Exception:
                    pass
            return {"collection_count": len(names), "total_vectors": total}
        except Exception:
            return {"collection_count": 0, "total_vectors": 0}


# =====================================================================
# Lightweight workflow engine
# =====================================================================

class _WorkflowEngine:
    """Thin workflow dispatcher for neurology workflows."""

    WORKFLOW_TYPES = [
        "acute_stroke_triage", "dementia_evaluation", "epilepsy_focus",
        "brain_tumor_grading", "ms_monitoring", "parkinsons_assessment",
        "headache_classification", "neuromuscular_evaluation", "general",
    ]

    def __init__(self, llm_client=None, rag_engine=None):
        self.llm_client = llm_client
        self.rag_engine = rag_engine

    def list_workflows(self) -> List[Dict]:
        """Return workflow definitions."""
        return [
            {"id": wf, "name": wf.replace("_", " ").title()}
            for wf in self.WORKFLOW_TYPES
        ]

    def execute(self, workflow_type: str, data: dict) -> dict:
        """Execute a workflow. Falls back to LLM if no dedicated engine."""
        if self.llm_client and self.rag_engine:
            context = ""
            try:
                results = self.rag_engine.search(
                    data.get("question", data.get("query", str(data))),
                    top_k=5,
                )
                context = "\n".join(
                    r.get("content", r.get("text", "")) for r in results
                )
            except Exception:
                pass

            prompt = (
                f"Neurology workflow: {workflow_type}\n\n"
                f"Input data:\n{data}\n\n"
                f"Relevant evidence:\n{context}\n\n"
                f"Provide a detailed clinical analysis and recommendations."
            )
            try:
                answer = self.llm_client.generate(
                    prompt, system_prompt=_NEURO_SYSTEM_PROMPT,
                )
                return {
                    "workflow_type": workflow_type,
                    "status": "completed",
                    "result": answer,
                    "evidence_used": bool(context),
                }
            except Exception as exc:
                logger.warning(f"LLM workflow execution failed: {exc}")

        return {
            "workflow_type": workflow_type,
            "status": "completed",
            "result": f"Workflow '{workflow_type}' executed with provided data.",
            "note": "LLM unavailable; returning placeholder result.",
        }


# =====================================================================
# Lifespan -- initialize engine on startup, disconnect on shutdown
# =====================================================================

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize the RAG engine, workflow engine, and Milvus on startup."""
    global _engine, _manager, _workflow_engine

    # -- Collection manager --
    try:
        _manager = _CollectionManager(
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        _manager.connect()
        logger.info("Collection manager connected to Milvus")
    except Exception as exc:
        logger.warning(f"Collection manager unavailable: {exc}")
        _manager = None

    # -- Embedder --
    try:
        from sentence_transformers import SentenceTransformer

        class _Embedder:
            def __init__(self):
                self.model = SentenceTransformer(settings.EMBEDDING_MODEL)

            def embed_text(self, text: str) -> List[float]:
                return self.model.encode(text).tolist()

        embedder = _Embedder()
        logger.info(f"Embedding model loaded: {settings.EMBEDDING_MODEL}")
    except ImportError:
        embedder = None
        logger.warning("sentence-transformers not available; embedder disabled")

    # -- LLM client --
    llm_client = None
    try:
        import anthropic

        class _LLMClient:
            def __init__(self):
                self.client = anthropic.Anthropic()

            def generate(
                self, prompt: str, system_prompt: str = "",
                max_tokens: int = 2048, temperature: float = 0.7,
            ) -> str:
                messages = [{"role": "user", "content": prompt}]
                resp = self.client.messages.create(
                    model=settings.LLM_MODEL,
                    max_tokens=max_tokens,
                    temperature=temperature,
                    system=system_prompt or _NEURO_SYSTEM_PROMPT,
                    messages=messages,
                )
                return resp.content[0].text

        llm_client = _LLMClient()
        logger.info("Anthropic LLM client initialized")
    except Exception as exc:
        logger.warning(f"LLM client unavailable: {exc}")

    # -- RAG engine --
    try:
        from src.rag_engine import NeuroRAGEngine
        _engine = NeuroRAGEngine(
            embedding_model=embedder,
            llm_client=llm_client,
            milvus_client=_manager,
        )
        logger.info("Neurology RAG engine initialized")
    except Exception as exc:
        logger.warning(f"RAG engine unavailable: {exc}")
        _engine = None

    # -- Workflow engine --
    _workflow_engine = _WorkflowEngine(
        llm_client=llm_client,
        rag_engine=_engine,
    )
    logger.info("Workflow engine initialized (8 workflows + GENERAL)")

    # -- Inject into app state so routes can access them --
    app.state.engine = _engine
    app.state.manager = _manager
    app.state.workflow_engine = _workflow_engine
    app.state.llm_client = llm_client
    app.state.metrics = _metrics
    app.state.metrics_lock = _metrics_lock

    yield  # -- Application runs here --

    # -- Shutdown --
    if _manager:
        try:
            _manager.disconnect()
            logger.info("Milvus disconnected")
        except Exception as exc:
            logger.debug("Shutdown disconnect ignored: %s", exc)
    logger.info("Neurology Intelligence Agent shut down")


# =====================================================================
# Application factory
# =====================================================================

app = FastAPI(
    title="Neurology Intelligence Agent API",
    description=(
        "RAG-powered neurology clinical decision support with "
        "acute stroke triage, dementia evaluation, epilepsy classification, "
        "brain tumor grading, MS monitoring, Parkinson's assessment, "
        "headache classification, neuromuscular evaluation, validated "
        "clinical scale calculators (NIHSS, GCS, MoCA, UPDRS, EDSS, "
        "mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr), and guideline-"
        "driven recommendations."
    ),
    version="1.0.0",
    lifespan=lifespan,
)

# -- CORS (use configured origins, not wildcard) --
_cors_origins = [o.strip() for o in settings.CORS_ORIGINS.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Content-Type", "Authorization", "X-API-Key", "Accept"],
)

# -- Include versioned routers --
app.include_router(clinical_router)
app.include_router(reports_router)
app.include_router(events_router)


# =====================================================================
# Middleware -- authentication, request limits, metrics
# =====================================================================

_AUTH_SKIP_PATHS = {"/health", "/healthz", "/metrics"}


@app.middleware("http")
async def check_api_key(request: Request, call_next):
    """Validate API key if API_KEY is configured in settings."""
    api_key = settings.API_KEY
    if not api_key:
        return await call_next(request)
    if request.url.path in _AUTH_SKIP_PATHS:
        return await call_next(request)
    provided = request.headers.get("X-API-Key") or request.query_params.get("api_key")
    if provided != api_key:
        return JSONResponse(
            status_code=401,
            content={"detail": "Invalid or missing API key"},
        )
    return await call_next(request)


@app.middleware("http")
async def limit_request_size(request: Request, call_next):
    """Reject request bodies that exceed the configured size limit."""
    content_length = request.headers.get("content-length")
    max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    if content_length:
        try:
            if int(content_length) > max_bytes:
                return JSONResponse(
                    status_code=413,
                    content={"detail": "Request too large"},
                )
        except ValueError:
            pass
    return await call_next(request)


_rate_limit_store: Dict[str, list] = defaultdict(list)
_RATE_LIMIT_MAX = 100  # requests per window
_RATE_LIMIT_WINDOW = 60  # seconds

_RATE_LIMIT_SKIP_PATHS = {"/health", "/healthz", "/metrics"}


@app.middleware("http")
async def rate_limit_middleware(request: Request, call_next):
    """Simple in-memory rate limiting by client IP."""
    if request.url.path in _RATE_LIMIT_SKIP_PATHS:
        return await call_next(request)
    client_ip = request.client.host if request.client else "unknown"
    now = time.time()
    # Clean old entries
    _rate_limit_store[client_ip] = [
        t for t in _rate_limit_store[client_ip] if now - t < _RATE_LIMIT_WINDOW
    ]
    if len(_rate_limit_store[client_ip]) >= _RATE_LIMIT_MAX:
        return JSONResponse(
            status_code=429,
            content={"detail": "Rate limit exceeded. Try again later."},
        )
    _rate_limit_store[client_ip].append(now)
    return await call_next(request)


@app.middleware("http")
async def metrics_middleware(request: Request, call_next):
    """Increment request counter for every inbound request."""
    with _metrics_lock:
        _metrics["requests_total"] += 1
    try:
        response = await call_next(request)
        return response
    except Exception:
        with _metrics_lock:
            _metrics["errors_total"] += 1
        raise


# =====================================================================
# Core endpoints
# =====================================================================

@app.get("/health", tags=["system"])
async def health_check():
    """Service health reflecting actual component readiness."""
    milvus_connected = False
    collection_count = 0
    vector_count = 0
    if _manager:
        try:
            stats = _manager.get_stats()
            collection_count = stats.get("collection_count", 0)
            vector_count = stats.get("total_vectors", 0)
            milvus_connected = collection_count > 0
        except Exception:
            pass

    engine_ready = _engine is not None
    workflow_ready = _workflow_engine is not None
    all_healthy = milvus_connected and engine_ready

    return {
        "status": "healthy" if all_healthy else "degraded",
        "agent": "neurology-intelligence-agent",
        "version": "1.0.0",
        "components": {
            "milvus": "connected" if milvus_connected else "unavailable",
            "rag_engine": "ready" if engine_ready else "unavailable",
            "workflow_engine": "ready" if workflow_ready else "unavailable",
        },
        "collections": collection_count,
        "total_vectors": vector_count,
        "workflows": 9,
        "scales": 10,
    }


@app.get("/collections", tags=["system"])
async def list_collections():
    """Return names and record counts for all loaded collections."""
    if _manager:
        try:
            return {"collections": _manager.list_collections()}
        except Exception as exc:
            logger.error(f"Failed to list collections: {exc}")
            raise HTTPException(status_code=503, detail="Service temporarily unavailable")

    raise HTTPException(
        status_code=503,
        detail="Service temporarily unavailable",
    )


@app.get("/workflows", tags=["system"])
async def list_workflows():
    """Return available neurology workflow definitions."""
    return {
        "workflows": [
            {
                "id": "acute_stroke_triage",
                "name": "Acute Stroke Triage",
                "description": "NIHSS-based stroke severity, ASPECTS scoring, thrombolysis/thrombectomy eligibility, and vascular territory localization",
            },
            {
                "id": "dementia_evaluation",
                "name": "Dementia Evaluation",
                "description": "MoCA-based cognitive screening, ATN biomarker staging, differential diagnosis (AD, FTD, LBD, VaD), and treatment planning",
            },
            {
                "id": "epilepsy_focus",
                "name": "Epilepsy Focus Localization",
                "description": "ILAE 2017 seizure classification, syndrome identification, EEG-MRI concordance analysis, and surgical candidacy assessment",
            },
            {
                "id": "brain_tumor_grading",
                "name": "Brain Tumor Grading",
                "description": "WHO 2021 CNS tumor classification with IDH, MGMT, 1p/19q molecular integration and NCCN treatment guidelines",
            },
            {
                "id": "ms_monitoring",
                "name": "MS Disease Monitoring",
                "description": "EDSS scoring, NEDA-3/4 status assessment, DMT escalation evaluation, and relapse risk stratification",
            },
            {
                "id": "parkinsons_assessment",
                "name": "Parkinson's Assessment",
                "description": "MDS-UPDRS Part III motor scoring, Hoehn-Yahr staging, motor subtype classification, and medication optimization",
            },
            {
                "id": "headache_classification",
                "name": "Headache Classification",
                "description": "ICHD-3 criteria-based classification, HIT-6 disability scoring, red flag screening, and preventive therapy selection",
            },
            {
                "id": "neuromuscular_evaluation",
                "name": "Neuromuscular Evaluation",
                "description": "ALSFRS-R functional scoring, EMG/NCS pattern analysis, NMJ localization, and genetic testing guidance",
            },
            {
                "id": "general",
                "name": "General Neurology Query",
                "description": "Free-form RAG-powered Q&A across all neurology knowledge collections",
            },
        ]
    }


@app.get("/metrics", tags=["system"], response_class=PlainTextResponse)
async def prometheus_metrics():
    """Prometheus-compatible metrics export."""
    try:
        from src.metrics import get_metrics_text
        text = get_metrics_text()
        if text and text.strip():
            return text
    except Exception:
        pass
    lines = []
    with _metrics_lock:
        for key, val in _metrics.items():
            lines.append(f"# HELP neuro_agent_{key} Neurology agent {key}")
            lines.append(f"# TYPE neuro_agent_{key} counter")
            lines.append(f"neuro_agent_{key} {val}")
    return "\n".join(lines) + "\n"


# =====================================================================
# Error handlers
# =====================================================================

@app.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    with _metrics_lock:
        _metrics["errors_total"] += 1
    return JSONResponse(
        status_code=exc.status_code,
        content={"error": exc.detail, "agent": "neurology-intelligence-agent"},
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    logger.error(f"Unhandled exception: {exc}")
    with _metrics_lock:
        _metrics["errors_total"] += 1
    return JSONResponse(
        status_code=500,
        content={"error": "Internal server error", "agent": "neurology-intelligence-agent"},
    )
