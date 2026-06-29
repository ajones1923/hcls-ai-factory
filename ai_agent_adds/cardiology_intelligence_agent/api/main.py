"""Cardiology Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine, clinical risk calculators, GDMT
optimizer, and cardiology workflow modules as a production-ready REST API
with CORS, health checks, Prometheus-compatible metrics, and Pydantic
request/response schemas.

Endpoints:
    GET  /health           -- Service health with collection and vector counts
    GET  /collections      -- Collection names and record counts
    GET  /workflows        -- Available clinical workflows
    GET  /metrics          -- Prometheus-compatible metrics (placeholder)

    Versioned routes (via api/routes/):
    POST /v1/cardio/query          -- RAG Q&A query
    POST /v1/cardio/search         -- Multi-collection search
    POST /v1/cardio/find-related   -- Related entity lookup
    POST /v1/cardio/risk/ascvd     -- ASCVD 10-year risk
    POST /v1/cardio/risk/heart-score -- HEART score for ACS
    POST /v1/cardio/risk/cha2ds2-vasc -- CHA2DS2-VASc stroke risk
    POST /v1/cardio/risk/has-bled  -- HAS-BLED bleeding risk
    POST /v1/cardio/risk/maggic    -- MAGGIC HF mortality
    POST /v1/cardio/risk/euroscore -- EuroSCORE II surgical risk
    POST /v1/cardio/gdmt/optimize  -- GDMT optimization
    POST /v1/cardio/workflow/*     -- Clinical workflows (8 types)
    GET  /v1/cardio/guidelines     -- Guideline library
    GET  /v1/cardio/conditions     -- Condition catalogue
    GET  /v1/cardio/biomarkers     -- Cardiac biomarker reference
    GET  /v1/cardio/drugs          -- Drug class reference
    GET  /v1/cardio/genes          -- Cardio-relevant genes
    POST /v1/reports/generate      -- Report generation
    GET  /v1/reports/formats       -- Supported export formats
    GET  /v1/events/stream         -- SSE event stream

Port: 8126 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8126 --reload

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
    or os.environ.get("CARDIO_ANTHROPIC_API_KEY")
)
if _api_key:
    os.environ["ANTHROPIC_API_KEY"] = _api_key

from config.settings import settings

# System prompt for LLM fallback
try:
    from src.agent import CARDIO_SYSTEM_PROMPT as _CARDIO_SYSTEM_PROMPT
except ImportError:
    _CARDIO_SYSTEM_PROMPT = (
        "You are a cardiology clinical decision support system. "
        "Provide evidence-based cardiovascular recommendations citing "
        "ACC/AHA/ESC guidelines with Class of Recommendation and Level of Evidence."
    )

# Route modules
from api.routes.cardio_clinical import router as clinical_router
from api.routes.reports import router as reports_router
from api.routes.events import router as events_router

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine = None          # CardioRAGEngine
_agent = None           # CardiologyIntelligenceAgent
_manager = None         # Collection manager
_risk_calculators = {}  # name -> calculator instance
_gdmt_optimizer = None  # GDMT optimizer module

# Simple request counters for /metrics
_metrics: Dict[str, int] = {
    "requests_total": 0,
    "query_requests_total": 0,
    "search_requests_total": 0,
    "risk_calc_requests_total": 0,
    "workflow_requests_total": 0,
    "gdmt_requests_total": 0,
    "report_requests_total": 0,
    "errors_total": 0,
}
_metrics_lock = threading.Lock()


# =====================================================================
# Lightweight Milvus collection manager (replaces missing CardioCollectionManager)
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
        except Exception:
            pass

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
# Lifespan -- initialize engine on startup, disconnect on shutdown
# =====================================================================

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize the RAG engine, risk calculators, and Milvus on startup."""
    global _engine, _agent, _manager, _risk_calculators, _gdmt_optimizer

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
                    system=system_prompt or _CARDIO_SYSTEM_PROMPT,
                    messages=messages,
                )
                return resp.content[0].text

        llm_client = _LLMClient()
        logger.info("Anthropic LLM client initialized")
    except Exception as exc:
        llm_client = None
        logger.warning(f"LLM client unavailable: {exc}")

    # -- RAG engine --
    try:
        from src.rag_engine import CardioRAGEngine
        _engine = CardioRAGEngine(
            embedding_model=embedder,
            llm_client=llm_client,
            milvus_client=_manager,
        )
        logger.info("Cardiology RAG engine initialized")
    except Exception as exc:
        logger.warning(f"RAG engine unavailable: {exc}")
        _engine = None

    # -- Risk calculators --
    try:
        from src.risk_calculators import RiskCalculatorEngine
        _calc_engine = RiskCalculatorEngine()
        _risk_calculators = {
            "engine": _calc_engine,
            "ascvd": _calc_engine,
            "heart": _calc_engine,
            "cha2ds2_vasc": _calc_engine,
            "has_bled": _calc_engine,
            "maggic": _calc_engine,
            "euroscore": _calc_engine,
        }
        logger.info("Risk calculator engine initialized (6 calculators)")
    except Exception as exc:
        logger.warning(f"Risk calculators unavailable: {exc}")

    # -- GDMT optimizer --
    try:
        from src.gdmt_optimizer import GDMTOptimizer
        _gdmt_optimizer = GDMTOptimizer()
        logger.info("GDMT optimizer initialized")
    except Exception as exc:
        logger.warning(f"GDMT optimizer unavailable: {exc}")

    # -- Inject into app state so routes can access them --
    app.state.engine = _engine
    app.state.manager = _manager
    app.state.risk_calculators = _risk_calculators
    app.state.gdmt_optimizer = _gdmt_optimizer
    app.state.metrics = _metrics
    app.state.metrics_lock = _metrics_lock

    yield  # ── Application runs here ──

    # -- Shutdown --
    if _manager:
        try:
            _manager.disconnect()
            logger.info("Milvus disconnected")
        except Exception:
            pass
    logger.info("Cardiology Intelligence Agent shut down")


# =====================================================================
# Application factory
# =====================================================================

app = FastAPI(
    title="Cardiology Intelligence Agent API",
    description=(
        "RAG-powered cardiovascular clinical decision support with "
        "risk calculators, GDMT optimization, and multi-workflow "
        "assessment for coronary, heart failure, arrhythmia, valvular, "
        "imaging, and cardio-oncology domains."
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
    gdmt_ready = _gdmt_optimizer is not None
    all_healthy = milvus_connected and engine_ready

    return {
        "status": "healthy" if all_healthy else "degraded",
        "agent": "cardiology-intelligence-agent",
        "version": "1.0.0",
        "components": {
            "milvus": "connected" if milvus_connected else "unavailable",
            "rag_engine": "ready" if engine_ready else "unavailable",
            "gdmt_optimizer": "ready" if gdmt_ready else "unavailable",
            "risk_calculators": "ready" if _risk_calculators else "unavailable",
        },
        "collections": collection_count,
        "total_vectors": vector_count,
        "workflows": 11,
        "risk_calculators": len(_risk_calculators),
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
    """Return available clinical workflow definitions."""
    return {
        "workflows": [
            {
                "id": "cad",
                "name": "Coronary Artery Disease Assessment",
                "description": "Calcium score, CAD-RADS, plaque characterization, revascularization planning",
                "risk_calculators": ["ascvd"],
            },
            {
                "id": "heart_failure",
                "name": "Heart Failure Management",
                "description": "HFrEF/HFpEF/HFmrEF classification, GDMT optimization, device evaluation",
                "risk_calculators": ["maggic"],
            },
            {
                "id": "valvular",
                "name": "Valvular Heart Disease",
                "description": "Severity grading, intervention timing, prosthetic valve follow-up",
                "risk_calculators": ["euroscore"],
            },
            {
                "id": "arrhythmia",
                "name": "Arrhythmia & EP Assessment",
                "description": "AF management, CHA2DS2-VASc, ablation candidacy, device programming",
                "risk_calculators": ["cha2ds2_vasc", "has_bled"],
            },
            {
                "id": "cardiac_mri",
                "name": "Cardiac MRI Interpretation",
                "description": "LGE patterns, T1/T2 mapping, strain analysis, tissue characterization",
                "risk_calculators": [],
            },
            {
                "id": "stress_test",
                "name": "Stress Testing Protocol",
                "description": "Exercise/pharmacologic stress, Duke treadmill score, perfusion defects",
                "risk_calculators": ["ascvd"],
            },
            {
                "id": "prevention",
                "name": "Cardiovascular Prevention",
                "description": "Lipid management, ASCVD risk, statin/PCSK9i selection, lifestyle Rx",
                "risk_calculators": ["ascvd"],
            },
            {
                "id": "cardio_oncology",
                "name": "Cardio-Oncology Surveillance",
                "description": "Chemotherapy cardiotoxicity, GLS tracking, biomarker monitoring",
                "risk_calculators": [],
            },
            {
                "id": "acute_decompensated_hf",
                "name": "Acute Decompensated Heart Failure",
                "description": "Hemodynamic profiling (warm-wet/cold-wet), IV diuretics, inotropes, MCS escalation",
                "risk_calculators": ["maggic"],
            },
            {
                "id": "post_mi",
                "name": "Post-MI Secondary Prevention",
                "description": "Reperfusion assessment, DAPT strategy, beta-blocker/statin/ACEi, cardiac rehab, ICD timing",
                "risk_calculators": ["ascvd", "heart"],
            },
            {
                "id": "myocarditis_pericarditis",
                "name": "Myocarditis & Pericarditis",
                "description": "Lake Louise CMR criteria, biopsy indications, NSAIDs/colchicine, activity restriction",
                "risk_calculators": [],
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
    # Fallback to simple counters
    lines = []
    with _metrics_lock:
        for key, val in _metrics.items():
            lines.append(f"# TYPE cardio_agent_{key} counter")
            lines.append(f"cardio_agent_{key} {val}")
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
        content={"error": exc.detail, "agent": "cardiology-intelligence-agent"},
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    logger.error(f"Unhandled exception: {exc}")
    with _metrics_lock:
        _metrics["errors_total"] += 1
    return JSONResponse(
        status_code=500,
        content={"error": "Internal server error", "agent": "cardiology-intelligence-agent"},
    )
