"""Single-Cell Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine, single-cell analysis workflows
(cell type annotation, TME profiling, drug response prediction,
subclonal analysis, spatial niche mapping, trajectory inference,
ligand-receptor interaction, biomarker discovery, CAR-T target
validation, treatment monitoring), and reference catalogues as a
production-ready REST API.

Endpoints:
    GET  /health           -- Service health with collection and vector counts
    GET  /collections      -- Collection names and record counts
    GET  /workflows        -- Available single-cell workflows
    GET  /metrics          -- Prometheus-compatible metrics (placeholder)

    Versioned routes (via api/routes/):
    POST /v1/sc/query                -- RAG Q&A query
    POST /v1/sc/search               -- Multi-collection search
    POST /v1/sc/annotate             -- Cell type annotation
    POST /v1/sc/tme-profile          -- TME profiling
    POST /v1/sc/drug-response        -- Drug response prediction
    POST /v1/sc/subclonal            -- Subclonal architecture analysis
    POST /v1/sc/spatial-niche        -- Spatial niche mapping
    POST /v1/sc/trajectory           -- Trajectory inference
    POST /v1/sc/ligand-receptor      -- Ligand-receptor interaction
    POST /v1/sc/biomarker            -- Biomarker discovery
    POST /v1/sc/cart-validate        -- CAR-T target validation
    POST /v1/sc/treatment-monitor    -- Treatment monitoring
    POST /v1/sc/workflow/{workflow_type} -- Generic workflow dispatch
    GET  /v1/sc/cell-types           -- Cell type catalogue
    GET  /v1/sc/markers              -- Marker gene reference
    GET  /v1/sc/tme-classes          -- TME classification reference
    GET  /v1/sc/spatial-platforms     -- Spatial platform reference
    GET  /v1/sc/knowledge-version    -- Version metadata
    POST /v1/reports/generate        -- Report generation
    GET  /v1/reports/formats         -- Supported export formats
    GET  /v1/events/stream           -- SSE event stream

Port: 8540 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8540 --reload

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
    or os.environ.get("SC_ANTHROPIC_API_KEY")
)
if _api_key:
    os.environ["ANTHROPIC_API_KEY"] = _api_key

from config.settings import settings

# System prompt for LLM fallback
_SC_SYSTEM_PROMPT = (
    "You are a single-cell genomics intelligence system. "
    "Provide evidence-based analysis for single-cell RNA-seq, spatial "
    "transcriptomics, tumor microenvironment profiling, cell type "
    "annotation, trajectory inference, and drug response prediction, "
    "citing Human Cell Atlas, CellxGene, and published single-cell "
    "research."
)

# Route modules
from api.routes.sc_clinical import router as clinical_router
from api.routes.reports import router as reports_router
from api.routes.events import router as events_router

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine = None          # SingleCellRAGEngine
_manager = None         # Collection manager
_workflow_engine = None  # Workflow engine

# Simple request counters for /metrics
_metrics: Dict[str, int] = {
    "requests_total": 0,
    "query_requests_total": 0,
    "search_requests_total": 0,
    "annotate_requests_total": 0,
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
    """Thin workflow dispatcher for single-cell analysis workflows."""

    WORKFLOW_TYPES = [
        "cell_type_annotation", "tme_profiling", "drug_response",
        "subclonal_analysis", "spatial_niche", "trajectory_inference",
        "ligand_receptor", "biomarker_discovery", "cart_validation",
        "treatment_monitoring",
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

    async def execute(self, workflow_type: str, data: dict) -> dict:
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
                f"Single-cell analysis workflow: {workflow_type}\n\n"
                f"Input data:\n{data}\n\n"
                f"Relevant evidence:\n{context}\n\n"
                f"Provide a detailed single-cell analysis and recommendations."
            )
            try:
                answer = self.llm_client.generate(
                    prompt, system_prompt=_SC_SYSTEM_PROMPT,
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
                    system=system_prompt or _SC_SYSTEM_PROMPT,
                    messages=messages,
                )
                return resp.content[0].text

        llm_client = _LLMClient()
        logger.info("Anthropic LLM client initialized")
    except Exception as exc:
        logger.warning(f"LLM client unavailable: {exc}")

    # -- RAG engine --
    try:
        from src.rag_engine import SingleCellRAGEngine
        _engine = SingleCellRAGEngine(
            embedding_model=embedder,
            llm_client=llm_client,
            milvus_client=_manager,
        )
        logger.info("Single-Cell RAG engine initialized")
    except Exception as exc:
        logger.warning(f"RAG engine unavailable: {exc}")
        _engine = None

    # -- Workflow engine --
    _workflow_engine = _WorkflowEngine(
        llm_client=llm_client,
        rag_engine=_engine,
    )
    logger.info("Workflow engine initialized (10 workflows)")

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
    logger.info("Single-Cell Intelligence Agent shut down")


# =====================================================================
# Application factory
# =====================================================================

app = FastAPI(
    title="Single-Cell Intelligence Agent API",
    description=(
        "RAG-powered single-cell genomics clinical decision support with "
        "cell type annotation, tumor microenvironment profiling, drug "
        "response prediction, subclonal architecture analysis, spatial "
        "niche mapping, trajectory inference, ligand-receptor interaction "
        "analysis, biomarker discovery, CAR-T target validation, and "
        "treatment monitoring."
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
        "agent": "single-cell-intelligence-agent",
        "version": "1.0.0",
        "components": {
            "milvus": "connected" if milvus_connected else "unavailable",
            "rag_engine": "ready" if engine_ready else "unavailable",
            "workflow_engine": "ready" if workflow_ready else "unavailable",
        },
        "collections": collection_count,
        "total_vectors": vector_count,
        "workflows": 10,
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
    """Return available single-cell workflow definitions."""
    return {
        "workflows": [
            {
                "id": "cell_type_annotation",
                "name": "Cell Type Annotation",
                "description": "Multi-strategy cell type annotation using reference-based, marker-based, and LLM-augmented approaches with consensus scoring",
            },
            {
                "id": "tme_profiling",
                "name": "TME Profiling",
                "description": "Tumor microenvironment classification (hot/cold/excluded/suppressed) with immune infiltrate quantification and therapy response prediction",
            },
            {
                "id": "drug_response",
                "name": "Drug Response Prediction",
                "description": "Cell-type-resolved drug sensitivity prediction using GDSC/DepMap signatures with resistance mechanism identification",
            },
            {
                "id": "subclonal_analysis",
                "name": "Subclonal Architecture",
                "description": "CNV-based subclonal detection with clone frequency, driver mutation, proliferation index, and antigen escape risk assessment",
            },
            {
                "id": "spatial_niche",
                "name": "Spatial Niche Mapping",
                "description": "Spatial transcriptomics niche identification with cell-cell proximity, spatial autocorrelation, and tissue architecture analysis",
            },
            {
                "id": "trajectory_inference",
                "name": "Trajectory Inference",
                "description": "Pseudotime trajectory analysis with branch point detection, driver gene identification, and fate probability estimation",
            },
            {
                "id": "ligand_receptor",
                "name": "Ligand-Receptor Interaction",
                "description": "Cell-cell communication network analysis using CellPhoneDB/NicheNet with pathway enrichment and therapeutic target prioritization",
            },
            {
                "id": "biomarker_discovery",
                "name": "Biomarker Discovery",
                "description": "Differential expression-based biomarker identification with cell-type specificity scoring and clinical correlation assessment",
            },
            {
                "id": "cart_validation",
                "name": "CAR-T Target Validation",
                "description": "Single-cell-level CAR-T target assessment with on-tumor/off-tumor expression, TME compatibility, and antigen escape risk analysis",
            },
            {
                "id": "treatment_monitoring",
                "name": "Treatment Monitoring",
                "description": "Longitudinal single-cell tracking of treatment response including resistance clone emergence and immune repertoire shifts",
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
            lines.append(f"# HELP sc_agent_{key} Single-cell agent {key}")
            lines.append(f"# TYPE sc_agent_{key} counter")
            lines.append(f"sc_agent_{key} {val}")
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
        content={"error": exc.detail, "agent": "single-cell-intelligence-agent"},
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    logger.error(f"Unhandled exception: {exc}")
    with _metrics_lock:
        _metrics["errors_total"] += 1
    return JSONResponse(
        status_code=500,
        content={"error": "Internal server error", "agent": "single-cell-intelligence-agent"},
    )
