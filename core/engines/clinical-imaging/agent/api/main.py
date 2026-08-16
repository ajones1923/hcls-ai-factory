"""FastAPI server for Imaging Intelligence Agent.

Multi-collection RAG engine API with NIM integration,
workflow execution, and Prometheus metrics.
Runs on port 8524.

Author: Adam Jones
Date: February 2026
"""

import time
from collections import defaultdict
from contextlib import asynccontextmanager
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from loguru import logger
from pydantic import BaseModel, Field
from prometheus_client import (
    CONTENT_TYPE_LATEST,
    Counter,
    Histogram,
    generate_latest,
)
from starlette.responses import Response

from api.routes.meta_agent import router as meta_agent_router
from api.routes.nim import router as nim_router
from api.routes.reports import router as reports_router
from api.routes.workflows import router as workflows_router
from api.routes.events import events_router
from api.routes.preview import router as preview_router
from api.routes.demo_cases import router as demo_cases_router
from api.routes.protocol import router as protocol_router
from api.routes.dose import router as dose_router
from api.routes.analytics import router as analytics_router
from api.routes.live_analysis import router as live_analysis_router

# =====================================================================
# Prometheus Metrics
# =====================================================================

QUERY_COUNT = Counter("imaging_agent_queries_total", "Total RAG queries", ["endpoint"])
QUERY_LATENCY = Histogram(
    "imaging_agent_query_duration_seconds",
    "Query latency in seconds",
    ["endpoint"],
    buckets=[0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0],
)
SEARCH_HITS = Histogram(
    "imaging_agent_search_hits",
    "Number of evidence hits per query",
    buckets=[0, 5, 10, 20, 50, 100],
)

# =====================================================================
# Request / Response Models
# =====================================================================


class QueryRequest(BaseModel):
    """RAG query request."""
    question: str = Field(..., min_length=3, max_length=2000)
    modality: Optional[str] = Field(None, description="Filter by imaging modality")
    body_region: Optional[str] = Field(None, description="Filter by body region")
    top_k: int = Field(5, ge=1, le=50, description="Results per collection")
    include_genomic: bool = Field(True, description="Include genomic_evidence collection")
    include_nim: bool = Field(True, description="Allow NIM service invocation")
    collections: Optional[List[str]] = Field(None, description="Specific collections to search")
    year_min: Optional[int] = Field(None, ge=1990, le=2030)
    year_max: Optional[int] = Field(None, ge=1990, le=2030)
    conversation_context: str = Field("", description="Prior conversation context")


class QueryResponse(BaseModel):
    """RAG query response."""
    question: str
    answer: str
    evidence_count: int
    collections_searched: int
    search_time_ms: float
    nim_services_used: List[str] = []


class SearchRequest(BaseModel):
    """Evidence-only search request (no LLM synthesis)."""
    question: str = Field(..., min_length=3, max_length=2000)
    modality: Optional[str] = None
    body_region: Optional[str] = None
    top_k: int = Field(5, ge=1, le=50)
    collections: Optional[List[str]] = None
    year_min: Optional[int] = None
    year_max: Optional[int] = None


class SearchHitResponse(BaseModel):
    """Single search hit in response."""
    collection: str
    id: str
    score: float
    text: str
    metadata: Dict[str, Any] = {}


class SearchResponse(BaseModel):
    """Evidence-only search response."""
    query: str
    hits: List[SearchHitResponse]
    total_hits: int
    collections_searched: int
    search_time_ms: float
    knowledge_context: str = ""


class FindRelatedRequest(BaseModel):
    """Cross-collection entity linking request."""
    entity: str = Field(..., min_length=2, max_length=500)
    top_k: int = Field(5, ge=1, le=20)


class FindRelatedResponse(BaseModel):
    """Cross-collection entity linking response."""
    entity: str
    collections: Dict[str, List[SearchHitResponse]]
    total_hits: int


class CollectionInfo(BaseModel):
    """Information about a single collection."""
    name: str
    count: int
    label: str


class HealthResponse(BaseModel):
    """Service health response."""
    status: str
    collections: Dict[str, int]
    total_vectors: int
    nim_services: Dict[str, str]


# =====================================================================
# Application State
# =====================================================================

_state: Dict[str, Any] = {}


# =====================================================================
# Lifespan
# =====================================================================


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize and tear down application resources."""
    logger.info("Starting Imaging Intelligence Agent API on port 8524")

    try:
        from config.settings import settings
        from sentence_transformers import SentenceTransformer
        from src.collections import ImagingCollectionManager
        from src.nim.service_manager import NIMServiceManager
        from src.rag_engine import ImagingRAGEngine

        # Collection manager
        manager = ImagingCollectionManager(
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        # Connect and create-collections are logged separately: folding them into one
        # try/except reported an AttributeError here (the method is create_all_collections,
        # not ensure_collections) as "Milvus connection deferred", which reads as a network
        # fault on a connection that had in fact already succeeded.
        try:
            manager.connect()
            logger.info("Milvus connected")
        except Exception as e:
            logger.warning(f"Milvus connection deferred: {e}")
        else:
            try:
                manager.create_all_collections()
                logger.info("Milvus collections ensured")
            except Exception as e:
                logger.warning(f"Milvus collection setup failed: {e}")

        # Embedding model
        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
        logger.info(f"Loaded embedding model: {settings.EMBEDDING_MODEL}")

        # NIM service manager
        # NVIDIA_API_KEY can come from settings (IMAGING_NVIDIA_API_KEY env)
        # or from the bare NVIDIA_API_KEY env var
        import os
        nvidia_api_key = getattr(settings, "NVIDIA_API_KEY", None) or os.environ.get("NVIDIA_API_KEY")
        nim_manager = NIMServiceManager(settings, nvidia_api_key=nvidia_api_key)

        # RAG engine
        engine = ImagingRAGEngine(
            collection_manager=manager,
            embedder=embedder,
            llm_client=nim_manager.llm,
            nim_service_manager=nim_manager,
        )

        _state["manager"] = manager
        _state["embedder"] = embedder
        _state["nim_manager"] = nim_manager
        _state["engine"] = engine
        _state["settings"] = settings

        # Cross-modal trigger (imaging -> genomics)
        from src.cross_modal import CrossModalTrigger

        cross_modal = CrossModalTrigger(
            collection_manager=manager,
            embedder=embedder,
            enabled=settings.CROSS_MODAL_ENABLED,
        )
        _state["cross_modal_trigger"] = cross_modal
        logger.info(
            f"Cross-modal trigger initialized (enabled={settings.CROSS_MODAL_ENABLED})"
        )

        logger.info("All components initialized successfully")

        # Start ingest scheduler if enabled
        if settings.INGEST_ENABLED:
            try:
                from src.scheduler import ImagingIngestScheduler

                ingest_scheduler = ImagingIngestScheduler(
                    collection_manager=manager,
                    embedder=embedder,
                )
                ingest_scheduler.start(
                    pubmed_interval_hours=settings.INGEST_SCHEDULE_HOURS,
                    trials_interval_hours=settings.INGEST_TRIALS_SCHEDULE_HOURS,
                )
                _state["ingest_scheduler"] = ingest_scheduler
                logger.info(
                    f"Ingest scheduler started — "
                    f"PubMed every {settings.INGEST_SCHEDULE_HOURS}h, "
                    f"Trials every {settings.INGEST_TRIALS_SCHEDULE_HOURS}h"
                )
            except Exception as e:
                logger.warning(f"Ingest scheduler failed to start: {e}")
        else:
            logger.info("Ingest scheduler disabled (INGEST_ENABLED=False)")

    except Exception as e:
        logger.error(f"Initialization failed: {e}")
        _state["error"] = str(e)

    yield

    # Stop ingest scheduler on shutdown
    ingest_scheduler = _state.get("ingest_scheduler")
    if ingest_scheduler and ingest_scheduler.is_running:
        ingest_scheduler.stop()
        logger.info("Ingest scheduler stopped")

    logger.info("Shutting down Imaging Intelligence Agent API")
    _state.clear()


# =====================================================================
# FastAPI App
# =====================================================================

app = FastAPI(
    title="Imaging Intelligence Agent API",
    description=(
        "Multi-collection RAG engine for medical imaging intelligence. "
        "Searches 10 imaging-specific Milvus collections plus genomic evidence, "
        "augmented with NVIDIA NIM services (VISTA-3D, MAISI, VILA-M3)."
    ),
    version="1.0.0",
    docs_url="/docs",
    openapi_url="/openapi.json",
    lifespan=lifespan,
)

# -- Platform security gate (audit 2026-08-15) -------------------------------
# This service returned clinical decision-support output with NO authentication. The gate is a
# no-op until HCLS_API_KEY (or the per-service HCLS_API_KEY_<SLUG>) is set, so the default
# trusted-network posture is unchanged; once set it FAILS CLOSED on every route except /health
# and /docs. install_governance adds a request id, timing, and X-HCLS-Governed listing only the
# gates a handler actually invoked.
try:  # pragma: no cover - platform lib is optional at import time
    from hcls_common.api_auth import install_api_key_auth
    from hcls_common.api_gate import install_governance
    install_api_key_auth(app, service="imaging")
    install_governance(app, service="imaging", capability_id="imaging-intelligence-agent")
except Exception:  # never take a running service down over the gate
    pass


from config.settings import settings as _settings

_cors_origins = [o.strip() for o in _settings.CORS_ORIGINS.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── Auth middleware (optional, based on API_KEY setting) ──
_AUTH_SKIP_PATHS = {"/health", "/healthz", "/metrics", "/docs", "/openapi.json"}


@app.middleware("http")
async def check_api_key(request: Request, call_next):
    if request.url.path in _AUTH_SKIP_PATHS:
        return await call_next(request)
    api_key = getattr(_settings, 'API_KEY', '') or ''
    if not api_key:
        return await call_next(request)
    provided = request.headers.get("X-API-Key", "")
    if provided != api_key:
        return JSONResponse(status_code=401, content={"detail": "Invalid or missing API key"})
    return await call_next(request)


# ── Rate limiting middleware ──
_rate_limit_store: dict = defaultdict(list)
_RATE_LIMIT_MAX = 100
_RATE_LIMIT_WINDOW = 60


@app.middleware("http")
async def rate_limit(request: Request, call_next):
    if request.url.path in {"/health", "/healthz", "/metrics", "/docs", "/openapi.json"}:
        return await call_next(request)
    client_ip = request.client.host if request.client else "unknown"
    now = time.time()
    _rate_limit_store[client_ip] = [t for t in _rate_limit_store.get(client_ip, []) if now - t < _RATE_LIMIT_WINDOW]
    if len(_rate_limit_store[client_ip]) >= _RATE_LIMIT_MAX:
        return JSONResponse(status_code=429, content={"detail": "Rate limit exceeded"})
    _rate_limit_store[client_ip].append(now)
    return await call_next(request)


@app.middleware("http")
async def _limit_request_size(request: Request, call_next):
    """Reject request bodies that exceed the configured size limit."""
    content_length = request.headers.get("content-length")
    max_bytes = _settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    if content_length and int(content_length) > max_bytes:
        return JSONResponse(status_code=413, content={"detail": "Request body too large"})
    return await call_next(request)

# Include routers
app.include_router(meta_agent_router, prefix="/api", tags=["Meta-Agent"])
app.include_router(nim_router, prefix="/nim", tags=["NIM Services"])
app.include_router(workflows_router, tags=["Workflows"])
app.include_router(reports_router, tags=["Reports"])
app.include_router(events_router, tags=["Events (Phase 2)"])
app.include_router(preview_router, tags=["Preview"])
app.include_router(demo_cases_router, tags=["Demo Cases"])
app.include_router(protocol_router, tags=["Protocol Optimizer"])
app.include_router(dose_router, tags=["Dose Intelligence"])
app.include_router(analytics_router, tags=["Analytics"])
app.include_router(live_analysis_router, tags=["Live Analysis"])

# =====================================================================
# Static Files — serve sample images for UI display
# =====================================================================

import os
_sample_images_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "sample_images")
if os.path.isdir(_sample_images_dir):
    app.mount("/images", StaticFiles(directory=_sample_images_dir), name="sample_images")
    logger.info(f"Serving sample images from {_sample_images_dir}")

_segmentation_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "demo", "segmentation")
if os.path.isdir(_segmentation_dir):
    app.mount("/segmentation", StaticFiles(directory=_segmentation_dir), name="segmentation")
    logger.info(f"Serving segmentation overlays from {_segmentation_dir}")

# Coronary artery meshes rendered from CoronariesNC6 (real vessel geometry with manual rater
# ground truth). These replace the cardiac panels that previously showed abdominal and spine
# anatomy captioned as coronary — see scripts/render_coronary_mesh.py.
_coronary_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "demo", "coronary")
if os.path.isdir(_coronary_dir):
    app.mount("/coronary", StaticFiles(directory=_coronary_dir), name="coronary")
    logger.info(f"Serving coronary mesh renders from {_coronary_dir}")


# =====================================================================
# Root endpoint
# =====================================================================

@app.get("/", include_in_schema=False)
def root():
    return {"service": "Imaging Intelligence Agent", "docs": "/docs", "health": "/health"}


# =====================================================================
# Helper
# =====================================================================


def get_engine():
    """Get the RAG engine from application state."""
    engine = _state.get("engine")
    if engine is None:
        raise HTTPException(status_code=503, detail="Engine not initialized. Check /health for details.")
    return engine


def get_manager():
    """Get the collection manager from application state."""
    manager = _state.get("manager")
    if manager is None:
        raise HTTPException(status_code=503, detail="Collection manager not initialized.")
    return manager


def get_nim_manager():
    """Get the NIM service manager from application state."""
    nim_manager = _state.get("nim_manager")
    if nim_manager is None:
        raise HTTPException(status_code=503, detail="NIM service manager not initialized.")
    return nim_manager


# =====================================================================
# Core Endpoints
# =====================================================================


@app.get("/health", response_model=HealthResponse, tags=["Health"])
async def health_check():
    """Service health check with collection stats and NIM status."""
    manager = _state.get("manager")
    nim_manager = _state.get("nim_manager")

    collections = {}
    total_vectors = 0
    try:
        if manager:
            collections = manager.get_collection_stats()
            total_vectors = sum(collections.values())
    except Exception as e:
        logger.warning(f"Failed to get collection stats: {e}")

    nim_services = {}
    try:
        if nim_manager:
            nim_services = nim_manager.check_all_services()
    except Exception as e:
        logger.warning(f"Failed to check NIM services: {e}")

    error = _state.get("error")
    status = "degraded" if error else ("healthy" if manager else "initializing")

    return HealthResponse(
        status=status,
        collections=collections,
        total_vectors=total_vectors,
        nim_services=nim_services,
    )


@app.get("/collections", response_model=List[CollectionInfo], tags=["Collections"])
async def list_collections():
    """List all collections with their record counts."""
    from src.rag_engine import COLLECTION_CONFIG

    manager = get_manager()

    try:
        stats = manager.get_collection_stats()
    except Exception as e:
        logger.error(f"Cannot fetch collection stats: {e}")
        raise HTTPException(status_code=503, detail="Service temporarily unavailable")

    result = []
    for coll_name, config in COLLECTION_CONFIG.items():
        result.append(CollectionInfo(
            name=coll_name,
            count=stats.get(coll_name, 0),
            label=config.get("label", coll_name),
        ))
    return result


@app.post("/query", response_model=QueryResponse, tags=["RAG"])
async def rag_query(request: QueryRequest):
    """Full RAG query: multi-collection retrieval + LLM synthesis."""
    start = time.time()
    QUERY_COUNT.labels(endpoint="query").inc()

    engine = get_engine()

    kwargs = {"top_k_per_collection": request.top_k}
    if request.collections:
        kwargs["collections_filter"] = request.collections
    if request.modality:
        kwargs["modality_filter"] = request.modality
    if request.body_region:
        kwargs["body_region_filter"] = request.body_region
    if request.year_min:
        kwargs["year_min"] = request.year_min
    if request.year_max:
        kwargs["year_max"] = request.year_max

    try:
        # Retrieve evidence
        evidence = engine.retrieve(request.question, **kwargs)
        SEARCH_HITS.observe(evidence.hit_count)

        # Synthesize answer
        answer = engine.query(
            request.question,
            conversation_context=request.conversation_context,
            **kwargs,
        )

        elapsed = time.time() - start
        QUERY_LATENCY.labels(endpoint="query").observe(elapsed)

        # Determine which NIM services were used
        nim_used = []
        nim_manager = _state.get("nim_manager")
        if nim_manager and request.include_nim:
            status = nim_manager.check_all_services()
            nim_used = [name for name, s in status.items() if s in ("available", "cloud", "anthropic", "mock")]

        return QueryResponse(
            question=request.question,
            answer=answer,
            evidence_count=evidence.hit_count,
            collections_searched=evidence.total_collections_searched,
            search_time_ms=evidence.search_time_ms,
            nim_services_used=nim_used,
        )
    except Exception as e:
        logger.error(f"Query failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/search", response_model=SearchResponse, tags=["RAG"])
async def search_evidence(request: SearchRequest):
    """Evidence-only search across collections (no LLM synthesis)."""
    start = time.time()
    QUERY_COUNT.labels(endpoint="search").inc()

    engine = get_engine()

    kwargs = {"top_k_per_collection": request.top_k}
    if request.collections:
        kwargs["collections_filter"] = request.collections
    if request.modality:
        kwargs["modality_filter"] = request.modality
    if request.body_region:
        kwargs["body_region_filter"] = request.body_region
    if request.year_min:
        kwargs["year_min"] = request.year_min
    if request.year_max:
        kwargs["year_max"] = request.year_max

    try:
        evidence = engine.retrieve(request.question, **kwargs)
        SEARCH_HITS.observe(evidence.hit_count)

        elapsed = time.time() - start
        QUERY_LATENCY.labels(endpoint="search").observe(elapsed)

        hits = [
            SearchHitResponse(
                collection=h.collection,
                id=h.id,
                score=h.score,
                text=h.text,
                metadata=h.metadata,
            )
            for h in evidence.hits
        ]

        return SearchResponse(
            query=request.question,
            hits=hits,
            total_hits=evidence.hit_count,
            collections_searched=evidence.total_collections_searched,
            search_time_ms=evidence.search_time_ms,
            knowledge_context=evidence.knowledge_context,
        )
    except Exception as e:
        logger.error(f"Search failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/find-related", response_model=FindRelatedResponse, tags=["RAG"])
async def find_related_entities(request: FindRelatedRequest):
    """Cross-collection entity linking -- find related evidence for an entity."""
    QUERY_COUNT.labels(endpoint="find_related").inc()

    engine = get_engine()

    try:
        related = engine.find_related(request.entity, top_k=request.top_k)

        collections_response = {}
        total = 0
        for coll, hits in related.items():
            collections_response[coll] = [
                SearchHitResponse(
                    collection=h.collection,
                    id=h.id,
                    score=h.score,
                    text=h.text,
                    metadata=h.metadata,
                )
                for h in hits
            ]
            total += len(hits)

        return FindRelatedResponse(
            entity=request.entity,
            collections=collections_response,
            total_hits=total,
        )
    except Exception as e:
        logger.error(f"Find-related failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.get("/knowledge/stats", tags=["Knowledge"])
async def knowledge_stats():
    """Return statistics about the imaging domain knowledge graph."""
    from src.knowledge import get_knowledge_stats

    return get_knowledge_stats()


@app.get("/metrics", tags=["Monitoring"])
async def prometheus_metrics():
    """Prometheus-compatible metrics endpoint."""
    return Response(
        content=generate_latest(),
        media_type=CONTENT_TYPE_LATEST,
    )


# =====================================================================
# Run
# =====================================================================

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "api.main:app",
        host="0.0.0.0",
        port=8524,
        reload=True,
        log_level="info",
    )
