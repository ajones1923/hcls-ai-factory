"""FastAPI server for Imaging Intelligence Agent.

Multi-collection RAG engine API with NIM integration,
workflow execution, and Prometheus metrics.
Runs on port 8524.

Author: Adam Jones
Date: February 2026
"""

import time
from contextlib import asynccontextmanager
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
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
        try:
            manager.connect()
            manager.ensure_collections()
            logger.info("Milvus connected and collections ensured")
        except Exception as e:
            logger.warning(f"Milvus connection deferred: {e}")

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
    except Exception as e:
        logger.error(f"Initialization failed: {e}")
        _state["error"] = str(e)

    yield

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
    lifespan=lifespan,
)

from config.settings import settings as _settings

_cors_origins = [o.strip() for o in _settings.CORS_ORIGINS.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


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
        raise HTTPException(status_code=503, detail=f"Cannot fetch collection stats: {e}")

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
        raise HTTPException(status_code=500, detail=f"Query processing failed: {e}")


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
        raise HTTPException(status_code=500, detail=f"Search failed: {e}")


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
        raise HTTPException(status_code=500, detail=f"Entity linking failed: {e}")


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
