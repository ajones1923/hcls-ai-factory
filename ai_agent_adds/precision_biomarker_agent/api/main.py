"""Biomarker Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine and analysis modules as a
production-ready REST API with CORS, health checks, Prometheus-compatible
metrics, and Pydantic request/response schemas.

Endpoints:
    GET  /health           -- Service health with collection and vector counts
    GET  /collections      -- Collection names and record counts
    GET  /knowledge/stats  -- Knowledge graph statistics
    GET  /metrics          -- Prometheus-compatible metrics (placeholder)

    Versioned routes (via api/routes/):
    POST /v1/analyze       -- Full patient analysis
    POST /v1/biological-age -- Biological age calculation
    POST /v1/disease-risk  -- Disease trajectory analysis
    POST /v1/pgx           -- Pharmacogenomic mapping
    POST /v1/query         -- RAG Q&A query
    GET  /v1/health        -- V1 health check

Port: 8529 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8529 --reload

Author: Adam Jones
Date: March 2026
"""

import os
import sys
import threading
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, PlainTextResponse
from loguru import logger
from pydantic import BaseModel, Field

# =====================================================================
# Path setup -- ensure project root is importable
# =====================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Load API key: prefer env vars, then rag-chat-pipeline .env file
_api_key = (
    os.environ.get("ANTHROPIC_API_KEY")
    or os.environ.get("BIOMARKER_ANTHROPIC_API_KEY")
)
if not _api_key:
    _env_path = PROJECT_ROOT.parent.parent / "rag-chat-pipeline" / ".env"
    if _env_path.exists():
        for _line in _env_path.read_text().splitlines():
            if _line.startswith("ANTHROPIC_API_KEY="):
                _api_key = _line.split("=", 1)[1].strip().strip('"')
                break
if _api_key:
    os.environ["ANTHROPIC_API_KEY"] = _api_key

# Route modules
from api.routes.analysis import router as analysis_router
from api.routes.events import router as events_router
from api.routes.reports import router as reports_router
from config.settings import settings
from src.biological_age import BiologicalAgeCalculator
from src.disease_trajectory import DiseaseTrajectoryAnalyzer
from src.genotype_adjustment import GenotypeAdjuster
from src.knowledge import get_knowledge_stats
from src.pharmacogenomics import PharmacogenomicMapper
from src.rag_engine import BiomarkerRAGEngine

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine: Optional[BiomarkerRAGEngine] = None
_agent = None  # PrecisionBiomarkerAgent
_manager = None  # Collection manager
_bio_age_calc: Optional[BiologicalAgeCalculator] = None
_trajectory_analyzer: Optional[DiseaseTrajectoryAnalyzer] = None
_pgx_mapper: Optional[PharmacogenomicMapper] = None
_genotype_adjuster: Optional[GenotypeAdjuster] = None

# Simple request counters for /metrics
_metrics: Dict[str, int] = {
    "requests_total": 0,
    "query_requests_total": 0,
    "search_requests_total": 0,
    "analyze_requests_total": 0,
    "bio_age_requests_total": 0,
    "errors_total": 0,
}
_metrics_lock = threading.Lock()


# =====================================================================
# Lifespan -- initialize engine on startup, disconnect on shutdown
# =====================================================================

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize the RAG engine, analysis modules, and Milvus on startup."""
    global _engine, _agent, _manager
    global _bio_age_calc, _trajectory_analyzer, _pgx_mapper, _genotype_adjuster

    # -- Collection manager --
    try:
        from src.collections import BiomarkerCollectionManager
        _manager = BiomarkerCollectionManager(
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        _manager.connect()
    except Exception:
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
    except ImportError:
        embedder = None

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
                msg = self.client.messages.create(
                    model=settings.LLM_MODEL,
                    max_tokens=max_tokens,
                    temperature=temperature,
                    system=system_prompt,
                    messages=[{"role": "user", "content": prompt}],
                )
                return msg.content[0].text

            def generate_stream(
                self, prompt: str, system_prompt: str = "",
                max_tokens: int = 2048, temperature: float = 0.7,
            ):
                with self.client.messages.stream(
                    model=settings.LLM_MODEL,
                    max_tokens=max_tokens,
                    temperature=temperature,
                    system=system_prompt,
                    messages=[{"role": "user", "content": prompt}],
                ) as stream:
                    for text in stream.text_stream:
                        yield text

        llm_client = _LLMClient()
    except (ImportError, Exception):
        llm_client = None

    # -- Knowledge module --
    from src import knowledge as kg

    # -- Analysis modules --
    _bio_age_calc = BiologicalAgeCalculator()
    _trajectory_analyzer = DiseaseTrajectoryAnalyzer()
    _pgx_mapper = PharmacogenomicMapper()
    _genotype_adjuster = GenotypeAdjuster()

    # -- Build RAG engine --
    _engine = BiomarkerRAGEngine(
        collection_manager=_manager,
        embedder=embedder,
        llm_client=llm_client,
        knowledge=kg,
    )

    # -- Build agent --
    from src.agent import PrecisionBiomarkerAgent
    _agent = PrecisionBiomarkerAgent(
        rag_engine=_engine,
        bio_age_calc=_bio_age_calc,
        trajectory_analyzer=_trajectory_analyzer,
        pgx_mapper=_pgx_mapper,
        genotype_adjuster=_genotype_adjuster,
    )

    # Store references on app.state for route access
    app.state.engine = _engine
    app.state.agent = _agent
    app.state.manager = _manager
    app.state.bio_age_calc = _bio_age_calc
    app.state.trajectory_analyzer = _trajectory_analyzer
    app.state.pgx_mapper = _pgx_mapper
    app.state.genotype_adjuster = _genotype_adjuster
    app.state.metrics = _metrics

    yield

    # -- Shutdown --
    if _manager:
        _manager.disconnect()


# =====================================================================
# FastAPI app
# =====================================================================

app = FastAPI(
    title="Biomarker Intelligence Agent API",
    description=(
        "Biological age, disease trajectories, pharmacogenomics, and "
        "RAG-powered evidence retrieval"
    ),
    version="1.0.0",
    docs_url="/docs",
    openapi_url="/openapi.json",
    lifespan=lifespan,
)

# -- CORS middleware --
_cors_origins = [o.strip() for o in settings.CORS_ORIGINS.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Content-Type", "Authorization", "X-API-Key", "Accept"],
)



# -- API key authentication middleware --
@app.middleware("http")
async def _check_api_key(request: Request, call_next):
    """Validate API key if API_KEY is configured in settings."""
    api_key = getattr(settings, "API_KEY", None)
    if api_key:
        # Skip auth for health and metrics only; docs require auth in production
        skip_paths = {"/health", "/metrics"}
        if request.url.path not in skip_paths:
            provided_key = request.headers.get("X-API-Key", "")
            if provided_key != api_key:
                return JSONResponse(status_code=401, content={"detail": "Invalid or missing API key"})
    return await call_next(request)


# -- Request size limit middleware --
@app.middleware("http")
async def _limit_request_size(request: Request, call_next):
    """Reject request bodies that exceed the configured size limit."""
    content_length = request.headers.get("content-length")
    max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    if content_length:
        try:
            if int(content_length) > max_bytes:
                return JSONResponse(status_code=413, content={"detail": "Request too large"})
        except ValueError:
            pass  # Malformed header, let request proceed
    return await call_next(request)


# -- Include route modules --
app.include_router(analysis_router)
app.include_router(reports_router)
app.include_router(events_router)


# =====================================================================
# Root endpoint
# =====================================================================

@app.get("/", include_in_schema=False)
def root():
    return {"service": "Biomarker Intelligence Agent", "docs": "/docs", "health": "/health"}


# =====================================================================
# Pydantic request / response schemas
# =====================================================================

class HealthResponse(BaseModel):
    """Response schema for GET /health."""
    status: str = "healthy"
    collections: int = Field(..., description="Number of active collections")
    total_vectors: int = Field(..., description="Total vectors across all collections")
    agent_ready: bool = Field(..., description="Whether the full agent is initialized")


class CollectionInfo(BaseModel):
    """Single collection metadata."""
    name: str
    record_count: int


class CollectionsResponse(BaseModel):
    """Response schema for GET /collections."""
    collections: List[CollectionInfo]
    total: int


class KnowledgeStatsResponse(BaseModel):
    """Response schema for GET /knowledge/stats."""
    disease_domains: int
    total_biomarkers: int
    total_genetic_modifiers: int
    pharmacogenes: int
    pgx_drug_interactions: int
    phenoage_markers: int
    cross_modal_links: int


# =====================================================================
# Endpoints
# =====================================================================

@app.get("/health", response_model=HealthResponse, tags=["status"])
async def health():
    """Return service health with collection count and total vector count."""
    with _metrics_lock:
        _metrics["requests_total"] += 1

    agent_ready = _engine is not None and _agent is not None

    if not _manager:
        return HealthResponse(
            status="degraded",
            collections=0,
            total_vectors=0,
            agent_ready=agent_ready,
        )

    try:
        stats = _manager.get_collection_stats()
        total_collections = sum(1 for v in stats.values() if v > 0)
        total_vectors = sum(stats.values())
        return HealthResponse(
            status="healthy",
            collections=total_collections,
            total_vectors=total_vectors,
            agent_ready=agent_ready,
        )
    except Exception as e:
        logger.exception(f"Health check failed -- Milvus unavailable: {e}")
        with _metrics_lock:
            _metrics["errors_total"] += 1
        raise HTTPException(status_code=503, detail="Service temporarily unavailable")


@app.get("/collections", response_model=CollectionsResponse, tags=["status"])
async def list_collections():
    """Return all collection names and their record counts."""
    with _metrics_lock:
        _metrics["requests_total"] += 1

    if not _manager:
        raise HTTPException(status_code=503, detail="Engine not initialized")

    try:
        stats = _manager.get_collection_stats()
        items = [
            CollectionInfo(name=name, record_count=count)
            for name, count in stats.items()
        ]
        return CollectionsResponse(
            collections=items,
            total=len(items),
        )
    except Exception as e:
        logger.exception(f"Failed to fetch collection stats: {e}")
        with _metrics_lock:
            _metrics["errors_total"] += 1
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.get("/knowledge/stats", response_model=KnowledgeStatsResponse, tags=["knowledge"])
async def knowledge_stats():
    """Return statistics about the biomarker knowledge graph."""
    with _metrics_lock:
        _metrics["requests_total"] += 1

    try:
        stats = get_knowledge_stats()
        return KnowledgeStatsResponse(**stats)
    except Exception as e:
        logger.exception(f"Knowledge stats failed: {e}")
        with _metrics_lock:
            _metrics["errors_total"] += 1
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.get("/metrics", response_class=PlainTextResponse, tags=["monitoring"])
async def metrics():
    """Prometheus-compatible metrics endpoint."""
    lines = [
        "# HELP biomarker_api_requests_total Total API requests",
        "# TYPE biomarker_api_requests_total counter",
        f'biomarker_api_requests_total {_metrics["requests_total"]}',
        "",
        "# HELP biomarker_api_query_requests_total Total /query requests",
        "# TYPE biomarker_api_query_requests_total counter",
        f'biomarker_api_query_requests_total {_metrics["query_requests_total"]}',
        "",
        "# HELP biomarker_api_search_requests_total Total /search requests",
        "# TYPE biomarker_api_search_requests_total counter",
        f'biomarker_api_search_requests_total {_metrics["search_requests_total"]}',
        "",
        "# HELP biomarker_api_analyze_requests_total Total /analyze requests",
        "# TYPE biomarker_api_analyze_requests_total counter",
        f'biomarker_api_analyze_requests_total {_metrics["analyze_requests_total"]}',
        "",
        "# HELP biomarker_api_bio_age_requests_total Total /biological-age requests",
        "# TYPE biomarker_api_bio_age_requests_total counter",
        f'biomarker_api_bio_age_requests_total {_metrics["bio_age_requests_total"]}',
        "",
        "# HELP biomarker_api_errors_total Total error responses",
        "# TYPE biomarker_api_errors_total counter",
        f'biomarker_api_errors_total {_metrics["errors_total"]}',
        "",
    ]

    # Add collection vector counts if available
    if _manager:
        try:
            stats = _manager.get_collection_stats()
            lines.append("# HELP biomarker_collection_vectors Number of vectors per collection")
            lines.append("# TYPE biomarker_collection_vectors gauge")
            for name, count in stats.items():
                lines.append(f'biomarker_collection_vectors{{collection="{name}"}} {count}')
            lines.append("")
        except Exception as e:
            logger.warning(f"Failed to fetch collection vectors for metrics: {e}")

    return "\n".join(lines) + "\n"


# =====================================================================
# Entrypoint for direct execution
# =====================================================================

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "api.main:app",
        host=settings.API_HOST,
        port=settings.API_PORT,
        reload=True,
    )
