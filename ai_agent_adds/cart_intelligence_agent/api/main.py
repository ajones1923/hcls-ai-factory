"""CAR-T Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine as a production-ready REST API with
CORS, health checks, Prometheus-compatible metrics, and Pydantic request /
response schemas.

Endpoints:
    GET  /health          -- Service health with collection and vector counts
    GET  /collections     -- Collection names and record counts
    POST /query           -- Full RAG query (retrieve + LLM synthesis)
    POST /search          -- Evidence-only retrieval (no LLM, fast)
    POST /find-related    -- Cross-collection entity linking
    GET  /knowledge/stats -- Knowledge graph statistics
    GET  /metrics         -- Prometheus-compatible metrics (placeholder)

Port: 8522 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8522 --reload

Author: Adam Jones
Date: February 2026
"""

import logging
import os
import sys
import time
from collections import defaultdict
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, PlainTextResponse
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

# =====================================================================
# Path setup -- ensure project root is importable
# =====================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Load API key from rag-chat-pipeline .env if not already set
if not os.environ.get("ANTHROPIC_API_KEY"):
    _env_path = Path(os.environ.get("CART_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")) / ".env"
    if _env_path.exists():
        for _line in _env_path.read_text().splitlines():
            if _line.startswith("ANTHROPIC_API_KEY="):
                os.environ["ANTHROPIC_API_KEY"] = _line.split("=", 1)[1].strip().strip('"')
                break

from config.settings import settings
from src.collections import CARTCollectionManager
from src.knowledge import get_knowledge_stats
from src.models import AgentQuery, CrossCollectionResult, SearchHit
from src.rag_engine import CARTRAGEngine

# Route modules (meta-agent, reports, events)
from api.routes.meta_agent import router as meta_agent_router
from api.routes.reports import router as reports_router
from api.routes.events import router as events_router

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine: Optional[CARTRAGEngine] = None
_manager: Optional[CARTCollectionManager] = None

# Simple request counters for /metrics
_metrics: Dict[str, int] = {
    "requests_total": 0,
    "query_requests_total": 0,
    "search_requests_total": 0,
    "find_related_requests_total": 0,
    "errors_total": 0,
}


# =====================================================================
# Lifespan -- initialize engine on startup, disconnect on shutdown
# =====================================================================

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize the RAG engine and Milvus connection on startup."""
    global _engine, _manager

    # ── Collection manager ──
    _manager = CARTCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    _manager.connect()

    # ── Embedder ──
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

    # ── LLM client ──
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

    # ── Knowledge + query expansion modules ──
    from src import knowledge as kg
    from src import query_expansion as qe

    # ── Build engine ──
    _engine = CARTRAGEngine(
        collection_manager=_manager,
        embedder=embedder,
        llm_client=llm_client,
        knowledge=kg,
        query_expander=qe,
    )

    yield

    # ── Shutdown ──
    if _manager:
        _manager.disconnect()


# =====================================================================
# FastAPI app
# =====================================================================

app = FastAPI(
    title="CAR-T Intelligence Agent API",
    description=(
        "REST API for the CAR-T Intelligence Agent -- multi-collection RAG "
        "engine spanning literature, clinical trials, constructs, assays, "
        "manufacturing, safety, biomarkers, regulatory, sequences, and "
        "real-world evidence."
    ),
    version="1.0.0",
    docs_url="/docs",
    openapi_url="/openapi.json",
    lifespan=lifespan,
)

# ── CORS middleware ──
_cors_origins = [o.strip() for o in settings.CORS_ORIGINS.split(",") if o.strip()]
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
    api_key = getattr(settings, 'API_KEY', '') or ''
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


# ── Request size limit middleware ──
@app.middleware("http")
async def _limit_request_size(request: Request, call_next):
    """Reject request bodies that exceed the configured size limit."""
    content_length = request.headers.get("content-length")
    max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    if content_length and int(content_length) > max_bytes:
        return JSONResponse(status_code=413, content={"detail": "Request body too large"})
    return await call_next(request)

# ── Include route modules ──
app.include_router(meta_agent_router)
app.include_router(reports_router)
app.include_router(events_router)


# =====================================================================
# Root endpoint
# =====================================================================

@app.get("/", include_in_schema=False)
def root():
    return {"service": "CAR-T Intelligence Agent", "docs": "/docs", "health": "/health"}


# =====================================================================
# Pydantic request / response schemas
# =====================================================================

class HealthResponse(BaseModel):
    """Response schema for GET /health."""
    status: str = "healthy"
    collections: int = Field(..., description="Number of active collections")
    total_vectors: int = Field(..., description="Total vectors across all collections")


class CollectionInfo(BaseModel):
    """Single collection metadata."""
    name: str
    record_count: int


class CollectionsResponse(BaseModel):
    """Response schema for GET /collections."""
    collections: List[CollectionInfo]
    total: int


class QueryRequest(BaseModel):
    """Request schema for POST /query and POST /search."""
    question: str = Field(..., min_length=1, description="Natural-language question")
    target_antigen: Optional[str] = Field(None, max_length=100, description="Filter by target antigen (e.g. CD19, BCMA)")
    collections: Optional[List[str]] = Field(None, description="Restrict search to specific collections")
    year_min: Optional[int] = Field(None, ge=1990, le=2030, description="Minimum publication year")
    year_max: Optional[int] = Field(None, ge=1990, le=2030, description="Maximum publication year")


class EvidenceItem(BaseModel):
    """A single piece of evidence returned to the client."""
    collection: str
    id: str
    score: float
    text: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


class QueryResponse(BaseModel):
    """Response schema for POST /query (RAG with LLM)."""
    question: str
    answer: str
    evidence: List[EvidenceItem]
    knowledge_context: str = ""
    collections_searched: int = 0
    search_time_ms: float = 0.0


class SearchResponse(BaseModel):
    """Response schema for POST /search (evidence only, no LLM)."""
    question: str
    evidence: List[EvidenceItem]
    knowledge_context: str = ""
    collections_searched: int = 0
    search_time_ms: float = 0.0


class FindRelatedRequest(BaseModel):
    """Request schema for POST /find-related."""
    entity: str = Field(..., min_length=1, description="Entity name (product, antigen, trial ID, etc.)")
    top_k: int = Field(5, ge=1, le=50, description="Max results per collection")


class FindRelatedResponse(BaseModel):
    """Response schema for POST /find-related."""
    entity: str
    results: Dict[str, List[EvidenceItem]]
    total_hits: int


class KnowledgeStatsResponse(BaseModel):
    """Response schema for GET /knowledge/stats."""
    target_antigens: int
    targets_with_approved_products: int
    toxicity_profiles: int
    manufacturing_processes: int
    biomarkers: int
    regulatory_products: int
    immunogenicity_topics: int


# =====================================================================
# Helper -- convert internal SearchHit to API EvidenceItem
# =====================================================================

def _hit_to_evidence(hit: SearchHit) -> EvidenceItem:
    """Convert an internal SearchHit to the API EvidenceItem schema."""
    return EvidenceItem(
        collection=hit.collection,
        id=hit.id,
        score=hit.score,
        text=hit.text,
        metadata=hit.metadata,
    )


# =====================================================================
# Endpoints
# =====================================================================

@app.get("/health", response_model=HealthResponse, tags=["status"])
async def health():
    """Return service health with collection count and total vector count.

    Returns 503 if the engine or Milvus connection is unavailable.
    """
    _metrics["requests_total"] += 1

    if not _manager:
        raise HTTPException(status_code=503, detail="Engine not initialized")

    try:
        stats = _manager.get_collection_stats()
        total_collections = sum(1 for v in stats.values() if v > 0)
        total_vectors = sum(stats.values())
        return HealthResponse(
            status="healthy",
            collections=total_collections,
            total_vectors=total_vectors,
        )
    except Exception as e:
        _metrics["errors_total"] += 1
        logger.error(f"Milvus unavailable: {e}")
        raise HTTPException(status_code=503, detail="Service temporarily unavailable")


@app.get("/collections", response_model=CollectionsResponse, tags=["status"])
async def list_collections():
    """Return all collection names and their record counts."""
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
        _metrics["errors_total"] += 1
        logger.error(f"Failed to fetch collection stats: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/query", response_model=QueryResponse, tags=["rag"])
async def query(request: QueryRequest):
    """Full RAG query: retrieve evidence from Milvus, augment with the
    knowledge graph, and synthesize an LLM response.

    Requires both the embedding model and LLM client to be available.
    """
    _metrics["requests_total"] += 1
    _metrics["query_requests_total"] += 1

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")
    if not _engine.llm:
        raise HTTPException(status_code=503, detail="LLM client not available")
    if not _engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        agent_query = AgentQuery(
            question=request.question,
            target_antigen=request.target_antigen,
        )

        # Retrieve evidence
        evidence: CrossCollectionResult = _engine.retrieve(
            query=agent_query,
            collections_filter=request.collections,
            year_min=request.year_min,
            year_max=request.year_max,
        )

        # Generate LLM response
        from src.rag_engine import CART_SYSTEM_PROMPT
        prompt_text = _engine._build_prompt(request.question, evidence)
        answer = _engine.llm.generate(
            prompt=prompt_text,
            system_prompt=CART_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return QueryResponse(
            question=request.question,
            answer=answer,
            evidence=[_hit_to_evidence(h) for h in evidence.hits],
            knowledge_context=evidence.knowledge_context,
            collections_searched=evidence.total_collections_searched,
            search_time_ms=evidence.search_time_ms,
        )

    except HTTPException:
        raise
    except Exception as e:
        _metrics["errors_total"] += 1
        logger.error(f"Query failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/search", response_model=SearchResponse, tags=["rag"])
async def search(request: QueryRequest):
    """Evidence-only retrieval (no LLM). Useful for fast retrieval when
    only evidence snippets are needed without synthesis.
    """
    _metrics["requests_total"] += 1
    _metrics["search_requests_total"] += 1

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")
    if not _engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        agent_query = AgentQuery(
            question=request.question,
            target_antigen=request.target_antigen,
        )

        evidence: CrossCollectionResult = _engine.retrieve(
            query=agent_query,
            collections_filter=request.collections,
            year_min=request.year_min,
            year_max=request.year_max,
        )

        return SearchResponse(
            question=request.question,
            evidence=[_hit_to_evidence(h) for h in evidence.hits],
            knowledge_context=evidence.knowledge_context,
            collections_searched=evidence.total_collections_searched,
            search_time_ms=evidence.search_time_ms,
        )

    except HTTPException:
        raise
    except Exception as e:
        _metrics["errors_total"] += 1
        logger.error(f"Search failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/find-related", response_model=FindRelatedResponse, tags=["rag"])
async def find_related(request: FindRelatedRequest):
    """Find all evidence related to an entity across all 10 collections.

    Enables queries like "show me everything about Yescarta" spanning
    literature, trials, constructs, safety, regulatory, and more.
    """
    _metrics["requests_total"] += 1
    _metrics["find_related_requests_total"] += 1

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")
    if not _engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        raw_results: Dict[str, List[SearchHit]] = _engine.find_related(
            entity=request.entity,
            top_k=request.top_k,
        )

        # Convert to API schema
        api_results: Dict[str, List[EvidenceItem]] = {}
        total = 0
        for coll_name, hits in raw_results.items():
            api_results[coll_name] = [_hit_to_evidence(h) for h in hits]
            total += len(hits)

        return FindRelatedResponse(
            entity=request.entity,
            results=api_results,
            total_hits=total,
        )

    except HTTPException:
        raise
    except Exception as e:
        _metrics["errors_total"] += 1
        logger.error(f"Find-related failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.get("/knowledge/stats", response_model=KnowledgeStatsResponse, tags=["knowledge"])
async def knowledge_stats():
    """Return statistics about the CAR-T knowledge graph.

    Includes counts of target antigens, approved products, toxicity
    profiles, manufacturing processes, biomarkers, and regulatory entries.
    """
    _metrics["requests_total"] += 1

    try:
        stats = get_knowledge_stats()
        return KnowledgeStatsResponse(**stats)
    except Exception as e:
        _metrics["errors_total"] += 1
        logger.error(f"Knowledge stats failed: {e}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.get("/metrics", response_class=PlainTextResponse, tags=["monitoring"])
async def metrics():
    """Prometheus-compatible metrics endpoint (placeholder).

    Returns basic request counters in Prometheus exposition format.
    A full implementation would integrate with prometheus_client.
    """
    lines = [
        "# HELP cart_api_requests_total Total API requests",
        "# TYPE cart_api_requests_total counter",
        f'cart_api_requests_total {_metrics["requests_total"]}',
        "",
        "# HELP cart_api_query_requests_total Total /query requests",
        "# TYPE cart_api_query_requests_total counter",
        f'cart_api_query_requests_total {_metrics["query_requests_total"]}',
        "",
        "# HELP cart_api_search_requests_total Total /search requests",
        "# TYPE cart_api_search_requests_total counter",
        f'cart_api_search_requests_total {_metrics["search_requests_total"]}',
        "",
        "# HELP cart_api_find_related_requests_total Total /find-related requests",
        "# TYPE cart_api_find_related_requests_total counter",
        f'cart_api_find_related_requests_total {_metrics["find_related_requests_total"]}',
        "",
        "# HELP cart_api_errors_total Total error responses",
        "# TYPE cart_api_errors_total counter",
        f'cart_api_errors_total {_metrics["errors_total"]}',
        "",
    ]

    # Add collection vector counts if available
    if _manager:
        try:
            stats = _manager.get_collection_stats()
            lines.append("# HELP cart_collection_vectors Number of vectors per collection")
            lines.append("# TYPE cart_collection_vectors gauge")
            for name, count in stats.items():
                lines.append(f'cart_collection_vectors{{collection="{name}"}} {count}')
            lines.append("")
        except Exception:
            pass

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
