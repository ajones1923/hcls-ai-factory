"""Pharmacogenomics Intelligence Agent -- FastAPI REST API.

Wraps the multi-collection RAG engine as a production-ready REST API with
CORS, health checks, Prometheus-compatible metrics, and Pydantic request /
response schemas.  Searches across 15 PGx collections (14 domain-specific +
genomic_evidence) for gene references, drug guidelines, HLA screening,
phenoconversion, dosing algorithms, clinical evidence, population data,
clinical trials, FDA labels, drug alternatives, patient profiles,
implementation protocols, and educational materials.

Endpoints:
    GET  /                 -- Service root
    GET  /health           -- Service health with collection and vector counts
    GET  /collections      -- Collection names and record counts
    POST /query            -- Full RAG query (retrieve + LLM synthesis)
    POST /search           -- Evidence-only retrieval (no LLM, fast)
    POST /find-related     -- Cross-collection entity linking
    GET  /knowledge/stats  -- Knowledge graph statistics
    GET  /metrics          -- Prometheus-compatible metrics

Port: 8107 (from config/settings.py)

Usage:
    uvicorn api.main:app --host 0.0.0.0 --port 8107 --reload

Author: Adam Jones
Date: March 2026
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
    _env_path = Path(os.environ.get("PGX_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")) / ".env"
    if _env_path.exists():
        for _line in _env_path.read_text().splitlines():
            if _line.startswith("ANTHROPIC_API_KEY="):
                os.environ["ANTHROPIC_API_KEY"] = _line.split("=", 1)[1].strip().strip('"')
                break

from config.settings import settings
from src.collections import PGxCollectionManager
from src.knowledge import get_knowledge_stats
from src.models import AgentQuery
from src.rag_engine import PGxRAGEngine

# Route modules (pgx_clinical, reports, events)
from api.routes.pgx_clinical import router as pgx_clinical_router
from api.routes.reports import router as reports_router
from api.routes.events import router as events_router

# =====================================================================
# Module-level state (populated during lifespan startup)
# =====================================================================

_engine: Optional[PGxRAGEngine] = None
_manager: Optional[PGxCollectionManager] = None

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

    # -- Collection manager --
    _manager = PGxCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    _manager.connect()

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

    # -- Knowledge + query expansion modules --
    from src import knowledge as kg
    from src import query_expansion as qe

    # -- Build engine --
    _engine = PGxRAGEngine(
        collection_manager=_manager,
        embedder=embedder,
        llm_client=llm_client,
        knowledge=kg,
        query_expander=qe,
    )

    yield

    # -- Shutdown --
    if _manager:
        _manager.disconnect()


# =====================================================================
# FastAPI app
# =====================================================================

app = FastAPI(
    title="Pharmacogenomics Intelligence Agent API",
    description=(
        "REST API for the Pharmacogenomics Intelligence Agent -- multi-collection "
        "RAG engine spanning gene references, drug guidelines, drug interactions, "
        "HLA hypersensitivity screening, phenoconversion, dosing algorithms, "
        "clinical evidence, population data, clinical trials, FDA labels, "
        "drug alternatives, patient profiles, implementation, and education."
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
    allow_methods=["*"],
    allow_headers=["*"],
)


# -- Auth middleware (optional, based on API_KEY setting) --
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


# -- Rate limiting middleware --
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


# -- Request size limit middleware --
@app.middleware("http")
async def _limit_request_size(request: Request, call_next):
    """Reject request bodies that exceed the configured size limit."""
    content_length = request.headers.get("content-length")
    max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    if content_length and int(content_length) > max_bytes:
        return JSONResponse(status_code=413, content={"detail": "Request body too large"})
    return await call_next(request)

# -- Include route modules --
app.include_router(pgx_clinical_router)
app.include_router(reports_router)
app.include_router(events_router)


# =====================================================================
# Root endpoint
# =====================================================================

@app.get("/", include_in_schema=False)
def root():
    return {"service": "Pharmacogenomics Intelligence Agent", "docs": "/docs", "health": "/health"}


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
    gene_filter: Optional[str] = Field(None, max_length=100, description="Filter by gene symbol (e.g. CYP2D6, CYP2C19)")
    drug_filter: Optional[str] = Field(None, max_length=200, description="Filter by drug name (e.g. codeine, warfarin)")
    collections: Optional[List[str]] = Field(None, description="Restrict search to specific collections")
    patient_id: Optional[str] = Field(None, max_length=100, description="Patient identifier for profile-based queries")
    medication_list: Optional[List[str]] = Field(None, description="Current medications for interaction checking")


class EvidenceItem(BaseModel):
    """A single piece of evidence returned to the client."""
    collection: str
    id: str
    score: float
    text: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


class AlertItem(BaseModel):
    """A clinical decision support alert."""
    alert_level: str
    gene: str = ""
    drug: str = ""
    phenotype: str = ""
    recommendation: str = ""
    evidence_pmids: str = ""


class QueryResponse(BaseModel):
    """Response schema for POST /query (RAG with LLM)."""
    question: str
    answer: str
    evidence: List[EvidenceItem]
    alerts: List[AlertItem] = Field(default_factory=list)
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
    entity: str = Field(..., min_length=1, description="Entity name (gene, drug, star allele, HLA allele, etc.)")
    top_k: int = Field(5, ge=1, le=50, description="Max results per collection")


class FindRelatedResponse(BaseModel):
    """Response schema for POST /find-related."""
    entity: str
    results: Dict[str, List[EvidenceItem]]
    total_hits: int


class KnowledgeStatsResponse(BaseModel):
    """Response schema for GET /knowledge/stats."""
    pharmacogenes: int
    metabolizer_phenotypes: int
    drug_categories: int
    drugs_tracked: int
    cyp_inhibitor_enzymes: int
    cyp_inducer_enzymes: int
    hla_drug_associations: int
    drug_alternative_mappings: int
    activity_score_tables: int
    entity_aliases: int


# =====================================================================
# Helper -- convert internal SearchHit to API EvidenceItem
# =====================================================================

def _hit_to_evidence(hit) -> EvidenceItem:
    """Convert an internal SearchHit to the API EvidenceItem schema."""
    return EvidenceItem(
        collection=hit.collection,
        id=hit.id,
        score=hit.score,
        text=hit.text,
        metadata=hit.metadata,
    )


def _alert_to_item(alert) -> AlertItem:
    """Convert an internal PGxAlert to the API AlertItem schema."""
    return AlertItem(
        alert_level=alert.alert_level.value if hasattr(alert.alert_level, "value") else str(alert.alert_level),
        gene=alert.gene,
        drug=alert.drug,
        phenotype=alert.phenotype,
        recommendation=alert.recommendation,
        evidence_pmids=alert.evidence_pmids,
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
    Supports gene/drug filtering, collection restriction, and patient
    profile-based queries.
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
        engine_query = AgentQuery(
            question=request.question,
            patient_id=request.patient_id,
            medication_list=request.medication_list,
            gene=request.gene_filter,
            drug=request.drug_filter,
        )

        # Retrieve evidence
        evidence = _engine.retrieve(
            query=engine_query,
            collections_filter=request.collections,
        )

        # Generate LLM response
        from src.rag_engine import PGX_SYSTEM_PROMPT
        prompt_text = _engine._build_prompt(request.question, evidence)
        answer = _engine.llm.generate(
            prompt=prompt_text,
            system_prompt=PGX_SYSTEM_PROMPT,
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

    Supports gene/drug filtering and collection restriction.
    """
    _metrics["requests_total"] += 1
    _metrics["search_requests_total"] += 1

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")
    if not _engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        engine_query = AgentQuery(
            question=request.question,
            gene=request.gene_filter,
            drug=request.drug_filter,
        )

        evidence = _engine.retrieve(
            query=engine_query,
            collections_filter=request.collections,
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
    """Find all evidence related to an entity across all 15 collections.

    Enables queries like "show me everything about CYP2D6" spanning
    gene references, drug guidelines, interactions, HLA screening,
    dosing algorithms, clinical evidence, and more.
    """
    _metrics["requests_total"] += 1
    _metrics["find_related_requests_total"] += 1

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")
    if not _engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        raw_results = _engine.find_related(
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
    """Return statistics about the pharmacogenomics knowledge graph.

    Includes counts of pharmacogenes, metabolizer phenotypes, drug
    categories, CYP inhibitors/inducers, HLA-drug associations,
    drug alternatives, activity score tables, and entity aliases.
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
    """Prometheus-compatible metrics endpoint.

    Returns basic request counters in Prometheus exposition format.
    A full implementation would integrate with prometheus_client.
    """
    lines = [
        "# HELP pgx_api_requests_total Total API requests",
        "# TYPE pgx_api_requests_total counter",
        f'pgx_api_requests_total {_metrics["requests_total"]}',
        "",
        "# HELP pgx_api_query_requests_total Total /query requests",
        "# TYPE pgx_api_query_requests_total counter",
        f'pgx_api_query_requests_total {_metrics["query_requests_total"]}',
        "",
        "# HELP pgx_api_search_requests_total Total /search requests",
        "# TYPE pgx_api_search_requests_total counter",
        f'pgx_api_search_requests_total {_metrics["search_requests_total"]}',
        "",
        "# HELP pgx_api_find_related_requests_total Total /find-related requests",
        "# TYPE pgx_api_find_related_requests_total counter",
        f'pgx_api_find_related_requests_total {_metrics["find_related_requests_total"]}',
        "",
        "# HELP pgx_api_errors_total Total error responses",
        "# TYPE pgx_api_errors_total counter",
        f'pgx_api_errors_total {_metrics["errors_total"]}',
        "",
    ]

    # Add collection vector counts if available
    if _manager:
        try:
            stats = _manager.get_collection_stats()
            lines.append("# HELP pgx_collection_vectors Number of vectors per collection")
            lines.append("# TYPE pgx_collection_vectors gauge")
            for name, count in stats.items():
                lines.append(f'pgx_collection_vectors{{collection="{name}"}} {count}')
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
