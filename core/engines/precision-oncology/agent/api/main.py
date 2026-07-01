"""
Oncology Intelligence Agent - FastAPI Application
===================================================
Phase 4: API + UI layer for the Oncology Intelligence Agent.
Provides RAG-powered clinical decision support for molecular tumor boards.

Author: Adam Jones
Date: February 2026
"""

import time
import logging
from collections import defaultdict
from contextlib import asynccontextmanager
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field
from sentence_transformers import SentenceTransformer

from config.settings import settings as _settings_instance
from src.collections import OncoCollectionManager
from src.rag_engine import OncoRAGEngine
from src.agent import OncoIntelligenceAgent
from src.case_manager import OncologyCaseManager
from src.trial_matcher import TrialMatcher
from src.therapy_ranker import TherapyRanker
from src.cross_modal import OncoCrossModalTrigger

from api.routes import meta_agent, cases, trials, reports, events

import src.knowledge as knowledge_module

logger = logging.getLogger(__name__)


class EmbedderWrapper:
    """Adapter giving SentenceTransformer both .encode() and .embed() APIs."""

    def __init__(self, model: SentenceTransformer):
        self._model = model

    def encode(self, texts):
        return self._model.encode(texts)

    def embed(self, text):
        """Return a single embedding vector (list of floats)."""
        if isinstance(text, str):
            return self._model.encode([text])[0].tolist()
        return self._model.encode(text).tolist()


# ---------------------------------------------------------------------------
# Global application state
# ---------------------------------------------------------------------------
_state: Dict[str, Any] = {}

VERSION = "0.1.0"


def get_state() -> Dict[str, Any]:
    """Return the shared application state dict."""
    return _state


# ---------------------------------------------------------------------------
# Lifespan
# ---------------------------------------------------------------------------
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown lifecycle for the Oncology Intelligence Agent."""
    logger.info("Oncology Intelligence Agent starting up ...")

    # -- Settings ----------------------------------------------------------
    settings = _settings_instance
    _state["settings"] = settings

    # -- Milvus ------------------------------------------------------------
    collection_manager = OncoCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    collection_manager.connect()
    _state["collection_manager"] = collection_manager

    # -- Embedder ----------------------------------------------------------
    embedder = EmbedderWrapper(SentenceTransformer(settings.EMBEDDING_MODEL))
    _state["embedder"] = embedder

    # -- LLM Client --------------------------------------------------------
    llm_client = None
    try:
        import anthropic

        class _OncoLLMClient:
            def __init__(self):
                self._client = anthropic.Anthropic()
                self._model = getattr(settings, 'LLM_MODEL', 'claude-sonnet-4-20250514')

            def chat(self, messages):
                system = ""
                user_msgs = []
                for m in messages:
                    if m["role"] == "system":
                        system = m["content"]
                    else:
                        user_msgs.append(m)
                resp = self._client.messages.create(
                    model=self._model, max_tokens=4096,
                    system=system, messages=user_msgs,
                )
                return resp.content[0].text

            def chat_stream(self, messages):
                system = ""
                user_msgs = []
                for m in messages:
                    if m["role"] == "system":
                        system = m["content"]
                    else:
                        user_msgs.append(m)
                with self._client.messages.stream(
                    model=self._model, max_tokens=4096,
                    system=system, messages=user_msgs,
                ) as stream:
                    for text in stream.text_stream:
                        yield text

        llm_client = _OncoLLMClient()
        logger.info("Anthropic LLM client initialized for Oncology RAG")
    except Exception as exc:
        logger.warning("LLM client unavailable: %s", exc)

    # -- RAG Engine --------------------------------------------------------
    rag_engine = OncoRAGEngine(
        collection_manager=collection_manager,
        embedder=embedder,
        settings=settings,
        llm_client=llm_client,
    )
    _state["rag_engine"] = rag_engine

    # -- Intelligence Agent ------------------------------------------------
    intelligence_agent = OncoIntelligenceAgent(
        rag_engine=rag_engine,
    )
    _state["intelligence_agent"] = intelligence_agent

    # -- Case Manager ------------------------------------------------------
    case_manager = OncologyCaseManager(
        collection_manager=collection_manager,
        embedder=embedder,
        knowledge=knowledge_module,
        rag_engine=rag_engine,
    )
    _state["case_manager"] = case_manager

    # -- Trial Matcher -----------------------------------------------------
    trial_matcher = TrialMatcher(
        collection_manager=collection_manager,
        embedder=embedder,
    )
    _state["trial_matcher"] = trial_matcher

    # -- Therapy Ranker ----------------------------------------------------
    therapy_ranker = TherapyRanker(
        collection_manager=collection_manager,
        embedder=embedder,
        knowledge=knowledge_module,
    )
    _state["therapy_ranker"] = therapy_ranker

    # -- Cross-Modal Trigger -----------------------------------------------
    cross_modal = OncoCrossModalTrigger(
        collection_manager=collection_manager,
        embedder=embedder,
        settings={
            "cross_modal_threshold": settings.CROSS_MODAL_THRESHOLD,
            "genomic_top_k": settings.GENOMIC_TOP_K,
            "imaging_top_k": settings.IMAGING_TOP_K,
        },
    )
    _state["cross_modal"] = cross_modal

    logger.info("All services initialised successfully.")

    yield  # --- application runs here ---

    # -- Shutdown ----------------------------------------------------------
    logger.info("Oncology Intelligence Agent shutting down ...")
    try:
        collection_manager.disconnect()
    except Exception as exc:
        logger.warning("Error disconnecting from Milvus: %s", exc)
    _state.clear()
    logger.info("Shutdown complete.")


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------
app = FastAPI(
    title="Oncology Intelligence MTB API",
    description=(
        "RAG-powered clinical decision support for molecular tumor boards. "
        "Part of the HCLS AI Factory pipeline."
    ),
    version=VERSION,
    docs_url="/docs",
    openapi_url="/openapi.json",
    lifespan=lifespan,
)

# -- CORS ------------------------------------------------------------------
_cors_origins = [
    o.strip() for o in _settings_instance.CORS_ORIGINS.split(",") if o.strip()
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# -- Auth middleware (optional, based on API_KEY setting) -------------------
_AUTH_SKIP_PATHS = {"/health", "/healthz", "/metrics", "/docs", "/openapi.json"}


@app.middleware("http")
async def check_api_key(request: Request, call_next):
    if request.url.path in _AUTH_SKIP_PATHS:
        return await call_next(request)
    api_key = getattr(_settings_instance, 'API_KEY', '') or ''
    if not api_key:
        return await call_next(request)
    provided = request.headers.get("X-API-Key", "")
    if provided != api_key:
        return JSONResponse(status_code=401, content={"detail": "Invalid or missing API key"})
    return await call_next(request)


# -- Rate limiting middleware ----------------------------------------------
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


# -- Request size limit ----------------------------------------------------
_MAX_BODY = _settings_instance.MAX_REQUEST_SIZE_MB * 1024 * 1024


@app.middleware("http")
async def limit_request_size(request: Request, call_next):
    """Reject requests exceeding the configured body size limit."""
    content_length = request.headers.get("content-length")
    if content_length and int(content_length) > _MAX_BODY:
        return JSONResponse(
            status_code=413,
            content={"detail": f"Request body exceeds {_settings_instance.MAX_REQUEST_SIZE_MB} MB limit"},
        )
    return await call_next(request)

# -- Include routers -------------------------------------------------------
app.include_router(meta_agent.router)
app.include_router(cases.router)
app.include_router(trials.router)
app.include_router(reports.router)
app.include_router(events.router)


# ---------------------------------------------------------------------------
# Root endpoint
# ---------------------------------------------------------------------------

@app.get("/", include_in_schema=False)
def root():
    return {"service": "Oncology Intelligence MTB", "docs": "/docs", "health": "/health"}


# ---------------------------------------------------------------------------
# Request / Response schemas
# ---------------------------------------------------------------------------
class QueryRequest(BaseModel):
    question: str
    cancer_type: Optional[str] = None
    gene: Optional[str] = None
    top_k: int = Field(default=10, ge=1, le=100)


class SearchRequest(BaseModel):
    question: str
    cancer_type: Optional[str] = None
    gene: Optional[str] = None
    top_k: int = Field(default=10, ge=1, le=100)


class FindRelatedRequest(BaseModel):
    entity: str
    entity_type: str = "gene"
    top_k: int = Field(default=10, ge=1, le=50)


# ---------------------------------------------------------------------------
# Core endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
async def health():
    """Service health check with collection statistics."""
    cm: OncoCollectionManager = _state.get("collection_manager")
    if cm is None:
        raise HTTPException(status_code=503, detail="Service not initialised")

    collection_info = {}
    total_vectors = 0
    try:
        for name in cm.list_collections():
            count = cm.get_collection_count(name)
            collection_info[name] = count
            total_vectors += count
    except Exception as exc:
        logger.warning("Failed to gather collection stats: %s", exc)

    services = {
        "milvus": cm.is_connected(),
        "embedder": _state.get("embedder") is not None,
        "rag_engine": _state.get("rag_engine") is not None,
        "intelligence_agent": _state.get("intelligence_agent") is not None,
        "case_manager": _state.get("case_manager") is not None,
        "trial_matcher": _state.get("trial_matcher") is not None,
        "therapy_ranker": _state.get("therapy_ranker") is not None,
    }

    return {
        "status": "healthy" if all(services.values()) else "degraded",
        "collections": collection_info,
        "total_vectors": total_vectors,
        "version": VERSION,
        "services": services,
    }


@app.get("/collections")
async def list_collections():
    """List all oncology knowledge collections with entity counts."""
    cm: OncoCollectionManager = _state.get("collection_manager")
    if cm is None:
        raise HTTPException(status_code=503, detail="Service not initialised")

    result = []
    for name in cm.list_collections():
        result.append({"name": name, "count": cm.get_collection_count(name)})
    return {"collections": result}


@app.post("/query")
async def query(req: QueryRequest):
    """Full RAG query — retrieves evidence and generates an LLM answer."""
    rag: OncoRAGEngine = _state.get("rag_engine")
    if rag is None:
        raise HTTPException(status_code=503, detail="RAG engine not initialised")

    t0 = time.time()
    try:
        answer = rag.query(question=req.question, top_k=req.top_k)
        elapsed_ms = round((time.time() - t0) * 1000, 1)
        return {
            "question": req.question,
            "answer": answer,
            "processing_time_ms": elapsed_ms,
        }
    except Exception as exc:
        logger.error("Query failed: %s", exc)
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/search")
async def search(req: SearchRequest):
    """Evidence-only vector search (no LLM generation)."""
    rag: OncoRAGEngine = _state.get("rag_engine")
    if rag is None:
        raise HTTPException(status_code=503, detail="RAG engine not initialised")

    t0 = time.time()
    try:
        hits = rag.search(question=req.question, top_k=req.top_k)
        elapsed_ms = round((time.time() - t0) * 1000, 1)
        return {"results": [h.model_dump() if hasattr(h, "model_dump") else h for h in hits],
                "count": len(hits), "processing_time_ms": elapsed_ms}
    except Exception as exc:
        logger.error("Search failed: %s", exc)
        raise HTTPException(status_code=500, detail="Internal processing error")


@app.post("/find-related")
async def find_related(req: FindRelatedRequest):
    """Cross-collection entity linking — find related knowledge across domains."""
    cross: OncoCrossModalTrigger = _state.get("cross_modal")
    if cross is None:
        raise HTTPException(status_code=503, detail="Cross-modal service not ready")

    t0 = time.time()
    try:
        related = cross.find_related(
            entity=req.entity,
            entity_type=req.entity_type,
            top_k=req.top_k,
        )
    except AttributeError:
        # find_related may not be implemented — fall back to evaluate
        related = []
    except Exception as exc:
        logger.error("Find-related failed: %s", exc)
        raise HTTPException(status_code=500, detail="Internal processing error")
    elapsed_ms = round((time.time() - t0) * 1000, 1)
    return {"entity": req.entity, "related": related, "processing_time_ms": elapsed_ms}


@app.get("/knowledge/stats")
async def knowledge_stats():
    """Aggregate knowledge-base statistics."""
    cm: OncoCollectionManager = _state.get("collection_manager")
    if cm is None:
        raise HTTPException(status_code=503, detail="Service not initialised")

    counts = {}
    for name in cm.list_collections():
        counts[name] = cm.get_collection_count(name)

    return {
        "target_count": counts.get("onco_variants", 0),
        "therapy_count": counts.get("onco_therapies", 0),
        "resistance_count": counts.get("onco_resistance", 0),
        "pathway_count": counts.get("onco_pathways", 0),
        "biomarker_count": counts.get("onco_biomarkers", 0),
        "collection_counts": counts,
    }


@app.get("/metrics")
async def metrics():
    """Prometheus-compatible metrics endpoint."""
    cm: OncoCollectionManager = _state.get("collection_manager")
    lines: List[str] = []

    lines.append("# HELP onco_agent_up Service availability.")
    lines.append("# TYPE onco_agent_up gauge")
    lines.append(f"onco_agent_up {1 if cm and cm.is_connected() else 0}")

    if cm:
        lines.append("# HELP onco_collection_vectors Vector count per collection.")
        lines.append("# TYPE onco_collection_vectors gauge")
        for name in cm.list_collections():
            count = cm.get_collection_count(name)
            lines.append(f'onco_collection_vectors{{collection="{name}"}} {count}')

    lines.append("")
    from fastapi.responses import PlainTextResponse

    return PlainTextResponse("\n".join(lines), media_type="text/plain; charset=utf-8")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import uvicorn

    logging.basicConfig(level=logging.INFO)
    uvicorn.run(
        "api.main:app",
        host=_settings_instance.API_HOST,
        port=_settings_instance.API_PORT,
        reload=False,
        log_level="info",
    )
