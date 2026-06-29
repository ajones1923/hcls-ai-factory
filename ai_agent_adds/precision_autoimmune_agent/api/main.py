"""
Precision Autoimmune Agent — FastAPI Server

Provides REST endpoints for:
- RAG-powered clinical queries
- Patient analysis pipeline
- Document ingestion
- Collection management
- Health monitoring

Ports: 8532 (API)
Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import json
import os
import sys
import time
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from loguru import logger
from pydantic import BaseModel

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config.logging import configure_logging
from config.settings import settings
from src.agent import AutoimmuneAgent
from src.collections import AutoimmuneCollectionManager
from src.diagnostic_engine import DiagnosticEngine
from src.document_processor import DocumentProcessor
from src.models import AutoimmunePatientProfile
from src.rag_engine import AutoimmuneRAGEngine
from src.timeline_builder import TimelineBuilder

# ── Global state ──────────────────────────────────────────────────────────
_state: Dict[str, Any] = {}


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown lifecycle with validation and graceful teardown."""

    # ── Centralized logging ──────────────────────────────────────────────
    configure_logging(
        log_level=settings.LOG_LEVEL,
        log_dir=settings.LOG_DIR or None,
        service_name="autoimmune-api",
    )

    logger.info("=" * 60)
    logger.info("  Precision Autoimmune Intelligence Agent — API")
    logger.info("=" * 60)

    # ── Startup validation ───────────────────────────────────────────────
    service_status: Dict[str, str] = {}

    # Check ANTHROPIC_API_KEY
    api_key = settings.ANTHROPIC_API_KEY or os.environ.get("ANTHROPIC_API_KEY", "")
    if not api_key:
        logger.warning("ANTHROPIC_API_KEY is not set — LLM features will be unavailable (demo mode)")
        service_status["llm"] = "DEGRADED (no API key)"
    else:
        service_status["llm"] = "PENDING"

    # Collection manager with retry
    cm = AutoimmuneCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
        embedding_dim=settings.EMBEDDING_DIM,
    )
    milvus_connected = False
    for attempt in range(1, 3):  # 2 attempts
        try:
            cm.connect()
            milvus_connected = True
            service_status["milvus"] = "OK"
            break
        except Exception as exc:
            logger.warning(f"Milvus connection attempt {attempt}/2 failed: {exc}")
            if attempt < 2:
                time.sleep(2)
    if not milvus_connected:
        logger.warning("Milvus unavailable — vector search will be degraded")
        service_status["milvus"] = "DEGRADED"

    # Embedder
    embedder = None
    try:
        from sentence_transformers import SentenceTransformer
        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
        logger.info(f"Loaded embedding model: {settings.EMBEDDING_MODEL}")
        service_status["embedder"] = "OK"
    except Exception as exc:
        logger.warning(f"Embedder load deferred: {exc}")
        service_status["embedder"] = "DEGRADED"

    # LLM client
    llm = None
    if api_key:
        try:
            import anthropic
            llm = anthropic.Anthropic(api_key=api_key)
            service_status["llm"] = "OK"
        except Exception as exc:
            logger.warning(f"LLM client deferred: {exc}")
            service_status["llm"] = "DEGRADED"

    # RAG engine
    from src.knowledge import KNOWLEDGE_VERSION
    rag = AutoimmuneRAGEngine(
        collection_manager=cm,
        embedder=embedder,
        llm_client=llm,
        settings=settings,
        knowledge=KNOWLEDGE_VERSION,
    )

    # Agent
    agent = AutoimmuneAgent()

    # Document processor
    doc_processor = DocumentProcessor(
        collection_manager=cm,
        embedder=embedder,
        max_chunk_size=settings.MAX_CHUNK_SIZE,
        chunk_overlap=settings.CHUNK_OVERLAP,
    )

    # Diagnostic engine
    diag = DiagnosticEngine(agent=agent, rag_engine=rag)

    # Timeline builder
    timeline = TimelineBuilder(rag_engine=rag)

    _state.update({
        "cm": cm,
        "embedder": embedder,
        "llm": llm,
        "rag": rag,
        "agent": agent,
        "doc_processor": doc_processor,
        "diagnostic": diag,
        "timeline": timeline,
        "start_time": time.time(),
    })

    # ── Startup banner ───────────────────────────────────────────────────
    logger.info("-" * 60)
    logger.info("  Service Status:")
    for svc, status in service_status.items():
        marker = "+" if status == "OK" else "!"
        logger.info(f"    [{marker}] {svc:.<20s} {status}")
    logger.info(f"  API port: {settings.API_PORT}")
    logger.info(f"  UI port:  {settings.STREAMLIT_PORT}")
    logger.info("-" * 60)
    logger.info("Autoimmune Agent API ready")

    yield

    # ── Graceful shutdown ────────────────────────────────────────────────
    logger.info("Shutting down Autoimmune Agent API...")
    try:
        cm.disconnect()
        logger.info("Milvus connection closed")
    except Exception as exc:
        logger.warning(f"Error disconnecting Milvus: {exc}")
    logger.info("Autoimmune Agent API shutdown complete")


# ── App ───────────────────────────────────────────────────────────────────

app = FastAPI(
    title="Precision Autoimmune Intelligence Agent",
    version="1.0.0",
    description="Multi-collection RAG agent for autoimmune disease analysis",
    lifespan=lifespan,
)

# CORS
origins = settings.CORS_ORIGINS.split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── API key authentication ──────────────────────────────────────────────
@app.middleware("http")
async def authenticate(request: Request, call_next):
    """Optional API key authentication (skip if AUTO_API_KEY not set)."""
    api_key = settings.API_KEY
    if not api_key:
        return await call_next(request)
    # Skip auth for health endpoints and root
    if request.url.path in ("/", "/health", "/healthz", "/metrics"):
        return await call_next(request)
    provided = request.headers.get("X-API-Key") or request.query_params.get("api_key")
    if provided != api_key:
        return JSONResponse(status_code=401, content={"detail": "Invalid or missing API key"})
    return await call_next(request)


@app.middleware("http")
async def add_timing(request: Request, call_next):
    """Add request timing header."""
    start = time.time()
    response = await call_next(request)
    response.headers["X-Process-Time-Ms"] = f"{(time.time() - start) * 1000:.1f}"
    return response


# ── Request size limiter ──────────────────────────────────────────────────
@app.middleware("http")
async def limit_request_size(request: Request, call_next):
    max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
    content_length = request.headers.get("content-length")
    if content_length and int(content_length) > max_bytes:
        return JSONResponse(status_code=413, content={"detail": "Request too large"})
    return await call_next(request)


# ── Request / Response models ─────────────────────────────────────────────

class QueryRequest(BaseModel):
    question: str
    patient_id: Optional[str] = None
    patient_context: Optional[str] = None
    collections_filter: Optional[List[str]] = None
    top_k: Optional[int] = None


class QueryResponse(BaseModel):
    answer: str
    evidence_count: int
    collections_searched: int
    search_time_ms: float


class SearchRequest(BaseModel):
    question: str
    patient_id: Optional[str] = None
    collections_filter: Optional[List[str]] = None
    top_k: Optional[int] = None


class DifferentialRequest(BaseModel):
    positive_antibodies: List[str]
    hla_alleles: Optional[List[str]] = None
    symptoms: Optional[List[str]] = None


# ── Core endpoints ────────────────────────────────────────────────────────

@app.get("/")
async def root():
    return {
        "service": "Precision Autoimmune Intelligence Agent",
        "version": "1.0.0",
        "status": "running",
        "ports": {"api": settings.API_PORT, "ui": settings.STREAMLIT_PORT},
    }


@app.get("/health")
async def health():
    cm = _state.get("cm")
    milvus_ok = False
    collections = []
    total_vectors = 0
    if cm:
        try:
            stats = cm.get_collection_stats()
            milvus_ok = True
            collections = list(stats.keys())
            total_vectors = sum(v for v in stats.values() if v > 0)
        except Exception:
            pass

    return {
        "status": "healthy",
        "service": "autoimmune-agent",
        "milvus_connected": milvus_ok,
        "collections": len(collections),
        "total_vectors": total_vectors,
        "embedder_loaded": _state.get("embedder") is not None,
        "llm_available": _state.get("llm") is not None,
        "uptime_seconds": round(time.time() - _state.get("start_time", time.time())),
    }


@app.get("/healthz")
async def healthz():
    """Lightweight health probe for landing page."""
    return {"status": "ok"}


# ── Cross-Agent Integration Endpoint ──────────────────────────────────────

@app.post("/v1/autoimmune/integrated-assessment")
async def integrated_assessment(request: dict):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, cardiology, neurology, pharmacogenomics, biomarker,
    and clinical trial agents for a comprehensive autoimmune assessment.
    """
    try:
        from src.cross_modal import (
            integrate_cross_agent_results,
            query_biomarker_agent,
            query_cardiology_agent,
            query_neurology_agent,
            query_oncology_agent,
            query_pgx_agent,
            query_trial_agent,
        )

        patient_profile = request.get("patient_profile", {})
        patient_id = request.get("patient_id", "")
        drug_list = request.get("drug_list", [])
        autoantibody_panel = request.get("autoantibody_panel", {})
        autoimmune_profile = request.get("autoimmune_profile", {})

        results = []

        # Query oncology agent for irAE assessment
        if patient_profile:
            results.append(query_oncology_agent(patient_profile))

        # Query cardiology agent for myocarditis risk
        if patient_id:
            results.append(query_cardiology_agent(patient_id))

        # Query neurology agent for autoimmune encephalitis risk
        if patient_id:
            results.append(query_neurology_agent(patient_id))

        # Query PGx agent for biologic drug-gene interactions
        if drug_list:
            results.append(query_pgx_agent(drug_list))

        # Query biomarker agent for autoantibody correlation
        if autoantibody_panel:
            results.append(query_biomarker_agent(autoantibody_panel))

        # Query trial agent for autoimmune trial matching
        if autoimmune_profile:
            results.append(query_trial_agent(autoimmune_profile))

        integrated = integrate_cross_agent_results(results)
        return {
            "status": "completed",
            "assessment": integrated,
            "agents_consulted": integrated.get("agents_consulted", []),
        }
    except Exception as exc:
        logger.error(f"Integrated assessment failed: {exc}")
        return {"status": "partial", "assessment": {}, "error": "Cross-agent integration unavailable"}


@app.post("/query", response_model=QueryResponse)
async def query(req: QueryRequest):
    """Full RAG query: retrieve evidence + synthesize with Claude."""
    rag = _state.get("rag")
    if not rag:
        raise HTTPException(503, "RAG engine not initialized")

    answer = rag.query(
        question=req.question,
        patient_context=req.patient_context,
        patient_id=req.patient_id,
        collections_filter=req.collections_filter,
        top_k_per_collection=req.top_k,
    )

    result = rag.retrieve(req.question, patient_id=req.patient_id)

    return QueryResponse(
        answer=answer,
        evidence_count=len(result.hits),
        collections_searched=result.total_collections_searched,
        search_time_ms=result.search_time_ms,
    )


@app.post("/query/stream")
async def query_stream(req: QueryRequest):
    """Streaming RAG query: retrieve evidence + stream Claude response via SSE."""
    rag = _state.get("rag")
    if not rag:
        raise HTTPException(503, "RAG engine not initialized")

    async def event_generator():
        try:
            for chunk in rag.query_stream(
                question=req.question,
                patient_context=req.patient_context,
                patient_id=req.patient_id,
                collections_filter=req.collections_filter,
            ):
                yield f"data: {json.dumps({'text': chunk})}\n\n"
            yield "data: [DONE]\n\n"
        except Exception as exc:
            logger.error(f"Streaming query failed: {exc}")
            yield f"data: {json.dumps({'error': 'Internal processing error'})}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")


@app.post("/search")
async def search(req: SearchRequest):
    """Evidence-only search without LLM synthesis."""
    rag = _state.get("rag")
    if not rag:
        raise HTTPException(503, "RAG engine not initialized")

    hits = rag.search(
        question=req.question,
        patient_id=req.patient_id,
        collections_filter=req.collections_filter,
        top_k_per_collection=req.top_k,
    )

    return {
        "hits": [
            {
                "collection": h.collection,
                "id": h.id,
                "score": h.score,
                "text": h.text[:1000],
                "relevance": h.relevance,
                "metadata": h.metadata,
            }
            for h in hits
        ],
        "total": len(hits),
    }


@app.post("/analyze")
async def analyze_patient(profile: AutoimmunePatientProfile):
    """Run full autoimmune analysis pipeline."""
    agent = _state.get("agent")
    if not agent:
        raise HTTPException(503, "Agent not initialized")

    result = agent.analyze_patient(profile)
    return result.model_dump()


@app.post("/differential")
async def differential_diagnosis(req: DifferentialRequest):
    """Generate differential diagnosis from antibodies and HLA data."""
    diag = _state.get("diagnostic")
    if not diag:
        raise HTTPException(503, "Diagnostic engine not initialized")

    differential = diag.generate_differential(
        positive_antibodies=req.positive_antibodies,
        hla_alleles=req.hla_alleles,
        symptoms=req.symptoms,
    )
    return {"differential": differential}


# ── Document ingestion endpoints ──────────────────────────────────────────

@app.post("/ingest/upload")
async def upload_document(
    file: UploadFile = File(...),
    patient_id: str = Form(...),
):
    """Upload and ingest a clinical PDF document."""
    doc_proc = _state.get("doc_processor")
    if not doc_proc:
        raise HTTPException(503, "Document processor not initialized")

    if not file.filename.lower().endswith(".pdf"):
        raise HTTPException(400, "Only PDF files are supported")

    # Save to temp location
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        content = await file.read()
        tmp.write(content)
        tmp_path = Path(tmp.name)

    try:
        records = doc_proc.process_pdf(tmp_path, patient_id)
        if not records:
            raise HTTPException(422, "No text extracted from PDF")

        records = doc_proc.embed_records(records)
        count = doc_proc.collection_manager.insert_batch(
            settings.COLL_CLINICAL_DOCUMENTS, records
        )

        return {
            "status": "success",
            "filename": file.filename,
            "patient_id": patient_id,
            "chunks_ingested": count,
            "doc_type": records[0].get("doc_type") if records else "unknown",
            "specialty": records[0].get("specialty") if records else "unknown",
        }
    finally:
        tmp_path.unlink(missing_ok=True)


@app.post("/ingest/demo-data")
async def ingest_demo_data():
    """Ingest all demo patient data from demo_data/ directory."""
    doc_proc = _state.get("doc_processor")
    if not doc_proc:
        raise HTTPException(503, "Document processor not initialized")

    demo_dir = settings.DEMO_DATA_DIR
    if not demo_dir.exists():
        raise HTTPException(404, f"Demo data directory not found: {demo_dir}")

    results = doc_proc.ingest_demo_data(demo_dir, settings.COLL_CLINICAL_DOCUMENTS)

    total = sum(results.values())
    return {
        "status": "success",
        "patients": results,
        "total_chunks": total,
    }


# ── Collection management endpoints ──────────────────────────────────────

@app.get("/collections")
async def list_collections():
    cm = _state.get("cm")
    if not cm:
        raise HTTPException(503, "Collection manager not initialized")

    stats = cm.get_collection_stats()
    return {
        "collections": [
            {"name": name, "count": count}
            for name, count in sorted(stats.items())
        ],
        "total_collections": len(stats),
        "total_vectors": sum(v for v in stats.values() if v > 0),
    }


@app.post("/collections/create")
async def create_collections(drop_existing: bool = False):
    cm = _state.get("cm")
    if not cm:
        raise HTTPException(503, "Collection manager not initialized")

    collections = cm.create_all_collections(drop_existing=drop_existing)
    return {
        "status": "success",
        "collections_created": list(collections.keys()),
        "count": len(collections),
    }


# ── Export endpoints ──────────────────────────────────────────────────────

class ExportRequest(BaseModel):
    patient_id: str
    format: str = "markdown"  # markdown, fhir, pdf
    query_answer: Optional[str] = None


@app.post("/export")
async def export_report(req: ExportRequest):
    """Export analysis report in specified format."""
    from src.export import AutoimmuneExporter
    from src.knowledge import KNOWLEDGE_VERSION

    exporter = AutoimmuneExporter(knowledge_version=KNOWLEDGE_VERSION)

    if req.format == "markdown":
        md = exporter.to_markdown(req.patient_id, query_answer=req.query_answer)
        return {"format": "markdown", "content": md}
    elif req.format == "fhir":
        fhir = exporter.to_fhir_r4(req.patient_id)
        return {"format": "fhir", "content": fhir}
    elif req.format == "pdf":
        pdf_bytes = exporter.to_pdf_bytes(req.patient_id, query_answer=req.query_answer)
        import base64
        return {"format": "pdf", "content_base64": base64.b64encode(pdf_bytes).decode()}
    else:
        raise HTTPException(400, f"Unsupported format: {req.format}")


# ── Metrics endpoint ─────────────────────────────────────────────────────

@app.get("/metrics")
async def metrics():
    """Prometheus-compatible metrics."""
    cm = _state.get("cm")
    lines = [
        "# HELP autoimmune_agent_up Whether the agent is running",
        "# TYPE autoimmune_agent_up gauge",
        "autoimmune_agent_up 1",
    ]
    if cm:
        try:
            stats = cm.get_collection_stats()
            lines.append("# HELP autoimmune_collection_vectors Number of vectors per collection")
            lines.append("# TYPE autoimmune_collection_vectors gauge")
            for name, count in stats.items():
                lines.append(f'autoimmune_collection_vectors{{collection="{name}"}} {max(count, 0)}')
        except Exception:
            pass

    uptime = time.time() - _state.get("start_time", time.time())
    lines.extend([
        "# HELP autoimmune_agent_uptime_seconds Agent uptime",
        "# TYPE autoimmune_agent_uptime_seconds gauge",
        f"autoimmune_agent_uptime_seconds {uptime:.0f}",
    ])
    from fastapi.responses import PlainTextResponse
    return PlainTextResponse("\n".join(lines) + "\n", media_type="text/plain; version=0.0.4")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "api.main:app",
        host="0.0.0.0",
        port=settings.API_PORT,
        workers=2,
        log_level="info",
    )
