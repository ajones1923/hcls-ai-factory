"""Protein embedding + sequence-search FastAPI service (B2). Default port 8571."""

from esmfold_service import SequenceError
from protein_search import ProteinSearchIndex


def create_app(index: ProteinSearchIndex | None = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel

    idx = index or ProteinSearchIndex()
    app = FastAPI(title="Protein Search Service", version="1.0")

    class EmbedReq(BaseModel):
        sequence: str

    class SearchReq(BaseModel):
        sequence: str
        top_k: int = 5

    class IndexReq(BaseModel):
        proteins: list[dict]   # [{"id": int, "name": str, "sequence": str}]

    @app.get("/health")
    def health():
        return {"status": "ok", "model": idx.embedder.model_name, "indexed": idx.count()}

    @app.post("/embed")
    def embed(req: EmbedReq):
        try:
            v = idx.embedder.embed(req.sequence)
        except SequenceError as e:
            raise HTTPException(422, str(e))
        return {"dim": len(v), "embedding": v}

    @app.post("/index")
    def index(req: IndexReq):
        try:
            n = idx.add(req.proteins)
        except SequenceError as e:
            raise HTTPException(422, str(e))
        return {"added": n, "total": idx.count()}

    @app.post("/search")
    def search(req: SearchReq):
        try:
            hits = idx.search(req.sequence, req.top_k)
        except SequenceError as e:
            raise HTTPException(422, str(e))
        return {"query_len": len(req.sequence), "hits": hits}

    return app


def _app_factory():
    return create_app()
