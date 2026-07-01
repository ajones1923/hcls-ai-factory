"""Developability + design FastAPI service (B3/B4). Port 8576. No future-annotations."""
from typing import Optional
from developability import DevelopabilityScorer, SequenceError


def create_app(scorer: Optional[DevelopabilityScorer] = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel
    sc = scorer or DevelopabilityScorer()
    app = FastAPI(title="Protein Developability Service", version="1.0")

    class SeqReq(BaseModel):
        sequence: str
        n: int = 5

    @app.get("/health")
    def health(): return {"status": "ok"}

    @app.post("/score")
    def score(req: SeqReq):
        try: return sc.score(req.sequence)
        except SequenceError as e: raise HTTPException(422, str(e))

    @app.post("/optimize")
    def optimize(req: SeqReq):
        try: return sc.optimize(req.sequence, n=req.n)
        except SequenceError as e: raise HTTPException(422, str(e))

    return app


def _app_factory(): return create_app()
