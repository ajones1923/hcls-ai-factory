"""MHCflurry immunogenicity service (B10). Port 8577. No future-annotations."""
from typing import List, Optional
import immunogenicity


def create_app():
    from fastapi import FastAPI
    from pydantic import BaseModel
    app = FastAPI(title="Immunogenicity (MHCflurry)", version="1.0")

    class ScanReq(BaseModel):
        sequence: str
        alleles: List[str]
        length: int = 9
        threshold: float = 0.5

    class PredReq(BaseModel):
        peptides: List[str]
        alleles: List[str]

    @app.get("/health")
    def health(): return {"status": "ok"}

    @app.post("/scan")
    def scan(req: ScanReq):
        return immunogenicity.scan_sequence(req.sequence, req.alleles, req.length, req.threshold)

    @app.post("/predict")
    def predict(req: PredReq):
        return {"predictions": immunogenicity.predict(req.peptides, req.alleles)}

    return app


def _app_factory(): return create_app()
