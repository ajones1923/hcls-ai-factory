"""Variant store FastAPI service + client (E2). Port 8575. No future-annotations."""
from typing import Optional, List
from variant_store import VariantStore


def create_app(store: Optional[VariantStore] = None):
    from fastapi import FastAPI
    from pydantic import BaseModel
    vs = store or VariantStore()
    app = FastAPI(title="Variant Store Service", version="1.0")

    class LoadReq(BaseModel):
        path: str
        sample: str = "HG002"
        limit: Optional[int] = None

    class RegionReq(BaseModel):
        chrom: str
        start: int
        end: int

    @app.get("/health")
    def health(): return {"status": "ok", "n_variants": vs.count()}

    @app.post("/load")
    def load(req: LoadReq): return {"loaded": vs.load_vcf(req.path, req.sample, req.limit), "total": vs.count()}

    @app.post("/query")
    def query(req: RegionReq): return {"variants": vs.query_region(req.chrom, req.start, req.end)}

    @app.get("/stats")
    def stats(): return vs.stats()

    return app


def _app_factory(): return create_app()


class VariantStoreClient:
    def __init__(self, endpoint: str = "http://localhost:8575", timeout: float = 600.0):
        self.endpoint = endpoint.rstrip("/"); self.timeout = timeout
    def stats(self) -> dict:
        import httpx; return httpx.get(f"{self.endpoint}/stats", timeout=self.timeout).json()
    def query(self, chrom, start, end) -> dict:
        import httpx
        return httpx.post(f"{self.endpoint}/query", json={"chrom":chrom,"start":start,"end":end}, timeout=self.timeout).json()
