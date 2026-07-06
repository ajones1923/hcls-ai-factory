"""Variant store FastAPI service + client (E1 Genomic Foundation). Port 8575. No future-annotations."""
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

    class MosaicReq(BaseModel):
        vaf_lo: float = 0.02
        vaf_hi: float = 0.35
        min_dp: int = 30

    @app.get("/health")
    def health(): return {"status": "ok", "n_variants": vs.count()}

    @app.post("/load")
    def load(req: LoadReq): return {"loaded": vs.load_vcf(req.path, req.sample, req.limit), "total": vs.count()}

    @app.post("/query")
    def query(req: RegionReq): return {"variants": vs.query_region(req.chrom, req.start, req.end)}

    @app.get("/stats")
    def stats(): return vs.stats()

    @app.get("/qc")
    def qc(): return vs.qc_report()

    @app.post("/mosaic")
    def mosaic(req: MosaicReq):
        # E1 F1: low-VAF mosaic CANDIDATES (evidence, not a classification)
        cands = vs.mosaic_candidates(req.vaf_lo, req.vaf_hi, req.min_dp)
        return {"candidates": cands,
                "params": {"vaf_lo": req.vaf_lo, "vaf_hi": req.vaf_hi, "min_dp": req.min_dp},
                "n": len(cands)}

    return app


def _app_factory(): return create_app()


class VariantStoreClient:
    def __init__(self, endpoint: str = "http://localhost:8575", timeout: float = 600.0):
        self.endpoint = endpoint.rstrip("/"); self.timeout = timeout
    def stats(self) -> dict:
        import httpx; return httpx.get(f"{self.endpoint}/stats", timeout=self.timeout).json()
    def qc(self) -> dict:
        import httpx; return httpx.get(f"{self.endpoint}/qc", timeout=self.timeout).json()
    def query(self, chrom, start, end) -> dict:
        import httpx
        return httpx.post(f"{self.endpoint}/query", json={"chrom":chrom,"start":start,"end":end}, timeout=self.timeout).json()
    def mosaic(self, vaf_lo=0.02, vaf_hi=0.35, min_dp=30) -> dict:
        import httpx
        return httpx.post(f"{self.endpoint}/mosaic", json={"vaf_lo":vaf_lo,"vaf_hi":vaf_hi,"min_dp":min_dp}, timeout=self.timeout).json()
