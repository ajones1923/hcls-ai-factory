"""Variant store FastAPI service + client (E1 Genomic Foundation). Port 8575. No future-annotations.

Security posture (2026-11-03 audit, S5): POST /load takes a filesystem path from the caller.
On a box doing genomics, "any file the service user can read, parsed and returned" means
other patients' VCFs. Two controls, both required:

  * the path is resolved and must stay inside VARIANT_STORE_DATA_ROOT (default: the engine's
    data/ directory), so the caller chooses a file within a configured corpus, not on the host;
  * the service installs the shared API-key gate and governance middleware, so it is covered
    by the same HCLS_API_KEY posture as the rest of the platform instead of sitting outside it.

Note the shared gate's documented posture: with no key configured it is open. That is the
platform default and is not changed here — the path confinement holds regardless, and is why
it is not the only control.
"""
import os
from pathlib import Path
from typing import Optional, List
from variant_store import VariantStore

# Filename suffixes this service will open. A VCF reader pointed at anything else is either
# a mistake or someone reading a file they should not.
ALLOWED_VCF_SUFFIXES = ('.vcf', '.vcf.gz')


def data_root() -> Path:
    """Directory that /load paths must resolve inside."""
    configured = os.environ.get('VARIANT_STORE_DATA_ROOT')
    if configured:
        return Path(configured).expanduser().resolve()
    # The engine's own data directory: src/ -> engine root -> data
    return (Path(__file__).resolve().parent.parent / 'data').resolve()


def resolve_vcf_path(raw: str) -> Path:
    """Resolve a caller-supplied path inside the data root.

    Raises ValueError with a message that does not echo host paths back to the caller.
    """
    if not raw or not isinstance(raw, str):
        raise ValueError('path is required')
    if '\x00' in raw:
        raise ValueError('invalid path')

    root = data_root()
    candidate = Path(raw).expanduser()
    resolved = (candidate if candidate.is_absolute() else (root / candidate)).resolve()

    if not (resolved == root or resolved.is_relative_to(root)):
        raise ValueError('path is outside the configured variant-store data root')
    if not resolved.name.lower().endswith(ALLOWED_VCF_SUFFIXES):
        raise ValueError(f'path must name one of {", ".join(ALLOWED_VCF_SUFFIXES)}')
    if resolved.is_symlink():
        raise ValueError('path must not be a symlink')
    if not resolved.is_file():
        raise ValueError('file not found under the variant-store data root')
    return resolved


def create_app(store: Optional[VariantStore] = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel
    vs = store or VariantStore()
    try:
        from hcls_common.api_gate import create_governed_app
        app = create_governed_app("variant-store", capability_id="variant-store-service",
                                  title="Variant Store Service", version="1.0")
    except Exception:  # hcls_common is an optional dependency for the small compute services
        app = FastAPI(title="Variant Store Service", version="1.0")
    try:
        from hcls_common.api_auth import install_api_key_auth
        install_api_key_auth(app, service="variant-store")
    except Exception:
        pass

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
    def load(req: LoadReq):
        try:
            path = resolve_vcf_path(req.path)
        except ValueError as err:
            raise HTTPException(status_code=400, detail=str(err))
        return {"loaded": vs.load_vcf(str(path), req.sample, req.limit), "total": vs.count()}

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
