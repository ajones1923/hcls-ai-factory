"""Single-cell compute FastAPI service + client (D1). Default port 8573."""
from typing import Any
from single_cell_compute import SingleCellAnalysis


def create_app(analysis: SingleCellAnalysis | None = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel
    sca = analysis or SingleCellAnalysis()
    app = FastAPI(title="Single-Cell Compute Service", version="1.0")

    class AnalyzeReq(BaseModel):
        h5ad_path: str | None = None
        demo: bool = False           # use the bundled PBMC 3k demo dataset
        resolution: float = 1.0

    @app.get("/health")
    def health():
        return {"status": "ok"}

    @app.post("/analyze")
    def analyze(req: AnalyzeReq):
        import scanpy as sc
        if req.demo:
            adata = sc.datasets.pbmc3k()
        elif req.h5ad_path:
            adata = sc.read_h5ad(req.h5ad_path)
        else:
            raise HTTPException(422, "provide h5ad_path or demo=true")
        return sca.run(adata, resolution=req.resolution)

    return app


def _app_factory():
    return create_app()


class SingleCellClient:
    def __init__(self, endpoint: str = "http://localhost:8573", timeout: float = 600.0) -> None:
        self.endpoint = endpoint.rstrip("/"); self.timeout = timeout

    def analyze(self, h5ad_path: str | None = None, demo: bool = False, resolution: float = 1.0) -> dict[str, Any]:
        import httpx
        r = httpx.post(f"{self.endpoint}/analyze",
                       json={"h5ad_path": h5ad_path, "demo": demo, "resolution": resolution},
                       timeout=self.timeout)
        r.raise_for_status(); return r.json()
