"""Single-cell compute FastAPI service + client (D1). Default port 8573."""
import os
from typing import Any
from single_cell_compute import SingleCellAnalysis

# scanpy writes downloaded datasets to `settings.datasetdir`, which defaults to ./data relative
# to the working directory. In the container that is /app/data, created by the Dockerfile and
# owned by the runtime user; set SC_DATASET_DIR to point somewhere else, such as a mounted
# directory already holding the file on an air-gapped host.
DATASET_DIR = os.environ.get("SC_DATASET_DIR", "")


def load_demo_dataset():
    """The PBMC 3k demo dataset.

    `sc.datasets.pbmc3k()` DOWNLOADS this on first use -- it is not bundled in the image, and
    the Dockerfile comment calling it "the bundled PBMC 3k demo dataset" was wrong. A request
    on a host with no egress, or with an unwritable dataset directory, used to surface as a
    bare 500 from a scanpy traceback. The failure now says what it actually was.
    """
    import scanpy as sc

    if DATASET_DIR:
        sc.settings.datasetdir = DATASET_DIR
    try:
        return sc.datasets.pbmc3k()
    except Exception as exc:
        raise RuntimeError(
            f"the PBMC 3k demo dataset could not be loaded ({type(exc).__name__}: {exc}). It is "
            f"downloaded on first use rather than shipped in the image; this host may have no "
            f"egress, or {sc.settings.datasetdir} may not be writable. Set SC_DATASET_DIR to a "
            f"directory that already holds it, or pass h5ad_path instead of demo=true."
        ) from exc


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
            try:
                adata = load_demo_dataset()
            except RuntimeError as exc:
                # 503, not 500: the service is fine, the dataset is not available here.
                raise HTTPException(503, str(exc)) from exc
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
