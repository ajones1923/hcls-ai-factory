"""ADMET prediction FastAPI service + client (C2). Default port 8572."""

from typing import Any

from admet_predictor import ADMETPredictor, SmilesError


def create_app(predictor: ADMETPredictor | None = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel

    pred = predictor or ADMETPredictor()
    app = FastAPI(title="ADMET Service", version="1.0")

    class AdmetReq(BaseModel):
        smiles: str
        full: bool = False

    @app.get("/health")
    def health():
        return {"status": "ok", "loaded": pred._loaded is not None}

    @app.post("/admet")
    def admet(req: AdmetReq):
        try:
            res = pred.predict(req.smiles)
        except SmilesError as e:
            raise HTTPException(422, str(e))
        return res.to_dict(full=req.full)

    return app


def _app_factory():
    return create_app()


class ADMETClient:
    """Calls a running ADMET service."""

    def __init__(self, endpoint: str = "http://localhost:8572", timeout: float = 120.0) -> None:
        self.endpoint = endpoint.rstrip("/")
        self.timeout = timeout

    def predict(self, smiles: str, full: bool = False) -> dict[str, Any]:
        import httpx
        r = httpx.post(f"{self.endpoint}/admet", json={"smiles": smiles, "full": full}, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def health(self) -> dict[str, Any]:
        import httpx
        return httpx.get(f"{self.endpoint}/health", timeout=10).json()
