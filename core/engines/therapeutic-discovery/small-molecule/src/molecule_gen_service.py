"""Molecule generation FastAPI service + client (C1). Port 8574. No future-annotations (FastAPI body fix)."""
from typing import Any, Optional, List
from molecule_generator_v2 import MoleculeGenerator


def create_app(generator: Optional[MoleculeGenerator] = None):
    from fastapi import FastAPI
    from pydantic import BaseModel
    gen = generator or MoleculeGenerator()
    app = FastAPI(title="Molecule Generator Service", version="1.0")

    class GenReq(BaseModel):
        seeds: List[str]
        n: int = 10
        scaffold: Optional[str] = None

    @app.get("/health")
    def health():
        return {"status": "ok", "backend": gen.backend()}

    @app.post("/generate")
    def generate(req: GenReq):
        return gen.generate(req.seeds, n=req.n, scaffold=req.scaffold)

    return app


def _app_factory():
    return create_app()


class MoleculeGenClient:
    def __init__(self, endpoint: str = "http://localhost:8574", timeout: float = 300.0):
        self.endpoint = endpoint.rstrip("/"); self.timeout = timeout

    def generate(self, seeds, n: int = 10, scaffold=None) -> dict:
        import httpx
        r = httpx.post(f"{self.endpoint}/generate", json={"seeds": seeds, "n": n, "scaffold": scaffold}, timeout=self.timeout)
        r.raise_for_status(); return r.json()
