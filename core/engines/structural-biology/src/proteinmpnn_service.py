"""ProteinMPNN inverse-folding service (B5). Port 8578. No future-annotations."""
from typing import Optional
import proteinmpnn_design


def create_app():
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel
    app = FastAPI(title="ProteinMPNN Inverse Folding", version="1.0")

    class DesignReq(BaseModel):
        pdb: str                      # PDB text or a path on the box
        num_seq: int = 3
        temperature: float = 0.2
        seed: int = 37

    @app.get("/health")
    def health(): return {"status": "ok"}

    @app.post("/design")
    def design(req: DesignReq):
        try:
            return proteinmpnn_design.design(req.pdb, num_seq=req.num_seq,
                                             temperature=req.temperature, seed=req.seed)
        except Exception as e:
            raise HTTPException(500, f"{type(e).__name__}: {e}")

    return app


def _app_factory(): return create_app()
