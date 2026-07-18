"""
ESMFold protein structure prediction (B1) — the entry point of the Large-Molecule engine.

Sequence -> 3D structure (PDB) using ESMFold (facebook/esmfold_v1) via transformers.
Real model, real output. GPU acceleration (ESMFold NIM or aarch64-CUDA torch) is a
drop-in swap behind the same interface; CPU is the default fallback on this box.

The module separates pure, testable logic (validation, PDB stats) from the heavy model
load, and accepts an injectable folder so the service can be unit-tested without the
2.8 GB model.
"""

import re
from dataclasses import dataclass, asdict
from typing import Any, Callable

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
# GPU-served on the GB10 → VCP-scale (806 aa) folds in ~1s; 1024 covers most single-chain proteins.
DEFAULT_MAX_LEN = 1024


class SequenceError(ValueError):
    """Raised for an invalid protein sequence."""


# --------------------------------------------------------------------------- #
# Pure logic (no model needed) — fully unit-testable
# --------------------------------------------------------------------------- #
def validate_sequence(seq: str, max_len: int = DEFAULT_MAX_LEN) -> str:
    """Normalize + validate a one-letter amino-acid sequence."""
    s = re.sub(r"\s+", "", (seq or "")).upper()
    if not s:
        raise SequenceError("empty sequence")
    bad = sorted(set(s) - AMINO_ACIDS)
    if bad:
        raise SequenceError(f"invalid residues: {bad} (expected 20 standard amino acids)")
    if len(s) > max_len:
        raise SequenceError(
            f"sequence length {len(s)} exceeds CPU limit {max_len}; fold on GPU (NIM) for larger proteins"
        )
    return s


def parse_pdb_stats(pdb: str) -> dict[str, Any]:
    """Extract structure stats from a PDB string. ESMFold writes per-residue pLDDT
    into the B-factor column, so mean B-factor == mean pLDDT confidence (0-100)."""
    atoms = [l for l in pdb.splitlines() if l.startswith("ATOM")]
    residues = {(l[21:22], l[22:26].strip()) for l in atoms}
    bfactors = []
    for l in atoms:
        col = l[60:66].strip()
        if col:
            try:
                bfactors.append(float(col))
            except ValueError:
                pass
    # transformers' ESMFold writes per-residue pLDDT in [0,1]; the AlphaFold/ESMFold
    # convention (and this function's contract) is 0-100 — normalize so a confident
    # fold reads 86, not 0.86. Values already on the 0-100 scale are left untouched.
    if bfactors and max(bfactors) <= 1.0:
        bfactors = [b * 100.0 for b in bfactors]
    plddt = round(sum(bfactors) / len(bfactors), 2) if bfactors else None
    return {"n_atoms": len(atoms), "n_residues": len(residues), "mean_plddt": plddt}


@dataclass
class FoldResult:
    sequence: str
    length: int
    pdb: str
    n_atoms: int
    n_residues: int
    mean_plddt: float | None
    served_by: str          # "esmfold" (real) | "mock" (labeled) | "test"

    def to_dict(self, include_pdb: bool = True) -> dict[str, Any]:
        d = asdict(self)
        if not include_pdb:
            d.pop("pdb")
        return d


# --------------------------------------------------------------------------- #
# Model wrapper (lazy heavy load; injectable for tests)
# --------------------------------------------------------------------------- #
class ESMFoldModel:
    def __init__(
        self,
        model_name: str = "facebook/esmfold_v1",
        max_len: int = DEFAULT_MAX_LEN,
        _folder: Callable[[str], str] | None = None,
    ) -> None:
        self.model_name = model_name
        self.max_len = max_len
        self._folder = _folder        # callable(seq)->pdb for tests; None => real model
        self._model = None
        self._torch = None
        self.device = "cpu"

    def _load(self) -> None:
        if self._model is None:
            import torch
            from transformers import EsmForProteinFolding
            self._torch = torch
            self.device = "cuda" if torch.cuda.is_available() else "cpu"   # GB10 when available
            m = EsmForProteinFolding.from_pretrained(self.model_name, low_cpu_mem_usage=True).eval()
            if self.device == "cuda":
                m = m.cuda()
                m.esm = m.esm.half()              # ESMFold's recommended GPU setup (half-precision trunk)
            self._model = m

    def fold(self, sequence: str) -> FoldResult:
        seq = validate_sequence(sequence, self.max_len)
        if self._folder is not None:                 # test path
            pdb, served = self._folder(seq), "test"
        else:                                        # real path
            self._load()
            with self._torch.no_grad():
                pdb = self._model.infer_pdb(seq)
            served = f"esmfold-{self.device}"
        stats = parse_pdb_stats(pdb)
        return FoldResult(
            sequence=seq, length=len(seq), pdb=pdb,
            n_atoms=stats["n_atoms"], n_residues=stats["n_residues"],
            mean_plddt=stats["mean_plddt"], served_by=served,
        )


# --------------------------------------------------------------------------- #
# FastAPI service
# --------------------------------------------------------------------------- #
def create_app(model: ESMFoldModel | None = None):
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel

    mdl = model or ESMFoldModel()
    app = FastAPI(title="ESMFold Service", version="1.0")

    class FoldRequest(BaseModel):
        sequence: str
        include_pdb: bool = True

    @app.get("/health")
    def health():
        return {"status": "ok", "model": mdl.model_name, "loaded": mdl._model is not None}

    @app.post("/fold")
    def fold(req: FoldRequest):
        try:
            res = mdl.fold(req.sequence)
        except SequenceError as e:
            raise HTTPException(status_code=422, detail=str(e))
        return res.to_dict(include_pdb=req.include_pdb)

    return app


# uvicorn entrypoint: `uvicorn esmfold_service:app` (lazy-loads the model on first /fold)
def _app_factory():
    return create_app()
