"""
ESMFold client — follows the drug-discovery `nim_clients` local/mock pattern so the
Large-Molecule engine plugs into the factory the same way MolMIM/DiffDock do.

Honesty: the mock client is for offline/demo only and labels every result
`served_by="mock"`; it must never be presented as a real prediction.
"""
from __future__ import annotations

import os
from typing import Any


class ESMFoldClient:
    """Calls a running ESMFold service (real model)."""

    def __init__(self, endpoint: str = "http://localhost:8570", timeout: float = 600.0) -> None:
        self.endpoint = endpoint.rstrip("/")
        self.timeout = timeout

    def fold(self, sequence: str, include_pdb: bool = True) -> dict[str, Any]:
        import httpx
        r = httpx.post(
            f"{self.endpoint}/fold",
            json={"sequence": sequence, "include_pdb": include_pdb},
            timeout=self.timeout,
        )
        r.raise_for_status()
        return r.json()

    def health(self) -> dict[str, Any]:
        import httpx
        return httpx.get(f"{self.endpoint}/health", timeout=10).json()


class MockESMFoldClient:
    """Offline placeholder — clearly labeled, NEVER a real structure."""

    def fold(self, sequence: str, include_pdb: bool = True) -> dict[str, Any]:
        seq = "".join(c for c in sequence.upper() if c.isalpha())
        # a single CA-only placeholder backbone so downstream code has *something* shaped
        # like a PDB, but it is explicitly mock and carries no real geometry.
        pdb = "".join(
            f"ATOM  {i+1:>5}  CA  {'GLY'} A{i+1:>4}    "
            f"{i*3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
            for i in range(len(seq))
        ) + "END\n"
        return {
            "sequence": seq, "length": len(seq),
            "pdb": pdb if include_pdb else None,
            "n_atoms": len(seq), "n_residues": len(seq), "mean_plddt": None,
            "served_by": "mock",
            "_warning": "MOCK structure — not a real prediction. For offline/demo only.",
        }


def create_esmfold_client(mode: str | None = None):
    """Factory mirroring NIMServiceManager: ESMFOLD_MODE = local | mock (default local)."""
    mode = (mode or os.getenv("ESMFOLD_MODE", "local")).lower()
    if mode == "mock":
        return MockESMFoldClient()
    return ESMFoldClient(endpoint=os.getenv("ESMFOLD_ENDPOINT", "http://localhost:8570"))


class ProteinSearchClient:
    """Calls a running protein embedding + search service (B2)."""

    def __init__(self, endpoint: str = "http://localhost:8571", timeout: float = 300.0) -> None:
        self.endpoint = endpoint.rstrip("/")
        self.timeout = timeout

    def _post(self, path: str, body: dict) -> Any:
        import httpx
        r = httpx.post(f"{self.endpoint}{path}", json=body, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def embed(self, sequence: str) -> dict[str, Any]:
        return self._post("/embed", {"sequence": sequence})

    def index(self, proteins: list[dict]) -> dict[str, Any]:
        return self._post("/index", {"proteins": proteins})

    def search(self, sequence: str, top_k: int = 5) -> dict[str, Any]:
        return self._post("/search", {"sequence": sequence, "top_k": top_k})

    def health(self) -> dict[str, Any]:
        import httpx
        return httpx.get(f"{self.endpoint}/health", timeout=10).json()
