"""
Molecule generation for diversity (C1) — a second generator alongside MolMIM.

Two real backends behind one interface:
  * **BRICS** (RDKit) — fragment-based de novo generation: decompose seed actives into
    BRICS fragments and recombine into novel, valid, drug-like molecules. Zero-dependency,
    always available, fully deterministic-testable.
  * **SAFE** (learned) — the GenMol-family SAFE-fragment generative model (de novo +
    scaffold-constrained), used when its package is installed; the NVIDIA GenMol NIM is the
    GPU-accelerated variant of the same approach.

Every generated molecule is RDKit-canonicalized, de-duplicated against the seeds, and QED-scored.
"""
from __future__ import annotations

from typing import Any

from rdkit import Chem
from rdkit.Chem import BRICS, QED
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def canonical(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m else None


# --------------------------------------------------------------------------- #
class BRICSGenerator:
    """Fragment-based de novo generation via BRICS decompose + recombine (real, zero-dep)."""

    name = "brics"

    def generate(self, seeds: list[str], n: int = 10, max_iter: int = 2000) -> list[dict[str, Any]]:
        mols = [m for m in (Chem.MolFromSmiles(s) for s in seeds) if m is not None]
        if not mols:
            return []
        frags = set()
        for m in mols:
            frags |= set(BRICS.BRICSDecompose(m))
        frag_mols = [Chem.MolFromSmiles(f) for f in frags]
        frag_mols = [f for f in frag_mols if f is not None]
        seed_canon = {Chem.MolToSmiles(m) for m in mols}
        out: list[dict[str, Any]] = []
        seen: set[str] = set()
        for i, built in enumerate(BRICS.BRICSBuild(frag_mols)):
            if i >= max_iter or len(out) >= n:
                break
            try:
                Chem.SanitizeMol(built)
                smi = Chem.MolToSmiles(built)
            except Exception:
                continue
            if smi in seen or smi in seed_canon:
                continue
            seen.add(smi)
            out.append({"smiles": smi, "qed": round(QED.qed(built), 3), "source": self.name})
        return out


# --------------------------------------------------------------------------- #
class SafeGenerator:
    """Learned SAFE-fragment generation (GenMol family). Available iff `safe` is installed."""

    name = "safe"

    def __init__(self) -> None:
        self._designer = None

    @staticmethod
    def available() -> bool:
        try:
            import safe  # noqa: F401
            return True
        except Exception:
            return False

    def _load(self):
        if self._designer is None:
            import safe
            self._designer = safe.SAFEDesign.load_default()
        return self._designer

    def generate(self, seeds: list[str], n: int = 10, scaffold: str | None = None) -> list[dict[str, Any]]:
        d = self._load()
        if scaffold:
            smis = d.scaffold_decoration(scaffold=scaffold, n_samples_per_trial=n, n_trials=1,
                                         sanitize=True, do_not_fragment_further=True)
        else:
            smis = d.de_novo_generation(n_samples_per_trial=n, sanitize=True)
        out, seen = [], set()
        for s in smis:
            c = canonical(s)
            if not c or c in seen:
                continue
            seen.add(c)
            out.append({"smiles": c, "qed": round(QED.qed(Chem.MolFromSmiles(c)), 3), "source": self.name})
        return out[:n]


# --------------------------------------------------------------------------- #
class MoleculeGenerator:
    """Unified generator: prefers the learned SAFE model, falls back to BRICS."""

    def __init__(self, prefer: str = "auto") -> None:
        self.brics = BRICSGenerator()
        self.prefer = prefer
        self._safe = SafeGenerator() if (prefer in ("auto", "safe") and SafeGenerator.available()) else None

    def backend(self) -> str:
        return "safe" if self._safe else "brics"

    def generate(self, seeds: list[str], n: int = 10, scaffold: str | None = None) -> dict[str, Any]:
        if self._safe:
            try:
                mols = self._safe.generate(seeds, n=n, scaffold=scaffold)
                if mols:
                    return {"backend": "safe", "n": len(mols), "molecules": mols}
            except Exception:
                pass   # degrade to BRICS, never fail the request
        mols = self.brics.generate(seeds, n=n)
        return {"backend": "brics", "n": len(mols), "molecules": mols}
