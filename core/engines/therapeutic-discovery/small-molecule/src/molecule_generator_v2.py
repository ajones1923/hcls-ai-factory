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
        """Generate from a scaffold, or from a seed-derived scaffold. Never ignores both.

        With a seed list and no scaffold this used to call
        `d.de_novo_generation(n_samples_per_trial=n)`, which takes no seed argument: the seeds
        were discarded and the molecules sampled from the model's prior bore no relation to the
        requested chemistry. The capability is registered `live`/`verified` with `seeds` among
        its inputs, so `POST /generate {"seeds": [CB-5083], "n": 10}` returned ten unrelated
        de-novo molecules under `{"backend": "safe", "n": 10}` -- and `GuidedDesignLoop.run`
        reseeds the next iteration from them, so the optimization drifted off the target series
        while reporting a rising reward trajectory.

        A seeded request is now honoured by decorating the seed's Murcko scaffold, and refused
        if that cannot be done, rather than silently answering a different question.
        """
        d = self._load()
        applied_scaffold = scaffold
        if not applied_scaffold and seeds:
            applied_scaffold = self._seed_scaffold(seeds)
            if not applied_scaffold:
                raise ValueError(
                    "SAFE generation was asked to honour seeds but no usable scaffold could be "
                    f"derived from them ({seeds[:3]}). Refusing to sample from the model's prior "
                    "and return molecules unrelated to the seeds; pass an explicit `scaffold`, "
                    "or call the BRICS backend, which recombines fragments of the seeds."
                )

        if applied_scaffold:
            smis = d.scaffold_decoration(scaffold=applied_scaffold, n_samples_per_trial=n,
                                         n_trials=1, sanitize=True,
                                         do_not_fragment_further=True)
        else:
            # No seeds and no scaffold: de novo is what was asked for.
            smis = d.de_novo_generation(n_samples_per_trial=n, sanitize=True)

        out, seen = [], set()
        for s in smis:
            c = canonical(s)
            if not c or c in seen:
                continue
            seen.add(c)
            out.append({
                "smiles": c,
                "qed": round(QED.qed(Chem.MolFromSmiles(c)), 3),
                "source": self.name,
                "scaffold_applied": applied_scaffold,
                "seeds_honoured": bool(applied_scaffold and seeds),
            })
        return out[:n]

    @staticmethod
    def _seed_scaffold(seeds: list[str]) -> str | None:
        """The Bemis-Murcko scaffold of the first parseable seed, or None."""
        from rdkit.Chem.Scaffolds import MurckoScaffold

        for seed in seeds:
            mol = Chem.MolFromSmiles(seed or "")
            if mol is None:
                continue
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            except Exception:
                continue
            if scaffold is None or scaffold.GetNumAtoms() == 0:
                continue
            return Chem.MolToSmiles(scaffold)
        return None


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
        """Generate molecules, and state which constraints were actually honoured.

        The BRICS fallback accepts `scaffold` at the API boundary and ignores it entirely, so a
        scaffold-constrained request returned unconstrained fragment recombinations labelled
        `{"backend": "brics"}` with nothing recording that the constraint had been dropped --
        which is how the registry's "scaffold-constrained generation" became unreachable in
        practice. Every response now reports `scaffold_applied` and `seeds_honoured`, and a
        dropped scaffold is named in `constraints_dropped`.
        """
        if self._safe:
            try:
                mols = self._safe.generate(seeds, n=n, scaffold=scaffold)
                if mols:
                    return {"backend": "safe", "n": len(mols), "molecules": mols,
                            "scaffold_applied": mols[0].get("scaffold_applied"),
                            "seeds_honoured": mols[0].get("seeds_honoured", False),
                            "constraints_dropped": []}
            except Exception as e:
                safe_error = f"{type(e).__name__}: {e}"
            else:
                safe_error = "the SAFE backend returned no molecules"
        else:
            safe_error = None

        mols = self.brics.generate(seeds, n=n)
        # BRICS decomposes and recombines the seeds, so seeds ARE honoured here; a scaffold
        # constraint is not supported and is reported as dropped rather than ignored.
        dropped = ["scaffold"] if scaffold else []
        result = {"backend": "brics", "n": len(mols), "molecules": mols,
                  "scaffold_applied": None,
                  "seeds_honoured": bool(seeds),
                  "constraints_dropped": dropped}
        if dropped:
            result["constraints_dropped_detail"] = (
                "the BRICS backend cannot constrain generation to a scaffold; the molecules "
                "below are unconstrained recombinations of the seeds' fragments"
            )
        if safe_error:
            result["degraded_from"] = f"safe ({safe_error})"
        return result
