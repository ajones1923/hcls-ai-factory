"""
Guided generate-score-reseed loop (C5).

The headline small-molecule gap: a closed optimization loop. Each iteration generates
candidates from the current seeds, scores them (QED + optional real ADMET), filters on hard
constraints, ranks by a composite reward, and reseeds the next iteration from the top-k.
Reuses the existing molecule generator and ADMET scorer — this is orchestration, not a new
model. Every run is logged to single-box MLOps with a reward trajectory.

Generator and scorer are injectable, so the loop is unit-testable without models, and the
real run uses BRICS/SAFE generation + the live ADMET service.
"""
from __future__ import annotations

from typing import Any, Callable

from rdkit import Chem
from rdkit.Chem import QED
from rdkit import RDLogger

from molecule_generator_v2 import MoleculeGenerator

RDLogger.DisableLog("rdApp.*")


class GuidedDesignLoop:
    def __init__(
        self,
        generator: Any = None,
        admet_scorer: Callable[[str], dict] | None = None,
        mlops: Any = None,
    ) -> None:
        self.generator = generator or MoleculeGenerator()
        self.admet_scorer = admet_scorer       # callable(smiles)->{verdict,flags,...}; optional
        self.mlops = mlops                      # MLOpsStore; optional

    def score(self, smiles: str) -> dict[str, Any]:
        mol = Chem.MolFromSmiles(smiles)
        qed = round(QED.qed(mol), 3) if mol else 0.0
        verdict, n_flags = None, 0
        if self.admet_scorer:
            a = self.admet_scorer(smiles) or {}
            verdict, n_flags = a.get("verdict"), len(a.get("flags", []))
        reward = qed / (1 + n_flags)            # composite: drug-likeness penalized by ADMET flags
        return {"smiles": smiles, "qed": qed, "admet_verdict": verdict,
                "n_flags": n_flags, "reward": round(reward, 4)}

    def run(self, seeds: list[str], iterations: int = 5, per_iter: int = 12,
            top_k: int = 3, min_qed: float = 0.5) -> dict[str, Any]:
        run_id = None
        if self.mlops:
            run_id = self.mlops.start_run("guided-design", capability="molecule-generator",
                                          params={"iterations": iterations, "seeds": seeds})
            self.mlops.set_status(run_id, "running")
        current = list(seeds)
        explored: dict[str, dict] = {}
        history = []
        for it in range(iterations):
            gen = self.generator.generate(current, n=per_iter)
            scored = [self.score(m["smiles"]) for m in gen.get("molecules", [])]
            for s in scored:
                explored[s["smiles"]] = s
            passed = [s for s in scored if s["qed"] >= min_qed and s["admet_verdict"] != "high-risk"]
            passed.sort(key=lambda s: -s["reward"])
            best = passed[0]["reward"] if passed else (max((s["reward"] for s in scored), default=0.0))
            history.append({"iteration": it, "generated": len(scored),
                            "passed": len(passed), "best_reward": best})
            if self.mlops and run_id:
                self.mlops.log_metric(run_id, "best_reward", best, step=it)
            current = [s["smiles"] for s in passed[:top_k]] or current   # reseed; keep seeds if none pass
        ranked = sorted(explored.values(), key=lambda s: -s["reward"])
        if self.mlops and run_id:
            self.mlops.set_status(run_id, "complete")
        return {
            "iterations": iterations,
            "n_explored": len(explored),
            "reward_trajectory": [h["best_reward"] for h in history],
            "history": history,
            "top": ranked[:top_k],
            "run_id": run_id,
        }
