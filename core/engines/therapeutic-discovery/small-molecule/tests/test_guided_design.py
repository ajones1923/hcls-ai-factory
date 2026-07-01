"""C5 guided generate-score-reseed loop. Generator + scorer injected (deterministic)."""
from guided_design import GuidedDesignLoop


class FakeGen:
    """Returns a fixed candidate set per call (ignores seeds/n)."""
    def __init__(self, mols): self.mols = mols; self.calls = 0
    def generate(self, seeds, n=10, scaffold=None):
        self.calls += 1
        return {"backend": "fake", "n": len(self.mols), "molecules": [{"smiles": s} for s in self.mols]}


def admet(smiles):
    # mark benzene 'high-risk', others clean
    return {"verdict": "high-risk", "flags": ["x", "y", "z"]} if smiles == "c1ccccc1" else {"verdict": "clean", "flags": []}


class TestScore:
    def test_reward_penalized_by_flags(self):
        loop = GuidedDesignLoop(generator=FakeGen([]), admet_scorer=admet)
        clean = loop.score("CCO"); risky = loop.score("c1ccccc1")
        assert clean["reward"] > risky["reward"]            # flags lower the reward
        assert clean["n_flags"] == 0 and risky["n_flags"] == 3


class TestLoop:
    def test_loop_explores_filters_and_reseeds(self):
        gen = FakeGen(["CCO", "CCN", "c1ccccc1", "CC(=O)O"])
        loop = GuidedDesignLoop(generator=gen, admet_scorer=admet)
        out = loop.run(["CCO"], iterations=3, per_iter=4, top_k=2)
        assert gen.calls == 3                               # one generation per iteration
        assert out["n_explored"] >= 3
        assert len(out["reward_trajectory"]) == 3
        # high-risk benzene is filtered out of the top set
        assert all(t["admet_verdict"] != "high-risk" for t in out["top"])

    def test_mlops_logging(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location("_hmlops",
            "./lib/hcls_common/mlops.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        s = m.MLOpsStore(":memory:")
        gen = FakeGen(["CCO", "CCN"])
        out = GuidedDesignLoop(generator=gen, admet_scorer=admet, mlops=s).run(["CCO"], iterations=2)
        run = s.get_run(out["run_id"])
        assert run["status"] == "complete" and "best_reward" in run["metrics"]
