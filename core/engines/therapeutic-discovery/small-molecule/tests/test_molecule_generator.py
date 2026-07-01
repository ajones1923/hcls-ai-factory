"""Tests for C1 molecule generation. BRICS backend is real, deterministic, no model needed."""
from rdkit import Chem
from molecule_generator_v2 import BRICSGenerator, MoleculeGenerator, canonical

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
IBUPROFEN = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"


class TestBRICS:
    def test_generates_valid_novel_molecules(self):
        out = BRICSGenerator().generate([ASPIRIN, IBUPROFEN], n=5)
        assert out, "BRICS should produce molecules"
        for m in out:
            assert Chem.MolFromSmiles(m["smiles"]) is not None      # valid
            assert 0.0 <= m["qed"] <= 1.0 and m["source"] == "brics"

    def test_excludes_the_seeds_themselves(self):
        seeds = [ASPIRIN, IBUPROFEN]
        seed_canon = {canonical(s) for s in seeds}
        out = BRICSGenerator().generate(seeds, n=10)
        assert all(m["smiles"] not in seed_canon for m in out)

    def test_results_are_unique(self):
        out = BRICSGenerator().generate([ASPIRIN, IBUPROFEN], n=10)
        smis = [m["smiles"] for m in out]
        assert len(smis) == len(set(smis))

    def test_empty_seed_returns_empty(self):
        assert BRICSGenerator().generate(["not_a_smiles"], n=5) == []


class TestUnified:
    def test_falls_back_to_brics_and_reports_backend(self):
        g = MoleculeGenerator(prefer="brics")
        res = g.generate([ASPIRIN, IBUPROFEN], n=4)
        assert res["backend"] == "brics" and res["n"] >= 1
        assert res["molecules"][0]["source"] == "brics"
