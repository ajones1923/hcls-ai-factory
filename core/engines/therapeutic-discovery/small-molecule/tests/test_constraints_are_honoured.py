"""Regression tests: the generator must not silently answer a different question.

Audit compute#43 / C18. `MoleculeGenerator.generate(seeds, n, scaffold)` is the unified entry
point, and capability `molecule-generator` is registered `status: live, maturity: verified` with
`seeds` among its inputs.

  * When the learned SAFE backend was active and no `scaffold` was given, `generate` called
    `d.de_novo_generation(n_samples_per_trial=n, sanitize=True)` -- which takes no seed argument.
    The seeds were discarded and the molecules sampled from the model's prior bore no relation to
    the requested chemistry. `POST /generate {"seeds": [CB-5083], "n": 10}` returned ten
    unrelated molecules under `{"backend": "safe", "n": 10}`, and `GuidedDesignLoop.run` reseeds
    the next iteration from them -- so the whole optimization drifted off the target series while
    reporting a rising reward trajectory.
  * Conversely the BRICS fallback accepts `scaffold` at the API boundary and ignores it
    entirely, returning unconstrained fragment recombinations labelled `{"backend": "brics"}`
    with nothing recording that the constraint had been dropped.
"""
import sys
from pathlib import Path

import pytest
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from molecule_generator_v2 import MoleculeGenerator, SafeGenerator  # noqa: E402

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
IBUPROFEN = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
CB5083 = "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5"


class _FakeDesigner:
    """Stands in for safe.SAFEDesign; records which call the generator made."""

    def __init__(self):
        self.calls = []

    def de_novo_generation(self, n_samples_per_trial=10, **kwargs):
        self.calls.append(("de_novo", kwargs))
        return [ASPIRIN] * n_samples_per_trial

    def scaffold_decoration(self, scaffold=None, n_samples_per_trial=10, **kwargs):
        self.calls.append(("scaffold_decoration", scaffold))
        return [IBUPROFEN] * n_samples_per_trial


@pytest.fixture
def safe_gen():
    gen = SafeGenerator()
    gen._designer = _FakeDesigner()
    return gen


class TestSeedsReachTheSafeBackend:
    def test_a_seeded_request_decorates_a_seed_derived_scaffold(self, safe_gen):
        out = safe_gen.generate([CB5083], n=3)
        kind, scaffold = safe_gen._designer.calls[0]
        assert kind == "scaffold_decoration", (
            "de_novo_generation takes no seeds; calling it discards them"
        )
        assert Chem.MolFromSmiles(scaffold) is not None
        assert all(m["seeds_honoured"] for m in out)
        assert all(m["scaffold_applied"] == scaffold for m in out)

    def test_an_explicit_scaffold_still_wins(self, safe_gen):
        safe_gen.generate([CB5083], n=2, scaffold="c1ccccc1")
        assert safe_gen._designer.calls[0] == ("scaffold_decoration", "c1ccccc1")

    def test_de_novo_is_only_used_when_nothing_was_requested(self, safe_gen):
        safe_gen.generate([], n=2)
        assert safe_gen._designer.calls[0][0] == "de_novo"

    def test_unusable_seeds_are_refused_not_ignored(self, safe_gen):
        with pytest.raises(ValueError, match="Refusing"):
            safe_gen.generate(["not_a_smiles"], n=5)
        assert safe_gen._designer.calls == [], "nothing may be generated for a refused request"

    def test_a_seed_with_no_ring_system_is_refused(self, safe_gen):
        # Murcko scaffold of an acyclic molecule is empty: there is nothing to decorate.
        with pytest.raises(ValueError):
            safe_gen.generate(["CCCCO"], n=5)


class TestDroppedConstraintsAreReported:
    def test_brics_reports_that_it_dropped_the_scaffold(self):
        res = MoleculeGenerator(prefer="brics").generate(
            [ASPIRIN, IBUPROFEN], n=4, scaffold="c1ccccc1")
        assert res["backend"] == "brics"
        assert res["constraints_dropped"] == ["scaffold"]
        assert res["scaffold_applied"] is None
        assert "cannot constrain" in res["constraints_dropped_detail"]

    def test_brics_without_a_scaffold_drops_nothing(self):
        res = MoleculeGenerator(prefer="brics").generate([ASPIRIN, IBUPROFEN], n=4)
        assert res["constraints_dropped"] == []
        assert res["seeds_honoured"] is True

    def test_every_response_declares_what_was_honoured(self):
        res = MoleculeGenerator(prefer="brics").generate([ASPIRIN, IBUPROFEN], n=2)
        for key in ("scaffold_applied", "seeds_honoured", "constraints_dropped"):
            assert key in res
