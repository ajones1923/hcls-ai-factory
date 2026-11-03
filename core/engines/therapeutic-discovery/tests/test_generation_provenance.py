"""Regression tests: MoleculeGenerator must not pass a hand-written library off as a model.

`_generate_bionemo()` used to print "BioNeMo generation for ..." and then silently return
`_generate_fallback()` -- a nine-entry hand-written analogue list. `_get_demo_molecules()`
additionally synthesized molecular properties from the list index (`450 + i * 10` Da,
`3.5 + i * 0.1` logP), numbers that described no molecule.
"""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.molecule_generator import MoleculeGenerator  # noqa: E402

SEED = "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5"

try:
    from rdkit import Chem  # noqa: F401
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@pytest.fixture
def gen(tmp_path):
    return MoleculeGenerator(use_bionemo=False, output_dir=tmp_path)


class TestBioNeMoIsNotFaked:
    def test_generate_bionemo_refuses_instead_of_falling_back(self, gen):
        with pytest.raises(NotImplementedError) as e:
            gen._generate_bionemo(SEED, "VCP", 5, 0.3)
        # The message must point at the real client, not just fail opaquely.
        assert "nim_clients" in str(e.value)

    def test_require_real_refuses_the_library_path(self, tmp_path):
        g = MoleculeGenerator(use_bionemo=False, output_dir=tmp_path, require_real=True)
        with pytest.raises(RuntimeError) as e:
            g.generate_from_seed(SEED, "VCP", num_molecules=5)
        assert "hand-written" in str(e.value)

    def test_fallback_still_available_when_not_required(self, gen):
        mols = gen.generate_from_seed(SEED, "VCP", num_molecules=5)
        assert mols, "the labeled fallback should still work when real is not required"


class TestProvenanceIsStamped:
    def test_library_molecules_declare_their_origin(self, gen):
        mols = gen._get_demo_molecules("VCP", 4)
        assert mols
        for m in mols:
            assert m.properties["provenance"] == "hand_written_analogue_library"
            assert m.properties["score_is_placeholder"] is True
            assert m.generation_method == "demo_library"

    def test_unknown_target_yields_nothing_rather_than_invented_molecules(self, gen):
        assert gen._get_demo_molecules("NOT_A_GENE", 10) == []

    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    def test_properties_are_computed_not_derived_from_the_index(self, gen):
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mols = gen._get_demo_molecules("VCP", 5)
        for i, m in enumerate(mols):
            expected_mw = Descriptors.MolWt(Chem.MolFromSmiles(m.smiles))
            assert abs(m.properties["molecular_weight"] - expected_mw) < 0.05
            # The old index formula; none of these molecules weigh 450 + i*10 Da.
            assert m.properties["molecular_weight"] != 450 + i * 10
            assert m.properties["logP"] != 3.5 + i * 0.1

    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    def test_seed_and_analogues_are_labeled_differently(self, gen):
        mols = gen._generate_rdkit_analogues(SEED, "VCP", 4)
        assert mols[0].generation_method == "seed"
        assert mols[0].properties["provenance"] == "input_seed"
        for m in mols[1:]:
            assert m.properties["provenance"] == "hand_written_analogue_library"
