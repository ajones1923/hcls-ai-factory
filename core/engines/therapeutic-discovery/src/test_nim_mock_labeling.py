"""F2: every mock NIM output must be labeled, and downstream guards must catch it."""
import pytest

from src.nim_clients import (
    MockMolMIMClient, MockDiffDockClient,
    is_mock_result, contains_mock, assert_real,
)

# a minimal protein PDB stub for the mock docker (it only hashes the string)
STUB_PDB = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"


class TestMockLabeling:
    def test_mock_molmim_labels_every_molecule(self):
        mols = MockMolMIMClient().generate("CC(=O)OC1=CC=CC=C1C(=O)O", num_molecules=3)
        assert mols, "mock should return molecules"
        assert all(m["_provenance"] == "mock" for m in mols)
        assert all("SIMULATED" in m["_warning"] for m in mols)

    def test_mock_diffdock_labels_every_pose(self):
        poses = MockDiffDockClient().dock(protein_pdb=STUB_PDB, ligand_smiles="CCO", num_poses=4)
        assert poses and all(p["_provenance"] == "mock" for p in poses)


class TestHonestyGuards:
    def test_is_mock_result(self):
        assert is_mock_result({"smiles": "CCO", "_provenance": "mock"})
        assert not is_mock_result({"smiles": "CCO"})            # real-shaped, unlabeled
        assert not is_mock_result("not a dict")

    def test_contains_mock_over_a_list(self):
        real = [{"smiles": "CCO"}, {"smiles": "CCN"}]
        mixed = real + [{"smiles": "CCC", "_provenance": "mock"}]
        assert not contains_mock(real)
        assert contains_mock(mixed)

    def test_assert_real_raises_on_mock(self):
        with pytest.raises(RuntimeError, match="MOCK"):
            assert_real([{"_provenance": "mock"}], context="final report")

    def test_assert_real_passes_on_real(self):
        assert_real([{"smiles": "CCO", "score": 0.9}])           # no raise
