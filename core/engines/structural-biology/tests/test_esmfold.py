"""Tests for the ESMFold service (B1). No model download needed — folder is injected."""
import pytest

from esmfold_service import (
    ESMFoldModel, FoldResult, SequenceError, parse_pdb_stats, validate_sequence,
)
from client import MockESMFoldClient, create_esmfold_client

# a minimal real-shaped PDB snippet (2 residues, with pLDDT in the B-factor column)
SAMPLE_PDB = (
    "ATOM      1  N   MET A   1      0.000   0.000   0.000  1.00 85.20           N\n"
    "ATOM      2  CA  MET A   1      1.458   0.000   0.000  1.00 88.40           C\n"
    "ATOM      3  N   ALA A   2      2.000   1.400   0.000  1.00 72.10           N\n"
    "ATOM      4  CA  ALA A   2      3.400   1.800   0.000  1.00 70.00           C\n"
    "END\n"
)


class TestValidate:
    def test_normalizes_and_uppercases(self):
        assert validate_sequence(" mkt v\nl ") == "MKTVL"

    def test_rejects_empty(self):
        with pytest.raises(SequenceError, match="empty"):
            validate_sequence("   ")

    def test_rejects_non_amino_acids(self):
        with pytest.raises(SequenceError, match="invalid residues"):
            validate_sequence("MKTXZ1")

    def test_rejects_too_long(self):
        with pytest.raises(SequenceError, match="exceeds"):
            validate_sequence("A" * 10, max_len=5)


class TestPdbStats:
    def test_counts_atoms_residues_and_plddt(self):
        s = parse_pdb_stats(SAMPLE_PDB)
        assert s["n_atoms"] == 4
        assert s["n_residues"] == 2
        assert s["mean_plddt"] == pytest.approx((85.2 + 88.4 + 72.1 + 70.0) / 4, abs=0.01)

    def test_plddt_normalized_from_0_1_scale(self):
        # transformers ESMFold writes pLDDT in [0,1]; parse_pdb_stats must report 0-100.
        # Reuse SAMPLE_PDB's exact column alignment, just swap the B-factor values.
        pdb01 = (SAMPLE_PDB.replace("85.20", " 0.85").replace("88.40", " 0.88")
                 .replace("72.10", " 0.72").replace("70.00", " 0.70"))
        s = parse_pdb_stats(pdb01)
        assert s["mean_plddt"] == pytest.approx((85 + 88 + 72 + 70) / 4, abs=0.5)   # ×100


class TestModelInjected:
    def test_fold_uses_injected_folder(self):
        m = ESMFoldModel(_folder=lambda seq: SAMPLE_PDB)
        res = m.fold("MA")
        assert isinstance(res, FoldResult)
        assert res.served_by == "test"
        assert res.n_residues == 2 and res.n_atoms == 4
        assert res.mean_plddt > 70

    def test_fold_validates_before_folding(self):
        m = ESMFoldModel(_folder=lambda seq: SAMPLE_PDB)
        with pytest.raises(SequenceError):
            m.fold("not-a-sequence-123")

    def test_to_dict_can_drop_pdb(self):
        m = ESMFoldModel(_folder=lambda seq: SAMPLE_PDB)
        d = m.fold("MA").to_dict(include_pdb=False)
        assert "pdb" not in d and d["n_atoms"] == 4


class TestClient:
    def test_mock_is_labeled(self):
        out = MockESMFoldClient().fold("MKTV")
        assert out["served_by"] == "mock" and "_warning" in out
        assert out["n_residues"] == 4

    def test_factory_mock_mode(self):
        assert isinstance(create_esmfold_client("mock"), MockESMFoldClient)
