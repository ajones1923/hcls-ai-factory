"""B5: ProteinMPNN inverse folding. Parse test (no subprocess) + a real-run test."""
import os
from pathlib import Path
import pytest
from proteinmpnn_design import _parse_fasta, design, MPNN_DIR

SAMPLE_FA = """>5L33, score=1.5720, global_score=1.5720, designed_chains=['A'], model_name=v_48_020
HMPEEEKAARLFIEALEKGD
>T=0.2, sample=1, score=0.9732, global_score=0.9732, seq_recovery=0.4340
PVDADTRVALDFIRALEAAD
>T=0.2, sample=2, score=0.8113, global_score=0.8113, seq_recovery=0.4623
SINEEEKKALDFVEALEKAD
"""


def test_parse_fasta(tmp_path):
    fa = tmp_path / "x.fa"; fa.write_text(SAMPLE_FA)
    native, designs = _parse_fasta(fa)
    assert native["native"] and native["sequence"].startswith("HMPEEE")
    assert len(designs) == 2
    assert designs[0]["seq_recovery"] == 0.434 and designs[1]["sample"] == 2


@pytest.mark.skipif(not (MPNN_DIR / "protein_mpnn_run.py").exists(), reason="ProteinMPNN not vendored")
def test_design_real():
    pdb = sorted((MPNN_DIR / "inputs/PDB_monomers/pdbs").glob("*.pdb"))[0]
    res = design(str(pdb), num_seq=2, temperature=0.2, seed=37)
    assert res["n_designs"] == 2
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    for d in res["designs"]:
        assert d["sequence"] and set(d["sequence"]) <= aa     # valid amino acids
        assert 0.0 < d["seq_recovery"] <= 1.0
    assert res["mean_seq_recovery"] > 0.2                      # ProteinMPNN recovers native ~40-50%
