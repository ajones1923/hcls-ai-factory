"""
ProteinMPNN inverse folding (B5).

Structure-conditioned sequence design: given a backbone (PDB), generate amino-acid sequences
predicted to fold into it. Real ProteinMPNN (MIT) — the vendored repo's `protein_mpnn_run.py`
is invoked as a subprocess and its FASTA output parsed into structured designs (with the model's
own score + native-sequence-recovery). Chains naturally to ESMFold (B1) for re-fold validation
and to RFdiffusion backbones (B6, when available).
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

MPNN_DIR = Path(__file__).resolve().parents[1] / "vendor_proteinmpnn"
_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _parse_fasta(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    recs = [(lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]
    native = {"sequence": recs[0][1], "native": True} if recs else {"sequence": "", "native": True}
    designs = []
    for header, seq in recs[1:]:
        fields = {}
        for part in header.lstrip(">").split(","):
            if "=" in part:
                k, _, v = part.strip().partition("=")
                fields[k.strip()] = v.strip()
        designs.append({
            "sequence": seq,
            "score": float(fields.get("score", 0.0)),
            "seq_recovery": float(fields.get("seq_recovery", 0.0)),
            "sample": int(float(fields.get("sample", 0))),
        })
    return native, designs


def design(pdb: str, num_seq: int = 3, temperature: float = 0.2, seed: int = 37,
           python_exe: str | None = None) -> dict[str, Any]:
    """Run ProteinMPNN on a backbone (PDB text or a path) → designed sequences."""
    py = python_exe or os.sys.executable
    runner = MPNN_DIR / "protein_mpnn_run.py"
    if not runner.exists():
        raise FileNotFoundError(f"ProteinMPNN not vendored at {MPNN_DIR}")
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        pdb_path = Path(pdb) if (len(pdb) < 1024 and Path(pdb).exists()) else None
        if pdb_path is None:                       # treat as PDB text
            pdb_path = tdp / "backbone.pdb"
            pdb_path.write_text(pdb)
        out = tdp / "out"
        out.mkdir()
        cmd = [py, str(runner), "--pdb_path", str(pdb_path), "--out_folder", str(out),
               "--num_seq_per_target", str(num_seq), "--sampling_temp", str(temperature),
               "--seed", str(seed), "--batch_size", "1"]
        proc = subprocess.run(cmd, cwd=str(MPNN_DIR), capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            raise RuntimeError(f"ProteinMPNN failed: {proc.stderr[-400:]}")
        fastas = list((out / "seqs").glob("*.fa"))
        if not fastas:
            raise RuntimeError("ProteinMPNN produced no sequences")
        native, designs = _parse_fasta(fastas[0])
        return {
            "backbone_name": fastas[0].stem,
            "native": native,
            "n_designs": len(designs),
            "designs": designs,
            "temperature": temperature,
            "mean_seq_recovery": round(sum(d["seq_recovery"] for d in designs) / max(len(designs), 1), 4),
        }
