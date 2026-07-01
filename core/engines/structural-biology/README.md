# Large-Molecule Engine

Protein structure prediction, embeddings, and design for the HCLS AI Factory.

## Modules
- **ESMFold** (`src/esmfold_service.py`) — sequence → 3D structure (PDB) with per-residue
  pLDDT confidence (B-factor column). Real model (`facebook/esmfold_v1`).

## Run
```bash
pip install -r requirements.txt
uvicorn esmfold_service:_app_factory --factory --port 8570   # from src/
# POST /fold {"sequence": "NLYIQWLKDGGPSSGRPPPS"} -> {pdb, n_atoms, n_residues, mean_plddt}
```

## Serving notes
- Defaults to CPU on this box (native torch is CPU-only). Folds small proteins (≤400 aa).
- GPU acceleration is a drop-in swap behind the same client interface: an ESMFold NIM
  container or an aarch64-CUDA torch build. Larger proteins (e.g. VCP, 806 aa) need GPU.

## Vendored dependencies (bootstrap)
Two third-party design tools are used as external clones under this engine and are
**not** committed here (they carry their own licenses and history). Hydrate them on a
fresh checkout:

```bash
cd core/engines/structural-biology
git clone https://github.com/dauparas/ProteinMPNN.git      vendor_proteinmpnn
git clone https://github.com/RosettaCommons/RFdiffusion.git vendor_rfdiffusion
```

- `vendor_proteinmpnn/` — used by `src/proteinmpnn_design.py` for inverse-folding design
  (the code resolves it at `../vendor_proteinmpnn`).
- `vendor_rfdiffusion/` — de novo backbone generation; runs on remote GPU (RunPod) per the
  factory's on-demand GPU path, not on the local ARM box.

Both directories are gitignored (`vendor_*/`), so they stay local and never publish.

## Roadmap (this engine)
ESMFold (B1) ✓ · ESM-2 embeddings + sequence search (B2) · Boltz-1 complexes · ProteinMPNN
design (B3) · developability suite (B4).
