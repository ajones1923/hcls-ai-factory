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

## Roadmap (this engine)
ESMFold (B1) ✓ · ESM-2 embeddings + sequence search (B2) · Boltz-1 complexes · ProteinMPNN
design (B3) · developability suite (B4).
