"""
Molecule Generation Module - BioNeMo integration for drug-like molecule generation.

Uses MolMIM or similar models to generate candidate molecules
from seed SMILES strings derived from known inhibitors.
"""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime
import subprocess
import os


@dataclass
class GeneratedMolecule:
    """A generated drug candidate molecule."""
    smiles: str
    name: Optional[str]
    source_seed: str
    generation_method: str
    target_gene: str
    properties: Dict[str, Any]
    score: float  # Overall druglikeness/relevance score
    generated_at: str = None

    def __post_init__(self):
        if self.generated_at is None:
            self.generated_at = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class MoleculeGenerator:
    """
    Generate drug-like molecules using BioNeMo or fallback methods.

    Supports:
    - BioNeMo MolMIM (when available)
    - RDKit-based analogue generation (fallback)
    - Pre-computed demo molecules
    """

    def __init__(self, use_bionemo: bool = True, output_dir: Path = None,
                 require_real: bool = False):
        """`require_real=True` refuses the non-model fallback instead of returning a
        hand-written library that would read as generated output."""
        self.use_bionemo = use_bionemo
        self.require_real = require_real
        self.output_dir = output_dir or Path(__file__).parent.parent / "outputs" / "molecules"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check if BioNeMo is available
        self.bionemo_available = self._check_bionemo_available()
        if use_bionemo and not self.bionemo_available:
            print("Warning: BioNeMo not available, using fallback generation")

    def _check_bionemo_available(self) -> bool:
        """Check if BioNeMo container/API is accessible."""
        # Check for BioNeMo environment
        try:
            # Try importing BioNeMo modules (if running inside container)
            import importlib
            spec = importlib.util.find_spec("bionemo")
            if spec is not None:
                return True
        except ImportError:
            pass

        # Check for BioNeMo container
        try:
            result = subprocess.run(
                ["docker", "images", "nvcr.io/nvidia/bionemo/bionemo-framework", "--format", "{{.Repository}}"],
                capture_output=True, text=True, timeout=5
            )
            if "bionemo" in result.stdout.lower():
                return True
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

        return False

    def generate_from_seed(
        self,
        seed_smiles: str,
        target_gene: str,
        num_molecules: int = 10,
        diversity: float = 0.3,
    ) -> List[GeneratedMolecule]:
        """
        Generate molecules from a seed SMILES string.

        Args:
            seed_smiles: Starting molecule SMILES
            target_gene: Target gene for context
            num_molecules: Number of molecules to generate
            diversity: How different from seed (0-1)

        Returns:
            List of generated molecules
        """
        if self.use_bionemo and self.bionemo_available:
            return self._generate_bionemo(seed_smiles, target_gene, num_molecules, diversity)

        if self.require_real:
            raise RuntimeError(
                "Real generative-model output was required, but no MolMIM/BioNeMo backend is "
                "reachable from MoleculeGenerator. This module's only other path is a "
                "hand-written analogue library, which is not model output. Use "
                "src.nim_clients.create_nim_clients() with NIM_MODE=cloud|local instead."
            )

        return self._generate_fallback(seed_smiles, target_gene, num_molecules, diversity)

    def _generate_bionemo(
        self,
        seed_smiles: str,
        target_gene: str,
        num_molecules: int,
        diversity: float,
    ) -> List[GeneratedMolecule]:
        """BioNeMo MolMIM generation -- NOT IMPLEMENTED in this module.

        This method never called BioNeMo: it printed a "BioNeMo generation ..." line and
        returned the local fallback, so a caller that asked for MolMIM got a hand-written
        analogue library labeled as if a model had produced it. There is a real MolMIM client
        (`nim_clients.MolMIMClient` / `CloudMolMIMClient`, used by `pipeline.py`); this module
        is the RDKit-only path. Rather than fake the model, say so.
        """
        raise NotImplementedError(
            "MoleculeGenerator has no BioNeMo/MolMIM implementation. Use "
            "src.nim_clients.create_nim_clients() (NIM_MODE=cloud|local) for real MolMIM "
            "generation, or construct MoleculeGenerator(use_bionemo=False) to use the "
            "explicitly-labeled RDKit/library fallback."
        )

    def _generate_fallback(
        self,
        seed_smiles: str,
        target_gene: str,
        num_molecules: int,
        diversity: float,
    ) -> List[GeneratedMolecule]:
        """Generate molecules using RDKit or pre-computed analogues."""
        molecules = []

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, Descriptors, Lipinski
            rdkit_available = True
        except ImportError:
            rdkit_available = False

        if rdkit_available and seed_smiles:
            molecules = self._generate_rdkit_analogues(seed_smiles, target_gene, num_molecules)
        else:
            # Use pre-computed demo molecules for VCP
            molecules = self._get_demo_molecules(target_gene, num_molecules)

        return molecules

    def _generate_rdkit_analogues(
        self,
        seed_smiles: str,
        target_gene: str,
        num_molecules: int,
    ) -> List[GeneratedMolecule]:
        """Generate simple analogues using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        molecules = []
        seed_mol = Chem.MolFromSmiles(seed_smiles)

        if seed_mol is None:
            return self._get_demo_molecules(target_gene, num_molecules)

        # Generate fingerprint-based similar molecules
        # In production, would use MMP or scaffold hopping
        seed_fp = AllChem.GetMorganFingerprintAsBitVect(seed_mol, 2, nBits=2048)

        # For demo, create simple modifications
        modifications = [
            ("F", "Cl"),  # Halogen swap
            ("C", "N"),   # C to N swap in rings
            ("O", "S"),   # O to S swap
        ]

        # Generate the seed as first molecule
        seed_props = self._calculate_properties(seed_mol)
        seed_props["provenance"] = "input_seed"
        molecules.append(GeneratedMolecule(
            smiles=seed_smiles,
            name=f"{target_gene}_seed",
            source_seed=seed_smiles,
            generation_method="seed",
            target_gene=target_gene,
            properties=seed_props,
            score=1.0,
        ))

        # Add pre-designed VCP inhibitor analogues. These are hand-written, NOT model output;
        # `score` is a fixed rank-order placeholder, not a prediction -- labeled as such so a
        # downstream surface cannot present it as a generative-model confidence.
        #
        # The library is VCP-specific and used to be appended for ANY target_gene, each molecule
        # stamped with the requested gene: ask the :8505 UI for EGFR and you were handed
        # VCP-Inh-A1..H1 labelled target_gene="EGFR". There is no analogue library for other
        # targets, so there is nothing honest to return for them beyond the caller's own seed.
        if target_gene.upper() != "VCP":
            seed_props["library_available"] = False
            seed_props["note"] = (
                f"no analogue library exists for {target_gene}; only the supplied seed is "
                "returned. The hand-written library in this module is VCP-specific."
            )
            return molecules

        vcp_analogues = self._get_vcp_inhibitor_analogues()
        for i, (smiles, name) in enumerate(vcp_analogues[:num_molecules-1]):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                props = self._calculate_properties(mol)
                props["provenance"] = "hand_written_analogue_library"
                props["score_is_placeholder"] = True
                molecules.append(GeneratedMolecule(
                    smiles=smiles,
                    name=name,
                    source_seed=seed_smiles,
                    generation_method="analogue_library",
                    target_gene=target_gene,
                    properties=props,
                    score=0.85 - (i * 0.05),
                ))

        return molecules[:num_molecules]

    def _calculate_properties(self, mol) -> Dict[str, Any]:
        """Calculate molecular properties using RDKit."""
        from rdkit.Chem import Descriptors, Lipinski

        return {
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "logP": round(Descriptors.MolLogP(mol), 2),
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "lipinski_violations": self._count_lipinski_violations(mol),
        }

    def _count_lipinski_violations(self, mol) -> int:
        """Count Lipinski Rule of 5 violations."""
        from rdkit.Chem import Descriptors, Lipinski

        violations = 0
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Lipinski.NumHDonors(mol) > 5:
            violations += 1
        if Lipinski.NumHAcceptors(mol) > 10:
            violations += 1
        return violations

    def _get_vcp_inhibitor_analogues(self) -> List[Tuple[str, str]]:
        """Pre-designed VCP inhibitor analogues for demo."""
        return [
            # CB-5083 analogues and related VCP inhibitors
            ("CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-A1"),
            ("CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCNCC5)cc4)c1", "VCP-Inh-A2"),
            ("Cc1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-B1"),
            ("CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)CC)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-C1"),
            ("CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4cccc(CN5CCOCC5)c4)c1", "VCP-Inh-D1"),
            # Additional scaffolds
            ("COc1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-E1"),
            ("CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(N5CCOCC5)cc4)c1", "VCP-Inh-F1"),
            ("Fc1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-G1"),
            ("CC(C)c1cc(F)c(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1", "VCP-Inh-H1"),
        ]

    def _get_demo_molecules(self, target_gene: str, num_molecules: int) -> List[GeneratedMolecule]:
        """Get the hand-written demo analogue library (no model involved).

        The properties used to be invented from the list index (`450 + i * 10` Da,
        `3.5 + i * 0.1` logP, hbd/hba/violations hard-coded) -- plausible-looking numbers that
        belonged to no molecule. They are now computed from the SMILES with RDKit; without
        RDKit they are reported as unavailable rather than guessed.
        """
        if target_gene.upper() != "VCP":
            # No library for other targets; say nothing rather than invent molecules.
            return []

        try:
            from rdkit import Chem
            has_rdkit = True
        except ImportError:
            has_rdkit = False

        analogues = self._get_vcp_inhibitor_analogues()
        molecules = []
        for i, (smiles, name) in enumerate(analogues[:num_molecules]):
            props: Dict[str, Any] = {}
            if has_rdkit:
                rdmol = Chem.MolFromSmiles(smiles)
                if rdmol is not None:
                    props = self._calculate_properties(rdmol)
            if not props:
                props = {"properties_unavailable": "RDKit not installed; not computed"}
            props["provenance"] = "hand_written_analogue_library"
            props["score_is_placeholder"] = True
            molecules.append(GeneratedMolecule(
                smiles=smiles,
                name=name,
                source_seed="CB-5083",
                generation_method="demo_library",
                target_gene=target_gene,
                properties=props,
                score=0.9 - (i * 0.05),
            ))
        return molecules

    def save_molecules(self, molecules: List[GeneratedMolecule], filename: str = None) -> Path:
        """Save generated molecules to JSON."""
        if filename is None:
            filename = f"generated_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

        output_file = self.output_dir / filename
        data = [m.to_dict() for m in molecules]

        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)

        return output_file


def generate_vcp_molecules(num_molecules: int = 10) -> List[GeneratedMolecule]:
    """Convenience function to generate VCP inhibitor candidates."""
    try:
        from .cryoem_evidence import get_vcp_inhibitor_seed
    except ImportError:
        from cryoem_evidence import get_vcp_inhibitor_seed

    # This module has no generative-model backend (see _generate_bionemo); ask for the
    # labeled RDKit/library path explicitly rather than appearing to try MolMIM first.
    generator = MoleculeGenerator(use_bionemo=False)
    seed = get_vcp_inhibitor_seed()

    if seed is None:
        # Use CB-5083 as default seed
        seed = "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5"

    return generator.generate_from_seed(
        seed_smiles=seed,
        target_gene="VCP",
        num_molecules=num_molecules,
        diversity=0.3,
    )
