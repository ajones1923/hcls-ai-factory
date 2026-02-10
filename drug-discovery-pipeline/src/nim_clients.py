"""
NIM Microservice Clients for Drug Discovery Pipeline.

Provides clients for:
- MolMIM: Molecule generation via masked modeling
- DiffDock: Molecular docking via diffusion models

Based on NVIDIA NIM API specifications from phase-5-6.pdf.
"""
import os
import json
import time
import requests
from pathlib import Path
from typing import List, Dict, Any, Optional, Generator
from dataclasses import dataclass
from loguru import logger

from .models import GeneratedMolecule, DockingResult, MoleculeProperties


@dataclass
class NIMServiceConfig:
    """Configuration for NIM service connection."""
    host: str = "localhost"
    port: int = 8000
    api_version: str = "v1"
    timeout: int = 300
    max_retries: int = 3

    @property
    def base_url(self) -> str:
        return f"http://{self.host}:{self.port}/{self.api_version}"


class MolMIMClient:
    """
    Client for MolMIM molecule generation service.

    MolMIM uses masked language modeling to generate novel molecules
    by sampling from a learned molecular distribution.
    """

    def __init__(self, config: NIMServiceConfig = None):
        self.config = config or NIMServiceConfig(port=8001)
        self.health_checked = False

    def check_health(self) -> bool:
        """Check if the MolMIM service is available."""
        try:
            response = requests.get(
                f"{self.config.base_url}/health",
                timeout=5
            )
            self.health_checked = response.status_code == 200
            return self.health_checked
        except requests.RequestException as e:
            logger.warning(f"MolMIM health check failed: {e}")
            return False

    def generate(
        self,
        seed_smiles: str,
        num_molecules: int = 10,
        temperature: float = 1.0,
        num_samples_per_token: int = 10,
        masked_ratio: float = 0.1,
    ) -> List[Dict[str, Any]]:
        """
        Generate molecules using MolMIM.

        Args:
            seed_smiles: Starting molecule SMILES
            num_molecules: Number of molecules to generate
            temperature: Sampling temperature
            num_samples_per_token: Samples per masked token
            masked_ratio: Ratio of tokens to mask

        Returns:
            List of generated molecule dictionaries

        Raises:
            RuntimeError: If generation fails after all retries
        """
        if not self.health_checked:
            self.check_health()

        payload = {
            "smiles": seed_smiles,
            "num_molecules": num_molecules,
            "temperature": temperature,
            "num_samples_per_token": num_samples_per_token,
            "masked_ratio": masked_ratio,
        }

        last_error = None
        for attempt in range(1, self.config.max_retries + 1):
            try:
                response = requests.post(
                    f"{self.config.base_url}/generate",
                    json=payload,
                    timeout=self.config.timeout,
                )
                response.raise_for_status()
                molecules = response.json().get("molecules", [])
                if not molecules:
                    logger.warning(f"MolMIM returned 0 molecules (attempt {attempt}/{self.config.max_retries})")
                return molecules

            except requests.RequestException as e:
                last_error = e
                if attempt < self.config.max_retries:
                    wait = 2 ** attempt  # Exponential backoff: 2s, 4s, 8s
                    logger.warning(f"MolMIM generation failed (attempt {attempt}/{self.config.max_retries}), retrying in {wait}s: {e}")
                    time.sleep(wait)
                else:
                    logger.error(f"MolMIM generation failed after {self.config.max_retries} attempts: {e}")

        raise RuntimeError(
            f"MolMIM generation failed after {self.config.max_retries} attempts. "
            f"Last error: {last_error}. "
            f"Check that MolMIM is running at {self.config.base_url}"
        )

    def generate_batch(
        self,
        seed_smiles_list: List[str],
        num_molecules_per_seed: int = 5,
        **kwargs
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Generate molecules from multiple seeds.

        Returns:
            Dictionary mapping seed SMILES to generated molecules
        """
        results = {}
        for seed in seed_smiles_list:
            results[seed] = self.generate(
                seed_smiles=seed,
                num_molecules=num_molecules_per_seed,
                **kwargs
            )
        return results


class DiffDockClient:
    """
    Client for DiffDock molecular docking service.

    DiffDock uses diffusion models to predict protein-ligand
    binding poses directly from sequence and structure.
    """

    def __init__(self, config: NIMServiceConfig = None):
        self.config = config or NIMServiceConfig(port=8002)
        self.health_checked = False

    def check_health(self) -> bool:
        """Check if the DiffDock service is available."""
        try:
            response = requests.get(
                f"{self.config.base_url}/health",
                timeout=5
            )
            self.health_checked = response.status_code == 200
            return self.health_checked
        except requests.RequestException as e:
            logger.warning(f"DiffDock health check failed: {e}")
            return False

    def dock(
        self,
        protein_pdb: str,
        ligand_smiles: str,
        num_poses: int = 10,
        confidence_threshold: float = 0.0,
    ) -> List[Dict[str, Any]]:
        """
        Dock a ligand to a protein structure.

        Args:
            protein_pdb: Path to protein PDB file or PDB string
            ligand_smiles: SMILES string of the ligand
            num_poses: Number of poses to generate
            confidence_threshold: Minimum confidence score

        Returns:
            List of docking pose dictionaries
        """
        if not self.health_checked:
            self.check_health()

        # Read PDB if path provided
        if Path(protein_pdb).exists():
            with open(protein_pdb, 'r') as f:
                protein_pdb = f.read()

        payload = {
            "protein": protein_pdb,
            "ligand": ligand_smiles,
            "num_poses": num_poses,
        }

        last_error = None
        for attempt in range(1, self.config.max_retries + 1):
            try:
                response = requests.post(
                    f"{self.config.base_url}/dock",
                    json=payload,
                    timeout=self.config.timeout,
                )
                response.raise_for_status()

                poses = response.json().get("poses", [])

                # Filter by confidence
                if confidence_threshold > 0:
                    poses = [p for p in poses if p.get("confidence", 0) >= confidence_threshold]

                return poses

            except requests.RequestException as e:
                last_error = e
                if attempt < self.config.max_retries:
                    wait = 2 ** attempt
                    logger.warning(f"DiffDock docking failed (attempt {attempt}/{self.config.max_retries}), retrying in {wait}s: {e}")
                    time.sleep(wait)
                else:
                    logger.error(f"DiffDock docking failed after {self.config.max_retries} attempts: {e}")

        raise RuntimeError(
            f"DiffDock docking failed after {self.config.max_retries} attempts. "
            f"Last error: {last_error}. "
            f"Check that DiffDock is running at {self.config.base_url}"
        )

    def dock_batch(
        self,
        protein_pdb: str,
        ligand_smiles_list: List[str],
        **kwargs
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Dock multiple ligands to a protein.

        Returns:
            Dictionary mapping SMILES to docking poses
        """
        results = {}
        for smiles in ligand_smiles_list:
            results[smiles] = self.dock(
                protein_pdb=protein_pdb,
                ligand_smiles=smiles,
                **kwargs
            )
        return results


class NIMServiceManager:
    """
    Manages NIM service connections and provides fallback behavior.
    """

    def __init__(
        self,
        molmim_config: NIMServiceConfig = None,
        diffdock_config: NIMServiceConfig = None,
    ):
        self.molmim_config = molmim_config or NIMServiceConfig(port=8001)
        self.diffdock_config = diffdock_config or NIMServiceConfig(port=8002)

        self._molmim: Optional[MolMIMClient] = None
        self._diffdock: Optional[DiffDockClient] = None

    @property
    def molmim(self) -> MolMIMClient:
        """Get or create MolMIM client."""
        if self._molmim is None:
            self._molmim = MolMIMClient(self.molmim_config)
        return self._molmim

    @property
    def diffdock(self) -> DiffDockClient:
        """Get or create DiffDock client."""
        if self._diffdock is None:
            self._diffdock = DiffDockClient(self.diffdock_config)
        return self._diffdock

    def check_services(self) -> Dict[str, bool]:
        """Check health of all NIM services."""
        return {
            "molmim": self.molmim.check_health(),
            "diffdock": self.diffdock.check_health(),
        }

    def get_available_services(self) -> List[str]:
        """Get list of available services."""
        status = self.check_services()
        return [name for name, available in status.items() if available]


class MockMolMIMClient(MolMIMClient):
    """
    Mock MolMIM client for testing and demo purposes.
    Uses RDKit-based generation as fallback.
    """

    def check_health(self) -> bool:
        """Mock always returns True."""
        self.health_checked = True
        return True

    def generate(
        self,
        seed_smiles: str,
        num_molecules: int = 10,
        **kwargs
    ) -> List[Dict[str, Any]]:
        """Generate molecules using RDKit fallback."""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, Descriptors
            HAS_RDKIT = True
        except ImportError:
            HAS_RDKIT = False
            logger.warning("RDKit not available for mock generation")
            return []

        molecules = []
        seed_mol = Chem.MolFromSmiles(seed_smiles)

        if seed_mol is None:
            logger.warning(f"Invalid seed SMILES: {seed_smiles}")
            return []

        # Generate the seed as first molecule
        molecules.append({
            "smiles": seed_smiles,
            "score": 1.0,
            "method": "seed",
        })

        # Generate simple analogues by modifying atoms
        # In production, MolMIM would do this properly
        atom_swaps = [
            ("F", "Cl"),
            ("C", "N"),
            ("O", "S"),
        ]

        for old, new in atom_swaps:
            modified = seed_smiles.replace(old, new, 1)
            mol = Chem.MolFromSmiles(modified)
            if mol and modified != seed_smiles:
                molecules.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "score": 0.85,
                    "method": "atom_swap",
                })
                if len(molecules) >= num_molecules:
                    break

        # Add some pre-designed VCP inhibitor analogues if targeting VCP
        vcp_analogues = [
            "CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1",
            "Cc1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1",
            "CC(C)c1ccc(NC2=NC3=C(C=N2)N(C=C3)CC)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1",
        ]

        for i, smiles in enumerate(vcp_analogues):
            if len(molecules) >= num_molecules:
                break
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                molecules.append({
                    "smiles": smiles,
                    "score": 0.8 - i * 0.05,
                    "method": "library",
                })

        return molecules[:num_molecules]


class MockDiffDockClient(DiffDockClient):
    """
    Mock DiffDock client for testing and demo purposes.
    Returns simulated docking scores.
    """

    def check_health(self) -> bool:
        """Mock always returns True."""
        self.health_checked = True
        return True

    def dock(
        self,
        protein_pdb: str,
        ligand_smiles: str,
        num_poses: int = 10,
        **kwargs
    ) -> List[Dict[str, Any]]:
        """Generate simulated docking results."""
        import random
        import hashlib

        # Use hash of inputs for reproducible "random" scores
        seed = int(hashlib.md5(f"{protein_pdb[:100]}{ligand_smiles}".encode()).hexdigest()[:8], 16)
        random.seed(seed)

        poses = []
        for i in range(num_poses):
            # Simulate realistic docking scores
            # Lower is better for docking scores
            base_score = -8.0 + random.gauss(0, 2)
            confidence = max(0, min(1, 0.7 + random.gauss(0, 0.2)))

            poses.append({
                "pose_id": i + 1,
                "docking_score": round(base_score + i * 0.5, 2),
                "confidence": round(confidence - i * 0.05, 3),
                "rmsd": round(random.uniform(0.5, 3.0), 2) if i > 0 else 0.0,
                "hydrogen_bonds": random.randint(0, 4),
                "contacts": [
                    f"ALA{random.randint(100, 500)}",
                    f"GLY{random.randint(100, 500)}",
                    f"ASP{random.randint(100, 500)}",
                ],
            })

        # Sort by score
        poses.sort(key=lambda x: x["docking_score"])
        return poses


def create_nim_clients(use_mock: bool = False) -> NIMServiceManager:
    """
    Factory function to create NIM clients.

    Args:
        use_mock: If True, use mock clients for testing

    Returns:
        NIMServiceManager instance
    """
    if use_mock:
        manager = NIMServiceManager()
        manager._molmim = MockMolMIMClient()
        manager._diffdock = MockDiffDockClient()
        return manager

    # Check environment for service URLs
    molmim_url = os.environ.get("MOLMIM_URL", "http://localhost:8001")
    diffdock_url = os.environ.get("DIFFDOCK_URL", "http://localhost:8002")

    # Parse URLs into configs
    def parse_url(url: str) -> NIMServiceConfig:
        from urllib.parse import urlparse
        parsed = urlparse(url)
        return NIMServiceConfig(
            host=parsed.hostname or "localhost",
            port=parsed.port or 8000,
        )

    manager = NIMServiceManager(
        molmim_config=parse_url(molmim_url),
        diffdock_config=parse_url(diffdock_url),
    )

    # Check if real services are available
    available = manager.get_available_services()
    if not available:
        # Check if we should allow mock fallback (for development/testing only)
        allow_mock = os.environ.get("NIM_ALLOW_MOCK_FALLBACK", "false").lower() == "true"

        if allow_mock:
            logger.warning("=" * 60)
            logger.warning("WARNING: No NIM services available!")
            logger.warning(f"  Expected MolMIM at: {molmim_url}")
            logger.warning(f"  Expected DiffDock at: {diffdock_url}")
            logger.warning("USING MOCK CLIENTS - Results are NOT scientifically valid!")
            logger.warning("=" * 60)
            return create_nim_clients(use_mock=True)
        else:
            logger.error("=" * 60)
            logger.error("CRITICAL: No NIM services available!")
            logger.error(f"  Expected MolMIM at: {molmim_url}")
            logger.error(f"  Expected DiffDock at: {diffdock_url}")
            logger.error("")
            logger.error("To fix this:")
            logger.error("  1. Start the NIM containers: docker-compose up -d")
            logger.error("  2. Ensure NGC_API_KEY is set")
            logger.error("")
            logger.error("To use mock data (dev only): export NIM_ALLOW_MOCK_FALLBACK=true")
            logger.error("=" * 60)
            raise RuntimeError(
                "NIM services unavailable. Set NIM_ALLOW_MOCK_FALLBACK=true to use mock data."
            )

    return manager
