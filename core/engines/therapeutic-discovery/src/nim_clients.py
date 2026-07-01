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
    max_retries: int = 5

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
        # Validate SMILES before sending to NIM service
        try:
            from rdkit import Chem
            if Chem.MolFromSmiles(seed_smiles) is None:
                raise ValueError(f"Invalid SMILES: {seed_smiles!r}")
        except ImportError:
            pass  # RDKit optional — skip validation if not installed

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
        # Validate SMILES before sending to NIM service
        try:
            from rdkit import Chem
            if Chem.MolFromSmiles(ligand_smiles) is None:
                raise ValueError(f"Invalid SMILES: {ligand_smiles!r}")
        except ImportError:
            pass  # RDKit optional — skip validation if not installed

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


class CloudMolMIMClient(MolMIMClient):
    """
    MolMIM client using NVIDIA's hosted cloud API at health.api.nvidia.com.
    """

    def __init__(self, api_key: str, endpoint_url: str = None):
        super().__init__()
        self.api_key = api_key
        self.endpoint_url = endpoint_url or "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"

    def check_health(self) -> bool:
        """Cloud service is assumed healthy if API key is set."""
        self.health_checked = bool(self.api_key)
        return self.health_checked

    def generate(
        self,
        seed_smiles: str,
        num_molecules: int = 10,
        temperature: float = 1.0,
        num_samples_per_token: int = 10,
        masked_ratio: float = 0.1,
        **kwargs,
    ) -> List[Dict[str, Any]]:
        """Generate molecules using NVIDIA cloud MolMIM API."""
        try:
            from rdkit import Chem
            if Chem.MolFromSmiles(seed_smiles) is None:
                raise ValueError(f"Invalid SMILES: {seed_smiles!r}")
        except ImportError:
            pass

        payload = {
            "smi": seed_smiles,
            "num_molecules": num_molecules,
            "algorithm": "CMA-ES",
            "property_name": "QED",
            "iterations": 10,
        }

        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

        last_error = None
        for attempt in range(1, self.config.max_retries + 1):
            try:
                response = requests.post(
                    self.endpoint_url,
                    json=payload,
                    headers=headers,
                    timeout=self.config.timeout,
                )
                response.raise_for_status()
                data = response.json()

                # Cloud API returns {"molecules": "[{...}]", "score_type": "QED"}
                molecules_raw = data.get("molecules", "[]")
                if isinstance(molecules_raw, str):
                    import json as _json
                    molecules_list = _json.loads(molecules_raw)
                else:
                    molecules_list = molecules_raw

                # Normalize to match local NIM format
                molecules = []
                for m in molecules_list:
                    molecules.append({
                        "smiles": m.get("sample", m.get("smiles", "")),
                        "score": m.get("score", 0.0),
                        "method": "cloud_molmim",
                    })

                logger.info(f"Cloud MolMIM generated {len(molecules)} molecules")
                return molecules

            except requests.RequestException as e:
                last_error = e
                if attempt < self.config.max_retries:
                    wait = 2 ** attempt
                    logger.warning(f"Cloud MolMIM failed (attempt {attempt}/{self.config.max_retries}), retrying in {wait}s: {e}")
                    time.sleep(wait)
                else:
                    logger.error(f"Cloud MolMIM failed after {self.config.max_retries} attempts: {e}")

        raise RuntimeError(
            f"Cloud MolMIM generation failed after {self.config.max_retries} attempts. "
            f"Last error: {last_error}"
        )


class CloudDiffDockClient(DiffDockClient):
    """
    DiffDock client using NVIDIA's hosted cloud API at health.api.nvidia.com.
    Uses NVCF asset staging for file uploads.
    """

    NVCF_ASSETS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

    def __init__(self, api_key: str, endpoint_url: str = None):
        super().__init__()
        self.api_key = api_key
        self.endpoint_url = endpoint_url or "https://health.api.nvidia.com/v1/biology/mit/diffdock"

    def check_health(self) -> bool:
        """Cloud service is assumed healthy if API key is set."""
        self.health_checked = bool(self.api_key)
        return self.health_checked

    def _upload_asset(self, content: str, content_type: str, description: str) -> str:
        """Upload a file to NVCF asset storage and return the asset ID."""
        # Step 1: Create asset and get upload URL
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }
        create_resp = requests.post(
            self.NVCF_ASSETS_URL,
            json={"contentType": content_type, "description": description},
            headers=headers,
            timeout=30,
        )
        create_resp.raise_for_status()
        asset_data = create_resp.json()
        asset_id = asset_data["assetId"]
        upload_url = asset_data["uploadUrl"]

        # Step 2: Upload content to pre-signed S3 URL
        upload_headers = {
            "Content-Type": content_type,
            "x-amz-meta-nvcf-asset-description": description,
        }
        upload_resp = requests.put(
            upload_url,
            data=content.encode("utf-8"),
            headers=upload_headers,
            timeout=60,
        )
        upload_resp.raise_for_status()

        logger.info(f"Uploaded NVCF asset: {asset_id} ({description})")
        return asset_id

    def _smiles_to_sdf(self, smiles: str) -> str:
        """Convert SMILES to SDF format using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles!r}")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        return Chem.MolToMolBlock(mol)

    def dock(
        self,
        protein_pdb: str,
        ligand_smiles: str,
        num_poses: int = 10,
        confidence_threshold: float = 0.0,
        **kwargs,
    ) -> List[Dict[str, Any]]:
        """Dock a ligand to a protein using NVIDIA cloud DiffDock API."""
        try:
            from rdkit import Chem
            if Chem.MolFromSmiles(ligand_smiles) is None:
                raise ValueError(f"Invalid SMILES: {ligand_smiles!r}")
        except ImportError:
            pass

        # Read PDB file if path provided (check length first to avoid OS error)
        if len(protein_pdb) < 4096 and Path(protein_pdb).exists():
            with open(protein_pdb, 'r') as f:
                protein_pdb = f.read()

        # Upload protein PDB as NVCF asset
        protein_asset_id = self._upload_asset(
            protein_pdb, "chemical/x-pdb", "protein structure"
        )

        # Convert SMILES to SDF and upload as NVCF asset
        sdf_content = self._smiles_to_sdf(ligand_smiles)
        ligand_asset_id = self._upload_asset(
            sdf_content, "chemical/x-mdl-sdfile", "ligand molecule"
        )

        payload = {
            "ligand": ligand_asset_id,
            "ligand_file_type": "sdf",
            "protein": protein_asset_id,
            "num_poses": min(num_poses, 100),
            "time_divisions": 20,
            "steps": 18,
            "save_trajectory": False,
            "is_staged": True,
        }

        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
            "NVCF-INPUT-ASSET-REFERENCES": f"{protein_asset_id},{ligand_asset_id}",
        }

        last_error = None
        for attempt in range(1, self.config.max_retries + 1):
            try:
                response = requests.post(
                    self.endpoint_url,
                    json=payload,
                    headers=headers,
                    timeout=self.config.timeout,
                )
                response.raise_for_status()

                # Parse cloud response into standard pose format
                data = response.json()
                poses = self._parse_cloud_response(data, ligand_smiles)

                if confidence_threshold > 0:
                    poses = [p for p in poses if p.get("confidence", 0) >= confidence_threshold]

                logger.info(f"Cloud DiffDock returned {len(poses)} poses")
                return poses

            except requests.RequestException as e:
                last_error = e
                if attempt < self.config.max_retries:
                    wait = 2 ** attempt
                    logger.warning(f"Cloud DiffDock failed (attempt {attempt}/{self.config.max_retries}), retrying in {wait}s: {e}")
                    time.sleep(wait)
                else:
                    logger.error(f"Cloud DiffDock failed after {self.config.max_retries} attempts: {e}")

        raise RuntimeError(
            f"Cloud DiffDock failed after {self.config.max_retries} attempts. "
            f"Last error: {last_error}"
        )

    def _parse_cloud_response(self, data: dict, ligand_smiles: str) -> List[Dict[str, Any]]:
        """Parse cloud DiffDock response into standard pose format.

        The NVIDIA cloud API returns a structure like:
        {
            "position_confidence": [-3.77, -4.46, -5.05],  # per-pose confidence scores
            "ligand_positions": ["SDF block 1", "SDF block 2", ...],
            ...
        }
        Each index corresponds to one docking pose.
        """
        poses = []

        # Extract confidence scores (may be a list for multiple poses)
        confidence_scores = data.get("position_confidence", data.get("confidence", []))
        ligand_positions = data.get("ligand_positions", data.get("poses", []))

        # Normalize to lists
        if isinstance(confidence_scores, (int, float)):
            confidence_scores = [confidence_scores]
        if isinstance(ligand_positions, str):
            ligand_positions = [ligand_positions]

        num_poses = max(len(confidence_scores), len(ligand_positions), 1)

        for i in range(num_poses):
            score = confidence_scores[i] if i < len(confidence_scores) else 0.0
            position = ligand_positions[i] if i < len(ligand_positions) else ""

            poses.append({
                "pose_id": i + 1,
                "docking_score": float(score),
                "confidence": float(score),
                "rmsd": 0.0,
                "hydrogen_bonds": 0,
                "contacts": [],
                "ligand_position": position,
            })

        poses.sort(key=lambda x: x["docking_score"])
        return poses


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

        out = molecules[:num_molecules]
        for m in out:                       # F2: stamp mock provenance on the data itself
            m["_provenance"] = "mock"
            m["_warning"] = "SIMULATED — not a real generation; offline/demo only"
        return out


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
        for p in poses:                     # F2: stamp mock provenance on the data itself
            p["_provenance"] = "mock"
            p["_warning"] = "SIMULATED — not a real docking; offline/demo only"
        return poses


# --------------------------------------------------------------------------- #
# F2: mock-provenance honesty helpers — so simulated data can never be presented
# as real. Every mock client stamps `_provenance="mock"` on its outputs; these
# helpers let any downstream surface (pipeline, report, MCP) detect and refuse it.
# --------------------------------------------------------------------------- #
def is_mock_result(item) -> bool:
    return isinstance(item, dict) and item.get("_provenance") == "mock"


def contains_mock(results) -> bool:
    """True if any item in the results is mock-derived."""
    if isinstance(results, dict):
        results = [results]
    return any(is_mock_result(r) for r in (results or []))


def assert_real(results, context: str = "") -> None:
    """Raise if any result is mock-derived — guards 'real, never mocked' surfaces."""
    if contains_mock(results):
        where = f" in {context}" if context else ""
        raise RuntimeError(f"refusing to present MOCK results as real{where}")


def create_nim_clients(use_mock: bool = False) -> NIMServiceManager:
    """
    Factory function to create NIM clients.

    Supports three modes via NIM_MODE env var:
        - "cloud": Use NVIDIA hosted APIs at health.api.nvidia.com (default if NVIDIA_API_KEY is set)
        - "local": Use local Docker NIM containers
        - "mock": Use mock clients for testing

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

    nim_mode = os.environ.get("NIM_MODE", "").lower()
    api_key = os.environ.get("NVIDIA_API_KEY", os.environ.get("NGC_API_KEY", ""))

    # Cloud mode: use NVIDIA hosted NIM APIs
    if nim_mode == "cloud" or (not nim_mode and api_key.startswith("nvapi-")):
        if not api_key:
            raise RuntimeError(
                "NIM_MODE=cloud requires NVIDIA_API_KEY to be set. "
                "Get one at https://build.nvidia.com"
            )

        molmim_url = os.environ.get(
            "MOLMIM_URL",
            "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate",
        )
        diffdock_url = os.environ.get(
            "DIFFDOCK_URL",
            "https://health.api.nvidia.com/v1/biology/mit/diffdock",
        )

        manager = NIMServiceManager()
        manager._molmim = CloudMolMIMClient(api_key=api_key, endpoint_url=molmim_url)
        manager._diffdock = CloudDiffDockClient(api_key=api_key, endpoint_url=diffdock_url)

        logger.info("=" * 60)
        logger.info("Using NVIDIA Cloud NIM APIs")
        logger.info(f"  MolMIM:   {molmim_url}")
        logger.info(f"  DiffDock: {diffdock_url}")
        logger.info("=" * 60)
        return manager

    # Local mode: use Docker NIM containers
    molmim_url = os.environ.get("MOLMIM_URL", "http://localhost:8001")
    diffdock_url = os.environ.get("DIFFDOCK_URL", "http://localhost:8002")

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
            logger.error("  1. Set NIM_MODE=cloud with NVIDIA_API_KEY for cloud APIs")
            logger.error("  2. Or start NIM containers: docker-compose up -d")
            logger.error("")
            logger.error("To use mock data (dev only): export NIM_ALLOW_MOCK_FALLBACK=true")
            logger.error("=" * 60)
            raise RuntimeError(
                "NIM services unavailable. Set NIM_MODE=cloud with NVIDIA_API_KEY, "
                "or set NIM_ALLOW_MOCK_FALLBACK=true for mock data."
            )

    return manager
