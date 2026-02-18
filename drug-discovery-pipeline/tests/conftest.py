"""
Shared pytest fixtures for Drug Discovery Pipeline tests.
"""
import sys
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import (
    GeneratedMolecule,
    MoleculeProperties,
    PipelineConfig,
    StructureInfo,
    TargetHypothesis,
)


@pytest.fixture
def sample_target():
    """Create a sample VCP target."""
    return TargetHypothesis(
        gene="VCP",
        protein="Valosin-containing protein (p97)",
        uniprot_id="P55072",
        rationale="VCP mutations cause FTD",
        therapeutic_area="Neurodegeneration",
        mechanism="AAA+ ATPase inhibition",
        confidence="high",
        priority=5,
        pdb_ids=["5FTK", "8OOI"],
        status="validated",
    )


@pytest.fixture
def sample_structure():
    """Create a sample structure."""
    return StructureInfo(
        pdb_id="5FTK",
        method="X-ray",
        resolution=2.3,
        chain="A",
        ligand_id="CB5",
        ligand_smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
        prepared=True,
    )


@pytest.fixture
def sample_molecule():
    """Create a sample molecule."""
    return GeneratedMolecule(
        id="test-mol-0001",
        smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
        name="VCP-Test-01",
        source_seed="CB-5083",
        generation_method="MolMIM",
        target_gene="VCP",
        generation_score=0.85,
        passed_qc=True,
        properties=MoleculeProperties(
            molecular_weight=489.6,
            logP=4.2,
            hbd=2,
            hba=6,
            tpsa=78.5,
            rotatable_bonds=7,
            lipinski_violations=0,
            qed=0.65,
        ),
    )


@pytest.fixture
def sample_config():
    """Create a sample pipeline config."""
    return PipelineConfig(
        target_gene="VCP",
        reference_compound_smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
        num_molecules=10,
        output_dir="test_outputs",
    )


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_nim_manager():
    """Create a mock NIM service manager."""
    manager = Mock()

    # Mock MolMIM
    manager.molmim = Mock()
    manager.molmim.check_health.return_value = True
    manager.molmim.generate.return_value = [
        {"smiles": "CCO", "score": 0.9, "method": "mock"},
        {"smiles": "CCCO", "score": 0.85, "method": "mock"},
    ]

    # Mock DiffDock
    manager.diffdock = Mock()
    manager.diffdock.check_health.return_value = True
    manager.diffdock.dock.return_value = [
        {"docking_score": -8.5, "confidence": 0.8, "pose_id": 1},
    ]

    manager.get_available_services.return_value = ["molmim", "diffdock"]

    return manager
