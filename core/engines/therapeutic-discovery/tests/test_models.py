"""
Tests for Drug Discovery Pipeline data models.
"""
import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import (
    TargetHypothesis,
    TargetStatus,
    Confidence,
    StructureInfo,
    StructureManifest,
    GeneratedMolecule,
    MoleculeProperties,
    DockingResult,
    RankedCandidate,
    PipelineConfig,
)


class TestTargetHypothesis:
    """Tests for TargetHypothesis model."""

    def test_create_basic(self):
        """Test basic target creation."""
        target = TargetHypothesis(gene="VCP")
        assert target.gene == "VCP"
        assert target.confidence == Confidence.MEDIUM
        assert target.priority == 3

    def test_create_full(self, sample_target):
        """Test full target creation."""
        assert sample_target.gene == "VCP"
        assert sample_target.protein == "Valosin-containing protein (p97)"
        assert sample_target.uniprot_id == "P55072"
        assert sample_target.priority == 5
        assert len(sample_target.pdb_ids) == 2

    def test_priority_validation(self):
        """Test priority must be 1-5."""
        with pytest.raises(ValueError):
            TargetHypothesis(gene="TEST", priority=0)

        with pytest.raises(ValueError):
            TargetHypothesis(gene="TEST", priority=6)

    def test_dict_conversion(self, sample_target):
        """Test dictionary conversion."""
        d = sample_target.dict()
        assert d["gene"] == "VCP"
        assert d["confidence"] == "high"
        assert d["status"] == "validated"


class TestStructureInfo:
    """Tests for StructureInfo model."""

    def test_create_basic(self):
        """Test basic structure creation."""
        struct = StructureInfo(pdb_id="5FTK")
        assert struct.pdb_id == "5FTK"
        assert struct.method == "X-ray"
        assert struct.chain == "A"

    def test_resolution_parsing(self):
        """Test resolution parsing from string."""
        struct = StructureInfo(pdb_id="5FTK", resolution="2.3 Å")
        assert struct.resolution == 2.3

    def test_with_ligand(self, sample_structure):
        """Test structure with ligand info."""
        assert sample_structure.ligand_id == "CB5"
        assert sample_structure.ligand_smiles is not None
        assert sample_structure.prepared == True


class TestStructureManifest:
    """Tests for StructureManifest model."""

    def test_create_manifest(self, sample_structure):
        """Test manifest creation."""
        manifest = StructureManifest(
            target_gene="VCP",
            structures=[sample_structure],
        )
        assert manifest.target_gene == "VCP"
        assert len(manifest.structures) == 1

    def test_get_best_structure(self):
        """Test getting best resolution structure."""
        structures = [
            StructureInfo(pdb_id="A", resolution=3.0),
            StructureInfo(pdb_id="B", resolution=2.0),
            StructureInfo(pdb_id="C", resolution=2.5),
        ]
        manifest = StructureManifest(
            target_gene="TEST",
            structures=structures,
        )

        best = manifest.get_best_structure()
        assert best.pdb_id == "B"
        assert best.resolution == 2.0


class TestGeneratedMolecule:
    """Tests for GeneratedMolecule model."""

    def test_create_basic(self):
        """Test basic molecule creation."""
        mol = GeneratedMolecule(
            id="mol-001",
            smiles="CCO",
            target_gene="VCP",
        )
        assert mol.id == "mol-001"
        assert mol.smiles == "CCO"
        assert mol.generation_method == "MolMIM"

    def test_with_properties(self, sample_molecule):
        """Test molecule with properties."""
        assert sample_molecule.properties is not None
        assert sample_molecule.properties.molecular_weight == 489.6
        assert sample_molecule.properties.lipinski_violations == 0
        assert sample_molecule.passed_qc == True


class TestMoleculeProperties:
    """Tests for MoleculeProperties model."""

    def test_create_properties(self):
        """Test properties creation."""
        props = MoleculeProperties(
            molecular_weight=450.0,
            logP=3.5,
            hbd=2,
            hba=6,
            tpsa=80.0,
            rotatable_bonds=5,
        )
        assert props.molecular_weight == 450.0
        assert props.lipinski_violations == 0

    def test_qed_validation(self):
        """Test QED must be 0-1."""
        props = MoleculeProperties(
            molecular_weight=450.0,
            logP=3.5,
            hbd=2,
            hba=6,
            tpsa=80.0,
            rotatable_bonds=5,
            qed=0.75,
        )
        assert props.qed == 0.75

        with pytest.raises(ValueError):
            MoleculeProperties(
                molecular_weight=450.0,
                logP=3.5,
                hbd=2,
                hba=6,
                tpsa=80.0,
                rotatable_bonds=5,
                qed=1.5,  # Invalid
            )


class TestDockingResult:
    """Tests for DockingResult model."""

    def test_create_result(self):
        """Test docking result creation."""
        result = DockingResult(
            molecule_id="mol-001",
            structure_id="5FTK",
            docking_score=-8.5,
        )
        assert result.molecule_id == "mol-001"
        assert result.docking_score == -8.5
        assert result.method == "DiffDock"

    def test_with_details(self):
        """Test result with full details."""
        result = DockingResult(
            molecule_id="mol-001",
            structure_id="5FTK",
            docking_score=-9.2,
            binding_energy=-10.5,
            confidence=0.85,
            hydrogen_bonds=3,
            contacts=["ALA123", "GLY456"],
        )
        assert result.confidence == 0.85
        assert len(result.contacts) == 2


class TestRankedCandidate:
    """Tests for RankedCandidate model."""

    def test_create_candidate(self):
        """Test candidate creation."""
        candidate = RankedCandidate(
            rank=1,
            molecule_id="mol-001",
            smiles="CCO",
            target_gene="VCP",
            docking_score=-8.5,
            generation_score=0.9,
            composite_score=0.85,
        )
        assert candidate.rank == 1
        assert candidate.composite_score == 0.85

    def test_rank_validation(self):
        """Test rank must be >= 1."""
        with pytest.raises(ValueError):
            RankedCandidate(
                rank=0,  # Invalid
                molecule_id="mol-001",
                smiles="CCO",
                target_gene="VCP",
                docking_score=-8.5,
                generation_score=0.9,
                composite_score=0.85,
            )


class TestPipelineConfig:
    """Tests for PipelineConfig model."""

    def test_create_basic(self):
        """Test basic config creation."""
        config = PipelineConfig(target_gene="VCP")
        assert config.target_gene == "VCP"
        assert config.num_molecules == 50
        assert config.generation_method == "molmim"

    def test_create_custom(self, sample_config):
        """Test custom config."""
        assert sample_config.num_molecules == 10
        assert sample_config.reference_compound_smiles is not None

    def test_weight_defaults(self):
        """Test scoring weights default correctly."""
        config = PipelineConfig(target_gene="VCP")
        assert config.docking_weight == 0.4
        assert config.generation_weight == 0.3
        assert config.qed_weight == 0.3
        # Weights should sum to 1
        total = config.docking_weight + config.generation_weight + config.qed_weight
        assert abs(total - 1.0) < 0.001
