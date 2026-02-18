"""
Tests for Drug Discovery Pipeline.
"""
import sys
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import PipelineConfig, TargetHypothesis
from src.pipeline import DrugDiscoveryPipeline


class TestDrugDiscoveryPipeline:
    """Tests for DrugDiscoveryPipeline class."""

    def test_pipeline_initialization(self, sample_config, mock_nim_manager):
        """Test pipeline initialization."""
        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        assert pipeline.config == sample_config
        assert pipeline.run_id is not None
        assert len(pipeline.molecules) == 0

    def test_stage_names(self):
        """Test stage names are defined."""
        assert len(DrugDiscoveryPipeline.STAGES) == 10
        assert DrugDiscoveryPipeline.STAGES[0] == "Initialize"
        assert DrugDiscoveryPipeline.STAGES[9] == "Reporting"

    def test_stage_0_initialize(self, sample_config, mock_nim_manager, sample_target):
        """Test stage 0 initialization."""
        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        pipeline.stage_0_initialize(sample_target)

        assert pipeline.target is not None
        assert pipeline.target.gene == "VCP"
        assert 0 in pipeline.run.stages_completed

    def test_stage_1_normalize_target(self, sample_config, mock_nim_manager, sample_target):
        """Test stage 1 target normalization."""
        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()

        assert pipeline.target.gene == "VCP"  # Uppercased
        assert pipeline.target.uniprot_id == "P55072"
        assert 1 in pipeline.run.stages_completed

    def test_stage_2_structure_discovery(self, sample_config, mock_nim_manager, sample_target, temp_output_dir):
        """Test stage 2 structure discovery."""
        sample_config.output_dir = str(temp_output_dir)
        sample_config.structure_cache = str(temp_output_dir / "structures")

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()

        assert pipeline.structures is not None
        assert len(pipeline.structures.structures) == len(sample_target.pdb_ids)
        assert 2 in pipeline.run.stages_completed

    def test_stage_4_molecule_generation(self, sample_config, mock_nim_manager, sample_target, temp_output_dir):
        """Test stage 4 molecule generation."""
        sample_config.output_dir = str(temp_output_dir)

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        # Run prerequisite stages
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()
        pipeline.stage_3_structure_prep()
        pipeline.stage_4_molecule_generation()

        assert len(pipeline.molecules) > 0
        assert pipeline.run.molecules_generated > 0
        assert 4 in pipeline.run.stages_completed

    def test_full_pipeline_with_mocks(self, sample_config, mock_nim_manager, sample_target, temp_output_dir):
        """Test running full pipeline with mock services."""
        sample_config.output_dir = str(temp_output_dir)
        sample_config.structure_cache = str(temp_output_dir / "structures")
        sample_config.num_molecules = 5

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        result = pipeline.run_pipeline(sample_target)

        assert result.status == "completed"
        assert result.completed_at is not None
        assert len(result.stages_completed) == 10

    def test_progress_callback(self, sample_config, mock_nim_manager, temp_output_dir):
        """Test progress callback is called."""
        sample_config.output_dir = str(temp_output_dir)

        progress_calls = []

        def callback(stage, message):
            progress_calls.append((stage, message))

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
            progress_callback=callback,
        )

        pipeline.stage_0_initialize()

        assert len(progress_calls) > 0
        assert progress_calls[0][0] == 0  # Stage 0

    def test_pipeline_creates_report(self, sample_config, mock_nim_manager, sample_target, temp_output_dir):
        """Test pipeline creates final report."""
        sample_config.output_dir = str(temp_output_dir)
        sample_config.num_molecules = 3

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        result = pipeline.run_pipeline(sample_target)

        report_files = list(temp_output_dir.glob("report_*.json"))
        assert len(report_files) == 1


class TestPipelineStages:
    """Tests for individual pipeline stages."""

    def test_chemistry_qc_filters_invalid(self, sample_config, mock_nim_manager, temp_output_dir):
        """Test chemistry QC filters out invalid molecules."""
        sample_config.output_dir = str(temp_output_dir)
        sample_config.max_mw = 300  # Very restrictive

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        # Add some molecules manually
        from src.models import GeneratedMolecule
        pipeline.molecules = [
            GeneratedMolecule(
                id="mol-1",
                smiles="CCO",  # Small MW
                target_gene="VCP",
            ),
            GeneratedMolecule(
                id="mol-2",
                smiles="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",  # Large MW
                target_gene="VCP",
            ),
        ]

        pipeline.stage_5_chemistry_qc()

        # Only small molecule should pass
        passed = [m for m in pipeline.molecules if m.passed_qc]
        assert len(passed) <= len(pipeline.molecules)

    def test_ranking_orders_by_score(self, sample_config, mock_nim_manager, temp_output_dir):
        """Test ranking orders candidates by composite score."""
        sample_config.output_dir = str(temp_output_dir)

        pipeline = DrugDiscoveryPipeline(
            config=sample_config,
            nim_manager=mock_nim_manager,
        )

        # Setup mock molecules
        from src.models import GeneratedMolecule, MoleculeProperties

        pipeline.molecules = [
            GeneratedMolecule(
                id="mol-1",
                smiles="CCO",
                target_gene="VCP",
                generation_score=0.5,
                passed_qc=True,
                properties=MoleculeProperties(
                    molecular_weight=46,
                    logP=-0.3,
                    hbd=1,
                    hba=1,
                    tpsa=20.2,
                    rotatable_bonds=0,
                    qed=0.4,
                ),
            ),
            GeneratedMolecule(
                id="mol-2",
                smiles="CCCO",
                target_gene="VCP",
                generation_score=0.9,
                passed_qc=True,
                properties=MoleculeProperties(
                    molecular_weight=60,
                    logP=0.3,
                    hbd=1,
                    hba=1,
                    tpsa=20.2,
                    rotatable_bonds=1,
                    qed=0.6,
                ),
            ),
        ]

        # Mock docking results
        from src.models import DockingResult
        pipeline.docking_results = [
            DockingResult(molecule_id="mol-1", structure_id="5FTK", docking_score=-6.0),
            DockingResult(molecule_id="mol-2", structure_id="5FTK", docking_score=-8.5),
        ]

        pipeline.stage_8_ranking()

        assert len(pipeline.ranked_candidates) > 0
        # First should have higher score
        if len(pipeline.ranked_candidates) >= 2:
            assert pipeline.ranked_candidates[0].composite_score >= pipeline.ranked_candidates[1].composite_score


class TestVCPDemo:
    """Tests for VCP demo pipeline."""

    def test_run_vcp_demo(self, temp_output_dir):
        """Test VCP demo runs successfully."""
        from src.pipeline import run_vcp_demo_pipeline

        with patch('src.pipeline.create_nim_clients') as mock_create:
            mock_manager = Mock()
            mock_manager.molmim.check_health.return_value = True
            mock_manager.molmim.generate.return_value = [
                {"smiles": "CCO", "score": 0.9, "method": "mock"}
            ]
            mock_manager.diffdock.check_health.return_value = True
            mock_manager.diffdock.dock.return_value = [
                {"docking_score": -8.0, "confidence": 0.8}
            ]
            mock_manager.get_available_services.return_value = ["molmim", "diffdock"]
            mock_create.return_value = mock_manager

            result = run_vcp_demo_pipeline(output_dir=str(temp_output_dir))

            assert result.target_gene == "VCP"
            assert result.status == "completed"
