"""
Tests for pipeline checkpoint/resume functionality.
"""
import json
import pytest
from pathlib import Path
from unittest.mock import Mock, patch, call

from src.checkpoint import CheckpointManager
from src.pipeline import DrugDiscoveryPipeline
from src.models import (
    TargetHypothesis,
    StructureInfo,
    StructureManifest,
    GeneratedMolecule,
    MoleculeProperties,
    DockingResult,
    RankedCandidate,
    PipelineConfig,
    PipelineRun,
)


class TestCheckpointManager:
    """Tests for CheckpointManager save/load operations."""

    def test_save_creates_file(self, temp_output_dir, sample_config, mock_nim_manager, sample_target):
        """Verify checkpoint file is created after save."""
        config = PipelineConfig(**{**sample_config.model_dump(), "output_dir": str(temp_output_dir)})
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)

        # Checkpoint should have been saved by the stage
        checkpoint_dir = temp_output_dir / "checkpoints"
        files = list(checkpoint_dir.glob("checkpoint_*_stage0.json"))
        assert len(files) == 1

        data = json.loads(files[0].read_text())
        assert data["stage"] == 0
        assert data["run_id"] == pipeline.run_id
        assert data["target"]["gene"] == "VCP"

    def test_load_restores_state(self, temp_output_dir, sample_config, mock_nim_manager, sample_target):
        """Save state, load into new pipeline, verify fields match."""
        config = PipelineConfig(**{**sample_config.model_dump(), "output_dir": str(temp_output_dir)})
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()

        run_id = pipeline.run_id

        # Create a new pipeline and restore
        pipeline2 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        mgr = CheckpointManager(temp_output_dir)
        checkpoint = mgr.load(run_id, 1)
        assert checkpoint is not None

        mgr.restore_pipeline_state(pipeline2, checkpoint)
        assert pipeline2.target.gene == "VCP"
        assert pipeline2.run.current_stage == 1
        assert 1 in pipeline2.run.stages_completed

    def test_find_latest_checkpoint(self, temp_output_dir, sample_config, mock_nim_manager, sample_target):
        """Save stages 0, 1, 2 — verify find_latest returns 2."""
        config = PipelineConfig(**{**sample_config.model_dump(), "output_dir": str(temp_output_dir)})
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()

        mgr = CheckpointManager(temp_output_dir)
        latest = mgr.find_latest_checkpoint(pipeline.run_id)
        assert latest == 2

    def test_find_latest_nonexistent(self, temp_output_dir):
        """No checkpoints returns None."""
        mgr = CheckpointManager(temp_output_dir)
        assert mgr.find_latest_checkpoint("nonexistent") is None

    def test_serialization_roundtrip(self, temp_output_dir):
        """Full state with molecules/docking results survives save/load."""
        mgr = CheckpointManager(temp_output_dir)

        # Create a mock pipeline with rich state
        pipeline = Mock()
        pipeline.run_id = "test123"
        pipeline.config = PipelineConfig(target_gene="VCP", output_dir=str(temp_output_dir))
        pipeline.run = PipelineRun(
            run_id="test123",
            target_gene="VCP",
            started_at="2025-01-01T00:00:00",
            current_stage=4,
            stages_completed=[0, 1, 2, 3, 4],
            molecules_generated=2,
        )
        pipeline.target = TargetHypothesis(
            gene="VCP",
            protein="p97",
            rationale="test",
            pdb_ids=["5FTK"],
        )
        pipeline.structures = StructureManifest(
            target_gene="VCP",
            structures=[StructureInfo(pdb_id="5FTK", method="X-ray", resolution=2.3)],
            primary_structure="5FTK",
        )
        pipeline.molecules = [
            GeneratedMolecule(
                id="mol-0001",
                smiles="CCO",
                target_gene="VCP",
                generation_score=0.9,
                properties=MoleculeProperties(
                    molecular_weight=46.07,
                    logP=-0.31,
                    hbd=1,
                    hba=1,
                    tpsa=20.23,
                    rotatable_bonds=0,
                    lipinski_violations=0,
                    qed=0.4,
                ),
            )
        ]
        pipeline.docking_results = [
            DockingResult(
                molecule_id="mol-0001",
                structure_id="5FTK",
                docking_score=-8.5,
                confidence=0.8,
                hydrogen_bonds=2,
            )
        ]
        pipeline.ranked_candidates = [
            RankedCandidate(
                rank=1,
                molecule_id="mol-0001",
                smiles="CCO",
                target_gene="VCP",
                docking_score=-8.5,
                generation_score=0.9,
                qed_score=0.4,
                composite_score=0.72,
            )
        ]

        # Save
        path = mgr.save(pipeline, 4)
        assert path.exists()

        # Load and restore
        checkpoint = mgr.load("test123", 4)
        assert checkpoint is not None

        restored = Mock()
        mgr.restore_pipeline_state(restored, checkpoint)

        assert restored.target.gene == "VCP"
        assert restored.structures.primary_structure == "5FTK"
        assert len(restored.molecules) == 1
        assert restored.molecules[0].smiles == "CCO"
        assert restored.molecules[0].properties.molecular_weight == 46.07
        assert len(restored.docking_results) == 1
        assert restored.docking_results[0].docking_score == -8.5
        assert len(restored.ranked_candidates) == 1
        assert restored.ranked_candidates[0].composite_score == 0.72
        assert restored.run.current_stage == 4


class TestPipelineResume:
    """Tests for pipeline resume-from-checkpoint functionality."""

    def test_resume_skips_completed_stages(self, temp_output_dir, mock_nim_manager, sample_target):
        """Resume from stage 2 should skip stages 0-2."""
        config = PipelineConfig(
            target_gene="VCP",
            num_molecules=2,
            output_dir=str(temp_output_dir),
        )

        # Run first pipeline through stage 2
        p1 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        p1.stage_0_initialize(sample_target)
        p1.stage_1_normalize_target()
        p1.stage_2_structure_discovery()
        run_id = p1.run_id

        # Resume with a new pipeline
        p2 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)

        with patch.object(p2, 'stage_0_initialize') as s0, \
             patch.object(p2, 'stage_1_normalize_target') as s1, \
             patch.object(p2, 'stage_2_structure_discovery') as s2, \
             patch.object(p2, 'stage_3_structure_prep') as s3, \
             patch.object(p2, 'stage_4_molecule_generation') as s4, \
             patch.object(p2, 'stage_5_chemistry_qc') as s5, \
             patch.object(p2, 'stage_6_conformers') as s6, \
             patch.object(p2, 'stage_7_docking') as s7, \
             patch.object(p2, 'stage_8_ranking') as s8, \
             patch.object(p2, 'stage_9_reporting') as s9:

            p2.run_pipeline(resume_from_run_id=run_id)

            # Stages 0-2 should NOT be called (already checkpointed)
            s0.assert_not_called()
            s1.assert_not_called()
            s2.assert_not_called()

            # Stages 3-9 SHOULD be called
            s3.assert_called_once()
            s4.assert_called_once()
            s5.assert_called_once()
            s6.assert_called_once()
            s7.assert_called_once()
            s8.assert_called_once()
            s9.assert_called_once()

    def test_resume_nonexistent_starts_fresh(self, temp_output_dir, mock_nim_manager, sample_target):
        """Invalid run_id should start a fresh run."""
        config = PipelineConfig(
            target_gene="VCP",
            num_molecules=2,
            output_dir=str(temp_output_dir),
        )
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)

        with patch.object(pipeline, 'stage_0_initialize') as s0:
            with patch.object(pipeline, 'stage_1_normalize_target'), \
                 patch.object(pipeline, 'stage_2_structure_discovery'), \
                 patch.object(pipeline, 'stage_3_structure_prep'), \
                 patch.object(pipeline, 'stage_4_molecule_generation'), \
                 patch.object(pipeline, 'stage_5_chemistry_qc'), \
                 patch.object(pipeline, 'stage_6_conformers'), \
                 patch.object(pipeline, 'stage_7_docking'), \
                 patch.object(pipeline, 'stage_8_ranking'), \
                 patch.object(pipeline, 'stage_9_reporting'):

                pipeline.run_pipeline(
                    target=sample_target,
                    resume_from_run_id="nonexistent-id",
                )

                # Stage 0 should be called (fresh start)
                s0.assert_called_once()
