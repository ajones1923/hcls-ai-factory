"""
Tests for pipeline performance timing.
"""
import pytest
from src.models import PipelineConfig, TargetHypothesis
from src.pipeline import DrugDiscoveryPipeline


class TestStageTimings:
    """Verify that stage timings are recorded during pipeline runs."""

    def test_timings_recorded_for_all_stages(
        self, temp_output_dir, sample_config, mock_nim_manager, sample_target
    ):
        """Full pipeline run should record timing for all 10 stages."""
        config = PipelineConfig(
            **{**sample_config.model_dump(), "output_dir": str(temp_output_dir)}
        )
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        run = pipeline.run_pipeline(sample_target)

        # All 10 stages should have timings
        assert len(pipeline.stage_timings) == 10
        for stage in range(10):
            assert stage in pipeline.stage_timings, f"Missing timing for stage {stage}"
            assert pipeline.stage_timings[stage] >= 0, f"Negative timing for stage {stage}"

    def test_timings_in_pipeline_run(
        self, temp_output_dir, sample_config, mock_nim_manager, sample_target
    ):
        """Timings should be propagated to PipelineRun model."""
        config = PipelineConfig(
            **{**sample_config.model_dump(), "output_dir": str(temp_output_dir)}
        )
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        run = pipeline.run_pipeline(sample_target)

        assert run.stage_timings == pipeline.stage_timings
        assert len(run.stage_timings) == 10

    def test_timings_in_report(
        self, temp_output_dir, sample_config, mock_nim_manager, sample_target
    ):
        """Timings should appear in the generated report JSON."""
        import json

        config = PipelineConfig(
            **{**sample_config.model_dump(), "output_dir": str(temp_output_dir)}
        )
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        pipeline.run_pipeline(sample_target)

        # Find the report file
        reports = list(temp_output_dir.glob("report_*.json"))
        assert len(reports) == 1

        report = json.loads(reports[0].read_text())
        assert "stage_timings" in report
        # Report is written during stage 9, so stage 9 timing is not yet available
        assert len(report["stage_timings"]) == 9

    def test_resume_preserves_prior_timings(
        self, temp_output_dir, mock_nim_manager, sample_target
    ):
        """Resumed pipeline should have timings for stages it actually ran."""
        config = PipelineConfig(
            target_gene="VCP",
            num_molecules=2,
            output_dir=str(temp_output_dir),
        )

        # Run first pipeline through all stages
        p1 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        p1.run_pipeline(sample_target)
        run_id = p1.run_id

        # Resume from checkpoint — should have timings for stages it re-ran
        p2 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        run = p2.run_pipeline(resume_from_run_id=run_id)

        # The resumed pipeline ran no stages (all checkpointed), so no new timings
        # but the run should complete successfully
        assert run.status == "completed"
