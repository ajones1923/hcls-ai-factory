"""Regression tests: a resume must not serve results computed under a different config as fresh.

Audit compute#30 / C13b. `CheckpointManager.save` stores `config`, and `restore_pipeline_state`
never read it back -- it restored target, structures, molecules, docking_results,
ranked_candidates and run, and nothing else. So a resume with a different `PipelineConfig`
silently mixed: stages at or below `resume_stage` kept results computed under the OLD config,
every later stage used the NEW one, and `stage_9_reporting` wrote
`"config": self.config.dict()` -- the new config -- into the report.

Run the VCP pipeline to completion with `num_molecules=20, max_mw=550, docking_weight=0.4`, then
resume the same run_id with `max_mw=400, docking_weight=0.1`: stages 0-9 are all at or below
`resume_stage=9`, so nothing re-executes, yet the report is written with the new parameters
beside candidates filtered at 550 Da and ranked with a 0.4 docking weight. The artifact documents
parameters that did not produce its own numbers.

`stage_timings` was likewise neither saved nor restored, so after a resume `run.stage_timings`
held only the stages that actually re-ran while `run.status` was `"completed"` -- any
"end-to-end in N minutes" figure read off that report understates the real cost by every skipped
stage. `tests/test_checkpoint.py` built the same config on both sides of every resume, so the
gap was untested.
"""
import json
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from src.checkpoint import ConfigMismatchError  # noqa: E402
from src.models import PipelineConfig  # noqa: E402
from src.pipeline import DrugDiscoveryPipeline  # noqa: E402


def _config(output_dir, structure_cache, **overrides):
    base = dict(target_gene="VCP", num_molecules=2, output_dir=str(output_dir),
                structure_cache=str(structure_cache))
    base.update(overrides)
    return PipelineConfig(**base)


class TestAChangedConfigIsRefused:
    def test_resuming_under_a_different_config_raises(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        p1 = DrugDiscoveryPipeline(_config(temp_output_dir, structure_cache, max_mw=550),
                                   nim_manager=mock_nim_manager)
        p1.run_pipeline(sample_target)

        p2 = DrugDiscoveryPipeline(_config(temp_output_dir, structure_cache, max_mw=400),
                                   nim_manager=mock_nim_manager)
        with pytest.raises(ConfigMismatchError) as e:
            p2.run_pipeline(resume_from_run_id=p1.run_id)
        assert "max_mw" in str(e.value)

    def test_the_same_config_still_resumes(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        config = _config(temp_output_dir, structure_cache)
        p1 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        p1.run_pipeline(sample_target)
        p2 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        run = p2.run_pipeline(resume_from_run_id=p1.run_id)
        assert run.status == "completed"

    def test_an_explicit_override_is_recorded_as_a_divergence(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        p1 = DrugDiscoveryPipeline(_config(temp_output_dir, structure_cache, max_mw=550),
                                   nim_manager=mock_nim_manager)
        p1.run_pipeline(sample_target)

        p2 = DrugDiscoveryPipeline(_config(temp_output_dir, structure_cache, max_mw=400),
                                   nim_manager=mock_nim_manager)
        p2.run_pipeline(resume_from_run_id=p1.run_id, allow_config_change=True)

        assert p2.config_divergence["resumed_under_different_config"] is True
        assert "max_mw" in p2.config_divergence["differing_keys"]


class TestTimingsCannotBeReadAsEndToEnd:
    @staticmethod
    def _resume_from_stage_5(output_dir, structure_cache, mock_nim_manager, sample_target):
        """Run to completion, then drop the stage 6-9 checkpoints so the resume re-runs 5-9."""
        config = _config(output_dir, structure_cache)
        p1 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        p1.run_pipeline(sample_target)
        for stage in range(6, 10):
            (Path(output_dir) / "checkpoints"
             / f"checkpoint_{p1.run_id}_stage{stage}.json").unlink()
        for report in Path(output_dir).glob("report_*.json"):
            report.unlink()

        p2 = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        p2.run_pipeline(resume_from_run_id=p1.run_id)
        return p1, p2

    def test_timings_from_the_earlier_invocation_survive_the_resume(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        p1, p2 = self._resume_from_stage_5(
            temp_output_dir, structure_cache, mock_nim_manager, sample_target)
        # The stage-5 checkpoint was written inside stage 5, so it carries timings for 0-4.
        assert set(p2.restored_stages) == {0, 1, 2, 3, 4}
        for stage in p2.restored_stages:
            assert p2.stage_timings[stage] == p1.stage_timings[stage]
        # The resume timed the stages it actually ran (6-9). Stage 5 has no timing in either
        # invocation: its checkpoint is written inside the stage, before run_pipeline records
        # the duration, and the resume skips it -- so the report names it under
        # `stages_with_no_timing` rather than letting the remaining timings be summed as a
        # complete total.
        assert set(p2.stage_timings) == {0, 1, 2, 3, 4, 6, 7, 8, 9}

    def test_a_resumed_report_says_it_is_not_end_to_end(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        self._resume_from_stage_5(
            temp_output_dir, structure_cache, mock_nim_manager, sample_target)
        report = json.loads(
            next(Path(temp_output_dir).glob("report_*.json")).read_text())
        prov = report["timing_provenance"]
        assert prov["is_end_to_end"] is False
        assert prov["stages_restored_from_checkpoint"] == [0, 1, 2, 3, 4]
        # Stage 9 writes this report, so its own duration is not known yet.
        assert prov["stages_executed_now"] == [6, 7, 8]
        # Stage 5 ran in the earlier invocation and was skipped by the resume: no timing exists.
        assert 5 in prov["stages_with_no_timing"]

    def test_a_fresh_report_says_it_is_end_to_end(
        self, temp_output_dir, structure_cache, mock_nim_manager, sample_target
    ):
        p = DrugDiscoveryPipeline(_config(temp_output_dir, structure_cache),
                                  nim_manager=mock_nim_manager)
        p.run_pipeline(sample_target)
        report = json.loads(
            next(Path(temp_output_dir).glob("report_*.json")).read_text())
        prov = report["timing_provenance"]
        assert prov["is_end_to_end"] is True
        assert prov["stages_restored_from_checkpoint"] == []
        # Nothing ran outside this invocation: every completed stage has a timing (stage 9 is
        # marked complete only after this report is written).
        assert prov["stages_with_no_timing"] == []
        assert report["config_divergence"] is None
