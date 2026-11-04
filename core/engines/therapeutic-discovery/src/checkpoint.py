"""
Pipeline Checkpoint Manager.

Saves and restores pipeline state between stages to enable
resume-after-failure capability.
"""
import json
from pathlib import Path
from typing import Optional, List, TYPE_CHECKING
from datetime import datetime
from loguru import logger

from .models import (
    TargetHypothesis,
    StructureManifest,
    GeneratedMolecule,
    DockingResult,
    RankedCandidate,
    PipelineRun,
)

if TYPE_CHECKING:
    from .pipeline import DrugDiscoveryPipeline


def _timings(pipeline) -> dict:
    """The pipeline's stage timings as a plain dict, or {} if it has none we can read."""
    timings = getattr(pipeline, "stage_timings", None)
    if isinstance(timings, dict):
        return dict(timings)
    return {}


class ConfigMismatchError(RuntimeError):
    """A resume would mix results computed under two different configurations."""


class CheckpointManager:
    """Manages checkpoint save/load for pipeline runs."""

    def __init__(self, output_dir: Path):
        self.checkpoint_dir = output_dir / "checkpoints"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def _checkpoint_path(self, run_id: str, stage: int) -> Path:
        return self.checkpoint_dir / f"checkpoint_{run_id}_stage{stage}.json"

    def save(self, pipeline: "DrugDiscoveryPipeline", stage: int) -> Path:
        """Save pipeline state after a stage completes."""
        checkpoint = {
            "run_id": pipeline.run_id,
            "stage": stage,
            "saved_at": datetime.now().isoformat(),
            "config": pipeline.config.model_dump(),
            "run": pipeline.run.model_dump(),
            "target": pipeline.target.model_dump() if pipeline.target else None,
            "structures": pipeline.structures.model_dump() if pipeline.structures else None,
            "molecules": [m.model_dump() for m in pipeline.molecules],
            "docking_results": [d.model_dump() for d in pipeline.docking_results],
            "ranked_candidates": [r.model_dump() for r in pipeline.ranked_candidates],
            # Timings were saved nowhere, so after a resume `run.stage_timings` held only the
            # stages that re-ran while `run.status` was set to "completed" -- an "end-to-end in
            # N minutes" figure read off that report understates the real cost by every skipped
            # stage.
            "stage_timings": _timings(pipeline),
        }

        path = self._checkpoint_path(pipeline.run_id, stage)
        with open(path, "w") as f:
            json.dump(checkpoint, f, indent=2, default=str)

        logger.info(f"Checkpoint saved: stage {stage} -> {path.name}")
        return path

    def load(self, run_id: str, stage: int) -> Optional[dict]:
        """Load a checkpoint file, returning raw dict."""
        path = self._checkpoint_path(run_id, stage)
        if not path.exists():
            return None

        with open(path) as f:
            return json.load(f)

    def find_latest_checkpoint(self, run_id: str) -> Optional[int]:
        """Find the highest completed stage for a run_id."""
        latest = -1
        for f in self.checkpoint_dir.glob(f"checkpoint_{run_id}_stage*.json"):
            try:
                stage_num = int(f.stem.split("stage")[-1])
                if stage_num > latest:
                    latest = stage_num
            except ValueError:
                continue
        return latest if latest >= 0 else None

    def restore_pipeline_state(
        self, pipeline: "DrugDiscoveryPipeline", checkpoint: dict,
        allow_config_change: bool = False,
    ):
        """Restore pipeline state from a checkpoint dict.

        `save()` stores `config`, and this method never read it back. So a resume with a
        different `PipelineConfig` silently mixed: stages at or below `resume_stage` kept results
        computed under the OLD config, every later stage used the NEW one, and
        `stage_9_reporting` wrote the new config into the report. Run to completion with
        `num_molecules=20, max_mw=550, docking_weight=0.4`, resume the same run_id with
        `max_mw=400, docking_weight=0.1`, and nothing re-executes -- yet the report documents
        parameters that did not produce its own numbers. Checkpoints are keyed on run_id and
        stage only, with no hash of inputs or config, so a stage-8 result computed days earlier
        is served as current.

        The mismatch is now refused. `allow_config_change=True` accepts it deliberately and the
        pipeline records the divergence so the report cannot present it as a single coherent run.
        """
        saved_config = checkpoint.get("config")
        current_config = None
        config = getattr(pipeline, "config", None)
        if config is not None:
            try:
                current_config = config.model_dump()
            except Exception:                       # not a PipelineConfig; nothing to compare
                current_config = None
        if not isinstance(saved_config, dict) or not isinstance(current_config, dict):
            saved_config = current_config = None
        if saved_config is not None and saved_config != current_config:
            differing = sorted(
                k for k in set(saved_config) | set(current_config)
                if saved_config.get(k) != current_config.get(k)
            )
            if not allow_config_change:
                raise ConfigMismatchError(
                    f"Refusing to resume run {checkpoint.get('run_id')} from stage "
                    f"{checkpoint.get('stage')}: the checkpoint was computed under a different "
                    f"configuration. Differing keys: {', '.join(differing)}. The restored stages "
                    "would keep results produced by the old settings while the report documented "
                    "the new ones. Re-run from scratch, or pass allow_config_change=True to "
                    "accept a mixed run (it will be labelled as one)."
                )
            pipeline.config_divergence = {
                "resumed_under_different_config": True,
                "differing_keys": differing,
                "checkpoint_config": saved_config,
            }
            logger.warning(
                f"Resuming under a CHANGED config; differing keys: {', '.join(differing)}. "
                "Stages restored from the checkpoint were computed under the old settings."
            )

        pipeline.stage_timings = {
            int(k): v for k, v in (checkpoint.get("stage_timings") or {}).items()
        }
        pipeline.restored_stages = sorted(pipeline.stage_timings)

        if checkpoint.get("target"):
            pipeline.target = TargetHypothesis.model_validate(checkpoint["target"])

        if checkpoint.get("structures"):
            pipeline.structures = StructureManifest.model_validate(
                checkpoint["structures"]
            )

        pipeline.molecules = [
            GeneratedMolecule.model_validate(m)
            for m in checkpoint.get("molecules", [])
        ]
        pipeline.docking_results = [
            DockingResult.model_validate(d)
            for d in checkpoint.get("docking_results", [])
        ]
        pipeline.ranked_candidates = [
            RankedCandidate.model_validate(r)
            for r in checkpoint.get("ranked_candidates", [])
        ]

        pipeline.run = PipelineRun.model_validate(checkpoint["run"])
