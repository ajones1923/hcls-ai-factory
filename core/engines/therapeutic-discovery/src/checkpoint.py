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
        self, pipeline: "DrugDiscoveryPipeline", checkpoint: dict
    ):
        """Restore pipeline state from a checkpoint dict."""
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
