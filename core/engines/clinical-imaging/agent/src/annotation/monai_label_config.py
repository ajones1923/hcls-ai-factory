"""MONAI Label integration for interactive annotation.

Connects MONAI Label's active learning annotation server with
the existing OHIF Viewer, VISTA-3D segmentation, and FLARE
federated learning pipeline.

Apache 2.0 Licensed. pip install monailabel.

Author: Adam Jones
Date: April 2026
"""

import json
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import requests
from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# Configuration Models
# ═══════════════════════════════════════════════════════════════════════


class AnnotationConfig(BaseModel):
    """MONAI Label server configuration."""

    host: str = "0.0.0.0"
    port: int = 8527
    app_dir: str = "data/monai_label_apps"
    studies_dir: str = "data/monai_label_studies"
    datastore_type: str = "orthanc"  # "orthanc" or "local"
    orthanc_url: str = "http://localhost:8042"
    models: List[str] = Field(
        default_factory=lambda: [
            "vista3d_interactive",
            "segmentation_spleen",
            "deepedit",
        ]
    )
    auto_segmentation: bool = True  # Auto-segment on study load
    active_learning: bool = True  # Active learning suggestions


class AnnotationSession(BaseModel):
    """Tracks an annotation session."""

    session_id: str = ""
    study_id: str = ""
    annotator: str = "anonymous"
    model_used: str = ""
    labels_created: int = 0
    labels_corrected: int = 0
    time_seconds: float = 0.0
    status: str = "active"  # "active", "completed", "abandoned"
    created_at: str = ""


class AnnotationStats(BaseModel):
    """Aggregate annotation statistics."""

    total_sessions: int = 0
    total_labels: int = 0
    total_corrections: int = 0
    correction_rate: float = 0.0  # corrections / total labels
    avg_session_time_seconds: float = 0.0
    models_used: Dict[str, int] = Field(default_factory=dict)  # model -> count
    annotators: Dict[str, int] = Field(default_factory=dict)  # annotator -> count
    studies_annotated: int = 0
    labels_for_flare: int = 0  # labels ready for federated training


# ═══════════════════════════════════════════════════════════════════════
# Default Model Catalog
# ═══════════════════════════════════════════════════════════════════════

_DEFAULT_MODELS: List[Dict] = [
    {
        "name": "vista3d_interactive",
        "type": "segmentation",
        "description": "VISTA-3D interactive segmentation (127+ classes)",
        "source": "nvidia",
    },
    {
        "name": "segmentation_spleen",
        "type": "segmentation",
        "description": "Spleen segmentation (MONAI bundle)",
        "source": "monai",
    },
    {
        "name": "deepedit",
        "type": "interactive",
        "description": "DeepEdit interactive segmentation",
        "source": "monai",
    },
    {
        "name": "segmentation_liver",
        "type": "segmentation",
        "description": "Liver and tumor segmentation",
        "source": "monai",
    },
]


# ═══════════════════════════════════════════════════════════════════════
# MONAILabelManager
# ═══════════════════════════════════════════════════════════════════════


class MONAILabelManager:
    """Manages MONAI Label annotation server and sessions.

    Data flywheel:
    clinical imaging -> AI inference -> radiologist review ->
    annotation correction -> FLARE federated training ->
    improved models -> better inference

    Works gracefully whether or not the ``monailabel`` package is
    installed or the MONAI Label server is reachable.
    """

    def __init__(self, config: Optional[AnnotationConfig] = None):
        self.config = config or AnnotationConfig()
        self._monailabel_available = self._check_monailabel()
        self._sessions: Dict[str, AnnotationSession] = {}
        self._completed_sessions: List[AnnotationSession] = []

    # ── Bootstrap checks ──────────────────────────────────────────────

    def _check_monailabel(self) -> bool:
        """Return True if the ``monailabel`` package is importable."""
        try:
            import monailabel  # noqa: F401

            logger.info("MONAI Label package detected")
            return True
        except ImportError:
            logger.info(
                "MONAI Label not installed. Install with: pip install monailabel"
            )
            return False

    # ── Server status ─────────────────────────────────────────────────

    def get_server_status(self) -> Dict:
        """Check if MONAI Label server is running.

        Returns:
            Dict with ``available`` (bool), ``version`` (str), and
            ``models`` (list) when the server is reachable.  Returns a
            degraded-status dict when it is not.
        """
        url = f"http://{self.config.host}:{self.config.port}/info"
        try:
            resp = requests.get(url, timeout=5)
            resp.raise_for_status()
            info = resp.json()
            return {
                "available": True,
                "version": info.get("version", "unknown"),
                "models": info.get("models", {}),
                "labels": info.get("labels", {}),
                "datastore": info.get("datastore", {}),
            }
        except Exception as exc:
            logger.debug(f"MONAI Label server not reachable: {exc}")
            return {
                "available": False,
                "version": None,
                "models": {},
                "labels": {},
                "datastore": {},
                "message": (
                    f"MONAI Label server not reachable at {url}. "
                    "Start with: docker compose up monai-label"
                ),
            }

    # ── Session management ────────────────────────────────────────────

    def start_session(
        self,
        study_id: str,
        annotator: str = "anonymous",
        model: str = "vista3d_interactive",
    ) -> AnnotationSession:
        """Start a new annotation session for a study.

        Args:
            study_id:  Orthanc study identifier.
            annotator: Name or ID of the person annotating.
            model:     MONAI Label model to use for pre-annotation.

        Returns:
            A new :class:`AnnotationSession` in ``active`` status.
        """
        session = AnnotationSession(
            session_id=str(uuid.uuid4()),
            study_id=study_id,
            annotator=annotator,
            model_used=model,
            labels_created=0,
            labels_corrected=0,
            time_seconds=0.0,
            status="active",
            created_at=datetime.now(timezone.utc).isoformat(),
        )
        self._sessions[session.session_id] = session
        logger.info(
            f"Annotation session started: {session.session_id} "
            f"(study={study_id}, annotator={annotator}, model={model})"
        )
        return session

    def end_session(
        self, session_id: str, status: str = "completed"
    ) -> AnnotationSession:
        """End an annotation session and record stats.

        Args:
            session_id: Session to close.
            status:     Terminal status (``completed`` or ``abandoned``).

        Returns:
            The updated :class:`AnnotationSession`.

        Raises:
            KeyError: If the session is not found.
        """
        session = self._sessions.get(session_id)
        if session is None:
            raise KeyError(f"Session not found: {session_id}")

        # Calculate elapsed time from created_at to now
        try:
            start = datetime.fromisoformat(session.created_at)
            elapsed = (datetime.now(timezone.utc) - start).total_seconds()
            session.time_seconds = round(elapsed, 2)
        except (ValueError, TypeError):
            session.time_seconds = 0.0

        session.status = status

        # Move to completed list
        self._completed_sessions.append(session)
        del self._sessions[session_id]

        logger.info(
            f"Annotation session ended: {session_id} "
            f"(status={status}, labels={session.labels_created}, "
            f"corrections={session.labels_corrected}, "
            f"time={session.time_seconds:.1f}s)"
        )
        return session

    # ── Label tracking ────────────────────────────────────────────────

    def record_label(self, session_id: str, label_type: str = "created") -> bool:
        """Record a label creation or correction in the active session.

        Args:
            session_id: Active session ID.
            label_type: ``"created"`` or ``"corrected"``.

        Returns:
            True on success, False if the session was not found or
            the label_type is invalid.
        """
        session = self._sessions.get(session_id)
        if session is None:
            logger.warning(f"Cannot record label: session {session_id} not found")
            return False

        if label_type == "created":
            session.labels_created += 1
        elif label_type == "corrected":
            session.labels_corrected += 1
        else:
            logger.warning(f"Invalid label_type: {label_type}")
            return False

        logger.debug(
            f"Label recorded ({label_type}) in session {session_id}: "
            f"created={session.labels_created}, corrected={session.labels_corrected}"
        )
        return True

    # ── Session queries ───────────────────────────────────────────────

    def get_session(self, session_id: str) -> Optional[AnnotationSession]:
        """Get session details by ID.

        Searches both active and completed sessions.

        Returns:
            The session, or ``None`` if not found.
        """
        session = self._sessions.get(session_id)
        if session is not None:
            return session

        # Check completed sessions
        for s in self._completed_sessions:
            if s.session_id == session_id:
                return s

        return None

    # ── Aggregate statistics ──────────────────────────────────────────

    def get_stats(self) -> AnnotationStats:
        """Get aggregate annotation statistics across all sessions."""
        all_sessions = list(self._sessions.values()) + self._completed_sessions

        if not all_sessions:
            return AnnotationStats()

        total_labels = sum(s.labels_created for s in all_sessions)
        total_corrections = sum(s.labels_corrected for s in all_sessions)
        completed = [s for s in all_sessions if s.status == "completed"]
        total_time = sum(s.time_seconds for s in completed)

        # Model usage counts
        models_used: Dict[str, int] = {}
        for s in all_sessions:
            if s.model_used:
                models_used[s.model_used] = models_used.get(s.model_used, 0) + 1

        # Annotator counts
        annotators: Dict[str, int] = {}
        for s in all_sessions:
            annotators[s.annotator] = annotators.get(s.annotator, 0) + 1

        # Unique studies
        studies_annotated = len({s.study_id for s in all_sessions})

        # Correction rate (avoid division by zero)
        total_all_labels = total_labels + total_corrections
        correction_rate = (
            total_corrections / total_all_labels if total_all_labels > 0 else 0.0
        )

        # Labels ready for FLARE: all labels from completed sessions
        labels_for_flare = sum(
            s.labels_created + s.labels_corrected
            for s in completed
        )

        return AnnotationStats(
            total_sessions=len(all_sessions),
            total_labels=total_labels,
            total_corrections=total_corrections,
            correction_rate=round(correction_rate, 4),
            avg_session_time_seconds=(
                round(total_time / len(completed), 2) if completed else 0.0
            ),
            models_used=models_used,
            annotators=annotators,
            studies_annotated=studies_annotated,
            labels_for_flare=labels_for_flare,
        )

    # ── Model catalog ─────────────────────────────────────────────────

    def get_available_models(self) -> List[Dict]:
        """List available MONAI Label models.

        Queries the running MONAI Label server first.  Falls back to a
        preset catalog when the server is unavailable.

        Returns:
            List of model descriptor dicts.
        """
        # Try live server first
        status = self.get_server_status()
        if status.get("available"):
            server_models = status.get("models", {})
            if server_models:
                return [
                    {
                        "name": name,
                        "type": info.get("type", "unknown"),
                        "description": info.get("description", ""),
                        "source": "server",
                    }
                    for name, info in server_models.items()
                ]

        # Fallback: return preset catalog
        logger.debug("Returning preset model catalog (server unavailable)")
        return list(_DEFAULT_MODELS)

    # ── FLARE export ──────────────────────────────────────────────────

    def export_for_flare(self, output_dir: str = "data/flare_training") -> Dict:
        """Export annotated data in format suitable for NVIDIA FLARE training.

        Generates a manifest file summarizing completed annotation sessions
        and their label counts, grouped by model.  The manifest can be
        consumed by the FLARE job configurations in ``flare/``.

        Args:
            output_dir: Directory where the manifest will be written.

        Returns:
            Dict with ``total_labels``, ``by_model`` breakdown, and
            ``manifest_path``.
        """
        completed = [
            s for s in self._completed_sessions if s.status == "completed"
        ]

        # Aggregate by model
        by_model: Dict[str, int] = {}
        total_labels = 0
        for s in completed:
            count = s.labels_created + s.labels_corrected
            total_labels += count
            by_model[s.model_used] = by_model.get(s.model_used, 0) + count

        # Build manifest
        manifest = {
            "export_timestamp": datetime.now(timezone.utc).isoformat(),
            "total_sessions": len(completed),
            "total_labels": total_labels,
            "by_model": by_model,
            "studies": list({s.study_id for s in completed}),
        }

        # Write manifest to disk
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        manifest_path = out_path / "annotation_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2))

        logger.info(
            f"FLARE export manifest written to {manifest_path}: "
            f"{total_labels} labels from {len(completed)} sessions"
        )

        return {
            "total_labels": total_labels,
            "by_model": by_model,
            "manifest_path": str(manifest_path),
        }

    # ── OHIF integration ──────────────────────────────────────────────

    def get_ohif_annotation_url(self, study_id: str) -> str:
        """Generate OHIF Viewer URL with MONAI Label extension enabled.

        Args:
            study_id: Orthanc study identifier.

        Returns:
            URL string that opens OHIF with the annotation panel active
            and the MONAI Label server pre-configured.
        """
        base = self.config.orthanc_url.rstrip("/")
        # OHIF v3 URL pattern with MONAI Label extension
        ohif_url = (
            f"http://localhost:8526"
            f"/viewer?StudyInstanceUIDs={study_id}"
            f"&hangingprotocolId=monaiLabel"
        )
        logger.debug(f"OHIF annotation URL: {ohif_url}")
        return ohif_url
