"""DICOM watcher service for Orthanc change tracking.

Polls the Orthanc ``/changes`` API at a configurable interval and
invokes a callback function whenever a ``StableStudy`` event is detected
(i.e. a study has finished receiving all instances and is stable).

Usage::

    from src.ingest.dicom_watcher import DicomWatcher

    def on_study(change: dict):
        print(f"New study: {change}")

    watcher = DicomWatcher(
        orthanc_url="http://localhost:8042",
        callback=on_study,
        poll_interval=5,
    )
    watcher.start()
    # ... later ...
    watcher.stop()

Author: Adam Jones
Date: February 2026
"""

import threading
import time
from typing import Callable, Optional

import httpx
from loguru import logger


class DicomWatcher:
    """Watches Orthanc for new stable studies via the /changes API.

    The Orthanc ``/changes`` endpoint returns an ordered stream of server
    events.  This watcher tracks the ``last_change`` sequence number so
    it only processes new events on each poll cycle.

    Attributes:
        orthanc_url: Base URL for the Orthanc REST API.
        callback: Callable invoked with each ``StableStudy`` change dict.
        poll_interval: Seconds between consecutive polls.
        username: HTTP basic-auth username (default ``admin``).
        password: HTTP basic-auth password (default ``orthanc``).
    """

    # Orthanc change types that represent a completed study
    STUDY_STABLE_TYPE = "StableStudy"

    def __init__(
        self,
        orthanc_url: str,
        callback: Callable[[dict], None],
        poll_interval: int = 5,
        username: str = "admin",
        password: str = "orthanc",
    ):
        self.orthanc_url = orthanc_url.rstrip("/")
        self.callback = callback
        self.poll_interval = poll_interval
        self.username = username
        self.password = password

        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._last_change: int = 0
        self._lock = threading.Lock()

    # ── Public API ────────────────────────────────────────────────

    @property
    def is_running(self) -> bool:
        """Return True if the watcher polling loop is active."""
        return self._running

    def start(self) -> None:
        """Start the background polling thread.

        Raises:
            RuntimeError: If the watcher is already running.
        """
        if self._running:
            raise RuntimeError("DicomWatcher is already running")

        self._running = True
        self._thread = threading.Thread(
            target=self._run_loop,
            name="dicom-watcher",
            daemon=True,
        )
        self._thread.start()
        logger.info(
            f"DicomWatcher started (url={self.orthanc_url}, "
            f"interval={self.poll_interval}s)"
        )

    def stop(self) -> None:
        """Signal the polling loop to stop and wait for the thread."""
        self._running = False
        if self._thread and self._thread.is_alive():
            self._thread.join(timeout=self.poll_interval + 2)
        logger.info("DicomWatcher stopped")

    # ── Internal ──────────────────────────────────────────────────

    def _run_loop(self) -> None:
        """Main polling loop executed in a daemon thread."""
        while self._running:
            try:
                self._poll_changes()
            except Exception as exc:
                logger.warning(f"DicomWatcher poll error: {exc}")
            time.sleep(self.poll_interval)

    def _poll_changes(self) -> None:
        """Poll Orthanc /changes API for new StableStudy events.

        The ``since`` parameter tells Orthanc to return only changes
        after the given sequence number.  After each successful poll
        the sequence cursor is advanced so the same events are never
        re-processed.
        """
        url = f"{self.orthanc_url}/changes"
        params = {"since": self._last_change, "limit": 100}
        auth = (self.username, self.password)

        with httpx.Client(timeout=15.0) as client:
            response = client.get(url, params=params, auth=auth)
            response.raise_for_status()
            data = response.json()

        changes = data.get("Changes", [])
        done = data.get("Done", True)
        last = data.get("Last", self._last_change)

        for change in changes:
            change_type = change.get("ChangeType", "")
            if change_type == self.STUDY_STABLE_TYPE:
                resource_id = change.get("ID", "")
                logger.info(f"StableStudy detected: {resource_id}")
                try:
                    self.callback(change)
                except Exception as exc:
                    logger.error(
                        f"DicomWatcher callback failed for {resource_id}: {exc}"
                    )

        with self._lock:
            self._last_change = last

        # If there are more changes, poll again immediately
        if not done:
            self._poll_changes()
