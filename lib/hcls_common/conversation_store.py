"""Session-memory persistence for the RAG engines.

Five subjects carried a byte-identical copy of `_save_conversation` / `_load_conversation` /
`_cleanup_expired_conversations` — 43 lines each, differing in nothing at all, and with no test
anywhere in the repo covering any of them. Session memory that silently stops persisting looks
exactly like a working agent that has forgotten the conversation, which is the hardest kind of
fault to notice from the outside.

The behaviour is preserved exactly, including that every operation swallows its exception and
warns: a conversation that fails to persist must never take down the answer the user asked for.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_TTL = timedelta(hours=24)


class ConversationStore:
    """JSON-per-session conversation history on disk, with a TTL."""

    def __init__(self, directory: Path | str, ttl: timedelta = DEFAULT_TTL):
        self.directory = Path(directory)
        self.ttl = ttl

    def _path(self, session_id: str) -> Path:
        return self.directory / f"{session_id}.json"

    def save(self, session_id: str, history: list) -> None:
        """Persist conversation to disk as JSON."""
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            self._path(session_id).write_text(json.dumps({
                "session_id": session_id,
                "updated": datetime.now(timezone.utc).isoformat(),
                "messages": history,
            }, indent=2))
        except Exception as exc:
            logger.warning("Failed to persist conversation %s: %s", session_id, exc)

    def load(self, session_id: str) -> list:
        """Load conversation from disk, respecting the TTL. Expired files are removed."""
        try:
            path = self._path(session_id)
            if path.exists():
                data = json.loads(path.read_text())
                updated = datetime.fromisoformat(data["updated"])
                if datetime.now(timezone.utc) - updated < self.ttl:
                    return data.get("messages", [])
                path.unlink(missing_ok=True)          # expired
        except Exception as exc:
            logger.warning("Failed to load conversation %s: %s", session_id, exc)
        return []

    def cleanup_expired(self) -> None:
        """Remove every conversation file older than the TTL."""
        try:
            if not self.directory.exists():
                return
            cutoff = datetime.now(timezone.utc) - self.ttl
            for path in self.directory.glob("*.json"):
                try:
                    if datetime.fromisoformat(json.loads(path.read_text())["updated"]) < cutoff:
                        path.unlink()
                except Exception:
                    pass                              # one corrupt file must not stop the sweep
        except Exception as exc:
            logger.warning("Conversation cleanup error: %s", exc)
