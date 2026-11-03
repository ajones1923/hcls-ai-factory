# Copyright 2026 Adam Jones
# SPDX-License-Identifier: Apache-2.0
# Part of the HCLS AI Factory — https://github.com/ajones1923/hcls-ai-factory
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
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

import hashlib
import json
import logging
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_TTL = timedelta(hours=24)

# `session_id` arrives in a request payload and used to be interpolated straight into a
# filename, so an id of "../../etc/cron.d/pwn" wrote JSON outside the conversation directory
# and an id of "../../../some/other.json" read a file that is not a conversation. Every
# operation here swallows its exception by design, so both failed silently.
_SAFE_SESSION_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def safe_session_name(session_id: str) -> str:
    """A filesystem-safe stem for a session id, preserving distinctness.

    An id that is already safe is used as-is, so existing conversation files keep their names.
    Anything else is replaced by a digest of the original rather than by a stripped version of
    it: two ids that sanitized to the same name would serve one session's history to another.
    """
    sid = "" if session_id is None else str(session_id)
    if _SAFE_SESSION_ID.match(sid):
        return sid
    digest = hashlib.sha256(sid.encode("utf-8", "surrogatepass")).hexdigest()[:32]
    return f"sid-{digest}"


class ConversationStore:
    """JSON-per-session conversation history on disk, with a TTL."""

    def __init__(self, directory: Path | str, ttl: timedelta = DEFAULT_TTL):
        self.directory = Path(directory)
        self.ttl = ttl

    def _path(self, session_id: str) -> Path:
        path = self.directory / f"{safe_session_name(session_id)}.json"
        # Belt and braces: the name is already constrained above, and this makes the
        # containment the actual invariant rather than a property of the regex.
        root = self.directory.resolve()
        resolved = path.resolve()
        if resolved.parent != root:
            raise ValueError(f"session id escapes the conversation directory: {session_id!r}")
        return path

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
