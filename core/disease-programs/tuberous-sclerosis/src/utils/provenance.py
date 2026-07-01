"""
Provenance envelope (PRD §2.5.1, NFR-R-1; master paper §16).

Every agent invocation produces a structured provenance record: agent id/version,
model id/version, prompt-template version, retrieved RAG source URIs, input hash,
and latency. Append-only and queryable; provenance is the substrate, not a bolt-on.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

from pydantic import BaseModel, Field


def stable_hash(obj) -> str:
    """Deterministic SHA-256 of any JSON-serializable object."""
    blob = json.dumps(obj, sort_keys=True, default=str, ensure_ascii=False)
    return "sha256:" + hashlib.sha256(blob.encode("utf-8")).hexdigest()


class Provenance(BaseModel):
    agent: str
    agent_version: str = "0.1.0"
    step: str | None = None
    tier: str | None = None                       # haiku | sonnet | opus | classical | deterministic
    model_id: str | None = None
    model_version: str | None = None
    prompt_template_version: str | None = None
    rag_source_uris: list[str] = Field(default_factory=list)
    input_hash: str | None = None
    latency_ms: float | None = None
    watermark: str = "SYNTHETIC"
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    def to_dict(self) -> dict:
        return self.model_dump(mode="json")
