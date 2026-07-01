"""
Event-sourced state (PRD §2.5). Append-only event log + materialized projections.

You can replay the log and reproduce any surface exactly; provenance is the
substrate. Two backends share one interface and one projection-folding function:
  - SqliteEventStore  — stdlib only; the default, runs on the Spark with no services.
  - PostgresEventStore — the production target (PRD §2.5.1), psycopg2.

`make_store()` selects Postgres when TSC_POSTGRES_DSN is set, else SQLite.
"""
from __future__ import annotations

import json
import sqlite3
import threading
from abc import ABC, abstractmethod
from pathlib import Path

from config.settings import settings
from src.orchestrator.events import AGENT_EMITS, Event, EventType

# event_type -> projection section it updates (PRD §2.5.2)
_SECTION_FOR_EVENT: dict[EventType, str] = {
    EventType.phenome_mapped: "hpo_profile",
    EventType.variant_curated: "variant_interp",
    EventType.trajectory_forecast: "trajectory",
    EventType.tand_surveyed: "tand_briefing",
    EventType.therapeutics_briefed: "therapeutics",
}
# section name -> agent name (for mapping agent_failed back to the projection section)
_AGENT_FOR_SECTION = {
    _SECTION_FOR_EVENT[ev]: agent
    for agent, ev in AGENT_EMITS.items()
    if ev in _SECTION_FOR_EVENT
}


def materialize_projection(events: list[Event]) -> dict:
    """Fold the event log for one patient into proj_patient_current (pure)."""
    proj: dict = {
        "patient_id": events[0].patient_id if events else None,
        "hpo_profile": None,
        "variant_interp": None,
        "trajectory": None,
        "tand_briefing": None,
        "therapeutics": None,
        "staleness": {},
        "last_event_id": None,
    }
    for ev in events:
        proj["last_event_id"] = ev.event_id
        section = _SECTION_FOR_EVENT.get(ev.event_type)
        if section is not None:
            proj[section] = ev.payload
            proj["staleness"][section] = {
                "status": "ok",
                "source_event_id": ev.event_id,
                "ts": ev.created_at.isoformat(),
            }
        elif ev.event_type == EventType.agent_failed:
            agent = ev.payload.get("agent")
            failed_section = next(
                (s for s, a in _AGENT_FOR_SECTION.items() if a == agent), None
            )
            if failed_section:
                # conservative failure: keep prior value, mark pending (PRD §2.5.2)
                proj["staleness"][failed_section] = {
                    "status": "pending",
                    "source_event_id": ev.event_id,
                    "ts": ev.created_at.isoformat(),
                    "error": ev.payload.get("error"),
                }
    return proj


class EventStore(ABC):
    @abstractmethod
    def append(self, event: Event) -> int: ...

    @abstractmethod
    def events_for(self, patient_id: str) -> list[Event]: ...

    def projection(self, patient_id: str) -> dict:
        return materialize_projection(self.events_for(patient_id))


class SqliteEventStore(EventStore):
    def __init__(self, path: str | Path | None = None) -> None:
        if path is None:
            settings.ensure_dirs()
            path = settings.SQLITE_PATH
        self.path = str(path)
        if self.path != ":memory:":
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._lock = threading.Lock()   # serialize access across FastAPI threadpool threads
        self._init()

    def _init(self) -> None:
        self._conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS engine_events (
                event_id     INTEGER PRIMARY KEY AUTOINCREMENT,
                patient_id   TEXT NOT NULL,
                event_type   TEXT NOT NULL,
                payload      TEXT NOT NULL,
                provenance   TEXT NOT NULL,
                created_at   TEXT NOT NULL,
                parent_event INTEGER REFERENCES engine_events(event_id)
            );
            CREATE INDEX IF NOT EXISTS ix_events_patient ON engine_events(patient_id, event_id);
            CREATE INDEX IF NOT EXISTS ix_events_type ON engine_events(event_type, created_at);
            """
        )
        self._conn.commit()

    def append(self, event: Event) -> int:
        with self._lock:
            cur = self._conn.execute(
                "INSERT INTO engine_events (patient_id, event_type, payload, provenance, created_at, parent_event)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (
                    event.patient_id,
                    event.event_type.value,
                    json.dumps(event.payload, default=str),
                    json.dumps(event.provenance, default=str),
                    event.created_at.isoformat(),
                    event.parent_event,
                ),
            )
            self._conn.commit()
            event.event_id = int(cur.lastrowid)
        return event.event_id

    def events_for(self, patient_id: str) -> list[Event]:
        with self._lock:
            rows = self._conn.execute(
                "SELECT * FROM engine_events WHERE patient_id = ? ORDER BY event_id",
                (patient_id,),
            ).fetchall()
        return [
            Event(
                event_id=r["event_id"],
                patient_id=r["patient_id"],
                event_type=EventType(r["event_type"]),
                payload=json.loads(r["payload"]),
                provenance=json.loads(r["provenance"]),
                created_at=r["created_at"],
                parent_event=r["parent_event"],
            )
            for r in rows
        ]


class PostgresEventStore(EventStore):  # pragma: no cover - prod target, needs PG
    """Production backend (PRD §2.5.1). Same interface; requires a running PostgreSQL."""

    def __init__(self, dsn: str | None = None) -> None:
        import psycopg2  # noqa: WPS433

        self._conn = psycopg2.connect(dsn or settings.POSTGRES_DSN)
        with self._conn, self._conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS engine_events (
                    event_id     BIGSERIAL PRIMARY KEY,
                    patient_id   TEXT NOT NULL,
                    event_type   TEXT NOT NULL,
                    payload      JSONB NOT NULL,
                    provenance   JSONB NOT NULL,
                    created_at   TIMESTAMPTZ NOT NULL DEFAULT now(),
                    parent_event BIGINT REFERENCES engine_events(event_id)
                );
                CREATE INDEX IF NOT EXISTS ix_events_patient ON engine_events(patient_id, event_id);
                CREATE INDEX IF NOT EXISTS ix_events_type ON engine_events(event_type, created_at);
                """
            )

    def append(self, event: Event) -> int:
        with self._conn, self._conn.cursor() as cur:
            cur.execute(
                "INSERT INTO engine_events (patient_id, event_type, payload, provenance, parent_event)"
                " VALUES (%s, %s, %s, %s, %s) RETURNING event_id",
                (
                    event.patient_id,
                    event.event_type.value,
                    json.dumps(event.payload, default=str),
                    json.dumps(event.provenance, default=str),
                    event.parent_event,
                ),
            )
            event.event_id = int(cur.fetchone()[0])
        return event.event_id

    def events_for(self, patient_id: str) -> list[Event]:
        with self._conn.cursor() as cur:
            cur.execute(
                "SELECT event_id, patient_id, event_type, payload, provenance, created_at, parent_event"
                " FROM engine_events WHERE patient_id = %s ORDER BY event_id",
                (patient_id,),
            )
            rows = cur.fetchall()
        return [
            Event(
                event_id=r[0], patient_id=r[1], event_type=EventType(r[2]),
                payload=r[3], provenance=r[4], created_at=r[5], parent_event=r[6],
            )
            for r in rows
        ]


def make_store() -> EventStore:
    if settings.POSTGRES_DSN:
        return PostgresEventStore()
    return SqliteEventStore()
