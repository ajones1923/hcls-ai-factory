"""
Single-box MLOps (A3) — experiment tracking + model registry + run lineage, with no
data warehouse and no external service. A SQLite-backed store that closes the platform's
biggest structural gap: every model run is recorded (params, metrics, artifacts, and the
reproducibility manifest), models are versioned and staged, and runs are queryable (by the
composer, the MCP surface, or a human).

Stdlib only. Use ``:memory:`` for tests, or a file path (e.g. data/mlops.db) for a durable
store. Timestamps are real wall-clock; ids are uuid4.
"""
from __future__ import annotations

import json
import sqlite3
import time
import uuid
from pathlib import Path
from typing import Any

_SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id     TEXT PRIMARY KEY,
    name       TEXT,
    capability TEXT,
    status     TEXT,
    started_at REAL,
    ended_at   REAL,
    manifest   TEXT
);
CREATE TABLE IF NOT EXISTS run_params  (run_id TEXT, key TEXT, value TEXT);
CREATE TABLE IF NOT EXISTS run_metrics (run_id TEXT, key TEXT, value REAL, step INTEGER);
CREATE TABLE IF NOT EXISTS model_versions (
    name       TEXT,
    version    TEXT,
    stage      TEXT,         -- none | staging | production | archived
    source     TEXT,
    run_id     TEXT,
    metrics    TEXT,
    created_at REAL,
    PRIMARY KEY (name, version)
);
"""

_STAGES = {"none", "staging", "production", "archived"}
# A6: mandatory batch-job status state-machine
RUN_STATES = ("submitted", "started", "running", "complete", "failed")


class MLOpsStore:
    def __init__(self, db_path: str | Path = ":memory:") -> None:
        if db_path != ":memory:":
            Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(db_path))
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(_SCHEMA)
        self.conn.commit()

    # -- experiment tracking ------------------------------------------------ #
    def start_run(self, name: str, capability: str | None = None,
                  params: dict[str, Any] | None = None,
                  manifest: dict[str, Any] | None = None) -> str:
        run_id = uuid.uuid4().hex
        self.conn.execute(
            "INSERT INTO runs (run_id,name,capability,status,started_at,ended_at,manifest) "
            "VALUES (?,?,?,?,?,?,?)",
            (run_id, name, capability, "running", time.time(), None,
             json.dumps(manifest) if manifest else None))
        for k, v in (params or {}).items():
            self.conn.execute("INSERT INTO run_params VALUES (?,?,?)", (run_id, k, json.dumps(v)))
        self.conn.commit()
        return run_id

    def log_metric(self, run_id: str, key: str, value: float, step: int = 0) -> None:
        self.conn.execute("INSERT INTO run_metrics VALUES (?,?,?,?)", (run_id, key, float(value), step))
        self.conn.commit()

    def log_param(self, run_id: str, key: str, value: Any) -> None:
        self.conn.execute("INSERT INTO run_params VALUES (?,?,?)", (run_id, key, json.dumps(value)))
        self.conn.commit()

    def end_run(self, run_id: str, status: str = "completed") -> None:
        self.conn.execute("UPDATE runs SET status=?, ended_at=? WHERE run_id=?",
                          (status, time.time(), run_id))
        self.conn.commit()

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        row = self.conn.execute("SELECT * FROM runs WHERE run_id=?", (run_id,)).fetchone()
        if not row:
            return None
        d = dict(row)
        d["manifest"] = json.loads(d["manifest"]) if d["manifest"] else None
        d["params"] = {r["key"]: json.loads(r["value"]) for r in
                       self.conn.execute("SELECT key,value FROM run_params WHERE run_id=?", (run_id,))}
        d["metrics"] = {r["key"]: r["value"] for r in
                        self.conn.execute("SELECT key,value FROM run_metrics WHERE run_id=? ORDER BY step", (run_id,))}
        return d

    def list_runs(self, capability: str | None = None, limit: int = 50) -> list[dict[str, Any]]:
        if capability:
            rows = self.conn.execute(
                "SELECT run_id,name,capability,status,started_at FROM runs WHERE capability=? "
                "ORDER BY started_at DESC LIMIT ?", (capability, limit))
        else:
            rows = self.conn.execute(
                "SELECT run_id,name,capability,status,started_at FROM runs "
                "ORDER BY started_at DESC LIMIT ?", (limit,))
        return [dict(r) for r in rows]

    # -- A6: status state-machine + run search -------------------------------- #
    def set_status(self, run_id: str, status: str) -> None:
        """Advance a run through submitted→started→running→complete/failed."""
        if status not in RUN_STATES:
            raise ValueError(f"invalid status {status!r}; one of {RUN_STATES}")
        ended = time.time() if status in ("complete", "failed") else None
        self.conn.execute("UPDATE runs SET status=?, ended_at=COALESCE(?, ended_at) WHERE run_id=?",
                          (status, ended, run_id))
        self.conn.commit()

    def search_runs(self, name_like: str | None = None, capability: str | None = None,
                    status: str | None = None, limit: int = 50) -> list[dict[str, Any]]:
        """'Search past runs' surface — filter by name substring / capability / status."""
        clauses, args = [], []
        if name_like:
            clauses.append("name LIKE ?"); args.append(f"%{name_like}%")
        if capability:
            clauses.append("capability = ?"); args.append(capability)
        if status:
            clauses.append("status = ?"); args.append(status)
        where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
        args.append(limit)
        rows = self.conn.execute(
            f"SELECT run_id,name,capability,status,started_at,ended_at FROM runs "
            f"{where} ORDER BY started_at DESC LIMIT ?", args)
        return [dict(r) for r in rows]

    def best_run(self, metric: str, capability: str | None = None, maximize: bool = True) -> dict[str, Any] | None:
        order = "DESC" if maximize else "ASC"
        q = ("SELECT m.run_id, m.value FROM run_metrics m JOIN runs r ON r.run_id=m.run_id "
             "WHERE m.key=? " + ("AND r.capability=? " if capability else "") +
             f"ORDER BY m.value {order} LIMIT 1")
        args = (metric, capability) if capability else (metric,)
        row = self.conn.execute(q, args).fetchone()
        return self.get_run(row["run_id"]) if row else None

    # -- model registry ----------------------------------------------------- #
    def register_model(self, name: str, version: str, source: str | None = None,
                       run_id: str | None = None, metrics: dict[str, Any] | None = None,
                       stage: str = "none") -> None:
        if stage not in _STAGES:
            raise ValueError(f"invalid stage {stage!r}; one of {sorted(_STAGES)}")
        self.conn.execute(
            "INSERT OR REPLACE INTO model_versions VALUES (?,?,?,?,?,?,?)",
            (name, version, stage, source, run_id, json.dumps(metrics or {}), time.time()))
        self.conn.commit()

    def transition_stage(self, name: str, version: str, stage: str) -> None:
        if stage not in _STAGES:
            raise ValueError(f"invalid stage {stage!r}")
        # only one production version at a time -> archive the incumbent
        if stage == "production":
            self.conn.execute("UPDATE model_versions SET stage='archived' "
                              "WHERE name=? AND stage='production'", (name,))
        self.conn.execute("UPDATE model_versions SET stage=? WHERE name=? AND version=?",
                          (stage, name, version))
        self.conn.commit()

    def get_model_version(self, name: str, version: str) -> dict[str, Any] | None:
        row = self.conn.execute("SELECT * FROM model_versions WHERE name=? AND version=?",
                               (name, version)).fetchone()
        if not row:
            return None
        d = dict(row); d["metrics"] = json.loads(d["metrics"]); return d

    def list_model_versions(self, name: str) -> list[dict[str, Any]]:
        rows = self.conn.execute("SELECT * FROM model_versions WHERE name=? ORDER BY created_at DESC", (name,))
        out = []
        for r in rows:
            d = dict(r); d["metrics"] = json.loads(d["metrics"]); out.append(d)
        return out

    def get_production(self, name: str) -> dict[str, Any] | None:
        row = self.conn.execute("SELECT * FROM model_versions WHERE name=? AND stage='production'",
                               (name,)).fetchone()
        if not row:
            return None
        d = dict(row); d["metrics"] = json.loads(d["metrics"]); return d


# process-wide durable store (under the lib's data dir by default)
_DEFAULT_DB = Path(__file__).resolve().parents[2] / "data" / "mlops.db"
_STORE: MLOpsStore | None = None


def get_mlops_store(db_path: str | Path | None = None) -> MLOpsStore:
    global _STORE
    if _STORE is None or db_path is not None:
        _STORE = MLOpsStore(db_path or _DEFAULT_DB)
    return _STORE
