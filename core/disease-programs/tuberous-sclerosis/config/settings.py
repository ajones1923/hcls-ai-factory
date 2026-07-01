"""
TSC Intelligence Engine — central configuration (PRD §2.4, §2.7).

Pydantic-settings; every value overridable via environment with the TSC_ prefix.
Defaults are chosen so the engine RUNS on the DGX Spark with no external services
(SQLite event store, stub model router) and upgrades to the full Compose stack
(PostgreSQL/Redis/Milvus/Ollama) purely via environment configuration.
"""
from __future__ import annotations

from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

PROJECT_ROOT = Path(__file__).resolve().parents[1]


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="TSC_", case_sensitive=False, env_file=".env", extra="ignore"
    )

    # ── Paths ───────────────────────────────────────────────────────────
    PROJECT_ROOT: Path = PROJECT_ROOT
    DATA_DIR: Path = PROJECT_ROOT / "data"
    COHORT_DIR: Path = PROJECT_ROOT / "data" / "cohort"
    STATE_DIR: Path = PROJECT_ROOT / "data" / "state"
    CONFIG_DIR: Path = PROJECT_ROOT / "config"

    # ── Event store (PRD §2.5) ──────────────────────────────────────────
    # If POSTGRES_DSN is set, the Postgres backend is used; otherwise SQLite.
    POSTGRES_DSN: str | None = None
    SQLITE_PATH: Path = PROJECT_ROOT / "data" / "state" / "engine.db"

    # ── Service ports (PRD §2.4 — 856x band) ────────────────────────────
    API_PORT: int = 8560
    BRIEFING_PORT: int = 8561
    DASHBOARD_PORT: int = 8562
    ALERTS_PORT: int = 8563
    REDIS_URL: str = "redis://localhost:6379/0"
    MILVUS_URI: str = "http://localhost:19530"
    OLLAMA_URL: str = "http://localhost:11434"

    # ── Models (PRD §2.2; policy in config/model_policy.yaml) ───────────
    ANTHROPIC_API_KEY: str | None = None
    MODEL_POLICY_PATH: Path = PROJECT_ROOT / "config" / "model_policy.yaml"
    DEMO_CONFIG_PATH: Path = PROJECT_ROOT / "config" / "demo.yaml"
    TRAJECTORY_CONFIG_PATH: Path = PROJECT_ROOT / "config" / "trajectory_config.yaml"
    TAND_SCORING_PATH: Path = PROJECT_ROOT / "config" / "tand_scoring.yaml"
    # When no API key is present (or TSC_OFFLINE=1), the model router returns
    # deterministic, watermarked stub responses so the engine runs end-to-end.
    OFFLINE: bool = False
    # P3-3 — dispatch enrollment through the compiled LangGraph StateGraph (PRD's named
    # production runtime). Byte-equivalent to the plain dispatcher (verified); off by
    # default so deployments without langgraph still run. Set TSC_USE_LANGGRAPH=1 to enable.
    USE_LANGGRAPH: bool = False

    # ── Access control (P1-7) ───────────────────────────────────────────
    # When TSC_API_KEY is set, the API requires an X-API-Key header on data routes
    # (fail-closed). Unset (default) = open, preserving the trusted-LAN synthetic-demo
    # posture. MUST be set before any non-LAN exposure or real PHI (see SEC-17/P3-5).
    API_KEY: str | None = None

    # ── Provenance / posture ────────────────────────────────────────────
    WATERMARK: str = "SYNTHETIC — TSC Intelligence Engine demonstration data"
    AGENT_VERSION: str = "0.1.0"

    # ── Alert discipline (PRD §2.4) ─────────────────────────────────────
    ALERTS_PER_CLINICIAN_PER_WEEK_MAX: int = 3

    # ── Cost governance (PRD §5.1.3 NFR-COST) ───────────────────────────
    MAX_RUN_USD: float = 0.0          # per-process billed-spend cap; 0 = no cap

    def ensure_dirs(self) -> None:
        for p in (self.DATA_DIR, self.COHORT_DIR, self.STATE_DIR):
            p.mkdir(parents=True, exist_ok=True)


settings = Settings()
