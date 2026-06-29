"""
Precision Autoimmune Agent — Centralized Logging Configuration

Sets up loguru with console + file output, configurable via environment.

Usage:
    from config.logging import configure_logging
    configure_logging()

Env vars (all use AUTO_ prefix via settings):
    AUTO_LOG_LEVEL  — DEBUG, INFO, WARNING, ERROR (default: INFO)
    AUTO_LOG_DIR    — Log directory path (default: <project_root>/logs)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from loguru import logger


def configure_logging(
    log_level: str | None = None,
    log_dir: str | Path | None = None,
    service_name: str = "autoimmune-agent",
) -> None:
    """Configure loguru with console and file sinks.

    Args:
        log_level: Override log level (default: AUTO_LOG_LEVEL env or INFO).
        log_dir: Override log directory (default: AUTO_LOG_DIR env or <project>/logs).
        service_name: Used in log file name.
    """
    # Resolve settings
    level = (log_level or os.environ.get("AUTO_LOG_LEVEL", "INFO")).upper()
    project_root = Path(__file__).resolve().parent.parent
    directory = Path(log_dir or os.environ.get("AUTO_LOG_DIR", str(project_root / "logs")))
    directory.mkdir(parents=True, exist_ok=True)

    # Structured format
    fmt_console = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
        "<level>{level: <8}</level> | "
        "<cyan>{module}</cyan> | "
        "{message}"
    )
    fmt_file = (
        "{time:YYYY-MM-DD HH:mm:ss.SSS} | "
        "{level: <8} | "
        "{module} | "
        "{message}"
    )

    # Remove default loguru handler
    logger.remove()

    # Console sink with color
    logger.add(
        sys.stderr,
        format=fmt_console,
        level=level,
        colorize=True,
    )

    # File sink with rotation
    log_file = directory / f"{service_name}.log"
    logger.add(
        str(log_file),
        format=fmt_file,
        level=level,
        rotation="10 MB",
        retention=5,
        encoding="utf-8",
        enqueue=True,  # thread-safe
    )

    logger.info(
        f"Logging configured: level={level}, file={log_file}"
    )
