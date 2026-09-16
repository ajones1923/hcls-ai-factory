"""Idempotent Prometheus collector registration.

Registering the same metric name twice raises `ValueError: Duplicated timeseries in
CollectorRegistry`. Because every collector in this library is created at MODULE level, that
exception kills the import of the module — and with it every caller.

That is reachable without anyone doing anything unusual. `mock.patch("src.rag_engine.Thing")`
resolves its target by importing `src.rag_engine`; if the test file already imported a bare
`rag_engine`, Python now holds two distinct module objects for the same file and executes its
imports a second time. A cardiology test did exactly that, and the moment a `hcls_common` import
was added to that engine the whole module failed to load.

A library module must not explode because it was imported twice. A duplicate registration reuses
the collector already registered under that name, which is the only sane outcome: the metric is
process-global anyway.
"""
from __future__ import annotations

from typing import Any, Sequence


def metric(cls: Any, name: str, doc: str, labels: Sequence[str] | None = None, **kwargs: Any):
    """Register a Prometheus collector, or return the one already registered as `name`.

    Args:
        cls: the collector class (`Counter`, `Histogram`, `Gauge`, …).
        name: metric name — the registry key.
        doc: help text.
        labels: label names, if the metric has any.
        **kwargs: passed through (e.g. `buckets=`).
    """
    args = (name, doc) + ((list(labels),) if labels else ())
    try:
        return cls(*args, **kwargs)
    except ValueError:
        try:
            from prometheus_client import REGISTRY
        except ImportError:
            raise
        existing = getattr(REGISTRY, "_names_to_collectors", {}).get(name)
        if existing is None:
            raise                       # a genuine ValueError, not a re-registration
        return existing
