"""
OpenTelemetry integration for the HCLS AI Factory.

Provides:
  - ``init_tracing(service_name)`` -- bootstrap the OTel SDK (gRPC exporter)
  - ``@traced`` decorator            -- auto-span for sync and async functions
  - ``inject_trace_context`` / ``extract_trace_context`` -- W3C traceparent
    propagation for inter-service calls (Nextflow, HTTP, etc.)

Graceful degradation: if ``opentelemetry-api`` is not installed every
function becomes a no-op so downstream code never needs to guard imports.
"""

from __future__ import annotations

import asyncio
import functools
import logging
from typing import Any, Callable, Dict, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Detect OpenTelemetry availability
# ---------------------------------------------------------------------------

try:
    from opentelemetry import trace
    from opentelemetry.context import attach, detach
    from opentelemetry.trace import (
        StatusCode,
        Tracer,
        set_span_in_context,
    )
    from opentelemetry.trace.propagation import get_current_span

    _HAS_OTEL = True
except ImportError:
    _HAS_OTEL = False

# ---------------------------------------------------------------------------
# No-op fallbacks
# ---------------------------------------------------------------------------


class _NoOpSpan:
    """Minimal span stand-in when OTel is absent."""

    def set_attribute(self, key: str, value: Any) -> None:
        pass

    def set_status(self, *args: Any, **kwargs: Any) -> None:
        pass

    def record_exception(self, exc: BaseException) -> None:
        pass

    def end(self) -> None:
        pass

    def __enter__(self) -> "_NoOpSpan":
        return self

    def __exit__(self, *exc: Any) -> None:
        pass

    async def __aenter__(self) -> "_NoOpSpan":
        return self

    async def __aexit__(self, *exc: Any) -> None:
        pass


class _NoOpTracer:
    """Minimal tracer stand-in when OTel is absent."""

    def start_span(self, name: str, **kwargs: Any) -> _NoOpSpan:
        return _NoOpSpan()

    def start_as_current_span(self, name: str, **kwargs: Any) -> _NoOpSpan:
        return _NoOpSpan()


_noop_tracer = _NoOpTracer()

# Module-level tracer -- set by init_tracing()
_tracer: Any = _noop_tracer

# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------


def init_tracing(
    service_name: str = "hcls-ai-factory",
    otlp_endpoint: Optional[str] = None,
) -> Any:
    """
    Bootstrap the OpenTelemetry SDK.

    Parameters
    ----------
    service_name : str
        Logical service name used in span metadata.
    otlp_endpoint : str or None
        gRPC endpoint for the OTel collector (e.g. ``http://localhost:4317``).
        If ``None``, reads ``OTEL_EXPORTER_OTLP_ENDPOINT`` from the env,
        and falls back to a console exporter when no endpoint is set.

    Returns
    -------
    opentelemetry.trace.Tracer or _NoOpTracer
        A tracer instance.  If OTel is not installed, returns ``_NoOpTracer``.
    """
    global _tracer

    if not _HAS_OTEL:
        logger.info(
            "OpenTelemetry SDK not installed -- tracing disabled. "
            "Install with: pip install opentelemetry-api opentelemetry-sdk "
            "opentelemetry-exporter-otlp-proto-grpc"
        )
        return _noop_tracer

    import os

    from opentelemetry import trace as _trace
    from opentelemetry.sdk.resources import Resource
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import (
        BatchSpanProcessor,
        ConsoleSpanExporter,
    )

    resource = Resource.create({"service.name": service_name})
    provider = TracerProvider(resource=resource)

    endpoint = otlp_endpoint or os.environ.get("OTEL_EXPORTER_OTLP_ENDPOINT")

    if endpoint:
        try:
            from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import (
                OTLPSpanExporter,
            )

            exporter = OTLPSpanExporter(endpoint=endpoint)
            provider.add_span_processor(BatchSpanProcessor(exporter))
            logger.info(
                "OTel tracing enabled: service=%s  endpoint=%s",
                service_name,
                endpoint,
            )
        except ImportError:
            logger.warning(
                "opentelemetry-exporter-otlp-proto-grpc not installed; "
                "falling back to console exporter"
            )
            provider.add_span_processor(
                BatchSpanProcessor(ConsoleSpanExporter())
            )
    else:
        # Console exporter for local dev
        provider.add_span_processor(
            BatchSpanProcessor(ConsoleSpanExporter())
        )
        logger.info(
            "OTel tracing enabled (console exporter): service=%s", service_name
        )

    _trace.set_tracer_provider(provider)
    _tracer = _trace.get_tracer(service_name)
    return _tracer


def get_tracer() -> Any:
    """Return the current module-level tracer (may be a no-op)."""
    return _tracer


# ---------------------------------------------------------------------------
# @traced decorator
# ---------------------------------------------------------------------------


def traced(
    name: Optional[str] = None,
    attributes: Optional[Dict[str, Any]] = None,
) -> Callable:
    """
    Decorator that wraps a function in an OTel span.

    Works with both sync and async functions.

    Parameters
    ----------
    name : str or None
        Span name.  Defaults to ``module.function``.
    attributes : dict or None
        Static span attributes set at span start.

    Example::

        @traced()
        def search_variants(gene: str):
            ...

        @traced(name="embed_batch", attributes={"component": "embedder"})
        async def embed_batch(texts):
            ...
    """

    def decorator(func: Callable) -> Callable:
        span_name = name or f"{func.__module__}.{func.__qualname__}"

        if asyncio.iscoroutinefunction(func):

            @functools.wraps(func)
            async def _async_wrapper(*args: Any, **kwargs: Any) -> Any:
                tracer = _tracer
                if isinstance(tracer, _NoOpTracer):
                    return await func(*args, **kwargs)

                with tracer.start_as_current_span(span_name) as span:
                    if attributes:
                        for k, v in attributes.items():
                            span.set_attribute(k, v)
                    try:
                        result = await func(*args, **kwargs)
                        return result
                    except Exception as exc:
                        span.set_status(
                            StatusCode.ERROR, str(exc)
                        )
                        span.record_exception(exc)
                        raise

            return _async_wrapper
        else:

            @functools.wraps(func)
            def _sync_wrapper(*args: Any, **kwargs: Any) -> Any:
                tracer = _tracer
                if isinstance(tracer, _NoOpTracer):
                    return func(*args, **kwargs)

                with tracer.start_as_current_span(span_name) as span:
                    if attributes:
                        for k, v in attributes.items():
                            span.set_attribute(k, v)
                    try:
                        result = func(*args, **kwargs)
                        return result
                    except Exception as exc:
                        span.set_status(
                            StatusCode.ERROR, str(exc)
                        )
                        span.record_exception(exc)
                        raise

            return _sync_wrapper

    return decorator


# ---------------------------------------------------------------------------
# W3C traceparent propagation helpers
# ---------------------------------------------------------------------------


def inject_trace_context() -> Dict[str, str]:
    """
    Inject the current trace context into a carrier dict (W3C traceparent
    format).  Pass the returned dict as HTTP headers or Nextflow metadata.

    Returns an empty dict if OTel is not available.
    """
    if not _HAS_OTEL:
        return {}

    from opentelemetry.propagate import inject

    carrier: Dict[str, str] = {}
    inject(carrier)
    return carrier


def extract_trace_context(carrier: Dict[str, str]) -> Any:
    """
    Extract trace context from a carrier dict and return an OTel ``Context``
    that can be used with ``attach()``.

    Returns ``None`` if OTel is not available.
    """
    if not _HAS_OTEL:
        return None

    from opentelemetry.propagate import extract

    return extract(carrier)
