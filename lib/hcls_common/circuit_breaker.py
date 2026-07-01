"""
General-purpose circuit breaker for the HCLS AI Factory.

Protects against cascading failures when calling external services (Milvus,
NIM microservices, LLM APIs).  Works as a decorator, sync context manager,
and async context manager.  Thread-safe.

States:
    CLOSED    -- normal operation; failures are counted
    OPEN      -- requests are rejected immediately for ``reset_timeout`` seconds
    HALF_OPEN -- a limited number of probe requests are allowed through

Example usage::

    cb = CircuitBreaker("milvus", failure_threshold=5, reset_timeout=30)

    # As decorator
    @cb
    def search_milvus(query):
        ...

    # As context manager
    with cb:
        client.search(...)

    # Async context manager
    async with cb:
        await client.search(...)
"""

from __future__ import annotations

import asyncio
import functools
import logging
import threading
import time
from enum import Enum, unique
from typing import Any, Callable, Dict, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Prometheus metrics (optional)
# ---------------------------------------------------------------------------

try:
    from prometheus_client import Counter, Gauge

    _cb_state_gauge = Gauge(
        "hcls_circuit_breaker_state",
        "Current circuit breaker state (0=closed, 1=open, 2=half_open)",
        ["name"],
    )
    _cb_rejected_total = Counter(
        "hcls_circuit_breaker_rejected_total",
        "Total requests rejected by open circuit breaker",
        ["name"],
    )
    _cb_failure_total = Counter(
        "hcls_circuit_breaker_failure_total",
        "Total failures recorded by circuit breaker",
        ["name"],
    )
    _cb_success_total = Counter(
        "hcls_circuit_breaker_success_total",
        "Total successes recorded by circuit breaker",
        ["name"],
    )
    _HAS_PROMETHEUS = True
except ImportError:
    _HAS_PROMETHEUS = False


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

@unique
class CircuitState(str, Enum):
    CLOSED = "closed"
    OPEN = "open"
    HALF_OPEN = "half_open"


class CircuitBreakerOpen(Exception):
    """Raised when the circuit breaker is open and rejects a call."""

    def __init__(self, name: str, remaining_seconds: float):
        self.name = name
        self.remaining_seconds = remaining_seconds
        super().__init__(
            f"Circuit breaker '{name}' is OPEN. "
            f"Retry after {remaining_seconds:.1f}s."
        )


# ---------------------------------------------------------------------------
# Registry -- allows global lookup of breakers by name
# ---------------------------------------------------------------------------

_registry: Dict[str, "CircuitBreaker"] = {}
_registry_lock = threading.Lock()


def get_breaker(name: str) -> Optional["CircuitBreaker"]:
    """Look up a circuit breaker by name. Returns ``None`` if not registered."""
    return _registry.get(name)


def all_breakers() -> Dict[str, "CircuitBreaker"]:
    """Return a snapshot of all registered circuit breakers."""
    return dict(_registry)


# ---------------------------------------------------------------------------
# CircuitBreaker
# ---------------------------------------------------------------------------

class CircuitBreaker:
    """
    Thread-safe circuit breaker with configurable thresholds.

    Parameters
    ----------
    name : str
        Unique identifier (used in metrics labels and the registry).
    failure_threshold : int
        Number of consecutive failures before opening the circuit (default 5).
    reset_timeout : float
        Seconds the circuit stays OPEN before moving to HALF_OPEN (default 30).
    success_threshold : int
        Consecutive successes in HALF_OPEN needed to close the circuit (default 2).
    excluded_exceptions : tuple
        Exception types that should **not** count as failures (e.g. ValueError).
    """

    def __init__(
        self,
        name: str,
        failure_threshold: int = 5,
        reset_timeout: float = 30.0,
        success_threshold: int = 2,
        excluded_exceptions: tuple = (),
    ) -> None:
        self.name = name
        self.failure_threshold = failure_threshold
        self.reset_timeout = reset_timeout
        self.success_threshold = success_threshold
        self.excluded_exceptions = excluded_exceptions

        self._lock = threading.Lock()
        self._state = CircuitState.CLOSED
        self._failure_count = 0
        self._success_count = 0
        self._last_failure_time: float = 0.0
        self._opened_at: float = 0.0

        # Register globally
        with _registry_lock:
            _registry[name] = self

        self._emit_state()

    # -- State inspection ---------------------------------------------------

    @property
    def state(self) -> CircuitState:
        with self._lock:
            self._maybe_transition_to_half_open()
            return self._state

    @property
    def failure_count(self) -> int:
        return self._failure_count

    @property
    def remaining_open_seconds(self) -> float:
        """Seconds remaining until the circuit transitions to HALF_OPEN."""
        with self._lock:
            if self._state != CircuitState.OPEN:
                return 0.0
            elapsed = time.monotonic() - self._opened_at
            return max(0.0, self.reset_timeout - elapsed)

    # -- Internal transitions -----------------------------------------------

    def _maybe_transition_to_half_open(self) -> None:
        """Must be called while holding ``_lock``."""
        if self._state == CircuitState.OPEN:
            if time.monotonic() - self._opened_at >= self.reset_timeout:
                logger.info(
                    "Circuit breaker '%s': OPEN -> HALF_OPEN", self.name
                )
                self._state = CircuitState.HALF_OPEN
                self._success_count = 0
                self._emit_state()

    def _record_success(self) -> None:
        with self._lock:
            if self._state == CircuitState.HALF_OPEN:
                self._success_count += 1
                if self._success_count >= self.success_threshold:
                    logger.info(
                        "Circuit breaker '%s': HALF_OPEN -> CLOSED", self.name
                    )
                    self._state = CircuitState.CLOSED
                    self._failure_count = 0
                    self._success_count = 0
                    self._emit_state()
            elif self._state == CircuitState.CLOSED:
                # Reset failure count on success
                self._failure_count = 0

        if _HAS_PROMETHEUS:
            _cb_success_total.labels(name=self.name).inc()

    def _record_failure(self) -> None:
        with self._lock:
            self._failure_count += 1
            self._last_failure_time = time.monotonic()

            if self._state == CircuitState.HALF_OPEN:
                # Any failure in HALF_OPEN immediately re-opens
                logger.warning(
                    "Circuit breaker '%s': HALF_OPEN -> OPEN (probe failed)",
                    self.name,
                )
                self._state = CircuitState.OPEN
                self._opened_at = time.monotonic()
                self._emit_state()
            elif (
                self._state == CircuitState.CLOSED
                and self._failure_count >= self.failure_threshold
            ):
                logger.warning(
                    "Circuit breaker '%s': CLOSED -> OPEN (%d failures)",
                    self.name,
                    self._failure_count,
                )
                self._state = CircuitState.OPEN
                self._opened_at = time.monotonic()
                self._emit_state()

        if _HAS_PROMETHEUS:
            _cb_failure_total.labels(name=self.name).inc()

    def _check_state(self) -> None:
        """Raise if the circuit is OPEN."""
        with self._lock:
            self._maybe_transition_to_half_open()
            if self._state == CircuitState.OPEN:
                remaining = max(
                    0.0,
                    self.reset_timeout
                    - (time.monotonic() - self._opened_at),
                )
                if _HAS_PROMETHEUS:
                    _cb_rejected_total.labels(name=self.name).inc()
                raise CircuitBreakerOpen(self.name, remaining)

    def _emit_state(self) -> None:
        if _HAS_PROMETHEUS:
            val = {
                CircuitState.CLOSED: 0,
                CircuitState.OPEN: 1,
                CircuitState.HALF_OPEN: 2,
            }[self._state]
            _cb_state_gauge.labels(name=self.name).set(val)

    # -- Manual controls ----------------------------------------------------

    def reset(self) -> None:
        """Manually reset the circuit breaker to CLOSED."""
        with self._lock:
            self._state = CircuitState.CLOSED
            self._failure_count = 0
            self._success_count = 0
            self._emit_state()
        logger.info("Circuit breaker '%s': manually reset to CLOSED", self.name)

    def trip(self) -> None:
        """Manually trip (open) the circuit breaker."""
        with self._lock:
            self._state = CircuitState.OPEN
            self._opened_at = time.monotonic()
            self._emit_state()
        logger.warning("Circuit breaker '%s': manually tripped to OPEN", self.name)

    # -- Decorator ----------------------------------------------------------

    def __call__(self, func: Callable) -> Callable:
        """Use as ``@circuit_breaker`` decorator."""
        if asyncio.iscoroutinefunction(func):

            @functools.wraps(func)
            async def _async_wrapper(*args: Any, **kwargs: Any) -> Any:
                self._check_state()
                try:
                    result = await func(*args, **kwargs)
                except self.excluded_exceptions:
                    raise
                except Exception:
                    self._record_failure()
                    raise
                else:
                    self._record_success()
                    return result

            return _async_wrapper
        else:

            @functools.wraps(func)
            def _sync_wrapper(*args: Any, **kwargs: Any) -> Any:
                self._check_state()
                try:
                    result = func(*args, **kwargs)
                except self.excluded_exceptions:
                    raise
                except Exception:
                    self._record_failure()
                    raise
                else:
                    self._record_success()
                    return result

            return _sync_wrapper

    # -- Sync context manager -----------------------------------------------

    def __enter__(self) -> "CircuitBreaker":
        self._check_state()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if exc_type is None:
            self._record_success()
        elif self.excluded_exceptions and issubclass(
            exc_type, self.excluded_exceptions
        ):
            pass  # excluded -- do not count
        else:
            self._record_failure()

    # -- Async context manager ----------------------------------------------

    async def __aenter__(self) -> "CircuitBreaker":
        self._check_state()
        return self

    async def __aexit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if exc_type is None:
            self._record_success()
        elif self.excluded_exceptions and issubclass(
            exc_type, self.excluded_exceptions
        ):
            pass
        else:
            self._record_failure()

    # -- repr ---------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"CircuitBreaker(name={self.name!r}, state={self._state.value}, "
            f"failures={self._failure_count}/{self.failure_threshold})"
        )
