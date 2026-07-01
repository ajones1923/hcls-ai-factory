"""Tests for hcls_common.circuit_breaker — CircuitBreaker state machine."""

import time
import pytest

from hcls_common.circuit_breaker import (
    CircuitBreaker,
    CircuitBreakerOpen,
    CircuitState,
    get_breaker,
    all_breakers,
)


class TestCircuitBreakerInit:
    def test_default_state_is_closed(self):
        cb = CircuitBreaker("test-init", failure_threshold=3, reset_timeout=1.0)
        assert cb.state == CircuitState.CLOSED

    def test_registered_in_global_registry(self):
        cb = CircuitBreaker("test-registry", failure_threshold=3)
        assert get_breaker("test-registry") is cb

    def test_all_breakers_returns_snapshot(self):
        cb = CircuitBreaker("test-all-breakers")
        breakers = all_breakers()
        assert "test-all-breakers" in breakers

    def test_initial_failure_count_is_zero(self):
        cb = CircuitBreaker("test-fc-init")
        assert cb.failure_count == 0


class TestCircuitBreakerTransitions:
    def test_closed_to_open_after_threshold(self):
        cb = CircuitBreaker("test-c2o", failure_threshold=3, reset_timeout=60)
        for _ in range(3):
            cb._record_failure()
        assert cb.state == CircuitState.OPEN

    def test_stays_closed_below_threshold(self):
        cb = CircuitBreaker("test-below", failure_threshold=5, reset_timeout=60)
        for _ in range(4):
            cb._record_failure()
        assert cb.state == CircuitState.CLOSED

    def test_success_resets_failure_count(self):
        cb = CircuitBreaker("test-reset-fc", failure_threshold=5, reset_timeout=60)
        cb._record_failure()
        cb._record_failure()
        cb._record_success()
        assert cb.failure_count == 0

    def test_open_to_half_open_after_timeout(self):
        cb = CircuitBreaker("test-o2ho", failure_threshold=1, reset_timeout=0.1)
        cb._record_failure()
        assert cb.state == CircuitState.OPEN
        time.sleep(0.15)
        assert cb.state == CircuitState.HALF_OPEN

    def test_half_open_to_closed_after_successes(self):
        cb = CircuitBreaker(
            "test-ho2c", failure_threshold=1, reset_timeout=0.1, success_threshold=2
        )
        cb._record_failure()
        time.sleep(0.15)
        assert cb.state == CircuitState.HALF_OPEN
        cb._record_success()
        cb._record_success()
        assert cb.state == CircuitState.CLOSED

    def test_half_open_to_open_on_failure(self):
        cb = CircuitBreaker("test-ho2o", failure_threshold=1, reset_timeout=0.1)
        cb._record_failure()
        time.sleep(0.15)
        assert cb.state == CircuitState.HALF_OPEN
        cb._record_failure()
        assert cb.state == CircuitState.OPEN


class TestCircuitBreakerRejection:
    def test_open_breaker_raises(self):
        cb = CircuitBreaker("test-reject", failure_threshold=1, reset_timeout=60)
        cb._record_failure()
        with pytest.raises(CircuitBreakerOpen) as exc_info:
            cb._check_state()
        assert exc_info.value.name == "test-reject"
        assert exc_info.value.remaining_seconds > 0

    def test_remaining_seconds_when_closed(self):
        cb = CircuitBreaker("test-remain-closed")
        assert cb.remaining_open_seconds == 0.0


class TestCircuitBreakerDecorator:
    def test_decorator_passes_on_success(self):
        cb = CircuitBreaker("test-dec-ok", failure_threshold=3)

        @cb
        def add(a, b):
            return a + b

        assert add(1, 2) == 3

    def test_decorator_records_failure(self):
        cb = CircuitBreaker("test-dec-fail", failure_threshold=3)

        @cb
        def fail():
            raise RuntimeError("boom")

        with pytest.raises(RuntimeError):
            fail()
        assert cb.failure_count == 1

    def test_decorator_rejects_when_open(self):
        cb = CircuitBreaker("test-dec-open", failure_threshold=1, reset_timeout=60)

        @cb
        def func():
            return "ok"

        # Trip the breaker
        cb._record_failure()
        with pytest.raises(CircuitBreakerOpen):
            func()

    def test_excluded_exceptions_not_counted(self):
        cb = CircuitBreaker(
            "test-excluded",
            failure_threshold=3,
            excluded_exceptions=(ValueError,),
        )

        @cb
        def validate(x):
            if x < 0:
                raise ValueError("negative")
            return x

        with pytest.raises(ValueError):
            validate(-1)
        assert cb.failure_count == 0


class TestCircuitBreakerContextManager:
    def test_sync_context_manager_success(self):
        cb = CircuitBreaker("test-ctx-ok", failure_threshold=3)
        with cb:
            result = 42
        assert result == 42
        assert cb.failure_count == 0

    def test_sync_context_manager_failure(self):
        cb = CircuitBreaker("test-ctx-fail", failure_threshold=3)
        with pytest.raises(RuntimeError):
            with cb:
                raise RuntimeError("boom")
        assert cb.failure_count == 1


class TestCircuitBreakerManualControls:
    def test_manual_reset(self):
        cb = CircuitBreaker("test-manual-reset", failure_threshold=1, reset_timeout=60)
        cb._record_failure()
        assert cb.state == CircuitState.OPEN
        cb.reset()
        assert cb.state == CircuitState.CLOSED
        assert cb.failure_count == 0

    def test_manual_trip(self):
        cb = CircuitBreaker("test-manual-trip")
        cb.trip()
        assert cb.state == CircuitState.OPEN

    def test_repr(self):
        cb = CircuitBreaker("test-repr", failure_threshold=5)
        r = repr(cb)
        assert "test-repr" in r
        assert "closed" in r


class TestCircuitBreakerAsync:
    @pytest.mark.asyncio
    async def test_async_context_manager_success(self):
        cb = CircuitBreaker("test-async-ok", failure_threshold=3)
        async with cb:
            result = 42
        assert result == 42

    @pytest.mark.asyncio
    async def test_async_context_manager_failure(self):
        cb = CircuitBreaker("test-async-fail", failure_threshold=3)
        with pytest.raises(RuntimeError):
            async with cb:
                raise RuntimeError("async boom")
        assert cb.failure_count == 1

    @pytest.mark.asyncio
    async def test_async_decorator(self):
        cb = CircuitBreaker("test-async-dec", failure_threshold=3)

        @cb
        async def async_add(a, b):
            return a + b

        result = await async_add(3, 4)
        assert result == 7
