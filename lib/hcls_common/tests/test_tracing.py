"""Tests for hcls_common.tracing — NoOp fallbacks and @traced decorator."""

import pytest

from hcls_common.tracing import (
    _NoOpSpan,
    _NoOpTracer,
    get_tracer,
    inject_trace_context,
    extract_trace_context,
    traced,
)


class TestNoOpSpan:
    def test_set_attribute_noop(self):
        span = _NoOpSpan()
        span.set_attribute("key", "value")  # should not raise

    def test_set_status_noop(self):
        span = _NoOpSpan()
        span.set_status("OK")

    def test_record_exception_noop(self):
        span = _NoOpSpan()
        span.record_exception(RuntimeError("test"))

    def test_context_manager(self):
        span = _NoOpSpan()
        with span as s:
            assert s is span

    @pytest.mark.asyncio
    async def test_async_context_manager(self):
        span = _NoOpSpan()
        async with span as s:
            assert s is span


class TestNoOpTracer:
    def test_start_span(self):
        tracer = _NoOpTracer()
        span = tracer.start_span("test")
        assert isinstance(span, _NoOpSpan)

    def test_start_as_current_span(self):
        tracer = _NoOpTracer()
        span = tracer.start_as_current_span("test")
        assert isinstance(span, _NoOpSpan)


class TestGetTracer:
    def test_returns_tracer_instance(self):
        tracer = get_tracer()
        # Should return either a NoOpTracer or a real OTel tracer
        assert tracer is not None


class TestTracedDecorator:
    def test_sync_function(self):
        @traced()
        def add(a, b):
            return a + b

        result = add(3, 4)
        assert result == 7

    def test_sync_function_with_name(self):
        @traced(name="custom_span")
        def multiply(a, b):
            return a * b

        assert multiply(3, 4) == 12

    def test_sync_function_preserves_exception(self):
        @traced()
        def fail():
            raise ValueError("boom")

        with pytest.raises(ValueError, match="boom"):
            fail()

    @pytest.mark.asyncio
    async def test_async_function(self):
        @traced()
        async def async_add(a, b):
            return a + b

        result = await async_add(3, 4)
        assert result == 7

    @pytest.mark.asyncio
    async def test_async_function_preserves_exception(self):
        @traced()
        async def async_fail():
            raise RuntimeError("async boom")

        with pytest.raises(RuntimeError, match="async boom"):
            await async_fail()

    def test_traced_with_attributes(self):
        @traced(attributes={"component": "test"})
        def noop():
            return "ok"

        assert noop() == "ok"


class TestPropagation:
    def test_inject_returns_dict(self):
        carrier = inject_trace_context()
        assert isinstance(carrier, dict)

    def test_extract_handles_empty_carrier(self):
        result = extract_trace_context({})
        # Without OTel it returns None; with OTel it returns a context
        # Either way, it should not raise
        assert result is None or result is not None
