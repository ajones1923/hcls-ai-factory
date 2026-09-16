"""Tests for idempotent collector registration.

Every collector in this library is created at MODULE level, so a duplicate registration raises
during import and takes the whole module — and every caller — down with it. A cardiology test
reached exactly that: `mock.patch("src.rag_engine.X")` imported a module under a second identity
and re-executed its imports.
"""
import pytest

prometheus_client = pytest.importorskip("prometheus_client")
from prometheus_client import CollectorRegistry, Counter, Gauge, Histogram  # noqa: E402

from hcls_common.metrics import metric  # noqa: E402


@pytest.fixture
def reg():
    """An isolated registry, so these tests never touch the process-global one."""
    return CollectorRegistry()


class TestFirstRegistration:
    def test_returns_a_working_collector(self, reg):
        c = metric(Counter, "t_calls_total", "calls", ["provider"], registry=reg)
        c.labels(provider="anthropic").inc()
        assert reg.get_sample_value("t_calls_total", {"provider": "anthropic"}) == 1.0

    def test_labels_are_optional(self, reg):
        g = metric(Gauge, "t_pool", "pool", registry=reg)
        g.set(3)
        assert reg.get_sample_value("t_pool") == 3.0

    def test_kwargs_pass_through(self, reg):
        metric(Histogram, "t_latency_seconds", "latency", ["c"], buckets=(0.1, 1.0), registry=reg)
        assert reg.get_sample_value(
            "t_latency_seconds_bucket", {"c": "x", "le": "0.1"}) is None  # not yet observed


class TestDuplicateRegistration:
    def test_second_registration_reuses_the_first(self):
        """The real case: process-global registry, module imported twice."""
        a = metric(Counter, "t_dupe_total", "dupe")
        b = metric(Counter, "t_dupe_total", "dupe")
        assert a is b, "a re-import must reuse the collector, not raise"

    def test_reused_collector_still_records(self):
        a = metric(Counter, "t_dupe2_total", "dupe2")
        a.inc()
        b = metric(Counter, "t_dupe2_total", "dupe2")
        b.inc()
        assert prometheus_client.REGISTRY.get_sample_value("t_dupe2_total") == 2.0

    def test_module_level_import_twice_does_not_raise(self):
        import importlib
        import sys
        for name in ("hcls_common.milvus_client", "hcls_common.embedder",
                     "hcls_common.llm_client", "hcls_common.circuit_breaker",
                     "hcls_common.event_bus", "hcls_common.query_router",
                     "hcls_common.bidirectional_triggers"):
            importlib.import_module(name)
            del sys.modules[name]
            importlib.import_module(name)          # the second identity must not explode


class TestGenuineErrorsStillRaise:
    def test_a_real_value_error_is_not_swallowed(self, reg):
        """Only a name collision may be absorbed — a malformed metric must still fail.

        `le` is reserved on a histogram, so this ValueError is not a re-registration and must
        reach the caller. Absorbing every ValueError would hide genuine misuse.
        """
        with pytest.raises(ValueError):
            metric(Histogram, "t_reserved_seconds", "reserved label", ["le"], registry=reg)

    def test_a_reserved_label_on_a_counter_also_raises(self, reg):
        with pytest.raises(ValueError):
            metric(Counter, "t_reserved_total", "reserved label", ["__private"], registry=reg)
