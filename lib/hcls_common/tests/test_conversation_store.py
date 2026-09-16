"""Tests for shared session-memory persistence.

Five subjects carried this code byte-for-byte and NONE of them tested it. Session memory that
stops persisting looks exactly like an agent that has forgotten the conversation — a fault with
no outward signal.
"""
import json
from datetime import datetime, timedelta, timezone

import pytest

from hcls_common.conversation_store import ConversationStore


@pytest.fixture
def store(tmp_path):
    return ConversationStore(tmp_path / "conversations")


class TestRoundTrip:
    def test_save_then_load(self, store):
        history = [{"role": "user", "content": "CD19 CAR-T options?"},
                   {"role": "assistant", "content": "Tisagenlecleucel…"}]
        store.save("s1", history)
        assert store.load("s1") == history

    def test_directory_is_created_on_demand(self, store):
        assert not store.directory.exists()
        store.save("s1", [])
        assert store.directory.is_dir()

    def test_unknown_session_is_empty_not_an_error(self, store):
        assert store.load("never-seen") == []

    def test_save_overwrites(self, store):
        store.save("s1", [{"role": "user", "content": "first"}])
        store.save("s1", [{"role": "user", "content": "second"}])
        assert store.load("s1") == [{"role": "user", "content": "second"}]


class TestTTL:
    def _age(self, store, session_id, delta):
        p = store._path(session_id)
        d = json.loads(p.read_text())
        d["updated"] = (datetime.now(timezone.utc) - delta).isoformat()
        p.write_text(json.dumps(d))

    def test_within_ttl_is_returned(self, store):
        store.save("s1", [{"role": "user", "content": "x"}])
        self._age(store, "s1", timedelta(hours=23))
        assert store.load("s1") != []

    def test_expired_is_dropped_and_deleted(self, store):
        store.save("s1", [{"role": "user", "content": "x"}])
        self._age(store, "s1", timedelta(hours=25))
        assert store.load("s1") == []
        assert not store._path("s1").exists(), "an expired conversation must not linger on disk"

    def test_ttl_is_configurable(self, tmp_path):
        s = ConversationStore(tmp_path, ttl=timedelta(seconds=0))
        s.save("s1", [{"role": "user", "content": "x"}])
        assert s.load("s1") == []


class TestCleanup:
    def test_removes_only_expired(self, store):
        store.save("fresh", [{"role": "user", "content": "a"}])
        store.save("stale", [{"role": "user", "content": "b"}])
        TestTTL()._age(store, "stale", timedelta(hours=48))
        store.cleanup_expired()
        assert store._path("fresh").exists()
        assert not store._path("stale").exists()

    def test_missing_directory_is_not_an_error(self, store):
        store.cleanup_expired()          # never created

    def test_one_corrupt_file_does_not_stop_the_sweep(self, store):
        store.save("stale", [{"role": "user", "content": "b"}])
        TestTTL()._age(store, "stale", timedelta(hours=48))
        (store.directory / "corrupt.json").write_text("{not json")
        store.cleanup_expired()
        assert not store._path("stale").exists(), "a corrupt neighbour blocked the sweep"


class TestFailuresDegrade:
    """A conversation that fails to persist must never take down the answer."""

    def test_save_failure_is_swallowed_and_warned(self, tmp_path, caplog):
        blocked = tmp_path / "afile"
        blocked.write_text("x")                       # a file where a directory is needed
        s = ConversationStore(blocked / "conversations")
        with caplog.at_level("WARNING"):
            s.save("s1", [{"role": "user", "content": "x"}])   # must not raise
        assert "s1" in caplog.text

    def test_corrupt_file_loads_as_empty(self, store):
        store.directory.mkdir(parents=True, exist_ok=True)
        store._path("s1").write_text("{not json")
        assert store.load("s1") == []
