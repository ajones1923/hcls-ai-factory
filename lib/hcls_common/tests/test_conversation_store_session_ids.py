"""A session id arrives in a request payload and must not choose a path on disk.

`_path` interpolated it straight into a filename:

    return self.directory / f"{session_id}.json"

so "../../etc/cron.d/pwn" wrote JSON outside the conversation directory, and a traversing id
passed to `load` read a file that is not a conversation. Every operation in this class swallows
its exception by design -- a conversation that fails to persist must never take down the answer
the user asked for -- so both failures were silent.
"""
import pytest

from hcls_common.conversation_store import ConversationStore, safe_session_name


@pytest.fixture
def store(tmp_path):
    return ConversationStore(tmp_path / "conversations")


HOSTILE = [
    "../../etc/cron.d/pwn",
    "../outside",
    "a/../../b",
    "/etc/passwd",
    "..",
    ".",
    ".hidden",
    "with spaces",
    "semi;colon",
    "nul\x00byte",
]


class TestNothingIsWrittenOutsideTheDirectory:
    @pytest.mark.parametrize("sid", HOSTILE)
    def test_the_file_lands_in_the_conversation_directory(self, store, sid):
        store.save(sid, [{"role": "user", "content": "hello"}])
        root = store.directory.resolve()
        written = list(root.glob("**/*.json"))
        assert written, f"{sid!r} wrote nothing at all"
        for path in written:
            assert path.resolve().parent == root, f"{sid!r} escaped to {path}"

    @pytest.mark.parametrize("sid", HOSTILE)
    def test_a_hostile_id_round_trips_to_its_own_history(self, store, sid):
        """Sanitizing must not break the feature it protects."""
        store.save(sid, [{"role": "user", "content": sid}])
        assert store.load(sid) == [{"role": "user", "content": sid}]

    def test_nothing_is_created_above_the_directory(self, store, tmp_path):
        store.save("../../escaped", [{"role": "user", "content": "x"}])
        assert not (tmp_path / "escaped.json").exists()
        assert not (tmp_path.parent / "escaped.json").exists()


class TestOrdinaryIdsAreUnchanged:
    @pytest.mark.parametrize("sid", ["abc123", "session-1", "a.b_c-d", "A1"])
    def test_a_safe_id_keeps_its_filename(self, store, sid):
        """Existing conversation files must keep their names."""
        assert safe_session_name(sid) == sid
        store.save(sid, [{"role": "user", "content": "hi"}])
        assert (store.directory / f"{sid}.json").exists()

    def test_a_safe_id_round_trips(self, store):
        store.save("abc123", [{"role": "assistant", "content": "answer"}])
        assert store.load("abc123")[0]["content"] == "answer"


class TestDistinctIdsStayDistinct:
    def test_two_hostile_ids_do_not_share_a_file(self, store):
        """A sanitizer that collapsed ids would serve one session's history to another."""
        assert safe_session_name("../a") != safe_session_name("../b")

    def test_two_hostile_ids_keep_separate_histories(self, store):
        store.save("../a", [{"role": "user", "content": "first"}])
        store.save("../b", [{"role": "user", "content": "second"}])
        assert store.load("../a")[0]["content"] == "first"
        assert store.load("../b")[0]["content"] == "second"

    def test_the_mapping_is_deterministic(self):
        assert safe_session_name("../x") == safe_session_name("../x")

    @pytest.mark.parametrize("sid", [None, ""])
    def test_an_empty_id_still_gets_a_name(self, sid):
        assert safe_session_name(sid).startswith("sid-")


class TestLoadCannotReadArbitraryFiles:
    def test_a_traversing_id_does_not_read_a_file_outside(self, store, tmp_path):
        decoy = tmp_path / "secret.json"
        decoy.write_text('{"updated": "2999-01-01T00:00:00+00:00", "messages": [{"x": 1}]}')
        assert store.load("../secret") == []
