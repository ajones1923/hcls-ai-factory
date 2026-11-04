"""The endpoint ratchet must fail for the reasons a ratchet exists, and not for other ones.

Check 1 of the claim-integrity rubric: a capability registered `live` answers on the port the
registry advertises. On 2026-11-03 nine did not, and nothing in the repository could notice,
because the static guard (`test_supervisor_covers_live_ports.py`) checks that a capability is
*supervised*, which a crash-looping service also is.

These tests drive the real script as CI drives it, rather than importing a function out of it:
the thing that can silently stop gating is the exit code, and a test that never runs the
command cannot see one. The precedent is in this repo -- the test-depth ratchet raised
`KeyError: 'orchestrator'` on every invocation for long enough that nobody noticed it had never
run, and its badge-guard sibling could not fail at all, because a regex never matched and the
no-match branch returned PASS.
"""
import json
import subprocess
import sys
import pathlib

import pytest

SCRIPT = (pathlib.Path(__file__).resolve().parents[3] / "scripts" / "check_endpoints.py")

pytestmark = pytest.mark.skipif(not SCRIPT.exists(),
                                reason="run from a source checkout; scripts/ not found")


def _registry(tmp_path, *ids, status="live", endpoint="localhost:8126"):
    p = tmp_path / "capabilities.json"
    p.write_text(json.dumps({"capabilities": [
        {"id": i, "status": status, "endpoint": endpoint} for i in ids]}))
    return p


def _snapshot(tmp_path, results, measured_at="2026-11-04T12:00:00+00:00"):
    p = tmp_path / "snapshot.json"
    p.write_text(json.dumps({"measured_at": measured_at, "host": "test", "results": results}))
    return p


def _baseline(tmp_path, known_down):
    p = tmp_path / "baseline.json"
    p.write_text(json.dumps({"known_down": known_down}))
    return p


def run(registry, snapshot, baseline, *extra):
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--enforce", "--max-age-days", "36500",
         "--registry", str(registry), "--snapshot", str(snapshot),
         "--baseline", str(baseline), *extra],
        capture_output=True, text=True)
    return r.returncode, r.stdout


UP = {"reachable": True, "how": "http", "status": 200, "endpoint": "localhost:8126"}
DOWN = {"reachable": False, "how": "tcp", "detail": "Connection refused",
        "endpoint": "localhost:8126"}


class TestItCatchesWhatItIsFor:
    def test_a_capability_that_stopped_answering_fails(self, tmp_path):
        code, out = run(_registry(tmp_path, "cardiology-intelligence-agent"),
                        _snapshot(tmp_path, {"cardiology-intelligence-agent": DOWN}),
                        _baseline(tmp_path, []))
        assert code == 1
        assert "cardiology-intelligence-agent" in out

    def test_a_live_capability_never_probed_fails(self, tmp_path):
        """An unprobed endpoint is unknown, not up. This is the registry-edit case."""
        code, out = run(_registry(tmp_path, "brand-new-agent"),
                        _snapshot(tmp_path, {"some-other-agent": UP}),
                        _baseline(tmp_path, []))
        assert code == 1
        assert "never been probed" in out
        assert "brand-new-agent" in out

    def test_a_stale_probe_fails(self, tmp_path):
        """A frozen green probe would certify an instance nobody has looked at."""
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--enforce", "--max-age-days", "30",
             "--registry", str(_registry(tmp_path, "a")),
             "--snapshot", str(_snapshot(tmp_path, {"a": UP}, "2020-01-01T00:00:00+00:00")),
             "--baseline", str(_baseline(tmp_path, []))],
            capture_output=True, text=True)
        assert r.returncode == 1
        assert "days old" in r.stdout

    def test_a_missing_snapshot_fails(self, tmp_path):
        code, out = run(_registry(tmp_path, "a"), tmp_path / "nope.json",
                        _baseline(tmp_path, []))
        assert code == 1
        assert "no endpoint snapshot" in out


class TestItDoesNotFailForOtherReasons:
    def test_a_baselined_outage_passes(self, tmp_path):
        """The nine that v2.2.0 has not been deployed to fix are accepted debt, not a block."""
        code, out = run(_registry(tmp_path, "cardiology-intelligence-agent"),
                        _snapshot(tmp_path, {"cardiology-intelligence-agent": DOWN}),
                        _baseline(tmp_path, ["cardiology-intelligence-agent"]))
        assert code == 0
        assert "no new capability stopped answering" in out

    def test_a_recovered_capability_passes_and_says_so(self, tmp_path):
        """Failing a build because something got fixed teaches people to disable the gate."""
        code, out = run(_registry(tmp_path, "a"), _snapshot(tmp_path, {"a": UP}),
                        _baseline(tmp_path, ["a"]))
        assert code == 0
        assert "now answer" in out
        assert "--update-baseline" in out

    def test_everything_up_and_no_baseline_passes(self, tmp_path):
        code, out = run(_registry(tmp_path, "a"), _snapshot(tmp_path, {"a": UP}),
                        _baseline(tmp_path, []))
        assert code == 0
        assert "answers on its advertised port" in out

    def test_a_planned_capability_is_not_required_to_answer(self, tmp_path):
        """`planned` is an honest status, not a missing service."""
        code, out = run(_registry(tmp_path, "future-agent", status="planned"),
                        _snapshot(tmp_path, {}), _baseline(tmp_path, []))
        assert code == 0, out

    def test_a_live_capability_with_no_endpoint_is_out_of_scope(self, tmp_path):
        """honesty-gate and event-bus are library surfaces; they serve no port by design."""
        reg = tmp_path / "capabilities.json"
        reg.write_text(json.dumps({"capabilities": [
            {"id": "honesty-gate", "status": "live", "endpoint": None}]}))
        code, out = run(reg, _snapshot(tmp_path, {"x": UP}), _baseline(tmp_path, []))
        assert code == 0, out


class TestTheBaselineCannotRotQuietly:
    def test_update_baseline_records_exactly_what_is_down(self, tmp_path):
        base = _baseline(tmp_path, [])
        subprocess.run(
            [sys.executable, str(SCRIPT), "--update-baseline", "--reason", "awaiting deploy",
             "--registry", str(_registry(tmp_path, "a", "b")),
             "--snapshot", str(_snapshot(tmp_path, {"a": DOWN, "b": UP})),
             "--baseline", str(base)], capture_output=True, text=True, check=True)
        recorded = json.loads(base.read_text())
        assert recorded["known_down"] == ["a"]
        assert recorded["reason"] == "awaiting deploy"

    def test_accepting_an_outage_without_a_reason_is_refused(self, tmp_path):
        """An unexplained baseline cannot be told from a gate somebody switched off."""
        base = tmp_path / "fresh.json"
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--update-baseline",
             "--registry", str(_registry(tmp_path, "a")),
             "--snapshot", str(_snapshot(tmp_path, {"a": DOWN})),
             "--baseline", str(base)], capture_output=True, text=True)
        assert r.returncode == 1
        assert "no reason recorded" in r.stdout
        assert not base.exists()

    def test_recording_zero_outages_needs_no_reason(self, tmp_path):
        """Nothing to excuse, nothing to explain."""
        base = tmp_path / "fresh.json"
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--update-baseline",
             "--registry", str(_registry(tmp_path, "a")),
             "--snapshot", str(_snapshot(tmp_path, {"a": UP})),
             "--baseline", str(base)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout
        assert json.loads(base.read_text())["known_down"] == []

    def test_the_baseline_does_not_excuse_an_unprobed_capability(self, tmp_path):
        """Listing an id as known-down must not stand in for measuring it."""
        code, out = run(_registry(tmp_path, "a"), _snapshot(tmp_path, {}),
                        _baseline(tmp_path, ["a"]))
        assert code == 1
        assert "never been probed" in out
