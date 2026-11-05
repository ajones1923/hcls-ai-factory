"""A capability marked `live` must have something that starts it.

Found 2026-11-04, while answering a question about how to deploy. Nine capabilities were marked
`live` in the registry, advertised a port as their `endpoint`, and were **absent from
`health-monitor.sh`'s service roster** -- so nothing on the box ever started them, and nothing
ever would. Probed on the live instance, all nine ports were closed:

    8126 cardiology  8507 pharmacogenomics  8521 cart        8523 imaging  8526 oncology
    8531 autoimmune  8538 clinical-trial    8540 single-cell 8544 rare-disease

Each one ships an `app/*_ui.py` and a venv with streamlit in it. They were not broken; they were
never supervised. The roster had eleven agents' APIs and only two of their UIs.

This is the same shape as the recurring finding in the claim-integrity review -- the honest words
were somewhere no generator reads -- except the two surfaces here are the registry and the
supervisor, and the thing they disagreed about was whether a service exists at all. It also
explains why two measurements of "is the factory up" disagreed: `health-monitor.sh status` said
31/32 healthy because it reports on its own roster, while a probe of the registry's advertised
endpoints answered 15 of 24. Both numbers were right about different questions.

LIMIT WORTH STATING: this gate reads the repo, so it catches the roster drifting away from the
registry. It cannot see the box. A correct repo deployed nowhere still serves nothing --
`scripts/check_endpoints.py` against a fresh `probe_live_endpoints.py` snapshot is what covers
that, and `docs/build/SPARK_DEPLOY.md` is the procedure.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
MONITOR = REPO / "health-monitor.sh"
MANIFEST = REPO / "lib" / "hcls_common" / "capabilities.json"

# `"id|port|name|health_path|type|dir|cmd"` -- one roster entry per line.
ROSTER_ENTRY = re.compile(r'^\s*"([a-z0-9-]+)\|(\d+)\|')
ENDPOINT_PORT = re.compile(r":(\d+)$")


def supervised_ports() -> dict[str, str]:
    """port -> service id, from the supervisor's own SERVICES array."""
    if not MONITOR.exists():                       # pragma: no cover - checkout sanity
        pytest.skip(f"no supervisor script at {MONITOR}")
    found = {}
    for line in MONITOR.read_text().splitlines():
        m = ROSTER_ENTRY.match(line)
        if m:
            found[m.group(2)] = m.group(1)
    return found


def live_capabilities_with_a_port() -> list[tuple[str, str]]:
    """(port, capability id) for every `live` capability advertising a host:port endpoint."""
    caps = json.loads(MANIFEST.read_text())["capabilities"]
    out = []
    for c in caps:
        if c.get("status") != "live":
            continue
        m = ENDPOINT_PORT.search(c.get("endpoint") or "")
        if m:
            out.append((m.group(1), c["id"]))
    return out


class TestEveryLiveEndpointIsSupervised:
    def test_the_roster_covers_them_all(self):
        """The finding in one assertion."""
        sup = supervised_ports()
        orphans = [(p, cid) for p, cid in live_capabilities_with_a_port() if p not in sup]
        assert not orphans, (
            "these capabilities are marked `live` and advertise a port that nothing in "
            "health-monitor.sh starts, so the endpoint can never answer:\n"
            + "\n".join(f"  :{p}  {cid}" for p, cid in sorted(orphans))
        )

    def test_the_check_has_something_to_check(self):
        """Guards against the assertion above passing because both sides came back empty.

        An unparseable roster or a renamed manifest key would make this file report green while
        checking nothing, which is the failure mode the whole review program is about.
        """
        assert len(supervised_ports()) >= 30, "roster did not parse; the entry format changed"
        assert len(live_capabilities_with_a_port()) >= 15, "no live endpoints read from manifest"

    def test_the_nine_agent_uis_are_on_the_roster(self):
        """Pinned explicitly, because these are the ones that were missing.

        Named rather than counted: a regression here would most likely drop one, and a count
        would still pass if something unrelated were added in the same change.
        """
        sup = supervised_ports()
        for port in ("8126", "8507", "8521", "8523", "8526", "8531", "8538", "8540", "8544"):
            assert port in sup, f":{port} fell off the supervisor roster again"


class TestTheRosterDoesNotInventServices:
    def test_every_streamlit_entry_points_at_a_file_that_exists(self):
        """A roster entry for a missing app is a start that fails every five minutes forever.

        `start_service` has preflights for a missing directory and an unbuilt venv, but not for
        a mistyped app path -- that one waits the full 60s and logs a generic failure.
        """
        missing = []
        for line in MONITOR.read_text().splitlines():
            if not ROSTER_ENTRY.match(line):
                continue
            parts = line.strip().strip('"').split("|")
            if len(parts) < 7 or parts[4] != "streamlit":
                continue
            svc_dir, cmd = parts[5], parts[6]
            m = re.search(r"streamlit run (\S+)", cmd)
            if not m:
                continue
            base = svc_dir.replace("${SCRIPT_DIR}", str(REPO))
            if not (Path(base) / m.group(1)).exists():
                missing.append(f"  {parts[0]}: {base}/{m.group(1)}")
        assert not missing, "roster entries naming an app file that is not there:\n" + "\n".join(missing)
