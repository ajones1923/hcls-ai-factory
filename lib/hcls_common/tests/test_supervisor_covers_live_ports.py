"""A `live` capability needs something that actually starts it.

`health-monitor.sh` is the deployment mechanism here -- a `*/5` cron that starts anything it
finds down. On 2026-11-03 all eleven agents shipped an `app/*_ui.py`, the registry advertised
each one's UI port as the capability `endpoint`, and only TWO were supervised. Probed on the
live box: :8126, :8507, :8521, :8523, :8526, :8531, :8538, :8540 and :8544 were all closed. So
nine capabilities marked `live` advertised an endpoint nothing served.

`scripts/validate_registry.py` could not catch it. It checked that every supervised port is one
the registry allocates, which is the other direction entirely; a port the registry allocated and
nobody supervised was invisible to it. These tests pin both directions, and they parse the same
two files the validator does so the guard cannot be satisfied by a comment.
"""
import json
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[3]
HEALTH_MONITOR = REPO / "health-monitor.sh"
REGISTRY = REPO / "lib" / "hcls_common" / "capabilities.json"

pytestmark = pytest.mark.skipif(
    not (HEALTH_MONITOR.exists() and REGISTRY.exists()),
    reason="run from a source checkout; health-monitor.sh / capabilities.json not found",
)

# Ports that are legitimately supervised without being registry capabilities: the landing page,
# the two shared Streamlit UIs of the intelligence/discovery engines, the orchestrator portal,
# and the monitoring stack.
INFRA_PORTS = {3000, 8080, 8501, 8510, 9099, 9100, 9400, 19530}


def supervised() -> dict:
    out = {}
    for line in HEALTH_MONITOR.read_text().splitlines():
        m = re.match(r'\s*"([a-z0-9-]+)\|(\d+)\|', line)
        if m:
            out[m.group(1)] = int(m.group(2))
    return out


def capabilities() -> list:
    data = json.loads(REGISTRY.read_text())
    return data["capabilities"] if isinstance(data, dict) else data


def port_of(endpoint) -> int | None:
    m = re.search(r":(\d+)$", str(endpoint or ""))
    return int(m.group(1)) if m else None


class TestEveryLiveEndpointIsSupervised:
    def test_no_live_capability_advertises_an_unsupervised_endpoint(self):
        ports = set(supervised().values())
        orphans = [
            (c["id"], field, c[field])
            for c in capabilities()
            if c.get("status") == "live"
            for field in ("endpoint", "api_endpoint")
            if c.get(field) and port_of(c[field]) not in ports
        ]
        assert not orphans, (
            "these live capabilities advertise a port nothing starts or restarts: " + repr(orphans)
        )

    def test_all_eleven_agent_uis_are_supervised(self):
        """The nine that were missing, plus the two that were not."""
        ports = set(supervised().values())
        for port in (8126, 8507, 8521, 8523, 8526, 8528, 8531, 8535, 8538, 8540, 8544):
            assert port in ports, f"agent UI :{port} is not supervised"

    def test_all_eleven_agent_apis_are_supervised(self):
        ports = set(supervised().values())
        for port in (8127, 8508, 8522, 8524, 8527, 8529, 8532, 8536, 8539, 8541, 8545):
            assert port in ports, f"agent API :{port} is not supervised"


class TestTheOtherDirectionStillHolds:
    def test_every_supervised_port_is_allocated_or_infrastructure(self):
        allocated = set()
        for c in capabilities():
            for field in ("endpoint", "api_endpoint"):
                p = port_of(c.get(field))
                if p:
                    allocated.add(p)
        strays = {name: p for name, p in supervised().items()
                  if p not in allocated and p not in INFRA_PORTS}
        assert not strays, f"supervised but not in the registry: {strays}"


class TestTheRosterIsWellFormed:
    def test_no_port_is_supervised_twice(self):
        """Two entries on one port cannot both bind it; one would restart-loop forever."""
        sup = supervised()
        seen: dict = {}
        for name, port in sup.items():
            seen.setdefault(port, []).append(name)
        dupes = {p: names for p, names in seen.items() if len(names) > 1}
        assert not dupes, f"two supervisor entries claim the same port: {dupes}"

    def test_every_ui_entry_runs_a_streamlit_app_that_exists(self):
        for line in HEALTH_MONITOR.read_text().splitlines():
            m = re.match(r'\s*"([a-z0-9-]+)\|\d+\|[^|]*\|[^|]*\|streamlit\|([^|]+)\|(.+)"', line)
            if not m:
                continue
            name, directory, command = m.groups()
            directory = directory.replace("${SCRIPT_DIR}/", "")
            app = re.search(r"run (\S+\.py)", command)
            assert app, f"{name}: streamlit entry with no app path"
            assert (REPO / directory / app.group(1)).exists(), \
                f"{name}: {directory}/{app.group(1)} does not exist"

    def test_every_ui_entry_binds_loopback_only(self):
        """A clinical UI must not be reachable off-box by default."""
        for line in HEALTH_MONITOR.read_text().splitlines():
            if "|streamlit|" not in line:
                continue
            assert "--server.address 127.0.0.1" in line, line.strip()
