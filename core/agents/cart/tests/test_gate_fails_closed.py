"""If an API key is configured, the service must not start without its gate.

The auth and governance install was wrapped in `except Exception: pass` in all eleven
services, with the comment "never take a running service down over the gate". That failed
OPEN: when hcls_common was missing or install raised, the service booted with no
authentication and no governance, silently, while the comment directly above promised it
fails closed on every route but /health.

The decision now rests on whether a key is configured. One is set -> the operator expects the
gate, and refusing to start beats serving clinical decision support unauthenticated. None set
-> the gate was a no-op anyway, so it logs loudly and carries on.
"""
import re
from pathlib import Path

import pytest

SERVICE_MAINS = sorted(
    p for p in Path(__file__).resolve().parents[3].glob("*/*/api/main.py")
    if "install_api_key_auth" in p.read_text()
) + sorted(
    p for p in Path(__file__).resolve().parents[3].glob("*/*/*/api/main.py")
    if "install_api_key_auth" in p.read_text()
)


def test_the_services_were_found():
    """Guard against this file silently testing nothing if the layout moves."""
    assert len(SERVICE_MAINS) >= 8, [str(p) for p in SERVICE_MAINS]


def _gate_block(src: str) -> str:
    """The region around the auth/governance install.

    Scoped deliberately: these files contain unrelated `except Exception: pass` handlers for
    optional features, and a whole-file search would fail on those instead.
    """
    start = src.index("install_api_key_auth(app")
    return src[start:start + 1400]


@pytest.mark.parametrize("main", SERVICE_MAINS, ids=lambda p: p.parts[-4])
def test_gate_install_does_not_swallow_silently(main):
    block = _gate_block(main.read_text())
    assert "except Exception:  # never take a running service down over the gate" not in block
    assert not re.search(r"except Exception:\s*\n\s*pass", block), (
        f"{main} still swallows a gate-install failure"
    )


@pytest.mark.parametrize("main", SERVICE_MAINS, ids=lambda p: p.parts[-4])
def test_gate_install_refuses_to_start_when_a_key_is_configured(main):
    src = _gate_block(main.read_text())
    assert 'os.environ.get("HCLS_API_KEY")' in src, (
        f"{main} does not check for a configured key before continuing without the gate"
    )
    assert "refusing to start" in src


@pytest.mark.parametrize("main", SERVICE_MAINS, ids=lambda p: p.parts[-4])
def test_governance_is_installed_exactly_once(main):
    """cart mounted the middleware twice, so every claim check ran twice per answer."""
    src = main.read_text()
    calls = len(re.findall(r"^\s*install_governance\(app", src, re.M))
    assert calls == 1, f"{main} calls install_governance {calls} times"
