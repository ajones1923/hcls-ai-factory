"""The clinician-facing surface is a Streamlit page, and it has to be gated too.

Two findings, both about the UI layer rather than the API layer:

1. Four UIs pointed their API base at a port that is not their API -- cardiology at :8126 (its
   own UI port), rare-disease at :8134 (a docker-compose host port), single-cell at :8540 and
   clinical-trial at :8538 (both UI ports). Every call from those pages got a refused
   connection, so they never reached their gated API at all; and once those UI ports are
   supervised, the same calls would get Streamlit HTML back instead, which is worse.

2. Five UIs import their engine and call it in-process, so the output-honesty gate -- response
   middleware -- never ran on them. Those now call `gate_text`.

These tests read the UI sources, which is crude, but the alternative is importing Streamlit
pages that execute at import time.
"""
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[3]

pytestmark = pytest.mark.skipif(not (REPO / "core").is_dir(),
                                reason="run from a source checkout")

# UI source -> the API port it must call. Verified against health-monitor.sh and the live box.
HTTP_UIS = {
    "core/engines/cardiology/app/cardio_ui.py": 8127,
    "core/agents/rare-disease-diagnostic/app/diagnostic_ui.py": 8545,
    "core/agents/single-cell/app/sc_ui.py": 8541,
    "core/agents/clinical-trial/app/trial_ui.py": 8539,
    "core/agents/neurology/app/neuro_ui.py": 8536,
    "core/engines/precision-oncology/agent/app/oncology_ui.py": 8527,
}

# UIs that call their engine in-process and must therefore gate their own output.
IN_PROCESS_UIS = [
    "core/agents/cart/app/cart_ui.py",
    "core/agents/pharmacogenomics/app/pgx_ui.py",
    "core/agents/precision-biomarker/app/biomarker_ui.py",
    "core/agents/precision-autoimmune/app/autoimmune_ui.py",
    "core/engines/precision-intelligence/app/chat_ui.py",
]

# The agent UI ports, which must never appear as an API base.
UI_PORTS = {8126, 8507, 8521, 8523, 8526, 8528, 8531, 8535, 8538, 8540, 8544}


class TestHttpUisCallTheirApi:
    @pytest.mark.parametrize("path,port", sorted(HTTP_UIS.items()))
    def test_the_default_api_base_is_the_api_port(self, path, port):
        src = (REPO / path).read_text()
        bases = re.findall(r'os\.environ\.get\(\s*"[A-Z_]+"\s*,\s*"http://localhost:(\d+)"', src)
        assert bases, f"{path}: no API base default found"
        for found in bases:
            assert int(found) == port, f"{path} points at :{found}, API is :{port}"

    @pytest.mark.parametrize("path", sorted(HTTP_UIS))
    def test_no_ui_calls_a_streamlit_port(self, path):
        """A UI that POSTs JSON at a Streamlit port gets HTML back and reports nonsense."""
        src = (REPO / path).read_text()
        for found in re.findall(r'"http://localhost:(\d+)"', src):
            assert int(found) not in UI_PORTS, f"{path} calls UI port :{found}"


class TestInProcessUisGateTheirOutput:
    @pytest.mark.parametrize("path", sorted(IN_PROCESS_UIS))
    def test_the_ui_imports_the_gate(self, path):
        src = (REPO / path).read_text()
        assert "gate_text" in src, (
            f"{path} calls its engine in-process, so response middleware never sees its "
            f"output; it must call gate_text itself"
        )

    @pytest.mark.parametrize("path", sorted(IN_PROCESS_UIS))
    def test_the_gate_is_called_and_not_merely_imported(self, path):
        src = (REPO / path).read_text()
        assert re.search(r"gate_text\(", src), f"{path} imports gate_text but never calls it"

    @pytest.mark.parametrize("path", sorted(IN_PROCESS_UIS))
    def test_every_stream_placeholder_is_rewritten_after_gating(self, path):
        """Tokens are already on screen; a withheld claim must not be left displayed."""
        src = (REPO / path).read_text()
        if "generate_stream" not in src and "query_stream" not in src:
            pytest.skip("no streaming path in this UI")
        gated = [m.start() for m in re.finditer(r"gate_text\(", src)]
        assert gated, f"{path} streams but never gates"
