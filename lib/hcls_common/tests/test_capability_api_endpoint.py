"""An invocation must be aimed at the API, not at the UI beside it.

`Capability.endpoint` advertises the port a human opens. For every agent in this factory that
is the Streamlit UI — the convention since 2026-08-15 is that the registry advertises the UI
port and the HTTP API is UI+1 — and `invoke_capability` built its POST URL from that field.
Every agent invocation through MCP or the composer was therefore aimed at Streamlit.

The demo transcripts are the corroboration: they exercise 8508, 8522, 8524, 8527, 8529, 8532,
8536, 8539, 8541, 8545 and 8127, each exactly one above the registry's advertised port.

The fix is a separate `api_endpoint` field rather than a +1 rule applied to `endpoint`, because
the non-agent capabilities (genomics :5000, Milvus :19530, the NIMs, the model services)
advertise their API directly and incrementing those would break what works.
"""
import pytest

from hcls_common.capability_registry import Capability, get_registry
from hcls_common.mcp_server import FactoryTools

AGENTS_WITH_SEPARATE_UI = {
    "precision-oncology-agent": ("localhost:8526", "localhost:8527"),
    "cart-intelligence-agent": ("localhost:8521", "localhost:8522"),
    "precision-biomarker-agent": ("localhost:8528", "localhost:8529"),
    "cardiology-intelligence-agent": ("localhost:8126", "localhost:8127"),
    "neurology-intelligence-agent": ("localhost:8535", "localhost:8536"),
    "precision-autoimmune-agent": ("localhost:8531", "localhost:8532"),
    "rare-disease-diagnostic-agent": ("localhost:8544", "localhost:8545"),
    "pharmacogenomics-intelligence-agent": ("localhost:8507", "localhost:8508"),
    "imaging-intelligence-agent": ("localhost:8523", "localhost:8524"),
    "single-cell-intelligence-agent": ("localhost:8540", "localhost:8541"),
    "clinical-trial-intelligence-agent": ("localhost:8538", "localhost:8539"),
}


class TestApiBaseResolution:
    def test_api_base_prefers_the_declared_api(self):
        cap = Capability.from_dict({
            "id": "x", "type": "agent", "name": "x", "description": "x",
            "endpoint": "localhost:8523", "api_endpoint": "localhost:8524",
        })
        assert cap.api_base == "localhost:8524"

    def test_api_base_falls_back_to_the_endpoint(self):
        """Capabilities that advertise their API directly must be left alone."""
        cap = Capability.from_dict({
            "id": "x", "type": "engine", "name": "x", "description": "x",
            "endpoint": "localhost:5000",
        })
        assert cap.api_base == "localhost:5000"

    def test_api_base_is_none_for_a_library(self):
        cap = Capability.from_dict({"id": "x", "type": "service", "name": "x", "description": "x"})
        assert cap.api_base is None

    def test_api_endpoint_survives_a_round_trip(self):
        cap = Capability.from_dict({
            "id": "x", "type": "agent", "name": "x", "description": "x",
            "endpoint": "localhost:8523", "api_endpoint": "localhost:8524",
        })
        assert Capability.from_dict(cap.to_dict()).api_endpoint == "localhost:8524"


class TestTheRegistryDeclaresTheSplit:
    @pytest.mark.parametrize("cap_id,ports", sorted(AGENTS_WITH_SEPARATE_UI.items()))
    def test_each_agent_declares_its_api_port(self, cap_id, ports):
        ui, api = ports
        cap = get_registry().get(cap_id)
        assert cap.endpoint == ui, "the advertised UI port changed"
        assert cap.api_endpoint == api, f"{cap_id} would be invoked on its Streamlit UI"

    def test_non_agent_capabilities_declare_no_separate_api(self):
        """Nothing was guessed: only the entries with a known UI/API split carry the field."""
        for cap in get_registry().all():
            if cap.api_endpoint:
                assert cap.id in AGENTS_WITH_SEPARATE_UI, cap.id

    def test_the_api_port_is_one_above_the_ui_port(self):
        for cap_id, (ui, api) in AGENTS_WITH_SEPARATE_UI.items():
            assert int(api.rsplit(":", 1)[1]) == int(ui.rsplit(":", 1)[1]) + 1, cap_id


class TestInvocationUsesTheApiPort:
    def test_the_posted_url_is_the_api_not_the_ui(self):
        tools = FactoryTools()
        posted = {}

        tools._port_open = lambda host, port: True
        def fake_post(url, payload):
            posted["url"] = url
            return {"status": "ok"}
        tools._http_post = fake_post

        tools.invoke_capability("imaging-intelligence-agent",
                                {"question": "describe the study"}, path="/v1/imaging/query")
        assert posted["url"].startswith("http://localhost:8524/")
        assert "8523" not in posted["url"]

    def test_an_unreachable_api_is_reported_by_its_api_address(self):
        tools = FactoryTools()
        tools._port_open = lambda host, port: False
        out = tools.invoke_capability("imaging-intelligence-agent", {})
        assert out["status"] == "down"
        assert "8524" in out["reason"]

    def test_health_probes_the_api_and_says_which_port_it_probed(self):
        tools = FactoryTools()
        probed = []
        tools._port_open = lambda host, port: (probed.append(port), True)[1]
        report = tools.check_factory_health()
        entry = report["services"]["imaging-intelligence-agent"]
        assert entry["endpoint"] == "localhost:8523"      # still advertises the UI
        assert entry["probed"] == "localhost:8524"        # but probes what it would invoke
        assert 8524 in probed
