"""Tests for the MCP tool-surface core (A2). Network calls are injected — no live services needed."""
from hcls_common.capability_registry import get_registry
from hcls_common.mcp_server import FactoryTools


def make_tools(up_ports=None, post_result=None):
    up = set(up_ports or [])
    calls = []

    def port_open(host, port):
        return port in up

    def http_post(url, payload):
        calls.append((url, payload))
        return post_result if post_result is not None else {"status": "ok", "result": {"echo": payload}}

    ft = FactoryTools(get_registry(reload=True), port_open=port_open, http_post=http_post)
    ft._calls = calls  # type: ignore[attr-defined]
    return ft


# ── discovery ────────────────────────────────────────────────────────────────
class TestDiscovery:
    def test_list_all(self):
        ft = make_tools()
        assert len(ft.list_capabilities()) >= 20

    def test_filter_by_type_and_status(self):
        ft = make_tools()
        agents = ft.list_capabilities(type="agent")
        assert len(agents) >= 10 and all(a["type"] == "agent" for a in agents)
        planned = ft.list_capabilities(status="planned")
        assert any(p["id"] == "genmol-nim" for p in planned)

    def test_describe_known_and_unknown(self):
        ft = make_tools()
        d = ft.describe_capability("genomics-engine")
        assert d["name"] == "Genomics Foundations Engine" and d["outputs"]
        assert "error" in ft.describe_capability("nope")


# ── health ───────────────────────────────────────────────────────────────────
class TestHealth:
    def test_health_reports_up_down(self):
        ft = make_tools(up_ports={5001, 19530})   # pretend only RAG API + Milvus are up
        rep = ft.check_factory_health()
        assert rep["up"] == 2 and rep["down"] >= 1
        assert rep["services"]["precision-intelligence-engine"]["up"] is True
        assert rep["services"]["genomics-engine"]["up"] is False


# ── invocation (honesty) ─────────────────────────────────────────────────────
class TestInvoke:
    def test_refuses_planned_capability(self):
        ft = make_tools(up_ports={9999})
        r = ft.invoke_capability("genmol-nim", {"seed_smiles": "CCO"})
        assert r["status"] == "unavailable" and "not live" in r["reason"]

    def test_library_has_no_endpoint(self):
        ft = make_tools()
        assert ft.invoke_capability("honesty-gate")["status"] == "no_endpoint"

    def test_down_endpoint_reported(self):
        ft = make_tools(up_ports=set())            # nothing up
        assert ft.invoke_capability("molmim-nim", {"seed_smiles": "CCO"})["status"] == "down"

    def test_live_up_endpoint_calls_http(self):
        ft = make_tools(up_ports={8001})           # MolMIM NIM "up"
        r = ft.invoke_capability("molmim-nim", {"seed_smiles": "CCO"}, path="generate")
        assert r["status"] == "ok"
        assert ft._calls and ft._calls[0][0].endswith(":8001/generate")


# ── planning (composer seam) ─────────────────────────────────────────────────
class TestPlan:
    def test_plan_surfaces_relevant_capabilities_and_wiring(self):
        ft = make_tools()
        plan = ft.plan_pipeline("predict a protein structure and dock a small molecule")
        ids = {c["id"] for c in plan["relevant_capabilities"]}
        assert ids & {"esmfold-model", "diffdock-nim", "therapeutic-discovery-engine"}
        # structure-shape edge esmfold -> diffdock should be proposed when both are relevant
        if {"esmfold-model", "diffdock-nim"} <= ids:
            assert any(e["shape"] == "structure" for e in plan["candidate_wiring"])
