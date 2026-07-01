"""A3 external-tool ingestion + A4 shared executor."""
from hcls_common.capability_registry import (
    CapabilityRegistry, Serving, validate_inputs, inputs_ok)
from hcls_common.external_tools import ingest_external_tools, external_tool_to_capability
from hcls_common.mcp_server import FactoryTools
from hcls_common import get_registry

FAKE_TOOLS = [
    {"name": "get_weather", "description": "Get weather", "inputSchema": {"type": "object",
        "properties": {"city": {"type": "string"},
                       "units": {"type": "string", "enum": ["c", "f"], "default": "c"},
                       "days": {"type": "integer", "minimum": 1, "maximum": 14, "default": 3}},
        "required": ["city"]}},
    {"name": "list_alerts", "description": "Alerts", "inputSchema": {"type": "object",
        "properties": {"region": {"type": "string"}}}},
]


class TestA3Ingestion:
    def test_ingest_creates_capabilities(self):
        r = CapabilityRegistry()
        ids = ingest_external_tools(r, "weathersvc", "weather.example:443", FAKE_TOOLS)
        assert ids == ["ext.weathersvc.get_weather", "ext.weathersvc.list_alerts"]
        cap = r.get("ext.weathersvc.get_weather")
        assert cap.serving is Serving.EXTERNAL and cap.domain == "external"
        assert {p.name for p in cap.inputs} == {"city", "units", "days"}

    def test_json_schema_populates_a1_contract(self):
        cap = external_tool_to_capability("s", "e:443", FAKE_TOOLS[0])
        units = next(p for p in cap.inputs if p.name == "units")
        days = next(p for p in cap.inputs if p.name == "days")
        city = next(p for p in cap.inputs if p.name == "city")
        assert units.enum == ["c", "f"] and units.default == "c"
        assert days.minimum == 1 and days.maximum == 14
        assert city.required

    def test_a2_gate_applies_to_external_tool(self):
        cap = external_tool_to_capability("s", "e:443", FAKE_TOOLS[0])
        _, bad = validate_inputs(cap, {"city": "NYC", "units": "kelvin"})
        assert not inputs_ok(bad)                              # bad enum rejected
        out, ok = validate_inputs(cap, {"city": "NYC", "days": 99})
        assert inputs_ok(ok) and out["days"] == 14 and out["units"] == "c"   # clamp + default


class TestA4SharedExecutor:
    def test_execute_aliases_invoke(self):
        ft = FactoryTools(get_registry(reload=True), port_open=lambda h, p: False)
        assert ft.execute_capability("esmfold-model", {}) == ft.invoke_capability("esmfold-model", {})

    def test_composer_dispatches_through_shared_executor(self):
        from hcls_common.workflow_composer import WorkflowComposer, Pipeline, Node, NodeInput
        calls = []

        class FakeTools:
            def invoke_capability(self, cap, payload, path=None):
                calls.append(cap); return {"status": "ok", "result": {}}
        comp = WorkflowComposer(get_registry(reload=True), tools=FakeTools())
        pipe = Pipeline(goal="t", nodes=[Node(id="n1", capability="esmfold-model",
                                              inputs=[NodeInput(name="sequence", value="MKT")])])
        comp.run(pipe)
        assert calls == ["esmfold-model"]                     # composer used the shared path
