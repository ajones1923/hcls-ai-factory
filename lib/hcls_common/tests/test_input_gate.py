"""A1/A2: param-contract + input-validation gate."""
from hcls_common.capability_registry import (
    Capability, CapabilityType, Port, ValueShape, Status,
    validate_inputs, inputs_ok, get_registry,
)
from hcls_common.mcp_server import FactoryTools


def _cap():
    return Capability(
        id="x", type=CapabilityType.MODEL, name="X", description="",
        inputs=[
            Port("mode", ValueShape.SCALAR, enum=["fast", "accurate"]),
            Port("n", ValueShape.SCALAR, required=False, minimum=1, maximum=100, default=10),
            Port("seq", ValueShape.SCALAR, required=True),
        ],
    )


class TestValidate:
    def test_default_applied(self):
        out, issues = validate_inputs(_cap(), {"mode": "fast", "seq": "M"})
        assert out["n"] == 10 and inputs_ok(issues)

    def test_required_missing_is_error(self):
        out, issues = validate_inputs(_cap(), {"mode": "fast"})
        assert not inputs_ok(issues) and any("required" in i and "seq" in i for i in issues)

    def test_enum_violation_is_error(self):
        out, issues = validate_inputs(_cap(), {"mode": "turbo", "seq": "M"})
        assert not inputs_ok(issues) and any("not in allowed" in i for i in issues)

    def test_numeric_clamped_with_warning(self):
        out, issues = validate_inputs(_cap(), {"mode": "fast", "seq": "M", "n": 999})
        assert out["n"] == 100 and inputs_ok(issues)      # clamp is a WARN, not blocking
        assert any(i.startswith("WARN") and "clamped" in i for i in issues)

    def test_real_manifest_contract(self):
        r = get_registry(reload=True)
        out, issues = validate_inputs(r.get("molecule-generator"), {"seeds": ["CCO"], "n": 5000})
        assert out["n"] == 100                            # clamped to the manifest maximum


class TestGate:
    def test_invoke_rejects_invalid_before_http(self):
        called = []
        ft = FactoryTools(get_registry(reload=True),
                          port_open=lambda h, p: True,
                          http_post=lambda u, pl: called.append(u) or {"status": "ok"})
        # singlecell-compute now requires expression_matrix; omit it -> gate blocks before HTTP
        r = ft.invoke_capability("singlecell-compute", {"resolution": 1.0})
        assert r["status"] == "invalid_input" and not called
