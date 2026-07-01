"""A7 license/pin gate + A8 deploy/destroy wizard."""
from hcls_common.license_gate import classify_license, is_pinned, check_requirements, audit
from hcls_common.deploy_wizard import teardown_order, deploy_order, status
from hcls_common.capability_registry import CapabilityType
from hcls_common import get_registry


class TestLicenseGate:
    def test_classify(self):
        assert classify_license("Apache-2.0") == "permissive"
        assert classify_license("MIT") == "permissive"
        assert classify_license("CC-BY-NC-SA-4.0") == "non-commercial"
        assert classify_license("GPL-3.0") == "copyleft"
        assert classify_license("CC-BY-SA-3.0") == "copyleft"
        assert classify_license(None) == "unknown"

    def test_pinning(self):
        assert is_pinned("torch==2.12.1") and not is_pinned("torch>=2.0") and not is_pinned("torch")
        assert check_requirements(["a==1", "b>=2", "# c", "d"]) == ["UNPINNED: b>=2", "UNPINNED: d"]

    def test_audit_blocks_noncommercial_and_unpinned(self):
        pkgs = [{"name": "esm", "license": "MIT"}, {"name": "nt", "license": "CC-BY-NC-SA"}]
        rep = audit(pkgs, requirements=["torch==2.12", "foo>=1"])
        assert not rep["ok"]
        assert "nt" in rep["blocking"] and any("foo" in b for b in rep["blocking"])

    def test_audit_passes_clean(self):
        rep = audit([{"name": "a", "license": "Apache-2.0"}], requirements=["a==1.0"])
        assert rep["ok"] and rep["blocking"] == []


class TestDeployWizard:
    def test_teardown_core_last(self):
        r = get_registry(reload=True)
        order = teardown_order(r)
        types = {c.id: c.type for c in r.all()}
        # the last torn-down item must be an engine/service, the first a model/agent
        assert types[order[0]] in (CapabilityType.MODEL, CapabilityType.AGENT, CapabilityType.NIM)
        assert types[order[-1]] in (CapabilityType.ENGINE, CapabilityType.SERVICE, CapabilityType.STAGE)
        assert deploy_order(r) == list(reversed(order))

    def test_status(self):
        class FakeTools:
            def check_factory_health(self):
                return {"up": 2, "down": 1, "services": {"a": {"up": True}, "b": {"up": True}, "c": {"up": False}}}
        s = status(FakeTools())
        assert s["up"] == 2 and s["down_services"] == ["c"] and not s["ready"]
