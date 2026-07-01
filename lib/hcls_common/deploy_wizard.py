"""
Guided deploy / destroy wizard (A8).

Operational helpers for a single-box deploy: a health-gated status view and an ordered
teardown/deploy plan that enforces **core-last teardown** (tear down models/agents first; bring
down platform services + engines last, so dependents never outlive their dependencies).
"""
from __future__ import annotations

from typing import Any

from hcls_common.capability_registry import CapabilityType, CapabilityRegistry

# lower number = torn down earlier (leaf capabilities first, core/platform last)
_TEARDOWN_PRIORITY = {
    CapabilityType.MODEL: 0,
    CapabilityType.AGENT: 0,
    CapabilityType.NIM: 1,
    CapabilityType.STAGE: 2,
    CapabilityType.SERVICE: 3,
    CapabilityType.ENGINE: 4,     # engines/core last
}


def teardown_order(registry: CapabilityRegistry) -> list[str]:
    """Live capability ids ordered for safe teardown — models/agents first, services+engines last."""
    live = [c for c in registry.all() if c.status.value == "live"]
    return [c.id for c in sorted(live, key=lambda c: _TEARDOWN_PRIORITY.get(c.type, 2))]


def deploy_order(registry: CapabilityRegistry) -> list[str]:
    """Bring-up order — the reverse of teardown (core/services first, then models/agents)."""
    return list(reversed(teardown_order(registry)))


def status(tools: Any) -> dict[str, Any]:
    """Health-gated status for the wizard (uses FactoryTools.check_factory_health)."""
    h = tools.check_factory_health()
    return {
        "up": h["up"],
        "down": h["down"],
        "down_services": [k for k, v in h["services"].items() if not v["up"]],
        "ready": h["down"] == 0,
    }
