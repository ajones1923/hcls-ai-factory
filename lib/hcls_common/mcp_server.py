"""
Assistant Tool-Surface (A2) — exposes the HCLS AI Factory's Capability Registry as
tools an AI assistant (Claude, Cursor, AI Playground, custom agents) can call to
*drive the factory*.

Two layers:
  * ``FactoryTools`` — transport-agnostic core (discovery, health, invocation,
    pipeline planning). Fully unit-testable without a running MCP client or live
    services (the network calls are injectable).
  * ``build_server()`` — wraps ``FactoryTools`` as a FastMCP server (stdio or
    streamable-HTTP transport).

Honesty rules carried through from the registry: a capability whose status is not
``live`` is reported as unavailable rather than invoked; mock-served capabilities are
labeled. Nothing here pretends a planned capability is real.
"""
from __future__ import annotations

import socket
from typing import Any, Callable

from hcls_common.capability_registry import (
    CapabilityRegistry,
    Status,
    ValueShape,
    get_registry,
    validate_inputs,
    inputs_ok,
)


# --------------------------------------------------------------------------- #
# Default network adapters (injectable for tests)
# --------------------------------------------------------------------------- #
def _default_port_open(host: str, port: int, timeout: float = 1.5) -> bool:
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except OSError:
        return False


def _default_http_post(url: str, payload: dict[str, Any], timeout: float = 60.0) -> dict[str, Any]:
    import httpx
    r = httpx.post(url, json=payload, timeout=timeout)
    r.raise_for_status()
    try:
        return {"status": "ok", "result": r.json()}
    except Exception:
        return {"status": "ok", "result": r.text}


def _split_endpoint(endpoint: str) -> tuple[str, int]:
    host, _, port = endpoint.partition(":")
    return host or "localhost", int(port or "80")


# --------------------------------------------------------------------------- #
# Core
# --------------------------------------------------------------------------- #
class FactoryTools:
    """Registry-driven tools for discovering, health-checking, and invoking the factory."""

    def __init__(
        self,
        registry: CapabilityRegistry | None = None,
        *,
        port_open: Callable[[str, int], bool] | None = None,
        http_post: Callable[[str, dict], dict] | None = None,
    ) -> None:
        self.registry = registry or get_registry()
        self._port_open = port_open or _default_port_open
        self._http_post = http_post or _default_http_post

    # -- discovery ---------------------------------------------------------- #
    def list_capabilities(
        self, type: str | None = None, domain: str | None = None, status: str | None = None
    ) -> list[dict[str, Any]]:
        """List capabilities, optionally filtered by type/domain/status."""
        caps = self.registry.find(type=type, domain=domain, status=status)
        return [
            {
                "id": c.id, "name": c.name, "type": c.type.value, "domain": c.domain,
                "status": c.status.value, "gpu": c.gpu, "description": c.description,
                "endpoint": c.endpoint,
                "inputs": [{"name": p.name, "shape": p.shape.value, "required": p.required} for p in c.inputs],
                "outputs": [{"name": p.name, "shape": p.shape.value} for p in c.outputs],
            }
            for c in caps
        ]

    def describe_capability(self, capability_id: str) -> dict[str, Any]:
        """Full detail for one capability, including its call contract."""
        try:
            return self.registry.get(capability_id).to_dict()
        except KeyError:
            return {"error": f"unknown capability: {capability_id!r}",
                    "known": self.registry.ids()}

    # -- health ------------------------------------------------------------- #
    def check_factory_health(self) -> dict[str, Any]:
        """Ping every live capability that exposes an endpoint; report up/down."""
        report: dict[str, Any] = {}
        for c in self.registry.find(status="live"):
            if not c.endpoint:
                continue
            host, port = _split_endpoint(c.endpoint)
            report[c.id] = {"endpoint": c.endpoint, "up": self._port_open(host, port)}
        up = sum(1 for v in report.values() if v["up"])
        return {"checked": len(report), "up": up, "down": len(report) - up, "services": report}

    # -- invocation --------------------------------------------------------- #
    def invoke_capability(
        self, capability_id: str, payload: dict[str, Any] | None = None, path: str | None = None
    ) -> dict[str, Any]:
        """Invoke a capability's HTTP endpoint. Honest about availability."""
        try:
            cap = self.registry.get(capability_id)
        except KeyError:
            return {"status": "error", "reason": f"unknown capability: {capability_id!r}"}
        if cap.status is not Status.LIVE:
            return {"status": "unavailable",
                    "reason": f"'{capability_id}' is {cap.status.value}, not live — cannot invoke a non-real capability"}
        if not cap.endpoint:
            return {"status": "no_endpoint",
                    "reason": f"'{capability_id}' has no callable endpoint (library/in-process)"}
        host, port = _split_endpoint(cap.endpoint)
        if not self._port_open(host, port):
            return {"status": "down", "reason": f"endpoint {cap.endpoint} is not reachable"}
        # A2 input-validation gate: reject bad enums / missing required, clamp numerics
        payload, issues = validate_inputs(cap, payload)
        if not inputs_ok(issues):
            return {"status": "invalid_input", "issues": issues}
        use_path = path or getattr(cap, "invoke_path", None) or "/"
        url = f"http://{cap.endpoint.rstrip('/')}/{use_path.lstrip('/')}"
        try:
            res = self._http_post(url, payload or {})
            warns = [i for i in issues if i.startswith("WARN")]
            return {**res, "_input_warnings": warns} if (warns and isinstance(res, dict)) else res
        except Exception as e:  # noqa: BLE001 - surface as data, never crash the tool
            return {"status": "error", "reason": f"{type(e).__name__}: {e}"}

    # -- A4: the single shared executor ------------------------------------- #
    def execute_capability(self, capability_id: str, payload: dict[str, Any] | None = None,
                           path: str | None = None) -> dict[str, Any]:
        """Canonical execute path — alias of invoke_capability so the tool-surface and the
        workflow composer never diverge (both inherit the A2 input-validation gate)."""
        return self.invoke_capability(capability_id, payload, path)

    # -- planning (precursor to the workflow composer) ---------------------- #
    def plan_pipeline(self, goal: str) -> dict[str, Any]:
        """Suggest capabilities + a shape-valid wiring for a natural-language goal.

        Deterministic and honest: it *suggests* from the registry's shape graph; it
        does not execute. This is the seam the AI Workflow Composer (G1) plugs into.
        """
        words = {w.lower().strip(".,") for w in goal.split()}
        scored = []
        for c in self.registry.all():
            hay = f"{c.name} {c.description} {c.domain} {' '.join(c.tags)}".lower()
            score = sum(1 for w in words if len(w) > 3 and w in hay)
            if score:
                scored.append((score, c))
        scored.sort(key=lambda t: -t[0])
        relevant = [c for _, c in scored[:6]]
        # propose shape-valid edges among the relevant capabilities
        edges = []
        for prod in relevant:
            for op in prod.outputs:
                for cons in relevant:
                    if cons.id == prod.id:
                        continue
                    for ip in cons.inputs:
                        if op.shape == ip.shape:
                            edges.append({
                                "from": f"{prod.id}.{op.name}",
                                "to": f"{cons.id}.{ip.name}",
                                "shape": op.shape.value,
                            })
        return {
            "goal": goal,
            "relevant_capabilities": [
                {"id": c.id, "name": c.name, "status": c.status.value} for c in relevant
            ],
            "candidate_wiring": edges,
            "note": "Suggestion only — the AI Workflow Composer validates, repairs, and runs.",
        }


# --------------------------------------------------------------------------- #
# FastMCP server wrapper
# --------------------------------------------------------------------------- #
def build_server(tools: FactoryTools | None = None):
    """Build a FastMCP server exposing the factory tools. Imports mcp lazily."""
    from mcp.server.fastmcp import FastMCP

    ft = tools or FactoryTools()
    mcp = FastMCP("hcls-ai-factory")

    @mcp.tool()
    def list_capabilities(type: str | None = None, domain: str | None = None,
                          status: str | None = None) -> list[dict]:
        """List the factory's capabilities (engines, agents, models, services), optionally filtered."""
        return ft.list_capabilities(type=type, domain=domain, status=status)

    @mcp.tool()
    def describe_capability(capability_id: str) -> dict:
        """Get the full spec + call contract for one capability."""
        return ft.describe_capability(capability_id)

    @mcp.tool()
    def check_factory_health() -> dict:
        """Ping every live capability with an endpoint and report which are up."""
        return ft.check_factory_health()

    @mcp.tool()
    def invoke_capability(capability_id: str, payload: dict | None = None, path: str = "/") -> dict:
        """Invoke a live capability's endpoint. Refuses to invoke planned/mock capabilities."""
        return ft.invoke_capability(capability_id, payload, path)

    @mcp.tool()
    def plan_pipeline(goal: str) -> dict:
        """Suggest capabilities + a shape-valid wiring for a natural-language goal."""
        return ft.plan_pipeline(goal)

    @mcp.tool()
    def compose_workflow(goal: str) -> dict:
        """Compose a validated, runnable pipeline from a natural-language goal
        (deterministic shape-based wiring + self-repair + a pre-run checklist)."""
        from hcls_common.workflow_composer import WorkflowComposer
        comp = WorkflowComposer(ft.registry, tools=ft)
        pipe, meta = comp.compose(goal)
        return {"pipeline": pipe.to_dict(), **meta}

    return mcp


def main() -> None:  # pragma: no cover - entrypoint
    build_server().run()  # stdio transport by default


if __name__ == "__main__":  # pragma: no cover
    main()
