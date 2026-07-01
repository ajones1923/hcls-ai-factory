"""
External tool-surface ingestion (A3).

Make the factory a *consumer* of other agents' tools, not only a producer: map an external
MCP server's ``tools/list`` JSON Schema directly into Capability registry rows. Each external
tool becomes a first-class capability — discoverable by the tool-surface, wireable by the
workflow composer, and governed by the same input-validation gate (A2) because the JSON Schema's
enum / minimum / maximum / default populate the A1 parameter contract automatically.
"""
from __future__ import annotations

from typing import Any

from hcls_common.capability_registry import (
    Capability, CapabilityRegistry, CapabilityType, Port, Serving, Status, ValueShape,
)


def json_type_to_shape(prop: dict[str, Any]) -> ValueShape:
    """Map a JSON Schema property to a registry value shape."""
    t = prop.get("type")
    if t == "array":
        items = prop.get("items", {}) or {}
        return ValueShape.LIST_OF_OBJECTS if items.get("type") == "object" else ValueShape.LIST
    if t == "object":
        return ValueShape.MAP
    return ValueShape.SCALAR        # string / number / integer / boolean


def external_tool_to_capability(server_name: str, endpoint: str, tool: dict[str, Any]) -> Capability:
    """Convert one MCP tool descriptor into a Capability (carrying its JSON-Schema param contract)."""
    schema = tool.get("inputSchema") or tool.get("input_schema") or {}
    props = schema.get("properties", {}) or {}
    required = set(schema.get("required", []) or [])
    inputs = [
        Port(
            name=name,
            shape=json_type_to_shape(p),
            description=p.get("description", ""),
            required=name in required,
            enum=p.get("enum"),
            minimum=p.get("minimum"),
            maximum=p.get("maximum"),
            default=p.get("default"),
        )
        for name, p in props.items()
    ]
    name = tool["name"]
    return Capability(
        id=f"ext.{server_name}.{name}",
        type=CapabilityType.SERVICE,
        name=name,
        description=(tool.get("description", "") + f" [external tool · {server_name}]").strip(),
        domain="external",
        inputs=inputs,
        outputs=[Port("result", ValueShape.MAP, required=False)],
        endpoint=endpoint,
        invoke_path=f"/tools/{name}",
        serving=Serving.EXTERNAL,
        status=Status.LIVE,
        tags=["external", server_name],
    )


def ingest_external_tools(registry: CapabilityRegistry, server_name: str, endpoint: str,
                          tools: list[dict[str, Any]]) -> list[str]:
    """Register an external server's tools/list into the registry. Returns the new capability ids."""
    added = []
    for tool in tools:
        cap = external_tool_to_capability(server_name, endpoint, tool)
        registry.register(cap, overwrite=True)
        added.append(cap.id)
    return added
