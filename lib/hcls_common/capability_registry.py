"""
Capability Registry — the single source of truth for everything the HCLS AI Factory
can do (engines, agents, accelerated inference services, pipeline stages, platform
services), with explicit input/output *shapes*.

This is the keystone for two higher-level features:
  * the assistant tool-surface (MCP) — auto-generates tools from these entries, and
  * the AI workflow composer — uses the declared value *shapes* to wire pipelines
    deterministically (derive extraction paths, insert bridge nodes, or reject
    impossible connections).

Design goals: stdlib-only (no heavy deps), declarative JSON manifest, self-registration,
and honest `status` (live / planned / mock) so nothing is advertised that isn't real.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Any, Iterable


# --------------------------------------------------------------------------- #
# Enumerations
# --------------------------------------------------------------------------- #
class CapabilityType(str, Enum):
    ENGINE = "engine"      # a full pipeline engine (genomics, intelligence, discovery)
    AGENT = "agent"        # an intelligence agent (clinical decision support)
    NIM = "nim"            # an accelerated inference microservice
    MODEL = "model"        # a served model endpoint (folding, ADMET, embeddings, ...)
    STAGE = "stage"        # a discrete pipeline stage
    SERVICE = "service"    # platform/infra service (vector db, monitoring, governance)


class ValueShape(str, Enum):
    """The shape of a port's value — used by the workflow composer to wire nodes."""
    SCALAR = "scalar"                 # str/int/float/bool
    LIST = "list"                     # list of scalars
    LIST_OF_OBJECTS = "list_of_objects"
    MAP = "map"                       # dict/object
    FILE = "file"                     # path/URI to a file (FASTQ/BAM/VCF/CSV/...)
    STRUCTURE = "structure"           # a 3D structure blob (PDB/mmCIF)


class Serving(str, Enum):
    NATIVE = "native"        # host process (FastAPI/Streamlit)
    CONTAINER = "container"  # docker service
    LOCAL_NIM = "local_nim"  # NIM container on this box
    CLOUD_NIM = "cloud_nim"  # hosted accelerated endpoint
    MOCK = "mock"            # simulated — never advertised as real
    NONE = "none"            # library/in-process, no endpoint
    EXTERNAL = "external"    # a third-party tool-surface ingested into the registry (A3)


class Status(str, Enum):
    LIVE = "live"          # real and running
    PLANNED = "planned"    # on the roadmap, not yet real
    MOCK = "mock"          # only a mock exists (must be labeled)


# --------------------------------------------------------------------------- #
# Schema
# --------------------------------------------------------------------------- #
@dataclass
class Port:
    name: str
    shape: ValueShape
    description: str = ""
    required: bool = True
    # A1: JSON-Schema-style parameter contract (optional, drives the input-validation gate)
    enum: list[Any] | None = None        # allowed values; anything else is rejected
    minimum: float | None = None         # numeric lower bound (clamped, logged)
    maximum: float | None = None         # numeric upper bound (clamped, logged)
    default: Any = None                  # applied when the input is absent

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Port":
        return cls(
            name=d["name"],
            shape=ValueShape(d["shape"]),
            description=d.get("description", ""),
            required=d.get("required", True),
            enum=d.get("enum"),
            minimum=d.get("minimum"),
            maximum=d.get("maximum"),
            default=d.get("default"),
        )


@dataclass
class Capability:
    id: str
    type: CapabilityType
    name: str
    description: str
    domain: str = ""                       # genomics / proteins / small_molecules / single_cell / clinical / platform
    inputs: list[Port] = field(default_factory=list)
    outputs: list[Port] = field(default_factory=list)
    endpoint: str | None = None            # url or host:port (None for libraries)
    invoke_path: str = "/"                 # POST path on the endpoint (e.g. /fold, /admet)
    serving: Serving = Serving.NATIVE
    gpu: bool = False
    cost_class: str = "low"                # low / medium / high
    status: Status = Status.LIVE
    tags: list[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Capability":
        return cls(
            id=d["id"],
            type=CapabilityType(d["type"]),
            name=d["name"],
            description=d["description"],
            domain=d.get("domain", ""),
            inputs=[Port.from_dict(p) for p in d.get("inputs", [])],
            outputs=[Port.from_dict(p) for p in d.get("outputs", [])],
            endpoint=d.get("endpoint"),
            invoke_path=d.get("invoke_path", "/"),
            serving=Serving(d.get("serving", "native")),
            gpu=d.get("gpu", False),
            cost_class=d.get("cost_class", "low"),
            status=Status(d.get("status", "live")),
            tags=d.get("tags", []),
        )

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # serialize enums to their string values
        d["type"] = self.type.value
        d["serving"] = self.serving.value
        d["status"] = self.status.value
        d["inputs"] = [{**asdict(p), "shape": p.shape.value} for p in self.inputs]
        d["outputs"] = [{**asdict(p), "shape": p.shape.value} for p in self.outputs]
        return d


# --------------------------------------------------------------------------- #
# A2: input-validation gate (input-side bracket; pairs with the output honesty gate)
# --------------------------------------------------------------------------- #
def validate_inputs(cap: Capability, payload: dict[str, Any] | None) -> tuple[dict[str, Any], list[str]]:
    """Validate a payload against a capability's input contract.

    Returns (cleaned_payload, issues). Issues prefixed ERROR (blocking) or WARN (auto-fixed):
      * required input missing (and no default)      -> ERROR
      * value not in the declared enum               -> ERROR (reject)
      * numeric value outside [minimum, maximum]      -> WARN (clamp + log)
      * absent input with a default                   -> default applied
    """
    out = dict(payload or {})
    issues: list[str] = []
    ports = {p.name: p for p in cap.inputs}
    for p in cap.inputs:
        if p.name not in out:
            if p.default is not None:
                out[p.name] = p.default
            elif p.required:
                issues.append(f"ERROR: required input '{p.name}' ({p.shape.value}) missing")
    for name, val in list(out.items()):
        p = ports.get(name)
        if p is None:
            continue
        if p.enum is not None and val not in p.enum:
            issues.append(f"ERROR: input '{name}'={val!r} not in allowed values {p.enum}")
        if isinstance(val, (int, float)) and not isinstance(val, bool):
            if p.minimum is not None and val < p.minimum:
                out[name] = p.minimum
                issues.append(f"WARN: input '{name}' {val} clamped up to minimum {p.minimum}")
            if p.maximum is not None and val > p.maximum:
                out[name] = p.maximum
                issues.append(f"WARN: input '{name}' {val} clamped down to maximum {p.maximum}")
    return out, issues


def inputs_ok(issues: list[str]) -> bool:
    return not any(i.startswith("ERROR") for i in issues)


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
DEFAULT_MANIFEST = Path(__file__).with_name("capabilities.json")


class CapabilityRegistry:
    """In-memory registry of Capabilities, loadable from a JSON manifest."""

    def __init__(self) -> None:
        self._caps: dict[str, Capability] = {}

    # -- population --------------------------------------------------------- #
    def register(self, cap: Capability, *, overwrite: bool = False) -> None:
        if cap.id in self._caps and not overwrite:
            raise ValueError(f"duplicate capability id: {cap.id!r}")
        self._validate(cap)
        self._caps[cap.id] = cap

    def load_manifest(self, path: str | Path | None = None) -> "CapabilityRegistry":
        path = Path(path) if path else DEFAULT_MANIFEST
        data = json.loads(Path(path).read_text())
        for entry in data.get("capabilities", []):
            self.register(Capability.from_dict(entry), overwrite=True)
        return self

    @staticmethod
    def _validate(cap: Capability) -> None:
        if not cap.id or not cap.name:
            raise ValueError("capability requires id and name")
        for port in (*cap.inputs, *cap.outputs):
            if not isinstance(port.shape, ValueShape):
                raise ValueError(f"{cap.id}: port {port.name} has invalid shape")
        # honesty rule: a 'live' capability must not be served by a mock
        if cap.status is Status.LIVE and cap.serving is Serving.MOCK:
            raise ValueError(f"{cap.id}: a LIVE capability cannot be MOCK-served (honesty rule)")

    # -- query -------------------------------------------------------------- #
    def get(self, cap_id: str) -> Capability:
        return self._caps[cap_id]

    def ids(self) -> list[str]:
        return sorted(self._caps)

    def all(self) -> list[Capability]:
        return [self._caps[i] for i in self.ids()]

    def find(
        self,
        *,
        type: CapabilityType | str | None = None,
        domain: str | None = None,
        status: Status | str | None = None,
        gpu: bool | None = None,
    ) -> list[Capability]:
        t = CapabilityType(type) if isinstance(type, str) else type
        s = Status(status) if isinstance(status, str) else status
        out = []
        for c in self.all():
            if t and c.type != t:
                continue
            if domain and c.domain != domain:
                continue
            if s and c.status != s:
                continue
            if gpu is not None and c.gpu != gpu:
                continue
            out.append(c)
        return out

    def live(self) -> list[Capability]:
        return self.find(status=Status.LIVE)

    # -- shape graph helpers (used by the composer) ------------------------- #
    def producers_of(self, shape: ValueShape) -> list[Capability]:
        return [c for c in self.all() if any(p.shape == shape for p in c.outputs)]

    def consumers_of(self, shape: ValueShape) -> list[Capability]:
        return [c for c in self.all() if any(p.shape == shape for p in c.inputs)]

    def can_connect(self, producer_id: str, output: str, consumer_id: str, input: str) -> bool:
        """True iff the named output port shape matches the named input port shape."""
        op = next(p for p in self.get(producer_id).outputs if p.name == output)
        ip = next(p for p in self.get(consumer_id).inputs if p.name == input)
        return op.shape == ip.shape

    # -- serialization ------------------------------------------------------ #
    def to_manifest(self) -> dict[str, Any]:
        return {"capabilities": [c.to_dict() for c in self.all()]}

    def __len__(self) -> int:
        return len(self._caps)


# --------------------------------------------------------------------------- #
# Singleton accessor
# --------------------------------------------------------------------------- #
_REGISTRY: CapabilityRegistry | None = None


def get_registry(reload: bool = False) -> CapabilityRegistry:
    """Process-wide registry, loaded from the default manifest on first use."""
    global _REGISTRY
    if _REGISTRY is None or reload:
        _REGISTRY = CapabilityRegistry()
        if DEFAULT_MANIFEST.exists():
            _REGISTRY.load_manifest()
    return _REGISTRY
