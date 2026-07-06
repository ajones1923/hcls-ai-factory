"""
AI Workflow Composer (G1) — natural-language goal -> a validated, executable, governed
pipeline over the factory's Capability Registry.

The composer turns a stated goal into a DAG of capability invocations and then:
  * **deterministic wiring** — every edge is checked against the registry's declared value
    *shapes*; incompatible connections are rejected with a reason,
  * **self-repair** — nodes referencing unknown/planned capabilities and type-incompatible
    edges are dropped, with a changelog,
  * **pre-run validation** — a checklist of every blocking issue; the pipeline is not
    runnable until clean,
  * **governed execution** — runs each node through the live tool-surface in topological
    order, passing typed outputs downstream,
  * **AI root-cause** — on failure, a structured cause + a data-vs-system verdict.

The LLM planner is injectable (any client with ``generate_json``); without one, a
deterministic shape-graph backward-chaining planner is used, so the composer works — and
is fully testable — offline with no model spend.
"""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from typing import Any, Callable

from hcls_common.capability_registry import (
    ArtifactShape, Capability, CapabilityRegistry, Status, ValueShape, get_registry,
)
from hcls_common.artifact import (
    Honesty, Maturity, derive_artifact, non_inflation_issues, weakest_maturity,
)


# --------------------------------------------------------------------------- #
# DAG model
# --------------------------------------------------------------------------- #
@dataclass
class NodeInput:
    name: str
    value: Any = None                 # a literal value, OR ...
    from_node: str | None = None      # ... a reference to another node's output
    from_port: str | None = None

    @property
    def is_ref(self) -> bool:
        return self.from_node is not None


@dataclass
class Node:
    id: str
    capability: str
    inputs: list[NodeInput] = field(default_factory=list)


@dataclass
class Pipeline:
    goal: str
    nodes: list[Node] = field(default_factory=list)

    def node(self, nid: str) -> Node | None:
        return next((n for n in self.nodes if n.id == nid), None)

    def to_dict(self) -> dict[str, Any]:
        return {"goal": self.goal, "nodes": [
            {"id": n.id, "capability": n.capability,
             "inputs": [{"name": i.name, **({"from": f"{i.from_node}.{i.from_port}"} if i.is_ref
                                            else {"value": i.value})} for i in n.inputs]}
            for n in self.nodes]}


@dataclass
class Issue:
    severity: str        # "error" | "warning"
    node: str | None
    message: str


# --------------------------------------------------------------------------- #
# Composer
# --------------------------------------------------------------------------- #
class WorkflowComposer:
    def __init__(self, registry: CapabilityRegistry | None = None,
                 llm: Any = None, tools: Any = None) -> None:
        self.registry = registry or get_registry()
        self.llm = llm          # NL->pipeline planner (optional; needs .generate_json)
        self.tools = tools      # execution backend (e.g. FactoryTools; needs .invoke_capability)

    # -- validation (deterministic) ---------------------------------------- #
    def validate(self, pipeline: Pipeline) -> list[Issue]:
        issues: list[Issue] = []
        ids = {n.id for n in pipeline.nodes}
        for n in pipeline.nodes:
            try:
                cap = self.registry.get(n.capability)
            except KeyError:
                issues.append(Issue("error", n.id, f"unknown capability '{n.capability}'"))
                continue
            if cap.status is not Status.LIVE:
                issues.append(Issue("error", n.id, f"'{cap.id}' is {cap.status.value}, not live"))
            in_ports = {p.name: p for p in cap.inputs}
            satisfied = {i.name for i in n.inputs}
            # required inputs present?
            for p in cap.inputs:
                if p.required and p.name not in satisfied:
                    issues.append(Issue("error", n.id, f"required input '{p.name}' ({p.shape.value}) unconnected"))
            # each input valid?
            for i in n.inputs:
                if i.name not in in_ports:
                    issues.append(Issue("warning", n.id, f"input '{i.name}' is not a port of '{cap.id}'"))
                    continue
                if i.is_ref:
                    if i.from_node not in ids:
                        issues.append(Issue("error", n.id, f"input '{i.name}' references missing node '{i.from_node}'"))
                        continue
                    src = pipeline.node(i.from_node)
                    src_cap = self._cap(src.capability)
                    out = next((o for o in (src_cap.outputs if src_cap else []) if o.name == i.from_port), None)
                    if out is None:
                        issues.append(Issue("error", n.id, f"'{i.from_node}' has no output '{i.from_port}'"))
                    elif out.shape != in_ports[i.name].shape:
                        issues.append(Issue("error", n.id,
                                            f"shape mismatch on '{i.name}': {out.shape.value} -> {in_ports[i.name].shape.value}"))
        if self._has_cycle(pipeline):
            issues.append(Issue("error", None, "pipeline has a cycle (must be a DAG)"))
        return issues

    def is_runnable(self, pipeline: Pipeline) -> bool:
        return not any(i.severity == "error" for i in self.validate(pipeline))

    def checklist(self, pipeline: Pipeline) -> dict[str, Any]:
        issues = self.validate(pipeline)
        return {"runnable": not any(i.severity == "error" for i in issues),
                "errors": [{"node": i.node, "message": i.message} for i in issues if i.severity == "error"],
                "warnings": [{"node": i.node, "message": i.message} for i in issues if i.severity == "warning"]}

    # -- self-repair -------------------------------------------------------- #
    def repair(self, pipeline: Pipeline) -> tuple[Pipeline, list[str]]:
        log: list[str] = []
        # 1. drop nodes whose capability is unknown or not live
        kept = []
        for n in pipeline.nodes:
            cap = self._cap(n.capability)
            if cap is None:
                log.append(f"dropped node '{n.id}': unknown capability '{n.capability}'")
            elif cap.status is not Status.LIVE:
                log.append(f"dropped node '{n.id}': '{cap.id}' is {cap.status.value}, not live")
            else:
                kept.append(n)
        kept_ids = {n.id for n in kept}
        # 2. strip edges that are type-incompatible or reference dropped/missing nodes
        for n in kept:
            cap = self._cap(n.capability)
            in_ports = {p.name: p for p in cap.inputs}
            new_inputs = []
            for i in n.inputs:
                if not i.is_ref:
                    new_inputs.append(i); continue
                if i.from_node not in kept_ids:
                    log.append(f"stripped edge {i.from_node}.{i.from_port} -> {n.id}.{i.name}: source removed")
                    continue
                src_cap = self._cap(pipeline.node(i.from_node).capability)
                out = next((o for o in src_cap.outputs if o.name == i.from_port), None)
                if out is None or (i.name in in_ports and out.shape != in_ports[i.name].shape):
                    log.append(f"stripped edge {i.from_node}.{i.from_port} -> {n.id}.{i.name}: shape-incompatible")
                    continue
                new_inputs.append(i)
            n.inputs = new_inputs
        return Pipeline(goal=pipeline.goal, nodes=kept), log

    # -- NL -> pipeline ----------------------------------------------------- #
    def compose(self, goal: str) -> tuple[Pipeline, dict[str, Any]]:
        """Draft -> repair -> validate, iterating until the graph stabilizes."""
        draft = self._llm_compose(goal) if self.llm else self._deterministic_compose(goal)
        history = []
        prev = None
        for _ in range(3):                       # converge: repair until unchanged
            repaired, log = self.repair(draft)
            history.append(log)
            sig = repaired.to_dict()
            if sig == prev:
                break
            prev, draft = sig, repaired
        return draft, {"repair_log": [x for x in history if x], "checklist": self.checklist(draft)}

    def _deterministic_compose(self, goal: str) -> Pipeline:
        """Backward-chain a DAG from the shape graph to satisfy a goal's target capability."""
        target = self._best_match(goal)
        nodes: dict[str, Node] = {}

        def ensure(cap_id: str, depth: int = 0) -> None:
            if cap_id in nodes or depth > 6:
                return
            cap = self.registry.get(cap_id)
            node = Node(id=cap_id, capability=cap_id, inputs=[])
            for ip in cap.inputs:
                if not ip.required:
                    continue
                prod = self._find_producer(ip.shape, exclude={cap_id, *nodes})
                if prod:
                    ensure(prod.id, depth + 1)
                    out_port = next(o.name for o in prod.outputs if o.shape == ip.shape)
                    node.inputs.append(NodeInput(name=ip.name, from_node=prod.id, from_port=out_port))
                else:
                    node.inputs.append(NodeInput(name=ip.name, value=f"<provide {ip.shape.value}>"))
            nodes[cap_id] = node

        ensure(target.id)
        return Pipeline(goal=goal, nodes=list(nodes.values()))

    def _llm_compose(self, goal: str) -> Pipeline:
        catalog = [{"id": c.id, "inputs": [(p.name, p.shape.value) for p in c.inputs],
                    "outputs": [(p.name, p.shape.value) for p in c.outputs]}
                   for c in self.registry.live()]
        prompt = (f"Goal: {goal}\nCapabilities (id, inputs[name,shape], outputs[name,shape]):\n"
                  f"{catalog}\nReturn JSON {{'nodes':[{{'id','capability','inputs':"
                  f"[{{'name','from':'node.port' OR 'value':...}}]}}]}}. Wire outputs to inputs only when shapes match.")
        data = self.llm.generate_json(prompt, system_prompt="You compose validated bioinformatics pipelines.")
        nodes = []
        for nd in data.get("nodes", []):
            inputs = []
            for i in nd.get("inputs", []):
                if "from" in i and "." in str(i["from"]):
                    fn, fp = str(i["from"]).rsplit(".", 1)
                    inputs.append(NodeInput(name=i["name"], from_node=fn, from_port=fp))
                else:
                    inputs.append(NodeInput(name=i["name"], value=i.get("value")))
            nodes.append(Node(id=nd["id"], capability=nd.get("capability", nd["id"]), inputs=inputs))
        return Pipeline(goal=goal, nodes=nodes)

    # -- execution ---------------------------------------------------------- #
    @staticmethod
    def _result_honesty(res: dict[str, Any]) -> Honesty:
        """A node's *declared* output honesty, read from its result if present (a capability may
        return ``honesty`` as a dict, or ``maturity`` + ``labels``/``requires``). Absent → live
        (no extra caution declared); the non-inflation combine then caps it to its inputs."""
        h = res.get("honesty")
        if isinstance(h, dict) and "maturity" in h:
            return Honesty.from_dict(h)
        m = res.get("maturity")
        if isinstance(m, str):
            try:
                return Honesty(Maturity(m), labels=list(res.get("labels", [])),
                               requires=list(res.get("requires", [])))
            except ValueError:
                pass
        return Honesty(Maturity.live)

    def run(self, pipeline: Pipeline, inputs: dict[str, dict] | None = None, *,
            patient_id: str | None = None, run_id: str | None = None,
            governed: bool = False) -> dict[str, Any]:
        """Execute a runnable pipeline in topological order via the tool-surface.

        With ``governed=True`` (PF-3/PF-4 wired in), each node's output is additionally wrapped in
        an ``Artifact``: provenance chained to the upstream nodes it consumed, semantic ``shape`` from
        the capability's output port, honesty combined **non-inflatingly** (a node can never emit a
        claim stronger than the weakest input it derives from — a tool that tries is capped and the
        attempt is flagged), and ``patient_id``/``run_id`` threaded through the whole chain. Default
        ``governed=False`` preserves the exact prior behavior (and every existing test)."""
        if not self.is_runnable(pipeline):
            return {"status": "blocked", "checklist": self.checklist(pipeline)}
        if self.tools is None:
            return {"status": "no_executor", "reason": "no tool-surface configured"}
        order = self._topo(pipeline)
        outputs: dict[str, Any] = {}
        trace = []
        node_artifacts: dict[str, Any] = {}          # nid -> Artifact (governed only)
        honesty_issues: list[dict[str, Any]] = []
        rid = run_id or (uuid.uuid4().hex if governed else "")
        for nid in order:
            n = pipeline.node(nid)
            payload = dict((inputs or {}).get(nid, {}))
            for i in n.inputs:
                payload[i.name] = (outputs.get(i.from_node, {}).get(i.from_port)
                                   if i.is_ref else i.value)
            res = self.tools.invoke_capability(n.capability, payload)
            trace.append({"node": nid, "result_status": res.get("status", "ok")})
            if res.get("status") in ("error", "down", "unavailable", "no_endpoint"):
                return {"status": "failed", "failed_node": nid, "trace": trace,
                        "root_cause": self.root_cause(n, res)}
            outputs[nid] = res.get("result", res)

            if governed:
                cap = self._cap(n.capability)
                in_arts = [node_artifacts[i.from_node] for i in n.inputs
                           if i.is_ref and i.from_node in node_artifacts]
                shape = cap.outputs[0].semantic if (cap and cap.outputs) else ArtifactShape.UNSPECIFIED
                own = self._result_honesty(res)
                # flag (do not silently swallow) a capability that claims more than its inputs allow
                if in_arts:
                    weakest = weakest_maturity([a.honesty.maturity for a in in_arts])
                    if own.maturity.rank < weakest.rank:
                        honesty_issues.append({"node": nid,
                            "issue": f"'{n.capability}' declared '{own.maturity.value}' but is capped "
                                     f"to '{weakest.value}' (weakest input) — non-inflation"})
                body = outputs[nid] if isinstance(outputs[nid], dict) else {"value": outputs[nid]}
                art = derive_artifact(shape, body, producer_id=n.capability, inputs=in_arts,
                                      own_maturity=own.maturity, extra_labels=own.labels,
                                      extra_requires=own.requires,
                                      serving=(cap.serving.value if cap else "native"),
                                      patient_id=patient_id, run_id=rid)
                node_artifacts[nid] = art
                honesty_issues.extend({"node": nid, "issue": m} for m in non_inflation_issues(art, in_arts))

        result: dict[str, Any] = {"status": "succeeded", "trace": trace, "outputs": outputs}
        if governed:
            result["run_id"] = rid
            result["artifacts"] = {nid: a.to_dict() for nid, a in node_artifacts.items()}
            result["honesty_issues"] = honesty_issues
        return result

    def root_cause(self, node: Node, result: dict[str, Any]) -> dict[str, Any]:
        st = result.get("status")
        verdict = "system" if st in ("down", "unavailable", "no_endpoint") else "data"
        fixes = {"down": "the capability's service is not running — start it",
                 "unavailable": "the capability is planned/mock, not live — pick a live one",
                 "no_endpoint": "this is a library capability with no callable endpoint",
                 "error": "the service rejected the inputs — check shapes/values"}
        return {"node": node.id, "capability": node.capability, "status": st,
                "verdict": verdict, "suggested_fix": fixes.get(st, "inspect the service logs"),
                "detail": result.get("reason")}

    # -- helpers ------------------------------------------------------------ #
    def _cap(self, cap_id: str) -> Capability | None:
        try:
            return self.registry.get(cap_id)
        except KeyError:
            return None

    def _best_match(self, goal: str) -> Capability:
        words = {w.lower().strip(".,") for w in goal.split() if len(w) > 3}
        scored = []
        for c in self.registry.live():
            hay = f"{c.name} {c.description} {c.domain} {' '.join(c.tags)}".lower()
            scored.append((sum(1 for w in words if w in hay), c))
        scored.sort(key=lambda t: (-t[0], 0 if t[1].outputs else 1))
        return scored[0][1]

    def _find_producer(self, shape: ValueShape, exclude: set[str]) -> Capability | None:
        for c in self.registry.live():
            if c.id not in exclude and any(o.shape == shape for o in c.outputs):
                return c
        return None

    def _edges(self, pipeline: Pipeline) -> list[tuple[str, str]]:
        return [(i.from_node, n.id) for n in pipeline.nodes for i in n.inputs if i.is_ref]

    def _has_cycle(self, pipeline: Pipeline) -> bool:
        try:
            self._topo(pipeline)
            return False
        except ValueError:
            return True

    def _topo(self, pipeline: Pipeline) -> list[str]:
        ids = [n.id for n in pipeline.nodes]
        deps = {nid: set() for nid in ids}
        for a, b in self._edges(pipeline):
            if a in deps and b in deps:
                deps[b].add(a)
        order, ready = [], [nid for nid in ids if not deps[nid]]
        seen = set(ready)
        while ready:
            nid = ready.pop(0)
            order.append(nid)
            for other in ids:
                if nid in deps[other]:
                    deps[other].discard(nid)
                    if not deps[other] and other not in seen:
                        ready.append(other); seen.add(other)
        if len(order) != len(ids):
            raise ValueError("cycle")
        return order
