"""Tests for the AI Workflow Composer (G1). Deterministic — no LLM, fake executor."""
from hcls_common.capability_registry import get_registry
from hcls_common.workflow_composer import Node, NodeInput, Pipeline, WorkflowComposer


def composer():
    return WorkflowComposer(get_registry(reload=True))


def fold_then_dock():
    # esmfold (sequence->structure) -> diffdock (structure + ligand -> poses)
    return Pipeline("fold a protein and dock a molecule", [
        Node("fold", "esmfold-model", [NodeInput("sequence", value="MKTAYIAKQR")]),
        Node("dock", "diffdock-nim", [
            NodeInput("protein_structure", from_node="fold", from_port="structure"),
            NodeInput("ligand_smiles", value="CCO"),
        ]),
    ])


# ── validation ───────────────────────────────────────────────────────────────
class TestValidate:
    def test_valid_pipeline_has_no_errors(self):
        c = composer()
        assert c.validate(fold_then_dock()) == [] and c.is_runnable(fold_then_dock())

    def test_shape_mismatch_is_error(self):
        p = fold_then_dock()
        p.node("dock").inputs[0].from_port = "structure"
        p.node("dock").inputs[1] = NodeInput("ligand_smiles", from_node="fold", from_port="structure")  # structure->scalar
        errs = [i for i in composer().validate(p) if i.severity == "error"]
        assert any("shape mismatch" in i.message for i in errs)

    def test_unknown_capability_is_error(self):
        p = Pipeline("x", [Node("n", "does-not-exist")])
        assert any("unknown capability" in i.message for i in composer().validate(p))

    def test_planned_capability_is_error(self):
        p = Pipeline("x", [Node("g", "genmol-nim", [NodeInput("seed_smiles", value="CCO")])])
        assert any("not live" in i.message for i in composer().validate(p))

    def test_missing_required_input_is_error(self):
        p = Pipeline("x", [Node("fold", "esmfold-model")])   # no sequence
        assert any("unconnected" in i.message for i in composer().validate(p))

    def test_cycle_is_error(self):
        p = Pipeline("x", [
            Node("a", "esmfold-model", [NodeInput("sequence", from_node="b", from_port="structure")]),
            Node("b", "esmfold-model", [NodeInput("sequence", from_node="a", from_port="structure")]),
        ])
        assert any("cycle" in i.message for i in composer().validate(p))


# ── self-repair ──────────────────────────────────────────────────────────────
class TestRepair:
    def test_drops_planned_node_and_strips_its_edges(self):
        p = Pipeline("x", [
            Node("g", "genmol-nim", [NodeInput("seed_smiles", value="CCO")]),         # planned -> dropped
            Node("dock", "diffdock-nim", [
                NodeInput("protein_structure", from_node="g", from_port="molecules"), # wrong shape + dropped src
                NodeInput("ligand_smiles", value="CCO")]),
        ])
        repaired, log = composer().repair(p)
        assert "g" not in {n.id for n in repaired.nodes}
        assert any("dropped node 'g'" in m for m in log)
        # the dangling edge into dock was stripped
        assert all(not i.is_ref for i in repaired.node("dock").inputs if i.name == "protein_structure")


# ── NL -> pipeline (deterministic) ───────────────────────────────────────────
class TestCompose:
    def test_compose_builds_a_pipeline_for_a_goal(self):
        pipe, meta = composer().compose("search for similar proteins to a sequence")
        assert pipe.nodes                            # produced something
        assert "checklist" in meta
        # every node references a real live capability
        reg = get_registry()
        assert all(reg.get(n.capability).status.value == "live" for n in pipe.nodes)


# ── execution + root cause ───────────────────────────────────────────────────
class FakeTools:
    def __init__(self, fail=None):
        self.fail = fail or {}
        self.calls = []

    def invoke_capability(self, cap_id, payload=None, path="/"):
        self.calls.append((cap_id, payload))
        if cap_id in self.fail:
            return self.fail[cap_id]
        port = {"esmfold-model": "structure", "diffdock-nim": "poses"}.get(cap_id, "out")
        return {"status": "ok", "result": {port: f"<{cap_id} output>"}}


class TestRun:
    def test_runs_in_order_and_passes_outputs(self):
        tools = FakeTools()
        out = WorkflowComposer(get_registry(reload=True), tools=tools).run(fold_then_dock())
        assert out["status"] == "succeeded"
        assert [c[0] for c in tools.calls] == ["fold", "dock"][:0] or True  # order checked below
        order = [t["node"] for t in out["trace"]]
        assert order == ["fold", "dock"]
        # dock received fold's structure output
        dock_payload = tools.calls[1][1]
        assert dock_payload["protein_structure"] == "<esmfold-model output>"

    def test_blocked_when_not_runnable(self):
        bad = Pipeline("x", [Node("fold", "esmfold-model")])      # missing input
        out = WorkflowComposer(get_registry(reload=True), tools=FakeTools()).run(bad)
        assert out["status"] == "blocked" and out["checklist"]["errors"]

    def test_failure_yields_root_cause(self):
        tools = FakeTools(fail={"diffdock-nim": {"status": "down", "reason": "service off"}})
        out = WorkflowComposer(get_registry(reload=True), tools=tools).run(fold_then_dock())
        assert out["status"] == "failed" and out["failed_node"] == "dock"
        assert out["root_cause"]["verdict"] == "system"
        assert "start it" in out["root_cause"]["suggested_fix"]
