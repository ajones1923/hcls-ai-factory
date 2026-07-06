"""Tests for the composed-chain lineage / reproducibility manifest (§8)."""
from hcls_common.artifact import new_artifact, derive_artifact, Maturity
from hcls_common.capability_registry import ArtifactShape
from hcls_common.lineage import chain_lineage


def _chain(dock_maturity=Maturity.live):
    fold = new_artifact(ArtifactShape.PROTEIN_STRUCTURE, {"s": 1}, producer_id="esmfold-model",
                        serving="native", run_id="R", patient_id="P0")
    dock = derive_artifact(ArtifactShape.DOCKING_POSES, {"p": 1}, producer_id="diffdock-nim",
                           inputs=[fold], own_maturity=dock_maturity, serving="cloud_nim",
                           run_id="R", patient_id="P0")
    return {"fold": fold, "dock": dock}


class TestChainLineage:
    def test_nodes_edges_and_serving_transparency(self):
        m = chain_lineage(_chain())
        assert m["n_hops"] == 2 and m["run_id"] == "R"
        ids = {n["producer"]: n["id"] for n in m["nodes"]}
        assert {n["producer"] for n in m["nodes"]} == {"esmfold-model", "diffdock-nim"}
        # the lineage edge dock <- fold
        assert {"from": ids["esmfold-model"], "to": ids["diffdock-nim"]} in m["edges"]
        # where each hop ACTUALLY ran is recorded (native vs cloud_nim — burst transparency)
        serving = {n["producer"]: n["serving"] for n in m["nodes"]}
        assert serving == {"esmfold-model": "native", "diffdock-nim": "cloud_nim"}

    def test_chain_honesty_is_the_floor(self):
        # dock declared hypothesis_only -> the whole chain's honesty is hypothesis_only
        m = chain_lineage(_chain(dock_maturity=Maturity.hypothesis_only))
        assert m["chain_honesty"]["maturity"] == "hypothesis_only"

    def test_replay_hash_stable_across_reruns_but_structure_sensitive(self):
        h1 = chain_lineage(_chain())["lineage_hash"]
        h2 = chain_lineage(_chain())["lineage_hash"]        # fresh uuids/timestamps
        assert h1 == h2                                     # same structure -> same hash
        single = {"only": new_artifact(ArtifactShape.VCF_VARIANTS, {}, producer_id="genomics-engine")}
        assert chain_lineage(single)["lineage_hash"] != h1  # different chain -> different hash

    def test_accepts_to_dict_forms(self):
        as_dicts = {k: v.to_dict() for k, v in _chain().items()}
        assert chain_lineage(as_dicts)["n_hops"] == 2

    def test_flows_from_governed_composer_run(self):
        from hcls_common.workflow_composer import WorkflowComposer, Pipeline, Node, NodeInput
        from hcls_common import get_registry

        class FakeTools:
            def invoke_capability(self, cap_id, payload=None, path="/"):
                port = {"esmfold-model": "structure", "diffdock-nim": "poses"}.get(cap_id, "out")
                return {"status": "ok", "result": {port: f"<{cap_id}>"}}

        p = Pipeline("x", [
            Node("fold", "esmfold-model", [NodeInput("sequence", value="MKT")]),
            Node("dock", "diffdock-nim", [
                NodeInput("protein_structure", from_node="fold", from_port="structure"),
                NodeInput("ligand_smiles", value="CCO")])])
        out = WorkflowComposer(get_registry(reload=True), tools=FakeTools()).run(
            p, governed=True, patient_id="P0")
        m = chain_lineage(out["artifacts"])
        assert m["n_hops"] == 2 and m["run_id"] == out["run_id"] and len(m["edges"]) == 1
