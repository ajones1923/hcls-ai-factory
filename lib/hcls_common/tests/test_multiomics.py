"""F1: multi-omics per-patient join + cross-omics convergence."""
from hcls_common.multiomics import MultiOmicsStore, MultiOmicsRecord


def _patient(store, pid="P1"):
    store.add_genomics(pid, {"variants": [{"gene": "TSC2", "clinical_significance": "Pathogenic"},
                                          {"gene": "BRCA1", "clinical_significance": "Benign"}],
                             "secondary_findings": [{"gene": "TSC2", "acmg_sf_condition": "Tuberous sclerosis"}]})
    store.add_transcriptomics(pid, {"cell_types": ["T cell", "Monocyte"], "marker_genes": ["TSC2", "CD3D"]})
    store.add_proteomics(pid, {"structures": [{"name": "tsc2_model"}], "targets": ["TSC2"]})


class TestJoin:
    def test_cross_omics_convergence(self):
        s = MultiOmicsStore(); _patient(s)
        rec = s.get("P1")
        links = rec.cross_omics_links()
        tsc2 = next(l for l in links if l["gene"] == "TSC2")
        assert tsc2["n_layers"] == 3                          # variant + marker + target all converge
        assert set(tsc2["layers"]) == {"genomics", "transcriptomics", "proteomics"}
        # BRCA1 only in genomics -> not a convergent link
        assert all(l["gene"] != "BRCA1" for l in links)

    def test_summary(self):
        s = MultiOmicsStore(); _patient(s)
        summ = s.get("P1").summary()
        assert summ["n_variants"] == 2 and summ["n_secondary_findings"] == 1
        assert summ["cell_types"] == ["T cell", "Monocyte"]
        assert summ["cross_omics_links"][0]["gene"] == "TSC2"   # ranked by breadth

    def test_empty_patient(self):
        assert MultiOmicsStore().get("none") is None
        assert MultiOmicsRecord("X").cross_omics_links() == []


# ── PF-1: PatientContext (the 13-layer longitudinal memory) ──────────────────
from hcls_common.multiomics import PatientContext, PatientContextStore  # noqa: E402


class TestPatientContext:
    def test_thirteen_layers_and_ledger(self):
        pc = PatientContext("P0")
        assert len(pc.LAYERS) == 13
        assert pc.artifacts == [] and pc.populated_layers() == []

    def test_same_gene_converges_across_omics(self):
        pc = PatientContext(
            "P0",
            genomics={"variants": [{"gene": "TSC2"}]},
            transcriptomics={"driver_genes": ["TSC2"], "marker_genes": ["CD3D"]},  # driver, not marker
            proteomics_structural={"targets": ["TSC2"]},
            pathway_activity={"genes": ["TSC2"]},
        )
        links = pc.cross_omics_links()
        tsc2 = next(l for l in links if l["gene"] == "TSC2")
        assert set(tsc2["layers"]) == {"genomics", "transcriptomics", "proteomics_structural", "pathway_activity"}
        assert tsc2["n_layers"] == 4

    def test_marker_genes_alone_do_not_converge(self):
        # the E8 fix: cell-type markers must not fire convergence, only DE driver genes
        pc = PatientContext(
            "P0",
            genomics={"variants": [{"gene": "CD3D"}]},
            transcriptomics={"marker_genes": ["CD3D"]},   # marker only, no driver_genes
        )
        assert pc.cross_omics_links() == []

    def test_add_artifact_dedups(self):
        pc = PatientContext("P0")
        pc.add_artifact("a1"); pc.add_artifact("a1"); pc.add_artifact("a2")
        assert pc.artifacts == ["a1", "a2"]

    def test_round_trip(self):
        pc = PatientContext("P0", genomics={"variants": [{"gene": "TSC2"}]}, artifacts=["a1"])
        assert PatientContext.from_dict(pc.to_dict()) == pc

    def test_bridge_from_legacy_record(self):
        store = MultiOmicsStore(); _patient(store, "P1")
        pc = PatientContext.from_multiomics_record(store.get("P1"))
        assert pc.patient_id == "P1"
        assert pc.genomics["variants"][0]["gene"] == "TSC2"
        assert pc.proteomics_structural["targets"] == ["TSC2"]


class TestPatientContextStore:
    def test_layer_writers_and_thread(self):
        s = PatientContextStore()
        s.add_genomics("P0", {"variants": [{"gene": "TSC2"}]})
        s.add_pathway_activity("P0", {"genes": ["TSC2"], "mTOR": "hyperactive"})
        s.add_therapeutics("P0", [{"drug": "everolimus", "line": 1}])
        s.record_artifact("P0", "art-1")
        pc = s.get("P0")
        assert pc.clinical["prior_therapy"][0]["drug"] == "everolimus"
        assert pc.pathway_activity["mTOR"] == "hyperactive"
        assert pc.artifacts == ["art-1"]
        assert {"genomics", "pathway_activity", "clinical"} <= set(pc.populated_layers())

    def test_unknown_layer_rejected(self):
        import pytest
        with pytest.raises(ValueError, match="unknown patient-context layer"):
            PatientContextStore().update_layer("P0", "not_a_layer", {})
