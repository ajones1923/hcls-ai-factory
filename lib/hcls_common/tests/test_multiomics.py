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
