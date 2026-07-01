"""E6: ACMG SF secondary-findings categorization (deterministic join, no model)."""
from acmg_sf import is_on_panel, is_reportable, secondary_findings, panel_summary, ACMG_SF_PANEL


class TestPanel:
    def test_membership(self):
        assert is_on_panel("BRCA1") and is_on_panel("brca1")      # case-insensitive
        assert not is_on_panel("MADEUPGENE")

    def test_reportable_requires_pathogenic_and_panel(self):
        assert is_reportable("BRCA1", "Pathogenic")
        assert is_reportable("MYBPC3", "Likely_pathogenic")
        assert not is_reportable("BRCA1", "Benign")               # on panel, benign -> no
        assert not is_reportable("EGFR", "Pathogenic")            # pathogenic, off panel -> no
        assert not is_reportable("BRCA1", None)

    def test_secondary_findings_filters_and_enriches(self):
        variants = [
            {"gene": "BRCA1", "clinical_significance": "Pathogenic", "pos": 100},
            {"gene": "BRCA1", "clinical_significance": "Benign", "pos": 200},
            {"gene": "TTN", "clinical_significance": "Pathogenic", "pos": 300},   # off panel
            {"gene": "TSC2", "clinical_significance": "Likely pathogenic", "pos": 400},
        ]
        sf = secondary_findings(variants)
        assert len(sf) == 2
        assert {v["gene"] for v in sf} == {"BRCA1", "TSC2"}
        assert all("acmg_sf_condition" in v for v in sf)
        assert sf[0]["acmg_sf_condition"]  # condition populated

    def test_summary(self):
        s = panel_summary()
        assert s["n_genes"] == len(ACMG_SF_PANEL) and s["n_conditions"] > 10
