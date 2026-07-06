"""Tests for the BioKey resolver (PF-5)."""
from hcls_common.biokey import BioKey, BioKeyResolver, EntityKind, resolve
from hcls_common.multiomics import PatientContext


class TestResolver:
    def test_seed_aliases_fold_to_canonical(self):
        r = BioKeyResolver()
        assert r.resolve_gene("N-myc") == "MYCN"
        assert r.resolve_gene("nmyc") == "MYCN"
        assert r.resolve_gene("tuberin") == "TSC2"
        assert r.resolve_gene("hamartin") == "TSC1"

    def test_unknown_gene_normalizes_to_upper(self):
        r = BioKeyResolver()
        assert r.resolve_gene("brca1") == "BRCA1"          # no alias → upper-cased passthrough
        assert r.resolve("  tp53 ").id == "TP53"

    def test_register_and_load_aliases(self):
        r = BioKeyResolver(gene_aliases={})                # empty seed
        assert r.resolve_gene("N-myc") == "N-MYC"          # not folded yet
        r.register_alias("N-myc", "MYCN")
        assert r.resolve_gene("N-myc") == "MYCN"
        r.load_aliases({"ERBB2": ["HER2", "NEU"]})
        assert r.resolve_gene("her2") == "ERBB2"

    def test_typed_biokey_and_verbatim_ids(self):
        r = BioKeyResolver()
        assert r.resolve("HLA-A*02:01", EntityKind.HLA) == BioKey(EntityKind.HLA, "HLA-A*02:01")
        # HPO / VRS ids are already canonical → kept verbatim (only trimmed)
        assert r.resolve(" HP:0001250 ", EntityKind.HPO).id == "HP:0001250"
        assert str(r.resolve("TSC2")) == "gene:TSC2"

    def test_module_level_resolve_uses_default(self):
        assert resolve("tuberin").id == "TSC2"


class TestConvergenceWithResolver:
    def test_aliases_converge_only_with_resolver(self):
        # TSC2 as a variant, but 'tuberin' as the structural target — same biology, different string
        pc = PatientContext(
            "P0",
            genomics={"variants": [{"gene": "TSC2"}]},
            proteomics_structural={"targets": ["tuberin"]},
        )
        # without the resolver: raw strings don't meet → no convergence
        assert pc.cross_omics_links() == []
        # with the resolver: tuberin folds to TSC2 → one convergent node across 2 layers
        links = pc.cross_omics_links(resolver=BioKeyResolver())
        assert len(links) == 1
        assert links[0]["gene"] == "TSC2" and links[0]["n_layers"] == 2

    def test_default_behavior_unchanged_without_resolver(self):
        pc = PatientContext(
            "P0",
            genomics={"variants": [{"gene": "TSC2"}]},
            transcriptomics={"driver_genes": ["TSC2"]},
        )
        assert pc.cross_omics_links()[0]["n_layers"] == 2   # identical to PF-1 behavior


# ── NF-9: ontology-version pinning + reconciliation ──────────────────────────
class TestOntologyVersionPinning:
    def test_pin_and_read_versions(self):
        r = BioKeyResolver(ontology_versions={"hpo": "2025-01"})
        r.pin_version("HGNC", "2025-03")
        assert r.versions() == {"HPO": "2025-01", "HGNC": "2025-03"}

    def test_reconcile_flags_version_mismatch(self):
        e1 = BioKeyResolver(ontology_versions={"HPO": "2024-10", "HGNC": "2025-01"})
        a7 = BioKeyResolver(ontology_versions={"HPO": "2025-01", "HGNC": "2025-01"})
        mism = e1.reconcile(a7)
        assert mism == ["ontology HPO: 2024-10 vs 2025-01"]      # HGNC agrees, HPO differs

    def test_reconcile_clean_when_versions_agree(self):
        a = BioKeyResolver(ontology_versions={"HPO": "2025-01"})
        assert a.reconcile({"HPO": "2025-01"}) == []
