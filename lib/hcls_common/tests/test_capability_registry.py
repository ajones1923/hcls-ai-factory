"""Tests for the Capability Registry (A1)."""
import pytest

from hcls_common.capability_registry import (
    Capability, CapabilityRegistry, CapabilityType, Port, Serving, Status, ValueShape,
    get_registry,
)


# ── manifest load ───────────────────────────────────────────────────────────
class TestManifest:
    def test_default_manifest_loads(self):
        r = get_registry(reload=True)
        assert len(r) >= 20, "expected the real factory capabilities to load"

    def test_engines_agents_nims_present(self):
        r = get_registry(reload=True)
        assert any(c.id == "genomics-engine" for c in r.find(type="engine"))
        assert len(r.find(type="agent")) >= 8   # 8 canonical intelligence agents
        assert {"molmim-nim", "diffdock-nim"} <= set(r.ids())

    def test_live_vs_planned_split(self):
        r = get_registry(reload=True)
        live = {c.id for c in r.find(status="live")}
        planned = {c.id for c in r.find(status="planned")}
        assert "genomics-engine" in live
        assert "genmol-nim" in planned               # roadmap item, honestly marked
        assert not (live & planned)


# ── honesty rule ────────────────────────────────────────────────────────────
class TestHonestyRule:
    def test_live_cannot_be_mock_served(self):
        reg = CapabilityRegistry()
        bad = Capability(
            id="x", type=CapabilityType.MODEL, name="X", description="",
            serving=Serving.MOCK, status=Status.LIVE,
        )
        with pytest.raises(ValueError, match="honesty"):
            reg.register(bad)

    def test_duplicate_id_rejected(self):
        reg = CapabilityRegistry()
        c = Capability(id="dup", type=CapabilityType.SERVICE, name="D", description="")
        reg.register(c)
        with pytest.raises(ValueError, match="duplicate"):
            reg.register(c)


# ── shape graph (composer foundation) ───────────────────────────────────────
class TestShapeWiring:
    def test_vcf_flows_genomics_to_intelligence(self):
        r = get_registry(reload=True)
        # genomics emits a VCF file; the intelligence engine consumes a VCF file
        assert r.can_connect("genomics-engine", "vcf", "precision-intelligence-engine", "vcf")

    def test_structure_flows_esmfold_to_diffdock(self):
        r = get_registry(reload=True)
        # planned ESMFold emits a structure; DiffDock consumes a protein structure
        assert r.can_connect("esmfold-model", "structure", "diffdock-nim", "protein_structure")

    def test_incompatible_shapes_do_not_connect(self):
        r = get_registry(reload=True)
        # a scalar query cannot feed a file input
        assert not r.can_connect("molmim-nim", "molecules", "genomics-engine", "fastq")

    def test_producers_and_consumers_by_shape(self):
        r = get_registry(reload=True)
        assert any(c.id == "genomics-engine" for c in r.producers_of(ValueShape.FILE))
        assert any(c.id == "diffdock-nim" for c in r.consumers_of(ValueShape.STRUCTURE))


# ── serialization round-trip ────────────────────────────────────────────────
class TestSerialization:
    def test_round_trip(self):
        r = get_registry(reload=True)
        manifest = r.to_manifest()
        r2 = CapabilityRegistry()
        for entry in manifest["capabilities"]:
            r2.register(Capability.from_dict(entry), overwrite=True)
        assert r2.ids() == r.ids()
