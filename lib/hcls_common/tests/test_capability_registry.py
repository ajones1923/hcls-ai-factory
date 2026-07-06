"""Tests for the Capability Registry (A1)."""
import pytest

from hcls_common.capability_registry import (
    Capability, CapabilityRegistry, CapabilityType, Port, Serving, Status, ValueShape,
    ArtifactShape, get_registry,
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


# ── PF-2: semantic ArtifactShape (additive on top of ValueShape) ─────────────
class TestSemanticShape:
    def _reg(self):
        """Two producers of a MAP: one semantically pgx_diplotypes, one hla_alleles;
        and a consumer that wants pgx_diplotypes."""
        reg = CapabilityRegistry()
        reg.register(Capability(
            id="pgx-producer", type=CapabilityType.STAGE, name="PGx", description="",
            outputs=[Port("out", ValueShape.MAP, semantic=ArtifactShape.PGX_DIPLOTYPES)],
        ))
        reg.register(Capability(
            id="hla-producer", type=CapabilityType.STAGE, name="HLA", description="",
            outputs=[Port("out", ValueShape.MAP, semantic=ArtifactShape.HLA_ALLELES)],
        ))
        reg.register(Capability(
            id="pgx-consumer", type=CapabilityType.AGENT, name="A3", description="",
            inputs=[Port("dip", ValueShape.MAP, semantic=ArtifactShape.PGX_DIPLOTYPES)],
        ))
        return reg

    def test_default_semantic_is_unspecified(self):
        p = Port("x", ValueShape.MAP)
        assert p.semantic is ArtifactShape.UNSPECIFIED

    def test_existing_manifest_ports_default_to_unspecified(self):
        # backward compatibility: the real manifest has no `semantic` key yet
        r = get_registry(reload=True)
        for c in r.all():
            for port in (*c.inputs, *c.outputs):
                assert port.semantic is ArtifactShape.UNSPECIFIED

    def test_semantic_survives_round_trip(self):
        reg = self._reg()
        r2 = CapabilityRegistry()
        for entry in reg.to_manifest()["capabilities"]:
            r2.register(Capability.from_dict(entry), overwrite=True)
        out = next(p for p in r2.get("pgx-producer").outputs if p.name == "out")
        assert out.semantic is ArtifactShape.PGX_DIPLOTYPES

    def test_semantic_producers_and_consumers(self):
        reg = self._reg()
        assert [c.id for c in reg.semantic_producers_of(ArtifactShape.PGX_DIPLOTYPES)] == ["pgx-producer"]
        assert [c.id for c in reg.semantic_consumers_of(ArtifactShape.PGX_DIPLOTYPES)] == ["pgx-consumer"]

    def test_coarse_wiring_unchanged_semantic_is_stricter(self):
        reg = self._reg()
        # coarse (default): both MAP producers connect to the MAP consumer
        assert reg.can_connect("pgx-producer", "out", "pgx-consumer", "dip")
        assert reg.can_connect("hla-producer", "out", "pgx-consumer", "dip")
        # semantic (strict): only the matching artifact connects
        assert reg.can_connect("pgx-producer", "out", "pgx-consumer", "dip", semantic=True)
        assert not reg.can_connect("hla-producer", "out", "pgx-consumer", "dip", semantic=True)

    def test_unspecified_never_semantically_connects(self):
        reg = CapabilityRegistry()
        reg.register(Capability(
            id="p", type=CapabilityType.STAGE, name="P", description="",
            outputs=[Port("o", ValueShape.MAP)],           # UNSPECIFIED
        ))
        reg.register(Capability(
            id="c", type=CapabilityType.AGENT, name="C", description="",
            inputs=[Port("i", ValueShape.MAP)],            # UNSPECIFIED
        ))
        assert reg.can_connect("p", "o", "c", "i")                    # coarse: yes
        assert not reg.can_connect("p", "o", "c", "i", semantic=True) # semantic: no (unspecified)
