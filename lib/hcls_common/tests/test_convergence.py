"""Tests for the multi-omics ConvergenceReasoner (§7.2)."""
from hcls_common.convergence import ConvergenceReasoner, ConvergenceSignal
from hcls_common.multiomics import PatientContext
from hcls_common.biokey import BioKeyResolver
from hcls_common.artifact import Maturity


def _tsc_patient():
    # TSC2 lights up across four layers (the flagship convergence)
    return PatientContext(
        "P0",
        genomics={"variants": [{"gene": "TSC2"}]},
        transcriptomics={"driver_genes": ["TSC2"]},
        proteomics_structural={"targets": ["tuberin"]},          # alias -> TSC2 via the resolver
        pathway_activity={"genes": ["TSC2"]},
    )


class TestConvergenceReasoner:
    def test_flagship_convergence_is_a_finding_when_all_layers_live(self):
        sigs = ConvergenceReasoner(resolver=BioKeyResolver()).signals(_tsc_patient())
        assert len(sigs) == 1
        s = sigs[0]
        assert s.entity == "TSC2" and s.breadth == 4
        assert set(s.layers) == {"genomics", "transcriptomics", "proteomics_structural", "pathway_activity"}
        assert s.honesty.maturity is Maturity.live
        assert s.is_finding and s.presentation == "finding"

    def test_honesty_floor_downgrades_to_hypothesis(self):
        # mark the transcriptomics layer hypothesis_only -> the whole convergence is a HYPOTHESIS
        r = ConvergenceReasoner(resolver=BioKeyResolver(),
                                layer_maturity={"transcriptomics": Maturity.hypothesis_only})
        s = r.signals(_tsc_patient())[0]
        assert s.honesty.maturity is Maturity.hypothesis_only    # computed, not copied
        assert not s.is_finding and s.presentation == "hypothesis"
        assert "clinician-review" in s.honesty.requires

    def test_floor_lowers_score(self):
        live = ConvergenceReasoner(resolver=BioKeyResolver()).signals(_tsc_patient())[0].score
        floored = ConvergenceReasoner(
            resolver=BioKeyResolver(),
            layer_maturity={"pathway_activity": Maturity.hypothesis_only}
        ).signals(_tsc_patient())[0].score
        assert floored < live                                    # a hypothesis weighs less

    def test_ranks_broader_convergence_first(self):
        pc = PatientContext(
            "P0",
            genomics={"variants": [{"gene": "TSC2"}, {"gene": "BRCA1"}]},
            transcriptomics={"driver_genes": ["TSC2", "BRCA1"]},
            proteomics_structural={"targets": ["TSC2"]},         # TSC2 in 3 layers, BRCA1 in 2
        )
        sigs = ConvergenceReasoner(resolver=BioKeyResolver()).signals(pc)
        assert [s.entity for s in sigs] == ["TSC2", "BRCA1"]
        assert sigs[0].breadth == 3 and sigs[1].breadth == 2

    def test_lr_hook_combines_when_supplied(self):
        r = ConvergenceReasoner(resolver=BioKeyResolver())
        lrs = {"TSC2": {"genomics": 8.0, "transcriptomics": 4.0,
                        "proteomics_structural": 2.0, "pathway_activity": 3.0}}
        s = r.signals(_tsc_patient(), layer_lrs=lrs)[0]
        # floor_weight(live)=1.0 × product(8*4*2*3)=192
        assert s.score == 192.0

    def test_no_convergence_yields_empty(self):
        pc = PatientContext("P0", genomics={"variants": [{"gene": "TSC2"}]})   # single layer only
        assert ConvergenceReasoner(resolver=BioKeyResolver()).signals(pc) == []
