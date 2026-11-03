"""PP4 claims the phenotype is highly specific for TSC, so it has to be.

The curator granted it on `bool(ctx.cohort.get("phenotypes"))` — any recorded phenotype at
all. A patient whose single entry was "Seizure" (HP:0001250) earned a supporting pathogenic
criterion from a finding that is not specific to anything, and seizures are not a TSC
diagnostic criterion in the first place.

PP4 now asks whether the recorded features meet the published TSC clinical diagnostic
criteria: two major features, or one major plus two or more minor (Northrup et al., Pediatr
Neurol 2021;123:50-66).
"""
from src.agents.variant_curator.annotation import (
    TSC_MAJOR_HPO, TSC_MINOR_HPO, assign_criteria, meets_tsc_clinical_criteria,
)


def _p(*hpo_ids):
    return [{"hpo_id": h, "label": h} for h in hpo_ids]


CORTICAL_TUBER = "HP:0009717"
HYPOMELANOTIC = "HP:0009719"
SEIZURE = "HP:0001250"
INFANTILE_SPASMS = "HP:0012469"
AUTISM = "HP:0000717"
INTELLECTUAL_DISABILITY = "HP:0001249"
ENAMEL_PITS = "HP:0006297"
RENAL_CYSTS = "HP:0000107"


class TestNonspecificFeaturesDoNotEarnPP4:
    def test_a_lone_seizure_does_not_meet_the_criteria(self):
        met, basis = meets_tsc_clinical_criteria(_p(SEIZURE))
        assert met is False
        assert "0 major, 0 minor" in basis

    def test_epilepsy_and_tand_features_are_not_diagnostic_criteria(self):
        """Common in TSC, and not features specific to it."""
        met, _ = meets_tsc_clinical_criteria(
            _p(SEIZURE, INFANTILE_SPASMS, AUTISM, INTELLECTUAL_DISABILITY))
        assert met is False

    def test_an_empty_or_missing_list_does_not_meet_the_criteria(self):
        assert meets_tsc_clinical_criteria([])[0] is False
        assert meets_tsc_clinical_criteria(None)[0] is False

    def test_one_major_feature_alone_is_not_enough(self):
        met, basis = meets_tsc_clinical_criteria(_p(CORTICAL_TUBER, SEIZURE))
        assert met is False
        assert "1 major, 0 minor" in basis


class TestTheCriteriaAreApplied:
    def test_two_major_features_meet_the_criteria(self):
        met, basis = meets_tsc_clinical_criteria(_p(CORTICAL_TUBER, HYPOMELANOTIC))
        assert met is True
        assert "cortical tubers" in basis and "hypomelanotic" in basis

    def test_one_major_plus_two_minor_meets_the_criteria(self):
        met, basis = meets_tsc_clinical_criteria(_p(CORTICAL_TUBER, ENAMEL_PITS, RENAL_CYSTS))
        assert met is True
        assert "minor:" in basis

    def test_the_same_feature_twice_is_counted_once(self):
        """Two notes recording cortical tubers is one major feature, not two."""
        assert meets_tsc_clinical_criteria(_p(CORTICAL_TUBER, CORTICAL_TUBER))[0] is False

    def test_major_and_minor_sets_do_not_overlap(self):
        assert not set(TSC_MAJOR_HPO) & set(TSC_MINOR_HPO)


class TestTheCriterionStatesItsBasis:
    def test_pp4_reason_names_the_features_and_the_source(self):
        ann = {"is_null": False, "lof_gene": True, "gene": "TSC2", "consequence": "missense",
               "gnomad_af": 0.0, "known_pathogenic": False, "known_source": None}
        met, basis = meets_tsc_clinical_criteria(_p(CORTICAL_TUBER, HYPOMELANOTIC))
        crits = assign_criteria(ann, mosaic=False, phenotype_consistent=met,
                                phenotype_basis=basis)
        pp4 = next(c for c in crits if c.code == "PP4")
        assert "Northrup 2021" in pp4.rationale
        assert "cortical tubers" in pp4.rationale

    def test_pp4_is_absent_when_the_criteria_are_not_met(self):
        ann = {"is_null": False, "lof_gene": True, "gene": "TSC2", "consequence": "missense",
               "gnomad_af": 0.0, "known_pathogenic": False, "known_source": None}
        met, basis = meets_tsc_clinical_criteria(_p(SEIZURE))
        crits = assign_criteria(ann, mosaic=False, phenotype_consistent=met,
                                phenotype_basis=basis)
        assert all(c.code != "PP4" for c in crits)
