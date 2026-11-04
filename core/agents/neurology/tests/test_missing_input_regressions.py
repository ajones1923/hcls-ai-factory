"""Where an input is absent, the answer must name it -- not substitute a benign value.

One root cause produced most of this agent's clinical defects: `inputs.get(key, <default>)`
with a default that happened to be the best or most reassuring value, and the default was then
printed as a finding. ALSFRS-R substituted 4 = normal function; NIHSS substituted 0 = no
deficit; HIT-6 padded with "never"; KPS became 100; tumour size became 0.0 cm; MGMT became
unmethylated; `progression` became "chronic"; `cognitive_status` became "normal"; NEDA-3
started True.

Each class below pins one of those, with the executed reproduction in its docstring.
"""
import asyncio
import types

import pytest

from src.clinical_workflows import (
    BrainTumorGradingWorkflow,
    DementiaEvaluationWorkflow,
    MSMonitoringWorkflow,
    NeuromuscularEvaluationWorkflow,
    ParkinsonsAssessmentWorkflow,
)
from src.models import SeverityLevel


class _FakeRequest:
    def __init__(self, llm=None):
        self.app = types.SimpleNamespace(
            state=types.SimpleNamespace(llm_client=llm, metrics=None, metrics_lock=None)
        )


def _findings(result) -> str:
    return " || ".join(result.findings)


def _recommendations(result) -> str:
    return " || ".join(result.recommendations)


# =====================================================================
# N7 -- "NEDA-3 ACHIEVED" on empty input
# =====================================================================

class TestN7NEDA3IsNotAssertedFromSilence:
    """`MSMonitoringWorkflow().run({})` yielded "EDSS: 0.0", "NEDA-3 STATUS: ACHIEVED -- no
    disease activity detected" and severity LOW: a patient with no MRI, no relapse history and
    no recorded EDSS was certified disease-free. `relapses_last_year`, `new_t2_lesions` and
    `gad_enhancing` defaulted to 0, `edss_current` to 0.0, and `neda3` started True.
    `edss_prior is None` also skipped the progression limb in silence while NEDA-3 was still
    asserted.
    """

    def test_an_empty_request_is_not_neda3_achieved(self):
        result = MSMonitoringWorkflow().run({})
        assert "NEDA-3 STATUS: ACHIEVED" not in _findings(result)
        assert "NEDA-3 STATUS: NOT ASSESSABLE" in _findings(result)

    def test_an_empty_request_is_not_severity_low(self):
        result = MSMonitoringWorkflow().run({})
        assert result.severity != SeverityLevel.LOW

    def test_an_unrecorded_edss_is_not_printed_as_0(self):
        result = MSMonitoringWorkflow().run({})
        assert "EDSS: 0.0" not in _findings(result)
        assert "EDSS: NOT RECORDED" in _findings(result)

    def test_the_missing_limbs_are_named(self):
        result = MSMonitoringWorkflow().run({})
        text = _findings(result)
        assert "relapse count" in text
        assert "MRI lesion activity" in text
        assert "EDSS progression" in text

    def test_a_genuinely_quiet_year_is_still_neda3(self):
        result = MSMonitoringWorkflow().run({
            "relapses_last_year": 0, "new_t2_lesions": 0, "gad_enhancing": 0,
            "edss_current": 2.0, "edss_prior": 2.0,
        })
        assert "NEDA-3 STATUS: ACHIEVED" in _findings(result)
        assert result.severity == SeverityLevel.LOW

    def test_genuine_activity_is_still_not_met(self):
        result = MSMonitoringWorkflow().run({
            "relapses_last_year": 1, "new_t2_lesions": 0, "gad_enhancing": 0,
            "edss_current": 2.0, "edss_prior": 2.0,
        })
        assert "NEDA-3 STATUS: NOT MET" in _findings(result)


class TestN7EDSSProgressionIsAChangeNotALevel:
    """`/ms/assess` read `if edss_score is not None and edss_score > 3.0: neda_3 = False;
    edss_progression = True`, and the endpoint had no prior-EDSS field at all. Wrong in both
    directions: a patient stable at EDSS 4.0 for ten years with no relapses and no new lesions
    got `no_edss_progression: false`, `treatment_response: "suboptimal"` and a DMT-escalation
    suggestion; while a patient who went 1.0 -> 3.0 (genuine progression) got
    `edss_progression: false` and `neda_3_achieved: true`.

    The >= 1.0 (baseline <= 5.5) / >= 0.5 (baseline >= 6.0) rule was already coded in
    clinical_workflows.py and is now applied here too.
    """

    def _assess(self, **kwargs):
        from api.routes.neuro_clinical import MSAssessRequest, ms_assess

        return asyncio.run(ms_assess(MSAssessRequest(**kwargs), _FakeRequest()))

    def test_a_patient_stable_at_edss_4_is_not_reported_as_progressing(self):
        response = self._assess(
            edss_score=4.0, edss_prior=4.0,
            relapse_count_2yr=0, new_t2_lesions=0, gad_enhancing_lesions=0,
        )
        assert response.neda_status["no_edss_progression"] is True
        assert response.neda_status["neda_3_achieved"] is True
        assert response.dmt_assessment["treatment_response"] == "adequate"

    def test_a_patient_who_went_1_to_3_is_reported_as_progressing(self):
        response = self._assess(
            edss_score=3.0, edss_prior=1.0,
            relapse_count_2yr=0, new_t2_lesions=0, gad_enhancing_lesions=0,
        )
        assert response.neda_status["no_edss_progression"] is False
        assert response.neda_status["neda_3_achieved"] is False

    def test_a_half_point_change_at_a_high_baseline_counts(self):
        response = self._assess(
            edss_score=6.5, edss_prior=6.0,
            relapse_count_2yr=0, new_t2_lesions=0, gad_enhancing_lesions=0,
        )
        assert response.neda_status["no_edss_progression"] is False

    def test_a_half_point_change_at_a_low_baseline_does_not(self):
        response = self._assess(
            edss_score=2.5, edss_prior=2.0,
            relapse_count_2yr=0, new_t2_lesions=0, gad_enhancing_lesions=0,
        )
        assert response.neda_status["no_edss_progression"] is True

    def test_neda3_is_not_assessable_with_no_prior_edss(self):
        response = self._assess(edss_score=4.0, relapse_count_2yr=0,
                                new_t2_lesions=0, gad_enhancing_lesions=0)
        assert response.neda_status["neda_3_achieved"] is None
        assert response.neda_status["neda_3_assessable"] is False
        assert "edss_prior" in response.neda_status["note"]

    def test_an_empty_request_is_not_neda3_achieved(self):
        response = self._assess()
        assert response.neda_status["neda_3_achieved"] is not True
        assert response.dmt_assessment["treatment_response"] == "not_assessable"


# =====================================================================
# N11 -- ATN printed unmeasured biomarkers as negative
# =====================================================================

class TestN11ATNAndAntiAmyloid:
    """`{amyloid_status: "positive", age: 72}` yielded "ATN biomarker profile: A+T-N-",
    "A+T- = Alzheimer's pathological change (amyloid-positive, pre-tangle stage)" -- a staging
    conclusion from a tau test that was never run -- and "ANTI-AMYLOID THERAPY: Patient may be
    eligible for lecanemab (Leqembi) or donanemab (Kisunla)", with "[INPUT WARNING] No
    cognitive test scores provided" in the same findings list.

    `eligible` started True and was only cleared by a cognitive score < 18, so with no MoCA
    and no MMSE the check was skipped entirely.
    """

    REPRODUCTION = {"amyloid_status": "positive", "age": 72}

    def test_unmeasured_analytes_are_not_printed_as_negative(self):
        result = DementiaEvaluationWorkflow().run(self.REPRODUCTION)
        assert "A+T-N-" not in _findings(result)
        assert "A+T?N?" in _findings(result)

    def test_the_unmeasured_analytes_are_named(self):
        result = DementiaEvaluationWorkflow().run(self.REPRODUCTION)
        assert "not measured: tau, neurodegeneration" in _findings(result)

    def test_no_staging_conclusion_is_drawn_from_an_unrun_test(self):
        result = DementiaEvaluationWorkflow().run(self.REPRODUCTION)
        assert "pre-tangle stage" not in _findings(result)
        assert "tau NOT MEASURED" in _findings(result)

    def test_anti_amyloid_therapy_is_not_offered_without_cognitive_staging(self):
        result = DementiaEvaluationWorkflow().run(self.REPRODUCTION)
        assert "Patient may be eligible for" not in _findings(result)
        assert "eligibility INDETERMINATE" in _findings(result)

    def test_the_inputs_the_safety_monitoring_turns_on_are_named(self):
        result = DementiaEvaluationWorkflow().run(self.REPRODUCTION)
        text = _findings(result)
        assert "cognitive stage" in text
        assert "anticoagulation status" in text
        assert "intracerebral haemorrhage" in text
        assert "microhaemorrhage" in text

    def test_a_measured_negative_tau_still_stages_as_a_plus_t_minus(self):
        result = DementiaEvaluationWorkflow().run({
            **self.REPRODUCTION, "tau_status": "negative",
            "neurodegeneration_status": "negative",
        })
        assert "A+T-N-" in _findings(result)
        assert "pre-tangle stage" in _findings(result)

    def test_a_fully_specified_eligible_case_is_still_offered_therapy(self):
        result = DementiaEvaluationWorkflow().run({
            "amyloid_status": "positive", "tau_status": "positive", "age": 72,
            "mmse_score": 25, "on_anticoagulation": False,
            "prior_ich": False, "microhemorrhage_count": 0,
        })
        assert "ANTI-AMYLOID THERAPY: Patient may be eligible" in _findings(result)

    def test_anticoagulation_blocks_anti_amyloid_therapy(self):
        result = DementiaEvaluationWorkflow().run({
            "amyloid_status": "positive", "tau_status": "positive", "age": 72,
            "mmse_score": 25, "on_anticoagulation": True,
            "prior_ich": False, "microhemorrhage_count": 0,
        })
        assert "Patient may be eligible" not in _findings(result)
        assert "On anticoagulation" in _findings(result)

    def test_an_incomplete_moca_is_not_a_cognitive_stage(self):
        result = DementiaEvaluationWorkflow().run({
            "amyloid_status": "positive", "tau_status": "positive", "age": 72,
            "moca_domain_scores": {"orientation": 6}, "education_years": 16,
            "on_anticoagulation": False, "prior_ich": False, "microhemorrhage_count": 0,
        })
        assert "eligibility INDETERMINATE" in _findings(result)


class TestN11CholinesteraseInhibitorForMCIContradictedItsOwnCitation:
    """`/dementia/evaluate` recommended a cholinesterase inhibitor for any MoCA < 26 -- i.e.
    for MCI -- while citing the AAN 2018 MCI Practice Guideline, whose MCI recommendation is
    that clinicians *may choose not to* prescribe them (Level B) and must discuss the absent
    evidence first (Level A).
    """

    def _evaluate(self, **kwargs):
        from api.routes.neuro_clinical import DementiaEvaluateRequest, dementia_evaluate

        return asyncio.run(
            dementia_evaluate(DementiaEvaluateRequest(**kwargs), _FakeRequest())
        )

    def test_an_mci_score_no_longer_gets_a_bare_cholinesterase_recommendation(self):
        response = self._evaluate(moca_total=22)
        joined = " ".join(response.treatment_options)
        assert "donepezil" not in joined
        assert "not approved for MCI" in joined

    def test_the_guideline_position_is_stated(self):
        response = self._evaluate(moca_total=22)
        joined = " ".join(response.treatment_options)
        assert "may choose not to prescribe" in joined
        assert "AAN 2018" in joined

    def test_a_dementia_range_score_still_gets_the_recommendation(self):
        response = self._evaluate(moca_total=14)
        joined = " ".join(response.treatment_options)
        assert "donepezil" in joined
        assert "memantine" in joined

    def test_education_years_can_no_longer_inflate_the_total(self):
        """The endpoint never passed its declared `education_years` to the calculator; it
        relied on the client smuggling it inside `moca_items`, which is N1."""
        response = self._evaluate(
            moca_items={
                "visuospatial": 2, "naming": 2, "attention": 3, "language": 1,
                "abstraction": 0, "delayed_recall": 0, "orientation": 4,
            },
            education_years=16,
        )
        assert response.moca_score == 12
        assert "Normal" not in (response.moca_interpretation or "")


# =====================================================================
# N12 -- KPS 100, size 0.0 cm, MGMT unmethylated, all printed as findings
# =====================================================================

class TestN12TumorDefaultsArePrintedAsFindings:
    """A GBM case with neither supplied reported "Location: not specified, Size: 0.0 cm" and
    "KPS: 100". The fabricated 100 also kept the `kps < 60` palliative-care recommendation
    from ever firing. For meningioma, `if size < 3.0 and not enhancing` gave an unmeasured
    tumour "Observation with serial MRI (6 months, then annually) for small, asymptomatic
    meningioma".
    """

    def test_an_unrecorded_kps_is_not_printed_as_100(self):
        result = BrainTumorGradingWorkflow().run({"histology": "glioblastoma"})
        assert "KPS: 100" not in _findings(result)
        assert "KPS: NOT RECORDED" in _findings(result)

    def test_an_unmeasured_tumour_is_not_printed_as_0_cm(self):
        result = BrainTumorGradingWorkflow().run({"histology": "glioblastoma"})
        assert "Size: 0.0 cm" not in _findings(result)
        assert "Size: NOT MEASURED" in _findings(result)

    def test_an_unmeasured_meningioma_is_not_routed_to_observation(self):
        result = BrainTumorGradingWorkflow().run({"histology": "meningioma"})
        assert "Observation with serial MRI" not in _recommendations(result)
        assert "was not measured" in _recommendations(result)

    def test_a_measured_small_meningioma_is_still_observed(self):
        result = BrainTumorGradingWorkflow().run({
            "histology": "meningioma", "size_cm": 1.5, "enhancing": False,
        })
        assert "Observation with serial MRI" in _recommendations(result)

    def test_an_untested_mgmt_is_not_reported_as_unmethylated(self):
        result = BrainTumorGradingWorkflow().run({"histology": "glioblastoma"})
        assert "UNMETHYLATED" not in _findings(result)
        assert "MGMT promoter methylation: NOT TESTED" in _findings(result)

    def test_a_low_kps_still_triggers_palliative_care(self):
        result = BrainTumorGradingWorkflow().run({
            "histology": "glioblastoma", "kps": 50, "size_cm": 4.0,
        })
        assert "palliative" in _recommendations(result).lower()

    def test_the_endpoint_does_not_report_untested_mgmt_as_unmethylated(self):
        from api.routes.neuro_clinical import TumorGradeRequest, tumor_grade

        response = asyncio.run(
            tumor_grade(
                TumorGradeRequest(histology="glioblastoma", molecular_markers=[]),
                _FakeRequest(),
            )
        )
        joined = " ".join(response.treatment_recommendations)
        assert "MGMT unmethylated -- consider clinical trial alternatives" not in joined
        assert "NOT TESTED" in joined

    def test_the_endpoint_still_reports_a_measured_unmethylated_result(self):
        from api.routes.neuro_clinical import TumorGradeRequest, tumor_grade

        response = asyncio.run(
            tumor_grade(
                TumorGradeRequest(
                    histology="glioblastoma", molecular_markers=["mgmt_unmethylated"]
                ),
                _FakeRequest(),
            )
        )
        joined = " ".join(response.treatment_recommendations)
        assert "MGMT unmethylated" in joined

    def test_an_unrecorded_kps_is_flagged_on_the_endpoint_too(self):
        from api.routes.neuro_clinical import TumorGradeRequest, tumor_grade

        response = asyncio.run(
            tumor_grade(TumorGradeRequest(histology="glioblastoma"), _FakeRequest())
        )
        assert any("not recorded" in r for r in response.recommendations)


# =====================================================================
# N14 -- progression defaulted to "chronic": acute GBS reported as CIDP
# =====================================================================

class TestN14AcuteGBSIsNotReportedAsCIDP:
    """`{emg_ncs: {conduction_block: true, f_wave_prolonged: true}, weakness_pattern:
    "generalized"}` with `progression` absent reported "Progression: chronic" -- a value
    nobody supplied -- a CIDP/MMN/anti-MAG/CMT1 differential, severity MODERATE, and an
    antibody panel plus a CMT gene panel. No airway monitoring for a demyelinating
    polyneuropathy of unknown tempo. The GBS limb carries the only urgent actions in the file:
    "URGENT: Monitor respiratory function (FVC q4h) -- intubate if FVC < 20 mL/kg",
    IVIg/plasmapheresis, and severity CRITICAL.
    """

    REPRODUCTION = {
        "emg_ncs": {"conduction_block": True, "f_wave_prolonged": True},
        "weakness_pattern": "generalized",
    }

    def test_an_unrecorded_tempo_is_not_printed_as_chronic(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        assert "Progression: chronic" not in _findings(result)
        assert "Progression: NOT RECORDED" in _findings(result)

    def test_gbs_stays_on_the_differential(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        assert "Guillain-Barre" in _findings(result)
        assert "CANNOT BE EXCLUDED" in _findings(result)

    def test_respiratory_monitoring_is_recommended(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        text = _recommendations(result)
        assert "FVC q4h" in text
        assert "20 mL/kg" in text

    def test_the_severity_is_escalated(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        assert result.severity == SeverityLevel.CRITICAL

    def test_the_tempo_is_asked_for(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        assert "Establish the tempo" in _recommendations(result)

    def test_no_cmt_gene_panel_is_ordered_on_an_unknown_tempo(self):
        result = NeuromuscularEvaluationWorkflow().run(self.REPRODUCTION)
        assert "CMT gene panel" not in _recommendations(result)

    def test_a_recorded_chronic_tempo_still_gets_the_cidp_pathway(self):
        result = NeuromuscularEvaluationWorkflow().run({
            **self.REPRODUCTION, "progression": "chronic",
        })
        assert "CIDP" in _findings(result)
        assert "FVC q4h" not in _recommendations(result)

    def test_a_recorded_acute_tempo_still_gets_the_gbs_pathway(self):
        result = NeuromuscularEvaluationWorkflow().run({
            **self.REPRODUCTION, "progression": "acute",
        })
        assert result.severity == SeverityLevel.CRITICAL
        assert "IVIg" in _recommendations(result)


class TestN14EndpointEscalatesOnAcuteNotOnlyRapid:
    """`/neuromuscular/evaluate` escalated only when `progression_rate` contained the literal
    "rapid", so a patient described as "acute" got no urgent recommendation; and
    `"GBS (if acute)"` matched none of the treatment branches, so a GBS differential yielded
    no treatment at all.
    """

    def _evaluate(self, **kwargs):
        from api.routes.neuro_clinical import (
            NeuromuscularEvaluateRequest,
            neuromuscular_evaluate,
        )

        base = dict(
            presentation="ascending weakness over four days",
            emg_findings="demyelinating with conduction block",
        )
        return asyncio.run(
            neuromuscular_evaluate(
                NeuromuscularEvaluateRequest(**{**base, **kwargs}), _FakeRequest()
            )
        )

    def test_an_acute_tempo_escalates(self):
        response = self._evaluate(progression_rate="acute")
        assert any("Urgent evaluation" in r for r in response.recommendations)

    def test_a_gbs_differential_now_has_a_treatment(self):
        response = self._evaluate(progression_rate="acute")
        names = [d["diagnosis"] for d in response.differential_diagnosis]
        assert "Guillain-Barre syndrome" in names
        assert "GBS (if acute)" not in names
        joined = " ".join(response.treatment_options)
        assert "IVIg" in joined
        assert "plasmapheresis" in joined

    def test_steroids_are_called_out_as_ineffective_in_gbs(self):
        response = self._evaluate(progression_rate="acute")
        assert "NOT effective in GBS" in " ".join(response.treatment_options)

    def test_respiratory_monitoring_is_the_first_recommendation(self):
        response = self._evaluate(progression_rate="acute")
        assert "URGENT" in response.recommendations[0]
        assert "FVC q4h" in response.recommendations[0]

    def test_an_unrecorded_tempo_keeps_gbs_and_its_monitoring(self):
        response = self._evaluate()
        names = [d["diagnosis"] for d in response.differential_diagnosis]
        assert "Guillain-Barre syndrome" in names
        assert any("FVC q4h" in r for r in response.recommendations)
        assert any("NOT RECORDED" in r for r in response.recommendations)

    def test_a_chronic_tempo_does_not_raise_a_gbs_alarm(self):
        response = self._evaluate(progression_rate="chronic, over two years")
        names = [d["diagnosis"] for d in response.differential_diagnosis]
        assert "Guillain-Barre syndrome" not in names
        assert "CIDP" in names


# =====================================================================
# N13 -- DBS cleared on an unexamined cognitive domain; PIGD by substring
# =====================================================================

class TestN13DBSAndMotorSubtype:
    """`cog_status = inputs.get("cognitive_status", "normal")` and `dbs_candidate = True`
    initially, so a patient with motor fluctuations, dyskinesias, 8-year duration, excellent
    levodopa response and UPDRS-III 66, with `cognitive_status` absent, yielded "DBS CANDIDATE:
    ... no cognitive contraindication" -- an assertion about a domain never assessed. A MoCA,
    when supplied, was computed and printed but never used in the decision.
    """

    CANDIDATE = {
        "motor_fluctuations": True, "dyskinesias": True,
        "symptom_duration_years": 8.0, "levodopa_response": "excellent",
        "age": 62,
    }

    def _updrs(self, per_item=2):
        from src.clinical_scales import UPDRSCalculator

        return {item: per_item for item in UPDRSCalculator.ITEMS}

    def test_dbs_is_not_cleared_on_an_unassessed_cognitive_domain(self):
        result = ParkinsonsAssessmentWorkflow().run({
            **self.CANDIDATE, "updrs_scores": self._updrs(),
        })
        assert "DBS CANDIDATE" not in _findings(result)
        assert "Cognitive status NOT ASSESSED" in _findings(result)

    def test_the_moca_that_was_supplied_is_actually_used(self):
        result = ParkinsonsAssessmentWorkflow().run({
            **self.CANDIDATE, "updrs_scores": self._updrs(),
            "moca_domain_scores": {
                "visuospatial": 2, "naming": 2, "attention": 3, "language": 1,
                "abstraction": 0, "delayed_recall": 0, "orientation": 4,
            },
            "education_years": 16,
        })
        assert "DBS CANDIDATE" not in _findings(result)
        assert "below the cognitive threshold" in _findings(result)

    def test_a_recorded_normal_cognition_still_clears_the_domain(self):
        result = ParkinsonsAssessmentWorkflow().run({
            **self.CANDIDATE, "updrs_scores": self._updrs(),
            "cognitive_status": "normal",
        })
        assert "DBS CANDIDATE" in _findings(result)

    def test_dementia_still_contraindicates(self):
        result = ParkinsonsAssessmentWorkflow().run({
            **self.CANDIDATE, "updrs_scores": self._updrs(),
            "cognitive_status": "dementia",
        })
        assert "DBS CANDIDATE" not in _findings(result)
        assert "dementia" in _findings(result).lower()

    def test_unrecorded_demographics_are_not_printed_as_zero(self):
        result = ParkinsonsAssessmentWorkflow().run({})
        assert "Age: 0, onset age: 0, duration: 0.0 years" not in _findings(result)
        assert "not recorded" in _findings(result)

    def test_postural_tremor_no_longer_counts_towards_pigd(self):
        """`"postural" in k.lower()` also matched 3.15a/3.15b postural TREMOR, so a
        pure-tremor patient contributed tremor points to the PIGD sum."""
        from api.routes.neuro_clinical import ParkinsonsAssessRequest, parkinsons_assess

        response = asyncio.run(parkinsons_assess(
            ParkinsonsAssessRequest(updrs_items={
                "3.15a_postural_tremor_r": 4,
                "3.17a_rest_tremor_rue": 4,
                "3.10_gait": 0,
                "3.12_postural_stability": 0,
            }),
            _FakeRequest(),
        ))
        assert response.motor_subtype == "tremor_dominant"
        assert response.motor_subtype != "postural_instability_gait_difficulty"

    def test_a_genuine_pigd_pattern_is_still_classified_as_pigd(self):
        from api.routes.neuro_clinical import ParkinsonsAssessRequest, parkinsons_assess

        response = asyncio.run(parkinsons_assess(
            ParkinsonsAssessRequest(updrs_items={
                "3.15a_postural_tremor_r": 0,
                "3.17a_rest_tremor_rue": 0,
                "3.10_gait": 3,
                "3.12_postural_stability": 3,
            }),
            _FakeRequest(),
        ))
        assert response.motor_subtype == "postural_instability_gait_difficulty"

    def test_subtyping_is_refused_when_one_item_group_was_never_scored(self):
        from api.routes.neuro_clinical import ParkinsonsAssessRequest, parkinsons_assess

        response = asyncio.run(parkinsons_assess(
            ParkinsonsAssessRequest(updrs_items={"3.10_gait": 3}),
            _FakeRequest(),
        ))
        assert "not_assessable" in response.motor_subtype


# =====================================================================
# N15 -- the MOH 10-day threshold collapsed to 15 days
# =====================================================================

class TestN15MedicationOveruseThresholdIsDrugClassSpecific:
    """A patient on sumatriptan 12 days/month for 6 months MEETS ICHD-3 8.2.2, but was
    reported as merely "approaching" the threshold, and the withdrawal-protocol recommendation
    (gated on "High" in moh_risk) never fired.

    ICHD-3, the issuing body's own criteria
    (https://ichd-3.org/8-headache-attributed-to-a-substance-or-its-withdrawal/
     8-2-medication-overuse-headache-moh/, sections 8.2.2, 8.2.4, 8.2.6):
    simple analgesics / paracetamol / NSAIDs >= 15 days/month for > 3 months; triptans,
    ergotamine, opioids >= 10 days/month; any combination of these classes >= 10 days/month.
    """

    def _classify(self, **kwargs):
        from api.routes.neuro_clinical import HeadacheClassifyRequest, headache_classify

        base = dict(headache_description="unilateral throbbing headaches")
        return asyncio.run(headache_classify(
            HeadacheClassifyRequest(**{**base, **kwargs}), _FakeRequest()
        ))

    def test_sumatriptan_12_days_meets_the_criteria(self):
        response = self._classify(
            analgesic_days_per_month=12, current_medications=["sumatriptan 100mg"]
        )
        assert "High risk" in response.medication_overuse_risk
        assert "10 days/month" in response.medication_overuse_risk

    def test_the_withdrawal_protocol_now_fires(self):
        response = self._classify(
            analgesic_days_per_month=12, current_medications=["sumatriptan 100mg"]
        )
        assert any("withdrawal" in r.lower() for r in response.recommendations)

    def test_ibuprofen_12_days_does_not_meet_the_criteria(self):
        response = self._classify(
            analgesic_days_per_month=12, current_medications=["ibuprofen 400mg"]
        )
        assert "High risk" not in response.medication_overuse_risk
        assert "15 days/month" in response.medication_overuse_risk

    def test_ibuprofen_16_days_does_meet_the_criteria(self):
        response = self._classify(
            analgesic_days_per_month=16, current_medications=["ibuprofen 400mg"]
        )
        assert "High risk" in response.medication_overuse_risk

    @pytest.mark.parametrize("drug", ["rizatriptan", "codeine", "dihydroergotamine",
                                      "tramadol", "Excedrin"])
    def test_every_ten_day_class_uses_the_ten_day_threshold(self, drug):
        response = self._classify(analgesic_days_per_month=11, current_medications=[drug])
        assert "High risk" in response.medication_overuse_risk

    def test_the_three_month_duration_criterion_is_mentioned(self):
        response = self._classify(
            analgesic_days_per_month=12, current_medications=["sumatriptan"]
        )
        assert "> 3 months" in response.medication_overuse_risk

    def test_an_unknown_drug_class_does_not_assume_a_threshold(self):
        response = self._classify(analgesic_days_per_month=12)
        assert "Not determinable" in response.medication_overuse_risk
        assert "class is required" in response.medication_overuse_risk


# =====================================================================
# N19 -- Dravet contraindication fired only from genetics
# =====================================================================

class TestN19DravetContraindicationDoesNotWaitOnGenetics:
    """`syndrome = "dravet"` was reachable only through a pathogenic SCN1A result, and the
    contraindication scan is gated on `syndrome in _CONTRAINDICATED_ASMS`. A clinically
    diagnosed Dravet patient taking carbamazepine, with no genetic testing in the payload, got
    no "DANGER: ... CONTRAINDICATED" finding and no severity escalation. Sodium-channel
    blockers worsen Dravet seizures. The `/epilepsy/classify` endpoint had no contraindication
    warning at all -- it only replaced the ASM list.
    """

    def _classify(self, **kwargs):
        from api.routes.neuro_clinical import EpilepsyClassifyRequest, epilepsy_classify

        base = dict(seizure_description="prolonged generalized convulsions with fever")
        return asyncio.run(epilepsy_classify(
            EpilepsyClassifyRequest(**{**base, **kwargs}), _FakeRequest()
        ))

    def test_a_clinical_diagnosis_triggers_the_contraindication_warning(self):
        response = self._classify(
            clinical_diagnosis="Dravet syndrome", current_aeds=["carbamazepine 400mg BID"]
        )
        assert response.contraindication_warnings
        joined = " ".join(response.contraindication_warnings)
        assert "DANGER" in joined
        assert "carbamazepine" in joined
        assert "CONTRAINDICATED" in joined

    def test_the_warning_reaches_the_recommendation_list(self):
        response = self._classify(
            clinical_diagnosis="Dravet syndrome", current_aeds=["carbamazepine"]
        )
        assert any("DANGER" in r for r in response.recommendations)

    def test_the_syndrome_specific_asm_list_is_still_replaced(self):
        response = self._classify(clinical_diagnosis="Dravet syndrome")
        assert any("Stiripentol" in a for a in response.aed_recommendations)

    def test_a_genetic_result_still_works(self):
        response = self._classify(
            genetic_results="SCN1A pathogenic variant", current_aeds=["phenytoin"]
        )
        assert any("DANGER" in w for w in response.contraindication_warnings)

    def test_the_workflow_warns_from_a_clinical_diagnosis_with_no_genetics(self):
        """This is the reproduction: the contraindication scan is gated on
        `syndrome in _CONTRAINDICATED_ASMS`, and `syndrome` was only ever set by SCN1A."""
        from src.clinical_workflows import EpilepsyFocusWorkflow

        result = EpilepsyFocusWorkflow().run({
            "clinical_syndrome_diagnosis": "Dravet syndrome",
            "current_asms": ["carbamazepine 400mg BID"],
        })
        text = _findings(result)
        assert "DANGER" in text
        assert "Carbamazepine is CONTRAINDICATED in dravet" in text
        assert result.severity == SeverityLevel.CRITICAL

    def test_the_workflow_issues_no_such_warning_without_the_diagnosis(self):
        """Guards the test above: it must be the diagnosis doing the work."""
        from src.clinical_workflows import EpilepsyFocusWorkflow

        result = EpilepsyFocusWorkflow().run({
            "current_asms": ["carbamazepine 400mg BID"],
        })
        assert "DANGER" not in _findings(result)


class TestN19DrugResistanceRequiresFailureNotConcurrency:
    """`surgical_candidacy` fired on `len(current_aeds) >= 2 and seizure_frequency`, so a
    patient on two concurrent ASMs was told "Drug-resistant epilepsy criteria may be met"
    under a cited drug-resistant-epilepsy guideline. The ILAE definition, and this workflow's
    own `num_failed_asms >= 2`, require FAILURE of two adequate trials.
    """

    def _classify(self, **kwargs):
        from api.routes.neuro_clinical import EpilepsyClassifyRequest, epilepsy_classify

        base = dict(seizure_description="focal seizures with impaired awareness")
        return asyncio.run(epilepsy_classify(
            EpilepsyClassifyRequest(**{**base, **kwargs}), _FakeRequest()
        ))

    def test_two_concurrent_asms_do_not_establish_drug_resistance(self):
        response = self._classify(
            current_aeds=["levetiracetam", "lamotrigine"], seizure_frequency="2/month"
        )
        assert "may be met" not in response.surgical_candidacy
        assert "Drug-resistance not established" in response.surgical_candidacy

    def test_two_failed_trials_do_establish_it(self):
        response = self._classify(
            current_aeds=["lacosamide"],
            failed_aeds=["levetiracetam", "lamotrigine"],
            seizure_frequency="2/month",
        )
        assert "may be met" in response.surgical_candidacy
        assert "2 ASMs failed" in response.surgical_candidacy


# =====================================================================
# The endpoints that surface a refused scale rather than a reassuring score
# =====================================================================

class TestRefusedScalesReachTheHTTPSurface:
    """A calculator that refuses to impute must not become a 500 or a silent best case."""

    def test_scale_calculate_returns_422_naming_the_problem(self):
        from fastapi import HTTPException

        from api.routes.neuro_clinical import ScaleCalculateRequest, scale_calculate

        with pytest.raises(HTTPException) as exc:
            asyncio.run(scale_calculate(
                ScaleCalculateRequest(scale_name="hit6", items={"a": 13, "b": 13}),
                _FakeRequest(),
            ))
        assert exc.value.status_code == 422
        assert "exactly 6" in exc.value.detail

    def test_scale_calculate_still_returns_400_for_an_unknown_scale(self):
        from fastapi import HTTPException

        from api.routes.neuro_clinical import ScaleCalculateRequest, scale_calculate

        with pytest.raises(HTTPException) as exc:
            asyncio.run(scale_calculate(
                ScaleCalculateRequest(scale_name="not_a_scale", items={}), _FakeRequest()
            ))
        assert exc.value.status_code == 400

    def test_scale_calculate_reports_completeness_on_the_wire(self):
        from api.routes.neuro_clinical import ScaleCalculateRequest, scale_calculate

        response = asyncio.run(scale_calculate(
            ScaleCalculateRequest(
                scale_name="alsfrs", items={"walking": 1, "climbing_stairs": 0}
            ),
            _FakeRequest(),
        ))
        assert response.data_complete is False
        assert response.items_scored == 2
        assert response.items_expected == 12
        assert response.data_warnings

    def test_headache_classify_refuses_a_short_hit6_rather_than_scoring_it(self):
        from fastapi import HTTPException

        from api.routes.neuro_clinical import HeadacheClassifyRequest, headache_classify

        with pytest.raises(HTTPException) as exc:
            asyncio.run(headache_classify(
                HeadacheClassifyRequest(
                    headache_description="unilateral throbbing",
                    hit6_items={"pain_severity": 13, "daily_activities": 13},
                ),
                _FakeRequest(),
            ))
        assert exc.value.status_code == 422
        assert "HIT-6" in exc.value.detail

    def test_neuromuscular_reports_an_incomplete_alsfrs(self):
        from api.routes.neuro_clinical import (
            NeuromuscularEvaluateRequest,
            neuromuscular_evaluate,
        )

        response = asyncio.run(neuromuscular_evaluate(
            NeuromuscularEvaluateRequest(
                presentation="progressive limb weakness",
                alsfrs_items={"walking": 1, "climbing_stairs": 0},
            ),
            _FakeRequest(),
        ))
        assert response.alsfrs_score == 1
        assert response.alsfrs_score != 41
        assert any("incomplete" in w.lower() for w in response.data_warnings)

    def test_no_bulbar_claim_is_made_when_no_bulbar_item_was_submitted(self):
        from api.routes.neuro_clinical import (
            NeuromuscularEvaluateRequest,
            neuromuscular_evaluate,
        )

        response = asyncio.run(neuromuscular_evaluate(
            NeuromuscularEvaluateRequest(
                presentation="progressive limb weakness",
                alsfrs_items={"walking": 1, "climbing_stairs": 0},
            ),
            _FakeRequest(),
        ))
        assert not any(
            "Prominent bulbar involvement" in r
            for r in response.recommended_workup + response.recommendations
        )
