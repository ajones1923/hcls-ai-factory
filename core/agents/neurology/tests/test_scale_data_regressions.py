"""Reference cases for the scale calculators, one per confirmed defect.

Each test here asserts that a specific wrong output is gone. The agent's existing
`test_clinical_scales.py` builds `ScaleResult` models by hand and asserts the number it just
typed, so it passes identically with all ten calculators deleted -- none of these defects
could have been caught by it.

`src.clinical_scales.calculate_scale` is the single shared implementation that the route-local
`_calculate_*` functions in `api/routes/neuro_clinical.py` now delegate to, so exercising it
exercises the path the shipped Streamlit UI actually calls.
"""
import pytest

from src.clinical_scales import (
    ALSFRSCalculator,
    ASPECTSCalculator,
    EDSSCalculator,
    HIT6Calculator,
    MoCACalculator,
    NIHSSCalculator,
    ScaleDataError,
    UPDRSCalculator,
    calculate_scale,
)


# =====================================================================
# N1 -- MoCA inflated by years of education, on the live UI path
# =====================================================================

class TestMoCAEducationIsNotAScore:
    """A true 12/30 exam reported 28/30 "Normal cognition" with zero recommendations.

    `_calculate_moca` summed every numeric value in `items`, and the UI puts
    `education_years` (default 16) in that same dict two lines above the +1 adjustment that
    reads it. The whole workup block is gated on `total < 26`, so the response carried no
    neuropsychological testing, no MRI and no B12/TSH/RPR screen for reversible causes: a
    moderate dementia reported as a normal cognitive screen.
    """

    # visuospatial 2, naming 2, attention 3, language 1, abstraction 0, recall 0,
    # orientation 4 = a true MoCA of 12/30.
    TRUE_12 = {
        "visuospatial": 2, "naming": 2, "attention": 3, "language": 1,
        "abstraction": 0, "delayed_recall": 0, "orientation": 4,
    }

    def test_education_16_no_longer_inflates_the_total(self):
        result = calculate_scale("moca", {**self.TRUE_12, "education_years": 16})
        assert result.score == 12
        assert result.score != 28

    def test_a_true_12_is_not_reported_as_normal_cognition(self):
        result = calculate_scale("moca", {**self.TRUE_12, "education_years": 16})
        assert "normal" not in result.severity_category.lower()
        assert result.severity_category == "Moderate cognitive impairment"

    def test_the_workup_recommendations_fire(self):
        """They were gated on total < 26 and so were empty at the inflated 28."""
        result = calculate_scale("moca", {**self.TRUE_12, "education_years": 16})
        joined = " ".join(result.recommendations).lower()
        assert result.recommendations
        assert "neuropsychological" in joined
        assert "mri" in joined
        assert "b12" in joined

    def test_education_12_adds_exactly_one_point_once(self):
        result = calculate_scale("moca", {**self.TRUE_12, "education_years": 12})
        assert result.score == 13

    def test_education_13_adds_nothing(self):
        result = calculate_scale("moca", {**self.TRUE_12, "education_years": 13})
        assert result.score == 12

    def test_the_adjustment_cannot_push_a_perfect_exam_past_30(self):
        perfect = dict(MoCACalculator.DOMAINS)
        result = MoCACalculator.calculate(perfect, education_years=8)
        assert result.score == 30

    def test_unknown_education_applies_no_adjustment_and_says_so(self):
        result = MoCACalculator.calculate(self.TRUE_12, education_years=None)
        assert result.score == 12
        assert any("not provided" in w for w in result.warnings)

    def test_the_total_cannot_exceed_30_from_stray_keys(self):
        result = calculate_scale(
            "moca", {**self.TRUE_12, "education_years": 16, "mmse": 29, "age": 71}
        )
        assert result.score == 12
        assert any("unrecognized item 'mmse'" in w for w in result.warnings)

    def test_an_incomplete_moca_is_not_reported_as_normal(self):
        """Unrecorded domains can only move the total down, so 'normal' is unsafe."""
        result = calculate_scale("moca", {"orientation": 6, "naming": 3})
        assert result.complete is False
        assert "Indeterminate" in result.severity_category


# =====================================================================
# N4 -- NIHSS interpretation recommended IV tPA, inside hemorrhage reports
# =====================================================================

class TestNIHSSCarriesNoTreatmentAdvice:
    """`NIHSSCalculator` baked "IV tPA indicated" into every band above minor.

    `AcuteStrokeTriageWorkflow` runs NIHSS for every stroke type including hemorrhagic, and
    `NeuroReportExporter.export_stroke_report` prints each scale's interpretation verbatim, so
    an exported ICH report endorsed alteplase three lines above "Intracerebral haemorrhage
    identified". Alteplase in ICH is lethal. NIHSS measures severity and has no stroke-type
    input, so all thrombolysis language belongs to the workflow.
    """

    FULL_18 = {
        "1a_loc": 1, "1b_loc_questions": 2, "1c_loc_commands": 2, "2_gaze": 1,
        "3_visual": 1, "4_facial": 2, "5a_left_arm": 2, "5b_right_arm": 1,
        "6a_left_leg": 2, "6b_right_leg": 1, "7_ataxia": 1, "8_sensory": 1,
        "9_language": 1, "10_dysarthria": 0, "11_extinction": 0,
    }

    def test_the_moderate_to_severe_band_no_longer_says_tpa_indicated(self):
        result = NIHSSCalculator.calculate(self.FULL_18)
        assert result.score == 18
        assert "tpa" not in result.interpretation.lower()
        assert "IV tPA indicated" not in result.interpretation

    @pytest.mark.parametrize("target", [3, 10, 18, 23, 30])
    def test_no_band_mentions_thrombolysis_or_thrombectomy(self, target):
        scores = {}
        remaining = target
        for item, max_val in NIHSSCalculator.ITEMS.items():
            take = min(max_val, remaining)
            scores[item] = take
            remaining -= take
        result = NIHSSCalculator.calculate(scores)
        text = (result.interpretation + " " + " ".join(result.recommendations)).lower()
        for banned in ("tpa", "alteplase", "thrombolysis", "thrombectomy"):
            assert banned not in text, f"NIHSS {result.score} still mentions {banned}"

    def test_the_severity_description_is_still_there(self):
        result = NIHSSCalculator.calculate(self.FULL_18)
        assert result.severity_category == "Moderate-to-severe stroke"
        assert "deficit" in result.interpretation.lower()

    def test_the_exported_hemorrhagic_report_endorses_no_thrombolytic(self):
        from src.clinical_workflows import AcuteStrokeTriageWorkflow
        from src.export import NeuroReportExporter

        result = AcuteStrokeTriageWorkflow().run({
            "stroke_type_hint": "hemorrhagic",
            "onset_hours": 1.0,
            "nihss_scores": self.FULL_18,
        })
        report = NeuroReportExporter().export_stroke_report(result).lower()
        assert "intracerebral haemorrhage identified" in report
        for banned in ("tpa", "alteplase", "thrombolysis"):
            assert banned not in report, f"hemorrhagic report still mentions {banned}"


# =====================================================================
# N20 -- NIHSS from a partial exam scored the unexamined items 0
# =====================================================================

class TestPartialNIHSSIsNotAScore:
    """`scores.get(item, 0)` scored unexamined items as "no deficit".

    `NIHSSCalculator.calculate({"1a_loc": 2, "9_language": 3})` returned 5/42, "Moderate
    stroke", with nothing saying 13 of 15 items had never been examined. NIHSS >= 4 gates
    alteplase and >= 6 gates thrombectomy in this agent, so a partial exam systematically
    under-triages.
    """

    PARTIAL = {"1a_loc": 2, "9_language": 3}

    def test_the_partial_exam_is_flagged_incomplete(self):
        result = NIHSSCalculator.calculate(self.PARTIAL)
        assert result.complete is False
        assert result.items_scored == 2
        assert result.items_expected == 15

    def test_it_is_not_reported_as_moderate_stroke(self):
        result = NIHSSCalculator.calculate(self.PARTIAL)
        assert result.severity_category != "Moderate stroke"
        assert "Indeterminate" in result.severity_category

    def test_the_interpretation_says_it_is_a_lower_bound(self):
        result = NIHSSCalculator.calculate(self.PARTIAL)
        assert "lower bound" in result.interpretation
        assert "2 of 15" in result.interpretation

    def test_a_complete_exam_is_still_scored_normally(self):
        full = {item: 0 for item in NIHSSCalculator.ITEMS}
        full["1a_loc"] = 2
        full["9_language"] = 3
        result = NIHSSCalculator.calculate(full)
        assert result.complete is True
        assert result.score == 5
        assert result.severity_category == "Moderate stroke"

    def test_a_mistyped_item_is_clamped_to_its_maximum(self):
        """The route twin applied no per-item maximum, so `4_facial: 30` passed straight in."""
        full = {item: 0 for item in NIHSSCalculator.ITEMS}
        full["4_facial"] = 30
        result = calculate_scale("nihss", full)
        assert result.score == 3
        assert any("outside the valid range" in w for w in result.warnings)

    def test_a_ui_item_name_is_resolved_rather_than_dropped(self):
        full = {item: 0 for item in NIHSSCalculator.ITEMS}
        del full["4_facial"]
        full["4_facial_palsy"] = 3
        result = calculate_scale("nihss", full)
        assert result.complete is True
        assert result.score == 3


# =====================================================================
# N5 -- ALSFRS-R scored every unrecorded item as 4 = normal function
# =====================================================================

class TestPartialALSFRSIsNotNormalFunction:
    """Two answered items returned 41/48, "Mild functional impairment. Early disease stage."

    `scores.get(item, 4)` meant the respiratory sub-score computed 12/12 out of defaults, so
    the PFT/BiPAP recommendation never fired, and the bulbar sub-score likewise. The decline
    rate, computed as (48 - total) / months, was understated by the same mechanism.
    """

    TWO_ITEMS = {"walking": 1, "climbing_stairs": 0}

    def test_two_items_no_longer_total_41(self):
        result = calculate_scale("alsfrs", self.TWO_ITEMS)
        assert result.score == 1
        assert result.score != 41

    def test_it_is_not_reported_as_mild_early_disease(self):
        result = calculate_scale("alsfrs", self.TWO_ITEMS)
        assert "Mild functional impairment" not in result.severity_category
        assert "Indeterminate" in result.severity_category
        assert "Early disease stage" not in result.interpretation

    def test_the_respiratory_subscore_is_not_computed_from_defaults(self):
        result = calculate_scale("alsfrs", self.TWO_ITEMS)
        assert any(
            "Respiratory sub-score not evaluated" in w for w in result.warnings
        )

    def test_the_bulbar_subscore_is_not_computed_from_defaults(self):
        result = calculate_scale("alsfrs", self.TWO_ITEMS)
        assert any("Bulbar sub-score not evaluated" in w for w in result.warnings)

    def test_the_decline_rate_is_not_computed_from_an_incomplete_instrument(self):
        result = ALSFRSCalculator.calculate(self.TWO_ITEMS, months_since_onset=12.0)
        assert "points/month" not in result.interpretation
        assert any("Decline rate not computed" in w for w in result.warnings)

    def test_the_ui_dressing_key_is_resolved_not_defaulted_to_4(self):
        """The UI posts `dressing`; the calculator expected `dressing_hygiene`."""
        items = {item: 0 for item in ALSFRSCalculator.ITEMS}
        del items["dressing_hygiene"]
        items["dressing"] = 2
        result = calculate_scale("alsfrs", items)
        assert result.complete is True
        assert result.score == 2

    def test_a_complete_instrument_still_bands_and_recommends(self):
        items = {item: 1 for item in ALSFRSCalculator.ITEMS}
        result = ALSFRSCalculator.calculate(items, months_since_onset=12.0)
        assert result.complete is True
        assert result.score == 12
        joined = " ".join(result.recommendations).lower()
        assert "pulmonary function" in joined
        assert "speech-language pathology" in joined


# =====================================================================
# N9 -- HIT-6 scored unanswered and unrecognized responses as "never"
# =====================================================================

class TestHIT6RefusesToImputeTheBestAnswer:
    """Short lists were padded with "never" (6 points, the best answer on every item).

    `HIT6Calculator.calculate(["always", "always"])` returned 50, "Some impact" -- below the
    >= 56 preventive-therapy trigger this agent uses, so a patient who answered both available
    items at maximum severity was denied preventive therapy. `ITEM_VALUES.get(r.lower(), 6)`
    mapped unrecognized strings to "never" too, so the typo "very often" (space instead of
    underscore) scored 6 instead of 11.
    """

    def test_two_maximum_severity_answers_are_refused_not_scored_50(self):
        with pytest.raises(ScaleDataError) as exc:
            HIT6Calculator.calculate(["always", "always"])
        assert "exactly 6" in str(exc.value)

    def test_a_separator_typo_scores_as_very_often_not_as_never(self):
        """"very often" with a space used to fall through `.get(..., 6)` and score as
        "never" -- 6 instead of 11, on every one of six items."""
        result = HIT6Calculator.calculate(["very often"] * 6)
        assert result.score == 66
        assert result.score != 36

    def test_a_response_that_is_not_a_category_at_all_is_refused(self):
        with pytest.raises(ScaleDataError) as exc:
            HIT6Calculator.calculate(["frequently"] * 6)
        assert "not a recognised response" in str(exc.value)

    def test_a_non_hit6_weight_is_refused(self):
        """The route twin summed six 1-5 Likert values to ~6-30 = "No significant impact"."""
        with pytest.raises(ScaleDataError) as exc:
            calculate_scale("hit6", {f"q{i}": 3 for i in range(1, 7)})
        assert "not a HIT-6 item weight" in str(exc.value)

    def test_six_maximum_answers_score_78_severe(self):
        result = HIT6Calculator.calculate(["always"] * 6)
        assert result.score == 78
        assert result.severity_category == "Severe impact"

    def test_the_ui_weight_payload_is_accepted(self):
        result = calculate_scale("hit6", {
            "pain_severity": 13, "daily_activities": 13, "lie_down": 11,
            "too_tired": 11, "fed_up": 10, "concentration": 10,
        })
        assert result.score == 68
        assert "preventive" in " ".join(result.recommendations).lower()


# =====================================================================
# N17 -- ASPECTS silently dropped unrecognized region names
# =====================================================================

class TestASPECTSDoesNotImproveOnUnknownRegions:
    """`[r.upper() for r in affected if r.upper() in REGIONS]` discarded the rest in silence.

    `calculate(["Caudate", "Insula", "M1 cortex"])` returned 10/10, "Favorable for
    intervention" -- a radiologist writing words instead of codes turned a 7/10 CT into a
    perfect one, and this score feeds the thrombectomy gate in the stroke workflow.
    """

    def test_long_region_names_are_scored_not_dropped(self):
        result = ASPECTSCalculator.calculate(["Caudate", "Insula", "M1 cortex"])
        assert result.score == 7
        assert result.score != 10
        assert result.severity_category != "Favorable for intervention"

    def test_a_genuinely_unknown_token_stops_the_calculation(self):
        with pytest.raises(ScaleDataError) as exc:
            ASPECTSCalculator.calculate(["C", "left frontal lobe"])
        assert "not recognised" in str(exc.value)
        assert "left frontal lobe" in str(exc.value)

    def test_short_codes_still_work(self):
        result = ASPECTSCalculator.calculate(["C", "I", "M1"])
        assert result.score == 7

    def test_a_region_flagged_true_is_counted_as_affected(self):
        """The route twin counted only values equal to the string "0"."""
        result = calculate_scale("aspects", {"M1": True, "M2": False})
        assert result.score == 9

    def test_the_ui_one_or_zero_payload_is_read_correctly(self):
        items = {r: 1 for r in ASPECTSCalculator.REGIONS}
        items["M1"] = 0
        items["I"] = 0
        result = calculate_scale("aspects", items)
        assert result.score == 8

    def test_a_twelve_key_payload_cannot_drive_the_score_negative(self):
        items = {r: 0 for r in ASPECTSCalculator.REGIONS}
        with pytest.raises(ScaleDataError):
            calculate_scale("aspects", {**items, "extra_region": 0, "another": 0})
        assert calculate_scale("aspects", items).score == 0


# =====================================================================
# N8 -- EDSS 9.5 reported as death; FS scores ignored
# =====================================================================

class TestEDSS:
    """Two defects.

    `elif edss <= 9.0 ... else "Death due to MS"` labelled a living EDSS 9.5 patient dead.
    This repo's own `src/knowledge.py:CLINICAL_SCALES["edss"]` maps "7.0-9.5" to "Wheelchair /
    bed-bound" and "10.0" to "Death due to MS", so the band contradicted the agent's own
    documented scale.

    And `calculate(fs_scores, edss_step)` never read `fs_scores` while its docstring claimed it
    "validates" the score: `calculate({"pyramidal": 6, "cerebellar": 5, "brainstem": 5}, 0.0)`
    returned 0.0 "No or minimal disability. Fully ambulatory."
    """

    def test_edss_9_5_is_not_death(self):
        result = EDSSCalculator.calculate({}, 9.5)
        assert "death" not in result.severity_category.lower()
        assert "death" not in result.interpretation.lower()
        assert result.severity_category == "Restricted to bed"

    def test_edss_10_is_still_death(self):
        result = EDSSCalculator.calculate({}, 10.0)
        assert result.severity_category == "Death due to MS"

    def test_the_repo_reference_table_and_the_calculator_now_agree(self):
        from src.knowledge import CLINICAL_SCALES

        bands = CLINICAL_SCALES["edss"]["interpretation"]
        assert "death" in bands["10.0"].lower()
        assert "death" not in bands["7.0-9.5"].lower()
        assert "death" not in EDSSCalculator.calculate({}, 9.5).interpretation.lower()

    def test_fs_scores_inconsistent_with_the_step_are_flagged(self):
        result = EDSSCalculator.calculate(
            {"pyramidal": 6, "cerebellar": 5, "brainstem": 5}, 0.0
        )
        assert result.warnings
        assert any("inconsistent" in w for w in result.warnings)

    def test_fs_scores_are_reported_rather_than_discarded(self):
        result = EDSSCalculator.calculate({"pyramidal": 6}, 6.0)
        assert any("pyramidal=6" in r for r in result.recommendations)

    def test_a_missing_step_is_indeterminate_not_derived_from_max_fs(self):
        """The route twin set `edss = min(float(max_fs), 9.5)`, so a blind but fully
        ambulatory patient with visual FS 6 was assigned EDSS 6.0 -- "requires a walking
        aid"."""
        result = calculate_scale("edss", {"visual": 6, "ambulation": 0.0})
        assert result.complete is False
        assert "INSUFFICIENT DATA" in result.interpretation
        assert result.severity_category == "Indeterminate -- EDSS step not recorded"

    def test_an_explicit_step_is_still_interpreted(self):
        result = calculate_scale("edss", {"visual": 6, "ambulation": 2.0})
        assert result.score == 2.0
        assert result.severity_category == "Mild disability"


# =====================================================================
# Shared-implementation guarantee (the structural finding)
# =====================================================================

class TestOneImplementationPerScale:
    """The route-local `_calculate_*` set was a second implementation, and the UI called it.

    It disagreed with `src/clinical_scales.py` on bandings, on defaults and on eligibility,
    which is why a tPA fix that landed in `clinical_workflows.py` did not change what a
    clinician saw. Both now resolve to `calculate_scale`.
    """

    CASES = [
        ("nihss", {item: 1 for item in NIHSSCalculator.ITEMS}),
        ("moca", dict(MoCACalculator.DOMAINS)),
        ("alsfrs", {item: 2 for item in ALSFRSCalculator.ITEMS}),
        ("aspects", {r: 1 for r in ASPECTSCalculator.REGIONS}),
        ("mrs", {"score": 3}),
        ("gcs", {"eye": 3, "verbal": 4, "motor": 5}),
        ("hoehn_yahr", {"stage": 3.0}),
    ]

    @pytest.mark.parametrize("scale,items", CASES)
    def test_route_and_library_return_the_same_answer(self, scale, items):
        import api.routes.neuro_clinical as routes

        route_result = routes._SCALE_CALCULATORS[scale](items)
        library_result = calculate_scale(scale, items)
        assert route_result["total_score"] == library_result.score
        assert route_result["interpretation"] == library_result.interpretation
        assert route_result["severity_category"] == library_result.severity_category

    def test_every_advertised_scale_has_an_adapter(self):
        import api.routes.neuro_clinical as routes

        advertised = {
            "nihss", "gcs", "moca", "updrs", "edss", "mrs",
            "hit6", "alsfrs", "aspects", "hoehn_yahr",
        }
        assert advertised <= set(routes._SCALE_CALCULATORS)


# =====================================================================
# N18 -- three mutually contradictory bandings for the same scale
# =====================================================================

class TestN18OneBandingPerScale:
    """`src/knowledge.py:CLINICAL_SCALES`, `src/clinical_scales.py` and the route-local
    `_calculate_*` set each carried their own severity bands.

    UPDRS-III: knowledge.py said 0-10 / 11-22 / 23-40 / 41-132, both calculators said
    0-10 / 11-32 / 33-58 / 59-80 / 81+, and the route said 0-10 / 11-32 / 33-58 / 59+. A score
    of 35 was "Severe" per the knowledge base -- which is what the RAG and `/scales` surfaces
    return -- and "Moderate" per both calculators.

    ALSFRS-R: knowledge.py said 39-48 / 25-38 / 13-24 / 0-12, both calculators said
    >= 40 / >= 30 / >= 20 / < 20. A score of 39 was "Mild" in one and "Moderate" in the other;
    22 was "Severe" against "Significant".

    NIHSS: `/stroke/triage`'s pre-computed-total fallback omitted the 16-20
    "moderate-to-severe" band present in both the calculator and knowledge.py, so the same
    NIHSS 18 was "moderate_severe" when submitted as items and "severe" when submitted as a
    total.

    The bands are reconciled to the calculators -- the table the clinical logic acts on -- and
    the uncited ones are labelled as local convention rather than resolved from memory.
    """

    def _band_for(self, scale, score):
        from src.knowledge import CLINICAL_SCALES

        for span, label in CLINICAL_SCALES[scale]["interpretation"].items():
            parts = span.split("-")
            low = float(parts[0])
            high = float(parts[-1])
            if low <= score <= high:
                return label
        return None

    def test_updrs_35_gets_one_answer(self):
        items = {item: 0 for item in UPDRSCalculator.ITEMS}
        remaining = 35
        for item in UPDRSCalculator.ITEMS:
            take = min(4, remaining)
            items[item] = take
            remaining -= take
        calculated = UPDRSCalculator.calculate(items)
        assert calculated.score == 35
        assert calculated.severity_category == "Moderate"
        assert "Moderate" in self._band_for("updrs", 35)

    def test_alsfrs_39_gets_one_answer(self):
        """39 was "Mild" in knowledge.py and "Moderate" in both calculators."""
        items = {item: 3 for item in ALSFRSCalculator.ITEMS}
        items["speech"] = 4
        items["salivation"] = 4
        items["swallowing"] = 4  # 36 + 3 = 39
        calculated = ALSFRSCalculator.calculate(items)
        assert calculated.score == 39
        assert calculated.severity_category == "Moderate functional impairment"
        assert self._band_for("alsfrs", 39) == "Moderate functional impairment"

    def test_alsfrs_40_is_mild_on_both_surfaces(self):
        items = {item: 3 for item in ALSFRSCalculator.ITEMS}
        for extra in ("speech", "salivation", "swallowing", "handwriting"):
            items[extra] = 4  # 36 + 4 = 40
        calculated = ALSFRSCalculator.calculate(items)
        assert calculated.score == 40
        assert calculated.severity_category == "Mild functional impairment"
        assert self._band_for("alsfrs", 40) == "Mild functional impairment"

    def test_alsfrs_22_gets_one_answer(self):
        assert "Significant" in self._band_for("alsfrs", 22)

    def test_nihss_18_gets_one_answer_from_both_input_shapes(self):
        import asyncio
        import types

        from api.routes.neuro_clinical import StrokeTriageRequest, stroke_triage

        class _Req:
            app = types.SimpleNamespace(
                state=types.SimpleNamespace(
                    llm_client=None, metrics=None, metrics_lock=None
                )
            )

        items = {item: 0 for item in NIHSSCalculator.ITEMS}
        remaining = 18
        for item, max_val in NIHSSCalculator.ITEMS.items():
            take = min(max_val, remaining)
            items[item] = take
            remaining -= take

        as_items = asyncio.run(stroke_triage(
            StrokeTriageRequest(nihss_items=items), _Req()
        ))
        as_total = asyncio.run(stroke_triage(
            StrokeTriageRequest(nihss_total=18), _Req()
        ))
        assert as_total.nihss_severity == "moderate_severe"
        assert "Moderate-to-severe" in as_items.nihss_severity
        assert "Moderate to severe" in self._band_for("nihss", 18)
