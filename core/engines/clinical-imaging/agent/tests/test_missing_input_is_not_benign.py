"""An absent or unrecognised input must not be scored as the most benign value.

One bug class, found in four of this engine's workflows by the 2026-11-03 pre-release audit.
In every instance the code reached for a default when an input was missing or unrecognised,
and in every instance the default available was the reassuring end of the scale: zero points
in a TI-RADS feature table, 0.0 mm for an unmeasured nodule, "solid" for an unstated nodule
density, category NORMAL and severity ROUTINE for a classifier that produced nothing. The
result each time is a report that looks complete and reads as benign about a patient who was
never assessed -- which is worse than an error, because an error gets investigated.

The cases here are the audit's own reproductions. Each one of them passed silently on the
release candidate with all 1,406 imaging tests green.

Fixes recorded alongside:
  A-5  Lung-RADS graded part-solid and ground-glass nodules with the solid table
  A-6  an unrecognised TI-RADS feature label scored 0 -- the benign value in all five tables
  A-7  an absent nodule diameter was scored as the smallest band
  A-14 the CXR workflow returned COMPLETED / ROUTINE / NORMAL on no classifier output
  A-15 a positive cardiomegaly was filed under FindingCategory.NORMAL
  A-16 the TI-RADS cross-modal genomic trigger was attached to the wrong nodule
  A-17 highest_tirads was tracked by severity rank, so a TR2-only gland reported TR1
  A-18 canned demo docstrings described categories their own code does not compute
"""
import inspect

import pytest

from src.models import (
    FindingCategory,
    FindingSeverity,
    LungRADS,
    TIRADS,
    WorkflowStatus,
)
from src.workflows.ct_chest_lung_nodule import (
    CTChestLungNoduleWorkflow,
    LungRADSGradingError,
    classify_nodule_lung_rads,
    classify_solid_nodule_lung_rads,
    resolve_nodule_type,
)
from src.workflows.cxr_rapid_findings import (
    CXR_CLASS_TO_CATEGORY,
    CXRRapidFindingsWorkflow,
)
from src.workflows.thyroid_tirads import (
    TIRADSFeatureError,
    ThyroidTIRADSWorkflow,
    calculate_tirads_points,
    normalize_feature_label,
    points_to_tirads,
    score_tirads_features,
    tirads_fna_recommendation,
)


# =====================================================================
# A-6 — an unrecognised TI-RADS feature label used to score zero
# =====================================================================

#: The audit's TR5 nodule, written the way a feature extractor or a parsed report plausibly
#: emits it: capitalised, hyphenated, spaced. ACR TI-RADS gives this 13 points -- solid 2 +
#: very hypoechoic 3 + taller-than-wide 3 + irregular 2 + punctate echogenic foci 3 -- which is
#: TR5, FNA indicated at >= 1 cm. The old code returned 0 points and
#: "TR1: Benign. No FNA indicated. No follow-up needed."
NON_CANONICAL_TR5 = {
    "id": "nodule_1",
    "max_diameter_mm": 12.0,
    "composition": "Solid",
    "echogenicity": "very hypoechoic",
    "shape": "taller-than-wide",
    "margin": "Irregular",
    "echogenic_foci": "punctate echogenic foci",
}

CANONICAL_TR5 = {
    "id": "nodule_1",
    "max_diameter_mm": 12.0,
    "composition": "solid",
    "echogenicity": "very_hypoechoic",
    "shape": "taller_than_wide",
    "margin": "irregular",
    "echogenic_foci": "punctate_echogenic_foci",
}


class TestNonCanonicalSpellingScoresTheSameThirteenPoints:
    def test_the_canonical_nodule_is_thirteen_points_tr5(self):
        assert calculate_tirads_points(CANONICAL_TR5) == 13
        assert points_to_tirads(13) is TIRADS.TR5

    def test_spelling_variants_score_identically_rather_than_zero(self):
        assert calculate_tirads_points(NON_CANONICAL_TR5) == 13

    @pytest.mark.parametrize("raw,expected", [
        ("Solid", "solid"),
        ("very hypoechoic", "very_hypoechoic"),
        ("taller-than-wide", "taller_than_wide"),
        ("  Mixed Cystic-Solid  ", "mixed_cystic_solid"),
        ("PUNCTATE  ECHOGENIC   FOCI", "punctate_echogenic_foci"),
    ])
    def test_normalisation_covers_case_spaces_and_hyphens(self, raw, expected):
        assert normalize_feature_label(raw) == expected

    def test_the_full_gland_is_reported_as_tr5_with_fna(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [NON_CANONICAL_TR5]})
        assert result.status is WorkflowStatus.COMPLETED
        assert result.classification == "TI-RADS TR5"
        assert result.severity is FindingSeverity.CRITICAL
        assert "FNA recommended" in result.findings[0]["recommendation"]


class TestAnUnrecognisedFeatureLabelIsNeverWorthZero:
    @pytest.mark.parametrize("feature", [
        "composition", "echogenicity", "shape", "margin", "echogenic_foci",
    ])
    def test_a_label_no_table_contains_is_reported_not_scored(self, feature):
        nodule = dict(CANONICAL_TR5, **{feature: "a vendor word for this"})
        score = score_tirads_features(nodule)
        assert not score.is_complete
        assert feature in score.unrecognised
        assert feature not in score.scored

    @pytest.mark.parametrize("feature", [
        "composition", "echogenicity", "shape", "margin", "echogenic_foci",
    ])
    def test_an_absent_feature_is_reported_not_scored(self, feature):
        nodule = {k: v for k, v in CANONICAL_TR5.items() if k != feature}
        score = score_tirads_features(nodule)
        assert not score.is_complete
        assert feature in score.missing

    def test_calculate_tirads_points_refuses_and_names_the_feature(self):
        with pytest.raises(TIRADSFeatureError) as exc:
            calculate_tirads_points(dict(CANONICAL_TR5, margin="frilly"))
        assert "margin" in str(exc.value)
        assert "frilly" in str(exc.value)

    def test_the_workflow_reports_indeterminate_rather_than_tr1(self):
        """The whole point: no TR level, no benign recommendation, and a failed status."""
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [dict(CANONICAL_TR5, margin="frilly")]})

        assert result.status is WorkflowStatus.FAILED
        assert result.severity is FindingSeverity.UNKNOWN
        assert "TR1" not in result.classification
        assert "margin" in result.classification

        finding = result.findings[0]
        assert finding["tirads"] is None
        assert finding["category"] == FindingCategory.INDETERMINATE.value
        assert finding["severity"] == FindingSeverity.UNKNOWN.value
        assert "No FNA indicated" not in finding["recommendation"]
        assert finding["unscored_features"]["unrecognised"] == {"margin": "frilly"}

    def test_an_unscored_feature_does_not_publish_a_point_total(self):
        """A `nodule_1_tirads_points` of 0.0 is indistinguishable from a benign nodule."""
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [dict(CANONICAL_TR5, margin="frilly")]})
        assert "nodule_1_tirads_points" not in result.measurements


# =====================================================================
# A-7 — an absent nodule size used to be scored as the smallest size
# =====================================================================


class TestAnUnmeasuredNoduleGetsNoSizeDependentRecommendation:
    def test_tirads_refuses_to_recommend_without_a_diameter(self):
        """This used to read "TR5 ... No FNA needed at current size (< 0.5cm)"."""
        rec = tirads_fna_recommendation(TIRADS.TR5, None)
        assert "CANNOT be determined" in rec
        assert "not measured" in rec
        assert "0.5cm" not in rec
        assert "not a recommendation against FNA" in rec

    @pytest.mark.parametrize("tr", [TIRADS.TR3, TIRADS.TR4, TIRADS.TR5])
    def test_every_size_dependent_level_refuses(self, tr):
        assert "CANNOT be determined" in tirads_fna_recommendation(tr, None)

    @pytest.mark.parametrize("tr", [TIRADS.TR1, TIRADS.TR2])
    def test_the_two_levels_that_do_not_depend_on_size_still_answer(self, tr):
        """At TR1 and TR2 the published management is the same at every size."""
        assert "No FNA indicated" in tirads_fna_recommendation(tr, None)

    def test_the_workflow_marks_an_unmeasured_tr5_gland_indeterminate(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        nodule = {k: v for k, v in CANONICAL_TR5.items() if k != "max_diameter_mm"}
        result = wf.postprocess({"nodules": [nodule]})

        assert result.status is WorkflowStatus.FAILED
        assert result.severity is FindingSeverity.UNKNOWN
        assert "max_diameter_mm" in result.classification
        # The TR level was computable from the features and is still reported honestly.
        assert result.findings[0]["tirads"] == TIRADS.TR5.value
        assert "CANNOT be determined" in result.findings[0]["recommendation"]
        assert "nodule_1_diameter_mm" not in result.measurements

    def test_lung_rads_refuses_to_grade_an_unmeasured_solid_nodule(self):
        """This used to be "category 2 ... Solid nodule < 6mm -- benign appearance"."""
        with pytest.raises(LungRADSGradingError) as exc:
            classify_nodule_lung_rads("solid", max_diameter_mm=None)
        assert "no measured diameter" in str(exc.value)


# =====================================================================
# A-5 — Lung-RADS graded every density with the solid-nodule table
# =====================================================================


class TestLungRadsWillNotGradeADensityItHasNoTableFor:
    """The replacement thresholds are deliberately absent, so the refusal is the whole fix.

    The ACR's Lung-RADS v2022 asset could not be retrieved (it redirects to a tenant page and
    serves no document), so no part-solid or non-solid band is implemented here. Inventing one
    would repeat the defect in a new direction. See CLINICAL_REVIEW.md Lung-RADS entry.
    """

    def test_the_solid_table_is_unchanged(self):
        assert classify_solid_nodule_lung_rads(4.0) is LungRADS.CAT_2
        assert classify_solid_nodule_lung_rads(7.0) is LungRADS.CAT_3
        assert classify_solid_nodule_lung_rads(12.0) is LungRADS.CAT_4A
        assert classify_solid_nodule_lung_rads(18.0) is LungRADS.CAT_4B

    def test_a_solid_nodule_still_grades(self):
        category, from_mm = classify_nodule_lung_rads("solid", max_diameter_mm=7.2)
        assert category is LungRADS.CAT_3
        assert from_mm == 7.2

    @pytest.mark.parametrize("density", [
        "part-solid", "Part Solid", "part_solid", "subsolid",
    ])
    def test_part_solid_is_refused_not_graded_on_its_solid_component(self, density):
        """The audit's case: 25 mm part-solid, 9.3 mm solid component, reported 4A."""
        with pytest.raises(LungRADSGradingError) as exc:
            classify_nodule_lung_rads(
                density, max_diameter_mm=25.0, solid_component_diameter_mm=9.3
            )
        assert "part-solid" in str(exc.value)
        assert "Christensen" in str(exc.value)

    @pytest.mark.parametrize("density", [
        "ground glass", "ground-glass", "non-solid", "GGN", "pure ground glass",
    ])
    def test_ground_glass_is_refused_not_called_4b_critical(self, density):
        """A 20 mm pure ground-glass nodule used to come back 4B, CRITICAL, "tissue sampling"."""
        with pytest.raises(LungRADSGradingError) as exc:
            classify_nodule_lung_rads(density, max_diameter_mm=20.0)
        assert "non-solid" in str(exc.value)

    def test_the_density_labels_this_module_understands(self):
        assert resolve_nodule_type("part-solid") == "part_solid"
        assert resolve_nodule_type("ground glass") == "non_solid"
        assert resolve_nodule_type("Solid") == "solid"

    def test_an_unstated_density_is_not_assumed_solid(self):
        for absent in (None, "", "   "):
            with pytest.raises(LungRADSGradingError) as exc:
                classify_nodule_lung_rads(absent, max_diameter_mm=9.0)
            assert "not reported" in str(exc.value)

    def test_an_unrecognised_density_is_refused(self):
        with pytest.raises(LungRADSGradingError) as exc:
            classify_nodule_lung_rads("spongy", max_diameter_mm=9.0)
        assert "spongy" in str(exc.value)

    def test_the_workflow_grades_the_solid_nodule_and_refuses_the_part_solid_one(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        result = wf.postprocess(wf._mock_inference())

        assert result.status is WorkflowStatus.FAILED
        assert result.severity is FindingSeverity.UNKNOWN
        assert "incomplete" in result.classification

        by_id = {f["nodule_id"]: f for f in result.findings}
        assert by_id["nodule_1"]["lung_rads"] == LungRADS.CAT_3.value
        assert by_id["nodule_2"]["lung_rads"] is None
        assert by_id["nodule_2"]["category"] == FindingCategory.INDETERMINATE.value
        assert by_id["nodule_2"]["severity"] == FindingSeverity.UNKNOWN.value
        assert "annual" not in by_id["nodule_2"]["recommendation"].lower()


# =====================================================================
# A-14 / A-15 — the CXR workflow
# =====================================================================


class TestAChestFilmThatWasNotReadIsNotANormalFilm:
    def test_no_classifier_output_is_failed_unknown_indeterminate(self):
        """Was: status COMPLETED, severity ROUTINE, one finding with category NORMAL."""
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.postprocess({"class_probabilities": {}})

        assert result.status is WorkflowStatus.FAILED
        assert result.severity is FindingSeverity.UNKNOWN
        assert result.classification == "not read: no classifier output"

        finding = result.findings[0]
        assert finding["category"] == FindingCategory.INDETERMINATE.value
        assert finding["severity"] == FindingSeverity.UNKNOWN.value
        assert finding["missing_input"] == "class_probabilities"
        assert "NOT READ" in finding["description"]

    def test_no_finding_on_that_path_claims_to_be_normal(self):
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.postprocess({"class_probabilities": {}})
        categories = [f["category"] for f in result.findings]
        assert FindingCategory.NORMAL.value not in categories

    def test_a_label_with_no_threshold_is_not_thresholded_at_half(self):
        """An unrecognised class used to take threshold 0.50 and category NORMAL."""
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.postprocess({"class_probabilities": {"aortic_dissection": 0.91}})

        assert result.status is WorkflowStatus.FAILED
        assert result.severity is FindingSeverity.UNKNOWN
        assert "NOT EVALUATED" in result.classification
        finding = result.findings[0]
        assert finding["category"] == FindingCategory.INDETERMINATE.value
        assert finding["above_threshold"] is None

    def test_a_genuinely_negative_film_still_reports_normal(self):
        """The fix must not turn every quiet film into a failure."""
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.postprocess({"class_probabilities": {
            "pneumothorax": 0.01, "consolidation": 0.02, "pleural_effusion": 0.03,
            "cardiomegaly": 0.04, "fracture": 0.05,
        }})
        assert result.status is WorkflowStatus.COMPLETED
        assert result.classification == "negative"
        assert result.findings[0]["category"] == FindingCategory.NORMAL.value


class TestAPositiveCardiomegalyIsNotCategorisedNormal:
    def test_the_class_maps_to_enlargement(self):
        assert CXR_CLASS_TO_CATEGORY["cardiomegaly"] is FindingCategory.ENLARGEMENT

    def test_a_positive_cardiomegaly_survives_a_category_filter(self):
        """Filtering on `category != normal` -- a worklist, an echo cohort -- used to drop it."""
        wf = CXRRapidFindingsWorkflow(mock_mode=True)
        result = wf.postprocess({"class_probabilities": {"cardiomegaly": 0.95}})

        assert result.classification == "positive: cardiomegaly"
        finding = result.findings[0]
        assert finding["category"] == FindingCategory.ENLARGEMENT.value
        assert finding["category"] != FindingCategory.NORMAL.value
        assert [f for f in result.findings
                if f["category"] != FindingCategory.NORMAL.value]


# =====================================================================
# A-16 / A-17 — which nodule, and which level
# =====================================================================

BENIGN_CYST = {
    "id": "benign_n1", "location": "right lobe", "laterality": "right",
    "max_diameter_mm": 9.0, "composition": "cystic", "echogenicity": "anechoic",
    "shape": "wider_than_tall", "margin": "smooth", "echogenic_foci": "none",
}

TR5_CANCER = dict(CANONICAL_TR5, id="cancer_n2", location="left lobe", laterality="left")


class TestTheGenomicTriggerIsAttachedToTheSuspiciousNodule:
    def test_the_trigger_follows_the_tr5_nodule_not_the_first_one(self):
        """Was: the BRAF/RAS/RET-PTC panel recorded against a 0-point cyst in the other lobe."""
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [BENIGN_CYST, TR5_CANCER]})
        by_id = {f["nodule_id"]: f for f in result.findings}

        assert by_id["cancer_n2"]["cross_modal_trigger"] is True
        assert "BRAF V600E" in by_id["cancer_n2"]["genomic_queries"]
        assert "cross_modal_trigger" not in by_id["benign_n1"]
        assert "genomic_queries" not in by_id["benign_n1"]

    def test_the_annotated_nodule_is_in_the_right_lobe(self):
        """An FNA ordered off this record has to point at the correct side of the neck."""
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [BENIGN_CYST, TR5_CANCER]})
        annotated = [f for f in result.findings if f.get("cross_modal_trigger")]
        assert len(annotated) == 1
        assert annotated[0]["laterality"] == "left"

    def test_every_tr4_or_tr5_nodule_is_annotated(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        tr4 = dict(BENIGN_CYST, id="tr4", composition="solid", echogenicity="hypoechoic")
        result = wf.postprocess({"nodules": [BENIGN_CYST, tr4, TR5_CANCER]})
        annotated = {f["nodule_id"] for f in result.findings if f.get("cross_modal_trigger")}
        assert annotated == {"tr4", "cancer_n2"}

    def test_a_benign_only_gland_triggers_nothing(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [BENIGN_CYST]})
        assert not [f for f in result.findings if f.get("cross_modal_trigger")]


class TestTheHeadlineTrLevelMatchesTheFindings:
    def test_a_tr2_only_gland_is_reported_as_tr2(self):
        """Was TR1: the maximum was tracked by severity rank, and TR1 and TR2 share one."""
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        tr2 = dict(BENIGN_CYST, composition="mixed_cystic_solid", echogenicity="isoechoic")
        result = wf.postprocess({"nodules": [tr2]})

        assert result.findings[0]["tirads"] == TIRADS.TR2.value
        assert result.classification == "TI-RADS TR2"

    @pytest.mark.parametrize("composition,echogenicity,expected", [
        ("cystic", "anechoic", TIRADS.TR1),
        ("mixed_cystic_solid", "isoechoic", TIRADS.TR2),
        ("solid", "isoechoic", TIRADS.TR3),
        ("solid", "hypoechoic", TIRADS.TR4),
    ])
    def test_the_headline_never_contradicts_its_own_body(
        self, composition, echogenicity, expected
    ):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        nodule = dict(BENIGN_CYST, composition=composition, echogenicity=echogenicity)
        result = wf.postprocess({"nodules": [nodule]})
        assert result.findings[0]["tirads"] == expected.value
        assert result.classification == f"TI-RADS {expected.value}"

    def test_the_highest_level_in_a_multi_nodule_gland_wins(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [BENIGN_CYST, TR5_CANCER]})
        assert result.classification == "TI-RADS TR5"


# =====================================================================
# A-18 — the demo docstrings described categories the code does not compute
# =====================================================================


class TestTheCannedDemoNarrativesMatchWhatTheCodeComputes:
    """A drifting docstring is what reviewers read instead of running the code.

    In the lung workflow the stale narrative ("Lung-RADS 4A") is what masked A-5 for a
    release, so these figures are asserted rather than trusted.
    """

    def test_the_thyroid_demo_nodules_are_tr5_twelve_points_and_tr2_two_points(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        nodules = wf._mock_inference()["nodules"]

        first = score_tirads_features(nodules[0])
        second = score_tirads_features(nodules[1])
        assert (first.points, points_to_tirads(first.points)) == (12, TIRADS.TR5)
        assert (second.points, points_to_tirads(second.points)) == (2, TIRADS.TR2)

    @staticmethod
    def _narrative(workflow_class):
        """The part of _mock_inference's docstring that DESCRIBES the canned nodules.

        Read off the class dict, because __init__ rebinds the instance attribute to the
        canned-payload marking wrapper. The audit-fix note after it is excluded: it quotes the
        wrong old figures on purpose, which is the point of the note.
        """
        doc = inspect.getdoc(workflow_class.__dict__["_mock_inference"]) or ""
        return doc.split("AUDIT FIX")[0]

    def test_the_thyroid_docstring_quotes_those_figures(self):
        doc = self._narrative(ThyroidTIRADSWorkflow)
        assert "TR5, 12 points" in doc
        assert "TR2, 2 points" in doc
        assert "TR3, 2 points" not in doc
        assert "TR5, 10 points" not in doc

    def test_the_lung_docstring_no_longer_calls_the_part_solid_nodule_4a(self):
        doc = self._narrative(CTChestLungNoduleWorkflow)
        assert "Lung-RADS 4A" not in doc
        assert "9.3mm solid" in doc
        assert "NOT GRADED" in doc

    def test_the_lung_demo_part_solid_nodule_really_has_a_93mm_solid_component(self):
        wf = CTChestLungNoduleWorkflow(mock_mode=True)
        second = wf._mock_inference()["nodules"][1]
        assert second["type"] == "part-solid"
        assert second["solid_component_diameter_mm"] == 9.3


# =====================================================================
# The exported report has to carry the same refusal
# =====================================================================


class TestTheExportedSrDoesNotCodeAnUngradedStudyAsRoutine:
    """The nine MAP services shared one `severity_map.get(severity, (..., "Routine"))`.

    A severity the map does not contain -- which is now every ungraded study, since the
    workflows report `unknown` -- was written into a DICOM SR document as a ROUTINE finding,
    leaving the service as a graded result. No SNOMED concept is invented for "not graded":
    the content item becomes TEXT and says so in words.

    The thyroid MAP is exercised here because it had a second, worse instance of the same
    class: every value in its report header was read from a measurement key the workflow does
    not publish, so each `.get` default applied, and EVERY thyroid study was exported as
    "TI-RADS points: 0, TI-RADS category: TR1, maximum nodule size 0.0 mm, FNA recommended:
    No" -- including the TR5 cancer in its own demo payload.
    """

    @staticmethod
    def _sr_text(result):
        import io

        import pydicom

        from maps.thyroid_tirads.app import ThyroidTIRADSApp

        raw = ThyroidTIRADSApp().build_dicom_sr(result, None)
        assert raw is not None
        dataset = pydicom.dcmread(io.BytesIO(raw), force=True)
        return dataset, dataset.ContentSequence[0].TextValue

    def test_the_demo_tr5_gland_exports_its_real_score(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        dataset, text = self._sr_text(wf.run())

        assert "TI-RADS points: 12" in text
        assert "TI-RADS category: TI-RADS TR5" in text
        assert "Maximum nodule size: 22.0 mm" in text
        assert "FNA recommended: Yes" in text
        assert "TI-RADS points: 0" not in text
        assert "category: TR1" not in text
        assert "Nodule 1: 0.0 mm" not in text

    def test_an_ungraded_gland_exports_no_score_and_no_routine_severity(self):
        wf = ThyroidTIRADSWorkflow(mock_mode=True)
        result = wf.postprocess({"nodules": [dict(CANONICAL_TR5, margin="frilly")]})
        dataset, text = self._sr_text(result)

        assert "TI-RADS points: not reported" in text
        assert "NOT ASSESSED" in text
        assert "SEVERITY: UNKNOWN" in text

        items = list(dataset.ContentSequence)[1:]
        values = [getattr(i, "TextValue", "") for i in items]
        assert any("severity not determined" in v for v in values)
        assert any("TI-RADS score not determined" in v for v in values)
        # Nothing was coded Routine, and no score of 0 was asserted as a measurement.
        for item in items:
            for coded in getattr(item, "ConceptCodeSequence", []):
                assert coded.CodeMeaning != "Routine"
            for measured in getattr(item, "MeasuredValueSequence", []):
                assert str(measured.NumericValue) != "0"

    @pytest.mark.parametrize("module_name", [
        "maps.breast_birads.app",
        "maps.ct_chest_lung_nodule.app",
        "maps.ct_coronary_angiography.app",
        "maps.ct_head_hemorrhage.app",
        "maps.cxr_rapid_findings.app",
        "maps.liver_lirads.app",
        "maps.mri_brain_ms_lesion.app",
        "maps.mri_prostate_pirads.app",
        "maps.thyroid_tirads.app",
    ])
    def test_no_map_service_still_falls_back_to_routine(self, module_name):
        """All nine carried the identical line, so all nine are checked for its absence."""
        import importlib

        import re

        source = inspect.getsource(importlib.import_module(module_name))
        # The comment left in place deliberately quotes the old line, so comments are stripped
        # before looking for it. The `"routine": ("17621005", "Routine")` dict entry is a
        # genuine mapping and is left alone -- what must not come back is that tuple used as
        # the DEFAULT for a severity the dict does not contain.
        code = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("#")
        )
        assert not re.search(r"severity_map\.get\(\s*severity\s*,", code)
        assert "coded_severity = severity_map.get(severity)" in code
        assert "severity not determined" in code
