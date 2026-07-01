"""
Tests for the Precision Autoimmune Agent diagnostic engine module.

Covers:
- ACR/EULAR classification criteria evaluation (RA, SLE, and others)
- Differential diagnosis from antibody panels and HLA alleles
- Diagnostic odyssey timeline analysis
- Overlap syndrome detection
- Edge cases: empty inputs, unsupported diseases, missing data

All tests run offline — no external services required.
"""

import pytest

from src.diagnostic_engine import (
    DiagnosticEngine,
)
from src.models import AutoimmuneDisease

# ── Fixtures ──────────────────────────────────────────────────────────


@pytest.fixture
def engine():
    return DiagnosticEngine(agent=None, rag_engine=None)


@pytest.fixture
def ra_clinical_data_meets():
    """Clinical data that meets RA classification criteria (>=6 points)."""
    return {
        # joint_involvement: 4-10_small = 3 points
        "4-10_small": True,
        # serology: high_positive_RF_or_CCP = 3 points
        "high_positive_RF_or_CCP": True,
        # acute_phase: abnormal_CRP_or_ESR = 1 point
        "abnormal_CRP_or_ESR": True,
        # duration: >=6_weeks = 1 point
        ">=6_weeks": True,
        # Total: 3 + 3 + 1 + 1 = 8 >= 6
    }


@pytest.fixture
def ra_clinical_data_below():
    """Clinical data that does NOT meet RA classification (< 6 points)."""
    return {
        "1_large": True,        # 0 points
        "negative_RF_negative_CCP": True,  # 0 points
        "normal_CRP_normal_ESR": True,     # 0 points
        "<6_weeks": True,       # 0 points
        # Total: 0
    }


@pytest.fixture
def sle_clinical_data_meets():
    """Clinical data that meets SLE classification criteria (>=10 points)."""
    return {
        # musculoskeletal: joint_involvement = 6 points
        "joint_involvement": True,
        # immunology: anti_dsDNA_or_anti_Smith = 6 points
        "anti_dsDNA_or_anti_Smith": True,
        # Total: 12 >= 10
    }


@pytest.fixture
def sle_clinical_data_below():
    """Clinical data that does NOT meet SLE classification."""
    return {
        "fever": True,          # 2 points
        "alopecia": True,       # 2 points
        # Total: 4 < 10
    }


@pytest.fixture
def timeline_events():
    """Timeline events for a diagnostic odyssey."""
    return [
        {
            "event_type": "symptom_onset",
            "event_date": "2021-03-15",
            "description": "Joint pain and morning stiffness",
            "specialty": "Primary Care",
            "provider": "Dr. Adams",
        },
        {
            "event_type": "misdiagnosis",
            "event_date": "2021-05-20",
            "description": "Initially treated as osteoarthritis",
            "specialty": "Orthopedics",
            "provider": "Dr. Brown",
        },
        {
            "event_type": "lab_result",
            "event_date": "2021-09-10",
            "description": "Positive RF and anti-CCP",
            "specialty": "Rheumatology",
            "provider": "Dr. Carter",
        },
        {
            "event_type": "imaging",
            "event_date": "2021-10-01",
            "description": "MRI shows synovitis in MCP joints",
            "specialty": "Radiology",
            "provider": "Dr. Davis",
        },
        {
            "event_type": "diagnosis",
            "event_date": "2021-10-15",
            "description": "Diagnosed with RA, meets ACR/EULAR criteria",
            "specialty": "Rheumatology",
            "provider": "Dr. Carter",
        },
    ]


# ── Classification criteria evaluation tests ──────────────────────────


class TestEvaluateClassificationCriteria:

    def test_ra_meets_criteria(self, engine, ra_clinical_data_meets):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            ra_clinical_data_meets,
        )
        assert result["meets_criteria"] is True
        assert result["total_points"] >= 6
        assert result["disease"] == "rheumatoid_arthritis"
        assert result["criteria_set"] == "2010 ACR/EULAR RA Classification"

    def test_ra_below_threshold(self, engine, ra_clinical_data_below):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            ra_clinical_data_below,
        )
        assert result["meets_criteria"] is False
        assert result["total_points"] < 6

    def test_ra_met_criteria_list(self, engine, ra_clinical_data_meets):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            ra_clinical_data_meets,
        )
        assert len(result["met_criteria"]) > 0
        # Each met criterion should show category, item, and points
        for item in result["met_criteria"]:
            assert "(+" in item

    def test_ra_unmet_criteria_list(self, engine, ra_clinical_data_below):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            ra_clinical_data_below,
        )
        assert len(result["unmet_criteria"]) > 0

    def test_sle_meets_criteria(self, engine, sle_clinical_data_meets):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.SYSTEMIC_LUPUS,
            sle_clinical_data_meets,
        )
        assert result["meets_criteria"] is True
        assert result["total_points"] >= 10
        assert result["disease"] == "systemic_lupus_erythematosus"

    def test_sle_below_threshold(self, engine, sle_clinical_data_below):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.SYSTEMIC_LUPUS,
            sle_clinical_data_below,
        )
        assert result["meets_criteria"] is False
        assert result["total_points"] < 10

    def test_sle_has_entry_criterion(self, engine, sle_clinical_data_meets):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.SYSTEMIC_LUPUS,
            sle_clinical_data_meets,
        )
        assert result["entry_criterion"] is not None
        assert "ANA" in result["entry_criterion"]

    def test_empty_clinical_data(self, engine):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            {},
        )
        assert result["meets_criteria"] is False
        assert result["total_points"] == 0
        assert len(result["met_criteria"]) == 0

    def test_unsupported_disease(self, engine):
        # GRAVES_DISEASE has no classification criteria defined
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.GRAVES_DISEASE,
            {"some_item": True},
        )
        assert result["supported"] is False
        assert "No criteria defined" in result["message"]

    def test_ankylosing_spondylitis(self, engine):
        data = {"sacroiliitis_on_imaging_plus_1_SpA_feature": True}
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.ANKYLOSING_SPONDYLITIS,
            data,
        )
        assert result["meets_criteria"] is True
        assert result["entry_criterion"] is not None

    def test_systemic_sclerosis(self, engine):
        data = {"skin_thickening_fingers_extending_proximal_to_MCP": True}
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.SYSTEMIC_SCLEROSIS,
            data,
        )
        assert result["meets_criteria"] is True
        assert result["total_points"] == 9

    def test_sjogrens_meets(self, engine):
        data = {
            "anti_SSA_Ro_positive": True,   # 3 points
            "schirmer_test_le_5mm_5min": True,  # 1 point
            # Total: 4 >= 4
        }
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.SJOGRENS_SYNDROME,
            data,
        )
        assert result["meets_criteria"] is True

    def test_multiple_sclerosis(self, engine):
        data = {
            "ge_1_T2_lesion_in_ge_2_MS_typical_CNS_regions": True,  # 1 point
            "CSF_oligoclonal_bands": True,  # 1 point
            # Total: 2 >= 2
        }
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.MULTIPLE_SCLEROSIS,
            data,
        )
        assert result["meets_criteria"] is True

    def test_threshold_value_returned(self, engine, ra_clinical_data_meets):
        result = engine.evaluate_classification_criteria(
            AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
            ra_clinical_data_meets,
        )
        assert result["threshold"] == 6


# ── Differential diagnosis tests ──────────────────────────────────────


class TestGenerateDifferential:

    def test_returns_ranked_list(self, engine):
        result = engine.generate_differential(["ANA", "anti-dsDNA"])
        assert isinstance(result, list)
        assert len(result) > 0

    def test_sle_ranked_high_with_specific_antibodies(self, engine):
        result = engine.generate_differential(["ANA", "anti-dsDNA", "anti-Smith"])
        # SLE should be top-ranked with these highly specific markers
        assert result[0]["disease"] == "systemic_lupus_erythematosus"

    def test_ra_ranked_high_with_rf_and_ccp(self, engine):
        result = engine.generate_differential(["RF", "anti-CCP"])
        diseases = [r["disease"] for r in result]
        assert "rheumatoid_arthritis" in diseases
        # RA should be ranked highly
        ra_entry = next(r for r in result if r["disease"] == "rheumatoid_arthritis")
        assert ra_entry["rank"] <= 2

    def test_evidence_includes_antibody_info(self, engine):
        result = engine.generate_differential(["anti-CCP"])
        ra_entry = next(r for r in result if r["disease"] == "rheumatoid_arthritis")
        assert any("anti-CCP" in e for e in ra_entry["evidence"])

    def test_score_is_positive(self, engine):
        result = engine.generate_differential(["ANA"])
        for entry in result:
            assert entry["score"] > 0

    def test_rank_order(self, engine):
        result = engine.generate_differential(["ANA", "RF", "anti-dsDNA"])
        for i, entry in enumerate(result):
            assert entry["rank"] == i + 1

    def test_with_hla_alleles(self, engine):
        result = engine.generate_differential(
            positive_antibodies=["ANA"],
            hla_alleles=["DRB1*04:01"],
        )
        diseases = [r["disease"] for r in result]
        assert "rheumatoid_arthritis" in diseases

    def test_hla_evidence_included(self, engine):
        result = engine.generate_differential(
            positive_antibodies=[],
            hla_alleles=["DRB1*04:01"],
        )
        if result:
            ra_entry = next((r for r in result if r["disease"] == "rheumatoid_arthritis"), None)
            if ra_entry:
                assert any("HLA" in e for e in ra_entry["evidence"])

    def test_empty_antibodies(self, engine):
        result = engine.generate_differential([])
        assert result == []

    def test_unknown_antibody_ignored(self, engine):
        result = engine.generate_differential(["nonexistent-antibody-xyz"])
        assert result == []

    def test_ssc_with_specific_antibodies(self, engine):
        result = engine.generate_differential(["anti-Scl-70", "ANA"])
        diseases = [r["disease"] for r in result]
        assert "systemic_sclerosis" in diseases


# ── Diagnostic odyssey analysis tests ─────────────────────────────────


class TestAnalyzeDiagnosticOdyssey:

    def test_returns_expected_keys(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert "first_symptom_date" in result
        assert "diagnosis_date" in result
        assert "diagnostic_delay" in result
        assert "num_specialists_seen" in result
        assert "num_misdiagnoses" in result
        assert "key_diagnostic_tests" in result
        assert "total_events" in result

    def test_first_symptom_detected(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert result["first_symptom_date"] == "2021-03-15"

    def test_diagnosis_date_detected(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert result["diagnosis_date"] == "2021-10-15"

    def test_diagnostic_delay_calculated(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        delay = result["diagnostic_delay"]
        assert delay is not None
        assert delay["days"] > 0
        # 2021-03-15 to 2021-10-15 = 214 days
        assert delay["days"] == 214
        assert delay["months"] > 0
        assert delay["years"] > 0

    def test_specialists_counted(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert result["num_specialists_seen"] == 4
        assert "Rheumatology" in result["specialists"]
        assert "Radiology" in result["specialists"]

    def test_misdiagnoses_detected(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert result["num_misdiagnoses"] == 1
        assert result["misdiagnoses"][0]["wrong_diagnosis"] == "Initially treated as osteoarthritis"

    def test_key_diagnostic_tests(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        tests = result["key_diagnostic_tests"]
        assert len(tests) == 2  # lab_result + imaging
        test_types = {t["type"] for t in tests}
        assert "lab_result" in test_types
        assert "imaging" in test_types

    def test_total_events(self, engine, timeline_events):
        result = engine.analyze_diagnostic_odyssey(timeline_events)
        assert result["total_events"] == 5

    def test_empty_events(self, engine):
        result = engine.analyze_diagnostic_odyssey([])
        assert result == {"status": "no_data"}

    def test_no_symptom_onset(self, engine):
        events = [
            {"event_type": "diagnosis", "event_date": "2023-01-01", "description": "RA"},
        ]
        result = engine.analyze_diagnostic_odyssey(events)
        assert result["first_symptom_date"] is None
        assert result["diagnostic_delay"] is None

    def test_no_diagnosis_event(self, engine):
        events = [
            {"event_type": "symptom_onset", "event_date": "2023-01-01", "description": "Pain"},
            {"event_type": "lab_result", "event_date": "2023-02-01", "description": "Positive ANA"},
        ]
        result = engine.analyze_diagnostic_odyssey(events)
        assert result["diagnosis_date"] is None
        assert result["diagnostic_delay"] is None

    def test_events_sorted_by_date(self, engine):
        """Events provided out of order should still be processed correctly."""
        events = [
            {"event_type": "diagnosis", "event_date": "2023-06-01", "description": "RA", "specialty": "Rheumatology"},
            {"event_type": "symptom_onset", "event_date": "2023-01-01", "description": "Pain", "specialty": "Primary Care"},
        ]
        result = engine.analyze_diagnostic_odyssey(events)
        assert result["first_symptom_date"] == "2023-01-01"
        assert result["diagnosis_date"] == "2023-06-01"

    def test_missing_specialty(self, engine):
        events = [
            {"event_type": "symptom_onset", "event_date": "2023-01-01", "description": "Pain"},
        ]
        result = engine.analyze_diagnostic_odyssey(events)
        assert result["num_specialists_seen"] == 0


# ── Overlap syndrome detection tests ──────────────────────────────────


class TestDetectOverlapSyndromes:

    def test_mctd_detected_with_rnp(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-RNP"])
        syndromes = [d["syndrome"] for d in detected]
        assert "mixed_connective_tissue_disease" in syndromes

    def test_sle_ra_overlap_with_shared_markers(self, engine):
        detected = engine.detect_overlap_syndromes(["ANA", "anti-CCP", "RF"])
        syndromes = [d["syndrome"] for d in detected]
        assert "overlap_syndrome_sle_ra" in syndromes

    def test_sle_sjogrens_overlap(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-SSA/Ro", "ANA"])
        syndromes = [d["syndrome"] for d in detected]
        assert "sjogrens_sle_overlap" in syndromes

    def test_ra_sjogrens_overlap(self, engine):
        detected = engine.detect_overlap_syndromes(["RF", "anti-SSA/Ro", "ANA"])
        syndromes = [d["syndrome"] for d in detected]
        assert "ra_sjogrens_overlap" in syndromes

    def test_confidence_high_with_multiple_markers(self, engine):
        detected = engine.detect_overlap_syndromes(["ANA", "anti-CCP", "RF"])
        sle_ra = next(d for d in detected if d["syndrome"] == "overlap_syndrome_sle_ra")
        assert sle_ra["confidence"] == "high"

    def test_confidence_moderate_with_single_marker(self, engine):
        detected = engine.detect_overlap_syndromes(["ANA"])
        # Syndromes with only 1 shared marker match should be "moderate"
        for d in detected:
            if len(d["matched_markers"]) < 2:
                assert d["confidence"] == "moderate"

    def test_diseases_involved_populated(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-RNP"])
        mctd = next(d for d in detected if d["syndrome"] == "mixed_connective_tissue_disease")
        assert len(mctd["diseases_involved"]) > 0
        assert "systemic_lupus_erythematosus" in mctd["diseases_involved"]

    def test_no_overlap_with_irrelevant_antibodies(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-mitochondrial"])
        assert detected == []

    def test_empty_antibodies(self, engine):
        detected = engine.detect_overlap_syndromes([])
        assert detected == []

    def test_ssc_myositis_overlap(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-Pm-Scl"])
        syndromes = [d["syndrome"] for d in detected]
        assert "ssc_myositis_overlap" in syndromes

    def test_multiple_overlaps_detected(self, engine):
        """A broad antibody panel should trigger multiple overlap syndromes."""
        detected = engine.detect_overlap_syndromes(
            ["ANA", "anti-CCP", "RF", "anti-SSA/Ro", "anti-RNP"]
        )
        assert len(detected) >= 3

    def test_t1d_celiac_overlap(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-tTG IgA", "HLA-DQ2"])
        syndromes = [d["syndrome"] for d in detected]
        assert "t1d_celiac_overlap" in syndromes

    def test_lupus_aps_overlap(self, engine):
        detected = engine.detect_overlap_syndromes(["anti-cardiolipin", "lupus anticoagulant"])
        syndromes = [d["syndrome"] for d in detected]
        assert "lupus_aps_overlap" in syndromes
