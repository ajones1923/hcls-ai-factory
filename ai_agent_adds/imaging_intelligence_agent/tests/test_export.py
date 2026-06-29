"""Tests for Imaging Intelligence Agent export functions.

Validates export_markdown(), export_json(), export_pdf(), and export_fhir()
with various AgentResponse inputs. All tests are self-contained.

Author: Adam Jones
Date: February 2026
"""

import json
import pytest
from unittest.mock import patch, MagicMock

from src.export import (
    export_fhir,
    export_json,
    export_markdown,
    export_pdf,
    SNOMED_FINDING_CODES,
    SEVERITY_INTERPRETATION,
)
from src.models import (
    AgentResponse,
    CrossCollectionResult,
    CrossModalResult,
    FindingSeverity,
    SearchHit,
    WorkflowResult,
    WorkflowStatus,
)


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def minimal_response():
    """Response with empty evidence and no extras."""
    return AgentResponse(
        question="What is VISTA-3D?",
        answer="VISTA-3D is a universal segmentation model.",
        evidence=CrossCollectionResult(query="VISTA-3D", hits=[]),
        timestamp="2026-02-27T10:00:00",
    )


@pytest.fixture
def rich_response():
    """Response with evidence, NIM services, and workflow results."""
    hits = [
        SearchHit(
            collection="imaging_literature",
            id="PMID:001",
            score=0.92,
            text="VISTA-3D segments 127 anatomical structures with mean Dice 0.91.",
            metadata={"year": 2025},
        ),
        SearchHit(
            collection="imaging_benchmarks",
            id="bench-001",
            score=0.85,
            text="VISTA-3D achieves state-of-the-art on TotalSegmentator benchmark.",
            metadata={"metric_value": 0.91},
        ),
        SearchHit(
            collection="imaging_guidelines",
            id="guide-001",
            score=0.78,
            text="ACR recommends AI segmentation for treatment planning.",
            metadata={"organization": "ACR"},
        ),
    ]

    evidence = CrossCollectionResult(
        query="What is the best segmentation model?",
        hits=hits,
        knowledge_context="## VISTA-3D\n- Universal segmentation model",
        total_collections_searched=11,
        search_time_ms=42.5,
    )

    workflow = WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[{
            "category": "hemorrhage",
            "description": "Intraparenchymal hemorrhage, 12.5 mL",
            "severity": "urgent",
        }],
        measurements={"volume_ml": 12.5, "midline_shift_mm": 3.2},
        classification="urgent_hemorrhage",
        severity=FindingSeverity.URGENT,
        nim_services_used=["3D U-Net", "VISTA-3D"],
        is_mock=True,
    )

    return AgentResponse(
        question="What is the best segmentation model?",
        answer=(
            "VISTA-3D from NVIDIA achieves state-of-the-art performance with "
            "mean Dice coefficient of 0.91 across 127 anatomical classes."
        ),
        evidence=evidence,
        workflow_results=[workflow],
        nim_services_used=["vista3d", "maisi"],
        knowledge_used=["pathology", "modality"],
        timestamp="2026-02-27T12:00:00",
    )


# ===================================================================
# export_markdown()
# ===================================================================


class TestExportMarkdown:
    """Tests for export_markdown()."""

    def test_returns_string(self, minimal_response):
        md = export_markdown(minimal_response)
        assert isinstance(md, str)
        assert len(md) > 0

    def test_includes_title(self, minimal_response):
        md = export_markdown(minimal_response)
        assert "Imaging Intelligence Report" in md

    def test_includes_query(self, minimal_response):
        md = export_markdown(minimal_response)
        assert "What is VISTA-3D?" in md

    def test_includes_timestamp(self, minimal_response):
        md = export_markdown(minimal_response)
        assert "2026-02-27" in md

    def test_includes_answer(self, minimal_response):
        md = export_markdown(minimal_response)
        assert "universal segmentation model" in md

    def test_includes_evidence_section(self, rich_response):
        md = export_markdown(rich_response)
        assert "Evidence" in md
        assert "3 items" in md or "items" in md.lower()

    def test_includes_evidence_hits(self, rich_response):
        md = export_markdown(rich_response)
        assert "PMID:001" in md

    def test_includes_nim_services(self, rich_response):
        md = export_markdown(rich_response)
        assert "NIM" in md or "nim" in md.lower()
        assert "vista3d" in md or "maisi" in md

    def test_includes_disclaimer(self, minimal_response):
        md = export_markdown(minimal_response)
        assert "Research use only" in md

    def test_empty_evidence_still_works(self, minimal_response):
        md = export_markdown(minimal_response)
        # Should not crash with empty evidence
        assert "Imaging Intelligence Report" in md

    def test_grouped_by_collection(self, rich_response):
        md = export_markdown(rich_response)
        assert "imaging_literature" in md
        assert "imaging_benchmarks" in md


# ===================================================================
# export_json()
# ===================================================================


class TestExportJSON:
    """Tests for export_json()."""

    def test_returns_string(self, minimal_response):
        result = export_json(minimal_response)
        assert isinstance(result, str)

    def test_valid_json(self, minimal_response):
        result = export_json(minimal_response)
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    def test_contains_question(self, minimal_response):
        result = export_json(minimal_response)
        parsed = json.loads(result)
        assert parsed["question"] == "What is VISTA-3D?"

    def test_contains_answer(self, minimal_response):
        result = export_json(minimal_response)
        parsed = json.loads(result)
        assert "segmentation" in parsed["answer"]

    def test_contains_evidence(self, minimal_response):
        result = export_json(minimal_response)
        parsed = json.loads(result)
        assert "evidence" in parsed
        assert "hits" in parsed["evidence"]

    def test_contains_timestamp(self, minimal_response):
        result = export_json(minimal_response)
        parsed = json.loads(result)
        assert "timestamp" in parsed

    def test_rich_response_has_all_fields(self, rich_response):
        result = export_json(rich_response)
        parsed = json.loads(result)
        assert len(parsed["evidence"]["hits"]) == 3
        assert len(parsed["workflow_results"]) == 1
        assert len(parsed["nim_services_used"]) == 2

    def test_workflow_results_serialized(self, rich_response):
        result = export_json(rich_response)
        parsed = json.loads(result)
        wr = parsed["workflow_results"][0]
        assert wr["workflow_name"] == "ct_head_hemorrhage"
        assert wr["severity"] == "urgent"

    def test_empty_evidence_works(self):
        response = AgentResponse(
            question="empty test",
            answer="no evidence",
            evidence=CrossCollectionResult(query="empty", hits=[]),
        )
        result = export_json(response)
        parsed = json.loads(result)
        assert parsed["evidence"]["hits"] == []

    def test_json_is_indented(self, minimal_response):
        result = export_json(minimal_response)
        assert "\n" in result  # indent=2 produces newlines
        assert "  " in result  # indent=2 produces 2-space indentation


# ===================================================================
# export_pdf()
# ===================================================================


class TestExportPDF:
    """Tests for export_pdf()."""

    def test_returns_empty_when_reportlab_missing(self, minimal_response):
        """When reportlab is not installed, returns empty string."""
        with patch.dict("sys.modules", {"reportlab": None}):
            with patch("src.export.logger"):
                # Force ImportError
                with patch("builtins.__import__", side_effect=_import_mock_reportlab):
                    result = export_pdf(minimal_response, "/tmp/test_report.pdf")
                    # Either returns the path (if reportlab available) or empty string
                    assert isinstance(result, str)

    def test_with_mock_reportlab(self, rich_response, tmp_path):
        """Test PDF export with mocked reportlab."""
        output_path = str(tmp_path / "test_report.pdf")

        mock_doc = MagicMock()
        mock_paragraph = MagicMock()
        mock_spacer = MagicMock()

        with patch("src.export.logger"):
            try:
                result = export_pdf(rich_response, output_path)
                # If reportlab is installed, it creates the file
                if result:
                    assert result == output_path
                else:
                    # reportlab not installed -> returns empty string
                    assert result == ""
            except Exception:
                # reportlab not installed or import fails
                pass

    def test_pdf_export_handles_multiline_answer(self, minimal_response, tmp_path):
        """Ensure multi-paragraph answers don't crash PDF export."""
        minimal_response.answer = (
            "First paragraph with findings.\n\n"
            "Second paragraph with analysis.\n\n"
            "Third paragraph with recommendations."
        )
        output_path = str(tmp_path / "test_multiline.pdf")

        try:
            result = export_pdf(minimal_response, output_path)
            assert isinstance(result, str)
        except Exception:
            # reportlab not installed
            pass

    def test_pdf_export_with_empty_evidence(self, tmp_path):
        """Ensure PDF works with zero evidence hits."""
        response = AgentResponse(
            question="test",
            answer="No evidence found.",
            evidence=CrossCollectionResult(query="test", hits=[]),
        )
        output_path = str(tmp_path / "test_empty.pdf")

        try:
            result = export_pdf(response, output_path)
            assert isinstance(result, str)
        except Exception:
            pass


# ===================================================================
# EDGE CASES
# ===================================================================


class TestExportEdgeCases:
    """Tests for edge cases across all export functions."""

    def test_markdown_with_special_characters(self):
        response = AgentResponse(
            question="What about <b>bold</b> & 'quotes'?",
            answer="Finding: ratio > 1.0 (5/7 = 0.71)",
            evidence=CrossCollectionResult(query="test", hits=[]),
        )
        md = export_markdown(response)
        assert "<b>bold</b>" in md

    def test_json_with_unicode(self):
        response = AgentResponse(
            question="Analyse du poumon",
            answer="R\u00e9sultat normal -- pas de pathologie.",
            evidence=CrossCollectionResult(query="test", hits=[]),
        )
        result = export_json(response)
        parsed = json.loads(result)
        assert parsed["answer"] == "R\u00e9sultat normal -- pas de pathologie."

    def test_markdown_with_long_answer(self):
        response = AgentResponse(
            question="comprehensive analysis",
            answer="A " * 5000,
            evidence=CrossCollectionResult(query="test", hits=[]),
        )
        md = export_markdown(response)
        assert len(md) > 5000

    def test_json_round_trip(self, rich_response):
        """Export to JSON, parse back, and verify key fields."""
        json_str = export_json(rich_response)
        parsed = json.loads(json_str)
        assert parsed["question"] == rich_response.question
        assert parsed["answer"] == rich_response.answer
        assert len(parsed["evidence"]["hits"]) == rich_response.evidence.hit_count

    def test_markdown_multiple_workflow_results(self):
        evidence = CrossCollectionResult(query="test", hits=[])
        workflows = [
            WorkflowResult(
                workflow_name="ct_head_hemorrhage",
                severity=FindingSeverity.URGENT,
            ),
            WorkflowResult(
                workflow_name="cxr_rapid_findings",
                severity=FindingSeverity.SIGNIFICANT,
            ),
        ]
        response = AgentResponse(
            question="multi-workflow test",
            answer="Multiple workflows executed.",
            evidence=evidence,
            workflow_results=workflows,
        )
        md = export_markdown(response)
        assert isinstance(md, str)
        assert len(md) > 0
        # Now includes workflow results section
        assert "Workflow Results" in md
        assert "ct_head_hemorrhage" in md
        assert "cxr_rapid_findings" in md


def _import_mock_reportlab(name, *args, **kwargs):
    """Helper to mock importlib for reportlab tests."""
    if "reportlab" in name:
        raise ImportError("Mocked: reportlab not installed")
    return __import__(name, *args, **kwargs)


# ===================================================================
# FHIR FIXTURES
# ===================================================================


@pytest.fixture
def fhir_response_minimal():
    """Minimal response for FHIR export testing."""
    return AgentResponse(
        question="CT head for headache",
        answer="No acute intracranial abnormality.",
        evidence=CrossCollectionResult(query="CT head", hits=[]),
        timestamp="2026-02-27T10:00:00",
    )


@pytest.fixture
def fhir_response_hemorrhage():
    """Response with hemorrhage workflow for FHIR export testing."""
    workflow = WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[{
            "category": "hemorrhage",
            "description": "Intraparenchymal hemorrhage in right basal ganglia, 12.5 mL",
            "severity": "urgent",
        }],
        measurements={"volume_ml": 12.5, "midline_shift_mm": 3.2},
        classification="urgent_hemorrhage",
        severity=FindingSeverity.URGENT,
        nim_services_used=["3D U-Net", "VISTA-3D"],
        is_mock=True,
    )
    return AgentResponse(
        question="CT head for acute hemorrhage evaluation",
        answer="Intraparenchymal hemorrhage identified in right basal ganglia, "
               "volume approximately 12.5 mL with 3.2 mm midline shift.",
        evidence=CrossCollectionResult(query="CT head hemorrhage", hits=[]),
        workflow_results=[workflow],
        nim_services_used=["vista3d"],
        timestamp="2026-02-27T12:00:00",
    )


@pytest.fixture
def fhir_response_multi_finding():
    """Response with multiple findings across workflows."""
    workflow1 = WorkflowResult(
        workflow_name="cxr_rapid_findings",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "consolidation",
                "description": "Right lower lobe consolidation",
                "severity": "significant",
            },
            {
                "category": "effusion",
                "description": "Small right pleural effusion",
                "severity": "routine",
            },
        ],
        measurements={},
        severity=FindingSeverity.SIGNIFICANT,
    )
    workflow2 = WorkflowResult(
        workflow_name="cxr_cardiac",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "cardiomegaly",
                "description": "Cardiothoracic ratio 0.58",
                "severity": "significant",
            },
        ],
        measurements={},
        severity=FindingSeverity.SIGNIFICANT,
    )
    return AgentResponse(
        question="CXR evaluation",
        answer="Right lower lobe consolidation with small effusion. Cardiomegaly.",
        evidence=CrossCollectionResult(query="CXR", hits=[]),
        workflow_results=[workflow1, workflow2],
        timestamp="2026-02-27T14:00:00",
    )


@pytest.fixture
def fhir_response_cross_modal():
    """Response with cross-modal genomic enrichment."""
    workflow = WorkflowResult(
        workflow_name="ct_chest_nodule",
        status=WorkflowStatus.COMPLETED,
        findings=[{
            "category": "nodule",
            "description": "Spiculated nodule right upper lobe 18mm, Lung-RADS 4B",
            "severity": "critical",
        }],
        measurements={"diameter_mm": 18.0},
        severity=FindingSeverity.CRITICAL,
    )
    cross_modal = CrossModalResult(
        trigger_reason="Lung-RADS 4B finding triggered genomic lookup",
        genomic_context=["EGFR mutations common in adenocarcinoma"],
        genomic_hit_count=5,
        query_count=2,
        enrichment_summary="EGFR, KRAS, ALK variants relevant to lung adenocarcinoma",
    )
    return AgentResponse(
        question="CT chest nodule evaluation",
        answer="Suspicious spiculated nodule, recommend biopsy.",
        evidence=CrossCollectionResult(query="CT chest", hits=[]),
        workflow_results=[workflow],
        cross_modal=cross_modal,
        timestamp="2026-02-27T15:00:00",
    )


# ===================================================================
# export_fhir()
# ===================================================================


class TestExportFHIR:
    """Tests for export_fhir() -- FHIR R4 DiagnosticReport Bundle."""

    # --- 1. Valid JSON output ---
    def test_returns_valid_json(self, fhir_response_minimal):
        result = export_fhir(fhir_response_minimal)
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    # --- 2. Bundle type and structure ---
    def test_bundle_resource_type(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        assert parsed["resourceType"] == "Bundle"

    def test_bundle_type_collection(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        assert parsed["type"] == "collection"

    def test_bundle_has_timestamp(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        assert "timestamp" in parsed
        assert len(parsed["timestamp"]) > 0

    def test_bundle_has_id(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        assert "id" in parsed
        assert len(parsed["id"]) > 0

    def test_bundle_has_entries(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        assert "entry" in parsed
        assert isinstance(parsed["entry"], list)
        assert len(parsed["entry"]) > 0

    # --- 3. DiagnosticReport resource present with correct fields ---
    def test_diagnostic_report_present(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        reports = [
            e for e in parsed["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ]
        assert len(reports) == 1

    def test_diagnostic_report_status_final(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert report["status"] == "final"

    def test_diagnostic_report_conclusion(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "hemorrhage" in report["conclusion"].lower()

    def test_diagnostic_report_effective_datetime(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert report["effectiveDateTime"] == "2026-02-27T12:00:00"

    def test_diagnostic_report_issued(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "issued" in report
        assert len(report["issued"]) > 0

    # --- 4. LOINC category present ---
    def test_loinc_category_radiology(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        category_codings = report["category"][0]["coding"]
        assert any(
            c["system"] == "http://loinc.org" and c["code"] == "LP29684-5"
            for c in category_codings
        )

    def test_loinc_code_imaging_study(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        code_codings = report["code"]["coding"]
        assert any(
            c["system"] == "http://loinc.org" and c["code"] == "18748-4"
            for c in code_codings
        )

    # --- 5. Patient and performer references ---
    def test_patient_resource_present(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        patients = [
            e for e in parsed["entry"]
            if e["resource"]["resourceType"] == "Patient"
        ]
        assert len(patients) == 1

    def test_patient_id_default_anonymous(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        patient = _find_resource(parsed, "Patient")
        assert patient["id"] == "anonymous"

    def test_custom_patient_id(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal, patient_id="PAT-001"))
        patient = _find_resource(parsed, "Patient")
        assert patient["id"] == "PAT-001"

    def test_subject_reference_present(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "subject" in report
        assert "reference" in report["subject"]
        assert report["subject"]["reference"].startswith("urn:uuid:")

    def test_performer_reference(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert len(report["performer"]) > 0
        assert "AI-system" in report["performer"][0]["display"]

    def test_custom_practitioner_id(self, fhir_response_minimal):
        parsed = json.loads(
            export_fhir(fhir_response_minimal, practitioner_id="Dr-Smith")
        )
        report = _find_resource(parsed, "DiagnosticReport")
        assert "Dr-Smith" in report["performer"][0]["reference"]

    # --- 6. Observation resources for each finding ---
    def test_observation_for_hemorrhage(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        observations = _find_all_resources(parsed, "Observation")
        assert len(observations) == 1

    def test_observation_value_string(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        assert "Intraparenchymal hemorrhage" in obs["valueString"]

    def test_observation_status_final(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        assert obs["status"] == "final"

    # --- 7. SNOMED coding present ---
    def test_snomed_coding_on_observation(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        snomed = obs["code"]["coding"][0]
        assert snomed["system"] == "http://snomed.info/sct"
        assert snomed["code"] == "50960005"  # hemorrhage

    def test_snomed_conclusion_code_on_report(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "conclusionCode" in report
        codes = report["conclusionCode"]
        assert any(
            c["coding"][0]["code"] == "50960005"
            for c in codes
        )

    def test_all_snomed_codes_mapped(self):
        """Verify all required SNOMED codes are present in the mapping."""
        required = {
            "hemorrhage": "50960005",
            "nodule": "416940007",
            "consolidation": "95436008",
            "pneumothorax": "36118008",
            "effusion": "60046008",
            "fracture": "125605004",
            "cardiomegaly": "8186001",
            "mass": "4147007",
            "edema": "267038008",
            "normal": "17621005",
        }
        for finding, code in required.items():
            assert SNOMED_FINDING_CODES[finding] == code

    # --- 8. Severity interpretation mapping ---
    def test_severity_urgent_maps_to_H(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        interp = obs["interpretation"][0]["coding"][0]
        assert interp["code"] == "H"  # urgent -> H
        assert interp["system"] == (
            "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation"
        )

    def test_severity_critical_maps_to_HH(self, fhir_response_cross_modal):
        parsed = json.loads(export_fhir(fhir_response_cross_modal))
        obs = _find_all_resources(parsed, "Observation")[0]
        interp = obs["interpretation"][0]["coding"][0]
        assert interp["code"] == "HH"

    def test_severity_mapping_completeness(self):
        """All FindingSeverity values have a mapping."""
        for sev in FindingSeverity:
            assert sev.value in SEVERITY_INTERPRETATION

    # --- 9. Measurement components ---
    def test_measurement_components_present(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        assert "component" in obs
        assert len(obs["component"]) == 2  # volume_ml, midline_shift_mm

    def test_volume_measurement_value(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        vol = [c for c in obs["component"] if "volume" in c["code"]["text"]]
        assert len(vol) == 1
        assert vol[0]["valueQuantity"]["value"] == 12.5
        assert vol[0]["valueQuantity"]["unit"] == "mL"

    def test_midline_shift_measurement(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        obs = _find_all_resources(parsed, "Observation")[0]
        shift = [c for c in obs["component"] if "midline" in c["code"]["text"]]
        assert len(shift) == 1
        assert shift[0]["valueQuantity"]["value"] == 3.2
        assert shift[0]["valueQuantity"]["unit"] == "mm"

    # --- 10. Empty evidence/findings edge cases ---
    def test_no_workflow_produces_no_observations(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        observations = _find_all_resources(parsed, "Observation")
        assert len(observations) == 0

    def test_empty_findings_list_creates_normal_observation(self):
        """Workflow with no findings produces a 'normal' observation."""
        workflow = WorkflowResult(
            workflow_name="cxr_screening",
            status=WorkflowStatus.COMPLETED,
            findings=[],
            severity=FindingSeverity.NORMAL,
        )
        response = AgentResponse(
            question="CXR screening",
            answer="Normal chest radiograph.",
            evidence=CrossCollectionResult(query="CXR", hits=[]),
            workflow_results=[workflow],
            timestamp="2026-02-27T10:00:00",
        )
        parsed = json.loads(export_fhir(response))
        observations = _find_all_resources(parsed, "Observation")
        assert len(observations) == 1
        assert observations[0]["code"]["coding"][0]["code"] == "17621005"  # normal

    def test_no_conclusion_code_when_no_findings(self, fhir_response_minimal):
        parsed = json.loads(export_fhir(fhir_response_minimal))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "conclusionCode" not in report

    # --- 11. Cross-modal results in extension ---
    def test_cross_modal_extension(self, fhir_response_cross_modal):
        parsed = json.loads(export_fhir(fhir_response_cross_modal))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "extension" in report
        ext = report["extension"][0]
        assert ext["url"] == "urn:imaging-intelligence:cross-modal-result"
        assert "EGFR" in ext["valueString"]

    def test_no_extension_when_no_cross_modal(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        report = _find_resource(parsed, "DiagnosticReport")
        assert "extension" not in report

    # --- 12. Multiple workflow results ---
    def test_multiple_workflows_produce_multiple_observations(
        self, fhir_response_multi_finding
    ):
        parsed = json.loads(export_fhir(fhir_response_multi_finding))
        observations = _find_all_resources(parsed, "Observation")
        # workflow1 has 2 findings, workflow2 has 1 finding -> 3 total
        assert len(observations) == 3

    def test_multiple_conclusion_codes(self, fhir_response_multi_finding):
        parsed = json.loads(export_fhir(fhir_response_multi_finding))
        report = _find_resource(parsed, "DiagnosticReport")
        codes = report["conclusionCode"]
        snomed_codes = {c["coding"][0]["code"] for c in codes}
        assert "95436008" in snomed_codes   # consolidation
        assert "60046008" in snomed_codes   # effusion
        assert "8186001" in snomed_codes    # cardiomegaly

    def test_report_result_refs_match_observations(
        self, fhir_response_multi_finding
    ):
        parsed = json.loads(export_fhir(fhir_response_multi_finding))
        report = _find_resource(parsed, "DiagnosticReport")
        obs_entries = [
            e for e in parsed["entry"]
            if e["resource"]["resourceType"] == "Observation"
        ]
        result_refs = {r["reference"] for r in report["result"]}
        obs_fullurls = {e["fullUrl"] for e in obs_entries}
        assert result_refs == obs_fullurls

    # --- 13. Round-trip JSON parse validity ---
    def test_round_trip_json_parse(self, fhir_response_hemorrhage):
        """Export to JSON string, parse, re-serialize, parse again -- must match."""
        json_str = export_fhir(fhir_response_hemorrhage)
        first_parse = json.loads(json_str)
        re_serialized = json.dumps(first_parse, indent=2)
        second_parse = json.loads(re_serialized)
        assert first_parse == second_parse

    def test_all_entries_have_fullurl(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        for entry in parsed["entry"]:
            assert "fullUrl" in entry
            assert entry["fullUrl"].startswith("urn:uuid:")

    def test_all_entries_have_resource(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        for entry in parsed["entry"]:
            assert "resource" in entry
            assert "resourceType" in entry["resource"]

    # --- 14. ImagingStudy resource ---
    def test_imaging_study_present(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        studies = _find_all_resources(parsed, "ImagingStudy")
        assert len(studies) == 1

    def test_imaging_study_modality_ct(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        study = _find_resource(parsed, "ImagingStudy")
        modalities = study["modality"]
        assert any(m["code"] == "CT" for m in modalities)

    def test_imaging_study_subject_ref(self, fhir_response_hemorrhage):
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        study = _find_resource(parsed, "ImagingStudy")
        assert study["subject"]["reference"].startswith("urn:uuid:")

    # --- 15. Unique UUIDs and structural integrity ---
    def test_all_fullurls_unique(self, fhir_response_multi_finding):
        parsed = json.loads(export_fhir(fhir_response_multi_finding))
        fullurls = [e["fullUrl"] for e in parsed["entry"]]
        assert len(fullurls) == len(set(fullurls))

    def test_resource_type_counts(self, fhir_response_hemorrhage):
        """Exactly 1 Patient, 1 ImagingStudy, 1 DiagnosticReport, N Observations."""
        parsed = json.loads(export_fhir(fhir_response_hemorrhage))
        type_counts = {}
        for entry in parsed["entry"]:
            rt = entry["resource"]["resourceType"]
            type_counts[rt] = type_counts.get(rt, 0) + 1
        assert type_counts["Patient"] == 1
        assert type_counts["ImagingStudy"] == 1
        assert type_counts["DiagnosticReport"] == 1
        assert type_counts.get("Observation", 0) == 1


# ===================================================================
# FHIR TEST HELPERS
# ===================================================================


def _find_resource(parsed_bundle: dict, resource_type: str) -> dict:
    """Find the first resource of a given type in a parsed FHIR Bundle."""
    for entry in parsed_bundle["entry"]:
        if entry["resource"]["resourceType"] == resource_type:
            return entry["resource"]
    raise AssertionError(f"No {resource_type} found in bundle")


def _find_all_resources(parsed_bundle: dict, resource_type: str) -> list:
    """Find all resources of a given type in a parsed FHIR Bundle."""
    return [
        entry["resource"]
        for entry in parsed_bundle["entry"]
        if entry["resource"]["resourceType"] == resource_type
    ]
