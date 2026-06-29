"""Tests for the Cardiology Intelligence Agent report export system.

Validates CardioReportExporter creation, Markdown/JSON/FHIR R4 export
formats, SEVERITY_COLORS mapping, REPORT_TEMPLATES dictionary, and
format helper functions.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import json
import unittest
from datetime import datetime

from src.export import (
    REPORT_TEMPLATES,
    SEVERITY_COLORS,
    VERSION,
    CardioReportExporter,
    _DISCLAIMER,
    _generate_filename,
    _now_display,
    _now_iso,
    _severity_badge,
    _severity_emoji,
    format_guideline_citation,
    format_risk_score_table,
)
from src.models import (
    SeverityLevel,
)


# =====================================================================
# CardioReportExporter creation
# =====================================================================


class TestCardioReportExporterCreation(unittest.TestCase):
    """Tests for CardioReportExporter instantiation."""

    def test_create_exporter(self):
        exporter = CardioReportExporter()
        self.assertIsInstance(exporter, CardioReportExporter)

    def test_exporter_has_export_markdown(self):
        exporter = CardioReportExporter()
        self.assertTrue(hasattr(exporter, "export_markdown"))
        self.assertTrue(callable(exporter.export_markdown))

    def test_exporter_has_export_json(self):
        exporter = CardioReportExporter()
        self.assertTrue(hasattr(exporter, "export_json"))
        self.assertTrue(callable(exporter.export_json))

    def test_exporter_has_export_fhir_r4(self):
        exporter = CardioReportExporter()
        self.assertTrue(hasattr(exporter, "export_fhir_r4"))
        self.assertTrue(callable(exporter.export_fhir_r4))

    def test_exporter_has_export_pdf(self):
        exporter = CardioReportExporter()
        self.assertTrue(hasattr(exporter, "export_pdf"))
        self.assertTrue(callable(exporter.export_pdf))


# =====================================================================
# VERSION constant
# =====================================================================


class TestVersion(unittest.TestCase):
    """Tests for the VERSION constant."""

    def test_version_is_string(self):
        self.assertIsInstance(VERSION, str)

    def test_version_is_semver(self):
        parts = VERSION.split(".")
        self.assertEqual(len(parts), 3)
        for part in parts:
            self.assertTrue(part.isdigit())

    def test_version_value(self):
        self.assertEqual(VERSION, "1.0.0")


# =====================================================================
# SEVERITY_COLORS
# =====================================================================


class TestSeverityColors(unittest.TestCase):
    """Tests for the SEVERITY_COLORS mapping."""

    def test_severity_colors_is_dict(self):
        self.assertIsInstance(SEVERITY_COLORS, dict)

    def test_severity_colors_not_empty(self):
        self.assertGreater(len(SEVERITY_COLORS), 0)

    def test_severity_colors_has_low(self):
        self.assertIn(SeverityLevel.LOW, SEVERITY_COLORS)

    def test_severity_colors_has_high(self):
        self.assertIn(SeverityLevel.HIGH, SEVERITY_COLORS)

    def test_severity_colors_has_critical(self):
        self.assertIn(SeverityLevel.CRITICAL, SEVERITY_COLORS)

    def test_severity_colors_has_moderate(self):
        self.assertIn(SeverityLevel.MODERATE, SEVERITY_COLORS)

    def test_all_colors_are_hex(self):
        for level, color in SEVERITY_COLORS.items():
            self.assertTrue(
                color.startswith("#"),
                f"Color for {level} does not start with #: {color}",
            )

    def test_all_colors_are_7_chars(self):
        for level, color in SEVERITY_COLORS.items():
            self.assertEqual(
                len(color), 7,
                f"Color for {level} is not 7 chars: {color}",
            )

    def test_low_color_is_green(self):
        self.assertEqual(SEVERITY_COLORS[SeverityLevel.LOW], "#2ECC71")

    def test_high_color_is_red(self):
        self.assertEqual(SEVERITY_COLORS[SeverityLevel.HIGH], "#E74C3C")

    def test_critical_color_is_purple(self):
        self.assertEqual(SEVERITY_COLORS[SeverityLevel.CRITICAL], "#8E44AD")

    def test_at_least_five_entries(self):
        self.assertGreaterEqual(len(SEVERITY_COLORS), 5)


# =====================================================================
# REPORT_TEMPLATES
# =====================================================================


class TestReportTemplates(unittest.TestCase):
    """Tests for REPORT_TEMPLATES dictionary."""

    def test_report_templates_is_dict(self):
        self.assertIsInstance(REPORT_TEMPLATES, dict)

    def test_report_templates_not_empty(self):
        self.assertGreater(len(REPORT_TEMPLATES), 0)

    def test_has_consultation_template(self):
        self.assertIn("consultation", REPORT_TEMPLATES)

    def test_has_risk_assessment_template(self):
        self.assertIn("risk_assessment", REPORT_TEMPLATES)

    def test_has_gdmt_optimization_template(self):
        self.assertIn("gdmt_optimization", REPORT_TEMPLATES)

    def test_has_workflow_template(self):
        self.assertIn("workflow", REPORT_TEMPLATES)

    def test_consultation_contains_patient_id(self):
        self.assertIn("{patient_id}", REPORT_TEMPLATES["consultation"])

    def test_consultation_contains_timestamp(self):
        self.assertIn("{timestamp}", REPORT_TEMPLATES["consultation"])

    def test_consultation_contains_version(self):
        self.assertIn("{version}", REPORT_TEMPLATES["consultation"])

    def test_risk_assessment_contains_patient_id(self):
        self.assertIn("{patient_id}", REPORT_TEMPLATES["risk_assessment"])

    def test_workflow_contains_workflow_type(self):
        self.assertIn("{workflow_type}", REPORT_TEMPLATES["workflow"])

    def test_all_templates_are_strings(self):
        for key, template in REPORT_TEMPLATES.items():
            self.assertIsInstance(
                template, str,
                f"Template '{key}' is not a string",
            )

    def test_all_templates_contain_header(self):
        for key, template in REPORT_TEMPLATES.items():
            self.assertIn("#", template, f"Template '{key}' has no markdown header")


# =====================================================================
# Helper functions
# =====================================================================


class TestNowIso(unittest.TestCase):
    """Tests for _now_iso helper."""

    def test_returns_string(self):
        result = _now_iso()
        self.assertIsInstance(result, str)

    def test_ends_with_z(self):
        result = _now_iso()
        self.assertTrue(result.endswith("Z"))

    def test_contains_t_separator(self):
        result = _now_iso()
        self.assertIn("T", result)

    def test_parseable_as_datetime(self):
        result = _now_iso()
        parsed = datetime.strptime(result, "%Y-%m-%dT%H:%M:%SZ")
        self.assertIsInstance(parsed, datetime)


class TestNowDisplay(unittest.TestCase):
    """Tests for _now_display helper."""

    def test_returns_string(self):
        result = _now_display()
        self.assertIsInstance(result, str)

    def test_contains_utc(self):
        result = _now_display()
        self.assertIn("UTC", result)


class TestSeverityBadge(unittest.TestCase):
    """Tests for _severity_badge helper."""

    def test_returns_string(self):
        result = _severity_badge(SeverityLevel.LOW)
        self.assertIsInstance(result, str)

    def test_contains_shield_url(self):
        result = _severity_badge(SeverityLevel.LOW)
        self.assertIn("img.shields.io", result)

    def test_contains_severity_label(self):
        result = _severity_badge(SeverityLevel.HIGH)
        self.assertIn("High", result)

    def test_low_badge_format(self):
        result = _severity_badge(SeverityLevel.LOW)
        self.assertIn("Severity", result)


class TestSeverityEmoji(unittest.TestCase):
    """Tests for _severity_emoji helper."""

    def test_returns_string(self):
        result = _severity_emoji(SeverityLevel.LOW)
        self.assertIsInstance(result, str)

    def test_low_returns_low_marker(self):
        result = _severity_emoji(SeverityLevel.LOW)
        self.assertEqual(result, "[LOW]")

    def test_high_returns_high_marker(self):
        result = _severity_emoji(SeverityLevel.HIGH)
        self.assertEqual(result, "[HIGH]")

    def test_critical_returns_critical_marker(self):
        result = _severity_emoji(SeverityLevel.CRITICAL)
        self.assertEqual(result, "[CRITICAL]")

    def test_moderate_returns_moderate_marker(self):
        result = _severity_emoji(SeverityLevel.MODERATE)
        self.assertEqual(result, "[MODERATE]")


class TestFormatGuidelineCitation(unittest.TestCase):
    """Tests for format_guideline_citation helper."""

    def test_returns_string(self):
        result = format_guideline_citation("ACC/AHA 2022")
        self.assertIsInstance(result, str)

    def test_contains_guideline_text(self):
        result = format_guideline_citation("ACC/AHA 2022 Heart Failure")
        self.assertIn("ACC/AHA 2022 Heart Failure", result)

    def test_is_blockquote(self):
        result = format_guideline_citation("Some guideline")
        self.assertTrue(result.startswith(">"))

    def test_empty_guideline_returns_empty(self):
        result = format_guideline_citation("")
        self.assertEqual(result, "")

    def test_none_like_empty_returns_empty(self):
        result = format_guideline_citation("")
        self.assertEqual(result, "")


class TestFormatRiskScoreTable(unittest.TestCase):
    """Tests for format_risk_score_table helper."""

    def test_no_scores_returns_placeholder(self):
        result = format_risk_score_table([])
        self.assertIn("No risk scores", result)

    def test_returns_string(self):
        result = format_risk_score_table([])
        self.assertIsInstance(result, str)


class TestGenerateFilename(unittest.TestCase):
    """Tests for _generate_filename helper."""

    def test_returns_string(self):
        result = _generate_filename("report", "md")
        self.assertIsInstance(result, str)

    def test_contains_prefix(self):
        result = _generate_filename("cardio_report", "pdf")
        self.assertTrue(result.startswith("cardio_report_"))

    def test_contains_extension(self):
        result = _generate_filename("report", "json")
        self.assertTrue(result.endswith(".json"))

    def test_contains_timestamp(self):
        result = _generate_filename("report", "md")
        self.assertIn("T", result)
        self.assertIn("Z", result)


class TestDisclaimer(unittest.TestCase):
    """Tests for the _DISCLAIMER constant."""

    def test_disclaimer_is_string(self):
        self.assertIsInstance(_DISCLAIMER, str)

    def test_disclaimer_contains_warning(self):
        self.assertIn("Disclaimer", _DISCLAIMER)

    def test_disclaimer_mentions_clinical(self):
        self.assertIn("clinical", _DISCLAIMER.lower())


# =====================================================================
# export_markdown
# =====================================================================


class TestExportMarkdown(unittest.TestCase):
    """Tests for CardioReportExporter.export_markdown."""

    def setUp(self):
        self.exporter = CardioReportExporter()

    def test_with_dict_response(self):
        response = {
            "severity": "low",
            "summary": "Normal cardiac assessment.",
            "findings": ["No abnormalities detected."],
            "recommendations": ["Continue current management."],
        }
        result = self.exporter.export_markdown(response, patient_id="P-100")
        self.assertIsInstance(result, str)

    def test_markdown_contains_patient_id(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response, patient_id="P-42")
        self.assertIn("P-42", result)

    def test_markdown_contains_date(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response)
        self.assertIn("UTC", result)

    def test_markdown_contains_headers(self):
        response = {"severity": "low", "summary": "Test summary"}
        result = self.exporter.export_markdown(response)
        self.assertIn("#", result)

    def test_markdown_contains_summary_section(self):
        response = {"severity": "low", "summary": "Test summary text"}
        result = self.exporter.export_markdown(response)
        self.assertIn("Clinical Summary", result)

    def test_markdown_contains_findings(self):
        response = {
            "severity": "low",
            "summary": "Test",
            "findings": ["Finding A", "Finding B"],
        }
        result = self.exporter.export_markdown(response)
        self.assertIn("Finding A", result)
        self.assertIn("Finding B", result)

    def test_markdown_contains_recommendations(self):
        response = {
            "severity": "low",
            "summary": "Test",
            "recommendations": ["Rec 1"],
        }
        result = self.exporter.export_markdown(response)
        self.assertIn("Rec 1", result)

    def test_markdown_contains_disclaimer(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response)
        self.assertIn("Disclaimer", result)

    def test_markdown_default_patient_id(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response)
        self.assertIn("DEMO", result)

    def test_markdown_report_type_consultation(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response, report_type="consultation")
        self.assertIn("Consultation", result)

    def test_markdown_report_type_risk_assessment(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response, report_type="risk_assessment")
        self.assertIn("Risk Assessment", result)

    def test_markdown_with_string_response(self):
        result = self.exporter.export_markdown("raw string response")
        self.assertIsInstance(result, str)

    def test_markdown_with_empty_dict(self):
        result = self.exporter.export_markdown({})
        self.assertIsInstance(result, str)

    def test_markdown_contains_version(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_markdown(response)
        self.assertIn(VERSION, result)


# =====================================================================
# export_json
# =====================================================================


class TestExportJson(unittest.TestCase):
    """Tests for CardioReportExporter.export_json."""

    def setUp(self):
        self.exporter = CardioReportExporter()

    def test_returns_dict(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_json(response)
        self.assertIsInstance(result, dict)

    def test_has_report_type_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("report_type", result)

    def test_has_version_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("version", result)

    def test_has_generated_at_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("generated_at", result)

    def test_has_patient_id_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("patient_id", result)

    def test_has_agent_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("agent", result)

    def test_has_data_key(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertIn("data", result)

    def test_patient_id_is_set(self):
        result = self.exporter.export_json({"summary": "Test"}, patient_id="P-99")
        self.assertEqual(result["patient_id"], "P-99")

    def test_default_patient_id(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertEqual(result["patient_id"], "DEMO")

    def test_version_matches(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertEqual(result["version"], VERSION)

    def test_agent_name(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertEqual(result["agent"], "cardiology_intelligence_agent")

    def test_report_type_value(self):
        result = self.exporter.export_json({"summary": "Test"})
        self.assertEqual(result["report_type"], "cardio_intelligence_report")

    def test_data_contains_input(self):
        result = self.exporter.export_json({"summary": "Test data"})
        self.assertIn("summary", result["data"])

    def test_json_serializable(self):
        result = self.exporter.export_json({"summary": "Test"})
        serialized = json.dumps(result)
        self.assertIsInstance(serialized, str)

    def test_with_string_response(self):
        result = self.exporter.export_json("raw string")
        self.assertIn("data", result)
        self.assertIn("raw", result["data"])


# =====================================================================
# export_fhir_r4
# =====================================================================


class TestExportFhirR4(unittest.TestCase):
    """Tests for CardioReportExporter.export_fhir_r4."""

    def setUp(self):
        self.exporter = CardioReportExporter()

    def test_has_export_fhir_r4_method(self):
        self.assertTrue(hasattr(self.exporter, "export_fhir_r4"))

    def test_export_fhir_r4_is_callable(self):
        self.assertTrue(callable(self.exporter.export_fhir_r4))

    def test_returns_dict(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_fhir_r4(response)
        self.assertIsInstance(result, dict)

    def test_has_resource_type(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_fhir_r4(response)
        self.assertIn("resourceType", result)

    def test_resource_type_is_diagnostic_report_or_bundle(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_fhir_r4(response)
        self.assertIn(
            result["resourceType"],
            ("DiagnosticReport", "Bundle"),
        )

    def test_fhir_with_patient_id(self):
        response = {"severity": "low", "summary": "Test"}
        result = self.exporter.export_fhir_r4(response, patient_id="P-123")
        self.assertIsInstance(result, dict)


# =====================================================================
# Export with mock CardioResponse data
# =====================================================================


class TestExportWithMockData(unittest.TestCase):
    """Tests for exports with realistic mock CardioResponse data."""

    def setUp(self):
        self.exporter = CardioReportExporter()
        self.mock_response = {
            "severity": "moderate",
            "summary": "Patient presents with moderate heart failure symptoms.",
            "findings": [
                "LVEF 35%",
                "Elevated NT-proBNP at 1200 pg/mL",
                "Mild mitral regurgitation",
            ],
            "recommendations": [
                "Initiate ARNI therapy",
                "Optimize beta-blocker dose",
                "Add SGLT2 inhibitor",
            ],
            "risk_scores": [],
            "cross_modal_triggers": [],
            "guideline_references": [
                "2022 AHA/ACC/HFSA Heart Failure Guideline",
            ],
            "confidence": 0.85,
        }

    def test_markdown_with_full_mock(self):
        result = self.exporter.export_markdown(
            self.mock_response, patient_id="P-HF-001"
        )
        self.assertIn("P-HF-001", result)
        self.assertIn("LVEF 35%", result)
        self.assertIn("ARNI", result)

    def test_json_with_full_mock(self):
        result = self.exporter.export_json(
            self.mock_response, patient_id="P-HF-001"
        )
        self.assertEqual(result["patient_id"], "P-HF-001")
        self.assertIn("summary", result["data"])

    def test_markdown_guideline_references(self):
        result = self.exporter.export_markdown(self.mock_response)
        self.assertIn("2022 AHA/ACC/HFSA", result)

    def test_json_data_integrity(self):
        result = self.exporter.export_json(self.mock_response)
        data = result["data"]
        self.assertEqual(len(data["findings"]), 3)
        self.assertEqual(len(data["recommendations"]), 3)


if __name__ == "__main__":
    unittest.main()
