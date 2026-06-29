"""Tests for the Cardiology Intelligence Agent ingest system.

Validates IngestRecord creation, BaseIngestParser abstract interface,
record validation, and domain-specific parsers (PubMed, ClinicalTrials,
Guideline, Device) including their reference data dictionaries.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import unittest
from abc import ABC

from src.ingest.base import BaseIngestParser, IngestRecord


# =====================================================================
# IngestRecord creation
# =====================================================================


class TestIngestRecordCreation(unittest.TestCase):
    """Tests for IngestRecord dataclass creation."""

    def test_create_record(self):
        record = IngestRecord(
            text="Test content",
            metadata={"key": "value"},
            collection="cardio_literature",
            source="PubMed",
        )
        self.assertIsInstance(record, IngestRecord)

    def test_record_text(self):
        record = IngestRecord(
            text="Heart failure treatment",
            metadata={},
            collection="cardio_literature",
            source="PubMed",
        )
        self.assertEqual(record.text, "Heart failure treatment")

    def test_record_metadata(self):
        meta = {"pmid": "12345", "year": 2024}
        record = IngestRecord(
            text="Test", metadata=meta, collection="test", source="test"
        )
        self.assertEqual(record.metadata["pmid"], "12345")
        self.assertEqual(record.metadata["year"], 2024)

    def test_record_collection(self):
        record = IngestRecord(
            text="Test", metadata={}, collection="cardio_trials", source="CT"
        )
        self.assertEqual(record.collection, "cardio_trials")

    def test_record_source(self):
        record = IngestRecord(
            text="Test", metadata={}, collection="test", source="PubMed"
        )
        self.assertEqual(record.source, "PubMed")

    def test_record_source_id_default_none(self):
        record = IngestRecord(
            text="Test", metadata={}, collection="test", source="test"
        )
        self.assertIsNone(record.source_id)

    def test_record_source_id_set(self):
        record = IngestRecord(
            text="Test",
            metadata={},
            collection="test",
            source="test",
            source_id="PMID12345",
        )
        self.assertEqual(record.source_id, "PMID12345")

    def test_record_empty_metadata(self):
        record = IngestRecord(
            text="Test", metadata={}, collection="test", source="test"
        )
        self.assertEqual(len(record.metadata), 0)

    def test_record_complex_metadata(self):
        meta = {
            "title": "Study Title",
            "authors": ["Author A", "Author B"],
            "year": 2024,
            "mesh_terms": ["Heart Failure"],
        }
        record = IngestRecord(
            text="Abstract text", metadata=meta, collection="test", source="PM"
        )
        self.assertEqual(len(record.metadata["authors"]), 2)


# =====================================================================
# BaseIngestParser abstract methods
# =====================================================================


class TestBaseIngestParserAbstract(unittest.TestCase):
    """Tests for BaseIngestParser abstract interface."""

    def test_is_abstract(self):
        self.assertTrue(issubclass(BaseIngestParser, ABC))

    def test_cannot_instantiate_directly(self):
        with self.assertRaises(TypeError):
            BaseIngestParser("test_collection")

    def test_has_fetch_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "fetch"))

    def test_has_parse_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "parse"))

    def test_has_run_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "run"))

    def test_has_validate_record_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "validate_record"))

    def test_has_filter_valid_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "filter_valid"))

    def test_has_truncate_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "truncate"))

    def test_has_safe_str_method(self):
        self.assertTrue(hasattr(BaseIngestParser, "safe_str"))


# =====================================================================
# Concrete parser for testing base class methods
# =====================================================================


class _TestParser(BaseIngestParser):
    """Concrete implementation for testing base class behaviour."""

    def fetch(self, **kwargs):
        return [{"text": "test"}]

    def parse(self, raw_data):
        return [
            IngestRecord(
                text=r.get("text", ""),
                metadata=r,
                collection=self.collection,
                source="test",
            )
            for r in raw_data
        ]


class TestBaseIngestParserValidation(unittest.TestCase):
    """Tests for validate_record and filter_valid."""

    def setUp(self):
        self.parser = _TestParser("test_collection")

    def test_validate_valid_record(self):
        record = IngestRecord(
            text="Valid text", metadata={}, collection="test_collection", source="test"
        )
        self.assertTrue(self.parser.validate_record(record))

    def test_validate_empty_text_returns_false(self):
        record = IngestRecord(
            text="", metadata={}, collection="test_collection", source="test"
        )
        self.assertFalse(self.parser.validate_record(record))

    def test_validate_whitespace_text_returns_false(self):
        record = IngestRecord(
            text="   ", metadata={}, collection="test_collection", source="test"
        )
        self.assertFalse(self.parser.validate_record(record))

    def test_validate_none_text_returns_false(self):
        record = IngestRecord(
            text=None, metadata={}, collection="test_collection", source="test"
        )
        self.assertFalse(self.parser.validate_record(record))

    def test_validate_no_collection_returns_false(self):
        record = IngestRecord(
            text="Valid", metadata={}, collection="", source="test"
        )
        self.assertFalse(self.parser.validate_record(record))

    def test_filter_valid_keeps_valid(self):
        records = [
            IngestRecord(
                text="Valid", metadata={}, collection="test_collection", source="test"
            ),
        ]
        result = self.parser.filter_valid(records)
        self.assertEqual(len(result), 1)

    def test_filter_valid_removes_invalid(self):
        records = [
            IngestRecord(
                text="", metadata={}, collection="test_collection", source="test"
            ),
        ]
        result = self.parser.filter_valid(records)
        self.assertEqual(len(result), 0)

    def test_filter_valid_mixed(self):
        records = [
            IngestRecord(
                text="Valid 1", metadata={}, collection="test_collection", source="test"
            ),
            IngestRecord(
                text="", metadata={}, collection="test_collection", source="test"
            ),
            IngestRecord(
                text="Valid 2", metadata={}, collection="test_collection", source="test"
            ),
        ]
        result = self.parser.filter_valid(records)
        self.assertEqual(len(result), 2)

    def test_filter_valid_empty_list(self):
        result = self.parser.filter_valid([])
        self.assertEqual(len(result), 0)


# =====================================================================
# BaseIngestParser utility helpers
# =====================================================================


class TestTruncate(unittest.TestCase):
    """Tests for BaseIngestParser.truncate static method."""

    def test_short_text_unchanged(self):
        result = BaseIngestParser.truncate("Hello", 100)
        self.assertEqual(result, "Hello")

    def test_exact_length_unchanged(self):
        result = BaseIngestParser.truncate("12345", 5)
        self.assertEqual(result, "12345")

    def test_long_text_truncated(self):
        result = BaseIngestParser.truncate("Hello, World!", 8)
        self.assertEqual(len(result), 8)
        self.assertTrue(result.endswith("..."))

    def test_truncated_has_ellipsis(self):
        result = BaseIngestParser.truncate("A very long string", 10)
        self.assertIn("...", result)


class TestSafeStr(unittest.TestCase):
    """Tests for BaseIngestParser.safe_str static method."""

    def test_none_returns_default(self):
        result = BaseIngestParser.safe_str(None)
        self.assertEqual(result, "")

    def test_none_with_custom_default(self):
        result = BaseIngestParser.safe_str(None, "N/A")
        self.assertEqual(result, "N/A")

    def test_string_passthrough(self):
        result = BaseIngestParser.safe_str("hello")
        self.assertEqual(result, "hello")

    def test_int_converted(self):
        result = BaseIngestParser.safe_str(42)
        self.assertEqual(result, "42")

    def test_float_converted(self):
        result = BaseIngestParser.safe_str(3.14)
        self.assertEqual(result, "3.14")


# =====================================================================
# BaseIngestParser run orchestration
# =====================================================================


class TestBaseIngestParserRun(unittest.TestCase):
    """Tests for BaseIngestParser.run orchestration."""

    def test_run_returns_list(self):
        parser = _TestParser("test_collection")
        result = parser.run()
        self.assertIsInstance(result, list)

    def test_run_returns_valid_records(self):
        parser = _TestParser("test_collection")
        result = parser.run()
        self.assertGreater(len(result), 0)
        for record in result:
            self.assertIsInstance(record, IngestRecord)

    def test_parser_collection_attribute(self):
        parser = _TestParser("my_collection")
        self.assertEqual(parser.collection, "my_collection")


# =====================================================================
# PubMedCardioParser
# =====================================================================


class TestPubMedCardioParser(unittest.TestCase):
    """Tests for PubMedCardioParser reference data."""

    @classmethod
    def setUpClass(cls):
        from src.ingest.pubmed_parser import PubMedCardioParser
        cls.ParserClass = PubMedCardioParser

    def test_has_cardio_mesh_terms(self):
        self.assertTrue(hasattr(self.ParserClass, "CARDIO_MESH_TERMS"))

    def test_mesh_terms_is_list(self):
        self.assertIsInstance(self.ParserClass.CARDIO_MESH_TERMS, list)

    def test_mesh_terms_not_empty(self):
        self.assertGreater(len(self.ParserClass.CARDIO_MESH_TERMS), 0)

    def test_mesh_terms_contains_cardiovascular(self):
        self.assertIn("Cardiovascular Diseases", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_contains_heart_failure(self):
        self.assertIn("Heart Failure", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_contains_cad(self):
        self.assertIn("Coronary Artery Disease", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_contains_af(self):
        self.assertIn("Atrial Fibrillation", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_contains_mi(self):
        self.assertIn("Myocardial Infarction", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_contains_cardiomyopathies(self):
        self.assertIn("Cardiomyopathies", self.ParserClass.CARDIO_MESH_TERMS)

    def test_mesh_terms_at_least_10(self):
        self.assertGreaterEqual(len(self.ParserClass.CARDIO_MESH_TERMS), 10)

    def test_parser_collection_name(self):
        parser = self.ParserClass()
        self.assertEqual(parser.collection, "cardio_literature")

    def test_is_base_ingest_parser(self):
        self.assertTrue(issubclass(self.ParserClass, BaseIngestParser))


# =====================================================================
# ClinicalTrialsCardioParser
# =====================================================================


class TestClinicalTrialsCardioParser(unittest.TestCase):
    """Tests for ClinicalTrialsCardioParser reference data."""

    @classmethod
    def setUpClass(cls):
        from src.ingest.clinical_trials_parser import ClinicalTrialsCardioParser
        cls.ParserClass = ClinicalTrialsCardioParser

    def test_has_landmark_trials(self):
        self.assertTrue(hasattr(self.ParserClass, "LANDMARK_TRIALS"))

    def test_landmark_trials_is_dict(self):
        self.assertIsInstance(self.ParserClass.LANDMARK_TRIALS, dict)

    def test_landmark_trials_not_empty(self):
        self.assertGreater(len(self.ParserClass.LANDMARK_TRIALS), 0)

    def test_landmark_trials_contains_paradigm_hf(self):
        self.assertIn("PARADIGM-HF", self.ParserClass.LANDMARK_TRIALS)

    def test_landmark_trials_contains_dapa_hf(self):
        self.assertIn("DAPA-HF", self.ParserClass.LANDMARK_TRIALS)

    def test_landmark_trials_contains_emperor_reduced(self):
        self.assertIn("EMPEROR-Reduced", self.ParserClass.LANDMARK_TRIALS)

    def test_landmark_trials_contains_emperor_preserved(self):
        self.assertIn("EMPEROR-Preserved", self.ParserClass.LANDMARK_TRIALS)

    def test_landmark_trials_contains_deliver(self):
        self.assertIn("DELIVER", self.ParserClass.LANDMARK_TRIALS)

    def test_landmark_trials_contains_rales(self):
        self.assertIn("RALES", self.ParserClass.LANDMARK_TRIALS)

    def test_paradigm_hf_has_nct(self):
        trial = self.ParserClass.LANDMARK_TRIALS["PARADIGM-HF"]
        self.assertIn("nct", trial)
        self.assertTrue(trial["nct"].startswith("NCT"))

    def test_paradigm_hf_has_condition(self):
        trial = self.ParserClass.LANDMARK_TRIALS["PARADIGM-HF"]
        self.assertIn("condition", trial)

    def test_paradigm_hf_has_drug(self):
        trial = self.ParserClass.LANDMARK_TRIALS["PARADIGM-HF"]
        self.assertIn("drug", trial)

    def test_paradigm_hf_has_key_findings(self):
        trial = self.ParserClass.LANDMARK_TRIALS["PARADIGM-HF"]
        self.assertIn("key_findings", trial)
        self.assertTrue(len(trial["key_findings"]) > 0)

    def test_dapa_hf_drug_is_dapagliflozin(self):
        trial = self.ParserClass.LANDMARK_TRIALS["DAPA-HF"]
        self.assertEqual(trial["drug"], "Dapagliflozin")

    def test_landmark_trials_at_least_5(self):
        self.assertGreaterEqual(len(self.ParserClass.LANDMARK_TRIALS), 5)

    def test_parser_collection_name(self):
        parser = self.ParserClass()
        self.assertEqual(parser.collection, "cardio_trials")

    def test_is_base_ingest_parser(self):
        self.assertTrue(issubclass(self.ParserClass, BaseIngestParser))

    def test_has_cardio_conditions(self):
        self.assertTrue(hasattr(self.ParserClass, "CARDIO_CONDITIONS"))
        self.assertIsInstance(self.ParserClass.CARDIO_CONDITIONS, list)


# =====================================================================
# GuidelineParser
# =====================================================================


class TestGuidelineParser(unittest.TestCase):
    """Tests for GuidelineParser reference data."""

    @classmethod
    def setUpClass(cls):
        from src.ingest.guideline_parser import GuidelineParser
        cls.ParserClass = GuidelineParser

    def test_has_guidelines(self):
        self.assertTrue(hasattr(self.ParserClass, "GUIDELINES"))

    def test_guidelines_is_dict(self):
        self.assertIsInstance(self.ParserClass.GUIDELINES, dict)

    def test_guidelines_not_empty(self):
        self.assertGreater(len(self.ParserClass.GUIDELINES), 0)

    def test_guidelines_contains_hf(self):
        hf_keys = [k for k in self.ParserClass.GUIDELINES if "Heart Failure" in k]
        self.assertGreater(len(hf_keys), 0)

    def test_guideline_entry_has_year(self):
        first_key = next(iter(self.ParserClass.GUIDELINES))
        entry = self.ParserClass.GUIDELINES[first_key]
        self.assertIn("year", entry)

    def test_guideline_entry_has_society(self):
        first_key = next(iter(self.ParserClass.GUIDELINES))
        entry = self.ParserClass.GUIDELINES[first_key]
        self.assertIn("society", entry)

    def test_guideline_entry_has_sections(self):
        first_key = next(iter(self.ParserClass.GUIDELINES))
        entry = self.ParserClass.GUIDELINES[first_key]
        self.assertIn("sections", entry)

    def test_parser_collection_name(self):
        parser = self.ParserClass()
        self.assertEqual(parser.collection, "cardio_guidelines")

    def test_is_base_ingest_parser(self):
        self.assertTrue(issubclass(self.ParserClass, BaseIngestParser))


# =====================================================================
# DeviceParser
# =====================================================================


class TestDeviceParser(unittest.TestCase):
    """Tests for DeviceParser reference data."""

    @classmethod
    def setUpClass(cls):
        from src.ingest.device_parser import DeviceParser
        cls.ParserClass = DeviceParser

    def test_has_fda_ai_cardiac_devices(self):
        self.assertTrue(hasattr(self.ParserClass, "FDA_AI_CARDIAC_DEVICES"))

    def test_fda_devices_is_list(self):
        self.assertIsInstance(self.ParserClass.FDA_AI_CARDIAC_DEVICES, list)

    def test_fda_devices_not_empty(self):
        self.assertGreater(len(self.ParserClass.FDA_AI_CARDIAC_DEVICES), 0)

    def test_fda_device_entry_has_name(self):
        first = self.ParserClass.FDA_AI_CARDIAC_DEVICES[0]
        self.assertIn("name", first)

    def test_fda_device_entry_has_type(self):
        first = self.ParserClass.FDA_AI_CARDIAC_DEVICES[0]
        self.assertIn("type", first)

    def test_fda_device_entry_has_manufacturer(self):
        first = self.ParserClass.FDA_AI_CARDIAC_DEVICES[0]
        self.assertIn("manufacturer", first)

    def test_fda_device_entry_has_application(self):
        first = self.ParserClass.FDA_AI_CARDIAC_DEVICES[0]
        self.assertIn("application", first)

    def test_parser_collection_name(self):
        parser = self.ParserClass()
        self.assertEqual(parser.collection, "cardio_devices")

    def test_is_base_ingest_parser(self):
        self.assertTrue(issubclass(self.ParserClass, BaseIngestParser))


if __name__ == "__main__":
    unittest.main()
