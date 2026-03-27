"""Tests for src/ingest/ — all PGx data ingest parsers.

Covers BaseIngestPipeline, CPICGuidelineParser, PharmVarParser,
PharmGKBParser, FDALabelParser, PopulationFrequencyParser,
PubMedPGxParser, and ClinicalTrialsPGxParser.

All HTTP requests are mocked — no network calls are made.

Author: Adam Jones
Date: March 2026
"""

import json
import xml.etree.ElementTree as ET
from typing import Any, Dict, List
from unittest.mock import MagicMock, patch

import pytest

from src.ingest.base import BaseIngestPipeline, _truncate_utf8
from src.ingest.cpic_parser import CPICGuidelineParser
from src.ingest.pharmvar_parser import PharmVarParser
from src.ingest.pharmgkb_parser import PharmGKBParser
from src.ingest.fda_label_parser import FDALabelParser
from src.ingest.population_parser import PopulationFrequencyParser
from src.ingest.pubmed_parser import PubMedPGxParser
from src.ingest.clinical_trials_parser import ClinicalTrialsPGxParser

from src.models import (
    AlertLevel,
    ClinicalAction,
    ClinicalEvidence,
    CPICLevel,
    DrugGuideline,
    DrugInteraction,
    EvidenceLevel,
    FDALabel,
    GeneReference,
    GuidelineBody,
    InteractionType,
    PGxClinicalTrial,
    PopulationData,
)


# ═══════════════════════════════════════════════════════════════════════
# Shared fixtures
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def mock_collection_manager():
    """Create a mock PGxCollectionManager."""
    mgr = MagicMock()
    mgr.insert_batch = MagicMock(return_value=None)
    return mgr


@pytest.fixture
def mock_embedder():
    """Create a mock embedder that returns 384-dim zero vectors."""
    embedder = MagicMock()
    embedder.encode = MagicMock(side_effect=lambda texts: [[0.0] * 384 for _ in texts])
    return embedder


# ═══════════════════════════════════════════════════════════════════════
# _truncate_utf8
# ═══════════════════════════════════════════════════════════════════════

class TestTruncateUtf8:
    """Tests for the byte-aware UTF-8 truncation helper."""

    def test_short_ascii_unchanged(self):
        assert _truncate_utf8("hello", 100) == "hello"

    def test_exact_length_unchanged(self):
        text = "abc"
        assert _truncate_utf8(text, 3) == "abc"

    def test_truncates_long_ascii(self):
        text = "a" * 200
        result = _truncate_utf8(text, 100)
        assert len(result.encode("utf-8")) <= 100

    def test_truncates_multibyte_safely(self):
        """Multi-byte characters should not be split mid-byte."""
        text = "\u00e9" * 100  # each e-acute is 2 bytes
        result = _truncate_utf8(text, 50)
        assert len(result.encode("utf-8")) <= 50
        # Should not raise on re-encode
        result.encode("utf-8")

    def test_empty_string(self):
        assert _truncate_utf8("", 100) == ""


# ═══════════════════════════════════════════════════════════════════════
# BaseIngestPipeline
# ═══════════════════════════════════════════════════════════════════════

class TestBaseIngestPipeline:
    """Tests for the abstract base class."""

    def test_cannot_instantiate_directly(self):
        """BaseIngestPipeline is abstract and should not be instantiated."""
        with pytest.raises(TypeError):
            BaseIngestPipeline(MagicMock(), MagicMock())

    def test_concrete_subclass(self, mock_collection_manager, mock_embedder):
        """A concrete subclass implementing fetch/parse can be instantiated."""

        class DummyPipeline(BaseIngestPipeline):
            COLLECTION_NAME = "test_collection"

            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = DummyPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.COLLECTION_NAME == "test_collection"
        assert pipeline.fetch() == []
        assert pipeline.parse([]) == []

    def test_run_raises_without_collection_name(self, mock_collection_manager, mock_embedder):
        """run() should raise ValueError when no collection_name is set."""

        class NoColl(BaseIngestPipeline):
            COLLECTION_NAME = ""

            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = NoColl(mock_collection_manager, mock_embedder)
        with pytest.raises(ValueError, match="collection_name must be provided"):
            pipeline.run()

    def test_embed_and_store_inserts_batches(self, mock_collection_manager, mock_embedder):
        """embed_and_store should batch insert records and return count."""

        class Dummy(BaseIngestPipeline):
            COLLECTION_NAME = "test"

            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = Dummy(mock_collection_manager, mock_embedder)

        # Create mock records with to_embedding_text and model_dump
        records = []
        for i in range(5):
            rec = MagicMock()
            rec.to_embedding_text.return_value = f"text {i}"
            rec.model_dump.return_value = {"id": f"rec_{i}", "text_chunk": f"text {i}"}
            records.append(rec)

        count = pipeline.embed_and_store(records, "test", batch_size=2)
        assert count == 5
        assert mock_collection_manager.insert_batch.call_count == 3  # 2+2+1

    def test_embed_and_store_continues_on_error(self, mock_collection_manager, mock_embedder):
        """embed_and_store should skip failed batches and continue."""

        class Dummy(BaseIngestPipeline):
            COLLECTION_NAME = "test"

            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = Dummy(mock_collection_manager, mock_embedder)
        mock_collection_manager.insert_batch.side_effect = [None, Exception("Milvus error"), None]

        records = []
        for i in range(6):
            rec = MagicMock()
            rec.to_embedding_text.return_value = f"text {i}"
            rec.model_dump.return_value = {"id": f"rec_{i}", "text_chunk": f"text {i}"}
            records.append(rec)

        count = pipeline.embed_and_store(records, "test", batch_size=2)
        # First batch (2) + second fails (0) + third batch (2) = 4
        assert count == 4


# ═══════════════════════════════════════════════════════════════════════
# CPICGuidelineParser
# ═══════════════════════════════════════════════════════════════════════

class TestCPICGuidelineParser:
    """Tests for the CPIC guideline ingest parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return CPICGuidelineParser(mock_collection_manager, mock_embedder)

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_drug_guidelines"

    @patch("src.ingest.cpic_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return dict with pairs, recommendations, guidelines."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.side_effect = [
            [{"genesymbol": "CYP2D6", "drugname": "codeine", "cpiclevel": "A"}],
            [{"genesymbol": "CYP2D6", "drugname": "codeine", "drugrecommendation": "Avoid"}],
            [{"id": 1, "version": "2024"}],
        ]
        mock_get.return_value = mock_resp

        result = parser.fetch()
        assert len(result["pairs"]) == 1
        assert len(result["recommendations"]) == 1
        assert len(result["guidelines"]) == 1

    @patch("src.ingest.cpic_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        """fetch() should return empty lists on HTTP error."""
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("timeout")

        result = parser.fetch()
        assert result["pairs"] == []
        assert result["recommendations"] == []
        assert result["guidelines"] == []

    def test_parse_valid_recommendations(self, parser):
        """parse() should produce DrugGuideline records from valid data."""
        raw_data = {
            "pairs": [
                {
                    "genesymbol": "CYP2D6",
                    "drugname": "codeine",
                    "cpiclevel": "A",
                    "guidelineid": 1,
                    "lastmodified": "2024-01-15T00:00:00Z",
                },
            ],
            "recommendations": [
                {
                    "genesymbol": "CYP2D6",
                    "drugname": "codeine",
                    "lookupkey": {"CYP2D6": "Poor Metabolizer"},
                    "drugrecommendation": "Avoid codeine. Use alternative analgesic not metabolized by CYP2D6.",
                    "classification": "Strong",
                    "alternatedrugavailable": "morphine, oxycodone",
                },
            ],
            "guidelines": [
                {"id": 1, "version": "2024.1"},
            ],
        }

        records = parser.parse(raw_data)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, DrugGuideline)
        assert rec.gene == "CYP2D6"
        assert rec.drug == "codeine"
        assert "Poor Metabolizer" in rec.phenotype
        assert rec.guideline_body == GuidelineBody.CPIC
        assert rec.cpic_level == CPICLevel.A
        assert rec.clinical_action == ClinicalAction.AVOID
        assert rec.alert_level == AlertLevel.CRITICAL

    def test_parse_empty_data(self, parser):
        """parse() should return empty list for empty input."""
        result = parser.parse({"pairs": [], "recommendations": [], "guidelines": []})
        assert result == []

    def test_parse_skips_malformed_records(self, parser):
        """parse() should skip records that cause exceptions."""
        raw_data = {
            "pairs": [],
            "recommendations": [
                {"genesymbol": None, "drugname": None, "drugrecommendation": None},
            ],
            "guidelines": [],
        }
        # Should not raise — malformed records are skipped
        records = parser.parse(raw_data)
        # May produce a record with empty fields or skip it
        assert isinstance(records, list)

    def test_classify_action_keywords(self):
        """_classify_action should match known keywords."""
        assert CPICGuidelineParser._classify_action("Avoid this drug") == ClinicalAction.AVOID
        assert CPICGuidelineParser._classify_action("Contraindicated in PM") == ClinicalAction.CONTRAINDICATED
        assert CPICGuidelineParser._classify_action("Reduce dose by 50%") == ClinicalAction.DOSE_ADJUST
        assert CPICGuidelineParser._classify_action("Use alternative agent") == ClinicalAction.ALTERNATIVE
        assert CPICGuidelineParser._classify_action("Standard dosing") == ClinicalAction.STANDARD

    def test_classify_alert_keywords(self):
        """_classify_alert should match known keywords."""
        assert CPICGuidelineParser._classify_alert("Do not use this drug") == AlertLevel.CRITICAL
        assert CPICGuidelineParser._classify_alert("Monitor closely") == AlertLevel.WARNING
        assert CPICGuidelineParser._classify_alert("No known concerns") == AlertLevel.INFO


# ═══════════════════════════════════════════════════════════════════════
# PharmVarParser
# ═══════════════════════════════════════════════════════════════════════

class TestPharmVarParser:
    """Tests for the PharmVar star allele definition parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return PharmVarParser(mock_collection_manager, mock_embedder, genes=["CYP2D6"])

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_gene_reference"

    @patch("src.ingest.pharmvar_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return allele data with _gene annotations."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = [
            {
                "alleleName": "*4",
                "function": "no function",
                "activityScore": 0.0,
                "pharmvarId": "PV00001",
                "variants": [{"rsid": "rs3892097"}],
            },
        ]
        mock_get.return_value = mock_resp

        result = parser.fetch()
        assert len(result) == 1
        assert result[0]["_gene"] == "CYP2D6"

    @patch("src.ingest.pharmvar_parser.requests.get")
    def test_fetch_handles_dict_response(self, mock_get, parser):
        """fetch() should handle responses wrapped in a 'data' key."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = {
            "data": [
                {"alleleName": "*1", "function": "normal function", "pharmvarId": "PV00002"},
            ],
        }
        mock_get.return_value = mock_resp

        result = parser.fetch()
        assert len(result) == 1
        assert result[0]["_gene"] == "CYP2D6"

    @patch("src.ingest.pharmvar_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        """fetch() should return empty list on HTTP error."""
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("refused")

        result = parser.fetch()
        assert result == []

    def test_parse_valid_alleles(self, parser):
        """parse() should produce GeneReference records from allele data."""
        raw_data = [
            {
                "_gene": "CYP2D6",
                "alleleName": "*4",
                "function": "no function",
                "activityScore": 0.0,
                "pharmvarId": "PV00001",
                "variants": [{"rsid": "rs3892097"}, {"rsid": "rs1065852"}],
                "frequencyGlobal": 0.12,
                "frequencyEuropean": 0.20,
                "frequencyAfrican": 0.06,
            },
        ]

        records = parser.parse(raw_data)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, GeneReference)
        assert rec.gene == "CYP2D6"
        assert rec.star_allele == "*4"
        assert rec.function_status == "no function"
        assert rec.activity_score == 0.0
        assert "rs3892097" in rec.defining_variants
        assert rec.allele_frequency_global == 0.12
        assert rec.source == "PharmVar"

    def test_parse_adds_star_prefix(self, parser):
        """parse() should add * prefix to allele names missing it."""
        raw_data = [{"_gene": "CYP2D6", "alleleName": "4", "pharmvarId": "PV1"}]
        records = parser.parse(raw_data)
        assert records[0].star_allele == "*4"

    def test_parse_empty_data(self, parser):
        records = parser.parse([])
        assert records == []

    def test_safe_float_valid(self):
        assert PharmVarParser._safe_float(0.5) == 0.5
        assert PharmVarParser._safe_float("0.3") == 0.3
        assert PharmVarParser._safe_float(0.0) == 0.0
        assert PharmVarParser._safe_float(1.0) == 1.0

    def test_safe_float_invalid(self):
        assert PharmVarParser._safe_float(None) is None
        assert PharmVarParser._safe_float("not_a_number") is None
        assert PharmVarParser._safe_float(1.5) is None  # Out of 0-1 range
        assert PharmVarParser._safe_float(-0.1) is None


# ═══════════════════════════════════════════════════════════════════════
# PharmGKBParser
# ═══════════════════════════════════════════════════════════════════════

class TestPharmGKBParser:
    """Tests for the PharmGKB clinical annotation parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return PharmGKBParser(mock_collection_manager, mock_embedder)

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_drug_interactions"

    @patch("src.ingest.pharmgkb_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return clinical and variant annotations."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.side_effect = [
            {
                "data": [
                    {
                        "id": "CA001",
                        "relatedGenes": [{"symbol": "CYP2D6"}],
                        "relatedChemicals": [{"name": "codeine"}],
                        "evidenceLevel": "1A",
                        "text": "CYP2D6 poor metabolizers have reduced morphine formation.",
                    },
                ],
            },
            {
                "data": [
                    {
                        "id": "VA001",
                        "variant": {"name": "rs3892097"},
                        "relatedGenes": [{"symbol": "CYP2D6"}],
                        "relatedChemicals": [{"name": "codeine"}],
                        "sentence": "rs3892097 is associated with reduced codeine metabolism.",
                    },
                ],
            },
        ]
        mock_get.return_value = mock_resp

        result = parser.fetch()
        assert len(result["clinical_annotations"]) == 1
        assert len(result["variant_annotations"]) == 1

    @patch("src.ingest.pharmgkb_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        """fetch() should return empty lists on HTTP errors."""
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("refused")

        result = parser.fetch()
        assert result["clinical_annotations"] == []
        assert result["variant_annotations"] == []

    def test_parse_clinical_annotation(self, parser):
        """parse() should produce DrugInteraction from clinical annotations."""
        raw_data = {
            "clinical_annotations": [
                {
                    "id": "CA12345",
                    "relatedGenes": [{"symbol": "CYP2D6"}],
                    "relatedChemicals": [{"name": "codeine"}],
                    "relatedVariants": [{"name": "rs3892097"}],
                    "evidenceLevel": "1A",
                    "text": "CYP2D6 poor metabolizers produce less morphine from codeine.",
                    "phenotypeCategory": "efficacy",
                    "significance": "Reduced analgesic effect",
                },
            ],
            "variant_annotations": [],
        }

        records = parser.parse(raw_data)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, DrugInteraction)
        assert rec.gene == "CYP2D6"
        assert rec.drug == "codeine"
        assert rec.variant_rsid == "rs3892097"
        assert rec.evidence_level == EvidenceLevel.LEVEL_1A
        assert rec.interaction_type == InteractionType.EFFICACY

    def test_parse_variant_annotation(self, parser):
        """parse() should produce ClinicalEvidence from variant annotations."""
        raw_data = {
            "clinical_annotations": [],
            "variant_annotations": [
                {
                    "id": "VA54321",
                    "variant": {"name": "rs3892097"},
                    "relatedGenes": [{"symbol": "CYP2D6"}],
                    "relatedChemicals": [{"name": "codeine"}],
                    "sentence": "rs3892097 is the defining variant of CYP2D6*4.",
                    "studyType": "cohort",
                    "literatures": [{"pmid": "12345678"}],
                    "significance": "functional variant",
                },
            ],
        }

        records = parser.parse(raw_data)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, ClinicalEvidence)
        assert rec.gene == "CYP2D6"
        assert rec.drug == "codeine"
        assert rec.pmid == "12345678"
        assert rec.source == "PharmGKB"

    def test_parse_skips_annotation_without_gene_or_drug(self, parser):
        """parse() should skip clinical annotations missing gene or drug."""
        raw_data = {
            "clinical_annotations": [
                {"id": "CA999", "relatedGenes": [], "relatedChemicals": [], "text": "orphan"},
            ],
            "variant_annotations": [],
        }

        records = parser.parse(raw_data)
        assert len(records) == 0

    def test_classify_interaction_type(self):
        assert PharmGKBParser._classify_interaction_type("toxicity observed", "") == InteractionType.TOXICITY
        assert PharmGKBParser._classify_interaction_type("response to drug", "efficacy") == InteractionType.EFFICACY
        assert PharmGKBParser._classify_interaction_type("pharmacodynamic effect", "") == InteractionType.PD
        assert PharmGKBParser._classify_interaction_type("no keywords here", "") == InteractionType.PK

    def test_parse_empty_data(self, parser):
        records = parser.parse({"clinical_annotations": [], "variant_annotations": []})
        assert records == []


# ═══════════════════════════════════════════════════════════════════════
# FDALabelParser
# ═══════════════════════════════════════════════════════════════════════

class TestFDALabelParser:
    """Tests for the FDA pharmacogenomic drug label parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return FDALabelParser(mock_collection_manager, mock_embedder)

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_fda_labels"

    @patch("src.ingest.fda_label_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return deduplicated FDA labels."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = {
            "results": [
                {
                    "id": "FDA001",
                    "openfda": {"brand_name": ["PLAVIX"], "generic_name": ["clopidogrel"]},
                    "warnings_and_precautions": ["CYP2C19 poor metabolizers..."],
                },
            ],
        }
        mock_get.return_value = mock_resp

        result = parser.fetch(genes=["CYP2C19"], max_per_gene=10)
        assert len(result) == 1
        assert result[0]["_query_gene"] == "CYP2C19"

    @patch("src.ingest.fda_label_parser.requests.get")
    def test_fetch_deduplicates(self, mock_get, parser):
        """fetch() should not include duplicate labels."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = {
            "results": [
                {"id": "FDA001", "openfda": {"brand_name": ["PLAVIX"]}},
                {"id": "FDA001", "openfda": {"brand_name": ["PLAVIX"]}},
            ],
        }
        mock_get.return_value = mock_resp

        result = parser.fetch(genes=["CYP2C19"], max_per_gene=10)
        assert len(result) == 1

    @patch("src.ingest.fda_label_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        """fetch() should return empty list on HTTP error."""
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("timeout")

        result = parser.fetch(genes=["CYP2D6"])
        assert result == []

    def test_parse_label_with_pgx_section(self, parser):
        """parse() should extract sections mentioning the query gene."""
        raw_data = [
            {
                "id": "FDA001",
                "_query_gene": "CYP2C19",
                "openfda": {"brand_name": ["PLAVIX"], "generic_name": ["clopidogrel"]},
                "warnings_and_precautions": [
                    "CYP2C19 poor metabolizers have reduced active metabolite formation."
                ],
                "dosage_and_administration": [
                    "Standard dose for normal CYP2C19 metabolizers."
                ],
                "clinical_pharmacology": [
                    "Unrelated pharmacology text without gene mention."
                ],
                "effective_time": "20240115",
            },
        ]

        records = parser.parse(raw_data)
        # Should find sections that mention CYP2C19
        assert len(records) >= 1
        rec = records[0]
        assert isinstance(rec, FDALabel)
        assert rec.drug == "PLAVIX"
        assert rec.gene == "CYP2C19"
        assert rec.last_updated == "2024-01-15"

    def test_parse_skips_sections_without_gene(self, parser):
        """parse() should skip label sections that don't mention the query gene."""
        raw_data = [
            {
                "id": "FDA002",
                "_query_gene": "CYP2D6",
                "openfda": {"brand_name": ["TYLENOL"]},
                "clinical_pharmacology": ["No pharmacogenomic information."],
            },
        ]

        records = parser.parse(raw_data)
        assert len(records) == 0

    def test_parse_empty_data(self, parser):
        records = parser.parse([])
        assert records == []

    def test_classify_label_type(self):
        assert FDALabelParser._classify_label_type("Testing required before use") == "testing required"
        assert FDALabelParser._classify_label_type("Consider dose adjustment") == "actionable PGx"
        assert FDALabelParser._classify_label_type("Polymorphism may affect levels") == "informative PGx"
        assert FDALabelParser._classify_label_type("No keywords") == "informative PGx"

    def test_classify_requirement(self):
        assert FDALabelParser._classify_requirement("Testing is required") == "required"
        assert FDALabelParser._classify_requirement("Testing is recommended") == "recommended"
        assert FDALabelParser._classify_requirement("May influence response") == "informational"
        assert FDALabelParser._classify_requirement("No keywords") == "informational"


# ═══════════════════════════════════════════════════════════════════════
# PopulationFrequencyParser
# ═══════════════════════════════════════════════════════════════════════

class TestPopulationFrequencyParser:
    """Tests for the population allele frequency parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return PopulationFrequencyParser(
            mock_collection_manager, mock_embedder, genes=["CYP2D6"]
        )

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_population_data"

    @patch("src.ingest.population_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return frequency data indexed by gene."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = {
            "data": [
                {
                    "alleleName": "*4",
                    "populations": [
                        {"population": "European", "frequency": 0.20, "sampleSize": 1000},
                        {"population": "African", "frequency": 0.06, "sampleSize": 500},
                    ],
                },
            ],
        }
        mock_get.return_value = mock_resp

        result = parser.fetch()
        assert "CYP2D6" in result["pharmgkb_frequencies"]

    @patch("src.ingest.population_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("refused")

        result = parser.fetch()
        assert result["pharmgkb_frequencies"] == {}

    def test_parse_list_populations(self, parser):
        """parse() should handle list-format population frequencies."""
        raw_data = {
            "pharmgkb_frequencies": {
                "CYP2D6": [
                    {
                        "alleleName": "*4",
                        "populations": [
                            {"population": "European", "frequency": 0.20, "sampleSize": 1000},
                            {"population": "East Asian", "frequency": 0.01, "sampleSize": 800},
                        ],
                        "source": "PharmGKB",
                    },
                ],
            },
        }

        records = parser.parse(raw_data)
        assert len(records) == 2
        for rec in records:
            assert isinstance(rec, PopulationData)
            assert rec.gene == "CYP2D6"
            assert rec.star_allele == "*4"
        # Check populations
        pops = {r.population for r in records}
        assert "European" in pops
        assert "East Asian" in pops

    def test_parse_dict_populations(self, parser):
        """parse() should handle dict-format population frequencies."""
        raw_data = {
            "pharmgkb_frequencies": {
                "CYP2D6": [
                    {
                        "alleleName": "*1",
                        "populations": {"European": 0.35, "African": 0.50},
                        "source": "gnomAD",
                    },
                ],
            },
        }

        records = parser.parse(raw_data)
        assert len(records) == 2

    def test_parse_global_frequency_fallback(self, parser):
        """parse() should fall back to global frequency when no population breakdown."""
        raw_data = {
            "pharmgkb_frequencies": {
                "CYP2D6": [
                    {
                        "alleleName": "*10",
                        "populations": "unavailable",
                        "frequency": 0.15,
                    },
                ],
            },
        }

        records = parser.parse(raw_data)
        assert len(records) == 1
        assert records[0].population == "Global"
        assert records[0].allele_frequency == 0.15

    def test_parse_rejects_invalid_frequencies(self, parser):
        """parse() should reject frequencies outside 0-1 range."""
        raw_data = {
            "pharmgkb_frequencies": {
                "CYP2D6": [
                    {
                        "alleleName": "*99",
                        "populations": {"European": 1.5},
                    },
                ],
            },
        }

        records = parser.parse(raw_data)
        assert len(records) == 0

    def test_parse_empty_data(self, parser):
        records = parser.parse({"pharmgkb_frequencies": {}})
        assert records == []


# ═══════════════════════════════════════════════════════════════════════
# PubMedPGxParser
# ═══════════════════════════════════════════════════════════════════════

class TestPubMedPGxParser:
    """Tests for the PubMed PGx literature parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return PubMedPGxParser(mock_collection_manager, mock_embedder, api_key="test_key")

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_clinical_evidence"

    def test_rate_limit_interval_with_api_key(self, parser):
        """With an API key, rate limit should be 10 req/s."""
        assert parser._min_interval == pytest.approx(0.1, abs=0.01)

    def test_rate_limit_interval_without_api_key(self, mock_collection_manager, mock_embedder):
        """Without an API key, rate limit should be 3 req/s."""
        p = PubMedPGxParser(mock_collection_manager, mock_embedder)
        assert p._min_interval == pytest.approx(1.0 / 3, abs=0.01)

    @patch("src.ingest.pubmed_parser.requests.get")
    @patch("src.ingest.pubmed_parser.time.sleep")
    def test_fetch_success(self, mock_sleep, mock_get, parser):
        """fetch() should retrieve PMIDs via esearch then articles via efetch."""
        esearch_resp = MagicMock()
        esearch_resp.raise_for_status = MagicMock()
        esearch_resp.json.return_value = {
            "esearchresult": {"count": "1", "idlist": ["39000001"]},
        }

        xml_content = b"""<?xml version="1.0"?>
        <PubmedArticleSet>
            <PubmedArticle>
                <MedlineCitation>
                    <PMID>39000001</PMID>
                    <Article>
                        <ArticleTitle>CYP2D6 and codeine pharmacogenomics</ArticleTitle>
                        <Abstract>
                            <AbstractText>This meta-analysis of CYP2D6 codeine studies...</AbstractText>
                        </Abstract>
                        <Journal><Title>Clin Pharmacol Ther</Title></Journal>
                    </Article>
                    <MeshHeadingList>
                        <MeshHeading><DescriptorName>Pharmacogenomics</DescriptorName></MeshHeading>
                    </MeshHeadingList>
                </MedlineCitation>
                <PubmedData>
                    <History>
                        <PubMedPubDate PubStatus="pubmed">
                            <Year>2025</Year>
                        </PubMedPubDate>
                    </History>
                </PubmedData>
            </PubmedArticle>
        </PubmedArticleSet>"""

        efetch_resp = MagicMock()
        efetch_resp.raise_for_status = MagicMock()
        efetch_resp.content = xml_content

        mock_get.side_effect = [esearch_resp, efetch_resp]

        articles = parser.fetch(query="CYP2D6 codeine", max_results=10)
        assert len(articles) == 1
        assert articles[0]["pmid"] == "39000001"
        assert "CYP2D6" in articles[0]["title"]

    @patch("src.ingest.pubmed_parser.requests.get")
    @patch("src.ingest.pubmed_parser.time.sleep")
    def test_fetch_returns_empty_on_esearch_error(self, mock_sleep, mock_get, parser):
        """fetch() should return empty list on esearch failure."""
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("timeout")

        result = parser.fetch()
        assert result == []

    def test_parse_articles(self, parser):
        """parse() should produce ClinicalEvidence records."""
        articles = [
            {
                "pmid": "39000001",
                "title": "CYP2D6 and codeine pharmacogenomics meta-analysis",
                "abstract": "This meta-analysis evaluates CYP2D6 genotype-guided codeine prescribing.",
                "journal": "Clin Pharmacol Ther",
                "year": "2025",
                "mesh_terms": ["Pharmacogenomics"],
            },
        ]

        records = parser.parse(articles)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, ClinicalEvidence)
        assert rec.pmid == "39000001"
        assert rec.gene == "CYP2D6"
        assert rec.drug == "codeine"
        assert rec.study_type == "meta-analysis"
        assert rec.year == 2025
        assert rec.source == "PubMed"

    def test_parse_empty_data(self, parser):
        records = parser.parse([])
        assert records == []

    def test_extract_gene(self):
        assert PubMedPGxParser._extract_gene("CYP2D6 poor metabolizers") == "CYP2D6"
        assert PubMedPGxParser._extract_gene("DPYD deficiency and 5-FU") == "DPYD"
        assert PubMedPGxParser._extract_gene("No gene mentioned") == ""
        assert PubMedPGxParser._extract_gene("") == ""

    def test_extract_drug(self):
        assert PubMedPGxParser._extract_drug("codeine metabolism") == "codeine"
        assert PubMedPGxParser._extract_drug("warfarin dosing") == "warfarin"
        assert PubMedPGxParser._extract_drug("No drug mentioned") == ""
        assert PubMedPGxParser._extract_drug("") == ""

    def test_classify_study_type(self):
        assert PubMedPGxParser._classify_study_type("meta-analysis of 20 studies") == "meta-analysis"
        assert PubMedPGxParser._classify_study_type("randomized controlled trial") == "RCT"
        assert PubMedPGxParser._classify_study_type("retrospective cohort study") == "cohort"
        assert PubMedPGxParser._classify_study_type("clinical implementation of PGx") == "implementation"
        assert PubMedPGxParser._classify_study_type("a standard study") == "primary research"
        assert PubMedPGxParser._classify_study_type("") == ""

    def test_parse_xml_article(self):
        """_parse_xml_article should extract fields from PubMed XML."""
        xml_str = """<PubmedArticle>
            <MedlineCitation>
                <PMID>12345</PMID>
                <Article>
                    <ArticleTitle>Test Title</ArticleTitle>
                    <Abstract>
                        <AbstractText Label="BACKGROUND">Background text.</AbstractText>
                        <AbstractText>Main text.</AbstractText>
                    </Abstract>
                    <Journal><Title>Test Journal</Title></Journal>
                </Article>
            </MedlineCitation>
            <PubmedData>
                <History>
                    <PubMedPubDate PubStatus="pubmed">
                        <Year>2024</Year>
                    </PubMedPubDate>
                </History>
            </PubmedData>
        </PubmedArticle>"""
        elem = ET.fromstring(xml_str)
        result = PubMedPGxParser._parse_xml_article(elem)
        assert result is not None
        assert result["pmid"] == "12345"
        assert result["title"] == "Test Title"
        assert "BACKGROUND: Background text." in result["abstract"]
        assert result["journal"] == "Test Journal"


# ═══════════════════════════════════════════════════════════════════════
# ClinicalTrialsPGxParser
# ═══════════════════════════════════════════════════════════════════════

class TestClinicalTrialsPGxParser:
    """Tests for the ClinicalTrials.gov PGx trial parser."""

    @pytest.fixture
    def parser(self, mock_collection_manager, mock_embedder):
        return ClinicalTrialsPGxParser(mock_collection_manager, mock_embedder)

    def test_collection_name(self, parser):
        assert parser.COLLECTION_NAME == "pgx_clinical_trials"

    @patch("src.ingest.clinical_trials_parser.requests.get")
    def test_fetch_success(self, mock_get, parser):
        """fetch() should return study dicts from ClinicalTrials.gov API."""
        mock_resp = MagicMock()
        mock_resp.raise_for_status = MagicMock()
        mock_resp.json.return_value = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00001234",
                            "briefTitle": "CYP2D6-guided codeine dosing trial",
                        },
                        "statusModule": {
                            "overallStatus": "RECRUITING",
                            "startDateStruct": {"date": "2024-06-01"},
                        },
                        "designModule": {
                            "phases": ["PHASE3"],
                            "enrollmentInfo": {"count": 200},
                        },
                        "descriptionModule": {
                            "briefSummary": "A randomized trial of CYP2D6 genotype-guided codeine prescribing.",
                        },
                        "conditionsModule": {},
                        "armsInterventionsModule": {
                            "interventions": [{"name": "CYP2D6 genotyping"}],
                        },
                        "outcomesModule": {
                            "primaryOutcomes": [{"measure": "Pain relief score"}],
                        },
                    },
                },
            ],
            "nextPageToken": None,
        }
        mock_get.return_value = mock_resp

        result = parser.fetch(query="CYP2D6", max_results=10)
        assert len(result) == 1

    @patch("src.ingest.clinical_trials_parser.requests.get")
    def test_fetch_pagination(self, mock_get, parser):
        """fetch() should paginate using nextPageToken."""
        resp1 = MagicMock()
        resp1.raise_for_status = MagicMock()
        resp1.json.return_value = {
            "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT001"}}}],
            "nextPageToken": "token2",
        }

        resp2 = MagicMock()
        resp2.raise_for_status = MagicMock()
        resp2.json.return_value = {
            "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT002"}}}],
            "nextPageToken": None,
        }

        mock_get.side_effect = [resp1, resp2]

        result = parser.fetch(max_results=200)
        assert len(result) == 2

    @patch("src.ingest.clinical_trials_parser.requests.get")
    def test_fetch_handles_http_error(self, mock_get, parser):
        import requests as req
        mock_get.side_effect = req.exceptions.ConnectionError("timeout")

        result = parser.fetch()
        assert result == []

    def test_parse_valid_study(self, parser):
        """parse() should produce PGxClinicalTrial records from valid study data."""
        raw_data = [
            {
                "protocolSection": {
                    "identificationModule": {
                        "nctId": "NCT00001234",
                        "briefTitle": "CYP2D6-guided codeine dosing trial",
                        "officialTitle": "Genotype-Guided Codeine Prescribing",
                    },
                    "statusModule": {
                        "overallStatus": "RECRUITING",
                        "startDateStruct": {"date": "2024-06-01"},
                    },
                    "designModule": {
                        "phases": ["PHASE3"],
                        "enrollmentInfo": {"count": 200},
                    },
                    "descriptionModule": {
                        "briefSummary": "Testing CYP2D6 genotype-guided codeine prescribing.",
                    },
                    "conditionsModule": {},
                    "armsInterventionsModule": {
                        "interventions": [{"name": "CYP2D6 genotyping"}],
                    },
                    "outcomesModule": {
                        "primaryOutcomes": [{"measure": "Pain control"}],
                    },
                },
            },
        ]

        records = parser.parse(raw_data)
        assert len(records) == 1
        rec = records[0]
        assert isinstance(rec, PGxClinicalTrial)
        assert rec.nct_id == "NCT00001234"
        assert rec.gene == "CYP2D6"
        assert rec.drug == "codeine"
        assert rec.phase == "Phase 3"
        assert rec.status == "RECRUITING"
        assert rec.enrollment == 200
        assert rec.start_year == 2024

    def test_parse_empty_data(self, parser):
        records = parser.parse([])
        assert records == []

    def test_extract_gene(self):
        assert ClinicalTrialsPGxParser._extract_gene("CYP2D6 genotype-guided") == "CYP2D6"
        assert ClinicalTrialsPGxParser._extract_gene("DPYD and fluorouracil") == "DPYD"
        assert ClinicalTrialsPGxParser._extract_gene("No genes mentioned") == ""
        assert ClinicalTrialsPGxParser._extract_gene("") == ""

    def test_extract_drug(self):
        assert ClinicalTrialsPGxParser._extract_drug("warfarin dosing study") == "warfarin"
        assert ClinicalTrialsPGxParser._extract_drug("codeine and CYP2D6") == "codeine"
        assert ClinicalTrialsPGxParser._extract_drug("No drugs") == ""
        assert ClinicalTrialsPGxParser._extract_drug("") == ""


# ═══════════════════════════════════════════════════════════════════════
# Record model validation (cross-parser)
# ═══════════════════════════════════════════════════════════════════════

class TestRecordFormats:
    """Verify that parsed records have all expected fields populated."""

    def test_drug_guideline_has_embedding_text(self):
        rec = DrugGuideline(
            id="test1", gene="CYP2D6", drug="codeine", phenotype="PM",
            recommendation="Avoid codeine", text_chunk="CYP2D6 codeine PM",
        )
        text = rec.to_embedding_text()
        assert "CYP2D6" in text
        assert "codeine" in text

    def test_gene_reference_has_embedding_text(self):
        rec = GeneReference(
            id="test2", gene="CYP2D6", star_allele="*4",
            text_chunk="CYP2D6 *4 no function", source="PharmVar",
        )
        text = rec.to_embedding_text()
        assert "CYP2D6" in text
        assert "*4" in text

    def test_drug_interaction_has_embedding_text(self):
        rec = DrugInteraction(
            id="test3", drug="codeine", gene="CYP2D6",
            effect_description="Reduced morphine formation",
            text_chunk="CYP2D6 codeine interaction",
        )
        text = rec.to_embedding_text()
        assert "codeine" in text

    def test_clinical_evidence_has_embedding_text(self):
        rec = ClinicalEvidence(
            id="test4", title="CYP2D6 and codeine",
            text_chunk="Study of CYP2D6 codeine",
        )
        text = rec.to_embedding_text()
        assert "CYP2D6" in text

    def test_fda_label_has_embedding_text(self):
        rec = FDALabel(
            id="test5", drug="clopidogrel", gene="CYP2C19",
            text_chunk="CYP2C19 clopidogrel FDA",
        )
        text = rec.to_embedding_text()
        assert "clopidogrel" in text

    def test_population_data_has_embedding_text(self):
        rec = PopulationData(
            id="test6", gene="CYP2D6", star_allele="*4",
            population="European", allele_frequency=0.20,
            text_chunk="CYP2D6 *4 European 0.20",
        )
        text = rec.to_embedding_text()
        assert "European" in text

    def test_clinical_trial_has_embedding_text(self):
        rec = PGxClinicalTrial(
            id="test7", title="CYP2D6 codeine trial",
            text_summary="CYP2D6 genotype-guided codeine study",
        )
        text = rec.to_embedding_text()
        assert "CYP2D6" in text
