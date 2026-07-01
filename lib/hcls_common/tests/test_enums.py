"""Tests for hcls_common.enums — domain enumerations."""

import pytest

from hcls_common.enums import (
    PipelineStage,
    ClinicalSignificance,
    VariantImpact,
    TherapeuticArea,
    DrugLikenessFilter,
    DockingMethod,
    LLMProvider,
    EmbeddingProvider,
    NIMService,
    AlertSeverity,
)


class TestPipelineStage:
    def test_all_stages_exist(self):
        assert PipelineStage.GENOMICS.value == "genomics"
        assert PipelineStage.RAG_CHAT.value == "rag_chat"
        assert PipelineStage.DRUG_DISCOVERY.value == "drug_discovery"

    def test_display_name(self):
        assert PipelineStage.GENOMICS.display_name == "Genomics"
        assert PipelineStage.RAG_CHAT.display_name == "Rag Chat"
        assert PipelineStage.DRUG_DISCOVERY.display_name == "Drug Discovery"

    def test_string_enum(self):
        assert isinstance(PipelineStage.GENOMICS, str)
        assert PipelineStage.GENOMICS == "genomics"


class TestClinicalSignificance:
    def test_actionable(self):
        assert ClinicalSignificance.PATHOGENIC.is_actionable is True
        assert ClinicalSignificance.LIKELY_PATHOGENIC.is_actionable is True

    def test_not_actionable(self):
        assert ClinicalSignificance.UNCERTAIN.is_actionable is False
        assert ClinicalSignificance.LIKELY_BENIGN.is_actionable is False
        assert ClinicalSignificance.BENIGN.is_actionable is False

    def test_display_name(self):
        assert ClinicalSignificance.PATHOGENIC.display_name == "Pathogenic"
        assert ClinicalSignificance.UNCERTAIN.display_name == "Uncertain Significance"


class TestVariantImpact:
    def test_rank_ordering(self):
        assert VariantImpact.HIGH.rank < VariantImpact.MODERATE.rank
        assert VariantImpact.MODERATE.rank < VariantImpact.LOW.rank
        assert VariantImpact.LOW.rank < VariantImpact.MODIFIER.rank

    def test_values(self):
        assert VariantImpact.HIGH.value == "HIGH"
        assert VariantImpact.MODIFIER.value == "MODIFIER"

    def test_sortable_by_rank(self):
        impacts = [VariantImpact.LOW, VariantImpact.HIGH, VariantImpact.MODERATE]
        sorted_impacts = sorted(impacts, key=lambda x: x.rank)
        assert sorted_impacts == [VariantImpact.HIGH, VariantImpact.MODERATE, VariantImpact.LOW]


class TestTherapeuticArea:
    def test_all_areas(self):
        areas = list(TherapeuticArea)
        assert len(areas) == 5
        assert TherapeuticArea.NEURODEGENERATION in areas
        assert TherapeuticArea.ONCOLOGY in areas

    def test_display_name(self):
        assert TherapeuticArea.RARE_DISEASE.display_name == "Rare Disease"


class TestDrugLikenessFilter:
    def test_descriptions_exist(self):
        for f in DrugLikenessFilter:
            assert len(f.description) > 10

    def test_lipinski_description(self):
        assert "Rule of Five" in DrugLikenessFilter.LIPINSKI.description


class TestAlertSeverity:
    def test_rank_ordering(self):
        assert AlertSeverity.CRITICAL.rank < AlertSeverity.HIGH.rank
        assert AlertSeverity.HIGH.rank < AlertSeverity.MEDIUM.rank
        assert AlertSeverity.INFO.rank == 4

    def test_all_severities_have_ranks(self):
        for sev in AlertSeverity:
            assert isinstance(sev.rank, int)


class TestLLMProvider:
    def test_providers(self):
        assert LLMProvider.ANTHROPIC.value == "anthropic"
        assert LLMProvider.OLLAMA.value == "ollama"


class TestEmbeddingProvider:
    def test_providers(self):
        assert EmbeddingProvider.LOCAL_BGE.value == "local_bge"


class TestNIMService:
    def test_services(self):
        assert NIMService.MOLMIM.value == "molmim"
        assert NIMService.DIFFDOCK.value == "diffdock"
