"""Tests for all diagnostic workflow types.

Covers all 13 DiagnosticWorkflowType members and their properties.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import DiagnosticWorkflowType


class TestAllWorkflowTypes:
    """Test all diagnostic workflow types."""

    @pytest.fixture
    def all_workflows(self):
        return list(DiagnosticWorkflowType)

    def test_workflow_count(self, all_workflows):
        assert len(all_workflows) == 21

    def test_phenotype_driven_workflow(self):
        wf = DiagnosticWorkflowType.PHENOTYPE_DRIVEN
        assert wf.value == "phenotype_driven"
        assert isinstance(wf, str)

    def test_variant_interpretation_workflow(self):
        wf = DiagnosticWorkflowType.VARIANT_INTERPRETATION
        assert wf.value == "variant_interpretation"

    def test_differential_diagnosis_workflow(self):
        wf = DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS
        assert wf.value == "differential_diagnosis"

    def test_gene_therapy_eligibility_workflow(self):
        wf = DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY
        assert wf.value == "gene_therapy_eligibility"

    def test_newborn_screening_workflow(self):
        wf = DiagnosticWorkflowType.NEWBORN_SCREENING
        assert wf.value == "newborn_screening"

    def test_metabolic_workup_workflow(self):
        wf = DiagnosticWorkflowType.METABOLIC_WORKUP
        assert wf.value == "metabolic_workup"

    def test_carrier_screening_workflow(self):
        wf = DiagnosticWorkflowType.CARRIER_SCREENING
        assert wf.value == "carrier_screening"

    def test_prenatal_diagnosis_workflow(self):
        wf = DiagnosticWorkflowType.PRENATAL_DIAGNOSIS
        assert wf.value == "prenatal_diagnosis"

    def test_natural_history_workflow(self):
        wf = DiagnosticWorkflowType.NATURAL_HISTORY
        assert wf.value == "natural_history"

    def test_therapy_selection_workflow(self):
        wf = DiagnosticWorkflowType.THERAPY_SELECTION
        assert wf.value == "therapy_selection"

    def test_clinical_trial_matching_workflow(self):
        wf = DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING
        assert wf.value == "clinical_trial_matching"

    def test_genetic_counseling_workflow(self):
        wf = DiagnosticWorkflowType.GENETIC_COUNSELING
        assert wf.value == "genetic_counseling"

    def test_general_workflow(self):
        wf = DiagnosticWorkflowType.GENERAL
        assert wf.value == "general"

    def test_workflow_roundtrip(self, all_workflows):
        """All workflow values should roundtrip correctly."""
        for wf in all_workflows:
            assert DiagnosticWorkflowType(wf.value) is wf

    def test_workflow_values_are_snake_case(self, all_workflows):
        """All workflow values should be lowercase snake_case."""
        for wf in all_workflows:
            assert wf.value == wf.value.lower()
            assert " " not in wf.value

    def test_workflow_coverage_for_rare_disease(self):
        """Verify key rare disease workflows are covered."""
        values = {wf.value for wf in DiagnosticWorkflowType}
        assert "phenotype_driven" in values
        assert "variant_interpretation" in values
        assert "gene_therapy_eligibility" in values
        assert "newborn_screening" in values
        assert "metabolic_workup" in values
        assert "carrier_screening" in values
        assert "prenatal_diagnosis" in values
