"""Tests for all enums and Pydantic models in src/models.py.

Covers:
  - Enum member counts and values
  - DiagnosticWorkflowType members

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import DiagnosticWorkflowType


# ===================================================================
# ENUM TESTS
# ===================================================================


class TestDiagnosticWorkflowType:
    """Tests for DiagnosticWorkflowType enum."""

    def test_member_count(self):
        """DiagnosticWorkflowType must have exactly 21 members."""
        assert len(DiagnosticWorkflowType) == 21

    def test_phenotype_driven_value(self):
        assert DiagnosticWorkflowType.PHENOTYPE_DRIVEN.value == "phenotype_driven"

    def test_variant_interpretation_value(self):
        assert DiagnosticWorkflowType.VARIANT_INTERPRETATION.value == "variant_interpretation"

    def test_differential_diagnosis_value(self):
        assert DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS.value == "differential_diagnosis"

    def test_gene_therapy_eligibility_value(self):
        assert DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY.value == "gene_therapy_eligibility"

    def test_newborn_screening_value(self):
        assert DiagnosticWorkflowType.NEWBORN_SCREENING.value == "newborn_screening"

    def test_metabolic_workup_value(self):
        assert DiagnosticWorkflowType.METABOLIC_WORKUP.value == "metabolic_workup"

    def test_carrier_screening_value(self):
        assert DiagnosticWorkflowType.CARRIER_SCREENING.value == "carrier_screening"

    def test_prenatal_diagnosis_value(self):
        assert DiagnosticWorkflowType.PRENATAL_DIAGNOSIS.value == "prenatal_diagnosis"

    def test_natural_history_value(self):
        assert DiagnosticWorkflowType.NATURAL_HISTORY.value == "natural_history"

    def test_therapy_selection_value(self):
        assert DiagnosticWorkflowType.THERAPY_SELECTION.value == "therapy_selection"

    def test_clinical_trial_matching_value(self):
        assert DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING.value == "clinical_trial_matching"

    def test_genetic_counseling_value(self):
        assert DiagnosticWorkflowType.GENETIC_COUNSELING.value == "genetic_counseling"

    def test_general_value(self):
        assert DiagnosticWorkflowType.GENERAL.value == "general"

    def test_all_members_are_strings(self):
        for member in DiagnosticWorkflowType:
            assert isinstance(member.value, str)

    def test_from_value_roundtrip(self):
        for member in DiagnosticWorkflowType:
            assert DiagnosticWorkflowType(member.value) == member

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            DiagnosticWorkflowType("nonexistent_workflow")

    def test_all_values_unique(self):
        values = [m.value for m in DiagnosticWorkflowType]
        assert len(values) == len(set(values))

    def test_all_names_unique(self):
        names = [m.name for m in DiagnosticWorkflowType]
        assert len(names) == len(set(names))

    def test_str_enum_behavior(self):
        """DiagnosticWorkflowType members should behave as strings."""
        assert DiagnosticWorkflowType.GENERAL == "general"
