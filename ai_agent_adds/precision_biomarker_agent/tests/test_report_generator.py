"""Tests for Biomarker Intelligence Agent report generation.

Validates that reports contain all 12 sections, critical alerts are
highlighted, PGx warnings are formatted, and markdown is valid.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    AnalysisResult,
    BiologicalAgeResult,
    DiseaseCategory,
    DiseaseTrajectoryResult,
    GenotypeAdjustmentResult,
    MetabolizerPhenotype,
    PGxResult,
    RiskLevel,
)
from src.report_generator import ReportGenerator

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def full_analysis(sample_patient_profile):
    """Return a complete AnalysisResult for report testing."""
    return AnalysisResult(
        patient_profile=sample_patient_profile,
        biological_age=BiologicalAgeResult(
            chronological_age=45,
            biological_age=47.2,
            age_acceleration=2.2,
            phenoage_score=0.003,
            mortality_risk=0.015,
            aging_drivers=[
                {"biomarker": "rdw", "value": 13.5, "contribution": 4.46, "direction": "aging"},
                {"biomarker": "glucose", "value": 105.0, "contribution": 2.05, "direction": "aging"},
                {"biomarker": "albumin", "value": 4.2, "contribution": -0.14, "direction": "protective"},
            ],
        ),
        disease_trajectories=[
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.DIABETES,
                risk_level=RiskLevel.HIGH,
                current_markers={"hba1c": 5.9, "fasting_glucose": 105},
                genetic_risk_factors=["TCF7L2 rs7903146 CT"],
                intervention_recommendations=["Dietary optimization", "Recheck in 3 months"],
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.CARDIOVASCULAR,
                risk_level=RiskLevel.MODERATE,
                current_markers={"ldl_c": 145, "lpa": 65},
                intervention_recommendations=["Optimize modifiable risk factors"],
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.LIVER,
                risk_level=RiskLevel.LOW,
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.THYROID,
                risk_level=RiskLevel.LOW,
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.IRON,
                risk_level=RiskLevel.LOW,
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.NUTRITIONAL,
                risk_level=RiskLevel.MODERATE,
                intervention_recommendations=["Vitamin D supplementation"],
            ),
        ],
        pgx_results=[
            PGxResult(
                gene="CYP2D6",
                star_alleles="*1/*4",
                phenotype=MetabolizerPhenotype.INTERMEDIATE,
                drugs_affected=[
                    {"drug": "Codeine", "recommendation": "Use with caution", "cpic_level": "1A"},
                ],
            ),
            PGxResult(
                gene="CYP2C19",
                star_alleles="*1/*2",
                phenotype=MetabolizerPhenotype.INTERMEDIATE,
                drugs_affected=[
                    {"drug": "Clopidogrel", "recommendation": "Consider alternative", "cpic_level": "1A"},
                ],
            ),
        ],
        genotype_adjustments=[
            GenotypeAdjustmentResult(
                biomarker="ALT",
                standard_range="0-56 U/L",
                adjusted_range="0-45 U/L",
                genotype="CG",
                gene="PNPLA3",
                rationale="PNPLA3 CG carrier has 2x NAFLD risk; tighter ALT threshold warranted.",
            ),
        ],
        critical_alerts=[
            "HbA1c 5.9% in pre-diabetic range with TCF7L2 CT carrier status",
            "Lp(a) 65 nmol/L elevated - genetically determined cardiovascular risk",
        ],
    )


@pytest.fixture
def generator():
    """Return a ReportGenerator instance."""
    return ReportGenerator()


# =====================================================================
# REPORT SECTIONS
# =====================================================================


class TestReportSections:
    """Tests for report section completeness."""

    def test_all_12_sections_present(self, generator, full_analysis):
        """Report should contain all 12 sections."""
        report = generator.generate(full_analysis)
        assert "# Biomarker Intelligence Report" in report
        assert "## 1. Biological Age Assessment" in report
        assert "## 2. Executive Findings" in report
        assert "## 3. Biomarker-Gene Correlation Map" in report
        assert "## 4. Disease Trajectory Analysis" in report
        assert "## 5. Pharmacogenomic Profile" in report
        assert "## 6. Nutritional Genomics Analysis" in report
        assert "## 7. Interconnected Pathways" in report
        assert "## 8. Prioritized Action Plan" in report
        assert "## 9. Monitoring Schedule" in report
        assert "## 10. Supplement Protocol Summary" in report
        assert "## 11. Clinical Summary for Healthcare Provider" in report
        assert "## 12. References" in report

    def test_patient_info_in_header(self, generator, full_analysis):
        """Report header should include patient ID, age, and sex."""
        report = generator.generate(full_analysis)
        assert "HG002" in report
        assert "45" in report
        assert "M" in report

    def test_biological_age_values_present(self, generator, full_analysis):
        """Report should contain biological age values."""
        report = generator.generate(full_analysis)
        assert "47.2" in report
        assert "+2.2" in report

    def test_disease_trajectories_all_listed(self, generator, full_analysis):
        """All 6 disease categories should appear in the report."""
        report = generator.generate(full_analysis)
        assert "Diabetes" in report
        assert "Cardiovascular" in report
        assert "Liver" in report
        assert "Thyroid" in report
        assert "Iron" in report
        assert "Nutritional" in report

    def test_pgx_genes_listed(self, generator, full_analysis):
        """PGx gene results should appear in the report."""
        report = generator.generate(full_analysis)
        assert "CYP2D6" in report
        assert "CYP2C19" in report
        assert "*1/*4" in report

    def test_genotype_adjustments_listed(self, generator, full_analysis):
        """Genotype adjustments should appear in the report."""
        report = generator.generate(full_analysis)
        assert "PNPLA3" in report
        assert "0-45 U/L" in report

    def test_aging_drivers_table(self, generator, full_analysis):
        """Report should contain an aging drivers table."""
        report = generator.generate(full_analysis)
        assert "Top Aging Drivers" in report
        assert "rdw" in report
        assert "glucose" in report

    def test_monitoring_schedule_present(self, generator, full_analysis):
        """Report should contain monitoring schedule with disease-specific panels."""
        report = generator.generate(full_analysis)
        assert "PhenoAge Panel" in report
        assert "Every 3 months" in report or "Every 6 months" in report

    def test_references_contain_citations(self, generator, full_analysis):
        """Report should contain key academic references."""
        report = generator.generate(full_analysis)
        assert "Levine ME" in report
        assert "CPIC" in report


# =====================================================================
# CRITICAL ALERTS
# =====================================================================


class TestCriticalAlertFormatting:
    """Tests for critical alert highlighting in reports."""

    def test_critical_alerts_highlighted(self, generator, full_analysis):
        """Critical alerts should appear in executive findings."""
        report = generator.generate(full_analysis)
        assert "CRITICAL" in report

    def test_alert_contains_hba1c_warning(self, generator, full_analysis):
        """HbA1c pre-diabetic alert should be in report."""
        report = generator.generate(full_analysis)
        assert "HbA1c" in report
        assert "pre-diabetic" in report

    def test_alert_contains_lpa_warning(self, generator, full_analysis):
        """Lp(a) elevation alert should be in report."""
        report = generator.generate(full_analysis)
        assert "Lp(a)" in report

    def test_no_alerts_when_none(self, generator, sample_patient_profile):
        """Report without critical alerts should still generate."""
        analysis = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=45,
                biological_age=44.0,
                age_acceleration=-1.0,
                phenoage_score=0.002,
            ),
        )
        report = generator.generate(analysis)
        assert "No critical or high-priority findings" in report


# =====================================================================
# PGx WARNINGS
# =====================================================================


class TestPGxWarningFormatting:
    """Tests for PGx warning formatting in reports."""

    def test_pgx_phenotype_displayed(self, generator, full_analysis):
        """PGx phenotype should be displayed in report."""
        report = generator.generate(full_analysis)
        assert "Intermediate" in report

    def test_pgx_drug_recommendations_for_poor_metabolizer(self, generator, sample_patient_profile):
        """Drug-specific recommendations should appear for poor metabolizers."""
        analysis = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=45,
                biological_age=45.0,
                age_acceleration=0.0,
                phenoage_score=0.002,
            ),
            pgx_results=[
                PGxResult(
                    gene="CYP2D6",
                    star_alleles="*4/*4",
                    phenotype=MetabolizerPhenotype.POOR,
                    drugs_affected=[
                        {"drug": "Codeine", "recommendation": "Avoid",
                         "cpic_level": "1A"},
                    ],
                ),
            ],
        )
        report = generator.generate(analysis)
        assert "Codeine" in report
        assert "Avoid" in report

    def test_poor_metabolizer_triggers_alert(self, generator, sample_patient_profile):
        """Poor metabolizers should get HIGH alert level."""
        analysis = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=45,
                biological_age=45.0,
                age_acceleration=0.0,
                phenoage_score=0.002,
            ),
            pgx_results=[
                PGxResult(
                    gene="CYP2D6",
                    star_alleles="*4/*4",
                    phenotype=MetabolizerPhenotype.POOR,
                    drugs_affected=[
                        {"drug": "Codeine", "recommendation": "Avoid",
                         "cpic_level": "1A"},
                    ],
                ),
            ],
        )
        report = generator.generate(analysis)
        assert "HIGH" in report
        assert "CYP2D6" in report


# =====================================================================
# MARKDOWN VALIDITY
# =====================================================================


class TestMarkdownValidity:
    """Tests for markdown format validity."""

    def test_report_is_string(self, generator, full_analysis):
        """Report should be a non-empty string."""
        report = generator.generate(full_analysis)
        assert isinstance(report, str)
        assert len(report) > 500

    def test_headers_have_proper_format(self, generator, full_analysis):
        """All headers should use proper markdown # syntax."""
        report = generator.generate(full_analysis)
        lines = report.split("\n")
        headers = [line for line in lines if line.startswith("#")]
        for h in headers:
            assert h.startswith("# ") or h.startswith("## ") or h.startswith("### ")

    def test_tables_have_header_and_separator(self, generator, full_analysis):
        """Markdown tables should have header rows and separator rows."""
        report = generator.generate(full_analysis)
        lines = report.split("\n")
        table_headers = [line for line in lines if line.startswith("|") and "---" not in line]
        table_seps = [line for line in lines if line.startswith("|") and "---" in line]
        assert len(table_headers) > 0
        assert len(table_seps) > 0

    def test_no_unclosed_bold_markers(self, generator, full_analysis):
        """Bold markers (**) should appear in pairs."""
        report = generator.generate(full_analysis)
        bold_count = report.count("**")
        assert bold_count % 2 == 0

    def test_footer_present(self, generator, full_analysis):
        """Report should end with a footer disclaimer."""
        report = generator.generate(full_analysis)
        assert "HCLS AI Factory" in report
        assert "healthcare provider" in report.lower()
