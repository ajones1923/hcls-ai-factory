"""Integration tests for the Biomarker Intelligence Agent.

These tests exercise the full analysis pipeline with real computation
modules (no mocking of core logic). Only external services (Milvus,
Claude API) are stubbed.

Run with: python3 -m pytest tests/test_integration.py -v
Skip with: python3 -m pytest tests/ -m "not integration"

Author: Adam Jones
Date: March 2026
"""

import json

import pytest

from src.biological_age import BiologicalAgeCalculator, validate_biomarker_ranges
from src.disease_trajectory import DiseaseTrajectoryAnalyzer
from src.export import export_csv, export_fhir_diagnostic_report, export_json, export_markdown, export_pdf
from src.genotype_adjustment import GenotypeAdjuster
from src.knowledge import BIOMARKER_PLAUSIBLE_RANGES, GENOTYPE_THRESHOLDS, KNOWLEDGE_VERSION
from src.models import (
    AnalysisResult,
    BiologicalAgeResult,
    DiseaseCategory,
    DiseaseTrajectoryResult,
    GenotypeAdjustmentResult,
    MetabolizerPhenotype,
    PatientProfile,
    PGxResult,
    RiskLevel,
)
from src.pharmacogenomics import PharmacogenomicMapper
from src.report_generator import ReportGenerator

# =====================================================================
# Realistic patient profiles for testing
# =====================================================================

HEALTHY_45M = PatientProfile(
    patient_id="INT-HEALTHY-001",
    age=45,
    sex="M",
    biomarkers={
        "albumin": 4.5, "creatinine": 0.9, "glucose": 92,
        "hs_crp": 0.8, "lymphocyte_pct": 32, "mcv": 88,
        "rdw": 12.5, "alkaline_phosphatase": 65, "wbc": 5.8,
        "hba1c": 5.2, "fasting_glucose": 88, "ldl": 110,
        "hdl": 55, "triglycerides": 120, "alt": 25, "ast": 22,
        "tsh": 2.1, "ferritin": 150, "hemoglobin": 15.2,
        "platelets": 250, "homocysteine": 8.5,
        "vitamin_d_25oh": 42, "vitamin_b12": 450, "folate_serum": 12,
    },
    genotypes={
        "TCF7L2_rs7903146": "CC",
        "APOE": "E3/E3",
        "PNPLA3_rs738409": "CC",
        "MTHFR_rs1801133": "CC",
    },
    star_alleles={
        "CYP2D6": "*1/*1",
        "CYP2C19": "*1/*1",
        "SLCO1B1": "*1a/*1a",
    },
)

PREDIABETIC_62F = PatientProfile(
    patient_id="INT-PREDIAB-002",
    age=62,
    sex="F",
    biomarkers={
        "albumin": 3.8, "creatinine": 1.1, "glucose": 118,
        "hs_crp": 4.2, "lymphocyte_pct": 24, "mcv": 94,
        "rdw": 14.8, "alkaline_phosphatase": 95, "wbc": 7.8,
        "hba1c": 6.2, "fasting_glucose": 115, "ldl": 155,
        "hdl": 42, "triglycerides": 195, "alt": 38, "ast": 32,
        "tsh": 3.8, "ferritin": 280, "hemoglobin": 12.8,
        "platelets": 210, "lpa": 85, "homocysteine": 14.5,
        "vitamin_d_25oh": 22, "vitamin_b12": 320, "folate_serum": 6,
    },
    genotypes={
        "TCF7L2_rs7903146": "CT",
        "APOE": "E3/E4",
        "PNPLA3_rs738409": "CG",
        "MTHFR_rs1801133": "CT",
        "HFE_rs1800562": "GG",
    },
    star_alleles={
        "CYP2D6": "*4/*4",
        "CYP2C19": "*1/*2",
        "SLCO1B1": "*1a/*5",
    },
)

HIGH_RISK_78M = PatientProfile(
    patient_id="INT-HIGHRISK-003",
    age=78,
    sex="M",
    biomarkers={
        "albumin": 3.2, "creatinine": 1.8, "glucose": 165,
        "hs_crp": 12.5, "lymphocyte_pct": 18, "mcv": 102,
        "rdw": 16.5, "alkaline_phosphatase": 140, "wbc": 10.2,
        "hba1c": 8.1, "fasting_glucose": 155, "ldl": 185,
        "hdl": 35, "triglycerides": 310, "alt": 55, "ast": 48,
        "tsh": 0.3, "ferritin": 520, "hemoglobin": 11.5,
        "platelets": 145, "lpa": 180, "homocysteine": 22,
    },
    genotypes={
        "TCF7L2_rs7903146": "TT",
        "APOE": "E4/E4",
        "PNPLA3_rs738409": "GG",
        "MTHFR_rs1801133": "TT",
        "HFE_rs1800562": "GA",
    },
    star_alleles={
        "CYP2D6": "*4/*4",
        "CYP2C19": "*2/*2",
        "CYP2C9": "*2/*3",
        "SLCO1B1": "*5/*5",
    },
)


# =====================================================================
# Helpers -- build AnalysisResult from raw module outputs
# =====================================================================

def _build_bio_age_result(
    chronological_age: int,
    ba_raw: dict,
) -> BiologicalAgeResult:
    """Convert raw BiologicalAgeCalculator output to BiologicalAgeResult model."""
    phenoage = ba_raw.get("phenoage", ba_raw)
    return BiologicalAgeResult(
        chronological_age=chronological_age,
        biological_age=ba_raw["biological_age"],
        age_acceleration=ba_raw["age_acceleration"],
        phenoage_score=ba_raw["biological_age"],
        grimage_score=(
            ba_raw["grimage"]["grimage_score"]
            if ba_raw.get("grimage") else None
        ),
        mortality_risk=phenoage.get("mortality_score", 0.0),
        aging_drivers=phenoage.get("top_aging_drivers", []),
        confidence_interval=phenoage.get("confidence_interval"),
        risk_confidence=phenoage.get("risk_confidence"),
    )


def _build_analysis_result(
    profile: PatientProfile,
    bio_age: BiologicalAgeCalculator,
    trajectory: DiseaseTrajectoryAnalyzer,
    pgx: PharmacogenomicMapper,
    adjuster: GenotypeAdjuster,
) -> AnalysisResult:
    """Run all analysis modules and build a full AnalysisResult."""
    ba_raw = bio_age.calculate(profile.age, profile.biomarkers)
    bio_age_result = _build_bio_age_result(profile.age, ba_raw)

    trajs_raw = trajectory.analyze_all(
        profile.biomarkers, profile.genotypes, profile.age, profile.sex,
    )

    pgx_raw = pgx.map_all(
        star_alleles=profile.star_alleles, genotypes=profile.genotypes,
    )

    sex_arg = "male" if profile.sex == "M" else "female"
    adjs_raw = adjuster.adjust_all(
        profile.biomarkers, profile.genotypes, sex=sex_arg,
    )

    # Convert trajectory dicts to model objects (simplified -- omit unmatched)
    _disease_str_map = {"type2_diabetes": "diabetes"}
    disease_results = []
    for td in trajs_raw:
        disease_str = td.get("disease", "")
        disease_str = _disease_str_map.get(disease_str, disease_str)
        if disease_str not in [e.value for e in DiseaseCategory]:
            continue
        risk_str = td.get("risk_level", "LOW").lower()
        if risk_str not in [e.value for e in RiskLevel]:
            risk_str = "normal"
        genetic_factors = [
            f"{grf.get('gene', '')} {grf.get('genotype', '')}"
            for grf in td.get("genetic_risk_factors", [])
        ]
        disease_results.append(DiseaseTrajectoryResult(
            disease=DiseaseCategory(disease_str),
            risk_level=RiskLevel(risk_str),
            current_markers=td.get("current_markers", {}),
            genetic_risk_factors=genetic_factors,
            years_to_onset_estimate=td.get("years_to_onset_estimate"),
            intervention_recommendations=td.get("recommendations", []),
        ))

    # Convert PGx results to model objects
    pgx_results = []
    for gr in pgx_raw.get("gene_results", []):
        phenotype_str = gr.get("phenotype", "normal_metabolizer")
        if phenotype_str is None:
            continue
        phenotype_val = phenotype_str.lower().replace(" ", "_").replace("-", "_")
        if phenotype_val not in [e.value for e in MetabolizerPhenotype]:
            phenotype_val = "normal"
        pgx_results.append(PGxResult(
            gene=gr.get("gene", ""),
            star_alleles=gr.get("star_alleles", "") or "",
            phenotype=MetabolizerPhenotype(phenotype_val),
            drugs_affected=gr.get("affected_drugs", []),
        ))

    # Convert adjustments to model objects
    adjustments = []
    for adj in adjs_raw.get("adjustments", []):
        std_range = adj.get("standard_range", {})
        adj_range = adj.get("adjusted_range", {})
        unit = adj.get("unit", "")
        adjustments.append(GenotypeAdjustmentResult(
            biomarker=adj.get("biomarker", ""),
            standard_range=f"{std_range.get('lower', '')}-{std_range.get('upper', '')} {unit}".strip(),
            adjusted_range=f"{adj_range.get('lower', '')}-{adj_range.get('upper', '')} {unit}".strip(),
            genotype=adj.get("genotype_value", ""),
            gene=adj.get("gene_display_name", ""),
            rationale=adj.get("rationale", ""),
        ))

    return AnalysisResult(
        patient_profile=profile,
        biological_age=bio_age_result,
        disease_trajectories=disease_results,
        pgx_results=pgx_results,
        genotype_adjustments=adjustments,
        critical_alerts=[],
    )


# =====================================================================
# Integration Tests -- Full Pipeline
# =====================================================================

@pytest.mark.integration
class TestFullPipelineIntegration:
    """Test complete analysis pipeline with real modules, no mocking."""

    def setup_method(self):
        self.bio_age = BiologicalAgeCalculator()
        self.trajectory = DiseaseTrajectoryAnalyzer()
        self.pgx = PharmacogenomicMapper()
        self.adjuster = GenotypeAdjuster()
        self.reporter = ReportGenerator()

    def test_healthy_patient_low_risk(self):
        """Healthy 45M should have normal biological age and low disease risk."""
        p = HEALTHY_45M

        # Biological age
        ba = self.bio_age.calculate(p.age, p.biomarkers)
        assert ba["biological_age"] < p.age + 5  # Not significantly accelerated
        assert ba["mortality_risk"] in ("NORMAL", "LOW")
        assert "confidence_interval" in ba["phenoage"]

        # Disease trajectories
        trajs = self.trajectory.analyze_all(p.biomarkers, p.genotypes, p.age, p.sex)
        high_risk = [t for t in trajs if t.get("risk_level") in ("CRITICAL", "HIGH")]
        assert len(high_risk) == 0, f"Healthy patient should have no high-risk trajectories: {high_risk}"

        # PGx
        pgx = self.pgx.map_all(star_alleles=p.star_alleles, genotypes=p.genotypes)
        assert pgx["drugs_to_avoid"] == []  # Normal metabolizer, no drugs to avoid
        assert "guideline_versions" in pgx
        assert "drug_interactions" in pgx

        # Genotype adjustments
        adjs = self.adjuster.adjust_all(p.biomarkers, p.genotypes, sex="male")
        assert isinstance(adjs, dict)
        assert "adjustments" in adjs

    def test_prediabetic_patient_moderate_risk(self):
        """Prediabetic 62F should show metabolic risk and PGx alerts."""
        p = PREDIABETIC_62F

        # Biological age should be accelerated
        ba = self.bio_age.calculate(p.age, p.biomarkers)
        assert ba["age_acceleration"] > 0  # Aging faster than chronological

        # Disease trajectories should flag diabetes risk
        trajs = self.trajectory.analyze_all(p.biomarkers, p.genotypes, p.age, p.sex)
        diabetes_traj = [t for t in trajs if "diabet" in t.get("disease", "").lower()]
        assert len(diabetes_traj) > 0, "Should detect diabetes trajectory"

        # CYP2D6 *4/*4 = poor metabolizer
        pgx = self.pgx.map_all(star_alleles=p.star_alleles, genotypes=p.genotypes)
        assert len(pgx["drugs_to_avoid"]) > 0  # PM should have drugs to avoid

        # Genotype adjustments for TCF7L2 CT
        adjs = self.adjuster.adjust_all(p.biomarkers, p.genotypes, sex="female")
        assert isinstance(adjs, dict)
        assert "adjustments" in adjs

    def test_high_risk_patient_critical_findings(self):
        """High-risk 78M should trigger multiple critical alerts."""
        p = HIGH_RISK_78M

        # Biological age significantly accelerated
        ba = self.bio_age.calculate(p.age, p.biomarkers)
        assert ba["age_acceleration"] > 2  # Significantly accelerated
        assert ba["mortality_risk"] in ("HIGH", "MODERATE")

        # Multiple disease trajectories should be elevated
        trajs = self.trajectory.analyze_all(p.biomarkers, p.genotypes, p.age, p.sex)
        elevated = [t for t in trajs if t.get("risk_level") in ("CRITICAL", "HIGH", "MODERATE")]
        assert len(elevated) >= 2, "High-risk patient should have multiple elevated trajectories"

        # Multiple PGx poor metabolizers
        pgx = self.pgx.map_all(star_alleles=p.star_alleles, genotypes=p.genotypes)
        assert len(pgx["drugs_to_avoid"]) >= 2  # Multiple PMs

    def test_pipeline_data_flows_to_report(self):
        """Verify data flows correctly from analysis through to report generation."""
        p = PREDIABETIC_62F

        analysis = _build_analysis_result(
            p, self.bio_age, self.trajectory, self.pgx, self.adjuster,
        )

        # Generate report
        report = self.reporter.generate(analysis, p)
        assert isinstance(report, str)
        assert len(report) > 1000
        assert p.patient_id in report
        assert "Biological Age" in report
        assert "Clinical Validation" in report
        assert "Knowledge base v" in report


@pytest.mark.integration
class TestExportPipelineIntegration:
    """Test full export pipeline with real data."""

    def setup_method(self):
        self.bio_age = BiologicalAgeCalculator()
        self.trajectory = DiseaseTrajectoryAnalyzer()
        self.pgx = PharmacogenomicMapper()
        self.adjuster = GenotypeAdjuster()
        self.reporter = ReportGenerator()
        self.analysis = _build_analysis_result(
            HEALTHY_45M, self.bio_age, self.trajectory, self.pgx, self.adjuster,
        )
        self.profile = HEALTHY_45M

    def test_markdown_export_complete(self):
        report = self.reporter.generate(self.analysis, self.profile)
        md = export_markdown(query="Integration test query", response_text=report)
        assert isinstance(md, str)
        assert "Biological Age" in md

    def test_json_export_roundtrip(self):
        json_str = export_json(analysis_result=self.analysis, query="Integration test")
        parsed = json.loads(json_str)
        assert parsed["report_type"] == "precision_biomarker_analysis"
        assert "analysis" in parsed
        assert parsed["analysis"]["patient_profile"]["patient_id"] == "INT-HEALTHY-001"

    def test_pdf_export_valid(self):
        report = self.reporter.generate(self.analysis, self.profile)
        pdf_bytes = export_pdf(report)
        assert isinstance(pdf_bytes, bytes)
        assert len(pdf_bytes) > 1000
        assert pdf_bytes[:5] == b"%PDF-"

    def test_csv_export_parseable(self):
        import csv
        import io
        csv_bytes = export_csv(self.analysis)
        reader = csv.reader(io.StringIO(csv_bytes.decode("utf-8")))
        rows = list(reader)
        assert len(rows) > 5  # Header + data rows
        # Should have section labels in first column (non-empty, non-blank)
        section_labels = set()
        for r in rows:
            if r and r[0].strip():
                section_labels.add(r[0].strip())
        # Expect at least Patient Info, Biological Age, Biomarker Values sections
        assert len(section_labels) >= 3, f"Expected section labels, got: {section_labels}"

    def test_fhir_bundle_valid(self):
        fhir_str = export_fhir_diagnostic_report(self.analysis, self.profile)
        fhir = json.loads(fhir_str)
        assert fhir["resourceType"] == "Bundle"
        assert fhir["type"] == "collection"
        assert len(fhir["entry"]) > 0
        # First entry should be DiagnosticReport
        dr = fhir["entry"][0]["resource"]
        assert dr["resourceType"] == "DiagnosticReport"
        assert dr["status"] == "final"


@pytest.mark.integration
class TestCrossModuleConsistency:
    """Verify data consistency across all analysis modules."""

    def test_biomarker_keys_consistent_across_modules(self):
        """All modules should accept the same biomarker key format."""
        bio_age = BiologicalAgeCalculator()
        trajectory = DiseaseTrajectoryAnalyzer()

        # Provide all 9 PhenoAge biomarkers so the result is reliable
        shared_biomarkers = {
            "albumin": 4.2, "creatinine": 0.9, "glucose": 95,
            "hs_crp": 1.5, "lymphocyte_pct": 30, "mcv": 89,
            "rdw": 13.0, "alkaline_phosphatase": 65, "wbc": 6.0,
            "hba1c": 5.4, "fasting_glucose": 90,
        }

        # Both should accept without error
        ba = bio_age.calculate(50, shared_biomarkers)
        assert ba["biological_age"] > 0

        trajs = trajectory.analyze_all(shared_biomarkers, {}, 50, "M")
        assert isinstance(trajs, list)

    def test_genotype_keys_consistent(self):
        """Genotype format should work across trajectory and adjustment modules."""
        trajectory = DiseaseTrajectoryAnalyzer()
        adjuster = GenotypeAdjuster()

        genotypes = {"TCF7L2_rs7903146": "CT", "PNPLA3_rs738409": "CG"}
        biomarkers = {"hba1c": 5.8, "fasting_glucose": 105, "alt": 35}

        trajs = trajectory.analyze_all(biomarkers, genotypes, 55, "M")
        adjs = adjuster.adjust_all(biomarkers, genotypes, sex="male")

        assert isinstance(trajs, list)
        assert isinstance(adjs, dict)
        assert "adjustments" in adjs

    def test_pgx_star_alleles_multiple_genes(self):
        """Multiple supported genes should be mappable simultaneously."""
        pgx = PharmacogenomicMapper()

        all_genes = {
            "CYP2D6": "*1/*1", "CYP2C19": "*1/*1", "CYP2C9": "*1/*1",
            "VKORC1": "*1/*1", "SLCO1B1": "*1a/*1a", "CYP3A5": "*1/*1",
            "TPMT": "*1/*1", "DPYD": "*1/*1",
        }
        genotypes = {
            "MTHFR_rs1801133": "CC",
        }

        result = pgx.map_all(star_alleles=all_genes, genotypes=genotypes)
        assert result["genes_analyzed"] >= 8  # At minimum the star allele genes

    def test_plausible_ranges_cover_phenoage_biomarkers(self):
        """All PhenoAge biomarkers should have plausible range validation."""
        from src.biological_age import PHENOAGE_COEFFICIENTS
        phenoage_markers = set(PHENOAGE_COEFFICIENTS.keys())
        plausible_markers = set(BIOMARKER_PLAUSIBLE_RANGES.keys())

        # ln_crp is derived from hs_crp, so check hs_crp instead
        check_markers = phenoage_markers - {"ln_crp"}
        check_markers.add("hs_crp")

        missing = check_markers - plausible_markers
        assert missing == set(), f"PhenoAge biomarkers missing plausible ranges: {missing}"

    def test_genotype_thresholds_used_by_trajectory(self):
        """GENOTYPE_THRESHOLDS should be importable and non-empty."""
        assert len(GENOTYPE_THRESHOLDS) >= 4
        assert "TCF7L2_hba1c" in GENOTYPE_THRESHOLDS
        assert "PNPLA3_alt_upper" in GENOTYPE_THRESHOLDS

    def test_knowledge_version_is_current(self):
        """Knowledge base version should be recent."""
        assert KNOWLEDGE_VERSION["version"] == "1.0.0"
        assert "2026" in KNOWLEDGE_VERSION["last_updated"]

    def test_ancestry_adjustments_available(self):
        """Ancestry adjustments should cover major populations."""
        from src.knowledge import ANCESTRY_ADJUSTMENTS
        assert "african" in ANCESTRY_ADJUSTMENTS
        assert "south_asian" in ANCESTRY_ADJUSTMENTS
        assert "east_asian" in ANCESTRY_ADJUSTMENTS
        assert "hispanic" in ANCESTRY_ADJUSTMENTS

    def test_age_sex_ranges_available(self):
        """Age/sex reference ranges should be available."""
        from src.knowledge import AGE_SEX_REFERENCE_RANGES
        assert "creatinine" in AGE_SEX_REFERENCE_RANGES
        assert "M" in AGE_SEX_REFERENCE_RANGES["creatinine"]
        assert "F" in AGE_SEX_REFERENCE_RANGES["creatinine"]


@pytest.mark.integration
class TestInputValidationIntegration:
    """Test that input validation works end-to-end."""

    def test_implausible_values_logged_but_analysis_completes(self):
        """Analysis should complete even with edge-case values, logging warnings."""
        bio_age = BiologicalAgeCalculator()
        # Very high glucose -- should warn but still compute
        result = bio_age.calculate(50, {
            "albumin": 4.2, "creatinine": 0.9, "glucose": 450,
            "hs_crp": 1.5, "lymphocyte_pct": 30, "mcv": 90,
            "rdw": 13, "alkaline_phosphatase": 60, "wbc": 6.0,
        })
        assert result["biological_age"] > 0
        assert result["phenoage"]["confidence_interval"] is not None

    def test_validate_biomarker_ranges_catches_errors(self):
        """Validation should catch clearly wrong values."""
        warnings = validate_biomarker_ranges({
            "albumin": 0.1,    # Too low (min 1.0)
            "glucose": 5000,   # Too high (max 600)
            "creatinine": 0.9, # Normal -- no warning
        })
        assert len(warnings) == 2
        flagged = {w.split("=")[0] for w in warnings}
        assert "albumin" in flagged
        assert "glucose" in flagged


@pytest.mark.integration
class TestAuditIntegration:
    """Test audit logging integrates correctly."""

    def test_audit_all_action_types(self):
        from src.audit import AuditAction, audit_log
        for action in AuditAction:
            event_id = audit_log(action, patient_id="AUDIT-TEST")
            assert len(event_id) == 16

    def test_audit_with_details(self):
        from src.audit import AuditAction, audit_log
        event_id = audit_log(
            AuditAction.PATIENT_ANALYSIS,
            patient_id="AUDIT-DETAIL",
            request_id="req-123",
            details={"modules": ["bio_age", "pgx"]},
            source_ip="127.0.0.1",
        )
        assert isinstance(event_id, str)
