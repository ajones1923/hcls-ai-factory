"""Tests for api/routes/ — PGx clinical endpoints, reports, and events.

Tests the Pydantic request/response models, endpoint logic functions,
and route module structure.  All external dependencies (Milvus, LLM,
embedding model) are mocked.

Author: Adam Jones
Date: March 2026
"""

import pytest
from pydantic import ValidationError

# Import route modules to verify they have valid router objects
from api.routes import pgx_clinical, reports, events

# Import Pydantic models from pgx_clinical
from api.routes.pgx_clinical import (
    DrugCheckRequest,
    DrugCheckAlert,
    DrugCheckResponse,
    MedicationReviewRequest,
    MedicationInteraction,
    MedicationReviewResponse,
    WarfarinDosingRequest,
    WarfarinDosingResponse,
    HLAScreenRequest,
    HLAScreenResult,
    HLAScreenResponse,
    PhenoconversionRequest,
    PhenoconversionResult,
    PhenoconversionResponse,
    GeneInfo,
    DrugInfo,
)

# Import Pydantic models from reports and events
from api.routes.reports import ReportMeta, _build_report_data, _render_markdown
from api.routes.events import (
    PipelineEvent,
    EventListResponse,
    emit_event,
    _event_store,
)

# Import endpoint functions for direct testing
from api.routes.pgx_clinical import (
    drug_check,
    medication_review,
    warfarin_dosing,
    hla_screen,
    phenoconversion_analysis,
    list_pharmacogenes,
    list_tracked_drugs,
)

# Try to import pytest-asyncio; skip async tests gracefully if unavailable
try:
    import pytest_asyncio  # noqa: F401
    HAS_ASYNCIO = True
except ImportError:
    HAS_ASYNCIO = False

skipif_no_asyncio = pytest.mark.skipif(
    not HAS_ASYNCIO,
    reason="pytest-asyncio not installed",
)


# ═══════════════════════════════════════════════════════════════════════
# Route module structure
# ═══════════════════════════════════════════════════════════════════════

class TestRouteModules:
    """Verify that each route module exports a valid APIRouter."""

    def test_pgx_clinical_has_router(self):
        assert hasattr(pgx_clinical, "router")
        assert pgx_clinical.router.prefix == "/v1/pgx"

    def test_reports_has_router(self):
        assert hasattr(reports, "router")
        assert reports.router.prefix == "/api"

    def test_events_has_router(self):
        assert hasattr(events, "router")
        assert events.router.prefix == "/api"

    def test_pgx_clinical_router_has_routes(self):
        route_paths = [r.path for r in pgx_clinical.router.routes]
        assert "/v1/pgx/drug-check" in route_paths
        assert "/v1/pgx/medication-review" in route_paths
        assert "/v1/pgx/dosing/warfarin" in route_paths
        assert "/v1/pgx/hla-screen" in route_paths
        assert "/v1/pgx/phenoconversion" in route_paths
        assert "/v1/pgx/genes" in route_paths
        assert "/v1/pgx/drugs" in route_paths

    def test_reports_router_has_routes(self):
        route_paths = [r.path for r in reports.router.routes]
        assert "/api/reports/{patient_id}" in route_paths
        assert "/api/reports/{patient_id}/{fmt}" in route_paths

    def test_events_router_has_routes(self):
        route_paths = [r.path for r in events.router.routes]
        assert "/api/events" in route_paths
        assert "/api/events/{event_id}" in route_paths


# ═══════════════════════════════════════════════════════════════════════
# Pydantic model validation — DrugCheck
# ═══════════════════════════════════════════════════════════════════════

class TestDrugCheckModels:
    """Tests for DrugCheck request/response Pydantic models."""

    def test_drug_check_request_valid(self):
        req = DrugCheckRequest(drug="codeine", gene="CYP2D6", phenotype="poor_metabolizer")
        assert req.drug == "codeine"
        assert req.gene == "CYP2D6"
        assert req.phenotype == "poor_metabolizer"

    def test_drug_check_request_drug_required(self):
        with pytest.raises(ValidationError):
            DrugCheckRequest(gene="CYP2D6")

    def test_drug_check_request_empty_drug_rejected(self):
        with pytest.raises(ValidationError):
            DrugCheckRequest(drug="")

    def test_drug_check_request_optional_fields(self):
        req = DrugCheckRequest(drug="warfarin")
        assert req.gene is None
        assert req.phenotype is None
        assert req.diplotype is None

    def test_drug_check_alert_creation(self):
        alert = DrugCheckAlert(
            alert_level="critical",
            gene="CYP2D6",
            drug="codeine",
            phenotype="poor_metabolizer",
            recommendation="Avoid codeine",
            guideline_body="CPIC",
            cpic_level="A",
            alternative_drugs=["morphine", "oxycodone"],
        )
        assert alert.alert_level == "critical"
        assert len(alert.alternative_drugs) == 2

    def test_drug_check_response_creation(self):
        resp = DrugCheckResponse(
            drug="codeine",
            gene="CYP2D6",
            phenotype="poor_metabolizer",
            alerts=[DrugCheckAlert(alert_level="critical")],
        )
        assert resp.drug == "codeine"
        assert len(resp.alerts) == 1

    def test_drug_check_response_serialization(self):
        resp = DrugCheckResponse(drug="codeine", alerts=[])
        d = resp.model_dump()
        assert d["drug"] == "codeine"
        assert d["alerts"] == []
        assert "processing_time_ms" in d


# ═══════════════════════════════════════════════════════════════════════
# Pydantic model validation — MedicationReview
# ═══════════════════════════════════════════════════════════════════════

class TestMedicationReviewModels:
    """Tests for MedicationReview request/response models."""

    def test_medication_review_request_valid(self):
        req = MedicationReviewRequest(
            medications=["codeine", "fluoxetine", "omeprazole"],
            gene="CYP2D6",
            phenotype="normal_metabolizer",
        )
        assert len(req.medications) == 3
        assert req.include_interactions is True

    def test_medication_review_request_medications_required(self):
        with pytest.raises(ValidationError):
            MedicationReviewRequest(medications=[])

    def test_medication_interaction_creation(self):
        interaction = MedicationInteraction(
            drug_a="fluoxetine",
            drug_b="codeine",
            enzyme="CYP2D6",
            interaction_type="strong inhibition",
            severity="high",
        )
        assert interaction.drug_a == "fluoxetine"

    def test_medication_review_response_creation(self):
        resp = MedicationReviewResponse(medications_reviewed=3)
        assert resp.medications_reviewed == 3
        assert resp.alerts == []
        assert resp.interactions == []


# ═══════════════════════════════════════════════════════════════════════
# Pydantic model validation — Warfarin Dosing
# ═══════════════════════════════════════════════════════════════════════

class TestWarfarinDosingModels:
    """Tests for warfarin dosing request/response models."""

    def test_warfarin_request_valid(self):
        req = WarfarinDosingRequest(
            age=65, height_cm=170, weight_kg=80,
            vkorc1_genotype="GA", cyp2c9_genotype="*1/*3",
            race="white", amiodarone=False, enzyme_inducer=False,
        )
        assert req.age == 65
        assert req.target_inr == 2.5  # default

    def test_warfarin_request_age_bounds(self):
        with pytest.raises(ValidationError):
            WarfarinDosingRequest(
                age=15, height_cm=170, weight_kg=80,
                vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            )

    def test_warfarin_request_height_bounds(self):
        with pytest.raises(ValidationError):
            WarfarinDosingRequest(
                age=50, height_cm=50, weight_kg=80,
                vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            )

    def test_warfarin_response_creation(self):
        resp = WarfarinDosingResponse(
            predicted_weekly_dose_mg=35.0,
            predicted_daily_dose_mg=5.0,
            dose_category="standard",
        )
        assert resp.algorithm == "IWPC"
        assert resp.predicted_weekly_dose_mg == 35.0


# ═══════════════════════════════════════════════════════════════════════
# Pydantic model validation — HLA Screen
# ═══════════════════════════════════════════════════════════════════════

class TestHLAScreenModels:
    """Tests for HLA screening request/response models."""

    def test_hla_screen_request_valid(self):
        req = HLAScreenRequest(drug="abacavir", hla_alleles=["HLA-B*57:01"])
        assert req.drug == "abacavir"

    def test_hla_screen_request_drug_required(self):
        with pytest.raises(ValidationError):
            HLAScreenRequest(hla_alleles=["HLA-B*57:01"])

    def test_hla_screen_result_creation(self):
        result = HLAScreenResult(
            hla_allele="HLA-B*57:01",
            drug="abacavir",
            reaction_type="Abacavir Hypersensitivity Syndrome",
            severity="potentially fatal",
            screening_mandatory=True,
        )
        assert result.screening_mandatory is True

    def test_hla_screen_response_creation(self):
        resp = HLAScreenResponse(drug="abacavir", associations_found=1)
        assert resp.drug == "abacavir"


# ═══════════════════════════════════════════════════════════════════════
# Pydantic model validation — Phenoconversion
# ═══════════════════════════════════════════════════════════════════════

class TestPhenoconversionModels:
    """Tests for phenoconversion analysis request/response models."""

    def test_phenoconversion_request_valid(self):
        req = PhenoconversionRequest(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["fluoxetine"],
        )
        assert req.enzyme == "CYP2D6"

    def test_phenoconversion_request_drugs_required(self):
        with pytest.raises(ValidationError):
            PhenoconversionRequest(
                enzyme="CYP2D6",
                baseline_phenotype="normal_metabolizer",
                concomitant_drugs=[],
            )

    def test_phenoconversion_result_creation(self):
        result = PhenoconversionResult(
            precipitant_drug="fluoxetine",
            inhibitor_strength="strong",
            converted_phenotype="poor_metabolizer",
            affected_substrates=["codeine", "tramadol"],
        )
        assert result.inhibitor_strength == "strong"
        assert len(result.affected_substrates) == 2

    def test_phenoconversion_response_creation(self):
        resp = PhenoconversionResponse(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            effective_phenotype="poor_metabolizer",
            phenoconversion_detected=True,
        )
        assert resp.phenoconversion_detected is True


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — drug_check
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestDrugCheckEndpoint:
    """Tests for the drug_check endpoint logic."""

    @pytest.mark.asyncio
    async def test_codeine_cyp2d6_poor_metabolizer_critical(self):
        """codeine + CYP2D6 + poor_metabolizer should return critical alert."""
        req = DrugCheckRequest(
            drug="codeine", gene="CYP2D6", phenotype="poor_metabolizer",
        )
        resp = await drug_check(req)
        assert resp.drug == "codeine"
        assert resp.gene == "CYP2D6"
        assert len(resp.alerts) >= 1

        # Find the CPIC guideline alert
        cpic_alerts = [a for a in resp.alerts if a.guideline_body == "CPIC"]
        assert len(cpic_alerts) >= 1
        assert cpic_alerts[0].alert_level == "critical"
        assert "poor_metabolizer" in cpic_alerts[0].phenotype

    @pytest.mark.asyncio
    async def test_codeine_cyp2d6_normal_metabolizer_warning(self):
        """codeine + CYP2D6 + normal_metabolizer should return warning (not critical)."""
        req = DrugCheckRequest(
            drug="codeine", gene="CYP2D6", phenotype="normal_metabolizer",
        )
        resp = await drug_check(req)
        cpic_alerts = [a for a in resp.alerts if a.guideline_body == "CPIC"]
        assert len(cpic_alerts) >= 1
        assert cpic_alerts[0].alert_level == "warning"

    @pytest.mark.asyncio
    async def test_codeine_cyp2d6_ultra_rapid(self):
        """codeine + CYP2D6 + ultra_rapid should mention dose reduction."""
        req = DrugCheckRequest(
            drug="codeine", gene="CYP2D6", phenotype="ultra_rapid_metabolizer",
        )
        resp = await drug_check(req)
        cpic_alerts = [a for a in resp.alerts if a.guideline_body == "CPIC"]
        assert len(cpic_alerts) >= 1
        assert "dose reduction" in cpic_alerts[0].recommendation.lower() or "alternative" in cpic_alerts[0].recommendation.lower()

    @pytest.mark.asyncio
    async def test_drug_without_gene(self):
        """Drug-only check should not crash."""
        req = DrugCheckRequest(drug="aspirin")
        resp = await drug_check(req)
        assert resp.drug == "aspirin"
        assert isinstance(resp.alerts, list)

    @pytest.mark.asyncio
    async def test_abacavir_hla_screening(self):
        """abacavir should trigger HLA-B*57:01 screening alert."""
        req = DrugCheckRequest(drug="abacavir")
        resp = await drug_check(req)
        hla_alerts = [a for a in resp.alerts if "HLA" in a.gene]
        assert len(hla_alerts) >= 1
        assert hla_alerts[0].alert_level == "critical"
        assert "HLA-B*57:01" in hla_alerts[0].gene

    @pytest.mark.asyncio
    async def test_carbamazepine_hla_screening(self):
        """carbamazepine should trigger HLA-B*15:02 screening alert."""
        req = DrugCheckRequest(drug="carbamazepine")
        resp = await drug_check(req)
        hla_alerts = [a for a in resp.alerts if "HLA" in a.gene]
        assert len(hla_alerts) >= 1

    @pytest.mark.asyncio
    async def test_unknown_gene_drug_pair(self):
        """Unknown gene-drug pair should return info-level alert."""
        req = DrugCheckRequest(drug="aspirin", gene="CYP2D6")
        resp = await drug_check(req)
        info_alerts = [a for a in resp.alerts if a.alert_level == "info"]
        assert len(info_alerts) >= 1

    @pytest.mark.asyncio
    async def test_response_has_processing_time(self):
        req = DrugCheckRequest(drug="codeine")
        resp = await drug_check(req)
        assert resp.processing_time_ms >= 0


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — medication_review
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestMedicationReviewEndpoint:
    """Tests for the medication_review endpoint logic."""

    @pytest.mark.asyncio
    async def test_basic_review(self):
        """Basic medication review should return count."""
        req = MedicationReviewRequest(
            medications=["codeine", "omeprazole"],
        )
        resp = await medication_review(req)
        assert resp.medications_reviewed == 2

    @pytest.mark.asyncio
    async def test_review_with_cyp_inhibitor_interaction(self):
        """fluoxetine + codeine should flag CYP2D6 inhibition interaction."""
        req = MedicationReviewRequest(
            medications=["codeine", "fluoxetine"],
            gene="CYP2D6",
            phenotype="normal_metabolizer",
            include_interactions=True,
        )
        resp = await medication_review(req)
        # fluoxetine is a strong CYP2D6 inhibitor, codeine is a CYP2D6 substrate
        assert len(resp.interactions) >= 1
        interaction = resp.interactions[0]
        assert "CYP2D6" in interaction.enzyme
        assert "inhibition" in interaction.interaction_type

    @pytest.mark.asyncio
    async def test_review_flags_phenoconversion_risk(self):
        """fluoxetine + codeine should flag phenoconversion risk."""
        req = MedicationReviewRequest(
            medications=["codeine", "fluoxetine"],
            include_interactions=True,
        )
        resp = await medication_review(req)
        assert len(resp.phenoconversion_risks) >= 1

    @pytest.mark.asyncio
    async def test_review_with_hla_drug(self):
        """abacavir in medication list should trigger HLA alert."""
        req = MedicationReviewRequest(medications=["abacavir", "metformin"])
        resp = await medication_review(req)
        hla_alerts = [a for a in resp.alerts if "HLA" in a.gene]
        assert len(hla_alerts) >= 1

    @pytest.mark.asyncio
    async def test_review_without_interactions(self):
        """Setting include_interactions=False should skip DDI analysis."""
        req = MedicationReviewRequest(
            medications=["codeine", "fluoxetine"],
            include_interactions=False,
        )
        resp = await medication_review(req)
        assert resp.interactions == []


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — warfarin_dosing
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestWarfarinDosingEndpoint:
    """Tests for the IWPC warfarin dosing algorithm endpoint."""

    @pytest.mark.asyncio
    async def test_standard_dose_reference_genotype(self):
        """GG + *1/*1 + standard demographics should give standard dose."""
        req = WarfarinDosingRequest(
            age=60, height_cm=170, weight_kg=80,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
        )
        resp = await warfarin_dosing(req)
        assert resp.predicted_weekly_dose_mg > 0
        assert resp.predicted_daily_dose_mg > 0
        assert resp.algorithm == "IWPC"
        assert resp.dose_category in ("low", "standard", "high")

    @pytest.mark.asyncio
    async def test_low_dose_aa_star3star3(self):
        """AA + *3/*3 should produce a low dose."""
        req = WarfarinDosingRequest(
            age=70, height_cm=160, weight_kg=60,
            vkorc1_genotype="AA", cyp2c9_genotype="*3/*3",
        )
        resp = await warfarin_dosing(req)
        assert resp.predicted_weekly_dose_mg < 21
        assert resp.dose_category == "low"

    @pytest.mark.asyncio
    async def test_amiodarone_warning(self):
        """Amiodarone use should produce a warning."""
        req = WarfarinDosingRequest(
            age=60, height_cm=170, weight_kg=80,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            amiodarone=True,
        )
        resp = await warfarin_dosing(req)
        assert any("Amiodarone" in w for w in resp.warnings)

    @pytest.mark.asyncio
    async def test_enzyme_inducer_increases_dose(self):
        """Enzyme inducer should increase predicted dose."""
        base_req = WarfarinDosingRequest(
            age=50, height_cm=175, weight_kg=85,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
        )
        inducer_req = WarfarinDosingRequest(
            age=50, height_cm=175, weight_kg=85,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            enzyme_inducer=True,
        )
        base_resp = await warfarin_dosing(base_req)
        inducer_resp = await warfarin_dosing(inducer_req)
        assert inducer_resp.predicted_weekly_dose_mg > base_resp.predicted_weekly_dose_mg

    @pytest.mark.asyncio
    async def test_unrecognized_vkorc1_warning(self):
        """Unrecognized VKORC1 genotype should produce a warning."""
        req = WarfarinDosingRequest(
            age=50, height_cm=170, weight_kg=80,
            vkorc1_genotype="XY", cyp2c9_genotype="*1/*1",
        )
        resp = await warfarin_dosing(req)
        assert any("VKORC1" in w for w in resp.warnings)

    @pytest.mark.asyncio
    async def test_race_asian_reduces_dose(self):
        """Asian race should reduce dose compared to white/reference."""
        white_req = WarfarinDosingRequest(
            age=50, height_cm=170, weight_kg=80,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            race="white",
        )
        asian_req = WarfarinDosingRequest(
            age=50, height_cm=170, weight_kg=80,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
            race="asian",
        )
        white_resp = await warfarin_dosing(white_req)
        asian_resp = await warfarin_dosing(asian_req)
        assert asian_resp.predicted_weekly_dose_mg < white_resp.predicted_weekly_dose_mg

    @pytest.mark.asyncio
    async def test_input_summary_populated(self):
        req = WarfarinDosingRequest(
            age=50, height_cm=170, weight_kg=80,
            vkorc1_genotype="GG", cyp2c9_genotype="*1/*1",
        )
        resp = await warfarin_dosing(req)
        assert resp.input_summary["age"] == 50
        assert resp.input_summary["vkorc1"] == "GG"

    @pytest.mark.asyncio
    async def test_dose_clamped_to_safety_bounds(self):
        """Extreme inputs should be clamped to 3-105 mg/week."""
        req = WarfarinDosingRequest(
            age=120, height_cm=100, weight_kg=30,
            vkorc1_genotype="AA", cyp2c9_genotype="*3/*3",
        )
        resp = await warfarin_dosing(req)
        assert resp.predicted_weekly_dose_mg >= 3.0


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — hla_screen
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestHLAScreenEndpoint:
    """Tests for the HLA hypersensitivity screening endpoint."""

    @pytest.mark.asyncio
    async def test_abacavir_finds_hla_b5701(self):
        """abacavir should return HLA-B*57:01 association."""
        req = HLAScreenRequest(drug="abacavir")
        resp = await hla_screen(req)
        assert resp.associations_found >= 1
        alleles = [r.hla_allele for r in resp.results]
        assert "HLA-B*57:01" in alleles

    @pytest.mark.asyncio
    async def test_abacavir_mandatory_screening(self):
        """abacavir HLA-B*57:01 screening should be mandatory."""
        req = HLAScreenRequest(drug="abacavir")
        resp = await hla_screen(req)
        b5701 = [r for r in resp.results if r.hla_allele == "HLA-B*57:01"]
        assert len(b5701) >= 1
        assert b5701[0].screening_mandatory is True

    @pytest.mark.asyncio
    async def test_patient_positive_contraindication(self):
        """Patient positive for HLA-B*57:01 should get CONTRAINDICATED message."""
        req = HLAScreenRequest(
            drug="abacavir", hla_alleles=["HLA-B*57:01"],
        )
        resp = await hla_screen(req)
        b5701 = [r for r in resp.results if r.hla_allele == "HLA-B*57:01"]
        assert len(b5701) >= 1
        assert "CONTRAINDICATED" in b5701[0].recommendation

    @pytest.mark.asyncio
    async def test_drug_without_hla_association(self):
        """A drug without HLA associations should return 0 results."""
        req = HLAScreenRequest(drug="metformin")
        resp = await hla_screen(req)
        assert resp.associations_found == 0
        assert resp.results == []

    @pytest.mark.asyncio
    async def test_carbamazepine_finds_hla_b1502(self):
        """carbamazepine should return HLA-B*15:02 association."""
        req = HLAScreenRequest(drug="carbamazepine")
        resp = await hla_screen(req)
        alleles = [r.hla_allele for r in resp.results]
        assert "HLA-B*15:02" in alleles


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — phenoconversion_analysis
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestPhenoconversionEndpoint:
    """Tests for the phenoconversion analysis endpoint."""

    @pytest.mark.asyncio
    async def test_strong_inhibitor_converts_nm_to_pm(self):
        """fluoxetine (strong CYP2D6 inhibitor) should convert NM to PM."""
        req = PhenoconversionRequest(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["fluoxetine"],
        )
        resp = await phenoconversion_analysis(req)
        assert resp.phenoconversion_detected is True
        assert resp.effective_phenotype == "poor_metabolizer"
        assert len(resp.results) >= 1
        assert resp.results[0].inhibitor_strength == "strong"

    @pytest.mark.asyncio
    async def test_moderate_inhibitor_converts_nm_to_im(self):
        """duloxetine (moderate CYP2D6 inhibitor) should convert NM to IM."""
        req = PhenoconversionRequest(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["duloxetine"],
        )
        resp = await phenoconversion_analysis(req)
        assert resp.phenoconversion_detected is True
        assert resp.effective_phenotype == "intermediate_metabolizer"

    @pytest.mark.asyncio
    async def test_no_inhibitor_no_conversion(self):
        """A non-inhibitor drug should not trigger phenoconversion."""
        req = PhenoconversionRequest(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["metformin"],
        )
        resp = await phenoconversion_analysis(req)
        assert resp.phenoconversion_detected is False
        assert resp.effective_phenotype == "normal_metabolizer"

    @pytest.mark.asyncio
    async def test_affected_substrates_listed(self):
        """Phenoconversion result should list affected substrate drugs."""
        req = PhenoconversionRequest(
            enzyme="CYP2D6",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["paroxetine"],
        )
        resp = await phenoconversion_analysis(req)
        assert resp.phenoconversion_detected is True
        assert len(resp.results[0].affected_substrates) >= 1

    @pytest.mark.asyncio
    async def test_enzyme_not_in_inhibitor_table(self):
        """An enzyme not in CYP_INHIBITORS should return no conversion."""
        req = PhenoconversionRequest(
            enzyme="CYP99Z9",
            baseline_phenotype="normal_metabolizer",
            concomitant_drugs=["fluoxetine"],
        )
        resp = await phenoconversion_analysis(req)
        assert resp.phenoconversion_detected is False


# ═══════════════════════════════════════════════════════════════════════
# Endpoint logic — list_pharmacogenes, list_tracked_drugs
# ═══════════════════════════════════════════════════════════════════════

@skipif_no_asyncio
class TestReferenceEndpoints:
    """Tests for the gene and drug reference list endpoints."""

    @pytest.mark.asyncio
    async def test_list_pharmacogenes_returns_genes(self):
        """list_pharmacogenes should return non-empty list of GeneInfo."""
        result = await list_pharmacogenes()
        assert len(result) >= 1
        assert isinstance(result[0], GeneInfo)
        gene_symbols = [g.gene for g in result]
        assert "CYP2D6" in gene_symbols
        assert "CYP2C19" in gene_symbols

    @pytest.mark.asyncio
    async def test_list_pharmacogenes_fields(self):
        """GeneInfo entries should have key fields populated."""
        result = await list_pharmacogenes()
        cyp2d6 = [g for g in result if g.gene == "CYP2D6"][0]
        assert cyp2d6.full_name != ""
        assert cyp2d6.chromosome != ""
        assert cyp2d6.substrates_count > 0
        assert len(cyp2d6.cpic_guidelines) > 0

    @pytest.mark.asyncio
    async def test_list_tracked_drugs_returns_drugs(self):
        """list_tracked_drugs should return non-empty list of DrugInfo."""
        result = await list_tracked_drugs()
        assert len(result) >= 1
        assert isinstance(result[0], DrugInfo)

    @pytest.mark.asyncio
    async def test_list_tracked_drugs_has_categories(self):
        """Each DrugInfo should have a non-empty category."""
        result = await list_tracked_drugs()
        for drug_info in result:
            assert drug_info.category != ""
            assert drug_info.drug != ""


# ═══════════════════════════════════════════════════════════════════════
# Reports module
# ═══════════════════════════════════════════════════════════════════════

class TestReports:
    """Tests for report generation helpers and models."""

    def test_report_meta_creation(self):
        meta = ReportMeta(
            patient_id="P001", format="json",
            generated_at="2026-03-12T00:00:00Z", sections=7,
        )
        assert meta.patient_id == "P001"
        assert meta.sections == 7

    def test_build_report_data(self):
        """_build_report_data should return structured report dict."""
        data = _build_report_data("P001")
        assert data["patient_id"] == "P001"
        assert "sections" in data
        assert len(data["sections"]) == 7
        # Each section should have heading and body
        for section in data["sections"]:
            assert "heading" in section
            assert "body" in section

    def test_render_markdown(self):
        """_render_markdown should produce valid markdown."""
        data = _build_report_data("P001")
        md = _render_markdown(data)
        assert md.startswith("# ")
        assert "## Diplotype" in md
        assert "## Drug Guideline" in md


# ═══════════════════════════════════════════════════════════════════════
# Events module
# ═══════════════════════════════════════════════════════════════════════

class TestEvents:
    """Tests for the event audit trail module."""

    @pytest.fixture(autouse=True)
    def clear_event_store(self):
        """Clear the in-memory event store before each test."""
        _event_store.clear()
        yield
        _event_store.clear()

    def test_emit_event_returns_id(self):
        eid = emit_event("query", source="test", summary="Test query")
        assert isinstance(eid, str)
        assert len(eid) == 12

    def test_emit_event_stores_record(self):
        emit_event("drug_check", source="pgx_clinical", summary="Checked codeine")
        assert len(_event_store) == 1
        record = _event_store[0]
        assert record["event_type"] == "drug_check"
        assert record["source"] == "pgx_clinical"
        assert "timestamp" in record

    def test_emit_event_with_metadata(self):
        emit_event("hla_screen", metadata={"drug": "abacavir", "hla": "HLA-B*57:01"})
        assert _event_store[0]["metadata"]["drug"] == "abacavir"

    def test_pipeline_event_model(self):
        event = PipelineEvent(
            event_id="abc123def456",
            event_type="query",
            timestamp="2026-03-12T00:00:00Z",
        )
        assert event.event_id == "abc123def456"

    def test_event_list_response_model(self):
        resp = EventListResponse(
            events=[], total=0, page=1, page_size=50,
        )
        assert resp.total == 0

    def test_multiple_events_stored(self):
        for i in range(5):
            emit_event("query", summary=f"Query {i}")
        assert len(_event_store) == 5
