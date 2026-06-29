"""Tests for all enums and Pydantic models in src/models.py.

Covers:
  - Enum member counts and values
  - SCWorkflowType members
  - SeverityLevel members
  - EvidenceLevel members
  - TMEClass, CellTypeConfidence, SpatialPlatform, etc.
  - Pydantic model validation (SCQuery, WorkflowResult, CellTypeAnnotation, etc.)

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    AssayType,
    BiomarkerCandidate,
    CARTTargetValidation,
    CellTypeAnnotation,
    CellTypeConfidence,
    ClusteringMethod,
    DrugResponsePrediction,
    EvidenceLevel,
    NormalizationMethod,
    ResistanceRisk,
    SCQuery,
    SCResponse,
    SCWorkflowType,
    SearchPlan,
    SeverityLevel,
    SpatialPlatform,
    TMEClass,
    TMEProfile,
    TrajectoryResult,
    TrajectoryType,
    WorkflowResult,
)


# ===================================================================
# ENUM TESTS
# ===================================================================


class TestSCWorkflowType:
    """Tests for SCWorkflowType enum."""

    def test_member_count(self):
        """SCWorkflowType must have exactly 11 members."""
        assert len(SCWorkflowType) == 11

    def test_cell_type_annotation_value(self):
        assert SCWorkflowType.CELL_TYPE_ANNOTATION.value == "cell_type_annotation"

    def test_tme_profiling_value(self):
        assert SCWorkflowType.TME_PROFILING.value == "tme_profiling"

    def test_drug_response_value(self):
        assert SCWorkflowType.DRUG_RESPONSE.value == "drug_response"

    def test_spatial_niche_value(self):
        assert SCWorkflowType.SPATIAL_NICHE.value == "spatial_niche"

    def test_trajectory_analysis_value(self):
        assert SCWorkflowType.TRAJECTORY_ANALYSIS.value == "trajectory_analysis"

    def test_ligand_receptor_value(self):
        assert SCWorkflowType.LIGAND_RECEPTOR.value == "ligand_receptor"

    def test_biomarker_discovery_value(self):
        assert SCWorkflowType.BIOMARKER_DISCOVERY.value == "biomarker_discovery"

    def test_cart_target_validation_value(self):
        assert SCWorkflowType.CART_TARGET_VALIDATION.value == "cart_target_validation"

    def test_treatment_monitoring_value(self):
        assert SCWorkflowType.TREATMENT_MONITORING.value == "treatment_monitoring"

    def test_general_value(self):
        assert SCWorkflowType.GENERAL.value == "general"


class TestSeverityLevel:
    """Tests for SeverityLevel enum."""

    def test_member_count(self):
        assert len(SeverityLevel) == 5

    def test_critical_value(self):
        assert SeverityLevel.CRITICAL.value == "critical"

    def test_informational_value(self):
        assert SeverityLevel.INFORMATIONAL.value == "informational"


class TestEvidenceLevel:
    """Tests for EvidenceLevel enum."""

    def test_member_count(self):
        assert len(EvidenceLevel) == 5

    def test_strong_value(self):
        assert EvidenceLevel.STRONG.value == "strong"

    def test_uncertain_value(self):
        assert EvidenceLevel.UNCERTAIN.value == "uncertain"


class TestTMEClass:
    """Tests for TMEClass enum."""

    def test_member_count(self):
        assert len(TMEClass) == 4

    def test_hot_inflamed_value(self):
        assert TMEClass.HOT_INFLAMED.value == "hot_inflamed"

    def test_cold_desert_value(self):
        assert TMEClass.COLD_DESERT.value == "cold_desert"


class TestSpatialPlatform:
    """Tests for SpatialPlatform enum."""

    def test_member_count(self):
        assert len(SpatialPlatform) == 4

    def test_visium_value(self):
        assert SpatialPlatform.VISIUM.value == "visium"

    def test_merfish_value(self):
        assert SpatialPlatform.MERFISH.value == "merfish"


class TestTrajectoryType:
    """Tests for TrajectoryType enum."""

    def test_member_count(self):
        assert len(TrajectoryType) == 6

    def test_differentiation_value(self):
        assert TrajectoryType.DIFFERENTIATION.value == "differentiation"

    def test_exhaustion_value(self):
        assert TrajectoryType.EXHAUSTION.value == "exhaustion"


class TestAssayType:
    """Tests for AssayType enum."""

    def test_member_count(self):
        assert len(AssayType) == 6


class TestResistanceRisk:
    """Tests for ResistanceRisk enum."""

    def test_member_count(self):
        assert len(ResistanceRisk) == 3


class TestClusteringMethod:
    """Tests for ClusteringMethod enum."""

    def test_member_count(self):
        assert len(ClusteringMethod) == 4

    def test_leiden_value(self):
        assert ClusteringMethod.LEIDEN.value == "leiden"


class TestNormalizationMethod:
    """Tests for NormalizationMethod enum."""

    def test_member_count(self):
        assert len(NormalizationMethod) == 4


# ===================================================================
# PYDANTIC MODEL TESTS
# ===================================================================


class TestSCQuery:
    """Tests for SCQuery model."""

    def test_minimal_query(self):
        q = SCQuery(query="What cell types are in PBMC?")
        assert q.query == "What cell types are in PBMC?"
        assert q.workflow_type is None

    def test_full_query(self):
        q = SCQuery(
            query="Annotate cell types in lung tumor",
            tissue_type="lung",
            disease_context="NSCLC",
            assay_type=AssayType.SCRNASEQ,
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            genes_of_interest=["CD3D", "CD8A"],
            top_k=20,
        )
        assert q.tissue_type == "lung"
        assert q.assay_type == AssayType.SCRNASEQ
        assert len(q.genes_of_interest) == 2


class TestCellTypeAnnotation:
    """Tests for CellTypeAnnotation model."""

    def test_basic_annotation(self):
        ann = CellTypeAnnotation(
            cluster_id="0",
            cell_type="CD8+ T cell",
            confidence=CellTypeConfidence.HIGH,
            confidence_score=0.95,
            marker_genes=["CD8A", "GZMB", "PRF1"],
            cell_count=500,
            fraction=0.15,
        )
        assert ann.cell_type == "CD8+ T cell"
        assert ann.confidence == CellTypeConfidence.HIGH
        assert ann.fraction == 0.15


class TestTMEProfile:
    """Tests for TMEProfile model."""

    def test_basic_profile(self):
        profile = TMEProfile(
            tme_class=TMEClass.HOT_INFLAMED,
            immune_score=0.8,
            stromal_score=0.3,
            exhaustion_signature=0.4,
            predicted_immunotherapy_response="responder",
        )
        assert profile.tme_class == TMEClass.HOT_INFLAMED
        assert profile.immune_score == 0.8


class TestDrugResponsePrediction:
    """Tests for DrugResponsePrediction model."""

    def test_basic_prediction(self):
        pred = DrugResponsePrediction(
            drug_name="Pembrolizumab",
            drug_class="checkpoint inhibitor",
            predicted_sensitivity=0.85,
            resistance_risk=ResistanceRisk.LOW,
        )
        assert pred.drug_name == "Pembrolizumab"
        assert pred.resistance_risk == ResistanceRisk.LOW


class TestWorkflowResult:
    """Tests for WorkflowResult model."""

    def test_cell_type_workflow(self):
        ann = CellTypeAnnotation(
            cluster_id="0",
            cell_type="T cell",
            confidence=CellTypeConfidence.HIGH,
            confidence_score=0.9,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[ann],
            severity=SeverityLevel.INFORMATIONAL,
        )
        assert result.workflow_type == SCWorkflowType.CELL_TYPE_ANNOTATION
        assert len(result.cell_annotations) == 1

    def test_model_dump(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
        )
        dumped = result.model_dump()
        assert dumped["workflow_type"] == "general"
        assert dumped["severity"] == "informational"


class TestSCResponse:
    """Tests for SCResponse model."""

    def test_basic_response(self):
        resp = SCResponse(
            query="What cell types?",
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            answer="The sample contains T cells and macrophages.",
            confidence=0.85,
        )
        assert resp.confidence == 0.85
        assert resp.workflow_type == SCWorkflowType.CELL_TYPE_ANNOTATION


class TestSearchPlan:
    """Tests for SearchPlan dataclass."""

    def test_default_plan(self):
        plan = SearchPlan()
        assert plan.workflow_type == SCWorkflowType.GENERAL
        assert plan.collections == []

    def test_configured_plan(self):
        plan = SearchPlan(
            workflow_type=SCWorkflowType.TME_PROFILING,
            collections=["sc_tme", "sc_cell_types"],
            weights={"sc_tme": 0.6, "sc_cell_types": 0.4},
        )
        assert len(plan.collections) == 2
        assert plan.weights["sc_tme"] == 0.6


class TestTrajectoryResult:
    """Tests for TrajectoryResult model."""

    def test_basic_trajectory(self):
        traj = TrajectoryResult(
            trajectory_id="traj_1",
            trajectory_type=TrajectoryType.EXHAUSTION,
            start_cell_type="Naive CD8+ T cell",
            end_cell_type="Exhausted CD8+ T cell",
            driver_genes=["TOX", "PDCD1", "HAVCR2"],
        )
        assert traj.trajectory_type == TrajectoryType.EXHAUSTION
        assert len(traj.driver_genes) == 3


class TestBiomarkerCandidate:
    """Tests for BiomarkerCandidate model."""

    def test_basic_biomarker(self):
        bm = BiomarkerCandidate(
            gene="CXCL13",
            biomarker_type="predictive",
            cell_type_specific="CD8+ T cell",
            fold_change=3.2,
            p_value_adjusted=0.001,
            specificity_score=0.92,
            surface_protein=False,
        )
        assert bm.gene == "CXCL13"
        assert bm.fold_change == 3.2


class TestCARTTargetValidation:
    """Tests for CARTTargetValidation model."""

    def test_basic_target(self):
        target = CARTTargetValidation(
            target_gene="CD19",
            tumor_expression_fraction=0.95,
            tumor_expression_level=5.2,
            on_target_off_tumor_risk="low",
            evidence_level=EvidenceLevel.STRONG,
        )
        assert target.target_gene == "CD19"
        assert target.tumor_expression_fraction == 0.95
