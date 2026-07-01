"""Tests for single-cell clinical workflows.

Tests all 11 workflow types and WorkflowResult construction.

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    CARTTargetValidation,
    CellTypeAnnotation,
    CellTypeConfidence,
    DrugResponsePrediction,
    EvidenceLevel,
    LigandReceptorInteraction,
    ResistanceRisk,
    SCWorkflowType,
    SeverityLevel,
    SpatialNiche,
    SpatialPlatform,
    SubclonalResult,
    TMEClass,
    TMEProfile,
    TrajectoryResult,
    TrajectoryType,
    TreatmentResponse,
    WorkflowResult,
)


class TestWorkflowTypes:
    """Tests for workflow type completeness."""

    def test_all_11_workflows_exist(self):
        expected = [
            "cell_type_annotation",
            "tme_profiling",
            "drug_response",
            "subclonal_architecture",
            "spatial_niche",
            "trajectory_analysis",
            "ligand_receptor",
            "biomarker_discovery",
            "cart_target_validation",
            "treatment_monitoring",
            "general",
        ]
        for wf_value in expected:
            assert SCWorkflowType(wf_value) is not None

    def test_general_workflow_exists(self):
        assert SCWorkflowType.GENERAL.value == "general"


class TestCellTypeAnnotationWorkflow:
    """Tests for cell type annotation workflow results."""

    def test_cell_type_result(self):
        ann = CellTypeAnnotation(
            cluster_id="0",
            cell_type="CD8+ exhausted T cell",
            confidence=CellTypeConfidence.HIGH,
            confidence_score=0.92,
            marker_genes=["CD8A", "PDCD1", "HAVCR2", "LAG3"],
            cell_count=1200,
            fraction=0.08,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[ann],
            severity=SeverityLevel.INFORMATIONAL,
        )
        assert result.workflow_type == SCWorkflowType.CELL_TYPE_ANNOTATION
        assert len(result.cell_annotations) == 1
        assert result.cell_annotations[0].cell_type == "CD8+ exhausted T cell"


class TestTMEProfilingWorkflow:
    """Tests for TME profiling workflow results."""

    def test_tme_result(self):
        profile = TMEProfile(
            tme_class=TMEClass.HOT_INFLAMED,
            immune_score=0.82,
            stromal_score=0.25,
            exhaustion_signature=0.45,
            predicted_immunotherapy_response="responder",
            evidence_level=EvidenceLevel.MODERATE,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.TME_PROFILING,
            tme_profile=profile,
            severity=SeverityLevel.MODERATE,
        )
        assert result.tme_profile.tme_class == TMEClass.HOT_INFLAMED
        assert result.tme_profile.immune_score == 0.82


class TestDrugResponseWorkflow:
    """Tests for drug response prediction workflow results."""

    def test_drug_response_result(self):
        pred = DrugResponsePrediction(
            drug_name="Pembrolizumab",
            drug_class="anti-PD-1",
            predicted_sensitivity=0.78,
            resistance_risk=ResistanceRisk.LOW,
            resistance_mechanisms=["TIM-3 upregulation"],
            synergy_candidates=["anti-LAG-3"],
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.DRUG_RESPONSE,
            drug_predictions=[pred],
        )
        assert len(result.drug_predictions) == 1
        assert result.drug_predictions[0].drug_name == "Pembrolizumab"


class TestSubclonalArchitectureWorkflow:
    """Tests for subclonal architecture workflow results."""

    def test_subclonal_result(self):
        clone = SubclonalResult(
            clone_id="clone_A",
            clone_fraction=0.6,
            cell_count=3000,
            driver_mutations=["KRAS G12D", "TP53 R175H"],
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.SUBCLONAL_ARCHITECTURE,
            subclones=[clone],
        )
        assert len(result.subclones) == 1
        assert result.subclones[0].clone_fraction == 0.6


class TestSpatialNicheWorkflow:
    """Tests for spatial niche workflow results."""

    def test_spatial_result(self):
        niche = SpatialNiche(
            niche_id="niche_1",
            niche_label="Tumor-immune interface",
            dominant_cell_types=["CD8+ T cell", "Macrophage", "Tumor cell"],
            spatial_platform=SpatialPlatform.VISIUM,
            area_fraction=0.15,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.SPATIAL_NICHE,
            spatial_niches=[niche],
        )
        assert len(result.spatial_niches) == 1
        assert result.spatial_niches[0].niche_label == "Tumor-immune interface"


class TestTrajectoryWorkflow:
    """Tests for trajectory analysis workflow results."""

    def test_trajectory_result(self):
        traj = TrajectoryResult(
            trajectory_id="traj_1",
            trajectory_type=TrajectoryType.EXHAUSTION,
            start_cell_type="Naive CD8+ T cell",
            end_cell_type="Exhausted CD8+ T cell",
            driver_genes=["TOX", "PDCD1", "HAVCR2"],
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.TRAJECTORY_ANALYSIS,
            trajectories=[traj],
        )
        assert len(result.trajectories) == 1


class TestLigandReceptorWorkflow:
    """Tests for ligand-receptor interaction workflow results."""

    def test_lr_result(self):
        interaction = LigandReceptorInteraction(
            ligand_gene="CXCL12",
            receptor_gene="CXCR4",
            source_cell_type="CAF",
            target_cell_type="T cell",
            interaction_score=0.9,
            p_value=0.001,
            pathway="CXCL12-CXCR4",
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.LIGAND_RECEPTOR,
            interactions=[interaction],
        )
        assert len(result.interactions) == 1


class TestCARTTargetWorkflow:
    """Tests for CAR-T target validation workflow results."""

    def test_cart_result(self):
        target = CARTTargetValidation(
            target_gene="CD19",
            tumor_expression_fraction=0.95,
            tumor_expression_level=5.2,
            on_target_off_tumor_risk="low",
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CART_TARGET_VALIDATION,
            cart_targets=[target],
        )
        assert len(result.cart_targets) == 1


class TestTreatmentMonitoringWorkflow:
    """Tests for treatment monitoring workflow results."""

    def test_treatment_result(self):
        response = TreatmentResponse(
            timepoint="week_4",
            treatment="Pembrolizumab",
            response_category="partial",
            compositional_shifts={"CD8+ T cell": 0.12, "Treg": -0.05},
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.TREATMENT_MONITORING,
            treatment_responses=[response],
        )
        assert len(result.treatment_responses) == 1
        assert result.treatment_responses[0].response_category == "partial"
