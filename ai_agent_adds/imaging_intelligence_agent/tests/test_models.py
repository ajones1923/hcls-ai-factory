"""Tests for Imaging Intelligence Agent Pydantic data models.

Validates all 10 collection models, 17 enums, NIM result models,
search result models, embedding text generation, and agent I/O models.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.models import (
    # Enums
    AITaskType,
    BiRADS,
    BodyRegion,
    DeviceRegulatory,
    EvidenceLevel,
    FindingCategory,
    FindingSeverity,
    ImagingModality,
    LungRADS,
    ModelArchitecture,
    NIMServiceStatus,
    SourceType,
    TrialPhase,
    TrialStatus,
    WorkflowStatus,
    # Collection models
    AnatomyRecord,
    BenchmarkRecord,
    DatasetRecord,
    GuidelineRecord,
    ImagingDevice,
    ImagingFinding,
    ImagingLiterature,
    ImagingProtocol,
    ImagingTrial,
    ReportTemplate,
    # NIM result models
    SegmentationResult,
    SyntheticCTResult,
    VLMResponse,
    WorkflowResult,
    # Search result models
    AgentQuery,
    AgentResponse,
    ComparativeResult,
    CrossCollectionResult,
    SearchHit,
)


# ===================================================================
# ENUM VALUES
# ===================================================================


class TestImagingModality:
    """Tests for ImagingModality enum."""

    def test_has_ct(self):
        assert ImagingModality.CT.value == "ct"

    def test_has_mri(self):
        assert ImagingModality.MRI.value == "mri"

    def test_has_xray(self):
        assert ImagingModality.XRAY.value == "xray"

    def test_has_cxr(self):
        assert ImagingModality.CXR.value == "cxr"

    def test_has_ultrasound(self):
        assert ImagingModality.ULTRASOUND.value == "ultrasound"

    def test_has_pet(self):
        assert ImagingModality.PET.value == "pet"

    def test_has_pet_ct(self):
        assert ImagingModality.PET_CT.value == "pet_ct"

    def test_has_mammography(self):
        assert ImagingModality.MAMMOGRAPHY.value == "mammography"

    def test_has_fluoroscopy(self):
        assert ImagingModality.FLUOROSCOPY.value == "fluoroscopy"

    def test_count(self):
        assert len(ImagingModality) == 9


class TestBodyRegion:
    """Tests for BodyRegion enum."""

    def test_has_head(self):
        assert BodyRegion.HEAD.value == "head"

    def test_has_brain(self):
        assert BodyRegion.BRAIN.value == "brain"

    def test_has_chest(self):
        assert BodyRegion.CHEST.value == "chest"

    def test_has_abdomen(self):
        assert BodyRegion.ABDOMEN.value == "abdomen"

    def test_has_cardiac(self):
        assert BodyRegion.CARDIAC.value == "cardiac"

    def test_has_breast(self):
        assert BodyRegion.BREAST.value == "breast"

    def test_count(self):
        assert len(BodyRegion) == 11


class TestFindingSeverity:
    """Tests for FindingSeverity enum."""

    def test_has_critical(self):
        assert FindingSeverity.CRITICAL.value == "critical"

    def test_has_urgent(self):
        assert FindingSeverity.URGENT.value == "urgent"

    def test_has_significant(self):
        assert FindingSeverity.SIGNIFICANT.value == "significant"

    def test_has_routine(self):
        assert FindingSeverity.ROUTINE.value == "routine"

    def test_has_normal(self):
        assert FindingSeverity.NORMAL.value == "normal"

    def test_count(self):
        assert len(FindingSeverity) == 5


class TestFindingCategory:
    """Tests for FindingCategory enum."""

    def test_has_hemorrhage(self):
        assert FindingCategory.HEMORRHAGE.value == "hemorrhage"

    def test_has_nodule(self):
        assert FindingCategory.NODULE.value == "nodule"

    def test_has_pneumothorax(self):
        assert FindingCategory.PNEUMOTHORAX.value == "pneumothorax"

    def test_has_fracture(self):
        assert FindingCategory.FRACTURE.value == "fracture"

    def test_has_normal(self):
        assert FindingCategory.NORMAL.value == "normal"

    def test_count(self):
        assert len(FindingCategory) == 14


class TestOtherEnums:
    """Tests for remaining enums: SourceType, TrialPhase, etc."""

    def test_source_type_values(self):
        assert SourceType.PUBMED.value == "pubmed"
        assert SourceType.PMC.value == "pmc"
        assert SourceType.GUIDELINE.value == "guideline"
        assert len(SourceType) == 5

    def test_trial_phase_values(self):
        assert TrialPhase.EARLY_1.value == "Early Phase 1"
        assert TrialPhase.PHASE_3.value == "Phase 3"
        assert TrialPhase.NA.value == "N/A"
        assert len(TrialPhase) == 8

    def test_trial_status_values(self):
        assert TrialStatus.RECRUITING.value == "Recruiting"
        assert TrialStatus.COMPLETED.value == "Completed"
        assert TrialStatus.TERMINATED.value == "Terminated"
        assert len(TrialStatus) == 8

    def test_device_regulatory_values(self):
        assert DeviceRegulatory.CLEARED_510K.value == "510k_cleared"
        assert DeviceRegulatory.DE_NOVO.value == "de_novo"
        assert len(DeviceRegulatory) == 6

    def test_ai_task_type_values(self):
        assert AITaskType.DETECTION.value == "detection"
        assert AITaskType.SEGMENTATION.value == "segmentation"
        assert AITaskType.TRIAGE.value == "triage"
        assert len(AITaskType) == 8

    def test_model_architecture_values(self):
        assert ModelArchitecture.UNET_3D.value == "3d_unet"
        assert ModelArchitecture.VISTA3D.value == "vista3d"
        assert ModelArchitecture.VILAM3.value == "vila_m3"
        assert len(ModelArchitecture) == 11

    def test_evidence_level_values(self):
        assert EvidenceLevel.VALIDATED.value == "validated"
        assert EvidenceLevel.EMERGING.value == "emerging"
        assert EvidenceLevel.EXPLORATORY.value == "exploratory"
        assert len(EvidenceLevel) == 3

    def test_lung_rads_values(self):
        assert LungRADS.CAT_0.value == "0"
        assert LungRADS.CAT_4A.value == "4A"
        assert LungRADS.CAT_4B.value == "4B"
        assert LungRADS.CAT_4X.value == "4X"
        assert len(LungRADS) == 8

    def test_bi_rads_values(self):
        assert BiRADS.CAT_0.value == "0"
        assert BiRADS.CAT_5.value == "5"
        assert BiRADS.CAT_6.value == "6"
        assert len(BiRADS) == 7

    def test_nim_service_status_values(self):
        assert NIMServiceStatus.AVAILABLE.value == "available"
        assert NIMServiceStatus.MOCK.value == "mock"
        assert len(NIMServiceStatus) == 3

    def test_workflow_status_values(self):
        assert WorkflowStatus.PENDING.value == "pending"
        assert WorkflowStatus.COMPLETED.value == "completed"
        assert WorkflowStatus.FAILED.value == "failed"
        assert len(WorkflowStatus) == 4


# ===================================================================
# COLLECTION MODEL CREATION
# ===================================================================


class TestImagingLiterature:
    """Tests for the ImagingLiterature collection model."""

    def test_create_with_valid_data(self):
        lit = ImagingLiterature(
            id="PMID:38001234",
            title="Deep Learning for ICH Detection on CT",
            text_chunk="We developed a 3D ResNet for hemorrhage classification.",
            year=2024,
        )
        assert lit.id == "PMID:38001234"
        assert lit.year == 2024
        assert lit.source_type == SourceType.PUBMED

    def test_to_embedding_text_non_empty(self):
        lit = ImagingLiterature(
            id="PMID:001",
            title="VISTA-3D Segmentation Study",
            text_chunk="Universal segmentation of 127 structures.",
            year=2025,
            modality=ImagingModality.CT,
            disease="lung cancer",
        )
        text = lit.to_embedding_text()
        assert len(text) > 0
        assert "VISTA-3D" in text
        assert "lung cancer" in text

    def test_defaults(self):
        lit = ImagingLiterature(
            id="1", title="T", text_chunk="C", year=2024,
        )
        assert lit.source_type == SourceType.PUBMED
        assert lit.modality == ImagingModality.CT
        assert lit.body_region == BodyRegion.CHEST


class TestImagingTrial:
    """Tests for the ImagingTrial collection model."""

    def test_create_with_valid_data(self):
        trial = ImagingTrial(
            id="NCT06001234",
            title="AI Lung Cancer Screening Trial",
            text_summary="Phase 3 trial evaluating AI triage.",
        )
        assert trial.id == "NCT06001234"
        assert trial.phase == TrialPhase.NA

    def test_to_embedding_text_includes_disease(self):
        trial = ImagingTrial(
            id="NCT06005678",
            title="CT Head AI Triage",
            text_summary="Evaluating AI for ICH detection.",
            disease="intracranial hemorrhage",
            outcome_summary="AUC 0.97",
        )
        text = trial.to_embedding_text()
        assert "intracranial hemorrhage" in text
        assert "AUC 0.97" in text


class TestImagingFinding:
    """Tests for the ImagingFinding collection model."""

    def test_create_with_valid_data(self):
        finding = ImagingFinding(
            id="finding-001",
            text_summary="Intraparenchymal hemorrhage in right basal ganglia.",
            finding_category=FindingCategory.HEMORRHAGE,
            severity=FindingSeverity.URGENT,
        )
        assert finding.finding_category == FindingCategory.HEMORRHAGE
        assert finding.severity == FindingSeverity.URGENT

    def test_to_embedding_text_non_empty(self):
        finding = ImagingFinding(
            id="f-002",
            text_summary="Solid nodule in RUL, 12mm.",
            clinical_significance="Suspicious for malignancy",
            differential_diagnosis="Primary lung cancer vs metastasis",
        )
        text = finding.to_embedding_text()
        assert "Suspicious" in text
        assert "DDx:" in text


class TestImagingProtocol:
    """Tests for the ImagingProtocol collection model."""

    def test_create_with_valid_data(self):
        proto = ImagingProtocol(
            id="proto-001",
            text_summary="Non-contrast CT head for hemorrhage evaluation.",
            protocol_name="CT Head Hemorrhage",
        )
        assert proto.protocol_name == "CT Head Hemorrhage"

    def test_to_embedding_text_includes_contrast(self):
        proto = ImagingProtocol(
            id="p-002",
            text_summary="CTPA for pulmonary embolism.",
            contrast_agent="iodinated",
            clinical_indication="rule out PE",
        )
        text = proto.to_embedding_text()
        assert "iodinated" in text
        assert "rule out PE" in text


class TestImagingDevice:
    """Tests for the ImagingDevice collection model."""

    def test_create_with_valid_data(self):
        device = ImagingDevice(
            id="K223456",
            text_summary="AI-based ICH detection on CT.",
            device_name="BrainScan AI",
            manufacturer="Acme Medical",
            regulatory_status=DeviceRegulatory.CLEARED_510K,
        )
        assert device.regulatory_status == DeviceRegulatory.CLEARED_510K

    def test_to_embedding_text_non_empty(self):
        device = ImagingDevice(
            id="K999",
            text_summary="Lung nodule CADe.",
            device_name="LungAI",
            intended_use="Assist in detecting pulmonary nodules",
            performance_summary="AUC 0.94",
        )
        text = device.to_embedding_text()
        assert "LungAI" in text
        assert "AUC 0.94" in text


class TestAnatomyRecord:
    """Tests for the AnatomyRecord collection model."""

    def test_create_with_valid_data(self):
        anatomy = AnatomyRecord(
            id="anat-brain",
            text_summary="Brain anatomy with cortex and deep structures.",
            structure_name="Brain",
            system="central nervous system",
        )
        assert anatomy.structure_name == "Brain"

    def test_to_embedding_text_non_empty(self):
        anatomy = AnatomyRecord(
            id="anat-liver",
            text_summary="Liver with 8 Couinaud segments.",
            structure_name="Liver",
            common_pathologies="HCC, metastasis, hemangioma",
        )
        text = anatomy.to_embedding_text()
        assert "Liver" in text
        assert "HCC" in text


class TestBenchmarkRecord:
    """Tests for the BenchmarkRecord collection model."""

    def test_create_with_valid_data(self):
        bench = BenchmarkRecord(
            id="bench-001",
            text_summary="VISTA-3D on TotalSegmentator.",
            model_name="VISTA-3D",
            metric_name="Dice",
            metric_value=0.91,
        )
        assert bench.metric_value == 0.91

    def test_to_embedding_text_includes_metric(self):
        bench = BenchmarkRecord(
            id="b-002",
            text_summary="nnU-Net lung nodule segmentation.",
            model_name="nnU-Net",
            metric_name="AUC",
            metric_value=0.95,
            hardware="DGX Spark",
        )
        text = bench.to_embedding_text()
        assert "nnU-Net" in text
        assert "0.95" in text
        assert "DGX Spark" in text


class TestGuidelineRecord:
    """Tests for the GuidelineRecord collection model."""

    def test_create_with_valid_data(self):
        guide = GuidelineRecord(
            id="guide-001",
            text_summary="ACR Lung-RADS v2022 assessment categories.",
            guideline_name="ACR Lung-RADS v2022",
            organization="ACR",
            year=2022,
        )
        assert guide.organization == "ACR"

    def test_to_embedding_text_non_empty(self):
        guide = GuidelineRecord(
            id="g-002",
            text_summary="Fleischner Society nodule management.",
            guideline_name="Fleischner 2017",
            key_recommendation="Follow-up CT in 6-12 months for 6-8mm nodule.",
        )
        text = guide.to_embedding_text()
        assert "Fleischner" in text
        assert "Follow-up" in text


class TestReportTemplate:
    """Tests for the ReportTemplate collection model."""

    def test_create_with_valid_data(self):
        tmpl = ReportTemplate(
            id="tmpl-001",
            text_summary="Structured chest CT report template.",
            template_name="Chest CT Normal",
        )
        assert tmpl.template_name == "Chest CT Normal"

    def test_to_embedding_text_non_empty(self):
        tmpl = ReportTemplate(
            id="t-002",
            text_summary="Pulmonary nodule structured report.",
            finding_type="pulmonary_nodule",
            coding_system="RadLex",
        )
        text = tmpl.to_embedding_text()
        assert "pulmonary_nodule" in text
        assert "RadLex" in text


class TestDatasetRecord:
    """Tests for the DatasetRecord collection model."""

    def test_create_with_valid_data(self):
        ds = DatasetRecord(
            id="ds-001",
            text_summary="RSNA intracranial hemorrhage detection dataset.",
            dataset_name="RSNA ICH Detection",
            source="RSNA",
            num_studies=25000,
        )
        assert ds.num_studies == 25000

    def test_to_embedding_text_non_empty(self):
        ds = DatasetRecord(
            id="d-002",
            text_summary="LIDC-IDRI lung nodule dataset.",
            dataset_name="LIDC-IDRI",
            disease_labels="lung nodule",
            annotation_type="segmentation_mask",
        )
        text = ds.to_embedding_text()
        assert "LIDC-IDRI" in text
        assert "lung nodule" in text


# ===================================================================
# PARAMETRIZED: ALL 10 MODELS PRODUCE NON-EMPTY EMBEDDING TEXT
# ===================================================================


@pytest.mark.parametrize(
    "model_cls,kwargs",
    [
        (ImagingLiterature, {"id": "1", "title": "T", "text_chunk": "C", "year": 2024}),
        (ImagingTrial, {"id": "NCT00000001", "title": "T", "text_summary": "S"}),
        (ImagingFinding, {"id": "f1", "text_summary": "S"}),
        (ImagingProtocol, {"id": "p1", "text_summary": "S"}),
        (ImagingDevice, {"id": "d1", "text_summary": "S"}),
        (AnatomyRecord, {"id": "a1", "text_summary": "S"}),
        (BenchmarkRecord, {"id": "b1", "text_summary": "S"}),
        (GuidelineRecord, {"id": "g1", "text_summary": "S"}),
        (ReportTemplate, {"id": "t1", "text_summary": "S"}),
        (DatasetRecord, {"id": "ds1", "text_summary": "S"}),
    ],
    ids=[
        "ImagingLiterature",
        "ImagingTrial",
        "ImagingFinding",
        "ImagingProtocol",
        "ImagingDevice",
        "AnatomyRecord",
        "BenchmarkRecord",
        "GuidelineRecord",
        "ReportTemplate",
        "DatasetRecord",
    ],
)
def test_all_models_embedding_text(model_cls, kwargs):
    """Every collection model's to_embedding_text() returns a non-empty string."""
    instance = model_cls(**kwargs)
    text = instance.to_embedding_text()
    assert isinstance(text, str)
    assert len(text) > 0


# ===================================================================
# NIM RESULT MODELS
# ===================================================================


class TestSegmentationResult:
    """Tests for SegmentationResult (VISTA-3D output)."""

    def test_create_default(self):
        result = SegmentationResult()
        assert result.classes_detected == []
        assert result.volumes == {}
        assert result.model == "vista3d"

    def test_create_with_data(self):
        result = SegmentationResult(
            classes_detected=["liver", "spleen"],
            volumes={"liver": 1500.0, "spleen": 200.0},
            num_classes=2,
            inference_time_ms=3000.0,
            is_mock=True,
        )
        assert len(result.classes_detected) == 2
        assert result.volumes["liver"] == 1500.0
        assert result.is_mock is True


class TestSyntheticCTResult:
    """Tests for SyntheticCTResult (MAISI output)."""

    def test_create_default(self):
        result = SyntheticCTResult()
        assert result.resolution == "512x512x512"
        assert result.model == "maisi"

    def test_create_with_data(self):
        result = SyntheticCTResult(
            resolution="256x256x256",
            body_region="chest",
            num_classes_annotated=127,
            generation_time_ms=45000.0,
            is_mock=True,
        )
        assert result.body_region == "chest"
        assert result.num_classes_annotated == 127


class TestVLMResponse:
    """Tests for VLMResponse (VILA-M3 output)."""

    def test_create_default(self):
        result = VLMResponse()
        assert result.answer == ""
        assert result.findings == []
        assert result.model == "vila_m3"

    def test_create_with_data(self):
        result = VLMResponse(
            answer="No acute abnormality.",
            findings=["clear lungs", "normal heart"],
            confidence=0.85,
            is_mock=True,
        )
        assert len(result.findings) == 2
        assert result.confidence == 0.85


class TestWorkflowResult:
    """Tests for WorkflowResult."""

    def test_create_with_required_only(self):
        result = WorkflowResult(workflow_name="test_workflow")
        assert result.workflow_name == "test_workflow"
        assert result.status == WorkflowStatus.COMPLETED
        assert result.severity == FindingSeverity.ROUTINE

    def test_create_with_all_fields(self):
        result = WorkflowResult(
            workflow_name="ct_head_hemorrhage",
            status=WorkflowStatus.COMPLETED,
            findings=[{"category": "hemorrhage", "severity": "critical"}],
            measurements={"volume_ml": 35.0, "midline_shift_mm": 6.0},
            classification="critical_hemorrhage",
            severity=FindingSeverity.CRITICAL,
            nim_services_used=["3D U-Net", "VISTA-3D"],
            is_mock=True,
        )
        assert result.severity == FindingSeverity.CRITICAL
        assert len(result.findings) == 1
        assert result.measurements["volume_ml"] == 35.0


# ===================================================================
# SEARCH RESULT MODELS
# ===================================================================


class TestSearchHit:
    """Tests for SearchHit creation and fields."""

    def test_create_search_hit(self):
        hit = SearchHit(
            collection="imaging_literature",
            id="PMID:001",
            score=0.85,
            text="AI hemorrhage detection study.",
            metadata={"year": 2024},
        )
        assert hit.collection == "imaging_literature"
        assert hit.score == 0.85
        assert hit.metadata["year"] == 2024

    def test_default_metadata_is_empty_dict(self):
        hit = SearchHit(
            collection="imaging_trials",
            id="NCT00000001",
            score=0.5,
            text="test",
        )
        assert hit.metadata == {}


class TestCrossCollectionResult:
    """Tests for CrossCollectionResult creation and properties."""

    def test_create_cross_collection_result(self, sample_search_hits):
        result = CrossCollectionResult(
            query="test query",
            hits=sample_search_hits,
            total_collections_searched=11,
            search_time_ms=50.0,
        )
        assert result.query == "test query"
        assert result.hit_count == 5
        assert result.total_collections_searched == 11

    def test_hits_by_collection(self, sample_search_hits):
        result = CrossCollectionResult(query="test", hits=sample_search_hits)
        grouped = result.hits_by_collection()
        assert "imaging_literature" in grouped
        assert "imaging_trials" in grouped
        assert len(grouped["imaging_literature"]) == 1

    def test_empty_result(self):
        result = CrossCollectionResult(query="empty")
        assert result.hit_count == 0
        assert result.hits_by_collection() == {}


class TestComparativeResult:
    """Tests for ComparativeResult creation and properties."""

    def test_create_comparative_result(self):
        ev_a = CrossCollectionResult(
            query="CT",
            hits=[SearchHit(collection="imaging_literature", id="1", score=0.9, text="A")],
        )
        ev_b = CrossCollectionResult(
            query="MRI",
            hits=[
                SearchHit(collection="imaging_literature", id="2", score=0.85, text="B"),
                SearchHit(collection="imaging_trials", id="NCT00000001", score=0.8, text="C"),
            ],
        )
        comp = ComparativeResult(
            query="Compare CT vs MRI for hemorrhage",
            entity_a="CT",
            entity_b="MRI",
            evidence_a=ev_a,
            evidence_b=ev_b,
        )
        assert comp.total_hits == 3
        assert comp.entity_a == "CT"
        assert comp.entity_b == "MRI"


# ===================================================================
# AGENT MODELS
# ===================================================================


class TestAgentQuery:
    """Tests for the AgentQuery input model."""

    def test_create_with_required_only(self):
        query = AgentQuery(question="What causes hemorrhage on CT?")
        assert query.question == "What causes hemorrhage on CT?"
        assert query.modality is None
        assert query.body_region is None
        assert query.include_genomic is True
        assert query.include_nim is True

    def test_create_with_all_fields(self):
        query = AgentQuery(
            question="Lung nodule classification",
            modality=ImagingModality.CT,
            body_region=BodyRegion.CHEST,
            include_genomic=False,
            include_nim=True,
        )
        assert query.modality == ImagingModality.CT
        assert query.body_region == BodyRegion.CHEST
        assert query.include_genomic is False


class TestAgentResponse:
    """Tests for the AgentResponse output model."""

    def test_create_with_required_fields(self):
        evidence = CrossCollectionResult(query="test", hits=[])
        response = AgentResponse(
            question="test question",
            answer="test answer",
            evidence=evidence,
        )
        assert response.question == "test question"
        assert response.workflow_results == []
        assert response.nim_services_used == []
        assert len(response.timestamp) > 0

    def test_create_with_all_fields(self, sample_search_hits):
        evidence = CrossCollectionResult(
            query="test",
            hits=sample_search_hits,
            total_collections_searched=11,
        )
        response = AgentResponse(
            question="Hemorrhage detection AI",
            answer="3D U-Net achieves 95% sensitivity.",
            evidence=evidence,
            workflow_results=[WorkflowResult(workflow_name="ct_head_hemorrhage")],
            nim_services_used=["vista3d"],
            knowledge_used=["pathology"],
        )
        assert len(response.workflow_results) == 1
        assert "vista3d" in response.nim_services_used
