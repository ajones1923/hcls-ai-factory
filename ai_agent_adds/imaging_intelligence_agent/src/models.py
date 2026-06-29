"""Pydantic data models for Imaging Intelligence Agent.

Maps to the 10 Milvus collections + NIM service results.
Follows the same dataclass/Pydantic pattern as:
  - rag-chat-pipeline/src/vcf_parser.py (VariantEvidence)
  - drug-discovery-pipeline/src/models.py (GeneratedMolecule, DockingResult)
  - ai_agent_adds/cart_intelligence_agent/src/models.py (CAR-T Agent)
"""

from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# ENUMS
# ═══════════════════════════════════════════════════════════════════════


class ImagingModality(str, Enum):
    """Supported medical imaging modalities."""
    CT = "ct"
    MRI = "mri"
    XRAY = "xray"
    CXR = "cxr"
    ULTRASOUND = "ultrasound"
    PET = "pet"
    PET_CT = "pet_ct"
    MAMMOGRAPHY = "mammography"
    FLUOROSCOPY = "fluoroscopy"


class BodyRegion(str, Enum):
    """Anatomical body regions for imaging studies."""
    HEAD = "head"
    NECK = "neck"
    CHEST = "chest"
    ABDOMEN = "abdomen"
    PELVIS = "pelvis"
    SPINE = "spine"
    EXTREMITY = "extremity"
    WHOLE_BODY = "whole_body"
    BRAIN = "brain"
    CARDIAC = "cardiac"
    BREAST = "breast"


class FindingSeverity(str, Enum):
    """Clinical severity of an imaging finding."""
    CRITICAL = "critical"
    URGENT = "urgent"
    SIGNIFICANT = "significant"
    ROUTINE = "routine"
    NORMAL = "normal"


class FindingCategory(str, Enum):
    """Categories of imaging findings."""
    HEMORRHAGE = "hemorrhage"
    NODULE = "nodule"
    MASS = "mass"
    FRACTURE = "fracture"
    CONSOLIDATION = "consolidation"
    EFFUSION = "effusion"
    PNEUMOTHORAX = "pneumothorax"
    EDEMA = "edema"
    ATELECTASIS = "atelectasis"
    LESION = "lesion"
    INFARCT = "infarct"
    STENOSIS = "stenosis"
    CALCIFICATION = "calcification"
    NORMAL = "normal"


class SourceType(str, Enum):
    PUBMED = "pubmed"
    PMC = "pmc"
    PREPRINT = "preprint"
    GUIDELINE = "guideline"
    MANUAL = "manual"


class TrialPhase(str, Enum):
    EARLY_1 = "Early Phase 1"
    PHASE_1 = "Phase 1"
    PHASE_1_2 = "Phase 1/Phase 2"
    PHASE_2 = "Phase 2"
    PHASE_2_3 = "Phase 2/Phase 3"
    PHASE_3 = "Phase 3"
    PHASE_4 = "Phase 4"
    NA = "N/A"


class TrialStatus(str, Enum):
    RECRUITING = "Recruiting"
    ACTIVE = "Active, not recruiting"
    COMPLETED = "Completed"
    TERMINATED = "Terminated"
    WITHDRAWN = "Withdrawn"
    SUSPENDED = "Suspended"
    NOT_YET = "Not yet recruiting"
    UNKNOWN = "Unknown status"


class DeviceRegulatory(str, Enum):
    """FDA regulatory pathway for AI/ML medical devices."""
    CLEARED_510K = "510k_cleared"
    DE_NOVO = "de_novo"
    PMA = "pma"
    BREAKTHROUGH = "breakthrough"
    PENDING = "pending"
    NOT_SUBMITTED = "not_submitted"


class AITaskType(str, Enum):
    """Types of AI tasks in medical imaging."""
    DETECTION = "detection"
    SEGMENTATION = "segmentation"
    CLASSIFICATION = "classification"
    REGISTRATION = "registration"
    RECONSTRUCTION = "reconstruction"
    TRIAGE = "triage"
    QUANTIFICATION = "quantification"
    REPORT_GENERATION = "report_generation"


class ModelArchitecture(str, Enum):
    """Neural network architectures for medical imaging AI."""
    UNET_3D = "3d_unet"
    SEGRESNET = "segresnet"
    RETINANET = "retinanet"
    DENSENET = "densenet"
    RESNET = "resnet"
    VISION_TRANSFORMER = "vit"
    SWIN_TRANSFORMER = "swin"
    NNUNET = "nnunet"
    VISTA3D = "vista3d"
    MONAI_BUNDLE = "monai_bundle"
    VILAM3 = "vila_m3"


class EvidenceLevel(str, Enum):
    VALIDATED = "validated"
    EMERGING = "emerging"
    EXPLORATORY = "exploratory"


class LungRADS(str, Enum):
    """ACR Lung-RADS v2022 categories."""
    CAT_0 = "0"
    CAT_1 = "1"
    CAT_2 = "2"
    CAT_3 = "3"
    CAT_4A = "4A"
    CAT_4B = "4B"
    CAT_4X = "4X"
    CAT_S = "S"


class BiRADS(str, Enum):
    """ACR BI-RADS categories."""
    CAT_0 = "0"
    CAT_1 = "1"
    CAT_2 = "2"
    CAT_3 = "3"
    CAT_4 = "4"
    CAT_5 = "5"
    CAT_6 = "6"


class NIMServiceStatus(str, Enum):
    """Status of a NVIDIA NIM microservice."""
    AVAILABLE = "available"
    UNAVAILABLE = "unavailable"
    MOCK = "mock"


class WorkflowStatus(str, Enum):
    """Status of an imaging analysis workflow."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MODELS (map to Milvus schemas)
# ═══════════════════════════════════════════════════════════════════════


class ImagingLiterature(BaseModel):
    """Published research paper — maps to imaging_literature collection."""
    id: str = Field(..., description="PMID or DOI")
    title: str = Field(..., max_length=500)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    source_type: SourceType = SourceType.PUBMED
    year: int = Field(..., ge=1990, le=2030)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    ai_task: str = Field("", max_length=100)
    disease: str = Field("", max_length=200)
    keywords: str = Field("", max_length=1000)
    journal: str = Field("", max_length=200)

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [self.title]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.modality:
            parts.append(f"Modality: {self.modality.value}")
        if self.body_region:
            parts.append(f"Region: {self.body_region.value}")
        if self.disease:
            parts.append(f"Disease: {self.disease}")
        return " ".join(parts)


class ImagingTrial(BaseModel):
    """ClinicalTrials.gov record — maps to imaging_trials collection."""
    id: str = Field(..., description="NCT number", pattern=r"^NCT\d{8}$")
    title: str = Field(..., max_length=500)
    text_summary: str = Field(..., max_length=3000)
    phase: TrialPhase = TrialPhase.NA
    status: TrialStatus = TrialStatus.UNKNOWN
    sponsor: str = Field("", max_length=200)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    ai_task: str = Field("", max_length=100)
    disease: str = Field("", max_length=200)
    enrollment: int = Field(0, ge=0)
    start_year: int = Field(0, ge=0, le=2030)
    outcome_summary: str = Field("", max_length=2000)

    def to_embedding_text(self) -> str:
        parts = [self.title, self.text_summary]
        if self.modality:
            parts.append(f"Modality: {self.modality.value}")
        if self.disease:
            parts.append(f"Indication: {self.disease}")
        if self.outcome_summary:
            parts.append(f"Outcome: {self.outcome_summary}")
        return " ".join(parts)


class ImagingFinding(BaseModel):
    """Imaging finding record — maps to imaging_findings collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=3000)
    finding_category: FindingCategory = FindingCategory.NORMAL
    severity: FindingSeverity = FindingSeverity.ROUTINE
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    clinical_significance: str = Field("", max_length=500)
    differential_diagnosis: str = Field("", max_length=500)
    recommended_followup: str = Field("", max_length=500)
    measurement_type: str = Field("", max_length=100, description="e.g., diameter, volume, HU")
    measurement_value: str = Field("", max_length=100, description="e.g., 12mm, 45 HU")
    classification_system: str = Field("", max_length=50, description="e.g., Lung-RADS, BI-RADS")
    classification_score: str = Field("", max_length=20, description="e.g., 4A, 3")

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.finding_category:
            parts.append(f"Finding: {self.finding_category.value}")
        if self.severity:
            parts.append(f"Severity: {self.severity.value}")
        if self.clinical_significance:
            parts.append(f"Significance: {self.clinical_significance}")
        if self.differential_diagnosis:
            parts.append(f"DDx: {self.differential_diagnosis}")
        return " ".join(parts)


class ImagingProtocol(BaseModel):
    """Imaging protocol / acquisition parameters — maps to imaging_protocols collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=2000)
    protocol_name: str = Field("", max_length=200)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    contrast_agent: str = Field("", max_length=100, description="e.g., iodinated, gadolinium, none")
    slice_thickness_mm: str = Field("", max_length=20, description="e.g., 1.0, 3.0")
    radiation_dose: str = Field("", max_length=100, description="e.g., CTDIvol 8 mGy")
    scan_duration: str = Field("", max_length=50, description="e.g., 30 seconds")
    clinical_indication: str = Field("", max_length=500)
    preprocessing_steps: str = Field("", max_length=500, description="e.g., windowing, resampling")

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.protocol_name:
            parts.append(f"Protocol: {self.protocol_name}")
        if self.modality:
            parts.append(f"Modality: {self.modality.value}")
        if self.clinical_indication:
            parts.append(f"Indication: {self.clinical_indication}")
        if self.contrast_agent:
            parts.append(f"Contrast: {self.contrast_agent}")
        return " ".join(parts)


class ImagingDevice(BaseModel):
    """FDA-cleared AI/ML medical device — maps to imaging_devices collection."""
    id: str = Field(..., description="510(k) number or device identifier")
    text_summary: str = Field(..., max_length=3000)
    device_name: str = Field("", max_length=200)
    manufacturer: str = Field("", max_length=200)
    regulatory_status: DeviceRegulatory = DeviceRegulatory.PENDING
    clearance_date: str = Field("", max_length=20, description="YYYY-MM-DD or YYYY-MM")
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    ai_task: str = Field("", max_length=100)
    intended_use: str = Field("", max_length=500)
    performance_summary: str = Field("", max_length=1000)
    model_architecture: str = Field("", max_length=100)

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.device_name:
            parts.append(f"Device: {self.device_name}")
        if self.manufacturer:
            parts.append(f"Manufacturer: {self.manufacturer}")
        if self.intended_use:
            parts.append(f"Use: {self.intended_use}")
        if self.performance_summary:
            parts.append(f"Performance: {self.performance_summary}")
        return " ".join(parts)


class AnatomyRecord(BaseModel):
    """Anatomical structure reference — maps to imaging_anatomy collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=2000)
    structure_name: str = Field("", max_length=200)
    body_region: BodyRegion = BodyRegion.CHEST
    system: str = Field("", max_length=100, description="e.g., cardiovascular, respiratory")
    snomed_code: str = Field("", max_length=50)
    fma_id: str = Field("", max_length=50, description="Foundational Model of Anatomy ID")
    imaging_characteristics: str = Field("", max_length=500)
    common_pathologies: str = Field("", max_length=500)
    segmentation_label_id: int = Field(0, ge=0, description="VISTA-3D segmentation label ID for Phase 2")

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.structure_name:
            parts.append(f"Structure: {self.structure_name}")
        if self.system:
            parts.append(f"System: {self.system}")
        if self.imaging_characteristics:
            parts.append(f"Characteristics: {self.imaging_characteristics}")
        if self.common_pathologies:
            parts.append(f"Pathologies: {self.common_pathologies}")
        return " ".join(parts)


class BenchmarkRecord(BaseModel):
    """Model benchmark / performance record — maps to imaging_benchmarks collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=2000)
    model_name: str = Field("", max_length=200)
    model_architecture: ModelArchitecture = ModelArchitecture.UNET_3D
    ai_task: str = Field("", max_length=100)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    dataset_name: str = Field("", max_length=200)
    metric_name: str = Field("", max_length=100, description="e.g., Dice, AUC, sensitivity")
    metric_value: float = Field(0.0)
    training_data_size: str = Field("", max_length=100, description="e.g., 10000 images")
    inference_time_ms: str = Field("", max_length=50)
    hardware: str = Field("", max_length=100, description="e.g., A100, DGX Spark")

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.model_name:
            parts.append(f"Model: {self.model_name}")
        if self.dataset_name:
            parts.append(f"Dataset: {self.dataset_name}")
        if self.metric_name and self.metric_value:
            parts.append(f"{self.metric_name}: {self.metric_value}")
        if self.hardware:
            parts.append(f"Hardware: {self.hardware}")
        return " ".join(parts)


class GuidelineRecord(BaseModel):
    """Clinical practice guideline — maps to imaging_guidelines collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=3000)
    guideline_name: str = Field("", max_length=300)
    organization: str = Field("", max_length=100, description="e.g., ACR, RSNA, NCCN")
    year: int = Field(0, ge=0, le=2030)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    clinical_indication: str = Field("", max_length=500)
    classification_system: str = Field("", max_length=100, description="e.g., Lung-RADS, BI-RADS, LI-RADS")
    key_recommendation: str = Field("", max_length=1000)
    evidence_level: EvidenceLevel = EvidenceLevel.EMERGING

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.guideline_name:
            parts.append(f"Guideline: {self.guideline_name}")
        if self.organization:
            parts.append(f"Organization: {self.organization}")
        if self.key_recommendation:
            parts.append(f"Recommendation: {self.key_recommendation}")
        if self.classification_system:
            parts.append(f"Classification: {self.classification_system}")
        return " ".join(parts)


class ReportTemplate(BaseModel):
    """Structured radiology report template — maps to imaging_report_templates collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=2000)
    template_name: str = Field("", max_length=200)
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    finding_type: str = Field("", max_length=100, description="e.g., pulmonary_nodule, fracture")
    structured_fields: str = Field("", max_length=1000, description="Comma-separated field names")
    example_report: str = Field("", max_length=2000)
    coding_system: str = Field("", max_length=50, description="e.g., RadLex, SNOMED-CT, ICD-10")

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.template_name:
            parts.append(f"Template: {self.template_name}")
        if self.finding_type:
            parts.append(f"Finding type: {self.finding_type}")
        if self.coding_system:
            parts.append(f"Coding: {self.coding_system}")
        if self.example_report:
            parts.append(f"Example: {self.example_report[:200]}")
        return " ".join(parts)


class DatasetRecord(BaseModel):
    """Public imaging dataset — maps to imaging_datasets collection."""
    id: str = Field(..., max_length=100)
    text_summary: str = Field(..., max_length=2000)
    dataset_name: str = Field("", max_length=200)
    source: str = Field("", max_length=100, description="e.g., TCIA, PhysioNet, Kaggle")
    modality: ImagingModality = ImagingModality.CT
    body_region: BodyRegion = BodyRegion.CHEST
    num_studies: int = Field(0, ge=0)
    num_images: int = Field(0, ge=0)
    disease_labels: str = Field("", max_length=500)
    annotation_type: str = Field("", max_length=200, description="e.g., bounding_box, segmentation_mask, report")
    license_type: str = Field("", max_length=100, description="e.g., CC-BY-4.0, TCIA restricted")
    download_url: str = Field("", max_length=500)

    def to_embedding_text(self) -> str:
        parts = [self.text_summary]
        if self.dataset_name:
            parts.append(f"Dataset: {self.dataset_name}")
        if self.source:
            parts.append(f"Source: {self.source}")
        if self.disease_labels:
            parts.append(f"Labels: {self.disease_labels}")
        if self.annotation_type:
            parts.append(f"Annotations: {self.annotation_type}")
        return " ".join(parts)


# ═══════════════════════════════════════════════════════════════════════
# NIM RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════


class SegmentationResult(BaseModel):
    """Result from VISTA-3D NIM segmentation."""
    classes_detected: List[str] = Field(default_factory=list)
    volumes: Dict[str, float] = Field(default_factory=dict)  # class -> volume_ml
    num_classes: int = 0
    inference_time_ms: float = 0.0
    segmentation_mask_path: Optional[str] = None
    model_name: str = "vista3d"
    model: str = "vista3d"
    is_mock: bool = False


class SyntheticCTResult(BaseModel):
    """Result from MAISI NIM synthetic CT generation."""
    volume_path: str = ""
    segmentation_mask_path: Optional[str] = None
    resolution: str = "512x512x512"
    body_region: str = ""
    num_classes_annotated: int = 0
    generation_time_ms: float = 0.0
    voxel_spacing_mm: List[float] = Field(default_factory=lambda: [1.0, 1.0, 1.0])
    model_name: str = "maisi"
    model: str = "maisi"
    is_mock: bool = False


class VLMResponse(BaseModel):
    """Result from VILA-M3 vision-language model."""
    answer: str = ""
    findings: List[str] = Field(default_factory=list)
    confidence: float = 0.0
    inference_time_ms: float = 0.0
    model: str = "vila_m3"
    is_mock: bool = False


class WorkflowResult(BaseModel):
    """Result from a reference imaging workflow."""
    workflow_name: str
    status: WorkflowStatus = WorkflowStatus.COMPLETED
    findings: List[Dict] = Field(default_factory=list)
    measurements: Dict[str, float] = Field(default_factory=dict)
    classification: str = ""
    severity: FindingSeverity = FindingSeverity.ROUTINE
    inference_time_ms: float = 0.0
    nim_services_used: List[str] = Field(default_factory=list)
    is_mock: bool = False


# ═══════════════════════════════════════════════════════════════════════
# SEARCH RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════


class SearchHit(BaseModel):
    """A single search result from any collection."""
    collection: str
    id: str
    score: float = Field(..., ge=0.0)
    text: str
    metadata: Dict = Field(default_factory=dict)


class CrossCollectionResult(BaseModel):
    """Merged results from multi-collection search."""
    query: str
    hits: List[SearchHit] = Field(default_factory=list)
    knowledge_context: str = ""
    total_collections_searched: int = 0
    search_time_ms: float = 0.0

    @property
    def hit_count(self) -> int:
        return len(self.hits)

    def hits_by_collection(self) -> Dict[str, List[SearchHit]]:
        grouped: Dict[str, List[SearchHit]] = {}
        for hit in self.hits:
            grouped.setdefault(hit.collection, []).append(hit)
        return grouped


class ComparativeResult(BaseModel):
    """Results from a comparative analysis query."""
    query: str
    entity_a: str
    entity_b: str
    evidence_a: CrossCollectionResult
    evidence_b: CrossCollectionResult
    comparison_context: str = ""
    total_search_time_ms: float = 0.0

    @property
    def total_hits(self) -> int:
        return self.evidence_a.hit_count + self.evidence_b.hit_count


# ═══════════════════════════════════════════════════════════════════════
# AGENT MODELS
# ═══════════════════════════════════════════════════════════════════════


class CrossModalResult(BaseModel):
    """Result from cross-modal genomics query triggered by imaging findings.

    When an imaging workflow produces a high-severity result (e.g.,
    Lung-RADS 4A+), the cross-modal trigger automatically queries the
    genomic_evidence Milvus collection for relevant cancer genomics
    context, enriching imaging findings with molecular insights.
    """
    trigger_reason: str = ""
    genomic_context: List[str] = Field(default_factory=list)
    genomic_hit_count: int = 0
    query_count: int = 0
    enrichment_summary: str = ""


class AgentQuery(BaseModel):
    """Input to the Imaging Intelligence Agent."""
    question: str
    modality: Optional[ImagingModality] = None
    body_region: Optional[BodyRegion] = None
    include_genomic: bool = True
    include_nim: bool = True  # Whether to invoke NIM services


class AgentResponse(BaseModel):
    """Output from the Imaging Intelligence Agent."""
    question: str
    answer: str
    evidence: CrossCollectionResult
    workflow_results: List[WorkflowResult] = Field(default_factory=list)
    nim_services_used: List[str] = Field(default_factory=list)
    knowledge_used: List[str] = Field(default_factory=list)
    cross_modal: Optional[CrossModalResult] = None
    timestamp: str = Field(default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
