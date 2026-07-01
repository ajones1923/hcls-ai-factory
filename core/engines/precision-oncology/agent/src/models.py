"""
Oncology Intelligence Agent - Pydantic Models
================================================
Domain models, search containers, and agent I/O types
for the Oncology Intelligence RAG agent.
"""

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, ConfigDict, Field


# ═══════════════════════════════════════════════════════════════════════════
#  Enums
# ═══════════════════════════════════════════════════════════════════════════


class CancerType(str, Enum):
    """Major cancer types supported by the agent."""
    NSCLC = "nsclc"
    SCLC = "sclc"
    BREAST = "breast"
    COLORECTAL = "colorectal"
    MELANOMA = "melanoma"
    PANCREATIC = "pancreatic"
    OVARIAN = "ovarian"
    PROSTATE = "prostate"
    RENAL = "renal"
    BLADDER = "bladder"
    HEAD_NECK = "head_neck"
    HEPATOCELLULAR = "hepatocellular"
    GASTRIC = "gastric"
    GLIOBLASTOMA = "glioblastoma"
    AML = "aml"
    CML = "cml"
    ALL = "all"
    CLL = "cll"
    DLBCL = "dlbcl"
    MULTIPLE_MYELOMA = "multiple_myeloma"
    OTHER = "other"
    CHOLANGIOCARCINOMA = "cholangiocarcinoma"
    ENDOMETRIAL = "endometrial"
    THYROID = "thyroid"
    MESOTHELIOMA = "mesothelioma"
    SARCOMA = "sarcoma"


class VariantType(str, Enum):
    """Genomic variant classification."""
    SNV = "snv"
    INDEL = "indel"
    CNV_AMP = "cnv_amplification"
    CNV_DEL = "cnv_deletion"
    FUSION = "fusion"
    REARRANGEMENT = "rearrangement"
    SV = "structural_variant"


class EvidenceLevel(str, Enum):
    """Tiered evidence levels for clinical actionability."""
    LEVEL_A = "A"   # Validated, FDA-approved companion diagnostic
    LEVEL_B = "B"   # Clinical evidence from well-powered studies
    LEVEL_C = "C"   # Case reports and small series
    LEVEL_D = "D"   # Preclinical / in-vitro data
    LEVEL_E = "E"   # Inferential / computational prediction


class TherapyCategory(str, Enum):
    """Therapeutic modality classification."""
    TARGETED = "targeted"
    IMMUNOTHERAPY = "immunotherapy"
    CHEMOTHERAPY = "chemotherapy"
    HORMONAL = "hormonal"
    COMBINATION = "combination"
    RADIOTHERAPY = "radiotherapy"
    CELL_THERAPY = "cell_therapy"
    ADC = "adc"
    BISPECIFIC = "bispecific"


class TrialPhase(str, Enum):
    """Clinical trial phase designation."""
    EARLY_PHASE_1 = "Early Phase 1"
    PHASE_1 = "Phase 1"
    PHASE_1_2 = "Phase 1/Phase 2"
    PHASE_2 = "Phase 2"
    PHASE_2_3 = "Phase 2/Phase 3"
    PHASE_3 = "Phase 3"
    PHASE_4 = "Phase 4"
    NA = "N/A"


class TrialStatus(str, Enum):
    """ClinicalTrials.gov recruitment status."""
    NOT_YET_RECRUITING = "Not yet recruiting"
    RECRUITING = "Recruiting"
    ENROLLING_BY_INVITATION = "Enrolling by invitation"
    ACTIVE_NOT_RECRUITING = "Active, not recruiting"
    SUSPENDED = "Suspended"
    TERMINATED = "Terminated"
    COMPLETED = "Completed"
    WITHDRAWN = "Withdrawn"
    UNKNOWN = "Unknown status"


class ResponseCategory(str, Enum):
    """RECIST-style tumor response classification."""
    CR = "complete_response"
    PR = "partial_response"
    SD = "stable_disease"
    PD = "progressive_disease"
    NE = "not_evaluable"


class BiomarkerType(str, Enum):
    """Biomarker functional classification."""
    PREDICTIVE = "predictive"
    PROGNOSTIC = "prognostic"
    DIAGNOSTIC = "diagnostic"
    MONITORING = "monitoring"
    RESISTANCE = "resistance"
    PHARMACODYNAMIC = "pharmacodynamic"
    SCREENING = "screening"
    THERAPEUTIC_SELECTION = "therapeutic_selection"


class PathwayName(str, Enum):
    """Major oncogenic signalling pathways."""
    MAPK = "mapk"
    PI3K_AKT_MTOR = "pi3k_akt_mtor"
    DDR = "dna_damage_repair"
    CELL_CYCLE = "cell_cycle"
    APOPTOSIS = "apoptosis"
    WNT = "wnt"
    NOTCH = "notch"
    HEDGEHOG = "hedgehog"
    JAK_STAT = "jak_stat"
    ANGIOGENESIS = "angiogenesis"
    HIPPO = "hippo"
    NF_KB = "nf_kb"
    TGF_BETA = "tgf_beta"


class GuidelineOrg(str, Enum):
    """Guideline-issuing organisations."""
    NCCN = "NCCN"
    ESMO = "ESMO"
    ASCO = "ASCO"
    WHO = "WHO"
    CAP_AMP = "CAP/AMP"
    FDA = "FDA"
    EMA = "EMA"
    AACR = "AACR"


class SourceType(str, Enum):
    """Literature source provenance."""
    PUBMED = "pubmed"
    PMC = "pmc"
    PREPRINT = "preprint"
    MANUAL = "manual"


# ═══════════════════════════════════════════════════════════════════════════
#  Domain Models
# ═══════════════════════════════════════════════════════════════════════════


class OncologyLiterature(BaseModel):
    """Chunk of oncology literature indexed for RAG retrieval."""
    id: str
    title: str
    text_chunk: str
    source_type: SourceType
    year: Optional[int] = None
    cancer_type: Optional[CancerType] = None
    gene: Optional[str] = None
    variant: Optional[str] = None
    keywords: List[str] = Field(default_factory=list)
    journal: Optional[str] = None

    def to_embedding_text(self) -> str:
        parts = [self.title, self.text_chunk]
        if self.gene:
            parts.append(f"Gene: {self.gene}")
        if self.variant:
            parts.append(f"Variant: {self.variant}")
        if self.cancer_type:
            parts.append(f"Cancer: {self.cancer_type.value}")
        if self.keywords:
            parts.append(f"Keywords: {', '.join(self.keywords)}")
        return " | ".join(parts)


class OncologyTrial(BaseModel):
    """Clinical trial record for oncology-specific matching."""
    id: str = Field(..., pattern=r"^NCT\d+$")
    title: str
    text_summary: str
    phase: TrialPhase
    status: TrialStatus
    sponsor: Optional[str] = None
    cancer_types: List[CancerType] = Field(default_factory=list)
    biomarker_criteria: List[str] = Field(default_factory=list)
    enrollment: Optional[int] = None
    start_year: Optional[int] = None
    outcome_summary: Optional[str] = None

    def to_embedding_text(self) -> str:
        parts = [self.title, self.text_summary]
        parts.append(f"Phase: {self.phase.value}")
        parts.append(f"Status: {self.status.value}")
        if self.cancer_types:
            parts.append(f"Cancers: {', '.join(c.value for c in self.cancer_types)}")
        if self.biomarker_criteria:
            parts.append(f"Biomarkers: {', '.join(self.biomarker_criteria)}")
        if self.outcome_summary:
            parts.append(f"Outcomes: {self.outcome_summary}")
        return " | ".join(parts)


class OncologyVariant(BaseModel):
    """Clinically annotated genomic variant."""
    id: str
    gene: str
    variant_name: str
    variant_type: VariantType
    cancer_type: Optional[CancerType] = None
    evidence_level: EvidenceLevel
    drugs: List[str] = Field(default_factory=list)
    civic_id: Optional[str] = None
    vrs_id: Optional[str] = None
    text_summary: str
    clinical_significance: Optional[str] = None
    allele_frequency: Optional[float] = None

    def to_embedding_text(self) -> str:
        parts = [
            f"{self.gene} {self.variant_name}",
            self.text_summary,
            f"Type: {self.variant_type.value}",
            f"Evidence: {self.evidence_level.value}",
        ]
        if self.cancer_type:
            parts.append(f"Cancer: {self.cancer_type.value}")
        if self.drugs:
            parts.append(f"Drugs: {', '.join(self.drugs)}")
        if self.clinical_significance:
            parts.append(f"Significance: {self.clinical_significance}")
        return " | ".join(parts)


class OncologyBiomarker(BaseModel):
    """Biomarker with clinical testing context."""
    id: str
    name: str
    biomarker_type: BiomarkerType
    cancer_types: List[CancerType] = Field(default_factory=list)
    predictive_value: Optional[str] = None
    testing_method: Optional[str] = None
    clinical_cutoff: Optional[str] = None
    text_summary: str
    evidence_level: EvidenceLevel

    def to_embedding_text(self) -> str:
        parts = [
            self.name,
            self.text_summary,
            f"Type: {self.biomarker_type.value}",
            f"Evidence: {self.evidence_level.value}",
        ]
        if self.cancer_types:
            parts.append(f"Cancers: {', '.join(c.value for c in self.cancer_types)}")
        if self.testing_method:
            parts.append(f"Method: {self.testing_method}")
        if self.clinical_cutoff:
            parts.append(f"Cutoff: {self.clinical_cutoff}")
        return " | ".join(parts)


class OncologyTherapy(BaseModel):
    """Therapeutic agent with mechanism and indication data."""
    id: str
    drug_name: str
    category: TherapyCategory
    targets: List[str] = Field(default_factory=list)
    approved_indications: List[str] = Field(default_factory=list)
    resistance_mechanisms: List[str] = Field(default_factory=list)
    evidence_level: EvidenceLevel
    text_summary: str
    mechanism_of_action: Optional[str] = None

    def to_embedding_text(self) -> str:
        parts = [
            self.drug_name,
            self.text_summary,
            f"Category: {self.category.value}",
            f"Evidence: {self.evidence_level.value}",
        ]
        if self.targets:
            parts.append(f"Targets: {', '.join(self.targets)}")
        if self.approved_indications:
            parts.append(f"Indications: {', '.join(self.approved_indications)}")
        if self.mechanism_of_action:
            parts.append(f"MoA: {self.mechanism_of_action}")
        return " | ".join(parts)


class OncologyPathway(BaseModel):
    """Oncogenic signalling pathway with therapeutic context."""
    id: str
    name: PathwayName
    key_genes: List[str] = Field(default_factory=list)
    therapeutic_targets: List[str] = Field(default_factory=list)
    cross_talk: List[str] = Field(default_factory=list)
    text_summary: str

    def to_embedding_text(self) -> str:
        parts = [
            f"Pathway: {self.name.value}",
            self.text_summary,
        ]
        if self.key_genes:
            parts.append(f"Genes: {', '.join(self.key_genes)}")
        if self.therapeutic_targets:
            parts.append(f"Targets: {', '.join(self.therapeutic_targets)}")
        if self.cross_talk:
            parts.append(f"Cross-talk: {', '.join(self.cross_talk)}")
        return " | ".join(parts)


class OncologyGuideline(BaseModel):
    """Clinical practice guideline recommendation."""
    id: str
    org: GuidelineOrg
    cancer_type: CancerType
    version: str
    year: int
    key_recommendations: List[str] = Field(default_factory=list)
    text_summary: str
    evidence_level: EvidenceLevel

    def to_embedding_text(self) -> str:
        parts = [
            f"{self.org.value} {self.cancer_type.value} v{self.version} ({self.year})",
            self.text_summary,
            f"Evidence: {self.evidence_level.value}",
        ]
        if self.key_recommendations:
            parts.append(f"Recommendations: {'; '.join(self.key_recommendations)}")
        return " | ".join(parts)


class ResistanceMechanism(BaseModel):
    """Documented mechanism of therapeutic resistance."""
    id: str
    primary_therapy: str
    gene: str
    mechanism: str
    bypass_pathway: Optional[str] = None
    alternative_therapies: List[str] = Field(default_factory=list)
    text_summary: str

    def to_embedding_text(self) -> str:
        parts = [
            f"Resistance to {self.primary_therapy}",
            f"Gene: {self.gene}",
            f"Mechanism: {self.mechanism}",
            self.text_summary,
        ]
        if self.bypass_pathway:
            parts.append(f"Bypass: {self.bypass_pathway}")
        if self.alternative_therapies:
            parts.append(f"Alternatives: {', '.join(self.alternative_therapies)}")
        return " | ".join(parts)


class OutcomeRecord(BaseModel):
    """Real-world or trial outcome observation."""
    id: str
    case_id: str
    therapy: str
    cancer_type: CancerType
    response: ResponseCategory
    duration_months: Optional[float] = None
    toxicities: List[str] = Field(default_factory=list)
    biomarkers_at_baseline: Dict[str, str] = Field(default_factory=dict)
    text_summary: str

    def to_embedding_text(self) -> str:
        parts = [
            f"Therapy: {self.therapy}",
            f"Cancer: {self.cancer_type.value}",
            f"Response: {self.response.value}",
            self.text_summary,
        ]
        if self.duration_months is not None:
            parts.append(f"Duration: {self.duration_months} months")
        if self.toxicities:
            parts.append(f"Toxicities: {', '.join(self.toxicities)}")
        if self.biomarkers_at_baseline:
            bm = "; ".join(f"{k}={v}" for k, v in self.biomarkers_at_baseline.items())
            parts.append(f"Biomarkers: {bm}")
        return " | ".join(parts)


class CaseSnapshot(BaseModel):
    """De-identified patient snapshot for MTB-style reasoning."""
    case_id: str = Field(..., alias="id", description="Unique case identifier (also accessible as 'id')")
    patient_id: str
    cancer_type: str = Field(..., description="Cancer type string (e.g. 'NSCLC', 'nsclc', 'breast')")
    stage: Optional[str] = None
    variants: List = Field(default_factory=list, description="Variants as list of strings or dicts")
    biomarkers: Dict[str, Any] = Field(default_factory=dict)
    prior_therapies: List[str] = Field(default_factory=list)
    text_summary: str = Field(default="", description="Free-text case summary")
    created_at: Optional[str] = None
    updated_at: Optional[str] = None

    model_config = ConfigDict(populate_by_name=True)

    def to_embedding_text(self) -> str:
        parts = [
            f"Cancer: {self.cancer_type}",
            self.text_summary,
        ]
        if self.stage:
            parts.append(f"Stage: {self.stage}")
        if self.variants:
            var_strs = [
                v if isinstance(v, str) else v.get("gene", str(v))
                for v in self.variants
            ]
            parts.append(f"Variants: {', '.join(var_strs)}")
        if self.biomarkers:
            bm = "; ".join(f"{k}={v}" for k, v in self.biomarkers.items())
            parts.append(f"Biomarkers: {bm}")
        if self.prior_therapies:
            parts.append(f"Prior therapies: {', '.join(self.prior_therapies)}")
        return " | ".join(parts)


# ═══════════════════════════════════════════════════════════════════════════
#  Search & Agent Models
# ═══════════════════════════════════════════════════════════════════════════


class SearchHit(BaseModel):
    """Single result from a Milvus collection search."""
    collection: str
    id: str
    score: float
    text: str
    metadata: Dict = Field(default_factory=dict)
    label: str = ""
    citation: str = ""
    relevance: str = ""

    @property
    def record_id(self) -> str:
        """Alias for id, used by the RAG engine for deduplication."""
        return self.id


class CrossCollectionResult(BaseModel):
    """Aggregated search results spanning multiple Milvus collections."""
    query: str
    hits: List[SearchHit] = Field(default_factory=list)
    knowledge_context: Optional[str] = None
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
    """Side-by-side evidence comparison for two entities."""
    query: str
    entity_a: str
    entity_b: str
    evidence_a: List[SearchHit] = Field(default_factory=list)
    evidence_b: List[SearchHit] = Field(default_factory=list)
    comparison_context: Optional[str] = None
    total_search_time_ms: float = 0.0

    @property
    def total_hits(self) -> int:
        return len(self.evidence_a) + len(self.evidence_b)


class MTBPacket(BaseModel):
    """Molecular Tumor Board decision-support packet."""
    case_id: str
    patient_summary: str = Field(default="", description="Patient summary text")
    patient_id: Optional[str] = None
    cancer_type: Optional[str] = None
    stage: Optional[str] = None
    variant_table: List[Dict] = Field(default_factory=list)
    evidence_table: List[Dict] = Field(default_factory=list)
    therapy_ranking: List[Dict] = Field(default_factory=list)
    trial_matches: List[Dict] = Field(default_factory=list)
    open_questions: List[str] = Field(default_factory=list)
    generated_at: Any = Field(default_factory=lambda: datetime.now(timezone.utc))
    citations: List[str] = Field(default_factory=list)


class AgentQuery(BaseModel):
    """Inbound query to the Oncology Intelligence Agent."""
    question: str
    cancer_type: Optional[CancerType] = None
    gene: Optional[str] = None
    include_genomic: bool = True

    @property
    def text(self) -> str:
        """Alias for question, used by the RAG engine."""
        return self.question


class AgentResponse(BaseModel):
    """Structured response from the Oncology Intelligence Agent."""
    question: str
    answer: str
    evidence: CrossCollectionResult
    knowledge_used: List[str] = Field(default_factory=list)
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    plan: Any = None
    report: Optional[str] = None
