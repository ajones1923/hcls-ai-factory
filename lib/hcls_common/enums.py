"""
Domain enumerations for the HCLS AI Factory precision medicine pipeline.

Provides strongly-typed enumerations used across all three pipeline stages
(Genomics, RAG/Chat, Drug Discovery) and the multi-agent constellation
(Biomarker, Oncology, Imaging, CAR-T, Autoimmune).
"""

from enum import Enum, unique


@unique
class PipelineStage(str, Enum):
    """Stages in the end-to-end Patient DNA -> Drug Candidates pipeline."""

    GENOMICS = "genomics"
    ANNOTATION = "annotation"
    RAG_INGEST = "rag_ingest"
    RAG_CHAT = "rag_chat"
    DRUG_DISCOVERY = "drug_discovery"
    CART_ANALYSIS = "cart_analysis"
    BIOMARKER_ANALYSIS = "biomarker_analysis"
    ONCOLOGY_ANALYSIS = "oncology_analysis"
    IMAGING_ANALYSIS = "imaging_analysis"
    AUTOIMMUNE_ANALYSIS = "autoimmune_analysis"
    PGX_CONTRAINDICATION = "pgx_contraindication"
    REPORTING = "reporting"

    @property
    def display_name(self) -> str:
        """Human-readable stage name."""
        return self.value.replace("_", " ").title()


@unique
class ClinicalSignificance(str, Enum):
    """ClinVar clinical significance classifications for genomic variants."""

    PATHOGENIC = "pathogenic"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    UNCERTAIN = "uncertain_significance"
    LIKELY_BENIGN = "likely_benign"
    BENIGN = "benign"

    @property
    def is_actionable(self) -> bool:
        """Whether this significance level warrants clinical follow-up."""
        return self in (
            ClinicalSignificance.PATHOGENIC,
            ClinicalSignificance.LIKELY_PATHOGENIC,
        )

    @property
    def display_name(self) -> str:
        return self.value.replace("_", " ").title()


@unique
class VariantImpact(str, Enum):
    """VEP / SnpEff variant impact categories."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"

    @property
    def rank(self) -> int:
        """Numeric rank for sorting (lower = more impactful)."""
        return {
            VariantImpact.HIGH: 0,
            VariantImpact.MODERATE: 1,
            VariantImpact.LOW: 2,
            VariantImpact.MODIFIER: 3,
        }[self]


@unique
class TherapeuticArea(str, Enum):
    """Therapeutic areas supported by the precision medicine platform."""

    NEURODEGENERATION = "neurodegeneration"
    ONCOLOGY = "oncology"
    IMMUNOLOGY = "immunology"
    CARDIOLOGY = "cardiology"
    RARE_DISEASE = "rare_disease"

    @property
    def display_name(self) -> str:
        return self.value.replace("_", " ").title()


@unique
class DrugLikenessFilter(str, Enum):
    """Drug-likeness rule sets used during molecule scoring in Stage 3."""

    LIPINSKI = "lipinski"
    VEBER = "veber"
    LEAD_LIKE = "lead_like"
    FRAGMENT = "fragment"

    @property
    def description(self) -> str:
        descs = {
            DrugLikenessFilter.LIPINSKI: "Lipinski Rule of Five (MW<=500, LogP<=5, HBD<=5, HBA<=10)",
            DrugLikenessFilter.VEBER: "Veber rules (RotBonds<=10, TPSA<=140)",
            DrugLikenessFilter.LEAD_LIKE: "Lead-like (MW 200-350, LogP -1..3.5, RotBonds<=7)",
            DrugLikenessFilter.FRAGMENT: "Fragment (MW<=300, LogP<=3, HBD<=3, HBA<=3, RotBonds<=3)",
        }
        return descs[self]


@unique
class DockingMethod(str, Enum):
    """Molecular docking methods supported or referenced by the platform."""

    DIFFDOCK = "diffdock"
    VINA = "vina"
    GLIDE = "glide"
    GOLD = "gold"


@unique
class LLMProvider(str, Enum):
    """Supported LLM providers for the RAG Chat pipeline."""

    ANTHROPIC = "anthropic"
    OPENAI = "openai"
    OLLAMA = "ollama"
    VLLM = "vllm"


@unique
class EmbeddingProvider(str, Enum):
    """Supported embedding providers."""

    LOCAL_BGE = "local_bge"
    TEI = "tei"
    OPENAI = "openai"


@unique
class NIMService(str, Enum):
    """NVIDIA NIM microservices used in the Drug Discovery pipeline."""

    MOLMIM = "molmim"
    DIFFDOCK = "diffdock"


@unique
class AlertSeverity(str, Enum):
    """Severity levels for observability alerts (Grafana/Prometheus)."""

    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    INFO = "info"

    @property
    def rank(self) -> int:
        """Numeric rank for sorting (lower = more severe)."""
        return {
            AlertSeverity.CRITICAL: 0,
            AlertSeverity.HIGH: 1,
            AlertSeverity.MEDIUM: 2,
            AlertSeverity.LOW: 3,
            AlertSeverity.INFO: 4,
        }[self]
