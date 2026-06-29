"""Pydantic data models for Pharmacogenomics Intelligence Agent.

Maps to the 15 Milvus collections (14 PGx-specific + genomic_evidence).
Follows the same dataclass/Pydantic pattern as:
  - rag-chat-pipeline/src/vcf_parser.py (VariantEvidence)
  - drug-discovery-pipeline/src/models.py (GeneratedMolecule, DockingResult)
  - cart_intelligence_agent/src/models.py (CARTLiterature, ClinicalTrial)
"""

from datetime import datetime, timezone
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ===================================================================
# ENUMS
# ===================================================================


class MetabolizerPhenotype(str, Enum):
    """CYP enzyme metabolizer status based on diplotype activity score."""
    ULTRA_RAPID = "ultra_rapid"
    RAPID = "rapid"
    NORMAL = "normal"
    INTERMEDIATE = "intermediate"
    POOR = "poor"


class TransporterFunction(str, Enum):
    """Transporter protein functional status (e.g., SLCO1B1, ABCB1)."""
    INCREASED = "increased"
    NORMAL = "normal"
    DECREASED = "decreased"
    POOR = "poor"


class HLAStatus(str, Enum):
    """HLA allele carrier status for hypersensitivity screening."""
    POSITIVE = "positive"
    NEGATIVE = "negative"


class EnzymeDeficiency(str, Enum):
    """Non-CYP enzyme deficiency status (e.g., G6PD, TPMT, DPYD)."""
    NORMAL = "normal"
    INTERMEDIATE = "intermediate"
    DEFICIENT = "deficient"


class GuidelineBody(str, Enum):
    """Pharmacogenomics guideline issuing organization."""
    CPIC = "CPIC"
    DPWG = "DPWG"
    FDA = "FDA"
    CPNDS = "CPNDS"


class CPICLevel(str, Enum):
    """CPIC level of evidence for gene-drug pairs."""
    A = "A"
    A_B = "A/B"
    B = "B"
    C = "C"
    D = "D"


class ClinicalAction(str, Enum):
    """Recommended clinical action based on PGx result."""
    STANDARD = "standard"
    DOSE_ADJUST = "dose_adjust"
    ALTERNATIVE = "alternative"
    AVOID = "avoid"
    CONTRAINDICATED = "contraindicated"


class AlertLevel(str, Enum):
    """Clinical decision support alert severity."""
    INFO = "info"
    WARNING = "warning"
    CRITICAL = "critical"


class InteractionType(str, Enum):
    """Type of pharmacogenomic drug interaction."""
    PK = "pharmacokinetic"
    PD = "pharmacodynamic"
    EFFICACY = "efficacy"
    TOXICITY = "toxicity"


class EvidenceLevel(str, Enum):
    """Strength of pharmacogenomic evidence."""
    LEVEL_1A = "1A"
    LEVEL_1B = "1B"
    LEVEL_2A = "2A"
    LEVEL_2B = "2B"
    LEVEL_3 = "3"
    LEVEL_4 = "4"


class ReactionSeverity(str, Enum):
    """Severity of adverse drug reaction."""
    LIFE_THREATENING = "life_threatening"
    SEVERE = "severe"
    MODERATE = "moderate"
    MILD = "mild"


class DrugCategory(str, Enum):
    """Therapeutic drug category with known PGx implications."""
    OPIOIDS = "opioids"
    ANTICOAGULANTS = "anticoagulants"
    ANTIDEPRESSANTS = "antidepressants"
    ANTIPSYCHOTICS = "antipsychotics"
    STATINS = "statins"
    CHEMOTHERAPY = "chemotherapy"
    ANTICONVULSANTS = "anticonvulsants"
    ANTIVIRALS = "antivirals"
    IMMUNOSUPPRESSANTS = "immunosuppressants"
    CARDIOVASCULAR = "cardiovascular"
    PPI = "proton_pump_inhibitors"
    ANTI_GOUT = "anti_gout"


class PGxWorkflowType(str, Enum):
    """Types of pharmacogenomics query workflows."""
    GENE_QUERY = "gene_query"
    DRUG_QUERY = "drug_query"
    PROFILE_QUERY = "profile_query"
    INTERACTION_QUERY = "interaction_query"
    DOSING_QUERY = "dosing_query"
    HLA_SCREEN = "hla_screen"


class InhibitorStrength(str, Enum):
    """CYP enzyme inhibitor potency classification."""
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"


# ===================================================================
# COLLECTION MODELS (map to Milvus schemas)
# ===================================================================


class GeneReference(BaseModel):
    """PGx gene and star allele reference data -- maps to pgx_gene_reference."""
    id: str = Field(..., description="Unique identifier (e.g., CYP2D6_star4)")
    gene: str = Field(..., max_length=50, description="Gene symbol (e.g., CYP2D6)")
    star_allele: str = Field(..., max_length=50, description="Star allele designation (e.g., *4)")
    defining_variants: str = Field("", max_length=1000, description="rsIDs defining this allele")
    activity_score: Optional[float] = Field(None, ge=0.0, description="CPIC activity score")
    function_status: str = Field("", max_length=100, description="e.g., no function, decreased function")
    allele_frequency_global: Optional[float] = Field(None, ge=0.0, le=1.0)
    allele_frequency_european: Optional[float] = Field(None, ge=0.0, le=1.0)
    allele_frequency_african: Optional[float] = Field(None, ge=0.0, le=1.0)
    allele_frequency_east_asian: Optional[float] = Field(None, ge=0.0, le=1.0)
    allele_frequency_south_asian: Optional[float] = Field(None, ge=0.0, le=1.0)
    allele_frequency_latino: Optional[float] = Field(None, ge=0.0, le=1.0)
    pharmvar_id: str = Field("", max_length=50, description="PharmVar accession")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    source: str = Field("", max_length=200, description="Data source (PharmVar, CPIC, etc.)")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.star_allele}"]
        if self.function_status:
            parts.append(f"Function: {self.function_status}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.defining_variants:
            parts.append(f"Variants: {self.defining_variants}")
        return " ".join(parts)


class DrugGuideline(BaseModel):
    """Clinical PGx guideline recommendation -- maps to pgx_drug_guidelines."""
    id: str = Field(..., description="Unique identifier")
    gene: str = Field(..., max_length=50)
    drug: str = Field(..., max_length=200)
    phenotype: str = Field(..., max_length=100, description="e.g., CYP2D6 poor metabolizer")
    guideline_body: GuidelineBody = GuidelineBody.CPIC
    cpic_level: CPICLevel = CPICLevel.A
    recommendation: str = Field(..., max_length=2000)
    clinical_action: ClinicalAction = ClinicalAction.STANDARD
    alert_level: AlertLevel = AlertLevel.INFO
    alternative_drugs: str = Field("", max_length=500, description="Comma-separated alternatives")
    dose_adjustment: str = Field("", max_length=500, description="Dosing modification details")
    evidence_pmids: str = Field("", max_length=500, description="Comma-separated PMIDs")
    guideline_version: str = Field("", max_length=50)
    last_updated: str = Field("", max_length=20, description="YYYY-MM-DD")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.drug} {self.phenotype}"]
        if self.recommendation:
            parts.append(self.recommendation)
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.dose_adjustment:
            parts.append(f"Dose adjustment: {self.dose_adjustment}")
        return " ".join(parts)


class DrugInteraction(BaseModel):
    """Gene-drug interaction record -- maps to pgx_drug_interactions."""
    id: str = Field(..., description="Unique identifier")
    drug: str = Field(..., max_length=200)
    gene: str = Field(..., max_length=50)
    variant_rsid: str = Field("", max_length=50, description="e.g., rs1045642")
    interaction_type: InteractionType = InteractionType.PK
    effect_description: str = Field(..., max_length=2000)
    evidence_level: EvidenceLevel = EvidenceLevel.LEVEL_2A
    clinical_significance: str = Field("", max_length=500)
    pharmgkb_id: str = Field("", max_length=50, description="PharmGKB clinical annotation ID")
    affected_phenotype: str = Field("", max_length=200)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.drug} {self.gene}"]
        if self.effect_description:
            parts.append(self.effect_description)
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.clinical_significance:
            parts.append(f"Significance: {self.clinical_significance}")
        return " ".join(parts)


class HLAHypersensitivity(BaseModel):
    """HLA-mediated drug hypersensitivity record -- maps to pgx_hla_hypersensitivity."""
    id: str = Field(..., description="Unique identifier")
    hla_allele: str = Field(..., max_length=50, description="e.g., HLA-B*57:01")
    drug: str = Field(..., max_length=200)
    reaction_type: str = Field(..., max_length=200, description="e.g., SJS/TEN, DRESS, HSR")
    risk_if_positive: str = Field(..., max_length=500)
    severity: ReactionSeverity = ReactionSeverity.SEVERE
    cpic_level: CPICLevel = CPICLevel.A
    recommendation: str = Field(..., max_length=2000)
    screening_mandatory: bool = Field(False, description="Whether pre-treatment screening is required")
    prevalence_european: Optional[float] = Field(None, ge=0.0, le=1.0)
    prevalence_african: Optional[float] = Field(None, ge=0.0, le=1.0)
    prevalence_east_asian: Optional[float] = Field(None, ge=0.0, le=1.0)
    prevalence_south_asian: Optional[float] = Field(None, ge=0.0, le=1.0)
    prevalence_latino: Optional[float] = Field(None, ge=0.0, le=1.0)
    alternative_drugs: str = Field("", max_length=500, description="Comma-separated alternatives")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.hla_allele} {self.drug} {self.reaction_type}"]
        if self.risk_if_positive:
            parts.append(self.risk_if_positive)
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.recommendation:
            parts.append(f"Recommendation: {self.recommendation}")
        return " ".join(parts)


class Phenoconversion(BaseModel):
    """Drug-induced phenoconversion record -- maps to pgx_phenoconversion."""
    id: str = Field(..., description="Unique identifier")
    affected_enzyme: str = Field(..., max_length=50, description="e.g., CYP2D6")
    precipitant_drug: str = Field(..., max_length=200, description="Drug causing the inhibition")
    interaction_type: InhibitorStrength = InhibitorStrength.STRONG
    effect_on_phenotype: str = Field(
        ..., max_length=500,
        description="e.g., converts NM to PM phenotype"
    )
    clinical_significance: str = Field("", max_length=500)
    affected_substrate_drugs: str = Field(
        "", max_length=1000,
        description="Comma-separated substrate drugs affected"
    )
    time_to_onset: str = Field("", max_length=100, description="e.g., 2-3 days")
    reversibility: str = Field("", max_length=200, description="Reversible, irreversible, time to recovery")
    evidence_level: EvidenceLevel = EvidenceLevel.LEVEL_2A
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.precipitant_drug} inhibits {self.affected_enzyme}"]
        if self.effect_on_phenotype:
            parts.append(self.effect_on_phenotype)
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.affected_substrate_drugs:
            parts.append(f"Affected drugs: {self.affected_substrate_drugs}")
        return " ".join(parts)


class DosingAlgorithm(BaseModel):
    """Genotype-guided dosing algorithm -- maps to pgx_dosing_algorithms."""
    id: str = Field(..., description="Unique identifier")
    drug: str = Field(..., max_length=200)
    genes_involved: str = Field(..., max_length=200, description="Comma-separated gene symbols")
    algorithm_name: str = Field(..., max_length=200)
    input_variables: str = Field(
        "", max_length=1000,
        description="Comma-separated input variables (genotype, age, weight, etc.)"
    )
    formula_description: str = Field("", max_length=2000, description="Algorithm or formula details")
    validation_cohort: str = Field("", max_length=500)
    accuracy_metrics: str = Field("", max_length=500, description="R-squared, MAE, etc.")
    clinical_context: str = Field("", max_length=500)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.algorithm_name} for {self.drug}"]
        if self.genes_involved:
            parts.append(f"Genes: {self.genes_involved}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.formula_description:
            parts.append(self.formula_description)
        return " ".join(parts)


class ClinicalEvidence(BaseModel):
    """Published PGx clinical study -- maps to pgx_clinical_evidence."""
    id: str = Field(..., description="Unique identifier (e.g., PMID)")
    title: str = Field(..., max_length=500)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    study_type: str = Field("", max_length=100, description="RCT, cohort, case-control, meta-analysis")
    gene: str = Field("", max_length=50)
    drug: str = Field("", max_length=200)
    phenotype: str = Field("", max_length=100)
    outcome_measure: str = Field("", max_length=200)
    outcome_value: str = Field("", max_length=200)
    sample_size: int = Field(0, ge=0)
    pmid: str = Field("", max_length=20)
    year: int = Field(0, ge=0, le=2030)
    source: str = Field("", max_length=200)

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [self.title]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.gene and self.drug:
            parts.append(f"Gene-drug: {self.gene}-{self.drug}")
        if self.outcome_measure and self.outcome_value:
            parts.append(f"Outcome: {self.outcome_measure} {self.outcome_value}")
        return " ".join(parts)


class PopulationData(BaseModel):
    """Population allele frequency data -- maps to pgx_population_data."""
    id: str = Field(..., description="Unique identifier")
    gene: str = Field(..., max_length=50)
    star_allele: str = Field(..., max_length=50)
    population: str = Field(..., max_length=100, description="e.g., East Asian, European, African")
    allele_frequency: float = Field(..., ge=0.0, le=1.0)
    sample_size: int = Field(0, ge=0)
    source_study: str = Field("", max_length=200)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.star_allele} in {self.population}"]
        parts.append(f"Frequency: {self.allele_frequency}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.source_study:
            parts.append(f"Source: {self.source_study}")
        return " ".join(parts)


class PGxClinicalTrial(BaseModel):
    """PGx-related clinical trial -- maps to pgx_clinical_trials."""
    id: str = Field(..., description="Unique identifier")
    title: str = Field(..., max_length=500)
    text_summary: str = Field(..., max_length=3000)
    nct_id: str = Field("", max_length=20, description="ClinicalTrials.gov NCT number")
    phase: str = Field("", max_length=50)
    status: str = Field("", max_length=50)
    gene: str = Field("", max_length=50)
    drug: str = Field("", max_length=200)
    enrollment: int = Field(0, ge=0)
    start_year: int = Field(0, ge=0, le=2030)
    outcome_summary: str = Field("", max_length=2000)

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [self.title, self.text_summary]
        if self.gene and self.drug:
            parts.append(f"Gene-drug: {self.gene}-{self.drug}")
        if self.outcome_summary:
            parts.append(f"Outcome: {self.outcome_summary}")
        return " ".join(parts)


class FDALabel(BaseModel):
    """FDA drug label PGx information -- maps to pgx_fda_labels."""
    id: str = Field(..., description="Unique identifier")
    drug: str = Field(..., max_length=200)
    gene: str = Field(..., max_length=50)
    labeling_section: str = Field(
        "", max_length=200,
        description="e.g., Dosage and Administration, Warnings, Clinical Pharmacology"
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    label_type: str = Field(
        "", max_length=100,
        description="e.g., actionable PGx, informative PGx, testing required"
    )
    requirement_level: str = Field(
        "", max_length=100,
        description="e.g., required, recommended, informational"
    )
    last_updated: str = Field("", max_length=20, description="YYYY-MM-DD")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.drug} {self.gene} FDA label"]
        if self.labeling_section:
            parts.append(f"Section: {self.labeling_section}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.label_type:
            parts.append(f"Type: {self.label_type}")
        return " ".join(parts)


class DrugAlternative(BaseModel):
    """PGx-informed drug alternative recommendation -- maps to pgx_drug_alternatives."""
    id: str = Field(..., description="Unique identifier")
    primary_drug: str = Field(..., max_length=200)
    gene: str = Field(..., max_length=50)
    phenotype: str = Field(..., max_length=100)
    alternative_drug: str = Field(..., max_length=200)
    rationale: str = Field(..., max_length=1000)
    evidence_level: EvidenceLevel = EvidenceLevel.LEVEL_2A
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [
            f"Alternative to {self.primary_drug}: {self.alternative_drug}",
            f"for {self.gene} {self.phenotype}",
        ]
        if self.rationale:
            parts.append(self.rationale)
        if self.text_chunk:
            parts.append(self.text_chunk)
        return " ".join(parts)


class PatientProfile(BaseModel):
    """Patient PGx genotype profile -- maps to pgx_patient_profiles."""
    id: str = Field(..., description="Unique identifier")
    patient_id: str = Field(..., max_length=100, description="De-identified patient ID")
    gene: str = Field(..., max_length=50)
    diplotype: str = Field(..., max_length=100, description="e.g., *1/*4")
    phenotype: str = Field("", max_length=100, description="e.g., intermediate metabolizer")
    activity_score: Optional[float] = Field(None, ge=0.0, description="CPIC activity score")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.diplotype}"]
        if self.phenotype:
            parts.append(f"Phenotype: {self.phenotype}")
        if self.activity_score is not None:
            parts.append(f"Activity score: {self.activity_score}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        return " ".join(parts)


class ImplementationProtocol(BaseModel):
    """PGx clinical implementation protocol -- maps to pgx_implementation."""
    id: str = Field(..., description="Unique identifier")
    institution: str = Field(..., max_length=200)
    title: str = Field(..., max_length=500)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    workflow_type: str = Field("", max_length=100, description="e.g., preemptive, reactive, hybrid")
    outcomes: str = Field("", max_length=1000, description="Implementation outcomes summary")
    year: int = Field(0, ge=0, le=2030)

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.institution}: {self.title}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.workflow_type:
            parts.append(f"Workflow: {self.workflow_type}")
        if self.outcomes:
            parts.append(f"Outcomes: {self.outcomes}")
        return " ".join(parts)


class EducationMaterial(BaseModel):
    """PGx educational content -- maps to pgx_education."""
    id: str = Field(..., description="Unique identifier")
    title: str = Field(..., max_length=500)
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    audience: str = Field("", max_length=100, description="e.g., clinician, pharmacist, patient")
    topic: str = Field("", max_length=200)
    complexity_level: str = Field("", max_length=50, description="e.g., basic, intermediate, advanced")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [self.title]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.audience:
            parts.append(f"Audience: {self.audience}")
        if self.topic:
            parts.append(f"Topic: {self.topic}")
        return " ".join(parts)


# ===================================================================
# SEARCH RESULT MODELS
# ===================================================================


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


# ===================================================================
# AGENT MODELS
# ===================================================================


class PGxAlert(BaseModel):
    """Clinical decision support alert for a PGx finding."""
    alert_level: AlertLevel = AlertLevel.INFO
    gene: str = Field("", max_length=50)
    drug: str = Field("", max_length=200)
    phenotype: str = Field("", max_length=100)
    recommendation: str = Field("", max_length=2000)
    evidence_pmids: str = Field("", max_length=500, description="Comma-separated PMIDs")


class PGxQuery(BaseModel):
    """Input to the Pharmacogenomics Intelligence Agent."""
    question: str
    patient_id: Optional[str] = None
    medication_list: Optional[List[str]] = Field(default=None, description="Current medications")
    gene_filter: Optional[str] = Field(default=None, description="Filter by gene symbol")
    drug_filter: Optional[str] = Field(default=None, description="Filter by drug name")


class AgentQuery(BaseModel):
    """Alias for RAG engine queries with gene/drug fields.

    Used by PGxRAGEngine for retrieve() and comparative analysis calls.
    Compatible with PGxQuery but uses ``gene`` / ``drug`` field names
    to match the RAG engine's filter-building logic.
    """
    question: str
    patient_id: Optional[str] = None
    medication_list: Optional[List[str]] = Field(default=None)
    gene: Optional[str] = Field(default=None, description="Gene symbol filter")
    drug: Optional[str] = Field(default=None, description="Drug name filter")


class PGxResponse(BaseModel):
    """Output from the Pharmacogenomics Intelligence Agent."""
    question: str
    answer: str
    evidence: CrossCollectionResult
    knowledge_used: List[str] = Field(default_factory=list)
    alerts: List[PGxAlert] = Field(default_factory=list)
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )
