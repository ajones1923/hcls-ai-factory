"""Pydantic data models for Biomarker Intelligence Agent.

Maps to the 10 biomarker Milvus collections + patient analysis models.
Follows the same Pydantic pattern as:
  - cart_intelligence_agent/src/models.py
  - rag-chat-pipeline/src/vcf_parser.py (VariantEvidence)

Author: Adam Jones
Date: March 2026
"""

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, model_validator

# ═══════════════════════════════════════════════════════════════════════
# ENUMS
# ═══════════════════════════════════════════════════════════════════════


class RiskLevel(str, Enum):
    """Patient risk classification levels."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    NORMAL = "normal"


class ClockType(str, Enum):
    """Epigenetic clock algorithms for biological age estimation."""
    PHENOAGE = "PhenoAge"
    GRIMAGE = "GrimAge"


class DiseaseCategory(str, Enum):
    """Major disease categories tracked by biomarker trajectories."""
    DIABETES = "diabetes"
    CARDIOVASCULAR = "cardiovascular"
    LIVER = "liver"
    THYROID = "thyroid"
    IRON = "iron"
    NUTRITIONAL = "nutritional"
    KIDNEY = "kidney"
    BONE_HEALTH = "bone_health"
    COGNITIVE = "cognitive"


class MetabolizerPhenotype(str, Enum):
    """CYP enzyme metabolizer phenotype classifications (CPIC)."""
    ULTRA_RAPID = "ultra_rapid"
    NORMAL = "normal"
    INTERMEDIATE = "intermediate"
    POOR = "poor"


class CPICLevel(str, Enum):
    """CPIC evidence levels for pharmacogenomic guidelines."""
    LEVEL_1A = "1A"
    LEVEL_1B = "1B"
    LEVEL_2A = "2A"
    LEVEL_2B = "2B"
    LEVEL_3 = "3"


class Zygosity(str, Enum):
    """Genotype zygosity classifications."""
    HOMOZYGOUS_REF = "homozygous_ref"
    HETEROZYGOUS = "heterozygous"
    HOMOZYGOUS_ALT = "homozygous_alt"


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MODELS (map to Milvus schemas)
# ═══════════════════════════════════════════════════════════════════════


class BiomarkerReference(BaseModel):
    """Reference biomarker definition — maps to biomarker_reference collection."""
    id: str = Field(..., max_length=100, description="Unique biomarker identifier")
    name: str = Field(..., max_length=100, description="Biomarker display name")
    unit: str = Field("", max_length=20, description="Measurement unit (e.g., mg/dL)")
    category: str = Field(
        "", max_length=30,
        description="CBC, CMP, Lipids, Thyroid, Inflammation, Nutrients",
    )
    ref_range_min: float = Field(0.0, description="Standard reference range lower bound")
    ref_range_max: float = Field(0.0, description="Standard reference range upper bound")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")
    clinical_significance: str = Field(
        "", max_length=2000,
        description="Clinical interpretation and significance",
    )
    epigenetic_clock: Optional[str] = Field(
        None, max_length=50,
        description="PhenoAge or GrimAge coefficient if applicable",
    )
    genetic_modifiers: str = Field(
        "", max_length=500,
        description="Comma-separated genes that modify this biomarker",
    )

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.name} ({self.unit})" if self.unit else self.name]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.clinical_significance:
            parts.append(f"Significance: {self.clinical_significance}")
        if self.category:
            parts.append(f"Category: {self.category}")
        if self.genetic_modifiers:
            parts.append(f"Genetic modifiers: {self.genetic_modifiers}")
        return " ".join(parts)

    @model_validator(mode="after")
    def _check_ref_range(self) -> "BiomarkerReference":
        if self.ref_range_max > 0 and self.ref_range_min > self.ref_range_max:
            raise ValueError(
                f"ref_range_min ({self.ref_range_min}) > ref_range_max ({self.ref_range_max})"
            )
        return self


class GeneticVariant(BaseModel):
    """Genetic variant affecting biomarker levels — maps to biomarker_genetic_variants."""
    id: str = Field(..., max_length=100, description="Unique variant identifier")
    gene: str = Field(..., max_length=50, description="Gene symbol (e.g., MTHFR)")
    rs_id: str = Field(..., max_length=20, description="dbSNP rsID (e.g., rs1801133)")
    risk_allele: str = Field("", max_length=5, description="Risk allele (e.g., T)")
    protective_allele: str = Field("", max_length=5, description="Protective allele (e.g., C)")
    effect_size: str = Field(
        "", max_length=100,
        description="Effect size description (e.g., 20-30% reduction in enzyme activity)",
    )
    mechanism: str = Field(
        "", max_length=2000,
        description="Molecular mechanism of the variant effect",
    )
    disease_associations: str = Field(
        "", max_length=1000,
        description="Comma-separated disease associations",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.rs_id}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.mechanism:
            parts.append(f"Mechanism: {self.mechanism}")
        if self.disease_associations:
            parts.append(f"Diseases: {self.disease_associations}")
        return " ".join(parts)


class PGxRule(BaseModel):
    """Pharmacogenomic dosing rule — maps to biomarker_pgx_rules collection."""
    id: str = Field(..., max_length=100, description="Unique PGx rule identifier")
    gene: str = Field(..., max_length=50, description="Pharmacogene (e.g., CYP2D6)")
    star_alleles: str = Field(
        ..., max_length=100,
        description="Star allele combination (e.g., *1/*2)",
    )
    drug: str = Field(..., max_length=100, description="Drug name")
    phenotype: MetabolizerPhenotype = MetabolizerPhenotype.NORMAL
    cpic_level: CPICLevel = CPICLevel.LEVEL_1A
    recommendation: str = Field(
        "", max_length=2000,
        description="CPIC dosing recommendation",
    )
    evidence_url: str = Field(
        "", max_length=500,
        description="URL to CPIC guideline or PharmGKB entry",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.gene} {self.star_alleles} — {self.drug}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.recommendation:
            parts.append(f"Recommendation: {self.recommendation}")
        parts.append(f"Phenotype: {self.phenotype.value}")
        parts.append(f"CPIC Level: {self.cpic_level.value}")
        return " ".join(parts)


class DiseaseTrajectory(BaseModel):
    """Disease progression trajectory — maps to biomarker_disease_trajectories."""
    id: str = Field(..., max_length=100, description="Unique trajectory identifier")
    disease: DiseaseCategory = DiseaseCategory.DIABETES
    stage: str = Field(
        "", max_length=30,
        description="Disease stage (e.g., pre-diabetes, early, advanced)",
    )
    biomarker_pattern: str = Field(
        "", max_length=2000,
        description="JSON string of biomarker thresholds defining this stage",
    )
    years_to_diagnosis: float = Field(
        0.0, description="Estimated years from this stage to clinical diagnosis",
    )
    intervention_window: str = Field(
        "", max_length=500,
        description="Description of the intervention opportunity",
    )
    risk_reduction_pct: float = Field(
        0.0, description="Potential risk reduction percentage with intervention",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.disease.value} — stage: {self.stage}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.intervention_window:
            parts.append(f"Intervention: {self.intervention_window}")
        if self.risk_reduction_pct:
            parts.append(f"Risk reduction: {self.risk_reduction_pct}%")
        return " ".join(parts)


class ClinicalEvidence(BaseModel):
    """Published clinical evidence — maps to biomarker_clinical_evidence."""
    id: str = Field(..., max_length=100, description="Unique evidence identifier")
    pmid: str = Field("", max_length=20, description="PubMed ID")
    title: str = Field(..., max_length=500, description="Publication title")
    finding: str = Field(
        "", max_length=3000,
        description="Key finding from the publication",
    )
    year: int = Field(0, ge=0, le=2030, description="Publication year")
    disease_area: str = Field("", max_length=100, description="Disease area or specialty")
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [self.title]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.finding:
            parts.append(f"Finding: {self.finding}")
        if self.disease_area:
            parts.append(f"Area: {self.disease_area}")
        return " ".join(parts)


class NutritionGuideline(BaseModel):
    """Genotype-aware nutrition guideline — maps to biomarker_nutrition."""
    id: str = Field(..., max_length=100, description="Unique guideline identifier")
    nutrient: str = Field(..., max_length=100, description="Nutrient name (e.g., Folate)")
    genetic_context: str = Field(
        "", max_length=200,
        description="Genetic context (e.g., MTHFR C677T heterozygous)",
    )
    recommended_form: str = Field(
        "", max_length=200,
        description="Recommended supplement form (e.g., methylfolate)",
    )
    dose_range: str = Field(
        "", max_length=100,
        description="Dosing range (e.g., 400-800 mcg/day)",
    )
    evidence_summary: str = Field(
        "", max_length=2000,
        description="Summary of evidence supporting this guideline",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.nutrient}"]
        if self.genetic_context:
            parts.append(f"for {self.genetic_context}")
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.recommended_form:
            parts.append(f"Form: {self.recommended_form}")
        if self.dose_range:
            parts.append(f"Dose: {self.dose_range}")
        return " ".join(parts)


class DrugInteraction(BaseModel):
    """Gene-drug interaction — maps to biomarker_drug_interactions."""
    id: str = Field(..., max_length=100, description="Unique interaction identifier")
    drug: str = Field(..., max_length=100, description="Drug name")
    gene: str = Field(..., max_length=50, description="Gene involved (e.g., CYP2C19)")
    interaction_type: str = Field(
        "", max_length=50,
        description="Interaction type (e.g., substrate, inhibitor, inducer)",
    )
    severity: str = Field(
        "", max_length=20,
        description="Severity (e.g., major, moderate, minor)",
    )
    alternative: str = Field(
        "", max_length=200,
        description="Alternative drug recommendation",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.drug} — {self.gene}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.interaction_type:
            parts.append(f"Type: {self.interaction_type}")
        if self.severity:
            parts.append(f"Severity: {self.severity}")
        if self.alternative:
            parts.append(f"Alternative: {self.alternative}")
        return " ".join(parts)


class AgingMarker(BaseModel):
    """Epigenetic aging marker — maps to biomarker_aging_markers."""
    id: str = Field(..., max_length=100, description="Unique marker identifier")
    marker_name: str = Field(..., max_length=100, description="Marker name (e.g., Albumin)")
    clock_type: ClockType = ClockType.PHENOAGE
    coefficient: float = Field(
        0.0, description="Coefficient weight in the aging clock algorithm",
    )
    unit: str = Field("", max_length=30, description="Measurement unit")
    interpretation: str = Field(
        "", max_length=2000,
        description="Clinical interpretation of this marker in aging context",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.marker_name} ({self.clock_type.value} clock)"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.interpretation:
            parts.append(f"Interpretation: {self.interpretation}")
        if self.coefficient:
            parts.append(f"Coefficient: {self.coefficient}")
        return " ".join(parts)


class GenotypeAdjustment(BaseModel):
    """Genotype-based reference range adjustment — maps to biomarker_genotype_adjustments."""
    id: str = Field(..., max_length=100, description="Unique adjustment identifier")
    biomarker: str = Field(..., max_length=100, description="Biomarker name")
    gene: str = Field(..., max_length=50, description="Gene symbol")
    rs_id: str = Field("", max_length=20, description="dbSNP rsID")
    genotype_ref: str = Field(
        "", max_length=10,
        description="Reference genotype (e.g., CC)",
    )
    genotype_het: str = Field(
        "", max_length=10,
        description="Heterozygous genotype (e.g., CT)",
    )
    genotype_hom: str = Field(
        "", max_length=10,
        description="Homozygous alternate genotype (e.g., TT)",
    )
    adjusted_min: float = Field(
        0.0, description="Adjusted reference range lower bound",
    )
    adjusted_max: float = Field(
        0.0, description="Adjusted reference range upper bound",
    )
    rationale: str = Field(
        "", max_length=2000,
        description="Rationale for the genotype-based adjustment",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"{self.biomarker} adjustment for {self.gene} {self.rs_id}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.rationale:
            parts.append(f"Rationale: {self.rationale}")
        return " ".join(parts)


class MonitoringProtocol(BaseModel):
    """Condition-specific monitoring protocol — maps to biomarker_monitoring."""
    id: str = Field(..., max_length=100, description="Unique protocol identifier")
    condition: str = Field(..., max_length=100, description="Medical condition")
    biomarker_panel: str = Field(
        "", max_length=500,
        description="Comma-separated biomarker panel",
    )
    frequency: str = Field(
        "", max_length=50,
        description="Monitoring frequency (e.g., every 3 months)",
    )
    trigger_values: str = Field(
        "", max_length=1000,
        description="Threshold values that trigger clinical action",
    )
    text_chunk: str = Field(..., max_length=3000, description="Text chunk for embedding")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small embedding."""
        parts = [f"Monitoring for {self.condition}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.biomarker_panel:
            parts.append(f"Panel: {self.biomarker_panel}")
        if self.frequency:
            parts.append(f"Frequency: {self.frequency}")
        if self.trigger_values:
            parts.append(f"Triggers: {self.trigger_values}")
        return " ".join(parts)


class CriticalValue(BaseModel):
    """Critical biomarker value threshold — maps to biomarker_critical_values collection."""
    id: str = Field(..., max_length=100)
    biomarker: str = Field(..., max_length=100)
    loinc_code: str = Field("", max_length=20)
    critical_low: float = Field(0.0)
    critical_high: float = Field(0.0)
    severity: str = Field("", max_length=20, description="critical/urgent/warning")
    escalation_target: str = Field("", max_length=100)
    clinical_action: str = Field("", max_length=2000)
    cross_checks: str = Field("", max_length=500, description="Comma-separated cross-check biomarkers")
    text_chunk: str = Field(..., max_length=3000)

    def to_embedding_text(self) -> str:
        parts = [f"Critical value: {self.biomarker} ({self.severity})"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.clinical_action:
            parts.append(f"Action: {self.clinical_action}")
        return " ".join(parts)


class DiscordanceRule(BaseModel):
    """Cross-biomarker discordance pattern — maps to biomarker_discordance_rules collection."""
    id: str = Field(..., max_length=100)
    name: str = Field(..., max_length=200)
    biomarker_a: str = Field(..., max_length=100)
    biomarker_b: str = Field(..., max_length=100)
    condition: str = Field("", max_length=500)
    expected_relationship: str = Field("", max_length=1000)
    differential_diagnosis: str = Field("", max_length=2000, description="Comma-separated differentials")
    agent_handoff: str = Field("", max_length=500, description="Comma-separated agent names")
    priority: str = Field("", max_length=10, description="high/medium/low")
    text_chunk: str = Field(..., max_length=3000)

    def to_embedding_text(self) -> str:
        parts = [f"Discordance: {self.name} — {self.biomarker_a} vs {self.biomarker_b}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.condition:
            parts.append(f"Pattern: {self.condition}")
        return " ".join(parts)


class AJCarrierScreeningEntry(BaseModel):
    """AJ carrier screening gene — maps to biomarker_aj_carrier_screening collection."""
    id: str = Field(..., max_length=100)
    gene: str = Field(..., max_length=50)
    disease: str = Field(..., max_length=200)
    inheritance: str = Field("", max_length=50)
    aj_carrier_frequency: str = Field("", max_length=30)
    general_carrier_frequency: str = Field("", max_length=30)
    common_mutations: str = Field("", max_length=500, description="Comma-separated mutations")
    loinc_code: str = Field("", max_length=20)
    method: str = Field("", max_length=500)
    clinical_significance: str = Field("", max_length=2000)
    reproductive_implications: str = Field("", max_length=2000)
    compound_risks: str = Field("", max_length=1000, description="Comma-separated compound risks")
    text_chunk: str = Field(..., max_length=3000)

    def to_embedding_text(self) -> str:
        parts = [f"AJ Carrier Screening: {self.gene} — {self.disease}"]
        if self.text_chunk:
            parts.append(self.text_chunk)
        if self.clinical_significance:
            parts.append(f"Significance: {self.clinical_significance}")
        return " ".join(parts)


# ═══════════════════════════════════════════════════════════════════════
# PATIENT & ANALYSIS MODELS
# ═══════════════════════════════════════════════════════════════════════


class PatientProfile(BaseModel):
    """Patient input profile with biomarker values and genotype data."""
    patient_id: str = Field(..., description="Unique patient identifier")
    age: int = Field(..., ge=0, le=150, description="Patient age in years")
    sex: str = Field(..., max_length=10, description="Patient sex (M/F)")
    biomarkers: Dict[str, float] = Field(
        default_factory=dict,
        description="Biomarker name -> measured value (e.g., {'HbA1c': 5.7})",
    )
    genotypes: Dict[str, str] = Field(
        default_factory=dict,
        description="rsID -> genotype mapping (e.g., {'rs1801133': 'CT'})",
    )
    star_alleles: Dict[str, str] = Field(
        default_factory=dict,
        description="Gene -> star allele mapping (e.g., {'CYP2D6': '*1/*2'})",
    )
    ancestry: Optional[str] = Field(
        None,
        description="Self-reported ancestry for population-specific adjustments (e.g., 'european', 'african', 'east_asian', 'south_asian', 'hispanic')",
        pattern="^(european|african|east_asian|south_asian|hispanic|mixed|other)$",
    )


class BiologicalAgeResult(BaseModel):
    """Result of biological age estimation using epigenetic clock biomarkers."""
    chronological_age: int = Field(..., description="Calendar age")
    biological_age: float = Field(..., description="Estimated biological age")
    age_acceleration: float = Field(
        ..., description="biological_age - chronological_age (positive = accelerated)",
    )
    phenoage_score: float = Field(..., description="PhenoAge composite score")
    grimage_score: Optional[float] = Field(
        None, description="GrimAge composite score (if available)",
    )
    grimage_data: Optional[Dict[str, Any]] = Field(
        None, description="Full GrimAge surrogate result dict (CI, confidence, validation)",
    )
    mortality_risk: float = Field(
        0.0, ge=0.0, le=1.0,
        description="Relative mortality risk multiplier",
    )
    aging_drivers: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="List of dicts with biomarker, value, contribution, direction",
    )
    confidence_interval: Optional[Dict[str, Any]] = Field(
        None,
        description="95% CI with lower, upper, confidence_level, standard_error, note",
    )
    risk_confidence: Optional[str] = Field(
        None,
        description="Confidence qualifier for risk classification: high, moderate, or low",
    )

    @model_validator(mode="after")
    def _check_age_acceleration(self) -> "BiologicalAgeResult":
        expected = round(self.biological_age - self.chronological_age, 1)
        if abs(self.age_acceleration - expected) > 0.2:
            import logging
            logging.getLogger(__name__).warning(
                f"age_acceleration ({self.age_acceleration}) doesn't match "
                f"biological_age - chronological_age ({expected})"
            )
        return self


class DiseaseTrajectoryResult(BaseModel):
    """Result of disease trajectory analysis for a patient."""
    disease: DiseaseCategory
    risk_level: RiskLevel = RiskLevel.NORMAL
    current_markers: Dict[str, float] = Field(
        default_factory=dict,
        description="Relevant biomarker name -> current value",
    )
    genetic_risk_factors: List[str] = Field(
        default_factory=list,
        description="Genetic variants contributing to disease risk",
    )
    years_to_onset_estimate: Optional[float] = Field(
        None, description="Estimated years to clinical onset (if applicable)",
    )
    intervention_recommendations: List[str] = Field(
        default_factory=list,
        description="Recommended interventions to reduce risk",
    )


class PGxResult(BaseModel):
    """Pharmacogenomic result for a single gene."""
    gene: str = Field(..., description="Pharmacogene (e.g., CYP2D6)")
    star_alleles: str = Field(..., description="Star allele result (e.g., *1/*2)")
    phenotype: MetabolizerPhenotype = MetabolizerPhenotype.NORMAL
    drugs_affected: List[Dict[str, Any]] = Field(
        default_factory=list,
        description="List of dicts with drug, recommendation, cpic_level",
    )


class GenotypeAdjustmentResult(BaseModel):
    """Result of genotype-based reference range adjustment for a biomarker."""
    biomarker: str = Field(..., description="Biomarker name")
    standard_range: str = Field(
        ..., description="Standard reference range (e.g., '4.0-5.6%')",
    )
    adjusted_range: str = Field(
        ..., description="Genotype-adjusted range (e.g., '4.2-5.8%')",
    )
    genotype: str = Field(..., description="Patient genotype at this locus")
    gene: str = Field(..., description="Gene symbol")
    rationale: str = Field(
        "", description="Explanation for the adjustment",
    )


# ═══════════════════════════════════════════════════════════════════════
# SEARCH RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════


class SearchHit(BaseModel):
    """A single search result from any collection."""
    collection: str
    id: str
    score: float = Field(..., ge=0.0)
    text: str
    metadata: Dict[str, Any] = Field(default_factory=dict)


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


# ═══════════════════════════════════════════════════════════════════════
# AGENT MODELS
# ═══════════════════════════════════════════════════════════════════════


class AgentQuery(BaseModel):
    """Input to the Biomarker Intelligence Agent."""
    question: str = Field(..., max_length=10000, description="Natural language biomarker question")
    patient_profile: Optional[PatientProfile] = None
    include_genomic: bool = True  # Also search genomic_evidence collection


class AgentResponse(BaseModel):
    """Output from the Biomarker Intelligence Agent."""
    question: str
    answer: str
    evidence: CrossCollectionResult
    biological_age: Optional[BiologicalAgeResult] = None
    disease_trajectories: Optional[List[DiseaseTrajectoryResult]] = None
    pgx_results: Optional[List[PGxResult]] = None
    genotype_adjustments: Optional[List[GenotypeAdjustmentResult]] = None
    critical_alerts: List[str] = Field(
        default_factory=list,
        description="Critical findings requiring immediate attention",
    )
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )


class AnalysisResult(BaseModel):
    """Comprehensive analysis result combining all sub-analyses."""
    patient_profile: PatientProfile
    biological_age: BiologicalAgeResult
    disease_trajectories: List[DiseaseTrajectoryResult] = Field(default_factory=list)
    pgx_results: List[PGxResult] = Field(default_factory=list)
    genotype_adjustments: List[GenotypeAdjustmentResult] = Field(default_factory=list)
    critical_alerts: List[str] = Field(
        default_factory=list,
        description="Critical findings requiring immediate attention",
    )
    timestamp: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )


# =====================================================================
# Longitudinal Tracking Models
# =====================================================================

class BiomarkerPanel(BaseModel):
    """A single time-point biomarker panel for longitudinal tracking."""
    date: str = Field(..., description="ISO-8601 date of the panel (e.g., '2025-06-15')")
    biomarkers: Dict[str, float] = Field(..., description="Biomarker name -> measured value")
    biological_age: Optional[float] = Field(None, description="Computed biological age at this time point")
    notes: Optional[str] = None


class PatientHistory(BaseModel):
    """Longitudinal patient data across multiple visits.

    Enables biological age trajectory analysis: is the patient aging
    faster or slower over time? This turns a snapshot into a story.
    """
    patient_id: str
    panels: List[BiomarkerPanel] = Field(default_factory=list, description="Time-ordered biomarker panels")

    @property
    def panel_count(self) -> int:
        return len(self.panels)

    @property
    def date_range(self) -> tuple:
        """Return (earliest_date, latest_date) or (None, None) if empty."""
        if not self.panels:
            return (None, None)
        dates = sorted(p.date for p in self.panels)
        return (dates[0], dates[-1])

    def biological_age_trajectory(self) -> List[Dict[str, Any]]:
        """Extract biological age trajectory for plotting."""
        return [
            {"date": p.date, "biological_age": p.biological_age}
            for p in self.panels
            if p.biological_age is not None
        ]

    def age_acceleration_trend(self) -> Optional[str]:
        """Determine if age acceleration is improving or worsening.

        Returns 'improving', 'stable', or 'worsening', or None if insufficient data.
        """
        trajectory = self.biological_age_trajectory()
        if len(trajectory) < 2:
            return None

        recent = trajectory[-1]["biological_age"]
        previous = trajectory[-2]["biological_age"]
        diff = recent - previous

        if diff < -0.5:
            return "improving"
        elif diff > 0.5:
            return "worsening"
        else:
            return "stable"


class WearableData(BaseModel):
    """Wearable device data for correlation with biomarker analysis.

    Captures physiological signals from consumer/clinical wearables
    that correlate with biological aging and disease trajectories.
    Schema designed for Apple Watch, Fitbit, Garmin, Oura Ring data.
    """
    device_type: Optional[str] = Field(None, description="Device manufacturer/model")
    measurement_date: str = Field(..., description="ISO-8601 date")

    # Heart rate metrics
    resting_heart_rate: Optional[float] = Field(None, ge=20, le=200, description="Resting HR in bpm")
    heart_rate_variability: Optional[float] = Field(None, ge=0, le=500, description="HRV (RMSSD) in ms")
    max_heart_rate: Optional[float] = Field(None, ge=40, le=250, description="Max HR in bpm")

    # Blood oxygen
    spo2_average: Optional[float] = Field(None, ge=70, le=100, description="Average SpO2 %")
    spo2_minimum: Optional[float] = Field(None, ge=50, le=100, description="Minimum SpO2 %")

    # Sleep
    sleep_duration_hours: Optional[float] = Field(None, ge=0, le=24, description="Total sleep duration")
    deep_sleep_pct: Optional[float] = Field(None, ge=0, le=100, description="Deep sleep percentage")
    rem_sleep_pct: Optional[float] = Field(None, ge=0, le=100, description="REM sleep percentage")
    sleep_score: Optional[float] = Field(None, ge=0, le=100, description="Composite sleep score")

    # Activity
    steps: Optional[int] = Field(None, ge=0, description="Daily step count")
    active_calories: Optional[float] = Field(None, ge=0, description="Active calories burned")
    vo2_max: Optional[float] = Field(None, ge=10, le=90, description="Estimated VO2 max (ml/kg/min)")

    # Stress / recovery
    stress_score: Optional[float] = Field(None, ge=0, le=100, description="Stress score (0=calm, 100=high)")
    recovery_score: Optional[float] = Field(None, ge=0, le=100, description="Recovery readiness score")
    body_temperature_delta: Optional[float] = Field(None, ge=-3, le=3, description="Temp deviation from baseline (°C)")
