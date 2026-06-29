"""Pydantic data models for the Rare Disease Diagnostic Agent.

Comprehensive enums and models for a rare disease RAG-based diagnostic
support system covering phenotype-driven diagnosis, WES/WGS interpretation,
metabolic screening, dysmorphology assessment, neurogenetic evaluation,
cardiac genetics, connective tissue disorders, inborn errors of metabolism,
gene therapy eligibility, and undiagnosed disease workup.

Follows the same dataclass/Pydantic pattern as:
  - clinical_trial_intelligence_agent/src/models.py
  - cardiology_intelligence_agent/src/models.py

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ===================================================================
# ENUMS
# ===================================================================


class DiagnosticWorkflowType(str, Enum):
    """Types of rare disease diagnostic query workflows."""
    PHENOTYPE_DRIVEN = "phenotype_driven"
    VARIANT_INTERPRETATION = "variant_interpretation"
    DIFFERENTIAL_DIAGNOSIS = "differential_diagnosis"
    GENE_THERAPY_ELIGIBILITY = "gene_therapy_eligibility"
    NEWBORN_SCREENING = "newborn_screening"
    METABOLIC_WORKUP = "metabolic_workup"
    CARRIER_SCREENING = "carrier_screening"
    PRENATAL_DIAGNOSIS = "prenatal_diagnosis"
    NATURAL_HISTORY = "natural_history"
    THERAPY_SELECTION = "therapy_selection"
    CLINICAL_TRIAL_MATCHING = "clinical_trial_matching"
    GENETIC_COUNSELING = "genetic_counseling"
    # Legacy aliases (kept for backward compatibility with models.py consumers)
    WES_WGS_INTERPRETATION = "wes_wgs_interpretation"
    METABOLIC_SCREENING = "metabolic_screening"
    DYSMORPHOLOGY = "dysmorphology"
    NEUROGENETIC = "neurogenetic"
    CARDIAC_GENETICS = "cardiac_genetics"
    CONNECTIVE_TISSUE = "connective_tissue"
    INBORN_ERRORS = "inborn_errors"
    UNDIAGNOSED_DISEASE = "undiagnosed_disease"
    GENERAL = "general"


class InheritancePattern(str, Enum):
    """Mendelian and non-Mendelian inheritance patterns."""
    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    X_LINKED_DOMINANT = "x_linked_dominant"
    X_LINKED_RECESSIVE = "x_linked_recessive"
    MITOCHONDRIAL = "mitochondrial"
    MULTIFACTORIAL = "multifactorial"
    DE_NOVO = "de_novo"


class ACMGClassification(str, Enum):
    """ACMG/AMP variant classification categories."""
    PATHOGENIC = "pathogenic"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    VUS = "vus"
    LIKELY_BENIGN = "likely_benign"
    BENIGN = "benign"


class VariantType(str, Enum):
    """Types of genetic variants."""
    SNV = "snv"
    INSERTION = "insertion"
    DELETION = "deletion"
    INDEL = "indel"
    CNV = "cnv"
    STRUCTURAL = "structural"
    REPEAT_EXPANSION = "repeat_expansion"


class DiseaseCategory(str, Enum):
    """Broad rare disease category classification."""
    METABOLIC = "metabolic"
    NEUROLOGICAL = "neurological"
    HEMATOLOGIC = "hematologic"
    CONNECTIVE_TISSUE = "connective_tissue"
    IMMUNOLOGIC = "immunologic"
    CARDIAC = "cardiac"
    CANCER_PREDISPOSITION = "cancer_predisposition"
    ENDOCRINE = "endocrine"
    SKELETAL = "skeletal"
    RENAL = "renal"
    PULMONARY = "pulmonary"
    DERMATOLOGIC = "dermatologic"
    OPHTHALMOLOGIC = "ophthalmologic"
    OTHER = "other"


class TherapyStatus(str, Enum):
    """Regulatory approval status for rare disease therapies."""
    APPROVED_FDA = "approved_fda"
    APPROVED_EMA = "approved_ema"
    INVESTIGATIONAL = "investigational"
    COMPASSIONATE_USE = "compassionate_use"
    EXPANDED_ACCESS = "expanded_access"


class Urgency(str, Enum):
    """Clinical urgency level for diagnostic workup."""
    ROUTINE = "routine"
    PRIORITY = "priority"
    EMERGENT = "emergent"


class SeverityLevel(str, Enum):
    """Clinical finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class EvidenceLevel(str, Enum):
    """Strength of evidence supporting a diagnostic finding or recommendation."""
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    CONFLICTING = "conflicting"
    UNCERTAIN = "uncertain"


# ===================================================================
# PYDANTIC MODELS - PATIENT QUERY
# ===================================================================


class PatientQuery(BaseModel):
    """Input query for rare disease diagnostic evaluation.

    Captures patient phenotype (HPO terms), clinical notes, optional
    genomic data (VCF path), and demographic context needed for
    phenotype-driven differential diagnosis.
    """
    patient_id: Optional[str] = Field(
        default=None,
        description="Unique patient identifier",
    )
    hpo_terms: List[str] = Field(
        default_factory=list,
        description="List of HPO term IDs (e.g., ['HP:0001250', 'HP:0001263'])",
    )
    clinical_notes: Optional[str] = Field(
        default=None,
        max_length=10000,
        description="Free-text clinical notes describing the patient presentation",
    )
    vcf_path: Optional[str] = Field(
        default=None,
        description="Path to VCF file for variant-level analysis",
    )
    age: Optional[str] = Field(
        default=None,
        description="Patient age (e.g., '3 years', '6 months', 'neonatal')",
    )
    sex: Optional[str] = Field(
        default=None,
        description="Patient sex (male, female, other, unknown)",
    )
    urgency: Urgency = Field(
        default=Urgency.ROUTINE,
        description="Clinical urgency level for the diagnostic workup",
    )
    workflow_type: Optional[DiagnosticWorkflowType] = Field(
        default=None,
        description="Specific diagnostic workflow; auto-detected if omitted",
    )
    family_history: Optional[str] = Field(
        default=None,
        max_length=5000,
        description="Family history details relevant to genetic diagnosis",
    )
    consanguinity: bool = Field(
        default=False,
        description="Whether consanguinity is reported in the family",
    )
    ethnicity: Optional[str] = Field(
        default=None,
        description="Patient ethnicity (relevant for carrier frequency estimates)",
    )
    prior_testing: Optional[List[str]] = Field(
        default=None,
        description="List of prior genetic/metabolic tests performed",
    )
    top_k: int = Field(
        default=10, ge=1, le=100,
        description="Number of results to return per collection",
    )
    include_guidelines: bool = Field(
        default=True,
        description="Whether to include guideline references in the response",
    )


# ===================================================================
# PYDANTIC MODELS - HPO PHENOTYPE
# ===================================================================


class HPOTerm(BaseModel):
    """Human Phenotype Ontology (HPO) term representation."""
    hpo_id: str = Field(
        ...,
        description="HPO identifier (e.g., 'HP:0001250')",
    )
    name: str = Field(
        ...,
        description="HPO term name (e.g., 'Seizure')",
    )
    definition: Optional[str] = Field(
        default=None,
        description="HPO term definition",
    )
    synonyms: List[str] = Field(
        default_factory=list,
        description="Alternative names for this phenotype",
    )
    information_content: float = Field(
        default=0.0, ge=0.0,
        description="Information content score (higher = more specific)",
    )


# ===================================================================
# PYDANTIC MODELS - DISEASE CANDIDATE
# ===================================================================


class DiseaseCandidate(BaseModel):
    """A candidate rare disease in the differential diagnosis.

    Ranked by phenotypic similarity to the patient presentation,
    with matched and unmatched HPO terms for clinical review.
    """
    disease_id: str = Field(
        ...,
        description="Disease identifier (OMIM or Orphanet)",
    )
    disease_name: str = Field(
        ...,
        description="Disease name",
    )
    rank: int = Field(
        ..., ge=1,
        description="Rank in the differential diagnosis list",
    )
    similarity_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Phenotypic similarity score (0.0 - 1.0)",
    )
    matched_phenotypes: List[str] = Field(
        default_factory=list,
        description="HPO terms matching the patient that are features of this disease",
    )
    unmatched_phenotypes: List[str] = Field(
        default_factory=list,
        description="Disease features not present in the patient (expected but absent)",
    )
    inheritance_pattern: Optional[InheritancePattern] = Field(
        default=None,
        description="Inheritance pattern of the disease",
    )
    prevalence: Optional[str] = Field(
        default=None,
        description="Disease prevalence estimate (e.g., '1/100,000')",
    )
    omim_id: Optional[str] = Field(
        default=None,
        description="OMIM disease identifier",
    )
    orpha_code: Optional[str] = Field(
        default=None,
        description="Orphanet disease code (ORPHA number)",
    )
    category: Optional[DiseaseCategory] = Field(
        default=None,
        description="Broad disease category",
    )
    causal_genes: List[str] = Field(
        default_factory=list,
        description="Known causal genes for this disease",
    )


# ===================================================================
# PYDANTIC MODELS - VARIANT CLASSIFICATION
# ===================================================================


class VariantClassification(BaseModel):
    """ACMG-classified genetic variant with evidence summary.

    Represents a variant found during WES/WGS interpretation
    with its pathogenicity classification and supporting evidence.
    """
    variant_id: str = Field(
        ...,
        description="Variant identifier (chr:pos:ref:alt or dbSNP rsID)",
    )
    gene: str = Field(
        ...,
        description="Gene symbol (HGNC)",
    )
    hgvs_c: Optional[str] = Field(
        default=None,
        description="HGVS coding DNA notation (e.g., 'c.1234A>G')",
    )
    hgvs_p: Optional[str] = Field(
        default=None,
        description="HGVS protein notation (e.g., 'p.Arg412Gly')",
    )
    variant_type: Optional[VariantType] = Field(
        default=None,
        description="Type of genetic variant",
    )
    classification: ACMGClassification = Field(
        ...,
        description="ACMG/AMP pathogenicity classification",
    )
    acmg_criteria: List[str] = Field(
        default_factory=list,
        description="Applied ACMG criteria (e.g., ['PVS1', 'PM2', 'PP3'])",
    )
    population_frequency: Optional[float] = Field(
        default=None, ge=0.0, le=1.0,
        description="Population allele frequency (gnomAD)",
    )
    associated_diseases: List[str] = Field(
        default_factory=list,
        description="Diseases associated with this variant in ClinVar/literature",
    )
    evidence_summary: Optional[str] = Field(
        default=None,
        max_length=2000,
        description="Summary of evidence supporting the classification",
    )
    zygosity: Optional[str] = Field(
        default=None,
        description="Zygosity (heterozygous, homozygous, hemizygous, compound_het)",
    )
    clinvar_stars: Optional[int] = Field(
        default=None, ge=0, le=4,
        description="ClinVar review star rating (0-4)",
    )


# ===================================================================
# PYDANTIC MODELS - THERAPY MATCH
# ===================================================================


class TherapyMatch(BaseModel):
    """A therapeutic option matched to a rare disease diagnosis.

    Includes FDA/EMA-approved orphan drugs, investigational therapies,
    gene therapies, and compassionate use programs.
    """
    therapy_name: str = Field(
        ...,
        description="Drug or therapy name",
    )
    indication: str = Field(
        ...,
        description="Approved or investigational indication",
    )
    status: TherapyStatus = Field(
        ...,
        description="Regulatory approval or access status",
    )
    trial_id: Optional[str] = Field(
        default=None,
        description="ClinicalTrials.gov NCT ID (for investigational therapies)",
    )
    eligibility_criteria: Optional[str] = Field(
        default=None,
        max_length=2000,
        description="Key eligibility criteria for the therapy or trial",
    )
    mechanism: Optional[str] = Field(
        default=None,
        max_length=1000,
        description="Mechanism of action",
    )
    gene_target: Optional[str] = Field(
        default=None,
        description="Target gene (for gene therapies)",
    )
    orphan_designation: bool = Field(
        default=False,
        description="Whether the therapy has orphan drug designation",
    )


# ===================================================================
# PYDANTIC MODELS - SEARCH RESULT
# ===================================================================


class DiagnosticSearchResult(BaseModel):
    """A single search result from any rare disease knowledge collection."""
    collection: str = Field(..., description="Source Milvus collection name")
    content: str = Field(..., description="Retrieved text content")
    score: float = Field(..., ge=0.0, description="Similarity score")
    metadata: Dict = Field(default_factory=dict, description="Source metadata")


# ===================================================================
# PYDANTIC MODELS - DIAGNOSTIC RESULT
# ===================================================================


class DiagnosticResult(BaseModel):
    """Complete diagnostic output combining phenotype matching, variant
    analysis, and therapy identification for a rare disease case.
    """
    candidate_diseases: List[DiseaseCandidate] = Field(
        default_factory=list,
        description="Ranked list of candidate diseases",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Overall diagnostic confidence score",
    )
    variants: List[VariantClassification] = Field(
        default_factory=list,
        description="Classified variants from genomic analysis",
    )
    therapies: List[TherapyMatch] = Field(
        default_factory=list,
        description="Matched therapeutic options",
    )
    evidence_summary: Optional[str] = Field(
        default=None,
        max_length=5000,
        description="Narrative evidence summary synthesizing findings",
    )
    workflow_type: Optional[DiagnosticWorkflowType] = Field(
        default=None,
        description="Diagnostic workflow that produced this result",
    )
    search_results: List[DiagnosticSearchResult] = Field(
        default_factory=list,
        description="Raw search results from knowledge collections",
    )


# ===================================================================
# PYDANTIC MODELS - WORKFLOW RESULT
# ===================================================================


class WorkflowResult(BaseModel):
    """Output from a specific diagnostic workflow execution.

    Captures structured findings, recommendations, guideline references,
    and cross-agent triggers for a single workflow run.
    """
    workflow_type: DiagnosticWorkflowType = Field(
        ...,
        description="Workflow that produced this result",
    )
    findings: List[str] = Field(
        default_factory=list,
        description="Key diagnostic findings",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Recommended next steps (testing, referrals, management)",
    )
    guideline_references: List[str] = Field(
        default_factory=list,
        description="Relevant clinical practice guideline citations",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Overall severity of the findings",
    )
    cross_agent_triggers: List[str] = Field(
        default_factory=list,
        description="Other intelligence agents to consult (e.g., 'cardiology', 'pgx')",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Confidence in the workflow conclusions",
    )
    diagnostic_result: Optional[DiagnosticResult] = Field(
        default=None,
        description="Full diagnostic result (if applicable)",
    )


# ===================================================================
# DATACLASS - SEARCH PLAN
# ===================================================================


@dataclass
class SearchPlan:
    """Plan for searching across multiple Milvus collections.

    Generated by the query router to specify which collections to
    search, how many results to retrieve from each, and the relative
    weights for result merging.
    """
    workflow_type: DiagnosticWorkflowType
    collections: List[str] = field(default_factory=list)
    weights: Dict[str, float] = field(default_factory=dict)
    top_k_per_collection: Dict[str, int] = field(default_factory=dict)
    query_text: str = ""
    hpo_terms: List[str] = field(default_factory=list)
    filters: Dict[str, str] = field(default_factory=dict)
    include_genomic: bool = True
    urgency: Urgency = Urgency.ROUTINE
