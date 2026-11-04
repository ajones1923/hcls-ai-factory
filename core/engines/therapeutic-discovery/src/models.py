"""
Data contracts for the Drug Discovery Pipeline.

Pydantic models defining the data structures passed between pipeline stages.
Based on phase-5-6.pdf specification.
"""
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Any, Optional
from datetime import datetime
from enum import Enum


class TargetStatus(str, Enum):
    """Status of a drug target."""
    HYPOTHESIS = "hypothesis"
    VALIDATED = "validated"
    SELECTED = "selected"
    REJECTED = "rejected"


class Confidence(str, Enum):
    """Confidence level."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class TargetHypothesis(BaseModel):
    """
    Target hypothesis imported from RAG Chat Pipeline.
    Stage 0 input.
    """
    gene: str = Field(..., description="Gene symbol (e.g., VCP)")
    protein: Optional[str] = Field(None, description="Protein name")
    uniprot_id: Optional[str] = Field(None, description="UniProt accession")
    rationale: str = Field("", description="Scientific rationale for targeting")
    therapeutic_area: Optional[str] = Field(None, description="Disease area")
    mechanism: Optional[str] = Field(None, description="Proposed mechanism of action")
    confidence: Confidence = Field(Confidence.MEDIUM)
    priority: int = Field(3, ge=1, le=5, description="Priority 1-5")
    pdb_ids: List[str] = Field(default_factory=list, description="Known PDB structures")
    variant_count: int = Field(0)
    status: TargetStatus = Field(TargetStatus.HYPOTHESIS)

    class Config:
        use_enum_values = True


class StructureInfo(BaseModel):
    """
    Protein structure information.
    Stage 2-3 output.
    """
    pdb_id: str = Field(..., description="PDB identifier")
    method: str = Field("X-ray", description="Experimental method (X-ray, Cryo-EM, NMR)")
    resolution: Optional[float] = Field(None, description="Resolution in Angstroms")
    chain: str = Field("A", description="Chain identifier")
    binding_site_residues: List[str] = Field(default_factory=list)
    ligand_id: Optional[str] = Field(None, description="Co-crystallized ligand")
    ligand_smiles: Optional[str] = Field(None)
    structure_file: Optional[str] = Field(None, description="Path to structure file")
    prepared: bool = Field(False, description="Whether structure has been prepared")
    preparation: Optional[str] = Field(
        None, description="What stage 3 actually did, so `prepared` cannot be read as a claim")

    @validator("resolution", pre=True)
    def parse_resolution(cls, v):
        if isinstance(v, str):
            # Handle "2.3 Å" format
            return float(v.replace("Å", "").strip())
        return v


class StructureManifest(BaseModel):
    """
    Collection of structures for a target.
    Stage 3 output.
    """
    target_gene: str
    structures: List[StructureInfo]
    primary_structure: Optional[str] = Field(None, description="Preferred PDB ID for docking")
    created_at: datetime = Field(default_factory=datetime.now)

    def get_best_structure(self) -> Optional[StructureInfo]:
        """Get the best resolution structure."""
        valid = [s for s in self.structures if s.resolution is not None]
        if not valid:
            return self.structures[0] if self.structures else None
        return min(valid, key=lambda s: s.resolution)


class MoleculeProperties(BaseModel):
    """
    Calculated molecular properties.
    """
    molecular_weight: float = Field(..., ge=0)
    logP: float = Field(..., description="Lipophilicity")
    hbd: int = Field(..., ge=0, description="H-bond donors")
    hba: int = Field(..., ge=0, description="H-bond acceptors")
    tpsa: float = Field(..., ge=0, description="Topological polar surface area")
    rotatable_bonds: int = Field(..., ge=0)
    lipinski_violations: int = Field(0, ge=0, le=4)
    qed: Optional[float] = Field(None, ge=0, le=1, description="Quantitative Estimate of Drug-likeness")
    sa_score: Optional[float] = Field(None, description="Synthetic accessibility score")


class GeneratedMolecule(BaseModel):
    """
    A generated drug candidate.
    Stage 4 output.
    """
    id: str = Field(..., description="Unique molecule identifier")
    smiles: str = Field(..., description="SMILES string")
    name: Optional[str] = Field(None)
    source_seed: Optional[str] = Field(None, description="Seed SMILES if applicable")
    generation_method: str = Field("MolMIM", description="Generation method used")
    target_gene: str
    properties: Optional[MoleculeProperties] = None
    generation_score: float = Field(0.0, ge=0, le=1)
    passed_qc: bool = Field(False)
    conformer_file: Optional[str] = Field(None, description="Path to 3D conformer")
    generated_at: datetime = Field(default_factory=datetime.now)


class DockingResult(BaseModel):
    """
    Result of molecular docking.
    Stage 7 output.
    """
    molecule_id: str
    structure_id: str
    # `docking_score` is a binding energy in kcal/mol (lower is better) and is None when the
    # docking method does not predict one. DiffDock does not: it predicts a POSE and a
    # confidence, so its results carry `confidence_score` instead. Keeping these apart stops
    # a confidence being read (and normalized) as an affinity.
    score_type: str = Field(
        "binding_energy",
        description="'binding_energy' (kcal/mol, lower better) or 'diffdock_confidence' (log-odds, higher better)",
    )
    docking_score: Optional[float] = Field(None, description="Binding energy kcal/mol (lower is better)")
    confidence_score: Optional[float] = Field(
        None, description="Raw pose confidence (DiffDock log-odds; higher is better)"
    )
    binding_energy: Optional[float] = Field(None, description="Binding energy kcal/mol")
    pose_file: Optional[str] = Field(None, description="Path to docked pose")
    rmsd: Optional[float] = Field(None, description="RMSD to reference pose")
    contacts: List[str] = Field(default_factory=list, description="Residues in contact")
    hydrogen_bonds: int = Field(0)
    method: str = Field("DiffDock", description="Docking method used")
    confidence: float = Field(0.0, ge=0, le=1, description="Calibrated P(pose within 2 A RMSD)")
    provenance: str = Field(
        "real",
        description="'real' or 'mock'. The mock client stamps _provenance on its poses; this "
                    "carries it into the record so a report cannot drop it.",
    )


class PediatricSafetyFlag(BaseModel):
    """A safety flag raised during pediatric assessment."""
    flag: str = Field(..., description="Flag identifier (e.g., HIGH_MW_BBB_RISK)")
    detail: str = Field(..., description="Human-readable explanation")
    severity: str = Field(..., description="critical, high, medium, or low")


class PediatricSafetyAssessment(BaseModel):
    """
    Pediatric safety assessment result for a drug candidate.

    `pediatric_safe` is three-valued. It was a bare bool, which forced the assessment to answer
    "safe" whenever no critical flag had been raised -- including when the only critical rule
    (cardiac liability) had no input to evaluate. `None` means "not assessed", and
    `cardiac_liability_assessed` records whether a hERG prediction was available at all.
    """
    pediatric_safe: Optional[bool] = Field(
        ..., description="True if no critical flags; None if the assessment could not be made")
    safety_flags: List[PediatricSafetyFlag] = Field(default_factory=list)
    recommendation: str = Field("", description="Summary recommendation")
    cardiac_liability_assessed: bool = Field(
        False, description="True only if a hERG prediction from a model was available")
    unassessed_reasons: List[str] = Field(
        default_factory=list, description="Why the verdict is None, when it is")


class RankedCandidate(BaseModel):
    """
    A ranked drug candidate with all scores.
    Stage 8 output.
    """
    rank: int = Field(..., ge=1)
    molecule_id: str
    smiles: str
    name: Optional[str] = None
    target_gene: str

    # Scores
    score_type: str = Field("binding_energy", description="Units of docking_score / confidence_score")
    docking_score: Optional[float] = None       # kcal/mol; None if the method predicts no affinity
    confidence_score: Optional[float] = None    # raw pose confidence (DiffDock log-odds)
    generation_score: float
    qed_score: Optional[float] = None
    composite_score: float = Field(..., description="Weighted composite score")

    # Provenance of the docking term in `composite_score`, carried so a report cannot present a
    # mock pose as a measurement by omission.
    docking_provenance: str = Field("none", description="'real', 'mock' or 'none' (not docked)")

    # Properties
    properties: Optional[MoleculeProperties] = None

    # Docking details
    best_pose_file: Optional[str] = None
    binding_residues: List[str] = Field(default_factory=list)

    # Pediatric safety
    pediatric_safety: Optional[PediatricSafetyAssessment] = Field(
        None, description="Pediatric safety assessment result"
    )

    # Flags
    passes_filters: bool = Field(True)
    alerts: List[str] = Field(default_factory=list, description="Any warnings")


class PipelineRun(BaseModel):
    """
    Metadata for a pipeline run.
    """
    run_id: str
    target_gene: str
    started_at: datetime
    completed_at: Optional[datetime] = None
    status: str = "running"  # running, completed, failed
    current_stage: int = 0
    stages_completed: List[int] = Field(default_factory=list)
    config: Dict[str, Any] = Field(default_factory=dict)
    error_message: Optional[str] = None

    # Results summary
    molecules_generated: int = 0
    molecules_docked: int = 0
    top_candidates: List[str] = Field(default_factory=list)

    # Performance
    stage_timings: Dict[int, float] = Field(default_factory=dict, description="Stage durations in seconds")


class PipelineConfig(BaseModel):
    """
    Pipeline configuration.
    """
    # Target settings
    target_gene: str
    reference_compound_smiles: Optional[str] = None

    # Generation settings
    num_molecules: int = Field(50, ge=1, le=1000)
    diversity: float = Field(0.3, ge=0, le=1)
    generation_method: str = Field("molmim")

    # Docking settings
    docking_method: str = Field("diffdock")
    num_poses: int = Field(10, ge=1, le=100)

    # Filtering
    max_mw: float = Field(550, description="Max molecular weight")
    max_logp: float = Field(5.0)
    max_lipinski_violations: int = Field(1)

    # Ranking
    docking_weight: float = Field(0.4)
    generation_weight: float = Field(0.3)
    qed_weight: float = Field(0.3)
    top_n_candidates: int = Field(10)

    # Paths
    output_dir: str = Field("outputs")
    structure_cache: str = Field("data/structures")

    # Service endpoints
    molmim_url: Optional[str] = Field(None)
    diffdock_url: Optional[str] = Field(None)
