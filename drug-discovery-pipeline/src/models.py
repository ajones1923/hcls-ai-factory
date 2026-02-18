"""
Data contracts for the Drug Discovery Pipeline.

Pydantic models defining the data structures passed between pipeline stages.
Based on phase-5-6.pdf specification.
"""
from datetime import datetime
from enum import Enum
from typing import Any

from pydantic import BaseModel, Field, validator


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
    protein: str | None = Field(None, description="Protein name")
    uniprot_id: str | None = Field(None, description="UniProt accession")
    rationale: str = Field("", description="Scientific rationale for targeting")
    therapeutic_area: str | None = Field(None, description="Disease area")
    mechanism: str | None = Field(None, description="Proposed mechanism of action")
    confidence: Confidence = Field(Confidence.MEDIUM)
    priority: int = Field(3, ge=1, le=5, description="Priority 1-5")
    pdb_ids: list[str] = Field(default_factory=list, description="Known PDB structures")
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
    resolution: float | None = Field(None, description="Resolution in Angstroms")
    chain: str = Field("A", description="Chain identifier")
    binding_site_residues: list[str] = Field(default_factory=list)
    ligand_id: str | None = Field(None, description="Co-crystallized ligand")
    ligand_smiles: str | None = Field(None)
    structure_file: str | None = Field(None, description="Path to structure file")
    prepared: bool = Field(False, description="Whether structure has been prepared")

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
    structures: list[StructureInfo]
    primary_structure: str | None = Field(None, description="Preferred PDB ID for docking")
    created_at: datetime = Field(default_factory=datetime.now)

    def get_best_structure(self) -> StructureInfo | None:
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
    qed: float | None = Field(None, ge=0, le=1, description="Quantitative Estimate of Drug-likeness")
    sa_score: float | None = Field(None, description="Synthetic accessibility score")


class GeneratedMolecule(BaseModel):
    """
    A generated drug candidate.
    Stage 4 output.
    """
    id: str = Field(..., description="Unique molecule identifier")
    smiles: str = Field(..., description="SMILES string")
    name: str | None = Field(None)
    source_seed: str | None = Field(None, description="Seed SMILES if applicable")
    generation_method: str = Field("MolMIM", description="Generation method used")
    target_gene: str
    properties: MoleculeProperties | None = None
    generation_score: float = Field(0.0, ge=0, le=1)
    passed_qc: bool = Field(False)
    conformer_file: str | None = Field(None, description="Path to 3D conformer")
    generated_at: datetime = Field(default_factory=datetime.now)


class DockingResult(BaseModel):
    """
    Result of molecular docking.
    Stage 7 output.
    """
    molecule_id: str
    structure_id: str
    docking_score: float = Field(..., description="Docking score (lower is better)")
    binding_energy: float | None = Field(None, description="Binding energy kcal/mol")
    pose_file: str | None = Field(None, description="Path to docked pose")
    rmsd: float | None = Field(None, description="RMSD to reference pose")
    contacts: list[str] = Field(default_factory=list, description="Residues in contact")
    hydrogen_bonds: int = Field(0)
    method: str = Field("DiffDock", description="Docking method used")
    confidence: float = Field(0.0, ge=0, le=1)


class RankedCandidate(BaseModel):
    """
    A ranked drug candidate with all scores.
    Stage 8 output.
    """
    rank: int = Field(..., ge=1)
    molecule_id: str
    smiles: str
    name: str | None = None
    target_gene: str

    # Scores
    docking_score: float
    generation_score: float
    qed_score: float | None = None
    composite_score: float = Field(..., description="Weighted composite score")

    # Properties
    properties: MoleculeProperties | None = None

    # Docking details
    best_pose_file: str | None = None
    binding_residues: list[str] = Field(default_factory=list)

    # Flags
    passes_filters: bool = Field(True)
    alerts: list[str] = Field(default_factory=list, description="Any warnings")


class PipelineRun(BaseModel):
    """
    Metadata for a pipeline run.
    """
    run_id: str
    target_gene: str
    started_at: datetime
    completed_at: datetime | None = None
    status: str = "running"  # running, completed, failed
    current_stage: int = 0
    stages_completed: list[int] = Field(default_factory=list)
    config: dict[str, Any] = Field(default_factory=dict)
    error_message: str | None = None

    # Results summary
    molecules_generated: int = 0
    molecules_docked: int = 0
    top_candidates: list[str] = Field(default_factory=list)

    # Performance
    stage_timings: dict[int, float] = Field(default_factory=dict, description="Stage durations in seconds")


class PipelineConfig(BaseModel):
    """
    Pipeline configuration.
    """
    # Target settings
    target_gene: str
    reference_compound_smiles: str | None = None

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
    molmim_url: str | None = Field(None)
    diffdock_url: str | None = Field(None)
