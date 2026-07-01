"""Single-cell clinical API routes.

Provides endpoints for RAG-powered single-cell queries, cell type
annotation, tumor microenvironment profiling, drug response prediction,
subclonal architecture analysis, spatial niche mapping, trajectory
inference, ligand-receptor interaction, biomarker discovery, CAR-T
target validation, treatment monitoring, and reference catalogues.

Author: Adam Jones
Date: March 2026
"""

from typing import List, Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

from src.knowledge import (
    KNOWLEDGE_VERSION,
    CELL_TYPE_ATLAS,
    TME_PROFILES,
    DRUG_SENSITIVITY_DATABASE,
    LIGAND_RECEPTOR_PAIRS,
    CANCER_TME_ATLAS,
    IMMUNE_SIGNATURES,
    MARKER_GENE_DATABASE,
)

router = APIRouter(prefix="/v1/sc", tags=["single-cell"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, CAR-T, biomarker, drug discovery, and imaging agents
    for a comprehensive single-cell-informed assessment.
    """
    try:
        from src.cross_modal import (
            query_oncology_agent,
            query_cart_agent,
            query_biomarker_agent,
            query_drug_discovery_agent,
            query_imaging_agent,
            integrate_cross_agent_results,
        )

        tumor_data = request.get("tumor_data", {})
        target_data = request.get("target_data", {})
        biomarker_data = request.get("biomarker_data", {})
        drug_data = request.get("drug_data", {})
        spatial_data = request.get("spatial_data", {})

        results = []

        # Query oncology agent for tumor profiling correlation
        if tumor_data:
            results.append(query_oncology_agent(tumor_data))

        # Query CAR-T agent for target validation
        if target_data:
            results.append(query_cart_agent(target_data))

        # Query biomarker agent for panel enrichment
        if biomarker_data:
            results.append(query_biomarker_agent(biomarker_data))

        # Query drug discovery for compound integration
        if drug_data:
            results.append(query_drug_discovery_agent(drug_data))

        # Query imaging agent for spatial-imaging correlation
        if spatial_data:
            results.append(query_imaging_agent(spatial_data))

        integrated = integrate_cross_agent_results(results)
        return {
            "status": "completed",
            "assessment": integrated,
            "agents_consulted": integrated.get("agents_consulted", []),
        }
    except Exception as exc:
        logger.error(f"Integrated assessment failed: {exc}")
        return {"status": "partial", "assessment": {}, "error": "Cross-agent integration unavailable"}


# =====================================================================
# Request / Response Schemas
# =====================================================================

# -- Query --

class QueryRequest(BaseModel):
    """Free-text RAG query with optional domain and patient context."""
    question: str = Field(..., min_length=3, description="Single-cell biology question")
    domain: Optional[str] = Field(
        None,
        description=(
            "Domain hint: transcriptomics | immunology | oncology | "
            "spatial | multiomics | pharmacology | cell_therapy | general"
        ),
    )
    patient_context: Optional[dict] = Field(None, description="Tumor type, sample info, clinical data")
    top_k: int = Field(5, ge=1, le=50, description="Number of evidence passages")
    include_guidelines: bool = Field(True, description="Include guideline citations")


class QueryResponse(BaseModel):
    answer: str
    evidence: List[dict]
    guidelines_cited: List[str] = []
    confidence: float
    domain_applied: Optional[str] = None


class SearchRequest(BaseModel):
    """Multi-collection semantic search."""
    question: str = Field(..., min_length=3)
    collections: Optional[List[str]] = None
    top_k: int = Field(5, ge=1, le=100)
    threshold: float = Field(0.0, ge=0.0, le=1.0)


class SearchResult(BaseModel):
    collection: str
    text: str
    score: float
    metadata: dict = {}


class SearchResponse(BaseModel):
    results: List[SearchResult]
    total: int
    collections_searched: List[str]


# -- Cell Type Annotation --

class AnnotateRequest(BaseModel):
    """Cell type annotation request."""
    gene_expression: Optional[dict] = Field(
        None,
        description="Gene expression profile (gene -> expression value)",
    )
    marker_genes: Optional[List[str]] = Field(
        None,
        description="Detected marker genes for annotation",
    )
    tissue_type: Optional[str] = Field(None, description="Source tissue (e.g., lung, breast, PBMC)")
    num_cells: Optional[int] = Field(None, ge=1, description="Number of cells in cluster")
    cluster_id: Optional[str] = Field(None, description="Cluster identifier")
    strategy: str = Field(
        "marker_based",
        description="Annotation strategy: reference_based | marker_based | llm_based",
    )
    top_n: int = Field(5, ge=1, le=20, description="Max cell types to return")
    clinical_notes: Optional[str] = Field(None)


class AnnotateResponse(BaseModel):
    cell_types: List[dict]
    cluster_id: Optional[str] = None
    tissue_type: Optional[str] = None
    strategy_used: str
    num_cells: Optional[int] = None
    confidence: float
    marker_evidence: List[dict] = []
    recommendations: List[str] = []


# -- TME Profiling --

class TMEProfileRequest(BaseModel):
    """Tumor microenvironment profiling request."""
    cell_type_proportions: Optional[dict] = Field(
        None,
        description="Cell type -> proportion mapping from deconvolution",
    )
    immune_markers: Optional[dict] = Field(None, description="Immune marker expression levels")
    tumor_type: Optional[str] = Field(None, description="Cancer type (e.g., NSCLC, melanoma)")
    spatial_data: Optional[dict] = Field(None, description="Spatial distribution data")
    gene_signatures: Optional[dict] = Field(None, description="Gene signature scores")
    clinical_notes: Optional[str] = Field(None)


class TMEProfileResponse(BaseModel):
    tme_classification: str
    immune_score: float
    stromal_score: float
    cell_composition: dict
    immune_phenotype: str
    checkpoint_expression: dict = {}
    therapy_prediction: dict = {}
    recommendations: List[str] = []
    evidence: List[dict] = []
    guidelines_cited: List[str] = []


# -- Drug Response --

class DrugResponseRequest(BaseModel):
    """Drug response prediction request."""
    gene_expression: Optional[dict] = Field(None, description="Gene expression profile")
    cell_type: Optional[str] = Field(None, description="Target cell type")
    drug_name: Optional[str] = Field(None, description="Drug/compound to evaluate")
    drug_class: Optional[str] = Field(None, description="Drug class (e.g., checkpoint inhibitor, TKI)")
    tumor_type: Optional[str] = Field(None, description="Cancer type")
    genomic_alterations: Optional[List[str]] = Field(None, description="Known mutations/CNVs")
    clinical_notes: Optional[str] = Field(None)


class DrugResponseResponse(BaseModel):
    predictions: List[dict]
    sensitivity_score: float
    resistance_mechanisms: List[str] = []
    biomarkers: List[str] = []
    recommendations: List[str] = []
    evidence: List[dict] = []
    guidelines_cited: List[str] = []


# -- Subclonal Analysis --

class SubclonalRequest(BaseModel):
    """Subclonal architecture analysis request."""
    cnv_profile: Optional[dict] = Field(None, description="Copy number variation data")
    mutation_data: Optional[dict] = Field(None, description="Somatic mutation data")
    cell_count: Optional[int] = Field(None, ge=1, description="Total cells analyzed")
    tumor_type: Optional[str] = Field(None)
    target_antigen: Optional[str] = Field(None, description="CAR-T/therapy target antigen")
    clinical_notes: Optional[str] = Field(None)


class SubclonalResponse(BaseModel):
    clones: List[dict]
    num_clones: int
    dominant_clone: Optional[str] = None
    heterogeneity_index: float
    escape_risk: str
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- Spatial Niche --

class SpatialNicheRequest(BaseModel):
    """Spatial niche mapping request."""
    spatial_coordinates: Optional[dict] = Field(None, description="Cell coordinate data")
    cell_types: Optional[dict] = Field(None, description="Cell type assignments per spot/cell")
    platform: str = Field("visium", description="Spatial platform: visium | merfish | cosmx | xenium | slide_seq | stereo_seq | seqfish")
    gene_expression: Optional[dict] = Field(None, description="Spatially-resolved expression")
    tissue_type: Optional[str] = Field(None)
    clinical_notes: Optional[str] = Field(None)


class SpatialNicheResponse(BaseModel):
    niches: List[dict]
    num_niches: int
    spatial_statistics: dict = {}
    cell_cell_interactions: List[dict] = []
    platform: str
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- Trajectory --

class TrajectoryRequest(BaseModel):
    """Trajectory inference request."""
    gene_expression: Optional[dict] = Field(None, description="Expression matrix summary")
    cell_types: Optional[List[str]] = Field(None, description="Cell types in trajectory")
    root_cell_type: Optional[str] = Field(None, description="Root/stem cell type")
    method: str = Field("monocle3", description="Method: monocle3 | paga | rna_velocity | scvelo | cellrank | palantir")
    tissue_type: Optional[str] = Field(None)
    clinical_notes: Optional[str] = Field(None)


class TrajectoryResponse(BaseModel):
    trajectory: dict
    branch_points: List[dict] = []
    driver_genes: List[dict] = []
    pseudotime_range: List[float] = []
    method_used: str
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- Ligand-Receptor --

class LigandReceptorRequest(BaseModel):
    """Ligand-receptor interaction analysis request."""
    cell_types: Optional[List[str]] = Field(None, description="Cell types to analyze")
    gene_expression: Optional[dict] = Field(None, description="Expression data")
    source_cell_type: Optional[str] = Field(None, description="Ligand source cell type")
    target_cell_type: Optional[str] = Field(None, description="Receptor target cell type")
    database: str = Field("cellphonedb", description="Database: cellphonedb | nichenet | cellchat | celltalkdb")
    tissue_type: Optional[str] = Field(None)
    clinical_notes: Optional[str] = Field(None)


class LigandReceptorResponse(BaseModel):
    interactions: List[dict]
    num_significant: int
    pathways_enriched: List[dict] = []
    network_summary: dict = {}
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- Biomarker --

class BiomarkerRequest(BaseModel):
    """Biomarker discovery request."""
    differential_expression: Optional[dict] = Field(None, description="DE results (gene -> stats)")
    cell_type: Optional[str] = Field(None, description="Cell type of interest")
    condition_a: Optional[str] = Field(None, description="Comparison group A")
    condition_b: Optional[str] = Field(None, description="Comparison group B")
    tumor_type: Optional[str] = Field(None)
    min_log2fc: float = Field(1.0, description="Minimum log2 fold change threshold")
    max_pval: float = Field(0.05, description="Maximum adjusted p-value threshold")
    clinical_notes: Optional[str] = Field(None)


class BiomarkerResponse(BaseModel):
    biomarkers: List[dict]
    num_candidates: int
    top_markers: List[str] = []
    validation_suggestions: List[str] = []
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- CAR-T Validation --

class CARTValidateRequest(BaseModel):
    """CAR-T target validation request."""
    target_gene: str = Field(..., min_length=1, description="Target antigen gene symbol")
    tumor_type: Optional[str] = Field(None, description="Cancer type")
    expression_data: Optional[dict] = Field(None, description="Expression across cell types")
    tme_data: Optional[dict] = Field(None, description="TME composition data")
    normal_tissue_expression: Optional[dict] = Field(None, description="Normal tissue expression")
    clinical_notes: Optional[str] = Field(None)


class CARTValidateResponse(BaseModel):
    target_gene: str
    on_tumor_pct: float
    off_tumor_risk: dict = {}
    tme_compatibility: float
    escape_risk: str
    therapeutic_index: float
    subclonal_heterogeneity: Optional[dict] = None
    recommendations: List[str] = []
    evidence: List[dict] = []
    guidelines_cited: List[str] = []


# -- Treatment Monitor --

class TreatmentMonitorRequest(BaseModel):
    """Treatment monitoring request."""
    timepoints: Optional[List[dict]] = Field(None, description="Longitudinal data per timepoint")
    treatment: Optional[str] = Field(None, description="Treatment regimen")
    tumor_type: Optional[str] = Field(None)
    baseline_composition: Optional[dict] = Field(None, description="Pre-treatment cell composition")
    current_composition: Optional[dict] = Field(None, description="Current cell composition")
    target_antigen: Optional[str] = Field(None, description="Therapy target antigen")
    clinical_notes: Optional[str] = Field(None)


class TreatmentMonitorResponse(BaseModel):
    response_assessment: str
    composition_changes: dict = {}
    resistance_indicators: List[dict] = []
    clone_dynamics: List[dict] = []
    immune_shift: dict = {}
    recommendations: List[str] = []
    evidence: List[dict] = []


# -- Workflow (generic) --

class WorkflowRequest(BaseModel):
    """Generic workflow execution request."""
    data: dict = Field(default={}, description="Workflow-specific input data")
    question: Optional[str] = Field(None, description="Free-text question for the workflow")


# =====================================================================
# Helper functions
# =====================================================================

def _get_engine(request: Request):
    """Get the RAG engine from app state."""
    engine = getattr(request.app.state, "engine", None)
    return engine


def _get_llm(request: Request):
    """Get the LLM client from app state."""
    return getattr(request.app.state, "llm_client", None)


def _get_workflow_engine(request: Request):
    """Get the workflow engine from app state."""
    return getattr(request.app.state, "workflow_engine", None)


def _increment_metric(request: Request, metric_name: str):
    """Increment a metrics counter."""
    metrics = getattr(request.app.state, "metrics", None)
    lock = getattr(request.app.state, "metrics_lock", None)
    if metrics and lock:
        with lock:
            metrics[metric_name] = metrics.get(metric_name, 0) + 1


async def _llm_fallback(request: Request, prompt: str, system_prompt: str = "") -> str:
    """Generate an LLM response as fallback when engine is unavailable."""
    llm = _get_llm(request)
    if not llm:
        return "LLM unavailable. Please ensure ANTHROPIC_API_KEY is configured."
    try:
        return llm.generate(prompt, system_prompt=system_prompt)
    except Exception as exc:
        logger.error(f"LLM fallback failed: {exc}")
        return f"LLM generation failed: {exc}"


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/query", response_model=QueryResponse)
async def query(request: QueryRequest, req: Request):
    """RAG-powered Q&A query across single-cell knowledge collections."""
    _increment_metric(req, "query_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.query(
                question=request.question,
                domain=request.domain,
                patient_context=request.patient_context,
                top_k=request.top_k,
            )
            return QueryResponse(**result)
        except Exception as exc:
            logger.error(f"RAG query failed: {exc}")

    # LLM fallback
    system_prompt = (
        "You are a single-cell genomics intelligence system. "
        "Provide evidence-based answers about single-cell analysis, "
        "cell types, gene expression, TME profiling, and spatial "
        "transcriptomics."
    )
    answer = await _llm_fallback(req, request.question, system_prompt)
    return QueryResponse(
        answer=answer,
        evidence=[],
        guidelines_cited=[],
        confidence=0.3,
        domain_applied=request.domain,
    )


@router.post("/search", response_model=SearchResponse)
async def search(request: SearchRequest, req: Request):
    """Multi-collection semantic search across single-cell collections."""
    _increment_metric(req, "search_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            results = engine.search(
                query=request.question,
                collections=request.collections,
                top_k=request.top_k,
                threshold=request.threshold,
            )
            return SearchResponse(
                results=[SearchResult(**r) for r in results],
                total=len(results),
                collections_searched=request.collections or ["all"],
            )
        except Exception as exc:
            logger.error(f"Search failed: {exc}")

    return SearchResponse(results=[], total=0, collections_searched=[])


@router.post("/annotate", response_model=AnnotateResponse)
async def annotate_cell_types(request: AnnotateRequest, req: Request):
    """Annotate cell types from marker genes or expression profiles."""
    _increment_metric(req, "annotate_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.annotate_cell_types(
                gene_expression=request.gene_expression,
                marker_genes=request.marker_genes,
                tissue_type=request.tissue_type,
                strategy=request.strategy,
                top_n=request.top_n,
            )
            return AnnotateResponse(**result)
        except Exception as exc:
            logger.error(f"Cell type annotation failed: {exc}")

    # LLM fallback for annotation
    markers_str = ", ".join(request.marker_genes) if request.marker_genes else "not provided"
    prompt = (
        f"Annotate cell types based on the following information:\n"
        f"Marker genes: {markers_str}\n"
        f"Tissue: {request.tissue_type or 'unknown'}\n"
        f"Strategy: {request.strategy}\n\n"
        f"Provide the most likely cell type annotations with confidence."
    )
    answer = await _llm_fallback(req, prompt)

    # Knowledge base fallback: match marker genes against CELL_TYPE_ATLAS
    kb_cell_types = []
    marker_evidence = []
    if request.marker_genes:
        input_markers = set(g.upper() for g in request.marker_genes)
        for ct_name, ct_data in CELL_TYPE_ATLAS.items():
            ct_markers = set(m.upper() for m in ct_data.get("markers", []))
            overlap = input_markers & ct_markers
            if overlap:
                score = len(overlap) / max(len(ct_markers), 1)
                kb_cell_types.append({
                    "cell_type": ct_name.replace("_", " "),
                    "confidence": round(score, 2),
                    "markers": list(overlap),
                    "source": "knowledge_base",
                })
                marker_evidence.append({
                    "cell_type": ct_name,
                    "matched_markers": list(overlap),
                    "total_markers": len(ct_markers),
                })
        kb_cell_types.sort(key=lambda x: x["confidence"], reverse=True)
        kb_cell_types = kb_cell_types[:request.top_n]

    if not kb_cell_types:
        kb_cell_types = [{"cell_type": "unknown", "confidence": 0.0, "markers": request.marker_genes or []}]

    return AnnotateResponse(
        cell_types=kb_cell_types,
        cluster_id=request.cluster_id,
        tissue_type=request.tissue_type,
        strategy_used=request.strategy,
        num_cells=request.num_cells,
        confidence=kb_cell_types[0]["confidence"] if kb_cell_types else 0.2,
        marker_evidence=marker_evidence,
        recommendations=[answer] if answer and "unavailable" not in answer.lower() else [],
    )


@router.post("/tme-profile", response_model=TMEProfileResponse)
async def tme_profile(request: TMEProfileRequest, req: Request):
    """Profile the tumor microenvironment from single-cell data."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.profile_tme(
                cell_type_proportions=request.cell_type_proportions,
                immune_markers=request.immune_markers,
                tumor_type=request.tumor_type,
            )
            return TMEProfileResponse(**result)
        except Exception as exc:
            logger.error(f"TME profiling failed: {exc}")

    prompt = (
        f"Classify the tumor microenvironment:\n"
        f"Cell proportions: {request.cell_type_proportions}\n"
        f"Tumor type: {request.tumor_type or 'unknown'}\n"
        f"Classify as hot/cold/excluded/immunosuppressive and provide therapy recommendations."
    )
    answer = await _llm_fallback(req, prompt)

    # Knowledge base fallback: match tumor type against TME_PROFILES and CANCER_TME_ATLAS
    kb_tme_class = "unknown"
    kb_therapy = {}
    kb_evidence = []
    if request.tumor_type:
        tumor_key = request.tumor_type.lower().replace(" ", "_")
        for atlas_key, atlas_data in CANCER_TME_ATLAS.items():
            if tumor_key in atlas_key.lower():
                kb_tme_class = atlas_data.get("dominant_tme", "unknown")
                kb_therapy = atlas_data.get("therapy_implications", {})
                kb_evidence.append({"source": "cancer_tme_atlas", "entry": atlas_key})
                break
    for tme_key, tme_data in TME_PROFILES.items():
        kb_evidence.append({"source": "tme_profiles", "class": tme_key, "description": tme_data.get("description", "")})

    return TMEProfileResponse(
        tme_classification=kb_tme_class,
        immune_score=0.0,
        stromal_score=0.0,
        cell_composition=request.cell_type_proportions or {},
        immune_phenotype="unclassified",
        therapy_prediction=kb_therapy,
        evidence=kb_evidence[:5],
        recommendations=[answer] if answer and "unavailable" not in answer.lower() else [],
    )


@router.post("/drug-response", response_model=DrugResponseResponse)
async def drug_response(request: DrugResponseRequest, req: Request):
    """Predict drug response from single-cell expression profiles."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.predict_drug_response(
                gene_expression=request.gene_expression,
                cell_type=request.cell_type,
                drug_name=request.drug_name,
                tumor_type=request.tumor_type,
            )
            return DrugResponseResponse(**result)
        except Exception as exc:
            logger.error(f"Drug response prediction failed: {exc}")

    prompt = (
        f"Predict drug response:\n"
        f"Drug: {request.drug_name or 'not specified'}\n"
        f"Cell type: {request.cell_type or 'not specified'}\n"
        f"Tumor type: {request.tumor_type or 'unknown'}\n"
        f"Genomic alterations: {request.genomic_alterations or 'none'}\n"
        f"Provide sensitivity prediction and resistance mechanisms."
    )
    answer = await _llm_fallback(req, prompt)

    # Knowledge base fallback: match drug against DRUG_SENSITIVITY_DATABASE
    kb_predictions = []
    kb_biomarkers = []
    kb_resistance = []
    if request.drug_name:
        drug_key = request.drug_name.lower().replace(" ", "_")
        for db_key, db_data in DRUG_SENSITIVITY_DATABASE.items():
            if drug_key in db_key.lower() or db_key.lower() in drug_key:
                kb_predictions.append({
                    "drug": db_key,
                    "drug_class": db_data.get("drug_class", ""),
                    "targets": db_data.get("targets", []),
                    "sensitive_cell_types": db_data.get("sensitive_cell_types", []),
                    "source": "knowledge_base",
                })
                kb_biomarkers.extend(db_data.get("biomarkers", []))
                kb_resistance.extend(db_data.get("resistance_mechanisms", []))
    if not kb_predictions and request.drug_class:
        for db_key, db_data in DRUG_SENSITIVITY_DATABASE.items():
            if request.drug_class.lower() in db_data.get("drug_class", "").lower():
                kb_predictions.append({
                    "drug": db_key,
                    "drug_class": db_data.get("drug_class", ""),
                    "targets": db_data.get("targets", []),
                    "source": "knowledge_base",
                })

    return DrugResponseResponse(
        predictions=kb_predictions[:5],
        sensitivity_score=0.0,
        biomarkers=list(set(kb_biomarkers))[:10],
        resistance_mechanisms=list(set(kb_resistance))[:5],
        recommendations=[answer] if answer and "unavailable" not in answer.lower() else [],
    )


@router.post("/subclonal", response_model=SubclonalResponse)
async def subclonal_analysis(request: SubclonalRequest, req: Request):
    """Analyze subclonal architecture from single-cell data."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.analyze_subclones(
                cnv_profile=request.cnv_profile,
                mutation_data=request.mutation_data,
                target_antigen=request.target_antigen,
            )
            return SubclonalResponse(**result)
        except Exception as exc:
            logger.error(f"Subclonal analysis failed: {exc}")

    # Knowledge base fallback: use CANCER_TME_ATLAS for tumor heterogeneity context
    kb_recs = ["Connect RAG engine for full subclonal analysis."]
    kb_evidence = []
    if request.tumor_type:
        tumor_key = request.tumor_type.lower().replace(" ", "_")
        for atlas_key, atlas_data in CANCER_TME_ATLAS.items():
            if tumor_key in atlas_key.lower():
                kb_recs.append(f"Reference TME profile: {atlas_data.get('dominant_tme', 'N/A')}")
                kb_evidence.append({"source": "cancer_tme_atlas", "entry": atlas_key})
                break

    return SubclonalResponse(
        clones=[],
        num_clones=0,
        heterogeneity_index=0.0,
        escape_risk="unknown",
        recommendations=kb_recs,
        evidence=kb_evidence,
    )


@router.post("/spatial-niche", response_model=SpatialNicheResponse)
async def spatial_niche(request: SpatialNicheRequest, req: Request):
    """Map spatial niches from spatial transcriptomics data."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.map_spatial_niches(
                spatial_coordinates=request.spatial_coordinates,
                cell_types=request.cell_types,
                platform=request.platform,
            )
            return SpatialNicheResponse(**result)
        except Exception as exc:
            logger.error(f"Spatial niche mapping failed: {exc}")

    # Knowledge base fallback: populate from LIGAND_RECEPTOR_PAIRS for cell interactions
    kb_interactions = []
    if request.cell_types:
        for lr_key, lr_data in LIGAND_RECEPTOR_PAIRS.items():
            source_types = set(ct.lower() for ct in lr_data.get("source_cell_types", []))
            target_types = set(ct.lower() for ct in lr_data.get("target_cell_types", []))
            input_types = set(ct.lower() for ct in (request.cell_types.values() if isinstance(request.cell_types, dict) else []))
            if input_types & (source_types | target_types):
                kb_interactions.append({
                    "ligand": lr_data.get("ligand", lr_key),
                    "receptor": lr_data.get("receptor", ""),
                    "source": "knowledge_base",
                })

    return SpatialNicheResponse(
        niches=[],
        num_niches=0,
        cell_cell_interactions=kb_interactions[:10],
        platform=request.platform,
        recommendations=["Connect RAG engine for full spatial niche analysis."],
    )


@router.post("/trajectory", response_model=TrajectoryResponse)
async def trajectory_inference(request: TrajectoryRequest, req: Request):
    """Infer cell trajectories and pseudotime ordering."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.infer_trajectory(
                cell_types=request.cell_types,
                root_cell_type=request.root_cell_type,
                method=request.method,
            )
            return TrajectoryResponse(**result)
        except Exception as exc:
            logger.error(f"Trajectory inference failed: {exc}")

    # Knowledge base fallback: provide driver genes from MARKER_GENE_DATABASE
    kb_driver_genes = []
    if request.cell_types:
        for ct in request.cell_types:
            ct_key = ct.lower().replace(" ", "_")
            for mk_key, mk_data in MARKER_GENE_DATABASE.items():
                if ct_key in mk_key.lower():
                    for gene in mk_data.get("markers", [])[:3]:
                        kb_driver_genes.append({
                            "gene": gene,
                            "cell_type": mk_key,
                            "source": "knowledge_base",
                        })

    return TrajectoryResponse(
        trajectory={},
        driver_genes=kb_driver_genes[:10],
        method_used=request.method,
        recommendations=["Connect RAG engine for full trajectory inference."],
    )


@router.post("/ligand-receptor", response_model=LigandReceptorResponse)
async def ligand_receptor(request: LigandReceptorRequest, req: Request):
    """Analyze ligand-receptor interactions between cell types."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.analyze_ligand_receptor(
                cell_types=request.cell_types,
                source_cell_type=request.source_cell_type,
                target_cell_type=request.target_cell_type,
                database=request.database,
            )
            return LigandReceptorResponse(**result)
        except Exception as exc:
            logger.error(f"Ligand-receptor analysis failed: {exc}")

    # Knowledge base fallback: match cell types against LIGAND_RECEPTOR_PAIRS
    kb_interactions = []
    if request.cell_types or request.source_cell_type or request.target_cell_type:
        query_types = set()
        if request.cell_types:
            query_types.update(ct.lower() for ct in request.cell_types)
        if request.source_cell_type:
            query_types.add(request.source_cell_type.lower())
        if request.target_cell_type:
            query_types.add(request.target_cell_type.lower())

        for lr_key, lr_data in LIGAND_RECEPTOR_PAIRS.items():
            source_types = set(ct.lower() for ct in lr_data.get("source_cell_types", []))
            target_types = set(ct.lower() for ct in lr_data.get("target_cell_types", []))
            if query_types & (source_types | target_types):
                kb_interactions.append({
                    "ligand": lr_data.get("ligand", lr_key),
                    "receptor": lr_data.get("receptor", ""),
                    "pathway": lr_data.get("pathway", ""),
                    "source_cell_types": lr_data.get("source_cell_types", []),
                    "target_cell_types": lr_data.get("target_cell_types", []),
                    "source": "knowledge_base",
                })

    return LigandReceptorResponse(
        interactions=kb_interactions[:15],
        num_significant=len(kb_interactions),
        recommendations=["Connect RAG engine for full ligand-receptor analysis."],
    )


@router.post("/biomarker", response_model=BiomarkerResponse)
async def biomarker_discovery(request: BiomarkerRequest, req: Request):
    """Discover biomarker candidates from differential expression."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.discover_biomarkers(
                differential_expression=request.differential_expression,
                cell_type=request.cell_type,
                min_log2fc=request.min_log2fc,
                max_pval=request.max_pval,
            )
            return BiomarkerResponse(**result)
        except Exception as exc:
            logger.error(f"Biomarker discovery failed: {exc}")

    # Knowledge base fallback: use MARKER_GENE_DATABASE and IMMUNE_SIGNATURES
    kb_biomarkers = []
    kb_top_markers = []
    if request.cell_type:
        ct_key = request.cell_type.lower().replace(" ", "_")
        for mk_key, mk_data in MARKER_GENE_DATABASE.items():
            if ct_key in mk_key.lower():
                for gene in mk_data.get("markers", []):
                    kb_biomarkers.append({
                        "gene": gene,
                        "cell_type": mk_key,
                        "source": "marker_gene_database",
                    })
                    kb_top_markers.append(gene)
        # Also check immune signatures
        for sig_key, sig_data in IMMUNE_SIGNATURES.items():
            sig_genes = sig_data.get("genes", [])
            if any(ct_key in g.lower() for g in sig_data.get("cell_types", [])):
                for gene in sig_genes[:3]:
                    kb_biomarkers.append({
                        "gene": gene,
                        "signature": sig_key,
                        "source": "immune_signatures",
                    })
                    kb_top_markers.append(gene)

    return BiomarkerResponse(
        biomarkers=kb_biomarkers[:15],
        num_candidates=len(kb_biomarkers),
        top_markers=list(dict.fromkeys(kb_top_markers))[:10],
        recommendations=["Connect RAG engine for full biomarker discovery."],
    )


@router.post("/cart-validate", response_model=CARTValidateResponse)
async def cart_validate(request: CARTValidateRequest, req: Request):
    """Validate a CAR-T therapy target using single-cell expression data."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.validate_cart_target(
                target_gene=request.target_gene,
                tumor_type=request.tumor_type,
                expression_data=request.expression_data,
                tme_data=request.tme_data,
            )
            return CARTValidateResponse(**result)
        except Exception as exc:
            logger.error(f"CAR-T validation failed: {exc}")

    prompt = (
        f"Evaluate CAR-T target {request.target_gene} for "
        f"{request.tumor_type or 'solid tumor'}.\n"
        f"Assess on-tumor expression, off-tumor toxicity risk, TME "
        f"compatibility, and antigen escape risk."
    )
    answer = await _llm_fallback(req, prompt)

    # Knowledge base fallback: use CELL_TYPE_ATLAS and CANCER_TME_ATLAS
    kb_evidence = []
    kb_recs = [answer] if answer and "unavailable" not in answer.lower() else []
    off_tumor_risk = {}
    for ct_name, ct_data in CELL_TYPE_ATLAS.items():
        markers = [m.upper() for m in ct_data.get("markers", [])]
        if request.target_gene.upper() in markers:
            kb_evidence.append({
                "cell_type": ct_name,
                "expresses_target": True,
                "source": "cell_type_atlas",
            })
            tissues = ct_data.get("tissues", [])
            if tissues:
                off_tumor_risk[ct_name] = tissues
    if request.tumor_type:
        tumor_key = request.tumor_type.lower().replace(" ", "_")
        for atlas_key, atlas_data in CANCER_TME_ATLAS.items():
            if tumor_key in atlas_key.lower():
                kb_evidence.append({"source": "cancer_tme_atlas", "entry": atlas_key})
                break

    return CARTValidateResponse(
        target_gene=request.target_gene,
        on_tumor_pct=0.0,
        off_tumor_risk=off_tumor_risk,
        tme_compatibility=0.0,
        escape_risk="unknown",
        therapeutic_index=0.0,
        evidence=kb_evidence[:10],
        recommendations=kb_recs,
    )


@router.post("/treatment-monitor", response_model=TreatmentMonitorResponse)
async def treatment_monitor(request: TreatmentMonitorRequest, req: Request):
    """Monitor treatment response using longitudinal single-cell data."""
    _increment_metric(req, "workflow_requests_total")

    engine = _get_engine(req)
    if engine:
        try:
            result = engine.monitor_treatment(
                timepoints=request.timepoints,
                treatment=request.treatment,
                baseline_composition=request.baseline_composition,
                current_composition=request.current_composition,
            )
            return TreatmentMonitorResponse(**result)
        except Exception as exc:
            logger.error(f"Treatment monitoring failed: {exc}")

    # Knowledge base fallback: use IMMUNE_SIGNATURES for monitoring context
    kb_resistance = []
    kb_immune_shift = {}
    if request.treatment:
        treatment_lower = request.treatment.lower()
        for drug_key, drug_data in DRUG_SENSITIVITY_DATABASE.items():
            if drug_key.lower() in treatment_lower or treatment_lower in drug_key.lower():
                kb_resistance.extend([
                    {"mechanism": m, "source": "drug_sensitivity_database"}
                    for m in drug_data.get("resistance_mechanisms", [])[:3]
                ])
    if request.baseline_composition and request.current_composition:
        for ct in request.baseline_composition:
            baseline_val = request.baseline_composition.get(ct, 0)
            current_val = request.current_composition.get(ct, 0)
            if baseline_val > 0:
                kb_immune_shift[ct] = {
                    "baseline": baseline_val,
                    "current": current_val,
                    "change_pct": round((current_val - baseline_val) / baseline_val * 100, 1),
                }

    return TreatmentMonitorResponse(
        response_assessment="insufficient_data",
        resistance_indicators=kb_resistance[:5],
        immune_shift=kb_immune_shift,
        recommendations=["Connect RAG engine for full treatment monitoring."],
    )


@router.post("/workflow/{workflow_type}")
async def run_workflow(workflow_type: str, request: WorkflowRequest, req: Request):
    """Generic workflow dispatch for any single-cell analysis workflow."""
    _increment_metric(req, "workflow_requests_total")

    valid_workflows = [
        "cell_type_annotation", "tme_profiling", "drug_response",
        "subclonal_analysis", "spatial_niche", "trajectory_inference",
        "ligand_receptor", "biomarker_discovery", "cart_validation",
        "treatment_monitoring",
    ]
    if workflow_type not in valid_workflows:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown workflow: '{workflow_type}'. Valid: {valid_workflows}",
        )

    wf_engine = _get_workflow_engine(req)
    if wf_engine:
        try:
            data = request.data or {}
            if request.question:
                data["question"] = request.question
            result = await wf_engine.execute(workflow_type, data)
            return result
        except Exception as exc:
            logger.error(f"Workflow '{workflow_type}' failed: {exc}")
            raise HTTPException(status_code=500, detail="Internal processing error")

    return {
        "workflow_type": workflow_type,
        "status": "completed",
        "result": f"Workflow '{workflow_type}' executed (engine unavailable).",
        "note": "Connect RAG engine for full workflow execution.",
    }


# =====================================================================
# Reference Catalogues
# =====================================================================

@router.get("/cell-types")
async def list_cell_types():
    """Return the cell type catalogue with canonical markers."""
    from src.knowledge import CELL_TYPE_ATLAS
    cell_types = []
    for ct_name, ct_data in CELL_TYPE_ATLAS.items():
        cell_types.append({
            "name": ct_name.replace("_", " "),
            "compartment": ct_data.get("cell_ontology_id", ""),
            "canonical_markers": ct_data.get("markers", []),
            "description": ct_data.get("description", ""),
            "tissues": ct_data.get("tissues", []),
            "subtypes": ct_data.get("subtypes", []),
        })
    return {"cell_types": cell_types, "total": len(cell_types)}


@router.get("/markers")
async def list_markers():
    """Return marker gene reference per cell type."""
    from src.knowledge import MARKER_GENE_DATABASE
    return {"markers": MARKER_GENE_DATABASE, "total": len(MARKER_GENE_DATABASE)}


@router.get("/tme-classes")
async def list_tme_classes():
    """Return TME classification reference."""
    from src.knowledge import TME_PROFILES
    return {"tme_classes": TME_PROFILES, "total": len(TME_PROFILES)}


@router.get("/spatial-platforms")
async def list_spatial_platforms():
    """Return spatial transcriptomics platform reference."""
    from src.knowledge import SPATIAL_PLATFORMS
    return {"platforms": SPATIAL_PLATFORMS, "total": len(SPATIAL_PLATFORMS)}


@router.get("/knowledge-version")
async def knowledge_version():
    """Return knowledge base version metadata."""
    return KNOWLEDGE_VERSION
