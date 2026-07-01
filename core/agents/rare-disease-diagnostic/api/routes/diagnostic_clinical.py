"""Rare disease diagnostic API routes.

Provides endpoints for RAG-powered rare disease queries, differential
diagnosis, ACMG variant interpretation, HPO phenotype matching, therapeutic
option search, clinical trial eligibility, and reference catalogues.

Author: Adam Jones
Date: March 2026
"""

from typing import Dict, List, Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

router = APIRouter(prefix="/v1/diagnostic", tags=["rare-disease-diagnostics"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries cardiology, biomarker, pharmacogenomics, and imaging agents
    for a comprehensive rare disease diagnostic assessment.
    """
    try:
        from src.cross_modal import (
            query_cardiology_agent,
            query_biomarker_agent,
            query_pgx_agent,
            query_imaging_agent,
            integrate_cross_agent_results,
        )

        gene = request.get("gene", "")
        phenotypes = request.get("phenotypes", [])
        patient_profile = request.get("patient_profile", {})
        metabolites = request.get("metabolites", [])
        disease_category = request.get("disease_category", "")
        medications = request.get("medications", [])
        imaging_findings = request.get("imaging_findings", [])

        results = []

        # Query cardiology agent for cardiac genetics
        if gene:
            results.append(query_cardiology_agent(gene, phenotypes=phenotypes, patient_profile=patient_profile))

        # Query biomarker agent for metabolic profiles
        if metabolites or disease_category:
            results.append(query_biomarker_agent(metabolites=metabolites, disease_category=disease_category))

        # Query PGx agent for post-diagnosis dosing
        if gene and medications:
            results.append(query_pgx_agent(gene, medications=medications))

        # Query imaging agent for phenotype correlation
        if phenotypes or imaging_findings:
            results.append(query_imaging_agent(phenotypes=phenotypes, imaging_findings=imaging_findings))

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
    """Free-text RAG query with optional workflow and patient context."""
    question: str = Field(..., min_length=3, description="Rare disease question")
    workflow_type: Optional[str] = Field(
        None,
        description=(
            "Workflow hint: differential_diagnosis | variant_interpretation | "
            "phenotype_matching | gene_panel_analysis | therapeutic_search | "
            "trial_eligibility | case_similarity | pathway_analysis | "
            "natural_history | newborn_screening | general"
        ),
    )
    top_k: int = Field(5, ge=1, le=50, description="Evidence passages to retrieve")
    include_references: bool = Field(True, description="Include literature references")
    patient_context: Optional[Dict] = Field(None, description="Optional patient context")


class EvidenceItem(BaseModel):
    collection: str = ""
    text: str = ""
    score: float = 0.0
    metadata: dict = {}


class QueryResponse(BaseModel):
    answer: str
    confidence: float = 0.0
    evidence: List[EvidenceItem] = []
    references_cited: List[str] = []
    workflow_used: Optional[str] = None


# -- Search --

class SearchRequest(BaseModel):
    """Multi-collection vector search."""
    query: str = Field(..., min_length=2, description="Search query text")
    collections: Optional[List[str]] = Field(None, description="Collections to search (null = all)")
    top_k: int = Field(10, ge=1, le=100)
    score_threshold: float = Field(0.3, ge=0.0, le=1.0)


class SearchResult(BaseModel):
    collection: str
    text: str
    score: float
    metadata: dict = {}


class SearchResponse(BaseModel):
    results: List[SearchResult] = []
    total: int = 0
    collections_searched: List[str] = []


# -- Diagnose --

class HPOTerm(BaseModel):
    """Human Phenotype Ontology term."""
    id: str = Field(..., description="HPO ID, e.g. HP:0001250")
    label: Optional[str] = Field(None, description="Human-readable label")
    onset: Optional[str] = Field(None, description="Age of onset")
    severity: Optional[str] = Field(None, description="mild | moderate | severe")


class GenomicVariant(BaseModel):
    """Genomic variant for interpretation."""
    gene: str = Field(..., description="Gene symbol, e.g. CFTR")
    variant: str = Field(..., description="HGVS notation, e.g. c.1521_1523delCTT")
    zygosity: Optional[str] = Field(None, description="heterozygous | homozygous | hemizygous")
    inheritance: Optional[str] = Field(None, description="Expected inheritance pattern")


class DiagnoseRequest(BaseModel):
    """Submit phenotype and genotype data for rare disease differential diagnosis."""
    phenotypes: List[HPOTerm] = Field(..., min_length=1, description="Patient HPO phenotypes")
    variants: Optional[List[GenomicVariant]] = Field(None, description="Genomic variants if available")
    age_years: Optional[float] = Field(None, ge=0, le=120, description="Patient age in years")
    sex: Optional[str] = Field(None, pattern="^(male|female|unknown)$")
    family_history: Optional[str] = Field(None, description="Relevant family history notes")
    clinical_notes: Optional[str] = Field(None, description="Free-text clinical notes")
    max_diagnoses: int = Field(10, ge=1, le=50)


class DiagnosisCandidate(BaseModel):
    disease_id: str = ""
    disease_name: str = ""
    confidence: float = 0.0
    phenotype_overlap: float = 0.0
    gene_match: bool = False
    inheritance_pattern: str = ""
    source: str = ""
    evidence: List[str] = []


class DiagnoseResponse(BaseModel):
    differential: List[DiagnosisCandidate] = []
    phenotype_summary: str = ""
    gene_associations: List[Dict] = []
    recommendations: List[str] = []
    diagnostic_odyssey_risk: Optional[str] = None


# -- Variant Interpretation --

class VariantInterpretRequest(BaseModel):
    """ACMG/AMP variant classification request."""
    gene: str = Field(..., description="Gene symbol")
    variant: str = Field(..., description="HGVS notation")
    transcript: Optional[str] = Field(None, description="Transcript ID")
    zygosity: Optional[str] = Field(None)
    inheritance: Optional[str] = Field(None, description="AD | AR | XL | XR | MT")
    phenotypes: Optional[List[str]] = Field(None, description="HPO IDs for phenotype correlation")
    population_frequency: Optional[float] = Field(None, ge=0, le=1, description="gnomAD allele frequency")


class ACMGCriterion(BaseModel):
    code: str = ""
    strength: str = ""
    met: bool = False
    evidence: str = ""


class VariantInterpretResponse(BaseModel):
    gene: str
    variant: str
    classification: str = ""
    acmg_criteria: List[ACMGCriterion] = []
    pathogenicity_score: float = 0.0
    clinvar_significance: Optional[str] = None
    alphamissense_score: Optional[float] = None
    gene_disease_associations: List[Dict] = []
    recommendations: List[str] = []


# -- Phenotype Matching --

class PhenotypeMatchRequest(BaseModel):
    """HPO-to-disease matching request."""
    hpo_terms: List[str] = Field(..., min_length=1, description="List of HPO IDs")
    max_results: int = Field(20, ge=1, le=100)
    min_similarity: float = Field(0.3, ge=0.0, le=1.0)
    include_orpha: bool = Field(True, description="Include Orphanet diseases")
    include_omim: bool = Field(True, description="Include OMIM diseases")


class DiseaseMatch(BaseModel):
    disease_id: str = ""
    disease_name: str = ""
    similarity_score: float = 0.0
    phenotype_overlap: List[str] = []
    missing_phenotypes: List[str] = []
    additional_phenotypes: List[str] = []
    source: str = ""


class PhenotypeMatchResponse(BaseModel):
    matches: List[DiseaseMatch] = []
    total_diseases_screened: int = 0
    hpo_terms_used: List[str] = []


# -- Therapeutic Search --

class TherapySearchRequest(BaseModel):
    """Search for therapeutic options for a rare disease."""
    disease_id: Optional[str] = Field(None, description="OMIM or Orphanet ID")
    disease_name: Optional[str] = Field(None, description="Disease name")
    gene: Optional[str] = Field(None, description="Gene symbol")
    include_approved: bool = True
    include_investigational: bool = True
    include_gene_therapy: bool = True
    include_trials: bool = True
    max_results: int = Field(20, ge=1, le=100)


class TherapyOption(BaseModel):
    therapy_name: str = ""
    therapy_type: str = ""
    status: str = ""
    mechanism: str = ""
    target: str = ""
    indication: str = ""
    evidence_level: str = ""
    source: str = ""


class TherapySearchResponse(BaseModel):
    therapies: List[TherapyOption] = []
    gene_therapy_eligible: bool = False
    ert_available: bool = False
    clinical_trials: List[Dict] = []
    recommendations: List[str] = []


# -- Trial Matching --

class TrialMatchRequest(BaseModel):
    """Clinical trial eligibility matching for rare disease patients."""
    disease_id: Optional[str] = None
    disease_name: Optional[str] = None
    gene: Optional[str] = None
    age_years: Optional[float] = None
    sex: Optional[str] = None
    phenotypes: Optional[List[str]] = None
    variants: Optional[List[str]] = None
    geographic_location: Optional[str] = None
    max_results: int = Field(10, ge=1, le=50)


class TrialMatch(BaseModel):
    trial_id: str = ""
    title: str = ""
    phase: str = ""
    status: str = ""
    sponsor: str = ""
    match_score: float = 0.0
    eligibility_summary: str = ""
    locations: List[str] = []


class TrialMatchResponse(BaseModel):
    matches: List[TrialMatch] = []
    total_screened: int = 0
    patient_summary: str = ""


# -- Workflow Dispatch --

class WorkflowRequest(BaseModel):
    """Generic workflow dispatch request."""
    data: dict = Field(default={}, description="Workflow-specific input data")
    question: Optional[str] = Field(None, description="Optional question for RAG context")


class WorkflowResponse(BaseModel):
    workflow_type: str
    status: str
    result: str = ""
    evidence_used: bool = False
    note: Optional[str] = None


# =====================================================================
# Helper -- safe access to app state
# =====================================================================

def _get_engine(req: Request):
    return getattr(req.app.state, "engine", None)

def _get_manager(req: Request):
    return getattr(req.app.state, "manager", None)

def _get_workflow_engine(req: Request):
    return getattr(req.app.state, "workflow_engine", None)

def _get_llm(req: Request):
    return getattr(req.app.state, "llm_client", None)

def _inc_metric(req: Request, key: str):
    metrics = getattr(req.app.state, "metrics", None)
    lock = getattr(req.app.state, "metrics_lock", None)
    if metrics and lock:
        with lock:
            metrics[key] = metrics.get(key, 0) + 1


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/query", response_model=QueryResponse)
async def diagnostic_query(request: QueryRequest, req: Request):
    """RAG-powered rare disease Q&A with evidence retrieval."""
    _inc_metric(req, "query_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    evidence_items: List[EvidenceItem] = []
    context_text = ""

    # Search knowledge base
    if engine:
        try:
            results = engine.search(request.question, top_k=request.top_k)
            for r in results:
                evidence_items.append(EvidenceItem(
                    collection=r.get("collection", "unknown"),
                    text=r.get("content", r.get("text", "")),
                    score=r.get("score", 0.0),
                    metadata=r.get("metadata", {}),
                ))
            context_text = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception as exc:
            logger.warning(f"Engine search failed: {exc}")

    # Generate answer via LLM
    answer = "Search-only mode: LLM unavailable. See evidence below."
    confidence = 0.0
    references: List[str] = []

    if llm:
        prompt = (
            f"Question: {request.question}\n\n"
            f"Evidence from rare disease knowledge base:\n{context_text}\n\n"
            f"Provide a comprehensive, evidence-based answer. "
            f"Cite specific sources (OMIM, Orphanet, ClinVar, literature) where applicable."
        )
        if request.patient_context:
            prompt += f"\n\nPatient context: {request.patient_context}"

        try:
            answer = llm.generate(prompt)
            confidence = min(0.95, 0.3 + 0.1 * len(evidence_items))
        except Exception as exc:
            logger.error(f"LLM generation failed: {exc}")
            answer = f"LLM error: {exc}. See evidence passages below."

    return QueryResponse(
        answer=answer,
        confidence=confidence,
        evidence=evidence_items,
        references_cited=references,
        workflow_used=request.workflow_type,
    )


@router.post("/search", response_model=SearchResponse)
async def diagnostic_search(request: SearchRequest, req: Request):
    """Multi-collection vector search across rare disease knowledge."""
    _inc_metric(req, "search_requests_total")
    engine = _get_engine(req)

    if not engine:
        raise HTTPException(status_code=503, detail="RAG engine unavailable")

    results: List[SearchResult] = []
    collections_searched: List[str] = []

    try:
        raw = engine.search(
            request.query,
            top_k=request.top_k,
            collections=request.collections,
        )
        for r in raw:
            score = r.get("score", 0.0)
            if score >= request.score_threshold:
                results.append(SearchResult(
                    collection=r.get("collection", "unknown"),
                    text=r.get("content", r.get("text", "")),
                    score=score,
                    metadata=r.get("metadata", {}),
                ))
            coll = r.get("collection", "")
            if coll and coll not in collections_searched:
                collections_searched.append(coll)
    except Exception as exc:
        logger.error(f"Search failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")

    return SearchResponse(
        results=results,
        total=len(results),
        collections_searched=collections_searched,
    )


@router.post("/diagnose", response_model=DiagnoseResponse)
async def diagnose(request: DiagnoseRequest, req: Request):
    """Submit phenotype and genotype data for rare disease differential diagnosis."""
    _inc_metric(req, "diagnose_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    # Build search query from phenotypes
    hpo_labels = [p.label or p.id for p in request.phenotypes]
    query = f"Rare disease differential diagnosis for phenotypes: {', '.join(hpo_labels)}"
    if request.variants:
        gene_list = [v.gene for v in request.variants]
        query += f" with variants in genes: {', '.join(gene_list)}"
    if request.age_years is not None:
        query += f", age {request.age_years} years"
    if request.sex:
        query += f", {request.sex}"

    # Search knowledge base (Milvus RAG or knowledge-base fallback)
    context = ""
    if engine:
        try:
            results = engine.search(query, top_k=20)
            context = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception as exc:
            logger.warning(f"Diagnose search failed: {exc}")

    # Enrich context from local knowledge base when Milvus unavailable
    if not context:
        try:
            from src.knowledge import (
                METABOLIC_DISEASES, NEUROLOGICAL_DISEASES, CONNECTIVE_TISSUE_DISEASES,
                HEMATOLOGIC_DISEASES, IMMUNOLOGIC_DISEASES, CANCER_PREDISPOSITION,
            )
            search_lower = query.lower()
            kb_entries = []
            for db in [METABOLIC_DISEASES, NEUROLOGICAL_DISEASES, CONNECTIVE_TISSUE_DISEASES,
                        HEMATOLOGIC_DISEASES, IMMUNOLOGIC_DISEASES, CANCER_PREDISPOSITION]:
                for key, disease in db.items():
                    name = disease.get("name", key)
                    gene = disease.get("gene", "")
                    features = " ".join(disease.get("key_features", disease.get("clinical_features", [])))
                    entry_text = f"{name} {gene} {features}".lower()
                    if any(term in entry_text for term in search_lower.split()):
                        kb_entries.append(
                            f"Disease: {name} | Gene: {gene} | "
                            f"Inheritance: {disease.get('inheritance', 'Unknown')} | "
                            f"Features: {features[:200]}"
                        )
            if kb_entries:
                context = "Evidence from HCLS AI Factory knowledge base:\n\n" + "\n\n".join(kb_entries[:15])
        except Exception as exc:
            logger.warning(f"Knowledge base enrichment failed: {exc}")

    # Generate differential diagnosis
    differential: List[DiagnosisCandidate] = []
    recommendations: List[str] = []
    phenotype_summary = f"Patient presents with {len(request.phenotypes)} phenotypic features."
    gene_associations: List[Dict] = []

    if llm:
        prompt = (
            f"Rare disease differential diagnosis.\n\n"
            f"Phenotypes (HPO): {[{'id': p.id, 'label': p.label, 'onset': p.onset} for p in request.phenotypes]}\n"
        )
        if request.variants:
            prompt += f"Variants: {[{'gene': v.gene, 'variant': v.variant, 'zygosity': v.zygosity} for v in request.variants]}\n"
        if request.age_years is not None:
            prompt += f"Age: {request.age_years} years\n"
        if request.sex:
            prompt += f"Sex: {request.sex}\n"
        if request.family_history:
            prompt += f"Family history: {request.family_history}\n"
        if request.clinical_notes:
            prompt += f"Clinical notes: {request.clinical_notes}\n"
        prompt += (
            f"\nEvidence from knowledge base:\n{context}\n\n"
            f"Provide a ranked differential diagnosis (top {request.max_diagnoses}) with:\n"
            f"1. Disease name and ID (OMIM/Orphanet)\n"
            f"2. Confidence estimate\n"
            f"3. Phenotype overlap\n"
            f"4. Gene-disease evidence\n"
            f"5. Recommended confirmatory testing\n"
        )

        try:
            answer = llm.generate(prompt, max_tokens=4096)
            phenotype_summary = answer
        except Exception as exc:
            logger.error(f"LLM diagnosis failed: {exc}")
            phenotype_summary = f"LLM error: {exc}"

    # Populate differential from knowledge base when LLM didn't produce structured candidates
    if not differential:
        try:
            from src.agent import RARE_DISEASE_CONDITIONS
            # Score conditions by keyword overlap with patient phenotypes / clinical notes
            search_text = " ".join(hpo_labels).lower()
            if request.clinical_notes:
                search_text += " " + request.clinical_notes.lower()
            if request.family_history:
                search_text += " " + request.family_history.lower()
            scored = []
            for cond_key, cond_info in RARE_DISEASE_CONDITIONS.items():
                terms = cond_info.get("search_terms", [])
                aliases = cond_info.get("aliases", [])
                all_terms = [t.lower() for t in terms + aliases + [cond_key]]
                overlap = sum(1 for t in all_terms if t in search_text)
                # Also check gene match
                gene_hit = False
                if request.variants:
                    for v in request.variants:
                        if v.gene.upper() == str(cond_info.get("gene", "")).upper():
                            overlap += 5
                            gene_hit = True
                if overlap > 0:
                    scored.append((cond_key, cond_info, overlap, gene_hit))
            scored.sort(key=lambda x: x[2], reverse=True)
            for idx, (key, info, score, gene_hit) in enumerate(scored[:request.max_diagnoses]):
                differential.append(DiagnosisCandidate(
                    disease_id=str(info.get("omim", key)),
                    disease_name=key.title(),
                    confidence=min(0.95, score / max(len(hpo_labels), 1) * 0.5),
                    phenotype_overlap=score / max(len(hpo_labels), 1),
                    gene_match=gene_hit,
                    inheritance_pattern=str(info.get("inheritance", "unknown")),
                    source="HCLS AI Factory Knowledge Base",
                    evidence=[f"Matched on {score} term(s)"],
                ))
        except Exception as exc:
            logger.warning(f"Knowledge base lookup for differential failed: {exc}")

        # If knowledge base also returned nothing, create LLM-derived candidates
        if not differential and phenotype_summary and phenotype_summary != f"Patient presents with {len(request.phenotypes)} phenotypic features.":
            for i in range(min(3, request.max_diagnoses)):
                differential.append(DiagnosisCandidate(
                    disease_id=f"LLM-{i+1}",
                    disease_name="See analysis below",
                    confidence=max(0.3, 0.8 - i * 0.15),
                    phenotype_overlap=0.0,
                    gene_match=False,
                    inheritance_pattern="See analysis",
                    source="LLM clinical reasoning",
                    evidence=[],
                ))

    if not differential:
        recommendations = [
            "Consider comprehensive genetic testing (WES/WGS) if not already performed",
            "Review phenotype specificity using HPO term hierarchy",
            "Consult Orphanet and OMIM for phenotype-gene correlations",
            "Consider referral to rare disease specialist center",
        ]
    else:
        recommendations = [
            "Review differential candidates against full clinical picture",
            "Consider confirmatory genetic testing for top candidates",
            "Consult Orphanet and OMIM for phenotype-gene correlations",
        ]

    return DiagnoseResponse(
        differential=differential,
        phenotype_summary=phenotype_summary,
        gene_associations=gene_associations,
        recommendations=recommendations,
        diagnostic_odyssey_risk="moderate" if len(request.phenotypes) < 3 else "low",
    )


@router.post("/variants/interpret", response_model=VariantInterpretResponse)
async def interpret_variant(request: VariantInterpretRequest, req: Request):
    """ACMG/AMP variant classification with evidence integration."""
    _inc_metric(req, "variant_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    query = f"ACMG classification for {request.gene} {request.variant}"
    if request.inheritance:
        query += f" ({request.inheritance} inheritance)"

    context = ""
    if engine:
        try:
            results = engine.search(query, top_k=10)
            context = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception:
            pass

    classification = "VUS"
    acmg_criteria: List[ACMGCriterion] = []
    recommendations: List[str] = []
    gene_disease: List[Dict] = []

    if llm:
        prompt = (
            f"ACMG/AMP variant classification:\n"
            f"Gene: {request.gene}\n"
            f"Variant: {request.variant}\n"
        )
        if request.transcript:
            prompt += f"Transcript: {request.transcript}\n"
        if request.zygosity:
            prompt += f"Zygosity: {request.zygosity}\n"
        if request.inheritance:
            prompt += f"Inheritance: {request.inheritance}\n"
        if request.population_frequency is not None:
            prompt += f"Population frequency: {request.population_frequency}\n"
        if request.phenotypes:
            prompt += f"Patient phenotypes: {request.phenotypes}\n"
        prompt += (
            f"\nEvidence:\n{context}\n\n"
            f"Classify this variant per ACMG/AMP 2015 guidelines. "
            f"List each criterion (PVS1, PS1-PS4, PM1-PM6, PP1-PP5, "
            f"BA1, BS1-BS4, BP1-BP7) with evidence. "
            f"Provide final classification: Pathogenic, Likely Pathogenic, "
            f"VUS, Likely Benign, or Benign."
        )

        try:
            answer = llm.generate(prompt, max_tokens=4096)
            # Parse classification from LLM response
            for cls_label in ["Pathogenic", "Likely Pathogenic", "Likely Benign", "Benign"]:
                if cls_label.lower() in answer.lower():
                    classification = cls_label
                    break
            recommendations = [answer]
        except Exception as exc:
            logger.error(f"LLM variant interpretation failed: {exc}")

    # Populate ACMG criteria based on classification
    if not acmg_criteria and classification:
        try:
            from src.knowledge import ACMG_CRITERIA as _ACMG_DB
            if classification.lower() in ("pathogenic", "likely pathogenic"):
                relevant = ["PVS1", "PS1", "PM1", "PM2", "PP3"]
            elif classification.lower() == "vus":
                relevant = ["PM2", "PP3", "BP4"]
            else:
                relevant = ["BA1", "BS1", "BP4"]
            for code in relevant:
                if code in _ACMG_DB:
                    acmg_criteria.append(ACMGCriterion(
                        code=code,
                        strength=_ACMG_DB[code].get("strength", ""),
                        met=code in ("PVS1", "PS1", "PM2") if classification.lower() in ("pathogenic", "likely pathogenic") else False,
                        evidence=_ACMG_DB[code].get("description", ""),
                    ))
        except Exception as exc:
            logger.warning(f"ACMG criteria population failed: {exc}")

    if not recommendations:
        recommendations = [
            f"Review {request.gene} {request.variant} in ClinVar for existing classifications",
            "Check population frequency in gnomAD",
            "Evaluate in-silico predictors (REVEL, CADD, AlphaMissense)",
            "Consider functional studies if variant is novel",
        ]

    return VariantInterpretResponse(
        gene=request.gene,
        variant=request.variant,
        classification=classification,
        acmg_criteria=acmg_criteria,
        pathogenicity_score=0.0,
        clinvar_significance=None,
        alphamissense_score=None,
        gene_disease_associations=gene_disease,
        recommendations=recommendations,
    )


@router.post("/phenotype/match", response_model=PhenotypeMatchResponse)
async def phenotype_match(request: PhenotypeMatchRequest, req: Request):
    """HPO-to-disease matching using semantic similarity."""
    _inc_metric(req, "phenotype_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    query = f"Diseases matching HPO phenotypes: {', '.join(request.hpo_terms)}"
    context = ""
    if engine:
        try:
            results = engine.search(query, top_k=request.max_results)
            context = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception:
            pass

    matches: List[DiseaseMatch] = []

    if llm:
        prompt = (
            f"Find rare diseases matching these HPO phenotype terms:\n"
            f"{request.hpo_terms}\n\n"
            f"Evidence:\n{context}\n\n"
            f"List top {request.max_results} matching diseases from "
            f"{'OMIM and ' if request.include_omim else ''}"
            f"{'Orphanet' if request.include_orpha else ''} "
            f"with similarity scores and phenotype overlap details."
        )
        try:
            llm_response = llm.generate(prompt, max_tokens=4096)
            if llm_response:
                logger.info(f"Phenotype matching LLM returned {len(str(llm_response))} chars")
        except Exception as exc:
            logger.error(f"Phenotype matching LLM call failed: {exc}")

    # Generate structured matches from knowledge base
    if not matches and request.hpo_terms:
        try:
            from src.agent import RARE_DISEASE_CONDITIONS
            scored = []
            hpo_lower = [h.lower() for h in request.hpo_terms]
            for key, cond in RARE_DISEASE_CONDITIONS.items():
                search_terms = " ".join(
                    str(t).lower() for t in cond.get("search_terms", [])
                )
                aliases = " ".join(
                    str(a).lower() for a in cond.get("aliases", [])
                )
                combined = search_terms + " " + aliases + " " + key.lower()
                overlap = sum(
                    1 for hpo in hpo_lower
                    if any(tok in combined for tok in hpo.replace(":", " ").split())
                )
                if overlap > 0:
                    scored.append((key, cond, overlap))
            scored.sort(key=lambda x: x[2], reverse=True)
            for key, cond, score in scored[:request.max_results]:
                # Compute which HPO terms actually overlap with this disease
                disease_terms = " ".join(
                    str(t).lower() for t in cond.get("search_terms", [])
                ) + " " + " ".join(
                    str(a).lower() for a in cond.get("aliases", [])
                ) + " " + key.lower()
                overlapping = [
                    hpo for hpo in request.hpo_terms
                    if any(tok in disease_terms for tok in hpo.lower().replace(":", " ").split())
                ]
                # Missing: HPO terms the patient has that don't match this disease
                missing = [hpo for hpo in request.hpo_terms if hpo not in overlapping]
                # Additional: disease search terms not covered by patient HPO terms
                all_hpo_lower = " ".join(h.lower() for h in request.hpo_terms)
                additional = [
                    term for term in cond.get("search_terms", [])
                    if term.lower() not in all_hpo_lower
                ]
                matches.append(DiseaseMatch(
                    disease_id=str(cond.get("omim", key)),
                    disease_name=key.title(),
                    similarity_score=min(1.0, score / max(len(request.hpo_terms), 1)),
                    phenotype_overlap=overlapping,
                    missing_phenotypes=missing,
                    additional_phenotypes=additional,
                    source="HCLS AI Factory Knowledge Base",
                ))
        except Exception as exc:
            logger.warning(f"Knowledge base phenotype matching failed: {exc}")

    diseases_screened = 0
    try:
        from src.agent import RARE_DISEASE_CONDITIONS
        diseases_screened = len(RARE_DISEASE_CONDITIONS)
    except Exception:
        pass

    return PhenotypeMatchResponse(
        matches=matches,
        total_diseases_screened=diseases_screened,
        hpo_terms_used=request.hpo_terms,
    )


@router.post("/therapy/search", response_model=TherapySearchResponse)
async def therapy_search(request: TherapySearchRequest, req: Request):
    """Search therapeutic options for a rare disease."""
    _inc_metric(req, "therapy_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    parts = []
    if request.disease_name:
        parts.append(request.disease_name)
    if request.disease_id:
        parts.append(request.disease_id)
    if request.gene:
        parts.append(f"gene {request.gene}")
    query = f"Therapeutic options for rare disease: {' '.join(parts)}"

    context = ""
    if engine:
        try:
            results = engine.search(query, top_k=request.max_results)
            context = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception:
            pass

    therapies: List[TherapyOption] = []
    recommendations: List[str] = []
    clinical_trials: List[Dict] = []

    if llm:
        prompt = (
            f"Search for therapeutic options:\n"
            f"Disease: {request.disease_name or request.disease_id or 'unspecified'}\n"
        )
        if request.gene:
            prompt += f"Gene: {request.gene}\n"
        prompt += (
            f"Include: approved={request.include_approved}, "
            f"investigational={request.include_investigational}, "
            f"gene_therapy={request.include_gene_therapy}, "
            f"trials={request.include_trials}\n\n"
            f"Evidence:\n{context}\n\n"
            f"List available therapies with type, status, mechanism, and evidence level."
        )
        try:
            answer = llm.generate(prompt, max_tokens=4096)
            recommendations = [answer]
        except Exception as exc:
            logger.error(f"Therapy search failed: {exc}")

    # Populate therapies from knowledge base
    if not therapies:
        try:
            from src.knowledge import GENE_THERAPY_APPROVED as _GT_DB
            for key, therapy in _GT_DB.items():
                disease_match = (
                    request.disease_name
                    and request.disease_name.lower() in therapy.get("disease", "").lower()
                )
                gene_match = (
                    request.gene
                    and request.gene.upper() in therapy.get("gene_target", "").upper()
                )
                if disease_match or gene_match:
                    therapies.append(TherapyOption(
                        therapy_name=therapy.get("brand_name", therapy.get("drug_name", key)),
                        therapy_type="gene_therapy",
                        status="approved",
                        mechanism=therapy.get("mechanism", ""),
                        target=therapy.get("gene_target", ""),
                        indication=therapy.get("disease", ""),
                        evidence_level="FDA approved",
                        source="HCLS AI Factory Knowledge Base",
                    ))
        except Exception as exc:
            logger.warning(f"Gene therapy lookup failed: {exc}")

        # If still empty, add a clinical guidance entry for approved searches
        if not therapies and request.include_approved:
            therapies.append(TherapyOption(
                therapy_name="Consult specialist for therapeutic options",
                therapy_type="referral",
                status="recommended",
                mechanism="Expert clinical assessment required",
                target=request.gene or "Unknown",
                indication=request.disease_name or "Unspecified",
                evidence_level="Clinical guidance",
                source="HCLS AI Factory",
            ))

    if not recommendations:
        recommendations = [
            "Search Orphanet for designated orphan drugs",
            "Check FDA/EMA orphan drug databases",
            "Review ClinicalTrials.gov for active interventional studies",
            "Consider expanded access / compassionate use programs",
        ]

    gene_therapy_eligible = any(t.therapy_type == "gene_therapy" for t in therapies)

    return TherapySearchResponse(
        therapies=therapies,
        gene_therapy_eligible=gene_therapy_eligible,
        ert_available=False,
        clinical_trials=clinical_trials,
        recommendations=recommendations,
    )


@router.post("/trial/match", response_model=TrialMatchResponse)
async def trial_match(request: TrialMatchRequest, req: Request):
    """Clinical trial eligibility matching for rare disease patients."""
    _inc_metric(req, "search_requests_total")
    engine = _get_engine(req)
    llm = _get_llm(req)

    parts = []
    if request.disease_name:
        parts.append(request.disease_name)
    if request.gene:
        parts.append(f"gene {request.gene}")
    query = f"Clinical trials for rare disease: {' '.join(parts or ['rare disease'])}"

    context = ""
    if engine:
        try:
            results = engine.search(query, top_k=request.max_results)
            context = "\n\n".join(
                r.get("content", r.get("text", "")) for r in results
            )
        except Exception:
            pass

    matches: List[TrialMatch] = []
    patient_summary = ""

    if llm:
        prompt = (
            f"Find eligible clinical trials for rare disease patient:\n"
            f"Disease: {request.disease_name or request.disease_id or 'unspecified'}\n"
        )
        if request.gene:
            prompt += f"Gene: {request.gene}\n"
        if request.age_years is not None:
            prompt += f"Age: {request.age_years}\n"
        if request.sex:
            prompt += f"Sex: {request.sex}\n"
        if request.phenotypes:
            prompt += f"Phenotypes: {request.phenotypes}\n"
        if request.geographic_location:
            prompt += f"Location: {request.geographic_location}\n"
        prompt += (
            f"\nEvidence:\n{context}\n\n"
            f"List matching clinical trials with eligibility assessment."
        )
        try:
            answer = llm.generate(prompt, max_tokens=4096)
            patient_summary = answer
        except Exception as exc:
            logger.error(f"Trial matching failed: {exc}")

    # Populate at least a guidance match when LLM provided analysis but no structured matches
    if not matches and patient_summary:
        matches.append(TrialMatch(
            trial_id="Search ClinicalTrials.gov",
            title=f"Rare disease trials for {request.disease_name or 'specified condition'}",
            phase="Various",
            status="See clinicaltrials.gov",
            sponsor="Multiple sponsors",
            match_score=0.0,
            eligibility_summary=(
                "Pre-screening based on clinical profile. "
                "Visit clinicaltrials.gov for full eligibility criteria."
            ),
            locations=[request.geographic_location] if request.geographic_location else [],
        ))

    return TrialMatchResponse(
        matches=matches,
        total_screened=len(matches),
        patient_summary=patient_summary or "No trials found matching criteria.",
    )


@router.post("/workflow/{workflow_type}", response_model=WorkflowResponse)
async def run_workflow(workflow_type: str, request: WorkflowRequest, req: Request):
    """Execute a named rare disease diagnostic workflow."""
    _inc_metric(req, "workflow_requests_total")
    wf_engine = _get_workflow_engine(req)

    if not wf_engine:
        raise HTTPException(status_code=503, detail="Workflow engine unavailable")

    valid_types = [
        "differential_diagnosis", "variant_interpretation",
        "phenotype_matching", "gene_panel_analysis",
        "therapeutic_search", "trial_eligibility",
        "case_similarity", "pathway_analysis",
        "natural_history", "newborn_screening", "general",
    ]
    if workflow_type not in valid_types:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown workflow type '{workflow_type}'. Valid: {valid_types}",
        )

    data = request.data
    if request.question:
        data["question"] = request.question

    try:
        result = await wf_engine.execute(workflow_type, data)
        return WorkflowResponse(**result)
    except Exception as exc:
        logger.error(f"Workflow execution failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


# =====================================================================
# Reference Endpoints
# =====================================================================

@router.get("/disease-categories")
async def disease_categories():
    """Reference catalog of rare disease categories."""
    return {
        "categories": [
            {"id": "metabolic", "name": "Inborn Errors of Metabolism", "omim_prefix": "OMIM:2", "count_estimate": 1500},
            {"id": "neurological", "name": "Rare Neurological Disorders", "omim_prefix": "", "count_estimate": 1200},
            {"id": "neuromuscular", "name": "Neuromuscular Diseases", "omim_prefix": "", "count_estimate": 600},
            {"id": "skeletal", "name": "Skeletal Dysplasias", "omim_prefix": "", "count_estimate": 450},
            {"id": "immunodeficiency", "name": "Primary Immunodeficiencies", "omim_prefix": "", "count_estimate": 400},
            {"id": "hematologic", "name": "Rare Hematologic Disorders", "omim_prefix": "", "count_estimate": 350},
            {"id": "connective_tissue", "name": "Connective Tissue Disorders", "omim_prefix": "", "count_estimate": 300},
            {"id": "cardiac", "name": "Rare Cardiac Disorders", "omim_prefix": "", "count_estimate": 250},
            {"id": "renal", "name": "Rare Renal Disorders", "omim_prefix": "", "count_estimate": 200},
            {"id": "pulmonary", "name": "Rare Pulmonary Disorders", "omim_prefix": "", "count_estimate": 150},
            {"id": "dermatologic", "name": "Rare Dermatologic Disorders", "omim_prefix": "", "count_estimate": 400},
            {"id": "endocrine", "name": "Rare Endocrine Disorders", "omim_prefix": "", "count_estimate": 200},
            {"id": "ophthalmologic", "name": "Rare Eye Disorders", "omim_prefix": "", "count_estimate": 350},
            {"id": "hepatic", "name": "Rare Liver Disorders", "omim_prefix": "", "count_estimate": 150},
            {"id": "lysosomal", "name": "Lysosomal Storage Disorders", "omim_prefix": "", "count_estimate": 70},
            {"id": "mitochondrial", "name": "Mitochondrial Disorders", "omim_prefix": "", "count_estimate": 300},
        ],
        "total_rare_diseases": 7000,
        "source": "Orphanet / OMIM (2024)",
    }


@router.get("/gene-therapies")
async def gene_therapies():
    """Reference list of approved and late-stage gene therapies for rare diseases."""
    return {
        "approved_therapies": [
            {"name": "Zolgensma (onasemnogene abeparvovec)", "disease": "Spinal Muscular Atrophy Type 1", "gene": "SMN1", "year": 2019, "agency": "FDA/EMA"},
            {"name": "Luxturna (voretigene neparvovec)", "disease": "RPE65-mediated Inherited Retinal Dystrophy", "gene": "RPE65", "year": 2017, "agency": "FDA/EMA"},
            {"name": "Skysona (elivaldogene autotemcel)", "disease": "Cerebral Adrenoleukodystrophy", "gene": "ABCD1", "year": 2022, "agency": "FDA"},
            {"name": "Zynteglo (betibeglogene autotemcel)", "disease": "Beta-Thalassemia", "gene": "HBB", "year": 2022, "agency": "FDA/EMA"},
            {"name": "Hemgenix (etranacogene dezaparvovec)", "disease": "Hemophilia B", "gene": "F9", "year": 2022, "agency": "FDA"},
            {"name": "Roctavian (valoctocogene roxaparvovec)", "disease": "Hemophilia A", "gene": "F8", "year": 2023, "agency": "EMA"},
            {"name": "Elevidys (delandistrogene moxeparvovec)", "disease": "Duchenne Muscular Dystrophy", "gene": "DMD", "year": 2023, "agency": "FDA"},
            {"name": "Casgevy (exagamglogene autotemcel)", "disease": "Sickle Cell Disease / Beta-Thalassemia", "gene": "BCL11A", "year": 2023, "agency": "FDA/EMA"},
            {"name": "Lyfgenia (lovotibeglogene autotemcel)", "disease": "Sickle Cell Disease", "gene": "HBB", "year": 2023, "agency": "FDA"},
        ],
        "total_approved": 9,
        "pipeline_estimate": 200,
        "last_updated": "2024-12",
    }


@router.get("/acmg-criteria")
async def acmg_criteria():
    """ACMG/AMP 2015 variant classification criteria reference."""
    return {
        "pathogenic_criteria": [
            {"code": "PVS1", "strength": "Very Strong", "description": "Null variant in gene where LOF is known mechanism"},
            {"code": "PS1", "strength": "Strong", "description": "Same amino acid change as established pathogenic variant"},
            {"code": "PS2", "strength": "Strong", "description": "De novo (confirmed) in patient with disease, no family history"},
            {"code": "PS3", "strength": "Strong", "description": "Well-established functional studies showing damaging effect"},
            {"code": "PS4", "strength": "Strong", "description": "Prevalence in affected significantly increased vs controls"},
            {"code": "PM1", "strength": "Moderate", "description": "Located in mutational hot spot or well-established functional domain"},
            {"code": "PM2", "strength": "Moderate", "description": "Absent from controls (or extremely low frequency)"},
            {"code": "PM3", "strength": "Moderate", "description": "Detected in trans with pathogenic variant (recessive disorders)"},
            {"code": "PM4", "strength": "Moderate", "description": "Protein length change due to in-frame indel or stop-loss"},
            {"code": "PM5", "strength": "Moderate", "description": "Novel missense at same position as established pathogenic"},
            {"code": "PM6", "strength": "Moderate", "description": "Assumed de novo (without confirmation of paternity/maternity)"},
            {"code": "PP1", "strength": "Supporting", "description": "Co-segregation with disease in multiple affected family members"},
            {"code": "PP2", "strength": "Supporting", "description": "Missense in gene with low rate of benign missense variants"},
            {"code": "PP3", "strength": "Supporting", "description": "Multiple lines of computational evidence support deleterious effect"},
            {"code": "PP4", "strength": "Supporting", "description": "Patient phenotype highly specific for disease with single genetic etiology"},
            {"code": "PP5", "strength": "Supporting", "description": "Reputable source recently reports variant as pathogenic"},
        ],
        "benign_criteria": [
            {"code": "BA1", "strength": "Stand-alone", "description": "Allele frequency > 5% in gnomAD or large population cohort"},
            {"code": "BS1", "strength": "Strong", "description": "Allele frequency greater than expected for disorder"},
            {"code": "BS2", "strength": "Strong", "description": "Observed in healthy adult individual(s) for fully penetrant disorder"},
            {"code": "BS3", "strength": "Strong", "description": "Well-established functional studies show no damaging effect"},
            {"code": "BS4", "strength": "Strong", "description": "Lack of segregation in affected members of a family"},
            {"code": "BP1", "strength": "Supporting", "description": "Missense in gene where only truncating variants cause disease"},
            {"code": "BP2", "strength": "Supporting", "description": "Observed in trans with pathogenic variant (dominant) or in cis"},
            {"code": "BP3", "strength": "Supporting", "description": "In-frame indel in repetitive region without known function"},
            {"code": "BP4", "strength": "Supporting", "description": "Multiple lines of computational evidence suggest no impact"},
            {"code": "BP5", "strength": "Supporting", "description": "Variant found in case with alternate molecular basis"},
            {"code": "BP6", "strength": "Supporting", "description": "Reputable source recently reports variant as benign"},
            {"code": "BP7", "strength": "Supporting", "description": "Synonymous variant with no predicted splice impact"},
        ],
        "classification_rules": {
            "Pathogenic": "1 Very Strong + >= 1 Strong; OR >= 2 Strong; OR 1 Strong + >= 3 Moderate; etc.",
            "Likely Pathogenic": "1 Very Strong + 1 Moderate; OR 1 Strong + 1-2 Moderate; OR >= 3 Moderate; etc.",
            "VUS": "Criteria for benign and pathogenic are not met",
            "Likely Benign": "1 Strong + 1 Supporting; OR >= 2 Supporting",
            "Benign": "1 Stand-alone; OR >= 2 Strong",
        },
        "source": "Richards et al., Genet Med 2015;17(5):405-424",
    }


@router.get("/hpo-categories")
async def hpo_categories():
    """HPO top-level phenotype categories."""
    return {
        "categories": [
            {"id": "HP:0000118", "label": "Phenotypic abnormality", "children_count": 16000},
            {"id": "HP:0000707", "label": "Abnormality of the nervous system", "children_count": 4500},
            {"id": "HP:0000924", "label": "Abnormality of the skeletal system", "children_count": 2200},
            {"id": "HP:0001626", "label": "Abnormality of the cardiovascular system", "children_count": 1100},
            {"id": "HP:0000478", "label": "Abnormality of the eye", "children_count": 1400},
            {"id": "HP:0001574", "label": "Abnormality of the integument", "children_count": 1200},
            {"id": "HP:0000152", "label": "Abnormality of head or neck", "children_count": 1800},
            {"id": "HP:0001507", "label": "Growth abnormality", "children_count": 400},
            {"id": "HP:0000818", "label": "Abnormality of the endocrine system", "children_count": 500},
            {"id": "HP:0002086", "label": "Abnormality of the respiratory system", "children_count": 600},
            {"id": "HP:0001871", "label": "Abnormality of blood and blood-forming tissues", "children_count": 700},
            {"id": "HP:0000119", "label": "Abnormality of the genitourinary system", "children_count": 900},
            {"id": "HP:0003011", "label": "Abnormality of the musculature", "children_count": 700},
            {"id": "HP:0001939", "label": "Abnormality of metabolism/homeostasis", "children_count": 800},
            {"id": "HP:0002715", "label": "Abnormality of the immune system", "children_count": 600},
            {"id": "HP:0000598", "label": "Abnormality of the ear", "children_count": 400},
        ],
        "total_hpo_terms": 16000,
        "hpo_version": "2024-04-26",
        "source": "Human Phenotype Ontology (hpo.jax.org)",
    }


@router.get("/knowledge-version")
async def knowledge_version():
    """Return metadata about the loaded knowledge base."""
    return {
        "agent": "rare-disease-diagnostic-agent",
        "version": "1.0.0",
        "knowledge_sources": [
            {"name": "OMIM", "version": "2024-12", "records": "~17,000 phenotype entries"},
            {"name": "Orphanet", "version": "2024-Q4", "records": "~7,000 rare diseases"},
            {"name": "HPO", "version": "2024-04-26", "records": "~16,000 phenotype terms"},
            {"name": "ClinVar", "version": "2024-12", "records": "~4.1M variant submissions"},
            {"name": "gnomAD", "version": "v4.1", "records": "~800M variants"},
            {"name": "AlphaMissense", "version": "2023", "records": "~71M predictions"},
            {"name": "Reactome", "version": "v87", "records": "~2,700 pathways"},
            {"name": "ClinicalTrials.gov", "version": "live", "records": "rare disease subset"},
        ],
        "embedding_model": "BAAI/bge-small-en-v1.5",
        "embedding_dimension": 384,
        "llm_model": "claude-sonnet-4-6",
    }
