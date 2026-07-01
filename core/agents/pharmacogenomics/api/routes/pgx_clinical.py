"""PGx-specific clinical decision support endpoints.

Provides pharmacogenomics-specific REST endpoints for clinical workflows
including single-drug PGx checks, polypharmacy medication reviews,
genotype-guided warfarin dosing (IWPC algorithm), HLA hypersensitivity
screening, phenoconversion analysis, and reference gene/drug lookups.

All endpoints use the /v1/pgx prefix and are designed for integration
with clinical decision support systems and EHR platforms.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException

logger = logging.getLogger(__name__)
from pydantic import BaseModel, Field

from src.knowledge import (
    PHARMACOGENES,
    DRUG_CATEGORIES,
    HLA_DRUG_ASSOCIATIONS,
    CYP_INHIBITORS,
    CYP_INDUCERS,
    DRUG_ALTERNATIVES,
    get_gene_context,
    get_drug_context,
    get_hla_context,
)
from src.metrics import record_pipeline_stage
from api.routes.events import emit_event

router = APIRouter(prefix="/v1/pgx", tags=["pgx-clinical"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, cardiology, neurology, and clinical trial agents
    for a comprehensive pharmacogenomic-guided assessment.
    """
    try:
        from src.cross_modal import (
            query_oncology_agent,
            query_cardiology_agent,
            query_neurology_agent,
            query_trial_agent,
            integrate_cross_agent_results,
        )

        patient_profile = request.get("patient_profile", {})
        drug_list = request.get("drug_list", {})
        pgx_profile = request.get("pgx_profile", {})

        results = []

        # Query oncology agent for planned therapy context
        if patient_profile:
            results.append(query_oncology_agent(patient_profile))

        # Query cardiology agent for cardiac drug interactions
        if drug_list:
            results.append(query_cardiology_agent(drug_list))

        # Query neurology agent for neurotoxic drug interactions
        if drug_list:
            results.append(query_neurology_agent(drug_list))

        # Query trial agent for PGx-guided trials
        if pgx_profile:
            results.append(query_trial_agent(pgx_profile))

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

class DrugCheckRequest(BaseModel):
    """Request for a single drug PGx check."""
    drug: str = Field(..., min_length=1, max_length=200, description="Drug name (e.g. codeine, clopidogrel)")
    gene: Optional[str] = Field(None, max_length=50, description="Gene symbol (e.g. CYP2D6)")
    phenotype: Optional[str] = Field(None, max_length=100, description="Metabolizer phenotype (e.g. poor_metabolizer)")
    diplotype: Optional[str] = Field(None, max_length=100, description="Diplotype (e.g. *1/*4)")


class DrugCheckAlert(BaseModel):
    """A single alert from a drug PGx check."""
    alert_level: str = Field(..., description="info, warning, or critical")
    gene: str = ""
    drug: str = ""
    phenotype: str = ""
    recommendation: str = ""
    guideline_body: str = ""
    cpic_level: str = ""
    alternative_drugs: List[str] = Field(default_factory=list)


class DrugCheckResponse(BaseModel):
    """Response for a single drug PGx check."""
    drug: str
    gene: Optional[str] = None
    phenotype: Optional[str] = None
    alerts: List[DrugCheckAlert] = Field(default_factory=list)
    knowledge_context: str = ""
    processing_time_ms: float = 0.0


class MedicationReviewRequest(BaseModel):
    """Request for a polypharmacy medication review."""
    medications: List[str] = Field(..., min_length=1, description="List of current medications")
    gene: Optional[str] = Field(None, max_length=50, description="Gene symbol to focus review on")
    phenotype: Optional[str] = Field(None, max_length=100, description="Metabolizer phenotype")
    include_interactions: bool = Field(True, description="Include drug-drug-gene interactions")


class MedicationInteraction(BaseModel):
    """A detected drug-drug-gene interaction."""
    drug_a: str
    drug_b: str
    enzyme: str = ""
    interaction_type: str = ""
    severity: str = ""
    description: str = ""


class MedicationReviewResponse(BaseModel):
    """Response for a polypharmacy medication review."""
    medications_reviewed: int
    alerts: List[DrugCheckAlert] = Field(default_factory=list)
    interactions: List[MedicationInteraction] = Field(default_factory=list)
    phenoconversion_risks: List[str] = Field(default_factory=list)
    processing_time_ms: float = 0.0


class WarfarinDosingRequest(BaseModel):
    """Request for IWPC warfarin dosing algorithm."""
    age: float = Field(..., ge=18, le=120, description="Patient age in years")
    height_cm: float = Field(..., ge=100, le=250, description="Height in centimeters")
    weight_kg: float = Field(..., ge=30, le=300, description="Weight in kilograms")
    vkorc1_genotype: str = Field(..., description="VKORC1 -1639G>A genotype: GG, GA, or AA")
    cyp2c9_genotype: str = Field(..., description="CYP2C9 genotype: *1/*1, *1/*2, *1/*3, *2/*2, *2/*3, *3/*3")
    race: Optional[str] = Field(None, description="Self-reported race: asian, black, white, other")
    amiodarone: bool = Field(False, description="Concomitant amiodarone use")
    enzyme_inducer: bool = Field(False, description="Concomitant enzyme inducer use (rifampin, carbamazepine, phenytoin)")
    target_inr: float = Field(2.5, ge=1.5, le=4.0, description="Target INR")


class WarfarinDosingResponse(BaseModel):
    """Response for IWPC warfarin dosing algorithm."""
    predicted_weekly_dose_mg: float
    predicted_daily_dose_mg: float
    dose_category: str = Field(..., description="low, standard, or high")
    algorithm: str = "IWPC"
    input_summary: Dict[str, Any] = Field(default_factory=dict)
    warnings: List[str] = Field(default_factory=list)
    processing_time_ms: float = 0.0


class HLAScreenRequest(BaseModel):
    """Request for HLA hypersensitivity screening."""
    drug: str = Field(..., min_length=1, max_length=200, description="Drug to screen")
    hla_alleles: Optional[List[str]] = Field(None, description="Known HLA alleles (e.g. HLA-B*57:01)")


class HLAScreenResult(BaseModel):
    """A single HLA screening result."""
    hla_allele: str
    drug: str
    reaction_type: str = ""
    severity: str = ""
    recommendation: str = ""
    screening_mandatory: bool = False
    prevalence_note: str = ""


class HLAScreenResponse(BaseModel):
    """Response for HLA hypersensitivity screening."""
    drug: str
    associations_found: int
    results: List[HLAScreenResult] = Field(default_factory=list)
    processing_time_ms: float = 0.0


class PhenoconversionRequest(BaseModel):
    """Request for phenoconversion analysis."""
    enzyme: str = Field(..., max_length=50, description="Enzyme to assess (e.g. CYP2D6)")
    baseline_phenotype: str = Field(..., max_length=100, description="Genotype-predicted phenotype (e.g. normal_metabolizer)")
    concomitant_drugs: List[str] = Field(..., min_length=1, description="Concomitant medications to evaluate")


class PhenoconversionResult(BaseModel):
    """A single phenoconversion finding."""
    precipitant_drug: str
    inhibitor_strength: str = ""
    converted_phenotype: str = ""
    clinical_impact: str = ""
    affected_substrates: List[str] = Field(default_factory=list)


class PhenoconversionResponse(BaseModel):
    """Response for phenoconversion analysis."""
    enzyme: str
    baseline_phenotype: str
    effective_phenotype: str
    phenoconversion_detected: bool = False
    results: List[PhenoconversionResult] = Field(default_factory=list)
    processing_time_ms: float = 0.0


class GeneInfo(BaseModel):
    """Summary information about a pharmacogene."""
    gene: str
    full_name: str = ""
    chromosome: str = ""
    function: str = ""
    substrates_count: int = 0
    star_alleles_defined: int = 0
    cpic_guidelines: List[str] = Field(default_factory=list)


class DrugInfo(BaseModel):
    """Summary information about a tracked drug."""
    drug: str
    category: str = ""


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/drug-check", response_model=DrugCheckResponse)
async def drug_check(request: DrugCheckRequest):
    """Perform a single-drug pharmacogenomic check.

    Queries the knowledge graph for gene-drug associations, guideline
    recommendations, and alternative drug suggestions based on the
    patient's phenotype.
    """
    t0 = time.perf_counter()

    try:
        return _drug_check_impl(request, t0)
    except Exception as exc:
        logger.error(f"Drug check failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


def _drug_check_impl(request: DrugCheckRequest, t0: float) -> DrugCheckResponse:
    drug_upper = request.drug.upper()
    alerts: List[DrugCheckAlert] = []
    context_parts = []

    # Look up drug context from knowledge graph
    drug_ctx = get_drug_context(request.drug.lower())
    if drug_ctx:
        context_parts.append(drug_ctx)

    # Look up gene context if provided
    if request.gene:
        gene_ctx = get_gene_context(request.gene.upper())
        if gene_ctx:
            context_parts.append(gene_ctx)

        # Check for gene-drug pair in PHARMACOGENES
        gene_data = PHARMACOGENES.get(request.gene.upper(), {})
        cpic_drugs = gene_data.get("cpic_guidelines", [])

        if request.drug.lower() in [d.lower() for d in cpic_drugs]:
            alert_level = "critical" if request.phenotype and "poor" in request.phenotype.lower() else "warning"
            recommendation = (
                f"CPIC guideline exists for {request.gene.upper()}-{request.drug}. "
            )
            if request.phenotype:
                recommendation += f"Patient phenotype: {request.phenotype}. "

            # Check for alternatives
            alts = []
            for alt_key, alt_data in DRUG_ALTERNATIVES.items():
                if isinstance(alt_data, dict):
                    if (alt_data.get("drug", "").lower() == request.drug.lower()
                            and alt_data.get("gene", "").upper() == request.gene.upper()):
                        for alt_entry in alt_data.get("alternatives", []):
                            if isinstance(alt_entry, dict):
                                alt_drug = alt_entry.get("drug", "")
                            else:
                                alt_drug = str(alt_entry)
                            if alt_drug:
                                alts.append(alt_drug)

            if request.phenotype and "poor" in request.phenotype.lower():
                recommendation += "Consider alternative therapy or significant dose reduction."
            elif request.phenotype and "ultra" in request.phenotype.lower():
                recommendation += "Consider dose reduction or alternative to avoid toxicity from rapid metabolism."

            alerts.append(DrugCheckAlert(
                alert_level=alert_level,
                gene=request.gene.upper(),
                drug=request.drug,
                phenotype=request.phenotype or "",
                recommendation=recommendation.strip(),
                guideline_body="CPIC",
                cpic_level="A",
                alternative_drugs=alts,
            ))
        elif gene_data:
            alerts.append(DrugCheckAlert(
                alert_level="info",
                gene=request.gene.upper(),
                drug=request.drug,
                phenotype=request.phenotype or "",
                recommendation=(
                    f"No CPIC guideline found for {request.gene.upper()}-{request.drug} pair. "
                    f"Gene {request.gene.upper()} has guidelines for: {', '.join(cpic_drugs[:5])}."
                ),
                guideline_body="",
                cpic_level="",
            ))

    # Check HLA associations for this drug
    for assoc_key, assoc_data in HLA_DRUG_ASSOCIATIONS.items():
        if isinstance(assoc_data, dict):
            assoc_drug = assoc_data.get("drug", "")
            if assoc_drug.upper() == drug_upper:
                hla_allele = assoc_data.get("hla_allele", assoc_key)
                alerts.append(DrugCheckAlert(
                    alert_level="critical",
                    gene=hla_allele,
                    drug=request.drug,
                    recommendation=(
                        f"HLA screening required: {hla_allele} testing before {request.drug}. "
                        f"Risk: {assoc_data.get('reaction_type', 'hypersensitivity reaction')}. "
                        f"Severity: {assoc_data.get('severity', 'severe')}."
                    ),
                    guideline_body="CPIC",
                    cpic_level=assoc_data.get("cpic_level", "A"),
                ))

    elapsed_ms = (time.perf_counter() - t0) * 1000
    record_pipeline_stage("drug_check", elapsed_ms / 1000)
    emit_event("drug_check", source="pgx_clinical", summary=f"Drug check: {request.drug}",
               metadata={"drug": request.drug, "gene": request.gene, "alerts": len(alerts)})

    return DrugCheckResponse(
        drug=request.drug,
        gene=request.gene,
        phenotype=request.phenotype,
        alerts=alerts,
        knowledge_context="\n\n".join(context_parts),
        processing_time_ms=round(elapsed_ms, 2),
    )


@router.post("/medication-review", response_model=MedicationReviewResponse)
async def medication_review(request: MedicationReviewRequest):
    """Perform a polypharmacy medication review.

    Checks each medication for PGx implications and identifies
    potential drug-drug-gene interactions and phenoconversion risks
    across the full medication list.
    """
    t0 = time.perf_counter()

    try:
        return _medication_review_impl(request, t0)
    except Exception as exc:
        logger.error(f"Medication review failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


def _medication_review_impl(request: MedicationReviewRequest, t0: float) -> MedicationReviewResponse:
    all_alerts: List[DrugCheckAlert] = []
    interactions: List[MedicationInteraction] = []
    phenoconversion_risks: List[str] = []

    # Check each medication individually
    for med in request.medications:
        med_upper = med.upper()

        # Check HLA associations
        for assoc_key, assoc_data in HLA_DRUG_ASSOCIATIONS.items():
            if isinstance(assoc_data, dict):
                if assoc_data.get("drug", "").upper() == med_upper:
                    hla_allele = assoc_data.get("hla_allele", assoc_key)
                    all_alerts.append(DrugCheckAlert(
                        alert_level="critical",
                        gene=hla_allele,
                        drug=med,
                        recommendation=(
                            f"HLA screening required: {hla_allele} before {med}."
                        ),
                        guideline_body="CPIC",
                        cpic_level=assoc_data.get("cpic_level", "A"),
                    ))

        # Check gene-specific alerts if gene provided
        if request.gene:
            gene_data = PHARMACOGENES.get(request.gene.upper(), {})
            cpic_drugs = gene_data.get("cpic_guidelines", [])
            if med.lower() in [d.lower() for d in cpic_drugs]:
                alert_level = "warning"
                if request.phenotype and "poor" in request.phenotype.lower():
                    alert_level = "critical"
                all_alerts.append(DrugCheckAlert(
                    alert_level=alert_level,
                    gene=request.gene.upper(),
                    drug=med,
                    phenotype=request.phenotype or "",
                    recommendation=(
                        f"CPIC guideline exists for {request.gene.upper()}-{med}. "
                        f"Review dosing based on phenotype."
                    ),
                    guideline_body="CPIC",
                    cpic_level="A",
                ))

    # Check for drug-drug interactions via CYP inhibitors/inducers
    if request.include_interactions:
        med_set_upper = {m.upper() for m in request.medications}

        for enzyme, inhibitor_data in CYP_INHIBITORS.items():
            if not isinstance(inhibitor_data, dict):
                continue
            for strength in ["strong", "moderate", "weak"]:
                inhibitor_list = inhibitor_data.get(strength, [])
                for inhib in inhibitor_list:
                    if inhib.upper() in med_set_upper:
                        # Check if any other medication is a substrate of this enzyme
                        gene_data = PHARMACOGENES.get(enzyme, {})
                        substrates = gene_data.get("cpic_guidelines", [])
                        for substrate in substrates:
                            if substrate.upper() in med_set_upper and substrate.upper() != inhib.upper():
                                interactions.append(MedicationInteraction(
                                    drug_a=inhib,
                                    drug_b=substrate,
                                    enzyme=enzyme,
                                    interaction_type=f"{strength} inhibition",
                                    severity="high" if strength == "strong" else "moderate",
                                    description=(
                                        f"{inhib} is a {strength} {enzyme} inhibitor; "
                                        f"{substrate} is a {enzyme} substrate. Risk of "
                                        f"increased {substrate} exposure."
                                    ),
                                ))

                                # Flag phenoconversion risk
                                phenoconversion_risks.append(
                                    f"{inhib} ({strength} {enzyme} inhibitor) may convert "
                                    f"genotype-predicted phenotype to a poorer metabolizer "
                                    f"status for {substrate}."
                                )

        for enzyme, inducer_data in CYP_INDUCERS.items():
            if not isinstance(inducer_data, dict):
                continue
            for strength in ["strong", "moderate"]:
                inducer_list = inducer_data.get(strength, [])
                for ind in inducer_list:
                    if ind.upper() in med_set_upper:
                        gene_data = PHARMACOGENES.get(enzyme, {})
                        substrates = gene_data.get("cpic_guidelines", [])
                        for substrate in substrates:
                            if substrate.upper() in med_set_upper and substrate.upper() != ind.upper():
                                interactions.append(MedicationInteraction(
                                    drug_a=ind,
                                    drug_b=substrate,
                                    enzyme=enzyme,
                                    interaction_type=f"{strength} induction",
                                    severity="high" if strength == "strong" else "moderate",
                                    description=(
                                        f"{ind} is a {strength} {enzyme} inducer; "
                                        f"{substrate} is a {enzyme} substrate. Risk of "
                                        f"decreased {substrate} efficacy."
                                    ),
                                ))

    elapsed_ms = (time.perf_counter() - t0) * 1000
    record_pipeline_stage("medication_review", elapsed_ms / 1000)
    emit_event("medication_review", source="pgx_clinical",
               summary=f"Reviewed {len(request.medications)} medications",
               metadata={"medications": request.medications, "alerts": len(all_alerts),
                         "interactions": len(interactions)})

    return MedicationReviewResponse(
        medications_reviewed=len(request.medications),
        alerts=all_alerts,
        interactions=interactions,
        phenoconversion_risks=phenoconversion_risks,
        processing_time_ms=round(elapsed_ms, 2),
    )


@router.post("/dosing/warfarin", response_model=WarfarinDosingResponse)
async def warfarin_dosing(request: WarfarinDosingRequest):
    """Calculate genotype-guided warfarin dose using the IWPC algorithm.

    Implements the International Warfarin Pharmacogenetics Consortium
    (IWPC) dosing algorithm incorporating VKORC1, CYP2C9 genotype,
    age, height, weight, race, amiodarone use, and enzyme inducer use.

    Reference: Klein TE et al. NEJM 2009;360:753-764 (PMID: 19228618)
    """
    t0 = time.perf_counter()

    try:
        return _warfarin_dosing_impl(request, t0)
    except Exception as exc:
        logger.error(f"Warfarin dosing calculation failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


def _warfarin_dosing_impl(request: WarfarinDosingRequest, t0: float) -> WarfarinDosingResponse:
    warnings: List[str] = []

    # IWPC algorithm: calculates square root of weekly dose
    # Base constant
    dose_sq = 5.6044

    # Age contribution (decades)
    age_decades = request.age / 10.0
    dose_sq -= 0.2614 * age_decades

    # Height contribution (cm)
    dose_sq += 0.0087 * request.height_cm

    # Weight contribution (kg)
    dose_sq += 0.0128 * request.weight_kg

    # VKORC1 -1639G>A contribution
    vkorc1 = request.vkorc1_genotype.upper()
    if vkorc1 == "GA" or vkorc1 == "AG":
        dose_sq -= 0.8677
    elif vkorc1 == "AA":
        dose_sq -= 1.6974
    elif vkorc1 == "GG":
        pass  # Reference genotype, no adjustment
    else:
        warnings.append(f"Unrecognized VKORC1 genotype '{request.vkorc1_genotype}'. Using GG (reference).")

    # CYP2C9 contribution
    cyp2c9 = request.cyp2c9_genotype.upper().replace(" ", "")
    cyp2c9_adjustments = {
        "*1/*2": -0.5211,
        "*1/*3": -0.9357,
        "*2/*2": -1.0616,
        "*2/*3": -1.9206,
        "*3/*3": -2.3312,
        "*1/*1": 0.0,
    }
    if cyp2c9 in cyp2c9_adjustments:
        dose_sq += cyp2c9_adjustments[cyp2c9]
    else:
        warnings.append(f"Unrecognized CYP2C9 genotype '{request.cyp2c9_genotype}'. Using *1/*1 (reference).")

    # Race contribution
    if request.race:
        race_lower = request.race.lower()
        if race_lower == "asian":
            dose_sq -= 0.1092
        elif race_lower == "black" or race_lower == "african_american":
            dose_sq -= 0.2760
        elif race_lower == "white" or race_lower == "caucasian":
            pass  # Reference
        # "other" or unrecognized: no adjustment

    # Amiodarone
    if request.amiodarone:
        dose_sq -= 0.5503
        warnings.append("Amiodarone is a potent CYP2C9 inhibitor; close INR monitoring recommended.")

    # Enzyme inducer (rifampin, carbamazepine, phenytoin)
    if request.enzyme_inducer:
        dose_sq += 1.2799
        warnings.append("Enzyme inducer detected; higher warfarin doses may be needed with close monitoring.")

    # Calculate weekly dose (square the result)
    weekly_dose = dose_sq ** 2 if dose_sq > 0 else 0.0
    weekly_dose = round(weekly_dose, 1)
    daily_dose = round(weekly_dose / 7.0, 1)

    # Categorize dose
    if weekly_dose < 21:
        dose_category = "low"
        warnings.append("Predicted dose is in the low range (<21 mg/week). Start conservatively.")
    elif weekly_dose > 49:
        dose_category = "high"
        warnings.append("Predicted dose is in the high range (>49 mg/week). Verify genotype and clinical factors.")
    else:
        dose_category = "standard"

    # Safety bounds
    if weekly_dose < 3:
        weekly_dose = 3.0
        daily_dose = round(3.0 / 7.0, 1)
        warnings.append("Predicted dose below clinical minimum; clamped to 3 mg/week.")
    elif weekly_dose > 105:
        weekly_dose = 105.0
        daily_dose = round(105.0 / 7.0, 1)
        warnings.append("Predicted dose above clinical maximum; clamped to 105 mg/week.")

    elapsed_ms = (time.perf_counter() - t0) * 1000
    record_pipeline_stage("warfarin_dosing", elapsed_ms / 1000)
    emit_event("dosing", source="pgx_clinical",
               summary=f"Warfarin IWPC dosing: {weekly_dose:.1f} mg/week ({dose_category})",
               metadata={"algorithm": "IWPC", "weekly_dose_mg": weekly_dose, "category": dose_category})

    return WarfarinDosingResponse(
        predicted_weekly_dose_mg=weekly_dose,
        predicted_daily_dose_mg=daily_dose,
        dose_category=dose_category,
        algorithm="IWPC",
        input_summary={
            "age": request.age,
            "height_cm": request.height_cm,
            "weight_kg": request.weight_kg,
            "vkorc1": request.vkorc1_genotype,
            "cyp2c9": request.cyp2c9_genotype,
            "race": request.race,
            "amiodarone": request.amiodarone,
            "enzyme_inducer": request.enzyme_inducer,
            "target_inr": request.target_inr,
        },
        warnings=warnings,
        processing_time_ms=round(elapsed_ms, 2),
    )


@router.post("/hla-screen", response_model=HLAScreenResponse)
async def hla_screen(request: HLAScreenRequest):
    """Screen a drug for HLA-mediated hypersensitivity associations.

    Queries the knowledge graph for known HLA-drug associations and
    returns screening recommendations, reaction types, severity levels,
    and prevalence information.
    """
    t0 = time.perf_counter()

    try:
        return _hla_screen_impl(request, t0)
    except Exception as exc:
        logger.error(f"HLA screening failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


def _hla_screen_impl(request: HLAScreenRequest, t0: float) -> HLAScreenResponse:
    drug_upper = request.drug.upper()
    results: List[HLAScreenResult] = []

    for assoc_key, assoc_data in HLA_DRUG_ASSOCIATIONS.items():
        if not isinstance(assoc_data, dict):
            continue
        assoc_drug = assoc_data.get("drug", "")
        if assoc_drug.upper() != drug_upper:
            continue

        hla_allele = assoc_data.get("hla_allele", assoc_key)

        # Check if patient has this allele (if alleles provided)
        patient_positive = False
        if request.hla_alleles:
            patient_positive = any(
                a.upper().replace(" ", "") == hla_allele.upper().replace(" ", "")
                for a in request.hla_alleles
            )

        prevalence_parts = []
        prev_by_pop = assoc_data.get("prevalence_by_population", {})
        if isinstance(prev_by_pop, dict):
            for pop, freq in prev_by_pop.items():
                prevalence_parts.append(f"{pop}: {freq}")
        prevalence_note = "; ".join(prevalence_parts) if prevalence_parts else ""

        recommendation = assoc_data.get("recommendation", "")
        if not recommendation:
            if assoc_data.get("screening_mandatory", False):
                recommendation = f"Mandatory {hla_allele} screening before {request.drug} initiation."
            else:
                recommendation = f"Consider {hla_allele} screening before {request.drug}."

        if patient_positive:
            recommendation = f"POSITIVE for {hla_allele}: {request.drug} is CONTRAINDICATED. " + recommendation

        results.append(HLAScreenResult(
            hla_allele=hla_allele,
            drug=request.drug,
            reaction_type=assoc_data.get("reaction_type", ""),
            severity=assoc_data.get("severity", ""),
            recommendation=recommendation,
            screening_mandatory=assoc_data.get("screening_mandatory", False),
            prevalence_note=prevalence_note,
        ))

    # Also check knowledge context
    get_hla_context(request.drug)

    elapsed_ms = (time.perf_counter() - t0) * 1000
    record_pipeline_stage("hla_screen", elapsed_ms / 1000)
    emit_event("hla_screen", source="pgx_clinical",
               summary=f"HLA screen: {request.drug} ({len(results)} associations)",
               metadata={"drug": request.drug, "associations": len(results)})

    return HLAScreenResponse(
        drug=request.drug,
        associations_found=len(results),
        results=results,
        processing_time_ms=round(elapsed_ms, 2),
    )


@router.post("/phenoconversion", response_model=PhenoconversionResponse)
async def phenoconversion_analysis(request: PhenoconversionRequest):
    """Analyze potential phenoconversion from concomitant medications.

    Evaluates whether concomitant CYP inhibitors or inducers may alter
    the genotype-predicted metabolizer phenotype for a given enzyme,
    and identifies affected substrate drugs.
    """
    t0 = time.perf_counter()

    try:
        return _phenoconversion_impl(request, t0)
    except Exception as exc:
        logger.error(f"Phenoconversion analysis failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


def _phenoconversion_impl(request: PhenoconversionRequest, t0: float) -> PhenoconversionResponse:
    enzyme_upper = request.enzyme.upper()
    results: List[PhenoconversionResult] = []
    phenoconversion_detected = False
    effective_phenotype = request.baseline_phenotype

    # Check each concomitant drug against CYP inhibitors
    inhibitor_data = CYP_INHIBITORS.get(enzyme_upper, {})
    if isinstance(inhibitor_data, dict):
        for drug in request.concomitant_drugs:
            drug_upper = drug.upper()
            for strength in ["strong", "moderate", "weak"]:
                inhibitor_list = inhibitor_data.get(strength, [])
                if drug_upper in [i.upper() for i in inhibitor_list]:
                    phenoconversion_detected = True

                    # Determine converted phenotype
                    baseline_lower = request.baseline_phenotype.lower()
                    if strength == "strong":
                        if "normal" in baseline_lower or "rapid" in baseline_lower:
                            converted = "poor_metabolizer"
                        elif "intermediate" in baseline_lower:
                            converted = "poor_metabolizer"
                        else:
                            converted = baseline_lower
                        effective_phenotype = converted
                    elif strength == "moderate":
                        if "normal" in baseline_lower or "rapid" in baseline_lower:
                            converted = "intermediate_metabolizer"
                        elif "intermediate" in baseline_lower:
                            converted = "poor_metabolizer"
                        else:
                            converted = baseline_lower
                        if "poor" in effective_phenotype.lower() and "poor" not in converted.lower():
                            pass  # Keep the more severe phenotype
                        else:
                            effective_phenotype = converted
                    else:  # weak
                        converted = baseline_lower
                        # Weak inhibitors generally do not cause clinically significant phenoconversion

                    # Get affected substrates
                    gene_data = PHARMACOGENES.get(enzyme_upper, {})
                    substrates = gene_data.get("cpic_guidelines", [])

                    clinical_impact = ""
                    if strength == "strong":
                        clinical_impact = (
                            f"Strong {enzyme_upper} inhibition by {drug} converts "
                            f"{request.baseline_phenotype} to {converted}. "
                            f"Treat patient as {converted} for {enzyme_upper} substrate dosing."
                        )
                    elif strength == "moderate":
                        clinical_impact = (
                            f"Moderate {enzyme_upper} inhibition by {drug} may convert "
                            f"{request.baseline_phenotype} to {converted}. "
                            f"Consider phenotype adjustment for {enzyme_upper} substrates."
                        )
                    else:
                        clinical_impact = (
                            f"Weak {enzyme_upper} inhibition by {drug}. Generally not "
                            f"clinically significant for phenoconversion."
                        )

                    results.append(PhenoconversionResult(
                        precipitant_drug=drug,
                        inhibitor_strength=strength,
                        converted_phenotype=converted,
                        clinical_impact=clinical_impact,
                        affected_substrates=substrates[:10],
                    ))
                    break  # Found strength for this drug, move to next drug

    # Check CYP inducers
    inducer_data = CYP_INDUCERS.get(enzyme_upper, {})
    if isinstance(inducer_data, dict):
        for drug in request.concomitant_drugs:
            drug_upper = drug.upper()
            for strength in ["strong", "moderate"]:
                inducer_list = inducer_data.get(strength, [])
                if drug_upper in [i.upper() for i in inducer_list]:
                    phenoconversion_detected = True

                    baseline_lower = request.baseline_phenotype.lower()
                    if "poor" in baseline_lower:
                        converted = "intermediate_metabolizer"
                    elif "intermediate" in baseline_lower:
                        converted = "normal_metabolizer"
                    elif "normal" in baseline_lower:
                        converted = "ultra_rapid_metabolizer"
                    else:
                        converted = baseline_lower

                    effective_phenotype = converted

                    gene_data = PHARMACOGENES.get(enzyme_upper, {})
                    substrates = gene_data.get("cpic_guidelines", [])

                    results.append(PhenoconversionResult(
                        precipitant_drug=drug,
                        inhibitor_strength=f"{strength}_inducer",
                        converted_phenotype=converted,
                        clinical_impact=(
                            f"{strength.capitalize()} {enzyme_upper} induction by {drug} "
                            f"may convert {request.baseline_phenotype} to {converted}. "
                            f"Substrate drugs may have reduced efficacy."
                        ),
                        affected_substrates=substrates[:10],
                    ))
                    break

    elapsed_ms = (time.perf_counter() - t0) * 1000
    record_pipeline_stage("phenoconversion", elapsed_ms / 1000)
    emit_event("phenoconversion", source="pgx_clinical",
               summary=f"Phenoconversion: {request.enzyme} {'detected' if phenoconversion_detected else 'none'}",
               metadata={"enzyme": request.enzyme, "detected": phenoconversion_detected,
                         "drugs_checked": len(request.concomitant_drugs)})

    return PhenoconversionResponse(
        enzyme=request.enzyme,
        baseline_phenotype=request.baseline_phenotype,
        effective_phenotype=effective_phenotype,
        phenoconversion_detected=phenoconversion_detected,
        results=results,
        processing_time_ms=round(elapsed_ms, 2),
    )


@router.get("/genes", response_model=List[GeneInfo])
async def list_pharmacogenes():
    """List all pharmacogenes in the knowledge graph.

    Returns summary information for each tracked pharmacogene including
    gene symbol, full name, chromosome location, function, number of
    substrates, star alleles, and CPIC guideline drugs.
    """
    genes = []
    for gene_symbol, data in sorted(PHARMACOGENES.items()):
        if not isinstance(data, dict):
            continue
        genes.append(GeneInfo(
            gene=gene_symbol,
            full_name=data.get("full_name", ""),
            chromosome=data.get("chromosome", ""),
            function=data.get("function", ""),
            substrates_count=data.get("substrates_count") or 0,
            star_alleles_defined=data.get("star_alleles_defined") or 0,
            cpic_guidelines=data.get("cpic_guidelines", []),
        ))
    return genes


@router.get("/drugs", response_model=List[DrugInfo])
async def list_tracked_drugs():
    """List all drugs tracked across PGx drug categories.

    Returns each drug with its therapeutic category (e.g. opioids,
    anticoagulants, antidepressants, statins, chemotherapy, etc.).
    """
    drugs = []
    for category, cat_data in sorted(DRUG_CATEGORIES.items()):
        if not isinstance(cat_data, dict):
            continue
        drug_list = cat_data.get("drugs", [])
        for drug in drug_list:
            drugs.append(DrugInfo(
                drug=drug,
                category=category,
            ))
    return drugs
