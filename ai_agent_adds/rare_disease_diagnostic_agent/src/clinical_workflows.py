"""Clinical workflows for the Rare Disease Diagnostic Agent.

Author: Adam Jones
Date: March 2026

Implements ten evidence-based rare disease diagnostic workflows covering
phenotype-driven diagnosis, WES/WGS interpretation, metabolic screening,
dysmorphology assessment, neurogenetic evaluation, cardiac genetics,
connective tissue disorders, inborn errors of metabolism, gene therapy
eligibility, and undiagnosed disease workup.

Each workflow follows the BaseRareDiseaseWorkflow contract
(preprocess -> execute -> postprocess) and is registered in the
WorkflowEngine for unified dispatch.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Dict, List, Optional

from src.models import (
    ACMGClassification,
    DiseaseCandidate,
    DiseaseCategory,
    DiagnosticResult,
    DiagnosticWorkflowType,
    InheritancePattern,
    SeverityLevel,
    TherapyMatch,
    TherapyStatus,
    VariantClassification,
    VariantType,
    WorkflowResult,
)

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

_SEVERITY_ORDER: List[SeverityLevel] = [
    SeverityLevel.INFORMATIONAL,
    SeverityLevel.LOW,
    SeverityLevel.MODERATE,
    SeverityLevel.HIGH,
    SeverityLevel.CRITICAL,
]


def _max_severity(*levels: SeverityLevel) -> SeverityLevel:
    """Return the highest severity among the given levels."""
    return max(levels, key=lambda s: _SEVERITY_ORDER.index(s))


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    """Clamp a numeric value to [lo, hi]."""
    return max(lo, min(hi, value))


def _trigger_string(trigger_type: str, detail: str) -> str:
    """Build a human-readable cross-agent trigger string."""
    return f"[{trigger_type}] {detail}"


# ═══════════════════════════════════════════════════════════════════════════
# BASE CLASS
# ═══════════════════════════════════════════════════════════════════════════


class BaseRareDiseaseWorkflow(ABC):
    """Abstract base for all rare disease diagnostic workflows."""

    workflow_type: DiagnosticWorkflowType

    # ── template-method orchestrator ──────────────────────────────────
    def run(self, inputs: dict) -> WorkflowResult:
        """Orchestrate preprocess -> execute -> postprocess."""
        logger.info("Running workflow %s", self.workflow_type.value)
        processed_inputs = self.preprocess(inputs)
        result = self.execute(processed_inputs)
        result = self.postprocess(result)
        # Inject any validation warnings as findings
        warnings = processed_inputs.get("_validation_warnings", [])
        if warnings:
            result.findings = [
                f"[INPUT WARNING] {w}" for w in warnings
            ] + result.findings
        return result

    def preprocess(self, inputs: dict) -> dict:
        """Validate and normalise raw inputs.  Override for workflow-specific logic."""
        return dict(inputs)

    @abstractmethod
    def execute(self, inputs: dict) -> WorkflowResult:
        """Core workflow logic.  Must be implemented by each workflow."""
        ...

    def postprocess(self, result: WorkflowResult) -> WorkflowResult:
        """Shared enrichment after execution."""
        try:
            from api.routes.events import publish_event
            publish_event("workflow_complete", {
                "workflow": result.workflow_type.value if hasattr(result.workflow_type, 'value') else str(result.workflow_type),
                "severity": result.severity.value if hasattr(result.severity, 'value') else str(result.severity),
                "findings_count": len(result.findings),
            })
        except Exception:
            pass  # Don't break workflow for event publishing failure
        return result

    @staticmethod
    def _init_warnings(inp: dict) -> list:
        """Initialise and return the validation warnings list on *inp*."""
        warnings: list = inp.setdefault("_validation_warnings", [])
        return warnings


# ═══════════════════════════════════════════════════════════════════════════
# HPO-DISEASE KNOWLEDGE BASE (inline for self-contained workflows)
# ═══════════════════════════════════════════════════════════════════════════

_HPO_DISEASE_DB: List[Dict] = [
    {"disease_id": "OMIM:219700", "name": "Cystic Fibrosis", "hpo": ["HP:0002205", "HP:0006538", "HP:0001508", "HP:0002110", "HP:0006536"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["CFTR"], "prevalence": "1:3,500", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:310200", "name": "Duchenne Muscular Dystrophy", "hpo": ["HP:0003202", "HP:0003236", "HP:0001644", "HP:0002093", "HP:0003560"],
     "inheritance": InheritancePattern.X_LINKED_RECESSIVE, "genes": ["DMD"], "prevalence": "1:3,500 males", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:141900", "name": "Huntington Disease", "hpo": ["HP:0002072", "HP:0000726", "HP:0001300", "HP:0001289", "HP:0002354"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["HTT"], "prevalence": "1:10,000", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:154700", "name": "Marfan Syndrome", "hpo": ["HP:0001166", "HP:0001519", "HP:0004382", "HP:0002816", "HP:0001083"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["FBN1"], "prevalence": "1:5,000", "category": DiseaseCategory.CONNECTIVE_TISSUE},
    {"disease_id": "OMIM:256700", "name": "Maple Syrup Urine Disease", "hpo": ["HP:0001250", "HP:0001252", "HP:0003217", "HP:0001263", "HP:0001508"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["BCKDHA", "BCKDHB", "DBT"], "prevalence": "1:185,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:261600", "name": "Phenylketonuria", "hpo": ["HP:0001249", "HP:0001250", "HP:0001252", "HP:0000964", "HP:0001263"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["PAH"], "prevalence": "1:10,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:232300", "name": "Gaucher Disease", "hpo": ["HP:0001744", "HP:0002240", "HP:0001882", "HP:0002653", "HP:0001376"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["GBA1"], "prevalence": "1:40,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:230800", "name": "Galactosemia", "hpo": ["HP:0001399", "HP:0002240", "HP:0000952", "HP:0001250", "HP:0001263"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["GALT"], "prevalence": "1:30,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:607014", "name": "Dravet Syndrome", "hpo": ["HP:0001250", "HP:0001263", "HP:0001249", "HP:0001252", "HP:0002373"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["SCN1A"], "prevalence": "1:15,700", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:192500", "name": "Long QT Syndrome", "hpo": ["HP:0001657", "HP:0004756", "HP:0001279", "HP:0001678", "HP:0001695"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["KCNQ1", "KCNH2", "SCN5A"], "prevalence": "1:2,500", "category": DiseaseCategory.CARDIAC},
    {"disease_id": "OMIM:115200", "name": "Hypertrophic Cardiomyopathy", "hpo": ["HP:0001639", "HP:0001635", "HP:0001663", "HP:0002875", "HP:0012764"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["MYH7", "MYBPC3"], "prevalence": "1:500", "category": DiseaseCategory.CARDIAC},
    {"disease_id": "OMIM:130050", "name": "Ehlers-Danlos Syndrome, Hypermobility Type", "hpo": ["HP:0001382", "HP:0000974", "HP:0001252", "HP:0002758", "HP:0001388"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["COL5A1", "COL5A2", "TNXB"], "prevalence": "1:5,000", "category": DiseaseCategory.CONNECTIVE_TISSUE},
    {"disease_id": "OMIM:253300", "name": "Spinal Muscular Atrophy", "hpo": ["HP:0003202", "HP:0001252", "HP:0001319", "HP:0002093", "HP:0001558"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["SMN1"], "prevalence": "1:10,000", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:300624", "name": "Rett Syndrome", "hpo": ["HP:0001263", "HP:0001249", "HP:0001250", "HP:0002167", "HP:0000733"],
     "inheritance": InheritancePattern.X_LINKED_DOMINANT, "genes": ["MECP2"], "prevalence": "1:10,000 females", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:176270", "name": "Prader-Willi Syndrome", "hpo": ["HP:0001252", "HP:0001263", "HP:0001513", "HP:0000046", "HP:0001249"],
     "inheritance": InheritancePattern.DE_NOVO, "genes": ["SNRPN"], "prevalence": "1:15,000", "category": DiseaseCategory.NEUROLOGICAL},
    {"disease_id": "OMIM:601419", "name": "Arrhythmogenic Right Ventricular Cardiomyopathy", "hpo": ["HP:0001695", "HP:0001663", "HP:0001279", "HP:0001714", "HP:0001635"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["PKP2", "DSP", "DSG2", "DSC2", "JUP"], "prevalence": "1:5,000", "category": DiseaseCategory.CARDIAC},
    {"disease_id": "OMIM:237300", "name": "Ornithine Transcarbamylase Deficiency", "hpo": ["HP:0001987", "HP:0001250", "HP:0001259", "HP:0002910", "HP:0003571"],
     "inheritance": InheritancePattern.X_LINKED_RECESSIVE, "genes": ["OTC"], "prevalence": "1:14,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:252010", "name": "Methylmalonic Acidemia", "hpo": ["HP:0001508", "HP:0001250", "HP:0001252", "HP:0001263", "HP:0001942"],
     "inheritance": InheritancePattern.AUTOSOMAL_RECESSIVE, "genes": ["MUT", "MMAA", "MMAB"], "prevalence": "1:50,000", "category": DiseaseCategory.METABOLIC},
    {"disease_id": "OMIM:200100", "name": "Achondroplasia", "hpo": ["HP:0004322", "HP:0000256", "HP:0003027", "HP:0002986", "HP:0000286"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["FGFR3"], "prevalence": "1:15,000", "category": DiseaseCategory.SKELETAL},
    {"disease_id": "OMIM:211980", "name": "Osteogenesis Imperfecta Type I", "hpo": ["HP:0002757", "HP:0000678", "HP:0000365", "HP:0000592", "HP:0001382"],
     "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "genes": ["COL1A1", "COL1A2"], "prevalence": "1:15,000", "category": DiseaseCategory.SKELETAL},
]


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 1 — Phenotype-Driven Diagnosis
# ═══════════════════════════════════════════════════════════════════════════


class PhenotypeDrivenWorkflow(BaseRareDiseaseWorkflow):
    """Phenotype-driven differential diagnosis using HPO terms.

    Accepts patient HPO terms, scores candidate diseases by phenotypic
    overlap (matched/unmatched terms) and inheritance pattern fit,
    and returns a ranked DiseaseCandidate list.

    Inputs
    ------
    hpo_terms : list[str]
        Patient HPO term IDs (e.g., ['HP:0001250', 'HP:0001263']).
    sex : str, optional
        Patient sex for X-linked inheritance filtering.
    age : str, optional
        Patient age for age-of-onset filtering.
    consanguinity : bool
        Whether parental consanguinity is reported.
    family_history : str, optional
        Free-text family history.
    top_k : int
        Number of candidates to return (default 10).
    """

    workflow_type = DiagnosticWorkflowType.PHENOTYPE_DRIVEN

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("hpo_terms"):
            warnings.append("No HPO terms provided — differential will be empty")
            inp["hpo_terms"] = []
        elif isinstance(inp["hpo_terms"], str):
            inp["hpo_terms"] = [t.strip() for t in inp["hpo_terms"].split(",")]

        inp.setdefault("sex", "unknown")
        inp.setdefault("age", None)
        inp.setdefault("consanguinity", False)
        inp.setdefault("family_history", None)
        inp.setdefault("top_k", 10)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        candidates: List[DiseaseCandidate] = []

        hpo_terms = set(inputs["hpo_terms"])
        sex = inputs.get("sex", "unknown")
        consanguinity = inputs.get("consanguinity", False)
        top_k = inputs.get("top_k", 10)

        if not hpo_terms:
            findings.append("No HPO terms provided; cannot generate differential")
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=findings,
                recommendations=["Provide HPO terms for phenotype-driven analysis"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        findings.append(f"Evaluating {len(hpo_terms)} patient HPO terms against {len(_HPO_DISEASE_DB)} diseases")

        # Score each disease
        scored: List[Dict] = []
        for disease in _HPO_DISEASE_DB:
            disease_hpo = set(disease["hpo"])
            matched = hpo_terms & disease_hpo
            unmatched = disease_hpo - hpo_terms

            if not matched:
                continue

            # Jaccard-like similarity
            union = hpo_terms | disease_hpo
            similarity = len(matched) / len(union) if union else 0.0

            # Inheritance fit bonus
            inheritance_bonus = 0.0
            inh = disease.get("inheritance")
            if consanguinity and inh == InheritancePattern.AUTOSOMAL_RECESSIVE:
                inheritance_bonus = 0.10
            if sex == "male" and inh in (InheritancePattern.X_LINKED_RECESSIVE, InheritancePattern.X_LINKED_DOMINANT):
                inheritance_bonus = 0.05
            if inh == InheritancePattern.DE_NOVO:
                inheritance_bonus = 0.02

            total_score = _clamp(similarity + inheritance_bonus)

            scored.append({
                "disease": disease,
                "matched": list(matched),
                "unmatched": list(unmatched),
                "score": total_score,
            })

        # Sort by score descending
        scored.sort(key=lambda x: x["score"], reverse=True)
        scored = scored[:top_k]

        for rank, entry in enumerate(scored, 1):
            d = entry["disease"]
            candidates.append(DiseaseCandidate(
                disease_id=d["disease_id"],
                disease_name=d["name"],
                rank=rank,
                similarity_score=round(entry["score"], 4),
                matched_phenotypes=entry["matched"],
                unmatched_phenotypes=entry["unmatched"],
                inheritance_pattern=d.get("inheritance"),
                prevalence=d.get("prevalence"),
                category=d.get("category"),
                causal_genes=d.get("genes", []),
            ))
            findings.append(
                f"  #{rank} {d['name']} (score={entry['score']:.3f}, "
                f"matched={len(entry['matched'])}, unmatched={len(entry['unmatched'])})"
            )

        if candidates:
            top = candidates[0]
            recommendations.append(
                f"Top candidate: {top.disease_name} — consider confirmatory "
                f"testing for gene(s): {', '.join(top.causal_genes)}"
            )
            if consanguinity:
                recommendations.append(
                    "Consanguinity reported — prioritise autosomal recessive candidates"
                )

        confidence = candidates[0].similarity_score if candidates else 0.0
        severity = SeverityLevel.MODERATE if confidence > 0.3 else SeverityLevel.LOW

        references.append("Kohler et al. (2021) The Human Phenotype Ontology in 2021. NAR")
        references.append("Robinson et al. (2014) Improved exome prioritization of disease genes. Genome Res")

        diagnostic = DiagnosticResult(
            candidate_diseases=candidates,
            confidence=confidence,
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            confidence=confidence,
            diagnostic_result=diagnostic,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 2 — WES/WGS Interpretation
# ═══════════════════════════════════════════════════════════════════════════


class WESWGSInterpretationWorkflow(BaseRareDiseaseWorkflow):
    """WES/WGS variant interpretation workflow.

    Filters VCF variants by population frequency (gnomAD <0.01 for dominant,
    <0.05 for recessive), applies simplified ACMG classification, and ranks
    candidate genes.

    Inputs
    ------
    variants : list[dict]
        Each dict: gene, variant_id, population_frequency, zygosity,
        variant_type, hgvs_c, hgvs_p, in_clinvar, clinvar_classification,
        computational_prediction.
    inheritance_filter : str, optional
        Expected inheritance pattern for filtering.
    hpo_terms : list[str], optional
        Patient HPO terms for phenotype-gene correlation.
    """

    workflow_type = DiagnosticWorkflowType.WES_WGS_INTERPRETATION

    # gnomAD frequency thresholds
    _FREQ_THRESHOLD_DOMINANT = 0.01
    _FREQ_THRESHOLD_RECESSIVE = 0.05

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("variants"):
            warnings.append("No variants provided — interpretation will be empty")
            inp["variants"] = []
        inp.setdefault("inheritance_filter", None)
        inp.setdefault("hpo_terms", [])
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        classified_variants: List[VariantClassification] = []

        variants = inputs.get("variants", [])
        hpo_terms = set(inputs.get("hpo_terms", []))

        if not variants:
            findings.append("No variants provided for interpretation")
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=findings,
                recommendations=["Submit VCF variants for WES/WGS interpretation"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        findings.append(f"Received {len(variants)} variants for interpretation")

        # Step 1: Frequency filter
        passing: List[Dict] = []
        for v in variants:
            freq = float(v.get("population_frequency", 0.0))
            zyg = v.get("zygosity", "heterozygous")
            threshold = self._FREQ_THRESHOLD_RECESSIVE if zyg in ("homozygous", "compound_het") else self._FREQ_THRESHOLD_DOMINANT
            if freq <= threshold:
                passing.append(v)

        findings.append(f"After frequency filter: {len(passing)}/{len(variants)} variants pass")

        # Step 2: ACMG classification (simplified)
        for v in passing:
            criteria: List[str] = []
            score = 0

            # PVS1: null variant in gene where LOF is known mechanism
            vtype = v.get("variant_type", "")
            if vtype in ("frameshift", "nonsense", "splice_site", "deletion"):
                criteria.append("PVS1")
                score += 8

            # PM2: absent from controls (very low frequency)
            freq = float(v.get("population_frequency", 0.0))
            if freq < 0.0001:
                criteria.append("PM2")
                score += 2

            # PP3: computational prediction supports deleterious
            if v.get("computational_prediction") in ("deleterious", "damaging", "pathogenic"):
                criteria.append("PP3")
                score += 1

            # PS1: same amino acid change as established pathogenic
            if v.get("in_clinvar") and v.get("clinvar_classification") in ("pathogenic", "likely_pathogenic"):
                criteria.append("PS1")
                score += 4

            # PP4: phenotype specificity
            gene = v.get("gene", "")
            gene_matches_phenotype = False
            for disease in _HPO_DISEASE_DB:
                if gene in disease.get("genes", []):
                    disease_hpo = set(disease["hpo"])
                    if hpo_terms & disease_hpo:
                        gene_matches_phenotype = True
                        break
            if gene_matches_phenotype:
                criteria.append("PP4")
                score += 1

            # BA1: population frequency too high -> benign
            if freq > 0.05:
                classification = ACMGClassification.BENIGN
            elif score >= 10:
                classification = ACMGClassification.PATHOGENIC
            elif score >= 6:
                classification = ACMGClassification.LIKELY_PATHOGENIC
            elif score <= 0:
                classification = ACMGClassification.LIKELY_BENIGN
            else:
                classification = ACMGClassification.VUS

            evidence_summary = (
                f"ACMG criteria met: {', '.join(criteria) if criteria else 'none'}. "
                f"Score: {score}. Frequency: {freq:.6f}."
            )

            classified_variants.append(VariantClassification(
                variant_id=v.get("variant_id", "unknown"),
                gene=gene,
                hgvs_c=v.get("hgvs_c"),
                hgvs_p=v.get("hgvs_p"),
                variant_type=VariantType(vtype) if vtype in [e.value for e in VariantType] else None,
                classification=classification,
                acmg_criteria=criteria,
                population_frequency=freq if freq <= 1.0 else None,
                evidence_summary=evidence_summary,
                zygosity=v.get("zygosity"),
            ))

        # Rank by pathogenicity
        path_order = {
            ACMGClassification.PATHOGENIC: 0,
            ACMGClassification.LIKELY_PATHOGENIC: 1,
            ACMGClassification.VUS: 2,
            ACMGClassification.LIKELY_BENIGN: 3,
            ACMGClassification.BENIGN: 4,
        }
        classified_variants.sort(key=lambda vc: path_order.get(vc.classification, 5))

        pathogenic_count = sum(
            1 for vc in classified_variants
            if vc.classification in (ACMGClassification.PATHOGENIC, ACMGClassification.LIKELY_PATHOGENIC)
        )
        vus_count = sum(1 for vc in classified_variants if vc.classification == ACMGClassification.VUS)
        findings.append(f"Classification: {pathogenic_count} pathogenic/likely pathogenic, {vus_count} VUS")

        if pathogenic_count > 0:
            recommendations.append("Review pathogenic/likely pathogenic variants with clinical team")
            recommendations.append("Consider segregation analysis in parents if available")
        if vus_count > 0:
            recommendations.append(f"{vus_count} VUS identified — consider functional studies or reanalysis in 12 months")

        confidence = _clamp(0.3 + 0.1 * pathogenic_count)
        severity = SeverityLevel.HIGH if pathogenic_count > 0 else SeverityLevel.MODERATE

        references.append("Richards et al. (2015) Standards and guidelines for interpretation of sequence variants. Genet Med")
        references.append("ClinGen Sequence Variant Interpretation (SVI) Working Group recommendations")

        diagnostic = DiagnosticResult(
            variants=classified_variants,
            confidence=confidence,
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            confidence=confidence,
            diagnostic_result=diagnostic,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 3 — Metabolic Screening
# ═══════════════════════════════════════════════════════════════════════════


class MetabolicScreeningWorkflow(BaseRareDiseaseWorkflow):
    """Metabolic screening workflow for inborn errors of metabolism.

    Accepts analyte values, performs pathway analysis, generates
    confirmatory test recommendations, and includes emergency protocols
    for metabolic crisis (hyperammonemia, metabolic acidosis).

    Inputs
    ------
    analytes : dict
        Analyte name -> value dict (e.g., {"phenylalanine": 1200, "ammonia": 350}).
    age : str, optional
        Patient age (neonatal presentation changes thresholds).
    clinical_presentation : str, optional
        Free-text description of clinical presentation.
    is_acute : bool
        Whether this is an acute metabolic crisis presentation.
    """

    workflow_type = DiagnosticWorkflowType.METABOLIC_SCREENING

    # Reference ranges and critical thresholds for key analytes (umol/L unless noted)
    _ANALYTE_THRESHOLDS: Dict[str, Dict[str, float]] = {
        "phenylalanine": {"upper_normal": 120, "action": 360, "critical": 1200},
        "ammonia": {"upper_normal": 50, "action": 100, "critical": 200},
        "lactate": {"upper_normal": 2.0, "action": 4.0, "critical": 10.0},
        "leucine": {"upper_normal": 200, "action": 500, "critical": 1500},
        "citrulline": {"upper_normal": 45, "action": 100, "critical": 1000},
        "tyrosine": {"upper_normal": 100, "action": 500, "critical": 1000},
        "methionine": {"upper_normal": 40, "action": 100, "critical": 500},
        "homocysteine": {"upper_normal": 15, "action": 50, "critical": 200},
        "propionylcarnitine": {"upper_normal": 5, "action": 10, "critical": 25},
        "octanoylcarnitine": {"upper_normal": 0.3, "action": 0.5, "critical": 2.0},
        "galactose_1_phosphate": {"upper_normal": 1.0, "action": 5.0, "critical": 10.0},
    }

    # Pathway mapping: elevated analyte -> suspected pathway
    _PATHWAY_MAP: Dict[str, Dict] = {
        "phenylalanine": {"pathway": "Phenylalanine hydroxylase", "disease": "Phenylketonuria (PKU)", "gene": "PAH",
                          "confirmatory": ["PAH gene sequencing", "BH4 loading test", "Urine pterin analysis"]},
        "ammonia": {"pathway": "Urea cycle", "disease": "Urea cycle defect", "gene": "OTC/ASS1/ASL/ARG1",
                    "confirmatory": ["Plasma amino acids", "Urine orotic acid", "Urea cycle gene panel"]},
        "lactate": {"pathway": "Mitochondrial/pyruvate metabolism", "disease": "Mitochondrial disorder", "gene": "Multiple",
                    "confirmatory": ["Pyruvate level", "Lactate/pyruvate ratio", "Mitochondrial DNA analysis", "Muscle biopsy"]},
        "leucine": {"pathway": "Branched-chain amino acid catabolism", "disease": "Maple Syrup Urine Disease", "gene": "BCKDHA/BCKDHB/DBT",
                    "confirmatory": ["Alloisoleucine level", "BCKD enzyme activity", "MSUD gene panel"]},
        "citrulline": {"pathway": "Urea cycle (ASS1)", "disease": "Citrullinemia", "gene": "ASS1",
                       "confirmatory": ["ASS1 gene sequencing", "Enzyme activity in fibroblasts"]},
        "tyrosine": {"pathway": "Tyrosine catabolism", "disease": "Tyrosinemia", "gene": "FAH/TAT/HPD",
                     "confirmatory": ["Succinylacetone (urine/blood)", "FAH gene sequencing", "Liver function tests"]},
        "propionylcarnitine": {"pathway": "Propionate metabolism", "disease": "Propionic/Methylmalonic Acidemia", "gene": "PCCA/PCCB/MUT",
                               "confirmatory": ["Urine organic acids", "Acylcarnitine profile", "Gene sequencing"]},
        "octanoylcarnitine": {"pathway": "Fatty acid beta-oxidation", "disease": "MCAD Deficiency", "gene": "ACADM",
                              "confirmatory": ["Urine organic acids", "ACADM gene sequencing", "Acylcarnitine profile"]},
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("analytes"):
            warnings.append("No analyte values provided — screening will be limited")
            inp["analytes"] = {}
        inp.setdefault("age", None)
        inp.setdefault("clinical_presentation", None)
        inp.setdefault("is_acute", False)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        analytes = inputs.get("analytes", {})
        is_acute = inputs.get("is_acute", False)
        severity = SeverityLevel.INFORMATIONAL
        confidence = 0.0

        if not analytes:
            findings.append("No analyte values provided for metabolic screening")
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=findings,
                recommendations=["Submit newborn screening or metabolic panel results"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        findings.append(f"Evaluating {len(analytes)} analyte values")

        # Check each analyte against thresholds
        elevated_analytes: List[Dict] = []
        critical_analytes: List[str] = []

        for analyte_name, value in analytes.items():
            name_lower = analyte_name.lower().replace(" ", "_")
            thresholds = self._ANALYTE_THRESHOLDS.get(name_lower)

            if thresholds is None:
                findings.append(f"  {analyte_name}: {value} (no reference range available)")
                continue

            if value > thresholds["critical"]:
                findings.append(f"  {analyte_name}: {value} — CRITICAL (>{thresholds['critical']})")
                critical_analytes.append(analyte_name)
                elevated_analytes.append({"name": name_lower, "value": value, "level": "critical"})
                severity = _max_severity(severity, SeverityLevel.CRITICAL)
            elif value > thresholds["action"]:
                findings.append(f"  {analyte_name}: {value} — ELEVATED (>{thresholds['action']})")
                elevated_analytes.append({"name": name_lower, "value": value, "level": "elevated"})
                severity = _max_severity(severity, SeverityLevel.HIGH)
            elif value > thresholds["upper_normal"]:
                findings.append(f"  {analyte_name}: {value} — borderline (>{thresholds['upper_normal']})")
                elevated_analytes.append({"name": name_lower, "value": value, "level": "borderline"})
                severity = _max_severity(severity, SeverityLevel.MODERATE)
            else:
                findings.append(f"  {analyte_name}: {value} — normal")

        # Pathway analysis for elevated analytes
        suspected_pathways: List[Dict] = []
        for ea in elevated_analytes:
            pathway_info = self._PATHWAY_MAP.get(ea["name"])
            if pathway_info:
                suspected_pathways.append(pathway_info)
                findings.append(
                    f"Pathway analysis: {ea['name']} elevation suggests "
                    f"{pathway_info['disease']} ({pathway_info['pathway']} pathway)"
                )
                recommendations.append(
                    f"Confirmatory testing for {pathway_info['disease']}: "
                    f"{', '.join(pathway_info['confirmatory'])}"
                )

        # Emergency protocols for metabolic crisis
        if is_acute or critical_analytes:
            ammonia_val = analytes.get("ammonia", analytes.get("Ammonia", 0))
            if isinstance(ammonia_val, (int, float)) and ammonia_val > 200:
                recommendations.insert(0, "EMERGENCY: Hyperammonemia protocol — stop protein intake immediately")
                recommendations.insert(1, "Start IV glucose (D10) at 1.5x maintenance rate")
                recommendations.insert(2, "Administer sodium benzoate/sodium phenylacetate (Ammonul)")
                recommendations.insert(3, "Consider hemodialysis if ammonia >500 umol/L")
                triggers.append(_trigger_string("EMERGENCY", "Hyperammonemia >200 umol/L — metabolic crisis protocol"))

            lactate_val = analytes.get("lactate", analytes.get("Lactate", 0))
            if isinstance(lactate_val, (int, float)) and lactate_val > 10:
                recommendations.insert(0, "EMERGENCY: Severe lactic acidosis — assess for tissue hypoperfusion")
                recommendations.insert(1, "Check blood gas, electrolytes, glucose")
                recommendations.insert(2, "Consider IV bicarbonate if pH <7.1")
                triggers.append(_trigger_string("EMERGENCY", "Lactic acidosis >10 mmol/L"))

        confidence = _clamp(0.2 + 0.15 * len(elevated_analytes) + 0.1 * len(suspected_pathways))

        references.append("ACMG ACT Sheets and Algorithms for newborn screening follow-up")
        references.append("Saudubray et al. (2016) Inborn Metabolic Diseases: Diagnosis and Treatment, 6th ed.")
        references.append("Blau et al. (2014) Physician's Guide to the Diagnosis, Treatment, and Follow-Up of Inherited Metabolic Diseases")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 4 — Dysmorphology Assessment
# ═══════════════════════════════════════════════════════════════════════════


class DysmorphologyWorkflow(BaseRareDiseaseWorkflow):
    """Dysmorphology assessment for syndromic pattern matching.

    Accepts facial features, growth parameters, and skeletal findings,
    then matches against known syndrome patterns to produce ranked
    syndrome candidates.

    Inputs
    ------
    facial_features : list[str]
        Observed facial dysmorphic features.
    growth_parameters : dict
        Growth data (height_percentile, weight_percentile, hc_percentile).
    skeletal_findings : list[str]
        Skeletal abnormalities observed.
    other_findings : list[str], optional
        Additional clinical findings.
    """

    workflow_type = DiagnosticWorkflowType.DYSMORPHOLOGY

    # Syndrome database with characteristic feature sets
    _SYNDROME_DB: List[Dict] = [
        {"name": "Down Syndrome (Trisomy 21)", "id": "OMIM:190685",
         "facial": ["upslanting palpebral fissures", "epicanthal folds", "flat nasal bridge", "protruding tongue", "small ears"],
         "growth": {"short_stature": True, "microcephaly": False},
         "skeletal": ["single palmar crease", "clinodactyly", "sandal gap", "atlantoaxial instability"],
         "other": ["intellectual disability", "congenital heart defect", "hypothyroidism", "hearing loss"],
         "genes": ["Trisomy 21"], "inheritance": InheritancePattern.DE_NOVO, "prevalence": "1:700"},
        {"name": "Turner Syndrome (45,X)", "id": "OMIM:300082",
         "facial": ["webbed neck", "low posterior hairline", "widely spaced nipples"],
         "growth": {"short_stature": True, "microcephaly": False},
         "skeletal": ["cubitus valgus", "short 4th metacarpal", "shield chest"],
         "other": ["primary amenorrhea", "coarctation of aorta", "lymphedema"],
         "genes": ["45,X"], "inheritance": InheritancePattern.DE_NOVO, "prevalence": "1:2,500 females"},
        {"name": "Noonan Syndrome", "id": "OMIM:163950",
         "facial": ["hypertelorism", "ptosis", "low-set ears", "webbed neck", "downslanting palpebral fissures"],
         "growth": {"short_stature": True, "microcephaly": False},
         "skeletal": ["pectus excavatum", "cubitus valgus"],
         "other": ["pulmonary stenosis", "bleeding diathesis", "cryptorchidism", "developmental delay"],
         "genes": ["PTPN11", "SOS1", "RAF1", "KRAS"], "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "prevalence": "1:2,500"},
        {"name": "Williams Syndrome", "id": "OMIM:194050",
         "facial": ["broad forehead", "stellate iris", "long philtrum", "wide mouth", "full lips"],
         "growth": {"short_stature": True, "microcephaly": False},
         "skeletal": ["joint laxity", "kyphoscoliosis"],
         "other": ["supravalvular aortic stenosis", "hypercalcemia", "intellectual disability", "overfriendly personality"],
         "genes": ["ELN", "7q11.23 deletion"], "inheritance": InheritancePattern.DE_NOVO, "prevalence": "1:7,500"},
        {"name": "Cornelia de Lange Syndrome", "id": "OMIM:122470",
         "facial": ["synophrys", "long eyelashes", "thin upper lip", "anteverted nares", "micrognathia"],
         "growth": {"short_stature": True, "microcephaly": True},
         "skeletal": ["upper limb reduction", "small hands", "clinodactyly", "proximally placed thumbs"],
         "other": ["intellectual disability", "hirsutism", "GERD", "hearing loss"],
         "genes": ["NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8"], "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "prevalence": "1:10,000"},
        {"name": "Treacher Collins Syndrome", "id": "OMIM:154500",
         "facial": ["malar hypoplasia", "downslanting palpebral fissures", "micrognathia", "ear anomalies", "coloboma of lower eyelid"],
         "growth": {"short_stature": False, "microcephaly": False},
         "skeletal": ["zygomatic arch hypoplasia", "mandibular hypoplasia"],
         "other": ["conductive hearing loss", "choanal atresia", "cleft palate"],
         "genes": ["TCOF1", "POLR1C", "POLR1D"], "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "prevalence": "1:50,000"},
        {"name": "Kabuki Syndrome", "id": "OMIM:147920",
         "facial": ["long palpebral fissures", "arched eyebrows", "eversion of lower eyelid", "broad nose", "prominent ears"],
         "growth": {"short_stature": True, "microcephaly": False},
         "skeletal": ["persistent fetal finger pads", "brachydactyly", "scoliosis", "joint laxity"],
         "other": ["intellectual disability", "congenital heart defect", "cleft palate", "renal anomaly"],
         "genes": ["KMT2D", "KDM6A"], "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "prevalence": "1:32,000"},
        {"name": "Sotos Syndrome", "id": "OMIM:117550",
         "facial": ["macrocephaly", "frontal bossing", "pointed chin", "downslanting palpebral fissures", "high hairline"],
         "growth": {"short_stature": False, "macrocephaly": True, "overgrowth": True},
         "skeletal": ["advanced bone age", "large hands and feet", "scoliosis"],
         "other": ["intellectual disability", "neonatal hypotonia", "seizures", "cardiac anomaly"],
         "genes": ["NSD1"], "inheritance": InheritancePattern.AUTOSOMAL_DOMINANT, "prevalence": "1:14,000"},
    ]

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("facial_features") and not inp.get("skeletal_findings"):
            warnings.append("No dysmorphic features provided — assessment will be limited")
        inp.setdefault("facial_features", [])
        inp.setdefault("growth_parameters", {})
        inp.setdefault("skeletal_findings", [])
        inp.setdefault("other_findings", [])
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        candidates: List[DiseaseCandidate] = []

        facial = [f.lower() for f in inputs.get("facial_features", [])]
        skeletal = [f.lower() for f in inputs.get("skeletal_findings", [])]
        other = [f.lower() for f in inputs.get("other_findings", [])]
        growth = inputs.get("growth_parameters", {})
        all_features = set(facial + skeletal + other)

        if not all_features:
            findings.append("No dysmorphic features provided for assessment")
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=findings,
                recommendations=["Provide facial features, skeletal findings, or other clinical features"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        findings.append(f"Evaluating {len(all_features)} clinical features against {len(self._SYNDROME_DB)} syndromes")

        # Score each syndrome
        scored: List[Dict] = []
        for syndrome in self._SYNDROME_DB:
            syndrome_features = set(
                [f.lower() for f in syndrome["facial"]]
                + [f.lower() for f in syndrome["skeletal"]]
                + [f.lower() for f in syndrome.get("other", [])]
            )

            # Partial string matching (feature substring in syndrome feature or vice versa)
            matched = []
            for patient_f in all_features:
                for syn_f in syndrome_features:
                    if patient_f in syn_f or syn_f in patient_f:
                        matched.append(patient_f)
                        break

            if not matched:
                continue

            unmatched = [f for f in syndrome_features if not any(
                pf in f or f in pf for pf in all_features
            )]

            # Score: matched fraction with weight for facial features
            total_syn_features = len(syndrome_features)
            match_score = len(matched) / total_syn_features if total_syn_features else 0.0

            # Growth parameter bonus
            growth_bonus = 0.0
            if growth.get("height_percentile", 50) < 3 and syndrome["growth"].get("short_stature"):
                growth_bonus = 0.05
            if growth.get("hc_percentile", 50) < 3 and syndrome["growth"].get("microcephaly"):
                growth_bonus = 0.05
            if growth.get("hc_percentile", 50) > 97 and syndrome["growth"].get("macrocephaly"):
                growth_bonus = 0.05

            total_score = _clamp(match_score + growth_bonus)

            scored.append({
                "syndrome": syndrome,
                "matched": matched,
                "unmatched": unmatched,
                "score": total_score,
            })

        scored.sort(key=lambda x: x["score"], reverse=True)
        scored = scored[:10]

        for rank, entry in enumerate(scored, 1):
            s = entry["syndrome"]
            genes = s.get("genes", [])
            candidates.append(DiseaseCandidate(
                disease_id=s["id"],
                disease_name=s["name"],
                rank=rank,
                similarity_score=round(entry["score"], 4),
                matched_phenotypes=entry["matched"],
                unmatched_phenotypes=entry["unmatched"],
                inheritance_pattern=s.get("inheritance"),
                prevalence=s.get("prevalence"),
                causal_genes=genes,
            ))
            findings.append(
                f"  #{rank} {s['name']} (score={entry['score']:.3f}, "
                f"matched={len(entry['matched'])})"
            )

        if candidates:
            top = candidates[0]
            recommendations.append(f"Top syndromic match: {top.disease_name}")
            if top.causal_genes and top.causal_genes[0] not in ("Trisomy 21", "45,X"):
                recommendations.append(
                    f"Recommend genetic testing: {', '.join(top.causal_genes)}"
                )
            else:
                recommendations.append("Recommend chromosomal analysis (karyotype/CMA)")

        confidence = candidates[0].similarity_score if candidates else 0.0
        severity = SeverityLevel.MODERATE if confidence > 0.2 else SeverityLevel.LOW

        references.append("Jones et al. (2021) Smith's Recognizable Patterns of Human Malformation, 8th ed.")
        references.append("GeneReviews - NCBI Bookshelf (https://www.ncbi.nlm.nih.gov/books/NBK1116/)")

        diagnostic = DiagnosticResult(
            candidate_diseases=candidates,
            confidence=confidence,
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            confidence=confidence,
            diagnostic_result=diagnostic,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 5 — Neurogenetic Evaluation
# ═══════════════════════════════════════════════════════════════════════════


class NeurogeneticWorkflow(BaseRareDiseaseWorkflow):
    """Neurogenetic evaluation for developmental delay, epilepsy, and
    movement disorders.

    Analyses DD/ID + seizure type + movement disorder presentation to
    recommend appropriate gene panels and testing algorithms
    (epilepsy: SCN1A/CDKL5/STXBP1, DD: CMA -> WES -> WGS).

    Inputs
    ------
    presentation : str
        Primary neurological presentation (e.g., 'epilepsy', 'developmental_delay',
        'movement_disorder', 'intellectual_disability').
    seizure_types : list[str], optional
        Types of seizures observed (focal, generalized, absence, myoclonic, etc.).
    age_of_onset : str, optional
        Age at symptom onset.
    movement_findings : list[str], optional
        Movement abnormalities (dystonia, ataxia, chorea, spasticity, etc.).
    mri_findings : str, optional
        Brain MRI findings.
    eeg_findings : str, optional
        EEG findings.
    prior_testing : list[str], optional
        Prior genetic tests already performed.
    """

    workflow_type = DiagnosticWorkflowType.NEUROGENETIC

    # Gene panels by presentation
    _EPILEPSY_PANELS: Dict[str, List[str]] = {
        "infantile_spasms": ["ARX", "CDKL5", "STXBP1", "KCNQ2", "SCN2A", "SCN8A", "KCNT1", "PIGA", "FOXG1", "TSC1", "TSC2"],
        "dravet_like": ["SCN1A", "SCN2A", "SCN8A", "PCDH19", "GABRA1", "GABRG2", "CHD2", "HCN1", "STXBP1"],
        "progressive_myoclonic": ["CSTB", "EPM2A", "NHLRC1", "CLN3", "CLN5", "CLN6", "CLN8", "PPT1", "TPP1", "KCTD7"],
        "absence": ["GABRA1", "GABRG2", "GABRB3", "CACNA1A", "SLC2A1"],
        "focal": ["DEPDC5", "NPRL2", "NPRL3", "LGI1", "KCNT1", "SCN1A", "CHRNA4", "CHRNB2"],
        "general": ["SCN1A", "CDKL5", "STXBP1", "KCNQ2", "SCN2A", "MECP2", "ARX", "FOXG1", "SLC2A1", "TSC1", "TSC2"],
    }

    _DD_ALGORITHM: List[Dict[str, str]] = [
        {"step": 1, "test": "Chromosomal microarray (CMA)", "yield": "15-20%",
         "rationale": "First-tier for DD/ID — detects CNVs missed by karyotype"},
        {"step": 2, "test": "Fragile X testing (FMR1)", "yield": "2-3%",
         "rationale": "Common single-gene cause of DD/ID, especially in males"},
        {"step": 3, "test": "Whole Exome Sequencing (WES)", "yield": "25-40%",
         "rationale": "Broad gene-level detection for Mendelian disorders"},
        {"step": 4, "test": "Whole Genome Sequencing (WGS)", "yield": "Additional 5-10%",
         "rationale": "Detects non-coding, deep intronic, and structural variants missed by WES"},
        {"step": 5, "test": "Metabolic workup", "yield": "1-5%",
         "rationale": "Treatable metabolic causes should not be missed"},
    ]

    _MOVEMENT_GENES: Dict[str, List[str]] = {
        "dystonia": ["TOR1A", "THAP1", "GCH1", "TH", "DDC", "SLC6A3", "SGCE", "ATP1A3", "KMT2B", "GNAO1"],
        "ataxia": ["ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "FXN", "APTX", "SETX", "ATM", "SACS"],
        "chorea": ["HTT", "NKX2-1", "ADCY5", "PDE10A", "GNAO1"],
        "spasticity": ["SPG4", "SPG3A", "SPG7", "SPG11", "SPG5", "ATL1", "REEP1", "KIF5A"],
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("presentation"):
            warnings.append("No neurological presentation specified — defaulting to general evaluation")
            inp["presentation"] = "general"
        inp.setdefault("seizure_types", [])
        inp.setdefault("age_of_onset", None)
        inp.setdefault("movement_findings", [])
        inp.setdefault("mri_findings", None)
        inp.setdefault("eeg_findings", None)
        inp.setdefault("prior_testing", [])
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        presentation = inputs["presentation"].lower()
        seizure_types = [s.lower() for s in inputs.get("seizure_types", [])]
        movement_findings = [m.lower() for m in inputs.get("movement_findings", [])]
        prior_testing = [t.lower() for t in inputs.get("prior_testing", [])]
        age_onset = inputs.get("age_of_onset", "")
        mri_findings = inputs.get("mri_findings", "")
        eeg_findings = inputs.get("eeg_findings", "")

        findings.append(f"Neurogenetic evaluation for: {presentation}")

        # Epilepsy pathway
        if "epilepsy" in presentation or "seizure" in presentation or seizure_types:
            findings.append("Epilepsy/seizure pathway activated")

            # Determine panel based on seizure type and age
            panel_key = "general"
            if any("spasm" in s for s in seizure_types) or (age_onset and "infant" in age_onset.lower()):
                panel_key = "infantile_spasms"
            elif any("myoclonic" in s for s in seizure_types):
                if age_onset and any(kw in age_onset.lower() for kw in ["teen", "adolesc", "juvenile"]):
                    panel_key = "progressive_myoclonic"
                else:
                    panel_key = "dravet_like"
            elif any("absence" in s for s in seizure_types):
                panel_key = "absence"
            elif any("focal" in s for s in seizure_types):
                panel_key = "focal"

            genes = self._EPILEPSY_PANELS.get(panel_key, self._EPILEPSY_PANELS["general"])
            findings.append(f"Recommended epilepsy gene panel ({panel_key}): {len(genes)} genes")
            findings.append(f"  Key genes: {', '.join(genes[:5])}")
            recommendations.append(f"Order {panel_key} epilepsy gene panel: {', '.join(genes)}")

            if age_onset and any(kw in age_onset.lower() for kw in ["neonat", "newborn", "first week"]):
                findings.append("Neonatal seizure onset — high priority for rapid WGS")
                recommendations.insert(0, "URGENT: Consider rapid WGS for neonatal seizures (turnaround <2 weeks)")

            if eeg_findings:
                findings.append(f"EEG findings noted: {eeg_findings}")

        # Developmental delay pathway
        if "delay" in presentation or "dd" in presentation or "id" in presentation or "intellectual" in presentation:
            findings.append("Developmental delay / intellectual disability pathway activated")

            # Recommend step-wise algorithm, skipping already-completed tests
            for step in self._DD_ALGORITHM:
                test_lower = step["test"].lower()
                already_done = any(t in test_lower or test_lower in t for t in prior_testing)
                status = "COMPLETED" if already_done else "RECOMMENDED"
                findings.append(f"  Step {step['step']}: {step['test']} (yield: {step['yield']}) — {status}")
                if not already_done:
                    recommendations.append(f"Step {step['step']}: {step['test']} (expected yield: {step['yield']})")

        # Movement disorder pathway
        if "movement" in presentation or movement_findings:
            findings.append("Movement disorder pathway activated")
            recommended_genes: List[str] = []
            for finding in movement_findings:
                for disorder_type, genes in self._MOVEMENT_GENES.items():
                    if disorder_type in finding or finding in disorder_type:
                        findings.append(f"  {disorder_type} genes: {', '.join(genes[:5])}")
                        recommended_genes.extend(genes)
                        break

            if recommended_genes:
                unique_genes = list(dict.fromkeys(recommended_genes))
                recommendations.append(f"Movement disorder gene panel: {', '.join(unique_genes[:15])}")
            else:
                recommendations.append("Consider comprehensive movement disorder gene panel")

        if mri_findings:
            findings.append(f"MRI findings: {mri_findings}")
            if any(kw in mri_findings.lower() for kw in ["white matter", "leukodystrophy", "demyelination"]):
                recommendations.append("Consider leukodystrophy gene panel (MLC1, GFAP, EIF2B1-5, GALC, ARSA)")
                triggers.append(_trigger_string("NEUROIMAGING", "White matter abnormalities — leukodystrophy workup"))
            if any(kw in mri_findings.lower() for kw in ["cerebellar", "atrophy", "hypoplasia"]):
                recommendations.append("Consider cerebellar ataxia gene panel")

        confidence = _clamp(0.3 + 0.05 * len(findings))
        severity = SeverityLevel.HIGH if "URGENT" in str(recommendations) else SeverityLevel.MODERATE

        references.append("Scheffer et al. (2017) ILAE classification of the epilepsies. Epilepsia")
        references.append("Michelson et al. (2011) Evidence report: Genetic and metabolic testing on children with GDD. Neurology")
        references.append("ACMG practice guideline: Evaluation of the child with DD/ID (2021)")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 6 — Cardiac Genetics
# ═══════════════════════════════════════════════════════════════════════════


class CardiacGeneticsWorkflow(BaseRareDiseaseWorkflow):
    """Cardiac genetics workflow for arrhythmia and cardiomyopathy phenotypes.

    Maps cardiac phenotype to specific gene panels (channelopathy vs
    sarcomeric vs desmosomal) and triggers cross-modal referral to
    Cardiology Intelligence Agent.

    Inputs
    ------
    cardiac_phenotype : str
        Primary cardiac diagnosis (e.g., 'long_qt', 'hcm', 'dcm', 'arvc', 'brugada').
    ecg_findings : str, optional
        ECG findings.
    echo_findings : str, optional
        Echocardiographic findings.
    family_history : str, optional
        Cardiac family history (SCD, cardiomyopathy, etc.).
    age : str, optional
        Patient age.
    symptoms : list[str], optional
        Cardiac symptoms (syncope, palpitations, chest pain, etc.).
    """

    workflow_type = DiagnosticWorkflowType.CARDIAC_GENETICS

    # Gene panels by cardiac phenotype
    _CARDIAC_PANELS: Dict[str, Dict] = {
        "long_qt": {
            "genes": ["KCNQ1", "KCNH2", "SCN5A", "KCNE1", "KCNE2", "KCNJ2", "CACNA1C", "CAV3", "SCN4B",
                       "AKAP9", "ANK2", "SNTA1", "CALM1", "CALM2", "CALM3", "TRDN"],
            "category": "channelopathy",
            "management": ["Beta-blocker therapy (nadolol preferred)", "Avoid QT-prolonging drugs",
                           "ICD consideration if high-risk (LQT3, cardiac arrest survivor)"],
        },
        "brugada": {
            "genes": ["SCN5A", "CACNA1C", "CACNB2", "CACNA2D1", "GPD1L", "KCNE3",
                       "KCND3", "HCN4", "SCN1B", "SCN2B", "SCN3B", "TRPM4"],
            "category": "channelopathy",
            "management": ["Avoid drugs on BrugadaDrugs.org", "Fever management (antipyretics)",
                           "ICD if symptomatic or inducible VF on EPS"],
        },
        "cpvt": {
            "genes": ["RYR2", "CASQ2", "TRDN", "CALM1", "CALM2", "CALM3", "TECRL"],
            "category": "channelopathy",
            "management": ["Flecainide + beta-blocker", "Exercise restriction",
                           "ICD if recurrent events on medical therapy"],
        },
        "hcm": {
            "genes": ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL2", "MYL3",
                       "CSRP3", "TNNC1", "MYOZ2", "JPH2", "PLN"],
            "category": "sarcomeric",
            "management": ["Serial echocardiography", "Holter monitoring for arrhythmia",
                           "SCD risk stratification (ESC/AHA)", "Mavacamten consideration for obstructive HCM"],
        },
        "dcm": {
            "genes": ["TTN", "LMNA", "MYH7", "TNNT2", "SCN5A", "RBM20", "BAG3", "DES",
                       "FLNC", "PLN", "TNNC1", "TPM1", "DSP"],
            "category": "sarcomeric_structural",
            "management": ["GDMT for heart failure", "ICD if LVEF <35%",
                           "LMNA carriers: early ICD consideration even with preserved EF"],
        },
        "arvc": {
            "genes": ["PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43", "PLN", "DES",
                       "CDH2", "CTNNA3", "FLNC"],
            "category": "desmosomal",
            "management": ["Exercise restriction (competitive sports prohibited)",
                           "Holter monitoring", "ICD if high-risk features",
                           "Annual cardiac MRI for disease progression"],
        },
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("cardiac_phenotype"):
            warnings.append("No cardiac phenotype specified — using general cardiac panel")
            inp["cardiac_phenotype"] = "general"
        inp.setdefault("ecg_findings", None)
        inp.setdefault("echo_findings", None)
        inp.setdefault("family_history", None)
        inp.setdefault("age", None)
        inp.setdefault("symptoms", [])
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        phenotype = inputs["cardiac_phenotype"].lower().replace(" ", "_").replace("-", "_")
        symptoms = inputs.get("symptoms", [])
        family_history = inputs.get("family_history", "")
        ecg_findings = inputs.get("ecg_findings", "")
        echo_findings = inputs.get("echo_findings", "")

        findings.append(f"Cardiac genetics evaluation for: {phenotype}")

        # Look up specific panel
        panel_info = self._CARDIAC_PANELS.get(phenotype)

        if panel_info:
            genes = panel_info["genes"]
            category = panel_info["category"]
            management = panel_info["management"]

            findings.append(f"Gene panel category: {category}")
            findings.append(f"Panel includes {len(genes)} genes: {', '.join(genes[:8])}...")
            recommendations.append(f"Order {category} gene panel ({len(genes)} genes) for {phenotype}")

            for mgmt in management:
                recommendations.append(f"Management: {mgmt}")
        else:
            # General comprehensive cardiac panel
            all_genes: List[str] = []
            for panel in self._CARDIAC_PANELS.values():
                all_genes.extend(panel["genes"])
            unique_genes = list(dict.fromkeys(all_genes))
            findings.append(f"No specific panel match — recommending comprehensive cardiac panel ({len(unique_genes)} genes)")
            recommendations.append(f"Order comprehensive cardiac gene panel ({len(unique_genes)} genes)")

        # Symptom analysis
        high_risk_symptoms = ["syncope", "cardiac arrest", "sudden death", "aborted scd"]
        if any(s.lower() in high_risk_symptoms for s in symptoms):
            findings.append("HIGH RISK: Syncope or cardiac arrest history")
            recommendations.insert(0, "URGENT: Immediate cardiology evaluation and risk stratification")
            recommendations.insert(1, "Consider ICD evaluation")

        # Family history
        if family_history:
            findings.append(f"Family history: {family_history}")
            if any(kw in family_history.lower() for kw in ["sudden death", "scd", "died suddenly"]):
                findings.append("Family history of sudden cardiac death identified")
                recommendations.append("Cascade screening of first-degree relatives recommended")

        if ecg_findings:
            findings.append(f"ECG findings: {ecg_findings}")
        if echo_findings:
            findings.append(f"Echocardiographic findings: {echo_findings}")

        # Cross-modal trigger to Cardiology Agent
        triggers.append(_trigger_string("CARDIOLOGY", f"Cardiac genetics referral for {phenotype} — panel testing initiated"))

        confidence = _clamp(0.4 + (0.1 if panel_info else 0.0) + 0.05 * len(findings))
        severity = SeverityLevel.HIGH if any(s.lower() in high_risk_symptoms for s in symptoms) else SeverityLevel.MODERATE

        references.append("Wilde et al. (2022) ESC Guidelines on ventricular arrhythmias and SCD. Eur Heart J")
        references.append("Ommen et al. (2020) AHA/ACC Guideline for HCM. Circulation")
        references.append("Towbin et al. (2019) ARVC diagnosis: HRS expert consensus. Heart Rhythm")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 7 — Connective Tissue Disorders
# ═══════════════════════════════════════════════════════════════════════════


class ConnectiveTissueWorkflow(BaseRareDiseaseWorkflow):
    """Connective tissue disorder assessment applying diagnostic criteria.

    Evaluates Beighton score for hypermobility, Ghent criteria for Marfan
    syndrome, and Villefranche/2017 criteria for Ehlers-Danlos syndrome.

    Inputs
    ------
    beighton_score : int, optional
        Beighton hypermobility score (0-9).
    ghent_criteria : dict, optional
        Ghent criteria features (aortic_root_z, ectopia_lentis, systemic_score, fbn1_mutation, family_history).
    eds_features : dict, optional
        EDS features (skin_hyperextensibility, joint_hypermobility, skin_fragility, etc.).
    eds_subtype : str, optional
        Suspected EDS subtype.
    skeletal_features : list[str], optional
        Skeletal features observed.
    vascular_features : list[str], optional
        Vascular features (aneurysm, dissection, etc.).
    """

    workflow_type = DiagnosticWorkflowType.CONNECTIVE_TISSUE

    # EDS gene panels by subtype
    _EDS_PANELS: Dict[str, Dict] = {
        "classical": {"genes": ["COL5A1", "COL5A2", "COL1A1"], "criteria": "Major: skin hyperextensibility + atrophic scarring, GJH (Beighton >= 5)"},
        "hypermobility": {"genes": ["TNXB"], "criteria": "2017 criteria: GJH + 2 of 3 features (systemic, family history, MSK complications)"},
        "vascular": {"genes": ["COL3A1", "COL1A1"], "criteria": "Major: arterial/organ rupture, thin skin, characteristic facies. URGENT genetic testing."},
        "kyphoscoliotic": {"genes": ["PLOD1", "FKBP14"], "criteria": "Major: congenital hypotonia, congenital/early kyphoscoliosis, GJH"},
        "arthrochalasia": {"genes": ["COL1A1", "COL1A2"], "criteria": "Major: congenital bilateral hip dislocation, severe GJH"},
        "dermatosparaxis": {"genes": ["ADAMTS2"], "criteria": "Major: extreme skin fragility, redundant skin"},
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        inp.setdefault("beighton_score", None)
        inp.setdefault("ghent_criteria", None)
        inp.setdefault("eds_features", None)
        inp.setdefault("eds_subtype", None)
        inp.setdefault("skeletal_features", [])
        inp.setdefault("vascular_features", [])

        if inp["beighton_score"] is not None:
            try:
                score = int(inp["beighton_score"])
                if score < 0 or score > 9:
                    warnings.append(f"Beighton score {score} out of range (0-9) — clamping")
                    inp["beighton_score"] = max(0, min(9, score))
            except (ValueError, TypeError):
                warnings.append(f"Invalid Beighton score '{inp['beighton_score']}' — ignoring")
                inp["beighton_score"] = None
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        beighton = inputs.get("beighton_score")
        ghent = inputs.get("ghent_criteria")
        eds_features = inputs.get("eds_features", {})
        eds_subtype = inputs.get("eds_subtype")
        vascular = inputs.get("vascular_features", [])

        # Beighton assessment
        if beighton is not None:
            findings.append(f"Beighton score: {beighton}/9")
            if beighton >= 5:
                findings.append("Generalised joint hypermobility (GJH) confirmed (Beighton >= 5)")
            elif beighton >= 4:
                findings.append("Borderline joint hypermobility (Beighton 4/9)")
            else:
                findings.append("Joint hypermobility not confirmed by Beighton criteria")

        # Ghent criteria for Marfan
        if ghent:
            findings.append("Evaluating revised Ghent criteria for Marfan syndrome")
            aortic_z = ghent.get("aortic_root_z", 0)
            ectopia_lentis = ghent.get("ectopia_lentis", False)
            systemic_score = ghent.get("systemic_score", 0)
            fbn1_mutation = ghent.get("fbn1_mutation", False)
            ghent.get("family_history", False)

            marfan_met = False

            # Criterion A: Aortic root Z >= 2 + ectopia lentis
            if aortic_z >= 2 and ectopia_lentis:
                findings.append("Ghent criteria MET: Aortic root dilation (Z >= 2) + ectopia lentis")
                marfan_met = True
            # Criterion B: Aortic root Z >= 2 + FBN1 mutation
            elif aortic_z >= 2 and fbn1_mutation:
                findings.append("Ghent criteria MET: Aortic root dilation (Z >= 2) + FBN1 mutation")
                marfan_met = True
            # Criterion C: Ectopia lentis + FBN1 known to cause aortic disease
            elif ectopia_lentis and fbn1_mutation:
                findings.append("Ghent criteria MET: Ectopia lentis + pathogenic FBN1 variant")
                marfan_met = True
            # Criterion D: Systemic score >= 7 + aortic root Z >= 2
            elif systemic_score >= 7 and aortic_z >= 2:
                findings.append("Ghent criteria MET: Systemic score >= 7 + aortic root dilation")
                marfan_met = True
            else:
                findings.append(
                    f"Ghent criteria NOT fully met (aortic Z={aortic_z}, "
                    f"ectopia_lentis={ectopia_lentis}, systemic={systemic_score}, "
                    f"FBN1={fbn1_mutation})"
                )

            if marfan_met:
                recommendations.append("Diagnosis of Marfan syndrome supported by Ghent criteria")
                recommendations.append("Management: annual aortic root imaging, beta-blocker/ARB therapy")
                recommendations.append("Ophthalmology referral for ectopia lentis monitoring")
                recommendations.append("Activity restriction: avoid contact sports, isometric exercise")
                if not fbn1_mutation:
                    recommendations.append("Confirm with FBN1 gene sequencing if not yet performed")
                triggers.append(_trigger_string("CARDIOLOGY", "Marfan syndrome — aortic root monitoring"))
            else:
                recommendations.append("Consider FBN1 gene sequencing to support/refute Marfan diagnosis")

        # EDS evaluation
        if eds_subtype or eds_features:
            subtype = (eds_subtype or "hypermobility").lower()
            panel_info = self._EDS_PANELS.get(subtype, self._EDS_PANELS.get("hypermobility"))

            findings.append(f"EDS subtype evaluation: {subtype}")
            if panel_info:
                findings.append(f"Diagnostic criteria: {panel_info['criteria']}")
                recommendations.append(f"Gene panel for EDS-{subtype}: {', '.join(panel_info['genes'])}")

            # Vascular EDS is urgent
            if subtype == "vascular" or vascular:
                findings.append("URGENT: Vascular EDS or vascular features identified")
                recommendations.insert(0, "URGENT: COL3A1 gene sequencing — vascular EDS carries high mortality risk")
                recommendations.insert(1, "Avoid invasive procedures, arteriography, colonoscopy unless essential")
                triggers.append(_trigger_string("VASCULAR", "Possible vascular EDS — avoid invasive procedures"))

        if not ghent and not eds_features and beighton is None and not eds_subtype:
            findings.append("Insufficient data for connective tissue disorder assessment")
            recommendations.append("Provide Beighton score, Ghent criteria, or EDS features for evaluation")

        confidence = _clamp(0.3 + 0.1 * len(findings))
        severity_level = SeverityLevel.CRITICAL if "vascular" in str(findings).lower() else SeverityLevel.MODERATE

        references.append("Loeys et al. (2010) Revised Ghent nosology for Marfan syndrome. J Med Genet")
        references.append("Malfait et al. (2017) International EDS classification. Am J Med Genet C")
        references.append("ACMG Practice Resource: Marfan syndrome and related disorders (2021)")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity_level,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 8 — Inborn Errors of Metabolism
# ═══════════════════════════════════════════════════════════════════════════


class InbornErrorsWorkflow(BaseRareDiseaseWorkflow):
    """Inborn errors of metabolism workflow for enzyme assay and biomarker
    interpretation with dietary management recommendations.

    Includes emergency protocols: stop protein, glucose infusion, carnitine
    supplementation for acute decompensation.

    Inputs
    ------
    enzyme_assays : dict, optional
        Enzyme name -> activity level dict.
    biomarkers : dict, optional
        Biomarker name -> value dict.
    suspected_pathway : str, optional
        Suspected metabolic pathway (amino_acid, organic_acid, fatty_acid, lysosomal, etc.).
    is_acute_decompensation : bool
        Whether patient is in acute metabolic decompensation.
    current_diet : str, optional
        Current dietary restrictions.
    """

    workflow_type = DiagnosticWorkflowType.INBORN_ERRORS

    # Enzyme deficiency patterns
    _ENZYME_DEFICIENCIES: Dict[str, Dict] = {
        "glucocerebrosidase": {"disease": "Gaucher Disease", "gene": "GBA1", "treatment": "Enzyme replacement therapy (imiglucerase, velaglucerase)"},
        "alpha_galactosidase": {"disease": "Fabry Disease", "gene": "GLA", "treatment": "ERT (agalsidase beta) or oral chaperone (migalastat)"},
        "acid_sphingomyelinase": {"disease": "Niemann-Pick A/B", "gene": "SMPD1", "treatment": "Olipudase alfa (for type B)"},
        "hexosaminidase_a": {"disease": "Tay-Sachs Disease", "gene": "HEXA", "treatment": "Supportive care only (no specific treatment)"},
        "galactose_1_phosphate_uridylyltransferase": {"disease": "Classic Galactosemia", "gene": "GALT", "treatment": "Galactose-free diet"},
        "phenylalanine_hydroxylase": {"disease": "Phenylketonuria", "gene": "PAH", "treatment": "Low-phenylalanine diet, sapropterin (BH4-responsive)"},
        "branched_chain_ketoacid_dehydrogenase": {"disease": "Maple Syrup Urine Disease", "gene": "BCKDHA/BCKDHB/DBT", "treatment": "BCAA-restricted diet, thiamine trial"},
        "biotinidase": {"disease": "Biotinidase Deficiency", "gene": "BTD", "treatment": "Oral biotin supplementation (5-20 mg/day)"},
        "medium_chain_acyl_coa_dehydrogenase": {"disease": "MCAD Deficiency", "gene": "ACADM", "treatment": "Avoid fasting >8-12h, emergency glucose protocol"},
    }

    # Dietary management by pathway
    _DIETARY_MANAGEMENT: Dict[str, List[str]] = {
        "amino_acid": [
            "Restrict specific amino acid(s) based on diagnosis",
            "Monitor plasma amino acid levels regularly",
            "Supplement with medical formula (e.g., PKU formula)",
            "Ensure adequate total protein from medical food",
        ],
        "organic_acid": [
            "Restrict precursor amino acids (isoleucine, valine, methionine, threonine)",
            "Supplement L-carnitine (50-100 mg/kg/day)",
            "Provide sick-day emergency protocol to family",
            "Monitor urine organic acids during illness",
        ],
        "fatty_acid": [
            "Avoid fasting (age-appropriate fasting tolerance)",
            "Provide cornstarch at bedtime for young children",
            "Medium-chain triglyceride (MCT) oil for LCHAD/VLCAD",
            "Emergency IV glucose protocol for illness",
        ],
        "urea_cycle": [
            "Restrict protein to prescribed amount",
            "Nitrogen scavenger therapy (sodium benzoate, sodium phenylbutyrate, glycerol phenylbutyrate)",
            "Essential amino acid supplementation",
            "Emergency ammonia protocol for family",
        ],
        "lysosomal": [
            "Enzyme replacement therapy (ERT) as indicated",
            "Consider substrate reduction therapy (SRT)",
            "Monitor for organ-specific complications",
            "Consider hematopoietic stem cell transplant for selected disorders",
        ],
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        inp.setdefault("enzyme_assays", {})
        inp.setdefault("biomarkers", {})
        inp.setdefault("suspected_pathway", None)
        inp.setdefault("is_acute_decompensation", False)
        inp.setdefault("current_diet", None)

        if not inp.get("enzyme_assays") and not inp.get("biomarkers") and not inp.get("suspected_pathway"):
            warnings.append("No enzyme assays, biomarkers, or suspected pathway provided")
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        enzyme_assays = inputs.get("enzyme_assays", {})
        biomarkers = inputs.get("biomarkers", {})
        pathway = inputs.get("suspected_pathway", "")
        is_acute = inputs.get("is_acute_decompensation", False)
        severity = SeverityLevel.INFORMATIONAL

        # EMERGENCY: Acute decompensation protocol
        if is_acute:
            severity = SeverityLevel.CRITICAL
            findings.append("ACUTE METABOLIC DECOMPENSATION — emergency protocol activated")
            recommendations.insert(0, "EMERGENCY: Stop all protein intake immediately")
            recommendations.insert(1, "Start IV D10 at 1.5x maintenance rate (prevent catabolism)")
            recommendations.insert(2, "Administer L-carnitine 100 mg/kg IV (max 6g)")
            recommendations.insert(3, "Draw STAT: ammonia, lactate, glucose, blood gas, electrolytes")
            recommendations.insert(4, "Draw: plasma amino acids, acylcarnitine profile, urine organic acids")
            recommendations.insert(5, "Contact metabolic specialist on-call immediately")
            triggers.append(_trigger_string("EMERGENCY", "Acute metabolic decompensation — PICU/NICU admission"))

        # Enzyme assay interpretation
        for enzyme_name, activity in enzyme_assays.items():
            enzyme_key = enzyme_name.lower().replace(" ", "_").replace("-", "_")
            deficiency_info = self._ENZYME_DEFICIENCIES.get(enzyme_key)

            if deficiency_info:
                findings.append(f"Enzyme: {enzyme_name} — activity: {activity}")
                if isinstance(activity, (int, float)) and activity < 10:
                    findings.append(f"  DEFICIENT: consistent with {deficiency_info['disease']}")
                    recommendations.append(f"Confirmed {deficiency_info['disease']} — {deficiency_info['treatment']}")
                    recommendations.append(f"Confirm with {deficiency_info['gene']} gene sequencing")
                    severity = _max_severity(severity, SeverityLevel.HIGH)
                elif isinstance(activity, (int, float)) and activity < 30:
                    findings.append("  REDUCED: possible carrier or partial deficiency")
                    recommendations.append(f"Consider {deficiency_info['gene']} gene sequencing")
                    severity = _max_severity(severity, SeverityLevel.MODERATE)
            else:
                findings.append(f"Enzyme: {enzyme_name} — activity: {activity} (no reference pattern)")

        # Biomarker interpretation
        for biomarker, value in biomarkers.items():
            findings.append(f"Biomarker: {biomarker} = {value}")

        # Dietary management recommendations
        if pathway:
            pathway_key = pathway.lower().replace(" ", "_")
            diet_recs = self._DIETARY_MANAGEMENT.get(pathway_key, [])
            if diet_recs:
                findings.append(f"Dietary management for {pathway_key} pathway disorder:")
                for rec in diet_recs:
                    recommendations.append(f"Diet: {rec}")
            else:
                recommendations.append(f"Consult metabolic dietitian for {pathway} pathway management")

        if not enzyme_assays and not biomarkers and not pathway and not is_acute:
            findings.append("No enzyme assays, biomarkers, or pathway specified")
            recommendations.append("Provide enzyme assay results, biomarkers, or suspected pathway")

        confidence = _clamp(0.2 + 0.1 * len(enzyme_assays) + 0.05 * len(biomarkers))

        references.append("Saudubray et al. (2016) Inborn Metabolic Diseases, 6th ed.")
        references.append("ACMG ACT Sheets for newborn screening conditions")
        references.append("Vockley et al. (2019) ACMG Guideline: Phenylalanine hydroxylase deficiency. Genet Med")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 9 — Gene Therapy Eligibility
# ═══════════════════════════════════════════════════════════════════════════


class GeneTherapyEligibilityWorkflow(BaseRareDiseaseWorkflow):
    """Gene therapy eligibility assessment.

    Matches patient diagnosis + genotype against 6 FDA-approved gene
    therapies and investigational programs, applying age/weight/genotype-
    specific eligibility checks.

    Inputs
    ------
    diagnosis : str
        Confirmed disease diagnosis.
    gene : str, optional
        Causative gene.
    genotype : str, optional
        Specific genotype/mutation (e.g., 'homozygous SMN1 deletion').
    age : str, optional
        Patient age.
    weight_kg : float, optional
        Patient weight in kg.
    prior_therapies : list[str], optional
        Prior treatments received.
    aav_antibodies : bool, optional
        Whether pre-existing AAV antibodies have been tested.
    liver_function : str, optional
        Liver function status.
    """

    workflow_type = DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY

    # Approved gene therapies with eligibility criteria
    _APPROVED_THERAPIES: List[Dict] = [
        {"name": "Zolgensma (onasemnogene abeparvovec)", "gene": "SMN1",
         "disease": "Spinal Muscular Atrophy",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": 24, "weight_max_kg": 13.5,
         "requirements": ["Bi-allelic SMN1 mutation", "Age <2 years", "Weight <13.5 kg",
                           "Anti-AAV9 antibody titer <1:50", "No active hepatitis"],
         "vector": "AAV9", "mechanism": "SMN1 gene replacement"},
        {"name": "Luxturna (voretigene neparvovec)", "gene": "RPE65",
         "disease": "RPE65-mediated Inherited Retinal Dystrophy",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": None, "weight_max_kg": None,
         "requirements": ["Bi-allelic RPE65 mutation", "Viable retinal cells on OCT",
                           "No active ocular infection"],
         "vector": "AAV2", "mechanism": "RPE65 gene replacement"},
        {"name": "Hemgenix (etranacogene dezaparvovec)", "gene": "F9",
         "disease": "Hemophilia B",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": None, "weight_max_kg": None,
         "requirements": ["F9 mutation confirmed", "Age >=18 years", "No F9 inhibitors",
                           "Anti-AAV5 antibody negative"],
         "vector": "AAV5", "mechanism": "F9 gene replacement"},
        {"name": "Roctavian (valoctocogene roxaparvovec)", "gene": "F8",
         "disease": "Hemophilia A",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": None, "weight_max_kg": None,
         "requirements": ["F8 mutation confirmed", "Severe hemophilia A (FVIII <1%)", "Age >=18 years",
                           "No F8 inhibitors", "Anti-AAV5 antibody negative"],
         "vector": "AAV5", "mechanism": "F8 gene replacement"},
        {"name": "Skysona (elivaldogene autotemcel)", "gene": "ABCD1",
         "disease": "Cerebral Adrenoleukodystrophy",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": 204, "weight_max_kg": None,
         "requirements": ["ABCD1 mutation confirmed", "Age 4-17 years", "Early cerebral disease (Loes <= 9)",
                           "Adequate performance IQ", "HLA-matched donor not available"],
         "vector": "Lentiviral", "mechanism": "ABCD1 gene addition in autologous HSCs"},
        {"name": "Casgevy (exagamglogene autotemcel)", "gene": "BCL11A",
         "disease": "Sickle Cell Disease / Transfusion-Dependent Beta-Thalassemia",
         "status": TherapyStatus.APPROVED_FDA,
         "age_max_months": None, "weight_max_kg": None,
         "requirements": ["Confirmed SCD or TDT diagnosis", "Age >=12 years",
                           "Adequate organ function for myeloablative conditioning"],
         "vector": "CRISPR-edited autologous HSCs", "mechanism": "BCL11A enhancer disruption"},
    ]

    _INVESTIGATIONAL: List[Dict] = [
        {"name": "Fidanacogene elaparvovec", "gene": "F9", "disease": "Hemophilia B",
         "status": TherapyStatus.INVESTIGATIONAL, "trial_id": "NCT03861273"},
        {"name": "Giroctocogene fitelparvovec", "gene": "F8", "disease": "Hemophilia A",
         "status": TherapyStatus.INVESTIGATIONAL, "trial_id": "NCT03370913"},
        {"name": "AVXS-101 IT", "gene": "SMN1", "disease": "SMA Type 2",
         "status": TherapyStatus.INVESTIGATIONAL, "trial_id": "NCT03381729"},
        {"name": "AT132 (resamirigene bilparvovec)", "gene": "MTM1", "disease": "X-linked Myotubular Myopathy",
         "status": TherapyStatus.INVESTIGATIONAL, "trial_id": "NCT03199469"},
        {"name": "BMN 307", "gene": "PAH", "disease": "Phenylketonuria",
         "status": TherapyStatus.INVESTIGATIONAL, "trial_id": "NCT04480567"},
    ]

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("diagnosis"):
            warnings.append("No diagnosis provided — eligibility assessment will be limited")
            inp["diagnosis"] = ""
        inp.setdefault("gene", None)
        inp.setdefault("genotype", None)
        inp.setdefault("age", None)
        inp.setdefault("weight_kg", None)
        inp.setdefault("prior_therapies", [])
        inp.setdefault("aav_antibodies", None)
        inp.setdefault("liver_function", None)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        therapies: List[TherapyMatch] = []

        diagnosis = inputs.get("diagnosis", "").lower()
        gene = inputs.get("gene", "") or ""
        genotype = inputs.get("genotype", "") or ""
        age_str = inputs.get("age", "") or ""
        weight_kg = inputs.get("weight_kg")

        findings.append(f"Gene therapy eligibility evaluation for: {diagnosis}")
        if gene:
            findings.append(f"Causative gene: {gene}")
        if genotype:
            findings.append(f"Genotype: {genotype}")

        # Check approved therapies
        for therapy in self._APPROVED_THERAPIES:
            # Match by disease name or gene
            disease_match = therapy["disease"].lower() in diagnosis or diagnosis in therapy["disease"].lower()
            gene_match = gene.upper() == therapy["gene"].upper() if gene else False

            if not disease_match and not gene_match:
                continue

            findings.append(f"Potential match: {therapy['name']} ({therapy['status'].value})")

            # Eligibility checks
            eligible = True
            issues: List[str] = []

            # Age check
            if therapy.get("age_max_months") and age_str:
                age_months = self._parse_age_months(age_str)
                if age_months is not None and age_months > therapy["age_max_months"]:
                    eligible = False
                    issues.append(f"Age exceeds maximum ({therapy['age_max_months']} months)")

            # Weight check
            if therapy.get("weight_max_kg") and weight_kg is not None:
                if weight_kg > therapy["weight_max_kg"]:
                    eligible = False
                    issues.append(f"Weight exceeds maximum ({therapy['weight_max_kg']} kg)")

            if eligible:
                findings.append(f"  ELIGIBLE: {therapy['name']} — meets initial criteria")
                findings.append(f"  Requirements: {'; '.join(therapy['requirements'])}")
                recommendations.append(f"Refer for {therapy['name']} evaluation")
            else:
                findings.append(f"  NOT ELIGIBLE: {'; '.join(issues)}")

            therapies.append(TherapyMatch(
                therapy_name=therapy["name"],
                indication=therapy["disease"],
                status=therapy["status"],
                eligibility_criteria="; ".join(therapy["requirements"]),
                mechanism=therapy["mechanism"],
                gene_target=therapy["gene"],
                orphan_designation=True,
            ))

        # Check investigational therapies
        for trial in self._INVESTIGATIONAL:
            disease_match = trial["disease"].lower() in diagnosis or diagnosis in trial["disease"].lower()
            gene_match = gene.upper() == trial["gene"].upper() if gene else False

            if disease_match or gene_match:
                findings.append(f"Investigational: {trial['name']} (NCT: {trial['trial_id']})")
                recommendations.append(f"Consider clinical trial: {trial['name']} ({trial['trial_id']})")
                therapies.append(TherapyMatch(
                    therapy_name=trial["name"],
                    indication=trial["disease"],
                    status=TherapyStatus.INVESTIGATIONAL,
                    trial_id=trial["trial_id"],
                    gene_target=trial["gene"],
                    orphan_designation=True,
                ))

        if not therapies:
            findings.append("No approved or investigational gene therapies matched")
            recommendations.append("Monitor ClinicalTrials.gov for emerging gene therapy trials")
            recommendations.append("Consider referral to NIH Undiagnosed Disease Network for novel approaches")

        confidence = _clamp(0.3 + 0.15 * len(therapies))
        severity = SeverityLevel.MODERATE if therapies else SeverityLevel.LOW

        references.append("FDA Approved Cellular and Gene Therapy Products (2024)")
        references.append("Mendell et al. (2021) Gene therapy for neuromuscular disorders. N Engl J Med")
        references.append("High & Roncarolo (2019) Gene therapy. N Engl J Med")

        diagnostic = DiagnosticResult(
            therapies=therapies,
            confidence=confidence,
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            confidence=confidence,
            diagnostic_result=diagnostic,
        )

    @staticmethod
    def _parse_age_months(age_str: str) -> Optional[int]:
        """Parse age string to months (approximate)."""
        age_lower = age_str.lower().strip()
        try:
            if "month" in age_lower:
                return int("".join(c for c in age_lower.split("month")[0] if c.isdigit()))
            elif "year" in age_lower:
                years = int("".join(c for c in age_lower.split("year")[0] if c.isdigit()))
                return years * 12
            elif "day" in age_lower or "week" in age_lower:
                return 0  # Neonatal
            else:
                # Try plain number as years
                num = int("".join(c for c in age_lower if c.isdigit()))
                return num * 12 if num < 100 else num
        except (ValueError, IndexError):
            return None


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 10 — Undiagnosed Disease Program
# ═══════════════════════════════════════════════════════════════════════════


class UndiagnosedDiseaseWorkflow(BaseRareDiseaseWorkflow):
    """Undiagnosed Disease Program (UDP) workflow combining multi-modal data.

    Integrates phenotype + genomic + metabolic + imaging data following the
    NIH UDP model for patients without a diagnosis after standard workup.
    Provides reanalysis recommendations and research referral options.

    Inputs
    ------
    hpo_terms : list[str]
        Patient HPO terms.
    genomic_data : dict, optional
        Prior genomic testing results (wes_done, wgs_done, variants_of_interest).
    metabolic_data : dict, optional
        Metabolic screening results.
    imaging_data : dict, optional
        Imaging findings (MRI, X-ray, etc.).
    prior_evaluations : list[str], optional
        List of prior specialist evaluations.
    years_undiagnosed : int, optional
        Years the patient has been without diagnosis.
    """

    workflow_type = DiagnosticWorkflowType.UNDIAGNOSED_DISEASE

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("hpo_terms") and not inp.get("genomic_data") and not inp.get("clinical_summary"):
            warnings.append("Minimal data provided — UDP evaluation will be limited")
        inp.setdefault("hpo_terms", [])
        inp.setdefault("genomic_data", {})
        inp.setdefault("metabolic_data", {})
        inp.setdefault("imaging_data", {})
        inp.setdefault("prior_evaluations", [])
        inp.setdefault("years_undiagnosed", None)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []
        triggers: List[str] = []

        hpo_terms = inputs.get("hpo_terms", [])
        genomic = inputs.get("genomic_data", {})
        metabolic = inputs.get("metabolic_data", {})
        imaging = inputs.get("imaging_data", {})
        prior_evals = inputs.get("prior_evaluations", [])
        years = inputs.get("years_undiagnosed")

        findings.append("Undiagnosed Disease Program (UDP) evaluation initiated")
        if years:
            findings.append(f"Patient has been undiagnosed for approximately {years} years")

        # Data completeness assessment
        completeness_score = 0
        data_summary: List[str] = []

        if hpo_terms:
            completeness_score += 1
            data_summary.append(f"Phenotype: {len(hpo_terms)} HPO terms")
        else:
            recommendations.append("Perform comprehensive phenotyping with HPO terms")

        if genomic:
            completeness_score += 1
            wes = genomic.get("wes_done", False)
            wgs = genomic.get("wgs_done", False)
            data_summary.append(f"Genomic: WES={'done' if wes else 'not done'}, WGS={'done' if wgs else 'not done'}")

            # Reanalysis recommendations
            if wes and not wgs:
                recommendations.append("Consider WGS — detects structural variants, deep intronic, and non-coding variants missed by WES")
            if wes:
                recommendations.append("Re-analyse WES data against updated gene-disease databases (yield: 10-15% on reanalysis)")

            voi = genomic.get("variants_of_interest", [])
            if voi:
                findings.append(f"Variants of interest from prior testing: {len(voi)}")
                for v in voi[:5]:
                    findings.append(f"  - {v}")
                recommendations.append("Review VUS with updated ACMG criteria and ClinVar submissions")
        else:
            recommendations.append("Genomic testing (WES/WGS) is essential for UDP evaluation")

        if metabolic:
            completeness_score += 1
            data_summary.append("Metabolic screening: available")
        else:
            recommendations.append("Consider comprehensive metabolic panel if not already performed")

        if imaging:
            completeness_score += 1
            data_summary.append("Imaging data: available")
            for modality, finding in imaging.items():
                findings.append(f"Imaging ({modality}): {finding}")

        findings.append(f"Data completeness: {completeness_score}/4 domains available")
        if data_summary:
            findings.append("Available data: " + "; ".join(data_summary))

        # Prior evaluations
        if prior_evals:
            findings.append(f"Prior specialist evaluations: {len(prior_evals)}")
            for eval_item in prior_evals:
                findings.append(f"  - {eval_item}")

        # Multi-modal integration recommendations
        if completeness_score >= 3:
            recommendations.append("Multi-modal data integration: combine phenotype-genomic-metabolic data for pattern analysis")
            recommendations.append("Consider computational phenotype matching (Exomiser, LIRICAL, Phen2Gene)")

        # Research program referral
        recommendations.append("Consider referral to NIH Undiagnosed Diseases Program (UDP)")
        recommendations.append("Consider Undiagnosed Diseases Network (UDN) application")
        recommendations.append("Consider SOLVE-RD (European) or equivalent international programs")

        # Additional testing modalities
        if not genomic.get("rna_seq"):
            recommendations.append("Consider RNA sequencing from affected tissue (diagnostic yield: 7-17%)")
        if not genomic.get("long_read"):
            recommendations.append("Consider long-read sequencing (PacBio/ONT) for repeat expansions and structural variants")
        if not metabolic.get("metabolomics"):
            recommendations.append("Consider untargeted metabolomics for novel biomarker discovery")

        # Phenotype-driven search if HPO terms available
        if hpo_terms:
            triggers.append(_trigger_string("PHENOTYPE", "Run phenotype-driven differential with HPO terms"))

        confidence = _clamp(0.15 + 0.05 * completeness_score)
        severity = SeverityLevel.MODERATE

        references.append("Splinter et al. (2018) Effect of genetic diagnosis on patients with previously undiagnosed disease. N Engl J Med")
        references.append("NIH Undiagnosed Diseases Program (https://undiagnosed.hms.harvard.edu/)")
        references.append("Taruscio et al. (2015) Undiagnosed diseases: international efforts. Orphanet J Rare Dis")
        references.append("Boycott et al. (2017) International cooperation to enable diagnosis of all rare genetic diseases. Am J Hum Genet")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=confidence,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW ENGINE
# ═══════════════════════════════════════════════════════════════════════════


class WorkflowEngine:
    """Central dispatcher that maps DiagnosticWorkflowType to the appropriate
    workflow implementation and handles query-based workflow detection."""

    _KEYWORD_MAP: Dict[str, DiagnosticWorkflowType] = {
        # Phenotype-Driven
        "phenotype": DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
        "hpo": DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
        "differential diagnosis": DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
        "phenotype driven": DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
        "phenotype matching": DiagnosticWorkflowType.PHENOTYPE_DRIVEN,
        # WES/WGS Interpretation
        "wes": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "wgs": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "exome": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "genome sequencing": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "variant interpretation": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "acmg": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        "variant classification": DiagnosticWorkflowType.WES_WGS_INTERPRETATION,
        # Metabolic Screening
        "metabolic screening": DiagnosticWorkflowType.METABOLIC_SCREENING,
        "newborn screening": DiagnosticWorkflowType.METABOLIC_SCREENING,
        "metabolic panel": DiagnosticWorkflowType.METABOLIC_SCREENING,
        "tandem mass": DiagnosticWorkflowType.METABOLIC_SCREENING,
        "amino acid": DiagnosticWorkflowType.METABOLIC_SCREENING,
        "acylcarnitine": DiagnosticWorkflowType.METABOLIC_SCREENING,
        # Dysmorphology
        "dysmorphology": DiagnosticWorkflowType.DYSMORPHOLOGY,
        "dysmorphic": DiagnosticWorkflowType.DYSMORPHOLOGY,
        "facial features": DiagnosticWorkflowType.DYSMORPHOLOGY,
        "syndrome": DiagnosticWorkflowType.DYSMORPHOLOGY,
        "syndromic": DiagnosticWorkflowType.DYSMORPHOLOGY,
        # Neurogenetic
        "neurogenetic": DiagnosticWorkflowType.NEUROGENETIC,
        "epilepsy": DiagnosticWorkflowType.NEUROGENETIC,
        "seizure": DiagnosticWorkflowType.NEUROGENETIC,
        "developmental delay": DiagnosticWorkflowType.NEUROGENETIC,
        "intellectual disability": DiagnosticWorkflowType.NEUROGENETIC,
        "movement disorder": DiagnosticWorkflowType.NEUROGENETIC,
        "dystonia": DiagnosticWorkflowType.NEUROGENETIC,
        "ataxia": DiagnosticWorkflowType.NEUROGENETIC,
        # Cardiac Genetics
        "cardiac genetics": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "cardiomyopathy": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "long qt": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "brugada": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "arrhythmia": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "channelopathy": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "sudden cardiac death": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "hcm": DiagnosticWorkflowType.CARDIAC_GENETICS,
        "arvc": DiagnosticWorkflowType.CARDIAC_GENETICS,
        # Connective Tissue
        "connective tissue": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "marfan": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "ehlers danlos": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "eds": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "hypermobility": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "beighton": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        "ghent": DiagnosticWorkflowType.CONNECTIVE_TISSUE,
        # Inborn Errors
        "inborn error": DiagnosticWorkflowType.INBORN_ERRORS,
        "enzyme deficiency": DiagnosticWorkflowType.INBORN_ERRORS,
        "lysosomal storage": DiagnosticWorkflowType.INBORN_ERRORS,
        "metabolic crisis": DiagnosticWorkflowType.INBORN_ERRORS,
        "urea cycle": DiagnosticWorkflowType.INBORN_ERRORS,
        "organic acidemia": DiagnosticWorkflowType.INBORN_ERRORS,
        # Gene Therapy Eligibility
        "gene therapy": DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
        "zolgensma": DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
        "luxturna": DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
        "gene replacement": DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
        "aav": DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY,
        # Undiagnosed Disease
        "undiagnosed": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
        "udp": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
        "udn": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
        "diagnostic odyssey": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
        "no diagnosis": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
        "reanalysis": DiagnosticWorkflowType.UNDIAGNOSED_DISEASE,
    }

    def __init__(self) -> None:
        workflow_instances: List[BaseRareDiseaseWorkflow] = [
            PhenotypeDrivenWorkflow(),
            WESWGSInterpretationWorkflow(),
            MetabolicScreeningWorkflow(),
            DysmorphologyWorkflow(),
            NeurogeneticWorkflow(),
            CardiacGeneticsWorkflow(),
            ConnectiveTissueWorkflow(),
            InbornErrorsWorkflow(),
            GeneTherapyEligibilityWorkflow(),
            UndiagnosedDiseaseWorkflow(),
        ]
        self._workflows: Dict[DiagnosticWorkflowType, BaseRareDiseaseWorkflow] = {
            wf.workflow_type: wf for wf in workflow_instances
        }

    @property
    def _WORKFLOWS(self) -> Dict[DiagnosticWorkflowType, BaseRareDiseaseWorkflow]:
        """Public access to the workflows registry."""
        return self._workflows

    # ── public API ────────────────────────────────────────────────────

    def run_workflow(
        self, workflow_type: DiagnosticWorkflowType, inputs: dict
    ) -> WorkflowResult:
        """Execute a specific workflow by type."""
        wf = self._workflows.get(workflow_type)
        if wf is None:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
        return wf.run(inputs)

    def detect_workflow(self, query: str) -> DiagnosticWorkflowType:
        """Detect the most appropriate workflow from a free-text query."""
        query_lower = query.lower()
        # Score each workflow type by keyword matches
        scores: Dict[DiagnosticWorkflowType, int] = {}
        for keyword, wf_type in self._KEYWORD_MAP.items():
            if keyword in query_lower:
                scores[wf_type] = scores.get(wf_type, 0) + len(keyword)

        if scores:
            return max(scores, key=scores.get)
        return DiagnosticWorkflowType.GENERAL

    def list_workflows(self) -> List[str]:
        """Return list of registered workflow type values."""
        return [wt.value for wt in self._workflows]
