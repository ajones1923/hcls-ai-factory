"""Decision support engines for the Single-Cell Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements four clinical decision support engines that provide
specialised analytical capabilities for single-cell data interpretation:

1. TMEClassifier -- immune phenotype classification with treatment recs
2. SubclonalRiskScorer -- escape risk assessment with timeline
3. TargetExpressionValidator -- CAR-T target safety evaluation
4. CellularDeconvolutionEngine -- bulk-to-single-cell deconvolution
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional

from src.models import (
    EvidenceLevel,
    ResistanceRisk,
    SeverityLevel,
    TMEClass,
)

logger = logging.getLogger(__name__)


# ===================================================================
# HELPERS
# ===================================================================

def _clamp(val: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, val))


def _safe_div(num: float, den: float, default: float = 0.0) -> float:
    return num / den if den != 0 else default


# ===================================================================
# ENGINE 1 -- TME Classifier
# ===================================================================


class TMEClassifier:
    """Classify the tumor microenvironment into one of four immune
    phenotypes and generate treatment recommendations.

    Scoring logic
    -------------
    - Immune infiltration: CD8+ T cells >15% of total => hot signal
    - Stromal density: fibroblast/CAF fraction
    - PD-L1 expression: gene-level or TPS score
    - Exhaustion signature: terminal exhaustion markers

    Output: TMEClass (HOT_INFLAMED, COLD_DESERT, EXCLUDED,
    IMMUNOSUPPRESSIVE) with treatment recommendations.
    """

    # Immune cell types that contribute to infiltration score
    _IMMUNE_TYPES = [
        "CD8_T", "CD4_T", "NK", "B_cell", "Macrophage_M1",
        "Dendritic", "Plasma", "Neutrophil",
    ]
    # Suppressive cell types
    _SUPPRESSIVE_TYPES = ["Treg", "MDSC", "Macrophage_M2"]
    # Stromal cell types
    _STROMAL_TYPES = ["Fibroblast", "CAF", "Pericyte", "Myofibroblast"]
    # Checkpoint genes
    _CHECKPOINT_GENES = [
        "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    ]
    # Immunosuppressive cytokines / enzymes
    _SUPPRESSIVE_GENES = ["IDO1", "TGFB1", "IL10", "VEGFA", "ARG1", "NOS2"]

    def classify(
        self,
        cell_type_proportions: Dict[str, float],
        gene_expression: Optional[Dict[str, float]] = None,
        pdl1_tps: Optional[float] = None,
        spatial_context: Optional[str] = None,
    ) -> Dict:
        """Classify TME and return structured result.

        Parameters
        ----------
        cell_type_proportions : dict
            cell_type -> proportion (0-1).
        gene_expression : dict | None
            gene -> mean expression value.
        pdl1_tps : float | None
            PD-L1 tumor proportion score (0-100).
        spatial_context : str | None
            'infiltrating', 'margin', 'absent'.

        Returns
        -------
        dict with keys:
            tme_class, immune_score, stromal_score, suppressive_score,
            pdl1_status, treatment_recommendations, evidence_level,
            severity, detail.
        """
        gene_expression = gene_expression or {}

        # --- Score computation ---
        immune_score = sum(
            cell_type_proportions.get(ct, 0.0) for ct in self._IMMUNE_TYPES
        )
        cd8_pct = cell_type_proportions.get("CD8_T", 0.0)
        suppressive_frac = sum(
            cell_type_proportions.get(ct, 0.0) for ct in self._SUPPRESSIVE_TYPES
        )
        stromal_frac = sum(
            cell_type_proportions.get(ct, 0.0) for ct in self._STROMAL_TYPES
        )

        # Checkpoint expression
        checkpoint_score = 0.0
        active_checkpoints: List[str] = []
        for gene in self._CHECKPOINT_GENES:
            val = gene_expression.get(gene, 0.0)
            if val > 1.0:
                checkpoint_score += min(val / 5.0, 1.0)
                active_checkpoints.append(gene)
        checkpoint_score = _clamp(
            checkpoint_score / max(len(self._CHECKPOINT_GENES), 1)
        )

        # Suppressive gene score
        suppressive_gene_score = 0.0
        for gene in self._SUPPRESSIVE_GENES:
            val = gene_expression.get(gene, 0.0)
            if val > 1.0:
                suppressive_gene_score += min(val / 5.0, 1.0)
        suppressive_gene_score = _clamp(
            suppressive_gene_score / max(len(self._SUPPRESSIVE_GENES), 1)
        )

        # Combined suppressive score
        total_suppressive = _clamp(
            0.5 * suppressive_frac / 0.2 + 0.5 * suppressive_gene_score
        )

        # PD-L1 status
        if pdl1_tps is not None:
            pdl1_high = pdl1_tps >= 50
            pdl1_positive = pdl1_tps >= 1
        else:
            pdl1_val = gene_expression.get("CD274", 0.0)
            pdl1_high = pdl1_val > 2.0
            pdl1_positive = pdl1_val > 0.5

        # --- Classification ---
        tme_class = self._classify(
            cd8_pct=cd8_pct,
            immune_score=immune_score,
            suppressive_score=total_suppressive,
            stromal_frac=stromal_frac,
            spatial=spatial_context,
            pdl1_high=pdl1_high,
        )

        # --- Treatment recommendations ---
        recommendations = self._get_recommendations(
            tme_class, cd8_pct, suppressive_frac, stromal_frac,
            pdl1_high, pdl1_positive, active_checkpoints,
            cell_type_proportions,
        )

        # Evidence level
        if spatial_context and pdl1_tps is not None:
            evidence = EvidenceLevel.STRONG
        elif pdl1_tps is not None or spatial_context:
            evidence = EvidenceLevel.MODERATE
        else:
            evidence = EvidenceLevel.LIMITED

        # Severity
        severity = SeverityLevel.INFORMATIONAL
        if tme_class == TMEClass.COLD_DESERT:
            severity = SeverityLevel.MODERATE
        elif tme_class == TMEClass.IMMUNOSUPPRESSIVE:
            severity = SeverityLevel.HIGH
        elif tme_class == TMEClass.EXCLUDED:
            severity = SeverityLevel.MODERATE

        return {
            "tme_class": tme_class.value,
            "immune_score": round(immune_score, 4),
            "stromal_score": round(stromal_frac, 4),
            "suppressive_score": round(total_suppressive, 4),
            "cd8_pct": round(cd8_pct, 4),
            "checkpoint_score": round(checkpoint_score, 4),
            "pdl1_status": "high" if pdl1_high else ("positive" if pdl1_positive else "negative"),
            "treatment_recommendations": recommendations,
            "evidence_level": evidence.value,
            "severity": severity.value,
            "detail": {
                "active_checkpoints": active_checkpoints,
                "suppressive_cell_fraction": round(suppressive_frac, 4),
                "spatial_context": spatial_context,
            },
        }

    @staticmethod
    def _classify(
        cd8_pct: float,
        immune_score: float,
        suppressive_score: float,
        stromal_frac: float,
        spatial: Optional[str],
        pdl1_high: bool,
    ) -> TMEClass:
        """Core classification logic."""
        # Spatial overrides
        if spatial == "absent" and immune_score < 0.05:
            return TMEClass.COLD_DESERT
        if spatial == "margin" and immune_score > 0.05:
            return TMEClass.EXCLUDED

        # Score-based
        if cd8_pct >= 0.15 and immune_score >= 0.25:
            if suppressive_score > 0.4:
                return TMEClass.IMMUNOSUPPRESSIVE
            return TMEClass.HOT_INFLAMED
        if immune_score >= 0.10 and stromal_frac > 0.20:
            return TMEClass.EXCLUDED
        if suppressive_score > 0.3 and immune_score >= 0.10:
            return TMEClass.IMMUNOSUPPRESSIVE
        if immune_score < 0.10:
            return TMEClass.COLD_DESERT
        if pdl1_high and cd8_pct >= 0.05:
            return TMEClass.HOT_INFLAMED
        return TMEClass.COLD_DESERT

    @staticmethod
    def _get_recommendations(
        tme_class: TMEClass,
        cd8_pct: float,
        suppressive_frac: float,
        stromal_frac: float,
        pdl1_high: bool,
        pdl1_positive: bool,
        active_checkpoints: List[str],
        proportions: Dict[str, float],
    ) -> List[str]:
        """Generate treatment recommendations based on TME classification."""
        recs: List[str] = []

        if tme_class == TMEClass.HOT_INFLAMED:
            recs.append(
                "Immune-hot TME: strong candidate for checkpoint inhibitor "
                "therapy (anti-PD-1/PD-L1)"
            )
            if pdl1_high:
                recs.append(
                    "PD-L1 TPS >= 50%: consider first-line pembrolizumab monotherapy"
                )
            elif pdl1_positive:
                recs.append(
                    "PD-L1 TPS 1-49%: consider pembrolizumab + chemotherapy combination"
                )
            if "LAG3" in active_checkpoints:
                recs.append(
                    "LAG-3 co-expressed: consider relatlimab + nivolumab (Opdualag)"
                )

        elif tme_class == TMEClass.COLD_DESERT:
            recs.append(
                "Immune-cold TME: low likelihood of single-agent checkpoint "
                "inhibitor benefit"
            )
            recs.append(
                "Consider priming strategies: oncolytic virus (T-VEC), "
                "STING agonist, or radiation to induce immunogenic cell death"
            )
            recs.append(
                "Evaluate bispecific T-cell engager (BiTE) or adoptive "
                "cell therapy to bypass recruitment defect"
            )

        elif tme_class == TMEClass.EXCLUDED:
            recs.append(
                "Immune-excluded TME: target stromal barrier to enable "
                "T-cell infiltration"
            )
            if stromal_frac > 0.3:
                recs.append(
                    f"High stromal density ({stromal_frac:.1%}): consider "
                    "anti-TGFb (bintrafusp alfa) or anti-VEGF combination"
                )
            recs.append(
                "Evaluate anti-CXCL12/CXCR4 axis to promote T-cell migration"
            )

        elif tme_class == TMEClass.IMMUNOSUPPRESSIVE:
            recs.append(
                "Immunosuppressive TME: consider dual checkpoint blockade "
                "(anti-PD-1 + anti-CTLA-4)"
            )
            treg_pct = proportions.get("Treg", 0.0)
            if treg_pct > 0.10:
                recs.append(
                    f"High Treg infiltration ({treg_pct:.1%}): consider anti-CCR8 "
                    "or low-dose cyclophosphamide for selective Treg depletion"
                )
            m2_pct = proportions.get("Macrophage_M2", 0.0)
            if m2_pct > 0.10:
                recs.append(
                    f"M2 macrophage enrichment ({m2_pct:.1%}): consider CSF1R "
                    "inhibitor for macrophage repolarisation"
                )
            if proportions.get("MDSC", 0.0) > 0.05:
                recs.append(
                    "MDSC enrichment: consider entinostat (HDAC inhibitor) "
                    "or all-trans retinoic acid for MDSC differentiation"
                )

        return recs


# ===================================================================
# ENGINE 2 -- Subclonal Risk Scorer
# ===================================================================


class SubclonalRiskScorer:
    """Score subclonal escape risk based on clone frequency,
    proliferation index, and antigen expression.

    Output: HIGH / MEDIUM / LOW risk with estimated timeline.
    """

    def score(
        self,
        clones: List[Dict],
        target_antigen: str = "CD19",
        total_cells: int = 0,
        doubling_time_days: float = 14.0,
    ) -> Dict:
        """Score subclonal escape risk.

        Parameters
        ----------
        clones : list[dict]
            Each clone: {clone_id, cell_count, proliferation_index,
                        antigen_expression, resistance_genes, is_expanding}.
        target_antigen : str
            Therapeutic target antigen gene.
        total_cells : int
            Total cells in sample.
        doubling_time_days : float
            Estimated tumour doubling time in days.

        Returns
        -------
        dict with keys:
            overall_risk, risk_timeline_days, clone_risks, recommendations,
            antigen_negative_fraction, dominant_clone_id.
        """
        if not clones:
            return {
                "overall_risk": ResistanceRisk.LOW.value,
                "risk_timeline_days": None,
                "clone_risks": [],
                "recommendations": ["No subclonal data available"],
                "antigen_negative_fraction": 0.0,
                "dominant_clone_id": None,
            }

        total = max(total_cells, sum(c.get("cell_count", 0) for c in clones))
        clone_risks: List[Dict] = []
        neg_cells = 0
        max_freq = 0.0
        dominant_id = None

        for clone in clones:
            cid = clone.get("clone_id", "unknown")
            cell_count = clone.get("cell_count", 0)
            freq = _safe_div(cell_count, total)
            prolif = clone.get("proliferation_index", 0.0)
            antigen = clone.get("antigen_expression", 1.0)
            res_genes = clone.get("resistance_genes", [])
            expanding = clone.get("is_expanding", False)

            if freq > max_freq:
                max_freq = freq
                dominant_id = cid

            # Antigen-negative check
            is_neg = antigen < 0.1
            if is_neg:
                neg_cells += cell_count

            # Per-clone risk score (0-1)
            risk_score = 0.0
            if is_neg:
                risk_score += 0.4
            if expanding:
                risk_score += 0.2
            risk_score += min(prolif * 0.2, 0.2)
            risk_score += min(len(res_genes) * 0.05, 0.2)
            risk_score = _clamp(risk_score)

            # Risk level
            if risk_score >= 0.6:
                level = ResistanceRisk.HIGH
            elif risk_score >= 0.3:
                level = ResistanceRisk.MEDIUM
            else:
                level = ResistanceRisk.LOW

            clone_risks.append({
                "clone_id": cid,
                "frequency": round(freq, 4),
                "antigen_negative": is_neg,
                "risk_score": round(risk_score, 3),
                "risk_level": level.value,
                "resistance_genes": res_genes[:5],
                "is_expanding": expanding,
            })

        neg_fraction = _safe_div(neg_cells, total)

        # Overall risk
        if neg_fraction > 0.10:
            overall_risk = ResistanceRisk.HIGH
        elif neg_fraction > 0.03:
            overall_risk = ResistanceRisk.MEDIUM
        elif any(cr["risk_level"] == ResistanceRisk.HIGH.value for cr in clone_risks):
            overall_risk = ResistanceRisk.MEDIUM
        else:
            overall_risk = ResistanceRisk.LOW

        # Timeline estimation
        # If antigen-negative clone is expanding, estimate time to dominance
        timeline_days = None
        if neg_fraction > 0.01 and neg_fraction < 0.5:
            # Exponential growth model: time to reach 50%
            # N(t) = N0 * 2^(t/Td)  =>  t = Td * log2(0.5/neg_fraction)
            doublings_needed = math.log2(0.5 / neg_fraction)
            timeline_days = round(doubling_time_days * doublings_needed)
            if timeline_days < 0:
                timeline_days = 0

        # Recommendations
        recommendations: List[str] = []
        if overall_risk == ResistanceRisk.HIGH:
            recommendations.append(
                f"HIGH escape risk: {target_antigen}-negative fraction at "
                f"{neg_fraction:.1%}. Consider dual-target CAR-T or "
                f"bispecific approach immediately"
            )
            if timeline_days is not None:
                recommendations.append(
                    f"Estimated time to {target_antigen}-negative dominance: "
                    f"~{timeline_days} days at current growth rate"
                )
            recommendations.append(
                "Urgent: initiate serial monitoring of antigen-negative "
                "clone expansion at 2-week intervals"
            )
        elif overall_risk == ResistanceRisk.MEDIUM:
            recommendations.append(
                f"MEDIUM escape risk: monitor {target_antigen} expression "
                "at monthly intervals"
            )
            recommendations.append(
                "Prepare contingency plan for antigen loss (alternative "
                "target identification)"
            )
        else:
            recommendations.append(
                f"LOW escape risk: continue standard monitoring for "
                f"{target_antigen} expression"
            )

        # Flag any resistance genes
        all_res = set()
        for cr in clone_risks:
            all_res.update(cr.get("resistance_genes", []))
        if all_res:
            recommendations.append(
                f"Resistance-associated genes detected: {', '.join(sorted(all_res)[:8])}. "
                "Consider pre-emptive combination therapy"
            )

        return {
            "overall_risk": overall_risk.value,
            "risk_timeline_days": timeline_days,
            "clone_risks": clone_risks,
            "recommendations": recommendations,
            "antigen_negative_fraction": round(neg_fraction, 4),
            "dominant_clone_id": dominant_id,
        }


# ===================================================================
# ENGINE 3 -- Target Expression Validator
# ===================================================================


class TargetExpressionValidator:
    """Validate a CAR-T / antibody-drug conjugate target by evaluating
    on-tumour expression, off-tumour vital organ safety, and therapeutic
    index.

    Logic
    -----
    - On-tumour: fraction of tumour cells expressing target > threshold.
    - Off-tumour: check expression in 8 vital organs. Flag if >0.5.
    - Therapeutic index: mean_on_tumour / max_off_tumour.
    """

    # Vital organs for safety check
    _VITAL_ORGANS = [
        "brain", "heart", "lung", "liver", "kidney",
        "pancreas", "bone_marrow", "intestine",
    ]

    # Therapeutic index thresholds
    _TI_FAVORABLE = 10.0
    _TI_ACCEPTABLE = 3.0

    def validate(
        self,
        target_gene: str,
        tumor_expression_values: List[float],
        normal_tissue_expression: Dict[str, float],
        expression_threshold: float = 0.1,
    ) -> Dict:
        """Validate target expression safety and efficacy.

        Parameters
        ----------
        target_gene : str
            Target gene symbol.
        tumor_expression_values : list[float]
            Per-cell expression values in tumour cells.
        normal_tissue_expression : dict[str, float]
            organ/tissue -> mean expression in normal tissue atlas.
        expression_threshold : float
            Minimum expression to count as "expressing" (default 0.1).

        Returns
        -------
        dict with keys:
            target_gene, on_tumor_pct, mean_tumor_expression,
            off_tumor_hits, max_off_tumor, therapeutic_index,
            safety_assessment, efficacy_assessment, overall_verdict,
            recommendations.
        """
        # On-tumour metrics
        if tumor_expression_values:
            expressing = [v for v in tumor_expression_values if v > expression_threshold]
            on_tumor_pct = len(expressing) / len(tumor_expression_values)
            mean_tumor = sum(tumor_expression_values) / len(tumor_expression_values)
        else:
            on_tumor_pct = 0.0
            mean_tumor = 0.0

        # Off-tumour vital organ check
        off_tumor_hits: Dict[str, float] = {}
        max_off_tumor = 0.0
        for organ in self._VITAL_ORGANS:
            organ_expr = normal_tissue_expression.get(organ, 0.0)
            if organ_expr > 0.5:
                off_tumor_hits[organ] = round(organ_expr, 3)
            max_off_tumor = max(max_off_tumor, organ_expr)

        # Also check non-vital tissues
        non_vital_hits: Dict[str, float] = {}
        for tissue, expr in normal_tissue_expression.items():
            if tissue not in self._VITAL_ORGANS and expr > 1.0:
                non_vital_hits[tissue] = round(expr, 3)

        # Therapeutic index
        therapeutic_index = _safe_div(mean_tumor, max_off_tumor + 0.01)

        # Safety assessment
        if not off_tumor_hits:
            safety = "favorable"
            safety_detail = f"No significant {target_gene} expression in vital organs"
        elif any(v > 2.0 for v in off_tumor_hits.values()):
            safety = "high_risk"
            high_organs = [o for o, v in off_tumor_hits.items() if v > 2.0]
            safety_detail = (
                f"HIGH expression in vital organs: {', '.join(high_organs)}. "
                "Requires safety switch (iCasp9) or affinity-tuned design"
            )
        else:
            safety = "moderate_risk"
            safety_detail = (
                f"Low-level expression in: {', '.join(off_tumor_hits.keys())}. "
                "Monitor for on-target off-tumour toxicity"
            )

        # Efficacy assessment
        if on_tumor_pct >= 0.90:
            efficacy = "excellent"
            efficacy_detail = f"{on_tumor_pct:.1%} tumour coverage -- minimal antigen-negative escape risk"
        elif on_tumor_pct >= 0.70:
            efficacy = "adequate"
            efficacy_detail = f"{on_tumor_pct:.1%} tumour coverage -- some antigen-negative cells present"
        elif on_tumor_pct >= 0.50:
            efficacy = "marginal"
            efficacy_detail = f"{on_tumor_pct:.1%} tumour coverage -- significant escape risk"
        else:
            efficacy = "insufficient"
            efficacy_detail = f"{on_tumor_pct:.1%} tumour coverage -- unlikely to achieve durable response"

        # Overall verdict
        if safety == "favorable" and efficacy in ("excellent", "adequate"):
            verdict = "FAVORABLE"
        elif safety == "high_risk":
            verdict = "UNFAVORABLE"
        elif efficacy == "insufficient":
            verdict = "UNFAVORABLE"
        elif safety == "moderate_risk" and efficacy in ("excellent", "adequate"):
            verdict = "CONDITIONAL"
        else:
            verdict = "CONDITIONAL"

        # Therapeutic index assessment
        if therapeutic_index >= self._TI_FAVORABLE:
            ti_assessment = "favorable"
        elif therapeutic_index >= self._TI_ACCEPTABLE:
            ti_assessment = "acceptable"
        else:
            ti_assessment = "unfavorable"

        # Recommendations
        recommendations: List[str] = []
        if verdict == "FAVORABLE":
            recommendations.append(
                f"{target_gene}: strong target candidate. Proceed with "
                "CAR-T construct development"
            )
        elif verdict == "CONDITIONAL":
            recommendations.append(
                f"{target_gene}: viable with risk mitigation. Consider "
                "dose escalation protocol and enhanced monitoring"
            )
            if safety != "favorable":
                recommendations.append(
                    "Incorporate safety switch (iCasp9 or EGFRt) in CAR design"
                )
            if efficacy in ("marginal", "adequate"):
                recommendations.append(
                    "Consider tandem or dual-target CAR to improve coverage"
                )
        else:
            recommendations.append(
                f"{target_gene}: not recommended as primary target. "
                "Explore alternative antigens"
            )
            if safety == "high_risk":
                recommendations.append(
                    f"Vital organ toxicity risk from {', '.join(off_tumor_hits.keys())}"
                )
            if efficacy == "insufficient":
                recommendations.append(
                    f"Only {on_tumor_pct:.1%} on-tumour expression -- "
                    "insufficient for therapeutic efficacy"
                )

        return {
            "target_gene": target_gene,
            "on_tumor_pct": round(on_tumor_pct, 4),
            "mean_tumor_expression": round(mean_tumor, 3),
            "off_tumor_hits": off_tumor_hits,
            "non_vital_hits": non_vital_hits,
            "max_off_tumor": round(max_off_tumor, 3),
            "therapeutic_index": round(therapeutic_index, 2),
            "ti_assessment": ti_assessment,
            "safety_assessment": safety,
            "safety_detail": safety_detail,
            "efficacy_assessment": efficacy,
            "efficacy_detail": efficacy_detail,
            "overall_verdict": verdict,
            "recommendations": recommendations,
        }


# ===================================================================
# ENGINE 4 -- Cellular Deconvolution Engine
# ===================================================================


class CellularDeconvolutionEngine:
    """Estimate cell type proportions from bulk RNA-seq expression data
    using a simplified CIBERSORTx-style approach.

    Uses a reference signature matrix (cell_type -> {gene: expression})
    and solves for proportions via non-negative least squares (simplified).
    """

    # Default reference signature (representative marker genes per cell type)
    _DEFAULT_SIGNATURE: Dict[str, Dict[str, float]] = {
        "CD8_T": {
            "CD8A": 8.5, "CD8B": 7.2, "GZMB": 6.8, "PRF1": 5.5,
            "CD3D": 7.0, "CD3E": 6.5, "NKG7": 4.5, "IFNG": 3.2,
        },
        "CD4_T": {
            "CD4": 7.0, "IL7R": 6.5, "CCR7": 5.8, "CD3D": 7.0,
            "CD3E": 6.5, "LEF1": 4.0, "TCF7": 4.5, "SELL": 3.8,
        },
        "Treg": {
            "FOXP3": 8.0, "IL2RA": 7.5, "CTLA4": 6.0, "IKZF2": 5.5,
            "CD4": 5.0, "TNFRSF18": 4.5, "TIGIT": 4.0, "CD3D": 5.0,
        },
        "NK": {
            "NCAM1": 7.5, "NKG7": 8.0, "GNLY": 7.8, "KLRD1": 6.5,
            "KLRB1": 6.0, "KLRC1": 5.5, "GZMB": 5.0, "PRF1": 4.5,
        },
        "B_cell": {
            "CD19": 8.0, "MS4A1": 7.5, "CD79A": 7.8, "CD79B": 7.0,
            "PAX5": 5.5, "CD22": 5.0, "BANK1": 4.5, "BLK": 4.0,
        },
        "Macrophage": {
            "CD68": 7.5, "CD163": 6.8, "CSF1R": 6.5, "MARCO": 5.5,
            "MRC1": 5.0, "MSR1": 4.5, "CD14": 6.0, "ITGAM": 5.5,
        },
        "Dendritic": {
            "ITGAX": 7.0, "CLEC9A": 6.5, "CD1C": 6.0, "FCER1A": 5.5,
            "HLA-DRA": 7.5, "HLA-DRB1": 7.0, "CD80": 4.5, "CD86": 4.0,
        },
        "Fibroblast": {
            "COL1A1": 8.5, "COL1A2": 8.0, "DCN": 7.5, "LUM": 7.0,
            "PDGFRA": 5.5, "FAP": 5.0, "ACTA2": 4.5, "VIM": 6.0,
        },
        "Endothelial": {
            "PECAM1": 7.5, "VWF": 7.0, "CDH5": 6.5, "FLT1": 6.0,
            "KDR": 5.5, "CLDN5": 5.0, "EMCN": 4.5, "ERG": 4.0,
        },
        "Epithelial": {
            "EPCAM": 8.0, "KRT18": 7.5, "KRT19": 7.0, "CDH1": 6.5,
            "KRT8": 6.0, "MUC1": 5.0, "CLDN4": 4.5, "TJP1": 4.0,
        },
    }

    def deconvolve(
        self,
        bulk_expression: Dict[str, float],
        signature_matrix: Optional[Dict[str, Dict[str, float]]] = None,
        max_iterations: int = 100,
        convergence_threshold: float = 1e-4,
    ) -> Dict:
        """Estimate cell type proportions from bulk expression.

        Parameters
        ----------
        bulk_expression : dict
            gene -> expression value from bulk RNA-seq.
        signature_matrix : dict | None
            cell_type -> {gene: reference_expression}.  None = use default.
        max_iterations : int
            Maximum NNLS iterations.
        convergence_threshold : float
            Convergence criterion for proportion changes.

        Returns
        -------
        dict with keys:
            proportions, confidence, residual, dominant_cell_type,
            quality_metrics, n_genes_used.
        """
        sig = signature_matrix or self._DEFAULT_SIGNATURE

        # Find overlapping genes
        sig_genes: set = set()
        for ct_genes in sig.values():
            sig_genes.update(ct_genes.keys())
        overlap_genes = sorted(sig_genes & set(bulk_expression.keys()))

        if len(overlap_genes) < 5:
            return {
                "proportions": {},
                "confidence": 0.0,
                "residual": 1.0,
                "dominant_cell_type": None,
                "quality_metrics": {
                    "n_genes_used": len(overlap_genes),
                    "warning": "Insufficient gene overlap for reliable deconvolution",
                },
                "n_genes_used": len(overlap_genes),
            }

        cell_types = list(sig.keys())
        n_types = len(cell_types)
        n_genes = len(overlap_genes)

        # Build signature matrix (genes x cell_types) and bulk vector
        S: List[List[float]] = []
        b: List[float] = []
        for gene in overlap_genes:
            row = [sig[ct].get(gene, 0.0) for ct in cell_types]
            S.append(row)
            b.append(bulk_expression[gene])

        # Simplified NNLS: iterative proportional fitting
        # Initialize uniform proportions
        props = [1.0 / n_types] * n_types

        for iteration in range(max_iterations):
            prev_props = list(props)

            # For each cell type, compute correlation with residual
            for t in range(n_types):
                # Predicted without this cell type
                pred_others = [
                    sum(S[g][j] * props[j] for j in range(n_types) if j != t)
                    for g in range(n_genes)
                ]
                # Residual to explain
                residuals = [b[g] - pred_others[g] for g in range(n_genes)]

                # Optimal proportion for this cell type
                sig_col = [S[g][t] for g in range(n_genes)]
                numerator = sum(residuals[g] * sig_col[g] for g in range(n_genes))
                denominator = sum(sig_col[g] ** 2 for g in range(n_genes))
                if denominator > 0:
                    optimal = numerator / denominator
                    props[t] = max(0.0, optimal)
                else:
                    props[t] = 0.0

            # Normalize to sum to 1
            total_prop = sum(props)
            if total_prop > 0:
                props = [p / total_prop for p in props]
            else:
                props = [1.0 / n_types] * n_types

            # Check convergence
            max_change = max(abs(props[i] - prev_props[i]) for i in range(n_types))
            if max_change < convergence_threshold:
                break

        # Compute residual (goodness of fit)
        predicted = [
            sum(S[g][t] * props[t] for t in range(n_types))
            for g in range(n_genes)
        ]
        ss_res = sum((b[g] - predicted[g]) ** 2 for g in range(n_genes))
        b_mean = sum(b) / max(n_genes, 1)
        ss_tot = sum((b[g] - b_mean) ** 2 for g in range(n_genes))
        r_squared = 1.0 - _safe_div(ss_res, ss_tot, default=1.0)
        r_squared = _clamp(r_squared)

        # Build results
        proportions = {
            ct: round(props[i], 4) for i, ct in enumerate(cell_types)
            if props[i] > 0.005  # Filter out noise below 0.5%
        }

        # Re-normalize filtered proportions
        total_filtered = sum(proportions.values())
        if total_filtered > 0:
            proportions = {
                ct: round(p / total_filtered, 4) for ct, p in proportions.items()
            }

        dominant = max(proportions, key=proportions.get) if proportions else None

        # Confidence based on R-squared and gene overlap
        gene_overlap_score = min(len(overlap_genes) / 50.0, 1.0)
        confidence = _clamp(0.5 * r_squared + 0.5 * gene_overlap_score)

        return {
            "proportions": proportions,
            "confidence": round(confidence, 3),
            "residual": round(1.0 - r_squared, 4),
            "r_squared": round(r_squared, 4),
            "dominant_cell_type": dominant,
            "quality_metrics": {
                "n_genes_used": len(overlap_genes),
                "n_cell_types": len(proportions),
                "iterations": min(iteration + 1, max_iterations),
                "convergence_achieved": max_change < convergence_threshold if 'max_change' in dir() else True,
            },
            "n_genes_used": len(overlap_genes),
        }
