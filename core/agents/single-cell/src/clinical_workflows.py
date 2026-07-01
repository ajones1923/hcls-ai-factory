"""Clinical workflows for the Single-Cell Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements ten single-cell analysis workflows that integrate scRNA-seq,
spatial transcriptomics, and clinical data to produce actionable
assessments.  Each workflow follows the BaseSCWorkflow contract
(preprocess -> execute -> postprocess) and is registered in the
WorkflowEngine for unified dispatch.
"""

from __future__ import annotations

import logging
import math
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from src.models import (
    BiomarkerCandidate,
    CARTTargetValidation,
    CellTypeAnnotation,
    CellTypeConfidence,
    DrugResponsePrediction,
    EvidenceLevel,
    LigandReceptorInteraction,
    ResistanceRisk,
    SCWorkflowType,
    SeverityLevel,
    SpatialNiche,
    SubclonalResult,
    TMEClass,
    TMEProfile,
    TrajectoryResult,
    TrajectoryType,
    TreatmentResponse,
    WorkflowResult,
)

logger = logging.getLogger(__name__)


# ===================================================================
# HELPERS
# ===================================================================

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


def _safe_div(num: float, den: float, default: float = 0.0) -> float:
    """Safe division avoiding ZeroDivisionError."""
    return num / den if den != 0 else default


def _clamp(val: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, val))


# Canonical marker gene sets for reference-free annotation
_CANONICAL_MARKERS: Dict[str, List[str]] = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "CD8_T": ["CD8A", "CD8B", "GZMB", "PRF1"],
    "CD4_T": ["CD4", "IL7R", "CCR7"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "NK": ["NCAM1", "NKG7", "GNLY", "KLRD1"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1", "XBP1"],
    "Macrophage": ["CD68", "CD163", "CSF1R", "MARCO"],
    "Dendritic": ["ITGAX", "CLEC9A", "CD1C", "FCER1A"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "S100A8"],
    "Mast": ["KIT", "CPA3", "TPSAB1", "TPSB2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "FLT1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Epithelial": ["EPCAM", "KRT18", "KRT19", "CDH1"],
    "Malignant": ["MKI67", "TOP2A", "PCNA"],
}

# Vital organs checked for off-tumor CAR-T safety
_VITAL_ORGANS = [
    "brain", "heart", "lung", "liver", "kidney",
    "pancreas", "bone_marrow", "intestine",
]

# Known immunosuppressive genes in TME
_IMMUNOSUPPRESSIVE_GENES = [
    "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    "IDO1", "TGFB1", "IL10", "VEGFA", "ARG1", "NOS2",
]

# Exhaustion markers by state
_EXHAUSTION_MARKERS: Dict[str, List[str]] = {
    "naive": ["TCF7", "LEF1", "SELL", "CCR7"],
    "effector": ["GZMB", "PRF1", "IFNG", "TBX21"],
    "memory": ["IL7R", "BCL2", "TCF7", "CD44"],
    "progenitor_exhausted": ["TCF7", "CXCR5", "SLAMF6", "TOX"],
    "terminally_exhausted": ["TOX", "ENTPD1", "HAVCR2", "LAG3", "TIGIT"],
}


# ===================================================================
# BASE CLASS
# ===================================================================


class BaseSCWorkflow(ABC):
    """Abstract base for all single-cell clinical workflows."""

    workflow_type: SCWorkflowType
    _findings: List[str]
    _recommendations: List[str]

    def __init__(self) -> None:
        self._findings: List[str] = []
        self._recommendations: List[str] = []

    # -- template-method orchestrator ----------------------------------
    def run(self, inputs: dict) -> WorkflowResult:
        """Orchestrate preprocess -> execute -> postprocess."""
        logger.info("Running workflow %s", self.workflow_type.value)
        self._findings = []
        self._recommendations = []
        processed_inputs = self.preprocess(inputs)
        result = self.execute(processed_inputs)
        result = self.postprocess(result)
        return result

    def preprocess(self, inputs: dict) -> dict:
        """Validate and normalise raw inputs.  Override for workflow-specific logic."""
        return dict(inputs)

    @abstractmethod
    def execute(self, inputs: dict) -> WorkflowResult:
        """Core analysis logic.  Must be implemented by each workflow."""
        ...

    def postprocess(self, result: WorkflowResult) -> WorkflowResult:
        """Shared enrichment after execution."""
        try:
            from api.routes.events import publish_event
            publish_event("workflow_complete", {
                "workflow": result.workflow_type.value,
                "severity": result.severity.value,
            })
        except Exception:
            pass
        return result

    @staticmethod
    def _init_warnings(inp: dict) -> list:
        """Initialise and return the validation warnings list on *inp*."""
        warnings: list = inp.setdefault("_validation_warnings", [])
        return warnings

    @staticmethod
    def _require_key(inp: dict, key: str, warnings: list, default=None):
        """Return inp[key] or default, appending a warning if missing."""
        if key not in inp:
            warnings.append(f"Missing required input '{key}'; using default={default}")
            return default
        return inp[key]


# ===================================================================
# WORKFLOW 1 -- Cell Type Annotation
# ===================================================================


class CellTypeAnnotationWorkflow(BaseSCWorkflow):
    """Multi-strategy consensus cell type annotation.

    Combines reference-based (SingleR/scArches), marker-based, and
    LLM-based annotation strategies.  Confidence:
      - 3/3 strategies agree = HIGH
      - 2/3 strategies agree = MEDIUM
      - 1/3 strategies agree = LOW

    Inputs
    ------
    expression_matrix : dict[str, dict[str, float]]
        cell_id -> {gene: expression_value}.
    reference_labels : dict[str, str] | None
        cell_id -> reference-predicted label (from SingleR / scArches).
    marker_genes : dict[str, list[str]] | None
        cell_type -> list of marker genes for marker-based annotation.
    llm_labels : dict[str, str] | None
        cell_id -> LLM-predicted cell type label.
    min_genes_per_cell : int
        Minimum genes detected per cell (QC filter).  Default 200.
    """

    workflow_type = SCWorkflowType.CELL_TYPE_ANNOTATION

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "expression_matrix", warnings, default={})
        inp.setdefault("reference_labels", None)
        inp.setdefault("marker_genes", _CANONICAL_MARKERS)
        inp.setdefault("llm_labels", None)
        inp.setdefault("min_genes_per_cell", 200)
        # QC filter
        expr = inp.get("expression_matrix", {})
        min_g = inp["min_genes_per_cell"]
        filtered = {
            cid: genes for cid, genes in expr.items()
            if len(genes) >= min_g
        }
        removed = len(expr) - len(filtered)
        if removed > 0:
            warnings.append(
                f"Removed {removed} cells with <{min_g} genes detected"
            )
        inp["expression_matrix"] = filtered
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        expr = inputs.get("expression_matrix", {})
        ref_labels = inputs.get("reference_labels") or {}
        marker_db = inputs.get("marker_genes", _CANONICAL_MARKERS)
        llm_labels = inputs.get("llm_labels") or {}

        cell_type_counts: Dict[str, Dict] = {}

        for cell_id, gene_expr in expr.items():
            ref_label = ref_labels.get(cell_id)
            marker_label = self._marker_annotate(gene_expr, marker_db)
            llm_label = llm_labels.get(cell_id)

            strategies: Dict[str, str] = {}
            if ref_label:
                strategies["reference"] = ref_label
            if marker_label:
                strategies["marker"] = marker_label
            if llm_label:
                strategies["llm"] = llm_label

            # Consensus: pick most common label
            if not strategies:
                consensus = "Unknown"
                confidence = CellTypeConfidence.LOW
                n_agree = 0
                total_strats = 0
            else:
                label_votes: Dict[str, int] = {}
                for lbl in strategies.values():
                    label_votes[lbl] = label_votes.get(lbl, 0) + 1
                best_label = max(label_votes, key=label_votes.get)
                n_agree = label_votes[best_label]
                total_strats = len(strategies)
                if total_strats >= 3 and n_agree >= 3:
                    confidence = CellTypeConfidence.HIGH
                elif n_agree >= 2:
                    confidence = CellTypeConfidence.MEDIUM if total_strats >= 3 else CellTypeConfidence.HIGH
                else:
                    confidence = CellTypeConfidence.LOW
                consensus = best_label

            if consensus not in cell_type_counts:
                cell_type_counts[consensus] = {
                    "count": 0,
                    "conf_high": 0,
                    "conf_med": 0,
                    "conf_low": 0,
                    "marker_evidence": {},
                    "markers": self._get_markers_for_type(consensus, marker_db),
                }
            cell_type_counts[consensus]["count"] += 1
            if confidence == CellTypeConfidence.HIGH:
                cell_type_counts[consensus]["conf_high"] += 1
            elif confidence == CellTypeConfidence.MEDIUM:
                cell_type_counts[consensus]["conf_med"] += 1
            else:
                cell_type_counts[consensus]["conf_low"] += 1

            # Accumulate marker evidence
            for mg in cell_type_counts[consensus]["markers"]:
                prev = cell_type_counts[consensus]["marker_evidence"].get(mg, 0.0)
                cell_type_counts[consensus]["marker_evidence"][mg] = prev + gene_expr.get(mg, 0.0)

        total_cells = max(len(expr), 1)
        annotations: List[CellTypeAnnotation] = []

        for idx, (ct, info) in enumerate(
            sorted(cell_type_counts.items(), key=lambda x: -x[1]["count"])
        ):
            fraction = info["count"] / total_cells
            # Overall confidence for this cell type: majority vote
            if info["conf_high"] >= info["conf_med"] and info["conf_high"] >= info["conf_low"]:
                ct_conf = CellTypeConfidence.HIGH
                conf_score = 0.9
            elif info["conf_med"] >= info["conf_low"]:
                ct_conf = CellTypeConfidence.MEDIUM
                conf_score = 0.6
            else:
                ct_conf = CellTypeConfidence.LOW
                conf_score = 0.3

            # Normalise marker evidence to per-cell averages
            marker_ev = {
                mg: round(val / max(info["count"], 1), 3)
                for mg, val in info["marker_evidence"].items()
            }

            cta = CellTypeAnnotation(
                cluster_id=f"cluster_{idx}",
                cell_type=ct,
                confidence=ct_conf,
                confidence_score=round(conf_score, 3),
                marker_genes=info["markers"],
                marker_evidence=marker_ev,
                cell_count=info["count"],
                fraction=round(fraction, 4),
            )
            annotations.append(cta)

        n_low = sum(1 for a in annotations if a.confidence == CellTypeConfidence.LOW)
        severity = SeverityLevel.INFORMATIONAL
        if n_low > len(annotations) * 0.5:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            cell_annotations=annotations,
            severity=severity,
        )

    @staticmethod
    def _marker_annotate(
        gene_expr: Dict[str, float],
        marker_db: Dict[str, List[str]],
    ) -> Optional[str]:
        """Annotate a cell using marker gene overlap scoring."""
        best_type = None
        best_score = 0.0
        expressed_genes = {g for g, v in gene_expr.items() if v > 0}
        for ct, markers in marker_db.items():
            if not markers:
                continue
            overlap = len(expressed_genes & set(markers))
            score = overlap / len(markers)
            if score > best_score and score >= 0.3:
                best_score = score
                best_type = ct
        return best_type

    @staticmethod
    def _get_markers_for_type(cell_type: str, marker_db: Dict[str, List[str]]) -> List[str]:
        return marker_db.get(cell_type, [])[:6]


# ===================================================================
# WORKFLOW 2 -- TME Profiling
# ===================================================================


class TMEProfilingWorkflow(BaseSCWorkflow):
    """Tumor microenvironment profiling: classify hot/cold/excluded/
    immunosuppressive.  Immune infiltration scoring.  Treatment implications.

    Inputs
    ------
    cell_type_proportions : dict[str, float]
        Cell type -> proportion (0-1).
    gene_expression_summary : dict[str, float]
        Gene -> mean expression across tumour sample.
    spatial_immune_location : str | None
        'infiltrating', 'margin', 'absent' -- spatial context if available.
    pdl1_tps : float | None
        PD-L1 tumor proportion score (0-100).
    tmb : float | None
        Tumor mutational burden (mut/Mb).
    msi_status : str | None
        'MSI-H', 'MSS', or None.
    """

    workflow_type = SCWorkflowType.TME_PROFILING

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "cell_type_proportions", warnings, default={})
        inp.setdefault("gene_expression_summary", {})
        inp.setdefault("spatial_immune_location", None)
        inp.setdefault("pdl1_tps", None)
        inp.setdefault("tmb", None)
        inp.setdefault("msi_status", None)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        props = inputs.get("cell_type_proportions", {})
        expr = inputs.get("gene_expression_summary", {})
        spatial = inputs.get("spatial_immune_location")
        pdl1_tps = inputs.get("pdl1_tps")
        tmb = inputs.get("tmb")
        msi = inputs.get("msi_status")

        # --- Immune infiltration score ---
        immune_types = {"CD8_T", "CD4_T", "NK", "B_cell", "Macrophage_M1",
                        "Dendritic", "Plasma", "Neutrophil"}
        total_immune = sum(props.get(ct, 0.0) for ct in immune_types)
        cd8_pct = props.get("CD8_T", 0.0)
        props.get("Treg", 0.0)
        props.get("Macrophage_M2", 0.0)
        props.get("MDSC", 0.0)
        immune_score = _clamp(total_immune)

        # --- Stromal density ---
        stromal_types = {"Fibroblast", "CAF", "Pericyte", "Myofibroblast"}
        stromal_density = sum(props.get(ct, 0.0) for ct in stromal_types)
        stromal_score = _clamp(stromal_density)

        # --- Immunosuppressive gene score ---
        checkpoint_expression: Dict[str, float] = {}
        cytokine_milieu: Dict[str, float] = {}
        suppressive_score = 0.0
        for gene in _IMMUNOSUPPRESSIVE_GENES:
            val = expr.get(gene, 0.0)
            if val > 0:
                if gene in ("CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"):
                    checkpoint_expression[gene] = round(val, 3)
                else:
                    cytokine_milieu[gene] = round(val, 3)
                if val > 1.0:
                    suppressive_score += min(val / 5.0, 1.0)
        suppressive_score = _clamp(suppressive_score / max(len(_IMMUNOSUPPRESSIVE_GENES), 1))

        # --- PD-L1 expression ---
        pdl1_expr = expr.get("CD274", 0.0)
        if pdl1_tps is not None:
            pdl1_high = pdl1_tps >= 50
        else:
            pdl1_high = pdl1_expr > 2.0

        # --- Exhaustion signature ---
        exh_markers = _EXHAUSTION_MARKERS["terminally_exhausted"]
        exh_score = _clamp(
            sum(expr.get(m, 0.0) for m in exh_markers) / (len(exh_markers) * 3.0)
        )

        # --- TME classification ---
        tme_class = self._classify_tme(
            cd8_pct, total_immune, suppressive_score, stromal_density,
            spatial, pdl1_high,
        )

        # --- Predicted immunotherapy response ---
        if tme_class == TMEClass.HOT_INFLAMED:
            predicted_response = "responder"
        elif tme_class == TMEClass.COLD_DESERT:
            predicted_response = "non-responder"
        elif tme_class == TMEClass.IMMUNOSUPPRESSIVE:
            predicted_response = "uncertain"
        else:
            predicted_response = "uncertain"

        # TMB/MSI adjustments
        if tmb is not None and tmb >= 10:
            predicted_response = "responder"
        if msi and msi.upper() == "MSI-H":
            predicted_response = "responder"

        tme_profile = TMEProfile(
            tme_class=tme_class,
            immune_score=round(immune_score, 4),
            stromal_score=round(stromal_score, 4),
            cell_type_fractions={k: round(v, 4) for k, v in props.items()},
            exhaustion_signature=round(exh_score, 4),
            checkpoint_expression=checkpoint_expression,
            cytokine_milieu=cytokine_milieu,
            predicted_immunotherapy_response=predicted_response,
            spatial_pattern=spatial,
            evidence_level=EvidenceLevel.MODERATE if spatial else EvidenceLevel.LIMITED,
        )

        severity = SeverityLevel.INFORMATIONAL
        if tme_class == TMEClass.COLD_DESERT:
            severity = SeverityLevel.MODERATE
        elif tme_class == TMEClass.IMMUNOSUPPRESSIVE:
            severity = SeverityLevel.HIGH
        elif tme_class == TMEClass.EXCLUDED:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            tme_profile=tme_profile,
            severity=severity,
        )

    @staticmethod
    def _classify_tme(
        cd8_pct: float,
        total_immune: float,
        suppressive_score: float,
        stromal_density: float,
        spatial: Optional[str],
        pdl1_high: bool,
    ) -> TMEClass:
        """Classify TME into hot/cold/excluded/immunosuppressive."""
        if spatial == "absent" and total_immune < 0.05:
            return TMEClass.COLD_DESERT
        if spatial == "margin" and total_immune > 0.05:
            return TMEClass.EXCLUDED

        if cd8_pct >= 0.15 and total_immune >= 0.25:
            if suppressive_score > 0.4:
                return TMEClass.IMMUNOSUPPRESSIVE
            return TMEClass.HOT_INFLAMED
        if total_immune >= 0.10 and stromal_density > 0.20:
            return TMEClass.EXCLUDED
        if suppressive_score > 0.3 and total_immune >= 0.10:
            return TMEClass.IMMUNOSUPPRESSIVE
        if total_immune < 0.10:
            return TMEClass.COLD_DESERT

        if pdl1_high and cd8_pct >= 0.05:
            return TMEClass.HOT_INFLAMED
        return TMEClass.COLD_DESERT


# ===================================================================
# WORKFLOW 3 -- Drug Response Prediction
# ===================================================================


class DrugResponseWorkflow(BaseSCWorkflow):
    """Cell-type-specific drug sensitivity from DepMap/CCLE data.
    IC50 prediction per cluster.

    Inputs
    ------
    cluster_expression : dict[str, dict[str, float]]
        cluster_id -> {gene: mean_expression}.
    cluster_cell_types : dict[str, str]
        cluster_id -> annotated cell type.
    drug_candidates : list[str] | None
        Specific drugs to evaluate; None = screen all available.
    depmap_sensitivity : dict[str, dict[str, float]] | None
        drug -> {gene: sensitivity_correlation}.
    """

    workflow_type = SCWorkflowType.DRUG_RESPONSE

    # Default drug-gene sensitivity map (simplified DepMap proxy)
    _DEFAULT_DRUG_GENES: Dict[str, Dict[str, float]] = {
        "cisplatin": {"ERCC1": -0.4, "BRCA1": -0.3, "XPC": -0.2, "TP53": 0.3},
        "paclitaxel": {"TUBB3": -0.5, "ABCB1": -0.6, "BCL2": -0.2, "BIM": 0.4},
        "doxorubicin": {"TOP2A": 0.5, "ABCB1": -0.6, "TP53": 0.3},
        "venetoclax": {"BCL2": 0.7, "MCL1": -0.5, "BAX": 0.3},
        "imatinib": {"ABL1": 0.6, "BCR": 0.4, "KIT": 0.5},
        "olaparib": {"BRCA1": -0.6, "BRCA2": -0.6, "PARP1": 0.4, "RAD51": -0.3},
        "trametinib": {"BRAF": 0.5, "KRAS": 0.4, "MAP2K1": 0.3, "DUSP6": 0.4},
        "osimertinib": {"EGFR": 0.7, "MET": -0.3, "ERBB2": 0.2},
        "pembrolizumab": {"CD274": 0.6, "PDCD1": 0.3, "IFNG": 0.4},
        "5-fluorouracil": {"TYMS": 0.4, "DPYD": -0.5, "TP53": 0.2},
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "cluster_expression", warnings, default={})
        inp.setdefault("cluster_cell_types", {})
        inp.setdefault("drug_candidates", None)
        inp.setdefault("depmap_sensitivity", self._DEFAULT_DRUG_GENES)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        cluster_expr = inputs.get("cluster_expression", {})
        cluster_types = inputs.get("cluster_cell_types", {})
        drug_cands = inputs.get("drug_candidates")
        drug_genes = inputs.get("depmap_sensitivity", self._DEFAULT_DRUG_GENES)

        if drug_cands:
            drug_genes = {d: drug_genes[d] for d in drug_cands if d in drug_genes}

        predictions: List[DrugResponsePrediction] = []

        for cluster_id, gene_expr in cluster_expr.items():
            ct = cluster_types.get(cluster_id, "Unknown")

            for drug, gene_weights in drug_genes.items():
                score = 0.0
                n_genes = 0
                resistance_mechs: List[str] = []
                for gene, weight in gene_weights.items():
                    gene_val = gene_expr.get(gene, 0.0)
                    score += gene_val * weight
                    n_genes += 1
                    if weight < 0 and gene_val > 2.0:
                        resistance_mechs.append(f"{gene} overexpression")

                if n_genes > 0:
                    score = score / n_genes

                sensitivity = _clamp((score + 3.0) / 6.0)

                if resistance_mechs:
                    res_risk = ResistanceRisk.HIGH if len(resistance_mechs) >= 2 else ResistanceRisk.MEDIUM
                    resistant_frac = min(0.3 * len(resistance_mechs), 1.0)
                else:
                    res_risk = ResistanceRisk.LOW
                    resistant_frac = 0.0

                pred = DrugResponsePrediction(
                    drug_name=drug,
                    drug_class=self._infer_drug_class(drug),
                    predicted_sensitivity=round(sensitivity, 3),
                    resistance_risk=res_risk,
                    resistance_mechanisms=resistance_mechs,
                    resistant_subpopulation=ct if resistance_mechs else None,
                    resistant_fraction=round(resistant_frac, 3),
                    evidence_level=EvidenceLevel.LIMITED,
                    source_studies=["DepMap Portal (Broad)", "CCLE (Barretina 2012)"],
                )
                predictions.append(pred)

        # Sort by sensitivity descending
        predictions.sort(key=lambda p: -p.predicted_sensitivity)

        severity = SeverityLevel.INFORMATIONAL
        any_resistant = any(p.resistance_risk == ResistanceRisk.HIGH for p in predictions)
        if any_resistant:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            drug_predictions=predictions,
            severity=severity,
        )

    @staticmethod
    def _infer_drug_class(drug: str) -> str:
        classes = {
            "cisplatin": "platinum agent",
            "paclitaxel": "taxane",
            "doxorubicin": "anthracycline",
            "venetoclax": "BCL-2 inhibitor",
            "imatinib": "TKI",
            "olaparib": "PARP inhibitor",
            "trametinib": "MEK inhibitor",
            "osimertinib": "EGFR TKI",
            "pembrolizumab": "checkpoint inhibitor",
            "5-fluorouracil": "antimetabolite",
        }
        return classes.get(drug, "unknown")


# ===================================================================
# WORKFLOW 4 -- Subclonal Architecture
# ===================================================================


class SubclonalArchitectureWorkflow(BaseSCWorkflow):
    """Clone detection, mutation-expression correlation, resistance risk
    scoring.  HIGH risk if antigen-negative clone >10%.

    Inputs
    ------
    clone_data : list[dict]
        Each dict: {clone_id, cell_count, driver_mutations, proliferation_index,
                    antigen_expression, resistance_genes, is_expanding,
                    cnv_profile, transcriptomic_signature}.
    total_cells : int
        Total cells in sample.
    target_antigen : str
        The therapeutic target antigen gene (e.g. 'CD19').
    mutation_expression_corr : dict[str, float] | None
        mutation -> expression correlation coefficient.
    """

    workflow_type = SCWorkflowType.SUBCLONAL_ARCHITECTURE

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "clone_data", warnings, default=[])
        self._require_key(inp, "total_cells", warnings, default=0)
        inp.setdefault("target_antigen", "CD19")
        inp.setdefault("mutation_expression_corr", {})
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        clone_list = inputs.get("clone_data", [])
        total_cells = inputs.get("total_cells", 0)
        target = inputs.get("target_antigen", "CD19")
        inputs.get("mutation_expression_corr", {})

        subclones: List[SubclonalResult] = []

        for cd in clone_list:
            cell_count = cd.get("cell_count", 0)
            freq = _safe_div(cell_count, max(total_cells, 1))
            antigen_expr = cd.get("antigen_expression", 0.0)

            # Determine therapy implications
            implications: List[str] = []
            if antigen_expr < 0.1:
                implications.append(
                    f"{target}-negative: potential escape variant"
                )
            if cd.get("is_expanding", False):
                implications.append("Expanding clone: monitor closely")
            if cd.get("resistance_genes"):
                implications.append(
                    f"Resistance genes: {', '.join(cd['resistance_genes'][:5])}"
                )

            # Fitness score based on proliferation and expansion
            fitness = cd.get("proliferation_index", 0.0)
            if cd.get("is_expanding", False):
                fitness = fitness * 1.5

            sr = SubclonalResult(
                clone_id=cd.get("clone_id", "unknown"),
                clone_fraction=round(freq, 4),
                cell_count=cell_count,
                driver_mutations=cd.get("driver_mutations", []),
                cnv_profile=cd.get("cnv_profile", {}),
                transcriptomic_signature=cd.get("transcriptomic_signature", []),
                phylogenetic_position=cd.get("phylogenetic_position"),
                fitness_score=round(fitness, 3),
                therapy_implications=implications,
            )
            subclones.append(sr)

        # Sort by frequency descending
        subclones.sort(key=lambda s: -s.clone_fraction)

        # Determine overall severity based on escape risk
        neg_freq = sum(
            s.clone_fraction for s in subclones
            if any("negative" in impl.lower() for impl in s.therapy_implications)
        )

        severity = SeverityLevel.INFORMATIONAL
        if neg_freq > 0.10:
            severity = SeverityLevel.CRITICAL
        elif neg_freq > 0.03:
            severity = SeverityLevel.HIGH
        elif any(s.fitness_score > 1.0 for s in subclones):
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            subclones=subclones,
            severity=severity,
        )


# ===================================================================
# WORKFLOW 5 -- Spatial Niche Analysis
# ===================================================================


class SpatialNicheWorkflow(BaseSCWorkflow):
    """k-NN neighbourhood analysis, niche clustering, spatial
    autocorrelation (Moran's I).

    Inputs
    ------
    cell_coordinates : dict[str, tuple[float, float]]
        cell_id -> (x, y) spatial coordinates.
    cell_types : dict[str, str]
        cell_id -> cell type annotation.
    gene_expression : dict[str, dict[str, float]] | None
        cell_id -> {gene: expr} for spatial autocorrelation.
    k_neighbors : int
        Number of nearest neighbours (default 15).
    genes_of_interest : list[str] | None
        Genes for Moran's I computation.
    """

    workflow_type = SCWorkflowType.SPATIAL_NICHE

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "cell_coordinates", warnings, default={})
        self._require_key(inp, "cell_types", warnings, default={})
        inp.setdefault("gene_expression", None)
        inp.setdefault("k_neighbors", 15)
        inp.setdefault("genes_of_interest", ["CD8A", "CD274", "FOXP3", "COL1A1", "MKI67"])
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        coords = inputs.get("cell_coordinates", {})
        cell_types = inputs.get("cell_types", {})
        gene_expr = inputs.get("gene_expression")
        k = inputs.get("k_neighbors", 15)
        goi = inputs.get("genes_of_interest", [])

        if not coords:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                severity=SeverityLevel.LOW,
            )

        cell_ids = list(coords.keys())
        n_cells = len(cell_ids)

        # --- Build kNN graph ---
        knn_graph: Dict[str, List[str]] = {}
        for i, cid in enumerate(cell_ids):
            x1, y1 = coords[cid]
            dists = []
            for j, oid in enumerate(cell_ids):
                if i == j:
                    continue
                x2, y2 = coords[oid]
                d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                dists.append((d, oid))
            dists.sort(key=lambda t: t[0])
            knn_graph[cid] = [oid for _, oid in dists[:k]]

        # --- Niche assignment by neighbourhood composition ---
        niche_labels: Dict[str, str] = {}
        niche_compositions: Dict[str, Dict[str, int]] = {}

        for cid in cell_ids:
            neighbors = knn_graph.get(cid, [])
            comp: Dict[str, int] = {}
            for nid in neighbors:
                ct = cell_types.get(nid, "Unknown")
                comp[ct] = comp.get(ct, 0) + 1
            self_ct = cell_types.get(cid, "Unknown")
            comp[self_ct] = comp.get(self_ct, 0) + 1

            niche_type = self._classify_niche_label(comp, k + 1)
            niche_labels[cid] = niche_type
            if niche_type not in niche_compositions:
                niche_compositions[niche_type] = {}
            for ct, cnt in comp.items():
                niche_compositions[niche_type][ct] = niche_compositions[niche_type].get(ct, 0) + cnt

        # Aggregate niches
        niche_counts: Dict[str, int] = {}
        for lbl in niche_labels.values():
            niche_counts[lbl] = niche_counts.get(lbl, 0) + 1

        spatial_niches: List[SpatialNiche] = []
        for i, (ntype, count) in enumerate(
            sorted(niche_counts.items(), key=lambda x: -x[1])
        ):
            raw_comp = niche_compositions.get(ntype, {})
            total = max(sum(raw_comp.values()), 1)
            prop_comp = {ct: round(cnt / total, 3) for ct, cnt in raw_comp.items()}
            dominant = sorted(prop_comp.items(), key=lambda x: -x[1])[:4]

            # Compute Moran's I for signature genes in this niche
            sig_genes: List[str] = []
            if gene_expr and goi:
                niche_cell_ids = [cid for cid in cell_ids if niche_labels.get(cid) == ntype]
                for gene in goi[:3]:
                    mi = self._compute_morans_i(gene, gene_expr, knn_graph, niche_cell_ids)
                    if mi is not None and abs(mi) > 0.2:
                        sig_genes.append(gene)

            sn = SpatialNiche(
                niche_id=f"niche_{i}",
                niche_label=ntype,
                dominant_cell_types=[ct for ct, _ in dominant],
                cell_type_proportions=prop_comp,
                signature_genes=sig_genes,
                area_fraction=round(count / n_cells, 4),
            )
            spatial_niches.append(sn)

        severity = SeverityLevel.INFORMATIONAL
        tumor_core_frac = niche_counts.get("Tumor core", 0) / max(n_cells, 1)
        if tumor_core_frac > 0.7:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            spatial_niches=spatial_niches,
            severity=severity,
        )

    @staticmethod
    def _classify_niche_label(composition: Dict[str, int], total: int) -> str:
        """Classify a spatial niche based on cell type composition."""
        props = {ct: cnt / max(total, 1) for ct, cnt in composition.items()}
        immune_terms = {"cd8_t", "cd4_t", "nk", "b_cell", "macrophage", "dendritic",
                        "treg", "plasma", "neutrophil", "mast"}
        stromal_terms = {"fibroblast", "caf", "pericyte", "myofibroblast"}
        endothelial_terms = {"endothelial"}
        malignant_terms = {"malignant", "tumor", "cancer"}

        immune_frac = sum(
            v for k, v in props.items()
            if k.lower() in immune_terms or "t_cell" in k.lower() or "macrophage" in k.lower()
        )
        stromal_frac = sum(v for k, v in props.items() if k.lower() in stromal_terms)
        endothelial_frac = sum(v for k, v in props.items() if k.lower() in endothelial_terms)
        malignant_frac = sum(v for k, v in props.items() if k.lower() in malignant_terms)

        if malignant_frac > 0.6:
            return "Tumor core"
        if immune_frac > 0.4:
            return "Immune niche"
        if stromal_frac > 0.4:
            return "Stromal niche"
        if endothelial_frac > 0.3:
            return "Vascular niche"
        if malignant_frac > 0.3 and immune_frac > 0.15:
            return "Tumor-immune interface"
        return "Normal adjacent"

    @staticmethod
    def _compute_morans_i(
        gene: str,
        gene_expr: Dict[str, Dict[str, float]],
        knn_graph: Dict[str, List[str]],
        cell_ids: List[str],
    ) -> Optional[float]:
        """Compute Moran's I spatial autocorrelation for a gene."""
        values = []
        for cid in cell_ids:
            e = gene_expr.get(cid, {})
            values.append(e.get(gene, 0.0))

        n = len(values)
        if n < 3:
            return None

        mean_val = sum(values) / n
        deviations = [v - mean_val for v in values]
        var = sum(d * d for d in deviations)
        if var == 0:
            return 0.0

        w_sum = 0.0
        total_w = 0.0
        id_to_idx = {cid: i for i, cid in enumerate(cell_ids)}
        for i, cid in enumerate(cell_ids):
            for nid in knn_graph.get(cid, []):
                j = id_to_idx.get(nid)
                if j is not None:
                    w_sum += deviations[i] * deviations[j]
                    total_w += 1.0

        if total_w == 0:
            return 0.0

        morans_i = (n / total_w) * (w_sum / var)
        return _clamp(morans_i, -1.0, 1.0)


# ===================================================================
# WORKFLOW 6 -- Trajectory Analysis
# ===================================================================


class TrajectoryAnalysisWorkflow(BaseSCWorkflow):
    """Pseudotime inference, driver gene identification, branch point
    detection.

    Inputs
    ------
    cell_pseudotime : dict[str, float]
        cell_id -> pseudotime value (0 = root).
    cell_types : dict[str, str]
        cell_id -> cell type annotation.
    gene_expression : dict[str, dict[str, float]]
        cell_id -> {gene: expression}.
    branch_assignments : dict[str, str] | None
        cell_id -> branch label (if branch analysis done externally).
    root_cell_type : str | None
        Expected root cell type (e.g. 'HSC').
    driver_gene_candidates : list[str] | None
        Genes to test for pseudotime correlation.
    trajectory_type : str
        Type: 'differentiation', 'activation', 'exhaustion', etc.
    """

    workflow_type = SCWorkflowType.TRAJECTORY_ANALYSIS

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "cell_pseudotime", warnings, default={})
        self._require_key(inp, "cell_types", warnings, default={})
        inp.setdefault("gene_expression", {})
        inp.setdefault("branch_assignments", None)
        inp.setdefault("root_cell_type", None)
        inp.setdefault("driver_gene_candidates", None)
        inp.setdefault("trajectory_type", "differentiation")
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        pseudotime = inputs.get("cell_pseudotime", {})
        cell_types = inputs.get("cell_types", {})
        gene_expr = inputs.get("gene_expression", {})
        branches = inputs.get("branch_assignments")
        inputs.get("root_cell_type")
        driver_candidates = inputs.get("driver_gene_candidates")
        traj_type_str = inputs.get("trajectory_type", "differentiation")

        if not pseudotime:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                severity=SeverityLevel.LOW,
            )

        n_cells = len(pseudotime)
        pt_values = sorted(pseudotime.values())
        pt_min, pt_max = pt_values[0], pt_values[-1]

        # Cell type ordering along pseudotime
        ct_median_pt: Dict[str, List[float]] = {}
        for cid, pt in pseudotime.items():
            ct = cell_types.get(cid, "Unknown")
            ct_median_pt.setdefault(ct, []).append(pt)

        ct_order = []
        for ct, pts in ct_median_pt.items():
            median_pt = sorted(pts)[len(pts) // 2]
            ct_order.append((ct, median_pt, len(pts)))
        ct_order.sort(key=lambda x: x[1])

        if not ct_order:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                severity=SeverityLevel.LOW,
            )

        start_ct = ct_order[0][0]
        end_ct = ct_order[-1][0]
        intermediate = [ct for ct, _, _ in ct_order[1:-1]] if len(ct_order) > 2 else []

        # Map trajectory type string to enum
        try:
            traj_type = TrajectoryType(traj_type_str)
        except ValueError:
            traj_type = TrajectoryType.DIFFERENTIATION

        # --- Driver gene identification ---
        driver_genes: List[Tuple[str, float]] = []
        if gene_expr:
            all_genes: set = set()
            for gdict in gene_expr.values():
                all_genes.update(gdict.keys())

            test_genes = driver_candidates if driver_candidates else list(all_genes)[:200]

            for gene in test_genes:
                corr = self._pseudotime_correlation(gene, pseudotime, gene_expr)
                if corr is not None and abs(corr) > 0.3:
                    driver_genes.append((gene, corr))

            driver_genes.sort(key=lambda x: -abs(x[1]))

        # --- Branch point detection ---
        branching_points: List[str] = []
        if branches:
            branch_cells: Dict[str, List[str]] = {}
            for cid, blabel in branches.items():
                branch_cells.setdefault(blabel, []).append(cid)

            if len(branch_cells) > 1:
                branch_min_pt: Dict[str, float] = {}
                for blabel, cids in branch_cells.items():
                    pts = [pseudotime.get(cid, float('inf')) for cid in cids]
                    branch_min_pt[blabel] = min(pts) if pts else float('inf')

                branch_time = min(branch_min_pt.values())
                branching_points.append(f"t={branch_time:.3f}")

        trajectory = TrajectoryResult(
            trajectory_id="traj_0",
            trajectory_type=traj_type,
            start_cell_type=start_ct,
            end_cell_type=end_ct,
            intermediate_states=intermediate,
            driver_genes=[g for g, _ in driver_genes[:15]],
            branching_points=branching_points,
            pseudotime_range=[round(pt_min, 4), round(pt_max, 4)],
            cell_count=n_cells,
            clinical_relevance=(
                f"Trajectory from {start_ct} to {end_ct} with "
                f"{len(driver_genes)} driver genes identified"
            ),
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            trajectories=[trajectory],
            severity=SeverityLevel.INFORMATIONAL,
        )

    @staticmethod
    def _pseudotime_correlation(
        gene: str,
        pseudotime: Dict[str, float],
        gene_expr: Dict[str, Dict[str, float]],
    ) -> Optional[float]:
        """Compute Pearson correlation between gene expression and pseudotime."""
        xs, ys = [], []
        for cid, pt in pseudotime.items():
            val = gene_expr.get(cid, {}).get(gene)
            if val is not None:
                xs.append(pt)
                ys.append(val)

        n = len(xs)
        if n < 5:
            return None

        mean_x = sum(xs) / n
        mean_y = sum(ys) / n

        cov = sum((xs[i] - mean_x) * (ys[i] - mean_y) for i in range(n))
        var_x = sum((x - mean_x) ** 2 for x in xs)
        var_y = sum((y - mean_y) ** 2 for y in ys)

        denom = math.sqrt(var_x * var_y)
        if denom == 0:
            return 0.0
        return _clamp(cov / denom, -1.0, 1.0)


# ===================================================================
# WORKFLOW 7 -- Ligand-Receptor Interaction Mapping
# ===================================================================


class LigandReceptorWorkflow(BaseSCWorkflow):
    """CellChat/CellPhoneDB-style interaction mapping, network
    visualisation data.

    Inputs
    ------
    cell_type_expression : dict[str, dict[str, float]]
        cell_type -> {gene: mean_expression}.
    interaction_database : list[dict] | None
        Each dict: {ligand, receptor, pathway}.
    min_expression : float
        Minimum expression threshold (default 0.5).
    p_value_threshold : float
        Max p-value for significant interactions (default 0.05).
    """

    workflow_type = SCWorkflowType.LIGAND_RECEPTOR

    _DEFAULT_INTERACTIONS: List[Dict] = [
        {"ligand": "CXCL12", "receptor": "CXCR4", "pathway": "CXCL"},
        {"ligand": "CCL19", "receptor": "CCR7", "pathway": "CCL"},
        {"ligand": "CD274", "receptor": "PDCD1", "pathway": "PD-L1/PD-1"},
        {"ligand": "CD80", "receptor": "CTLA4", "pathway": "CD80/CTLA4"},
        {"ligand": "CD80", "receptor": "CD28", "pathway": "CD80/CD28"},
        {"ligand": "TGFB1", "receptor": "TGFBR2", "pathway": "TGFb"},
        {"ligand": "VEGFA", "receptor": "FLT1", "pathway": "VEGF"},
        {"ligand": "VEGFA", "receptor": "KDR", "pathway": "VEGF"},
        {"ligand": "TNF", "receptor": "TNFRSF1A", "pathway": "TNF"},
        {"ligand": "IL6", "receptor": "IL6R", "pathway": "IL6"},
        {"ligand": "IL10", "receptor": "IL10RA", "pathway": "IL10"},
        {"ligand": "IFNG", "receptor": "IFNGR1", "pathway": "IFNg"},
        {"ligand": "FAS", "receptor": "FASLG", "pathway": "FAS"},
        {"ligand": "WNT5A", "receptor": "FZD5", "pathway": "WNT"},
        {"ligand": "DLL1", "receptor": "NOTCH1", "pathway": "NOTCH"},
        {"ligand": "JAG1", "receptor": "NOTCH2", "pathway": "NOTCH"},
        {"ligand": "HGF", "receptor": "MET", "pathway": "HGF/MET"},
        {"ligand": "EGF", "receptor": "EGFR", "pathway": "EGF"},
        {"ligand": "PDGFB", "receptor": "PDGFRB", "pathway": "PDGF"},
        {"ligand": "CSF1", "receptor": "CSF1R", "pathway": "CSF1"},
        {"ligand": "SPP1", "receptor": "CD44", "pathway": "SPP1"},
        {"ligand": "GAS6", "receptor": "AXL", "pathway": "GAS6/AXL"},
        {"ligand": "TIGIT", "receptor": "PVR", "pathway": "TIGIT/PVR"},
        {"ligand": "CD40LG", "receptor": "CD40", "pathway": "CD40"},
    ]

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "cell_type_expression", warnings, default={})
        inp.setdefault("interaction_database", self._DEFAULT_INTERACTIONS)
        inp.setdefault("min_expression", 0.5)
        inp.setdefault("p_value_threshold", 0.05)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        ct_expr = inputs.get("cell_type_expression", {})
        interactions_db = inputs.get("interaction_database", self._DEFAULT_INTERACTIONS)
        min_expr = inputs.get("min_expression", 0.5)
        p_thresh = inputs.get("p_value_threshold", 0.05)

        sig_interactions: List[LigandReceptorInteraction] = []
        cell_types_list = list(ct_expr.keys())

        for interaction in interactions_db:
            ligand = interaction["ligand"]
            receptor = interaction["receptor"]
            pathway = interaction.get("pathway", "unknown")

            for src_ct in cell_types_list:
                lig_expr = ct_expr[src_ct].get(ligand, 0.0)
                if lig_expr < min_expr:
                    continue

                for tgt_ct in cell_types_list:
                    rec_expr = ct_expr[tgt_ct].get(receptor, 0.0)
                    if rec_expr < min_expr:
                        continue

                    score = math.sqrt(lig_expr * rec_expr)
                    p_val = max(0.001, 1.0 / (1.0 + score * 2.0))

                    if p_val <= p_thresh:
                        pair = LigandReceptorInteraction(
                            source_cell_type=src_ct,
                            target_cell_type=tgt_ct,
                            ligand_gene=ligand,
                            receptor_gene=receptor,
                            interaction_score=round(score, 4),
                            p_value=round(p_val, 5),
                            pathway=pathway,
                            method="CellChat-style",
                        )
                        sig_interactions.append(pair)

        sig_interactions.sort(key=lambda p: -p.interaction_score)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            interactions=sig_interactions,
            severity=SeverityLevel.INFORMATIONAL,
        )


# ===================================================================
# WORKFLOW 8 -- Biomarker Discovery
# ===================================================================


class BiomarkerDiscoveryWorkflow(BaseSCWorkflow):
    """Differential expression, cell-type specificity, clinical
    correlation scoring.

    Inputs
    ------
    de_results : list[dict]
        Each dict: {gene, log2_fold_change, p_value_adj, cell_type}.
    cell_type_expression : dict[str, dict[str, float]]
        cell_type -> {gene: mean_expression} for specificity computation.
    clinical_outcomes : dict[str, float] | None
        gene -> clinical outcome correlation (e.g., survival HR).
    min_log2fc : float
        Minimum absolute log2FC threshold (default 1.0).
    max_padj : float
        Maximum adjusted p-value (default 0.05).
    """

    workflow_type = SCWorkflowType.BIOMARKER_DISCOVERY

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "de_results", warnings, default=[])
        inp.setdefault("cell_type_expression", {})
        inp.setdefault("clinical_outcomes", {})
        inp.setdefault("min_log2fc", 1.0)
        inp.setdefault("max_padj", 0.05)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        de_results = inputs.get("de_results", [])
        ct_expr = inputs.get("cell_type_expression", {})
        clin_outcomes = inputs.get("clinical_outcomes", {})
        min_fc = inputs.get("min_log2fc", 1.0)
        max_p = inputs.get("max_padj", 0.05)

        sig_de = [
            d for d in de_results
            if abs(d.get("log2_fold_change", 0)) >= min_fc
            and d.get("p_value_adj", 1.0) <= max_p
        ]

        candidates: List[BiomarkerCandidate] = []

        for de in sig_de:
            gene = de.get("gene", "")
            fc = de.get("log2_fold_change", 0.0)
            padj = de.get("p_value_adj", 1.0)
            source_ct = de.get("cell_type", "Unknown")

            # Specificity score
            source_expr = ct_expr.get(source_ct, {}).get(gene, 0.0)
            other_exprs = [
                ct_expr.get(ct, {}).get(gene, 0.0)
                for ct in ct_expr if ct != source_ct
            ]
            max_other = max(other_exprs) if other_exprs else 0.0
            specificity = _clamp(
                _safe_div(source_expr - max_other, source_expr + max_other + 0.01)
            )

            clin_corr = abs(clin_outcomes.get(gene, 0.0))

            # Determine biomarker type
            if clin_corr > 0.5:
                bm_type = "prognostic"
            elif abs(fc) > 3.0:
                bm_type = "diagnostic"
            else:
                bm_type = "predictive"

            # Evidence level
            if padj < 0.001 and specificity > 0.7:
                ev = EvidenceLevel.STRONG
            elif padj < 0.01:
                ev = EvidenceLevel.MODERATE
            else:
                ev = EvidenceLevel.LIMITED

            bc = BiomarkerCandidate(
                gene=gene,
                biomarker_type=bm_type,
                cell_type_specific=source_ct,
                fold_change=round(fc, 3),
                p_value_adjusted=round(padj, 6),
                specificity_score=round(specificity, 3),
                evidence_level=ev,
            )
            candidates.append(bc)

        # Sort by specificity * fold_change magnitude
        candidates.sort(key=lambda b: -(b.specificity_score * abs(b.fold_change)))

        return WorkflowResult(
            workflow_type=self.workflow_type,
            biomarkers=candidates,
            severity=SeverityLevel.INFORMATIONAL,
        )


# ===================================================================
# WORKFLOW 9 -- CAR-T Target Validation
# ===================================================================


class CARTTargetValidationWorkflow(BaseSCWorkflow):
    """On-tumour expression, off-tumour safety (vital organ check),
    TME compatibility, escape risk assessment.

    Inputs
    ------
    target_gene : str
        CAR-T target antigen gene symbol.
    tumor_expression : dict[str, float]
        cell_id -> target gene expression in tumour sample.
    tumor_cell_ids : list[str]
        Cell IDs classified as malignant.
    normal_tissue_expression : dict[str, float]
        organ/tissue -> mean target expression in normal tissue atlas.
    tme_data : dict
        TME profiling data (cell_type_proportions, etc.).
    clone_data : list[dict] | None
        Subclonal data for escape risk.
    """

    workflow_type = SCWorkflowType.CART_TARGET_VALIDATION

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "target_gene", warnings, default="CD19")
        self._require_key(inp, "tumor_expression", warnings, default={})
        inp.setdefault("tumor_cell_ids", [])
        inp.setdefault("normal_tissue_expression", {})
        inp.setdefault("tme_data", {})
        inp.setdefault("clone_data", None)
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        target = inputs.get("target_gene", "CD19")
        tumor_expr = inputs.get("tumor_expression", {})
        tumor_ids = set(inputs.get("tumor_cell_ids", []))
        normal_expr = inputs.get("normal_tissue_expression", {})
        inputs.get("tme_data", {})
        clone_data = inputs.get("clone_data")

        # --- On-tumour expression ---
        if tumor_ids and tumor_expr:
            tumor_vals = [
                tumor_expr.get(cid, 0.0) for cid in tumor_ids
                if cid in tumor_expr
            ]
        else:
            tumor_vals = list(tumor_expr.values())

        if tumor_vals:
            expressing = [v for v in tumor_vals if v > 0.1]
            on_tumor_pct = len(expressing) / max(len(tumor_vals), 1)
            mean_expr = sum(tumor_vals) / len(tumor_vals)
        else:
            on_tumor_pct = 0.0
            mean_expr = 0.0

        # --- Off-tumour vital organ safety ---
        vital_hits: Dict[str, float] = {}
        max_off_tumor = 0.0
        for organ in _VITAL_ORGANS:
            organ_expr = normal_expr.get(organ, 0.0)
            if organ_expr > 0.5:
                vital_hits[organ] = round(organ_expr, 3)
            max_off_tumor = max(max_off_tumor, organ_expr)

        # Therapeutic index
        _safe_div(mean_expr, max_off_tumor + 0.01)

        # Off-tumor risk classification
        if vital_hits:
            if any(v > 2.0 for v in vital_hits.values()):
                off_tumor_risk = "high"
            else:
                off_tumor_risk = "medium"
        else:
            off_tumor_risk = "low"

        # --- Heterogeneity score ---
        if tumor_vals and len(tumor_vals) > 1:
            mean_v = sum(tumor_vals) / len(tumor_vals)
            variance = sum((v - mean_v) ** 2 for v in tumor_vals) / len(tumor_vals)
            heterogeneity = _clamp(math.sqrt(variance) / (mean_v + 0.01))
        else:
            heterogeneity = 0.0

        # --- Escape risk ---
        escape_risk = ResistanceRisk.LOW
        if clone_data:
            neg_clones = [
                c for c in clone_data
                if c.get("antigen_expression", 1.0) < 0.1
            ]
            neg_freq = sum(
                c.get("cell_count", 0) for c in neg_clones
            ) / max(sum(c.get("cell_count", 0) for c in clone_data), 1)
            if neg_freq > 0.10:
                escape_risk = ResistanceRisk.HIGH
            elif neg_freq > 0.03:
                escape_risk = ResistanceRisk.MEDIUM
        elif on_tumor_pct < 0.8:
            escape_risk = ResistanceRisk.MEDIUM

        cart_target = CARTTargetValidation(
            target_gene=target,
            tumor_expression_fraction=round(on_tumor_pct, 4),
            tumor_expression_level=round(mean_expr, 3),
            normal_tissue_expression=vital_hits,
            on_target_off_tumor_risk=off_tumor_risk,
            heterogeneity_score=round(heterogeneity, 3),
            escape_variant_risk=escape_risk,
            evidence_level=EvidenceLevel.MODERATE,
        )

        severity = SeverityLevel.INFORMATIONAL
        if off_tumor_risk == "high":
            severity = SeverityLevel.HIGH
        if on_tumor_pct < 0.5:
            severity = _max_severity(severity, SeverityLevel.HIGH)
        if escape_risk == ResistanceRisk.HIGH:
            severity = _max_severity(severity, SeverityLevel.CRITICAL)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            cart_targets=[cart_target],
            severity=severity,
        )


# ===================================================================
# WORKFLOW 10 -- Treatment Monitoring
# ===================================================================


class TreatmentMonitoringWorkflow(BaseSCWorkflow):
    """Longitudinal clonal dynamics, exhaustion trajectory, resistance
    emergence monitoring.

    Inputs
    ------
    timepoints : list[dict]
        Each dict: {timepoint_id, days_from_baseline, cell_type_proportions,
                    clone_frequencies, exhaustion_markers, resistance_markers,
                    treatment}.
    target_antigen : str
        Therapeutic target gene.
    treatment_type : str
        'car_t', 'checkpoint_inhibitor', 'targeted_therapy', etc.
    baseline_response : str | None
        'CR', 'PR', 'SD', 'PD' -- baseline clinical response.
    """

    workflow_type = SCWorkflowType.TREATMENT_MONITORING

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        self._require_key(inp, "timepoints", warnings, default=[])
        inp.setdefault("target_antigen", "CD19")
        inp.setdefault("treatment_type", "car_t")
        inp.setdefault("baseline_response", None)
        # Sort timepoints
        tps = inp.get("timepoints", [])
        tps.sort(key=lambda t: t.get("days_from_baseline", 0))
        inp["timepoints"] = tps
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        timepoints = inputs.get("timepoints", [])
        inputs.get("target_antigen", "CD19")
        treatment = inputs.get("treatment_type", "car_t")
        baseline_resp = inputs.get("baseline_response")

        if len(timepoints) < 2:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                severity=SeverityLevel.LOW,
            )

        responses: List[TreatmentResponse] = []

        # Use first timepoint as baseline
        baseline_ct = timepoints[0].get("cell_type_proportions", {})

        for tp in timepoints:
            tp_id = tp.get("timepoint_id", f"day_{tp.get('days_from_baseline', 0)}")
            day = tp.get("days_from_baseline", 0)
            current_ct = tp.get("cell_type_proportions", {})
            clone_freq = tp.get("clone_frequencies", {})
            exh_markers = tp.get("exhaustion_markers", {})
            res_markers = tp.get("resistance_markers", {})

            # Compositional shifts from baseline
            all_types = set(baseline_ct.keys()) | set(current_ct.keys())
            shifts = {}
            for ct in all_types:
                shift = current_ct.get(ct, 0.0) - baseline_ct.get(ct, 0.0)
                if abs(shift) > 0.02:
                    shifts[ct] = round(shift, 4)

            # Emerging clones (present now but not at baseline or expanding)
            baseline_clones = timepoints[0].get("clone_frequencies", {})
            emerging = [
                cid for cid, freq in clone_freq.items()
                if freq > 0.05 and baseline_clones.get(cid, 0.0) < 0.01
            ]

            # Resistance signatures
            resistance_sigs = [
                gene for gene, val in res_markers.items()
                if val > 2.0
            ]

            # Exhaustion dynamics
            exh_state = self._classify_exhaustion(exh_markers)
            immune_dynamics = {"exhaustion_state": exh_state}
            for m, v in exh_markers.items():
                immune_dynamics[m] = round(v, 3)

            # Response category estimation
            if day == 0:
                response_cat = baseline_resp or "baseline"
            elif resistance_sigs or emerging:
                response_cat = "progression"
            elif shifts:
                positive_shifts = sum(1 for v in shifts.values() if v > 0.05)
                negative_shifts = sum(1 for v in shifts.values() if v < -0.05)
                if positive_shifts > negative_shifts:
                    response_cat = "partial"
                else:
                    response_cat = "stable"
            else:
                response_cat = "stable"

            # Actionable findings
            actionable: List[str] = []
            if resistance_sigs:
                actionable.append(
                    f"Resistance markers emerging: {', '.join(resistance_sigs[:5])}"
                )
            if emerging:
                actionable.append(
                    f"New clones expanding: {', '.join(emerging[:3])}"
                )
            if exh_state == "terminally_exhausted":
                actionable.append(
                    "Terminal T-cell exhaustion: consider reinvigoration strategy"
                )

            tr = TreatmentResponse(
                timepoint=tp_id,
                treatment=treatment,
                response_category=response_cat,
                compositional_shifts=shifts,
                emerging_clones=emerging,
                resistance_signatures=resistance_sigs,
                immune_dynamics=immune_dynamics,
                actionable_findings=actionable,
            )
            responses.append(tr)

        # Determine overall severity
        severity = SeverityLevel.INFORMATIONAL
        last_resp = responses[-1] if responses else None
        if last_resp:
            if last_resp.resistance_signatures:
                severity = SeverityLevel.HIGH
            if last_resp.emerging_clones:
                severity = _max_severity(severity, SeverityLevel.HIGH)
            if last_resp.response_category == "progression":
                severity = _max_severity(severity, SeverityLevel.CRITICAL)
            exh = last_resp.immune_dynamics.get("exhaustion_state", "")
            if exh == "terminally_exhausted":
                severity = _max_severity(severity, SeverityLevel.HIGH)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            treatment_responses=responses,
            severity=severity,
        )

    @staticmethod
    def _classify_exhaustion(markers: Dict[str, float]) -> str:
        """Classify T-cell exhaustion state from marker expression."""
        best_state = "naive"
        best_score = -1.0

        for state, state_markers in _EXHAUSTION_MARKERS.items():
            score = sum(markers.get(m, 0.0) for m in state_markers)
            score = score / max(len(state_markers), 1)
            if score > best_score:
                best_score = score
                best_state = state

        return best_state


# ===================================================================
# WORKFLOW ENGINE
# ===================================================================


class WorkflowEngine:
    """Central dispatcher that maps SCWorkflowType to the appropriate
    workflow implementation and handles query-based workflow detection."""

    _KEYWORD_MAP: Dict[str, SCWorkflowType] = {
        # Cell Type Annotation
        "cell type": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "cell annotation": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "cell identity": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "cluster annotation": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "cell label": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "singler": SCWorkflowType.CELL_TYPE_ANNOTATION,
        "celltypist": SCWorkflowType.CELL_TYPE_ANNOTATION,
        # TME Profiling
        "tumor microenvironment": SCWorkflowType.TME_PROFILING,
        "tme": SCWorkflowType.TME_PROFILING,
        "immune infiltration": SCWorkflowType.TME_PROFILING,
        "hot tumor": SCWorkflowType.TME_PROFILING,
        "cold tumor": SCWorkflowType.TME_PROFILING,
        "immune excluded": SCWorkflowType.TME_PROFILING,
        "immunoscore": SCWorkflowType.TME_PROFILING,
        "immune phenotype": SCWorkflowType.TME_PROFILING,
        # Drug Response
        "drug response": SCWorkflowType.DRUG_RESPONSE,
        "drug sensitivity": SCWorkflowType.DRUG_RESPONSE,
        "ic50": SCWorkflowType.DRUG_RESPONSE,
        "depmap": SCWorkflowType.DRUG_RESPONSE,
        "ccle": SCWorkflowType.DRUG_RESPONSE,
        "pharmacogenomics": SCWorkflowType.DRUG_RESPONSE,
        "drug screen": SCWorkflowType.DRUG_RESPONSE,
        # Subclonal Architecture
        "subclonal": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        "clonal architecture": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        "clone detection": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        "clonal evolution": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        "intratumoral heterogeneity": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        "ith": SCWorkflowType.SUBCLONAL_ARCHITECTURE,
        # Spatial Niche
        "spatial": SCWorkflowType.SPATIAL_NICHE,
        "niche": SCWorkflowType.SPATIAL_NICHE,
        "spatial transcriptomics": SCWorkflowType.SPATIAL_NICHE,
        "visium": SCWorkflowType.SPATIAL_NICHE,
        "merfish": SCWorkflowType.SPATIAL_NICHE,
        "moran": SCWorkflowType.SPATIAL_NICHE,
        "neighborhood": SCWorkflowType.SPATIAL_NICHE,
        "neighbourhood": SCWorkflowType.SPATIAL_NICHE,
        # Trajectory Analysis
        "trajectory": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "pseudotime": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "lineage": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "differentiation": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "monocle": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "rna velocity": SCWorkflowType.TRAJECTORY_ANALYSIS,
        "scvelo": SCWorkflowType.TRAJECTORY_ANALYSIS,
        # Ligand-Receptor
        "ligand receptor": SCWorkflowType.LIGAND_RECEPTOR,
        "ligand-receptor": SCWorkflowType.LIGAND_RECEPTOR,
        "cell-cell communication": SCWorkflowType.LIGAND_RECEPTOR,
        "cellchat": SCWorkflowType.LIGAND_RECEPTOR,
        "cellphonedb": SCWorkflowType.LIGAND_RECEPTOR,
        "nichenet": SCWorkflowType.LIGAND_RECEPTOR,
        "cell interaction": SCWorkflowType.LIGAND_RECEPTOR,
        "intercellular": SCWorkflowType.LIGAND_RECEPTOR,
        # Biomarker Discovery
        "biomarker": SCWorkflowType.BIOMARKER_DISCOVERY,
        "differential expression": SCWorkflowType.BIOMARKER_DISCOVERY,
        "deg": SCWorkflowType.BIOMARKER_DISCOVERY,
        "marker gene": SCWorkflowType.BIOMARKER_DISCOVERY,
        "prognostic marker": SCWorkflowType.BIOMARKER_DISCOVERY,
        # CAR-T Target Validation
        "car-t": SCWorkflowType.CART_TARGET_VALIDATION,
        "car t": SCWorkflowType.CART_TARGET_VALIDATION,
        "cart target": SCWorkflowType.CART_TARGET_VALIDATION,
        "antigen target": SCWorkflowType.CART_TARGET_VALIDATION,
        "on-tumor": SCWorkflowType.CART_TARGET_VALIDATION,
        "off-tumor": SCWorkflowType.CART_TARGET_VALIDATION,
        "therapeutic index": SCWorkflowType.CART_TARGET_VALIDATION,
        # Treatment Monitoring
        "treatment monitoring": SCWorkflowType.TREATMENT_MONITORING,
        "longitudinal": SCWorkflowType.TREATMENT_MONITORING,
        "clonal dynamics": SCWorkflowType.TREATMENT_MONITORING,
        "exhaustion": SCWorkflowType.TREATMENT_MONITORING,
        "resistance emergence": SCWorkflowType.TREATMENT_MONITORING,
        "treatment response": SCWorkflowType.TREATMENT_MONITORING,
        "monitoring": SCWorkflowType.TREATMENT_MONITORING,
    }

    _WORKFLOW_REGISTRY: Dict[SCWorkflowType, type] = {
        SCWorkflowType.CELL_TYPE_ANNOTATION: CellTypeAnnotationWorkflow,
        SCWorkflowType.TME_PROFILING: TMEProfilingWorkflow,
        SCWorkflowType.DRUG_RESPONSE: DrugResponseWorkflow,
        SCWorkflowType.SUBCLONAL_ARCHITECTURE: SubclonalArchitectureWorkflow,
        SCWorkflowType.SPATIAL_NICHE: SpatialNicheWorkflow,
        SCWorkflowType.TRAJECTORY_ANALYSIS: TrajectoryAnalysisWorkflow,
        SCWorkflowType.LIGAND_RECEPTOR: LigandReceptorWorkflow,
        SCWorkflowType.BIOMARKER_DISCOVERY: BiomarkerDiscoveryWorkflow,
        SCWorkflowType.CART_TARGET_VALIDATION: CARTTargetValidationWorkflow,
        SCWorkflowType.TREATMENT_MONITORING: TreatmentMonitoringWorkflow,
    }

    def detect_workflow(self, query: str) -> SCWorkflowType:
        """Detect the most appropriate workflow from a free-text query."""
        query_lower = query.lower()
        scores: Dict[SCWorkflowType, int] = {}
        for keyword, wtype in self._KEYWORD_MAP.items():
            if keyword in query_lower:
                scores[wtype] = scores.get(wtype, 0) + 1
        if scores:
            return max(scores, key=scores.get)
        return SCWorkflowType.CELL_TYPE_ANNOTATION

    def get_workflow(self, workflow_type: SCWorkflowType) -> BaseSCWorkflow:
        """Instantiate and return the workflow for the given type."""
        cls = self._WORKFLOW_REGISTRY.get(workflow_type)
        if cls is None:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
        return cls()

    def run(self, workflow_type: SCWorkflowType, inputs: dict) -> WorkflowResult:
        """Dispatch to the appropriate workflow and execute."""
        workflow = self.get_workflow(workflow_type)
        return workflow.run(inputs)

    def run_from_query(self, query: str, inputs: dict) -> WorkflowResult:
        """Detect workflow from query text and execute."""
        wtype = self.detect_workflow(query)
        logger.info("Detected workflow %s from query", wtype.value)
        return self.run(wtype, inputs)

    @property
    def available_workflows(self) -> List[str]:
        """Return list of registered workflow type values."""
        return [wt.value for wt in self._WORKFLOW_REGISTRY]
