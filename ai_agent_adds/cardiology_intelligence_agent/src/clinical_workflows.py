"""Clinical workflows for the Cardiology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements eight evidence-based clinical workflows that integrate imaging,
biomarker, electrophysiology, and genomic data to produce actionable
cardiology assessments.  Each workflow follows the BaseCardioWorkflow
contract (preprocess -> execute -> postprocess) and is registered in the
WorkflowEngine for unified dispatch.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Dict, List, Optional

from src.models import (
    AnticoagulationRecommendation,
    CardiotoxicityRisk,
    CardioWorkflowType,
    EjectionFractionCategory,
    HeartFailureClass,
    HeartFailureStage,
    RiskScoreResult,
    RiskScoreType,
    SeverityLevel,
    ValveSeverity,
    WorkflowResult,
)

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

# Ordered severity levels for max() comparisons (lowest -> highest).
_SEVERITY_ORDER: List[SeverityLevel] = [
    SeverityLevel.INFORMATIONAL,
    SeverityLevel.BORDERLINE,
    SeverityLevel.LOW,
    SeverityLevel.INTERMEDIATE,
    SeverityLevel.MODERATE,
    SeverityLevel.HIGH,
    SeverityLevel.VERY_HIGH,
    SeverityLevel.CRITICAL,
]


def _max_severity(*levels: SeverityLevel) -> SeverityLevel:
    """Return the highest severity among the given levels."""
    return max(levels, key=lambda s: _SEVERITY_ORDER.index(s))


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _trigger_string(trigger_type: str, genes: List[str], reason: str) -> str:
    """Build a human-readable cross-modal trigger string."""
    gene_str = ", ".join(genes[:8])
    if len(genes) > 8:
        gene_str += f" (+{len(genes) - 8} more)"
    return f"[{trigger_type}] Genes: {gene_str} — {reason}"


# ═══════════════════════════════════════════════════════════════════════════
# BASE CLASS
# ═══════════════════════════════════════════════════════════════════════════


class BaseCardioWorkflow(ABC):
    """Abstract base for all cardiology clinical workflows."""

    workflow_type: CardioWorkflowType

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
        """Core clinical logic.  Must be implemented by each workflow."""
        ...

    def postprocess(self, result: WorkflowResult) -> WorkflowResult:
        """Shared enrichment after execution.  Override to add workflow-specific post-steps."""
        return result

    @staticmethod
    def _init_warnings(inp: dict) -> list:
        """Initialise and return the validation warnings list on *inp*."""
        warnings: list = inp.setdefault("_validation_warnings", [])
        return warnings


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 1 — CAD Assessment
# ═══════════════════════════════════════════════════════════════════════════


class CADAssessmentWorkflow(BaseCardioWorkflow):
    """Coronary artery disease assessment integrating calcium scoring,
    CAD-RADS grading, and high-risk plaque feature analysis.

    Inputs
    ------
    calcium_score : int
        Agatston coronary artery calcium score.
    cad_rads : str
        CAD-RADS grade (0, 1, 2, 3, 4A, 4B, 5).
    plaque_features : list[str]
        High-risk plaque descriptors (e.g. 'low_attenuation',
        'positive_remodeling', 'napkin_ring_sign', 'spotty_calcification').
    symptoms : list[str]
        Current cardiac symptoms.
    risk_factors : list[str]
        Patient cardiovascular risk factors.
    """

    workflow_type = CardioWorkflowType.CAD_ASSESSMENT

    # CAD-RADS management map
    _CADRADS_MGMT: Dict[str, dict] = {
        "0": {
            "stenosis": "0% — no plaque or stenosis",
            "management": "No further testing; preventive care",
            "severity": SeverityLevel.LOW,
        },
        "1": {
            "stenosis": "1-24% — minimal non-obstructive",
            "management": "Preventive therapy; risk factor modification",
            "severity": SeverityLevel.LOW,
        },
        "2": {
            "stenosis": "25-49% — mild non-obstructive",
            "management": "Preventive therapy; consider functional testing if symptomatic",
            "severity": SeverityLevel.LOW,
        },
        "3": {
            "stenosis": "50-69% — moderate stenosis",
            "management": "Consider functional testing; anti-ischaemic therapy",
            "severity": SeverityLevel.MODERATE,
        },
        "4A": {
            "stenosis": "70-99% — severe stenosis (one- or two-vessel)",
            "management": "Functional testing or invasive angiography; consider revascularisation",
            "severity": SeverityLevel.HIGH,
        },
        "4B": {
            "stenosis": "70-99% — severe stenosis (left main or three-vessel)",
            "management": "Invasive angiography with intent to revascularise",
            "severity": SeverityLevel.CRITICAL,
        },
        "5": {
            "stenosis": "100% — total occlusion",
            "management": "Invasive angiography; consider CTO PCI if clinically indicated",
            "severity": SeverityLevel.CRITICAL,
        },
    }

    _HIGH_RISK_PLAQUE = {
        "low_attenuation": (
            "Low-attenuation plaque (<30 HU) — associated with vulnerable "
            "plaque and future ACS"
        ),
        "positive_remodeling": (
            "Positive remodelling (RI >1.1) — marker of plaque vulnerability"
        ),
        "napkin_ring_sign": (
            "Napkin-ring sign — strong predictor of plaque rupture"
        ),
        "spotty_calcification": (
            "Spotty calcification — associated with plaque instability"
        ),
    }

    _VALID_CADRADS = {"0", "1", "2", "3", "4A", "4B", "5"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("calcium_score", 0)
        inp.setdefault("cad_rads", "0")
        inp.setdefault("plaque_features", [])
        inp.setdefault("symptoms", [])
        inp.setdefault("risk_factors", [])
        # Validate calcium_score
        try:
            cs = float(inp["calcium_score"])
            if cs < 0:
                warnings.append(
                    f"Calcium score {inp['calcium_score']} is negative — must be >= 0"
                )
                inp["calcium_score"] = 0
        except (TypeError, ValueError):
            warnings.append(
                f"Calcium score '{inp['calcium_score']}' is not numeric — defaulting to 0"
            )
            inp["calcium_score"] = 0
        # Normalise CAD-RADS string
        inp["cad_rads"] = str(inp["cad_rads"]).upper().replace("CAD-RADS ", "")
        # Validate cad_rads
        if inp["cad_rads"] not in self._VALID_CADRADS:
            warnings.append(
                f"CAD-RADS grade '{inp['cad_rads']}' is not valid "
                f"(expected one of {', '.join(sorted(self._VALID_CADRADS))}) — defaulting to 0"
            )
            inp["cad_rads"] = "0"
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        # --- Calcium score risk stratification ---
        cac = int(inputs["calcium_score"])
        if cac == 0:
            cac_sev = SeverityLevel.LOW
            findings.append(
                f"Calcium score {cac}: very low risk (zero plaque burden)"
            )
        elif cac < 100:
            cac_sev = SeverityLevel.LOW
            findings.append(
                f"Calcium score {cac}: low risk (mild plaque burden)"
            )
        elif cac < 300:
            cac_sev = SeverityLevel.MODERATE
            findings.append(
                f"Calcium score {cac}: moderate risk (significant plaque burden)"
            )
            recommendations.append(
                "Intensify statin therapy; target LDL <70 mg/dL"
            )
        elif cac < 1000:
            cac_sev = SeverityLevel.HIGH
            findings.append(
                f"Calcium score {cac}: high risk (extensive plaque burden)"
            )
            recommendations.append(
                "High-intensity statin; consider PCSK9i if LDL not at goal"
            )
            recommendations.append(
                "Consider functional stress testing to assess ischaemia"
            )
        else:
            cac_sev = SeverityLevel.CRITICAL
            findings.append(
                f"Calcium score {cac}: very high risk (>=1000; severe, "
                "diffuse atherosclerotic plaque burden)"
            )
            recommendations.append(
                "High-intensity statin; strongly consider PCSK9i. "
                "Functional stress testing recommended."
            )
            recommendations.append(
                "Aggressive risk factor modification; LDL target <55 mg/dL"
            )

        # --- CAD-RADS grading ---
        cr = inputs["cad_rads"]
        cr_info = self._CADRADS_MGMT.get(cr, self._CADRADS_MGMT["0"])
        findings.append(f"CAD-RADS {cr}: {cr_info['stenosis']}")
        recommendations.append(
            f"CAD-RADS {cr} management: {cr_info['management']}"
        )

        # --- High-risk plaque features ---
        plaque = [
            p.lower().replace(" ", "_").replace("-", "_")
            for p in inputs["plaque_features"]
        ]
        plaque_findings: List[str] = []
        for feat in plaque:
            desc = self._HIGH_RISK_PLAQUE.get(feat)
            if desc:
                plaque_findings.append(desc)
        if plaque_findings:
            findings.append(
                f"High-risk plaque features identified ({len(plaque_findings)}): "
                + "; ".join(plaque_findings)
            )
            recommendations.append(
                "Aggressive medical therapy warranted given high-risk plaque morphology"
            )

        # --- Overall severity ---
        severity = _max_severity(
            cac_sev,
            cr_info["severity"],
            SeverityLevel.HIGH if plaque_findings else SeverityLevel.INFORMATIONAL,
        )

        # --- Cross-modal genomics trigger ---
        cr_numeric = {
            "0": 0, "1": 1, "2": 2, "3": 3, "4A": 4, "4B": 4, "5": 5,
        }.get(cr, 0)
        if cac >= 400 or cr_numeric >= 3:
            genes = ["LDLR", "PCSK9", "APOB"]
            triggers.append(
                _trigger_string(
                    "familial_hypercholesterolemia_screen",
                    genes,
                    f"CAC score {cac} and/or CAD-RADS {cr} indicate significant "
                    "atherosclerotic burden — FH genetic screening recommended",
                )
            )
            recommendations.append(
                "Recommend genetic testing for familial hypercholesterolemia "
                "(LDLR, PCSK9, APOB)"
            )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "SCCT 2022 CAD-RADS 2.0 (Cury et al., Radiology: Cardiothoracic Imaging 2022)",
                "2019 ACC/AHA Guideline on Primary Prevention (Arnett et al., Circulation 2019)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 2 — Heart Failure
# ═══════════════════════════════════════════════════════════════════════════


class HeartFailureWorkflow(BaseCardioWorkflow):
    """Heart failure classification, GDMT assessment, and device
    eligibility evaluation.

    Inputs
    ------
    lvef : float            — Left ventricular ejection fraction (%).
    nyha_class : str        — NYHA functional class (I-IV or 1-4).
    bnp : float             — BNP (pg/mL).
    nt_probnp : float       — NT-proBNP (pg/mL).
    current_meds : list     — Current medication names.
    qrs_duration : int      — QRS duration (ms).
    qrs_morphology : str    — e.g. 'normal', 'lbbb', 'rbbb'.
    age : int               — Patient age.
    race : str              — Patient race.
    """

    workflow_type = CardioWorkflowType.HEART_FAILURE

    _CARDIOMYOPATHY_GENES = [
        "TTN", "LMNA", "MYH7", "MYBPC3", "TNNT2", "TNNI3",
        "TPM1", "ACTC1", "SCN5A", "PLN", "DES", "FLNC",
        "RBM20", "BAG3", "DSP", "PKP2", "DSG2", "DSC2",
    ]

    _VALID_NYHA = {"I", "II", "III", "IV", "1", "2", "3", "4"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("lvef", 55.0)
        inp.setdefault("nyha_class", "I")
        inp.setdefault("bnp", 0.0)
        inp.setdefault("nt_probnp", 0.0)
        inp.setdefault("current_meds", [])
        inp.setdefault("qrs_duration", 100)
        inp.setdefault("qrs_morphology", "normal")
        inp.setdefault("age", 65)
        inp.setdefault("race", "other")
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(
                        f"LVEF {lv}% outside valid range (5-90%)"
                    )
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 55.0"
                )
                inp["lvef"] = 55.0
        else:
            warnings.append("Missing LVEF — results may be incomplete")
            inp["lvef"] = 55.0
        # Validate NYHA class
        nyha_str = str(inp["nyha_class"]).strip().upper()
        if nyha_str not in self._VALID_NYHA:
            warnings.append(
                f"NYHA class '{inp['nyha_class']}' is not valid "
                f"(expected I-IV or 1-4) — defaulting to I"
            )
            inp["nyha_class"] = "I"
        return inp

    # ── helpers ───────────────────────────────────────────────────────
    @staticmethod
    def _classify_ef(lvef: float) -> EjectionFractionCategory:
        if lvef <= 40:
            return EjectionFractionCategory.HFrEF
        if lvef <= 49:
            return EjectionFractionCategory.HFmrEF
        return EjectionFractionCategory.HFpEF

    @staticmethod
    def _classify_nyha(nyha_str: str) -> HeartFailureClass:
        mapping = {
            "I": HeartFailureClass.NYHA_I,
            "II": HeartFailureClass.NYHA_II,
            "III": HeartFailureClass.NYHA_III,
            "IV": HeartFailureClass.NYHA_IV,
            "1": HeartFailureClass.NYHA_I,
            "2": HeartFailureClass.NYHA_II,
            "3": HeartFailureClass.NYHA_III,
            "4": HeartFailureClass.NYHA_IV,
        }
        return mapping.get(str(nyha_str).strip().upper(), HeartFailureClass.NYHA_II)

    @staticmethod
    def _classify_stage(
        ef_cat: EjectionFractionCategory,
        nyha: HeartFailureClass,
        bnp: float,
        nt_probnp: float,
    ) -> HeartFailureStage:
        symptomatic = nyha in (HeartFailureClass.NYHA_III, HeartFailureClass.NYHA_IV)
        has_structural = (
            ef_cat != EjectionFractionCategory.HFpEF
            or bnp > 100
            or nt_probnp > 300
        )
        if symptomatic and (bnp > 600 or nt_probnp > 5000):
            return HeartFailureStage.STAGE_D
        if symptomatic or (has_structural and nyha != HeartFailureClass.NYHA_I):
            return HeartFailureStage.STAGE_C
        if has_structural:
            return HeartFailureStage.STAGE_B
        return HeartFailureStage.STAGE_A

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        lvef = float(inputs["lvef"])
        nyha = self._classify_nyha(inputs["nyha_class"])
        ef_cat = self._classify_ef(lvef)
        bnp = float(inputs["bnp"])
        nt_probnp = float(inputs["nt_probnp"])
        stage = self._classify_stage(ef_cat, nyha, bnp, nt_probnp)
        current_meds = [m.lower() for m in inputs["current_meds"]]
        qrs_dur = int(inputs["qrs_duration"])
        qrs_morph = str(inputs["qrs_morphology"]).lower()
        age = int(inputs["age"])

        findings.append(f"LVEF {lvef}% — classified as {ef_cat.value}")
        findings.append(f"NYHA class {nyha.value}")
        findings.append(f"ACC/AHA stage {stage.value}")

        if bnp > 0:
            findings.append(
                f"BNP {bnp} pg/mL" + (" (elevated)" if bnp > 100 else "")
            )
        if nt_probnp > 0:
            findings.append(
                f"NT-proBNP {nt_probnp} pg/mL"
                + (" (elevated)" if nt_probnp > 300 else "")
            )

        # --- GDMT assessment (simplified GDMTOptimizer delegate) ------
        gdmt_pillars = {
            "beta_blocker": [
                "carvedilol", "metoprolol succinate", "bisoprolol",
            ],
            "arni_acei_arb": [
                "sacubitril/valsartan", "sacubitril-valsartan",
                "enalapril", "lisinopril", "ramipril",
                "losartan", "valsartan", "candesartan",
            ],
            "mra": ["spironolactone", "eplerenone"],
            "sglt2i": ["dapagliflozin", "empagliflozin", "sotagliflozin"],
        }

        if ef_cat == EjectionFractionCategory.HFrEF:
            findings.append("HFrEF — eligible for all four GDMT pillars")
            gaps: List[str] = []
            for pillar_name, drug_list in gdmt_pillars.items():
                on_pillar = any(d in current_meds for d in drug_list)
                if not on_pillar:
                    gaps.append(pillar_name.replace("_", " ").title())
            if gaps:
                recommendations.append(
                    f"GDMT gaps identified: {', '.join(gaps)} — "
                    "initiate missing therapies"
                )
            else:
                findings.append(
                    "All four GDMT pillars are represented in current medications"
                )
                recommendations.append(
                    "Ensure all GDMT medications are up-titrated to target doses"
                )
        elif ef_cat == EjectionFractionCategory.HFmrEF:
            findings.append(
                "HFmrEF — consider GDMT (SGLT2i has Class I evidence)"
            )
            recommendations.append(
                "Initiate SGLT2 inhibitor (Class I, Level A)"
            )
            recommendations.append(
                "Consider ARNI, beta-blocker, and MRA (Class IIb)"
            )
        else:
            findings.append("HFpEF — limited disease-modifying therapies")
            recommendations.append(
                "Initiate SGLT2 inhibitor (Class IIa, Level B-R)"
            )
            recommendations.append("Diuretics for volume management")
            recommendations.append(
                "Treat comorbidities: hypertension, AF, obesity, CAD"
            )

        # --- Device eligibility ----------------------------------------
        if ef_cat == EjectionFractionCategory.HFrEF and lvef <= 35:
            # ICD
            if nyha in (HeartFailureClass.NYHA_II, HeartFailureClass.NYHA_III):
                recommendations.append(
                    "ICD recommended for primary prevention (LVEF <=35%, "
                    "NYHA II-III, on optimal GDMT >= 3 months; Class I, Level A)"
                )
            # CRT
            if qrs_dur >= 150 and "lbbb" in qrs_morph:
                recommendations.append(
                    "CRT strongly recommended (LVEF <=35%, QRS >=150 ms "
                    "with LBBB; Class I, Level A)"
                )
            elif qrs_dur >= 130 and "lbbb" in qrs_morph:
                recommendations.append(
                    "CRT recommended (LVEF <=35%, QRS 130-149 ms with LBBB; "
                    "Class IIa, Level B)"
                )
            elif qrs_dur >= 150 and "lbbb" not in qrs_morph:
                recommendations.append(
                    "CRT may be considered (LVEF <=35%, QRS >=150 ms, "
                    "non-LBBB; Class IIa, Level B)"
                )
            elif qrs_dur >= 130 and "lbbb" not in qrs_morph:
                recommendations.append(
                    "CRT may be considered (LVEF <=35%, QRS 130-149 ms, "
                    "non-LBBB; Class IIb, Level B)"
                )

        # --- Cross-modal genomics trigger ------------------------------
        if age < 50 and ef_cat == EjectionFractionCategory.HFrEF:
            triggers.append(
                _trigger_string(
                    "cardiomyopathy_gene_panel",
                    self._CARDIOMYOPATHY_GENES,
                    f"Age {age} with unexplained dilated cardiomyopathy "
                    f"(LVEF {lvef}%) — genetic testing recommended per "
                    "2022 ACC/AHA HF guidelines",
                )
            )
            recommendations.append(
                "Recommend cardiomyopathy gene panel (TTN, LMNA, MYH7, "
                "MYBPC3, etc.) and referral to cardiovascular genetics"
            )

        # --- Severity --------------------------------------------------
        if stage == HeartFailureStage.STAGE_D:
            severity = SeverityLevel.CRITICAL
        elif stage == HeartFailureStage.STAGE_C and nyha in (
            HeartFailureClass.NYHA_III,
            HeartFailureClass.NYHA_IV,
        ):
            severity = SeverityLevel.HIGH
        elif stage == HeartFailureStage.STAGE_C:
            severity = SeverityLevel.MODERATE
        elif stage == HeartFailureStage.STAGE_B:
            severity = SeverityLevel.LOW
        else:
            severity = SeverityLevel.INFORMATIONAL

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2022 AHA/ACC/HFSA Guideline for Management of Heart Failure "
                "(Heidenreich et al., Circulation 2022)",
                "2021 ESC Guidelines for Diagnosis and Treatment of HF "
                "(McDonagh et al., Eur Heart J 2021)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 3 — Valvular Disease
# ═══════════════════════════════════════════════════════════════════════════


class ValvularDiseaseWorkflow(BaseCardioWorkflow):
    """Valvular heart disease severity grading and intervention criteria
    per ASE/ACC/AHA guidelines.

    Inputs
    ------
    valve : str             — Valve name (aortic, mitral, tricuspid, pulmonic).
    pathology : str         — Lesion type (stenosis, regurgitation).
    measurements : dict     — Haemodynamic values (peak_velocity, mean_gradient,
                              valve_area, dimensionless_index, regurgitant_volume,
                              ero, vena_contracta, etc.).
    symptoms : list[str]    — Current symptoms.
    lvef : float            — LVEF (%).
    sts_score : float|None  — STS predicted risk of mortality (%).
    """

    workflow_type = CardioWorkflowType.VALVULAR_DISEASE

    # ASE criteria tables
    _AS_CRITERIA = {
        ValveSeverity.MILD: {
            "vmax": (2.0, 2.9), "mean_grad": (0, 19), "ava": (1.5, 999),
        },
        ValveSeverity.MODERATE: {
            "vmax": (3.0, 3.9), "mean_grad": (20, 39), "ava": (1.0, 1.5),
        },
        ValveSeverity.SEVERE: {
            "vmax": (4.0, 999), "mean_grad": (40, 999), "ava": (0, 0.99),
        },
    }

    _MR_CRITERIA = {
        ValveSeverity.MILD: {
            "regurg_vol": (0, 29), "ero": (0, 0.19),
            "vena_contracta": (0, 0.29),
        },
        ValveSeverity.MODERATE: {
            "regurg_vol": (30, 59), "ero": (0.20, 0.39),
            "vena_contracta": (0.30, 0.69),
        },
        ValveSeverity.SEVERE: {
            "regurg_vol": (60, 999), "ero": (0.40, 999),
            "vena_contracta": (0.70, 999),
        },
    }

    _VALID_VALVES = {"aortic", "mitral", "tricuspid", "pulmonic"}
    _VALID_PATHOLOGIES = {"stenosis", "regurgitation"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("valve", "aortic")
        inp.setdefault("pathology", "stenosis")
        inp.setdefault("measurements", {})
        inp.setdefault("symptoms", [])
        inp.setdefault("lvef", 60.0)
        inp.setdefault("sts_score", None)
        # Validate valve
        valve_lower = str(inp["valve"]).lower().strip()
        if valve_lower not in self._VALID_VALVES:
            warnings.append(
                f"Valve '{inp['valve']}' is not valid "
                f"(expected one of {', '.join(sorted(self._VALID_VALVES))}) — defaulting to aortic"
            )
            inp["valve"] = "aortic"
        # Validate pathology
        path_lower = str(inp["pathology"]).lower().strip()
        if path_lower not in self._VALID_PATHOLOGIES:
            warnings.append(
                f"Pathology '{inp['pathology']}' is not valid "
                f"(expected stenosis or regurgitation) — defaulting to stenosis"
            )
            inp["pathology"] = "stenosis"
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(f"LVEF {lv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 60.0"
                )
                inp["lvef"] = 60.0
        return inp

    # ── severity grading helpers ──────────────────────────────────────
    def _grade_aortic_stenosis(self, m: dict) -> tuple:
        vmax = m.get("peak_velocity", 0.0)
        mg = m.get("mean_gradient", 0.0)
        ava = m.get("valve_area", 999.0)
        di = m.get("dimensionless_index", None)

        findings: List[str] = []
        severity = ValveSeverity.MILD

        for sev, crit in self._AS_CRITERIA.items():
            vmax_match = crit["vmax"][0] <= vmax <= crit["vmax"][1]
            mg_match = crit["mean_grad"][0] <= mg <= crit["mean_grad"][1]
            ava_match = crit["ava"][0] <= ava <= crit["ava"][1]
            if sum([vmax_match, mg_match, ava_match]) >= 2:
                severity = sev

        findings.append(f"Peak velocity: {vmax} m/s")
        findings.append(f"Mean gradient: {mg} mmHg")
        findings.append(f"Aortic valve area: {ava} cm2")
        if di is not None:
            findings.append(f"Dimensionless index: {di}")
            if di < 0.25:
                findings.append("DI <0.25 supports severe AS")
        findings.append(f"Aortic stenosis severity: {severity.value}")
        return severity, findings

    def _grade_mitral_regurgitation(self, m: dict) -> tuple:
        rv = m.get("regurgitant_volume", 0.0)
        ero = m.get("ero", 0.0)
        vc = m.get("vena_contracta", 0.0)

        findings: List[str] = []
        severity = ValveSeverity.MILD

        for sev, crit in self._MR_CRITERIA.items():
            rv_match = crit["regurg_vol"][0] <= rv <= crit["regurg_vol"][1]
            ero_match = crit["ero"][0] <= ero <= crit["ero"][1]
            vc_match = crit["vena_contracta"][0] <= vc <= crit["vena_contracta"][1]
            if sum([rv_match, ero_match, vc_match]) >= 2:
                severity = sev

        findings.append(f"Regurgitant volume: {rv} mL")
        findings.append(f"ERO: {ero} cm2")
        findings.append(f"Vena contracta: {vc} cm")
        findings.append(f"Mitral regurgitation severity: {severity.value}")
        return severity, findings

    def _grade_generic_valve(
        self, valve: str, pathology: str, m: dict
    ) -> tuple:
        findings = [
            f"{valve.title()} {pathology}: quantitative grading "
            "deferred to specialist review"
        ]
        mg = m.get("mean_gradient", 0.0)
        va = m.get("valve_area", 999.0)
        rv = m.get("regurgitant_volume", 0.0)
        if mg > 0:
            findings.append(f"Mean gradient: {mg} mmHg")
        if va < 999:
            findings.append(f"Valve area: {va} cm2")
        if rv > 0:
            findings.append(f"Regurgitant volume: {rv} mL")
        return ValveSeverity.MODERATE, findings

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        valve = str(inputs["valve"]).lower()
        pathology = str(inputs["pathology"]).lower()
        measurements = inputs["measurements"]
        symptoms = inputs.get("symptoms", [])
        lvef = float(inputs.get("lvef", 60.0))
        sts_score = inputs.get("sts_score")

        # --- Severity grading -----------------------------------------
        if valve == "aortic" and pathology == "stenosis":
            severity, grade_findings = self._grade_aortic_stenosis(measurements)
        elif valve == "mitral" and pathology == "regurgitation":
            severity, grade_findings = self._grade_mitral_regurgitation(
                measurements
            )
        else:
            severity, grade_findings = self._grade_generic_valve(
                valve, pathology, measurements
            )
        findings.extend(grade_findings)

        # --- Intervention criteria (aortic stenosis) -------------------
        if valve == "aortic" and pathology == "stenosis":
            symptomatic = len(symptoms) > 0
            if severity == ValveSeverity.SEVERE and symptomatic:
                recommendations.append(
                    "Symptomatic severe aortic stenosis — AVR indicated "
                    "(Class I, Level A)"
                )
            elif severity == ValveSeverity.SEVERE and lvef < 50:
                recommendations.append(
                    "Asymptomatic severe AS with LVEF <50% — AVR indicated "
                    "(Class I, Level B)"
                )
            elif severity == ValveSeverity.SEVERE:
                recommendations.append(
                    "Asymptomatic severe AS with preserved EF — close "
                    "surveillance; consider AVR if very severe (Vmax >=5 m/s) "
                    "or rapid progression"
                )
            elif severity == ValveSeverity.MODERATE:
                recommendations.append(
                    "Moderate AS — serial echocardiography every 1-2 years"
                )

            # TAVR vs SAVR decision support
            if severity == ValveSeverity.SEVERE and (symptomatic or lvef < 50):
                if sts_score is not None:
                    if sts_score < 3:
                        recommendations.append(
                            f"STS score {sts_score}% (low risk) — SAVR or TAVR "
                            "per shared decision-making (for age >=65)"
                        )
                    elif sts_score < 8:
                        recommendations.append(
                            f"STS score {sts_score}% (intermediate risk) — "
                            "TAVR or SAVR per Heart Team discussion"
                        )
                    else:
                        recommendations.append(
                            f"STS score {sts_score}% (high/prohibitive risk) — "
                            "TAVR preferred; consider palliative care if "
                            "prohibitive"
                        )
                else:
                    recommendations.append(
                        "Calculate STS score for TAVR vs SAVR decision support"
                    )

        # --- Intervention criteria (mitral regurgitation) --------------
        elif valve == "mitral" and pathology == "regurgitation":
            symptomatic = len(symptoms) > 0
            if severity == ValveSeverity.SEVERE and symptomatic and lvef > 30:
                recommendations.append(
                    "Symptomatic severe MR — mitral valve surgery indicated "
                    "(Class I)"
                )
            elif severity == ValveSeverity.SEVERE and lvef <= 60:
                recommendations.append(
                    "Severe MR with LV dysfunction (LVEF <=60%) — surgery "
                    "indicated (Class I)"
                )
            elif severity == ValveSeverity.SEVERE:
                recommendations.append(
                    "Asymptomatic severe MR with preserved LV function — "
                    "consider surgery if high repair likelihood (Class IIa)"
                )
            if severity == ValveSeverity.MODERATE:
                recommendations.append(
                    "Moderate MR — serial echocardiography every 1-2 years; "
                    "medical therapy"
                )

        # --- Map valve severity to clinical SeverityLevel --------------
        _sev_map = {
            ValveSeverity.MILD: SeverityLevel.LOW,
            ValveSeverity.MODERATE: SeverityLevel.MODERATE,
            ValveSeverity.SEVERE: SeverityLevel.HIGH,
            ValveSeverity.CRITICAL: SeverityLevel.CRITICAL,
        }
        clinical_severity = _sev_map.get(severity, SeverityLevel.MODERATE)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=clinical_severity,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "2020 ACC/AHA Guideline for Management of Valvular Heart "
                "Disease (Otto et al., Circulation 2021)",
                "ASE Guidelines for Echocardiographic Assessment of Valvular "
                "Stenosis and Regurgitation",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 4 — Arrhythmia
# ═══════════════════════════════════════════════════════════════════════════


class ArrhythmiaWorkflow(BaseCardioWorkflow):
    """Arrhythmia evaluation including AF stroke risk, QTc assessment,
    bradycardia evaluation, and wide-complex tachycardia DDx.

    Inputs
    ------
    rhythm : str            — Rhythm description.
    heart_rate : int        — Ventricular rate (bpm).
    pr_interval : int       — PR interval (ms).
    qrs_duration : int      — QRS duration (ms).
    qtc : int               — Corrected QT interval (ms).
    findings : list[str]    — ECG findings.
    age : int               — Patient age.
    sex : str               — 'male' or 'female'.
    comorbidities : dict    — Boolean flags for CHF, hypertension,
                              diabetes, stroke_tia, vascular_disease,
                              family_hx_scd, etc.
    """

    workflow_type = CardioWorkflowType.ARRHYTHMIA

    _CHANNELOPATHY_GENES = [
        "SCN5A", "KCNQ1", "KCNH2", "KCNE1", "KCNE2", "KCNJ2",
        "CACNA1C", "ANK2", "RYR2", "CASQ2", "CALM1", "CALM2", "CALM3",
    ]

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("rhythm", "sinus")
        inp.setdefault("heart_rate", 72)
        inp.setdefault("pr_interval", 160)
        inp.setdefault("qrs_duration", 90)
        inp.setdefault("qtc", 420)
        inp.setdefault("findings", [])
        inp.setdefault("age", 65)
        inp.setdefault("sex", "male")
        inp.setdefault("comorbidities", {})
        # Validate heart_rate
        try:
            hr = int(inp["heart_rate"])
            if not (20 <= hr <= 300):
                warnings.append(
                    f"Heart rate {hr} bpm outside plausible range (20-300)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"Heart rate '{inp['heart_rate']}' is not numeric — defaulting to 72"
            )
            inp["heart_rate"] = 72
        # Validate QTc
        try:
            qtc = int(inp["qtc"])
            if not (200 <= qtc <= 700):
                warnings.append(
                    f"QTc {qtc} ms outside plausible range (200-700)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"QTc '{inp['qtc']}' is not numeric — defaulting to 420"
            )
            inp["qtc"] = 420
        # Validate QRS duration
        try:
            qrs = int(inp["qrs_duration"])
            if not (40 <= qrs <= 300):
                warnings.append(
                    f"QRS duration {qrs} ms outside plausible range (40-300)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"QRS duration '{inp['qrs_duration']}' is not numeric — defaulting to 90"
            )
            inp["qrs_duration"] = 90
        # Validate rhythm is not empty
        if not inp["rhythm"] or not str(inp["rhythm"]).strip():
            warnings.append("Missing rhythm — defaulting to sinus")
            inp["rhythm"] = "sinus"
        return inp

    # ── CHA2DS2-VASc ─────────────────────────────────────────────────
    @staticmethod
    def _cha2ds2_vasc(age: int, sex: str, comorbidities: dict) -> tuple:
        score = 0
        breakdown: Dict[str, int] = {}

        if comorbidities.get("chf", False):
            score += 1; breakdown["CHF"] = 1  # noqa: E702
        if comorbidities.get("hypertension", False):
            score += 1; breakdown["Hypertension"] = 1  # noqa: E702
        if age >= 75:
            score += 2; breakdown["Age_>=75"] = 2  # noqa: E702
        if comorbidities.get("diabetes", False):
            score += 1; breakdown["Diabetes"] = 1  # noqa: E702
        if comorbidities.get("stroke_tia", False):
            score += 2; breakdown["Stroke/TIA"] = 2  # noqa: E702
        if comorbidities.get("vascular_disease", False):
            score += 1; breakdown["Vascular_disease"] = 1  # noqa: E702
        if 65 <= age < 75:
            score += 1; breakdown["Age_65-74"] = 1  # noqa: E702
        if sex.lower() == "female":
            score += 1; breakdown["Sex_female"] = 1  # noqa: E702

        return score, breakdown

    @staticmethod
    def _anticoag_recommendation(
        score: int, sex: str
    ) -> tuple:
        effective = score
        if sex.lower() == "female" and score == 1:
            effective = 0

        if effective == 0:
            return (
                AnticoagulationRecommendation.NO_ANTICOAG,
                "CHA2DS2-VASc 0 (male) or 1 (female, sex-point only) — "
                "anticoagulation not indicated",
            )
        if effective == 1:
            return (
                AnticoagulationRecommendation.CONSIDER_ANTICOAG,
                f"CHA2DS2-VASc {score} — consider anticoagulation; "
                "shared decision-making recommended (Class IIa)",
            )
        return (
            AnticoagulationRecommendation.ANTICOAG_RECOMMENDED,
            f"CHA2DS2-VASc {score} — oral anticoagulation recommended "
            "(Class I, Level A). DOACs preferred over warfarin.",
        )

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []
        risk_scores: List[RiskScoreResult] = []

        rhythm = str(inputs["rhythm"]).lower()
        hr = int(inputs["heart_rate"])
        pr = int(inputs["pr_interval"])
        qrs = int(inputs["qrs_duration"])
        qtc = int(inputs["qtc"])
        age = int(inputs["age"])
        sex = str(inputs["sex"]).lower()
        comorbidities = inputs["comorbidities"]
        ecg_findings = inputs["findings"]

        findings.append(f"Rhythm: {rhythm}, HR {hr} bpm")

        # --- Atrial fibrillation pathway -------------------------------
        is_af = (
            "atrial fibrillation" in rhythm
            or "afib" in rhythm
            or "af" == rhythm or rhythm.startswith("af ")
        )
        if is_af:
            findings.append("Atrial fibrillation identified")
            score, breakdown = self._cha2ds2_vasc(age, sex, comorbidities)
            rec_enum, rec_text = self._anticoag_recommendation(score, sex)
            findings.append(f"CHA2DS2-VASc score: {score}")
            recommendations.append(rec_text)

            risk_scores.append(
                RiskScoreResult(
                    score_type=RiskScoreType.CHA2DS2_VASC,
                    score_value=float(score),
                    risk_category=(
                        "Low" if score <= 1
                        else ("Moderate" if score <= 3 else "High")
                    ),
                    interpretation=rec_text,
                    recommendations=[rec_text],
                    guideline_reference=(
                        "2023 ACC/AHA/ACCP/HRS Guideline for AF "
                        "(Joglar et al., Circulation 2024)"
                    ),
                )
            )

            if hr > 110:
                recommendations.append(
                    "Inadequate rate control (HR >110) — uptitrate "
                    "rate-control agents (beta-blocker or CCB)"
                )
            recommendations.append(
                "Assess candidacy for rhythm control (catheter ablation "
                "or antiarrhythmic drugs)"
            )

        # --- QTc assessment --------------------------------------------
        qtc_upper = 440 if sex == "male" else 460
        if qtc <= qtc_upper:
            findings.append(
                f"QTc {qtc} ms — normal (upper limit {qtc_upper} ms for {sex})"
            )
        elif qtc <= 500:
            findings.append(f"QTc {qtc} ms — borderline prolonged")
            recommendations.append(
                "Review QT-prolonging medications; monitor electrolytes "
                "(K+, Mg2+, Ca2+)"
            )
        else:
            findings.append(
                f"QTc {qtc} ms — significantly prolonged (>500 ms)"
            )
            recommendations.append(
                "URGENT: Discontinue all QT-prolonging medications; correct "
                "electrolytes; consider telemetry monitoring"
            )

        # QTc genetic trigger
        family_hx_scd = comorbidities.get("family_hx_scd", False)
        if qtc > 500 or family_hx_scd:
            triggers.append(
                _trigger_string(
                    "channelopathy_gene_panel",
                    self._CHANNELOPATHY_GENES,
                    f"QTc {qtc} ms"
                    + (
                        " and family history of sudden cardiac death"
                        if family_hx_scd else ""
                    )
                    + " — genetic testing for channelopathies recommended",
                )
            )
            recommendations.append(
                "Recommend channelopathy gene panel (SCN5A, KCNQ1, KCNH2, "
                "etc.) and referral to cardiac electrophysiology/genetics"
            )

        # --- Bradycardia assessment ------------------------------------
        if hr < 50:
            findings.append(f"Bradycardia (HR {hr} bpm)")
            if pr > 200:
                findings.append(
                    f"PR interval {pr} ms — first-degree AV block"
                )
            symptomatic_brady = any(
                s in str(ecg_findings).lower()
                for s in ["syncope", "presyncope", "dizziness", "fatigue"]
            )
            high_grade = any(
                f.lower() in [
                    "complete_heart_block", "third_degree_av_block",
                    "mobitz_ii", "complete heart block",
                    "third degree av block", "mobitz ii",
                    "high grade av block", "high_grade_av_block",
                ]
                for f in ecg_findings
            )
            if high_grade:
                recommendations.append(
                    "High-grade AV block — permanent pacemaker indicated "
                    "(Class I, Level B)"
                )
            elif symptomatic_brady:
                recommendations.append(
                    "Symptomatic bradycardia — pacemaker evaluation "
                    "indicated (Class I, Level C)"
                )
            else:
                recommendations.append(
                    "Asymptomatic bradycardia — monitor; withhold "
                    "rate-slowing agents if present"
                )

        # --- Wide-complex tachycardia ----------------------------------
        if qrs >= 120 and hr >= 100:
            findings.append(
                f"Wide-complex tachycardia (QRS {qrs} ms, HR {hr} bpm)"
            )
            vt_features: List[str] = []
            svt_features: List[str] = []

            if qrs >= 160:
                vt_features.append("Very wide QRS (>=160 ms)")
            if any(
                "av_dissociation" in str(f).lower()
                or "av dissociation" in str(f).lower()
                for f in ecg_findings
            ):
                vt_features.append("AV dissociation present")
            if any("fusion" in str(f).lower() for f in ecg_findings):
                vt_features.append("Fusion/capture beats present")
            if any("concordance" in str(f).lower() for f in ecg_findings):
                vt_features.append("Precordial concordance")

            if any(
                "bundle_branch" in str(f).lower()
                or "bundle branch" in str(f).lower()
                for f in ecg_findings
            ):
                svt_features.append("Bundle branch block pattern")
            if any(
                "baseline_bbb" in str(f).lower()
                or "baseline bbb" in str(f).lower()
                for f in ecg_findings
            ):
                svt_features.append("Known baseline BBB")

            if vt_features:
                findings.append(
                    f"Features favouring VT: {', '.join(vt_features)}"
                )
                recommendations.append(
                    "Wide-complex tachycardia with features of VT — "
                    "treat as VT until proven otherwise"
                )
            else:
                findings.append(
                    "No classic VT morphology criteria identified; "
                    "consider SVT with aberrancy"
                )
                recommendations.append(
                    "If haemodynamically unstable: electrical cardioversion. "
                    "If stable: adenosine trial may help differentiate SVT "
                    "vs VT."
                )

        # --- Severity --------------------------------------------------
        severity = SeverityLevel.INFORMATIONAL
        if qtc > 500 or (qrs >= 120 and hr >= 100):
            severity = SeverityLevel.HIGH
        elif is_af:
            severity = SeverityLevel.MODERATE
        elif qtc > qtc_upper or hr < 50:
            severity = SeverityLevel.LOW

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            risk_scores=risk_scores,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2023 ACC/AHA/ACCP/HRS Guideline for AF "
                "(Joglar et al., Circulation 2024)",
                "2018 ACC/AHA/HRS Guideline on Bradycardia and Cardiac "
                "Conduction Delay",
                "2017 AHA/ACC/HRS Guideline for Management of Ventricular "
                "Arrhythmias and SCD",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 5 — Cardiac MRI
# ═══════════════════════════════════════════════════════════════════════════


class CardiacMRIWorkflow(BaseCardioWorkflow):
    """Cardiac MRI interpretation including LGE pattern analysis,
    T1/T2 mapping, and ECV quantification.

    Inputs
    ------
    lvef : float              — LV ejection fraction (%).
    rvef : float              — RV ejection fraction (%).
    lv_edv : float            — LV end-diastolic volume (mL).
    lv_esv : float            — LV end-systolic volume (mL).
    lge_pattern : str         — LGE pattern (subendocardial, mid_wall,
                                epicardial, rv_insertion, diffuse, etc.).
    lge_extent_percent : float — LGE extent as percentage of LV mass.
    t1_native : int           — Native T1 value (ms).
    t2_value : int            — T2 mapping value (ms).
    ecv_percent : float       — Extracellular volume fraction (%).
    """

    workflow_type = CardioWorkflowType.CARDIAC_MRI

    _LGE_INTERPRETATION: Dict[str, dict] = {
        "subendocardial": {
            "etiology": "Ischaemic cardiomyopathy",
            "description": (
                "Subendocardial LGE follows a coronary territory distribution "
                "— consistent with prior myocardial infarction"
            ),
            "gene_panel": None,
        },
        "transmural": {
            "etiology": "Ischaemic cardiomyopathy (transmural MI)",
            "description": (
                "Transmural LGE indicates full-thickness infarction in a "
                "coronary distribution"
            ),
            "gene_panel": None,
        },
        "mid_wall": {
            "etiology": "Non-ischaemic: DCM, myocarditis, sarcoidosis",
            "description": (
                "Mid-wall LGE is characteristic of non-ischaemic "
                "cardiomyopathy — differential includes idiopathic DCM, "
                "resolved myocarditis, and cardiac sarcoidosis"
            ),
            "gene_panel": [
                "TTN", "LMNA", "MYH7", "MYBPC3", "TNNT2", "SCN5A",
                "FLNC", "RBM20", "BAG3", "DES", "PLN",
            ],
        },
        "epicardial": {
            "etiology": "Myocarditis",
            "description": (
                "Epicardial LGE pattern is characteristic of active or "
                "resolved myocarditis"
            ),
            "gene_panel": None,
        },
        "rv_insertion": {
            "etiology": "HCM or pulmonary hypertension",
            "description": (
                "RV insertion point LGE — seen in hypertrophic "
                "cardiomyopathy and pulmonary hypertension with RV "
                "pressure overload"
            ),
            "gene_panel": [
                "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1",
                "MYL2", "MYL3",
            ],
        },
        "diffuse": {
            "etiology": "Cardiac amyloidosis",
            "description": (
                "Diffuse LGE with difficulty nulling myocardium — highly "
                "suggestive of cardiac amyloidosis (AL or ATTR)"
            ),
            "gene_panel": ["TTR"],
        },
        "patchy": {
            "etiology": "Sarcoidosis, HCM, or myocarditis",
            "description": (
                "Patchy LGE has a broad differential including cardiac "
                "sarcoidosis, HCM, and myocarditis"
            ),
            "gene_panel": None,
        },
    }

    _VALID_LGE_PATTERNS = {
        "none", "subendocardial", "transmural", "mid_wall", "epicardial",
        "rv_insertion", "diffuse", "patchy",
    }

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("lvef", 55.0)
        inp.setdefault("rvef", 55.0)
        inp.setdefault("lv_edv", 150.0)
        inp.setdefault("lv_esv", 60.0)
        inp.setdefault("lge_pattern", "none")
        inp.setdefault("lge_extent_percent", 0.0)
        inp.setdefault("t1_native", 1000)
        inp.setdefault("t2_value", 45)
        inp.setdefault("ecv_percent", 27.0)
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(f"LVEF {lv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 55.0"
                )
                inp["lvef"] = 55.0
        # Validate RVEF
        if inp["rvef"] is not None:
            try:
                rv = float(inp["rvef"])
                if not (5 <= rv <= 90):
                    warnings.append(f"RVEF {rv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"RVEF '{inp['rvef']}' is not numeric — defaulting to 55.0"
                )
                inp["rvef"] = 55.0
        # Validate LGE pattern
        lge_norm = str(inp["lge_pattern"]).lower().replace(" ", "_").replace("-", "_")
        if lge_norm not in self._VALID_LGE_PATTERNS:
            warnings.append(
                f"LGE pattern '{inp['lge_pattern']}' is not a standard pattern — "
                "may require manual review"
            )
        # Validate ECV percent
        try:
            ecv = float(inp["ecv_percent"])
            if not (0 <= ecv <= 100):
                warnings.append(
                    f"ECV {ecv}% outside valid range (0-100%)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"ECV '{inp['ecv_percent']}' is not numeric — defaulting to 27.0"
            )
            inp["ecv_percent"] = 27.0
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        lvef = float(inputs["lvef"])
        rvef = float(inputs["rvef"])
        lv_edv = float(inputs["lv_edv"])
        lv_esv = float(inputs["lv_esv"])
        lge_pattern = (
            str(inputs["lge_pattern"]).lower().replace(" ", "_").replace("-", "_")
        )
        lge_extent = float(inputs["lge_extent_percent"])
        t1 = int(inputs["t1_native"])
        t2 = int(inputs["t2_value"])
        ecv = float(inputs["ecv_percent"])

        # --- Ventricular function --------------------------------------
        findings.append(f"LV ejection fraction: {lvef}%")
        findings.append(f"RV ejection fraction: {rvef}%")
        findings.append(f"LV EDV: {lv_edv} mL, LV ESV: {lv_esv} mL")

        if lvef < 50:
            findings.append("LV systolic dysfunction identified")
        if rvef < 45:
            findings.append("RV systolic dysfunction identified")

        # --- LGE interpretation ----------------------------------------
        lge_info = self._LGE_INTERPRETATION.get(lge_pattern)
        if lge_info:
            findings.append(
                f"LGE pattern: {lge_pattern.replace('_', ' ')}"
            )
            findings.append(f"LGE extent: {lge_extent}% of LV mass")
            findings.append(f"LGE interpretation: {lge_info['description']}")
            findings.append(f"Most likely etiology: {lge_info['etiology']}")

            if lge_extent > 25:
                recommendations.append(
                    f"Extensive LGE ({lge_extent}%) — associated with "
                    "increased arrhythmic risk and adverse prognosis"
                )
            if lge_extent > 15 and lvef <= 35:
                recommendations.append(
                    "Consider ICD for primary prevention given significant "
                    "scar burden with reduced LVEF"
                )

            # Cross-modal triggers
            if lge_info["gene_panel"]:
                trigger_reason_map = {
                    "mid_wall": (
                        "Non-ischaemic DCM pattern on CMR — cardiomyopathy "
                        "gene panel recommended"
                    ),
                    "rv_insertion": (
                        "HCM pattern on CMR — sarcomeric gene panel "
                        "recommended"
                    ),
                    "diffuse": (
                        "Amyloid pattern on CMR — TTR gene testing "
                        "recommended; also obtain serum/urine immunofixation "
                        "for AL amyloid"
                    ),
                }
                reason = trigger_reason_map.get(
                    lge_pattern,
                    f"{lge_pattern} LGE pattern — genetic evaluation "
                    "may be indicated",
                )
                triggers.append(
                    _trigger_string(
                        "cmr_genetic_panel",
                        lge_info["gene_panel"],
                        reason,
                    )
                )
                gene_str = ", ".join(lge_info["gene_panel"][:6])
                if len(lge_info["gene_panel"]) > 6:
                    gene_str += "..."
                recommendations.append(
                    f"Genetic testing recommended: {gene_str}"
                )
        elif lge_pattern == "none":
            findings.append("No late gadolinium enhancement identified")
        else:
            findings.append(
                f"LGE pattern: {lge_pattern} (atypical — manual review "
                "recommended)"
            )

        # --- T1 mapping ------------------------------------------------
        if t1 > 1100:
            findings.append(
                f"Native T1: {t1} ms — elevated (normal ~950-1050 ms at "
                "1.5T). Differential: myocardial fibrosis, amyloid "
                "infiltration, or oedema."
            )
            recommendations.append(
                "Correlate elevated T1 with clinical context and LGE pattern"
            )
        elif t1 < 900:
            findings.append(
                f"Native T1: {t1} ms — low. Differential: iron overload "
                "(Anderson-Fabry disease, haemochromatosis), fat "
                "infiltration."
            )
            recommendations.append(
                "Check serum ferritin and transferrin saturation; consider "
                "GLA gene testing"
            )
        else:
            findings.append(f"Native T1: {t1} ms — within normal range")

        # --- T2 mapping ------------------------------------------------
        if t2 > 55:
            findings.append(
                f"T2: {t2} ms — elevated (normal <55 ms at 1.5T). "
                "Suggests myocardial oedema/active inflammation."
            )
            recommendations.append(
                "Elevated T2 indicates active inflammation — consider "
                "endomyocardial biopsy if clinical picture supports "
                "myocarditis or sarcoidosis"
            )
        else:
            findings.append(
                f"T2: {t2} ms — within normal range (no evidence of "
                "active oedema)"
            )

        # --- ECV -------------------------------------------------------
        if ecv > 30:
            findings.append(
                f"ECV: {ecv}% — elevated (normal 25-30%). Indicates "
                "diffuse myocardial fibrosis or infiltrative process."
            )
            if ecv > 40:
                findings.append(
                    f"ECV {ecv}% is markedly elevated — strongly suggestive "
                    "of infiltrative cardiomyopathy (amyloid, sarcoid) or "
                    "severe diffuse fibrosis"
                )
                recommendations.append(
                    "Markedly elevated ECV — pursue tissue diagnosis; "
                    "consider bone scintigraphy (Tc99m-PYP) for ATTR amyloid"
                )
        else:
            findings.append(f"ECV: {ecv}% — within normal range")

        # --- Severity --------------------------------------------------
        if lvef < 35 or lge_extent > 25 or ecv > 40:
            severity = SeverityLevel.HIGH
        elif lvef < 50 or lge_extent > 10 or ecv > 30 or t2 > 55:
            severity = SeverityLevel.MODERATE
        elif lge_pattern != "none" and lge_info:
            severity = SeverityLevel.LOW
        else:
            severity = SeverityLevel.INFORMATIONAL

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "SCMR 2020 Recommendations for CMR in Clinical Practice",
                "2022 AHA Scientific Statement on CMR in Non-Ischaemic "
                "Myocardial Inflammation",
                "Society for CMR Consensus Statement on Parametric Mapping",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 6 — Stress Test
# ═══════════════════════════════════════════════════════════════════════════


class StressTestWorkflow(BaseCardioWorkflow):
    """Stress test interpretation including Duke Treadmill Score,
    nuclear perfusion analysis, and heart rate recovery.

    Inputs
    ------
    test_type : str               — 'exercise' or 'pharm'.
    exercise_time_min : float     — Exercise duration (minutes, Bruce protocol).
    max_st_deviation_mm : float   — Maximum ST segment deviation (mm).
    angina_index : int            — 0 = none, 1 = non-limiting, 2 = limiting.
    peak_hr : int                 — Peak heart rate (bpm).
    percent_mphr : float          — Percent of maximum predicted HR achieved.
    bp_response : str             — 'normal', 'hypotensive', 'hypertensive'.
    perfusion_defect : str        — 'none', 'fixed', 'reversible', 'mixed'.
    hr_recovery_1min : int|None   — HR drop in the first minute of recovery.
    """

    workflow_type = CardioWorkflowType.STRESS_TEST

    _VALID_TEST_TYPES = {"exercise", "pharm"}
    _VALID_BP_RESPONSES = {"normal", "hypotensive", "hypertensive", "drop", "decrease", "exaggerated"}
    _VALID_PERFUSION = {"none", "normal", "fixed", "reversible", "mixed"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("test_type", "exercise")
        inp.setdefault("exercise_time_min", 0.0)
        inp.setdefault("max_st_deviation_mm", 0.0)
        inp.setdefault("angina_index", 0)
        inp.setdefault("peak_hr", 150)
        inp.setdefault("percent_mphr", 85.0)
        inp.setdefault("bp_response", "normal")
        inp.setdefault("perfusion_defect", "none")
        inp.setdefault("hr_recovery_1min", None)
        inp.setdefault("age", 55)
        # Validate test_type
        tt = str(inp["test_type"]).lower().strip()
        if tt not in self._VALID_TEST_TYPES:
            warnings.append(
                f"Test type '{inp['test_type']}' is not valid "
                f"(expected exercise or pharm) — defaulting to exercise"
            )
            inp["test_type"] = "exercise"
        # Validate exercise_time_min
        try:
            et = float(inp["exercise_time_min"])
            if et < 0:
                warnings.append(
                    f"Exercise time {et} min is negative — must be >= 0"
                )
                inp["exercise_time_min"] = 0.0
        except (TypeError, ValueError):
            warnings.append(
                f"Exercise time '{inp['exercise_time_min']}' is not numeric — defaulting to 0"
            )
            inp["exercise_time_min"] = 0.0
        # Validate angina_index
        try:
            ai = int(inp["angina_index"])
            if ai not in (0, 1, 2):
                warnings.append(
                    f"Angina index {ai} is not valid (expected 0, 1, or 2) — defaulting to 0"
                )
                inp["angina_index"] = 0
        except (TypeError, ValueError):
            warnings.append(
                f"Angina index '{inp['angina_index']}' is not numeric — defaulting to 0"
            )
            inp["angina_index"] = 0
        # Validate bp_response
        bp = str(inp["bp_response"]).lower().strip()
        if bp not in self._VALID_BP_RESPONSES:
            warnings.append(
                f"BP response '{inp['bp_response']}' is not valid — defaulting to normal"
            )
            inp["bp_response"] = "normal"
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        test_type = str(inputs["test_type"]).lower()
        ex_time = float(inputs["exercise_time_min"])
        st_dev = float(inputs["max_st_deviation_mm"])
        angina = int(inputs["angina_index"])
        peak_hr = int(inputs["peak_hr"])
        pmphr = float(inputs["percent_mphr"])
        bp_resp = str(inputs["bp_response"]).lower()
        perfusion = str(inputs["perfusion_defect"]).lower()
        hr_rec = inputs.get("hr_recovery_1min")

        dts_sev = SeverityLevel.INFORMATIONAL
        findings.append(f"Stress test type: {test_type}")

        # --- Exercise parameters ---------------------------------------
        if test_type == "exercise":
            # Bruce protocol MET approximation (ACSM formula)
            mets = round(2.94 + ex_time * 1.06, 1) if ex_time <= 21 else round(23.0, 1)
            findings.append(
                f"Exercise time: {ex_time} min (estimated {mets} METs)"
            )
            if mets < 5:
                findings.append(
                    "Poor exercise capacity (<5 METs) — associated with "
                    "worse prognosis"
                )
            elif mets >= 10:
                findings.append(
                    "Excellent exercise capacity (>=10 METs) — favourable "
                    "prognosis"
                )

        # --- Duke Treadmill Score --------------------------------------
        if test_type == "exercise":
            dts = ex_time - (5.0 * st_dev) - (4.0 * angina)

            if dts >= 5:
                dts_sev = SeverityLevel.LOW
                findings.append(
                    f"Duke Treadmill Score: {dts:.1f} — Low risk "
                    "(annual mortality <1%)"
                )
            elif dts >= -10:
                dts_sev = SeverityLevel.MODERATE
                findings.append(
                    f"Duke Treadmill Score: {dts:.1f} — Moderate risk "
                    "(annual mortality 2-3%)"
                )
                recommendations.append(
                    "Moderate DTS — consider additional imaging (stress echo "
                    "or nuclear perfusion)"
                )
            else:
                dts_sev = SeverityLevel.HIGH
                findings.append(
                    f"Duke Treadmill Score: {dts:.1f} — High risk "
                    "(annual mortality >=5%)"
                )
                recommendations.append(
                    "High-risk DTS — recommend coronary angiography"
                )

        # --- ST segment analysis ---------------------------------------
        findings.append(f"Max ST deviation: {st_dev} mm")
        if st_dev >= 2.0:
            findings.append(
                "Significant ST depression (>=2 mm) — high probability "
                "of ischaemia"
            )
            recommendations.append(
                "Significant ischaemic ECG changes — refer for angiography"
            )
        elif st_dev >= 1.0:
            findings.append(
                "ST depression >=1 mm — positive for ischaemia"
            )
        else:
            findings.append("No significant ST changes")

        # --- Angina ----------------------------------------------------
        if angina == 2:
            findings.append("Exercise-limiting angina (angina index 2)")
        elif angina == 1:
            findings.append(
                "Non-limiting angina during exercise (angina index 1)"
            )
        else:
            findings.append("No angina during testing")

        # --- Chronotropic response -------------------------------------
        findings.append(f"Peak HR: {peak_hr} bpm ({pmphr}% MPHR)")
        if pmphr < 85:
            findings.append(
                "Chronotropic incompetence (<85% MPHR) — limits test "
                "sensitivity"
            )
            recommendations.append(
                "Sub-maximal heart rate response — consider pharmacological "
                "stress imaging for definitive ischaemia assessment"
            )

        # --- BP response -----------------------------------------------
        if bp_resp in ("hypotensive", "drop", "decrease"):
            findings.append(
                "Hypotensive BP response during exercise — concerning for "
                "severe CAD or LV dysfunction"
            )
            recommendations.append(
                "Abnormal BP response — evaluate LV function; consider "
                "angiography"
            )
        elif bp_resp in ("hypertensive", "exaggerated"):
            findings.append(
                "Exaggerated hypertensive BP response during exercise"
            )
        else:
            findings.append("Normal blood pressure response to exercise")

        # --- Heart rate recovery ---------------------------------------
        if hr_rec is not None:
            hr_rec = int(hr_rec)
            if hr_rec <= 12:
                findings.append(
                    f"Abnormal heart rate recovery: {hr_rec} bpm drop in "
                    "1st minute (normal >=12 bpm). Associated with increased "
                    "all-cause mortality."
                )
                recommendations.append(
                    "Abnormal HR recovery — independent predictor of "
                    "mortality; assess autonomic function and overall CV risk"
                )
            else:
                findings.append(
                    f"Normal heart rate recovery: {hr_rec} bpm drop in "
                    "1st minute"
                )

        # --- Nuclear perfusion -----------------------------------------
        if perfusion in ("normal", "none"):
            findings.append(
                "Nuclear perfusion: normal — no perfusion defects"
            )
        elif "reversible" in perfusion and "fixed" in perfusion:
            findings.append(
                "Nuclear perfusion: mixed defect — both reversible "
                "(ischaemia) and fixed (scar) components"
            )
            recommendations.append(
                "Mixed perfusion defect — coronary angiography recommended "
                "to assess viability and guide revascularisation"
            )
        elif "reversible" in perfusion:
            findings.append(
                "Nuclear perfusion: reversible defect — indicates "
                "inducible ischaemia"
            )
            recommendations.append(
                "Reversible perfusion defect (ischaemia) — consider "
                "coronary angiography"
            )
        elif "fixed" in perfusion:
            findings.append(
                "Nuclear perfusion: fixed defect — indicates myocardial scar"
            )
            recommendations.append(
                "Fixed perfusion defect (scar) — assess viability if "
                "revascularisation considered"
            )
        else:
            findings.append(f"Nuclear perfusion: {perfusion}")

        # --- Severity --------------------------------------------------
        severity = dts_sev
        if "reversible" in perfusion:
            severity = _max_severity(severity, SeverityLevel.MODERATE)
        if st_dev >= 2.0 or bp_resp in ("hypotensive", "drop"):
            severity = _max_severity(severity, SeverityLevel.HIGH)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "2021 AHA/ACC Guideline for Chest Pain Evaluation "
                "(Gulati et al., Circulation 2021)",
                "ASNC Guidelines for Nuclear Cardiology Procedures",
                "Mark DB et al. Duke Treadmill Score (NEJM 1991)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 7 — Preventive Risk
# ═══════════════════════════════════════════════════════════════════════════


class PreventiveRiskWorkflow(BaseCardioWorkflow):
    """ASCVD risk assessment, statin decision framework, and
    risk-enhancing factor evaluation.

    Inputs
    ------
    age : int                        — Patient age.
    sex : str                        — 'male' or 'female'.
    race : str                       — 'white', 'african_american', etc.
    total_cholesterol : float        — Total cholesterol (mg/dL).
    hdl : float                      — HDL cholesterol (mg/dL).
    ldl : float                      — LDL cholesterol (mg/dL).
    systolic_bp : float              — Systolic blood pressure (mmHg).
    bp_treatment : bool              — On antihypertensive therapy.
    diabetes : bool                  — Diabetes mellitus.
    smoker : bool                    — Current smoker.
    family_hx_premature_ascvd : bool — Family history of premature ASCVD.
    lpa : float|None                 — Lipoprotein(a) (mg/dL).
    hscrp : float|None               — High-sensitivity CRP (mg/L).
    cac_score : int|None             — Coronary artery calcium score.
    abi : float|None                 — Ankle-brachial index.
    """

    workflow_type = CardioWorkflowType.PREVENTIVE_RISK

    _FH_GENES = ["LDLR", "PCSK9", "APOB"]

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("age", 55)
        inp.setdefault("sex", "male")
        inp.setdefault("race", "white")
        inp.setdefault("total_cholesterol", 200)
        inp.setdefault("hdl", 50)
        inp.setdefault("ldl", 130)
        inp.setdefault("systolic_bp", 130)
        inp.setdefault("bp_treatment", False)
        inp.setdefault("diabetes", False)
        inp.setdefault("smoker", False)
        inp.setdefault("family_hx_premature_ascvd", False)
        inp.setdefault("lpa", None)
        inp.setdefault("hscrp", None)
        inp.setdefault("cac_score", None)
        inp.setdefault("abi", None)
        # Validate age
        try:
            age = int(inp["age"])
            if not (18 <= age <= 120):
                warnings.append(
                    f"Age {age} outside plausible range (18-120)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"Age '{inp['age']}' is not numeric — defaulting to 55"
            )
            inp["age"] = 55
        # Validate lipid values
        for field, label, lo, hi in [
            ("total_cholesterol", "Total cholesterol", 50, 500),
            ("hdl", "HDL", 5, 150),
            ("ldl", "LDL", 10, 400),
        ]:
            try:
                val = float(inp[field])
                if not (lo <= val <= hi):
                    warnings.append(
                        f"{label} {val} mg/dL outside plausible range ({lo}-{hi})"
                    )
            except (TypeError, ValueError):
                warnings.append(
                    f"{label} '{inp[field]}' is not numeric — using default"
                )
        # Validate systolic BP
        try:
            sbp = float(inp["systolic_bp"])
            if not (60 <= sbp <= 250):
                warnings.append(
                    f"Systolic BP {sbp} mmHg outside plausible range (60-250)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"Systolic BP '{inp['systolic_bp']}' is not numeric — defaulting to 130"
            )
            inp["systolic_bp"] = 130
        return inp

    # ── simplified ASCVD risk (Pooled Cohort Equations approx.) ──────
    @staticmethod
    def _estimate_ascvd_risk(
        age: int, sex: str, race: str, tc: float, hdl: float,
        sbp: float, bp_tx: bool, dm: bool, smoker: bool,
    ) -> float:
        """Simplified 10-year ASCVD risk estimation.

        Point-based approximation of the Pooled Cohort Equations for
        clinical decision support — not a replacement for the full
        validated calculator (see risk_calculators.py for the full
        implementation).
        """
        if age < 40 or age > 79:
            return -1.0  # outside validated range

        # Simplified point-based model calibrated to PCE outputs
        points = 0.0

        # Age contribution (exponential increase)
        if sex.lower() == "male":
            base_risk = 0.5 + (age - 40) * 0.25
        else:
            base_risk = 0.2 + (age - 40) * 0.15

        # Race adjustment
        if race.lower() in ("african_american", "black"):
            base_risk *= 1.5

        # Cholesterol ratio
        tc_hdl_ratio = tc / max(hdl, 1)
        if tc_hdl_ratio > 5:
            points += (tc_hdl_ratio - 5) * 1.5
        elif tc_hdl_ratio < 3.5:
            points -= 1.0

        # Blood pressure
        if sbp >= 180:
            points += 4.0
        elif sbp >= 160:
            points += 3.0
        elif sbp >= 140:
            points += 2.0
        elif sbp >= 130:
            points += 1.0
        elif sbp < 120:
            points -= 0.5

        if bp_tx:
            points += 1.0

        # Risk factors
        if dm:
            points += 3.0
        if smoker:
            points += 3.0

        risk = base_risk + points
        return round(max(0.0, min(risk, 100.0)), 1)

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []
        risk_scores: List[RiskScoreResult] = []

        age = int(inputs["age"])
        sex = str(inputs["sex"])
        race = str(inputs["race"])
        tc = float(inputs["total_cholesterol"])
        hdl = float(inputs["hdl"])
        ldl = float(inputs["ldl"])
        sbp = float(inputs["systolic_bp"])
        bp_tx = bool(inputs["bp_treatment"])
        dm = bool(inputs["diabetes"])
        smoker = bool(inputs["smoker"])
        fhx = bool(inputs["family_hx_premature_ascvd"])
        lpa = inputs.get("lpa")
        hscrp = inputs.get("hscrp")
        cac = inputs.get("cac_score")
        abi = inputs.get("abi")

        # --- ASCVD risk ------------------------------------------------
        ascvd_risk = self._estimate_ascvd_risk(
            age, sex, race, tc, hdl, sbp, bp_tx, dm, smoker
        )

        if ascvd_risk < 0:
            findings.append(
                f"Age {age} is outside the validated range (40-79) for PCE"
            )
            risk_cat = "Not calculable"
        else:
            if ascvd_risk < 5:
                risk_cat = "Low"
            elif ascvd_risk < 7.5:
                risk_cat = "Borderline"
            elif ascvd_risk < 20:
                risk_cat = "Intermediate"
            else:
                risk_cat = "High"

            findings.append(
                f"10-year ASCVD risk: {ascvd_risk}% — {risk_cat} risk "
                "category"
            )
            risk_scores.append(
                RiskScoreResult(
                    score_type=RiskScoreType.ASCVD,
                    score_value=ascvd_risk,
                    risk_category=risk_cat,
                    interpretation=(
                        f"Estimated 10-year ASCVD risk: {ascvd_risk}% "
                        f"({risk_cat})"
                    ),
                    recommendations=[],
                    guideline_reference=(
                        "2019 ACC/AHA Primary Prevention Guideline "
                        "(Arnett et al., Circulation 2019)"
                    ),
                )
            )

        # --- Risk-enhancing factors ------------------------------------
        enhancers: List[str] = []
        if fhx:
            enhancers.append("Family history of premature ASCVD")
        if ldl >= 160:
            enhancers.append(
                f"Persistently elevated LDL (LDL {ldl} mg/dL >= 160)"
            )
        if dm:
            enhancers.append("Diabetes mellitus")
        if lpa is not None and lpa >= 50:
            enhancers.append(f"Elevated Lp(a): {lpa} mg/dL (>=50)")
        if hscrp is not None and hscrp >= 2.0:
            enhancers.append(f"Elevated hs-CRP: {hscrp} mg/L (>=2.0)")
        if abi is not None and abi < 0.9:
            enhancers.append(
                f"Low ABI: {abi} (<0.9) — peripheral arterial disease"
            )

        if enhancers:
            findings.append(
                f"Risk-enhancing factors ({len(enhancers)}): "
                + "; ".join(enhancers)
            )

        # --- CAC-based reclassification --------------------------------
        if cac is not None:
            cac = int(cac)
            if cac == 0:
                findings.append(
                    f"CAC score {cac} — may reclassify risk downward; "
                    "consider deferring statin if no other high-risk features"
                )
            elif cac < 100:
                findings.append(f"CAC score {cac} — mild plaque burden")
            elif cac < 300:
                findings.append(
                    f"CAC score {cac} — moderate plaque burden; favours "
                    "statin therapy"
                )
            else:
                findings.append(
                    f"CAC score {cac} — extensive plaque burden "
                    "(>=75th percentile likely); statin strongly recommended"
                )

        # --- Statin decision framework ---------------------------------
        if ldl >= 190:
            recommendations.append(
                "LDL >=190 mg/dL — high-intensity statin indicated "
                "(Class I, Level A). If not at goal, add ezetimibe "
                "and/or PCSK9 inhibitor."
            )
            findings.append(
                "Very high LDL (>=190 mg/dL) — suspect familial "
                "hypercholesterolemia"
            )
        elif dm and 40 <= age <= 75:
            if ascvd_risk >= 20 or len(enhancers) >= 2:
                recommendations.append(
                    "Diabetes with multiple risk factors or high ASCVD "
                    "risk — high-intensity statin (Class I, Level A)"
                )
            else:
                recommendations.append(
                    "Diabetes, age 40-75 — moderate-intensity statin "
                    "(Class I, Level A)"
                )
        elif ascvd_risk >= 20:
            recommendations.append(
                "High ASCVD risk (>=20%) — high-intensity statin to reduce "
                "LDL by >=50% (Class I, Level A)"
            )
        elif 7.5 <= ascvd_risk < 20:
            recommendations.append(
                "Intermediate ASCVD risk (7.5-20%) — moderate-to-high "
                "intensity statin (Class I, Level A). Use risk-enhancing "
                "factors and CAC to guide shared decision-making."
            )
        elif 5.0 <= ascvd_risk < 7.5:
            if enhancers:
                recommendations.append(
                    "Borderline ASCVD risk with risk-enhancing factors — "
                    "consider moderate-intensity statin (Class IIb, "
                    "Level B-R)"
                )
            else:
                recommendations.append(
                    "Borderline risk without enhancers — lifestyle "
                    "modification; consider CAC scoring for further risk "
                    "stratification"
                )
        elif ascvd_risk >= 0:
            recommendations.append(
                "Low ASCVD risk (<5%) — emphasise lifestyle modification; "
                "statin generally not recommended"
            )

        # --- FH screening criteria -------------------------------------
        fh_suspected = False
        if ldl >= 190:
            fh_suspected = True
            findings.append(
                "LDL >=190 mg/dL — meets Dutch Lipid Clinic Network "
                "screening threshold for FH"
            )
        if fhx and ldl >= 160:
            fh_suspected = True
            findings.append(
                "Family history of premature ASCVD + LDL >=160 — clinical "
                "criteria suggestive of familial hypercholesterolemia"
            )

        if fh_suspected:
            triggers.append(
                _trigger_string(
                    "familial_hypercholesterolemia_screen",
                    self._FH_GENES,
                    f"LDL {ldl} mg/dL"
                    + (
                        " with family history of premature ASCVD"
                        if fhx else ""
                    )
                    + " — FH genetic testing recommended",
                )
            )
            recommendations.append(
                "Recommend genetic testing for FH (LDLR, PCSK9, APOB) and "
                "cascade screening of first-degree relatives"
            )

        # --- Severity --------------------------------------------------
        if ascvd_risk >= 20 or ldl >= 190:
            severity = SeverityLevel.HIGH
        elif ascvd_risk >= 7.5 or (dm and enhancers):
            severity = SeverityLevel.MODERATE
        elif ascvd_risk >= 5.0 or enhancers:
            severity = SeverityLevel.LOW
        else:
            severity = SeverityLevel.INFORMATIONAL

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            risk_scores=risk_scores,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2019 ACC/AHA Guideline on Primary Prevention of CVD "
                "(Arnett et al., Circulation 2019)",
                "2018 AHA/ACC Multi-Society Cholesterol Guideline "
                "(Grundy et al., Circulation 2019)",
                "2022 ACC Expert Consensus on the Role of Coronary "
                "Artery Calcium",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 8 — Cardio-Oncology
# ═══════════════════════════════════════════════════════════════════════════


class CardioOncologyWorkflow(BaseCardioWorkflow):
    """Cancer therapy-related cardiac dysfunction (CTRCD) monitoring,
    risk stratification, and cardioprotective management.

    Inputs
    ------
    chemotherapy_agent : str     — Agent name (anthracycline, trastuzumab,
                                   immune_checkpoint_inhibitor, vegf_inhibitor).
    cumulative_dose : float      — Cumulative dose (mg/m2).
    baseline_lvef : float        — Baseline LVEF (%).
    current_lvef : float         — Current LVEF (%).
    baseline_gls : float         — Baseline global longitudinal strain (%).
    current_gls : float          — Current GLS (%).
    troponin : float             — Troponin (ng/mL).
    nt_probnp : float            — NT-proBNP (pg/mL).
    risk_factors : list[str]     — Patient risk factors.
    """

    workflow_type = CardioWorkflowType.CARDIO_ONCOLOGY

    _AGENT_PROTOCOLS: Dict[str, dict] = {
        "anthracycline": {
            "mechanism": (
                "Dose-dependent irreversible Type I cardiotoxicity "
                "(oxidative stress, topoisomerase IIb inhibition)"
            ),
            "monitoring": (
                "Echocardiogram (LVEF + GLS) at baseline, 250 mg/m2, "
                "before each cycle after 300 mg/m2, 3/6/12 months "
                "post-treatment"
            ),
            "cumulative_threshold_mg_m2": 400,
            "risk_factors": [
                "age >65", "pre-existing CVD", "prior mediastinal RT",
                "concurrent trastuzumab", "hypertension", "diabetes",
            ],
            "cardioprotection": [
                "Dexrazoxane if cumulative dose >300 mg/m2",
                "ACE inhibitor or ARB",
                "Beta-blocker (carvedilol)",
            ],
        },
        "trastuzumab": {
            "mechanism": (
                "Dose-independent reversible Type II cardiotoxicity "
                "(HER2 pathway inhibition in cardiomyocytes)"
            ),
            "monitoring": (
                "LVEF every 3 months during treatment; GLS at baseline "
                "and every 3 months"
            ),
            "cumulative_threshold_mg_m2": None,
            "risk_factors": [
                "prior anthracycline", "age >65", "BMI >30",
                "baseline LVEF <55%", "hypertension",
            ],
            "cardioprotection": [
                "ACE inhibitor or ARB if LVEF decline",
                "Beta-blocker if tachycardia or LVEF decline",
            ],
        },
        "immune_checkpoint_inhibitor": {
            "mechanism": (
                "Immune-mediated myocarditis (T-cell infiltration); "
                "rare but potentially fulminant"
            ),
            "monitoring": (
                "Baseline ECG, troponin, BNP. Troponin with each cycle "
                "for first 3 months. ECG if symptoms."
            ),
            "cumulative_threshold_mg_m2": None,
            "risk_factors": [
                "combination ICI therapy", "pre-existing autoimmune disease",
                "prior myocarditis",
            ],
            "cardioprotection": [
                "High-dose corticosteroids for ICI myocarditis",
                "Withhold ICI if confirmed myocarditis",
            ],
        },
        "vegf_inhibitor": {
            "mechanism": (
                "Hypertension, arterial thromboembolism, heart failure "
                "(endothelial dysfunction)"
            ),
            "monitoring": (
                "BP every 2 weeks for first cycle, then every cycle. "
                "LVEF if symptoms or BP uncontrolled."
            ),
            "cumulative_threshold_mg_m2": None,
            "risk_factors": [
                "pre-existing hypertension", "prior CVD", "age >65",
                "renal impairment",
            ],
            "cardioprotection": [
                "Aggressive BP management (ACEi/ARB preferred)",
                "Low-dose aspirin if thrombotic risk",
            ],
        },
    }

    _VALID_AGENTS = {"anthracycline", "trastuzumab", "immune_checkpoint_inhibitor", "vegf_inhibitor"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("chemotherapy_agent", "anthracycline")
        inp.setdefault("cumulative_dose", 0.0)
        inp.setdefault("baseline_lvef", 60.0)
        inp.setdefault("current_lvef", 58.0)
        inp.setdefault("baseline_gls", -20.0)
        inp.setdefault("current_gls", -19.0)
        inp.setdefault("troponin", 0.0)
        inp.setdefault("nt_probnp", 0.0)
        inp.setdefault("risk_factors", [])
        # Validate chemotherapy_agent
        agent_norm = (
            str(inp["chemotherapy_agent"])
            .lower().replace(" ", "_").replace("-", "_")
        )
        if agent_norm not in self._VALID_AGENTS:
            # Check partial match (the execute method also does this)
            partial = any(k in agent_norm or agent_norm in k for k in self._VALID_AGENTS)
            if not partial:
                warnings.append(
                    f"Chemotherapy agent '{inp['chemotherapy_agent']}' not recognised — "
                    "will use anthracycline protocol as fallback"
                )
        # Validate baseline LVEF
        if inp["baseline_lvef"] is not None:
            try:
                bl = float(inp["baseline_lvef"])
                if not (5 <= bl <= 90):
                    warnings.append(
                        f"Baseline LVEF {bl}% outside valid range (5-90%)"
                    )
            except (TypeError, ValueError):
                warnings.append(
                    f"Baseline LVEF '{inp['baseline_lvef']}' is not numeric — defaulting to 60.0"
                )
                inp["baseline_lvef"] = 60.0
        else:
            warnings.append("Missing baseline LVEF — results may be incomplete")
            inp["baseline_lvef"] = 60.0
        # Validate current LVEF
        if inp["current_lvef"] is not None:
            try:
                cl = float(inp["current_lvef"])
                if not (5 <= cl <= 90):
                    warnings.append(
                        f"Current LVEF {cl}% outside valid range (5-90%)"
                    )
            except (TypeError, ValueError):
                warnings.append(
                    f"Current LVEF '{inp['current_lvef']}' is not numeric — defaulting to 58.0"
                )
                inp["current_lvef"] = 58.0
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        agent = (
            str(inputs["chemotherapy_agent"])
            .lower()
            .replace(" ", "_")
            .replace("-", "_")
        )
        cum_dose = float(inputs["cumulative_dose"])
        bl_lvef = float(inputs["baseline_lvef"])
        cur_lvef = float(inputs["current_lvef"])
        bl_gls = float(inputs["baseline_gls"])
        cur_gls = float(inputs["current_gls"])
        troponin = float(inputs["troponin"])
        nt_probnp = float(inputs["nt_probnp"])
        patient_rfs = inputs.get("risk_factors", [])

        # Look up agent protocol
        protocol = self._AGENT_PROTOCOLS.get(agent)
        if protocol is None:
            for key, val in self._AGENT_PROTOCOLS.items():
                if key in agent or agent in key:
                    protocol = val
                    agent = key
                    break
        if protocol is None:
            protocol = self._AGENT_PROTOCOLS["anthracycline"]
            findings.append(
                f"Agent '{inputs['chemotherapy_agent']}' not in database — "
                "using anthracycline protocol as default"
            )

        findings.append(
            f"Chemotherapy agent: {agent.replace('_', ' ').title()}"
        )
        findings.append(f"Mechanism: {protocol['mechanism']}")
        findings.append(f"Monitoring protocol: {protocol['monitoring']}")

        # --- GLS assessment --------------------------------------------
        if bl_gls != 0:
            gls_change_pct = abs((cur_gls - bl_gls) / bl_gls) * 100
            gls_worsened = cur_gls > bl_gls  # less negative = worse
        else:
            gls_change_pct = 0.0
            gls_worsened = False

        findings.append(
            f"Baseline GLS: {bl_gls}%, Current GLS: {cur_gls}%"
        )
        subclinical_cardiotox = False
        if gls_worsened and gls_change_pct > 15:
            findings.append(
                f"GLS relative decline: {gls_change_pct:.1f}% (>15%) — "
                "SUBCLINICAL CARDIOTOXICITY detected"
            )
            recommendations.append(
                "Subclinical cardiotoxicity by GLS criteria — initiate "
                "cardioprotective therapy (ACE inhibitor and/or beta-"
                "blocker) and consider oncology consultation regarding "
                "treatment modification"
            )
            subclinical_cardiotox = True
        elif gls_worsened and gls_change_pct > 8:
            findings.append(
                f"GLS relative decline: {gls_change_pct:.1f}% (8-15%) — "
                "borderline; repeat in 2-3 weeks"
            )
        else:
            findings.append("GLS within acceptable range")

        # --- LVEF assessment -------------------------------------------
        lvef_decline = bl_lvef - cur_lvef
        findings.append(
            f"Baseline LVEF: {bl_lvef}%, Current LVEF: {cur_lvef}%"
        )
        ctrcd = False
        if lvef_decline > 10 and cur_lvef < 50:
            findings.append(
                f"LVEF decline of {lvef_decline:.1f}% to below 50% — "
                "CANCER THERAPY-RELATED CARDIAC DYSFUNCTION (CTRCD) "
                "confirmed"
            )
            ctrcd = True
            recommendations.append(
                "CTRCD confirmed — hold cardiotoxic therapy; initiate GDMT "
                "(ACE inhibitor/ARB + beta-blocker). Reassess LVEF in "
                "2-3 weeks."
            )
        elif lvef_decline > 10:
            findings.append(
                f"LVEF decline of {lvef_decline:.1f}% but remains >=50% "
                "— significant; monitor closely"
            )
            recommendations.append(
                "Significant LVEF decline — intensify monitoring; initiate "
                "cardioprotection"
            )
        elif cur_lvef < 50:
            findings.append(
                f"Current LVEF {cur_lvef}% is below 50% — evaluate for "
                "CTRCD"
            )
            ctrcd = True
        else:
            findings.append("LVEF remains within acceptable range")

        # --- Biomarkers ------------------------------------------------
        if troponin > 0:
            if troponin > 0.04:
                findings.append(
                    f"Troponin {troponin} ng/mL — elevated (above URL)"
                )
                recommendations.append(
                    "Elevated troponin during cardiotoxic therapy — suggests "
                    "myocardial injury; correlate with imaging"
                )
            else:
                findings.append(
                    f"Troponin {troponin} ng/mL — within normal range"
                )
        if nt_probnp > 0:
            if nt_probnp > 300:
                findings.append(
                    f"NT-proBNP {nt_probnp} pg/mL — elevated "
                    "(myocardial stress)"
                )
            else:
                findings.append(
                    f"NT-proBNP {nt_probnp} pg/mL — within normal range"
                )

        # --- Cumulative dose -------------------------------------------
        threshold = protocol.get("cumulative_threshold_mg_m2")
        if threshold and cum_dose > 0:
            pct_threshold = (cum_dose / threshold) * 100
            if pct_threshold >= 100:
                findings.append(
                    f"Cumulative dose {cum_dose} mg/m2 exceeds threshold "
                    f"({threshold} mg/m2) — high cardiotoxicity risk"
                )
                recommendations.append(
                    f"Cumulative dose exceeds {threshold} mg/m2 — strongly "
                    "consider dose reduction or alternative agent"
                )
            elif pct_threshold >= 75:
                findings.append(
                    f"Cumulative dose {cum_dose} mg/m2 approaching threshold "
                    f"({pct_threshold:.0f}% of {threshold} mg/m2)"
                )
                recommendations.append(
                    "Approaching cumulative dose limit — ensure dexrazoxane "
                    "and cardioprotective agents are in place"
                )

        # --- Risk stratification ---------------------------------------
        agent_rfs = protocol.get("risk_factors", [])
        matched_rfs = [
            rf for rf in patient_rfs
            if any(r.lower() in rf.lower() for r in agent_rfs)
        ]
        n_rfs = len(matched_rfs) + len(
            [rf for rf in agent_rfs
             if any(rf.lower() in p.lower() for p in patient_rfs)]
        )

        if ctrcd or (threshold and cum_dose > threshold):
            tox_risk = CardiotoxicityRisk.VERY_HIGH
        elif subclinical_cardiotox or n_rfs >= 3:
            tox_risk = CardiotoxicityRisk.HIGH
        elif n_rfs >= 1 or (gls_worsened and gls_change_pct > 8):
            tox_risk = CardiotoxicityRisk.MODERATE
        else:
            tox_risk = CardiotoxicityRisk.LOW

        findings.append(f"Overall cardiotoxicity risk: {tox_risk.value}")

        # --- Cardioprotective recommendations --------------------------
        if tox_risk in (CardiotoxicityRisk.HIGH, CardiotoxicityRisk.VERY_HIGH):
            for cp in protocol.get("cardioprotection", []):
                recommendations.append(f"Cardioprotection: {cp}")

        # --- Severity --------------------------------------------------
        _tox_sev = {
            CardiotoxicityRisk.LOW: SeverityLevel.LOW,
            CardiotoxicityRisk.MODERATE: SeverityLevel.MODERATE,
            CardiotoxicityRisk.HIGH: SeverityLevel.HIGH,
            CardiotoxicityRisk.VERY_HIGH: SeverityLevel.CRITICAL,
        }
        severity = _tox_sev.get(tox_risk, SeverityLevel.MODERATE)

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "2022 ESC Guidelines on Cardio-Oncology "
                "(Lyon et al., Eur Heart J 2022)",
                "2023 ACC Expert Consensus on Cardiovascular Complications "
                "of Cancer Therapy",
                "ASCO Clinical Practice Guideline on Prevention of Cardiac "
                "Dysfunction (Armenian et al., JCO 2017)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 9 — Acute Decompensated Heart Failure
# ═══════════════════════════════════════════════════════════════════════════


class AcuteDecompensatedHFWorkflow(BaseCardioWorkflow):
    """Acute decompensated heart failure (ADHF) assessment with hemodynamic
    profiling, IV diuretic escalation, and inotrope/MCS decision support.

    Inputs
    ------
    systolic_bp : float         — Systolic blood pressure (mmHg).
    heart_rate : int            — Heart rate (bpm).
    creatinine : float          — Serum creatinine (mg/dL).
    potassium : float           — Serum potassium (mEq/L).
    bnp_or_ntprobnp : float     — BNP or NT-proBNP (pg/mL).
    lvef : float                — Left ventricular ejection fraction (%).
    congestion_signs : list     — Signs of congestion (e.g. 'jvd', 'rales',
                                  'peripheral_edema', 'ascites').
    perfusion_status : str      — 'warm' or 'cold'.
    congestion_status : str     — 'wet' or 'dry'.
    current_iv_diuretics : bool — Currently on IV loop diuretics.
    current_inotropes : bool    — Currently on inotropes.
    """

    workflow_type = CardioWorkflowType.ACUTE_DECOMPENSATED_HF

    _CARDIOMYOPATHY_GENES = [
        "TTN", "LMNA", "MYH7", "MYBPC3", "TNNT2", "TNNI3",
        "TPM1", "ACTC1", "SCN5A", "PLN", "DES", "FLNC",
        "RBM20", "BAG3", "DSP", "PKP2", "DSG2", "DSC2",
    ]

    _VALID_PERFUSION = {"warm", "cold"}
    _VALID_CONGESTION = {"wet", "dry"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("systolic_bp", 110.0)
        inp.setdefault("heart_rate", 90)
        inp.setdefault("creatinine", 1.0)
        inp.setdefault("potassium", 4.0)
        inp.setdefault("bnp_or_ntprobnp", 0.0)
        inp.setdefault("lvef", 35.0)
        inp.setdefault("congestion_signs", [])
        inp.setdefault("perfusion_status", "warm")
        inp.setdefault("congestion_status", "wet")
        inp.setdefault("current_iv_diuretics", False)
        inp.setdefault("current_inotropes", False)
        # Normalise
        inp["perfusion_status"] = str(inp["perfusion_status"]).lower().strip()
        inp["congestion_status"] = str(inp["congestion_status"]).lower().strip()
        # Validate perfusion_status
        if inp["perfusion_status"] not in self._VALID_PERFUSION:
            warnings.append(
                f"Perfusion status '{inp['perfusion_status']}' is not valid "
                "(expected warm or cold) — defaulting to warm"
            )
            inp["perfusion_status"] = "warm"
        # Validate congestion_status
        if inp["congestion_status"] not in self._VALID_CONGESTION:
            warnings.append(
                f"Congestion status '{inp['congestion_status']}' is not valid "
                "(expected wet or dry) — defaulting to wet"
            )
            inp["congestion_status"] = "wet"
        # Validate systolic_bp
        try:
            sbp = float(inp["systolic_bp"])
            if not (30 <= sbp <= 300):
                warnings.append(
                    f"Systolic BP {sbp} mmHg outside plausible range (30-300)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"Systolic BP '{inp['systolic_bp']}' is not numeric — defaulting to 110"
            )
            inp["systolic_bp"] = 110.0
        # Validate potassium
        try:
            k = float(inp["potassium"])
            if not (1.0 <= k <= 9.0):
                warnings.append(
                    f"Potassium {k} mEq/L outside plausible range (1.0-9.0)"
                )
        except (TypeError, ValueError):
            warnings.append(
                f"Potassium '{inp['potassium']}' is not numeric — defaulting to 4.0"
            )
            inp["potassium"] = 4.0
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(f"LVEF {lv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 35.0"
                )
                inp["lvef"] = 35.0
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        sbp = float(inputs["systolic_bp"])
        hr = int(inputs["heart_rate"])
        creatinine = float(inputs["creatinine"])
        potassium = float(inputs["potassium"])
        natriuretic = float(inputs["bnp_or_ntprobnp"])
        lvef = float(inputs["lvef"])
        congestion_signs = inputs["congestion_signs"]
        perfusion = inputs["perfusion_status"]
        congestion = inputs["congestion_status"]
        on_iv_diuretics = bool(inputs["current_iv_diuretics"])
        on_inotropes = bool(inputs["current_inotropes"])

        # --- Hemodynamic profile classification ---
        profile = f"{perfusion}-{congestion}"
        profile_labels = {
            "warm-wet": "Warm-Wet (most common ADHF profile — congestion predominant)",
            "cold-wet": "Cold-Wet (low-output with congestion — cardiogenic shock spectrum)",
            "warm-dry": "Warm-Dry (euvolemic — compensated, may not require acute intervention)",
            "cold-dry": "Cold-Dry (low-output without congestion — hypovolemic / over-diuresed)",
        }
        findings.append(
            f"Hemodynamic profile: {profile_labels.get(profile, profile)}"
        )
        findings.append(f"SBP {sbp} mmHg, HR {hr} bpm, LVEF {lvef}%")

        if congestion_signs:
            findings.append(
                f"Congestion signs: {', '.join(congestion_signs)}"
            )
        if natriuretic > 0:
            findings.append(
                f"BNP/NT-proBNP {natriuretic} pg/mL"
                + (" (markedly elevated)" if natriuretic > 1000 else "")
            )

        # --- Profile-specific management ---
        if profile == "warm-wet":
            recommendations.append(
                "IV loop diuretic escalation: furosemide 40-80 mg IV; "
                "double dose if inadequate response at 2 hours"
            )
            recommendations.append(
                "Target net negative 1-2 L/day fluid balance"
            )
            if on_iv_diuretics:
                recommendations.append(
                    "If diuretic-resistant: consider adding metolazone 2.5-5 mg PO "
                    "or chlorothiazide 500 mg IV (sequential nephron blockade)"
                )
            if sbp > 100:
                recommendations.append(
                    "SGLT2i can be initiated in-hospital if SBP >100 "
                    "(EMPULSE trial: empagliflozin safe and beneficial in ADHF)"
                )
        elif profile == "cold-wet":
            findings.append("WARNING: Cardiogenic shock spectrum identified")
            recommendations.append(
                "Consider inotrope initiation: dobutamine 2-5 mcg/kg/min "
                "or milrinone 0.125-0.375 mcg/kg/min"
            )
            if sbp < 90:
                recommendations.append(
                    "SBP <90 mmHg — evaluate for mechanical circulatory support "
                    "(Impella/IABP); consider vasopressors"
                )
            recommendations.append(
                "Cautious diuresis with hemodynamic monitoring; "
                "consider PA catheter-guided management"
            )
        elif profile == "cold-dry":
            recommendations.append(
                "Cautious volume challenge: 250 mL normal saline bolus; "
                "reassess hemodynamics before further fluid administration"
            )
            recommendations.append(
                "Re-evaluate diuretic dosing — may be over-diuresed"
            )
        else:  # warm-dry
            recommendations.append(
                "Euvolemic and well-perfused — optimise oral GDMT; "
                "may not require acute inpatient intervention"
            )

        # --- Electrolyte and renal monitoring ---
        if potassium < 3.5:
            findings.append(f"Hypokalaemia (K+ {potassium} mEq/L)")
            recommendations.append("Replete potassium to target >4.0 mEq/L")
        elif potassium > 5.5:
            findings.append(f"Hyperkalaemia (K+ {potassium} mEq/L)")
            recommendations.append(
                "Hold MRA/ACEi; consider potassium-lowering therapy"
            )
        if creatinine > 2.0:
            findings.append(f"Elevated creatinine {creatinine} mg/dL — monitor for cardiorenal syndrome")

        # --- SGLT2i in-hospital initiation ---
        if sbp > 100 and profile != "cold-wet":
            recommendations.append(
                "Initiate SGLT2 inhibitor in-hospital per EMPULSE trial data "
                "(clinical benefit with empagliflozin 10 mg daily in ADHF)"
            )

        # --- Severity ---
        if sbp < 90 or on_inotropes or profile == "cold-wet":
            severity = SeverityLevel.CRITICAL
        elif on_iv_diuretics or profile == "warm-wet":
            severity = SeverityLevel.HIGH
        else:
            severity = SeverityLevel.MODERATE

        # --- Cross-modal genomics trigger ---
        if lvef < 40:
            triggers.append(
                _trigger_string(
                    "new_onset_hf_genomic_workup",
                    self._CARDIOMYOPATHY_GENES,
                    f"New-onset HF with LVEF {lvef}% — consider genomic workup "
                    "for DCM/HCM genes if etiology unclear",
                )
            )
            recommendations.append(
                "If new-onset HF without clear etiology, recommend "
                "cardiomyopathy gene panel (TTN, LMNA, MYH7, MYBPC3, etc.)"
            )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2022 AHA/ACC/HFSA Guideline for Management of Heart Failure "
                "(Heidenreich et al., Circulation 2022)",
                "2021 ESC Guidelines for Diagnosis and Treatment of Acute and "
                "Chronic HF (McDonagh et al., Eur Heart J 2021)",
                "EMPULSE Trial (Voors et al., Nat Med 2022)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 10 — Post-MI Management
# ═══════════════════════════════════════════════════════════════════════════


class PostMIWorkflow(BaseCardioWorkflow):
    """Post-myocardial infarction assessment covering reperfusion evaluation,
    DAPT strategy, secondary prevention, and device eligibility.

    Inputs
    ------
    mi_type : str               — 'STEMI' or 'NSTEMI'.
    time_from_onset_hours : float — Hours from symptom onset.
    reperfusion_strategy : str  — 'PCI', 'thrombolysis', or 'none'.
    lvef : float                — LVEF (%).
    culprit_vessel : str        — Culprit vessel (e.g. 'LAD', 'RCA', 'LCx').
    troponin_peak : float       — Peak troponin (ng/mL).
    complications : list        — Post-MI complications (e.g. 'cardiogenic_shock',
                                  'vsd', 'papillary_rupture', 'lv_thrombus',
                                  'pericarditis').
    current_medications : list  — Current medication names.
    """

    workflow_type = CardioWorkflowType.POST_MI

    _VALID_MI_TYPES = {"STEMI", "NSTEMI"}
    _VALID_REPERFUSION = {"pci", "thrombolysis", "none"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("mi_type", "NSTEMI")
        inp.setdefault("time_from_onset_hours", 6.0)
        inp.setdefault("reperfusion_strategy", "PCI")
        inp.setdefault("lvef", 50.0)
        inp.setdefault("culprit_vessel", "LAD")
        inp.setdefault("troponin_peak", 0.0)
        inp.setdefault("complications", [])
        inp.setdefault("current_medications", [])
        inp["mi_type"] = str(inp["mi_type"]).upper().strip()
        inp["reperfusion_strategy"] = str(inp["reperfusion_strategy"]).lower().strip()
        inp["culprit_vessel"] = str(inp["culprit_vessel"]).upper().strip()
        # Validate MI type
        if inp["mi_type"] not in self._VALID_MI_TYPES:
            warnings.append(
                f"MI type '{inp['mi_type']}' is not valid "
                "(expected STEMI or NSTEMI) — defaulting to NSTEMI"
            )
            inp["mi_type"] = "NSTEMI"
        # Validate reperfusion strategy
        if inp["reperfusion_strategy"] not in self._VALID_REPERFUSION:
            warnings.append(
                f"Reperfusion strategy '{inp['reperfusion_strategy']}' is not valid "
                "(expected PCI, thrombolysis, or none) — defaulting to PCI"
            )
            inp["reperfusion_strategy"] = "pci"
        # Validate time_from_onset_hours
        try:
            t = float(inp["time_from_onset_hours"])
            if t < 0:
                warnings.append(
                    f"Time from onset {t} hours is negative — must be >= 0"
                )
                inp["time_from_onset_hours"] = 0.0
        except (TypeError, ValueError):
            warnings.append(
                f"Time from onset '{inp['time_from_onset_hours']}' is not numeric — defaulting to 6.0"
            )
            inp["time_from_onset_hours"] = 6.0
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(f"LVEF {lv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 50.0"
                )
                inp["lvef"] = 50.0
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        mi_type = inputs["mi_type"]
        time_hrs = float(inputs["time_from_onset_hours"])
        reperfusion = inputs["reperfusion_strategy"]
        lvef = float(inputs["lvef"])
        culprit = inputs["culprit_vessel"]
        troponin_peak = float(inputs["troponin_peak"])
        complications = [c.lower() for c in inputs["complications"]]
        current_meds = [m.lower() for m in inputs["current_medications"]]

        findings.append(f"{mi_type} — culprit vessel: {culprit}")
        findings.append(f"Time from symptom onset: {time_hrs} hours")
        findings.append(f"LVEF: {lvef}%")
        if troponin_peak > 0:
            findings.append(f"Peak troponin: {troponin_peak} ng/mL")

        # --- Reperfusion assessment ---
        if mi_type == "STEMI":
            if reperfusion == "pci":
                if time_hrs <= 1.5:
                    findings.append("Door-to-balloon time within 90-minute target")
                else:
                    findings.append(
                        f"Door-to-balloon time {time_hrs}h — exceeds 90-minute target"
                    )
            elif reperfusion == "thrombolysis":
                if time_hrs <= 0.5:
                    findings.append("Door-to-needle time within 30-minute target")
                else:
                    findings.append(
                        f"Door-to-needle time {time_hrs}h — exceeds 30-minute target"
                    )
            else:
                findings.append(
                    "No reperfusion strategy employed — evaluate for rescue PCI"
                )
                recommendations.append(
                    "Urgent consideration for rescue PCI if within window"
                )
        else:
            findings.append(f"NSTEMI — reperfusion strategy: {reperfusion}")

        # --- DAPT ---
        has_aspirin = any("aspirin" in m or "asa" in m for m in current_meds)
        has_p2y12 = any(
            d in " ".join(current_meds)
            for d in ["ticagrelor", "clopidogrel", "prasugrel"]
        )
        recommendations.append(
            "DAPT: aspirin 81 mg daily + P2Y12 inhibitor "
            "(ticagrelor 90 mg BID preferred over clopidogrel per PLATO trial)"
        )
        recommendations.append(
            "DAPT duration: 12 months standard; consider 6 months if high "
            "bleeding risk; extend 18-36 months if high ischemic/low bleeding risk"
        )
        if not has_aspirin:
            recommendations.append("Initiate aspirin 81 mg daily")
        if not has_p2y12:
            recommendations.append("Initiate ticagrelor 90 mg BID (preferred ACS P2Y12)")

        # --- Beta-blocker ---
        any(
            d in " ".join(current_meds)
            for d in ["metoprolol", "carvedilol", "bisoprolol", "atenolol"]
        )
        if "cardiogenic_shock" not in complications:
            recommendations.append(
                "Beta-blocker within 24 hours if no cardiogenic shock "
                "(metoprolol succinate or carvedilol)"
            )
        else:
            findings.append("Cardiogenic shock present — defer beta-blocker initiation")

        # --- High-intensity statin ---
        recommendations.append(
            "High-intensity statin: atorvastatin 80 mg or rosuvastatin 40 mg "
            "(Class I, Level A)"
        )

        # --- ACEi/ARB ---
        if culprit == "LAD" or lvef < 40:
            recommendations.append(
                "ACEi/ARB within 24 hours (anterior MI or LVEF <40%)"
            )

        # --- MRA ---
        if lvef <= 40 and (
            "heart_failure" in complications
            or "diabetes" in complications
            or any("diabet" in c for c in complications)
        ):
            recommendations.append(
                "Initiate MRA (eplerenone 25-50 mg daily) — LVEF <=40% "
                "with HF symptoms or diabetes (EPHESUS trial)"
            )

        # --- Cardiac rehab ---
        recommendations.append(
            "Cardiac rehabilitation referral before discharge "
            "(Class I recommendation)"
        )

        # --- ICD assessment ---
        if lvef <= 35:
            recommendations.append(
                "Reassess LVEF at 40 days post-MI; if LVEF remains <=35% "
                "on optimal GDMT, ICD indicated for primary prevention "
                "(wait for GDMT optimisation before implant)"
            )
        elif lvef <= 40:
            recommendations.append(
                "LVEF mildly reduced — reassess at 40 days post-MI after "
                "GDMT optimisation to determine ICD eligibility"
            )

        # --- Complications ---
        if complications:
            findings.append(f"Post-MI complications: {', '.join(complications)}")

        # --- Severity ---
        has_ongoing_ischemia = "ongoing_ischemia" in complications
        has_shock = "cardiogenic_shock" in complications
        if has_ongoing_ischemia or has_shock:
            severity = SeverityLevel.CRITICAL
        elif lvef < 40 or len(complications) > 0:
            severity = SeverityLevel.HIGH
        else:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2023 ACC/AHA Guideline for Management of Patients with "
                "Chronic Coronary Disease (Virani et al., Circulation 2023)",
                "2021 ACC/AHA/SCAI Guideline for Coronary Artery "
                "Revascularisation (Lawton et al., Circulation 2022)",
                "2021 ACC/AHA Chest Pain Guideline "
                "(Gulati et al., Circulation 2021)",
                "PLATO Trial (Wallentin et al., NEJM 2009)",
                "EPHESUS Trial (Pitt et al., NEJM 2003)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 11 — Myocarditis / Pericarditis
# ═══════════════════════════════════════════════════════════════════════════


class MyocarditisPericarditisWorkflow(BaseCardioWorkflow):
    """Myocarditis and pericarditis diagnosis, CMR-based assessment, and
    evidence-based treatment including return-to-play guidance.

    Inputs
    ------
    presentation : str          — 'myocarditis', 'pericarditis', or 'myopericarditis'.
    symptom_onset_days : int    — Days since symptom onset.
    troponin : float            — Troponin level (ng/mL).
    crp : float                 — C-reactive protein (mg/L).
    lvef : float                — LVEF (%).
    cmr_findings : dict         — CMR results with keys: lge_pattern (str),
                                  t1_elevated (bool), t2_elevated (bool),
                                  pericardial_effusion (bool).
    ecg_findings : list         — ECG findings (e.g. 'diffuse_st_elevation',
                                  'pr_depression', 'st_changes').
    biopsy_result : str|None    — Endomyocardial biopsy result if performed
                                  (e.g. 'viral', 'giant_cell', 'eosinophilic').
    suspected_etiology : str    — Suspected cause (e.g. 'viral', 'autoimmune',
                                  'drug_reaction', 'idiopathic').
    """

    workflow_type = CardioWorkflowType.MYOCARDITIS_PERICARDITIS

    _VALID_PRESENTATIONS = {"myocarditis", "pericarditis", "myopericarditis"}

    # ── preprocess ────────────────────────────────────────────────────
    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)
        inp.setdefault("presentation", "pericarditis")
        inp.setdefault("symptom_onset_days", 7)
        inp.setdefault("troponin", 0.0)
        inp.setdefault("crp", 0.0)
        inp.setdefault("lvef", 55.0)
        inp.setdefault("cmr_findings", {})
        inp.setdefault("ecg_findings", [])
        inp.setdefault("biopsy_result", None)
        inp.setdefault("suspected_etiology", "viral")
        inp["presentation"] = str(inp["presentation"]).lower().strip()
        inp["suspected_etiology"] = str(inp["suspected_etiology"]).lower().strip()
        # Validate presentation type
        if inp["presentation"] not in self._VALID_PRESENTATIONS:
            warnings.append(
                f"Presentation type '{inp['presentation']}' is not valid "
                "(expected myocarditis, pericarditis, or myopericarditis) "
                "— defaulting to pericarditis"
            )
            inp["presentation"] = "pericarditis"
        # Validate CMR findings is a dict
        if not isinstance(inp["cmr_findings"], dict):
            warnings.append(
                "CMR findings should be a dict — defaulting to empty"
            )
            inp["cmr_findings"] = {}
        # Validate LVEF
        if inp["lvef"] is not None:
            try:
                lv = float(inp["lvef"])
                if not (5 <= lv <= 90):
                    warnings.append(f"LVEF {lv}% outside valid range (5-90%)")
            except (TypeError, ValueError):
                warnings.append(
                    f"LVEF '{inp['lvef']}' is not numeric — defaulting to 55.0"
                )
                inp["lvef"] = 55.0
        # Validate symptom_onset_days
        try:
            days = int(inp["symptom_onset_days"])
            if days < 0:
                warnings.append(
                    f"Symptom onset {days} days is negative — must be >= 0"
                )
                inp["symptom_onset_days"] = 0
        except (TypeError, ValueError):
            warnings.append(
                f"Symptom onset days '{inp['symptom_onset_days']}' is not numeric — defaulting to 7"
            )
            inp["symptom_onset_days"] = 7
        return inp

    # ── execute ───────────────────────────────────────────────────────
    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        presentation = inputs["presentation"]
        onset_days = int(inputs["symptom_onset_days"])
        troponin = float(inputs["troponin"])
        crp = float(inputs["crp"])
        lvef = float(inputs["lvef"])
        cmr = inputs["cmr_findings"]
        ecg_findings = inputs["ecg_findings"]
        biopsy = inputs["biopsy_result"]
        etiology = inputs["suspected_etiology"]

        findings.append(f"Presentation: {presentation}")
        findings.append(f"Symptom onset: {onset_days} days ago")
        findings.append(f"LVEF: {lvef}%")

        if troponin > 0:
            findings.append(
                f"Troponin: {troponin} ng/mL"
                + (" (elevated)" if troponin > 0.04 else "")
            )
        if crp > 0:
            findings.append(
                f"CRP: {crp} mg/L"
                + (" (elevated)" if crp > 3.0 else "")
            )

        # --- CMR-based diagnosis (Modified Lake Louise Criteria) ---
        lge_pattern = cmr.get("lge_pattern", "none")
        t1_elevated = cmr.get("t1_elevated", False)
        t2_elevated = cmr.get("t2_elevated", False)
        pericardial_effusion = cmr.get("pericardial_effusion", False)

        if presentation in ("myocarditis", "myopericarditis"):
            lake_louise_criteria = 0
            if t2_elevated:
                lake_louise_criteria += 1
                findings.append("CMR: T2-weighted edema present (Lake Louise criterion 1)")
            if t1_elevated or (lge_pattern != "none" and lge_pattern):
                lake_louise_criteria += 1
                findings.append(
                    "CMR: T1/ECV elevation or non-ischemic LGE present "
                    "(Lake Louise criterion 2)"
                )
            if lge_pattern not in ("none", "", None) and lge_pattern != "subendocardial":
                findings.append(f"LGE pattern: {lge_pattern} (non-ischemic)")

            if lake_louise_criteria >= 2:
                findings.append(
                    "Modified Lake Louise CMR criteria MET (>=2/3 criteria) "
                    "— supports myocarditis diagnosis"
                )
            else:
                findings.append(
                    f"Modified Lake Louise CMR criteria: {lake_louise_criteria}/3 met "
                    "— consider clinical correlation"
                )

        # --- Pericarditis diagnosis (2/4 criteria) ---
        if presentation in ("pericarditis", "myopericarditis"):
            pericarditis_criteria = 0
            pericarditis_details: List[str] = []

            ecg_lower = [e.lower() for e in ecg_findings]
            if any("pleuritic" in s or "chest_pain" in s for s in ecg_lower):
                pericarditis_criteria += 1
                pericarditis_details.append("pleuritic chest pain")
            if any("rub" in s for s in ecg_lower):
                pericarditis_criteria += 1
                pericarditis_details.append("pericardial rub")
            if any("st_elevation" in s or "pr_depression" in s for s in ecg_lower):
                pericarditis_criteria += 1
                pericarditis_details.append("diffuse ST elevation / PR depression")
            if pericardial_effusion:
                pericarditis_criteria += 1
                pericarditis_details.append("pericardial effusion")

            findings.append(
                f"Pericarditis diagnostic criteria: {pericarditis_criteria}/4 met"
                + (f" ({', '.join(pericarditis_details)})" if pericarditis_details else "")
            )

        # --- Biopsy indications ---
        if biopsy:
            findings.append(f"Endomyocardial biopsy result: {biopsy}")
        else:
            biopsy_reasons: List[str] = []
            if lvef < 40:
                biopsy_reasons.append("new-onset HF with hemodynamic compromise")
            if etiology in ("giant_cell", "eosinophilic"):
                biopsy_reasons.append(f"suspected {etiology} myocarditis")
            if biopsy_reasons:
                recommendations.append(
                    "Endomyocardial biopsy indicated: "
                    + "; ".join(biopsy_reasons)
                )

        # --- Treatment: Myocarditis ---
        if presentation in ("myocarditis", "myopericarditis"):
            recommendations.append(
                "Activity restriction: avoid exercise for 3-6 months "
                "(ACC/AHA guidelines)"
            )
            if lvef < 50:
                recommendations.append(
                    "Initiate HF therapy for reduced LVEF: ACEi/ARB, "
                    "beta-blocker, and consider SGLT2i"
                )
            # Immunosuppression only for specific etiologies
            if biopsy in ("giant_cell", "eosinophilic") or etiology in (
                "giant_cell", "eosinophilic"
            ):
                recommendations.append(
                    f"Immunosuppression indicated for {biopsy or etiology} "
                    "myocarditis (corticosteroids + additional immunosuppressive agent)"
                )
            else:
                recommendations.append(
                    "Immunosuppression NOT recommended for viral myocarditis "
                    "— supportive care only"
                )

        # --- Treatment: Pericarditis ---
        if presentation in ("pericarditis", "myopericarditis"):
            recommendations.append(
                "NSAIDs: ibuprofen 600 mg TID for 1-2 weeks with taper "
                "+ colchicine 0.5 mg BID for 3 months "
                "(COPE/ICAP trials — reduces recurrence by ~50%)"
            )
            recommendations.append(
                "Corticosteroids ONLY if contraindication to NSAIDs "
                "(promotes recurrence — avoid as first-line)"
            )
            if onset_days > 90 or "recurrent" in etiology:
                recommendations.append(
                    "Recurrent pericarditis: long-term colchicine; "
                    "if refractory, consider IL-1 blocker "
                    "(rilonacept per RHAPSODY trial or anakinra)"
                )

        # --- Return-to-play ---
        if presentation in ("myocarditis", "myopericarditis"):
            recommendations.append(
                "Return-to-play criteria: no competitive sports for 3-6 months; "
                "require normal LVEF, no LGE progression, and no arrhythmias "
                "on exercise testing before clearance (ACC/AHA)"
            )

        # --- Severity ---
        is_giant_cell = (biopsy == "giant_cell" or etiology == "giant_cell")
        hemodynamic_instability = lvef < 30
        large_effusion = pericardial_effusion and lvef < 50

        if hemodynamic_instability or is_giant_cell:
            severity = SeverityLevel.CRITICAL
        elif lvef < 50 or large_effusion:
            severity = SeverityLevel.HIGH
        else:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            severity=severity,
            findings=findings,
            recommendations=recommendations,
            cross_modal_triggers=triggers,
            guideline_references=[
                "2024 ACC Expert Consensus Decision Pathway on Myocarditis",
                "2015 ESC Guidelines for Diagnosis and Management of "
                "Pericardial Diseases (Adler et al., Eur Heart J 2015)",
                "Modified Lake Louise CMR Criteria "
                "(Ferreira et al., JACC Cardiovasc Imaging 2018)",
                "COPE Trial (Imazio et al., Circulation 2005)",
                "ICAP Trial (Imazio et al., NEJM 2013)",
                "RHAPSODY Trial (Klein et al., NEJM 2021)",
            ],
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW ENGINE
# ═══════════════════════════════════════════════════════════════════════════


class WorkflowEngine:
    """Central dispatcher that maps CardioWorkflowType to the appropriate
    workflow implementation and handles query-based workflow detection."""

    _KEYWORD_MAP: Dict[str, CardioWorkflowType] = {
        # CAD Assessment
        "coronary": CardioWorkflowType.CAD_ASSESSMENT,
        "cad": CardioWorkflowType.CAD_ASSESSMENT,
        "calcium score": CardioWorkflowType.CAD_ASSESSMENT,
        "cac": CardioWorkflowType.CAD_ASSESSMENT,
        "cadrads": CardioWorkflowType.CAD_ASSESSMENT,
        "cad-rads": CardioWorkflowType.CAD_ASSESSMENT,
        "plaque": CardioWorkflowType.CAD_ASSESSMENT,
        "coronary stenosis": CardioWorkflowType.CAD_ASSESSMENT,
        # Heart Failure
        "heart failure": CardioWorkflowType.HEART_FAILURE,
        "hfref": CardioWorkflowType.HEART_FAILURE,
        "hfpef": CardioWorkflowType.HEART_FAILURE,
        "hfmref": CardioWorkflowType.HEART_FAILURE,
        "ejection fraction": CardioWorkflowType.HEART_FAILURE,
        "gdmt": CardioWorkflowType.HEART_FAILURE,
        "nyha": CardioWorkflowType.HEART_FAILURE,
        "bnp": CardioWorkflowType.HEART_FAILURE,
        "cardiomyopathy": CardioWorkflowType.HEART_FAILURE,
        # Valvular Disease
        "valve": CardioWorkflowType.VALVULAR_DISEASE,
        "valvular": CardioWorkflowType.VALVULAR_DISEASE,
        "aortic stenosis": CardioWorkflowType.VALVULAR_DISEASE,
        "mitral regurgitation": CardioWorkflowType.VALVULAR_DISEASE,
        "tavr": CardioWorkflowType.VALVULAR_DISEASE,
        "savr": CardioWorkflowType.VALVULAR_DISEASE,
        "avr": CardioWorkflowType.VALVULAR_DISEASE,
        # Arrhythmia
        "arrhythmia": CardioWorkflowType.ARRHYTHMIA,
        "atrial fibrillation": CardioWorkflowType.ARRHYTHMIA,
        "afib": CardioWorkflowType.ARRHYTHMIA,
        "a-fib": CardioWorkflowType.ARRHYTHMIA,
        "qtc": CardioWorkflowType.ARRHYTHMIA,
        "qt prolongation": CardioWorkflowType.ARRHYTHMIA,
        "bradycardia": CardioWorkflowType.ARRHYTHMIA,
        "tachycardia": CardioWorkflowType.ARRHYTHMIA,
        "cha2ds2": CardioWorkflowType.ARRHYTHMIA,
        "anticoagulation": CardioWorkflowType.ARRHYTHMIA,
        "pacemaker": CardioWorkflowType.ARRHYTHMIA,
        "ventricular tachycardia": CardioWorkflowType.ARRHYTHMIA,
        "atrial flutter": CardioWorkflowType.ARRHYTHMIA,
        "brugada": CardioWorkflowType.ARRHYTHMIA,
        "long qt": CardioWorkflowType.ARRHYTHMIA,
        "cpvt": CardioWorkflowType.ARRHYTHMIA,
        "channelopathy": CardioWorkflowType.ARRHYTHMIA,
        "wolff-parkinson-white": CardioWorkflowType.ARRHYTHMIA,
        "wpw": CardioWorkflowType.ARRHYTHMIA,
        # Cardiac MRI
        "cardiac mri": CardioWorkflowType.CARDIAC_MRI,
        "cmr": CardioWorkflowType.CARDIAC_MRI,
        "lge": CardioWorkflowType.CARDIAC_MRI,
        "late gadolinium": CardioWorkflowType.CARDIAC_MRI,
        "t1 mapping": CardioWorkflowType.CARDIAC_MRI,
        "t2 mapping": CardioWorkflowType.CARDIAC_MRI,
        "ecv": CardioWorkflowType.CARDIAC_MRI,
        "myocarditis": CardioWorkflowType.CARDIAC_MRI,
        "amyloid": CardioWorkflowType.CARDIAC_MRI,
        # Stress Test
        "stress test": CardioWorkflowType.STRESS_TEST,
        "treadmill": CardioWorkflowType.STRESS_TEST,
        "duke": CardioWorkflowType.STRESS_TEST,
        "exercise test": CardioWorkflowType.STRESS_TEST,
        "perfusion": CardioWorkflowType.STRESS_TEST,
        "nuclear": CardioWorkflowType.STRESS_TEST,
        "ischaemia": CardioWorkflowType.STRESS_TEST,
        "ischemia": CardioWorkflowType.STRESS_TEST,
        # Preventive Risk
        "ascvd": CardioWorkflowType.PREVENTIVE_RISK,
        "prevention": CardioWorkflowType.PREVENTIVE_RISK,
        "preventive": CardioWorkflowType.PREVENTIVE_RISK,
        "statin": CardioWorkflowType.PREVENTIVE_RISK,
        "cholesterol": CardioWorkflowType.PREVENTIVE_RISK,
        "ldl": CardioWorkflowType.PREVENTIVE_RISK,
        "lipid": CardioWorkflowType.PREVENTIVE_RISK,
        "familial hypercholesterolemia": CardioWorkflowType.PREVENTIVE_RISK,
        "risk score": CardioWorkflowType.PREVENTIVE_RISK,
        # Cardio-Oncology
        "cardio-oncology": CardioWorkflowType.CARDIO_ONCOLOGY,
        "cardio oncology": CardioWorkflowType.CARDIO_ONCOLOGY,
        "cardiotoxicity": CardioWorkflowType.CARDIO_ONCOLOGY,
        "ctrcd": CardioWorkflowType.CARDIO_ONCOLOGY,
        "anthracycline": CardioWorkflowType.CARDIO_ONCOLOGY,
        "trastuzumab": CardioWorkflowType.CARDIO_ONCOLOGY,
        "checkpoint inhibitor": CardioWorkflowType.CARDIO_ONCOLOGY,
        "chemotherapy heart": CardioWorkflowType.CARDIO_ONCOLOGY,
        "gls": CardioWorkflowType.CARDIO_ONCOLOGY,
        "dexrazoxane": CardioWorkflowType.CARDIO_ONCOLOGY,
        # Acute Decompensated HF
        "acute decompensated": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "adhf": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "acute heart failure": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "iv diuretic": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "cardiogenic shock": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "inotrope": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "warm-wet": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "cold-wet": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        "impella": CardioWorkflowType.ACUTE_DECOMPENSATED_HF,
        # Post-MI
        "post-mi": CardioWorkflowType.POST_MI,
        "post mi": CardioWorkflowType.POST_MI,
        "myocardial infarction": CardioWorkflowType.POST_MI,
        "dapt": CardioWorkflowType.POST_MI,
        "dual antiplatelet": CardioWorkflowType.POST_MI,
        "door-to-balloon": CardioWorkflowType.POST_MI,
        "ticagrelor": CardioWorkflowType.POST_MI,
        "cardiac rehab": CardioWorkflowType.POST_MI,
        # Myocarditis / Pericarditis
        "myocarditis": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,  # noqa: F601
        "pericarditis": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "myopericarditis": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "pericardial effusion": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "lake louise": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "colchicine": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "giant cell": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "return-to-play": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
        "rilonacept": CardioWorkflowType.MYOCARDITIS_PERICARDITIS,
    }

    def __init__(self) -> None:
        workflow_instances: List[BaseCardioWorkflow] = [
            CADAssessmentWorkflow(),
            HeartFailureWorkflow(),
            ValvularDiseaseWorkflow(),
            ArrhythmiaWorkflow(),
            CardiacMRIWorkflow(),
            StressTestWorkflow(),
            PreventiveRiskWorkflow(),
            CardioOncologyWorkflow(),
            AcuteDecompensatedHFWorkflow(),
            PostMIWorkflow(),
            MyocarditisPericarditisWorkflow(),
        ]
        self._workflows: Dict[CardioWorkflowType, BaseCardioWorkflow] = {
            wf.workflow_type: wf for wf in workflow_instances
        }

    # ── public API ────────────────────────────────────────────────────

    def run_workflow(
        self, workflow_type: CardioWorkflowType, inputs: dict
    ) -> WorkflowResult:
        """Execute a specific workflow by type."""
        wf = self._workflows.get(workflow_type)
        if wf is None:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
        return wf.run(inputs)

    def detect_workflow(self, query: str) -> Optional[CardioWorkflowType]:
        """Detect the most relevant workflow from a free-text query.

        Returns the workflow type with the most keyword matches, or None
        if no keywords match.
        """
        query_lower = query.lower()
        scores: Dict[CardioWorkflowType, int] = {}
        for keyword, wf_type in self._KEYWORD_MAP.items():
            if keyword in query_lower:
                scores[wf_type] = scores.get(wf_type, 0) + 1
        if not scores:
            return None
        return max(scores, key=scores.get)  # type: ignore[arg-type]

    def get_available_workflows(self) -> List[CardioWorkflowType]:
        """Return all registered workflow types."""
        return list(self._workflows.keys())
