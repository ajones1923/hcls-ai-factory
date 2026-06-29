"""Clinical workflows for the Clinical Trial Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements ten evidence-based clinical trial workflows that integrate protocol
design, patient matching, site selection, eligibility optimization, adaptive
design, safety signal detection, regulatory document generation, competitive
intelligence, diversity assessment, and decentralized trial planning.

Each workflow follows the BaseTrialWorkflow contract
(preprocess -> execute -> postprocess) and is registered in the
WorkflowEngine for unified dispatch.
"""

from __future__ import annotations

import logging
import math
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from typing import Dict, List

from src.models import (
    CriterionType,
    DCTComponent,
    MatchScore,
    SeverityLevel,
    TrialWorkflowType,
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


class BaseTrialWorkflow(ABC):
    """Abstract base for all clinical trial workflows."""

    workflow_type: TrialWorkflowType

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
        """Shared enrichment after execution.  Override to add workflow-specific post-steps."""
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
# WORKFLOW 1 — Protocol Design
# ═══════════════════════════════════════════════════════════════════════════


class ProtocolDesignWorkflow(BaseTrialWorkflow):
    """Protocol design workflow that generates evidence-based protocol
    blueprints by analysing historical protocols, success patterns, and
    endpoint benchmarks for a given indication/phase combination.

    Inputs
    ------
    indication : str
        Target disease or condition.
    phase : str
        Trial phase (phase_i, phase_ii, phase_iii, etc.).
    comparator : str
        Comparator arm description (e.g. 'placebo', 'standard_of_care').
    target_population : str
        Description of the target patient population.
    mechanism_of_action : str
        Drug mechanism of action.
    """

    workflow_type = TrialWorkflowType.PROTOCOL_DESIGN

    # Historical success rates by indication category and phase
    _HISTORICAL_SUCCESS: Dict[str, Dict[str, float]] = {
        "oncology": {
            "phase_i": 0.63, "phase_ii": 0.29, "phase_iii": 0.36, "phase_iv": 0.85,
        },
        "cardiology": {
            "phase_i": 0.72, "phase_ii": 0.38, "phase_iii": 0.52, "phase_iv": 0.88,
        },
        "neurology": {
            "phase_i": 0.68, "phase_ii": 0.25, "phase_iii": 0.30, "phase_iv": 0.82,
        },
        "immunology": {
            "phase_i": 0.70, "phase_ii": 0.35, "phase_iii": 0.48, "phase_iv": 0.86,
        },
        "rare_disease": {
            "phase_i": 0.74, "phase_ii": 0.42, "phase_iii": 0.55, "phase_iv": 0.90,
        },
        "infectious_disease": {
            "phase_i": 0.71, "phase_ii": 0.40, "phase_iii": 0.58, "phase_iv": 0.87,
        },
        "default": {
            "phase_i": 0.67, "phase_ii": 0.33, "phase_iii": 0.45, "phase_iv": 0.85,
        },
    }

    # Recommended primary endpoints by indication
    _ENDPOINT_BENCHMARKS: Dict[str, Dict[str, str]] = {
        "oncology": {
            "primary": "Overall Survival (OS) or Progression-Free Survival (PFS)",
            "secondary": "Objective Response Rate (ORR), Duration of Response (DoR)",
            "safety": "CTCAE v5.0 graded adverse events",
        },
        "cardiology": {
            "primary": "Major Adverse Cardiovascular Events (MACE)",
            "secondary": "Change in LVEF, NT-proBNP levels, 6MWD",
            "safety": "Bleeding events (BARC criteria), arrhythmias",
        },
        "neurology": {
            "primary": "Disease-specific functional scale (e.g. EDSS, ADAS-Cog)",
            "secondary": "MRI lesion volume, biomarker levels",
            "safety": "PML screening, hepatotoxicity monitoring",
        },
        "default": {
            "primary": "Disease-specific validated outcome measure",
            "secondary": "Quality of life (EQ-5D), biomarker response",
            "safety": "Treatment-emergent adverse events (TEAEs)",
        },
    }

    # Design complexity factors
    _DESIGN_COMPLEXITY: Dict[str, float] = {
        "open_label": 0.3, "single_blind": 0.5, "double_blind": 0.7,
        "adaptive": 0.85, "basket": 0.8, "umbrella": 0.85, "platform": 0.95,
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        # Required fields with defaults
        if not inp.get("indication"):
            warnings.append("No indication provided — using 'unspecified'")
            inp["indication"] = "unspecified"
        if not inp.get("phase"):
            warnings.append("No phase provided — defaulting to phase_ii")
            inp["phase"] = "phase_ii"
        inp.setdefault("comparator", "placebo")
        inp.setdefault("target_population", "adult patients")
        inp.setdefault("mechanism_of_action", "unspecified")
        inp.setdefault("design_type", "double_blind")
        inp.setdefault("num_arms", 2)

        # Normalise phase
        phase = str(inp["phase"]).lower().replace(" ", "_").replace("-", "_")
        if not phase.startswith("phase_"):
            phase = f"phase_{phase}"
        valid_phases = {"phase_i", "phase_i_ii", "phase_ii", "phase_ii_iii",
                        "phase_iii", "phase_iv"}
        if phase not in valid_phases:
            warnings.append(
                f"Phase '{inp['phase']}' not recognised — defaulting to phase_ii"
            )
            phase = "phase_ii"
        inp["phase"] = phase

        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        references: List[str] = []

        indication = inputs["indication"].lower()
        phase = inputs["phase"]
        comparator = inputs["comparator"]
        moa = inputs["mechanism_of_action"]
        design_type = inputs.get("design_type", "double_blind")
        num_arms = inputs.get("num_arms", 2)

        # --- Determine therapeutic category ---
        category = "default"
        for cat in self._HISTORICAL_SUCCESS:
            if cat != "default" and cat in indication:
                category = cat
                break

        # --- Success probability ---
        base_success = self._HISTORICAL_SUCCESS.get(
            category, self._HISTORICAL_SUCCESS["default"]
        ).get(phase, 0.40)

        # Adjust for design factors
        design_factor = self._DESIGN_COMPLEXITY.get(design_type, 0.5)
        # More complex designs have slightly lower success but better outcomes
        adjusted_success = base_success * (1.0 - 0.1 * design_factor)

        # Adjust for biomarker-driven enrichment
        if any(kw in moa.lower() for kw in ["targeted", "biomarker", "precision"]):
            adjusted_success *= 1.15
            findings.append(
                "Biomarker-driven enrichment strategy detected — "
                f"success probability elevated by 15% (base: {base_success:.0%})"
            )

        adjusted_success = _clamp(adjusted_success, 0.05, 0.95)

        findings.append(
            f"Historical success rate for {category}/{phase}: {base_success:.0%}"
        )
        findings.append(
            f"Adjusted success probability: {adjusted_success:.0%} "
            f"(design: {design_type}, arms: {num_arms})"
        )

        # --- Endpoint benchmarks ---
        endpoints = self._ENDPOINT_BENCHMARKS.get(
            category, self._ENDPOINT_BENCHMARKS["default"]
        )
        findings.append(f"Recommended primary endpoint: {endpoints['primary']}")
        findings.append(f"Recommended secondary endpoints: {endpoints['secondary']}")
        findings.append(f"Safety monitoring: {endpoints['safety']}")

        # --- Sample size estimation (simplified) ---
        phase_sample_ranges = {
            "phase_i": (20, 80), "phase_i_ii": (50, 150),
            "phase_ii": (100, 300), "phase_ii_iii": (200, 600),
            "phase_iii": (300, 3000), "phase_iv": (500, 10000),
        }
        lo, hi = phase_sample_ranges.get(phase, (100, 500))
        estimated_n = int((lo + hi) / 2 * num_arms / 2)
        findings.append(
            f"Estimated sample size: {estimated_n} "
            f"(range: {lo * num_arms // 2}-{hi * num_arms // 2} for {num_arms}-arm design)"
        )

        # --- Protocol complexity ---
        complexity_score = _clamp(
            design_factor * 0.4 + (num_arms / 5) * 0.3 + (estimated_n / 3000) * 0.3
        )
        findings.append(f"Protocol complexity score: {complexity_score:.2f}")

        # --- Recommendations ---
        recommendations.append(
            f"Design a {design_type} {phase.replace('_', ' ')} trial with "
            f"{num_arms} arms ({comparator} comparator)"
        )
        if adjusted_success < 0.35:
            recommendations.append(
                "Consider adaptive design with interim futility analysis "
                "to mitigate low expected success rate"
            )
        if complexity_score > 0.7:
            recommendations.append(
                "High protocol complexity — consider simplifying visit "
                "schedules and reducing secondary endpoints"
            )
        if "placebo" in comparator.lower() and phase in ("phase_iii", "phase_iv"):
            recommendations.append(
                "Evaluate ethical appropriateness of placebo comparator "
                "if standard-of-care exists (Helsinki Declaration)"
            )

        recommendations.append(
            f"Primary endpoint: {endpoints['primary']}"
        )

        # --- References ---
        references.extend([
            "ICH E6(R3) Good Clinical Practice (2023)",
            "ICH E8(R1) General Considerations for Clinical Studies (2021)",
            "ICH E9(R1) Estimands and Sensitivity Analysis",
            f"FDA Guidance: {category.title()} Trial Design (current)",
        ])

        severity = SeverityLevel.INFORMATIONAL
        if adjusted_success < 0.30:
            severity = SeverityLevel.HIGH
        elif adjusted_success < 0.45:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=severity,
            confidence=_clamp(0.5 + adjusted_success * 0.3),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 2 — Patient Matching
# ═══════════════════════════════════════════════════════════════════════════


class PatientMatchingWorkflow(BaseTrialWorkflow):
    """Patient-trial matching workflow that evaluates a patient profile
    against trial eligibility criteria and ranks trials by match score.

    Inputs
    ------
    patient : dict
        Patient demographics, biomarkers, genomic data, medications,
        comorbidities, and location (PatientProfile fields).
    trials : list[dict]
        List of trial dicts, each containing: trial_id, title, phase,
        status, inclusion_criteria (list[str]), exclusion_criteria (list[str]).
    """

    workflow_type = TrialWorkflowType.PATIENT_MATCHING

    # Age-related keyword patterns
    _AGE_PATTERNS = {
        "adult": (18, 120), "pediatric": (0, 17), "elderly": (65, 120),
        "geriatric": (65, 120), "adolescent": (12, 17), "child": (0, 11),
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("patient"):
            warnings.append("No patient profile provided — using empty profile")
            inp["patient"] = {}
        if not inp.get("trials"):
            warnings.append("No trials provided — no matching will be performed")
            inp["trials"] = []

        # Normalise patient
        p = inp["patient"]
        p.setdefault("age", None)
        p.setdefault("sex", None)
        p.setdefault("diagnosis", None)
        p.setdefault("biomarkers", [])
        p.setdefault("medications", [])
        p.setdefault("genomic_variants", [])
        p.setdefault("comorbidities", [])
        p.setdefault("geographic_location", None)

        if p["age"] is not None:
            try:
                p["age"] = int(p["age"])
                if not (0 <= p["age"] <= 120):
                    warnings.append(
                        f"Patient age {p['age']} out of range (0-120) — setting to None"
                    )
                    p["age"] = None
            except (TypeError, ValueError):
                warnings.append(
                    f"Patient age '{p['age']}' is not numeric — setting to None"
                )
                p["age"] = None

        return inp

    def _evaluate_criterion(
        self, criterion: str, patient: dict, is_exclusion: bool
    ) -> MatchScore:
        """Evaluate a single eligibility criterion against the patient profile."""
        crit_lower = criterion.lower()
        met = True
        confidence = 0.5  # Default uncertain
        evidence = ""

        # --- Age checks ---
        age = patient.get("age")
        if age is not None:
            if "age" in crit_lower or "years" in crit_lower or "≥" in criterion:
                # Try to extract age bounds
                for pattern, (lo, hi) in self._AGE_PATTERNS.items():
                    if pattern in crit_lower:
                        met = lo <= age <= hi
                        confidence = 0.9
                        evidence = f"Patient age {age}, criterion requires {pattern} ({lo}-{hi})"
                        break
                else:
                    # Check for numeric bounds
                    import re
                    nums = re.findall(r"(\d+)\s*(?:years|yrs|y\.o\.)", crit_lower)
                    if len(nums) >= 2:
                        lo, hi = int(nums[0]), int(nums[1])
                        met = lo <= age <= hi
                        confidence = 0.9
                        evidence = f"Patient age {age}, criterion range {lo}-{hi}"
                    elif nums:
                        bound = int(nums[0])
                        if any(kw in crit_lower for kw in ["≥", ">=", "at least", "older"]):
                            met = age >= bound
                            confidence = 0.9
                            evidence = f"Patient age {age} vs minimum {bound}"
                        elif any(kw in crit_lower for kw in ["≤", "<=", "no more", "younger"]):
                            met = age <= bound
                            confidence = 0.9
                            evidence = f"Patient age {age} vs maximum {bound}"

        # --- Sex checks ---
        sex = patient.get("sex")
        if sex and ("male" in crit_lower or "female" in crit_lower):
            if "female" in crit_lower and "only" in crit_lower:
                met = sex.lower() == "female"
                confidence = 0.95
                evidence = f"Patient sex: {sex}, criterion requires female"
            elif "male" in crit_lower and "only" in crit_lower:
                met = sex.lower() == "male"
                confidence = 0.95
                evidence = f"Patient sex: {sex}, criterion requires male"

        # --- Diagnosis check ---
        diagnosis = (patient.get("diagnosis") or "").lower()
        if diagnosis and any(
            kw in crit_lower
            for kw in ["diagnosis", "confirmed", "histologically", "documented"]
        ):
            if diagnosis and any(
                word in crit_lower for word in diagnosis.split()
                if len(word) > 3
            ):
                met = True
                confidence = 0.7
                evidence = f"Patient diagnosis '{diagnosis}' matches criterion keywords"
            else:
                met = False
                confidence = 0.5
                evidence = f"Patient diagnosis '{diagnosis}' may not match criterion"

        # --- Biomarker checks ---
        biomarkers = [b.lower() for b in patient.get("biomarkers", [])]
        if biomarkers:
            for bm in biomarkers:
                if bm in crit_lower:
                    met = True
                    confidence = 0.85
                    evidence = f"Biomarker '{bm}' found in criterion"
                    break

        # --- Comorbidity exclusion checks ---
        comorbidities = [c.lower() for c in patient.get("comorbidities", [])]
        if is_exclusion and comorbidities:
            for comor in comorbidities:
                if comor in crit_lower or any(
                    word in crit_lower for word in comor.split() if len(word) > 4
                ):
                    met = True  # exclusion criterion triggered
                    confidence = 0.75
                    evidence = f"Comorbidity '{comor}' matches exclusion criterion"
                    break

        # --- Medication exclusion checks ---
        medications = [m.lower() for m in patient.get("medications", [])]
        if is_exclusion and medications:
            for med in medications:
                if med in crit_lower:
                    met = True  # exclusion criterion triggered
                    confidence = 0.8
                    evidence = f"Medication '{med}' matches exclusion criterion"
                    break

        crit_type = CriterionType.EXCLUSION if is_exclusion else CriterionType.INCLUSION
        return MatchScore(
            criterion_text=criterion,
            criterion_type=crit_type,
            met=met,
            confidence=confidence,
            evidence=evidence,
        )

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        patient = inputs["patient"]
        trials = inputs.get("trials", [])

        if not trials:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=["No trials provided for matching"],
                recommendations=["Provide trial eligibility criteria for matching"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        matches: List[dict] = []
        for trial in trials:
            trial_id = trial.get("trial_id", "UNKNOWN")
            title = trial.get("title", "Untitled")
            phase = trial.get("phase", "phase_ii")
            status = trial.get("status", "recruiting")
            inclusion = trial.get("inclusion_criteria", [])
            exclusion = trial.get("exclusion_criteria", [])

            # Evaluate inclusion criteria
            inc_scores = [
                self._evaluate_criterion(c, patient, is_exclusion=False)
                for c in inclusion
            ]
            inc_met = sum(1 for s in inc_scores if s.met)

            # Evaluate exclusion criteria
            exc_scores = [
                self._evaluate_criterion(c, patient, is_exclusion=True)
                for c in exclusion
            ]
            # For exclusion: "clear" means NOT triggered
            exc_clear = sum(1 for s in exc_scores if not s.met)

            # Overall score
            inc_ratio = inc_met / max(len(inclusion), 1)
            exc_ratio = exc_clear / max(len(exclusion), 1)
            overall = 0.6 * inc_ratio + 0.4 * exc_ratio

            avg_confidence = 0.5
            all_scores = inc_scores + exc_scores
            if all_scores:
                avg_confidence = sum(s.confidence for s in all_scores) / len(all_scores)

            matches.append({
                "trial_id": trial_id,
                "title": title,
                "phase": phase,
                "status": status,
                "inclusion_met": inc_met,
                "inclusion_total": len(inclusion),
                "exclusion_clear": exc_clear,
                "exclusion_total": len(exclusion),
                "overall_score": round(overall, 3),
                "confidence": round(avg_confidence, 3),
            })

        # Sort by overall score descending
        matches.sort(key=lambda m: m["overall_score"], reverse=True)

        findings.append(f"Evaluated {len(trials)} trials for patient matching")
        for i, m in enumerate(matches[:5], 1):
            findings.append(
                f"  #{i}: {m['trial_id']} — score {m['overall_score']:.2f} "
                f"(inclusion: {m['inclusion_met']}/{m['inclusion_total']}, "
                f"exclusion clear: {m['exclusion_clear']}/{m['exclusion_total']})"
            )

        if matches and matches[0]["overall_score"] >= 0.7:
            recommendations.append(
                f"Top match: {matches[0]['trial_id']} — "
                f"consider immediate pre-screening"
            )
        elif matches:
            recommendations.append(
                "No high-confidence matches found — consider broadening search "
                "criteria or reviewing waiver options"
            )

        # Determine overall severity
        top_score = matches[0]["overall_score"] if matches else 0
        severity = SeverityLevel.INFORMATIONAL
        if top_score >= 0.8:
            severity = SeverityLevel.LOW
        elif top_score >= 0.5:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "FDA Guidance: Enrichment Strategies for Clinical Trials (2019)",
                "ICH E6(R3) GCP — Eligibility Assessment",
            ],
            severity=severity,
            confidence=_clamp(top_score * 0.8 + 0.1),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 3 — Site Selection
# ═══════════════════════════════════════════════════════════════════════════


class SiteSelectionWorkflow(BaseTrialWorkflow):
    """Site selection and ranking workflow based on weighted scoring
    across seven performance dimensions.

    Inputs
    ------
    sites : list[dict]
        Candidate site data with fields: site_id, facility_name, city,
        country, enrollment_rate, screen_failure_rate, investigator_h_index,
        therapeutic_experience (years), population_access (int),
        diversity_index (0-1), regulatory_readiness (0-1).
    target_enrollment : int
        Target total enrollment.
    therapeutic_area : str
        Trial therapeutic area.
    """

    workflow_type = TrialWorkflowType.SITE_SELECTION

    # Weight configuration for scoring dimensions
    _WEIGHTS = {
        "enrollment_rate": 0.25,
        "screen_failure_rate": 0.15,
        "investigator_h_index": 0.10,
        "therapeutic_experience": 0.15,
        "population_access": 0.15,
        "diversity_index": 0.10,
        "regulatory_readiness": 0.10,
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("sites"):
            warnings.append("No sites provided — returning empty rankings")
            inp["sites"] = []
        inp.setdefault("target_enrollment", 300)
        inp.setdefault("therapeutic_area", "general")

        try:
            inp["target_enrollment"] = int(inp["target_enrollment"])
            if inp["target_enrollment"] <= 0:
                warnings.append("Target enrollment must be positive — defaulting to 300")
                inp["target_enrollment"] = 300
        except (TypeError, ValueError):
            warnings.append("Invalid target enrollment — defaulting to 300")
            inp["target_enrollment"] = 300

        return inp

    def _score_site(self, site: dict) -> Dict[str, float]:
        """Score a site across all dimensions (each normalised to 0-1)."""
        scores: Dict[str, float] = {}

        # Enrollment rate: normalise to 0-1 (assume max 10 pts/month is excellent)
        er = float(site.get("enrollment_rate", 0))
        scores["enrollment_rate"] = _clamp(er / 10.0)

        # Screen failure: lower is better
        sfr = float(site.get("screen_failure_rate", 0.5))
        scores["screen_failure_rate"] = _clamp(1.0 - sfr)

        # Investigator H-index: normalise (assume max 80)
        h_idx = float(site.get("investigator_h_index", 10))
        scores["investigator_h_index"] = _clamp(h_idx / 80.0)

        # Therapeutic experience: years, normalise (assume max 30)
        te = float(site.get("therapeutic_experience", 5))
        scores["therapeutic_experience"] = _clamp(te / 30.0)

        # Population access: normalise (assume max 5M catchment)
        pa = float(site.get("population_access", 100000))
        scores["population_access"] = _clamp(pa / 5_000_000)

        # Diversity index: already 0-1
        di = float(site.get("diversity_index", 0.3))
        scores["diversity_index"] = _clamp(di)

        # Regulatory readiness: already 0-1
        rr = float(site.get("regulatory_readiness", 0.5))
        scores["regulatory_readiness"] = _clamp(rr)

        return scores

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        sites = inputs.get("sites", [])
        target = inputs["target_enrollment"]

        if not sites:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=["No candidate sites provided"],
                recommendations=["Provide site data for selection analysis"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        scored_sites: List[dict] = []
        for site in sites:
            dim_scores = self._score_site(site)
            weighted = sum(
                dim_scores[dim] * self._WEIGHTS[dim] for dim in self._WEIGHTS
            )
            scored_sites.append({
                "site_id": site.get("site_id", "UNKNOWN"),
                "facility_name": site.get("facility_name", "Unknown Facility"),
                "city": site.get("city", ""),
                "country": site.get("country", ""),
                "enrollment_rate": float(site.get("enrollment_rate", 0)),
                "screen_failure_rate": float(site.get("screen_failure_rate", 0.5)),
                "diversity_index": float(site.get("diversity_index", 0.3)),
                "overall_score": round(weighted, 4),
                "dimension_scores": dim_scores,
            })

        # Sort by overall score
        scored_sites.sort(key=lambda s: s["overall_score"], reverse=True)

        findings.append(
            f"Evaluated {len(scored_sites)} sites for trial "
            f"(target enrollment: {target})"
        )

        # Estimate required number of sites
        avg_rate = sum(s["enrollment_rate"] for s in scored_sites) / len(scored_sites)
        months = 18  # typical enrollment window
        if avg_rate > 0:
            needed = math.ceil(target / (avg_rate * months))
        else:
            needed = len(scored_sites)
        findings.append(
            f"Estimated sites needed: {needed} "
            f"(avg rate: {avg_rate:.1f} pts/mo, {months}-month window)"
        )

        for i, s in enumerate(scored_sites[:10], 1):
            findings.append(
                f"  #{i}: {s['site_id']} ({s['facility_name']}, {s['city']}, "
                f"{s['country']}) — score {s['overall_score']:.3f}"
            )

        # Recommendations
        if needed > len(scored_sites):
            recommendations.append(
                f"Insufficient sites: need ~{needed} but only {len(scored_sites)} "
                f"candidates — expand site network"
            )
        low_diversity = [
            s for s in scored_sites[:needed]
            if s["diversity_index"] < 0.3
        ]
        if low_diversity:
            recommendations.append(
                f"{len(low_diversity)} top-ranked sites have low diversity index "
                f"(<0.3) — supplement with demographically diverse sites"
            )

        recommendations.append(
            f"Select top {min(needed, len(scored_sites))} sites for activation"
        )

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "ICH E6(R3) GCP — Investigator Qualifications",
                "FDA Guidance: Enhancing the Diversity of Clinical Trial Populations (2020)",
                "TransCelerate Clinical Operations Best Practices",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=_clamp(0.6 + len(sites) * 0.01),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 4 — Eligibility Optimization
# ═══════════════════════════════════════════════════════════════════════════


class EligibilityOptimizationWorkflow(BaseTrialWorkflow):
    """Optimise eligibility criteria by identifying overly restrictive
    criteria and recommending broadening strategies.

    Inputs
    ------
    eligibility_criteria : list[dict]
        List of dicts with keys: text, type (inclusion/exclusion),
        scientific_justification (optional).
    indication : str
        Target indication for competitor comparison.
    """

    workflow_type = TrialWorkflowType.ELIGIBILITY_OPTIMIZATION

    # Known overly restrictive patterns and their typical population impact
    _RESTRICTIVE_PATTERNS: Dict[str, Dict[str, float]] = {
        "ecog 0-1": {"impact": 0.25, "note": "Excludes ECOG 2 patients who may benefit"},
        "no prior therapy": {"impact": 0.35, "note": "Excludes pretreated population"},
        "no organ dysfunction": {"impact": 0.20, "note": "Excludes common comorbidities"},
        "hemoglobin": {"impact": 0.12, "note": "May exclude anemic patients unnecessarily"},
        "creatinine clearance": {"impact": 0.15, "note": "May exclude mild renal impairment"},
        "hepatic function": {"impact": 0.10, "note": "May over-restrict liver function requirements"},
        "no autoimmune": {"impact": 0.08, "note": "Excludes controlled autoimmune conditions"},
        "no brain metastases": {"impact": 0.15, "note": "Excludes treated/stable CNS disease"},
        "washout period": {"impact": 0.10, "note": "Extended washout may delay enrollment"},
        "cardiac ejection fraction": {"impact": 0.08, "note": "May exclude subclinical dysfunction"},
        "bmi": {"impact": 0.12, "note": "BMI restrictions exclude significant populations"},
        "age": {"impact": 0.18, "note": "Narrow age ranges limit generalisability"},
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("eligibility_criteria"):
            warnings.append(
                "No eligibility criteria provided — returning empty analysis"
            )
            inp["eligibility_criteria"] = []
        inp.setdefault("indication", "unspecified")
        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        criteria = inputs["eligibility_criteria"]

        if not criteria:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=["No eligibility criteria provided for analysis"],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        total_impact = 0.0
        restrictive_criteria: List[dict] = []

        for crit in criteria:
            text = crit.get("text", "") if isinstance(crit, dict) else str(crit)
            crit_type = crit.get("type", "inclusion") if isinstance(crit, dict) else "inclusion"
            justification = crit.get("scientific_justification", "") if isinstance(crit, dict) else ""
            text_lower = text.lower()

            # Check against known restrictive patterns
            population_impact = 0.0
            matched_pattern = None
            for pattern, info in self._RESTRICTIVE_PATTERNS.items():
                if pattern in text_lower:
                    population_impact = info["impact"]
                    matched_pattern = pattern
                    break

            # If no known pattern, estimate based on criterion type
            if not matched_pattern:
                population_impact = 0.05 if crit_type == "inclusion" else 0.03

            # Assess scientific justification strength
            justification_score = 0.5  # Default moderate
            if justification:
                just_lower = justification.lower()
                if any(kw in just_lower for kw in ["rct", "meta-analysis", "phase iii"]):
                    justification_score = 0.9
                elif any(kw in just_lower for kw in ["observational", "cohort", "registry"]):
                    justification_score = 0.6
                elif any(kw in just_lower for kw in ["expert", "opinion", "convention"]):
                    justification_score = 0.3
            else:
                justification_score = 0.2

            total_impact += population_impact

            if population_impact >= 0.10 and justification_score < 0.6:
                restrictive_criteria.append({
                    "text": text,
                    "impact": population_impact,
                    "justification_score": justification_score,
                    "pattern": matched_pattern,
                })

        # Compound impact (assuming independence)
        compound_exclusion = 1.0 - math.prod(
            1.0 - rc["impact"] for rc in restrictive_criteria
        ) if restrictive_criteria else 0.0

        findings.append(
            f"Analysed {len(criteria)} eligibility criteria"
        )
        findings.append(
            f"Identified {len(restrictive_criteria)} potentially "
            f"overly restrictive criteria"
        )
        findings.append(
            f"Estimated compound population exclusion: {compound_exclusion:.0%}"
        )

        for rc in sorted(restrictive_criteria, key=lambda x: x["impact"], reverse=True):
            findings.append(
                f"  - '{rc['text'][:80]}...' — impact: {rc['impact']:.0%}, "
                f"justification: {rc['justification_score']:.0%}"
            )
            note = (
                self._RESTRICTIVE_PATTERNS.get(rc["pattern"], {}).get("note", "")
                if rc["pattern"]
                else ""
            )
            if note:
                recommendations.append(
                    f"Consider broadening: {note} "
                    f"(predicted population gain: {rc['impact']:.0%})"
                )

        if compound_exclusion > 0.40:
            recommendations.append(
                "ALERT: Cumulative eligibility restrictions may exclude >40% "
                "of the target population — review for FDA Broadening Eligibility Criteria guidance"
            )

        severity = SeverityLevel.INFORMATIONAL
        if compound_exclusion > 0.40:
            severity = SeverityLevel.HIGH
        elif compound_exclusion > 0.20:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "FDA Guidance: Enhancing the Diversity of Clinical Trial Populations (2020)",
                "FDA Guidance: Broadening Eligibility Criteria (2019)",
                "ASCO-Friends Broadening Eligibility Criteria Project",
            ],
            severity=severity,
            confidence=_clamp(0.55 + len(criteria) * 0.01),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 5 — Adaptive Design
# ═══════════════════════════════════════════════════════════════════════════


class AdaptiveDesignWorkflow(BaseTrialWorkflow):
    """Adaptive trial design workflow that evaluates interim data against
    pre-specified decision rules and recommends design adaptations.

    Inputs
    ------
    design_type : str
        One of: group_sequential, sample_size_reestimation,
        response_adaptive, biomarker_adaptive, platform, seamless.
    interim_data : dict
        Interim analysis results (effect_size, p_value, n_enrolled,
        n_events, response_rates per arm).
    enrollment_status : dict
        Current enrollment (enrolled, target, rate_per_month).
    safety_data : dict
        Accumulated safety (sae_count, discontinuation_rate,
        dose_modifications).
    decision_boundaries : dict
        Pre-specified boundaries (futility_boundary, efficacy_boundary,
        alpha_spent).
    """

    workflow_type = TrialWorkflowType.ADAPTIVE_DESIGN

    _VALID_DESIGNS = {
        "group_sequential", "sample_size_reestimation",
        "response_adaptive", "biomarker_adaptive", "platform", "seamless",
    }

    # Default alpha spending functions (O'Brien-Fleming-like boundaries)
    _DEFAULT_BOUNDARIES = {
        "futility_boundary": 0.50,
        "efficacy_boundary": 0.001,
        "alpha_spent": 0.005,
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        design = inp.get("design_type", "group_sequential")
        if design not in self._VALID_DESIGNS:
            warnings.append(
                f"Design type '{design}' not recognised — "
                f"defaulting to group_sequential"
            )
            design = "group_sequential"
        inp["design_type"] = design

        inp.setdefault("interim_data", {})
        inp.setdefault("enrollment_status", {})
        inp.setdefault("safety_data", {})
        inp.setdefault("decision_boundaries", dict(self._DEFAULT_BOUNDARIES))

        # Validate interim data
        interim = inp["interim_data"]
        interim.setdefault("effect_size", None)
        interim.setdefault("p_value", None)
        interim.setdefault("n_enrolled", 0)
        interim.setdefault("n_events", 0)
        interim.setdefault("response_rates", {})

        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        design = inputs["design_type"]
        interim = inputs["interim_data"]
        enrollment = inputs["enrollment_status"]
        safety = inputs["safety_data"]
        boundaries = inputs["decision_boundaries"]

        findings.append(f"Adaptive design type: {design.replace('_', ' ').title()}")

        # --- Enrollment assessment ---
        enrolled = enrollment.get("enrolled", 0)
        target = enrollment.get("target", 100)
        rate = enrollment.get("rate_per_month", 0)
        info_fraction = enrolled / max(target, 1)
        findings.append(
            f"Information fraction: {info_fraction:.0%} "
            f"({enrolled}/{target} enrolled)"
        )

        if rate > 0 and enrolled < target:
            months_remaining = (target - enrolled) / rate
            findings.append(
                f"Estimated time to full enrollment: {months_remaining:.1f} months "
                f"(rate: {rate:.1f}/month)"
            )

        # --- Futility assessment ---
        p_value = interim.get("p_value")
        effect_size = interim.get("effect_size")
        futility_bound = boundaries.get(
            "futility_boundary", self._DEFAULT_BOUNDARIES["futility_boundary"]
        )
        efficacy_bound = boundaries.get(
            "efficacy_boundary", self._DEFAULT_BOUNDARIES["efficacy_boundary"]
        )

        if p_value is not None:
            findings.append(f"Interim p-value: {p_value:.4f}")
            if info_fraction >= 0.5 and p_value > futility_bound:
                findings.append(
                    f"FUTILITY BOUNDARY CROSSED: p={p_value:.4f} > {futility_bound} "
                    f"at {info_fraction:.0%} information fraction"
                )
                recommendations.append(
                    "Consider stopping for futility — conditional power is low"
                )
            elif p_value < efficacy_bound:
                findings.append(
                    f"EFFICACY BOUNDARY CROSSED: p={p_value:.4f} < {efficacy_bound}"
                )
                recommendations.append(
                    "Consider early stopping for overwhelming efficacy"
                )

        if effect_size is not None:
            findings.append(f"Observed effect size: {effect_size:.3f}")

        # --- Safety assessment ---
        sae_count = safety.get("sae_count", 0)
        disc_rate = safety.get("discontinuation_rate", 0)
        safety.get("dose_modifications", 0)

        if sae_count > 0:
            sae_rate = sae_count / max(enrolled, 1)
            findings.append(
                f"SAE rate: {sae_rate:.1%} ({sae_count}/{enrolled})"
            )
            if sae_rate > 0.10:
                recommendations.append(
                    f"Elevated SAE rate ({sae_rate:.1%}) — recommend "
                    f"DSMB review and potential protocol amendment"
                )
                triggers.append(
                    _trigger_string("SAFETY", "Elevated SAE rate requires DSMB review")
                )

        if disc_rate > 0.15:
            findings.append(f"Discontinuation rate: {disc_rate:.0%}")
            recommendations.append(
                "High discontinuation rate — assess tolerability "
                "and consider dose modification strategy"
            )

        # --- Design-specific recommendations ---
        if design == "sample_size_reestimation":
            if effect_size is not None and effect_size > 0:
                # Conditional power calculation (simplified)
                cp = _clamp(1.0 - p_value if p_value else 0.5)
                if cp < 0.5 and info_fraction < 0.75:
                    new_n = int(target * (0.8 / max(cp, 0.1)))
                    recommendations.append(
                        f"Sample size re-estimation: increase N from {target} "
                        f"to {min(new_n, target * 3)} based on conditional power {cp:.0%}"
                    )
                    findings.append(f"Conditional power: {cp:.0%}")

        elif design == "response_adaptive":
            rates = interim.get("response_rates", {})
            if rates:
                best_arm = max(rates, key=rates.get)
                findings.append(
                    f"Response rates by arm: {rates}"
                )
                recommendations.append(
                    f"Shift allocation ratio to favour '{best_arm}' "
                    f"arm (response rate: {rates[best_arm]:.0%})"
                )

        elif design == "biomarker_adaptive":
            recommendations.append(
                "Assess biomarker-defined subgroup treatment effects "
                "for potential enrichment or stratification refinement"
            )

        elif design == "platform":
            recommendations.append(
                "Evaluate graduated arms for futility; "
                "assess new arm additions based on emerging data"
            )

        elif design == "seamless":
            if info_fraction >= 0.5:
                recommendations.append(
                    "Phase II/III seamless transition point reached — "
                    "confirm endpoint alignment and regulatory acceptance"
                )

        # General recommendations
        if not recommendations:
            recommendations.append(
                "Continue as planned — no adaptation triggers at this interim analysis"
            )

        severity = SeverityLevel.INFORMATIONAL
        if any("FUTILITY" in f for f in findings):
            severity = SeverityLevel.HIGH
        elif any("EFFICACY" in f for f in findings):
            severity = SeverityLevel.MODERATE
        elif sae_count > 0 and sae_count / max(enrolled, 1) > 0.10:
            severity = SeverityLevel.HIGH

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "FDA Guidance: Adaptive Designs for Clinical Trials (2019)",
                "EMA Reflection Paper on Methodological Issues (2007, updated 2019)",
                "ICH E9(R1) Estimands Framework",
                "Mehta & Pocock: Adaptive Increase in Sample Size (NEJM 2011)",
            ],
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=_clamp(0.4 + info_fraction * 0.4),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 6 — Safety Signal Detection
# ═══════════════════════════════════════════════════════════════════════════


class SafetySignalWorkflow(BaseTrialWorkflow):
    """Safety signal detection workflow using disproportionality metrics,
    WHO-UMC causality assessment, and pharmacogenomic risk factors.

    Inputs
    ------
    adverse_event : str
        MedDRA preferred term for the adverse event.
    drug : str
        Drug name or identifier.
    event_count_drug : int
        Number of events in the drug arm.
    total_drug : int
        Total patients in the drug arm.
    event_count_comparator : int
        Number of events in the comparator arm.
    total_comparator : int
        Total patients in the comparator arm.
    patient_profile : dict
        Patient demographics and PGx data.
    dechallenge_positive : bool
        Whether dechallenge was positive (event resolved on drug withdrawal).
    rechallenge_positive : bool
        Whether rechallenge was positive (event recurred on re-exposure).
    time_to_onset : str
        Temporal relationship (e.g. 'days', 'weeks', 'months').
    literature_reports : int
        Number of published reports of this AE with this drug.
    """

    workflow_type = TrialWorkflowType.SAFETY_SIGNAL

    # PGx risk alleles associated with increased AE risk
    _PGX_RISK_FACTORS: Dict[str, Dict[str, str]] = {
        "CYP2D6": {
            "poor_metabolizer": "Increased drug exposure — higher AE risk",
            "ultra_rapid": "Reduced efficacy, potential toxicity from active metabolites",
        },
        "CYP2C19": {
            "poor_metabolizer": "Altered metabolism — dose adjustment needed",
        },
        "CYP3A4": {
            "poor_metabolizer": "Increased systemic exposure",
        },
        "HLA-B*57:01": {
            "positive": "Abacavir hypersensitivity risk",
        },
        "HLA-B*15:02": {
            "positive": "Carbamazepine SJS/TEN risk (Asian populations)",
        },
        "DPYD": {
            "deficient": "Fluoropyrimidine severe toxicity risk",
        },
        "UGT1A1": {
            "poor_metabolizer": "Irinotecan toxicity risk (*28/*28)",
        },
        "TPMT": {
            "deficient": "Thiopurine myelosuppression risk",
        },
    }

    # WHO-UMC causality categories and scoring
    _CAUSALITY_CRITERIA = {
        "certain": {"min_score": 8, "description": "Definitive causal relationship"},
        "probable": {"min_score": 5, "description": "Likely causal relationship"},
        "possible": {"min_score": 3, "description": "Cannot rule out causal relationship"},
        "unlikely": {"min_score": 1, "description": "Doubtful causal relationship"},
        "unassessable": {"min_score": 0, "description": "Insufficient information"},
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        if not inp.get("adverse_event"):
            warnings.append("No adverse event specified — using 'unspecified AE'")
            inp["adverse_event"] = "unspecified AE"
        if not inp.get("drug"):
            warnings.append("No drug specified")
            inp["drug"] = "unspecified"

        inp.setdefault("event_count_drug", 0)
        inp.setdefault("total_drug", 1)
        inp.setdefault("event_count_comparator", 0)
        inp.setdefault("total_comparator", 1)
        inp.setdefault("patient_profile", {})
        inp.setdefault("dechallenge_positive", False)
        inp.setdefault("rechallenge_positive", False)
        inp.setdefault("time_to_onset", "unknown")
        inp.setdefault("literature_reports", 0)

        # Prevent division by zero
        for key in ("total_drug", "total_comparator"):
            try:
                v = int(inp[key])
                if v <= 0:
                    inp[key] = 1
                else:
                    inp[key] = v
            except (TypeError, ValueError):
                inp[key] = 1

        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []
        triggers: List[str] = []

        ae = inputs["adverse_event"]
        drug = inputs["drug"]
        a = int(inputs["event_count_drug"])     # AE + drug
        b = int(inputs["total_drug"]) - a       # no AE + drug
        c = int(inputs["event_count_comparator"])  # AE + comparator
        d = int(inputs["total_comparator"]) - c   # no AE + comparator
        b = max(b, 0)
        d = max(d, 0)

        findings.append(
            f"Safety signal analysis: '{ae}' with '{drug}'"
        )

        # --- Disproportionality metrics ---
        # PRR = (a/(a+b)) / (c/(c+d))
        rate_drug = a / max(a + b, 1)
        rate_comp = c / max(c + d, 1)
        prr = rate_drug / max(rate_comp, 0.001)
        findings.append(
            f"Drug arm AE rate: {rate_drug:.1%} ({a}/{a + b})"
        )
        findings.append(
            f"Comparator arm AE rate: {rate_comp:.1%} ({c}/{c + d})"
        )
        findings.append(f"Proportional Reporting Ratio (PRR): {prr:.2f}")

        # ROR = (a*d) / (b*c)
        ror = (a * d) / max(b * c, 1)
        findings.append(f"Reporting Odds Ratio (ROR): {ror:.2f}")

        # Signal threshold (PRR >= 2 AND N >= 3)
        signal_detected = prr >= 2.0 and a >= 3
        if signal_detected:
            findings.append(
                "SIGNAL DETECTED: PRR >= 2.0 with >= 3 cases"
            )

        # --- WHO-UMC Causality Assessment ---
        causality_score = 0

        # Temporal relationship
        onset = inputs["time_to_onset"].lower()
        if onset in ("days", "hours", "immediate"):
            causality_score += 2
            findings.append("Temporal relationship: compatible (days/hours)")
        elif onset in ("weeks",):
            causality_score += 1
            findings.append("Temporal relationship: plausible (weeks)")
        else:
            findings.append(f"Temporal relationship: {onset}")

        # Dechallenge
        if inputs["dechallenge_positive"]:
            causality_score += 2
            findings.append("Positive dechallenge: event resolved on withdrawal")

        # Rechallenge
        if inputs["rechallenge_positive"]:
            causality_score += 3
            findings.append("Positive rechallenge: event recurred on re-exposure")

        # Literature support
        lit = inputs["literature_reports"]
        if lit > 10:
            causality_score += 2
            findings.append(f"Strong literature support: {lit} published reports")
        elif lit > 0:
            causality_score += 1
            findings.append(f"Some literature support: {lit} published reports")

        # Dose-response (if PRR > 3, suggestive)
        if prr > 3.0:
            causality_score += 1

        # Determine causality category
        causality = "unassessable"
        for cat in ["certain", "probable", "possible", "unlikely"]:
            if causality_score >= self._CAUSALITY_CRITERIA[cat]["min_score"]:
                causality = cat
                break

        findings.append(
            f"WHO-UMC causality assessment: {causality.upper()} "
            f"(score: {causality_score})"
        )

        # --- PGx risk factor check ---
        patient = inputs["patient_profile"]
        pgx_variants = patient.get("genomic_variants", [])
        pgx_risks: List[str] = []
        for gene, alleles in self._PGX_RISK_FACTORS.items():
            for allele, risk_desc in alleles.items():
                for variant in pgx_variants:
                    if gene.lower() in variant.lower() and allele.lower() in variant.lower():
                        pgx_risks.append(f"{gene} {allele}: {risk_desc}")

        if pgx_risks:
            findings.append("PGx risk factors identified:")
            for risk in pgx_risks:
                findings.append(f"  - {risk}")
            triggers.append(
                _trigger_string("PGx", "Pharmacogenomic risk factors detected")
            )

        # --- Recommendations ---
        if causality in ("certain", "probable"):
            recommendations.append(
                f"Immediate DSMB notification for '{ae}' — "
                f"causality assessed as {causality}"
            )
            recommendations.append(
                "Update Investigator Brochure and informed consent"
            )
        elif causality == "possible":
            recommendations.append(
                f"Enhanced monitoring for '{ae}' — add to AE of special interest list"
            )
        if signal_detected:
            recommendations.append(
                "Signal detected — conduct formal benefit-risk assessment"
            )
        if prr > 5.0:
            recommendations.append(
                f"Strong signal (PRR={prr:.1f}) — consider regulatory notification "
                f"(IND Safety Report / SUSAR)"
            )
        if pgx_risks:
            recommendations.append(
                "Consider pharmacogenomic screening for at-risk populations"
            )

        # Severity
        severity = SeverityLevel.INFORMATIONAL
        if causality in ("certain", "probable"):
            severity = SeverityLevel.CRITICAL
        elif signal_detected:
            severity = SeverityLevel.HIGH
        elif causality == "possible":
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "WHO-UMC Causality Assessment System",
                "ICH E2B(R3) Individual Case Safety Report",
                "FDA Guidance: Safety Reporting Requirements (2015)",
                "Evans & Waller: Use of PRR for Signal Detection (Pharmacoepidemiology 2001)",
                "CIOMS Working Group Reports",
            ],
            severity=severity,
            cross_agent_triggers=triggers,
            confidence=_clamp(0.3 + causality_score * 0.08),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 7 — Regulatory Document Generation
# ═══════════════════════════════════════════════════════════════════════════


class RegulatoryDocumentWorkflow(BaseTrialWorkflow):
    """Regulatory document framework generation workflow.

    Inputs
    ------
    document_type : str
        One of: IND, CSR, Briefing, PSP, RMP, DSUR.
    trial_data : dict
        Trial metadata (title, phase, indication, sponsor, nct_id,
        enrollment, primary_endpoint, results_summary).
    regulatory_agency : str
        Target agency (fda, ema, pmda, etc.).
    """

    workflow_type = TrialWorkflowType.REGULATORY_DOCS

    # Document section templates by type
    _DOCUMENT_TEMPLATES: Dict[str, Dict[str, List[str]]] = {
        "ind": {
            "title": "Investigational New Drug Application",
            "sections": [
                "1. Cover Sheet (Form FDA 1571)",
                "2. Table of Contents",
                "3. Introductory Statement and General Investigational Plan",
                "4. Investigator's Brochure",
                "5. Clinical Protocol(s)",
                "  5.1 Study Objectives",
                "  5.2 Study Design and Methods",
                "  5.3 Selection of Study Population",
                "  5.4 Dosing Regimen",
                "  5.5 Efficacy Assessments",
                "  5.6 Safety Assessments",
                "  5.7 Statistical Analysis Plan",
                "6. Chemistry, Manufacturing, and Controls (CMC)",
                "7. Pharmacology and Toxicology Data",
                "8. Previous Human Experience",
                "9. Additional Information",
            ],
        },
        "csr": {
            "title": "Clinical Study Report",
            "sections": [
                "1. Title Page",
                "2. Synopsis",
                "3. Table of Contents and List of Abbreviations",
                "4. Ethics — IRB/IEC Approval, Informed Consent",
                "5. Investigators and Study Administrative Structure",
                "6. Introduction",
                "7. Study Objectives",
                "8. Investigational Plan",
                "  8.1 Overall Study Design",
                "  8.2 Discussion of Study Design",
                "  8.3 Selection of Study Population",
                "  8.4 Treatments Administered",
                "  8.5 Efficacy and Safety Variables",
                "  8.6 Data Quality Assurance",
                "  8.7 Statistical Methods",
                "9. Study Patients",
                "10. Efficacy Evaluation",
                "11. Safety Evaluation",
                "12. Discussion and Overall Conclusions",
                "13. Tables, Figures, and Graphs",
                "14. Reference List",
                "15. Appendices",
            ],
        },
        "briefing": {
            "title": "Briefing Document for Regulatory Meeting",
            "sections": [
                "1. Executive Summary",
                "2. Product Overview and Development Rationale",
                "3. Nonclinical Summary",
                "4. Clinical Development Program Overview",
                "5. Key Efficacy Results",
                "6. Key Safety Results",
                "7. Benefit-Risk Assessment",
                "8. Questions for Discussion",
                "9. Proposed Regulatory Strategy",
                "10. Appendices — Supporting Data",
            ],
        },
        "psp": {
            "title": "Pediatric Study Plan",
            "sections": [
                "1. Executive Summary",
                "2. Product Information",
                "3. Condition or Indication",
                "4. Analysis of Condition in Pediatric Population",
                "5. Existing Data Relevant to Pediatric Populations",
                "6. Planned Pediatric Development Strategy",
                "  6.1 Study Design(s)",
                "  6.2 Age Groups and Extrapolation",
                "  6.3 Formulation Development",
                "  6.4 Timeline and Milestones",
                "7. Request for Deferral or Waiver (if applicable)",
                "8. References",
            ],
        },
        "rmp": {
            "title": "Risk Management Plan",
            "sections": [
                "1. Product Overview",
                "2. Safety Specification",
                "  2.1 Epidemiology of the Indication",
                "  2.2 Non-clinical Part of the Safety Specification",
                "  2.3 Clinical Trial Exposure",
                "  2.4 Populations Not Studied in Clinical Trials",
                "  2.5 Post-authorisation Experience",
                "  2.6 Identified Risks",
                "  2.7 Potential Risks",
                "  2.8 Missing Information",
                "3. Pharmacovigilance Plan",
                "  3.1 Routine Pharmacovigilance",
                "  3.2 Additional Pharmacovigilance Activities",
                "4. Risk Minimisation Measures",
                "5. Summary of the Risk Management Plan",
            ],
        },
        "dsur": {
            "title": "Development Safety Update Report",
            "sections": [
                "1. Executive Summary",
                "2. Introduction",
                "3. Worldwide Marketing Approval Status",
                "4. Update on Actions Taken for Safety Reasons",
                "5. Changes to Reference Safety Information",
                "6. Inventory of Clinical Trials",
                "7. Estimated Cumulative Exposure",
                "8. Presentation of Data in Line Listings and Summary Tables",
                "9. Significant Findings from Clinical Trials",
                "10. Safety Findings from Non-Interventional Studies",
                "11. Other Safety Information",
                "12. Late-Breaking Information",
                "13. Overall Safety Assessment",
                "14. Summary of Important Risks",
                "15. Benefit-Risk Analysis",
                "16. Conclusions",
                "17. Appendices",
            ],
        },
    }

    _AGENCY_REQUIREMENTS: Dict[str, Dict[str, str]] = {
        "fda": {
            "format": "eCTD v4.0",
            "regulation": "21 CFR Part 312 (IND), 21 CFR Part 314 (NDA)",
            "timeline": "30-day IND review; 12-month NDA PDUFA",
        },
        "ema": {
            "format": "eCTD v4.0",
            "regulation": "Regulation (EC) No 726/2004, Directive 2001/83/EC",
            "timeline": "210-day centralised procedure",
        },
        "pmda": {
            "format": "eCTD v4.0 (Japan Module)",
            "regulation": "Pharmaceutical and Medical Device Act (PMD Act)",
            "timeline": "12-month standard review",
        },
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        doc_type = str(inp.get("document_type", "csr")).lower()
        valid_types = set(self._DOCUMENT_TEMPLATES.keys())
        if doc_type not in valid_types:
            warnings.append(
                f"Document type '{doc_type}' not recognised — "
                f"defaulting to CSR. Valid: {', '.join(sorted(valid_types))}"
            )
            doc_type = "csr"
        inp["document_type"] = doc_type

        inp.setdefault("trial_data", {})
        agency = str(inp.get("regulatory_agency", "fda")).lower()
        if agency not in self._AGENCY_REQUIREMENTS:
            warnings.append(f"Agency '{agency}' not in database — defaulting to FDA")
            agency = "fda"
        inp["regulatory_agency"] = agency

        return inp

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        doc_type = inputs["document_type"]
        trial = inputs["trial_data"]
        agency = inputs["regulatory_agency"]

        template = self._DOCUMENT_TEMPLATES[doc_type]
        agency_info = self._AGENCY_REQUIREMENTS.get(agency, {})

        findings.append(
            f"Document: {template['title']} ({doc_type.upper()})"
        )
        findings.append(
            f"Target agency: {agency.upper()} — "
            f"format: {agency_info.get('format', 'eCTD')}"
        )
        findings.append(
            f"Regulatory framework: {agency_info.get('regulation', 'N/A')}"
        )
        findings.append(
            f"Expected timeline: {agency_info.get('timeline', 'N/A')}"
        )

        # Populate sections with trial data
        findings.append(f"\nDocument outline ({len(template['sections'])} sections):")
        for section in template["sections"]:
            findings.append(f"  {section}")

        # Trial-specific population
        if trial:
            findings.append("\nTrial context:")
            if trial.get("title"):
                findings.append(f"  Title: {trial['title']}")
            if trial.get("nct_id"):
                findings.append(f"  NCT: {trial['nct_id']}")
            if trial.get("phase"):
                findings.append(f"  Phase: {trial['phase']}")
            if trial.get("indication"):
                findings.append(f"  Indication: {trial['indication']}")
            if trial.get("enrollment"):
                findings.append(f"  Enrollment: {trial['enrollment']}")

        # Recommendations
        recommendations.append(
            f"Prepare {doc_type.upper()} following {agency_info.get('format', 'eCTD')} "
            f"format requirements"
        )
        if doc_type == "ind":
            recommendations.append(
                "Schedule pre-IND meeting to align on CMC and nonclinical requirements"
            )
        elif doc_type == "csr":
            recommendations.append(
                "Follow ICH E3 structure; ensure all protocol amendments are documented"
            )
        elif doc_type == "rmp":
            recommendations.append(
                "Include all identified and potential risks from the safety database"
            )

        references = [
            "ICH M4 Common Technical Document (eCTD)",
        ]
        if doc_type == "csr":
            references.append("ICH E3 Structure and Content of Clinical Study Reports")
        elif doc_type == "ind":
            references.append("21 CFR Part 312 — Investigational New Drug Application")
        elif doc_type == "dsur":
            references.append("ICH E2F Development Safety Update Report")
        elif doc_type == "rmp":
            references.append("EMA Guideline on Good Pharmacovigilance Practices Module V")

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=references,
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.75,
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 8 — Competitive Intelligence
# ═══════════════════════════════════════════════════════════════════════════


class CompetitiveIntelligenceWorkflow(BaseTrialWorkflow):
    """Competitive landscape analysis workflow that assesses competing
    trials, enrollment progress, design differentiation, and threat scores.

    Inputs
    ------
    therapeutic_area : str
        Therapeutic area to search.
    mechanism_of_action : str
        Drug mechanism of action.
    indication : str
        Target indication.
    competitors : list[dict]
        Known competitor trials with fields: trial_id, sponsor, phase,
        indication, mechanism, enrollment_target, enrollment_actual,
        estimated_completion, design_type, primary_endpoint.
    own_trial : dict
        Own trial details for differentiation comparison.
    """

    workflow_type = TrialWorkflowType.COMPETITIVE_INTEL

    # Threat weights
    _THREAT_WEIGHTS = {
        "phase_advancement": 0.30,
        "enrollment_progress": 0.25,
        "sponsor_resources": 0.20,
        "differentiation": 0.25,
    }

    # Sponsor tier classification
    _SPONSOR_TIERS: Dict[str, float] = {
        "large_pharma": 0.9, "mid_pharma": 0.7, "biotech": 0.5,
        "academic": 0.3, "unknown": 0.4,
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        inp.setdefault("therapeutic_area", "unspecified")
        inp.setdefault("mechanism_of_action", "unspecified")
        inp.setdefault("indication", "unspecified")
        inp.setdefault("competitors", [])
        inp.setdefault("own_trial", {})

        if not inp["competitors"]:
            warnings.append(
                "No competitor data provided — analysis will be limited"
            )

        return inp

    def _calculate_threat(self, competitor: dict, own_trial: dict) -> float:
        """Calculate composite threat score for a competitor."""
        # Phase advancement (higher phase = greater threat)
        phase_scores = {
            "phase_i": 0.2, "phase_i_ii": 0.3, "phase_ii": 0.4,
            "phase_ii_iii": 0.6, "phase_iii": 0.8, "phase_iv": 0.95,
        }
        phase = competitor.get("phase", "phase_ii").lower().replace(" ", "_")
        phase_score = phase_scores.get(phase, 0.4)

        # Enrollment progress
        target = competitor.get("enrollment_target", 100)
        actual = competitor.get("enrollment_actual", 0)
        enrollment_score = _clamp(actual / max(target, 1))

        # Sponsor resources
        sponsor_type = competitor.get("sponsor_type", "unknown")
        sponsor_score = self._SPONSOR_TIERS.get(sponsor_type, 0.4)

        # Differentiation (lower = more threatening, i.e. more similar)
        own_moa = own_trial.get("mechanism", "").lower()
        comp_moa = competitor.get("mechanism", "").lower()
        if own_moa and comp_moa:
            # Simple word overlap similarity
            own_words = set(own_moa.split())
            comp_words = set(comp_moa.split())
            overlap = len(own_words & comp_words) / max(len(own_words | comp_words), 1)
            differentiation_score = overlap  # Higher overlap = more threat
        else:
            differentiation_score = 0.5

        threat = (
            self._THREAT_WEIGHTS["phase_advancement"] * phase_score
            + self._THREAT_WEIGHTS["enrollment_progress"] * enrollment_score
            + self._THREAT_WEIGHTS["sponsor_resources"] * sponsor_score
            + self._THREAT_WEIGHTS["differentiation"] * differentiation_score
        )
        return _clamp(threat)

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        ta = inputs["therapeutic_area"]
        moa = inputs["mechanism_of_action"]
        indication = inputs["indication"]
        competitors = inputs["competitors"]
        own_trial = inputs["own_trial"]

        findings.append(
            f"Competitive landscape: {indication} / {ta} / {moa}"
        )

        if not competitors:
            findings.append("No competitor data available")
            recommendations.append(
                "Perform ClinicalTrials.gov search for competing trials "
                "in the target indication"
            )
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=findings,
                recommendations=recommendations,
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.3,
            )

        scored: List[dict] = []
        for comp in competitors:
            threat = self._calculate_threat(comp, own_trial)
            scored.append({**comp, "_threat": threat})

        scored.sort(key=lambda c: c["_threat"], reverse=True)

        findings.append(f"Analysed {len(scored)} competing trials:")
        high_threats = 0
        for i, comp in enumerate(scored[:10], 1):
            threat = comp["_threat"]
            level = "HIGH" if threat > 0.7 else "MODERATE" if threat > 0.4 else "LOW"
            if threat > 0.7:
                high_threats += 1
            findings.append(
                f"  #{i}: {comp.get('trial_id', 'N/A')} "
                f"({comp.get('sponsor', 'Unknown')}) — "
                f"Phase {comp.get('phase', 'N/A')}, "
                f"threat: {threat:.2f} [{level}]"
            )

        # Landscape summary
        phases = [c.get("phase", "") for c in scored]
        phase_counts = {}
        for p in phases:
            phase_counts[p] = phase_counts.get(p, 0) + 1
        findings.append(
            f"Phase distribution: {phase_counts}"
        )

        # Recommendations
        if high_threats > 0:
            recommendations.append(
                f"{high_threats} high-threat competitor(s) identified — "
                f"prioritise differentiation strategy"
            )
        if high_threats > 3:
            recommendations.append(
                "Crowded competitive landscape — consider niche "
                "subpopulation or combination strategy"
            )
        recommendations.append(
            "Monitor competitor enrollment rates and regulatory milestones quarterly"
        )

        severity = SeverityLevel.INFORMATIONAL
        if high_threats > 3:
            severity = SeverityLevel.HIGH
        elif high_threats > 0:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "ClinicalTrials.gov Data Standards",
                "Tufts CSDD Competitive Intelligence Benchmarks",
            ],
            severity=severity,
            confidence=_clamp(0.4 + len(competitors) * 0.03),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 9 — Diversity Assessment
# ═══════════════════════════════════════════════════════════════════════════


class DiversityAssessmentWorkflow(BaseTrialWorkflow):
    """Diversity and representativeness assessment workflow that calculates
    Population Representativeness Index (PRI), geographic diversity, and
    genetic ancestry coverage to support FDA Diversity Action Plans.

    Inputs
    ------
    enrollment_data : dict
        Keys: total_enrolled, demographics (dict of race/ethnicity counts),
        sex_distribution (dict), age_distribution (dict of age band counts).
    site_demographics : list[dict]
        Per-site demographic data.
    indication_epidemiology : dict
        Disease prevalence by demographic group (race_prevalence,
        sex_prevalence, age_prevalence).
    cyp450_coverage : dict
        CYP450 polymorphism coverage by genetic ancestry group.
    """

    workflow_type = TrialWorkflowType.DIVERSITY_ASSESSMENT

    # US Census 2020 approximate proportions for reference
    _US_CENSUS: Dict[str, float] = {
        "white": 0.576, "black": 0.134, "hispanic": 0.189,
        "asian": 0.060, "native_american": 0.013,
        "pacific_islander": 0.003, "multiracial": 0.025,
    }

    # Known CYP450 polymorphism frequencies by ancestry
    _CYP_ANCESTRY: Dict[str, Dict[str, float]] = {
        "CYP2D6_poor_metabolizer": {
            "european": 0.07, "african": 0.03, "east_asian": 0.01,
            "south_asian": 0.04, "hispanic": 0.05,
        },
        "CYP2C19_poor_metabolizer": {
            "european": 0.02, "african": 0.04, "east_asian": 0.15,
            "south_asian": 0.12, "hispanic": 0.03,
        },
        "CYP3A5_expresser": {
            "european": 0.15, "african": 0.75, "east_asian": 0.30,
            "south_asian": 0.40, "hispanic": 0.25,
        },
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        warnings = self._init_warnings(inp)

        inp.setdefault("enrollment_data", {})
        inp.setdefault("site_demographics", [])
        inp.setdefault("indication_epidemiology", {})
        inp.setdefault("cyp450_coverage", {})

        enrollment = inp["enrollment_data"]
        enrollment.setdefault("total_enrolled", 0)
        enrollment.setdefault("demographics", {})
        enrollment.setdefault("sex_distribution", {})
        enrollment.setdefault("age_distribution", {})

        if enrollment["total_enrolled"] == 0:
            warnings.append(
                "No enrollment data — diversity metrics cannot be calculated"
            )

        return inp

    def _calculate_pri(
        self, enrollment_pct: Dict[str, float],
        population_pct: Dict[str, float],
    ) -> float:
        """Calculate Population Representativeness Index.

        PRI = 1 - (sum |enrolled% - population%|) / 2
        Range: 0 (no representation) to 1 (perfect representation).
        """
        all_groups = set(enrollment_pct.keys()) | set(population_pct.keys())
        total_diff = sum(
            abs(enrollment_pct.get(g, 0) - population_pct.get(g, 0))
            for g in all_groups
        )
        return _clamp(1.0 - total_diff / 2.0)

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        enrollment = inputs["enrollment_data"]
        epi = inputs["indication_epidemiology"]
        total = enrollment["total_enrolled"]
        demographics = enrollment["demographics"]

        findings.append(f"Diversity assessment — {total} participants enrolled")

        if total == 0:
            return WorkflowResult(
                workflow_type=self.workflow_type,
                findings=["No enrollment data available for diversity analysis"],
                recommendations=[
                    "Establish diversity enrollment targets aligned with "
                    "FDA Diversity Action Plan requirements"
                ],
                severity=SeverityLevel.INFORMATIONAL,
                confidence=0.0,
            )

        # --- Race/ethnicity PRI ---
        enrolled_pct = {
            group: count / max(total, 1)
            for group, count in demographics.items()
        }

        # Compare against disease epidemiology if available, else census
        reference_pct = epi.get("race_prevalence", self._US_CENSUS)
        pri = self._calculate_pri(enrolled_pct, reference_pct)
        findings.append(
            f"Population Representativeness Index (PRI): {pri:.2f}"
        )

        if pri < 0.5:
            findings.append(
                "UNDERREPRESENTATION DETECTED: PRI < 0.50 indicates "
                "significant demographic disparity"
            )
        elif pri < 0.7:
            findings.append(
                "Moderate representation gap: PRI between 0.50-0.70"
            )
        else:
            findings.append("Good demographic representation: PRI >= 0.70")

        # Per-group analysis
        findings.append("\nDemographic representation vs reference:")
        for group in sorted(set(enrolled_pct.keys()) | set(reference_pct.keys())):
            e_pct = enrolled_pct.get(group, 0)
            r_pct = reference_pct.get(group, 0)
            gap = e_pct - r_pct
            status = "OVER" if gap > 0.05 else "UNDER" if gap < -0.05 else "OK"
            findings.append(
                f"  {group}: enrolled {e_pct:.1%} vs reference {r_pct:.1%} "
                f"[{status}]"
            )

        # --- Sex distribution ---
        sex_dist = enrollment.get("sex_distribution", {})
        if sex_dist:
            findings.append(f"\nSex distribution: {sex_dist}")
            sex_ref = epi.get("sex_prevalence", {"male": 0.5, "female": 0.5})
            sex_pri = self._calculate_pri(
                {k: v / max(total, 1) if isinstance(v, (int, float)) and v > 1 else v
                 for k, v in sex_dist.items()},
                sex_ref,
            )
            findings.append(f"Sex PRI: {sex_pri:.2f}")

        # --- Geographic diversity ---
        sites = inputs.get("site_demographics", [])
        if sites:
            countries = set(s.get("country", "") for s in sites)
            findings.append(
                f"\nGeographic coverage: {len(sites)} sites across "
                f"{len(countries)} countries"
            )
            if len(countries) < 3:
                recommendations.append(
                    "Limited geographic diversity — consider adding sites "
                    "in underrepresented regions"
                )

        # --- CYP450 genetic ancestry coverage ---
        cyp_coverage = inputs.get("cyp450_coverage", {})
        if cyp_coverage:
            findings.append("\nCYP450 polymorphism coverage by ancestry:")
            for gene, ancestries in cyp_coverage.items():
                findings.append(f"  {gene}: {ancestries}")
        else:
            # Use reference data for recommendations
            findings.append(
                "\nCYP450 polymorphism consideration: frequency variation "
                "by ancestry affects drug metabolism"
            )
            for gene, freqs in self._CYP_ANCESTRY.items():
                max_pop = max(freqs, key=freqs.get)
                findings.append(
                    f"  {gene}: highest in {max_pop} ({freqs[max_pop]:.0%}) — "
                    f"ensure adequate representation"
                )

        # --- FDA Diversity Action Plan recommendations ---
        recommendations.append(
            "FDA Diversity Action Plan components:"
        )
        if pri < 0.7:
            underrepresented = [
                g for g in reference_pct
                if enrolled_pct.get(g, 0) < reference_pct[g] * 0.5
            ]
            if underrepresented:
                recommendations.append(
                    f"  1. Enrich enrollment of: {', '.join(underrepresented)}"
                )
            recommendations.append(
                "  2. Expand site network to areas with higher minority populations"
            )
            recommendations.append(
                "  3. Implement community engagement and outreach programs"
            )
            recommendations.append(
                "  4. Consider decentralised trial elements to improve access"
            )
        else:
            recommendations.append(
                "  Current enrollment meets diversity benchmarks — "
                "maintain monitoring and enrollment targets"
            )

        severity = SeverityLevel.INFORMATIONAL
        if pri < 0.5:
            severity = SeverityLevel.HIGH
        elif pri < 0.7:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "FDA Guidance: Diversity Plans to Improve Enrollment (2024)",
                "FDA Guidance: Enhancing the Diversity of Clinical Trial Populations (2020)",
                "NIH Policy on Inclusion of Women and Minorities",
                "ICH E5(R1) Ethnic Factors in Acceptability of Foreign Clinical Data",
            ],
            severity=severity,
            confidence=_clamp(0.5 + pri * 0.3),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW 10 — Decentralized Trial Planning
# ═══════════════════════════════════════════════════════════════════════════


class DecentralizedPlanningWorkflow(BaseTrialWorkflow):
    """Decentralized clinical trial (DCT) feasibility assessment and
    implementation planning workflow.

    Inputs
    ------
    protocol : dict
        Protocol details (phase, indication, visit_schedule, procedures,
        endpoints, duration_months).
    dct_requirements : dict
        DCT component requirements and preferences. Keys: components
        (list of DCTComponent strings), patient_population (str),
        geographic_scope (str), budget_constraints (str).
    regulatory_context : str
        Target regulatory agency context.
    """

    workflow_type = TrialWorkflowType.DECENTRALIZED_PLANNING

    # DCT component assessment dimensions
    _ASSESSMENT_DIMENSIONS = [
        "clinical_necessity",
        "regulatory_acceptance",
        "patient_preference",
        "data_quality",
        "technology_readiness",
    ]

    # Component-level feasibility baselines
    _COMPONENT_BASELINES: Dict[str, Dict[str, float]] = {
        "econsent": {
            "clinical_necessity": 0.6,
            "regulatory_acceptance": 0.85,
            "patient_preference": 0.80,
            "data_quality": 0.90,
            "technology_readiness": 0.90,
        },
        "telemedicine": {
            "clinical_necessity": 0.70,
            "regulatory_acceptance": 0.75,
            "patient_preference": 0.85,
            "data_quality": 0.65,
            "technology_readiness": 0.85,
        },
        "home_health": {
            "clinical_necessity": 0.65,
            "regulatory_acceptance": 0.70,
            "patient_preference": 0.90,
            "data_quality": 0.60,
            "technology_readiness": 0.70,
        },
        "local_labs": {
            "clinical_necessity": 0.75,
            "regulatory_acceptance": 0.80,
            "patient_preference": 0.85,
            "data_quality": 0.75,
            "technology_readiness": 0.80,
        },
        "wearables": {
            "clinical_necessity": 0.60,
            "regulatory_acceptance": 0.55,
            "patient_preference": 0.70,
            "data_quality": 0.55,
            "technology_readiness": 0.65,
        },
        "epro_ecoa": {
            "clinical_necessity": 0.70,
            "regulatory_acceptance": 0.80,
            "patient_preference": 0.75,
            "data_quality": 0.80,
            "technology_readiness": 0.85,
        },
        "direct_to_patient": {
            "clinical_necessity": 0.55,
            "regulatory_acceptance": 0.65,
            "patient_preference": 0.90,
            "data_quality": 0.70,
            "technology_readiness": 0.60,
        },
    }

    # Phase-specific adjustments
    _PHASE_ADJUSTMENTS: Dict[str, float] = {
        "phase_i": -0.15,  # More caution in early phase
        "phase_ii": -0.05,
        "phase_iii": 0.0,
        "phase_iv": 0.10,  # More flexibility post-approval
    }

    def preprocess(self, inputs: dict) -> dict:
        inp = dict(inputs)
        self._init_warnings(inp)

        inp.setdefault("protocol", {})
        inp.setdefault("dct_requirements", {})
        inp.setdefault("regulatory_context", "fda")

        protocol = inp["protocol"]
        protocol.setdefault("phase", "phase_iii")
        protocol.setdefault("indication", "unspecified")
        protocol.setdefault("visit_schedule", [])
        protocol.setdefault("procedures", [])
        protocol.setdefault("endpoints", [])
        protocol.setdefault("duration_months", 12)

        dct = inp["dct_requirements"]
        if not dct.get("components"):
            # Default: assess all components
            dct["components"] = [c.value for c in DCTComponent]
        dct.setdefault("patient_population", "general")
        dct.setdefault("geographic_scope", "us")

        return inp

    def _assess_component(
        self, component: str, phase: str, protocol: dict
    ) -> Dict[str, float]:
        """Assess a DCT component across all dimensions."""
        baselines = self._COMPONENT_BASELINES.get(
            component, {dim: 0.5 for dim in self._ASSESSMENT_DIMENSIONS}
        )
        phase_adj = self._PHASE_ADJUSTMENTS.get(phase, 0.0)

        scores: Dict[str, float] = {}
        for dim in self._ASSESSMENT_DIMENSIONS:
            base = baselines.get(dim, 0.5)
            adjusted = _clamp(base + phase_adj)

            # Procedure-specific adjustments
            procedures = protocol.get("procedures", [])
            if component == "home_health" and any(
                p in str(procedures).lower()
                for p in ["infusion", "injection", "blood draw"]
            ):
                if dim == "clinical_necessity":
                    adjusted = _clamp(adjusted + 0.15)

            if component == "wearables" and any(
                ep in str(protocol.get("endpoints", [])).lower()
                for ep in ["activity", "heart rate", "sleep", "steps"]
            ):
                if dim == "clinical_necessity":
                    adjusted = _clamp(adjusted + 0.20)
                if dim == "regulatory_acceptance":
                    adjusted = _clamp(adjusted + 0.10)

            scores[dim] = round(adjusted, 3)

        return scores

    def execute(self, inputs: dict) -> WorkflowResult:
        findings: List[str] = []
        recommendations: List[str] = []

        protocol = inputs["protocol"]
        dct = inputs["dct_requirements"]
        phase = protocol["phase"].lower().replace(" ", "_")
        components = dct["components"]

        findings.append(
            f"DCT feasibility assessment — {phase.replace('_', ' ').title()}, "
            f"{protocol['indication']}"
        )
        findings.append(f"Components assessed: {len(components)}")

        component_results: List[dict] = []
        for comp in components:
            comp_key = comp.lower().replace(" ", "_")
            dim_scores = self._assess_component(comp_key, phase, protocol)
            overall = sum(dim_scores.values()) / max(len(dim_scores), 1)
            component_results.append({
                "component": comp,
                "dimensions": dim_scores,
                "overall": round(overall, 3),
            })

        # Sort by overall feasibility
        component_results.sort(key=lambda c: c["overall"], reverse=True)

        findings.append("\nComponent feasibility scores:")
        feasible_count = 0
        for cr in component_results:
            status = (
                "RECOMMENDED" if cr["overall"] >= 0.70
                else "CONDITIONAL" if cr["overall"] >= 0.55
                else "NOT RECOMMENDED"
            )
            if cr["overall"] >= 0.70:
                feasible_count += 1
            findings.append(
                f"  {cr['component']}: {cr['overall']:.2f} [{status}]"
            )
            for dim, score in cr["dimensions"].items():
                findings.append(
                    f"    - {dim}: {score:.2f}"
                )

        # Overall DCT feasibility
        if component_results:
            avg_feasibility = sum(
                c["overall"] for c in component_results
            ) / len(component_results)
        else:
            avg_feasibility = 0.0

        findings.append(
            f"\nOverall DCT feasibility: {avg_feasibility:.2f}"
        )
        findings.append(
            f"Recommended components: {feasible_count}/{len(component_results)}"
        )

        # Implementation plan
        recommendations.append("DCT Implementation Plan:")
        recommended = [c for c in component_results if c["overall"] >= 0.70]
        conditional = [c for c in component_results if 0.55 <= c["overall"] < 0.70]
        not_rec = [c for c in component_results if c["overall"] < 0.55]

        if recommended:
            recommendations.append(
                f"  Phase 1 (deploy immediately): "
                f"{', '.join(c['component'] for c in recommended)}"
            )
        if conditional:
            recommendations.append(
                f"  Phase 2 (pilot with monitoring): "
                f"{', '.join(c['component'] for c in conditional)}"
            )
        if not_rec:
            recommendations.append(
                f"  Not recommended at this time: "
                f"{', '.join(c['component'] for c in not_rec)}"
            )

        # Regulatory considerations
        recommendations.append(
            "Regulatory considerations:"
        )
        recommendations.append(
            "  - Submit DCT protocol elements in pre-submission meeting package"
        )
        recommendations.append(
            "  - Ensure GCP compliance for all remote data collection"
        )
        if "wearables" in [c.lower() for c in components]:
            recommendations.append(
                "  - Wearable devices: validate against reference standards "
                "per FDA Digital Health guidance"
            )

        severity = SeverityLevel.INFORMATIONAL
        if avg_feasibility < 0.5:
            severity = SeverityLevel.HIGH
        elif avg_feasibility < 0.65:
            severity = SeverityLevel.MODERATE

        return WorkflowResult(
            workflow_type=self.workflow_type,
            findings=findings,
            recommendations=recommendations,
            guideline_references=[
                "FDA Guidance: Decentralized Clinical Trials for Drugs, Biologics, "
                "and Devices (2024)",
                "ICH E6(R3) GCP — Technology in Clinical Trials",
                "EMA Recommendation on DCT Elements (2022)",
                "CTTI Decentralized Clinical Trials Recommendations",
            ],
            severity=severity,
            confidence=_clamp(0.5 + avg_feasibility * 0.3),
        )


# ═══════════════════════════════════════════════════════════════════════════
# WORKFLOW ENGINE
# ═══════════════════════════════════════════════════════════════════════════


class WorkflowEngine:
    """Central dispatcher that maps TrialWorkflowType to the appropriate
    workflow implementation and handles query-based workflow detection."""

    _KEYWORD_MAP: Dict[str, TrialWorkflowType] = {
        # Protocol Design
        "protocol design": TrialWorkflowType.PROTOCOL_DESIGN,
        "study design": TrialWorkflowType.PROTOCOL_DESIGN,
        "trial design": TrialWorkflowType.PROTOCOL_DESIGN,
        "protocol": TrialWorkflowType.PROTOCOL_DESIGN,
        "sample size": TrialWorkflowType.PROTOCOL_DESIGN,
        "endpoints": TrialWorkflowType.PROTOCOL_DESIGN,
        "comparator": TrialWorkflowType.PROTOCOL_DESIGN,
        "randomization": TrialWorkflowType.PROTOCOL_DESIGN,
        # Patient Matching
        "patient matching": TrialWorkflowType.PATIENT_MATCHING,
        "eligibility": TrialWorkflowType.PATIENT_MATCHING,
        "patient match": TrialWorkflowType.PATIENT_MATCHING,
        "trial matching": TrialWorkflowType.PATIENT_MATCHING,
        "enroll patient": TrialWorkflowType.PATIENT_MATCHING,
        "inclusion criteria": TrialWorkflowType.PATIENT_MATCHING,
        "exclusion criteria": TrialWorkflowType.PATIENT_MATCHING,
        # Site Selection
        "site selection": TrialWorkflowType.SITE_SELECTION,
        "site ranking": TrialWorkflowType.SITE_SELECTION,
        "investigator": TrialWorkflowType.SITE_SELECTION,
        "site feasibility": TrialWorkflowType.SITE_SELECTION,
        "enrollment rate": TrialWorkflowType.SITE_SELECTION,
        # Eligibility Optimization
        "eligibility optimization": TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
        "broaden eligibility": TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
        "restrictive criteria": TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
        "population impact": TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
        # Adaptive Design
        "adaptive design": TrialWorkflowType.ADAPTIVE_DESIGN,
        "interim analysis": TrialWorkflowType.ADAPTIVE_DESIGN,
        "futility": TrialWorkflowType.ADAPTIVE_DESIGN,
        "sample size reestimation": TrialWorkflowType.ADAPTIVE_DESIGN,
        "group sequential": TrialWorkflowType.ADAPTIVE_DESIGN,
        "platform trial": TrialWorkflowType.ADAPTIVE_DESIGN,
        "basket trial": TrialWorkflowType.ADAPTIVE_DESIGN,
        "seamless design": TrialWorkflowType.ADAPTIVE_DESIGN,
        # Safety Signal
        "safety signal": TrialWorkflowType.SAFETY_SIGNAL,
        "adverse event": TrialWorkflowType.SAFETY_SIGNAL,
        "pharmacovigilance": TrialWorkflowType.SAFETY_SIGNAL,
        "dsmb": TrialWorkflowType.SAFETY_SIGNAL,
        "causality": TrialWorkflowType.SAFETY_SIGNAL,
        "disproportionality": TrialWorkflowType.SAFETY_SIGNAL,
        "sae": TrialWorkflowType.SAFETY_SIGNAL,
        "susar": TrialWorkflowType.SAFETY_SIGNAL,
        # Regulatory Documents
        "regulatory": TrialWorkflowType.REGULATORY_DOCS,
        "ind application": TrialWorkflowType.REGULATORY_DOCS,
        "clinical study report": TrialWorkflowType.REGULATORY_DOCS,
        "csr": TrialWorkflowType.REGULATORY_DOCS,
        "briefing document": TrialWorkflowType.REGULATORY_DOCS,
        "risk management plan": TrialWorkflowType.REGULATORY_DOCS,
        "dsur": TrialWorkflowType.REGULATORY_DOCS,
        # Competitive Intelligence
        "competitive": TrialWorkflowType.COMPETITIVE_INTEL,
        "competitor": TrialWorkflowType.COMPETITIVE_INTEL,
        "landscape": TrialWorkflowType.COMPETITIVE_INTEL,
        "competing trials": TrialWorkflowType.COMPETITIVE_INTEL,
        "threat assessment": TrialWorkflowType.COMPETITIVE_INTEL,
        # Diversity Assessment
        "diversity": TrialWorkflowType.DIVERSITY_ASSESSMENT,
        "representativeness": TrialWorkflowType.DIVERSITY_ASSESSMENT,
        "underrepresentation": TrialWorkflowType.DIVERSITY_ASSESSMENT,
        "minority enrollment": TrialWorkflowType.DIVERSITY_ASSESSMENT,
        "diversity action plan": TrialWorkflowType.DIVERSITY_ASSESSMENT,
        # Decentralized Planning
        "decentralized": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "dct": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "remote monitoring": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "telemedicine": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "wearable": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "econsent": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "epro": TrialWorkflowType.DECENTRALIZED_PLANNING,
        "home health": TrialWorkflowType.DECENTRALIZED_PLANNING,
    }

    def __init__(self) -> None:
        workflow_instances: List[BaseTrialWorkflow] = [
            ProtocolDesignWorkflow(),
            PatientMatchingWorkflow(),
            SiteSelectionWorkflow(),
            EligibilityOptimizationWorkflow(),
            AdaptiveDesignWorkflow(),
            SafetySignalWorkflow(),
            RegulatoryDocumentWorkflow(),
            CompetitiveIntelligenceWorkflow(),
            DiversityAssessmentWorkflow(),
            DecentralizedPlanningWorkflow(),
        ]
        self._workflows: Dict[TrialWorkflowType, BaseTrialWorkflow] = {
            wf.workflow_type: wf for wf in workflow_instances
        }

    @property
    def _WORKFLOWS(self) -> Dict[TrialWorkflowType, BaseTrialWorkflow]:
        """Public access to the workflows registry."""
        return self._workflows

    # ── public API ────────────────────────────────────────────────────

    def run_workflow(
        self, workflow_type: TrialWorkflowType, inputs: dict
    ) -> WorkflowResult:
        """Execute a specific workflow by type."""
        wf = self._workflows.get(workflow_type)
        if wf is None:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
        return wf.run(inputs)

    def detect_workflow(self, query: str) -> TrialWorkflowType:
        """Detect the most appropriate workflow from a free-text query."""
        query_lower = query.lower()
        # Score each workflow type by keyword matches
        scores: Dict[TrialWorkflowType, int] = {}
        for keyword, wf_type in self._KEYWORD_MAP.items():
            if keyword in query_lower:
                scores[wf_type] = scores.get(wf_type, 0) + len(keyword)

        if scores:
            return max(scores, key=scores.get)
        return TrialWorkflowType.GENERAL

    def list_workflows(self) -> List[str]:
        """Return list of registered workflow type values."""
        return [wt.value for wt in self._workflows]
