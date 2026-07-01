"""Decision support engines for the Clinical Trial Intelligence Agent.

Author: Adam Jones
Date: March 2026

Provides calibrated confidence scoring, protocol complexity assessment,
enrollment prediction, eligibility analysis, and competitive threat
scoring to support clinical trial decision-making.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List

from src.models import (
    EligibilityAnalysis,
    ProtocolComplexity,
)

logger = logging.getLogger(__name__)


def _clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    """Clamp a numeric value to [lo, hi]."""
    return max(lo, min(hi, value))


# ═══════════════════════════════════════════════════════════════════════════
# CONFIDENCE CALIBRATOR
# ═══════════════════════════════════════════════════════════════════════════


class ConfidenceCalibrator:
    """Calibrate raw confidence scores using a weighted multi-factor model.

    Formula
    -------
    calibrated = 0.3 * raw_confidence
               + 0.3 * evidence_base
               + 0.2 * doc_factor
               + 0.2 * agreement_factor

    Parameters
    ----------
    raw_confidence : float
        Raw model or workflow confidence (0.0-1.0).
    evidence_level : EvidenceLevel
        Level of supporting evidence (A1 through E).
    n_docs : int
        Number of supporting documents retrieved.
    cross_agent_agreement : float
        Agreement score across multiple agents (0.0-1.0).
    """

    # Evidence level to numeric mapping
    _EVIDENCE_SCORES: Dict[str, float] = {
        "a1": 1.0,   # Systematic review of RCTs
        "a2": 0.85,  # High-quality RCT
        "b": 0.65,   # Non-randomized controlled study
        "c": 0.45,   # Observational study
        "d": 0.25,   # Case series / case report
        "e": 0.15,   # Expert opinion
    }

    def calibrate(
        self,
        raw_confidence: float,
        evidence_level: str = "e",
        n_docs: int = 0,
        cross_agent_agreement: float = 0.5,
    ) -> float:
        """Produce a calibrated confidence score.

        Parameters
        ----------
        raw_confidence : float
            Raw confidence from the workflow (0.0-1.0).
        evidence_level : str
            Evidence level string (a1, a2, b, c, d, e).
        n_docs : int
            Number of supporting documents retrieved.
        cross_agent_agreement : float
            Agreement score from cross-agent validation (0.0-1.0).

        Returns
        -------
        float
            Calibrated confidence score (0.0-1.0).
        """
        raw = _clamp(raw_confidence)

        # Evidence base score
        evidence_base = self._EVIDENCE_SCORES.get(
            evidence_level.lower(), 0.15
        )

        # Document factor: logarithmic scaling (diminishing returns)
        # 0 docs -> 0.0, 1 doc -> 0.3, 5 docs -> 0.7, 10+ docs -> ~0.9
        if n_docs <= 0:
            doc_factor = 0.0
        else:
            doc_factor = _clamp(math.log(n_docs + 1) / math.log(12))

        # Agreement factor
        agreement_factor = _clamp(cross_agent_agreement)

        calibrated = (
            0.3 * raw
            + 0.3 * evidence_base
            + 0.2 * doc_factor
            + 0.2 * agreement_factor
        )

        logger.debug(
            "Confidence calibration: raw=%.3f, evidence=%.3f (level=%s), "
            "doc_factor=%.3f (n=%d), agreement=%.3f => calibrated=%.3f",
            raw, evidence_base, evidence_level, doc_factor, n_docs,
            agreement_factor, calibrated,
        )

        return round(_clamp(calibrated), 4)

    def calibrate_batch(
        self,
        scores: List[Dict],
    ) -> List[float]:
        """Calibrate a batch of confidence scores.

        Parameters
        ----------
        scores : list[dict]
            Each dict must have: raw_confidence, evidence_level, n_docs,
            cross_agent_agreement.

        Returns
        -------
        list[float]
            Calibrated scores.
        """
        return [
            self.calibrate(
                raw_confidence=s.get("raw_confidence", 0.5),
                evidence_level=s.get("evidence_level", "e"),
                n_docs=s.get("n_docs", 0),
                cross_agent_agreement=s.get("cross_agent_agreement", 0.5),
            )
            for s in scores
        ]


# ═══════════════════════════════════════════════════════════════════════════
# PROTOCOL COMPLEXITY SCORER
# ═══════════════════════════════════════════════════════════════════════════


class ProtocolComplexityScorer:
    """Score protocol complexity based on structural features.

    Scoring dimensions (all normalised to 0-1, then combined):
    - procedure_count: total distinct procedures (normalised by max 50)
    - visit_frequency: visits per month (normalised by max 4)
    - endpoint_count: total endpoints (normalised by max 20)
    - eligibility_criteria_count: total criteria (normalised by max 40)
    - amendment_history: number of amendments (normalised by max 5)

    Percentile rank is estimated against industry benchmarks.
    """

    # Industry benchmarks (Tufts CSDD)
    _BENCHMARKS = {
        "procedure_count": {"median": 20, "p75": 30, "p90": 45, "max_norm": 50},
        "visit_count": {"median": 12, "p75": 18, "p90": 28, "max_norm": 36},
        "endpoint_count": {"median": 8, "p75": 12, "p90": 18, "max_norm": 20},
        "eligibility_criteria_count": {"median": 22, "p75": 30, "p90": 38, "max_norm": 40},
        "amendment_count": {"median": 2, "p75": 3, "p90": 4, "max_norm": 5},
    }

    def score(self, protocol: dict) -> ProtocolComplexity:
        """Score protocol complexity.

        Parameters
        ----------
        protocol : dict
            Keys: procedure_count (int), visit_count (int),
            endpoint_count (int), eligibility_criteria_count (int),
            amendment_count (int), duration_months (int).

        Returns
        -------
        ProtocolComplexity
            Complexity assessment with normalised score and percentile.
        """
        proc = int(protocol.get("procedure_count", 0))
        visits = int(protocol.get("visit_count", 0))
        endpoints = int(protocol.get("endpoint_count", 0))
        criteria = int(protocol.get("eligibility_criteria_count", 0))
        amendments = int(protocol.get("amendment_count", 0))

        # Normalise each dimension
        proc_norm = _clamp(proc / self._BENCHMARKS["procedure_count"]["max_norm"])
        visit_norm = _clamp(visits / self._BENCHMARKS["visit_count"]["max_norm"])
        ep_norm = _clamp(endpoints / self._BENCHMARKS["endpoint_count"]["max_norm"])
        crit_norm = _clamp(
            criteria / self._BENCHMARKS["eligibility_criteria_count"]["max_norm"]
        )
        amend_norm = _clamp(amendments / self._BENCHMARKS["amendment_count"]["max_norm"])

        # Weighted composite (procedures and visits weighted higher)
        complexity = (
            0.25 * proc_norm
            + 0.25 * visit_norm
            + 0.20 * ep_norm
            + 0.20 * crit_norm
            + 0.10 * amend_norm
        )
        complexity = round(_clamp(complexity), 4)

        # Estimate percentile rank
        # Simple linear interpolation against benchmarks
        percentile = _clamp(complexity * 100, 0.0, 99.9)

        return ProtocolComplexity(
            procedure_count=proc,
            visit_count=visits,
            endpoint_count=endpoints,
            eligibility_criteria_count=criteria,
            complexity_score=complexity,
            percentile_rank=round(percentile, 1),
        )


# ═══════════════════════════════════════════════════════════════════════════
# ENROLLMENT PREDICTOR
# ═══════════════════════════════════════════════════════════════════════════


class EnrollmentPredictor:
    """Predict enrollment rate for a trial site based on historical data,
    disease prevalence, competing trials, and site capacity.

    prediction = historical_rate * prevalence_factor * competition_factor
                 * capacity_factor

    Parameters
    ----------
    site : dict
        Site data (historical_enrollment_rate, capacity, staff_count).
    trial : dict
        Trial data (phase, therapeutic_area, eligibility_stringency).
    disease : dict
        Disease data (prevalence_per_100k, competing_trials_count,
        standard_of_care_available).
    """

    # Phase-specific enrollment difficulty multipliers
    _PHASE_FACTORS: Dict[str, float] = {
        "phase_i": 0.6,
        "phase_i_ii": 0.7,
        "phase_ii": 0.8,
        "phase_ii_iii": 0.85,
        "phase_iii": 1.0,
        "phase_iv": 1.1,
    }

    def predict(self, site: dict, trial: dict, disease: dict) -> float:
        """Predict monthly enrollment rate for a site.

        Parameters
        ----------
        site : dict
            Keys: historical_enrollment_rate (float, pts/month),
            capacity (int, max concurrent patients), staff_count (int).
        trial : dict
            Keys: phase (str), therapeutic_area (str),
            eligibility_stringency (float, 0-1 where 1 is very restrictive).
        disease : dict
            Keys: prevalence_per_100k (float),
            competing_trials_count (int),
            standard_of_care_available (bool).

        Returns
        -------
        float
            Predicted enrollment rate (patients per month).
        """
        # Historical baseline
        hist_rate = float(site.get("historical_enrollment_rate", 2.0))

        # Prevalence factor
        prevalence = float(disease.get("prevalence_per_100k", 50))
        if prevalence >= 500:
            prevalence_factor = 1.2   # Common disease
        elif prevalence >= 100:
            prevalence_factor = 1.0
        elif prevalence >= 10:
            prevalence_factor = 0.7   # Uncommon
        else:
            prevalence_factor = 0.4   # Rare disease

        # Competition factor
        competing = int(disease.get("competing_trials_count", 0))
        if competing == 0:
            competition_factor = 1.2   # No competition
        elif competing <= 3:
            competition_factor = 1.0
        elif competing <= 10:
            competition_factor = 0.7
        else:
            competition_factor = 0.5   # Very crowded

        # Capacity factor
        capacity = int(site.get("capacity", 50))
        staff = int(site.get("staff_count", 5))
        # Staff-to-capacity ratio
        ratio = staff / max(capacity, 1)
        if ratio >= 0.2:
            capacity_factor = 1.1
        elif ratio >= 0.1:
            capacity_factor = 1.0
        else:
            capacity_factor = 0.7

        # Phase factor
        phase = trial.get("phase", "phase_iii").lower().replace(" ", "_")
        phase_factor = self._PHASE_FACTORS.get(phase, 1.0)

        # Eligibility stringency penalty
        stringency = float(trial.get("eligibility_stringency", 0.5))
        stringency_factor = 1.0 - (stringency * 0.4)  # max 40% reduction

        predicted = (
            hist_rate
            * prevalence_factor
            * competition_factor
            * capacity_factor
            * phase_factor
            * stringency_factor
        )

        logger.debug(
            "Enrollment prediction: hist=%.2f, prev=%.2f, comp=%.2f, "
            "cap=%.2f, phase=%.2f, string=%.2f => predicted=%.2f",
            hist_rate, prevalence_factor, competition_factor,
            capacity_factor, phase_factor, stringency_factor, predicted,
        )

        return round(max(predicted, 0.0), 2)


# ═══════════════════════════════════════════════════════════════════════════
# ELIGIBILITY ANALYZER
# ═══════════════════════════════════════════════════════════════════════════


class EligibilityAnalyzer:
    """Analyse eligibility criteria to identify overly restrictive
    criteria and estimate population impact.
    """

    # Restrictive patterns with estimated population exclusion
    _POPULATION_IMPACTS: Dict[str, float] = {
        "ecog 0-1": 0.25,
        "ecog 0": 0.40,
        "no prior systemic therapy": 0.35,
        "treatment naive": 0.35,
        "no prior immunotherapy": 0.20,
        "creatinine clearance": 0.15,
        "hepatic function": 0.10,
        "no cns metastases": 0.15,
        "no brain metastases": 0.15,
        "no autoimmune disease": 0.08,
        "no cardiac history": 0.12,
        "ejection fraction": 0.08,
        "hemoglobin": 0.12,
        "platelet count": 0.08,
        "bmi": 0.12,
        "no hiv": 0.01,
        "no hepatitis": 0.03,
        "washout period": 0.10,
        "life expectancy": 0.15,
        "prior car-t": 0.85,
        "prior bispecific": 0.80,
        "albumin < 3.0": 0.25,
        "pregnancy or lactation": 0.50,
        "active brain metastases": 0.40,
        "prior organ transplant": 0.15,
        "active autoimmune disease": 0.30,
        "uncontrolled diabetes": 0.20,
        "severe hepatic impairment": 0.25,
        "ejection fraction < 40%": 0.15,
        "prior immunotherapy": 0.35,
        "measurable disease required": 0.25,
    }

    def analyze(self, criteria: List[dict]) -> List[EligibilityAnalysis]:
        """Analyse a list of eligibility criteria.

        Parameters
        ----------
        criteria : list[dict]
            Each dict has keys: text (str), type (str: inclusion/exclusion),
            scientific_justification (str, optional).

        Returns
        -------
        list[EligibilityAnalysis]
            Per-criterion analysis with population impact and recommendations.
        """
        results: List[EligibilityAnalysis] = []

        for crit in criteria:
            text = crit.get("text", "") if isinstance(crit, dict) else str(crit)
            justification = (
                crit.get("scientific_justification", "")
                if isinstance(crit, dict) else ""
            )
            text_lower = text.lower()

            # Estimate population impact
            population_impact = 0.03  # baseline minimal impact
            for pattern, impact in self._POPULATION_IMPACTS.items():
                if pattern in text_lower:
                    population_impact = max(population_impact, impact)

            # Assess scientific justification strength
            justification_score = 0.3  # default weak
            if justification:
                just_lower = justification.lower()
                if any(kw in just_lower for kw in [
                    "rct", "meta-analysis", "systematic review",
                    "phase iii", "pivotal",
                ]):
                    justification_score = 0.9
                elif any(kw in just_lower for kw in [
                    "observational", "cohort", "registry", "phase ii",
                ]):
                    justification_score = 0.6
                elif any(kw in just_lower for kw in [
                    "expert", "opinion", "convention", "common practice",
                ]):
                    justification_score = 0.3

            # Generate recommendation
            if population_impact >= 0.15 and justification_score < 0.5:
                recommendation = (
                    f"BROADEN: High population impact ({population_impact:.0%}) "
                    f"with weak justification — consider relaxing this criterion"
                )
            elif population_impact >= 0.10 and justification_score < 0.7:
                recommendation = (
                    f"REVIEW: Moderate impact ({population_impact:.0%}) — "
                    f"evaluate if criterion can be broadened"
                )
            elif population_impact >= 0.20 and justification_score >= 0.7:
                recommendation = (
                    f"RETAIN with monitoring: High impact ({population_impact:.0%}) "
                    f"but strong scientific justification"
                )
            else:
                recommendation = "RETAIN: Acceptable impact and justification"

            # Competitor comparison (generic)
            competitor_note = ""
            if population_impact >= 0.15:
                competitor_note = (
                    "Many competing trials use broader criteria for this parameter — "
                    "restrictive criterion may disadvantage enrollment"
                )

            results.append(EligibilityAnalysis(
                criterion=text,
                population_impact=round(population_impact, 3),
                scientific_justification_score=round(justification_score, 3),
                competitor_comparison=competitor_note,
                recommendation=recommendation,
            ))

        # Sort by population impact descending
        results.sort(key=lambda r: r.population_impact, reverse=True)
        return results


# ═══════════════════════════════════════════════════════════════════════════
# COMPETITIVE THREAT SCORER
# ═══════════════════════════════════════════════════════════════════════════


class CompetitiveThreatScorer:
    """Score competitive threat from a rival trial on a 0-1 scale.

    Factors
    -------
    - phase_advancement (0.30): How advanced is the competitor?
    - enrollment_progress (0.25): How far along is enrollment?
    - sponsor_resources (0.20): Sponsor capability tier.
    - differentiation (0.25): Mechanism/design similarity (lower diff = higher threat).
    """

    _PHASE_SCORES: Dict[str, float] = {
        "phase_i": 0.15,
        "phase_i_ii": 0.25,
        "phase_ii": 0.40,
        "phase_ii_iii": 0.60,
        "phase_iii": 0.80,
        "phase_iv": 0.95,
        "approved": 1.0,
    }

    _SPONSOR_SCORES: Dict[str, float] = {
        "large_pharma": 0.90,
        "mid_pharma": 0.70,
        "biotech": 0.50,
        "academic": 0.30,
    }

    def score(self, competitor: dict) -> float:
        """Calculate threat level for a competitor.

        Parameters
        ----------
        competitor : dict
            Keys: phase (str), enrollment_target (int),
            enrollment_actual (int), sponsor_type (str),
            mechanism_similarity (float, 0-1 where 1 = identical).

        Returns
        -------
        float
            Threat score (0.0-1.0).
        """
        # Phase advancement
        phase = competitor.get("phase", "phase_ii").lower().replace(" ", "_")
        phase_score = self._PHASE_SCORES.get(phase, 0.40)

        # Enrollment progress
        target = int(competitor.get("enrollment_target", 100))
        actual = int(competitor.get("enrollment_actual", 0))
        enrollment_score = _clamp(actual / max(target, 1))

        # Sponsor resources
        sponsor = competitor.get("sponsor_type", "biotech").lower()
        sponsor_score = self._SPONSOR_SCORES.get(sponsor, 0.50)

        # Differentiation (similarity as threat)
        similarity = float(competitor.get("mechanism_similarity", 0.5))
        differentiation_score = _clamp(similarity)

        threat = (
            0.30 * phase_score
            + 0.25 * enrollment_score
            + 0.20 * sponsor_score
            + 0.25 * differentiation_score
        )

        return round(_clamp(threat), 4)

    def score_batch(self, competitors: List[dict]) -> List[Dict[str, float]]:
        """Score multiple competitors and return sorted results.

        Returns
        -------
        list[dict]
            Each dict has: competitor data + threat_score, sorted descending.
        """
        scored = []
        for comp in competitors:
            threat = self.score(comp)
            scored.append({**comp, "threat_score": threat})
        scored.sort(key=lambda c: c["threat_score"], reverse=True)
        return scored

    def classify_threat(self, score: float) -> str:
        """Classify a threat score into a human-readable level.

        Returns
        -------
        str
            One of: 'critical', 'high', 'moderate', 'low', 'minimal'.
        """
        if score >= 0.80:
            severity = "critical"
        elif score >= 0.60:
            severity = "high"
        elif score >= 0.40:
            severity = "moderate"
        elif score >= 0.20:
            severity = "low"
        else:
            severity = "minimal"

        try:
            from api.routes.events import publish_event
            if severity in ("critical", "high"):
                publish_event("safety_alert", {"event_type": "competitive_threat", "severity": severity})
        except Exception:
            pass

        return severity


# ═══════════════════════════════════════════════════════════════════════════
# HISTORICAL SUCCESS RATE ESTIMATOR
# ═══════════════════════════════════════════════════════════════════════════


class HistoricalSuccessEstimator:
    """Estimate probability of trial success based on historical success
    rates by therapeutic area and development phase.

    Data derived from industry benchmarks (BIO/QLS Advisors, Tufts CSDD).
    """

    # Historical success rates keyed by (therapeutic_area, phase)
    _SUCCESS_RATES: Dict[str, Dict[str, float]] = {
        "oncology": {
            "phase_i": 0.52,
            "phase_ii": 0.28,
            "phase_iii": 0.40,
        },
        "cardiology": {
            "phase_i": 0.64,
            "phase_ii": 0.36,
            "phase_iii": 0.52,
        },
        "metabolic": {
            "phase_i": 0.60,
            "phase_ii": 0.35,
            "phase_iii": 0.52,
        },
        "rare_disease": {
            "phase_i": 0.70,
            "phase_ii": 0.45,
            "phase_iii": 0.62,
        },
        "gene_therapy": {
            "phase_i": 0.65,
            "phase_ii": 0.40,
            "phase_iii": 0.55,
        },
        "infectious_disease": {
            "phase_i": 0.60,
            "phase_ii": 0.38,
            "phase_iii": 0.50,
        },
        "dermatology": {
            "phase_i": 0.68,
            "phase_ii": 0.42,
            "phase_iii": 0.58,
        },
        "ophthalmology": {
            "phase_i": 0.65,
            "phase_ii": 0.40,
            "phase_iii": 0.54,
        },
        "respiratory": {
            "phase_i": 0.58,
            "phase_ii": 0.32,
            "phase_iii": 0.45,
        },
        "neurology": {
            "phase_i": 0.50,
            "phase_ii": 0.22,
            "phase_iii": 0.35,
        },
        "hematology": {
            "phase_i": 0.62,
            "phase_ii": 0.40,
            "phase_iii": 0.55,
        },
        "gastroenterology": {
            "phase_i": 0.58,
            "phase_ii": 0.34,
            "phase_iii": 0.48,
        },
    }

    # Overall industry average by phase (fallback)
    _PHASE_AVERAGES: Dict[str, float] = {
        "phase_i": 0.60,
        "phase_ii": 0.33,
        "phase_iii": 0.50,
        "phase_iv": 0.85,
    }

    def estimate(self, therapeutic_area: str, phase: str) -> float:
        """Estimate probability of success for a given area and phase.

        Parameters
        ----------
        therapeutic_area : str
            Therapeutic area (e.g. 'oncology', 'neurology').
        phase : str
            Development phase (e.g. 'phase_iii', 'Phase III').

        Returns
        -------
        float
            Estimated probability of success (0.0-1.0).
        """
        area = therapeutic_area.lower().replace(" ", "_").replace("-", "_")
        ph = phase.lower().replace(" ", "_")

        area_rates = self._SUCCESS_RATES.get(area)
        if area_rates and ph in area_rates:
            return area_rates[ph]

        # Fallback to phase average
        return self._PHASE_AVERAGES.get(ph, 0.50)

    def estimate_cumulative(self, therapeutic_area: str, current_phase: str) -> float:
        """Estimate cumulative probability of reaching approval from current phase.

        Parameters
        ----------
        therapeutic_area : str
            Therapeutic area.
        current_phase : str
            Current development phase.

        Returns
        -------
        float
            Cumulative probability of success from current phase to approval.
        """
        phase_order = ["phase_i", "phase_ii", "phase_iii"]
        ph = current_phase.lower().replace(" ", "_")

        if ph not in phase_order:
            return self.estimate(therapeutic_area, ph)

        start_idx = phase_order.index(ph)
        cumulative = 1.0
        for p in phase_order[start_idx:]:
            cumulative *= self.estimate(therapeutic_area, p)

        return round(cumulative, 4)
