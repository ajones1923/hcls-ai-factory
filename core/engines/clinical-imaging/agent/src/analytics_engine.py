"""NVIDIA RAPIDS GPU-accelerated analytics for medical imaging.

Provides population-level imaging analytics using cuDF (GPU DataFrame)
with transparent fallback to pandas when RAPIDS is unavailable.

Apache 2.0 Licensed.

Author: Adam Jones
Date: April 2026
"""

import random
import time
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# ANALYTICS DATA MODELS
# ═══════════════════════════════════════════════════════════════════════


class PopulationSummary(BaseModel):
    """Aggregate statistics across all registered imaging studies."""

    total_studies: int
    modality_distribution: Dict[str, int]
    body_region_distribution: Dict[str, int]
    severity_distribution: Dict[str, int]
    finding_prevalence: Dict[str, float]  # finding -> percentage
    mean_processing_time_ms: float
    date_range: Optional[Dict[str, str]] = None  # {"start": ..., "end": ...}
    studies_with_critical_findings: int
    critical_finding_rate: float


class CohortCriteria(BaseModel):
    """Criteria for filtering imaging study cohorts."""

    modality: Optional[str] = None
    body_region: Optional[str] = None
    min_finding_size_mm: Optional[float] = None
    severity: Optional[str] = None  # "critical", "urgent", etc.
    start_date: Optional[str] = None
    end_date: Optional[str] = None
    finding_type: Optional[str] = None  # "nodule", "hemorrhage", etc.
    age_min: Optional[int] = None
    age_max: Optional[int] = None


class CohortResult(BaseModel):
    """Result of a cohort query against the study registry."""

    criteria: CohortCriteria
    matching_studies: int
    total_studies: int
    match_rate: float
    studies: List[Dict[str, Any]]  # summary of each matching study


class TrendPoint(BaseModel):
    """A single data point in a temporal trend series."""

    period: str  # "2025-Q1", "2025-W12", etc.
    count: int
    rate: float


class TrendResult(BaseModel):
    """Temporal trend analysis result."""

    metric: str  # "critical_findings", "ct_volume", etc.
    granularity: str  # "weekly", "monthly", "quarterly"
    data_points: List[TrendPoint]
    trend_direction: str  # "increasing", "decreasing", "stable"


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

_MODALITIES = ["ct", "mri", "xray", "cxr", "ultrasound", "pet", "pet_ct", "mammography"]
_BODY_REGIONS = [
    "head", "neck", "chest", "abdomen", "pelvis", "spine",
    "extremity", "brain", "cardiac", "breast",
]
_SEVERITIES = ["critical", "urgent", "significant", "routine", "normal"]
_FINDING_TYPES = [
    "nodule", "hemorrhage", "mass", "fracture", "consolidation",
    "effusion", "pneumothorax", "edema", "atelectasis", "lesion",
    "infarct", "stenosis", "calcification", "normal",
]


# ═══════════════════════════════════════════════════════════════════════
# ANALYTICS ENGINE
# ═══════════════════════════════════════════════════════════════════════


class ImagingAnalyticsEngine:
    """GPU-accelerated imaging analytics via RAPIDS cuDF.

    Falls back to pandas when RAPIDS is not available.
    """

    def __init__(self) -> None:
        self._rapids_available: bool = self._check_rapids()
        self._df_lib = self._get_df_library()
        self._study_registry: List[Dict[str, Any]] = []  # In-memory store for demo
        backend = "cuDF (GPU)" if self._rapids_available else "pandas (CPU)"
        logger.info(f"ImagingAnalyticsEngine initialized — backend: {backend}")

    # ── Library detection ──────────────────────────────────────────

    def _check_rapids(self) -> bool:
        """Check whether NVIDIA RAPIDS cuDF is importable."""
        try:
            import cudf  # noqa: F401

            return True
        except ImportError:
            logger.debug("RAPIDS cuDF not available, falling back to pandas")
            return False

    def _get_df_library(self):
        """Return cudf if available, otherwise pandas."""
        if self._rapids_available:
            import cudf

            return cudf
        import pandas

        return pandas

    @property
    def backend_name(self) -> str:
        """Return human-readable backend name."""
        return "cudf" if self._rapids_available else "pandas"

    # ── Study registration ─────────────────────────────────────────

    def register_study(self, study: Dict[str, Any]) -> None:
        """Register a completed imaging study for analytics."""
        self._study_registry.append(study)

    def register_studies(self, studies: List[Dict[str, Any]]) -> None:
        """Batch register studies."""
        self._study_registry.extend(studies)

    @property
    def study_count(self) -> int:
        """Number of registered studies."""
        return len(self._study_registry)

    # ── Population summary ─────────────────────────────────────────

    def population_summary(self) -> PopulationSummary:
        """Aggregate statistics across all registered studies."""
        total = len(self._study_registry)
        if total == 0:
            return PopulationSummary(
                total_studies=0,
                modality_distribution={},
                body_region_distribution={},
                severity_distribution={},
                finding_prevalence={},
                mean_processing_time_ms=0.0,
                date_range=None,
                studies_with_critical_findings=0,
                critical_finding_rate=0.0,
            )

        df = self._df_lib.DataFrame(self._study_registry)

        # Modality distribution
        modality_dist: Dict[str, int] = {}
        if "modality" in df.columns:
            vc = df["modality"].value_counts()
            modality_dist = {str(k): int(v) for k, v in zip(vc.index, vc.values)}

        # Body region distribution
        region_dist: Dict[str, int] = {}
        if "body_region" in df.columns:
            vc = df["body_region"].value_counts()
            region_dist = {str(k): int(v) for k, v in zip(vc.index, vc.values)}

        # Severity distribution
        severity_dist: Dict[str, int] = {}
        if "severity" in df.columns:
            vc = df["severity"].value_counts()
            severity_dist = {str(k): int(v) for k, v in zip(vc.index, vc.values)}

        # Finding prevalence (percentage)
        finding_prev: Dict[str, float] = {}
        if "finding_type" in df.columns:
            vc = df["finding_type"].value_counts()
            for k, v in zip(vc.index, vc.values):
                finding_prev[str(k)] = round(float(v) / total * 100, 2)

        # Mean processing time
        mean_time = 0.0
        if "processing_time_ms" in df.columns:
            mean_time = float(df["processing_time_ms"].mean())

        # Date range
        date_range = None
        if "study_date" in df.columns:
            dates = sorted([s["study_date"] for s in self._study_registry if s.get("study_date")])
            if dates:
                date_range = {"start": dates[0], "end": dates[-1]}

        # Critical findings
        critical_count = severity_dist.get("critical", 0)
        critical_rate = round(critical_count / total, 4) if total > 0 else 0.0

        return PopulationSummary(
            total_studies=total,
            modality_distribution=modality_dist,
            body_region_distribution=region_dist,
            severity_distribution=severity_dist,
            finding_prevalence=finding_prev,
            mean_processing_time_ms=round(mean_time, 2),
            date_range=date_range,
            studies_with_critical_findings=critical_count,
            critical_finding_rate=critical_rate,
        )

    # ── Cohort query ───────────────────────────────────────────────

    def cohort_query(self, criteria: CohortCriteria) -> CohortResult:
        """Find studies matching specific criteria."""
        total = len(self._study_registry)
        if total == 0:
            return CohortResult(
                criteria=criteria,
                matching_studies=0,
                total_studies=0,
                match_rate=0.0,
                studies=[],
            )

        df = self._df_lib.DataFrame(self._study_registry)
        mask = self._df_lib.Series([True] * len(df))

        if criteria.modality and "modality" in df.columns:
            mask = mask & (df["modality"] == criteria.modality)
        if criteria.body_region and "body_region" in df.columns:
            mask = mask & (df["body_region"] == criteria.body_region)
        if criteria.severity and "severity" in df.columns:
            mask = mask & (df["severity"] == criteria.severity)
        if criteria.finding_type and "finding_type" in df.columns:
            mask = mask & (df["finding_type"] == criteria.finding_type)
        if criteria.min_finding_size_mm is not None and "finding_size_mm" in df.columns:
            mask = mask & (df["finding_size_mm"] >= criteria.min_finding_size_mm)
        if criteria.age_min is not None and "patient_age" in df.columns:
            mask = mask & (df["patient_age"] >= criteria.age_min)
        if criteria.age_max is not None and "patient_age" in df.columns:
            mask = mask & (df["patient_age"] <= criteria.age_max)
        if criteria.start_date and "study_date" in df.columns:
            mask = mask & (df["study_date"] >= criteria.start_date)
        if criteria.end_date and "study_date" in df.columns:
            mask = mask & (df["study_date"] <= criteria.end_date)

        filtered = df[mask]
        n_match = int(len(filtered))

        # Convert matched rows to list of dicts
        if self._rapids_available:
            matched_studies = filtered.to_pandas().to_dict(orient="records")
        else:
            matched_studies = filtered.to_dict(orient="records")

        match_rate = round(n_match / total, 4) if total > 0 else 0.0

        return CohortResult(
            criteria=criteria,
            matching_studies=n_match,
            total_studies=total,
            match_rate=match_rate,
            studies=matched_studies,
        )

    # ── Temporal trends ────────────────────────────────────────────

    def temporal_trends(
        self, metric: str, granularity: str = "monthly"
    ) -> TrendResult:
        """Compute temporal trends for a given metric.

        Supported metrics:
          - "critical_findings": count of critical-severity studies
          - "ct_volume": count of CT studies
          - "total_volume": count of all studies

        Supported granularity: "weekly", "monthly", "quarterly"
        """
        if not self._study_registry:
            return TrendResult(
                metric=metric,
                granularity=granularity,
                data_points=[],
                trend_direction="stable",
            )

        df = self._df_lib.DataFrame(self._study_registry)

        # Filter by metric
        if metric == "critical_findings" and "severity" in df.columns:
            if self._rapids_available:
                subset = df[df["severity"] == "critical"].to_pandas()
            else:
                subset = df[df["severity"] == "critical"]
        elif metric == "ct_volume" and "modality" in df.columns:
            if self._rapids_available:
                subset = df[df["modality"] == "ct"].to_pandas()
            else:
                subset = df[df["modality"] == "ct"]
        else:
            # total_volume or unknown metric — use all studies
            subset = df.to_pandas() if self._rapids_available else df

        all_df = df.to_pandas() if self._rapids_available else df

        # Group by period
        data_points: List[TrendPoint] = []
        period_counts: Dict[str, int] = {}
        period_totals: Dict[str, int] = {}

        for _, row in subset.iterrows():
            date_str = row.get("study_date", "")
            if not date_str:
                continue
            period_key = self._date_to_period(date_str, granularity)
            period_counts[period_key] = period_counts.get(period_key, 0) + 1

        for _, row in all_df.iterrows():
            date_str = row.get("study_date", "")
            if not date_str:
                continue
            period_key = self._date_to_period(date_str, granularity)
            period_totals[period_key] = period_totals.get(period_key, 0) + 1

        for period in sorted(period_counts.keys()):
            count = period_counts[period]
            total_in_period = period_totals.get(period, 1)
            rate = round(count / total_in_period, 4) if total_in_period > 0 else 0.0
            data_points.append(TrendPoint(period=period, count=count, rate=rate))

        # Determine trend direction
        trend_direction = self._compute_trend_direction(data_points)

        return TrendResult(
            metric=metric,
            granularity=granularity,
            data_points=data_points,
            trend_direction=trend_direction,
        )

    @staticmethod
    def _date_to_period(date_str: str, granularity: str) -> str:
        """Convert a date string (YYYY-MM-DD) to a period key."""
        try:
            dt = datetime.strptime(date_str[:10], "%Y-%m-%d")
        except (ValueError, TypeError):
            return "unknown"

        if granularity == "weekly":
            iso = dt.isocalendar()
            return f"{iso[0]}-W{iso[1]:02d}"
        elif granularity == "quarterly":
            quarter = (dt.month - 1) // 3 + 1
            return f"{dt.year}-Q{quarter}"
        else:  # monthly
            return f"{dt.year}-{dt.month:02d}"

    @staticmethod
    def _compute_trend_direction(data_points: List[TrendPoint]) -> str:
        """Determine if trend is increasing, decreasing, or stable."""
        if len(data_points) < 2:
            return "stable"

        counts = [dp.count for dp in data_points]
        # Simple linear comparison: first half vs second half
        mid = len(counts) // 2
        first_half_avg = sum(counts[:mid]) / max(mid, 1)
        second_half_avg = sum(counts[mid:]) / max(len(counts) - mid, 1)

        diff_pct = (second_half_avg - first_half_avg) / max(first_half_avg, 1)
        if diff_pct > 0.10:
            return "increasing"
        elif diff_pct < -0.10:
            return "decreasing"
        return "stable"

    # ── Finding correlation ────────────────────────────────────────

    def finding_correlation(self, finding_a: str, finding_b: str) -> Dict[str, Any]:
        """Compute co-occurrence rate between two finding types.

        Uses patient-level co-occurrence: how often a patient with
        finding_a also has finding_b.
        """
        total = len(self._study_registry)
        if total == 0:
            return {
                "finding_a": finding_a,
                "finding_b": finding_b,
                "co_occurrence_count": 0,
                "finding_a_count": 0,
                "finding_b_count": 0,
                "co_occurrence_rate": 0.0,
                "total_studies": 0,
            }

        # Group findings by patient_id
        patient_findings: Dict[str, set] = {}
        for study in self._study_registry:
            pid = study.get("patient_id", study.get("study_id", "unknown"))
            ft = study.get("finding_type", "")
            patient_findings.setdefault(pid, set()).add(ft)

        a_count = sum(1 for fts in patient_findings.values() if finding_a in fts)
        b_count = sum(1 for fts in patient_findings.values() if finding_b in fts)
        co_count = sum(
            1 for fts in patient_findings.values()
            if finding_a in fts and finding_b in fts
        )
        co_rate = round(co_count / max(a_count, 1), 4)

        return {
            "finding_a": finding_a,
            "finding_b": finding_b,
            "co_occurrence_count": co_count,
            "finding_a_count": a_count,
            "finding_b_count": b_count,
            "co_occurrence_rate": co_rate,
            "total_studies": total,
        }

    # ── Severity by modality ───────────────────────────────────────

    def severity_by_modality(self) -> Dict[str, Dict[str, int]]:
        """Cross-tabulation of severity by modality.

        Returns: {"ct": {"critical": 5, "urgent": 10, ...}, ...}
        """
        if not self._study_registry:
            return {}

        result: Dict[str, Dict[str, int]] = {}
        for study in self._study_registry:
            mod = study.get("modality", "unknown")
            sev = study.get("severity", "unknown")
            result.setdefault(mod, {})
            result[mod][sev] = result[mod].get(sev, 0) + 1

        return result

    # ── AI concordance ─────────────────────────────────────────────

    def ai_concordance(self) -> Dict[str, Any]:
        """Compare AI workflow results vs radiologist reports.

        Looks for studies that have both 'ai_severity' and
        'radiologist_severity' fields and computes agreement rate.
        """
        paired = [
            s for s in self._study_registry
            if s.get("ai_severity") and s.get("radiologist_severity")
        ]
        total_paired = len(paired)
        if total_paired == 0:
            return {
                "total_paired": 0,
                "agreement_count": 0,
                "agreement_rate": 0.0,
                "disagreements": [],
            }

        agreements = 0
        disagreements: List[Dict[str, Any]] = []
        for s in paired:
            if s["ai_severity"] == s["radiologist_severity"]:
                agreements += 1
            else:
                disagreements.append({
                    "study_id": s.get("study_id", "unknown"),
                    "ai_severity": s["ai_severity"],
                    "radiologist_severity": s["radiologist_severity"],
                })

        return {
            "total_paired": total_paired,
            "agreement_count": agreements,
            "agreement_rate": round(agreements / total_paired, 4),
            "disagreements": disagreements[:50],  # cap for response size
        }

    # ── Demo data generation ───────────────────────────────────────

    def generate_demo_data(self, n_studies: int = 500, seed: int = 42) -> int:
        """Generate synthetic imaging study data for demonstration.

        Returns the number of studies generated.
        """
        rng = random.Random(seed)
        base_date = datetime(2025, 1, 1, tzinfo=timezone.utc)
        studies: List[Dict[str, Any]] = []

        for i in range(n_studies):
            modality = rng.choice(_MODALITIES)
            body_region = rng.choice(_BODY_REGIONS)
            severity = rng.choices(
                _SEVERITIES,
                weights=[5, 10, 20, 40, 25],
                k=1,
            )[0]
            finding_type = rng.choice(_FINDING_TYPES)
            patient_age = rng.randint(18, 95)
            study_date = (base_date + timedelta(days=rng.randint(0, 365))).strftime(
                "%Y-%m-%d"
            )
            processing_time = round(rng.uniform(50.0, 5000.0), 1)
            finding_size = round(rng.uniform(1.0, 80.0), 1) if finding_type != "normal" else 0.0

            study = {
                "study_id": f"STUDY-{i:06d}",
                "patient_id": f"PAT-{rng.randint(1, n_studies // 3):05d}",
                "modality": modality,
                "body_region": body_region,
                "severity": severity,
                "finding_type": finding_type,
                "finding_size_mm": finding_size,
                "patient_age": patient_age,
                "study_date": study_date,
                "processing_time_ms": processing_time,
            }

            # 20% of studies have AI + radiologist severity for concordance
            if rng.random() < 0.20:
                study["ai_severity"] = severity
                # Radiologist agrees 80% of the time
                if rng.random() < 0.80:
                    study["radiologist_severity"] = severity
                else:
                    study["radiologist_severity"] = rng.choice(
                        [s for s in _SEVERITIES if s != severity]
                    )

            studies.append(study)

        self.register_studies(studies)
        logger.info(
            f"Generated {n_studies} synthetic imaging studies "
            f"(backend={self.backend_name})"
        )
        return n_studies
