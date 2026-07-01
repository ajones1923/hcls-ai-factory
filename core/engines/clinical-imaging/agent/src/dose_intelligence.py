"""Radiation dose tracking and optimization.

Tracks cumulative radiation exposure, compares to Diagnostic Reference
Levels (DRLs), generates alerts, and provides dose optimization
recommendations following ALARA principles.

Pure Python. Apache 2.0 compatible.

Author: Adam Jones
Date: April 2026
"""

import hashlib
import random
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional, Tuple

from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# DATA MODELS
# ═══════════════════════════════════════════════════════════════════════


class DoseRecord(BaseModel):
    """A single radiation dose record for a patient study."""

    patient_id: str
    study_date: str  # ISO format YYYY-MM-DD
    modality: str
    protocol: str
    body_region: str
    effective_dose_msv: float
    dlp_mgy_cm: Optional[float] = None  # Dose-Length Product (CT)
    ctdi_vol_mgy: Optional[float] = None  # CT Dose Index
    dap_mgy_cm2: Optional[float] = None  # Dose-Area Product (fluoroscopy)
    num_exposures: int = 1
    scanner_model: Optional[str] = None
    pediatric: bool = False


class CumulativeDose(BaseModel):
    """Cumulative dose summary for a patient over a time period."""

    patient_id: str
    total_effective_dose_msv: float
    study_count: int
    date_range: Dict[str, str] = Field(default_factory=dict)  # {"first": ..., "last": ...}
    by_modality: Dict[str, float] = Field(default_factory=dict)  # modality -> cumulative mSv
    by_body_region: Dict[str, float] = Field(default_factory=dict)
    alert_level: str = "normal"  # "normal", "elevated", "high", "critical"
    alert_message: Optional[str] = None


class DRLComparison(BaseModel):
    """Comparison of a patient dose to national Diagnostic Reference Levels."""

    protocol: str
    patient_dose_msv: float
    drl_msv: float  # national DRL
    achievable_dose_msv: float  # achievable dose level
    ratio: float  # patient/DRL
    status: str  # "below_achievable", "below_drl", "above_drl", "significantly_above"
    optimization_suggestions: List[str] = Field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════
# DOSE INTELLIGENCE ENGINE
# ═══════════════════════════════════════════════════════════════════════


class DoseIntelligenceEngine:
    """Radiation dose tracking, DRL comparison, and optimization.

    Maintains an in-memory registry of dose records, supports cumulative
    dose calculations, DRL comparisons, and population-level analytics.
    """

    # National Diagnostic Reference Levels (mSv effective dose)
    NATIONAL_DRLS: Dict[str, Dict[str, float]] = {
        "CT Head without contrast": {"drl": 2.5, "achievable": 1.8},
        "CT Head": {"drl": 2.5, "achievable": 1.8},
        "CTA Head/Neck": {"drl": 7.0, "achievable": 4.5},
        "CT Chest": {"drl": 8.0, "achievable": 5.5},
        "CT Chest with contrast": {"drl": 8.0, "achievable": 5.5},
        "Low-dose CT Chest": {"drl": 2.0, "achievable": 1.2},
        "CT Abdomen/Pelvis": {"drl": 15.0, "achievable": 10.0},
        "CT Abdomen/Pelvis with contrast": {"drl": 15.0, "achievable": 10.0},
        "CT Abdomen/Pelvis without contrast": {"drl": 13.0, "achievable": 9.0},
        "CT Abdomen multiphase": {"drl": 25.0, "achievable": 18.0},
        "CTA Coronary": {"drl": 7.0, "achievable": 3.5},
        "CT Chest PE protocol": {"drl": 10.0, "achievable": 6.0},
        "CTA Chest PE protocol": {"drl": 10.0, "achievable": 6.0},
        "CT Spine Lumbar": {"drl": 8.0, "achievable": 5.5},
        "CT Spine Cervical": {"drl": 4.0, "achievable": 2.5},
    }

    # Cumulative dose thresholds (mSv per year)
    ANNUAL_THRESHOLDS: Dict[str, float] = {
        "normal": 20.0,    # < 20 mSv/year
        "elevated": 50.0,  # 20-50 mSv/year
        "high": 100.0,     # 50-100 mSv/year
        # > 100 mSv/year = critical
    }

    PEDIATRIC_MULTIPLIER: float = 0.5  # pediatric thresholds are half of adult

    # Protocols for demo data generation
    _DEMO_PROTOCOLS: List[Dict[str, Any]] = [
        {"protocol": "CT Head without contrast", "modality": "CT", "region": "head", "dose_range": (1.5, 3.0)},
        {"protocol": "CT Chest", "modality": "CT", "region": "chest", "dose_range": (5.0, 10.0)},
        {"protocol": "Low-dose CT Chest", "modality": "CT", "region": "chest", "dose_range": (1.0, 2.0)},
        {"protocol": "CT Abdomen/Pelvis with contrast", "modality": "CT", "region": "abdomen", "dose_range": (8.0, 18.0)},
        {"protocol": "CTA Coronary", "modality": "CT", "region": "cardiac", "dose_range": (2.5, 8.0)},
        {"protocol": "CT Chest PE protocol", "modality": "CT", "region": "chest", "dose_range": (5.0, 12.0)},
        {"protocol": "Chest X-Ray", "modality": "XR", "region": "chest", "dose_range": (0.01, 0.05)},
        {"protocol": "Abdominal X-Ray", "modality": "XR", "region": "abdomen", "dose_range": (0.5, 1.0)},
        {"protocol": "CT Abdomen multiphase", "modality": "CT", "region": "abdomen", "dose_range": (12.0, 30.0)},
        {"protocol": "Mammography", "modality": "MAMMO", "region": "breast", "dose_range": (0.3, 0.7)},
    ]

    def __init__(self) -> None:
        self._dose_registry: List[DoseRecord] = []
        logger.info("DoseIntelligenceEngine initialized")

    # ───────────────────────────────────────────────────────────────────
    # DOSE RECORDING
    # ───────────────────────────────────────────────────────────────────

    def record_dose(self, record: DoseRecord) -> None:
        """Record a dose entry in the registry."""
        self._dose_registry.append(record)
        logger.debug(
            "Recorded dose: patient={}, protocol={}, dose={:.1f} mSv",
            record.patient_id,
            record.protocol,
            record.effective_dose_msv,
        )

    @property
    def registry_size(self) -> int:
        """Return the number of dose records in the registry."""
        return len(self._dose_registry)

    # ───────────────────────────────────────────────────────────────────
    # CUMULATIVE DOSE
    # ───────────────────────────────────────────────────────────────────

    def get_cumulative_dose(
        self,
        patient_id: str,
        period_days: int = 365,
    ) -> CumulativeDose:
        """Calculate cumulative dose for a patient over a given period.

        Parameters
        ----------
        patient_id : str
            Patient identifier.
        period_days : int
            Lookback period in days (default 365 = 1 year).

        Returns
        -------
        CumulativeDose
            Cumulative dose summary with breakdown and alert level.
        """
        cutoff = (datetime.now() - timedelta(days=period_days)).strftime("%Y-%m-%d")
        records = [
            r for r in self._dose_registry
            if r.patient_id == patient_id and r.study_date >= cutoff
        ]

        if not records:
            return CumulativeDose(
                patient_id=patient_id,
                total_effective_dose_msv=0.0,
                study_count=0,
                alert_level="normal",
                alert_message=None,
            )

        total_msv = sum(r.effective_dose_msv for r in records)
        by_modality: Dict[str, float] = {}
        by_region: Dict[str, float] = {}
        dates = []

        for r in records:
            by_modality[r.modality] = by_modality.get(r.modality, 0.0) + r.effective_dose_msv
            by_region[r.body_region] = by_region.get(r.body_region, 0.0) + r.effective_dose_msv
            dates.append(r.study_date)

        dates.sort()
        pediatric = any(r.pediatric for r in records)
        alert_level, alert_message = self._classify_alert_level(total_msv, pediatric=pediatric)

        return CumulativeDose(
            patient_id=patient_id,
            total_effective_dose_msv=round(total_msv, 2),
            study_count=len(records),
            date_range={"first": dates[0], "last": dates[-1]},
            by_modality={k: round(v, 2) for k, v in by_modality.items()},
            by_body_region={k: round(v, 2) for k, v in by_region.items()},
            alert_level=alert_level,
            alert_message=alert_message,
        )

    # ───────────────────────────────────────────────────────────────────
    # DRL COMPARISON
    # ───────────────────────────────────────────────────────────────────

    def compare_to_drl(self, record: DoseRecord) -> DRLComparison:
        """Compare a dose record to national Diagnostic Reference Levels.

        Parameters
        ----------
        record : DoseRecord
            The dose record to compare.

        Returns
        -------
        DRLComparison
            Comparison result with status and optimization suggestions.
        """
        drl_entry = self._find_drl(record.protocol)

        if drl_entry is None:
            logger.debug("No DRL found for protocol '{}'", record.protocol)
            return DRLComparison(
                protocol=record.protocol,
                patient_dose_msv=record.effective_dose_msv,
                drl_msv=0.0,
                achievable_dose_msv=0.0,
                ratio=0.0,
                status="no_drl_available",
                optimization_suggestions=[],
            )

        drl_msv = drl_entry["drl"]
        achievable = drl_entry["achievable"]
        ratio = record.effective_dose_msv / drl_msv if drl_msv > 0 else 0.0

        if record.effective_dose_msv <= achievable:
            status = "below_achievable"
        elif record.effective_dose_msv <= drl_msv:
            status = "below_drl"
        elif ratio <= 1.5:
            status = "above_drl"
        else:
            status = "significantly_above"

        suggestions = self._get_optimization_suggestions(record.protocol, ratio, record.pediatric)

        return DRLComparison(
            protocol=record.protocol,
            patient_dose_msv=round(record.effective_dose_msv, 2),
            drl_msv=drl_msv,
            achievable_dose_msv=achievable,
            ratio=round(ratio, 2),
            status=status,
            optimization_suggestions=suggestions,
        )

    def get_optimization_suggestions(self, record: DoseRecord) -> List[str]:
        """Get optimization suggestions for a given dose record."""
        drl_entry = self._find_drl(record.protocol)
        if drl_entry is None:
            return ["No DRL available for this protocol — review local benchmarks."]
        ratio = record.effective_dose_msv / drl_entry["drl"] if drl_entry["drl"] > 0 else 0.0
        return self._get_optimization_suggestions(record.protocol, ratio, record.pediatric)

    # ───────────────────────────────────────────────────────────────────
    # POPULATION ANALYTICS
    # ───────────────────────────────────────────────────────────────────

    def get_population_dose_summary(self) -> Dict[str, Any]:
        """Generate population-level dose statistics from the registry.

        Returns
        -------
        dict
            Summary statistics including total records, unique patients,
            per-protocol and per-modality breakdowns, and alert distribution.
        """
        if not self._dose_registry:
            return {
                "total_records": 0,
                "unique_patients": 0,
                "mean_dose_msv": 0.0,
                "median_dose_msv": 0.0,
                "by_protocol": {},
                "by_modality": {},
                "alert_distribution": {},
            }

        doses = [r.effective_dose_msv for r in self._dose_registry]
        sorted_doses = sorted(doses)
        n = len(sorted_doses)
        median = sorted_doses[n // 2] if n % 2 == 1 else (sorted_doses[n // 2 - 1] + sorted_doses[n // 2]) / 2

        # Per-protocol stats
        by_protocol: Dict[str, List[float]] = {}
        for r in self._dose_registry:
            by_protocol.setdefault(r.protocol, []).append(r.effective_dose_msv)

        protocol_summary = {}
        for proto, proto_doses in by_protocol.items():
            protocol_summary[proto] = {
                "count": len(proto_doses),
                "mean_msv": round(sum(proto_doses) / len(proto_doses), 2),
                "min_msv": round(min(proto_doses), 2),
                "max_msv": round(max(proto_doses), 2),
            }

        # Per-modality stats
        by_modality: Dict[str, List[float]] = {}
        for r in self._dose_registry:
            by_modality.setdefault(r.modality, []).append(r.effective_dose_msv)

        modality_summary = {}
        for mod, mod_doses in by_modality.items():
            modality_summary[mod] = {
                "count": len(mod_doses),
                "mean_msv": round(sum(mod_doses) / len(mod_doses), 2),
                "total_msv": round(sum(mod_doses), 2),
            }

        # Alert distribution across unique patients
        patient_ids = set(r.patient_id for r in self._dose_registry)
        alert_counts: Dict[str, int] = {"normal": 0, "elevated": 0, "high": 0, "critical": 0}
        for pid in patient_ids:
            cumulative = self.get_cumulative_dose(pid)
            alert_counts[cumulative.alert_level] = alert_counts.get(cumulative.alert_level, 0) + 1

        return {
            "total_records": n,
            "unique_patients": len(patient_ids),
            "mean_dose_msv": round(sum(doses) / n, 2),
            "median_dose_msv": round(median, 2),
            "by_protocol": protocol_summary,
            "by_modality": modality_summary,
            "alert_distribution": alert_counts,
        }

    # ───────────────────────────────────────────────────────────────────
    # DEMO DATA GENERATION
    # ───────────────────────────────────────────────────────────────────

    def generate_demo_data(self, n_records: int = 200) -> None:
        """Generate synthetic dose records for demonstration purposes.

        Parameters
        ----------
        n_records : int
            Number of synthetic records to generate (default 200).
        """
        rng = random.Random(42)
        n_patients = max(n_records // 5, 10)
        patient_ids = [f"PAT-{i:05d}" for i in range(1, n_patients + 1)]
        scanner_models = [
            "Siemens SOMATOM Force",
            "GE Revolution CT",
            "Philips iQon Spectral CT",
            "Canon Aquilion ONE",
            "Siemens SOMATOM X.cite",
        ]

        base_date = datetime.now() - timedelta(days=365)

        for _ in range(n_records):
            proto = rng.choice(self._DEMO_PROTOCOLS)
            patient_id = rng.choice(patient_ids)
            study_date = base_date + timedelta(days=rng.randint(0, 365))
            dose = round(rng.uniform(*proto["dose_range"]), 2)
            pediatric = rng.random() < 0.1  # 10% pediatric

            if pediatric:
                dose *= 0.5

            # Generate plausible DLP and CTDIvol for CT
            dlp = None
            ctdi = None
            if proto["modality"] == "CT":
                # Rough conversion: effective dose (mSv) ~ DLP * k, where k ~ 0.014-0.017
                dlp = round(dose / 0.015, 1)
                ctdi = round(dlp / rng.uniform(25, 40), 1)

            record = DoseRecord(
                patient_id=patient_id,
                study_date=study_date.strftime("%Y-%m-%d"),
                modality=proto["modality"],
                protocol=proto["protocol"],
                body_region=proto["region"],
                effective_dose_msv=dose,
                dlp_mgy_cm=dlp,
                ctdi_vol_mgy=ctdi,
                num_exposures=rng.randint(1, 3) if proto["modality"] == "XR" else 1,
                scanner_model=rng.choice(scanner_models) if proto["modality"] == "CT" else None,
                pediatric=pediatric,
            )
            self.record_dose(record)

        logger.info(
            "Generated {} demo dose records for {} patients",
            n_records,
            n_patients,
        )

    # ───────────────────────────────────────────────────────────────────
    # INTERNAL HELPERS
    # ───────────────────────────────────────────────────────────────────

    def _classify_alert_level(
        self,
        total_msv: float,
        pediatric: bool = False,
    ) -> Tuple[str, Optional[str]]:
        """Classify cumulative dose into alert levels.

        Parameters
        ----------
        total_msv : float
            Total effective dose in mSv.
        pediatric : bool
            If True, apply pediatric multiplier (halve thresholds).

        Returns
        -------
        tuple[str, str | None]
            (alert_level, alert_message)
        """
        multiplier = self.PEDIATRIC_MULTIPLIER if pediatric else 1.0
        thresholds = {k: v * multiplier for k, v in self.ANNUAL_THRESHOLDS.items()}

        if total_msv >= thresholds["high"]:
            return (
                "critical",
                f"CRITICAL: Cumulative dose {total_msv:.1f} mSv exceeds {thresholds['high']:.0f} mSv threshold. "
                "Immediate review of imaging utilization required. "
                "Consider alternative non-ionizing modalities."
                + (" Pediatric thresholds applied." if pediatric else ""),
            )
        elif total_msv >= thresholds["elevated"]:
            return (
                "high",
                f"HIGH: Cumulative dose {total_msv:.1f} mSv exceeds {thresholds['elevated']:.0f} mSv threshold. "
                "Review imaging frequency and consider dose optimization."
                + (" Pediatric thresholds applied." if pediatric else ""),
            )
        elif total_msv >= thresholds["normal"]:
            return (
                "elevated",
                f"ELEVATED: Cumulative dose {total_msv:.1f} mSv approaching annual limit. "
                "Monitor ongoing imaging utilization."
                + (" Pediatric thresholds applied." if pediatric else ""),
            )
        else:
            return ("normal", None)

    def _find_drl(self, protocol: str) -> Optional[Dict[str, float]]:
        """Find the matching DRL entry for a protocol name."""
        # Exact match first
        if protocol in self.NATIONAL_DRLS:
            return self.NATIONAL_DRLS[protocol]

        # Partial match
        protocol_lower = protocol.lower()
        for key, value in self.NATIONAL_DRLS.items():
            if key.lower() in protocol_lower or protocol_lower in key.lower():
                return value

        return None

    def _get_optimization_suggestions(
        self,
        protocol: str,
        ratio: float,
        pediatric: bool = False,
    ) -> List[str]:
        """Generate dose optimization suggestions based on DRL ratio.

        Parameters
        ----------
        protocol : str
            Protocol name.
        ratio : float
            Patient dose / DRL ratio.
        pediatric : bool
            Whether the patient is pediatric.

        Returns
        -------
        list[str]
            Practical suggestions for dose reduction.
        """
        suggestions: List[str] = []

        if ratio > 1.5:
            suggestions.extend([
                "Consider reducing kV from 120 to 100 kV (up to 40% dose reduction).",
                "Enable iterative reconstruction (ADMIRE, ASIR-V, IMR, or AiCE).",
                "Review scan range — reduce to clinically necessary anatomy only.",
                "Verify automatic exposure control (AEC) is active and properly calibrated.",
                "Consider reducing number of scan phases if multiphase protocol.",
            ])
        elif ratio > 1.0:
            suggestions.extend([
                "Review mAs settings — current dose exceeds DRL.",
                "Consider automatic exposure control (AEC) optimization.",
                "Evaluate if iterative or deep-learning reconstruction could maintain "
                "image quality at lower dose.",
            ])
        elif ratio > 0.75:
            suggestions.append(
                "Dose is within acceptable range but above achievable level — "
                "consider incremental optimization."
            )
        else:
            suggestions.append(
                "Dose is at or below achievable level — current protocol is well optimized."
            )

        if pediatric:
            suggestions.extend([
                "PEDIATRIC: Use weight-based or age-based pediatric protocol.",
                "PEDIATRIC: Reduce kV to 80-100 kV for smaller body habitus.",
                "PEDIATRIC: Apply size-specific dose estimates (SSDE) rather than CTDIvol.",
                "PEDIATRIC: Consider alternative non-ionizing modalities (US, MRI) where appropriate.",
            ])

        return suggestions
