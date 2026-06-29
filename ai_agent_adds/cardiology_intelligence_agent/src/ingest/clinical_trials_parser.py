"""ClinicalTrials.gov cardiovascular trials ingest parser.

Fetches cardiovascular clinical trials via the ClinicalTrials.gov API v2,
parses JSON responses into IngestRecord instances, and targets the
``cardio_trials`` Milvus collection.

Includes a curated dictionary of 25+ landmark cardiovascular trials
(PARADIGM-HF, DAPA-HF, EMPEROR-Reduced, ISCHEMIA, etc.) that are
seeded regardless of API availability.

API v2 docs: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: March 2026
"""

import re
import time
from typing import Any, Dict, List, Optional

import requests

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

CT_GOV_BASE_URL = "https://clinicaltrials.gov/api/v2"


# ═══════════════════════════════════════════════════════════════════════
# CLINICAL TRIALS CARDIO PARSER
# ═══════════════════════════════════════════════════════════════════════


class ClinicalTrialsCardioParser(BaseIngestParser):
    """Ingest parser for ClinicalTrials.gov cardiovascular trials.

    Searches across nine cardiovascular condition categories, retrieves
    trial metadata via the ClinicalTrials.gov API v2, and produces
    IngestRecord instances for the ``cardio_trials`` collection.

    Also provides :meth:`seed_landmark_trials` to inject curated records
    for 25+ landmark cardiovascular trials with known outcomes.

    Usage::

        parser = ClinicalTrialsCardioParser()
        records = parser.run(max_results=500)
        landmark = parser.seed_landmark_trials()
    """

    CARDIO_CONDITIONS: List[str] = [
        "Heart Failure",
        "Coronary Artery Disease",
        "Atrial Fibrillation",
        "Myocardial Infarction",
        "Cardiomyopathy",
        "Valvular Heart Disease",
        "Pulmonary Hypertension",
        "Aortic Stenosis",
        "Heart Transplant",
    ]
    """Cardiovascular conditions used for the API search query."""

    LANDMARK_TRIALS: Dict[str, Dict[str, Any]] = {
        # ── Heart Failure ─────────────────────────────────────────
        "PARADIGM-HF": {
            "nct": "NCT01035255",
            "condition": "Heart Failure",
            "drug": "Sacubitril/Valsartan",
            "phase": "Phase 3",
            "enrollment": 8442,
            "status": "Completed",
            "key_findings": (
                "Sacubitril/valsartan reduced cardiovascular death or "
                "HF hospitalization by 20% vs enalapril (HR 0.80, "
                "p<0.001). Established ARNI as foundational HFrEF therapy."
            ),
        },
        "DAPA-HF": {
            "nct": "NCT03036124",
            "condition": "Heart Failure",
            "drug": "Dapagliflozin",
            "phase": "Phase 3",
            "enrollment": 4744,
            "status": "Completed",
            "key_findings": (
                "Dapagliflozin reduced composite of worsening HF or CV "
                "death by 26% vs placebo (HR 0.74, p<0.001) in HFrEF, "
                "regardless of diabetes status."
            ),
        },
        "EMPEROR-Reduced": {
            "nct": "NCT03057977",
            "condition": "Heart Failure",
            "drug": "Empagliflozin",
            "phase": "Phase 3",
            "enrollment": 3730,
            "status": "Completed",
            "key_findings": (
                "Empagliflozin reduced CV death or HF hospitalization "
                "by 25% (HR 0.75, p<0.001) in HFrEF. Confirmed SGLT2i "
                "class effect in heart failure."
            ),
        },
        "EMPEROR-Preserved": {
            "nct": "NCT03057951",
            "condition": "Heart Failure",
            "drug": "Empagliflozin",
            "phase": "Phase 3",
            "enrollment": 5988,
            "status": "Completed",
            "key_findings": (
                "Empagliflozin reduced CV death or HF hospitalization "
                "by 21% (HR 0.79, p<0.001) in HFpEF. First drug to "
                "show benefit in HFpEF."
            ),
        },
        "DELIVER": {
            "nct": "NCT03619213",
            "condition": "Heart Failure",
            "drug": "Dapagliflozin",
            "phase": "Phase 3",
            "enrollment": 6263,
            "status": "Completed",
            "key_findings": (
                "Dapagliflozin reduced worsening HF or CV death by 18% "
                "(HR 0.82, p<0.001) in HFmrEF/HFpEF."
            ),
        },
        "RALES": {
            "nct": "NCT00000536",
            "condition": "Heart Failure",
            "drug": "Spironolactone",
            "phase": "Phase 3",
            "enrollment": 1663,
            "status": "Completed",
            "key_findings": (
                "Spironolactone reduced all-cause mortality by 30% "
                "(RR 0.70, p<0.001) in severe HFrEF. Established MRA "
                "as core GDMT pillar."
            ),
        },
        "COPERNICUS": {
            "nct": "NCT00000613",
            "condition": "Heart Failure",
            "drug": "Carvedilol",
            "phase": "Phase 3",
            "enrollment": 2289,
            "status": "Completed",
            "key_findings": (
                "Carvedilol reduced all-cause mortality by 35% "
                "(p=0.0014) in severe HFrEF (LVEF <25%)."
            ),
        },
        # ── Coronary Artery Disease ───────────────────────────────
        "ISCHEMIA": {
            "nct": "NCT01471522",
            "condition": "Stable Coronary Artery Disease",
            "drug": "Invasive strategy vs conservative",
            "phase": "Phase 4",
            "enrollment": 5179,
            "status": "Completed",
            "key_findings": (
                "Initial invasive strategy did not reduce death, MI, "
                "or hospitalization for angina vs conservative medical "
                "therapy in stable CAD with moderate-severe ischemia."
            ),
        },
        "COURAGE": {
            "nct": "NCT00007657",
            "condition": "Stable Coronary Artery Disease",
            "drug": "PCI + OMT vs OMT alone",
            "phase": "Phase 4",
            "enrollment": 2287,
            "status": "Completed",
            "key_findings": (
                "PCI added to optimal medical therapy did not reduce "
                "death, MI, or major CV events vs OMT alone in stable CAD."
            ),
        },
        "SYNTAX": {
            "nct": "NCT00114972",
            "condition": "Complex Coronary Artery Disease",
            "drug": "PCI vs CABG",
            "phase": "Phase 4",
            "enrollment": 1800,
            "status": "Completed",
            "key_findings": (
                "In three-vessel or left main CAD, CABG had lower "
                "MACCE rates than PCI at 5 years. SYNTAX score guides "
                "revascularization strategy selection."
            ),
        },
        # ── Atrial Fibrillation ───────────────────────────────────
        "AFFIRM": {
            "nct": "NCT00000556",
            "condition": "Atrial Fibrillation",
            "drug": "Rate control vs rhythm control",
            "phase": "Phase 3",
            "enrollment": 4060,
            "status": "Completed",
            "key_findings": (
                "Rate control was non-inferior to rhythm control for "
                "mortality in AF. Rhythm control did not reduce stroke "
                "or mortality."
            ),
        },
        "CASTLE-AF": {
            "nct": "NCT00643188",
            "condition": "Atrial Fibrillation with Heart Failure",
            "drug": "Catheter ablation",
            "phase": "Phase 4",
            "enrollment": 363,
            "status": "Completed",
            "key_findings": (
                "AF catheter ablation reduced composite of death or HF "
                "hospitalization by 38% (HR 0.62, p=0.007) in patients "
                "with AF and HFrEF."
            ),
        },
        "EAST-AFNET 4": {
            "nct": "NCT01288352",
            "condition": "Atrial Fibrillation",
            "drug": "Early rhythm control",
            "phase": "Phase 4",
            "enrollment": 2789,
            "status": "Completed",
            "key_findings": (
                "Early systematic rhythm control reduced CV death, "
                "stroke, and HF hospitalization by 21% (HR 0.79, "
                "p=0.005) in recently diagnosed AF."
            ),
        },
        # ── Valvular Heart Disease ────────────────────────────────
        "PARTNER 3": {
            "nct": "NCT02675114",
            "condition": "Aortic Stenosis",
            "drug": "TAVR (SAPIEN 3) vs SAVR",
            "phase": "Phase 3",
            "enrollment": 1000,
            "status": "Completed",
            "key_findings": (
                "TAVR was superior to SAVR for the composite of death, "
                "stroke, or rehospitalization at 1 year in low-risk "
                "patients with severe aortic stenosis."
            ),
        },
        "Evolut Low Risk": {
            "nct": "NCT02701283",
            "condition": "Aortic Stenosis",
            "drug": "TAVR (CoreValve Evolut) vs SAVR",
            "phase": "Phase 3",
            "enrollment": 1468,
            "status": "Completed",
            "key_findings": (
                "Self-expanding TAVR was non-inferior to surgery for "
                "death or disabling stroke at 24 months in low-risk "
                "severe AS patients."
            ),
        },
        "COAPT": {
            "nct": "NCT01626079",
            "condition": "Mitral Regurgitation with Heart Failure",
            "drug": "MitraClip (TEER)",
            "phase": "Phase 4",
            "enrollment": 614,
            "status": "Completed",
            "key_findings": (
                "MitraClip TEER reduced HF hospitalization by 47% and "
                "all-cause mortality by 38% at 2 years in severe "
                "secondary MR with HF."
            ),
        },
        # ── Prevention ────────────────────────────────────────────
        "FOURIER": {
            "nct": "NCT01764633",
            "condition": "Atherosclerotic CVD",
            "drug": "Evolocumab (PCSK9 inhibitor)",
            "phase": "Phase 3",
            "enrollment": 27564,
            "status": "Completed",
            "key_findings": (
                "Evolocumab reduced LDL by 59% and major CV events by "
                "15% (HR 0.85, p<0.001) in patients with established "
                "atherosclerotic disease."
            ),
        },
        "REDUCE-IT": {
            "nct": "NCT01492361",
            "condition": "Hypertriglyceridemia + ASCVD",
            "drug": "Icosapent ethyl",
            "phase": "Phase 3",
            "enrollment": 8179,
            "status": "Completed",
            "key_findings": (
                "Icosapent ethyl reduced major adverse CV events by "
                "25% (HR 0.75, p<0.001) in statin-treated patients "
                "with elevated triglycerides."
            ),
        },
        "CLEAR Outcomes": {
            "nct": "NCT02993406",
            "condition": "Statin-intolerant hypercholesterolemia",
            "drug": "Bempedoic acid",
            "phase": "Phase 3",
            "enrollment": 13970,
            "status": "Completed",
            "key_findings": (
                "Bempedoic acid reduced major adverse CV events by 13% "
                "(HR 0.87, p=0.004) in statin-intolerant patients."
            ),
        },
        # ── Pulmonary Hypertension ────────────────────────────────
        "GRIPHON": {
            "nct": "NCT01106014",
            "condition": "Pulmonary Arterial Hypertension",
            "drug": "Selexipag",
            "phase": "Phase 3",
            "enrollment": 1156,
            "status": "Completed",
            "key_findings": (
                "Selexipag reduced the composite of death or PAH "
                "complications by 40% (HR 0.60, p<0.001)."
            ),
        },
        # ── Myocardial Infarction ─────────────────────────────────
        "PLATO": {
            "nct": "NCT00391872",
            "condition": "Acute Coronary Syndromes",
            "drug": "Ticagrelor vs clopidogrel",
            "phase": "Phase 3",
            "enrollment": 18624,
            "status": "Completed",
            "key_findings": (
                "Ticagrelor reduced CV death, MI, or stroke by 16% "
                "(HR 0.84, p<0.001) vs clopidogrel in ACS without "
                "increase in overall major bleeding."
            ),
        },
        "COMPLETE": {
            "nct": "NCT01740479",
            "condition": "STEMI with Multivessel Disease",
            "drug": "Complete revascularization vs culprit-only PCI",
            "phase": "Phase 4",
            "enrollment": 4041,
            "status": "Completed",
            "key_findings": (
                "Routine complete revascularization reduced CV death "
                "or MI by 26% (HR 0.74, p=0.004) vs culprit-only PCI "
                "in STEMI with multivessel CAD."
            ),
        },
        # ── Cardiomyopathy ────────────────────────────────────────
        "MAVERICK-HCM": {
            "nct": "NCT03470545",
            "condition": "Hypertrophic Cardiomyopathy",
            "drug": "Mavacamten",
            "phase": "Phase 2",
            "enrollment": 59,
            "status": "Completed",
            "key_findings": (
                "Mavacamten reduced NT-proBNP and cardiac troponin I "
                "in obstructive HCM, supporting Phase 3 development."
            ),
        },
        "EXPLORER-HCM": {
            "nct": "NCT03470545",
            "condition": "Obstructive Hypertrophic Cardiomyopathy",
            "drug": "Mavacamten",
            "phase": "Phase 3",
            "enrollment": 251,
            "status": "Completed",
            "key_findings": (
                "Mavacamten improved exercise capacity and NYHA class "
                "in obstructive HCM. 37% achieved complete response "
                "vs 17% placebo (p=0.0005). FDA-approved as Camzyos."
            ),
        },
        # ── Heart Failure (continued) ────────────────────────────
        "STRONG-HF": {
            "nct": "NCT03412201",
            "condition": "Heart Failure",
            "drug": "Rapid GDMT uptitration (intensive vs usual care)",
            "phase": "Phase 4",
            "enrollment": 1078,
            "status": "Completed",
            "key_findings": (
                "Rapid uptitration of GDMT (within 2 weeks of discharge) "
                "reduced 180-day HF death or readmission by 34% (HR 0.66, "
                "p=0.002). Supports aggressive, early optimization of "
                "neurohormonal therapy after HF hospitalization."
            ),
        },
        "EMPULSE": {
            "nct": "NCT04157751",
            "condition": "Acute Decompensated Heart Failure",
            "drug": "Empagliflozin",
            "phase": "Phase 3",
            "enrollment": 530,
            "status": "Completed",
            "key_findings": (
                "Empagliflozin initiated in-hospital for acute decompensated "
                "HF showed significant clinical benefit at 90 days (win ratio "
                "1.36, p=0.0054) including death, HF events, change in KCCQ, "
                "and weight loss. Supports in-hospital SGLT2i initiation."
            ),
        },
        "ADVOR": {
            "nct": "NCT03505788",
            "condition": "Acute Decompensated Heart Failure",
            "drug": "Acetazolamide + loop diuretics",
            "phase": "Phase 3",
            "enrollment": 519,
            "status": "Completed",
            "key_findings": (
                "Acetazolamide 500 mg IV added to standardized loop diuretic "
                "therapy improved successful decongestion at 3 days (42.2% "
                "vs 30.5%, RR 1.46, p<0.001) in acute decompensated HF. "
                "Supports combination diuretic strategy."
            ),
        },
        "DAPA-CKD": {
            "nct": "NCT03036150",
            "condition": "Chronic Kidney Disease",
            "drug": "Dapagliflozin",
            "phase": "Phase 3",
            "enrollment": 4304,
            "status": "Completed",
            "key_findings": (
                "Dapagliflozin reduced the kidney composite endpoint "
                "(sustained eGFR decline >=50%, ESKD, renal/CV death) by "
                "39% (HR 0.61, p<0.001) in CKD with or without diabetes. "
                "Also reduced CV death and HF hospitalization."
            ),
        },
        # ── Prevention / Lipid ──────────────────────────────────
        "SELECT": {
            "nct": "NCT03574597",
            "condition": "Cardiovascular Disease in Obesity",
            "drug": "Semaglutide 2.4 mg weekly",
            "phase": "Phase 3",
            "enrollment": 17604,
            "status": "Completed",
            "key_findings": (
                "Semaglutide 2.4 mg weekly reduced MACE (CV death, non-fatal "
                "MI, non-fatal stroke) by 20% (HR 0.80, p<0.001) in "
                "overweight/obese patients with established CVD but without "
                "diabetes. Practice-changing for cardiometabolic risk reduction."
            ),
        },
        "PROMINENT": {
            "nct": "NCT03071692",
            "condition": "Hypertriglyceridemia with Type 2 Diabetes",
            "drug": "Pemafibrate",
            "phase": "Phase 3",
            "enrollment": 10497,
            "status": "Completed",
            "key_findings": (
                "Pemafibrate reduced triglycerides by 26% but did NOT reduce "
                "cardiovascular events (HR 1.03, p=0.67) in patients with "
                "T2DM and elevated triglycerides. Negative trial; demonstrates "
                "that TG lowering alone does not improve CV outcomes."
            ),
        },
        # ── Coronary (continued) ────────────────────────────────
        "ORBITA-2": {
            "nct": "NCT03742050",
            "condition": "Stable Angina",
            "drug": "PCI vs placebo procedure",
            "phase": "Phase 4",
            "enrollment": 301,
            "status": "Completed",
            "key_findings": (
                "PCI provided significant symptom improvement over placebo "
                "procedure in patients with stable angina (angina symptom "
                "score difference -2.9, p<0.001) when antianginal medications "
                "were withdrawn. Clarifies that PCI provides genuine "
                "anti-ischemic benefit beyond placebo effect."
            ),
        },
        "REVIVED-BCIS2": {
            "nct": "NCT01920048",
            "condition": "Ischemic Cardiomyopathy",
            "drug": "PCI + OMT vs OMT alone",
            "phase": "Phase 4",
            "enrollment": 700,
            "status": "Completed",
            "key_findings": (
                "PCI plus optimal medical therapy did NOT reduce all-cause "
                "death or HF hospitalization (HR 0.99, p=0.96) compared "
                "with OMT alone in patients with severe ischemic "
                "cardiomyopathy (LVEF <=35%) and extensive CAD. Negative, "
                "practice-changing trial supporting OMT-first approach."
            ),
        },
        # ── Arrhythmia (continued) ──────────────────────────────
        "ADVENT": {
            "nct": "NCT04612244",
            "condition": "Paroxysmal Atrial Fibrillation",
            "drug": "Pulsed field ablation (PFA)",
            "phase": "Phase 3",
            "enrollment": 607,
            "status": "Completed",
            "key_findings": (
                "Pulsed field ablation was non-inferior to conventional "
                "thermal ablation (RF or cryo) for the primary effectiveness "
                "endpoint (freedom from AF/AT/AFL at 1 year: 73.3% vs "
                "71.3%). PFA showed tissue selectivity with lower risk of "
                "collateral damage to esophagus and phrenic nerve."
            ),
        },
        "OPTION": {
            "nct": "NCT03795298",
            "condition": "Atrial Fibrillation",
            "drug": "LAA occlusion (Watchman FLX) vs DOAC",
            "phase": "Phase 4",
            "enrollment": 1600,
            "status": "Completed",
            "key_findings": (
                "Left atrial appendage occlusion with Watchman FLX was "
                "non-inferior to DOAC therapy for the composite of stroke, "
                "systemic embolism, or CV/unexplained death in patients "
                "with AF undergoing catheter ablation. Supports LAAO as "
                "alternative stroke prevention strategy."
            ),
        },
        # ── Structural (continued) ──────────────────────────────
        "TRILUMINATE Pivotal": {
            "nct": "NCT04482907",
            "condition": "Severe Tricuspid Regurgitation",
            "drug": "Tricuspid TEER (TriClip)",
            "phase": "Phase 3",
            "enrollment": 350,
            "status": "Completed",
            "key_findings": (
                "Tricuspid TEER (TriClip) met safety endpoint and reduced "
                "TR severity (87% achieved moderate or less TR vs 4.8% "
                "control). Quality of life (KCCQ) improved significantly. "
                "Did not meet composite primary endpoint of death/tricuspid "
                "surgery/HF hospitalization, but FDA approved based on TR "
                "reduction and symptom improvement."
            ),
        },
        "CLASP IID": {
            "nct": "NCT03706833",
            "condition": "Degenerative Mitral Regurgitation",
            "drug": "Mitral TEER (PASCAL) vs MitraClip",
            "phase": "Phase 3",
            "enrollment": 180,
            "status": "Completed",
            "key_findings": (
                "PASCAL mitral TEER was non-inferior to MitraClip for the "
                "composite of CV mortality, mitral valve surgery, or MR "
                ">=3+ at 2 years in degenerative MR (DMR). Comparable "
                "safety and effectiveness profile. FDA approved as "
                "alternative TEER device."
            ),
        },
        # ── Cardiomyopathy (continued) ──────────────────────────
        "VALOR-HCM": {
            "nct": "NCT04349072",
            "condition": "Obstructive Hypertrophic Cardiomyopathy",
            "drug": "Mavacamten",
            "phase": "Phase 3",
            "enrollment": 112,
            "status": "Completed",
            "key_findings": (
                "Mavacamten markedly reduced the proportion of patients "
                "meeting guideline criteria for septal reduction therapy "
                "(17.9% vs 76.8% placebo at 16 weeks, p<0.001). Established "
                "mavacamten as medical alternative to surgical myectomy/"
                "alcohol septal ablation in eligible patients."
            ),
        },
        "ATTR-ACT": {
            "nct": "NCT01994889",
            "condition": "Transthyretin Amyloid Cardiomyopathy",
            "drug": "Tafamidis",
            "phase": "Phase 3",
            "enrollment": 441,
            "status": "Completed",
            "key_findings": (
                "Tafamidis reduced all-cause mortality by 30% (HR 0.70, "
                "p=0.026) and cardiovascular hospitalizations by 32% in "
                "ATTR cardiomyopathy (wild-type and hereditary). First "
                "and only disease-modifying therapy approved for ATTR-CM."
            ),
        },
        # ── Heart Transplant ──────────────────────────────────────
        "MOMENTUM 3": {
            "nct": "NCT02224755",
            "condition": "Advanced Heart Failure",
            "drug": "HeartMate 3 LVAD",
            "phase": "Phase 3",
            "enrollment": 1028,
            "status": "Completed",
            "key_findings": (
                "HeartMate 3 magnetically levitated centrifugal LVAD "
                "was superior to HeartMate II for survival free of "
                "disabling stroke or reoperation at 2 years (77.9% "
                "vs 56.4%, p<0.001)."
            ),
        },
    }
    """Curated landmark cardiovascular trials with known outcomes."""

    def __init__(self, base_url: str = CT_GOV_BASE_URL):
        """Initialize the ClinicalTrials.gov cardiology parser.

        Args:
            base_url: ClinicalTrials.gov API v2 base URL.
        """
        super().__init__("cardio_trials")
        self.base_url = base_url

    # ── Fetch ─────────────────────────────────────────────────────────

    def fetch(self, max_results: int = 500, page_size: int = 100) -> List[dict]:
        """Fetch cardiovascular trials from ClinicalTrials.gov API v2.

        Searches across all nine :attr:`CARDIO_CONDITIONS` using the
        ``query.cond`` parameter with pagination.

        Args:
            max_results: Maximum total trials to retrieve.
            page_size: Trials per API request (max 1000).

        Returns:
            List of study JSON objects from the API response.
        """
        query = " OR ".join(
            f'"{cond}"' for cond in self.CARDIO_CONDITIONS
        )
        url = f"{self.base_url}/studies"
        all_studies: List[dict] = []
        page_token: Optional[str] = None
        page_num = 0

        while len(all_studies) < max_results:
            params: Dict[str, Any] = {
                "query.cond": query,
                "pageSize": min(page_size, max_results - len(all_studies)),
                "sort": "LastUpdatePostDate",
            }
            if page_token:
                params["pageToken"] = page_token

            try:
                response = requests.get(
                    url,
                    params=params,
                    headers={"Accept": "application/json"},
                    timeout=30,
                )
                response.raise_for_status()
                data = response.json()
            except requests.RequestException as exc:
                self.logger.error(
                    f"ClinicalTrials.gov API error on page {page_num + 1}: {exc}"
                )
                break

            studies = data.get("studies", [])
            all_studies.extend(studies)
            page_num += 1
            self.logger.info(
                f"Fetched page {page_num}: {len(studies)} studies "
                f"(total {len(all_studies)})"
            )

            page_token = data.get("nextPageToken")
            if not page_token:
                break

            # Rate limit: 1 request/second
            time.sleep(1)

        return all_studies[:max_results]

    # ── Parse ─────────────────────────────────────────────────────────

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Parse ClinicalTrials.gov JSON studies into IngestRecord instances.

        Extracts fields from the API v2 JSON structure and maps them
        to the ``cardio_trials`` collection schema.

        Args:
            raw_data: List of study JSON objects from the API.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []

        for study in raw_data:
            try:
                record = self._parse_trial(study)
                if record:
                    records.append(record)
            except Exception as exc:
                nct = (
                    study.get("protocolSection", {})
                    .get("identificationModule", {})
                    .get("nctId", "unknown")
                )
                self.logger.warning(f"Failed to parse trial {nct}: {exc}")
                continue

        self.logger.info(
            f"Parsed {len(records)} trial records from "
            f"{len(raw_data)} API studies"
        )
        return records

    # ── Landmark trial seeding ────────────────────────────────────────

    def seed_landmark_trials(self) -> List[IngestRecord]:
        """Generate IngestRecord instances for curated landmark trials.

        Creates records from the :attr:`LANDMARK_TRIALS` dictionary,
        providing pre-populated outcomes data for the most important
        cardiovascular trials that every cardiologist should know.

        Returns:
            List of :class:`IngestRecord` instances for landmark trials.
        """
        records: List[IngestRecord] = []

        for trial_name, info in self.LANDMARK_TRIALS.items():
            condition = info.get("condition", "")
            drug = info.get("drug", "")
            findings = info.get("key_findings", "")

            text = (
                f"{trial_name}: {condition}. "
                f"Intervention: {drug}. "
                f"{findings}"
            )

            metadata = {
                "trial_name": trial_name,
                "nct_id": info.get("nct", ""),
                "phase": info.get("phase", ""),
                "condition": condition,
                "intervention": drug,
                "primary_outcome": "",
                "enrollment": info.get("enrollment", 0),
                "status": info.get("status", "Completed"),
                "start_date": "",
                "completion_date": "",
                "sponsor": "",
                "key_findings": self.truncate(findings, 4094),
            }

            records.append(
                IngestRecord(
                    text=self.truncate(text, 4096),
                    metadata=metadata,
                    collection=self.collection,
                    source="Landmark Trial Seed",
                    source_id=info.get("nct"),
                )
            )

        self.logger.info(
            f"Seeded {len(records)} landmark cardiovascular trials"
        )
        return records

    # ── Private helpers ───────────────────────────────────────────────

    def _parse_trial(self, trial_json: dict) -> Optional[IngestRecord]:
        """Parse a single ClinicalTrials.gov study JSON into an IngestRecord.

        Navigates the API v2 JSON structure to extract:
          - protocolSection.identificationModule (NCT ID, title)
          - protocolSection.descriptionModule (brief summary)
          - protocolSection.designModule (phase, enrollment)
          - protocolSection.statusModule (status, dates)
          - protocolSection.sponsorCollaboratorsModule (sponsor)
          - protocolSection.conditionsModule (conditions)
          - protocolSection.armsInterventionsModule (interventions)
          - protocolSection.outcomesModule (primary outcomes)

        Args:
            trial_json: Single study JSON object from the API.

        Returns:
            An :class:`IngestRecord`, or ``None`` if required fields
            are missing.
        """
        protocol = trial_json.get("protocolSection", {})
        id_mod = protocol.get("identificationModule", {})
        desc_mod = protocol.get("descriptionModule", {})
        design_mod = protocol.get("designModule", {})
        status_mod = protocol.get("statusModule", {})
        sponsor_mod = protocol.get("sponsorCollaboratorsModule", {})
        cond_mod = protocol.get("conditionsModule", {})
        arms_mod = protocol.get("armsInterventionsModule", {})
        outcomes_mod = protocol.get("outcomesModule", {})

        # Required: NCT ID
        nct_id = id_mod.get("nctId", "")
        if not nct_id:
            return None

        # Title
        title = (
            id_mod.get("officialTitle")
            or id_mod.get("briefTitle")
            or "Untitled"
        )
        brief_summary = desc_mod.get("briefSummary", "")

        # Phase
        phases = design_mod.get("phases", [])
        phase = ", ".join(phases) if phases else ""

        # Enrollment
        enrollment_info = design_mod.get("enrollmentInfo", {})
        enrollment = enrollment_info.get("count", 0) or 0

        # Status
        status = status_mod.get("overallStatus", "")

        # Dates
        start_struct = status_mod.get("startDateStruct", {})
        start_date = start_struct.get("date", "")
        completion_struct = status_mod.get("completionDateStruct", {})
        completion_date = completion_struct.get("date", "")

        # Sponsor
        lead_sponsor = sponsor_mod.get("leadSponsor", {})
        sponsor = lead_sponsor.get("name", "")

        # Conditions
        conditions = cond_mod.get("conditions", [])
        condition_str = "; ".join(conditions[:5])

        # Interventions
        interventions = arms_mod.get("interventions", [])
        intervention_names = [
            intv.get("name", "") for intv in interventions if intv.get("name")
        ]
        intervention_str = "; ".join(intervention_names[:3])

        # Primary outcomes
        primary_outcomes = outcomes_mod.get("primaryOutcomes", [])
        outcome_measures = [
            po.get("measure", "") for po in primary_outcomes if po.get("measure")
        ]
        primary_outcome_str = "; ".join(outcome_measures[:3])

        # Build embedding text
        text = (
            f"{title}. {brief_summary}. "
            f"Conditions: {condition_str}. "
            f"Interventions: {intervention_str}."
        )

        # Try to extract trial acronym from title
        trial_name = self._extract_trial_name(title)

        metadata = {
            "trial_name": self.truncate(trial_name or title, 254),
            "nct_id": nct_id,
            "phase": self.truncate(phase, 14),
            "condition": self.truncate(condition_str, 510),
            "intervention": self.truncate(intervention_str, 510),
            "primary_outcome": self.truncate(primary_outcome_str, 1022),
            "enrollment": enrollment,
            "status": self.truncate(status, 30),
            "start_date": self.truncate(start_date, 30),
            "completion_date": self.truncate(completion_date, 30),
            "sponsor": self.truncate(sponsor, 254),
            "key_findings": "",
        }

        return IngestRecord(
            text=self.truncate(text, 4096),
            metadata=metadata,
            collection=self.collection,
            source="ClinicalTrials.gov",
            source_id=nct_id,
        )

    @staticmethod
    def _extract_trial_name(title: str) -> str:
        """Extract trial acronym from title if present.

        Looks for patterns like ``PARADIGM-HF:``, ``(DAPA-HF)``, or
        ``The ISCHEMIA Trial``.

        Args:
            title: Full official trial title.

        Returns:
            Extracted acronym, or empty string if none found.
        """
        # Match acronym in parentheses: (PARADIGM-HF)
        paren_match = re.search(r"\(([A-Z][A-Z0-9\-]{2,20})\)", title)
        if paren_match:
            return paren_match.group(1)

        # Match acronym before colon: DAPA-HF: ...
        colon_match = re.match(r"^([A-Z][A-Z0-9\-]{2,20}):", title)
        if colon_match:
            return colon_match.group(1)

        # Match "The ACRONYM Trial"
        the_match = re.search(
            r"\bThe\s+([A-Z][A-Z0-9\-]{2,20})\s+(?:Trial|Study)\b", title
        )
        if the_match:
            return the_match.group(1)

        return ""
