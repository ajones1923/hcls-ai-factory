"""ACC/AHA 2022 Heart Failure Guidelines GDMT Optimization Engine.

Implements guideline-directed medical therapy (GDMT) optimization for heart
failure patients based on the 2022 AHA/ACC/HFSA Guideline for the Management
of Heart Failure (Heidenreich et al., Circulation 2022;145:e895-e1032).

Covers:
  - Four-pillar HFrEF GDMT (beta-blocker, ARNI/ACEi/ARB, MRA, SGLT2i)
  - Additional therapies (hydralazine-ISDN, ivabradine, digoxin, vericiguat)
  - Device therapy eligibility (ICD, CRT, CRT-D)
  - HFpEF-specific therapies (SGLT2i, GLP-1 RA, diuretics)
  - Titration planning, contraindication checking, drug interaction screening
  - Monitoring plan generation

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import re
from copy import deepcopy
from typing import Any, Dict, List, Optional, Tuple

from src.models import (
    EjectionFractionCategory,
    EvidenceLevel,
    GDMTMedication,
    GDMTPillar,
    GDMTRecommendation,
    GDMTStatus,
    GuidelineClass,
    HeartFailureClass,
)

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# GDMT DRUG DATABASE
# ═══════════════════════════════════════════════════════════════════════════════

GDMT_DRUG_DATABASE: Dict[GDMTPillar, Dict[str, Dict[str, Any]]] = {
    # ── Pillar 1: Beta-Blockers ──────────────────────────────────────────
    GDMTPillar.BETA_BLOCKER: {
        "carvedilol": {
            "starting_dose": "3.125 mg BID",
            "target_dose": "25 mg BID (50 mg BID if >85 kg)",
            "titration_interval": "2 weeks",
            "titration_steps": [
                "3.125 mg BID",
                "6.25 mg BID",
                "12.5 mg BID",
                "25 mg BID",
                "50 mg BID (if >85 kg)",
            ],
            "contraindications": [
                "severe bradycardia (HR <50 bpm)",
                "2nd or 3rd degree AV block without pacemaker",
                "cardiogenic shock",
                "decompensated heart failure requiring inotropes",
                "severe reactive airway disease",
                "sick sinus syndrome without pacemaker",
            ],
            "cautions": [
                "diabetes (may mask hypoglycemia symptoms)",
                "peripheral vascular disease",
                "hepatic impairment",
            ],
            "monitoring": [
                "heart rate (target resting HR 60-70 bpm)",
                "blood pressure (hold if SBP <90 mmHg)",
                "symptoms of fatigue or dizziness",
                "weight (fluid retention)",
                "blood glucose in diabetics",
            ],
            "key_trial": "COPERNICUS",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "35% reduction in all-cause mortality (COPERNICUS)",
            "mechanism": "Non-selective beta-blocker with alpha-1 blocking activity",
            "preferred_in": ["diabetes", "peripheral vascular disease"],
        },
        "metoprolol_succinate": {
            "starting_dose": "12.5-25 mg daily",
            "target_dose": "200 mg daily",
            "titration_interval": "2 weeks",
            "titration_steps": [
                "12.5 mg daily",
                "25 mg daily",
                "50 mg daily",
                "100 mg daily",
                "200 mg daily",
            ],
            "contraindications": [
                "severe bradycardia (HR <50 bpm)",
                "2nd or 3rd degree AV block without pacemaker",
                "cardiogenic shock",
                "decompensated heart failure requiring inotropes",
                "sick sinus syndrome without pacemaker",
            ],
            "cautions": [
                "reactive airway disease (beta-1 selective but caution at high doses)",
                "peripheral vascular disease",
                "abrupt withdrawal may exacerbate angina",
            ],
            "monitoring": [
                "heart rate (target resting HR 60-70 bpm)",
                "blood pressure (hold if SBP <90 mmHg)",
                "symptoms of fatigue or dizziness",
                "weight (fluid retention)",
            ],
            "key_trial": "MERIT-HF",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "34% reduction in all-cause mortality (MERIT-HF)",
            "mechanism": "Selective beta-1 adrenergic receptor blocker (extended release)",
            "note": "Must use succinate (extended-release) formulation; tartrate NOT indicated for HF",
        },
        "bisoprolol": {
            "starting_dose": "1.25 mg daily",
            "target_dose": "10 mg daily",
            "titration_interval": "2 weeks",
            "titration_steps": [
                "1.25 mg daily",
                "2.5 mg daily",
                "3.75 mg daily",
                "5 mg daily",
                "7.5 mg daily",
                "10 mg daily",
            ],
            "contraindications": [
                "severe bradycardia (HR <50 bpm)",
                "2nd or 3rd degree AV block without pacemaker",
                "cardiogenic shock",
                "decompensated heart failure requiring inotropes",
                "sick sinus syndrome without pacemaker",
                "severe reactive airway disease",
            ],
            "cautions": [
                "renal impairment (dose adjustment may be needed)",
                "hepatic impairment",
            ],
            "monitoring": [
                "heart rate (target resting HR 60-70 bpm)",
                "blood pressure (hold if SBP <90 mmHg)",
                "symptoms of fatigue or dizziness",
            ],
            "key_trial": "CIBIS-II",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "34% reduction in all-cause mortality (CIBIS-II)",
            "mechanism": "Highly selective beta-1 adrenergic receptor blocker",
        },
    },
    # ── Pillar 2: ARNI / ACEi / ARB ─────────────────────────────────────
    GDMTPillar.ARNI_ACEI_ARB: {
        "sacubitril_valsartan": {
            "starting_dose": "24/26 mg BID (49/51 mg BID if tolerating ACEi/ARB)",
            "target_dose": "97/103 mg BID",
            "titration_interval": "2-4 weeks",
            "titration_steps": [
                "24/26 mg BID",
                "49/51 mg BID",
                "97/103 mg BID",
            ],
            "contraindications": [
                "history of angioedema with ACEi or ARB",
                "concurrent ACEi use (36-hour washout required)",
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to sacubitril or valsartan",
                "severe hepatic impairment (Child-Pugh C)",
                "concomitant aliskiren use in diabetic patients",
            ],
            "cautions": [
                "SBP <100 mmHg (start at lower dose)",
                "eGFR <30 mL/min/1.73m2 (start at lower dose, monitor closely)",
                "moderate hepatic impairment (Child-Pugh B)",
                "hyperkalemia (K+ >5.0 mEq/L)",
                "volume depletion",
            ],
            "monitoring": [
                "blood pressure (SBP >=100 mmHg preferred for titration)",
                "serum potassium (maintain <5.0 mEq/L)",
                "serum creatinine / eGFR",
                "symptoms of hypotension or dizziness",
                "angioedema symptoms",
            ],
            "key_trial": "PARADIGM-HF",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "20% reduction in CV death (PARADIGM-HF vs enalapril)",
            "mechanism": "Neprilysin inhibitor + ARB (angiotensin receptor-neprilysin inhibitor)",
            "preferred": True,
            "note": "Preferred over ACEi/ARB. 36-hour washout from ACEi before starting.",
        },
        "enalapril": {
            "starting_dose": "2.5 mg BID",
            "target_dose": "10-20 mg BID",
            "titration_interval": "1-2 weeks",
            "titration_steps": [
                "2.5 mg BID",
                "5 mg BID",
                "10 mg BID",
                "20 mg BID",
            ],
            "contraindications": [
                "history of angioedema",
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ACE inhibitors",
                "concomitant ARNI use",
            ],
            "cautions": [
                "renal impairment (eGFR <30 mL/min/1.73m2)",
                "hyperkalemia (K+ >5.0 mEq/L)",
                "hypotension (SBP <90 mmHg)",
                "volume depletion",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR (check 1-2 weeks after initiation/titration)",
                "cough (class effect, switch to ARB if intolerable)",
            ],
            "key_trial": "SOLVD-Treatment",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "16% reduction in all-cause mortality (SOLVD-Treatment)",
            "mechanism": "Angiotensin-converting enzyme inhibitor",
        },
        "lisinopril": {
            "starting_dose": "2.5-5 mg daily",
            "target_dose": "20-40 mg daily",
            "titration_interval": "1-2 weeks",
            "titration_steps": [
                "2.5 mg daily",
                "5 mg daily",
                "10 mg daily",
                "20 mg daily",
                "40 mg daily",
            ],
            "contraindications": [
                "history of angioedema",
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ACE inhibitors",
                "concomitant ARNI use",
            ],
            "cautions": [
                "renal impairment",
                "hyperkalemia",
                "hypotension",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR",
                "cough",
            ],
            "key_trial": "ATLAS",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mechanism": "Angiotensin-converting enzyme inhibitor",
        },
        "ramipril": {
            "starting_dose": "1.25-2.5 mg daily",
            "target_dose": "10 mg daily",
            "titration_interval": "1-2 weeks",
            "titration_steps": [
                "1.25 mg daily",
                "2.5 mg daily",
                "5 mg daily",
                "10 mg daily",
            ],
            "contraindications": [
                "history of angioedema",
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ACE inhibitors",
                "concomitant ARNI use",
            ],
            "cautions": [
                "renal impairment",
                "hyperkalemia",
                "hypotension",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR",
                "cough",
            ],
            "key_trial": "AIRE",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "27% reduction in all-cause mortality post-MI (AIRE)",
            "mechanism": "Angiotensin-converting enzyme inhibitor",
        },
        "losartan": {
            "starting_dose": "25-50 mg daily",
            "target_dose": "150 mg daily",
            "titration_interval": "1-2 weeks",
            "titration_steps": [
                "25 mg daily",
                "50 mg daily",
                "100 mg daily",
                "150 mg daily",
            ],
            "contraindications": [
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ARBs",
                "concomitant ARNI use",
                "concomitant ACEi use",
            ],
            "cautions": [
                "renal impairment",
                "hyperkalemia",
                "hypotension",
                "hepatic impairment",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR",
            ],
            "key_trial": "HEAAL",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mechanism": "Angiotensin II receptor blocker",
            "note": "Use if ACEi not tolerated (cough). Transition to ARNI when possible.",
        },
        "valsartan": {
            "starting_dose": "20-40 mg BID",
            "target_dose": "160 mg BID",
            "titration_interval": "2 weeks",
            "titration_steps": [
                "20 mg BID",
                "40 mg BID",
                "80 mg BID",
                "160 mg BID",
            ],
            "contraindications": [
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ARBs",
                "concomitant ARNI use",
                "concomitant ACEi use",
            ],
            "cautions": [
                "renal impairment",
                "hyperkalemia",
                "hypotension",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR",
            ],
            "key_trial": "Val-HeFT",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "13% reduction in combined mortality/morbidity (Val-HeFT)",
            "mechanism": "Angiotensin II receptor blocker",
            "note": "Use if ACEi not tolerated. Transition to ARNI when possible.",
        },
        "candesartan": {
            "starting_dose": "4-8 mg daily",
            "target_dose": "32 mg daily",
            "titration_interval": "2 weeks",
            "titration_steps": [
                "4 mg daily",
                "8 mg daily",
                "16 mg daily",
                "32 mg daily",
            ],
            "contraindications": [
                "bilateral renal artery stenosis",
                "pregnancy",
                "hypersensitivity to ARBs",
                "concomitant ARNI use",
                "concomitant ACEi use",
            ],
            "cautions": [
                "renal impairment",
                "hyperkalemia",
                "hypotension",
            ],
            "monitoring": [
                "blood pressure",
                "serum potassium",
                "serum creatinine / eGFR",
            ],
            "key_trial": "CHARM-Alternative",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_B_R,
            },
            "mortality_reduction": "23% reduction in CV death or HF hospitalization (CHARM-Alternative)",
            "mechanism": "Angiotensin II receptor blocker",
            "note": "Use if ACEi not tolerated. Transition to ARNI when possible.",
        },
    },
    # ── Pillar 3: Mineralocorticoid Receptor Antagonists (MRA) ───────────
    GDMTPillar.MRA: {
        "spironolactone": {
            "starting_dose": "12.5-25 mg daily",
            "target_dose": "25-50 mg daily",
            "titration_interval": "4-8 weeks",
            "titration_steps": [
                "12.5 mg daily",
                "25 mg daily",
                "50 mg daily",
            ],
            "contraindications": [
                "serum potassium >5.0 mEq/L",
                "eGFR <30 mL/min/1.73m2",
                "concomitant potassium-sparing diuretic",
                "Addison's disease",
                "hypersensitivity to spironolactone",
            ],
            "cautions": [
                "eGFR 30-45 mL/min/1.73m2 (use lower dose, monitor frequently)",
                "concomitant ACEi/ARB/ARNI (increased hyperkalemia risk)",
                "gynecomastia or breast tenderness (consider switch to eplerenone)",
                "hepatic impairment",
            ],
            "monitoring": [
                "serum potassium (check within 1 week of initiation, then at 1 month, "
                "then every 3 months)",
                "serum creatinine / eGFR",
                "signs of gynecomastia or breast pain",
                "blood pressure",
            ],
            "key_trial": "RALES",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "30% reduction in all-cause mortality (RALES)",
            "mechanism": "Non-selective mineralocorticoid receptor antagonist",
            "note": "Preferred first-line MRA due to cost. Switch to eplerenone if gynecomastia.",
        },
        "eplerenone": {
            "starting_dose": "25 mg daily",
            "target_dose": "50 mg daily",
            "titration_interval": "4 weeks",
            "titration_steps": [
                "25 mg daily",
                "50 mg daily",
            ],
            "contraindications": [
                "serum potassium >5.0 mEq/L",
                "eGFR <30 mL/min/1.73m2",
                "concomitant strong CYP3A4 inhibitors (ketoconazole, itraconazole, "
                "nefazodone, ritonavir, nelfinavir, clarithromycin)",
                "concomitant potassium-sparing diuretic",
                "hypersensitivity to eplerenone",
            ],
            "cautions": [
                "eGFR 30-49 mL/min/1.73m2 (monitor potassium closely)",
                "concomitant ACEi/ARB/ARNI (increased hyperkalemia risk)",
                "moderate CYP3A4 inhibitors (erythromycin, fluconazole, verapamil)",
                "hepatic impairment",
            ],
            "monitoring": [
                "serum potassium (check within 1 week, then 1 month, then every 3 months)",
                "serum creatinine / eGFR",
                "blood pressure",
            ],
            "key_trial": "EMPHASIS-HF",
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "24% reduction in CV death or HF hospitalization (EMPHASIS-HF)",
            "mechanism": "Selective mineralocorticoid receptor antagonist",
            "note": "Preferred if spironolactone causes gynecomastia. Also studied post-MI (EPHESUS).",
        },
    },
    # ── Pillar 4: SGLT2 Inhibitors ──────────────────────────────────────
    GDMTPillar.SGLT2I: {
        "dapagliflozin": {
            "starting_dose": "10 mg daily",
            "target_dose": "10 mg daily",
            "titration_interval": "N/A (no titration required)",
            "titration_steps": [
                "10 mg daily",
            ],
            "contraindications": [
                "type 1 diabetes (risk of DKA)",
                "hypersensitivity to dapagliflozin",
                "dialysis",
            ],
            "cautions": [
                "eGFR <20 mL/min/1.73m2 (limited data; benefit for HF outcomes "
                "persists at lower eGFR but glycemic efficacy reduced)",
                "recurrent urogenital infections",
                "volume depletion or hypotension",
                "history of diabetic ketoacidosis",
                "concomitant loop diuretics (may need dose reduction)",
            ],
            "monitoring": [
                "eGFR (expected initial dip of 10-15%, stabilizes within weeks)",
                "blood glucose in diabetics",
                "volume status and blood pressure",
                "signs of genital mycotic infections",
                "ketone monitoring if intercurrent illness",
            ],
            "key_trial": "DAPA-HF",
            "additional_trials": ["DELIVER (HFpEF)", "DAPA-CKD"],
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "17% reduction in CV death or worsening HF (DAPA-HF)",
            "mechanism": "Sodium-glucose cotransporter 2 inhibitor",
            "note": "Benefits independent of diabetes status. Works across EF spectrum.",
        },
        "empagliflozin": {
            "starting_dose": "10 mg daily",
            "target_dose": "10 mg daily",
            "titration_interval": "N/A (no titration required)",
            "titration_steps": [
                "10 mg daily",
            ],
            "contraindications": [
                "type 1 diabetes (risk of DKA)",
                "hypersensitivity to empagliflozin",
                "dialysis",
            ],
            "cautions": [
                "eGFR <20 mL/min/1.73m2 (limited data)",
                "recurrent urogenital infections",
                "volume depletion or hypotension",
                "history of diabetic ketoacidosis",
                "concomitant loop diuretics (may need dose reduction)",
            ],
            "monitoring": [
                "eGFR (expected initial dip, stabilizes)",
                "blood glucose in diabetics",
                "volume status and blood pressure",
                "signs of genital mycotic infections",
                "ketone monitoring if intercurrent illness",
            ],
            "key_trial": "EMPEROR-Reduced",
            "additional_trials": ["EMPEROR-Preserved (HFpEF)", "EMPA-REG OUTCOME"],
            "evidence": {
                "class": GuidelineClass.CLASS_I,
                "level": EvidenceLevel.LEVEL_A,
            },
            "mortality_reduction": "25% reduction in CV death or HF hospitalization (EMPEROR-Reduced)",
            "mechanism": "Sodium-glucose cotransporter 2 inhibitor",
            "note": "Benefits independent of diabetes status. Also proven in HFpEF.",
        },
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# ADDITIONAL THERAPIES (beyond the 4 pillars)
# ═══════════════════════════════════════════════════════════════════════════════

ADDITIONAL_THERAPIES: Dict[str, Dict[str, Any]] = {
    "hydralazine_isosorbide_dinitrate": {
        "combination": "Hydralazine 37.5 mg + Isosorbide dinitrate 20 mg TID",
        "starting_dose": "Hydralazine 25 mg + ISDN 20 mg TID",
        "target_dose": "Hydralazine 75 mg + ISDN 40 mg TID",
        "titration_interval": "2-4 weeks",
        "titration_steps": [
            "Hydralazine 25 mg + ISDN 20 mg TID",
            "Hydralazine 50 mg + ISDN 20 mg TID",
            "Hydralazine 50 mg + ISDN 40 mg TID",
            "Hydralazine 75 mg + ISDN 40 mg TID",
        ],
        "indication": (
            "Class I recommendation for self-identified African American patients "
            "with NYHA III-IV HFrEF on optimal GDMT. Class IIa for patients who "
            "cannot tolerate ACEi/ARB/ARNI."
        ),
        "contraindications": [
            "concomitant PDE-5 inhibitors (sildenafil, tadalafil)",
            "severe hypotension (SBP <90 mmHg)",
            "drug-induced lupus from hydralazine",
            "right ventricular infarction",
        ],
        "monitoring": [
            "blood pressure",
            "heart rate (reflex tachycardia)",
            "ANA titer if prolonged hydralazine use (drug-induced lupus)",
            "headache (common; often improves with continued use)",
        ],
        "key_trial": "A-HeFT",
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_A,
        },
        "mortality_reduction": "43% reduction in all-cause mortality in African Americans (A-HeFT)",
        "note": (
            "Fixed-dose combination (BiDil) available. Also consider if intolerant to "
            "ACEi/ARB/ARNI regardless of race (Class IIa, LOE B-R)."
        ),
    },
    "ivabradine": {
        "starting_dose": "2.5-5 mg BID",
        "target_dose": "7.5 mg BID",
        "titration_interval": "2 weeks",
        "titration_steps": [
            "2.5 mg BID",
            "5 mg BID",
            "7.5 mg BID",
        ],
        "indication": (
            "Resting HR >=70 bpm in sinus rhythm despite maximally tolerated "
            "beta-blocker dose. LVEF <=35%. NYHA II-III."
        ),
        "contraindications": [
            "acute decompensated heart failure",
            "SBP <90 mmHg",
            "sick sinus syndrome, SA block, or 3rd degree AV block without pacemaker",
            "resting HR <60 bpm before initiation",
            "severe hepatic impairment",
            "atrial fibrillation or flutter",
            "concomitant strong CYP3A4 inhibitors",
            "pregnancy",
        ],
        "monitoring": [
            "heart rate (reduce dose if HR <50 bpm or symptomatic bradycardia)",
            "visual disturbances (phosphenes; transient luminous phenomena)",
            "atrial fibrillation (discontinue if develops)",
            "blood pressure",
        ],
        "key_trial": "SHIFT",
        "evidence": {
            "class": GuidelineClass.CLASS_IIA,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "outcome": "18% reduction in CV death or HF hospitalization (SHIFT)",
        "mechanism": "Selective If channel inhibitor; reduces HR without negative inotropy",
    },
    "digoxin": {
        "starting_dose": "0.125 mg daily",
        "target_dose": "0.125-0.25 mg daily (target level 0.5-0.9 ng/mL)",
        "indication": (
            "Persistent symptoms despite optimal GDMT. May reduce HF hospitalizations "
            "but does NOT reduce mortality."
        ),
        "contraindications": [
            "hypokalemia (K+ <3.5 mEq/L increases toxicity risk)",
            "hypomagnesemia",
            "hypercalcemia",
            "severe sinus node or AV nodal disease without pacemaker",
            "ventricular tachycardia/fibrillation",
            "hypertrophic obstructive cardiomyopathy (HOCM)",
            "Wolff-Parkinson-White (WPW) syndrome",
            "constrictive pericarditis",
        ],
        "cautions": [
            "renal impairment (reduce dose; digoxin is renally cleared)",
            "elderly patients (reduce dose)",
            "hypothyroidism (increased sensitivity)",
            "concomitant amiodarone (reduce digoxin dose by 50%)",
            "concomitant verapamil/diltiazem (increase digoxin levels)",
        ],
        "monitoring": [
            "serum digoxin level (target 0.5-0.9 ng/mL; check 6+ hours post-dose)",
            "serum potassium (hypokalemia potentiates toxicity)",
            "serum magnesium",
            "serum creatinine / eGFR",
            "ECG (signs of toxicity: ST scooping, PAT with block, bidirectional VT)",
            "GI symptoms (nausea, anorexia - early toxicity signs)",
        ],
        "key_trial": "DIG",
        "evidence": {
            "class": GuidelineClass.CLASS_IIB,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "outcome": "No mortality benefit; 28% reduction in HF hospitalization (DIG trial)",
        "mechanism": "Cardiac glycoside; positive inotrope via Na/K-ATPase inhibition",
        "note": (
            "Keep levels 0.5-0.9 ng/mL. Higher levels (>1.0) associated with increased "
            "mortality. Post-hoc analyses suggest harm at levels >1.2 ng/mL."
        ),
    },
    "vericiguat": {
        "starting_dose": "2.5 mg daily",
        "target_dose": "10 mg daily",
        "titration_interval": "2 weeks",
        "titration_steps": [
            "2.5 mg daily",
            "5 mg daily",
            "10 mg daily",
        ],
        "indication": (
            "Worsening HFrEF (LVEF <45%) despite optimal GDMT. Recent HF hospitalization "
            "or need for outpatient IV diuretics within prior 6 months."
        ),
        "contraindications": [
            "concomitant nitrates or nitric oxide donors",
            "concomitant PDE-5 inhibitors",
            "concomitant soluble guanylate cyclase stimulators (riociguat)",
            "pregnancy",
            "severe hypotension (SBP <100 mmHg)",
        ],
        "cautions": [
            "symptomatic hypotension",
            "anemia (may worsen)",
        ],
        "monitoring": [
            "blood pressure (SBP; hold if <100 mmHg)",
            "hemoglobin (anemia monitoring)",
            "signs and symptoms of hypotension",
        ],
        "key_trial": "VICTORIA",
        "evidence": {
            "class": GuidelineClass.CLASS_IIB,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "outcome": "10% reduction in CV death or first HF hospitalization (VICTORIA)",
        "mechanism": "Soluble guanylate cyclase (sGC) stimulator; enhances NO-sGC-cGMP pathway",
    },
    "finerenone": {
        "starting_dose": "10 mg daily (eGFR >=60) or 10 mg daily (eGFR 25-59)",
        "target_dose": "20 mg daily",
        "indication": (
            "CKD with type 2 diabetes and albuminuria, emerging HFrEF data"
        ),
        "contraindications": [
            "eGFR <25",
            "potassium >5.0",
            "severe hepatic impairment",
        ],
        "monitoring": [
            "potassium (1 week, then monthly)",
            "eGFR (baseline, 1 month, then quarterly)",
            "serum creatinine",
        ],
        "key_trials": ["FIDELIO-DKD", "FIGARO-DKD", "FINEARTS-HF"],
        "evidence": {
            "class": GuidelineClass.CLASS_IIA,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "mechanism": (
            "Non-steroidal selective mineralocorticoid receptor antagonist "
            "— less hyperkalemia than spironolactone/eplerenone"
        ),
    },
    "omecamtiv_mecarbil": {
        "starting_dose": "25 mg BID",
        "target_dose": "50 mg BID",
        "indication": (
            "HFrEF LVEF <=35% despite optimized GDMT, NYHA II-IV"
        ),
        "contraindications": [
            "severe hepatic impairment",
            "concomitant strong CYP3A4 inhibitors",
        ],
        "monitoring": [
            "LVEF (serial echo)",
            "troponin (baseline, 2 weeks)",
            "BP (modest reduction expected)",
        ],
        "key_trial": "GALACTIC-HF",
        "evidence": {
            "class": GuidelineClass.CLASS_IIB,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "outcome": "8% reduction in CV death/HF events (GALACTIC-HF)",
        "mechanism": (
            "Selective cardiac myosin activator; increases systolic ejection time "
            "and stroke volume without increasing intracellular calcium"
        ),
    },
    "sotagliflozin": {
        "starting_dose": "200 mg daily",
        "target_dose": "400 mg daily",
        "indication": (
            "HF (HFrEF and HFpEF) with type 2 diabetes; initiated before or "
            "shortly after discharge for acute decompensated HF"
        ),
        "contraindications": [
            "type 1 diabetes",
            "eGFR <25",
            "dialysis",
        ],
        "monitoring": [
            "blood glucose",
            "ketone awareness (DKA risk)",
            "volume status",
            "eGFR",
        ],
        "key_trials": ["SOLOIST-WHF", "SCORED"],
        "evidence": {
            "class": GuidelineClass.CLASS_IIA,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "mechanism": (
            "Dual SGLT1/SGLT2 inhibitor — additional SGLT1 inhibition reduces "
            "postprandial glucose absorption and delays intestinal glucose uptake"
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# DEVICE THERAPY CRITERIA
# ═══════════════════════════════════════════════════════════════════════════════

DEVICE_THERAPY_CRITERIA: Dict[str, Dict[str, Any]] = {
    "ICD": {
        "description": "Implantable Cardioverter-Defibrillator for primary prevention of SCD",
        "criteria": {
            "lvef": "<=35%",
            "nyha_class": ["II", "III"],
            "gdmt_duration": ">=90 days on optimal GDMT (or >=40 days post-MI)",
            "post_mi_waiting": ">=40 days post-MI",
            "expected_survival": ">1 year with meaningful functional status",
            "exclusion": [
                "NYHA IV refractory to medical therapy (unless CRT or LVAD candidate)",
                "within 40 days of acute MI",
                "within 90 days of revascularization (CABG or PCI)",
                "incessant VT/VF (ablation first)",
            ],
        },
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_A,
        },
        "key_trials": ["SCD-HeFT", "MADIT-II"],
        "mortality_reduction": "23% reduction in all-cause mortality (SCD-HeFT)",
        "note": (
            "Reassess LVEF after 3-6 months of optimized GDMT before implantation. "
            "Some patients improve LVEF above threshold after medical optimization."
        ),
    },
    "CRT": {
        "description": "Cardiac Resynchronization Therapy (biventricular pacing)",
        "criteria": {
            "lvef": "<=35%",
            "qrs_duration": ">=150 ms with LBBB morphology (Class I); "
                            "120-149 ms with LBBB (Class IIa); "
                            ">=150 ms non-LBBB (Class IIa); "
                            "120-149 ms non-LBBB (Class IIb)",
            "rhythm": "sinus rhythm preferred (AF with rate control may benefit)",
            "nyha_class": ["II", "III", "IV (ambulatory)"],
            "on_gdmt": True,
        },
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_A,
        },
        "key_trials": ["COMPANION", "CARE-HF", "MADIT-CRT", "RAFT"],
        "outcome": (
            "Reduces mortality, HF hospitalization, and improves quality of life. "
            "CARE-HF: 36% reduction in all-cause mortality."
        ),
        "note": (
            "Greatest benefit with LBBB QRS >=150 ms. Response rate ~65-70%. "
            "Consider CRT even in NYHA I if LVEF <=30% and LBBB QRS >=150 ms."
        ),
    },
    "CRT_D": {
        "description": "CRT with Defibrillator (CRT-D) - combined resynchronization and ICD",
        "criteria": {
            "requirement": "Meets BOTH ICD and CRT criteria",
            "lvef": "<=35%",
            "qrs_duration": ">=150 ms with LBBB (Class I)",
            "nyha_class": ["II", "III"],
            "expected_survival": ">1 year",
        },
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_A,
        },
        "key_trials": ["COMPANION", "RAFT"],
        "note": (
            "Preferred over CRT-P (pacemaker only) when patient also meets ICD criteria. "
            "COMPANION: 36% reduction in all-cause mortality with CRT-D vs. optimal medical therapy."
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# HFpEF-SPECIFIC THERAPIES
# ═══════════════════════════════════════════════════════════════════════════════

HFPEF_THERAPIES: Dict[str, Dict[str, Any]] = {
    "sglt2i": {
        "drugs": ["dapagliflozin 10 mg daily", "empagliflozin 10 mg daily"],
        "indication": "HFpEF (LVEF >40%) - reduces HF hospitalization and CV death",
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_A,
        },
        "key_trials": ["EMPEROR-Preserved", "DELIVER"],
        "outcome": (
            "EMPEROR-Preserved: 21% reduction in CV death or HF hospitalization. "
            "DELIVER: 18% reduction in worsening HF or CV death."
        ),
        "note": "First Class I therapy for HFpEF with robust evidence.",
    },
    "diuretics": {
        "drugs": [
            "furosemide 20-600 mg daily (divided doses)",
            "bumetanide 0.5-10 mg daily",
            "torsemide 10-200 mg daily",
        ],
        "indication": "Volume overload / congestion in HFpEF",
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_C_EO,
        },
        "note": (
            "Use lowest effective dose. Torsemide may have advantages over furosemide "
            "(longer half-life, anti-aldosterone effect). Monitor electrolytes and renal function."
        ),
    },
    "glp1_ra": {
        "drugs": ["semaglutide 2.4 mg weekly"],
        "indication": "HFpEF with obesity (BMI >=30 kg/m2)",
        "evidence": {
            "class": GuidelineClass.CLASS_IIA,
            "level": EvidenceLevel.LEVEL_B_R,
        },
        "key_trials": ["STEP-HFpEF", "STEP-HFpEF DM"],
        "outcome": (
            "STEP-HFpEF: Significant improvement in symptoms (KCCQ), 6MWD, "
            "body weight, and CRP. Mean weight loss ~13%."
        ),
        "note": (
            "Emerging therapy for obese HFpEF phenotype. Addresses the obesity-driven "
            "inflammatory and metabolic contributors to HFpEF."
        ),
    },
    "treat_comorbidities": {
        "conditions": {
            "hypertension": {
                "target": "BP <130/80 mmHg",
                "agents": "ACEi/ARB, thiazide diuretics, amlodipine; avoid non-DHP CCBs",
            },
            "atrial_fibrillation": {
                "approach": "Rate vs rhythm control; anticoagulation per CHA2DS2-VASc",
                "note": "AF very common in HFpEF; contributes to symptoms",
            },
            "coronary_artery_disease": {
                "approach": "Standard CAD management; revascularize if ischemia-driven symptoms",
            },
            "obesity": {
                "approach": "Weight loss (GLP-1 RA, bariatric surgery if BMI >=40 or >=35 with comorbidities)",
            },
            "diabetes": {
                "approach": "SGLT2i preferred; metformin safe; avoid TZDs (fluid retention)",
            },
            "iron_deficiency": {
                "criteria": "Ferritin <100 or ferritin 100-300 with TSAT <20%",
                "treatment": "IV iron (ferric carboxymaltose or ferric derisomaltose)",
                "evidence": {
                    "class": GuidelineClass.CLASS_IIA,
                    "level": EvidenceLevel.LEVEL_B_R,
                },
                "key_trials": ["FAIR-HF", "AFFIRM-AHF", "IRONMAN"],
            },
            "sleep_disordered_breathing": {
                "approach": "CPAP for obstructive sleep apnea; avoid ASV in HFrEF with CSA",
                "note": "SERVE-HF showed harm with ASV in HFrEF with central sleep apnea",
            },
        },
        "evidence": {
            "class": GuidelineClass.CLASS_I,
            "level": EvidenceLevel.LEVEL_C_EO,
        },
        "note": (
            "HFpEF management is largely comorbidity-driven. Addressing each comorbidity "
            "individually often provides greater symptom improvement than any single therapy."
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# DRUG INTERACTION DATABASE (clinically significant interactions in HF)
# ═══════════════════════════════════════════════════════════════════════════════

_DRUG_INTERACTIONS: Dict[str, List[Dict[str, str]]] = {
    "sacubitril_valsartan": [
        {
            "interacting_drug": "enalapril",
            "severity": "contraindicated",
            "effect": "Increased risk of angioedema. 36-hour washout required.",
        },
        {
            "interacting_drug": "lisinopril",
            "severity": "contraindicated",
            "effect": "Increased risk of angioedema. 36-hour washout required.",
        },
        {
            "interacting_drug": "ramipril",
            "severity": "contraindicated",
            "effect": "Increased risk of angioedema. 36-hour washout required.",
        },
        {
            "interacting_drug": "aliskiren",
            "severity": "contraindicated",
            "effect": "Dual RAAS blockade contraindicated in diabetic patients.",
        },
        {
            "interacting_drug": "spironolactone",
            "severity": "monitor",
            "effect": "Increased hyperkalemia risk. Monitor K+ closely.",
        },
        {
            "interacting_drug": "eplerenone",
            "severity": "monitor",
            "effect": "Increased hyperkalemia risk. Monitor K+ closely.",
        },
        {
            "interacting_drug": "potassium_supplements",
            "severity": "monitor",
            "effect": "Increased hyperkalemia risk.",
        },
    ],
    "spironolactone": [
        {
            "interacting_drug": "enalapril",
            "severity": "monitor",
            "effect": "Hyperkalemia risk; monitor K+ within 3-7 days of initiation.",
        },
        {
            "interacting_drug": "sacubitril_valsartan",
            "severity": "monitor",
            "effect": "Hyperkalemia risk; monitor K+ within 3-7 days.",
        },
        {
            "interacting_drug": "trimethoprim",
            "severity": "major",
            "effect": "Marked hyperkalemia risk. Avoid combination or monitor K+ very closely.",
        },
        {
            "interacting_drug": "digoxin",
            "severity": "monitor",
            "effect": "Spironolactone may increase digoxin levels. Monitor digoxin level.",
        },
        {
            "interacting_drug": "nsaids",
            "severity": "major",
            "effect": "NSAIDs reduce diuretic efficacy and increase hyperkalemia/AKI risk.",
        },
        {
            "interacting_drug": "potassium_supplements",
            "severity": "major",
            "effect": "Life-threatening hyperkalemia risk. Avoid concurrent K+ supplements or monitor K+ within 48 hours.",
        },
    ],
    "digoxin": [
        {
            "interacting_drug": "amiodarone",
            "severity": "major",
            "effect": "Amiodarone increases digoxin levels ~70%. Reduce digoxin dose by 50%.",
        },
        {
            "interacting_drug": "verapamil",
            "severity": "major",
            "effect": "Verapamil increases digoxin levels. Also additive AV nodal block.",
        },
        {
            "interacting_drug": "diltiazem",
            "severity": "major",
            "effect": "Diltiazem increases digoxin levels. Additive AV nodal block.",
        },
        {
            "interacting_drug": "carvedilol",
            "severity": "monitor",
            "effect": "Carvedilol may increase digoxin levels by ~15%. Monitor level.",
        },
        {
            "interacting_drug": "spironolactone",
            "severity": "monitor",
            "effect": "May interfere with digoxin assays. Monitor clinical response.",
        },
    ],
    "ivabradine": [
        {
            "interacting_drug": "verapamil",
            "severity": "contraindicated",
            "effect": "Excessive bradycardia and CYP3A4 inhibition.",
        },
        {
            "interacting_drug": "diltiazem",
            "severity": "contraindicated",
            "effect": "Excessive bradycardia and CYP3A4 inhibition.",
        },
        {
            "interacting_drug": "ketoconazole",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor markedly increases ivabradine levels.",
        },
        {
            "interacting_drug": "carvedilol",
            "severity": "monitor",
            "effect": "Additive bradycardia. Monitor HR closely.",
        },
        {
            "interacting_drug": "metoprolol_succinate",
            "severity": "monitor",
            "effect": "Additive bradycardia. Expected when used together; monitor HR.",
        },
        {
            "interacting_drug": "bisoprolol",
            "severity": "monitor",
            "effect": "Additive bradycardia. Monitor HR.",
        },
    ],
    "vericiguat": [
        {
            "interacting_drug": "nitrates",
            "severity": "contraindicated",
            "effect": "Risk of severe hypotension. Do not combine.",
        },
        {
            "interacting_drug": "sildenafil",
            "severity": "contraindicated",
            "effect": "PDE-5 inhibitor increases cGMP; severe hypotension risk.",
        },
        {
            "interacting_drug": "tadalafil",
            "severity": "contraindicated",
            "effect": "PDE-5 inhibitor increases cGMP; severe hypotension risk.",
        },
        {
            "interacting_drug": "riociguat",
            "severity": "contraindicated",
            "effect": "Dual sGC stimulation; excessive hypotension.",
        },
    ],
    "hydralazine_isosorbide_dinitrate": [
        {
            "interacting_drug": "sildenafil",
            "severity": "contraindicated",
            "effect": "Nitrate + PDE-5 inhibitor: severe hypotension risk.",
        },
        {
            "interacting_drug": "tadalafil",
            "severity": "contraindicated",
            "effect": "Nitrate + PDE-5 inhibitor: severe hypotension risk.",
        },
        {
            "interacting_drug": "vericiguat",
            "severity": "contraindicated",
            "effect": "Nitrate + sGC stimulator: severe hypotension risk.",
        },
    ],
    "finerenone": [
        {
            "interacting_drug": "itraconazole",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases finerenone exposure 3-fold.",
        },
        {
            "interacting_drug": "ketoconazole",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases finerenone exposure 3-fold.",
        },
        {
            "interacting_drug": "ritonavir",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases finerenone exposure 3-fold.",
        },
        {
            "interacting_drug": "clarithromycin",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases finerenone exposure 3-fold.",
        },
        {
            "interacting_drug": "potassium_supplements",
            "severity": "contraindicated",
            "effect": "Concomitant potassium supplements or potassium-sparing diuretics — hyperkalemia risk.",
        },
        {
            "interacting_drug": "erythromycin",
            "severity": "monitor",
            "effect": "Moderate CYP3A4 inhibitor — dose adjustment may be needed.",
        },
        {
            "interacting_drug": "diltiazem",
            "severity": "monitor",
            "effect": "Moderate CYP3A4 inhibitor — dose adjustment may be needed.",
        },
        {
            "interacting_drug": "verapamil",
            "severity": "monitor",
            "effect": "Moderate CYP3A4 inhibitor — dose adjustment may be needed.",
        },
        {
            "interacting_drug": "fluconazole",
            "severity": "monitor",
            "effect": "Moderate CYP3A4 inhibitor — dose adjustment may be needed.",
        },
        {
            "interacting_drug": "enalapril",
            "severity": "monitor",
            "effect": "ACEi/ARB/ARNI — additive hyperkalemia risk, check K+ within 1 week.",
        },
        {
            "interacting_drug": "sacubitril_valsartan",
            "severity": "monitor",
            "effect": "ACEi/ARB/ARNI — additive hyperkalemia risk, check K+ within 1 week.",
        },
        {
            "interacting_drug": "losartan",
            "severity": "monitor",
            "effect": "ACEi/ARB/ARNI — additive hyperkalemia risk, check K+ within 1 week.",
        },
        {
            "interacting_drug": "valsartan",
            "severity": "monitor",
            "effect": "ACEi/ARB/ARNI — additive hyperkalemia risk, check K+ within 1 week.",
        },
    ],
    "omecamtiv_mecarbil": [
        {
            "interacting_drug": "itraconazole",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases omecamtiv mecarbil exposure.",
        },
        {
            "interacting_drug": "ketoconazole",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases omecamtiv mecarbil exposure.",
        },
        {
            "interacting_drug": "ritonavir",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases omecamtiv mecarbil exposure.",
        },
        {
            "interacting_drug": "clarithromycin",
            "severity": "contraindicated",
            "effect": "Strong CYP3A4 inhibitor increases omecamtiv mecarbil exposure.",
        },
        {
            "interacting_drug": "digoxin",
            "severity": "monitor",
            "effect": "May enhance inotropic effects. Monitor clinical response.",
        },
        {
            "interacting_drug": "warfarin",
            "severity": "monitor",
            "effect": "No significant interaction but monitor INR during initiation.",
        },
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# KNOWN DRUG NAME ALIASES
# ═══════════════════════════════════════════════════════════════════════════════

_DRUG_ALIASES: Dict[str, Tuple[GDMTPillar, str]] = {
    # Beta-blockers
    "coreg": (GDMTPillar.BETA_BLOCKER, "carvedilol"),
    "carvedilol": (GDMTPillar.BETA_BLOCKER, "carvedilol"),
    "toprol xl": (GDMTPillar.BETA_BLOCKER, "metoprolol_succinate"),
    "toprol-xl": (GDMTPillar.BETA_BLOCKER, "metoprolol_succinate"),
    "metoprolol succinate": (GDMTPillar.BETA_BLOCKER, "metoprolol_succinate"),
    "metoprolol er": (GDMTPillar.BETA_BLOCKER, "metoprolol_succinate"),
    "metoprolol xl": (GDMTPillar.BETA_BLOCKER, "metoprolol_succinate"),
    "bisoprolol": (GDMTPillar.BETA_BLOCKER, "bisoprolol"),
    "zebeta": (GDMTPillar.BETA_BLOCKER, "bisoprolol"),
    # ARNI / ACEi / ARB
    "entresto": (GDMTPillar.ARNI_ACEI_ARB, "sacubitril_valsartan"),
    "sacubitril/valsartan": (GDMTPillar.ARNI_ACEI_ARB, "sacubitril_valsartan"),
    "sacubitril_valsartan": (GDMTPillar.ARNI_ACEI_ARB, "sacubitril_valsartan"),
    "sacubitril valsartan": (GDMTPillar.ARNI_ACEI_ARB, "sacubitril_valsartan"),
    "enalapril": (GDMTPillar.ARNI_ACEI_ARB, "enalapril"),
    "vasotec": (GDMTPillar.ARNI_ACEI_ARB, "enalapril"),
    "lisinopril": (GDMTPillar.ARNI_ACEI_ARB, "lisinopril"),
    "zestril": (GDMTPillar.ARNI_ACEI_ARB, "lisinopril"),
    "prinivil": (GDMTPillar.ARNI_ACEI_ARB, "lisinopril"),
    "ramipril": (GDMTPillar.ARNI_ACEI_ARB, "ramipril"),
    "altace": (GDMTPillar.ARNI_ACEI_ARB, "ramipril"),
    "losartan": (GDMTPillar.ARNI_ACEI_ARB, "losartan"),
    "cozaar": (GDMTPillar.ARNI_ACEI_ARB, "losartan"),
    "valsartan": (GDMTPillar.ARNI_ACEI_ARB, "valsartan"),
    "diovan": (GDMTPillar.ARNI_ACEI_ARB, "valsartan"),
    "candesartan": (GDMTPillar.ARNI_ACEI_ARB, "candesartan"),
    "atacand": (GDMTPillar.ARNI_ACEI_ARB, "candesartan"),
    # MRA
    "spironolactone": (GDMTPillar.MRA, "spironolactone"),
    "aldactone": (GDMTPillar.MRA, "spironolactone"),
    "eplerenone": (GDMTPillar.MRA, "eplerenone"),
    "inspra": (GDMTPillar.MRA, "eplerenone"),
    # SGLT2i
    "dapagliflozin": (GDMTPillar.SGLT2I, "dapagliflozin"),
    "farxiga": (GDMTPillar.SGLT2I, "dapagliflozin"),
    "forxiga": (GDMTPillar.SGLT2I, "dapagliflozin"),
    "empagliflozin": (GDMTPillar.SGLT2I, "empagliflozin"),
    "jardiance": (GDMTPillar.SGLT2I, "empagliflozin"),
}

# Aliases for additional therapies (not tied to a specific GDMT pillar)
_ADDITIONAL_THERAPY_ALIASES: Dict[str, str] = {
    "kerendia": "finerenone",
    "omecamtiv": "omecamtiv_mecarbil",
    "xatmep": "sotagliflozin",
    "inpefa": "sotagliflozin",
}


# ═══════════════════════════════════════════════════════════════════════════════
# COMORBIDITY ADJUSTMENTS
# ═══════════════════════════════════════════════════════════════════════════════

COMORBIDITY_ADJUSTMENTS: Dict[str, Dict[str, Any]] = {
    "ckd_stage_4": {
        "description": "CKD Stage 4 (eGFR 15-29)",
        "adjustments": [
            "Prefer finerenone over spironolactone if K+ <5.0",
            "SGLT2i can be continued if started with eGFR >20",
        ],
    },
    "liver_disease": {
        "description": "Liver disease (Child-Pugh B)",
        "adjustments": [
            "Reduce carvedilol dose",
            "Avoid spironolactone if ascites",
            "Eplerenone preferred MRA",
        ],
    },
    "frailty_elderly": {
        "description": "Frailty/Elderly (>80 years)",
        "adjustments": [
            "Consider lower target doses",
            "Slower titration (4-week intervals)",
            "Prioritize tolerability over target dose",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: DOSE PARSING
# ═══════════════════════════════════════════════════════════════════════════════

_DOSE_PATTERN = re.compile(
    r"(\d+(?:\.\d+)?)\s*(?:/\s*(\d+(?:\.\d+)?))?\s*mg",
    re.IGNORECASE,
)


def _parse_dose_mg(dose_str: str) -> float:
    """Extract the primary numeric mg value from a dose string.

    For combination drugs like sacubitril/valsartan "97/103 mg BID",
    returns the first number (97).  Returns 0.0 if unparseable.
    """
    if not dose_str:
        return 0.0
    match = _DOSE_PATTERN.search(dose_str)
    if match:
        return float(match.group(1))
    return 0.0


def _normalise_drug_name(name: str) -> str:
    """Normalise a drug name to its canonical database key."""
    key = name.strip().lower().replace("-", " ").replace("/", " ")
    # Remove dosage info that might be appended
    key = re.sub(r"\d+(\.\d+)?\s*mg.*", "", key).strip()
    return key


# ═══════════════════════════════════════════════════════════════════════════════
# GDMT OPTIMIZER
# ═══════════════════════════════════════════════════════════════════════════════


class GDMTOptimizer:
    """GDMT optimization engine based on ACC/AHA 2022 HF Guidelines.

    Analyses a patient's current heart failure medications, lab values, and
    vitals against guideline targets and generates structured recommendations
    for initiation, titration, additional therapies, device therapy, and
    monitoring.

    Clinical principles encoded
    --------------------------
    * HFrEF: start all 4 pillars simultaneously (do NOT wait for one to
      reach target before starting the next).
    * Prefer ARNI over ACEi/ARB in eligible patients.
    * Check K+ before MRA (must be <5.0 mEq/L); use caution if eGFR <30.
    * HR must be adequate for beta-blocker titration (target resting 60-70).
    * SBP must support ARNI/ACEi/ARB titration (SBP >=100 preferred).
    * After optimization: reassess EF at 3-6 months.  If EF improves >40%,
      reclassify as HFimpEF and continue all GDMT.
    """

    GUIDELINE_REFERENCE = (
        "2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure. "
        "Heidenreich PA, et al. Circulation. 2022;145:e895-e1032."
    )

    # Thresholds
    _HFREF_EF_CUTOFF = 40.0
    _HFMREF_UPPER = 49.0
    _HFPEF_EF_CUTOFF = 50.0

    _MIN_SBP_FOR_RAAS = 100.0  # mmHg
    _MIN_SBP_ABSOLUTE = 90.0
    _MIN_HR_FOR_BB = 55.0  # bpm
    _TARGET_HR_LOW = 60.0
    _TARGET_HR_HIGH = 70.0
    _MAX_K_FOR_MRA = 5.0  # mEq/L
    _MIN_EGFR_FOR_MRA = 30.0  # mL/min/1.73m2
    _MIN_EGFR_CAUTION_ARNI = 20.0
    _ICD_EF_CUTOFF = 35.0
    _CRT_EF_CUTOFF = 35.0
    _CRT_QRS_CLASS_I = 150  # ms
    _CRT_QRS_CLASS_IIA = 120  # ms
    _IVABRADINE_HR_THRESHOLD = 70  # bpm
    _IVABRADINE_EF_CUTOFF = 35.0

    def __init__(self) -> None:
        self._drug_db = deepcopy(GDMT_DRUG_DATABASE)
        self._additional = deepcopy(ADDITIONAL_THERAPIES)
        self._device = deepcopy(DEVICE_THERAPY_CRITERIA)
        self._hfpef = deepcopy(HFPEF_THERAPIES)
        self._interactions = deepcopy(_DRUG_INTERACTIONS)

    # ──────────────────────────────────────────────────────────────────────
    # PUBLIC API
    # ──────────────────────────────────────────────────────────────────────

    def optimize(
        self,
        lvef: float,
        nyha_class: str,
        current_medications: List[Dict[str, Any]],
        patient_data: Dict[str, Any],
    ) -> GDMTRecommendation:
        """Generate comprehensive GDMT optimization recommendations.

        Parameters
        ----------
        lvef : float
            Left ventricular ejection fraction (%).
        nyha_class : str
            NYHA functional class as string ("I", "II", "III", "IV").
        current_medications : list of dict
            Each dict should have at least ``{"name": str}`` and optionally
            ``"dose"``, ``"frequency"``, ``"status"``.
        patient_data : dict
            Patient context including labs, vitals, demographics.
            Expected keys (all optional): ``sbp``, ``hr``, ``potassium``,
            ``egfr``, ``creatinine``, ``age``, ``sex``, ``race``,
            ``weight_kg``, ``diabetes``, ``qrs_duration``, ``qrs_morphology``,
            ``rhythm``, ``days_post_mi``, ``days_on_gdmt``,
            ``contraindications`` (list of str), ``previous_ef`` (float).

        Returns
        -------
        GDMTRecommendation
            Structured recommendation object.
        """
        logger.info(
            "GDMT optimization requested: LVEF=%.1f%%, NYHA %s, %d current meds",
            lvef, nyha_class, len(current_medications),
        )

        # Ensure lvef is available in patient_data for downstream methods
        if "lvef" not in patient_data and "ef" not in patient_data:
            patient_data["lvef"] = lvef

        # 1. Classify EF
        ef_category = self._classify_ef(lvef, patient_data)

        # 2. Map NYHA class
        nyha = self._parse_nyha(nyha_class)

        # 3. Assess current GDMT
        current_gdmt = self._assess_current_gdmt(current_medications)

        # 4. Gather patient contraindications
        patient_contraindications: List[str] = patient_data.get("contraindications", [])

        # 5. Build recommendation object
        rec = GDMTRecommendation(
            ef_category=ef_category,
            current_meds=[med for med in current_gdmt.values()],
            recommendations=[],
            next_steps=[],
            guideline_references=[self.GUIDELINE_REFERENCE],
        )

        # 6. Route by EF category
        if ef_category in (EjectionFractionCategory.HFrEF, EjectionFractionCategory.HFimpEF):
            self._optimize_hfref(
                rec, ef_category, nyha, current_gdmt, patient_data, patient_contraindications,
            )
        elif ef_category == EjectionFractionCategory.HFmrEF:
            self._optimize_hfmref(
                rec, nyha, current_gdmt, patient_data, patient_contraindications,
            )
        else:
            self._optimize_hfpef(
                rec, nyha, current_gdmt, patient_data, patient_contraindications,
            )

        # 7. Assess additional therapies (applicable across EF categories)
        additional = self._assess_additional_therapies(
            ef_category, nyha, patient_data, current_medications,
        )
        if additional:
            rec.recommendations.extend(additional)

        # 8. Assess device therapy
        device_recs = self._assess_device_therapy(
            lvef,
            nyha,
            patient_data.get("qrs_duration"),
            patient_data.get("qrs_morphology"),
            patient_data,
        )
        if device_recs:
            rec.recommendations.extend(device_recs)

        # 9. Drug interaction check
        all_med_names = [
            _normalise_drug_name(m.get("name", "")) for m in current_medications
        ]
        interaction_warnings = self._check_all_interactions(all_med_names, rec)
        if interaction_warnings:
            rec.next_steps.extend(interaction_warnings)

        # 10. Generate monitoring plan
        monitoring = self._generate_monitoring_plan(rec, patient_data)
        rec.next_steps.extend(monitoring)

        # 11. Reassessment
        rec.next_steps.append(
            "Reassess LVEF with echocardiography at 3-6 months after GDMT optimisation."
        )
        if ef_category == EjectionFractionCategory.HFrEF:
            rec.next_steps.append(
                "If LVEF improves to >40%, reclassify as HFimpEF and CONTINUE all GDMT "
                "(do not withdraw therapies that led to improvement)."
            )

        logger.info(
            "GDMT optimization complete: %d recommendations, %d next steps",
            len(rec.recommendations), len(rec.next_steps),
        )
        return rec

    # ──────────────────────────────────────────────────────────────────────
    # EF CLASSIFICATION
    # ──────────────────────────────────────────────────────────────────────

    def _classify_ef(
        self, lvef: float, patient_data: Dict[str, Any],
    ) -> EjectionFractionCategory:
        """Classify EF into HFrEF / HFmrEF / HFpEF / HFimpEF.

        If previous_ef is provided in patient_data and was <=40% but current
        LVEF is >40%, the patient is classified as HFimpEF.
        """
        previous_ef = patient_data.get("previous_ef")

        # Check for HFimpEF: previously <=40%, now >40%
        if previous_ef is not None and previous_ef <= self._HFREF_EF_CUTOFF and lvef > self._HFREF_EF_CUTOFF:
            logger.info(
                "EF improved from %.1f%% to %.1f%% -> HFimpEF", previous_ef, lvef,
            )
            return EjectionFractionCategory.HFimpEF

        if lvef <= self._HFREF_EF_CUTOFF:
            return EjectionFractionCategory.HFrEF
        elif lvef <= self._HFMREF_UPPER:
            return EjectionFractionCategory.HFmrEF
        else:
            return EjectionFractionCategory.HFpEF

    # ──────────────────────────────────────────────────────────────────────
    # NYHA PARSING
    # ──────────────────────────────────────────────────────────────────────

    @staticmethod
    def _parse_nyha(nyha_str: str) -> HeartFailureClass:
        """Parse a NYHA class string into the enum."""
        mapping = {
            "I": HeartFailureClass.NYHA_I,
            "1": HeartFailureClass.NYHA_I,
            "II": HeartFailureClass.NYHA_II,
            "2": HeartFailureClass.NYHA_II,
            "III": HeartFailureClass.NYHA_III,
            "3": HeartFailureClass.NYHA_III,
            "IV": HeartFailureClass.NYHA_IV,
            "4": HeartFailureClass.NYHA_IV,
        }
        normalised = nyha_str.strip().upper().replace("NYHA", "").replace("CLASS", "").strip()
        result = mapping.get(normalised)
        if result is None:
            logger.warning("Unrecognised NYHA class '%s', defaulting to NYHA II", nyha_str)
            return HeartFailureClass.NYHA_II
        return result

    # ──────────────────────────────────────────────────────────────────────
    # CURRENT GDMT ASSESSMENT
    # ──────────────────────────────────────────────────────────────────────

    def _assess_current_gdmt(
        self, current_medications: List[Dict[str, Any]],
    ) -> Dict[GDMTPillar, GDMTMedication]:
        """Map current medications to GDMT pillars and assess status.

        Returns a dict keyed by pillar with a GDMTMedication for each pillar
        that the patient is currently on.  Pillars not represented in the
        patient's medication list are omitted (indicating a gap).
        """
        result: Dict[GDMTPillar, GDMTMedication] = {}

        for med in current_medications:
            raw_name = med.get("name", "")
            normalised = _normalise_drug_name(raw_name)
            alias_match = _DRUG_ALIASES.get(normalised)

            if alias_match is None:
                # Try partial matching
                for alias_key, alias_val in _DRUG_ALIASES.items():
                    if alias_key in normalised or normalised in alias_key:
                        alias_match = alias_val
                        break

            if alias_match is None:
                logger.debug("Medication '%s' not in GDMT database; skipping", raw_name)
                continue

            pillar, canonical_name = alias_match
            db_entry = self._drug_db.get(pillar, {}).get(canonical_name, {})
            target_dose = db_entry.get("target_dose", "unknown")
            current_dose = med.get("dose", "")

            # Determine status
            status = self._determine_med_status(
                current_dose, target_dose, med.get("status", ""),
            )

            gdmt_med = GDMTMedication(
                pillar=pillar,
                drug_name=canonical_name,
                current_dose=current_dose or None,
                target_dose=target_dose,
                status=status,
                contraindications=med.get("contraindications", []),
            )
            # Only keep one med per pillar (first match or prefer the one already present)
            if pillar not in result:
                result[pillar] = gdmt_med
            else:
                # If both present on same pillar, keep the one with higher status
                existing = result[pillar]
                if self._status_rank(status) > self._status_rank(existing.status):
                    result[pillar] = gdmt_med

        return result

    @staticmethod
    def _status_rank(status: GDMTStatus) -> int:
        """Rank GDMT statuses for comparison (higher = more optimised)."""
        ranking = {
            GDMTStatus.NOT_STARTED: 0,
            GDMTStatus.CONTRAINDICATED: 1,
            GDMTStatus.INTOLERANT: 1,
            GDMTStatus.INITIATED: 2,
            GDMTStatus.UPTITRATING: 3,
            GDMTStatus.AT_TARGET: 5,
        }
        return ranking.get(status, 0)

    def _determine_med_status(
        self, current_dose: str, target_dose: str, explicit_status: str,
    ) -> GDMTStatus:
        """Determine medication status from dose comparison or explicit status."""
        # Honour explicit status if provided
        explicit_lower = explicit_status.strip().lower()
        status_map = {
            "contraindicated": GDMTStatus.CONTRAINDICATED,
            "intolerant": GDMTStatus.INTOLERANT,
            "not started": GDMTStatus.NOT_STARTED,
            "not_started": GDMTStatus.NOT_STARTED,
            "at target": GDMTStatus.AT_TARGET,
            "at_target": GDMTStatus.AT_TARGET,
            "uptitrating": GDMTStatus.UPTITRATING,
            "initiated": GDMTStatus.INITIATED,
        }
        if explicit_lower in status_map:
            return status_map[explicit_lower]

        # Infer from dose comparison
        if not current_dose:
            return GDMTStatus.INITIATED

        current_mg = _parse_dose_mg(current_dose)
        target_mg = _parse_dose_mg(target_dose)

        if current_mg <= 0 or target_mg <= 0:
            return GDMTStatus.INITIATED

        ratio = current_mg / target_mg
        if ratio >= 0.95:
            return GDMTStatus.AT_TARGET
        elif ratio >= 0.25:
            return GDMTStatus.UPTITRATING
        else:
            return GDMTStatus.INITIATED

    # ──────────────────────────────────────────────────────────────────────
    # GAP IDENTIFICATION
    # ──────────────────────────────────────────────────────────────────────

    def _identify_gaps(
        self,
        ef_category: EjectionFractionCategory,
        current_gdmt: Dict[GDMTPillar, GDMTMedication],
        contraindications: List[str],
    ) -> List[str]:
        """Identify missing or sub-optimally treated GDMT pillars.

        For HFrEF, all 4 pillars should be on board.  For HFmrEF, SGLT2i
        is Class I; others are Class IIb.  For HFpEF, only SGLT2i is Class I.
        """
        gaps: List[str] = []

        if ef_category in (EjectionFractionCategory.HFrEF, EjectionFractionCategory.HFimpEF):
            required_pillars = [
                GDMTPillar.BETA_BLOCKER,
                GDMTPillar.ARNI_ACEI_ARB,
                GDMTPillar.MRA,
                GDMTPillar.SGLT2I,
            ]
        elif ef_category == EjectionFractionCategory.HFmrEF:
            required_pillars = [GDMTPillar.SGLT2I]
        else:
            required_pillars = [GDMTPillar.SGLT2I]

        for pillar in required_pillars:
            if pillar not in current_gdmt:
                gaps.append(
                    f"GAP: {pillar.value} not initiated - "
                    f"Class I indication per ACC/AHA 2022 guidelines"
                )
            else:
                med = current_gdmt[pillar]
                if med.status == GDMTStatus.CONTRAINDICATED:
                    gaps.append(
                        f"NOTE: {pillar.value} contraindicated ({med.drug_name}); "
                        f"consider alternative agent within pillar if available"
                    )
                elif med.status == GDMTStatus.INTOLERANT:
                    gaps.append(
                        f"NOTE: {pillar.value} intolerant ({med.drug_name}); "
                        f"consider alternative agent within pillar"
                    )
                elif med.status in (GDMTStatus.INITIATED, GDMTStatus.UPTITRATING):
                    gaps.append(
                        f"UPTITRATE: {pillar.value} ({med.drug_name}) currently "
                        f"{med.current_dose or 'sub-target'} -> target {med.target_dose}"
                    )

        # For HFrEF, also flag if using ACEi/ARB but not ARNI
        if ef_category in (EjectionFractionCategory.HFrEF, EjectionFractionCategory.HFimpEF):
            arni_pillar = current_gdmt.get(GDMTPillar.ARNI_ACEI_ARB)
            if arni_pillar and arni_pillar.drug_name != "sacubitril_valsartan":
                is_arni_contraindicated = any(
                    "angioedema" in c.lower() for c in contraindications
                )
                if not is_arni_contraindicated:
                    gaps.append(
                        f"SWITCH: Currently on {arni_pillar.drug_name}; guideline preference "
                        f"is to switch to sacubitril/valsartan (ARNI) per PARADIGM-HF "
                        f"(Class I, LOE A). 36-hour ACEi washout required before starting ARNI."
                    )

        return gaps

    # ──────────────────────────────────────────────────────────────────────
    # TITRATION PLANNING
    # ──────────────────────────────────────────────────────────────────────

    def _generate_titration_plan(
        self,
        medication: GDMTMedication,
        patient_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Generate a step-by-step titration plan for a GDMT medication.

        Takes into account current dose, target dose, titration interval,
        and patient haemodynamic parameters.
        """
        db_entry = self._drug_db.get(medication.pillar, {}).get(medication.drug_name, {})
        if not db_entry:
            return {
                "drug": medication.drug_name,
                "plan": "Unable to generate titration plan: drug not in database.",
            }

        titration_steps = db_entry.get("titration_steps", [])
        interval = db_entry.get("titration_interval", "2 weeks")
        monitoring = db_entry.get("monitoring", [])

        # Find current step
        current_dose_mg = _parse_dose_mg(medication.current_dose or "")
        current_step_idx = 0
        for i, step in enumerate(titration_steps):
            step_mg = _parse_dose_mg(step)
            if step_mg > 0 and current_dose_mg >= step_mg * 0.9:
                current_step_idx = i

        remaining_steps = titration_steps[current_step_idx + 1:]

        # Add haemodynamic precautions
        precautions: List[str] = []
        sbp = patient_data.get("sbp")
        hr = patient_data.get("hr")

        if medication.pillar == GDMTPillar.ARNI_ACEI_ARB:
            if sbp is not None and sbp < self._MIN_SBP_FOR_RAAS:
                precautions.append(
                    f"SBP is {sbp} mmHg (below 100). Titrate cautiously; consider "
                    f"reducing diuretic dose if euvolemic. Hold if SBP <90."
                )
        elif medication.pillar == GDMTPillar.BETA_BLOCKER:
            if hr is not None and hr < self._TARGET_HR_LOW:
                precautions.append(
                    f"Resting HR is {hr} bpm (below 60). Do not uptitrate. "
                    f"Reduce dose if HR <55 or symptomatic."
                )
            if sbp is not None and sbp < self._MIN_SBP_ABSOLUTE:
                precautions.append(
                    f"SBP is {sbp} mmHg (below 90). Hold beta-blocker titration."
                )

        return {
            "drug": medication.drug_name,
            "pillar": medication.pillar.value,
            "current_dose": medication.current_dose,
            "target_dose": db_entry.get("target_dose", "unknown"),
            "current_step": titration_steps[current_step_idx] if titration_steps else "unknown",
            "remaining_steps": remaining_steps,
            "titration_interval": interval,
            "total_remaining_titrations": len(remaining_steps),
            "estimated_time_to_target": f"{len(remaining_steps) * 2} weeks (minimum)" if remaining_steps else "At or near target",
            "monitoring_at_each_step": monitoring,
            "precautions": precautions,
            "key_trial": db_entry.get("key_trial", ""),
        }

    # ──────────────────────────────────────────────────────────────────────
    # CONTRAINDICATION CHECKING
    # ──────────────────────────────────────────────────────────────────────

    def _check_contraindications(
        self, drug: str, patient_data: Dict[str, Any],
    ) -> List[str]:
        """Check drug-specific contraindications against patient data.

        Returns a list of active contraindications found.
        """
        active: List[str] = []
        normalised = _normalise_drug_name(drug)
        alias = _DRUG_ALIASES.get(normalised)

        if alias is None:
            return active

        pillar, canonical = alias
        db_entry = self._drug_db.get(pillar, {}).get(canonical, {})
        if not db_entry:
            # Check additional therapies
            db_entry = self._additional.get(normalised, {})

        contraindications = db_entry.get("contraindications", [])
        patient_conditions: List[str] = patient_data.get("contraindications", [])
        patient_conditions_lower = [c.lower() for c in patient_conditions]

        for ci in contraindications:
            ci_lower = ci.lower()
            # Check against explicit patient contraindications
            for pc in patient_conditions_lower:
                if pc in ci_lower or ci_lower in pc:
                    active.append(f"{canonical}: {ci}")
                    break

        # Lab-based contraindication checks
        k = patient_data.get("potassium")
        egfr = patient_data.get("egfr")
        hr = patient_data.get("hr")
        sbp = patient_data.get("sbp")

        if pillar == GDMTPillar.MRA:
            if k is not None and k > self._MAX_K_FOR_MRA:
                active.append(
                    f"{canonical}: Potassium {k} mEq/L > {self._MAX_K_FOR_MRA}. "
                    f"MRA contraindicated until K+ normalised."
                )
            if egfr is not None and egfr < self._MIN_EGFR_FOR_MRA:
                active.append(
                    f"{canonical}: eGFR {egfr} mL/min/1.73m2 < {self._MIN_EGFR_FOR_MRA}. "
                    f"MRA contraindicated."
                )

        if pillar == GDMTPillar.ARNI_ACEI_ARB:
            if sbp is not None and sbp < self._MIN_SBP_ABSOLUTE:
                active.append(
                    f"{canonical}: SBP {sbp} mmHg < {self._MIN_SBP_ABSOLUTE}. "
                    f"ARNI/ACEi/ARB contraindicated due to hypotension."
                )
            if egfr is not None and egfr < self._MIN_EGFR_CAUTION_ARNI:
                active.append(
                    f"{canonical}: eGFR {egfr} mL/min/1.73m2 < {self._MIN_EGFR_CAUTION_ARNI}. "
                    f"Use with extreme caution or consider alternative."
                )

        if pillar == GDMTPillar.BETA_BLOCKER:
            if hr is not None and hr < self._MIN_HR_FOR_BB:
                active.append(
                    f"{canonical}: HR {hr} bpm < {self._MIN_HR_FOR_BB}. "
                    f"Beta-blocker contraindicated due to bradycardia."
                )
            if sbp is not None and sbp < self._MIN_SBP_ABSOLUTE:
                active.append(
                    f"{canonical}: SBP {sbp} mmHg < {self._MIN_SBP_ABSOLUTE}. "
                    f"Hold beta-blocker."
                )

        return active

    # ──────────────────────────────────────────────────────────────────────
    # DRUG INTERACTION CHECKING
    # ──────────────────────────────────────────────────────────────────────

    def _check_drug_interactions(
        self, drug: str, current_meds: List[str],
    ) -> List[str]:
        """Check for clinically significant drug interactions.

        Returns warnings for each identified interaction.
        """
        warnings: List[str] = []
        normalised = _normalise_drug_name(drug)

        # Resolve to canonical name via aliases
        alias = _DRUG_ALIASES.get(normalised)
        canonical = alias[1] if alias else normalised

        interactions = self._interactions.get(canonical, [])
        current_normalised = {_normalise_drug_name(m) for m in current_meds}

        # Also resolve current med names to canonical
        current_canonical: set[str] = set()
        for m in current_normalised:
            a = _DRUG_ALIASES.get(m)
            if a:
                current_canonical.add(a[1])
            else:
                current_canonical.add(m)

        for interaction in interactions:
            interacting = interaction["interacting_drug"]
            if interacting in current_canonical or interacting in current_normalised:
                severity = interaction["severity"].upper()
                effect = interaction["effect"]
                warnings.append(
                    f"INTERACTION [{severity}]: {canonical} + {interacting} - {effect}"
                )

        return warnings

    def _check_all_interactions(
        self,
        all_med_names: List[str],
        rec: GDMTRecommendation,
    ) -> List[str]:
        """Check interactions for all current + newly recommended medications."""
        warnings: List[str] = []
        seen: set[str] = set()

        # Extract newly recommended drug names from recommendations text
        new_drug_names: List[str] = []
        for r in rec.recommendations:
            for alias_key in _DRUG_ALIASES:
                if alias_key in r.lower():
                    new_drug_names.append(alias_key)

        all_names = all_med_names + new_drug_names

        for drug in all_names:
            interactions = self._check_drug_interactions(drug, all_names)
            for w in interactions:
                if w not in seen:
                    seen.add(w)
                    warnings.append(w)

        return warnings

    # ──────────────────────────────────────────────────────────────────────
    # ADDITIONAL THERAPIES ASSESSMENT
    # ──────────────────────────────────────────────────────────────────────

    def _assess_additional_therapies(
        self,
        ef_category: EjectionFractionCategory,
        nyha: HeartFailureClass,
        patient_data: Dict[str, Any],
        current_medications: List[Dict[str, Any]],
    ) -> List[str]:
        """Assess eligibility for therapies beyond the 4 pillars.

        Checks hydralazine-ISDN, ivabradine, digoxin, and vericiguat
        eligibility based on patient characteristics.
        """
        recs: List[str] = []
        current_names = {_normalise_drug_name(m.get("name", "")) for m in current_medications}
        current_canonical: set[str] = set()
        for n in current_names:
            a = _DRUG_ALIASES.get(n)
            current_canonical.add(a[1] if a else n)

        hr = patient_data.get("hr")
        race = patient_data.get("race", "").lower()
        rhythm = patient_data.get("rhythm", "").lower()

        is_hfref = ef_category in (
            EjectionFractionCategory.HFrEF, EjectionFractionCategory.HFimpEF,
        )

        # ── Hydralazine-ISDN ────────────────────────────────────────────
        if is_hfref:
            already_on_h_isdn = any(
                "hydralazine" in n or "isosorbide" in n for n in current_canonical
            )
            if not already_on_h_isdn:
                # Class I for self-identified African American with NYHA III-IV
                if race in ("african_american", "african american", "black") and nyha in (
                    HeartFailureClass.NYHA_III, HeartFailureClass.NYHA_IV,
                ):
                    recs.append(
                        "ADDITIONAL THERAPY: Hydralazine + Isosorbide dinitrate (Class I, LOE A). "
                        "A-HeFT trial demonstrated 43% mortality reduction in self-identified "
                        "African American patients with NYHA III-IV HFrEF. Start hydralazine "
                        "25 mg + ISDN 20 mg TID, titrate to hydralazine 75 mg + ISDN 40 mg TID."
                    )

                # Class IIa if intolerant to ACEi/ARB/ARNI
                arni_pillar_med = None
                for med in current_medications:
                    n = _normalise_drug_name(med.get("name", ""))
                    a = _DRUG_ALIASES.get(n)
                    if a and a[0] == GDMTPillar.ARNI_ACEI_ARB:
                        arni_pillar_med = med
                        break

                if arni_pillar_med is None:
                    status = ""
                    for ci_str in patient_data.get("contraindications", []):
                        ci_lower = ci_str.lower()
                        if "acei" in ci_lower or "arb" in ci_lower or "arni" in ci_lower or "angioedema" in ci_lower:
                            status = "intolerant"
                            break
                    if status == "intolerant":
                        recs.append(
                            "ADDITIONAL THERAPY: Consider hydralazine + isosorbide dinitrate "
                            "as ACEi/ARB/ARNI alternative (Class IIa, LOE B-R). Patient appears "
                            "intolerant to RAAS inhibitors."
                        )

        # ── Ivabradine ──────────────────────────────────────────────────
        ef_val = float(patient_data.get("lvef", patient_data.get("ef", 50)))
        if (is_hfref and ef_val <= self._IVABRADINE_EF_CUTOFF
                and "ivabradine" not in current_canonical):
            if hr is not None and hr >= self._IVABRADINE_HR_THRESHOLD:
                if "sinus" in rhythm:
                    if nyha in (HeartFailureClass.NYHA_II, HeartFailureClass.NYHA_III):
                        # Check if already on maximally tolerated beta-blocker
                        bb_on_board = any(
                            _DRUG_ALIASES.get(n, (None, None))[0] == GDMTPillar.BETA_BLOCKER
                            for n in current_canonical
                        )
                        recs.append(
                            f"ADDITIONAL THERAPY: Ivabradine (Class IIa, LOE B-R). "
                            f"Resting HR is {hr} bpm (>=70) "
                            f"{'on beta-blocker' if bb_on_board else '(ensure beta-blocker maximised first)'}. "
                            f"SHIFT trial: 18% reduction in CV death or HF hospitalization. "
                            f"Start 2.5-5 mg BID, titrate to 7.5 mg BID. "
                            f"Requires sinus rhythm; contraindicated in AF."
                        )

        # ── Vericiguat ──────────────────────────────────────────────────
        if is_hfref and "vericiguat" not in current_canonical:
            recent_hf_hosp = patient_data.get("recent_hf_hospitalization", False)
            if recent_hf_hosp and nyha in (
                HeartFailureClass.NYHA_II, HeartFailureClass.NYHA_III, HeartFailureClass.NYHA_IV,
            ):
                recs.append(
                    "ADDITIONAL THERAPY: Consider vericiguat (Class IIb, LOE B-R). "
                    "Patient has recent HF hospitalization despite optimized GDMT. "
                    "VICTORIA trial: 10% reduction in CV death or first HF hospitalization. "
                    "Start 2.5 mg daily, titrate to 10 mg daily every 2 weeks. "
                    "Contraindicated with concomitant nitrates or PDE-5 inhibitors."
                )

        # ── Digoxin ─────────────────────────────────────────────────────
        if is_hfref and "digoxin" not in current_canonical:
            if nyha in (HeartFailureClass.NYHA_III, HeartFailureClass.NYHA_IV):
                recs.append(
                    "ADDITIONAL THERAPY: Consider digoxin (Class IIb, LOE B-R). "
                    "For persistent symptoms despite optimized GDMT. DIG trial: no mortality "
                    "benefit but 28% reduction in HF hospitalizations. Target serum level "
                    "0.5-0.9 ng/mL. Monitor K+, Mg2+, renal function. "
                    "Caution with amiodarone (reduce dose 50%) and verapamil."
                )

        return recs

    # ──────────────────────────────────────────────────────────────────────
    # DEVICE THERAPY ASSESSMENT
    # ──────────────────────────────────────────────────────────────────────

    def _assess_device_therapy(
        self,
        lvef: float,
        nyha: HeartFailureClass,
        qrs_duration: Optional[int],
        qrs_morphology: Optional[str],
        patient_data: Dict[str, Any],
    ) -> List[str]:
        """Assess eligibility for ICD, CRT, or CRT-D based on guidelines.

        Parameters
        ----------
        lvef : float
            Current LVEF (%).
        nyha : HeartFailureClass
            Current NYHA functional class.
        qrs_duration : int or None
            QRS duration in ms (from ECG).
        qrs_morphology : str or None
            "LBBB", "RBBB", "IVCD", or "normal".
        patient_data : dict
            Additional patient context.
        """
        recs: List[str] = []
        days_post_mi = patient_data.get("days_post_mi")
        days_on_gdmt = patient_data.get("days_on_gdmt")
        expected_survival = patient_data.get("expected_survival_years", 2)

        nyha_eligible_icd = nyha in (HeartFailureClass.NYHA_II, HeartFailureClass.NYHA_III)
        nyha_eligible_crt = nyha in (
            HeartFailureClass.NYHA_II, HeartFailureClass.NYHA_III, HeartFailureClass.NYHA_IV,
        )

        # ── ICD assessment ──────────────────────────────────────────────
        icd_eligible = True
        icd_reasons: List[str] = []
        icd_barriers: List[str] = []

        if lvef > self._ICD_EF_CUTOFF:
            icd_eligible = False
            icd_barriers.append(f"LVEF {lvef}% > {self._ICD_EF_CUTOFF}%")

        if not nyha_eligible_icd:
            icd_eligible = False
            icd_barriers.append(f"NYHA class {nyha.value} not in II-III range")

        if days_post_mi is not None and days_post_mi < 40:
            icd_eligible = False
            icd_barriers.append(
                f"Only {days_post_mi} days post-MI (must wait >=40 days)"
            )

        if days_on_gdmt is not None and days_on_gdmt < 90:
            icd_eligible = False
            icd_barriers.append(
                f"Only {days_on_gdmt} days on optimal GDMT (must wait >=90 days)"
            )

        if expected_survival is not None and expected_survival < 1:
            icd_eligible = False
            icd_barriers.append("Expected survival <1 year")

        if icd_eligible:
            icd_reasons.append(
                f"LVEF {lvef}% <=35%, NYHA {nyha.value}"
            )
            recs.append(
                f"DEVICE THERAPY - ICD (Class I, LOE A): Patient meets primary prevention "
                f"ICD criteria. LVEF {lvef}% <=35%, NYHA {nyha.value}. "
                f"SCD-HeFT: 23% mortality reduction. "
                f"Ensure >=90 days on optimal GDMT and reassess EF before implant."
            )
        elif lvef <= self._ICD_EF_CUTOFF:
            barrier_str = "; ".join(icd_barriers)
            recs.append(
                f"DEVICE THERAPY - ICD: Not currently eligible ({barrier_str}). "
                f"Reassess when barriers resolved."
            )

        # ── CRT assessment ──────────────────────────────────────────────
        crt_eligible = False
        crt_class = ""

        if lvef <= self._CRT_EF_CUTOFF and nyha_eligible_crt and qrs_duration is not None:
            morphology = (qrs_morphology or "").upper()
            is_lbbb = "LBBB" in morphology or "LEFT BUNDLE" in morphology

            if qrs_duration >= self._CRT_QRS_CLASS_I and is_lbbb:
                crt_eligible = True
                crt_class = "Class I, LOE A"
            elif qrs_duration >= self._CRT_QRS_CLASS_IIA and is_lbbb:
                crt_eligible = True
                crt_class = "Class IIa, LOE B-R"
            elif qrs_duration >= self._CRT_QRS_CLASS_I and not is_lbbb:
                crt_eligible = True
                crt_class = "Class IIa, LOE B-NR"
            elif qrs_duration >= self._CRT_QRS_CLASS_IIA and not is_lbbb:
                crt_eligible = True
                crt_class = "Class IIb, LOE B-NR"

        if crt_eligible:
            morphology_str = qrs_morphology or "not specified"
            recs.append(
                f"DEVICE THERAPY - CRT ({crt_class}): Patient meets CRT criteria. "
                f"LVEF {lvef}%, QRS {qrs_duration} ms ({morphology_str}), NYHA {nyha.value}. "
                f"CARE-HF: 36% mortality reduction. Greatest benefit with LBBB QRS >=150 ms."
            )

            # CRT-D if also ICD eligible
            if icd_eligible:
                recs.append(
                    f"DEVICE THERAPY - CRT-D ({crt_class}): Meets both ICD and CRT criteria. "
                    f"CRT-D preferred over CRT-P when ICD indication coexists. "
                    f"COMPANION: 36% all-cause mortality reduction with CRT-D."
                )
        elif lvef <= self._CRT_EF_CUTOFF and qrs_duration is not None and qrs_duration < self._CRT_QRS_CLASS_IIA:
            recs.append(
                f"DEVICE THERAPY - CRT: QRS duration {qrs_duration} ms < {self._CRT_QRS_CLASS_IIA} ms. "
                f"CRT not indicated. Continue optimal GDMT."
            )

        return recs

    # ──────────────────────────────────────────────────────────────────────
    # MONITORING PLAN
    # ──────────────────────────────────────────────────────────────────────

    def _generate_monitoring_plan(
        self, rec: GDMTRecommendation, patient_data: Dict[str, Any],
    ) -> List[str]:
        """Generate a monitoring plan based on active and recommended therapies.

        Returns a list of monitoring action items.
        """
        plan: List[str] = []

        # General monitoring for all HF patients
        plan.append(
            "MONITORING: Daily weights. Contact provider if >2 lb gain in 24 hours "
            "or >5 lb gain in 1 week."
        )
        plan.append(
            "MONITORING: Sodium restriction (<2 g/day) and fluid restriction if "
            "hyponatraemic (Na <130 mEq/L)."
        )

        # Medication-specific monitoring
        active_pillars: set[str] = set()
        for med in rec.current_meds:
            active_pillars.add(med.pillar.value)

        has_raas = GDMTPillar.ARNI_ACEI_ARB.value in active_pillars
        has_mra = GDMTPillar.MRA.value in active_pillars
        has_bb = GDMTPillar.BETA_BLOCKER.value in active_pillars

        if has_raas or has_mra:
            plan.append(
                "MONITORING - LABS: Check BMP (potassium, creatinine, eGFR) within "
                "1-2 weeks of initiating or titrating RAAS inhibitor or MRA, then "
                "at 1 month, then every 3-6 months."
            )

        if has_mra:
            plan.append(
                "MONITORING - MRA: Serum K+ must remain <5.0 mEq/L. Hold MRA if "
                "K+ >=5.5 mEq/L. Recheck K+ within 3-7 days of any dose change."
            )

        if has_bb:
            plan.append(
                "MONITORING - BETA-BLOCKER: Check HR and BP before each titration. "
                "Target resting HR 60-70 bpm. Hold if HR <55 bpm or SBP <90 mmHg."
            )

        if has_raas:
            plan.append(
                "MONITORING - ARNI/ACEi/ARB: Check BP before each titration. "
                "Preferred SBP >=100 mmHg for uptitration. Hold if SBP <90 mmHg. "
                "Monitor for angioedema, cough (ACEi), dizziness."
            )

        # Check for newly recommended therapies and add monitoring
        rec_text = " ".join(rec.recommendations).lower()
        if "sglt2" in rec_text or GDMTPillar.SGLT2I.value in active_pillars:
            plan.append(
                "MONITORING - SGLT2i: Check eGFR at baseline and periodically. "
                "Expected initial eGFR dip of 10-15% (stabilises). Monitor for "
                "genital mycotic infections. Check blood glucose in diabetics."
            )

        if "ivabradine" in rec_text:
            plan.append(
                "MONITORING - IVABRADINE: Check HR before each titration. "
                "Reduce dose if HR <50 bpm. Discontinue if AF develops. "
                "Ask about visual disturbances (phosphenes)."
            )

        if "digoxin" in rec_text:
            plan.append(
                "MONITORING - DIGOXIN: Check serum digoxin level 6+ hours post-dose "
                "(target 0.5-0.9 ng/mL). Monitor K+, Mg2+, renal function. "
                "Toxicity signs: nausea, visual changes, arrhythmias."
            )

        if "vericiguat" in rec_text:
            plan.append(
                "MONITORING - VERICIGUAT: Check BP before each titration. "
                "Monitor haemoglobin (may worsen anaemia). "
                "Hold if SBP <100 mmHg."
            )

        # NT-proBNP / BNP trending
        plan.append(
            "MONITORING - BIOMARKERS: Obtain NT-proBNP or BNP at baseline and with "
            "clinical changes. Trending biomarkers can guide therapy intensification "
            "(GUIDE-IT approach)."
        )

        return plan

    # ──────────────────────────────────────────────────────────────────────
    # HFrEF OPTIMIZATION
    # ──────────────────────────────────────────────────────────────────────

    def _optimize_hfref(
        self,
        rec: GDMTRecommendation,
        ef_category: EjectionFractionCategory,
        nyha: HeartFailureClass,
        current_gdmt: Dict[GDMTPillar, GDMTMedication],
        patient_data: Dict[str, Any],
        contraindications: List[str],
    ) -> None:
        """Generate HFrEF-specific GDMT optimization recommendations.

        Key principle: initiate all 4 pillars simultaneously.  Do NOT wait
        for one pillar to reach target dose before starting the next.
        """
        rec.recommendations.append(
            "HFrEF MANAGEMENT: ACC/AHA 2022 recommends initiation of all 4 GDMT pillars "
            "(beta-blocker, ARNI/ACEi/ARB, MRA, SGLT2i) simultaneously in eligible patients. "
            "Do NOT wait for one to reach target dose before starting others."
        )

        if ef_category == EjectionFractionCategory.HFimpEF:
            rec.recommendations.append(
                "HFimpEF: LVEF has improved from prior measurement. CONTINUE all GDMT. "
                "Withdrawal of therapies has been associated with recurrent LV dysfunction "
                "and clinical deterioration (TRED-HF trial)."
            )

        # Identify gaps
        gaps = self._identify_gaps(ef_category, current_gdmt, contraindications)
        rec.recommendations.extend(gaps)

        # For each missing pillar, recommend initiation
        for pillar in (GDMTPillar.BETA_BLOCKER, GDMTPillar.ARNI_ACEI_ARB, GDMTPillar.MRA, GDMTPillar.SGLT2I):
            if pillar not in current_gdmt:
                self._recommend_pillar_initiation(
                    rec, pillar, patient_data, contraindications,
                )

        # For each existing sub-target pillar, generate titration plan
        for pillar, med in current_gdmt.items():
            if med.status in (GDMTStatus.INITIATED, GDMTStatus.UPTITRATING):
                # Check for contraindications before titrating
                ci = self._check_contraindications(med.drug_name, patient_data)
                if ci:
                    rec.next_steps.append(
                        f"CAUTION before titrating {med.drug_name}: " + "; ".join(ci)
                    )
                else:
                    plan = self._generate_titration_plan(med, patient_data)
                    if plan.get("remaining_steps"):
                        rec.next_steps.append(
                            f"TITRATE {med.drug_name}: {med.current_dose} -> "
                            f"{plan['target_dose']}. Next step: "
                            f"{plan['remaining_steps'][0]}. Interval: "
                            f"{plan['titration_interval']}."
                        )

    def _recommend_pillar_initiation(
        self,
        rec: GDMTRecommendation,
        pillar: GDMTPillar,
        patient_data: Dict[str, Any],
        contraindications: List[str],
    ) -> None:
        """Recommend initiating a specific GDMT pillar with preferred agent."""
        sbp = patient_data.get("sbp")
        hr = patient_data.get("hr")
        k = patient_data.get("potassium")
        egfr = patient_data.get("egfr")
        weight = patient_data.get("weight_kg")

        if pillar == GDMTPillar.BETA_BLOCKER:
            # Check haemodynamic eligibility
            if hr is not None and hr < self._MIN_HR_FOR_BB:
                rec.recommendations.append(
                    f"BETA-BLOCKER: HR {hr} bpm is below {self._MIN_HR_FOR_BB}. "
                    f"Defer initiation until HR supports it."
                )
                return
            if sbp is not None and sbp < self._MIN_SBP_ABSOLUTE:
                rec.recommendations.append(
                    f"BETA-BLOCKER: SBP {sbp} mmHg too low. Defer initiation."
                )
                return

            rec.recommendations.append(
                f"INITIATE BETA-BLOCKER: Start carvedilol 3.125 mg BID "
                f"(Class I, LOE A, COPERNICUS trial). Titrate every 2 weeks to "
                f"target 25 mg BID{' (50 mg BID if >85 kg)' if weight and weight > 85 else ''}. "
                f"Monitor HR (target 60-70 bpm) and BP at each visit. "
                f"Alternatives: metoprolol succinate (MERIT-HF), bisoprolol (CIBIS-II)."
            )

        elif pillar == GDMTPillar.ARNI_ACEI_ARB:
            if sbp is not None and sbp < self._MIN_SBP_FOR_RAAS:
                rec.recommendations.append(
                    f"ARNI/ACEi/ARB: SBP {sbp} mmHg below {self._MIN_SBP_FOR_RAAS}. "
                    f"Use low starting dose with careful monitoring. Consider reducing "
                    f"diuretic if patient is euvolemic."
                )

            # Check for angioedema history
            angioedema_hx = any("angioedema" in c.lower() for c in contraindications)
            if angioedema_hx:
                rec.recommendations.append(
                    "ARNI/ACEi/ARB: History of angioedema. ARNI and ACEi are "
                    "contraindicated. Consider ARB with close monitoring or "
                    "hydralazine-ISDN alternative."
                )
                return

            # Check if currently on ACEi (need washout before ARNI)
            current_on_acei = False
            for med in rec.current_meds:
                if med.drug_name in ("enalapril", "lisinopril", "ramipril"):
                    current_on_acei = True
                    break

            if current_on_acei:
                rec.recommendations.append(
                    "SWITCH TO ARNI: Currently on ACEi. Discontinue ACEi and wait "
                    "36 hours before starting sacubitril/valsartan 24/26 mg BID "
                    "(PARADIGM-HF: 20% CV death reduction vs enalapril). "
                    "36-hour washout is mandatory to reduce angioedema risk."
                )
            else:
                starting = "49/51 mg BID" if (sbp and sbp >= 110) else "24/26 mg BID"
                rec.recommendations.append(
                    f"INITIATE ARNI: Start sacubitril/valsartan {starting} "
                    f"(Class I, LOE A, PARADIGM-HF). Preferred over ACEi/ARB. "
                    f"Titrate every 2-4 weeks to target 97/103 mg BID. "
                    f"Monitor BP, K+, creatinine. "
                    f"If ARNI not available: start enalapril 2.5 mg BID (SOLVD-Treatment)."
                )

        elif pillar == GDMTPillar.MRA:
            if k is not None and k > self._MAX_K_FOR_MRA:
                rec.recommendations.append(
                    f"MRA: K+ is {k} mEq/L (>={self._MAX_K_FOR_MRA}). "
                    f"Defer MRA initiation until K+ <5.0. Address causes of hyperkalemia."
                )
                return

            if egfr is not None and egfr < self._MIN_EGFR_FOR_MRA:
                rec.recommendations.append(
                    f"MRA: eGFR {egfr} mL/min/1.73m2 < {self._MIN_EGFR_FOR_MRA}. "
                    f"MRA contraindicated due to high hyperkalemia risk."
                )
                return

            rec.recommendations.append(
                "INITIATE MRA: Start spironolactone 12.5-25 mg daily "
                "(Class I, LOE A, RALES trial: 30% mortality reduction). "
                "Check K+ within 3-7 days, then at 1 month, then every 3 months. "
                "Must maintain K+ <5.0 mEq/L and eGFR >=30. "
                "If gynecomastia occurs, switch to eplerenone 25 mg daily (EMPHASIS-HF)."
            )

        elif pillar == GDMTPillar.SGLT2I:
            rec.recommendations.append(
                "INITIATE SGLT2i: Start dapagliflozin 10 mg daily (DAPA-HF) or "
                "empagliflozin 10 mg daily (EMPEROR-Reduced). Class I, LOE A. "
                "No titration required. Benefits independent of diabetes status. "
                "Monitor eGFR (expect initial dip), volume status, and genital infections. "
                "May reduce diuretic dose if patient becomes volume-depleted."
            )

    # ──────────────────────────────────────────────────────────────────────
    # HFmrEF OPTIMIZATION
    # ──────────────────────────────────────────────────────────────────────

    def _optimize_hfmref(
        self,
        rec: GDMTRecommendation,
        nyha: HeartFailureClass,
        current_gdmt: Dict[GDMTPillar, GDMTMedication],
        patient_data: Dict[str, Any],
        contraindications: List[str],
    ) -> None:
        """Generate HFmrEF-specific recommendations.

        SGLT2i is Class I.  Other GDMT pillars may be reasonable (Class IIb)
        based on emerging evidence and clinical judgement.
        """
        rec.recommendations.append(
            "HFmrEF MANAGEMENT (LVEF 41-49%): SGLT2i is Class I (DELIVER, "
            "EMPEROR-Preserved included HFmrEF patients). Other GDMT pillars "
            "(beta-blocker, ARNI/ACEi/ARB, MRA) may be reasonable (Class IIb) "
            "based on clinical context - many patients with HFmrEF had prior HFrEF."
        )

        gaps = self._identify_gaps(EjectionFractionCategory.HFmrEF, current_gdmt, contraindications)
        rec.recommendations.extend(gaps)

        # Always recommend SGLT2i if not on it
        if GDMTPillar.SGLT2I not in current_gdmt:
            self._recommend_pillar_initiation(
                rec, GDMTPillar.SGLT2I, patient_data, contraindications,
            )

        # Manage congestion
        rec.recommendations.append(
            "DIURETICS: Use loop diuretics as needed for volume management. "
            "Target euvolemia. Adjust doses based on daily weights."
        )

        # Comorbidity management
        rec.recommendations.append(
            "COMORBIDITY MANAGEMENT: Address hypertension (target <130/80), AF, "
            "CAD, diabetes, and obesity. These drive much of HFmrEF pathology."
        )

    # ──────────────────────────────────────────────────────────────────────
    # HFpEF OPTIMIZATION
    # ──────────────────────────────────────────────────────────────────────

    def _optimize_hfpef(
        self,
        rec: GDMTRecommendation,
        nyha: HeartFailureClass,
        current_gdmt: Dict[GDMTPillar, GDMTMedication],
        patient_data: Dict[str, Any],
        contraindications: List[str],
    ) -> None:
        """Generate HFpEF-specific recommendations.

        SGLT2i is Class I.  Management is otherwise heavily comorbidity-driven.
        """
        rec.recommendations.append(
            "HFpEF MANAGEMENT (LVEF >=50%): SGLT2i is Class I (EMPEROR-Preserved, "
            "DELIVER). Management otherwise focuses on diuretics for congestion "
            "and aggressive comorbidity treatment."
        )

        # SGLT2i
        if GDMTPillar.SGLT2I not in current_gdmt:
            self._recommend_pillar_initiation(
                rec, GDMTPillar.SGLT2I, patient_data, contraindications,
            )

        # Diuretics for congestion
        rec.recommendations.append(
            "DIURETICS: Loop diuretics for volume management (Class I, LOE C-EO). "
            "Use lowest effective dose. Consider torsemide over furosemide "
            "(longer half-life, anti-aldosterone properties)."
        )

        # GLP-1 RA for obesity
        bmi = patient_data.get("bmi")
        if bmi is not None and bmi >= 30:
            rec.recommendations.append(
                f"GLP-1 RA: BMI {bmi} kg/m2 (>=30). Consider semaglutide for "
                f"obese HFpEF (Class IIa, LOE B-R). STEP-HFpEF: significant "
                f"improvement in symptoms, exercise capacity, and body weight. "
                f"Addresses obesity-driven inflammation and metabolic HFpEF phenotype."
            )

        # Comorbidity management
        comorbidity_recs = self._hfpef.get("treat_comorbidities", {})
        conditions = comorbidity_recs.get("conditions", {})

        if patient_data.get("hypertension"):
            ht = conditions.get("hypertension", {})
            rec.recommendations.append(
                f"HYPERTENSION: Target {ht.get('target', 'BP <130/80 mmHg')}. "
                f"{ht.get('agents', 'ACEi/ARB, thiazide, amlodipine')}."
            )

        if patient_data.get("atrial_fibrillation"):
            rec.recommendations.append(
                "ATRIAL FIBRILLATION: Rate vs rhythm control per clinical context. "
                "Anticoagulation per CHA2DS2-VASc score. AF very common in HFpEF "
                "and contributes significantly to symptom burden."
            )

        if patient_data.get("diabetes"):
            rec.recommendations.append(
                "DIABETES: SGLT2i preferred (dual benefit). Metformin is safe. "
                "AVOID thiazolidinediones (fluid retention worsens HF). "
                "GLP-1 RA has additional HFpEF benefit in obesity."
            )

        # Iron deficiency
        ferritin = patient_data.get("ferritin")
        tsat = patient_data.get("tsat")
        if ferritin is not None:
            iron_deficient = ferritin < 100 or (ferritin <= 300 and tsat is not None and tsat < 20)
            if iron_deficient:
                rec.recommendations.append(
                    f"IRON DEFICIENCY: Ferritin {ferritin}"
                    f"{f', TSAT {tsat}%' if tsat is not None else ''}. "
                    f"IV iron replacement recommended (Class IIa, LOE B-R). "
                    f"FAIR-HF and AFFIRM-AHF support IV ferric carboxymaltose. "
                    f"Improves symptoms, functional capacity, and quality of life."
                )

    # ──────────────────────────────────────────────────────────────────────
    # UTILITY: SUMMARY
    # ──────────────────────────────────────────────────────────────────────

    def get_drug_info(self, drug_name: str) -> Optional[Dict[str, Any]]:
        """Retrieve full database entry for a drug by name or alias.

        Returns None if not found in GDMT database or additional therapies.
        """
        normalised = _normalise_drug_name(drug_name)
        alias = _DRUG_ALIASES.get(normalised)

        if alias:
            pillar, canonical = alias
            entry = self._drug_db.get(pillar, {}).get(canonical)
            if entry:
                return {"pillar": pillar.value, "drug": canonical, **entry}

        # Check additional therapy aliases first
        resolved = _ADDITIONAL_THERAPY_ALIASES.get(normalised, normalised)

        # Check additional therapies
        for therapy_key, therapy_data in self._additional.items():
            if resolved in therapy_key or therapy_key in resolved:
                return {"category": "additional_therapy", "drug": therapy_key, **therapy_data}

        return None

    def get_all_pillars_status(
        self, current_medications: List[Dict[str, Any]],
    ) -> Dict[str, Dict[str, str]]:
        """Return a summary of all 4 pillars and their current status.

        Useful for dashboard display or quick clinical overview.
        """
        current_gdmt = self._assess_current_gdmt(current_medications)

        summary: Dict[str, Dict[str, str]] = {}
        for pillar in GDMTPillar:
            if pillar in current_gdmt:
                med = current_gdmt[pillar]
                summary[pillar.value] = {
                    "status": med.status.value,
                    "drug": med.drug_name,
                    "current_dose": med.current_dose or "unknown",
                    "target_dose": med.target_dose,
                }
            else:
                summary[pillar.value] = {
                    "status": "not_started",
                    "drug": "none",
                    "current_dose": "N/A",
                    "target_dose": "N/A",
                }

        return summary

    def get_guideline_reference(self) -> str:
        """Return the primary guideline citation."""
        return self.GUIDELINE_REFERENCE
