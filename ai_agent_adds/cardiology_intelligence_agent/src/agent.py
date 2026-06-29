"""Cardiology Intelligence Agent -- autonomous reasoning across cardiovascular data silos.

Implements the plan -> search -> evaluate -> synthesize -> report pattern from the
VAST AI OS AgentEngine model. The agent can:

1. Parse complex multi-part questions about cardiovascular medicine
2. Plan a search strategy across 12 domain-specific collections
3. Execute multi-collection retrieval via the CardioRAGEngine
4. Evaluate evidence quality and completeness
5. Synthesize cross-functional insights with clinical alerts
6. Generate structured reports with cardiology-specific formatting

Mapping to VAST AI OS:
  - AgentEngine entry point: CardioIntelligenceAgent.run()
  - Plan -> search_plan()
  - Execute -> rag_engine.query()
  - Reflect -> evaluate_evidence()
  - Report -> generate_report()

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, List

from .models import (
    CardioWorkflowType,
    WorkflowResult,
)


# =====================================================================
# CARDIOLOGY SYSTEM PROMPT
# =====================================================================

CARDIO_SYSTEM_PROMPT = """\
You are a cardiology clinical decision support system embedded within the HCLS AI Factory \
precision medicine platform. You have deep expertise across the full spectrum of \
cardiovascular medicine including coronary artery disease, heart failure, valvular \
heart disease, arrhythmias, cardiac imaging, interventional cardiology, preventive \
cardiology, and cardio-oncology.

Your responses must adhere to the following clinical standards:

1. **Guideline Citations** -- Always cite ACC/AHA/ESC/HRS/SCAI guidelines with the \
   Class of Recommendation (I, IIa, IIb, III) and Level of Evidence (A, B-R, B-NR, \
   C-LD, C-EO). Use the format: [Class I, LOE A] or [Class IIa, LOE B-R]. Include \
   the guideline year and first author or writing committee when available.

2. **Evidence References** -- Include PubMed identifiers (PMID) for landmark trials \
   and key evidence supporting your recommendations. Format as clickable links: \
   [PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). Reference trial \
   acronyms (e.g., PARADIGM-HF, DAPA-HF, ISCHEMIA, PARTNER) where applicable.

3. **CRITICAL Findings** -- Flag the following as CRITICAL with prominent visual \
   markers and immediate action recommendations:
   - Acute ST-elevation myocardial infarction (STEMI)
   - Cardiac tamponade or pericardial effusion with hemodynamic compromise
   - Acute aortic dissection (Stanford Type A or Type B with complications)
   - Sustained ventricular tachycardia (VT) or ventricular fibrillation (VF)
   - Complete heart block with hemodynamic instability
   - Acute decompensated heart failure with cardiogenic shock
   - Massive pulmonary embolism with RV dysfunction
   - Acute severe aortic or mitral regurgitation
   - Hypertensive emergency with end-organ damage
   - Infective endocarditis with embolic complications

4. **Severity Badges** -- Classify all findings using standardised severity levels: \
   [CRITICAL], [HIGH], [MODERATE], [LOW], [BORDERLINE]. Place the badge at the \
   start of each finding or recommendation line.

5. **Measurements and Units** -- Include relevant measurements with SI and \
   conventional units where appropriate. Always provide the normal reference range \
   in parentheses. Example: LVEF 35% (normal >= 52% male, >= 54% female).

6. **Structured Formatting** -- Organise responses with clear sections: \
   Assessment, Risk Stratification, Recommendations, Monitoring Plan, and \
   Follow-up Timeline. Use bullet points and numbered lists for actionable items.

7. **Genomic Cross-Reference** -- When genetic or pharmacogenomic implications \
   are relevant (e.g., familial cardiomyopathy, channelopathies, PCSK9 variants, \
   CYP2C19 and clopidogrel), cross-reference with the genomic_evidence collection \
   and recommend appropriate genetic testing or genotype-guided therapy.

8. **Risk Calculators** -- Reference and interpret validated risk scores: \
   ASCVD Pooled Cohort Equations, HEART score, CHA2DS2-VASc, HAS-BLED, \
   MAGGIC, EuroSCORE II, SYNTAX score, STS risk calculator.

9. **GDMT Optimisation** -- For heart failure queries, systematically address \
   all four pillars of guideline-directed medical therapy (ARNI/ACEi/ARB, \
   beta-blocker, MRA, SGLT2i) with titration status and target doses.

10. **Limitations** -- You are a clinical decision support tool. You do NOT \
    replace physician judgment. All recommendations require clinical correlation \
    with the patient's complete history, physical examination, and diagnostic \
    workup. Explicitly state when evidence is limited or when specialist \
    referral is recommended."""


# =====================================================================
# WORKFLOW-SPECIFIC COLLECTION BOOST WEIGHTS
# =====================================================================
# Maps each CardioWorkflowType to collection weight overrides (multipliers).
# Collections not listed retain their base weight (1.0x). Values > 1.0
# boost the collection; values < 1.0 would suppress it.

WORKFLOW_COLLECTION_BOOST: Dict[CardioWorkflowType, Dict[str, float]] = {

    # ── Coronary Artery Disease Assessment ────────────────────────────
    CardioWorkflowType.CAD_ASSESSMENT: {
        "cardio_imaging": 2.0,
        "cardio_interventional": 1.8,
        "cardio_guidelines": 1.5,
        "cardio_trials": 1.3,
        "cardio_hemodynamics": 1.2,
        "cardio_prevention": 1.2,
        "cardio_literature": 1.1,
        "cardio_electrophysiology": 1.0,
    },

    # ── Heart Failure Evaluation ──────────────────────────────────────
    CardioWorkflowType.HEART_FAILURE: {
        "cardio_heart_failure": 2.5,
        "cardio_guidelines": 1.5,
        "cardio_trials": 1.5,
        "cardio_imaging": 1.3,
        "cardio_hemodynamics": 1.3,
        "cardio_devices": 1.2,
        "cardio_literature": 1.1,
        "cardio_electrophysiology": 1.0,
    },

    # ── Valvular Heart Disease ────────────────────────────────────────
    CardioWorkflowType.VALVULAR_DISEASE: {
        "cardio_valvular": 2.5,
        "cardio_imaging": 2.0,
        "cardio_interventional": 1.8,
        "cardio_guidelines": 1.5,
        "cardio_hemodynamics": 1.5,
        "cardio_trials": 1.2,
        "cardio_literature": 1.1,
    },

    # ── Arrhythmia / Electrophysiology ────────────────────────────────
    CardioWorkflowType.ARRHYTHMIA: {
        "cardio_electrophysiology": 2.5,
        "cardio_guidelines": 1.5,
        "cardio_devices": 1.5,
        "cardio_trials": 1.3,
        "cardio_literature": 1.2,
        "cardio_imaging": 1.1,
        "cardio_interventional": 1.0,
    },

    # ── Cardiac MRI Interpretation ────────────────────────────────────
    CardioWorkflowType.CARDIAC_MRI: {
        "cardio_imaging": 2.5,
        "cardio_heart_failure": 1.5,
        "cardio_guidelines": 1.3,
        "cardio_valvular": 1.2,
        "cardio_hemodynamics": 1.2,
        "cardio_literature": 1.1,
    },

    # ── Stress Testing ────────────────────────────────────────────────
    CardioWorkflowType.STRESS_TEST: {
        "cardio_imaging": 2.0,
        "cardio_guidelines": 1.8,
        "cardio_interventional": 1.5,
        "cardio_prevention": 1.3,
        "cardio_electrophysiology": 1.2,
        "cardio_trials": 1.1,
        "cardio_literature": 1.0,
    },

    # ── Preventive Cardiology / Risk Assessment ───────────────────────
    CardioWorkflowType.PREVENTIVE_RISK: {
        "cardio_prevention": 2.5,
        "cardio_guidelines": 1.8,
        "cardio_trials": 1.5,
        "cardio_literature": 1.3,
        "cardio_imaging": 1.0,
        "cardio_devices": 0.8,
    },

    # ── Cardio-Oncology Consultation ──────────────────────────────────
    CardioWorkflowType.CARDIO_ONCOLOGY: {
        "cardio_oncology": 2.5,
        "cardio_imaging": 1.5,
        "cardio_guidelines": 1.3,
        "cardio_trials": 1.3,
        "cardio_heart_failure": 1.2,
        "cardio_literature": 1.2,
        "cardio_electrophysiology": 1.0,
    },

    # ── General cardiology (balanced across all collections) ─────────
    CardioWorkflowType.GENERAL: {
        "cardio_literature": 1.2,
        "cardio_guidelines": 1.2,
        "cardio_trials": 1.1,
        "cardio_imaging": 1.0,
        "cardio_heart_failure": 1.0,
        "cardio_electrophysiology": 1.0,
        "cardio_prevention": 1.0,
        "cardio_valvular": 1.0,
        "cardio_interventional": 1.0,
        "cardio_hemodynamics": 1.0,
        "cardio_oncology": 0.9,
        "cardio_devices": 0.9,
    },

    # ── Acute Decompensated Heart Failure ──────────────────────────────
    CardioWorkflowType.ACUTE_DECOMPENSATED_HF: {
        "cardio_heart_failure": 2.5,
        "cardio_hemodynamics": 2.0,
        "cardio_guidelines": 1.5,
        "cardio_imaging": 1.3,
        "cardio_devices": 1.2,
        "cardio_trials": 1.2,
        "cardio_literature": 1.1,
        "cardio_electrophysiology": 1.0,
    },

    # ── Post-MI Secondary Prevention ──────────────────────────────────
    CardioWorkflowType.POST_MI: {
        "cardio_guidelines": 2.0,
        "cardio_interventional": 2.0,
        "cardio_imaging": 1.5,
        "cardio_trials": 1.5,
        "cardio_prevention": 1.3,
        "cardio_heart_failure": 1.2,
        "cardio_electrophysiology": 1.1,
        "cardio_literature": 1.0,
        "cardio_devices": 1.0,
    },

    # ── Myocarditis / Pericarditis ────────────────────────────────────
    CardioWorkflowType.MYOCARDITIS_PERICARDITIS: {
        "cardio_imaging": 2.5,
        "cardio_guidelines": 1.8,
        "cardio_literature": 1.5,
        "cardio_heart_failure": 1.3,
        "cardio_trials": 1.2,
        "cardio_hemodynamics": 1.1,
        "cardio_electrophysiology": 1.0,
    },
}


# =====================================================================
# KNOWLEDGE DOMAIN DICTIONARIES
# =====================================================================
# Comprehensive clinical knowledge for entity detection and context
# enrichment. These are used by the agent's search_plan() to identify
# clinical entities in user queries and map them to workflows.

CARDIO_CONDITIONS: Dict[str, Dict[str, object]] = {
    # ── Coronary Artery Disease ───────────────────────────────────────
    "coronary artery disease": {
        "aliases": ["cad", "ischemic heart disease", "ihd", "coronary disease",
                    "atherosclerotic heart disease", "ashd"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
        "collections": ["cardio_imaging", "cardio_interventional", "cardio_guidelines"],
    },
    "acute coronary syndrome": {
        "aliases": ["acs", "nstemi", "nste-acs", "unstable angina", "ua"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
        "collections": ["cardio_interventional", "cardio_guidelines"],
    },
    "stemi": {
        "aliases": ["st elevation myocardial infarction", "st-elevation mi",
                    "acute mi", "myocardial infarction"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
        "collections": ["cardio_interventional", "cardio_guidelines"],
        "critical": True,
    },
    "stable angina": {
        "aliases": ["chronic stable angina", "exertional angina", "angina pectoris",
                    "stable ischemic heart disease", "sihd"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT, CardioWorkflowType.STRESS_TEST],
        "collections": ["cardio_imaging", "cardio_guidelines"],
    },

    # ── Heart Failure ─────────────────────────────────────────────────
    "heart failure": {
        "aliases": ["hf", "chf", "congestive heart failure", "cardiac failure"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_heart_failure", "cardio_guidelines"],
    },
    "hfref": {
        "aliases": ["heart failure with reduced ejection fraction",
                    "systolic heart failure", "hf with reduced ef"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_heart_failure", "cardio_guidelines"],
    },
    "hfpef": {
        "aliases": ["heart failure with preserved ejection fraction",
                    "diastolic heart failure", "hf with preserved ef"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_heart_failure", "cardio_guidelines"],
    },
    "hfmref": {
        "aliases": ["heart failure with mildly reduced ejection fraction",
                    "hf mid-range", "hf with mildly reduced ef"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_heart_failure", "cardio_guidelines"],
    },
    "cardiogenic shock": {
        "aliases": ["cs", "cardiogenic shock syndrome"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_heart_failure", "cardio_hemodynamics"],
        "critical": True,
    },
    "cardiomyopathy": {
        "aliases": ["dcm", "dilated cardiomyopathy", "hcm",
                    "hypertrophic cardiomyopathy", "rcm",
                    "restrictive cardiomyopathy", "arvc", "arvd",
                    "arrhythmogenic right ventricular cardiomyopathy",
                    "takotsubo", "stress cardiomyopathy",
                    "peripartum cardiomyopathy"],
        "workflows": [CardioWorkflowType.HEART_FAILURE, CardioWorkflowType.CARDIAC_MRI],
        "collections": ["cardio_heart_failure", "cardio_imaging"],
    },

    # ── Valvular Disease ──────────────────────────────────────────────
    "aortic stenosis": {
        "aliases": ["as", "aortic valve stenosis", "avs", "severe as",
                    "calcific aortic stenosis"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging"],
    },
    "aortic regurgitation": {
        "aliases": ["ar", "aortic insufficiency", "ai", "aortic valve regurgitation"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging"],
    },
    "mitral regurgitation": {
        "aliases": ["mr", "mitral insufficiency", "mi", "mitral valve regurgitation",
                    "mitral valve prolapse", "mvp"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging"],
    },
    "mitral stenosis": {
        "aliases": ["ms", "mitral valve stenosis", "rheumatic mitral stenosis"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging"],
    },
    "tricuspid regurgitation": {
        "aliases": ["tr", "tricuspid insufficiency", "tricuspid valve regurgitation"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging"],
    },
    "endocarditis": {
        "aliases": ["infective endocarditis", "ie", "bacterial endocarditis",
                    "prosthetic valve endocarditis"],
        "workflows": [CardioWorkflowType.VALVULAR_DISEASE],
        "collections": ["cardio_valvular", "cardio_imaging", "cardio_guidelines"],
        "critical": True,
    },

    # ── Arrhythmias ───────────────────────────────────────────────────
    "atrial fibrillation": {
        "aliases": ["afib", "af", "a-fib", "atrial fib", "paroxysmal af",
                    "persistent af", "permanent af"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
    },
    "atrial flutter": {
        "aliases": ["aflutter", "typical flutter", "atypical flutter"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
    },
    "ventricular tachycardia": {
        "aliases": ["vt", "v-tach", "vtach", "sustained vt", "monomorphic vt",
                    "polymorphic vt"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
        "critical": True,
    },
    "ventricular fibrillation": {
        "aliases": ["vf", "v-fib", "vfib"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
        "critical": True,
    },
    "supraventricular tachycardia": {
        "aliases": ["svt", "avnrt", "avrt", "wpw",
                    "wolff-parkinson-white", "accessory pathway"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
    },
    "bradycardia": {
        "aliases": ["sinus bradycardia", "sick sinus syndrome", "sss",
                    "heart block", "av block", "complete heart block",
                    "third degree av block", "mobitz", "wenckebach"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_devices"],
    },
    "long qt syndrome": {
        "aliases": ["lqts", "long qt", "torsades de pointes", "tdp",
                    "qt prolongation", "brugada", "brugada syndrome"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_guidelines"],
    },
    "sudden cardiac death": {
        "aliases": ["scd", "sudden cardiac arrest", "sca",
                    "aborted sudden death"],
        "workflows": [CardioWorkflowType.ARRHYTHMIA],
        "collections": ["cardio_electrophysiology", "cardio_devices"],
        "critical": True,
    },

    # ── Vascular / Emergencies ────────────────────────────────────────
    "aortic dissection": {
        "aliases": ["dissection", "stanford type a", "stanford type b",
                    "debakey", "acute aortic syndrome", "intramural hematoma",
                    "penetrating aortic ulcer"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
        "collections": ["cardio_imaging", "cardio_interventional"],
        "critical": True,
    },
    "pulmonary embolism": {
        "aliases": ["pe", "pulmonary embolus", "massive pe", "submassive pe"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
        "collections": ["cardio_imaging", "cardio_guidelines"],
        "critical": True,
    },
    "cardiac tamponade": {
        "aliases": ["tamponade", "pericardial tamponade", "pericardial effusion"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_imaging", "cardio_hemodynamics"],
        "critical": True,
    },
    "pericarditis": {
        "aliases": ["acute pericarditis", "constrictive pericarditis",
                    "pericardial disease"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_imaging", "cardio_guidelines"],
    },
    "pulmonary hypertension": {
        "aliases": ["pah", "ph", "pulmonary arterial hypertension",
                    "group 1 ph", "group 2 ph"],
        "workflows": [CardioWorkflowType.HEART_FAILURE],
        "collections": ["cardio_hemodynamics", "cardio_guidelines"],
    },

    # ── Prevention / Metabolic ────────────────────────────────────────
    "hypertension": {
        "aliases": ["htn", "high blood pressure", "resistant hypertension",
                    "hypertensive crisis", "hypertensive emergency"],
        "workflows": [CardioWorkflowType.PREVENTIVE_RISK],
        "collections": ["cardio_prevention", "cardio_guidelines"],
    },
    "dyslipidemia": {
        "aliases": ["hyperlipidemia", "hypercholesterolemia",
                    "high cholesterol", "elevated ldl",
                    "familial hypercholesterolemia", "fh"],
        "workflows": [CardioWorkflowType.PREVENTIVE_RISK],
        "collections": ["cardio_prevention", "cardio_guidelines"],
    },

    # ── Cardio-Oncology ──────────────────────────────────────────────
    "cardiotoxicity": {
        "aliases": ["cancer therapy cardiotoxicity", "chemotherapy cardiotoxicity",
                    "anthracycline cardiotoxicity", "trastuzumab cardiotoxicity",
                    "immune checkpoint inhibitor myocarditis",
                    "radiation heart disease"],
        "workflows": [CardioWorkflowType.CARDIO_ONCOLOGY],
        "collections": ["cardio_oncology", "cardio_imaging"],
    },
}


CARDIO_BIOMARKERS: Dict[str, Dict[str, str]] = {
    "troponin": {
        "full_name": "High-sensitivity cardiac troponin (hs-cTnI / hs-cTnT)",
        "normal_range": "< 14 ng/L (hs-cTnI), < 22 ng/L (hs-cTnT, male) / < 14 ng/L (female)",
        "significance": "Myocardial injury; ACS rule-in/rule-out with serial measurements",
        "workflows": ["cad_assessment", "heart_failure"],
    },
    "bnp": {
        "full_name": "B-type natriuretic peptide",
        "normal_range": "< 100 pg/mL (HF unlikely), > 400 pg/mL (HF likely)",
        "significance": "Heart failure diagnosis and prognosis; volume status marker",
        "workflows": ["heart_failure"],
    },
    "nt-probnp": {
        "full_name": "N-terminal pro-B-type natriuretic peptide",
        "normal_range": "< 125 pg/mL (< 75 yrs), < 450 pg/mL (> 75 yrs)",
        "significance": "Heart failure diagnosis, prognosis, and treatment monitoring",
        "workflows": ["heart_failure"],
    },
    "crp": {
        "full_name": "High-sensitivity C-reactive protein (hs-CRP)",
        "normal_range": "< 2.0 mg/L (low risk), > 3.0 mg/L (high risk)",
        "significance": "Cardiovascular inflammation; residual risk assessment",
        "workflows": ["preventive_risk"],
    },
    "d-dimer": {
        "full_name": "D-dimer",
        "normal_range": "< 500 ng/mL FEU (age-adjusted: age x 10 for > 50 yrs)",
        "significance": "VTE exclusion; pulmonary embolism rule-out",
        "workflows": ["cad_assessment"],
    },
    "ldl-c": {
        "full_name": "Low-density lipoprotein cholesterol",
        "normal_range": "< 70 mg/dL (very high risk), < 100 mg/dL (high risk)",
        "significance": "Primary target for ASCVD risk reduction; statin therapy guide",
        "workflows": ["preventive_risk"],
    },
    "lp(a)": {
        "full_name": "Lipoprotein(a)",
        "normal_range": "< 50 nmol/L (< 30 mg/dL)",
        "significance": "Independent ASCVD risk factor; risk-enhancing factor per ACC/AHA",
        "workflows": ["preventive_risk"],
    },
    "hemoglobin a1c": {
        "full_name": "Glycated hemoglobin",
        "normal_range": "< 5.7% (normal), 5.7-6.4% (prediabetes), >= 6.5% (diabetes)",
        "significance": "Diabetes diagnosis; cardiovascular risk modifier",
        "workflows": ["preventive_risk"],
    },
    "creatinine": {
        "full_name": "Serum creatinine / eGFR",
        "normal_range": "0.7-1.3 mg/dL (male), 0.6-1.1 mg/dL (female); eGFR > 60 mL/min/1.73m2",
        "significance": "Renal function; cardiorenal syndrome; contrast nephropathy risk",
        "workflows": ["heart_failure", "cad_assessment"],
    },
    "potassium": {
        "full_name": "Serum potassium",
        "normal_range": "3.5-5.0 mEq/L",
        "significance": "Arrhythmia risk; MRA/ACEi/ARB monitoring; digoxin toxicity",
        "workflows": ["heart_failure", "arrhythmia"],
    },
}


CARDIO_DRUGS: Dict[str, Dict[str, object]] = {
    # ── HF GDMT ──────────────────────────────────────────────────────
    "sacubitril/valsartan": {
        "aliases": ["entresto", "arni", "sacubitril-valsartan"],
        "class": "ARNI",
        "indications": ["HFrEF", "HFmrEF"],
        "workflows": ["heart_failure"],
    },
    "carvedilol": {
        "aliases": ["coreg"],
        "class": "Beta-blocker",
        "indications": ["HFrEF", "hypertension"],
        "workflows": ["heart_failure"],
    },
    "metoprolol succinate": {
        "aliases": ["toprol xl", "metoprolol er"],
        "class": "Beta-blocker",
        "indications": ["HFrEF", "hypertension", "rate control"],
        "workflows": ["heart_failure", "arrhythmia"],
    },
    "bisoprolol": {
        "aliases": ["zebeta"],
        "class": "Beta-blocker",
        "indications": ["HFrEF"],
        "workflows": ["heart_failure"],
    },
    "spironolactone": {
        "aliases": ["aldactone"],
        "class": "MRA",
        "indications": ["HFrEF", "resistant hypertension"],
        "workflows": ["heart_failure"],
    },
    "eplerenone": {
        "aliases": ["inspra"],
        "class": "MRA",
        "indications": ["HFrEF", "post-MI"],
        "workflows": ["heart_failure"],
    },
    "dapagliflozin": {
        "aliases": ["farxiga"],
        "class": "SGLT2i",
        "indications": ["HFrEF", "HFpEF", "diabetes"],
        "workflows": ["heart_failure"],
    },
    "empagliflozin": {
        "aliases": ["jardiance"],
        "class": "SGLT2i",
        "indications": ["HFrEF", "HFpEF", "diabetes"],
        "workflows": ["heart_failure"],
    },

    # ── Antiplatelet / Anticoagulant ──────────────────────────────────
    "aspirin": {
        "aliases": ["asa", "acetylsalicylic acid"],
        "class": "Antiplatelet",
        "indications": ["ACS", "post-PCI", "secondary prevention"],
        "workflows": ["cad_assessment", "preventive_risk"],
    },
    "clopidogrel": {
        "aliases": ["plavix"],
        "class": "P2Y12 inhibitor",
        "indications": ["ACS", "post-PCI", "stroke prevention"],
        "workflows": ["cad_assessment"],
        "pgx_relevant": True,
        "pgx_gene": "CYP2C19",
    },
    "ticagrelor": {
        "aliases": ["brilinta"],
        "class": "P2Y12 inhibitor",
        "indications": ["ACS", "post-PCI"],
        "workflows": ["cad_assessment"],
    },
    "prasugrel": {
        "aliases": ["effient"],
        "class": "P2Y12 inhibitor",
        "indications": ["ACS with PCI"],
        "workflows": ["cad_assessment"],
    },
    "apixaban": {
        "aliases": ["eliquis"],
        "class": "DOAC",
        "indications": ["AF anticoagulation", "VTE treatment"],
        "workflows": ["arrhythmia"],
    },
    "rivaroxaban": {
        "aliases": ["xarelto"],
        "class": "DOAC",
        "indications": ["AF anticoagulation", "VTE treatment", "CAD"],
        "workflows": ["arrhythmia", "cad_assessment"],
    },
    "warfarin": {
        "aliases": ["coumadin"],
        "class": "Vitamin K antagonist",
        "indications": ["AF", "mechanical valve", "VTE"],
        "workflows": ["arrhythmia", "valvular_disease"],
        "pgx_relevant": True,
        "pgx_gene": "CYP2C9/VKORC1",
    },

    # ── Statins ──────────────────────────────────────────────────────
    "atorvastatin": {
        "aliases": ["lipitor"],
        "class": "Statin",
        "indications": ["Dyslipidemia", "ASCVD prevention"],
        "workflows": ["preventive_risk"],
    },
    "rosuvastatin": {
        "aliases": ["crestor"],
        "class": "Statin",
        "indications": ["Dyslipidemia", "ASCVD prevention"],
        "workflows": ["preventive_risk"],
    },

    # ── Antiarrhythmics ──────────────────────────────────────────────
    "amiodarone": {
        "aliases": ["cordarone", "pacerone"],
        "class": "Class III antiarrhythmic",
        "indications": ["AF", "VT", "cardiac arrest"],
        "workflows": ["arrhythmia"],
    },
    "sotalol": {
        "aliases": ["betapace"],
        "class": "Class III antiarrhythmic / Beta-blocker",
        "indications": ["AF", "VT"],
        "workflows": ["arrhythmia"],
    },
    "flecainide": {
        "aliases": ["tambocor"],
        "class": "Class IC antiarrhythmic",
        "indications": ["SVT", "AF (structurally normal heart)"],
        "workflows": ["arrhythmia"],
    },
    "digoxin": {
        "aliases": ["lanoxin"],
        "class": "Cardiac glycoside",
        "indications": ["HFrEF", "AF rate control"],
        "workflows": ["heart_failure", "arrhythmia"],
    },

    # ── Vasodilators / Others ────────────────────────────────────────
    "nitroglycerin": {
        "aliases": ["ntg", "glyceryl trinitrate", "nitro"],
        "class": "Nitrate",
        "indications": ["Angina", "ACS", "acute HF"],
        "workflows": ["cad_assessment", "heart_failure"],
    },
    "hydralazine/isosorbide dinitrate": {
        "aliases": ["bidil", "h-isdn"],
        "class": "Vasodilator combination",
        "indications": ["HFrEF (ACEi/ARB intolerant)", "HFrEF in Black patients"],
        "workflows": ["heart_failure"],
    },
    "ivabradine": {
        "aliases": ["corlanor"],
        "class": "If channel inhibitor",
        "indications": ["HFrEF with HR >= 70 bpm on max beta-blocker"],
        "workflows": ["heart_failure"],
    },
}


CARDIO_GENES: Dict[str, Dict[str, object]] = {
    "MYH7": {
        "full_name": "Myosin heavy chain 7",
        "conditions": ["Hypertrophic cardiomyopathy", "Dilated cardiomyopathy"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Family history of HCM, unexplained LVH in young patient",
    },
    "MYBPC3": {
        "full_name": "Myosin-binding protein C3",
        "conditions": ["Hypertrophic cardiomyopathy"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Most common HCM gene; family screening",
    },
    "TNNT2": {
        "full_name": "Troponin T type 2",
        "conditions": ["Hypertrophic cardiomyopathy", "Dilated cardiomyopathy"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "HCM with high SCD risk despite mild hypertrophy",
    },
    "LMNA": {
        "full_name": "Lamin A/C",
        "conditions": ["Dilated cardiomyopathy", "Conduction disease", "Emery-Dreifuss MD"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "DCM with conduction disease; high SCD risk",
    },
    "TTN": {
        "full_name": "Titin",
        "conditions": ["Dilated cardiomyopathy"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Most common DCM gene (~25% of familial DCM)",
    },
    "SCN5A": {
        "full_name": "Sodium channel Nav1.5",
        "conditions": ["Brugada syndrome", "Long QT syndrome type 3",
                       "Conduction disease", "Sick sinus syndrome"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Brugada pattern on ECG, unexplained syncope in young",
    },
    "KCNQ1": {
        "full_name": "Potassium channel Kv7.1",
        "conditions": ["Long QT syndrome type 1", "Short QT syndrome"],
        "inheritance": "Autosomal dominant (Romano-Ward) / Autosomal recessive (Jervell-Lange-Nielsen)",
        "testing_indication": "QTc prolongation, exercise-triggered syncope",
    },
    "KCNH2": {
        "full_name": "Potassium channel hERG (Kv11.1)",
        "conditions": ["Long QT syndrome type 2"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "QTc prolongation, auditory-triggered arrhythmia",
    },
    "RYR2": {
        "full_name": "Ryanodine receptor 2",
        "conditions": ["Catecholaminergic polymorphic VT (CPVT)"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Exercise/stress-induced bidirectional VT, SCA in young",
    },
    "PKP2": {
        "full_name": "Plakophilin-2",
        "conditions": ["Arrhythmogenic right ventricular cardiomyopathy"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "RV dilation/dysfunction, epsilon waves, T-wave inversion V1-V3",
    },
    "DSP": {
        "full_name": "Desmoplakin",
        "conditions": ["ARVC", "Dilated cardiomyopathy", "Carvajal syndrome"],
        "inheritance": "Autosomal dominant/recessive",
        "testing_indication": "LV-dominant arrhythmogenic cardiomyopathy, LGE on CMR",
    },
    "PCSK9": {
        "full_name": "Proprotein convertase subtilisin/kexin type 9",
        "conditions": ["Familial hypercholesterolemia"],
        "inheritance": "Autosomal dominant (gain-of-function)",
        "testing_indication": "Severely elevated LDL-C, family history of premature ASCVD",
    },
    "LDLR": {
        "full_name": "Low-density lipoprotein receptor",
        "conditions": ["Familial hypercholesterolemia"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "LDL-C > 190 mg/dL, tendon xanthomas, premature ASCVD",
    },
    "APOB": {
        "full_name": "Apolipoprotein B",
        "conditions": ["Familial hypercholesterolemia (FDB)",
                       "Familial hypobetalipoproteinemia"],
        "inheritance": "Autosomal dominant",
        "testing_indication": "Elevated LDL-C not explained by LDLR mutations",
    },
    "CYP2C19": {
        "full_name": "Cytochrome P450 2C19",
        "conditions": ["Clopidogrel resistance"],
        "inheritance": "Autosomal codominant",
        "testing_indication": "Pre-PCI clopidogrel prescribing; poor metabolizer screening",
        "pgx_relevant": True,
    },
    "CYP2C9": {
        "full_name": "Cytochrome P450 2C9",
        "conditions": ["Warfarin sensitivity"],
        "inheritance": "Autosomal codominant",
        "testing_indication": "Warfarin dose guidance with VKORC1",
        "pgx_relevant": True,
    },
    "VKORC1": {
        "full_name": "Vitamin K epoxide reductase complex subunit 1",
        "conditions": ["Warfarin sensitivity"],
        "inheritance": "Autosomal codominant",
        "testing_indication": "Warfarin dose guidance with CYP2C9",
        "pgx_relevant": True,
    },
}


CARDIO_IMAGING_MODALITIES: Dict[str, Dict[str, object]] = {
    "echocardiography": {
        "aliases": ["echo", "tte", "transthoracic echo", "tee",
                    "transesophageal echo", "transthoracic echocardiogram",
                    "strain imaging", "gls", "speckle tracking"],
        "workflows": [CardioWorkflowType.HEART_FAILURE, CardioWorkflowType.VALVULAR_DISEASE],
    },
    "cardiac ct": {
        "aliases": ["ccta", "coronary ct angiography", "coronary cta",
                    "calcium score", "cac score", "ct coronary"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT, CardioWorkflowType.PREVENTIVE_RISK],
    },
    "cardiac mri": {
        "aliases": ["cmr", "cardiac magnetic resonance", "cardiac mr",
                    "lge", "late gadolinium enhancement", "t1 mapping",
                    "t2 mapping", "myocardial perfusion mri"],
        "workflows": [CardioWorkflowType.CARDIAC_MRI],
    },
    "nuclear perfusion": {
        "aliases": ["myocardial perfusion imaging", "mpi", "spect",
                    "nuclear stress test", "sestamibi", "thallium"],
        "workflows": [CardioWorkflowType.STRESS_TEST, CardioWorkflowType.CAD_ASSESSMENT],
    },
    "pet": {
        "aliases": ["cardiac pet", "rubidium pet", "ammonia pet",
                    "fdg pet", "pet-ct cardiac", "myocardial viability pet"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT, CardioWorkflowType.CARDIAC_MRI],
    },
    "angiography": {
        "aliases": ["cardiac catheterization", "cath", "coronary angiogram",
                    "cardiac cath", "left heart cath", "right heart cath",
                    "rhc", "lhc", "diagnostic angiography"],
        "workflows": [CardioWorkflowType.CAD_ASSESSMENT],
    },
    "stress test": {
        "aliases": ["exercise stress test", "treadmill test", "bruce protocol",
                    "stress echo", "dobutamine stress echo", "dse",
                    "pharmacologic stress test", "adenosine stress",
                    "regadenoson stress", "exercise ecg"],
        "workflows": [CardioWorkflowType.STRESS_TEST],
    },
}


CARDIO_GUIDELINES: Dict[str, Dict[str, str]] = {
    "acc_aha_hf_2022": {
        "title": "2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure",
        "citation": "Heidenreich PA, et al. Circulation. 2022;145(18):e895-e1032.",
        "pmid": "35363499",
        "year": "2022",
    },
    "acc_aha_chest_pain_2021": {
        "title": "2021 AHA/ACC Guideline for the Evaluation and Diagnosis of Chest Pain",
        "citation": "Gulati M, et al. Circulation. 2021;144(22):e368-e454.",
        "pmid": "34709928",
        "year": "2021",
    },
    "acc_aha_vhd_2020": {
        "title": "2020 ACC/AHA Guideline for the Management of Patients With Valvular Heart Disease",
        "citation": "Otto CM, et al. Circulation. 2021;143(5):e72-e227.",
        "pmid": "33332150",
        "year": "2020",
    },
    "acc_aha_af_2023": {
        "title": "2023 ACC/AHA/ACCP/HRS Guideline for Diagnosis and Management of Atrial Fibrillation",
        "citation": "Joglar JA, et al. Circulation. 2024;149(1):e1-e156.",
        "pmid": "38033089",
        "year": "2023",
    },
    "esc_af_2024": {
        "title": "2024 ESC Guidelines for the Management of Atrial Fibrillation",
        "citation": "Van Gelder IC, et al. Eur Heart J. 2024;45(36):3314-3414.",
        "pmid": "39210723",
        "year": "2024",
    },
    "acc_aha_cholesterol_2018": {
        "title": "2018 AHA/ACC Guideline on the Management of Blood Cholesterol",
        "citation": "Grundy SM, et al. Circulation. 2019;139(25):e1082-e1143.",
        "pmid": "30586774",
        "year": "2018",
    },
    "acc_aha_prevention_2019": {
        "title": "2019 ACC/AHA Guideline on the Primary Prevention of Cardiovascular Disease",
        "citation": "Arnett DK, et al. Circulation. 2019;140(11):e596-e646.",
        "pmid": "30879355",
        "year": "2019",
    },
    "esc_cardio_oncology_2022": {
        "title": "2022 ESC Guidelines on Cardio-Oncology",
        "citation": "Lyon AR, et al. Eur Heart J. 2022;43(41):4229-4361.",
        "pmid": "36017568",
        "year": "2022",
    },
    "acc_aha_hcm_2024": {
        "title": "2024 AHA/ACC Guideline for the Management of Hypertrophic Cardiomyopathy",
        "citation": "Ommen SR, et al. Circulation. 2024;149(23):e1239-e1311.",
        "pmid": "38718139",
        "year": "2024",
    },
    "acc_aha_hypertension_2017": {
        "title": "2017 ACC/AHA Guideline for the Prevention, Detection, Evaluation, and Management of High Blood Pressure in Adults",
        "citation": "Whelton PK, et al. J Am Coll Cardiol. 2018;71(19):e127-e248.",
        "pmid": "29146535",
        "year": "2017",
    },
}


# =====================================================================
# SEARCH PLAN DATACLASS
# =====================================================================

@dataclass
class SearchPlan:
    """Agent's plan for answering a cardiovascular medicine question.

    The search plan captures all entities detected in the user's question
    and the strategy the agent will use to retrieve evidence from the
    12 cardiology-specific Milvus collections.
    """
    question: str
    conditions: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    imaging_modalities: List[str] = field(default_factory=list)
    genes: List[str] = field(default_factory=list)
    relevant_workflows: List[CardioWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad, targeted, comparative, clinical
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# =====================================================================
# CARDIOLOGY INTELLIGENCE AGENT
# =====================================================================

class CardioIntelligenceAgent:
    """Autonomous Cardiology Intelligence Agent.

    Wraps the multi-collection CardioRAGEngine with planning and reasoning
    capabilities. Designed to answer complex cross-functional questions
    about cardiovascular disease, treatment, imaging, and risk assessment.

    Example queries this agent handles:
    - "What is the optimal GDMT regimen for a patient with HFrEF and LVEF 30%?"
    - "Interpret a cardiac MRI showing mid-wall LGE in a patient with DCM"
    - "Should I start anticoagulation for a patient with AF and CHA2DS2-VASc of 2?"
    - "Compare TAVR vs SAVR for severe aortic stenosis in an 80-year-old"
    - "Is clopidogrel safe for a CYP2C19 poor metabolizer post-PCI?"
    - "What is the work-up for a young patient with aborted sudden cardiac death?"
    - "Evaluate cardiotoxicity risk for a breast cancer patient starting doxorubicin"
    - "What stress test modality is best for a patient with LBBB?"

    Usage:
        agent = CardioIntelligenceAgent(rag_engine)
        plan = agent.search_plan("What is GDMT for HFrEF?")
        response = agent.run("What is GDMT for HFrEF?")
    """

    def __init__(self, rag_engine):
        """Initialize agent with a configured RAG engine.

        Args:
            rag_engine: CardioRAGEngine instance with Milvus collections connected.
        """
        self.rag = rag_engine
        self.knowledge = {
            "conditions": CARDIO_CONDITIONS,
            "biomarkers": CARDIO_BIOMARKERS,
            "drugs": CARDIO_DRUGS,
            "genes": CARDIO_GENES,
            "imaging": CARDIO_IMAGING_MODALITIES,
            "guidelines": CARDIO_GUIDELINES,
        }

    # ── Public API ──────────────────────────────────────────────────

    def run(self, question: str, **kwargs) -> WorkflowResult:
        """Execute the full agent pipeline: plan -> search -> evaluate -> synthesize.

        Args:
            question: Natural language question about cardiovascular medicine.
            **kwargs: Additional query parameters (patient_context, workflow_override,
                      top_k, collection_filter).

        Returns:
            WorkflowResult with findings, recommendations, and metadata.
        """
        # Phase 1: Plan
        plan = self.search_plan(question)

        # Phase 2: Determine workflow (allow override)
        workflow_override = kwargs.get("workflow_override")
        workflow = workflow_override or (
            plan.relevant_workflows[0] if plan.relevant_workflows else None
        )

        # Phase 3: Search via RAG engine
        patient_context = kwargs.get("patient_context")
        top_k = kwargs.get("top_k", 5)

        response = self.rag.query(
            question=question,
            workflow=workflow,
            top_k=top_k,
            patient_context=patient_context,
        )

        # Phase 4: Evaluate and potentially expand
        if hasattr(response, "results") and response.results is not None:
            quality = self._evaluate_evidence_quality(response.results)
            if quality == "insufficient" and plan.sub_questions:
                for sub_q in plan.sub_questions[:2]:
                    sub_response = self.rag.search(sub_q, top_k=top_k)
                    if hasattr(response, "results") and sub_response:
                        response.results.extend(sub_response)

        return response

    def search_plan(self, question: str) -> SearchPlan:
        """Analyze a question and create an optimised search plan.

        Detects cardiovascular conditions, drugs, genes, and imaging modalities
        in the question text. Determines relevant clinical workflows, chooses
        a search strategy, and generates sub-questions for comprehensive
        retrieval across collections.

        Args:
            question: The user's natural language question.

        Returns:
            SearchPlan with all detected entities and retrieval strategy.
        """
        plan = SearchPlan(question=question)

        # Step 1: Detect clinical entities
        plan.conditions = self._detect_conditions(question)
        plan.drugs = self._detect_drugs(question)
        plan.genes = self._detect_genes(question)
        plan.imaging_modalities = self._detect_imaging(question)

        # Step 2: Determine relevant workflows
        plan.relevant_workflows = self._determine_workflows(
            plan.conditions, plan.drugs, plan.genes,
            plan.imaging_modalities, question,
        )

        # Step 3: Choose search strategy
        plan.search_strategy = self._choose_strategy(
            question, plan.conditions, plan.drugs,
        )

        # Step 4: Generate sub-questions
        plan.sub_questions = self._generate_sub_questions(plan)

        # Step 5: Compile identified topics
        plan.identified_topics = (
            plan.conditions + plan.drugs + plan.genes + plan.imaging_modalities
        )

        return plan

    # ── Entity Detection ────────────────────────────────────────────

    def _detect_conditions(self, text: str) -> List[str]:
        """Detect cardiovascular conditions mentioned in text.

        Performs case-insensitive matching against the CARDIO_CONDITIONS
        knowledge dictionary, checking both canonical names and aliases.

        Args:
            text: Input text to scan for condition mentions.

        Returns:
            List of canonical condition names found in the text.
        """
        detected: List[str] = []
        text_lower = text.lower()

        for condition, info in CARDIO_CONDITIONS.items():
            # Check canonical name
            if condition in text_lower:
                if condition not in detected:
                    detected.append(condition)
                continue

            # Check aliases
            aliases = info.get("aliases", [])
            for alias in aliases:
                # For very short aliases (2-3 chars), require word boundary
                if len(alias) <= 3:
                    # Check as isolated word
                    import re
                    pattern = r'\b' + re.escape(alias) + r'\b'
                    if re.search(pattern, text_lower):
                        if condition not in detected:
                            detected.append(condition)
                        break
                else:
                    if alias.lower() in text_lower:
                        if condition not in detected:
                            detected.append(condition)
                        break

        return detected

    def _detect_drugs(self, text: str) -> List[str]:
        """Detect cardiac drugs mentioned in text.

        Performs case-insensitive matching against the CARDIO_DRUGS
        knowledge dictionary, checking both generic names and brand aliases.

        Args:
            text: Input text to scan for drug mentions.

        Returns:
            List of canonical drug names found in the text.
        """
        detected: List[str] = []
        text_lower = text.lower()

        for drug, info in CARDIO_DRUGS.items():
            if drug.lower() in text_lower:
                if drug not in detected:
                    detected.append(drug)
                continue

            aliases = info.get("aliases", [])
            for alias in aliases:
                if alias.lower() in text_lower:
                    if drug not in detected:
                        detected.append(drug)
                    break

        # Also check drug class names
        drug_class_keywords = {
            "arni": "sacubitril/valsartan",
            "beta-blocker": None,  # multiple drugs
            "beta blocker": None,
            "mra": None,
            "sglt2i": None,
            "sglt2 inhibitor": None,
            "statin": None,
            "doac": None,
            "antiplatelet": None,
            "antiarrhythmic": None,
            "nitrate": "nitroglycerin",
            "ace inhibitor": None,
            "acei": None,
            "arb": None,
        }

        for keyword, mapped_drug in drug_class_keywords.items():
            if keyword in text_lower:
                if mapped_drug and mapped_drug not in detected:
                    detected.append(mapped_drug)
                if keyword not in [d.lower() for d in detected]:
                    # Add the class as a topic for search expansion
                    pass

        return detected

    def _detect_genes(self, text: str) -> List[str]:
        """Detect cardiovascular genes mentioned in text.

        Performs case-insensitive matching against the CARDIO_GENES
        knowledge dictionary. Gene symbols are typically uppercase
        (MYH7, SCN5A, etc.) but the search is case-insensitive.

        Args:
            text: Input text to scan for gene mentions.

        Returns:
            List of gene symbols found in the text.
        """
        detected: List[str] = []
        text_upper = text.upper()

        for gene in CARDIO_GENES:
            if gene.upper() in text_upper:
                if gene not in detected:
                    detected.append(gene)

        return detected

    def _detect_imaging(self, text: str) -> List[str]:
        """Detect imaging modalities mentioned in text.

        Performs case-insensitive matching against the CARDIO_IMAGING_MODALITIES
        knowledge dictionary, checking both canonical names and aliases.

        Args:
            text: Input text to scan for imaging modality mentions.

        Returns:
            List of canonical imaging modality names found in the text.
        """
        detected: List[str] = []
        text_lower = text.lower()

        for modality, info in CARDIO_IMAGING_MODALITIES.items():
            if modality in text_lower:
                if modality not in detected:
                    detected.append(modality)
                continue

            aliases = info.get("aliases", [])
            for alias in aliases:
                if len(alias) <= 3:
                    import re
                    pattern = r'\b' + re.escape(alias) + r'\b'
                    if re.search(pattern, text_lower):
                        if modality not in detected:
                            detected.append(modality)
                        break
                else:
                    if alias.lower() in text_lower:
                        if modality not in detected:
                            detected.append(modality)
                        break

        return detected

    # ── Workflow Determination ──────────────────────────────────────

    def _determine_workflows(
        self,
        conditions: List[str],
        drugs: List[str],
        genes: List[str],
        imaging: List[str],
        text: str,
    ) -> List[CardioWorkflowType]:
        """Determine relevant clinical workflows based on detected entities.

        Uses a combination of entity-workflow mappings from the knowledge
        dictionaries and keyword-based heuristics to identify which of
        the 8 clinical workflows are most relevant to the query.

        Args:
            conditions: Detected cardiovascular conditions.
            drugs: Detected cardiac drugs.
            genes: Detected cardiovascular genes.
            imaging: Detected imaging modalities.
            text: Original query text for keyword fallback.

        Returns:
            Ordered list of CardioWorkflowType values, most relevant first.
        """
        workflow_scores: Dict[CardioWorkflowType, float] = {}

        # Score from conditions
        for condition in conditions:
            info = CARDIO_CONDITIONS.get(condition, {})
            for wf in info.get("workflows", []):
                workflow_scores[wf] = workflow_scores.get(wf, 0) + 2.0
                # Boost critical conditions
                if info.get("critical"):
                    workflow_scores[wf] += 1.0

        # Score from drugs
        for drug in drugs:
            info = CARDIO_DRUGS.get(drug, {})
            for wf_name in info.get("workflows", []):
                try:
                    wf = CardioWorkflowType(wf_name)
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.5
                except ValueError:
                    pass

        # Score from genes
        for gene in genes:
            info = CARDIO_GENES.get(gene, {})
            gene_conditions = info.get("conditions", [])
            for cond_name in gene_conditions:
                # Map gene conditions back to workflows
                for cond_key, cond_info in CARDIO_CONDITIONS.items():
                    cond_aliases = [cond_key] + cond_info.get("aliases", [])
                    if any(cond_name.lower() in alias.lower() for alias in cond_aliases):
                        for wf in cond_info.get("workflows", []):
                            workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.0
                        break

        # Score from imaging modalities
        for modality in imaging:
            info = CARDIO_IMAGING_MODALITIES.get(modality, {})
            for wf in info.get("workflows", []):
                workflow_scores[wf] = workflow_scores.get(wf, 0) + 1.5

        # Keyword-based fallback scoring
        text_upper = text.upper()
        keyword_workflow_map = {
            CardioWorkflowType.CAD_ASSESSMENT: [
                "CHEST PAIN", "ANGINA", "CORONARY", "STENT", "PCI", "CABG",
                "ISCHEMIA", "INFARCTION", "ACS", "TROPONIN", "CAD",
            ],
            CardioWorkflowType.HEART_FAILURE: [
                "HEART FAILURE", "HF", "EJECTION FRACTION", "LVEF", "GDMT",
                "NYHA", "CONGESTION", "DIURETIC", "CARDIOMYOPATHY",
                "DECOMPENSATED", "BNP", "NT-PROBNP",
            ],
            CardioWorkflowType.VALVULAR_DISEASE: [
                "VALVE", "VALVULAR", "STENOSIS", "REGURGITATION", "TAVR",
                "SAVR", "MITRACLIP", "PROLAPSE", "ENDOCARDITIS",
            ],
            CardioWorkflowType.ARRHYTHMIA: [
                "ARRHYTHMIA", "RHYTHM", "ATRIAL FIBRILLATION", "AFIB",
                "FLUTTER", "TACHYCARDIA", "BRADYCARDIA", "ABLATION",
                "PACEMAKER", "ICD", "DEFIBRILLATOR", "QT", "ECG", "EKG",
            ],
            CardioWorkflowType.CARDIAC_MRI: [
                "CMR", "CARDIAC MRI", "CARDIAC MR", "LGE",
                "LATE GADOLINIUM", "T1 MAPPING", "T2 MAPPING",
                "MYOCARDIAL TISSUE", "FIBROSIS",
            ],
            CardioWorkflowType.STRESS_TEST: [
                "STRESS TEST", "EXERCISE TEST", "TREADMILL", "BRUCE",
                "DOBUTAMINE", "STRESS ECHO", "NUCLEAR STRESS",
                "PHARMACOLOGIC STRESS", "REGADENOSON", "ADENOSINE",
            ],
            CardioWorkflowType.PREVENTIVE_RISK: [
                "PREVENTION", "RISK FACTOR", "ASCVD", "CHOLESTEROL",
                "LDL", "STATIN", "BLOOD PRESSURE", "HYPERTENSION",
                "SCREENING", "RISK SCORE", "RISK CALCULATOR",
                "LIFESTYLE", "SMOKING CESSATION",
            ],
            CardioWorkflowType.CARDIO_ONCOLOGY: [
                "CARDIOTOXICITY", "CHEMOTHERAPY", "ANTHRACYCLINE",
                "TRASTUZUMAB", "IMMUNOTHERAPY", "CHECKPOINT INHIBITOR",
                "RADIATION HEART", "CANCER", "ONCOLOGY",
            ],
        }

        for wf, keywords in keyword_workflow_map.items():
            for kw in keywords:
                if kw in text_upper:
                    workflow_scores[wf] = workflow_scores.get(wf, 0) + 0.5

        # Sort by score descending and return
        if not workflow_scores:
            return [CardioWorkflowType.CAD_ASSESSMENT]

        sorted_workflows = sorted(
            workflow_scores.items(), key=lambda x: x[1], reverse=True,
        )

        return [wf for wf, _score in sorted_workflows]

    # ── Strategy Selection ─────────────────────────────────────────

    def _choose_strategy(
        self,
        text: str,
        conditions: List[str],
        drugs: List[str],
    ) -> str:
        """Choose search strategy: broad, targeted, comparative, or clinical.

        Args:
            text: Original query text.
            conditions: Detected conditions.
            drugs: Detected drugs.

        Returns:
            Strategy name string.
        """
        text_upper = text.upper()

        # Comparative queries
        if ("COMPARE" in text_upper or " VS " in text_upper
                or "VERSUS" in text_upper or "DIFFERENCE BETWEEN" in text_upper
                or "COMPARING" in text_upper):
            return "comparative"

        # Clinical case / patient-specific queries
        clinical_keywords = [
            "PATIENT", "YEAR OLD", "YO ", "Y/O", "PRESENTS WITH",
            "ADMITTED", "HISTORY OF", "PMH", "LABS SHOW",
            "ECG SHOWS", "ECHO SHOWS", "CT SHOWS", "MRI SHOWS",
        ]
        if any(kw in text_upper for kw in clinical_keywords):
            return "clinical"

        # Targeted: specific condition + drug or single focused entity
        if (len(conditions) == 1 and len(drugs) <= 1) or (
            len(conditions) <= 1 and len(drugs) == 1
        ):
            if conditions or drugs:
                return "targeted"

        # Default to broad
        return "broad"

    # ── Sub-Question Generation ────────────────────────────────────

    def _generate_sub_questions(self, plan: SearchPlan) -> List[str]:
        """Generate sub-questions for comprehensive retrieval.

        Decomposes the main question into focused sub-queries based on
        the detected entities and workflow type. This enables multi-hop
        retrieval across different aspects of the clinical question.

        Args:
            plan: SearchPlan with detected entities and workflows.

        Returns:
            List of sub-question strings (typically 2-4 questions).
        """
        sub_questions: List[str] = []
        plan.question.upper()

        condition_label = plan.conditions[0] if plan.conditions else "the condition"
        drug_label = plan.drugs[0] if plan.drugs else "the medication"
        gene_label = plan.genes[0] if plan.genes else "the gene"

        # ── Pattern 1: Heart Failure / GDMT ───────────────────────
        if any(wf == CardioWorkflowType.HEART_FAILURE
               for wf in plan.relevant_workflows):
            sub_questions = [
                f"What are the current ACC/AHA guideline recommendations for {condition_label}?",
                f"What is the optimal GDMT titration strategy for {condition_label}?",
                f"What are the key landmark trials supporting treatment of {condition_label}?",
                f"What device therapy or advanced therapies are indicated for {condition_label}?",
            ]

        # ── Pattern 2: CAD / Ischemia ─────────────────────────────
        elif any(wf == CardioWorkflowType.CAD_ASSESSMENT
                 for wf in plan.relevant_workflows):
            sub_questions = [
                f"What diagnostic approach is recommended for {condition_label}?",
                f"What are the interventional options for {condition_label}?",
                f"What is the optimal medical therapy for {condition_label}?",
                f"What risk stratification tools apply to {condition_label}?",
            ]

        # ── Pattern 3: Arrhythmia ─────────────────────────────────
        elif any(wf == CardioWorkflowType.ARRHYTHMIA
                 for wf in plan.relevant_workflows):
            sub_questions = [
                f"What are the diagnostic criteria for {condition_label}?",
                f"What pharmacologic treatment options exist for {condition_label}?",
                f"When is ablation or device therapy indicated for {condition_label}?",
                f"What is the anticoagulation strategy for {condition_label}?",
            ]

        # ── Pattern 4: Valvular Disease ───────────────────────────
        elif any(wf == CardioWorkflowType.VALVULAR_DISEASE
                 for wf in plan.relevant_workflows):
            sub_questions = [
                f"What are the severity criteria for {condition_label}?",
                f"When is intervention indicated for {condition_label}?",
                f"What are the surgical vs transcatheter options for {condition_label}?",
                f"What serial monitoring is recommended for {condition_label}?",
            ]

        # ── Pattern 5: Imaging Interpretation ─────────────────────
        elif any(wf in (CardioWorkflowType.CARDIAC_MRI, CardioWorkflowType.STRESS_TEST)
                 for wf in plan.relevant_workflows):
            modality = plan.imaging_modalities[0] if plan.imaging_modalities else "imaging"
            sub_questions = [
                f"What are the normal reference ranges for {modality} measurements?",
                f"What findings on {modality} suggest pathology?",
                f"What are the clinical implications of abnormal {modality} findings?",
                f"When should additional imaging be performed after {modality}?",
            ]

        # ── Pattern 6: Prevention / Risk ──────────────────────────
        elif any(wf == CardioWorkflowType.PREVENTIVE_RISK
                 for wf in plan.relevant_workflows):
            sub_questions = [
                f"What is the recommended risk assessment approach for {condition_label}?",
                f"What are the treatment targets for {condition_label}?",
                f"What lifestyle modifications are recommended for {condition_label}?",
                f"What pharmacotherapy is indicated for {condition_label}?",
            ]

        # ── Pattern 7: Cardio-Oncology ────────────────────────────
        elif any(wf == CardioWorkflowType.CARDIO_ONCOLOGY
                 for wf in plan.relevant_workflows):
            sub_questions = [
                f"What is the baseline cardiac assessment for {drug_label} therapy?",
                f"What are the cardiotoxicity monitoring protocols for {drug_label}?",
                f"What biomarker thresholds trigger cardiac intervention during {drug_label}?",
                f"What cardioprotective strategies are available during {drug_label}?",
            ]

        # ── Pattern 8: Comparative ────────────────────────────────
        elif plan.search_strategy == "comparative":
            entities = plan.conditions + plan.drugs
            if len(entities) >= 2:
                sub_questions = [
                    f"What is the evidence for {entities[0]} in cardiovascular disease?",
                    f"What is the evidence for {entities[1]} in cardiovascular disease?",
                    f"What head-to-head trials compare {entities[0]} and {entities[1]}?",
                ]
            else:
                sub_questions = [
                    f"What are the treatment options for {condition_label}?",
                    f"What comparative evidence exists for {condition_label} treatments?",
                ]

        # ── Pattern 9: Genomic / PGx ─────────────────────────────
        elif plan.genes:
            sub_questions = [
                f"What cardiovascular conditions are associated with {gene_label} variants?",
                f"What is the recommended genetic testing approach for {gene_label}?",
                f"What are the management implications of {gene_label} variants?",
                f"What family screening is recommended for {gene_label} carriers?",
            ]

        # ── Default ──────────────────────────────────────────────
        else:
            sub_questions = [
                f"What are the current guidelines for {condition_label}?",
                f"What is the recommended management approach for {condition_label}?",
                f"What is the prognosis and follow-up for {condition_label}?",
            ]

        return sub_questions

    # ── Evidence Evaluation ────────────────────────────────────────

    def _evaluate_evidence_quality(self, results) -> str:
        """Evaluate the quality and coverage of retrieved evidence.

        Uses collection diversity and hit count to assess whether
        the retrieved evidence is sufficient for a comprehensive answer.

        Args:
            results: List of search results from the RAG engine.

        Returns:
            "sufficient", "partial", or "insufficient".
        """
        if not results:
            return "insufficient"

        total_hits = len(results)
        collections_seen = set()

        for result in results:
            if hasattr(result, "collection"):
                collections_seen.add(result.collection)
            elif isinstance(result, dict):
                collections_seen.add(result.get("collection", "unknown"))

        num_collections = len(collections_seen)

        if num_collections >= 3 and total_hits >= 10:
            return "sufficient"
        elif num_collections >= 2 and total_hits >= 5:
            return "partial"
        else:
            return "insufficient"

    # ── Report Generation ──────────────────────────────────────────

    def generate_report(self, question: str, response) -> str:
        """Generate a structured cardiology analysis report.

        Args:
            question: Original query text.
            response: Response object from run() or rag.query().

        Returns:
            Formatted markdown report string.
        """
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        plan = self.search_plan(question)

        report_lines = [
            "# Cardiology Intelligence Report",
            f"**Query:** {question}",
            f"**Generated:** {timestamp}",
            f"**Workflows:** {', '.join(wf.value for wf in plan.relevant_workflows)}",
            f"**Strategy:** {plan.search_strategy}",
            "",
        ]

        # Detected entities
        if plan.conditions or plan.drugs or plan.genes:
            report_lines.extend([
                "---",
                "",
                "## Detected Clinical Entities",
                "",
            ])
            if plan.conditions:
                report_lines.append(
                    f"- **Conditions:** {', '.join(plan.conditions)}"
                )
            if plan.drugs:
                report_lines.append(
                    f"- **Medications:** {', '.join(plan.drugs)}"
                )
            if plan.genes:
                report_lines.append(
                    f"- **Genes:** {', '.join(plan.genes)}"
                )
            if plan.imaging_modalities:
                report_lines.append(
                    f"- **Imaging:** {', '.join(plan.imaging_modalities)}"
                )
            report_lines.append("")

        # Critical findings check
        critical_conditions = [
            c for c in plan.conditions
            if CARDIO_CONDITIONS.get(c, {}).get("critical")
        ]
        if critical_conditions:
            report_lines.extend([
                "---",
                "",
                "## [CRITICAL] Urgent Findings",
                "",
            ])
            for cond in critical_conditions:
                report_lines.append(
                    f"- **[CRITICAL]** {cond.title()} detected -- "
                    f"immediate clinical evaluation required."
                )
            report_lines.append("")

        # Analysis section
        report_lines.extend([
            "---",
            "",
            "## Analysis",
            "",
        ])

        if hasattr(response, "answer"):
            report_lines.append(response.answer)
        elif hasattr(response, "summary"):
            report_lines.append(response.summary)
        elif isinstance(response, str):
            report_lines.append(response)
        else:
            report_lines.append("No analysis generated.")

        report_lines.append("")

        # Guideline references
        if plan.relevant_workflows:
            report_lines.extend([
                "---",
                "",
                "## Applicable Guidelines",
                "",
            ])
            # Map workflows to relevant guidelines
            workflow_guideline_map = {
                CardioWorkflowType.HEART_FAILURE: ["acc_aha_hf_2022"],
                CardioWorkflowType.CAD_ASSESSMENT: ["acc_aha_chest_pain_2021"],
                CardioWorkflowType.VALVULAR_DISEASE: ["acc_aha_vhd_2020"],
                CardioWorkflowType.ARRHYTHMIA: ["acc_aha_af_2023", "esc_af_2024"],
                CardioWorkflowType.PREVENTIVE_RISK: [
                    "acc_aha_cholesterol_2018", "acc_aha_prevention_2019",
                    "acc_aha_hypertension_2017",
                ],
                CardioWorkflowType.CARDIO_ONCOLOGY: ["esc_cardio_oncology_2022"],
                CardioWorkflowType.CARDIAC_MRI: [],
                CardioWorkflowType.STRESS_TEST: ["acc_aha_chest_pain_2021"],
            }
            referenced = set()
            for wf in plan.relevant_workflows:
                for gl_key in workflow_guideline_map.get(wf, []):
                    if gl_key not in referenced:
                        referenced.add(gl_key)
                        gl = CARDIO_GUIDELINES.get(gl_key, {})
                        if gl:
                            report_lines.append(
                                f"- {gl['title']} "
                                f"([PMID:{gl['pmid']}]"
                                f"(https://pubmed.ncbi.nlm.nih.gov/{gl['pmid']}/))"
                            )
            report_lines.append("")

        # Disclaimer
        report_lines.extend([
            "---",
            "",
            "*This report is generated by the Cardiology Intelligence Agent "
            "for clinical decision support only. All recommendations require "
            "clinical correlation and physician judgment.*",
        ])

        return "\n".join(report_lines)
