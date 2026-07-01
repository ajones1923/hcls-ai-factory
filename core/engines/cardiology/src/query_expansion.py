"""Cardiology-specific query expansion for the Cardiology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Maps cardiovascular terminology to related terms, enabling comprehensive
retrieval across the 12 Milvus collections.  Each expansion map covers a
clinical domain (heart failure, coronary artery disease, arrhythmia, etc.)
and widens the semantic net so that a user searching for, e.g., "CHF" also
pulls back documents mentioning LVEF, BNP, ARNI, SGLT2i, and so on.

The QueryExpander class orchestrates alias resolution, entity detection,
comparative-query parsing, workflow-aware term boosting, and collection
weight suggestions.
"""

import re
from typing import Dict, List, Optional

from src.models import CardioWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# 1. HEART FAILURE MAP
# ═══════════════════════════════════════════════════════════════════════

HEART_FAILURE_MAP: Dict[str, List[str]] = {
    "hf": [
        "heart failure", "CHF", "congestive heart failure", "LVEF",
        "ejection fraction", "BNP", "NT-proBNP", "ARNI", "sacubitril",
        "valsartan", "entresto", "SGLT2i", "dapagliflozin", "empagliflozin",
        "diuretic", "furosemide", "bumetanide", "carvedilol", "metoprolol",
        "bisoprolol", "spironolactone", "eplerenone", "GDMT",
        "guideline-directed medical therapy", "NYHA", "Stage C", "Stage D",
        "congestion", "volume overload", "decompensation", "pulmonary edema",
        "jugular venous distension", "S3 gallop", "cardiomyopathy",
        "vericiguat", "verquvo",
    ],
    "chf": [
        "congestive heart failure", "heart failure", "HF", "LVEF",
        "ejection fraction", "BNP", "NT-proBNP", "volume overload",
        "congestion", "pulmonary edema", "decompensation", "diuretic",
        "furosemide", "ARNI", "SGLT2i", "GDMT", "NYHA",
    ],
    "heart failure": [
        "HF", "CHF", "LVEF", "ejection fraction", "BNP", "NT-proBNP",
        "ARNI", "sacubitril", "valsartan", "SGLT2i", "dapagliflozin",
        "empagliflozin", "diuretic", "furosemide", "carvedilol",
        "metoprolol", "spironolactone", "eplerenone", "GDMT", "NYHA",
        "Stage C", "congestion", "volume overload", "decompensation",
        "cardiomyopathy", "HFrEF", "HFpEF", "HFmrEF",
    ],
    "systolic dysfunction": [
        "HFrEF", "heart failure with reduced ejection fraction", "LVEF",
        "ejection fraction", "reduced EF", "EF less than 40", "cardiomyopathy",
        "dilated cardiomyopathy", "GDMT", "ARNI", "beta blocker", "MRA",
        "SGLT2i", "ICD", "CRT", "wall motion abnormality",
    ],
    "diastolic dysfunction": [
        "HFpEF", "heart failure with preserved ejection fraction",
        "E/A ratio", "E/e' ratio", "LA enlargement", "left atrial volume",
        "diastolic filling", "relaxation abnormality", "restrictive filling",
        "pseudonormal", "grade I", "grade II", "grade III", "SGLT2i",
        "diuretic", "volume management", "spironolactone",
    ],
    "hfref": [
        "heart failure with reduced ejection fraction", "systolic dysfunction",
        "LVEF", "ejection fraction less than 40", "GDMT", "ARNI",
        "sacubitril", "beta blocker", "carvedilol", "metoprolol",
        "MRA", "spironolactone", "SGLT2i", "ICD", "CRT",
        "decompensation", "BNP", "NT-proBNP",
    ],
    "hfpef": [
        "heart failure with preserved ejection fraction",
        "diastolic dysfunction", "LVEF greater than 50", "SGLT2i",
        "empagliflozin", "dapagliflozin", "diuretic", "volume management",
        "E/e' ratio", "LA enlargement", "exercise intolerance",
        "congestion", "spironolactone", "BNP", "NT-proBNP",
    ],
    "bnp": [
        "brain natriuretic peptide", "NT-proBNP", "heart failure",
        "congestion", "volume overload", "decompensation", "biomarker",
        "HF", "CHF", "LVEF", "dyspnea", "prognosis",
    ],
    "gdmt": [
        "guideline-directed medical therapy", "ARNI", "sacubitril",
        "valsartan", "beta blocker", "carvedilol", "metoprolol",
        "bisoprolol", "MRA", "spironolactone", "eplerenone", "SGLT2i",
        "dapagliflozin", "empagliflozin", "ACE inhibitor", "ARB",
        "four pillars", "heart failure", "HFrEF", "titration",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 2. CORONARY ARTERY DISEASE MAP
# ═══════════════════════════════════════════════════════════════════════

CORONARY_ARTERY_MAP: Dict[str, List[str]] = {
    "cad": [
        "coronary artery disease", "atherosclerosis", "coronary",
        "ischemia", "MI", "myocardial infarction", "ACS",
        "acute coronary syndrome", "angina", "stable angina",
        "troponin", "angiography", "PCI", "percutaneous coronary intervention",
        "CABG", "coronary artery bypass graft", "stent", "DES",
        "drug-eluting stent", "FFR", "fractional flow reserve", "iFR",
        "instantaneous wave-free ratio", "calcium score", "CTA",
        "coronary CT angiography", "plaque", "stenosis", "LAD", "RCA",
        "LCx", "left main", "multivessel", "revascularization",
    ],
    "coronary": [
        "coronary artery disease", "CAD", "atherosclerosis", "stenosis",
        "plaque", "ischemia", "angina", "troponin", "angiography", "PCI",
        "CABG", "stent", "FFR", "iFR", "CTA", "calcium score",
        "LAD", "RCA", "LCx", "revascularization",
    ],
    "mi": [
        "myocardial infarction", "heart attack", "STEMI",
        "ST-elevation myocardial infarction", "NSTEMI",
        "non-ST-elevation myocardial infarction", "ACS",
        "acute coronary syndrome", "troponin", "troponin I", "troponin T",
        "CK-MB", "angiography", "PCI", "thrombolysis", "aspirin",
        "P2Y12 inhibitor", "clopidogrel", "ticagrelor", "prasugrel",
        "heparin", "DAPT", "dual antiplatelet therapy", "beta blocker",
        "ACE inhibitor", "statin", "wall motion abnormality",
        "cardiac rehabilitation",
    ],
    "myocardial infarction": [
        "MI", "heart attack", "STEMI", "NSTEMI", "ACS",
        "acute coronary syndrome", "troponin", "angiography", "PCI",
        "thrombolysis", "CABG", "reperfusion", "ischemia", "necrosis",
        "DAPT", "beta blocker", "ACE inhibitor", "statin",
        "cardiac rehabilitation", "LV dysfunction",
    ],
    "acs": [
        "acute coronary syndrome", "unstable angina", "STEMI", "NSTEMI",
        "MI", "myocardial infarction", "troponin", "chest pain",
        "ECG changes", "ST elevation", "ST depression", "T wave inversion",
        "angiography", "PCI", "heparin", "aspirin", "P2Y12 inhibitor",
        "ticagrelor", "prasugrel", "clopidogrel", "DAPT",
        "risk stratification", "GRACE score", "TIMI score",
    ],
    "stemi": [
        "ST-elevation myocardial infarction", "MI", "acute MI",
        "primary PCI", "door-to-balloon time", "thrombolysis",
        "fibrinolysis", "reperfusion", "troponin", "ST elevation",
        "reciprocal changes", "culprit lesion", "total occlusion",
        "LAD", "RCA", "LCx", "cardiogenic shock", "Killip class",
    ],
    "nstemi": [
        "non-ST-elevation myocardial infarction", "NSTE-ACS",
        "unstable angina", "troponin elevation", "ST depression",
        "T wave inversion", "early invasive strategy",
        "ischemia-guided strategy", "angiography", "PCI",
        "anticoagulation", "DAPT", "GRACE score", "TIMI score",
        "risk stratification",
    ],
    "pci": [
        "percutaneous coronary intervention", "angioplasty", "stent",
        "DES", "drug-eluting stent", "BMS", "bare-metal stent",
        "balloon angioplasty", "FFR", "iFR", "IVUS",
        "intravascular ultrasound", "OCT", "optical coherence tomography",
        "DAPT", "dual antiplatelet therapy", "restenosis", "thrombosis",
        "stent thrombosis", "revascularization",
    ],
    "cabg": [
        "coronary artery bypass graft", "bypass surgery",
        "LIMA", "left internal mammary artery", "RIMA",
        "saphenous vein graft", "SVG", "off-pump", "on-pump",
        "cardiopulmonary bypass", "sternotomy", "MIDCAB",
        "revascularization", "graft patency", "multivessel disease",
        "left main disease",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 3. ARRHYTHMIA MAP
# ═══════════════════════════════════════════════════════════════════════

ARRHYTHMIA_MAP: Dict[str, List[str]] = {
    "atrial fibrillation": [
        "AFib", "AF", "irregular rhythm", "CHA2DS2-VASc",
        "anticoagulation", "DOAC", "warfarin", "apixaban", "rivaroxaban",
        "edoxaban", "dabigatran", "ablation", "pulmonary vein isolation",
        "PVI", "cardioversion", "rate control", "rhythm control",
        "amiodarone", "flecainide", "dronedarone", "sotalol",
        "left atrial appendage", "LAA", "Watchman", "stroke risk",
        "paroxysmal", "persistent", "permanent", "long-standing persistent",
    ],
    "afib": [
        "atrial fibrillation", "AF", "CHA2DS2-VASc", "anticoagulation",
        "DOAC", "warfarin", "ablation", "cardioversion", "rate control",
        "rhythm control", "amiodarone", "stroke risk", "irregular rhythm",
        "rapid ventricular response", "RVR", "LAA closure",
    ],
    "af": [
        "atrial fibrillation", "AFib", "CHA2DS2-VASc", "anticoagulation",
        "DOAC", "ablation", "cardioversion", "rate control",
        "rhythm control", "stroke risk",
    ],
    "vt": [
        "ventricular tachycardia", "VT", "wide complex tachycardia",
        "sustained VT", "non-sustained VT", "NSVT", "monomorphic VT",
        "polymorphic VT", "ICD", "implantable cardioverter defibrillator",
        "ablation", "amiodarone", "lidocaine", "procainamide",
        "electrical storm", "VF", "ventricular fibrillation",
        "sudden cardiac death", "SCD", "antitachycardia pacing", "ATP",
        "cardiomyopathy", "scar-mediated",
    ],
    "svt": [
        "supraventricular tachycardia", "SVT", "narrow complex tachycardia",
        "AVNRT", "AVRT", "atrial tachycardia", "adenosine",
        "vagal maneuvers", "Valsalva", "carotid sinus massage",
        "ablation", "verapamil", "diltiazem", "beta blocker",
        "accessory pathway", "WPW", "Wolff-Parkinson-White",
    ],
    "bradycardia": [
        "slow heart rate", "sinus bradycardia", "sinus node dysfunction",
        "sick sinus syndrome", "SSS", "heart block", "AV block",
        "first degree", "second degree", "Mobitz", "Wenckebach",
        "third degree", "complete heart block", "pacemaker",
        "chronotropic incompetence", "symptomatic bradycardia",
        "atropine", "isoproterenol", "temporary pacing",
    ],
    "tachycardia": [
        "fast heart rate", "sinus tachycardia", "SVT",
        "supraventricular tachycardia", "VT", "ventricular tachycardia",
        "atrial fibrillation", "atrial flutter", "wide complex",
        "narrow complex", "rate control", "rhythm control",
        "cardioversion", "adenosine", "amiodarone",
    ],
    "long qt": [
        "QT prolongation", "long QT syndrome", "LQTS", "torsades de pointes",
        "TdP", "QTc", "corrected QT", "drug-induced QT prolongation",
        "SCN5A", "KCNQ1", "KCNH2", "hERG", "channelopathy",
        "beta blocker", "ICD", "magnesium", "isoproterenol",
        "congenital long QT",
    ],
    "brugada": [
        "Brugada syndrome", "Brugada pattern", "SCN5A", "type 1 Brugada",
        "coved ST elevation", "right precordial", "V1", "V2",
        "sodium channel", "ICD", "sudden cardiac death", "SCD",
        "quinidine", "isoproterenol", "fever-induced", "ajmaline",
        "procainamide", "flecainide challenge",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 4. VALVULAR DISEASE MAP
# ═══════════════════════════════════════════════════════════════════════

VALVULAR_MAP: Dict[str, List[str]] = {
    "aortic stenosis": [
        "AS", "aortic valve", "calcific aortic stenosis",
        "mean gradient", "peak gradient", "AVA", "aortic valve area",
        "EOA", "effective orifice area", "dimensionless index",
        "TAVR", "transcatheter aortic valve replacement", "SAVR",
        "surgical aortic valve replacement", "AVR", "aortic valve replacement",
        "bicuspid aortic valve", "BAV", "low-flow low-gradient",
        "paradoxical low-flow", "dobutamine stress echo", "severe AS",
        "critical AS", "valve-in-valve",
    ],
    "as": [
        "aortic stenosis", "aortic valve", "gradient", "AVA",
        "TAVR", "SAVR", "calcific", "bicuspid", "valve replacement",
        "EOA", "severe", "critical", "low-flow low-gradient",
    ],
    "mitral regurgitation": [
        "MR", "mitral valve", "regurgitant volume", "regurgitant fraction",
        "EROA", "effective regurgitant orifice area", "vena contracta",
        "proximal isovelocity surface area", "PISA", "mitral valve repair",
        "MVR", "mitral valve replacement", "MitraClip", "TEER",
        "transcatheter edge-to-edge repair", "primary MR", "secondary MR",
        "functional MR", "degenerative MR", "flail leaflet",
        "mitral valve prolapse", "MVP", "chordal rupture",
        "annuloplasty", "ring",
    ],
    "mr": [
        "mitral regurgitation", "mitral valve", "regurgitant volume",
        "vena contracta", "EROA", "repair", "replacement", "MitraClip",
        "TEER", "primary MR", "secondary MR", "MVP",
    ],
    "tavr": [
        "transcatheter aortic valve replacement", "TAVI",
        "transcatheter aortic valve implantation", "aortic stenosis",
        "Edwards SAPIEN", "Medtronic CoreValve", "Evolut",
        "valve-in-valve", "paravalvular leak", "PVL",
        "conduction disturbance", "pacemaker", "vascular access",
        "transfemoral", "transapical", "CT planning", "annulus sizing",
        "STS score", "frailty", "heart team",
    ],
    "valve": [
        "valvular disease", "aortic valve", "mitral valve",
        "tricuspid valve", "pulmonic valve", "stenosis", "regurgitation",
        "repair", "replacement", "prosthetic valve", "mechanical valve",
        "bioprosthetic valve", "endocarditis", "vegetation",
        "anticoagulation", "INR", "valve-in-valve",
    ],
    "ross procedure": [
        "Ross operation", "pulmonary autograft", "aortic valve replacement",
        "homograft", "bicuspid aortic valve", "young patient",
        "autograft", "neoaortic root", "pulmonary valve",
        "durability", "reintervention",
    ],
    "tricuspid regurgitation": [
        "TR", "tricuspid valve", "severe TR", "functional TR",
        "annular dilation", "right heart failure", "RV dilation",
        "tricuspid annuloplasty", "TriClip", "TEER",
        "hepatic vein flow reversal", "peripheral edema", "ascites",
        "pulmonary hypertension", "right ventricular dysfunction",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 5. ECHOCARDIOGRAPHY / IMAGING - ECHO MAP
# ═══════════════════════════════════════════════════════════════════════

IMAGING_ECHO_MAP: Dict[str, List[str]] = {
    "echo": [
        "echocardiogram", "echocardiography", "TTE",
        "transthoracic echocardiogram", "TEE",
        "transesophageal echocardiogram", "LVEF",
        "left ventricular ejection fraction", "GLS",
        "global longitudinal strain", "diastolic function",
        "E/A ratio", "E/e' ratio", "TAPSE",
        "tricuspid annular plane systolic excursion",
        "chamber dimensions", "wall thickness", "IVSd", "LVIDd",
        "LVIDs", "PWd", "wall motion", "wall motion abnormality",
        "WMA", "speckle tracking", "strain imaging", "RV function",
        "RVSP", "right ventricular systolic pressure",
        "pulmonary artery systolic pressure", "PASP",
        "pericardial effusion", "Doppler", "color Doppler",
        "tissue Doppler", "contrast echo",
    ],
    "echocardiogram": [
        "echo", "TTE", "TEE", "LVEF", "GLS", "diastolic function",
        "E/A", "E/e'", "TAPSE", "wall motion", "strain",
        "chamber dimensions", "Doppler", "valvular assessment",
    ],
    "tte": [
        "transthoracic echocardiogram", "echo", "echocardiography",
        "LVEF", "GLS", "diastolic function", "wall motion",
        "valvular assessment", "chamber dimensions", "Doppler",
    ],
    "tee": [
        "transesophageal echocardiogram", "echo", "echocardiography",
        "intracardiac", "valvular assessment", "endocarditis",
        "vegetation", "atrial appendage", "thrombus", "PFO",
        "ASD", "prosthetic valve", "mitral valve", "aortic valve",
        "intraoperative", "procedural guidance",
    ],
    "strain": [
        "speckle tracking", "GLS", "global longitudinal strain",
        "myocardial deformation", "strain imaging", "strain rate",
        "radial strain", "circumferential strain", "longitudinal strain",
        "segmental strain", "bull's eye plot", "cardiotoxicity",
        "subclinical dysfunction", "CTRCD", "echocardiography",
    ],
    "lvef": [
        "left ventricular ejection fraction", "ejection fraction", "EF",
        "systolic function", "LV function", "biplane Simpson",
        "3D ejection fraction", "visual estimate", "wall motion",
        "contractility", "HFrEF", "HFpEF", "HFmrEF",
        "cardiomyopathy", "recovery", "improvement",
    ],
    "diastolic function": [
        "diastolic dysfunction", "E/A ratio", "E/e' ratio",
        "e' velocity", "deceleration time", "DT", "IVRT",
        "isovolumic relaxation time", "LA volume", "LA volume index",
        "TR velocity", "pulmonary vein flow", "grade I", "grade II",
        "grade III", "relaxation abnormality", "pseudonormal",
        "restrictive filling", "elevated filling pressures",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 6. CARDIAC CT MAP
# ═══════════════════════════════════════════════════════════════════════

IMAGING_CT_MAP: Dict[str, List[str]] = {
    "cardiac ct": [
        "coronary CT angiography", "CTA", "calcium score",
        "Agatston score", "CAD-RADS", "plaque", "stenosis",
        "coronary artery", "LAD", "left anterior descending",
        "RCA", "right coronary artery", "LCx", "left circumflex",
        "left main", "calcified plaque", "non-calcified plaque",
        "mixed plaque", "positive remodeling", "napkin-ring sign",
        "FFR-CT", "CT-derived FFR", "TAVR planning", "annulus sizing",
        "structural CT", "gating", "prospective", "retrospective",
    ],
    "cta": [
        "coronary CT angiography", "cardiac CT", "stenosis", "plaque",
        "CAD-RADS", "coronary artery", "LAD", "RCA", "LCx",
        "FFR-CT", "non-invasive", "anatomic assessment",
        "contrast", "gating",
    ],
    "calcium score": [
        "coronary artery calcium", "CAC", "Agatston score",
        "calcium scoring", "cardiac CT", "risk stratification",
        "ASCVD risk", "zero calcium score", "high calcium score",
        "percentile", "age-sex-race percentile",
        "CAC greater than 100", "CAC greater than 400",
        "statin indication", "preventive cardiology",
    ],
    "cad-rads": [
        "Coronary Artery Disease Reporting and Data System",
        "CAD-RADS 0", "CAD-RADS 1", "CAD-RADS 2", "CAD-RADS 3",
        "CAD-RADS 4A", "CAD-RADS 4B", "CAD-RADS 5", "CAD-RADS N",
        "stenosis grade", "plaque burden", "high-risk plaque",
        "coronary CTA", "reporting", "standardized reporting",
    ],
    "ffr-ct": [
        "CT-derived fractional flow reserve", "FFR-CT", "HeartFlow",
        "coronary CTA", "functional assessment", "non-invasive FFR",
        "hemodynamic significance", "ischemia", "stenosis",
        "computational fluid dynamics",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 7. CARDIAC MRI MAP
# ═══════════════════════════════════════════════════════════════════════

IMAGING_MRI_MAP: Dict[str, List[str]] = {
    "cardiac mri": [
        "CMR", "cardiac magnetic resonance", "cine imaging",
        "T1 mapping", "T2 mapping", "ECV", "extracellular volume",
        "LGE", "late gadolinium enhancement", "perfusion",
        "stress perfusion", "rest perfusion", "adenosine stress",
        "regadenoson", "tissue characterization", "fibrosis",
        "edema", "myocarditis", "cardiomyopathy", "iron overload",
        "T2*", "amyloid", "sarcoid", "ARVC",
        "arrhythmogenic right ventricular cardiomyopathy",
        "strain", "feature tracking", "phase contrast",
        "flow quantification", "4D flow",
    ],
    "cmr": [
        "cardiac MRI", "cardiac magnetic resonance", "LGE",
        "T1 mapping", "T2 mapping", "ECV", "tissue characterization",
        "perfusion", "fibrosis", "edema", "cardiomyopathy",
        "myocarditis", "strain", "feature tracking",
    ],
    "lge": [
        "late gadolinium enhancement", "LGE", "gadolinium",
        "delayed enhancement", "scar", "fibrosis",
        "subendocardial", "transmural", "mid-wall", "epicardial",
        "RV insertion point", "diffuse", "patchy",
        "ischemic pattern", "non-ischemic pattern",
        "myocardial infarction", "cardiomyopathy", "myocarditis",
        "HCM", "DCM", "sarcoid", "amyloid",
    ],
    "t1 mapping": [
        "T1 map", "native T1", "post-contrast T1", "ECV",
        "extracellular volume fraction", "myocardial characterization",
        "diffuse fibrosis", "amyloid", "Anderson-Fabry",
        "iron overload", "edema", "MOLLI", "ShMOLLI",
        "parametric mapping",
    ],
    "t2 mapping": [
        "T2 map", "myocardial edema", "inflammation", "myocarditis",
        "acute MI", "area at risk", "tissue characterization",
        "T2-weighted imaging", "STIR", "parametric mapping",
    ],
    "tissue characterization": [
        "T1 mapping", "T2 mapping", "ECV", "LGE",
        "late gadolinium enhancement", "fibrosis", "edema",
        "inflammation", "iron overload", "T2*", "amyloid",
        "Anderson-Fabry", "sarcoid", "myocarditis",
        "parametric mapping", "native T1", "native T2",
    ],
    "perfusion": [
        "stress perfusion", "rest perfusion", "adenosine",
        "regadenoson", "dipyridamole", "vasodilator stress",
        "myocardial blood flow", "MBF", "ischemia", "perfusion defect",
        "subendocardial", "transmural", "coronary reserve",
        "microvascular disease", "CMR", "cardiac MRI",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 8. PREVENTIVE CARDIOLOGY MAP
# ═══════════════════════════════════════════════════════════════════════

PREVENTIVE_MAP: Dict[str, List[str]] = {
    "prevention": [
        "primary prevention", "secondary prevention", "ASCVD",
        "atherosclerotic cardiovascular disease", "risk assessment",
        "risk factor", "lifestyle modification", "diet", "exercise",
        "Mediterranean diet", "DASH diet", "weight management",
        "smoking cessation", "blood pressure control",
        "lipid management", "statin", "aspirin",
    ],
    "risk": [
        "cardiovascular risk", "ASCVD risk", "PCE",
        "pooled cohort equations", "10-year risk", "lifetime risk",
        "risk factors", "risk stratification", "risk reduction",
        "risk enhancers", "coronary artery calcium",
        "Framingham risk score", "SCORE", "SCORE2",
    ],
    "cholesterol": [
        "LDL", "low-density lipoprotein", "LDL-C",
        "HDL", "high-density lipoprotein", "HDL-C",
        "total cholesterol", "triglycerides", "VLDL",
        "non-HDL cholesterol", "Lp(a)", "lipoprotein(a)",
        "ApoB", "apolipoprotein B", "lipid panel",
        "lipid profile", "dyslipidemia", "hyperlipidemia",
        "familial hypercholesterolemia", "FH",
    ],
    "lipids": [
        "cholesterol", "LDL", "HDL", "triglycerides", "Lp(a)",
        "ApoB", "lipid panel", "dyslipidemia", "statin",
        "ezetimibe", "PCSK9 inhibitor", "bempedoic acid",
        "inclisiran", "icosapent ethyl", "omega-3",
        "lipid management", "lipid-lowering therapy",
    ],
    "statin": [
        "HMG-CoA reductase inhibitor", "atorvastatin", "rosuvastatin",
        "simvastatin", "pravastatin", "pitavastatin", "lovastatin",
        "fluvastatin", "high-intensity statin", "moderate-intensity statin",
        "statin benefit group", "statin intolerance",
        "myopathy", "rhabdomyolysis", "LDL reduction",
        "ASCVD risk reduction", "primary prevention",
        "secondary prevention", "ezetimibe", "PCSK9i",
    ],
    "ascvd": [
        "atherosclerotic cardiovascular disease", "10-year ASCVD risk",
        "PCE", "pooled cohort equations", "risk assessment",
        "risk enhancers", "coronary artery calcium", "statin therapy",
        "LDL", "cholesterol", "lifestyle", "prevention",
        "risk factor modification",
    ],
    "pcsk9": [
        "PCSK9 inhibitor", "PCSK9i", "evolocumab", "alirocumab",
        "Repatha", "Praluent", "LDL lowering", "familial hypercholesterolemia",
        "FH", "statin intolerance", "ASCVD", "secondary prevention",
        "LDL goal", "very high risk",
    ],
    "lpa": [
        "lipoprotein(a)", "Lp(a)", "lipoprotein little a",
        "elevated Lp(a)", "cardiovascular risk", "aortic stenosis",
        "ASCVD", "risk enhancer", "pelacarsen", "olpasiran",
        "antisense oligonucleotide", "siRNA", "inherited risk",
        "familial", "LPA gene",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 9. HEMODYNAMICS MAP
# ═══════════════════════════════════════════════════════════════════════

HEMODYNAMICS_MAP: Dict[str, List[str]] = {
    "catheterization": [
        "cardiac catheterization", "cath", "heart cath",
        "right heart catheterization", "RHC", "left heart catheterization",
        "LHC", "coronary angiography", "hemodynamics", "pressure",
        "LVEDP", "left ventricular end-diastolic pressure",
        "PCWP", "pulmonary capillary wedge pressure",
        "PA pressure", "pulmonary artery pressure",
        "cardiac output", "cardiac index", "PVR",
        "pulmonary vascular resistance", "SVR",
        "systemic vascular resistance", "Fick", "thermodilution",
        "FFR", "fractional flow reserve", "shunt calculation",
        "oximetry run",
    ],
    "cath": [
        "catheterization", "cardiac catheterization", "angiography",
        "hemodynamics", "pressure", "LVEDP", "PCWP", "PA pressure",
        "cardiac output", "FFR", "PCI", "coronary angiography",
    ],
    "hemodynamics": [
        "cardiac catheterization", "pressure", "LVEDP", "PCWP",
        "PA pressure", "RA pressure", "RV pressure",
        "cardiac output", "CO", "cardiac index", "CI",
        "PVR", "pulmonary vascular resistance", "SVR",
        "systemic vascular resistance", "Fick method",
        "thermodilution", "stroke volume", "SV",
        "stroke volume index", "SVI", "mixed venous saturation",
        "SvO2", "transpulmonary gradient", "TPG",
        "diastolic pressure gradient", "DPG",
    ],
    "pressure": [
        "blood pressure", "hemodynamics", "LVEDP", "PCWP",
        "PA pressure", "systolic pressure", "diastolic pressure",
        "mean pressure", "RA pressure", "RV pressure",
        "aortic pressure", "gradient", "catheterization",
    ],
    "lvedp": [
        "left ventricular end-diastolic pressure", "filling pressure",
        "hemodynamics", "catheterization", "diastolic function",
        "heart failure", "elevated filling pressures",
        "congestion", "PCWP", "wedge pressure",
    ],
    "pcwp": [
        "pulmonary capillary wedge pressure", "wedge pressure",
        "left atrial pressure", "filling pressure", "hemodynamics",
        "right heart catheterization", "heart failure", "congestion",
        "pulmonary edema", "volume status", "LVEDP",
    ],
    "cardiac output": [
        "CO", "cardiac index", "CI", "stroke volume", "SV",
        "heart rate", "Fick method", "thermodilution",
        "hemodynamics", "catheterization", "low cardiac output",
        "cardiogenic shock", "inotropes", "dobutamine",
        "milrinone", "vasopressors",
    ],
    "ffr": [
        "fractional flow reserve", "FFR", "coronary physiology",
        "iFR", "instantaneous wave-free ratio", "pressure wire",
        "hyperemia", "adenosine", "hemodynamic significance",
        "ischemia", "stenosis severity", "FFR less than 0.80",
        "revascularization decision", "PCI", "deferral",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 10. ELECTROPHYSIOLOGY MAP
# ═══════════════════════════════════════════════════════════════════════

ELECTROPHYSIOLOGY_MAP: Dict[str, List[str]] = {
    "ecg": [
        "electrocardiogram", "EKG", "12-lead ECG", "rhythm",
        "heart rhythm", "sinus rhythm", "axis", "P wave", "PR interval",
        "QRS complex", "QRS duration", "QTc", "corrected QT interval",
        "ST segment", "ST elevation", "ST depression", "T wave",
        "T wave inversion", "LBBB", "left bundle branch block",
        "RBBB", "right bundle branch block", "LAFB",
        "left anterior fascicular block", "LPFB",
        "left posterior fascicular block", "bifascicular block",
        "LVH", "left ventricular hypertrophy", "RVH",
        "right ventricular hypertrophy", "Q waves", "R wave progression",
        "poor R wave progression", "U wave",
    ],
    "ekg": [
        "ECG", "electrocardiogram", "12-lead", "rhythm", "axis",
        "intervals", "PR", "QRS", "QTc", "ST segment", "T wave",
        "LBBB", "RBBB", "LVH", "RVH",
    ],
    "ep study": [
        "electrophysiology study", "EP", "intracardiac electrograms",
        "mapping", "3D mapping", "CARTO", "EnSite",
        "programmed stimulation", "inducibility", "ablation",
        "radiofrequency ablation", "RF ablation", "cryoablation",
        "pulsed field ablation", "PFA", "reentrant circuit",
        "accessory pathway", "slow pathway", "fast pathway",
        "His bundle", "AV node", "sinus node",
    ],
    "electrophysiology": [
        "EP study", "ablation", "mapping", "arrhythmia",
        "rhythm", "conduction", "intracardiac", "programmed stimulation",
        "radiofrequency", "cryoablation", "pulsed field ablation",
        "catheter ablation", "3D mapping",
    ],
    "ablation": [
        "catheter ablation", "radiofrequency ablation", "RF ablation",
        "cryoablation", "pulsed field ablation", "PFA",
        "pulmonary vein isolation", "PVI", "atrial fibrillation ablation",
        "VT ablation", "SVT ablation", "AVNRT ablation",
        "accessory pathway ablation", "mapping", "3D mapping",
        "substrate ablation", "pace mapping", "entrainment",
        "success rate", "recurrence", "complication",
    ],
    "qrs": [
        "QRS complex", "QRS duration", "wide QRS", "narrow QRS",
        "bundle branch block", "LBBB", "RBBB", "intraventricular conduction delay",
        "IVCD", "ventricular depolarization", "CRT criteria",
        "QRS morphology", "delta wave", "pre-excitation",
    ],
    "qtc": [
        "corrected QT interval", "QT prolongation", "long QT syndrome",
        "LQTS", "short QT syndrome", "SQTS", "torsades de pointes",
        "TdP", "drug-induced QT prolongation", "Bazett formula",
        "Fridericia formula", "QTc greater than 500",
        "QTc greater than 470", "risk of arrhythmia",
    ],
    "lbbb": [
        "left bundle branch block", "QRS greater than 120",
        "broad QRS", "CRT indication", "cardiac resynchronization",
        "new LBBB", "acute MI", "rate-related LBBB",
        "conduction abnormality", "His bundle pacing",
        "left bundle branch area pacing", "LBBAP",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 11. CARDIO-ONCOLOGY MAP
# ═══════════════════════════════════════════════════════════════════════

CARDIO_ONCOLOGY_MAP: Dict[str, List[str]] = {
    "cardiotoxicity": [
        "cancer therapy-related cardiac dysfunction", "CTRCD",
        "anthracycline", "doxorubicin", "trastuzumab", "Herceptin",
        "checkpoint inhibitor", "immune checkpoint inhibitor", "ICI",
        "myocarditis", "ICI myocarditis", "tyrosine kinase inhibitor",
        "TKI", "VEGF inhibitor", "proteasome inhibitor",
        "carfilzomib", "radiation-induced heart disease",
        "GLS", "global longitudinal strain", "LVEF monitoring",
        "baseline LVEF", "dexrazoxane", "cardioprotection",
        "surveillance", "serial echocardiography",
    ],
    "chemotherapy": [
        "anthracycline", "doxorubicin", "epirubicin", "daunorubicin",
        "trastuzumab", "pertuzumab", "cardiotoxicity", "CTRCD",
        "LVEF decline", "GLS decline", "heart failure",
        "troponin", "BNP", "cardioprotection", "dexrazoxane",
        "monitoring", "cumulative dose", "lifetime dose",
    ],
    "cardio-oncology": [
        "cardiotoxicity", "cancer therapy", "CTRCD",
        "cardio-oncology assessment", "baseline risk",
        "cardiovascular risk", "anthracycline", "trastuzumab",
        "checkpoint inhibitor", "radiation therapy",
        "radiation-induced heart disease", "CAD", "valvular disease",
        "pericardial disease", "arrhythmia", "GLS monitoring",
        "LVEF surveillance", "cardioprotection", "survivorship",
    ],
    "checkpoint inhibitor": [
        "immune checkpoint inhibitor", "ICI", "pembrolizumab",
        "nivolumab", "ipilimumab", "atezolizumab", "durvalumab",
        "ICI myocarditis", "immune-mediated myocarditis",
        "fulminant myocarditis", "troponin", "ECG changes",
        "conduction abnormality", "high-dose steroids",
        "methylprednisolone", "immunosuppression", "cardiotoxicity",
    ],
    "anthracycline": [
        "doxorubicin", "epirubicin", "daunorubicin", "idarubicin",
        "anthracycline cardiotoxicity", "cumulative dose",
        "dose-dependent", "CTRCD", "LVEF decline",
        "GLS decline", "dexrazoxane", "cardioprotection",
        "dilated cardiomyopathy", "heart failure", "irreversible",
        "monitoring", "serial echocardiography",
    ],
    "radiation heart disease": [
        "radiation-induced heart disease", "RIHD",
        "radiation-associated cardiovascular disease",
        "mediastinal radiation", "chest radiation",
        "CAD", "valvular disease", "pericarditis",
        "constrictive pericarditis", "restrictive cardiomyopathy",
        "conduction disease", "carotid stenosis",
        "surveillance", "screening", "late effects",
    ],
    "gls": [
        "global longitudinal strain", "GLS", "speckle tracking",
        "myocardial strain", "GLS decline", "relative GLS drop",
        "15 percent drop", "subclinical dysfunction", "CTRCD",
        "cardiotoxicity", "anthracycline", "trastuzumab",
        "surveillance", "echocardiography", "baseline GLS",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 12. CONGENITAL HEART DISEASE MAP
# ═══════════════════════════════════════════════════════════════════════

CONGENITAL_MAP: Dict[str, List[str]] = {
    "congenital": [
        "congenital heart disease", "CHD", "structural heart disease",
        "pediatric cardiology", "adult congenital heart disease", "ACHD",
        "ASD", "atrial septal defect", "VSD", "ventricular septal defect",
        "PDA", "patent ductus arteriosus", "PFO", "patent foramen ovale",
        "coarctation", "aortic coarctation", "Tetralogy of Fallot", "TOF",
        "transposition of the great arteries", "TGA", "Ebstein anomaly",
        "bicuspid aortic valve", "BAV", "shunt", "cyanotic", "acyanotic",
    ],
    "chd": [
        "congenital heart disease", "structural heart disease",
        "adult congenital heart disease", "ACHD", "ASD", "VSD", "PDA",
        "TOF", "TGA", "shunt", "cyanotic", "acyanotic",
        "surgical repair", "palliation",
    ],
    "asd": [
        "atrial septal defect", "secundum ASD", "primum ASD",
        "sinus venosus defect", "left-to-right shunt", "Qp:Qs",
        "right heart volume overload", "RV dilation",
        "paradoxical embolism", "device closure", "Amplatzer",
        "surgical closure", "Eisenmenger syndrome",
    ],
    "vsd": [
        "ventricular septal defect", "perimembranous VSD",
        "muscular VSD", "inlet VSD", "outlet VSD",
        "left-to-right shunt", "Qp:Qs", "LV volume overload",
        "Eisenmenger syndrome", "device closure", "surgical closure",
        "spontaneous closure",
    ],
    "pfo": [
        "patent foramen ovale", "PFO closure", "cryptogenic stroke",
        "paradoxical embolism", "right-to-left shunt",
        "bubble study", "agitated saline", "contrast echo",
        "device closure", "Amplatzer PFO Occluder",
        "GORE CARDIOFORM", "migraine",
    ],
    "tetralogy": [
        "Tetralogy of Fallot", "TOF", "VSD",
        "right ventricular outflow tract obstruction", "RVOTO",
        "overriding aorta", "right ventricular hypertrophy",
        "cyanosis", "tet spell", "Blalock-Taussig shunt", "BT shunt",
        "complete repair", "pulmonary regurgitation",
        "pulmonary valve replacement", "RV dilation",
        "RV dysfunction", "arrhythmia",
    ],
    "eisenmenger": [
        "Eisenmenger syndrome", "pulmonary arterial hypertension",
        "PAH", "shunt reversal", "right-to-left shunt", "cyanosis",
        "irreversible pulmonary vascular disease",
        "pulmonary vascular resistance", "PVR",
        "bosentan", "sildenafil", "prostacyclin",
        "contraindication to closure",
    ],
    "fontan": [
        "Fontan circulation", "Fontan procedure", "total cavopulmonary connection",
        "TCPC", "Glenn", "bidirectional Glenn", "single ventricle",
        "hypoplastic left heart syndrome", "HLHS",
        "protein-losing enteropathy", "PLE", "plastic bronchitis",
        "Fontan-associated liver disease", "FALD",
        "failing Fontan", "heart transplant",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 13. CARDIAC DEVICE MAP
# ═══════════════════════════════════════════════════════════════════════

DEVICE_MAP: Dict[str, List[str]] = {
    "pacemaker": [
        "PPM", "permanent pacemaker", "dual chamber pacemaker", "DDD",
        "single chamber pacemaker", "VVI", "AAI", "leadless pacemaker",
        "Micra", "lead", "atrial lead", "ventricular lead",
        "sensing", "pacing", "threshold", "impedance", "battery",
        "battery longevity", "ERI", "elective replacement indicator",
        "EOL", "end of life", "interrogation", "device check",
        "rate response", "mode switch", "His bundle pacing",
        "left bundle branch area pacing", "LBBAP",
        "conduction system pacing", "CSP",
    ],
    "icd": [
        "implantable cardioverter defibrillator", "ICD",
        "primary prevention ICD", "secondary prevention ICD",
        "shock", "appropriate shock", "inappropriate shock",
        "antitachycardia pacing", "ATP", "VT detection",
        "VF detection", "sensing", "lead", "impedance", "threshold",
        "battery", "interrogation", "S-ICD", "subcutaneous ICD",
        "transvenous ICD", "programming", "shock zone", "VT zone",
    ],
    "crt": [
        "cardiac resynchronization therapy", "CRT",
        "CRT-D", "CRT-P", "biventricular pacing", "BiV",
        "left ventricular lead", "coronary sinus lead", "CS lead",
        "LBBB", "QRS greater than 150", "super-responder",
        "non-responder", "optimization", "AV delay", "VV delay",
        "heart failure", "LVEF less than 35", "CRT response",
    ],
    "lvad": [
        "left ventricular assist device", "LVAD", "mechanical circulatory support",
        "MCS", "HeartMate", "HeartMate 3", "HeartWare", "HVAD",
        "destination therapy", "DT", "bridge to transplant", "BTT",
        "bridge to decision", "BTD", "driveline", "pump speed",
        "pump flow", "power", "pulsatility index", "PI",
        "suction event", "pump thrombosis", "driveline infection",
        "right heart failure", "GI bleeding", "stroke",
        "hemocompatibility", "anticoagulation",
    ],
    "lead": [
        "pacemaker lead", "ICD lead", "CRT lead", "atrial lead",
        "ventricular lead", "coronary sinus lead", "lead extraction",
        "lead malfunction", "lead fracture", "lead dislodgement",
        "insulation breach", "sensing", "impedance", "threshold",
        "capture", "undersensing", "oversensing",
    ],
    "interrogation": [
        "device interrogation", "device check", "remote monitoring",
        "Medtronic CareLink", "Abbott Merlin", "Boston Scientific LATITUDE",
        "Biotronik Home Monitoring", "battery status", "lead parameters",
        "arrhythmia episodes", "therapy history", "shock history",
        "pacing percentage", "mode", "programming",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 14. VASCULAR MAP
# ═══════════════════════════════════════════════════════════════════════

VASCULAR_MAP: Dict[str, List[str]] = {
    "aorta": [
        "aortic disease", "aortic aneurysm", "AAA",
        "abdominal aortic aneurysm", "TAA",
        "thoracic aortic aneurysm", "aortic dissection",
        "Stanford type A", "Stanford type B", "DeBakey",
        "aortic root", "ascending aorta", "aortic arch",
        "descending aorta", "abdominal aorta",
        "aortic root dilation", "annuloaortic ectasia",
        "Marfan syndrome", "Loeys-Dietz syndrome",
        "Ehlers-Danlos vascular type", "bicuspid aortopathy",
        "surveillance imaging", "surgical threshold",
    ],
    "pad": [
        "peripheral artery disease", "peripheral arterial disease",
        "peripheral vascular disease", "PVD", "claudication",
        "intermittent claudication", "critical limb ischemia", "CLI",
        "chronic limb-threatening ischemia", "CLTI",
        "ABI", "ankle-brachial index", "exercise ABI",
        "arterial duplex", "angiography", "CTA", "MRA",
        "revascularization", "angioplasty", "stent",
        "bypass surgery", "cilostazol", "antiplatelet",
        "supervised exercise", "wound healing",
    ],
    "peripheral": [
        "peripheral artery disease", "PAD", "PVD", "claudication",
        "limb ischemia", "ABI", "duplex ultrasound", "angiography",
        "revascularization", "bypass", "stent", "antiplatelet",
    ],
    "carotid": [
        "carotid artery disease", "carotid stenosis",
        "carotid artery stenosis", "internal carotid artery", "ICA",
        "carotid endarterectomy", "CEA", "carotid stenting", "CAS",
        "carotid duplex", "carotid ultrasound", "IMT",
        "intima-media thickness", "carotid plaque",
        "stroke prevention", "TIA", "transient ischemic attack",
        "symptomatic carotid stenosis", "asymptomatic carotid stenosis",
        "NASCET criteria",
    ],
    "aneurysm": [
        "aortic aneurysm", "AAA", "abdominal aortic aneurysm",
        "TAA", "thoracic aortic aneurysm", "TAAA",
        "thoracoabdominal aortic aneurysm", "saccular", "fusiform",
        "dilation", "rupture risk", "surgical repair",
        "open repair", "EVAR", "endovascular aneurysm repair",
        "TEVAR", "thoracic endovascular aortic repair",
        "stent-graft", "endovascular", "surveillance",
        "growth rate", "5.5 cm threshold",
    ],
    "dissection": [
        "aortic dissection", "Stanford type A", "Stanford type B",
        "DeBakey classification", "intimal tear", "false lumen",
        "true lumen", "malperfusion", "aortic root",
        "ascending aorta", "descending aorta",
        "emergency surgery", "TEVAR", "medical management",
        "blood pressure control", "heart rate control",
        "esmolol", "labetalol", "nicardipine",
        "complicated dissection", "uncomplicated dissection",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 15. CARDIOVASCULAR GENOMICS MAP
# ═══════════════════════════════════════════════════════════════════════

GENOMICS_CARDIO_MAP: Dict[str, List[str]] = {
    "cardiomyopathy": [
        "dilated cardiomyopathy", "DCM", "hypertrophic cardiomyopathy",
        "HCM", "restrictive cardiomyopathy", "RCM",
        "arrhythmogenic cardiomyopathy", "ACM", "ARVC",
        "left ventricular non-compaction", "LVNC",
        "MYH7", "MYBPC3", "TTN", "LMNA", "TNNT2", "TNNI3",
        "DSP", "PKP2", "DSG2", "DSC2",
        "genetic testing", "gene panel", "cascade screening",
        "family screening", "variant", "pathogenic variant",
        "likely pathogenic", "VUS", "variant of uncertain significance",
        "penetrance", "expressivity", "autosomal dominant",
    ],
    "channelopathy": [
        "ion channel disease", "long QT syndrome", "LQTS",
        "Brugada syndrome", "catecholaminergic polymorphic VT", "CPVT",
        "short QT syndrome", "SQTS", "early repolarization syndrome",
        "SCN5A", "KCNQ1", "KCNH2", "hERG", "KCNJ2", "CACNA1C",
        "RYR2", "CASQ2", "sodium channel", "potassium channel",
        "calcium channel", "genetic testing", "gene panel",
        "cascade screening", "sudden cardiac death", "SCD",
        "ICD", "beta blocker", "flecainide",
    ],
    "marfan": [
        "Marfan syndrome", "FBN1", "fibrillin-1", "connective tissue disorder",
        "aortic root dilation", "aortic dissection",
        "mitral valve prolapse", "MVP", "lens subluxation",
        "tall stature", "arm span", "Ghent criteria",
        "aortic root replacement", "Bentall procedure",
        "valve-sparing root replacement", "David procedure",
        "losartan", "beta blocker", "exercise restriction",
        "surveillance imaging", "Loeys-Dietz syndrome",
        "Ehlers-Danlos vascular type", "TGFBR1", "TGFBR2",
    ],
    "genetic": [
        "genetic testing", "gene panel", "whole exome sequencing",
        "WES", "whole genome sequencing", "WGS",
        "variant classification", "pathogenic", "likely pathogenic",
        "VUS", "likely benign", "benign", "ACMG criteria",
        "ClinGen", "ClinVar", "gnomAD",
        "cascade screening", "family screening", "genetic counseling",
        "autosomal dominant", "autosomal recessive",
        "X-linked", "de novo", "penetrance", "expressivity",
        "genotype-phenotype correlation",
    ],
    "myh7": [
        "MYH7", "beta-myosin heavy chain", "hypertrophic cardiomyopathy",
        "HCM", "dilated cardiomyopathy", "DCM",
        "sarcomeric gene", "missense variant", "pathogenic variant",
        "genetic testing", "cascade screening", "family screening",
        "phenotype", "penetrance", "risk stratification",
        "sudden cardiac death", "SCD", "ICD",
    ],
    "mybpc3": [
        "MYBPC3", "myosin-binding protein C", "cardiac myosin-binding protein C",
        "hypertrophic cardiomyopathy", "HCM",
        "sarcomeric gene", "truncating variant", "frameshift",
        "nonsense", "pathogenic variant", "genetic testing",
        "cascade screening", "founder mutation",
        "late-onset HCM", "variable penetrance",
    ],
    "ttn": [
        "TTN", "titin", "dilated cardiomyopathy", "DCM",
        "truncating variant", "TTNtv", "titin-truncating variant",
        "A-band", "sarcomere", "pathogenic variant",
        "genetic testing", "peripartum cardiomyopathy",
        "alcohol-related cardiomyopathy", "prognosis",
        "recovery", "LVEF recovery",
    ],
    "scn5a": [
        "SCN5A", "sodium channel", "Nav1.5",
        "Brugada syndrome", "long QT syndrome type 3", "LQT3",
        "progressive cardiac conduction disease", "PCCD",
        "Lenegre disease", "sick sinus syndrome",
        "dilated cardiomyopathy", "overlap syndrome",
        "gain of function", "loss of function",
        "genetic testing", "variant", "pathogenic variant",
    ],
    "kcnq1": [
        "KCNQ1", "potassium channel", "Kv7.1",
        "long QT syndrome type 1", "LQT1",
        "Jervell and Lange-Nielsen syndrome", "JLNS",
        "Romano-Ward syndrome", "exercise-triggered arrhythmia",
        "swimming-triggered", "beta blocker", "ICD",
        "genetic testing", "cascade screening",
    ],
    "ldlr": [
        "LDLR", "LDL receptor", "familial hypercholesterolemia", "FH",
        "heterozygous FH", "HeFH", "homozygous FH", "HoFH",
        "elevated LDL", "premature CAD", "xanthoma", "tendon xanthoma",
        "corneal arcus", "statin therapy", "PCSK9 inhibitor",
        "ezetimibe", "LDL apheresis", "genetic testing",
        "cascade screening", "Dutch Lipid Network criteria",
        "Simon Broome criteria",
    ],
    "fbn1": [
        "FBN1", "fibrillin-1", "Marfan syndrome",
        "aortic root dilation", "aortic dissection",
        "connective tissue disorder", "Ghent criteria",
        "missense variant", "haploinsufficiency",
        "dominant negative", "genetic testing",
        "cascade screening", "neonatal Marfan",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 16. SPECIALTY SYNONYM MAP
# ═══════════════════════════════════════════════════════════════════════

SPECIALTY_SYNONYM_MAP: Dict[str, List[str]] = {
    "scad": [
        "spontaneous coronary artery dissection", "SCAD",
        "coronary dissection",
    ],
    "cmd": [
        "coronary microvascular disease", "microvascular angina",
        "syndrome X cardiac",
    ],
    "ppcm": [
        "peripartum cardiomyopathy", "postpartum cardiomyopathy",
        "pregnancy cardiomyopathy",
    ],
    "lvnc": [
        "left ventricular non-compaction", "noncompaction cardiomyopathy",
        "spongy myocardium",
    ],
    "attr": [
        "transthyretin amyloidosis", "ATTR cardiomyopathy", "ATTR-CM",
        "wild-type ATTR", "hereditary ATTR", "tafamidis",
    ],
    "pfo": [
        "patent foramen ovale", "PFO", "right-to-left shunt",
        "cryptogenic stroke cardiac",
    ],
    "bav": [
        "bicuspid aortic valve", "BAV", "bicuspid aortopathy",
        "congenital aortic valve",
    ],
    "finerenone": [
        "finerenone", "Kerendia", "non-steroidal MRA", "nsMRA",
    ],
    "semaglutide": [
        "semaglutide", "Wegovy", "Ozempic", "GLP-1 cardiovascular",
    ],
    "inclisiran": [
        "inclisiran", "Leqvio", "siRNA PCSK9",
    ],
    "tafamidis": [
        "tafamidis", "Vyndamax", "Vyndaqel", "TTR stabilizer",
    ],
    "pfa": [
        "pulsed field ablation", "PFA", "electroporation ablation",
    ],
    "ecpr": [
        "ECMO CPR", "extracorporeal CPR", "ECPR",
        "VA-ECMO resuscitation",
    ],
    "teer": [
        "transcatheter edge-to-edge repair", "TEER", "MitraClip",
        "PASCAL", "TriClip",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 17. RARE DISEASE MAP
# ═══════════════════════════════════════════════════════════════════════

RARE_DISEASE_MAP: Dict[str, List[str]] = {
    "fabry": [
        "Fabry disease", "Anderson-Fabry",
        "alpha-galactosidase A deficiency", "GLA mutation", "agalsidase",
    ],
    "chagas": [
        "Chagas disease", "Chagas cardiomyopathy",
        "Trypanosoma cruzi", "American trypanosomiasis",
    ],
    "noonan": [
        "Noonan syndrome", "RASopathy", "RAF1", "PTPN11",
        "congenital heart disease genetic",
    ],
    "marfan": [
        "Marfan syndrome", "FBN1", "fibrillin-1",
        "aortic root dilation", "ectopia lentis",
    ],
    "loeys_dietz": [
        "Loeys-Dietz syndrome", "TGFBR1", "TGFBR2", "SMAD3",
        "arterial tortuosity",
    ],
    "williams": [
        "Williams syndrome", "ELN deletion",
        "supravalvar aortic stenosis", "SVAS",
    ],
    "turner": [
        "Turner syndrome", "45X", "bicuspid aortic valve Turner",
        "coarctation Turner",
    ],
    "barth": [
        "Barth syndrome", "TAZ", "tafazzin",
        "3-methylglutaconic aciduria", "LVNC X-linked",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 18. PROCEDURE MAP
# ═══════════════════════════════════════════════════════════════════════

PROCEDURE_MAP: Dict[str, List[str]] = {
    "tavr": [
        "transcatheter aortic valve replacement", "TAVR", "TAVI",
        "Edwards SAPIEN", "Medtronic Evolut", "valve-in-valve",
    ],
    "pci": [
        "percutaneous coronary intervention", "PCI", "angioplasty",
        "stenting", "DES", "drug-eluting stent", "balloon angioplasty",
    ],
    "cabg": [
        "coronary artery bypass grafting", "CABG", "bypass surgery",
        "LIMA", "SVG", "arterial graft",
    ],
    "ablation": [
        "catheter ablation", "radiofrequency ablation", "RFA",
        "cryoablation", "PFA", "pulmonary vein isolation", "PVI",
    ],
    "icd_implant": [
        "ICD implantation", "defibrillator implant", "S-ICD",
        "subcutaneous ICD", "transvenous ICD",
    ],
    "crt_implant": [
        "CRT implantation", "biventricular pacemaker", "CRT-D",
        "CRT-P", "cardiac resynchronization",
    ],
    "laa_occlusion": [
        "LAA occlusion", "left atrial appendage closure", "WATCHMAN",
        "Amulet", "LARIAT",
    ],
    "lvad": [
        "left ventricular assist device", "LVAD", "HeartMate 3",
        "mechanical circulatory support", "MCS", "destination therapy",
        "bridge to transplant",
    ],
    "transplant": [
        "heart transplant", "cardiac transplantation", "OHT",
        "orthotopic heart transplant", "donor heart",
    ],
    "alcohol_septal": [
        "alcohol septal ablation", "ASA",
        "septal myectomy alternative", "HCM intervention",
    ],
    "endomyocardial_biopsy": [
        "EMB", "endomyocardial biopsy", "RV biopsy",
        "transplant rejection biopsy",
    ],
    "pericardiocentesis": [
        "pericardiocentesis", "pericardial drainage",
        "pericardial window", "tamponade drainage",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# ENTITY ALIASES (abbreviation → canonical name)
# ═══════════════════════════════════════════════════════════════════════

ENTITY_ALIASES: Dict[str, str] = {
    # Acute coronary syndromes
    "MI": "myocardial infarction",
    "STEMI": "ST-elevation myocardial infarction",
    "NSTEMI": "non-ST-elevation myocardial infarction",
    "ACS": "acute coronary syndrome",
    "UA": "unstable angina",

    # Heart failure
    "HF": "heart failure",
    "CHF": "congestive heart failure",
    "HFrEF": "heart failure with reduced ejection fraction",
    "HFpEF": "heart failure with preserved ejection fraction",
    "HFmrEF": "heart failure with mildly reduced ejection fraction",
    "HFimpEF": "heart failure with improved ejection fraction",
    "LVEF": "left ventricular ejection fraction",
    "EF": "ejection fraction",
    "BNP": "brain natriuretic peptide",
    "NT-proBNP": "N-terminal pro-brain natriuretic peptide",
    "GDMT": "guideline-directed medical therapy",
    "ARNI": "angiotensin receptor-neprilysin inhibitor",
    "SGLT2i": "sodium-glucose cotransporter 2 inhibitor",
    "MRA": "mineralocorticoid receptor antagonist",
    "ACEi": "angiotensin-converting enzyme inhibitor",
    "ARB": "angiotensin receptor blocker",

    # Arrhythmia
    "AFib": "atrial fibrillation",
    "AF": "atrial fibrillation",
    "AFL": "atrial flutter",
    "VT": "ventricular tachycardia",
    "VF": "ventricular fibrillation",
    "SVT": "supraventricular tachycardia",
    "AVNRT": "atrioventricular nodal reentrant tachycardia",
    "AVRT": "atrioventricular reentrant tachycardia",
    "WPW": "Wolff-Parkinson-White syndrome",
    "LQTS": "long QT syndrome",
    "SQTS": "short QT syndrome",
    "CPVT": "catecholaminergic polymorphic ventricular tachycardia",
    "TdP": "torsades de pointes",
    "SCD": "sudden cardiac death",
    "NSVT": "non-sustained ventricular tachycardia",
    "PVC": "premature ventricular complex",
    "PAC": "premature atrial complex",

    # Devices
    "ICD": "implantable cardioverter defibrillator",
    "CRT": "cardiac resynchronization therapy",
    "CRT-D": "cardiac resynchronization therapy defibrillator",
    "CRT-P": "cardiac resynchronization therapy pacemaker",
    "PPM": "permanent pacemaker",
    "LVAD": "left ventricular assist device",
    "MCS": "mechanical circulatory support",
    "S-ICD": "subcutaneous implantable cardioverter defibrillator",
    "ATP": "antitachycardia pacing",
    "CSP": "conduction system pacing",
    "LBBAP": "left bundle branch area pacing",

    # Valvular
    "AS": "aortic stenosis",
    "AR": "aortic regurgitation",
    "MS": "mitral stenosis",
    "MR": "mitral regurgitation",
    "TR": "tricuspid regurgitation",
    "TS": "tricuspid stenosis",
    "PR": "pulmonic regurgitation",
    "PS": "pulmonic stenosis",
    "TAVR": "transcatheter aortic valve replacement",
    "TAVI": "transcatheter aortic valve implantation",
    "SAVR": "surgical aortic valve replacement",
    "AVR": "aortic valve replacement",
    "MVR": "mitral valve replacement",
    "TEER": "transcatheter edge-to-edge repair",
    "AVA": "aortic valve area",
    "EOA": "effective orifice area",
    "EROA": "effective regurgitant orifice area",
    "MVP": "mitral valve prolapse",
    "BAV": "bicuspid aortic valve",

    # Imaging
    "TTE": "transthoracic echocardiogram",
    "TEE": "transesophageal echocardiogram",
    "CMR": "cardiac magnetic resonance",
    "CTA": "coronary computed tomography angiography",
    "CAC": "coronary artery calcium",
    "LGE": "late gadolinium enhancement",
    "ECV": "extracellular volume fraction",
    "GLS": "global longitudinal strain",
    "TAPSE": "tricuspid annular plane systolic excursion",
    "RVSP": "right ventricular systolic pressure",
    "PASP": "pulmonary artery systolic pressure",
    "FFR-CT": "CT-derived fractional flow reserve",
    "IVUS": "intravascular ultrasound",
    "OCT": "optical coherence tomography",

    # Coronary / Hemodynamics
    "PCI": "percutaneous coronary intervention",
    "CABG": "coronary artery bypass graft",
    "DES": "drug-eluting stent",
    "BMS": "bare-metal stent",
    "FFR": "fractional flow reserve",
    "iFR": "instantaneous wave-free ratio",
    "DAPT": "dual antiplatelet therapy",
    "LAD": "left anterior descending artery",
    "LCx": "left circumflex artery",
    "RCA": "right coronary artery",
    "LVEDP": "left ventricular end-diastolic pressure",
    "PCWP": "pulmonary capillary wedge pressure",
    "CO": "cardiac output",
    "CI": "cardiac index",
    "PVR": "pulmonary vascular resistance",
    "SVR": "systemic vascular resistance",
    "SV": "stroke volume",
    "SVI": "stroke volume index",
    "RHC": "right heart catheterization",

    # ECG / EP
    "ECG": "electrocardiogram",
    "EKG": "electrocardiogram",
    "LBBB": "left bundle branch block",
    "RBBB": "right bundle branch block",
    "LAFB": "left anterior fascicular block",
    "LPFB": "left posterior fascicular block",
    "LVH": "left ventricular hypertrophy",
    "RVH": "right ventricular hypertrophy",
    "PFA": "pulsed field ablation",
    "PVI": "pulmonary vein isolation",

    # Congenital
    "ASD": "atrial septal defect",
    "VSD": "ventricular septal defect",
    "PDA": "patent ductus arteriosus",
    "PFO": "patent foramen ovale",
    "TOF": "Tetralogy of Fallot",
    "TGA": "transposition of the great arteries",
    "ACHD": "adult congenital heart disease",
    "CHD": "congenital heart disease",
    "HLHS": "hypoplastic left heart syndrome",
    "TCPC": "total cavopulmonary connection",
    "PAH": "pulmonary arterial hypertension",

    # Coronary / Specialty
    "SCAD": "spontaneous coronary artery dissection",
    "MINOCA": "myocardial infarction with non-obstructive coronary arteries",
    "INOCA": "ischemia with non-obstructive coronary arteries",

    # Vascular
    "AAA": "abdominal aortic aneurysm",
    "TAA": "thoracic aortic aneurysm",
    "TAAA": "thoracoabdominal aortic aneurysm",
    "EVAR": "endovascular aneurysm repair",
    "TEVAR": "thoracic endovascular aortic repair",
    "PAD": "peripheral artery disease",
    "PVD": "peripheral vascular disease",
    "ABI": "ankle-brachial index",
    "CEA": "carotid endarterectomy",
    "CAS": "carotid artery stenting",
    "CLI": "critical limb ischemia",
    "CLTI": "chronic limb-threatening ischemia",

    # Prevention / Lipids
    "ASCVD": "atherosclerotic cardiovascular disease",
    "PCE": "pooled cohort equations",
    "FH": "familial hypercholesterolemia",
    "HeFH": "heterozygous familial hypercholesterolemia",
    "HoFH": "homozygous familial hypercholesterolemia",
    "PCSK9i": "PCSK9 inhibitor",
    "LDL": "low-density lipoprotein",
    "HDL": "high-density lipoprotein",
    "TG": "triglycerides",
    "Lp(a)": "lipoprotein(a)",
    "ApoB": "apolipoprotein B",

    # Cardio-oncology
    "CTRCD": "cancer therapy-related cardiac dysfunction",
    "ICI": "immune checkpoint inhibitor",
    "TKI": "tyrosine kinase inhibitor",
    "RIHD": "radiation-induced heart disease",

    # Genomics
    "DCM": "dilated cardiomyopathy",
    "HCM": "hypertrophic cardiomyopathy",
    "RCM": "restrictive cardiomyopathy",
    "ACM": "arrhythmogenic cardiomyopathy",
    "ARVC": "arrhythmogenic right ventricular cardiomyopathy",
    "LVNC": "left ventricular non-compaction",
    "VUS": "variant of uncertain significance",
    "WES": "whole exome sequencing",
    "WGS": "whole genome sequencing",

    # Specialty conditions
    "CMD": "coronary microvascular disease",
    "PPCM": "peripartum cardiomyopathy",
    "ATTR": "transthyretin amyloidosis",
    "ATTR-CM": "ATTR cardiomyopathy",

    # Newer therapeutics
    "nsMRA": "non-steroidal mineralocorticoid receptor antagonist",

    # Newer procedures / modalities
    "ECPR": "extracorporeal CPR",
}


# ═══════════════════════════════════════════════════════════════════════
# COMPARATIVE QUERY PATTERNS
# ═══════════════════════════════════════════════════════════════════════

COMPARATIVE_PATTERNS: List[str] = [
    r"\bvs\.?\b",
    r"\bversus\b",
    r"\bcompared\s+to\b",
    r"\bcompared\s+with\b",
    r"\bcomparison\s+of\b",
    r"\bcomparison\s+between\b",
    r"\bbetter\s+than\b",
    r"\bsuperior\s+to\b",
    r"\binferior\s+to\b",
    r"\bnon-inferior\b",
    r"\bnoninferi",
    r"\bprefer(?:red|able)?\s+(?:to|over)\b",
    r"\bdifference\s+between\b",
    r"\bdifferences\s+between\b",
    r"\badvantage(?:s)?\s+(?:of|over)\b",
    r"\bdisadvantage(?:s)?\s+(?:of|over)\b",
    r"\bhead[- ]to[- ]head\b",
    r"\brandomized\s+(?:to|comparison)\b",
    r"\bswitch(?:ing)?\s+from\b",
    r"\btransition(?:ing)?\s+from\b",
]


# ═══════════════════════════════════════════════════════════════════════
# COMMON CARDIOLOGY COMPARISONS
# ═══════════════════════════════════════════════════════════════════════

CARDIO_COMPARISONS: Dict[str, Dict[str, List[str]]] = {
    "tavr vs savr": {
        "side_a": [
            "TAVR", "transcatheter aortic valve replacement",
            "transfemoral", "less invasive", "shorter recovery",
            "paravalvular leak", "pacemaker", "conduction disturbance",
            "PARTNER trial", "Evolut trial", "valve durability",
        ],
        "side_b": [
            "SAVR", "surgical aortic valve replacement", "sternotomy",
            "cardiopulmonary bypass", "open heart surgery",
            "surgical risk", "STS score", "mechanical valve option",
            "concomitant procedures", "longer durability",
        ],
        "shared": [
            "aortic stenosis", "severe AS", "aortic valve replacement",
            "heart team", "mortality", "stroke", "quality of life",
        ],
    },
    "doac vs warfarin": {
        "side_a": [
            "DOAC", "direct oral anticoagulant", "apixaban", "rivaroxaban",
            "edoxaban", "dabigatran", "fixed dose", "no INR monitoring",
            "fewer drug interactions", "fewer food interactions",
            "lower intracranial bleeding", "RE-LY", "ROCKET-AF",
            "ARISTOTLE", "ENGAGE AF-TIMI 48",
        ],
        "side_b": [
            "warfarin", "vitamin K antagonist", "VKA", "INR monitoring",
            "variable dosing", "bridging", "reversal with vitamin K",
            "mechanical valve", "moderate-severe mitral stenosis",
            "antiphospholipid syndrome", "cost", "generic",
        ],
        "shared": [
            "anticoagulation", "atrial fibrillation", "stroke prevention",
            "bleeding risk", "HAS-BLED", "CHA2DS2-VASc",
            "thromboembolism", "reversal agent",
        ],
    },
    "pci vs cabg": {
        "side_a": [
            "PCI", "percutaneous coronary intervention", "stent",
            "drug-eluting stent", "DES", "less invasive",
            "shorter recovery", "catheterization", "DAPT",
            "restenosis", "repeat revascularization",
        ],
        "side_b": [
            "CABG", "coronary artery bypass graft", "bypass surgery",
            "LIMA", "saphenous vein graft", "SVG", "sternotomy",
            "longer durability", "complete revascularization",
            "survival benefit", "diabetic patients",
        ],
        "shared": [
            "coronary artery disease", "CAD", "revascularization",
            "multivessel disease", "left main disease", "SYNTAX score",
            "heart team", "mortality", "MI", "stroke",
            "FREEDOM trial", "EXCEL trial", "NOBLE trial",
        ],
    },
    "arni vs acei": {
        "side_a": [
            "ARNI", "sacubitril-valsartan", "Entresto",
            "angiotensin receptor-neprilysin inhibitor",
            "neprilysin inhibition", "natriuretic peptide enhancement",
            "PARADIGM-HF", "superior mortality reduction",
            "reduced HF hospitalization",
        ],
        "side_b": [
            "ACEi", "ACE inhibitor", "angiotensin-converting enzyme inhibitor",
            "enalapril", "lisinopril", "ramipril",
            "established therapy", "lower cost", "generic",
            "cough", "angioedema",
        ],
        "shared": [
            "heart failure", "HFrEF", "RAAS inhibition",
            "GDMT", "mortality", "hospitalization",
            "hyperkalemia", "hypotension", "renal function",
            "washout period",
        ],
    },
    "ablation vs rate control": {
        "side_a": [
            "ablation", "catheter ablation", "pulmonary vein isolation",
            "PVI", "rhythm control", "sinus rhythm restoration",
            "CASTLE-AF", "CABANA", "quality of life",
            "symptom improvement", "reduced heart failure",
        ],
        "side_b": [
            "rate control", "beta blocker", "calcium channel blocker",
            "diltiazem", "verapamil", "metoprolol", "digoxin",
            "heart rate target", "less invasive", "fewer complications",
            "AFFIRM", "RACE", "AF-CHF",
        ],
        "shared": [
            "atrial fibrillation", "AFib", "AF",
            "anticoagulation", "CHA2DS2-VASc", "stroke risk",
            "symptom management", "quality of life", "heart failure",
        ],
    },
    "apixaban vs rivaroxaban": {
        "side_a": [
            "apixaban", "Eliquis", "twice daily", "ARISTOTLE trial",
            "lower bleeding rates", "lower dose available",
            "renal dosing", "dose adjustment criteria",
        ],
        "side_b": [
            "rivaroxaban", "Xarelto", "once daily", "ROCKET-AF trial",
            "with food", "single daily dose", "VTE treatment",
            "PAD indication", "COMPASS trial",
        ],
        "shared": [
            "DOAC", "direct oral anticoagulant", "factor Xa inhibitor",
            "anticoagulation", "atrial fibrillation", "DVT", "PE",
            "stroke prevention", "bleeding risk",
        ],
    },
    "stress echo vs nuclear": {
        "side_a": [
            "stress echocardiography", "exercise echo", "dobutamine echo",
            "wall motion abnormality", "no radiation", "lower cost",
            "valvular assessment", "diastolic function",
            "exercise capacity", "hemodynamics",
        ],
        "side_b": [
            "nuclear stress test", "myocardial perfusion imaging", "MPI",
            "SPECT", "PET", "technetium", "thallium", "rubidium",
            "pharmacologic stress", "adenosine", "regadenoson",
            "dipyridamole", "perfusion defect", "quantitative flow",
        ],
        "shared": [
            "stress testing", "ischemia", "CAD evaluation",
            "functional assessment", "risk stratification",
            "sensitivity", "specificity", "prognostic value",
        ],
    },
    "mechanical vs bioprosthetic valve": {
        "side_a": [
            "mechanical valve", "St. Jude", "On-X", "tilting disc",
            "bileaflet", "lifelong anticoagulation", "warfarin",
            "INR monitoring", "durability", "younger patients",
            "longer lifespan", "thrombosis risk",
        ],
        "side_b": [
            "bioprosthetic valve", "tissue valve", "porcine", "bovine",
            "pericardial", "Edwards", "Medtronic", "limited anticoagulation",
            "structural valve deterioration", "SVD", "reintervention",
            "valve-in-valve", "older patients",
        ],
        "shared": [
            "valve replacement", "prosthetic valve", "endocarditis",
            "paravalvular leak", "hemolysis", "patient-prosthesis mismatch",
            "heart team", "shared decision-making",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════
# WORKFLOW-SPECIFIC SEARCH TERMS
# ═══════════════════════════════════════════════════════════════════════

_WORKFLOW_TERMS: Dict[CardioWorkflowType, List[str]] = {
    CardioWorkflowType.CAD_ASSESSMENT: [
        "coronary artery disease", "CAD", "ischemia", "angina",
        "troponin", "stress test", "angiography", "PCI", "CABG",
        "statin", "antiplatelet", "revascularization", "calcium score",
        "CTA", "FFR", "risk stratification",
        "endocarditis", "pericarditis",
    ],
    CardioWorkflowType.HEART_FAILURE: [
        "heart failure", "HF", "LVEF", "BNP", "NT-proBNP",
        "GDMT", "ARNI", "SGLT2i", "beta blocker", "MRA",
        "diuretic", "NYHA", "decompensation", "congestion",
        "cardiomyopathy", "ICD", "CRT", "transplant", "LVAD",
        "pulmonary hypertension", "pericarditis",
    ],
    CardioWorkflowType.VALVULAR_DISEASE: [
        "valvular disease", "aortic stenosis", "mitral regurgitation",
        "TAVR", "valve replacement", "valve repair", "echo",
        "gradient", "valve area", "regurgitant volume",
        "endocarditis", "prosthetic valve", "anticoagulation",
    ],
    CardioWorkflowType.ARRHYTHMIA: [
        "arrhythmia", "atrial fibrillation", "ventricular tachycardia",
        "SVT", "bradycardia", "heart block", "ablation",
        "anticoagulation", "CHA2DS2-VASc", "ICD", "pacemaker",
        "ECG", "EP study", "antiarrhythmic", "cardioversion",
    ],
    CardioWorkflowType.CARDIAC_MRI: [
        "cardiac MRI", "CMR", "LGE", "T1 mapping", "T2 mapping",
        "ECV", "fibrosis", "edema", "perfusion", "strain",
        "tissue characterization", "cardiomyopathy", "myocarditis",
        "sarcoid", "amyloid", "iron overload",
        "pericarditis", "endocarditis", "pulmonary hypertension",
    ],
    CardioWorkflowType.STRESS_TEST: [
        "stress test", "exercise stress", "pharmacologic stress",
        "dobutamine", "adenosine", "regadenoson",
        "stress echocardiography", "nuclear perfusion", "SPECT",
        "PET", "ischemia", "wall motion", "perfusion defect",
        "exercise capacity", "METs", "Duke treadmill score",
    ],
    CardioWorkflowType.PREVENTIVE_RISK: [
        "prevention", "ASCVD risk", "statin", "LDL", "cholesterol",
        "lipid", "lifestyle", "exercise", "diet", "blood pressure",
        "smoking cessation", "calcium score", "risk assessment",
        "Lp(a)", "PCSK9i", "risk factors",
    ],
    CardioWorkflowType.CARDIO_ONCOLOGY: [
        "cardio-oncology", "cardiotoxicity", "CTRCD",
        "anthracycline", "trastuzumab", "checkpoint inhibitor",
        "GLS", "LVEF monitoring", "dexrazoxane", "cardioprotection",
        "surveillance", "radiation heart disease", "cancer therapy",
    ],
    CardioWorkflowType.GENERAL: [
        "cardiology", "cardiac", "cardiovascular", "heart",
        "general cardiology", "consultation",
        "endocarditis", "pericarditis", "pulmonary hypertension",
    ],
    CardioWorkflowType.ACUTE_DECOMPENSATED_HF: [
        "acute heart failure", "decompensated", "ADHF", "acute HF",
        "IV diuretics", "furosemide IV", "inotropes", "dobutamine",
        "milrinone", "cardiogenic shock", "Impella", "IABP",
        "warm wet", "cold wet", "congestion", "volume overload",
        "decongestion", "diuretic resistance", "EMPULSE",
    ],
    CardioWorkflowType.POST_MI: [
        "myocardial infarction", "MI", "STEMI", "NSTEMI",
        "post-MI", "ACS", "reperfusion", "PCI", "thrombolysis",
        "DAPT", "ticagrelor", "clopidogrel", "aspirin",
        "cardiac rehabilitation", "secondary prevention",
        "door-to-balloon", "troponin", "culprit vessel",
    ],
    CardioWorkflowType.MYOCARDITIS_PERICARDITIS: [
        "myocarditis", "pericarditis", "myopericarditis",
        "Lake Louise", "troponin elevation", "pericardial effusion",
        "colchicine", "NSAID", "ibuprofen", "rilonacept",
        "giant cell myocarditis", "eosinophilic",
        "endomyocardial biopsy", "activity restriction",
        "return to play", "recurrent pericarditis",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# ENTITY CATEGORIES (for detect_entities)
# ═══════════════════════════════════════════════════════════════════════

_ENTITY_CATEGORIES: Dict[str, List[str]] = {
    "conditions": [
        "heart failure", "coronary artery disease", "myocardial infarction",
        "atrial fibrillation", "ventricular tachycardia", "aortic stenosis",
        "mitral regurgitation", "cardiomyopathy", "hypertrophic cardiomyopathy",
        "dilated cardiomyopathy", "myocarditis", "pericarditis",
        "endocarditis", "pulmonary hypertension", "aortic dissection",
        "aortic aneurysm", "peripheral artery disease", "deep vein thrombosis",
        "pulmonary embolism", "stroke", "heart block",
        "long QT syndrome", "Brugada syndrome", "Marfan syndrome",
        "familial hypercholesterolemia", "Eisenmenger syndrome",
        "Tetralogy of Fallot", "amyloidosis", "sarcoidosis",
        "cardiotoxicity", "congenital heart disease",
    ],
    "drugs": [
        "sacubitril", "valsartan", "entresto", "carvedilol", "metoprolol",
        "bisoprolol", "spironolactone", "eplerenone", "dapagliflozin",
        "empagliflozin", "furosemide", "bumetanide", "torsemide",
        "enalapril", "lisinopril", "ramipril", "losartan", "candesartan",
        "amiodarone", "flecainide", "dronedarone", "sotalol", "dofetilide",
        "apixaban", "rivaroxaban", "dabigatran", "edoxaban", "warfarin",
        "aspirin", "clopidogrel", "ticagrelor", "prasugrel",
        "atorvastatin", "rosuvastatin", "simvastatin", "pravastatin",
        "ezetimibe", "evolocumab", "alirocumab", "inclisiran",
        "bempedoic acid", "icosapent ethyl",
        "dobutamine", "milrinone", "dopamine", "norepinephrine",
        "adenosine", "verapamil", "diltiazem", "digoxin",
        "nitroglycerin", "nitroprusside", "hydralazine",
        "doxorubicin", "trastuzumab", "pembrolizumab", "nivolumab",
        "dexrazoxane", "heparin", "enoxaparin", "bivalirudin",
    ],
    "genes": [
        "MYH7", "MYBPC3", "TTN", "LMNA", "TNNT2", "TNNI3",
        "SCN5A", "KCNQ1", "KCNH2", "KCNJ2", "CACNA1C",
        "RYR2", "CASQ2", "DSP", "PKP2", "DSG2", "DSC2",
        "FBN1", "TGFBR1", "TGFBR2", "LDLR", "PCSK9", "APOB",
        "LPA", "ACTC1", "MYL2", "MYL3", "TPM1", "FLNC",
        "PLN", "BAG3", "RBM20", "ANKRD1",
    ],
    "imaging_modalities": [
        "echocardiography", "echocardiogram", "echo",
        "transthoracic echocardiogram", "TTE",
        "transesophageal echocardiogram", "TEE",
        "cardiac MRI", "CMR", "cardiac magnetic resonance",
        "cardiac CT", "coronary CT angiography", "CTA",
        "nuclear perfusion", "SPECT", "PET",
        "angiography", "coronary angiography",
        "intravascular ultrasound", "IVUS",
        "optical coherence tomography", "OCT",
    ],
    "procedures": [
        "PCI", "percutaneous coronary intervention",
        "CABG", "coronary artery bypass graft",
        "TAVR", "transcatheter aortic valve replacement",
        "SAVR", "surgical aortic valve replacement",
        "ablation", "catheter ablation", "cardioversion",
        "pacemaker implant", "ICD implant", "CRT implant",
        "LVAD implant", "heart transplant",
        "valve repair", "valve replacement",
        "lead extraction", "PFO closure", "ASD closure",
        "EVAR", "TEVAR", "carotid endarterectomy",
    ],
    "biomarkers": [
        "troponin", "troponin I", "troponin T", "high-sensitivity troponin",
        "BNP", "NT-proBNP", "CK-MB", "D-dimer",
        "CRP", "hs-CRP", "Lp(a)", "ApoB", "LDL", "HDL",
        "triglycerides", "HbA1c", "creatinine", "GFR", "eGFR",
        "INR", "PTT", "anti-Xa level",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION BOOST RULES
# ═══════════════════════════════════════════════════════════════════════

_COLLECTION_BOOST_RULES: Dict[str, Dict[str, float]] = {
    "heart_failure": {
        "clinical_guidelines": 1.5,
        "drug_interactions": 1.3,
        "biomarker_reference": 1.2,
        "genomic_variants": 1.1,
    },
    "coronary": {
        "clinical_guidelines": 1.4,
        "imaging_protocols": 1.3,
        "drug_interactions": 1.2,
        "procedural_reference": 1.5,
    },
    "arrhythmia": {
        "clinical_guidelines": 1.4,
        "drug_interactions": 1.4,
        "device_reference": 1.5,
        "electrophysiology": 1.5,
    },
    "valvular": {
        "clinical_guidelines": 1.4,
        "imaging_protocols": 1.5,
        "procedural_reference": 1.5,
        "hemodynamic_reference": 1.3,
    },
    "imaging": {
        "imaging_protocols": 1.5,
        "clinical_guidelines": 1.2,
        "hemodynamic_reference": 1.2,
    },
    "prevention": {
        "clinical_guidelines": 1.5,
        "drug_interactions": 1.3,
        "biomarker_reference": 1.3,
        "risk_calculators": 1.5,
    },
    "hemodynamics": {
        "hemodynamic_reference": 1.5,
        "procedural_reference": 1.4,
        "clinical_guidelines": 1.2,
    },
    "electrophysiology": {
        "electrophysiology": 1.5,
        "device_reference": 1.4,
        "drug_interactions": 1.3,
        "clinical_guidelines": 1.3,
    },
    "cardio_oncology": {
        "clinical_guidelines": 1.4,
        "drug_interactions": 1.5,
        "biomarker_reference": 1.3,
        "imaging_protocols": 1.3,
    },
    "congenital": {
        "clinical_guidelines": 1.4,
        "imaging_protocols": 1.3,
        "procedural_reference": 1.3,
        "hemodynamic_reference": 1.2,
    },
    "device": {
        "device_reference": 1.5,
        "clinical_guidelines": 1.3,
        "electrophysiology": 1.3,
        "procedural_reference": 1.2,
    },
    "vascular": {
        "clinical_guidelines": 1.3,
        "imaging_protocols": 1.4,
        "procedural_reference": 1.5,
        "hemodynamic_reference": 1.2,
    },
    "genomics": {
        "genomic_variants": 1.5,
        "clinical_guidelines": 1.3,
        "biomarker_reference": 1.2,
        "family_screening": 1.4,
    },
}


# ═══════════════════════════════════════════════════════════════════════
# QueryExpander CLASS
# ═══════════════════════════════════════════════════════════════════════


class QueryExpander:
    """Expand cardiology queries with related clinical terms.

    Orchestrates 15 domain-specific expansion maps, entity alias
    resolution, comparative-query detection, and collection-weight
    boosting to maximise recall across the 12 Milvus collections.
    """

    def __init__(self) -> None:
        self.expansion_maps: List[Dict[str, List[str]]] = [
            HEART_FAILURE_MAP,
            CORONARY_ARTERY_MAP,
            ARRHYTHMIA_MAP,
            VALVULAR_MAP,
            IMAGING_ECHO_MAP,
            IMAGING_CT_MAP,
            IMAGING_MRI_MAP,
            PREVENTIVE_MAP,
            HEMODYNAMICS_MAP,
            ELECTROPHYSIOLOGY_MAP,
            CARDIO_ONCOLOGY_MAP,
            CONGENITAL_MAP,
            DEVICE_MAP,
            VASCULAR_MAP,
            GENOMICS_CARDIO_MAP,
            SPECIALTY_SYNONYM_MAP,
            RARE_DISEASE_MAP,
            PROCEDURE_MAP,
        ]
        self.entity_aliases: Dict[str, str] = ENTITY_ALIASES
        self._comparative_re = re.compile(
            "|".join(f"({p})" for p in COMPARATIVE_PATTERNS),
            re.IGNORECASE,
        )

    # ───────────────────────────────────────────────────────────────
    # Public API
    # ───────────────────────────────────────────────────────────────

    def expand(
        self,
        query: str,
        workflow: Optional[CardioWorkflowType] = None,
    ) -> dict:
        """Expand a query with related cardiovascular terms.

        Parameters
        ----------
        query : str
            The raw user query.
        workflow : CardioWorkflowType, optional
            If provided, additional workflow-specific terms are appended.

        Returns
        -------
        dict
            ``original``            – the raw query text
            ``expanded_terms``      – list of additional search terms
            ``detected_entities``   – dict of categorised entities found
            ``is_comparative``      – whether this is a comparison query
            ``comparison``          – parsed comparison details (or None)
            ``workflow_hint``       – workflow type used for expansion
        """
        query_lower = query.lower().strip()

        # Collect expanded terms from all maps
        expanded_terms: List[str] = []
        for exp_map in self.expansion_maps:
            for trigger, terms in exp_map.items():
                if self._trigger_matches(trigger, query_lower):
                    for term in terms:
                        if term.lower() not in query_lower and term not in expanded_terms:
                            expanded_terms.append(term)

        # Add workflow-specific terms
        workflow_terms: List[str] = []
        if workflow is not None:
            workflow_terms = self.get_workflow_terms(workflow)
            for term in workflow_terms:
                if term.lower() not in query_lower and term not in expanded_terms:
                    expanded_terms.append(term)

        # Detect entities
        detected_entities = self.detect_entities(query)

        # Detect comparative queries
        comparison = self.detect_comparative(query)
        is_comparative = comparison is not None

        return {
            "original": query,
            "expanded_terms": expanded_terms,
            "detected_entities": detected_entities,
            "is_comparative": is_comparative,
            "comparison": comparison,
            "workflow_hint": workflow.value if workflow else None,
        }

    def resolve_aliases(self, text: str) -> str:
        """Replace known abbreviations with their canonical names.

        Parameters
        ----------
        text : str
            Input text that may contain cardiology abbreviations.

        Returns
        -------
        str
            Text with abbreviations expanded to full canonical names.
        """
        result = text
        for abbr, canonical in sorted(
            self.entity_aliases.items(), key=lambda x: -len(x[0])
        ):
            # Only replace whole-word matches (case-sensitive for
            # abbreviations since they are typically uppercase)
            pattern = re.compile(r"\b" + re.escape(abbr) + r"\b")
            if pattern.search(result):
                result = pattern.sub(f"{abbr} ({canonical})", result)
        return result

    def detect_entities(self, text: str) -> dict:
        """Detect cardiology entities in text by category.

        Parameters
        ----------
        text : str
            Input text to scan.

        Returns
        -------
        dict
            Keys are category names (``conditions``, ``drugs``,
            ``genes``, ``imaging_modalities``, ``procedures``,
            ``biomarkers``); values are lists of matched entity strings.
        """
        text.lower()
        detected: Dict[str, List[str]] = {}

        for category, entities in _ENTITY_CATEGORIES.items():
            matches: List[str] = []
            for entity in entities:
                # Check for whole-word match (case-insensitive)
                pattern = re.compile(
                    r"\b" + re.escape(entity) + r"\b", re.IGNORECASE
                )
                if pattern.search(text):
                    matches.append(entity)
                # Also check if an alias resolves to this entity
                elif entity.upper() in self.entity_aliases:
                    alias_pattern = re.compile(
                        r"\b" + re.escape(entity.upper()) + r"\b"
                    )
                    if alias_pattern.search(text):
                        matches.append(entity)

            if matches:
                detected[category] = matches

        # Also detect aliases that appear directly
        alias_matches: List[str] = []
        for abbr in self.entity_aliases:
            pattern = re.compile(r"\b" + re.escape(abbr) + r"\b")
            if pattern.search(text):
                canonical = self.entity_aliases[abbr]
                entry = f"{abbr} ({canonical})"
                if entry not in alias_matches:
                    alias_matches.append(entry)

        if alias_matches:
            detected["resolved_aliases"] = alias_matches

        return detected

    def detect_comparative(self, text: str) -> Optional[dict]:
        """Detect and parse comparative queries.

        Parameters
        ----------
        text : str
            Query text to analyse.

        Returns
        -------
        dict or None
            If a comparative pattern is found, returns a dict with:
            ``pattern_matched`` – the regex pattern that triggered
            ``known_comparison`` – matching entry from CARDIO_COMPARISONS
                                  (or None if not a known pair)
            ``raw_sides``       – tuple of (side_a_text, side_b_text) if
                                  a ``vs``-style split was possible
        """
        match = self._comparative_re.search(text)
        if match is None:
            return None

        result: dict = {
            "pattern_matched": match.group(0),
            "known_comparison": None,
            "raw_sides": None,
        }

        # Try to identify a known comparison
        text_lower = text.lower()
        for comp_key, comp_data in CARDIO_COMPARISONS.items():
            # Build variants of the comparison key
            parts = comp_key.split(" vs ")
            if len(parts) == 2:
                a, b = parts
                if (a in text_lower and b in text_lower) or (
                    b in text_lower and a in text_lower
                ):
                    result["known_comparison"] = {
                        "key": comp_key,
                        **comp_data,
                    }
                    break

        # Try to extract raw sides from a "vs"-style split
        vs_pattern = re.compile(
            r"(.+?)\s+(?:vs\.?|versus|compared\s+to|compared\s+with)\s+(.+)",
            re.IGNORECASE,
        )
        vs_match = vs_pattern.search(text)
        if vs_match:
            result["raw_sides"] = (
                vs_match.group(1).strip(),
                vs_match.group(2).strip(),
            )

        return result

    def get_workflow_terms(self, workflow: CardioWorkflowType) -> List[str]:
        """Return additional search terms for a specific workflow.

        Parameters
        ----------
        workflow : CardioWorkflowType
            The clinical workflow.

        Returns
        -------
        list of str
            Terms relevant to the workflow.
        """
        return list(_WORKFLOW_TERMS.get(workflow, []))

    def boost_collections(self, expanded: dict) -> Dict[str, float]:
        """Suggest Milvus collection weight boosts based on expanded query.

        Analyses the expanded terms and detected entities to determine
        which domain(s) the query falls into, then returns a mapping
        of collection names to multiplicative boost factors.

        Parameters
        ----------
        expanded : dict
            Output from :meth:`expand`.

        Returns
        -------
        dict
            Collection name → boost factor (float >= 1.0).
        """
        boosts: Dict[str, float] = {}
        all_terms = " ".join(expanded.get("expanded_terms", [])).lower()
        original = expanded.get("original", "").lower()
        combined = f"{original} {all_terms}"

        # Determine which domain rules apply
        domain_scores: Dict[str, float] = {}

        domain_keywords: Dict[str, List[str]] = {
            "heart_failure": [
                "heart failure", "hf", "chf", "lvef", "bnp", "gdmt",
                "arni", "sglt2i", "cardiomyopathy", "decompensation",
            ],
            "coronary": [
                "coronary", "cad", "mi", "myocardial infarction", "acs",
                "stemi", "nstemi", "pci", "cabg", "angina", "troponin",
            ],
            "arrhythmia": [
                "arrhythmia", "atrial fibrillation", "afib", "vt", "svt",
                "bradycardia", "tachycardia", "ablation", "antiarrhythmic",
            ],
            "valvular": [
                "aortic stenosis", "mitral regurgitation", "tavr", "valve",
                "gradient", "regurgitant", "vena contracta",
            ],
            "imaging": [
                "echo", "echocardiogram", "cardiac ct", "cta", "cardiac mri",
                "cmr", "strain", "lge", "perfusion", "calcium score",
            ],
            "prevention": [
                "prevention", "ascvd", "statin", "cholesterol", "lipid",
                "ldl", "hdl", "pcsk9", "lifestyle", "risk",
            ],
            "hemodynamics": [
                "catheterization", "cath", "hemodynamics", "pressure",
                "lvedp", "pcwp", "cardiac output", "ffr",
            ],
            "electrophysiology": [
                "ecg", "ekg", "ep study", "electrophysiology", "ablation",
                "qrs", "qtc", "lbbb", "rbbb", "mapping",
            ],
            "cardio_oncology": [
                "cardiotoxicity", "cardio-oncology", "anthracycline",
                "trastuzumab", "checkpoint inhibitor", "ctrcd", "gls",
            ],
            "congenital": [
                "congenital", "chd", "asd", "vsd", "pfo", "tetralogy",
                "fontan", "eisenmenger", "shunt",
            ],
            "device": [
                "pacemaker", "icd", "crt", "lvad", "lead", "interrogation",
                "device", "shock", "pacing",
            ],
            "vascular": [
                "aorta", "pad", "peripheral", "carotid", "aneurysm",
                "dissection", "endovascular", "stent-graft",
            ],
            "genomics": [
                "cardiomyopathy", "channelopathy", "marfan", "genetic",
                "myh7", "mybpc3", "ttn", "scn5a", "kcnq1", "ldlr",
                "gene panel", "variant", "pathogenic",
            ],
        }

        for domain, keywords in domain_keywords.items():
            score = sum(1.0 for kw in keywords if kw in combined)
            if score > 0:
                domain_scores[domain] = score

        if not domain_scores:
            return boosts

        # Apply boost rules for all matched domains, weighted by score
        max_score = max(domain_scores.values())
        for domain, score in domain_scores.items():
            weight = score / max_score  # normalise to [0, 1]
            if domain in _COLLECTION_BOOST_RULES:
                for collection, boost_val in _COLLECTION_BOOST_RULES[domain].items():
                    # Scale the boost by domain relevance
                    scaled_boost = 1.0 + (boost_val - 1.0) * weight
                    if collection in boosts:
                        boosts[collection] = max(boosts[collection], scaled_boost)
                    else:
                        boosts[collection] = scaled_boost

        # Round for readability
        boosts = {k: round(v, 2) for k, v in boosts.items() if v > 1.0}

        return boosts

    # ───────────────────────────────────────────────────────────────
    # Private helpers
    # ───────────────────────────────────────────────────────────────

    @staticmethod
    def _trigger_matches(trigger: str, query_lower: str) -> bool:
        """Check whether *trigger* appears as a whole word in *query_lower*."""
        pattern = re.compile(r"\b" + re.escape(trigger) + r"\b", re.IGNORECASE)
        return bool(pattern.search(query_lower))
