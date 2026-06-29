"""ACC/AHA/ESC guideline ingest parser for Cardiology Intelligence Agent.

Provides a curated registry of major cardiovascular clinical practice
guidelines from the American College of Cardiology (ACC), American Heart
Association (AHA), European Society of Cardiology (ESC), Heart Rhythm
Society (HRS), and Society for Cardiovascular Angiography and
Interventions (SCAI).

Each guideline entry includes structured metadata (society, year, DOI,
condition, key sections) and representative recommendation text suitable
for vector embedding and RAG retrieval.

Targets the ``cardio_guidelines`` Milvus collection.

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# GUIDELINE PARSER
# ═══════════════════════════════════════════════════════════════════════


class GuidelineParser(BaseIngestParser):
    """Ingest parser for ACC/AHA/ESC clinical practice guidelines.

    Maintains a curated dictionary of all major cardiovascular guidelines
    with their DOIs, publication years, issuing societies, and key
    recommendation sections.  Generates IngestRecord instances for the
    ``cardio_guidelines`` collection.

    Guidelines are seeded statically because full-text PDFs are
    copyright-protected.  The seed data provides enough context for
    RAG to surface guideline-concordant recommendations.

    Usage::

        parser = GuidelineParser()
        records = parser.seed_guidelines()

        # Or via the standard run() interface (fetch returns the
        # static guideline list):
        records = parser.run()
    """

    GUIDELINES: Dict[str, Dict[str, Any]] = {
        # ── Heart Failure ─────────────────────────────────────────
        "ACC/AHA 2022 Heart Failure": {
            "year": 2022,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000001063",
            "condition": "Heart Failure",
            "sections": [
                {
                    "section": "GDMT for HFrEF",
                    "recommendation": (
                        "In patients with HFrEF (LVEF <=40%), the use of "
                        "ARNI (or ACEi/ARB), evidence-based beta-blocker, "
                        "MRA, and SGLT2i is recommended to reduce morbidity "
                        "and mortality."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "SGLT2i in Heart Failure",
                    "recommendation": (
                        "In patients with HFrEF or HFpEF, SGLT2 inhibitors "
                        "are recommended to reduce HF hospitalizations and "
                        "cardiovascular mortality."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Device Therapy for HFrEF",
                    "recommendation": (
                        "ICD is recommended for primary prevention of SCD in "
                        "patients with LVEF <=35% despite >=3 months of GDMT, "
                        "with NYHA class II-III symptoms and life expectancy >1 year."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "CRT in Heart Failure",
                    "recommendation": (
                        "CRT is recommended in HFrEF patients with LVEF <=35%, "
                        "sinus rhythm, LBBB with QRS >=150 ms, and NYHA II-IV "
                        "symptoms despite GDMT."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Advanced Heart Failure",
                    "recommendation": (
                        "Referral to an advanced HF specialist is recommended "
                        "for patients with Stage D HF to evaluate candidacy for "
                        "MCS, heart transplant, or palliative care."
                    ),
                    "class": "I",
                    "evidence": "C-EO",
                },
            ],
        },
        # ── Valvular Heart Disease ────────────────────────────────
        "ACC/AHA 2020 Valvular Heart Disease": {
            "year": 2020,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000000923",
            "condition": "Valvular Heart Disease",
            "sections": [
                {
                    "section": "Severe Aortic Stenosis — Intervention",
                    "recommendation": (
                        "AVR (SAVR or TAVR) is recommended in symptomatic "
                        "patients with severe high-gradient aortic stenosis "
                        "(mean gradient >=40 mmHg or Vmax >=4.0 m/s)."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "TAVR in Low Surgical Risk",
                    "recommendation": (
                        "In patients >=65 years with severe symptomatic AS "
                        "and low surgical risk, either SAVR or TAVR is "
                        "recommended based on shared decision-making."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Mitral Regurgitation — TEER",
                    "recommendation": (
                        "Transcatheter edge-to-edge repair (MitraClip) is "
                        "recommended for severely symptomatic patients with "
                        "chronic severe secondary MR who remain symptomatic "
                        "despite GDMT and meet anatomic criteria."
                    ),
                    "class": "IIa",
                    "evidence": "B-R",
                },
            ],
        },
        # ── Prevention ────────────────────────────────────────────
        "ACC/AHA 2019 Prevention": {
            "year": 2019,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000000625",
            "condition": "Cardiovascular Prevention",
            "sections": [
                {
                    "section": "Statin Therapy — Primary Prevention",
                    "recommendation": (
                        "In adults 40-75 years with LDL-C >=70 mg/dL and "
                        "10-year ASCVD risk >=7.5%, moderate-intensity statin "
                        "therapy is recommended."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "ASCVD Risk Assessment",
                    "recommendation": (
                        "For adults 40-75 years, the 10-year ASCVD risk using "
                        "the Pooled Cohort Equations should be estimated to "
                        "guide primary prevention decisions."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Coronary Artery Calcium Score",
                    "recommendation": (
                        "When risk-based decisions are uncertain, CAC scoring "
                        "is reasonable to guide statin therapy decisions, "
                        "particularly for borderline-risk patients."
                    ),
                    "class": "IIa",
                    "evidence": "B-NR",
                },
            ],
        },
        # ── Atrial Fibrillation ───────────────────────────────────
        "ACC/AHA/HRS 2023 Atrial Fibrillation": {
            "year": 2023,
            "society": "ACC/AHA/HRS",
            "doi": "10.1161/CIR.0000000000001193",
            "condition": "Atrial Fibrillation",
            "sections": [
                {
                    "section": "Anticoagulation in AF",
                    "recommendation": (
                        "In patients with AF and CHA2DS2-VASc score >=2 in men "
                        "or >=3 in women, oral anticoagulation is recommended "
                        "to reduce stroke risk. DOACs are preferred over warfarin."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Catheter Ablation for AF",
                    "recommendation": (
                        "Catheter ablation of AF is recommended as first-line "
                        "rhythm control in select symptomatic patients with "
                        "paroxysmal AF."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Rate Control in AF",
                    "recommendation": (
                        "A lenient rate control strategy (resting HR <110 bpm) "
                        "is reasonable in patients with AF without significant "
                        "symptoms."
                    ),
                    "class": "IIa",
                    "evidence": "B-R",
                },
            ],
        },
        # ── Coronary Artery Disease ───────────────────────────────
        "ACC/AHA 2021 Chest Pain": {
            "year": 2021,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000001029",
            "condition": "Chest Pain / Stable CAD",
            "sections": [
                {
                    "section": "High-Sensitivity Troponin",
                    "recommendation": (
                        "High-sensitivity cardiac troponin I or T should be "
                        "the preferred biomarker for evaluating patients with "
                        "suspected acute coronary syndromes."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "CCTA for Chest Pain",
                    "recommendation": (
                        "Coronary CTA is an effective first-line test for "
                        "evaluation of stable chest pain in intermediate-risk "
                        "patients with no known CAD."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── ESC Heart Failure ─────────────────────────────────────
        "ESC 2021 Heart Failure": {
            "year": 2021,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehab368",
            "condition": "Heart Failure",
            "sections": [
                {
                    "section": "Four Pillars of GDMT",
                    "recommendation": (
                        "ACEi/ARNI, beta-blocker, MRA, and SGLT2i are "
                        "recommended as first-line therapies in all patients "
                        "with HFrEF to reduce hospitalization and death."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Iron Deficiency in HF",
                    "recommendation": (
                        "IV iron supplementation with ferric carboxymaltose "
                        "should be considered in symptomatic HFrEF patients "
                        "with iron deficiency to improve symptoms and QoL."
                    ),
                    "class": "IIa",
                    "evidence": "A",
                },
            ],
        },
        # ── ESC Atrial Fibrillation ───────────────────────────────
        "ESC 2024 Atrial Fibrillation": {
            "year": 2024,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehae176",
            "condition": "Atrial Fibrillation",
            "sections": [
                {
                    "section": "AF-CARE Pathway",
                    "recommendation": (
                        "The AF-CARE (Comorbidity/risk, Avoid stroke, Reduce "
                        "symptoms, Evaluate) pathway is recommended for "
                        "holistic, structured management of AF patients."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── Hypertrophic Cardiomyopathy ───────────────────────────
        "ACC/AHA 2020 Hypertrophic Cardiomyopathy": {
            "year": 2020,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000000937",
            "condition": "Hypertrophic Cardiomyopathy",
            "sections": [
                {
                    "section": "SCD Risk Stratification in HCM",
                    "recommendation": (
                        "ICD implantation is recommended in HCM patients with "
                        "prior cardiac arrest or sustained ventricular "
                        "tachycardia."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Septal Reduction in HCM",
                    "recommendation": (
                        "Septal reduction therapy (surgical myectomy preferred) "
                        "is recommended for drug-refractory symptomatic "
                        "obstructive HCM with LVOT gradient >=50 mmHg."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
            ],
        },
        # ── Pulmonary Hypertension ────────────────────────────────
        "ESC/ERS 2022 Pulmonary Hypertension": {
            "year": 2022,
            "society": "ESC/ERS",
            "doi": "10.1093/eurheartj/ehac237",
            "condition": "Pulmonary Hypertension",
            "sections": [
                {
                    "section": "Hemodynamic Definition of PH",
                    "recommendation": (
                        "Pulmonary hypertension is defined as mean pulmonary "
                        "arterial pressure >20 mmHg at rest by right heart "
                        "catheterization (revised from >25 mmHg threshold)."
                    ),
                    "class": "I",
                    "evidence": "C-LD",
                },
                {
                    "section": "Initial Combination Therapy in PAH",
                    "recommendation": (
                        "Initial combination therapy with an ERA and PDE5i "
                        "is recommended in treatment-naive PAH patients at "
                        "low or intermediate risk."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── Cardiac Implantable Devices ───────────────────────────
        "HRS 2018 Bradycardia / CRT": {
            "year": 2018,
            "society": "ACC/AHA/HRS",
            "doi": "10.1016/j.hrthm.2018.10.037",
            "condition": "Bradycardia and Cardiac Resynchronization",
            "sections": [
                {
                    "section": "Permanent Pacemaker — Sinus Node Dysfunction",
                    "recommendation": (
                        "Permanent pacing is recommended in patients with "
                        "symptomatic sinus node dysfunction including "
                        "documented symptomatic bradycardia or sinus pauses."
                    ),
                    "class": "I",
                    "evidence": "B-R",
                },
            ],
        },
        # ── Ventricular Arrhythmias / SCD ─────────────────────────
        "ACC/AHA/HRS 2017 Ventricular Arrhythmias": {
            "year": 2017,
            "society": "ACC/AHA/HRS",
            "doi": "10.1016/j.hrthm.2017.10.036",
            "condition": "Ventricular Arrhythmias and SCD Prevention",
            "sections": [
                {
                    "section": "ICD for Primary Prevention",
                    "recommendation": (
                        "ICD implantation is recommended for primary "
                        "prevention of SCD in patients with LVEF <=35% due "
                        "to prior MI (>=40 days) or non-ischemic "
                        "cardiomyopathy with NYHA II-III despite GDMT."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── Peripheral Artery Disease ─────────────────────────────
        "ACC/AHA 2016 PAD": {
            "year": 2016,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000000471",
            "condition": "Peripheral Artery Disease",
            "sections": [
                {
                    "section": "ABI Screening",
                    "recommendation": (
                        "Resting ankle-brachial index (ABI) is recommended "
                        "as a first-line noninvasive test for diagnosis of "
                        "PAD in patients with suspected lower extremity PAD."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
            ],
        },
        # ── Chronic Coronary Disease ─────────────────────────────
        "ACC/AHA 2023 Chronic Coronary Disease": {
            "year": 2023,
            "society": "ACC/AHA",
            "doi": "10.1161/CIR.0000000000001168",
            "condition": "Chronic Coronary Disease",
            "sections": [
                {
                    "section": "DAPT Duration Post-PCI",
                    "recommendation": (
                        "In patients with chronic coronary disease who have "
                        "undergone PCI, DAPT with aspirin and a P2Y12 inhibitor "
                        "should be given for 6 months. Shorter DAPT (1-3 months) "
                        "may be reasonable in patients at high bleeding risk."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Ischemia-Guided Revascularization",
                    "recommendation": (
                        "In patients with chronic coronary disease, an ischemia-guided "
                        "strategy using stress testing or FFR/iFR is recommended to "
                        "determine whether revascularization provides benefit beyond "
                        "optimal medical therapy alone."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Cardiac Rehabilitation",
                    "recommendation": (
                        "Cardiac rehabilitation is recommended for all eligible "
                        "patients with chronic coronary disease, including post-MI, "
                        "post-PCI, and post-CABG, to reduce cardiovascular mortality "
                        "and improve functional capacity."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "SGLT2i for CKD with Diabetes and CAD",
                    "recommendation": (
                        "In patients with chronic coronary disease and comorbid "
                        "type 2 diabetes and chronic kidney disease (eGFR 20-60), "
                        "SGLT2 inhibitors are recommended to reduce cardiovascular "
                        "events and kidney disease progression."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Long-Term Antithrombotic Therapy",
                    "recommendation": (
                        "Low-dose aspirin (75-100 mg daily) is recommended for "
                        "long-term secondary prevention in patients with chronic "
                        "coronary disease. Addition of low-dose rivaroxaban 2.5 mg "
                        "BID may be considered in high-risk patients without high "
                        "bleeding risk."
                    ),
                    "class": "IIa",
                    "evidence": "A",
                },
            ],
        },
        # ── ACS / NSTE-ACS ──────────────────────────────────────
        "ESC 2023 ACS/NSTE-ACS": {
            "year": 2023,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehad191",
            "condition": "Acute Coronary Syndromes / NSTE-ACS",
            "sections": [
                {
                    "section": "Very High-Risk Criteria for Immediate Catheterization",
                    "recommendation": (
                        "Immediate invasive strategy (<2 hours) is recommended in "
                        "NSTE-ACS patients with very high-risk criteria: "
                        "hemodynamic instability or cardiogenic shock, recurrent "
                        "or refractory chest pain, life-threatening arrhythmias, "
                        "mechanical complications, or acute heart failure clearly "
                        "related to ACS."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "High-Sensitivity Troponin 0/1h Algorithm",
                    "recommendation": (
                        "The 0h/1h high-sensitivity cardiac troponin algorithm is "
                        "recommended as the preferred rapid rule-in/rule-out "
                        "strategy for NSTE-ACS, with rule-out if baseline hs-cTn "
                        "is very low and 1h delta is below the assay-specific "
                        "threshold."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Early Invasive vs Conservative Strategy",
                    "recommendation": (
                        "An early invasive strategy (coronary angiography within "
                        "24 hours) is recommended in NSTE-ACS patients with at "
                        "least one high-risk criterion including significant "
                        "troponin rise/fall, dynamic ST/T-wave changes, or GRACE "
                        "risk score >140."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Bleeding Risk Assessment",
                    "recommendation": (
                        "Systematic assessment of bleeding risk using validated "
                        "scores (e.g., ARC-HBR, PRECISE-DAPT) is recommended "
                        "at the time of PCI to guide the intensity and duration "
                        "of antithrombotic therapy."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "P2Y12 Inhibitor Selection in ACS",
                    "recommendation": (
                        "Prasugrel or ticagrelor is recommended over clopidogrel "
                        "in patients with NSTE-ACS proceeding to PCI, unless "
                        "contraindicated or at very high bleeding risk."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── Endocarditis ────────────────────────────────────────
        "ESC 2023 Endocarditis": {
            "year": 2023,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehad193",
            "condition": "Infective Endocarditis",
            "sections": [
                {
                    "section": "Modified Duke Criteria for Diagnosis",
                    "recommendation": (
                        "The modified Duke criteria incorporating imaging findings "
                        "(echocardiography, cardiac CT, 18F-FDG PET/CT, and WBC "
                        "SPECT/CT) are recommended for the diagnosis of infective "
                        "endocarditis."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Empiric Antibiotic Therapy",
                    "recommendation": (
                        "Empiric antibiotic therapy should be started promptly "
                        "after obtaining at least 3 sets of blood cultures from "
                        "separate venipuncture sites. In native valve endocarditis, "
                        "ampicillin-sulbactam plus gentamicin; in prosthetic valve, "
                        "vancomycin plus gentamicin plus rifampin."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Surgical Indications — Heart Failure",
                    "recommendation": (
                        "Early surgery is recommended in patients with IE "
                        "complicated by heart failure due to severe aortic or "
                        "mitral regurgitation, intracardiac fistula, or valve "
                        "obstruction."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Surgical Indications — Uncontrolled Infection and Abscess",
                    "recommendation": (
                        "Surgery is recommended for uncontrolled infection "
                        "including persistent bacteremia >7 days despite "
                        "appropriate antibiotics, perivalvular abscess or "
                        "pseudoaneurysm, or infection caused by fungi or "
                        "multi-resistant organisms."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Surgical Indications — Recurrent Emboli",
                    "recommendation": (
                        "Surgery should be considered in patients with IE and "
                        "recurrent embolic events despite appropriate antibiotic "
                        "therapy, particularly when vegetations >10 mm persist "
                        "after one or more embolic episodes."
                    ),
                    "class": "IIa",
                    "evidence": "B-NR",
                },
            ],
        },
        # ── Cardiac Arrest / CPR ────────────────────────────────
        "AHA 2023 Cardiac Arrest/CPR": {
            "year": 2023,
            "society": "AHA",
            "doi": "10.1161/CIR.0000000000001179",
            "condition": "Cardiac Arrest and Resuscitation",
            "sections": [
                {
                    "section": "High-Quality CPR Parameters",
                    "recommendation": (
                        "High-quality CPR is recommended with chest compression "
                        "rate of 100-120/min, depth of at least 2 inches (5 cm), "
                        "full chest recoil, minimizing interruptions (<10 sec), "
                        "and avoiding excessive ventilation."
                    ),
                    "class": "I",
                    "evidence": "B-R",
                },
                {
                    "section": "ECMO-Facilitated Resuscitation (ECPR)",
                    "recommendation": (
                        "ECMO-facilitated resuscitation (ECPR) may be considered "
                        "for select patients with refractory cardiac arrest when "
                        "performed by experienced teams in systems with established "
                        "ECPR protocols."
                    ),
                    "class": "IIb",
                    "evidence": "C-LD",
                },
                {
                    "section": "Dual Sequential Defibrillation",
                    "recommendation": (
                        "Dual sequential defibrillation (applying two sets of "
                        "defibrillation pads and delivering near-simultaneous "
                        "shocks) may be considered in refractory ventricular "
                        "fibrillation unresponsive to standard defibrillation."
                    ),
                    "class": "IIb",
                    "evidence": "C-LD",
                },
                {
                    "section": "Post-Cardiac Arrest Care Bundle",
                    "recommendation": (
                        "A comprehensive post-cardiac arrest care bundle is "
                        "recommended including targeted temperature management "
                        "(32-36 C for >=24 hours), coronary angiography for "
                        "STEMI, hemodynamic optimization, multimodal neuroprognostication "
                        "(>=72 hours after ROSC), and seizure management."
                    ),
                    "class": "I",
                    "evidence": "B-R",
                },
                {
                    "section": "Epinephrine Timing in Cardiac Arrest",
                    "recommendation": (
                        "Epinephrine 1 mg IV/IO should be administered every "
                        "3-5 minutes during cardiac arrest. For non-shockable "
                        "rhythms, early epinephrine administration is recommended."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
            ],
        },
        # ── Myocarditis ─────────────────────────────────────────
        "ACC 2024 Myocarditis": {
            "year": 2024,
            "society": "ACC",
            "doi": "10.1016/j.jacc.2024.03.411",
            "condition": "Myocarditis",
            "sections": [
                {
                    "section": "Endomyocardial Biopsy Indications",
                    "recommendation": (
                        "Endomyocardial biopsy is recommended in patients with "
                        "clinically suspected myocarditis presenting with new-onset "
                        "heart failure with hemodynamic compromise, life-threatening "
                        "arrhythmias, or failure to respond to standard therapy "
                        "within 1-2 weeks."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "CMR Lake Louise Criteria",
                    "recommendation": (
                        "Cardiac magnetic resonance (CMR) using the updated Lake "
                        "Louise criteria (T1 mapping, T2 mapping, and LGE) is "
                        "recommended as the primary noninvasive diagnostic tool "
                        "for clinically suspected myocarditis in hemodynamically "
                        "stable patients."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Immunosuppression for Giant Cell Myocarditis",
                    "recommendation": (
                        "Immunosuppressive therapy with corticosteroids and "
                        "calcineurin inhibitors (cyclosporine or tacrolimus) is "
                        "recommended in biopsy-proven giant cell myocarditis."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Immunosuppression for Eosinophilic Myocarditis",
                    "recommendation": (
                        "High-dose corticosteroids are recommended as first-line "
                        "therapy in biopsy-proven eosinophilic myocarditis after "
                        "exclusion of parasitic and hypersensitivity etiologies."
                    ),
                    "class": "I",
                    "evidence": "C-LD",
                },
                {
                    "section": "Activity Restriction After Myocarditis",
                    "recommendation": (
                        "Competitive athletes and highly active individuals with "
                        "myocarditis should abstain from competitive sports and "
                        "intense exercise for 3-6 months from symptom onset, with "
                        "return to activity guided by resolution of symptoms, "
                        "biomarkers, arrhythmia, and CMR findings."
                    ),
                    "class": "I",
                    "evidence": "C-EO",
                },
            ],
        },
        # ── ATTR Cardiac Amyloidosis ────────────────────────────
        "ESC 2024 ATTR Cardiac Amyloidosis": {
            "year": 2024,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehae404",
            "condition": "ATTR Cardiac Amyloidosis",
            "sections": [
                {
                    "section": "Non-Biopsy Diagnosis of ATTR-CM",
                    "recommendation": (
                        "Non-biopsy diagnosis of ATTR cardiac amyloidosis is "
                        "recommended when bone scintigraphy (99mTc-PYP, DPD, or "
                        "HMDP) shows grade >=2 myocardial uptake AND monoclonal "
                        "protein is excluded by serum/urine immunofixation and "
                        "serum free light chain assay."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Tafamidis for ATTR-CM",
                    "recommendation": (
                        "Tafamidis 80 mg (meglumine) or 61 mg (free acid) daily "
                        "is recommended for patients with ATTR cardiomyopathy "
                        "(wild-type or hereditary) and NYHA class I-III symptoms "
                        "to reduce mortality and cardiovascular hospitalizations."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Diuretic Management in ATTR-CM",
                    "recommendation": (
                        "Loop diuretics are recommended for volume management in "
                        "patients with ATTR-CM and signs of congestion, with "
                        "careful dose titration to avoid hypotension given the "
                        "restrictive physiology and preload dependence."
                    ),
                    "class": "I",
                    "evidence": "C-EO",
                },
                {
                    "section": "Anticoagulation in AF with ATTR-CM",
                    "recommendation": (
                        "Anticoagulation is recommended in all patients with "
                        "cardiac amyloidosis and atrial fibrillation, regardless "
                        "of CHA2DS2-VASc score, given the high thromboembolic "
                        "risk associated with atrial amyloid infiltration."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Avoidance of Specific Medications",
                    "recommendation": (
                        "Calcium channel blockers, digoxin, and beta-blockers "
                        "should generally be avoided or used with extreme "
                        "caution in cardiac amyloidosis due to potential binding "
                        "to amyloid fibrils and worsening of heart failure."
                    ),
                    "class": "III",
                    "evidence": "C-LD",
                },
            ],
        },
        # ── PAD Update ──────────────────────────────────────────
        "ACC/AHA 2024 PAD Update": {
            "year": 2024,
            "society": "ACC/AHA",
            "doi": "10.1016/j.jacc.2024.02.013",
            "condition": "Peripheral Artery Disease",
            "sections": [
                {
                    "section": "ABI Screening for PAD",
                    "recommendation": (
                        "Ankle-brachial index (ABI) measurement is recommended "
                        "as the first-line diagnostic test for PAD in patients "
                        "with exertional leg symptoms, nonhealing wounds, age "
                        ">=65 years, or age >=50 years with risk factors "
                        "(diabetes, smoking, hypertension, dyslipidemia)."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Supervised Exercise Therapy",
                    "recommendation": (
                        "Supervised exercise therapy (SET) is recommended as "
                        "first-line treatment for claudication, consisting of "
                        ">=12 weeks of sessions (>=3 per week, >=30 minutes "
                        "per session) using intermittent walking to near-maximal "
                        "claudication pain."
                    ),
                    "class": "I",
                    "evidence": "A",
                },
                {
                    "section": "Cilostazol for Claudication",
                    "recommendation": (
                        "Cilostazol 100 mg twice daily is recommended to improve "
                        "walking distance and quality of life in patients with "
                        "lifestyle-limiting claudication who do not have heart "
                        "failure."
                    ),
                    "class": "IIa",
                    "evidence": "A",
                },
                {
                    "section": "Revascularization Indications",
                    "recommendation": (
                        "Revascularization (endovascular or surgical) is "
                        "recommended for patients with claudication that is "
                        "lifestyle-limiting despite exercise and pharmacotherapy, "
                        "with adequate anatomy for intervention."
                    ),
                    "class": "IIa",
                    "evidence": "B-R",
                },
                {
                    "section": "Critical Limb-Threatening Ischemia Management",
                    "recommendation": (
                        "In patients with chronic limb-threatening ischemia "
                        "(CLTI), urgent vascular assessment and revascularization "
                        "are recommended to relieve rest pain, heal wounds, and "
                        "prevent limb loss. A multidisciplinary team approach "
                        "including wound care is recommended."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
            ],
        },
        # ── Cardio-Oncology ──────────────────────────────────────
        "ESC 2022 Cardio-Oncology": {
            "year": 2022,
            "society": "ESC",
            "doi": "10.1093/eurheartj/ehac244",
            "condition": "Cardio-Oncology",
            "sections": [
                {
                    "section": "Baseline CV Assessment Before Cardiotoxic Therapy",
                    "recommendation": (
                        "Baseline cardiovascular assessment including "
                        "echocardiography (LVEF, GLS), ECG, and cardiac "
                        "biomarkers (troponin, BNP) is recommended before "
                        "starting potentially cardiotoxic cancer therapy."
                    ),
                    "class": "I",
                    "evidence": "B-NR",
                },
                {
                    "section": "Cardioprotection During Anthracyclines",
                    "recommendation": (
                        "ACEi/ARB and/or beta-blocker should be considered "
                        "for primary prevention of anthracycline cardiotoxicity "
                        "in high-risk patients."
                    ),
                    "class": "IIa",
                    "evidence": "B-R",
                },
            ],
        },
    }
    """Curated registry of major cardiovascular clinical practice guidelines."""

    def __init__(self):
        """Initialize the guideline parser."""
        super().__init__("cardio_guidelines")

    # ── Fetch ─────────────────────────────────────────────────────────

    def fetch(self, **kwargs) -> List[dict]:
        """Return the static guideline registry as a list of dicts.

        Since guideline full text is copyright-protected, this parser
        operates on curated seed data rather than fetching from an API.

        Returns:
            List of guideline dictionaries from :attr:`GUIDELINES`.
        """
        guidelines = []
        for name, info in self.GUIDELINES.items():
            entry = dict(info)
            entry["guideline_name"] = name
            guidelines.append(entry)
        self.logger.info(f"Loaded {len(guidelines)} guideline entries")
        return guidelines

    # ── Parse ─────────────────────────────────────────────────────────

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Parse guideline entries into IngestRecord instances.

        Each guideline may have multiple sections; each section becomes
        a separate IngestRecord to enable fine-grained retrieval.

        Args:
            raw_data: List of guideline dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            guideline_name = entry.get("guideline_name", "")
            society = entry.get("society", "")
            year = entry.get("year", 0)
            doi = entry.get("doi", "")
            condition = entry.get("condition", "")
            sections = entry.get("sections", [])

            for sec in sections:
                section_name = sec.get("section", "")
                recommendation = sec.get("recommendation", "")
                class_of_rec = sec.get("class", "")
                evidence_level = sec.get("evidence", "")

                if not recommendation:
                    continue

                text = (
                    f"{guideline_name} ({year}) — {section_name}: "
                    f"{recommendation}"
                )

                metadata = {
                    "society": society,
                    "guideline_title": self.truncate(guideline_name, 510),
                    "year": year,
                    "recommendation": self.truncate(recommendation, 2046),
                    "class_of_rec": class_of_rec,
                    "evidence_level": evidence_level,
                    "condition": self.truncate(condition, 126),
                    "section": self.truncate(section_name, 254),
                }

                records.append(
                    IngestRecord(
                        text=self.truncate(text, 2048),
                        metadata=metadata,
                        collection=self.collection,
                        source=f"{society} {year}",
                        source_id=doi,
                    )
                )

        self.logger.info(
            f"Parsed {len(records)} guideline recommendation records "
            f"from {len(raw_data)} guidelines"
        )
        return records

    # ── Convenience seed method ───────────────────────────────────────

    def seed_guidelines(self) -> List[IngestRecord]:
        """Convenience wrapper: fetch + parse + filter in one call.

        Equivalent to calling ``run()`` but with a clearer name for
        the static-seeding use case.

        Returns:
            List of valid :class:`IngestRecord` instances.
        """
        return self.run()
