"""Patient: Emma Williams -- 34F, Multiple Sclerosis.

Timeline: March 2024 -> March 2026
Specialists seen: PCP, Ophthalmology, Neurology, Radiology
Key pattern: Optic neuritis -> periventricular MRI lesions -> oligoclonal bands
           -> McDonald criteria met -> Ocrelizumab treatment

The system should detect: Relapsing-remitting MS, HLA-DRB1*15:01 risk.
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
    generate_lab_report,
    generate_progress_note,
)

PATIENT = {
    "name": "Emma Williams",
    "dob": "1991-09-22",
    "mrn": "EWI-2024-55912",
    "sex": "F",
    "age_at_start": 32,
}


def generate(output_dir: str):
    """Generate all clinical documents for Emma Williams."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: PCP -- March 2024 (blurred vision, tingling)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-03-11_PCP_Progress_Note.pdf"),
        PATIENT, "03/11/2024", "pcp",
        "OFFICE VISIT -- NEW VISUAL SYMPTOMS AND PARESTHESIAS",
        {
            "chief_complaint": "Blurred vision in right eye x 5 days. Tingling in left hand and forearm x 2 weeks.",
            "hpi": (
                "Ms. Williams is a 32-year-old woman presenting with subacute onset of blurred "
                "vision in her right eye over the past 5 days. She reports the vision is 'foggy' "
                "centrally and worsened over the first 3 days, now stable. She also notes pain "
                "behind the right eye, worsened with eye movement. Separately, she has experienced "
                "intermittent tingling and numbness in her left hand and forearm for approximately "
                "2 weeks. She attributes the hand tingling to carpal tunnel from computer work "
                "(she is a graphic designer). No recent illness, trauma, or vaccination. No fever. "
                "No headache. No weakness. No bowel or bladder changes. No prior neurological symptoms."
            ),
            "vitals": "BP 116/72  HR 68  Temp 98.6F  SpO2 99%  Wt 138 lbs  Ht 5'6\"",
            "ros": (
                "Constitutional: Mild fatigue over past month. No fever, weight change.\n"
                "HEENT: Right eye blurred vision and retro-orbital pain as above. No diplopia.\n"
                "Neurological: Left hand/forearm tingling. No weakness. No gait problems.\n"
                "Musculoskeletal: No joint pain or stiffness.\n"
                "Psychiatric: Mild stress from work deadlines. Sleep adequate."
            ),
            "exam": (
                "General: Well-appearing woman in no acute distress.\n"
                "HEENT: Right pupil mildly sluggish to direct light, brisk consensual. "
                "Left pupil normal. Possible relative afferent pupillary defect (RAPD) on right.\n"
                "Visual acuity: OD 20/50 (baseline 20/20 per patient). OS 20/20.\n"
                "Color vision: OD impaired (Ishihara 8/14 plates). OS normal (14/14).\n"
                "Neurological: Sensation decreased to light touch left C6-C8 dermatome. "
                "Strength 5/5 all extremities. DTRs 2+ symmetric. Babinski equivocal on right."
            ),
            "assessment": (
                "1. Acute right optic neuritis -- subacute monocular visual loss with pain on "
                "eye movement, RAPD, and dyschromatopsia. URGENT ophthalmology referral.\n"
                "2. Left upper extremity paresthesias -- dermatomal pattern, concerning in "
                "context of optic neuritis for demyelinating disease vs. cervical radiculopathy.\n"
                "3. Need to rule out multiple sclerosis -- optic neuritis in a young woman with "
                "concurrent paresthesias warrants neurological evaluation."
            ),
            "plan": (
                "1. URGENT ophthalmology referral -- same day if possible\n"
                "2. Neurology referral -- expedited, within 1-2 weeks\n"
                "3. MRI brain with and without contrast -- ordered STAT\n"
                "4. CBC, CMP, ESR, CRP, vitamin B12, TSH, ANA\n"
                "5. Patient counseled on warning signs (vision worsening, weakness, "
                "bladder dysfunction) -- ER if develops\n"
                "6. Return after specialist evaluations"
            ),
        },
    )

    # Lab: March 2024 -- initial workup
    generate_lab_report(
        os.path.join(output_dir, "2024-03-11_Lab_Initial.pdf"),
        PATIENT, "03/11/2024", PROVIDERS["pcp"],
        [{
            "panel_name": "COMPLETE BLOOD COUNT (CBC)",
            "results": [
                {"test": "WBC", "value": "6.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "242", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }, {
            "panel_name": "COMPREHENSIVE METABOLIC PANEL",
            "results": [
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "ALT", "value": "19", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "22", "unit": "U/L", "ref_range": "10-40", "flag": ""},
            ],
        }, {
            "panel_name": "ADDITIONAL STUDIES",
            "results": [
                {"test": "ESR", "value": "14", "unit": "mm/hr", "ref_range": "0-20", "flag": ""},
                {"test": "CRP", "value": "2.1", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "TSH", "value": "2.4", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Vitamin B12", "value": "520", "unit": "pg/mL", "ref_range": "200-900", "flag": ""},
                {"test": "ANA", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: Ophthalmology -- March 2024 (optic neuritis)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-03-13_Ophth_Optic_Neuritis.pdf"),
        PATIENT, "03/13/2024", "ophthalmology",
        "OPHTHALMOLOGY CONSULTATION -- ACUTE OPTIC NEURITIS",
        {
            "chief_complaint": "Referred for evaluation of right eye visual loss with pain on eye movement.",
            "hpi": (
                "Ms. Williams is a 32-year-old woman referred urgently by Dr. Martinez for subacute "
                "right visual loss over 5-7 days associated with retro-orbital pain exacerbated by "
                "eye movement. No prior episodes. She also reports new left hand paresthesias. "
                "No history of autoimmune disease. Family history notable for maternal aunt with "
                "lupus but no MS."
            ),
            "exam": (
                "Visual Acuity: OD 20/60 (pinhole: no improvement). OS 20/20.\n"
                "Color Vision (Ishihara): OD 6/14 plates. OS 14/14.\n"
                "Pupil: Right RAPD confirmed (0.6 log unit).\n"
                "Confrontation Fields: OD central scotoma. OS full.\n"
                "Slit Lamp: Anterior segment quiet OU. No uveitis. No keratic precipitates.\n"
                "IOP: OD 15, OS 14.\n"
                "Fundoscopy: OD optic disc mildly edematous with blurred margins. Peripapillary "
                "flame hemorrhage. No pallor. OS optic disc normal, well-defined margins.\n"
                "OCT RNFL: OD thickened (118 mcm, norm 80-100) consistent with acute edema. OS normal (92 mcm).\n"
                "Visual evoked potential: OD prolonged P100 latency 128 ms (norm <102). OS normal 98 ms."
            ),
            "assessment": (
                "1. Acute demyelinating optic neuritis, right eye -- classic presentation: "
                "young woman, subacute visual loss, pain with eye movement, RAPD, "
                "dyschromatopsia, disc edema, prolonged VEP.\n"
                "2. Concurrent left upper extremity paresthesias raise strong concern for "
                "multiple sclerosis -- awaiting MRI brain."
            ),
            "plan": (
                "1. IV methylprednisolone 1g daily x 3 days (ONTT protocol) followed by oral "
                "prednisone taper -- scheduled for infusion center\n"
                "2. Expedite MRI brain and cervical spine with contrast -- critical for MS workup\n"
                "3. Neurology referral confirmed -- Dr. Kim, appointment in 5 days\n"
                "4. Follow-up in 2 weeks to reassess vision\n"
                "5. Patient counseled: 50% of patients with optic neuritis develop MS within 15 "
                "years; risk higher if MRI shows lesions"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: MRI Brain -- March 2024
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2024-03-14_MRI_Brain.pdf"),
        PATIENT, "03/14/2024", "neurology",
        "MRI", "Brain with and without Gadolinium",
        "32-year-old female with acute right optic neuritis and left upper extremity paresthesias. "
        "Evaluate for demyelinating disease.",
        "MRI brain performed on 3T scanner with sagittal T1, axial T2, FLAIR, diffusion-weighted, "
        "and post-gadolinium T1-weighted sequences.",
        (
            "White matter:\n"
            "- Multiple (>9) ovoid T2/FLAIR hyperintense lesions identified in the periventricular "
            "white matter, oriented perpendicular to the lateral ventricles (Dawson fingers pattern).\n"
            "- Lesion sizes range from 3 mm to 12 mm in longest dimension.\n"
            "- Three juxtacortical lesions identified: left frontal, right parietal, left temporal.\n"
            "- Two infratentorial lesions: right middle cerebellar peduncle (6 mm) and dorsal pons (4 mm).\n"
            "- Post-gadolinium: TWO lesions show ring enhancement -- left periventricular (8 mm) "
            "and right frontal juxtacortical (5 mm), indicating active demyelination.\n"
            "- Right optic nerve: abnormal T2 signal and enhancement in the intraorbital segment, "
            "consistent with acute optic neuritis.\n\n"
            "Other findings:\n"
            "- No mass lesion, midline shift, or hydrocephalus.\n"
            "- Normal gray matter signal. No restricted diffusion to suggest acute infarct.\n"
            "- Normal posterior fossa structures apart from the noted lesions."
        ),
        (
            "1. Multiple periventricular, juxtacortical, and infratentorial T2/FLAIR hyperintense "
            "lesions with Dawson fingers morphology -- HIGHLY SUSPICIOUS FOR DEMYELINATING DISEASE "
            "(multiple sclerosis).\n"
            "2. Two enhancing lesions indicating active demyelination, dissemination in time.\n"
            "3. Right optic nerve enhancement consistent with acute optic neuritis.\n"
            "4. Lesion distribution satisfies 3 of 4 spatial criteria for McDonald 2017 "
            "(periventricular, juxtacortical, infratentorial). Spinal cord imaging recommended.\n\n"
            "CLINICAL CORRELATION: These findings, combined with clinical optic neuritis, "
            "fulfill McDonald 2017 criteria for MS (dissemination in space AND time)."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: MRI Spine -- March 2024
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2024-03-14_MRI_Cervical_Spine.pdf"),
        PATIENT, "03/14/2024", "neurology",
        "MRI", "Cervical and Thoracic Spine with and without Gadolinium",
        "32-year-old female with acute optic neuritis and brain MRI showing multiple "
        "demyelinating lesions. Evaluate for spinal cord involvement.",
        "MRI cervical and thoracic spine with sagittal T1, T2, STIR, and post-gadolinium "
        "T1-weighted sequences on 3T scanner.",
        (
            "Cervical spine:\n"
            "- Two intramedullary T2 hyperintense lesions identified:\n"
            "  * C3-C4 level: 8 mm lesion, posterolateral cord, non-enhancing (chronic)\n"
            "  * C5 level: 5 mm lesion, lateral cord on left, no enhancement\n"
            "- No cord expansion or syrinx.\n"
            "- Normal disc heights and vertebral body signal.\n\n"
            "Thoracic spine:\n"
            "- One T2 hyperintense lesion at T6, 6 mm, dorsal cord, non-enhancing.\n"
            "- No additional lesions identified.\n"
            "- Normal conus medullaris."
        ),
        (
            "1. Three spinal cord lesions (two cervical, one thoracic) consistent with "
            "demyelinating disease. Non-enhancing, suggesting chronic/subacute lesions.\n"
            "2. Combined with brain MRI findings, this completes 4 of 4 spatial criteria "
            "for McDonald 2017 (periventricular, juxtacortical, infratentorial, spinal cord).\n"
            "3. C3-C4 lesion location correlates with patient's reported left upper extremity "
            "paresthesias (posterior column involvement)."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: Neurology Consult -- March 2024 (MS diagnosis)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-03-18_Neurology_Consult.pdf"),
        PATIENT, "03/18/2024", "neurology",
        "NEUROLOGY CONSULTATION -- MULTIPLE SCLEROSIS EVALUATION",
        {
            "chief_complaint": (
                "Referred for evaluation of acute optic neuritis with MRI findings "
                "concerning for multiple sclerosis."
            ),
            "hpi": (
                "Ms. Williams is a 32-year-old woman presenting with acute right optic neuritis "
                "treated with IV methylprednisolone (completed 3-day course). Vision is improving. "
                "MRI brain shows >9 T2/FLAIR lesions in periventricular, juxtacortical, and "
                "infratentorial locations with 2 enhancing lesions. MRI spine shows 3 additional "
                "spinal cord lesions. She also has left hand paresthesias correlating with C3-C4 "
                "cord lesion. No prior neurological episodes, though she recalls a 2-week episode "
                "of right leg heaviness 6 months ago that resolved spontaneously (possible prior "
                "relapse). No family history of MS. Maternal aunt with lupus."
            ),
            "exam": (
                "Mental Status: Alert, oriented x4. Normal affect. MMSE 30/30.\n"
                "Cranial Nerves: Right RAPD (improving). Visual acuity OD 20/40 (improving). "
                "OS 20/20. Color vision OD 10/14 (improved from 6/14). Remaining CN intact.\n"
                "Motor: 5/5 all extremities. No drift. Normal tone.\n"
                "Sensory: Decreased light touch and pinprick left C6-C8. Decreased vibration "
                "left fingers (128 Hz tuning fork: 4 sec vs 12 sec right). Lhermitte sign POSITIVE.\n"
                "Reflexes: 3+ bilateral lower extremities (hyperreflexia). Right Babinski positive. "
                "Left Babinski equivocal. Hoffman sign positive bilaterally.\n"
                "Cerebellar: No dysmetria. No ataxia.\n"
                "Gait: Normal. Tandem gait mildly unsteady."
            ),
            "assessment": (
                "DIAGNOSIS: RELAPSING-REMITTING MULTIPLE SCLEROSIS (RRMS)\n\n"
                "McDonald 2017 criteria fulfilled:\n"
                "- Dissemination in space: periventricular (>3 lesions), juxtacortical (3 lesions), "
                "infratentorial (2 lesions), spinal cord (3 lesions) -- 4/4 criteria met\n"
                "- Dissemination in time: simultaneous enhancing and non-enhancing lesions on MRI "
                "plus clinical optic neuritis (current) and possible prior relapse (right leg "
                "heaviness 6 months ago)\n\n"
                "EDSS Score: 2.0 (minimal disability -- visual function and sensory findings)\n\n"
                "This is a young woman with high lesion burden at presentation, enhancing lesions, "
                "and spinal cord involvement -- features indicating ACTIVE disease warranting "
                "high-efficacy therapy."
            ),
            "plan": (
                "1. Lumbar puncture -- CSF for oligoclonal bands, IgG index, cell count, protein, "
                "glucose, cytology. Supports diagnosis and rules out mimics.\n"
                "2. JC virus antibody testing (anti-JCV Ab index) -- baseline before treatment decision\n"
                "3. Additional labs: Vitamin D 25-OH, NMO-IgG (aquaporin-4 Ab), MOG antibodies "
                "(rule out NMOSD/MOGAD), quantitative immunoglobulins, hepatitis B/C serologies\n"
                "4. Neurofilament light chain (NfL) -- serum biomarker of neuroaxonal injury, "
                "baseline for monitoring\n"
                "5. Treatment discussion: Given high lesion burden and active disease, recommend "
                "HIGH-EFFICACY disease-modifying therapy. Top options:\n"
                "   a. Ocrelizumab (Ocrevus) -- anti-CD20, preferred given efficacy profile\n"
                "   b. Natalizumab (Tysabri) -- if JCV negative, highly effective\n"
                "6. EDSS baseline documented: 2.0\n"
                "7. MS specialist nurse education appointment scheduled\n"
                "8. Follow-up in 2 weeks with LP and lab results for treatment decision"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Lumbar Puncture Results -- March 2024
    # ─────────────────────────────────────────────────────────────
    generate_lab_report(
        os.path.join(output_dir, "2024-03-20_Lab_Lumbar_Puncture.pdf"),
        PATIENT, "03/20/2024", PROVIDERS["neurology"],
        [{
            "panel_name": "CEREBROSPINAL FLUID ANALYSIS",
            "results": [
                {"test": "CSF Appearance", "value": "Clear, colorless", "unit": "", "ref_range": "Clear, colorless", "flag": ""},
                {"test": "CSF WBC", "value": "8", "unit": "cells/mcL", "ref_range": "0-5", "flag": "H"},
                {"test": "CSF RBC", "value": "0", "unit": "cells/mcL", "ref_range": "0", "flag": ""},
                {"test": "CSF Protein", "value": "52", "unit": "mg/dL", "ref_range": "15-45", "flag": "H"},
                {"test": "CSF Glucose", "value": "62", "unit": "mg/dL", "ref_range": "40-70", "flag": ""},
                {"test": "CSF IgG", "value": "6.8", "unit": "mg/dL", "ref_range": "0.0-3.4", "flag": "HH"},
                {"test": "IgG Index", "value": "0.92", "unit": "", "ref_range": "<0.70", "flag": "H"},
                {"test": "IgG Synthesis Rate", "value": "14.2", "unit": "mg/day", "ref_range": "<3.3", "flag": "HH"},
                {"test": "Oligoclonal Bands", "value": "6 bands", "unit": "", "ref_range": "0-1", "flag": "A"},
                {"test": "Oligoclonal Bands (CSF-specific)", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "CSF Cytology", "value": "Negative for malignancy", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }, {
            "panel_name": "PAIRED SERUM",
            "results": [
                {"test": "Serum IgG", "value": "1040", "unit": "mg/dL", "ref_range": "700-1600", "flag": ""},
                {"test": "Serum Albumin", "value": "4.2", "unit": "g/dL", "ref_range": "3.5-5.0", "flag": ""},
            ],
        }],
    )

    # Lab: JCV and additional MS workup
    generate_lab_report(
        os.path.join(output_dir, "2024-03-20_Lab_MS_Workup.pdf"),
        PATIENT, "03/20/2024", PROVIDERS["neurology"],
        [{
            "panel_name": "MS BIOMARKERS AND SAFETY SCREENING",
            "results": [
                {"test": "JC Virus Antibody", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "JCV Ab Index", "value": "1.8", "unit": "", "ref_range": "<0.9 low risk", "flag": "H"},
                {"test": "Neurofilament Light Chain (NfL)", "value": "32.4", "unit": "pg/mL", "ref_range": "<10.0", "flag": "HH"},
                {"test": "NMO-IgG (Aquaporin-4 Ab)", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "MOG Antibody", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "18", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
            ],
        }, {
            "panel_name": "HEPATITIS AND IMMUNOGLOBULIN SCREENING",
            "results": [
                {"test": "Hepatitis B Surface Antigen", "value": "Non-reactive", "unit": "", "ref_range": "Non-reactive", "flag": ""},
                {"test": "Hepatitis B Core Antibody", "value": "Non-reactive", "unit": "", "ref_range": "Non-reactive", "flag": ""},
                {"test": "Hepatitis C Antibody", "value": "Non-reactive", "unit": "", "ref_range": "Non-reactive", "flag": ""},
                {"test": "IgG", "value": "1040", "unit": "mg/dL", "ref_range": "700-1600", "flag": ""},
                {"test": "IgM", "value": "125", "unit": "mg/dL", "ref_range": "40-230", "flag": ""},
                {"test": "IgA", "value": "198", "unit": "mg/dL", "ref_range": "70-400", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 7: Neurology Follow-up -- April 2024 (treatment decision)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-04-01_Neurology_Treatment_Decision.pdf"),
        PATIENT, "04/01/2024", "neurology",
        "NEUROLOGY FOLLOW-UP -- MS TREATMENT DECISION",
        {
            "chief_complaint": "Follow-up RRMS. Discuss treatment options with LP and lab results.",
            "hpi": (
                "Ms. Williams returns 2 weeks after MS diagnosis. Right eye vision has improved "
                "significantly (now 20/25). Left hand tingling persists but is less bothersome. "
                "She has been processing the diagnosis with family support. LP confirmed CSF "
                "oligoclonal bands (6 bands) and elevated IgG index (0.92). NfL significantly "
                "elevated at 32.4 pg/mL indicating active neuroaxonal injury. NMO and MOG "
                "antibodies negative, confirming MS. JCV antibody POSITIVE with index 1.8 -- "
                "this influences treatment selection (natalizumab carries PML risk with JCV+)."
            ),
            "assessment": (
                "Relapsing-remitting multiple sclerosis -- CONFIRMED\n\n"
                "Favorable prognostic features: young age, female, optic neuritis as presenting event\n"
                "Unfavorable features: high lesion burden (>9 brain, 3 spinal), enhancing lesions, "
                "elevated NfL (32.4), low vitamin D\n\n"
                "TREATMENT SELECTION: Ocrelizumab (Ocrevus) 600 mg IV q6m\n"
                "Rationale: High-efficacy therapy needed given active disease. JCV+ status "
                "makes natalizumab higher risk for PML. Ocrelizumab has robust phase III data "
                "(OPERA I/II) showing 46-47% reduction in annualized relapse rate vs IFN beta-1a."
            ),
            "plan": (
                "1. Start Ocrelizumab (Ocrevus):\n"
                "   - Initial dose: 300 mg IV x 2 infusions, 2 weeks apart\n"
                "   - Subsequent doses: 600 mg IV every 6 months\n"
                "   - Pre-medication: methylprednisolone 100 mg IV, diphenhydramine, acetaminophen\n"
                "2. Pre-infusion labs: CBC, hepatitis B (completed), quantitative immunoglobulins\n"
                "3. Vitamin D supplementation: cholecalciferol 5000 IU daily (level 18, goal >40)\n"
                "4. Baseline EDSS: 2.0 -- reassess at 6 months\n"
                "5. MRI brain and spine in 6 months (new baseline on treatment)\n"
                "6. NfL recheck in 6 months -- goal: normalization (<10 pg/mL)\n"
                "7. MS specialist nurse for injection training, symptom management education\n"
                "8. Referral to MS support group\n"
                "9. First infusion scheduled: April 15, 2024"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 8: Ocrelizumab Infusion Note -- April 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-04-15_Infusion_Ocrelizumab_Dose1.pdf"),
        PATIENT, "04/15/2024", "neurology",
        "INFUSION NOTE -- OCRELIZUMAB (OCREVUS) DOSE 1 OF 2",
        {
            "chief_complaint": "First Ocrelizumab infusion for relapsing-remitting MS.",
            "hpi": (
                "Ms. Williams presents for first Ocrelizumab infusion (loading dose 1 of 2). "
                "She has completed pre-infusion labs which are satisfactory. Hepatitis B negative. "
                "Immunoglobulins within normal limits. She reports continued improvement in right eye "
                "vision (now 20/25). Left hand tingling is mild and intermittent."
            ),
            "vitals": "BP 118/74  HR 72  Temp 98.4F  SpO2 99%  Wt 137 lbs",
            "medications": (
                "Pre-medications administered:\n"
                "- Methylprednisolone 100 mg IV\n"
                "- Diphenhydramine 50 mg IV\n"
                "- Acetaminophen 650 mg PO\n\n"
                "Infusion:\n"
                "- Ocrelizumab 300 mg IV over 2.5 hours\n"
                "- Start rate 30 mL/hr, increased to 60 mL/hr at 30 min, then 180 mL/hr at 60 min"
            ),
            "assessment": (
                "RRMS, first Ocrelizumab infusion (loading dose 1/2). Infusion completed without "
                "infusion-related reaction. Patient monitored for 1 hour post-infusion. Vital signs "
                "stable throughout. No throat tightness, pruritus, flushing, or hypotension."
            ),
            "plan": (
                "1. Second loading dose: Ocrelizumab 300 mg IV -- scheduled April 29, 2024\n"
                "2. Monitor for delayed reactions: rash, fever, malaise in next 24-48 hours\n"
                "3. Continue vitamin D 5000 IU daily\n"
                "4. Call clinic if any infusion-related symptoms develop\n"
                "5. Avoid live vaccines"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 9: 6-month MRI -- October 2024
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2024-10-14_MRI_Brain_6Month.pdf"),
        PATIENT, "10/14/2024", "neurology",
        "MRI", "Brain with and without Gadolinium",
        "33-year-old female with RRMS on Ocrelizumab x 6 months. Surveillance scan.",
        "MRI brain 3T with sagittal T1, axial T2, FLAIR, and post-gadolinium T1. "
        "Comparison: March 14, 2024.",
        (
            "White matter lesions:\n"
            "- Previously identified periventricular, juxtacortical, and infratentorial "
            "T2/FLAIR lesions are STABLE in number and size.\n"
            "- No new T2/FLAIR lesions identified.\n"
            "- Previously enhancing left periventricular and right frontal lesions are NO LONGER "
            "enhancing -- resolved active inflammation.\n"
            "- No new enhancing lesions.\n"
            "- Right optic nerve: resolution of previously seen T2 signal abnormality and "
            "enhancement. Normal optic nerve appearance.\n\n"
            "Brain volume: Within normal limits. No interval volume loss appreciated on "
            "qualitative assessment."
        ),
        (
            "1. STABLE MS lesion burden. No new or enlarging T2/FLAIR lesions.\n"
            "2. RESOLVED gadolinium enhancement -- no active inflammation.\n"
            "3. Resolution of right optic neuritis on imaging.\n"
            "4. No evidence of disease activity (NEDA criteria met on imaging).\n\n"
            "IMPRESSION: Favorable treatment response to Ocrelizumab. No Evidence of "
            "Disease Activity on imaging."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 10: Neurology 6-month follow-up -- October 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-10-21_Neurology_6Month.pdf"),
        PATIENT, "10/21/2024", "neurology",
        "NEUROLOGY FOLLOW-UP -- 6 MONTHS ON OCRELIZUMAB",
        {
            "chief_complaint": "Routine MS follow-up. No new symptoms.",
            "hpi": (
                "Ms. Williams returns for 6-month follow-up on Ocrelizumab. She has received "
                "the loading doses and first maintenance dose without infusion reactions. She reports "
                "complete resolution of right eye visual symptoms. Left hand tingling has improved "
                "significantly, now only noticeable when fatigued. No new neurological symptoms. "
                "Energy level good. Working full-time without limitations. Exercising regularly "
                "(yoga 3x/week, walking). Taking vitamin D 5000 IU daily."
            ),
            "exam": (
                "Visual Acuity: OD 20/20 (recovered). OS 20/20.\n"
                "Pupil: No RAPD (resolved).\n"
                "Color Vision: OD 13/14 (near-normal). OS 14/14.\n"
                "Sensory: Minimal decreased vibration left fingers (improved).\n"
                "Reflexes: 2+ symmetric (normalized). Babinski negative bilaterally.\n"
                "Gait: Normal, including tandem."
            ),
            "assessment": (
                "RRMS -- EXCELLENT TREATMENT RESPONSE\n\n"
                "EDSS: 1.0 (improved from 2.0 at diagnosis)\n\n"
                "NEDA-3 criteria at 6 months: MET\n"
                "  1. No relapses -- MET\n"
                "  2. No new/enlarging MRI lesions -- MET\n"
                "  3. No confirmed disability worsening -- MET (EDSS improved)\n\n"
                "NfL decreased from 32.4 to 8.6 pg/mL -- normalizing, excellent biomarker response.\n"
                "Vitamin D improved from 18 to 42 ng/mL on supplementation."
            ),
            "plan": (
                "1. Continue Ocrelizumab 600 mg IV q6m -- next infusion April 2025\n"
                "2. Continue vitamin D 5000 IU daily\n"
                "3. Annual MRI brain and spine\n"
                "4. NfL and immunoglobulins q6m with each infusion\n"
                "5. Annual EDSS assessment\n"
                "6. Ophthalmology OCT annually to monitor RNFL\n"
                "7. Discuss pregnancy planning timeline (Ocrelizumab washout needed)\n"
                "8. Follow-up in 6 months"
            ),
        },
    )

    # Lab: 6-month monitoring
    generate_lab_report(
        os.path.join(output_dir, "2024-10-21_Lab_MS_Monitoring.pdf"),
        PATIENT, "10/21/2024", PROVIDERS["neurology"],
        [{
            "panel_name": "MS TREATMENT MONITORING",
            "results": [
                {"test": "Neurofilament Light Chain (NfL)", "value": "8.6", "unit": "pg/mL", "ref_range": "<10.0", "flag": ""},
                {"test": "IgG", "value": "860", "unit": "mg/dL", "ref_range": "700-1600", "flag": ""},
                {"test": "IgM", "value": "68", "unit": "mg/dL", "ref_range": "40-230", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "42", "unit": "ng/mL", "ref_range": "30-100", "flag": ""},
            ],
        }, {
            "panel_name": "CBC WITH DIFFERENTIAL",
            "results": [
                {"test": "WBC", "value": "4.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Lymphocytes", "value": "0.8", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": "L"},
                {"test": "CD19 B Cells", "value": "0", "unit": "cells/mcL", "ref_range": "100-500", "flag": "L"},
                {"test": "Hemoglobin", "value": "13.6", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "235", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient Emma Williams: {count} documents generated")


if __name__ == "__main__":
    generate(os.path.join(os.path.dirname(__file__), "..", "demo_data", "emma_williams"))
