"""Patient 1: Sarah Mitchell -- 34F, 3.5-year lupus diagnostic odyssey.

Timeline: June 2022 -> December 2025
Specialists seen: PCP, Dermatology, Rheumatology, Nephrology, Cardiology
Key pattern: Scattered ANA, rising anti-dsDNA, falling C3/C4, proteinuria
across 5 specialists who each missed the unified picture.

The system should detect: SLE with lupus nephritis, HLA-DRB1*03:01 risk.
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_genetic_report,
    generate_imaging_report,
    generate_lab_report,
    generate_pathology_report,
    generate_progress_note,
)

PATIENT = {
    "name": "Sarah Mitchell",
    "dob": "1991-08-14",
    "mrn": "SMI-2022-44871",
    "sex": "F",
    "age_at_start": 30,
}


def generate(output_dir: str):
    """Generate all clinical documents for Sarah Mitchell."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: PCP -- June 2022 (fatigue, joint pain)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-06-15_PCP_Progress_Note.pdf"),
        PATIENT, "06/15/2022", "pcp",
        "OFFICE VISIT -- ANNUAL PHYSICAL WITH NEW COMPLAINTS",
        {
            "chief_complaint": "Annual wellness exam. New complaints of fatigue and intermittent joint pain x 3 months.",
            "hpi": (
                "Ms. Mitchell is a 30-year-old woman presenting for annual physical. She reports progressive "
                "fatigue over the past 3 months that is not relieved by rest. She has noticed intermittent "
                "pain and stiffness in her fingers and wrists, worse in the morning, typically lasting 30-45 "
                "minutes. She denies any recent illness, travel, or new medications. She attributes fatigue "
                "to work stress (software engineer, long hours). No fever, weight loss, or night sweats. "
                "No rash, hair loss, or oral ulcers. No chest pain or shortness of breath."
            ),
            "vitals": "BP 118/74  HR 78  Temp 98.4F  SpO2 99%  Wt 132 lbs  Ht 5'5\"",
            "ros": (
                "Constitutional: Fatigue as above. No fever, chills, or weight change.\n"
                "HEENT: No oral ulcers, dry eyes, or dry mouth.\n"
                "Cardiovascular: No chest pain, palpitations, or edema.\n"
                "Pulmonary: No cough, dyspnea, or pleurisy.\n"
                "Musculoskeletal: Bilateral hand/wrist stiffness and pain as above. No swelling noted by patient.\n"
                "Skin: No rash. No photosensitivity.\n"
                "Neurological: Occasional brain fog. No headaches, seizures, or focal deficits.\n"
                "Psychiatric: Mild anxiety related to work. Sleep adequate."
            ),
            "exam": (
                "General: Well-appearing woman in no acute distress.\n"
                "HEENT: PERRL, oral mucosa clear, no ulcers, no lymphadenopathy.\n"
                "Cardiovascular: RRR, no murmurs, no JVD, no peripheral edema.\n"
                "Pulmonary: Clear to auscultation bilaterally.\n"
                "Musculoskeletal: Mild tenderness to palpation at MCP joints bilaterally. No frank synovitis. "
                "Full ROM all joints. No joint effusions.\n"
                "Skin: No rashes. No alopecia.\n"
                "Neurological: Alert, oriented x4. CN II-XII intact."
            ),
            "assessment": (
                "1. Fatigue -- likely multifactorial (work stress, possible iron deficiency)\n"
                "2. Arthralgias, bilateral hands -- likely repetitive strain from computer work, "
                "though morning stiffness duration is somewhat longer than typical mechanical causes\n"
                "3. Annual wellness exam -- up to date on preventive care"
            ),
            "plan": (
                "1. CBC, CMP, TSH, iron studies, vitamin D to evaluate fatigue\n"
                "2. Empiric trial of ergonomic wrist supports and OTC NSAIDs for hand symptoms\n"
                "3. If joint symptoms persist > 4-6 weeks, will check inflammatory markers (CRP, ESR)\n"
                "4. Continue current medications (OCP)\n"
                "5. Return in 3 months or sooner if symptoms worsen"
            ),
        },
    )

    # Lab: June 2022 -- basic workup (mostly normal, mild anemia)
    generate_lab_report(
        os.path.join(output_dir, "2022-06-15_Lab_CBC_CMP.pdf"),
        PATIENT, "06/15/2022", PROVIDERS["pcp"],
        [{
            "panel_name": "COMPLETE BLOOD COUNT (CBC)",
            "results": [
                {"test": "WBC", "value": "5.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "RBC", "value": "4.1", "unit": "M/uL", "ref_range": "3.8-5.2", "flag": ""},
                {"test": "Hemoglobin", "value": "11.8", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": "L"},
                {"test": "Hematocrit", "value": "35.2", "unit": "%", "ref_range": "36.0-46.0", "flag": "L"},
                {"test": "MCV", "value": "86", "unit": "fL", "ref_range": "80-100", "flag": ""},
                {"test": "Platelets", "value": "198", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Lymphocytes", "value": "1.4", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": ""},
            ],
        }, {
            "panel_name": "COMPREHENSIVE METABOLIC PANEL",
            "results": [
                {"test": "Glucose", "value": "88", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "BUN", "value": "14", "unit": "mg/dL", "ref_range": "7-20", "flag": ""},
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "eGFR", "value": ">90", "unit": "mL/min", "ref_range": ">60", "flag": ""},
                {"test": "Albumin", "value": "4.0", "unit": "g/dL", "ref_range": "3.5-5.5", "flag": ""},
            ],
        }, {
            "panel_name": "THYROID / IRON / VITAMIN D",
            "results": [
                {"test": "TSH", "value": "2.1", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Iron", "value": "48", "unit": "ug/dL", "ref_range": "60-170", "flag": "L"},
                {"test": "Ferritin", "value": "18", "unit": "ng/mL", "ref_range": "12-150", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "22", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: PCP -- October 2022 (persistent symptoms, ANA ordered)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-10-03_PCP_Follow_Up.pdf"),
        PATIENT, "10/03/2022", "pcp",
        "OFFICE VISIT -- FOLLOW-UP FATIGUE AND JOINT PAIN",
        {
            "chief_complaint": "Follow-up fatigue and joint symptoms. Symptoms persisting despite iron supplementation.",
            "hpi": (
                "Ms. Mitchell returns for follow-up. She has been taking iron supplements and vitamin D "
                "as directed. Fatigue has improved slightly but remains significant. Joint pain in hands "
                "and wrists continues, with morning stiffness now lasting approximately 60 minutes. "
                "She reports a new symptom: intermittent facial flushing and mild rash across her cheeks "
                "after sun exposure last month, which she attributed to sunburn. The rash resolved after "
                "a few days indoors. She also notes her hair seems to be thinning, though she is unsure "
                "if this is related. No fevers, weight loss, or new symptoms otherwise."
            ),
            "vitals": "BP 122/76  HR 82  Temp 98.6F  SpO2 98%  Wt 129 lbs",
            "exam": (
                "General: Appears fatigued but in no acute distress.\n"
                "HEENT: No active facial rash today. Mild diffuse hair thinning noted. No oral ulcers.\n"
                "Musculoskeletal: Bilateral MCP and PIP joints mildly tender. Possible mild swelling at "
                "right 2nd and 3rd MCP joints. Full ROM maintained.\n"
                "Skin: No active rash. No discoid lesions."
            ),
            "assessment": (
                "1. Persistent arthralgias with morning stiffness >30 min -- concerning for inflammatory etiology\n"
                "2. Photosensitive facial rash -- resolved, will monitor\n"
                "3. Mild alopecia -- differential includes telogen effluvium vs autoimmune\n"
                "4. Iron deficiency -- improving with supplementation"
            ),
            "plan": (
                "1. Check CRP, ESR, ANA with reflex, CBC\n"
                "2. Continue iron supplementation\n"
                "3. Avoid excessive sun exposure, use SPF 50+\n"
                "4. If ANA positive or inflammatory markers elevated, will refer to Rheumatology\n"
                "5. Return in 4-6 weeks for results review"
            ),
        },
    )

    # Lab: October 2022 -- ANA weakly positive, CRP mildly elevated
    generate_lab_report(
        os.path.join(output_dir, "2022-10-03_Lab_ANA_Inflammatory.pdf"),
        PATIENT, "10/03/2022", PROVIDERS["pcp"],
        [{
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP, High Sensitivity", "value": "8.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "28", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }, {
            "panel_name": "AUTOIMMUNE SCREENING",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:160", "unit": "", "ref_range": "<1:40", "flag": "A"},
                {"test": "ANA Pattern", "value": "Homogeneous", "unit": "", "ref_range": "N/A", "flag": ""},
            ],
        }, {
            "panel_name": "COMPLETE BLOOD COUNT",
            "results": [
                {"test": "WBC", "value": "4.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "12.1", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Lymphocytes", "value": "1.2", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": ""},
                {"test": "Platelets", "value": "185", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: Dermatology -- January 2023 (facial rash recurrence)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-01-18_Derm_Consult.pdf"),
        PATIENT, "01/18/2023", "dermatology",
        "DERMATOLOGY CONSULTATION -- FACIAL RASH",
        {
            "chief_complaint": "Recurrent facial rash, referred by PCP.",
            "hpi": (
                "Ms. Mitchell is a 31-year-old woman referred for evaluation of recurrent facial rash. "
                "She reports 3 episodes over the past 6 months of erythematous rash across both cheeks "
                "and the bridge of the nose, typically following sun exposure. Each episode lasts 3-7 days "
                "and resolves with sun avoidance. She uses SPF 50 sunscreen inconsistently. She also "
                "reports ongoing hair thinning. No oral ulcers, no discoid lesions elsewhere. "
                "PCP noted ANA positive at 1:160 homogeneous pattern in October 2022. "
                "No prior dermatologic history. No family history of lupus or psoriasis."
            ),
            "exam": (
                "Face: Faint erythema across bilateral malar eminences with sparing of nasolabral folds. "
                "No scaling. No discoid lesions.\n"
                "Scalp: Mild diffuse thinning, no scarring alopecia, no discoid plaques.\n"
                "Oral: No ulcers.\n"
                "Trunk/Extremities: No rashes, no livedo reticularis.\n"
                "Nails: Normal. No splinter hemorrhages."
            ),
            "assessment": (
                "1. Malar erythema -- differential includes photosensitive dermatitis, rosacea, "
                "vs malar rash of SLE. Distribution pattern with nasolabial fold sparing is notable.\n"
                "2. Diffuse non-scarring alopecia -- telogen effluvium vs lupus hair loss\n"
                "3. ANA positive 1:160 -- noted, but low titer, can be seen in healthy individuals"
            ),
            "plan": (
                "1. Strict photoprotection: SPF 50+ daily, UPF clothing, avoid midday sun\n"
                "2. Topical tacrolimus 0.1% ointment to affected areas BID for flares\n"
                "3. Consider skin biopsy if rash persists or worsens\n"
                "4. I note the positive ANA but the titer is low and nonspecific. I will not pursue "
                "further autoimmune workup from dermatology at this time. If the PCP or rheumatology "
                "wishes to evaluate further, that would be appropriate.\n"
                "5. Follow-up in 3 months"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: PCP -- April 2023 (worsening, refers to rheumatology)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-04-10_PCP_Follow_Up.pdf"),
        PATIENT, "04/10/2023", "pcp",
        "OFFICE VISIT -- FOLLOW-UP JOINT PAIN, FATIGUE",
        {
            "chief_complaint": "Worsening joint pain and fatigue. New pleuritic chest pain.",
            "hpi": (
                "Ms. Mitchell returns with worsening symptoms. Joint pain now involves wrists, MCPs, "
                "and knees bilaterally. Morning stiffness lasts 1-2 hours. She reports a new symptom: "
                "sharp chest pain that worsens with deep breathing, present intermittently for 2 weeks. "
                "No fever. Fatigue remains debilitating -- she has reduced her work hours. Hair loss "
                "continues. No rash currently. She saw Dermatology in January who prescribed topical "
                "tacrolimus for the facial rash. She reports 2 small oral ulcers last month that resolved "
                "spontaneously."
            ),
            "vitals": "BP 126/78  HR 88  Temp 99.1F  SpO2 97%  Wt 126 lbs",
            "exam": (
                "General: Fatigued-appearing. Low-grade temperature.\n"
                "HEENT: No active facial rash. No oral ulcers today. Diffuse hair thinning.\n"
                "Cardiovascular: Tachycardic. Faint pericardial friction rub vs. pleural rub heard.\n"
                "Pulmonary: Diminished breath sounds at left base. Pleuritic pain with deep inspiration.\n"
                "Musculoskeletal: Bilateral MCP, PIP, wrist tenderness with mild synovitis. "
                "Bilateral knee small effusions. No deformity."
            ),
            "assessment": (
                "1. Polyarthritis with synovitis -- now clearly inflammatory pattern\n"
                "2. Pleuritis/possible pericarditis -- concerning for serositis\n"
                "3. ANA positive (1:160, homogeneous) from October\n"
                "4. Oral ulcers, photosensitive rash, alopecia, fatigue\n"
                "5. This constellation is concerning. Differential includes SLE, RA, MCTD."
            ),
            "plan": (
                "1. URGENT rheumatology referral -- called Dr. Chen's office, appointment in 2 weeks\n"
                "2. Repeat ANA, add anti-dsDNA, anti-Smith, complement C3/C4, CBC, CMP, urinalysis\n"
                "3. Chest X-ray today to evaluate pleural/pericardial process\n"
                "4. Prednisone 20 mg daily x 7 days for acute symptom relief while awaiting rheum eval\n"
                "5. If chest pain worsens or fever develops, go to ED"
            ),
        },
    )

    # Lab: April 2023 -- anti-dsDNA borderline, complement starting to fall
    generate_lab_report(
        os.path.join(output_dir, "2023-04-10_Lab_Autoimmune_Panel.pdf"),
        PATIENT, "04/10/2023", PROVIDERS["pcp"],
        [{
            "panel_name": "AUTOIMMUNE PANEL",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:320", "unit": "", "ref_range": "<1:40", "flag": "A"},
                {"test": "ANA Pattern", "value": "Homogeneous", "unit": "", "ref_range": "N/A", "flag": ""},
                {"test": "Anti-dsDNA Ab (ELISA)", "value": "28", "unit": "IU/mL", "ref_range": "<25", "flag": "H"},
                {"test": "Anti-Smith Ab", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-SSA/Ro", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-SSB/La", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }, {
            "panel_name": "COMPLEMENT LEVELS",
            "results": [
                {"test": "Complement C3", "value": "82", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "16", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
            ],
        }, {
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP", "value": "14.6", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "42", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }, {
            "panel_name": "URINALYSIS",
            "results": [
                {"test": "Protein", "value": "Trace", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Blood", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "WBC", "value": "0-2", "unit": "/HPF", "ref_range": "0-5", "flag": ""},
                {"test": "RBC", "value": "0-2", "unit": "/HPF", "ref_range": "0-3", "flag": ""},
            ],
        }],
    )

    # Imaging: Chest X-ray -- April 2023
    generate_imaging_report(
        os.path.join(output_dir, "2023-04-10_CXR.pdf"),
        PATIENT, "04/10/2023", "pcp",
        "X-Ray", "Chest PA and Lateral",
        "Pleuritic chest pain, evaluate for effusion or pericardial disease.",
        "PA and lateral views of the chest were obtained.",
        (
            "Heart size is at the upper limits of normal. No focal consolidation. "
            "Small left-sided pleural effusion is noted. No pneumothorax. "
            "Osseous structures are unremarkable. No mediastinal widening."
        ),
        (
            "1. Small left pleural effusion.\n"
            "2. Borderline cardiomegaly -- recommend echocardiogram if clinically indicated.\n"
            "3. No acute pulmonary parenchymal disease."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: Rheumatology -- May 2023 (first rheum visit)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-05-01_Rheum_New_Patient.pdf"),
        PATIENT, "05/01/2023", "rheumatology",
        "RHEUMATOLOGY NEW PATIENT CONSULTATION",
        {
            "chief_complaint": "Referred by PCP for polyarthritis, positive ANA, pleuritis, suspected autoimmune disease.",
            "hpi": (
                "Ms. Mitchell is a 31-year-old woman referred for evaluation of possible systemic autoimmune "
                "disease. History is notable for: (1) polyarthritis involving MCPs, PIPs, wrists, and knees "
                "with inflammatory morning stiffness >1 hour x 10 months; (2) photosensitive malar rash x 3 "
                "episodes; (3) pleuritic chest pain with small pleural effusion on CXR; (4) ANA positive 1:320 "
                "homogeneous; (5) anti-dsDNA borderline positive at 28 IU/mL; (6) C3 slightly low at 82; "
                "(7) oral ulcers x1 episode; (8) diffuse alopecia; (9) fatigue. Symptoms responded partially "
                "to prednisone 20 mg taper from PCP. Currently on prednisone 10 mg daily.\n\n"
                "PMH: Iron deficiency (treated). No prior autoimmune diagnosis.\n"
                "FH: Mother with hypothyroidism. Maternal aunt with 'arthritis' (type unknown).\n"
                "Social: Software engineer. Non-smoker. Occasional alcohol. No recreational drugs."
            ),
            "vitals": "BP 124/76  HR 80  Temp 98.8F  SpO2 98%  Wt 127 lbs",
            "exam": (
                "General: Well-appearing, alert. No acute distress.\n"
                "HEENT: No active malar rash today (on prednisone). Mild diffuse alopecia. No oral ulcers. "
                "No lymphadenopathy.\n"
                "Cardiovascular: RRR. No friction rub. No murmurs.\n"
                "Pulmonary: Clear bilaterally. No pleural rub.\n"
                "Musculoskeletal: Mild synovitis at bilateral 2nd/3rd MCP joints. Wrists tender bilaterally. "
                "No knee effusions currently (improved on prednisone). No deformities. Grip strength mildly "
                "reduced bilaterally.\n"
                "Skin: No active rash. No livedo reticularis. No Raynaud's changes.\n"
                "Neurological: Intact. No focal deficits."
            ),
            "assessment": (
                "31-year-old woman with multisystem inflammatory disease. Features include inflammatory "
                "polyarthritis, photosensitive malar rash, serositis (pleuritis with effusion), oral ulcers, "
                "non-scarring alopecia, ANA 1:320 homogeneous, borderline anti-dsDNA, and low C3.\n\n"
                "Differential diagnosis:\n"
                "1. Systemic lupus erythematosus -- highest on differential given constellation\n"
                "2. Mixed connective tissue disease -- would need anti-U1-RNP\n"
                "3. Undifferentiated connective tissue disease\n"
                "4. Rheumatoid arthritis -- less likely given extra-articular features\n\n"
                "2019 ACR/EULAR SLE classification score (preliminary):\n"
                "- ANA >= 1:80: YES (entry criterion met)\n"
                "- Arthritis (synovitis in 2+ joints): 6 points\n"
                "- Pleural effusion (serositis): 5 points\n"
                "- Oral ulcers: 2 points\n"
                "- Non-scarring alopecia: 2 points\n"
                "- Anti-dsDNA: borderline, need to confirm -- if positive: 6 points\n"
                "- Low C3: 3 points\n"
                "Total (without anti-dsDNA): 18 points -- WELL ABOVE 10 point threshold\n"
                "Assessment: Highly consistent with SLE"
            ),
            "plan": (
                "1. Confirm anti-dsDNA with Farr assay (more specific than ELISA)\n"
                "2. Complete serologic workup: anti-Smith, anti-U1-RNP, anti-SSA/Ro, anti-SSB/La, "
                "antiphospholipid antibodies (anticardiolipin, lupus anticoagulant, anti-beta2-glycoprotein I)\n"
                "3. Complement C3, C4 repeat\n"
                "4. CBC with differential, CMP, urinalysis with microscopy, spot urine protein/creatinine ratio\n"
                "5. Consider HLA typing if genetic risk assessment desired\n"
                "6. Start hydroxychloroquine 200 mg BID (after ophthalmology baseline exam)\n"
                "7. Continue prednisone 10 mg, plan to taper once hydroxychloroquine therapeutic (6-8 weeks)\n"
                "8. Referral to Ophthalmology for baseline retinal exam (pre-HCQ)\n"
                "9. Discuss sun protection, fatigue management, disease education\n"
                "10. Follow-up in 6 weeks with labs"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Lab: May 2023 -- rheumatology comprehensive panel
    # ─────────────────────────────────────────────────────────────
    generate_lab_report(
        os.path.join(output_dir, "2023-05-01_Lab_Rheum_Comprehensive.pdf"),
        PATIENT, "05/01/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "AUTOIMMUNE COMPREHENSIVE PANEL",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:640", "unit": "", "ref_range": "<1:40", "flag": "A"},
                {"test": "ANA Pattern", "value": "Homogeneous", "unit": "", "ref_range": "N/A", "flag": ""},
                {"test": "Anti-dsDNA (Farr)", "value": "45", "unit": "IU/mL", "ref_range": "<10", "flag": "HH"},
                {"test": "Anti-Smith Ab", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Anti-U1-RNP", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-SSA/Ro", "value": "Weak Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Anti-SSB/La", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "RF", "value": "12", "unit": "IU/mL", "ref_range": "<14", "flag": ""},
                {"test": "Anti-CCP", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }, {
            "panel_name": "COMPLEMENT AND INFLAMMATORY",
            "results": [
                {"test": "Complement C3", "value": "72", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "11", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
                {"test": "CRP", "value": "11.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "38", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }, {
            "panel_name": "ANTIPHOSPHOLIPID ANTIBODIES",
            "results": [
                {"test": "Anticardiolipin IgG", "value": "8", "unit": "GPL", "ref_range": "<15", "flag": ""},
                {"test": "Anticardiolipin IgM", "value": "6", "unit": "MPL", "ref_range": "<15", "flag": ""},
                {"test": "Lupus Anticoagulant", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-B2-Glycoprotein I IgG", "value": "4", "unit": "SGU", "ref_range": "<20", "flag": ""},
            ],
        }, {
            "panel_name": "RENAL FUNCTION / URINALYSIS",
            "results": [
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "eGFR", "value": ">90", "unit": "mL/min", "ref_range": ">60", "flag": ""},
                {"test": "Urine Protein/Creatinine", "value": "0.18", "unit": "mg/mg", "ref_range": "<0.15", "flag": "H"},
                {"test": "Urinalysis Protein", "value": "Trace", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Urinalysis RBC", "value": "3-5", "unit": "/HPF", "ref_range": "0-3", "flag": "A"},
            ],
        }, {
            "panel_name": "HEMATOLOGY",
            "results": [
                {"test": "WBC", "value": "4.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Lymphocytes", "value": "0.9", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": "L"},
                {"test": "Hemoglobin", "value": "11.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": "L"},
                {"test": "Platelets", "value": "162", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Nephrology referral -- September 2023 (proteinuria worsening)
    # ─────────────────────────────────────────────────────────────
    generate_lab_report(
        os.path.join(output_dir, "2023-08-14_Lab_Follow_Up.pdf"),
        PATIENT, "08/14/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "LUPUS MONITORING PANEL",
            "results": [
                {"test": "Anti-dsDNA (Farr)", "value": "78", "unit": "IU/mL", "ref_range": "<10", "flag": "HH"},
                {"test": "Complement C3", "value": "58", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "8", "unit": "mg/dL", "ref_range": "10-40", "flag": "L"},
                {"test": "CRP", "value": "18.4", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "52", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }, {
            "panel_name": "RENAL FUNCTION",
            "results": [
                {"test": "Creatinine", "value": "1.0", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "eGFR", "value": "85", "unit": "mL/min", "ref_range": ">60", "flag": ""},
                {"test": "Urine Protein/Creatinine", "value": "0.62", "unit": "mg/mg", "ref_range": "<0.15", "flag": "HH"},
                {"test": "Urinalysis RBC", "value": "8-12", "unit": "/HPF", "ref_range": "0-3", "flag": "H"},
                {"test": "Urinalysis WBC", "value": "5-8", "unit": "/HPF", "ref_range": "0-5", "flag": "H"},
                {"test": "Urine RBC Casts", "value": "Present", "unit": "", "ref_range": "Absent", "flag": "A"},
            ],
        }, {
            "panel_name": "HEMATOLOGY",
            "results": [
                {"test": "WBC", "value": "3.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": "L"},
                {"test": "Lymphocytes", "value": "0.7", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": "L"},
                {"test": "Hemoglobin", "value": "10.8", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": "L"},
                {"test": "Platelets", "value": "138", "unit": "K/uL", "ref_range": "150-400", "flag": "L"},
            ],
        }],
    )

    generate_progress_note(
        os.path.join(output_dir, "2023-09-18_Nephrology_Consult.pdf"),
        PATIENT, "09/18/2023", "nephrology",
        "NEPHROLOGY CONSULTATION -- SUSPECTED LUPUS NEPHRITIS",
        {
            "chief_complaint": "Referred by Rheumatology for worsening proteinuria and active urine sediment in setting of SLE.",
            "hpi": (
                "Ms. Mitchell is a 32-year-old woman with recently diagnosed SLE (diagnosed May 2023 by "
                "Dr. Chen, Rheumatology) now with worsening renal parameters. Current labs show urine "
                "protein/creatinine ratio 0.62 mg/mg (up from 0.18 in May), active urine sediment with "
                "RBC casts, rising anti-dsDNA (78 IU/mL, up from 45), and falling complements (C3 58, C4 8). "
                "Creatinine remains 1.0 mg/dL (eGFR 85). She is currently on hydroxychloroquine 200 mg BID "
                "and prednisone 10 mg daily.\n\n"
                "She reports mild bilateral lower extremity edema over the past 2 weeks and foamy urine. "
                "No dysuria, no gross hematuria. Blood pressure has been running higher (130s/80s at home)."
            ),
            "vitals": "BP 138/84  HR 82  Temp 98.6F  SpO2 98%  Wt 131 lbs (up 4 lbs from baseline)",
            "exam": (
                "General: No acute distress. Mild facial puffiness.\n"
                "Cardiovascular: RRR, no murmurs. No JVD.\n"
                "Extremities: 1+ bilateral pretibial edema.\n"
                "Abdomen: Soft, non-tender, no organomegaly."
            ),
            "assessment": (
                "SLE with probable lupus nephritis. Active serologic markers (rising anti-dsDNA, "
                "hypocomplementemia) with progressive proteinuria and active urine sediment including "
                "RBC casts are highly suggestive of proliferative lupus nephritis (ISN/RPS Class III or IV).\n\n"
                "Renal biopsy is indicated to confirm diagnosis and guide therapy."
            ),
            "plan": (
                "1. Schedule renal biopsy within 2 weeks\n"
                "2. 24-hour urine collection for protein quantification\n"
                "3. Increase prednisone to 40 mg daily pending biopsy results\n"
                "4. Start lisinopril 10 mg daily for proteinuria and BP control\n"
                "5. Continue hydroxychloroquine\n"
                "6. Renal biopsy results will guide immunosuppressive therapy "
                "(mycophenolate vs cyclophosphamide per ACR 2024 lupus nephritis guidelines)\n"
                "7. Follow-up 1 week post-biopsy with results"
            ),
        },
    )

    # Pathology: Renal biopsy -- October 2023
    generate_pathology_report(
        os.path.join(output_dir, "2023-10-02_Path_Renal_Biopsy.pdf"),
        PATIENT, "10/02/2023",
        "Kidney, left, percutaneous needle biopsy",
        "32-year-old woman with SLE, worsening proteinuria, active urine sediment, rising anti-dsDNA, low complements. "
        "Rule out lupus nephritis.",
        (
            "Received in formalin: two cores of renal tissue measuring 1.2 cm and 0.8 cm in greatest "
            "dimension. Tissue submitted entirely for light microscopy, immunofluorescence, and electron "
            "microscopy."
        ),
        (
            "Light Microscopy: Twenty-two glomeruli are present, none globally sclerosed. Fourteen of 22 "
            "glomeruli (64%) show diffuse endocapillary hypercellularity with segmental areas of "
            "mesangial proliferation. Wire loop deposits are identified in 6 glomeruli. Focal fibrinoid "
            "necrosis is present in 2 glomeruli. No crescents are identified. Tubular atrophy and "
            "interstitial fibrosis are minimal (<5%). Arteries and arterioles show no vasculitis.\n\n"
            "Immunofluorescence: Full house pattern with granular staining for IgG (3+), IgA (2+), "
            "IgM (2+), C3 (3+), C4 (2+), C1q (3+), kappa (2+), lambda (2+) along glomerular "
            "capillary walls and mesangium.\n\n"
            "Electron Microscopy: Subendothelial, subepithelial, and mesangial electron-dense deposits. "
            "Tubuloreticular inclusions present in endothelial cells. Foot process effacement "
            "approximately 40%."
        ),
        (
            "LUPUS NEPHRITIS, ISN/RPS CLASS IV-G (A)\n"
            "- Diffuse global proliferative lupus nephritis, active\n"
            "- Activity Index: 8/24 (endocapillary proliferation, wire loops, fibrinoid necrosis)\n"
            "- Chronicity Index: 1/12 (minimal tubular atrophy/interstitial fibrosis)\n"
            "- Full house immunofluorescence pattern consistent with lupus nephritis\n"
            "- Tubuloreticular inclusions (interferon signature)"
        ),
        comment=(
            "The biopsy findings are diagnostic of Class IV-G (A) diffuse proliferative lupus nephritis "
            "per the ISN/RPS 2003 classification. The relatively low chronicity index is favorable. "
            "Induction immunosuppressive therapy is recommended per current guidelines."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Genetic testing -- November 2023
    # ─────────────────────────────────────────────────────────────
    generate_genetic_report(
        os.path.join(output_dir, "2023-11-15_Genetics_HLA_Report.pdf"),
        PATIENT, "11/15/2023",
        "HLA CLASS I AND CLASS II TYPING",
        "Systemic lupus erythematosus with Class IV lupus nephritis. HLA typing for disease association "
        "and pharmacogenomic assessment.",
        "Next-generation sequencing (NGS) of HLA Class I (A, B, C) and Class II (DRB1, DQB1, DPB1) loci. "
        "Resolution: high-resolution (2-field).",
        [
            {"Locus": "HLA-A", "Allele 1": "A*02:01", "Allele 2": "A*24:02"},
            {"Locus": "HLA-B", "Allele 1": "B*08:01", "Allele 2": "B*44:02"},
            {"Locus": "HLA-C", "Allele 1": "C*05:01", "Allele 2": "C*07:01"},
            {"Locus": "HLA-DRB1", "Allele 1": "DRB1*03:01", "Allele 2": "DRB1*07:01"},
            {"Locus": "HLA-DQB1", "Allele 1": "DQB1*02:01", "Allele 2": "DQB1*02:02"},
        ],
        (
            "HLA-DRB1*03:01 is associated with increased susceptibility to systemic lupus erythematosus "
            "(OR 2.4, PMID:19864127), Sjogren's syndrome (OR 3.1), type 1 diabetes (OR 3.6), and "
            "celiac disease (OR 7.0). This allele is one of the strongest HLA risk factors for SLE in "
            "European-descent populations.\n\n"
            "HLA-DQB1*02:01 (part of the HLA-DQ2 heterodimer with DQA1*05:01) is associated with "
            "celiac disease (OR 7.0) and type 1 diabetes (OR 3.0).\n\n"
            "HLA-B*08:01 is part of the A1-B8-DR3 ancestral haplotype (8.1 haplotype), which carries "
            "broad autoimmune susceptibility.\n\n"
            "Pharmacogenomic note: HLA-DRB1*03:01 has been associated with increased risk of "
            "anti-drug antibody formation to adalimumab and other TNF inhibitors."
        ),
        recommendations=(
            "1. The HLA profile confirms strong genetic susceptibility to SLE, consistent with clinical diagnosis.\n"
            "2. Given HLA-DQ2 positivity, consider screening for celiac disease (anti-tTG IgA) if GI symptoms develop.\n"
            "3. Monitor for other autoimmune conditions given broad autoimmune susceptibility haplotype.\n"
            "4. If TNF inhibitor therapy is considered for any indication, monitor for anti-drug antibodies given DRB1*03:01."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Later visits: Dec 2023 - Dec 2025 (treatment, monitoring)
    # ─────────────────────────────────────────────────────────────

    # Lab: December 2023 -- post-induction (improving)
    generate_lab_report(
        os.path.join(output_dir, "2023-12-11_Lab_Post_Induction.pdf"),
        PATIENT, "12/11/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "LUPUS MONITORING",
            "results": [
                {"test": "Anti-dsDNA (Farr)", "value": "52", "unit": "IU/mL", "ref_range": "<10", "flag": "HH"},
                {"test": "Complement C3", "value": "68", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "10", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
                {"test": "CRP", "value": "6.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "Urine Protein/Creatinine", "value": "0.38", "unit": "mg/mg", "ref_range": "<0.15", "flag": "H"},
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
            ],
        }],
    )

    # Lab: June 2024 -- partial remission
    generate_lab_report(
        os.path.join(output_dir, "2024-06-17_Lab_Maintenance.pdf"),
        PATIENT, "06/17/2024", PROVIDERS["rheumatology"],
        [{
            "panel_name": "LUPUS MONITORING",
            "results": [
                {"test": "Anti-dsDNA (Farr)", "value": "22", "unit": "IU/mL", "ref_range": "<10", "flag": "H"},
                {"test": "Complement C3", "value": "88", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "15", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
                {"test": "CRP", "value": "3.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "22", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "Urine Protein/Creatinine", "value": "0.14", "unit": "mg/mg", "ref_range": "<0.15", "flag": ""},
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
            ],
        }, {
            "panel_name": "HEMATOLOGY",
            "results": [
                {"test": "WBC", "value": "5.4", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Lymphocytes", "value": "1.1", "unit": "K/uL", "ref_range": "1.0-4.0", "flag": ""},
                {"test": "Hemoglobin", "value": "12.2", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "178", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # Rheum follow-up: December 2024
    generate_progress_note(
        os.path.join(output_dir, "2024-12-09_Rheum_Follow_Up.pdf"),
        PATIENT, "12/09/2024", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- SLE / LUPUS NEPHRITIS",
        {
            "chief_complaint": "Routine follow-up SLE with lupus nephritis Class IV, on maintenance therapy.",
            "hpi": (
                "Ms. Mitchell is doing well overall. She completed induction therapy with mycophenolate "
                "mofetil 3g/day and prednisone taper. Currently on maintenance mycophenolate 2g/day, "
                "hydroxychloroquine 200 mg BID, and prednisone 5 mg daily. Proteinuria has normalized. "
                "Anti-dsDNA trending down. Complements improving. No arthritis flares. No rash. Energy "
                "much improved. She has returned to full-time work.\n\n"
                "Only concern: mild hair thinning persists, likely medication-related. No new symptoms."
            ),
            "vitals": "BP 118/72  HR 74  Temp 98.4F  SpO2 99%  Wt 130 lbs",
            "exam": (
                "General: Well-appearing. Good energy.\n"
                "HEENT: No rash, no oral ulcers, mild diffuse hair thinning (stable).\n"
                "Musculoskeletal: No synovitis. Full ROM. No tenderness.\n"
                "Extremities: No edema.\n"
                "Skin: No rash. No livedo."
            ),
            "assessment": (
                "SLE with Class IV lupus nephritis -- currently in partial clinical remission on maintenance "
                "immunosuppression. Serologic markers continue to improve. Proteinuria resolved."
            ),
            "plan": (
                "1. Continue current regimen: MMF 2g/day, HCQ 200 mg BID, prednisone 5 mg daily\n"
                "2. Consider adding belimumab if unable to taper prednisone below 5 mg\n"
                "3. Labs today: anti-dsDNA, C3/C4, CBC, CMP, urine P/C ratio\n"
                "4. Annual ophthalmology exam for HCQ monitoring\n"
                "5. Follow-up in 3 months\n"
                "6. Discussed pregnancy planning -- will need to switch from MMF to azathioprine if planning conception"
            ),
        },
    )

    # Lab: December 2024 -- near remission
    generate_lab_report(
        os.path.join(output_dir, "2024-12-09_Lab_Monitoring.pdf"),
        PATIENT, "12/09/2024", PROVIDERS["rheumatology"],
        [{
            "panel_name": "LUPUS MONITORING",
            "results": [
                {"test": "Anti-dsDNA (Farr)", "value": "14", "unit": "IU/mL", "ref_range": "<10", "flag": "H"},
                {"test": "Complement C3", "value": "95", "unit": "mg/dL", "ref_range": "90-180", "flag": ""},
                {"test": "Complement C4", "value": "18", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
                {"test": "CRP", "value": "2.1", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "14", "unit": "mm/hr", "ref_range": "0-20", "flag": ""},
                {"test": "Urine Protein/Creatinine", "value": "0.08", "unit": "mg/mg", "ref_range": "<0.15", "flag": ""},
            ],
        }],
    )

    # Final lab: flare warning -- late 2025
    generate_lab_report(
        os.path.join(output_dir, "2025-09-22_Lab_Flare_Warning.pdf"),
        PATIENT, "09/22/2025", PROVIDERS["rheumatology"],
        [{
            "panel_name": "LUPUS MONITORING -- ROUTINE",
            "results": [
                {"test": "Anti-dsDNA (Farr)", "value": "38", "unit": "IU/mL", "ref_range": "<10", "flag": "HH"},
                {"test": "Complement C3", "value": "74", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "Complement C4", "value": "9", "unit": "mg/dL", "ref_range": "10-40", "flag": "L"},
                {"test": "CRP", "value": "8.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "Urine Protein/Creatinine", "value": "0.22", "unit": "mg/mg", "ref_range": "<0.15", "flag": "H"},
            ],
        }],
    )

    print(f"  Patient Sarah Mitchell: {len(os.listdir(output_dir))} documents generated")
