"""Patient: Karen Foster -- 48F, Systemic Sclerosis (Scleroderma).

Timeline: January 2023 -> January 2026
Specialists seen: PCP, Rheumatology, Dermatology, Pulmonology, Cardiology
Key pattern: Raynaud's -> skin thickening -> anti-Scl-70 positive ->
             ILD on PFTs -> PAH screening -> mycophenolate + nintedanib

The system should detect: Diffuse cutaneous systemic sclerosis (dcSSc),
                           anti-Scl-70 positive, ILD, early PAH risk.
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
    generate_lab_report,
    generate_pathology_report,
    generate_progress_note,
)

PATIENT = {
    "name": "Karen Foster",
    "dob": "1977-06-10",
    "mrn": "KFO-2023-38216",
    "sex": "F",
    "age_at_start": 45,
}


def generate(output_dir: str):
    """Generate all clinical documents for Karen Foster."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: PCP -- January 2023 (Raynaud's, puffy fingers)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-01-18_PCP_Progress_Note.pdf"),
        PATIENT, "01/18/2023", "pcp",
        "OFFICE VISIT -- COLD-SENSITIVITY AND FINGER SWELLING",
        {
            "chief_complaint": "Color changes in fingers with cold exposure x 6 months. Swollen, puffy fingers x 3 months. Fatigue.",
            "hpi": (
                "Mrs. Foster is a 45-year-old woman presenting with 6-month history of episodic "
                "color changes in her fingers triggered by cold temperatures or emotional stress. "
                "She describes a triphasic pattern: fingers turn white, then blue/purple, then red "
                "with painful throbbing upon rewarming. Episodes last 15-30 minutes and occur "
                "multiple times daily in winter. She has also noticed bilateral finger swelling "
                "(puffy fingers) over the past 3 months, making it difficult to remove her rings. "
                "She reports progressive fatigue, mild arthralgia in her hands, and new difficulty "
                "swallowing solid foods (food feels like it 'sticks' in her mid-chest). No rash "
                "or skin changes noticed by patient. No shortness of breath. No muscle weakness."
            ),
            "vitals": "BP 128/82  HR 76  Temp 98.2F  SpO2 98%  Wt 148 lbs  Ht 5'4\"",
            "ros": (
                "Constitutional: Fatigue x 3 months. No fever or weight loss.\n"
                "HEENT: Dry mouth. No oral ulcers.\n"
                "Cardiovascular: No chest pain, palpitations, or edema.\n"
                "Pulmonary: No cough or dyspnea.\n"
                "GI: Dysphagia to solids as above. Occasional heartburn. No diarrhea.\n"
                "Musculoskeletal: Bilateral hand stiffness and swelling. Mild arthralgias.\n"
                "Skin: Color changes in digits as above. No rash reported.\n"
                "Neurological: No numbness or weakness."
            ),
            "exam": (
                "General: Well-appearing woman, no acute distress.\n"
                "HEENT: Mild perioral skin tightening (subtle). Decreased oral aperture (3.5 cm, "
                "borderline reduced). Dry oral mucosa.\n"
                "Cardiovascular: RRR, no murmurs. No JVD. No peripheral edema.\n"
                "Pulmonary: Clear bilaterally. No crackles.\n"
                "Extremities: Bilateral diffuse finger swelling ('puffy fingers'). Fingers appear "
                "sausage-like. Mild skin induration over dorsum of hands extending to wrists. "
                "No digital ulcers. No calcinosis. Nail fold examination: dilated capillary loops "
                "visible with ophthalmoscope (abnormal).\n"
                "Skin: Subtle skin tightening over forearms. No telangiectases.\n"
                "Musculoskeletal: Mild tenderness MCPs bilaterally. Full ROM but limited finger "
                "flexion due to swelling."
            ),
            "assessment": (
                "1. Raynaud phenomenon, severe -- triphasic color changes, multiple daily episodes. "
                "New onset at age 45 with associated puffy fingers is concerning for SECONDARY "
                "Raynaud's (connective tissue disease).\n"
                "2. Puffy fingers (sclerodactyly precursor) -- high concern for early systemic sclerosis\n"
                "3. Dysphagia -- could indicate esophageal dysmotility (scleroderma GI involvement)\n"
                "4. Subtle perioral tightening and forearm skin changes\n"
                "5. Abnormal nailfold capillaries on bedside exam\n\n"
                "CLINICAL CONCERN: This presentation is highly suspicious for EARLY SYSTEMIC "
                "SCLEROSIS. Urgent rheumatology referral indicated."
            ),
            "plan": (
                "1. URGENT rheumatology referral\n"
                "2. ANA panel, anti-Scl-70, anti-centromere, anti-RNA Polymerase III\n"
                "3. CBC, CMP, ESR, CRP\n"
                "4. Nifedipine ER 30 mg daily for Raynaud's\n"
                "5. Omeprazole 20 mg BID for dysphagia/reflux\n"
                "6. Avoid cold exposure, use heated gloves\n"
                "7. Return after rheumatology evaluation"
            ),
        },
    )

    # Lab: January 2023
    generate_lab_report(
        os.path.join(output_dir, "2023-01-18_Lab_Initial.pdf"),
        PATIENT, "01/18/2023", PROVIDERS["pcp"],
        [{
            "panel_name": "AUTOIMMUNE PANEL",
            "results": [
                {"test": "ANA", "value": "Positive (1:640)", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Pattern", "value": "Nucleolar", "unit": "", "ref_range": "N/A", "flag": ""},
                {"test": "Anti-Scl-70 (Anti-Topoisomerase I)", "value": "Positive (>200)", "unit": "AU/mL", "ref_range": "<20", "flag": "HH"},
                {"test": "Anti-Centromere Antibody", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-RNA Polymerase III", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }, {
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "ESR", "value": "32", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "CRP", "value": "8.4", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC AND CMP",
            "results": [
                {"test": "WBC", "value": "7.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "12.6", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "268", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "ALT", "value": "24", "unit": "U/L", "ref_range": "7-56", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: Rheumatology Consult -- February 2023
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-02-14_Rheumatology_Consult.pdf"),
        PATIENT, "02/14/2023", "rheumatology",
        "RHEUMATOLOGY CONSULTATION -- SUSPECTED SYSTEMIC SCLEROSIS",
        {
            "chief_complaint": (
                "Referred for evaluation of Raynaud's, puffy fingers, dysphagia, positive "
                "ANA 1:640 nucleolar pattern, and anti-Scl-70 positive."
            ),
            "hpi": (
                "Mrs. Foster is a 45-year-old woman referred by Dr. Martinez with a constellation "
                "of findings highly suspicious for systemic sclerosis: 6-month history of severe "
                "Raynaud phenomenon (triphasic), 3 months of bilateral puffy fingers with early "
                "skin induration over hands and forearms, new dysphagia to solids, and ANA 1:640 "
                "with nucleolar pattern and strongly positive anti-Scl-70 (anti-topoisomerase I). "
                "She was started on nifedipine and omeprazole by PCP."
            ),
            "exam": (
                "Skin Assessment (modified Rodnan Skin Score):\n"
                "- Fingers: 2/3 bilateral (moderate thickening)\n"
                "- Dorsum of hands: 2/3 bilateral\n"
                "- Forearms: 1/3 bilateral (early involvement)\n"
                "- Upper arms: 0/3\n"
                "- Face: 1/3 (perioral tightening, reduced oral aperture 3.2 cm)\n"
                "- Chest: 0/3\n"
                "- Abdomen: 0/3\n"
                "- Thighs/legs/feet: 0/3\n"
                "mRSS Total: 10/51\n\n"
                "Nailfold Capillaroscopy: Performed with dermatoscope.\n"
                "- Dilated giant capillaries\n"
                "- Multiple areas of capillary dropout (avascular zones)\n"
                "- Microhemorrhages\n"
                "- Pattern: ACTIVE SCLERODERMA PATTERN\n\n"
                "Hands: Puffy fingers, sclerodactyly developing. No digital pitting scars or ulcers. "
                "No calcinosis. Tendon friction rubs bilaterally at wrists (SIGNIFICANT).\n"
                "Chest: Clear auscultation. No bibasilar crackles.\n"
                "Heart: RRR. Loud P2 (note for future monitoring)."
            ),
            "assessment": (
                "DIAGNOSIS: SYSTEMIC SCLEROSIS -- Diffuse Cutaneous (dcSSc)\n\n"
                "ACR/EULAR 2013 Classification Criteria: MET (score >=9)\n"
                "- Skin thickening of fingers extending proximal to MCPs: 9 points\n"
                "- Puffy fingers: 2 points (not counted, superseded by above)\n"
                "- Raynaud phenomenon: 3 points\n"
                "- SSc-related autoantibody (anti-Scl-70): 3 points\n"
                "- Abnormal nailfold capillaries: 2 points\n"
                "Total: >=9 -- DIAGNOSTIC\n\n"
                "ANTI-SCL-70 POSITIVE dcSSc carries HIGH RISK for:\n"
                "1. Interstitial lung disease (ILD) -- 70-80% develop ILD\n"
                "2. Renal crisis (lower risk than anti-RNA Pol III, but monitor)\n"
                "3. Cardiac involvement\n"
                "4. GI dysmotility\n\n"
                "Tendon friction rubs are a poor prognostic sign in dcSSc -- predictor of "
                "rapid skin progression and internal organ involvement."
            ),
            "plan": (
                "1. BASELINE ORGAN ASSESSMENT (urgent):\n"
                "   a. High-resolution CT chest (HRCT) -- ILD screening\n"
                "   b. Pulmonary function tests with DLCO\n"
                "   c. Echocardiogram -- PAH screening and cardiac assessment\n"
                "   d. BNP/NT-proBNP\n"
                "2. Additional labs: CBC, CMP, urinalysis (renal crisis screening), CK (myopathy), "
                "BNP, anti-Pm-Scl, anti-Th/To\n"
                "3. Dermatology referral for skin biopsy (histologic confirmation)\n"
                "4. Start mycophenolate mofetil (CellCept) 500 mg BID, increase to 1500 mg BID "
                "over 4 weeks -- immunosuppression for skin and ILD prevention\n"
                "5. Continue nifedipine for Raynaud's\n"
                "6. BP monitoring at home -- scleroderma renal crisis can present as sudden HTN\n"
                "7. Follow-up in 4 weeks with imaging results"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: Dermatology Biopsy -- February 2023
    # ─────────────────────────────────────────────────────────────
    generate_pathology_report(
        os.path.join(output_dir, "2023-02-22_Pathology_Skin_Biopsy.pdf"),
        PATIENT, "02/22/2023",
        "Punch biopsy, left forearm (4 mm)",
        (
            "45-year-old female with clinically suspected diffuse cutaneous systemic sclerosis. "
            "Anti-Scl-70 positive. Skin induration over forearms and hands."
        ),
        (
            "Received in formalin: 4 mm punch biopsy of skin, 0.4 x 0.4 x 0.5 cm. "
            "Skin surface is smooth and slightly shiny. Bisected and entirely submitted."
        ),
        (
            "Sections show skin with marked dermal fibrosis. The reticular dermis is expanded "
            "with dense, homogenized collagen bundles replacing normal loose connective tissue. "
            "Collagen extends into subcutaneous fat with entrapment of fat lobules (fat trapping). "
            "Eccrine glands are atrophic and surrounded by dense collagen. Hair follicles are "
            "diminished. The epidermis is mildly atrophic with loss of rete ridges. A mild "
            "perivascular lymphocytic infiltrate is present in the superficial dermis, predominantly "
            "CD4+ T cells. No vasculitis. Trichrome stain confirms excessive collagen deposition "
            "throughout the dermis and superficial subcutis."
        ),
        (
            "SKIN, LEFT FOREARM (PUNCH BIOPSY):\n"
            "- Dermal fibrosis with collagen deposition consistent with SCLERODERMA\n"
            "- Dense homogenized collagen replacing normal dermal architecture\n"
            "- Fat trapping (collagen extension into subcutis)\n"
            "- Eccrine gland atrophy\n"
            "- Mild perivascular lymphocytic inflammation\n\n"
            "CONSISTENT WITH SYSTEMIC SCLEROSIS (SCLERODERMA)"
        ),
        (
            "Histologic findings are characteristic of scleroderma and correlate with the clinical "
            "presentation and positive anti-Scl-70 serology. The presence of fat trapping and "
            "perivascular inflammation suggests an active fibrotic process."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: HRCT Chest -- March 2023 (ILD detected)
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2023-03-06_HRCT_Chest.pdf"),
        PATIENT, "03/06/2023", "pulmonology",
        "HRCT", "Chest (High-Resolution CT)",
        "45-year-old female with newly diagnosed diffuse cutaneous systemic sclerosis, "
        "anti-Scl-70 positive. Screening for interstitial lung disease.",
        "Non-contrast high-resolution CT of the chest with 1 mm axial slices, "
        "prone and supine acquisitions.",
        (
            "Lungs:\n"
            "- Bilateral ground-glass opacities in the dependent portions of the lower lobes, "
            "predominantly basilar and posterior, persisting on prone images (not positional).\n"
            "- Early peripheral reticulation in bilateral lower lobes with subtle traction "
            "bronchiolectasis.\n"
            "- No honeycombing.\n"
            "- Pattern is NONSPECIFIC INTERSTITIAL PNEUMONIA (NSIP) pattern, typical for "
            "SSc-associated ILD.\n"
            "- Extent: approximately 15-20% of total lung volume involved.\n\n"
            "Mediastinum:\n"
            "- Mild esophageal dilation with air-fluid level, consistent with esophageal "
            "dysmotility.\n"
            "- No significant lymphadenopathy.\n"
            "- Main pulmonary artery diameter 28 mm (upper normal limit, monitor for PAH).\n\n"
            "Other:\n"
            "- No pleural effusion.\n"
            "- No pericardial effusion."
        ),
        (
            "1. INTERSTITIAL LUNG DISEASE consistent with NSIP pattern -- EXPECTED IN CONTEXT "
            "OF ANTI-SCL-70 POSITIVE SYSTEMIC SCLEROSIS.\n"
            "   - Bilateral lower lobe ground-glass opacities and early reticulation\n"
            "   - Approximately 15-20% lung involvement (moderate)\n"
            "   - No honeycombing (favorable)\n"
            "2. Esophageal dilation with air-fluid level -- esophageal dysmotility (SSc-related).\n"
            "3. Main pulmonary artery at upper normal limits -- baseline for PAH monitoring.\n\n"
            "RECOMMENDATION: PFTs with DLCO, echocardiogram for PAH screening, and immunosuppressive "
            "therapy per rheumatology."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: PFTs -- March 2023
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-03-10_Pulmonology_PFTs.pdf"),
        PATIENT, "03/10/2023", "pulmonology",
        "PULMONARY FUNCTION TEST REPORT",
        {
            "chief_complaint": "Baseline PFTs for SSc-ILD. HRCT shows NSIP pattern with 15-20% involvement.",
            "hpi": (
                "Mrs. Foster is a 45-year-old woman with newly diagnosed diffuse cutaneous "
                "systemic sclerosis (anti-Scl-70 positive). HRCT chest shows NSIP pattern ILD "
                "with bilateral lower lobe ground-glass opacities. She currently reports no "
                "significant dyspnea at rest but notes mild exercise intolerance (gets winded "
                "climbing 2 flights of stairs, previously could do 4-5 without difficulty). "
                "No cough. Non-smoker."
            ),
            "exam": (
                "SPIROMETRY:\n"
                "- FVC: 2.72 L (74% predicted) -- REDUCED\n"
                "- FEV1: 2.28 L (78% predicted)\n"
                "- FEV1/FVC ratio: 0.84 (normal -- no obstruction)\n\n"
                "LUNG VOLUMES:\n"
                "- TLC: 3.86 L (72% predicted) -- REDUCED (restrictive pattern)\n"
                "- RV: 1.14 L (84% predicted)\n\n"
                "DIFFUSION CAPACITY:\n"
                "- DLCO: 16.8 mL/min/mmHg (62% predicted) -- REDUCED\n"
                "- DLCO/VA (KCO): 3.8 mL/min/mmHg/L (78% predicted)\n\n"
                "6-MINUTE WALK TEST:\n"
                "- Distance: 420 meters (reduced for age, predicted >500 m)\n"
                "- SpO2 baseline: 97%, nadir: 92% at 4 minutes (desaturation)\n"
                "- HR baseline: 78, peak: 112\n"
                "- Borg dyspnea scale: 5/10 at completion"
            ),
            "assessment": (
                "1. RESTRICTIVE VENTILATORY DEFECT -- FVC 74%, TLC 72%\n"
                "2. REDUCED DIFFUSION CAPACITY -- DLCO 62%, disproportionate to restriction\n"
                "   - Disproportionately low DLCO raises concern for:\n"
                "     a. ILD component (confirmed on HRCT)\n"
                "     b. Possible early pulmonary vascular disease / PAH\n"
                "3. EXERTIONAL DESATURATION to 92% on 6MWT\n"
                "4. Reduced exercise capacity (420 m)\n\n"
                "In context of SSc-ILD with DLCO <70% predicted, patient should be monitored "
                "closely. DLCO <60% or FVC decline >10% absolute would warrant treatment escalation."
            ),
            "plan": (
                "1. Confirm mycophenolate started by rheumatology for ILD\n"
                "2. Consider adding nintedanib (SENSCIS trial data) if FVC declines\n"
                "3. Echocardiogram -- DLCO disproportionately low, evaluate for PAH\n"
                "4. Repeat PFTs in 3-6 months to assess trajectory\n"
                "5. Pneumococcal and influenza vaccination\n"
                "6. Follow-up with results of echocardiogram"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Echocardiogram -- March 2023
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2023-03-15_Echo_PAH_Screening.pdf"),
        PATIENT, "03/15/2023", "cardiology",
        "Echocardiogram", "Transthoracic",
        "45-year-old female with diffuse cutaneous SSc, anti-Scl-70 positive, ILD on HRCT. "
        "DLCO disproportionately reduced (62%). PAH screening.",
        "Complete 2D, M-mode, and Doppler transthoracic echocardiogram.",
        (
            "Left ventricle:\n"
            "- Normal size. LVEDD 4.6 cm. Normal wall thickness.\n"
            "- Ejection fraction: 60% (normal)\n"
            "- No regional wall motion abnormalities.\n"
            "- Grade I diastolic dysfunction (impaired relaxation).\n\n"
            "Right ventricle:\n"
            "- Mildly dilated. RVEDD 3.8 cm (upper normal 3.5 cm).\n"
            "- RV systolic function: TAPSE 2.0 cm (borderline, normal >1.7 cm).\n"
            "- Estimated RVSP: 38 mmHg (mildly elevated, normal <35 mmHg).\n\n"
            "Tricuspid valve: Mild regurgitation. TR jet velocity 2.9 m/s.\n\n"
            "Pulmonary valve: Normal. Pulmonic acceleration time 98 ms (borderline reduced).\n\n"
            "Pericardium: Trivial pericardial effusion (5 mm).\n\n"
            "IVC: Normal caliber (1.8 cm) with >50% inspiratory collapse. Estimated RAP 5 mmHg."
        ),
        (
            "1. MILDLY ELEVATED RVSP (38 mmHg) -- borderline pulmonary hypertension. In context "
            "of SSc with ILD, this could represent:\n"
            "   a. Group 1 PAH (SSc-associated) -- most concerning\n"
            "   b. Group 3 PH secondary to ILD\n"
            "2. Mild RV dilation with borderline function -- early right heart strain.\n"
            "3. Trivial pericardial effusion -- may be SSc-related.\n"
            "4. LV function preserved.\n\n"
            "RECOMMENDATION: NT-proBNP levels. Consider right heart catheterization if RVSP "
            "progresses. Annual echo screening per DETECT algorithm for SSc-PAH."
        ),
    )

    # Lab: Cardiac biomarkers
    generate_lab_report(
        os.path.join(output_dir, "2023-03-15_Lab_Cardiac.pdf"),
        PATIENT, "03/15/2023", PROVIDERS["cardiology"],
        [{
            "panel_name": "CARDIAC BIOMARKERS",
            "results": [
                {"test": "NT-proBNP", "value": "285", "unit": "pg/mL", "ref_range": "<125", "flag": "H"},
                {"test": "Troponin I", "value": "<0.01", "unit": "ng/mL", "ref_range": "<0.04", "flag": ""},
                {"test": "BNP", "value": "92", "unit": "pg/mL", "ref_range": "<100", "flag": ""},
            ],
        }, {
            "panel_name": "ADDITIONAL LABS",
            "results": [
                {"test": "CK", "value": "142", "unit": "U/L", "ref_range": "30-135", "flag": "H"},
                {"test": "Aldolase", "value": "8.2", "unit": "U/L", "ref_range": "1.0-7.5", "flag": "H"},
                {"test": "Uric Acid", "value": "7.8", "unit": "mg/dL", "ref_range": "2.4-6.0", "flag": "H"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 7: Rheumatology Follow-up -- April 2023
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-04-18_Rheumatology_Follow_Up.pdf"),
        PATIENT, "04/18/2023", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- SSc BASELINE ASSESSMENT COMPLETE",
        {
            "chief_complaint": "Follow-up dcSSc. Review baseline organ assessment results.",
            "hpi": (
                "Mrs. Foster returns after completing baseline organ assessment. She has been on "
                "mycophenolate 1000 mg BID (uptitrated from 500 mg BID over 4 weeks). Tolerating "
                "well with mild GI upset. Raynaud episodes less frequent on nifedipine. Skin "
                "thickening seems stable but she reports new tightening over anterior chest."
            ),
            "exam": (
                "Skin Assessment (mRSS re-evaluation):\n"
                "- Fingers: 2/3 bilateral\n"
                "- Hands: 2/3 bilateral\n"
                "- Forearms: 2/3 bilateral (INCREASED from 1/3)\n"
                "- Upper arms: 1/3 bilateral (NEW involvement)\n"
                "- Face: 1/3\n"
                "- Chest: 1/3 (NEW involvement)\n"
                "- Abdomen: 0/3\n"
                "- Thighs: 0/3\n"
                "- Legs/Feet: 0/3\n"
                "mRSS Total: 18/51 (INCREASED from 10 -- PROGRESSING)\n\n"
                "Digital exam: No ulcers. No calcinosis.\n"
                "Chest: Bibasilar fine crackles (new finding).\n"
                "Tendon friction rubs: Still present bilateral wrists."
            ),
            "assessment": (
                "DIFFUSE CUTANEOUS SSc -- PROGRESSING despite 8 weeks of mycophenolate\n\n"
                "1. Skin: mRSS increased 10 -> 18 (new forearm, upper arm, chest involvement)\n"
                "   - This is concerning and within the first 3 years of disease -- peak risk "
                "for progression\n"
                "2. ILD: NSIP pattern on HRCT, FVC 74%, DLCO 62%\n"
                "   - New bibasilar crackles on exam\n"
                "   - Will need repeat PFTs in 2 months to assess trajectory\n"
                "3. PAH concern: RVSP 38 mmHg, NT-proBNP 285\n"
                "   - Borderline -- continue monitoring\n"
                "4. Esophageal dysmotility: stable symptoms on PPI\n"
                "5. Elevated CK/aldolase -- possible inflammatory myopathy overlap\n"
                "6. Tendon friction rubs persist -- poor prognostic marker"
            ),
            "plan": (
                "1. Increase mycophenolate to 1500 mg BID (maximum dose)\n"
                "2. ADD nintedanib (Ofev) 150 mg BID -- SENSCIS trial showed 44% reduction in "
                "FVC decline in SSc-ILD. Start given progressing ILD on exam.\n"
                "3. Repeat PFTs with DLCO in 2 months (May 2023)\n"
                "4. Repeat HRCT chest in 6 months\n"
                "5. EMG/NCS and MRI thighs if CK rises further (myopathy workup)\n"
                "6. Continue nifedipine, omeprazole\n"
                "7. BP log -- bring to each visit (renal crisis surveillance)\n"
                "8. Annual echocardiogram for PAH surveillance\n"
                "9. Follow-up in 8 weeks"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 8: Nailfold Capillaroscopy Report -- April 2023
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-04-18_Capillaroscopy.pdf"),
        PATIENT, "04/18/2023", "rheumatology",
        "NAILFOLD CAPILLAROSCOPY REPORT",
        {
            "chief_complaint": "Formal nailfold capillaroscopy for systemic sclerosis staging.",
            "exam": (
                "NAILFOLD CAPILLAROSCOPY -- All 10 digits examined with video capillaroscope\n\n"
                "Pattern: ACTIVE SCLERODERMA PATTERN (Cutolo classification)\n\n"
                "Findings by finger:\n"
                "- Right D2: Giant capillaries (>50 mcm), 2 areas of dropout, microhemorrhage\n"
                "- Right D3: Multiple giant capillaries, bushy formations (neoangiogenesis)\n"
                "- Right D4: Capillary dropout (2 avascular areas), ramified capillaries\n"
                "- Right D5: Moderate dropout, microhemorrhages\n"
                "- Left D2: Giant capillaries, dropout, bushy neoangiogenesis\n"
                "- Left D3: Extensive dropout (3 avascular areas), ramified capillaries\n"
                "- Left D4: Giant capillaries, moderate dropout\n"
                "- Left D5: Dropout, microhemorrhages\n\n"
                "Mean capillary density: 4.2 per mm (reduced; normal >7 per mm)\n"
                "Giant capillary count: 12 per digit average (elevated; normal 0)"
            ),
            "assessment": (
                "ACTIVE SCLERODERMA PATTERN on nailfold capillaroscopy:\n"
                "- Extensive capillary dropout with avascular areas\n"
                "- Giant capillaries\n"
                "- Neoangiogenesis (bushy and ramified capillaries)\n"
                "- Microhemorrhages\n\n"
                "This pattern indicates ACTIVE microvascular disease and correlates with higher "
                "risk of digital ulceration and internal organ involvement (ILD, PAH). The "
                "progression from 'early' to 'active' pattern typically occurs within the first "
                "2-3 years of disease."
            ),
            "plan": (
                "1. Document for longitudinal monitoring\n"
                "2. Repeat capillaroscopy in 6-12 months\n"
                "3. Supportive care for digital vascular protection\n"
                "4. If digital ulcers develop, consider bosentan or iloprost"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 9: PFTs Follow-up -- September 2023
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-09-12_Pulmonology_PFTs_Follow_Up.pdf"),
        PATIENT, "09/12/2023", "pulmonology",
        "PULMONARY FUNCTION TESTS -- 6-MONTH FOLLOW-UP",
        {
            "chief_complaint": "Follow-up PFTs for SSc-ILD on mycophenolate + nintedanib x 5 months.",
            "hpi": (
                "Mrs. Foster returns for follow-up PFTs. She has been on mycophenolate 1500 mg BID "
                "and nintedanib 150 mg BID for 5 months. She reports stable exercise tolerance -- "
                "can climb 2 flights without significant dyspnea. Mild diarrhea from nintedanib "
                "managed with loperamide. No cough. No hemoptysis."
            ),
            "exam": (
                "SPIROMETRY:\n"
                "- FVC: 2.64 L (72% predicted) -- stable (prior 74%, within variability)\n"
                "- FEV1: 2.22 L (76% predicted)\n"
                "- FEV1/FVC: 0.84 (normal)\n\n"
                "DIFFUSION CAPACITY:\n"
                "- DLCO: 16.2 mL/min/mmHg (60% predicted) -- stable (prior 62%)\n\n"
                "6-MINUTE WALK TEST:\n"
                "- Distance: 412 meters (stable from 420 m)\n"
                "- SpO2 nadir: 93% (improved from 92%)\n"
                "- Borg dyspnea: 4/10 (improved from 5/10)"
            ),
            "assessment": (
                "SSc-ILD -- STABLE on dual therapy (mycophenolate + nintedanib)\n\n"
                "- FVC: 72% (stable from 74% -- within test variability, no clinically "
                "significant decline)\n"
                "- DLCO: 60% (stable from 62%)\n"
                "- 6MWT stable with mild improvement in desaturation and dyspnea\n\n"
                "This stability is a TREATMENT SUCCESS -- untreated SSc-ILD with anti-Scl-70 "
                "typically shows 5-10% FVC decline per year. The combination of mycophenolate "
                "and nintedanib appears to be halting progression."
            ),
            "plan": (
                "1. Continue current regimen: mycophenolate 1500 mg BID + nintedanib 150 mg BID\n"
                "2. Repeat PFTs in 6 months\n"
                "3. Annual HRCT (due December 2023)\n"
                "4. Annual echocardiogram for PAH surveillance\n"
                "5. Manage nintedanib-related diarrhea with loperamide PRN\n"
                "6. Follow-up in 6 months"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 10: Rheumatology 1-year follow-up -- January 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-01-22_Rheumatology_1Year.pdf"),
        PATIENT, "01/22/2024", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- 1 YEAR SINCE DIAGNOSIS",
        {
            "chief_complaint": "Annual SSc review. 1 year on immunosuppression.",
            "hpi": (
                "Mrs. Foster returns for 1-year review of dcSSc. She reports improved quality of "
                "life compared to diagnosis. Raynaud episodes reduced with nifedipine. Skin "
                "tightening has not progressed further and she feels mild softening over forearms. "
                "Exercise tolerance stable. Dysphagia improved with PPI. She developed one small "
                "digital pit on right index finger last month (healed). No renal crisis symptoms "
                "(BP well-controlled at home, average 122/78)."
            ),
            "exam": (
                "mRSS reassessment:\n"
                "- Fingers: 2/3 bilateral\n"
                "- Hands: 2/3 bilateral\n"
                "- Forearms: 1/3 bilateral (IMPROVED from 2/3)\n"
                "- Upper arms: 1/3 bilateral\n"
                "- Face: 1/3\n"
                "- Chest: 1/3\n"
                "mRSS Total: 16/51 (improved from 18 at peak -- stabilizing)\n\n"
                "Digital: Healed pit scar right D2 tip. No active ulcers.\n"
                "Chest: Bibasilar fine crackles (unchanged).\n"
                "Tendon friction rubs: Resolved at wrists."
            ),
            "assessment": (
                "dcSSc -- STABILIZED on mycophenolate + nintedanib\n\n"
                "1. Skin: mRSS 16 (improved from peak 18) -- skin progression halted, mild regression\n"
                "   - Resolution of tendon friction rubs is a positive sign\n"
                "2. ILD: Stable (FVC 72%, DLCO 60%)\n"
                "3. Pulmonary hypertension: Pending annual echo\n"
                "4. Digital vasculopathy: One healed pit -- monitoring\n"
                "5. GI: Esophageal dysmotility managed with PPI\n"
                "6. No renal crisis (BP normal, creatinine stable)\n"
                "7. CK normalizing (98 U/L) -- no clinical myopathy"
            ),
            "plan": (
                "1. Continue mycophenolate 1500 mg BID and nintedanib 150 mg BID\n"
                "2. Annual echocardiogram scheduled\n"
                "3. PFTs with DLCO in 3 months\n"
                "4. Labs: CBC, CMP, CRP, CK, urinalysis, NT-proBNP\n"
                "5. Continue nifedipine, omeprazole\n"
                "6. If digital ulcers recur, consider adding bosentan\n"
                "7. Influenza and COVID vaccination up to date\n"
                "8. Follow-up in 4 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2024-01-22_Lab_Annual.pdf"),
        PATIENT, "01/22/2024", PROVIDERS["rheumatology"],
        [{
            "panel_name": "SSc MONITORING",
            "results": [
                {"test": "CRP", "value": "4.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "22", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "CK", "value": "98", "unit": "U/L", "ref_range": "30-135", "flag": ""},
                {"test": "NT-proBNP", "value": "198", "unit": "pg/mL", "ref_range": "<125", "flag": "H"},
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "Hemoglobin", "value": "12.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "WBC", "value": "5.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "ALT", "value": "28", "unit": "U/L", "ref_range": "7-56", "flag": ""},
            ],
        }, {
            "panel_name": "URINALYSIS",
            "results": [
                {"test": "Protein", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Blood", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Specific Gravity", "value": "1.018", "unit": "", "ref_range": "1.005-1.030", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient Karen Foster: {count} documents generated")


if __name__ == "__main__":
    generate(os.path.join(os.path.dirname(__file__), "..", "demo_data", "karen_foster"))
