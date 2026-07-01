"""Patient: Michael Torres -- 41M, Graves' Disease.

Timeline: May 2024 -> May 2026
Specialists seen: PCP, Endocrinology, Cardiology, Ophthalmology, Nuclear Medicine
Key pattern: Weight loss + palpitations -> suppressed TSH / elevated fT4/fT3 ->
             TRAb positive -> thyroid uptake scan -> radioactive iodine ->
             hypothyroid management

The system should detect: Graves' disease, TRAb/TSI positive, thyroid eye disease,
                           atrial fibrillation.
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
    generate_lab_report,
    generate_progress_note,
)

PATIENT = {
    "name": "Michael Torres",
    "dob": "1984-02-15",
    "mrn": "MTO-2024-71843",
    "sex": "M",
    "age_at_start": 40,
}


def generate(output_dir: str):
    """Generate all clinical documents for Michael Torres."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: PCP -- May 2024 (weight loss, palpitations, tremor)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-05-06_PCP_Progress_Note.pdf"),
        PATIENT, "05/06/2024", "pcp",
        "OFFICE VISIT -- UNINTENTIONAL WEIGHT LOSS AND PALPITATIONS",
        {
            "chief_complaint": "Unintentional weight loss of 18 lbs over 3 months. Heart racing. Hand tremor. Heat intolerance.",
            "hpi": (
                "Mr. Torres is a 40-year-old male engineer presenting with a 3-month history of "
                "progressive unintentional weight loss (182 lbs down to 164 lbs) despite increased "
                "appetite. He reports daily palpitations described as his heart 'racing and pounding' "
                "at rest, lasting 30-60 minutes. He has developed a fine hand tremor that impairs "
                "his ability to write and use a mouse at work. He feels excessively hot and sweats "
                "easily -- his wife has noticed he is 'radiating heat' at night. He reports increased "
                "frequency of bowel movements (3-4 per day, soft, non-bloody). He feels anxious and "
                "irritable, with disrupted sleep (wakes at 3-4 AM unable to return to sleep). He "
                "noticed his eyes look 'wider' or 'more prominent' in photos from the past month. "
                "No neck pain or tenderness. No recent illness. Mother has hypothyroidism."
            ),
            "vitals": (
                "BP 148/72 (widened pulse pressure)  HR 112 (resting tachycardia, irregular)  "
                "Temp 99.2F  SpO2 99%  Wt 164 lbs  Ht 5'10\"  BMI 23.5"
            ),
            "ros": (
                "Constitutional: Weight loss, heat intolerance, excessive sweating, fatigue.\n"
                "HEENT: Eyes feel gritty and dry. Tearing. Eye prominence noticed by patient.\n"
                "Cardiovascular: Palpitations, racing heart. No chest pain. No syncope. No edema.\n"
                "Pulmonary: Mild dyspnea with exertion. No cough.\n"
                "GI: Increased stool frequency (3-4/day). No blood. Good appetite.\n"
                "Musculoskeletal: Proximal muscle weakness. Difficulty climbing stairs.\n"
                "Neurological: Bilateral hand tremor. Anxiety. Insomnia.\n"
                "Skin: Warm, moist. No rash."
            ),
            "exam": (
                "General: Thin, anxious-appearing man. Diaphoretic. Fidgety.\n"
                "HEENT: Eyes -- bilateral proptosis (mild exophthalmos). Lid retraction with "
                "visible sclera above iris (Dalrymple sign). Mild conjunctival injection. Lid lag "
                "on downgaze (von Graefe sign). EOMs intact, no diplopia. No periorbital edema.\n"
                "Neck: Thyroid diffusely enlarged (approximately 2x normal), firm, non-tender, "
                "no nodules. Thyroid bruit audible bilaterally. No lymphadenopathy.\n"
                "Cardiovascular: Tachycardic, IRREGULARLY IRREGULAR rhythm consistent with "
                "atrial fibrillation. Hyperdynamic precordium. Systolic flow murmur II/VI. "
                "No JVD. No peripheral edema.\n"
                "Pulmonary: Clear bilaterally.\n"
                "Extremities: Warm, moist palms. Fine bilateral hand tremor on extension. "
                "Hyperreflexia (3+) throughout.\n"
                "Skin: Warm, moist, velvety. No pretibial myxedema."
            ),
            "assessment": (
                "1. THYROTOXICOSIS -- classic presentation: weight loss despite increased appetite, "
                "tachycardia, tremor, heat intolerance, hyperdefecation, widened pulse pressure, "
                "diffuse goiter with bruit, and lid signs\n"
                "2. ATRIAL FIBRILLATION -- likely thyrotoxicosis-induced. New onset. CHA2DS2-VASc "
                "score 0, but needs rate control\n"
                "3. Probable GRAVES' DISEASE given: diffuse goiter, thyroid bruit, and clinical "
                "ophthalmopathy (proptosis, lid retraction, lid lag)\n"
                "4. Thyroid eye disease (Graves ophthalmopathy) -- mild, needs ophthalmology eval"
            ),
            "plan": (
                "1. STAT thyroid labs: TSH, free T4, free T3, total T3\n"
                "2. Thyroid antibodies: TRAb (TSH receptor antibody), anti-TPO, thyroglobulin\n"
                "3. CBC, CMP, LFTs (baseline before anti-thyroid meds)\n"
                "4. ECG -- confirm atrial fibrillation\n"
                "5. Start metoprolol tartrate 25 mg BID for rate control and symptom relief\n"
                "6. URGENT endocrinology referral\n"
                "7. Ophthalmology referral for thyroid eye disease assessment\n"
                "8. Cardiology referral for new atrial fibrillation\n"
                "9. CHA2DS2-VASc = 0 -- no anticoagulation needed at this time\n"
                "10. Patient counseled: no heavy exercise until heart rate controlled"
            ),
        },
    )

    # Lab: May 2024 -- initial thyroid labs
    generate_lab_report(
        os.path.join(output_dir, "2024-05-06_Lab_Thyroid_Initial.pdf"),
        PATIENT, "05/06/2024", PROVIDERS["pcp"],
        [{
            "panel_name": "THYROID FUNCTION",
            "results": [
                {"test": "TSH", "value": "<0.01", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": "LL"},
                {"test": "Free T4", "value": "4.8", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": "HH"},
                {"test": "Free T3", "value": "18.2", "unit": "pg/mL", "ref_range": "2.3-4.2", "flag": "HH"},
                {"test": "Total T3", "value": "412", "unit": "ng/dL", "ref_range": "80-200", "flag": "HH"},
            ],
        }, {
            "panel_name": "THYROID ANTIBODIES",
            "results": [
                {"test": "TRAb (TSH Receptor Antibody)", "value": "12.8", "unit": "IU/L", "ref_range": "<1.75", "flag": "HH"},
                {"test": "TSI (Thyroid Stimulating Immunoglobulin)", "value": "485", "unit": "%", "ref_range": "<140", "flag": "HH"},
                {"test": "Anti-TPO Antibody", "value": "218", "unit": "IU/mL", "ref_range": "<35", "flag": "HH"},
                {"test": "Thyroglobulin", "value": "68", "unit": "ng/mL", "ref_range": "3-40", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC AND LFTs",
            "results": [
                {"test": "WBC", "value": "5.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "14.2", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "Platelets", "value": "262", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "ALT", "value": "48", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "42", "unit": "U/L", "ref_range": "10-40", "flag": "H"},
                {"test": "ALP", "value": "128", "unit": "U/L", "ref_range": "44-147", "flag": ""},
                {"test": "Total Bilirubin", "value": "0.8", "unit": "mg/dL", "ref_range": "0.1-1.2", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: Endocrinology Consult -- May 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-05-13_Endocrinology_Consult.pdf"),
        PATIENT, "05/13/2024", "endocrinology",
        "ENDOCRINOLOGY CONSULTATION -- GRAVES' DISEASE",
        {
            "chief_complaint": (
                "Referred for management of thyrotoxicosis. TSH <0.01, FT4 4.8, FT3 18.2. "
                "TRAb 12.8, TSI 485%. New atrial fibrillation."
            ),
            "hpi": (
                "Mr. Torres is a 40-year-old man with severe biochemical and clinical "
                "thyrotoxicosis. Labs confirm Graves' disease: TSH suppressed (<0.01), "
                "markedly elevated FT4 (4.8) and FT3 (18.2), and strongly positive TRAb (12.8) "
                "and TSI (485%). He has a diffuse goiter with bruit, mild clinical ophthalmopathy "
                "(proptosis, lid retraction), and new atrial fibrillation. He was started on "
                "metoprolol 25 mg BID by PCP. Heart rate currently 88, regular (converted to "
                "sinus rhythm on beta-blocker). Symptoms mildly improved on beta-blocker but "
                "still has tremor, heat intolerance, and 2-lb weight loss since last visit."
            ),
            "exam": (
                "Thyroid: Diffuse goiter, approximately 40 grams (2x normal). Firm, smooth, "
                "non-tender. Bilateral bruit. No nodules palpable.\n"
                "Eyes: Mild bilateral proptosis. Hertel exophthalmometry: OD 21 mm, OS 20 mm "
                "(upper normal 20 mm). Lid retraction bilateral. No periorbital edema. "
                "EOMs full, no diplopia or pain.\n"
                "Cardiovascular: Now in sinus rhythm (SR). Rate 88. No murmur.\n"
                "Tremor: Fine bilateral, persistent."
            ),
            "assessment": (
                "GRAVES' DISEASE -- confirmed\n\n"
                "Diagnostic certainty: definitive\n"
                "- Suppressed TSH with elevated FT4/FT3\n"
                "- Strongly positive TRAb (12.8) and TSI (485%)\n"
                "- Diffuse goiter with bruit\n"
                "- Clinical ophthalmopathy\n\n"
                "Severity: SEVERE thyrotoxicosis\n"
                "- FT4 nearly 3x upper limit\n"
                "- FT3 >4x upper limit (T3-predominant production)\n"
                "- Atrial fibrillation (converted with beta-blocker)\n\n"
                "Graves Ophthalmopathy: Mild (CAS 2/7)\n"
                "- Lid retraction, proptosis\n"
                "- No inflammatory signs (edema, chemosis)\n"
                "- EUGOGO classification: mild, inactive\n\n"
                "TREATMENT OPTIONS discussed:\n"
                "1. Anti-thyroid drugs (methimazole) -- first-line, trial x 12-18 months\n"
                "2. Radioactive iodine (RAI) -- definitive, but contraindicated with active TED\n"
                "3. Thyroidectomy -- surgical option\n\n"
                "Patient prefers trial of methimazole first. RAI discussed as definitive option "
                "if relapse occurs."
            ),
            "plan": (
                "1. START methimazole 20 mg daily (high dose given severity)\n"
                "2. Continue metoprolol 25 mg BID\n"
                "3. Baseline ANC (agranulocytosis risk with methimazole)\n"
                "4. LFTs monitoring q4-6 weeks\n"
                "5. CBC q4 weeks x 3 months (agranulocytosis screening)\n"
                "6. Thyroid function (TSH, FT4) in 4 weeks -- goal: normalize within 4-8 weeks\n"
                "7. Thyroid uptake and scan -- ordered to document baseline and confirm diagnosis\n"
                "8. Ophthalmology evaluation for TED staging\n"
                "9. Smoking cessation counseling (smoking worsens TED) -- patient smokes 3-4 "
                "cigarettes/day socially\n"
                "10. Patient education: sore throat with fever = STOP methimazole and get "
                "STAT CBC (agranulocytosis warning)\n"
                "11. Follow-up in 4 weeks"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: Thyroid Uptake and Scan -- May 2024
    # ─────────────────────────────────────────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2024-05-20_Thyroid_Uptake_Scan.pdf"),
        PATIENT, "05/20/2024", "nuclear_medicine",
        "Nuclear Medicine", "Thyroid Uptake and Scan (I-123)",
        "40-year-old male with Graves' disease. Baseline uptake and scan prior to treatment.",
        "I-123 (0.2 mCi) administered orally. Uptake measurements at 4 hours and 24 hours. "
        "Pinhole gamma camera images of thyroid in anterior and oblique projections.",
        (
            "UPTAKE:\n"
            "- 4-hour uptake: 38% (normal 5-15%) -- ELEVATED\n"
            "- 24-hour uptake: 62% (normal 15-35%) -- MARKEDLY ELEVATED\n\n"
            "SCAN:\n"
            "- Thyroid gland is diffusely enlarged.\n"
            "- Homogeneous, intensely increased radiotracer uptake throughout the entire gland.\n"
            "- No focal cold or hot nodules.\n"
            "- Both lobes symmetrically enlarged.\n"
            "- Estimated thyroid weight: 40 grams (normal 15-20 g).\n"
            "- Isthmus shows proportional uptake.\n"
            "- No substernal extension.\n"
            "- Pyramidal lobe visualized (consistent with thyroid stimulation)."
        ),
        (
            "1. DIFFUSELY ELEVATED RADIOIODINE UPTAKE consistent with GRAVES' DISEASE.\n"
            "   - 24-hour uptake 62% (nearly 2x upper normal limit)\n"
            "   - Homogeneous, intense uptake without focal lesions\n"
            "2. Diffuse thyromegaly (estimated 40 grams).\n"
            "3. No cold nodules to suggest malignancy.\n"
            "4. If RAI therapy planned: estimated therapeutic dose would be approximately "
            "10-15 mCi I-131 based on gland size and uptake."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: Ophthalmology -- June 2024 (TED evaluation)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-06-03_Ophth_TED_Evaluation.pdf"),
        PATIENT, "06/03/2024", "ophthalmology",
        "OPHTHALMOLOGY CONSULTATION -- THYROID EYE DISEASE EVALUATION",
        {
            "chief_complaint": "Referred for evaluation of thyroid eye disease in Graves' disease.",
            "hpi": (
                "Mr. Torres is a 40-year-old male with newly diagnosed Graves' disease (TRAb 12.8, "
                "TSI 485%) referred for TED assessment. He reports bilateral dry, gritty eyes, "
                "excessive tearing, and his wife notes his eyes appear more prominent. No diplopia. "
                "No color vision changes. No pain with eye movement. He is a light smoker "
                "(3-4 cigarettes/day). Currently on methimazole and metoprolol."
            ),
            "exam": (
                "Visual Acuity: OD 20/20. OS 20/20.\n"
                "Color Vision: Normal OU.\n"
                "Pupil: No RAPD. Equal and reactive.\n"
                "Hertel Exophthalmometry (base 100 mm): OD 22 mm, OS 21 mm. (Upper normal 20 mm "
                "for males.)\n"
                "Lid Measurements:\n"
                "- Upper lid retraction: OD 3 mm, OS 2 mm (normal <2 mm)\n"
                "- Lagophthalmos: 1 mm bilateral on gentle closure\n"
                "- Lid lag present on downgaze bilateral\n"
                "Slit Lamp: Mild conjunctival injection. No chemosis. Corneas clear, no SPK. "
                "Tear film adequate but unstable (TBUT 6 sec OD, 7 sec OS).\n"
                "IOP: OD 18, OS 17 (primary). OD 22, OS 21 (upgaze -- elevation on upgaze suggests "
                "mild inferior rectus enlargement).\n"
                "Fundoscopy: Discs sharp, pink, well-perfused OU.\n"
                "EOMs: Full OU. No restriction. No diplopia in any gaze.\n\n"
                "Clinical Activity Score (CAS): 2/7\n"
                "- Spontaneous retrobulbar pain: No (0)\n"
                "- Pain on eye movement: No (0)\n"
                "- Eyelid erythema: No (0)\n"
                "- Conjunctival injection: Yes (1)\n"
                "- Chemosis: No (0)\n"
                "- Caruncle swelling: No (0)\n"
                "- Eyelid edema: Yes, minimal (1)"
            ),
            "assessment": (
                "THYROID EYE DISEASE -- MILD, INACTIVE\n\n"
                "EUGOGO Severity: MILD\n"
                "- Lid retraction <3 mm\n"
                "- Proptosis <3 mm above normal\n"
                "- No diplopia\n"
                "- No optic neuropathy\n\n"
                "CAS 2/7 -- INACTIVE disease (CAS <3 = inactive)\n\n"
                "Risk factors for progression:\n"
                "- Active smoking (strongest modifiable risk factor)\n"
                "- High TRAb levels\n"
                "- Male sex (males have worse TED outcomes)"
            ),
            "plan": (
                "1. Lubricating eye drops (artificial tears) QID\n"
                "2. Nighttime ointment for lagophthalmos\n"
                "3. CRITICAL: smoking cessation -- strongest modifiable risk factor for TED "
                "progression. Referred to smoking cessation program.\n"
                "4. Selenium supplementation 200 mcg daily (EUGOGO trial -- reduces TED "
                "progression and improves QOL in mild TED)\n"
                "5. Sleep with head elevated 30 degrees (reduces periorbital edema)\n"
                "6. UV protection with wrap-around sunglasses\n"
                "7. Follow-up in 3 months or sooner if symptoms worsen\n"
                "8. If RAI therapy chosen for Graves: MUST have steroid prophylaxis to prevent "
                "TED worsening (prednisone 0.3-0.5 mg/kg x 3 months)"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: Cardiology Consult -- June 2024 (atrial fibrillation)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-06-10_Cardiology_Consult.pdf"),
        PATIENT, "06/10/2024", "cardiology",
        "CARDIOLOGY CONSULTATION -- THYROTOXICOSIS-INDUCED ATRIAL FIBRILLATION",
        {
            "chief_complaint": "Referred for new-onset atrial fibrillation in setting of Graves' disease.",
            "hpi": (
                "Mr. Torres is a 40-year-old male with Graves' disease diagnosed 1 month ago, "
                "presenting with new-onset atrial fibrillation detected at PCP visit. He converted "
                "to sinus rhythm on metoprolol 25 mg BID. Currently on methimazole 20 mg daily. "
                "He still has intermittent palpitations lasting 5-10 minutes, 2-3x per week, "
                "though less severe than before beta-blocker. No syncope. No chest pain. "
                "No lower extremity edema."
            ),
            "exam": (
                "Cardiovascular: Sinus rhythm, rate 76. Regular. No murmurs. No JVD. "
                "No peripheral edema. S1/S2 normal. No S3 or S4.\n"
                "ECG: Normal sinus rhythm, rate 78. Normal axis. PR interval 160 ms. "
                "QRS 88 ms. No ST changes. No LVH criteria."
            ),
            "assessment": (
                "1. Paroxysmal atrial fibrillation -- thyrotoxicosis-induced\n"
                "   - AF occurs in 10-15% of hyperthyroid patients\n"
                "   - Expected to resolve with achievement of euthyroid state\n"
                "   - Currently in sinus rhythm on metoprolol\n"
                "2. CHA2DS2-VASc: 0 -- no indication for anticoagulation\n"
                "3. No structural heart disease on echocardiogram (ordered today)"
            ),
            "plan": (
                "1. Echocardiogram -- evaluate LV function and chamber sizes\n"
                "2. Continue metoprolol 25 mg BID -- increase to 50 mg BID if palpitations persist\n"
                "3. Holter monitor x 48 hours -- document AF burden\n"
                "4. Primary management: achieve euthyroid state (endocrinology managing)\n"
                "5. If AF persists despite euthyroid state x 3 months, consider cardioversion\n"
                "6. No anticoagulation needed (CHA2DS2-VASc 0)\n"
                "7. Follow-up in 3 months or sooner if symptoms worsen"
            ),
        },
    )

    # Echocardiogram
    generate_imaging_report(
        os.path.join(output_dir, "2024-06-10_Echo_AF_Evaluation.pdf"),
        PATIENT, "06/10/2024", "cardiology",
        "Echocardiogram", "Transthoracic",
        "40-year-old male with Graves' disease and paroxysmal atrial fibrillation. "
        "Evaluate cardiac structure and function.",
        "Complete 2D, M-mode, and Doppler transthoracic echocardiogram.",
        (
            "Left ventricle:\n"
            "- Normal size. LVEDD 5.0 cm. Normal wall thickness.\n"
            "- Ejection fraction: 58% (low-normal)\n"
            "- No regional wall motion abnormalities.\n"
            "- Hyperdynamic LV function pattern.\n\n"
            "Left atrium: Mildly dilated (4.2 cm, normal <4.0 cm).\n\n"
            "Right ventricle: Normal size and function. TAPSE 2.2 cm.\n\n"
            "Valves:\n"
            "- Mitral: Trace regurgitation. No stenosis.\n"
            "- Tricuspid: Trace regurgitation. RVSP 28 mmHg (normal).\n"
            "- Aortic: No stenosis or regurgitation.\n\n"
            "Pericardium: No effusion.\n"
            "IVC: Normal caliber with normal respiratory variation."
        ),
        (
            "1. MILDLY DILATED LEFT ATRIUM (4.2 cm) -- consistent with thyrotoxic state "
            "and paroxysmal AF.\n"
            "2. LOW-NORMAL LV ejection fraction (58%) -- hyperthyroid cardiomyopathy "
            "should be considered. May improve with euthyroid state.\n"
            "3. No structural valve disease.\n"
            "4. No pulmonary hypertension.\n\n"
            "RECOMMENDATION: Repeat echo after achieving euthyroid state to assess for "
            "reversibility of LA dilation and LV function."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Endocrinology 4-week follow-up -- June 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-06-17_Endocrinology_Follow_Up.pdf"),
        PATIENT, "06/17/2024", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- 4 WEEKS ON METHIMAZOLE",
        {
            "chief_complaint": "Follow-up Graves' disease. 4 weeks on methimazole 20 mg.",
            "hpi": (
                "Mr. Torres returns 4 weeks after starting methimazole. He reports moderate "
                "improvement: tremor less prominent, heat intolerance improving, weight stable "
                "(164 lbs -- no further loss). Palpitations reduced to 1-2x/week. Still has "
                "mild anxiety and sleep disruption. Eye symptoms stable. He has reduced smoking "
                "to 1-2 cigarettes/day (from 3-4). No sore throat or fever. No rash."
            ),
            "assessment": (
                "Graves' Disease -- IMPROVING on methimazole\n\n"
                "Thyroid function trending toward normal:\n"
                "- TSH: still suppressed (<0.01) -- expected, lag of 4-8 weeks\n"
                "- FT4: 2.1 (down from 4.8 -- 56% reduction, excellent response)\n"
                "- FT3: 8.4 (down from 18.2 -- 54% reduction)\n\n"
                "CBC: ANC 2.8 K/uL (normal -- no agranulocytosis)\n"
                "LFTs: ALT 32 (improving from 48 -- no hepatotoxicity)\n\n"
                "Clinical improvement correlates with biochemical improvement."
            ),
            "plan": (
                "1. Reduce methimazole to 10 mg daily (block-replace strategy not needed)\n"
                "2. Continue metoprolol 25 mg BID\n"
                "3. Recheck TSH, FT4, FT3, CBC in 4 weeks\n"
                "4. Goal: FT4 in mid-normal range\n"
                "5. Continue smoking reduction -- goal: cessation\n"
                "6. If euthyroid for 12-18 months, may trial methimazole discontinuation\n"
                "7. Discussed: 50% relapse rate after methimazole -- RAI remains definitive option\n"
                "8. Follow-up in 4 weeks"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2024-06-17_Lab_Thyroid_Follow_Up.pdf"),
        PATIENT, "06/17/2024", PROVIDERS["endocrinology"],
        [{
            "panel_name": "THYROID FUNCTION -- 4 WEEKS ON METHIMAZOLE",
            "results": [
                {"test": "TSH", "value": "<0.01", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": "LL"},
                {"test": "Free T4", "value": "2.1", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": "H"},
                {"test": "Free T3", "value": "8.4", "unit": "pg/mL", "ref_range": "2.3-4.2", "flag": "H"},
            ],
        }, {
            "panel_name": "SAFETY MONITORING",
            "results": [
                {"test": "WBC", "value": "6.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "ANC", "value": "2.8", "unit": "K/uL", "ref_range": ">1.5", "flag": ""},
                {"test": "ALT", "value": "32", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "28", "unit": "U/L", "ref_range": "10-40", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 7: Endocrinology -- November 2024 (relapse, RAI decision)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-11-04_Endocrinology_Relapse.pdf"),
        PATIENT, "11/04/2024", "endocrinology",
        "ENDOCRINOLOGY VISIT -- GRAVES' DISEASE RELAPSE ON METHIMAZOLE",
        {
            "chief_complaint": "Recurrent hyperthyroid symptoms. Palpitations, tremor, and 8-lb weight loss over 6 weeks.",
            "hpi": (
                "Mr. Torres was doing well on methimazole 10 mg daily, achieved near-euthyroid "
                "state by August 2024 (FT4 1.4, FT3 3.8, TSH 0.08). Methimazole was reduced to "
                "5 mg daily in September as a trial taper. Over the past 6 weeks, he has had "
                "recurrence of palpitations, hand tremor, heat intolerance, and weight loss (172 "
                "to 164 lbs). He admits to intermittent medication non-adherence (missed doses "
                "3-4x/week). TRAb remains strongly positive (10.2). This represents treatment "
                "failure / relapse and indicates need for definitive therapy."
            ),
            "assessment": (
                "Graves' Disease -- RELAPSE during methimazole taper\n\n"
                "Current labs:\n"
                "- TSH: <0.01 (re-suppressed)\n"
                "- FT4: 3.2 (re-elevated)\n"
                "- FT3: 12.6 (re-elevated)\n"
                "- TRAb: 10.2 (persistently elevated -- poor prognostic factor for remission)\n\n"
                "Relapse predictors present (3 of 4 = 90% relapse risk):\n"
                "1. Male sex\n"
                "2. Large goiter (40 g)\n"
                "3. Persistently elevated TRAb (>3.5)\n"
                "4. T3-predominant thyrotoxicosis\n\n"
                "DEFINITIVE THERAPY RECOMMENDED: Radioactive iodine (RAI)\n"
                "- Ophthalmology confirms TED is mild and INACTIVE (CAS 2/7)\n"
                "- Can proceed with RAI with steroid prophylaxis\n"
                "- Patient agrees to RAI after discussing options"
            ),
            "plan": (
                "1. Restart methimazole 20 mg daily -- bridge to RAI\n"
                "2. Schedule RAI (I-131) in 4-6 weeks\n"
                "3. Stop methimazole 3-5 days before RAI\n"
                "4. Steroid prophylaxis for TED: prednisone 0.3 mg/kg/day starting day of RAI, "
                "taper over 3 months (per EUGOGO guidelines for mild TED)\n"
                "5. Pre-RAI: pregnancy test (N/A male), thyroid uptake (already documented 62%)\n"
                "6. Post-RAI: expect hypothyroidism within 2-6 months -- will start levothyroxine\n"
                "7. Increase metoprolol to 50 mg BID for symptom control\n"
                "8. Follow-up in 2 weeks, then RAI day"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2024-11-04_Lab_Relapse.pdf"),
        PATIENT, "11/04/2024", PROVIDERS["endocrinology"],
        [{
            "panel_name": "THYROID FUNCTION -- RELAPSE",
            "results": [
                {"test": "TSH", "value": "<0.01", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": "LL"},
                {"test": "Free T4", "value": "3.2", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": "HH"},
                {"test": "Free T3", "value": "12.6", "unit": "pg/mL", "ref_range": "2.3-4.2", "flag": "HH"},
                {"test": "TRAb", "value": "10.2", "unit": "IU/L", "ref_range": "<1.75", "flag": "HH"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 8: RAI Treatment -- December 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-12-09_Nuclear_Medicine_RAI.pdf"),
        PATIENT, "12/09/2024", "nuclear_medicine",
        "RADIOACTIVE IODINE THERAPY -- GRAVES' DISEASE",
        {
            "chief_complaint": "Scheduled I-131 ablation for Graves' disease refractory to methimazole.",
            "hpi": (
                "Mr. Torres presents for radioactive iodine (I-131) therapy. Methimazole was "
                "discontinued 5 days ago. He has mild recurrent symptoms (tremor, palpitations) "
                "managed with metoprolol 50 mg BID. Pre-treatment labs confirm persistent "
                "thyrotoxicosis. Ophthalmology cleared for RAI with steroid prophylaxis. "
                "Prednisone 0.3 mg/kg/day (15 mg) started today."
            ),
            "assessment": (
                "Graves' disease, refractory to methimazole. Proceeding with definitive RAI therapy.\n\n"
                "Dose calculation:\n"
                "- Estimated thyroid weight: 40 grams\n"
                "- 24-hour RAIU: 62%\n"
                "- Target dose: 150 mcCi/gram\n"
                "- Calculated activity: (40 x 150) / 0.62 = 9.7 mCi, rounded to 12 mCi\n"
                "  (higher dose chosen given large gland and high TRAb)\n\n"
                "I-131 dose administered: 12 mCi (444 MBq) orally"
            ),
            "plan": (
                "1. Radiation safety precautions reviewed with patient:\n"
                "   - Avoid close contact with children and pregnant women x 5 days\n"
                "   - Sleep alone x 3 days\n"
                "   - Use separate utensils, flush twice after bathroom use\n"
                "   - Return to work after 3 days if desk job\n"
                "2. Prednisone 15 mg daily x 4 weeks, then taper by 5 mg q2w (total 3 months)\n"
                "3. Continue metoprolol 50 mg BID\n"
                "4. DO NOT restart methimazole\n"
                "5. Thyroid labs (TSH, FT4) in 4 weeks, then q4-6 weeks\n"
                "6. Expect hypothyroidism within 2-6 months -- start levothyroxine when TSH rises\n"
                "7. Endocrinology follow-up in 4 weeks"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 9: Endocrinology Follow-up -- February 2025 (now hypothyroid)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2025-02-10_Endocrinology_Post_RAI.pdf"),
        PATIENT, "02/10/2025", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- POST-RAI HYPOTHYROIDISM",
        {
            "chief_complaint": "Follow-up 2 months post-RAI. New symptoms of fatigue, cold intolerance, constipation.",
            "hpi": (
                "Mr. Torres returns 2 months after I-131 therapy. Initially he felt well as "
                "hyperthyroid symptoms resolved. Over the past 3 weeks, he has developed new "
                "symptoms: fatigue, cold intolerance (reversal from heat intolerance), constipation, "
                "dry skin, and weight gain of 6 lbs (164 to 170 lbs). Palpitations have completely "
                "resolved. Tremor resolved. Eye symptoms stable -- no worsening on prednisone taper "
                "(currently 5 mg daily, last 2 weeks of taper)."
            ),
            "assessment": (
                "POST-RAI HYPOTHYROIDISM -- EXPECTED outcome\n\n"
                "Labs:\n"
                "- TSH: 48.2 (elevated -- hypothyroid)\n"
                "- FT4: 0.4 (below normal)\n"
                "- FT3: 1.8 (below normal)\n"
                "- TRAb: 4.2 (declining from 10.2 -- good response to RAI)\n\n"
                "RAI was successful in ablating thyroid function. Patient is now hypothyroid "
                "as expected. Needs lifelong levothyroxine replacement.\n\n"
                "Graves Ophthalmopathy: STABLE on steroid taper. No worsening post-RAI."
            ),
            "plan": (
                "1. START levothyroxine 88 mcg daily (1.6 mcg/kg based on weight 170 lbs)\n"
                "   - Take on empty stomach, 30-60 min before breakfast\n"
                "   - Separate from calcium, iron, PPIs by 4 hours\n"
                "2. Complete prednisone taper (5 mg x 2 more weeks, then stop)\n"
                "3. Discontinue metoprolol (no more tachycardia)\n"
                "4. Recheck TSH, FT4 in 6 weeks -- adjust levothyroxine to TSH goal 0.5-2.0\n"
                "5. Annual thyroid labs once stable\n"
                "6. Ophthalmology follow-up in 3 months (post-steroid taper)\n"
                "7. Follow-up in 6 weeks"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2025-02-10_Lab_Post_RAI.pdf"),
        PATIENT, "02/10/2025", PROVIDERS["endocrinology"],
        [{
            "panel_name": "THYROID FUNCTION -- POST-RAI",
            "results": [
                {"test": "TSH", "value": "48.2", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": "HH"},
                {"test": "Free T4", "value": "0.4", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": "L"},
                {"test": "Free T3", "value": "1.8", "unit": "pg/mL", "ref_range": "2.3-4.2", "flag": "L"},
                {"test": "TRAb", "value": "4.2", "unit": "IU/L", "ref_range": "<1.75", "flag": "H"},
                {"test": "Thyroglobulin", "value": "12", "unit": "ng/mL", "ref_range": "3-40", "flag": ""},
            ],
        }, {
            "panel_name": "METABOLIC",
            "results": [
                {"test": "Total Cholesterol", "value": "248", "unit": "mg/dL", "ref_range": "<200", "flag": "H"},
                {"test": "LDL", "value": "158", "unit": "mg/dL", "ref_range": "<130", "flag": "H"},
                {"test": "CK", "value": "210", "unit": "U/L", "ref_range": "30-200", "flag": "H"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 10: Endocrinology Stable -- May 2025
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2025-05-12_Endocrinology_Stable.pdf"),
        PATIENT, "05/12/2025", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- EUTHYROID ON LEVOTHYROXINE",
        {
            "chief_complaint": "Follow-up post-RAI Graves' disease. Now on levothyroxine replacement.",
            "hpi": (
                "Mr. Torres returns for follow-up. He has been on levothyroxine 100 mcg daily "
                "(increased from 88 mcg at last visit). He reports feeling well -- energy normalized, "
                "weight stable at 172 lbs, no cold intolerance, no palpitations. Bowel habits "
                "normal. Eye symptoms stable -- discontinued prednisone in March without TED flare. "
                "He quit smoking completely 2 months ago. AF has not recurred (confirmed by recent "
                "Holter monitor showing sinus rhythm throughout)."
            ),
            "assessment": (
                "POST-RAI GRAVES' DISEASE -- EUTHYROID on levothyroxine\n\n"
                "Labs:\n"
                "- TSH: 1.8 (TARGET range 0.5-2.0 -- OPTIMAL)\n"
                "- FT4: 1.2 (normal)\n"
                "- TRAb: 2.1 (declining, approaching normal -- RAI effectiveness)\n\n"
                "Outcomes:\n"
                "1. Thyroid: Euthyroid on levothyroxine 100 mcg -- stable\n"
                "2. Cardiac: No recurrence of AF. LVEF likely normalized (repeat echo recommended)\n"
                "3. TED: Stable, no worsening post-RAI. Steroid prophylaxis successful.\n"
                "4. Smoking: Quit -- major win for TED prognosis\n"
                "5. Lipids normalizing with euthyroid state (cholesterol 202)"
            ),
            "plan": (
                "1. Continue levothyroxine 100 mcg daily\n"
                "2. Annual TSH, FT4\n"
                "3. Ophthalmology annual follow-up for TED monitoring\n"
                "4. Cardiology: repeat echocardiogram to document LV function recovery\n"
                "5. Lipid panel in 3 months on euthyroid state\n"
                "6. TRAb annual monitoring -- if normalizes, very low risk of TED worsening\n"
                "7. Lifelong levothyroxine replacement\n"
                "8. Follow-up in 6 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2025-05-12_Lab_Euthyroid.pdf"),
        PATIENT, "05/12/2025", PROVIDERS["endocrinology"],
        [{
            "panel_name": "THYROID FUNCTION -- ON LEVOTHYROXINE",
            "results": [
                {"test": "TSH", "value": "1.8", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Free T4", "value": "1.2", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": ""},
                {"test": "TRAb", "value": "2.1", "unit": "IU/L", "ref_range": "<1.75", "flag": "H"},
                {"test": "Anti-TPO", "value": "142", "unit": "IU/mL", "ref_range": "<35", "flag": "H"},
                {"test": "Thyroglobulin", "value": "3.2", "unit": "ng/mL", "ref_range": "3-40", "flag": ""},
            ],
        }, {
            "panel_name": "METABOLIC",
            "results": [
                {"test": "Total Cholesterol", "value": "202", "unit": "mg/dL", "ref_range": "<200", "flag": "H"},
                {"test": "LDL", "value": "128", "unit": "mg/dL", "ref_range": "<130", "flag": ""},
                {"test": "CK", "value": "118", "unit": "U/L", "ref_range": "30-200", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient Michael Torres: {count} documents generated")


if __name__ == "__main__":
    generate(os.path.join(os.path.dirname(__file__), "..", "demo_data", "michael_torres"))
