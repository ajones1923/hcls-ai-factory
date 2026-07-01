"""Dismissal documents -- where real symptoms are attributed to benign causes.

These are the highest-impact demo documents. When the AI system ingests the
full patient record, it flags these visits as missed diagnostic opportunities:
"Patient was told this was anxiety, but orthostatic HR data shows POTS."

Documents:
- Maya Rodriguez: Psychiatry evaluation diagnosing anxiety/panic disorder (actually POTS)
- Maya Rodriguez: Dermatology consult diagnosing rosacea (actually MCAS flushing)
- Sarah Mitchell: ER visit diagnosing "viral syndrome" (actually lupus flare)
- Sarah Mitchell: Urgent care visit for pleurisy dismissed as costochondritis
- Linda Chen: PCP visit attributing sicca to perimenopause (actually Sjogren's)
- David Park: Sleep medicine consult for insomnia (actually AS night pain)
- Rachel Thompson: PCP dismissing ANA as "could be nothing" (already in existing
  notes as text, but adding explicit dismissal documentation)

Author: Adam Jones
Date: March 2026
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_lab_report,
    generate_progress_note,
    generate_referral_letter,
)

# Patient references (same as other files)
MAYA = {
    "name": "Maya Rodriguez",
    "dob": "1996-08-14",
    "mrn": "MRD-2021-45678",
    "sex": "F",
    "age_at_start": 24,
}

SARAH = {
    "name": "Sarah Mitchell",
    "dob": "1991-03-22",
    "mrn": "SMI-2022-12345",
    "sex": "F",
    "age_at_start": 31,
}

LINDA = {
    "name": "Linda Chen",
    "dob": "1973-05-19",
    "mrn": "LCH-2022-67890",
    "sex": "F",
    "age_at_start": 47,
}

DAVID = {
    "name": "David Park",
    "dob": "1980-11-03",
    "mrn": "DPA-2019-33102",
    "sex": "M",
    "age_at_start": 38,
}


def generate_maya_dismissals(output_dir: str):
    """Generate dismissal documents for Maya Rodriguez."""
    os.makedirs(output_dir, exist_ok=True)

    # ──────────────────────────────────────────────────────────────
    # PCP Referral to Psychiatry (2021-08-30)
    # Context: After first ER syncope visit (July 2021), PCP refers to psych.
    # ──────────────────────────────────────────────────────────────
    generate_referral_letter(
        os.path.join(output_dir, "2021-08-30_Referral_PCP_to_Psychiatry.pdf"),
        MAYA, "08/30/2021",
        "pcp", "psychiatry",
        PROVIDERS["psychiatry"],
        (
            "Evaluation for anxiety disorder. This 25-year-old woman has recurrent episodes "
            "of tachycardia, dizziness, lightheadedness, and near-syncope. She was seen in the "
            "Emergency Department in July 2021 for syncope -- workup was unrevealing. Episodes "
            "are accompanied by significant anxiety and fear of passing out."
        ),
        (
            "Ms. Rodriguez is a 25-year-old graduate student in generally good health. Over "
            "the past 5 months she has developed episodes of rapid heartbeat, dizziness, and "
            "lightheadedness occurring 3-4 times per week, often when standing from seated "
            "position or during prolonged standing. She had a syncopal episode in July and was "
            "evaluated in the ED where basic labs, EKG, and CT head were normal.\n\n"
            "She reports significant distress about the episodes and has begun avoiding certain "
            "activities (standing in lines, crowded places). She has reduced her exercise due to "
            "fear of passing out. She also endorses intermittent nausea, abdominal discomfort, "
            "and brain fog.\n\n"
            "Relevant vitals from today: BP 108/62 seated, HR 72 seated. Standing: BP 100/68, "
            "HR 108 (increase of 36 bpm). [NOTE: These positional vital signs were recorded "
            "by the medical assistant but were not reviewed by the referring physician before "
            "sending this referral.]\n\n"
            "PMH: Unremarkable. No prior psychiatric history.\n"
            "FH: Mother with fibromyalgia. Maternal aunt with lupus.\n"
            "Meds: None currently.\n"
            "Allergies: NKDA"
        ),
        (
            "1. Does this patient have generalized anxiety disorder or panic disorder?\n"
            "2. Would she benefit from an SSRI or anxiolytic?\n"
            "3. Any concerns for somatization disorder?"
        ),
        urgency="Routine",
    )

    # ──────────────────────────────────────────────────────────────
    # Psychiatry Evaluation (2021-09-13)
    # THE KEY DISMISSAL: POTS symptoms diagnosed as GAD + Panic Disorder
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-09-13_Psychiatry_Evaluation.pdf"),
        MAYA, "09/13/2021", "psychiatry",
        "PSYCHIATRIC EVALUATION -- NEW PATIENT",
        {
            "chief_complaint": (
                "Referred by PCP for anxiety evaluation. Recurrent episodes of tachycardia, "
                "dizziness, and near-syncope."
            ),
            "hpi": (
                "Ms. Rodriguez is a 25-year-old Hispanic female graduate student (neuroscience) "
                "referred by Dr. Martinez for evaluation of anxiety symptoms. She describes "
                "recurrent episodes of rapid heartbeat (she can feel her heart pounding), "
                "dizziness, lightheadedness, and a sensation of being about to faint. Episodes "
                "occur predominantly when standing up, standing in lines, or in warm environments. "
                "She had one syncopal episode in a grocery store in July 2021.\n\n"
                "The episodes are accompanied by:\n"
                "- Palpitations and tachycardia\n"
                "- Dizziness and tunnel vision\n"
                "- Nausea\n"
                "- Tremulousness\n"
                "- Feeling of impending doom\n"
                "- Difficulty concentrating ('brain fog')\n\n"
                "She has begun avoiding grocery stores, long lines, and standing events. She "
                "stopped going to the gym. She is increasingly worried about 'what's wrong with "
                "me' and has difficulty sleeping due to worry.\n\n"
                "She denies depressed mood. She denies suicidal or homicidal ideation. She denies "
                "substance use. No prior psychiatric treatment. No childhood trauma.\n\n"
                "She mentions that her symptoms are worse in the morning, with heat, after meals, "
                "and during her menstrual period. She feels better when lying down. She has noticed "
                "she needs to increase her salt and fluid intake to feel better. [NOTE: These are "
                "classic POTS features that were documented but not recognized as such.]\n\n"
                "She also mentions intermittent facial flushing, hives on her chest, and episodes "
                "of diarrhea alternating with constipation that she attributes to stress."
            ),
            "exam": (
                "MENTAL STATUS EXAMINATION:\n"
                "Appearance: Well-groomed, appropriately dressed 25-year-old woman.\n"
                "Behavior: Cooperative, good eye contact, mildly anxious.\n"
                "Speech: Normal rate, rhythm, and volume.\n"
                "Mood: 'Anxious and frustrated'\n"
                "Affect: Anxious but appropriate range, congruent with stated mood.\n"
                "Thought Process: Linear, goal-directed. No tangentiality or circumstantiality.\n"
                "Thought Content: No SI/HI. No delusions. No obsessions. Preoccupied with health "
                "concerns.\n"
                "Perceptions: No hallucinations.\n"
                "Cognition: Alert, oriented x4. Concentration intact.\n"
                "Insight: Good. Judgment: Good.\n\n"
                "PHYSICAL OBSERVATIONS (limited, not a full medical exam):\n"
                "Patient appeared flushed during interview. She mentioned feeling lightheaded when "
                "she stood from the waiting room chair. Hands appeared slightly swollen."
            ),
            "assessment": (
                "DSM-5 Diagnoses:\n\n"
                "1. Generalized Anxiety Disorder (F41.1)\n"
                "   - Excessive worry about health and daily activities > 6 months\n"
                "   - Difficulty controlling worry\n"
                "   - Associated symptoms: sleep disturbance, difficulty concentrating, "
                "muscle tension, fatigue\n\n"
                "2. Panic Disorder (F41.0)\n"
                "   - Recurrent unexpected panic attacks with prominent somatic symptoms "
                "(palpitations, dizziness, nausea, feeling of impending doom)\n"
                "   - Persistent concern about additional attacks\n"
                "   - Significant behavioral avoidance\n\n"
                "3. Agoraphobia (F40.00) -- mild\n"
                "   - Avoidance of situations where panic attacks have occurred "
                "(grocery stores, lines, gym)\n\n"
                "DIFFERENTIAL CONSIDERATIONS:\n"
                "- Medical causes of tachycardia have reportedly been ruled out per ER workup "
                "(labs, EKG normal). The positional nature of symptoms could suggest an "
                "autonomic component, but this is beyond the scope of psychiatric evaluation. "
                "If symptoms do not improve with psychiatric treatment, recommend medical "
                "re-evaluation with possible autonomic testing.\n\n"
                "NOTE: The somatic symptom pattern (GI complaints, flushing, brain fog) is "
                "consistent with autonomic hyperarousal in anxiety. However, if an underlying "
                "medical condition is later identified, these symptoms may need to be "
                "reattributed."
            ),
            "plan": (
                "1. Start sertraline 25 mg daily x 1 week, then increase to 50 mg daily\n"
                "   (SSRI for GAD + panic disorder)\n"
                "2. Lorazepam 0.5 mg PRN (max 2x/day) for acute panic episodes -- short-term "
                "bridge until SSRI takes effect\n"
                "3. CBT referral: cognitive-behavioral therapy for panic disorder\n"
                "4. Psychoeducation provided: nature of anxiety, fight-or-flight response, "
                "connection between physical symptoms and anxiety\n"
                "5. Sleep hygiene counseling\n"
                "6. Exercise encouraged (patient is avoiding due to fear -- graded exposure)\n"
                "7. Follow-up in 4 weeks to assess SSRI response\n"
                "8. If symptoms persist despite adequate SSRI trial (8-12 weeks), recommend "
                "medical re-evaluation for possible autonomic dysfunction"
            ),
        },
    )

    # ──────────────────────────────────────────────────────────────
    # Psychiatry Follow-up (2021-10-11)
    # Sertraline makes things WORSE -- classic for POTS/MCAS patients
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-10-11_Psychiatry_Follow_Up.pdf"),
        MAYA, "10/11/2021", "psychiatry",
        "PSYCHIATRIC FOLLOW-UP",
        {
            "chief_complaint": "Sertraline causing worsened dizziness and nausea. Panic episodes unchanged.",
            "hpi": (
                "Ms. Rodriguez returns 4 weeks after starting sertraline. She titrated to 50 mg "
                "as directed. She reports:\n"
                "- Dizziness has WORSENED since starting sertraline\n"
                "- Increased nausea, especially in the morning\n"
                "- One near-syncopal episode while standing in line at pharmacy\n"
                "- Facial flushing has increased in frequency\n"
                "- She developed a widespread itchy rash on her trunk after the first week, which "
                "she treated with Benadryl (resolved)\n"
                "- Sleep has improved slightly\n"
                "- Mood is more frustrated -- 'I keep being told it's anxiety but something "
                "doesn't feel right'\n"
                "- She used lorazepam 3 times for acute episodes -- it helped with the fear but "
                "NOT with the dizziness or heart racing\n\n"
                "She mentions she has been tracking her heart rate with a smartwatch. Standing "
                "heart rate regularly goes to 120-140 bpm. Resting heart rate is 65-75 bpm. "
                "She shows me screenshots. [NOTE: This data was acknowledged but not acted upon.]"
            ),
            "assessment": (
                "1. GAD / Panic Disorder -- partial response to sertraline. GI side effects "
                "are common with SSRI initiation and often resolve by week 6-8.\n"
                "2. The rash is likely a non-specific drug reaction or coincidental.\n"
                "3. Patient's smartwatch data is noted. HR variability with position change is "
                "not uncommon in anxiety states due to hyperadrenergic response.\n"
                "4. Patient is articulate and health-literate (neuroscience student) -- she has "
                "been researching her symptoms online and has expressed concern about 'POTS.' "
                "While I appreciate her engagement, internet-driven health anxiety can reinforce "
                "somatic hypervigilance. Discussed the importance of not 'doctor shopping' and "
                "allowing current treatment adequate time.\n\n"
                "NOTE: If symptoms persist after 8-12 weeks of adequate SSRI, I do recommend "
                "cardiology evaluation with tilt table testing to formally evaluate for "
                "autonomic dysfunction. The positional HR changes she reports are notable."
            ),
            "plan": (
                "1. Continue sertraline 50 mg -- GI side effects should improve\n"
                "2. Add hydroxyzine 25 mg PRN for acute anxiety (replacing lorazepam to avoid "
                "benzodiazepine dependence)\n"
                "3. Continue CBT referral (not yet started -- waitlist)\n"
                "4. Increase hydration and salt intake (general wellness)\n"
                "5. Follow-up in 6 weeks\n"
                "6. If no improvement by 12 weeks, discontinue SSRI and refer to Cardiology "
                "for tilt table testing"
            ),
        },
    )

    # ──────────────────────────────────────────────────────────────
    # Dermatology Consult (2022-08-15)
    # MCAS flushing diagnosed as rosacea
    # Between GI consult (June 2022) and Allergy consult (Oct 2022)
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-08-15_Derm_Consult.pdf"),
        MAYA, "08/15/2022", "dermatology",
        "DERMATOLOGY CONSULTATION -- FACIAL FLUSHING AND RASH",
        {
            "chief_complaint": "Recurrent facial flushing, hives on chest and neck x 18 months.",
            "hpi": (
                "Ms. Rodriguez is a 26-year-old woman referred by PCP for evaluation of recurrent "
                "facial flushing and urticarial rash. She reports episodic bright red flushing "
                "of face, neck, and upper chest occurring 4-5 times per week x 18 months. "
                "Episodes are triggered by:\n"
                "- Alcohol (even small amounts)\n"
                "- Spicy foods\n"
                "- Heat exposure / hot showers\n"
                "- Emotional stress\n"
                "- Exercise\n"
                "- Certain fragrances and cleaning products\n\n"
                "During flushing episodes, she also develops scattered urticarial wheals on her "
                "chest and neck that resolve within 2-4 hours. She sometimes feels her throat "
                "tighten slightly during severe episodes. She denies frank anaphylaxis.\n\n"
                "She has a history of POTS (diagnosed February 2022) and is being evaluated "
                "for GI complaints. She mentions her tryptase was elevated at 18.4 during an ER "
                "visit, but she's not sure what that means.\n\n"
                "No history of atopy. No childhood eczema. No food allergies diagnosed."
            ),
            "exam": (
                "Face: Diffuse erythema bilateral cheeks and nose. No papules or pustules. "
                "No telangiectasia. No rhinophyma.\n"
                "Chest: 3 urticarial wheals on upper chest (patient reports these appeared "
                "20 minutes ago in the warm waiting room). Faint dermographism noted.\n"
                "Neck: Mild erythema. No wheals currently.\n"
                "Hands: Hypermobile joints noted incidentally (fingers hyperextend past 90 degrees).\n"
                "No malar rash per se -- flushing is diffuse, not butterfly-pattern."
            ),
            "assessment": (
                "1. Rosacea, erythematotelangiectatic subtype (ETR) -- facial flushing with "
                "triggers including alcohol, heat, and spicy food is classic\n"
                "2. Chronic spontaneous urticaria -- intermittent wheals without clear allergic "
                "etiology\n\n"
                "NOTE: Patient mentions elevated tryptase level. This could be relevant to the "
                "flushing but is outside the scope of dermatologic evaluation. If flushing does "
                "not respond to rosacea treatment, consider Allergy/Immunology referral for "
                "mast cell evaluation."
            ),
            "plan": (
                "1. Metronidazole 0.75% cream to face BID\n"
                "2. Gentle cleanser, SPF 30+ daily\n"
                "3. Avoid known triggers (alcohol, heat, spicy food)\n"
                "4. Cetirizine 10 mg daily for urticaria\n"
                "5. If flushing persists despite treatment, consider Allergy referral\n"
                "6. Follow-up in 8 weeks"
            ),
        },
    )

    count_new = len([f for f in os.listdir(output_dir)
                     if f.endswith(".pdf") and any(x in f for x in
                     ["Psychiatry", "Referral_PCP_to_Psych", "Derm_Consult"])])
    print(f"    + Maya dismissals: {count_new} documents added")


def generate_sarah_dismissals(output_dir: str):
    """Generate dismissal documents for Sarah Mitchell."""
    os.makedirs(output_dir, exist_ok=True)

    # ──────────────────────────────────────────────────────────────
    # ER Visit (2023-02-14) -- Lupus flare dismissed as viral syndrome
    # Between derm consult (Jan 2023) and PCP follow-up (Apr 2023)
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-02-14_ER_Viral_Syndrome.pdf"),
        SARAH, "02/14/2023", "er",
        "EMERGENCY DEPARTMENT VISIT",
        {
            "chief_complaint": "Joint pain, facial rash, low-grade fever x 5 days.",
            "hpi": (
                "Ms. Mitchell is a 31-year-old woman presenting to the ED with worsening bilateral "
                "joint pain (hands, wrists, knees), a facial rash, and low-grade fever (100.4F at "
                "home) x 5 days. She reports extreme fatigue, unable to get out of bed today. "
                "She went to a dermatologist last month for a facial rash who 'took pictures and "
                "said it could be several things.' She had lab work recently showing 'some positive "
                "test' but has not followed up with her PCP yet.\n\n"
                "She came to the ER because the joint pain became severe overnight and she "
                "could not open a water bottle this morning. She also reports pleuritic chest "
                "pain (sharp, worse with deep breathing) x 2 days.\n\n"
                "No recent travel. No sick contacts. No cough, no URI symptoms. No GI symptoms. "
                "No vaginal discharge. No urinary symptoms."
            ),
            "vitals": (
                "BP 128/82  HR 96  Temp 100.8F  RR 18  SpO2 98% RA\n"
                "Pain: 8/10"
            ),
            "exam": (
                "General: Young woman in moderate distress due to joint pain. Appears fatigued.\n"
                "HEENT: Erythematous rash across bilateral malar eminences and bridge of nose, "
                "sparing nasolabial folds. Oral mucosa with 2 small painless ulcers on hard palate.\n"
                "Neck: Supple. Small bilateral cervical lymphadenopathy.\n"
                "Chest: Mild bilateral basilar crackles. No friction rub. Tenderness to palpation "
                "over left costochondral junction.\n"
                "CV: Tachycardic, regular rhythm. No murmur.\n"
                "Abdomen: Soft, non-tender.\n"
                "MSK: Bilateral MCP, PIP swelling and tenderness. Bilateral wrist swelling. "
                "Bilateral knee effusions, small. Reduced grip strength.\n"
                "Skin: Malar rash as above. Mottled erythema on extensor forearms (photosensitive "
                "distribution). Alopecia noted at temples.\n"
                "Neuro: Intact."
            ),
            "labs_reviewed": (
                "CBC: WBC 3.2 (L), Hgb 11.4 (L), Plt 132 (L), Lymph 0.8 (L)\n"
                "CMP: Cr 0.9, BUN 18, otherwise normal\n"
                "CRP: 24.6 (H)\n"
                "Urinalysis: Trace protein, no RBC, no casts [NOTE: Dipstick only, no microscopy]\n"
                "Chest X-ray: Small bilateral pleural effusions. No infiltrate.\n"
                "EKG: Sinus tachycardia, no ST changes."
            ),
            "assessment": (
                "1. Polyarthritis with fever -- likely viral arthritis vs. reactive arthritis vs. "
                "early inflammatory arthritis. Given her age and acute presentation with viral "
                "prodrome pattern, most likely VIRAL SYNDROME with arthritis.\n"
                "2. Facial rash -- likely viral exanthem vs. photosensitive reaction. Could be "
                "consistent with lupus but no autoimmune workup available tonight.\n"
                "3. Small bilateral pleural effusions -- likely reactive/viral.\n"
                "4. Leukopenia and lymphopenia -- consistent with viral illness.\n"
                "5. Pleuritic chest pain -- likely pleuritis from effusions.\n"
                "6. Oral ulcers -- noted, likely aphthous.\n\n"
                "NOTE: The combination of polyarthritis, malar rash, oral ulcers, pleuritis, "
                "leukopenia, and lymphopenia DOES raise the possibility of SLE. However, she has "
                "no prior history of autoimmune disease. ANA was not ordered in the ED as it is "
                "not part of the acute workup algorithm."
            ),
            "plan": (
                "1. Prednisone 20 mg x 5 days for symptomatic relief of polyarthritis\n"
                "2. Ibuprofen 600 mg TID with food PRN\n"
                "3. Diagnosis: VIRAL SYNDROME WITH POLYARTHRITIS (B34.9, M13.80)\n"
                "4. STRONGLY recommend PCP follow-up within 1 week for ANA testing and "
                "autoimmune workup given clinical features\n"
                "5. Return precautions: worsening fever, chest pain, shortness of breath, "
                "new rash, blood in urine\n"
                "6. Patient discharged in stable condition"
            ),
        },
    )

    # ER labs for that visit
    generate_lab_report(
        os.path.join(output_dir, "2023-02-14_Lab_ER.pdf"),
        SARAH, "02/14/2023", PROVIDERS["er"],
        [{
            "panel_name": "CBC WITH DIFFERENTIAL",
            "results": [
                {"test": "WBC", "value": "3.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": "L"},
                {"test": "Hemoglobin", "value": "11.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": "L"},
                {"test": "Hematocrit", "value": "34.2", "unit": "%", "ref_range": "36.0-46.0", "flag": "L"},
                {"test": "Platelets", "value": "132", "unit": "K/uL", "ref_range": "150-400", "flag": "L"},
                {"test": "Neutrophils, Abs", "value": "2.0", "unit": "K/uL", "ref_range": "1.8-7.7", "flag": ""},
                {"test": "Lymphocytes, Abs", "value": "0.8", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": "L"},
                {"test": "Monocytes, Abs", "value": "0.3", "unit": "K/uL", "ref_range": "0.2-0.8", "flag": ""},
            ],
        }, {
            "panel_name": "COMPREHENSIVE METABOLIC PANEL",
            "results": [
                {"test": "Glucose", "value": "94", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "BUN", "value": "18", "unit": "mg/dL", "ref_range": "7-20", "flag": ""},
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
                {"test": "Sodium", "value": "138", "unit": "mEq/L", "ref_range": "136-145", "flag": ""},
                {"test": "Potassium", "value": "4.0", "unit": "mEq/L", "ref_range": "3.5-5.1", "flag": ""},
                {"test": "CO2", "value": "24", "unit": "mEq/L", "ref_range": "22-29", "flag": ""},
                {"test": "Albumin", "value": "3.4", "unit": "g/dL", "ref_range": "3.5-5.5", "flag": "L"},
                {"test": "Total Protein", "value": "8.2", "unit": "g/dL", "ref_range": "6.0-8.3", "flag": ""},
            ],
        }, {
            "panel_name": "INFLAMMATORY / ACUTE",
            "results": [
                {"test": "CRP", "value": "24.6", "unit": "mg/L", "ref_range": "<3.0", "flag": "HH"},
                {"test": "Procalcitonin", "value": "0.08", "unit": "ng/mL", "ref_range": "<0.10", "flag": ""},
            ],
        }, {
            "panel_name": "URINALYSIS",
            "results": [
                {"test": "Appearance", "value": "Slightly cloudy", "unit": "", "ref_range": "Clear", "flag": ""},
                {"test": "Protein", "value": "Trace", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Blood", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Glucose", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Leukocyte Esterase", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }],
    )

    # ──────────────────────────────────────────────────────────────
    # Urgent Care Visit (2022-12-05) -- pleuritis dismissed as costochondritis
    # Before the derm consult (Jan 2023)
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-12-05_UC_Chest_Pain.pdf"),
        SARAH, "12/05/2022", "urgent_care",
        "URGENT CARE VISIT -- CHEST PAIN",
        {
            "chief_complaint": "Sharp left-sided chest pain x 3 days, worse with deep breathing.",
            "hpi": (
                "Ms. Mitchell is a 31-year-old woman presenting with sharp left-sided chest pain "
                "that started 3 days ago. Pain is worse with deep breathing, coughing, and lying "
                "flat. Improves slightly when sitting forward. She also reports joint pain in her "
                "hands and fatigue over the past few months. She has had a 'rash on her face' that "
                "comes and goes, worse with sun. She has no PCP appointment for another 3 weeks.\n\n"
                "No fever, no cough, no shortness of breath at rest. No recent travel. No leg "
                "swelling. No recent immobilization."
            ),
            "vitals": "BP 118/74  HR 88  Temp 99.2F  RR 16  SpO2 99% RA",
            "exam": (
                "Chest: Clear to auscultation. No friction rub appreciated. Reproducible "
                "tenderness over left parasternal region at 3rd-4th costal cartilage.\n"
                "CV: Regular rate and rhythm. No murmur.\n"
                "MSK: Mild swelling bilateral PIP joints. Mild tenderness bilateral wrists."
            ),
            "assessment": (
                "1. Costochondritis -- reproducible chest wall tenderness, pleuritic quality. "
                "Low probability for PE given low risk score.\n"
                "2. Arthralgias -- nonspecific. May be early OA vs. inflammatory. "
                "Recommend PCP follow-up."
            ),
            "plan": (
                "1. Naproxen 500 mg BID x 7 days\n"
                "2. If not improving, see PCP for further evaluation\n"
                "3. Return to ER for fever >101, worsening SOB, or hemoptysis"
            ),
        },
    )

    count_new = len([f for f in os.listdir(output_dir)
                     if f.endswith(".pdf") and any(x in f for x in
                     ["ER_Viral", "Lab_ER", "UC_Chest"])])
    print(f"    + Sarah dismissals: {count_new} documents added")


def generate_linda_dismissals(output_dir: str):
    """Generate dismissal documents for Linda Chen."""
    os.makedirs(output_dir, exist_ok=True)

    # ──────────────────────────────────────────────────────────────
    # PCP Visit (2021-03-08) -- Dry eyes/mouth attributed to perimenopause
    # Between ophth follow-up (Jan 2021) and next ophth (Aug 2021)
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-03-08_PCP_Perimenopause.pdf"),
        LINDA, "03/08/2021", "pcp",
        "OFFICE VISIT -- FATIGUE, DRY EYES, DRY MOUTH, JOINT STIFFNESS",
        {
            "chief_complaint": (
                "Fatigue x 3 months. Worsening dry eyes (under ophthalmology care). "
                "New dry mouth. Mild joint stiffness in the mornings."
            ),
            "hpi": (
                "Ms. Chen is a 47-year-old schoolteacher presenting with several complaints:\n\n"
                "1. FATIGUE: Progressive fatigue x 3 months. She sleeps 8 hours but wakes "
                "unrefreshed. Difficulty getting through the school day. Previously very active.\n\n"
                "2. DRY EYES: Under care of Ophthalmology (Dr. Wells) since June 2020. "
                "Currently on artificial tears and Restasis. Schirmer tests have been trending "
                "down (8 mm -> 6 mm over 6 months).\n\n"
                "3. DRY MOUTH: New onset x 4-5 months. Needs to sip water constantly while "
                "teaching. Difficulty swallowing dry foods. Wakes at night with dry mouth.\n\n"
                "4. JOINT STIFFNESS: Mild morning stiffness in hands and knees x 2 months, "
                "lasting 15-20 minutes, improves with activity.\n\n"
                "5. She also notes irregular menstrual periods over the past year and occasional "
                "hot flashes.\n\n"
                "PMH: Hypothyroidism (levothyroxine 75 mcg), mild dry eye\n"
                "Meds: Levothyroxine 75 mcg, Restasis, artificial tears, calcium + vitamin D\n"
                "Social: Teacher, married, 2 children (ages 18, 15). Non-smoker."
            ),
            "vitals": "BP 122/74  HR 72  Temp 98.4F  Wt 134 lbs",
            "exam": (
                "General: Well-appearing, no acute distress.\n"
                "HEENT: Oral mucosa appears mildly dry. No parotid enlargement appreciated "
                "(not specifically examined). Dentition intact.\n"
                "Thyroid: Non-enlarged, no nodules.\n"
                "MSK: Full ROM all joints. No synovitis. No tenderness.\n"
                "Skin: Normal turgor. No rash."
            ),
            "assessment": (
                "1. PERIMENOPAUSE -- age 47 with irregular menses, hot flashes, fatigue, and "
                "new-onset mucosal dryness. These are classic perimenopausal symptoms. "
                "Declining estrogen causes decreased secretory function including reduced "
                "tear and saliva production.\n"
                "2. Dry eyes -- multifactorial: perimenopause + screen time + environmental. "
                "Continue ophthalmology management.\n"
                "3. Dry mouth -- likely perimenopausal. Not on any anticholinergic medications.\n"
                "4. Morning joint stiffness -- common in perimenopause (estrogen decline "
                "affects joint lubrication). Brief duration (<30 min) and no synovitis "
                "argue against inflammatory arthritis.\n"
                "5. Fatigue -- perimenopausal. R/O anemia, thyroid dysfunction.\n\n"
                "NOTE: The constellation of dry eyes + dry mouth + fatigue + arthralgias "
                "COULD suggest Sjogren's syndrome. However, given her age and perimenopausal "
                "status, hormonal changes are the more likely explanation. Will check thyroid "
                "and basic labs. Autoimmune workup not indicated at this time unless symptoms "
                "progress significantly."
            ),
            "plan": (
                "1. TSH, Free T4 (recheck -- may need levothyroxine dose adjustment)\n"
                "2. CBC, CMP, vitamin D level\n"
                "3. Consider FSH if menstrual pattern clarification needed\n"
                "4. Biotene mouth rinse for xerostomia\n"
                "5. Continue ophthalmology follow-up for dry eyes\n"
                "6. OTC glucosamine/chondroitin for joint stiffness\n"
                "7. Exercise and stress management for fatigue\n"
                "8. Return in 3 months. If dry mouth/eyes worsen significantly, may "
                "consider ANA/SSA testing\n"
                "9. Discussed perimenopause -- patient may benefit from HRT evaluation if "
                "symptoms become bothersome"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2021-03-08_Lab_Perimenopause.pdf"),
        LINDA, "03/08/2021", PROVIDERS["pcp"],
        [{
            "panel_name": "THYROID / GENERAL",
            "results": [
                {"test": "TSH", "value": "2.8", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Free T4", "value": "1.2", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "28", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
                {"test": "FSH", "value": "28.4", "unit": "mIU/mL", "ref_range": "3.5-12.5 (follicular)", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "5.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.0", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "210", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Lymphocytes, Abs", "value": "1.2", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": ""},
            ],
        }],
    )

    count_new = len([f for f in os.listdir(output_dir)
                     if f.endswith(".pdf") and any(x in f for x in
                     ["Perimenopause"])])
    print(f"    + Linda dismissals: {count_new} documents added")


def generate_david_dismissals(output_dir: str):
    """Generate dismissal documents for David Park."""
    os.makedirs(output_dir, exist_ok=True)

    # ──────────────────────────────────────────────────────────────
    # Sleep Medicine Consult (2020-06-22) -- Insomnia workup for AS night pain
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2020-06-22_Sleep_Consult.pdf"),
        DAVID, "06/22/2020", "sleep_medicine",
        "SLEEP MEDICINE CONSULTATION",
        {
            "chief_complaint": "Insomnia and poor sleep quality x 1 year. Wakes multiple times at night.",
            "hpi": (
                "Mr. Park is a 39-year-old construction supervisor referred by PCP for evaluation "
                "of chronic insomnia. He reports difficulty maintaining sleep -- wakes 3-4 times "
                "per night. He falls asleep within 15-20 minutes but wakes after 2-3 hours and "
                "has difficulty returning to sleep. He rates sleep quality 3/10.\n\n"
                "He notes that waking is often associated with back pain and stiffness. He needs "
                "to move around or take a hot shower to get comfortable. Pain is worst in the "
                "second half of the night and early morning. He attributes his insomnia to 'chronic "
                "back pain from work.'\n\n"
                "He has tried melatonin 5 mg (minimal benefit), sleep hygiene changes, and a new "
                "mattress (no improvement). He does not snore per his wife. No witnessed apneas. "
                "BMI 27.4. No excessive daytime sleepiness (Epworth = 8/24, normal).\n\n"
                "He consumes 3-4 cups of coffee daily, stops by 2 PM. No alcohol before bed. "
                "No shift work."
            ),
            "exam": (
                "General: Alert, well-appearing male. BMI 27.4.\n"
                "HEENT: Mallampati II. Normal oropharynx. No retrognathia.\n"
                "Neck circumference: 16 inches (low risk for OSA).\n"
                "Chest: Clear. Normal expansion.\n"
                "Neuro: Alert, oriented. No focal deficits."
            ),
            "assessment": (
                "1. Chronic insomnia, sleep maintenance type -- likely secondary to chronic "
                "pain syndrome (low back pain). His sleep fragmentation is driven by pain "
                "awakening rather than primary insomnia.\n"
                "2. Low probability for obstructive sleep apnea (STOP-BANG score 1, low risk; "
                "no snoring, no witnessed apneas, not obese, normal Epworth).\n"
                "3. The pattern of waking in the second half of the night with stiffness that "
                "requires movement is somewhat atypical for mechanical back pain, which usually "
                "hurts MORE with movement. However, this is not my area of expertise.\n\n"
                "NOTE: Sleep study is not strongly indicated given low OSA risk, but patient "
                "requests one for reassurance. Will accommodate."
            ),
            "plan": (
                "1. Home sleep study (HST) to rule out OSA -- low pre-test probability\n"
                "2. Trazodone 50 mg at bedtime for sleep initiation and maintenance\n"
                "3. Continue addressing underlying back pain with PCP\n"
                "4. Sleep hygiene optimization reviewed\n"
                "5. CBT-I referral if insomnia persists after pain management\n"
                "6. Follow-up with sleep study results"
            ),
        },
    )

    # ──────────────────────────────────────────────────────────────
    # Sleep Study Results (2020-07-13) -- Normal, as expected
    # ──────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2020-07-13_Sleep_Study_Results.pdf"),
        DAVID, "07/13/2020", "sleep_medicine",
        "SLEEP STUDY RESULTS -- HOME SLEEP TEST",
        {
            "chief_complaint": "Follow-up home sleep study results.",
            "labs_reviewed": (
                "HOME SLEEP TEST (HST) -- 07/06/2020:\n\n"
                "Recording time: 7.2 hours\n"
                "Total sleep time (estimated): 5.8 hours\n"
                "Sleep efficiency: 80.6% (reduced)\n\n"
                "RESPIRATORY EVENTS:\n"
                "AHI (Apnea-Hypopnea Index): 2.1 events/hour (NORMAL; <5 = normal)\n"
                "Obstructive apneas: 4\n"
                "Hypopneas: 8\n"
                "Central apneas: 0\n"
                "Lowest SpO2: 93%\n"
                "Mean SpO2: 96%\n"
                "Time SpO2 <90%: 0%\n\n"
                "BODY POSITION:\n"
                "Supine: 45% of recording time\n"
                "Non-supine: 55%\n"
                "Supine AHI: 3.2  |  Non-supine AHI: 1.4\n\n"
                "HEART RATE:\n"
                "Mean: 68 bpm  |  Min: 52 bpm  |  Max: 94 bpm\n"
                "Multiple position changes noted during recording (12 major position changes) -- "
                "suggests significant sleep fragmentation from discomfort."
            ),
            "assessment": (
                "1. NO OBSTRUCTIVE SLEEP APNEA -- AHI 2.1 (normal)\n"
                "2. Reduced sleep efficiency (80.6%) with frequent position changes -- "
                "consistent with pain-related sleep fragmentation\n"
                "3. The underlying issue is his chronic back pain, not a primary sleep disorder"
            ),
            "plan": (
                "1. No CPAP indicated\n"
                "2. Continue trazodone 50 mg for sleep maintenance\n"
                "3. Primary management should focus on back pain treatment\n"
                "4. Recommend he continue following up with PCP regarding back pain -- "
                "the pain pattern (worse with rest, better with movement) is unusual and "
                "may warrant further evaluation\n"
                "5. Discharged from sleep medicine"
            ),
        },
    )

    count_new = len([f for f in os.listdir(output_dir)
                     if f.endswith(".pdf") and any(x in f for x in
                     ["Sleep"])])
    print(f"    + David dismissals: {count_new} documents added")
