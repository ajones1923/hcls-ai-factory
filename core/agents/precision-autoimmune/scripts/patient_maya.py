"""Patient 2: Maya Rodriguez -- 28F, 4-year POTS/hEDS/MCAS diagnostic odyssey.

Timeline: March 2021 -> January 2025
Specialists seen: PCP, ER (x3), Cardiology, Neurology, GI, Allergy
Key pattern: Orthostatic HR spikes buried in ER vitals, joint hypermobility
noted but never scored, elevated tryptase dismissed as "lab error."

The system should detect: POTS + hEDS + MCAS triad, TPSAB1 duplication.
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
    "name": "Maya Rodriguez",
    "dob": "1997-03-22",
    "mrn": "MRO-2021-55203",
    "sex": "F",
    "age_at_start": 23,
}


def generate(output_dir: str):
    """Generate all clinical documents for Maya Rodriguez."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: PCP -- March 2021 (fatigue, dizziness)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-03-08_PCP_Visit.pdf"),
        PATIENT, "03/08/2021", "pcp",
        "OFFICE VISIT -- FATIGUE AND DIZZINESS",
        {
            "chief_complaint": "Fatigue, lightheadedness, and intermittent palpitations x 2 months.",
            "hpi": (
                "Ms. Rodriguez is a 23-year-old graduate student presenting with progressive fatigue "
                "and lightheadedness over the past 2 months. She describes feeling dizzy when standing "
                "up from sitting or lying down, sometimes seeing 'stars' or experiencing tunnel vision. "
                "She has had several near-syncopal episodes but no frank syncope. She reports palpitations -- "
                "a racing heart -- that occur with standing and sometimes at rest. Symptoms are worse in "
                "hot weather, after meals, and with prolonged standing. She has tried increasing fluid "
                "intake without significant improvement.\n\n"
                "She also reports chronic joint pain, especially in her fingers, wrists, and knees, which "
                "she has had 'for years' and attributed to being 'double-jointed.' She was a dancer in "
                "high school and had multiple ankle sprains and one shoulder subluxation.\n\n"
                "GI symptoms: intermittent abdominal cramping, bloating, alternating diarrhea/constipation "
                "x 6 months. She has noticed flushing episodes and occasional hives without clear trigger.\n\n"
                "No fever, weight loss, chest pain, or dyspnea. No family history of heart disease."
            ),
            "vitals": (
                "Sitting: BP 112/68  HR 72\n"
                "Standing (1 min): BP 108/70  HR 98\n"
                "Standing (3 min): BP 106/72  HR 108\n"
                "Temp 98.2F  SpO2 99%  Wt 118 lbs  Ht 5'4\""
            ),
            "ros": (
                "Constitutional: Fatigue, lightheadedness as above. No fever or weight loss.\n"
                "Cardiovascular: Palpitations with standing. No chest pain. No edema.\n"
                "GI: Bloating, cramping, alternating bowel habits. No blood in stool.\n"
                "Musculoskeletal: Chronic joint hypermobility, frequent subluxations. Chronic pain.\n"
                "Skin: Intermittent flushing and urticaria. Easy bruising.\n"
                "Neurological: Brain fog, difficulty concentrating. No headaches, no seizures."
            ),
            "exam": (
                "General: Thin young woman, appears anxious. No acute distress.\n"
                "Cardiovascular: Tachycardic when standing (confirmed). Regular rhythm. No murmurs.\n"
                "Pulmonary: Clear bilaterally.\n"
                "Abdomen: Soft, mild diffuse tenderness. No organomegaly.\n"
                "Musculoskeletal: Notable joint hypermobility -- hyperextension at elbows, knees, "
                "and MCPs. Can place palms flat on floor with knees straight. Thumbs touch forearms. "
                "Skin: Soft, velvety texture. Mild striae on thighs. Easy bruising on arms.\n"
                "Neurological: Alert, oriented. Cranial nerves intact."
            ),
            "assessment": (
                "1. Orthostatic tachycardia -- HR increase of 36 bpm on standing. Differential "
                "includes dehydration, anxiety, anemia, or less likely POTS.\n"
                "2. Fatigue and brain fog -- likely related to orthostatic symptoms\n"
                "3. Joint hypermobility -- lifelong, likely benign hypermobility spectrum\n"
                "4. GI symptoms -- IBS suspected\n"
                "5. Intermittent urticaria/flushing -- allergic vs. idiopathic"
            ),
            "plan": (
                "1. CBC, TSH, CMP, iron studies to rule out anemia/thyroid\n"
                "2. Increase fluid intake to 2-3 L/day, increase salt intake\n"
                "3. If orthostatic symptoms persist, will refer to Cardiology for tilt table test\n"
                "4. Trial of omeprazole 20 mg for GI symptoms\n"
                "5. OTC cetirizine 10 mg daily for urticaria\n"
                "6. Return in 4-6 weeks"
            ),
        },
    )

    # Lab: March 2021 -- normal basic workup
    generate_lab_report(
        os.path.join(output_dir, "2021-03-08_Lab_Basic.pdf"),
        PATIENT, "03/08/2021", PROVIDERS["pcp"],
        [{
            "panel_name": "COMPLETE BLOOD COUNT",
            "results": [
                {"test": "WBC", "value": "6.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.2", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "245", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }, {
            "panel_name": "METABOLIC / THYROID",
            "results": [
                {"test": "TSH", "value": "1.8", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Glucose", "value": "82", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "Creatinine", "value": "0.7", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "Iron", "value": "72", "unit": "ug/dL", "ref_range": "60-170", "flag": ""},
                {"test": "Ferritin", "value": "38", "unit": "ng/mL", "ref_range": "12-150", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: ER -- July 2021 (syncope at work)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-07-14_ER_Syncope.pdf"),
        PATIENT, "07/14/2021", "er",
        "EMERGENCY DEPARTMENT -- SYNCOPE",
        {
            "chief_complaint": "Syncopal episode at work. Found on floor by coworker.",
            "hpi": (
                "24-year-old woman brought in by coworker after witnessed syncopal episode at work. "
                "Patient was standing in line at cafeteria, felt lightheaded and warm, then lost "
                "consciousness for approximately 10-15 seconds. No head strike. No seizure activity "
                "witnessed. No tongue biting or incontinence. Patient reports she had only a small "
                "breakfast and was in a warm, crowded cafeteria.\n\n"
                "She reports ongoing symptoms of lightheadedness with standing, palpitations, and "
                "fatigue x 6 months. She was told by PCP to drink more water. She states she drinks "
                "'about 1.5 liters a day.' She denies alcohol, drugs, or new medications."
            ),
            "vitals": (
                "Triage: BP 108/62  HR 68 (supine)\n"
                "Repeat standing: BP 104/64  HR 112 (standing 2 min)  **HR increase 44 bpm**\n"
                "Temp 98.4F  SpO2 99%"
            ),
            "exam": (
                "General: Alert, oriented x4. Appears well. No injuries from fall.\n"
                "HEENT: No head laceration. PERRL. No tongue bite.\n"
                "Cardiovascular: Tachycardic when upright, regular when supine. No murmurs.\n"
                "Neurological: Alert, oriented. GCS 15. No focal deficits. Normal gait."
            ),
            "assessment": (
                "1. Syncope -- likely vasovagal given prodrome and triggers (warm environment, "
                "standing, inadequate hydration). Orthostatic vitals notable for significant "
                "tachycardia (HR increase 44 bpm) without significant BP drop.\n"
                "2. ECG: normal sinus rhythm, no prolonged QT, no pre-excitation, no Brugada pattern."
            ),
            "plan": (
                "1. ECG -- normal (see results)\n"
                "2. BMP -- normal\n"
                "3. Diagnosis: vasovagal syncope\n"
                "4. Discharge with instructions: increase fluids, avoid prolonged standing, "
                "counter-pressure maneuvers\n"
                "5. Follow-up with PCP within 1 week\n"
                "6. If recurrent syncope, consider Holter monitor or cardiology referral"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: ER -- November 2021 (tachycardia, flushing, abdominal pain)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-11-22_ER_Tachycardia.pdf"),
        PATIENT, "11/22/2021", "er",
        "EMERGENCY DEPARTMENT -- TACHYCARDIA AND FLUSHING",
        {
            "chief_complaint": "Acute episode of tachycardia, facial flushing, abdominal cramping, and diffuse hives.",
            "hpi": (
                "24-year-old woman presents with acute onset tachycardia (patient reports HR 130s on "
                "home monitor), diffuse flushing, urticaria (hives) on trunk and extremities, and severe "
                "abdominal cramping with watery diarrhea. Episode started approximately 2 hours ago "
                "without clear trigger. She ate a normal dinner. No new exposures, no insect stings, "
                "no new medications. She reports these episodes have happened 'maybe 5-6 times in the "
                "past year' but this is the worst one. Previous episodes typically resolved in 1-2 hours. "
                "She takes cetirizine 10 mg daily which helps somewhat.\n\n"
                "No dyspnea, no throat swelling, no wheezing. No hypotension."
            ),
            "vitals": "BP 118/72  HR 118  Temp 99.0F  SpO2 98%  RR 18",
            "exam": (
                "General: Anxious-appearing, flushed. Diffuse urticaria on trunk and arms.\n"
                "HEENT: No angioedema. No throat swelling. Facial flushing.\n"
                "Cardiovascular: Tachycardic, regular. No murmurs.\n"
                "Pulmonary: Clear. No wheezing.\n"
                "Abdomen: Diffuse tenderness, hyperactive bowel sounds. No guarding.\n"
                "Skin: Diffuse urticarial wheals, 1-3 cm, trunk and extremities. No petechiae."
            ),
            "assessment": (
                "1. Acute urticaria with systemic symptoms (flushing, tachycardia, GI distress) -- "
                "differential includes allergic reaction vs. mast cell activation episode\n"
                "2. No anaphylaxis criteria met (no hypotension, no respiratory compromise)\n"
                "3. Recurrent episodes concerning for chronic urticaria vs. MCAS vs. carcinoid (less likely)"
            ),
            "plan": (
                "1. Diphenhydramine 50 mg IV, famotidine 20 mg IV\n"
                "2. CBC, CMP, tryptase level (drawn during episode)\n"
                "3. Symptoms improving after antihistamines\n"
                "4. Discharge on famotidine 20 mg BID + cetirizine 10 mg BID (increase from daily)\n"
                "5. Follow-up with Allergy/Immunology for recurrent urticaria workup\n"
                "6. If future episodes involve throat swelling or breathing difficulty, use EpiPen and call 911"
            ),
        },
    )

    # Lab: ER tryptase -- elevated but dismissed
    generate_lab_report(
        os.path.join(output_dir, "2021-11-22_Lab_ER_Tryptase.pdf"),
        PATIENT, "11/22/2021", PROVIDERS["er"],
        [{
            "panel_name": "EMERGENCY LABS",
            "results": [
                {"test": "Tryptase, Serum", "value": "18.4", "unit": "ng/mL", "ref_range": "<11.4", "flag": "H"},
                {"test": "WBC", "value": "7.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Glucose", "value": "92", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "Creatinine", "value": "0.7", "unit": "mg/dL", "ref_range": "0.6-1.2", "flag": ""},
                {"test": "Potassium", "value": "3.8", "unit": "mEq/L", "ref_range": "3.5-5.0", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: Cardiology -- February 2022
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-02-07_Cardiology_Consult.pdf"),
        PATIENT, "02/07/2022", "cardiology",
        "CARDIOLOGY CONSULTATION -- RECURRENT SYNCOPE AND TACHYCARDIA",
        {
            "chief_complaint": "Referred for recurrent syncope and palpitations. One syncopal episode (July 2021), "
                              "ongoing pre-syncopal episodes.",
            "hpi": (
                "Ms. Rodriguez is a 24-year-old woman referred for evaluation of recurrent near-syncope "
                "and tachycardia. She has had one witnessed syncopal episode (July 2021, ER visit) and "
                "reports daily lightheadedness with standing, racing heart, and exercise intolerance. "
                "Symptoms are worst in the morning, after meals, and in warm environments. She has tried "
                "increasing fluids and salt per PCP recommendations with minimal improvement.\n\n"
                "She also reports an ER visit in November 2021 for flushing, hives, and tachycardia -- "
                "treated with antihistamines. She takes cetirizine and famotidine daily.\n\n"
                "No family history of sudden cardiac death, long QT, or cardiomyopathy."
            ),
            "vitals": (
                "Supine: BP 114/68  HR 66\n"
                "Seated (5 min): BP 112/70  HR 78\n"
                "Standing (1 min): BP 110/72  HR 100\n"
                "Standing (5 min): BP 108/74  HR 112\n"
                "Standing (10 min): BP 106/70  HR 118  **HR increase 52 bpm**"
            ),
            "exam": (
                "General: Thin, well-appearing. Becomes visibly uncomfortable after standing 5+ minutes.\n"
                "Cardiovascular: Regular rate and rhythm when supine. Tachycardic when standing. "
                "No murmurs. No S3/S4. JVP normal.\n"
                "Pulmonary: Clear.\n"
                "Extremities: No edema. Mild dependent acrocyanosis of feet when standing.\n"
                "Skin: No rash today."
            ),
            "assessment": (
                "1. Orthostatic tachycardia -- HR increase of 52 bpm on 10-minute active stand test. "
                "This exceeds the 30 bpm threshold. No orthostatic hypotension. Symptoms have been "
                "present >6 months.\n"
                "2. This is consistent with POTS (postural orthostatic tachycardia syndrome), though "
                "I would like to confirm with formal tilt table testing.\n"
                "3. Dependent acrocyanosis is commonly seen with POTS.\n"
                "4. Echocardiogram recommended to rule out structural heart disease."
            ),
            "plan": (
                "1. Formal tilt table test -- scheduled for March 2022\n"
                "2. Echocardiogram today\n"
                "3. Holter monitor x 48 hours\n"
                "4. Continue increased fluids (2.5-3L/day) and salt (3-5g added/day)\n"
                "5. Consider compression stockings (waist-high, 20-30 mmHg)\n"
                "6. If tilt table confirms POTS, will discuss pharmacotherapy (fludrocortisone, midodrine, "
                "or ivabradine)\n"
                "7. Follow-up after tilt table test"
            ),
        },
    )

    # Echo -- normal
    generate_imaging_report(
        os.path.join(output_dir, "2022-02-07_Echo.pdf"),
        PATIENT, "02/07/2022", "cardiology",
        "Echocardiogram", "Transthoracic",
        "Recurrent syncope, orthostatic tachycardia, rule out structural heart disease.",
        "Standard 2D, M-mode, and Doppler echocardiogram performed.",
        (
            "Left ventricle: Normal size and systolic function. LVEF 62% by biplane Simpson's method. "
            "No wall motion abnormalities. Normal diastolic function.\n"
            "Right ventricle: Normal size and function. TAPSE 2.2 cm.\n"
            "Valves: All valves structurally normal. Trivial TR. No significant regurgitation or stenosis.\n"
            "Atria: Normal left and right atrial size.\n"
            "Aortic root: Normal.\n"
            "Pericardium: Normal. No effusion.\n"
            "IVC: Normal size with normal respiratory variation."
        ),
        "Normal transthoracic echocardiogram. No structural heart disease identified.",
        reading_radiologist="Michael Thompson, MD, FACC",
    )

    # Tilt table -- March 2022
    generate_progress_note(
        os.path.join(output_dir, "2022-03-14_Tilt_Table_Test.pdf"),
        PATIENT, "03/14/2022", "cardiology",
        "TILT TABLE TEST REPORT",
        {
            "chief_complaint": "Formal tilt table test for suspected POTS.",
            "hpi": "See previous cardiology consultation (02/07/2022). Symptoms unchanged.",
            "vitals": (
                "TILT TABLE TEST RESULTS:\n"
                "Supine baseline (10 min): HR 64, BP 116/72\n"
                "Tilt 70 degrees -- Minute 1: HR 88, BP 114/74\n"
                "Tilt 70 degrees -- Minute 3: HR 102, BP 112/76\n"
                "Tilt 70 degrees -- Minute 5: HR 114, BP 110/74\n"
                "Tilt 70 degrees -- Minute 8: HR 122, BP 108/72\n"
                "Tilt 70 degrees -- Minute 10: HR 128, BP 106/70  **HR increase 64 bpm**\n"
                "Tilt 70 degrees -- Minute 12: HR 132, BP 104/68\n"
                "Patient reported lightheadedness, palpitations, and nausea at minute 8.\n"
                "Near-syncopal at minute 12 -- test terminated. Returned to supine with rapid HR recovery."
            ),
            "assessment": (
                "POSITIVE TILT TABLE TEST FOR POTS\n"
                "Heart rate increase of 64 bpm within 10 minutes of head-up tilt, without orthostatic "
                "hypotension (BP drop <20/10). Symptoms reproduced during test (lightheadedness, "
                "palpitations, nausea). This meets diagnostic criteria for Postural Orthostatic "
                "Tachycardia Syndrome (POTS).\n\n"
                "No vasovagal response (no cardioinhibitory or vasodepressor response).\n"
                "Classification: Likely neuropathic POTS given the pattern."
            ),
            "plan": (
                "1. Diagnosis confirmed: POTS\n"
                "2. Start fludrocortisone 0.1 mg daily for volume expansion\n"
                "3. Continue increased salt and fluid intake\n"
                "4. Waist-high compression stockings 20-30 mmHg\n"
                "5. Graded exercise program (recumbent exercise initially -- rowing, swimming, recumbent bike)\n"
                "6. Avoid prolonged standing, hot environments, and large carbohydrate-heavy meals\n"
                "7. Follow-up in 6 weeks to assess response\n"
                "8. Note: I am not pursuing the cause of her POTS at this time. Common etiologies include "
                "deconditioning, autoimmune, and hypermobility-associated. Her joint hypermobility is "
                "notable but I will defer further evaluation to her PCP or Rheumatology if warranted."
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: GI -- June 2022 (IBS diagnosis)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-06-20_GI_Consult.pdf"),
        PATIENT, "06/20/2022", "gi",
        "GASTROENTEROLOGY CONSULTATION -- CHRONIC GI SYMPTOMS",
        {
            "chief_complaint": "Chronic bloating, abdominal cramping, alternating diarrhea/constipation x 1 year.",
            "hpi": (
                "Ms. Rodriguez is a 25-year-old woman with recently diagnosed POTS (March 2022) referred "
                "for chronic GI symptoms. She reports daily bloating, intermittent cramping abdominal pain, "
                "and alternating loose stools and constipation. Symptoms are worse after meals, especially "
                "large meals. She also reports intermittent nausea, early satiety, and occasional reflux. "
                "She takes omeprazole 20 mg daily and famotidine 20 mg BID (for MCAS symptoms per ER). "
                "No blood in stool, no weight loss, no dysphagia. No family history of IBD, celiac, or "
                "colon cancer."
            ),
            "exam": (
                "Abdomen: Soft, mild diffuse tenderness without guarding or rebound. Mildly distended. "
                "Active bowel sounds. No masses. No hepatosplenomegaly."
            ),
            "assessment": (
                "1. Chronic GI symptoms consistent with irritable bowel syndrome (IBS), mixed type\n"
                "2. Gastroparesis should be considered given POTS diagnosis (frequently comorbid)\n"
                "3. Screen for celiac disease given symptoms"
            ),
            "plan": (
                "1. Anti-tTG IgA and total IgA to screen for celiac disease\n"
                "2. CBC, CMP, CRP, fecal calprotectin\n"
                "3. If celiac serologies negative, empiric low-FODMAP diet trial x 6 weeks\n"
                "4. If gastroparesis suspected, consider gastric emptying study\n"
                "5. Continue PPI and H2 blocker\n"
                "6. Follow-up in 6 weeks with results"
            ),
        },
    )

    # GI labs -- celiac negative, calprotectin normal
    generate_lab_report(
        os.path.join(output_dir, "2022-06-20_Lab_GI.pdf"),
        PATIENT, "06/20/2022", PROVIDERS["gi"],
        [{
            "panel_name": "CELIAC SCREENING / GI PANEL",
            "results": [
                {"test": "Anti-tTG IgA", "value": "3", "unit": "U/mL", "ref_range": "<20", "flag": ""},
                {"test": "Total IgA", "value": "185", "unit": "mg/dL", "ref_range": "70-400", "flag": ""},
                {"test": "Fecal Calprotectin", "value": "28", "unit": "ug/g", "ref_range": "<50", "flag": ""},
                {"test": "CRP", "value": "1.8", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Allergy -- October 2022 (tryptase follow-up)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-10-17_Allergy_Consult.pdf"),
        PATIENT, "10/17/2022", "allergy",
        "ALLERGY/IMMUNOLOGY CONSULTATION -- RECURRENT URTICARIA AND FLUSHING",
        {
            "chief_complaint": "Recurrent episodes of flushing, urticaria, and GI distress. Elevated tryptase "
                              "noted in ER (18.4 ng/mL, November 2021).",
            "hpi": (
                "Ms. Rodriguez is a 25-year-old woman with POTS (diagnosed March 2022) referred for "
                "evaluation of recurrent flushing, urticaria, and GI symptoms. She reports approximately "
                "monthly episodes of facial/upper body flushing, diffuse hives, abdominal cramping, and "
                "diarrhea. Episodes last 1-4 hours and partially respond to dual antihistamines "
                "(cetirizine + famotidine). A tryptase level drawn during an acute episode in the ER "
                "(November 2021) was elevated at 18.4 ng/mL (ref <11.4). No anaphylaxis history. "
                "No known specific allergies by skin testing.\n\n"
                "She reports the POTS is partially controlled with fludrocortisone, salt, and fluids. "
                "She continues to have joint hypermobility and chronic pain."
            ),
            "exam": (
                "General: Well-appearing. No acute flushing or urticaria today.\n"
                "Skin: No active hives. Mild dermatographism (positive dermatographic test).\n"
                "No mastocytosis-type skin lesions (no urticaria pigmentosa/maculopapular cutaneous mastocytosis)."
            ),
            "assessment": (
                "1. Recurrent episodic flushing, urticaria, and GI symptoms with elevated acute tryptase -- "
                "raises concern for mast cell activation syndrome (MCAS)\n"
                "2. However, the single elevated tryptase could represent an acute allergic reaction. "
                "Need to establish baseline tryptase and confirm with additional mediators.\n"
                "3. No skin lesions to suggest systemic mastocytosis\n"
                "4. The combination of POTS + joint hypermobility + mast cell symptoms is recognized -- "
                "some literature describes a POTS/hEDS/MCAS triad"
            ),
            "plan": (
                "1. Baseline tryptase (when well, not during episode)\n"
                "2. 24-hour urine for N-methylhistamine, prostaglandin D2 metabolite (11-beta-PGF2alpha), "
                "and leukotriene E4\n"
                "3. If baseline tryptase elevated, consider bone marrow biopsy to rule out systemic mastocytosis\n"
                "4. Continue H1 (cetirizine 10 mg BID) + H2 (famotidine 20 mg BID) antihistamines\n"
                "5. Add cromolyn sodium 200 mg QID before meals for GI mast cell stabilization\n"
                "6. Epinephrine auto-injector prescribed as precaution\n"
                "7. Keep symptom diary with triggers\n"
                "8. Follow-up in 6 weeks with results"
            ),
        },
    )

    # Allergy labs -- baseline tryptase still elevated
    generate_lab_report(
        os.path.join(output_dir, "2022-10-17_Lab_MCAS_Workup.pdf"),
        PATIENT, "10/17/2022", PROVIDERS["allergy"],
        [{
            "panel_name": "MAST CELL MEDIATORS (BASELINE, ASYMPTOMATIC)",
            "results": [
                {"test": "Tryptase, Serum (baseline)", "value": "14.2", "unit": "ng/mL", "ref_range": "<11.4", "flag": "H"},
                {"test": "N-Methylhistamine, 24hr Urine", "value": "285", "unit": "mcg/g Cr", "ref_range": "30-200", "flag": "H"},
                {"test": "Prostaglandin D2 (11-beta-PGF2a)", "value": "1850", "unit": "ng/24hr", "ref_range": "100-1000", "flag": "H"},
                {"test": "Leukotriene E4, 24hr Urine", "value": "124", "unit": "pg/mg Cr", "ref_range": "<104", "flag": "H"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 7: Neurology -- March 2023 (brain fog, autonomic testing)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-03-06_Neuro_Consult.pdf"),
        PATIENT, "03/06/2023", "neurology",
        "NEUROLOGY CONSULTATION -- AUTONOMIC DYSFUNCTION, BRAIN FOG",
        {
            "chief_complaint": "Evaluation of POTS with brain fog and suspected small fiber neuropathy.",
            "hpi": (
                "Ms. Rodriguez is a 25-year-old woman with confirmed POTS (tilt table March 2022), "
                "suspected MCAS (elevated tryptase, elevated urine mediators), and joint hypermobility. "
                "Referred for autonomic evaluation and assessment of brain fog and possible small fiber "
                "neuropathy. She reports burning pain in her feet, numbness in fingertips, temperature "
                "dysregulation (sweating abnormalities), and significant cognitive difficulty ('brain fog') "
                "that interferes with her graduate studies. She was previously an honors student.\n\n"
                "Current medications: fludrocortisone 0.1 mg daily, cetirizine 10 mg BID, famotidine "
                "20 mg BID, cromolyn sodium 200 mg QID, omeprazole 20 mg daily."
            ),
            "exam": (
                "Neurological:\n"
                "Mental status: Alert, oriented x4. Difficulty with serial 7s (self-reported cognitive fog).\n"
                "Cranial nerves: II-XII intact.\n"
                "Motor: 5/5 throughout. No atrophy.\n"
                "Sensory: Decreased pinprick sensation in stocking distribution (feet to mid-calf bilaterally). "
                "Decreased temperature discrimination in feet. Vibration and proprioception intact.\n"
                "Reflexes: 2+ throughout. No pathologic reflexes.\n"
                "Coordination: Normal finger-to-nose and heel-to-shin.\n"
                "Gait: Normal. Tandem gait intact.\n\n"
                "Autonomic: Dependent acrocyanosis when standing. Reduced sweating on feet by exam."
            ),
            "assessment": (
                "1. POTS with possible autoimmune and/or neuropathic etiology\n"
                "2. Small fiber neuropathy suspected -- burning pain, sensory changes in stocking "
                "distribution. Consider skin biopsy (IENFD) for confirmation.\n"
                "3. Autonomic neuropathy -- QSART and autonomic reflex testing recommended\n"
                "4. Brain fog -- likely secondary to cerebral hypoperfusion from POTS\n"
                "5. Consider autoimmune etiologies: ganglionic AChR antibodies, Sjogren's panel"
            ),
            "plan": (
                "1. Autonomic reflex testing (QSART, Valsalva, deep breathing, tilt)\n"
                "2. Skin biopsy for intraepidermal nerve fiber density (IENFD) -- bilateral distal leg and proximal thigh\n"
                "3. Ganglionic AChR antibodies to evaluate for autoimmune autonomic ganglionopathy\n"
                "4. Sjogren's panel (anti-SSA/Ro, anti-SSB/La) -- Sjogren's is a treatable cause of autonomic neuropathy\n"
                "5. NCS/EMG to evaluate for large fiber neuropathy\n"
                "6. B12, methylmalonic acid, A1c to rule out metabolic neuropathy\n"
                "7. Follow-up with results in 6 weeks"
            ),
        },
    )

    # Neurology labs -- ganglionic AChR negative, Sjogren's negative, B12 normal
    generate_lab_report(
        os.path.join(output_dir, "2023-03-06_Lab_Neuro.pdf"),
        PATIENT, "03/06/2023", PROVIDERS["neurology"],
        [{
            "panel_name": "AUTOIMMUNE / NEUROPATHY PANEL",
            "results": [
                {"test": "Ganglionic AChR Ab", "value": "<0.02", "unit": "nmol/L", "ref_range": "<0.02", "flag": ""},
                {"test": "Anti-SSA/Ro", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-SSB/La", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "ANA", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "B12", "value": "485", "unit": "pg/mL", "ref_range": "200-900", "flag": ""},
                {"test": "Methylmalonic Acid", "value": "0.12", "unit": "umol/L", "ref_range": "<0.40", "flag": ""},
                {"test": "Hemoglobin A1c", "value": "5.1", "unit": "%", "ref_range": "<5.7", "flag": ""},
            ],
        }],
    )

    # IENFD biopsy result -- reduced
    generate_pathology_report(
        os.path.join(output_dir, "2023-04-10_Path_Skin_Biopsy_IENFD.pdf"),
        PATIENT, "04/10/2023",
        "Skin, bilateral distal leg and proximal thigh, punch biopsy for IENFD",
        "26-year-old woman with POTS, burning pain in feet, suspected small fiber neuropathy.",
        "Four 3mm punch biopsies received: two from distal leg (10 cm above lateral malleolus, bilateral) "
        "and two from proximal thigh (20 cm below iliac spine, bilateral).",
        (
            "Specimens processed for PGP9.5 immunohistochemistry per standardized IENFD protocol.\n\n"
            "Right distal leg: IENFD = 2.8 fibers/mm (reference: >7.8 fibers/mm for age/sex) -- REDUCED\n"
            "Left distal leg: IENFD = 3.1 fibers/mm -- REDUCED\n"
            "Right proximal thigh: IENFD = 8.2 fibers/mm (reference: >11.2 fibers/mm) -- REDUCED\n"
            "Left proximal thigh: IENFD = 7.9 fibers/mm -- REDUCED\n\n"
            "Length-dependent pattern of nerve fiber loss confirmed."
        ),
        (
            "SMALL FIBER NEUROPATHY -- CONFIRMED\n"
            "Reduced intraepidermal nerve fiber density at all four sites with length-dependent "
            "gradient, confirming small fiber neuropathy.\n"
            "IENFD values are significantly below age/sex-matched normative values (5th percentile cutoffs)."
        ),
        comment=(
            "Small fiber neuropathy is confirmed by IENFD testing. This finding is consistent with "
            "the clinical presentation of burning pain, autonomic dysfunction (POTS), and temperature "
            "dysregulation. Etiological workup should include evaluation for autoimmune, metabolic, "
            "and genetic causes."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 8: PCP -- July 2023 (Beighton score, hEDS evaluation)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-07-24_PCP_hEDS_Eval.pdf"),
        PATIENT, "07/24/2023", "pcp",
        "OFFICE VISIT -- JOINT HYPERMOBILITY ASSESSMENT / hEDS EVALUATION",
        {
            "chief_complaint": "Follow-up multiple conditions. Patient requests formal hypermobility assessment.",
            "hpi": (
                "Ms. Rodriguez returns requesting formal evaluation for Ehlers-Danlos syndrome. She has been "
                "reading about the POTS/hEDS/MCAS triad and notes she has all three conditions. She has "
                "lifelong joint hypermobility, chronic pain, frequent subluxations (shoulders, fingers), "
                "easy bruising, soft/velvety skin, and striae. Her mother also has 'loose joints' and GERD.\n\n"
                "Current diagnoses: POTS (confirmed tilt table), suspected MCAS (elevated mediators), "
                "small fiber neuropathy (confirmed IENFD). POTS partially controlled with fludrocortisone."
            ),
            "exam": (
                "BEIGHTON SCORE ASSESSMENT:\n"
                "1. Passive dorsiflexion of 5th MCP >90 degrees: LEFT yes, RIGHT yes (2 points)\n"
                "2. Passive apposition of thumb to forearm: LEFT yes, RIGHT yes (2 points)\n"
                "3. Hyperextension of elbow >10 degrees: LEFT yes, RIGHT yes (2 points)\n"
                "4. Hyperextension of knee >10 degrees: LEFT no, RIGHT yes (1 point)\n"
                "5. Forward flexion, palms flat on floor with knees straight: yes (1 point)\n"
                "TOTAL BEIGHTON SCORE: 8/9\n\n"
                "Skin: Soft, velvety texture. Mild skin hyperextensibility (extends >1.5 cm at forearm). "
                "Atrophic scars on knees. Striae on thighs and hips (non-pregnancy related). Easy bruising.\n\n"
                "2017 International hEDS Criteria Assessment:\n"
                "Criterion 1 (Generalized joint hypermobility): Beighton >= 5 -- MET\n"
                "Criterion 2 (must meet 2 of 3):\n"
                "  Feature A (5+ of 12 systemic features): skin hyperextensibility, atrophic scarring, "
                "striae, bilateral piezogenic papules on heels, recurrent subluxations, dental crowding, "
                "arachnodactyly -- 7/12 -- MET\n"
                "  Feature B (positive family history): mother with hypermobility -- MET\n"
                "  Feature C (MSK complications): chronic widespread pain >3 months, recurrent subluxations -- MET\n"
                "  2 of 3 features met: YES\n"
                "Criterion 3 (exclusions): No Marfan (normal echo), no other heritable CTD -- MET\n\n"
                "ALL THREE CRITERIA MET FOR hEDS DIAGNOSIS"
            ),
            "assessment": (
                "1. Hypermobile Ehlers-Danlos syndrome (hEDS) -- meets 2017 International Criteria\n"
                "2. POTS -- likely hEDS-associated (venous pooling from connective tissue laxity)\n"
                "3. Suspected MCAS -- elevated mediators, clinical episodes\n"
                "4. Small fiber neuropathy -- confirmed\n"
                "5. The POTS/hEDS/MCAS triad is now established in this patient"
            ),
            "plan": (
                "1. Formal diagnosis: hypermobile Ehlers-Danlos syndrome (hEDS)\n"
                "2. Referral to Genetics for comprehensive evaluation and genetic testing\n"
                "3. Physical therapy referral for joint stabilization program\n"
                "4. Continue current medications for POTS and MCAS\n"
                "5. Discuss genetic testing: consider TPSAB1 copy number (hereditary alpha-tryptasemia "
                "associated with elevated baseline tryptase and MCAS symptoms)\n"
                "6. Follow-up in 3 months"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 9: Genetics -- January 2024
    # ─────────────────────────────────────────────────────────────
    generate_genetic_report(
        os.path.join(output_dir, "2024-01-15_Genetics_Panel.pdf"),
        PATIENT, "01/15/2024",
        "HEREDITARY CONNECTIVE TISSUE DISORDER AND TRYPTASEMIA PANEL",
        "Hypermobile EDS, POTS, MCAS with elevated baseline tryptase. Evaluate for hereditary "
        "alpha-tryptasemia (HaT) and vascular EDS exclusion.",
        (
            "Targeted next-generation sequencing panel for connective tissue disorders (COL3A1, COL5A1, "
            "COL5A2, COL1A1, COL1A2, TNXB, FLNA, PLOD1, ADAMTS2, B4GALT7) plus TPSAB1 copy number "
            "analysis by digital droplet PCR."
        ),
        [
            {"Gene": "COL3A1", "Result": "No pathogenic variants", "Classification": "Normal"},
            {"Gene": "COL5A1", "Result": "No pathogenic variants", "Classification": "Normal"},
            {"Gene": "COL5A2", "Result": "No pathogenic variants", "Classification": "Normal"},
            {"Gene": "TNXB", "Result": "Heterozygous VUS (c.8092G>A)", "Classification": "VUS"},
            {"Gene": "TPSAB1", "Result": "3 copies (duplication)", "Classification": "Positive"},
        ],
        (
            "1. TPSAB1 DUPLICATION DETECTED (3 copies, normal = 2 copies)\n"
            "This is diagnostic of Hereditary Alpha-Tryptasemia (HaT). HaT is present in 4-6% of the "
            "general population and is strongly associated with elevated baseline tryptase levels, mast "
            "cell activation symptoms, and the POTS/hEDS/MCAS triad. This finding explains the patient's "
            "persistently elevated baseline tryptase (14.2 ng/mL) and recurrent mast cell activation "
            "episodes.\n\n"
            "2. No pathogenic variants identified in vascular EDS genes (COL3A1) -- vascular EDS excluded.\n\n"
            "3. TNXB heterozygous VUS (c.8092G>A, p.Val2698Met) -- variant of uncertain significance in "
            "tenascin-X gene. Heterozygous TNXB variants have been described in some hEDS-like phenotypes "
            "but clinical significance of this specific variant is uncertain."
        ),
        recommendations=(
            "1. Hereditary alpha-tryptasemia confirmed -- supports MCAS diagnosis and management\n"
            "2. Genetic counseling for autosomal dominant HaT (50% inheritance risk)\n"
            "3. Family screening recommended: mother should be tested for HaT given hypermobility history\n"
            "4. Vascular EDS excluded -- no arterial/organ rupture surveillance needed\n"
            "5. TNXB VUS: no clinical action at this time; reclassification may occur with future evidence"
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 10: Allergy follow-up -- May 2024 (MCAS confirmed)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-05-13_Allergy_Follow_Up.pdf"),
        PATIENT, "05/13/2024", "allergy",
        "ALLERGY/IMMUNOLOGY FOLLOW-UP -- MCAS MANAGEMENT",
        {
            "chief_complaint": "Follow-up MCAS management. Hereditary alpha-tryptasemia confirmed.",
            "hpi": (
                "Ms. Rodriguez returns for follow-up. Genetics confirmed TPSAB1 duplication (hereditary "
                "alpha-tryptasemia). She reports mast cell activation episodes have improved significantly "
                "with dual antihistamines and cromolyn sodium -- frequency reduced from monthly to "
                "approximately every 2-3 months, and severity reduced. She has identified some triggers: "
                "heat, alcohol, NSAIDs, and strong fragrances. She avoids these when possible.\n\n"
                "POTS remains partially controlled. hEDS management with PT ongoing."
            ),
            "assessment": (
                "1. Mast cell activation syndrome (MCAS) confirmed: recurrent episodic symptoms, "
                "elevated tryptase (both baseline and acute), elevated 24-hour urine mediators "
                "(N-methylhistamine, PGD2, LTE4), response to mast cell-directed therapy, AND "
                "hereditary alpha-tryptasemia (TPSAB1 duplication) as underlying genetic basis.\n"
                "2. This meets consensus diagnostic criteria for MCAS (Valent et al., 2019).\n"
                "3. Continue current successful regimen."
            ),
            "plan": (
                "1. Continue: cetirizine 10 mg BID, famotidine 20 mg BID, cromolyn sodium 200 mg QID\n"
                "2. Add: montelukast 10 mg daily (leukotriene receptor antagonist)\n"
                "3. Maintain trigger avoidance (heat, alcohol, NSAIDs, fragrances)\n"
                "4. Epinephrine auto-injector: ensure not expired, review technique\n"
                "5. Follow-up in 6 months\n"
                "6. Consider low-dose aspirin 81 mg for PGD2 blockade if symptoms recur despite current regimen"
            ),
        },
    )

    # Final summary note from PCP -- January 2025
    generate_progress_note(
        os.path.join(output_dir, "2025-01-13_PCP_Annual.pdf"),
        PATIENT, "01/13/2025", "pcp",
        "OFFICE VISIT -- ANNUAL WELLNESS / MULTI-CONDITION MANAGEMENT",
        {
            "chief_complaint": "Annual wellness exam. Multi-condition management review.",
            "hpi": (
                "Ms. Rodriguez is now 27 years old with established diagnoses of: (1) POTS -- confirmed "
                "by tilt table March 2022; (2) Hypermobile Ehlers-Danlos syndrome (hEDS) -- diagnosed "
                "July 2023 per 2017 criteria; (3) Mast cell activation syndrome (MCAS) -- confirmed with "
                "elevated mediators and hereditary alpha-tryptasemia (TPSAB1 duplication); (4) Small "
                "fiber neuropathy -- confirmed by IENFD biopsy April 2023.\n\n"
                "Overall she reports significant improvement with multi-modal management. POTS symptoms "
                "are 60-70% improved with fludrocortisone, salt/fluids, compression, and exercise. MCAS "
                "episodes have decreased to every 3-4 months with antihistamines, cromolyn, and montelukast. "
                "Physical therapy has helped joint stability. She completed her graduate degree and is "
                "working part-time.\n\n"
                "Time from initial symptoms to complete diagnosis: approximately 4 years."
            ),
            "vitals": (
                "Sitting: BP 116/72  HR 70\n"
                "Standing (3 min): BP 112/74  HR 94  (HR increase 24 bpm -- improved from 36+ bpm)\n"
                "Temp 98.4F  Wt 122 lbs"
            ),
            "medications": (
                "1. Fludrocortisone 0.1 mg daily\n"
                "2. Cetirizine 10 mg BID\n"
                "3. Famotidine 20 mg BID\n"
                "4. Cromolyn sodium 200 mg QID before meals\n"
                "5. Montelukast 10 mg daily\n"
                "6. Omeprazole 20 mg daily"
            ),
            "assessment": (
                "27F with POTS/hEDS/MCAS/SFN tetrad. Significant improvement with multi-modal therapy. "
                "Diagnostic odyssey of 4 years from initial symptoms (March 2021) to complete diagnosis "
                "(January 2024 with genetic confirmation). This case illustrates the typical multi-year, "
                "multi-specialist journey that patients with this symptom complex endure."
            ),
            "plan": (
                "1. Continue all current medications\n"
                "2. Annual labs: CBC, CMP, tryptase, vitamin D\n"
                "3. Continue PT for joint stabilization\n"
                "4. Recumbent exercise program ongoing\n"
                "5. Follow-up routine in 6 months\n"
                "6. Specialty follow-ups: Cardiology (annual), Allergy (every 6 months), Neurology (annual)"
            ),
        },
    )

    print(f"  Patient Maya Rodriguez: {len(os.listdir(output_dir))} documents generated")
