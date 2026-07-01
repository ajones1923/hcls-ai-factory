"""Additional specialist documents and subtle diagnostic labs.

Contains:
- Rachel Thompson: Genetics/HLA report (missing from original)
- Maya Rodriguez: Thyroid panel (normal -- rules out thyroid, misses POTS)
- Maya Rodriguez: Holter monitor (shows positional tachycardia but read as "sinus")
- Linda Chen: ANA from 2021 that was "weakly positive" and dismissed
- Sarah Mitchell: Early CBC showing subtle lymphopenia trend
- David Park: CRP trending data (standalone lab showing chronic inflammation)

Author: Adam Jones
Date: March 2026
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_genetic_report,
    generate_lab_report,
    generate_progress_note,
)

RACHEL = {
    "name": "Rachel Thompson",
    "dob": "1987-09-28",
    "mrn": "RTH-2020-78234",
    "sex": "F",
    "age_at_start": 32,
}

MAYA = {
    "name": "Maya Rodriguez",
    "dob": "1996-08-14",
    "mrn": "MRD-2021-45678",
    "sex": "F",
    "age_at_start": 24,
}

LINDA = {
    "name": "Linda Chen",
    "dob": "1973-05-19",
    "mrn": "LCH-2022-67890",
    "sex": "F",
    "age_at_start": 47,
}

SARAH = {
    "name": "Sarah Mitchell",
    "dob": "1991-03-22",
    "mrn": "SMI-2022-12345",
    "sex": "F",
    "age_at_start": 31,
}

DAVID = {
    "name": "David Park",
    "dob": "1980-11-03",
    "mrn": "DPA-2019-33102",
    "sex": "M",
    "age_at_start": 38,
}


def generate_rachel_genetics(output_dir: str):
    """Rachel Thompson -- HLA typing and autoimmune genetic panel."""
    os.makedirs(output_dir, exist_ok=True)

    generate_genetic_report(
        os.path.join(output_dir, "2021-01-18_Genetics_HLA.pdf"),
        RACHEL, "01/18/2021",
        "HLA HIGH-RESOLUTION TYPING AND AUTOIMMUNE RISK PANEL",
        (
            "Mixed connective tissue disease (MCTD) with strongly positive anti-U1-RNP. "
            "HLA typing for prognostic and pharmacogenomic purposes."
        ),
        (
            "Next-generation sequencing, high-resolution (2-field). Autoimmune susceptibility "
            "SNP panel for selected immune-regulatory genes."
        ),
        [
            {"Locus": "HLA-A", "Allele 1": "A*02:01", "Allele 2": "A*24:02"},
            {"Locus": "HLA-B", "Allele 1": "B*08:01", "Allele 2": "B*35:01"},
            {"Locus": "HLA-C", "Allele 1": "C*07:01", "Allele 2": "C*04:01"},
            {"Locus": "HLA-DRB1", "Allele 1": "DRB1*04:01", "Allele 2": "DRB1*03:01"},
            {"Locus": "HLA-DQB1", "Allele 1": "DQB1*03:01", "Allele 2": "DQB1*02:01"},
        ],
        (
            "HLA-DRB1*04:01 DETECTED -- This allele carries the 'shared epitope' associated with:\n"
            "- Rheumatoid arthritis susceptibility (OR 3.5-5.0)\n"
            "- Anti-U1-RNP-associated MCTD (established association)\n"
            "- More aggressive inflammatory joint disease\n"
            "- Higher risk of extra-articular manifestations including ILD\n\n"
            "HLA-DRB1*03:01 DETECTED -- This allele is associated with:\n"
            "- SLE susceptibility (OR 2.0-3.0)\n"
            "- Anti-Ro/SSA positivity (patient currently negative)\n"
            "- Sjogren's syndrome risk\n"
            "- Type 1 diabetes (note: no current evidence)\n\n"
            "HLA-B*08:01 DETECTED -- Part of the autoimmune-associated '8.1 ancestral haplotype' "
            "(A1-B8-DR3), which confers increased susceptibility to multiple autoimmune conditions "
            "including SLE, Sjogren's, myasthenia gravis, and celiac disease.\n\n"
            "AUTOIMMUNE SNP PANEL:\n"
            "- STAT4 rs7574865: T/T (risk allele homozygous) -- associated with SLE and RA "
            "susceptibility (OR 1.5-2.0)\n"
            "- IRF5 rs2004640: T/G (heterozygous) -- associated with SLE susceptibility\n"
            "- PTPN22 rs2476601: C/T (heterozygous) -- associated with RA and other autoimmune "
            "diseases (OR 1.5-2.0)\n\n"
            "SUMMARY: Ms. Thompson carries multiple genetic risk factors for autoimmune "
            "connective tissue disease. The combination of HLA-DRB1*04:01, DRB1*03:01, "
            "B*08:01, and risk alleles in STAT4/IRF5/PTPN22 places her at substantially "
            "elevated genetic risk. Her current MCTD diagnosis with anti-U1-RNP positivity "
            "is consistent with this genetic background. The HLA-DRB1*03:01 allele suggests "
            "she should be monitored for evolution toward SLE or Sjogren's features over time."
        ),
        recommendations=(
            "1. Monitor for evolution of MCTD toward a more defined CTD, particularly SLE "
            "(given DRB1*03:01 and STAT4 risk alleles)\n"
            "2. Annual anti-SSA/SSB testing given DRB1*03:01-associated Sjogren's risk\n"
            "3. The shared epitope (DRB1*04:01) suggests more aggressive joint disease -- "
            "early DMARD therapy is appropriate\n"
            "4. Genetic counseling offered -- patient interested in family planning discussion "
            "regarding autoimmune risk in offspring"
        ),
    )
    print("    + Rachel genetics report: 1 document added")


def generate_maya_additional(output_dir: str):
    """Maya Rodriguez -- thyroid panel and Holter monitor results."""
    os.makedirs(output_dir, exist_ok=True)

    # ── Thyroid Panel (2021-08-02) -- Normal, rules out thyroid ───
    # PCP orders thyroid before sending to psychiatry. Normal. Reinforces
    # the "it must be anxiety" narrative.
    generate_lab_report(
        os.path.join(output_dir, "2021-08-02_Lab_Thyroid.pdf"),
        MAYA, "08/02/2021", PROVIDERS["pcp"],
        [{
            "panel_name": "THYROID FUNCTION",
            "results": [
                {"test": "TSH", "value": "2.1", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Free T4", "value": "1.3", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": ""},
                {"test": "Free T3", "value": "3.2", "unit": "pg/mL", "ref_range": "2.3-4.2", "flag": ""},
            ],
        }, {
            "panel_name": "METABOLIC / ANEMIA SCREENING",
            "results": [
                {"test": "Glucose, Fasting", "value": "82", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "Hemoglobin", "value": "13.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Iron, Serum", "value": "68", "unit": "mcg/dL", "ref_range": "60-170", "flag": ""},
                {"test": "Ferritin", "value": "32", "unit": "ng/mL", "ref_range": "12-150", "flag": ""},
                {"test": "Vitamin B12", "value": "412", "unit": "pg/mL", "ref_range": "200-900", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "24", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
            ],
        }],
    )

    # ── Holter Monitor (2021-12-20) -- 24-hr recording ────────────
    # Ordered after 2nd ER visit. Shows positional tachycardia but reported
    # as "normal sinus rhythm" because the standard Holter interpretation
    # doesn't flag positional changes.
    generate_progress_note(
        os.path.join(output_dir, "2021-12-20_Holter_Results.pdf"),
        MAYA, "12/20/2021", "cardiology",
        "HOLTER MONITOR -- 24-HOUR AMBULATORY ECG RESULTS",
        {
            "chief_complaint": "24-hour Holter monitor ordered after ER visit for tachycardia. Rule out arrhythmia.",
            "labs_reviewed": (
                "24-HOUR HOLTER MONITOR -- Recorded 12/16/2021 - 12/17/2021\n\n"
                "HEART RATE SUMMARY:\n"
                "Minimum HR: 52 bpm (03:42 AM, sleeping)\n"
                "Maximum HR: 148 bpm (11:14 AM, patient diary: 'standing in kitchen making lunch')\n"
                "Mean HR: 84 bpm\n"
                "Mean daytime HR: 92 bpm\n"
                "Mean nighttime HR: 62 bpm\n\n"
                "RHYTHM:\n"
                "- Predominant rhythm: Normal sinus rhythm throughout\n"
                "- No supraventricular tachycardia (SVT)\n"
                "- No ventricular tachycardia (VT)\n"
                "- No atrial fibrillation or flutter\n"
                "- Isolated PACs: 42 (within normal limits)\n"
                "- Isolated PVCs: 18 (within normal limits)\n"
                "- No pauses >2.0 seconds\n\n"
                "PATIENT DIARY CORRELATION:\n"
                "- 08:15 AM: 'Got out of bed, felt dizzy' -- HR 62 -> 128 bpm over 2 minutes\n"
                "- 11:14 AM: 'Standing cooking' -- HR 138-148 bpm sustained for 22 minutes\n"
                "- 12:30 PM: 'Sat down to eat' -- HR 148 -> 78 bpm within 4 minutes\n"
                "- 02:45 PM: 'Walking in Target' -- HR 124 bpm, diary: 'felt heart pounding, "
                "had to sit down'\n"
                "- 05:30 PM: 'Hot shower' -- HR 142 bpm\n"
                "- 09:00 PM: 'Lying on couch watching TV' -- HR 68 bpm, diary: 'feel fine'\n\n"
                "ST SEGMENT ANALYSIS:\n"
                "No significant ST depression or elevation."
            ),
            "assessment": (
                "1. NORMAL 24-HOUR HOLTER MONITOR -- no arrhythmia detected\n"
                "2. Sinus tachycardia episodes noted, correlated with upright activity per "
                "patient diary. Maximum HR 148 bpm during routine standing activity (making "
                "lunch) is NOTABLE.\n"
                "3. Significant heart rate variability with position changes: resting/supine "
                "62-78 bpm, standing 124-148 bpm (delta 50-86 bpm)\n\n"
                "INTERPRETATION NOTE: While no primary arrhythmia was identified, the heart "
                "rate pattern shows marked positional dependence. Heart rate rises > 30 bpm "
                "with standing on multiple occasions during this recording. This pattern is "
                "consistent with postural orthostatic tachycardia syndrome (POTS) and warrants "
                "formal tilt table evaluation."
            ),
            "plan": (
                "1. No anti-arrhythmic therapy needed\n"
                "2. Recommend formal tilt table testing to evaluate for POTS\n"
                "3. Referring to Dr. Thompson (Cardiology) for follow-up\n"
                "4. Continue metoprolol ER 25 mg daily for symptomatic relief"
            ),
        },
    )

    print("    + Maya additional reports: 2 documents added")


def generate_linda_additional(output_dir: str):
    """Linda Chen -- early dismissed ANA test."""
    os.makedirs(output_dir, exist_ok=True)

    # ── ANA from June 2022 that PCP ordered but downplayed ───────
    # PCP finally orders ANA after ophth and dentist both urged it, but
    # the initial screen is "weakly positive" (1:80) and dismissed.
    # The FULL panel isn't ordered until September 2022.
    generate_lab_report(
        os.path.join(output_dir, "2022-06-27_Lab_ANA_Screen.pdf"),
        LINDA, "06/27/2022", PROVIDERS["pcp"],
        [{
            "panel_name": "AUTOIMMUNE SCREENING",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:80", "unit": "", "ref_range": "<1:40", "flag": "A"},
                {"test": "ANA Pattern", "value": "Speckled", "unit": "", "ref_range": "N/A", "flag": ""},
            ],
        }],
    )

    # PCP follow-up where 1:80 ANA is dismissed
    generate_progress_note(
        os.path.join(output_dir, "2022-07-11_PCP_ANA_Review.pdf"),
        LINDA, "07/11/2022", "pcp",
        "TELEPHONE ENCOUNTER -- LAB RESULTS",
        {
            "chief_complaint": "Lab result review. ANA test result.",
            "hpi": (
                "Telephone call with Ms. Chen to review ANA results ordered 06/27/2022 after "
                "letters from both Ophthalmology (Dr. Wells, 09/2021) and Dentistry (Dr. Rivera, "
                "07/2022) recommending Sjogren's evaluation.\n\n"
                "ANA result: Positive at 1:80, speckled pattern."
            ),
            "assessment": (
                "ANA positive at 1:80 (low titer). A titer of 1:80 can be found in up to 15% of "
                "healthy individuals and increases with age. Given that the patient is 49 years old "
                "and perimenopausal, a low-titer ANA is of uncertain clinical significance.\n\n"
                "However, both Ophthalmology and Dentistry have raised concern for Sjogren's "
                "syndrome. The anti-SSA/Ro antibody is more specific for Sjogren's than ANA alone. "
                "A negative ANA at low titer does NOT rule out Sjogren's -- up to 30% of Sjogren's "
                "patients are ANA-negative, and anti-SSA can be positive with a negative or low-titer "
                "ANA.\n\n"
                "Decision: Will monitor clinically. If symptoms worsen, will order anti-SSA/SSB panel."
            ),
            "plan": (
                "1. Discussed with patient: ANA is weakly positive. This is likely not significant "
                "but will continue to monitor.\n"
                "2. If dry eyes/mouth worsen or new symptoms develop (joint pain, rash, fatigue), "
                "return for anti-SSA/SSB testing\n"
                "3. Continue current medications\n"
                "4. Routine follow-up in 3 months\n\n"
                "NOTE: Patient asked 'could this be Sjogren's like my eye doctor mentioned?' "
                "Advised that a 1:80 ANA is very common and not diagnostic. Reassured that her "
                "symptoms are most likely perimenopausal."
            ),
        },
    )

    print("    + Linda additional reports: 2 documents added")


def generate_sarah_additional(output_dir: str):
    """Sarah Mitchell -- early subtle labs showing lupus signatures."""
    os.makedirs(output_dir, exist_ok=True)

    # ── Pre-clinical CBC (2022-03-07) -- Annual physical labs ─────
    # 3 months BEFORE she presents with symptoms. Already showing subtle
    # lymphopenia that a human wouldn't flag but the AI system can track.
    generate_lab_report(
        os.path.join(output_dir, "2022-03-07_Lab_Annual_Physical.pdf"),
        SARAH, "03/07/2022", PROVIDERS["pcp"],
        [{
            "panel_name": "CBC WITH DIFFERENTIAL",
            "results": [
                {"test": "WBC", "value": "4.4", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "12.8", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Hematocrit", "value": "38.4", "unit": "%", "ref_range": "36.0-46.0", "flag": ""},
                {"test": "Platelets", "value": "168", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Neutrophils, Abs", "value": "2.8", "unit": "K/uL", "ref_range": "1.8-7.7", "flag": ""},
                {"test": "Lymphocytes, Abs", "value": "1.1", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": ""},
                {"test": "Monocytes, Abs", "value": "0.4", "unit": "K/uL", "ref_range": "0.2-0.8", "flag": ""},
            ],
        }, {
            "panel_name": "COMPREHENSIVE METABOLIC PANEL",
            "results": [
                {"test": "Glucose, Fasting", "value": "88", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "BUN", "value": "12", "unit": "mg/dL", "ref_range": "7-20", "flag": ""},
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
                {"test": "Sodium", "value": "140", "unit": "mEq/L", "ref_range": "136-145", "flag": ""},
                {"test": "Potassium", "value": "4.2", "unit": "mEq/L", "ref_range": "3.5-5.1", "flag": ""},
                {"test": "Albumin", "value": "4.0", "unit": "g/dL", "ref_range": "3.5-5.5", "flag": ""},
                {"test": "ALT", "value": "18", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "22", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
            ],
        }, {
            "panel_name": "LIPID PANEL",
            "results": [
                {"test": "Total Cholesterol", "value": "186", "unit": "mg/dL", "ref_range": "<200", "flag": ""},
                {"test": "LDL, Calculated", "value": "108", "unit": "mg/dL", "ref_range": "<130", "flag": ""},
                {"test": "HDL", "value": "58", "unit": "mg/dL", "ref_range": ">40", "flag": ""},
                {"test": "Triglycerides", "value": "100", "unit": "mg/dL", "ref_range": "<150", "flag": ""},
            ],
        }],
    )

    # ── Annual Physical Progress Note (2022-03-07) ────────────────
    # Completely normal visit. But in hindsight, the lymphocyte count
    # of 1.1 (low-normal) and low vitamin D are early lupus signals.
    generate_progress_note(
        os.path.join(output_dir, "2022-03-07_PCP_Annual_Physical.pdf"),
        SARAH, "03/07/2022", "pcp",
        "ANNUAL PHYSICAL EXAM",
        {
            "chief_complaint": "Annual physical exam. No acute complaints.",
            "hpi": (
                "Ms. Mitchell is a 30-year-old woman presenting for her annual physical exam. "
                "She reports generally feeling well. She mentions occasional fatigue that she "
                "attributes to work stress (she recently started a new position as a marketing "
                "manager). She exercises 3-4 times per week. Diet is balanced. No significant "
                "medical history. No prior surgeries. No medications."
            ),
            "vitals": "BP 112/68  HR 72  Temp 98.4F  Wt 138 lbs  Ht 5'6\"  BMI 22.3",
            "ros": (
                "Constitutional: Mild fatigue (attributes to stress). HEENT: Occasional headaches. "
                "CV: No chest pain. Resp: No SOB. GI: No issues. MSK: No joint pain. "
                "Skin: Notes she gets 'red cheeks' after sun exposure but this has always been "
                "the case. Neuro: Occasional brain fog. Psych: Mild work stress."
            ),
            "exam": (
                "General: Well-appearing, well-nourished female.\n"
                "HEENT: Normal. No oral ulcers. No lymphadenopathy.\n"
                "CV: RRR, no murmur. Lungs: CTA bilateral.\n"
                "Abdomen: Soft, non-tender. Skin: Clear, no rash today.\n"
                "MSK: Full ROM all joints. No tenderness."
            ),
            "assessment": (
                "1. Healthy 30-year-old woman. Annual exam unremarkable.\n"
                "2. Mild fatigue -- likely stress-related. Sleep hygiene discussed.\n"
                "3. Low vitamin D (22 ng/mL) -- supplement recommended.\n"
                "4. Labs otherwise within normal limits."
            ),
            "plan": (
                "1. Vitamin D3 2000 IU daily\n"
                "2. Continue exercise, healthy diet\n"
                "3. Age-appropriate cancer screening: up to date\n"
                "4. Return in 1 year or sooner if concerns"
            ),
        },
    )

    print("    + Sarah additional reports: 2 documents added")


def generate_david_additional(output_dir: str):
    """David Park -- CRP trending lab showing chronic inflammation over years."""
    os.makedirs(output_dir, exist_ok=True)

    # ── Annual Physical Labs (2021-01-11) -- CRP still high, unaddressed ──
    generate_lab_report(
        os.path.join(output_dir, "2021-01-11_Lab_Annual.pdf"),
        DAVID, "01/11/2021", PROVIDERS["pcp"],
        [{
            "panel_name": "ANNUAL LABS",
            "results": [
                {"test": "CRP", "value": "20.6", "unit": "mg/L", "ref_range": "<3.0", "flag": "HH"},
                {"test": "Glucose, Fasting", "value": "96", "unit": "mg/dL", "ref_range": "70-100", "flag": ""},
                {"test": "Hemoglobin", "value": "14.0", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "Creatinine", "value": "1.0", "unit": "mg/dL", "ref_range": "0.7-1.3", "flag": ""},
                {"test": "Total Cholesterol", "value": "204", "unit": "mg/dL", "ref_range": "<200", "flag": "H"},
                {"test": "LDL, Calculated", "value": "128", "unit": "mg/dL", "ref_range": "<130", "flag": ""},
                {"test": "HDL", "value": "44", "unit": "mg/dL", "ref_range": ">40", "flag": ""},
                {"test": "Triglycerides", "value": "160", "unit": "mg/dL", "ref_range": "<150", "flag": "H"},
                {"test": "ALT", "value": "32", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "HbA1c", "value": "5.4", "unit": "%", "ref_range": "<5.7", "flag": ""},
            ],
        }],
    )

    print("    + David additional reports: 1 document added")
