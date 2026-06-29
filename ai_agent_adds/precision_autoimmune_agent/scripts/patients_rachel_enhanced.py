"""Enhanced documents for Rachel Thompson -- expand MCTD follow-up through 2023.

Adds ~9 documents:
- Rheumatology follow-ups with steroid taper monitoring (2021)
- Echocardiogram for PAH screening (2021)
- PFTs for ILD monitoring (2021)
- Ongoing rheumatology management (2022, 2023)
- Serial labs showing disease trajectory (2021, 2022, 2023)
- Follow-up HRCT chest (2022)
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
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


def generate_rachel_enhanced(output_dir: str):
    """Generate additional documents for Rachel Thompson's MCTD journey."""
    os.makedirs(output_dir, exist_ok=True)

    # ── 2021-05-17  Rheum Follow-up -- Steroid Taper ──────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-05-17_Rheum_Follow_Up.pdf"),
        RACHEL, "05/17/2021", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- MCTD STEROID TAPER",
        {
            "chief_complaint": "Follow-up MCTD. Steroid taper in progress. Myositis improving.",
            "hpi": (
                "Ms. Thompson returns 3 months after last labs. She has tapered prednisone from "
                "40 mg to 10 mg daily over the past 3 months per protocol. Myositis symptoms much "
                "improved -- able to climb stairs without difficulty, lifting overhead improved. "
                "Raynaud's unchanged (2-3 episodes/week). Puffy fingers slightly improved. "
                "Hydroxychloroquine 200 mg BID ongoing. She reports new symptom: dysphagia with "
                "solid foods, occasional food getting 'stuck.' No weight loss."
            ),
            "vitals": "BP 120/74  HR 72  Temp 98.6F  Wt 142 lbs",
            "exam": (
                "Motor: Proximal strength improved -- hip flexors 4+/5 bilateral, deltoids 5-/5 "
                "bilateral. No atrophy.\n"
                "Hands: Puffy fingers mildly improved. No sclerodactyly. Mild Raynaud's changes.\n"
                "Joints: MCP tenderness minimal. No active synovitis.\n"
                "Lungs: Clear bilateral."
            ),
            "assessment": (
                "1. MCTD -- myositis component improving on steroid taper. CK trending down.\n"
                "2. New dysphagia -- esophageal dysmotility is common in MCTD (70% of patients). "
                "The prior HRCT noted mildly dilated esophagus. Need barium swallow or manometry.\n"
                "3. Raynaud's -- persistent but manageable\n"
                "4. Early NSIP on HRCT -- will need PFTs for monitoring"
            ),
            "plan": (
                "1. Continue prednisone taper: 10 mg x 2 weeks, then 7.5 mg, then 5 mg maintenance\n"
                "2. Continue hydroxychloroquine 200 mg BID\n"
                "3. If myositis flares on taper below 10 mg, add methotrexate as steroid-sparing agent\n"
                "4. Barium swallow for dysphagia evaluation\n"
                "5. Labs: CK, aldolase, CRP, CBC, CMP\n"
                "6. PFTs ordered (baseline for ILD monitoring)\n"
                "7. Echocardiogram ordered (PAH screening -- leading cause of death in MCTD)\n"
                "8. Follow-up in 3 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2021-05-17_Lab_Follow_Up.pdf"),
        RACHEL, "05/17/2021", PROVIDERS["rheumatology"],
        [{
            "panel_name": "MCTD MONITORING",
            "results": [
                {"test": "CK, Total", "value": "142", "unit": "U/L", "ref_range": "30-135", "flag": "H"},
                {"test": "Aldolase", "value": "7.8", "unit": "U/L", "ref_range": "1.0-7.5", "flag": "H"},
                {"test": "CRP", "value": "2.8", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "18", "unit": "mm/hr", "ref_range": "0-20", "flag": ""},
                {"test": "Anti-U1-RNP", "value": "Positive (5.8)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC / CMP",
            "results": [
                {"test": "WBC", "value": "8.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.2", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "278", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Glucose, Fasting", "value": "108", "unit": "mg/dL", "ref_range": "70-100", "flag": "H"},
                {"test": "Creatinine", "value": "0.7", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
            ],
        }],
    )

    # ── 2021-07-12  Echocardiogram -- PAH Screening ───────────────
    generate_imaging_report(
        os.path.join(output_dir, "2021-07-12_Echo_PAH_Screening.pdf"),
        RACHEL, "07/12/2021", "cardiology",
        "Echocardiogram", "Transthoracic",
        "MCTD with ILD screening. Evaluate for pulmonary arterial hypertension.",
        "Complete 2D and Doppler transthoracic echocardiogram with color flow mapping.",
        (
            "Left ventricle: Normal size (LVEDD 4.6 cm). Normal systolic function (EF 60-65%). "
            "No wall motion abnormalities. Normal diastolic function.\n"
            "Right ventricle: Normal size and function. TAPSE 2.2 cm (normal >1.7 cm).\n"
            "Atria: Normal size bilaterally.\n"
            "Valves: Trivial tricuspid regurgitation.\n"
            "RVSP: Estimated 28 mmHg (normal <35 mmHg). Unable to fully estimate due to "
            "trivial TR jet -- may be underestimated.\n"
            "Pericardium: No effusion.\n"
            "IVC: Normal size (1.8 cm), >50% collapse with inspiration."
        ),
        (
            "1. Normal biventricular size and function.\n"
            "2. RVSP estimated at 28 mmHg -- NORMAL, but note: TR jet was trivial and RVSP "
            "may be underestimated. If clinical suspicion for PAH is high, consider right "
            "heart catheterization.\n"
            "3. No pericardial effusion.\n\n"
            "RECOMMENDATION: Annual echocardiographic screening recommended for MCTD patients "
            "given 10-45% lifetime risk of PAH."
        ),
        reading_radiologist="Michael Thompson, MD, FACC",
    )

    # ── 2021-08-09  PFTs ──────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-08-09_PFT_Results.pdf"),
        RACHEL, "08/09/2021", "pcp",
        "PULMONARY FUNCTION TEST RESULTS",
        {
            "chief_complaint": "Baseline PFTs for ILD monitoring in MCTD.",
            "hpi": (
                "Ms. Thompson is a 33-year-old woman with MCTD (anti-U1-RNP positive) with early "
                "NSIP on HRCT (November 2020). Baseline PFTs ordered by Rheumatology for monitoring."
            ),
            "labs_reviewed": (
                "SPIROMETRY:\n"
                "FVC: 3.42 L (96% predicted) -- Normal\n"
                "FEV1: 2.88 L (98% predicted) -- Normal\n"
                "FEV1/FVC: 84% -- Normal\n\n"
                "LUNG VOLUMES:\n"
                "TLC: 5.10 L (93% predicted) -- Normal\n"
                "RV: 1.68 L (98% predicted) -- Normal\n\n"
                "DIFFUSING CAPACITY:\n"
                "DLCO: 21.8 mL/min/mmHg (82% predicted) -- LOW-NORMAL\n"
                "DLCO/VA (KCO): 4.4 mL/min/mmHg/L (88% predicted) -- Normal"
            ),
            "assessment": (
                "1. Normal spirometry and lung volumes\n"
                "2. DLCO at low end of normal (82% predicted) -- in context of MCTD with "
                "early NSIP on HRCT, this requires close monitoring. DLCO decline is often "
                "the earliest functional change in ILD.\n"
                "3. Baseline values established for longitudinal comparison"
            ),
            "plan": (
                "1. Results forwarded to Rheumatology\n"
                "2. Repeat PFTs in 6-12 months\n"
                "3. If DLCO declines >10% absolute or FVC declines >10% predicted, "
                "consider treatment escalation for ILD"
            ),
        },
    )

    # ── 2021-11-15  Rheum Follow-up ───────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-11-15_Rheum_Follow_Up.pdf"),
        RACHEL, "11/15/2021", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- MCTD MAINTENANCE",
        {
            "chief_complaint": "6-month follow-up. Steroid taper complete. Mild flare of myositis.",
            "hpi": (
                "Ms. Thompson is now 1 year post-MCTD diagnosis. She tapered prednisone to 5 mg "
                "but noted return of proximal muscle weakness over the past 3 weeks. Difficulty "
                "with overhead activities again. Also reports intermittent pleuritic chest pain "
                "(sharp, worse with deep breathing) x 2 weeks. Raynaud's stable. Dysphagia "
                "evaluated -- barium swallow showed mild esophageal dysmotility. Echo was normal "
                "(RVSP 28 mmHg). PFTs showed low-normal DLCO (82%)."
            ),
            "vitals": "BP 118/72  HR 76  Temp 98.8F  Wt 140 lbs",
            "exam": (
                "Motor: Proximal weakness returned -- deltoids 4/5 bilateral, hip flexors 4/5.\n"
                "Chest: No pleural rub. Clear to auscultation.\n"
                "Hands: Puffy fingers persist. Mild Raynaud's changes."
            ),
            "assessment": (
                "1. MCTD -- myositis FLARE on low-dose steroids. Need steroid-sparing agent.\n"
                "2. Pleuritic chest pain -- serositis (pleuritis) is common in MCTD. "
                "May need chest imaging if persists.\n"
                "3. Esophageal dysmotility -- confirmed, mild.\n"
                "4. DLCO 82% -- borderline, will monitor.\n"
                "5. PAH screening normal."
            ),
            "plan": (
                "1. Increase prednisone to 15 mg daily temporarily\n"
                "2. START methotrexate 15 mg weekly (steroid-sparing for myositis)\n"
                "3. Folic acid 1 mg daily (MTX supplement)\n"
                "4. Pre-MTX labs: CBC, CMP, hepatitis B/C, chest X-ray\n"
                "5. If pleuritic pain persists >2 weeks, chest CT to rule out effusion/pericarditis\n"
                "6. PPI for esophageal dysmotility symptoms\n"
                "7. Labs in 4 weeks (MTX monitoring), follow-up in 3 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2021-11-15_Lab_Flare.pdf"),
        RACHEL, "11/15/2021", PROVIDERS["rheumatology"],
        [{
            "panel_name": "MCTD FLARE ASSESSMENT",
            "results": [
                {"test": "CK, Total", "value": "524", "unit": "U/L", "ref_range": "30-135", "flag": "HH"},
                {"test": "Aldolase", "value": "10.2", "unit": "U/L", "ref_range": "1.0-7.5", "flag": "H"},
                {"test": "CRP", "value": "8.4", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "42", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "Anti-U1-RNP", "value": "Positive (7.1)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
                {"test": "C3", "value": "78", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "C4", "value": "14", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
            ],
        }],
    )

    # ── 2022-03-14  Rheum Follow-up -- MTX Response ───────────────
    generate_progress_note(
        os.path.join(output_dir, "2022-03-14_Rheum_Follow_Up.pdf"),
        RACHEL, "03/14/2022", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- MCTD ON METHOTREXATE",
        {
            "chief_complaint": "4-month follow-up on methotrexate. Myositis improving.",
            "hpi": (
                "Ms. Thompson is 4 months into methotrexate 15 mg weekly. Prednisone tapered back "
                "to 7.5 mg daily. Proximal weakness significantly improved -- back to normal activities. "
                "Pleuritic chest pain resolved after 3 weeks (was self-limited serositis). Raynaud's "
                "same. Dysphagia mild, managed with PPI and dietary modification. No MTX side effects "
                "(no nausea, no mouth sores, no hair loss). Folic acid 1 mg daily."
            ),
            "vitals": "BP 116/72  HR 70  Temp 98.4F  Wt 138 lbs",
            "exam": (
                "Motor: Proximal strength 5/5 throughout -- normalized.\n"
                "Hands: Puffy fingers improving. No sclerodactyly.\n"
                "Joints: Non-tender. No synovitis.\n"
                "Lungs: Clear.\n"
                "Skin: No rash."
            ),
            "assessment": (
                "1. MCTD -- EXCELLENT RESPONSE to methotrexate. Myositis resolved. CK normalizing.\n"
                "2. Steroid taper progressing well -- target 5 mg then off.\n"
                "3. Serositis resolved.\n"
                "4. Continue monitoring for ILD progression and PAH development."
            ),
            "plan": (
                "1. Continue methotrexate 15 mg weekly\n"
                "2. Continue prednisone taper: 7.5 mg x 4 weeks -> 5 mg x 4 weeks -> consider d/c\n"
                "3. Continue hydroxychloroquine 200 mg BID\n"
                "4. Labs: CK, CBC, CMP, LFTs (MTX monitoring)\n"
                "5. Repeat PFTs in 3 months (due at 12-month mark from baseline)\n"
                "6. Annual echocardiogram (next July 2022)\n"
                "7. Follow-up in 4 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2022-03-14_Lab_MTX_Monitoring.pdf"),
        RACHEL, "03/14/2022", PROVIDERS["rheumatology"],
        [{
            "panel_name": "MCTD / MTX MONITORING",
            "results": [
                {"test": "CK, Total", "value": "118", "unit": "U/L", "ref_range": "30-135", "flag": ""},
                {"test": "Aldolase", "value": "6.4", "unit": "U/L", "ref_range": "1.0-7.5", "flag": ""},
                {"test": "CRP", "value": "2.1", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "16", "unit": "mm/hr", "ref_range": "0-20", "flag": ""},
                {"test": "Anti-U1-RNP", "value": "Positive (4.8)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
            ],
        }, {
            "panel_name": "MTX SAFETY MONITORING",
            "results": [
                {"test": "WBC", "value": "6.4", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.6", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "298", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "ALT", "value": "32", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "28", "unit": "U/L", "ref_range": "10-40", "flag": ""},
                {"test": "Creatinine", "value": "0.7", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
                {"test": "Albumin", "value": "4.2", "unit": "g/dL", "ref_range": "3.5-5.5", "flag": ""},
            ],
        }],
    )

    # ── 2022-09-19  Follow-up HRCT Chest ──────────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2022-09-19_HRCT_Chest_Follow_Up.pdf"),
        RACHEL, "09/19/2022", "rheumatology",
        "CT", "Chest HRCT",
        "MCTD with early NSIP (November 2020). Follow-up for ILD progression.",
        "High-resolution CT of the chest without contrast. Comparison: November 2020.",
        (
            "Lungs: Persistent minimal ground-glass opacity in bilateral lower lobes posteriorly, "
            "UNCHANGED compared to November 2020. No new areas of ground-glass opacity. "
            "No honeycombing. No traction bronchiectasis. No consolidation.\n"
            "Airways: Normal. No air trapping.\n"
            "Pleura: No effusion or thickening.\n"
            "Mediastinum: No lymphadenopathy.\n"
            "Heart: Normal size.\n"
            "Esophagus: Mildly dilated, unchanged. Consistent with known dysmotility."
        ),
        (
            "1. STABLE minimal bilateral lower lobe ground-glass opacities compared to November "
            "2020 (22 months prior). No progression of NSIP pattern.\n"
            "2. No new ILD findings.\n"
            "3. Stable esophageal dilatation.\n\n"
            "IMPRESSION: Stable early NSIP with no evidence of progression on current immunosuppressive "
            "therapy (methotrexate + hydroxychloroquine)."
        ),
    )

    # ── 2023-03-13  Rheum Annual -- Disease Status ────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-03-13_Rheum_Annual.pdf"),
        RACHEL, "03/13/2023", "rheumatology",
        "RHEUMATOLOGY ANNUAL REVIEW -- MCTD",
        {
            "chief_complaint": "Annual comprehensive MCTD review. Overall stable.",
            "hpi": (
                "Ms. Thompson is now 2.5 years post-MCTD diagnosis. Current medications: methotrexate "
                "15 mg weekly, hydroxychloroquine 200 mg BID, nifedipine 30 mg ER daily, folic acid "
                "1 mg daily. Prednisone discontinued 6 months ago -- no flare. Myositis in complete "
                "remission (CK normal x 9 months). Raynaud's reduced to 1-2 episodes/week on "
                "nifedipine. Puffy fingers resolved. Dysphagia improved on PPI. No pleuritic pain "
                "recurrence. Follow-up HRCT (September 2022) showed stable early NSIP.\n\n"
                "New concern: She is planning pregnancy in the next 1-2 years and wants to discuss "
                "medication safety."
            ),
            "vitals": "BP 114/68  HR 68  Temp 98.2F  Wt 136 lbs",
            "exam": (
                "General: Well-appearing.\n"
                "Motor: Full proximal and distal strength 5/5.\n"
                "Hands: No puffy fingers. No sclerodactyly. Nailfold capillaroscopy: mild "
                "architectural changes (dilated capillaries, no avascular areas).\n"
                "Joints: Non-tender. No synovitis.\n"
                "Lungs: Clear bilateral.\n"
                "Skin: No rash. No calcinosis."
            ),
            "assessment": (
                "MCTD -- IN REMISSION on methotrexate + hydroxychloroquine.\n\n"
                "Disease Summary at 2.5 years:\n"
                "- Myositis: REMISSION (CK normalized, strength 5/5)\n"
                "- Arthritis: REMISSION (no synovitis)\n"
                "- Raynaud's: CONTROLLED on nifedipine\n"
                "- ILD (NSIP): STABLE, no progression\n"
                "- Esophageal dysmotility: STABLE, mild\n"
                "- Serositis: RESOLVED, no recurrence\n"
                "- PAH: No evidence (RVSP 28, normal)\n\n"
                "Anti-U1-RNP remains positive but declining titer.\n\n"
                "Re: pregnancy planning -- MTX is TERATOGENIC (Category X). Must be discontinued "
                "3 months prior to conception. HCQ is safe in pregnancy and should be continued."
            ),
            "plan": (
                "1. Continue current regimen for now\n"
                "2. PREGNANCY PLANNING:\n"
                "   - When ready: d/c MTX 3 months before conception\n"
                "   - Continue HCQ throughout pregnancy (protective against neonatal lupus)\n"
                "   - Add azathioprine if steroid-sparing agent needed during pregnancy\n"
                "   - Anti-SSA/SSB testing (for neonatal lupus risk) -- currently negative, recheck\n"
                "   - High-risk OB referral when pregnant\n"
                "3. Labs: CK, CBC, CMP, LFTs, anti-SSA/SSB, complement C3/C4\n"
                "4. Annual echocardiogram (PAH screening)\n"
                "5. PFTs in 6 months\n"
                "6. Follow-up in 6 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2023-03-13_Lab_Annual.pdf"),
        RACHEL, "03/13/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "MCTD ANNUAL MONITORING",
            "results": [
                {"test": "CK, Total", "value": "98", "unit": "U/L", "ref_range": "30-135", "flag": ""},
                {"test": "CRP", "value": "1.4", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "12", "unit": "mm/hr", "ref_range": "0-20", "flag": ""},
                {"test": "Anti-U1-RNP", "value": "Positive (3.6)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
                {"test": "Anti-SSA/Ro", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-SSB/La", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "C3", "value": "102", "unit": "mg/dL", "ref_range": "90-180", "flag": ""},
                {"test": "C4", "value": "22", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
            ],
        }, {
            "panel_name": "MTX SAFETY + CBC",
            "results": [
                {"test": "WBC", "value": "5.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.8", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "312", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "ALT", "value": "26", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "22", "unit": "U/L", "ref_range": "10-40", "flag": ""},
                {"test": "Creatinine", "value": "0.7", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient Rachel Thompson: {count} documents total (enhanced)")
