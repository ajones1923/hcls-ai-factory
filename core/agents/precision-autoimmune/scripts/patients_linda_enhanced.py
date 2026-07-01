"""Enhanced documents for Linda Chen -- expand Sjogren's to true 3-year odyssey.

Adds ~10 documents:
- Earlier ophthalmology visits showing progressive dry eye (2020, 2021)
- Earlier dentistry visit (2021)
- Lab monitoring (2023, 2024)
- Rheumatology follow-ups (2023, 2024)
- Pulmonary function test screening (2023)
- Chest X-ray ILD screening (2023)
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
    generate_lab_report,
    generate_progress_note,
)

LINDA = {
    "name": "Linda Chen",
    "dob": "1973-05-19",
    "mrn": "LCH-2022-67890",
    "sex": "F",
    "age_at_start": 47,
}


def generate_linda_enhanced(output_dir: str):
    """Generate additional documents for Linda Chen's Sjogren's odyssey."""
    os.makedirs(output_dir, exist_ok=True)

    # ── 2020-06-15  Ophthalmology -- First Dry Eye Visit ──────────
    generate_progress_note(
        os.path.join(output_dir, "2020-06-15_Ophth_Dry_Eyes_Initial.pdf"),
        LINDA, "06/15/2020", "ophthalmology",
        "OPHTHALMOLOGY VISIT -- DRY EYE EVALUATION",
        {
            "chief_complaint": "Eyes feel dry and irritated x 3 months. Worse with reading and computer work.",
            "hpi": (
                "Ms. Chen is a 47-year-old schoolteacher presenting with bilateral dry eye symptoms. "
                "She describes intermittent grittiness, burning, and foreign body sensation x 3 months, "
                "worse toward end of day and with screen time. She spends 6+ hours daily on computer "
                "for lesson planning. Denies any other systemic symptoms. Tried OTC Systane drops "
                "with modest relief."
            ),
            "exam": (
                "OD: VA 20/20. Schirmer test (without anesthesia): 8 mm/5 min (borderline low).\n"
                "OS: VA 20/20. Schirmer test: 9 mm/5 min (borderline low).\n"
                "Slit lamp: Mild punctate epithelial erosions inferior cornea bilateral on "
                "fluorescein staining. Tear breakup time 6 seconds (reduced, normal >10s). "
                "Meibomian glands: mild inspissation."
            ),
            "assessment": (
                "1. Mild-moderate dry eye disease, bilateral\n"
                "2. Evaporative component (meibomian gland dysfunction) and possible "
                "aqueous-deficient component (borderline Schirmer)"
            ),
            "plan": (
                "1. Preservative-free artificial tears QID\n"
                "2. Warm compresses and lid hygiene BID\n"
                "3. Omega-3 supplements 2g daily\n"
                "4. 20-20-20 rule for screen time\n"
                "5. Return in 6 months -- if worsening, consider Restasis\n"
                "6. Schirmer scores are borderline -- will monitor trend"
            ),
        },
    )

    # ── 2021-01-11  Ophthalmology Follow-up ───────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-01-11_Ophth_Follow_Up.pdf"),
        LINDA, "01/11/2021", "ophthalmology",
        "OPHTHALMOLOGY FOLLOW-UP -- DRY EYES",
        {
            "chief_complaint": "Dry eyes worsening despite OTC drops and lid hygiene.",
            "hpi": (
                "Ms. Chen returns 6 months later with progressive dry eye symptoms. Artificial tears "
                "help briefly but she needs them every 1-2 hours. Morning eye crusting. She also "
                "reports new symptom: her mouth has been feeling dry, though she attributes this to "
                "increased caffeine intake. No joint pain, no rash."
            ),
            "exam": (
                "OD: VA 20/25. Schirmer test: 6 mm/5 min (decreased from 8 mm).\n"
                "OS: VA 20/25. Schirmer test: 7 mm/5 min (decreased from 9 mm).\n"
                "Slit lamp: Moderate punctate epithelial erosions bilateral. Tear breakup time "
                "4 seconds (worsened). Reduced tear meniscus height."
            ),
            "assessment": (
                "1. Progressive aqueous-deficient dry eye disease, bilateral. Schirmer trending "
                "down (8/9 -> 6/7 in 6 months).\n"
                "2. Mention of dry mouth -- if persistent, could suggest systemic process "
                "(e.g., Sjogren's syndrome). Currently no indication for workup."
            ),
            "plan": (
                "1. Start Restasis (cyclosporine 0.05%) BID both eyes\n"
                "2. Continue preservative-free artificial tears q2h\n"
                "3. Omega-3 supplements continue\n"
                "4. Follow-up in 6 months\n"
                "5. If dry mouth persists, recommend discussing with PCP"
            ),
        },
    )

    # ── 2021-08-23  Ophthalmology Follow-up -- worsening ──────────
    generate_progress_note(
        os.path.join(output_dir, "2021-08-23_Ophth_Follow_Up.pdf"),
        LINDA, "08/23/2021", "ophthalmology",
        "OPHTHALMOLOGY FOLLOW-UP -- PROGRESSIVE DRY EYES",
        {
            "chief_complaint": "Dry eyes not controlled with Restasis. Needs drops every hour.",
            "hpi": (
                "Ms. Chen returns after 7 months on Restasis with limited improvement. She is "
                "using artificial tears almost hourly. Eyes constantly gritty. She now also reports "
                "dry mouth is persistent and worsening -- she carries a water bottle everywhere and "
                "wakes at night to drink. She has had 2 dental cavities in the past year (unusual "
                "for her). No joint symptoms yet."
            ),
            "exam": (
                "OD: VA 20/25. Schirmer test: 4 mm/5 min (worsened from 6 mm).\n"
                "OS: VA 20/25. Schirmer test: 5 mm/5 min (worsened from 7 mm).\n"
                "Slit lamp: Moderate-severe punctate erosions. Filamentary keratitis early changes. "
                "TBUT 3 seconds. Tear meniscus barely visible."
            ),
            "assessment": (
                "1. SEVERE aqueous-deficient dry eye disease, bilateral. Schirmer now <5 mm "
                "(consistent with severe keratoconjunctivitis sicca).\n"
                "2. Progressive decline over 14 months: Schirmer 8-9 -> 6-7 -> 4-5 mm.\n"
                "3. Combined sicca symptoms (dry eyes + dry mouth) with accelerated dental caries "
                "raise concern for SJOGREN'S SYNDROME. This is no longer simple dry eye.\n"
                "4. Early filamentary keratitis -- may need more aggressive treatment."
            ),
            "plan": (
                "1. Add punctal plugs, lower lids bilateral (office procedure today)\n"
                "2. Continue Restasis\n"
                "3. Serum eye drops if not improving in 6 weeks\n"
                "4. RECOMMEND: ANA, anti-SSA/Ro, anti-SSB/La through PCP\n"
                "5. RECOMMEND: Discuss with PCP re: Sjogren's syndrome evaluation\n"
                "6. Follow-up in 3 months\n\n"
                "NOTE: Letter sent to PCP Dr. Martinez recommending Sjogren's workup."
            ),
        },
    )

    # ── 2021-11-01  Dentistry -- Earlier Visit ────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2021-11-01_Dentistry_Caries.pdf"),
        LINDA, "11/01/2021", "dentistry",
        "DENTAL VISIT -- ACCELERATED CARIES",
        {
            "chief_complaint": "Three new cavities found on routine exam. Dry mouth worsening.",
            "hpi": (
                "Ms. Chen presents for semi-annual dental exam. She reports persistent dry mouth "
                "x 9+ months. Difficulty swallowing crackers or bread without water. She chews "
                "sugar-free gum frequently. She has excellent oral hygiene habits (brushes 2x daily, "
                "flosses daily) yet is developing cavities at an unprecedented rate -- 3 new carious "
                "lesions found today. She had no cavities in the 10 years prior to this."
            ),
            "exam": (
                "Oral exam: Dry oral mucosa. Scant saliva. Tongue slightly dry. "
                "Caries at #5 (cervical), #12 (root surface), #29 (cervical). "
                "All at cervical/root surface -- classic xerostomia pattern. "
                "Mild bilateral parotid fullness on palpation."
            ),
            "assessment": (
                "1. Accelerated cervical/root surface caries secondary to xerostomia\n"
                "2. Xerostomia of unclear etiology -- not medication-related (patient on no "
                "anticholinergic medications)\n"
                "3. Bilateral parotid fullness\n"
                "4. Pattern concerning for possible Sjogren's syndrome"
            ),
            "plan": (
                "1. Restore carious teeth (#5, #12, #29)\n"
                "2. Prescription fluoride toothpaste (PreviDent 5000)\n"
                "3. Biotene mouth rinse\n"
                "4. Saliva stimulant: sugar-free xylitol lozenges\n"
                "5. Increase dental visit frequency to every 3 months\n"
                "6. Recommend PCP evaluation for dry mouth cause -- possible autoimmune etiology"
            ),
        },
    )

    # ── 2023-02-20  Rheumatology Follow-up (3 months post-dx) ────
    generate_progress_note(
        os.path.join(output_dir, "2023-02-20_Rheum_Follow_Up.pdf"),
        LINDA, "02/20/2023", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- SJOGREN'S SYNDROME",
        {
            "chief_complaint": "3-month follow-up on hydroxychloroquine and pilocarpine.",
            "hpi": (
                "Ms. Chen returns 3 months after starting hydroxychloroquine 200 mg BID and "
                "pilocarpine 5 mg TID. She reports modest improvement in fatigue (60% of baseline) "
                "and slight increase in saliva production. Arthralgia improved but not resolved. "
                "Dry eyes stable with Restasis and punctal plugs. She has had no parotid swelling "
                "episodes since starting treatment. No new symptoms."
            ),
            "vitals": "BP 118/72  HR 68  Temp 98.4F  Wt 132 lbs",
            "exam": (
                "HEENT: Oral mucosa slightly improved from prior. No parotid enlargement today.\n"
                "Musculoskeletal: Minimal MCP tenderness. No synovitis.\n"
                "Skin: No purpura. No vasculitic changes.\n"
                "Lymph nodes: No lymphadenopathy."
            ),
            "assessment": (
                "Primary Sjogren's Syndrome -- PARTIAL RESPONSE to hydroxychloroquine + pilocarpine.\n"
                "- Fatigue: improved (ESSPRI fatigue domain 6 -> 4)\n"
                "- Dryness: mildly improved (ESSPRI dryness domain 7 -> 6)\n"
                "- Pain: improved (ESSPRI pain domain 5 -> 3)\n"
                "- ESSPRI total: 6.0 -> 4.3 (clinically meaningful improvement >1 point)"
            ),
            "plan": (
                "1. Continue hydroxychloroquine 200 mg BID\n"
                "2. Continue pilocarpine 5 mg TID\n"
                "3. Ophthalmology annual retinal exam (HCQ monitoring)\n"
                "4. Labs: CBC, CMP, CRP, ESR, immunoglobulins, complement C3/C4\n"
                "5. Baseline chest X-ray for ILD screening (ordered but not yet completed)\n"
                "6. Discuss lymphoma surveillance -- patient counseled on warning signs "
                "(persistent lymphadenopathy, B symptoms, parotid mass)\n"
                "7. Follow-up in 4 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2023-02-20_Lab_Monitoring.pdf"),
        LINDA, "02/20/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "SJOGREN'S MONITORING",
            "results": [
                {"test": "CRP", "value": "3.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "28", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "IgG, Total", "value": "1842", "unit": "mg/dL", "ref_range": "700-1600", "flag": "H"},
                {"test": "IgA, Total", "value": "380", "unit": "mg/dL", "ref_range": "70-400", "flag": ""},
                {"test": "IgM, Total", "value": "245", "unit": "mg/dL", "ref_range": "40-230", "flag": "H"},
                {"test": "C3", "value": "98", "unit": "mg/dL", "ref_range": "90-180", "flag": ""},
                {"test": "C4", "value": "18", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "4.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "12.8", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "188", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "Lymphocytes, Abs", "value": "0.9", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": "L"},
            ],
        }, {
            "panel_name": "RENAL",
            "results": [
                {"test": "Creatinine", "value": "0.8", "unit": "mg/dL", "ref_range": "0.6-1.1", "flag": ""},
                {"test": "BUN", "value": "14", "unit": "mg/dL", "ref_range": "7-20", "flag": ""},
                {"test": "CO2 (Bicarb)", "value": "23", "unit": "mEq/L", "ref_range": "22-29", "flag": ""},
                {"test": "Potassium", "value": "4.1", "unit": "mEq/L", "ref_range": "3.5-5.1", "flag": ""},
            ],
        }],
    )

    # ── 2023-04-03  Chest X-ray -- ILD Screening ──────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2023-04-03_CXR_ILD_Screening.pdf"),
        LINDA, "04/03/2023", "rheumatology",
        "X-Ray", "Chest PA and Lateral",
        "Primary Sjogren's syndrome. Screening for interstitial lung disease.",
        "PA and lateral chest radiographs obtained.",
        (
            "Heart size normal. Mediastinal contours normal. No lymphadenopathy.\n"
            "Lungs: Clear bilaterally. No focal consolidation, mass, or nodule. "
            "No interstitial markings suggesting fibrosis. No pleural effusion.\n"
            "Osseous structures unremarkable."
        ),
        (
            "1. Normal chest radiograph.\n"
            "2. No evidence of interstitial lung disease.\n"
            "3. Note: Chest X-ray has limited sensitivity for early ILD. If clinical suspicion "
            "increases (new cough, dyspnea on exertion), HRCT recommended."
        ),
    )

    # ── 2023-06-26  PFT -- Pulmonary Function Testing ─────────────
    # Using progress note format since we don't have a dedicated PFT generator
    generate_progress_note(
        os.path.join(output_dir, "2023-06-26_PFT_Results.pdf"),
        LINDA, "06/26/2023", "pcp",
        "PULMONARY FUNCTION TEST RESULTS",
        {
            "chief_complaint": "Routine PFT screening for interstitial lung disease in Sjogren's syndrome.",
            "hpi": (
                "Ms. Chen is a 50-year-old woman with Primary Sjogren's Syndrome (anti-SSA/SSB+, "
                "focus score 3) undergoing baseline PFT screening. ILD occurs in 9-20% of pSS "
                "patients. She denies cough, dyspnea, or exercise intolerance."
            ),
            "labs_reviewed": (
                "SPIROMETRY:\n"
                "FVC: 3.21 L (92% predicted) -- Normal\n"
                "FEV1: 2.68 L (94% predicted) -- Normal\n"
                "FEV1/FVC: 83% -- Normal\n\n"
                "LUNG VOLUMES:\n"
                "TLC: 4.88 L (95% predicted) -- Normal\n"
                "RV: 1.67 L (102% predicted) -- Normal\n\n"
                "DIFFUSING CAPACITY:\n"
                "DLCO: 19.2 mL/min/mmHg (78% predicted) -- MILDLY REDUCED\n"
                "DLCO/VA (KCO): 4.1 mL/min/mmHg/L (84% predicted) -- Low-normal"
            ),
            "assessment": (
                "1. Normal spirometry and lung volumes\n"
                "2. MILDLY REDUCED DLCO (78% predicted) -- may represent early subclinical "
                "interstitial involvement or small airway disease in Sjogren's. This is a "
                "common early finding even with normal chest X-ray.\n"
                "3. No symptoms currently"
            ),
            "plan": (
                "1. Results forwarded to Rheumatology (Dr. Chen)\n"
                "2. Repeat PFTs in 12 months to monitor trend\n"
                "3. If DLCO declines further (<70% predicted) or symptoms develop, "
                "proceed to HRCT chest\n"
                "4. No treatment change needed at this time"
            ),
        },
    )

    # ── 2023-11-06  Rheumatology Follow-up ────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-11-06_Rheum_Follow_Up.pdf"),
        LINDA, "11/06/2023", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- SJOGREN'S MAINTENANCE",
        {
            "chief_complaint": "Routine follow-up. Overall stable.",
            "hpi": (
                "Ms. Chen returns for routine Sjogren's follow-up. She is 1 year post-diagnosis. "
                "Hydroxychloroquine and pilocarpine continue. Fatigue is manageable. Arthralgia mild. "
                "Dry eyes stable on Restasis + punctal plugs. Dry mouth improved with pilocarpine. "
                "No dental cavities at last check (3-month interval visits). PFTs showed mildly "
                "reduced DLCO (78%) -- asymptomatic from pulmonary standpoint.\n\n"
                "She reports one new symptom: intermittent purpuric lesions on lower legs x 2 months, "
                "non-pruritic, resolve spontaneously."
            ),
            "vitals": "BP 116/70  HR 64  Temp 98.2F  Wt 130 lbs",
            "exam": (
                "HEENT: Oral mucosa moist (improved on pilocarpine). No parotid enlargement.\n"
                "Skin: 3-4 small non-palpable purpuric macules bilateral shins. Non-blanching.\n"
                "Musculoskeletal: No synovitis. Minimal MCP tenderness.\n"
                "Lymph nodes: No lymphadenopathy."
            ),
            "assessment": (
                "1. Primary Sjogren's Syndrome -- overall stable on current regimen\n"
                "2. Non-palpable purpura lower extremities -- hypergammaglobulinemic purpura "
                "is common in Sjogren's (polyclonal B-cell activation). Monitor for evolution "
                "to palpable purpura (vasculitis) or cryoglobulinemia.\n"
                "3. Mildly reduced DLCO -- monitor annually\n"
                "4. Lymphopenia (0.9 K/uL) -- common in Sjogren's, monitor"
            ),
            "plan": (
                "1. Continue hydroxychloroquine 200 mg BID, pilocarpine 5 mg TID\n"
                "2. Labs: CBC, CMP, CRP, immunoglobulins, cryoglobulins, free light chains, "
                "beta-2 microglobulin\n"
                "3. If cryoglobulins positive, add hepatitis C screening\n"
                "4. Annual ophthalmology retinal screening (HCQ toxicity monitoring)\n"
                "5. Repeat PFTs in 6 months\n"
                "6. Follow-up in 6 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2023-11-06_Lab_Surveillance.pdf"),
        LINDA, "11/06/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "SJOGREN'S SURVEILLANCE",
            "results": [
                {"test": "IgG, Total", "value": "2104", "unit": "mg/dL", "ref_range": "700-1600", "flag": "H"},
                {"test": "Kappa Free Light Chain", "value": "32.4", "unit": "mg/L", "ref_range": "3.3-19.4", "flag": "H"},
                {"test": "Lambda Free Light Chain", "value": "18.8", "unit": "mg/L", "ref_range": "5.7-26.3", "flag": ""},
                {"test": "Kappa/Lambda Ratio", "value": "1.72", "unit": "", "ref_range": "0.26-1.65", "flag": "H"},
                {"test": "Beta-2 Microglobulin", "value": "2.8", "unit": "mg/L", "ref_range": "0.7-1.8", "flag": "H"},
                {"test": "Cryoglobulins", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "CRP", "value": "3.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "30", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "3.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": "L"},
                {"test": "Lymphocytes, Abs", "value": "0.8", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": "L"},
                {"test": "Hemoglobin", "value": "12.4", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "172", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # ── 2024-05-13  Annual Labs + Lymphoma Surveillance ───────────
    generate_lab_report(
        os.path.join(output_dir, "2024-05-13_Lab_Annual.pdf"),
        LINDA, "05/13/2024", PROVIDERS["rheumatology"],
        [{
            "panel_name": "SJOGREN'S ANNUAL MONITORING",
            "results": [
                {"test": "CRP", "value": "2.4", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "26", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "IgG, Total", "value": "2280", "unit": "mg/dL", "ref_range": "700-1600", "flag": "H"},
                {"test": "IgA, Total", "value": "402", "unit": "mg/dL", "ref_range": "70-400", "flag": "H"},
                {"test": "C3", "value": "88", "unit": "mg/dL", "ref_range": "90-180", "flag": "L"},
                {"test": "C4", "value": "12", "unit": "mg/dL", "ref_range": "10-40", "flag": ""},
                {"test": "RF", "value": "62", "unit": "IU/mL", "ref_range": "<14", "flag": "H"},
                {"test": "Kappa/Lambda Ratio", "value": "1.84", "unit": "", "ref_range": "0.26-1.65", "flag": "H"},
                {"test": "Beta-2 Microglobulin", "value": "3.1", "unit": "mg/L", "ref_range": "0.7-1.8", "flag": "H"},
                {"test": "LDH", "value": "198", "unit": "U/L", "ref_range": "120-246", "flag": ""},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "3.6", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": "L"},
                {"test": "Lymphocytes, Abs", "value": "0.7", "unit": "K/uL", "ref_range": "1.0-4.8", "flag": "L"},
                {"test": "Hemoglobin", "value": "12.2", "unit": "g/dL", "ref_range": "12.0-16.0", "flag": ""},
                {"test": "Platelets", "value": "164", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient Linda Chen: {count} documents total (enhanced)")
