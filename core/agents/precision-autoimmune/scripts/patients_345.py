"""Patients 3-5: David (AS), Linda (Sjogren's), Rachel (MCTD).

Shorter timelines focused on key diagnostic patterns.
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

# =====================================================================
# Patient 3: David Park -- 45M, 6-year AS diagnostic odyssey
# =====================================================================

DAVID = {
    "name": "David Park",
    "dob": "1980-11-03",
    "mrn": "DPA-2019-33102",
    "sex": "M",
    "age_at_start": 38,
}


def generate_david(output_dir: str):
    """Generate documents for David Park -- ankylosing spondylitis."""
    os.makedirs(output_dir, exist_ok=True)

    generate_progress_note(
        os.path.join(output_dir, "2019-02-11_PCP_Back_Pain.pdf"),
        DAVID, "02/11/2019", "pcp",
        "OFFICE VISIT -- LOW BACK PAIN",
        {
            "chief_complaint": "Low back pain x 6 months, worse in the morning.",
            "hpi": (
                "Mr. Park is a 38-year-old construction supervisor with chronic low back pain. Pain "
                "started insidiously about 6 months ago without specific injury. He reports significant "
                "morning stiffness lasting 45-60 minutes that improves with activity. Pain is worse with "
                "rest and at night -- sometimes wakes him. Improves with movement and NSAIDs. He rates "
                "pain 6/10 at worst. No radiculopathy, no bowel/bladder changes. No history of psoriasis "
                "or IBD. He reports occasional achilles tendon pain bilaterally."
            ),
            "vitals": "BP 128/82  HR 72  Temp 98.6F  Wt 185 lbs  Ht 5'10\"",
            "exam": (
                "Lumbar spine: Reduced forward flexion (modified Schober: 3 cm increase, normal >5 cm). "
                "Tenderness over SI joints bilaterally. Negative straight leg raise. No neurologic deficits.\n"
                "Achilles: Mild tenderness at bilateral insertions (enthesitis).\n"
                "Chest expansion: 3.5 cm (borderline reduced, normal >5 cm)."
            ),
            "assessment": (
                "1. Chronic low back pain -- inflammatory features (morning stiffness >30 min, improvement "
                "with activity, night pain, insidious onset). However, given his occupation, mechanical "
                "causes are also possible.\n"
                "2. Bilateral achilles enthesitis"
            ),
            "plan": (
                "1. Lumbar spine X-ray\n"
                "2. Naproxen 500 mg BID with meals\n"
                "3. Physical therapy referral for core strengthening\n"
                "4. If not improving in 6 weeks, consider MRI and inflammatory markers\n"
                "5. Return in 6 weeks"
            ),
        },
    )

    generate_imaging_report(
        os.path.join(output_dir, "2019-02-11_XR_Lumbar.pdf"),
        DAVID, "02/11/2019", "pcp",
        "X-Ray", "Lumbar Spine AP and Lateral",
        "Chronic low back pain, evaluate for degenerative changes.",
        "AP and lateral views of the lumbar spine obtained.",
        (
            "Vertebral body heights are maintained. Mild disc space narrowing at L4-L5. "
            "No spondylolisthesis. No compression fractures. Sacroiliac joints are partially "
            "obscured by overlying bowel gas but appear grossly normal on this limited view. "
            "No significant osteophyte formation."
        ),
        (
            "1. Mild L4-L5 degenerative disc disease.\n"
            "2. No acute fracture or malalignment.\n"
            "3. SI joints not well visualized on this study -- dedicated SI joint views or MRI "
            "recommended if clinical concern for sacroiliitis."
        ),
    )

    # PCP attributes to mechanical -- treats with PT for 2 years
    generate_progress_note(
        os.path.join(output_dir, "2019-08-19_PCP_Follow_Up.pdf"),
        DAVID, "08/19/2019", "pcp",
        "OFFICE VISIT -- FOLLOW-UP BACK PAIN",
        {
            "chief_complaint": "Ongoing back pain despite PT. Morning stiffness persists.",
            "hpi": (
                "Mr. Park completed 12 weeks of PT with modest improvement. Morning stiffness "
                "persists (45 min daily). He notes the pain is now radiating into the buttocks "
                "bilaterally (alternating). NSAIDs help but he cannot take them continuously. "
                "He has had one episode of red, painful right eye (saw Ophthalmology -- diagnosed "
                "with anterior uveitis, treated with steroid drops, resolved)."
            ),
            "assessment": (
                "1. Persistent inflammatory back pain -- not responding well to mechanical treatments\n"
                "2. History of anterior uveitis -- noted\n"
                "3. Bilateral achilles enthesitis -- persistent"
            ),
            "plan": (
                "1. Continue NSAIDs PRN\n"
                "2. CRP, ESR to check for systemic inflammation\n"
                "3. If inflammatory markers elevated, will consider rheumatology referral\n"
                "4. Return in 3 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2019-08-19_Lab_Inflammatory.pdf"),
        DAVID, "08/19/2019", PROVIDERS["pcp"],
        [{
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP", "value": "18.4", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "34", "unit": "mm/hr", "ref_range": "0-15", "flag": "H"},
            ],
        }],
    )

    # Ophthalmology note for uveitis
    generate_progress_note(
        os.path.join(output_dir, "2019-07-02_Ophth_Uveitis.pdf"),
        DAVID, "07/02/2019", "ophthalmology",
        "OPHTHALMOLOGY VISIT -- ACUTE ANTERIOR UVEITIS",
        {
            "chief_complaint": "Red, painful right eye x 3 days. Photophobia.",
            "hpi": (
                "Mr. Park presents with acute onset right eye pain, redness, and photophobia x 3 days. "
                "No trauma. No prior episodes. No contact lens use."
            ),
            "exam": (
                "OD: VA 20/30 (baseline 20/20). Ciliary injection. 2+ anterior chamber cells. "
                "Fine KP on endothelium. No posterior synechiae. IOP 12.\n"
                "OS: Normal. VA 20/20. No cells."
            ),
            "assessment": "Acute anterior uveitis, right eye. Non-granulomatous.",
            "plan": (
                "1. Prednisolone acetate 1% q2h x 1 week, then taper\n"
                "2. Cyclopentolate 1% TID for 5 days\n"
                "3. If recurrent, recommend HLA-B27 testing and rheumatology referral\n"
                "4. Follow-up in 1 week"
            ),
        },
    )

    # Years pass -- finally referred to rheumatology in 2022
    generate_progress_note(
        os.path.join(output_dir, "2022-03-07_Orthopedics_Consult.pdf"),
        DAVID, "03/07/2022", "orthopedics",
        "ORTHOPEDIC CONSULTATION -- CHRONIC LOW BACK PAIN",
        {
            "chief_complaint": "Persistent low back pain x 3 years, failing conservative management.",
            "hpi": (
                "Mr. Park is a 41-year-old man with 3+ years of chronic inflammatory-type low back pain "
                "not adequately managed with PT and NSAIDs. Referred to evaluate for surgical candidacy. "
                "Morning stiffness 60+ minutes daily. CRP persistently elevated (18-24 mg/L). History of "
                "anterior uveitis. Bilateral alternating buttock pain. Achilles enthesitis."
            ),
            "exam": (
                "Severely reduced lumbar flexion. Modified Schober: 2 cm increase. Chest expansion 2.8 cm. "
                "Bilateral SI joint tenderness. FABER test positive bilaterally."
            ),
            "assessment": (
                "This patient's presentation is more consistent with inflammatory spondyloarthritis "
                "than mechanical low back pain. The combination of inflammatory back pain, reduced "
                "spinal mobility, enthesitis, anterior uveitis, and elevated CRP strongly suggests "
                "ankylosing spondylitis. I do NOT recommend surgical intervention.\n\n"
                "I recommend urgent rheumatology referral and MRI of sacroiliac joints."
            ),
            "plan": (
                "1. MRI sacroiliac joints with STIR sequences\n"
                "2. URGENT rheumatology referral\n"
                "3. HLA-B27 testing\n"
                "4. Continue NSAIDs\n"
                "5. No surgical intervention indicated"
            ),
        },
    )

    generate_imaging_report(
        os.path.join(output_dir, "2022-03-14_MRI_SI_Joints.pdf"),
        DAVID, "03/14/2022", "orthopedics",
        "MRI", "Sacroiliac Joints",
        "Suspected inflammatory spondyloarthritis. Chronic back pain with elevated CRP and uveitis history.",
        "MRI of the sacroiliac joints with T1, T2, STIR, and post-gadolinium sequences.",
        (
            "Bilateral sacroiliac joints demonstrate:\n"
            "- Bilateral subchondral bone marrow edema on STIR sequences (left > right), "
            "consistent with active sacroiliitis\n"
            "- Erosions at bilateral iliac surfaces, more prominent on the left\n"
            "- Subchondral sclerosis bilaterally\n"
            "- Partial ankylosis of the inferior portion of the left SI joint\n"
            "- Post-gadolinium enhancement along bilateral SI joint surfaces\n\n"
            "Lumbar spine: mild disc desiccation at L4-L5 and L5-S1. No significant "
            "spinal stenosis. Early squaring of L2 and L3 vertebral bodies noted."
        ),
        (
            "1. BILATERAL SACROILIITIS with active inflammation (bone marrow edema) and chronic "
            "structural changes (erosions, sclerosis, partial ankylosis) -- findings consistent "
            "with ankylosing spondylitis.\n"
            "2. Early vertebral squaring at L2-L3 suggesting early spinal involvement.\n"
            "3. Modified New York criteria: Grade III bilateral sacroiliitis."
        ),
    )

    generate_lab_report(
        os.path.join(output_dir, "2022-03-14_Lab_HLA_B27.pdf"),
        DAVID, "03/14/2022", PROVIDERS["orthopedics"],
        [{
            "panel_name": "HLA AND INFLAMMATORY",
            "results": [
                {"test": "HLA-B27", "value": "POSITIVE", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "CRP", "value": "24.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "HH"},
                {"test": "ESR", "value": "42", "unit": "mm/hr", "ref_range": "0-15", "flag": "H"},
            ],
        }],
    )

    generate_progress_note(
        os.path.join(output_dir, "2022-04-18_Rheum_AS_Diagnosis.pdf"),
        DAVID, "04/18/2022", "rheumatology",
        "RHEUMATOLOGY NEW PATIENT -- ANKYLOSING SPONDYLITIS DIAGNOSIS",
        {
            "chief_complaint": "Referred for suspected ankylosing spondylitis.",
            "hpi": (
                "Mr. Park is a 41-year-old man with 3+ years of inflammatory back pain, bilateral "
                "sacroiliitis on MRI, positive HLA-B27, history of anterior uveitis, bilateral achilles "
                "enthesitis, and persistently elevated CRP. Referred by Orthopedics."
            ),
            "assessment": (
                "DIAGNOSIS: Ankylosing Spondylitis\n\n"
                "Modified New York Criteria: MET\n"
                "- Clinical: inflammatory back pain >3 months, limited lumbar motion, reduced chest expansion\n"
                "- Radiographic: bilateral sacroiliitis Grade III\n"
                "- HLA-B27: positive\n\n"
                "ASDAS-CRP = 3.8 (very high disease activity; >3.5 = very high)\n"
                "BASDAI = 6.4 (active disease; >=4 = active)\n\n"
                "Time from symptom onset to diagnosis: approximately 6 years.\n"
                "This delay is unfortunately typical -- average AS diagnostic delay is 6-8 years."
            ),
            "plan": (
                "1. Start biologic therapy: adalimumab 40 mg SC every 2 weeks (TNF inhibitor)\n"
                "2. Pre-biologic screening: TB (QuantiFERON), Hepatitis B/C, CBC, CMP\n"
                "3. Continue naproxen 500 mg BID as bridge\n"
                "4. PT referral: AS-specific exercise program (emphasis on extension, posture, flexibility)\n"
                "5. Ophthalmology: annual screening for uveitis recurrence\n"
                "6. Follow-up in 12 weeks to assess biologic response\n"
                "7. Target: ASDAS <1.3 (inactive disease) or ASDAS improvement >=1.1"
            ),
        },
    )

    # HLA typing -- confirm B*27:05
    generate_genetic_report(
        os.path.join(output_dir, "2022-05-02_Genetics_HLA.pdf"),
        DAVID, "05/02/2022",
        "HLA HIGH-RESOLUTION TYPING",
        "Ankylosing spondylitis, positive HLA-B27 screening. High-resolution typing for subtype.",
        "Next-generation sequencing, high-resolution (2-field).",
        [
            {"Locus": "HLA-A", "Allele 1": "A*02:01", "Allele 2": "A*33:03"},
            {"Locus": "HLA-B", "Allele 1": "B*27:05", "Allele 2": "B*44:03"},
            {"Locus": "HLA-C", "Allele 1": "C*01:02", "Allele 2": "C*16:01"},
            {"Locus": "HLA-DRB1", "Allele 1": "DRB1*04:05", "Allele 2": "DRB1*13:02"},
        ],
        (
            "HLA-B*27:05 confirmed. This is the subtype most strongly associated with ankylosing "
            "spondylitis (OR 87.4, the strongest known HLA-disease association). The arthritogenic "
            "peptide hypothesis proposes that B27 presents self-peptides to autoreactive CD8+ T-cells; "
            "additionally, B27 misfolding triggers the unfolded protein response (UPR) and IL-23 "
            "production, driving Th17-mediated inflammation at entheseal sites.\n\n"
            "HLA-DRB1*04:05 is also noted -- this allele is associated with RA susceptibility (shared "
            "epitope), though no RA features are present. It may be relevant if inflammatory arthritis "
            "of peripheral joints develops."
        ),
    )

    print(f"  Patient David Park: {len(os.listdir(output_dir))} documents generated")


# =====================================================================
# Patient 4: Linda Chen -- 52F, 3-year Sjogren's odyssey
# =====================================================================

LINDA = {
    "name": "Linda Chen",
    "dob": "1973-05-19",
    "mrn": "LCH-2022-67890",
    "sex": "F",
    "age_at_start": 49,
}


def generate_linda(output_dir: str):
    """Generate documents for Linda Chen -- Sjogren's syndrome."""
    os.makedirs(output_dir, exist_ok=True)

    generate_progress_note(
        os.path.join(output_dir, "2022-04-04_Ophth_Dry_Eyes.pdf"),
        LINDA, "04/04/2022", "ophthalmology",
        "OPHTHALMOLOGY VISIT -- CHRONIC DRY EYES",
        {
            "chief_complaint": "Dry, gritty eyes x 1 year. OTC artificial tears not helping.",
            "hpi": (
                "Ms. Chen is a 49-year-old woman with progressive dry eye symptoms. She describes "
                "persistent gritty, burning sensation in both eyes. Worse with computer use and in "
                "air-conditioned environments. Has tried multiple OTC artificial tears with limited relief. "
                "She also notes her eyes tire easily when reading."
            ),
            "exam": (
                "OD: VA 20/25. Schirmer test (without anesthesia): 3 mm/5 min (abnormal, <5 mm).\n"
                "OS: VA 20/25. Schirmer test: 4 mm/5 min (abnormal).\n"
                "Slit lamp: Bilateral punctate epithelial erosions on fluorescein staining. "
                "Reduced tear meniscus. No anterior chamber inflammation."
            ),
            "assessment": (
                "1. Severe aqueous-deficient dry eye disease, bilateral\n"
                "2. Schirmer test markedly reduced bilaterally (<5 mm)"
            ),
            "plan": (
                "1. Restasis (cyclosporine 0.05%) BID both eyes\n"
                "2. Preservative-free artificial tears q2h\n"
                "3. Punctal plugs lower lids bilaterally\n"
                "4. Omega-3 supplement 2g daily\n"
                "5. If not improving, consider Sjogren's workup (ANA, anti-SSA/Ro)\n"
                "6. Follow-up in 3 months"
            ),
        },
    )

    generate_progress_note(
        os.path.join(output_dir, "2022-06-13_Dentistry_Dry_Mouth.pdf"),
        LINDA, "06/13/2022", "dentistry",
        "DENTAL VISIT -- INCREASED CARIES, DRY MOUTH",
        {
            "chief_complaint": "Multiple new dental caries. Persistent dry mouth.",
            "hpi": (
                "Ms. Chen presents for routine dental exam with 4 new carious lesions despite good oral "
                "hygiene. She reports persistent dry mouth (xerostomia) x 1 year, difficulty swallowing "
                "dry foods, need to sip water constantly. She uses sugar-free gum and lozenges."
            ),
            "exam": (
                "Oral exam: Dry oral mucosa. Scant saliva pooling. Tongue appears dry and fissured. "
                "Caries at #3, #14, #19, #30 -- all cervical/root surface caries. Mild bilateral "
                "parotid gland enlargement on palpation."
            ),
            "assessment": (
                "1. Rampant cervical caries secondary to xerostomia\n"
                "2. Bilateral parotid enlargement\n"
                "3. Pattern of dry mouth + dry eyes + parotid enlargement raises concern for Sjogren's syndrome"
            ),
            "plan": (
                "1. Restore carious teeth\n"
                "2. Prescription fluoride toothpaste (5000 ppm)\n"
                "3. Saliva substitute (Biotene)\n"
                "4. Recommend medical evaluation for Sjogren's syndrome -- communicate with PCP\n"
                "5. Return in 3 months for caries check"
            ),
        },
    )

    generate_progress_note(
        os.path.join(output_dir, "2022-09-19_PCP_Fatigue.pdf"),
        LINDA, "09/19/2022", "pcp",
        "OFFICE VISIT -- FATIGUE, DRY EYES, DRY MOUTH, JOINT PAIN",
        {
            "chief_complaint": "Progressive fatigue, dry eyes (under ophthalmology care), dry mouth (dentist noted concern for Sjogren's), new joint pain.",
            "hpi": (
                "Ms. Chen is a 49-year-old teacher with escalating complaints: (1) dry eyes treated by "
                "Ophthalmology with Restasis and punctal plugs; (2) dry mouth leading to dental caries, "
                "dentist raised concern for Sjogren's; (3) new bilateral hand and knee pain x 3 months "
                "with morning stiffness 20-30 minutes; (4) debilitating fatigue; (5) intermittent "
                "bilateral parotid swelling. No rash, no fever, no weight loss."
            ),
            "assessment": (
                "1. Sicca symptoms (dry eyes + dry mouth) with bilateral parotid enlargement, fatigue, "
                "and arthralgia -- Sjogren's syndrome is a strong consideration\n"
                "2. Arthralgias -- inflammatory vs. osteoarthritis"
            ),
            "plan": (
                "1. ANA, anti-SSA/Ro, anti-SSB/La, RF, CBC, CMP, CRP, ESR\n"
                "2. Rheumatology referral if serologies positive\n"
                "3. Continue ophthalmology and dental management\n"
                "4. Return with results"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2022-09-19_Lab_Sjogren_Screen.pdf"),
        LINDA, "09/19/2022", PROVIDERS["pcp"],
        [{
            "panel_name": "AUTOIMMUNE PANEL",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:320", "unit": "", "ref_range": "<1:40", "flag": "A"},
                {"test": "ANA Pattern", "value": "Speckled", "unit": "", "ref_range": "N/A", "flag": ""},
                {"test": "Anti-SSA/Ro", "value": "Positive (>8.0)", "unit": "AI", "ref_range": "<1.0", "flag": "HH"},
                {"test": "Anti-SSB/La", "value": "Positive (3.2)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
                {"test": "RF", "value": "48", "unit": "IU/mL", "ref_range": "<14", "flag": "H"},
                {"test": "CRP", "value": "4.2", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "32", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }],
    )

    generate_progress_note(
        os.path.join(output_dir, "2022-11-14_Rheum_Sjogrens_Dx.pdf"),
        LINDA, "11/14/2022", "rheumatology",
        "RHEUMATOLOGY NEW PATIENT -- SJOGREN'S SYNDROME DIAGNOSIS",
        {
            "chief_complaint": "Dry eyes, dry mouth, positive anti-SSA/SSB. Suspected Sjogren's syndrome.",
            "hpi": (
                "Ms. Chen referred with sicca symptoms, bilateral parotid enlargement, fatigue, arthralgias, "
                "ANA 1:320 speckled, anti-SSA positive (>8.0), anti-SSB positive (3.2), RF 48."
            ),
            "exam": (
                "HEENT: Dry oral mucosa. Bilateral parotid enlargement. Schirmer test <5 mm bilateral "
                "(per ophthalmology records).\n"
                "Musculoskeletal: Mild tenderness MCP and PIP joints bilaterally. No synovitis.\n"
                "Skin: No purpura. No vasculitic lesions."
            ),
            "assessment": (
                "2016 ACR/EULAR Classification Criteria for Sjogren's Syndrome:\n"
                "- Anti-SSA/Ro positive: 3 points\n"
                "- Schirmer test <= 5mm: 1 point\n"
                "- Lip biopsy: not yet performed (would add 3 points if focus score >= 1)\n"
                "Current score: 4 points -- MEETS CRITERIA (>=4 required)\n\n"
                "Diagnosis: Primary Sjogren's Syndrome\n"
                "Time from symptom onset to diagnosis: approximately 3 years (dry eyes since ~2020)."
            ),
            "plan": (
                "1. Lip biopsy for histopathologic confirmation (optional, criteria already met)\n"
                "2. Hydroxychloroquine 200 mg BID for fatigue and arthralgias\n"
                "3. Pilocarpine 5 mg TID for xerostomia\n"
                "4. Continue Restasis for dry eyes\n"
                "5. Baseline chest X-ray (ILD screening)\n"
                "6. Renal function monitoring (renal tubular acidosis risk)\n"
                "7. Monitor for lymphoma risk (4-44x increased in Sjogren's)\n"
                "8. Follow-up in 3 months"
            ),
        },
    )

    # Lip biopsy
    generate_pathology_report(
        os.path.join(output_dir, "2022-12-05_Path_Lip_Biopsy.pdf"),
        LINDA, "12/05/2022",
        "Lip, lower, minor salivary gland biopsy",
        "52-year-old woman with Sjogren's syndrome (anti-SSA+, anti-SSB+, sicca symptoms). Confirmatory biopsy.",
        "Received: lip tissue measuring 0.8 x 0.5 x 0.3 cm with visible minor salivary glands.",
        (
            "Four minor salivary gland lobules are present. Three of four lobules show dense "
            "periductal and parenchymal lymphocytic infiltrates. Focus score: 3 foci per 4 mm^2 "
            "(abnormal; >=1 focus per 4 mm^2 is positive). Foci consist predominantly of CD4+ "
            "T lymphocytes with admixed CD20+ B lymphocytes. No granulomas. No evidence of "
            "lymphoma. Acinar atrophy is present in affected lobules."
        ),
        (
            "FOCAL LYMPHOCYTIC SIALADENITIS -- CONSISTENT WITH SJOGREN'S SYNDROME\n"
            "Focus score: 3 (strongly positive; threshold >=1)\n"
            "No evidence of mucosa-associated lymphoid tissue (MALT) lymphoma."
        ),
    )

    print(f"  Patient Linda Chen: {len(os.listdir(output_dir))} documents generated")


# =====================================================================
# Patient 5: Rachel Thompson -- 38F, 5-year MCTD odyssey
# =====================================================================

RACHEL = {
    "name": "Rachel Thompson",
    "dob": "1987-09-28",
    "mrn": "RTH-2020-78234",
    "sex": "F",
    "age_at_start": 32,
}


def generate_rachel(output_dir: str):
    """Generate documents for Rachel Thompson -- MCTD."""
    os.makedirs(output_dir, exist_ok=True)

    generate_progress_note(
        os.path.join(output_dir, "2020-01-20_PCP_Raynauds.pdf"),
        RACHEL, "01/20/2020", "pcp",
        "OFFICE VISIT -- RAYNAUD'S PHENOMENON, FATIGUE",
        {
            "chief_complaint": "Color changes in fingers with cold exposure x 6 months. Fatigue.",
            "hpi": (
                "Ms. Thompson is a 32-year-old accountant presenting with episodic color changes in "
                "fingers (white -> blue -> red) triggered by cold exposure and stress x 6 months. "
                "Episodes last 15-30 minutes, sometimes painful. She also reports increasing fatigue, "
                "puffy fingers in the morning, and mild joint stiffness in hands lasting 20 minutes."
            ),
            "exam": (
                "Hands: No active color changes today. Mild non-pitting edema of fingers ('puffy fingers'). "
                "No sclerodactyly. Nailfold capillaries: not formally assessed.\n"
                "Joints: Mild tenderness MCP and PIP joints bilaterally. No synovitis."
            ),
            "assessment": (
                "1. Raynaud's phenomenon -- primary vs. secondary (autoimmune)\n"
                "2. Puffy fingers -- nonspecific but concerning in context of Raynaud's\n"
                "3. Arthralgias\n"
                "4. Fatigue"
            ),
            "plan": (
                "1. ANA, CBC, CMP, CRP, ESR\n"
                "2. Nifedipine 30 mg ER daily for Raynaud's\n"
                "3. Warm gloves, avoid cold triggers\n"
                "4. Return with results in 2-3 weeks"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2020-01-20_Lab_ANA.pdf"),
        RACHEL, "01/20/2020", PROVIDERS["pcp"],
        [{
            "panel_name": "AUTOIMMUNE SCREENING",
            "results": [
                {"test": "ANA by IIF", "value": "Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "ANA Titer", "value": "1:1280", "unit": "", "ref_range": "<1:40", "flag": "HH"},
                {"test": "ANA Pattern", "value": "Speckled", "unit": "", "ref_range": "N/A", "flag": ""},
                {"test": "CRP", "value": "6.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "38", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
            ],
        }],
    )

    # Years of scattered visits -- muscle pain, swelling, different specialists
    generate_progress_note(
        os.path.join(output_dir, "2020-09-14_PCP_Muscle_Pain.pdf"),
        RACHEL, "09/14/2020", "pcp",
        "OFFICE VISIT -- MUSCLE PAIN, WORSENING FATIGUE",
        {
            "chief_complaint": "Proximal muscle weakness and pain x 2 months. Worsening fatigue.",
            "hpi": (
                "Ms. Thompson returns with new complaints of proximal muscle pain and weakness. She "
                "has difficulty climbing stairs and lifting objects overhead. Raynaud's continues. "
                "Puffy fingers persist. ANA was positive 1:1280 speckled in January -- she was told "
                "'it could be nothing' and to follow up if symptoms worsened."
            ),
            "exam": (
                "Motor: Proximal weakness -- hip flexors 4/5 bilateral, deltoids 4/5 bilateral. "
                "Distal strength 5/5. No atrophy.\n"
                "Hands: Persistent puffy fingers. Mild Raynaud's changes (cool, pale fingertips).\n"
                "Joints: Mild synovitis bilateral MCPs."
            ),
            "assessment": (
                "1. Proximal myopathy -- inflammatory vs. metabolic. CK, aldolase needed.\n"
                "2. ANA 1:1280 speckled with Raynaud's, puffy fingers, arthritis, and now myositis "
                "features -- this pattern raises concern for MCTD or overlap syndrome.\n"
                "3. URGENT rheumatology referral needed."
            ),
            "plan": (
                "1. CK, aldolase, LDH, AST, ALT\n"
                "2. Anti-U1-RNP antibody (key test for MCTD)\n"
                "3. Anti-dsDNA, anti-Smith, anti-Scl-70, anti-Jo-1\n"
                "4. URGENT rheumatology referral\n"
                "5. Prednisone 20 mg daily if CK markedly elevated"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2020-09-14_Lab_Myositis.pdf"),
        RACHEL, "09/14/2020", PROVIDERS["pcp"],
        [{
            "panel_name": "MUSCLE ENZYMES",
            "results": [
                {"test": "CK, Total", "value": "842", "unit": "U/L", "ref_range": "30-135", "flag": "HH"},
                {"test": "Aldolase", "value": "12.8", "unit": "U/L", "ref_range": "1.0-7.5", "flag": "H"},
                {"test": "LDH", "value": "310", "unit": "U/L", "ref_range": "120-246", "flag": "H"},
                {"test": "AST", "value": "68", "unit": "U/L", "ref_range": "10-40", "flag": "H"},
                {"test": "ALT", "value": "42", "unit": "U/L", "ref_range": "7-56", "flag": ""},
            ],
        }, {
            "panel_name": "AUTOIMMUNE PANEL",
            "results": [
                {"test": "Anti-U1-RNP", "value": "Positive (>8.0)", "unit": "AI", "ref_range": "<1.0", "flag": "HH"},
                {"test": "ANA Titer", "value": "1:2560", "unit": "", "ref_range": "<1:40", "flag": "HH"},
                {"test": "Anti-dsDNA", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-Smith", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-Scl-70", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "Anti-Jo-1", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
            ],
        }],
    )

    generate_progress_note(
        os.path.join(output_dir, "2020-11-02_Rheum_MCTD_Diagnosis.pdf"),
        RACHEL, "11/02/2020", "rheumatology",
        "RHEUMATOLOGY NEW PATIENT -- MIXED CONNECTIVE TISSUE DISEASE",
        {
            "chief_complaint": "Raynaud's, puffy fingers, polyarthritis, proximal myopathy, anti-U1-RNP strongly positive.",
            "hpi": (
                "Ms. Thompson is a 33-year-old woman with 10-month history of progressive multisystem "
                "inflammatory disease: Raynaud's phenomenon, puffy fingers (sclerodactyly absent), "
                "inflammatory polyarthritis, proximal myopathy with CK 842, and strongly positive "
                "anti-U1-RNP (>8.0 AI) with ANA 1:2560 speckled. Anti-dsDNA, anti-Smith, anti-Scl-70, "
                "and anti-Jo-1 are all negative."
            ),
            "assessment": (
                "Alarcon-Segovia Criteria for MCTD:\n"
                "Serologic: Anti-U1-RNP >= 1:1600 -- YES\n"
                "Clinical (need 3 of 5):\n"
                "1. Edema of hands (puffy fingers) -- YES\n"
                "2. Synovitis -- YES (bilateral MCP synovitis)\n"
                "3. Myositis (clinical + CK elevation) -- YES\n"
                "4. Raynaud's phenomenon -- YES\n"
                "5. Acrosclerosis -- NO\n"
                "4 of 5 clinical criteria met -- DIAGNOSTIC FOR MCTD\n\n"
                "Diagnosis: Mixed Connective Tissue Disease (MCTD)\n\n"
                "Note: MCTD may evolve over time into a more defined CTD (SLE, SSc, or PM). "
                "Requires longitudinal monitoring."
            ),
            "plan": (
                "1. Prednisone 40 mg daily taper over 8 weeks (for active myositis)\n"
                "2. Hydroxychloroquine 200 mg BID (for arthritis, Raynaud's, long-term immunomodulation)\n"
                "3. Continue nifedipine for Raynaud's\n"
                "4. Baseline PFTs and HRCT chest (ILD screening -- significant risk in MCTD)\n"
                "5. Echocardiogram (pulmonary hypertension screening -- leading cause of mortality in MCTD)\n"
                "6. If myositis does not respond to steroids, add methotrexate or azathioprine\n"
                "7. Follow-up in 6 weeks with CK, inflammatory markers"
            ),
        },
    )

    # PFTs and HRCT
    generate_imaging_report(
        os.path.join(output_dir, "2020-11-16_HRCT_Chest.pdf"),
        RACHEL, "11/16/2020", "rheumatology",
        "CT", "Chest HRCT",
        "MCTD with ILD screening. Anti-U1-RNP positive.",
        "High-resolution CT of the chest without contrast. Inspiratory and expiratory images.",
        (
            "Lungs: Minimal ground-glass opacity in bilateral lower lobes posteriorly. No honeycombing. "
            "No traction bronchiectasis. No consolidation. No nodules.\n"
            "Airways: Normal. No air trapping on expiratory images.\n"
            "Pleura: No effusion. No thickening.\n"
            "Mediastinum: No lymphadenopathy. Heart size normal.\n"
            "Esophagus: Mildly dilated, may suggest esophageal dysmotility."
        ),
        (
            "1. Minimal bilateral lower lobe ground-glass opacities -- may represent early NSIP "
            "(nonspecific interstitial pneumonia), the most common ILD pattern in MCTD. "
            "Clinical correlation and PFTs recommended. Repeat imaging in 6-12 months.\n"
            "2. Mildly dilated esophagus -- consider esophageal dysmotility evaluation if symptomatic.\n"
            "3. No evidence of pulmonary hypertension by imaging."
        ),
    )

    # Follow-up labs showing improvement
    generate_lab_report(
        os.path.join(output_dir, "2021-02-08_Lab_Follow_Up.pdf"),
        RACHEL, "02/08/2021", PROVIDERS["rheumatology"],
        [{
            "panel_name": "MCTD MONITORING",
            "results": [
                {"test": "CK, Total", "value": "188", "unit": "U/L", "ref_range": "30-135", "flag": "H"},
                {"test": "Aldolase", "value": "8.2", "unit": "U/L", "ref_range": "1.0-7.5", "flag": "H"},
                {"test": "CRP", "value": "3.4", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "24", "unit": "mm/hr", "ref_range": "0-20", "flag": "H"},
                {"test": "Anti-U1-RNP", "value": "Positive (6.2)", "unit": "AI", "ref_range": "<1.0", "flag": "H"},
            ],
        }],
    )

    print(f"  Patient Rachel Thompson: {len(os.listdir(output_dir))} documents generated")
