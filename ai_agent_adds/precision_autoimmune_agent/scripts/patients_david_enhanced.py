"""Enhanced documents for David Park -- fill 6-year AS odyssey gaps.

Adds ~12 documents spanning the 2019-2024 gap:
- PT evaluation (2019)
- PCP follow-ups with trending symptoms (2020, 2021)
- Trending inflammatory labs (2020, 2021)
- Second uveitis episode (2021)
- GI consult for IBD screening (2023)
- Post-diagnosis rheumatology follow-ups (2022, 2023, 2024)
- Follow-up MRI showing biologic response (2024)
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_imaging_report,
    generate_lab_report,
    generate_progress_note,
)

DAVID = {
    "name": "David Park",
    "dob": "1980-11-03",
    "mrn": "DPA-2019-33102",
    "sex": "M",
    "age_at_start": 38,
}


def generate_david_enhanced(output_dir: str):
    """Generate additional documents for David Park's AS odyssey."""
    os.makedirs(output_dir, exist_ok=True)

    # ── 2019-04-29  PT Evaluation ──────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2019-04-29_PT_Evaluation.pdf"),
        DAVID, "04/29/2019", "pcp",
        "PHYSICAL THERAPY INITIAL EVALUATION",
        {
            "chief_complaint": "Low back pain referred for core strengthening program.",
            "hpi": (
                "Mr. Park is a 38-year-old male referred by Dr. Martinez for PT evaluation of "
                "chronic low back pain. Onset ~8 months ago without trauma. He describes significant "
                "morning stiffness lasting 45-60 minutes daily. Pain improves with activity and NSAIDs. "
                "X-rays show mild L4-L5 DDD. He works as a construction supervisor with moderate "
                "physical demands."
            ),
            "exam": (
                "Posture: Mild loss of lumbar lordosis. Forward head posture.\n"
                "ROM: Lumbar flexion 40 degrees (norm 60). Extension 15 degrees (norm 25). "
                "Lateral flexion reduced bilaterally.\n"
                "Special tests: Modified Schober 3 cm (reduced). Positive sacroiliac compression "
                "and distraction tests bilaterally. Gaenslen test positive bilaterally.\n"
                "Strength: Core stability 3/5. Hip extensors 4/5 bilat.\n"
                "Gait: Mildly stiff, reduced trunk rotation."
            ),
            "assessment": (
                "Chronic low back pain with inflammatory features. Significant morning stiffness "
                "and positive SI joint provocation tests suggest possible sacroiliac pathology "
                "rather than purely mechanical LBP. Noted for primary care but proceeding with "
                "stabilization program as ordered."
            ),
            "plan": (
                "1. Core stabilization program 2x/week x 12 weeks\n"
                "2. Lumbar extension exercises, hip flexor stretching\n"
                "3. Aquatic therapy 1x/week for pain management\n"
                "4. Home exercise program with emphasis on posture correction\n"
                "5. Note to Dr. Martinez: SI joint provocation tests positive bilaterally -- "
                "may warrant further workup if not improving with PT"
            ),
        },
    )

    # ── 2020-03-16  PCP Follow-up (year 2 of symptoms) ────────────
    generate_progress_note(
        os.path.join(output_dir, "2020-03-16_PCP_Follow_Up.pdf"),
        DAVID, "03/16/2020", "pcp",
        "OFFICE VISIT -- BACK PAIN FOLLOW-UP",
        {
            "chief_complaint": "Persistent back pain. PT completed, modest benefit. Morning stiffness worsening.",
            "hpi": (
                "Mr. Park returns for follow-up of chronic low back pain. He completed 12 weeks of PT "
                "with modest improvement in strength but no change in morning stiffness (now 60-90 minutes). "
                "He is taking naproxen 500 mg BID daily -- reports good relief but concerned about long-term "
                "use. New symptom: bilateral anterior chest wall pain (costochondral junction), worse with "
                "deep breathing. He had another episode of right eye redness 2 months ago -- resolved "
                "spontaneously in 1 week without treatment. No bowel changes. No psoriasis. No urethritis."
            ),
            "vitals": "BP 130/84  HR 74  Temp 98.4F  Wt 188 lbs",
            "exam": (
                "Spine: Continued reduced lumbar flexion. Modified Schober 2.5 cm. "
                "Bilateral SI joint tenderness. Reduced chest expansion 3.0 cm.\n"
                "Chest wall: Tenderness at bilateral costochondral junctions (enthesitis).\n"
                "Achilles: Bilateral insertional tenderness persists."
            ),
            "assessment": (
                "1. Chronic inflammatory back pain -- persistent despite PT and NSAIDs. Morning stiffness "
                "now 60-90 min. Multiple enthesitis sites. Recurrent eye redness (possible uveitis).\n"
                "2. Costochondral enthesitis\n"
                "3. Possible recurrent anterior uveitis -- undocumented episode\n\n"
                "Given persistent NSAID-dependent symptoms, considering MRI but will defer given current "
                "COVID-19 situation. Reassess at next visit."
            ),
            "plan": (
                "1. Continue naproxen 500 mg BID\n"
                "2. Add omeprazole 20 mg daily for GI protection\n"
                "3. Inflammatory markers: CRP, ESR\n"
                "4. Defer MRI SI joints per current COVID-19 imaging restrictions\n"
                "5. If eye redness recurs, urgent ophthalmology evaluation\n"
                "6. Return in 6 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2020-03-16_Lab_Inflammatory.pdf"),
        DAVID, "03/16/2020", PROVIDERS["pcp"],
        [{
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP", "value": "22.1", "unit": "mg/L", "ref_range": "<3.0", "flag": "HH"},
                {"test": "ESR", "value": "38", "unit": "mm/hr", "ref_range": "0-15", "flag": "H"},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "8.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "13.8", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "Hematocrit", "value": "41.2", "unit": "%", "ref_range": "38.3-48.6", "flag": ""},
                {"test": "Platelets", "value": "312", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # ── 2020-11-09  PCP Annual ─────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2020-11-09_PCP_Annual.pdf"),
        DAVID, "11/09/2020", "pcp",
        "ANNUAL PHYSICAL EXAM",
        {
            "chief_complaint": "Annual exam. Ongoing back pain.",
            "hpi": (
                "Mr. Park here for annual physical. Chief ongoing concern is chronic low back pain "
                "with morning stiffness now consistently 60-90 minutes. He has been on continuous "
                "naproxen for over a year. Pain now occasionally wakes him in the second half of the "
                "night. He endorses occasional buttock pain alternating sides. Chest wall pain "
                "comes and goes. Energy level declining -- he attributes to poor sleep from pain.\n\n"
                "PMH: Chronic low back pain, anterior uveitis (2019, possible recurrence 2020), "
                "bilateral achilles enthesitis, costochondral enthesitis\n\n"
                "Social: Construction supervisor. Married, 2 children. Non-smoker. Occasional beer."
            ),
            "vitals": "BP 132/86  HR 70  Temp 98.6F  Wt 191 lbs  BMI 27.4",
            "ros": (
                "Constitutional: Fatigue, poor sleep. HEENT: No current eye redness. "
                "MSK: Back pain, morning stiffness >60 min, chest wall pain. "
                "GI: No changes. Neuro: No deficits."
            ),
            "assessment": (
                "1. Chronic inflammatory back pain -- I am increasingly concerned this is not "
                "simple mechanical LBP. The pattern of morning stiffness >60 min, night pain, "
                "alternating buttock pain, multiple enthesitis sites, and history of anterior "
                "uveitis is suspicious for spondyloarthropathy.\n"
                "2. Hypertension, stage 1 -- partially related to chronic NSAID use\n"
                "3. Weight gain -- likely inactivity from pain"
            ),
            "plan": (
                "1. Strongly considering rheumatology referral -- will reassess after imaging\n"
                "2. MRI SI joints when available (still limited access due to COVID-19 backlog)\n"
                "3. HLA-B27 testing -- will order if imaging positive\n"
                "4. BP monitoring, consider ACEI if persistent (may be NSAID-related)\n"
                "5. Return in 6 months or sooner if worsening"
            ),
        },
    )

    # ── 2021-06-14  Ophthalmology -- Second Uveitis Episode ───────
    generate_progress_note(
        os.path.join(output_dir, "2021-06-14_Ophth_Uveitis_Recurrence.pdf"),
        DAVID, "06/14/2021", "ophthalmology",
        "OPHTHALMOLOGY VISIT -- RECURRENT ANTERIOR UVEITIS",
        {
            "chief_complaint": "Red, painful left eye x 4 days. Photophobia. History of right eye uveitis in 2019.",
            "hpi": (
                "Mr. Park presents with acute left eye pain, redness, photophobia, and tearing x 4 days. "
                "He had a prior episode of anterior uveitis in the RIGHT eye in July 2019 treated with "
                "topical steroids. He also reports a possible self-resolving episode of right eye redness "
                "in early 2020 for which he did not seek care. He has chronic low back pain and achilles "
                "tendon issues."
            ),
            "exam": (
                "OS: VA 20/40 (baseline 20/20). Conjunctival injection. 3+ anterior chamber cells, "
                "1+ flare. Fine KP. No posterior synechiae. IOP 14.\n"
                "OD: Quiet. VA 20/20. No cells. Old fine KP remnants."
            ),
            "assessment": (
                "Recurrent acute anterior uveitis, left eye. Non-granulomatous. This is now his "
                "THIRD episode (2019 OD, 2020 OD probable, 2021 OS). Alternating eyes is classic "
                "for HLA-B27-associated uveitis."
            ),
            "plan": (
                "1. Prednisolone acetate 1% q1h while awake x 3 days, then taper over 6 weeks\n"
                "2. Cyclopentolate 1% TID x 7 days\n"
                "3. HLA-B27 testing -- STRONGLY recommended given recurrent AAU and history of "
                "inflammatory back pain\n"
                "4. Rheumatology referral -- recurrent uveitis with inflammatory back pain warrants "
                "evaluation for spondyloarthropathy\n"
                "5. Follow-up in 1 week\n\n"
                "Letter sent to PCP Dr. Martinez recommending HLA-B27 and rheumatology referral."
            ),
        },
    )

    # ── 2021-07-12  PCP Follow-up -- After recurrent uveitis ─────
    generate_progress_note(
        os.path.join(output_dir, "2021-07-12_PCP_Follow_Up.pdf"),
        DAVID, "07/12/2021", "pcp",
        "OFFICE VISIT -- RECURRENT UVEITIS, BACK PAIN UPDATE",
        {
            "chief_complaint": "Follow-up after third uveitis episode. Ophthalmology letter recommending HLA-B27 and rheumatology referral.",
            "hpi": (
                "Mr. Park seen after third episode of anterior uveitis (now left eye). "
                "Dr. Wells (Ophthalmology) has sent letter strongly recommending HLA-B27 testing "
                "and rheumatology referral given recurrent uveitis + inflammatory back pain. "
                "Back pain continues, morning stiffness 60-90 min. He is now having difficulty "
                "turning his head -- new onset neck stiffness."
            ),
            "exam": (
                "Cervical spine: Mildly reduced rotation bilateral. Tenderness over cervical facets.\n"
                "Lumbar spine: Further reduced flexion. Modified Schober 2 cm.\n"
                "Chest expansion: 2.5 cm (reduced from 3.5 cm in 2019)."
            ),
            "assessment": (
                "1. Recurrent anterior uveitis (3 episodes in 2 years)\n"
                "2. Progressive inflammatory back pain with now CERVICAL involvement\n"
                "3. Progressive reduction in spinal mobility and chest expansion\n"
                "4. Multiple enthesitis sites\n\n"
                "In retrospect, this presentation has been suspicious for ankylosing spondylitis "
                "for some time. The recurrent uveitis is the final prompt for referral."
            ),
            "plan": (
                "1. CRP, ESR\n"
                "2. Rheumatology referral -- scheduled, but wait time is 6-8 months\n"
                "3. In the interim, referral to Orthopedics for MRI SI joints (shorter wait)\n"
                "4. Continue naproxen, omeprazole\n"
                "5. Return after specialist evaluations"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2021-07-12_Lab_Inflammatory.pdf"),
        DAVID, "07/12/2021", PROVIDERS["pcp"],
        [{
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP", "value": "26.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "HH"},
                {"test": "ESR", "value": "44", "unit": "mm/hr", "ref_range": "0-15", "flag": "H"},
            ],
        }],
    )

    # ── 2022-07-25  Rheum Follow-up -- 12 weeks on adalimumab ─────
    generate_progress_note(
        os.path.join(output_dir, "2022-07-25_Rheum_Follow_Up.pdf"),
        DAVID, "07/25/2022", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- 12 WEEKS ON ADALIMUMAB",
        {
            "chief_complaint": "Follow-up AS on biologic therapy. Subjective improvement.",
            "hpi": (
                "Mr. Park returns 12 weeks after starting adalimumab 40 mg SC q2w. He reports "
                "significant improvement: morning stiffness now 15-20 minutes (down from 60-90). "
                "Night pain resolved. Able to reduce naproxen to PRN (1-2x per week). Achilles "
                "pain improved. No uveitis flares since starting biologic. Energy improving. "
                "Returned to full work duties. No injection site reactions. No infections."
            ),
            "vitals": "BP 124/78  HR 68  Temp 98.4F  Wt 186 lbs",
            "exam": (
                "Cervical: Improved rotation. Still mildly reduced.\n"
                "Lumbar: Modified Schober 3 cm (improved from 2 cm pre-biologic).\n"
                "Chest expansion: 3.2 cm (improved from 2.5 cm).\n"
                "Achilles: Minimal tenderness bilateral. Costochondral: non-tender.\n"
                "No peripheral joint synovitis."
            ),
            "assessment": (
                "Ankylosing Spondylitis -- GOOD RESPONSE to TNF inhibition.\n\n"
                "BASDAI: 2.8 (down from 6.4 -- clinically significant improvement)\n"
                "ASDAS-CRP: pending lab results today\n\n"
                "Excellent clinical response at 12 weeks. Target is ASDAS <1.3 (inactive) or "
                "ASDAS improvement >=1.1 units."
            ),
            "plan": (
                "1. Continue adalimumab 40 mg SC q2w\n"
                "2. Labs: CRP, ESR, CBC, CMP, LFTs (biologic monitoring)\n"
                "3. Discontinue omeprazole (NSAID use now minimal)\n"
                "4. Continue AS-specific exercise program\n"
                "5. Annual ophthalmology screening (uveitis)\n"
                "6. Follow-up in 4 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2022-07-25_Lab_Biologic_Monitoring.pdf"),
        DAVID, "07/25/2022", PROVIDERS["rheumatology"],
        [{
            "panel_name": "INFLAMMATORY MARKERS",
            "results": [
                {"test": "CRP", "value": "4.8", "unit": "mg/L", "ref_range": "<3.0", "flag": "H"},
                {"test": "ESR", "value": "14", "unit": "mm/hr", "ref_range": "0-15", "flag": ""},
            ],
        }, {
            "panel_name": "BIOLOGIC THERAPY MONITORING",
            "results": [
                {"test": "WBC", "value": "6.8", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "14.6", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "Platelets", "value": "268", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
                {"test": "ALT", "value": "28", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "AST", "value": "24", "unit": "U/L", "ref_range": "10-40", "flag": ""},
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.7-1.3", "flag": ""},
            ],
        }],
    )

    # ── 2023-01-16  GI Consult -- IBD screening ───────────────────
    generate_progress_note(
        os.path.join(output_dir, "2023-01-16_GI_Consult.pdf"),
        DAVID, "01/16/2023", "gi",
        "GI CONSULTATION -- IBD SCREENING IN ANKYLOSING SPONDYLITIS",
        {
            "chief_complaint": "Rheumatology referral for IBD screening. AS patient on TNF inhibitor.",
            "hpi": (
                "Mr. Park is a 42-year-old man with ankylosing spondylitis (HLA-B*27:05 positive) "
                "on adalimumab, referred for IBD screening. He reports occasional loose stools "
                "(2-3x per month) and mild abdominal cramping. No bloody stools, no weight loss, "
                "no nocturnal diarrhea. He denies tenesmus or urgency. IBD is present in up to 10% "
                "of AS patients and may be subclinical."
            ),
            "exam": (
                "Abdomen: Soft, non-distended. Mild tenderness RLQ. No masses. No hepatosplenomegaly. "
                "BS active. No perianal disease."
            ),
            "assessment": (
                "1. AS patient with mild GI symptoms -- subclinical IBD is common (5-10% of AS "
                "patients; up to 60% have subclinical gut inflammation on ileocolonoscopy)\n"
                "2. Currently on adalimumab which treats both AS and IBD\n"
                "3. Mild symptoms may represent IBS vs. subclinical IBD"
            ),
            "plan": (
                "1. Fecal calprotectin (non-invasive IBD marker)\n"
                "2. CBC, CMP, iron studies, B12, folate\n"
                "3. If calprotectin elevated (>150 mcg/g), recommend colonoscopy with ileal biopsies\n"
                "4. If calprotectin normal, reassurance and symptom monitoring\n"
                "5. Follow-up with results"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2023-01-16_Lab_GI_Screening.pdf"),
        DAVID, "01/16/2023", PROVIDERS["gi"],
        [{
            "panel_name": "GI / IBD SCREENING",
            "results": [
                {"test": "Fecal Calprotectin", "value": "82", "unit": "mcg/g", "ref_range": "<50", "flag": "H"},
                {"test": "Iron, Serum", "value": "78", "unit": "mcg/dL", "ref_range": "60-170", "flag": ""},
                {"test": "Ferritin", "value": "42", "unit": "ng/mL", "ref_range": "20-250", "flag": ""},
                {"test": "Vitamin B12", "value": "488", "unit": "pg/mL", "ref_range": "200-900", "flag": ""},
                {"test": "Folate", "value": "12.4", "unit": "ng/mL", "ref_range": ">3.0", "flag": ""},
            ],
        }],
    )

    # ── 2023-07-10  Rheum Follow-up -- 15 months on biologic ─────
    generate_progress_note(
        os.path.join(output_dir, "2023-07-10_Rheum_Follow_Up.pdf"),
        DAVID, "07/10/2023", "rheumatology",
        "RHEUMATOLOGY FOLLOW-UP -- AS MAINTENANCE",
        {
            "chief_complaint": "Routine AS follow-up. Doing well on adalimumab.",
            "hpi": (
                "Mr. Park continues on adalimumab 40 mg SC q2w with sustained excellent response. "
                "Morning stiffness <15 min. Minimal back pain (2/10). No uveitis flares in over a "
                "year. Achilles and chest wall enthesitis resolved. Exercising regularly -- swimming "
                "and yoga 3x/week. GI workup showed mildly elevated calprotectin (82) but no "
                "colonoscopy needed per GI, monitoring symptomatically. No infections."
            ),
            "vitals": "BP 122/76  HR 66  Temp 98.4F  Wt 183 lbs",
            "exam": (
                "Cervical: Near-normal rotation.\n"
                "Lumbar: Modified Schober 3.5 cm (continued improvement).\n"
                "Chest expansion: 3.8 cm (near-normal).\n"
                "No enthesitis. No peripheral arthritis."
            ),
            "assessment": (
                "Ankylosing Spondylitis -- SUSTAINED REMISSION on TNF inhibition.\n\n"
                "BASDAI: 1.4 (inactive disease; <4 = inactive)\n"
                "ASDAS-CRP: 1.1 (inactive disease; <1.3 = inactive) -- TARGET MET\n\n"
                "Remarkable improvement from pre-treatment ASDAS 3.8 to current 1.1. "
                "This case illustrates the impact of diagnostic delay: 6 years of preventable "
                "disease activity and structural damage."
            ),
            "plan": (
                "1. Continue adalimumab 40 mg SC q2w -- indefinite\n"
                "2. Labs: CRP, ESR, CBC, CMP\n"
                "3. Annual ophthalmology exam\n"
                "4. Consider repeat MRI SI joints in 6 months to document structural stabilization\n"
                "5. Continue exercise program\n"
                "6. Follow-up in 6 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2023-07-10_Lab_Monitoring.pdf"),
        DAVID, "07/10/2023", PROVIDERS["rheumatology"],
        [{
            "panel_name": "AS MONITORING",
            "results": [
                {"test": "CRP", "value": "1.8", "unit": "mg/L", "ref_range": "<3.0", "flag": ""},
                {"test": "ESR", "value": "8", "unit": "mm/hr", "ref_range": "0-15", "flag": ""},
                {"test": "WBC", "value": "7.1", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": ""},
                {"test": "Hemoglobin", "value": "15.2", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "ALT", "value": "22", "unit": "U/L", "ref_range": "7-56", "flag": ""},
                {"test": "Creatinine", "value": "0.9", "unit": "mg/dL", "ref_range": "0.7-1.3", "flag": ""},
            ],
        }],
    )

    # ── 2024-01-22  Follow-up MRI SI Joints ───────────────────────
    generate_imaging_report(
        os.path.join(output_dir, "2024-01-22_MRI_SI_Joints_Follow_Up.pdf"),
        DAVID, "01/22/2024", "rheumatology",
        "MRI", "Sacroiliac Joints",
        "AS on TNF inhibitor therapy x 20 months. Follow-up to assess treatment response.",
        "MRI of the sacroiliac joints with T1, T2, STIR, and post-gadolinium sequences. "
        "Comparison: March 2022 MRI.",
        (
            "Bilateral sacroiliac joints:\n"
            "- STIR sequences: RESOLVED bone marrow edema bilaterally (previously extensive). "
            "No active inflammation identified.\n"
            "- Chronic structural changes UNCHANGED: erosions at bilateral iliac surfaces, "
            "subchondral sclerosis, partial ankylosis inferior left SI joint\n"
            "- No NEW erosions compared to prior\n"
            "- Post-gadolinium: No pathologic enhancement (previously present bilaterally)\n\n"
            "Lumbar spine: Stable vertebral squaring at L2-L3. No new syndesmophytes. "
            "Stable disc desiccation L4-L5, L5-S1."
        ),
        (
            "1. RESOLVED ACTIVE SACROILIITIS -- no bone marrow edema or pathologic enhancement "
            "(previously extensive bilateral active inflammation on March 2022 study).\n"
            "2. Stable chronic structural changes (erosions, sclerosis, partial ankylosis) -- "
            "these are irreversible but no progression.\n"
            "3. No new syndesmophyte formation.\n\n"
            "IMPRESSION: Excellent imaging response to TNF inhibitor therapy with resolution of "
            "active inflammation and no structural progression."
        ),
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient David Park: {count} documents total (enhanced)")
