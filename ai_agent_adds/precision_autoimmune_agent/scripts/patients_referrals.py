"""Cross-specialist referral letters for all patients.

These show the bouncing between providers and the lack of cross-provider
communication that defines the diagnostic odyssey. Each referral letter
contains clinical data that SHOULD have triggered the correct diagnosis
if the receiving provider had access to the full record.

Author: Adam Jones
Date: March 2026
"""

import os

from pdf_engine import PROVIDERS, generate_referral_letter

SARAH = {
    "name": "Sarah Mitchell",
    "dob": "1991-03-22",
    "mrn": "SMI-2022-12345",
    "sex": "F",
    "age_at_start": 31,
}

MAYA = {
    "name": "Maya Rodriguez",
    "dob": "1996-08-14",
    "mrn": "MRD-2021-45678",
    "sex": "F",
    "age_at_start": 24,
}

DAVID = {
    "name": "David Park",
    "dob": "1980-11-03",
    "mrn": "DPA-2019-33102",
    "sex": "M",
    "age_at_start": 38,
}

LINDA = {
    "name": "Linda Chen",
    "dob": "1973-05-19",
    "mrn": "LCH-2022-67890",
    "sex": "F",
    "age_at_start": 47,
}

RACHEL = {
    "name": "Rachel Thompson",
    "dob": "1987-09-28",
    "mrn": "RTH-2020-78234",
    "sex": "F",
    "age_at_start": 32,
}


def generate_sarah_referrals(output_dir: str):
    """Sarah Mitchell referral letters."""
    os.makedirs(output_dir, exist_ok=True)

    # PCP -> Dermatology (2022-12-15)
    generate_referral_letter(
        os.path.join(output_dir, "2022-12-15_Referral_PCP_to_Derm.pdf"),
        SARAH, "12/15/2022",
        "pcp", "dermatology",
        PROVIDERS["dermatology"],
        "Evaluation of persistent facial rash, photosensitive distribution.",
        (
            "Ms. Mitchell is a 31-year-old woman with a 6-month history of intermittent "
            "erythematous rash over bilateral cheeks and nose. The rash is photosensitive -- "
            "worsens with sun exposure. She also reports fatigue, joint pain in her hands and "
            "knees, and a recent episode of hair thinning.\n\n"
            "Recent labs (10/03/2022): ANA positive at 1:160 (speckled pattern). CRP mildly "
            "elevated at 8.2 mg/L. ESR 28 mm/hr. CBC shows mild leukopenia (WBC 3.8) and "
            "borderline anemia (Hgb 11.8).\n\n"
            "I am referring primarily for dermatologic evaluation of the facial rash. The "
            "positive ANA is noted but at 1:160 this titer can be seen in up to 5% of "
            "healthy individuals. The clinical picture may warrant further autoimmune workup "
            "depending on your assessment.\n\n"
            "PMH: No significant medical history\n"
            "Meds: Naproxen 500 mg BID, omeprazole 20 mg daily\n"
            "Allergies: Sulfa drugs (rash)"
        ),
        (
            "1. Is this photosensitive dermatitis, rosacea, or lupus-related malar rash?\n"
            "2. Would you recommend skin biopsy?\n"
            "3. Should rheumatology referral be expedited given the ANA?"
        ),
        urgency="Routine",
    )

    # PCP -> Rheumatology (2023-03-28) -- AFTER the ER visit
    generate_referral_letter(
        os.path.join(output_dir, "2023-03-28_Referral_PCP_to_Rheum.pdf"),
        SARAH, "03/28/2023",
        "pcp", "rheumatology",
        PROVIDERS["rheumatology"],
        "URGENT evaluation for suspected systemic lupus erythematosus.",
        (
            "Ms. Mitchell is a 31-year-old woman with a progressively concerning clinical "
            "picture that I believe may represent systemic lupus erythematosus. I am requesting "
            "urgent evaluation.\n\n"
            "CLINICAL TIMELINE:\n"
            "- 06/2022: Presented with fatigue, joint pain. Labs: mild leukopenia, borderline anemia\n"
            "- 10/2022: ANA positive 1:160 speckled, CRP 8.2, ESR 28. Facial rash developing\n"
            "- 12/2022: Urgent care visit for pleuritic chest pain (diagnosed costochondritis)\n"
            "- 01/2023: Dermatology evaluated malar rash -- photosensitive, recommended biopsy "
            "(pending)\n"
            "- 02/2023: ED visit with polyarthritis, malar rash, fever, bilateral pleural "
            "effusions, leukopenia (WBC 3.2), lymphopenia (0.8), oral ulcers. Diagnosed as "
            "'viral syndrome' and discharged on short prednisone burst.\n"
            "- 03/2023: ANA now 1:320 speckled (rising titer). Anti-dsDNA pending.\n\n"
            "She has responded dramatically to the prednisone burst but symptoms are returning "
            "as she tapers. In retrospect, I believe the ER presentation in February was a "
            "lupus flare, not a viral syndrome.\n\n"
            "CURRENT MEDICATIONS: Naproxen, omeprazole, hydroxyzine, triamcinolone cream, "
            "vitamin D, iron supplement\n"
            "ALLERGIES: Sulfa drugs (rash)"
        ),
        (
            "1. Does this patient meet ACR/EULAR criteria for SLE?\n"
            "2. Is anti-dsDNA, complement, and renal assessment needed urgently?\n"
            "3. Should she remain on corticosteroids pending your evaluation?"
        ),
        urgency="URGENT",
    )

    print("    + Sarah referral letters: 2 documents added")


def generate_maya_referrals(output_dir: str):
    """Maya Rodriguez referral letters."""
    os.makedirs(output_dir, exist_ok=True)

    # Psychiatry consultation letter BACK to PCP (2021-10-18)
    # This is the letter that acknowledges "maybe it's not just anxiety"
    generate_referral_letter(
        os.path.join(output_dir, "2021-10-18_Letter_Psychiatry_to_PCP.pdf"),
        MAYA, "10/18/2021",
        "psychiatry", "pcp",
        PROVIDERS["pcp"],
        (
            "Consultation follow-up for Maya Rodriguez. Partial response to psychiatric "
            "treatment. Possible organic component to symptoms."
        ),
        (
            "Dear Dr. Martinez,\n\n"
            "Thank you for referring Ms. Rodriguez for psychiatric evaluation. I have now seen "
            "her twice (09/13/2021 and 10/11/2021).\n\n"
            "DIAGNOSES: Generalized Anxiety Disorder (F41.1), Panic Disorder (F41.0)\n\n"
            "TREATMENT: Sertraline 50 mg daily, hydroxyzine 25 mg PRN, trazodone 50 mg QHS PRN\n\n"
            "RESPONSE: Partial. Her anxiety about the episodes has improved modestly, but the "
            "SOMATIC SYMPTOMS (tachycardia, dizziness, near-syncope) have NOT improved with "
            "SSRI treatment. In fact, she reports worsened dizziness since starting sertraline. "
            "She also developed a transient rash on sertraline.\n\n"
            "CONCERNS:\n"
            "While she meets criteria for GAD and panic disorder, I have several concerns that "
            "suggest an organic component may be contributing:\n"
            "1. Her tachycardia is consistently POSITIONAL (standing) -- she has been tracking "
            "with a smartwatch showing resting HR 65-75, standing HR 120-140\n"
            "2. Symptoms are worse with heat, after meals, and in the morning -- not typical "
            "anxiety triggers\n"
            "3. She feels better lying down -- anxiety does not typically improve with position\n"
            "4. The SSRI made somatic symptoms WORSE, which is unusual for panic disorder "
            "(SSRIs typically reduce somatic anxiety symptoms by week 4-6)\n"
            "5. She reports associated GI symptoms, flushing, and hives that predate anxiety\n"
            "6. She has notably hypermobile joints (I observed finger hyperextension)\n\n"
            "RECOMMENDATION:\n"
            "I recommend CARDIOLOGY REFERRAL with tilt table testing to evaluate for postural "
            "orthostatic tachycardia syndrome (POTS). If POTS is confirmed, many of her 'anxiety' "
            "symptoms would be reattributed to autonomic dysfunction, and treatment would be "
            "fundamentally different (volume expansion, compression, exercise rather than SSRI).\n\n"
            "I will continue psychiatric management for the secondary anxiety symptoms, but "
            "I believe the primary driver is medical, not psychiatric.\n\n"
            "If POTS is confirmed, I will reassess psychiatric diagnoses. The anxiety and "
            "avoidance behaviors may be SECONDARY to a frightening medical condition rather "
            "than primary psychiatric disorders."
        ),
        urgency="Routine",
    )

    # PCP -> Cardiology referral (2022-01-24) -- finally
    generate_referral_letter(
        os.path.join(output_dir, "2022-01-24_Referral_PCP_to_Cardiology.pdf"),
        MAYA, "01/24/2022",
        "pcp", "cardiology",
        PROVIDERS["cardiology"],
        (
            "Evaluation of positional tachycardia and recurrent presyncope. Possible POTS. "
            "Psychiatry recommends tilt table testing."
        ),
        (
            "Ms. Rodriguez is a 25-year-old graduate student with 10 months of recurrent "
            "episodes of tachycardia, dizziness, presyncope, and one syncopal episode. She has "
            "been under psychiatric care for anxiety/panic disorder since September 2021, but "
            "her psychiatrist (Dr. Crawford) now believes there may be an underlying autonomic "
            "component and specifically recommends tilt table testing.\n\n"
            "KEY CLINICAL FEATURES:\n"
            "- Tachycardia is POSITIONAL: resting HR 65-75, standing HR 120-140 (per smartwatch)\n"
            "- Symptoms improve with recumbency\n"
            "- Worse with heat, post-prandial, menstruation\n"
            "- SSRI (sertraline) worsened symptoms\n"
            "- ER visits x2 (07/2021 syncope, 11/2021 tachycardia + flushing + tryptase 18.4)\n"
            "- No structural heart disease on prior ER EKG\n\n"
            "ASSOCIATED SYMPTOMS:\n"
            "- GI: nausea, alternating diarrhea/constipation\n"
            "- Skin: facial flushing, intermittent urticaria\n"
            "- Neuro: brain fog, difficulty concentrating\n"
            "- MSK: joint hypermobility, chronic joint pain\n\n"
            "VITALS FROM TODAY: Seated BP 110/64, HR 68. Standing (3 min): BP 102/70, HR 104 "
            "(delta +36 bpm)\n\n"
            "CURRENT MEDS: Sertraline 50 mg, metoprolol ER 25 mg, omeprazole 20 mg, "
            "ondansetron PRN, hydroxyzine PRN\n\n"
            "ALLERGIES: NKDA. Reports flushing with alcohol (intolerance)."
        ),
        (
            "1. Please evaluate for postural orthostatic tachycardia syndrome (POTS)\n"
            "2. Tilt table testing requested per psychiatry recommendation\n"
            "3. Echocardiogram to rule out structural cardiac disease\n"
            "4. Should metoprolol be continued, adjusted, or changed?"
        ),
        urgency="URGENT",
    )

    # Ophthalmology letter to PCP about uveitis + back pain (David)
    # This is critical: ophthalmologist sees the connection no one else has
    generate_referral_letter(
        os.path.join(output_dir, "2022-09-12_Letter_Allergy_to_PCP.pdf"),
        MAYA, "09/12/2022",
        "allergy", "pcp",
        PROVIDERS["pcp"],
        (
            "Consultation update: Evaluation of elevated tryptase and recurrent flushing/urticaria. "
            "Mast cell activation syndrome is being considered."
        ),
        (
            "Dear Dr. Martinez,\n\n"
            "I am writing to update you on Ms. Rodriguez's evaluation in our clinic.\n\n"
            "Ms. Rodriguez was seen on 10/17/2022 for evaluation of recurrent flushing, "
            "urticaria, and an elevated serum tryptase of 18.4 ng/mL measured during an ER "
            "visit in November 2021 (normal <11.5). She was previously seen by Dermatology "
            "and diagnosed with rosacea and chronic urticaria.\n\n"
            "Our workup reveals:\n"
            "- Baseline tryptase: 14.2 ng/mL (elevated, suggesting persistent mast cell "
            "activation or hereditary alpha-tryptasemia)\n"
            "- 24-hour urine histamine: 86 mcg (elevated, normal <64)\n"
            "- 24-hour urine prostaglandin D2: 1,240 ng (elevated, normal <100-280)\n"
            "- Urine N-methylhistamine: 312 mcg/g Cr (elevated, normal <200)\n\n"
            "These findings are consistent with MAST CELL ACTIVATION SYNDROME (MCAS). The "
            "elevated baseline tryptase also raises the possibility of hereditary alpha-"
            "tryptasemia (HaT) from TPSAB1 gene duplication, which should be confirmed with "
            "genetic testing.\n\n"
            "IMPORTANTLY: MCAS frequently co-occurs with POTS (which she has, diagnosed by "
            "Dr. Thompson in February 2022) and hypermobile Ehlers-Danlos syndrome (hEDS). "
            "This triad (POTS + hEDS + MCAS) is an increasingly recognized clinical entity. "
            "I recommend evaluation for hEDS using the 2017 diagnostic criteria.\n\n"
            "I have started her on:\n"
            "- Cetirizine 10 mg BID (H1 blocker)\n"
            "- Famotidine 20 mg BID (H2 blocker)\n"
            "- Cromolyn sodium 200 mg QID (mast cell stabilizer)\n"
            "- Avoid known triggers (alcohol, heat, NSAIDs, certain foods)\n\n"
            "I will see her again in 3 months to assess response."
        ),
        urgency="Routine",
    )

    print("    + Maya referral letters: 3 documents added")


def generate_david_referrals(output_dir: str):
    """David Park referral letters."""
    os.makedirs(output_dir, exist_ok=True)

    # Ophthalmology -> PCP letter (2021-06-21) -- the critical missed connection
    generate_referral_letter(
        os.path.join(output_dir, "2021-06-21_Letter_Ophth_to_PCP.pdf"),
        DAVID, "06/21/2021",
        "ophthalmology", "pcp",
        PROVIDERS["pcp"],
        (
            "IMPORTANT: Recurrent anterior uveitis -- third episode. "
            "HLA-B27 testing and rheumatology referral strongly recommended."
        ),
        (
            "Dear Dr. Martinez,\n\n"
            "I am writing regarding your patient Mr. Park, whom I saw on 06/14/2021 for his "
            "THIRD episode of acute anterior uveitis (left eye this time; prior episodes: "
            "right eye July 2019, right eye probable self-resolving episode early 2020).\n\n"
            "CLINICAL SIGNIFICANCE:\n"
            "Recurrent acute anterior uveitis (AAU) with ALTERNATING eyes is strongly associated "
            "with HLA-B27 positivity. HLA-B27-associated AAU is the most common form of anterior "
            "uveitis and is frequently associated with spondyloarthropathies, particularly "
            "ankylosing spondylitis.\n\n"
            "I note from his history that Mr. Park has had chronic inflammatory-type back pain "
            "for over 2 years with morning stiffness >60 minutes, bilateral achilles enthesitis, "
            "and now three episodes of AAU. This combination is HIGHLY SUSPICIOUS for ankylosing "
            "spondylitis.\n\n"
            "RECOMMENDATIONS:\n"
            "1. HLA-B27 blood test -- STRONGLY recommended\n"
            "2. Rheumatology referral -- STRONGLY recommended for SpA evaluation\n"
            "3. MRI sacroiliac joints if HLA-B27 positive or if rheumatology deems appropriate\n\n"
            "If left undiagnosed, AS leads to progressive spinal fusion and irreversible "
            "structural damage. Early biologic therapy (TNF inhibitors) can prevent structural "
            "progression. The diagnostic delay for AS averages 6-8 years -- this patient is "
            "already at 2+ years.\n\n"
            "I have treated his current uveitis episode with prednisolone acetate taper and "
            "cyclopentolate. I will see him for follow-up in 1 week.\n\n"
            "Please do not hesitate to contact me if you have questions."
        ),
        urgency="URGENT",
    )

    print("    + David referral letters: 1 document added")


def generate_linda_referrals(output_dir: str):
    """Linda Chen referral letters."""
    os.makedirs(output_dir, exist_ok=True)

    # Ophthalmology -> PCP letter (2021-09-01) -- recommending Sjogren's workup
    generate_referral_letter(
        os.path.join(output_dir, "2021-09-01_Letter_Ophth_to_PCP.pdf"),
        LINDA, "09/01/2021",
        "ophthalmology", "pcp",
        PROVIDERS["pcp"],
        (
            "Progressive severe dry eye disease with sicca symptoms. "
            "Recommending Sjogren's syndrome workup."
        ),
        (
            "Dear Dr. Martinez,\n\n"
            "I am writing regarding Ms. Chen, whom I have followed since June 2020 for "
            "progressive dry eye disease. I am concerned that her ocular findings may be "
            "part of a systemic autoimmune process.\n\n"
            "CLINICAL COURSE IN MY CARE:\n"
            "- 06/2020: Initial visit. Schirmer 8/9 mm (borderline). Started artificial tears.\n"
            "- 01/2021: Worsening. Schirmer 6/7 mm. Started Restasis. Patient mentioned dry "
            "mouth beginning.\n"
            "- 08/2021: Significantly worse. Schirmer now 4/5 mm (severe). Placed punctal plugs. "
            "Early filamentary keratitis developing. Patient reports persistent xerostomia with "
            "accelerated dental caries (per her dentist).\n\n"
            "CLINICAL CONCERN:\n"
            "The progressive decline in lacrimal function (Schirmer 8-9 -> 6-7 -> 4-5 over 14 "
            "months) combined with progressive xerostomia and accelerated dental caries is "
            "classic for SJOGREN'S SYNDROME. Her dentist (Dr. Rivera) has also independently "
            "raised this concern and noted bilateral parotid enlargement.\n\n"
            "RECOMMENDED WORKUP:\n"
            "1. ANA\n"
            "2. Anti-SSA/Ro antibody (most sensitive for Sjogren's)\n"
            "3. Anti-SSB/La antibody\n"
            "4. RF (rheumatoid factor)\n"
            "5. ESR, CRP\n"
            "6. If positive, rheumatology referral for formal diagnosis and treatment\n\n"
            "I understand that you evaluated Ms. Chen in March 2021 and attributed her symptoms "
            "to perimenopause. While hormonal changes can contribute to dryness, the SEVERITY "
            "and PROGRESSIVE NATURE of her presentation, plus bilateral parotid enlargement, "
            "is more consistent with an autoimmune etiology. I would recommend reconsidering "
            "the autoimmune workup.\n\n"
            "I will continue to manage her dry eye disease from the ophthalmology perspective."
        ),
        urgency="Routine",
    )

    # Dentistry -> PCP letter (2022-07-05) -- SECOND nudge to check for Sjogren's
    generate_referral_letter(
        os.path.join(output_dir, "2022-07-05_Letter_Dentist_to_PCP.pdf"),
        LINDA, "07/05/2022",
        "dentistry", "pcp",
        PROVIDERS["pcp"],
        (
            "Accelerated dental caries and severe xerostomia. Requesting medical evaluation "
            "for underlying cause -- Sjogren's syndrome suspected."
        ),
        (
            "Dear Dr. Martinez,\n\n"
            "I am writing regarding Ms. Chen, whom I have been following with increased "
            "frequency (every 3 months) due to an alarming pattern of accelerated dental caries.\n\n"
            "DENTAL HISTORY:\n"
            "- Prior to 2021: No caries in 10 years. Excellent oral hygiene.\n"
            "- 11/2021: 3 new cervical caries (#5, #12, #29). Restored.\n"
            "- 06/2022: 4 new cervical/root caries (#3, #14, #19, #30). Restored.\n"
            "- Total: 7 new caries in 8 months despite prescription fluoride toothpaste, "
            "Biotene rinse, and excellent compliance.\n\n"
            "All caries are at the CERVICAL (gumline) and ROOT SURFACE -- this is the classic "
            "pattern of XEROSTOMIA-INDUCED CARIES, not typical coronal caries. This pattern is "
            "highly associated with Sjogren's syndrome.\n\n"
            "PHYSICAL FINDINGS:\n"
            "- Severely dry oral mucosa despite Biotene and frequent water sipping\n"
            "- Bilateral parotid gland enlargement (progressive)\n"
            "- Fissured, dry tongue\n"
            "- Reduced saliva pooling (nearly absent)\n\n"
            "CONCERN:\n"
            "Both Dr. Wells (Ophthalmology) and I have recommended evaluation for Sjogren's "
            "syndrome. The patient tells me she was seen in March 2021 and told this was "
            "perimenopause. With respect, the severity of her presentation -- progressive "
            "destruction of dentition in a patient with excellent hygiene, bilateral parotid "
            "enlargement, and concurrent severe dry eye disease -- is NOT consistent with "
            "perimenopausal dryness alone.\n\n"
            "I strongly recommend ANA and anti-SSA/Ro testing at minimum. If positive, "
            "rheumatology referral for formal diagnosis and systemic treatment.\n\n"
            "I will continue to manage her dental care on a 3-month cycle."
        ),
        urgency="Routine",
    )

    print("    + Linda referral letters: 2 documents added")


def generate_rachel_referrals(output_dir: str):
    """Rachel Thompson referral letters."""
    os.makedirs(output_dir, exist_ok=True)

    # PCP -> Rheumatology URGENT referral (2020-09-21)
    generate_referral_letter(
        os.path.join(output_dir, "2020-09-21_Referral_PCP_to_Rheum.pdf"),
        RACHEL, "09/21/2020",
        "pcp", "rheumatology",
        PROVIDERS["rheumatology"],
        (
            "URGENT: 33-year-old woman with ANA 1:2560 speckled, strongly positive anti-U1-RNP, "
            "Raynaud's, polyarthritis, and proximal myopathy with CK 842. Suspected MCTD."
        ),
        (
            "Dear Dr. Chen,\n\n"
            "I am urgently referring Ms. Thompson, a 33-year-old accountant with a rapidly "
            "evolving multisystem inflammatory presentation.\n\n"
            "TIMELINE:\n"
            "01/2020: Presented with Raynaud's phenomenon and puffy fingers. ANA 1:1280 speckled. "
            "Started nifedipine. I initially told her the ANA 'could be nothing' and to follow up "
            "if symptoms worsened. In retrospect, I should have referred then.\n\n"
            "09/2020: Returns with NEW proximal muscle weakness (difficulty climbing stairs, "
            "overhead lifting), worsening fatigue, and now bilateral MCP synovitis.\n\n"
            "LABS FROM 09/14/2020:\n"
            "- ANA: 1:2560 speckled (RISING from 1:1280)\n"
            "- Anti-U1-RNP: Strongly positive (>8.0 AI)\n"
            "- Anti-dsDNA: Negative\n"
            "- Anti-Smith: Negative\n"
            "- Anti-Scl-70: Negative\n"
            "- Anti-Jo-1: Negative\n"
            "- CK: 842 U/L (markedly elevated)\n"
            "- Aldolase: 12.8 U/L (elevated)\n"
            "- CRP: 6.8, ESR: 38\n\n"
            "CLINICAL FEATURES:\n"
            "1. Raynaud's phenomenon (9 months)\n"
            "2. Puffy fingers / hand edema (9 months)\n"
            "3. Inflammatory polyarthritis (bilateral MCP synovitis)\n"
            "4. Proximal myopathy with CK 842\n"
            "5. Strongly positive anti-U1-RNP\n\n"
            "This pattern strongly suggests Mixed Connective Tissue Disease. I have started "
            "prednisone 20 mg daily given the active myositis and am hoping you can see her "
            "within 2-4 weeks.\n\n"
            "CURRENT MEDS: Nifedipine 30 mg ER daily, prednisone 20 mg daily (just started), "
            "ibuprofen 600 mg TID PRN\n"
            "ALLERGIES: NKDA"
        ),
        (
            "1. Does this patient have MCTD vs. overlap syndrome?\n"
            "2. Is the myositis severe enough to warrant higher-dose steroids?\n"
            "3. What steroid-sparing agent would you recommend?\n"
            "4. Does she need screening for ILD or PAH?"
        ),
        urgency="URGENT",
    )

    print("    + Rachel referral letters: 1 document added")
