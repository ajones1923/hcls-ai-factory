"""Medication reconciliation documents for all 5 patients.

These show the accumulating symptomatic treatments without a unifying diagnosis.
When the AI ingests these alongside clinical notes, it can flag:
"Patient is on 9 medications that map to a single underlying condition."

The medication pile IS the diagnostic signal that no single provider saw.

Author: Adam Jones
Date: March 2026
"""

import os

from pdf_engine import generate_medication_list

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


def generate_sarah_med_list(output_dir: str):
    """Sarah Mitchell -- medication pile BEFORE lupus diagnosis (April 2023).

    9 medications treating lupus symptoms individually. No provider had connected
    the dots: joint pain + rash + pleurisy + fatigue + cytopenias = SLE.
    """
    os.makedirs(output_dir, exist_ok=True)

    generate_medication_list(
        os.path.join(output_dir, "2023-03-15_Medication_Reconciliation.pdf"),
        SARAH, "03/15/2023", "pcp",
        [
            {
                "name": "Naproxen",
                "dose": "500 mg PO",
                "frequency": "BID",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "06/2022",
                "indication": "Joint pain, hands and knees",
            },
            {
                "name": "Omeprazole",
                "dose": "20 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "06/2022",
                "indication": "GI protection (NSAID use)",
            },
            {
                "name": "Hydroxyzine",
                "dose": "25 mg PO",
                "frequency": "BID PRN",
                "prescriber": "Dr. Foster (Derm)",
                "start_date": "01/2023",
                "indication": "Facial rash / itching",
            },
            {
                "name": "Triamcinolone 0.1% cream",
                "dose": "Topical",
                "frequency": "BID to face",
                "prescriber": "Dr. Foster (Derm)",
                "start_date": "01/2023",
                "indication": "Malar rash",
            },
            {
                "name": "Prednisone",
                "dose": "20 mg PO taper",
                "frequency": "Daily",
                "prescriber": "ER physician",
                "start_date": "02/2023",
                "indication": "Polyarthritis (viral syndrome dx)",
            },
            {
                "name": "Sunscreen SPF 50",
                "dose": "Topical",
                "frequency": "Daily",
                "prescriber": "Dr. Foster (Derm)",
                "start_date": "01/2023",
                "indication": "Photosensitive rash",
            },
            {
                "name": "Vitamin D3",
                "dose": "2000 IU PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "10/2022",
                "indication": "Low vitamin D (22 ng/mL)",
            },
            {
                "name": "Iron supplement",
                "dose": "325 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "10/2022",
                "indication": "Borderline anemia (Hgb 11.8)",
            },
            {
                "name": "Melatonin",
                "dose": "5 mg PO",
                "frequency": "QHS",
                "prescriber": "Self / OTC",
                "start_date": "09/2022",
                "indication": "Fatigue / poor sleep quality",
            },
        ],
        allergies="Sulfa drugs (rash)",
        pharmacy_notes=(
            "Patient is on 9 active medications/supplements from 4 different prescribers. "
            "Multiple medications are for symptomatic management of what may be related "
            "conditions. Patient has upcoming rheumatology referral.\n\n"
            "Reconciliation reviewed with patient on 03/15/2023. Patient confirms taking all "
            "medications as prescribed. Prednisone taper from ER visit: completed 5-day course."
        ),
    )
    print("    + Sarah medication reconciliation: 1 document added")


def generate_maya_med_list(output_dir: str):
    """Maya Rodriguez -- medication pile BEFORE POTS diagnosis (January 2022).

    10 medications from 5 prescribers treating POTS/MCAS/hEDS symptoms as separate
    conditions: anxiety, tachycardia, nausea, GI issues, flushing, insomnia.
    The medication list itself IS the unifying diagnostic clue.
    """
    os.makedirs(output_dir, exist_ok=True)

    generate_medication_list(
        os.path.join(output_dir, "2022-01-10_Medication_Reconciliation.pdf"),
        MAYA, "01/10/2022", "pcp",
        [
            {
                "name": "Sertraline",
                "dose": "50 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Crawford (Psych)",
                "start_date": "09/2021",
                "indication": "GAD / Panic disorder",
            },
            {
                "name": "Hydroxyzine",
                "dose": "25 mg PO",
                "frequency": "PRN (2x/day max)",
                "prescriber": "Dr. Crawford (Psych)",
                "start_date": "10/2021",
                "indication": "Acute anxiety episodes",
            },
            {
                "name": "Metoprolol succinate ER",
                "dose": "25 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "08/2021",
                "indication": "Tachycardia episodes",
            },
            {
                "name": "Ondansetron",
                "dose": "4 mg PO/ODT",
                "frequency": "PRN (3x/day max)",
                "prescriber": "ER physician",
                "start_date": "11/2021",
                "indication": "Nausea with tachycardia",
            },
            {
                "name": "Omeprazole",
                "dose": "20 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "11/2021",
                "indication": "Acid reflux / nausea",
            },
            {
                "name": "Diphenhydramine",
                "dose": "25 mg PO",
                "frequency": "PRN",
                "prescriber": "Self / OTC",
                "start_date": "11/2021",
                "indication": "Flushing and hives episodes",
            },
            {
                "name": "Electrolyte packets (Liquid IV)",
                "dose": "1 packet in 16oz water",
                "frequency": "1-2 daily",
                "prescriber": "Self / OTC",
                "start_date": "08/2021",
                "indication": "Dizziness (self-managed)",
            },
            {
                "name": "Melatonin",
                "dose": "3 mg PO",
                "frequency": "QHS",
                "prescriber": "Self / OTC",
                "start_date": "07/2021",
                "indication": "Insomnia / sleep fragmentation",
            },
            {
                "name": "Ibuprofen",
                "dose": "400 mg PO",
                "frequency": "PRN",
                "prescriber": "Self / OTC",
                "start_date": "06/2021",
                "indication": "Joint pain (knees, hips)",
            },
            {
                "name": "Trazodone",
                "dose": "50 mg PO",
                "frequency": "QHS PRN",
                "prescriber": "Dr. Crawford (Psych)",
                "start_date": "10/2021",
                "indication": "Insomnia",
            },
        ],
        allergies="NKDA (note: patient reports flushing with alcohol -- classified as intolerance, not allergy)",
        pharmacy_notes=(
            "Patient is on 10 active medications/supplements from 5 different prescribers "
            "(PCP, Psychiatry, ER, 3 OTC self-managed). Significant polypharmacy for a "
            "25-year-old. Multiple medications target symptoms that could be related to a "
            "single underlying autonomic or systemic condition.\n\n"
            "NOTE: Diphenhydramine (anticholinergic) may worsen tachycardia. Metoprolol may "
            "worsen orthostatic symptoms. These potential interactions were not addressed.\n\n"
            "Patient also reports that sertraline has worsened her dizziness and nausea but "
            "was told to 'give it more time.'"
        ),
    )
    print("    + Maya medication reconciliation: 1 document added")


def generate_david_med_list(output_dir: str):
    """David Park -- medication pile BEFORE AS diagnosis (early 2022).

    7 medications from 3 prescribers treating ankylosing spondylitis symptoms
    as mechanical back pain + insomnia + eye problems. 3 years of symptom chasing.
    """
    os.makedirs(output_dir, exist_ok=True)

    generate_medication_list(
        os.path.join(output_dir, "2021-12-13_Medication_Reconciliation.pdf"),
        DAVID, "12/13/2021", "pcp",
        [
            {
                "name": "Naproxen",
                "dose": "500 mg PO",
                "frequency": "BID (continuous)",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "02/2019",
                "indication": "Chronic low back pain",
            },
            {
                "name": "Omeprazole",
                "dose": "20 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "03/2020",
                "indication": "GI protection (chronic NSAID use)",
            },
            {
                "name": "Cyclobenzaprine",
                "dose": "10 mg PO",
                "frequency": "QHS PRN",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "11/2019",
                "indication": "Back muscle spasms / stiffness",
            },
            {
                "name": "Trazodone",
                "dose": "50 mg PO",
                "frequency": "QHS",
                "prescriber": "Dr. Wu (Sleep Med)",
                "start_date": "06/2020",
                "indication": "Insomnia (pain-related sleep fragmentation)",
            },
            {
                "name": "Prednisolone acetate 1% eye drops",
                "dose": "1 drop affected eye",
                "frequency": "Taper per Ophth",
                "prescriber": "Dr. Wells (Ophth)",
                "start_date": "06/2021 (3rd Rx)",
                "indication": "Recurrent anterior uveitis",
            },
            {
                "name": "Acetaminophen",
                "dose": "1000 mg PO",
                "frequency": "PRN (max 3g/day)",
                "prescriber": "Self / OTC",
                "start_date": "2019",
                "indication": "Back pain (NSAID augmentation)",
            },
            {
                "name": "Capsaicin cream 0.075%",
                "dose": "Topical",
                "frequency": "TID to low back",
                "prescriber": "PT recommendation",
                "start_date": "05/2019",
                "indication": "Localized back pain",
            },
        ],
        allergies="NKDA",
        pharmacy_notes=(
            "Patient has been on continuous naproxen for 34 months. Chronic daily NSAID use "
            "raises GI and CV risk. Omeprazole added for gastroprotection. Three separate "
            "courses of prednisolone eye drops for recurrent uveitis episodes (2019, 2020, "
            "2021) from different Ophthalmology visits.\n\n"
            "Reconciliation note: Patient's chronic back pain has not been adequately "
            "diagnosed. He has been treated symptomatically for 3 years. The combination of "
            "chronic NSAID-dependent back pain + recurrent uveitis + insomnia from night pain "
            "involves 3 separate providers who are not communicating with each other.\n\n"
            "Annual physical scheduled for 2022 -- consider consolidating care and reassessing "
            "underlying diagnosis."
        ),
    )
    print("    + David medication reconciliation: 1 document added")


def generate_linda_med_list(output_dir: str):
    """Linda Chen -- medication pile BEFORE Sjogren's diagnosis (September 2022).

    11 medications/products treating Sjogren's symptoms as dry eyes + dry mouth +
    menopause + joint pain, from 4 different providers. The sheer volume of
    moisture-replacement products IS the diagnostic clue.
    """
    os.makedirs(output_dir, exist_ok=True)

    generate_medication_list(
        os.path.join(output_dir, "2022-08-22_Medication_Reconciliation.pdf"),
        LINDA, "08/22/2022", "pcp",
        [
            {
                "name": "Levothyroxine",
                "dose": "75 mcg PO",
                "frequency": "Daily (AM, empty stomach)",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "2018",
                "indication": "Hypothyroidism",
            },
            {
                "name": "Restasis (cyclosporine 0.05%)",
                "dose": "1 drop each eye",
                "frequency": "BID",
                "prescriber": "Dr. Wells (Ophth)",
                "start_date": "01/2021",
                "indication": "Severe dry eye disease",
            },
            {
                "name": "Artificial tears (PF)",
                "dose": "1-2 drops each eye",
                "frequency": "q1-2h PRN",
                "prescriber": "Dr. Wells (Ophth)",
                "start_date": "06/2020",
                "indication": "Dry eyes -- supplemental lubrication",
            },
            {
                "name": "Punctal plugs (lower lids)",
                "dose": "In situ",
                "frequency": "Permanent",
                "prescriber": "Dr. Wells (Ophth)",
                "start_date": "08/2021",
                "indication": "Severe aqueous-deficient dry eye",
            },
            {
                "name": "Biotene Dry Mouth Rinse",
                "dose": "15 mL oral rinse",
                "frequency": "QID",
                "prescriber": "Dr. Rivera (Dentist)",
                "start_date": "11/2021",
                "indication": "Xerostomia",
            },
            {
                "name": "PreviDent 5000 (Rx fluoride)",
                "dose": "Thin ribbon",
                "frequency": "BID (brush)",
                "prescriber": "Dr. Rivera (Dentist)",
                "start_date": "11/2021",
                "indication": "Caries prevention (xerostomia-related)",
            },
            {
                "name": "Xylitol lozenges",
                "dose": "1 lozenge",
                "frequency": "6-8x/day",
                "prescriber": "Dr. Rivera (Dentist)",
                "start_date": "11/2021",
                "indication": "Saliva stimulation",
            },
            {
                "name": "Omega-3 Fish Oil",
                "dose": "2000 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Wells (Ophth)",
                "start_date": "06/2020",
                "indication": "Dry eye -- anti-inflammatory support",
            },
            {
                "name": "Calcium + Vitamin D",
                "dose": "600/400 IU PO",
                "frequency": "BID",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "2019",
                "indication": "Bone health (perimenopausal)",
            },
            {
                "name": "Ibuprofen",
                "dose": "400 mg PO",
                "frequency": "PRN",
                "prescriber": "Self / OTC",
                "start_date": "2022",
                "indication": "Joint stiffness (hands, knees)",
            },
            {
                "name": "Glucosamine/Chondroitin",
                "dose": "1500/1200 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "03/2021",
                "indication": "Joint stiffness (perimenopause)",
            },
        ],
        allergies="Penicillin (childhood rash -- may be outgrown)",
        pharmacy_notes=(
            "Patient is on 11 active medications/products from 4 providers (PCP, Ophthalmology, "
            "Dentistry, OTC self-managed). NOTABLE: 7 of 11 items are for moisture replacement "
            "or dryness management (eyes and mouth). This pattern of multi-system exocrine "
            "gland failure in a middle-aged woman with arthralgias is highly suggestive of "
            "Sjogren's syndrome, though autoimmune workup has not been performed.\n\n"
            "Patient's ophthalmologist and dentist have BOTH independently suggested Sjogren's "
            "evaluation. PCP attributed symptoms to perimenopause in March 2021. Patient is "
            "seeking PCP re-evaluation for comprehensive workup.\n\n"
            "Pharmacy flag: Patient purchases artificial tears in bulk (12+ boxes/month). "
            "This level of consumption is unusual for age-related dry eye."
        ),
    )
    print("    + Linda medication reconciliation: 1 document added")


def generate_rachel_med_list(output_dir: str):
    """Rachel Thompson -- medication pile at diagnosis (November 2020).

    Shorter list but shows 3 prescribers treating MCTD fragments separately.
    """
    os.makedirs(output_dir, exist_ok=True)

    generate_medication_list(
        os.path.join(output_dir, "2020-10-19_Medication_Reconciliation.pdf"),
        RACHEL, "10/19/2020", "pcp",
        [
            {
                "name": "Nifedipine ER",
                "dose": "30 mg PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "01/2020",
                "indication": "Raynaud's phenomenon",
            },
            {
                "name": "Ibuprofen",
                "dose": "600 mg PO",
                "frequency": "TID PRN",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "01/2020",
                "indication": "Joint pain / swelling (hands)",
            },
            {
                "name": "Vitamin D3",
                "dose": "2000 IU PO",
                "frequency": "Daily",
                "prescriber": "Dr. Martinez (PCP)",
                "start_date": "03/2020",
                "indication": "Low vitamin D level",
            },
            {
                "name": "CBD oil (OTC)",
                "dose": "25 mg sublingual",
                "frequency": "BID",
                "prescriber": "Self / OTC",
                "start_date": "07/2020",
                "indication": "Muscle pain / weakness",
            },
            {
                "name": "Heated gloves (medical grade)",
                "dose": "N/A",
                "frequency": "As needed",
                "prescriber": "Self-purchased",
                "start_date": "02/2020",
                "indication": "Raynaud's prevention (cold weather)",
            },
            {
                "name": "Acetaminophen",
                "dose": "1000 mg PO",
                "frequency": "PRN (max 3g/day)",
                "prescriber": "Self / OTC",
                "start_date": "08/2020",
                "indication": "Muscle pain augmentation",
            },
        ],
        allergies="NKDA",
        pharmacy_notes=(
            "Patient is a 33-year-old woman on 6 medications/products for what appear to be "
            "separate musculoskeletal and vascular complaints. The combination of Raynaud's "
            "phenomenon requiring daily calcium channel blocker + inflammatory joint pain + "
            "progressive proximal muscle weakness in a young woman with ANA 1:1280 is highly "
            "suspicious for a connective tissue disease.\n\n"
            "The patient has been self-managing muscle pain with CBD oil and acetaminophen "
            "for 3 months without improvement. URGENT rheumatology referral is pending.\n\n"
            "NOTE: CK drawn at last PCP visit was 842 U/L (markedly elevated). Patient may "
            "need systemic immunosuppression. Ibuprofen alone is inadequate."
        ),
    )
    print("    + Rachel medication reconciliation: 1 document added")
