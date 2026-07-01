"""Patient: James Cooper -- 19M, Type 1 Diabetes + Celiac Disease (overlap syndrome).

Timeline: September 2024 -> March 2026
Specialists seen: ER, Endocrinology, Gastroenterology, Genetics, Dietitian
Key pattern: DKA presentation -> T1D diagnosis -> celiac screening positive ->
             duodenal biopsy Marsh 3b -> gluten-free diet

The system should detect: T1D + Celiac autoimmune overlap, HLA-DQ2/DQ8 risk,
                           anti-GAD65, anti-IA-2, anti-ZnT8, anti-tTG IgA.
"""

import os

from pdf_engine import (
    PROVIDERS,
    generate_genetic_report,
    generate_lab_report,
    generate_pathology_report,
    generate_progress_note,
)

PATIENT = {
    "name": "James Cooper",
    "dob": "2006-04-18",
    "mrn": "JCO-2024-67234",
    "sex": "M",
    "age_at_start": 18,
}


def generate(output_dir: str):
    """Generate all clinical documents for James Cooper."""
    os.makedirs(output_dir, exist_ok=True)

    # ─────────────────────────────────────────────────────────────
    # Visit 1: ER -- September 2024 (DKA presentation)
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-09-08_ER_DKA_Presentation.pdf"),
        PATIENT, "09/08/2024", "er",
        "EMERGENCY DEPARTMENT NOTE -- DIABETIC KETOACIDOSIS",
        {
            "chief_complaint": "18-year-old male with nausea, vomiting, abdominal pain, polyuria, and polydipsia x 3 weeks. Found lethargic by roommate.",
            "hpi": (
                "Mr. Cooper is an 18-year-old male college freshman brought in by ambulance after "
                "his roommate found him lethargic and confused. Per roommate, patient has had "
                "progressive polyuria, polydipsia, and unintentional weight loss of approximately "
                "15 lbs over the past month. He has been drinking excessive amounts of water and "
                "Gatorade. Today he developed nausea, vomiting (x4 episodes), and diffuse abdominal "
                "pain. He has also noted intermittent loose stools and bloating for the past 2-3 "
                "months which he attributed to dining hall food. No prior medical history. No family "
                "history of diabetes per patient (adopted, limited family history). Non-smoker, "
                "denies alcohol or drug use."
            ),
            "vitals": (
                "BP 98/62  HR 118  Temp 97.8F  RR 28 (Kussmaul)  SpO2 98% RA  "
                "Wt 142 lbs (reported prior weight 157 lbs)  Ht 5'11\"  BMI 19.8"
            ),
            "exam": (
                "General: Thin young man, lethargic but arousable. Fruity breath odor. Dry mucous membranes.\n"
                "HEENT: Sunken eyes. Dry lips and oral mucosa. No thrush.\n"
                "Cardiovascular: Tachycardic, regular rhythm. No murmurs. Capillary refill 3 seconds.\n"
                "Pulmonary: Kussmaul respirations. Clear to auscultation.\n"
                "Abdomen: Diffuse tenderness, no guarding or rigidity. Hyperactive bowel sounds.\n"
                "Skin: Decreased turgor. No acanthosis nigricans. No rash.\n"
                "Neurological: GCS 13 (E3V4M6). Oriented to person, not place or time. No focal deficits."
            ),
            "labs_reviewed": (
                "Point-of-care glucose: 486 mg/dL (CRITICAL HIGH)\n"
                "Venous blood gas: pH 7.14, pCO2 18, HCO3 8.2 -- SEVERE METABOLIC ACIDOSIS\n"
                "Serum ketones: LARGE (beta-hydroxybutyrate 6.8 mmol/L, ref <0.6)\n"
                "Anion gap: 28 (elevated)\n"
                "Na 131 (corrected 138), K 5.4, Cl 92, BUN 32, Cr 1.4 (pre-renal)\n"
                "Serum osmolality: 312 mOsm/kg\n"
                "UA: glucose 4+, ketones 3+, specific gravity 1.030"
            ),
            "assessment": (
                "1. DIABETIC KETOACIDOSIS (DKA), severe -- new-onset diabetes mellitus, likely "
                "Type 1 given age, body habitus, and presentation\n"
                "2. Severe dehydration with pre-renal azotemia\n"
                "3. Hyponatremia (corrected Na 138 -- pseudohyponatremia from hyperglycemia)\n"
                "4. Weight loss -- 15 lbs over 1 month, consistent with undiagnosed T1D\n"
                "5. GI symptoms (loose stools, bloating) -- may be DKA-related vs. need further eval"
            ),
            "plan": (
                "1. DKA protocol:\n"
                "   - NS bolus 20 mL/kg, then 250 mL/hr\n"
                "   - Insulin drip 0.1 units/kg/hr (14 units/hr)\n"
                "   - Potassium replacement per protocol\n"
                "   - Hourly BMP, glucose q1h\n"
                "   - Transition to subQ insulin when gap closes and tolerating PO\n"
                "2. Endocrinology consult for new-onset T1D management\n"
                "3. Diabetes autoantibodies: anti-GAD65, anti-IA-2, anti-ZnT8, insulin antibodies\n"
                "4. HbA1c, C-peptide, fasting lipid panel\n"
                "5. Admit to ICU for DKA monitoring\n"
                "6. Social work consult (college student, newly diagnosed)"
            ),
        },
    )

    # Lab: ER admission labs
    generate_lab_report(
        os.path.join(output_dir, "2024-09-08_Lab_DKA_Admission.pdf"),
        PATIENT, "09/08/2024", PROVIDERS["er"],
        [{
            "panel_name": "METABOLIC PANEL -- DKA",
            "results": [
                {"test": "Glucose", "value": "486", "unit": "mg/dL", "ref_range": "70-100", "flag": "HH"},
                {"test": "Sodium", "value": "131", "unit": "mEq/L", "ref_range": "136-145", "flag": "L"},
                {"test": "Potassium", "value": "5.4", "unit": "mEq/L", "ref_range": "3.5-5.0", "flag": "H"},
                {"test": "Chloride", "value": "92", "unit": "mEq/L", "ref_range": "98-106", "flag": "L"},
                {"test": "CO2", "value": "8", "unit": "mEq/L", "ref_range": "22-28", "flag": "LL"},
                {"test": "BUN", "value": "32", "unit": "mg/dL", "ref_range": "7-20", "flag": "H"},
                {"test": "Creatinine", "value": "1.4", "unit": "mg/dL", "ref_range": "0.7-1.3", "flag": "H"},
                {"test": "Anion Gap", "value": "28", "unit": "mEq/L", "ref_range": "8-12", "flag": "HH"},
                {"test": "Beta-Hydroxybutyrate", "value": "6.8", "unit": "mmol/L", "ref_range": "<0.6", "flag": "HH"},
            ],
        }, {
            "panel_name": "BLOOD GAS",
            "results": [
                {"test": "pH (venous)", "value": "7.14", "unit": "", "ref_range": "7.32-7.42", "flag": "LL"},
                {"test": "pCO2", "value": "18", "unit": "mmHg", "ref_range": "38-42", "flag": "L"},
                {"test": "HCO3", "value": "8.2", "unit": "mEq/L", "ref_range": "22-26", "flag": "LL"},
            ],
        }, {
            "panel_name": "CBC",
            "results": [
                {"test": "WBC", "value": "14.2", "unit": "K/uL", "ref_range": "4.0-11.0", "flag": "H"},
                {"test": "Hemoglobin", "value": "14.8", "unit": "g/dL", "ref_range": "13.5-17.5", "flag": ""},
                {"test": "Platelets", "value": "298", "unit": "K/uL", "ref_range": "150-400", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 2: Endocrinology Consult -- September 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-09-09_Endocrinology_Consult.pdf"),
        PATIENT, "09/09/2024", "endocrinology",
        "ENDOCRINOLOGY CONSULTATION -- NEW-ONSET TYPE 1 DIABETES",
        {
            "chief_complaint": "Inpatient endocrinology consult for new-onset DKA, suspected T1D.",
            "hpi": (
                "Mr. Cooper is an 18-year-old male admitted last night with severe DKA "
                "(pH 7.14, glucose 486, BHB 6.8). Insulin drip initiated per protocol. He is now "
                "on hospital day 1. DKA is resolving: glucose 188, anion gap closing (14), pH "
                "normalizing (7.32). He is tolerating clear liquids. We are consulted for T1D "
                "confirmation and insulin initiation. He has been having polyuria, polydipsia, "
                "and 15-lb weight loss over the past month. He also reports 2-3 months of loose "
                "stools, bloating, and occasional abdominal cramping -- unclear if related to DKA "
                "or a separate process. He is adopted and has limited family medical history."
            ),
            "assessment": (
                "1. New-onset Type 1 Diabetes Mellitus with DKA\n"
                "   - Classic presentation: young male, lean, severe DKA\n"
                "   - Autoantibodies pending (GAD65, IA-2, ZnT8, IAA)\n"
                "   - HbA1c 12.8% indicating months of undiagnosed hyperglycemia\n"
                "   - C-peptide 0.3 ng/mL (LOW -- confirms minimal residual beta cell function)\n\n"
                "2. GI symptoms (loose stools, bloating x 2-3 months)\n"
                "   - In context of new T1D, must screen for CELIAC DISEASE -- prevalence 5-10% "
                "in T1D patients due to shared HLA-DQ2/DQ8 susceptibility\n"
                "   - ADA guidelines recommend celiac screening in ALL new T1D patients\n\n"
                "3. Adopted -- recommend HLA typing for diabetes risk stratification"
            ),
            "plan": (
                "1. Transition to basal-bolus insulin:\n"
                "   - Lantus (glargine) 18 units at bedtime\n"
                "   - Humalog (lispro) 1:10 carb ratio + correction factor 1:40 above 150\n"
                "2. Diabetes autoantibody panel: anti-GAD65, anti-IA-2, anti-ZnT8, insulin antibodies\n"
                "3. CELIAC SCREENING (ADA guidelines): tTG-IgA, total IgA, anti-endomysial antibody\n"
                "4. HLA-DQ2/DQ8 typing\n"
                "5. Thyroid screening: TSH, free T4, anti-TPO (autoimmune thyroiditis comorbidity)\n"
                "6. CGM (Dexcom G7) and insulin pump discussion -- defer until stable on MDI\n"
                "7. Diabetes educator, dietitian, social work referrals\n"
                "8. Discharge when tolerating full diet on subQ insulin with stable glucose\n"
                "9. Outpatient endocrine follow-up in 1 week"
            ),
        },
    )

    # Lab: T1D autoantibodies and celiac screening
    generate_lab_report(
        os.path.join(output_dir, "2024-09-10_Lab_Autoantibodies.pdf"),
        PATIENT, "09/10/2024", PROVIDERS["endocrinology"],
        [{
            "panel_name": "DIABETES AUTOANTIBODY PANEL",
            "results": [
                {"test": "Anti-GAD65 Antibody", "value": "248", "unit": "U/mL", "ref_range": "<5.0", "flag": "HH"},
                {"test": "Anti-IA-2 Antibody", "value": "156", "unit": "U/mL", "ref_range": "<8.0", "flag": "HH"},
                {"test": "Anti-ZnT8 Antibody", "value": "Positive (>200)", "unit": "U/mL", "ref_range": "<15.0", "flag": "HH"},
                {"test": "Insulin Autoantibody (IAA)", "value": "Negative", "unit": "", "ref_range": "Negative", "flag": ""},
                {"test": "C-Peptide (fasting)", "value": "0.3", "unit": "ng/mL", "ref_range": "0.8-3.1", "flag": "LL"},
                {"test": "HbA1c", "value": "12.8", "unit": "%", "ref_range": "<5.7", "flag": "HH"},
            ],
        }, {
            "panel_name": "CELIAC DISEASE SCREENING",
            "results": [
                {"test": "tTG-IgA (tissue transglutaminase)", "value": "142", "unit": "U/mL", "ref_range": "<4.0", "flag": "HH"},
                {"test": "Total IgA", "value": "186", "unit": "mg/dL", "ref_range": "70-400", "flag": ""},
                {"test": "Anti-Endomysial Antibody (EMA)", "value": "Positive (1:160)", "unit": "", "ref_range": "Negative", "flag": "A"},
                {"test": "Deamidated Gliadin IgG", "value": "88", "unit": "U/mL", "ref_range": "<20", "flag": "HH"},
            ],
        }, {
            "panel_name": "THYROID SCREENING",
            "results": [
                {"test": "TSH", "value": "3.2", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Free T4", "value": "1.1", "unit": "ng/dL", "ref_range": "0.8-1.8", "flag": ""},
                {"test": "Anti-TPO Antibody", "value": "12", "unit": "IU/mL", "ref_range": "<35", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 3: Endocrinology Follow-up -- September 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-09-18_Endocrinology_Follow_Up.pdf"),
        PATIENT, "09/18/2024", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- T1D + CELIAC SCREENING POSITIVE",
        {
            "chief_complaint": "1-week follow-up after DKA hospitalization. Autoantibody and celiac results to review.",
            "hpi": (
                "Mr. Cooper returns 10 days post-discharge. He is managing insulin injections "
                "well with assistance from his mother (who traveled to be with him). Blood glucose "
                "averaging 180 mg/dL, range 90-310. He is carb counting using an app. He continues "
                "to have loose stools (2-3 per day), bloating, and mild abdominal cramping. He also "
                "notes fatigue despite improving glucose control.\n\n"
                "Lab results received:\n"
                "- Strongly positive diabetes autoantibodies: GAD65 248, IA-2 156, ZnT8 >200\n"
                "- CELIAC SCREENING POSITIVE: tTG-IgA 142 (>10x ULN), EMA positive 1:160\n"
                "- Thyroid normal, TPO negative"
            ),
            "assessment": (
                "1. Type 1 Diabetes Mellitus -- CONFIRMED\n"
                "   - Triple autoantibody positive (GAD65, IA-2, ZnT8)\n"
                "   - C-peptide 0.3 (minimal beta cell reserve)\n"
                "   - Glucose improving but not yet at target\n\n"
                "2. CELIAC DISEASE -- screening strongly positive\n"
                "   - tTG-IgA 142 (>10x ULN -- per ESPGHAN guidelines, tTG >10x ULN with "
                "positive EMA may be sufficient for diagnosis without biopsy in some protocols)\n"
                "   - However, given T1D overlap and importance of definitive diagnosis, "
                "recommend confirmatory duodenal biopsy\n"
                "   - GI symptoms (loose stools, bloating) now clearly explained\n"
                "   - Iron and nutritional panel needed\n\n"
                "3. Autoimmune overlap syndrome: T1D + Celiac Disease\n"
                "   - Shared HLA-DQ2/DQ8 susceptibility\n"
                "   - 5-10% of T1D patients have celiac disease\n"
                "   - Ongoing monitoring for additional autoimmune conditions (thyroiditis, Addison's)"
            ),
            "plan": (
                "1. Insulin adjustment: increase Lantus to 20 units nightly, keep carb ratio 1:10\n"
                "2. CGM (Dexcom G7) -- insurance authorization submitted\n"
                "3. GI referral for EGD with duodenal biopsy -- celiac confirmation\n"
                "4. Iron studies, ferritin, folate, B12, vitamin D, calcium (celiac nutrient panel)\n"
                "5. HLA-DQ2/DQ8 typing ordered\n"
                "6. DO NOT start gluten-free diet until AFTER biopsy (need active exposure for dx)\n"
                "7. Dietitian referral for celiac + T1D combined management\n"
                "8. Follow-up in 2 weeks with biopsy results"
            ),
        },
    )

    # Lab: Nutritional panel
    generate_lab_report(
        os.path.join(output_dir, "2024-09-18_Lab_Nutritional.pdf"),
        PATIENT, "09/18/2024", PROVIDERS["endocrinology"],
        [{
            "panel_name": "NUTRITIONAL PANEL (CELIAC SCREENING)",
            "results": [
                {"test": "Iron, Serum", "value": "38", "unit": "mcg/dL", "ref_range": "60-170", "flag": "L"},
                {"test": "Ferritin", "value": "12", "unit": "ng/mL", "ref_range": "20-250", "flag": "L"},
                {"test": "TIBC", "value": "420", "unit": "mcg/dL", "ref_range": "250-370", "flag": "H"},
                {"test": "Folate", "value": "4.8", "unit": "ng/mL", "ref_range": ">3.0", "flag": ""},
                {"test": "Vitamin B12", "value": "310", "unit": "pg/mL", "ref_range": "200-900", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "16", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
                {"test": "Calcium", "value": "8.8", "unit": "mg/dL", "ref_range": "8.6-10.2", "flag": ""},
                {"test": "Albumin", "value": "3.4", "unit": "g/dL", "ref_range": "3.5-5.0", "flag": "L"},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 4: GI Consult + EGD -- October 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-10-02_GI_Consult.pdf"),
        PATIENT, "10/02/2024", "gi",
        "GI CONSULTATION -- CELIAC DISEASE EVALUATION",
        {
            "chief_complaint": "Referred by endocrinology for EGD with duodenal biopsy. New T1D patient with strongly positive celiac serology.",
            "hpi": (
                "Mr. Cooper is an 18-year-old male with newly diagnosed T1D (DKA presentation "
                "September 2024) found to have strongly positive celiac serology: tTG-IgA 142 "
                "(>10x ULN), EMA positive 1:160. He has been having chronic GI symptoms: loose "
                "stools 2-3x/day, bloating, abdominal cramping x 2-3 months. Iron deficiency "
                "identified (serum iron 38, ferritin 12). Low vitamin D. Mild hypoalbuminemia. "
                "He is still on a gluten-containing diet as instructed for biopsy accuracy."
            ),
            "exam": (
                "Abdomen: Soft, mildly distended. Mild diffuse tenderness, no focal tenderness. "
                "No hepatosplenomegaly. BS active. No perianal disease.\n"
                "General: Thin. BMI 19.8. No dermatitis herpetiformis. No angular cheilitis. "
                "No glossitis."
            ),
            "assessment": (
                "1. Probable celiac disease in setting of T1D -- strongly positive serology\n"
                "2. Iron deficiency (pre-anemia stage) -- consistent with celiac malabsorption\n"
                "3. Vitamin D deficiency -- multifactorial (celiac, limited sun exposure)\n"
                "4. Recommend EGD with 4-6 duodenal biopsies for histologic confirmation"
            ),
            "plan": (
                "1. EGD scheduled October 7, 2024\n"
                "   - 4-6 duodenal biopsies (D1 and D2) per ACG guidelines\n"
                "   - Patient to continue gluten-containing diet until biopsy\n"
                "2. If biopsy confirms celiac: strict gluten-free diet, dietitian education\n"
                "3. Repeat tTG-IgA at 6 months and 12 months on GFD to confirm serologic response\n"
                "4. Iron supplementation: ferrous sulfate 325 mg daily\n"
                "5. Follow-up with biopsy results"
            ),
        },
    )

    # Pathology: Duodenal biopsy
    generate_pathology_report(
        os.path.join(output_dir, "2024-10-07_Pathology_Duodenal_Biopsy.pdf"),
        PATIENT, "10/07/2024",
        "Duodenal biopsies (4 fragments from D1 and D2)",
        (
            "18-year-old male with T1D and strongly positive celiac serology "
            "(tTG-IgA 142, EMA positive 1:160). Chronic diarrhea, bloating, iron deficiency."
        ),
        (
            "Received in formalin: four fragments of tan-pink mucosa measuring 0.2-0.4 cm each. "
            "All submitted in one cassette."
        ),
        (
            "Sections show duodenal mucosa with SEVERE villous blunting (subtotal villous atrophy). "
            "Crypt hyperplasia is prominent with elongated, branching crypts. The villous:crypt "
            "ratio is markedly reduced (<1:1, normal >3:1). Marked increase in intraepithelial "
            "lymphocytes (IELs), counted at >40 IELs per 100 enterocytes (normal <25). Surface "
            "enterocytes show degenerative changes with loss of brush border. Lamina propria shows "
            "increased lymphoplasmacytic infiltrate.\n\n"
            "Immunohistochemistry: CD3 highlights markedly increased intraepithelial T lymphocytes. "
            "CD8 positive predominance."
        ),
        (
            "DUODENAL BIOPSIES (D1 AND D2):\n"
            "- CELIAC DISEASE, Marsh Classification Type 3b (subtotal villous atrophy)\n"
            "  * Marked villous blunting with subtotal atrophy\n"
            "  * Crypt hyperplasia\n"
            "  * >40 IELs per 100 enterocytes\n\n"
            "Consistent with CELIAC DISEASE in clinical context of positive tTG-IgA (142) "
            "and anti-endomysial antibody."
        ),
        (
            "Marsh 3b histology indicates moderate-to-severe celiac disease. The degree of villous "
            "atrophy correlates with clinical symptoms and nutritional deficiencies. Strict gluten-free "
            "diet is the cornerstone of treatment. Follow-up biopsy may be considered at 12-24 months "
            "to confirm mucosal recovery."
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 5: HLA Typing -- October 2024
    # ─────────────────────────────────────────────────────────────
    generate_genetic_report(
        os.path.join(output_dir, "2024-10-10_Genetic_HLA_Typing.pdf"),
        PATIENT, "10/10/2024",
        "HLA-DQ Typing for Celiac Disease and Type 1 Diabetes",
        (
            "18-year-old male with newly diagnosed T1D (GAD65+, IA-2+, ZnT8+) and confirmed "
            "celiac disease (Marsh 3b). HLA typing for autoimmune risk stratification."
        ),
        "DNA extracted from peripheral blood. HLA typing by sequence-specific oligonucleotide "
        "probe (SSOP) hybridization with high-resolution confirmation by Sanger sequencing.",
        [
            {"Locus": "HLA-DQA1", "Allele 1": "*05:01", "Allele 2": "*03:01"},
            {"Locus": "HLA-DQB1", "Allele 1": "*02:01", "Allele 2": "*03:02"},
            {"Locus": "HLA-DRB1", "Allele 1": "*03:01", "Allele 2": "*04:01"},
        ],
        (
            "HETEROZYGOUS HLA-DQ2/DQ8:\n\n"
            "1. HLA-DQ2.5 haplotype detected (DQA1*05:01/DQB1*02:01 with DRB1*03:01)\n"
            "   - Present in 90-95% of celiac disease patients\n"
            "   - Also associated with T1D, SLE, and Addison's disease\n\n"
            "2. HLA-DQ8 haplotype detected (DQA1*03:01/DQB1*03:02 with DRB1*04:01)\n"
            "   - The highest-risk HLA for T1D (DRB1*03:01/DRB1*04:01 heterozygote)\n"
            "   - DQ8 found in 5-10% of celiac patients who are DQ2-negative\n\n"
            "3. DRB1*03:01/DRB1*04:01 heterozygosity confers the HIGHEST genetic risk for "
            "T1D (OR ~15-20)\n\n"
            "INTERPRETATION: This patient carries both HLA-DQ2.5 and HLA-DQ8 haplotypes, "
            "explaining susceptibility to BOTH Type 1 Diabetes and Celiac Disease. The "
            "DRB1*03:01/04:01 genotype places him at the highest genetic risk category for T1D. "
            "Recommend screening for additional autoimmune conditions (thyroiditis, Addison's disease) "
            "given polyautoimmune risk."
        ),
        (
            "1. Annual thyroid function and anti-TPO screening (autoimmune thyroiditis risk)\n"
            "2. Consider adrenal antibody screening if symptoms develop (Addison's disease risk)\n"
            "3. Genetic counseling if family planning is desired in the future\n"
            "4. First-degree relatives should be screened for celiac disease and T1D autoantibodies"
        ),
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 6: Endocrinology Follow-up -- November 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-11-12_Endocrinology_Follow_Up.pdf"),
        PATIENT, "11/12/2024", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- T1D + CELIAC DISEASE MANAGEMENT",
        {
            "chief_complaint": "Follow-up T1D and confirmed celiac disease. 5 weeks on gluten-free diet.",
            "hpi": (
                "Mr. Cooper returns for follow-up. He started a strict gluten-free diet immediately "
                "after biopsy confirmed Marsh 3b celiac disease. His GI symptoms have improved "
                "significantly -- stools now formed, 1-2x daily, bloating resolved. He has gained "
                "4 lbs (now 146 lbs). Dexcom G7 placed 2 weeks ago; time-in-range 62% (goal >70%). "
                "He reports difficulty with carb counting on a gluten-free diet as portion sizes "
                "and carb content differ from standard products. He met with the dietitian twice. "
                "No hypoglycemic episodes. HLA typing confirmed DQ2.5/DQ8 -- highest risk genotype "
                "for both T1D and celiac."
            ),
            "assessment": (
                "1. T1D -- improving control, CGM in use, TIR 62% (target >70%)\n"
                "   - HbA1c today will reflect partial improvement (only 2 months since dx)\n"
                "2. Celiac Disease (Marsh 3b) -- GI symptoms improving on GFD\n"
                "   - Weight gain is reassuring\n"
                "   - Will recheck tTG-IgA at 6 months to confirm serologic response\n"
                "3. Iron deficiency -- on supplementation\n"
                "4. Vitamin D deficiency -- on supplementation"
            ),
            "plan": (
                "1. Insulin adjustment: increase Lantus to 22 units, adjust carb ratio to 1:8 "
                "for breakfast (dawn phenomenon)\n"
                "2. Continue CGM, review download at next visit\n"
                "3. Discuss insulin pump (Tandem t:slim X2) when TIR >65% for 2 weeks\n"
                "4. Continue strict GFD\n"
                "5. Recheck: HbA1c, tTG-IgA, iron/ferritin, vitamin D at next visit\n"
                "6. Continue iron 325 mg daily and vitamin D 2000 IU daily\n"
                "7. Annual TSH and anti-TPO (autoimmune thyroiditis screening)\n"
                "8. Follow-up in 3 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2024-11-12_Lab_Follow_Up.pdf"),
        PATIENT, "11/12/2024", PROVIDERS["endocrinology"],
        [{
            "panel_name": "DIABETES MONITORING",
            "results": [
                {"test": "HbA1c", "value": "10.2", "unit": "%", "ref_range": "<5.7", "flag": "HH"},
                {"test": "C-Peptide (fasting)", "value": "0.2", "unit": "ng/mL", "ref_range": "0.8-3.1", "flag": "LL"},
                {"test": "Fasting Glucose", "value": "142", "unit": "mg/dL", "ref_range": "70-100", "flag": "H"},
            ],
        }, {
            "panel_name": "NUTRITIONAL FOLLOW-UP",
            "results": [
                {"test": "Iron, Serum", "value": "52", "unit": "mcg/dL", "ref_range": "60-170", "flag": "L"},
                {"test": "Ferritin", "value": "18", "unit": "ng/mL", "ref_range": "20-250", "flag": "L"},
                {"test": "Vitamin D, 25-OH", "value": "22", "unit": "ng/mL", "ref_range": "30-100", "flag": "L"},
                {"test": "Albumin", "value": "3.8", "unit": "g/dL", "ref_range": "3.5-5.0", "flag": ""},
            ],
        }],
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 7: Dietitian Note -- December 2024
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2024-12-04_Dietitian_GFD_T1D.pdf"),
        PATIENT, "12/04/2024", "dietitian",
        "NUTRITION CONSULTATION -- GLUTEN-FREE DIET WITH T1D",
        {
            "chief_complaint": "Follow-up nutrition counseling for concurrent T1D carb counting and celiac disease gluten-free diet.",
            "hpi": (
                "Mr. Cooper is an 18-year-old college student managing both T1D (on MDI with CGM) "
                "and celiac disease (Marsh 3b, 3 months on GFD). He has been struggling with "
                "accurate carb counting for gluten-free products, as GF breads and pastas have "
                "different glycemic indices than standard products. He reports difficulty eating "
                "in the college dining hall -- limited GF options, risk of cross-contamination. "
                "Weight: 148 lbs (gained 6 lbs since diagnosis). GI symptoms resolved on GFD."
            ),
            "assessment": (
                "1. T1D + Celiac Disease -- dual dietary management\n"
                "   - GF products often have higher glycemic index -- contributing to post-meal spikes\n"
                "   - Need to focus on naturally GF whole grains (rice, quinoa, millet) vs. processed GF products\n"
                "   - College dining hall challenges identified -- will work with university dining services\n"
                "2. Nutritional recovery underway -- weight gain, improving iron and vitamin D\n"
                "3. Adequate protein intake at 0.9 g/kg/day (goal >1.0)"
            ),
            "plan": (
                "1. Updated GF carb counting guide -- with GF-specific glycemic index adjustments\n"
                "2. Meal plans using naturally GF whole foods (less reliance on processed GF products)\n"
                "3. Letter to university dining services requesting GF accommodations\n"
                "4. Taught label reading for hidden gluten (malt, modified food starch, soy sauce)\n"
                "5. Snack strategies for hypoglycemia using GF options\n"
                "6. Follow-up in 2 months"
            ),
        },
    )

    # ─────────────────────────────────────────────────────────────
    # Visit 8: Endocrinology 6-month -- March 2025
    # ─────────────────────────────────────────────────────────────
    generate_progress_note(
        os.path.join(output_dir, "2025-03-10_Endocrinology_6Month.pdf"),
        PATIENT, "03/10/2025", "endocrinology",
        "ENDOCRINOLOGY FOLLOW-UP -- 6 MONTHS T1D AND CELIAC",
        {
            "chief_complaint": "6-month follow-up T1D and celiac disease. On insulin pump and GFD.",
            "hpi": (
                "Mr. Cooper returns for 6-month review. He transitioned to Tandem t:slim X2 insulin "
                "pump with Control-IQ technology 2 months ago. Time-in-range now 72% (improved from "
                "62%). HbA1c improving. He has been strictly adherent to GFD -- no GI symptoms. "
                "GF dining accommodations arranged at university. Weight stable at 152 lbs. He "
                "reports good energy, mood improved, and academic performance back to baseline. "
                "No hypoglycemic events requiring assistance. Annual thyroid screening due."
            ),
            "assessment": (
                "1. T1D -- IMPROVING\n"
                "   - HbA1c 7.4% (down from 12.8% at diagnosis, 10.2% at 2 months)\n"
                "   - TIR 72% -- meeting ADA target (>70%)\n"
                "   - Control-IQ hybrid closed loop performing well\n"
                "   - C-peptide undetectable (<0.1) -- honeymoon period has ended\n\n"
                "2. Celiac Disease -- GOOD DIETARY RESPONSE\n"
                "   - tTG-IgA 18 (down from 142 at diagnosis) -- trending toward normal\n"
                "   - GI symptoms fully resolved\n"
                "   - Iron and ferritin improving\n"
                "   - Vitamin D repleted\n\n"
                "3. Thyroid function normal, anti-TPO negative -- continue annual screening\n\n"
                "4. Autoimmune overlap syndrome well-managed\n"
                "   - T1D + Celiac on shared HLA-DQ2/DQ8 background\n"
                "   - Excellent adherence and self-management"
            ),
            "plan": (
                "1. Continue current insulin pump settings\n"
                "2. Continue strict GFD\n"
                "3. Labs in 3 months: HbA1c, lipid panel\n"
                "4. Recheck tTG-IgA at 12 months (goal: normalization)\n"
                "5. Annual TSH, anti-TPO\n"
                "6. Continue vitamin D 2000 IU, discontinue iron (ferritin normalized)\n"
                "7. Transition care planning: will need adult endocrinologist when ages out of pediatric system\n"
                "8. Follow-up in 3 months"
            ),
        },
    )

    generate_lab_report(
        os.path.join(output_dir, "2025-03-10_Lab_6Month.pdf"),
        PATIENT, "03/10/2025", PROVIDERS["endocrinology"],
        [{
            "panel_name": "DIABETES MONITORING",
            "results": [
                {"test": "HbA1c", "value": "7.4", "unit": "%", "ref_range": "<7.0", "flag": "H"},
                {"test": "C-Peptide (fasting)", "value": "<0.1", "unit": "ng/mL", "ref_range": "0.8-3.1", "flag": "LL"},
                {"test": "Fasting Glucose", "value": "118", "unit": "mg/dL", "ref_range": "70-100", "flag": "H"},
            ],
        }, {
            "panel_name": "CELIAC MONITORING",
            "results": [
                {"test": "tTG-IgA", "value": "18", "unit": "U/mL", "ref_range": "<4.0", "flag": "H"},
                {"test": "Anti-Endomysial Antibody", "value": "Weakly Positive", "unit": "", "ref_range": "Negative", "flag": "A"},
            ],
        }, {
            "panel_name": "NUTRITIONAL AND THYROID",
            "results": [
                {"test": "Iron, Serum", "value": "82", "unit": "mcg/dL", "ref_range": "60-170", "flag": ""},
                {"test": "Ferritin", "value": "38", "unit": "ng/mL", "ref_range": "20-250", "flag": ""},
                {"test": "Vitamin D, 25-OH", "value": "38", "unit": "ng/mL", "ref_range": "30-100", "flag": ""},
                {"test": "TSH", "value": "2.8", "unit": "mIU/L", "ref_range": "0.4-4.0", "flag": ""},
                {"test": "Anti-TPO Antibody", "value": "14", "unit": "IU/mL", "ref_range": "<35", "flag": ""},
            ],
        }],
    )

    count = len([f for f in os.listdir(output_dir) if f.endswith(".pdf")])
    print(f"  Patient James Cooper: {count} documents generated")


if __name__ == "__main__":
    generate(os.path.join(os.path.dirname(__file__), "..", "demo_data", "james_cooper"))
