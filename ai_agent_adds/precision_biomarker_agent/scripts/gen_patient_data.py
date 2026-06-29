#!/usr/bin/env python3
"""Generate sample patient JSON files for the Precision Biomarker Agent."""

import json
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "reference")


def build_patient_1():
    return {
        "id": "HCLS-BIO-2026-00001",
        "demographics": {
            "sex": "Male",
            "age": 45,
            "height_cm": 183,
            "weight_kg": 79.4,
            "bmi": 23.7,
            "ethnicity": "Ashkenazi Jewish",
            "genome": "HG002 (NA24385)"
        },
        "clinical_context": {
            "medications": [
                {"name": "Atorvastatin", "dose": "10mg daily", "indication": "Lipid management (ApoE E4 carrier)", "biomarker_effects": ["LDL", "Total Cholesterol", "CoQ10", "GGT"]},
                {"name": "Lisinopril", "dose": "10mg daily", "indication": "Hypertension", "biomarker_effects": ["Potassium", "Creatinine", "BUN"]},
                {"name": "L-Methylfolate", "dose": "15mg daily", "indication": "MTHFR C677T support", "biomarker_effects": ["Homocysteine", "Folate"]},
                {"name": "Fish Oil", "dose": "2g EPA/DHA daily", "indication": "CV protection", "biomarker_effects": ["Triglycerides", "INR", "Omega-3 Index"]},
                {"name": "Vitamin D3", "dose": "5000 IU daily", "indication": "Vitamin D optimization", "biomarker_effects": ["Vitamin D", "Calcium", "PTH"]},
                {"name": "Methylcobalamin", "dose": "1000 mcg daily", "indication": "Active B12 for MTHFR", "biomarker_effects": ["Vitamin B12", "MMA"]},
                {"name": "CoQ10", "dose": "200mg daily", "indication": "Statin depletion prevention", "biomarker_effects": ["CoQ10"]}
            ],
            "family_history": [
                {"relation": "Father", "condition": "MI at age 58", "genetic_note": "ApoE E4 carrier", "relevance": "Familial hyperlipidemia"},
                {"relation": "Mother", "condition": "Type 2 DM at age 52", "genetic_note": None, "relevance": "Metabolic syndrome"},
                {"relation": "Paternal Grandmother", "condition": "Alzheimer disease at age 74", "genetic_note": "ApoE E4/E4 homozygous (post-mortem)", "relevance": "Multi-generational AD"},
                {"relation": "Maternal Uncle", "condition": "Colorectal cancer at age 61", "genetic_note": None, "relevance": "AJ elevated CRC risk"},
                {"relation": "Paternal Grandfather", "condition": "CAD at age 55", "genetic_note": None, "relevance": "Premature CAD pattern"},
                {"relation": "Sister", "condition": "Gaucher Disease carrier", "genetic_note": "GBA N370S heterozygous", "relevance": "AJ 1/15 carrier rate"},
                {"relation": "Brother", "condition": "No known conditions", "genetic_note": None, "relevance": None}
            ],
            "symptoms": [
                {"description": "Occasional fatigue", "duration": "6 months", "severity": 3, "related_biomarkers": ["Ferritin", "Vitamin D", "Free T3"]},
                {"description": "Mild morning joint stiffness", "duration": "3 months", "severity": 2, "related_biomarkers": ["hs-CRP", "RF", "ESR"]},
                {"description": "Intermittent brain fog", "duration": "4 months", "severity": 2, "related_biomarkers": ["ApoE Genotype", "Homocysteine", "Omega-3 Index"]},
                {"description": "Occasional easy bruising", "duration": "2 months", "severity": 2, "related_biomarkers": ["Platelet Count", "PT/INR", "Fibrinogen"]},
                {"description": "Cold extremities", "duration": "12 months", "severity": 2, "related_biomarkers": ["Free T3", "Ferritin", "Hemoglobin"]}
            ],
            "lifestyle": {
                "diet": "Mediterranean-leaning, AJ traditions",
                "exercise": "3x/week resistance + walking",
                "alcohol": "3-5 glasses red wine/week",
                "sleep_hours": 7.25,
                "smoking": "Never",
                "stress_level": 6,
                "occupation": "Software engineer"
            }
        },
        "genetic_markers": {
            "apoe": "E3/E4",
            "mthfr_c677t": "CT (Heterozygous)",
            "mthfr_a1298c": "AA (Wild-type)",
            "gba": "N370S (Heterozygous carrier)"
        },
        "panels_included": [
            "Annual Membership",
            "Extended Heart & Metabolic",
            "Extended Thyroid",
            "MTHFR Gene",
            "Galleri (GRAIL)",
            "Extended Vitamins/Minerals",
            "Extended Autoimmunity",
            "ApoE Genotype"
        ],
        "total_biomarkers": 257,
        "collection_date": "2026-03-02",
        "fasting": True,
        "lab_partner": "Quest Diagnostics"
    }


def build_patient_2():
    return {
        "id": "HCLS-BIO-2026-00002",
        "demographics": {
            "sex": "Female",
            "age": 38,
            "height_cm": 168,
            "weight_kg": 63.5,
            "bmi": 22.6,
            "ethnicity": "Ashkenazi Jewish",
            "genome": None
        },
        "clinical_context": {
            "medications": [
                {"name": "OCP (ethinyl estradiol/norgestimate)", "dose": "0.035/0.25mg daily", "indication": "Contraception", "biomarker_effects": ["SHBG", "TBG", "Total T4", "Free Testosterone", "Triglycerides", "Fibrinogen"]},
                {"name": "L-Methylfolate", "dose": "15mg daily", "indication": "MTHFR support", "biomarker_effects": ["Homocysteine", "Folate"]},
                {"name": "Vitamin D3", "dose": "2000 IU daily", "indication": "Vitamin D optimization", "biomarker_effects": ["Vitamin D", "Calcium", "PTH"]},
                {"name": "Iron bisglycinate", "dose": "25mg EOD", "indication": "Iron repletion", "biomarker_effects": ["Ferritin", "Iron", "TIBC", "Hemoglobin"]},
                {"name": "Prenatal DHA", "dose": "600mg daily", "indication": "Preconception omega-3 support", "biomarker_effects": ["Omega-3 Index", "DHA"]}
            ],
            "family_history": [
                {"relation": "Mother", "condition": "Breast cancer ER+ at age 48", "genetic_note": "BRCA1 185delAG confirmed", "relevance": "Hereditary breast/ovarian cancer syndrome"},
                {"relation": "Maternal Aunt", "condition": "Ovarian cancer at age 55", "genetic_note": "BRCA1-related", "relevance": "Hereditary breast/ovarian cancer syndrome"},
                {"relation": "Father", "condition": "Type 2 DM at age 58", "genetic_note": None, "relevance": "Metabolic syndrome risk"},
                {"relation": "Paternal Grandmother", "condition": "Gaucher Disease Type 1 at age 35", "genetic_note": "GBA homozygous", "relevance": "AJ Gaucher prevalence"},
                {"relation": "Sister", "condition": "Healthy", "genetic_note": "BRCA1 negative", "relevance": None}
            ],
            "symptoms": [
                {"description": "Heavy menstrual periods", "duration": "2 years", "severity": 4, "related_biomarkers": ["Ferritin", "Iron", "Hemoglobin", "TIBC"]},
                {"description": "Breast tenderness cyclical", "duration": "6 months", "severity": 2, "related_biomarkers": ["Estradiol", "Progesterone", "Prolactin"]},
                {"description": "Anxiety mild", "duration": "6 months", "severity": 3, "related_biomarkers": ["TSH", "Free T4", "Cortisol", "GABA"]},
                {"description": "Fatigue mid-cycle", "duration": "1 year", "severity": 3, "related_biomarkers": ["Ferritin", "Vitamin D", "Free T3", "Hemoglobin"]}
            ],
            "ob_gyn": {
                "menarche": 12,
                "cycle_days": "28-30 on OCP",
                "gravida_para": "G0P0",
                "preconception_planning": "ACTIVE - 12-18 months",
                "brca1_status": "NOT YET TESTED - URGENT",
                "mammography": "Not yet (age <40)"
            },
            "lifestyle": {
                "diet": "Mediterranean-leaning, AJ traditions",
                "exercise": "Yoga 3x/week, walking daily",
                "alcohol": "Rare, social only",
                "sleep_hours": 7.5,
                "smoking": "Never",
                "stress_level": 5,
                "occupation": "Genetic counselor"
            }
        },
        "genetic_markers": {
            "apoe": "Unknown",
            "mthfr_c677t": "Unknown",
            "mthfr_a1298c": "Unknown",
            "gba": "50% carrier risk (Father obligate carrier)"
        },
        "panels_included": [
            "Annual Membership",
            "Extended Heart & Metabolic",
            "Extended Thyroid",
            "MTHFR Gene",
            "Extended Vitamins/Minerals",
            "Extended Autoimmunity",
            "Women's Health Panel"
        ],
        "total_biomarkers": 215,
        "collection_date": "2026-03-02",
        "fasting": True,
        "lab_partner": "Quest Diagnostics"
    }


def build_longitudinal():
    visits = [
        {"visit_number": 1, "date": "2025-03-01", "label": "Baseline"},
        {"visit_number": 2, "date": "2025-09-01", "label": "6-Month Follow-up"},
        {"visit_number": 3, "date": "2026-03-02", "label": "12-Month Follow-up"}
    ]

    # Helper to build a biomarker entry
    def bm(name, unit, loinc, vals, trend, flag_v1, flag_v3, note):
        v1, v3 = vals[0], vals[2]
        pct = round((v3 - v1) / v1 * 100, 1) if v1 != 0 else 0.0
        return {
            "biomarker": name,
            "unit": unit,
            "loinc_code": loinc,
            "values": vals,
            "trend": trend,
            "pct_change_total": pct,
            "flag_v1": flag_v1,
            "flag_v3": flag_v3,
            "clinical_note": note
        }

    tracked = [
        # --- IMPROVING (18) ---
        bm("LDL Cholesterol", "mg/dL", "13457-7", [118, 95, 82], "improving", "Borderline High", "Optimal",
           "Statin + lifestyle approach effective; ApoE E4 carrier target <100 achieved"),
        bm("ApoB", "mg/dL", "1884-6", [95, 80, 72], "improving", "Borderline High", "Optimal",
           "Better discriminator than LDL for E4 carriers; target <80 achieved"),
        bm("Omega-3 Index", "%", "82810-0", [5.8, 7.2, 8.5], "improving", "Suboptimal", "Optimal",
           "Fish oil 2g/day achieving cardioprotective range >8%"),
        bm("AA/EPA Ratio", "ratio", "LP35579-3", [7.1, 4.5, 2.8], "improving", "High", "Optimal",
           "Anti-inflammatory balance improving with omega-3 supplementation"),
        bm("hs-CRP", "mg/L", "30522-7", [0.9, 0.5, 0.3], "improving", "Low Risk", "Optimal",
           "Systemic inflammation well controlled; CV risk minimal"),
        bm("Homocysteine", "umol/L", "2160-0", [9.2, 7.5, 6.5], "improving", "Borderline", "Optimal",
           "L-Methylfolate 15mg normalizing despite MTHFR C677T heterozygous"),
        bm("Vitamin D, 25-OH", "ng/mL", "1989-3", [42, 55, 65], "improving", "Adequate", "Optimal",
           "5000 IU D3 daily achieving functional medicine target 60-80"),
        bm("Insulin, Fasting", "uIU/mL", "1979-3", [7.2, 5.8, 4.5], "improving", "Normal", "Optimal",
           "Insulin sensitivity improving with resistance training and Mediterranean diet"),
        bm("HOMA-IR", "index", "LP73414-6", [1.64, 1.2, 0.95], "improving", "Normal", "Optimal",
           "Insulin resistance resolving; target <1.0 nearly achieved"),
        bm("HbA1c", "%", "4548-4", [5.4, 5.2, 5.1], "improving", "Normal", "Optimal",
           "Glycemic control excellent; maternal T2D risk being mitigated"),
        bm("Non-HDL Cholesterol", "mg/dL", "43396-1", [140, 110, 90], "improving", "Borderline High", "Optimal",
           "Comprehensive atherogenic particle reduction on track"),
        bm("HDL Cholesterol", "mg/dL", "2085-9", [55, 58, 62], "improving", "Normal", "Optimal",
           "HDL rising with exercise and omega-3; cardioprotective trend"),
        bm("LDL Particle Number", "nmol/L", "49132-4", [1180, 950, 780], "improving", "Borderline High", "Optimal",
           "Advanced lipid particle count normalizing; concordant with ApoB improvement"),
        bm("Total Cholesterol", "mg/dL", "2093-3", [195, 178, 165], "improving", "Borderline High", "Desirable",
           "Overall lipid burden reducing with multi-modal intervention"),
        bm("Triglycerides", "mg/dL", "2571-8", [110, 90, 75], "improving", "Normal", "Optimal",
           "Fish oil and carb moderation driving TG reduction"),
        bm("ApoE Protein", "mg/dL", "1871-3", [5.2, 5.0, 4.8], "improving", "Normal", "Normal",
           "ApoE protein stable-low; genotype E3/E4 confirmed separately"),
        bm("Omega-6/3 Ratio", "ratio", "LP35578-5", [4.2, 3.5, 2.8], "improving", "Elevated", "Optimal",
           "Inflammatory balance improving; target <3:1 achieved"),
        bm("Copper/Zinc Ratio", "ratio", "LP73415-3", [1.34, 1.1, 0.95], "improving", "Elevated", "Optimal",
           "Oxidative stress marker normalizing; zinc supplementation effective"),

        # --- CRISIS at V3 (9) ---
        bm("Platelet Count", "x10^3/uL", "777-3", [245, 238, 42], "crisis", "Normal", "CRITICAL LOW",
           "CRITICAL: Severe thrombocytopenia; rule out ITP, TTP, HIT, bone marrow suppression"),
        bm("INR", "ratio", "6301-6", [1.0, 1.0, 4.8], "crisis", "Normal", "CRITICAL HIGH",
           "CRITICAL: Severely elevated INR without anticoagulation; evaluate hepatic synthetic function, DIC"),
        bm("Glucose, Fasting", "mg/dL", "1558-6", [92, 88, 285], "crisis", "Normal", "CRITICAL HIGH",
           "CRITICAL: New-onset hyperglycemia; rule out DKA, steroid-induced, pancreatic pathology"),
        bm("Free T3", "pg/mL", "3051-0", [3.1, 3.2, 5.8], "crisis", "Normal", "CRITICAL HIGH",
           "CRITICAL: Thyrotoxicosis level; rule out thyroid storm, exogenous T3, Graves disease"),
        bm("Fibrinogen", "mg/dL", "3255-7", [285, 270, 480], "crisis", "Normal", "CRITICAL HIGH",
           "CRITICAL: Acute phase reactant markedly elevated; correlate with DIC panel, infection"),
        bm("Ferritin", "ng/mL", "2276-4", [125, 115, 12], "crisis", "Normal", "CRITICAL LOW",
           "CRITICAL: Severe iron depletion; evaluate acute blood loss, malabsorption, occult GI bleeding"),
        bm("Cystatin C", "mg/L", "33863-2", [0.82, 0.80, 1.15], "crisis", "Normal", "HIGH",
           "CRITICAL: Acute kidney injury marker; more reliable than creatinine; evaluate nephrotoxic exposure"),
        bm("GGT", "U/L", "2324-2", [22, 28, 85], "crisis", "Normal", "HIGH",
           "CRITICAL: Hepatobiliary injury; evaluate drug-induced (statin), alcohol, biliary obstruction"),
        bm("Potassium", "mEq/L", "2823-3", [4.3, 4.2, 6.2], "crisis", "Normal", "CRITICAL HIGH",
           "CRITICAL: Life-threatening hyperkalemia; evaluate renal function, ACEi effect, hemolysis artifact"),

        # --- STABLE (19) ---
        bm("WBC", "x10^3/uL", "6690-2", [5.8, 5.5, 5.8], "stable", "Normal", "Normal",
           "White cell count stable within reference range"),
        bm("RBC", "x10^6/uL", "789-8", [4.95, 4.90, 4.95], "stable", "Normal", "Normal",
           "Red cell mass stable; no evidence of polycythemia or anemia"),
        bm("Hemoglobin", "g/dL", "718-7", [15.2, 15.0, 15.2], "stable", "Normal", "Normal",
           "Hemoglobin stable despite ferritin crisis; may reflect early iron depletion before anemia"),
        bm("Hematocrit", "%", "4544-3", [44.8, 44.5, 44.8], "stable", "Normal", "Normal",
           "Hematocrit tracking with hemoglobin; stable"),
        bm("TSH", "mIU/L", "3016-3", [2.1, 2.0, 2.1], "stable", "Normal", "Normal",
           "TSH paradoxically normal despite elevated Free T3; evaluate for TSH-secreting adenoma"),
        bm("Free T4", "ng/dL", "3024-7", [1.25, 1.22, 1.25], "stable", "Normal", "Normal",
           "Free T4 normal despite Free T3 elevation; possible T3 thyrotoxicosis pattern"),
        bm("Sodium", "mEq/L", "2951-2", [140, 141, 140], "stable", "Normal", "Normal",
           "Sodium stable; no dilutional or depletional abnormalities"),
        bm("Creatinine", "mg/dL", "2160-0", [1.02, 1.00, 1.02], "stable", "Normal", "Normal",
           "Creatinine stable but may lag behind cystatin C in detecting AKI"),
        bm("eGFR", "mL/min/1.73m2", "33914-3", [95, 96, 95], "stable", "Normal", "Normal",
           "eGFR stable but cystatin C-based eGFR may show decline; recommend cystatin C-based calculation"),
        bm("Albumin", "g/dL", "1751-7", [4.5, 4.5, 4.5], "stable", "Normal", "Normal",
           "Albumin stable; hepatic synthetic function preserved despite GGT elevation"),
        bm("Calcium", "mg/dL", "17861-6", [9.6, 9.5, 9.6], "stable", "Normal", "Normal",
           "Calcium stable on vitamin D supplementation; no hypercalcemia"),
        bm("Total Protein", "g/dL", "2885-2", [7.1, 7.0, 7.1], "stable", "Normal", "Normal",
           "Total protein stable; no evidence of paraprotein or depletion"),
        bm("AST", "U/L", "1920-8", [24, 22, 24], "stable", "Normal", "Normal",
           "AST stable; hepatocellular function preserved despite GGT rise"),
        bm("ALT", "U/L", "1742-6", [28, 26, 28], "stable", "Normal", "Normal",
           "ALT stable; isolated GGT elevation suggests biliary rather than hepatocellular process"),
        bm("Magnesium, RBC", "mg/dL", "19123-9", [5.2, 5.3, 5.2], "stable", "Normal", "Normal",
           "Intracellular magnesium stable; superior to serum magnesium for status assessment"),
        bm("Zinc", "ug/dL", "2601-3", [82, 85, 82], "stable", "Normal", "Normal",
           "Zinc stable; copper/zinc ratio improvement driven by copper reduction"),
        bm("Vitamin B12", "pg/mL", "2132-9", [620, 650, 620], "stable", "Normal", "Normal",
           "B12 well-maintained on methylcobalamin 1000 mcg; active form for MTHFR support"),
        bm("Folate", "ng/mL", "2284-8", [15.8, 16.5, 15.8], "stable", "Normal", "Normal",
           "Folate optimal on L-methylfolate 15mg; bypassing MTHFR enzyme deficiency"),
        bm("Selenium", "ug/L", "2823-3", [135, 138, 135], "stable", "Normal", "Normal",
           "Selenium adequate; supports glutathione peroxidase and thyroid function"),
    ]

    # Genetic markers (fixed values across visits)
    genetic = [
        {
            "biomarker": "ApoE Genotype",
            "unit": "genotype",
            "loinc_code": "33640-4",
            "values": ["E3/E4", "E3/E4", "E3/E4"],
            "trend": "stable",
            "pct_change_total": None,
            "flag_v1": "Risk Allele Present",
            "flag_v3": "Risk Allele Present",
            "clinical_note": "E4 allele confirmed; increased AD and CVD risk; guides statin and lifestyle intensity"
        },
        {
            "biomarker": "MTHFR C677T",
            "unit": "genotype",
            "loinc_code": "53576-8",
            "values": ["CT", "CT", "CT"],
            "trend": "stable",
            "pct_change_total": None,
            "flag_v1": "Heterozygous",
            "flag_v3": "Heterozygous",
            "clinical_note": "Heterozygous C677T; ~35% reduced enzyme activity; methylfolate supplementation indicated"
        },
        {
            "biomarker": "GBA N370S",
            "unit": "genotype",
            "loinc_code": "55197-1",
            "values": [None, None, "Heterozygous"],
            "trend": "stable",
            "pct_change_total": None,
            "flag_v1": "Not Tested",
            "flag_v3": "Carrier Identified",
            "clinical_note": "NEW at V3: GBA N370S carrier confirmed; AJ 1/15 carrier frequency; sister also carrier; Gaucher and Parkinson risk counseling indicated"
        }
    ]

    tracked.extend(genetic)

    return {
        "patient_id": "HCLS-BIO-2026-00001",
        "schema_version": "1.0",
        "visits": visits,
        "tracked_biomarkers": tracked
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # File 1: Sample patients
    patients = [build_patient_1(), build_patient_2()]
    patients_path = os.path.join(OUTPUT_DIR, "biomarker_sample_patients.json")
    with open(patients_path, "w") as f:
        json.dump(patients, f, indent=2)
    print(f"Wrote {patients_path} ({len(patients)} patients)")

    # File 2: Longitudinal tracking
    longitudinal = build_longitudinal()
    longitudinal_path = os.path.join(OUTPUT_DIR, "biomarker_longitudinal_tracking.json")
    with open(longitudinal_path, "w") as f:
        json.dump(longitudinal, f, indent=2)
    n_biomarkers = len(longitudinal["tracked_biomarkers"])
    print(f"Wrote {longitudinal_path} ({n_biomarkers} tracked biomarkers)")


if __name__ == "__main__":
    main()
