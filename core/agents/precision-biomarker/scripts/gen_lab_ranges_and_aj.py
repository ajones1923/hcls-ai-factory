#!/usr/bin/env python3
"""
Generate two reference JSON files:
  1. biomarker_lab_ranges.json   -- multi-lab reference range comparison
  2. biomarker_aj_carrier_screening.json -- Ashkenazi Jewish carrier screening panel
"""

import json
import os
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "reference"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
def _r(min_val, max_val, unit):
    """Return a range dict; use None for open-ended bounds."""
    return {"min": min_val, "max": max_val, "unit": unit}


# ---------------------------------------------------------------------------
# FILE 1 -- biomarker_lab_ranges.json
# ---------------------------------------------------------------------------
def build_lab_ranges():
    # fmt: off
    ranges = {
        "glucose_fasting": {
            "name": "Glucose, Fasting",
            "quest":    _r(65, 99, "mg/dL"),
            "labcorp":  _r(70, 100, "mg/dL"),
            "optimal":  _r(72, 90, "mg/dL"),
        },
        "hemoglobin_a1c": {
            "name": "Hemoglobin A1c (HbA1c)",
            "quest":    _r(4.8, 5.6, "%"),
            "labcorp":  _r(None, 5.7, "%"),
            "optimal":  _r(None, 5.2, "%"),
        },
        "insulin_fasting": {
            "name": "Insulin, Fasting",
            "quest":    _r(2.6, 24.9, "uIU/mL"),
            "labcorp":  _r(2.6, 24.9, "uIU/mL"),
            "optimal":  _r(2.0, 5.0, "uIU/mL"),
        },
        "homa_ir": {
            "name": "HOMA-IR",
            "quest":    _r(None, 2.5, "index"),
            "labcorp":  _r(None, 2.5, "index"),
            "optimal":  _r(None, 1.0, "index"),
        },
        "ferritin": {
            "name": "Ferritin",
            "quest":    _r(30, 400, "ng/mL"),
            "labcorp":  _r(22, 322, "ng/mL"),
            "optimal":  _r(40, 200, "ng/mL"),
        },
        "transferrin_saturation": {
            "name": "Transferrin Saturation",
            "quest":    _r(15, 55, "%"),
            "labcorp":  _r(15, 55, "%"),
            "optimal":  _r(20, 45, "%"),
        },
        "iron_total": {
            "name": "Iron, Total",
            "quest":    _r(38, 169, "ug/dL"),
            "labcorp":  _r(38, 169, "ug/dL"),
            "optimal":  _r(60, 150, "ug/dL"),
        },
        "tibc": {
            "name": "Total Iron Binding Capacity (TIBC)",
            "quest":    _r(250, 370, "ug/dL"),
            "labcorp":  _r(250, 370, "ug/dL"),
            "optimal":  _r(250, 350, "ug/dL"),
        },
        "free_t3": {
            "name": "Free T3",
            "quest":    _r(2.0, 4.4, "pg/mL"),
            "labcorp":  _r(2.0, 4.4, "pg/mL"),
            "optimal":  _r(3.0, 4.0, "pg/mL"),
        },
        "free_t4": {
            "name": "Free T4",
            "quest":    _r(0.82, 1.77, "ng/dL"),
            "labcorp":  _r(0.82, 1.77, "ng/dL"),
            "optimal":  _r(1.0, 1.5, "ng/dL"),
        },
        "tsh": {
            "name": "TSH",
            "quest":    _r(0.45, 4.5, "uIU/mL"),
            "labcorp":  _r(0.45, 4.5, "uIU/mL"),
            "optimal":  _r(1.0, 2.5, "uIU/mL"),
        },
        "ldl_cholesterol": {
            "name": "LDL Cholesterol",
            "quest":    _r(None, 99, "mg/dL"),
            "labcorp":  _r(None, 100, "mg/dL"),
            "optimal":  _r(None, 70, "mg/dL"),
        },
        "hdl_cholesterol": {
            "name": "HDL Cholesterol",
            "quest":    _r(39, None, "mg/dL"),
            "labcorp":  _r(40, None, "mg/dL"),
            "optimal":  _r(60, None, "mg/dL"),
        },
        "total_cholesterol": {
            "name": "Total Cholesterol",
            "quest":    _r(100, 199, "mg/dL"),
            "labcorp":  _r(100, 199, "mg/dL"),
            "optimal":  _r(100, 180, "mg/dL"),
        },
        "triglycerides": {
            "name": "Triglycerides",
            "quest":    _r(None, 149, "mg/dL"),
            "labcorp":  _r(None, 150, "mg/dL"),
            "optimal":  _r(None, 80, "mg/dL"),
        },
        "apolipoprotein_b": {
            "name": "Apolipoprotein B (ApoB)",
            "quest":    _r(None, 90, "mg/dL"),
            "labcorp":  _r(None, 90, "mg/dL"),
            "optimal":  _r(None, 80, "mg/dL"),
        },
        "lipoprotein_a": {
            "name": "Lipoprotein(a)",
            "quest":    _r(None, 75, "nmol/L"),
            "labcorp":  _r(None, 75, "nmol/L"),
            "optimal":  _r(None, 30, "nmol/L"),
        },
        "hs_crp": {
            "name": "High-Sensitivity C-Reactive Protein (hs-CRP)",
            "quest":    _r(None, 3.0, "mg/L"),
            "labcorp":  _r(None, 3.0, "mg/L"),
            "optimal":  _r(None, 0.5, "mg/L"),
        },
        "homocysteine": {
            "name": "Homocysteine",
            "quest":    _r(None, 10.4, "umol/L"),
            "labcorp":  _r(None, 11.4, "umol/L"),
            "optimal":  _r(None, 7.0, "umol/L"),
        },
        "cystatin_c": {
            "name": "Cystatin C",
            "quest":    _r(0.53, 0.95, "mg/L"),
            "labcorp":  _r(0.53, 0.95, "mg/L"),
            "optimal":  _r(0.55, 0.85, "mg/L"),
        },
        "vitamin_d_25oh": {
            "name": "Vitamin D, 25-Hydroxy",
            "quest":    _r(30, 100, "ng/mL"),
            "labcorp":  _r(30, 100, "ng/mL"),
            "optimal":  _r(50, 80, "ng/mL"),
        },
        "vitamin_b12": {
            "name": "Vitamin B12",
            "quest":    _r(232, 1245, "pg/mL"),
            "labcorp":  _r(232, 1245, "pg/mL"),
            "optimal":  _r(500, 1000, "pg/mL"),
        },
        "folate": {
            "name": "Folate",
            "quest":    _r(2.7, 17.0, "ng/mL"),
            "labcorp":  _r(2.7, 17.0, "ng/mL"),
            "optimal":  _r(10.0, 25.0, "ng/mL"),
        },
        "omega_3_index": {
            "name": "Omega-3 Index",
            "quest":    _r(4.0, None, "%"),
            "labcorp":  _r(4.0, None, "%"),
            "optimal":  _r(8.0, None, "%"),
        },
        "magnesium_rbc": {
            "name": "Magnesium, RBC",
            "quest":    _r(4.2, 6.8, "mg/dL"),
            "labcorp":  _r(4.2, 6.8, "mg/dL"),
            "optimal":  _r(5.0, 6.5, "mg/dL"),
        },
        "ggt": {
            "name": "Gamma-Glutamyl Transferase (GGT)",
            "quest":    _r(None, 65, "U/L"),
            "labcorp":  _r(None, 65, "U/L"),
            "optimal":  _r(10, 40, "U/L"),
        },
        "alt": {
            "name": "Alanine Aminotransferase (ALT)",
            "quest":    _r(None, 44, "U/L"),
            "labcorp":  _r(None, 44, "U/L"),
            "optimal":  _r(None, 25, "U/L"),
        },
        "ast": {
            "name": "Aspartate Aminotransferase (AST)",
            "quest":    _r(None, 40, "U/L"),
            "labcorp":  _r(None, 40, "U/L"),
            "optimal":  _r(None, 25, "U/L"),
        },
        "uric_acid": {
            "name": "Uric Acid",
            "quest":    _r(3.7, 8.6, "mg/dL"),
            "labcorp":  _r(3.7, 8.6, "mg/dL"),
            "optimal":  _r(3.5, 5.5, "mg/dL"),
        },
        "creatinine": {
            "name": "Creatinine",
            "quest":    _r(0.76, 1.27, "mg/dL"),
            "labcorp":  _r(0.76, 1.27, "mg/dL"),
            "optimal":  _r(0.7, 1.1, "mg/dL"),
        },
        "bun": {
            "name": "Blood Urea Nitrogen (BUN)",
            "quest":    _r(6, 24, "mg/dL"),
            "labcorp":  _r(6, 24, "mg/dL"),
            "optimal":  _r(10, 20, "mg/dL"),
        },
        "egfr": {
            "name": "Estimated GFR",
            "quest":    _r(60, None, "mL/min/1.73m2"),
            "labcorp":  _r(60, None, "mL/min/1.73m2"),
            "optimal":  _r(90, None, "mL/min/1.73m2"),
        },
        "albumin": {
            "name": "Albumin",
            "quest":    _r(3.5, 5.5, "g/dL"),
            "labcorp":  _r(3.5, 5.5, "g/dL"),
            "optimal":  _r(4.2, 5.0, "g/dL"),
        },
        "total_protein": {
            "name": "Total Protein",
            "quest":    _r(6.0, 8.5, "g/dL"),
            "labcorp":  _r(6.0, 8.5, "g/dL"),
            "optimal":  _r(6.5, 7.5, "g/dL"),
        },
        "white_blood_cells": {
            "name": "White Blood Cells (WBC)",
            "quest":    _r(3.4, 10.8, "x10E3/uL"),
            "labcorp":  _r(3.4, 10.8, "x10E3/uL"),
            "optimal":  _r(4.0, 7.0, "x10E3/uL"),
        },
        "hemoglobin": {
            "name": "Hemoglobin",
            "quest":    _r(12.6, 17.7, "g/dL"),
            "labcorp":  _r(12.6, 17.7, "g/dL"),
            "optimal":  _r(13.5, 16.0, "g/dL"),
        },
        "platelets": {
            "name": "Platelets",
            "quest":    _r(150, 379, "x10E3/uL"),
            "labcorp":  _r(150, 379, "x10E3/uL"),
            "optimal":  _r(175, 300, "x10E3/uL"),
        },
        "dhea_s": {
            "name": "DHEA-Sulfate",
            "quest":    _r(44, 332, "ug/dL"),
            "labcorp":  _r(44, 332, "ug/dL"),
            "optimal":  _r(150, 300, "ug/dL"),
        },
        "cortisol_am": {
            "name": "Cortisol, AM",
            "quest":    _r(6.2, 19.4, "ug/dL"),
            "labcorp":  _r(6.2, 19.4, "ug/dL"),
            "optimal":  _r(10, 15, "ug/dL"),
        },
        "testosterone_total": {
            "name": "Testosterone, Total",
            "quest":    _r(264, 916, "ng/dL"),
            "labcorp":  _r(264, 916, "ng/dL"),
            "optimal":  _r(500, 900, "ng/dL"),
        },
        "estradiol": {
            "name": "Estradiol",
            "quest":    _r(7.6, 42.6, "pg/mL"),
            "labcorp":  _r(7.6, 42.6, "pg/mL"),
            "optimal":  _r(20, 35, "pg/mL"),
        },
    }
    # fmt: on

    # Restructure into per-lab format
    quest_ranges = {}
    labcorp_ranges = {}
    optimal_ranges = {}

    for key, entry in ranges.items():
        quest_ranges[key] = {
            "name": entry["name"],
            **entry["quest"],
        }
        labcorp_ranges[key] = {
            "name": entry["name"],
            **entry["labcorp"],
        }
        optimal_ranges[key] = {
            "name": entry["name"],
            **entry["optimal"],
        }

    return {
        "schema_version": "1.0",
        "description": "Multi-lab reference range comparison for biomarker interpretation across lab providers",
        "labs": {
            "quest_diagnostics": {
                "name": "Quest Diagnostics",
                "ranges": quest_ranges,
            },
            "labcorp": {
                "name": "LabCorp",
                "ranges": labcorp_ranges,
            },
            "function_health_optimal": {
                "name": "Function Health Optimal",
                "ranges": optimal_ranges,
            },
        },
    }


# ---------------------------------------------------------------------------
# FILE 2 -- biomarker_aj_carrier_screening.json
# ---------------------------------------------------------------------------
def build_aj_carrier_screening():
    return [
        {
            "id": "ajcs_brca1",
            "gene": "BRCA1",
            "disease": "Hereditary Breast and Ovarian Cancer Syndrome",
            "inheritance": "autosomal_dominant",
            "aj_carrier_frequency": "1 in 40",
            "general_carrier_frequency": "1 in 400",
            "common_mutations": ["185delAG"],
            "loinc_code": "55233-2",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Pathogenic variants in BRCA1 confer a 55-72% lifetime risk of breast cancer and 39-44% lifetime risk of ovarian cancer. Male carriers have increased risks of breast and prostate cancer.",
            "reproductive_implications": "Carriers should be offered genetic counseling. Family planning may include preimplantation genetic testing (PGT) to avoid transmitting the variant. Children of a carrier have a 50% chance of inheriting the mutation.",
            "compound_risks": [
                "Combined BRCA1 + CHEK2 variants may further elevate breast cancer risk",
                "Interaction with hormone replacement therapy increases risk profile"
            ],
            "text_chunk": "BRCA1 185delAG is one of three founder mutations highly prevalent in Ashkenazi Jewish populations. The carrier frequency is approximately 1 in 40 among AJ individuals compared to 1 in 400 in the general population. Pathogenic variants confer a 55-72% lifetime risk of breast cancer and a 39-44% lifetime risk of ovarian cancer. PARP inhibitors such as olaparib have demonstrated significant efficacy in BRCA1-mutated cancers. Screening recommendations include annual breast MRI starting at age 25 and consideration of risk-reducing salpingo-oophorectomy by age 35-40."
        },
        {
            "id": "ajcs_brca2",
            "gene": "BRCA2",
            "disease": "Hereditary Breast, Ovarian, and Prostate Cancer Syndrome",
            "inheritance": "autosomal_dominant",
            "aj_carrier_frequency": "1 in 80",
            "general_carrier_frequency": "1 in 400",
            "common_mutations": ["6174delT"],
            "loinc_code": "55234-0",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Pathogenic variants in BRCA2 confer a 45-69% lifetime risk of breast cancer, 11-17% risk of ovarian cancer, and significantly elevated prostate cancer risk in male carriers (up to 20% lifetime risk).",
            "reproductive_implications": "Carriers should receive genetic counseling prior to family planning. Preimplantation genetic testing is available. Each child of a carrier has a 50% chance of inheriting the variant.",
            "compound_risks": [
                "BRCA2 + PALB2 variants may further elevate pancreatic cancer risk",
                "Male carriers face elevated prostate and breast cancer risk"
            ],
            "text_chunk": "BRCA2 6174delT is an Ashkenazi Jewish founder mutation with a carrier frequency of approximately 1 in 80. Unlike BRCA1, BRCA2 mutations carry a somewhat lower ovarian cancer risk (11-17%) but still confer a 45-69% lifetime breast cancer risk. Male carriers face notably increased prostate cancer risk (up to 20% lifetime). BRCA2-mutated cancers also respond to PARP inhibitor therapy. Pancreatic cancer risk is also elevated. Screening follows similar protocols to BRCA1 with annual breast MRI and consideration of risk-reducing surgeries."
        },
        {
            "id": "ajcs_gba",
            "gene": "GBA",
            "disease": "Gaucher Disease",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 15",
            "general_carrier_frequency": "1 in 100",
            "common_mutations": ["N370S", "L444P", "84GG", "IVS2+1"],
            "loinc_code": "49874-2",
            "method": "PCR with targeted mutation analysis and enzyme assay",
            "clinical_significance": "Gaucher Disease Type 1 (non-neuronopathic) is the most common form in AJ populations. It causes hepatosplenomegaly, anemia, thrombocytopenia, and bone disease. Enzyme replacement therapy (ERT) and substrate reduction therapy (SRT) are available treatments.",
            "reproductive_implications": "With a carrier frequency of 1 in 15, AJ couples have approximately a 1 in 900 chance of having an affected child. Carrier testing and genetic counseling are strongly recommended before conception.",
            "compound_risks": [
                "GBA carriers (heterozygotes) have a 5-fold increased risk of Parkinson's disease",
                "Compound risk with APOE E4 for Parkinson's disease and Lewy body dementia",
                "GBA + LRRK2 variants may further elevate Parkinson's risk"
            ],
            "text_chunk": "GBA mutations are the most common genetic risk factor for Gaucher Disease in Ashkenazi Jewish populations, with a carrier frequency of approximately 1 in 15. The N370S mutation is the most prevalent, associated with Type 1 (non-neuronopathic) Gaucher Disease. Importantly, even heterozygous GBA carriers face a 5-fold increased risk of developing Parkinson's disease. The compound interaction between GBA carrier status and APOE E4 significantly elevates the risk of Parkinson's disease and Lewy body dementia. Enzyme replacement therapy with imiglucerase or velaglucerase alfa is effective for symptomatic Gaucher Disease. Substrate reduction therapy with eliglustat is an oral alternative."
        },
        {
            "id": "ajcs_hexa",
            "gene": "HEXA",
            "disease": "Tay-Sachs Disease",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 30",
            "general_carrier_frequency": "1 in 300",
            "common_mutations": ["1278insTATC", "IVS12+1G>C"],
            "loinc_code": "49875-9",
            "method": "PCR with targeted mutation analysis and hexosaminidase A enzyme assay",
            "clinical_significance": "Tay-Sachs Disease is a fatal neurodegenerative disorder caused by deficiency of hexosaminidase A enzyme, leading to accumulation of GM2 ganglioside in neurons. Infantile-onset form causes progressive neurodegeneration, seizures, blindness, and death typically by age 4-5.",
            "reproductive_implications": "Carrier screening is one of the most successful genetic screening programs in history, having reduced Tay-Sachs incidence by over 90% in AJ populations. Couples where both partners are carriers have a 25% chance of an affected child per pregnancy.",
            "compound_risks": [
                "Late-onset Tay-Sachs may present with psychiatric symptoms and motor neuron disease"
            ],
            "text_chunk": "HEXA mutations causing Tay-Sachs Disease have a carrier frequency of approximately 1 in 30 among Ashkenazi Jews. The 1278insTATC insertion and IVS12+1G>C splice-site mutation are the two most common AJ founder mutations. The carrier screening program for Tay-Sachs, initiated in 1971, is considered one of the most successful genetic disease prevention programs in history, reducing disease incidence by over 90%. Infantile Tay-Sachs presents with developmental regression at 3-6 months, followed by progressive neurodegeneration. Late-onset forms exist with milder phenotypes including psychiatric manifestations and motor neuron disease."
        },
        {
            "id": "ajcs_fancc",
            "gene": "FANCC",
            "disease": "Fanconi Anemia Type C",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 89",
            "general_carrier_frequency": "1 in 500",
            "common_mutations": ["IVS4+4A>T"],
            "loinc_code": "49876-7",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Fanconi Anemia Type C causes bone marrow failure, congenital malformations (short stature, abnormal thumbs, skin pigmentation), and markedly increased cancer susceptibility, particularly acute myeloid leukemia and squamous cell carcinomas.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child per pregnancy. Prenatal diagnosis and preimplantation genetic testing are available. Genetic counseling is recommended for all AJ couples.",
            "compound_risks": [
                "FANCC-associated bone marrow failure may require hematopoietic stem cell transplant",
                "Elevated risk of head and neck squamous cell carcinomas and gynecologic cancers"
            ],
            "text_chunk": "FANCC mutations causing Fanconi Anemia Type C have a carrier frequency of approximately 1 in 89 in Ashkenazi Jewish populations. The IVS4+4A>T splice-site mutation is the predominant AJ founder mutation. Fanconi Anemia causes progressive bone marrow failure typically presenting in the first decade of life, along with congenital anomalies and a dramatically elevated cancer risk. Hematopoietic stem cell transplant is the definitive treatment for bone marrow failure. Affected individuals require lifelong cancer surveillance due to a 50-fold increased risk of acute myeloid leukemia and elevated risk of solid tumors."
        },
        {
            "id": "ajcs_aspa",
            "gene": "ASPA",
            "disease": "Canavan Disease",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 40",
            "general_carrier_frequency": "rare",
            "common_mutations": ["E285A", "Y231X"],
            "loinc_code": "49877-5",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Canavan Disease is a severe leukodystrophy caused by aspartoacylase deficiency, leading to spongy degeneration of the white matter. It presents with macrocephaly, hypotonia, and progressive loss of motor and cognitive milestones, typically fatal in early childhood.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child. Prenatal diagnosis via N-acetylaspartic acid levels in amniotic fluid or molecular testing is available. Carrier screening is recommended for all AJ individuals.",
            "compound_risks": [
                "No significant compound genetic risks identified; disease severity is primarily driven by residual enzyme activity"
            ],
            "text_chunk": "ASPA mutations causing Canavan Disease have a carrier frequency of approximately 1 in 40 among Ashkenazi Jews. The E285A missense mutation accounts for approximately 84% of AJ carrier alleles, while Y231X accounts for most of the remainder. Canavan Disease presents in infancy with macrocephaly, severe hypotonia, and failure to achieve motor milestones. MRI shows diffuse white matter involvement with elevated N-acetylaspartic acid on MR spectroscopy. There is currently no cure, though gene therapy trials are ongoing. Supportive care focuses on seizure management, nutrition, and physical therapy."
        },
        {
            "id": "ajcs_blm",
            "gene": "BLM",
            "disease": "Bloom Syndrome",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 100",
            "general_carrier_frequency": "rare",
            "common_mutations": ["blmAsh"],
            "loinc_code": "49878-3",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Bloom Syndrome is characterized by proportional pre- and postnatal growth deficiency, sun-sensitive facial erythema, immunodeficiency, and a dramatically elevated risk of a wide range of cancers occurring at early ages. Mean age of cancer diagnosis is approximately 25 years.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child per pregnancy. The blmAsh founder mutation accounts for nearly all AJ carrier alleles, making targeted screening highly effective.",
            "compound_risks": [
                "Elevated risk of multiple cancer types including leukemia, lymphoma, and carcinomas at young ages",
                "Immunodeficiency increases susceptibility to infections"
            ],
            "text_chunk": "BLM mutations causing Bloom Syndrome have a carrier frequency of approximately 1 in 100 in Ashkenazi Jewish populations. The blmAsh (2281del6ins7) founder mutation accounts for virtually all AJ carrier alleles. Bloom Syndrome is caused by deficiency of the BLM RecQ helicase, leading to genomic instability and a dramatically elevated cancer predisposition. Affected individuals develop proportional dwarfism, a characteristic sun-sensitive facial rash, and immunodeficiency. The cancer risk is striking, with approximately 50% of affected individuals developing cancer by age 25. The Bloom Syndrome Registry has documented over 300 cancers in approximately 270 registered individuals."
        },
        {
            "id": "ajcs_smpd1",
            "gene": "SMPD1",
            "disease": "Niemann-Pick Disease Type A",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 90",
            "general_carrier_frequency": "rare",
            "common_mutations": ["R496L", "L302P"],
            "loinc_code": "49879-1",
            "method": "PCR with targeted mutation analysis and acid sphingomyelinase enzyme assay",
            "clinical_significance": "Niemann-Pick Disease Type A is a fatal neurovisceral lipid storage disorder caused by acid sphingomyelinase deficiency. It presents with hepatosplenomegaly in infancy, progressive neurodegeneration, and a characteristic cherry-red macular spot. Death typically occurs by age 3.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child per pregnancy. Enzyme assay and molecular testing can distinguish Type A (severe) from Type B (milder) SMPD1 mutations. Prenatal and preimplantation genetic diagnosis are available.",
            "compound_risks": [
                "SMPD1 Type B mutations cause a milder chronic visceral form without neurodegeneration",
                "Compound heterozygosity for Type A and Type B alleles results in intermediate phenotypes"
            ],
            "text_chunk": "SMPD1 mutations causing Niemann-Pick Disease Type A have a carrier frequency of approximately 1 in 90 among Ashkenazi Jews. The R496L and L302P mutations are the most common AJ alleles. Type A Niemann-Pick is a devastating infantile-onset neurodegenerative disorder caused by severe deficiency of acid sphingomyelinase. Affected infants develop hepatosplenomegaly within the first months of life, followed by progressive neurological deterioration. The characteristic cherry-red macular spot is present in about 50% of cases. There is currently no effective treatment, and death typically occurs by age 2-3 years."
        },
        {
            "id": "ajcs_ikbkap",
            "gene": "IKBKAP/ELP1",
            "disease": "Familial Dysautonomia (Riley-Day Syndrome)",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 30",
            "general_carrier_frequency": "extremely rare",
            "common_mutations": ["IVS20+6T>C"],
            "loinc_code": "49880-9",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Familial Dysautonomia is a severe disorder of the autonomic and sensory nervous system. It causes episodes of autonomic crises (hypertension, tachycardia, vomiting), insensitivity to pain, absent overflow tears, labile blood pressure, and progressive neurodegeneration.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child per pregnancy. The single founder mutation (IVS20+6T>C) accounts for over 99.5% of AJ disease alleles, making carrier screening highly sensitive in this population.",
            "compound_risks": [
                "Progressive renal disease is a significant long-term complication",
                "Autonomic crises can be life-threatening and require specialized management"
            ],
            "text_chunk": "IKBKAP (ELP1) mutations causing Familial Dysautonomia have a carrier frequency of approximately 1 in 30 in Ashkenazi Jewish populations. The IVS20+6T>C splice-site mutation accounts for over 99.5% of AJ disease alleles, causing tissue-specific exon skipping and reduced IKAP protein levels. Familial Dysautonomia affects the development and survival of sensory and autonomic neurons. Clinical features include absence of overflow tears, labile blood pressure, episodic vomiting crises, insensitivity to pain, and progressive gait ataxia. Median survival has improved to approximately 40 years with specialized care. Management includes prevention of aspiration pneumonia, blood pressure management, and treatment of autonomic crises."
        },
        {
            "id": "ajcs_mcoln1",
            "gene": "MCOLN1",
            "disease": "Mucolipidosis Type IV",
            "inheritance": "autosomal_recessive",
            "aj_carrier_frequency": "1 in 100",
            "general_carrier_frequency": "rare",
            "common_mutations": ["IVS3-2A>G"],
            "loinc_code": "49881-7",
            "method": "PCR with targeted mutation analysis",
            "clinical_significance": "Mucolipidosis IV is a lysosomal storage disorder caused by mucolipin-1 deficiency. It presents with psychomotor retardation evident by the first year of life, progressive visual impairment from corneal clouding and retinal degeneration, and achlorhydria with iron-deficiency anemia.",
            "reproductive_implications": "Carrier couples have a 25% chance of an affected child per pregnancy. The AJ founder mutation IVS3-2A>G accounts for approximately 72% of disease alleles in this population. Prenatal diagnosis is available through molecular testing.",
            "compound_risks": [
                "Iron-deficiency anemia from achlorhydria requires monitoring and supplementation",
                "Progressive visual loss may lead to blindness in adolescence"
            ],
            "text_chunk": "MCOLN1 mutations causing Mucolipidosis Type IV have a carrier frequency of approximately 1 in 100 in Ashkenazi Jewish populations. The IVS3-2A>G splice-site mutation is the most common AJ founder mutation, accounting for approximately 72% of disease alleles. Mucolipidosis IV is caused by deficiency of the mucolipin-1 channel protein, leading to impaired lysosomal calcium signaling and accumulation of lipids and mucopolysaccharides. Most affected individuals achieve limited motor milestones but do not progress beyond an early developmental level. Corneal clouding and retinal degeneration cause progressive visual impairment. Achlorhydria leads to constitutive iron-deficiency anemia requiring supplementation."
        },
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    # File 1
    lab_ranges = build_lab_ranges()
    path1 = DATA_DIR / "biomarker_lab_ranges.json"
    with open(path1, "w") as f:
        json.dump(lab_ranges, f, indent=2)
    n_biomarkers = len(next(iter(lab_ranges["labs"].values()))["ranges"])
    print(f"Wrote {path1}  ({n_biomarkers} biomarkers, 3 labs)")

    # File 2
    aj_data = build_aj_carrier_screening()
    path2 = DATA_DIR / "biomarker_aj_carrier_screening.json"
    with open(path2, "w") as f:
        json.dump(aj_data, f, indent=2)
    print(f"Wrote {path2}  ({len(aj_data)} genes)")


if __name__ == "__main__":
    main()
