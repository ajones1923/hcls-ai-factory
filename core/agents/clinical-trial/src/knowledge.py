"""Clinical Trial Intelligence Agent — Domain Knowledge Base.

Comprehensive clinical trial knowledge covering therapeutic areas, trial phases,
regulatory agencies, endpoint types, adaptive designs, biomarker strategies,
decentralized trial components, landmark trials, and safety signal metrics.

Author: Adam Jones
Date: March 2026
"""

from typing import Dict, Any


# ═══════════════════════════════════════════════════════════════════════════════
# KNOWLEDGE BASE VERSION
# ═══════════════════════════════════════════════════════════════════════════════

KNOWLEDGE_VERSION: Dict[str, Any] = {
    "version": "2.0.0",
    "last_updated": "2026-03-22",
    "revision_notes": "Expanded release — 13 therapeutic areas, 7 trial phases, "
                      "9 regulatory agencies, 9 endpoint types, 9 adaptive designs, "
                      "9 biomarker strategies, 9 DCT components, 40 landmark trials, "
                      "6 safety signal metrics.",
    "sources": [
        "ICH E6(R2) / E6(R3) Good Clinical Practice",
        "ICH E8(R1) General Considerations for Clinical Studies",
        "ICH E9(R1) Statistical Principles — Estimands",
        "FDA Guidance: Adaptive Designs for Clinical Trials (2019)",
        "FDA Guidance: Decentralized Clinical Trials (2023)",
        "EMA Guideline on Multiplicity Issues in Clinical Trials",
        "ClinicalTrials.gov registry data",
        "FDA CDER Drug Approval Reports 2015-2025",
        "PhRMA Clinical Trial Success Rates 2011-2020",
        "BIO/QLS Advisors Clinical Development Success Rates 2011-2020",
    ],
    "counts": {
        "therapeutic_areas": 13,
        "trial_phases": 7,
        "regulatory_agencies": 9,
        "endpoint_types": 9,
        "adaptive_designs": 9,
        "biomarker_strategies": 9,
        "dct_components": 9,
        "landmark_trials": 40,
        "safety_signal_metrics": 6,
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. THERAPEUTIC AREAS
# ═══════════════════════════════════════════════════════════════════════════════

THERAPEUTIC_AREAS: Dict[str, Dict[str, Any]] = {
    "oncology": {
        "name": "Oncology",
        "description": "Cancer therapeutics including solid tumors, hematologic malignancies, "
                       "and immuno-oncology approaches.",
        "key_indications": [
            "non-small cell lung cancer", "breast cancer", "colorectal cancer",
            "melanoma", "acute myeloid leukemia", "multiple myeloma",
            "renal cell carcinoma", "hepatocellular carcinoma",
        ],
        "avg_trial_duration": "3-5 years",
        "success_rate_phase3": 0.36,
        "key_biomarkers": ["PD-L1", "TMB", "MSI-H/dMMR", "HER2", "EGFR", "ALK", "BRCA1/2", "KRAS G12C"],
        "regulatory_considerations": [
            "Accelerated approval based on surrogate endpoints (ORR, PFS)",
            "Breakthrough therapy designation common",
            "Companion diagnostics often co-developed",
            "RTOR (Real-Time Oncology Review) program",
            "Project Orbis for simultaneous multi-country review",
        ],
    },
    "cardiovascular": {
        "name": "Cardiovascular",
        "description": "Heart and vascular diseases including heart failure, coronary artery "
                       "disease, arrhythmias, and cardiomyopathies.",
        "key_indications": [
            "heart failure", "acute coronary syndrome", "atrial fibrillation",
            "hypertension", "pulmonary arterial hypertension", "dyslipidemia",
        ],
        "avg_trial_duration": "3-6 years",
        "success_rate_phase3": 0.35,
        "key_biomarkers": ["NT-proBNP", "troponin", "LDL-C", "Lp(a)", "hsCRP", "PCSK9"],
        "regulatory_considerations": [
            "Cardiovascular outcome trials (CVOTs) often required",
            "Large enrollment for MACE endpoints",
            "FDA CV Safety Guidance for diabetes drugs",
            "DSMB oversight for mortality endpoints",
        ],
    },
    "neuroscience": {
        "name": "Neuroscience",
        "description": "CNS disorders including neurodegenerative diseases, psychiatric "
                       "conditions, and neurological conditions.",
        "key_indications": [
            "Alzheimer's disease", "Parkinson's disease", "multiple sclerosis",
            "major depressive disorder", "schizophrenia", "epilepsy", "migraine",
        ],
        "avg_trial_duration": "4-8 years",
        "success_rate_phase3": 0.28,
        "key_biomarkers": ["amyloid PET", "tau PET", "neurofilament light chain",
                           "p-tau217", "CSF biomarkers", "alpha-synuclein"],
        "regulatory_considerations": [
            "Accelerated approval based on amyloid reduction (AD)",
            "Long placebo-controlled periods challenging",
            "Subjective endpoints require validated PRO instruments",
            "CNS penetration pharmacokinetics critical",
        ],
    },
    "immunology": {
        "name": "Immunology",
        "description": "Autoimmune and inflammatory diseases spanning multiple organ systems.",
        "key_indications": [
            "rheumatoid arthritis", "psoriasis", "inflammatory bowel disease",
            "systemic lupus erythematosus", "atopic dermatitis", "ankylosing spondylitis",
        ],
        "avg_trial_duration": "3-5 years",
        "success_rate_phase3": 0.40,
        "key_biomarkers": ["CRP", "ESR", "anti-CCP", "RF", "ANA", "IL-6", "TNF-alpha"],
        "regulatory_considerations": [
            "Biosimilar competition for established biologics",
            "Treat-to-target trial designs common",
            "Long-term safety registries required",
            "Immunogenicity assessment (ADA) mandatory",
        ],
    },
    "infectious_disease": {
        "name": "Infectious Disease",
        "description": "Anti-infective therapies including antivirals, antibiotics, antifungals, "
                       "and vaccines.",
        "key_indications": [
            "HIV", "hepatitis B", "hepatitis C", "influenza",
            "COVID-19", "RSV", "antimicrobial resistance", "tuberculosis",
        ],
        "avg_trial_duration": "2-4 years",
        "success_rate_phase3": 0.58,
        "key_biomarkers": ["viral load", "CD4 count", "seroconversion rate",
                           "PCR positivity", "MIC values"],
        "regulatory_considerations": [
            "LPAD pathway for antibiotics (limited population)",
            "Emergency Use Authorization pathway",
            "Surrogate virologic endpoints often accepted",
            "Challenge trial designs for some pathogens",
        ],
    },
    "rare_diseases": {
        "name": "Rare Diseases",
        "description": "Orphan drugs targeting conditions affecting fewer than 200,000 patients "
                       "in the US.",
        "key_indications": [
            "spinal muscular atrophy", "Duchenne muscular dystrophy",
            "cystic fibrosis", "sickle cell disease", "hemophilia",
            "Fabry disease", "Huntington's disease",
        ],
        "avg_trial_duration": "3-6 years",
        "success_rate_phase3": 0.45,
        "key_biomarkers": ["genetic mutation status", "enzyme activity levels",
                           "functional biomarkers specific to condition"],
        "regulatory_considerations": [
            "Orphan Drug Designation — 7-year market exclusivity (US)",
            "Smaller trial sizes acceptable",
            "Natural history studies as external controls",
            "Single-arm trials may be sufficient",
            "Pediatric-first development common",
        ],
    },
    "metabolic": {
        "name": "Metabolic / Endocrinology",
        "description": "Metabolic disorders including diabetes, obesity, and lipid disorders.",
        "key_indications": [
            "type 2 diabetes", "type 1 diabetes", "obesity",
            "NASH/MASH", "dyslipidemia", "gout", "osteoporosis",
        ],
        "avg_trial_duration": "3-5 years",
        "success_rate_phase3": 0.42,
        "key_biomarkers": ["HbA1c", "fasting glucose", "BMI", "ALT/AST",
                           "liver fibrosis score", "LDL-C", "uric acid"],
        "regulatory_considerations": [
            "Cardiovascular outcome trials required for diabetes drugs",
            "Histologic endpoints for NASH (liver biopsy)",
            "Long-term weight maintenance data expected",
            "Combination therapy trials common",
        ],
    },
    "respiratory": {
        "name": "Respiratory",
        "description": "Pulmonary and airway diseases including asthma, COPD, and fibrotic conditions.",
        "key_indications": [
            "asthma", "COPD", "idiopathic pulmonary fibrosis",
            "cystic fibrosis", "pulmonary arterial hypertension",
        ],
        "avg_trial_duration": "2-4 years",
        "success_rate_phase3": 0.39,
        "key_biomarkers": ["FEV1", "eosinophil count", "FeNO", "IgE",
                           "FVC", "6MWD", "DLCO"],
        "regulatory_considerations": [
            "Lung function (FEV1) as primary endpoint well established",
            "Exacerbation rate reduction as co-primary endpoint",
            "Inhaled drug device combination products (505(b)(2))",
            "Biomarker-guided therapy (eosinophilic phenotyping)",
        ],
    },
    "hematology": {
        "name": "Hematology",
        "description": "Blood disorders including malignancies, coagulation disorders, and anemias.",
        "key_indications": [
            "acute myeloid leukemia", "chronic lymphocytic leukemia",
            "multiple myeloma", "myelodysplastic syndromes",
            "sickle cell disease", "hemophilia", "immune thrombocytopenia",
        ],
        "avg_trial_duration": "3-5 years",
        "success_rate_phase3": 0.41,
        "key_biomarkers": ["MRD status", "complete response rate", "transfusion independence",
                           "factor levels", "hemoglobin", "platelet count"],
        "regulatory_considerations": [
            "MRD as potential surrogate endpoint",
            "Accelerated approval based on response rates",
            "Gene therapy regulatory pathway (BLA)",
            "CAR-T cell therapy REMS requirements",
        ],
    },
    "gastroenterology": {
        "name": "Gastroenterology",
        "description": "Gastrointestinal and hepatic diseases.",
        "key_indications": [
            "Crohn's disease", "ulcerative colitis", "GERD",
            "celiac disease", "NASH/MASH", "primary biliary cholangitis",
        ],
        "avg_trial_duration": "3-5 years",
        "success_rate_phase3": 0.38,
        "key_biomarkers": ["fecal calprotectin", "CRP", "endoscopic scores",
                           "histologic remission", "ALT/AST", "liver stiffness"],
        "regulatory_considerations": [
            "Endoscopic endpoints increasingly required",
            "Histologic-endoscopic remission as composite",
            "Liver biopsy for NASH — non-invasive alternatives emerging",
            "Biosimilar competition in IBD space",
        ],
    },
    "dermatology": {
        "name": "Dermatology",
        "description": "Skin and dermatologic conditions.",
        "key_indications": [
            "atopic dermatitis", "psoriasis", "hidradenitis suppurativa",
            "alopecia areata", "vitiligo", "acne",
        ],
        "avg_trial_duration": "2-4 years",
        "success_rate_phase3": 0.44,
        "key_biomarkers": ["EASI score", "PASI score", "IGA response",
                           "BSA involvement", "DLQI"],
        "regulatory_considerations": [
            "IGA 0/1 response as standard primary endpoint",
            "Photographic documentation requirements",
            "Topical vs systemic regulatory pathways differ",
            "Patient-reported outcomes (DLQI, itch NRS) as key secondaries",
        ],
    },
    "ophthalmology": {
        "name": "Ophthalmology",
        "description": "Eye diseases including retinal conditions, glaucoma, and dry eye.",
        "key_indications": [
            "age-related macular degeneration", "diabetic macular edema",
            "glaucoma", "dry eye disease", "geographic atrophy",
        ],
        "avg_trial_duration": "2-4 years",
        "success_rate_phase3": 0.43,
        "key_biomarkers": ["BCVA", "central retinal thickness", "IOP",
                           "geographic atrophy lesion area", "Schirmer's test"],
        "regulatory_considerations": [
            "Visual acuity (ETDRS letters) as primary endpoint",
            "Anti-VEGF injection burden as key differentiator",
            "Bilateral disease — fellow eye as control considerations",
            "Intravitreal delivery device considerations",
        ],
    },
    "gene_cell_therapy": {
        "name": "Gene and Cell Therapy",
        "description": "Advanced therapy medicinal products including gene therapies, "
                       "cell therapies, and tissue-engineered products.",
        "key_indications": [
            "spinal muscular atrophy", "hemophilia A/B", "sickle cell disease",
            "beta-thalassemia", "CAR-T for lymphoma/leukemia",
            "inherited retinal dystrophies", "Duchenne muscular dystrophy",
        ],
        "avg_trial_duration": "4-7 years",
        "success_rate_phase3": 0.40,
        "key_biomarkers": ["transgene expression levels", "vector copy number",
                           "functional protein levels", "T-cell expansion",
                           "cytokine levels", "integration site analysis"],
        "regulatory_considerations": [
            "BLA pathway (not NDA) for biologics",
            "Long-term follow-up (LTFU) — 5-15 years post-treatment",
            "Manufacturing consistency (CMC) challenges",
            "REMS programs for CAR-T (CRS, ICANS monitoring)",
            "Expedited programs: RMAT designation",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. TRIAL PHASES
# ═══════════════════════════════════════════════════════════════════════════════

TRIAL_PHASES: Dict[str, Dict[str, Any]] = {
    "preclinical": {
        "name": "Preclinical",
        "description": "Laboratory and animal studies to evaluate safety and biological activity "
                       "before human testing.",
        "typical_duration": "1-3 years",
        "typical_enrollment": "N/A (animal studies)",
        "primary_objectives": [
            "Pharmacology and mechanism of action",
            "ADME (absorption, distribution, metabolism, excretion)",
            "Acute and chronic toxicology",
            "Genotoxicity and carcinogenicity assessment",
        ],
        "key_endpoints": ["NOAEL", "therapeutic index", "PK parameters", "off-target effects"],
        "regulatory_requirements": [
            "IND-enabling GLP toxicology studies",
            "ICH M3(R2) nonclinical safety studies",
            "Species selection rationale",
        ],
    },
    "phase_0": {
        "name": "Phase 0 (Exploratory IND)",
        "description": "First-in-human microdose studies to evaluate pharmacokinetics and "
                       "pharmacodynamics at sub-therapeutic doses.",
        "typical_duration": "1-2 months",
        "typical_enrollment": "10-15 participants",
        "primary_objectives": [
            "Confirm mechanism of action in humans",
            "Evaluate PK at sub-therapeutic doses",
            "Assess target engagement",
        ],
        "key_endpoints": ["PK parameters", "target occupancy", "biodistribution"],
        "regulatory_requirements": [
            "Exploratory IND (reduced nonclinical package)",
            "FDA Guidance on Exploratory IND Studies (2006)",
            "Microdose (<100 mcg or 1/100th of NOAEL)",
        ],
    },
    "phase_1": {
        "name": "Phase 1",
        "description": "First major testing in humans. Evaluates safety, tolerability, "
                       "pharmacokinetics, and pharmacodynamics.",
        "typical_duration": "6-12 months",
        "typical_enrollment": "20-80 participants",
        "primary_objectives": [
            "Determine maximum tolerated dose (MTD)",
            "Characterize dose-limiting toxicities (DLTs)",
            "Establish pharmacokinetic profile",
            "Identify recommended phase 2 dose (RP2D)",
        ],
        "key_endpoints": ["MTD", "DLT rate", "PK parameters (Cmax, AUC, t1/2)",
                          "food effect", "drug-drug interactions"],
        "regulatory_requirements": [
            "IND application (FDA) or CTA (EMA)",
            "IRB/IEC approval",
            "Informed consent",
            "ICH E6(R2) GCP compliance",
            "Safety monitoring and expedited SAE reporting",
        ],
    },
    "phase_2": {
        "name": "Phase 2",
        "description": "Evaluate efficacy in target patient population and refine dosing. "
                       "Phase 2a (proof of concept) and Phase 2b (dose-finding).",
        "typical_duration": "1-2 years",
        "typical_enrollment": "100-300 participants",
        "primary_objectives": [
            "Establish proof of concept (Phase 2a)",
            "Determine optimal dose and regimen (Phase 2b)",
            "Preliminary efficacy assessment",
            "Further safety characterization",
        ],
        "key_endpoints": ["disease-specific efficacy endpoints", "dose-response relationship",
                          "biomarker changes", "adverse event profile"],
        "regulatory_requirements": [
            "End-of-Phase-1 meeting with FDA (recommended)",
            "Protocol amendments for dose adjustments",
            "DSMB for controlled trials",
            "CRF/EDC data capture per GCP",
        ],
    },
    "phase_3": {
        "name": "Phase 3 (Pivotal)",
        "description": "Confirmatory trials to demonstrate efficacy and safety for regulatory "
                       "approval. Usually randomized, controlled, and adequately powered.",
        "typical_duration": "2-4 years",
        "typical_enrollment": "300-3000+ participants",
        "primary_objectives": [
            "Confirm efficacy vs. standard of care or placebo",
            "Establish benefit-risk profile",
            "Support labeling claims",
            "Characterize safety in larger population",
        ],
        "key_endpoints": ["primary efficacy endpoint (per indication)",
                          "key secondary endpoints", "safety/tolerability",
                          "health-related quality of life"],
        "regulatory_requirements": [
            "End-of-Phase-2 / Pre-Phase-3 meeting with FDA",
            "Special Protocol Assessment (SPA) available",
            "Two adequate and well-controlled studies (or one large)",
            "Statistical Analysis Plan (SAP) finalized before unblinding",
            "GMP manufacturing at commercial scale",
        ],
    },
    "phase_4": {
        "name": "Phase 4 (Post-Marketing)",
        "description": "Post-approval studies to monitor long-term safety, effectiveness in "
                       "broader populations, and new indications.",
        "typical_duration": "Ongoing (2-10+ years)",
        "typical_enrollment": "1000-100000+ participants",
        "primary_objectives": [
            "Long-term safety surveillance",
            "Real-world effectiveness",
            "New population or indication expansion",
            "Fulfill post-marketing requirements (PMRs/PMCs)",
        ],
        "key_endpoints": ["long-term safety outcomes", "real-world effectiveness",
                          "rare adverse event detection", "comparative effectiveness"],
        "regulatory_requirements": [
            "Post-Marketing Requirements (PMR) / Commitments (PMC)",
            "REMS compliance if applicable",
            "Periodic Safety Update Reports (PSURs)",
            "Pharmacovigilance system maintenance",
        ],
    },
    "expanded_access": {
        "name": "Expanded Access / Compassionate Use",
        "description": "Pre-approval access for seriously ill patients with no alternative "
                       "treatment options.",
        "typical_duration": "Variable",
        "typical_enrollment": "Individual to cohort-based",
        "primary_objectives": [
            "Provide access to investigational drug for unmet need",
            "Collect additional safety data",
        ],
        "key_endpoints": ["safety and tolerability", "clinical benefit assessment"],
        "regulatory_requirements": [
            "FDA Form 3926 (individual patients)",
            "IRB/IEC review and approval",
            "Informed consent with additional disclosures",
            "Manufacturer must agree to supply drug",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. REGULATORY AGENCIES
# ═══════════════════════════════════════════════════════════════════════════════

REGULATORY_AGENCIES: Dict[str, Dict[str, Any]] = {
    "FDA": {
        "name": "U.S. Food and Drug Administration",
        "country": "United States",
        "approval_pathway": [
            "NDA (New Drug Application) — small molecules",
            "BLA (Biologics License Application) — biologics",
            "505(b)(2) — modified new drugs referencing existing data",
            "ANDA (Abbreviated NDA) — generics",
        ],
        "avg_review_time": "10-12 months (standard), 6-8 months (priority)",
        "key_requirements": [
            "IND submission before clinical trials",
            "Two adequate and well-controlled studies",
            "GMP-compliant manufacturing",
            "Complete safety database (ICH E1)",
            "Pediatric study plan (PSP) submission",
        ],
        "expedited_programs": [
            "Breakthrough Therapy Designation (BTD)",
            "Fast Track Designation",
            "Accelerated Approval (subpart H/I)",
            "Priority Review",
            "RTOR (Real-Time Oncology Review)",
            "RMAT (Regenerative Medicine Advanced Therapy)",
        ],
    },
    "EMA": {
        "name": "European Medicines Agency",
        "country": "European Union",
        "approval_pathway": [
            "Centralised Procedure (mandatory for biotech, orphan, HIV, cancer, etc.)",
            "Mutual Recognition Procedure (MRP)",
            "Decentralised Procedure (DCP)",
            "National Procedure",
        ],
        "avg_review_time": "12-15 months (centralised)",
        "key_requirements": [
            "Clinical Trial Regulation EU 536/2014 (CTIS portal)",
            "EU Common Technical Document (CTD) format",
            "Risk Management Plan (RMP) mandatory",
            "Paediatric Investigation Plan (PIP) or waiver",
            "GMP certification by EU authority",
        ],
        "expedited_programs": [
            "PRIME (PRIority MEdicines)",
            "Accelerated Assessment (150-day timeline)",
            "Conditional Marketing Authorisation",
            "Marketing Authorisation under Exceptional Circumstances",
        ],
    },
    "PMDA": {
        "name": "Pharmaceuticals and Medical Devices Agency",
        "country": "Japan",
        "approval_pathway": [
            "New Drug Application (J-NDA)",
            "SAKIGAKE Designation (pioneer review)",
            "Conditional Early Approval System",
        ],
        "avg_review_time": "12 months (standard), 6 months (SAKIGAKE priority)",
        "key_requirements": [
            "Japanese bridging studies may be required",
            "CTN (Clinical Trial Notification) before trials",
            "J-GCP compliance",
            "Japanese patient data (ethnic sensitivity per ICH E5)",
        ],
        "expedited_programs": [
            "SAKIGAKE Designation System",
            "Conditional Early Approval",
            "Priority Review",
            "Orphan Drug Designation",
        ],
    },
    "Health_Canada": {
        "name": "Health Canada",
        "country": "Canada",
        "approval_pathway": [
            "New Drug Submission (NDS)",
            "Abbreviated New Drug Submission (ANDS)",
            "Supplement to NDS (SNDS)",
        ],
        "avg_review_time": "12 months (standard), 6 months (priority)",
        "key_requirements": [
            "Clinical Trial Application (CTA) with 30-day default",
            "ICH CTD format",
            "Canadian reference product for generics",
            "Drug Identification Number (DIN) required",
        ],
        "expedited_programs": [
            "Priority Review",
            "Notice of Compliance with Conditions (NOC/c)",
            "Agile Licensing",
        ],
    },
    "TGA": {
        "name": "Therapeutic Goods Administration",
        "country": "Australia",
        "approval_pathway": [
            "Category 1 (new chemical entity evaluation)",
            "Category 3 (bioequivalence / biosimilar)",
            "Provisional Approval (2+ years conditional)",
        ],
        "avg_review_time": "11 months (standard), 5-6 months (priority)",
        "key_requirements": [
            "CTN (Clinical Trial Notification) scheme",
            "CTX (Clinical Trial Exemption) for higher-risk trials",
            "Australian Register of Therapeutic Goods (ARTG) listing",
            "Alignment with ICH guidelines",
        ],
        "expedited_programs": [
            "Priority Review",
            "Provisional Approval Pathway",
            "Orphan Drug Program",
        ],
    },
    "MHRA": {
        "name": "Medicines and Healthcare products Regulatory Agency",
        "country": "United Kingdom",
        "approval_pathway": [
            "National Procedure (post-Brexit)",
            "International Recognition Procedure",
            "Innovative Licensing and Access Pathway (ILAP)",
        ],
        "avg_review_time": "12 months (standard), 150 days (accelerated)",
        "key_requirements": [
            "Clinical Trial Authorisation (CTA)",
            "Combined review with HRA (Health Research Authority)",
            "ICH CTD format",
            "GMP compliance with MHRA inspection",
        ],
        "expedited_programs": [
            "ILAP (Innovative Licensing and Access Pathway)",
            "Promising Innovative Medicine (PIM) designation",
            "EAMS (Early Access to Medicines Scheme)",
            "Rolling Review",
        ],
    },
    "NMPA": {
        "name": "China National Medical Products Administration",
        "country": "China",
        "approval_pathway": [
            "Category 1 New Drug Application (innovative drug)",
            "Category 2 Modified New Drug",
            "Category 5 Imported Drug Application",
            "Generic Drug Application (ANDA equivalent)",
        ],
        "avg_review_time": "~300 days (standard), 130 days (breakthrough therapy)",
        "key_requirements": [
            "Clinical Trial Application (CTA) with 60-day default approval",
            "ICH CTD format adopted since 2020",
            "Chinese patient data may be required for ethnic sensitivity",
            "GMP compliance with NMPA inspection",
            "Drug Master File (DMF) system",
        ],
        "expedited_programs": [
            "Breakthrough Therapy Designation",
            "Conditional Approval Pathway",
            "Priority Review (for innovative drugs, rare diseases)",
            "Special Approval Procedure (public health emergencies)",
        ],
    },
    "Swissmedic": {
        "name": "Swissmedic — Swiss Agency for Therapeutic Products",
        "country": "Switzerland",
        "approval_pathway": [
            "National Authorisation Procedure",
            "Authorisation by Notification (known active substances)",
            "Fast-track Authorisation Procedure",
            "Temporary Authorisation (pandemic/emergency)",
        ],
        "avg_review_time": "330 days (standard), 140 days (fast-track)",
        "key_requirements": [
            "Clinical Trial Application via BASEC portal",
            "ICH CTD format",
            "GMP compliance with Swissmedic inspection",
            "Periodic Safety Update Reports (PSURs)",
        ],
        "expedited_programs": [
            "Fast-track Authorisation (serious conditions, unmet need)",
            "Temporary Authorisation for Pandemic Use",
            "International cooperation agreements (FDA, EMA, TGA, Health Canada)",
            "ACSS consortium (Australia, Canada, Singapore, Switzerland)",
        ],
    },
    "ANVISA": {
        "name": "Agencia Nacional de Vigilancia Sanitaria",
        "country": "Brazil",
        "approval_pathway": [
            "New Drug Registration (Registro de Medicamento Novo)",
            "Generic Drug Registration",
            "Similar Drug Registration (branded generic)",
            "Biological Product Registration",
        ],
        "avg_review_time": "365 days (standard), 120 days (priority)",
        "key_requirements": [
            "Clinical Trial approval via CEP/CONEP ethics system",
            "ICH CTD format accepted",
            "GMP certification by ANVISA",
            "Stability studies under Zone IVb climatic conditions",
            "Local representative required",
        ],
        "expedited_programs": [
            "ANVISA Priority Review (rare diseases, neglected diseases, public health emergencies)",
            "Temporary Emergency Use Authorisation",
            "Simplified registration for WHO prequalified products",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. ENDPOINT TYPES
# ═══════════════════════════════════════════════════════════════════════════════

ENDPOINT_TYPES: Dict[str, Dict[str, Any]] = {
    "primary": {
        "description": "The main outcome measure used to determine whether a trial "
                       "meets its primary objective. Pre-specified in the SAP.",
        "examples": [
            "Overall survival (OS) in oncology",
            "HbA1c reduction in diabetes",
            "ACR20 response in rheumatoid arthritis",
            "MACE (major adverse cardiovascular events)",
            "FEV1 change from baseline in respiratory",
        ],
        "statistical_methods": [
            "Intent-to-treat (ITT) analysis",
            "Alpha-controlled hypothesis testing",
            "Pre-specified estimand framework (ICH E9(R1))",
            "Sample size powered for primary endpoint",
        ],
        "regulatory_acceptance": "Required for approval; must demonstrate clinical meaningfulness.",
    },
    "secondary": {
        "description": "Additional endpoints that support the primary endpoint or evaluate "
                       "other important treatment effects.",
        "examples": [
            "Progression-free survival (when OS is primary)",
            "Quality of life measures (EQ-5D, SF-36)",
            "Time to first exacerbation",
            "Change in biomarker from baseline",
        ],
        "statistical_methods": [
            "Hierarchical (fixed-sequence) testing",
            "Hochberg or Holm step-down procedures",
            "Gatekeeping strategies",
            "Multiplicity adjustment required for key secondaries",
        ],
        "regulatory_acceptance": "Supports labeling claims if multiplicity-controlled.",
    },
    "exploratory": {
        "description": "Hypothesis-generating endpoints not powered for statistical significance. "
                       "Used for future trial design.",
        "examples": [
            "Subgroup efficacy analyses",
            "Biomarker correlations",
            "Dose-exposure-response relationships",
            "Novel digital health metrics",
        ],
        "statistical_methods": [
            "Descriptive statistics",
            "No formal multiplicity adjustment",
            "Confidence intervals (not p-values)",
            "Bayesian exploratory models",
        ],
        "regulatory_acceptance": "Not used for approval decisions; hypothesis-generating only.",
    },
    "safety": {
        "description": "Endpoints evaluating the adverse event profile, tolerability, "
                       "and safety signals of the investigational agent.",
        "examples": [
            "Treatment-emergent adverse events (TEAEs)",
            "Serious adverse events (SAEs)",
            "Adverse events of special interest (AESIs)",
            "Laboratory abnormalities (liver, renal, hematologic)",
            "ECG changes (QTc prolongation)",
        ],
        "statistical_methods": [
            "Exposure-adjusted incidence rates",
            "Kaplan-Meier for time-to-event safety",
            "MedDRA coding (SOC → PT)",
            "Hepatotoxicity assessment (Hy's Law)",
        ],
        "regulatory_acceptance": "Integral to benefit-risk assessment; comprehensive safety "
                                 "database required (ICH E1).",
    },
    "patient_reported": {
        "description": "Outcomes reported directly by patients without clinician interpretation. "
                       "Captures treatment impact on daily life.",
        "examples": [
            "Patient Global Impression of Change (PGI-C)",
            "EORTC QLQ-C30 (oncology QoL)",
            "PHQ-9 (depression severity)",
            "DLQI (dermatology QoL)",
            "KOOS (knee injury and osteoarthritis)",
        ],
        "statistical_methods": [
            "Mixed-effects models for repeated measures (MMRM)",
            "Responder analysis (MCID threshold)",
            "Pattern mixture models for missing data",
            "FDA PRO Guidance (2009) validation requirements",
        ],
        "regulatory_acceptance": "Accepted for labeling claims if validated instrument used per "
                                 "FDA PRO Guidance; increasingly important.",
    },
    "digital": {
        "description": "Endpoints derived from digital health technologies (wearables, sensors, "
                       "smartphone apps) providing continuous real-world data.",
        "examples": [
            "Step count / physical activity (accelerometry)",
            "Continuous glucose monitoring (CGM) metrics",
            "Digital cognitive assessments",
            "Sleep quality metrics (actigraphy)",
            "Remote spirometry (FEV1)",
        ],
        "statistical_methods": [
            "Time-series analysis",
            "Mixed models for longitudinal digital data",
            "Algorithm-derived composite scores",
            "Verification & analytical validation (V3 framework)",
        ],
        "regulatory_acceptance": "Emerging acceptance; FDA Digital Health Technologies guidance. "
                                 "Requires analytical and clinical validation.",
    },
    "composite": {
        "description": "Pre-specified combination of multiple individual endpoints into a single "
                       "outcome measure, typically used when single endpoints are insufficiently "
                       "frequent or clinically incomplete.",
        "examples": [
            "MACE (major adverse cardiovascular events: CV death, MI, stroke)",
            "PFS or OS co-primary in oncology",
            "DAS28 composite in rheumatoid arthritis",
            "Composite of sustained eGFR decline, ESKD, or renal death in nephrology",
        ],
        "statistical_methods": [
            "Win ratio (Pocock et al.) for prioritized components",
            "Hierarchical testing of individual components",
            "Time-to-first-event analysis (standard Cox model)",
            "Finkelstein-Schoenfeld method for prioritized composites",
        ],
        "regulatory_acceptance": "Widely accepted in cardiovascular (MACE) and nephrology; "
                                 "FDA expects clinical meaningfulness of each component.",
    },
    "minimal_residual_disease": {
        "description": "Detection of residual cancer cells below the threshold of conventional "
                       "morphologic assessment, used as a sensitive molecular response endpoint "
                       "in hematologic malignancies.",
        "examples": [
            "MRD negativity at 10^-5 or 10^-6 in multiple myeloma (NGS/MFC)",
            "MRD negativity in CLL (iwCLL criteria)",
            "MRD negativity in ALL (flow cytometry or PCR)",
            "MRD-guided treatment discontinuation in CML",
        ],
        "statistical_methods": [
            "MRD negativity rate as binary endpoint",
            "Landmark MRD analysis at pre-specified timepoints",
            "MRD-free survival (time from MRD negativity to progression/death)",
            "Concordance between MRD assays (NGS vs flow cytometry)",
        ],
        "regulatory_acceptance": "FDA draft guidance (2020) on MRD as endpoint in hematologic "
                                 "malignancies; accepted as reasonably likely surrogate for "
                                 "accelerated approval. Key trials: MASTER (myeloma), FLAIR (CLL).",
    },
    "ctDNA_clearance": {
        "description": "Molecular response endpoint based on clearance of circulating tumor DNA "
                       "from plasma, serving as a real-time, non-invasive measure of treatment "
                       "response and residual disease.",
        "examples": [
            "ctDNA clearance at landmark timepoint post-treatment",
            "Molecular response (ctDNA reduction > 50%)",
            "ctDNA dynamics as early predictor of radiographic response",
            "Serial ctDNA monitoring for recurrence detection",
        ],
        "statistical_methods": [
            "Landmark analysis of ctDNA clearance vs non-clearance",
            "ctDNA-stratified survival analysis",
            "Longitudinal mixed models for ctDNA kinetics",
            "Sensitivity/specificity for recurrence prediction",
        ],
        "regulatory_acceptance": "Emerging; FDA has accepted ctDNA-guided endpoints in "
                                 "exploratory/secondary contexts. circularDNA trial (DYNAMIC, "
                                 "BESPOKE CRC) providing evidence for ctDNA-guided adjuvant therapy.",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ADAPTIVE DESIGNS
# ═══════════════════════════════════════════════════════════════════════════════

ADAPTIVE_DESIGNS: Dict[str, Dict[str, Any]] = {
    "group_sequential": {
        "description": "Pre-planned interim analyses allowing early stopping for efficacy "
                       "or futility while controlling type I error.",
        "when_to_use": [
            "Large confirmatory trials with long enrollment periods",
            "Ethical imperative to stop early if treatment clearly superior",
            "High cost per patient necessitating efficiency",
        ],
        "advantages": [
            "Reduces expected sample size under alternative hypothesis",
            "Ethical — stops early when efficacy proven",
            "Well-established regulatory acceptance",
            "Preserves type I error control",
        ],
        "disadvantages": [
            "Requires pre-specified analysis schedule",
            "Operational complexity of interim unblinding",
            "Alpha spending limits power at final analysis",
        ],
        "regulatory_guidance": "FDA Guidance on Adaptive Designs (2019); ICH E9(R1)",
        "key_trials": [
            "RECOVERY trial (COVID-19) — group sequential with interim analyses",
            "DAPA-HF (dapagliflozin) — stopped early for efficacy",
        ],
    },
    "sample_size_reestimation": {
        "description": "Allows re-calculation of sample size at an interim analysis based on "
                       "observed data (effect size, variance).",
        "when_to_use": [
            "Uncertain effect size from Phase 2 data",
            "High variability in the primary endpoint",
            "Novel endpoints with limited historical data",
        ],
        "advantages": [
            "Protects against underpowered studies",
            "Uses blinded or unblinded approaches",
            "Can maintain type I error under blinded SSR",
        ],
        "disadvantages": [
            "Unblinded SSR requires careful type I error control",
            "Budget and timeline uncertainty for sponsors",
            "Operational complexity in implementation",
        ],
        "regulatory_guidance": "FDA Guidance on Adaptive Designs (2019); "
                               "CHMP Reflection Paper on Methodological Issues",
        "key_trials": [
            "FOCUS trial (stroke) — blinded SSR for primary outcome variance",
        ],
    },
    "response_adaptive": {
        "description": "Randomization probabilities change during the trial based on "
                       "accumulating outcome data, allocating more patients to better-performing arms.",
        "when_to_use": [
            "Multi-arm dose-finding studies",
            "Rare diseases with limited patient populations",
            "Platform trials evaluating multiple treatments",
        ],
        "advantages": [
            "More patients receive better-performing treatments",
            "Ethical benefit in randomization allocation",
            "Efficient exploration of multiple doses/treatments",
        ],
        "disadvantages": [
            "Operational complexity in real-time data analysis",
            "Potential for time-trend bias",
            "Regulatory concerns about type I error inflation",
            "More complex statistical analysis",
        ],
        "regulatory_guidance": "FDA Guidance on Adaptive Designs (2019); "
                               "Bayesian framework often used (FDA Bayesian Guidance 2010)",
        "key_trials": [
            "I-SPY 2 (breast cancer) — Bayesian adaptive randomization",
        ],
    },
    "biomarker_adaptive": {
        "description": "Trial design adapts based on biomarker results, such as enriching "
                       "enrollment for biomarker-positive patients or modifying treatment assignment.",
        "when_to_use": [
            "Targeted therapies with companion diagnostic development",
            "When biomarker-treatment interaction is hypothesized but uncertain",
            "Precision medicine approaches",
        ],
        "advantages": [
            "Identifies optimal target population",
            "Supports companion diagnostic co-development",
            "More efficient than separate biomarker-defined trials",
        ],
        "disadvantages": [
            "Requires validated biomarker assay early in development",
            "Complex statistical design and analysis",
            "Potential for false enrichment",
        ],
        "regulatory_guidance": "FDA Guidance on Enrichment Strategies (2019); "
                               "FDA Guidance on Co-Development of Drugs and Diagnostics",
        "key_trials": [
            "KEYNOTE-024 (pembrolizumab) — PD-L1 enrichment",
            "SHIVA trial — biomarker-matched targeted therapy",
        ],
    },
    "platform_trial": {
        "description": "Perpetual trial infrastructure testing multiple treatments against a "
                       "shared control arm. Arms can be added or dropped based on interim results.",
        "when_to_use": [
            "Multiple candidate treatments for same disease",
            "Pandemic response requiring rapid evaluation",
            "Disease areas with many emerging therapies",
        ],
        "advantages": [
            "Shared control arm increases efficiency",
            "Perpetual enrollment reduces startup costs",
            "Arms can be added without new trial infrastructure",
            "Faster evaluation of multiple treatments",
        ],
        "disadvantages": [
            "Complex governance and data sharing agreements",
            "Heterogeneity in control arm over time",
            "Logistical complexity in multi-arm management",
            "Regulatory challenges with changing comparators",
        ],
        "regulatory_guidance": "FDA Guidance: Master Protocols (2022); "
                               "EMA Reflection Paper on Complex Clinical Trials (2023)",
        "key_trials": [
            "RECOVERY (COVID-19) — largest platform trial (40,000+ patients)",
            "GBM AGILE (glioblastoma) — Bayesian platform",
            "I-SPY COVID — adaptive platform for COVID treatments",
        ],
    },
    "seamless_phase": {
        "description": "Combines two phases (e.g., Phase 2/3) into a single trial, using "
                       "data from the learning phase to inform the confirmatory phase.",
        "when_to_use": [
            "Accelerated development timelines (breakthrough therapies)",
            "Rare diseases with limited patient populations",
            "When dose selection and confirmation can share a framework",
        ],
        "advantages": [
            "Shorter overall development timeline",
            "Reduced patient numbers vs. separate trials",
            "No inter-trial gap for enrollment continuity",
        ],
        "disadvantages": [
            "Complex statistical design for type I error control",
            "Manufacturing at scale may not be ready for Phase 3 portion",
            "Regulatory engagement required early and often",
            "Difficulty in protocol amendments across phases",
        ],
        "regulatory_guidance": "FDA Guidance on Adaptive Designs (2019); "
                               "EMA Qualification Advice for Novel Methodologies",
        "key_trials": [
            "DETERMINE (abemaciclib) — seamless Phase 2/3 in breast cancer",
        ],
    },
    "master_protocol": {
        "description": "Overarching protocol framework that encompasses basket trials (one drug, "
                       "multiple diseases), umbrella trials (one disease, multiple drugs), and "
                       "platform trials (perpetual infrastructure) under a single governance structure.",
        "when_to_use": [
            "Precision medicine initiatives with multiple biomarker-drug pairings",
            "Disease areas with molecularly defined subtypes",
            "Multi-stakeholder collaboration across sponsors",
        ],
        "advantages": [
            "Shared infrastructure (sites, IRBs, data management) across sub-studies",
            "Efficient screening and biomarker testing",
            "Central governance and DSMB oversight",
            "Accelerated enrollment through common protocol",
        ],
        "disadvantages": [
            "Complex governance with multiple sponsors",
            "IP and data sharing legal challenges",
            "Heterogeneous patient populations across sub-studies",
            "Logistical burden of multi-drug supply chains",
        ],
        "regulatory_guidance": "FDA Guidance: Master Protocols — Efficient Clinical Trial Design "
                               "Strategies (2022); EMA Reflection Paper on Complex Clinical Trials (2023)",
        "key_trials": [
            "NCI-MATCH (basket) — molecular-matched therapy across tumor types",
            "Lung-MAP (umbrella) — biomarker-driven treatment in squamous NSCLC",
            "RECOVERY (platform) — COVID-19 treatments",
        ],
    },
    "enrichment_adaptive": {
        "description": "Adaptive enrichment design that allows modification of the target population "
                       "during the trial based on interim biomarker response data, restricting "
                       "enrollment to biomarker-positive subgroups if differential benefit is observed.",
        "when_to_use": [
            "Uncertain biomarker-treatment interaction at trial start",
            "When broad enrollment is desired initially but enrichment may be needed",
            "Targeted therapies with evolving companion diagnostic data",
        ],
        "advantages": [
            "Begins with broad population, enriches adaptively",
            "Preserves type I error while increasing power in target subgroup",
            "Avoids separate all-comers and enriched trials",
        ],
        "disadvantages": [
            "Complex interim decision rules and statistical methods",
            "Requires real-time biomarker data availability",
            "Potential for bias if enrichment criteria are data-driven",
            "Regulatory scrutiny of adaptive enrichment decisions",
        ],
        "regulatory_guidance": "FDA Guidance on Enrichment Strategies for Clinical Trials (2019); "
                               "FDA Guidance on Adaptive Designs (2019)",
        "key_trials": [
            "TAPPAS (trebananib) — adaptive enrichment by ascites status in ovarian cancer",
        ],
    },
    "dose_finding": {
        "description": "Model-based adaptive designs for Phase I/II dose optimization, using "
                       "statistical models to guide dose escalation and identify the optimal "
                       "biological dose (OBD) or maximum tolerated dose (MTD).",
        "when_to_use": [
            "First-in-human dose escalation studies",
            "Dose optimization for combination therapies",
            "Immunotherapy and targeted agents where MTD may not be optimal",
        ],
        "advantages": [
            "More patients treated at doses near optimal",
            "Model-guided decisions reduce trial size vs. 3+3 design",
            "Incorporates both toxicity and efficacy data (OBD designs)",
            "Better characterization of dose-response relationship",
        ],
        "disadvantages": [
            "Requires statistical modeling expertise",
            "Model misspecification risk",
            "Operational complexity in real-time dose decisions",
            "May require Safety Review Committee for each cohort",
        ],
        "regulatory_guidance": "FDA Project Optimus (dose optimization initiative); "
                               "FDA Guidance on Adaptive Designs (2019); ICH E4",
        "key_trials": [
            "Multiple Phase I oncology trials using CRM (continual reassessment method)",
            "BOIN (Bayesian optimal interval) designs in immuno-oncology",
            "BLRM (Bayesian logistic regression model) dose escalation",
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6. BIOMARKER STRATEGIES
# ═══════════════════════════════════════════════════════════════════════════════

BIOMARKER_STRATEGIES: Dict[str, Dict[str, Any]] = {
    "enrichment": {
        "description": "Restricts enrollment to patients who are biomarker-positive, "
                       "increasing the likelihood of observing a treatment effect.",
        "examples": [
            "HER2-positive breast cancer for trastuzumab",
            "PD-L1 >= 50% for first-line pembrolizumab monotherapy",
            "EGFR mutation-positive NSCLC for osimertinib",
        ],
        "regulatory_pathway": "FDA Guidance on Enrichment Strategies (2019); "
                              "co-development with companion diagnostic (CDx) required.",
    },
    "stratification": {
        "description": "Randomization is stratified by biomarker status to ensure balance "
                       "and allow pre-planned subgroup analyses.",
        "examples": [
            "PD-L1 expression levels (TPS <1%, 1-49%, >=50%)",
            "MSI-H vs MSS in colorectal cancer",
            "BRCA mutation status in ovarian cancer",
        ],
        "regulatory_pathway": "Standard practice; stratification factors pre-specified in SAP.",
    },
    "prognostic": {
        "description": "Biomarker predicts disease outcome regardless of treatment. "
                       "Used for risk stratification, not treatment selection.",
        "examples": [
            "Oncotype DX recurrence score in breast cancer",
            "CRP in cardiovascular risk assessment",
            "Stage and grade in oncology",
        ],
        "regulatory_pathway": "Supports patient selection and risk-benefit assessment; "
                              "does not require companion diagnostic.",
    },
    "predictive": {
        "description": "Biomarker predicts differential treatment response between "
                       "biomarker-positive and biomarker-negative patients.",
        "examples": [
            "EGFR mutations predicting response to TKIs",
            "ALK rearrangements predicting response to crizotinib/alectinib",
            "KRAS G12C predicting response to sotorasib/adagrasib",
        ],
        "regulatory_pathway": "Companion diagnostic co-development (PMA for CDx); "
                              "FDA Guidance on In Vitro Companion Diagnostic Devices.",
    },
    "pharmacodynamic": {
        "description": "Biomarker measures the biological effect of the drug on its target. "
                       "Used for dose selection and proof of mechanism.",
        "examples": [
            "Phosphorylated target protein levels post-treatment",
            "Receptor occupancy by PET imaging",
            "Cytokine levels post-immunotherapy",
        ],
        "regulatory_pathway": "Supports dose selection rationale; typically used in Phase 1/2 "
                              "and included in Clinical Pharmacology sections of NDA.",
    },
    "surrogate": {
        "description": "Biomarker used as a substitute for a clinical endpoint. Must be "
                       "reasonably likely to predict clinical benefit.",
        "examples": [
            "HbA1c as surrogate for diabetic complications",
            "Viral load suppression as surrogate for clinical outcomes (HIV)",
            "Tumor response rate (ORR) as surrogate for overall survival",
            "Blood pressure reduction as surrogate for stroke/MI prevention",
        ],
        "regulatory_pathway": "FDA Table of Surrogate Endpoints; Accelerated Approval pathway "
                              "requires confirmatory trial. Subpart H (drugs) / Subpart E (biologics).",
    },
    "companion_diagnostic": {
        "description": "Co-development of a diagnostic test (CDx) with a therapeutic product, "
                       "where the CDx is essential for the safe and effective use of the drug "
                       "by identifying the appropriate patient population.",
        "examples": [
            "FoundationOne CDx for multiple targeted therapies (tissue NGS panel)",
            "Ventana PD-L1 (SP142) for atezolizumab",
            "cobas EGFR Mutation Test v2 for osimertinib",
            "Guardant360 CDx (liquid biopsy) for osimertinib in NSCLC",
        ],
        "regulatory_pathway": "FDA PMA (Premarket Approval) pathway for CDx device; "
                              "co-development with drug per FDA Guidance on In Vitro Companion "
                              "Diagnostic Devices (2014, updated 2023). Contemporaneous approval required.",
    },
    "liquid_biopsy": {
        "description": "Non-invasive biomarker strategy using circulating tumor DNA (ctDNA), "
                       "cell-free DNA (cfDNA), circulating tumor cells (CTCs), or exosomes "
                       "from blood samples for genotyping, monitoring, and MRD detection.",
        "examples": [
            "ctDNA genotyping for EGFR T790M resistance in NSCLC (Guardant360)",
            "Serial ctDNA monitoring for treatment response in CRC (DYNAMIC trial)",
            "cfDNA methylation patterns for early cancer detection (GRAIL Galleri)",
            "MRD detection via ctDNA in post-surgical monitoring",
        ],
        "regulatory_pathway": "FDA-approved liquid biopsy CDx platforms available (Guardant360 CDx, "
                              "FoundationOne Liquid CDx). FDA guidance evolving for ctDNA-based "
                              "MRD and early detection applications.",
    },
    "digital_biomarker": {
        "description": "Biomarkers derived from wearable devices, smartphones, and digital sensors "
                       "that capture continuous physiological and behavioral data. Enables "
                       "objective, real-world measurement of disease burden and treatment effect.",
        "examples": [
            "Actigraphy-derived sleep and activity metrics in depression/CNS trials",
            "Digital phenotyping via smartphone passive sensing (typing speed, gait)",
            "Continuous glucose monitoring (CGM) as pharmacodynamic biomarker",
            "Digital spirometry for remote pulmonary function assessment",
        ],
        "regulatory_pathway": "FDA Digital Health Center of Excellence; FDA DHT guidance for "
                              "drug development (2023). Requires V3 framework validation "
                              "(Verification, Analytical Validation, Clinical Validation). "
                              "CTTI recommendations for sensor-based endpoints.",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 7. DECENTRALIZED TRIAL COMPONENTS (DCT)
# ═══════════════════════════════════════════════════════════════════════════════

DECENTRALIZED_COMPONENTS: Dict[str, Dict[str, Any]] = {
    "econsent": {
        "description": "Electronic informed consent using digital platforms, allowing "
                       "remote review, multimedia education, and digital signatures.",
        "regulatory_status": "FDA Guidance on Use of Electronic Informed Consent (2016). "
                             "21 CFR Part 11 compliance for electronic records/signatures.",
        "data_quality_considerations": [
            "Audit trail for consent versioning",
            "Identity verification for digital signatures",
            "Accessibility compliance (Section 508)",
            "Comprehension assessment integration",
        ],
        "patient_preference": "High — convenience, ability to review at own pace, "
                              "share with family before signing.",
    },
    "telemedicine": {
        "description": "Remote clinical visits via video or phone, replacing some "
                       "in-person site visits for assessments and monitoring.",
        "regulatory_status": "FDA Guidance on Conduct of Clinical Trials During COVID-19 "
                             "(2020, updated 2023). State licensure requirements apply.",
        "data_quality_considerations": [
            "Standardized virtual visit checklists",
            "Video quality requirements for visual assessments",
            "Documentation of visit modality in CRF",
            "Limited physical examination capabilities",
        ],
        "patient_preference": "Very high — eliminates travel, reduces burden; "
                              "preferred for routine follow-up visits.",
    },
    "home_health": {
        "description": "Qualified nurses or phlebotomists visit patient homes for "
                       "study procedures, sample collection, and drug administration.",
        "regulatory_status": "Permitted under most regulatory frameworks. FDA Guidance on "
                             "DCTs (2023) addresses home-based procedures.",
        "data_quality_considerations": [
            "Training and certification of home health providers",
            "Sample handling and cold chain logistics",
            "Protocol deviation reporting for home visit deviations",
            "IP accountability and chain of custody",
        ],
        "patient_preference": "High — especially for elderly, disabled, or rural patients.",
    },
    "local_labs": {
        "description": "Use of local/community laboratories instead of central labs "
                       "for routine safety monitoring samples.",
        "regulatory_status": "Acceptable with standardization. CLIA certification "
                             "required in US. Cross-lab harmonization needed.",
        "data_quality_considerations": [
            "Reference range harmonization across labs",
            "Central lab for primary endpoint biomarkers (recommended)",
            "CLIA/ISO 15189 certification requirements",
            "Data integration from multiple LIMS systems",
        ],
        "patient_preference": "Moderate — proximity benefit; some prefer consistency of central lab.",
    },
    "wearables": {
        "description": "Continuous data collection via wearable devices (smartwatches, "
                       "patches, biosensors) for objective physiological monitoring.",
        "regulatory_status": "FDA Digital Health Technologies guidance (2023). "
                             "Device must be fit-for-purpose validated.",
        "data_quality_considerations": [
            "Device calibration and wear-time compliance",
            "Algorithm validation for derived endpoints",
            "Data volume management and storage",
            "Artifact detection and data cleaning pipelines",
            "V3 framework (Verification, Analytical Validation, Clinical Validation)",
        ],
        "patient_preference": "Moderate to high — depends on device burden and comfort.",
    },
    "epro_ecoa": {
        "description": "Electronic patient-reported outcomes and clinician-assessed outcomes "
                       "collected via apps, tablets, or web portals.",
        "regulatory_status": "FDA PRO Guidance (2009). ISPOR best practices for ePRO. "
                             "Instrument migration from paper requires equivalence testing.",
        "data_quality_considerations": [
            "Compliance rates and reminder systems",
            "Response window enforcement",
            "Device provisioning vs BYOD considerations",
            "Linguistic and cultural validation of instruments",
        ],
        "patient_preference": "High — real-time capture is preferred over recall-based reporting.",
    },
    "direct_to_patient": {
        "description": "Investigational product shipped directly to the patient's home, "
                       "bypassing site-based dispensing.",
        "regulatory_status": "Requires appropriate IP labeling, temperature monitoring, "
                             "and chain of custody documentation per FDA/EMA requirements.",
        "data_quality_considerations": [
            "Temperature-controlled shipping and monitoring",
            "IP accountability and reconciliation",
            "Compliance tracking (smart packaging, MEMS caps)",
            "State/country pharmacy regulations for DTP",
        ],
        "patient_preference": "Very high — eliminates site visits for drug pickup.",
    },
    "remote_monitoring": {
        "description": "Continuous remote patient monitoring using connected devices (blood pressure "
                       "cuffs, pulse oximeters, ECG patches, glucometers) that transmit data in "
                       "real-time to clinical trial sites and safety monitoring teams.",
        "regulatory_status": "FDA guidance on use of DHTs in clinical investigations (2023). "
                             "Requires device fit-for-purpose validation and data integration.",
        "data_quality_considerations": [
            "Device calibration and accuracy validation",
            "Real-time safety alerting thresholds and escalation pathways",
            "Data transmission reliability and connectivity requirements",
            "Integration with EDC systems and safety databases",
            "Patient training on device use and troubleshooting",
        ],
        "patient_preference": "High — reduces site visits while improving safety monitoring; "
                              "may increase sense of security for patients in remote areas.",
    },
    "digital_informed_consent": {
        "description": "Video-based informed consent process with interactive comprehension checks, "
                       "multi-language support, and dynamic consent updates, going beyond simple "
                       "e-signature to enhance participant understanding.",
        "regulatory_status": "FDA Guidance on Use of Electronic Informed Consent (2016). "
                             "21 CFR Part 11 compliance. IRB review of multimedia consent materials.",
        "data_quality_considerations": [
            "Comprehension assessment quiz with minimum passing threshold",
            "Multi-language video and text with certified translations",
            "Version control and re-consent tracking for protocol amendments",
            "Accessibility features (closed captions, text-to-speech, large font)",
            "Audit trail for consent interactions and comprehension scores",
        ],
        "patient_preference": "High — video explanations improve understanding, especially for "
                              "complex trials; ability to re-review consent materials valued highly.",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 8. LANDMARK TRIALS
# ═══════════════════════════════════════════════════════════════════════════════

LANDMARK_TRIALS: Dict[str, Dict[str, Any]] = {
    "KEYNOTE-024": {
        "name": "KEYNOTE-024",
        "nct_id": "NCT02142738",
        "phase": "Phase 3",
        "indication": "NSCLC (PD-L1 TPS >= 50%)",
        "intervention": "Pembrolizumab vs platinum-based chemotherapy",
        "primary_endpoint": "Progression-free survival (PFS)",
        "key_finding": "Pembrolizumab significantly improved PFS (HR 0.50) and OS; "
                       "established first-line IO monotherapy for high PD-L1 NSCLC.",
        "enrollment": 305,
        "impact": "Changed standard of care for PD-L1-high NSCLC; established biomarker-driven IO.",
    },
    "EMPEROR-Reduced": {
        "name": "EMPEROR-Reduced",
        "nct_id": "NCT03057977",
        "phase": "Phase 3",
        "indication": "Heart failure with reduced ejection fraction (HFrEF)",
        "intervention": "Empagliflozin vs placebo (on top of GDMT)",
        "primary_endpoint": "Composite of CV death or HF hospitalization",
        "key_finding": "Empagliflozin reduced primary endpoint by 25% (HR 0.75); "
                       "SGLT2i benefit extended beyond diabetes to all HFrEF.",
        "enrollment": 3730,
        "impact": "Established SGLT2i as fourth pillar of HF therapy regardless of diabetes status.",
    },
    "RECOVERY": {
        "name": "RECOVERY (Randomised Evaluation of COVid-19 thERapY)",
        "nct_id": "NCT04381936",
        "phase": "Phase 2/3 platform trial",
        "indication": "Hospitalized COVID-19",
        "intervention": "Multiple arms (dexamethasone, tocilizumab, baricitinib, etc.)",
        "primary_endpoint": "28-day all-cause mortality",
        "key_finding": "Dexamethasone reduced mortality by one-third in ventilated patients; "
                       "largest adaptive platform trial in history.",
        "enrollment": 47000,
        "impact": "Demonstrated power of platform trials; dexamethasone became global standard of care.",
    },
    "PARADIGM-HF": {
        "name": "PARADIGM-HF",
        "nct_id": "NCT01035255",
        "phase": "Phase 3",
        "indication": "Heart failure with reduced ejection fraction",
        "intervention": "Sacubitril/valsartan (ARNI) vs enalapril",
        "primary_endpoint": "Composite of CV death or HF hospitalization",
        "key_finding": "ARNI reduced primary endpoint by 20% (HR 0.80) vs ACE inhibitor; "
                       "stopped early for efficacy.",
        "enrollment": 8442,
        "impact": "Established ARNI as standard of care for HFrEF, replacing ACE inhibitors.",
    },
    "CheckMate-067": {
        "name": "CheckMate-067",
        "nct_id": "NCT01844505",
        "phase": "Phase 3",
        "indication": "Advanced melanoma",
        "intervention": "Nivolumab + ipilimumab vs nivolumab vs ipilimumab",
        "primary_endpoint": "PFS and OS (co-primary)",
        "key_finding": "Combination IO achieved 52% 5-year OS (vs 44% nivo, 26% ipi); "
                       "demonstrated durable responses with checkpoint combination.",
        "enrollment": 945,
        "impact": "Established combination checkpoint immunotherapy as transformative approach in oncology.",
    },
    "SPRINT": {
        "name": "SPRINT (Systolic Blood Pressure Intervention Trial)",
        "nct_id": "NCT01206062",
        "phase": "Phase 3",
        "indication": "Hypertension (high cardiovascular risk)",
        "intervention": "Intensive BP control (<120 mmHg) vs standard (<140 mmHg)",
        "primary_endpoint": "Composite MACE (MI, ACS, stroke, HF, CV death)",
        "key_finding": "Intensive control reduced primary endpoint by 25% and all-cause mortality "
                       "by 27%; stopped early for efficacy.",
        "enrollment": 9361,
        "impact": "Changed BP targets in guidelines; influenced ACC/AHA 2017 HTN guidelines.",
    },
    "EMPA-REG_OUTCOME": {
        "name": "EMPA-REG OUTCOME",
        "nct_id": "NCT01131676",
        "phase": "Phase 3",
        "indication": "Type 2 diabetes with established CV disease",
        "intervention": "Empagliflozin vs placebo",
        "primary_endpoint": "3-point MACE",
        "key_finding": "Empagliflozin reduced CV death by 38% and HF hospitalization by 35%; "
                       "first diabetes drug to show CV mortality benefit.",
        "enrollment": 7020,
        "impact": "Launched SGLT2i revolution in cardio-renal-metabolic medicine.",
    },
    "DAPA-CKD": {
        "name": "DAPA-CKD",
        "nct_id": "NCT03036150",
        "phase": "Phase 3",
        "indication": "Chronic kidney disease (with or without diabetes)",
        "intervention": "Dapagliflozin vs placebo",
        "primary_endpoint": "Composite of sustained eGFR decline, ESKD, renal or CV death",
        "key_finding": "Dapagliflozin reduced primary endpoint by 39% (HR 0.61); "
                       "benefit consistent regardless of diabetes status.",
        "enrollment": 4304,
        "impact": "Extended SGLT2i to non-diabetic CKD; changed nephrology practice.",
    },
    "HIMALAYA": {
        "name": "HIMALAYA",
        "nct_id": "NCT03298451",
        "phase": "Phase 3",
        "indication": "Unresectable hepatocellular carcinoma",
        "intervention": "Tremelimumab + durvalumab (STRIDE) vs sorafenib",
        "primary_endpoint": "Overall survival",
        "key_finding": "STRIDE regimen improved OS vs sorafenib (HR 0.78); "
                       "single priming dose of anti-CTLA-4 was sufficient.",
        "enrollment": 1171,
        "impact": "Introduced novel single-tremelimumab priming dose concept in IO combinations.",
    },
    "DESTINY-Breast04": {
        "name": "DESTINY-Breast04",
        "nct_id": "NCT03734029",
        "phase": "Phase 3",
        "indication": "HER2-low metastatic breast cancer",
        "intervention": "Trastuzumab deruxtecan (T-DXd) vs physician's choice chemotherapy",
        "primary_endpoint": "PFS in HR-positive cohort",
        "key_finding": "T-DXd improved PFS (HR 0.51) and OS in HER2-low breast cancer; "
                       "created a new targetable population.",
        "enrollment": 557,
        "impact": "Redefined HER2 as a continuous biomarker; expanded ADC treatment paradigm.",
    },
    "ADVANCE": {
        "name": "ADVANCE (Aducanumab)",
        "nct_id": "NCT02477800",
        "phase": "Phase 3",
        "indication": "Early Alzheimer's disease",
        "intervention": "Aducanumab vs placebo",
        "primary_endpoint": "CDR-SB (Clinical Dementia Rating Sum of Boxes)",
        "key_finding": "Controversial approval based on amyloid reduction as surrogate; "
                       "one trial met primary endpoint, one did not.",
        "enrollment": 3285,
        "impact": "Landmark for accelerated approval debate; spurred amyloid biomarker surrogate discussion.",
    },
    "CLARITY-AD": {
        "name": "CLARITY-AD (Lecanemab)",
        "nct_id": "NCT03887455",
        "phase": "Phase 3",
        "indication": "Early Alzheimer's disease",
        "intervention": "Lecanemab vs placebo",
        "primary_endpoint": "CDR-SB change from baseline at 18 months",
        "key_finding": "Lecanemab reduced clinical decline by 27% on CDR-SB; "
                       "first anti-amyloid to show clear clinical + biomarker benefit.",
        "enrollment": 1795,
        "impact": "First fully approved anti-amyloid therapy; validated amyloid hypothesis clinically.",
    },
    "FOURIER": {
        "name": "FOURIER",
        "nct_id": "NCT01764633",
        "phase": "Phase 3",
        "indication": "Atherosclerotic cardiovascular disease",
        "intervention": "Evolocumab (PCSK9 inhibitor) vs placebo",
        "primary_endpoint": "Composite MACE (CV death, MI, stroke, hospitalization for UA, "
                            "coronary revascularization)",
        "key_finding": "Evolocumab reduced LDL-C by 59% and MACE by 15% (HR 0.85); "
                       "proved PCSK9 inhibition reduces CV events.",
        "enrollment": 27564,
        "impact": "Validated PCSK9 as therapeutic target; established 'lower is better' for LDL-C.",
    },
    "FLAURA": {
        "name": "FLAURA",
        "nct_id": "NCT02296125",
        "phase": "Phase 3",
        "indication": "EGFR-mutated advanced NSCLC (first-line)",
        "intervention": "Osimertinib vs standard EGFR-TKI (gefitinib or erlotinib)",
        "primary_endpoint": "Progression-free survival",
        "key_finding": "Osimertinib nearly doubled PFS (18.9 vs 10.2 months, HR 0.46); "
                       "also improved OS (38.6 vs 31.8 months).",
        "enrollment": 556,
        "impact": "Established osimertinib as first-line standard for EGFR-mutated NSCLC.",
    },
    "ADAURA": {
        "name": "ADAURA",
        "nct_id": "NCT02511106",
        "phase": "Phase 3",
        "indication": "Resected EGFR-mutated NSCLC (adjuvant)",
        "intervention": "Osimertinib vs placebo (adjuvant setting)",
        "primary_endpoint": "Disease-free survival (DFS) in Stage II-IIIA",
        "key_finding": "Osimertinib reduced recurrence risk by 80% (HR 0.20) in Stage II-IIIA; "
                       "unprecedented DFS benefit in adjuvant NSCLC.",
        "enrollment": 682,
        "impact": "First targeted adjuvant therapy for NSCLC; changed post-surgical standard of care.",
    },
    "MAGELLAN": {
        "name": "MAGELLAN (Zolbetuximab)",
        "nct_id": "NCT03504397",
        "phase": "Phase 3",
        "indication": "CLDN18.2-positive gastric/GEJ cancer",
        "intervention": "Zolbetuximab + mFOLFOX6 vs placebo + mFOLFOX6",
        "primary_endpoint": "Progression-free survival and overall survival (co-primary)",
        "key_finding": "Zolbetuximab improved PFS (HR 0.75) and OS (HR 0.77) in "
                       "CLDN18.2-positive gastric cancer.",
        "enrollment": 868,
        "impact": "First CLDN18.2-targeted therapy; new biomarker-selected approach in gastric cancer.",
    },
    "I-SPY_2": {
        "name": "I-SPY 2 (Investigation of Serial Studies to Predict Your Therapeutic Response)",
        "nct_id": "NCT01042379",
        "phase": "Phase 2 adaptive platform",
        "indication": "High-risk early-stage breast cancer (neoadjuvant)",
        "intervention": "Multiple investigational agents + standard chemotherapy backbone",
        "primary_endpoint": "Pathological complete response (pCR)",
        "key_finding": "Graduated 10+ drugs to Phase 3 via Bayesian adaptive design; "
                       "identified pembrolizumab benefit in TNBC neoadjuvant setting.",
        "enrollment": "2000+ (ongoing)",
        "impact": "Gold standard for adaptive platform trial design; accelerated drug development.",
    },
    "VICTORIA": {
        "name": "VICTORIA",
        "nct_id": "NCT02861534",
        "phase": "Phase 3",
        "indication": "Worsening chronic heart failure (HFrEF)",
        "intervention": "Vericiguat vs placebo",
        "primary_endpoint": "Composite of CV death or HF hospitalization",
        "key_finding": "Vericiguat reduced primary endpoint by 10% (HR 0.90); "
                       "modest but significant benefit in high-risk HF population.",
        "enrollment": 5050,
        "impact": "Introduced soluble guanylate cyclase stimulation as new HF mechanism.",
    },
    "TOPAZ-1": {
        "name": "TOPAZ-1",
        "nct_id": "NCT03875235",
        "phase": "Phase 3",
        "indication": "Advanced biliary tract cancer (first-line)",
        "intervention": "Durvalumab + gemcitabine/cisplatin vs placebo + gemcitabine/cisplatin",
        "primary_endpoint": "Overall survival",
        "key_finding": "Durvalumab added to chemotherapy improved OS (HR 0.80, 12.8 vs 11.5 mo); "
                       "first IO to show OS benefit in biliary cancer.",
        "enrollment": 685,
        "impact": "Changed first-line standard of care in biliary tract cancer to include IO.",
    },
    "KRYSTAL-1": {
        "name": "KRYSTAL-1 (Adagrasib)",
        "nct_id": "NCT03785249",
        "phase": "Phase 1/2",
        "indication": "KRAS G12C-mutated NSCLC",
        "intervention": "Adagrasib monotherapy",
        "primary_endpoint": "Objective response rate (ORR)",
        "key_finding": "Adagrasib achieved 42.9% ORR with median DOR 8.5 months; "
                       "confirmed druggability of KRAS G12C.",
        "enrollment": 116,
        "impact": "Second KRAS G12C inhibitor to gain approval; validated target class.",
    },
    "ELARA": {
        "name": "ELARA (Tisagenlecleucel in FL)",
        "nct_id": "NCT03568461",
        "phase": "Phase 2",
        "indication": "Relapsed/refractory follicular lymphoma",
        "intervention": "Tisagenlecleucel (CAR-T, CD19-directed)",
        "primary_endpoint": "Complete response rate (CRR)",
        "key_finding": "CRR of 68% with tisagenlecleucel in r/r FL; "
                       "durable responses with manageable CRS.",
        "enrollment": 97,
        "impact": "Expanded CAR-T therapy to indolent lymphomas.",
    },
    "CREST": {
        "name": "CREST (Semaglutide in Obesity)",
        "nct_id": "NCT03693430",
        "phase": "Phase 3",
        "indication": "Obesity / Overweight with comorbidities",
        "intervention": "Semaglutide 2.4 mg weekly vs placebo",
        "primary_endpoint": "Percent change in body weight at 68 weeks",
        "key_finding": "Semaglutide achieved ~15% body weight reduction vs ~2.4% placebo; "
                       "transformed obesity pharmacotherapy.",
        "enrollment": 1961,
        "impact": "Established GLP-1 RA as highly effective obesity treatment; reshaped metabolic medicine.",
    },
    "CheckMate-227": {
        "name": "CheckMate-227",
        "nct_id": "NCT02477826",
        "phase": "Phase 3",
        "indication": "First-line advanced NSCLC",
        "intervention": "Nivolumab + ipilimumab vs chemotherapy",
        "primary_endpoint": "Overall survival",
        "key_finding": "Nivo+ipi demonstrated OS benefit regardless of PD-L1 expression level; "
                       "4-year OS rate 29% vs 18% with chemotherapy.",
        "enrollment": 1739,
        "impact": "Established dual IO combination as first-line option for NSCLC irrespective of PD-L1 status.",
    },
    "KEYNOTE-522": {
        "name": "KEYNOTE-522",
        "nct_id": "NCT03036488",
        "phase": "Phase 3",
        "indication": "Early-stage triple-negative breast cancer (TNBC)",
        "intervention": "Pembrolizumab + chemotherapy (neoadjuvant) followed by pembrolizumab (adjuvant) vs placebo + chemo",
        "primary_endpoint": "Pathological complete response (pCR) and event-free survival (EFS)",
        "key_finding": "Pembro+chemo achieved pCR 64.8% vs 51.2% placebo; significant EFS improvement (HR 0.63); "
                       "established immunotherapy in early TNBC.",
        "enrollment": 1174,
        "impact": "Changed standard of care for early-stage TNBC to include neoadjuvant/adjuvant pembrolizumab.",
    },
    "DESTINY-Lung02": {
        "name": "DESTINY-Lung02",
        "nct_id": "NCT04644237",
        "phase": "Phase 2",
        "indication": "HER2-mutant advanced NSCLC",
        "intervention": "Trastuzumab deruxtecan (T-DXd) 5.4 mg/kg vs 6.4 mg/kg",
        "primary_endpoint": "Objective response rate (ORR)",
        "key_finding": "T-DXd 5.4 mg/kg achieved ORR 49% in HER2-mutant NSCLC with manageable safety; "
                       "confirmed HER2 mutations as actionable target in lung cancer.",
        "enrollment": 152,
        "impact": "Led to accelerated approval of T-DXd for HER2-mutant NSCLC; expanded ADC paradigm beyond breast cancer.",
    },
    "CLEAR_Outcomes": {
        "name": "CLEAR Outcomes",
        "nct_id": "NCT02993406",
        "phase": "Phase 3",
        "indication": "Statin-intolerant patients at high cardiovascular risk",
        "intervention": "Bempedoic acid vs placebo",
        "primary_endpoint": "4-point MACE (CV death, MI, stroke, coronary revascularization)",
        "key_finding": "Bempedoic acid reduced 4-point MACE by 13% (HR 0.87); first non-statin "
                       "oral therapy to demonstrate CV outcomes benefit.",
        "enrollment": 13970,
        "impact": "Provided evidence-based option for statin-intolerant patients; validated ACL inhibition as CV target.",
    },
    "SELECT": {
        "name": "SELECT (Semaglutide Effects on Cardiovascular Outcomes in People with Overweight or Obesity)",
        "nct_id": "NCT03574597",
        "phase": "Phase 3",
        "indication": "Obesity/overweight with established CV disease (without diabetes)",
        "intervention": "Semaglutide 2.4 mg weekly vs placebo",
        "primary_endpoint": "3-point MACE (CV death, non-fatal MI, non-fatal stroke)",
        "key_finding": "Semaglutide reduced MACE by 20% (HR 0.80) in obese patients without diabetes; "
                       "demonstrated cardiometabolic benefit of GLP-1 RA beyond glucose control.",
        "enrollment": 17604,
        "impact": "Established obesity treatment as cardiovascular risk reduction strategy; blurred line between metabolic and CV medicine.",
    },
    "STEP_HFpEF": {
        "name": "STEP-HFpEF",
        "nct_id": "NCT04788511",
        "phase": "Phase 3",
        "indication": "Heart failure with preserved ejection fraction (HFpEF) and obesity",
        "intervention": "Semaglutide 2.4 mg weekly vs placebo",
        "primary_endpoint": "Dual primary: KCCQ-CSS and body weight at 52 weeks",
        "key_finding": "Semaglutide improved KCCQ-CSS by 7.8 points vs placebo and reduced body weight "
                       "by 13.3% vs 2.6%; reduced HF-related symptoms and physical limitations.",
        "enrollment": 529,
        "impact": "First pharmacotherapy to meaningfully improve HFpEF symptoms via obesity treatment; "
                  "paradigm shift in HFpEF management.",
    },
    "TRAILBLAZER-ALZ_2": {
        "name": "TRAILBLAZER-ALZ 2",
        "nct_id": "NCT04437511",
        "phase": "Phase 3",
        "indication": "Early symptomatic Alzheimer's disease (amyloid and tau positive)",
        "intervention": "Donanemab vs placebo",
        "primary_endpoint": "iADRS (integrated Alzheimer's Disease Rating Scale) change",
        "key_finding": "Donanemab slowed cognitive decline by 35% on iADRS in low/medium tau population; "
                       "76% of patients achieved amyloid negativity; treatment could be discontinued.",
        "enrollment": 1736,
        "impact": "Demonstrated tau-staging approach for patient selection in AD; introduced concept of "
                  "treatment discontinuation upon amyloid clearance.",
    },
    "EMERGE_ENGAGE": {
        "name": "EMERGE / ENGAGE (Aducanumab Phase 3 twins)",
        "nct_id": "NCT02484547",
        "phase": "Phase 3",
        "indication": "Early Alzheimer's disease (MCI or mild dementia due to AD)",
        "intervention": "Aducanumab (high-dose 10 mg/kg) vs placebo",
        "primary_endpoint": "CDR-SB (Clinical Dementia Rating Sum of Boxes) at 78 weeks",
        "key_finding": "EMERGE met primary endpoint (22% slowing on CDR-SB) but ENGAGE did not; "
                       "discordant results led to controversial FDA accelerated approval in 2021.",
        "enrollment": 3285,
        "impact": "Most controversial FDA accelerated approval; three advisory committee members resigned; "
                  "reshaped regulatory debate on surrogate-based approval in neuroscience.",
    },
    "SPRINT_SMA": {
        "name": "SPRINT (SPR1NT) — Pre-symptomatic SMA Gene Therapy",
        "nct_id": "NCT03505099",
        "phase": "Phase 3",
        "indication": "Pre-symptomatic spinal muscular atrophy (SMA) identified by newborn screening",
        "intervention": "Onasemnogene abeparvovec (Zolgensma) single IV infusion",
        "primary_endpoint": "Achievement of developmental motor milestones (sitting, standing, walking)",
        "key_finding": "Pre-symptomatic treatment enabled children with 2 SMN2 copies to sit, stand, "
                       "and walk within normal developmental windows; transformative outcomes.",
        "enrollment": 30,
        "impact": "Landmark gene therapy demonstrating curative potential when administered before symptom onset; "
                  "supported newborn screening programs for SMA worldwide.",
    },
    "SUNFISH": {
        "name": "SUNFISH",
        "nct_id": "NCT02908685",
        "phase": "Phase 2/3",
        "indication": "Spinal muscular atrophy Types 2 and 3 (aged 2-25 years)",
        "intervention": "Risdiplam (oral SMN2 splicing modifier) vs placebo",
        "primary_endpoint": "MFM32 (Motor Function Measure 32) change from baseline at Month 12",
        "key_finding": "Risdiplam significantly improved motor function (MFM32 change 1.36 vs -0.19 placebo); "
                       "first oral treatment for SMA across broad age range.",
        "enrollment": 231,
        "impact": "Established oral therapy option for SMA; expanded access beyond infusion-based treatments; "
                  "demonstrated benefit across SMA types and ages.",
    },
    "CASGEVY": {
        "name": "CASGEVY (Exa-cel / CTX001 — CRISPR Gene Editing)",
        "nct_id": "NCT03655678",
        "phase": "Phase 2/3",
        "indication": "Sickle cell disease and transfusion-dependent beta-thalassemia",
        "intervention": "Exagamglogene autotemcel (exa-cel) — CRISPR/Cas9 edited autologous CD34+ cells",
        "primary_endpoint": "Freedom from vaso-occlusive crises (SCD) or transfusion independence (beta-thal)",
        "key_finding": "97% of SCD patients were free of vaso-occlusive crises for 12+ months; "
                       "93% of beta-thalassemia patients achieved transfusion independence.",
        "enrollment": 97,
        "impact": "First CRISPR-based therapy approved (UK MHRA Dec 2023, FDA Dec 2023); "
                  "historic milestone for gene editing as therapeutic modality.",
    },
    "RINVOQ_SELECT": {
        "name": "SELECT-COMPARE (Upadacitinib)",
        "nct_id": "NCT02629159",
        "phase": "Phase 3",
        "indication": "Moderate-to-severe rheumatoid arthritis (MTX-IR)",
        "intervention": "Upadacitinib 15 mg vs adalimumab 40 mg vs placebo (on MTX background)",
        "primary_endpoint": "ACR20 at Week 12 (vs placebo); DAS28-CRP < 2.6 at Week 26 (vs adalimumab)",
        "key_finding": "Upadacitinib demonstrated superiority over adalimumab on DAS28-CRP remission "
                       "(29% vs 18%) and non-inferiority/superiority on ACR50/70 responses.",
        "enrollment": 1629,
        "impact": "First JAK inhibitor to demonstrate head-to-head superiority over adalimumab in RA; "
                  "challenged biologic-first treatment paradigms.",
    },
    "SURMOUNT-1": {
        "name": "SURMOUNT-1",
        "nct_id": "NCT04184622",
        "phase": "Phase 3",
        "indication": "Obesity (BMI >= 30) or overweight (BMI >= 27) with comorbidity, without diabetes",
        "intervention": "Tirzepatide (5, 10, or 15 mg) vs placebo",
        "primary_endpoint": "Percent change in body weight and proportion achieving >= 5% weight loss at 72 weeks",
        "key_finding": "Tirzepatide 15 mg achieved 22.5% mean weight loss (vs 2.4% placebo); "
                       "63% of patients lost >= 20% body weight; unprecedented efficacy.",
        "enrollment": 2539,
        "impact": "Set new benchmark for pharmacological weight loss; dual GIP/GLP-1 agonism emerged as "
                  "most effective anti-obesity mechanism.",
    },
    "PURPOSE_1": {
        "name": "PURPOSE 1 (Lenacapavir for HIV PrEP)",
        "nct_id": "NCT04994509",
        "phase": "Phase 3",
        "indication": "HIV pre-exposure prophylaxis (PrEP) in cisgender women",
        "intervention": "Lenacapavir (twice-yearly subcutaneous injection) vs daily oral TDF/FTC vs daily oral TAF/FTC",
        "primary_endpoint": "HIV incidence",
        "key_finding": "Lenacapavir achieved 100% efficacy — zero HIV infections among 2134 participants; "
                       "daily oral PrEP showed 83-89% risk reduction vs background incidence.",
        "enrollment": 5338,
        "impact": "Twice-yearly injection for HIV prevention is transformative for adherence; "
                  "potential to reshape global HIV prevention strategy.",
    },
    "PANORAMIC": {
        "name": "PANORAMIC (Platform Adaptive trial of NOvel antiviRals for eArly treatMent of COVID-19 In the Community)",
        "nct_id": "NCT04998396",
        "phase": "Phase 3 (platform trial)",
        "indication": "COVID-19 in the community (high-risk, vaccinated population)",
        "intervention": "Molnupiravir vs usual care",
        "primary_endpoint": "COVID-19 related hospitalization or death within 28 days",
        "key_finding": "Molnupiravir reduced time to first recovery by 4.2 days in vaccinated population; "
                       "hospitalization/death benefit was modest and not statistically significant in primary analysis.",
        "enrollment": 26411,
        "impact": "Demonstrated challenges of antiviral trials in vaccinated populations; "
                  "highlighted importance of community-based platform trial design.",
    },
    "EPIC-HR": {
        "name": "EPIC-HR (Evaluation of Protease Inhibition for COVID-19 in High-Risk patients)",
        "nct_id": "NCT04960202",
        "phase": "Phase 2/3",
        "indication": "Symptomatic COVID-19 in non-hospitalized high-risk adults",
        "intervention": "Nirmatrelvir/ritonavir (Paxlovid) vs placebo",
        "primary_endpoint": "COVID-19 related hospitalization or death through Day 28",
        "key_finding": "Nirmatrelvir/ritonavir reduced hospitalization or death by 89% (HR 0.11) "
                       "when given within 5 days of symptom onset; NNT = 18.",
        "enrollment": 2246,
        "impact": "Established oral protease inhibitor as first-line outpatient COVID-19 treatment; "
                  "received FDA EUA in Dec 2021, full approval in May 2023.",
    },
    "IMpower110": {
        "name": "IMpower110",
        "nct_id": "NCT02409342",
        "phase": "Phase 3",
        "indication": "First-line NSCLC with high PD-L1 expression (TC3/IC3)",
        "intervention": "Atezolizumab monotherapy vs platinum-based chemotherapy",
        "primary_endpoint": "Overall survival in PD-L1-high population",
        "key_finding": "Atezolizumab improved OS by 7.1 months vs chemotherapy (20.2 vs 13.1 mo, HR 0.59) "
                       "in PD-L1-high NSCLC; confirmed IO monotherapy benefit.",
        "enrollment": 572,
        "impact": "Expanded first-line IO monotherapy options beyond pembrolizumab in PD-L1-high NSCLC.",
    },
    "DAPA-HF": {
        "name": "DAPA-HF",
        "nct_id": "NCT03036124",
        "phase": "Phase 3",
        "indication": "Heart failure with reduced ejection fraction (with or without diabetes)",
        "intervention": "Dapagliflozin 10 mg vs placebo (on top of GDMT)",
        "primary_endpoint": "Composite of worsening HF (hospitalization or urgent HF visit) or CV death",
        "key_finding": "Dapagliflozin reduced primary endpoint by 26% (HR 0.74); benefit consistent "
                       "regardless of diabetes status; stopped early for overwhelming efficacy.",
        "enrollment": 4744,
        "impact": "First SGLT2 inhibitor to show HF benefit independent of diabetes; established new pillar of HF therapy.",
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# 9. SAFETY SIGNAL METRICS
# ═══════════════════════════════════════════════════════════════════════════════

SAFETY_SIGNAL_METRICS: Dict[str, Dict[str, Any]] = {
    "PRR": {
        "name": "Proportional Reporting Ratio",
        "description": "Compares the proportion of a specific adverse event for a drug "
                       "against the proportion of that event for all other drugs in the database.",
        "formula": "PRR = (a / (a + b)) / (c / (c + d)), where a = target drug-event pair, "
                   "b = target drug all other events, c = all other drugs same event, "
                   "d = all other drugs all other events.",
        "thresholds": {
            "signal": "PRR >= 2",
            "chi_squared": ">= 4",
            "n_cases": ">= 3",
        },
        "interpretation": "PRR >= 2 with chi-squared >= 4 and N >= 3 suggests a signal. "
                          "Higher PRR indicates stronger disproportionality.",
    },
    "ROR": {
        "name": "Reporting Odds Ratio",
        "description": "Odds ratio comparing the odds of a specific event being reported "
                       "for the drug vs. all other drugs.",
        "formula": "ROR = (a * d) / (b * c), where a/b/c/d as in PRR 2x2 table.",
        "thresholds": {
            "signal": "Lower bound of 95% CI > 1",
            "n_cases": ">= 3",
        },
        "interpretation": "ROR with lower 95% CI > 1 is considered a signal. "
                          "More conservative than PRR as it uses confidence intervals.",
    },
    "IC": {
        "name": "Information Component (Bayesian)",
        "description": "Bayesian measure of disproportionality based on the observed-to-expected "
                       "ratio of drug-event combinations. Used by WHO-UMC (Uppsala).",
        "formula": "IC = log2(observed / expected), where expected = (Ndrug * Nevent) / Ntotal. "
                   "IC025 is the lower bound of the 95% credibility interval.",
        "thresholds": {
            "signal": "IC025 > 0",
            "strong_signal": "IC025 > 1.5",
        },
        "interpretation": "IC025 > 0 indicates the drug-event combination is reported more "
                          "frequently than expected. Bayesian shrinkage reduces false positives.",
    },
    "EBGM": {
        "name": "Empirical Bayes Geometric Mean",
        "description": "Multi-item gamma-Poisson shrinkage (MGPS) estimator used by FDA FAERS. "
                       "Shrinks observed/expected ratios toward the null to reduce noise.",
        "formula": "EBGM = geometric mean of the posterior distribution of lambda (true RR). "
                   "EB05 = lower 5th percentile of the posterior.",
        "thresholds": {
            "signal": "EB05 >= 2",
            "n_cases": ">= 3",
        },
        "interpretation": "EB05 >= 2 is the standard FDA signal threshold. "
                          "EBGM provides point estimate; EB05 provides conservative lower bound.",
    },
    "MGPS": {
        "name": "Multi-item Gamma-Poisson Shrinker",
        "description": "Bayesian shrinkage estimator for pharmacovigilance signal detection that "
                       "models the expected count of drug-event pairs using a gamma-Poisson "
                       "mixture model. Shrinks observed-to-expected ratios toward the overall "
                       "mean, reducing false signals from rare drug-event combinations.",
        "formula": "Lambda_ij ~ Gamma(alpha, beta) mixture; EBGM = geometric mean of posterior. "
                   "Shrinkage is strongest for rare combinations with few reports.",
        "thresholds": {
            "signal": "EB05 >= 2 (lower 5th percentile of posterior)",
            "n_cases": ">= 1 (but shrinkage handles rare events)",
        },
        "interpretation": "MGPS is the method underlying the EBGM metric used by FDA FAERS. "
                          "Provides robust signal detection by borrowing strength across the "
                          "entire database, reducing false positives from random variation.",
    },
    "TreeScan": {
        "name": "Tree-based Scan Statistic",
        "description": "Hierarchical signal detection method that evaluates adverse event signals "
                       "across the MedDRA hierarchy (SOC -> HLGT -> HLT -> PT) simultaneously, "
                       "adjusting for multiple testing at all levels of the tree structure.",
        "formula": "Unconditional Bernoulli or Poisson model; likelihood ratio test statistic "
                   "at each node in the MedDRA tree. P-values via Monte Carlo simulation.",
        "thresholds": {
            "signal": "P-value < 0.05 after tree-adjusted multiple testing correction",
            "n_cases": ">= 2 at the node level",
        },
        "interpretation": "TreeScan detects signals at the most informative level of the MedDRA "
                          "hierarchy — from broad SOC-level patterns to specific PT-level signals. "
                          "Adjusts for multiple comparisons across the entire tree, reducing "
                          "false positives while maintaining power to detect hierarchical clusters.",
    },
}


# ===================================================================
# PEDIATRIC CLINICAL TRIAL KNOWLEDGE
# ===================================================================

PEDIATRIC_TRIAL_KNOWLEDGE: Dict[str, Dict[str, Any]] = {
    "age_filters": {
        "description": "ClinicalTrials.gov standard age group filters for pediatric trials",
        "categories": {
            "child": "Birth to 17 years",
            "adolescent": "13 to 17 years",
            "adult": "18 to 64 years",
        },
        "note": "Pediatric trials may span multiple age groups; neonates (0-27 days) and "
               "infants (28 days-23 months) are sometimes separated.",
    },
    "cooperative_groups": {
        "COG": {
            "name": "Children's Oncology Group",
            "description": "Primary pediatric cooperative group in North America; >200 institutions",
            "landmark_trials": {
                "AALL0932": {
                    "indication": "Standard-risk B-cell ALL",
                    "design": "Risk-stratified therapy based on MRD and cytogenetics",
                    "outcome": "5-year EFS >90%",
                },
                "AALL1231": {
                    "indication": "High-risk B-cell ALL",
                    "design": "MRD-directed therapy with augmented intensification",
                    "outcome": "Improved EFS with MRD-guided therapy escalation",
                },
                "ANBL1232": {
                    "indication": "High-risk neuroblastoma",
                    "design": "Immunotherapy (anti-GD2 + cytokines) in consolidation",
                    "outcome": "Improved EFS with immunotherapy backbone",
                },
            },
        },
    },
    "precision_oncology_programs": {
        "Pediatric_MATCH": {
            "name": "NCI-COG Pediatric MATCH",
            "description": "Molecularly targeted therapy matching for relapsed/refractory "
                          "pediatric solid tumors and lymphomas",
            "design": "Histology-agnostic, biomarker-driven basket trial",
            "key_feature": "Centralized molecular profiling with multiple treatment arms",
        },
        "INFORM": {
            "name": "INdividualized therapy FOr Relapsed Malignancies in childhood",
            "description": "International pediatric precision oncology registry "
                          "(European-led, multicenter)",
            "design": "Molecular profiling registry with treatment recommendations",
            "key_feature": "WGS, RNA-seq, and methylation profiling for each patient",
        },
    },
    "pediatric_trial_design": {
        "phase1_design": {
            "preferred_method": "Rolling-6 design (not traditional 3+3)",
            "rationale": "Faster enrollment in small patient populations; allows up to 6 "
                        "patients per dose level with continuous enrollment",
        },
        "endpoints": {
            "preferred_primary": "Event-free survival (EFS)",
            "rationale": "EFS preferred over OS in pediatric oncology due to long natural "
                        "survival post-relapse with salvage therapy, and ethical concerns "
                        "about waiting for OS maturity in children",
        },
        "consent_requirements": {
            "parental_consent": "Required for all children <18",
            "child_assent": "Required for children age 7 and older (per AAP/FDA guidance)",
            "adolescent_considerations": "Adolescents 14-17 may provide enhanced assent or "
                                        "co-consent depending on jurisdiction",
        },
        "long_term_follow_up": {
            "requirement": "Required for all pediatric clinical trials",
            "guideline": "COG Long-Term Follow-Up (LTFU) Guidelines",
            "duration": "Minimum 10 years; ideally lifetime monitoring for late effects",
        },
    },
}
