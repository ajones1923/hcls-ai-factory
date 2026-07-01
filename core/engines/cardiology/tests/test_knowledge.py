"""Unit tests for Cardiology Intelligence Agent knowledge graph.

Tests all knowledge dictionaries in src/knowledge.py including
CARDIOVASCULAR_CONDITIONS, CARDIAC_BIOMARKERS, CARDIAC_DRUG_CLASSES,
CARDIOVASCULAR_GENES, IMAGING_MODALITIES, GUIDELINE_RECOMMENDATIONS,
ENTITY_ALIASES, CONDITION_DRUG_MAP, and IMAGING_TRIGGER_MAP.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.knowledge import (
    CARDIAC_BIOMARKERS,
    CARDIAC_DRUG_CLASSES,
    CARDIOVASCULAR_CONDITIONS,
    CARDIOVASCULAR_GENES,
    CONDITION_DRUG_MAP,
    ENTITY_ALIASES,
    GUIDELINE_RECOMMENDATIONS,
    IMAGING_MODALITIES,
    IMAGING_TRIGGER_MAP,
)


# ===================================================================
# CARDIOVASCULAR_CONDITIONS (32 entries)
# ===================================================================


class TestCardiovascularConditions:
    """Tests for the CARDIOVASCULAR_CONDITIONS dictionary."""

    def test_total_count(self):
        assert len(CARDIOVASCULAR_CONDITIONS) == 45

    REQUIRED_KEYS = ["aliases", "icd10", "diagnostic_criteria", "treatment", "genes", "guidelines"]

    @pytest.mark.parametrize("key", REQUIRED_KEYS)
    def test_all_conditions_have_required_key(self, key):
        for condition_name, data in CARDIOVASCULAR_CONDITIONS.items():
            assert key in data, f"'{condition_name}' missing required key '{key}'"

    @pytest.mark.parametrize(
        "condition",
        [
            "Hypertrophic Cardiomyopathy",
            "Dilated Cardiomyopathy",
            "Arrhythmogenic Cardiomyopathy",
            "Coronary Artery Disease",
            "Acute Coronary Syndrome",
            "Stable Angina",
            "Heart Failure with Reduced EF",
            "Heart Failure with Preserved EF",
            "Atrial Fibrillation",
            "Atrial Flutter",
            "Ventricular Tachycardia",
            "Long QT Syndrome",
            "Brugada Syndrome",
            "CPVT",
            "Aortic Stenosis",
            "Aortic Regurgitation",
            "Mitral Regurgitation",
            "Mitral Stenosis",
            "Tricuspid Regurgitation",
            "Pulmonary Hypertension",
            "Aortic Dissection",
            "Thoracic Aortic Aneurysm",
            "Familial Hypercholesterolemia",
            "Cardiac Amyloidosis",
            "Myocarditis",
            "Pericarditis",
            "Infective Endocarditis",
            "Cardiac Sarcoidosis",
            "Takotsubo Cardiomyopathy",
            "Wolff-Parkinson-White Syndrome",
            "Cardiac Tamponade",
            "Restrictive Cardiomyopathy",
        ],
    )
    def test_condition_exists(self, condition):
        assert condition in CARDIOVASCULAR_CONDITIONS

    # ---- Specific condition data tests ----

    def test_hcm_genes_include_myh7(self):
        genes = CARDIOVASCULAR_CONDITIONS["Hypertrophic Cardiomyopathy"]["genes"]
        assert "MYH7" in genes

    def test_hcm_genes_include_mybpc3(self):
        genes = CARDIOVASCULAR_CONDITIONS["Hypertrophic Cardiomyopathy"]["genes"]
        assert "MYBPC3" in genes

    def test_af_icd10(self):
        assert CARDIOVASCULAR_CONDITIONS["Atrial Fibrillation"]["icd10"] == "I48.91"

    def test_af_aliases_contain_afib(self):
        aliases = CARDIOVASCULAR_CONDITIONS["Atrial Fibrillation"]["aliases"]
        assert "AFib" in aliases

    def test_cad_icd10(self):
        assert CARDIOVASCULAR_CONDITIONS["Coronary Artery Disease"]["icd10"] == "I25.1"

    def test_hfref_icd10(self):
        assert CARDIOVASCULAR_CONDITIONS["Heart Failure with Reduced EF"]["icd10"] == "I50.2"

    def test_hfpef_icd10(self):
        assert CARDIOVASCULAR_CONDITIONS["Heart Failure with Preserved EF"]["icd10"] == "I50.3"

    def test_dcm_genes_include_ttn(self):
        genes = CARDIOVASCULAR_CONDITIONS["Dilated Cardiomyopathy"]["genes"]
        assert "TTN" in genes

    def test_dcm_genes_include_lmna(self):
        genes = CARDIOVASCULAR_CONDITIONS["Dilated Cardiomyopathy"]["genes"]
        assert "LMNA" in genes

    def test_arvc_genes_include_pkp2(self):
        genes = CARDIOVASCULAR_CONDITIONS["Arrhythmogenic Cardiomyopathy"]["genes"]
        assert "PKP2" in genes

    def test_lqts_genes_include_kcnq1(self):
        genes = CARDIOVASCULAR_CONDITIONS["Long QT Syndrome"]["genes"]
        assert "KCNQ1" in genes

    def test_lqts_genes_include_scn5a(self):
        genes = CARDIOVASCULAR_CONDITIONS["Long QT Syndrome"]["genes"]
        assert "SCN5A" in genes

    def test_brugada_genes_include_scn5a(self):
        genes = CARDIOVASCULAR_CONDITIONS["Brugada Syndrome"]["genes"]
        assert "SCN5A" in genes

    def test_fh_genes_include_ldlr(self):
        genes = CARDIOVASCULAR_CONDITIONS["Familial Hypercholesterolemia"]["genes"]
        assert "LDLR" in genes

    def test_fh_genes_include_pcsk9(self):
        genes = CARDIOVASCULAR_CONDITIONS["Familial Hypercholesterolemia"]["genes"]
        assert "PCSK9" in genes

    def test_cardiac_amyloidosis_gene_ttr(self):
        genes = CARDIOVASCULAR_CONDITIONS["Cardiac Amyloidosis"]["genes"]
        assert "TTR" in genes

    def test_hcm_has_cross_modal_trigger(self):
        assert CARDIOVASCULAR_CONDITIONS["Hypertrophic Cardiomyopathy"]["cross_modal_trigger"] is True

    def test_af_no_cross_modal_trigger(self):
        assert CARDIOVASCULAR_CONDITIONS["Atrial Fibrillation"]["cross_modal_trigger"] is False

    def test_all_conditions_have_aliases_list(self):
        for name, data in CARDIOVASCULAR_CONDITIONS.items():
            assert isinstance(data["aliases"], list), f"'{name}' aliases is not a list"

    def test_all_conditions_have_genes_list(self):
        for name, data in CARDIOVASCULAR_CONDITIONS.items():
            assert isinstance(data["genes"], list), f"'{name}' genes is not a list"

    def test_all_conditions_have_guidelines_list(self):
        for name, data in CARDIOVASCULAR_CONDITIONS.items():
            assert isinstance(data["guidelines"], list), f"'{name}' guidelines is not a list"

    def test_aortic_stenosis_icd10(self):
        assert CARDIOVASCULAR_CONDITIONS["Aortic Stenosis"]["icd10"] == "I35.0"

    def test_pulmonary_hypertension_genes(self):
        genes = CARDIOVASCULAR_CONDITIONS["Pulmonary Hypertension"]["genes"]
        assert "BMPR2" in genes

    def test_aortic_dissection_has_cross_modal(self):
        assert CARDIOVASCULAR_CONDITIONS["Aortic Dissection"]["cross_modal_trigger"] is True


# ===================================================================
# CARDIAC_BIOMARKERS (21 entries)
# ===================================================================


class TestCardiacBiomarkers:
    """Tests for the CARDIAC_BIOMARKERS dictionary."""

    def test_total_count(self):
        assert len(CARDIAC_BIOMARKERS) == 29

    @pytest.mark.parametrize(
        "biomarker",
        [
            "hs-cTnI", "hs-cTnT", "NT-proBNP", "BNP", "CK-MB",
            "D-dimer", "hsCRP", "LDL-C", "HDL-C", "Total Cholesterol",
            "Triglycerides", "Lp(a)", "ApoB", "HbA1c", "Creatinine/eGFR",
            "Potassium", "Magnesium", "Ferritin/Iron/TSAT", "Lactate",
            "Procalcitonin", "IL-6",
        ],
    )
    def test_biomarker_exists(self, biomarker):
        assert biomarker in CARDIAC_BIOMARKERS

    def test_troponin_i_has_full_name(self):
        assert CARDIAC_BIOMARKERS["hs-cTnI"]["full_name"] == "High-sensitivity cardiac troponin I"

    def test_bnp_has_clinical_use(self):
        assert "clinical_use" in CARDIAC_BIOMARKERS["BNP"]

    def test_nt_probnp_has_reference_ranges(self):
        assert "reference_ranges" in CARDIAC_BIOMARKERS["NT-proBNP"]

    def test_all_biomarkers_have_full_name(self):
        for name, data in CARDIAC_BIOMARKERS.items():
            assert "full_name" in data, f"'{name}' missing 'full_name'"

    def test_all_biomarkers_have_clinical_use(self):
        for name, data in CARDIAC_BIOMARKERS.items():
            assert "clinical_use" in data, f"'{name}' missing 'clinical_use'"

    def test_potassium_hyperkalemia_threshold(self):
        ranges = CARDIAC_BIOMARKERS["Potassium"]["reference_ranges"]
        assert "hyperkalemia" in ranges

    def test_ldl_c_has_guideline_thresholds(self):
        assert "guideline_thresholds" in CARDIAC_BIOMARKERS["LDL-C"]


# ===================================================================
# CARDIAC_DRUG_CLASSES (26 entries)
# ===================================================================


class TestCardiacDrugClasses:
    """Tests for the CARDIAC_DRUG_CLASSES dictionary."""

    def test_total_count(self):
        assert len(CARDIAC_DRUG_CLASSES) == 32

    @pytest.mark.parametrize(
        "drug_class",
        [
            "beta_blockers", "ace_inhibitors", "arbs", "arni",
            "sglt2_inhibitors", "mra", "loop_diuretics", "thiazide_diuretics",
            "ccb_dihydropyridine", "ccb_nondihydropyridine", "statins",
            "pcsk9_inhibitors", "ezetimibe", "antiplatelets", "doacs",
            "warfarin", "heparin_lmwh", "antiarrhythmics_class_i",
            "antiarrhythmics_class_iii", "digoxin", "nitrates",
            "hydralazine_isdn", "ivabradine", "mavacamten",
            "glp1_receptor_agonists", "inotropes",
        ],
    )
    def test_drug_class_exists(self, drug_class):
        assert drug_class in CARDIAC_DRUG_CLASSES

    def test_all_classes_have_drugs_list(self):
        for name, data in CARDIAC_DRUG_CLASSES.items():
            assert "drugs" in data, f"'{name}' missing 'drugs'"
            assert isinstance(data["drugs"], list)

    def test_all_classes_have_mechanism(self):
        for name, data in CARDIAC_DRUG_CLASSES.items():
            assert "mechanism" in data, f"'{name}' missing 'mechanism'"

    def test_all_classes_have_indications(self):
        for name, data in CARDIAC_DRUG_CLASSES.items():
            assert "indications" in data, f"'{name}' missing 'indications'"

    def test_beta_blockers_include_carvedilol(self):
        assert "carvedilol" in CARDIAC_DRUG_CLASSES["beta_blockers"]["drugs"]

    def test_beta_blockers_include_metoprolol(self):
        assert "metoprolol succinate" in CARDIAC_DRUG_CLASSES["beta_blockers"]["drugs"]

    def test_sglt2_inhibitors_include_dapagliflozin(self):
        assert "dapagliflozin" in CARDIAC_DRUG_CLASSES["sglt2_inhibitors"]["drugs"]

    def test_arni_includes_sacubitril_valsartan(self):
        assert "sacubitril/valsartan" in CARDIAC_DRUG_CLASSES["arni"]["drugs"]

    def test_statins_include_atorvastatin(self):
        assert "atorvastatin" in CARDIAC_DRUG_CLASSES["statins"]["drugs"]

    def test_doacs_include_apixaban(self):
        assert "apixaban" in CARDIAC_DRUG_CLASSES["doacs"]["drugs"]

    def test_mra_includes_spironolactone(self):
        assert "spironolactone" in CARDIAC_DRUG_CLASSES["mra"]["drugs"]

    def test_beta_blockers_target_dose_carvedilol(self):
        target_doses = CARDIAC_DRUG_CLASSES["beta_blockers"]["target_doses"]
        assert "carvedilol" in target_doses
        assert "25 mg BID" in target_doses["carvedilol"]

    def test_sglt2_target_dose_dapagliflozin(self):
        target_doses = CARDIAC_DRUG_CLASSES["sglt2_inhibitors"]["target_doses"]
        assert target_doses["dapagliflozin"] == "10 mg daily"


# ===================================================================
# CARDIOVASCULAR_GENES (42 entries)
# ===================================================================


class TestCardiovascularGenes:
    """Tests for the CARDIOVASCULAR_GENES dictionary."""

    def test_total_count(self):
        assert len(CARDIOVASCULAR_GENES) >= 41

    @pytest.mark.parametrize(
        "gene",
        [
            "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1",
            "MYL2", "MYL3", "TTN", "LMNA", "RBM20",
            "PKP2", "DSP", "DSG2", "DSC2", "JUP", "TMEM43",
            "SCN5A", "KCNQ1", "KCNH2",
        ],
    )
    def test_gene_exists(self, gene):
        assert gene in CARDIOVASCULAR_GENES

    def test_all_genes_have_full_name(self):
        for gene, data in CARDIOVASCULAR_GENES.items():
            assert "full_name" in data, f"'{gene}' missing 'full_name'"

    def test_all_genes_have_chromosome(self):
        for gene, data in CARDIOVASCULAR_GENES.items():
            assert "chromosome" in data, f"'{gene}' missing 'chromosome'"

    def test_all_genes_have_associated_conditions(self):
        for gene, data in CARDIOVASCULAR_GENES.items():
            assert "associated_conditions" in data, f"'{gene}' missing 'associated_conditions'"

    def test_myh7_conditions_include_hcm(self):
        assert "HCM" in CARDIOVASCULAR_GENES["MYH7"]["associated_conditions"]

    def test_ttn_conditions_include_dcm(self):
        assert "DCM" in CARDIOVASCULAR_GENES["TTN"]["associated_conditions"]

    def test_pkp2_conditions_include_arvc(self):
        assert "ARVC" in CARDIOVASCULAR_GENES["PKP2"]["associated_conditions"]

    def test_scn5a_conditions_include_brugada(self):
        conditions = CARDIOVASCULAR_GENES["SCN5A"]["associated_conditions"]
        assert "Brugada syndrome" in conditions

    def test_kcnq1_conditions_include_lqts(self):
        conditions = CARDIOVASCULAR_GENES["KCNQ1"]["associated_conditions"]
        assert "LQTS type 1" in conditions

    def test_myh7_chromosome(self):
        assert CARDIOVASCULAR_GENES["MYH7"]["chromosome"] == "14q11.2"

    def test_ttn_chromosome(self):
        assert CARDIOVASCULAR_GENES["TTN"]["chromosome"] == "2q31.2"

    def test_lmna_inheritance(self):
        assert CARDIOVASCULAR_GENES["LMNA"]["inheritance"] == "AD"

    def test_genes_have_genetic_testing_indication(self):
        for gene, data in CARDIOVASCULAR_GENES.items():
            assert "genetic_testing_indication" in data, f"'{gene}' missing 'genetic_testing_indication'"


# ===================================================================
# IMAGING_MODALITIES (15 entries)
# ===================================================================


class TestImagingModalities:
    """Tests for the IMAGING_MODALITIES dictionary."""

    def test_total_count(self):
        assert len(IMAGING_MODALITIES) == 15

    @pytest.mark.parametrize(
        "modality",
        [
            "TTE", "TEE", "Stress Echo", "Cardiac CT Calcium Score",
            "Coronary CTA", "Cardiac MRI", "Cardiac PET", "SPECT MPI",
            "MUGA", "Right Heart Catheterization", "Coronary Angiography",
            "12-Lead ECG", "Holter Monitor", "Event Monitor",
            "Device Interrogation",
        ],
    )
    def test_modality_exists(self, modality):
        assert modality in IMAGING_MODALITIES

    def test_all_modalities_have_full_name(self):
        for name, data in IMAGING_MODALITIES.items():
            assert "full_name" in data, f"'{name}' missing 'full_name'"

    def test_all_modalities_have_indications(self):
        for name, data in IMAGING_MODALITIES.items():
            assert "indications" in data, f"'{name}' missing 'indications'"

    def test_all_modalities_have_guideline_society(self):
        for name, data in IMAGING_MODALITIES.items():
            assert "guideline_society" in data, f"'{name}' missing 'guideline_society'"

    def test_tte_full_name(self):
        assert IMAGING_MODALITIES["TTE"]["full_name"] == "Transthoracic Echocardiography"

    def test_cardiac_mri_protocols(self):
        protocols = IMAGING_MODALITIES["Cardiac MRI"]["protocols"]
        assert isinstance(protocols, list)
        assert len(protocols) > 0

    def test_ecg_guideline_society(self):
        assert IMAGING_MODALITIES["12-Lead ECG"]["guideline_society"] == "ACC/AHA/HRS"

    def test_rhc_key_measurements(self):
        measurements = IMAGING_MODALITIES["Right Heart Catheterization"]["key_measurements"]
        assert "PCWP" in measurements
        assert "CO" in measurements


# ===================================================================
# GUIDELINE_RECOMMENDATIONS (52+ entries)
# ===================================================================


class TestGuidelineRecommendations:
    """Tests for the GUIDELINE_RECOMMENDATIONS list."""

    def test_minimum_count(self):
        assert len(GUIDELINE_RECOMMENDATIONS) >= 50

    def test_is_list(self):
        assert isinstance(GUIDELINE_RECOMMENDATIONS, list)

    REQUIRED_KEYS = ["society", "year", "recommendation", "class", "level", "condition", "source"]

    @pytest.mark.parametrize("key", REQUIRED_KEYS)
    def test_all_entries_have_required_key(self, key):
        for idx, rec in enumerate(GUIDELINE_RECOMMENDATIONS):
            assert key in rec, f"Entry {idx} missing required key '{key}'"

    def test_hfref_arni_recommendation_exists(self):
        found = any(
            "ARNI" in r["recommendation"] and r["condition"] == "HFrEF"
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found, "No ARNI recommendation for HFrEF found"

    def test_hfref_sglt2i_recommendation_exists(self):
        found = any(
            "SGLT2" in r["recommendation"] and r["condition"] == "HFrEF"
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found, "No SGLT2i recommendation for HFrEF found"

    def test_hfref_beta_blocker_recommendation_exists(self):
        found = any(
            "beta-blocker" in r["recommendation"].lower() and r["condition"] == "HFrEF"
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found, "No beta-blocker recommendation for HFrEF found"

    def test_af_anticoagulation_recommendation(self):
        found = any(
            "anticoagulation" in r["recommendation"].lower() and "Atrial Fibrillation" in r["condition"]
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found

    def test_doac_over_warfarin_recommendation(self):
        found = any(
            "DOAC" in r["recommendation"] and "warfarin" in r["recommendation"].lower()
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found

    def test_stemi_pci_recommendation(self):
        found = any(
            r["condition"] == "STEMI" and "angiography" in r["recommendation"].lower()
            for r in GUIDELINE_RECOMMENDATIONS
        )
        assert found

    def test_valid_class_values(self):
        valid_classes = {"I", "IIa", "IIb", "III", "IIa-B"}
        for idx, rec in enumerate(GUIDELINE_RECOMMENDATIONS):
            assert rec["class"] in valid_classes, (
                f"Entry {idx} has invalid class '{rec['class']}'"
            )

    def test_valid_level_values(self):
        valid_levels = {"A", "B", "B-R", "B-NR", "C", "C-LD", "C-EO"}
        for idx, rec in enumerate(GUIDELINE_RECOMMENDATIONS):
            assert rec["level"] in valid_levels, (
                f"Entry {idx} has invalid level '{rec['level']}'"
            )

    def test_societies_include_acc_aha(self):
        societies = {r["society"] for r in GUIDELINE_RECOMMENDATIONS}
        assert "ACC/AHA" in societies

    def test_societies_include_esc(self):
        societies = {r["society"] for r in GUIDELINE_RECOMMENDATIONS}
        assert "ESC" in societies


# ===================================================================
# ENTITY_ALIASES (120+ entries)
# ===================================================================


class TestEntityAliases:
    """Tests for the ENTITY_ALIASES dictionary."""

    def test_minimum_count(self):
        assert len(ENTITY_ALIASES) >= 120

    def test_is_dict(self):
        assert isinstance(ENTITY_ALIASES, dict)

    @pytest.mark.parametrize(
        "abbrev, expansion",
        [
            ("MI", "myocardial infarction"),
            ("STEMI", "ST-elevation myocardial infarction"),
            ("ACS", "acute coronary syndrome"),
            ("CAD", "coronary artery disease"),
            ("HF", "heart failure"),
            ("HFrEF", "heart failure with reduced ejection fraction"),
            ("HFpEF", "heart failure with preserved ejection fraction"),
            ("AF", "atrial fibrillation"),
            ("VT", "ventricular tachycardia"),
            ("LQTS", "long QT syndrome"),
            ("HCM", "hypertrophic cardiomyopathy"),
            ("DCM", "dilated cardiomyopathy"),
            ("ARVC", "arrhythmogenic right ventricular cardiomyopathy"),
            ("GDMT", "guideline-directed medical therapy"),
            ("ARNI", "angiotensin receptor-neprilysin inhibitor"),
            ("MRA", "mineralocorticoid receptor antagonist"),
            ("SGLT2i", "SGLT2 inhibitor"),
            ("PCI", "percutaneous coronary intervention"),
            ("CABG", "coronary artery bypass graft"),
            ("ICD", "implantable cardioverter-defibrillator"),
            ("CRT", "cardiac resynchronization therapy"),
            ("ECG", "electrocardiogram"),
            ("GLS", "global longitudinal strain"),
            ("MACE", "major adverse cardiovascular events"),
        ],
    )
    def test_alias_mapping(self, abbrev, expansion):
        assert ENTITY_ALIASES[abbrev] == expansion

    def test_all_values_are_strings(self):
        for key, val in ENTITY_ALIASES.items():
            assert isinstance(val, str), f"Alias '{key}' value is not a string"

    def test_all_keys_are_strings(self):
        for key in ENTITY_ALIASES:
            assert isinstance(key, str)


# ===================================================================
# CONDITION_DRUG_MAP
# ===================================================================


class TestConditionDrugMap:
    """Tests for the CONDITION_DRUG_MAP dictionary."""

    def test_is_dict(self):
        assert isinstance(CONDITION_DRUG_MAP, dict)

    def test_hfref_has_four_pillars(self):
        hfref_drugs = CONDITION_DRUG_MAP["HFrEF"]
        assert "beta_blockers" in hfref_drugs
        assert "arni" in hfref_drugs
        assert "mra" in hfref_drugs
        assert "sglt2_inhibitors" in hfref_drugs

    def test_hfpef_has_sglt2i(self):
        assert "sglt2_inhibitors" in CONDITION_DRUG_MAP["HFpEF"]

    def test_af_has_doacs(self):
        assert "doacs" in CONDITION_DRUG_MAP["Atrial Fibrillation"]

    def test_cad_has_statins(self):
        assert "statins" in CONDITION_DRUG_MAP["Coronary Artery Disease"]

    def test_acs_has_antiplatelets(self):
        assert "antiplatelets" in CONDITION_DRUG_MAP["Acute Coronary Syndrome"]

    def test_hcm_has_beta_blockers(self):
        assert "beta_blockers" in CONDITION_DRUG_MAP["Hypertrophic Cardiomyopathy"]

    def test_fh_has_pcsk9_inhibitors(self):
        assert "pcsk9_inhibitors" in CONDITION_DRUG_MAP["Familial Hypercholesterolemia"]

    def test_all_drug_class_refs_valid(self):
        valid_classes = set(CARDIAC_DRUG_CLASSES.keys())
        for condition, drugs in CONDITION_DRUG_MAP.items():
            for drug in drugs:
                assert drug in valid_classes, (
                    f"CONDITION_DRUG_MAP['{condition}'] references unknown drug class '{drug}'"
                )


# ===================================================================
# IMAGING_TRIGGER_MAP
# ===================================================================


class TestImagingTriggerMap:
    """Tests for the IMAGING_TRIGGER_MAP dictionary."""

    def test_is_dict(self):
        assert isinstance(IMAGING_TRIGGER_MAP, dict)

    def test_has_entries(self):
        assert len(IMAGING_TRIGGER_MAP) >= 6

    @pytest.mark.parametrize(
        "trigger_key",
        [
            "unexplained_lvh_15mm",
            "non_ischemic_lge_dcm",
            "rv_dilation_with_arrhythmia",
            "diffuse_lge_low_voltage",
            "aortic_root_dilation",
            "premature_coronary_calcification",
        ],
    )
    def test_trigger_exists(self, trigger_key):
        assert trigger_key in IMAGING_TRIGGER_MAP

    def test_all_triggers_have_gene_panel(self):
        for key, data in IMAGING_TRIGGER_MAP.items():
            assert "gene_panel" in data, f"'{key}' missing 'gene_panel'"
            assert isinstance(data["gene_panel"], list)

    def test_all_triggers_have_conditions(self):
        for key, data in IMAGING_TRIGGER_MAP.items():
            assert "conditions" in data, f"'{key}' missing 'conditions'"

    def test_all_triggers_have_imaging_finding(self):
        for key, data in IMAGING_TRIGGER_MAP.items():
            assert "imaging_finding" in data, f"'{key}' missing 'imaging_finding'"

    def test_lvh_trigger_genes_include_myh7(self):
        genes = IMAGING_TRIGGER_MAP["unexplained_lvh_15mm"]["gene_panel"]
        assert "MYH7" in genes

    def test_lvh_trigger_genes_include_mybpc3(self):
        genes = IMAGING_TRIGGER_MAP["unexplained_lvh_15mm"]["gene_panel"]
        assert "MYBPC3" in genes

    def test_dcm_trigger_genes_include_ttn(self):
        genes = IMAGING_TRIGGER_MAP["non_ischemic_lge_dcm"]["gene_panel"]
        assert "TTN" in genes

    def test_arvc_trigger_genes_include_pkp2(self):
        genes = IMAGING_TRIGGER_MAP["rv_dilation_with_arrhythmia"]["gene_panel"]
        assert "PKP2" in genes

    def test_amyloid_trigger_genes_include_ttr(self):
        genes = IMAGING_TRIGGER_MAP["diffuse_lge_low_voltage"]["gene_panel"]
        assert "TTR" in genes

    def test_aortic_trigger_genes_include_fbn1(self):
        genes = IMAGING_TRIGGER_MAP["aortic_root_dilation"]["gene_panel"]
        assert "FBN1" in genes

    def test_fh_trigger_conditions(self):
        conditions = IMAGING_TRIGGER_MAP["premature_coronary_calcification"]["conditions"]
        assert "Familial Hypercholesterolemia" in conditions
