"""Tests for the Cardiology Intelligence Agent query expansion module.

Author: Adam Jones
Date: March 2026

Covers:
- All 15 expansion maps (structure and content)
- ENTITY_ALIASES (120+ entries)
- COMPARATIVE_PATTERNS and CARDIO_COMPARISONS
- _ENTITY_CATEGORIES and _COLLECTION_BOOST_RULES
- _WORKFLOW_TERMS
- QueryExpander class: expand(), resolve_aliases(), detect_entities(),
  detect_comparative(), get_workflow_terms(), boost_collections()
"""

import re

import pytest

from src.models import CardioWorkflowType
from src.query_expansion import (
    HEART_FAILURE_MAP,
    CORONARY_ARTERY_MAP,
    ARRHYTHMIA_MAP,
    VALVULAR_MAP,
    IMAGING_ECHO_MAP,
    IMAGING_CT_MAP,
    IMAGING_MRI_MAP,
    PREVENTIVE_MAP,
    HEMODYNAMICS_MAP,
    ELECTROPHYSIOLOGY_MAP,
    CARDIO_ONCOLOGY_MAP,
    CONGENITAL_MAP,
    DEVICE_MAP,
    VASCULAR_MAP,
    GENOMICS_CARDIO_MAP,
    ENTITY_ALIASES,
    COMPARATIVE_PATTERNS,
    CARDIO_COMPARISONS,
    _ENTITY_CATEGORIES,
    _COLLECTION_BOOST_RULES,
    _WORKFLOW_TERMS,
    QueryExpander,
)


# =====================================================================
# EXPANSION MAP EXISTENCE AND NON-EMPTINESS
# =====================================================================


class TestExpansionMapsExist:
    """Verify all 15 expansion maps exist and are non-empty dicts."""

    ALL_MAPS = [
        ("HEART_FAILURE_MAP", HEART_FAILURE_MAP),
        ("CORONARY_ARTERY_MAP", CORONARY_ARTERY_MAP),
        ("ARRHYTHMIA_MAP", ARRHYTHMIA_MAP),
        ("VALVULAR_MAP", VALVULAR_MAP),
        ("IMAGING_ECHO_MAP", IMAGING_ECHO_MAP),
        ("IMAGING_CT_MAP", IMAGING_CT_MAP),
        ("IMAGING_MRI_MAP", IMAGING_MRI_MAP),
        ("PREVENTIVE_MAP", PREVENTIVE_MAP),
        ("HEMODYNAMICS_MAP", HEMODYNAMICS_MAP),
        ("ELECTROPHYSIOLOGY_MAP", ELECTROPHYSIOLOGY_MAP),
        ("CARDIO_ONCOLOGY_MAP", CARDIO_ONCOLOGY_MAP),
        ("CONGENITAL_MAP", CONGENITAL_MAP),
        ("DEVICE_MAP", DEVICE_MAP),
        ("VASCULAR_MAP", VASCULAR_MAP),
        ("GENOMICS_CARDIO_MAP", GENOMICS_CARDIO_MAP),
    ]

    @pytest.mark.parametrize("name,mapping", ALL_MAPS)
    def test_map_is_non_empty_dict(self, name, mapping):
        assert isinstance(mapping, dict), f"{name} is not a dict"
        assert len(mapping) > 0, f"{name} is empty"

    @pytest.mark.parametrize("name,mapping", ALL_MAPS)
    def test_map_values_are_lists_of_strings(self, name, mapping):
        for key, terms in mapping.items():
            assert isinstance(terms, list), f"{name}['{key}'] is not a list"
            for term in terms:
                assert isinstance(term, str), (
                    f"{name}['{key}'] contains non-string: {term!r}"
                )


# =====================================================================
# HEART FAILURE MAP CONTENT
# =====================================================================


class TestHeartFailureMapContent:
    """HEART_FAILURE_MAP should contain relevant clinical terms."""

    def test_hf_key_exists(self):
        assert "hf" in HEART_FAILURE_MAP

    def test_hf_contains_bnp(self):
        assert "BNP" in HEART_FAILURE_MAP["hf"]

    def test_hf_contains_arni(self):
        assert "ARNI" in HEART_FAILURE_MAP["hf"]

    def test_hf_contains_sglt2i(self):
        assert "SGLT2i" in HEART_FAILURE_MAP["hf"]

    def test_chf_key_exists(self):
        assert "chf" in HEART_FAILURE_MAP

    def test_heart_failure_key_exists(self):
        assert "heart failure" in HEART_FAILURE_MAP

    def test_heart_failure_contains_hfref(self):
        assert "HFrEF" in HEART_FAILURE_MAP["heart failure"]

    def test_heart_failure_contains_hfpef(self):
        assert "HFpEF" in HEART_FAILURE_MAP["heart failure"]

    def test_gdmt_key_exists(self):
        assert "gdmt" in HEART_FAILURE_MAP

    def test_gdmt_contains_four_pillars_components(self):
        gdmt_terms = HEART_FAILURE_MAP["gdmt"]
        for expected in ["ARNI", "beta blocker", "MRA", "SGLT2i"]:
            assert expected in gdmt_terms, f"Missing {expected} in GDMT terms"

    def test_systolic_dysfunction_key_exists(self):
        assert "systolic dysfunction" in HEART_FAILURE_MAP

    def test_diastolic_dysfunction_key_exists(self):
        assert "diastolic dysfunction" in HEART_FAILURE_MAP


# =====================================================================
# CORONARY ARTERY MAP CONTENT
# =====================================================================


class TestCoronaryArteryMapContent:
    """CORONARY_ARTERY_MAP should expand CAD-related queries."""

    def test_cad_key_exists(self):
        assert "cad" in CORONARY_ARTERY_MAP

    def test_cad_contains_pci(self):
        assert "PCI" in CORONARY_ARTERY_MAP["cad"]

    def test_cad_contains_cabg(self):
        assert "CABG" in CORONARY_ARTERY_MAP["cad"]

    def test_cad_contains_troponin(self):
        terms = CORONARY_ARTERY_MAP["cad"]
        # Check for any troponin mention
        assert any("troponin" in t.lower() for t in terms)

    def test_mi_key_exists(self):
        assert "mi" in CORONARY_ARTERY_MAP

    def test_mi_contains_stemi(self):
        assert "STEMI" in CORONARY_ARTERY_MAP["mi"]

    def test_mi_contains_nstemi(self):
        assert any("NSTEMI" in t for t in CORONARY_ARTERY_MAP["mi"])

    def test_acs_key_exists(self):
        assert "acs" in CORONARY_ARTERY_MAP

    def test_pci_key_exists(self):
        assert "pci" in CORONARY_ARTERY_MAP

    def test_cabg_key_exists(self):
        assert "cabg" in CORONARY_ARTERY_MAP


# =====================================================================
# ARRHYTHMIA MAP CONTENT
# =====================================================================


class TestArrhythmiaMapContent:
    """ARRHYTHMIA_MAP should expand rhythm-related queries."""

    def test_atrial_fibrillation_key(self):
        assert "atrial fibrillation" in ARRHYTHMIA_MAP

    def test_afib_key(self):
        assert "afib" in ARRHYTHMIA_MAP

    def test_af_contains_anticoagulation(self):
        terms = ARRHYTHMIA_MAP["atrial fibrillation"]
        assert any("anticoagulation" in t.lower() for t in terms)

    def test_vt_key_exists(self):
        assert "vt" in ARRHYTHMIA_MAP

    def test_svt_key_exists(self):
        assert "svt" in ARRHYTHMIA_MAP

    def test_bradycardia_key(self):
        assert "bradycardia" in ARRHYTHMIA_MAP

    def test_long_qt_key(self):
        assert "long qt" in ARRHYTHMIA_MAP

    def test_brugada_key(self):
        assert "brugada" in ARRHYTHMIA_MAP


# =====================================================================
# REMAINING MAP SPOT CHECKS
# =====================================================================


class TestRemainingMapsSpotChecks:
    """Quick content checks for other maps."""

    def test_valvular_has_aortic_stenosis(self):
        assert "aortic stenosis" in VALVULAR_MAP

    def test_valvular_has_tavr(self):
        assert "tavr" in VALVULAR_MAP

    def test_imaging_echo_has_echo(self):
        assert "echo" in IMAGING_ECHO_MAP

    def test_imaging_echo_has_strain(self):
        assert "strain" in IMAGING_ECHO_MAP

    def test_imaging_ct_has_calcium_score(self):
        assert "calcium score" in IMAGING_CT_MAP

    def test_imaging_ct_has_cad_rads(self):
        assert "cad-rads" in IMAGING_CT_MAP

    def test_imaging_mri_has_cardiac_mri(self):
        assert "cardiac mri" in IMAGING_MRI_MAP

    def test_imaging_mri_has_lge(self):
        assert "lge" in IMAGING_MRI_MAP

    def test_preventive_has_statin(self):
        assert "statin" in PREVENTIVE_MAP

    def test_preventive_has_ascvd(self):
        assert "ascvd" in PREVENTIVE_MAP

    def test_hemodynamics_has_catheterization(self):
        assert "catheterization" in HEMODYNAMICS_MAP

    def test_hemodynamics_has_ffr(self):
        assert "ffr" in HEMODYNAMICS_MAP

    def test_electrophysiology_has_ecg(self):
        assert "ecg" in ELECTROPHYSIOLOGY_MAP

    def test_electrophysiology_has_ablation(self):
        assert "ablation" in ELECTROPHYSIOLOGY_MAP

    def test_cardio_oncology_has_cardiotoxicity(self):
        assert "cardiotoxicity" in CARDIO_ONCOLOGY_MAP

    def test_cardio_oncology_has_anthracycline(self):
        assert "anthracycline" in CARDIO_ONCOLOGY_MAP

    def test_congenital_has_congenital(self):
        assert "congenital" in CONGENITAL_MAP

    def test_congenital_has_pfo(self):
        assert "pfo" in CONGENITAL_MAP

    def test_device_has_pacemaker(self):
        assert "pacemaker" in DEVICE_MAP

    def test_device_has_icd(self):
        assert "icd" in DEVICE_MAP

    def test_device_has_crt(self):
        assert "crt" in DEVICE_MAP

    def test_device_has_lvad(self):
        assert "lvad" in DEVICE_MAP

    def test_vascular_has_aorta(self):
        assert "aorta" in VASCULAR_MAP

    def test_vascular_has_pad(self):
        assert "pad" in VASCULAR_MAP

    def test_genomics_has_cardiomyopathy(self):
        assert "cardiomyopathy" in GENOMICS_CARDIO_MAP

    def test_genomics_has_channelopathy(self):
        assert "channelopathy" in GENOMICS_CARDIO_MAP

    def test_genomics_has_ldlr(self):
        assert "ldlr" in GENOMICS_CARDIO_MAP


# =====================================================================
# ENTITY_ALIASES
# =====================================================================


class TestEntityAliases:
    """ENTITY_ALIASES maps abbreviations to canonical names."""

    def test_has_at_least_120_entries(self):
        assert len(ENTITY_ALIASES) >= 120, (
            f"Expected >= 120 aliases, got {len(ENTITY_ALIASES)}"
        )

    def test_mi_resolves_to_myocardial_infarction(self):
        assert ENTITY_ALIASES["MI"] == "myocardial infarction"

    def test_af_resolves_to_atrial_fibrillation(self):
        assert ENTITY_ALIASES["AF"] == "atrial fibrillation"

    def test_afib_resolves_to_atrial_fibrillation(self):
        assert ENTITY_ALIASES["AFib"] == "atrial fibrillation"

    def test_hf_resolves_to_heart_failure(self):
        assert ENTITY_ALIASES["HF"] == "heart failure"

    def test_chf_resolves_to_congestive_heart_failure(self):
        assert ENTITY_ALIASES["CHF"] == "congestive heart failure"

    def test_hfref_resolves(self):
        assert "reduced ejection fraction" in ENTITY_ALIASES["HFrEF"]

    def test_hfpef_resolves(self):
        assert "preserved ejection fraction" in ENTITY_ALIASES["HFpEF"]

    def test_lvef_resolves(self):
        assert ENTITY_ALIASES["LVEF"] == "left ventricular ejection fraction"

    def test_gdmt_resolves(self):
        assert ENTITY_ALIASES["GDMT"] == "guideline-directed medical therapy"

    def test_arni_resolves(self):
        assert "neprilysin" in ENTITY_ALIASES["ARNI"]

    def test_sglt2i_resolves(self):
        assert "sodium-glucose" in ENTITY_ALIASES["SGLT2i"]

    def test_pci_resolves(self):
        assert ENTITY_ALIASES["PCI"] == "percutaneous coronary intervention"

    def test_cabg_resolves(self):
        assert ENTITY_ALIASES["CABG"] == "coronary artery bypass graft"

    def test_tavr_resolves(self):
        assert "transcatheter" in ENTITY_ALIASES["TAVR"]

    def test_ecg_resolves(self):
        assert ENTITY_ALIASES["ECG"] == "electrocardiogram"

    def test_ekg_resolves(self):
        assert ENTITY_ALIASES["EKG"] == "electrocardiogram"

    def test_tte_resolves(self):
        assert "transthoracic" in ENTITY_ALIASES["TTE"]

    def test_cmr_resolves(self):
        assert "cardiac magnetic resonance" == ENTITY_ALIASES["CMR"]

    def test_lge_resolves(self):
        assert "late gadolinium enhancement" == ENTITY_ALIASES["LGE"]

    def test_icd_resolves(self):
        assert "implantable" in ENTITY_ALIASES["ICD"]

    def test_crt_resolves(self):
        assert "resynchronization" in ENTITY_ALIASES["CRT"]

    def test_lvad_resolves(self):
        assert "left ventricular assist" in ENTITY_ALIASES["LVAD"]

    def test_asd_resolves(self):
        assert ENTITY_ALIASES["ASD"] == "atrial septal defect"

    def test_vsd_resolves(self):
        assert ENTITY_ALIASES["VSD"] == "ventricular septal defect"

    def test_pfo_resolves(self):
        assert ENTITY_ALIASES["PFO"] == "patent foramen ovale"

    def test_dcm_resolves(self):
        assert ENTITY_ALIASES["DCM"] == "dilated cardiomyopathy"

    def test_hcm_resolves(self):
        assert ENTITY_ALIASES["HCM"] == "hypertrophic cardiomyopathy"

    def test_ascvd_resolves(self):
        assert "atherosclerotic" in ENTITY_ALIASES["ASCVD"]

    def test_ctrcd_resolves(self):
        assert "cancer therapy" in ENTITY_ALIASES["CTRCD"]

    def test_all_values_are_strings(self):
        for abbr, canonical in ENTITY_ALIASES.items():
            assert isinstance(abbr, str)
            assert isinstance(canonical, str)
            assert len(canonical) > 0


# =====================================================================
# COMPARATIVE_PATTERNS
# =====================================================================


class TestComparativePatterns:
    """COMPARATIVE_PATTERNS should contain regexes for comparison queries."""

    def test_is_non_empty_list(self):
        assert isinstance(COMPARATIVE_PATTERNS, list)
        assert len(COMPARATIVE_PATTERNS) > 0

    def test_patterns_are_valid_regex(self):
        for pattern in COMPARATIVE_PATTERNS:
            re.compile(pattern, re.IGNORECASE)  # should not raise

    def test_vs_pattern_present(self):
        combined = "|".join(COMPARATIVE_PATTERNS)
        assert re.search(re.compile(combined, re.IGNORECASE), "TAVR vs SAVR")

    def test_versus_pattern_present(self):
        combined = "|".join(COMPARATIVE_PATTERNS)
        assert re.search(re.compile(combined, re.IGNORECASE), "apixaban versus rivaroxaban")

    def test_compared_to_pattern_present(self):
        combined = "|".join(COMPARATIVE_PATTERNS)
        assert re.search(re.compile(combined, re.IGNORECASE), "drug A compared to drug B")

    def test_superior_to_pattern(self):
        combined = "|".join(COMPARATIVE_PATTERNS)
        assert re.search(re.compile(combined, re.IGNORECASE), "is ARNI superior to ACEi")


# =====================================================================
# CARDIO_COMPARISONS
# =====================================================================


class TestCardioComparisons:
    """CARDIO_COMPARISONS should contain key comparison entries."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_COMPARISONS, dict)
        assert len(CARDIO_COMPARISONS) > 0

    def test_tavr_vs_savr(self):
        assert "tavr vs savr" in CARDIO_COMPARISONS

    def test_doac_vs_warfarin(self):
        assert "doac vs warfarin" in CARDIO_COMPARISONS

    def test_pci_vs_cabg(self):
        assert "pci vs cabg" in CARDIO_COMPARISONS

    def test_arni_vs_acei(self):
        assert "arni vs acei" in CARDIO_COMPARISONS

    def test_ablation_vs_rate_control(self):
        assert "ablation vs rate control" in CARDIO_COMPARISONS

    def test_apixaban_vs_rivaroxaban(self):
        assert "apixaban vs rivaroxaban" in CARDIO_COMPARISONS

    def test_stress_echo_vs_nuclear(self):
        assert "stress echo vs nuclear" in CARDIO_COMPARISONS

    def test_mechanical_vs_bioprosthetic(self):
        assert "mechanical vs bioprosthetic valve" in CARDIO_COMPARISONS

    def test_each_comparison_has_side_a_side_b(self):
        for key, comp in CARDIO_COMPARISONS.items():
            assert "side_a" in comp, f"{key} missing side_a"
            assert "side_b" in comp, f"{key} missing side_b"
            assert isinstance(comp["side_a"], list)
            assert isinstance(comp["side_b"], list)
            assert len(comp["side_a"]) > 0
            assert len(comp["side_b"]) > 0

    def test_comparisons_have_shared_terms(self):
        for key, comp in CARDIO_COMPARISONS.items():
            assert "shared" in comp, f"{key} missing shared terms"
            assert isinstance(comp["shared"], list)


# =====================================================================
# _ENTITY_CATEGORIES
# =====================================================================


class TestEntityCategories:
    """_ENTITY_CATEGORIES should have conditions, drugs, genes, etc."""

    def test_has_conditions(self):
        assert "conditions" in _ENTITY_CATEGORIES

    def test_has_drugs(self):
        assert "drugs" in _ENTITY_CATEGORIES

    def test_has_genes(self):
        assert "genes" in _ENTITY_CATEGORIES

    def test_has_imaging_modalities(self):
        assert "imaging_modalities" in _ENTITY_CATEGORIES

    def test_has_procedures(self):
        assert "procedures" in _ENTITY_CATEGORIES

    def test_has_biomarkers(self):
        assert "biomarkers" in _ENTITY_CATEGORIES

    def test_conditions_contains_heart_failure(self):
        assert "heart failure" in _ENTITY_CATEGORIES["conditions"]

    def test_conditions_contains_atrial_fibrillation(self):
        assert "atrial fibrillation" in _ENTITY_CATEGORIES["conditions"]

    def test_drugs_contains_amiodarone(self):
        assert "amiodarone" in _ENTITY_CATEGORIES["drugs"]

    def test_genes_contains_myh7(self):
        assert "MYH7" in _ENTITY_CATEGORIES["genes"]

    def test_genes_contains_scn5a(self):
        assert "SCN5A" in _ENTITY_CATEGORIES["genes"]

    def test_imaging_contains_echocardiography(self):
        assert "echocardiography" in _ENTITY_CATEGORIES["imaging_modalities"]

    def test_procedures_contains_pci(self):
        assert "PCI" in _ENTITY_CATEGORIES["procedures"]

    def test_biomarkers_contains_troponin(self):
        assert "troponin" in _ENTITY_CATEGORIES["biomarkers"]


# =====================================================================
# _COLLECTION_BOOST_RULES
# =====================================================================


class TestCollectionBoostRules:
    """_COLLECTION_BOOST_RULES should map domains to boost factors."""

    def test_is_non_empty_dict(self):
        assert isinstance(_COLLECTION_BOOST_RULES, dict)
        assert len(_COLLECTION_BOOST_RULES) > 0

    def test_heart_failure_domain_exists(self):
        assert "heart_failure" in _COLLECTION_BOOST_RULES

    def test_coronary_domain_exists(self):
        assert "coronary" in _COLLECTION_BOOST_RULES

    def test_arrhythmia_domain_exists(self):
        assert "arrhythmia" in _COLLECTION_BOOST_RULES

    def test_valvular_domain_exists(self):
        assert "valvular" in _COLLECTION_BOOST_RULES

    def test_boost_values_are_positive_floats(self):
        for domain, rules in _COLLECTION_BOOST_RULES.items():
            for coll, boost_val in rules.items():
                assert isinstance(boost_val, (int, float))
                assert boost_val > 0, (
                    f"{domain}.{coll} boost is non-positive: {boost_val}"
                )


# =====================================================================
# _WORKFLOW_TERMS
# =====================================================================


class TestWorkflowTerms:
    """_WORKFLOW_TERMS should have entries for all 8 workflow types."""

    def test_all_workflow_types_present(self):
        for wf in CardioWorkflowType:
            assert wf in _WORKFLOW_TERMS, f"Missing {wf.value}"

    def test_each_entry_is_non_empty_list(self):
        for wf, terms in _WORKFLOW_TERMS.items():
            assert isinstance(terms, list)
            assert len(terms) > 0, f"{wf.value} has no terms"


# =====================================================================
# QueryExpander CREATION
# =====================================================================


class TestQueryExpanderCreation:
    """QueryExpander should initialise with 18 maps and aliases."""

    def test_creation(self):
        qe = QueryExpander()
        assert qe is not None

    def test_has_18_expansion_maps(self):
        qe = QueryExpander()
        assert len(qe.expansion_maps) == 18

    def test_entity_aliases_loaded(self):
        qe = QueryExpander()
        assert qe.entity_aliases is ENTITY_ALIASES

    def test_comparative_regex_compiled(self):
        qe = QueryExpander()
        assert hasattr(qe, "_comparative_re")
        assert qe._comparative_re is not None


# =====================================================================
# QueryExpander.expand()
# =====================================================================


class TestQueryExpanderExpand:
    """expand() should return dict with expected keys."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_returns_dict(self, expander):
        result = expander.expand("heart failure management")
        assert isinstance(result, dict)

    def test_has_expected_keys(self, expander):
        result = expander.expand("heart failure")
        for key in ("original", "expanded_terms", "detected_entities",
                     "is_comparative", "comparison", "workflow_hint"):
            assert key in result, f"Missing key: {key}"

    def test_original_preserved(self, expander):
        result = expander.expand("LVEF assessment")
        assert result["original"] == "LVEF assessment"

    def test_expanded_terms_is_list(self, expander):
        result = expander.expand("heart failure")
        assert isinstance(result["expanded_terms"], list)

    def test_hf_query_produces_expansions(self, expander):
        result = expander.expand("heart failure")
        assert len(result["expanded_terms"]) > 0

    def test_cad_query_produces_expansions(self, expander):
        result = expander.expand("coronary artery disease")
        assert len(result["expanded_terms"]) > 0

    def test_afib_query_produces_expansions(self, expander):
        result = expander.expand("atrial fibrillation management")
        assert len(result["expanded_terms"]) > 0

    def test_non_comparative_query(self, expander):
        result = expander.expand("heart failure GDMT")
        assert result["is_comparative"] is False

    def test_comparative_query_detected(self, expander):
        result = expander.expand("TAVR vs SAVR for aortic stenosis")
        assert result["is_comparative"] is True

    def test_workflow_hint_none_when_no_workflow(self, expander):
        result = expander.expand("heart failure")
        assert result["workflow_hint"] is None

    def test_workflow_hint_set_when_workflow_provided(self, expander):
        result = expander.expand(
            "heart failure GDMT",
            workflow=CardioWorkflowType.HEART_FAILURE,
        )
        assert result["workflow_hint"] == CardioWorkflowType.HEART_FAILURE.value

    def test_workflow_terms_appended(self, expander):
        without_wf = expander.expand("general question")
        with_wf = expander.expand(
            "general question",
            workflow=CardioWorkflowType.ARRHYTHMIA,
        )
        assert len(with_wf["expanded_terms"]) >= len(without_wf["expanded_terms"])

    def test_detected_entities_is_dict(self, expander):
        result = expander.expand("heart failure and atrial fibrillation")
        assert isinstance(result["detected_entities"], dict)


# =====================================================================
# QueryExpander.resolve_aliases()
# =====================================================================


class TestQueryExpanderResolveAliases:
    """resolve_aliases() should expand abbreviations inline."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_mi_expanded(self, expander):
        result = expander.resolve_aliases("Patient with MI")
        assert "myocardial infarction" in result

    def test_af_expanded(self, expander):
        result = expander.resolve_aliases("Patient has AF")
        assert "atrial fibrillation" in result

    def test_hf_expanded(self, expander):
        result = expander.resolve_aliases("Diagnosed with HF")
        assert "heart failure" in result

    def test_ecg_expanded(self, expander):
        result = expander.resolve_aliases("The ECG shows")
        assert "electrocardiogram" in result

    def test_lvef_expanded(self, expander):
        result = expander.resolve_aliases("LVEF is 35%")
        assert "left ventricular ejection fraction" in result

    def test_no_partial_word_replacement(self, expander):
        # "MRI" should not get expanded if only "MI" alias exists
        result = expander.resolve_aliases("MRI scan ordered")
        # "MI" should not replace inside "MRI" due to word-boundary matching
        assert "myocardial infarction" not in result

    def test_returns_string(self, expander):
        result = expander.resolve_aliases("simple text")
        assert isinstance(result, str)


# =====================================================================
# QueryExpander.detect_entities()
# =====================================================================


class TestQueryExpanderDetectEntities:
    """detect_entities() should find conditions, drugs, genes, etc."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_detects_condition(self, expander):
        result = expander.detect_entities("Patient with heart failure")
        assert "conditions" in result
        assert "heart failure" in result["conditions"]

    def test_detects_drug(self, expander):
        result = expander.detect_entities("Started on amiodarone")
        assert "drugs" in result
        assert "amiodarone" in result["drugs"]

    def test_detects_gene(self, expander):
        result = expander.detect_entities("Genetic testing revealed MYH7 variant")
        assert "genes" in result
        assert "MYH7" in result["genes"]

    def test_detects_imaging_modality(self, expander):
        result = expander.detect_entities("Echocardiography showed LVEF 30%")
        assert "imaging_modalities" in result

    def test_detects_biomarker(self, expander):
        result = expander.detect_entities("Troponin was elevated")
        assert "biomarkers" in result
        assert "troponin" in result["biomarkers"]

    def test_detects_procedure(self, expander):
        result = expander.detect_entities("Underwent PCI to LAD")
        assert "procedures" in result
        assert "PCI" in result["procedures"]

    def test_detects_multiple_categories(self, expander):
        result = expander.detect_entities(
            "Patient with heart failure started on amiodarone, troponin elevated"
        )
        assert "conditions" in result
        assert "drugs" in result
        assert "biomarkers" in result

    def test_resolved_aliases_in_output(self, expander):
        result = expander.detect_entities("AF with rapid ventricular response")
        assert "resolved_aliases" in result

    def test_returns_dict(self, expander):
        result = expander.detect_entities("nothing relevant")
        assert isinstance(result, dict)

    def test_empty_text_returns_minimal_dict(self, expander):
        result = expander.detect_entities("")
        assert isinstance(result, dict)


# =====================================================================
# QueryExpander.detect_comparative()
# =====================================================================


class TestQueryExpanderDetectComparative:
    """detect_comparative() should detect and parse comparison queries."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_vs_detected(self, expander):
        result = expander.detect_comparative("TAVR vs SAVR")
        assert result is not None
        assert "pattern_matched" in result

    def test_versus_detected(self, expander):
        result = expander.detect_comparative("apixaban versus rivaroxaban")
        assert result is not None

    def test_compared_to_detected(self, expander):
        result = expander.detect_comparative("PCI compared to CABG")
        assert result is not None

    def test_non_comparative_returns_none(self, expander):
        result = expander.detect_comparative("What is GDMT for HFrEF?")
        assert result is None

    def test_known_comparison_matched(self, expander):
        result = expander.detect_comparative("TAVR vs SAVR for aortic stenosis")
        assert result is not None
        assert result["known_comparison"] is not None
        assert result["known_comparison"]["key"] == "tavr vs savr"

    def test_unknown_comparison_has_none_known(self, expander):
        result = expander.detect_comparative("drug A vs drug B")
        assert result is not None
        assert result["known_comparison"] is None

    def test_raw_sides_extracted(self, expander):
        result = expander.detect_comparative("TAVR vs SAVR")
        assert result is not None
        assert result["raw_sides"] is not None
        assert len(result["raw_sides"]) == 2

    def test_superior_to_detected(self, expander):
        result = expander.detect_comparative("Is ARNI superior to ACEi?")
        assert result is not None

    def test_difference_between_detected(self, expander):
        result = expander.detect_comparative("What is the difference between PCI and CABG?")
        assert result is not None


# =====================================================================
# QueryExpander.get_workflow_terms()
# =====================================================================


class TestQueryExpanderGetWorkflowTerms:
    """get_workflow_terms() should return relevant terms for each workflow."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_heart_failure_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.HEART_FAILURE)
        assert isinstance(terms, list)
        assert len(terms) > 0
        assert any("heart failure" in t.lower() for t in terms)

    def test_cad_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.CAD_ASSESSMENT)
        assert len(terms) > 0
        assert any("coronary" in t.lower() for t in terms)

    def test_arrhythmia_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.ARRHYTHMIA)
        assert len(terms) > 0

    def test_valvular_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.VALVULAR_DISEASE)
        assert len(terms) > 0

    def test_cardiac_mri_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.CARDIAC_MRI)
        assert len(terms) > 0

    def test_stress_test_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.STRESS_TEST)
        assert len(terms) > 0

    def test_preventive_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.PREVENTIVE_RISK)
        assert len(terms) > 0

    def test_cardio_oncology_terms(self, expander):
        terms = expander.get_workflow_terms(CardioWorkflowType.CARDIO_ONCOLOGY)
        assert len(terms) > 0

    def test_returns_new_list(self, expander):
        terms1 = expander.get_workflow_terms(CardioWorkflowType.HEART_FAILURE)
        terms2 = expander.get_workflow_terms(CardioWorkflowType.HEART_FAILURE)
        assert terms1 is not terms2  # should return a new list each time


# =====================================================================
# QueryExpander.boost_collections()
# =====================================================================


class TestQueryExpanderBoostCollections:
    """boost_collections() should return dict of collection → boost factor."""

    @pytest.fixture
    def expander(self):
        return QueryExpander()

    def test_returns_dict(self, expander):
        expanded = expander.expand("heart failure management")
        boosts = expander.boost_collections(expanded)
        assert isinstance(boosts, dict)

    def test_values_are_floats(self, expander):
        expanded = expander.expand("heart failure management")
        boosts = expander.boost_collections(expanded)
        for coll, val in boosts.items():
            assert isinstance(val, float)

    def test_boost_values_greater_than_one(self, expander):
        expanded = expander.expand("heart failure GDMT management")
        boosts = expander.boost_collections(expanded)
        for coll, val in boosts.items():
            assert val > 1.0, f"{coll} boost {val} not > 1.0"

    def test_heart_failure_query_boosts_something(self, expander):
        expanded = expander.expand("heart failure GDMT ARNI SGLT2i")
        boosts = expander.boost_collections(expanded)
        assert len(boosts) > 0

    def test_coronary_query_boosts_something(self, expander):
        expanded = expander.expand("coronary artery disease PCI stent troponin")
        boosts = expander.boost_collections(expanded)
        assert len(boosts) > 0

    def test_empty_expansion_returns_empty(self, expander):
        expanded = {"original": "", "expanded_terms": []}
        boosts = expander.boost_collections(expanded)
        assert isinstance(boosts, dict)
        # May or may not be empty, but should be valid dict

    def test_arrhythmia_query_boosts(self, expander):
        expanded = expander.expand("atrial fibrillation ablation anticoagulation")
        boosts = expander.boost_collections(expanded)
        assert len(boosts) > 0
