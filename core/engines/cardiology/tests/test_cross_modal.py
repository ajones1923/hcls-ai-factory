"""Tests for cross_modal.py -- Cardiology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Tests the GENOMIC_TRIGGER_MAP structure, CrossModalEngine trigger evaluation,
ECG/clinical/imaging trigger pathways, trigger prioritisation, and the
PHENOTYPE_GENE_MAP reverse index.
"""

import pytest

from src.models import CrossModalTrigger, SeverityLevel
from src.cross_modal import (
    GENOMIC_TRIGGER_MAP,
    GENOMIC_COLLECTION,
    TOP_K_PER_QUERY,
    SCORE_THRESHOLD,
    EMBEDDING_DIM,
    FAMILY_SCREENING_PROTOCOLS,
    GENE_PANEL_COSTS,
    PHENOTYPE_GENE_MAP,
    CrossModalEngine,
    prioritize_triggers,
    prioritize_trigger_models,
    genes_for_phenotype,
    phenotypes_for_gene,
    estimate_panel_cost,
    _severity_from_urgency,
    _build_rationale,
    _TriggerContext,
)


# =====================================================================
# GENOMIC_TRIGGER_MAP STRUCTURE
# =====================================================================


class TestGenomicTriggerMap:
    """Verify GENOMIC_TRIGGER_MAP has >= 12 entries with correct structure."""

    def test_trigger_map_has_at_least_12_entries(self):
        assert len(GENOMIC_TRIGGER_MAP) >= 12

    def test_each_trigger_has_gene_panel(self):
        for key, entry in GENOMIC_TRIGGER_MAP.items():
            assert "gene_panel" in entry, f"Missing gene_panel in {key}"
            assert isinstance(entry["gene_panel"], list)
            assert len(entry["gene_panel"]) > 0

    def test_each_trigger_has_conditions(self):
        for key, entry in GENOMIC_TRIGGER_MAP.items():
            assert "conditions" in entry, f"Missing conditions in {key}"
            assert isinstance(entry["conditions"], list)

    def test_each_trigger_has_criteria(self):
        for key, entry in GENOMIC_TRIGGER_MAP.items():
            assert "criteria" in entry, f"Missing criteria in {key}"
            assert isinstance(entry["criteria"], str)
            assert len(entry["criteria"]) > 0

    def test_each_trigger_has_urgency(self):
        valid_urgencies = {"critical", "high", "moderate", "low"}
        for key, entry in GENOMIC_TRIGGER_MAP.items():
            assert "urgency" in entry, f"Missing urgency in {key}"
            assert entry["urgency"] in valid_urgencies, (
                f"Invalid urgency '{entry['urgency']}' in {key}"
            )

    # -- Specific trigger content --

    def test_unexplained_lvh_includes_myh7(self):
        panel = GENOMIC_TRIGGER_MAP["unexplained_lvh"]["gene_panel"]
        assert "MYH7" in panel

    def test_unexplained_lvh_includes_mybpc3(self):
        panel = GENOMIC_TRIGGER_MAP["unexplained_lvh"]["gene_panel"]
        assert "MYBPC3" in panel

    def test_unexplained_lvh_urgency_high(self):
        assert GENOMIC_TRIGGER_MAP["unexplained_lvh"]["urgency"] == "high"

    def test_long_qt_includes_kcnq1(self):
        panel = GENOMIC_TRIGGER_MAP["long_qt"]["gene_panel"]
        assert "KCNQ1" in panel

    def test_long_qt_includes_kcnh2(self):
        panel = GENOMIC_TRIGGER_MAP["long_qt"]["gene_panel"]
        assert "KCNH2" in panel

    def test_long_qt_includes_scn5a(self):
        panel = GENOMIC_TRIGGER_MAP["long_qt"]["gene_panel"]
        assert "SCN5A" in panel

    def test_long_qt_urgency_critical(self):
        assert GENOMIC_TRIGGER_MAP["long_qt"]["urgency"] == "critical"

    def test_premature_cad_includes_ldlr(self):
        panel = GENOMIC_TRIGGER_MAP["premature_cad"]["gene_panel"]
        assert "LDLR" in panel

    def test_premature_cad_includes_pcsk9(self):
        panel = GENOMIC_TRIGGER_MAP["premature_cad"]["gene_panel"]
        assert "PCSK9" in panel

    def test_premature_cad_includes_apob(self):
        panel = GENOMIC_TRIGGER_MAP["premature_cad"]["gene_panel"]
        assert "APOB" in panel

    def test_aortic_dilation_includes_fbn1(self):
        panel = GENOMIC_TRIGGER_MAP["aortic_dilation"]["gene_panel"]
        assert "FBN1" in panel

    def test_aortic_dilation_includes_tgfbr1(self):
        panel = GENOMIC_TRIGGER_MAP["aortic_dilation"]["gene_panel"]
        assert "TGFBR1" in panel

    def test_brugada_pattern_urgency_critical(self):
        assert GENOMIC_TRIGGER_MAP["brugada_pattern"]["urgency"] == "critical"

    def test_cpvt_includes_ryr2(self):
        panel = GENOMIC_TRIGGER_MAP["cpvt_suspected"]["gene_panel"]
        assert "RYR2" in panel

    def test_unexplained_dcm_includes_ttn(self):
        panel = GENOMIC_TRIGGER_MAP["unexplained_dcm"]["gene_panel"]
        assert "TTN" in panel

    def test_cardiac_amyloid_includes_ttr(self):
        panel = GENOMIC_TRIGGER_MAP["cardiac_amyloid"]["gene_panel"]
        assert "TTR" in panel

    def test_fabry_disease_includes_gla(self):
        panel = GENOMIC_TRIGGER_MAP["fabry_disease"]["gene_panel"]
        assert "GLA" in panel

    def test_scd_family_history_key_exists(self):
        assert "scd_family_history" in GENOMIC_TRIGGER_MAP

    def test_arrhythmogenic_cm_key_exists(self):
        assert "arrhythmogenic_cm" in GENOMIC_TRIGGER_MAP

    def test_lvnc_key_exists(self):
        assert "lvnc" in GENOMIC_TRIGGER_MAP


# =====================================================================
# CONSTANTS
# =====================================================================


class TestConstants:
    """Test module-level constants."""

    def test_genomic_collection_name(self):
        assert GENOMIC_COLLECTION == "genomic_evidence"

    def test_top_k_per_query(self):
        assert TOP_K_PER_QUERY > 0

    def test_score_threshold(self):
        assert 0 < SCORE_THRESHOLD < 1

    def test_embedding_dim(self):
        assert EMBEDDING_DIM == 384


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================


class TestHelperFunctions:
    """Test module-level helper functions."""

    def test_severity_from_urgency_critical(self):
        assert _severity_from_urgency("critical") == SeverityLevel.CRITICAL

    def test_severity_from_urgency_high(self):
        assert _severity_from_urgency("high") == SeverityLevel.HIGH

    def test_severity_from_urgency_moderate(self):
        assert _severity_from_urgency("moderate") == SeverityLevel.MODERATE

    def test_severity_from_urgency_low(self):
        assert _severity_from_urgency("low") == SeverityLevel.LOW

    def test_severity_from_urgency_default(self):
        assert _severity_from_urgency("unknown") == SeverityLevel.MODERATE

    def test_build_rationale_includes_criteria(self):
        entry = {"criteria": "QTc > 500ms", "guideline": "ACC/AHA", "urgency": "critical"}
        result = _build_rationale(entry, "ecg", "QTc 520", 1.0)
        assert "QTc > 500ms" in result

    def test_build_rationale_includes_urgency(self):
        entry = {"criteria": "test", "urgency": "high"}
        result = _build_rationale(entry, "imaging", "", 1.0)
        assert "high" in result.lower()

    def test_build_rationale_includes_low_confidence(self):
        entry = {"criteria": "test", "urgency": "moderate"}
        result = _build_rationale(entry, "imaging", "", 0.6)
        assert "60%" in result


# =====================================================================
# GENES / PHENOTYPE LOOKUPS
# =====================================================================


class TestGenePhenotypeLookups:
    """Test genes_for_phenotype and phenotypes_for_gene functions."""

    def test_genes_for_unexplained_lvh(self):
        genes = genes_for_phenotype("unexplained_lvh")
        assert "MYH7" in genes

    def test_genes_for_unknown(self):
        genes = genes_for_phenotype("nonexistent_key")
        assert genes == []

    def test_phenotypes_for_myh7(self):
        phenos = phenotypes_for_gene("MYH7")
        assert len(phenos) > 0
        assert "unexplained_lvh" in phenos

    def test_phenotypes_for_unknown_gene(self):
        phenos = phenotypes_for_gene("FAKEGENE")
        assert phenos == []


# =====================================================================
# PHENOTYPE_GENE_MAP
# =====================================================================


class TestPhenotypeGeneMap:
    """Test PHENOTYPE_GENE_MAP reverse mapping."""

    def test_phenotype_gene_map_is_dict(self):
        assert isinstance(PHENOTYPE_GENE_MAP, dict)

    def test_phenotype_gene_map_has_entries(self):
        assert len(PHENOTYPE_GENE_MAP) > 0

    def test_myh7_in_reverse_map(self):
        assert "MYH7" in PHENOTYPE_GENE_MAP

    def test_ttn_in_reverse_map(self):
        assert "TTN" in PHENOTYPE_GENE_MAP

    def test_scn5a_in_reverse_map(self):
        assert "SCN5A" in PHENOTYPE_GENE_MAP

    def test_ldlr_in_reverse_map(self):
        assert "LDLR" in PHENOTYPE_GENE_MAP

    def test_fbn1_in_reverse_map(self):
        assert "FBN1" in PHENOTYPE_GENE_MAP

    def test_reverse_map_values_are_sorted_lists(self):
        for gene, phenos in PHENOTYPE_GENE_MAP.items():
            assert isinstance(phenos, list)
            assert phenos == sorted(phenos), f"Phenotypes for {gene} are not sorted"


# =====================================================================
# FAMILY SCREENING PROTOCOLS
# =====================================================================


class TestFamilyScreeningProtocols:
    """Test FAMILY_SCREENING_PROTOCOLS structure."""

    def test_protocols_exist(self):
        assert len(FAMILY_SCREENING_PROTOCOLS) > 0

    def test_hcm_protocol_exists(self):
        assert "Hypertrophic Cardiomyopathy" in FAMILY_SCREENING_PROTOCOLS

    def test_dcm_protocol_exists(self):
        assert "Dilated Cardiomyopathy" in FAMILY_SCREENING_PROTOCOLS

    def test_lqts_protocol_exists(self):
        assert "Long QT Syndrome" in FAMILY_SCREENING_PROTOCOLS

    def test_protocol_has_who(self):
        for name, proto in FAMILY_SCREENING_PROTOCOLS.items():
            assert "who" in proto, f"Missing 'who' in {name}"

    def test_protocol_has_modalities(self):
        for name, proto in FAMILY_SCREENING_PROTOCOLS.items():
            assert "modalities" in proto, f"Missing 'modalities' in {name}"
            assert isinstance(proto["modalities"], list)


# =====================================================================
# GENE PANEL COSTS
# =====================================================================


class TestGenePanelCosts:
    """Test GENE_PANEL_COSTS structure."""

    def test_panel_costs_exist(self):
        assert len(GENE_PANEL_COSTS) > 0

    def test_each_panel_has_cost(self):
        for name, panel in GENE_PANEL_COSTS.items():
            assert "estimated_cost_usd" in panel, f"Missing cost in {name}"
            assert panel["estimated_cost_usd"] > 0

    def test_each_panel_has_turnaround(self):
        for name, panel in GENE_PANEL_COSTS.items():
            assert "turnaround_days" in panel, f"Missing turnaround in {name}"

    def test_hcm_panel_exists(self):
        assert "HCM_panel" in GENE_PANEL_COSTS

    def test_fh_panel_exists(self):
        assert "FH_panel" in GENE_PANEL_COSTS


# =====================================================================
# CROSS-MODAL ENGINE
# =====================================================================


class TestCrossModalEngine:
    """Test CrossModalEngine creation and trigger evaluation."""

    @pytest.fixture
    def engine(self):
        return CrossModalEngine()

    def test_engine_creation(self):
        engine = CrossModalEngine()
        assert engine is not None
        assert engine.cascade_enabled is True

    def test_engine_creation_no_cascade(self):
        engine = CrossModalEngine(cascade_enabled=False)
        assert engine.cascade_enabled is False

    def test_engine_has_trigger_map(self, engine):
        assert engine.trigger_map is GENOMIC_TRIGGER_MAP

    # -- evaluate_triggers --

    def test_evaluate_triggers_empty_findings(self, engine):
        result = engine.evaluate_triggers({})
        assert isinstance(result, list)

    def test_evaluate_triggers_with_lvh(self, engine):
        findings = {
            "modality": "echocardiography",
            "measurements": {"lv_wall_thickness_mm": 18},
            "imaging_findings": [],
        }
        result = engine.evaluate_triggers(findings)
        assert isinstance(result, list)

    def test_evaluate_triggers_returns_crossmodal_triggers(self, engine):
        findings = {
            "modality": "echocardiography",
            "measurements": {"lv_wall_thickness_mm": 18},
            "imaging_findings": [],
        }
        result = engine.evaluate_triggers(findings)
        for trigger in result:
            assert isinstance(trigger, CrossModalTrigger)

    # -- check_ecg_triggers --

    def test_check_ecg_triggers_high_qtc(self, engine):
        result = engine.check_ecg_triggers(
            qtc=520, rhythm="sinus", morphology="normal", family_hx={},
        )
        assert isinstance(result, list)

    def test_check_ecg_triggers_brugada(self, engine):
        result = engine.check_ecg_triggers(
            qtc=400, rhythm="sinus", morphology="brugada_type1", family_hx={},
        )
        assert isinstance(result, list)

    def test_check_ecg_triggers_normal(self, engine):
        result = engine.check_ecg_triggers(
            qtc=420, rhythm="sinus", morphology="normal", family_hx={},
        )
        assert isinstance(result, list)

    # -- check_clinical_triggers --

    def test_check_clinical_triggers_high_ldl(self, engine):
        result = engine.check_clinical_triggers(
            age=50, ldl=200, family_hx={}, conditions=[],
        )
        assert isinstance(result, list)

    def test_check_clinical_triggers_young_scd_family(self, engine):
        result = engine.check_clinical_triggers(
            age=30, ldl=100,
            family_hx={"scd": True},
            conditions=[],
        )
        assert isinstance(result, list)

    def test_check_clinical_triggers_normal(self, engine):
        result = engine.check_clinical_triggers(
            age=50, ldl=100, family_hx={}, conditions=[],
        )
        assert isinstance(result, list)

    # -- check_imaging_triggers --

    def test_check_imaging_triggers_lvh(self, engine):
        result = engine.check_imaging_triggers(
            modality="echocardiography",
            measurements={"lv_wall_thickness_mm": 18},
            findings=[],
        )
        assert isinstance(result, list)

    def test_check_imaging_triggers_normal(self, engine):
        result = engine.check_imaging_triggers(
            modality="echocardiography",
            measurements={"lv_wall_thickness_mm": 10},
            findings=[],
        )
        assert isinstance(result, list)

    # -- build_genomic_query --

    def test_build_genomic_query(self, engine):
        trigger = CrossModalTrigger(
            trigger_source="imaging",
            finding="unexplained_lvh",
            gene_panel=["MYH7", "MYBPC3"],
            conditions=["Hypertrophic Cardiomyopathy"],
            rationale="Test",
        )
        query = engine.build_genomic_query(trigger)
        assert query["collection"] == GENOMIC_COLLECTION
        assert len(query["query_texts"]) > 0
        assert query["top_k"] == TOP_K_PER_QUERY

    # -- format_trigger_report --

    def test_format_trigger_report_empty(self, engine):
        report = engine.format_trigger_report([])
        assert "No cross-modal" in report

    def test_format_trigger_report_with_triggers(self, engine):
        trigger = CrossModalTrigger(
            trigger_source="ecg",
            finding="long_qt",
            gene_panel=["KCNQ1"],
            conditions=["Long QT Syndrome"],
            rationale="Urgency: critical",
        )
        report = engine.format_trigger_report([trigger])
        assert "CROSS-MODAL" in report
        assert "CRITICAL" in report


# =====================================================================
# PRIORITIZE TRIGGERS
# =====================================================================


class TestPrioritizeTriggers:
    """Test trigger prioritisation by urgency."""

    def test_prioritize_sorts_by_urgency(self):
        ctx_low = _TriggerContext("test1", SeverityLevel.LOW)
        ctx_critical = _TriggerContext("test2", SeverityLevel.CRITICAL)
        ctx_high = _TriggerContext("test3", SeverityLevel.HIGH)

        trigger1 = CrossModalTrigger(
            trigger_source="t", finding="f1", gene_panel=["G1"],
            conditions=["C1"], rationale="Urgency: low",
        )
        trigger2 = CrossModalTrigger(
            trigger_source="t", finding="f2", gene_panel=["G2"],
            conditions=["C2"], rationale="Urgency: critical",
        )
        trigger3 = CrossModalTrigger(
            trigger_source="t", finding="f3", gene_panel=["G3"],
            conditions=["C3"], rationale="Urgency: high",
        )

        pairs = [(trigger1, ctx_low), (trigger2, ctx_critical), (trigger3, ctx_high)]
        sorted_pairs = prioritize_triggers(pairs)
        assert sorted_pairs[0][1].urgency == SeverityLevel.CRITICAL
        assert sorted_pairs[1][1].urgency == SeverityLevel.HIGH
        assert sorted_pairs[2][1].urgency == SeverityLevel.LOW

    def test_prioritize_trigger_models(self):
        t1 = CrossModalTrigger(
            trigger_source="t", finding="f1", gene_panel=["G"],
            conditions=["C"], rationale="Urgency: low",
        )
        t2 = CrossModalTrigger(
            trigger_source="t", finding="f2", gene_panel=["G"],
            conditions=["C"], rationale="Urgency: critical",
        )
        sorted_triggers = prioritize_trigger_models([t1, t2])
        assert "critical" in sorted_triggers[0].rationale.lower()

    def test_estimate_panel_cost_hcm(self):
        trigger = CrossModalTrigger(
            trigger_source="imaging", finding="lvh",
            gene_panel=["MYH7", "MYBPC3", "TNNT2"],
            conditions=["HCM"], rationale="test",
        )
        cost = estimate_panel_cost(trigger)
        assert cost is not None
        assert "estimated_cost_usd" in cost

    def test_estimate_panel_cost_no_match(self):
        trigger = CrossModalTrigger(
            trigger_source="imaging", finding="x",
            gene_panel=["FAKEGENE1"],
            conditions=["Unknown"], rationale="test",
        )
        cost = estimate_panel_cost(trigger)
        # May or may not match -- test it doesn't crash
        assert cost is None or isinstance(cost, dict)
