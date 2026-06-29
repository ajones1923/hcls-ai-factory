"""Tests for the Cardiology Intelligence Agent main agent module.

Author: Adam Jones
Date: March 2026

Covers:
- CARDIO_SYSTEM_PROMPT content and clinical standards
- WORKFLOW_COLLECTION_BOOST structure (all 8 workflow types)
- CARDIO_CONDITIONS, CARDIO_BIOMARKERS, CARDIO_DRUGS, CARDIO_GENES,
  CARDIO_IMAGING_MODALITIES, CARDIO_GUIDELINES knowledge dictionaries
- SearchPlan dataclass
- CardioIntelligenceAgent: creation, run(), search_plan(),
  _detect_conditions(), _detect_drugs(), _detect_genes(), _detect_imaging(),
  _determine_workflows(), _choose_strategy(), _generate_sub_questions(),
  _evaluate_evidence_quality(), generate_report()
"""

from dataclasses import fields as dc_fields
from unittest.mock import MagicMock

import pytest

from src.models import CardioWorkflowType
from src.agent import (
    CARDIO_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    CARDIO_CONDITIONS,
    CARDIO_BIOMARKERS,
    CARDIO_DRUGS,
    CARDIO_GENES,
    CARDIO_IMAGING_MODALITIES,
    CARDIO_GUIDELINES,
    SearchPlan,
    CardioIntelligenceAgent,
)


# =====================================================================
# CARDIO_SYSTEM_PROMPT CONTENT
# =====================================================================


class TestCardioSystemPrompt:
    """CARDIO_SYSTEM_PROMPT must mention key clinical standards."""

    def test_is_non_empty_string(self):
        assert isinstance(CARDIO_SYSTEM_PROMPT, str)
        assert len(CARDIO_SYSTEM_PROMPT) > 200

    def test_mentions_acc_aha(self):
        assert "ACC/AHA" in CARDIO_SYSTEM_PROMPT

    def test_mentions_esc(self):
        assert "ESC" in CARDIO_SYSTEM_PROMPT

    def test_mentions_hrs(self):
        assert "HRS" in CARDIO_SYSTEM_PROMPT

    def test_mentions_scai(self):
        assert "SCAI" in CARDIO_SYSTEM_PROMPT

    def test_mentions_critical_findings(self):
        assert "CRITICAL" in CARDIO_SYSTEM_PROMPT

    def test_mentions_pmid_format(self):
        assert "PMID" in CARDIO_SYSTEM_PROMPT

    def test_mentions_class_of_recommendation(self):
        assert "Class" in CARDIO_SYSTEM_PROMPT

    def test_mentions_level_of_evidence(self):
        assert "Level of Evidence" in CARDIO_SYSTEM_PROMPT

    def test_mentions_stemi(self):
        assert "STEMI" in CARDIO_SYSTEM_PROMPT

    def test_mentions_cardiac_tamponade(self):
        assert "tamponade" in CARDIO_SYSTEM_PROMPT.lower()

    def test_mentions_ventricular_tachycardia(self):
        assert "ventricular tachycardia" in CARDIO_SYSTEM_PROMPT.lower()

    def test_mentions_gdmt(self):
        assert "GDMT" in CARDIO_SYSTEM_PROMPT

    def test_mentions_risk_calculators(self):
        assert "CHA2DS2-VASc" in CARDIO_SYSTEM_PROMPT

    def test_mentions_limitations_disclaimer(self):
        prompt_lower = CARDIO_SYSTEM_PROMPT.lower()
        assert "clinical decision support" in prompt_lower
        assert "physician judgment" in prompt_lower or "physician judgement" in prompt_lower


# =====================================================================
# WORKFLOW_COLLECTION_BOOST
# =====================================================================


class TestWorkflowCollectionBoost:
    """WORKFLOW_COLLECTION_BOOST should cover all 8 workflow types."""

    def test_is_dict(self):
        assert isinstance(WORKFLOW_COLLECTION_BOOST, dict)

    def test_all_workflow_types_present(self):
        for wf in CardioWorkflowType:
            assert wf in WORKFLOW_COLLECTION_BOOST, (
                f"Missing workflow: {wf.value}"
            )

    def test_each_entry_is_dict_of_floats(self):
        for wf, boosts in WORKFLOW_COLLECTION_BOOST.items():
            assert isinstance(boosts, dict), f"{wf.value} boosts not a dict"
            for coll, weight in boosts.items():
                assert isinstance(coll, str)
                assert isinstance(weight, (int, float))
                assert weight > 0

    def test_cad_assessment_boosts_imaging(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.CAD_ASSESSMENT]
        assert "cardio_imaging" in boosts
        assert boosts["cardio_imaging"] > 1.0

    def test_heart_failure_boosts_heart_failure_collection(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.HEART_FAILURE]
        assert "cardio_heart_failure" in boosts
        assert boosts["cardio_heart_failure"] > 1.0

    def test_valvular_boosts_valvular_collection(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.VALVULAR_DISEASE]
        assert "cardio_valvular" in boosts
        assert boosts["cardio_valvular"] > 1.0

    def test_arrhythmia_boosts_ep_collection(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.ARRHYTHMIA]
        assert "cardio_electrophysiology" in boosts
        assert boosts["cardio_electrophysiology"] > 1.0

    def test_cardio_oncology_boosts_oncology_collection(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.CARDIO_ONCOLOGY]
        assert "cardio_oncology" in boosts
        assert boosts["cardio_oncology"] > 1.0

    def test_preventive_boosts_prevention_collection(self):
        boosts = WORKFLOW_COLLECTION_BOOST[CardioWorkflowType.PREVENTIVE_RISK]
        assert "cardio_prevention" in boosts
        assert boosts["cardio_prevention"] > 1.0


# =====================================================================
# KNOWLEDGE DICTIONARIES
# =====================================================================


class TestCardioConditions:
    """CARDIO_CONDITIONS should have clinical conditions with aliases."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_CONDITIONS, dict)
        assert len(CARDIO_CONDITIONS) > 10

    def test_heart_failure_entry(self):
        assert "heart failure" in CARDIO_CONDITIONS
        entry = CARDIO_CONDITIONS["heart failure"]
        assert "aliases" in entry
        assert "workflows" in entry

    def test_atrial_fibrillation_entry(self):
        assert "atrial fibrillation" in CARDIO_CONDITIONS

    def test_aortic_stenosis_entry(self):
        assert "aortic stenosis" in CARDIO_CONDITIONS

    def test_stemi_is_critical(self):
        assert CARDIO_CONDITIONS["stemi"].get("critical") is True

    def test_ventricular_tachycardia_is_critical(self):
        assert CARDIO_CONDITIONS["ventricular tachycardia"].get("critical") is True

    def test_aortic_dissection_is_critical(self):
        assert CARDIO_CONDITIONS["aortic dissection"].get("critical") is True

    def test_cardiac_tamponade_is_critical(self):
        assert CARDIO_CONDITIONS["cardiac tamponade"].get("critical") is True


class TestCardioBiomarkers:
    """CARDIO_BIOMARKERS should contain key cardiac biomarkers."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_BIOMARKERS, dict)
        assert len(CARDIO_BIOMARKERS) > 5

    def test_troponin_entry(self):
        assert "troponin" in CARDIO_BIOMARKERS
        assert "full_name" in CARDIO_BIOMARKERS["troponin"]
        assert "normal_range" in CARDIO_BIOMARKERS["troponin"]

    def test_bnp_entry(self):
        assert "bnp" in CARDIO_BIOMARKERS

    def test_nt_probnp_entry(self):
        assert "nt-probnp" in CARDIO_BIOMARKERS

    def test_ldl_c_entry(self):
        assert "ldl-c" in CARDIO_BIOMARKERS


class TestCardioDrugs:
    """CARDIO_DRUGS should contain cardiac medications."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_DRUGS, dict)
        assert len(CARDIO_DRUGS) > 10

    def test_sacubitril_valsartan(self):
        assert "sacubitril/valsartan" in CARDIO_DRUGS
        entry = CARDIO_DRUGS["sacubitril/valsartan"]
        assert "aliases" in entry
        assert "entresto" in [a.lower() for a in entry["aliases"]]

    def test_amiodarone(self):
        assert "amiodarone" in CARDIO_DRUGS

    def test_apixaban(self):
        assert "apixaban" in CARDIO_DRUGS

    def test_clopidogrel_has_pgx(self):
        entry = CARDIO_DRUGS["clopidogrel"]
        assert entry.get("pgx_relevant") is True
        assert entry.get("pgx_gene") == "CYP2C19"


class TestCardioGenes:
    """CARDIO_GENES should contain cardiovascular gene entries."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_GENES, dict)
        assert len(CARDIO_GENES) > 10

    def test_myh7(self):
        assert "MYH7" in CARDIO_GENES
        assert "conditions" in CARDIO_GENES["MYH7"]

    def test_scn5a(self):
        assert "SCN5A" in CARDIO_GENES
        conditions = CARDIO_GENES["SCN5A"]["conditions"]
        assert any("Brugada" in c for c in conditions)

    def test_ttn(self):
        assert "TTN" in CARDIO_GENES

    def test_lmna(self):
        assert "LMNA" in CARDIO_GENES

    def test_cyp2c19_is_pgx_relevant(self):
        assert CARDIO_GENES["CYP2C19"].get("pgx_relevant") is True


class TestCardioImagingModalities:
    """CARDIO_IMAGING_MODALITIES should map imaging types."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_IMAGING_MODALITIES, dict)
        assert len(CARDIO_IMAGING_MODALITIES) > 3

    def test_echocardiography(self):
        assert "echocardiography" in CARDIO_IMAGING_MODALITIES
        entry = CARDIO_IMAGING_MODALITIES["echocardiography"]
        assert "aliases" in entry
        assert "workflows" in entry

    def test_cardiac_ct(self):
        assert "cardiac ct" in CARDIO_IMAGING_MODALITIES

    def test_cardiac_mri(self):
        assert "cardiac mri" in CARDIO_IMAGING_MODALITIES


class TestCardioGuidelines:
    """CARDIO_GUIDELINES should contain guideline references."""

    def test_is_non_empty_dict(self):
        assert isinstance(CARDIO_GUIDELINES, dict)
        assert len(CARDIO_GUIDELINES) > 5

    def test_hf_guideline(self):
        assert "acc_aha_hf_2022" in CARDIO_GUIDELINES
        entry = CARDIO_GUIDELINES["acc_aha_hf_2022"]
        assert "title" in entry
        assert "pmid" in entry
        assert "year" in entry

    def test_af_guideline(self):
        assert "acc_aha_af_2023" in CARDIO_GUIDELINES

    def test_cholesterol_guideline(self):
        assert "acc_aha_cholesterol_2018" in CARDIO_GUIDELINES


# =====================================================================
# SearchPlan DATACLASS
# =====================================================================


class TestSearchPlan:
    """SearchPlan should hold all detected entities and strategy."""

    def test_creation_with_question(self):
        plan = SearchPlan(question="What is GDMT?")
        assert plan.question == "What is GDMT?"

    def test_defaults(self):
        plan = SearchPlan(question="test")
        assert plan.conditions == []
        assert plan.drugs == []
        assert plan.genes == []
        assert plan.imaging_modalities == []
        assert plan.relevant_workflows == []
        assert plan.search_strategy == "broad"
        assert plan.sub_questions == []
        assert plan.identified_topics == []

    def test_has_expected_fields(self):
        field_names = {f.name for f in dc_fields(SearchPlan)}
        for expected in ("question", "conditions", "drugs", "genes",
                         "imaging_modalities", "relevant_workflows",
                         "search_strategy", "sub_questions", "identified_topics"):
            assert expected in field_names


# =====================================================================
# CardioIntelligenceAgent CREATION
# =====================================================================


class TestCardioIntelligenceAgentCreation:
    """CardioIntelligenceAgent should initialise with a RAG engine."""

    def test_creation(self):
        mock_rag = MagicMock()
        agent = CardioIntelligenceAgent(mock_rag)
        assert agent.rag is mock_rag

    def test_knowledge_loaded(self):
        mock_rag = MagicMock()
        agent = CardioIntelligenceAgent(mock_rag)
        assert "conditions" in agent.knowledge
        assert "biomarkers" in agent.knowledge
        assert "drugs" in agent.knowledge
        assert "genes" in agent.knowledge
        assert "imaging" in agent.knowledge
        assert "guidelines" in agent.knowledge


# =====================================================================
# search_plan() ENTITY DETECTION
# =====================================================================


class TestSearchPlanEntityDetection:
    """search_plan() should detect conditions, drugs, genes, imaging."""

    @pytest.fixture
    def agent(self):
        mock_rag = MagicMock()
        return CardioIntelligenceAgent(mock_rag)

    # -- Conditions --

    def test_detects_heart_failure(self, agent):
        plan = agent.search_plan("Patient with heart failure and LVEF 25%")
        assert "heart failure" in plan.conditions

    def test_detects_atrial_fibrillation(self, agent):
        plan = agent.search_plan("New onset atrial fibrillation management")
        assert "atrial fibrillation" in plan.conditions

    def test_detects_aortic_stenosis(self, agent):
        plan = agent.search_plan("Severe aortic stenosis evaluation")
        assert "aortic stenosis" in plan.conditions

    def test_detects_stemi(self, agent):
        plan = agent.search_plan("Acute STEMI management")
        assert "stemi" in plan.conditions

    def test_detects_cardiomyopathy(self, agent):
        plan = agent.search_plan("Dilated cardiomyopathy workup")
        assert "cardiomyopathy" in plan.conditions

    def test_detects_coronary_artery_disease(self, agent):
        plan = agent.search_plan("Coronary artery disease risk assessment")
        assert "coronary artery disease" in plan.conditions

    # -- Drugs --

    def test_detects_sacubitril_valsartan(self, agent):
        plan = agent.search_plan("Start sacubitril/valsartan for HFrEF")
        assert "sacubitril/valsartan" in plan.drugs

    def test_detects_amiodarone(self, agent):
        plan = agent.search_plan("Amiodarone for ventricular tachycardia")
        assert "amiodarone" in plan.drugs

    def test_detects_apixaban(self, agent):
        plan = agent.search_plan("Apixaban dosing for AF")
        assert "apixaban" in plan.drugs

    def test_detects_entresto_alias(self, agent):
        plan = agent.search_plan("Patient on Entresto 97/103 mg BID")
        assert "sacubitril/valsartan" in plan.drugs

    def test_detects_clopidogrel(self, agent):
        plan = agent.search_plan("Clopidogrel resistance post-PCI")
        assert "clopidogrel" in plan.drugs

    # -- Genes --

    def test_detects_myh7(self, agent):
        plan = agent.search_plan("MYH7 variant identified in HCM screening")
        assert "MYH7" in plan.genes

    def test_detects_scn5a(self, agent):
        plan = agent.search_plan("SCN5A pathogenic variant and Brugada")
        assert "SCN5A" in plan.genes

    def test_detects_ttn(self, agent):
        plan = agent.search_plan("TTN truncating variant in DCM")
        assert "TTN" in plan.genes

    def test_detects_cyp2c19(self, agent):
        plan = agent.search_plan("CYP2C19 poor metabolizer and clopidogrel")
        assert "CYP2C19" in plan.genes

    # -- Imaging --

    def test_detects_echocardiography(self, agent):
        plan = agent.search_plan("Echocardiography for valvular assessment")
        assert "echocardiography" in plan.imaging_modalities

    def test_detects_cardiac_mri(self, agent):
        plan = agent.search_plan("Cardiac MRI for myocarditis")
        assert "cardiac mri" in plan.imaging_modalities

    def test_detects_cardiac_ct_via_alias(self, agent):
        plan = agent.search_plan("Coronary CTA for chest pain evaluation")
        assert "cardiac ct" in plan.imaging_modalities


# =====================================================================
# search_plan() WORKFLOW DETERMINATION
# =====================================================================


class TestSearchPlanWorkflowDetermination:
    """search_plan() should determine relevant workflows."""

    @pytest.fixture
    def agent(self):
        return CardioIntelligenceAgent(MagicMock())

    def test_heart_failure_workflow(self, agent):
        plan = agent.search_plan("Heart failure with LVEF 30%")
        assert CardioWorkflowType.HEART_FAILURE in plan.relevant_workflows

    def test_cad_workflow(self, agent):
        plan = agent.search_plan("Coronary artery disease assessment")
        assert CardioWorkflowType.CAD_ASSESSMENT in plan.relevant_workflows

    def test_valvular_workflow(self, agent):
        plan = agent.search_plan("Severe aortic stenosis TAVR evaluation")
        assert CardioWorkflowType.VALVULAR_DISEASE in plan.relevant_workflows

    def test_arrhythmia_workflow(self, agent):
        plan = agent.search_plan("Atrial fibrillation ablation candidacy")
        assert CardioWorkflowType.ARRHYTHMIA in plan.relevant_workflows

    def test_cardio_oncology_workflow(self, agent):
        plan = agent.search_plan("Cardiotoxicity from doxorubicin chemotherapy")
        assert CardioWorkflowType.CARDIO_ONCOLOGY in plan.relevant_workflows

    def test_preventive_workflow(self, agent):
        plan = agent.search_plan("ASCVD risk assessment and statin therapy")
        assert CardioWorkflowType.PREVENTIVE_RISK in plan.relevant_workflows


# =====================================================================
# search_plan() COMPARATIVE QUERIES & STRATEGY
# =====================================================================


class TestSearchPlanComparativeAndStrategy:
    """search_plan() should handle comparative queries and strategy."""

    @pytest.fixture
    def agent(self):
        return CardioIntelligenceAgent(MagicMock())

    def test_comparative_strategy(self, agent):
        plan = agent.search_plan("Compare TAVR vs SAVR for aortic stenosis")
        assert plan.search_strategy == "comparative"

    def test_clinical_strategy(self, agent):
        plan = agent.search_plan(
            "65 year old patient with history of MI presents with chest pain"
        )
        assert plan.search_strategy == "clinical"

    def test_targeted_strategy_single_condition(self, agent):
        plan = agent.search_plan("Heart failure GDMT optimisation")
        assert plan.search_strategy == "targeted"

    def test_broad_strategy_vague_query(self, agent):
        plan = agent.search_plan("cardiovascular")
        assert plan.search_strategy == "broad"


# =====================================================================
# search_plan() SUB-QUESTIONS & TOPICS
# =====================================================================


class TestSearchPlanSubQuestions:
    """search_plan() should generate sub-questions."""

    @pytest.fixture
    def agent(self):
        return CardioIntelligenceAgent(MagicMock())

    def test_sub_questions_generated(self, agent):
        plan = agent.search_plan("Heart failure management")
        assert isinstance(plan.sub_questions, list)
        assert len(plan.sub_questions) > 0

    def test_identified_topics_populated(self, agent):
        plan = agent.search_plan("Heart failure with amiodarone and MYH7 variant")
        assert isinstance(plan.identified_topics, list)
        assert len(plan.identified_topics) > 0


# =====================================================================
# CardioIntelligenceAgent.run() WITH MOCK RAG
# =====================================================================


class TestCardioIntelligenceAgentRun:
    """run() should call plan -> search -> evaluate -> synthesize."""

    @pytest.fixture
    def agent(self):
        mock_rag = MagicMock()
        mock_response = MagicMock()
        mock_response.results = []
        mock_response.answer = "Test answer"
        mock_rag.query.return_value = mock_response
        mock_rag.search.return_value = []
        return CardioIntelligenceAgent(mock_rag)

    def test_run_returns_response(self, agent):
        result = agent.run("What is GDMT for HFrEF?")
        assert result is not None

    def test_run_calls_rag_query(self, agent):
        agent.run("Heart failure management")
        agent.rag.query.assert_called_once()

    def test_run_with_workflow_override(self, agent):
        agent.run(
            "General question",
            workflow_override=CardioWorkflowType.HEART_FAILURE,
        )
        call_kwargs = agent.rag.query.call_args
        assert call_kwargs[1]["workflow"] == CardioWorkflowType.HEART_FAILURE

    def test_run_with_patient_context(self, agent):
        ctx = {"age": 65, "sex": "male", "lvef": 30}
        agent.run("HFrEF management", patient_context=ctx)
        call_kwargs = agent.rag.query.call_args
        assert call_kwargs[1]["patient_context"] == ctx


# =====================================================================
# _evaluate_evidence_quality()
# =====================================================================


class TestEvaluateEvidenceQuality:
    """_evaluate_evidence_quality should classify evidence quality."""

    @pytest.fixture
    def agent(self):
        return CardioIntelligenceAgent(MagicMock())

    def test_no_results_is_insufficient(self, agent):
        assert agent._evaluate_evidence_quality([]) == "insufficient"

    def test_none_results_is_insufficient(self, agent):
        assert agent._evaluate_evidence_quality(None) == "insufficient"

    def test_sufficient_with_many_collections(self, agent):
        results = []
        for i in range(12):
            r = MagicMock()
            r.collection = f"collection_{i % 4}"
            results.append(r)
        assert agent._evaluate_evidence_quality(results) == "sufficient"

    def test_partial_with_moderate_results(self, agent):
        results = []
        for i in range(6):
            r = MagicMock()
            r.collection = f"collection_{i % 2}"
            results.append(r)
        assert agent._evaluate_evidence_quality(results) == "partial"

    def test_insufficient_with_few_results(self, agent):
        results = [MagicMock(collection="one")]
        assert agent._evaluate_evidence_quality(results) == "insufficient"


# =====================================================================
# generate_report()
# =====================================================================


class TestGenerateReport:
    """generate_report() should produce a markdown report string."""

    @pytest.fixture
    def agent(self):
        mock_rag = MagicMock()
        mock_response = MagicMock()
        mock_response.answer = "Sample analysis"
        mock_response.results = []
        mock_rag.query.return_value = mock_response
        return CardioIntelligenceAgent(mock_rag)

    def _mock_response(self):
        resp = MagicMock()
        resp.answer = "Sample analysis text"
        return resp

    def test_returns_string(self, agent):
        report = agent.generate_report("Heart failure?", self._mock_response())
        assert isinstance(report, str)

    def test_contains_report_header(self, agent):
        report = agent.generate_report("Heart failure?", self._mock_response())
        assert "Cardiology Intelligence Report" in report

    def test_contains_query(self, agent):
        report = agent.generate_report("Heart failure management?", self._mock_response())
        assert "Heart failure management?" in report

    def test_critical_condition_flagged(self, agent):
        report = agent.generate_report("Acute STEMI management", self._mock_response())
        assert "CRITICAL" in report
