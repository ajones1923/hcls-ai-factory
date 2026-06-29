"""
Integration tests for the CAR-T Intelligence Agent.
=====================================================
Exercises the full agent pipeline without external dependencies (no Milvus,
no LLM API). Uses realistic CAR-T therapy scenarios and validates
cross-module consistency.

These tests verify:
  - Full plan -> search -> evaluate -> synthesize pipeline
  - Search planning with real target antigens and development stages
  - Evidence evaluation across various hit profiles
  - Export functions (Markdown, JSON)
  - Cross-module consistency (agent -> export round-trip)

Author: Adam Jones
Date: March 2026
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.agent import CARTIntelligenceAgent
from src.export import export_json, export_markdown, generate_filename
from src.models import (
    AgentResponse,
    AssayResult,
    AssayType,
    BiomarkerRecord,
    BiomarkerType,
    CARConstruct,
    CARGeneration,
    CARTLiterature,
    CARTStage,
    ClinicalTrial,
    CrossCollectionResult,
    EvidenceLevel,
    FDAStatus,
    ManufacturingRecord,
    ProcessStep,
    SafetyEventType,
    SafetyRecord,
    SearchHit,
    SourceType,
    TrialPhase,
    TrialStatus,
)


# ═══════════════════════════════════════════════════════════════════════════
# Mock RAG engine for integration tests
# ═══════════════════════════════════════════════════════════════════════════


def _make_mock_rag_engine(target: str = "CD19", disease: str = "B-ALL"):
    """Create a mock RAG engine with realistic CAR-T evidence."""
    rag = MagicMock()

    hits = [
        SearchHit(
            collection="cart_literature",
            id="PMID:30275568",
            score=0.93,
            text=f"{target} CAR-T therapy achieves high complete response rates in {disease}.",
            metadata={
                "title": f"{target} CAR-T in {disease}",
                "year": 2024,
                "target_antigen": target,
            },
        ),
        SearchHit(
            collection="cart_trials",
            id="NCT03958656",
            score=0.89,
            text=f"Phase 2 study of {target}-directed CAR-T in relapsed/refractory {disease}.",
            metadata={
                "phase": "Phase 2",
                "status": "Completed",
                "target_antigen": target,
            },
        ),
        SearchHit(
            collection="cart_constructs",
            id=f"construct-{target.lower()}-001",
            score=0.86,
            text=f"2nd-gen CAR with 4-1BB costimulatory domain targeting {target}.",
            metadata={
                "generation": "2nd",
                "costimulatory_domain": "4-1BB",
                "target_antigen": target,
            },
        ),
        SearchHit(
            collection="cart_safety",
            id="safety-crs-001",
            score=0.81,
            text=f"Grade 3+ CRS in 22% of patients receiving {target} CAR-T.",
            metadata={
                "event_type": "CRS",
                "severity_grade": "Grade 3-4",
                "product": f"{target}-CAR-T",
            },
        ),
        SearchHit(
            collection="cart_manufacturing",
            id="mfg-lenti-001",
            score=0.76,
            text="Lentiviral transduction efficiency 45% at MOI 5.",
            metadata={
                "process_step": "transduction",
                "parameter": "MOI",
                "parameter_value": "5",
            },
        ),
        SearchHit(
            collection="cart_biomarkers",
            id="biomarker-peak-expansion",
            score=0.74,
            text="Peak CAR-T expansion at day 10-14 correlates with response.",
            metadata={
                "biomarker_name": "CAR-T peak expansion",
                "biomarker_type": "predictive",
            },
        ),
    ]

    evidence = CrossCollectionResult(
        query=f"{target} CAR-T {disease}",
        hits=hits,
        knowledge_context=(
            f"## Target Antigen: {target}\n"
            f"- **Protein:** {target} surface antigen\n"
            f"- **Diseases:** {disease}\n"
        ),
        total_collections_searched=10,
        search_time_ms=55.2,
    )

    rag.retrieve.return_value = evidence

    def _build_prompt(question, evidence):
        return f"Question: {question}\nEvidence: {evidence.hit_count} items"

    rag._build_prompt.side_effect = _build_prompt

    # LLM generate mock
    rag.llm = MagicMock()
    rag.llm.generate.return_value = (
        f"{target} CAR-T therapy is effective in {disease} with high response rates. "
        f"Key considerations include CRS management and long-term monitoring."
    )

    return rag


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Full Agent Pipeline
# ═══════════════════════════════════════════════════════════════════════════


class TestAgentPipelineIntegration:
    """Test the full plan -> search -> evaluate -> synthesize pipeline."""

    def test_cd19_ball_full_pipeline(self):
        """Full pipeline for CD19 CAR-T in B-ALL."""
        rag = _make_mock_rag_engine("CD19", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("What is the efficacy of CD19 CAR-T in B-ALL?")

        assert isinstance(response, AgentResponse)
        assert "CD19" in response.answer
        assert response.evidence.hit_count > 0
        assert len(response.knowledge_used) > 0

    def test_bcma_myeloma_pipeline(self):
        """Full pipeline for BCMA CAR-T in multiple myeloma."""
        rag = _make_mock_rag_engine("BCMA", "Multiple Myeloma")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("BCMA CAR-T therapy outcomes in multiple myeloma")

        assert isinstance(response, AgentResponse)
        assert response.evidence.hit_count >= 3

    def test_cd22_pipeline(self):
        """Full pipeline for CD22 CAR-T."""
        rag = _make_mock_rag_engine("CD22", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("CD22 targeted CAR-T for relapsed B-ALL")

        assert isinstance(response, AgentResponse)
        assert "CD22" in response.answer

    def test_comparative_query(self):
        """Full pipeline for a comparative query."""
        rag = _make_mock_rag_engine("CD19", "DLBCL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("Compare 4-1BB vs CD28 costimulatory domains for CD19 CAR-T")

        assert isinstance(response, AgentResponse)
        # Comparative queries should have been detected in planning

    def test_failure_analysis_query(self):
        """Full pipeline for a failure analysis question."""
        rag = _make_mock_rag_engine("CD19", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("Why do CD19 CAR-T therapies fail in some patients?")

        assert isinstance(response, AgentResponse)
        assert response.evidence.hit_count >= 0

    def test_manufacturing_query(self):
        """Full pipeline for a manufacturing-focused question."""
        rag = _make_mock_rag_engine("CD19", "DLBCL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run(
            "What manufacturing parameters predict clinical response "
            "for CD19 CAR-T in DLBCL?"
        )

        assert isinstance(response, AgentResponse)

    def test_safety_query(self):
        """Full pipeline for a safety-focused question."""
        rag = _make_mock_rag_engine("CD19", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("CRS management protocols for CD19 CAR-T therapy")

        assert isinstance(response, AgentResponse)


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Search Planning
# ═══════════════════════════════════════════════════════════════════════════


class TestSearchPlanIntegration:
    """Test search planning with realistic CAR-T queries."""

    @pytest.fixture
    def agent(self):
        return CARTIntelligenceAgent(_make_mock_rag_engine())

    def test_cd19_target_detected(self, agent):
        """CD19 is identified as a target antigen."""
        plan = agent.search_plan("CD19 CAR-T efficacy in DLBCL")
        assert "CD19" in plan.target_antigens

    def test_bcma_target_detected(self, agent):
        """BCMA is identified as a target antigen."""
        plan = agent.search_plan("BCMA directed CAR-T for myeloma")
        assert "BCMA" in plan.target_antigens

    def test_clinical_stage_detected(self, agent):
        """Clinical stage keywords are detected."""
        plan = agent.search_plan("Clinical trial results for CAR-T patient survival")
        assert CARTStage.CLINICAL in plan.relevant_stages

    def test_manufacturing_stage_detected(self, agent):
        """Manufacturing keywords trigger vector_eng stage."""
        plan = agent.search_plan("Lentiviral vector manufacturing for CAR-T")
        assert CARTStage.VECTOR_ENG in plan.relevant_stages

    def test_testing_stage_detected(self, agent):
        """In vitro testing keywords trigger testing stage."""
        plan = agent.search_plan("In vitro cytotoxicity assay for CAR-T constructs")
        assert CARTStage.TESTING in plan.relevant_stages

    def test_car_design_stage_detected(self, agent):
        """CAR design keywords trigger car_design stage."""
        plan = agent.search_plan("4-1BB vs CD28 costimulatory domain design")
        assert CARTStage.CAR_DESIGN in plan.relevant_stages

    def test_target_id_stage_detected(self, agent):
        """Antigen expression keywords trigger target_id stage."""
        plan = agent.search_plan("CD19 antigen expression on B-cells")
        assert CARTStage.TARGET_ID in plan.relevant_stages

    def test_comparative_strategy(self, agent):
        """'Compare' keyword triggers comparative strategy."""
        plan = agent.search_plan("Compare CD19 vs BCMA CAR-T therapies")
        assert plan.search_strategy == "comparative"

    def test_targeted_strategy(self, agent):
        """Known antigen + single stage -> targeted strategy."""
        plan = agent.search_plan("CD19 CAR-T clinical outcomes")
        assert plan.search_strategy == "targeted"

    def test_failure_decomposition(self, agent):
        """'Why ... fail' triggers sub-question decomposition."""
        plan = agent.search_plan("Why do CD19 CAR-T therapies fail?")
        assert len(plan.sub_questions) >= 2


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Evidence Evaluation
# ═══════════════════════════════════════════════════════════════════════════


class TestEvidenceEvaluationIntegration:
    """Test evidence evaluation with various hit profiles."""

    @pytest.fixture
    def agent(self):
        return CARTIntelligenceAgent(_make_mock_rag_engine())

    def test_empty_evidence_insufficient(self, agent):
        """Empty evidence is insufficient."""
        evidence = CrossCollectionResult(query="test", hits=[])
        assert agent.evaluate_evidence(evidence) == "insufficient"

    def test_rich_evidence_sufficient(self, agent):
        """Many hits from many collections is sufficient."""
        hits = []
        for col in ["cart_literature", "cart_trials", "cart_constructs",
                     "cart_safety", "cart_manufacturing"]:
            for i in range(3):
                hits.append(SearchHit(
                    collection=col,
                    id=f"{col}-{i}",
                    score=0.8 + i * 0.05,
                    text=f"Evidence from {col}",
                ))
        evidence = CrossCollectionResult(
            query="CD19 CAR-T efficacy",
            hits=hits,
            total_collections_searched=10,
        )
        assert agent.evaluate_evidence(evidence) == "sufficient"

    def test_sparse_evidence_partial_or_insufficient(self, agent):
        """Few hits from one collection is partial or insufficient."""
        evidence = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(
                    collection="cart_literature",
                    id="1",
                    score=0.8,
                    text="Single evidence item",
                ),
            ],
            total_collections_searched=10,
        )
        result = agent.evaluate_evidence(evidence)
        assert result in ("partial", "insufficient")


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Report Generation
# ═══════════════════════════════════════════════════════════════════════════


class TestReportGenerationIntegration:
    """Test report generation from agent response."""

    def test_report_structure(self):
        """Agent generate_report produces well-structured Markdown."""
        rag = _make_mock_rag_engine("CD19", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("CD19 CAR-T efficacy in B-ALL")
        report = agent.generate_report(response)

        assert isinstance(report, str)
        assert "# CAR-T Intelligence Report" in report
        assert "CD19" in report
        assert "Analysis" in report
        assert "Evidence Sources" in report

    def test_report_includes_evidence_collections(self):
        """Report lists evidence from different collections."""
        rag = _make_mock_rag_engine("BCMA", "Myeloma")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("BCMA CAR-T outcomes")
        report = agent.generate_report(response)

        # Report should mention evidence collection names
        assert "results)" in report or "result)" in report

    def test_report_includes_knowledge_graph(self):
        """Report includes knowledge graph section if context available."""
        rag = _make_mock_rag_engine("CD19", "DLBCL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("CD19 CAR-T in DLBCL")
        report = agent.generate_report(response)

        if response.knowledge_used:
            assert "Knowledge Graph" in report


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Export Functions
# ═══════════════════════════════════════════════════════════════════════════


class TestExportIntegration:
    """Test export functions with agent-generated data."""

    @pytest.fixture
    def agent_evidence(self):
        """Evidence from a simulated agent run."""
        return CrossCollectionResult(
            query="CD19 CAR-T efficacy in B-ALL",
            hits=[
                SearchHit(
                    collection="cart_literature",
                    id="PMID:30275568",
                    score=0.93,
                    text="CD19 CAR-T achieves 81% CR in r/r B-ALL.",
                    metadata={"year": 2024, "target_antigen": "CD19"},
                ),
                SearchHit(
                    collection="cart_trials",
                    id="NCT03958656",
                    score=0.89,
                    text="Phase 2 tisagenlecleucel pediatric B-ALL.",
                    metadata={"phase": "Phase 2", "status": "Completed"},
                ),
                SearchHit(
                    collection="cart_safety",
                    id="safety-crs-kymriah",
                    score=0.85,
                    text="Grade 3+ CRS 22%, ICANS 12%.",
                    metadata={"event_type": "CRS", "product": "Kymriah"},
                ),
            ],
            knowledge_context="## Target: CD19\nB-cell antigen",
            total_collections_searched=10,
            search_time_ms=42.5,
        )

    def test_markdown_export(self, agent_evidence):
        """Markdown export produces non-empty report with all sections."""
        md = export_markdown(
            query="CD19 CAR-T efficacy in B-ALL",
            response_text="CD19 CAR-T achieves high CR rates in B-ALL.",
            evidence=agent_evidence,
        )
        assert isinstance(md, str)
        assert len(md) > 100
        assert "CD19" in md
        assert "Evidence Sources" in md
        assert "CAR-T Intelligence" in md

    def test_json_export(self, agent_evidence):
        """JSON export produces valid, parseable JSON."""
        json_str = export_json(
            query="CD19 CAR-T efficacy",
            response_text="CD19 CAR-T is effective.",
            evidence=agent_evidence,
        )
        data = json.loads(json_str)
        assert isinstance(data, dict)
        assert data["query"] == "CD19 CAR-T efficacy"
        assert "evidence" in data
        assert "search_metrics" in data
        assert data["search_metrics"]["total_results"] == 3

    def test_filename_generation(self):
        """Filename generation produces valid filenames."""
        md_name = generate_filename("md")
        json_name = generate_filename("json")
        pdf_name = generate_filename("pdf")

        assert md_name.endswith(".md")
        assert json_name.endswith(".json")
        assert pdf_name.endswith(".pdf")
        assert md_name.startswith("cart_query_")

    def test_markdown_no_evidence(self):
        """Markdown export handles None evidence gracefully."""
        md = export_markdown(
            query="test question",
            response_text="No evidence found.",
            evidence=None,
        )
        assert isinstance(md, str)
        assert "test question" in md


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Model Construction with Realistic Data
# ═══════════════════════════════════════════════════════════════════════════


class TestModelIntegration:
    """Test Pydantic models with realistic CAR-T data."""

    def test_cart_literature_embedding(self):
        """CARTLiterature.to_embedding_text includes key fields."""
        lit = CARTLiterature(
            id="PMID:30275568",
            title="Long-term outcomes of tisagenlecleucel in pediatric B-ALL",
            text_chunk=(
                "Tisagenlecleucel achieved 81% overall remission rate in "
                "pediatric and young adult patients with relapsed/refractory B-ALL."
            ),
            source_type=SourceType.PUBMED,
            year=2024,
            cart_stage=CARTStage.CLINICAL,
            target_antigen="CD19",
            disease="B-ALL",
            keywords="CAR-T, CD19, B-ALL, tisagenlecleucel, pediatric",
        )
        text = lit.to_embedding_text()
        assert "CD19" in text
        assert "tisagenlecleucel" in text

    def test_clinical_trial_model(self):
        """ClinicalTrial model with realistic data."""
        trial = ClinicalTrial(
            id="NCT03958656",
            title="Tisagenlecleucel in Pediatric B-ALL (ELIANA follow-up)",
            text_summary="Long-term follow-up of tisagenlecleucel in pediatric B-ALL.",
            phase=TrialPhase.PHASE_2,
            status=TrialStatus.COMPLETED,
            sponsor="Novartis",
            target_antigen="CD19",
            car_generation=CARGeneration.SECOND,
            costimulatory="4-1BB",
            disease="B-ALL",
            enrollment=75,
            start_year=2019,
            outcome_summary="81% ORR, 3-year OS 63%",
        )
        text = trial.to_embedding_text()
        assert "CD19" in text
        assert "B-ALL" in text
        assert "81%" in text

    def test_car_construct_model(self):
        """CARConstruct model with real product data."""
        construct = CARConstruct(
            id="construct-kymriah",
            name="Tisagenlecleucel (Kymriah)",
            text_summary=(
                "2nd-generation anti-CD19 CAR with 4-1BB costimulatory domain "
                "and CD3-zeta signaling domain. FMC63 scFv, lentiviral vector."
            ),
            target_antigen="CD19",
            scfv_origin="FMC63 (murine)",
            costimulatory_domain="4-1BB",
            signaling_domain="CD3-zeta",
            generation=CARGeneration.SECOND,
            hinge_tm="CD8a hinge and transmembrane",
            vector_type="lentiviral",
            fda_status=FDAStatus.APPROVED,
            known_toxicities="CRS 58%, ICANS 18%, cytopenias",
        )
        text = construct.to_embedding_text()
        assert "CD19" in text
        assert "4-1BB" in text
        assert construct.fda_status == FDAStatus.APPROVED

    def test_safety_record_model(self):
        """SafetyRecord with CRS data."""
        safety = SafetyRecord(
            id="safety-crs-kymriah-001",
            text_summary=(
                "In the ELIANA trial, Grade 3+ CRS was reported in 47% of patients. "
                "Median onset day 3, managed with tocilizumab and corticosteroids."
            ),
            product="Kymriah",
            event_type=SafetyEventType.CRS,
            severity_grade="Grade 3-4",
            onset_timing="Median day 3 post-infusion",
            incidence_rate="47% grade 3+; 77% any grade",
            management_protocol="Tocilizumab 8mg/kg, dexamethasone for refractory",
            outcome="Resolved in >95% of cases",
            reporting_source="ELIANA trial",
            year=2024,
        )
        text = safety.to_embedding_text()
        assert "CRS" in text
        assert "Kymriah" in text

    def test_manufacturing_record_model(self):
        """ManufacturingRecord with transduction data."""
        mfg = ManufacturingRecord(
            id="mfg-lenti-transduction-001",
            text_summary="Lentiviral transduction of T-cells at MOI 5.",
            process_step=ProcessStep.TRANSDUCTION,
            vector_type="lentiviral",
            parameter="MOI",
            parameter_value="5",
            target_spec="VCN 1-5 copies/cell",
            met_spec="yes",
            batch_id="BATCH-2024-001",
        )
        text = mfg.to_embedding_text()
        assert "MOI" in text
        assert "5" in text

    def test_biomarker_record_model(self):
        """BiomarkerRecord with predictive biomarker."""
        biomarker = BiomarkerRecord(
            id="biomarker-peak-expansion",
            text_summary=(
                "Peak CAR-T expansion (Cmax) at day 10-14 post-infusion "
                "is predictive of durable response in CD19 CAR-T therapy."
            ),
            biomarker_name="Peak CAR-T expansion (Cmax)",
            biomarker_type=BiomarkerType.PREDICTIVE,
            assay_method="Flow cytometry (qPCR confirmation)",
            clinical_cutoff=">50 copies/ug genomic DNA",
            predictive_value="Higher Cmax correlates with CR and durability",
            associated_outcome="Complete remission duration",
            target_antigen="CD19",
            disease="B-ALL",
            evidence_level=EvidenceLevel.VALIDATED,
        )
        text = biomarker.to_embedding_text()
        assert "expansion" in text.lower()
        assert biomarker.evidence_level == EvidenceLevel.VALIDATED

    def test_assay_result_model(self):
        """AssayResult with cytotoxicity data."""
        assay = AssayResult(
            id="assay-cytotox-cd19-001",
            text_summary=(
                "CD19 CAR-T cells showed 92% specific lysis of Nalm-6 cells "
                "at 5:1 E:T ratio after 4-hour co-culture."
            ),
            assay_type=AssayType.CYTOTOXICITY,
            construct_id="construct-kymriah",
            target_antigen="CD19",
            cell_line="Nalm-6",
            effector_ratio="5:1",
            key_metric="specific_lysis",
            metric_value=92.0,
            outcome="success",
        )
        text = assay.to_embedding_text()
        assert "Nalm-6" in text
        assert "92.0" in text


# ═══════════════════════════════════════════════════════════════════════════
# Integration: Cross-Module Consistency
# ═══════════════════════════════════════════════════════════════════════════


class TestCrossModuleConsistency:
    """Verify data flows consistently across agent modules."""

    def test_agent_response_to_report(self):
        """Agent response produces a valid report."""
        rag = _make_mock_rag_engine("CD19", "B-ALL")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("CD19 CAR-T efficacy in B-ALL")
        report = agent.generate_report(response)

        # Report should contain evidence from the response
        assert response.evidence.hit_count > 0
        assert isinstance(report, str)
        assert len(report) > 100

    def test_agent_response_to_export(self):
        """Agent response evidence can be exported consistently."""
        rag = _make_mock_rag_engine("BCMA", "Myeloma")
        agent = CARTIntelligenceAgent(rag)

        response = agent.run("BCMA CAR-T outcomes")

        # Export to markdown
        md = export_markdown(
            query=response.question,
            response_text=response.answer,
            evidence=response.evidence,
        )
        assert "BCMA" in md

        # Export to JSON
        json_str = export_json(
            query=response.question,
            response_text=response.answer,
            evidence=response.evidence,
        )
        data = json.loads(json_str)
        assert data["query"] == response.question
        assert data["search_metrics"]["total_results"] == response.evidence.hit_count

    def test_evidence_grouping_consistency(self):
        """hits_by_collection output is consistent with hit_count."""
        evidence = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(collection="cart_literature", id="1", score=0.9, text="A"),
                SearchHit(collection="cart_literature", id="2", score=0.8, text="B"),
                SearchHit(collection="cart_trials", id="3", score=0.7, text="C"),
                SearchHit(collection="cart_safety", id="4", score=0.6, text="D"),
                SearchHit(collection="cart_constructs", id="5", score=0.5, text="E"),
            ],
        )
        grouped = evidence.hits_by_collection()
        total_grouped = sum(len(v) for v in grouped.values())
        assert total_grouped == evidence.hit_count
        assert len(grouped["cart_literature"]) == 2
        assert len(grouped) == 4  # 4 unique collections
