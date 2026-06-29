"""Tests for RAG engine functionality.

Since the clinical trial agent does not yet have a full rag_engine.py,
these tests validate the data flow patterns that the RAG engine will use:
  - IngestRecord -> parse -> validate pipeline
  - Search result model creation
  - Response model composition

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    TrialQuery,
    TrialSearchResult,
    TrialResponse,
    WorkflowResult,
    TrialWorkflowType,
    SearchPlan,
)


class TestRAGSearchResultCreation:
    """Test search result creation for RAG pipeline."""

    def test_single_result(self):
        r = TrialSearchResult(
            collection="trial_protocols",
            content="KEYNOTE-024 Phase III trial of pembrolizumab",
            score=0.92,
            metadata={"nct_id": "NCT02477436"},
        )
        assert r.score == 0.92
        assert r.metadata["nct_id"] == "NCT02477436"

    def test_multiple_results_sorted(self):
        results = [
            TrialSearchResult(collection="trial_protocols", content="Trial A", score=0.7),
            TrialSearchResult(collection="trial_literature", content="Paper B", score=0.9),
            TrialSearchResult(collection="trial_safety", content="Signal C", score=0.6),
        ]
        sorted_results = sorted(results, key=lambda x: x.score, reverse=True)
        assert sorted_results[0].score == 0.9
        assert sorted_results[0].collection == "trial_literature"

    def test_result_from_different_collections(self):
        collections = [
            "trial_protocols", "trial_eligibility", "trial_endpoints",
            "trial_regulatory", "trial_literature",
        ]
        for coll in collections:
            r = TrialSearchResult(
                collection=coll,
                content=f"Content from {coll}",
                score=0.5,
            )
            assert r.collection == coll


class TestRAGResponseComposition:
    """Test composing TrialResponse from search results."""

    def test_simple_response(self):
        response = TrialResponse(
            answer="Found 3 relevant trials for NSCLC with EGFR mutation.",
            citations=[
                {"collection": "trial_protocols", "score": 0.92, "nct_id": "NCT02477436"},
                {"collection": "trial_protocols", "score": 0.88, "nct_id": "NCT03723655"},
            ],
            confidence=0.85,
        )
        assert len(response.citations) == 2
        assert response.confidence == 0.85

    def test_response_with_workflows(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=["3 trials match patient profile"],
            confidence=0.8,
        )
        response = TrialResponse(
            answer="Patient matching complete.",
            workflow_results=[wr],
            confidence=0.8,
        )
        assert len(response.workflow_results) == 1
        assert response.workflow_results[0].workflow_type == TrialWorkflowType.PATIENT_MATCHING

    def test_empty_response(self):
        response = TrialResponse(
            answer="No matching trials found.",
            confidence=0.1,
        )
        assert response.citations == []
        assert response.matches == []
        assert response.workflow_results == []


class TestSearchPlanCreation:
    """Test search plan creation for query analysis."""

    def test_broad_search(self):
        plan = SearchPlan(
            question="What clinical trials are available?",
            search_strategy="broad",
        )
        assert plan.search_strategy == "broad"

    def test_targeted_search(self):
        plan = SearchPlan(
            question="Find EGFR-mutant NSCLC trials with osimertinib",
            therapeutic_areas=["oncology"],
            drugs=["osimertinib"],
            biomarkers=["EGFR T790M"],
            relevant_workflows=[TrialWorkflowType.PATIENT_MATCHING],
            search_strategy="targeted",
            sub_questions=[
                "What is the EGFR mutation status?",
                "Which trials use osimertinib?",
            ],
        )
        assert plan.search_strategy == "targeted"
        assert len(plan.sub_questions) == 2
        assert "osimertinib" in plan.drugs

    def test_comparative_search(self):
        plan = SearchPlan(
            question="Compare checkpoint inhibitor trials in NSCLC",
            therapeutic_areas=["oncology"],
            drugs=["pembrolizumab", "nivolumab", "atezolizumab"],
            relevant_workflows=[TrialWorkflowType.COMPETITIVE_INTEL],
            search_strategy="comparative",
        )
        assert plan.search_strategy == "comparative"
        assert len(plan.drugs) == 3

    def test_regulatory_search(self):
        plan = SearchPlan(
            question="FDA approval requirements for KRAS G12C inhibitors",
            relevant_workflows=[TrialWorkflowType.REGULATORY_DOCS],
            search_strategy="regulatory",
            identified_topics=["KRAS G12C", "FDA approval", "accelerated approval"],
        )
        assert plan.search_strategy == "regulatory"
        assert len(plan.identified_topics) == 3


class TestConversationMemoryPattern:
    """Test conversation memory pattern for multi-turn queries."""

    def test_query_sequence(self):
        """Simulate a multi-turn conversation."""
        queries = [
            TrialQuery(question="What NSCLC trials are recruiting?"),
            TrialQuery(
                question="Which of those accept EGFR-mutant patients?",
                workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            ),
            TrialQuery(
                question="What are the eligibility criteria for the best match?",
                workflow_type=TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
            ),
        ]

        # Verify each query is valid
        for q in queries:
            assert len(q.question) > 0

        # Second query should have explicit workflow type
        assert queries[1].workflow_type == TrialWorkflowType.PATIENT_MATCHING

    def test_context_preservation(self):
        """Test that patient context can be preserved across queries."""
        patient_ctx = {
            "age": 62,
            "sex": "male",
            "diagnosis": "NSCLC",
            "biomarkers": ["EGFR T790M", "PD-L1 80%"],
        }

        q1 = TrialQuery(
            question="Match this patient to trials",
            patient_context=patient_ctx,
        )
        q2 = TrialQuery(
            question="What about cardiac safety?",
            patient_context=patient_ctx,  # Same context
        )

        assert q1.patient_context == q2.patient_context
