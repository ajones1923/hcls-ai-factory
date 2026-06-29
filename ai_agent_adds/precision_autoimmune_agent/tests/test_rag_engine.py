"""
Tests for src/rag_engine.py — AutoimmuneRAGEngine

All external dependencies (pymilvus, anthropic, sentence_transformers) are mocked
so tests run without any external services.
"""

from __future__ import annotations

import threading
from unittest.mock import MagicMock, patch

import pytest

from src.rag_engine import (
    _COMPARATIVE_RE,
    _SAFE_FILTER_RE,
    SYSTEM_PROMPT,
    AutoimmuneRAGEngine,
    CrossCollectionResult,
    SearchHit,
    _sanitize_filter_value,
)

# ── Fixtures ─────────────────────────────────────────────────────────────


class FakeSettings:
    """Minimal settings stand-in for tests (no pydantic dependency)."""

    MILVUS_HOST = "localhost"
    MILVUS_PORT = 19530
    EMBEDDING_MODEL = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIM = 384
    BGE_INSTRUCTION = "Represent this sentence for searching relevant passages: "
    LLM_MODEL = "claude-sonnet-4-20250514"
    LLM_MAX_TOKENS = 4096
    LLM_TEMPERATURE = 0.2
    TOP_K_PER_COLLECTION = 5
    SCORE_THRESHOLD = 0.40
    MAX_EVIDENCE_ITEMS = 30
    CONVERSATION_MEMORY_SIZE = 3
    CITATION_HIGH = 0.80
    CITATION_MEDIUM = 0.60
    MAX_EVIDENCE_TEXT_LENGTH = 1500
    MAX_KNOWLEDGE_CONTEXT_ITEMS = 25

    COLL_CLINICAL_DOCUMENTS = "autoimmune_clinical_documents"
    COLL_PATIENT_LABS = "autoimmune_patient_labs"
    COLL_PATIENT_TIMELINES = "autoimmune_patient_timelines"
    COLL_LITERATURE = "autoimmune_literature"

    @property
    def collection_config(self):
        return {
            "autoimmune_clinical_documents": {"weight": 0.18, "label": "Clinical Document"},
            "autoimmune_patient_labs": {"weight": 0.14, "label": "Lab Result"},
            "autoimmune_autoantibody_panels": {"weight": 0.12, "label": "Autoantibody"},
            "autoimmune_hla_associations": {"weight": 0.08, "label": "HLA Association"},
            "autoimmune_disease_criteria": {"weight": 0.08, "label": "Classification Criteria"},
            "autoimmune_disease_activity": {"weight": 0.07, "label": "Disease Activity"},
            "autoimmune_flare_patterns": {"weight": 0.06, "label": "Flare Pattern"},
            "autoimmune_biologic_therapies": {"weight": 0.06, "label": "Biologic Therapy"},
            "autoimmune_pgx_rules": {"weight": 0.04, "label": "PGx Rule"},
            "autoimmune_clinical_trials": {"weight": 0.05, "label": "Clinical Trial"},
            "autoimmune_literature": {"weight": 0.05, "label": "Literature"},
            "autoimmune_patient_timelines": {"weight": 0.03, "label": "Timeline"},
            "autoimmune_cross_disease": {"weight": 0.02, "label": "Cross-Disease"},
            "genomic_evidence": {"weight": 0.02, "label": "Genomic Evidence"},
        }

    @property
    def all_collection_names(self):
        return list(self.collection_config.keys())


@pytest.fixture
def fake_settings():
    return FakeSettings()


@pytest.fixture
def mock_embedder():
    emb = MagicMock()
    emb.encode.return_value = MagicMock(tolist=MagicMock(return_value=[0.1] * 384))
    return emb


@pytest.fixture
def mock_llm():
    llm = MagicMock()
    # .messages.create(...) returns a response with .content[0].text
    mock_response = MagicMock()
    mock_response.content = [MagicMock(text="LLM synthesized answer.")]
    llm.messages.create.return_value = mock_response
    return llm


@pytest.fixture
def mock_cm():
    cm = MagicMock()
    cm.search_all.return_value = {}
    return cm


@pytest.fixture
def engine(mock_cm, mock_embedder, mock_llm, fake_settings):
    return AutoimmuneRAGEngine(
        collection_manager=mock_cm,
        embedder=mock_embedder,
        llm_client=mock_llm,
        settings=fake_settings,
        knowledge=True,  # Truthy so _build_knowledge_context doesn't short-circuit
    )


# ── Sanitize filter value ────────────────────────────────────────────────


class TestSanitizeFilterValue:
    def test_valid_alphanumeric(self):
        assert _sanitize_filter_value("P001") == "P001"

    def test_valid_with_hyphens_underscores(self):
        assert _sanitize_filter_value("patient-001_v2") == "patient-001_v2"

    def test_valid_with_dots(self):
        assert _sanitize_filter_value("test.value") == "test.value"

    def test_empty_string_returns_none(self):
        assert _sanitize_filter_value("") is None

    def test_too_long_returns_none(self):
        assert _sanitize_filter_value("x" * 100, max_len=64) is None

    def test_exactly_max_len_passes(self):
        assert _sanitize_filter_value("x" * 64, max_len=64) == "x" * 64

    def test_sql_injection_rejected(self):
        assert _sanitize_filter_value('"; DROP TABLE') is None

    def test_special_chars_rejected(self):
        assert _sanitize_filter_value("value'OR 1=1") is None

    def test_unicode_rejected(self):
        assert _sanitize_filter_value("test\u00e9value") is None

    def test_newline_rejected(self):
        assert _sanitize_filter_value("hello\nworld") is None

    def test_custom_max_len(self):
        assert _sanitize_filter_value("short", max_len=3) is None
        assert _sanitize_filter_value("hi", max_len=3) == "hi"


# ── Disease area detection ───────────────────────────────────────────────


class TestDetectDiseaseAreas:
    def test_rheumatoid_arthritis_keywords(self, engine):
        assert "rheumatoid_arthritis" in engine._detect_disease_areas("patient with rheumatoid arthritis")
        assert "rheumatoid_arthritis" in engine._detect_disease_areas("anti-CCP results are elevated")
        assert "rheumatoid_arthritis" in engine._detect_disease_areas("DAS28 score is 5.3")
        assert "rheumatoid_arthritis" in engine._detect_disease_areas("morning stiffness lasting 2 hours")

    def test_systemic_lupus_keywords(self, engine):
        assert "systemic_lupus" in engine._detect_disease_areas("patient diagnosed with lupus")
        assert "systemic_lupus" in engine._detect_disease_areas("anti-dsDNA is elevated")
        assert "systemic_lupus" in engine._detect_disease_areas("SLEDAI score 12")
        assert "systemic_lupus" in engine._detect_disease_areas("butterfly rash on face")

    def test_multiple_sclerosis_keywords(self, engine):
        assert "multiple_sclerosis" in engine._detect_disease_areas("multiple sclerosis relapse")
        assert "multiple_sclerosis" in engine._detect_disease_areas("oligoclonal bands detected")

    def test_ankylosing_spondylitis_keywords(self, engine):
        assert "ankylosing_spondylitis" in engine._detect_disease_areas("ankylosing spondylitis")
        assert "ankylosing_spondylitis" in engine._detect_disease_areas("HLA-B27 positive")
        assert "ankylosing_spondylitis" in engine._detect_disease_areas("sacroiliitis on imaging")

    def test_sjogrens_keywords(self, engine):
        assert "sjogrens" in engine._detect_disease_areas("sjogren syndrome")
        assert "sjogrens" in engine._detect_disease_areas("anti-SSA positive")
        assert "sjogrens" in engine._detect_disease_areas("sicca symptoms")

    def test_systemic_sclerosis_keywords(self, engine):
        assert "systemic_sclerosis" in engine._detect_disease_areas("scleroderma diagnosis")
        assert "systemic_sclerosis" in engine._detect_disease_areas("anti-Scl-70 positive")
        assert "systemic_sclerosis" in engine._detect_disease_areas("raynaud phenomenon")

    def test_ibd_keywords(self, engine):
        assert "inflammatory_bowel" in engine._detect_disease_areas("crohn disease")
        assert "inflammatory_bowel" in engine._detect_disease_areas("ulcerative colitis")
        assert "inflammatory_bowel" in engine._detect_disease_areas("calprotectin elevated")

    def test_psoriasis_keywords(self, engine):
        assert "psoriasis" in engine._detect_disease_areas("psoriasis treatment")
        assert "psoriasis" in engine._detect_disease_areas("IL-17 pathway")

    def test_myasthenia_gravis_keywords(self, engine):
        assert "myasthenia_gravis" in engine._detect_disease_areas("myasthenia gravis")
        assert "myasthenia_gravis" in engine._detect_disease_areas("AChR antibody test")

    def test_celiac_keywords(self, engine):
        assert "celiac" in engine._detect_disease_areas("celiac disease screening")
        assert "celiac" in engine._detect_disease_areas("anti-tTG IgA elevated")
        assert "celiac" in engine._detect_disease_areas("villous atrophy found")

    def test_thyroid_keywords(self, engine):
        assert "thyroid_autoimmune" in engine._detect_disease_areas("graves disease")
        assert "thyroid_autoimmune" in engine._detect_disease_areas("hashimoto thyroiditis")
        assert "thyroid_autoimmune" in engine._detect_disease_areas("anti-TPO positive")

    def test_pots_eds_mcas_keywords(self, engine):
        assert "pots_eds_mcas" in engine._detect_disease_areas("POTS diagnosis")
        assert "pots_eds_mcas" in engine._detect_disease_areas("hypermobility type EDS")
        assert "pots_eds_mcas" in engine._detect_disease_areas("mast cell activation")

    def test_multiple_diseases_detected(self, engine):
        areas = engine._detect_disease_areas("patient with lupus and sjogren syndrome overlap")
        assert "systemic_lupus" in areas
        assert "sjogrens" in areas

    def test_no_disease_detected(self, engine):
        areas = engine._detect_disease_areas("general health checkup")
        assert areas == []

    def test_case_insensitive(self, engine):
        assert "rheumatoid_arthritis" in engine._detect_disease_areas("RHEUMATOID ARTHRITIS")


# ── Comparative query detection ──────────────────────────────────────────


class TestIsComparative:
    def test_compare_keyword(self, engine):
        assert engine._is_comparative("Compare adalimumab vs etanercept") is True

    def test_versus_keyword(self, engine):
        assert engine._is_comparative("adalimumab versus rituximab") is True

    def test_vs_abbreviation(self, engine):
        assert engine._is_comparative("TNF inhibitor vs JAK inhibitor") is True

    def test_vs_with_period(self, engine):
        assert engine._is_comparative("adalimumab vs. etanercept") is True

    def test_difference_between(self, engine):
        assert engine._is_comparative("What is the difference between RA and SLE?") is True

    def test_head_to_head(self, engine):
        assert engine._is_comparative("head-to-head trial data") is True
        assert engine._is_comparative("head to head comparison") is True

    def test_which_is_better(self, engine):
        assert engine._is_comparative("Which is better for RA?") is True

    def test_pros_and_cons(self, engine):
        assert engine._is_comparative("pros and cons of biologics") is True

    def test_not_comparative(self, engine):
        assert engine._is_comparative("What is the DAS28 score for this patient?") is False

    def test_empty_string(self, engine):
        assert engine._is_comparative("") is False


# ── Score relevance ──────────────────────────────────────────────────────


class TestScoreRelevance:
    def test_high_relevance(self, engine):
        assert engine._score_relevance(0.95) == "high"
        assert engine._score_relevance(0.80) == "high"

    def test_medium_relevance(self, engine):
        assert engine._score_relevance(0.75) == "medium"
        assert engine._score_relevance(0.60) == "medium"

    def test_low_relevance(self, engine):
        assert engine._score_relevance(0.50) == "low"
        assert engine._score_relevance(0.10) == "low"

    def test_boundary_values(self, engine):
        # 0.80 is high (>=), 0.60 is medium (>=)
        assert engine._score_relevance(0.80) == "high"
        assert engine._score_relevance(0.799) == "medium"
        assert engine._score_relevance(0.60) == "medium"
        assert engine._score_relevance(0.599) == "low"


# ── Build evidence block ────────────────────────────────────────────────


class TestBuildEvidenceBlock:
    def test_empty_hits(self, engine):
        result = CrossCollectionResult(query="test")
        block = engine._build_evidence_block(result)
        assert block == ""

    def test_single_hit_formatting(self, engine):
        hit = SearchHit(
            collection="autoimmune_patient_labs",
            id="lab_001",
            score=0.92,
            text="CRP elevated at 45 mg/L",
            metadata={},
            relevance="high",
        )
        result = CrossCollectionResult(query="test", hits=[hit])
        block = engine._build_evidence_block(result)

        assert "Evidence 1" in block
        assert "Lab Result" in block  # label from config
        assert "lab_001" in block
        assert "0.920" in block
        assert "[high relevance]" in block
        assert "CRP elevated at 45 mg/L" in block

    def test_hit_with_pmid_includes_pubmed_link(self, engine):
        hit = SearchHit(
            collection="autoimmune_literature",
            id="lit_001",
            score=0.85,
            text="Study results",
            metadata={"pmid": "12345678"},
            relevance="high",
        )
        result = CrossCollectionResult(query="test", hits=[hit])
        block = engine._build_evidence_block(result)
        assert "https://pubmed.ncbi.nlm.nih.gov/12345678/" in block

    def test_hit_with_nct_id_includes_clinicaltrials_link(self, engine):
        hit = SearchHit(
            collection="autoimmune_clinical_trials",
            id="trial_001",
            score=0.78,
            text="Phase 3 trial",
            metadata={"nct_id": "NCT12345678"},
            relevance="medium",
        )
        result = CrossCollectionResult(query="test", hits=[hit])
        block = engine._build_evidence_block(result)
        assert "https://clinicaltrials.gov/study/NCT12345678" in block

    def test_text_truncation(self, engine):
        hit = SearchHit(
            collection="autoimmune_patient_labs",
            id="lab_001",
            score=0.9,
            text="X" * 3000,
            metadata={},
            relevance="high",
        )
        result = CrossCollectionResult(query="test", hits=[hit])
        block = engine._build_evidence_block(result)
        # Text should be truncated to 1500 chars
        x_count = block.count("X")
        assert x_count == 1500

    def test_multiple_hits_numbered(self, engine):
        hits = [
            SearchHit(collection="autoimmune_patient_labs", id=f"id_{i}",
                      score=0.9 - i * 0.1, text=f"text {i}", metadata={}, relevance="high")
            for i in range(3)
        ]
        result = CrossCollectionResult(query="test", hits=hits)
        block = engine._build_evidence_block(result)
        assert "Evidence 1" in block
        assert "Evidence 2" in block
        assert "Evidence 3" in block


# ── Build knowledge context ─────────────────────────────────────────────


class TestBuildKnowledgeContext:
    def test_hla_allele_in_query(self, engine):
        with patch("src.knowledge.HLA_DISEASE_ASSOCIATIONS", {
            "HLA-B*27:05": [{"disease": "ankylosing_spondylitis", "odds_ratio": 87.4, "pmid": "25603694", "note": "Strongest"}],
        }), patch("src.knowledge.AUTOANTIBODY_DISEASE_MAP", {}), \
             patch("src.knowledge.BIOLOGIC_THERAPIES", []), \
             patch("src.knowledge.FLARE_BIOMARKER_PATTERNS", {}):
            ctx = engine._build_knowledge_context("HLA-B*27:05 positive patient", ["ankylosing_spondylitis"])
        assert "HLA-B*27:05" in ctx
        assert "Ankylosing Spondylitis" in ctx
        assert "87.4" in ctx

    def test_autoantibody_in_query(self, engine):
        with patch("src.knowledge.HLA_DISEASE_ASSOCIATIONS", {}), \
             patch("src.knowledge.AUTOANTIBODY_DISEASE_MAP", {
                 "anti-CCP": [{"disease": "rheumatoid_arthritis", "sensitivity": 0.67, "specificity": 0.95, "note": ""}],
             }), \
             patch("src.knowledge.BIOLOGIC_THERAPIES", []), \
             patch("src.knowledge.FLARE_BIOMARKER_PATTERNS", {}):
            ctx = engine._build_knowledge_context("anti-CCP results elevated", ["rheumatoid_arthritis"])
        assert "anti-CCP" in ctx
        assert "Rheumatoid Arthritis" in ctx
        assert "0.67" in ctx

    def test_therapy_in_query(self, engine):
        with patch("src.knowledge.HLA_DISEASE_ASSOCIATIONS", {}), \
             patch("src.knowledge.AUTOANTIBODY_DISEASE_MAP", {}), \
             patch("src.knowledge.BIOLOGIC_THERAPIES", [
                 {
                     "drug_name": "Adalimumab",
                     "drug_class": "TNF inhibitor",
                     "mechanism": "Anti-TNF-alpha",
                     "indicated_diseases": ["rheumatoid_arthritis"],
                     "pgx_considerations": ["HLA-DRB1*03:01 associated with anti-drug antibody formation"],
                 },
             ]), \
             patch("src.knowledge.FLARE_BIOMARKER_PATTERNS", {}):
            ctx = engine._build_knowledge_context("adalimumab dosing for RA", ["rheumatoid_arthritis"])
        assert "Adalimumab" in ctx
        assert "TNF inhibitor" in ctx
        assert "Rheumatoid Arthritis" in ctx

    def test_flare_pattern_in_query(self, engine):
        with patch("src.knowledge.HLA_DISEASE_ASSOCIATIONS", {}), \
             patch("src.knowledge.AUTOANTIBODY_DISEASE_MAP", {}), \
             patch("src.knowledge.BIOLOGIC_THERAPIES", []), \
             patch("src.knowledge.FLARE_BIOMARKER_PATTERNS", {
                 "rheumatoid_arthritis": {
                     "early_warning_biomarkers": ["CRP", "ESR", "IL-6"],
                     "protective_signals": ["stable_RF_titer", "normal_CRP_trend"],
                 },
             }):
            ctx = engine._build_knowledge_context("rheumatoid arthritis flare", ["rheumatoid_arthritis"])
        assert "Rheumatoid Arthritis" in ctx
        assert "CRP" in ctx

    def test_no_knowledge_returns_empty(self, mock_cm, mock_embedder, mock_llm, fake_settings):
        engine = AutoimmuneRAGEngine(
            collection_manager=mock_cm,
            embedder=mock_embedder,
            llm_client=mock_llm,
            settings=fake_settings,
            knowledge=None,
        )
        result = engine._build_knowledge_context("any query", ["rheumatoid_arthritis"])
        assert result == ""


# ── Embedding cache ──────────────────────────────────────────────────────


class TestEmbeddingCache:
    def test_cache_stores_result(self, engine):
        vec = engine._embed("test query")
        assert vec == [0.1] * 384
        assert "test query" in engine._embed_cache

    def test_cache_returns_cached_result(self, engine):
        engine._embed("test query")
        engine.embedder.encode.reset_mock()
        vec = engine._embed("test query")
        assert vec == [0.1] * 384
        engine.embedder.encode.assert_not_called()

    def test_cache_eviction_when_full(self, engine):
        engine._embed_cache_max = 3
        engine._embed("query1")
        engine._embed("query2")
        engine._embed("query3")
        assert len(engine._embed_cache) == 3

        engine._embed("query4")
        assert len(engine._embed_cache) == 3
        # First entry should have been evicted
        assert "query1" not in engine._embed_cache
        assert "query4" in engine._embed_cache

    def test_cache_key_truncated_to_512(self, engine):
        long_query = "A" * 1000
        engine._embed(long_query)
        cache_key = long_query[:512]
        assert cache_key in engine._embed_cache
        assert long_query not in engine._embed_cache

    def test_embed_without_embedder_raises(self, mock_cm, mock_llm, fake_settings):
        engine = AutoimmuneRAGEngine(
            collection_manager=mock_cm,
            embedder=None,
            llm_client=mock_llm,
            settings=fake_settings,
        )
        with pytest.raises(RuntimeError, match="Embedder not initialized"):
            engine._embed("test")


# ── Conversation history (thread-safe) ───────────────────────────────────


class TestConversationHistory:
    def test_clear_history(self, engine):
        with engine._conversation_lock:
            engine._conversation_history.append({"question": "q1", "answer": "a1"})
        assert len(engine._conversation_history) == 1
        engine.clear_history()
        assert len(engine._conversation_history) == 0

    def test_history_max_size(self, engine):
        # maxlen is CONVERSATION_MEMORY_SIZE = 3
        for i in range(5):
            with engine._conversation_lock:
                engine._conversation_history.append({"question": f"q{i}", "answer": f"a{i}"})
        assert len(engine._conversation_history) == 3
        # Only last 3 should remain
        assert engine._conversation_history[0]["question"] == "q2"

    def test_thread_safe_concurrent_access(self, engine):
        """Verify concurrent writes don't crash (basic thread safety test)."""
        errors = []

        def writer(n):
            try:
                for i in range(50):
                    with engine._conversation_lock:
                        engine._conversation_history.append(
                            {"question": f"thread{n}_q{i}", "answer": f"a{i}"}
                        )
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=writer, args=(t,)) for t in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert errors == []
        # Should be capped at maxlen
        assert len(engine._conversation_history) <= engine.settings.CONVERSATION_MEMORY_SIZE


# ── Service readiness properties ─────────────────────────────────────────


class TestServiceReadiness:
    def test_is_ready_true(self, engine):
        assert engine.is_ready is True

    def test_is_ready_false_no_embedder(self, mock_cm, mock_llm, fake_settings):
        e = AutoimmuneRAGEngine(mock_cm, None, mock_llm, fake_settings)
        assert e.is_ready is False

    def test_is_ready_false_no_llm(self, mock_cm, mock_embedder, fake_settings):
        e = AutoimmuneRAGEngine(mock_cm, mock_embedder, None, fake_settings)
        assert e.is_ready is False

    def test_can_search_true(self, engine):
        assert engine.can_search is True

    def test_can_search_false(self, mock_cm, mock_llm, fake_settings):
        e = AutoimmuneRAGEngine(mock_cm, None, mock_llm, fake_settings)
        assert e.can_search is False


# ── Query (end-to-end with mocks) ───────────────────────────────────────


class TestQuery:
    def test_query_not_ready_returns_error(self, mock_cm, mock_embedder, fake_settings):
        engine = AutoimmuneRAGEngine(mock_cm, mock_embedder, None, fake_settings)
        result = engine.query("What is DAS28?")
        assert "Error" in result
        assert "not fully initialized" in result

    def test_query_calls_llm(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {
            "autoimmune_patient_labs": [
                {"id": "lab1", "score": 0.9, "text_chunk": "CRP 45 mg/L"},
            ],
        }
        result = engine.query("What is the CRP level?")
        assert result == "LLM synthesized answer."
        mock_llm.messages.create.assert_called_once()

    def test_query_stores_conversation_history(self, engine, mock_cm):
        mock_cm.search_all.return_value = {}
        engine.query("test question")
        assert len(engine._conversation_history) == 1
        assert engine._conversation_history[0]["question"] == "test question"
        assert engine._conversation_history[0]["answer"] == "LLM synthesized answer."

    def test_query_with_patient_context(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}
        engine.query("DAS28 score?", patient_context="45-year-old female with RA")
        call_kwargs = mock_llm.messages.create.call_args
        messages = call_kwargs[1]["messages"] if "messages" in call_kwargs[1] else call_kwargs[0][0]
        # Find the user message that contains the patient context
        user_msgs = [m for m in messages if m["role"] == "user"]
        assert any("45-year-old female with RA" in m["content"] for m in user_msgs)

    def test_query_llm_error_fallback_with_hits(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {
            "autoimmune_patient_labs": [
                {"id": "lab1", "score": 0.9, "text_chunk": "CRP elevated"},
            ],
        }
        mock_llm.messages.create.side_effect = Exception("API error")
        result = engine.query("CRP levels?")
        assert "LLM synthesis unavailable" in result
        assert "CRP elevated" in result

    def test_query_llm_error_no_hits(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}
        mock_llm.messages.create.side_effect = Exception("API error")
        result = engine.query("anything")
        assert "Error generating response" in result

    def test_query_passes_system_prompt(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}
        engine.query("test")
        call_kwargs = mock_llm.messages.create.call_args[1]
        assert call_kwargs["system"] == SYSTEM_PROMPT

    def test_query_passes_model_and_max_tokens(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}
        engine.query("test")
        call_kwargs = mock_llm.messages.create.call_args[1]
        assert call_kwargs["model"] == "claude-sonnet-4-20250514"
        assert call_kwargs["max_tokens"] == 4096

    def test_query_includes_conversation_history_in_messages(self, engine, mock_cm, mock_llm):
        # Add prior conversation
        with engine._conversation_lock:
            engine._conversation_history.append({"question": "prior Q", "answer": "prior A"})

        mock_cm.search_all.return_value = {}
        engine.query("follow up question")

        call_kwargs = mock_llm.messages.create.call_args[1]
        messages = call_kwargs["messages"]
        # Should have history + current question
        assert len(messages) >= 3  # prior user, prior assistant, current user
        assert messages[0]["role"] == "user"
        assert "prior Q" in messages[0]["content"]
        assert messages[1]["role"] == "assistant"
        assert "prior A" in messages[1]["content"]


# ── Query stream ─────────────────────────────────────────────────────────


class TestQueryStream:
    def test_stream_not_ready_yields_error(self, mock_cm, mock_embedder, fake_settings):
        engine = AutoimmuneRAGEngine(mock_cm, mock_embedder, None, fake_settings)
        tokens = list(engine.query_stream("test"))
        assert any("Error" in t for t in tokens)

    def test_stream_yields_tokens(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}

        # Mock the streaming context manager
        mock_stream = MagicMock()
        mock_stream.text_stream = iter(["Hello", " world", "!"])
        mock_stream.__enter__ = MagicMock(return_value=mock_stream)
        mock_stream.__exit__ = MagicMock(return_value=False)
        mock_llm.messages.stream.return_value = mock_stream

        tokens = list(engine.query_stream("test question"))
        assert tokens == ["Hello", " world", "!"]

    def test_stream_stores_conversation_history(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}

        mock_stream = MagicMock()
        mock_stream.text_stream = iter(["Hello", " world"])
        mock_stream.__enter__ = MagicMock(return_value=mock_stream)
        mock_stream.__exit__ = MagicMock(return_value=False)
        mock_llm.messages.stream.return_value = mock_stream

        list(engine.query_stream("streamed question"))
        assert len(engine._conversation_history) == 1
        assert engine._conversation_history[0]["answer"] == "Hello world"

    def test_stream_error_yields_error_message(self, engine, mock_cm, mock_llm):
        mock_cm.search_all.return_value = {}

        mock_stream = MagicMock()
        mock_stream.text_stream = iter([])
        mock_stream.__enter__ = MagicMock(side_effect=ConnectionError("lost connection"))
        mock_stream.__exit__ = MagicMock(return_value=False)
        mock_llm.messages.stream.return_value = mock_stream

        tokens = list(engine.query_stream("test"))
        combined = "".join(tokens)
        assert "Streaming interrupted" in combined
        assert "ConnectionError" in combined


# ── Retrieve ─────────────────────────────────────────────────────────────


class TestRetrieve:
    def test_retrieve_without_embedder(self, mock_cm, mock_llm, fake_settings):
        engine = AutoimmuneRAGEngine(mock_cm, None, mock_llm, fake_settings)
        result = engine.retrieve("test query")
        assert isinstance(result, CrossCollectionResult)
        assert result.hits == []

    def test_retrieve_with_results(self, engine, mock_cm):
        mock_cm.search_all.return_value = {
            "autoimmune_patient_labs": [
                {"id": "lab1", "score": 0.9, "text_chunk": "CRP 45"},
                {"id": "lab2", "score": 0.7, "text_chunk": "ESR 30"},
            ],
        }
        result = engine.retrieve("CRP and ESR levels")
        assert len(result.hits) == 2
        assert result.hits[0].score >= result.hits[1].score  # sorted desc
        assert result.total_collections_searched == 1
        assert result.search_time_ms >= 0  # may round to 0 on fast machines

    def test_retrieve_deduplicates_by_id(self, engine, mock_cm):
        mock_cm.search_all.return_value = {
            "coll_a": [{"id": "dup1", "score": 0.9, "text_chunk": "text A"}],
            "coll_b": [{"id": "dup1", "score": 0.8, "text_chunk": "text B"}],
        }
        result = engine.retrieve("test")
        ids = [h.id for h in result.hits]
        assert ids.count("dup1") == 1

    def test_retrieve_deduplicates_by_text_hash(self, engine, mock_cm):
        mock_cm.search_all.return_value = {
            "coll_a": [{"id": "id1", "score": 0.9, "text_chunk": "same text content"}],
            "coll_b": [{"id": "id2", "score": 0.8, "text_chunk": "same text content"}],
        }
        result = engine.retrieve("test")
        assert len(result.hits) == 1

    def test_retrieve_caps_at_max_evidence(self, engine, mock_cm):
        engine.settings.MAX_EVIDENCE_ITEMS = 2
        mock_cm.search_all.return_value = {
            "coll_a": [
                {"id": f"id{i}", "score": 0.9 - i * 0.01, "text_chunk": f"unique text {i}"}
                for i in range(5)
            ],
        }
        result = engine.retrieve("test")
        assert len(result.hits) <= 2

    def test_retrieve_with_patient_id_filter(self, engine, mock_cm):
        engine.retrieve("test query", patient_id="P001")
        call_kwargs = mock_cm.search_all.call_args[1]
        filter_exprs = call_kwargs.get("filter_exprs", {})
        assert "autoimmune_clinical_documents" in filter_exprs
        assert 'patient_id == "P001"' in filter_exprs["autoimmune_clinical_documents"]

    def test_retrieve_rejects_unsafe_patient_id(self, engine, mock_cm):
        engine.retrieve("test", patient_id='"; DROP TABLE')
        call_kwargs = mock_cm.search_all.call_args[1]
        filter_exprs = call_kwargs.get("filter_exprs", {})
        # Unsafe patient_id should not produce any filter
        assert not any("DROP" in v for v in filter_exprs.values())

    def test_retrieve_embedding_failure(self, engine, mock_cm):
        engine.embedder.encode.side_effect = RuntimeError("model crash")
        result = engine.retrieve("test")
        assert result.hits == []

    def test_retrieve_milvus_search_failure(self, engine, mock_cm):
        mock_cm.search_all.side_effect = Exception("Milvus down")
        result = engine.retrieve("test")
        assert result.hits == []


# ── Search (evidence only) ───────────────────────────────────────────────


class TestSearch:
    def test_search_returns_hits_list(self, engine, mock_cm):
        mock_cm.search_all.return_value = {
            "autoimmune_patient_labs": [
                {"id": "lab1", "score": 0.9, "text_chunk": "CRP elevated"},
            ],
        }
        hits = engine.search("CRP levels")
        assert isinstance(hits, list)
        assert len(hits) == 1
        assert isinstance(hits[0], SearchHit)


# ── find_related ─────────────────────────────────────────────────────────


class TestFindRelated:
    def test_find_related_groups_by_collection(self, engine, mock_cm):
        mock_cm.search_all.return_value = {
            "autoimmune_patient_labs": [
                {"id": "lab1", "score": 0.9, "text_chunk": "CRP"},
            ],
            "autoimmune_literature": [
                {"id": "lit1", "score": 0.85, "text_chunk": "Study on CRP"},
            ],
        }
        grouped = engine.find_related("CRP")
        assert "autoimmune_patient_labs" in grouped
        assert "autoimmune_literature" in grouped
        assert len(grouped["autoimmune_patient_labs"]) == 1


# ── Data classes ─────────────────────────────────────────────────────────


class TestDataClasses:
    def test_search_hit_defaults(self):
        hit = SearchHit(collection="test", id="1", score=0.5, text="hello")
        assert hit.metadata == {}
        assert hit.relevance == "low"

    def test_cross_collection_result_defaults(self):
        result = CrossCollectionResult(query="q")
        assert result.hits == []
        assert result.knowledge_context == ""
        assert result.total_collections_searched == 0
        assert result.search_time_ms == 0.0


# ── Regex patterns ───────────────────────────────────────────────────────


class TestRegexPatterns:
    def test_safe_filter_accepts_valid(self):
        assert _SAFE_FILTER_RE.match("Patient_001")
        assert _SAFE_FILTER_RE.match("abc-def.ghi")
        assert _SAFE_FILTER_RE.match("123")

    def test_safe_filter_rejects_injection(self):
        assert _SAFE_FILTER_RE.match('"; DROP') is None
        assert _SAFE_FILTER_RE.match("a'b") is None
        assert _SAFE_FILTER_RE.match("a;b") is None

    def test_comparative_regex(self):
        assert _COMPARATIVE_RE.search("compare these drugs")
        assert _COMPARATIVE_RE.search("drug A vs drug B")
        assert _COMPARATIVE_RE.search("difference between RA and SLE")
        assert _COMPARATIVE_RE.search("which is better for psoriasis")
        assert _COMPARATIVE_RE.search("head-to-head trial")
        assert _COMPARATIVE_RE.search("pros and cons")
        assert not _COMPARATIVE_RE.search("just a normal question")
