"""Tests for CAR-T Intelligence Agent RAG engine module.

Validates initialization, embedding prefixing, comparative detection,
citation formatting, merge-and-rank logic, knowledge context extraction,
retrieval, and COLLECTION_CONFIG structure.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.models import AgentQuery, CARTStage, CrossCollectionResult, SearchHit
from src.rag_engine import COLLECTION_CONFIG, STAGE_COLLECTION_BOOST, CARTRAGEngine


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def rag_engine(mock_embedder, mock_llm_client, mock_collection_manager):
    """Return a CARTRAGEngine with all mock components and knowledge enabled."""
    return CARTRAGEngine(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        llm_client=mock_llm_client,
        knowledge=True,
        query_expander=None,  # Disable expansion for isolated tests
    )


@pytest.fixture
def rag_engine_no_knowledge(mock_embedder, mock_llm_client, mock_collection_manager):
    """Return a CARTRAGEngine with knowledge disabled."""
    return CARTRAGEngine(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        llm_client=mock_llm_client,
        knowledge=None,
        query_expander=None,
    )


# ═══════════════════════════════════════════════════════════════════════
# INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════


class TestCARTRAGEngineInit:
    """Tests for CARTRAGEngine initialization."""

    def test_stores_all_components(
        self, mock_embedder, mock_llm_client, mock_collection_manager
    ):
        """__init__ stores all components as instance attributes."""
        engine = CARTRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge="knowledge_module",
            query_expander="expander_module",
        )
        assert engine.collections is mock_collection_manager
        assert engine.embedder is mock_embedder
        assert engine.llm is mock_llm_client
        assert engine.knowledge == "knowledge_module"
        assert engine.expander == "expander_module"

    def test_optional_components_default_to_none(
        self, mock_embedder, mock_llm_client, mock_collection_manager
    ):
        """Knowledge and query_expander default to None when not provided."""
        engine = CARTRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )
        assert engine.knowledge is None
        assert engine.expander is None


# ═══════════════════════════════════════════════════════════════════════
# EMBEDDING
# ═══════════════════════════════════════════════════════════════════════


class TestEmbedQuery:
    """Tests for _embed_query BGE prefix logic."""

    def test_adds_bge_prefix(self, rag_engine, mock_embedder):
        """_embed_query prepends the BGE instruction prefix to the query text."""
        rag_engine._embed_query("CD19 therapy")
        call_args = mock_embedder.embed_text.call_args[0][0]
        assert call_args.startswith("Represent this sentence for searching relevant passages: ")
        assert "CD19 therapy" in call_args

    def test_returns_embedding_vector(self, rag_engine):
        """_embed_query returns the embedder's output (384-dim vector)."""
        result = rag_engine._embed_query("test query")
        assert isinstance(result, list)
        assert len(result) == 384


# ═══════════════════════════════════════════════════════════════════════
# COMPARATIVE DETECTION
# ═══════════════════════════════════════════════════════════════════════


class TestIsComparative:
    """Tests for _is_comparative() pattern detection."""

    @pytest.mark.parametrize(
        "question",
        [
            "Compare CD19 and BCMA",
            "Kymriah vs Yescarta efficacy",
            "CD28 versus 4-1BB costimulation",
            "Comparing lentiviral and retroviral vectors",
        ],
    )
    def test_detects_comparative_queries(self, rag_engine, question):
        """_is_comparative returns True for comparative patterns."""
        assert rag_engine._is_comparative(question) is True

    @pytest.mark.parametrize(
        "question",
        [
            "What is the efficacy of CD19 CAR-T?",
            "How does CRS develop?",
            "What are the manufacturing steps?",
        ],
    )
    def test_non_comparative_queries(self, rag_engine, question):
        """_is_comparative returns False for non-comparative questions."""
        assert rag_engine._is_comparative(question) is False


# ═══════════════════════════════════════════════════════════════════════
# CITATION FORMATTING
# ═══════════════════════════════════════════════════════════════════════


class TestFormatCitation:
    """Tests for _format_citation() URL generation."""

    def test_pubmed_citation_creates_url(self):
        """PMIDs get PubMed URLs."""
        citation = CARTRAGEngine._format_citation("Literature", "12345678")
        assert "pubmed.ncbi.nlm.nih.gov/12345678/" in citation
        assert "PMID" in citation

    def test_nct_citation_creates_url(self):
        """NCT IDs get ClinicalTrials.gov URLs."""
        citation = CARTRAGEngine._format_citation("Trial", "NCT03958656")
        assert "clinicaltrials.gov/study/NCT03958656" in citation

    def test_generic_citation_uses_brackets(self):
        """Non-PubMed/non-NCT citations use [Collection:ID] format."""
        citation = CARTRAGEngine._format_citation("Construct", "construct-001")
        assert citation == "[Construct:construct-001]"

    def test_literature_non_digit_id(self):
        """Literature IDs that are not all digits use bracket format."""
        citation = CARTRAGEngine._format_citation("Literature", "patent-US123")
        assert citation == "[Literature:patent-US123]"

    def test_trial_non_nct_id(self):
        """Trial IDs not starting with NCT use bracket format."""
        citation = CARTRAGEngine._format_citation("Trial", "EUDRACT-12345")
        assert citation == "[Trial:EUDRACT-12345]"


# ═══════════════════════════════════════════════════════════════════════
# MERGE AND RANK
# ═══════════════════════════════════════════════════════════════════════


class TestMergeAndRank:
    """Tests for _merge_and_rank() deduplication and sorting."""

    def test_deduplicates_by_id(self, rag_engine):
        """Duplicate IDs are collapsed to a single hit (first seen)."""
        hits = [
            SearchHit(collection="Literature", id="1", score=0.9, text="A"),
            SearchHit(collection="Literature", id="1", score=0.7, text="A duplicate"),
            SearchHit(collection="Trial", id="2", score=0.8, text="B"),
        ]
        result = rag_engine._merge_and_rank(hits)
        assert len(result) == 2
        ids = [h.id for h in result]
        assert ids.count("1") == 1

    def test_sorts_by_score_descending(self, rag_engine):
        """Results are sorted highest score first."""
        hits = [
            SearchHit(collection="A", id="low", score=0.3, text="low"),
            SearchHit(collection="B", id="high", score=0.95, text="high"),
            SearchHit(collection="C", id="mid", score=0.6, text="mid"),
        ]
        result = rag_engine._merge_and_rank(hits)
        scores = [h.score for h in result]
        assert scores == sorted(scores, reverse=True)

    def test_caps_at_30_results(self, rag_engine):
        """_merge_and_rank never returns more than 30 results."""
        hits = [
            SearchHit(collection="Lit", id=str(i), score=0.5, text=f"hit {i}")
            for i in range(50)
        ]
        result = rag_engine._merge_and_rank(hits)
        assert len(result) <= 30

    def test_empty_input(self, rag_engine):
        """An empty list returns an empty list."""
        result = rag_engine._merge_and_rank([])
        assert result == []


# ═══════════════════════════════════════════════════════════════════════
# KNOWLEDGE CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetKnowledgeContext:
    """Tests for _get_knowledge_context()."""

    def test_cd19_crs_returns_context(self, rag_engine):
        """A query about CD19 and CRS returns non-empty knowledge context."""
        ctx = rag_engine._get_knowledge_context("What is CD19 CRS rate?")
        assert len(ctx) > 0
        assert "CD19" in ctx

    def test_no_knowledge_module_returns_empty(self, rag_engine_no_knowledge):
        """When knowledge is None, _get_knowledge_context returns empty string."""
        ctx = rag_engine_no_knowledge._get_knowledge_context("CD19 CRS")
        assert ctx == ""

    def test_manufacturing_query(self, rag_engine):
        """A manufacturing query returns manufacturing context."""
        ctx = rag_engine._get_knowledge_context("lentiviral transduction efficiency")
        assert len(ctx) > 0

    def test_biomarker_query(self, rag_engine):
        """A biomarker query returns biomarker context."""
        ctx = rag_engine._get_knowledge_context("Does ferritin predict CRS?")
        assert len(ctx) > 0

    def test_regulatory_query(self, rag_engine):
        """A regulatory query about Kymriah returns regulatory context."""
        ctx = rag_engine._get_knowledge_context("When was Kymriah approved?")
        assert len(ctx) > 0
        assert "Kymriah" in ctx or "tisagenlecleucel" in ctx


# ═══════════════════════════════════════════════════════════════════════
# RETRIEVE
# ═══════════════════════════════════════════════════════════════════════


class TestRetrieve:
    """Tests for the retrieve() method."""

    def test_returns_cross_collection_result(self, rag_engine):
        """retrieve() returns a CrossCollectionResult instance."""
        query = AgentQuery(question="CD19 CAR-T therapy efficacy")
        result = rag_engine.retrieve(query)
        assert isinstance(result, CrossCollectionResult)
        assert result.query == "CD19 CAR-T therapy efficacy"

    def test_search_time_is_positive(self, rag_engine):
        """retrieve() records a positive search_time_ms."""
        query = AgentQuery(question="test query")
        result = rag_engine.retrieve(query)
        assert result.search_time_ms >= 0

    def test_collections_searched_count(self, rag_engine):
        """retrieve() reports the correct number of collections searched."""
        query = AgentQuery(question="test")
        result = rag_engine.retrieve(query)
        assert result.total_collections_searched == 11

    def test_with_target_antigen_filter(self, rag_engine):
        """retrieve() accepts a target_antigen filter without error."""
        query = AgentQuery(question="CD19 trials", target_antigen="CD19")
        result = rag_engine.retrieve(query)
        assert isinstance(result, CrossCollectionResult)

    def test_with_collections_filter(self, rag_engine):
        """retrieve() filters to a subset of collections."""
        query = AgentQuery(question="test")
        result = rag_engine.retrieve(
            query, collections_filter=["cart_literature", "cart_trials"]
        )
        assert result.total_collections_searched == 2

    def test_with_year_range(self, rag_engine):
        """retrieve() accepts year_min and year_max without error."""
        query = AgentQuery(question="recent CD19 trials")
        result = rag_engine.retrieve(query, year_min=2020, year_max=2025)
        assert isinstance(result, CrossCollectionResult)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION_CONFIG
# ═══════════════════════════════════════════════════════════════════════


class TestCollectionConfig:
    """Tests for the COLLECTION_CONFIG dictionary."""

    def test_has_11_entries(self):
        """COLLECTION_CONFIG contains exactly 11 collection entries."""
        assert len(COLLECTION_CONFIG) == 11

    def test_expected_collection_names(self):
        """COLLECTION_CONFIG contains all expected collection names."""
        expected = {
            "cart_literature",
            "cart_trials",
            "cart_constructs",
            "cart_assays",
            "cart_manufacturing",
            "cart_safety",
            "cart_biomarkers",
            "cart_regulatory",
            "cart_sequences",
            "cart_realworld",
            "genomic_evidence",
        }
        assert set(COLLECTION_CONFIG.keys()) == expected

    def test_each_entry_has_required_keys(self):
        """Each config entry has weight, label, has_target_antigen, and year_field."""
        required_keys = {"weight", "label", "has_target_antigen", "year_field"}
        for name, config in COLLECTION_CONFIG.items():
            for key in required_keys:
                assert key in config, f"Missing '{key}' in COLLECTION_CONFIG['{name}']"

    def test_weights_are_positive(self):
        """All collection weights are positive floats."""
        for name, config in COLLECTION_CONFIG.items():
            assert config["weight"] > 0, f"Weight for '{name}' is not positive"

    def test_labels_are_non_empty(self):
        """All collection labels are non-empty strings."""
        for name, config in COLLECTION_CONFIG.items():
            assert isinstance(config["label"], str)
            assert len(config["label"]) > 0


# ═══════════════════════════════════════════════════════════════════════
# DYNAMIC WEIGHT BOOSTING
# ═══════════════════════════════════════════════════════════════════════


class TestComputeBoostedWeights:
    """Tests for _compute_boosted_weights() dynamic weight adjustment."""

    def test_clinical_stage_boosts_expected_collections(self, rag_engine):
        """CLINICAL stage boosts cart_trials, cart_safety, cart_realworld."""
        weights = rag_engine._compute_boosted_weights([CARTStage.CLINICAL])
        base = {name: cfg["weight"] for name, cfg in COLLECTION_CONFIG.items()}
        # Boosted collections should have higher relative weights
        for coll in STAGE_COLLECTION_BOOST[CARTStage.CLINICAL]:
            base_ratio = base[coll] / sum(base.values())
            assert weights[coll] > base_ratio

    def test_car_design_stage_boosts_constructs_and_sequences(self, rag_engine):
        """CAR_DESIGN stage boosts cart_constructs and cart_sequences."""
        weights = rag_engine._compute_boosted_weights([CARTStage.CAR_DESIGN])
        base = {name: cfg["weight"] for name, cfg in COLLECTION_CONFIG.items()}
        for coll in STAGE_COLLECTION_BOOST[CARTStage.CAR_DESIGN]:
            base_ratio = base[coll] / sum(base.values())
            assert weights[coll] > base_ratio

    def test_weights_sum_to_one(self, rag_engine):
        """Boosted weights always sum to approximately 1.0."""
        weights = rag_engine._compute_boosted_weights([CARTStage.CLINICAL])
        assert abs(sum(weights.values()) - 1.0) < 1e-9

    def test_multiple_stages_sum_to_one(self, rag_engine):
        """Weights still sum to ~1.0 when multiple stages are boosted."""
        weights = rag_engine._compute_boosted_weights(
            [CARTStage.CLINICAL, CARTStage.TESTING]
        )
        assert abs(sum(weights.values()) - 1.0) < 1e-9

    def test_all_stages_sum_to_one(self, rag_engine):
        """Weights sum to ~1.0 even when all stages are boosted."""
        all_stages = list(CARTStage)
        weights = rag_engine._compute_boosted_weights(all_stages)
        assert abs(sum(weights.values()) - 1.0) < 1e-9

    def test_empty_stages_returns_normalized_base(self, rag_engine):
        """Empty stages list returns base weights normalized to 1.0."""
        weights = rag_engine._compute_boosted_weights([])
        assert abs(sum(weights.values()) - 1.0) < 1e-9
        # Should be proportional to base weights
        base = {name: cfg["weight"] for name, cfg in COLLECTION_CONFIG.items()}
        total = sum(base.values())
        for name in weights:
            expected = base[name] / total
            assert abs(weights[name] - expected) < 1e-9

    def test_all_collections_present_in_output(self, rag_engine):
        """Boosted weights contain entries for all collections in COLLECTION_CONFIG."""
        weights = rag_engine._compute_boosted_weights([CARTStage.CLINICAL])
        assert set(weights.keys()) == set(COLLECTION_CONFIG.keys())

    def test_boost_factor_is_1_5x(self, rag_engine):
        """Boosted collections get exactly 1.5x their base weight before normalization."""
        stages = [CARTStage.VECTOR_ENG]  # boosts only cart_manufacturing
        weights = rag_engine._compute_boosted_weights(stages)
        # Verify by computing manually
        base = {name: cfg["weight"] for name, cfg in COLLECTION_CONFIG.items()}
        base["cart_manufacturing"] *= 1.5
        total = sum(base.values())
        expected = base["cart_manufacturing"] / total
        assert abs(weights["cart_manufacturing"] - expected) < 1e-9

    def test_retrieve_accepts_stages_parameter(self, rag_engine):
        """retrieve() accepts an optional stages parameter without error."""
        query = AgentQuery(question="CD19 clinical trials")
        result = rag_engine.retrieve(query, stages=[CARTStage.CLINICAL])
        assert isinstance(result, CrossCollectionResult)
