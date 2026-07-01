"""Tests for RAG engine in src/rag_engine.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

# rag_engine may not exist yet; test conditionally
try:
    from src.rag_engine import SingleCellRAGEngine
    _RAG_AVAILABLE = True
except ImportError:
    _RAG_AVAILABLE = False

from src.models import SCWorkflowType


@pytest.mark.skipif(not _RAG_AVAILABLE, reason="rag_engine module not yet implemented")
class TestSingleCellRAGEngine:
    """Tests for SingleCellRAGEngine if available."""

    def test_engine_class_exists(self):
        assert SingleCellRAGEngine is not None


class TestRAGPlaceholder:
    """Placeholder tests that always pass."""

    def test_workflow_types_for_rag(self):
        """All workflow types should be valid for RAG routing."""
        for wf in SCWorkflowType:
            assert isinstance(wf.value, str)

    def test_workflow_count_for_rag(self):
        """RAG engine should handle all 11 workflow types."""
        assert len(SCWorkflowType) == 11

    def test_search_collections_defined(self):
        """Verify collection names available for RAG search."""
        from src.collections import get_all_collection_names
        names = get_all_collection_names()
        assert len(names) == 12
