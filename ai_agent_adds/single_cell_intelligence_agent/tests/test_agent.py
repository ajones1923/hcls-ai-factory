"""Tests for the Single-Cell Intelligence Agent in src/agent.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

# agent module may have heavy dependencies; test conditionally
try:
    from src.agent import SingleCellAgent
    _AGENT_AVAILABLE = True
except ImportError:
    _AGENT_AVAILABLE = False

from src.models import SCWorkflowType, SeverityLevel


@pytest.mark.skipif(not _AGENT_AVAILABLE, reason="agent module not yet fully available")
class TestSingleCellAgent:
    """Tests for SingleCellAgent if available."""

    def test_agent_class_exists(self):
        assert SingleCellAgent is not None


class TestAgentPlaceholder:
    """Placeholder tests validating agent-related models."""

    def test_workflow_types_complete(self):
        """Agent should support all 11 workflow types."""
        expected = {
            "cell_type_annotation", "tme_profiling", "drug_response",
            "subclonal_architecture", "spatial_niche", "trajectory_analysis",
            "ligand_receptor", "biomarker_discovery", "cart_target_validation",
            "treatment_monitoring", "general",
        }
        actual = {wf.value for wf in SCWorkflowType}
        assert expected == actual

    def test_severity_levels_complete(self):
        """Agent should have severity levels for clinical grading."""
        assert len(SeverityLevel) == 5

    def test_critical_severity_exists(self):
        assert SeverityLevel.CRITICAL.value == "critical"

    def test_informational_severity_exists(self):
        assert SeverityLevel.INFORMATIONAL.value == "informational"
