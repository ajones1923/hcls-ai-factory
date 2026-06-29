"""Tests for the Neurology Intelligence Agent in src/agent.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

# agent module may not exist yet; test conditionally
try:
    from src.agent import NeuroAgent
    _AGENT_AVAILABLE = True
except ImportError:
    _AGENT_AVAILABLE = False

from src.models import NeuroWorkflowType, SeverityLevel


@pytest.mark.skipif(not _AGENT_AVAILABLE, reason="agent module not yet implemented")
class TestNeuroAgent:
    """Tests for NeuroAgent if available."""

    def test_agent_instantiation(self):
        agent = NeuroAgent()
        assert agent is not None


class TestAgentPlaceholder:
    """Placeholder tests validating agent-related models."""

    def test_workflow_types_complete(self):
        """Agent should support all 9 workflow types."""
        expected = {
            "acute_stroke", "dementia_evaluation", "epilepsy_focus",
            "brain_tumor", "ms_monitoring", "parkinsons_assessment",
            "headache_classification", "neuromuscular_evaluation", "general",
        }
        actual = {wf.value for wf in NeuroWorkflowType}
        assert expected == actual

    def test_severity_levels_complete(self):
        """Agent should have severity levels for clinical grading."""
        assert len(SeverityLevel) == 5

    def test_critical_severity_exists(self):
        assert SeverityLevel.CRITICAL.value == "critical"

    def test_informational_severity_exists(self):
        assert SeverityLevel.INFORMATIONAL.value == "informational"
