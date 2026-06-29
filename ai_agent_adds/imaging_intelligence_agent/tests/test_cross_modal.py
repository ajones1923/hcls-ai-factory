"""Tests for cross-modal imaging -> genomics trigger.

Validates CrossModalTrigger.evaluate() for lung nodule, CXR, and
brain lesion workflows.  Confirms that:
  - Lung-RADS 4A, 4B, 4X trigger genomic queries
  - Lung-RADS 1-3 do NOT trigger
  - CXR critical consolidation triggers
  - CXR normal/routine does NOT trigger
  - Brain highly_active triggers
  - Brain stable/active does NOT trigger
  - Disabled trigger always returns None
  - Genomic hits are de-duplicated
  - Unknown workflows return None

Author: Adam Jones
Date: February 2026
"""

import pytest
from unittest.mock import MagicMock

from src.cross_modal import (
    CrossModalTrigger,
    LUNG_CANCER_QUERIES,
    INFECTION_GENOMICS_QUERIES,
    NEURO_GENOMICS_QUERIES,
)
from src.models import (
    CrossModalResult,
    FindingSeverity,
    WorkflowResult,
    WorkflowStatus,
)


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def mock_genomic_hits():
    """Return sample genomic evidence search hits."""
    return [
        {
            "id": "clinvar-EGFR-L858R",
            "score": 0.89,
            "collection": "genomic_evidence",
            "text_chunk": (
                "EGFR L858R is the most common activating mutation in "
                "non-small cell lung cancer, present in ~40% of EGFR-mutant "
                "NSCLC. Responds to erlotinib, gefitinib, and osimertinib."
            ),
        },
        {
            "id": "clinvar-ALK-EML4",
            "score": 0.82,
            "collection": "genomic_evidence",
            "text_chunk": (
                "EML4-ALK fusion is found in ~5% of NSCLC patients. "
                "First-line treatment with alectinib or lorlatinib."
            ),
        },
        {
            "id": "clinvar-KRAS-G12C",
            "score": 0.78,
            "collection": "genomic_evidence",
            "text_chunk": (
                "KRAS G12C mutation occurs in ~13% of NSCLC. "
                "Sotorasib (Lumakras) is a targeted covalent inhibitor."
            ),
        },
    ]


@pytest.fixture
def cross_modal_trigger(mock_collection_manager, mock_embedder):
    """Create a CrossModalTrigger with mock dependencies."""
    return CrossModalTrigger(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        enabled=True,
    )


@pytest.fixture
def disabled_trigger(mock_collection_manager, mock_embedder):
    """Create a disabled CrossModalTrigger."""
    return CrossModalTrigger(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        enabled=False,
    )


def _lung_nodule_result(classification: str, severity: FindingSeverity) -> WorkflowResult:
    """Helper to build a ct_chest_lung_nodule WorkflowResult."""
    return WorkflowResult(
        workflow_name="ct_chest_lung_nodule",
        status=WorkflowStatus.COMPLETED,
        classification=classification,
        severity=severity,
        findings=[{
            "category": "nodule",
            "description": f"Solid nodule, {classification}",
            "severity": severity.value,
        }],
        measurements={"nodule_count": 1.0},
        is_mock=True,
    )


def _cxr_result(
    severity: FindingSeverity, has_consolidation: bool = False
) -> WorkflowResult:
    """Helper to build a cxr_rapid_findings WorkflowResult."""
    findings = []
    if has_consolidation:
        findings.append({
            "category": "consolidation",
            "description": "Pulmonary consolidation",
            "severity": severity.value,
            "class_name": "consolidation",
            "confidence": 0.87,
        })
    else:
        findings.append({
            "category": "normal",
            "description": "No acute findings",
            "severity": "normal",
        })
    classification = "positive: consolidation" if has_consolidation else "negative"
    return WorkflowResult(
        workflow_name="cxr_rapid_findings",
        status=WorkflowStatus.COMPLETED,
        classification=classification,
        severity=severity,
        findings=findings,
        is_mock=True,
    )


def _brain_result(classification: str, severity: FindingSeverity) -> WorkflowResult:
    """Helper to build a mri_brain_ms_lesion WorkflowResult."""
    return WorkflowResult(
        workflow_name="mri_brain_ms_lesion",
        status=WorkflowStatus.COMPLETED,
        classification=classification,
        severity=severity,
        findings=[{
            "category": "lesion",
            "description": f"MS lesion activity: {classification}",
            "severity": severity.value,
        }],
        measurements={"total_lesion_count": 12.0, "new_lesion_count": 5.0},
        is_mock=True,
    )


# ===================================================================
# LUNG NODULE TRIGGERS
# ===================================================================


class TestLungNoduleTrigger:
    """Tests for Lung-RADS cross-modal triggers."""

    def test_lung_rads_4b_triggers(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert isinstance(cm, CrossModalResult)
        assert "4B" in cm.trigger_reason
        assert cm.query_count == len(LUNG_CANCER_QUERIES)

    def test_lung_rads_4a_triggers(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 4A", FindingSeverity.URGENT)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert "4A" in cm.trigger_reason

    def test_lung_rads_4x_triggers(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 4X", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert "4X" in cm.trigger_reason

    def test_lung_rads_1_does_not_trigger(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 1", FindingSeverity.NORMAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_lung_rads_2_does_not_trigger(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 2", FindingSeverity.ROUTINE)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_lung_rads_3_does_not_trigger(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 3", FindingSeverity.SIGNIFICANT)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_genomic_hits_returned(self, cross_modal_trigger, mock_genomic_hits):
        """Verify genomic hits from the mock collection manager are captured."""
        cross_modal_trigger.collection_manager.search.return_value = mock_genomic_hits
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        # Each query returns 3 unique hits; de-dup means max = 3
        assert cm.genomic_hit_count == 3
        assert len(cm.genomic_context) == 3
        assert "EGFR" in cm.genomic_context[0]

    def test_enrichment_summary_populated(self, cross_modal_trigger):
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert len(cm.enrichment_summary) > 0
        assert "Cross-modal" in cm.enrichment_summary


# ===================================================================
# CXR TRIGGERS
# ===================================================================


class TestCXRTrigger:
    """Tests for CXR cross-modal triggers."""

    def test_urgent_consolidation_triggers(self, cross_modal_trigger):
        result = _cxr_result(FindingSeverity.URGENT, has_consolidation=True)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert "consolidation" in cm.trigger_reason.lower()
        assert cm.query_count == len(INFECTION_GENOMICS_QUERIES)

    def test_critical_consolidation_triggers(self, cross_modal_trigger):
        result = _cxr_result(FindingSeverity.CRITICAL, has_consolidation=True)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None

    def test_routine_does_not_trigger(self, cross_modal_trigger):
        result = _cxr_result(FindingSeverity.ROUTINE, has_consolidation=True)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_normal_does_not_trigger(self, cross_modal_trigger):
        result = _cxr_result(FindingSeverity.NORMAL, has_consolidation=False)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_urgent_without_consolidation_does_not_trigger(self, cross_modal_trigger):
        """Urgent CXR without consolidation should not trigger infection genomics."""
        result = _cxr_result(FindingSeverity.URGENT, has_consolidation=False)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None


# ===================================================================
# BRAIN LESION TRIGGERS
# ===================================================================


class TestBrainLesionTrigger:
    """Tests for brain lesion cross-modal triggers."""

    def test_highly_active_triggers(self, cross_modal_trigger):
        result = _brain_result("ms_highly_active", FindingSeverity.URGENT)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert "highly active" in cm.trigger_reason.lower()
        assert cm.query_count == len(NEURO_GENOMICS_QUERIES)

    def test_active_does_not_trigger(self, cross_modal_trigger):
        result = _brain_result("ms_active", FindingSeverity.SIGNIFICANT)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None

    def test_stable_does_not_trigger(self, cross_modal_trigger):
        result = _brain_result("ms_stable", FindingSeverity.ROUTINE)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None


# ===================================================================
# DISABLED TRIGGER
# ===================================================================


class TestDisabledTrigger:
    """Tests that disabled trigger always returns None."""

    def test_disabled_lung_rads_4b(self, disabled_trigger):
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = disabled_trigger.evaluate(result)
        assert cm is None

    def test_disabled_cxr_urgent(self, disabled_trigger):
        result = _cxr_result(FindingSeverity.URGENT, has_consolidation=True)
        cm = disabled_trigger.evaluate(result)
        assert cm is None

    def test_disabled_brain_highly_active(self, disabled_trigger):
        result = _brain_result("ms_highly_active", FindingSeverity.URGENT)
        cm = disabled_trigger.evaluate(result)
        assert cm is None


# ===================================================================
# UNKNOWN WORKFLOW
# ===================================================================


class TestUnknownWorkflow:
    """Tests that unknown workflows do not trigger."""

    def test_unknown_workflow_returns_none(self, cross_modal_trigger):
        result = WorkflowResult(
            workflow_name="pet_ct_oncology",
            status=WorkflowStatus.COMPLETED,
            severity=FindingSeverity.CRITICAL,
            classification="SUV_max_15.2",
        )
        cm = cross_modal_trigger.evaluate(result)
        assert cm is None


# ===================================================================
# DE-DUPLICATION
# ===================================================================


class TestDeduplication:
    """Tests that genomic hits are de-duplicated across queries."""

    def test_duplicate_hits_are_removed(self, cross_modal_trigger):
        """If the same hit ID appears in multiple query results, count once."""
        duplicate_hits = [
            {"id": "clinvar-EGFR-L858R", "score": 0.89, "text_chunk": "EGFR L858R info"},
            {"id": "clinvar-EGFR-L858R", "score": 0.85, "text_chunk": "EGFR L858R info"},
            {"id": "clinvar-ALK-EML4", "score": 0.82, "text_chunk": "ALK info"},
        ]
        cross_modal_trigger.collection_manager.search.return_value = duplicate_hits
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        # 3 queries each return 3 hits, but IDs are only 2 unique
        assert cm.genomic_hit_count == 2


# ===================================================================
# SEARCH FAILURE RESILIENCE
# ===================================================================


class TestSearchFailureResilience:
    """Tests that the trigger handles search failures gracefully."""

    def test_search_exception_returns_empty_result(self, cross_modal_trigger):
        """If Milvus search raises, trigger should still return a result with 0 hits."""
        cross_modal_trigger.collection_manager.search.side_effect = Exception(
            "Milvus connection refused"
        )
        result = _lung_nodule_result("Lung-RADS 4B", FindingSeverity.CRITICAL)
        cm = cross_modal_trigger.evaluate(result)
        assert cm is not None
        assert cm.genomic_hit_count == 0
        assert cm.query_count == len(LUNG_CANCER_QUERIES)
        assert "threshold" in cm.enrichment_summary.lower()
