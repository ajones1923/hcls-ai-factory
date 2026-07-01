"""
Tests for Target Hypothesis module.
"""
import pytest
import json
from pathlib import Path
from datetime import datetime

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.target_hypothesis import (
    TargetHypothesis,
    TargetHypothesisManager,
    create_hypothesis_from_chat,
)


class TestTargetHypothesis:
    """Tests for TargetHypothesis dataclass."""

    def test_hypothesis_creation(self):
        """Test creating a target hypothesis."""
        hyp = TargetHypothesis(
            gene="VCP",
            protein="Valosin-containing protein",
            rationale="VCP mutations cause IBMPFD with FTD",
        )

        assert hyp.gene == "VCP"
        assert hyp.protein == "Valosin-containing protein"
        assert hyp.id is not None
        assert hyp.id.startswith("TH-")
        assert "VCP" in hyp.id

    def test_hypothesis_auto_id(self):
        """Test automatic ID generation."""
        hyp = TargetHypothesis(gene="EGFR")
        assert hyp.id is not None
        assert hyp.created_at is not None
        assert hyp.updated_at is not None

    def test_hypothesis_with_id(self):
        """Test hypothesis with provided ID."""
        hyp = TargetHypothesis(gene="BRAF", id="TH-CUSTOM-001")
        assert hyp.id == "TH-CUSTOM-001"

    def test_hypothesis_variant_count(self):
        """Test variant count is calculated."""
        variants = [
            {"chrom": "chr9", "pos": 123},
            {"chrom": "chr9", "pos": 456},
        ]
        hyp = TargetHypothesis(gene="VCP", variants=variants)
        assert hyp.variant_count == 2

    def test_add_variant(self):
        """Test adding variants."""
        hyp = TargetHypothesis(gene="VCP")
        assert hyp.variant_count == 0

        hyp.add_variant({"chrom": "chr9", "pos": 123})
        assert hyp.variant_count == 1

        hyp.add_variant({"chrom": "chr9", "pos": 456})
        assert hyp.variant_count == 2

    def test_to_dict(self):
        """Test conversion to dictionary."""
        hyp = TargetHypothesis(
            gene="VCP",
            protein="p97",
            therapeutic_area="Neurodegeneration",
            confidence="high",
            priority=5,
        )
        d = hyp.to_dict()

        assert d["gene"] == "VCP"
        assert d["protein"] == "p97"
        assert d["therapeutic_area"] == "Neurodegeneration"
        assert d["confidence"] == "high"
        assert d["priority"] == 5

    def test_from_dict(self):
        """Test creation from dictionary."""
        data = {
            "gene": "EGFR",
            "protein": "Epidermal growth factor receptor",
            "rationale": "Key oncology target",
            "confidence": "high",
            "priority": 5,
            "status": "validated",
            "id": "TH-TEST-001",
            "created_at": "2024-01-01T00:00:00",
            "updated_at": "2024-01-01T00:00:00",
            "variants": [],
            "variant_count": 0,
            "pdb_ids": ["7SYE"],
        }
        hyp = TargetHypothesis.from_dict(data)

        assert hyp.gene == "EGFR"
        assert hyp.confidence == "high"
        assert hyp.status == "validated"
        assert "7SYE" in hyp.pdb_ids

    def test_update_rationale(self):
        """Test updating rationale."""
        hyp = TargetHypothesis(gene="VCP", rationale="Initial rationale")
        old_updated = hyp.updated_at

        hyp.update_rationale("Updated rationale with new evidence")

        assert hyp.rationale == "Updated rationale with new evidence"
        assert hyp.updated_at != old_updated


class TestTargetHypothesisManager:
    """Tests for TargetHypothesisManager class."""

    def test_manager_initialization(self, temp_hypotheses_dir):
        """Test manager initialization."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        assert len(manager.hypotheses) == 0

    def test_add_hypothesis(self, temp_hypotheses_dir):
        """Test adding a hypothesis."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        hyp = TargetHypothesis(gene="VCP", rationale="Test")

        hyp_id = manager.add(hyp)

        assert hyp_id == hyp.id
        assert hyp_id in manager.hypotheses

    def test_get_hypothesis(self, temp_hypotheses_dir):
        """Test retrieving a hypothesis by ID."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        hyp = TargetHypothesis(gene="VCP", rationale="Test")
        hyp_id = manager.add(hyp)

        retrieved = manager.get(hyp_id)

        assert retrieved is not None
        assert retrieved.gene == "VCP"

    def test_get_by_gene(self, temp_hypotheses_dir):
        """Test retrieving a hypothesis by gene name."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager.add(TargetHypothesis(gene="VCP", rationale="Test VCP"))
        manager.add(TargetHypothesis(gene="EGFR", rationale="Test EGFR"))

        vcp = manager.get_by_gene("VCP")
        egfr = manager.get_by_gene("egfr")  # Case insensitive

        assert vcp is not None
        assert vcp.gene == "VCP"
        assert egfr is not None
        assert egfr.gene == "EGFR"

    def test_list_all(self, temp_hypotheses_dir):
        """Test listing all hypotheses."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager.add(TargetHypothesis(gene="VCP"))
        manager.add(TargetHypothesis(gene="EGFR"))
        manager.add(TargetHypothesis(gene="BRAF"))

        all_hyps = manager.list_all()

        assert len(all_hyps) == 3

    def test_list_by_status(self, temp_hypotheses_dir):
        """Test filtering by status."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager.add(TargetHypothesis(gene="VCP", status="validated"))
        manager.add(TargetHypothesis(gene="EGFR", status="validated"))
        manager.add(TargetHypothesis(gene="BRAF", status="hypothesis"))

        validated = manager.list_by_status("validated")
        hypotheses = manager.list_by_status("hypothesis")

        assert len(validated) == 2
        assert len(hypotheses) == 1

    def test_list_by_priority(self, temp_hypotheses_dir):
        """Test filtering by priority."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager.add(TargetHypothesis(gene="VCP", priority=5))
        manager.add(TargetHypothesis(gene="EGFR", priority=4))
        manager.add(TargetHypothesis(gene="BRAF", priority=2))

        high_priority = manager.list_by_priority(min_priority=4)

        assert len(high_priority) == 2
        # Should be sorted by priority descending
        assert high_priority[0].priority >= high_priority[1].priority

    def test_update_hypothesis(self, temp_hypotheses_dir):
        """Test updating a hypothesis."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        hyp = TargetHypothesis(gene="VCP", status="hypothesis", priority=3)
        hyp_id = manager.add(hyp)

        success = manager.update(hyp_id, {
            "status": "validated",
            "priority": 5,
            "rationale": "Updated rationale",
        })

        assert success
        updated = manager.get(hyp_id)
        assert updated.status == "validated"
        assert updated.priority == 5
        assert updated.rationale == "Updated rationale"

    def test_delete_hypothesis(self, temp_hypotheses_dir):
        """Test deleting a hypothesis."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        hyp = TargetHypothesis(gene="VCP")
        hyp_id = manager.add(hyp)

        success = manager.delete(hyp_id)

        assert success
        assert manager.get(hyp_id) is None

    def test_persistence(self, temp_hypotheses_dir):
        """Test hypotheses persist across manager instances."""
        # Create and save
        manager1 = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager1.add(TargetHypothesis(gene="VCP", rationale="Persistent"))

        # Create new manager pointing to same storage
        manager2 = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)

        assert len(manager2.hypotheses) == 1
        vcp = manager2.get_by_gene("VCP")
        assert vcp is not None
        assert vcp.rationale == "Persistent"

    def test_export_for_phase5(self, temp_hypotheses_dir):
        """Test export for Phase 5 (drug discovery)."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)

        # Add some hypotheses with different statuses
        manager.add(TargetHypothesis(
            gene="VCP",
            protein="p97",
            status="validated",
            priority=5,
            pdb_ids=["5FTK"],
        ))
        manager.add(TargetHypothesis(
            gene="EGFR",
            status="selected",
            priority=4,
        ))
        manager.add(TargetHypothesis(
            gene="BRAF",
            status="hypothesis",
            priority=2,  # Low priority, should not be exported
        ))

        output_file = manager.export_for_phase5()

        assert output_file.exists()

        with open(output_file) as f:
            data = json.load(f)

        # Only validated/selected/high-priority should be exported
        assert data["count"] == 2
        genes = [t["gene"] for t in data["targets"]]
        assert "VCP" in genes
        assert "EGFR" in genes
        assert "BRAF" not in genes

    def test_get_summary(self, temp_hypotheses_dir):
        """Test summary statistics."""
        manager = TargetHypothesisManager(storage_dir=temp_hypotheses_dir)
        manager.add(TargetHypothesis(gene="VCP", status="validated", confidence="high", priority=5))
        manager.add(TargetHypothesis(gene="EGFR", status="validated", confidence="medium", priority=4))
        manager.add(TargetHypothesis(gene="BRAF", status="hypothesis", confidence="low", priority=2))

        summary = manager.get_summary()

        assert summary["total"] == 3
        assert summary["by_status"]["validated"] == 2
        assert summary["by_status"]["hypothesis"] == 1
        assert summary["high_priority"] == 2  # Priority >= 4


class TestCreateHypothesisFromChat:
    """Tests for the helper function."""

    def test_create_from_chat(self):
        """Test creating hypothesis from chat interaction."""
        hyp = create_hypothesis_from_chat(
            gene="VCP",
            rationale="User asked about dementia targets, VCP is associated with FTD",
            source_query="What genes are associated with dementia?",
            therapeutic_area="Neurodegeneration",
        )

        assert hyp.gene == "VCP"
        assert hyp.rationale == "User asked about dementia targets, VCP is associated with FTD"
        assert hyp.source_query == "What genes are associated with dementia?"
        assert hyp.therapeutic_area == "Neurodegeneration"

    def test_create_with_variants(self):
        """Test creating hypothesis with variant evidence."""
        variants = [
            {"id": "chr9_123_A_G", "gene": "VCP", "impact": "HIGH"},
            {"id": "chr9_456_C_T", "gene": "VCP", "impact": "MODERATE"},
        ]
        hyp = create_hypothesis_from_chat(
            gene="VCP",
            rationale="Multiple pathogenic variants found",
            variants=variants,
        )

        assert hyp.variant_count == 2
        assert len(hyp.variants) == 2
