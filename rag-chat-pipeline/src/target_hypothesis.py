"""
Target Hypothesis Manager for RAG Chat Pipeline.
Manages drug target hypotheses generated from variant analysis.
"""
import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any
from loguru import logger


@dataclass
class TargetHypothesis:
    """Represents a potential drug target hypothesis."""

    # Target identification
    gene: str
    protein: Optional[str] = None
    uniprot_id: Optional[str] = None

    # Variant evidence
    variants: List[Dict[str, Any]] = field(default_factory=list)
    variant_count: int = 0

    # Hypothesis details
    rationale: str = ""
    therapeutic_area: Optional[str] = None
    mechanism: Optional[str] = None

    # Confidence and priority
    confidence: str = "medium"  # low, medium, high
    priority: int = 0  # 1-5 scale

    # Metadata
    id: Optional[str] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    source_query: Optional[str] = None
    notes: str = ""

    # Status for Phase 5 handoff
    status: str = "hypothesis"  # hypothesis, validated, selected, rejected
    pdb_ids: List[str] = field(default_factory=list)  # For Cryo-EM phase

    def __post_init__(self):
        if not self.id:
            self.id = f"TH-{datetime.now().strftime('%Y%m%d%H%M%S')}-{self.gene}"
        if not self.created_at:
            self.created_at = datetime.now().isoformat()
        self.updated_at = datetime.now().isoformat()
        self.variant_count = len(self.variants)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "TargetHypothesis":
        """Create from dictionary."""
        return cls(**data)

    def add_variant(self, variant: Dict[str, Any]):
        """Add a variant to the evidence."""
        self.variants.append(variant)
        self.variant_count = len(self.variants)
        self.updated_at = datetime.now().isoformat()

    def update_rationale(self, rationale: str):
        """Update the hypothesis rationale."""
        self.rationale = rationale
        self.updated_at = datetime.now().isoformat()


class TargetHypothesisManager:
    """Manages collection of target hypotheses."""

    def __init__(self, storage_dir: Optional[Path] = None):
        self.storage_dir = storage_dir or Path("data/targets")
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.hypotheses: Dict[str, TargetHypothesis] = {}
        self._load_existing()

    def _load_existing(self):
        """Load existing hypotheses from storage."""
        hypotheses_file = self.storage_dir / "hypotheses.json"
        if hypotheses_file.exists():
            try:
                with open(hypotheses_file, 'r') as f:
                    data = json.load(f)
                    for item in data.get('hypotheses', []):
                        hyp = TargetHypothesis.from_dict(item)
                        self.hypotheses[hyp.id] = hyp
                logger.info(f"Loaded {len(self.hypotheses)} existing hypotheses")
            except Exception as e:
                logger.error(f"Error loading hypotheses: {e}")

    def save(self):
        """Save all hypotheses to storage."""
        hypotheses_file = self.storage_dir / "hypotheses.json"
        data = {
            'version': '1.0',
            'updated_at': datetime.now().isoformat(),
            'count': len(self.hypotheses),
            'hypotheses': [h.to_dict() for h in self.hypotheses.values()]
        }
        with open(hypotheses_file, 'w') as f:
            json.dump(data, f, indent=2)
        logger.info(f"Saved {len(self.hypotheses)} hypotheses")

    def add(self, hypothesis: TargetHypothesis) -> str:
        """Add a new hypothesis."""
        self.hypotheses[hypothesis.id] = hypothesis
        self.save()
        return hypothesis.id

    def get(self, hypothesis_id: str) -> Optional[TargetHypothesis]:
        """Get a hypothesis by ID."""
        return self.hypotheses.get(hypothesis_id)

    def get_by_gene(self, gene: str) -> Optional[TargetHypothesis]:
        """Get hypothesis by gene name."""
        for hyp in self.hypotheses.values():
            if hyp.gene.upper() == gene.upper():
                return hyp
        return None

    def list_all(self) -> List[TargetHypothesis]:
        """Get all hypotheses."""
        return list(self.hypotheses.values())

    def list_by_status(self, status: str) -> List[TargetHypothesis]:
        """Get hypotheses by status."""
        return [h for h in self.hypotheses.values() if h.status == status]

    def list_by_priority(self, min_priority: int = 1) -> List[TargetHypothesis]:
        """Get hypotheses with minimum priority."""
        return sorted(
            [h for h in self.hypotheses.values() if h.priority >= min_priority],
            key=lambda x: x.priority,
            reverse=True
        )

    def update(self, hypothesis_id: str, updates: Dict[str, Any]) -> bool:
        """Update a hypothesis."""
        if hypothesis_id not in self.hypotheses:
            return False

        hyp = self.hypotheses[hypothesis_id]
        for key, value in updates.items():
            if hasattr(hyp, key):
                setattr(hyp, key, value)
        hyp.updated_at = datetime.now().isoformat()
        self.save()
        return True

    def delete(self, hypothesis_id: str) -> bool:
        """Delete a hypothesis."""
        if hypothesis_id in self.hypotheses:
            del self.hypotheses[hypothesis_id]
            self.save()
            return True
        return False

    def export_for_phase5(self, output_file: Optional[Path] = None) -> Path:
        """Export selected targets for Phase 5 (Cryo-EM)."""
        output_file = output_file or self.storage_dir / "targets_for_phase5.json"

        # Filter to selected/validated targets
        selected = [
            h for h in self.hypotheses.values()
            if h.status in ('validated', 'selected') or h.priority >= 4
        ]

        export_data = {
            'version': '1.0',
            'exported_at': datetime.now().isoformat(),
            'source': 'rag-chat-pipeline',
            'next_phase': 'target-selection-pipeline',
            'count': len(selected),
            'targets': []
        }

        for hyp in selected:
            export_data['targets'].append({
                'id': hyp.id,
                'gene': hyp.gene,
                'protein': hyp.protein,
                'uniprot_id': hyp.uniprot_id,
                'rationale': hyp.rationale,
                'therapeutic_area': hyp.therapeutic_area,
                'mechanism': hyp.mechanism,
                'confidence': hyp.confidence,
                'priority': hyp.priority,
                'variant_count': hyp.variant_count,
                'pdb_ids': hyp.pdb_ids,
                'status': hyp.status,
            })

        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)

        logger.info(f"Exported {len(selected)} targets to {output_file}")
        return output_file

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        total = len(self.hypotheses)
        by_status = {}
        by_confidence = {}

        for hyp in self.hypotheses.values():
            by_status[hyp.status] = by_status.get(hyp.status, 0) + 1
            by_confidence[hyp.confidence] = by_confidence.get(hyp.confidence, 0) + 1

        return {
            'total': total,
            'by_status': by_status,
            'by_confidence': by_confidence,
            'high_priority': len([h for h in self.hypotheses.values() if h.priority >= 4]),
        }


def create_hypothesis_from_chat(
    gene: str,
    rationale: str,
    variants: List[Dict[str, Any]] = None,
    source_query: str = None,
    **kwargs
) -> TargetHypothesis:
    """
    Helper to create a hypothesis from chat interaction.
    """
    return TargetHypothesis(
        gene=gene,
        rationale=rationale,
        variants=variants or [],
        source_query=source_query,
        **kwargs
    )
