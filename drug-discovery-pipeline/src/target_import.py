"""
Target Import Module - Imports target hypotheses from RAG Chat pipeline.
"""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
from datetime import datetime


@dataclass
class ImportedTarget:
    """Target hypothesis imported from RAG Chat."""
    gene: str
    protein: Optional[str]
    rationale: str
    confidence: str  # low, medium, high
    priority: int  # 1-5
    therapeutic_area: Optional[str]
    mechanism: Optional[str]
    pdb_ids: List[str]
    variant_count: int
    variants: List[Dict[str, Any]]
    source_query: Optional[str]
    status: str
    imported_at: str = None

    def __post_init__(self):
        if self.imported_at is None:
            self.imported_at = datetime.now().isoformat()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class TargetImporter:
    """Import targets from RAG Chat target hypothesis exports."""

    def __init__(self, rag_data_dir: Path = None):
        """
        Args:
            rag_data_dir: Path to RAG Chat data directory containing targets
        """
        self.rag_data_dir = rag_data_dir or Path("/home/adam/transfer/rag-chat-pipeline/data/targets")

    def import_from_export(self, export_file: Path) -> List[ImportedTarget]:
        """Import targets from a Phase 5 export JSON file."""
        if not export_file.exists():
            raise FileNotFoundError(f"Export file not found: {export_file}")

        with open(export_file, 'r') as f:
            data = json.load(f)

        targets = []
        for item in data.get('targets', []):
            target = ImportedTarget(
                gene=item.get('gene', 'Unknown'),
                protein=item.get('protein'),
                rationale=item.get('rationale', ''),
                confidence=item.get('confidence', 'medium'),
                priority=item.get('priority', 3),
                therapeutic_area=item.get('therapeutic_area'),
                mechanism=item.get('mechanism'),
                pdb_ids=item.get('pdb_ids', []),
                variant_count=item.get('variant_count', 0),
                variants=item.get('variants', []),
                source_query=item.get('source_query'),
                status=item.get('status', 'imported'),
            )
            targets.append(target)

        return targets

    def import_all(self) -> List[ImportedTarget]:
        """Import all targets from RAG Chat data directory."""
        targets = []

        # Look for individual target JSON files
        if self.rag_data_dir.exists():
            for json_file in self.rag_data_dir.glob("*.json"):
                if json_file.name.startswith("target_"):
                    try:
                        with open(json_file, 'r') as f:
                            data = json.load(f)
                        target = ImportedTarget(
                            gene=data.get('gene', 'Unknown'),
                            protein=data.get('protein'),
                            rationale=data.get('rationale', ''),
                            confidence=data.get('confidence', 'medium'),
                            priority=data.get('priority', 3),
                            therapeutic_area=data.get('therapeutic_area'),
                            mechanism=data.get('mechanism'),
                            pdb_ids=data.get('pdb_ids', []),
                            variant_count=data.get('variant_count', 0),
                            variants=data.get('variants', []),
                            source_query=data.get('source_query'),
                            status=data.get('status', 'imported'),
                        )
                        targets.append(target)
                    except (json.JSONDecodeError, KeyError) as e:
                        print(f"Warning: Could not parse {json_file}: {e}")

        return targets

    def create_vcp_demo_target(self) -> ImportedTarget:
        """Create a pre-configured VCP target for the FTD demo."""
        return ImportedTarget(
            gene="VCP",
            protein="Valosin-containing protein (p97/VCP)",
            rationale="""VCP/p97 is a critical AAA+ ATPase involved in protein quality control,
            autophagy, and endoplasmic reticulum-associated degradation (ERAD).
            Mutations in VCP cause frontotemporal dementia (FTD) with inclusion body myopathy
            and Paget's disease of bone (IBMPFD). The R159H variant found in the HG002 genome
            is located at the N-D1 domain interface, a region known to harbor disease-causing mutations.
            VCP is a validated drug target with existing small-molecule inhibitors in development.""",
            confidence="high",
            priority=5,
            therapeutic_area="Neurology / Neurodegeneration",
            mechanism="AAA+ ATPase modulation / Protein quality control restoration",
            pdb_ids=["8OOI", "9DIL", "7K56", "5FTK"],
            variant_count=10,
            variants=[
                {
                    "rsid": "rs188935092",
                    "chrom": "chr9",
                    "pos": 35056077,
                    "ref": "A",
                    "alt": "G",
                    "consequence": "missense_variant",
                    "impact": "HIGH",
                    "clinical_significance": "Conflicting_interpretations_of_pathogenicity",
                    "disease_association": "Frontotemporal dementia and/or amyotrophic lateral sclerosis 6"
                }
            ],
            source_query="What variants are found in the VCP gene? Are there any associated with frontotemporal dementia or ALS?",
            status="validated",
        )


def get_vcp_target() -> ImportedTarget:
    """Convenience function to get the VCP demo target."""
    importer = TargetImporter()
    return importer.create_vcp_demo_target()
