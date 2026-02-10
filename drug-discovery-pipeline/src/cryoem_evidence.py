"""
Cryo-EM Evidence Module - Structural ground truth for drug discovery.

"Genomics tells us what changed, Cryo-EM shows us how it changed,
and generative AI helps us design what could fix it."
"""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict


@dataclass
class CryoEMStructure:
    """Cryo-EM structural evidence object."""
    structure_id: str
    protein: str
    gene: str
    method: str
    resolution: str
    conformation: str
    variant_context: str
    binding_sites: List[str]
    druggable_pockets: List[str]
    source: str
    summary_text: str
    pdb_url: Optional[str] = None
    full_name: Optional[str] = None
    organism: Optional[str] = "Homo sapiens"
    inhibitor: Optional[str] = None
    inhibitor_smiles: Optional[str] = None
    image_url: Optional[str] = None
    emdb_id: Optional[str] = None
    emdb_image_url: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def get_display_name(self) -> str:
        """Get a display-friendly name."""
        return f"{self.structure_id} - {self.protein} ({self.resolution})"


class CryoEMEvidenceManager:
    """Manager for Cryo-EM structural evidence."""

    def __init__(self, structures_dir: Path = None):
        self.structures_dir = structures_dir or Path(__file__).parent.parent / "data" / "structures"
        self._structures_cache: Dict[str, List[CryoEMStructure]] = {}

    def load_structures_for_gene(self, gene: str) -> List[CryoEMStructure]:
        """Load Cryo-EM structures for a specific gene."""
        gene_upper = gene.upper()

        if gene_upper in self._structures_cache:
            return self._structures_cache[gene_upper]

        structures = []

        # Look for gene-specific structure file
        structure_file = self.structures_dir / f"{gene.lower()}_structures.json"
        if structure_file.exists():
            with open(structure_file, 'r') as f:
                data = json.load(f)
            for item in data:
                structure = CryoEMStructure(
                    structure_id=item.get('structure_id', ''),
                    protein=item.get('protein', ''),
                    gene=item.get('gene', gene_upper),
                    method=item.get('method', 'Cryo-EM'),
                    resolution=item.get('resolution', 'N/A'),
                    conformation=item.get('conformation', ''),
                    variant_context=item.get('variant_context', ''),
                    binding_sites=item.get('binding_sites', []),
                    druggable_pockets=item.get('druggable_pockets', []),
                    source=item.get('source', 'RCSB PDB'),
                    summary_text=item.get('summary_text', ''),
                    pdb_url=item.get('pdb_url'),
                    full_name=item.get('full_name'),
                    organism=item.get('organism', 'Homo sapiens'),
                    inhibitor=item.get('inhibitor'),
                    inhibitor_smiles=item.get('inhibitor_smiles'),
                    image_url=item.get('image_url'),
                    emdb_id=item.get('emdb_id'),
                    emdb_image_url=item.get('emdb_image_url'),
                )
                structures.append(structure)

        self._structures_cache[gene_upper] = structures
        return structures

    def get_best_structure_for_drug_design(self, gene: str) -> Optional[CryoEMStructure]:
        """
        Get the best structure for drug design purposes.
        Prefers: highest resolution, inhibitor-bound, druggable pockets defined.
        """
        structures = self.load_structures_for_gene(gene)
        if not structures:
            return None

        # Score structures for drug design suitability
        def score_structure(s: CryoEMStructure) -> float:
            score = 0.0

            # Higher resolution is better (lower number)
            try:
                res = float(s.resolution.replace('Å', '').replace(' ', ''))
                score += max(0, 5 - res)  # Higher score for better resolution
            except ValueError:
                pass

            # Inhibitor-bound structures are valuable
            if s.inhibitor:
                score += 3.0

            # More druggable pockets is better
            score += len(s.druggable_pockets) * 0.5

            # Cryo-EM preferred for flexibility information
            if 'Cryo-EM' in s.method:
                score += 0.5

            return score

        return max(structures, key=score_structure)

    def get_inhibitor_template(self, gene: str) -> Optional[str]:
        """Get a known inhibitor SMILES to use as a seed for molecule generation."""
        structures = self.load_structures_for_gene(gene)
        for s in structures:
            if s.inhibitor_smiles:
                return s.inhibitor_smiles
        return None

    def auto_load_structures(self, gene: str, pdb_ids: List[str]) -> List[CryoEMStructure]:
        """
        Load structures from cache, falling back to creating entries from PDB IDs.

        First checks the local JSON cache. For any PDB IDs not found in cache,
        creates basic CryoEMStructure entries so the pipeline can proceed.

        Args:
            gene: Gene symbol (e.g. "VCP")
            pdb_ids: List of PDB IDs to load

        Returns:
            List of CryoEMStructure objects
        """
        # Start with any cached structures
        cached = self.load_structures_for_gene(gene)
        cached_ids = {s.structure_id.replace("PDB:", "").upper() for s in cached}

        # Add entries for any PDB IDs not already cached
        for pdb_id in pdb_ids:
            pdb_upper = pdb_id.upper().strip()
            if pdb_upper not in cached_ids:
                cached.append(CryoEMStructure(
                    structure_id=f"PDB:{pdb_upper}",
                    protein=gene,
                    gene=gene.upper(),
                    method="Unknown",
                    resolution="N/A",
                    conformation="Unknown",
                    variant_context="",
                    binding_sites=[],
                    druggable_pockets=[],
                    source="RCSB PDB",
                    summary_text=f"Structure {pdb_upper} for {gene}",
                    pdb_url=f"https://www.rcsb.org/structure/{pdb_upper}",
                ))

        # Update cache
        self._structures_cache[gene.upper()] = cached
        return cached

    def format_evidence_for_rag(self, structure: CryoEMStructure) -> Dict[str, Any]:
        """Format structure as evidence object compatible with RAG system."""
        return {
            "type": "cryo_em_structure",
            "id": structure.structure_id,
            "gene": structure.gene,
            "text_summary": structure.summary_text,
            "score": 0.95,  # High score for structural evidence
            "retrieval_method": "structural_evidence",
            "metadata": {
                "resolution": structure.resolution,
                "conformation": structure.conformation,
                "binding_sites": structure.binding_sites,
                "druggable_pockets": structure.druggable_pockets,
                "pdb_url": structure.pdb_url,
            }
        }


def get_vcp_structures() -> List[CryoEMStructure]:
    """Convenience function to get VCP structures for the FTD demo."""
    manager = CryoEMEvidenceManager()
    return manager.load_structures_for_gene("VCP")


def get_vcp_inhibitor_seed() -> Optional[str]:
    """Get the CB-5083 inhibitor SMILES for VCP molecule generation."""
    manager = CryoEMEvidenceManager()
    return manager.get_inhibitor_template("VCP")
