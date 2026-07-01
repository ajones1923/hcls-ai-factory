"""
Multi-omics per-patient join (F1).

Unifies a patient's genomics (variants + ACMG secondary findings), transcriptomics (single-cell
cell types + marker genes), and proteomics (folded structures + developability + targets) into one
record the intelligence agents can reason over. The payoff is **cross-omics convergence**: genes
that surface in more than one layer (e.g. a pathogenic variant whose gene is also a single-cell
marker and a protein target) are flagged as convergent signals — exactly the multi-modal evidence a
clinical agent should weight highest.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


def _genes(items: list[dict], key: str = "gene") -> set[str]:
    return {i[key].upper() for i in items if i.get(key)}


@dataclass
class MultiOmicsRecord:
    patient_id: str
    genomics: dict[str, Any] = field(default_factory=dict)        # {variants, secondary_findings}
    transcriptomics: dict[str, Any] = field(default_factory=dict)  # {cell_types, marker_genes}
    proteomics: dict[str, Any] = field(default_factory=dict)       # {structures, developability, targets}

    def _layers(self) -> dict[str, set[str]]:
        return {
            "genomics": _genes(self.genomics.get("variants", []))
            | _genes(self.genomics.get("secondary_findings", [])),
            "transcriptomics": {m.upper() for m in self.transcriptomics.get("marker_genes", [])},
            "proteomics": {t.upper() for t in self.proteomics.get("targets", [])},
        }

    def genes(self) -> set[str]:
        out: set[str] = set()
        for s in self._layers().values():
            out |= s
        return out

    def cross_omics_links(self) -> list[dict[str, Any]]:
        """Genes present in >1 omics layer — convergent multi-modal signals, ranked by breadth."""
        layers = self._layers()
        links = []
        for gene in sorted(self.genes()):
            present = [name for name, genes in layers.items() if gene in genes]
            if len(present) > 1:
                links.append({"gene": gene, "layers": present, "n_layers": len(present)})
        links.sort(key=lambda l: -l["n_layers"])
        return links

    def summary(self) -> dict[str, Any]:
        return {
            "patient_id": self.patient_id,
            "n_variants": len(self.genomics.get("variants", [])),
            "n_secondary_findings": len(self.genomics.get("secondary_findings", [])),
            "cell_types": self.transcriptomics.get("cell_types", []),
            "n_protein_findings": len(self.proteomics.get("structures", [])),
            "cross_omics_links": self.cross_omics_links(),
        }


class MultiOmicsStore:
    def __init__(self) -> None:
        self._recs: dict[str, MultiOmicsRecord] = {}

    def _rec(self, patient_id: str) -> MultiOmicsRecord:
        return self._recs.setdefault(patient_id, MultiOmicsRecord(patient_id))

    def add_genomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self._rec(patient_id).genomics.update(data)

    def add_transcriptomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self._rec(patient_id).transcriptomics.update(data)

    def add_proteomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self._rec(patient_id).proteomics.update(data)

    def get(self, patient_id: str) -> MultiOmicsRecord | None:
        return self._recs.get(patient_id)

    def patients(self) -> list[str]:
        return sorted(self._recs)
