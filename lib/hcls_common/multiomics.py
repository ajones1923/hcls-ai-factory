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
from typing import Any, ClassVar


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


# --------------------------------------------------------------------------- #
# PF-1 — PatientContext: the longitudinal per-patient memory (upgrade of the
# 3-layer MultiOmicsRecord to the eleven-plus layers the seventeen surfaces need).
#
# Additive: MultiOmicsRecord/Store above are unchanged. Every engine/agent that
# touches a patient reads from and writes to this object; the `artifacts` ledger
# makes it the patient's durable memory. Cross-omics convergence here is over the
# gene-bearing layers (the flagship "same gene across omics" — e.g. TSC2 as a
# variant, a DE *driver* gene, a structural target, and a pathway member); the
# entity-resolution key (semantic concepts, cell-types, HPO) is upgraded by the
# BioKey resolver (PF-5) and the convergence engine (Wave 4).
# --------------------------------------------------------------------------- #
@dataclass
class PatientContext:
    patient_id: str
    genomics: dict[str, Any] = field(default_factory=dict)              # SNV/indel, mosaic VAF, CNV/SV, STR, mtDNA, splice, secondary findings (E1,A7,A4,A5)
    pharmacogenomics: dict[str, Any] = field(default_factory=dict)      # diplotypes, phenotypes, effective/measured phenotype (A3; E1 produces)
    immunogenetics: dict[str, Any] = field(default_factory=dict)        # HLA class I/II (the missing producer) (E1 HLA-node -> A3/A4/A5/E7/A1)
    transcriptomics: dict[str, Any] = field(default_factory=dict)       # cell types, DE-driver genes, TME, subclones, purity (E8/A8)
    surface_proteome: dict[str, Any] = field(default_factory=dict)      # CITE-seq/ADT antigen density (E8/A8)
    immune_repertoire: dict[str, Any] = field(default_factory=dict)     # TCR/BCR clonotypes, expansion, clonotype<->phenotype (E8/A8,A1,A4)
    epigenomics: dict[str, Any] = field(default_factory=dict)           # cell-type-resolved ATAC accessibility, episignatures (E8,A7)
    proteomics_structural: dict[str, Any] = field(default_factory=dict) # structures, pockets, ddG, PAE, PTM, epitope, developability (E7,E3)
    imaging: dict[str, Any] = field(default_factory=dict)               # image_read, segmentation, RECIST, CAC, volumetrics, imaging_hpo_terms (E4/E6)
    clinical: dict[str, Any] = field(default_factory=dict)              # HPO (present+excluded), staging, prior therapy, labs, activity indices (A7,A4,A5,A6)
    pathway_activity: dict[str, Any] = field(default_factory=dict)      # functional readouts (mTOR pS6/p-4E-BP1, IFN signature) (TSC,A4,A8)
    trajectory: dict[str, Any] = field(default_factory=dict)            # time-ordered serial signals + slopes (EF drift, NfL, EDSS, UPDRS) (E6,A5,A2)
    biological_age: dict[str, Any] = field(default_factory=dict)        # clocks, aging velocity, frailty, convergence signals (A2)
    artifacts: list[str] = field(default_factory=list)                  # every Artifact id ever written for this patient (the ledger)

    # the layer names, in order (the 13 data layers) — a class constant, not a field
    LAYERS: ClassVar[tuple[str, ...]] = (
        "genomics", "pharmacogenomics", "immunogenetics", "transcriptomics",
        "surface_proteome", "immune_repertoire", "epigenomics", "proteomics_structural",
        "imaging", "clinical", "pathway_activity", "trajectory", "biological_age",
    )

    def _gene_layers(self) -> dict[str, set[str]]:
        """Gene sets per gene-bearing layer (for cross-omics convergence).

        The transcriptomics contribution is DE **driver** genes (the E8-flagged fix), not
        cell-type markers — markers alone must not fire the flagship convergence.
        """
        g: dict[str, set[str]] = {
            "genomics": (
                _genes(self.genomics.get("variants", []))
                | _genes(self.genomics.get("secondary_findings", []))
                | _genes(self.genomics.get("cnv", []))
                | _genes(self.genomics.get("str_repeats", []))
            ),
            "pharmacogenomics": {k.upper() for k in self.pharmacogenomics.get("diplotypes", {})}
            | {x.upper() for x in self.pharmacogenomics.get("genes", [])},
            "transcriptomics": {x.upper() for x in self.transcriptomics.get("driver_genes", [])},
            "epigenomics": {x.upper() for x in self.epigenomics.get("genes", [])},
            "proteomics_structural": {t.upper() for t in self.proteomics_structural.get("targets", [])},
            "pathway_activity": {x.upper() for x in self.pathway_activity.get("genes", [])},
            "clinical": {x.upper() for x in self.clinical.get("genes", [])},
        }
        return {k: v for k, v in g.items() if v}

    def genes(self) -> set[str]:
        out: set[str] = set()
        for s in self._gene_layers().values():
            out |= s
        return out

    def cross_omics_links(self) -> list[dict[str, Any]]:
        """Genes present in >1 layer — convergent multi-modal signals, ranked by breadth."""
        layers = self._gene_layers()
        links = []
        for gene in sorted(self.genes()):
            present = [name for name, genes in layers.items() if gene in genes]
            if len(present) > 1:
                links.append({"gene": gene, "layers": present, "n_layers": len(present)})
        links.sort(key=lambda l: -l["n_layers"])
        return links

    def add_artifact(self, artifact_id: str) -> None:
        """Append an Artifact id to the patient's ledger (dedup, order-preserving)."""
        if artifact_id not in self.artifacts:
            self.artifacts.append(artifact_id)

    def populated_layers(self) -> list[str]:
        return [name for name in self.LAYERS if getattr(self, name)]

    def summary(self) -> dict[str, Any]:
        return {
            "patient_id": self.patient_id,
            "populated_layers": self.populated_layers(),
            "n_variants": len(self.genomics.get("variants", [])),
            "n_artifacts": len(self.artifacts),
            "cross_omics_links": self.cross_omics_links(),
        }

    def to_dict(self) -> dict[str, Any]:
        return {"patient_id": self.patient_id, "artifacts": list(self.artifacts),
                **{name: getattr(self, name) for name in self.LAYERS}}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "PatientContext":
        return cls(
            patient_id=d["patient_id"],
            artifacts=list(d.get("artifacts", [])),
            **{name: dict(d.get(name, {})) for name in cls.LAYERS},
        )

    @classmethod
    def from_multiomics_record(cls, rec: MultiOmicsRecord) -> "PatientContext":
        """Bridge the legacy 3-layer record into the new context (genomics, transcriptomics,
        and proteomics -> proteomics_structural). Legacy transcriptomics carries marker genes,
        which by design do NOT contribute to convergence (only DE driver genes do)."""
        return cls(
            patient_id=rec.patient_id,
            genomics=dict(rec.genomics),
            transcriptomics=dict(rec.transcriptomics),
            proteomics_structural=dict(rec.proteomics),
        )


class PatientContextStore:
    """In-memory store of PatientContext records — the factory's shared patient memory.

    Every capability reads via ``get`` and writes via ``update_layer`` (or a named helper).
    ``add_therapeutics`` / ``add_pathway_activity`` are the context-thread writers the TSC/agent
    flows use; ``record_artifact`` chains an emitted Artifact into the patient's ledger.
    """

    def __init__(self) -> None:
        self._recs: dict[str, PatientContext] = {}

    def _rec(self, patient_id: str) -> PatientContext:
        return self._recs.setdefault(patient_id, PatientContext(patient_id))

    def update_layer(self, patient_id: str, layer: str, data: dict[str, Any]) -> None:
        if layer not in PatientContext.LAYERS:
            raise ValueError(f"unknown patient-context layer: {layer!r}")
        getattr(self._rec(patient_id), layer).update(data)

    # named convenience writers (the common context-thread entry points)
    def add_genomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "genomics", data)

    def add_pharmacogenomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "pharmacogenomics", data)

    def add_transcriptomics(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "transcriptomics", data)

    def add_imaging(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "imaging", data)

    def add_clinical(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "clinical", data)

    def add_pathway_activity(self, patient_id: str, data: dict[str, Any]) -> None:
        self.update_layer(patient_id, "pathway_activity", data)

    def add_therapeutics(self, patient_id: str, therapies: list[dict[str, Any]]) -> None:
        """Append prior/active therapy into the clinical layer's `prior_therapy` (context thread)."""
        clinical = self._rec(patient_id).clinical
        clinical.setdefault("prior_therapy", []).extend(therapies)

    def record_artifact(self, patient_id: str, artifact_id: str) -> None:
        self._rec(patient_id).add_artifact(artifact_id)

    def get(self, patient_id: str) -> PatientContext | None:
        return self._recs.get(patient_id)

    def patients(self) -> list[str]:
        return sorted(self._recs)
