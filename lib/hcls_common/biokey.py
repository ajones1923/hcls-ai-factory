"""
BioKey resolver (PF-5) — the entity-resolution join primitive.

Different identifiers that point to the SAME biology — a gene symbol and its previous symbols /
synonyms / protein-product names, a druggable target, a cell-type marker, a pathway node, an HLA
allele, a VRS variant id — resolve to ONE canonical ``BioKey`` so cross-omics convergence unifies
them instead of failing on a raw string mismatch (the E8-flagged failure: ``N-myc`` never meeting
``MYCN``, ``tuberin`` never meeting ``TSC2``).

This ships a **small, honest seed** of well-known aliases plus ``register_alias`` / ``load_aliases``
so the real tables (HGNC previous symbols, UniProt, cell-ontology, HPO) are loaded by the producing
engines later — no network, no heavy deps, consistent with hcls_common. The full ``ConvergenceReasoner``
(breadth × clinical weight × honesty floor over resolved entities) is Wave 4; this is the join key it
stands on.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class EntityKind(str, Enum):
    GENE = "gene"
    PROTEIN = "protein"
    CELL_TYPE = "cell_type"
    PATHWAY = "pathway"
    HPO = "hpo"
    HLA = "hla"
    VARIANT = "variant"      # a VRS id (E1's genomic join key)


@dataclass(frozen=True)
class BioKey:
    """A canonical, typed entity key. Two references to the same biology share one BioKey."""
    kind: EntityKind
    id: str

    def __str__(self) -> str:
        return f"{self.kind.value}:{self.id}"


# A small, deliberately-conservative seed (canonical HGNC symbol -> known aliases / product names).
# Extend with register_alias()/load_aliases() from real HGNC/UniProt tables — do not grow this by hand.
_SEED_GENE_ALIASES: dict[str, list[str]] = {
    "MYCN": ["N-MYC", "NMYC", "N-myc"],
    "TSC1": ["hamartin"],
    "TSC2": ["tuberin"],
}


def _norm(s: str) -> str:
    return s.strip().upper()


class BioKeyResolver:
    """Resolves raw identifiers to canonical ``BioKey``s.

    Gene resolution folds known aliases/previous-symbols to the canonical symbol; other kinds are
    normalized (upper-cased, trimmed) and passed through unless an alias is registered. HPO/VRS ids
    are kept verbatim (they are already canonical identifiers).
    """

    def __init__(self, gene_aliases: dict[str, list[str]] | None = None) -> None:
        # reverse index: normalized alias -> canonical (per kind)
        self._alias: dict[EntityKind, dict[str, str]] = {k: {} for k in EntityKind}
        self.load_aliases(gene_aliases if gene_aliases is not None else _SEED_GENE_ALIASES)

    def register_alias(self, alias: str, canonical: str, kind: EntityKind = EntityKind.GENE) -> None:
        self._alias[kind][_norm(alias)] = _norm(canonical)

    def load_aliases(self, mapping: dict[str, list[str]], kind: EntityKind = EntityKind.GENE) -> None:
        """Load a canonical -> [aliases] table (e.g. HGNC previous symbols)."""
        for canonical, aliases in mapping.items():
            for a in aliases:
                self.register_alias(a, canonical, kind)

    def resolve(self, entity: str, kind: EntityKind = EntityKind.GENE) -> BioKey:
        """Resolve a raw identifier to its canonical BioKey."""
        if kind in (EntityKind.HPO, EntityKind.VARIANT):
            return BioKey(kind, entity.strip())          # already-canonical ids, kept verbatim
        key = _norm(entity)
        canonical = self._alias[kind].get(key, key)
        return BioKey(kind, canonical)

    def resolve_gene(self, symbol: str) -> str:
        """Convenience: the canonical gene symbol (upper-cased, alias-folded)."""
        return self.resolve(symbol, EntityKind.GENE).id


# a shared default resolver (seeded); engines extend it via register_alias/load_aliases
DEFAULT_RESOLVER = BioKeyResolver()


def resolve(entity: str, kind: EntityKind = EntityKind.GENE) -> BioKey:
    """Resolve against the shared default resolver."""
    return DEFAULT_RESOLVER.resolve(entity, kind)
