"""Orphanet parser for the Rare Disease Diagnostic Agent.

Parses Orphanet disease data and seeds 30+ rare diseases with structured
data including ORPHA code, name, prevalence, inheritance pattern, and
age of onset information.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: 35 ORPHANET RARE DISEASES
# ===================================================================

ORPHANET_DISEASES: List[Dict[str, Any]] = [
    {"orpha_code": "ORPHA:586", "name": "Cystic Fibrosis", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy, Neonatal", "gene": "CFTR"},
    {"orpha_code": "ORPHA:98896", "name": "Duchenne Muscular Dystrophy", "prevalence": "1-5 / 10 000", "inheritance": "X-linked recessive", "age_onset": "Childhood", "gene": "DMD"},
    {"orpha_code": "ORPHA:248", "name": "Huntington Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal dominant", "age_onset": "Adult", "gene": "HTT"},
    {"orpha_code": "ORPHA:355", "name": "Gaucher Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Childhood, Adult", "gene": "GBA"},
    {"orpha_code": "ORPHA:365", "name": "Pompe Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "All ages", "gene": "GAA"},
    {"orpha_code": "ORPHA:716", "name": "Phenylketonuria", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal recessive", "age_onset": "Neonatal", "gene": "PAH"},
    {"orpha_code": "ORPHA:666", "name": "Osteogenesis Imperfecta", "prevalence": "6-9 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Neonatal, Infancy", "gene": "COL1A1"},
    {"orpha_code": "ORPHA:636", "name": "Neurofibromatosis Type 1", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Childhood", "gene": "NF1"},
    {"orpha_code": "ORPHA:805", "name": "Tuberous Sclerosis Complex", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Infancy, Childhood", "gene": "TSC1"},
    {"orpha_code": "ORPHA:70", "name": "Spinal Muscular Atrophy", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy", "gene": "SMN1"},
    {"orpha_code": "ORPHA:778", "name": "Rett Syndrome", "prevalence": "1-9 / 100 000", "inheritance": "X-linked dominant", "age_onset": "Infancy, Childhood", "gene": "MECP2"},
    {"orpha_code": "ORPHA:511", "name": "Maple Syrup Urine Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Neonatal", "gene": "BCKDHA"},
    {"orpha_code": "ORPHA:558", "name": "Marfan Syndrome", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Childhood, Adult", "gene": "FBN1"},
    {"orpha_code": "ORPHA:730", "name": "Polycystic Kidney Disease (AD)", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Adult", "gene": "PKD1"},
    {"orpha_code": "ORPHA:324", "name": "Fabry Disease", "prevalence": "1-5 / 10 000", "inheritance": "X-linked recessive", "age_onset": "Childhood, Adult", "gene": "GLA"},
    {"orpha_code": "ORPHA:91378", "name": "Hereditary Angioedema", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal dominant", "age_onset": "Childhood", "gene": "SERPING1"},
    {"orpha_code": "ORPHA:56", "name": "Alkaptonuria", "prevalence": "<1 / 1 000 000", "inheritance": "Autosomal recessive", "age_onset": "Childhood, Adult", "gene": "HGD"},
    {"orpha_code": "ORPHA:43", "name": "X-linked Adrenoleukodystrophy", "prevalence": "1-9 / 100 000", "inheritance": "X-linked recessive", "age_onset": "Childhood, Adult", "gene": "ABCD1"},
    {"orpha_code": "ORPHA:105400", "name": "Hereditary ATTR Amyloidosis", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal dominant", "age_onset": "Adult", "gene": "TTR"},
    {"orpha_code": "ORPHA:169805", "name": "Hemophilia A", "prevalence": "1-5 / 10 000", "inheritance": "X-linked recessive", "age_onset": "Infancy, Childhood", "gene": "F8"},
    {"orpha_code": "ORPHA:33069", "name": "Dravet Syndrome", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal dominant", "age_onset": "Infancy", "gene": "SCN1A"},
    {"orpha_code": "ORPHA:910", "name": "Xeroderma Pigmentosum", "prevalence": "<1 / 1 000 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy, Childhood", "gene": "XPA"},
    {"orpha_code": "ORPHA:534", "name": "MPS I (Hurler Syndrome)", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy", "gene": "IDUA"},
    {"orpha_code": "ORPHA:580", "name": "MPS II (Hunter Syndrome)", "prevalence": "1-9 / 100 000", "inheritance": "X-linked recessive", "age_onset": "Childhood", "gene": "IDS"},
    {"orpha_code": "ORPHA:908", "name": "Fragile X Syndrome", "prevalence": "1-5 / 10 000", "inheritance": "X-linked dominant", "age_onset": "Childhood", "gene": "FMR1"},
    {"orpha_code": "ORPHA:803", "name": "Tay-Sachs Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy", "gene": "HEXA"},
    {"orpha_code": "ORPHA:79269", "name": "Niemann-Pick Disease Type C", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Variable", "gene": "NPC1"},
    {"orpha_code": "ORPHA:79257", "name": "Krabbe Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy", "gene": "GALC"},
    {"orpha_code": "ORPHA:512", "name": "Medium-chain Acyl-CoA Dehydrogenase Deficiency", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal recessive", "age_onset": "Infancy", "gene": "ACADM"},
    {"orpha_code": "ORPHA:79259", "name": "Metachromatic Leukodystrophy", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Childhood", "gene": "ARSA"},
    {"orpha_code": "ORPHA:540", "name": "MPS IIIA (Sanfilippo A)", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Childhood", "gene": "SGSH"},
    {"orpha_code": "ORPHA:905", "name": "Wilson Disease", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Childhood, Adult", "gene": "ATP7B"},
    {"orpha_code": "ORPHA:79233", "name": "Propionic Acidemia", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Neonatal", "gene": "PCCA"},
    {"orpha_code": "ORPHA:79312", "name": "Methylmalonic Acidemia", "prevalence": "1-9 / 100 000", "inheritance": "Autosomal recessive", "age_onset": "Neonatal", "gene": "MUT"},
    {"orpha_code": "ORPHA:308", "name": "Ehlers-Danlos Syndrome (Hypermobile)", "prevalence": "1-5 / 10 000", "inheritance": "Autosomal dominant", "age_onset": "Childhood", "gene": "Unknown"},
]


# ===================================================================
# ORPHANET PARSER IMPLEMENTATION
# ===================================================================


class OrphanetParser(BaseIngestParser):
    """Parse Orphanet disease data for the Rare Disease Diagnostic Agent.

    In offline/seed mode, returns the curated ORPHANET_DISEASES list.

    Usage::

        parser = OrphanetParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        api_key: str | None = None,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="orphanet",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch Orphanet disease data.

        In seed mode (no API key), returns the curated ORPHANET_DISEASES list.

        Returns:
            List of raw Orphanet disease dictionaries.
        """
        if self.api_key:
            self.logger.info("Orphanet API key provided but live fetch not implemented; using seed data")

        self.logger.info("Using curated Orphanet seed data (%d diseases)", len(ORPHANET_DISEASES))
        return list(ORPHANET_DISEASES)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw Orphanet disease data into IngestRecord objects.

        Args:
            raw_data: List of Orphanet disease dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            orpha_code = entry.get("orpha_code", "")
            name = entry.get("name", "")
            prevalence = entry.get("prevalence", "")
            inheritance = entry.get("inheritance", "")
            age_onset = entry.get("age_onset", "")
            gene = entry.get("gene", "")

            text = (
                f"Orphanet Disease: {name} ({orpha_code}). "
                f"Gene: {gene}. "
                f"Inheritance: {inheritance}. "
                f"Age of onset: {age_onset}. "
                f"Prevalence: {prevalence}."
            )

            record = IngestRecord(
                text=text,
                metadata={
                    "orpha_code": orpha_code,
                    "name": name,
                    "prevalence": prevalence,
                    "inheritance": inheritance,
                    "age_onset": age_onset,
                    "gene": gene,
                    "source_db": "Orphanet",
                },
                collection_name="rd_diseases",
                record_id=f"ORPHA_{orpha_code}",
                source="orphanet",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate an Orphanet IngestRecord.

        Requirements:
            - text must be non-empty
            - must have orpha_code in metadata
            - must have name in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("orpha_code"):
            return False
        if not meta.get("name"):
            return False

        return True


def get_orphanet_disease_count() -> int:
    """Return the number of curated Orphanet diseases."""
    return len(ORPHANET_DISEASES)


def get_orphanet_genes() -> List[str]:
    """Return a deduplicated list of genes from Orphanet seed data."""
    genes = list({d["gene"] for d in ORPHANET_DISEASES if d["gene"] != "Unknown"})
    genes.sort()
    return genes
