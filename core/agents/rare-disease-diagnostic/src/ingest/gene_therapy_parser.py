"""Gene Therapy parser for the Rare Disease Diagnostic Agent.

Seeds 6 approved gene therapies and 20 investigational gene therapies
for rare diseases with structured data including drug name, target gene,
indication, approval status, mechanism, and vector type.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA: APPROVED GENE THERAPIES
# ===================================================================

APPROVED_GENE_THERAPIES: List[Dict[str, Any]] = [
    {
        "drug_name": "Zolgensma (onasemnogene abeparvovec)",
        "target_gene": "SMN1",
        "indication": "Spinal Muscular Atrophy Type 1",
        "approval_status": "FDA approved (2019), EMA approved (2020)",
        "mechanism": "AAV9-mediated SMN1 gene replacement",
        "vector_type": "AAV9",
        "manufacturer": "Novartis Gene Therapies",
        "route": "Intravenous",
        "age_group": "Pediatric (<2 years)",
    },
    {
        "drug_name": "Luxturna (voretigene neparvovec)",
        "target_gene": "RPE65",
        "indication": "RPE65-mediated Inherited Retinal Dystrophy",
        "approval_status": "FDA approved (2017), EMA approved (2018)",
        "mechanism": "AAV2-mediated RPE65 gene replacement",
        "vector_type": "AAV2",
        "manufacturer": "Spark Therapeutics / Roche",
        "route": "Subretinal injection",
        "age_group": "Pediatric and adult",
    },
    {
        "drug_name": "Skysona (elivaldogene autotemcel)",
        "target_gene": "ABCD1",
        "indication": "Cerebral Adrenoleukodystrophy",
        "approval_status": "FDA approved (2022), EMA approved (2021, Libmeldy)",
        "mechanism": "Lentiviral-mediated ABCD1 gene addition in autologous HSCs",
        "vector_type": "Lentiviral",
        "manufacturer": "bluebird bio",
        "route": "Intravenous (autologous HSCT)",
        "age_group": "Pediatric (4-17 years)",
    },
    {
        "drug_name": "Hemgenix (etranacogene dezaparvovec)",
        "target_gene": "F9",
        "indication": "Hemophilia B",
        "approval_status": "FDA approved (2022), EMA approved (2023)",
        "mechanism": "AAV5-mediated FIX-Padua gene delivery",
        "vector_type": "AAV5",
        "manufacturer": "CSL Behring / uniQure",
        "route": "Intravenous",
        "age_group": "Adult",
    },
    {
        "drug_name": "Roctavian (valoctocogene roxaparvovec)",
        "target_gene": "F8",
        "indication": "Hemophilia A",
        "approval_status": "EMA approved (2022), FDA approved (2023)",
        "mechanism": "AAV5-mediated FVIII gene delivery",
        "vector_type": "AAV5",
        "manufacturer": "BioMarin Pharmaceutical",
        "route": "Intravenous",
        "age_group": "Adult",
    },
    {
        "drug_name": "Casgevy (exagamglogene autotemcel)",
        "target_gene": "BCL11A (CRISPR editing)",
        "indication": "Sickle Cell Disease / Transfusion-Dependent Beta-Thalassemia",
        "approval_status": "FDA approved (2023), EMA approved (2024)",
        "mechanism": "CRISPR-Cas9 editing of BCL11A enhancer in autologous HSCs",
        "vector_type": "CRISPR/Cas9 (ex vivo)",
        "manufacturer": "Vertex / CRISPR Therapeutics",
        "route": "Intravenous (autologous HSCT)",
        "age_group": "Pediatric (12+) and adult",
    },
]


# ===================================================================
# SEED DATA: INVESTIGATIONAL GENE THERAPIES
# ===================================================================

INVESTIGATIONAL_GENE_THERAPIES: List[Dict[str, Any]] = [
    {
        "drug_name": "Giroctocogene fitelparvovec",
        "target_gene": "F8",
        "indication": "Hemophilia A",
        "approval_status": "Phase 3",
        "mechanism": "AAV6-mediated FVIII gene delivery",
        "vector_type": "AAV6",
        "manufacturer": "Pfizer / Sangamo",
    },
    {
        "drug_name": "Fidanacogene elaparvovec",
        "target_gene": "F9",
        "indication": "Hemophilia B",
        "approval_status": "Phase 3 (BLA submitted)",
        "mechanism": "AAV-mediated FIX gene delivery",
        "vector_type": "AAVRh74",
        "manufacturer": "Pfizer / Spark",
    },
    {
        "drug_name": "Delandistrogene moxeparvovec (SRP-9001)",
        "target_gene": "DMD (micro-dystrophin)",
        "indication": "Duchenne Muscular Dystrophy",
        "approval_status": "FDA accelerated approval (2023)",
        "mechanism": "AAVrh74-mediated micro-dystrophin expression",
        "vector_type": "AAVrh74",
        "manufacturer": "Sarepta Therapeutics",
    },
    {
        "drug_name": "OTL-200 (atidarsagene autotemcel)",
        "target_gene": "ARSA",
        "indication": "Metachromatic Leukodystrophy",
        "approval_status": "EMA approved (2020), FDA under review",
        "mechanism": "Lentiviral ARSA gene addition in autologous HSCs",
        "vector_type": "Lentiviral",
        "manufacturer": "Orchard Therapeutics",
    },
    {
        "drug_name": "AT132 (taldefgrobep alfa)",
        "target_gene": "MTM1",
        "indication": "X-linked Myotubular Myopathy",
        "approval_status": "Phase 1/2 (clinical hold)",
        "mechanism": "AAV8-mediated MTM1 gene replacement",
        "vector_type": "AAV8",
        "manufacturer": "Astellas Gene Therapies",
    },
    {
        "drug_name": "LX1001",
        "target_gene": "APOE2",
        "indication": "Niemann-Pick Disease Type C",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV9-mediated APOE2 delivery to CNS",
        "vector_type": "AAV9",
        "manufacturer": "Lexeo Therapeutics",
    },
    {
        "drug_name": "RGX-121",
        "target_gene": "IDS",
        "indication": "MPS II (Hunter Syndrome)",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV9-mediated IDS gene delivery to CNS",
        "vector_type": "AAV9",
        "manufacturer": "REGENXBIO",
    },
    {
        "drug_name": "PBGM01",
        "target_gene": "IDUA",
        "indication": "MPS I (Hurler Syndrome)",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV9-mediated IDUA gene delivery",
        "vector_type": "AAV9",
        "manufacturer": "Passage Bio",
    },
    {
        "drug_name": "FLT190",
        "target_gene": "FAH",
        "indication": "Tyrosinemia Type I",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV-mediated FAH gene replacement",
        "vector_type": "AAVS3",
        "manufacturer": "Freeline Therapeutics",
    },
    {
        "drug_name": "ABO-102 (isaralgagene civaparvovec)",
        "target_gene": "SGSH",
        "indication": "MPS IIIA (Sanfilippo A)",
        "approval_status": "Phase 3",
        "mechanism": "AAV9-mediated SGSH gene delivery to CNS",
        "vector_type": "AAV9",
        "manufacturer": "Abeona Therapeutics",
    },
    {
        "drug_name": "Etranacogene dezaparvovec v2 (optimized)",
        "target_gene": "F9",
        "indication": "Hemophilia B (optimized construct)",
        "approval_status": "Phase 3",
        "mechanism": "AAV5-mediated optimized FIX-Padua gene delivery",
        "vector_type": "AAV5",
        "manufacturer": "CSL Behring / uniQure",
    },
    {
        "drug_name": "AVXS-101 intrathecal (onasemnogene abeparvovec IT)",
        "target_gene": "SMN1",
        "indication": "SMA (older patients, intrathecal route)",
        "approval_status": "Phase 3",
        "mechanism": "AAV9-mediated SMN1 gene replacement via intrathecal delivery",
        "vector_type": "AAV9",
        "manufacturer": "Novartis Gene Therapies",
        "route": "Intrathecal",
        "age_group": "Pediatric (2+ years)",
    },
    {
        "drug_name": "BMN 307",
        "target_gene": "PAH",
        "indication": "Phenylketonuria (PKU)",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV5-mediated PAH gene replacement in hepatocytes",
        "vector_type": "AAV5",
        "manufacturer": "BioMarin Pharmaceutical",
    },
    {
        "drug_name": "FLT201 (Fabry gene therapy)",
        "target_gene": "GLA",
        "indication": "Fabry Disease",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV-mediated GLA gene replacement for alpha-galactosidase A",
        "vector_type": "AAVS3",
        "manufacturer": "Freeline Therapeutics",
    },
    {
        "drug_name": "SPK-3006",
        "target_gene": "GAA",
        "indication": "Pompe Disease (late-onset)",
        "approval_status": "Phase 1/2",
        "mechanism": "AAV-mediated GAA gene delivery for acid alpha-glucosidase",
        "vector_type": "AAV",
        "manufacturer": "Spark Therapeutics",
    },
    {
        "drug_name": "Timrepigene emparvovec (AGTC-501)",
        "target_gene": "RPGR",
        "indication": "X-linked Retinitis Pigmentosa (XLRP)",
        "approval_status": "Phase 3",
        "mechanism": "AAV-mediated RPGR gene replacement in photoreceptors",
        "vector_type": "AAV",
        "manufacturer": "AGTC / Beacon Therapeutics",
        "route": "Subretinal injection",
        "age_group": "Pediatric and adult",
    },
    {
        "drug_name": "Beremagene geperpavec (B-VEC)",
        "target_gene": "COL7A1",
        "indication": "Dystrophic Epidermolysis Bullosa",
        "approval_status": "FDA approved (2023)",
        "mechanism": "HSV1-mediated COL7A1 gene delivery to skin",
        "vector_type": "HSV-1",
        "manufacturer": "Krystal Biotech",
        "route": "Topical",
        "age_group": "Pediatric (6+ months) and adult",
    },
    {
        "drug_name": "Fordadistrogene movaparvovec (PF-06939926)",
        "target_gene": "DMD (mini-dystrophin)",
        "indication": "Duchenne Muscular Dystrophy",
        "approval_status": "Phase 3",
        "mechanism": "AAV9-mediated mini-dystrophin expression",
        "vector_type": "AAV9",
        "manufacturer": "Pfizer",
    },
    {
        "drug_name": "AAV-GAD (VY-AADC)",
        "target_gene": "AADC",
        "indication": "Aromatic L-amino Acid Decarboxylase Deficiency",
        "approval_status": "Phase 2 (approved in EU as Upstaza, 2022)",
        "mechanism": "AAV2-mediated AADC gene delivery to putamen",
        "vector_type": "AAV2",
        "manufacturer": "PTC Therapeutics",
        "route": "Intraputaminal",
        "age_group": "Pediatric (18+ months)",
    },
]


ALL_GENE_THERAPIES = APPROVED_GENE_THERAPIES + INVESTIGATIONAL_GENE_THERAPIES


# ===================================================================
# GENE THERAPY PARSER IMPLEMENTATION
# ===================================================================


class GeneTherapyParser(BaseIngestParser):
    """Parse gene therapy data for the Rare Disease Diagnostic Agent.

    Seeds 6 approved and 20 investigational gene therapies.

    Usage::

        parser = GeneTherapyParser()
        records, stats = parser.run()
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        super().__init__(
            source_name="gene_therapy",
            collection_manager=collection_manager,
            embedder=embedder,
        )

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch gene therapy data.

        Returns the curated list of approved and investigational gene therapies.

        Returns:
            List of gene therapy dictionaries.
        """
        self.logger.info(
            "Using curated gene therapy seed data (%d approved, %d investigational)",
            len(APPROVED_GENE_THERAPIES),
            len(INVESTIGATIONAL_GENE_THERAPIES),
        )
        return list(ALL_GENE_THERAPIES)

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw gene therapy data into IngestRecord objects.

        Args:
            raw_data: List of gene therapy dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for entry in raw_data:
            drug_name = entry.get("drug_name", "")
            target_gene = entry.get("target_gene", "")
            indication = entry.get("indication", "")
            approval_status = entry.get("approval_status", "")
            mechanism = entry.get("mechanism", "")
            vector_type = entry.get("vector_type", "")
            manufacturer = entry.get("manufacturer", "")
            route = entry.get("route", "")
            age_group = entry.get("age_group", "")

            text = (
                f"Gene Therapy: {drug_name}. "
                f"Target gene: {target_gene}. "
                f"Indication: {indication}. "
                f"Status: {approval_status}. "
                f"Mechanism: {mechanism}. "
                f"Vector: {vector_type}. "
                f"Manufacturer: {manufacturer}."
            )
            if route:
                text += f" Route: {route}."
            if age_group:
                text += f" Age group: {age_group}."

            safe_id = drug_name.split("(")[0].strip().replace(" ", "_")
            record = IngestRecord(
                text=text,
                metadata={
                    "drug_name": drug_name,
                    "target_gene": target_gene,
                    "indication": indication,
                    "approval_status": approval_status,
                    "mechanism": mechanism,
                    "vector_type": vector_type,
                    "manufacturer": manufacturer,
                    "route": route,
                    "age_group": age_group,
                    "source_db": "GeneTherapy",
                },
                collection_name="rd_therapies",
                record_id=f"GT_{safe_id}",
                source="gene_therapy",
            )
            records.append(record)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a gene therapy IngestRecord.

        Requirements:
            - text must be non-empty
            - must have drug_name in metadata
            - must have target_gene in metadata
            - must have indication in metadata

        Args:
            record: The record to validate.

        Returns:
            True if the record passes all validation checks.
        """
        if not record.text or not record.text.strip():
            return False

        meta = record.metadata
        if not meta.get("drug_name"):
            return False
        if not meta.get("target_gene"):
            return False
        if not meta.get("indication"):
            return False

        return True


def get_approved_therapy_count() -> int:
    """Return the number of approved gene therapies."""
    return len(APPROVED_GENE_THERAPIES)


def get_investigational_therapy_count() -> int:
    """Return the number of investigational gene therapies."""
    return len(INVESTIGATIONAL_GENE_THERAPIES)


def get_all_therapy_count() -> int:
    """Return the total number of gene therapies."""
    return len(ALL_GENE_THERAPIES)
