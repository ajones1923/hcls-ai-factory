"""Population allele frequency ingest pipeline for PGx Intelligence Agent.

Parses population-specific allele frequency data from gnomAD and PharmGKB
sources into PopulationData models.  Population frequency data is critical
for understanding the prevalence of pharmacogenomic variants across
different ethnic groups and for guiding population-specific PGx testing.

gnomAD API: https://gnomad.broadinstitute.org/api
PharmGKB frequencies: https://api.pharmgkb.org/v1/data/alleleFrequency

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import PGxCollectionManager
from src.models import PopulationData

from .base import BaseIngestPipeline, _truncate_utf8


# API endpoints
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
PHARMGKB_FREQ_URL = "https://api.pharmgkb.org/v1/data/alleleFrequency"

# Population labels used in gnomAD v4
GNOMAD_POPULATIONS = [
    {"code": "afr", "label": "African/African American"},
    {"code": "amr", "label": "Latino/Admixed American"},
    {"code": "asj", "label": "Ashkenazi Jewish"},
    {"code": "eas", "label": "East Asian"},
    {"code": "fin", "label": "Finnish"},
    {"code": "nfe", "label": "Non-Finnish European"},
    {"code": "sas", "label": "South Asian"},
    {"code": "mid", "label": "Middle Eastern"},
]

# Core pharmacogenes for frequency retrieval
FREQUENCY_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6",
    "DPYD", "TPMT", "NUDT15", "UGT1A1", "SLCO1B1",
    "VKORC1", "CYP4F2", "NAT2", "G6PD", "IFNL3",
]


class PopulationFrequencyParser(BaseIngestPipeline):
    """Ingest pipeline for population-specific allele frequency data.

    Fetches allele frequencies from PharmGKB and gnomAD for core
    pharmacogenes, stratified by population group.  Stores results
    in the pgx_population_data Milvus collection.

    Usage:
        parser = PopulationFrequencyParser(collection_manager, embedder)
        count = parser.run()
    """

    COLLECTION_NAME = "pgx_population_data"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
        genes: Optional[List[str]] = None,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout
        self.genes = genes or FREQUENCY_GENES

    def fetch(self, **kwargs) -> Dict[str, Any]:
        """Fetch allele frequency data from PharmGKB.

        Returns:
            Dict with keys: pharmgkb_frequencies, indexed by gene.
        """
        genes = kwargs.get("genes", self.genes)
        raw: Dict[str, Any] = {"pharmgkb_frequencies": {}}

        for gene in genes:
            try:
                logger.info(f"Fetching allele frequencies for {gene}...")
                params = {"gene": gene, "view": "max"}
                resp = requests.get(
                    PHARMGKB_FREQ_URL,
                    params=params,
                    timeout=self.request_timeout,
                )
                resp.raise_for_status()
                data = resp.json()
                freq_data = data.get("data", data) if isinstance(data, dict) else data
                raw["pharmgkb_frequencies"][gene] = freq_data
                count = len(freq_data) if isinstance(freq_data, list) else 0
                logger.info(f"Fetched {count} frequency entries for {gene}")

            except requests.RequestException as e:
                logger.error(f"Failed to fetch frequencies for {gene}: {e}")
                continue

        return raw

    def parse(self, raw_data: Dict[str, Any]) -> List[PopulationData]:
        """Parse allele frequency data into PopulationData models."""
        records: List[PopulationData] = []

        pharmgkb_freqs = raw_data.get("pharmgkb_frequencies", {})
        for gene, freq_entries in pharmgkb_freqs.items():
            if not isinstance(freq_entries, list):
                continue

            for entry in freq_entries:
                try:
                    allele_name = entry.get("alleleName", "") or entry.get("name", "")
                    star_allele = allele_name if allele_name.startswith("*") else f"*{allele_name}"

                    # Process population-specific frequencies
                    populations = entry.get("populations", entry.get("frequencies", []))
                    if isinstance(populations, dict):
                        # Format: {"European": 0.15, "African": 0.05, ...}
                        for pop_name, freq in populations.items():
                            record = self._build_record(
                                gene, star_allele, pop_name, freq, entry
                            )
                            if record:
                                records.append(record)
                    elif isinstance(populations, list):
                        # Format: [{"population": "European", "frequency": 0.15}, ...]
                        for pop_entry in populations:
                            if isinstance(pop_entry, dict):
                                pop_name = pop_entry.get("population", "") or pop_entry.get("ethnicity", "")
                                freq = pop_entry.get("frequency", None) or pop_entry.get("alleleFrequency", None)
                                sample_size = pop_entry.get("sampleSize", 0) or pop_entry.get("subjectCount", 0)

                                record = self._build_record(
                                    gene, star_allele, pop_name, freq, entry, sample_size
                                )
                                if record:
                                    records.append(record)
                    else:
                        # Single global frequency
                        global_freq = entry.get("frequency", entry.get("alleleFrequency"))
                        if global_freq is not None:
                            record = self._build_record(
                                gene, star_allele, "Global", global_freq, entry
                            )
                            if record:
                                records.append(record)

                except Exception as e:
                    logger.warning(f"Failed to parse frequency entry for {gene}: {e}")
                    continue

        logger.info(f"Parsed {len(records)} population frequency records")
        return records

    def _build_record(
        self,
        gene: str,
        star_allele: str,
        population: str,
        frequency: Any,
        source_entry: Dict,
        sample_size: int = 0,
    ) -> Optional[PopulationData]:
        """Build a PopulationData record from parsed fields."""
        if not population:
            return None

        try:
            freq_val = float(frequency)
        except (TypeError, ValueError):
            return None

        if not (0.0 <= freq_val <= 1.0):
            return None

        source_study = source_entry.get("source", "") or source_entry.get("citation", "")
        pop_clean = population.replace(" ", "_")
        rec_id = _truncate_utf8(
            f"{gene}_{star_allele}_{pop_clean}".replace("*", "star")[:95], 95
        )

        text_chunk = (
            f"{gene} {star_allele} allele frequency in {population} population: "
            f"{freq_val:.4f}. "
            f"Sample size: {sample_size}. "
            f"Source: {source_study}"
        ).strip()

        return PopulationData(
            id=rec_id,
            gene=_truncate_utf8(gene, 48),
            star_allele=_truncate_utf8(star_allele, 48),
            population=_truncate_utf8(population, 95),
            allele_frequency=freq_val,
            sample_size=sample_size,
            source_study=_truncate_utf8(source_study, 195),
            text_chunk=_truncate_utf8(text_chunk, 2990),
        )
