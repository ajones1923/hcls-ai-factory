"""PharmVar allele definition ingest pipeline for PGx Intelligence Agent.

Parses PharmVar star allele definition tables, haplotype definitions,
and functional annotations into GeneReference models.  PharmVar is the
central repository for pharmacogene variation, providing standardized
star allele nomenclature for CYP enzymes and other pharmacogenes.

PharmVar API: https://www.pharmvar.org/api-service

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional

import requests
from loguru import logger
from pydantic import BaseModel

from src.collections import PGxCollectionManager
from src.models import GeneReference

from .base import BaseIngestPipeline, _truncate_utf8


# PharmVar API endpoints
PHARMVAR_BASE_URL = "https://www.pharmvar.org/api-service"
PHARMVAR_GENES_URL = f"{PHARMVAR_BASE_URL}/genes"
PHARMVAR_ALLELES_URL = f"{PHARMVAR_BASE_URL}/alleles"

# Core pharmacogenes tracked by PharmVar
CORE_PHARMACOGENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6",
    "CYP3A4", "CYP1A2", "CYP2A6", "CYP2E1", "CYP4F2",
    "DPYD", "TPMT", "NUDT15", "UGT1A1", "NAT2",
    "SLCO1B1", "ABCG2", "VKORC1", "IFNL3",
]


class PharmVarParser(BaseIngestPipeline):
    """Ingest pipeline for PharmVar star allele definitions.

    Fetches gene and allele data from the PharmVar API, parses star
    allele nomenclature, defining variants, and functional annotations
    into GeneReference models for the pgx_gene_reference collection.

    Usage:
        parser = PharmVarParser(collection_manager, embedder)
        count = parser.run()
    """

    COLLECTION_NAME = "pgx_gene_reference"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
        genes: Optional[List[str]] = None,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout
        self.genes = genes or CORE_PHARMACOGENES

    def fetch(self, **kwargs) -> List[Dict[str, Any]]:
        """Fetch allele definitions from PharmVar for each pharmacogene.

        Returns:
            List of dicts, each containing gene and allele data from the API.
        """
        genes = kwargs.get("genes", self.genes)
        all_allele_data: List[Dict[str, Any]] = []

        for gene in genes:
            try:
                url = f"{PHARMVAR_ALLELES_URL}/{gene}"
                logger.info(f"Fetching PharmVar alleles for {gene}...")
                resp = requests.get(
                    url,
                    headers={"Accept": "application/json"},
                    timeout=self.request_timeout,
                )
                resp.raise_for_status()
                alleles = resp.json()

                if isinstance(alleles, list):
                    for allele in alleles:
                        allele["_gene"] = gene
                    all_allele_data.extend(alleles)
                    logger.info(f"Fetched {len(alleles)} alleles for {gene}")
                elif isinstance(alleles, dict) and "data" in alleles:
                    for allele in alleles["data"]:
                        allele["_gene"] = gene
                    all_allele_data.extend(alleles["data"])
                    logger.info(f"Fetched {len(alleles['data'])} alleles for {gene}")

            except requests.RequestException as e:
                logger.error(f"Failed to fetch PharmVar alleles for {gene}: {e}")
                continue

        logger.info(f"Fetched {len(all_allele_data)} total alleles across {len(genes)} genes")
        return all_allele_data

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[GeneReference]:
        """Parse PharmVar allele data into GeneReference models."""
        records: List[GeneReference] = []

        for allele in raw_data:
            try:
                gene = allele.get("_gene", "")
                allele_name = allele.get("alleleName", "") or allele.get("name", "")
                star_allele = allele_name if allele_name.startswith("*") else f"*{allele_name}"

                # Extract defining variant rsIDs
                variants = allele.get("variants", [])
                if isinstance(variants, list):
                    rsids = [v.get("rsid", "") for v in variants if isinstance(v, dict) and v.get("rsid")]
                    defining_variants = ", ".join(rsids)
                elif isinstance(variants, str):
                    defining_variants = variants
                else:
                    defining_variants = ""

                # Activity score and function
                activity_score = allele.get("activityScore")
                if activity_score is not None:
                    try:
                        activity_score = float(activity_score)
                    except (TypeError, ValueError):
                        activity_score = None

                function_status = allele.get("function", "") or allele.get("functionalStatus", "") or ""
                pharmvar_id = str(allele.get("pharmvarId", "") or allele.get("id", ""))

                # Allele frequencies (may not be present)
                freq_global = self._safe_float(allele.get("frequencyGlobal"))
                freq_european = self._safe_float(allele.get("frequencyEuropean"))
                freq_african = self._safe_float(allele.get("frequencyAfrican"))
                freq_east_asian = self._safe_float(allele.get("frequencyEastAsian"))
                freq_south_asian = self._safe_float(allele.get("frequencySouthAsian"))
                freq_latino = self._safe_float(allele.get("frequencyLatino"))

                # Build text chunk for embedding
                text_parts = [f"{gene} {star_allele}"]
                if function_status:
                    text_parts.append(f"Function: {function_status}")
                if defining_variants:
                    text_parts.append(f"Defining variants: {defining_variants}")
                if activity_score is not None:
                    text_parts.append(f"Activity score: {activity_score}")
                text_chunk = ". ".join(text_parts)

                rec_id = _truncate_utf8(f"{gene}_{star_allele}".replace("*", "star"), 95)

                record = GeneReference(
                    id=rec_id,
                    gene=_truncate_utf8(gene, 48),
                    star_allele=_truncate_utf8(star_allele, 48),
                    defining_variants=_truncate_utf8(defining_variants, 990),
                    activity_score=activity_score,
                    function_status=_truncate_utf8(function_status, 95),
                    allele_frequency_global=freq_global,
                    allele_frequency_european=freq_european,
                    allele_frequency_african=freq_african,
                    allele_frequency_east_asian=freq_east_asian,
                    allele_frequency_south_asian=freq_south_asian,
                    allele_frequency_latino=freq_latino,
                    pharmvar_id=_truncate_utf8(pharmvar_id, 48),
                    text_chunk=_truncate_utf8(text_chunk, 2990),
                    source="PharmVar",
                )
                records.append(record)

            except Exception as e:
                logger.warning(
                    f"Failed to parse PharmVar allele "
                    f"{allele.get('_gene', '?')} {allele.get('alleleName', '?')}: {e}"
                )
                continue

        logger.info(f"Parsed {len(records)} PharmVar gene reference records")
        return records

    @staticmethod
    def _safe_float(val: Any) -> Optional[float]:
        """Safely convert a value to float, returning None on failure."""
        if val is None:
            return None
        try:
            f = float(val)
            return f if 0.0 <= f <= 1.0 else None
        except (TypeError, ValueError):
            return None
