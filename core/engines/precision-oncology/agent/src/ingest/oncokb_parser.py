"""
OncoKB ingest pipeline STUB.

OncoKB is a precision oncology knowledge base maintained by Memorial Sloan
Kettering Cancer Center.  Access to the OncoKB API requires a licensed API
key (free for academic / non-commercial use, commercial license required
otherwise).

This module provides the class skeleton so that the pipeline can be wired
up once a valid API key is available.  All fetch operations raise
``NotImplementedError`` until credentials are configured.

Reference: https://www.oncokb.org/
License info: https://www.oncokb.org/apiAccess

Author: Adam Jones
Date: February 2026
"""

import logging
from typing import Any, Dict, List, Optional

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)


class OncoKBIngestPipeline(BaseIngestPipeline):
    """
    Stub ingest pipeline for OncoKB variant and evidence data.

    Populates the ``onco_variants`` Milvus collection when a licensed
    API key is provided.  Until then, ``fetch`` raises
    ``NotImplementedError``.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance.
    embedder : object
        Embedding model with ``encode`` method.
    """

    def __init__(self, collection_manager: Any, embedder: Any) -> None:
        super().__init__(
            collection_manager=collection_manager,
            embedder=embedder,
            collection_name="onco_variants",
        )
        logger.warning(
            "OncoKBIngestPipeline initialised — note that OncoKB requires a "
            "licensed API key.  Visit https://www.oncokb.org/apiAccess to "
            "register for access before running this pipeline."
        )

    # ------------------------------------------------------------------
    # Fetch (not implemented)
    # ------------------------------------------------------------------

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict]:
        """
        Fetch variant / evidence data from OncoKB.

        Raises
        ------
        NotImplementedError
            Always — OncoKB requires a licensed API key.
        """
        raise NotImplementedError(
            "OncoKB requires a licensed API key. "
            "Visit https://www.oncokb.org/apiAccess to register for access. "
            "Once you have a key, set the ONCOKB_API_KEY environment variable "
            "and implement the fetch logic in this class."
        )

    # ------------------------------------------------------------------
    # Parse (pass-through)
    # ------------------------------------------------------------------

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse OncoKB records into ``onco_variants`` schema.

        Currently a pass-through since ``fetch`` is not yet implemented.
        When implemented, this should map OncoKB fields to:
            - id, gene, variant_name, variant_type, cancer_type
            - evidence_level (OncoKB levels 1-4, R1, R2)
            - drugs, oncokb_id, text_summary, source_type="oncokb"

        Parameters
        ----------
        raw_data : list of dict
            Raw OncoKB records.

        Returns
        -------
        list of dict
            Normalised records (pass-through for now).
        """
        return raw_data
