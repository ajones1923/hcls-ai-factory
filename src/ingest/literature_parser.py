"""PubMed literature ingest pipeline for CAR-T Intelligence Agent.

Fetches CAR-T research papers via NCBI E-utilities (esearch + efetch),
parses PubMed XML into CARTLiterature models, and stores embeddings
in the cart_literature Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import re
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import CARTLiterature, CARTStage, SourceType
from src.utils.pubmed_client import PubMedClient

from .base import BaseIngestPipeline


# ── Keyword sets for CAR-T stage classification ─────────────────────

_STAGE_KEYWORDS: Dict[CARTStage, List[str]] = {
    CARTStage.TARGET_ID: [
        "target identification",
        "target discovery",
        "antigen discovery",
        "tumor-associated antigen",
        "surface marker",
        "antigen screening",
        "neoantigen",
        "expression profiling",
        "single-cell RNA",
        "scRNA-seq",
        "proteomics",
        "immunopeptidomics",
    ],
    CARTStage.CAR_DESIGN: [
        "scfv",
        "single-chain variable fragment",
        "costimulatory domain",
        "4-1bb",
        "cd28",
        "car construct",
        "chimeric antigen receptor design",
        "nanobody",
        "vhh",
        "affinity maturation",
        "hinge region",
        "spacer",
        "transmembrane domain",
        "intracellular domain",
        "signaling domain",
    ],
    CARTStage.VECTOR_ENG: [
        "lentiviral vector",
        "retroviral vector",
        "viral vector",
        "vector production",
        "transduction efficiency",
        "crispr",
        "gene editing",
        "knock-in",
        "non-viral",
        "transposon",
        "sleeping beauty",
        "piggybac",
        "aav",
        "mrna delivery",
        "lipid nanoparticle",
    ],
    CARTStage.TESTING: [
        "cytotoxicity assay",
        "killing assay",
        "chromium release",
        "flow cytometry",
        "in vivo efficacy",
        "xenograft",
        "mouse model",
        "nsg mice",
        "tumor regression",
        "persistence",
        "exhaustion",
        "cytokine release",
        "ifn-gamma",
        "il-2",
        "t cell expansion",
    ],
    CARTStage.CLINICAL: [
        "clinical trial",
        "phase 1",
        "phase 2",
        "phase 3",
        "patient response",
        "complete remission",
        "partial remission",
        "overall survival",
        "progression-free",
        "cytokine release syndrome",
        "crs",
        "icans",
        "neurotoxicity",
        "dose escalation",
        "bridging therapy",
        "lymphodepletion",
        "fda approval",
        "bla",
    ],
}

# Keywords for extracting the target antigen from text
_ANTIGEN_PATTERNS: List[str] = [
    r"\bCD(?:19|20|22|30|33|38|123|138|269|276)\b",
    r"\bBCMA\b",
    r"\bGPC3\b",
    r"\bGD2\b",
    r"\bHER2\b",
    r"\bEGFR(?:vIII)?\b",
    r"\bMesothelin\b",
    r"\bCLAUDIN[\s-]?18\.2\b",
    r"\bPSMA\b",
    r"\bMUC1\b",
    r"\bROR1\b",
    r"\bCS1\b",
    r"\bSLAMF7\b",
    r"\bDLL3\b",
    r"\bFOLR1\b",
    r"\bTROP2\b",
    r"\bB7-H3\b",
    r"\bGPRC5D\b",
    r"\bFcRH5\b",
]

_ANTIGEN_REGEX = re.compile("|".join(_ANTIGEN_PATTERNS), re.IGNORECASE)

# Default PubMed query for CAR-T literature
DEFAULT_QUERY = (
    '"chimeric antigen receptor" OR "CAR-T" OR "CAR T cell"'
)


def _truncate_utf8(text: str, max_bytes: int) -> str:
    """Truncate a string to fit within max_bytes when UTF-8 encoded.

    Milvus VARCHAR max_length counts bytes, not characters, so we need
    byte-aware truncation for text containing multi-byte characters
    (Greek letters, mathematical symbols, CJK characters, etc.).
    """
    encoded = text.encode("utf-8")
    if len(encoded) <= max_bytes:
        return text
    # Truncate bytes and decode safely (ignoring partial multi-byte chars)
    return encoded[:max_bytes].decode("utf-8", errors="ignore")


class PubMedIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for PubMed CAR-T literature.

    Fetches abstracts from PubMed using NCBI E-utilities, classifies each
    paper by CAR-T development stage, extracts target antigens, and stores
    the results in the cart_literature Milvus collection.

    Usage:
        client = PubMedClient(api_key="...")
        pipeline = PubMedIngestPipeline(collection_manager, embedder, client)
        count = pipeline.run(query="CD19 CAR-T", max_results=500)
    """

    COLLECTION_NAME = "cart_literature"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        pubmed_client: Optional[PubMedClient] = None,
    ):
        """Initialize the PubMed ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            pubmed_client: PubMedClient instance.  If None, a default client
                is created.
        """
        super().__init__(collection_manager, embedder)
        self.pubmed_client = pubmed_client or PubMedClient()

    def fetch(
        self,
        query: str = DEFAULT_QUERY,
        max_results: int = 5000,
    ) -> List[Dict[str, Any]]:
        """Fetch abstracts from PubMed via NCBI E-utilities.

        Performs a two-step retrieval:
          1. esearch — get PMIDs matching the query
          2. efetch  — retrieve full abstract records for those PMIDs

        Args:
            query: PubMed search query string.
            max_results: Maximum number of articles to retrieve.

        Returns:
            List of dicts with keys: pmid, title, abstract, authors,
            journal, year, mesh_terms.

        """
        pmids = self.pubmed_client.search(query, max_results)
        logger.info(f"Found {len(pmids)} PMIDs")
        articles = self.pubmed_client.fetch_abstracts(pmids)
        return articles

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[CARTLiterature]:
        """Parse PubMed article dicts into CARTLiterature models.

        For each article, classifies the CAR-T development stage based on
        keyword analysis and extracts the target antigen if mentioned.

        Args:
            raw_data: List of dicts from fetch(), each containing:
                pmid, title, abstract, authors, journal, year, mesh_terms.

        Returns:
            List of validated CARTLiterature model instances.

        """
        records = []
        for article in raw_data:
            try:
                pmid = article.get("pmid", "")
                title = article.get("title", "")
                abstract = article.get("abstract", "")
                article.get("authors", [])
                journal = article.get("journal", "")
                year = article.get("year", None)
                mesh_terms = article.get("mesh_terms", [])

                # Combine title + abstract as the text chunk
                # Truncate to 2990 UTF-8 bytes (Milvus VARCHAR counts bytes)
                text_chunk = _truncate_utf8(f"{title} {abstract}".strip(), 2990)

                # Classify the CAR-T development stage
                cart_stage = self._classify_cart_stage(text_chunk)

                # Extract target antigen from text
                target_antigen = self._extract_target_antigen(text_chunk)

                # Convert year to int, default to 2020 if missing/invalid
                try:
                    year = int(year)
                except (TypeError, ValueError):
                    year = 2020

                # Join mesh terms list into a semicolon-separated string
                keywords = "; ".join(mesh_terms) if mesh_terms else ""

                record = CARTLiterature(
                    id=_truncate_utf8(pmid, 95),
                    title=_truncate_utf8(title, 490),
                    text_chunk=text_chunk,
                    source_type=SourceType.PUBMED,
                    year=year,
                    cart_stage=cart_stage,
                    target_antigen=_truncate_utf8(target_antigen, 95),
                    journal=_truncate_utf8(journal, 190),
                    keywords=_truncate_utf8(keywords, 950),
                )
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse article {article.get('pmid', '?')}: {e}")
                continue

        return records

    @staticmethod
    def _classify_cart_stage(abstract_text: str) -> CARTStage:
        """Classify a paper into a CAR-T development stage by keyword matching.

        Counts keyword hits for each CARTStage and returns the stage with
        the highest count.  Falls back to CARTStage.CLINICAL if no
        keywords match.

        Args:
            abstract_text: The abstract (or title + abstract) text to classify.

        Returns:
            The CARTStage with the most keyword matches.
        """
        if not abstract_text:
            return CARTStage.CLINICAL

        text_lower = abstract_text.lower()
        scores: Dict[CARTStage, int] = {}

        for stage, keywords in _STAGE_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw.lower() in text_lower)
            scores[stage] = count

        best_stage = max(scores, key=scores.get)  # type: ignore[arg-type]

        # If no keywords matched at all, default to CLINICAL
        if scores[best_stage] == 0:
            return CARTStage.CLINICAL

        return best_stage

    @staticmethod
    def _extract_target_antigen(abstract_text: str) -> str:
        """Extract the target antigen mentioned in abstract text.

        Uses regex patterns to find known CAR-T target antigens.
        Returns the first match found, or empty string if none.

        Args:
            abstract_text: The abstract text to search.

        Returns:
            The target antigen string (e.g. "CD19", "BCMA"), or "".
        """
        if not abstract_text:
            return ""

        match = _ANTIGEN_REGEX.search(abstract_text)
        if match:
            return match.group(0).upper()
        return ""

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full PubMed ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'cart_literature').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (query, max_results).

        Returns:
            Total number of records ingested.

        """
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
