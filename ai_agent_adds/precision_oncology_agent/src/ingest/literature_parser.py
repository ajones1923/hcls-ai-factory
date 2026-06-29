"""
PubMed oncology literature ingest pipeline.

Fetches precision oncology publications from PubMed via the E-utilities API
and normalises them into the ``onco_literature`` collection schema for the
Oncology Intelligence Agent's RAG knowledge base.

Follows the same structural pattern as the CAR-T Agent's literature
ingest for consistency across the HCLS AI Factory project.

Author: Adam Jones
Date: February 2026
"""

import logging
import re
from typing import Any, Dict, List

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Cancer-type and gene extraction patterns
# ---------------------------------------------------------------------------

CANCER_KEYWORDS = [
    "lung cancer", "nsclc", "sclc", "breast cancer", "colorectal cancer",
    "colon cancer", "melanoma", "pancreatic cancer", "ovarian cancer",
    "prostate cancer", "glioblastoma", "glioma", "leukemia", "leukaemia",
    "lymphoma", "myeloma", "hepatocellular", "liver cancer", "renal cell",
    "kidney cancer", "bladder cancer", "gastric cancer", "stomach cancer",
    "esophageal cancer", "thyroid cancer", "endometrial cancer",
    "cervical cancer", "head and neck cancer", "sarcoma", "mesothelioma",
    "cholangiocarcinoma", "neuroblastoma", "medulloblastoma",
    "acute myeloid leukemia", "aml", "chronic myeloid leukemia", "cml",
    "acute lymphoblastic leukemia", "all",
    "non-hodgkin lymphoma", "hodgkin lymphoma",
    "triple-negative breast cancer", "tnbc", "her2-positive",
]

COMMON_ONCOGENES = [
    "EGFR", "BRAF", "KRAS", "NRAS", "PIK3CA", "TP53", "ALK", "ROS1",
    "MET", "HER2", "ERBB2", "RET", "NTRK1", "NTRK2", "NTRK3", "FGFR1",
    "FGFR2", "FGFR3", "IDH1", "IDH2", "BRCA1", "BRCA2", "ATM", "CDK4",
    "CDK6", "PTEN", "APC", "VHL", "RB1", "NF1", "NF2", "KIT", "PDGFRA",
    "FLT3", "NPM1", "JAK2", "MPL", "CALR", "BCR-ABL", "ABL1", "NOTCH1",
    "FBXW7", "STK11", "KEAP1", "AKT1", "MTOR", "MAP2K1", "ARID1A",
    "SMAD4", "CTNNB1", "DDR2", "ESR1", "AR", "PALB2", "RAD51",
    "POLE", "MSH2", "MSH6", "MLH1", "PMS2",
]

# Pre-compile regex patterns
_cancer_pattern = re.compile(
    "|".join(re.escape(k) for k in CANCER_KEYWORDS),
    re.IGNORECASE,
)
_gene_pattern = re.compile(
    r"\b(" + "|".join(re.escape(g) for g in COMMON_ONCOGENES) + r")\b"
)


class PubMedIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for PubMed oncology literature.

    Populates the ``onco_literature`` Milvus collection with article
    abstracts, metadata, and extracted annotations (cancer type, gene
    mentions, keywords).
    """

    def __init__(self, collection_manager: Any, embedder: Any) -> None:
        super().__init__(
            collection_manager=collection_manager,
            embedder=embedder,
            collection_name="onco_literature",
        )

    # ------------------------------------------------------------------
    # Fetch
    # ------------------------------------------------------------------

    def fetch(
        self,
        query: str = "precision oncology targeted therapy biomarker",
        max_results: int = 5000,
    ) -> List[Dict]:
        """
        Fetch oncology articles from PubMed via the E-utilities API.

        Uses the ``src.utils.pubmed_client`` module for HTTP operations.

        Parameters
        ----------
        query : str
            PubMed search query (default covers broad precision oncology).
        max_results : int
            Maximum number of articles to retrieve (default 5000).

        Returns
        -------
        list of dict
            Raw article metadata dicts with keys: pmid, title, abstract,
            authors, journal, year, keywords, mesh_terms, etc.
        """
        try:
            from src.utils.pubmed_client import search_pubmed, fetch_articles
        except ImportError:
            logger.error(
                "Could not import src.utils.pubmed_client — "
                "ensure pubmed_client.py is available."
            )
            return []

        logger.info("Searching PubMed: query=%r, max_results=%d", query, max_results)

        try:
            pmids = search_pubmed(query=query, max_results=max_results)
        except Exception as exc:
            logger.error("PubMed search failed: %s", exc)
            return []

        if not pmids:
            logger.warning("No PMIDs returned for query: %s", query)
            return []

        logger.info("Retrieved %d PMIDs, fetching article metadata...", len(pmids))

        try:
            articles = fetch_articles(pmids=pmids)
        except Exception as exc:
            logger.error("PubMed article fetch failed: %s", exc)
            return []

        logger.info("Fetched metadata for %d articles", len(articles))
        return articles

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse PubMed article records into ``onco_literature`` schema.

        Schema fields:
            - id: PMID (string)
            - title: Article title
            - text_chunk: Abstract text (primary embedding source)
            - text: Same as text_chunk (alias for base class)
            - source_type: Always "pubmed"
            - year: Publication year
            - cancer_type: Extracted from title/abstract
            - gene: Extracted gene mentions from title/abstract
            - keywords: MeSH terms + author keywords
            - authors: Author list (semicolon-separated)
            - journal: Journal name

        Parameters
        ----------
        raw_data : list of dict
            Raw PubMed article records from ``fetch``.

        Returns
        -------
        list of dict
            Normalised records ready for embedding and insertion.
        """
        records: List[Dict] = []

        for article in raw_data:
            pmid = str(article.get("pmid", article.get("id", "")))
            title = article.get("title", "")
            abstract = article.get("abstract", "")

            # Skip articles without meaningful text
            if not abstract and not title:
                continue

            # Build the text chunk: title + abstract
            text_chunk = f"{title}. {abstract}" if abstract else title

            # Extract year
            year = article.get("year", "")
            if not year:
                pub_date = article.get("pub_date", "")
                if pub_date and len(str(pub_date)) >= 4:
                    year = str(pub_date)[:4]

            # Extract cancer type from title + abstract
            cancer_type = self._extract_cancer_type(text_chunk)

            # Extract gene mentions
            gene = self._extract_genes(text_chunk)

            # Combine keywords from various sources
            keywords_list = article.get("keywords", [])
            mesh_terms = article.get("mesh_terms", [])
            if isinstance(keywords_list, str):
                keywords_list = [keywords_list]
            if isinstance(mesh_terms, str):
                mesh_terms = [mesh_terms]
            all_keywords = list(set(keywords_list + mesh_terms))
            keywords = "; ".join(all_keywords) if all_keywords else ""

            # Authors
            authors = article.get("authors", "")
            if isinstance(authors, list):
                authors = "; ".join(authors)

            journal = article.get("journal", "")

            records.append({
                "id": f"pmid_{pmid}",
                "title": title,
                "text_chunk": text_chunk,
                "text": text_chunk,
                "source_type": "pubmed",
                "year": str(year),
                "cancer_type": cancer_type,
                "gene": gene,
                "keywords": keywords,
                "authors": authors,
                "journal": journal,
            })

        logger.info(
            "Parsed %d literature records from %d raw articles",
            len(records), len(raw_data),
        )
        return records

    # ------------------------------------------------------------------
    # Extraction helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_cancer_type(text: str) -> str:
        """Extract the first matching cancer type from text."""
        match = _cancer_pattern.search(text)
        return match.group(0).lower() if match else ""

    @staticmethod
    def _extract_genes(text: str) -> str:
        """Extract unique gene symbols mentioned in text."""
        matches = _gene_pattern.findall(text)
        unique_genes = list(dict.fromkeys(matches))  # preserve order, deduplicate
        return ", ".join(unique_genes) if unique_genes else ""
