"""PubMed clinical trial publication parser for the Clinical Trial Intelligence Agent.

Fetches and parses clinical trial publications from NCBI's E-utilities API,
extracting structured publication metadata for embedding and storage in Milvus.

Uses the Entrez E-utilities:
  - esearch: Find PMIDs matching a query
  - efetch: Retrieve article metadata in XML

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)

# ===================================================================
# CONSTANTS
# ===================================================================

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_MAX_RESULTS = 500
RATE_LIMIT_DELAY = 0.34  # ~3 requests/second without API key
BATCH_SIZE = 50  # PMIDs per efetch request

# Default query for clinical trial publications
DEFAULT_QUERY = (
    '"clinical trial"[Publication Type] AND '
    '("2024"[Date - Publication] : "3000"[Date - Publication])'
)

# Clinical trial MeSH terms for targeted searches
TRIAL_MESH_TERMS = [
    "Clinical Trial",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Randomized Controlled Trial",
    "Controlled Clinical Trial",
    "Multicenter Study",
    "Pragmatic Clinical Trial",
    "Adaptive Clinical Trial",
]


# ===================================================================
# PARSER IMPLEMENTATION
# ===================================================================


class PubMedTrialParser(BaseIngestParser):
    """Ingest parser for PubMed clinical trial publications via E-utilities.

    Searches PubMed for clinical trial publications, fetches article
    metadata (PMID, title, abstract, MeSH terms, journal, year),
    and produces IngestRecord objects for embedding and Milvus storage.

    Usage::

        parser = PubMedTrialParser()
        records, stats = parser.run(
            query='breast cancer AND "clinical trial"[pt]',
            max_results=100,
        )
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
        api_key: Optional[str] = None,
    ) -> None:
        """Initialize the PubMed parser.

        Args:
            collection_manager: Optional Milvus collection manager.
            embedder: Optional embedding model.
            api_key: Optional NCBI API key (increases rate limit to 10 req/s).
        """
        super().__init__(
            source_name="pubmed",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key
        self._last_request_time: float = 0.0

    def _rate_limit(self) -> None:
        """Enforce rate limiting between API requests."""
        delay = 0.1 if self.api_key else RATE_LIMIT_DELAY
        elapsed = time.time() - self._last_request_time
        if elapsed < delay:
            time.sleep(delay - elapsed)
        self._last_request_time = time.time()

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch clinical trial publications from PubMed.

        Args:
            query: PubMed search query string.
            max_results: Maximum number of articles to fetch (default 500).
            mesh_terms: Additional MeSH terms to include in the query.

        Returns:
            List of raw article dictionaries.
        """
        try:
            import requests
        except ImportError:
            self.logger.warning(
                "requests library not available -- returning empty results"
            )
            return []

        query = kwargs.get("query", DEFAULT_QUERY)
        max_results = kwargs.get("max_results", DEFAULT_MAX_RESULTS)
        mesh_terms = kwargs.get("mesh_terms", [])

        # Append MeSH terms to query if provided
        if mesh_terms:
            mesh_query = " OR ".join(
                f'"{term}"[MeSH Terms]' for term in mesh_terms
            )
            query = f"({query}) AND ({mesh_query})"

        # Step 1: Search for PMIDs
        search_params: Dict[str, Any] = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "date",
        }
        if self.api_key:
            search_params["api_key"] = self.api_key

        try:
            self._rate_limit()
            response = requests.get(ESEARCH_URL, params=search_params, timeout=30)
            response.raise_for_status()
            search_data = response.json()
        except Exception as exc:
            self.logger.error("PubMed search failed: %s", exc)
            return []

        pmids = search_data.get("esearchresult", {}).get("idlist", [])
        if not pmids:
            self.logger.info("No PubMed results for query: %s", query[:100])
            return []

        self.logger.info("Found %d PMIDs matching query", len(pmids))

        # Step 2: Fetch article details in batches
        all_articles: List[Dict[str, Any]] = []
        for i in range(0, len(pmids), BATCH_SIZE):
            batch = pmids[i : i + BATCH_SIZE]
            articles = self._fetch_batch(batch, requests)
            all_articles.extend(articles)

        return all_articles

    def _fetch_batch(
        self, pmids: List[str], requests_module: Any
    ) -> List[Dict[str, Any]]:
        """Fetch article metadata for a batch of PMIDs.

        Args:
            pmids: List of PubMed IDs.
            requests_module: The requests library module.

        Returns:
            List of parsed article dictionaries.
        """
        fetch_params: Dict[str, Any] = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "rettype": "abstract",
        }
        if self.api_key:
            fetch_params["api_key"] = self.api_key

        try:
            self._rate_limit()
            response = requests_module.get(
                EFETCH_URL, params=fetch_params, timeout=30
            )
            response.raise_for_status()
            return self._parse_xml_response(response.text)
        except Exception as exc:
            self.logger.error("PubMed fetch batch failed: %s", exc)
            return []

    def _parse_xml_response(self, xml_text: str) -> List[Dict[str, Any]]:
        """Parse PubMed XML response into article dictionaries."""
        articles: List[Dict[str, Any]] = []

        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as exc:
            self.logger.error("Failed to parse PubMed XML: %s", exc)
            return []

        for article_elem in root.findall(".//PubmedArticle"):
            try:
                article = self._extract_article(article_elem)
                if article:
                    articles.append(article)
            except Exception as exc:
                self.logger.warning("Failed to extract article: %s", exc)

        return articles

    def _extract_article(self, elem: ET.Element) -> Optional[Dict[str, Any]]:
        """Extract structured data from a PubmedArticle XML element."""
        medline = elem.find("MedlineCitation")
        if medline is None:
            return None

        pmid_elem = medline.find("PMID")
        pmid = pmid_elem.text if pmid_elem is not None else ""
        if not pmid:
            return None

        article = medline.find("Article")
        if article is None:
            return None

        # Title
        title_elem = article.find("ArticleTitle")
        title = title_elem.text if title_elem is not None else ""

        # Abstract
        abstract_parts: List[str] = []
        abstract_elem = article.find("Abstract")
        if abstract_elem is not None:
            for text_elem in abstract_elem.findall("AbstractText"):
                label = text_elem.get("Label", "")
                text = text_elem.text or ""
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
        abstract = " ".join(abstract_parts)

        # Journal
        journal_elem = article.find("Journal/Title")
        journal = journal_elem.text if journal_elem is not None else ""

        # Year
        year_elem = article.find("Journal/JournalIssue/PubDate/Year")
        year = year_elem.text if year_elem is not None else ""

        # MeSH terms
        mesh_terms: List[str] = []
        mesh_list = medline.find("MeshHeadingList")
        if mesh_list is not None:
            for heading in mesh_list.findall("MeshHeading/DescriptorName"):
                term = heading.text
                if term:
                    mesh_terms.append(term)

        # Publication types
        pub_types: List[str] = []
        for pt in article.findall("PublicationTypeList/PublicationType"):
            if pt.text:
                pub_types.append(pt.text)

        # Authors
        authors: List[str] = []
        author_list = article.find("AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                last = author.find("LastName")
                initials = author.find("Initials")
                if last is not None and last.text:
                    name = last.text
                    if initials is not None and initials.text:
                        name += f" {initials.text}"
                    authors.append(name)

        return {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "mesh_terms": mesh_terms,
            "publication_types": pub_types,
            "authors": authors[:10],  # Cap at 10 authors
        }

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw PubMed article data into IngestRecord objects.

        Args:
            raw_data: List of article dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for article in raw_data:
            try:
                pmid = article.get("pmid", "")
                title = article.get("title", "")
                abstract = article.get("abstract", "")
                journal = article.get("journal", "")
                year = article.get("year", "")
                mesh_terms = article.get("mesh_terms", [])
                pub_types = article.get("publication_types", [])
                authors = article.get("authors", [])

                # Compose text for embedding
                text_parts = [f"Title: {title}"]
                if abstract:
                    text_parts.append(f"Abstract: {abstract}")
                if mesh_terms:
                    text_parts.append(f"MeSH Terms: {', '.join(mesh_terms)}")

                text = "\n".join(text_parts)

                metadata = {
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract[:2000] if abstract else "",
                    "journal": journal,
                    "year": year,
                    "mesh_terms": mesh_terms,
                    "publication_types": pub_types,
                    "authors": authors,
                    "source": "pubmed",
                }

                records.append(
                    IngestRecord(
                        text=text,
                        metadata=metadata,
                        collection_name="trial_literature",
                        record_id=f"PMID:{pmid}",
                        source="pubmed",
                    )
                )
            except Exception as exc:
                self.logger.warning("Failed to parse article PMID %s: %s",
                                    article.get("pmid", "?"), exc)

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a PubMed IngestRecord.

        Args:
            record: The IngestRecord to validate.

        Returns:
            True if valid.
        """
        if not record.text or len(record.text.strip()) < 20:
            return False
        if not record.record_id:
            return False
        if not record.metadata.get("title"):
            return False
        return True
