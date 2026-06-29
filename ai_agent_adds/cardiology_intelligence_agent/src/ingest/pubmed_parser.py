"""PubMed cardiovascular literature ingest parser.

Fetches cardiology research articles from the NCBI E-Utilities API
(esearch + efetch), parses XML responses into IngestRecord instances,
and targets the ``cardio_literature`` Milvus collection.

Covers 12 cardiovascular MeSH term categories with configurable
date-range filtering.  Respects NCBI rate limits (3 req/s without
API key, 10 req/s with key).

API docs: https://www.ncbi.nlm.nih.gov/books/NBK25497/

Author: Adam Jones
Date: March 2026
"""

import os
import re
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

import requests

from .base import BaseIngestParser, IngestRecord


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"
EFETCH_URL = f"{EUTILS_BASE}/efetch.fcgi"

# NCBI rate limit: 3 requests/second without API key
RATE_LIMIT_DELAY = 0.34


# ═══════════════════════════════════════════════════════════════════════
# PUBMED CARDIOLOGY PARSER
# ═══════════════════════════════════════════════════════════════════════


class PubMedCardioParser(BaseIngestParser):
    """Ingest parser for PubMed cardiovascular literature.

    Searches PubMed using 12 cardiovascular MeSH terms, retrieves
    article metadata and abstracts, and produces IngestRecord instances
    for the ``cardio_literature`` collection.

    MeSH terms cover the full breadth of cardiovascular medicine:
    heart failure, coronary artery disease, arrhythmia, valvular
    disease, imaging, prevention, cardio-oncology, and interventional
    cardiology.

    Usage::

        parser = PubMedCardioParser()
        records = parser.run(max_results=3000, years_back=5)
    """

    CARDIO_MESH_TERMS: List[str] = [
        "Cardiovascular Diseases",
        "Heart Failure",
        "Coronary Artery Disease",
        "Atrial Fibrillation",
        "Myocardial Infarction",
        "Cardiomyopathies",
        "Heart Valve Diseases",
        "Cardiac Imaging Techniques",
        "Electrocardiography",
        "Cardio-Oncology",
        "Preventive Cardiology",
        "Interventional Cardiology",
    ]
    """MeSH terms used to build the PubMed search query."""

    # Map MeSH terms to cardiology subspecialties for metadata tagging
    _SUBSPECIALTY_MAP: Dict[str, str] = {
        "Cardiovascular Diseases": "general",
        "Heart Failure": "heart_failure",
        "Coronary Artery Disease": "interventional",
        "Atrial Fibrillation": "electrophysiology",
        "Myocardial Infarction": "acute_coronary",
        "Cardiomyopathies": "heart_failure",
        "Heart Valve Diseases": "valvular",
        "Cardiac Imaging Techniques": "imaging",
        "Electrocardiography": "electrophysiology",
        "Cardio-Oncology": "cardio_oncology",
        "Preventive Cardiology": "prevention",
        "Interventional Cardiology": "interventional",
    }

    # Study type detection patterns
    _STUDY_TYPE_PATTERNS: List[tuple] = [
        (re.compile(r"\bmeta[- ]?analysis\b", re.IGNORECASE), "meta-analysis"),
        (re.compile(r"\brandomized\b.*\btrial\b", re.IGNORECASE), "RCT"),
        (re.compile(r"\bcohort\b", re.IGNORECASE), "cohort"),
        (re.compile(r"\bcase[- ]?control\b", re.IGNORECASE), "case-control"),
        (re.compile(r"\bcross[- ]?sectional\b", re.IGNORECASE), "cross-sectional"),
        (re.compile(r"\bsystematic review\b", re.IGNORECASE), "systematic_review"),
        (re.compile(r"\bcase report\b", re.IGNORECASE), "case_report"),
        (re.compile(r"\bguideline\b", re.IGNORECASE), "guideline"),
        (re.compile(r"\breview\b", re.IGNORECASE), "review"),
    ]

    def __init__(self, api_key: Optional[str] = None):
        """Initialize the PubMed cardiology parser.

        Args:
            api_key: Optional NCBI API key for higher rate limits.
                Falls back to the ``NCBI_API_KEY`` environment variable.
        """
        super().__init__("cardio_literature")
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")

    # ── Fetch ─────────────────────────────────────────────────────────

    def fetch(
        self,
        max_results: int = 3000,
        years_back: int = 5,
        batch_size: int = 200,
    ) -> List[dict]:
        """Fetch cardiovascular articles from PubMed via E-Utilities.

        Performs a two-step retrieval:
          1. ``esearch`` to get PMIDs matching the cardiovascular query.
          2. ``efetch`` in batches to retrieve full article XML.

        Args:
            max_results: Maximum number of articles to retrieve.
            years_back: Number of years back from today to search.
            batch_size: Number of PMIDs per efetch request (max 500).

        Returns:
            List of article dictionaries extracted from PubMed XML.
        """
        query = self._build_search_query(self.CARDIO_MESH_TERMS, years_back)
        self.logger.info(f"PubMed search query: {query}")

        # Step 1: esearch to get PMIDs
        pmids = self._esearch(query, max_results)
        if not pmids:
            self.logger.warning("No PMIDs returned from esearch")
            return []
        self.logger.info(f"esearch returned {len(pmids)} PMIDs")

        # Step 2: efetch in batches
        articles: List[dict] = []
        for i in range(0, len(pmids), batch_size):
            batch_ids = pmids[i : i + batch_size]
            batch_articles = self._efetch(batch_ids)
            articles.extend(batch_articles)
            self.logger.info(
                f"efetch batch {i // batch_size + 1}: "
                f"retrieved {len(batch_articles)} articles "
                f"(total {len(articles)})"
            )
            time.sleep(RATE_LIMIT_DELAY)

        return articles[:max_results]

    # ── Parse ─────────────────────────────────────────────────────────

    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Parse PubMed article dictionaries into IngestRecord instances.

        Each article dictionary (produced by :meth:`_extract_article`)
        is converted to an IngestRecord with text = title + abstract and
        metadata matching the ``cardio_literature`` collection schema.

        Args:
            raw_data: List of article dictionaries from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances.
        """
        records: List[IngestRecord] = []

        for article in raw_data:
            try:
                title = article.get("title", "").strip()
                abstract = article.get("abstract", "").strip()
                if not title and not abstract:
                    continue

                # Build embedding text: title + abstract
                text = f"{title}. {abstract}" if abstract else title

                # Detect study type from title + abstract
                study_type = self._detect_study_type(text)

                # Detect subspecialty from MeSH terms
                mesh_terms = article.get("mesh_terms", [])
                subspecialty = self._detect_subspecialty(mesh_terms)

                metadata = {
                    "title": self.truncate(title, 510),
                    "abstract": self.truncate(abstract, 8190),
                    "authors": self.truncate(
                        article.get("authors", ""), 1020
                    ),
                    "journal": self.truncate(
                        article.get("journal", ""), 254
                    ),
                    "year": article.get("year", 0),
                    "pmid": article.get("pmid", ""),
                    "doi": article.get("doi", ""),
                    "mesh_terms": self.truncate(
                        "; ".join(mesh_terms), 2046
                    ),
                    "study_type": study_type,
                    "subspecialty": subspecialty,
                }

                records.append(
                    IngestRecord(
                        text=self.truncate(text, 8192),
                        metadata=metadata,
                        collection=self.collection,
                        source="PubMed",
                        source_id=article.get("pmid"),
                    )
                )

            except Exception as exc:
                pmid = article.get("pmid", "unknown")
                self.logger.warning(f"Failed to parse article {pmid}: {exc}")
                continue

        self.logger.info(
            f"Parsed {len(records)} IngestRecords from "
            f"{len(raw_data)} PubMed articles"
        )
        return records

    # ── Private helpers ───────────────────────────────────────────────

    def _build_search_query(
        self, mesh_terms: List[str], years_back: int
    ) -> str:
        """Build a PubMed search query from MeSH terms and date filter.

        Constructs a query of the form:
          ``("Term1"[MeSH] OR "Term2"[MeSH] OR ...) AND
          ("YYYY/01/01"[PDAT] : "3000"[PDAT])``

        Args:
            mesh_terms: List of MeSH terms to OR together.
            years_back: Number of years back from the current year.

        Returns:
            Formatted PubMed query string.
        """
        from datetime import datetime

        start_year = datetime.now().year - years_back
        mesh_clause = " OR ".join(
            f'"{term}"[MeSH]' for term in mesh_terms
        )
        date_clause = f'"{start_year}/01/01"[PDAT] : "3000"[PDAT]'
        return f"({mesh_clause}) AND ({date_clause})"

    def _esearch(self, query: str, max_results: int) -> List[str]:
        """Execute an esearch request and return PMIDs.

        Args:
            query: PubMed query string.
            max_results: Maximum number of PMIDs to retrieve.

        Returns:
            List of PMID strings.
        """
        params: Dict[str, Any] = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        response = requests.get(ESEARCH_URL, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        return data.get("esearchresult", {}).get("idlist", [])

    def _efetch(self, pmids: List[str]) -> List[dict]:
        """Fetch full article XML for a batch of PMIDs.

        Args:
            pmids: List of PMID strings to retrieve.

        Returns:
            List of article dictionaries extracted from the XML.
        """
        params: Dict[str, Any] = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "rettype": "abstract",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        response = requests.get(EFETCH_URL, params=params, timeout=60)
        response.raise_for_status()

        return self._parse_efetch_xml(response.text)

    def _parse_efetch_xml(self, xml_text: str) -> List[dict]:
        """Parse efetch XML response into article dictionaries.

        Args:
            xml_text: Raw XML string from efetch.

        Returns:
            List of article dictionaries with keys: title, abstract,
            authors, journal, year, pmid, doi, mesh_terms.
        """
        articles: List[dict] = []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as exc:
            self.logger.error(f"XML parse error: {exc}")
            return articles

        for article_elem in root.findall(".//PubmedArticle"):
            extracted = self._extract_article(article_elem)
            if extracted:
                articles.append(extracted)

        return articles

    def _extract_article(self, article_xml: ET.Element) -> Optional[dict]:
        """Extract structured fields from a single PubmedArticle element.

        Navigates the PubMed XML structure to pull out title, abstract,
        author list, journal, year, PMID, DOI, and MeSH terms.

        Args:
            article_xml: An ``<PubmedArticle>`` XML element.

        Returns:
            Dictionary with article fields, or ``None`` if the article
            lacks essential data (PMID + title).
        """
        medline = article_xml.find("MedlineCitation")
        if medline is None:
            return None

        # PMID
        pmid_elem = medline.find("PMID")
        pmid = pmid_elem.text if pmid_elem is not None else ""
        if not pmid:
            return None

        # Article element
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
            for abs_text in abstract_elem.findall("AbstractText"):
                label = abs_text.get("Label", "")
                text = abs_text.text or ""
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
        abstract = " ".join(abstract_parts)

        # Authors
        author_list = article.find("AuthorList")
        authors: List[str] = []
        if author_list is not None:
            for author in author_list.findall("Author"):
                last = author.findtext("LastName", "")
                initials = author.findtext("Initials", "")
                if last:
                    authors.append(f"{last} {initials}".strip())

        # Journal
        journal_elem = article.find("Journal/Title")
        journal = journal_elem.text if journal_elem is not None else ""

        # Year
        year = 0
        pub_date = article.find("Journal/JournalIssue/PubDate")
        if pub_date is not None:
            year_elem = pub_date.find("Year")
            if year_elem is not None and year_elem.text:
                try:
                    year = int(year_elem.text)
                except ValueError:
                    pass

        # DOI
        doi = ""
        article_data = article_xml.find("PubmedData")
        if article_data is not None:
            for id_elem in article_data.findall(
                "ArticleIdList/ArticleId"
            ):
                if id_elem.get("IdType") == "doi":
                    doi = id_elem.text or ""
                    break

        # MeSH terms
        mesh_terms: List[str] = []
        mesh_list = medline.find("MeshHeadingList")
        if mesh_list is not None:
            for heading in mesh_list.findall("MeshHeading"):
                desc = heading.find("DescriptorName")
                if desc is not None and desc.text:
                    mesh_terms.append(desc.text)

        return {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "authors": ", ".join(authors),
            "journal": journal,
            "year": year,
            "doi": doi,
            "mesh_terms": mesh_terms,
        }

    def _detect_study_type(self, text: str) -> str:
        """Detect study design type from article text.

        Scans the title + abstract against ordered regex patterns
        and returns the first match.  Order matters: more specific
        types (meta-analysis, RCT) are tested before generic (review).

        Args:
            text: Combined title and abstract text.

        Returns:
            Study type string (e.g. ``"RCT"``, ``"cohort"``), or
            ``"other"`` if no pattern matches.
        """
        for pattern, study_type in self._STUDY_TYPE_PATTERNS:
            if pattern.search(text):
                return study_type
        return "other"

    def _detect_subspecialty(self, mesh_terms: List[str]) -> str:
        """Detect cardiology subspecialty from MeSH terms.

        Maps each MeSH term through :attr:`_SUBSPECIALTY_MAP` and
        returns the first non-general match.

        Args:
            mesh_terms: List of MeSH descriptor strings.

        Returns:
            Subspecialty string (e.g. ``"heart_failure"``), or
            ``"general"`` if no specific match is found.
        """
        for term in mesh_terms:
            subspecialty = self._SUBSPECIALTY_MAP.get(term)
            if subspecialty and subspecialty != "general":
                return subspecialty
        return "general"
