"""NCBI E-utilities API wrapper for PubMed literature retrieval.

Provides a client for searching PubMed and fetching abstract records
via the NCBI E-utilities API:
  - esearch: search PubMed and return PMIDs
  - efetch:  fetch full records (XML) for a list of PMIDs

Rate limits:
  - Without API key: 3 requests/second
  - With API key:   10 requests/second

API docs: https://www.ncbi.nlm.nih.gov/books/NBK25500/

Author: Adam Jones
Date: February 2026
"""

import os
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional
from urllib.parse import urlencode

import requests
from loguru import logger


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{EUTILS_BASE_URL}/esearch.fcgi"
EFETCH_URL = f"{EUTILS_BASE_URL}/efetch.fcgi"

# Rate limits (requests per second)
RATE_LIMIT_NO_KEY = 3
RATE_LIMIT_WITH_KEY = 10

# Default batch size for efetch (max 10,000 per NCBI docs)
EFETCH_BATCH_SIZE = 200


class PubMedClient:
    """Client for NCBI E-utilities PubMed API.

    Handles search (esearch) and abstract retrieval (efetch) with
    automatic rate limiting and pagination.

    Usage:
        client = PubMedClient(api_key="your_ncbi_key")
        pmids = client.search("medical imaging AI", max_results=100)
        articles = client.fetch_abstracts(pmids)
        for article in articles:
            print(article["title"], article["year"])
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        tool: str = "imaging_intelligence_agent",
    ):
        """Initialize the PubMed client.

        Args:
            api_key: NCBI API key for higher rate limits (10 req/sec).
                If None, falls back to NCBI_API_KEY environment variable.
                Without a key, rate limit is 3 req/sec.
            email: Contact email for NCBI (recommended by their policy).
                Falls back to NCBI_EMAIL environment variable.
            tool: Tool name sent in requests (NCBI tracking).
        """
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.email = email or os.environ.get("NCBI_EMAIL", "")
        self.tool = tool

        # Set rate limit based on API key presence
        if self.api_key:
            self._min_interval = 1.0 / RATE_LIMIT_WITH_KEY
            logger.info("PubMed client initialized with API key (10 req/sec)")
        else:
            self._min_interval = 1.0 / RATE_LIMIT_NO_KEY
            logger.info("PubMed client initialized without API key (3 req/sec)")

        self._last_request_time: float = 0.0

    def _build_base_params(self) -> Dict[str, str]:
        """Build the base query parameters shared by all E-utilities requests.

        Returns:
            Dict of common query parameters (tool, email, api_key if set).
        """
        params: Dict[str, str] = {
            "tool": self.tool,
        }
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def _rate_limit(self) -> None:
        """Enforce rate limiting between API requests.

        Sleeps if necessary to maintain the required minimum interval
        between consecutive requests (3/sec without key, 10/sec with key).
        """
        now = time.time()
        elapsed = now - self._last_request_time
        if elapsed < self._min_interval:
            sleep_time = self._min_interval - elapsed
            time.sleep(sleep_time)
        self._last_request_time = time.time()

    def search(
        self,
        query: str,
        max_results: int = 5000,
    ) -> List[str]:
        """Search PubMed and return matching PMIDs.

        Uses the esearch E-utility to find articles matching the query.
        Automatically paginates if max_results exceeds the per-request
        limit (10,000).

        Args:
            query: PubMed search query string.  Supports full PubMed
                query syntax including field tags, boolean operators,
                and MeSH terms.
                Examples:
                  - '"medical imaging" AND "deep learning"'
                  - '"chest x-ray classification"[Title/Abstract]'
                  - '"CT segmentation" AND "clinical trial"[pt]'
            max_results: Maximum number of PMIDs to retrieve.

        Returns:
            List of PMID strings (e.g. ["12345678", "23456789"]).
        """
        all_pmids: List[str] = []
        retmax = 10000
        retstart = 0

        try:
            # First request to get total count and initial batch
            params = self._build_base_params()
            params.update({
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retstart": str(retstart),
                "retmax": str(retmax),
            })
            url = f"{ESEARCH_URL}?{urlencode(params)}"

            self._rate_limit()
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            result = response.json()

            id_list = result["esearchresult"]["idlist"]
            total_count = int(result["esearchresult"]["count"])
            all_pmids.extend(id_list)

            logger.info(
                f"PubMed search '{query[:80]}' returned {total_count} total results"
            )

            # Paginate if there are more results than the first batch
            while len(all_pmids) < total_count and len(all_pmids) < max_results:
                retstart += retmax
                params = self._build_base_params()
                params.update({
                    "db": "pubmed",
                    "term": query,
                    "retmode": "json",
                    "retstart": str(retstart),
                    "retmax": str(retmax),
                })
                url = f"{ESEARCH_URL}?{urlencode(params)}"

                self._rate_limit()
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                result = response.json()

                id_list = result["esearchresult"]["idlist"]
                if not id_list:
                    break
                all_pmids.extend(id_list)

            # Trim to max_results
            all_pmids = all_pmids[:max_results]
            logger.info(f"Returning {len(all_pmids)} PMIDs")
            return all_pmids

        except requests.RequestException as e:
            logger.error(f"HTTP error during PubMed search: {e}")
            return all_pmids
        except (KeyError, ValueError) as e:
            logger.error(f"Error parsing PubMed search response: {e}")
            return all_pmids

    def fetch_abstracts(
        self,
        pmids: List[str],
        batch_size: int = EFETCH_BATCH_SIZE,
    ) -> List[Dict[str, Any]]:
        """Fetch detailed abstract records for a list of PMIDs.

        Uses the efetch E-utility to retrieve PubMed XML records,
        then parses them into structured dicts.  Processes PMIDs
        in batches to respect API limits.

        Args:
            pmids: List of PubMed IDs to fetch records for.
            batch_size: Number of PMIDs per efetch request (max ~10,000,
                default 200 for reliability).

        Returns:
            List of article dicts, each containing:
                - pmid (str): PubMed ID
                - title (str): Article title
                - abstract (str): Abstract text
                - authors (List[str]): Author names
                - journal (str): Journal title
                - year (str): Publication year
                - mesh_terms (List[str]): MeSH descriptor terms
        """
        if not pmids:
            return []

        articles: List[Dict[str, Any]] = []

        # Split PMIDs into batches
        batches = [
            pmids[i : i + batch_size] for i in range(0, len(pmids), batch_size)
        ]
        total_batches = len(batches)

        for batch_idx, batch in enumerate(batches, start=1):
            try:
                params = self._build_base_params()
                params.update({
                    "db": "pubmed",
                    "id": ",".join(batch),
                    "rettype": "xml",
                    "retmode": "xml",
                })
                url = f"{EFETCH_URL}?{urlencode(params)}"

                self._rate_limit()
                response = requests.get(url, timeout=60)
                response.raise_for_status()

                root = ET.fromstring(response.content)

                for article_elem in root.findall(".//PubmedArticle"):
                    article: Dict[str, Any] = {}

                    # PMID
                    pmid_elem = article_elem.find(".//PMID")
                    article["pmid"] = (
                        pmid_elem.text if pmid_elem is not None and pmid_elem.text else ""
                    )

                    # Title
                    title_elem = article_elem.find(".//ArticleTitle")
                    article["title"] = (
                        title_elem.text if title_elem is not None and title_elem.text else ""
                    )

                    # Abstract - join multiple AbstractText elements
                    abstract_parts = []
                    for abs_elem in article_elem.findall(".//AbstractText"):
                        if abs_elem.text:
                            label = abs_elem.get("Label")
                            if label:
                                abstract_parts.append(f"{label}: {abs_elem.text}")
                            else:
                                abstract_parts.append(abs_elem.text)
                    article["abstract"] = " ".join(abstract_parts)

                    # Authors
                    authors = []
                    for author_elem in article_elem.findall(".//Author"):
                        last_name_elem = author_elem.find("LastName")
                        fore_name_elem = author_elem.find("ForeName")
                        last_name = (
                            last_name_elem.text
                            if last_name_elem is not None and last_name_elem.text
                            else ""
                        )
                        fore_name = (
                            fore_name_elem.text
                            if fore_name_elem is not None and fore_name_elem.text
                            else ""
                        )
                        if last_name and fore_name:
                            authors.append(f"{last_name} {fore_name}")
                        elif last_name:
                            authors.append(last_name)
                    article["authors"] = authors

                    # Journal
                    journal_elem = article_elem.find(".//Journal/Title")
                    article["journal"] = (
                        journal_elem.text
                        if journal_elem is not None and journal_elem.text
                        else ""
                    )

                    # Year - try PubDate/Year first, fallback to MedlineDate
                    year_elem = article_elem.find(".//PubDate/Year")
                    if year_elem is not None and year_elem.text:
                        article["year"] = year_elem.text
                    else:
                        medline_date_elem = article_elem.find(".//PubDate/MedlineDate")
                        if (
                            medline_date_elem is not None
                            and medline_date_elem.text
                        ):
                            date_text = medline_date_elem.text
                            article["year"] = date_text[:4] if len(date_text) >= 4 else ""
                        else:
                            article["year"] = ""

                    # MeSH terms
                    mesh_terms = []
                    for mesh_elem in article_elem.findall(
                        ".//MeshHeading/DescriptorName"
                    ):
                        if mesh_elem.text:
                            mesh_terms.append(mesh_elem.text)
                    article["mesh_terms"] = mesh_terms

                    articles.append(article)

                logger.info(f"Fetched batch {batch_idx}/{total_batches}")

            except requests.RequestException as e:
                logger.error(
                    f"HTTP error fetching batch {batch_idx}/{total_batches}: {e}"
                )
            except ET.ParseError as e:
                logger.error(
                    f"XML parse error for batch {batch_idx}/{total_batches}: {e}"
                )

        logger.info(f"Fetched {len(articles)} articles total from {len(pmids)} PMIDs")
        return articles
