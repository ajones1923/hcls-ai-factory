"""
Simple PubMed E-utilities client.

Provides functions for searching PubMed via the NCBI E-utilities API
(esearch + efetch) and parsing the XML responses into structured article
metadata dictionaries.

NCBI API guidelines:
    - Without an API key: max 3 requests/second
    - With an API key: max 10 requests/second
    - Set api_key or NCBI_API_KEY environment variable

Reference: https://www.ncbi.nlm.nih.gov/books/NBK25497/

Author: Adam Jones
Date: February 2026
"""

import logging
import os
import time
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"
EFETCH_URL = f"{EUTILS_BASE}/efetch.fcgi"

BATCH_SIZE = 200  # PMIDs per efetch request
RATE_LIMIT_DELAY = 0.35  # seconds between requests (stay under 3/sec)


# ===================================================================
# Public API
# ===================================================================

def search_pubmed(
    query: str,
    max_results: int = 5000,
    api_key: Optional[str] = None,
) -> List[str]:
    """
    Search PubMed and return a list of PMIDs.

    Parameters
    ----------
    query : str
        PubMed search query string.
    max_results : int
        Maximum number of PMIDs to return (default 5000).
    api_key : str, optional
        NCBI API key. Falls back to ``NCBI_API_KEY`` environment variable.

    Returns
    -------
    list of str
        List of PMID strings.
    """
    api_key = api_key or os.environ.get("NCBI_API_KEY")

    all_pmids: List[str] = []
    retstart = 0

    logger.info("Searching PubMed: query=%r, max_results=%d", query, max_results)

    while retstart < max_results:
        params: Dict[str, str] = {
            "db": "pubmed",
            "term": query,
            "retmax": str(min(BATCH_SIZE, max_results - retstart)),
            "retstart": str(retstart),
            "retmode": "json",
            "sort": "relevance",
        }
        if api_key:
            params["api_key"] = api_key

        try:
            response = requests.get(ESEARCH_URL, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
        except requests.RequestException as exc:
            logger.error("PubMed esearch request failed: %s", exc)
            break

        result = data.get("esearchresult", {})
        id_list = result.get("idlist", [])

        if not id_list:
            break

        all_pmids.extend(id_list)
        retstart += len(id_list)

        # Stop if we've exhausted available results
        total_count = int(result.get("count", 0))
        if retstart >= total_count:
            break

        time.sleep(RATE_LIMIT_DELAY)

    logger.info("PubMed search returned %d PMIDs", len(all_pmids))
    return all_pmids[:max_results]


def fetch_articles(
    pmids: List[str],
    api_key: Optional[str] = None,
) -> List[Dict]:
    """
    Fetch article metadata for a list of PMIDs using efetch.

    Parameters
    ----------
    pmids : list of str
        PMID strings to fetch.
    api_key : str, optional
        NCBI API key. Falls back to ``NCBI_API_KEY`` environment variable.

    Returns
    -------
    list of dict
        Article metadata dictionaries with keys: pmid, title, abstract,
        authors, journal, year, pub_date, keywords, mesh_terms.
    """
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    all_articles: List[Dict] = []

    logger.info("Fetching metadata for %d PMIDs", len(pmids))

    for batch_start in range(0, len(pmids), BATCH_SIZE):
        batch = pmids[batch_start: batch_start + BATCH_SIZE]

        params: Dict[str, str] = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
            "rettype": "abstract",
        }
        if api_key:
            params["api_key"] = api_key

        try:
            response = requests.get(EFETCH_URL, params=params, timeout=60)
            response.raise_for_status()
            articles = _parse_article_xml(response.text)
            all_articles.extend(articles)
        except requests.RequestException as exc:
            logger.error("PubMed efetch request failed for batch at %d: %s", batch_start, exc)
        except ET.ParseError as exc:
            logger.error("XML parsing failed for batch at %d: %s", batch_start, exc)

        time.sleep(RATE_LIMIT_DELAY)

    logger.info("Fetched metadata for %d articles", len(all_articles))
    return all_articles


# ===================================================================
# XML Parsing
# ===================================================================

def _parse_article_xml(xml_text: str) -> List[Dict]:
    """
    Parse PubMed XML efetch response into article metadata dicts.

    Parameters
    ----------
    xml_text : str
        Raw XML response from PubMed efetch.

    Returns
    -------
    list of dict
        Parsed article records.
    """
    articles: List[Dict] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        logger.error("Failed to parse PubMed XML: %s", exc)
        return []

    for article_elem in root.findall(".//PubmedArticle"):
        try:
            article = _parse_single_article(article_elem)
            if article:
                articles.append(article)
        except Exception as exc:
            logger.warning("Failed to parse article element: %s", exc)

    return articles


def _parse_single_article(article_elem: ET.Element) -> Optional[Dict]:
    """Parse a single PubmedArticle XML element."""
    medline = article_elem.find("MedlineCitation")
    if medline is None:
        return None

    # PMID
    pmid_elem = medline.find("PMID")
    pmid = pmid_elem.text if pmid_elem is not None else ""

    article_node = medline.find("Article")
    if article_node is None:
        return None

    # Title
    title_elem = article_node.find("ArticleTitle")
    title = title_elem.text if title_elem is not None else ""

    # Abstract
    abstract_parts: List[str] = []
    abstract_node = article_node.find("Abstract")
    if abstract_node is not None:
        for abs_text in abstract_node.findall("AbstractText"):
            label = abs_text.get("Label", "")
            text = abs_text.text or ""
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)
    abstract = " ".join(abstract_parts)

    # Authors
    authors: List[str] = []
    author_list = article_node.find("AuthorList")
    if author_list is not None:
        for author in author_list.findall("Author"):
            last = author.find("LastName")
            first = author.find("ForeName")
            if last is not None and first is not None:
                authors.append(f"{last.text} {first.text}")
            elif last is not None:
                authors.append(last.text)

    # Journal
    journal_node = article_node.find("Journal")
    journal = ""
    year = ""
    pub_date = ""
    if journal_node is not None:
        title_node = journal_node.find("Title")
        if title_node is not None:
            journal = title_node.text or ""

        ji_node = journal_node.find("JournalIssue")
        if ji_node is not None:
            pd_node = ji_node.find("PubDate")
            if pd_node is not None:
                year_elem = pd_node.find("Year")
                month_elem = pd_node.find("Month")
                if year_elem is not None:
                    year = year_elem.text or ""
                    pub_date = year
                    if month_elem is not None:
                        pub_date = f"{year} {month_elem.text}"

    # Keywords
    keywords: List[str] = []
    keyword_list = medline.find("KeywordList")
    if keyword_list is not None:
        for kw in keyword_list.findall("Keyword"):
            if kw.text:
                keywords.append(kw.text)

    # MeSH terms
    mesh_terms: List[str] = []
    mesh_list = medline.find("MeshHeadingList")
    if mesh_list is not None:
        for mesh_heading in mesh_list.findall("MeshHeading"):
            descriptor = mesh_heading.find("DescriptorName")
            if descriptor is not None and descriptor.text:
                mesh_terms.append(descriptor.text)

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "authors": authors,
        "journal": journal,
        "year": year,
        "pub_date": pub_date,
        "keywords": keywords,
        "mesh_terms": mesh_terms,
    }
