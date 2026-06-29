"""PubMed PGx literature ingest pipeline for Pharmacogenomics Intelligence Agent.

Fetches pharmacogenomics research papers via NCBI E-utilities (esearch +
efetch), parses PubMed XML into ClinicalEvidence models, and stores
embeddings in the pgx_clinical_evidence Milvus collection.

Follows the same pattern as the CAR-T literature_parser but with
PGx-specific search terms, study type classification, and gene/drug
extraction.

Author: Adam Jones
Date: March 2026
"""

import re
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import PGxCollectionManager
from src.models import ClinicalEvidence

from .base import BaseIngestPipeline, _truncate_utf8


# NCBI E-utilities endpoints
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{EUTILS_BASE_URL}/esearch.fcgi"
EFETCH_URL = f"{EUTILS_BASE_URL}/efetch.fcgi"

# Default PGx PubMed query
DEFAULT_PGX_QUERY = (
    '"pharmacogenomics" OR "pharmacogenetics" OR "pharmacogenomic" '
    'OR "CYP2D6" OR "CYP2C19" OR "CYP2C9" OR "DPYD" OR "TPMT" '
    'OR "genotype-guided dosing" OR "genotype-guided therapy" '
    'OR "precision prescribing"'
)

# Core PGx genes for extraction from abstracts
_PGX_GENE_PATTERNS = [
    r"\bCYP2D6\b", r"\bCYP2C19\b", r"\bCYP2C9\b", r"\bCYP3A5\b",
    r"\bCYP3A4\b", r"\bCYP2B6\b", r"\bCYP1A2\b", r"\bCYP2A6\b",
    r"\bCYP4F2\b", r"\bDPYD\b", r"\bTPMT\b", r"\bNUDT15\b",
    r"\bUGT1A1\b", r"\bSLCO1B1\b", r"\bABCG2\b", r"\bVKORC1\b",
    r"\bIFNL3\b", r"\bNAT2\b", r"\bG6PD\b", r"\bHLA-[AB]\b",
    r"\bABCB1\b", r"\bCOMT\b", r"\bOPRM1\b",
]
_PGX_GENE_REGEX = re.compile("|".join(_PGX_GENE_PATTERNS), re.IGNORECASE)

# High-interest PGx drug patterns
_PGX_DRUG_PATTERNS = [
    r"\bwarfarin\b", r"\bclopidogrel\b", r"\btamoxifen\b", r"\bcodeine\b",
    r"\btramadol\b", r"\bsimvastatin\b", r"\batorvastatin\b", r"\babacavir\b",
    r"\bcarbamazepine\b", r"\bphenytoin\b", r"\bfluorouracil\b", r"\b5-FU\b",
    r"\bcapecitabine\b", r"\bmercaptopurine\b", r"\bazathioprine\b",
    r"\btacrolimus\b", r"\bvoriconazole\b", r"\birinotecan\b",
    r"\ballopurinol\b", r"\bomeprazole\b", r"\besomeprazole\b",
    r"\bsertraline\b", r"\bescitalopram\b", r"\bamitriptyline\b",
    r"\bnortriptyline\b", r"\batomoxetine\b", r"\bondansetron\b",
]
_PGX_DRUG_REGEX = re.compile("|".join(_PGX_DRUG_PATTERNS), re.IGNORECASE)

# Study type classification keywords
_STUDY_TYPE_KEYWORDS = {
    "meta-analysis": ["meta-analysis", "meta analysis", "systematic review"],
    "RCT": ["randomized", "randomised", "rct", "controlled trial"],
    "cohort": ["cohort", "prospective", "retrospective", "observational"],
    "case-control": ["case-control", "case control"],
    "implementation": ["implementation", "clinical implementation", "pgx program"],
    "review": ["review", "overview", "narrative review"],
}


class PubMedPGxParser(BaseIngestPipeline):
    """Ingest pipeline for PubMed PGx literature.

    Fetches pharmacogenomics abstracts from PubMed using NCBI E-utilities,
    extracts gene and drug mentions, classifies study type, and stores
    results in the pgx_clinical_evidence Milvus collection.

    Usage:
        parser = PubMedPGxParser(collection_manager, embedder)
        count = parser.run(query="CYP2D6 codeine", max_results=1000)
    """

    COLLECTION_NAME = "pgx_clinical_evidence"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
    ):
        super().__init__(collection_manager, embedder)
        self.api_key = api_key
        self.email = email or ""
        self._min_interval = 1.0 / 10 if api_key else 1.0 / 3
        self._last_request_time: float = 0.0

    def _rate_limit(self) -> None:
        """Enforce NCBI rate limiting between requests."""
        now = time.time()
        elapsed = now - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

    def _build_params(self) -> Dict[str, str]:
        """Build base query parameters for E-utilities requests."""
        params: Dict[str, str] = {"tool": "pgx_intelligence_agent"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def fetch(
        self,
        query: str = DEFAULT_PGX_QUERY,
        max_results: int = 5000,
        **kwargs,
    ) -> List[Dict[str, Any]]:
        """Fetch PGx abstracts from PubMed via NCBI E-utilities.

        Performs esearch to get PMIDs, then efetch to retrieve abstracts.

        Returns:
            List of article dicts with pmid, title, abstract, authors,
            journal, year, mesh_terms.
        """
        # Step 1: esearch to get PMIDs
        pmids: List[str] = []
        try:
            params = self._build_params()
            params.update({
                "db": "pubmed", "term": query,
                "retmode": "json", "retmax": str(min(max_results, 10000)),
            })
            self._rate_limit()
            resp = requests.get(ESEARCH_URL, params=params, timeout=30)
            resp.raise_for_status()
            result = resp.json()
            pmids = result["esearchresult"]["idlist"]
            total = int(result["esearchresult"]["count"])
            logger.info(f"PubMed PGx search found {total} results, retrieving {len(pmids)}")
        except (requests.RequestException, KeyError, ValueError) as e:
            logger.error(f"PubMed esearch failed: {e}")
            return []

        if not pmids:
            return []

        # Step 2: efetch abstracts in batches
        articles: List[Dict[str, Any]] = []
        batch_size = 200

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            try:
                params = self._build_params()
                params.update({
                    "db": "pubmed", "id": ",".join(batch),
                    "rettype": "xml", "retmode": "xml",
                })
                self._rate_limit()
                resp = requests.get(EFETCH_URL, params=params, timeout=60)
                resp.raise_for_status()
                root = ET.fromstring(resp.content)

                for article_elem in root.findall(".//PubmedArticle"):
                    article = self._parse_xml_article(article_elem)
                    if article:
                        articles.append(article)

                logger.info(f"Fetched batch {i // batch_size + 1} ({len(articles)} articles so far)")

            except requests.RequestException as e:
                logger.error(f"efetch HTTP error for batch starting at {i}: {e}")
            except ET.ParseError as e:
                logger.error(f"XML parse error for batch starting at {i}: {e}")

        logger.info(f"Fetched {len(articles)} PGx articles from PubMed")
        return articles

    @staticmethod
    def _parse_xml_article(article_elem: ET.Element) -> Optional[Dict[str, Any]]:
        """Parse a single PubmedArticle XML element into a dict."""
        pmid_elem = article_elem.find(".//PMID")
        pmid = pmid_elem.text if pmid_elem is not None and pmid_elem.text else ""

        title_elem = article_elem.find(".//ArticleTitle")
        title = title_elem.text if title_elem is not None and title_elem.text else ""

        abstract_parts = []
        for abs_elem in article_elem.findall(".//AbstractText"):
            if abs_elem.text:
                label = abs_elem.get("Label")
                if label:
                    abstract_parts.append(f"{label}: {abs_elem.text}")
                else:
                    abstract_parts.append(abs_elem.text)
        abstract = " ".join(abstract_parts)

        journal_elem = article_elem.find(".//Journal/Title")
        journal = journal_elem.text if journal_elem is not None and journal_elem.text else ""

        year_elem = article_elem.find(".//PubDate/Year")
        if year_elem is not None and year_elem.text:
            year = year_elem.text
        else:
            medline_elem = article_elem.find(".//PubDate/MedlineDate")
            year = medline_elem.text[:4] if medline_elem is not None and medline_elem.text and len(medline_elem.text) >= 4 else ""

        mesh_terms = [
            m.text for m in article_elem.findall(".//MeshHeading/DescriptorName")
            if m.text
        ]

        return {
            "pmid": pmid, "title": title, "abstract": abstract,
            "journal": journal, "year": year, "mesh_terms": mesh_terms,
        }

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[ClinicalEvidence]:
        """Parse PubMed article dicts into ClinicalEvidence models."""
        records: List[ClinicalEvidence] = []

        for article in raw_data:
            try:
                pmid = article.get("pmid", "")
                title = article.get("title", "")
                abstract = article.get("abstract", "")
                full_text = f"{title} {abstract}".strip()

                # Extract gene and drug mentions
                gene = self._extract_gene(full_text)
                drug = self._extract_drug(full_text)

                # Classify study type
                study_type = self._classify_study_type(full_text)

                # Parse year
                try:
                    year = int(article.get("year", 0))
                except (TypeError, ValueError):
                    year = 0

                record = ClinicalEvidence(
                    id=_truncate_utf8(f"pubmed_{pmid}", 95),
                    title=_truncate_utf8(title, 490),
                    text_chunk=_truncate_utf8(full_text, 2990),
                    study_type=_truncate_utf8(study_type, 95),
                    gene=_truncate_utf8(gene, 48),
                    drug=_truncate_utf8(drug, 195),
                    phenotype="",
                    outcome_measure="",
                    outcome_value="",
                    sample_size=0,
                    pmid=_truncate_utf8(pmid, 18),
                    year=year,
                    source="PubMed",
                )
                records.append(record)

            except Exception as e:
                logger.warning(f"Failed to parse article {article.get('pmid', '?')}: {e}")
                continue

        logger.info(f"Parsed {len(records)} PGx literature records")
        return records

    @staticmethod
    def _extract_gene(text: str) -> str:
        """Extract the first PGx gene mentioned in text."""
        if not text:
            return ""
        match = _PGX_GENE_REGEX.search(text)
        return match.group(0).upper() if match else ""

    @staticmethod
    def _extract_drug(text: str) -> str:
        """Extract the first PGx-relevant drug mentioned in text."""
        if not text:
            return ""
        match = _PGX_DRUG_REGEX.search(text)
        return match.group(0).lower() if match else ""

    @staticmethod
    def _classify_study_type(text: str) -> str:
        """Classify study type from abstract text using keyword matching."""
        if not text:
            return ""
        text_lower = text.lower()
        for study_type, keywords in _STUDY_TYPE_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return study_type
        return "primary research"
