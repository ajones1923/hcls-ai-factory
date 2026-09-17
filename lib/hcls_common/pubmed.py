"""One PubMed client for the whole factory.

Three subjects ship their own PubMed ingester (cart, precision-oncology, clinical-imaging) and
eight do not — which is most of why 44 Milvus collections are empty and eight subjects hold fewer
than 500 vectors. The literature is public and free; the only thing missing was a way to fetch it
that was not welded to one agent's package layout.

Deliberately small and dependency-free (stdlib urllib, no requests): this has to run from any
subject's venv without adding to its install.

NCBI etiquette is enforced here rather than left to each caller — 3 requests/second without an
API key, 10 with one (`NCBI_API_KEY`). Exceeding it gets an IP blocked, which would break every
ingester at once.
"""
from __future__ import annotations

import json
import logging
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_UA = "hcls-ai-factory (https://github.com/ajones1923/hcls-ai-factory)"


@dataclass
class Article:
    pmid: str
    title: str = ""
    abstract: str = ""
    journal: str = ""
    year: str = ""
    doi: str = ""
    authors: list[str] = field(default_factory=list)

    @property
    def text(self) -> str:
        """What gets embedded: title plus abstract, which is how a reader would search it."""
        return f"{self.title}\n\n{self.abstract}".strip()

    def as_record(self, **extra) -> dict:
        """The shape `hcls_common.ingest_persist.persist_records` expects."""
        meta = {
            "pmid": self.pmid, "title": self.title, "journal": self.journal,
            "year": self.year, "doi": self.doi,
            "authors": ", ".join(self.authors[:8]),
            "source": "pubmed", "source_type": "literature",
            "abstract_summary": self.abstract[:900],
            "abstract_text": self.abstract,
            "text_chunk": self.text,
        }
        meta.update(extra)
        return {"text": self.text, "metadata": meta}


class PubMed:
    def __init__(self, api_key: str | None = None, timeout: int = 30):
        self.api_key = api_key or os.getenv("NCBI_API_KEY") or ""
        self.timeout = timeout
        # 10/s with a key, 3/s without. Stay just inside it.
        self._gap = 0.11 if self.api_key else 0.35
        self._last = 0.0

    def _wait(self) -> None:
        delta = time.monotonic() - self._last
        if delta < self._gap:
            time.sleep(self._gap - delta)
        self._last = time.monotonic()

    def _get(self, path: str, params: dict) -> bytes:
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{EUTILS}/{path}?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"User-Agent": _UA})
        # NCBI returns 502/503 under load — seen on the first live call from this box. A single
        # transient 5xx must not abort an ingest that is 40 batches deep, so back off and retry.
        last: Exception | None = None
        for attempt in range(4):
            self._wait()
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    return r.read()
            except urllib.error.HTTPError as exc:
                last = exc
                if exc.code not in (429, 500, 502, 503, 504):
                    raise
            except Exception as exc:                 # transport-level blips get the same treatment
                last = exc
            backoff = 1.5 * (2 ** attempt)
            logger.warning("NCBI %s failed (%s); retry %d/3 in %.1fs", path, last, attempt + 1, backoff)
            time.sleep(backoff)
        raise RuntimeError(f"NCBI {path} failed after 4 attempts: {last}")

    def search(self, query: str, max_results: int = 200) -> list[str]:
        """PMIDs for a query, newest first."""
        out: list[str] = []
        page = 200                                   # NCBI's comfortable retmax
        while len(out) < max_results:
            want = min(page, max_results - len(out))
            raw = self._get("esearch.fcgi", {
                "db": "pubmed", "term": query, "retmax": want, "retstart": len(out),
                "retmode": "json", "sort": "date",
            })
            ids = json.loads(raw).get("esearchresult", {}).get("idlist", [])
            if not ids:
                break
            out.extend(ids)
            if len(ids) < want:
                break
        return out[:max_results]

    def fetch(self, pmids: list[str], batch: int = 100) -> list[Article]:
        """Full records for PMIDs. Skips anything unparseable rather than aborting the batch."""
        arts: list[Article] = []
        for i in range(0, len(pmids), batch):
            chunk = pmids[i:i + batch]
            try:
                raw = self._get("efetch.fcgi", {
                    "db": "pubmed", "id": ",".join(chunk), "retmode": "xml",
                })
                arts.extend(_parse(raw))
            except Exception as exc:
                logger.warning("efetch failed for %d PMIDs (%s) — skipping that batch",
                               len(chunk), exc)
        return arts


def _text(node, path: str, default: str = "") -> str:
    el = node.find(path)
    return "".join(el.itertext()).strip() if el is not None else default


def _parse(raw: bytes) -> list[Article]:
    arts: list[Article] = []
    root = ET.fromstring(raw)
    for art in root.findall(".//PubmedArticle"):
        pmid = _text(art, ".//PMID")
        if not pmid:
            continue
        # An abstract can be split into labelled sections; join them in document order.
        abstract = " ".join(
            "".join(seg.itertext()).strip()
            for seg in art.findall(".//Abstract/AbstractText")
        ).strip()
        authors = []
        for a in art.findall(".//Author"):
            last, init = _text(a, "LastName"), _text(a, "Initials")
            if last:
                authors.append(f"{last} {init}".strip())
        doi = ""
        for aid in art.findall(".//ArticleId"):
            if aid.get("IdType") == "doi":
                doi = (aid.text or "").strip()
        arts.append(Article(
            pmid=pmid,
            title=_text(art, ".//ArticleTitle"),
            abstract=abstract,
            journal=_text(art, ".//Journal/Title"),
            year=_text(art, ".//PubDate/Year") or _text(art, ".//PubDate/MedlineDate")[:4],
            doi=doi,
            authors=authors,
        ))
    return arts
