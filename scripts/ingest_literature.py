#!/usr/bin/env python3
"""Ingest real PubMed literature into any subject's Milvus collection.

Three subjects shipped their own PubMed ingester and eight did not, which is most of why 44
collections are empty and eight subjects hold fewer than 500 vectors. The literature is public
and free — what was missing was a fetcher not welded to one agent's package layout, and a way to
write it into collections whose schemas all differ.

Both pieces already existed separately:

  hcls_common.pubmed          fetch + parse, with NCBI rate limits and 5xx retry
  hcls_common.ingest_persist  project a record onto whatever a collection declares, then upsert

This joins them. Every literature collection in this platform has a different schema
(`text` vs `text_chunk` vs `abstract_text`, `source` vs `source_type`, …), and dynamic fields
are off, so the projection step is not optional — an unexpected key aborts the whole batch.

    .venv/bin/python scripts/ingest_literature.py --collection rd_literature \\
        --query '"tuberous sclerosis"[Title/Abstract]' --max 200
    .venv/bin/python scripts/ingest_literature.py --collection rd_literature --query '…' --dry-run

Set NCBI_API_KEY to lift the rate limit from 3/s to 10/s.
"""
from __future__ import annotations

import argparse
import logging
import os
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("ingest_literature")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--collection", required=True, help="target Milvus collection")
    ap.add_argument("--query", required=True, help="PubMed query (E-utilities syntax)")
    ap.add_argument("--max", type=int, default=200, help="max articles (default 200)")
    ap.add_argument("--dry-run", action="store_true", help="fetch and show, write nothing")
    ap.add_argument("--min-abstract", type=int, default=200,
                    help="skip records with a shorter abstract (default 200 chars)")
    a = ap.parse_args()

    from hcls_common.pubmed import PubMed

    pm = PubMed()
    logger.info("searching PubMed: %s", a.query)
    pmids = pm.search(a.query, max_results=a.max)
    logger.info("%d PMIDs", len(pmids))
    if not pmids:
        return 0

    arts = pm.fetch(pmids)
    logger.info("%d articles fetched", len(arts))

    # An abstract-less record embeds to little more than its title and pollutes retrieval with
    # near-empty passages that still score. Skip rather than store.
    keep = [x for x in arts if len(x.abstract) >= a.min_abstract]
    skipped = len(arts) - len(keep)
    if skipped:
        logger.info("%d skipped (abstract under %d chars)", skipped, a.min_abstract)
    if not keep:
        return 0

    if a.dry_run:
        print(f"\n  would write {len(keep)} records to {a.collection}:\n")
        for x in keep[:5]:
            print(f"    {x.pmid} ({x.year}) {x.journal[:44]}")
            print(f"      {x.title[:96]}")
        print(f"    … and {max(0, len(keep) - 5)} more")
        return 0

    from hcls_common.ingest_persist import persist_records

    host = os.getenv("MILVUS_HOST", "localhost")
    port = os.getenv("MILVUS_PORT", "19530")
    written = persist_records(
        [x.as_record() for x in keep], a.collection,
        milvus_uri=f"http://{host}:{port}",
    )
    logger.info("wrote %d rows to %s", written, a.collection)
    # persist_records never raises; a zero here means the projection or the client failed, and
    # reporting success on zero rows is the exact failure this platform keeps finding.
    if written == 0:
        logger.error("nothing was written — check the collection schema and the Milvus connection")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
