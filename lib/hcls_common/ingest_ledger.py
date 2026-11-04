"""An append-only record of what was ingested into each collection, and with which model.

Why this exists rather than provenance columns on the records themselves: every collection in
this platform declares a fixed schema with **dynamic fields off and nothing nullable**
(`ingest_persist`), so an undeclared key aborts the whole batch. Adding `source_query`,
`embed_model` and `retrieved_at` to the records would mean altering 114 schemas. The ledger
gets the same auditability without touching any of them.

What it is for
--------------
**Provenance.** This repo's thesis is that a clinical claim should be traceable to a primary
source. A corpus without provenance is the same defect one layer down: a retrieved passage that
influences an answer, with no record of where it came from, which query selected it, or when.

**Catching a model change, which is the silent one.** Milvus validates vector dimension on
insert, so embedding with a 768-dim model into a 384-dim collection fails loudly. Embedding with
a *different 384-dim model* does not fail at all — the insert succeeds, the vectors are
meaningless relative to their neighbours, and retrieval quietly returns worse matches with no
error anywhere. Dimension checking cannot see that. A recorded model per collection can.

The ledger holds metadata only — collection, source, query, counts, model, timestamp. No
abstracts, because PubMed abstracts stay on the box under NCBI's terms while code and docs
publish.
"""
from __future__ import annotations

import datetime
import json
import pathlib
from dataclasses import asdict, dataclass, field
from typing import Iterator

#: Default location. Committed, on the same reasoning as docs/build/corpus_snapshot.json:
#: it is a measurement of the box that a reviewer needs in order to judge a change.
DEFAULT_LEDGER = "docs/build/corpus_provenance.jsonl"


@dataclass(frozen=True)
class IngestEvent:
    """One ingest run against one collection."""

    collection: str
    source: str
    """Where the records came from: "pubmed", "clinicaltrials.gov", "clinvar", ..."""
    query: str
    embed_model: str
    fetched: int
    written: int
    skipped: int = 0
    at: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc)
                    .replace(microsecond=0).isoformat())
    note: str = ""

    def as_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def read(path: str | pathlib.Path = DEFAULT_LEDGER) -> list[IngestEvent]:
    """Every recorded event, oldest first. A missing ledger is empty, not an error."""
    p = pathlib.Path(path)
    if not p.exists():
        return []
    out: list[IngestEvent] = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            # One corrupt line must not hide the rest of the history.
            continue
        known = {k: d[k] for k in d if k in IngestEvent.__dataclass_fields__}
        try:
            out.append(IngestEvent(**known))
        except TypeError:
            continue
    return out


def append(event: IngestEvent, path: str | pathlib.Path = DEFAULT_LEDGER) -> None:
    """Append one event. Creates the file and its parent if needed."""
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a") as fh:
        fh.write(event.as_json() + "\n")


def models_used(collection: str,
                path: str | pathlib.Path = DEFAULT_LEDGER) -> list[str]:
    """Distinct embedding models previously written into `collection`, oldest first."""
    seen: list[str] = []
    for e in read(path):
        if e.collection == collection and e.written > 0 and e.embed_model not in seen:
            seen.append(e.embed_model)
    return seen


def model_conflict(collection: str, embed_model: str,
                   path: str | pathlib.Path = DEFAULT_LEDGER) -> str | None:
    """Return an explanation if writing `embed_model` would mix models, else None.

    Mixing is not an error the database will raise. Two models of the same dimension insert
    equally happily, and the collection then holds vectors from two different spaces: a
    nearest-neighbour search across them returns confident nonsense. This is the check that
    would notice.
    """
    prior = models_used(collection, path)
    if not prior or embed_model in prior:
        return None
    return (f"{collection} already holds vectors embedded with {', '.join(prior)}; "
            f"writing {embed_model} would mix embedding spaces in one collection. "
            f"Nearest-neighbour search across two spaces returns confident nonsense and "
            f"nothing raises. Re-embed the collection, or write to a new one.")


def history(path: str | pathlib.Path = DEFAULT_LEDGER) -> Iterator[tuple[str, int, str]]:
    """(collection, rows written, model) per collection, for an operator summary."""
    totals: dict[str, tuple[int, str]] = {}
    for e in read(path):
        rows, model = totals.get(e.collection, (0, e.embed_model))
        totals[e.collection] = (rows + e.written, e.embed_model)
    for coll in sorted(totals):
        rows, model = totals[coll]
        yield coll, rows, model
