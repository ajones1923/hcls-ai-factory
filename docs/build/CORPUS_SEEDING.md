# Seeding the agent corpora

**Status:** operational note, repo-only. **Written:** 2026-09-15, taking the fleet from 440 vectors
(imaging only) to 3,570 across eight agents.

An agent with an empty corpus answers HTTP 200, passes its whole test suite, and returns nothing.
That is the failure mode this page exists to prevent recurring.

## The order that works

```bash
cd core/agents/<agent>
./venv/bin/python scripts/setup_collections.py   # create + index + load
./venv/bin/python scripts/seed_knowledge.py      # embed + insert
# then RESTART the service — engines cache collection handles at startup
```

Three steps, and **all three are required**. Seeding before the collections exist inserts nothing;
seeding without restarting leaves the running engine holding handles from before.

**Milvus collections must be `load()`ed into memory before search.** Nothing did this, so even a
correctly seeded collection answered `collection not loaded`. `setup_collections.py` now loads on
create; for pre-existing collections:

```python
from pymilvus import MilvusClient
c = MilvusClient(uri="http://localhost:19530")
for n in c.list_collections():
    c.load_collection(collection_name=n)
```

## Current state

| Agent | Vectors | Note |
|---|---:|---|
| precision-biomarker | 1,244 | `seed_all.py` |
| cart | 879 | 9 seeders (regulatory, safety, assays, biomarkers, patents, sequences, …) |
| clinical-imaging | 440 | pre-existing |
| rare-disease | 440 | seed + OMIM/HPO/Orphanet/gene-therapy ingest |
| single-cell | 279 | seed + ingest |
| pharmacogenomics | 240 | its pipelines always persisted |
| clinical-trial | 178 | seed + live ClinicalTrials.gov ingest |
| neurology | 168 | seed + PubMed ingest |
| precision-autoimmune | 116 | had **no seeder at all** until 2026-09-15 |

**~4,000 vectors across 8 agents**, up from 440 (imaging only).

## Five traps, each of which fails silently

1. **A setup script that only logs.** neurology's and single-cell's `setup_collections.py` printed
   `[create] <name>` and `Collection setup complete` while creating nothing. If a seeder reports
   "can't find collection" right after a successful-looking setup, this is why.
2. **Fields live inside `IngestRecord.metadata`.** Parsers emit
   `IngestRecord(text=…, metadata={pmid, title, …})`; the schemas declare `pmid`/`title` as
   top-level columns. Reading the record's own attributes yields
   `{text, metadata, collection_name, record_id, source}` — which shares no column with the
   schema, so nothing inserts.
3. **`str(DataType.FLOAT)` is `'10'`, not `'FLOAT'`.** Detecting numeric columns by string match
   classifies every float as text. Use `.name`. And keep int64 and float apart — coercing both to
   float fails an int64 insert.
4. **Dynamic fields are off and nothing is nullable.** One unexpected key aborts the entire batch;
   one missing declared column does too. Project each row onto the schema and fill the remainder.
5. **Display labels vs enum slugs.** CAR-T's regulatory corpus says
   `"Breakthrough Therapy Designation"`; the enum holds `breakthrough_therapy`. All 40 records
   failed validation and were logged away, leaving the collection empty.

## Ingest: fetch is not persist

Several `run_ingest.py` scripts fetched, parsed, validated — and wrote nothing. They logged
"N records validated" and left the corpus untouched, which is why agents sat at a few dozen rows
while their ingest "worked". Fixed in `lib/hcls_common/ingest_persist.py`, **one shared
implementation** used by clinical-trial, neurology, single-cell and rare-disease.

`pharmacogenomics` was never broken — its pipelines always persisted. Absence of the function
name is not absence of the behaviour; check what a script *does*, not what it is called.

```bash
cd core/agents/<agent>
./venv/bin/python scripts/run_ingest.py --source all            # fetch + persist
./venv/bin/python scripts/run_ingest.py --source all --dry-run  # the old behaviour
```

### Why the projection rules live in one file

Each rule below came from a failed insert. Milvus aborts the **whole batch** on a single bad key
and the caller logs a warning and moves on, so a per-agent copy that misses one fails silently:

- dynamic fields are OFF and nothing is nullable — every declared column must be present, and no
  undeclared key may appear
- `str(DataType.FLOAT)` is its numeric **code** (`"10"`), not `"FLOAT"` — use `.name`
- int64 and float are different; coercing both to float fails an int64 insert
- only genuine ARRAY columns may take a list; a list in a VARCHAR column aborts the batch
- **BOOL** is neither numeric nor text — defaulting it to `""` fails the insert
- the primary key differs per collection (VARCHAR here, int64 there) and is not `auto_id`
- a content-derived id makes re-ingest idempotent via upsert instead of duplicating the corpus

Two collections previously recorded here as unresolvable parser/schema drift — `sc_markers` and
`rd_diseases` — load correctly through this projection. The drift was real; it was a projection
problem, not a modelling one.

**What is still dropped:** PubMed records routed to a domain collection (`neuro_oncology`,
`neuro_headache`, ...) lose `pmid`/`title`, because those schemas do not declare them. The title
survives inside the embedded text, so retrieval works, but the structured citation id does not.
The loader reports exactly which fields it ignored, once per collection.

## Verify

```bash
.venv/bin/python scripts/run_demo.py A1   # cart
.venv/bin/python scripts/run_demo.py A6   # clinical-trial
```

Both assert on retrieved content, so an empty corpus fails them rather than passing quietly.
