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
| cart | 879 | 9 seeders (regulatory, safety, assays, biomarkers, patents, sequences, …) |
| precision-biomarker | 622 | `seed_all.py` |
| clinical-imaging | 440 | pre-existing |
| neurology | 308 | `neuro_electrophysiology` still empty — see below |
| pharmacogenomics | 240 | |
| rare-disease | 176 | `rd_diseases` fails on a required `disease_id` |
| single-cell | 164 | `sc_markers` still empty — see below |
| clinical-trial | 59 | |

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

## Known gaps — parser/schema drift, not bugs

`neuro_electrophysiology` and `sc_markers` genuinely share no fields with their seed records: the
parsers and the collection schemas were designed independently. The seeder now reports this
explicitly ("No rows survived schema projection") rather than inserting empty rows. Reconciling
them means deciding which side is canonical — a modelling decision, not a fix.

`rd_diseases` needs a required `disease_id` the parser does not emit.

**Where only a column or two overlaps**, the seeder writes the record's source text into whatever
prose column the schema offers (`description`, `text`, `summary`, `abstract`, …). Without that the
rows embed correctly — the vector is built from the real text — but return nothing readable.

## Verify

```bash
.venv/bin/python scripts/run_demo.py A1   # cart
.venv/bin/python scripts/run_demo.py A6   # clinical-trial
```

Both assert on retrieved content, so an empty corpus fails them rather than passing quietly.
