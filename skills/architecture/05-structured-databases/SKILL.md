---
name: 05-structured-databases
description: >-
  Best-practice standards for Pillar 05 (Structured Databases) of the HCLS AI Factory. Use when designing,
  building, operating, or reviewing tabular/relational state — variant stores, service databases, and
  platform metadata. Triggers: adding a SQL-backed store, a schema migration, a region/QC query, MLOps
  run tracking, or any "which database do I reach for?" decision on the single box.
---

# Pillar 05 — Structured Databases

The tabular/relational tier of the factory: columnar analytics over variants, small relational stores
per service, and the platform's own metadata (MLOps runs, model registry, lineage, governance state).
On a single DGX Spark there is no data warehouse — you pick the right embedded database for the job and
keep state small and local.

## In the HCLS AI Factory
- **DuckDB — the queryable variant store** (`variant-store`, :8575, `serving: native`, tags `duckdb`,
  `E1/E2`). VCF → columnar store; SQL/region queries by `chrom/pos/filter/sample`, PASS-vs-RefCall
  aware, with Ts/Tv QC (verified 2.40 over PASS SNVs on real HG002). `POST /query /load`, `GET /stats`.
  GWAS association testing (`gwas-association`) and the ACMG SF panel (`acmg-secondary-findings`) read
  the same columnar variant layer at single-box scale — no Spark.
- **PostgreSQL + pgvector** where a service genuinely needs a persistent relational + vector store on
  one engine (e.g. clinical-imaging, the tuberous-sclerosis disease program's orchestrator state) —
  not a platform-wide default.
- **SQLite — the platform metadata spine.** Single-box MLOps (`local-mlops`, `mlops.py`) stores runs,
  params, metrics, the model registry (`model_versions` with `none/staging/production/archived`
  stages), and run lineage in one file (`data/mlops.db`), stdlib-only. Governance/repro state is
  likewise small and local. `hcls_common.get_mlops_store()`.
- Registry facts: MLOps and governance capabilities are `serving: none` (library/in-process, no port);
  the variant store is a real `native` FastAPI service.

## Best-practice standards
- **Right database for the job at single-box scale.** Columnar analytics (variant scans, GWAS,
  region/QC aggregates) → DuckDB. Small relational service state → SQLite or PostgreSQL. Platform
  metadata (runs/registry/lineage) → the shared SQLite MLOps store — don't stand up a second one.
- **Embedded first.** Prefer in-process/embedded engines (DuckDB, SQLite) over a server process; only
  run PostgreSQL when a service needs concurrent writers or pgvector co-located with rows.
- **Keep state small and local.** The 700 GB+ working tree stays on disk, not in git; databases hold
  indexes, metadata, and query-ready columns — not raw FASTQ/BAM/weights.
- **Migrations are versioned and forward-only.** Every schema change ships as an idempotent migration
  (`CREATE TABLE IF NOT EXISTS`, additive columns) in the same commit as the code that needs it; never
  hand-edit a live schema.
- **Query by contract, parameterized.** Region queries take typed `{chrom, pos, filter, sample}`;
  bind parameters — never string-format SQL from request input.
- **QC is a first-class column path.** Expose the QC signals the science needs (Ts/Tv, PASS vs RefCall,
  MAF filters) as queryable outputs, not buried print statements.
- **One writer per file.** DuckDB/SQLite files are single-writer; serialize writes behind the owning
  service and let readers open read-only.
- **Back up what you can't recompute.** MLOps/lineage/governance `.db` files are durable state — put
  them on a named, backed-up volume; a recomputable variant store can be rebuilt from the VCF.

## Do / Don't
**Do:** reach for DuckDB for variant/region/QC analytics; keep MLOps in the one shared SQLite store;
ship schema changes as idempotent migrations; parameterize every query; register any SQL-backed
service in the capability registry with its real `serving`/port.
**Don't:** stand up a data warehouse or a second metadata database; hardcode `/home/<user>` DB paths
(derive from repo root or `$HCLS_ROOT`); commit `.db`/`.duckdb` files or query outputs to git;
interpolate untrusted input into SQL; open the same embedded file from two concurrent writers.

## Wiring it in
- Variant store service lives behind `variant-store` (:8575); load then query:
  ```bash
  curl -s localhost:8575/stats
  curl -s -X POST localhost:8575/query -d '{"region":{"chrom":"chr17","start":41196312,"end":41277500}}'
  ```
- Platform metadata via the shared store — don't reinvent tracking:
  ```python
  from hcls_common import get_mlops_store
  store = get_mlops_store()                       # data/mlops.db (or ":memory:" in tests)
  run = store.start_run(name="admet-scan", capability="chemprop-admet")
  store.log_metric(run, "bbb", 0.98); store.end_run(run, "complete")
  ```
- Register any new SQL-backed service in `lib/hcls_common/capabilities.json` (`type: service`/`model`,
  real `endpoint`, `serving: native|container`, `status: live`) and map its directory in
  `scripts/validate_registry.py` → `python scripts/validate_registry.py` must print `OK`.
- Config (DB paths, ports) comes from `hcls_common/config.py` / env vars, never literals.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **No warehouse, no distributed SQL.** Spark/Snowflake-shaped designs don't fit one box — keep it
  embedded. GWAS is pure-Python `statsmodels`, deliberately not Spark.
- **Unified 128 GB memory is shared** with GPU inference; a DuckDB scan that materializes a huge
  result can starve ESMFold/DiffDock. Stream/aggregate in SQL, don't `SELECT *` a whole chromosome.
- **ARM wheels.** Prefer engines with clean aarch64 wheels (DuckDB, SQLite ship them); a psycopg/Postgres
  client stack that assumes x86 build tools can bite — pin and test on the GB10.
- **Single-writer files corrupt under concurrency** — a crashed writer can leave a `.duckdb`/`.db`
  mid-write; back up durable metadata and prefer WAL where available.
- **Named volumes must be explicitly backed** — an unbacked docker volume holding `mlops.db` is lost
  on `down -v`.

## Related
- Pillars: 06-vector-databases-and-embeddings, 07-message-bus-and-async, 08-inference-serving
- build-housekeeping-standards
