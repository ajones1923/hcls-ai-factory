---
name: 02-storage-and-data-layer
description: >-
  Best-practice standards for Pillar 02 (Storage & the Data Layer) of the HCLS AI Factory. Use when
  designing, building, operating, or reviewing where bytes live — reference data, model weights,
  variant stores, object storage, and volumes. Concrete triggers: adding a dataset/model, touching
  `.gitignore`, sizing a Docker volume, shipping a fixture, or deciding what may enter git.
---

# Pillar 02 — Storage & the Data Layer

How the factory holds and serves data: a large local working tree (reference genomes, ClinVar,
AlphaMissense, PDB, model weights) that stays on the box, a columnar variant store, object storage
backing the vector database, and tiny fixtures that ship so anyone can reproduce a run.

## In the HCLS AI Factory
- **~711 GB working tree, code-and-docs-only in git.** `data/` and `hcls-ai-factory-core-data/` plus
  all weights/artifacts are gitignored and **never** published; the tree lives on the DGX Spark.
- **Reference data package** (local, out of git): reference genome + indices, ClinVar (~2.7M records),
  AlphaMissense (~71M predictions), PDB structures — downloaded/recreated, not tracked.
- **Columnar variant store (DuckDB):** VCF → DuckDB with Ts/Tv QC (Engine 1, "Variant store", :8575) —
  a queryable analytical layer over called variants, an embedded file database (no server).
- **Object storage backs the vector database:** `milvus-minio` (S3-compatible object storage) holds
  Milvus segment/blob data; `milvus-etcd` holds metadata — see `docker-compose.dgx-spark.yml`. Named
  volumes `minio_data`, `milvus_data`, `etcd_data` persist them.
- **Fixtures ship, data doesn't:** a tiny demo VCF (~28 KB) is tracked for reproducibility; the largest
  tracked file stays well under 1 MB. `.gitignore` is the single source of truth guarding the boundary.

## Best-practice standards
- **`.gitignore` is the guard — extend it before adding a big format.** New weights/data extensions get
  a rule in the same change; never `git add -f` a `data/`, `outputs/`, `results/`, or `*.bak` path.
- **Reference data is recreated, not committed.** Ship a downloader/manifest (URL + checksum + version),
  not the bytes; record the version so a run is reproducible.
- **Ship a tiny, real fixture for every data-consuming capability** (like the ~28 KB demo VCF) so tests
  and demos run without the 711 GB tree.
- **Keep the analytical store embedded and file-based** (DuckDB) — no extra server to operate on the box.
- **Persist stateful services to named volumes** (Milvus/MinIO/etcd, Prometheus, Grafana); don't rely on
  the container filesystem for anything you need after `down`.
- **Largest tracked file well under 1 MB** — the pre-commit oversize guard blocks the rest; if something
  large "needs" to be in git, it belongs in the data package instead.
- **No secrets or machine-specific paths in tracked data files** — derive paths from repo root / env,
  not `/home/<user>/...`.

## Do / Don't
**Do:** gitignore weights + `data/` + `hcls-ai-factory-core-data/`; ship a downloader + checksum for
reference data; commit a <1 MB fixture; back the vector DB with the MinIO object store on a named volume;
query variants from the DuckDB store.
**Don't:** commit BAM/FASTQ/`*.pt`/`*.safetensors`/`*.onnx`; `git add -f` generated outputs; hardcode an
absolute home path to a dataset; put a full ClinVar/AlphaMissense dump in git; run a stateful service on
an anonymous/ephemeral volume.

## Wiring it in
```gitignore
# .gitignore already guards the boundary — extend it, don't work around it
data/
hcls-ai-factory-core-data/
*.bam  *.fastq  *.fastq.gz          # sequence data
*.pt   *.pth    *.safetensors  *.onnx  *.ckpt   # model weights
```
- Object-storage-backed vector DB (from `docker-compose.dgx-spark.yml`): `milvus-minio` +
  `milvus-etcd` + `milvus`, volumes `minio_data` / `etcd_data` / `milvus_data`. MinIO creds come from
  `${HCLS_MINIO_USER:?}` / `${HCLS_MINIO_PASSWORD:?}` — never defaulted.
- Variant store: VCF → DuckDB (Engine 1, :8575); keep it a local `.duckdb` file, gitignored.

## Pitfalls (single-box DGX Spark / ARM / this factory)
- **The 711 GB tree is single-box** — there is no replica; losing the box loses the data unless the
  reference package can be re-downloaded and outputs are backed up separately.
- **Unified memory is shared with storage caches** — a big DuckDB scan or Milvus load competes with
  models for the same 128 GB; size queries accordingly.
- **A stray `git add -f` of weights bloats history irreversibly** — the oversize/neutrality hooks exist
  precisely to stop this; don't `--no-verify` past them for data.
- **Anonymous volumes vanish on recreate** — always name the volume for anything stateful.

## Related
- Pillars: 05-structured-databases, 06-vector-databases-and-embeddings, 11-security-and-secrets
- build-housekeeping-standards
