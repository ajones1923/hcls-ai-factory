"""
Precision Autoimmune Agent — Milvus Collection Manager

Manages 14 autoimmune-domain vector collections:
  1. autoimmune_clinical_documents    — Ingested patient records (PDFs)
  2. autoimmune_patient_labs          — Lab results with flag analysis
  3. autoimmune_autoantibody_panels   — Autoantibody test result panels
  4. autoimmune_hla_associations      — HLA allele → disease risk mapping
  5. autoimmune_disease_criteria      — ACR/EULAR classification criteria
  6. autoimmune_disease_activity      — Activity scoring reference (DAS28, SLEDAI, etc.)
  7. autoimmune_flare_patterns        — Flare prediction biomarker patterns
  8. autoimmune_biologic_therapies    — Biologic drug database with PGx
  9. autoimmune_pgx_rules             — Pharmacogenomic dosing rules
 10. autoimmune_clinical_trials       — Autoimmune clinical trials
 11. autoimmune_literature            — Published literature
 12. autoimmune_patient_timelines     — Patient diagnostic timelines
 13. autoimmune_cross_disease         — Cross-disease / overlap syndromes
 14. genomic_evidence                 — Shared read-only collection

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

from loguru import logger
from pymilvus import (
    Collection,
    CollectionSchema,
    DataType,
    FieldSchema,
    MilvusException,
    connections,
    utility,
)

# ── Index / search parameters ────────────────────────────────────────────
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}

_DIM = 384  # BGE-small-en-v1.5


# ── Schema helpers ────────────────────────────────────────────────────────
def _pk(max_len: int = 128) -> FieldSchema:
    return FieldSchema("id", DataType.VARCHAR, is_primary=True, max_length=max_len)


def _embedding(dim: int = _DIM) -> FieldSchema:
    return FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=dim)


def _varchar(name: str, max_length: int = 3000) -> FieldSchema:
    return FieldSchema(name, DataType.VARCHAR, max_length=max_length)


def _int(name: str) -> FieldSchema:
    return FieldSchema(name, DataType.INT64)


def _float(name: str) -> FieldSchema:
    return FieldSchema(name, DataType.FLOAT)


# ── Collection schemas ────────────────────────────────────────────────────

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {}


def _register(name: str, fields: list, description: str = ""):
    schema = CollectionSchema(fields=fields, description=description)
    COLLECTION_SCHEMAS[name] = schema


# 1. Clinical documents (ingested PDFs)
_register(
    "autoimmune_clinical_documents",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("patient_id", 64),
        _varchar("doc_type", 128),       # progress_note, lab_report, imaging, pathology, etc.
        _varchar("specialty", 128),       # rheumatology, neurology, etc.
        _varchar("provider", 256),
        _varchar("visit_date", 32),
        _varchar("source_file", 512),
        _int("page_number"),
        _int("chunk_index"),
    ],
    "Ingested clinical documents for autoimmune patients",
)

# 2. Patient labs
_register(
    "autoimmune_patient_labs",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("patient_id", 64),
        _varchar("test_name", 256),
        _float("value"),
        _varchar("unit", 64),
        _varchar("reference_range", 128),
        _varchar("flag", 32),             # normal, high, low, critical
        _varchar("collection_date", 32),
        _varchar("panel_name", 256),
    ],
    "Laboratory results with reference range analysis",
)

# 3. Autoantibody panels
_register(
    "autoimmune_autoantibody_panels",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("antibody_name", 128),
        _varchar("associated_diseases", 1024),
        _float("sensitivity"),
        _float("specificity"),
        _varchar("pattern", 128),         # homogeneous, speckled, etc.
        _varchar("clinical_significance", 2000),
        _varchar("interpretation_guide", 2000),
    ],
    "Autoantibody reference panels with disease associations",
)

# 4. HLA associations
_register(
    "autoimmune_hla_associations",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("allele", 64),
        _varchar("disease", 256),
        _float("odds_ratio"),
        _varchar("population", 128),
        _varchar("pmid", 32),
        _varchar("mechanism", 1024),
        _varchar("clinical_implication", 2000),
    ],
    "HLA allele to autoimmune disease associations with odds ratios",
)

# 5. Disease criteria
_register(
    "autoimmune_disease_criteria",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("disease", 256),
        _varchar("criteria_set", 256),    # ACR/EULAR 2010, SLICC 2012, etc.
        _varchar("criteria_type", 64),    # classification, diagnostic
        _int("year"),
        _varchar("required_score", 128),
        _varchar("criteria_items", 3000),
        _varchar("sensitivity_specificity", 256),
    ],
    "ACR/EULAR classification and diagnostic criteria",
)

# 6. Disease activity scoring
_register(
    "autoimmune_disease_activity",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("score_name", 128),      # DAS28-CRP, SLEDAI-2K, BASDAI, etc.
        _varchar("disease", 256),
        _varchar("components", 2000),
        _varchar("thresholds", 1024),     # JSON: remission/low/moderate/high cutoffs
        _varchar("interpretation", 2000),
        _varchar("monitoring_frequency", 512),
    ],
    "Disease activity scoring systems and interpretation guides",
)

# 7. Flare patterns
_register(
    "autoimmune_flare_patterns",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("disease", 256),
        _varchar("biomarker_pattern", 2000),
        _varchar("early_warning_signs", 2000),
        _varchar("typical_timeline", 512),
        _varchar("protective_factors", 1024),
        _varchar("intervention_triggers", 1024),
    ],
    "Flare prediction biomarker patterns and early warning signs",
)

# 8. Biologic therapies
_register(
    "autoimmune_biologic_therapies",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("drug_name", 128),
        _varchar("drug_class", 128),
        _varchar("mechanism", 512),
        _varchar("indicated_diseases", 1024),
        _varchar("pgx_considerations", 2000),
        _varchar("contraindications", 1024),
        _varchar("monitoring", 2000),
        _varchar("dosing", 512),
        _varchar("evidence_level", 64),
    ],
    "Biologic therapy database with pharmacogenomic considerations",
)

# 9. PGx rules
_register(
    "autoimmune_pgx_rules",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("gene", 64),
        _varchar("variant", 128),
        _varchar("drug", 128),
        _varchar("phenotype", 128),       # poor/intermediate/normal/rapid metabolizer
        _varchar("recommendation", 2000),
        _varchar("evidence_level", 64),   # CPIC Level A/B/C
        _varchar("pmid", 256),
    ],
    "Pharmacogenomic dosing rules for autoimmune therapies",
)

# 10. Clinical trials
_register(
    "autoimmune_clinical_trials",
    [
        _pk(256), _embedding(),
        _varchar("text_chunk"),
        _varchar("title", 1024),
        _varchar("nct_id", 32),
        _varchar("phase", 32),
        _varchar("status", 64),
        _varchar("disease", 256),
        _varchar("intervention", 512),
        _varchar("biomarker_criteria", 1024),
        _int("enrollment"),
        _int("start_year"),
        _varchar("sponsor", 256),
    ],
    "Autoimmune disease clinical trials",
)

# 11. Literature
_register(
    "autoimmune_literature",
    [
        _pk(256), _embedding(),
        _varchar("text_chunk"),
        _varchar("title", 1024),
        _varchar("authors", 1024),
        _varchar("journal", 256),
        _int("year"),
        _varchar("pmid", 32),
        _varchar("disease_focus", 256),
        _varchar("keywords", 1024),
        _varchar("abstract_summary", 3000),
    ],
    "Published autoimmune literature and research",
)

# 12. Patient timelines
_register(
    "autoimmune_patient_timelines",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("patient_id", 64),
        _varchar("event_type", 128),      # symptom_onset, diagnosis, treatment_start, flare, etc.
        _varchar("event_date", 32),
        _varchar("description", 2000),
        _varchar("provider", 256),
        _varchar("specialty", 128),
        _int("days_from_first_symptom"),
    ],
    "Patient diagnostic timeline events for odyssey analysis",
)

# 13. Cross-disease patterns
_register(
    "autoimmune_cross_disease",
    [
        _pk(), _embedding(),
        _varchar("text_chunk"),
        _varchar("primary_disease", 256),
        _varchar("associated_conditions", 1024),
        _varchar("shared_pathways", 1024),
        _varchar("shared_biomarkers", 1024),
        _varchar("overlap_criteria", 2000),
        _float("co_occurrence_rate"),
    ],
    "Cross-disease overlap syndromes and shared pathogenic mechanisms",
)


# ── Collection Manager ────────────────────────────────────────────────────

class AutoimmuneCollectionManager:
    """Manages Milvus collections for the Autoimmune Agent."""

    def __init__(self, host: str = "localhost", port: int = 19530, embedding_dim: int = _DIM):
        self.host = host
        self.port = port
        self.embedding_dim = embedding_dim
        self._alias = "autoimmune_agent"
        self._connected = False

    # ── Connection ────────────────────────────────────────────────────
    def connect(self) -> None:
        if self._connected:
            return
        try:
            connections.connect(alias=self._alias, host=self.host, port=self.port)
            self._connected = True
            logger.info(f"Connected to Milvus at {self.host}:{self.port}")
        except MilvusException as exc:
            logger.error(f"Milvus connection failed: {exc}")
            raise

    def disconnect(self) -> None:
        if self._connected:
            connections.disconnect(self._alias)
            self._connected = False
            logger.info("Disconnected from Milvus")

    def _ensure_connected(self) -> None:
        if not self._connected:
            self.connect()

    # ── Collection lifecycle ──────────────────────────────────────────
    def create_collection(self, name: str, drop_existing: bool = False) -> Collection:
        self._ensure_connected()
        schema = COLLECTION_SCHEMAS.get(name)
        if schema is None:
            raise ValueError(f"Unknown collection: {name}")

        if utility.has_collection(name, using=self._alias):
            if drop_existing:
                utility.drop_collection(name, using=self._alias)
                logger.info(f"Dropped existing collection: {name}")
            else:
                coll = Collection(name, using=self._alias)
                coll.load()
                return coll

        coll = Collection(name, schema=schema, using=self._alias)
        coll.create_index("embedding", INDEX_PARAMS)
        coll.load()
        logger.info(f"Created collection: {name} ({len(schema.fields)} fields)")
        return coll

    def create_all_collections(self, drop_existing: bool = False) -> Dict[str, Collection]:
        self._ensure_connected()
        collections = {}
        for name in COLLECTION_SCHEMAS:
            # Skip genomic_evidence — shared, not ours to create
            if name == "genomic_evidence":
                if utility.has_collection(name, using=self._alias):
                    collections[name] = Collection(name, using=self._alias)
                    collections[name].load()
                continue
            collections[name] = self.create_collection(name, drop_existing)
        logger.info(f"Created/loaded {len(collections)} collections")
        return collections

    def list_collections(self) -> List[str]:
        self._ensure_connected()
        all_colls = utility.list_collections(using=self._alias)
        return [c for c in all_colls if c.startswith("autoimmune_") or c == "genomic_evidence"]

    @property
    def is_connected(self) -> bool:
        """Check if Milvus connection is active."""
        return self._connected

    def get_collection_stats(self) -> Dict[str, int]:
        self._ensure_connected()
        stats = {}
        for name in self.list_collections():
            try:
                coll = Collection(name, using=self._alias)
                stats[name] = coll.num_entities
            except Exception as exc:
                logger.warning(f"Failed to get stats for {name}: {exc}")
                stats[name] = -1
        return stats

    def get_collection_count(self, name: str) -> int:
        self._ensure_connected()
        try:
            coll = Collection(name, using=self._alias)
            return coll.num_entities
        except Exception as exc:
            logger.warning(f"Failed to count {name}: {exc}")
            return 0

    # ── Insert ────────────────────────────────────────────────────────
    def insert(self, collection_name: str, records: List[Dict[str, Any]]) -> int:
        """Insert records into a collection. Each record must include 'embedding'."""
        self._ensure_connected()
        if not records:
            return 0

        coll = Collection(collection_name, using=self._alias)
        schema_fields = coll.schema.fields
        field_defaults = {}
        for f in schema_fields:
            if f.name == "embedding":
                field_defaults[f.name] = None  # must be provided
            elif f.dtype == DataType.VARCHAR:
                field_defaults[f.name] = ""
            elif f.dtype in (DataType.INT64, DataType.INT32, DataType.INT16, DataType.INT8):
                field_defaults[f.name] = 0
            elif f.dtype in (DataType.FLOAT, DataType.DOUBLE):
                field_defaults[f.name] = 0.0
            else:
                field_defaults[f.name] = ""

        fields = [f.name for f in schema_fields]
        data = {f: [] for f in fields}
        skipped = 0

        for rec in records:
            # Skip records without a real embedding (don't insert zero vectors)
            emb = rec.get("embedding")
            if emb is None or (isinstance(emb, list) and len(emb) != self.embedding_dim):
                skipped += 1
                continue

            for f in fields:
                if f == "embedding":
                    data[f].append(emb)
                else:
                    val = rec.get(f, field_defaults.get(f, ""))
                    # Truncate strings to field max_length to prevent Milvus errors
                    if isinstance(val, str):
                        for sf in schema_fields:
                            if sf.name == f and hasattr(sf, "max_length") and sf.max_length:
                                val = val[:sf.max_length]
                                break
                    data[f].append(val)

        if not data[fields[0]]:
            return 0

        if skipped:
            logger.warning(f"Skipped {skipped} records without valid embeddings in {collection_name}")

        try:
            insert_data = [data[f] for f in fields]
            result = coll.insert(insert_data)
            coll.flush()
            return result.insert_count
        except MilvusException as exc:
            logger.error(f"Insert failed for {collection_name}: {exc}")
            raise

    def insert_batch(self, collection_name: str, records: List[Dict[str, Any]], batch_size: int = 100) -> int:
        """Insert records in batches."""
        total = 0
        for i in range(0, len(records), batch_size):
            batch = records[i : i + batch_size]
            total += self.insert(collection_name, batch)
        return total

    # ── Search ────────────────────────────────────────────────────────
    def search(
        self,
        collection_name: str,
        query_embedding: List[float],
        top_k: int = 5,
        filter_expr: Optional[str] = None,
        output_fields: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Search a single collection."""
        self._ensure_connected()
        try:
            coll = Collection(collection_name, using=self._alias)
            if output_fields is None:
                output_fields = [
                    f.name for f in coll.schema.fields
                    if f.name not in ("embedding",)
                ]
            results = coll.search(
                data=[query_embedding],
                anns_field="embedding",
                param=SEARCH_PARAMS,
                limit=top_k,
                expr=filter_expr,
                output_fields=output_fields,
            )
            hits = []
            for hit in results[0]:
                entry = {"id": hit.id, "score": hit.score}
                for field in output_fields:
                    if field != "id":
                        entry[field] = hit.entity.get(field, "")
                hits.append(entry)
            return hits
        except Exception as exc:
            logger.warning(f"Search failed on {collection_name}: {exc}")
            return []

    def search_all(
        self,
        query_embedding: List[float],
        top_k_per_collection: int = 5,
        collections: Optional[List[str]] = None,
        filter_exprs: Optional[Dict[str, str]] = None,
        score_threshold: float = 0.0,
        max_workers: int = 6,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Search multiple collections in parallel."""
        self._ensure_connected()
        if collections is None:
            collections = self.list_collections()
        if filter_exprs is None:
            filter_exprs = {}

        results: Dict[str, List[Dict[str, Any]]] = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self.search,
                    coll_name,
                    query_embedding,
                    top_k_per_collection,
                    filter_exprs.get(coll_name),
                ): coll_name
                for coll_name in collections
            }
            for future in as_completed(futures):
                coll_name = futures[future]
                try:
                    hits = future.result()
                    if score_threshold > 0:
                        hits = [h for h in hits if h["score"] >= score_threshold]
                    if hits:
                        results[coll_name] = hits
                except Exception as exc:
                    logger.warning(f"Parallel search failed for {coll_name}: {exc}")

        return results
