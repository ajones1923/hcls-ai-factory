"""Create all 14 Milvus collections for the Rare Disease Diagnostic Agent.

Creates collections with proper schemas and IVF_FLAT + COSINE indexes.

Collections:
  1. rd_phenotypes         -- HPO phenotype terms
  2. rd_diseases           -- OMIM/Orphanet disease entries
  3. rd_genes              -- Gene-disease associations
  4. rd_variants           -- ClinVar/ACMG variant data
  5. rd_literature         -- PubMed rare disease publications
  6. rd_trials             -- Clinical trials for rare diseases
  7. rd_therapies          -- Gene therapies (approved + investigational)
  8. rd_case_reports       -- Published case reports
  9. rd_guidelines         -- ACMG/ACOG clinical guidelines
  10. rd_pathways           -- Biochemical/metabolic pathways
  11. rd_registries         -- Patient registries and natural history studies
  12. rd_natural_history    -- Disease natural history data
  13. rd_newborn_screening  -- Newborn screening panel data
  14. genomic_evidence      -- Shared genomic evidence collection

Usage:
    python scripts/setup_collections.py

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logger = logging.getLogger(__name__)

# ===================================================================
# CONSTANTS
# ===================================================================

EMBEDDING_DIM = 384  # BGE-small-en-v1.5
INDEX_TYPE = "IVF_FLAT"
METRIC_TYPE = "COSINE"
NLIST = 128

# ===================================================================
# COLLECTION SCHEMAS
# ===================================================================

COLLECTION_SCHEMAS = {
    "rd_phenotypes": {
        "description": "HPO phenotype terms with definitions and IC scores",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "hpo_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "name", "dtype": "VARCHAR", "max_length": 256},
            {"name": "ic_score", "dtype": "FLOAT"},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 20000,
        "search_weight": 0.12,
    },
    "rd_diseases": {
        "description": "OMIM/Orphanet rare disease entries",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "disease_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "disease_name", "dtype": "VARCHAR", "max_length": 512},
            {"name": "gene", "dtype": "VARCHAR", "max_length": 64},
            {"name": "inheritance", "dtype": "VARCHAR", "max_length": 64},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 10000,
        "search_weight": 0.11,
    },
    "rd_genes": {
        "description": "Gene-disease associations and functional data",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "gene_symbol", "dtype": "VARCHAR", "max_length": 64},
            {"name": "chromosome", "dtype": "VARCHAR", "max_length": 32},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 5000,
        "search_weight": 0.10,
    },
    "rd_variants": {
        "description": "ClinVar/ACMG variant classifications",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "gene", "dtype": "VARCHAR", "max_length": 64},
            {"name": "variant", "dtype": "VARCHAR", "max_length": 128},
            {"name": "classification", "dtype": "VARCHAR", "max_length": 64},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 500000,
        "search_weight": 0.10,
    },
    "rd_literature": {
        "description": "PubMed rare disease publications",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "pmid", "dtype": "VARCHAR", "max_length": 32},
            {"name": "title", "dtype": "VARCHAR", "max_length": 512},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 50000,
        "search_weight": 0.08,
    },
    "rd_trials": {
        "description": "Clinical trials for rare diseases",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "disease", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 10000,
        "search_weight": 0.06,
    },
    "rd_therapies": {
        "description": "Gene therapies (approved + investigational)",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "drug_name", "dtype": "VARCHAR", "max_length": 256},
            {"name": "target_gene", "dtype": "VARCHAR", "max_length": 64},
            {"name": "approval_status", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 500,
        "search_weight": 0.07,
    },
    "rd_case_reports": {
        "description": "Published rare disease case reports",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "disease", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 20000,
        "search_weight": 0.07,
    },
    "rd_guidelines": {
        "description": "ACMG/ACOG clinical genetics guidelines",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "guideline_id", "dtype": "VARCHAR", "max_length": 64},
            {"name": "organization", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 2000,
        "search_weight": 0.06,
    },
    "rd_pathways": {
        "description": "Biochemical and metabolic pathways",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "pathway_id", "dtype": "VARCHAR", "max_length": 64},
            {"name": "pathway_name", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 1000,
        "search_weight": 0.06,
    },
    "rd_registries": {
        "description": "Patient registries and natural history studies",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "registry_name", "dtype": "VARCHAR", "max_length": 256},
            {"name": "disease", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 500,
        "search_weight": 0.04,
    },
    "rd_natural_history": {
        "description": "Disease natural history and progression data",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "disease", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 2000,
        "search_weight": 0.05,
    },
    "rd_newborn_screening": {
        "description": "Newborn screening panel conditions and protocols",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "condition", "dtype": "VARCHAR", "max_length": 256},
            {"name": "screening_method", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 200,
        "search_weight": 0.05,
    },
    "genomic_evidence": {
        "description": "Shared genomic evidence collection",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "gene", "dtype": "VARCHAR", "max_length": 64},
            {"name": "variant", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 100000,
        "search_weight": 0.03,
    },
}


def get_collection_names() -> list:
    """Return the list of all collection names."""
    return list(COLLECTION_SCHEMAS.keys())


def get_collection_schema(name: str) -> dict:
    """Return the schema for a specific collection."""
    return COLLECTION_SCHEMAS.get(name, {})


def setup_all_collections(milvus_host: str = "localhost", milvus_port: int = 19530) -> None:
    """Create all 14 Milvus collections with proper schemas and indexes.

    Args:
        milvus_host: Milvus server hostname.
        milvus_port: Milvus server port.
    """
    try:
        from pymilvus import connections, Collection, CollectionSchema, FieldSchema, DataType, utility
    except ImportError:
        logger.error("pymilvus is not installed. Run: pip install pymilvus")
        return

    logger.info("Connecting to Milvus at %s:%d", milvus_host, milvus_port)
    connections.connect(host=milvus_host, port=milvus_port)

    dtype_map = {
        "INT64": DataType.INT64,
        "VARCHAR": DataType.VARCHAR,
        "FLOAT_VECTOR": DataType.FLOAT_VECTOR,
        "FLOAT": DataType.FLOAT,
    }

    for coll_name, config in COLLECTION_SCHEMAS.items():
        if utility.has_collection(coll_name):
            logger.info("Collection '%s' already exists, skipping.", coll_name)
            continue

        logger.info("Creating collection '%s': %s", coll_name, config["description"])

        fields = []
        for field_def in config["fields"]:
            kwargs = {
                "name": field_def["name"],
                "dtype": dtype_map[field_def["dtype"]],
            }
            if field_def.get("is_primary"):
                kwargs["is_primary"] = True
            if field_def.get("auto_id"):
                kwargs["auto_id"] = True
            if "max_length" in field_def:
                kwargs["max_length"] = field_def["max_length"]
            if "dim" in field_def:
                kwargs["dim"] = field_def["dim"]

            fields.append(FieldSchema(**kwargs))

        schema = CollectionSchema(
            fields=fields,
            description=config["description"],
        )

        collection = Collection(name=coll_name, schema=schema)

        index_params = {
            "metric_type": METRIC_TYPE,
            "index_type": INDEX_TYPE,
            "params": {"nlist": NLIST},
        }
        collection.create_index(field_name="embedding", index_params=index_params)
        logger.info("Collection '%s' created with index.", coll_name)

    logger.info("All %d collections set up successfully.", len(COLLECTION_SCHEMAS))
    connections.disconnect("default")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    setup_all_collections()
