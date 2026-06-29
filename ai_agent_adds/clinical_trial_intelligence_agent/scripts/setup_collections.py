"""Create all 14 Milvus collections for the Clinical Trial Intelligence Agent.

Creates collections with proper schemas and IVF_FLAT + COSINE indexes.

Collections:
  1. trial_protocols      -- Study protocols from ClinicalTrials.gov
  2. trial_eligibility    -- Eligibility criteria (inclusion/exclusion)
  3. trial_endpoints      -- Primary and secondary endpoints
  4. trial_sites          -- Trial sites and investigator facilities
  5. trial_investigators  -- Principal investigators and study teams
  6. trial_results        -- Published trial results and outcomes
  7. trial_regulatory     -- FDA/EMA regulatory documents and decisions
  8. trial_literature     -- PubMed clinical trial publications
  9. trial_biomarkers     -- Biomarker data and enrichment strategies
  10. trial_safety        -- Adverse event and safety signal data
  11. trial_rwe           -- Real-world evidence and observational data
  12. trial_adaptive      -- Adaptive design and interim analysis data
  13. trial_guidelines    -- ICH/FDA/EMA guidance documents
  14. genomic_evidence    -- Shared genomic evidence collection

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
    "trial_protocols": {
        "description": "Clinical trial protocols from ClinicalTrials.gov",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "title", "dtype": "VARCHAR", "max_length": 1024},
            {"name": "phase", "dtype": "VARCHAR", "max_length": 64},
            {"name": "status", "dtype": "VARCHAR", "max_length": 64},
            {"name": "sponsor", "dtype": "VARCHAR", "max_length": 256},
            {"name": "therapeutic_area", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 50000,
        "search_weight": 0.10,
    },
    "trial_eligibility": {
        "description": "Trial eligibility criteria (inclusion/exclusion)",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "criterion_type", "dtype": "VARCHAR", "max_length": 32},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 200000,
        "search_weight": 0.09,
    },
    "trial_endpoints": {
        "description": "Clinical trial primary and secondary endpoints",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "endpoint_type", "dtype": "VARCHAR", "max_length": 32},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 100000,
        "search_weight": 0.08,
    },
    "trial_sites": {
        "description": "Clinical trial sites and locations",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "facility", "dtype": "VARCHAR", "max_length": 256},
            {"name": "city", "dtype": "VARCHAR", "max_length": 128},
            {"name": "country", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 500000,
        "search_weight": 0.07,
    },
    "trial_investigators": {
        "description": "Principal investigators and study teams",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "name", "dtype": "VARCHAR", "max_length": 256},
            {"name": "institution", "dtype": "VARCHAR", "max_length": 256},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 100000,
        "search_weight": 0.05,
    },
    "trial_results": {
        "description": "Published trial results and outcomes",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "outcome_type", "dtype": "VARCHAR", "max_length": 64},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 50000,
        "search_weight": 0.09,
    },
    "trial_regulatory": {
        "description": "FDA/EMA regulatory documents and decisions",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "agency", "dtype": "VARCHAR", "max_length": 32},
            {"name": "drug", "dtype": "VARCHAR", "max_length": 256},
            {"name": "decision", "dtype": "VARCHAR", "max_length": 128},
            {"name": "date", "dtype": "VARCHAR", "max_length": 32},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 10000,
        "search_weight": 0.07,
    },
    "trial_literature": {
        "description": "PubMed clinical trial publications",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "pmid", "dtype": "VARCHAR", "max_length": 32},
            {"name": "title", "dtype": "VARCHAR", "max_length": 1024},
            {"name": "journal", "dtype": "VARCHAR", "max_length": 256},
            {"name": "year", "dtype": "VARCHAR", "max_length": 8},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 100000,
        "search_weight": 0.08,
    },
    "trial_biomarkers": {
        "description": "Biomarker data and enrichment strategies",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "biomarker_name", "dtype": "VARCHAR", "max_length": 128},
            {"name": "therapeutic_area", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 20000,
        "search_weight": 0.07,
    },
    "trial_safety": {
        "description": "Adverse event and safety signal data",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "event_type", "dtype": "VARCHAR", "max_length": 256},
            {"name": "severity", "dtype": "VARCHAR", "max_length": 32},
            {"name": "nct_id", "dtype": "VARCHAR", "max_length": 32},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 200000,
        "search_weight": 0.08,
    },
    "trial_rwe": {
        "description": "Real-world evidence and observational data",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "data_source", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 50000,
        "search_weight": 0.06,
    },
    "trial_adaptive": {
        "description": "Adaptive design and interim analysis data",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 4096},
            {"name": "design_type", "dtype": "VARCHAR", "max_length": 128},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 10000,
        "search_weight": 0.05,
    },
    "trial_guidelines": {
        "description": "ICH/FDA/EMA guidance documents",
        "fields": [
            {"name": "id", "dtype": "INT64", "is_primary": True, "auto_id": True},
            {"name": "embedding", "dtype": "FLOAT_VECTOR", "dim": EMBEDDING_DIM},
            {"name": "text", "dtype": "VARCHAR", "max_length": 8192},
            {"name": "guideline_id", "dtype": "VARCHAR", "max_length": 64},
            {"name": "agency", "dtype": "VARCHAR", "max_length": 32},
            {"name": "title", "dtype": "VARCHAR", "max_length": 512},
            {"name": "source", "dtype": "VARCHAR", "max_length": 64},
        ],
        "estimated_records": 5000,
        "search_weight": 0.08,
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

    for coll_name, config in COLLECTION_SCHEMAS.items():
        if utility.has_collection(coll_name):
            logger.info("Collection '%s' already exists, skipping.", coll_name)
            continue

        logger.info("Creating collection '%s': %s", coll_name, config["description"])

        dtype_map = {
            "INT64": DataType.INT64,
            "VARCHAR": DataType.VARCHAR,
            "FLOAT_VECTOR": DataType.FLOAT_VECTOR,
        }

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

        # Create index on embedding field
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
