"""One Milvus similarity search for the RAG engines.

Five subjects carried `_search_collection` — 88 lines each, identical in every line of code and
differing only in a docstring example and one comment. It is the single hottest path in the
platform: every clinical answer in the factory passes through it, once per collection searched.

Two traps are baked into the body and are the reason this must not be re-derived per subject:

  * `MilvusClient.search` takes `search_params`, NOT the ORM API's `param`. Passing the wrong one
    raises inside the try block, which returns [] — a silent empty result set, not an error.
  * The two client APIs return different shapes. `MilvusClient` yields plain dicts
    ({id, distance, entity}); the ORM API yields Hit objects with .id/.score. Both are handled so
    the engines work whichever client is injected.

The failure mode throughout is an empty list plus a warning: a retrieval failure must degrade the
answer, never break the request. That is also why it is worth having in one tested place — a
search that silently returns nothing looks identical to a corpus with no match.
"""
from __future__ import annotations

import logging
from typing import Any, List, Optional

logger = logging.getLogger(__name__)


def search_collection(
    client: Any,
    collection_name: str,
    query_vector: List[float],
    top_k: int,
    filter_expr: Optional[str] = None,
    *,
    metric_type: str = "COSINE",
    nprobe: int = 16,
) -> List[dict]:
    """Search one Milvus collection and flatten the hits.

    Args:
        client: a `MilvusClient` (or an ORM-style client exposing `.search`).
        collection_name: Milvus collection name.
        query_vector: query embedding.
        top_k: maximum number of results.
        filter_expr: optional Milvus boolean filter expression
            (e.g. 'phase == "Phase 3"', 'modality == "echocardiography"').

    Returns:
        List of result dicts carrying `id`, `score`, the entity fields, and `metadata`.
        An empty list on any failure.
    """
    try:
        search_kwargs = {
            "collection_name": collection_name,
            "data": [query_vector],
            "anns_field": "embedding",
            # MilvusClient.search takes `search_params`, not the ORM API's `param`.
            "search_params": {"metric_type": metric_type, "params": {"nprobe": nprobe}},
            "limit": top_k,
            "output_fields": ["*"],
        }
        if filter_expr:
            search_kwargs["filter"] = filter_expr

        results = client.search(**search_kwargs)

        flat: List[dict] = []
        if results and len(results) > 0:
            for hit in results[0]:
                if isinstance(hit, dict):
                    record = {
                        "id": str(hit.get("id", "")),
                        "score": float(hit.get("distance", hit.get("score", 0.0)) or 0.0),
                    }
                    ent = hit.get("entity") or {}
                    if isinstance(ent, dict):
                        for k, v in ent.items():
                            if k != "embedding":
                                record[k] = v
                    record["metadata"] = {k: v for k, v in record.items()
                                          if k not in ("id", "score", "metadata")}
                    flat.append(record)
                    continue
                record = {
                    "id": str(hit.id),
                    "score": float(hit.score) if hasattr(hit, "score") else 0.0,
                }
                if hasattr(hit, "entity"):
                    entity = hit.entity
                    if hasattr(entity, "fields"):
                        for name, value in entity.fields.items():
                            if name != "embedding":
                                record[name] = value
                    elif isinstance(entity, dict):
                        for k, v in entity.items():
                            if k != "embedding":
                                record[k] = v
                flat.append(record)
        return flat

    except Exception as exc:
        logger.warning("Search failed for collection '%s': %s", collection_name, exc)
        return []
