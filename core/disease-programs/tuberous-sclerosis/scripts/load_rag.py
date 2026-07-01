"""
Ingest the TSC literature corpus into the RAG store (PRD §2.5.4; scripts/load_rag.py).
Targets Milvus when TSC_USE_MILVUS=1 (prod), else the in-memory store. Uses real
BGE/BiomedBERT embeddings when TSC_USE_REAL_EMBEDDINGS=1, else the deterministic embedder.

  python3 scripts/load_rag.py
  TSC_USE_MILVUS=1 TSC_USE_REAL_EMBEDDINGS=1 python3 scripts/load_rag.py   # prod
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.rag.corpus import SEED_CORPUS          # noqa: E402
from src.rag.embeddings import mode             # noqa: E402
from src.rag.retriever import make_vector_store  # noqa: E402


def main() -> None:
    store = make_vector_store()
    store.upsert(SEED_CORPUS)
    print(f"Ingested {len(SEED_CORPUS)} chunks into {type(store).__name__} "
          f"(embeddings: {mode()})")
    print("Smoke query — 'everolimus reduces seizures in TSC':")
    for h in store.search("everolimus reduces seizures in TSC", k=3):
        print(f"  {h.get('score')}  [{h.get('partition')}]  {h.get('source_uri')}")


if __name__ == "__main__":
    main()
