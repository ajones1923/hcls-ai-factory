# Precision Autoimmune Intelligence Agent
#
# Core modules:
#   agent.py             — Main orchestrator (AutoimmuneAgent)
#   models.py            — Pydantic data models
#   knowledge.py         — Static knowledge base (HLA, autoantibodies, biologics)
#   collections.py       — Milvus collection manager (14 collections)
#   rag_engine.py        — Multi-collection RAG engine
#   document_processor.py — Clinical PDF ingestion pipeline
#   diagnostic_engine.py — Classification criteria, differential diagnosis
#   timeline_builder.py  — Diagnostic odyssey timeline construction
