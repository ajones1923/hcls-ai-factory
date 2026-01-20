# Phase 4: RAG Chat Pipeline - Summary

## Overview

**Purpose:** Query genomic variants via natural language and identify drug targets
**Input:** VCF from Phase 1 (Genomics Pipeline)
**Output:** Target hypotheses with PDB IDs for Phase 5 (Cryo-EM)
**Hardware:** NVIDIA DGX Spark (128GB unified memory)

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     PHASE 4: RAG CHAT PIPELINE                  │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│   VCF File (4M+ variants)                                       │
│        │                                                        │
│        ▼                                                        │
│   ┌─────────────┐    ┌─────────────┐    ┌─────────────┐        │
│   │ VCF Parser  │───▶│ VEP Annotate│───▶│  Embedder   │        │
│   │  (cyvcf2)   │    │  (REST API) │    │ (bge-small) │        │
│   └─────────────┘    └─────────────┘    └─────────────┘        │
│                                                │                │
│                                                ▼                │
│                                         ┌─────────────┐        │
│                                         │   Milvus    │        │
│                                         │ Vector DB   │        │
│                                         └─────────────┘        │
│                                                │                │
│   User Question ──────────────────────────────▼                │
│        │                              ┌─────────────┐          │
│        │         Retrieve Top 10      │    RAG      │          │
│        └─────────────────────────────▶│   Engine    │          │
│                                       └─────────────┘          │
│                                                │                │
│                                                ▼                │
│                                         ┌─────────────┐        │
│                                         │    vLLM     │        │
│                                         │ Llama 70B   │        │
│                                         └─────────────┘        │
│                                                │                │
│                                                ▼                │
│                                         ┌─────────────┐        │
│                                         │   Answer    │        │
│                                         │ + Evidence  │        │
│                                         └─────────────┘        │
│                                                │                │
│                                                ▼                │
│                                         ┌─────────────┐        │
│                                         │   Target    │        │
│                                         │ Hypotheses  │        │
│                                         └─────────────┘        │
│                                                │                │
│                                                ▼                │
│                                         ┌─────────────┐        │
│                                         │ Export JSON │───────▶ Phase 5
│                                         │  (PDB IDs)  │        │
│                                         └─────────────┘        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Components Built

### Docker Services
| Service | Port | Purpose |
|---------|------|---------|
| Milvus | 19530 | Vector database |
| Attu | 8000 | Milvus web UI |
| vLLM | 8080 | Local LLM server (Llama 70B) |

### Python Modules
| Module | Purpose |
|--------|---------|
| `src/vcf_parser.py` | Parse VCF files (cyvcf2 + fallback) |
| `src/annotator.py` | VEP annotation (REST + local) |
| `src/embedder.py` | Text embeddings (sentence-transformers) |
| `src/milvus_client.py` | Vector database operations |
| `src/llm_client.py` | LLM interface (vLLM/Anthropic/OpenAI) |
| `src/rag_engine.py` | RAG orchestration |
| `src/target_hypothesis.py` | Target management & export |

### Interfaces
| Interface | Port | Purpose |
|-----------|------|---------|
| Portal | 5001 | Management & monitoring |
| Streamlit Chat | 8501 | RAG chat interface |

### Scripts
| Script | Purpose |
|--------|---------|
| `scripts/ingest_vcf.py` | Load VCF into Milvus |
| `scripts/run_chat.py` | Launch Streamlit |
| `scripts/setup_demo_targets.py` | GTC demo prep |

---

## Steps to Complete

### Step 1: Environment Setup
```bash
cd /home/adam/transfer/rag-chat-pipeline
cp .env.example .env
# Edit .env: Add HF_TOKEN
./run.sh setup
```
- [ ] .env configured with HF_TOKEN
- [ ] Llama 3.1 70B license accepted on Hugging Face
- [ ] Python dependencies installed

### Step 2: Start Services
```bash
./run.sh start
```
- [ ] Milvus running (check: http://localhost:8000)
- [ ] vLLM running with 70B model loaded
- [ ] Health checks passing

### Step 3: Ingest VCF
```bash
# Test with small batch first
./run.sh ingest --limit 1000

# Full ingestion
./run.sh ingest
```
- [ ] Test ingestion successful (1000 variants)
- [ ] Full ingestion complete (4M+ variants)
- [ ] Verify in Attu UI

### Step 4: Test Chat Interface
```bash
./run.sh chat
```
- [ ] Streamlit opens at http://localhost:8501
- [ ] Query: "What variants are in EGFR?" returns results
- [ ] Query: "What are high impact variants?" works
- [ ] Evidence section shows retrieved variants
- [ ] Performance metrics display (TTFT, tokens/sec)

### Step 5: Test Target Hypotheses
- [ ] Navigate to Target Hypotheses tab
- [ ] Add target (EGFR, PDB: 7SYE)
- [ ] Verify saved to `data/targets/hypotheses.json`
- [ ] Change status to "validated"
- [ ] Export for Phase 5
- [ ] Verify `data/targets/targets_for_phase5.json`

### Step 6: Test Portal
```bash
./run.sh portal
```
- [ ] Portal opens at http://localhost:5001
- [ ] Service status shows green
- [ ] VCF Preview tab works
- [ ] Targets tab shows saved targets
- [ ] Export button works

### Step 7: Demo Preparation
```bash
./run.sh demo-setup
```
- [ ] Demo targets created (EGFR, BRAF, PIK3CA, KRAS, JAK2)
- [ ] Export ready for Phase 5

### Step 8: End-to-End Verification
- [ ] Full demo flow works smoothly
- [ ] Response times acceptable (<30 sec)
- [ ] Export JSON valid for Phase 5 handoff

---

## CLI Commands

```bash
./run.sh setup        # Install dependencies
./run.sh start        # Start Milvus + vLLM
./run.sh stop         # Stop all services
./run.sh status       # Check service status
./run.sh ingest       # Ingest VCF into Milvus
./run.sh chat         # Start Streamlit (port 8501)
./run.sh portal       # Start Portal (port 5001)
./run.sh demo-setup   # Setup GTC demo targets
./run.sh help         # Show help
```

---

## Service URLs

| Service | URL | Purpose |
|---------|-----|---------|
| Portal | http://localhost:5001 | Management interface |
| Chat | http://localhost:8501 | RAG chat interface |
| Attu | http://localhost:8000 | Milvus admin UI |
| vLLM | http://localhost:8080 | LLM API endpoint |
| Milvus | localhost:19530 | Vector DB (gRPC) |

---

## Data Flow

```
Phase 1 Output                    Phase 4                           Phase 5 Input
─────────────────────────────────────────────────────────────────────────────────

HG002.genome.vcf.gz ──▶ Ingest ──▶ Milvus ──▶ RAG Chat ──▶ targets_for_phase5.json
     (4M+ variants)               (vectors)   (LLM 70B)      (PDB IDs, priorities)
```

---

## Key Files

```
/home/adam/transfer/rag-chat-pipeline/
├── .env                              # Configuration (HF_TOKEN, model)
├── docker-compose.yml                # Milvus, Attu, vLLM
├── run.sh                            # CLI interface
├── PHASE4_BUILD_PLAN.md              # Detailed build plan
├── PHASE4_SUMMARY.md                 # This file
│
├── config/
│   └── settings.py                   # Application settings
│
├── src/
│   ├── vcf_parser.py                 # VCF parsing
│   ├── annotator.py                  # VEP annotation
│   ├── embedder.py                   # Text embeddings
│   ├── milvus_client.py              # Vector database
│   ├── llm_client.py                 # LLM interface
│   ├── rag_engine.py                 # RAG orchestration
│   └── target_hypothesis.py          # Target management
│
├── app/
│   └── chat_ui.py                    # Streamlit chat
│
├── portal/
│   ├── app/server.py                 # Flask backend
│   ├── templates/index.html          # Portal UI
│   └── static/{css,js}/              # Assets
│
├── scripts/
│   ├── ingest_vcf.py                 # VCF ingestion
│   ├── run_chat.py                   # Chat launcher
│   └── setup_demo_targets.py         # Demo setup
│
└── data/
    └── targets/
        ├── hypotheses.json           # Saved targets
        └── targets_for_phase5.json   # Export for Phase 5
```

---

## Success Criteria

Phase 4 is **COMPLETE** when:

1. ✅ All services running (Milvus, vLLM, Portal, Chat)
2. ✅ VCF fully ingested into Milvus
3. ✅ Chat answers variant questions with evidence
4. ✅ Target hypotheses can be created and saved
5. ✅ Export generates valid JSON for Phase 5
6. ✅ Demo runs smoothly end-to-end

---

## Handoff to Phase 5

**Output File:** `data/targets/targets_for_phase5.json`

**Contains:**
```json
{
  "targets": [
    {
      "gene": "EGFR",
      "protein": "Epidermal growth factor receptor",
      "pdb_ids": ["7SYE", "1M17"],
      "rationale": "Kinase domain variant...",
      "priority": 5,
      "status": "validated"
    }
  ]
}
```

**Phase 5 (Cryo-EM) uses:**
- PDB IDs to fetch protein structures
- Priority to order analysis
- Rationale for context

---

## Timeline Estimate

| Step | Duration |
|------|----------|
| Setup & Config | 10 min |
| First vLLM start (model download) | 30-60 min |
| Subsequent vLLM starts | 2-5 min |
| VCF ingestion (4M variants) | 30-60 min |
| Testing & verification | 30 min |
| **Total first-time setup** | **~2 hours** |
| **Subsequent starts** | **~5 min** |
