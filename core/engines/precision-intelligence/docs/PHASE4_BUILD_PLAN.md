# Phase 4: RAG Chat Pipeline - Build & Test Plan

## Overview

Complete build and test plan for Phase 4 (RAG Chat Pipeline) to be demo-ready for GTC.

**Target Hardware:** NVIDIA DGX Spark (128GB unified memory)
**LLM Model:** Llama 3.1 70B Instruct
**VCF Source:** HG002 WGS from Phase 1 (Genomics Pipeline)

---

## Prerequisites Checklist

- [ ] Phase 1 (Genomics Pipeline) complete with VCF output
- [ ] Docker with NVIDIA Container Runtime installed
- [ ] Hugging Face account with token
- [ ] Llama 3.1 70B license accepted on Hugging Face
- [ ] VCF file exists: `genomics-pipeline/data/output/HG002.genome.vcf.gz`

---

## Step 1: Environment Setup

### 1.1 Configure Environment Variables
```bash
cd rag-chat-pipeline
cp .env.example .env
```

Edit `.env` and set:
```bash
HF_TOKEN=hf_xxxxxxxxxxxxxxxxxxxx  # Your Hugging Face token
VLLM_MODEL=meta-llama/Llama-3.1-70B-Instruct
VCF_INPUT_PATH=genomics-pipeline/data/output/HG002.genome.vcf.gz
```

### 1.2 Install Python Dependencies
```bash
./run.sh setup
```

**Expected output:** Virtual environment created, all packages installed.

### 1.3 Verify Setup
```bash
source venv/bin/activate
python -c "import pymilvus; import streamlit; import torch; print('All imports OK')"
```

---

## Step 2: Start Services

### 2.1 Start Milvus (Vector Database)
```bash
./run.sh start
```

Or manually:
```bash
docker-compose up -d milvus attu
```

### 2.2 Verify Milvus
```bash
# Check container is running
docker-compose ps

# Check Milvus health
curl http://localhost:9091/healthz

# Open Attu UI
# Browser: http://localhost:8000
```

**Expected:** Milvus container healthy, Attu UI accessible.

### 2.3 Start vLLM (Local LLM)
```bash
docker-compose up -d vllm
```

**First run:** Downloads ~140GB model. Monitor progress:
```bash
docker-compose logs -f vllm
```

**Expected time:** 20-60 minutes for download (depends on network).

### 2.4 Verify vLLM
```bash
# Check health endpoint
curl http://localhost:8080/health

# Check model loaded
curl http://localhost:8080/v1/models

# Test inference
curl http://localhost:8080/v1/completions \
  -H "Content-Type: application/json" \
  -d '{"model": "meta-llama/Llama-3.1-70B-Instruct", "prompt": "Hello", "max_tokens": 10}'
```

**Expected:** Model name returned, test completion works.

---

## Step 3: Ingest VCF Data

### 3.1 Test Ingestion (Small Batch)
```bash
./run.sh ingest --limit 1000
```

**Expected:** 1000 variants ingested in ~1-2 minutes.

### 3.2 Verify in Milvus
```bash
# Via Python
source venv/bin/activate
python -c "
from pymilvus import connections, Collection
connections.connect('default', host='localhost', port='19530')
coll = Collection('genomic_evidence')
print(f'Entities: {coll.num_entities}')
"
```

Or check via Attu UI: http://localhost:8000

### 3.3 Full Ingestion (All Variants)
```bash
./run.sh ingest
```

**Expected time:** 30-60 minutes for ~4M variants
**Monitor progress:** Script shows progress bar

### 3.4 Verify Full Ingestion
```bash
python -c "
from pymilvus import connections, Collection
connections.connect('default', host='localhost', port='19530')
coll = Collection('genomic_evidence')
print(f'Total entities: {coll.num_entities:,}')
"
```

**Expected:** Millions of variants indexed.

---

## Step 4: Test RAG Chat Interface

### 4.1 Start Streamlit
```bash
./run.sh chat
```

**Expected:** Opens browser to http://localhost:8501

### 4.2 Test Basic Queries

| Query | Expected Result |
|-------|-----------------|
| "What variants are in EGFR?" | List of EGFR variants with positions |
| "What are the high impact variants?" | Variants with HIGH impact level |
| "How many variants are on chromosome 17?" | Count and examples |

### 4.3 Test Filters
- [ ] Filter by gene (e.g., "BRCA1")
- [ ] Filter by chromosome
- [ ] Filter by impact level

### 4.4 Test Performance Metrics
- [ ] TTFT (Time to First Token) displays
- [ ] Tokens/sec displays
- [ ] Metrics update after each query

### 4.5 Test VCF Preview
- [ ] Navigate to VCF Preview tab in sidebar
- [ ] Click "Load Preview"
- [ ] Verify variant table displays

---

## Step 5: Test Target Hypothesis Workflow

### 5.1 Create Target Hypothesis
1. Ask: "What kinase genes have high-impact variants that could be drug targets?"
2. Navigate to "Target Hypotheses" tab
3. Click "Add Target"
4. Fill in:
   - Gene: EGFR
   - Protein: Epidermal growth factor receptor
   - Rationale: [from LLM response]
   - Confidence: high
   - Priority: 5
   - PDB IDs: 7SYE
5. Save

### 5.2 Verify Target Saved
```bash
cat data/targets/hypotheses.json
```

### 5.3 Update Target Status
- [ ] Change status from "hypothesis" to "validated"
- [ ] Verify change persists

### 5.4 Export for Phase 5
1. Click "Export for Phase 5"
2. Download JSON
3. Verify content:
```bash
cat data/targets/targets_for_phase5.json
```

**Expected JSON structure:**
```json
{
  "targets": [
    {
      "gene": "EGFR",
      "pdb_ids": ["7SYE"],
      "priority": 5,
      ...
    }
  ]
}
```

---

## Step 6: Test Portal Interface

### 6.1 Start Portal
```bash
./run.sh portal
```

**Expected:** Opens http://localhost:5001

### 6.2 Verify Portal Functions
- [ ] Service status shows (Milvus, vLLM, Attu)
- [ ] System metrics display (GPU, Memory)
- [ ] VCF Preview tab works
- [ ] Targets tab shows saved targets
- [ ] Export button works
- [ ] Console output streams during operations

---

## Step 7: End-to-End Demo Test

### 7.1 Full Demo Run
1. Start fresh (stop all services)
   ```bash
   ./run.sh stop
   ```

2. Start services
   ```bash
   ./run.sh start
   ```

3. Open Portal (http://localhost:5001)

4. Open Chat (http://localhost:8501)

5. Run demo queries:
   ```
   "What are the most promising drug targets in this sample based on variant impact?"
   "Tell me about EGFR variants for structure-based drug design"
   "What kinase mutations could be targeted with small molecules?"
   ```

6. Add 2-3 target hypotheses with PDB IDs

7. Export for Phase 5

8. Verify export file ready for handoff

### 7.2 Timing Verification

| Step | Target Time |
|------|-------------|
| Query response (TTFT) | < 2 seconds |
| Full response | < 30 seconds |
| VCF preview load | < 5 seconds |
| Target save | < 1 second |
| Export | < 1 second |

---

## Step 8: Prepare Demo Scenarios

### 8.1 Pre-load Demo Data
Create script to pre-populate targets:
```bash
python scripts/setup_demo_targets.py
```

### 8.2 Prepare Question Bank
Save known-good queries that produce compelling answers.

### 8.3 Pre-fetch Structures
Download PDB files for quick Phase 5 transition.

---

## Troubleshooting

### vLLM Won't Start
```bash
# Check logs
docker-compose logs vllm

# Common issues:
# - HF_TOKEN not set or invalid
# - Model license not accepted
# - Insufficient GPU memory (need ~80GB for 70B)
```

### Milvus Connection Failed
```bash
# Check container
docker-compose ps milvus

# Check logs
docker-compose logs milvus

# Restart
docker-compose restart milvus
```

### Slow Ingestion
```bash
# Use smaller batch size
python scripts/ingest_vcf.py --batch-size 500

# Check GPU memory during embedding
nvidia-smi
```

### Chat Errors
```bash
# Check all services running
./run.sh status

# Verify collection loaded
python -c "
from pymilvus import connections, Collection
connections.connect('default', host='localhost', port='19530')
coll = Collection('genomic_evidence')
coll.load()
print('Collection loaded')
"
```

---

## Success Criteria

Phase 4 is complete when:

- [ ] All services start successfully
- [ ] VCF fully ingested into Milvus (4M+ variants)
- [ ] Chat responds to queries with relevant variants
- [ ] Performance metrics display correctly
- [ ] Target hypotheses can be created/saved
- [ ] Export generates valid JSON for Phase 5
- [ ] Portal shows all status correctly
- [ ] End-to-end demo runs smoothly
- [ ] Response times meet targets

---

## Files Created in Phase 4

```
rag-chat-pipeline/
├── .env                          # Configuration
├── docker-compose.yml            # Services
├── requirements.txt              # Python deps
├── run.sh                        # CLI
├── config/
│   └── settings.py               # App settings
├── src/
│   ├── vcf_parser.py             # VCF parsing
│   ├── annotator.py              # VEP annotation
│   ├── embedder.py               # Text embeddings
│   ├── milvus_client.py          # Vector DB
│   ├── llm_client.py             # LLM interface
│   ├── rag_engine.py             # RAG orchestration
│   └── target_hypothesis.py      # Target management
├── app/
│   └── chat_ui.py                # Streamlit UI
├── portal/
│   ├── app/server.py             # Flask backend
│   ├── templates/index.html      # Portal UI
│   └── static/                   # CSS/JS
├── scripts/
│   ├── ingest_vcf.py             # Ingestion script
│   └── run_chat.py               # Chat launcher
└── data/
    └── targets/                  # Saved hypotheses
        ├── hypotheses.json
        └── targets_for_phase5.json
```

---

## Next Phase Handoff

When Phase 4 is complete, Phase 5 (Cryo-EM) receives:
- `targets_for_phase5.json` with prioritized targets and PDB IDs
- Ready for structure fetching and analysis
