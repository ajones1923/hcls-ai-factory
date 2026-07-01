# Precision Autoimmune Intelligence Agent -- API Reference

## Base URL

```
http://localhost:8532
```

The API port is configurable via the `AUTO_API_PORT` environment variable (default: 8532).

## Authentication

Authentication is optional. When `AUTO_API_KEY` is set to a non-empty value, all endpoints except `/`, `/health`, `/healthz`, and `/metrics` require the API key to be provided via one of:

- **Header**: `X-API-Key: <your-api-key>`
- **Query parameter**: `?api_key=<your-api-key>`

When `AUTO_API_KEY` is empty or unset, no authentication is required.

## Response Headers

All responses include the `X-Process-Time-Ms` header indicating the server-side processing time in milliseconds.

## Request Size Limit

The maximum request body size is controlled by `AUTO_MAX_REQUEST_SIZE_MB` (default: 50 MB). Requests exceeding this limit receive a `413` response.

---

## Endpoints

### GET /

Root endpoint. Returns service identity and status.

**Authentication**: Not required.

**Response**:

```json
{
  "service": "Precision Autoimmune Intelligence Agent",
  "version": "1.0.0",
  "status": "running",
  "ports": {
    "api": 8532,
    "ui": 8531
  }
}
```

**cURL**:

```bash
curl http://localhost:8532/
```

---

### GET /health

Returns detailed health status including Milvus connectivity, collection counts, embedder status, and LLM availability.

**Authentication**: Not required.

**Response**:

```json
{
  "status": "healthy",
  "service": "autoimmune-agent",
  "milvus_connected": true,
  "collections": 14,
  "total_vectors": 12450,
  "embedder_loaded": true,
  "llm_available": true,
  "uptime_seconds": 3621
}
```

**cURL**:

```bash
curl http://localhost:8532/health
```

---

### GET /healthz

Lightweight health probe for load balancers and the HCLS AI Factory landing page.

**Authentication**: Not required.

**Response**:

```json
{
  "status": "ok"
}
```

**cURL**:

```bash
curl http://localhost:8532/healthz
```

---

### GET /metrics

Prometheus-compatible metrics in text exposition format.

**Authentication**: Not required.

**Response** (Content-Type: `text/plain; version=0.0.4`):

```
# HELP autoimmune_agent_up Whether the agent is running
# TYPE autoimmune_agent_up gauge
autoimmune_agent_up 1
# HELP autoimmune_collection_vectors Number of vectors per collection
# TYPE autoimmune_collection_vectors gauge
autoimmune_collection_vectors{collection="autoimmune_clinical_documents"} 2340
autoimmune_collection_vectors{collection="autoimmune_patient_labs"} 1850
# HELP autoimmune_agent_uptime_seconds Agent uptime
# TYPE autoimmune_agent_uptime_seconds gauge
autoimmune_agent_uptime_seconds 3621
```

**cURL**:

```bash
curl http://localhost:8532/metrics
```

---

### POST /query

Full RAG query: retrieves evidence from Milvus collections and synthesizes a response using Claude.

**Authentication**: Required (when enabled).

**Request Body**:

| Field | Type | Required | Description |
|---|---|---|---|
| `question` | string | Yes | The clinical question to answer |
| `patient_id` | string | No | Filter results to a specific patient |
| `patient_context` | string | No | Additional clinical context to include in the prompt |
| `collections_filter` | string[] | No | Limit search to specific collection names |
| `top_k` | integer | No | Override the number of results per collection (default: 5) |

**Request Example**:

```json
{
  "question": "What is the significance of anti-dsDNA antibodies with falling complement C3/C4 levels?",
  "patient_id": "SMI-2022-44871",
  "patient_context": "34-year-old female with ANA 1:640 homogeneous pattern, proteinuria on urinalysis",
  "top_k": 8
}
```

**Response**:

```json
{
  "answer": "The combination of positive anti-dsDNA antibodies with declining complement C3/C4 levels is highly specific for active systemic lupus erythematosus (SLE) and strongly suggests immune-complex-mediated disease activity...",
  "evidence_count": 18,
  "collections_searched": 14,
  "search_time_ms": 142.3
}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "question": "What is the significance of anti-dsDNA antibodies with falling complement C3/C4 levels?",
    "patient_id": "SMI-2022-44871",
    "patient_context": "34-year-old female with ANA 1:640 homogeneous pattern",
    "top_k": 8
  }'
```

---

### POST /query/stream

Streaming RAG query. Retrieves evidence and streams the Claude response as Server-Sent Events (SSE).

**Authentication**: Required (when enabled).

**Request Body**: Same schema as `POST /query`.

**Request Example**:

```json
{
  "question": "Compare TNF inhibitors vs IL-6 receptor antagonists for rheumatoid arthritis",
  "collections_filter": ["autoimmune_biologic_therapies", "autoimmune_clinical_trials"]
}
```

**Response** (Content-Type: `text/event-stream`):

Each event contains a JSON object with a `text` field. The stream terminates with `[DONE]`.

```
data: {"text": "When comparing TNF inhibitors"}

data: {"text": " and IL-6 receptor antagonists"}

data: {"text": " for rheumatoid arthritis..."}

data: [DONE]
```

If an error occurs during streaming:

```
data: {"error": "LLM streaming failed: connection timeout"}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/query/stream \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -N \
  -d '{
    "question": "Compare TNF inhibitors vs IL-6 receptor antagonists for rheumatoid arthritis"
  }'
```

---

### POST /search

Evidence-only search across Milvus collections without LLM synthesis. Returns raw search hits with scores and metadata.

**Authentication**: Required (when enabled).

**Request Body**:

| Field | Type | Required | Description |
|---|---|---|---|
| `question` | string | Yes | The search query |
| `patient_id` | string | No | Filter results to a specific patient |
| `collections_filter` | string[] | No | Limit search to specific collection names |
| `top_k` | integer | No | Override the number of results per collection (default: 5) |

**Request Example**:

```json
{
  "question": "HLA-B27 ankylosing spondylitis sacroiliitis",
  "collections_filter": ["autoimmune_hla_associations", "autoimmune_disease_criteria"],
  "top_k": 10
}
```

**Response**:

```json
{
  "hits": [
    {
      "collection": "autoimmune_hla_associations",
      "id": "hla-b27-as-001",
      "score": 0.923,
      "text": "HLA-B*27:05 is the strongest genetic risk factor for ankylosing spondylitis with an odds ratio of 87.4...",
      "relevance": "high",
      "metadata": {
        "allele": "B*27:05",
        "disease": "ankylosing_spondylitis",
        "odds_ratio": 87.4,
        "pmid": "29123456"
      }
    }
  ],
  "total": 12
}
```

The `text` field is truncated to 1000 characters. The `relevance` field is one of `high`, `medium`, or `low` based on the configurable score thresholds.

**cURL**:

```bash
curl -X POST http://localhost:8532/search \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "question": "HLA-B27 ankylosing spondylitis sacroiliitis",
    "top_k": 10
  }'
```

---

### POST /analyze

Runs the full autoimmune analysis pipeline on a patient profile. Performs autoantibody interpretation, HLA association analysis, disease activity scoring, flare prediction, and biologic therapy recommendation.

**Authentication**: Required (when enabled).

**Request Body**: An `AutoimmunePatientProfile` object. At least one of `autoantibody_panel`, `hla_profile`, or `biomarkers` must be provided.

| Field | Type | Required | Description |
|---|---|---|---|
| `patient_id` | string | Yes | Unique patient identifier |
| `age` | integer | Yes | Patient age (0-150) |
| `sex` | string | Yes | `M` or `F` |
| `autoantibody_panel` | object | No | Autoantibody panel results (see below) |
| `hla_profile` | object | No | HLA typing results (see below) |
| `biomarkers` | object | No | Key-value map of biomarker names to numeric values |
| `genotypes` | object | No | Key-value map of gene names to genotype strings |
| `diagnosed_conditions` | string[] | No | List of diagnosed autoimmune diseases (enum values) |
| `current_medications` | string[] | No | Current medication list |
| `symptom_duration_months` | integer | No | Duration of symptoms in months |
| `family_history` | string[] | No | Family history entries |

**Autoantibody Panel** (`autoantibody_panel`):

| Field | Type | Required | Description |
|---|---|---|---|
| `patient_id` | string | Yes | Patient identifier |
| `collection_date` | string | Yes | ISO-8601 date |
| `results` | array | No | List of autoantibody results |

Each result in `results`:

| Field | Type | Required | Description |
|---|---|---|---|
| `antibody` | string | Yes | Autoantibody name (e.g., `ANA`, `anti-dsDNA`, `RF`) |
| `value` | float | Yes | Measured value |
| `unit` | string | No | Unit of measurement |
| `reference_range` | string | No | Normal reference range |
| `positive` | boolean | No | Whether the result is positive |
| `titer` | string | No | Titer if applicable (e.g., `1:320`) |
| `pattern` | string | No | ANA staining pattern (e.g., `homogeneous`, `speckled`) |

**HLA Profile** (`hla_profile`):

| Field | Type | Required | Description |
|---|---|---|---|
| `hla_a` | string[] | No | HLA-A alleles |
| `hla_b` | string[] | No | HLA-B alleles |
| `hla_c` | string[] | No | HLA-C alleles |
| `hla_drb1` | string[] | No | HLA-DRB1 alleles |
| `hla_dqb1` | string[] | No | HLA-DQB1 alleles |

**Request Example**:

```json
{
  "patient_id": "DPA-2019-33102",
  "age": 45,
  "sex": "M",
  "hla_profile": {
    "hla_b": ["B*27:05"]
  },
  "biomarkers": {
    "CRP": 28.5,
    "ESR": 48
  },
  "diagnosed_conditions": ["ankylosing_spondylitis"],
  "current_medications": ["naproxen", "physical therapy"],
  "symptom_duration_months": 72
}
```

**Response** (an `AutoimmuneAnalysisResult`):

```json
{
  "patient_id": "DPA-2019-33102",
  "disease_activity_scores": [
    {
      "disease": "ankylosing_spondylitis",
      "score_name": "BASDAI",
      "score_value": 28.5,
      "level": "high",
      "components": {"CRP": 28.5},
      "thresholds": {"remission": 1.0, "low": 2.0, "moderate": 4.0, "high": 6.0}
    }
  ],
  "flare_predictions": [
    {
      "disease": "ankylosing_spondylitis",
      "current_activity": "moderate",
      "predicted_risk": "high",
      "risk_score": 0.6,
      "contributing_factors": ["Elevated CRP: 28.5"],
      "protective_factors": [],
      "recommended_monitoring": ["Repeat CRP in 2-4 weeks", "Repeat ESR in 2-4 weeks"],
      "time_horizon_days": 90
    }
  ],
  "hla_associations": [
    {
      "allele": "B*27:05",
      "disease": "ankylosing_spondylitis",
      "odds_ratio": 87.4,
      "pmid": "29123456",
      "note": "Strongest genetic risk factor for AS"
    }
  ],
  "biologic_recommendations": [
    {
      "drug_name": "Adalimumab",
      "drug_class": "TNF inhibitor",
      "mechanism": "Binds and neutralizes TNF-alpha",
      "indicated_diseases": ["ankylosing_spondylitis"],
      "pgx_considerations": ["FCGR3A V/F polymorphism may affect response"],
      "contraindications": ["Active TB", "Severe heart failure"],
      "monitoring_requirements": ["TB screening", "Hepatitis B/C", "CBC every 3 months"],
      "efficacy_evidence": ""
    }
  ],
  "critical_alerts": [
    "HIGH DISEASE ACTIVITY: ankylosing_spondylitis (BASDAI = 28.5)",
    "STRONG HLA ASSOCIATION: B*27:05 -> ankylosing_spondylitis (OR=87.4)"
  ],
  "cross_agent_findings": []
}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/analyze \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "patient_id": "DPA-2019-33102",
    "age": 45,
    "sex": "M",
    "hla_profile": {"hla_b": ["B*27:05"]},
    "biomarkers": {"CRP": 28.5, "ESR": 48},
    "diagnosed_conditions": ["ankylosing_spondylitis"]
  }'
```

---

### POST /differential

Generates a differential diagnosis from autoantibody and HLA data. Uses the diagnostic engine to evaluate positive antibodies against known disease associations.

**Authentication**: Required (when enabled).

**Request Body**:

| Field | Type | Required | Description |
|---|---|---|---|
| `positive_antibodies` | string[] | Yes | List of positive autoantibody names |
| `hla_alleles` | string[] | No | List of HLA alleles |
| `symptoms` | string[] | No | List of presenting symptoms |

**Request Example**:

```json
{
  "positive_antibodies": ["ANA", "anti-dsDNA", "anti-Smith"],
  "hla_alleles": ["DRB1*03:01", "DRB1*15:01"],
  "symptoms": ["malar rash", "arthralgia", "proteinuria", "photosensitivity"]
}
```

**Response**:

```json
{
  "differential": [
    {
      "disease": "systemic_lupus_erythematosus",
      "confidence": 0.92,
      "supporting_antibodies": ["ANA", "anti-dsDNA", "anti-Smith"],
      "supporting_hla": ["DRB1*03:01"],
      "matching_symptoms": ["malar rash", "arthralgia", "proteinuria", "photosensitivity"],
      "classification_criteria_met": 8,
      "notes": "Meets 2019 EULAR/ACR SLE classification criteria"
    }
  ]
}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/differential \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "positive_antibodies": ["ANA", "anti-dsDNA", "anti-Smith"],
    "hla_alleles": ["DRB1*03:01"],
    "symptoms": ["malar rash", "arthralgia", "proteinuria"]
  }'
```

---

### POST /ingest/upload

Uploads and ingests a clinical PDF document into the `autoimmune_clinical_documents` collection. The document is chunked, embedded, and stored in Milvus.

**Authentication**: Required (when enabled).

**Content-Type**: `multipart/form-data`

**Form Fields**:

| Field | Type | Required | Description |
|---|---|---|---|
| `file` | file | Yes | PDF file to ingest |
| `patient_id` | string | Yes | Patient identifier to associate with the document |

Only PDF files are accepted.

**Response**:

```json
{
  "status": "success",
  "filename": "2024-03-15_Rheumatology_Progress_Note.pdf",
  "patient_id": "SMI-2022-44871",
  "chunks_ingested": 12,
  "doc_type": "progress_note",
  "specialty": "rheumatology"
}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/ingest/upload \
  -H "X-API-Key: your-api-key" \
  -F "file=@/path/to/clinical_note.pdf" \
  -F "patient_id=SMI-2022-44871"
```

---

### POST /ingest/demo-data

Ingests all demo patient data from the `demo_data/` directory. Processes all PDF files across all patient subdirectories and stores the chunks in Milvus.

**Authentication**: Required (when enabled).

**Request Body**: None.

**Response**:

```json
{
  "status": "success",
  "patients": {
    "sarah_mitchell": 156,
    "maya_rodriguez": 134,
    "david_park": 98,
    "linda_chen": 112,
    "rachel_thompson": 105
  },
  "total_chunks": 605
}
```

**cURL**:

```bash
curl -X POST http://localhost:8532/ingest/demo-data \
  -H "X-API-Key: your-api-key"
```

---

### GET /collections

Lists all Milvus collections with their vector counts.

**Authentication**: Required (when enabled).

**Response**:

```json
{
  "collections": [
    {"name": "autoimmune_autoantibody_panels", "count": 245},
    {"name": "autoimmune_biologic_therapies", "count": 180},
    {"name": "autoimmune_clinical_documents", "count": 2340},
    {"name": "autoimmune_clinical_trials", "count": 420},
    {"name": "autoimmune_cross_disease", "count": 85},
    {"name": "autoimmune_disease_activity", "count": 310},
    {"name": "autoimmune_disease_criteria", "count": 156},
    {"name": "autoimmune_flare_patterns", "count": 198},
    {"name": "autoimmune_hla_associations", "count": 275},
    {"name": "autoimmune_literature", "count": 1850},
    {"name": "autoimmune_patient_labs", "count": 1920},
    {"name": "autoimmune_patient_timelines", "count": 450},
    {"name": "autoimmune_pgx_rules", "count": 120},
    {"name": "genomic_evidence", "count": 3500}
  ],
  "total_collections": 14,
  "total_vectors": 12049
}
```

**cURL**:

```bash
curl http://localhost:8532/collections \
  -H "X-API-Key: your-api-key"
```

---

### POST /collections/create

Creates all 14 Milvus collections with the configured schema and embedding dimension. Optionally drops existing collections first.

**Authentication**: Required (when enabled).

**Query Parameters**:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `drop_existing` | boolean | `false` | If `true`, drops and recreates existing collections |

**Response**:

```json
{
  "status": "success",
  "collections_created": [
    "autoimmune_clinical_documents",
    "autoimmune_patient_labs",
    "autoimmune_autoantibody_panels",
    "autoimmune_hla_associations",
    "autoimmune_disease_criteria",
    "autoimmune_disease_activity",
    "autoimmune_flare_patterns",
    "autoimmune_biologic_therapies",
    "autoimmune_pgx_rules",
    "autoimmune_clinical_trials",
    "autoimmune_literature",
    "autoimmune_patient_timelines",
    "autoimmune_cross_disease",
    "genomic_evidence"
  ],
  "count": 14
}
```

**cURL**:

```bash
# Create collections (skip existing)
curl -X POST http://localhost:8532/collections/create \
  -H "X-API-Key: your-api-key"

# Drop and recreate all collections
curl -X POST "http://localhost:8532/collections/create?drop_existing=true" \
  -H "X-API-Key: your-api-key"
```

---

### POST /export

Exports an analysis report for a patient in the specified format.

**Authentication**: Required (when enabled).

**Request Body**:

| Field | Type | Required | Default | Description |
|---|---|---|---|---|
| `patient_id` | string | Yes | -- | Patient identifier |
| `format` | string | No | `markdown` | Export format: `markdown`, `fhir`, or `pdf` |
| `query_answer` | string | No | -- | Optional query answer to include in the report |

**Request Example** (Markdown):

```json
{
  "patient_id": "SMI-2022-44871",
  "format": "markdown",
  "query_answer": "Based on the analysis, Sarah Mitchell meets 2019 EULAR/ACR criteria for SLE..."
}
```

**Response** (Markdown):

```json
{
  "format": "markdown",
  "content": "# Autoimmune Analysis Report\n## Patient: SMI-2022-44871\n..."
}
```

**Request Example** (FHIR R4):

```json
{
  "patient_id": "SMI-2022-44871",
  "format": "fhir"
}
```

**Response** (FHIR R4):

```json
{
  "format": "fhir",
  "content": {
    "resourceType": "DiagnosticReport",
    "status": "final",
    "code": {
      "coding": [{
        "system": "http://loinc.org",
        "code": "11526-1",
        "display": "Pathology study"
      }]
    },
    "subject": {
      "reference": "Patient/SMI-2022-44871"
    }
  }
}
```

**Request Example** (PDF):

```json
{
  "patient_id": "SMI-2022-44871",
  "format": "pdf"
}
```

**Response** (PDF):

```json
{
  "format": "pdf",
  "content_base64": "JVBERi0xLjQKMSAwIG9iago8PAovVHlwZSAvQ2F0YWxvZw..."
}
```

**cURL**:

```bash
# Markdown export
curl -X POST http://localhost:8532/export \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{"patient_id": "SMI-2022-44871", "format": "markdown"}'

# FHIR R4 export
curl -X POST http://localhost:8532/export \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{"patient_id": "SMI-2022-44871", "format": "fhir"}'

# PDF export (base64-encoded in response)
curl -X POST http://localhost:8532/export \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{"patient_id": "SMI-2022-44871", "format": "pdf"}'
```

---

## Autoimmune Disease Enum Values

The following string values are valid for the `diagnosed_conditions` field and other disease-related parameters:

| Enum Value | Disease |
|---|---|
| `rheumatoid_arthritis` | Rheumatoid Arthritis |
| `systemic_lupus_erythematosus` | Systemic Lupus Erythematosus |
| `multiple_sclerosis` | Multiple Sclerosis |
| `type_1_diabetes` | Type 1 Diabetes |
| `inflammatory_bowel_disease` | Inflammatory Bowel Disease |
| `psoriasis` | Psoriasis |
| `ankylosing_spondylitis` | Ankylosing Spondylitis |
| `sjogrens_syndrome` | Sjogren's Syndrome |
| `systemic_sclerosis` | Systemic Sclerosis |
| `myasthenia_gravis` | Myasthenia Gravis |
| `celiac_disease` | Celiac Disease |
| `graves_disease` | Graves' Disease |
| `hashimoto_thyroiditis` | Hashimoto's Thyroiditis |

---

## Error Responses

### 400 Bad Request

Returned when the request contains invalid parameters (e.g., unsupported export format, non-PDF file upload).

```json
{
  "detail": "Only PDF files are supported"
}
```

### 401 Unauthorized

Returned when `AUTO_API_KEY` is configured and the request does not include a valid API key.

```json
{
  "detail": "Invalid or missing API key"
}
```

### 413 Request Entity Too Large

Returned when the request body exceeds `AUTO_MAX_REQUEST_SIZE_MB` (default: 50 MB).

```json
{
  "detail": "Request too large"
}
```

### 422 Unprocessable Entity

Returned when the request body fails Pydantic validation or when document processing yields no content.

```json
{
  "detail": [
    {
      "loc": ["body", "patient_id"],
      "msg": "Field required",
      "type": "missing"
    }
  ]
}
```

Also returned when a PDF upload yields no extractable text:

```json
{
  "detail": "No text extracted from PDF"
}
```

The `POST /analyze` endpoint returns a 422 if the patient profile is missing all of `autoantibody_panel`, `hla_profile`, and `biomarkers`:

```json
{
  "detail": [
    {
      "loc": ["body"],
      "msg": "Value error, At least one of autoantibody_panel, hla_profile, or biomarkers must be provided",
      "type": "value_error"
    }
  ]
}
```

### 503 Service Unavailable

Returned when a required internal service (RAG engine, agent, document processor, diagnostic engine, collection manager) has not been initialized.

```json
{
  "detail": "RAG engine not initialized"
}
```

### 500 Internal Server Error

Returned for unexpected server-side errors.

```json
{
  "detail": "Internal server error"
}
```
