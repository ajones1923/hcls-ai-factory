#!/bin/bash
# ============================================================================
# HCLS AI Factory — Full Demo Script
#
# Demonstrates the complete precision medicine platform across all 11
# intelligence agents using a pediatric oncology persona:
#   8-year-old female with B-cell Acute Lymphoblastic Leukemia (ALL)
#
# Usage:
#   ./scripts/run_full_demo.sh              # Run full demo
#   ./scripts/run_full_demo.sh --health     # Health checks only
#   ./scripts/run_full_demo.sh --quick      # Skip LLM-dependent queries
#
# Requires: curl, jq (optional, for pretty output)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DEMO_DIR="${PROJECT_DIR}/demo"
REQUESTS_DIR="${DEMO_DIR}/requests"

# ── Colors ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
BOLD='\033[1m'
DIM='\033[2m'
RESET='\033[0m'

# ── Configuration ───────────────────────────────────────────────────────────
CURL_TIMEOUT=15
HEALTH_ONLY=false
QUICK_MODE=false
PASSED=0
FAILED=0
SKIPPED=0

for arg in "$@"; do
    case "$arg" in
        --health) HEALTH_ONLY=true ;;
        --quick)  QUICK_MODE=true ;;
    esac
done

# ── Agent Registry ──────────────────────────────────────────────────────────
# Format: ID|PORT|NAME
declare -a AGENTS=(
    "biomarker|8529|Precision Biomarker"
    "oncology|8527|Precision Oncology"
    "cart|8522|CAR-T Intelligence"
    "imaging|8524|Imaging Intelligence"
    "autoimmune|8532|Precision Autoimmune"
    "pharmacogenomics|8107|Pharmacogenomics"
    "cardiology|8126|Cardiology Intelligence"
    "clinical-trial|8538|Clinical Trial Intelligence"
    "rare-disease|8134|Rare Disease Diagnostic"
    "neurology|8528|Neurology Intelligence"
    "single-cell|8540|Single-Cell Intelligence"
)

# ── Utility Functions ───────────────────────────────────────────────────────

banner() {
    echo ""
    echo -e "${BOLD}${BLUE}╔══════════════════════════════════════════════════════════════╗${RESET}"
    echo -e "${BOLD}${BLUE}║${RESET}  ${BOLD}$1${RESET}"
    echo -e "${BOLD}${BLUE}╚══════════════════════════════════════════════════════════════╝${RESET}"
}

section() {
    echo ""
    echo -e "${CYAN}┌──────────────────────────────────────────────────────────────┐${RESET}"
    echo -e "${CYAN}│${RESET}  ${BOLD}$1${RESET}"
    echo -e "${CYAN}└──────────────────────────────────────────────────────────────┘${RESET}"
}

step() {
    echo -e "  ${MAGENTA}▸${RESET} $1"
}

ok() {
    echo -e "    ${GREEN}✓${RESET} $1"
    PASSED=$((PASSED + 1))
}

fail() {
    echo -e "    ${RED}✗${RESET} $1"
    FAILED=$((FAILED + 1))
}

skip() {
    echo -e "    ${YELLOW}⊘${RESET} $1 ${DIM}(skipped)${RESET}"
    SKIPPED=$((SKIPPED + 1))
}

pretty_json() {
    if command -v jq &>/dev/null; then
        jq '.' 2>/dev/null || cat
    else
        cat
    fi
}

# Perform a curl call; return 0 on success (HTTP 2xx/3xx), 1 on failure.
# Prints the response body.
api_call() {
    local method="$1"
    local url="$2"
    local data="${3:-}"
    local response code body

    if [ "$method" = "GET" ]; then
        response=$(curl -s -w "\n%{http_code}" --connect-timeout "$CURL_TIMEOUT" \
                        --max-time "$CURL_TIMEOUT" "$url" 2>/dev/null) || { fail "Connection refused: $url"; return 1; }
    else
        response=$(curl -s -w "\n%{http_code}" --connect-timeout "$CURL_TIMEOUT" \
                        --max-time "$CURL_TIMEOUT" -X POST \
                        -H "Content-Type: application/json" \
                        -d "$data" "$url" 2>/dev/null) || { fail "Connection refused: $url"; return 1; }
    fi

    code=$(echo "$response" | tail -1)
    body=$(echo "$response" | sed '$d')

    if [ -n "$code" ] && [ "$code" -ge 200 ] 2>/dev/null && [ "$code" -lt 400 ] 2>/dev/null; then
        echo "$body"
        return 0
    else
        echo "$body"
        return 1
    fi
}

# ============================================================================
# PHASE 0: Banner
# ============================================================================

clear 2>/dev/null || true
echo ""
echo -e "${BOLD}${BLUE}"
echo "    ╦ ╦╔═╗╦  ╔═╗  ╔═╗╦  ╔═╗╔═╗╔═╗╔╦╗╔═╗╦═╗╦ ╦"
echo "    ╠═╣║  ║  ╚═╗  ╠═╣║  ╠╣ ╠═╣║   ║ ║ ║╠╦╝╚╦╝"
echo "    ╩ ╩╚═╝╩═╝╚═╝  ╩ ╩╩  ╚  ╩ ╩╚═╝ ╩ ╚═╝╩╚═ ╩ "
echo -e "${RESET}"
echo -e "  ${BOLD}Precision Medicine Platform — Full Demo${RESET}"
echo -e "  ${DIM}Patient DNA to Drug Candidates in <5 hours${RESET}"
echo ""
echo -e "  ${YELLOW}Demo Persona:${RESET} ${BOLD}Pediatric Oncology${RESET}"
echo -e "  ${DIM}  8-year-old female | B-cell ALL | PEDS-ALL-001${RESET}"
echo -e "  ${DIM}  $(date '+%Y-%m-%d %H:%M:%S')${RESET}"

# ============================================================================
# PHASE 1: Health Checks
# ============================================================================

banner "PHASE 1: Service Health Checks"

echo ""
echo -e "  ${DIM}Checking all 11 intelligence agents...${RESET}"
echo ""

HEALTHY_COUNT=0
UNHEALTHY_LIST=""

for entry in "${AGENTS[@]}"; do
    id=$(echo "$entry" | cut -d'|' -f1)
    port=$(echo "$entry" | cut -d'|' -f2)
    name=$(echo "$entry" | cut -d'|' -f3)

    if result=$(api_call GET "http://localhost:${port}/health" 2>/dev/null); then
        printf "  ${GREEN}●${RESET} %-30s :%-6s ${GREEN}HEALTHY${RESET}\n" "$name" "$port"
        HEALTHY_COUNT=$((HEALTHY_COUNT + 1))
    else
        printf "  ${RED}○${RESET} %-30s :%-6s ${RED}DOWN${RESET}\n" "$name" "$port"
        UNHEALTHY_LIST="${UNHEALTHY_LIST} ${id}"
    fi
done

echo ""
echo -e "  ────────────────────────────────────────────"
echo -e "  ${BOLD}${HEALTHY_COUNT}/${#AGENTS[@]} agents healthy${RESET}"

if [ -n "$UNHEALTHY_LIST" ]; then
    echo -e "  ${RED}Down:${RESET}${UNHEALTHY_LIST}"
fi

# Also check core infrastructure
section "Core Infrastructure"

for svc_info in "Milvus|19530|TCP" "Landing Page|8080|/health" "Genomics Portal|5000|/health" "RAG/Chat API|5001|/health"; do
    name=$(echo "$svc_info" | cut -d'|' -f1)
    port=$(echo "$svc_info" | cut -d'|' -f2)
    path=$(echo "$svc_info" | cut -d'|' -f3)

    if [ "$path" = "TCP" ]; then
        if timeout 2 bash -c "echo > /dev/tcp/localhost/$port" 2>/dev/null; then
            printf "  ${GREEN}●${RESET} %-30s :%-6s ${GREEN}HEALTHY${RESET}\n" "$name" "$port"
        else
            printf "  ${RED}○${RESET} %-30s :%-6s ${RED}DOWN${RESET}\n" "$name" "$port"
        fi
    else
        if api_call GET "http://localhost:${port}${path}" &>/dev/null; then
            printf "  ${GREEN}●${RESET} %-30s :%-6s ${GREEN}HEALTHY${RESET}\n" "$name" "$port"
        else
            printf "  ${RED}○${RESET} %-30s :%-6s ${RED}DOWN${RESET}\n" "$name" "$port"
        fi
    fi
done

if [ "$HEALTH_ONLY" = true ]; then
    echo ""
    echo -e "  ${DIM}--health flag set; skipping demo flow.${RESET}"
    exit 0
fi

# ============================================================================
# PHASE 2: Clinical Demo Flow
# ============================================================================

banner "PHASE 2: Clinical Intelligence Demo"
echo -e "  ${DIM}Executing cross-agent workflow for patient PEDS-ALL-001${RESET}"

# ── 2a. Cardiology: ASCVD Risk (baseline cardiac assessment) ────────────────
section "Step 1: Cardiology — Baseline Cardiac Risk (ASCVD)"
step "Calculating 10-year ASCVD risk for chemotherapy cardiotoxicity monitoring"
step "POST http://localhost:8126/v1/cardio/risk/ascvd"

ASCVD_PAYLOAD='{
    "age": 45,
    "sex": "female",
    "race": "white",
    "total_cholesterol": 210,
    "hdl_cholesterol": 55,
    "systolic_bp": 128,
    "bp_treatment": false,
    "diabetes": false,
    "smoker": false
}'

if [ -f "${REQUESTS_DIR}/cardiology_risk.json" ]; then
    ASCVD_PAYLOAD=$(cat "${REQUESTS_DIR}/cardiology_risk.json")
fi

echo ""
if result=$(api_call POST "http://localhost:8126/v1/cardio/risk/ascvd" "$ASCVD_PAYLOAD"); then
    echo "$result" | pretty_json
    ok "ASCVD risk calculated (cardio-oncology baseline for patient's parent/guardian demo)"
else
    fail "Cardiology ASCVD endpoint unavailable"
fi

# ── 2b. Neurology: NIHSS Scale Calculation ──────────────────────────────────
section "Step 2: Neurology — NIHSS Scale Calculation"
step "Demonstrating validated neurological scale calculator"
step "POST http://localhost:8528/v1/neuro/scale/calculate"

NIHSS_PAYLOAD='{
    "scale_name": "nihss",
    "items": {
        "1a_loc": 0,
        "1b_loc_questions": 0,
        "1c_loc_commands": 0,
        "2_gaze": 0,
        "3_visual": 0,
        "4_facial_palsy": 0,
        "5a_motor_left_arm": 0,
        "5b_motor_right_arm": 0,
        "6a_motor_left_leg": 0,
        "6b_motor_right_leg": 0,
        "7_limb_ataxia": 0,
        "8_sensory": 0,
        "9_language": 0,
        "10_dysarthria": 0,
        "11_extinction": 0
    }
}'

if [ -f "${REQUESTS_DIR}/neurology_nihss.json" ]; then
    NIHSS_PAYLOAD=$(cat "${REQUESTS_DIR}/neurology_nihss.json")
fi

echo ""
if result=$(api_call POST "http://localhost:8528/v1/neuro/scale/calculate" "$NIHSS_PAYLOAD"); then
    echo "$result" | pretty_json
    ok "NIHSS scale calculated — baseline for CNS leukemia monitoring"
else
    fail "Neurology scale endpoint unavailable"
fi

# ── 2c. Rare Disease: Phenotype-Driven Diagnosis ───────────────────────────
section "Step 3: Rare Disease — Phenotype-Driven Differential Diagnosis"
step "Submitting HPO phenotypes + genomic variants for diagnostic workup"
step "POST http://localhost:8134/v1/diagnostic/diagnose"

DIAGNOSE_PAYLOAD='{
    "phenotypes": [
        {"id": "HP:0001875", "label": "Neutropenia", "onset": "infantile", "severity": "severe"},
        {"id": "HP:0001903", "label": "Anemia", "onset": "infantile"},
        {"id": "HP:0001744", "label": "Splenomegaly"},
        {"id": "HP:0002240", "label": "Hepatomegaly"},
        {"id": "HP:0001973", "label": "Autoimmune thrombocytopenia"}
    ],
    "variants": [
        {"gene": "PAX5", "variant": "c.547C>T", "zygosity": "heterozygous"},
        {"gene": "IKZF1", "variant": "del(7p12.2)", "zygosity": "heterozygous"}
    ],
    "age_years": 8,
    "sex": "female",
    "family_history": "No family history of hematologic malignancies",
    "clinical_notes": "Pediatric B-cell ALL with ETV6-RUNX1 fusion. Presenting with fatigue, pallor, petechiae, and bone pain.",
    "max_diagnoses": 5
}'

if [ -f "${REQUESTS_DIR}/rare_disease_diagnose.json" ]; then
    DIAGNOSE_PAYLOAD=$(cat "${REQUESTS_DIR}/rare_disease_diagnose.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Rare disease diagnosis (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8134/v1/diagnostic/diagnose" "$DIAGNOSE_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Rare disease differential generated"
    else
        fail "Rare disease diagnostic endpoint unavailable"
    fi
fi

# ── 2d. Clinical Trial: Patient-Trial Matching ─────────────────────────────
section "Step 4: Clinical Trial Intelligence — Patient-Trial Matching"
step "Matching pediatric ALL patient to eligible clinical trials"
step "POST http://localhost:8538/v1/trial/match"

MATCH_PAYLOAD='{
    "age": 8,
    "sex": "female",
    "diagnosis": "B-cell Acute Lymphoblastic Leukemia",
    "biomarkers": ["CD19+", "CD22+", "ETV6-RUNX1 fusion", "MRD negative"],
    "medications": ["Daunorubicin", "Vincristine", "PEG-asparaginase", "Dexamethasone"],
    "genomic_variants": ["PAX5 c.547C>T", "IKZF1 del(7p12.2)"],
    "comorbidities": [],
    "therapeutic_area": "pediatric_oncology",
    "max_results": 5
}'

if [ -f "${REQUESTS_DIR}/trial_match.json" ]; then
    MATCH_PAYLOAD=$(cat "${REQUESTS_DIR}/trial_match.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Trial matching (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8538/v1/trial/match" "$MATCH_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Patient-trial matching complete"
    else
        fail "Clinical trial match endpoint unavailable"
    fi
fi

# ============================================================================
# PHASE 3: Per-Agent RAG Queries
# ============================================================================

banner "PHASE 3: Intelligence Agent Queries"
echo -e "  ${DIM}Querying each agent with pediatric ALL context${RESET}"

# ── Precision Biomarker ────────────────────────────────────────────────────
section "Precision Biomarker Agent (:8529)"
step "POST http://localhost:8529/v1/query"

BIOMARKER_PAYLOAD='{
    "query": "What biomarkers are critical for monitoring minimal residual disease in pediatric B-cell ALL with ETV6-RUNX1 fusion?",
    "collections": null,
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/biomarker_query.json" ]; then
    BIOMARKER_PAYLOAD=$(cat "${REQUESTS_DIR}/biomarker_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Biomarker query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8529/v1/query" "$BIOMARKER_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Biomarker query complete"
    else
        fail "Biomarker agent query unavailable"
    fi
fi

# ── Precision Oncology ─────────────────────────────────────────────────────
section "Precision Oncology Agent (:8527)"
step "POST http://localhost:8527/query"

ONCOLOGY_PAYLOAD='{
    "query": "Treatment options for pediatric B-cell ALL with ETV6-RUNX1 fusion, PAX5 mutation, and IKZF1 deletion. Standard risk vs high risk stratification.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/oncology_query.json" ]; then
    ONCOLOGY_PAYLOAD=$(cat "${REQUESTS_DIR}/oncology_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Oncology query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8527/query" "$ONCOLOGY_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Oncology query complete"
    else
        fail "Oncology agent query unavailable"
    fi
fi

# ── CAR-T Intelligence ─────────────────────────────────────────────────────
section "CAR-T Intelligence Agent (:8522)"
step "POST http://localhost:8522/query"

CART_PAYLOAD='{
    "query": "CAR-T therapy options for relapsed/refractory pediatric B-cell ALL targeting CD19. Compare tisagenlecleucel vs brexucabtagene autoleucel efficacy and CRS rates.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/cart_query.json" ]; then
    CART_PAYLOAD=$(cat "${REQUESTS_DIR}/cart_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "CAR-T query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8522/query" "$CART_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "CAR-T query complete"
    else
        fail "CAR-T agent query unavailable"
    fi
fi

# ── Imaging Intelligence ───────────────────────────────────────────────────
section "Imaging Intelligence Agent (:8524)"
step "POST http://localhost:8524/query"

IMAGING_PAYLOAD='{
    "query": "Imaging protocol for initial staging and CNS evaluation in pediatric B-cell ALL. Brain MRI findings and chest CT considerations.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/imaging_query.json" ]; then
    IMAGING_PAYLOAD=$(cat "${REQUESTS_DIR}/imaging_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Imaging query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8524/query" "$IMAGING_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Imaging query complete"
    else
        fail "Imaging agent query unavailable"
    fi
fi

# ── Precision Autoimmune ───────────────────────────────────────────────────
section "Precision Autoimmune Agent (:8532)"
step "POST http://localhost:8532/query"

AUTOIMMUNE_PAYLOAD='{
    "query": "Autoimmune complications in pediatric ALL treatment: immune thrombocytopenia, autoimmune hemolytic anemia, and post-transplant autoimmune syndromes.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/autoimmune_query.json" ]; then
    AUTOIMMUNE_PAYLOAD=$(cat "${REQUESTS_DIR}/autoimmune_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Autoimmune query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8532/query" "$AUTOIMMUNE_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Autoimmune query complete"
    else
        fail "Autoimmune agent query unavailable"
    fi
fi

# ── Pharmacogenomics ───────────────────────────────────────────────────────
section "Pharmacogenomics Agent (:8107)"
step "POST http://localhost:8107/query"

PGX_PAYLOAD='{
    "query": "Pharmacogenomic considerations for 6-mercaptopurine and methotrexate in pediatric ALL. TPMT and NUDT15 genotype impact on thiopurine dosing.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/pharmacogenomics_query.json" ]; then
    PGX_PAYLOAD=$(cat "${REQUESTS_DIR}/pharmacogenomics_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Pharmacogenomics query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8107/query" "$PGX_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Pharmacogenomics query complete"
    else
        fail "Pharmacogenomics agent query unavailable"
    fi
fi

# ── Cardiology: RAG Query ──────────────────────────────────────────────────
section "Cardiology Intelligence Agent (:8126) — RAG Query"
step "POST http://localhost:8126/v1/cardio/query"

CARDIO_QUERY_PAYLOAD='{
    "query": "Anthracycline cardiotoxicity monitoring in pediatric oncology. Daunorubicin cumulative dose limits and echocardiographic surveillance guidelines.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/cardiology_query.json" ]; then
    CARDIO_QUERY_PAYLOAD=$(cat "${REQUESTS_DIR}/cardiology_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Cardiology RAG query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8126/v1/cardio/query" "$CARDIO_QUERY_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Cardiology RAG query complete"
    else
        fail "Cardiology agent query unavailable"
    fi
fi

# ── Clinical Trial: RAG Query ──────────────────────────────────────────────
section "Clinical Trial Intelligence Agent (:8538) — RAG Query"
step "POST http://localhost:8538/v1/trial/query"

TRIAL_QUERY_PAYLOAD='{
    "query": "Current phase III clinical trials for pediatric B-cell ALL with novel immunotherapy combinations. Blinatumomab and inotuzumab ozogamicin frontline integration.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/trial_query.json" ]; then
    TRIAL_QUERY_PAYLOAD=$(cat "${REQUESTS_DIR}/trial_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Clinical trial RAG query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8538/v1/trial/query" "$TRIAL_QUERY_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Clinical trial RAG query complete"
    else
        fail "Clinical trial agent query unavailable"
    fi
fi

# ── Neurology: RAG Query ──────────────────────────────────────────────────
section "Neurology Intelligence Agent (:8528) — RAG Query"
step "POST http://localhost:8528/v1/neuro/query"

NEURO_QUERY_PAYLOAD='{
    "query": "Neurotoxicity from intrathecal methotrexate in pediatric ALL. CNS prophylaxis protocols and leukoencephalopathy monitoring.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/neurology_query.json" ]; then
    NEURO_QUERY_PAYLOAD=$(cat "${REQUESTS_DIR}/neurology_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Neurology RAG query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8528/v1/neuro/query" "$NEURO_QUERY_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Neurology RAG query complete"
    else
        fail "Neurology agent query unavailable"
    fi
fi

# ── Single-Cell Intelligence ───────────────────────────────────────────────
section "Single-Cell Intelligence Agent (:8540)"
step "POST http://localhost:8540/v1/sc/query"

SC_PAYLOAD='{
    "query": "Single-cell RNA sequencing of bone marrow in pediatric B-cell ALL. Leukemic blast heterogeneity, immune microenvironment, and MRD clone tracking.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/single_cell_query.json" ]; then
    SC_PAYLOAD=$(cat "${REQUESTS_DIR}/single_cell_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Single-cell query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8540/v1/sc/query" "$SC_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Single-cell query complete"
    else
        fail "Single-cell agent query unavailable"
    fi
fi

# ── Rare Disease: RAG Query ───────────────────────────────────────────────
section "Rare Disease Diagnostic Agent (:8134) — RAG Query"
step "POST http://localhost:8134/v1/diagnostic/query"

RD_QUERY_PAYLOAD='{
    "query": "Germline predisposition syndromes associated with pediatric ALL. PAX5 and IKZF1 constitutional variants and Li-Fraumeni overlap.",
    "top_k": 5
}'

if [ -f "${REQUESTS_DIR}/rare_disease_query.json" ]; then
    RD_QUERY_PAYLOAD=$(cat "${REQUESTS_DIR}/rare_disease_query.json")
fi

echo ""
if [ "$QUICK_MODE" = true ]; then
    skip "Rare disease RAG query (requires LLM)"
else
    if result=$(api_call POST "http://localhost:8134/v1/diagnostic/query" "$RD_QUERY_PAYLOAD"); then
        echo "$result" | pretty_json
        ok "Rare disease RAG query complete"
    else
        fail "Rare disease agent query unavailable"
    fi
fi

# ============================================================================
# PHASE 4: Summary
# ============================================================================

banner "Demo Summary"

echo ""
echo -e "  ${BOLD}Patient:${RESET}  PEDS-ALL-001 | 8F | B-cell ALL with ETV6-RUNX1"
echo ""
echo -e "  ${GREEN}Passed:${RESET}   ${PASSED}"
echo -e "  ${RED}Failed:${RESET}   ${FAILED}"
echo -e "  ${YELLOW}Skipped:${RESET}  ${SKIPPED}"
echo -e "  ${DIM}Total:    $((PASSED + FAILED + SKIPPED))${RESET}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "  ${GREEN}${BOLD}All demo steps completed successfully.${RESET}"
else
    echo -e "  ${YELLOW}${BOLD}Some steps failed — check service health above.${RESET}"
fi

echo ""
echo -e "  ${DIM}Demo completed at $(date '+%Y-%m-%d %H:%M:%S')${RESET}"
echo -e "  ${DIM}Run with --health for health checks only, --quick to skip LLM queries${RESET}"
echo ""
