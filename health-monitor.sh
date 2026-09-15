#!/bin/bash
# ============================================================================
# HCLS AI Factory — Health Monitor & Auto-Recovery
#
# Monitors all pipeline services, detects failures, auto-restarts, and logs.
#
# Usage:
#   ./health-monitor.sh status          # One-shot health check (human-readable)
#   ./health-monitor.sh status --json   # One-shot health check (JSON output)
#   ./health-monitor.sh fix             # Check + auto-restart any failed services
#   ./health-monitor.sh watch           # Continuous watchdog (runs every 60s)
#   ./health-monitor.sh restart <svc>   # Restart a specific service
#   ./health-monitor.sh restart all     # Restart all services
#   ./health-monitor.sh stop <svc>      # Stop a specific service
#   ./health-monitor.sh log             # Show recent health monitor log
#   ./health-monitor.sh install         # Install as cron job (every 5 min)
#   ./health-monitor.sh uninstall       # Remove cron job
# ============================================================================

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="${SCRIPT_DIR}/logs/health-monitor.log"
LOG_DIR="${SCRIPT_DIR}/logs"
WATCH_INTERVAL=${WATCH_INTERVAL:-60}
MAX_LOG_SIZE=10485760  # 10MB

mkdir -p "${LOG_DIR}"

# ============================================================================
# SINGLE-INSTANCE GUARD  +  SELF-TIMEOUT
# ============================================================================
# A run that has to give up on N dead services costs N x 60s (30 attempts x 2s).
# On a partly-built box that exceeds the 5-minute cron cadence, so invocations
# pile up: 11 concurrent monitors -- three of them 51 days old -- were observed
# on 2026-09-15, racing each other to start the same services and appending
# 536 MB to logs/. flock makes a late run skip instead of stack.
#
# The self-timeout is the other half: without it a single hung run would hold
# the lock forever and silently end all supervision.
#
# Read-only commands (status / log / help) are deliberately NOT guarded.
LOCK_FILE="${LOG_DIR}/health-monitor.lock"
MAX_RUN_SECONDS=${MAX_RUN_SECONDS:-600}

acquire_lock() {
    exec 9>"${LOCK_FILE}" || return 0
    if ! flock -n 9; then
        echo "Another health-monitor run still holds ${LOCK_FILE} — skipping this tick."
        log "INFO" "Skipped tick: previous run still in progress"
        exit 0
    fi
    # Watchdog: SIGTERM this run if it outlives MAX_RUN_SECONDS, so the lock is
    # always released. `exec 9>&-` is load-bearing -- without it the watchdog
    # inherits the lock fd and keeps the lock held after the parent is killed,
    # silently blocking every later tick for the full MAX_RUN_SECONDS.
    ( exec 9>&-; sleep "${MAX_RUN_SECONDS}"; kill -TERM $$ 2>/dev/null ) &
    WATCHDOG_PID=$!
    trap 'kill "${WATCHDOG_PID}" 2>/dev/null' EXIT
}

# Cap the cron append-log, which nothing else rotates.
cap_cron_log() {
    local f="${LOG_DIR}/cron-health.log"
    [ -f "$f" ] || return 0
    local size
    size=$(stat --format="%s" "$f" 2>/dev/null || echo 0)
    if [ "$size" -gt "$MAX_LOG_SIZE" ]; then
        tail -c "$((MAX_LOG_SIZE / 2))" "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
        log "INFO" "cron-health.log truncated at ${size} bytes"
    fi
}

# ============================================================================
# SERVICE REGISTRY
# ============================================================================
# Format: ID|PORT|NAME|HEALTH_PATH|TYPE|START_DIR|START_CMD
#
# TYPE: python = Python process, streamlit = Streamlit app, docker = Docker container, external = no restart
#
declare -a SERVICES=(
    "genomics|5000|Genomics Portal|/health|python|${SCRIPT_DIR}/core/engines/genomic-foundation/web-portal|./venv/bin/python app/server.py"
    "rag-api|5001|RAG/Chat API|/health|python|${SCRIPT_DIR}/core/engines/precision-intelligence|./venv/bin/python3 portal/app/server.py"
    "rag-chat|8501|RAG Chat UI|/healthz|streamlit|${SCRIPT_DIR}/core/engines/precision-intelligence|./venv/bin/streamlit run app/chat_ui.py --server.port 8501 --server.address 0.0.0.0 --server.headless true"
    "drug-discovery|8505|Drug Discovery|/healthz|streamlit|${SCRIPT_DIR}/core/engines/therapeutic-discovery|./venv/bin/streamlit run app/discovery_ui.py --server.port 8505 --server.headless true"
    "drug-portal|8510|Discovery Portal|/healthz|streamlit|${SCRIPT_DIR}/hcls-orchestrator|source ${SCRIPT_DIR}/core/engines/therapeutic-discovery/venv/bin/activate && streamlit run portal/app.py --server.port 8510 --server.headless true"
    "landing|8080|Landing Page|/health|python|${SCRIPT_DIR}/landing-page|./venv/bin/python server.py"
    "milvus|19530|Milvus|TCP|docker|${SCRIPT_DIR}/core/engines/precision-intelligence|docker compose up -d milvus"
    "grafana|3000|Grafana|/api/health|docker|${SCRIPT_DIR}/core/engines/therapeutic-discovery/monitoring|docker compose up -d grafana"
    "prometheus|9099|Prometheus|/-/healthy|docker|${SCRIPT_DIR}/core/engines/therapeutic-discovery/monitoring|docker compose up -d prometheus"
    "node-exporter|9100|Node Exporter|/metrics|docker|${SCRIPT_DIR}/core/engines/therapeutic-discovery/monitoring|docker compose up -d node-exporter"
    "dcgm|9400|DCGM Exporter|/metrics|docker|${SCRIPT_DIR}/core/engines/therapeutic-discovery/monitoring|docker compose up -d dcgm-exporter"
    # ── Intelligence Agents ──
    "biomarker|8529|Precision Biomarker|/health|python|${SCRIPT_DIR}/core/agents/precision-biomarker|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8529"
    "biomarker-ui|8528|Biomarker UI|/healthz|streamlit|${SCRIPT_DIR}/core/agents/precision-biomarker|./venv/bin/streamlit run app/biomarker_ui.py --server.port 8528 --server.address 0.0.0.0 --server.headless true"
    "oncology|8527|Precision Oncology|/health|python|${SCRIPT_DIR}/core/engines/precision-oncology/agent|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8527"
    "cart|8522|CAR-T Intelligence|/health|python|${SCRIPT_DIR}/core/agents/cart|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8522"
    "imaging|8524|Imaging Intelligence|/health|python|${SCRIPT_DIR}/core/engines/clinical-imaging/agent|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8524"
    "autoimmune|8532|Precision Autoimmune|/health|python|${SCRIPT_DIR}/core/agents/precision-autoimmune|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8532"
    "pharmacogenomics|8508|Pharmacogenomics|/health|python|${SCRIPT_DIR}/core/agents/pharmacogenomics|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8508"
    "cardiology|8127|Cardiology Intelligence|/health|python|${SCRIPT_DIR}/core/engines/cardiology|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8127"
    "clinical-trial|8539|Clinical Trial Intelligence|/health|python|${SCRIPT_DIR}/core/agents/clinical-trial|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8539"
    "rare-disease|8545|Rare Disease Diagnostic|/health|python|${SCRIPT_DIR}/core/agents/rare-disease-diagnostic|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8545"
    "neurology|8536|Neurology Intelligence|/health|python|${SCRIPT_DIR}/core/agents/neurology|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8536"
    "neurology-ui|8535|Neurology UI|/healthz|streamlit|${SCRIPT_DIR}/core/agents/neurology|./venv/bin/streamlit run app/neuro_ui.py --server.port 8535 --server.address 0.0.0.0 --server.headless true"
    "single-cell|8541|Single-Cell Intelligence|/health|python|${SCRIPT_DIR}/core/agents/single-cell|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8541"
    # ── Model / compute services and the disease program ──
    # These were registered `live` while nothing ran them: they expose a create_app() factory
    # rather than a module-level `app`, so they need uvicorn --factory, and nothing supervised
    # them. Added so they survive a reboot (PRD R10/A3) and so --probe stays green.
    "variant-store|8575|Variant Store|/stats|python|${SCRIPT_DIR}/core/engines/genomic-foundation/src|${SCRIPT_DIR}/.venv/bin/python -m uvicorn --factory variant_store_service:create_app --host 0.0.0.0 --port 8575"
    "singlecell-compute|8573|Single-Cell Compute|/docs|python|${SCRIPT_DIR}/core/engines/single-cell/src|${SCRIPT_DIR}/.venv/bin/python -m uvicorn --factory single_cell_service:create_app --host 0.0.0.0 --port 8573"
    "proteinmpnn|8578|ProteinMPNN Design|/docs|python|${SCRIPT_DIR}/core/engines/structural-biology/src|${SCRIPT_DIR}/.venv/bin/python -m uvicorn --factory proteinmpnn_service:create_app --host 0.0.0.0 --port 8578"
    "molecule-generator|8574|Molecule Generator|/docs|python|${SCRIPT_DIR}/core/engines/therapeutic-discovery/small-molecule/src|${SCRIPT_DIR}/.venv/bin/python -m uvicorn --factory molecule_gen_service:create_app --host 0.0.0.0 --port 8574"
    "tuberous-sclerosis|8560|TSC Intelligence Engine|/health|python|${SCRIPT_DIR}/core/disease-programs/tuberous-sclerosis|./venv/bin/python -m uvicorn api.main:app --host 0.0.0.0 --port 8560"
)

# ============================================================================
# LOGGING
# ============================================================================
log() {
    local level="$1"
    shift
    local msg="$*"
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "${timestamp} [${level}] ${msg}" >> "${LOG_FILE}"

    # Rotate log if too large
    if [ -f "${LOG_FILE}" ]; then
        local size
        size=$(stat --format="%s" "${LOG_FILE}" 2>/dev/null || echo "0")
        if [ "$size" -gt "$MAX_LOG_SIZE" ]; then
            mv "${LOG_FILE}" "${LOG_FILE}.1"
            echo "${timestamp} [INFO] Log rotated" > "${LOG_FILE}"
        fi
    fi
}

# ============================================================================
# SERVICE FIELD EXTRACTORS
# ============================================================================
get_field() {
    local service="$1"
    local field="$2"
    echo "$service" | cut -d'|' -f"$field"
}

get_id()       { get_field "$1" 1; }
get_port()     { get_field "$1" 2; }
get_name()     { get_field "$1" 3; }
get_path()     { get_field "$1" 4; }
get_type()     { get_field "$1" 5; }
get_dir()      { get_field "$1" 6; }
get_cmd()      { get_field "$1" 7; }

# Find a service entry by ID
find_service() {
    local target_id="$1"
    for svc in "${SERVICES[@]}"; do
        if [ "$(get_id "$svc")" = "$target_id" ]; then
            echo "$svc"
            return 0
        fi
    done
    return 1
}

# ============================================================================
# HEALTH CHECKS
# ============================================================================
check_port_open() {
    local port="$1"
    local timeout="${2:-2}"
    timeout "$timeout" bash -c "echo > /dev/tcp/localhost/$port" 2>/dev/null
    return $?
}

check_http_health() {
    local port="$1"
    local path="$2"
    local timeout="${3:-3}"
    local code
    code=$(curl -s -o /dev/null -w "%{http_code}" --connect-timeout "$timeout" --max-time "$timeout" "http://localhost:${port}${path}" 2>/dev/null)
    if [ -n "$code" ] && [ "$code" -ge 200 ] && [ "$code" -lt 500 ]; then
        return 0
    fi
    return 1
}

check_service() {
    local svc="$1"
    local port
    port=$(get_port "$svc")
    local path
    path=$(get_path "$svc")

    if [ "$path" = "TCP" ]; then
        check_port_open "$port" 2
    else
        check_http_health "$port" "$path" 3
    fi
    return $?
}

# Check GPU health
check_gpu() {
    if nvidia-smi > /dev/null 2>&1; then
        return 0
    fi
    return 1
}

# GB10 is a UNIFIED-memory part: the GPU allocates from the same pool as the host, and
# nvidia-smi reports [N/A] for memory here, so MemFree is the only usable signal. It is
# MemFree and deliberately not MemAvailable -- the page cache counts as "available" to the
# host but is NOT reclaimable for a CUDA allocation, which is the whole reason the workbook
# prescribes dropping caches before a GPU run (PRD R5).
#
# Without this, `check_gpu` reported HEALTHY whenever nvidia-smi merely ran -- liveness, not
# usable capacity -- while only 5 GiB of 119 GiB was actually allocatable and anything the
# size of ESMFold or Parabricks would fail.
GPU_FREE_WARN_GIB=${GPU_FREE_WARN_GIB:-16}

gpu_free_gib() {
    awk '/^MemFree:/ {printf "%.1f", $2/1048576}' /proc/meminfo 2>/dev/null
}

gpu_line() {
    if ! check_gpu; then
        echo "  GPU                  :--     NOT RESPONDING"
        return
    fi
    local free
    free=$(gpu_free_gib)
    if [ -z "$free" ]; then
        echo "  GPU                  :--     HEALTHY"
    elif awk "BEGIN{exit !($free < $GPU_FREE_WARN_GIB)}"; then
        echo "  GPU                  :--     HEALTHY  (only ${free} GiB allocatable — page cache holds the rest; drop caches before a large model)"
    else
        echo "  GPU                  :--     HEALTHY  (${free} GiB allocatable)"
    fi
}

# ============================================================================
# SERVICE MANAGEMENT
# ============================================================================
get_pid_on_port() {
    local port="$1"
    ss -tlnp 2>/dev/null | grep ":${port} " | grep -oP 'pid=\K[0-9]+' | head -1
}

stop_service() {
    local svc="$1"
    local port
    port=$(get_port "$svc")
    local name
    name=$(get_name "$svc")
    local svc_type
    svc_type=$(get_type "$svc")

    if [ "$svc_type" = "docker" ]; then
        echo "  Stopping Docker service ${name}..."
        local svc_dir
        svc_dir=$(get_dir "$svc")
        local svc_id
        svc_id=$(get_id "$svc")
        # Docker compose service names map from our IDs
        cd "$svc_dir" 2>/dev/null && docker compose stop 2>/dev/null
        log "INFO" "Stopped Docker service: ${name}"
        return
    fi

    local pid
    pid=$(get_pid_on_port "$port")
    if [ -n "$pid" ]; then
        echo "  Stopping ${name} (PID ${pid})..."
        kill "$pid" 2>/dev/null
        # Wait up to 10 seconds for graceful shutdown
        local count=0
        while [ $count -lt 10 ] && kill -0 "$pid" 2>/dev/null; do
            sleep 1
            count=$((count + 1))
        done
        # Force kill if still running
        if kill -0 "$pid" 2>/dev/null; then
            kill -9 "$pid" 2>/dev/null
        fi
        log "INFO" "Stopped ${name} (PID ${pid})"
    fi
}

start_service() {
    local svc="$1"
    local port
    port=$(get_port "$svc")
    local name
    name=$(get_name "$svc")
    local svc_type
    svc_type=$(get_type "$svc")
    local svc_dir
    svc_dir=$(get_dir "$svc")
    local svc_cmd
    svc_cmd=$(get_cmd "$svc")

    # Check if already running
    if check_service "$svc"; then
        echo "  ${name} already healthy on port ${port}"
        return 0
    fi

    echo "  Starting ${name} on port ${port}..."
    local log_name
    log_name=$(get_id "$svc")

    # Preflight 1 — the working directory must exist. Previously `cd ... 2>/dev/null &&`
    # swallowed this, started nothing, then waited the full 60s to report a generic failure.
    if [ ! -d "$svc_dir" ]; then
        echo "  SKIP: ${name} — directory not found: ${svc_dir}"
        log "WARN" "Cannot start ${name}: missing directory ${svc_dir}"
        return 1
    fi

    # Preflight 2 — most services start via a per-service ./venv. 13 of 21 of those venvs
    # do not exist on this box, and waiting 60s to discover that (per service, every 5
    # minutes) is what produced the pileup. Detect it in O(1) and say so plainly.
    # Match an absolute venv path first (some commands `source` another service's venv);
    # only fall back to the relative form, which is resolved against $svc_dir.
    local venv_path
    venv_path=$(printf '%s\n' "$svc_cmd" | grep -oE '(^|[[:space:]])/[^[:space:]]*/[.]?venv/bin/[a-zA-Z0-9_.-]+' | head -1 | tr -d '[:space:]')
    if [ -z "$venv_path" ]; then
        local venv_rel
        venv_rel=$(printf '%s\n' "$svc_cmd" | grep -oE '(^|[[:space:]])\./[.]?venv/bin/[a-zA-Z0-9_.-]+' | head -1 | tr -d '[:space:]')
        [ -n "$venv_rel" ] && venv_path="${svc_dir}/${venv_rel#./}"
    fi
    # -e, not -x: `activate` is sourced and is not executable. The question this
    # preflight answers is "has the venv been built", for which existence is the test.
    if [ -n "$venv_path" ] && [ ! -e "$venv_path" ]; then
        echo "  SKIP: ${name} — no interpreter at ${venv_path} (venv not built)"
        log "WARN" "Cannot start ${name}: ${venv_path} missing — build the service venv first"
        return 1
    fi

    # `exec 9>&-` is load-bearing in BOTH branches. fd 9 is the flock held by this run; a
    # started service inherits every open descriptor, so without closing it the service keeps
    # the lock for its entire lifetime -- and since services are meant to run forever, the
    # first successful start would permanently wedge the lock and silently end all
    # supervision. That is not hypothetical: it happened here, and every tick for the next
    # five hours logged "Skipped tick: previous run still in progress" while a Streamlit
    # process sat holding the lock.
    local child_pid=""
    if [ "$svc_type" = "docker" ]; then
        ( exec 9>&-; cd "$svc_dir" && eval "$svc_cmd" ) >> "${LOG_DIR}/${log_name}.log" 2>&1
    else
        ( exec 9>&-; cd "$svc_dir" && nohup bash -c "$svc_cmd" >> "${LOG_DIR}/${log_name}.log" 2>&1 ) &
        child_pid=$!
    fi

    # Wait for the service to come up -- but stop early if the process we launched
    # has already died, instead of burning the full 60s on a corpse.
    local attempts=0
    local max_attempts=30
    while [ $attempts -lt $max_attempts ]; do
        sleep 2
        if check_service "$svc"; then
            echo "  ${name} started successfully (port ${port})"
            log "INFO" "Started ${name} on port ${port}"
            return 0
        fi
        if [ -n "$child_pid" ] && ! kill -0 "$child_pid" 2>/dev/null; then
            echo "  FAILED: ${name} — process exited immediately; see ${LOG_DIR}/${log_name}.log"
            log "WARN" "Start of ${name} exited immediately (port ${port})"
            return 1
        fi
        attempts=$((attempts + 1))
    done

    echo "  WARNING: ${name} did not become healthy after ${max_attempts} attempts"
    log "WARN" "Failed to start ${name} on port ${port}"
    return 1
}

restart_service() {
    local svc="$1"
    local name
    name=$(get_name "$svc")
    echo "  Restarting ${name}..."
    stop_service "$svc"
    sleep 2
    start_service "$svc"
}

# ============================================================================
# COMMANDS
# ============================================================================

cmd_status() {
    local json_mode=false
    [ "${1}" = "--json" ] && json_mode=true

    local total=0
    local healthy=0
    local unhealthy_list=""

    if [ "$json_mode" = true ]; then
        echo "{"
        echo '  "timestamp": "'$(date -Iseconds)'",'
        echo '  "services": {'
    else
        echo ""
        echo "  HCLS AI Factory — Service Health"
        echo "  $(date '+%Y-%m-%d %H:%M:%S')"
        echo "  ────────────────────────────────────────────────────────"
    fi

    local first=true
    for svc in "${SERVICES[@]}"; do
        local id
        id=$(get_id "$svc")
        local port
        port=$(get_port "$svc")
        local name
        name=$(get_name "$svc")
        local path
        path=$(get_path "$svc")

        total=$((total + 1))

        if check_service "$svc"; then
            healthy=$((healthy + 1))
            if [ "$json_mode" = true ]; then
                [ "$first" = true ] || echo ","
                printf '    "%s": {"port": %s, "name": "%s", "status": "healthy"}' "$id" "$port" "$name"
                first=false
            else
                printf "  %-22s :%-6s %s\n" "$name" "$port" "HEALTHY"
            fi
        else
            unhealthy_list="${unhealthy_list} ${id}"
            if [ "$json_mode" = true ]; then
                [ "$first" = true ] || echo ","
                printf '    "%s": {"port": %s, "name": "%s", "status": "unhealthy"}' "$id" "$port" "$name"
                first=false
            else
                printf "  %-22s :%-6s %s\n" "$name" "$port" "DOWN"
            fi
        fi
    done

    # GPU check
    local gpu_status="healthy"
    if ! check_gpu; then
        gpu_status="unhealthy"
    fi

    if [ "$json_mode" = true ]; then
        echo ""
        echo "  },"
        echo "  \"gpu\": \"${gpu_status}\","
        echo "  \"gpu_free_gib\": $(gpu_free_gib || echo null),"
        echo "  \"summary\": {\"healthy\": ${healthy}, \"total\": ${total}, \"all_healthy\": $([ $healthy -eq $total ] && echo true || echo false)}"
        echo "}"
    else
        echo "  ────────────────────────────────────────────────────────"
        gpu_line
        echo "  ────────────────────────────────────────────────────────"
        echo "  Total: ${healthy}/${total} services healthy"
        if [ $healthy -eq $total ]; then
            echo "  All systems operational."
        else
            echo "  DEGRADED:${unhealthy_list}"
        fi
        echo ""
    fi

    log "INFO" "Health check: ${healthy}/${total} healthy${unhealthy_list:+ | DOWN:${unhealthy_list}}"
    [ $healthy -eq $total ]
    return $?
}

cmd_fix() {
    echo ""
    echo "  HCLS AI Factory — Auto-Recovery"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  ────────────────────────────────────────────────────────"

    local fixed=0
    local failed=0

    for svc in "${SERVICES[@]}"; do
        local id
        id=$(get_id "$svc")
        local name
        name=$(get_name "$svc")
        local port
        port=$(get_port "$svc")

        if ! check_service "$svc"; then
            echo ""
            echo "  [DOWN] ${name} (port ${port}) — attempting recovery..."
            log "WARN" "Service down: ${name} (port ${port}) — attempting restart"

            # Kill any zombie processes on the port
            local pid
            pid=$(get_pid_on_port "$port")
            if [ -n "$pid" ]; then
                echo "  Killing zombie process (PID ${pid})..."
                kill "$pid" 2>/dev/null
                sleep 2
            fi

            if start_service "$svc"; then
                fixed=$((fixed + 1))
                log "INFO" "Auto-recovered: ${name} (port ${port})"
            else
                failed=$((failed + 1))
                log "ERROR" "Failed to recover: ${name} (port ${port})"
            fi
        fi
    done

    echo ""
    echo "  ────────────────────────────────────────────────────────"

    if [ $fixed -eq 0 ] && [ $failed -eq 0 ]; then
        echo "  All services already healthy — nothing to fix."
    else
        echo "  Recovered: ${fixed}  |  Failed: ${failed}"
    fi
    echo ""

    # Run a final status check
    cmd_status
}

cmd_watch() {
    echo "HCLS AI Factory — Health Watchdog"
    echo "Checking every ${WATCH_INTERVAL}s (Ctrl+C to stop)"
    echo ""
    log "INFO" "Watchdog started (interval: ${WATCH_INTERVAL}s)"

    while true; do
        local unhealthy=false

        for svc in "${SERVICES[@]}"; do
            if ! check_service "$svc"; then
                unhealthy=true
                local name
                name=$(get_name "$svc")
                local port
                port=$(get_port "$svc")
                local id
                id=$(get_id "$svc")

                log "WARN" "Watchdog detected down service: ${name} (port ${port})"
                echo "[$(date '+%H:%M:%S')] DOWN: ${name} (port ${port}) — restarting..."

                # Kill zombie
                local pid
                pid=$(get_pid_on_port "$port")
                if [ -n "$pid" ]; then
                    kill "$pid" 2>/dev/null
                    sleep 2
                fi

                if start_service "$svc"; then
                    log "INFO" "Watchdog recovered: ${name} (port ${port})"
                    echo "[$(date '+%H:%M:%S')] RECOVERED: ${name}"
                else
                    log "ERROR" "Watchdog failed to recover: ${name} (port ${port})"
                    echo "[$(date '+%H:%M:%S')] FAILED: ${name} — manual intervention needed"
                fi
            fi
        done

        if [ "$unhealthy" = false ]; then
            echo "[$(date '+%H:%M:%S')] All ${#SERVICES[@]} services healthy"
        fi

        sleep "$WATCH_INTERVAL"
    done
}

cmd_restart() {
    local target="$1"

    if [ -z "$target" ]; then
        echo "Usage: $0 restart <service-id|all>"
        echo ""
        echo "Available services:"
        for svc in "${SERVICES[@]}"; do
            printf "  %-20s (port %s)\n" "$(get_id "$svc")" "$(get_port "$svc")"
        done
        return 1
    fi

    if [ "$target" = "all" ]; then
        echo ""
        echo "  Restarting all services..."
        echo "  ────────────────────────────────────────────────────────"
        for svc in "${SERVICES[@]}"; do
            restart_service "$svc"
        done
        echo ""
        cmd_status
        return
    fi

    local svc
    svc=$(find_service "$target")
    if [ -z "$svc" ]; then
        echo "Unknown service: ${target}"
        echo "Run '$0 restart' to see available services."
        return 1
    fi

    restart_service "$svc"
}

cmd_stop() {
    local target="$1"

    if [ -z "$target" ]; then
        echo "Usage: $0 stop <service-id>"
        return 1
    fi

    local svc
    svc=$(find_service "$target")
    if [ -z "$svc" ]; then
        echo "Unknown service: ${target}"
        return 1
    fi

    stop_service "$svc"
}

cmd_log() {
    if [ -f "${LOG_FILE}" ]; then
        tail -50 "${LOG_FILE}"
    else
        echo "No log file found at ${LOG_FILE}"
    fi
}

cmd_install() {
    local script_path
    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

    # Check if already installed
    if crontab -l 2>/dev/null | grep -q "health-monitor.sh"; then
        echo "Cron job already installed. Use 'uninstall' first to replace."
        crontab -l 2>/dev/null | grep "health-monitor.sh"
        return 0
    fi

    # Add cron job: every 5 minutes, run fix (auto-restart failed services)
    (crontab -l 2>/dev/null; echo "*/5 * * * * ${script_path} fix >> ${LOG_DIR}/cron-health.log 2>&1") | crontab -
    echo "Installed cron job: every 5 minutes"
    echo "  Command: ${script_path} fix"
    echo "  Log: ${LOG_DIR}/cron-health.log"
    log "INFO" "Cron job installed"
}

cmd_uninstall() {
    crontab -l 2>/dev/null | grep -v "health-monitor.sh" | crontab -
    echo "Cron job removed."
    log "INFO" "Cron job removed"
}

# ============================================================================
# MAIN
# ============================================================================
case "${1:-status}" in
    status)
        cmd_status "$2"
        ;;
    fix)
        acquire_lock; cap_cron_log
        cmd_fix
        ;;
    watch)
        acquire_lock
        cmd_watch
        ;;
    restart)
        acquire_lock
        cmd_restart "$2"
        ;;
    stop)
        cmd_stop "$2"
        ;;
    log|logs)
        cmd_log
        ;;
    install)
        cmd_install
        ;;
    uninstall)
        cmd_uninstall
        ;;
    help|--help|-h)
        echo "HCLS AI Factory — Health Monitor"
        echo ""
        echo "Usage: $0 <command>"
        echo ""
        echo "Commands:"
        echo "  status [--json]     Show health of all services"
        echo "  fix                 Auto-restart any failed services"
        echo "  watch               Continuous watchdog (every ${WATCH_INTERVAL}s)"
        echo "  restart <svc|all>   Restart a specific service or all"
        echo "  stop <svc>          Stop a specific service"
        echo "  log                 Show recent health monitor log"
        echo "  install             Install as cron job (every 5 min)"
        echo "  uninstall           Remove cron job"
        echo ""
        echo "Services:"
        for svc in "${SERVICES[@]}"; do
            printf "  %-20s port %-6s %s\n" "$(get_id "$svc")" "$(get_port "$svc")" "$(get_name "$svc")"
        done
        echo ""
        echo "Environment:"
        echo "  WATCH_INTERVAL      Seconds between watchdog checks (default: 60)"
        echo ""
        ;;
    *)
        echo "Unknown command: $1"
        echo "Run '$0 help' for usage."
        exit 1
        ;;
esac
