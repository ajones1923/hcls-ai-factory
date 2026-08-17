# =============================================================================
# Structural Biology Engine — Dockerfile (multi-purpose image)
#
# This engine is NOT one service. It is five FastAPI services that share a codebase:
#
#   esmfold_service         :8570   ESMFold structure prediction        (GPU)
#   protein_search_service  :8571   ESM-2 embeddings + sequence search  (CPU-feasible)
#   developability_service  :8576   antibody developability             (planned)
#   immunogenicity_service  :8577   MHCflurry immunogenicity            (planned)
#   proteinmpnn_service     :8578   ProteinMPNN inverse folding         (GPU)
#
# One image, CMD selects the service — the same pattern precision-oncology uses. Compose declares
# one entry per service.
#
# NOTE ON :8581 — the registry lists `structural-biology-engine` at localhost:8581, but NO process
# in this tree binds 8581. There is no aggregator. Until one exists the engine capability is marked
# `planned`; the model-level capabilities (8570/8571/8578) are the real, reachable surface.
#
# NOTE ON GPU — ESMFold and ProteinMPNN need CUDA. This image installs the CPU wheel, so only
# protein_search is practically servable from it. (The box itself does have a working CUDA build —
# torch 2.12.1+cu130, GB10 at sm_121 — in core/engines/therapeutic-discovery/venv.) Do not mark
# the GPU services `live` until PRD R6 lands.
#
# `vendor_rfdiffusion/` and `vendor_proteinmpnn/` are deliberately NOT copied: vendored third-party
# trees needing gated GPU packages (dgl, rfdiffusion). They are excluded from tests for the same
# reason (scripts/run_all_tests.py).
# =============================================================================

FROM python:3.11-slim AS builder
WORKDIR /build
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ gfortran libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir -r requirements.txt

FROM python:3.11-slim
LABEL maintainer="Adam Jones"
LABEL description="Structural Biology Engine — HCLS AI Factory"
LABEL version="1.0.0"
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl libgomp1 libopenblas0 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY src/ /app/src/

# src/ on the path is safe here — this engine ships no module that shadows the standard library.
# It is NOT safe for the subjects carrying src/collections.py; see scripts/run_all_tests.py.
ENV PYTHONPATH="/app:/app/src"
ENV PYTHONUNBUFFERED=1

RUN useradd -r -s /bin/false sbuser && chown -R sbuser:sbuser /app
USER sbuser

EXPOSE 8570 8571 8576 8577 8578

# Probes the DEFAULT service's port. Compose entries that override CMD to a different service must
# override this too, or the healthcheck probes a port that image instance never binds.
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8571/health || exit 1

# Default: the CPU-feasible service. Override CMD in compose for the others.
CMD ["uvicorn", "protein_search_service:_app_factory", "--factory", \
     "--host", "0.0.0.0", "--port", "8571"]
