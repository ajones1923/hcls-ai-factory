# =============================================================================
# Single-Cell Analysis Engine (Engine 8 — compute layer) — Dockerfile
#
# Serves the scanpy scRNA-seq workflow behind FastAPI on :8573, the endpoint the capability
# registry already advertises (`singlecell-compute`, tags: D1 / scanpy / cpu-served).
#
# Before this file existed the capability was registered `live` with `endpoint: localhost:8573`
# but had no Dockerfile, no requirements and no compose entry, so nothing could ever bind that
# port. The code itself was always real -- a genuine scanpy pipeline verified on PBMC 3k
# (2,700 cells -> 9 clusters -> CD4 T / B / NK / CD14+ & FCGR3A+ monocytes / DC / megakaryocytes).
#
# CPU-served by design: scanpy + leidenalg run on the Grace cores. No CUDA needed.
# =============================================================================

# Stage 1: builder — compile the scientific wheels once
FROM python:3.11-slim AS builder
WORKDIR /build
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ gfortran libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Stage 2: runtime
FROM python:3.11-slim
LABEL maintainer="Adam Jones"
LABEL description="Single-Cell Analysis Engine — HCLS AI Factory"
LABEL version="1.0.0"
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl libgomp1 libopenblas0 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY src/ /app/src/
COPY data/ /app/data/

# src/ is on the path so `from single_cell_compute import ...` resolves. This engine ships NO
# module that shadows the standard library, so putting src/ on PYTHONPATH is safe here -- it is
# NOT safe for the subjects that carry src/collections.py (see scripts/run_all_tests.py).
ENV PYTHONPATH="/app:/app/src"
ENV PYTHONUNBUFFERED=1

RUN useradd -r -s /bin/false scuser && chown -R scuser:scuser /app
USER scuser

EXPOSE 8573
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8573/health || exit 1

# create_app() is a factory, so uvicorn needs --factory
CMD ["uvicorn", "single_cell_service:_app_factory", "--factory", \
     "--host", "0.0.0.0", "--port", "8573"]
