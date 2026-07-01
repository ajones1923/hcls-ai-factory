# =============================================================================
# Pharmacogenomics Intelligence Agent — Dockerfile
# HCLS AI Factory / ai_agent_adds / pharmacogenomics_intelligence_agent
#
# Multi-purpose image: runs Streamlit UI (8507), FastAPI server (8107),
# or one-shot setup/seed scripts depending on CMD override.
#
# Author: Adam Jones
# Date:   March 2026
# =============================================================================

# ── Stage 1: Builder ─────────────────────────────────────────────────────────
FROM python:3.10-slim AS builder

WORKDIR /build

# System dependencies required by sentence-transformers / numpy / lxml / cyvcf2
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        libxml2-dev \
        libxslt1-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# ── Stage 2: Runtime ─────────────────────────────────────────────────────────
FROM python:3.10-slim

LABEL maintainer="Adam Jones"
LABEL description="Pharmacogenomics Intelligence Agent — HCLS AI Factory"
LABEL version="1.0.0"

WORKDIR /app

# Minimal runtime libs (libgomp needed by torch/sentence-transformers)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        libgomp1 \
        libxml2 \
        libxslt1.1 \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy application source
COPY config/   /app/config/
COPY src/       /app/src/
COPY app/       /app/app/
COPY api/       /app/api/
COPY scripts/   /app/scripts/
COPY data/      /app/data/
COPY .streamlit/ /app/.streamlit/

# Ensure Python can find the project root
ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1

# Create non-root user
RUN useradd -r -s /bin/false pgxuser \
    && mkdir -p /app/data/cache /app/data/reference /app/data/events \
    && chown -R pgxuser:pgxuser /app
USER pgxuser

# Expose Streamlit and FastAPI ports
EXPOSE 8507
EXPOSE 8107

# Healthcheck against Streamlit (default service)
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:8507/_stcore/health || exit 1

# Default: launch Streamlit UI
CMD ["streamlit", "run", "app/pgx_ui.py", \
     "--server.port=8507", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--browser.gatherUsageStats=false"]
