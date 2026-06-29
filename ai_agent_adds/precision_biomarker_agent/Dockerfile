# =============================================================================
# Biomarker Intelligence Agent — Dockerfile
# HCLS AI Factory / ai_agent_adds / precision_biomarker_agent
#
# Multi-purpose image: runs Streamlit UI (8528), FastAPI server (8529),
# or one-shot setup/seed scripts depending on CMD override.
#
# Author: Adam Jones
# Date:   March 2026
# =============================================================================

# -- Stage 1: Builder --------------------------------------------------------
FROM python:3.10-slim AS builder

WORKDIR /build

# System dependencies required by sentence-transformers / numpy
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        libxml2-dev \
        libxslt1-dev \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# -- Stage 2: Runtime --------------------------------------------------------
FROM python:3.10-slim

LABEL maintainer="Adam Jones"
LABEL description="Biomarker Intelligence Agent — HCLS AI Factory"
LABEL version="1.0.0"

WORKDIR /app

# Minimal runtime libs (libgomp needed by torch/sentence-transformers)
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
        libxml2 \
        libxslt1.1 \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy application source
COPY config/       /app/config/
COPY src/           /app/src/
COPY api/           /app/api/
COPY app/           /app/app/
COPY scripts/       /app/scripts/
COPY data/          /app/data/
COPY .streamlit/    /app/.streamlit/

# Ensure Python can find the project root
ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1

# Create non-root user
RUN useradd -r -s /bin/false biomarkeruser \
    && mkdir -p /app/data/cache /app/data/reference \
    && chown -R biomarkeruser:biomarkeruser /app
USER biomarkeruser

# Expose Streamlit and FastAPI ports
EXPOSE 8528
EXPOSE 8529

# Health check: uses Python stdlib (no curl dependency required)
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8528/health')" || exit 1

# Default: launch Streamlit UI
CMD ["streamlit", "run", "app/biomarker_ui.py", \
     "--server.port=8528", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--browser.gatherUsageStats=false"]
