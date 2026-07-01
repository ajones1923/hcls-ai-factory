# =============================================================================
# Imaging Intelligence Agent — Dockerfile
# HCLS AI Factory / ai_agent_adds / imaging_intelligence_agent
#
# Multi-purpose image: runs Streamlit UI (8525), FastAPI server (8524),
# or one-shot setup/seed scripts depending on CMD override.
#
# Author: Adam Jones
# Date:   February 2026
# =============================================================================

# ── Stage 1: Builder ─────────────────────────────────────────────────────────
FROM python:3.10-slim AS builder

WORKDIR /build

# System dependencies required by sentence-transformers / numpy / lxml / SimpleITK
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

# ── Stage 2: Runtime ─────────────────────────────────────────────────────────
FROM python:3.10-slim

LABEL maintainer="Adam Jones"
LABEL description="Imaging Intelligence Agent — HCLS AI Factory"
LABEL version="1.0.0"

WORKDIR /app

# Minimal runtime libs (libgomp needed by torch/sentence-transformers)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        libgomp1 \
        libxml2 \
        libxslt1.1 \
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

# Ensure Python can find the project root
ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1

# Create non-root user with home directory for HuggingFace cache
RUN useradd -r -m -s /bin/false imaginguser \
    && mkdir -p /app/data/cache /app/data/reference /app/data/dicom \
    && mkdir -p /home/imaginguser/.cache/huggingface \
    && chown -R imaginguser:imaginguser /app /home/imaginguser
USER imaginguser

# Expose Streamlit and FastAPI ports
EXPOSE 8524
EXPOSE 8525

# Healthcheck against Streamlit (default service)
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:8525/_stcore/health || exit 1

# Default: launch Streamlit UI
CMD ["streamlit", "run", "app/imaging_ui.py", \
     "--server.port=8525", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--browser.gatherUsageStats=false"]
