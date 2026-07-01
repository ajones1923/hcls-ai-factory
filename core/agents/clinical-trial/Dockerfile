# =============================================================================
# Clinical Trial Intelligence Agent -- Dockerfile
# HCLS AI Factory / ai_agent_adds / clinical_trial_intelligence_agent
#
# Multi-purpose image: runs Streamlit UI (8128), FastAPI server (8538),
# or one-shot setup/seed scripts depending on CMD override.
#
# Author: Adam Jones
# Date:   March 2026
# =============================================================================

# -- Stage 1: Builder --------------------------------------------------------
FROM python:3.10-slim AS builder

WORKDIR /build

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

# -- Stage 2: Runtime --------------------------------------------------------
FROM python:3.10-slim

LABEL maintainer="Adam Jones"
LABEL description="Clinical Trial Intelligence Agent -- HCLS AI Factory"
LABEL version="1.0.0"

WORKDIR /app

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

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY config/   /app/config/
COPY src/       /app/src/
COPY app/       /app/app/
COPY api/       /app/api/
COPY scripts/   /app/scripts/
COPY data/      /app/data/
COPY .streamlit/ /app/.streamlit/

ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1

RUN useradd -r -s /bin/false trialuser \
    && mkdir -p /app/data/cache /app/data/reference /app/data/events \
    && chown -R trialuser:trialuser /app
USER trialuser

EXPOSE 8128
EXPOSE 8538

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:8538/health || exit 1

CMD ["streamlit", "run", "app/trial_ui.py", \
     "--server.port=8128", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--browser.gatherUsageStats=false"]
