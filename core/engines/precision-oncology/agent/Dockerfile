# =============================================================================
# Oncology Intelligence Agent — Dockerfile
# HCLS AI Factory / ai_agent_adds / precision_oncology_agent
#
# Multi-purpose image: runs Streamlit UI (8526), FastAPI server (8527),
# or one-shot setup/seed scripts depending on CMD override.
#
# Author: Adam Jones
# Date:   February 2026
# =============================================================================

# Stage 1: Builder
FROM python:3.10-slim AS builder
WORKDIR /build
RUN apt-get update && apt-get install -y --no-install-recommends build-essential gcc g++ libxml2-dev libxslt1-dev && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir -r requirements.txt

# Stage 2: Runtime
FROM python:3.10-slim
LABEL maintainer="Adam Jones"
LABEL description="Oncology Intelligence Agent — HCLS AI Factory"
LABEL version="1.0.0"
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends curl libgomp1 libxml2 libxslt1.1 && rm -rf /var/lib/apt/lists/*
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
COPY config/ /app/config/
COPY src/ /app/src/
COPY api/ /app/api/
COPY app/ /app/app/
COPY scripts/ /app/scripts/
COPY data/ /app/data/
ENV PYTHONPATH="/app"
ENV PYTHONUNBUFFERED=1
RUN useradd -r -s /bin/false oncouser && mkdir -p /app/data/cache /app/data/reference && chown -R oncouser:oncouser /app
USER oncouser
EXPOSE 8526 8527
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 CMD curl -f http://localhost:8526/_stcore/health || exit 1
CMD ["streamlit", "run", "app/oncology_ui.py", "--server.port=8526", "--server.address=0.0.0.0", "--server.headless=true", "--browser.gatherUsageStats=false"]
