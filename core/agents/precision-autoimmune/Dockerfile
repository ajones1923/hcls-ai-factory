###############################################################################
# Precision Autoimmune Intelligence Agent — Dockerfile
#
# Multi-stage build: builder → runtime
# Default CMD: Streamlit UI on port 8531
# Override for API: uvicorn api.main:app --host 0.0.0.0 --port 8532
###############################################################################

# ── Stage 1: Builder ────────────────────────────────────────────────────────
FROM python:3.10-slim AS builder

WORKDIR /build

COPY requirements.txt .
RUN python -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir --upgrade pip && \
    /opt/venv/bin/pip install --no-cache-dir -r requirements.txt

# ── Stage 2: Runtime ────────────────────────────────────────────────────────
FROM python:3.10-slim AS runtime

# System deps
RUN apt-get update && \
    apt-get install -y --no-install-recommends tini && \
    rm -rf /var/lib/apt/lists/*

# Copy virtualenv
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONUNBUFFERED=1

# Create non-root user
RUN useradd -m -s /bin/bash autouser
WORKDIR /app

# Copy application
COPY config/ config/
COPY src/ src/
COPY api/ api/
COPY app/ app/
COPY scripts/ scripts/
COPY data/ data/
COPY .streamlit/ .streamlit/
COPY requirements.txt .

# Fix permissions
RUN chown -R autouser:autouser /app

USER autouser

# Expose ports
EXPOSE 8531 8532

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8531/_stcore/health')" || exit 1

# Default: Streamlit UI
ENTRYPOINT ["tini", "--"]
CMD ["streamlit", "run", "app/autoimmune_ui.py", \
     "--server.port=8531", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
