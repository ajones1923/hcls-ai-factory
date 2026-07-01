# Drug Discovery Pipeline UI
# Streamlit-based interface for VCP/FTD drug discovery demo

# ---- Stage 1: Builder ----
FROM python:3.11-slim@sha256:0b23cfb36dbe81e1cfaabab14e0a0d19e8c1d0c543a2900cbc1b8a0a12e7a8f5 AS builder

WORKDIR /build

# Install system dependencies needed for building Python packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
        curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir -r requirements.txt

# ---- Stage 2: Runtime ----
FROM python:3.11-slim@sha256:0b23cfb36dbe81e1cfaabab14e0a0d19e8c1d0c543a2900cbc1b8a0a12e7a8f5

WORKDIR /app

# Install only runtime system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy application code
COPY src/ /app/src/
COPY app/ /app/app/
COPY data/ /app/data/
COPY generate_vcp_report.py /app/
COPY generate_vcp_report_enhanced.py /app/

# Create output directories and non-root user
RUN mkdir -p /app/outputs /app/outputs/molecules \
    && useradd --create-home --shell /bin/bash appuser \
    && chown -R appuser:appuser /app
USER appuser

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Run Streamlit
ENTRYPOINT ["streamlit", "run", "app/discovery_ui.py", \
    "--server.port=8501", \
    "--server.address=0.0.0.0", \
    "--server.headless=true"]
