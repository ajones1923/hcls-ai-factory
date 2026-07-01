# RAG Chat Pipeline Portal
# Flask-based API server for ClinVar/AlphaMissense RAG pipeline

# ---- Stage 1: Builder ----
FROM python:3.11-slim AS builder

WORKDIR /build

COPY requirements.txt .

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install dependencies, skipping GPU/large packages that are pre-installed on DGX
RUN pip install --no-cache-dir \
    python-dotenv==1.2.1 \
    pydantic==2.12.5 \
    pydantic-settings==2.12.0 \
    loguru==0.7.3 \
    pymilvus==2.6.6 \
    openai==2.15.0 \
    anthropic==0.75.0 \
    streamlit==1.52.2 \
    streamlit-chat==0.1.1 \
    pandas==2.3.3 \
    numpy==2.4.0 \
    tqdm==4.67.1 \
    fastapi==0.128.0 \
    uvicorn==0.40.0 \
    flask==3.1.2 \
    flask-cors==6.0.2 \
    psutil==7.2.1 \
    pynvml==13.0.1

# ---- Stage 2: Runtime ----
FROM python:3.11-slim

WORKDIR /app

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy portal application code
COPY portal/ /app/portal/

# Create non-root user
RUN useradd -r -s /bin/false appuser && chown -R appuser:appuser /app
USER appuser

EXPOSE 5001

HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:5001/health')" || exit 1

CMD ["python", "portal/app/server.py"]
