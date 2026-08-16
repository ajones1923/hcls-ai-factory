# TSC Intelligence Engine (Engine 7) — Apache 2.0
FROM python:3.12-slim

WORKDIR /app
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends curl && rm -rf /var/lib/apt/lists/*

COPY requirements.txt requirements-prod.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt -r requirements-prod.txt

COPY . .

# Run unprivileged. This is the flagship disease program and it served as root until 2026-08-15 --
# a container escape would have landed with full root in the namespace. `COPY . .` also ingests the
# whole build context, so .dockerignore is load-bearing here (secrets, logs, data must not enter
# the image).
RUN useradd -r -s /bin/false tscuser && chown -R tscuser:tscuser /app
USER tscuser

EXPOSE 8560 8561 8562 8563

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:8560/health || exit 1

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8560"]
