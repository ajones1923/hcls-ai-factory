# TSC Intelligence Engine (Engine 7) — Apache 2.0
FROM python:3.12-slim

WORKDIR /app
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends curl && rm -rf /var/lib/apt/lists/*

COPY requirements.txt requirements-prod.txt ./
RUN pip install --upgrade pip && pip install -r requirements.txt -r requirements-prod.txt

COPY . .

EXPOSE 8560 8561 8562 8563
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8560"]
