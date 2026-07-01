# ═══════════════════════════════════════════════════════════════════
# CXR Rapid Findings — MONAI Application Package Container
# ═══════════════════════════════════════════════════════════════════
#
# Build:
#   docker build -t cxr_rapid_findings_map:1.0.0 -f maps/cxr_rapid_findings/Dockerfile .
#
# Run (standalone):
#   docker run --gpus all -v /path/to/dicom:/input -v /path/to/output:/output \
#     cxr_rapid_findings_map:1.0.0 --input /input --output /output
#
# Package as MAP:
#   monai-deploy package maps/cxr_rapid_findings/app.py \
#     -o cxr_rapid_findings_map.map \
#     -t cxr_rapid_findings_map:1.0.0
#
# License: Apache 2.0

FROM nvcr.io/nvidia/pytorch:24.01-py3

LABEL maintainer="HCLS AI Factory"
LABEL version="1.0.0"
LABEL description="CXR Rapid Findings MAP — DenseNet-121 multi-label classification"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    dcmtk \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    monai-deploy-app-sdk>=0.6.0 \
    monai>=1.3.0 \
    torchxrayvision>=1.2.0 \
    torchvision>=0.15.0 \
    pydicom>=2.4.0 \
    highdicom>=0.22.0 \
    SimpleITK>=2.3.0 \
    nibabel>=5.0.0 \
    pydantic>=2.0 \
    loguru>=0.7.0 \
    Pillow>=10.0.0

# Set up application directory
WORKDIR /app

# Copy common MAP infrastructure
COPY maps/common/ /app/maps/common/

# Copy the CXR rapid findings MAP
COPY maps/cxr_rapid_findings/ /app/maps/cxr_rapid_findings/

# Copy the underlying workflow and models
COPY src/workflows/base.py /app/src/workflows/base.py
COPY src/workflows/cxr_rapid_findings.py /app/src/workflows/cxr_rapid_findings.py
COPY src/workflows/__init__.py /app/src/workflows/__init__.py
COPY src/models.py /app/src/models.py
COPY src/__init__.py /app/src/__init__.py

# Create input/output directories
RUN mkdir -p /input /output

# Default environment
ENV PYTHONPATH=/app

ENTRYPOINT ["python", "-m", "maps.cxr_rapid_findings"]
CMD ["--input", "/input", "--output", "/output", "--mock"]
