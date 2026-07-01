# ═══════════════════════════════════════════════════════════════════
# CT Head Hemorrhage Triage — MONAI Application Package Container
# ═══════════════════════════════════════════════════════════════════
#
# Build:
#   docker build -t ct_head_hemorrhage_map:1.0.0 -f maps/ct_head_hemorrhage/Dockerfile .
#
# Run (standalone):
#   docker run --gpus all -v /path/to/dicom:/input -v /path/to/output:/output \
#     ct_head_hemorrhage_map:1.0.0 --input /input --output /output
#
# Package as MAP:
#   monai-deploy package maps/ct_head_hemorrhage/app.py \
#     -o ct_head_hemorrhage_map.map \
#     -t ct_head_hemorrhage_map:1.0.0
#
# License: Apache 2.0

FROM nvcr.io/nvidia/pytorch:24.01-py3

LABEL maintainer="HCLS AI Factory"
LABEL version="1.0.0"
LABEL description="CT Head Hemorrhage Triage MAP — SegResNet + Brain Trauma Foundation"

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
    pydicom>=2.4.0 \
    highdicom>=0.22.0 \
    SimpleITK>=2.3.0 \
    nibabel>=5.0.0 \
    pydantic>=2.0 \
    loguru>=0.7.0

# Set up application directory
WORKDIR /app

# Copy common MAP infrastructure
COPY maps/common/ /app/maps/common/

# Copy the CT head hemorrhage MAP
COPY maps/ct_head_hemorrhage/ /app/maps/ct_head_hemorrhage/

# Copy the underlying workflow and models
COPY src/workflows/base.py /app/src/workflows/base.py
COPY src/workflows/ct_head_hemorrhage.py /app/src/workflows/ct_head_hemorrhage.py
COPY src/workflows/__init__.py /app/src/workflows/__init__.py
COPY src/models.py /app/src/models.py
COPY src/__init__.py /app/src/__init__.py

# Create input/output directories
RUN mkdir -p /input /output

# Default environment
ENV PYTHONPATH=/app
ENV MONAI_BUNDLE_DIR=/tmp/monai_bundles

ENTRYPOINT ["python", "-m", "maps.ct_head_hemorrhage"]
CMD ["--input", "/input", "--output", "/output", "--mock"]
