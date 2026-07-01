# ═══════════════════════════════════════════════════════════════════
# Breast BI-RADS Classification — MONAI Application Package Container
# ═══════════════════════════════════════════════════════════════════
#
# Build:
#   docker build -t breast_birads_map:1.0.0 -f maps/breast_birads/Dockerfile .
#
# Run (standalone):
#   docker run --gpus all -v /path/to/dicom:/input -v /path/to/output:/output \
#     breast_birads_map:1.0.0 --input /input --output /output
#
# Package as MAP:
#   monai-deploy package maps/breast_birads/app.py \
#     -o breast_birads_map.map \
#     -t breast_birads_map:1.0.0
#
# License: Apache 2.0

FROM nvcr.io/nvidia/pytorch:24.01-py3

LABEL maintainer="HCLS AI Factory"
LABEL version="1.0.0"
LABEL description="Breast BI-RADS Classification MAP — BI-RADS 5th Edition"

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
    loguru>=0.7.0 \
    Pillow>=10.0.0

# Set up application directory
WORKDIR /app

# Copy common MAP infrastructure
COPY maps/common/ /app/maps/common/

# Copy the breast BI-RADS MAP
COPY maps/breast_birads/ /app/maps/breast_birads/

# Copy the underlying workflow and models
COPY src/workflows/base.py /app/src/workflows/base.py
COPY src/workflows/breast_birads.py /app/src/workflows/breast_birads.py
COPY src/workflows/__init__.py /app/src/workflows/__init__.py
COPY src/models.py /app/src/models.py
COPY src/__init__.py /app/src/__init__.py

# Create input/output directories
RUN mkdir -p /input /output

# Default environment
ENV PYTHONPATH=/app

ENTRYPOINT ["python", "-m", "maps.breast_birads"]
CMD ["--input", "/input", "--output", "/output", "--mock"]
