# ═══════════════════════════════════════════════════════════════════
# MRI Brain MS Lesion Detection — MONAI Application Package Container
# ═══════════════════════════════════════════════════════════════════
#
# Build:
#   docker build -t mri_brain_ms_lesion_map:1.0.0 -f maps/mri_brain_ms_lesion/Dockerfile .
#
# Run (standalone):
#   docker run --gpus all -v /path/to/dicom:/input -v /path/to/output:/output \
#     mri_brain_ms_lesion_map:1.0.0 --input /input --output /output
#
# Package as MAP:
#   monai-deploy package maps/mri_brain_ms_lesion/app.py \
#     -o mri_brain_ms_lesion_map.map \
#     -t mri_brain_ms_lesion_map:1.0.0
#
# License: Apache 2.0

FROM nvcr.io/nvidia/pytorch:24.01-py3

LABEL maintainer="HCLS AI Factory"
LABEL version="1.0.0"
LABEL description="MRI Brain MS Lesion Detection MAP — FLAIR/T2 lesion volumetry"

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

# Copy the MRI brain MS lesion MAP
COPY maps/mri_brain_ms_lesion/ /app/maps/mri_brain_ms_lesion/

# Copy the underlying workflow and models
COPY src/workflows/base.py /app/src/workflows/base.py
COPY src/workflows/mri_brain_ms_lesion.py /app/src/workflows/mri_brain_ms_lesion.py
COPY src/workflows/__init__.py /app/src/workflows/__init__.py
COPY src/models.py /app/src/models.py
COPY src/__init__.py /app/src/__init__.py

# Create input/output directories
RUN mkdir -p /input /output

# Default environment
ENV PYTHONPATH=/app
ENV MONAI_BUNDLE_DIR=/tmp/monai_bundles

ENTRYPOINT ["python", "-m", "maps.mri_brain_ms_lesion"]
CMD ["--input", "/input", "--output", "/output", "--mock"]
