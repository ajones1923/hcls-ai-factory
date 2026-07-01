"""Tests for deployment tier configuration.

Validates that DeploymentTier enum, tier-aware properties, and
environment variable overrides work correctly across community,
enterprise, and research tiers.

Author: Adam Jones
Date: April 2026
"""

import os
from unittest.mock import patch

import pytest

from config.settings import DeploymentTier, ImagingSettings


# ===================================================================
# TIER DEFAULTS
# ===================================================================


def test_default_is_community():
    """DeploymentTier defaults to 'community'."""
    s = ImagingSettings(ANTHROPIC_API_KEY=None)
    assert s.DEPLOYMENT_TIER == "community"


def test_deployment_tier_enum_values():
    """DeploymentTier enum exposes community, enterprise, research."""
    assert DeploymentTier.COMMUNITY.value == "community"
    assert DeploymentTier.ENTERPRISE.value == "enterprise"
    assert DeploymentTier.RESEARCH.value == "research"


# ===================================================================
# SEGMENTATION BACKEND
# ===================================================================


def test_community_segmentation_backend():
    """Community tier uses NV-Segment-CT."""
    s = ImagingSettings(DEPLOYMENT_TIER="community", ANTHROPIC_API_KEY=None)
    assert s.segmentation_backend == "nv-segment-ct"


def test_enterprise_segmentation_backend():
    """Enterprise tier uses VISTA-3D NIM."""
    s = ImagingSettings(DEPLOYMENT_TIER="enterprise", ANTHROPIC_API_KEY=None)
    assert s.segmentation_backend == "vista3d-nim"


def test_research_segmentation_backend():
    """Research tier uses NV-Segment-CTMR."""
    s = ImagingSettings(DEPLOYMENT_TIER="research", ANTHROPIC_API_KEY=None)
    assert s.segmentation_backend == "nv-segment-ctmr"


# ===================================================================
# CXR BACKEND
# ===================================================================


def test_community_cxr_backend():
    """Community tier uses DenseNet-121 for CXR classification."""
    s = ImagingSettings(DEPLOYMENT_TIER="community", ANTHROPIC_API_KEY=None)
    assert s.cxr_backend == "densenet121"


def test_enterprise_cxr_backend():
    """Enterprise tier also uses DenseNet-121 (open model)."""
    s = ImagingSettings(DEPLOYMENT_TIER="enterprise", ANTHROPIC_API_KEY=None)
    assert s.cxr_backend == "densenet121"


def test_research_cxr_backend():
    """Research tier uses NV-Reason-CXR-3B (noncommercial)."""
    s = ImagingSettings(DEPLOYMENT_TIER="research", ANTHROPIC_API_KEY=None)
    assert s.cxr_backend == "nv-reason-cxr-3b"


# ===================================================================
# EMBEDDING & STORAGE BACKENDS
# ===================================================================


def test_community_embedding_backend():
    """Community tier uses BGE-small-en-v1.5."""
    s = ImagingSettings(DEPLOYMENT_TIER="community", ANTHROPIC_API_KEY=None)
    assert s.embedding_backend == "bge-small-en-v1.5"


def test_enterprise_embedding_backend():
    """Enterprise tier uses NeMo Retriever."""
    s = ImagingSettings(DEPLOYMENT_TIER="enterprise", ANTHROPIC_API_KEY=None)
    assert s.embedding_backend == "nemo-retriever"


def test_community_storage_backend():
    """Community tier uses local storage."""
    s = ImagingSettings(DEPLOYMENT_TIER="community", ANTHROPIC_API_KEY=None)
    assert s.storage_backend == "local"


def test_enterprise_storage_backend():
    """Enterprise tier uses the AI platform."""
    s = ImagingSettings(DEPLOYMENT_TIER="enterprise", ANTHROPIC_API_KEY=None)
    assert s.storage_backend == "vast-aios"


# ===================================================================
# WEIGHT SUM
# ===================================================================


def test_weight_sum():
    """All collection weights (including radiomics + reports) sum to ~1.0."""
    s = ImagingSettings(ANTHROPIC_API_KEY=None)
    total = (
        s.WEIGHT_LITERATURE
        + s.WEIGHT_TRIALS
        + s.WEIGHT_FINDINGS
        + s.WEIGHT_PROTOCOLS
        + s.WEIGHT_DEVICES
        + s.WEIGHT_ANATOMY
        + s.WEIGHT_BENCHMARKS
        + s.WEIGHT_GUIDELINES
        + s.WEIGHT_REPORT_TEMPLATES
        + s.WEIGHT_DATASETS
        + s.WEIGHT_GENOMIC
        + s.WEIGHT_RADIOMICS
        + s.WEIGHT_REPORTS
    )
    assert abs(total - 1.0) < 0.02, f"Weight sum {total:.4f} deviates from 1.0"


# ===================================================================
# ENV VAR OVERRIDE
# ===================================================================


def test_tier_from_env():
    """IMAGING_DEPLOYMENT_TIER env var overrides the default."""
    with patch.dict(os.environ, {"IMAGING_DEPLOYMENT_TIER": "enterprise"}):
        s = ImagingSettings(ANTHROPIC_API_KEY=None)
        assert s.DEPLOYMENT_TIER == "enterprise"
        assert s.segmentation_backend == "vista3d-nim"
