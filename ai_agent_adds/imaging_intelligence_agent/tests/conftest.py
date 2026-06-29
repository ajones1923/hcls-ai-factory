"""Shared pytest fixtures for Imaging Intelligence Agent test suite.

Provides mock embedder, LLM client, collection manager, NIM service manager,
and sample search results so that tests run without Milvus, NIM services,
or any network access.

Author: Adam Jones
Date: February 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, PropertyMock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so ``from src.…`` imports work
# regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.models import (  # noqa: E402
    AgentResponse,
    CrossCollectionResult,
    FindingSeverity,
    SearchHit,
    WorkflowResult,
    WorkflowStatus,
)


# ===================================================================
# MOCK EMBEDDER
# ===================================================================


@pytest.fixture
def mock_embedder():
    """Return a mock SentenceTransformer that produces 384-dim numpy vectors.

    Supports both single-string and batch encoding:
      - embedder.encode("text")           -> np.ndarray shape (384,)
      - embedder.encode(["a", "b"])       -> np.ndarray shape (2, 384)
    """
    embedder = MagicMock()

    def _encode(texts, normalize_embeddings=True, **kwargs):
        if isinstance(texts, str):
            return np.random.randn(384).astype(np.float32)
        return np.random.randn(len(texts), 384).astype(np.float32)

    embedder.encode = MagicMock(side_effect=_encode)
    return embedder


# ===================================================================
# MOCK LLM CLIENT
# ===================================================================


@pytest.fixture
def mock_llm_client():
    """Return a mock LlamaLLMClient with generate() returning a template answer."""
    client = MagicMock()
    client.generate.return_value = (
        "Based on the available imaging evidence, the findings are consistent "
        "with normal anatomy without acute pathology. Clinical correlation "
        "recommended."
    )
    client.generate_stream.return_value = iter([
        "Based on ", "the available ", "imaging evidence, ",
        "the findings ", "are consistent ", "with normal anatomy.",
    ])
    client.is_available.return_value = False
    client.mock_enabled = True
    return client


# ===================================================================
# MOCK COLLECTION MANAGER
# ===================================================================


IMAGING_COLLECTION_NAMES = [
    "imaging_literature",
    "imaging_trials",
    "imaging_findings",
    "imaging_protocols",
    "imaging_devices",
    "imaging_anatomy",
    "imaging_benchmarks",
    "imaging_guidelines",
    "imaging_report_templates",
    "imaging_datasets",
    "genomic_evidence",
]


@pytest.fixture
def mock_collection_manager():
    """Return a MagicMock collection manager with sane defaults.

    - search()                -> list of 3 sample result dicts
    - get_collection_stats()  -> dummy counts for all 11 collections
    - connect() / disconnect() -> no-ops
    """
    manager = MagicMock()

    manager.search.return_value = [
        {
            "id": "lit-001",
            "score": 0.85,
            "text_chunk": "AI-based hemorrhage detection achieves 95% sensitivity.",
            "modality": "ct",
            "body_region": "head",
        },
        {
            "id": "lit-002",
            "score": 0.78,
            "text_chunk": "VISTA-3D segments 127 anatomical structures from CT.",
            "modality": "ct",
            "body_region": "chest",
        },
        {
            "id": "lit-003",
            "score": 0.72,
            "text_chunk": "Lung-RADS v2022 classifies solid nodules by diameter.",
            "modality": "ct",
            "body_region": "chest",
        },
    ]

    manager.get_collection_stats.return_value = {
        name: 100 for name in IMAGING_COLLECTION_NAMES
    }

    manager.connect.return_value = None
    manager.disconnect.return_value = None

    return manager


# ===================================================================
# MOCK NIM SERVICES
# ===================================================================


@pytest.fixture
def mock_nim_services():
    """Return a mock NIMServiceManager with all clients returning mock results."""
    manager = MagicMock()

    # VISTA-3D
    vista3d = MagicMock()
    vista3d.is_available.return_value = False
    vista3d.mock_enabled = True
    vista3d.get_status.return_value = "mock"
    vista3d.segment.return_value = MagicMock(
        classes_detected=["liver", "spleen", "right_kidney"],
        volumes={"liver": 1450.0, "spleen": 180.0, "right_kidney": 155.0},
        num_classes=3,
        inference_time_ms=3500.0,
        is_mock=True,
    )
    vista3d.get_supported_classes.return_value = [
        "liver", "spleen", "right_kidney", "left_kidney", "pancreas",
    ]
    type(manager).vista3d = PropertyMock(return_value=vista3d)

    # MAISI
    maisi = MagicMock()
    maisi.is_available.return_value = False
    maisi.mock_enabled = True
    maisi.get_status.return_value = "mock"
    maisi.generate.return_value = MagicMock(
        resolution="512x512x512",
        body_region="chest",
        num_classes_annotated=127,
        generation_time_ms=45000.0,
        is_mock=True,
    )
    type(manager).maisi = PropertyMock(return_value=maisi)

    # VILA-M3
    vilam3 = MagicMock()
    vilam3.is_available.return_value = False
    vilam3.mock_enabled = True
    vilam3.get_status.return_value = "mock"
    vilam3.analyze_image.return_value = MagicMock(
        answer="No acute cardiopulmonary abnormality.",
        findings=["clear lungs", "normal cardiac silhouette"],
        confidence=0.85,
        is_mock=True,
    )
    vilam3.generate_report.return_value = (
        "TECHNIQUE: CT of the chest.\n"
        "FINDINGS: No acute abnormality.\n"
        "IMPRESSION: Normal chest CT."
    )
    type(manager).vilam3 = PropertyMock(return_value=vilam3)

    # LLM
    llm = MagicMock()
    llm.is_available.return_value = False
    llm.mock_enabled = True
    llm.get_status.return_value = "mock"
    llm.generate.return_value = "Mock clinical assessment based on imaging data."
    type(manager).llm = PropertyMock(return_value=llm)

    # Service management
    manager.check_all_services.return_value = {
        "vista3d": "mock",
        "maisi": "mock",
        "vila_m3": "mock",
        "llm": "mock",
    }
    manager.get_available_services.return_value = []
    manager.get_client.side_effect = lambda name: {
        "vista3d": vista3d,
        "maisi": maisi,
        "vila_m3": vilam3,
        "llm": llm,
    }[name]

    return manager


# ===================================================================
# SAMPLE SEARCH DATA
# ===================================================================


@pytest.fixture
def sample_search_hits():
    """Return a list of 5 SearchHit objects spanning different collections."""
    return [
        SearchHit(
            collection="imaging_literature",
            id="PMID:38001234",
            score=0.92,
            text="AI hemorrhage detection achieves 95% sensitivity on RSNA dataset.",
            metadata={
                "title": "Deep Learning for ICH Detection",
                "year": 2024,
                "modality": "ct",
            },
        ),
        SearchHit(
            collection="imaging_trials",
            id="NCT06001234",
            score=0.87,
            text="Phase 3 trial of AI-assisted lung cancer screening.",
            metadata={
                "phase": "Phase 3",
                "status": "Recruiting",
                "sponsor": "RSNA",
            },
        ),
        SearchHit(
            collection="imaging_findings",
            id="finding-ich-001",
            score=0.83,
            text="Intraparenchymal hemorrhage in right basal ganglia, 12.5 mL.",
            metadata={
                "finding_category": "hemorrhage",
                "severity": "urgent",
            },
        ),
        SearchHit(
            collection="imaging_guidelines",
            id="guideline-lungrads-001",
            score=0.78,
            text="ACR Lung-RADS v2022: solid nodule 8-15mm is category 4A.",
            metadata={
                "organization": "ACR",
                "year": 2022,
            },
        ),
        SearchHit(
            collection="imaging_benchmarks",
            id="bench-vista3d-001",
            score=0.71,
            text="VISTA-3D achieves mean Dice 0.91 across 127 classes on TotalSegmentator.",
            metadata={
                "model_name": "VISTA-3D",
                "metric_name": "Dice",
                "metric_value": 0.91,
            },
        ),
    ]


@pytest.fixture
def sample_evidence(sample_search_hits):
    """Return a CrossCollectionResult populated with 5 sample hits."""
    return CrossCollectionResult(
        query="What is the best AI model for hemorrhage detection?",
        hits=sample_search_hits,
        knowledge_context=(
            "## Pathology: Intracranial Hemorrhage\n"
            "- **ICD-10:** I62.9\n"
            "- **Body Region:** head\n"
            "- **Modalities:** ct, mri"
        ),
        total_collections_searched=11,
        search_time_ms=35.2,
    )


@pytest.fixture
def sample_workflow_result():
    """Return a WorkflowResult from CT head hemorrhage mock run."""
    return WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "hemorrhage",
                "description": (
                    "Intraparenchymal hemorrhage in right basal ganglia, "
                    "volume 12.5 mL, midline shift 3.2 mm, max thickness 8.1 mm"
                ),
                "severity": "urgent",
                "hemorrhage_type": "intraparenchymal",
                "location": "right basal ganglia",
            },
        ],
        measurements={
            "volume_ml": 12.5,
            "midline_shift_mm": 3.2,
            "max_thickness_mm": 8.1,
            "hounsfield_mean": 62.0,
            "hounsfield_max": 78.0,
        },
        classification="urgent_hemorrhage",
        severity=FindingSeverity.URGENT,
        nim_services_used=["3D U-Net (MONAI)", "VISTA-3D (optional)"],
        is_mock=True,
    )


@pytest.fixture
def sample_agent_response(sample_evidence, sample_workflow_result):
    """Return a complete AgentResponse for export testing."""
    return AgentResponse(
        question="What is the best AI model for hemorrhage detection?",
        answer=(
            "Based on the available evidence, 3D U-Net models trained with MONAI "
            "achieve 95% sensitivity for intracranial hemorrhage detection on CT."
        ),
        evidence=sample_evidence,
        workflow_results=[sample_workflow_result],
        nim_services_used=["vista3d"],
        knowledge_used=["pathology", "modality"],
        timestamp="2026-02-27T10:00:00",
    )


# ===================================================================
# SETTINGS
# ===================================================================


@pytest.fixture
def settings():
    """Return an ImagingSettings instance (no Milvus or NIM needed)."""
    from config.settings import ImagingSettings

    return ImagingSettings(
        MILVUS_HOST="localhost",
        MILVUS_PORT=19530,
        NIM_ALLOW_MOCK_FALLBACK=True,
        ANTHROPIC_API_KEY=None,
    )
