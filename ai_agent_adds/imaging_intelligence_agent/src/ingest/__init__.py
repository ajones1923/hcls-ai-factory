"""Data ingest pipelines for Imaging Intelligence Agent.

Provides the BaseIngestPipeline abstract class and all 10 concrete
ingest pipeline implementations for populating the Milvus collections.

Author: Adam Jones
Date: February 2026
"""

from .base import BaseIngestPipeline
from .literature_parser import PubMedImagingIngestPipeline
from .clinical_trials_parser import ImagingTrialsIngestPipeline
from .finding_parser import FindingIngestPipeline
from .protocol_parser import ProtocolIngestPipeline
from .device_parser import DeviceIngestPipeline
from .anatomy_parser import AnatomyIngestPipeline
from .benchmark_parser import BenchmarkIngestPipeline
from .guideline_parser import GuidelineIngestPipeline
from .report_template_parser import ReportTemplateIngestPipeline
from .dataset_parser import DatasetIngestPipeline
from .dicom_watcher import DicomWatcher

__all__ = [
    "BaseIngestPipeline",
    "PubMedImagingIngestPipeline",
    "ImagingTrialsIngestPipeline",
    "FindingIngestPipeline",
    "ProtocolIngestPipeline",
    "DeviceIngestPipeline",
    "AnatomyIngestPipeline",
    "BenchmarkIngestPipeline",
    "GuidelineIngestPipeline",
    "ReportTemplateIngestPipeline",
    "DatasetIngestPipeline",
    "DicomWatcher",
]
