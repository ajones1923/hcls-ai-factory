"""Data ingest parsers for Cardiology Intelligence Agent.

Provides seven domain-specific parsers that fetch, parse, and produce
IngestRecord instances for the 12 cardiology Milvus collections:

  - BaseIngestParser / IngestRecord  — shared base class and record type
  - PubMedCardioParser               — PubMed cardiovascular literature
  - ClinicalTrialsCardioParser       — ClinicalTrials.gov cardiovascular trials
  - GuidelineParser                  — ACC / AHA / ESC clinical practice guidelines
  - ImagingParser                    — Cardiac imaging (echo, CT, CMR, nuclear)
  - ECGParser                        — ECG / electrophysiology reference data
  - DeviceParser                     — FDA AI devices and implantable devices
  - HemodynamicsParser               — Hemodynamic parameters and cath lab protocols

Author: Adam Jones
Date: March 2026
"""

from .base import BaseIngestParser, IngestRecord
from .clinical_trials_parser import ClinicalTrialsCardioParser
from .device_parser import DeviceParser
from .ecg_parser import ECGParser
from .guideline_parser import GuidelineParser
from .hemodynamics_parser import HemodynamicsParser
from .imaging_parser import ImagingParser
from .pubmed_parser import PubMedCardioParser

__all__ = [
    "BaseIngestParser",
    "IngestRecord",
    "PubMedCardioParser",
    "ClinicalTrialsCardioParser",
    "GuidelineParser",
    "ImagingParser",
    "ECGParser",
    "DeviceParser",
    "HemodynamicsParser",
]
