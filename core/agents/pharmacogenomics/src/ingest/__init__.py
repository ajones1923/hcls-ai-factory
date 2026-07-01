"""Pharmacogenomics Intelligence Agent -- data ingest pipelines.

Provides parsers for ingesting PGx data from multiple authoritative
sources into Milvus vector collections:

  - CPICGuidelineParser       -- CPIC clinical prescribing guidelines
  - PharmVarParser            -- PharmVar star allele definitions
  - PharmGKBParser            -- PharmGKB clinical/variant annotations
  - FDALabelParser            -- FDA pharmacogenomic drug labeling
  - PopulationFrequencyParser -- Population allele frequency data
  - PubMedPGxParser           -- PubMed PGx literature
  - ClinicalTrialsPGxParser   -- ClinicalTrials.gov PGx trials

Author: Adam Jones
Date: March 2026
"""

from .base import BaseIngestPipeline
from .clinical_trials_parser import ClinicalTrialsPGxParser
from .cpic_parser import CPICGuidelineParser
from .fda_label_parser import FDALabelParser
from .pharmgkb_parser import PharmGKBParser
from .pharmvar_parser import PharmVarParser
from .population_parser import PopulationFrequencyParser
from .pubmed_parser import PubMedPGxParser

__all__ = [
    "BaseIngestPipeline",
    "CPICGuidelineParser",
    "PharmVarParser",
    "PharmGKBParser",
    "FDALabelParser",
    "PopulationFrequencyParser",
    "PubMedPGxParser",
    "ClinicalTrialsPGxParser",
]
