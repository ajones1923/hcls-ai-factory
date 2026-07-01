"""
Knowledge Base Versioning for the HCLS AI Factory.

Provides semantic versioning for all knowledge bases across all agents,
enabling reproducibility and staleness detection.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Any


@dataclass
class KnowledgeSource:
    """A single knowledge source with version tracking."""
    name: str
    version: str
    source_url: str = ""
    last_updated: str = ""
    record_count: int = 0
    description: str = ""


@dataclass
class KnowledgeManifest:
    """Version manifest for an agent's knowledge base.

    Each agent maintains a manifest documenting all knowledge sources,
    their versions, and last update dates. This enables:
    - Reproducibility: pin exact knowledge versions in reports
    - Staleness detection: alert when guidelines are outdated
    - Audit: trace every recommendation to its source data
    """
    agent_name: str
    manifest_version: str = "1.0.0"
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    sources: list[KnowledgeSource] = field(default_factory=list)

    def add_source(self, name: str, version: str, **kwargs: Any) -> None:
        self.sources.append(KnowledgeSource(name=name, version=version, **kwargs))

    def get_source(self, name: str) -> KnowledgeSource | None:
        for s in self.sources:
            if s.name == name:
                return s
        return None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def summary(self) -> str:
        """Human-readable summary for report footers."""
        lines = [f"Knowledge Base: {self.agent_name} v{self.manifest_version}"]
        for src in self.sources:
            lines.append(f"  - {src.name} v{src.version} ({src.last_updated})")
        return "\n".join(lines)


# Pre-built manifests for each agent
BIOMARKER_KNOWLEDGE = KnowledgeManifest(agent_name="Precision Biomarker Agent")
BIOMARKER_KNOWLEDGE.add_source("CPIC Guidelines", "1.22", source_url="https://cpicpgx.org", last_updated="2024-01", record_count=11, description="Clinical Pharmacogenetics Implementation Consortium guidelines for 11 genes")
BIOMARKER_KNOWLEDGE.add_source("PhenoAge Algorithm", "1.0", source_url="https://doi.org/10.18632/aging.101414", last_updated="2018-02", description="Levine et al. 2018 biological age calculation (PMID:29676998)")
BIOMARKER_KNOWLEDGE.add_source("GrimAge Surrogate", "1.0", source_url="https://doi.org/10.1186/s13148-020-00824-y", last_updated="2020-03", description="Hillary et al. 2020 DNA methylation proxy (PMID:32245506)")
BIOMARKER_KNOWLEDGE.add_source("ClinVar", "2024-12", source_url="https://www.ncbi.nlm.nih.gov/clinvar/", last_updated="2024-12", record_count=2700000, description="NCBI clinical variant database")
BIOMARKER_KNOWLEDGE.add_source("Ancestry Reference Ranges", "1.0", last_updated="2024-01", description="Population-specific biomarker ranges for 4 ancestries")

ONCOLOGY_KNOWLEDGE = KnowledgeManifest(agent_name="Precision Oncology Agent")
ONCOLOGY_KNOWLEDGE.add_source("OncoKB", "4.15", source_url="https://www.oncokb.org", last_updated="2024-06", description="MSK precision oncology knowledgebase")
ONCOLOGY_KNOWLEDGE.add_source("CIViC", "2024-03", source_url="https://civicdb.org", last_updated="2024-03", description="Clinical Interpretation of Variants in Cancer")
ONCOLOGY_KNOWLEDGE.add_source("NCCN Guidelines", "2024", source_url="https://www.nccn.org", last_updated="2024-01", description="National Comprehensive Cancer Network clinical practice guidelines")

CART_KNOWLEDGE = KnowledgeManifest(agent_name="CAR-T Intelligence Agent")
CART_KNOWLEDGE.add_source("CIBMTR Registry", "2024-Q2", source_url="https://cibmtr.org", last_updated="2024-06", description="Center for International Blood and Marrow Transplant Research")
CART_KNOWLEDGE.add_source("FDA Approved CAR-T Products", "2024", last_updated="2024-06", record_count=6, description="Kymriah, Yescarta, Tecartus, Breyanzi, Abecma, Carvykti")
CART_KNOWLEDGE.add_source("CAR-T Safety (FAERS)", "2024-Q1", source_url="https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html", last_updated="2024-03", description="FDA Adverse Event Reporting System for cell therapies")

IMAGING_KNOWLEDGE = KnowledgeManifest(agent_name="Imaging Intelligence Agent")
IMAGING_KNOWLEDGE.add_source("ACR Appropriateness Criteria", "2024", source_url="https://acsearch.acr.org/list", last_updated="2024-01", description="American College of Radiology imaging guidelines")
IMAGING_KNOWLEDGE.add_source("MONAI Medical Imaging Models", "1.3", source_url="https://monai.io", last_updated="2024-06", description="NVIDIA MONAI segmentation and classification models")
IMAGING_KNOWLEDGE.add_source("RadLex Ontology", "4.1", source_url="https://radlex.org", last_updated="2023-06", description="Radiology lexicon for structured reporting")

AUTOIMMUNE_KNOWLEDGE = KnowledgeManifest(agent_name="Precision Autoimmune Agent")
AUTOIMMUNE_KNOWLEDGE.add_source("ACR/EULAR Classification Criteria", "2019", source_url="https://www.rheumatology.org", last_updated="2019-09", description="SLE (2019), RA (2010), SSc (2013) classification criteria")
AUTOIMMUNE_KNOWLEDGE.add_source("HLA Disease Associations", "2024", last_updated="2024-01", description="HLA allele-disease associations across 50+ autoimmune conditions")
AUTOIMMUNE_KNOWLEDGE.add_source("Biologics PGx", "1.0", last_updated="2024-01", description="Pharmacogenomic predictors of biologic therapy response")
