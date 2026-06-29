"""FDA regulatory document parser for the Clinical Trial Intelligence Agent.

Parses FDA approval packages, guidance documents, and regulatory decisions
for embedding and storage in Milvus.

Data sources:
  - FDA Drugs@FDA API
  - FDA openFDA drug label API
  - Curated regulatory milestone database

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, List, Optional

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)

# ===================================================================
# CONSTANTS
# ===================================================================

OPENFDA_DRUG_URL = "https://api.fda.gov/drug/drugsfda.json"
OPENFDA_LABEL_URL = "https://api.fda.gov/drug/label.json"
RATE_LIMIT_DELAY = 0.25  # 4 requests/second

# Curated FDA regulatory milestones for seeding
REGULATORY_MILESTONES: List[Dict[str, Any]] = [
    {
        "drug": "Pembrolizumab (Keytruda)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Metastatic melanoma (unresectable)",
        "date": "2014-09-04",
        "mechanism": "PD-1 checkpoint inhibitor",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review",
        "significance": "First PD-1 inhibitor approved",
    },
    {
        "drug": "Nivolumab (Opdivo)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Metastatic melanoma",
        "date": "2014-12-22",
        "mechanism": "PD-1 checkpoint inhibitor",
        "regulatory_pathway": "Breakthrough Therapy",
        "significance": "Second PD-1 inhibitor approved, expanded to multiple indications",
    },
    {
        "drug": "Sacubitril/Valsartan (Entresto)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Heart failure with reduced ejection fraction",
        "date": "2015-07-07",
        "mechanism": "Angiotensin receptor-neprilysin inhibitor (ARNI)",
        "regulatory_pathway": "Priority Review",
        "significance": "First ARNI approved, paradigm shift in HFrEF treatment",
    },
    {
        "drug": "Tisagenlecleucel (Kymriah)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Pediatric/young adult relapsed/refractory B-cell ALL",
        "date": "2017-08-30",
        "mechanism": "CAR-T cell therapy",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review",
        "significance": "First FDA-approved CAR-T cell therapy",
    },
    {
        "drug": "Dapagliflozin (Farxiga)",
        "agency": "FDA",
        "decision": "sNDA Approval",
        "indication": "Heart failure with reduced ejection fraction",
        "date": "2020-05-05",
        "mechanism": "SGLT2 inhibitor",
        "regulatory_pathway": "Priority Review",
        "significance": "First SGLT2i approved for HF regardless of diabetes status",
    },
    {
        "drug": "Empagliflozin (Jardiance)",
        "agency": "FDA",
        "decision": "sNDA Approval",
        "indication": "Heart failure",
        "date": "2022-02-24",
        "mechanism": "SGLT2 inhibitor",
        "regulatory_pathway": "Standard Review",
        "significance": "Expanded to include HFpEF in addition to HFrEF",
    },
    {
        "drug": "Trastuzumab Deruxtecan (Enhertu)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "HER2-positive metastatic breast cancer",
        "date": "2019-12-20",
        "mechanism": "Antibody-drug conjugate (ADC)",
        "regulatory_pathway": "Breakthrough Therapy, Accelerated Approval",
        "significance": "Next-generation ADC with novel payload technology",
    },
    {
        "drug": "Osimertinib (Tagrisso)",
        "agency": "FDA",
        "decision": "sNDA Approval",
        "indication": "Adjuvant EGFR-mutant NSCLC",
        "date": "2020-12-18",
        "mechanism": "Third-generation EGFR TKI",
        "regulatory_pathway": "Priority Review",
        "significance": "First adjuvant targeted therapy in early-stage NSCLC",
    },
    {
        "drug": "Sotorasib (Lumakras)",
        "agency": "FDA",
        "decision": "Accelerated Approval",
        "indication": "KRAS G12C-mutant NSCLC",
        "date": "2021-05-28",
        "mechanism": "KRAS G12C inhibitor",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Accelerated Approval",
        "significance": "First KRAS-targeted therapy approved",
    },
    {
        "drug": "Mavacamten (Camzyos)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Symptomatic obstructive hypertrophic cardiomyopathy",
        "date": "2022-04-28",
        "mechanism": "Cardiac myosin inhibitor",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review",
        "significance": "First cardiac myosin inhibitor, novel mechanism for HCM",
    },
    {
        "drug": "Tafamidis (Vyndaqel/Vyndamax)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Transthyretin amyloid cardiomyopathy",
        "date": "2019-05-03",
        "mechanism": "Transthyretin stabilizer",
        "regulatory_pathway": "Priority Review",
        "significance": "First disease-modifying therapy for ATTR cardiomyopathy",
    },
    {
        "drug": "Inclisiran (Leqvio)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Heterozygous familial hypercholesterolemia, ASCVD",
        "date": "2021-12-22",
        "mechanism": "PCSK9 siRNA",
        "regulatory_pathway": "Priority Review",
        "significance": "First siRNA-based LDL-lowering therapy, twice-yearly dosing",
    },
    {
        "drug": "Resmetirom (Rezdiffra)",
        "agency": "FDA",
        "decision": "Accelerated Approval",
        "indication": "Non-cirrhotic NASH with moderate-to-advanced liver fibrosis",
        "date": "2024-03-14",
        "mechanism": "Thyroid hormone receptor-beta agonist",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Accelerated Approval",
        "significance": "First FDA-approved drug specifically for NASH",
    },
    {
        "drug": "Exagamglogene autotemcel (Casgevy)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Sickle cell disease and transfusion-dependent beta thalassemia",
        "date": "2023-12-08",
        "mechanism": "CRISPR/Cas9 gene-edited cell therapy",
        "regulatory_pathway": "Priority Review, Orphan Drug",
        "significance": "First CRISPR-based gene therapy approved by any regulatory agency",
    },
    {
        "drug": "Lecanemab (Leqembi)",
        "agency": "FDA",
        "decision": "Traditional Approval",
        "indication": "Early Alzheimer's disease",
        "date": "2023-07-06",
        "mechanism": "Anti-amyloid beta antibody",
        "regulatory_pathway": "Priority Review (converted from Accelerated Approval)",
        "significance": "First anti-amyloid antibody to receive traditional approval for AD",
    },
    {
        "drug": "Glofitamab (Columvi)",
        "agency": "FDA",
        "decision": "Accelerated Approval",
        "indication": "Relapsed/refractory diffuse large B-cell lymphoma",
        "date": "2023-06-15",
        "mechanism": "CD20xCD3 bispecific antibody",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Accelerated Approval",
        "significance": "First CD20xCD3 bispecific antibody approved for DLBCL",
    },
    {
        "drug": "Donanemab (Kisunla)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Early symptomatic Alzheimer's disease",
        "date": "2024-07-02",
        "mechanism": "Anti-amyloid beta (N3pG) antibody",
        "regulatory_pathway": "Priority Review",
        "significance": "Second anti-amyloid antibody approved for AD, targets pyroglutamate A-beta",
    },
    {
        "drug": "Teplizumab (Tzield)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Delay onset of Stage 3 type 1 diabetes in adults and pediatric patients >= 8 years",
        "date": "2022-11-17",
        "mechanism": "Anti-CD3 monoclonal antibody",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review",
        "significance": "First therapy approved to delay onset of type 1 diabetes",
    },
    {
        "drug": "Repotrectinib (Augtyro)",
        "agency": "FDA",
        "decision": "Accelerated Approval",
        "indication": "ROS1-positive non-small cell lung cancer",
        "date": "2023-11-15",
        "mechanism": "Next-generation ROS1/TRK tyrosine kinase inhibitor",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Accelerated Approval",
        "significance": "Next-generation ROS1/TRK inhibitor designed to overcome resistance mutations",
    },
    {
        "drug": "Tovorafenib (Ojemda)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Relapsed/refractory pediatric low-grade glioma with BRAF alteration",
        "date": "2024-04-22",
        "mechanism": "Pan-RAF kinase inhibitor",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Rare Pediatric Disease",
        "significance": "First RAF inhibitor approved for pediatric low-grade glioma",
    },
    {
        "drug": "Retifanlimab (Zynyz)",
        "agency": "FDA",
        "decision": "Approval",
        "indication": "Metastatic or recurrent locally advanced Merkel cell carcinoma",
        "date": "2023-03-22",
        "mechanism": "PD-1 checkpoint inhibitor",
        "regulatory_pathway": "Priority Review",
        "significance": "First PD-1 inhibitor approved specifically for Merkel cell carcinoma",
    },
    {
        "drug": "Tarlatamab (Imdelltra)",
        "agency": "FDA",
        "decision": "Accelerated Approval",
        "indication": "Extensive-stage small cell lung cancer after prior platinum-based chemotherapy",
        "date": "2024-05-16",
        "mechanism": "DLL3-targeted bispecific T-cell engager",
        "regulatory_pathway": "Breakthrough Therapy, Priority Review, Accelerated Approval",
        "significance": "First DLL3-targeted bispecific T-cell engager approved for SCLC",
    },
]


# ===================================================================
# PARSER IMPLEMENTATION
# ===================================================================


class RegulatoryParser(BaseIngestParser):
    """Ingest parser for FDA regulatory documents and approval packages.

    Extracts agency, decision type, drug name, indication, date, and
    regulatory pathway information from FDA data sources and curated
    regulatory milestone databases.

    Usage::

        parser = RegulatoryParser()
        records, stats = parser.run(
            drug_names=["pembrolizumab", "osimertinib"],
            max_results=50,
        )
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
        api_key: Optional[str] = None,
    ) -> None:
        """Initialize the regulatory parser.

        Args:
            collection_manager: Optional Milvus collection manager.
            embedder: Optional embedding model.
            api_key: Optional openFDA API key.
        """
        super().__init__(
            source_name="regulatory",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key
        self._last_request_time: float = 0.0

    def _rate_limit(self) -> None:
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < RATE_LIMIT_DELAY:
            time.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = time.time()

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch regulatory data from FDA APIs and curated milestones.

        Args:
            drug_names: List of drug names to search for.
            max_results: Maximum number of results (default 50).
            include_milestones: Whether to include curated milestones (default True).

        Returns:
            List of raw regulatory record dictionaries.
        """
        drug_names = kwargs.get("drug_names", [])
        max_results = kwargs.get("max_results", 50)
        include_milestones = kwargs.get("include_milestones", True)

        all_records: List[Dict[str, Any]] = []

        # Include curated milestones
        if include_milestones:
            milestones = REGULATORY_MILESTONES.copy()
            if drug_names:
                # Filter milestones by drug names
                milestones = [
                    m for m in milestones
                    if any(
                        dn.lower() in m["drug"].lower()
                        for dn in drug_names
                    )
                ]
            all_records.extend(milestones)

        # Try to fetch from openFDA API
        if drug_names:
            try:
                import requests
                for drug_name in drug_names:
                    if len(all_records) >= max_results:
                        break
                    fda_records = self._fetch_openfda(drug_name, requests)
                    all_records.extend(fda_records)
            except ImportError:
                self.logger.warning(
                    "requests library not available -- using curated milestones only"
                )

        return all_records[:max_results]

    def _fetch_openfda(
        self, drug_name: str, requests_module: Any
    ) -> List[Dict[str, Any]]:
        """Fetch drug approval data from openFDA.

        Args:
            drug_name: Drug name to search for.
            requests_module: The requests library module.

        Returns:
            List of regulatory record dictionaries.
        """
        params: Dict[str, Any] = {
            "search": f'openfda.brand_name:"{drug_name}"',
            "limit": 5,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            self._rate_limit()
            response = requests_module.get(
                OPENFDA_DRUG_URL, params=params, timeout=30
            )
            if response.status_code == 404:
                return []
            response.raise_for_status()
            data = response.json()
        except Exception as exc:
            self.logger.warning(
                "openFDA fetch failed for '%s': %s", drug_name, exc
            )
            return []

        records: List[Dict[str, Any]] = []
        for result in data.get("results", []):
            products = result.get("products", [])
            submissions = result.get("submissions", [])

            for submission in submissions:
                sub_type = submission.get("submission_type", "")
                sub_status = submission.get("submission_status", "")
                sub_date = submission.get("submission_status_date", "")

                if sub_status.lower() in ("ap", "approved"):
                    record = {
                        "drug": drug_name,
                        "agency": "FDA",
                        "decision": f"{sub_type} Approval",
                        "indication": "",
                        "date": sub_date,
                        "mechanism": "",
                        "regulatory_pathway": sub_type,
                        "significance": "",
                        "source": "openfda",
                    }
                    if products:
                        record["indication"] = products[0].get(
                            "active_ingredients", ""
                        )
                    records.append(record)

        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw regulatory data into IngestRecord objects.

        Args:
            raw_data: List of regulatory record dictionaries.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for item in raw_data:
            try:
                drug = item.get("drug", "")
                agency = item.get("agency", "FDA")
                decision = item.get("decision", "")
                indication = item.get("indication", "")
                date = item.get("date", "")
                mechanism = item.get("mechanism", "")
                pathway = item.get("regulatory_pathway", "")
                significance = item.get("significance", "")

                # Compose text for embedding
                text_parts = [
                    f"Regulatory Decision: {agency} {decision}",
                    f"Drug: {drug}",
                    f"Indication: {indication}",
                    f"Date: {date}",
                ]
                if mechanism:
                    text_parts.append(f"Mechanism: {mechanism}")
                if pathway:
                    text_parts.append(f"Regulatory Pathway: {pathway}")
                if significance:
                    text_parts.append(f"Significance: {significance}")

                text = "\n".join(text_parts)

                # Generate a unique record ID
                record_id = f"REG:{agency}:{drug}:{date}".replace(" ", "_")

                metadata = {
                    "drug": drug,
                    "agency": agency,
                    "decision": decision,
                    "indication": indication,
                    "date": date,
                    "mechanism": mechanism,
                    "regulatory_pathway": pathway,
                    "significance": significance,
                    "source": item.get("source", "curated"),
                }

                records.append(
                    IngestRecord(
                        text=text,
                        metadata=metadata,
                        collection_name="trial_regulatory",
                        record_id=record_id,
                        source="regulatory",
                    )
                )
            except Exception as exc:
                self.logger.warning(
                    "Failed to parse regulatory record for '%s': %s",
                    item.get("drug", "?"), exc,
                )

        return records

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a regulatory IngestRecord.

        Args:
            record: The IngestRecord to validate.

        Returns:
            True if valid.
        """
        if not record.text or len(record.text.strip()) < 20:
            return False
        if not record.record_id:
            return False
        if not record.metadata.get("drug"):
            return False
        if not record.metadata.get("agency"):
            return False
        return True
