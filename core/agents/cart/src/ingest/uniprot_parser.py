"""UniProt protein data ingest pipeline for CAR-T Intelligence Agent.

Fetches protein records for CAR-T target antigens via the UniProt REST API,
parses JSON responses into SequenceRecord models, and stores embeddings
in the cart_sequences Milvus collection.

UniProt REST API docs: https://www.uniprot.org/help/api_queries

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import CARTCollectionManager
from src.models import SequenceRecord

from .base import BaseIngestPipeline


# UniProt REST API search endpoint
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# Default query for CAR-T target proteins
DEFAULT_QUERY = "CD19 OR BCMA OR TNFRSF17 OR MS4A1"

# Mapping from UniProt gene/protein names to CAR-T target antigen names
_GENE_TO_ANTIGEN: Dict[str, str] = {
    "CD19": "CD19",
    "MS4A1": "CD20",
    "CD20": "CD20",
    "TNFRSF17": "BCMA",
    "BCMA": "BCMA",
    "CD22": "CD22",
    "CD33": "CD33",
    "CD38": "CD38",
    "CD123": "CD123",
    "IL3RA": "CD123",
    "CD269": "BCMA",
    "SDC1": "CD138",
    "CD138": "CD138",
    "GPRC5D": "GPRC5D",
    "FCRL5": "FcRH5",
    "FCRH5": "FcRH5",
    "SLAMF7": "CS1/SLAMF7",
    "CS1": "CS1/SLAMF7",
}

# Maximum retry attempts for network requests
MAX_RETRIES = 3

# Delay between requests to be respectful to UniProt servers
REQUEST_DELAY_SEC = 0.5


class UniProtIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for UniProt protein data relevant to CAR-T targets.

    Fetches protein records for key CAR-T target antigens (CD19, BCMA/TNFRSF17,
    CD20/MS4A1, etc.) from the UniProt REST API, parses the JSON responses into
    SequenceRecord models, and stores embeddings in the cart_sequences Milvus
    collection.

    Provides structural and functional protein data that enriches the agent's
    understanding of CAR-T binding targets: protein domains, binding sites,
    species information, and sequence annotations.

    Usage:
        pipeline = UniProtIngestPipeline(collection_manager, embedder)
        count = pipeline.run(max_results=50)
    """

    COLLECTION_NAME = "cart_sequences"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
    ):
        """Initialize the UniProt ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
        """
        super().__init__(collection_manager, embedder)

    def fetch(
        self,
        query: str = DEFAULT_QUERY,
        max_results: int = 50,
        format: str = "json",
    ) -> List[Dict[str, Any]]:
        """Fetch protein records from the UniProt REST API.

        Uses the UniProt search endpoint with a query targeting CAR-T-relevant
        proteins. Retrieves results in JSON format with full protein annotations.

        API endpoint: https://rest.uniprot.org/uniprotkb/search
        Query parameters:
            query  — UniProt query string (supports field searches)
            format — response format (json)
            size   — maximum number of results

        Args:
            query: UniProt search query string for CAR-T target proteins.
            max_results: Maximum number of protein entries to retrieve.
            format: Response format (default: json).

        Returns:
            List of protein entry JSON objects from the UniProt API.
        """
        params: Dict[str, Any] = {
            "query": query,
            "format": format,
            "size": min(max_results, 500),  # UniProt caps at 500 per page
        }

        all_entries: List[Dict[str, Any]] = []
        next_url: Optional[str] = None
        page_num = 0

        # First request
        response = self._request_with_retry(UNIPROT_SEARCH_URL, params)
        if response is None:
            logger.warning("UniProt API unavailable; returning empty result set")
            return []

        data = response.json()
        results = data.get("results", [])
        all_entries.extend(results)
        page_num += 1

        logger.info(
            f"UniProt: fetched page {page_num} ({len(results)} entries), "
            f"total {len(all_entries)} so far"
        )

        # Handle pagination via Link header
        while len(all_entries) < max_results:
            link_header = response.headers.get("Link", "")
            next_url = self._extract_next_link(link_header)
            if not next_url:
                break

            time.sleep(REQUEST_DELAY_SEC)

            response = self._request_with_retry(next_url, params=None)
            if response is None:
                logger.warning("UniProt pagination request failed; stopping")
                break

            data = response.json()
            results = data.get("results", [])
            if not results:
                break

            all_entries.extend(results)
            page_num += 1

            logger.info(
                f"UniProt: fetched page {page_num} ({len(results)} entries), "
                f"total {len(all_entries)} so far"
            )

        trimmed = all_entries[:max_results]
        logger.info(f"UniProt fetch complete: {len(trimmed)} protein entries")
        return trimmed

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[SequenceRecord]:
        """Parse UniProt protein entries into SequenceRecord models.

        Maps the UniProt JSON response structure to SequenceRecord fields:
            - proteinDescription.recommendedName.fullName -> construct_name
            - organism.scientificName -> species_origin
            - sequence.value -> structural_notes (summary)
            - features[] -> binding domain info, scfv_clone context
            - genes[].geneName.value -> target_antigen mapping

        Args:
            raw_data: List of UniProt protein entry JSON objects.

        Returns:
            List of validated SequenceRecord model instances.
        """
        records: List[SequenceRecord] = []

        for idx, entry in enumerate(raw_data):
            try:
                record = self._parse_single_entry(entry, idx)
                if record is not None:
                    records.append(record)
            except Exception as exc:
                accession = entry.get("primaryAccession", f"index-{idx}")
                logger.warning(f"Failed to parse UniProt entry {accession}: {exc}")
                continue

        logger.info(
            f"Parsed {len(records)} SequenceRecord instances "
            f"from {len(raw_data)} UniProt entries"
        )
        return records

    def _parse_single_entry(
        self, entry: Dict[str, Any], idx: int
    ) -> Optional[SequenceRecord]:
        """Parse a single UniProt protein entry into a SequenceRecord.

        Args:
            entry: Single UniProt protein entry JSON object.
            idx: Index in the batch for fallback ID generation.

        Returns:
            SequenceRecord instance, or None if the entry lacks required data.
        """
        # --- Accession (primary ID) ---
        accession = entry.get("primaryAccession", "")
        if not accession:
            return None
        record_id = f"UNIPROT-{accession}"

        # --- Protein description / construct name ---
        protein_desc = entry.get("proteinDescription", {})
        construct_name = self._extract_protein_name(protein_desc)

        # --- Gene name -> target antigen mapping ---
        genes = entry.get("genes", [])
        gene_name = ""
        target_antigen = ""
        for gene in genes:
            gn = gene.get("geneName", {})
            gene_name = gn.get("value", "")
            if gene_name:
                target_antigen = _GENE_TO_ANTIGEN.get(
                    gene_name.upper(), gene_name
                )
                break

        # --- Organism / species ---
        organism = entry.get("organism", {})
        scientific_name = organism.get("scientificName", "")
        species_origin = self._map_species_origin(scientific_name)

        # --- Sequence info ---
        sequence_info = entry.get("sequence", {})
        seq_length = sequence_info.get("length", 0)
        seq_mass = sequence_info.get("molWeight", 0)
        seq_value = sequence_info.get("value", "")

        # --- Features: extract binding domains, signal peptides, etc. ---
        features = entry.get("features", [])
        binding_info = self._extract_binding_domains(features)
        domain_info = self._extract_domain_info(features)

        # --- Function annotations (from comments) ---
        function_text = self._extract_function(entry.get("comments", []))

        # --- Build text summary ---
        text_parts = [
            f"UniProt protein record for {construct_name} ({accession}).",
            f"Gene: {gene_name}." if gene_name else "",
            f"Organism: {scientific_name}." if scientific_name else "",
            f"Length: {seq_length} aa, Mass: {seq_mass} Da." if seq_length else "",
        ]
        if function_text:
            text_parts.append(f"Function: {function_text}")
        if binding_info:
            text_parts.append(f"Binding domains: {binding_info}")
        if domain_info:
            text_parts.append(f"Domains: {domain_info}")

        text_summary = " ".join(p for p in text_parts if p)
        if len(text_summary) > 2900:
            text_summary = text_summary[:2897] + "..."

        # --- Structural notes: sequence excerpt + domain summary ---
        structural_parts = []
        if seq_value:
            # Include first 100 residues as a sequence excerpt
            excerpt = seq_value[:100]
            structural_parts.append(f"N-terminal sequence (first 100 aa): {excerpt}")
        if domain_info:
            structural_parts.append(f"Domain architecture: {domain_info}")
        if binding_info:
            structural_parts.append(f"Binding regions: {binding_info}")

        structural_notes = "; ".join(structural_parts)
        if len(structural_notes) > 990:
            structural_notes = structural_notes[:987] + "..."

        # --- Framework from protein type ---
        framework = self._extract_framework(protein_desc, features)

        # --- Immunogenicity risk (based on species) ---
        immunogenicity_risk = ""
        if "homo sapiens" in scientific_name.lower():
            immunogenicity_risk = "low"
        elif "mus musculus" in scientific_name.lower():
            immunogenicity_risk = "high"

        return SequenceRecord(
            id=record_id[:100],
            text_summary=text_summary,
            construct_name=construct_name[:200],
            target_antigen=target_antigen[:100],
            scfv_clone="",  # scFv clone info not directly in UniProt target entries
            binding_affinity_kd="",  # Not typically in UniProt
            heavy_chain_vregion="",
            light_chain_vregion="",
            framework=framework[:100],
            species_origin=species_origin[:30],
            immunogenicity_risk=immunogenicity_risk[:20],
            structural_notes=structural_notes,
        )

    @staticmethod
    def _extract_protein_name(protein_desc: Dict[str, Any]) -> str:
        """Extract the recommended or submitted protein name.

        Args:
            protein_desc: proteinDescription object from the UniProt entry.

        Returns:
            Protein full name string.
        """
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        if isinstance(full_name, dict):
            name = full_name.get("value", "")
        elif isinstance(full_name, str):
            name = full_name
        else:
            name = ""

        if name:
            return name

        # Fallback to submittedName
        submitted = protein_desc.get("submissionNames", [])
        if submitted:
            sub_name = submitted[0].get("fullName", {})
            if isinstance(sub_name, dict):
                return sub_name.get("value", "Unknown protein")
            return str(sub_name) if sub_name else "Unknown protein"

        # Fallback to alternativeNames
        alt_names = protein_desc.get("alternativeNames", [])
        if alt_names:
            alt_name = alt_names[0].get("fullName", {})
            if isinstance(alt_name, dict):
                return alt_name.get("value", "Unknown protein")
            return str(alt_name) if alt_name else "Unknown protein"

        return "Unknown protein"

    @staticmethod
    def _map_species_origin(scientific_name: str) -> str:
        """Map organism scientific name to a species origin category.

        Args:
            scientific_name: e.g., "Homo sapiens", "Mus musculus".

        Returns:
            Species origin string compatible with SequenceRecord.
        """
        name_lower = scientific_name.lower()
        if "homo sapiens" in name_lower:
            return "fully_human"
        if "mus musculus" in name_lower:
            return "murine"
        if "rattus" in name_lower:
            return "murine"
        if "cricetulus" in name_lower or "hamster" in name_lower:
            return "murine"
        return scientific_name[:30] if scientific_name else "unknown"

    @staticmethod
    def _extract_binding_domains(features: List[Dict[str, Any]]) -> str:
        """Extract binding site and domain information from UniProt features.

        Focuses on features relevant to CAR-T targeting: binding sites,
        active sites, and regions involved in receptor-ligand interactions.

        Args:
            features: List of feature objects from the UniProt entry.

        Returns:
            Semicolon-separated string of binding domain descriptions.
        """
        binding_features = []
        relevant_types = {"Binding site", "Active site", "Region", "Domain"}

        for feature in features:
            feat_type = feature.get("type", "")
            if feat_type not in relevant_types:
                continue

            description = feature.get("description", "")
            location = feature.get("location", {})

            # Extract position range
            start = location.get("start", {}).get("value", "?")
            end = location.get("end", {}).get("value", "?")
            position = f"{start}-{end}"

            if description:
                binding_features.append(f"{feat_type} ({position}): {description}")
            else:
                binding_features.append(f"{feat_type} ({position})")

            if len(binding_features) >= 5:
                break

        return "; ".join(binding_features)

    @staticmethod
    def _extract_domain_info(features: List[Dict[str, Any]]) -> str:
        """Extract protein domain annotations from UniProt features.

        Args:
            features: List of feature objects from the UniProt entry.

        Returns:
            Comma-separated string of domain names with positions.
        """
        domains = []
        for feature in features:
            if feature.get("type") not in ("Domain", "Topological domain"):
                continue

            description = feature.get("description", "")
            location = feature.get("location", {})
            start = location.get("start", {}).get("value", "?")
            end = location.get("end", {}).get("value", "?")

            if description:
                domains.append(f"{description} ({start}-{end})")
            if len(domains) >= 8:
                break

        return ", ".join(domains)

    @staticmethod
    def _extract_function(comments: List[Dict[str, Any]]) -> str:
        """Extract protein function description from UniProt comments.

        Args:
            comments: List of comment objects from the UniProt entry.

        Returns:
            Function description text, truncated to 500 characters.
        """
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    value = texts[0].get("value", "")
                    if len(value) > 500:
                        value = value[:497] + "..."
                    return value
        return ""

    @staticmethod
    def _extract_framework(
        protein_desc: Dict[str, Any],
        features: List[Dict[str, Any]],
    ) -> str:
        """Infer protein framework type from description and features.

        Args:
            protein_desc: proteinDescription object.
            features: List of feature objects.

        Returns:
            Framework string (e.g., "Type I transmembrane", "Immunoglobulin").
        """
        # Check for Ig-like domains in features
        ig_domains = []
        for feature in features:
            desc = feature.get("description", "").lower()
            if "ig-like" in desc or "immunoglobulin" in desc:
                ig_domains.append(feature.get("description", ""))

        if ig_domains:
            return f"Immunoglobulin superfamily ({len(ig_domains)} Ig-like domains)"

        # Check protein name for type hints
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        name_str = full_name.get("value", "") if isinstance(full_name, dict) else str(full_name)

        if "receptor" in name_str.lower():
            return "Transmembrane receptor"
        if "antigen" in name_str.lower():
            return "Cell surface antigen"

        return ""

    @staticmethod
    def _extract_next_link(link_header: str) -> Optional[str]:
        """Extract the 'next' pagination URL from the Link header.

        UniProt uses RFC 5988 Link headers for pagination:
            <url>; rel="next"

        Args:
            link_header: Link header string from the HTTP response.

        Returns:
            Next page URL, or None if not found.
        """
        if not link_header:
            return None

        # Parse Link header format: <url>; rel="next"
        parts = link_header.split(",")
        for part in parts:
            segments = part.strip().split(";")
            if len(segments) >= 2:
                rel = segments[1].strip()
                if 'rel="next"' in rel:
                    url = segments[0].strip().strip("<>")
                    return url
        return None

    @staticmethod
    def _request_with_retry(
        url: str,
        params: Optional[Dict[str, Any]],
        max_retries: int = MAX_RETRIES,
    ) -> Optional[requests.Response]:
        """Make an HTTP GET request with retry logic.

        Args:
            url: Request URL.
            params: Query parameters (None for paginated follow-up requests
                where the URL already contains parameters).
            max_retries: Maximum number of retry attempts.

        Returns:
            requests.Response on success, None on failure.
        """
        for attempt in range(1, max_retries + 1):
            try:
                kwargs: Dict[str, Any] = {
                    "headers": {"Accept": "application/json"},
                    "timeout": 30,
                }
                if params is not None:
                    kwargs["params"] = params

                response = requests.get(url, **kwargs)
                response.raise_for_status()
                return response
            except requests.exceptions.HTTPError as exc:
                status_code = exc.response.status_code if exc.response else 0
                if status_code == 429:
                    wait = 2 ** attempt
                    logger.warning(
                        f"UniProt rate limited (429); retrying in {wait}s "
                        f"(attempt {attempt}/{max_retries})"
                    )
                    time.sleep(wait)
                    continue
                logger.warning(
                    f"UniProt HTTP error {status_code} on attempt "
                    f"{attempt}/{max_retries}: {exc}"
                )
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return None
            except requests.exceptions.RequestException as exc:
                logger.warning(
                    f"UniProt request failed on attempt "
                    f"{attempt}/{max_retries}: {exc}"
                )
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return None

        return None

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        max_results: int = 50,
        **fetch_kwargs,
    ) -> int:
        """Execute the full UniProt ingest pipeline: fetch -> parse -> embed_and_store.

        Args:
            collection_name: Target Milvus collection (defaults to 'cart_sequences').
            batch_size: Batch size for embedding and insertion.
            max_results: Maximum number of protein entries to fetch.
            **fetch_kwargs: Additional keyword arguments passed to fetch().

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting UniProt ingest pipeline -> {target}")

        raw = self.fetch(max_results=max_results, **fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw protein entries from UniProt")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} SequenceRecord instances")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"UniProt ingest complete: {count} records into {target}")
        return count
