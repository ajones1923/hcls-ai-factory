"""CAR construct data ingest pipeline for CAR-T Intelligence Agent.

Parses CAR construct design data from reference files (JSON/CSV) and
FDA-approved product records into CARConstruct models, storing embeddings
in the cart_constructs Milvus collection.

Includes seed data for all 6 FDA-approved CAR-T products as of Feb 2026.

Author: Adam Jones
Date: February 2026
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import (
    CARConstruct,
    CARGeneration,
    FDAStatus,
)

from .base import BaseIngestPipeline


class ConstructIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR construct design data.

    Reads construct definitions from a reference JSON/CSV file and/or
    the built-in FDA-approved seed data, converts them into CARConstruct
    models, and stores embeddings in the cart_constructs collection.

    Usage:
        pipeline = ConstructIngestPipeline(collection_manager, embedder)
        # Ingest FDA-approved products only:
        count = pipeline.run(include_fda_seed=True)
        # Ingest from a reference file:
        count = pipeline.run(reference_file="/path/to/constructs.json")
    """

    COLLECTION_NAME = "cart_constructs"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        data_dir: Optional[Path] = None,
    ):
        """Initialize the construct ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            data_dir: Directory containing reference construct files.
                Defaults to the project data/reference directory.
        """
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data" / "reference"

    def fetch(
        self,
        reference_file: Optional[str] = None,
        include_fda_seed: bool = True,
    ) -> List[Dict[str, Any]]:
        """Fetch CAR construct data from reference files and/or seed data.

        Args:
            reference_file: Path to a JSON or CSV file containing additional
                construct definitions.  If None, only FDA seed data is used.
            include_fda_seed: If True, include the 6 FDA-approved CAR-T
                products as seed data.

        Returns:
            List of construct data dicts.

        """
        records = []

        if include_fda_seed:
            fda_constructs = self._get_fda_approved_constructs()
            records.extend(c.model_dump() for c in fda_constructs)

        if reference_file:
            ref_path = Path(reference_file)
            if ref_path.suffix.lower() == ".json":
                with open(ref_path, "r") as f:
                    file_data = json.load(f)
                records.extend(file_data)
            elif ref_path.suffix.lower() == ".csv":
                with open(ref_path, "r", newline="") as f:
                    reader = csv.DictReader(f)
                    records.extend(list(reader))

        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[CARConstruct]:
        """Parse construct data dicts into CARConstruct models.

        Args:
            raw_data: List of dicts from fetch(), each containing construct
                fields matching the CARConstruct model.

        Returns:
            List of validated CARConstruct model instances.

        """
        records = []
        for data in raw_data:
            try:
                record = CARConstruct(**data)
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse construct {data.get('id', '?')}: {e}")
                continue
        return records

    @staticmethod
    def _get_fda_approved_constructs() -> List[CARConstruct]:
        """Return the 6 FDA-approved CAR-T products as CARConstruct objects.

        Data is current as of February 2026.  All 6 products target
        hematologic malignancies (B-ALL, DLBCL, Multiple Myeloma, FL).

        Returns:
            List of 6 CARConstruct instances with real clinical data.
        """
        return [
            CARConstruct(
                id="fda-tisagenlecleucel",
                name="Tisagenlecleucel (Kymriah)",
                text_summary=(
                    "Tisagenlecleucel (Kymriah) is a CD19-directed CAR-T cell therapy "
                    "developed by Novartis. It was the first CAR-T product approved by "
                    "the FDA (August 2017) for pediatric and young adult patients with "
                    "relapsed/refractory B-cell acute lymphoblastic leukemia (B-ALL). "
                    "Subsequently approved for adult relapsed/refractory diffuse large "
                    "B-cell lymphoma (DLBCL) and follicular lymphoma (FL). Uses a "
                    "lentiviral vector with 4-1BB costimulatory domain. The ELIANA trial "
                    "showed 82% overall remission rate in pediatric B-ALL."
                ),
                target_antigen="CD19",
                scfv_origin="FMC63 murine anti-CD19",
                costimulatory_domain="4-1BB (CD137)",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="CD8-alpha hinge and transmembrane",
                vector_type="lentiviral",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (58% any grade, 22% grade 3+), neurological events (39%), B-cell aplasia, hypogammaglobulinemia",
            ),
            CARConstruct(
                id="fda-axicabtagene-ciloleucel",
                name="Axicabtagene ciloleucel (Yescarta)",
                text_summary=(
                    "Axicabtagene ciloleucel (Yescarta) is a CD19-directed CAR-T cell "
                    "therapy developed by Kite Pharma / Gilead Sciences. FDA approved "
                    "October 2017 for adult relapsed/refractory large B-cell lymphoma "
                    "after two or more lines of systemic therapy. Also approved for "
                    "follicular lymphoma (FL). Uses a retroviral vector with CD28 "
                    "costimulatory domain. The ZUMA-1 trial demonstrated 83% overall "
                    "response rate and 58% complete response rate in DLBCL."
                ),
                target_antigen="CD19",
                scfv_origin="FMC63 murine anti-CD19",
                costimulatory_domain="CD28",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="CD28 hinge and transmembrane",
                vector_type="retroviral (gamma)",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (93% any grade, 13% grade 3+), neurological events (64%, 28% grade 3+), B-cell aplasia, cytopenias",
            ),
            CARConstruct(
                id="fda-brexucabtagene-autoleucel",
                name="Brexucabtagene autoleucel (Tecartus)",
                text_summary=(
                    "Brexucabtagene autoleucel (Tecartus) is a CD19-directed CAR-T cell "
                    "therapy developed by Kite Pharma / Gilead Sciences. FDA approved "
                    "July 2020 for adult relapsed/refractory mantle cell lymphoma (MCL). "
                    "Also approved for adult B-cell acute lymphoblastic leukemia (B-ALL) "
                    "in October 2021. Uses the same CAR construct as Yescarta with CD28 "
                    "costimulatory domain but a different manufacturing process that "
                    "includes T cell enrichment to remove circulating tumor cells. "
                    "The ZUMA-2 trial showed 87% overall response rate in MCL."
                ),
                target_antigen="CD19",
                scfv_origin="FMC63 murine anti-CD19",
                costimulatory_domain="CD28",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="CD28 hinge and transmembrane",
                vector_type="retroviral (gamma)",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (91% any grade, 15% grade 3+), neurological events (63%, 31% grade 3+), cytopenias, infections",
            ),
            CARConstruct(
                id="fda-lisocabtagene-maraleucel",
                name="Lisocabtagene maraleucel (Breyanzi)",
                text_summary=(
                    "Lisocabtagene maraleucel (Breyanzi) is a CD19-directed CAR-T cell "
                    "therapy developed by Juno Therapeutics / Bristol Myers Squibb. "
                    "FDA approved February 2021 for adult relapsed/refractory large "
                    "B-cell lymphoma after two or more lines of systemic therapy. "
                    "Uniquely administered as a defined composition of CD8+ and CD4+ "
                    "CAR-T cells at a 1:1 ratio. Uses a lentiviral vector with 4-1BB "
                    "costimulatory domain. The TRANSCEND NHL 001 trial demonstrated "
                    "73% overall response rate and 53% complete response rate."
                ),
                target_antigen="CD19",
                scfv_origin="FMC63 murine anti-CD19",
                costimulatory_domain="4-1BB (CD137)",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="IgG4 hinge, CD28 transmembrane",
                vector_type="lentiviral",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (46% any grade, 4% grade 3+), neurological events (35%, 12% grade 3+), cytopenias, hypogammaglobulinemia",
            ),
            CARConstruct(
                id="fda-idecabtagene-vicleucel",
                name="Idecabtagene vicleucel (Abecma)",
                text_summary=(
                    "Idecabtagene vicleucel (Abecma) is a BCMA-directed CAR-T cell "
                    "therapy developed by Celgene / Bristol Myers Squibb and bluebird "
                    "bio. FDA approved March 2021 for adult relapsed/refractory multiple "
                    "myeloma after four or more prior lines of therapy including an "
                    "immunomodulatory agent, a proteasome inhibitor, and an anti-CD38 "
                    "antibody. First BCMA-targeting CAR-T product approved. Uses a "
                    "lentiviral vector with 4-1BB costimulatory domain. The KarMMa trial "
                    "showed 73% overall response rate and 33% complete response rate."
                ),
                target_antigen="BCMA",
                scfv_origin="Murine anti-BCMA",
                costimulatory_domain="4-1BB (CD137)",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="CD8-alpha hinge and transmembrane",
                vector_type="lentiviral",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (84% any grade, 5% grade 3+), neurological events (18%, 3% grade 3+), cytopenias (97%), infections",
            ),
            CARConstruct(
                id="fda-ciltacabtagene-autoleucel",
                name="Ciltacabtagene autoleucel (Carvykti)",
                text_summary=(
                    "Ciltacabtagene autoleucel (Carvykti) is a BCMA-directed CAR-T cell "
                    "therapy developed by Janssen / Legend Biotech. FDA approved February "
                    "2022 for adult relapsed/refractory multiple myeloma after four or "
                    "more prior lines of therapy. Features a unique dual-epitope binding "
                    "domain consisting of two BCMA-targeting single-domain antibodies "
                    "(VHH nanobodies) that confer high avidity binding. Uses a lentiviral "
                    "vector with 4-1BB costimulatory domain. The CARTITUDE-1 trial "
                    "demonstrated 98% overall response rate and 83% stringent complete "
                    "response rate, the highest response rates among approved CAR-T products."
                ),
                target_antigen="BCMA",
                scfv_origin="Dual VHH nanobody (llama-derived) anti-BCMA",
                costimulatory_domain="4-1BB (CD137)",
                signaling_domain="CD3-zeta",
                generation=CARGeneration.SECOND,
                hinge_tm="CD8-alpha hinge and transmembrane",
                vector_type="lentiviral",
                fda_status=FDAStatus.APPROVED,
                known_toxicities="CRS (95% any grade, 4% grade 3+), neurotoxicity (21%, 9% grade 3+ including parkinsonian-like movement disorders), cytopenias, infections",
            ),
        ]

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full construct ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'cart_constructs').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (reference_file, include_fda_seed).

        Returns:
            Total number of records ingested.

        """
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
